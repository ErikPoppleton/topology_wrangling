#!/usr/bin/env python3
"""Re-order the atoms of a GROMACS .itp so every hydrogen follows its heavy atom.

GROMACS is noticeably faster when hydrogens are listed immediately after the
heavy atom they are bonded to instead of in blocks.  This script

  1. parses the [ atoms ] and [ bonds ] directives,
  2. builds the connectivity graph from [ bonds ],
  3. re-orders [ atoms ] so that each heavy atom is followed by its hydrogens
     (heavy atoms keep their original relative order, as do the hydrogens of a
     given heavy atom),
  4. rewrites every directive that refers to atom indices (bonds, pairs,
     angles, dihedrals, exclusions, constraints, restraints, virtual sites, ...)
     with the new numbering, and
  5. optionally applies the same permutation to every molecule of that residue
     in a coordinate (.gro) file, so the coordinates still match the topology.

Only the numbers change: functional forms, parameters and the order of the
atoms *inside* a bonded interaction are left untouched, so improper dihedrals
and other order-sensitive terms keep their meaning.

Usage:
    # topology only
    python reorder_itp.py DOPC.itp                        # -> DOPC_reordered.itp

    # topology + coordinates in one go
    python reorder_itp.py DOPC.itp --gro md.gro           # -> ..._reordered.{itp,gro}

    # re-order a coordinate file with a map written by an earlier run
    python reorder_itp.py --read-map map.txt --gro md.gro

Options:
    -o/--output, --gro-out   output paths
    --write-map FILE         save the new/old index table
    --resname NAME           residue to permute in the .gro (default: from the .itp)
    --cgnr {renumber,keep,per-atom}
    --no-sort                keep the original row order in the bonded directives
    --no-check               skip the atom-name check against the .gro
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict


# Directives whose leading columns are atom indices, and how many of them.
# Everything after those columns (function type, parameters) is copied verbatim.
ATOM_INDEX_COLUMNS = {
    "bonds": 2,
    "pairs": 2,
    "pairs_nb": 2,
    "angles": 3,
    "dihedrals": 4,
    "constraints": 2,
    "settles": 1,
    "position_restraints": 1,
    "distance_restraints": 2,
    "dihedral_restraints": 4,
    "orientation_restraints": 2,
    "angle_restraints": 4,
    "angle_restraints_z": 2,
    "virtual_sites1": 1,
    "virtual_sites2": 3,
    "virtual_sites3": 4,
    "virtual_sites4": 5,
    "virtual_sitesn": None,  # site + funct + variable-length list
    "exclusions": None,      # whole line is atom indices
}

HYDROGEN_MASS_MAX = 3.5  # covers H and deuterium-ish masses


class ItpError(Exception):
    """Anything wrong with the input that the user needs to fix."""


# --------------------------------------------------------------------------
# .itp parsing
# --------------------------------------------------------------------------

class Line:
    """One physical line of the file, tagged with the directive it belongs to."""

    __slots__ = ("raw", "directive", "data", "comment")

    def __init__(self, raw, directive):
        self.raw = raw
        self.directive = directive
        self.data = None      # list of tokens, for data lines only
        self.comment = ""     # trailing ';...' comment, if any

    def is_data(self):
        return self.data is not None


def split_comment(line):
    """Return (code, comment) where comment keeps its leading ';'."""
    idx = line.find(";")
    if idx == -1:
        return line, ""
    return line[:idx], line[idx:]


def parse_itp(path):
    """Return the file as a list of Line objects."""
    lines = []
    directive = None
    with open(path) as fh:
        for raw in fh:
            raw = raw.rstrip("\n")
            code, comment = split_comment(raw)
            stripped = code.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                directive = stripped[1:-1].strip().lower()
                lines.append(Line(raw, directive))
                continue
            entry = Line(raw, directive)
            if stripped and not stripped.startswith("#"):
                entry.data = stripped.split()
                entry.comment = comment
            lines.append(entry)
    return lines


def is_hydrogen(atom):
    """Decide whether an [ atoms ] record is a hydrogen.

    Mass is the reliable signal; when the mass column is absent we fall back to
    the atom type / atom name, which in the Slipids and CHARMM naming schemes
    always starts with an H for hydrogens.
    """
    mass = atom["mass"]
    if mass is not None:
        return mass <= HYDROGEN_MASS_MAX
    for field in ("type", "name"):
        value = atom[field]
        if value:
            return value[0].upper() == "H"
    return False


def read_atoms(lines):
    """Collect the [ atoms ] records in file order."""
    atoms = []
    for entry in lines:
        if entry.directive != "atoms" or not entry.is_data():
            continue
        tok = entry.data
        if len(tok) < 6:
            raise ItpError("malformed [ atoms ] line: %r" % entry.raw)
        try:
            mass = float(tok[7]) if len(tok) > 7 else None
        except ValueError:
            mass = None
        atoms.append(
            {
                "nr": int(tok[0]),
                "type": tok[1],
                "resnr": tok[2],
                "residue": tok[3],
                "name": tok[4],
                "cgnr": tok[5],
                "charge": tok[6] if len(tok) > 6 else "",
                "mass": mass,
                "extra": tok[8:],
                "comment": entry.comment,
            }
        )
    if not atoms:
        raise ItpError("no [ atoms ] directive found")
    return atoms


def read_bond_graph(lines, natoms):
    """Adjacency lists built from [ bonds ] (1-based atom indices)."""
    graph = defaultdict(list)
    for entry in lines:
        if entry.directive != "bonds" or not entry.is_data():
            continue
        ai, aj = int(entry.data[0]), int(entry.data[1])
        for a in (ai, aj):
            if not 1 <= a <= natoms:
                raise ItpError(
                    "bond refers to atom %d, which is outside [ atoms ] (1-%d)"
                    % (a, natoms)
                )
        graph[ai].append(aj)
        graph[aj].append(ai)
    return graph


# --------------------------------------------------------------------------
# the permutation
# --------------------------------------------------------------------------

def build_order(atoms, graph):
    """Return (order, hydrogen_flags); order[new-1] == old 1-based index.

    Heavy atoms keep their original relative order; each is immediately
    followed by the hydrogens bonded to it, in their original order.  A
    hydrogen bonded to several heavy atoms is attached to the first one, and a
    hydrogen with no heavy neighbour (a bare or unbonded H) is left in place.
    """
    hydrogen = [is_hydrogen(a) for a in atoms]

    attached = defaultdict(list)   # heavy atom (old index) -> its hydrogens
    placed = set()
    for idx in range(1, len(atoms) + 1):
        if not hydrogen[idx - 1]:
            continue
        heavies = [n for n in graph.get(idx, ()) if not hydrogen[n - 1]]
        if heavies:
            attached[min(heavies)].append(idx)
            placed.add(idx)

    order = []
    for idx in range(1, len(atoms) + 1):
        if hydrogen[idx - 1]:
            if idx not in placed:  # orphan hydrogen: keep its original slot
                order.append(idx)
            continue
        order.append(idx)
        order.extend(sorted(attached[idx]))

    if sorted(order) != list(range(1, len(atoms) + 1)):
        raise AssertionError("re-ordering did not produce a permutation")
    return order, hydrogen


# --------------------------------------------------------------------------
# the index map file
# --------------------------------------------------------------------------

def write_map(path, order, atoms, resname):
    with open(path, "w") as fh:
        fh.write("; atom index map written by reorder_itp.py\n")
        fh.write("; resname %s\n" % resname)
        fh.write("; natoms %d\n" % len(order))
        fh.write(";  new    old  name\n")
        for pos, old in enumerate(order, start=1):
            fh.write("%6d %6d  %s\n" % (pos, old, atoms[old - 1]["name"]))


def read_map(path):
    """Return (order, names, resname); names is in the ORIGINAL atom order."""
    order, pairs, resname = [], [], None
    for raw in open(path):
        line = raw.strip()
        if line.startswith(";"):
            tok = line.lstrip(";").split()
            if len(tok) >= 2 and tok[0] == "resname":
                resname = tok[1]
            continue
        if not line:
            continue
        tok = line.split()
        if len(tok) < 2:
            raise ItpError("malformed map line: %r" % raw.rstrip())
        new, old = int(tok[0]), int(tok[1])
        name = tok[2] if len(tok) > 2 else ""
        order.append((new, old))
        pairs.append((old, name))
    order.sort()
    if [n for n, _ in order] != list(range(1, len(order) + 1)):
        raise ItpError("map file does not list new indices 1..N exactly once")
    old_order = [o for _, o in order]
    if sorted(old_order) != list(range(1, len(old_order) + 1)):
        raise ItpError("map file does not list old indices 1..N exactly once")
    names = [""] * len(old_order)
    for old, name in pairs:
        names[old - 1] = name
    return old_order, names, resname


# --------------------------------------------------------------------------
# .itp output
# --------------------------------------------------------------------------

def format_atom(atom, new_nr, new_cgnr):
    fields = [
        "%6d" % new_nr,
        "%11s" % atom["type"],
        "%7s" % atom["resnr"],
        "%7s" % atom["residue"],
        "%7s" % atom["name"],
        "%7s" % new_cgnr,
        "%11s" % atom["charge"],
    ]
    if atom["mass"] is not None:
        fields.append("%11s" % ("%.4f" % atom["mass"]).rstrip("0").rstrip("."))
    fields.extend("%9s" % t for t in atom["extra"])
    line = "".join(fields)
    if atom["comment"]:
        line += "   " + atom["comment"]
    return line


def remap_indices(tokens, directive, new_index):
    """Return tokens with the atom-index columns replaced by the new numbers."""
    ncols = ATOM_INDEX_COLUMNS.get(directive, 0)
    if ncols is None:
        if directive == "exclusions":
            ncols = len(tokens)
        elif directive == "virtual_sitesn":
            # site funct a1 a2 ... : index columns are 0 and 2..end
            out = list(tokens)
            out[0] = str(new_index[int(tokens[0])])
            for i in range(2, len(out)):
                out[i] = str(new_index[int(tokens[i])])
            return out
        else:
            ncols = 0
    out = list(tokens)
    for i in range(min(ncols, len(out))):
        out[i] = str(new_index[int(out[i])])
    return out


def format_bonded(tokens, ncols):
    """Right-align the atom columns, keep the rest loosely spaced."""
    parts = ["%6s" % t for t in tokens[:ncols]]
    parts.extend("%6s" % t for t in tokens[ncols : ncols + 1])   # funct
    parts.extend("%14s" % t for t in tokens[ncols + 1 :])        # parameters
    return "".join(parts) + " "


def rewrite_itp(lines, atoms, order, new_index, cgnr_mode, do_sort):
    """Produce the output lines of the re-ordered topology."""
    out = []
    i = 0
    n = len(lines)

    # Remap the charge-group numbers so they stay ascending in the new order.
    cgnr_map = {}
    if cgnr_mode == "renumber":
        for old in order:
            cg = atoms[old - 1]["cgnr"]
            if cg not in cgnr_map:
                cgnr_map[cg] = str(len(cgnr_map) + 1)

    def cgnr_for(old, pos):
        if cgnr_mode == "keep":
            return atoms[old - 1]["cgnr"]
        if cgnr_mode == "per-atom":
            return str(pos)
        return cgnr_map[atoms[old - 1]["cgnr"]]

    while i < n:
        entry = lines[i]
        directive = entry.directive

        # --- [ atoms ]: emit the whole block in the new order -------------
        if directive == "atoms" and entry.is_data():
            while i < n and lines[i].directive == "atoms" and lines[i].is_data():
                i += 1
            for pos, old in enumerate(order, start=1):
                out.append(format_atom(atoms[old - 1], pos, cgnr_for(old, pos)))
            continue

        # --- bonded / index-carrying directives --------------------------
        if directive in ATOM_INDEX_COLUMNS and entry.is_data():
            block = []
            while i < n and lines[i].directive == directive and lines[i].is_data():
                block.append(lines[i])
                i += 1
            ncols = ATOM_INDEX_COLUMNS[directive]
            rows = [
                (remap_indices(item.data, directive, new_index), item.comment)
                for item in block
            ]
            if do_sort and ncols:
                rows.sort(key=lambda r: tuple(int(t) for t in r[0][:ncols]))
            for tokens, comment in rows:
                if ncols:
                    text = format_bonded(tokens, ncols)
                else:
                    text = "".join("%6s" % t for t in tokens) + " "
                if comment:
                    text += "  " + comment
                out.append(text)
            continue

        out.append(entry.raw)
        i += 1

    return out


# --------------------------------------------------------------------------
# .gro handling
# --------------------------------------------------------------------------

def read_gro(path):
    """Return (title, atom_lines, box, trailing).  Lines keep their exact text."""
    with open(path) as fh:
        title = fh.readline().rstrip("\n")
        count_line = fh.readline()
        if not count_line:
            raise ItpError("%s: truncated, no atom count" % path)
        try:
            natoms = int(count_line.split()[0])
        except (ValueError, IndexError):
            raise ItpError("%s: second line is not an atom count: %r"
                           % (path, count_line.rstrip()))
        atom_lines = []
        for _ in range(natoms):
            raw = fh.readline()
            if not raw:
                raise ItpError(
                    "%s: file ends after %d of %d atoms"
                    % (path, len(atom_lines), natoms)
                )
            atom_lines.append(raw.rstrip("\n"))
        box = fh.readline().rstrip("\n")
        trailing = fh.read()
    if not box.strip():
        raise ItpError("%s: missing box line" % path)
    if trailing.strip():
        raise ItpError(
            "%s: extra content after the box line - multi-frame .gro files are "
            "not supported" % path
        )
    return title, atom_lines, box


def gro_resname(line):
    return line[5:10].strip()


def gro_resid(line):
    return line[0:5]


def gro_atomname(line):
    return line[10:15].strip()


def renumber_gro(line, nr):
    """Replace the atom-number column; .gro wraps it at 100000."""
    return line[:15] + "%5d" % (nr % 100000) + line[20:]


def reorder_gro(atom_lines, order, names, resname, check=True):
    """Permute the atoms of every `resname` molecule; return (lines, nmol)."""
    nat = len(order)
    out = []
    nmol = 0
    warned_resid = False
    i = 0
    total = len(atom_lines)

    while i < total:
        line = atom_lines[i]
        if len(line) < 20:
            raise ItpError("short .gro atom line at atom %d: %r" % (i + 1, line))
        if gro_resname(line) != resname:
            out.append(line)
            i += 1
            continue

        block = atom_lines[i : i + nat]
        if len(block) < nat:
            raise ItpError(
                "%s block starting at .gro atom %d has only %d atoms, expected %d"
                % (resname, i + 1, len(block), nat)
            )
        if check:
            got = [gro_atomname(l) for l in block]
            if got != names:
                bad = next(k for k in range(nat) if got[k] != names[k])
                raise ItpError(
                    "atom names in the %s molecule starting at .gro atom %d do "
                    "not match the topology: position %d is %r, expected %r.\n"
                    "The .gro must list this residue in the same order as the "
                    "original .itp.  Re-run with --no-check to override."
                    % (resname, i + 1, bad + 1, got[bad], names[bad])
                )
        if not warned_resid and len({gro_resid(l) for l in block}) != 1:
            sys.stderr.write(
                "warning: the %d-atom %s molecule at .gro atom %d spans more "
                "than one residue number\n" % (nat, resname, i + 1)
            )
            warned_resid = True

        out.extend(block[old - 1] for old in order)
        nmol += 1
        i += nat

    if nmol == 0:
        raise ItpError("no %s residues found in the .gro file" % resname)
    return [renumber_gro(l, k) for k, l in enumerate(out, start=1)], nmol


def write_gro(path, title, atom_lines, box):
    with open(path, "w") as fh:
        fh.write(title + "\n")
        fh.write("%5d\n" % len(atom_lines))
        fh.write("\n".join(atom_lines))
        fh.write("\n" + box + "\n")


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def default_out(path, suffix):
    base = path[: -len(suffix)] if path.endswith(suffix) else path
    return base + "_reordered" + suffix


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="\n".join(__doc__.splitlines()[19:]),
    )
    ap.add_argument("itp", nargs="?", help="input .itp file")
    ap.add_argument("-o", "--output", help="output .itp (default: <name>_reordered.itp)")
    ap.add_argument("--gro", help="coordinate file to permute with the same map")
    ap.add_argument("--gro-out", help="output .gro (default: <name>_reordered.gro)")
    ap.add_argument(
        "--write-map",
        "--map",
        dest="write_map",
        metavar="FILE",
        help="write the new/old index table (reusable with --read-map)",
    )
    ap.add_argument(
        "--read-map",
        metavar="FILE",
        help="take the permutation from a map file instead of an .itp; with "
        "--gro this only re-orders coordinates",
    )
    ap.add_argument(
        "--resname",
        help="residue name to permute in the .gro (default: from the .itp or map)",
    )
    ap.add_argument(
        "--cgnr",
        choices=("renumber", "keep", "per-atom"),
        default="renumber",
        help="charge-group numbers: renumber the existing groups so they ascend "
        "(default), keep them verbatim, or give every atom its own group",
    )
    ap.add_argument(
        "--no-sort",
        action="store_true",
        help="do not sort the rows of the bonded directives after re-indexing",
    )
    ap.add_argument(
        "--no-check",
        action="store_true",
        help="skip the atom-name consistency check against the .gro",
    )
    args = ap.parse_args(argv)

    if not args.itp and not args.read_map:
        ap.error("give an .itp file, or --read-map to reuse an existing map")
    if args.itp and args.read_map:
        ap.error("--read-map replaces the .itp; give one or the other")
    if args.read_map and not args.gro:
        ap.error("--read-map is only useful together with --gro")

    atoms = None

    # ---- get the permutation ------------------------------------------
    if args.read_map:
        order, names, resname = read_map(args.read_map)
        resname = args.resname or resname
        if not resname:
            ap.error("the map file records no residue name; pass --resname")
        if not args.no_check and not all(names):
            ap.error(
                "the map file has no atom-name column, so the .gro cannot be "
                "checked; pass --no-check to proceed anyway"
            )
        sys.stderr.write(
            "map: %d atoms, residue %s (from %s)\n"
            % (len(order), resname, args.read_map)
        )
    else:
        lines = parse_itp(args.itp)
        atoms = read_atoms(lines)
        if [a["nr"] for a in atoms] != list(range(1, len(atoms) + 1)):
            raise ItpError("[ atoms ] is not numbered 1..N consecutively")
        graph = read_bond_graph(lines, len(atoms))
        order, hydrogen = build_order(atoms, graph)
        names = [a["name"] for a in atoms]
        resnames = {a["residue"] for a in atoms}
        resname = args.resname or sorted(resnames)[0]
        if args.gro and not args.resname and len(resnames) > 1:
            raise ItpError(
                "the .itp covers several residue names (%s); pass --resname to "
                "say which one to permute in the .gro" % ", ".join(sorted(resnames))
            )

        # ---- write the topology ---------------------------------------
        new_index = {old: pos for pos, old in enumerate(order, start=1)}
        out_itp = args.output or default_out(args.itp, ".itp")
        out_lines = rewrite_itp(
            lines, atoms, order, new_index, args.cgnr, not args.no_sort
        )
        with open(out_itp, "w") as fh:
            fh.write("\n".join(out_lines) + "\n")

        nh = sum(hydrogen)
        moved = sum(1 for pos, old in enumerate(order, start=1) if pos != old)
        sys.stderr.write(
            "%s: %d atoms (%d hydrogens, %d heavy); %d atoms changed index\n"
            % (args.itp, len(atoms), nh, len(atoms) - nh, moved)
        )
        sys.stderr.write("wrote %s\n" % out_itp)

    if args.write_map:
        write_map(
            args.write_map,
            order,
            atoms or [{"name": n} for n in names],
            resname,
        )
        sys.stderr.write("wrote %s\n" % args.write_map)

    # ---- permute the coordinates --------------------------------------
    if args.gro:
        title, atom_lines, box = read_gro(args.gro)
        new_atoms, nmol = reorder_gro(
            atom_lines, order, names, resname, check=not args.no_check
        )
        out_gro = args.gro_out or default_out(args.gro, ".gro")
        write_gro(out_gro, title, new_atoms, box)
        sys.stderr.write(
            "%s: %d atoms, permuted %d %s molecules\n"
            % (args.gro, len(atom_lines), nmol, resname)
        )
        sys.stderr.write("wrote %s\n" % out_gro)

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ItpError as err:
        sys.stderr.write("error: %s\n" % err)
        sys.exit(1)
