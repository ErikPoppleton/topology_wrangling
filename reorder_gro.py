#!/usr/bin/env python3
"""
Reorder the atoms of every molecule of one residue in a GROMACS .gro file so
that the per-molecule atom order matches the order given by an .itp topology,
rewrite the residue and atom names to the topology's, then renumber the global
atom-index column sequentially. The output residue is then an exact match for
the .itp in both order and naming, which is what grompp wants -- the point of
the script is to make a .gro written by one tool usable with force-field files
that follow different naming conventions.

Matching is by atom name, which is unique within a molecule, so the mapping
is unambiguous. Only the residue the topology describes is touched: every
other residue stays exactly where it is, and the title, box and coordinates
(and velocities) are carried through verbatim.

The residue name comes from the .itp's [ atoms ] residue column, so nothing
about a particular molecule is built in. Where the two files disagree on a
name rather than only on the order:

    --resname SOD              # the gro calls this residue SOD, the .itp Na+
    --alias O=OH2              # ITP_NAME=GRO_NAME, repeatable

--resname names the residue as it appears *in the .gro*; the output takes the
.itp's residue name (SOD becomes Na+). An --alias likewise declares which gro
atom is which topology atom, and that atom takes the .itp's name in the output
(C2T becomes C26). Where the .itp's own residue name is unusable -- DCLE.itp
calls its residue 303 -- --out-resname sets what is written instead. Pass
--keep-names to reorder only, leaving both names alone. The atom-name check
lists the unmatched names on each side, so the aliases a file needs can be
read straight off the error.

GRO fixed-column format (per line):
    cols  0:5   residue number
    cols  5:10  residue name (rewritten to the topology's name)
    cols 10:15  atom name    (rewritten to the topology's name)
    cols 15:20  atom number  (renumbered here)
    cols 20:    coordinates (+ optional velocities) -- preserved verbatim

Usage:
    reorder_gro.py <input.gro> <topology.itp> [-o output.gro]
                   [--resname GRO_RESNAME] [--out-resname NAME]
                   [--alias ITP_NAME=GRO_NAME ...] [--keep-names]

With no -o, the output is written beside the input as <input>_reordered.gro.
"""
import os
import sys


def parse_itp_atoms(itp_path):
    """Return (residue name, atom names in [ atoms ] order) for one topology."""
    names, resnames = [], []
    in_atoms = False
    with open(itp_path) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.startswith("["):
                in_atoms = stripped.strip("[] \t").startswith("atoms")
                continue
            if not in_atoms or not stripped or stripped.startswith(";"):
                continue
            # columns: nr type resnr residue atom cgnr charge mass
            fields = stripped.split()
            if len(fields) < 5:
                continue
            resnames.append(fields[3])
            names.append(fields[4])
    if not names:
        sys.exit("%s has no [ atoms ] entries" % itp_path)
    distinct = set(resnames)
    if len(distinct) != 1:
        sys.exit("%s spans several residue names (%s); this script reorders "
                 "one residue" % (itp_path, ", ".join(sorted(distinct))))
    if len(names) != len(set(names)):
        sys.exit("%s repeats an atom name, so the order cannot be matched "
                 "by name" % itp_path)
    return resnames[0], names


def parse_args(argv):
    positional, out, resname, aliases, rename = [], None, None, {}, True
    out_resname = None
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg in ("-h", "--help"):
            sys.exit(__doc__)
        elif arg in ("-o", "--output"):
            out, i = argv[i + 1], i + 2
        elif arg == "--resname":
            resname, i = argv[i + 1], i + 2
        elif arg == "--out-resname":
            out_resname, i = argv[i + 1], i + 2
        elif arg == "--keep-names":
            rename, i = False, i + 1
        elif arg == "--alias":
            spec = argv[i + 1]
            if "=" not in spec:
                sys.exit("--alias wants ITP_NAME=GRO_NAME, got %r" % spec)
            itp_name, gro_name = spec.split("=", 1)
            aliases[itp_name] = gro_name
            i += 2
        else:
            positional.append(arg)
            i += 1
    if len(positional) != 2:
        sys.exit(__doc__)
    gro, itp = positional
    if out is None:
        root, ext = os.path.splitext(gro)
        out = root + "_reordered" + (ext or ".gro")
    return gro, itp, out, resname, out_resname, aliases, rename


def build_order(itp_names, aliases):
    """{gro atom name: position within the molecule}."""
    unknown = sorted(set(aliases) - set(itp_names))
    if unknown:
        sys.exit("--alias names not present in the topology: %s"
                 % ", ".join(unknown))
    order = {}
    for position, itp_name in enumerate(itp_names):
        order[aliases.get(itp_name, itp_name)] = position
    if len(order) != len(itp_names):
        sys.exit("the --alias mapping collides with an existing atom name")
    return order


def name_field(name, justify_left, what="atom"):
    """Format a name into one of the .gro's 5-column name fields."""
    if len(name) > 5:
        sys.exit("%s name %r does not fit the .gro's 5-column field"
                 % (what, name))
    return ("%-5s" if justify_left else "%5s") % name


def justifies_left(body, resname, lo=10, hi=15):
    """Whether this .gro left-justifies the name in columns lo:hi."""
    left = right = 0
    for line in body:
        if line[5:10].strip() != resname:
            continue
        field = line[lo:hi]
        if len(field.strip()) == 5:
            continue
        if field[0] != " ":
            left += 1
        elif field[-1] != " ":
            right += 1
    return left > right


def reorder(body, resname, size, order, itp_names, justify_left,
            out_resname=None, res_justify_left=True):
    """Return the reordered atom lines and the number of molecules touched."""
    out, count, i = [], 0, 0
    while i < len(body):
        if body[i][5:10].strip() != resname:
            out.append(body[i])
            i += 1
            continue

        block = body[i:i + size]
        resnum = block[0][0:5]
        if len(block) != size or any(line[5:10].strip() != resname
                                     or line[0:5] != resnum for line in block):
            sys.exit("the %s residue starting at atom %d is not %d contiguous "
                     "atoms of a single residue" % (resname, i + 1, size))

        names = set(line[10:15].strip() for line in block)
        if names != set(order):
            sys.exit("the %s residue starting at atom %d does not match the "
                     "topology.\n  only in .itp: %s\n  only in .gro: %s\n"
                     "Map differing names with --alias ITP_NAME=GRO_NAME."
                     % (resname, i + 1,
                        ", ".join(sorted(set(order) - names)) or "-",
                        ", ".join(sorted(names - set(order))) or "-"))

        block = sorted(block, key=lambda line: order[line[10:15].strip()])
        if itp_names is not None:
            block = [line[:10] + name_field(itp_names[position], justify_left)
                     + line[15:] for position, line in enumerate(block)]
        if out_resname is not None:
            field = name_field(out_resname, res_justify_left, "residue")
            block = [line[:5] + field + line[10:] for line in block]
        out.extend(block)
        count += 1
        i += size
    return out, count


def main():
    (gro, itp, out_path, resname, out_resname, aliases,
     rename) = parse_args(sys.argv[1:])

    itp_resname, itp_names = parse_itp_atoms(itp)
    order = build_order(itp_names, aliases)
    if resname is None:
        resname = itp_resname
    if not rename:
        out_resname = None
    elif out_resname is None:
        out_resname = itp_resname
    print("%s: %d atoms per molecule%s%s"
          % (resname, len(itp_names),
             "" if out_resname in (None, resname)
             else " -> residue %s" % out_resname,
             "" if not aliases else ", %s %s" % (
                 "renaming" if rename else "aliases",
                 ", ".join("%s->%s" % (gro_name, itp_name) if rename
                           else "%s->%s" % (itp_name, gro_name)
                           for itp_name, gro_name in sorted(aliases.items())))))

    with open(gro) as handle:
        lines = handle.read().splitlines(True)
    title, natoms = lines[0], int(lines[1])
    body, box = lines[2:2 + natoms], lines[2 + natoms:]
    if len(body) != natoms:
        sys.exit("%s claims %d atoms but holds %d" % (gro, natoms, len(body)))

    body, count = reorder(body, resname, len(itp_names), order,
                          itp_names if rename else None,
                          justifies_left(body, resname, 10, 15),
                          out_resname if out_resname != resname else None,
                          justifies_left(body, resname, 5, 10))
    if not count:
        sys.exit("no %s residues in %s; pass --resname if the gro names the "
                 "residue differently" % (resname, gro))

    # renumber the global atom-index column, wrapping at 100000 as GROMACS does
    for i, line in enumerate(body):
        body[i] = line[:15] + "%5d" % ((i + 1) % 100000) + line[20:]

    with open(out_path, "w") as handle:
        handle.write(title)
        handle.write("%d\n" % natoms)
        handle.writelines(body)
        handle.writelines(box)

    print("reordered%s %d %s molecules; wrote %s (%d atoms)"
          % ("" if not rename else " and renamed", count, resname,
             out_path, natoms))


if __name__ == "__main__":
    main()
