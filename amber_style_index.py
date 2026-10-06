#!/usr/bin/env python3
"""Write AMBER (Lipid17)-style head/tail index groups for CHARMM/SLipids lipids.

AMBER's Lipid17 splits each phospholipid into three residues -- tail, head,
tail -- while CHARMM36 and SLipids keep the whole lipid in a single residue.
This script reads a CHARMM- or SLipids-style ``.gro`` and writes a GROMACS
``.ndx`` whose groups hold exactly the atoms the corresponding AMBER residues
would, named after them:

    lipid   head group   tail group
    DOPC    PC           OL
    DMPC    PC           MY
    DMPG    PGR          MY

The boundary follows Lipid17, where the ester carbonyls belong to the *head*:
the head is phosphate + glycerol + both ``C=O`` groups (CHARMM/SLipids
``C21``/``O21``/``O22`` and ``C31``/``O31``/``O32``), and each tail starts at
the first methylene past the carbonyl (``C22``/``C32``) and runs to the
terminal methyl, hydrogens included.  Per lipid that gives 38 head atoms
(PC) or 31 (PGR), and 50 (OL) or 40 (MY) atoms per tail -- the same counts as
the AMBER residues.

Tail carbons are picked by name (``C2{n}``/``C3{n}``, n >= 2); hydrogens are
assigned to the nearest heavy atom of the same lipid (minimum image), so the
result does not depend on the force field's hydrogen naming.  A system holding
several of these lipids gets one group per AMBER residue name (e.g. ``PC``,
``PGR`` and ``MY`` for DMPC + DMPG).  Everything that is not one of the three
lipids is ignored.

Index numbers are 1-based positions in the ``.gro`` (i.e. GROMACS atom
indices), not the wrapped atom-number column.

Usage:
    python3 amber_style_index.py system.gro                 # -> system_amber.ndx
    python3 amber_style_index.py system.gro -o lipids.ndx
    python3 amber_style_index.py system.gro -o index.ndx --append

Standard library only.
"""

import argparse
import sys
from collections import OrderedDict
from pathlib import Path

# lipid resname -> (number of carbons per acyl chain, AMBER head, AMBER tail)
LIPIDS = {
    "DOPC": (18, "PC", "OL"),
    "DMPC": (14, "PC", "MY"),
    "DMPG": (14, "PGR", "MY"),
}

# Expected atoms per lipid in each AMBER residue (one tail residue per chain).
AMBER_SIZE = {"PC": 38, "PGR": 31, "OL": 50, "MY": 40}


def read_gro(path):
    """Return ``(atoms, box)``; atoms are ``(resid, resname, name, xyz)``."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[1])
    atom_lines = lines[2:2 + n]
    box = [float(v) for v in lines[2 + n].split()[:3]]

    # Coordinate field width from the decimal-point spacing, so files written
    # with extra precision parse too (and touching negative numbers don't).
    first = atom_lines[0]
    dot1 = first.index(".", 20)
    width = first.index(".", dot1 + 1) - dot1

    atoms = []
    for line in atom_lines:
        xyz = [float(line[20 + i * width:20 + (i + 1) * width]) for i in range(3)]
        atoms.append((line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), xyz))
    return atoms, box


def split_molecules(atoms):
    """Group consecutive atom indices into residues.

    A new residue starts when the resid or resname changes, or when an atom
    name repeats (guards against resid wrap-around and duplicated numbering).
    """
    residues, current, seen, key = [], [], set(), None
    for i, (resid, resname, name, _) in enumerate(atoms):
        if current and ((resid, resname) != key or name in seen):
            residues.append(current)
            current, seen = [], set()
        current.append(i)
        seen.add(name)
        key = (resid, resname)
    if current:
        residues.append(current)
    return residues


def min_image_dist2(a, b, box):
    d2 = 0.0
    for k in range(3):
        d = a[k] - b[k]
        if box[k] > 0:
            d -= box[k] * round(d / box[k])
        d2 += d * d
    return d2


def split_lipid(atoms, residue, box, n_carbons):
    """Return ``(head, tail)`` atom-index lists for one CHARMM/SLipids lipid."""
    names = {atoms[i][2]: i for i in residue}
    tail_carbons = ["C%d%d" % (chain, c) for chain in (2, 3) for c in range(2, n_carbons + 1)]
    missing = [nm for nm in tail_carbons if nm not in names]
    if missing:
        raise ValueError("lipid at atom %d lacks tail carbons %s"
                         % (residue[0] + 1, ", ".join(missing)))
    tail_heavy = set(names[nm] for nm in tail_carbons)

    heavy = [i for i in residue if not atoms[i][2].startswith("H")]
    in_tail = set(tail_heavy)
    for i in residue:
        if atoms[i][2].startswith("H"):
            parent = min(heavy, key=lambda j: min_image_dist2(atoms[i][3], atoms[j][3], box))
            if parent in tail_heavy:
                in_tail.add(i)

    head = [i for i in residue if i not in in_tail]
    tail = [i for i in residue if i in in_tail]
    return head, tail


def write_ndx(path, groups, append=False):
    with open(path, "a" if append else "w") as fh:
        for name, indices in groups.items():
            fh.write("[ %s ]\n" % name)
            for k in range(0, len(indices), 15):
                fh.write(" ".join("%d" % (i + 1) for i in indices[k:k + 15]) + "\n")
            fh.write("\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gro", help="CHARMM- or SLipids-style .gro")
    ap.add_argument("-o", "--output", help="output .ndx (default: <gro stem>_amber.ndx)")
    ap.add_argument("--append", action="store_true",
                    help="append the groups to an existing .ndx instead of overwriting")
    args = ap.parse_args()

    atoms, box = read_gro(args.gro)
    groups = OrderedDict()
    counts = OrderedDict()
    for residue in split_molecules(atoms):
        resname = atoms[residue[0]][1]
        if resname not in LIPIDS:
            continue
        n_carbons, head_name, tail_name = LIPIDS[resname]
        head, tail = split_lipid(atoms, residue, box, n_carbons)
        if len(head) != AMBER_SIZE[head_name] or len(tail) != 2 * AMBER_SIZE[tail_name]:
            sys.exit("%s at atom %d split %d/%d, expected %d/%d -- unexpected topology"
                     % (resname, residue[0] + 1, len(head), len(tail),
                        AMBER_SIZE[head_name], 2 * AMBER_SIZE[tail_name]))
        groups.setdefault(head_name, []).extend(head)
        groups.setdefault(tail_name, []).extend(tail)
        counts[resname] = counts.get(resname, 0) + 1

    if not groups:
        sys.exit("no %s residues found in %s" % ("/".join(LIPIDS), args.gro))

    out = args.output or str(Path(args.gro).with_suffix("")) + "_amber.ndx"
    write_ndx(out, groups, append=args.append)
    for resname, n in counts.items():
        _, head_name, tail_name = LIPIDS[resname]
        print("%-4s x %d -> %s (%d atoms/lipid) + %s (2 x %d atoms/lipid)"
              % (resname, n, head_name, AMBER_SIZE[head_name],
                 tail_name, AMBER_SIZE[tail_name]))
    print("wrote %s: %s" % (out, ", ".join("%s (%d)" % (g, len(ix)) for g, ix in groups.items())))


if __name__ == "__main__":
    main()
