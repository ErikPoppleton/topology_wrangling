"""Write AMBER (Lipid21)-style head/tail index groups for CHARMM/SLipids lipids.

AMBER's Lipid21 splits each phospholipid into three residues -- tail, head,
tail -- while CHARMM36 and SLipids keep the whole lipid in a single residue.
This reads a CHARMM- or SLipids-style .gro and writes a GROMACS .ndx whose
groups hold exactly the atoms the corresponding AMBER residues would, named
after them:

    lipid   head group   tail group
    DOPC    PC   (38)    OL   (2 x 50)
    DMPC    PC   (38)    MY   (2 x 40)
    DMPG    PGR  (31)    MY   (2 x 40)

The ester carbonyls belong to the head, as in Lipid21; each tail runs from
C22/C32 to the terminal methyl, hydrogens included.  A system holding several
of these lipids gets one group per AMBER residue name (PC, PGR and MY for
DMPC + DMPG), and everything else in the file is ignored.

Each hydrogen goes with the heavy atom it is bonded to.  Give the lipid's .itp
and that is read from the topology's bonds; without one it is the nearest
heavy atom, by minimum-image distance in the (rectangular or triclinic) box.
Atoms are matched to a topology by name, with the same --residue/--alias
options as reorder_gro.py.

Index numbers are 1-based positions in the .gro (GROMACS atom indices), not
the wrapped atom-number column.

Examples:
    amber_style_index.py system.gro                    # -> system_amber.ndx
    amber_style_index.py system.gro DOPC.itp -o lipids.ndx
    amber_style_index.py system.gro -o index.ndx --append
"""

from __future__ import annotations

import argparse

from ..gro import GroFrame
from ..itp import ItpFile
from ..lipid21 import AMBER_SIZE, LIPIDS, index_groups
from ..ndx import write_ndx
from ..residues import ResidueTopology
from .common import add_alias_options, alias_table, default_out, report, run


def build_parser():
    parser = argparse.ArgumentParser(
        prog="amber_style_index.py",
        description=__doc__.splitlines()[0],
        epilog="\n".join(__doc__.splitlines()[1:]),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("gro", help="CHARMM- or SLipids-style .gro")
    parser.add_argument("itp", nargs="*",
                        help="lipid topologies to take bonds from (optional; "
                             "lipids without one are split by distance)")
    parser.add_argument("-o", "--output",
                        help="output .ndx (default: <name>_amber.ndx)")
    parser.add_argument("--append", action="store_true",
                        help="append the groups to an existing .ndx instead "
                             "of overwriting it")
    add_alias_options(parser)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    itps = [ItpFile.read(path) for path in args.itp]
    table = alias_table(args, itps)
    topologies = [ResidueTopology.from_itp(itp, table) for itp in itps]

    frame = GroFrame.read(args.gro)
    groups, counts = index_groups(frame, topologies)

    bonded = {top.gro_resname: top.itp.path for top in topologies}
    for top in topologies:
        if top.gro_resname not in counts:
            report("%s: no %s in %s; not used"
                   % (top.itp.path, top.gro_resname, args.gro))
    for resname, n in counts.items():
        _, head, tail = LIPIDS[resname]
        report("%-4s x %d -> %s (%d atoms/lipid) + %s (2 x %d atoms/lipid), %s"
               % (resname, n, head, AMBER_SIZE[head], tail, AMBER_SIZE[tail],
                  "bonds from %s" % bonded[resname] if resname in bonded
                  else "hydrogens by distance"))

    out_path = args.output or default_out(args.gro, ".ndx", tag="_amber")
    write_ndx(out_path, groups, append=args.append)
    report("%s %s: %s" % ("appended to" if args.append else "wrote", out_path,
                          ", ".join("%s (%d)" % (name, len(ix))
                                    for name, ix in groups.items())))
    return 0


def cli():
    run(main)


if __name__ == "__main__":
    cli()
