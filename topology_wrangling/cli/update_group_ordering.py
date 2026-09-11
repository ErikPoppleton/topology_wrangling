"""Re-order the atoms of a GROMACS .itp so every hydrogen follows its heavy atom.

GROMACS update groups require that the hydrogens of a heavy atom be listed
immediately after it rather than in blocks.  This tool

  1. parses the [ atoms ] and [ bonds ] directives,
  2. builds the connectivity graph from [ bonds ],
  3. re-orders [ atoms ] so each heavy atom is followed by its hydrogens (heavy
     atoms keep their original relative order, as do the hydrogens of a given
     heavy atom),
  4. rewrites every directive that refers to atom indices (bonds, pairs,
     angles, dihedrals, exclusions, constraints, restraints, virtual sites...)
     with the new numbering, and
  5. optionally applies the same permutation to a coordinate file.

Only the numbers change: functional forms, parameters and the order of the
atoms *inside* a bonded interaction are left untouched, so improper dihedrals
and other order-sensitive terms keep their meaning.

Examples:
    update_group_ordering.py DOPC.itp                     # -> DOPC_reordered.itp
    update_group_ordering.py DOPC.itp --gro md.gro        # topology and coordinates
    update_group_ordering.py DOPC.itp --write-map map.txt
    update_group_ordering.py --read-map map.txt --gro md.gro
"""

from __future__ import annotations

import argparse

from ..errors import TopologyError
from ..gro import GroFrame
from ..itp import ItpFile
from ..mapping import AtomOrder
from ..ordering import update_group_order
from .common import default_out, report, run


def build_parser():
    parser = argparse.ArgumentParser(
        prog="update_group_ordering.py",
        description=__doc__.splitlines()[0],
        epilog="\n".join(__doc__.splitlines()[1:]),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("itp", nargs="?", help="input .itp file")
    parser.add_argument("-o", "--output",
                        help="output .itp (default: <name>_reordered.itp)")
    parser.add_argument("--gro", help="coordinate file to permute the same way")
    parser.add_argument("--gro-out",
                        help="output .gro (default: <name>_reordered.gro)")
    parser.add_argument("--write-map", "--map", dest="write_map", metavar="FILE",
                        help="write the new/old index table (reusable with "
                             "--read-map)")
    parser.add_argument("--read-map", metavar="FILE",
                        help="take the permutation from a map file instead of "
                             "an .itp; with --gro this only re-orders "
                             "coordinates")
    parser.add_argument("--resname",
                        help="residue to permute in the .gro (default: the "
                             "one the topology or map names)")
    parser.add_argument("--cgnr", choices=("renumber", "keep", "per-atom"),
                        default="renumber",
                        help="charge-group numbers: renumber the existing "
                             "groups so they ascend (default), keep them "
                             "verbatim, or give every atom its own group")
    parser.add_argument("--no-sort", action="store_true",
                        help="do not sort the rows of the bonded directives "
                             "after re-indexing")
    parser.add_argument("--no-check", action="store_true",
                        help="skip the atom-name check against the .gro")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.itp and not args.read_map:
        parser.error("give an .itp file, or --read-map to reuse an existing map")
    if args.itp and args.read_map:
        parser.error("--read-map replaces the .itp; give one or the other")
    if args.read_map and not args.gro:
        parser.error("--read-map is only useful together with --gro")

    if args.read_map:
        order = AtomOrder.read_map(args.read_map)
        if args.resname:
            order.resname = args.resname
        if not order.resname:
            parser.error("the map file records no residue name; pass --resname")
        if not args.no_check and not all(order.names):
            parser.error("the map file has no atom-name column, so the .gro "
                         "cannot be checked; pass --no-check to proceed anyway")
        report("map: %d atoms, residue %s (from %s)"
               % (order.natoms, order.resname, args.read_map))
    else:
        order = reorder_topology(args)

    if args.write_map:
        order.write_map(args.write_map, source=args.itp or args.read_map)
        report("wrote %s" % args.write_map)

    if args.gro:
        permute_coordinates(args, order)
    return 0


def reorder_topology(args):
    """Re-order the topology, write it, and return the permutation used."""
    itp = ItpFile.read(args.itp)
    resname = itp.residue_name(args.resname) if args.resname or args.gro \
        else itp.residue_names[0]
    order, hydrogen = update_group_order(itp.atoms, itp.bond_graph(), resname)

    out_path = args.output or default_out(args.itp, ".itp")
    itp.permute(order, args.cgnr, not args.no_sort).write(out_path)

    nh = sum(hydrogen)
    report("%s: %d atoms (%d hydrogens, %d heavy); %d atoms changed index"
           % (args.itp, itp.natoms, nh, itp.natoms - nh, order.moved()))
    report("wrote %s" % out_path)
    return order


def permute_coordinates(args, order):
    """Apply the permutation to every molecule of the residue in the .gro."""
    frame = GroFrame.read(args.gro)
    resname = args.resname or order.resname
    count = 0
    out = []
    i, total = 0, len(frame.atoms)
    while i < total:
        if frame.atoms[i].resname != resname:
            out.append(frame.atoms[i])
            i += 1
            continue
        start, block = i, frame.atoms[i:i + order.natoms]
        if len(block) != order.natoms:
            raise TopologyError(
                "%s: the %s molecule starting at atom %d has only %d atoms, "
                "expected %d"
                % (args.gro, resname, start + 1, len(block), order.natoms)
            )
        if not args.no_check:
            check_names(block, order, resname, start, args.gro)
        out.extend(order.apply(block))
        count += 1
        i += order.natoms

    if not count:
        raise TopologyError("%s: no %s residues found" % (args.gro, resname))

    frame.atoms = out
    frame.renumber()
    out_path = args.gro_out or default_out(args.gro, ".gro")
    frame.write(out_path)
    report("%s: %d atoms, permuted %d %s molecules"
           % (args.gro, total, count, resname))
    report("wrote %s" % out_path)


def check_names(block, order, resname, start, path):
    """The .gro must list the residue in the same order as the original .itp."""
    got = [a.name for a in block]
    if got == list(order.names):
        return
    bad = next(k for k in range(len(got)) if got[k] != order.names[k])
    raise TopologyError(
        "%s: the atom names of the %s molecule starting at atom %d do not "
        "match the topology: position %d is %r, expected %r.\n"
        "The .gro must list this residue in the same order as the original "
        ".itp.  Re-run with --no-check to override, or use reorder_gro.py, "
        "which matches atoms by name instead."
        % (path, resname, start + 1, bad + 1, got[bad], order.names[bad])
    )


def cli():
    run(main)


if __name__ == "__main__":
    cli()
