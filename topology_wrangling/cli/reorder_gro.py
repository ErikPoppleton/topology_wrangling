"""Make a .gro match the atom order and naming of one or more .itp topologies.

For every residue a topology describes, the atoms of each molecule are
re-ordered into the topology's order, the residue and atom names are rewritten
to the topology's, and the global atom-index column is renumbered.  The output
is then an exact match for the topologies in both order and naming, which is
what grompp wants -- the point of the tool is to make a .gro written by one
program usable with force-field files that follow different naming conventions.

Matching is by atom name, which is unique within a molecule, so the mapping is
unambiguous and the input order does not matter.  Residues that no topology
covers stay exactly where they are, and the title, box, coordinates and
velocities are carried through verbatim.

Where the two files disagree on a name rather than only on the order:

    --residue TIP3=TP3      # the .gro calls the TP3 residue TIP3
    --alias TP3:O=OH2       # within TP3, the .gro calls the atom O "OH2"
    --alias-file charmm.map # the same, read from a file

The output takes the topology's names (TIP3 becomes TP3, OH2 becomes O).  Where
a topology's own residue name is unusable -- a file that calls its residue 303 --
--out-resname sets what is written instead.  Pass --keep-names to re-order only.
The atom-name check lists the unmatched names on each side, so the aliases a
file needs can be read straight off the error.

Examples:
    reorder_gro.py system.gro DOPC.itp
    reorder_gro.py system.gro *.itp --alias-file charmm.map
"""

from __future__ import annotations

import argparse

from ..errors import TopologyError
from ..gro import GroFrame
from ..itp import ItpFile
from ..residues import AliasTable, ResidueTopology, apply_topologies
from .common import add_alias_options, default_out, report, run


def build_parser():
    parser = argparse.ArgumentParser(
        prog="reorder_gro.py",
        description=__doc__.splitlines()[0],
        epilog="\n".join(__doc__.splitlines()[1:]),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("gro", help="input .gro file")
    parser.add_argument("itp", nargs="+", help="one or more .itp topologies")
    parser.add_argument("-o", "--output",
                        help="output .gro (default: <name>_reordered.gro)")
    parser.add_argument("--resname", metavar="NAME",
                        help="with a single topology, the residue name the "
                             "coordinates use for it (same as --residue "
                             "NAME=<topology residue>)")
    parser.add_argument("--out-resname", action="append", default=[],
                        metavar="[ITP=]NAME",
                        help="write this residue name instead of the "
                             "topology's; repeatable")
    parser.add_argument("--keep-names", action="store_true",
                        help="re-order only, leaving residue and atom names "
                             "alone")
    parser.add_argument("--strict", action="store_true",
                        help="fail instead of passing through residues that "
                             "no topology covers")
    add_alias_options(parser)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    table = AliasTable()
    for path in args.alias_file:
        table.load(path)
    for spec in args.residue:
        table.add_residue_spec(spec)
    default_resname = None
    if len(args.itp) == 1:
        # With one topology an --alias needs no residue prefix.
        default_resname = ItpFile.read(args.itp[0]).residue_names[0]
    for spec in args.alias:
        table.add_atom_spec(spec, default_resname)

    if args.resname:
        if len(args.itp) != 1:
            parser.error("--resname needs a single topology; with several use "
                         "--residue GRO_NAME=ITP_NAME")
        table.add_residue(args.resname, default_resname)

    topologies = [
        ResidueTopology.from_itp(ItpFile.read(path), table)
        for path in args.itp
    ]
    out_resnames = parse_out_resnames(args.out_resname, topologies, parser)

    for top in topologies:
        report("%s: %d atoms per molecule, residue %s%s%s"
               % (top.itp.path, top.size, top.itp_resname,
                  "" if top.gro_resname == top.itp_resname
                  else " (%s in the coordinates)" % top.gro_resname,
                  "" if not top.aliases else ", aliases %s" % ", ".join(
                      "%s=%s" % pair for pair in sorted(top.aliases.items()))))

    frame = GroFrame.read(args.gro)
    counts = apply_topologies(
        frame, topologies,
        rename=not args.keep_names,
        strict=args.strict,
        out_resnames=out_resnames,
        report=report,
    )

    out_path = args.output or default_out(args.gro, ".gro")
    frame.write(out_path)
    report("%s: %d atoms; %s"
           % (args.gro, len(frame.atoms),
              ", ".join("%d %s" % (counts[name], name)
                        for name in sorted(counts))))
    report("wrote %s" % out_path)
    return 0


def parse_out_resnames(specs, topologies, parser):
    """`--out-resname NAME` or `--out-resname ITP=NAME` -> {itp resname: name}."""
    out = {}
    for spec in specs:
        if "=" in spec:
            itp_resname, name = spec.split("=", 1)
            out[itp_resname.strip()] = name.strip()
        elif len(topologies) == 1:
            out[topologies[0].itp_resname] = spec.strip()
        else:
            parser.error("--out-resname needs ITP_NAME=NAME when there is more "
                         "than one topology")
    known = {top.itp_resname for top in topologies}
    unknown = sorted(set(out) - known)
    if unknown:
        raise TopologyError(
            "--out-resname names no loaded topology: %s" % ", ".join(unknown)
        )
    return out


def cli():
    run(main)


if __name__ == "__main__":
    cli()
