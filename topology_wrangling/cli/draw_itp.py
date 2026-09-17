"""Draw a GROMACS topology as a labelled 2D graph.

There are plenty of 3D molecular viewers; this makes the flat picture a
topology file calls for.  Atoms become circles carrying their unique name from
the .itp, bonds become edges, and a force-directed layout untangles the graph
into something that reads like a structural formula.

Atoms are filled with the Jmol element colours by default -- the same scheme a
3D viewer uses -- and each label is drawn black or white, whichever contrasts
better with the fill.  They can be filled from a force field parameter instead:

    --color-by charge     the molecule's own per-atom partial charge, on a
                          blue-red diverging map centred on neutral
    --color-by epsilon    the Lennard-Jones well depth of each atom's type,
                          on viridis; needs --ff to say where the types are

Partial charge is read from the molecule's [ atoms ] records, not from the
force field: a charge set there overrides the force field's default for the
type, and the two routinely differ.  Epsilon has no per-atom form, so it is
looked up by atom type in the [ atomtypes ] directive of the file(s) given
with --ff.

Output is .svg (which needs nothing but the standard library) or .png (which
needs matplotlib).

Examples:
    draw_itp.py DOPC.itp                            # -> DOPC.svg
    draw_itp.py DOPC.itp -o dopc.png --hide-hydrogens
    draw_itp.py DOPC.itp --label type
    draw_itp.py DOPC.itp --color-by charge
    draw_itp.py DOPC.itp --color-by epsilon --ff ffnonbonded.itp
"""

from __future__ import annotations

import argparse
import os

from ..draw import draw_topology
from ..itp import ItpFile
from .common import report, run


def build_parser():
    parser = argparse.ArgumentParser(
        prog="draw_itp.py",
        description=__doc__.splitlines()[0],
        epilog="\n".join(__doc__.splitlines()[1:]),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("itp", help="topology to draw")
    parser.add_argument("-o", "--output",
                        help="output .svg or .png (default: <name>.svg)")
    parser.add_argument("--label", choices=("name", "type", "nr", "none"),
                        default="name",
                        help="what to write in each atom: its unique name from "
                             "the .itp (default), its atom type, its number, "
                             "or nothing")
    parser.add_argument("--color-by", choices=("element", "charge", "epsilon"),
                        default="element",
                        help="fill atoms with the Jmol element colours "
                             "(default), with the molecule's per-atom partial "
                             "charge, or with the LJ epsilon of each atom's "
                             "type (which needs --ff)")
    parser.add_argument("--ff", action="append", default=[], metavar="FILE",
                        help="force field file holding [ atomtypes ]; "
                             "repeatable, for a force field split over several")
    parser.add_argument("--cmap", metavar="NAME",
                        help="matplotlib colour map, overriding the one that "
                             "suits the property (coolwarm for charge, "
                             "viridis for epsilon)")
    parser.add_argument("--extend", choices=("auto", "neither", "min", "max",
                                            "both"),
                        default="auto",
                        help="mark the ends of the colour bar that the scale "
                             "does not reach, so an atom drawn in the end "
                             "colour is not read as sitting at that value "
                             "(default: worked out from the atoms drawn)")
    parser.add_argument("--range", dest="value_range", nargs=2, type=float,
                        metavar=("MIN", "MAX"),
                        help="fix the colour scale (default: symmetric about "
                             "zero for charge, the molecule's own span for "
                             "epsilon)")
    parser.add_argument("--layout", choices=("skeletal", "force"),
                        default="skeletal",
                        help="skeletal draws it like a structural formula, "
                             "with 120 degree chains and polygon rings "
                             "(default); force is the old spring model")
    parser.add_argument("--branch", choices=("compact", "spread", "auto"),
                        default="compact",
                        help="how a skeletal layout aims branches: compact "
                             "packs a lipid's tails roughly parallel "
                             "(default), spread puts them in a line, auto "
                             "measures both and keeps the tighter")
    parser.add_argument("--hide-hydrogens", action="store_true",
                        help="draw only the heavy atoms, which in an all-atom "
                             "topology is most of the picture")
    parser.add_argument("--node-radius", type=float, default=24.0,
                        metavar="PX",
                        help="drawn radius of a heavy atom; the canvas is "
                             "sized to suit (default: 24)")
    parser.add_argument("--max-size", type=float, default=5000.0, metavar="PX",
                        help="longest edge of the canvas, past which the "
                             "drawing is scaled down (default: 5000)")
    parser.add_argument("--iterations", type=int, default=400, metavar="N",
                        help="layout iterations (default: 400)")
    parser.add_argument("--seed", type=int, default=0,
                        help="layout seed; the same seed always gives the same "
                             "picture (default: 0)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="resolution for .png output (default: 200)")
    parser.add_argument("--no-title", action="store_true",
                        help="leave off the molecule name and atom count")
    parser.add_argument("--no-colorbar", action="store_true",
                        help="leave off the colour bar")
    parser.add_argument("--bare", action="store_true",
                        help="just the structure: no title, no colour bar. "
                             "For panels of one figure that share a single "
                             "bar, which draw_colorbar.py writes on its own")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    itp = ItpFile.read(args.itp)

    atomtypes = {}
    for path in args.ff:
        atomtypes.update(ItpFile.read(path).atomtypes)
    if atomtypes:
        report("%s: %d atom types" % (", ".join(args.ff), len(atomtypes)))

    scene = draw_topology(
        itp,
        label=args.label,
        color_by=args.color_by,
        cmap=args.cmap,
        value_range=tuple(args.value_range) if args.value_range
        else (None, None),
        atomtypes=atomtypes,
        extend=args.extend,
        hide_hydrogens=args.hide_hydrogens,
        layout=args.layout,
        branch=args.branch,
        iterations=args.iterations,
        seed=args.seed,
        node_radius=args.node_radius,
        max_size=args.max_size,
        title=not (args.no_title or args.bare),
        colorbar=not (args.no_colorbar or args.bare),
    )

    out_path = args.output or os.path.splitext(args.itp)[0] + ".svg"
    scene.save(out_path, dpi=args.dpi)

    drawn = len(scene.circles)
    report("%s: %s, %d of %d atoms drawn, coloured by %s"
           % (args.itp, itp.molecule_name or itp.residue_names[0], drawn,
              itp.natoms, args.color_by))
    report("wrote %s (%.0f x %.0f)" % (out_path, scene.width, scene.height))
    return 0


def cli():
    run(main)


if __name__ == "__main__":
    cli()
