"""Write a colour bar on its own, with no molecule attached.

Several molecules going into one figure should share one colour bar rather
than each carrying its own.  Draw the panels with `draw_itp.py --bare` and the
bar with this, on the same --range.

The range is required: there is no molecule here to take it from, and the
whole point is that it matches the panels.

Examples:
    draw_colorbar.py --color-by charge --range -1.58 1.58 -o charge_bar.png
    draw_colorbar.py --color-by epsilon --range 0.029 2.448 -o eps_bar.svg \
        --orientation vertical
"""

from __future__ import annotations

import argparse

from ..draw import PROPERTIES, colorbar_scene
from .common import report, run


def build_parser():
    parser = argparse.ArgumentParser(
        prog="draw_colorbar.py",
        description=__doc__.splitlines()[0],
        epilog="\n".join(__doc__.splitlines()[1:]),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--color-by", choices=sorted(PROPERTIES),
                        required=True,
                        help="which quantity the bar is for; this sets its "
                             "caption, units and default colour map")
    parser.add_argument("--range", dest="value_range", nargs=2, type=float,
                        required=True, metavar=("MIN", "MAX"),
                        help="the scale, which must match the figures it "
                             "accompanies")
    parser.add_argument("-o", "--output", required=True,
                        help="output .svg or .png")
    parser.add_argument("--cmap", metavar="NAME",
                        help="matplotlib colour map, overriding the one that "
                             "suits the property")
    parser.add_argument("--extend", choices=("neither", "min", "max", "both"),
                        default="neither",
                        help="mark an end the scale does not reach: --extend "
                             "max labels it \">value\", for a scale deliberately "
                             "topped below the largest atom (default: neither)")
    parser.add_argument("--orientation", choices=("horizontal", "vertical"),
                        default="horizontal",
                        help="default: horizontal")
    parser.add_argument("--length", type=float, default=360.0, metavar="PX",
                        help="length of the bar itself (default: 360)")
    parser.add_argument("--thickness", type=float, default=18.0, metavar="PX",
                        help="thickness of the bar (default: 18)")
    parser.add_argument("--font", type=float, default=14.0, metavar="PX",
                        help="label size (default: 14)")
    parser.add_argument("--dpi", type=int, default=200,
                        help="resolution for .png output (default: 200)")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    scene = colorbar_scene(
        args.color_by,
        cmap=args.cmap,
        value_range=tuple(args.value_range),
        length=args.length,
        thickness=args.thickness,
        orientation=args.orientation,
        font=args.font,
        extend=args.extend,
    )
    scene.save(args.output, dpi=args.dpi)
    report("%s colour bar, %g to %g" % (args.color_by, *args.value_range))
    report("wrote %s (%.0f x %.0f)" % (args.output, scene.width, scene.height))
    return 0


def cli():
    run(main)


if __name__ == "__main__":
    cli()
