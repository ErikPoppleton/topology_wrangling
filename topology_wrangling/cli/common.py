"""Plumbing shared by the command line front-ends.

Library code raises TopologyError; this is the only place that prints and the
only place that exits, so a new tool gets consistent diagnostics for free.
"""

from __future__ import annotations

import os
import sys

from ..errors import TopologyError
from ..residues import AliasTable


def report(message):
    """Progress and summary lines, on stderr so stdout stays usable in a pipe."""
    sys.stderr.write(message + "\n")


def default_out(path, suffix=None, tag="_reordered"):
    """`md.gro` -> `md_reordered.gro`."""
    root, ext = os.path.splitext(path)
    return root + tag + (suffix or ext or "")


def add_alias_options(parser):
    """The residue and atom naming options, shared by the coordinate tools."""
    parser.add_argument(
        "--residue",
        action="append",
        default=[],
        metavar="GRO=ITP",
        help="the coordinates call this topology's residue something else "
        "(e.g. --residue TIP3=TP3); repeatable",
    )
    parser.add_argument(
        "--alias",
        action="append",
        default=[],
        metavar="[RES:]ITP=GRO",
        help="the coordinates call this atom something else "
        "(e.g. --alias TP3:O=OH2); repeatable",
    )
    parser.add_argument(
        "--alias-file",
        action="append",
        default=[],
        metavar="FILE",
        help="read 'residue' and 'atom' alias lines from a file; repeatable",
    )
    return parser


def alias_table(args, itps):
    """The AliasTable that add_alias_options' arguments describe.

    With a single topology an --alias needs no residue prefix.
    """
    table = AliasTable()
    for path in args.alias_file:
        table.load(path)
    for spec in args.residue:
        table.add_residue_spec(spec)
    default_resname = itps[0].residue_names[0] if len(itps) == 1 else None
    for spec in args.alias:
        table.add_atom_spec(spec, default_resname)
    return table


def run(main):
    """Call a main() that may raise TopologyError, and exit accordingly."""
    try:
        sys.exit(main() or 0)
    except TopologyError as err:
        report("error: %s" % err)
        sys.exit(1)
    except BrokenPipeError:
        sys.exit(0)
