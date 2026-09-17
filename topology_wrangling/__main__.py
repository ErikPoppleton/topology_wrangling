"""`python -m topology_wrangling <tool> ...` -- one entry point for every tool."""

from __future__ import annotations

import sys

TOOLS = {
    "update-group-ordering": "topology_wrangling.cli.update_group_ordering",
    "draw-colorbar": "topology_wrangling.cli.draw_colorbar",
    "draw-itp": "topology_wrangling.cli.draw_itp",
    "reorder-gro": "topology_wrangling.cli.reorder_gro",
}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help") or argv[0] not in TOOLS:
        sys.stderr.write(
            "usage: python -m topology_wrangling <tool> [options]\n\ntools:\n"
            + "".join("  %s\n" % name for name in sorted(TOOLS))
        )
        return 0 if argv[:1] in ([], ["-h"], ["--help"]) else 2
    module = __import__(TOOLS[argv[0]], fromlist=["cli"])
    sys.argv = [argv[0]] + argv[1:]
    module.cli()


if __name__ == "__main__":
    sys.exit(main())
