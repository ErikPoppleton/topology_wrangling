#!/usr/bin/env python3
"""Entry point for `draw_colorbar.py`; the code lives in the library.

See topology_wrangling/cli/draw_colorbar.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from topology_wrangling.cli.draw_colorbar import cli

if __name__ == "__main__":
    cli()
