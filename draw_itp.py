#!/usr/bin/env python3
"""Entry point for `draw_itp.py`; the code lives in the library.

See topology_wrangling/cli/draw_itp.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from topology_wrangling.cli.draw_itp import cli

if __name__ == "__main__":
    cli()
