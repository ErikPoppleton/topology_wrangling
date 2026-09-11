#!/usr/bin/env python3
"""Entry point for `reorder_gro.py`; the code lives in the library.

See topology_wrangling/cli/reorder_gro.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from topology_wrangling.cli.reorder_gro import cli

if __name__ == "__main__":
    cli()
