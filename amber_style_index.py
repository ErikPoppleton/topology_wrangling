#!/usr/bin/env python3
"""Entry point for `amber_style_index.py`; the code lives in the library.

See topology_wrangling/cli/amber_style_index.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from topology_wrangling.cli.amber_style_index import cli

if __name__ == "__main__":
    cli()
