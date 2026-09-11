#!/usr/bin/env python3
"""Entry point for `update_group_ordering.py`; the code lives in the library.

See topology_wrangling/cli/update_group_ordering.py.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from topology_wrangling.cli.update_group_ordering import cli

if __name__ == "__main__":
    cli()
