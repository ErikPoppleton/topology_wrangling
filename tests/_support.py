"""Shared paths and imports for the tests."""

import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def example(name):
    return os.path.join(EXAMPLES, name)


# The topologies that describe the example system, and the residue name the
# example .gro uses for each.
SYSTEM = [
    ("DOPC_reordered.itp", "DOPC", 138, 128),
    ("TP3.itp", "TIP3", 3, 5909),
    ("Na+.itp", "SOD", 1, 14),
    ("Cl-.itp", "CLA", 1, 26),
    ("MG.itp", "MG", 1, 4),
    ("CA.itp", "CAL", 1, 2),
]


def read(path):
    with open(path) as handle:
        return handle.read()
