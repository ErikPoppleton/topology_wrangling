"""Element and connectivity helpers shared by every tool.

Nothing here knows about file formats: it works on the Atom records produced
by topology_wrangling.itp, so a re-ordering script, a visualiser and anything
added later all agree on what counts as a hydrogen and which element an atom
is.
"""

from __future__ import annotations

from collections import defaultdict

from .errors import TopologyError

# Anything at or below this mass is treated as a hydrogen.  The bound is loose
# enough to cover deuterium and the heavy-hydrogen (mass-repartitioned) atoms
# some force fields use.
HYDROGEN_MASS_MAX = 3.5

# Standard atomic masses for the elements that show up in biomolecular force
# fields, used to name an element from the mass column.
_MASSES = {
    "H": 1.008, "D": 2.014, "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990,
    "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Fe": 55.845,
    "Cu": 63.546, "Zn": 65.38, "Br": 79.904, "I": 126.904, "Cs": 132.905,
}

# Two-letter symbols we are willing to read off an atom type or name.  One
# letter symbols are the fallback, so a CHARMM type like "CTL2" gives C and
# "CLA" gives Cl.
_TWO_LETTER = ("Cl", "Br", "Na", "Mg", "Ca", "Zn", "Fe", "Cu", "Li", "Si")


def is_hydrogen(atom):
    """Whether an Atom record is a hydrogen.

    The mass column is the reliable signal.  When it is missing we fall back to
    the atom type and then the atom name, which start with an H for hydrogens
    in the CHARMM, Slipids and AMBER naming schemes alike.
    """
    if atom.mass is not None:
        return atom.mass <= HYDROGEN_MASS_MAX
    for value in (atom.type, atom.name):
        if value:
            return value[0].upper() == "H"
    return False


def element_of(atom):
    """Best guess at the element symbol of an Atom record.

    Mass first, because it is unambiguous; then the atom type and name, where
    a leading two-letter symbol wins over a one-letter one so CLA is chlorine
    rather than carbon.  Returns "X" when nothing matches.
    """
    if atom.mass is not None:
        best, best_diff = None, None
        for symbol, mass in _MASSES.items():
            diff = abs(mass - atom.mass)
            if best_diff is None or diff < best_diff:
                best, best_diff = symbol, diff
        # 0.5 u is wide enough for rounded force-field masses and narrow
        # enough that it never bridges two real elements.
        if best_diff is not None and best_diff <= 0.5:
            return best
    for value in (atom.type, atom.name):
        symbol = _element_from_label(value)
        if symbol:
            return symbol
    return "X"


def _element_from_label(label):
    if not label:
        return None
    text = label.strip().lstrip("0123456789")
    if not text:
        return None
    head = text[:2].capitalize()
    if head in _TWO_LETTER:
        return head
    head = text[:1].upper()
    return head if head in _MASSES else None


def bond_graph(bonds, natoms):
    """Adjacency lists (1-based) from an iterable of (i, j) atom index pairs."""
    graph = defaultdict(list)
    for ai, aj in bonds:
        for a in (ai, aj):
            if not 1 <= a <= natoms:
                raise TopologyError(
                    "bond refers to atom %d, which is outside [ atoms ] (1-%d)"
                    % (a, natoms)
                )
        graph[ai].append(aj)
        graph[aj].append(ai)
    return graph
