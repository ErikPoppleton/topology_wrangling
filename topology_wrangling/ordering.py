"""Deciding a new atom order.

Two sources, both returning an AtomOrder:

  update_group_order  -- put every hydrogen straight after its heavy atom, so
                         the topology works with GROMACS update groups
  order_from_names    -- match the atoms of a coordinate file to a topology by
                         atom name
"""

from __future__ import annotations

from collections import defaultdict

from .chem import heavy_parents, is_hydrogen
from .errors import TopologyError
from .mapping import AtomOrder


def update_group_order(atoms, graph, resname=None):
    """Order the atoms so each heavy atom is followed by its hydrogens.

    Heavy atoms keep their original relative order, as do the hydrogens of a
    given heavy atom.  A hydrogen bonded to several heavy atoms is attached to
    the first of them; a hydrogen with no heavy neighbour keeps its own slot.
    """
    hydrogen = [is_hydrogen(a) for a in atoms]
    natoms = len(atoms)

    parents = heavy_parents(atoms, graph)
    attached = defaultdict(list)   # heavy atom (old index) -> its hydrogens
    for idx, parent in parents.items():
        attached[parent].append(idx)
    placed = set(parents)

    old_for_new = []
    for idx in range(1, natoms + 1):
        if hydrogen[idx - 1]:
            if idx not in placed:   # orphan hydrogen: keep its original slot
                old_for_new.append(idx)
            continue
        old_for_new.append(idx)
        old_for_new.extend(sorted(attached[idx]))

    order = AtomOrder(old_for_new, [a.name for a in atoms], resname)
    return order, hydrogen


def order_from_names(itp_names, gro_names, aliases=None, where=""):
    """Order a molecule's coordinate atoms to match a topology's atom order.

    ``aliases`` maps a topology atom name to the name the coordinate file uses
    for the same atom.  Matching is by name, which is unique within a molecule,
    so the mapping is unambiguous.
    """
    aliases = aliases or {}
    unknown = sorted(set(aliases) - set(itp_names))
    if unknown:
        raise TopologyError(
            "%salias names not present in the topology: %s"
            % (where, ", ".join(unknown))
        )
    if len(set(itp_names)) != len(itp_names):
        raise TopologyError(
            "%sthe topology repeats an atom name, so the order cannot be "
            "matched by name" % where
        )

    wanted = [aliases.get(name, name) for name in itp_names]
    if len(set(wanted)) != len(wanted):
        raise TopologyError("%sthe alias mapping collides with an existing "
                            "atom name" % where)

    if len(gro_names) != len(wanted):
        raise TopologyError(
            "%sthe molecule has %d atoms but the topology has %d"
            % (where, len(gro_names), len(wanted))
        )
    if set(gro_names) != set(wanted):
        raise TopologyError(
            "%sthe molecule does not match the topology.\n"
            "  only in the topology: %s\n"
            "  only in the coordinates: %s\n"
            "Map differing names with an alias."
            % (where,
               ", ".join(sorted(set(wanted) - set(gro_names))) or "-",
               ", ".join(sorted(set(gro_names) - set(wanted))) or "-")
        )

    position = {name: i for i, name in enumerate(gro_names)}
    old_for_new = [position[name] + 1 for name in wanted]
    return AtomOrder(old_for_new, list(gro_names))
