"""Splitting CHARMM36 / SLipids phospholipids into AMBER Lipid21 residues.

Lipid21 splits each phospholipid into three residues -- tail, head, tail --
where CHARMM36 and SLipids keep the whole lipid in one.  The boundary follows
Lipid21, where the ester carbonyls belong to the *head*: the head is phosphate
+ glycerol + both C=O groups (``C21``/``O21``/``O22`` and ``C31``/``O31``/
``O32``), and each tail starts at the first methylene past the carbonyl
(``C22``/``C32``) and runs to the terminal methyl, hydrogens included.

Tail carbons are picked by name (``C2{n}``/``C3{n}``, n >= 2).  Each hydrogen
goes with the heavy atom it belongs to, found one of two ways:

  parents_by_bonds     -- from a topology's bonds (exact)
  parents_by_distance  -- the nearest heavy atom, minimum image (needs only
                          the coordinates)

so the result does not depend on the force field's hydrogen naming.
"""

from __future__ import annotations

import numpy as np

from .chem import heavy_parents, is_hydrogen
from .errors import TopologyError
from .ordering import order_from_names

# lipid resname -> (carbons per acyl chain, AMBER head residue, AMBER tail residue)
LIPIDS = {
    "DOPC": (18, "PC", "OL"),
    "DMPC": (14, "PC", "MY"),
    "DMPG": (14, "PGR", "MY"),
}

# Atoms per lipid in each AMBER residue (a lipid has one head and two tails).
AMBER_SIZE = {"PC": 38, "PGR": 31, "OL": 50, "MY": 40}


def lipid_size(resname):
    """Atoms in one whole lipid, which is the head plus both tails."""
    _, head, tail = LIPIDS[resname]
    return AMBER_SIZE[head] + 2 * AMBER_SIZE[tail]


def tail_carbon_names(n_carbons):
    return ["C%d%d" % (chain, c)
            for chain in (2, 3) for c in range(2, n_carbons + 1)]


def parents_by_distance(block, box):
    """For each atom of a molecule, the position of the heavy atom it goes with.

    A heavy atom goes with itself and a hydrogen (named H...) with the nearest
    heavy atom of the same molecule.  ``box`` is a pbc.Box, so a molecule split
    across the periodic boundary is handled.
    """
    positions = np.array([a.position for a in block])
    is_h = np.array([a.name.startswith("H") for a in block])
    heavy = np.flatnonzero(~is_h)
    hydrogens = np.flatnonzero(is_h)
    parents = np.arange(len(block))
    if len(hydrogens):
        d2 = box.distance2(positions[hydrogens, None, :],
                           positions[None, heavy, :])
        parents[hydrogens] = heavy[d2.argmin(axis=1)]
    return parents.tolist()


def parents_by_bonds(itp, order, where=""):
    """The same as parents_by_distance, from the bonds of a topology.

    ``order`` matches the molecule's coordinates to the topology, as
    ordering.order_from_names returns it: topology atom k is coordinate atom
    ``order.old_for_new[k - 1]``.
    """
    atoms = itp.atoms
    position = [old - 1 for old in order.old_for_new]
    parents = list(position)          # a heavy atom goes with itself
    of_hydrogen = heavy_parents(atoms, itp.bond_graph())
    for h, heavy in of_hydrogen.items():
        parents[position[h - 1]] = position[heavy - 1]
    orphans = [a.name for k, a in enumerate(atoms, start=1)
               if is_hydrogen(a) and k not in of_hydrogen]
    if orphans:
        raise TopologyError(
            "%s%s: hydrogen(s) bonded to no heavy atom: %s"
            % (where, itp.path or "topology", ", ".join(orphans))
        )
    return parents


def split_lipid(names, parents, n_carbons, where=""):
    """``(head, tail)``: 0-based positions within one lipid.

    ``parents`` is what parents_by_bonds or parents_by_distance returned.
    """
    index = {name: i for i, name in enumerate(names)}
    wanted = tail_carbon_names(n_carbons)
    missing = [name for name in wanted if name not in index]
    if missing:
        raise TopologyError("%slacks the tail carbons %s"
                            % (where, ", ".join(missing)))
    tail_heavy = {index[name] for name in wanted}
    head = [i for i, p in enumerate(parents) if p not in tail_heavy]
    tail = [i for i, p in enumerate(parents) if p in tail_heavy]
    return head, tail


def index_groups(frame, topologies=()):
    """The AMBER-residue index groups of every Lipid21 lipid in a frame.

    ``topologies`` are residues.ResidueTopology objects; a lipid with one is
    split by its bonds, any other by distance.  Returns ``(groups, counts)``:
    ``{AMBER residue name: sorted 1-based atom indices}`` and ``{lipid resname:
    lipids split}``.  Residues that are not one of LIPIDS are ignored.
    """
    by_resname = {}
    for top in topologies:
        if top.gro_resname not in LIPIDS:
            raise TopologyError(
                "%s describes residue %s, which is not one of %s"
                % (top.itp.path, top.gro_resname, "/".join(LIPIDS))
            )
        if top.gro_resname in by_resname:
            raise TopologyError("two topologies both describe %s"
                                % top.gro_resname)
        by_resname[top.gro_resname] = top

    box = None
    groups, counts = {}, {}
    for resname, (n_carbons, head_name, tail_name) in LIPIDS.items():
        top = by_resname.get(resname)
        solved = {}     # atom-name sequence -> parents, when there are bonds
        for start, block in frame.iter_residue_blocks(resname,
                                                      lipid_size(resname)):
            where = "%s: %s at atom %d: " % (frame._where(), resname, start + 1)
            names = [a.name for a in block]
            if top is None:
                if box is None:
                    box = frame.periodic_box
                parents = parents_by_distance(block, box)
            else:
                # Every molecule of a residue lists its atoms the same way in
                # practice, so match each distinct sequence only once.
                key = tuple(names)
                parents = solved.get(key)
                if parents is None:
                    order = order_from_names(top.names, names, top.aliases,
                                             where)
                    parents = solved[key] = parents_by_bonds(top.itp, order,
                                                             where)
            head, tail = split_lipid(names, parents, n_carbons, where)
            if (len(head) != AMBER_SIZE[head_name]
                    or len(tail) != 2 * AMBER_SIZE[tail_name]):
                raise TopologyError(
                    "%ssplit %d/%d, expected %d/%d -- unexpected topology"
                    % (where, len(head), len(tail), AMBER_SIZE[head_name],
                       2 * AMBER_SIZE[tail_name])
                )
            groups.setdefault(head_name, []).extend(start + 1 + i for i in head)
            groups.setdefault(tail_name, []).extend(start + 1 + i for i in tail)
            counts[resname] = counts.get(resname, 0) + 1

    if not counts:
        raise TopologyError("%s: no %s residues"
                            % (frame._where(), "/".join(LIPIDS)))
    return {name: sorted(ix) for name, ix in groups.items()}, counts
