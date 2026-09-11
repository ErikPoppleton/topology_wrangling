"""Tools for reading, re-ordering and checking GROMACS topologies.

    from topology_wrangling import ItpFile, GroFrame, chem

    top   = ItpFile.read("DOPC.itp")
    frame = GroFrame.read("system.gro")
    for start, block in frame.iter_residue_blocks(top.residue_name(), top.natoms):
        xyz = [atom.position for atom in block]

Everything here raises TopologyError on bad input and never exits, so it can be
driven from a script as well as from the command line.
"""

from .errors import TopologyError
from .gro import GroAtom, GroFrame
from .itp import Atom, ItpFile, DIRECTIVES
from .mapping import AtomOrder
from .ordering import order_from_names, update_group_order
from .residues import AliasTable, ResidueTopology, apply_topologies

__all__ = [
    "TopologyError",
    "GroAtom", "GroFrame",
    "Atom", "ItpFile", "DIRECTIVES",
    "AtomOrder",
    "update_group_order", "order_from_names",
    "AliasTable", "ResidueTopology", "apply_topologies",
]
