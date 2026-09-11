# topology_wrangling

Tools for reading, re-ordering and renaming GROMACS topologies and coordinate
files, built on a small library so that new tools do not have to re-implement
the parsing.

## The tools

### `update_group_ordering.py` — make a topology work with update groups

GROMACS update groups require that the hydrogens of a heavy atom be listed
immediately after it rather than in blocks.  This re-orders `[ atoms ]` and
rewrites every directive that refers to atom indices with the new numbering.
Only the numbers change: functional forms, parameters and the order of the
atoms *inside* a bonded interaction are untouched, so impropers keep their
meaning.

```sh
./update_group_ordering.py examples/update_group_incompatible.itp
# -> examples/update_group_incompatible_reordered.itp
```

### `reorder_gro.py` — make coordinates match a set of topologies

For every residue a topology describes, the atoms of each molecule are
re-ordered into the topology's order, the residue and atom names are rewritten
to the topology's, and the atom-index column is renumbered.  Matching is by
atom name, so the input order does not matter.  Residues no topology covers
are left exactly where they are.

```sh
./reorder_gro.py examples/charmm_DOPC.gro examples/DOPC_reordered.itp \
    examples/TP3.itp examples/Na+.itp examples/Cl-.itp examples/MG.itp \
    examples/CA.itp --alias-file examples/charmm.map -o system.gro
```

Either tool can also be reached as `python -m topology_wrangling update-group-ordering` or `... reorder-gro`.

## The worked example

`examples/` holds a CHARMM-GUI membrane system and the topologies that describe
it.  The DOPC topology is from Slipids and lists its hydrogens in blocks, and
the coordinate file uses CHARMM names (`TIP3`/`OH2`) where the topologies use
the converter's (`TP3`/`O`).  The two steps together fix both:

```sh
./update_group_ordering.py examples/update_group_incompatible.itp -o DOPC.itp
./reorder_gro.py examples/charmm_DOPC.gro DOPC.itp examples/TP3.itp \
    examples/Na+.itp examples/Cl-.itp examples/MG.itp examples/CA.itp \
    --alias-file examples/charmm.map -o system.gro
```

`examples/charmm.map` records the name differences; without it the tool reports
the unmatched names on each side, so the aliases a file needs can be read
straight off the error.  Name differences can also be given inline with
`--residue TIP3=TP3` and `--alias TP3:O=OH2`.

## The library

```python
from topology_wrangling import ItpFile, GroFrame, chem

top   = ItpFile.read("DOPC.itp")
frame = GroFrame.read("system.gro")

bonds = top.bond_graph()                      # {index: [neighbour, ...]}
elements = [chem.element_of(a) for a in top.atoms]

for start, molecule in frame.iter_residue_blocks(top.residue_name(), top.natoms):
    xyz = [atom.position for atom in molecule]
```

| module | holds |
| --- | --- |
| `itp` | `ItpFile`, `Atom`, the directive table; parses, edits and writes `.itp` |
| `gro` | `GroFrame`, `GroAtom`; fixed-column parsing, coordinates on demand |
| `mapping` | `AtomOrder`, the permutation type, and its map-file form |
| `ordering` | building an order from a bond graph or by atom-name matching |
| `residues` | `AliasTable`, matching a frame's residues to topologies |
| `chem` | hydrogen and element identification, bond graphs |
| `cli` | the command line front-ends |

Two properties the design depends on:

* an unmodified file written back out is **byte-identical** to the input --
  directives the library does not know about survive a round trip, and a
  `GroAtom` keeps the text of its coordinate columns verbatim, so re-ordering
  and renaming can never perturb a number;
* library code raises `TopologyError` and never exits, so every module can be
  driven from a script.

## Tests

```sh
python3 -m unittest discover -s tests
```

The suite runs both tools over the example system and checks the properties
that matter: the re-ordered topology still describes the same interactions,
every hydrogen follows a heavy atom it is bonded to, and the re-ordered
coordinate file conserves every coordinate of the input.
