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

### `amber_style_index.py` — AMBER Lipid21 head/tail groups for CHARMM/SLipids lipids

Lipid21 splits each phospholipid into three residues (tail, head, tail) where
CHARMM36 and SLipids keep one.  This writes a `.ndx` whose groups hold exactly
the atoms of the AMBER residues (`PC`/`PGR` heads, `OL`/`MY` tails) for every
DOPC, DMPC and DMPG in a `.gro`.  Hydrogens go with the heavy atom they are
bonded to: read from the lipid's `.itp` when one is given, otherwise the
nearest heavy atom by minimum-image distance (rectangular or triclinic box).

```sh
./amber_style_index.py examples/charmm_DOPC.gro                 # -> examples/charmm_DOPC_amber.ndx
./amber_style_index.py examples/charmm_DOPC.gro examples/DOPC_reordered.itp -o lipids.ndx
```

### `draw_itp.py` — see a topology as a labelled 2D graph

There are plenty of 3D molecular viewers; this makes the flat picture a
topology file calls for. Atoms become circles carrying their unique name from
the `.itp`, bonds become edges, and the layout follows the conventions of a
structural formula rather than a graph drawing:

* a chain zigzags, consecutive bonds meeting at a 120° interior angle,
* a branching atom fans its bonds evenly (three 120° apart, four 90° apart),
* a ring is a regular polygon, and a ring fused to one already placed is
  reflected across the bond they share,
* hydrogens hang off their heavy atom on short spokes, in whatever directions
  the skeleton leaves free, instead of joining the chain.

Every bond comes out the same length, so a long lipid reads as a long lipid.

Which way a bond turns is decided as the walk passes through it, but whether
that choice clashes is only known once the rest of the molecule is down — a
lipid head group folds back over the glycerol it hangs off, and nothing at the
moment of placing it could have said so. So after the first pass the layout
finds the atoms that are too close, looks at the bonds on the path between
them, and reflects the branch beyond one of those bonds across its own axis.
A reflection is a rigid motion, so the branch swings to the other side with
every bond length and every angle exactly as it was.

```sh
./draw_itp.py examples/DOPC_reordered.itp -o dopc.svg
./draw_itp.py examples/DOPC_reordered.itp -o dopc.png --hide-hydrogens
./draw_itp.py examples/DOPC_reordered.itp -o charges.svg --color-by charge --cmap bwr
```

Atoms are filled with the [Jmol element colours](https://jmol.sourceforge.net/jscolors/)
by default — the same scheme a 3D viewer uses — and each label is drawn black or
white, whichever gives the better contrast against its fill. `--label type`
writes atom types instead of names; `--label nr` writes indices.

### Force field parameters as a colour overlay

```sh
./draw_itp.py DOPC.itp --color-by charge
./draw_itp.py DOPC.itp --color-by epsilon --ff examples/slipids_ffnonbonded.itp
```

| overlay | source | colour map |
| --- | --- | --- |
| `charge` | the molecule's own `[ atoms ]` records | `coolwarm`, centred on neutral |
| `epsilon` | `[ atomtypes ]` of the file(s) given with `--ff` | `viridis`, over the molecule's range |

A signed quantity gets a diverging map centred on zero, so the sign reads
straight off the colour; one with no zero crossing gets a sequential map over
its own range. `--cmap` overrides either and a colour bar is drawn with the units.

`--range MIN MAX` fixes the scale, which is what puts several molecules on one
scale — take a first pass over every file to find the global range, then render
them all with it:

```sh
python3 - <<'EOF'
from topology_wrangling import ItpFile
# ... min/max of float(a.charge) and of types[a.type].epsilon over every file
EOF
./draw_itp.py mol.itp --color-by epsilon --ff mol_ff.itp --range 0.0293 2.4476
```

### Several molecules in one figure

Draw each panel bare and give them one shared colour bar:

```sh
for m in charmm_DOPC amber_DOPC slipids_DOPC; do
    ./draw_itp.py $m.itp --color-by epsilon --ff ${m}_ff.itp \
        --range 0.0293 2.4476 --bare -o ${m}.png --dpi 150
done
./draw_colorbar.py --color-by epsilon --range 0.0293 2.4476 -o bar.png --dpi 150
```

`--bare` drops the title and the colour bar, leaving just the structure
(`--no-title` and `--no-colorbar` do one each). `draw_colorbar.py` writes a bar
with no molecule attached — same caption, units and colour map, `--orientation
vertical` for a bar down the side of a panel grid. Its `--range` is required,
since there is no molecule to take one from and the whole point is that it
matches the panels.

### Capping a scale that one atom dominates

A single outlier flattens everything else: phosphorus at 2.45 kJ/mol leaves
every carbon in a lipid within the bottom tenth of the epsilon scale and so
effectively one colour. Topping the scale lower spreads the rest out, and
`--extend max` labels that end `>value` so an atom drawn in the end colour is
not read as sitting at it:

```sh
./draw_itp.py DOPC.itp --color-by epsilon --ff DOPC_ff.itp --range 0.0293 0.9205
./draw_colorbar.py --color-by epsilon --range 0.0293 0.9205 --extend max -o bar.png
```

Atoms past either end are drawn in that end's colour. `draw_itp.py` works out
`--extend` from the atoms it drew unless told otherwise; the standalone bar has
no molecule to work it out from, so it must be given.

**Partial charge is read from the molecule, not the force field.** A charge set
in a molecule's `[ atoms ]` record overrides the force field's default for that
atom type, and the two routinely differ — Slipids gives the type `NTL` a charge
of −0.60, while the nitrogen of DOPC carries +0.20. Epsilon has no per-atom
form, so it is looked up by atom type; an atom whose type is not in the files
given is an error naming the type, not a silently blank atom.

`[ settles ]` and `[ constraints ]` count as bonds for the drawing, so rigid
water comes out connected rather than as three loose dots.

After unfolding, the layout keeps flipping branches for as long as the drawing
gets smaller, scoring clashes first, then crossing bonds, then the area of page
used per atom. That last term is what packs a lipid's two tails roughly
parallel instead of splaying them: parallel chains fit in a fraction of the box
a V needs. There is deliberately no preference for a square page — a straight
alkane is thin, and thin is how one is drawn.

`--branch` sets where a branch aims: `compact` (the default) lets it carry on
forwards, which packs tails parallel; `spread` sends the longest branch away
from everything drawn so far, putting the tails in a line; `auto` measures both
and keeps the tighter. `--layout force` brings back the old spring embedding,
worth a try for a molecule the skeletal rules do badly on.

The canvas is sized from the atoms rather than the other way round — a long
lipid gets more room than a water, so labels stay readable at any molecule
size. `--node-radius` sets the scale and `--max-size` caps the result.

Any of the tools can also be reached as `python -m topology_wrangling <tool>`
(`update-group-ordering`, `reorder-gro`, `draw-itp`, `draw-colorbar`).

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
| `pbc` | `Box`, rectangular and triclinic, and the minimum-image convention |
| `ndx` | writing GROMACS index files |
| `lipid21` | splitting CHARMM/SLipids lipids into AMBER Lipid21 head and tail groups |
| `colors` | the Jmol element table, contrast-picked label colours, value maps |
| `skeleton` | the structural-formula layout: chain angles, branches, rings |
| `layout` | the force-directed alternative, and overlap relaxation |
| `draw` | topology to a Scene, and the SVG and PNG backends |
| `cli` | the command line front-ends |

Two properties the design depends on:

* an unmodified file written back out is **byte-identical** to the input --
  directives the library does not know about survive a round trip, and a
  `GroAtom` keeps the text of its coordinate columns verbatim, so re-ordering
  and renaming can never perturb a number;
* library code raises `TopologyError` and never exits, so every module can be
  driven from a script.

## Dependencies

**numpy** is required: distances between atoms (periodic boxes, the
`amber_style_index.py` split) and the force-directed drawing layout use it.
**matplotlib** is needed only for `.png` output and `--color-by charge`; an
`.svg` with element colours does without it.

## Tests

```sh
python3 -m unittest discover -s tests
```

The suite runs both tools over the example system and checks the properties
that matter: the re-ordered topology still describes the same interactions,
every hydrogen follows a heavy atom it is bonded to, and the re-ordered
coordinate file conserves every coordinate of the input.
