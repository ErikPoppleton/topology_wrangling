# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Python 3 with numpy, no build step and no installer. numpy may be used in any
module, and **distance calculations on coordinates use it** — vectorised over
whole arrays through `pbc.Box`, not pair by pair in Python. matplotlib is only
for `.png` output and charge colour maps; keep it out of everything else.

```sh
# all tests
python3 -m unittest discover -s tests

# one test, by name pattern, from the repo root
python3 -m unittest discover -s tests -k test_permutation_is_complete

# one module or class (must run from inside tests/, which puts _support.py on sys.path)
cd tests && python3 -m unittest test_itp.UpdateGroupOrderTest

# the three tools, on the example system
./update_group_ordering.py examples/update_group_incompatible.itp -o DOPC.itp
./reorder_gro.py examples/charmm_DOPC.gro DOPC.itp examples/TP3.itp examples/Na+.itp \
    examples/Cl-.itp examples/MG.itp examples/CA.itp --alias-file examples/charmm.map -o system.gro
```

`python3 -m unittest tests.test_itp` does **not** work — the tests import `_support`
by bare name, so the `tests/` directory itself has to be on `sys.path`.

## Architecture

The two re-ordering tools are the same pipeline: **read a topology → derive an
`AtomOrder` → apply it to something**. They differ only in where the order comes
from and what it is applied to. Keep new re-ordering tools on that spine.

The drawing tool is the other consumer of the same front half: **read a
topology → lay its graph out → build a `Scene` → write it**.

```
                       ordering.update_group_order   (from the bond graph)
itp.ItpFile  ───────>                                ───> mapping.AtomOrder ──┬─> ItpFile.permute()
gro.GroFrame ───────>  ordering.order_from_names     (by atom name)           └─> residues.apply_topologies()
```

`mapping.AtomOrder` is the currency between the halves, and its map-file form is
what lets `update_group_ordering --write-map` hand a permutation to a later
`--read-map` run. `residues.AliasTable` carries the name disagreements (a
CHARMM `.gro` says `TIP3`/`OH2` where a converter `.itp` says `TP3`/`O`) that
name matching needs.

```
                              skeleton.skeletal   (structural formula, default)
itp.ItpFile.connectivity() ─>                     ─> draw.Scene ─┬─> to_svg() (stdlib)
                              layout.force_directed (spring model)└─> to_png() (matplotlib)
```

`draw.Scene` is a bag of plain primitives, deliberately backend-neutral: the
SVG and PNG writers both consume it, which is what stops the two outputs from
drifting apart. Add a primitive there, not in one backend.

The `cli/` front-ends are thin drivers; `cli/common.py` holds the shared
argparse groups, `report()` and the `run()` wrapper.

### Three invariants the design rests on

1. **Byte-identical round-trip.** `ItpFile.read(p).text() == open(p).read()` and
   the same for `GroFrame`. `ItpFile` keeps every raw line tagged with its
   directive, so directives the library does not model survive untouched; a
   `GroAtom` keeps everything from column 20 on (`tail`) as verbatim text, so
   re-ordering and renaming can never perturb a coordinate. Tests assert this —
   do not "normalise" formatting on the way out.
2. **Library code raises `TopologyError` and never exits.** Exiting happens only
   under `cli/` (`common.run()` turns the exception into a message and a status;
   argparse handles usage errors). This is what makes every module importable
   from a script.
3. **Re-indexing must not reorder atoms *within* a bonded interaction.**
   Impropers and other order-sensitive terms depend on it. `ItpFile.permute`
   may sort the rows of a block, never the tokens inside a row.

### Extension points

- **A new `.itp` directive** that carries atom indices: add one entry to
  `itp.DIRECTIVES`. `DirectiveSpec` handles the three shapes (`fixed` N leading
  index columns, `all` of them for `[ exclusions ]`, `vsiten` for
  `[ virtual_sitesn ]`). Do not branch on directive names elsewhere.
- Values outside `--range` clamp to the end colour (matplotlib does this), so a
  capped scale needs `extend` to label that end `>value` — otherwise a
  phosphorus sitting on a scale topped at the nitrogens reads as though it were
  0.92. `draw_itp` defaults to `extend="auto"` and derives it from the drawn
  values (`draw._extend_for`); `colorbar_scene` refuses `"auto"`, having no
  molecule to derive it from.
- `draw.colorbar_scene` builds a colour bar with no molecule attached, for
  figures whose panels are drawn `--bare` and share one bar. It samples the map
  through the same `draw._bar_fills` the inline bar uses, which is what keeps
  the shared bar meaning the same thing as the panels — a test asserts the two
  produce identical blocks. Vertical bars put the caption above and the end
  labels beside the bar rather than rotating text, so no backend needs to grow
  a rotation primitive.
- **A new per-atom property to colour by**: add one `draw.Property` to
  `draw.PROPERTIES` (label, unit, default colour map, whether the scale is
  centred on zero, whether it needs a force field) and, if it is not read off
  the `Atom`, a branch in `draw._values`. Signed quantities get a diverging map
  centred on zero; unsigned ones a sequential map over their own range.
- **Partial charge comes from the molecule's `[ atoms ]`, never from the force
  field's `[ atomtypes ]`.** Both files have a charge column and they routinely
  disagree (Slipids type `NTL` is −0.60; DOPC's N is +0.20) — the molecule wins.
  `AtomType.charge` is parsed and kept, but nothing colours with it.
- **The element colour table** in `colors.py` is transcribed from Jmol's own
  constants (https://jmol.sourceforge.net/jscolors/) — it is data, not a
  judgement call, so correct it against that source rather than by eye. Label
  colours come from the WCAG contrast ratio, not a brightness threshold.
- **Distances between atoms** go through `pbc.Box` (`GroFrame.periodic_box`),
  which handles triclinic boxes; do not reduce by the box diagonal alone.
  `minimum_image` and `distance2` take numpy arrays of shape `(..., 3)` and
  broadcast, so a distance matrix is one call. Only vectors longer than half
  the smallest box height (never bonded neighbours) get the 26-image search.
- **A new tool**: read via `ItpFile` / `GroFrame`, get numbers from
  `GroAtom.position` / `.velocity` (parsed on demand), elements from
  `chem.element_of`, connectivity from `ItpFile.bond_graph()`, and drive it from
  a `cli/` module wired into `topology_wrangling/__main__.py`.

### Drawing gotchas

- `[ atomtypes ]` column layouts vary (the bonded type and atomic number are
  both optional), so `ItpFile._read_atomtypes` locates the lone `A`/`S`/`V`/`D`
  particle-type column and places everything relative to it. Do not count
  columns from the left.

- Sizing runs backwards from the usual: pick a node radius, derive the canvas.
  A fixed canvas shrinks a 138-atom lipid's labels past reading.
- `layout.relax_overlaps` runs only when `draw._has_overlap` finds a real
  collision — it moves atoms, which would spoil the skeletal angles, and a
  skeletal layout usually has nothing to fix. `draw._compose` re-measures the
  extent afterwards because relaxing changes it.
- Hydrogen spokes are deliberately shorter than skeleton bonds, so
  `_radius_units` measures only heavy-atom bonds; sizing from a spoke shrinks
  every atom in an all-atom topology.
- `skeleton.skeletal` is deterministic and takes no seed. `skeleton._score` /
  `_quality` is the single objective, in priority order: **clashes, then
  crossings, then area of page per atom**. Both the hill-climbing and the
  `branch="auto"` choice use it, so they cannot disagree.
- There is **no aspect-ratio term**, on purpose. A straight alkane is thin and
  thin is correct; penalising it folds the chain in half to square the page,
  which costs area and looks wrong. Sprawl is caught by area per atom, not by
  shape. Do not add one back without re-checking `test_a_chain_runs_along_one_axis`.
- `skeleton._unfold` is what stops a branch folding back over another (the head
  group over the glycerol). It cannot be decided while placing — the clash is
  global — so it runs afterwards: find clashing pairs, take the bonds on the
  path between them, and reflect the subtree beyond one across the bond axis.
  **Reflection is a rigid motion**, which is the whole point: bond lengths and
  angles are preserved exactly, so the skeletal geometry survives. Bonds inside
  a ring are never flipped (it would tear the ring open), and a bond that is
  not a bridge has no "beyond" to flip.
- Which branch strategy comes out squarer is not an invariant — unfolding
  changes it. Do not write tests that assume one is.
- A PNG's pixel dimensions do not depend on font size, so comparing them across
  dpi does NOT catch text scaling with resolution; the test that does compares
  the fraction of each image that is ink.

### Performance

`GroFrame` holds one `__slots__` object per atom for the whole frame; real
systems run to hundreds of thousands of atoms. Parse lazily and stay
line-oriented — do not add per-atom float parsing to the read path.

## Test data and golden files

`examples/` is a CHARMM-GUI membrane system plus the topologies describing it,
chosen so the two-step workflow is exercised end to end: the Slipids DOPC
topology lists hydrogens in blocks (fixed by step 1), and the `.gro` uses CHARMM
names where the topologies use the converter's (fixed by step 2, via
`examples/charmm.map`).

`examples/DOPC_reordered.itp` does double duty: it is the **golden reference**
for step 1 and an **input** to step 2. A deliberate change to the ordering
algorithm means regenerating it.

Tests check properties, not just output: every hydrogen follows a heavy atom it
is bonded to, every bonded term still joins the same atoms after renumbering,
re-ordering is idempotent, and the re-ordered `.gro` conserves every coordinate
string of the input.
