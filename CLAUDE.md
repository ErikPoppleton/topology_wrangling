# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

Pure standard library, Python 3. No build step, no dependencies, no installer.

```sh
# all tests
python3 -m unittest discover -s tests

# one test, by name pattern, from the repo root
python3 -m unittest discover -s tests -k test_permutation_is_complete

# one module or class (must run from inside tests/, which puts _support.py on sys.path)
cd tests && python3 -m unittest test_itp.UpdateGroupOrderTest

# the two tools, on the example system
./update_group_ordering.py examples/update_group_incompatible.itp -o DOPC.itp
./reorder_gro.py examples/charmm_DOPC.gro DOPC.itp examples/TP3.itp examples/Na+.itp \
    examples/Cl-.itp examples/MG.itp examples/CA.itp --alias-file examples/charmm.map -o system.gro
```

`python3 -m unittest tests.test_itp` does **not** work — the tests import `_support`
by bare name, so the `tests/` directory itself has to be on `sys.path`.

## Architecture

Both tools are the same pipeline: **read a topology → derive an `AtomOrder` →
apply it to something**. They differ only in where the order comes from and
what it is applied to. Keep new tools on that spine.

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
- **A new tool**: read via `ItpFile` / `GroFrame`, get numbers from
  `GroAtom.position` / `.velocity` (parsed on demand), elements from
  `chem.element_of`, connectivity from `ItpFile.bond_graph()`, and drive it from
  a `cli/` module wired into `topology_wrangling/__main__.py`.

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
