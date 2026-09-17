"""Force-field [ atomtypes ] parsing, and painting a parameter onto a drawing."""

import unittest

from _support import example, read
from topology_wrangling import ItpFile, TopologyError
from topology_wrangling import colors, draw
from topology_wrangling.draw import colorbar_scene

FF = example("slipids_ffnonbonded.itp")
MOLECULE = example("DOPC_reordered.itp")


class AtomTypeTest(unittest.TestCase):
    def setUp(self):
        self.types = ItpFile.read(FF).atomtypes

    def test_every_type_is_read(self):
        self.assertEqual(len(self.types), 32)
        self.assertIn("PL", self.types)

    def test_values_match_the_file(self):
        """Spot rows, including one written in exponent form."""
        phosphorus = self.types["PL"]
        self.assertEqual(phosphorus.atomic_number, 15)
        self.assertAlmostEqual(phosphorus.mass, 30.974)
        self.assertAlmostEqual(phosphorus.sigma, 0.3830864488)
        self.assertAlmostEqual(phosphorus.epsilon, 2.44764)
        self.assertEqual(phosphorus.ptype, "A")

        water_oxygen = self.types["OWT3"]        # written as 3.15058e-01
        self.assertAlmostEqual(water_oxygen.sigma, 0.315058)
        self.assertAlmostEqual(water_oxygen.epsilon, 0.636386)

    def test_a_zero_epsilon_is_kept_not_dropped(self):
        """TIP3P hydrogens have no LJ term at all; that is a value, not a gap."""
        self.assertEqual(self.types["HWT3"].epsilon, 0.0)
        self.assertEqual(self.types["HWT3"].sigma, 0.0)

    def test_trailing_comments_do_not_become_columns(self):
        self.assertAlmostEqual(self.types["CTL3"].epsilon, 0.34)
        self.assertIn("New", self.types["CTL3"].comment)

    def test_the_file_still_round_trips(self):
        self.assertEqual(ItpFile.read(FF).text(), read(FF))

    def test_column_layouts(self):
        """GROMACS allows the bonded type and atomic number to be absent."""
        layouts = {
            "name at.num mass charge ptype": " CT 6 12.011 -0.1 A 0.35 0.276",
            "no atomic number": " CT 12.011 -0.1 A 0.35 0.276",
            "bonded type and atomic number": " CT C 6 12.011 -0.1 A 0.35 0.276",
        }
        for why, row in layouts.items():
            with self.subTest(why):
                entry = ItpFile.from_text("[ atomtypes ]\n%s\n" % row).atomtypes
                self.assertAlmostEqual(entry["CT"].epsilon, 0.276)
                self.assertAlmostEqual(entry["CT"].sigma, 0.35)
                self.assertAlmostEqual(entry["CT"].charge, -0.1)
                self.assertAlmostEqual(entry["CT"].mass, 12.011)

    def test_a_file_without_the_directive_is_reported(self):
        with self.assertRaises(TopologyError):
            ItpFile.read(MOLECULE).atomtypes

    def test_a_malformed_row_is_reported(self):
        for row in (" CT 6 12.011", " CT 6 12.011 -0.1 0.35 0.276"):
            with self.subTest(row), self.assertRaises(TopologyError):
                ItpFile.from_text("[ atomtypes ]\n%s\n" % row).atomtypes


class ScaleTest(unittest.TestCase):
    def test_a_signed_quantity_is_centred_on_zero(self):
        _, vmin, vmax = colors.value_colors([-0.8, 0.2], symmetric=True)
        self.assertEqual((vmin, vmax), (-0.8, 0.8))

    def test_an_unsigned_quantity_spans_its_own_range(self):
        _, vmin, vmax = colors.value_colors([0.09, 2.45], symmetric=False)
        self.assertEqual((vmin, vmax), (0.09, 2.45))

    def test_a_constant_property_still_gives_a_usable_bar(self):
        _, vmin, vmax = colors.value_colors([0.3, 0.3])
        self.assertLess(vmin, vmax)


class OverlayTest(unittest.TestCase):
    def setUp(self):
        self.molecule = ItpFile.read(MOLECULE)
        self.types = ItpFile.read(FF).atomtypes

    def draw(self, **kwargs):
        kwargs.setdefault("hide_hydrogens", True)
        return draw.draw_topology(self.molecule, **kwargs)

    def test_each_property_uses_the_colour_map_that_suits_it(self):
        self.assertEqual(draw.PROPERTIES["charge"].cmap, "coolwarm")
        self.assertTrue(draw.PROPERTIES["charge"].symmetric)
        self.assertEqual(draw.PROPERTIES["epsilon"].cmap, "viridis")
        self.assertFalse(draw.PROPERTIES["epsilon"].symmetric)

    def test_epsilon_spans_the_force_field_values(self):
        scene = self.draw(color_by="epsilon", atomtypes=self.types)
        caption, vmin, vmax = _legend(scene)
        self.assertIn("epsilon", caption)
        self.assertIn("kJ/mol", caption)
        # phosphorus has the deepest well in this force field
        self.assertAlmostEqual(vmax, self.types["PL"].epsilon,
                               places=LEGEND_PLACES)
        self.assertGreater(vmin, 0.0)

    def test_the_colour_bar_names_the_quantity_not_the_colour_map(self):
        scene = self.draw(color_by="epsilon", atomtypes=self.types)
        captions = [t[2] for t in scene.texts]
        self.assertIn("LJ epsilon (kJ/mol)", captions)
        for name in ("viridis", "coolwarm", "plasma"):
            self.assertNotIn(name, " ".join(captions))

    def test_charge_comes_from_the_molecule_not_the_force_field(self):
        """A charge in [ atoms ] overrides the force field's type default."""
        nitrogen = next(a for a in self.molecule.atoms if a.name == "N")
        self.assertEqual(nitrogen.type, "NTL")
        self.assertNotAlmostEqual(float(nitrogen.charge),
                                  self.types["NTL"].charge)
        _, vmin, vmax = _legend(self.draw(color_by="charge",
                                          atomtypes=self.types))
        molecule_extreme = max(abs(float(a.charge))
                               for a in self.molecule.atoms
                               if a.name in {b.name for b in self.molecule.atoms})
        self.assertAlmostEqual(vmax, molecule_extreme, places=LEGEND_PLACES)
        self.assertAlmostEqual(vmin, -molecule_extreme, places=LEGEND_PLACES)

    def test_the_deepest_well_and_the_shallowest_differ_in_colour(self):
        scene = self.draw(color_by="epsilon", atomtypes=self.types)
        fills = [c[3] for c in scene.circles]
        self.assertGreater(len(set(fills)), 3)

    def test_an_explicit_colour_map_overrides_the_default(self):
        """Checked by the colours themselves, since the bar no longer names it."""
        default = self.draw(color_by="epsilon", atomtypes=self.types)
        plasma = self.draw(color_by="epsilon", atomtypes=self.types,
                           cmap="plasma")
        self.assertNotEqual([c[3] for c in default.circles],
                            [c[3] for c in plasma.circles])
        self.assertEqual(_legend(default)[1:], _legend(plasma)[1:])

    def test_an_explicit_range_is_honoured(self):
        """Fixing the range is what puts two molecules on one scale."""
        scene = self.draw(color_by="epsilon", atomtypes=self.types,
                          value_range=(0.0, 5.0))
        self.assertEqual(_legend(scene)[1:], (0.0, 5.0))

    def test_epsilon_without_a_force_field_is_reported(self):
        with self.assertRaises(TopologyError) as caught:
            self.draw(color_by="epsilon")
        self.assertIn("--ff", str(caught.exception))

    def test_an_atom_type_missing_from_the_force_field_is_named(self):
        partial = {k: v for k, v in self.types.items() if k != "PL"}
        with self.assertRaises(TopologyError) as caught:
            self.draw(color_by="epsilon", atomtypes=partial)
        self.assertIn("PL", str(caught.exception))

    def test_a_bare_drawing_carries_no_annotation(self):
        """For panels of one figure that share a single colour bar."""
        full = self.draw(color_by="epsilon", atomtypes=self.types)
        bare = self.draw(color_by="epsilon", atomtypes=self.types,
                         title=False, colorbar=False)
        self.assertFalse(bare.rects, "the colour bar is still drawn")
        labels = [t[2] for t in bare.texts]
        self.assertNotIn("LJ epsilon (kJ/mol)", labels)
        self.assertFalse([t for t in labels if "atoms" in t])
        # the structure itself is untouched
        self.assertEqual(len(bare.circles), len(full.circles))
        self.assertEqual([c[3] for c in bare.circles],
                         [c[3] for c in full.circles])
        self.assertLess(bare.height, full.height)

    def test_title_and_colorbar_can_be_dropped_separately(self):
        atomtypes = self.types
        no_bar = self.draw(color_by="epsilon", atomtypes=atomtypes,
                           colorbar=False)
        self.assertFalse(no_bar.rects)
        self.assertTrue([t for t in no_bar.texts if "atoms" in t[2]])

        no_title = self.draw(color_by="epsilon", atomtypes=atomtypes,
                             title=False)
        self.assertTrue(no_title.rects)
        self.assertFalse([t for t in no_title.texts if "atoms" in t[2]])

    def test_a_value_past_the_range_is_drawn_in_the_end_colour(self):
        """Capping the scale below the phosphorus is the point of extend."""
        capped = self.draw(color_by="epsilon", atomtypes=self.types,
                           value_range=(0.0293, 0.9205))
        fills = {c[3] for c in capped.circles}
        top, _, _ = colors.value_colors([0.9205], "viridis", 0.0293, 0.9205)
        self.assertIn(top[0], fills)          # phosphorus, clamped to the top

    def test_extend_is_worked_out_from_the_atoms_drawn(self):
        # phosphorus (2.45) runs past a scale topped at the nitrogens
        capped = self.draw(color_by="epsilon", atomtypes=self.types,
                           value_range=(0.0293, 0.9205))
        self.assertIn(">0.92", [t[2] for t in capped.texts])
        # a scale that covers everything says nothing
        full = self.draw(color_by="epsilon", atomtypes=self.types)
        self.assertFalse([t for t in full.texts if t[2].startswith((">", "<"))])

    def test_extend_can_be_set_by_hand(self):
        scene = self.draw(color_by="epsilon", atomtypes=self.types,
                          value_range=(0.0293, 0.9205), extend="neither")
        self.assertNotIn(">0.92", [t[2] for t in scene.texts])
        self.assertIn("0.92", [t[2] for t in scene.texts])

    def test_an_unknown_property_is_rejected(self):
        with self.assertRaises(TopologyError):
            self.draw(color_by="polarisability")


# The colour bar prints its ends to three significant figures, so values read
# back out of a finished scene carry only that much precision.
LEGEND_PLACES = 2


class StandaloneColorbarTest(unittest.TestCase):
    def test_it_matches_the_bar_drawn_under_a_molecule(self):
        """The shared bar has to mean the same thing as the panels it labels."""
        molecule = ItpFile.read(MOLECULE)
        types = ItpFile.read(FF).atomtypes
        inline = draw.draw_topology(molecule, color_by="epsilon",
                                    atomtypes=types, hide_hydrogens=True,
                                    value_range=(0.05, 2.0))
        alone = colorbar_scene("epsilon", value_range=(0.05, 2.0))
        # The two are drawn in different numbers of blocks (each scales its
        # smoothness to its own length), so compare the colours themselves at
        # matched points along the bar rather than block for block.
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            with self.subTest(fraction=fraction):
                self.assertEqual(_at(alone.rects, fraction),
                                 _at(inline.rects, fraction))
        self.assertIn("LJ epsilon (kJ/mol)", [t[2] for t in alone.texts])

    def test_a_longer_bar_is_drawn_in_more_blocks(self):
        """Smoothness has to scale, or a large bar bands visibly."""
        short = colorbar_scene("charge", value_range=(-1.0, 1.0), length=120)
        long = colorbar_scene("charge", value_range=(-1.0, 1.0), length=1200)
        self.assertGreater(len(long.rects), len(short.rects))
        # blocks stay under a scene unit wide, so they cannot band
        self.assertLess(1200.0 / len(long.rects), 1.0)

    def test_both_orientations(self):
        wide = colorbar_scene("charge", value_range=(-1.5, 1.5))
        tall = colorbar_scene("charge", value_range=(-1.5, 1.5),
                              orientation="vertical")
        self.assertGreater(wide.width, wide.height)
        self.assertGreater(tall.height, tall.width)
        for scene in (wide, tall):
            labels = [t[2] for t in scene.texts]
            self.assertIn("-1.5", labels)
            self.assertIn("1.5", labels)
            self.assertIn("partial charge (e)", labels)

    def test_the_low_end_is_at_the_bottom_when_vertical(self):
        scene = colorbar_scene("epsilon", value_range=(0.0, 1.0),
                               orientation="vertical")
        lowest = max(scene.rects, key=lambda r: r[1])      # largest y = bottom
        highest = min(scene.rects, key=lambda r: r[1])
        ends, _, _ = colors.value_colors([0.0, 1.0], "viridis", 0.0, 1.0)
        self.assertEqual(lowest[4], ends[0])
        self.assertEqual(highest[4], ends[1])

    def test_extend_marks_the_end_the_scale_does_not_reach(self):
        for extend, low, high in (("neither", "0", "1"),
                                  ("max", "0", ">1"),
                                  ("min", "<0", "1"),
                                  ("both", "<0", ">1")):
            with self.subTest(extend):
                scene = colorbar_scene("epsilon", value_range=(0.0, 1.0),
                                       extend=extend)
                labels = [t[2] for t in scene.texts]
                self.assertIn(low, labels)
                self.assertIn(high, labels)

    def test_a_standalone_bar_will_not_guess_extend(self):
        """There is no molecule here to work it out from."""
        with self.assertRaises(TopologyError):
            colorbar_scene("charge", value_range=(0.0, 1.0), extend="auto")

    def test_a_range_is_required_and_must_increase(self):
        for kwargs in ({"value_range": (None, None)},
                       {"value_range": (1.0, None)},
                       {"value_range": (1.0, 1.0)},
                       {"value_range": (2.0, 1.0)}):
            with self.subTest(kwargs), self.assertRaises(TopologyError):
                colorbar_scene("charge", **kwargs)

    def test_bad_inputs_are_rejected(self):
        with self.assertRaises(TopologyError):
            colorbar_scene("polarisability", value_range=(0.0, 1.0))
        with self.assertRaises(TopologyError):
            colorbar_scene("charge", value_range=(0.0, 1.0),
                           orientation="sideways")

    def test_an_explicit_colour_map_is_used(self):
        default = colorbar_scene("charge", value_range=(-1.0, 1.0))
        plasma = colorbar_scene("charge", value_range=(-1.0, 1.0),
                                cmap="plasma")
        self.assertNotEqual([r[4] for r in default.rects],
                            [r[4] for r in plasma.rects])


def _at(rects, fraction):
    """The colour a fraction of the way along a bar's blocks."""
    index = min(int(round(fraction * (len(rects) - 1))), len(rects) - 1)
    return rects[index][4]


def _legend(scene):
    """Pull the colour-bar caption and range back out of a finished scene."""
    labels = [t[2] for t in scene.texts]
    caption = next(t for t in labels if "(" in t and ")" in t
                   and not t.startswith("DOPC"))
    numbers = []
    for text in labels:
        try:
            numbers.append(float(text))
        except ValueError:
            continue
    return caption, min(numbers), max(numbers)


if __name__ == "__main__":
    unittest.main()
