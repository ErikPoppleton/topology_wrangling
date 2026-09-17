"""Colours, graph layout and the 2D drawing."""

import os
import shutil
import tempfile
import unittest
import xml.etree.ElementTree as ElementTree

from _support import example
from topology_wrangling import ItpFile, TopologyError
from topology_wrangling import colors, draw
from topology_wrangling.layout import force_directed, relax_overlaps

DOPC = example("DOPC_reordered.itp")
WATER = example("TP3.itp")


class JmolColorTest(unittest.TestCase):
    def test_table_matches_jmol(self):
        """Spot values from Jmol's own element colour table."""
        for symbol, expected in (("H", "#FFFFFF"), ("C", "#909090"),
                                 ("N", "#3050F8"), ("O", "#FF0D0D"),
                                 ("P", "#FF8000"), ("S", "#FFFF30"),
                                 ("Cl", "#1FF01F"), ("Na", "#AB5CF2"),
                                 ("Fe", "#E06633")):
            with self.subTest(symbol):
                self.assertEqual(colors.element_color(symbol), expected)

    def test_every_entry_is_a_hex_colour(self):
        self.assertGreaterEqual(len(colors.JMOL_COLORS), 100)
        for symbol, value in colors.JMOL_COLORS.items():
            with self.subTest(symbol):
                self.assertRegex(value, r"^#[0-9A-F]{6}$")

    def test_an_unknown_element_still_gets_a_colour(self):
        self.assertEqual(colors.element_color("Xx"), colors.UNKNOWN_COLOR)


class ContrastTest(unittest.TestCase):
    def test_label_colour_follows_the_background(self):
        """Dark fills take white text, light fills black."""
        for background, expected in (("#FFFFFF", "#000000"),   # hydrogen
                                     ("#FFFF30", "#000000"),   # sulphur
                                     ("#3050F8", "#ffffff"),   # nitrogen
                                     ("#000000", "#ffffff"),
                                     ("#1FF01F", "#000000")):  # chlorine
            with self.subTest(background):
                self.assertEqual(colors.text_color(background), expected)

    def test_luminance_bounds(self):
        self.assertAlmostEqual(colors.relative_luminance("#000000"), 0.0)
        self.assertAlmostEqual(colors.relative_luminance("#FFFFFF"), 1.0)

    def test_a_bad_colour_is_rejected(self):
        with self.assertRaises(TopologyError):
            colors.to_rgb("#12345")


class ChargeColorTest(unittest.TestCase):
    def test_scale_is_symmetric_about_zero(self):
        """So the middle of a diverging map is a neutral atom."""
        _, vmin, vmax = colors.charge_colors([-0.8, 0.2, 0.1])
        self.assertEqual((vmin, vmax), (-0.8, 0.8))

    def test_sign_is_visible_in_the_colour(self):
        (low, mid, high), _, _ = colors.charge_colors([-1.0, 0.0, 1.0],
                                                      "coolwarm")
        self.assertLess(colors.relative_luminance(mid),
                        1.0)                       # not white-out
        self.assertNotEqual(low, high)
        self.assertLess(colors.to_rgb(low)[0], colors.to_rgb(high)[0])  # blue->red

    def test_an_explicit_range_is_honoured(self):
        _, vmin, vmax = colors.charge_colors([0.1], vmin=-2.0, vmax=2.0)
        self.assertEqual((vmin, vmax), (-2.0, 2.0))

    def test_an_unknown_colour_map_is_rejected(self):
        with self.assertRaises(TopologyError):
            colors.charge_colors([0.0], cmap="not-a-colour-map")


class LayoutTest(unittest.TestCase):
    def setUp(self):
        self.itp = ItpFile.read(DOPC)
        self.edges = self.itp.connectivity()

    def test_is_deterministic_for_a_seed(self):
        a = force_directed(self.itp.natoms, self.edges, seed=7, iterations=50)
        b = force_directed(self.itp.natoms, self.edges, seed=7, iterations=50)
        self.assertEqual(a, b)

    def test_bonded_atoms_end_up_closer_than_average(self):
        """The point of a spring embedding: bonds are short."""
        import math
        pos = force_directed(self.itp.natoms, self.edges, seed=0)
        bonded = [math.dist(pos[i - 1], pos[j - 1]) for i, j in self.edges]
        allpairs = [math.dist(pos[i], pos[j])
                    for i in range(0, len(pos), 7)
                    for j in range(i + 1, len(pos), 7)]
        self.assertLess(sum(bonded) / len(bonded),
                        sum(allpairs) / len(allpairs))

    def test_coordinates_are_finite_and_bounded(self):
        pos = force_directed(self.itp.natoms, self.edges, seed=0)
        self.assertEqual(len(pos), self.itp.natoms)
        for x, y in pos:
            self.assertTrue(-1.001 <= x <= 1.001 and -1.001 <= y <= 1.001)

    def test_degenerate_molecules(self):
        self.assertEqual(force_directed(1, []), [(0.0, 0.0)])
        self.assertEqual(len(force_directed(3, [], iterations=20)), 3)
        with self.assertRaises(TopologyError):
            force_directed(0, [])

    def test_a_bond_outside_the_molecule_is_rejected(self):
        with self.assertRaises(TopologyError):
            force_directed(3, [(1, 9)])

    def test_overlaps_are_relaxed(self):
        import math
        crowded = [(0.0, 0.0), (0.01, 0.0), (0.0, 0.01)]
        spread = relax_overlaps(crowded, [0.1, 0.1, 0.1])
        for i in range(3):
            for j in range(i + 1, 3):
                self.assertGreater(math.dist(spread[i], spread[j]), 0.15)


class DrawTest(unittest.TestCase):
    def setUp(self):
        self.itp = ItpFile.read(DOPC)

    def test_one_circle_per_atom(self):
        scene = draw.draw_topology(self.itp, iterations=50)
        self.assertEqual(len(scene.circles), self.itp.natoms)
        self.assertEqual(len(scene.lines), len(self.itp.connectivity()))

    def test_labels_come_from_the_topology(self):
        names = {text for _, _, text, _, _, _
                 in draw.draw_topology(self.itp, iterations=50).texts}
        self.assertTrue(set(self.itp.atom_names) <= names)

        types = {text for _, _, text, _, _, _
                 in draw.draw_topology(self.itp, label="type",
                                       iterations=50).texts}
        self.assertTrue({a.type for a in self.itp.atoms} <= types)

    def test_label_none_writes_only_the_title(self):
        scene = draw.draw_topology(self.itp, label="none", iterations=50)
        self.assertEqual(len(scene.texts), 1)
        self.assertIn("DOPC", scene.texts[0][2])

    def test_hydrogens_can_be_hidden(self):
        scene = draw.draw_topology(self.itp, hide_hydrogens=True, iterations=50)
        self.assertEqual(len(scene.circles), 54)

    def test_element_colouring_uses_jmol(self):
        scene = draw.draw_topology(self.itp, iterations=50)
        fills = {fill for _, _, _, fill, _, _ in scene.circles}
        self.assertEqual(fills, {colors.element_color(e)
                                 for e in ("C", "H", "N", "O", "P")})

    def test_charge_colouring_adds_a_legend(self):
        plain = draw.draw_topology(self.itp, iterations=50)
        charged = draw.draw_topology(self.itp, color_by="charge",
                                     iterations=50)
        self.assertFalse(plain.rects)            # no colour bar
        self.assertTrue(charged.rects)           # a colour bar
        self.assertGreater(len({f for _, _, _, f, _, _ in charged.circles}), 10)

    def test_water_is_connected_through_settles(self):
        """Rigid water has no [ bonds ]; it must not draw as loose dots."""
        scene = draw.draw_topology(ItpFile.read(WATER), iterations=50)
        self.assertEqual(len(scene.circles), 3)
        self.assertEqual(len(scene.lines), 2)

    def test_bad_options_are_rejected(self):
        for kwargs in ({"label": "nonsense"}, {"color_by": "nonsense"}):
            with self.subTest(kwargs), self.assertRaises(TopologyError):
                draw.draw_topology(self.itp, iterations=10, **kwargs)

    def test_hiding_every_atom_is_an_error(self):
        water = ItpFile.from_text(
            "[ atoms ]\n 1 HW 1 W H1 1 0.4 1.008\n 2 HW 1 W H2 2 0.4 1.008\n"
        )
        with self.assertRaises(TopologyError):
            draw.draw_topology(water, hide_hydrogens=True)


class OutputTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="topowrangle-draw-")
        self.scene = draw.draw_topology(ItpFile.read(WATER), iterations=50)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_svg_is_well_formed(self):
        path = os.path.join(self.dir, "water.svg")
        self.scene.save(path)
        root = ElementTree.parse(path).getroot()
        self.assertTrue(root.tag.endswith("svg"))
        circles = root.iter("{http://www.w3.org/2000/svg}circle")
        self.assertEqual(len(list(circles)), 3)

    def test_svg_escapes_label_text(self):
        itp = ItpFile.from_text("[ atoms ]\n 1 C 1 R A<&B 1 0.0 12.011\n")
        scene = draw.draw_topology(itp, iterations=5)
        self.assertIn("A&lt;&amp;B", scene.to_svg())
        ElementTree.fromstring(scene.to_svg())      # still parses

    def test_png_is_written(self):
        path = os.path.join(self.dir, "water.png")
        self.scene.save(path, dpi=72)
        with open(path, "rb") as handle:
            self.assertEqual(handle.read(8), b"\x89PNG\r\n\x1a\n")

    def test_png_content_does_not_grow_with_dpi(self):
        """Text and line widths are in points and must not scale with dpi.

        Sizing them by dpi made labels grow with resolution and disagree with
        the SVG, which is the one thing the shared Scene exists to prevent.
        Comparing how much of each image is ink catches that; comparing the
        pixel dimensions does not, since those depend only on the figure size.
        """
        from matplotlib.image import imread

        low = os.path.join(self.dir, "low.png")
        high = os.path.join(self.dir, "high.png")
        self.scene.save(low, dpi=60)
        self.scene.save(high, dpi=180)
        self.assertAlmostEqual(_ink_fraction(imread(low)),
                               _ink_fraction(imread(high)), delta=0.02)

    def test_png_pixel_size_follows_dpi(self):
        low = os.path.join(self.dir, "low.png")
        high = os.path.join(self.dir, "high.png")
        self.scene.save(low, dpi=72)
        self.scene.save(high, dpi=144)
        wide, tall = _png_size(low)
        # a half-pixel of the canvas can round either way at each resolution,
        # in either dimension
        self.assertAlmostEqual(_png_size(high)[0], 2 * wide, delta=1)
        self.assertAlmostEqual(_png_size(high)[1], 2 * tall, delta=1)

    def test_an_unknown_extension_is_rejected(self):
        with self.assertRaises(TopologyError):
            self.scene.save(os.path.join(self.dir, "water.pdf"))


def _ink_fraction(pixels):
    """How much of a rendered image is not background."""
    dark = 0
    for row in pixels:
        for pixel in row:
            if float(pixel[0]) + float(pixel[1]) + float(pixel[2]) < 2.7:
                dark += 1
    return dark / float(len(pixels) * len(pixels[0]))


def _png_size(path):
    """(width, height) from a PNG's IHDR, without needing an image library."""
    import struct
    with open(path, "rb") as handle:
        header = handle.read(24)
    return struct.unpack(">II", header[16:24])


if __name__ == "__main__":
    unittest.main()
