"""The two command line tools, driven end to end on the example system.

This is the workflow the examples directory exists for:

  1. reorder the Slipids DOPC topology for update groups, then
  2. reorder and rename the CHARMM .gro to match it and the other topologies.
"""

import contextlib
import io
import os
import shutil
import tempfile
import unittest

from _support import example, read, SYSTEM
from topology_wrangling import AtomOrder, GroFrame, ItpFile, TopologyError
from topology_wrangling.cli import (amber_style_index, draw_colorbar, draw_itp,
                                    reorder_gro, update_group_ordering)
from topology_wrangling.cli.common import run


@contextlib.contextmanager
def quiet():
    """The tools report progress on stderr; keep it out of the test output."""
    with contextlib.redirect_stderr(io.StringIO()) as captured:
        yield captured


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="topowrangle-")
        self.itp = os.path.join(self.dir, "DOPC.itp")
        shutil.copy(example("update_group_incompatible.itp"), self.itp)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def path(self, name):
        return os.path.join(self.dir, name)

    def test_step_one_reorders_the_topology(self):
        with quiet():
            self.assertEqual(update_group_ordering.main([self.itp]), 0)
        out = self.path("DOPC_reordered.itp")
        self.assertTrue(os.path.exists(out), "the default output name changed")
        self.assertEqual(read(out), read(example("DOPC_reordered.itp")))

    def test_step_two_reorders_and_renames_the_coordinates(self):
        out = self.path("system.gro")
        with quiet():
            update_group_ordering.main([self.itp, "-o", self.path("top.itp")])
            code = reorder_gro.main(
                [example("charmm_DOPC.gro"), self.path("top.itp")]
                + [example(name) for name, _, _, _ in SYSTEM
                   if name != "DOPC_reordered.itp"]
                + ["--alias-file", example("charmm.map"), "-o", out]
            )
        self.assertEqual(code, 0)

        frame = GroFrame.read(out)
        source = GroFrame.read(example("charmm_DOPC.gro"))
        self.assertEqual(len(frame), len(source))
        self.assertEqual(sorted(a.tail for a in frame),
                         sorted(a.tail for a in source))
        # Every molecule now matches its topology in order and in naming.
        for name, _, size, count in SYSTEM:
            path = self.path("top.itp") if name == "DOPC_reordered.itp" \
                else example(name)
            top = ItpFile.read(path)
            blocks = list(frame.iter_residue_blocks(top.residue_name(), size))
            with self.subTest(name):
                self.assertEqual(len(blocks), count)
                for _, block in blocks:
                    self.assertEqual([a.name for a in block], top.atom_names)

    def test_the_map_carries_the_permutation_to_a_second_run(self):
        """--write-map then --read-map must reproduce the re-ordering."""
        with quiet():
            update_group_ordering.main([self.itp, "-o", self.path("top.itp"),
                              "--write-map", self.path("map.txt")])
        text = read(self.path("map.txt"))
        self.assertIn("resname DOPC", text)
        self.assertIn("natoms 138", text)

        # --read-map expects coordinates in the ORIGINAL topology order, so
        # build some by undoing the permutation on the example's DOPC atoms.
        order = AtomOrder.read_map(self.path("map.txt"))
        frame = GroFrame.read(example("charmm_DOPC.gro"))
        frame.atoms = [a for a in frame if a.resname == "DOPC"]
        inverse = order.inverse()
        frame.atoms = [a for start, block in
                       frame.iter_residue_blocks("DOPC", 138)
                       for a in inverse.apply(block)]
        frame.renumber()
        frame.write(self.path("original_order.gro"))

        with quiet():
            update_group_ordering.main(["--read-map", self.path("map.txt"),
                              "--gro", self.path("original_order.gro"),
                              "--gro-out", self.path("back.gro")])

        # Re-ordering by the map must give back the example's own order.
        back = GroFrame.read(self.path("back.gro"))
        expected = [a for a in GroFrame.read(example("charmm_DOPC.gro"))
                    if a.resname == "DOPC"]
        self.assertEqual([a.name for a in back], [a.name for a in expected])
        self.assertEqual([a.tail for a in back], [a.tail for a in expected])

    def test_default_output_names(self):
        with quiet():
            update_group_ordering.main([self.itp])
        self.assertTrue(os.path.exists(self.path("DOPC_reordered.itp")))

    def test_charge_group_option(self):
        with quiet():
            update_group_ordering.main([self.itp, "--cgnr", "per-atom",
                              "-o", self.path("per_atom.itp")])
        atoms = ItpFile.read(self.path("per_atom.itp")).atoms
        self.assertEqual([a.cgnr for a in atoms], [str(i) for i in range(1, 139)])


class DrawToolTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="topowrangle-draw-")

    def tearDown(self):
        shutil.rmtree(self.dir)

    def test_default_output_is_svg_beside_the_input(self):
        itp = os.path.join(self.dir, "TP3.itp")
        shutil.copy(example("TP3.itp"), itp)
        with quiet():
            self.assertEqual(draw_itp.main([itp, "--iterations", "20"]), 0)
        self.assertTrue(os.path.exists(os.path.join(self.dir, "TP3.svg")))

    def test_png_output(self):
        out = os.path.join(self.dir, "water.png")
        with quiet():
            draw_itp.main([example("TP3.itp"), "-o", out, "--iterations", "20",
                           "--dpi", "72"])
        self.assertGreater(os.path.getsize(out), 0)

    def test_options_reach_the_drawing(self):
        out = os.path.join(self.dir, "dopc.svg")
        with quiet():
            draw_itp.main([example("DOPC_reordered.itp"), "-o", out,
                           "--color-by", "charge", "--cmap", "bwr",
                           "--hide-hydrogens", "--label", "type",
                           "--iterations", "20"])
        text = read(out)
        self.assertIn("<circle", text)
        types = {a.type for a in ItpFile.read(example("DOPC_reordered.itp")).atoms}
        self.assertTrue(any(">%s<" % t in text for t in types))

    def test_a_force_field_parameter_overlay(self):
        out = os.path.join(self.dir, "eps.svg")
        with quiet():
            code = draw_itp.main([
                example("DOPC_reordered.itp"), "-o", out,
                "--color-by", "epsilon",
                "--ff", example("slipids_ffnonbonded.itp"),
                "--hide-hydrogens",
            ])
        self.assertEqual(code, 0)
        text = read(out)
        self.assertIn("LJ epsilon (kJ/mol)", text)
        self.assertNotIn("viridis", text)

    def test_epsilon_without_a_force_field_is_refused(self):
        with self.assertRaises(TopologyError), quiet():
            draw_itp.main([example("DOPC_reordered.itp"), "-o",
                           os.path.join(self.dir, "x.svg"),
                           "--color-by", "epsilon"])

    def test_bare_panels_and_a_shared_colour_bar(self):
        """The workflow for putting several molecules in one figure."""
        panel = os.path.join(self.dir, "panel.svg")
        bar = os.path.join(self.dir, "bar.svg")
        with quiet():
            draw_itp.main([example("DOPC_reordered.itp"), "-o", panel,
                           "--color-by", "charge", "--range", "-1.58", "1.58",
                           "--bare", "--hide-hydrogens"])
            draw_colorbar.main(["--color-by", "charge", "--range",
                                "-1.58", "1.58", "-o", bar])
        panel_text = read(panel)
        self.assertNotIn("partial charge", panel_text)
        self.assertNotIn("atoms", panel_text)
        self.assertIn("<circle", panel_text)
        self.assertIn("partial charge (e)", read(bar))

    def test_a_capped_scale_is_marked_on_the_shared_bar(self):
        bar = os.path.join(self.dir, "capped.svg")
        with quiet():
            draw_colorbar.main(["--color-by", "epsilon", "--range",
                                "0.0293", "0.9205", "--extend", "max",
                                "-o", bar])
        self.assertIn("&gt;0.92", read(bar))

    def test_a_colour_bar_needs_its_range(self):
        with self.assertRaises(SystemExit), quiet():
            draw_colorbar.main(["--color-by", "charge", "-o",
                                os.path.join(self.dir, "bar.svg")])

    def test_an_unwritable_format_is_reported(self):
        with self.assertRaises(TopologyError), quiet():
            draw_itp.main([example("TP3.itp"), "-o",
                           os.path.join(self.dir, "x.tiff")])


def read_ndx(path):
    groups, name = {}, None
    for line in read(path).splitlines():
        if line.startswith("["):
            name = line.strip("[] ")
            groups[name] = []
        elif line.strip():
            groups[name].extend(int(i) for i in line.split())
    return groups


class AmberIndexTest(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="topowrangle-ndx-")
        self.gro = os.path.join(self.dir, "system.gro")
        shutil.copy(example("charmm_DOPC.gro"), self.gro)

    def tearDown(self):
        shutil.rmtree(self.dir)

    def path(self, name):
        return os.path.join(self.dir, name)

    def test_default_output_is_beside_the_input(self):
        with quiet():
            self.assertEqual(amber_style_index.main([self.gro]), 0)
        groups = read_ndx(self.path("system_amber.ndx"))
        self.assertEqual(list(groups), ["PC", "OL"])
        self.assertEqual(len(groups["PC"]), 128 * 38)
        self.assertEqual(len(groups["OL"]), 128 * 100)
        self.assertEqual(groups["PC"][:3], [1, 2, 3])

    def test_a_topology_gives_the_same_groups(self):
        with quiet():
            amber_style_index.main([self.gro, "-o", self.path("dist.ndx")])
            amber_style_index.main([self.gro, example("DOPC_reordered.itp"),
                                    "-o", self.path("bond.ndx")])
        self.assertEqual(read(self.path("dist.ndx")),
                         read(self.path("bond.ndx")))

    def test_append_keeps_the_existing_groups(self):
        out = self.path("index.ndx")
        with open(out, "w") as handle:
            handle.write("[ System ]\n1 2 3\n\n")
        with quiet():
            amber_style_index.main([self.gro, "-o", out, "--append"])
        self.assertEqual(list(read_ndx(out)), ["System", "PC", "OL"])

    def test_a_non_lipid_topology_is_an_error(self):
        with self.assertRaises(SystemExit) as caught, quiet() as err:
            run(lambda: amber_style_index.main([self.gro,
                                                example("TP3.itp")]))
        self.assertEqual(caught.exception.code, 1)
        self.assertIn("not one of DOPC/DMPC/DMPG", err.getvalue())


class ArgumentTest(unittest.TestCase):
    def test_update_group_ordering_needs_an_input(self):
        with self.assertRaises(SystemExit), quiet():
            update_group_ordering.main([])

    def test_read_map_without_gro_is_refused(self):
        with self.assertRaises(SystemExit), quiet():
            update_group_ordering.main(["--read-map", "map.txt"])

    def test_resname_with_several_topologies_is_refused(self):
        with self.assertRaises(SystemExit), quiet():
            reorder_gro.main([example("charmm_DOPC.gro"), example("TP3.itp"),
                              example("MG.itp"), "--resname", "TIP3"])

    def test_resname_is_a_shorthand_for_residue(self):
        out = tempfile.mkdtemp(prefix="topowrangle-")
        try:
            with quiet():
                code = reorder_gro.main([
                    example("charmm_DOPC.gro"), example("MG.itp"),
                    "--resname", "MG", "-o", os.path.join(out, "x.gro"),
                ])
            self.assertEqual(code, 0)
        finally:
            shutil.rmtree(out)

    def test_unknown_out_resname_is_rejected(self):
        with self.assertRaises((SystemExit, TopologyError)), quiet():
            reorder_gro.main([example("charmm_DOPC.gro"), example("MG.itp"),
                              "--out-resname", "NOPE=X"])


if __name__ == "__main__":
    unittest.main()
