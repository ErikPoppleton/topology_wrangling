"""Coordinate file reading, writing and residue matching."""

import os
import tempfile
import unittest

from _support import example, read, SYSTEM
from topology_wrangling import (AliasTable, GroFrame, ItpFile, ResidueTopology,
                                TopologyError, apply_topologies)

GRO = example("charmm_DOPC.gro")
MAP = example("charmm.map")


def load_system(table):
    return [ResidueTopology.from_itp(ItpFile.read(example(name)), table)
            for name, _, _, _ in SYSTEM]


class GroReadWriteTest(unittest.TestCase):
    def setUp(self):
        self.frame = GroFrame.read(GRO)

    def test_frame_contents(self):
        self.assertEqual(len(self.frame), 35437)
        self.assertEqual(self.frame.title, "Title")
        self.assertEqual(self.frame.resnames,
                         ["DOPC", "SOD", "CAL", "MG", "CLA", "TIP3"])

    def test_round_trip_is_byte_identical(self):
        """Reading and writing an untouched frame must not perturb anything."""
        path = _tmp()
        try:
            self.frame.write(path)
            self.assertEqual(read(path), read(GRO))
        finally:
            os.remove(path)

    def test_columns_are_parsed(self):
        atom = self.frame.atoms[0]
        self.assertEqual((atom.resid, atom.resname, atom.name, atom.nr),
                         ("1", "DOPC", "N", "1"))

    def test_coordinates_and_velocities(self):
        atom = self.frame.atoms[0]
        self.assertEqual(atom.position, (3.436, 1.390, 5.388))
        self.assertEqual(atom.velocity, (-0.0230, 0.3964, -1.0131))

    def test_a_frame_without_velocities(self):
        frame = GroFrame.from_text_for_test = None   # placeholder, see below
        text = ("t\n2\n"
                "    1AAA      C    1   1.000   2.000   3.000\n"
                "    1AAA      H    2   1.100   2.000   3.000\n"
                "   1.0   1.0   1.0\n")
        path = _tmp()
        try:
            with open(path, "w") as handle:
                handle.write(text)
            frame = GroFrame.read(path)
        finally:
            os.remove(path)
        self.assertEqual(frame.atoms[0].position, (1.0, 2.0, 3.0))
        self.assertIsNone(frame.atoms[0].velocity)

    def test_truncated_and_multi_frame_files_are_rejected(self):
        for text, why in (
            ("t\n5\n    1AAA      C    1   1.0   2.0   3.0\n", "truncated"),
            ("t\n1\n    1AAA      C    1   1.0   2.0   3.0\n  1 1 1\nt\n1\n",
             "second frame"),
            ("t\n", "no atom count"),
        ):
            with self.subTest(why):
                path = _tmp()
                try:
                    with open(path, "w") as handle:
                        handle.write(text)
                    with self.assertRaises(TopologyError):
                        GroFrame.read(path)
                finally:
                    os.remove(path)

    def test_residue_blocks(self):
        blocks = list(self.frame.iter_residue_blocks("DOPC", 138))
        self.assertEqual(len(blocks), 128)
        self.assertEqual(blocks[0][0], 0)
        self.assertTrue(all(a.resname == "DOPC" for _, b in blocks for a in b))

    def test_a_wrong_block_size_is_an_error(self):
        with self.assertRaises(TopologyError):
            list(self.frame.iter_residue_blocks("TIP3", 4))

    def test_renumber_wraps_like_gromacs(self):
        frame = GroFrame.read(GRO)
        frame.renumber(start=99999)
        self.assertEqual([a.nr for a in frame.atoms[:3]], ["99999", "0", "1"])


class AliasTableTest(unittest.TestCase):
    def test_alias_file(self):
        table = AliasTable().load(MAP)
        self.assertEqual(table.residues["TIP3"], "TP3")
        self.assertEqual(table.residues["CAL"], "CA")
        self.assertEqual(table.atoms_for("TP3"), {"O": "OH2"})
        self.assertEqual(table.gro_name_for("TP3"), "TIP3")
        self.assertEqual(table.gro_name_for("DOPC"), "DOPC")

    def test_command_line_specs(self):
        table = AliasTable()
        table.add_residue_spec("TIP3=TP3")
        table.add_atom_spec("TP3:O=OH2")
        table.add_atom_spec("O=OH2", default_resname="TP3")
        self.assertEqual(table.residues, {"TIP3": "TP3"})
        self.assertEqual(table.atoms_for("TP3"), {"O": "OH2"})

    def test_malformed_specs_are_rejected(self):
        table = AliasTable()
        for bad in ("TIP3", "TP3:O"):
            with self.subTest(bad), self.assertRaises(TopologyError):
                table.add_residue_spec(bad) if ":" not in bad \
                    else table.add_atom_spec(bad)
        with self.assertRaises(TopologyError):
            table.add_atom_spec("O=OH2")     # no residue given


class ApplyTopologiesTest(unittest.TestCase):
    """The end-to-end job: make the example .gro match the example topologies."""

    @classmethod
    def setUpClass(cls):
        cls.table = AliasTable().load(MAP)
        cls.source = GroFrame.read(GRO)
        cls.frame = GroFrame.read(GRO)
        cls.counts = apply_topologies(cls.frame, load_system(cls.table))

    def test_every_molecule_was_matched(self):
        expected = {gro_name: count for _, gro_name, _, count in SYSTEM}
        self.assertEqual(self.counts, expected)

    def test_atom_order_and_names_now_match_the_topologies(self):
        for name, _, size, count in SYSTEM:
            top = ItpFile.read(example(name))
            blocks = list(self.frame.iter_residue_blocks(top.residue_name(), size))
            with self.subTest(name):
                self.assertEqual(len(blocks), count)
                for _, block in blocks:
                    self.assertEqual([a.name for a in block], top.atom_names)

    def test_no_coordinate_was_changed_lost_or_duplicated(self):
        self.assertEqual(sorted(a.tail for a in self.frame),
                         sorted(a.tail for a in self.source))

    def test_coordinates_travel_with_their_atom(self):
        """Renaming must not detach an atom from its position."""
        water_in = {a.name: a.tail for a in self.source.atoms[17710:17713]}
        first = [a for a in self.frame if a.resname == "TP3"][:3]
        self.assertEqual([a.name for a in first], ["O", "H1", "H2"])
        self.assertEqual(first[0].tail, water_in["OH2"])
        self.assertEqual(first[1].tail, water_in["H1"])

    def test_the_whole_system_is_covered(self):
        self.assertEqual(sorted(self.frame.resnames),
                         sorted(["DOPC", "TP3", "Na+", "Cl-", "MG", "CA"]))

    def test_atoms_are_renumbered_consecutively(self):
        self.assertEqual([a.nr for a in self.frame.atoms[:3]], ["1", "2", "3"])
        self.assertEqual(self.frame.atoms[-1].nr, str(len(self.frame) % 100000))

    def test_title_atom_count_and_box_are_preserved(self):
        self.assertEqual(len(self.frame), len(self.source))
        self.assertEqual(self.frame.box, self.source.box)
        self.assertEqual(self.frame.title, self.source.title)


class ApplyTopologiesFailureTest(unittest.TestCase):
    def test_a_missing_alias_names_both_sides(self):
        """The error must say which aliases the file needs."""
        top = ResidueTopology.from_itp(ItpFile.read(example("TP3.itp")),
                                       AliasTable(), gro_resname="TIP3")
        with self.assertRaises(TopologyError) as caught:
            apply_topologies(GroFrame.read(GRO), [top])
        message = str(caught.exception)
        self.assertIn("only in the topology: O", message)
        self.assertIn("only in the coordinates: OH2", message)

    def test_a_residue_with_no_topology_passes_through(self):
        """Everything the topologies do not cover must be left exactly alone."""
        table = AliasTable().load(MAP)
        partial = [top for top in load_system(table) if top.itp_resname != "CA"]
        frame = GroFrame.read(GRO)
        messages = []
        apply_topologies(frame, partial, report=messages.append)
        cal = [a for a in frame if a.resname == "CAL"]
        self.assertEqual([a.name for a in cal], ["CAL", "CAL"])
        self.assertTrue(any("CAL" in m for m in messages))

    def test_strict_rejects_an_uncovered_residue(self):
        table = AliasTable().load(MAP)
        partial = [top for top in load_system(table) if top.itp_resname != "CA"]
        with self.assertRaises(TopologyError) as caught:
            apply_topologies(GroFrame.read(GRO), partial, strict=True)
        self.assertIn("CAL", str(caught.exception))

    def test_two_topologies_cannot_claim_one_residue(self):
        table = AliasTable()
        tops = [ResidueTopology.from_itp(ItpFile.read(example("Na+.itp")),
                                         table, gro_resname="SOD")] * 2
        with self.assertRaises(TopologyError):
            apply_topologies(GroFrame.read(GRO), tops)

    def test_no_match_at_all_is_an_error(self):
        top = ResidueTopology.from_itp(ItpFile.read(example("Na+.itp")),
                                       AliasTable(), gro_resname="NOPE")
        with self.assertRaises(TopologyError):
            apply_topologies(GroFrame.read(GRO), [top])

    def test_keep_names_reorders_without_renaming(self):
        table = AliasTable().load(MAP)
        frame = GroFrame.read(GRO)
        apply_topologies(frame, load_system(table), rename=False)
        self.assertIn("TIP3", frame.resnames)
        water = [a for a in frame if a.resname == "TIP3"][:3]
        # re-ordered into the topology's order, but still under the .gro's names
        self.assertEqual([a.name for a in water], ["OH2", "H1", "H2"])


def _tmp():
    handle, path = tempfile.mkstemp(suffix=".gro")
    os.close(handle)
    return path


if __name__ == "__main__":
    unittest.main()
