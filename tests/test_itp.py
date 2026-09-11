"""Topology parsing, re-ordering and the index map."""

import unittest

from _support import example, read
from topology_wrangling import ItpFile, AtomOrder, TopologyError
from topology_wrangling.chem import element_of, is_hydrogen
from topology_wrangling.ordering import update_group_order

SLIPIDS = example("update_group_incompatible.itp")
GOLDEN = example("DOPC_reordered.itp")


class ParsingTest(unittest.TestCase):
    def test_round_trip_is_byte_identical(self):
        """An unmodified topology must come back out exactly as it went in."""
        for name in ("update_group_incompatible.itp", "TP3.itp", "Na+.itp",
                     "Cl-.itp", "MG.itp", "CA.itp"):
            with self.subTest(name):
                path = example(name)
                self.assertEqual(ItpFile.read(path).text(), read(path))

    def test_atoms(self):
        itp = ItpFile.read(SLIPIDS)
        self.assertEqual(itp.molecule_name, "DOPC")
        self.assertEqual(itp.natoms, 138)
        self.assertEqual(itp.residue_names, ["DOPC"])
        self.assertEqual(itp.residue_name(), "DOPC")
        first = itp.atoms[0]
        self.assertEqual((first.nr, first.name, first.type), (1, "N", "NTL"))
        self.assertAlmostEqual(first.mass, 14.007)

    def test_bond_graph(self):
        itp = ItpFile.read(SLIPIDS)
        self.assertEqual(len(itp.bonds()), 137)
        graph = itp.bond_graph()
        # Atom 1 is the choline nitrogen, bonded to four carbons.
        self.assertEqual(sorted(graph[1]), [2, 3, 4, 5])

    def test_bad_input_raises_topology_error(self):
        itp = ItpFile.from_text("[ atoms ]\n 1 C 1 RES C 1\n[ bonds ]\n 1 9 1\n")
        with self.assertRaises(TopologyError):
            itp.bond_graph()
        with self.assertRaises(TopologyError):
            ItpFile.from_text("[ moleculetype ]\nX 3\n").atoms

    def test_several_residues_must_be_disambiguated(self):
        itp = ItpFile.from_text(
            "[ atoms ]\n 1 C 1 AAA C 1\n 2 C 1 BBB C 1\n"
        )
        with self.assertRaises(TopologyError):
            itp.residue_name()
        self.assertEqual(itp.residue_name("BBB"), "BBB")


class ChemTest(unittest.TestCase):
    def test_hydrogens_are_recognised(self):
        atoms = ItpFile.read(SLIPIDS).atoms
        self.assertEqual(sum(is_hydrogen(a) for a in atoms), 84)

    def test_element_from_mass_and_from_name(self):
        atoms = {a.name: a for a in ItpFile.read(SLIPIDS).atoms}
        self.assertEqual(element_of(atoms["N"]), "N")
        self.assertEqual(element_of(atoms["C12"]), "C")
        self.assertEqual(element_of(atoms["H12A"]), "H")
        self.assertEqual(element_of(atoms["P"]), "P")
        # The ion topologies exercise the two-letter symbols.
        self.assertEqual(element_of(ItpFile.read(example("Cl-.itp")).atoms[0]), "Cl")
        self.assertEqual(element_of(ItpFile.read(example("MG.itp")).atoms[0]), "Mg")
        self.assertEqual(element_of(ItpFile.read(example("CA.itp")).atoms[0]), "Ca")


class UpdateGroupOrderTest(unittest.TestCase):
    def setUp(self):
        self.itp = ItpFile.read(SLIPIDS)
        self.order, self.hydrogen = update_group_order(
            self.itp.atoms, self.itp.bond_graph(), "DOPC"
        )

    def test_every_hydrogen_follows_its_heavy_atom(self):
        """The property the whole tool exists to establish."""
        graph = self.itp.bond_graph()
        names = self.order.new_names
        heavy_at = None
        for position, old in enumerate(self.order.old_for_new):
            if not self.hydrogen[old - 1]:
                heavy_at = old
                continue
            heavies = [n for n in graph.get(old, ()) if not self.hydrogen[n - 1]]
            if not heavies:
                continue    # an unbonded hydrogen keeps its own slot
            self.assertIn(
                heavy_at, heavies,
                "hydrogen %s at position %d does not follow a heavy atom it is "
                "bonded to" % (names[position], position + 1),
            )

    def test_heavy_atoms_keep_their_relative_order(self):
        heavy = [old for old in self.order.old_for_new
                 if not self.hydrogen[old - 1]]
        self.assertEqual(heavy, sorted(heavy))

    def test_permutation_is_complete(self):
        self.assertEqual(sorted(self.order.old_for_new), list(range(1, 139)))
        self.assertEqual(self.order.moved(), 11)

    def test_interactions_survive_the_renumbering(self):
        """Every bonded term must still join the same atoms."""
        new = self.itp.permute(self.order)
        mapping = self.order.new_for_old
        for directive, ncols in (("bonds", 2), ("angles", 3), ("dihedrals", 4),
                                 ("pairs", 2)):
            with self.subTest(directive):
                before = sorted(
                    tuple(mapping[int(t)] for t in tokens[:ncols])
                    for tokens, _ in self.itp.rows(directive)
                )
                after = sorted(
                    tuple(int(t) for t in tokens[:ncols])
                    for tokens, _ in new.rows(directive)
                )
                self.assertEqual(before, after)
                self.assertTrue(before)

    def test_atom_order_inside_an_interaction_is_untouched(self):
        """Impropers are order-sensitive, so rows must not be re-arranged."""
        new = self.itp.permute(self.order)
        mapping = self.order.new_for_old
        before = {tuple(mapping[int(t)] for t in tokens[:4])
                  for tokens, _ in self.itp.rows("dihedrals")
                  if tokens[4] == "2"}
        after = {tuple(int(t) for t in tokens[:4])
                 for tokens, _ in new.rows("dihedrals") if tokens[4] == "2"}
        self.assertEqual(before, after)

    def test_reordered_output_matches_the_reference(self):
        self.assertEqual(self.itp.permute(self.order).text(), read(GOLDEN))

    def test_reordering_is_idempotent(self):
        new = self.itp.permute(self.order)
        again, _ = update_group_order(new.atoms, new.bond_graph())
        self.assertEqual(again.moved(), 0)

    def test_charge_group_modes(self):
        per_atom = self.itp.permute(self.order, cgnr_mode="per-atom")
        self.assertEqual([a.cgnr for a in per_atom.atoms],
                         [str(i) for i in range(1, 139)])
        kept = self.itp.permute(self.order, cgnr_mode="keep")
        self.assertEqual([a.cgnr for a in kept.atoms],
                         [self.itp.atoms[old - 1].cgnr
                          for old in self.order.old_for_new])
        with self.assertRaises(TopologyError):
            self.itp.permute(self.order, cgnr_mode="nonsense")


class AtomOrderTest(unittest.TestCase):
    def test_map_file_round_trip(self):
        import tempfile, os
        itp = ItpFile.read(SLIPIDS)
        order, _ = update_group_order(itp.atoms, itp.bond_graph(), "DOPC")
        handle, path = tempfile.mkstemp(suffix=".txt")
        os.close(handle)
        try:
            order.write_map(path)
            back = AtomOrder.read_map(path)
        finally:
            os.remove(path)
        self.assertEqual(back.old_for_new, order.old_for_new)
        self.assertEqual(back.names, order.names)
        self.assertEqual(back.resname, "DOPC")

    def test_inverse_undoes_the_order(self):
        itp = ItpFile.read(SLIPIDS)
        order, _ = update_group_order(itp.atoms, itp.bond_graph())
        names = itp.atom_names
        self.assertEqual(order.inverse().apply(order.apply(names)), names)

    def test_a_non_permutation_is_rejected(self):
        with self.assertRaises(TopologyError):
            AtomOrder([1, 2, 2])


if __name__ == "__main__":
    unittest.main()
