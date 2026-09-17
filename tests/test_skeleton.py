"""The ChemDraw-style skeletal layout: chain angles, branches, rings."""

import math
import unittest

from _support import example
from topology_wrangling import ItpFile, TopologyError
from topology_wrangling.chem import is_hydrogen
from topology_wrangling.skeleton import (BOND, find_rings, skeletal,
                                         _count_crossings, _score)

# toluene: a six-ring carrying a methyl, enough to exercise rings and branches
TOLUENE_BONDS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 1), (6, 7)]
# naphthalene: two six-rings sharing one bond
NAPHTHALENE_BONDS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8),
                     (8, 9), (9, 10), (10, 5), (10, 1)]


def angle_between(a, b, c):
    """The interior angle at b, in degrees."""
    v1 = (a[0] - b[0], a[1] - b[1])
    v2 = (c[0] - b[0], c[1] - b[1])
    dot = v1[0] * v2[0] + v1[1] * v2[1]
    size = math.hypot(*v1) * math.hypot(*v2)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / size))))


class ChainTest(unittest.TestCase):
    def setUp(self):
        self.bonds = [(i, i + 1) for i in range(1, 12)]
        self.points = skeletal(12, self.bonds)

    def test_all_bonds_are_the_same_length(self):
        """Equal bond lengths are what make the drawing read as a structure."""
        lengths = [math.dist(self.points[i - 1], self.points[j - 1])
                   for i, j in self.bonds]
        for length in lengths:
            self.assertAlmostEqual(length, lengths[0], places=6)

    def test_a_chain_zigzags_at_120_degrees(self):
        for i in range(2, 12):
            with self.subTest(atom=i):
                self.assertAlmostEqual(
                    angle_between(self.points[i - 2], self.points[i - 1],
                                  self.points[i]),
                    120.0, places=4)

    def test_a_chain_runs_along_one_axis(self):
        """The zigzag must advance, not curl: a 12-carbon chain is long."""
        xs = [p[0] for p in self.points]
        ys = [p[1] for p in self.points]
        self.assertGreater((max(xs) - min(xs)) / (max(ys) - min(ys)), 5.0)

    def test_no_two_atoms_land_on_each_other(self):
        # coordinates are normalised, so measure against this drawing's bond
        bond = math.dist(self.points[0], self.points[1])
        for i in range(len(self.points)):
            for j in range(i + 1, len(self.points)):
                with self.subTest(pair=(i + 1, j + 1)):
                    self.assertGreater(
                        math.dist(self.points[i], self.points[j]), bond * 0.95)


class BranchTest(unittest.TestCase):
    def test_three_bonds_fan_120_degrees_apart(self):
        # a central atom with three chains hanging off it
        bonds = [(1, 2), (1, 3), (1, 4), (2, 5), (3, 6), (4, 7)]
        p = skeletal(7, bonds)
        for a, b in ((2, 3), (3, 4), (2, 4)):
            with self.subTest(pair=(a, b)):
                self.assertAlmostEqual(
                    angle_between(p[a - 1], p[0], p[b - 1]), 120.0, places=3)

    def test_four_bonds_fan_90_degrees_apart(self):
        bonds = [(1, 2), (1, 3), (1, 4), (1, 5),
                 (2, 6), (3, 7), (4, 8), (5, 9)]
        p = skeletal(9, bonds)
        angles = sorted(angle_between(p[n - 1], p[0], p[m - 1])
                        for n, m in ((2, 3), (3, 4), (4, 5), (2, 5)))
        for value in angles[:2]:
            self.assertAlmostEqual(value, 90.0, places=3)


class HydrogenTest(unittest.TestCase):
    def test_hydrogens_hang_off_and_do_not_join_the_chain(self):
        # C-C-C with two hydrogens on the middle carbon
        bonds = [(1, 2), (2, 3), (2, 4), (2, 5)]
        p = skeletal(5, bonds, pendant={4, 5})
        # the carbon chain keeps its 120 degrees despite the hydrogens
        self.assertAlmostEqual(angle_between(p[0], p[1], p[2]), 120.0, places=3)
        # and the hydrogens sit closer in than a skeleton bond
        for h in (4, 5):
            self.assertLess(math.dist(p[1], p[h - 1]),
                            math.dist(p[0], p[1]) * 0.95)

    def test_hydrogens_avoid_the_directions_the_skeleton_uses(self):
        bonds = [(1, 2), (2, 3), (2, 4), (2, 5)]
        p = skeletal(5, bonds, pendant={4, 5})
        for h in (4, 5):
            for heavy in (1, 3):
                with self.subTest(h=h, heavy=heavy):
                    self.assertGreater(
                        angle_between(p[h - 1], p[1], p[heavy - 1]), 30.0)


class RingTest(unittest.TestCase):
    def test_rings_are_found(self):
        adjacency = {}
        for i, j in TOLUENE_BONDS:
            adjacency.setdefault(i, []).append(j)
            adjacency.setdefault(j, []).append(i)
        rings = find_rings(adjacency)
        self.assertEqual(len(rings), 1)
        self.assertEqual(set(rings[0]), {1, 2, 3, 4, 5, 6})

    def test_a_six_ring_is_a_regular_hexagon(self):
        p = skeletal(7, TOLUENE_BONDS)
        ring = [p[i - 1] for i in (1, 2, 3, 4, 5, 6)]
        centre = (sum(x for x, _ in ring) / 6, sum(y for _, y in ring) / 6)
        radii = [math.dist(centre, point) for point in ring]
        for radius in radii:
            self.assertAlmostEqual(radius, radii[0], places=6)
        for i in range(6):
            self.assertAlmostEqual(
                angle_between(ring[i - 1], ring[i], ring[(i + 1) % 6]),
                120.0, places=4)

    def test_ring_bonds_are_the_normal_bond_length(self):
        p = skeletal(7, TOLUENE_BONDS)
        for i, j in TOLUENE_BONDS:
            with self.subTest(bond=(i, j)):
                self.assertAlmostEqual(math.dist(p[i - 1], p[j - 1]),
                                       math.dist(p[0], p[1]), places=6)

    def test_fused_rings_share_their_bond_and_do_not_overlap(self):
        p = skeletal(10, NAPHTHALENE_BONDS)
        # every atom distinct
        bond = math.dist(p[0], p[1])
        for i in range(10):
            for j in range(i + 1, 10):
                with self.subTest(pair=(i + 1, j + 1)):
                    self.assertGreater(math.dist(p[i], p[j]), bond * 0.95)
        # both rings still regular hexagons
        for ring in ((1, 2, 3, 4, 5, 10), (5, 6, 7, 8, 9, 10)):
            points = [p[i - 1] for i in ring]
            centre = (sum(x for x, _ in points) / 6,
                      sum(y for _, y in points) / 6)
            radii = [math.dist(centre, q) for q in points]
            with self.subTest(ring=ring):
                for radius in radii:
                    self.assertAlmostEqual(radius, radii[0], places=5)

    def test_a_ring_drawing_has_no_crossing_bonds(self):
        p = skeletal(10, NAPHTHALENE_BONDS)
        self.assertEqual(_count_crossings(p, NAPHTHALENE_BONDS), 0)


class RealMoleculeTest(unittest.TestCase):
    def setUp(self):
        self.itp = ItpFile.read(example("DOPC_reordered.itp"))
        self.hydrogens = {i for i, a in enumerate(self.itp.atoms, start=1)
                          if is_hydrogen(a)}
        self.bonds = self.itp.connectivity()

    def test_a_lipid_is_laid_out_without_crossings(self):
        points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        heavy = [(i, j) for i, j in self.bonds
                 if i not in self.hydrogens and j not in self.hydrogens]
        self.assertEqual(_count_crossings(points, heavy), 0)

    def test_it_is_far_more_compact_than_a_spring_layout(self):
        """The complaint that started this: the spring layout sprawled."""
        from topology_wrangling.layout import force_directed
        skeleton_points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        spring_points = force_directed(self.itp.natoms, self.bonds, seed=0)

        def bond_to_extent(points):
            lengths = sorted(math.dist(points[i - 1], points[j - 1])
                             for i, j in self.bonds)
            typical = lengths[len(lengths) // 2]
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            return max(max(xs) - min(xs), max(ys) - min(ys)) / typical

        self.assertLess(bond_to_extent(skeleton_points),
                        bond_to_extent(spring_points))

    def test_no_branch_folds_back_onto_another(self):
        """The reported defect: the head group looped back over the glycerol.

        Nothing at the moment of placing the P-O12 bond can know it will
        clash, so the layout has to detect the clash afterwards and flip the
        branch.  This asserts the result: no two non-bonded heavy atoms end up
        closer than a bond length.
        """
        points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        heavy = [i for i in range(1, self.itp.natoms + 1)
                 if i not in self.hydrogens]
        bonded = {(min(i, j), max(i, j)) for i, j in self.bonds}
        unit = sorted(math.dist(points[i - 1], points[j - 1])
                      for i, j in self.bonds
                      if i not in self.hydrogens and j not in self.hydrogens)
        unit = unit[len(unit) // 2]
        for a in heavy:
            for b in heavy:
                if b <= a or (a, b) in bonded:
                    continue
                gap = math.dist(points[a - 1], points[b - 1])
                self.assertGreater(
                    gap, unit * 0.9,
                    "atoms %s and %s sit %.2f apart, under one bond"
                    % (self.itp.atoms[a - 1].name,
                       self.itp.atoms[b - 1].name, gap))

    def test_flipping_a_branch_keeps_every_bond_length(self):
        """A reflection is a rigid motion, so unfolding must not stretch bonds."""
        points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        heavy = [math.dist(points[i - 1], points[j - 1])
                 for i, j in self.bonds
                 if i not in self.hydrogens and j not in self.hydrogens]
        for length in heavy:
            self.assertAlmostEqual(length, heavy[0], places=6)

    def test_the_drawing_does_not_waste_its_page(self):
        """Whitespace is the second thing the layout optimises, after clashes.

        A splayed lipid took about 8 square bond lengths per atom; packing the
        two tails roughly parallel brings that near 1.  Anything much above 2
        means the drawing has sprawled again.
        """
        points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        heavy = [(i, j) for i, j in self.bonds
                 if i not in self.hydrogens and j not in self.hydrogens]
        lengths = sorted(math.dist(points[i - 1], points[j - 1])
                         for i, j in heavy)
        unit = lengths[len(lengths) // 2]
        xs = [p[0] / unit for p in points]
        ys = [p[1] / unit for p in points]
        area = (max(xs) - min(xs)) * (max(ys) - min(ys))
        self.assertLess(area / self.itp.natoms, 2.0)

    def test_the_two_tails_end_up_roughly_parallel(self):
        """What packing a lipid tightly actually means."""
        points = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        names = {a.name: i for i, a in enumerate(self.itp.atoms, start=1)}

        def axis(first, last):
            a = points[names[first] - 1]
            b = points[names[last] - 1]
            return math.atan2(b[1] - a[1], b[0] - a[0])

        sn1 = axis("C32", "C318")
        sn2 = axis("C22", "C218")
        between = abs(math.degrees(_wrap(sn1 - sn2)))
        self.assertLess(min(between, 180.0 - between), 25.0,
                        "the tails run %.0f degrees apart" % between)

    def test_it_is_deterministic(self):
        first = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        second = skeletal(self.itp.natoms, self.bonds, self.hydrogens)
        self.assertEqual(first, second)

    def test_branch_strategies_differ_and_both_are_usable(self):
        """Which one is squarer is not fixed; that both are clean is."""
        heavy = [(i, j) for i, j in self.bonds
                 if i not in self.hydrogens and j not in self.hydrogens]
        layouts = {}
        for strategy in ("spread", "compact"):
            layouts[strategy] = skeletal(self.itp.natoms, self.bonds,
                                         self.hydrogens, branch=strategy)
            with self.subTest(strategy):
                self.assertEqual(
                    _count_crossings(layouts[strategy], heavy), 0)
        self.assertNotEqual(layouts["spread"], layouts["compact"])

    def test_auto_is_no_worse_than_either_strategy(self):
        chosen = skeletal(self.itp.natoms, self.bonds, self.hydrogens,
                          branch="auto")
        scores = [_score(skeletal(self.itp.natoms, self.bonds, self.hydrogens,
                                  branch=s), self.bonds)
                  for s in ("spread", "compact")]
        self.assertLessEqual(_score(chosen, self.bonds), min(scores) + 1e-9)

    def test_an_unknown_strategy_is_rejected(self):
        with self.assertRaises(TopologyError):
            skeletal(4, [(1, 2)], branch="sideways")


class DegenerateTest(unittest.TestCase):
    def test_a_single_atom(self):
        self.assertEqual(len(skeletal(1, [])), 1)

    def test_an_empty_molecule_is_rejected(self):
        with self.assertRaises(TopologyError):
            skeletal(0, [])

    def test_disconnected_molecules_are_placed_side_by_side(self):
        # two separate two-atom molecules, neither sitting on the other
        points = skeletal(4, [(1, 2), (3, 4)])
        bond = math.dist(points[0], points[1])
        for a, b in ((0, 2), (0, 3), (1, 2), (1, 3)):
            with self.subTest(pair=(a, b)):
                self.assertGreater(math.dist(points[a], points[b]), bond * 0.9)

    def test_a_bond_outside_the_molecule_is_rejected(self):
        with self.assertRaises(TopologyError):
            skeletal(3, [(1, 9)])

    def test_water_from_settles(self):
        itp = ItpFile.read(example("TP3.itp"))
        points = skeletal(itp.natoms, itp.connectivity(), {2, 3})
        self.assertEqual(len(points), 3)


def _wrap(angle):
    return (angle + math.pi) % (2 * math.pi) - math.pi


def _aspect(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    width = max(max(xs) - min(xs), 1e-9)
    height = max(max(ys) - min(ys), 1e-9)
    return max(width / height, height / width)


if __name__ == "__main__":
    unittest.main()
