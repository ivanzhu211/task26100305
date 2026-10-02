"""Tests for mannwhitney: reference enumeration and confidence-set checks.

The reference implementation below is an independent brute-force label
enumeration written straight from the problem statement; the library is
checked against it on small samples with ties and unequal group sizes.
"""

import itertools
import unittest
from fractions import Fraction

from mannwhitney import exact_test, invert_shift


# --------------------------------------------------------------------------
# Independent reference implementation (label enumeration over indices).
# --------------------------------------------------------------------------

def ref_u(a, b):
    u = Fraction(0)
    for x in a:
        for y in b:
            if x > y:
                u += 1
            elif x == y:
                u += Fraction(1, 2)
    return u


def ref_p(a, b):
    n, m = len(a), len(b)
    values = list(a) + list(b)
    center = Fraction(n * m, 2)
    d_obs = abs(ref_u(a, b) - center)
    extreme = 0
    total = 0
    for idx in itertools.combinations(range(n + m), n):
        keep = frozenset(idx)
        ga = [v for i, v in enumerate(values) if i in keep]
        gb = [v for i, v in enumerate(values) if i not in keep]
        if abs(ref_u(ga, gb) - center) >= d_obs:
            extreme += 1
        total += 1
    return Fraction(extreme, total)


def parse(group):
    return [Fraction(s) for s in group]


def ref_member(a, b, theta, alpha):
    shifted = [x - theta for x in a]
    return ref_p(shifted, b) > alpha


# Samples with duplicated values and unequal group sizes.
SAMPLE_TIES_A = ["1.5", "2", "2", "3.25"]
SAMPLE_TIES_B = ["2", "2.5"]

SAMPLE_HEAVY_TIES_A = ["-0.5", "0", "0", "1"]
SAMPLE_HEAVY_TIES_B = ["0", "0.5", "1", "1"]


class ExactTestTests(unittest.TestCase):
    def test_hand_computed_singleton(self):
        # A=['1'], B=['2']: U = 0; both label assignments give |U - 1/2| = 1/2.
        r = exact_test(["1"], ["2"])
        self.assertEqual(r.u, Fraction(0))
        self.assertEqual(r.p_value, Fraction(1))

    def test_hand_computed_separated(self):
        # A=['3','4','5'], B=['0','1','2']: U = 9, maximal deviation.
        # Only 2 of C(6,3)=20 assignments reach deviation 9/2.
        r = exact_test(["3", "4", "5"], ["0", "1", "2"])
        self.assertEqual(r.u, Fraction(9))
        self.assertEqual(r.p_value, Fraction(2, 20))

    def test_ties_unequal_sizes_against_reference(self):
        r = exact_test(SAMPLE_TIES_A, SAMPLE_TIES_B)
        self.assertEqual(r.u, ref_u(parse(SAMPLE_TIES_A), parse(SAMPLE_TIES_B)))
        self.assertEqual(r.p_value, ref_p(parse(SAMPLE_TIES_A), parse(SAMPLE_TIES_B)))

    def test_heavy_ties_against_reference(self):
        r = exact_test(SAMPLE_HEAVY_TIES_A, SAMPLE_HEAVY_TIES_B)
        self.assertEqual(
            r.u, ref_u(parse(SAMPLE_HEAVY_TIES_A), parse(SAMPLE_HEAVY_TIES_B))
        )
        self.assertEqual(
            r.p_value, ref_p(parse(SAMPLE_HEAVY_TIES_A), parse(SAMPLE_HEAVY_TIES_B))
        )

    def test_u_is_exact_half_integer(self):
        # 1.5 vs 2 (loss), 1.5 vs 2.5 (loss), 2 vs 2 (tie), 2 vs 2.5 (loss),
        # 2 vs 2 (tie), 2 vs 2.5 (loss), 3.25 vs 2 (win), 3.25 vs 2.5 (win)
        r = exact_test(SAMPLE_TIES_A, SAMPLE_TIES_B)
        self.assertEqual(r.u, Fraction(3))

    def test_invalid_input_rejected(self):
        with self.assertRaises(ValueError):
            exact_test([], ["1"])
        with self.assertRaises(ValueError):
            exact_test(["1"], [])
        with self.assertRaises(ValueError):
            exact_test(["abc"], ["1"])


class InversionTests(unittest.TestCase):
    def candidate_thetas(self, a, b):
        av, bv = parse(a), parse(b)
        diffs = sorted({x - y for x in av for y in bv})
        thetas = list(diffs)  # cross differences themselves
        for lo, hi in zip(diffs, diffs[1:]):
            thetas.append((lo + hi) / 2)  # between adjacent differences
        thetas.append(diffs[0] - 1)  # outside on the left
        thetas.append(diffs[-1] + 1)  # outside on the right
        thetas.append(diffs[0] - Fraction(1, 3))
        thetas.append(diffs[-1] + Fraction(1, 3))
        return thetas

    def check_membership(self, a, b, confidence, alpha):
        cs = invert_shift(a, b, confidence)
        av, bv = parse(a), parse(b)
        for theta in self.candidate_thetas(a, b):
            expected = ref_member(av, bv, theta, alpha)
            actual = theta in cs
            self.assertEqual(
                actual,
                expected,
                f"theta={theta} ({float(theta):.4f}): set says {actual}, "
                f"reference says {expected}",
            )
        return cs

    def test_membership_ties_unequal_90(self):
        self.check_membership(SAMPLE_TIES_A, SAMPLE_TIES_B, 0.90, Fraction(1, 10))

    def test_membership_ties_unequal_95(self):
        self.check_membership(SAMPLE_TIES_A, SAMPLE_TIES_B, "95%", Fraction(1, 20))

    def test_membership_heavy_ties_99(self):
        self.check_membership(
            SAMPLE_HEAVY_TIES_A, SAMPLE_HEAVY_TIES_B, 99, Fraction(1, 100)
        )

    def test_p_equal_alpha_is_excluded(self):
        # Fully separated samples: at theta outside [1, 5] the shifted groups
        # stay fully separated, giving p = 2/20 = 0.1 exactly.  At 90%
        # confidence (alpha = 0.1) those theta must NOT be accepted.
        a, b = ["3", "4", "5"], ["0", "1", "2"]
        cs = invert_shift(a, b, 0.90)
        av, bv = parse(a), parse(b)

        boundary_theta = Fraction(1, 2)  # shifted groups fully separated
        p_here = ref_p([x - boundary_theta for x in av], bv)
        self.assertEqual(p_here, Fraction(1, 10))  # p == alpha exactly
        self.assertNotIn(boundary_theta, cs)

        # The endpoints theta = 1 and theta = 5 have p > alpha: included.
        self.assertIn(Fraction(1), cs)
        self.assertIn(Fraction(5), cs)
        # Just outside the endpoints the set stops (open exterior).
        self.assertNotIn(Fraction(999, 1000), cs)
        self.assertNotIn(Fraction(5001, 1000), cs)

    def test_confidence_spellings(self):
        a, b = SAMPLE_TIES_A, SAMPLE_TIES_B
        base = invert_shift(a, b, 0.95)
        for spelling in ("95%", "0.95", 95, Fraction(19, 20)):
            self.assertEqual(invert_shift(a, b, spelling).intervals, base.intervals)

    def test_unsupported_confidence_rejected(self):
        with self.assertRaises(ValueError):
            invert_shift(SAMPLE_TIES_A, SAMPLE_TIES_B, 0.80)

    def test_degenerate_single_point_set(self):
        # Constant groups: only theta = -1 makes all pairs tie (U = 9/2,
        # deviation 0, p = 1); every other theta keeps the groups fully
        # separated (p = 0.1 == alpha).  The confidence set is exactly {-1}.
        cs = invert_shift(["0", "0", "0"], ["1", "1", "1"], 0.90)
        self.assertEqual(len(cs.intervals), 1)
        iv = cs.intervals[0]
        self.assertEqual(iv.lower, Fraction(-1))
        self.assertEqual(iv.upper, Fraction(-1))
        self.assertTrue(iv.lower_closed)
        self.assertTrue(iv.upper_closed)
        self.assertIn(Fraction(-1), cs)
        self.assertNotIn(Fraction(-1) + Fraction(1, 1000), cs)
        self.assertNotIn(Fraction(-2), cs)


if __name__ == "__main__":
    unittest.main()
