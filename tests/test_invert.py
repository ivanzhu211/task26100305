"""invert_shift checked against independent enumeration of p(theta)."""

from __future__ import annotations

import unittest
from fractions import Fraction

from mwexact import invert_shift
from mwexact.exact import parse_values
from tests.reference import reference_shift_p


def candidate_thetas(a, b):
    """Cross differences, gap midpoints, and points beyond both ends."""
    differences = sorted({x - y for x in a for y in b})
    points = [differences[0] - 2, differences[0] - 1]
    for index, point in enumerate(differences):
        points.append(point)
        if index + 1 < len(differences):
            points.append((point + differences[index + 1]) / 2)
    points += [differences[-1] + 1, differences[-1] + 2]
    return points


class MembershipAgainstReferenceTests(unittest.TestCase):
    DATASETS = [
        (["1", "2", "4", "5"], ["0", "1", "3"]),
        (["0.5", "0.5", "1.5"], ["0.5", "2.5"]),
        (["-1.25", "0", "3.5"], ["0", "-1.25", "2"]),
        (["2", "2", "2"], ["2", "2"]),
        (["0.1", "0.2"], ["0.15", "0.25", "0.05"]),
    ]
    LEVELS = [0.9, 0.95, 0.99]

    def test_membership_matches_reference(self):
        for a_raw, b_raw in self.DATASETS:
            a = parse_values(a_raw, "a")
            b = parse_values(b_raw, "b")
            for level in self.LEVELS:
                alpha = 1 - Fraction(str(level))
                cs = invert_shift(a_raw, b_raw, level)
                for theta in candidate_thetas(a, b):
                    with self.subTest(a=a_raw, b=b_raw, level=level, theta=theta):
                        expected = reference_shift_p(a, b, theta) > alpha
                        self.assertEqual(cs.contains(theta), expected)


class BoundaryEqualityTests(unittest.TestCase):
    def test_unbounded_whole_line(self):
        # n = m = 1: every relabelling is equally extreme, so p(theta) = 1
        # everywhere and the confidence set is the whole real line.
        cs = invert_shift(["1"], ["2"], 0.9)
        self.assertEqual(len(cs.intervals), 1)
        only = cs.intervals[0]
        self.assertIsNone(only.low)
        self.assertIsNone(only.high)
        self.assertTrue(cs.contains("-1000000"))
        self.assertTrue(cs.contains("1000000"))

    def test_p_equal_to_alpha_is_excluded(self):
        # All values tied: p(theta) = 1/10 exactly for every theta != 0,
        # so at 90% the confidence set collapses to the single point {0}.
        cs = invert_shift(["0", "0"], ["0", "0", "0"], 0.9)
        self.assertEqual(len(cs.intervals), 1)
        only = cs.intervals[0]
        self.assertEqual(only.low, Fraction(0))
        self.assertEqual(only.high, Fraction(0))
        self.assertTrue(only.low_closed and only.high_closed)
        self.assertTrue(cs.contains("0"))
        self.assertFalse(cs.contains("1"))
        self.assertFalse(cs.contains("-1"))

    def test_exact_endpoint_membership(self):
        # Confirm via the reference that the excluded endpoints really do
        # sit exactly at p == alpha, not merely below it.
        a = parse_values(["0", "0"], "a")
        b = parse_values(["0", "0", "0"], "b")
        self.assertEqual(reference_shift_p(a, b, Fraction(1)), Fraction(1, 10))
        self.assertEqual(reference_shift_p(a, b, Fraction(-1)), Fraction(1, 10))
        self.assertEqual(reference_shift_p(a, b, Fraction(0)), Fraction(1))


class ConfidenceParsingTests(unittest.TestCase):
    A = ["1", "2", "4", "5"]
    B = ["0", "1", "3"]

    def test_equivalent_forms_agree(self):
        expected = invert_shift(self.A, self.B, 0.95)
        for form in ("95%", "0.95", Fraction(19, 20), 95):
            with self.subTest(form=form):
                self.assertEqual(invert_shift(self.A, self.B, form), expected)

    def test_unsupported_level_rejected(self):
        for bad in (0.5, 0.999, "80%", Fraction(4, 5)):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    invert_shift(self.A, self.B, bad)


if __name__ == "__main__":
    unittest.main()
