"""exact_test checked against hand counts and independent enumeration."""

from __future__ import annotations

import unittest
from fractions import Fraction

from mwexact import exact_test
from mwexact.exact import parse_values, u_statistic
from tests.reference import reference_p, reference_u


class HandComputedTests(unittest.TestCase):
    def test_tied_pairs_u_and_p(self):
        # A={1,2}, B={1,3}: U = 1/2+0+1+0 = 3/2; all 6 relabellings land
        # at distance >= 1/2 from the centre, so p = 1.
        result = exact_test(["1", "2"], ["1", "3"])
        self.assertEqual(result.u, Fraction(3, 2))
        self.assertEqual(result.p, Fraction(1))
        self.assertEqual((result.n, result.m), (2, 2))

    def test_separated_single_observation(self):
        # A={3}, B={1,2}: U = 2; relabellings give distances 1, 1, 0,
        # so p = 2/3.
        result = exact_test(["3"], ["1", "2"])
        self.assertEqual(result.u, Fraction(2))
        self.assertEqual(result.p, Fraction(2, 3))

    def test_decimal_strings_parsed_exactly(self):
        a = parse_values(["0.1", "0.2"], "a")
        b = parse_values(["0.15"], "b")
        self.assertEqual(a[0], Fraction(1, 10))
        self.assertEqual(a[1], Fraction(1, 5))
        self.assertEqual(u_statistic(a, b), Fraction(1))

    def test_rejects_empty_group(self):
        with self.assertRaises(ValueError):
            exact_test([], ["1"])
        with self.assertRaises(ValueError):
            exact_test(["1"], [])

    def test_rejects_bad_decimal(self):
        with self.assertRaises(ValueError):
            exact_test(["1", "abc"], ["2"])


class EnumerationAgreementTests(unittest.TestCase):
    DATASETS = [
        (["1", "2", "2"], ["2", "3"]),          # ties, unequal sizes
        (["0.5", "0.5", "1.5"], ["0.5", "2.5"]),  # duplicate decimals
        (["-1.25", "0", "0", "3.5"], ["0", "-1.25"]),
        (["2", "2", "2"], ["2", "2"]),          # everything tied
        (["10", "20", "30", "40"], ["5", "15"]),  # separated, unequal
        (["0.1", "0.2"], ["0.15", "0.25", "0.05"]),
    ]

    def test_u_matches_reference(self):
        for a_raw, b_raw in self.DATASETS:
            with self.subTest(a=a_raw, b=b_raw):
                a = parse_values(a_raw, "a")
                b = parse_values(b_raw, "b")
                self.assertEqual(
                    exact_test(a_raw, b_raw).u, reference_u(a, b)
                )

    def test_p_matches_reference(self):
        for a_raw, b_raw in self.DATASETS:
            with self.subTest(a=a_raw, b=b_raw):
                a = parse_values(a_raw, "a")
                b = parse_values(b_raw, "b")
                self.assertEqual(
                    exact_test(a_raw, b_raw).p, reference_p(a, b)
                )


if __name__ == "__main__":
    unittest.main()
