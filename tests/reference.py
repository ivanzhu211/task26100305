"""Independent brute-force reference used only by the tests.

Enumerates every size-preserving relabelling of the pooled sample and
scores U straight from the pair-counting definition.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Sequence


def reference_u(a: Sequence[Fraction], b: Sequence[Fraction]) -> Fraction:
    total = Fraction(0)
    for x in a:
        for y in b:
            if x > y:
                total += 1
            elif x == y:
                total += Fraction(1, 2)
    return total


def reference_p(a: Sequence[Fraction], b: Sequence[Fraction]) -> Fraction:
    n, m = len(a), len(b)
    pooled = list(a) + list(b)
    centre = Fraction(n * m, 2)
    observed = abs(reference_u(a, b) - centre)
    extreme = 0
    total = 0
    for keep in combinations(range(n + m), n):
        labels = ["B"] * (n + m)
        for index in keep:
            labels[index] = "A"
        new_a = [value for value, label in zip(pooled, labels) if label == "A"]
        new_b = [value for value, label in zip(pooled, labels) if label == "B"]
        if abs(reference_u(new_a, new_b) - centre) >= observed:
            extreme += 1
        total += 1
    return Fraction(extreme, total)


def reference_shift_p(
    a: Sequence[Fraction], b: Sequence[Fraction], theta: Fraction
) -> Fraction:
    return reference_p([x - theta for x in a], list(b))
