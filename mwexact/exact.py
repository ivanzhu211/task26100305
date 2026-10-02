"""Exact Mann-Whitney U test under full label reallocation."""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence


@dataclass(frozen=True)
class ExactTestResult:
    """Observed U statistic and its exact two-sided p-value."""

    u: Fraction
    p: Fraction
    n: int
    m: int


def parse_values(values: Sequence[str], name: str) -> tuple[Fraction, ...]:
    """Parse non-empty finite decimal strings into exact Fractions."""
    if not values:
        raise ValueError(f"{name} must be a non-empty sequence")
    parsed: list[Fraction] = []
    for raw in values:
        try:
            value = Fraction(raw)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} contains invalid decimal string: {raw!r}") from exc
        parsed.append(value)
    return tuple(parsed)


def u_statistic(a: Sequence[Fraction], b: Sequence[Fraction]) -> Fraction:
    """U = sum over cross pairs of 1 (a>b), 1/2 (a==b), 0 (a<b)."""
    total = Fraction(0)
    for x in a:
        for y in b:
            if x > y:
                total += 1
            elif x == y:
                total += Fraction(1, 2)
    return total


def exact_p_value(a: Sequence[Fraction], b: Sequence[Fraction]) -> Fraction:
    """Exact two-sided p-value over all size-preserving reallocations."""
    n = len(a)
    m = len(b)
    pooled = tuple(a) + tuple(b)
    center = Fraction(n * m, 2)
    observed_distance = abs(u_statistic(a, b) - center)
    extreme = 0
    total = 0
    indices = range(n + m)
    for a_indices in itertools.combinations(indices, n):
        a_set = frozenset(a_indices)
        new_a = [pooled[i] for i in a_indices]
        new_b = [pooled[j] for j in indices if j not in a_set]
        distance = abs(u_statistic(new_a, new_b) - center)
        if distance >= observed_distance:
            extreme += 1
        total += 1
    return Fraction(extreme, total)


def exact_test(a: Sequence[str], b: Sequence[str]) -> ExactTestResult:
    """Exact Mann-Whitney test for two groups of decimal strings."""
    parsed_a = parse_values(a, "a")
    parsed_b = parse_values(b, "b")
    return ExactTestResult(
        u=u_statistic(parsed_a, parsed_b),
        p=exact_p_value(parsed_a, parsed_b),
        n=len(parsed_a),
        m=len(parsed_b),
    )
