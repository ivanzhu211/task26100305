"""Exact Mann-Whitney U test with ties and shift confidence-set inversion.

All inputs are non-empty iterables of finite decimal strings, interpreted as
exact numeric values (via fractions.Fraction).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Iterator, NamedTuple, Optional, Sequence, Tuple, Union

Number = Union[Fraction, int]
MaybeInfinite = Optional[Fraction]  # None means unbounded (-inf / +inf)


class ExactTestResult(NamedTuple):
    """Result of the exact Mann-Whitney test."""

    u: Fraction
    p_value: Fraction


class Interval(NamedTuple):
    """One connected piece of a confidence set.

    ``lower``/``upper`` are exact Fraction endpoints; ``None`` means the
    interval is unbounded on that side.  The *_closed flags tell whether a
    finite endpoint belongs to the set.
    """

    lower: MaybeInfinite
    upper: MaybeInfinite
    lower_closed: bool
    upper_closed: bool


class ConfidenceSet:
    """Full confidence set for the location shift (A - B)."""

    def __init__(self, intervals: Sequence[Interval], confidence: Fraction):
        self.intervals: Tuple[Interval, ...] = tuple(intervals)
        self.confidence = confidence

    def __contains__(self, theta: object) -> bool:
        t = Fraction(theta)  # type: ignore[arg-type]
        for iv in self.intervals:
            if iv.lower is not None:
                if t < iv.lower or (t == iv.lower and not iv.lower_closed):
                    continue
            if iv.upper is not None:
                if t > iv.upper or (t == iv.upper and not iv.upper_closed):
                    continue
            return True
        return False

    def __iter__(self) -> Iterator[Interval]:
        return iter(self.intervals)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return f"ConfidenceSet(intervals={self.intervals!r}, confidence={self.confidence!s})"


def exact_test(a: Sequence[str], b: Sequence[str]) -> ExactTestResult:
    """Observed U statistic and exact two-sided p-value."""
    av = _parse_group(a, "a")
    bv = _parse_group(b, "b")
    u_obs = _u_statistic(av, bv)
    p = _exact_p_value(av, bv, u_obs)
    return ExactTestResult(u=u_obs, p_value=p)


def _parse_group(values: Sequence[str], name: str) -> Tuple[Fraction, ...]:
    if values is None or len(values) == 0:
        raise ValueError(f"group {name!r} must be non-empty")
    parsed = []
    for v in values:
        try:
            parsed.append(Fraction(v))
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"group {name!r}: invalid decimal string {v!r}") from exc
    return tuple(parsed)


def _u_statistic(av: Sequence[Fraction], bv: Sequence[Fraction]) -> Fraction:
    u = Fraction(0)
    half = Fraction(1, 2)
    for x in av:
        for y in bv:
            if x > y:
                u += 1
            elif x == y:
                u += half
    return u


def _exact_p_value(av: Sequence[Fraction], bv: Sequence[Fraction], u_obs: Fraction) -> Fraction:
    """Exact two-sided p-value by full label re-enumeration.

    Every choice of n positions (out of n + m combined observations) to form
    group A is equally likely; equal values are still distinct observations.
    """
    n, m = len(av), len(bv)
    combined = tuple(av) + tuple(bv)
    total = n + m
    center = Fraction(n * m, 2)
    d_obs = abs(u_obs - center)

    # comp[i][j]: score of observation i versus observation j.
    half = Fraction(1, 2)
    comp = []
    for x in combined:
        row = []
        for y in combined:
            if x > y:
                row.append(1)
            elif x == y:
                row.append(half)
            else:
                row.append(0)
        comp.append(row)

    assignments = 0
    extreme = 0
    for a_idx in combinations(range(total), n):
        in_a = [False] * total
        for i in a_idx:
            in_a[i] = True
        u = Fraction(0)
        for i in a_idx:
            row = comp[i]
            for j in range(total):
                if not in_a[j]:
                    u += row[j]
        if abs(u - center) >= d_obs:
            extreme += 1
        assignments += 1
    return Fraction(extreme, assignments)


def invert_shift(a: Sequence[str], b: Sequence[str], confidence: object) -> ConfidenceSet:
    """Full confidence set for the shift theta (A - B) at ``confidence``.

    ``confidence`` may be given as 0.90/0.95/0.99 or 90/95/99.
    """
    av = _parse_group(a, "a")
    bv = _parse_group(b, "b")
    conf, alpha = _parse_confidence(confidence)

    # p(theta) only changes where some a_i - theta crosses some b_j,
    # i.e. at theta equal to a cross difference a_i - b_j.
    diffs = sorted({x - y for x in av for y in bv})

    # Membership atoms, in increasing theta order:
    #   gap_0, point_1, gap_1, ..., point_k, gap_k
    # gap_i is the open interval between two consecutive critical points
    # (unbounded on the outer sides); point_i is the critical value itself.
    def p_at(theta: Fraction) -> Fraction:
        shifted = [x - theta for x in av]
        u = _u_statistic(shifted, bv)
        return _exact_p_value(shifted, bv, u)

    def accepted(theta: Fraction) -> bool:
        return p_at(theta) > alpha

    atoms: list[tuple[str, MaybeInfinite, MaybeInfinite, bool]] = []
    # (kind, lo_bound, hi_bound, included); points have lo == hi.
    prev: Optional[Fraction] = None
    for d in diffs:
        if prev is None:
            atoms.append(("gap", None, d, accepted(d - 1)))
        else:
            mid = (prev + d) / 2
            atoms.append(("gap", prev, d, accepted(mid)))
        atoms.append(("point", d, d, accepted(d)))
        prev = d
    atoms.append(("gap", prev, None, accepted(prev + 1)))

    # Merge runs of included atoms into maximal intervals.  The upper end of
    # a run is the hi bound of its last included atom: closed for an included
    # critical point, open for a gap (the boundary point is excluded).
    intervals: list[Interval] = []
    cur_lo: MaybeInfinite = None
    cur_lo_closed = False
    last_hi: MaybeInfinite = None
    last_hi_closed = False
    active = False
    for kind, lo, hi, included in atoms:
        if included:
            if not active:
                cur_lo, cur_lo_closed = lo, kind == "point"
                active = True
            last_hi, last_hi_closed = hi, kind == "point"
        else:
            if active:
                intervals.append(Interval(cur_lo, last_hi, cur_lo_closed, last_hi_closed))
                active = False
    if active:
        intervals.append(Interval(cur_lo, last_hi, cur_lo_closed, last_hi_closed))

    return ConfidenceSet(intervals, conf)


_CONFIDENCE_TABLE = {
    Fraction(9, 10): Fraction(1, 10),
    Fraction(19, 20): Fraction(1, 20),
    Fraction(99, 100): Fraction(1, 100),
}


def _parse_confidence(confidence: object) -> Tuple[Fraction, Fraction]:
    if isinstance(confidence, str):
        text = confidence.strip()
        if text.endswith("%"):
            text = text[:-1].strip()
    elif isinstance(confidence, float):
        text = str(confidence)
    else:
        text = confidence
    try:
        conf = Fraction(text)  # type: ignore[arg-type]
    except (ValueError, TypeError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid confidence level {confidence!r}") from exc
    if conf in (90, 95, 99):
        conf = Fraction(conf, 100)
    if conf not in _CONFIDENCE_TABLE:
        raise ValueError(
            f"unsupported confidence level {confidence!r}; use 90%, 95% or 99%"
        )
    return conf, _CONFIDENCE_TABLE[conf]
