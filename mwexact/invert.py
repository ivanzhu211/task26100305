"""Inversion of the exact test into a shift (A - B) confidence set."""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from mwexact.exact import exact_p_value, parse_values
from mwexact.intervals import ConfidenceSet, Interval

_SUPPORTED_CONFIDENCES = frozenset(
    {Fraction(9, 10), Fraction(19, 20), Fraction(99, 100)}
)


def _parse_confidence(confidence: float | str | Fraction) -> Fraction:
    if isinstance(confidence, Fraction):
        key = confidence
    elif isinstance(confidence, float):
        key = Fraction(str(confidence))
    elif isinstance(confidence, int):
        key = Fraction(confidence, 100)
    elif isinstance(confidence, str):
        text = confidence.strip()
        if text.endswith("%"):
            key = Fraction(text[:-1].strip()) / 100
        else:
            key = Fraction(text)
    else:
        raise TypeError(f"unsupported confidence type: {type(confidence)!r}")
    if key not in _SUPPORTED_CONFIDENCES:
        allowed = ", ".join(str(c) for c in sorted(_SUPPORTED_CONFIDENCES))
        raise ValueError(f"confidence must be one of {allowed}; got {key}")
    return key


def _shifted(values: Sequence[Fraction], theta: Fraction) -> list[Fraction]:
    return [value - theta for value in values]


def _segments(
    a: Sequence[Fraction], b: Sequence[Fraction], alpha: Fraction
) -> list[Interval]:
    """Accepted regions between consecutive distinct cross differences."""
    differences = sorted({x - y for x in a for y in b})
    first, last = differences[0], differences[-1]

    candidates: list[Interval] = [Interval(None, first, False, False)]
    for index, point in enumerate(differences):
        candidates.append(Interval(point, point, True, True))
        if index + 1 < len(differences):
            candidates.append(Interval(point, differences[index + 1], False, False))
    candidates.append(Interval(last, None, False, False))

    def representative(interval: Interval) -> Fraction:
        if interval.low is None:
            return interval.high - 1
        if interval.high is None:
            return interval.low + 1
        return (interval.low + interval.high) / 2

    accepted: list[Interval] = []
    for interval in candidates:
        theta = representative(interval)
        p_value = exact_p_value(_shifted(a, theta), b)
        if p_value > alpha:
            accepted.append(interval)
    return accepted


def _merge(segments: Sequence[Interval]) -> tuple[Interval, ...]:
    """Fuse adjacent accepted segments into maximal intervals."""
    merged: list[Interval] = []
    for segment in segments:
        if merged and merged[-1].high == segment.low:
            previous = merged[-1]
            merged[-1] = Interval(
                previous.low, segment.high, previous.low_closed, segment.high_closed
            )
        else:
            merged.append(segment)
    return tuple(merged)


def invert_shift(
    a: Sequence[str], b: Sequence[str], confidence: float | str | Fraction
) -> ConfidenceSet:
    """Exact confidence set {theta : p(theta) > 1 - confidence} for A - B."""
    level = _parse_confidence(confidence)
    parsed_a = parse_values(a, "a")
    parsed_b = parse_values(b, "b")
    alpha = 1 - level
    return ConfidenceSet(_merge(_segments(parsed_a, parsed_b, alpha)), level)
