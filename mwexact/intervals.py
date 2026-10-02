"""Exact interval and confidence-set representations on the real line."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Interval:
    """Single exact interval; None bounds denote unbounded sides."""

    low: Fraction | None
    high: Fraction | None
    low_closed: bool
    high_closed: bool

    def contains(self, theta: Fraction) -> bool:
        if self.low is not None:
            if theta < self.low or (theta == self.low and not self.low_closed):
                return False
        if self.high is not None:
            if theta > self.high or (theta == self.high and not self.high_closed):
                return False
        return True

    def __str__(self) -> str:
        left = "[" if self.low_closed else "("
        right = "]" if self.high_closed else ")"
        low = "-inf" if self.low is None else str(self.low)
        high = "+inf" if self.high is None else str(self.high)
        return f"{left}{low}, {high}{right}"


@dataclass(frozen=True)
class ConfidenceSet:
    """Finite union of disjoint exact intervals."""

    intervals: tuple[Interval, ...]
    confidence: Fraction

    def contains(self, theta: Fraction | str | int | float) -> bool:
        point = theta if isinstance(theta, Fraction) else Fraction(theta)
        return any(interval.contains(point) for interval in self.intervals)

    def __str__(self) -> str:
        body = " U ".join(str(interval) for interval in self.intervals)
        return f"ConfidenceSet(confidence={self.confidence}, {body or 'empty'})"
