"""
Addon methods for Peak1D class.
"""

from __future__ import annotations

from . import addon


@addon("Peak1D")
def to_tuple(self) -> tuple:
    """Return peak as (mz, intensity) tuple."""
    return (self.getMZ(), self.getIntensity())
