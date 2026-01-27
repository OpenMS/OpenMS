"""
Addon methods for Peak1D class.
"""

from __future__ import annotations

from . import addon


@addon("Peak1D")
def __repr__(self) -> str:
    """Return string representation of Peak1D."""
    return f"Peak1D(mz={self.getMZ():.4f}, intensity={self.getIntensity():.2f})"


@addon("Peak1D")
def __str__(self) -> str:
    """Return user-friendly string of Peak1D."""
    return f"({self.getMZ():.4f}, {self.getIntensity():.2f})"


@addon("Peak1D")
def to_tuple(self) -> tuple:
    """Return peak as (mz, intensity) tuple."""
    return (self.getMZ(), self.getIntensity())
