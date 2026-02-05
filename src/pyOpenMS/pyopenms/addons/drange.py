"""Addon methods for DRange classes."""

from __future__ import annotations
from . import addon


@addon("DRange1")
def __repr__(self):
    return f"DRange1({float(self.minX())}, {float(self.maxX())})"


@addon("DRange1")
def __str__(self):
    return f"DRange1({float(self.minX())}, {float(self.maxX())})"


@addon("DRange2")
def __repr__(self):
    return f"DRange2(({float(self.minX())}, {float(self.minY())}), ({float(self.maxX())}, {float(self.maxY())}))"


@addon("DRange2")
def __str__(self):
    return f"DRange2(({float(self.minX())}, {float(self.minY())}), ({float(self.maxX())}, {float(self.maxY())}))"
