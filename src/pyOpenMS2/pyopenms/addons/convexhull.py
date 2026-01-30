"""Addon methods for ConvexHull2D to return numpy arrays."""

from __future__ import annotations
import numpy as np


def _setup_convexhull(cls):
    """Wrap ConvexHull2D methods to return/accept numpy arrays."""
    _orig_getHullPoints = cls.getHullPoints
    _orig_setHullPoints = cls.setHullPoints
    _orig_addPoints = cls.addPoints

    def getHullPoints(self):
        pts = _orig_getHullPoints(self)
        if len(pts) == 0:
            return np.empty((0, 2), dtype=np.float64)
        return np.array(pts, dtype=np.float64)

    def setHullPoints(self, points):
        if hasattr(points, 'tolist'):
            points = points.tolist()
        _orig_setHullPoints(self, [[float(x), float(y)] for x, y in points])

    def addPoints(self, points):
        if hasattr(points, 'tolist'):
            points = points.tolist()
        _orig_addPoints(self, [[float(x), float(y)] for x, y in points])

    cls.getHullPoints = getHullPoints
    cls.setHullPoints = setHullPoints
    cls.addPoints = addPoints


# Registration: this will be called from apply_addons
WRAP_CLASSES = {"ConvexHull2D": _setup_convexhull}
