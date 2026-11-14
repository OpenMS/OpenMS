"""Addon methods for MSImagingExperiment class."""

from __future__ import annotations
from . import register_element_views


register_element_views("MSImagingExperiment", "spectrum", "getNumberOfPixels", "pixel spectra")
