"""Addons for MSImagingExperiment: the plural/iterator forms of the spectrum view family.

``spectrum_view(i)`` itself is a C++ binding; ``spectrum_views()`` and
``iter_spectrum_views()`` are generated here so the imaging class shares the
view contract of MSExperiment (see OWNERSHIP.md).
"""
from __future__ import annotations

from . import register_element_views

register_element_views("MSImagingExperiment", "spectrum", "getNrSpectra", "spectra")
