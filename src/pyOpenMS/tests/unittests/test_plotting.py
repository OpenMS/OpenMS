"""
``plot_spectrum()`` and ``mirror_plot_spectrum()`` label each peak with the entry
of the spectrum's first StringDataArray. The bindings return those entries as
``str``; pyOpenMS 3.5.0 returned ``bytes``. Both must be drawn (OpenMS issue #10260).
"""

import pytest

import pyopenms

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pyopenms.plotting import mirror_plot_spectrum, plot_spectrum  # noqa: E402

ANNOTATIONS = ["b1", "y2", "b3"]


def _annotated_spectrum():
    spec = pyopenms.MSSpectrum()
    spec.set_peaks(([100.0, 200.0, 300.0], [10.0, 30.0, 20.0]))
    annotations = pyopenms.StringDataArray()
    for annotation in ANNOTATIONS:
        annotations.push_back(annotation)
    spec.setStringDataArrays([annotations])
    return spec


def test_plot_spectrum_draws_annotations():
    fig, ax = plt.subplots()
    try:
        plot_spectrum(_annotated_spectrum(), ax=ax)
        assert sorted(text.get_text() for text in ax.texts) == sorted(ANNOTATIONS)
    finally:
        plt.close(fig)


def test_mirror_plot_spectrum_draws_annotations():
    fig, ax = plt.subplots()
    try:
        mirror_plot_spectrum(_annotated_spectrum(), _annotated_spectrum(), ax=ax)
        assert sorted(text.get_text() for text in ax.texts) == sorted(ANNOTATIONS * 2)
    finally:
        plt.close(fig)
