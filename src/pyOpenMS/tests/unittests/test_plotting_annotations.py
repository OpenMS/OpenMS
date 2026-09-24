"""
Regression tests for pyopenms.plotting on annotated spectra (OpenMS issue #10260).

``plot_spectrum()`` reads peak annotations from the first ``StringDataArray``.
The bindings return these as ``str`` (they used to be ``bytes``), and the
unconditional ``.decode()`` crashed with
``AttributeError: 'str' object has no attribute 'decode'``, which also broke
``mirror_plot_spectrum()``.
"""

import pytest

import pyopenms

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pyopenms.plotting import plot_spectrum, mirror_plot_spectrum  # noqa: E402

ANNOTATIONS = ["b2+", "y3+", "y4+"]


def _make_annotated_spectrum():
    """An MS2 spectrum whose first StringDataArray annotates every peak."""
    spec = pyopenms.MSSpectrum()
    spec.setMSLevel(2)
    spec.set_peaks(([200.1, 350.2, 480.3], [100.0, 250.0, 50.0]))
    annot = pyopenms.StringDataArray()
    annot.setName("IonNames")
    for a in ANNOTATIONS:
        annot.push_back(a)
    spec.setStringDataArrays([annot])
    return spec


def _texts(ax):
    return [t.get_text() for t in ax.texts]


def test_plot_spectrum_with_annotations():
    spec = _make_annotated_spectrum()
    fig, ax = plt.subplots()
    try:
        plot_spectrum(spec, ax=ax)
        assert _texts(ax) == ANNOTATIONS
    finally:
        plt.close(fig)


def test_mirror_plot_spectrum_with_annotations():
    spec_top = _make_annotated_spectrum()
    spec_bottom = _make_annotated_spectrum()
    fig, ax = plt.subplots()
    try:
        mirror_plot_spectrum(spec_top, spec_bottom, ax=ax)
        assert _texts(ax) == ANNOTATIONS + ANNOTATIONS
    finally:
        plt.close(fig)


if __name__ == "__main__":
    test_plot_spectrum_with_annotations()
    test_mirror_plot_spectrum_with_annotations()
    print("All plotting annotation tests passed!")
