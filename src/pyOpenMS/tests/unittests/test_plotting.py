"""Tests for pyopenms.plotting."""

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt  # noqa: E402

import pyopenms  # noqa: E402
from pyopenms.plotting import plot_spectrum  # noqa: E402


def _annotated_spectrum():
    """An MS2 spectrum whose first StringDataArray holds one ion name per peak."""
    spec = pyopenms.MSSpectrum()
    spec.setMSLevel(2)
    spec.setNativeID('scan=1')
    spec.set_peaks([np.array([100.0, 200.0, 300.0]),
                    np.array([1.0, 2.0, 3.0], dtype=np.float32)])
    sda = pyopenms.StringDataArray()
    sda.setName('IonNames')
    for ion in ['b2+', 'y3+', 'b4+']:
        sda.push_back(ion)
    spec.setStringDataArrays([sda])
    return spec


class TestPlotSpectrum:
    def teardown_method(self):
        plt.close('all')

    def test_annotated_spectrum_plots(self):
        """Annotations are read from the StringDataArray, which yields str. This used to
        raise AttributeError because the annotations were .decode()'d as if they were
        bytes, a leftover from the Cython bindings."""
        ax = plot_spectrum(_annotated_spectrum(), annotate_ions=True)

        drawn = [t.get_text() for t in ax.texts]
        assert drawn == ['b2+', 'y3+', 'b4+']

    def test_unannotated_spectrum_plots(self):
        """A spectrum without a StringDataArray draws empty annotations, not ion names."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(2)
        spec.set_peaks([np.array([100.0, 200.0]),
                        np.array([1.0, 2.0], dtype=np.float32)])

        ax = plot_spectrum(spec, annotate_ions=True)
        assert [t.get_text() for t in ax.texts] == ['', '']

    def test_string_data_array_length_mismatch_is_ignored(self):
        """A StringDataArray that does not have one entry per peak is not used."""
        spec = _annotated_spectrum()
        sda = pyopenms.StringDataArray()
        sda.setName('IonNames')
        sda.push_back('b2+')  # one entry, three peaks
        spec.setStringDataArrays([sda])

        ax = plot_spectrum(spec, annotate_ions=True)
        assert [t.get_text() for t in ax.texts] == ['', '', '']
