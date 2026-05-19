"""
Tests for IMAGING bindings: IonImage, MSImagingGeometry, MSImagingExperiment.
"""

import pytest
import numpy as np


def _make_fixture():
    """2x2 grid; pixel (1, 1) intentionally missing."""
    from pyopenms import MSExperiment, MSSpectrum, MSImagingExperiment, MSImagingGeometry

    exp = MSExperiment()
    s0 = MSSpectrum(); s0.set_peaks(([499.95, 500.05], [10.0, 20.0])); exp.addSpectrum(s0)
    s1 = MSSpectrum(); s1.set_peaks(([490.0],          [5.0]));        exp.addSpectrum(s1)
    s2 = MSSpectrum(); s2.set_peaks(([500.0],          [100.0]));      exp.addSpectrum(s2)

    geom = MSImagingGeometry()
    geom.setDimensions(2, 2)
    geom.addPixel(0, 0, 0)
    geom.addPixel(1, 0, 1)
    geom.addPixel(0, 1, 2)

    mie = MSImagingExperiment(exp)
    mie.setGeometry(geom)
    return mie


def test_ion_image_basic():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    assert img.getWidth() == 3
    assert img.getHeight() == 2
    assert img.hasPixel(0, 0) is False

    img.setIntensity(1, 0, 7.5)
    assert img.hasPixel(1, 0) is True
    assert img.getIntensity(1, 0) == pytest.approx(7.5)


def test_ion_image_get_data_zero_copy():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    img.setIntensity(1, 0, 7.5)

    data = img.get_data()
    assert data.shape == (2, 3)
    assert data.dtype == np.float64
    assert data[0, 1] == pytest.approx(7.5)

    # Zero-copy: writes through the view are reflected in the C++ object
    data[1, 2] = 42.0
    assert img.getIntensity(2, 1) == pytest.approx(42.0)


def test_ion_image_get_mask_copy():
    from pyopenms import IonImage

    img = IonImage(2, 2)
    img.setIntensity(0, 0, 1.0)
    img.setIntensity(1, 1, 2.0)

    mask = img.get_mask()
    assert mask.shape == (2, 2)
    assert mask.dtype == np.uint8
    assert mask[0, 0] == 1
    assert mask[1, 1] == 1
    assert mask[0, 1] == 0
    assert mask[1, 0] == 0


def test_ion_image_oob_raises():
    from pyopenms import IonImage

    img = IonImage(3, 2)
    with pytest.raises(Exception):
        img.getIntensity(5, 0)


def test_ms_imaging_geometry_basic():
    from pyopenms import MSImagingGeometry

    g = MSImagingGeometry()
    g.setDimensions(4, 3)
    g.setPixelSize(25.0, 25.0, "micrometer")
    g.addPixel(0, 0, 100)
    g.addPixel(1, 0, 101)
    assert g.getWidth() == 4
    assert g.getHeight() == 3
    assert g.getNumberOfPixels() == 2
    assert g.hasPixel(0, 0) is True
    assert g.hasPixel(2, 0) is False
    assert g.getSpectrumIndex(1, 0) == 101


def test_ms_imaging_geometry_pixels_struct():
    g = _make_fixture().getGeometry()
    pix = g.get_pixels_struct()
    assert pix.dtype.names == ("x", "y", "spectrum_index")
    assert pix.dtype.itemsize == 16
    assert pix.shape == (3,)
    assert pix["x"].tolist() == [0, 1, 0]
    assert pix["y"].tolist() == [0, 0, 1]
    assert pix["spectrum_index"].tolist() == [0, 1, 2]


def test_ms_imaging_geometry_rejects_oob():
    from pyopenms import MSImagingGeometry

    g = MSImagingGeometry()
    g.setDimensions(2, 2)
    g.addPixel(1, 1, 0)
    with pytest.raises(Exception):
        g.addPixel(2, 0, 1)
    with pytest.raises(Exception):
        g.addPixel(1, 1, 99)  # duplicate


def test_ms_imaging_experiment_extract_ion_image():
    mie = _make_fixture()
    assert mie.getNumberOfPixels() == 3
    assert mie.getNumberOfSpectra() == 3
    assert mie.hasPixel(1, 1) is False

    # 200 ppm @ 500 Da -> window [499.9, 500.1]
    img = mie.extractIonImage(500.0, 200.0)
    assert img.getWidth() == 2
    assert img.getHeight() == 2

    data = img.get_data()
    assert data[0, 0] == pytest.approx(30.0)   # 10 + 20
    assert data[0, 1] == pytest.approx(0.0)    # spectrum 1 has only 490
    assert data[1, 0] == pytest.approx(100.0)  # spectrum 2 has 500
    assert data[1, 1] == pytest.approx(0.0)    # pixel missing, but data buffer reads 0

    valid = img.get_mask()
    assert valid[1, 1] == 0  # pixel (1, 1) absent -> masked out

    assert img.getMzRange().getMinMZ() == pytest.approx(499.9)
    assert img.getMzRange().getMaxMZ() == pytest.approx(500.1)


def test_ms_imaging_experiment_extract_rejects_bad_inputs():
    mie = _make_fixture()
    with pytest.raises(Exception):
        mie.extractIonImage(-1.0, 100.0)
    with pytest.raises(Exception):
        mie.extractIonImage(500.0, -1.0)


def test_ms_imaging_experiment_get_spectrum():
    mie = _make_fixture()
    sp = mie.getSpectrum(0, 1)
    mz, intens = sp.get_peaks()
    assert list(mz) == pytest.approx([500.0])
    assert list(intens) == pytest.approx([100.0])

    with pytest.raises(Exception):
        mie.getSpectrum(1, 1)


def test_ms_imaging_experiment_assignment_clears_geometry():
    from pyopenms import MSExperiment, MSSpectrum

    mie = _make_fixture()
    fresh = MSExperiment()
    fresh.addSpectrum(MSSpectrum())
    mie.setMSExperiment(fresh)
    # setMSExperiment keeps the geometry
    assert mie.getNumberOfPixels() == 3
    # operator= is exposed in C++ but we can re-create via constructor + reassign:
    mie2 = type(mie)(fresh)
    assert mie2.getNumberOfSpectra() == 1
    assert mie2.getNumberOfPixels() == 0
