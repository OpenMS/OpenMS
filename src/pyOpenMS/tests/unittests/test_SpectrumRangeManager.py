"""Tests for SpectrumRangeManager.byMSLevel() (issue #8751)."""

import pytest
import numpy as np
import pyopenms


def _make_experiment():
    """Create an MSExperiment with MS1 and MS2 spectra at known ranges."""
    exp = pyopenms.MSExperiment()

    # MS1 spectrum at RT=1.0, MZ=[100, 200], intensity=[1000, 2000]
    s1 = pyopenms.MSSpectrum()
    s1.setRT(1.0)
    s1.setMSLevel(1)
    s1.set_peaks((np.array([100.0, 200.0]), np.array([1000.0, 2000.0])))
    exp.addSpectrum(s1)

    # MS1 spectrum at RT=3.0, MZ=[150, 250], intensity=[500, 1500]
    s2 = pyopenms.MSSpectrum()
    s2.setRT(3.0)
    s2.setMSLevel(1)
    s2.set_peaks((np.array([150.0, 250.0]), np.array([500.0, 1500.0])))
    exp.addSpectrum(s2)

    # MS2 spectrum at RT=2.0, MZ=[50, 500], intensity=[100, 300]
    s3 = pyopenms.MSSpectrum()
    s3.setRT(2.0)
    s3.setMSLevel(2)
    s3.set_peaks((np.array([50.0, 500.0]), np.array([100.0, 300.0])))
    exp.addSpectrum(s3)

    exp.updateRanges()
    return exp


def test_byMSLevel_returns_SpectrumRanges():
    exp = _make_experiment()
    sr = exp.spectrumRanges()
    ms1 = sr.byMSLevel(1)
    assert type(ms1).__name__ == "SpectrumRanges"


def test_byMSLevel_ms1_ranges():
    exp = _make_experiment()
    sr = exp.spectrumRanges()
    ms1 = sr.byMSLevel(1)

    assert ms1.getMinRT() == pytest.approx(1.0)
    assert ms1.getMaxRT() == pytest.approx(3.0)
    assert ms1.getMinMZ() == pytest.approx(100.0)
    assert ms1.getMaxMZ() == pytest.approx(250.0)
    assert ms1.getMinIntensity() == pytest.approx(500.0)
    assert ms1.getMaxIntensity() == pytest.approx(2000.0)


def test_byMSLevel_ms2_ranges():
    exp = _make_experiment()
    sr = exp.spectrumRanges()
    ms2 = sr.byMSLevel(2)

    assert ms2.getMinRT() == pytest.approx(2.0)
    assert ms2.getMaxRT() == pytest.approx(2.0)
    assert ms2.getMinMZ() == pytest.approx(50.0)
    assert ms2.getMaxMZ() == pytest.approx(500.0)
    assert ms2.getMinIntensity() == pytest.approx(100.0)
    assert ms2.getMaxIntensity() == pytest.approx(300.0)


def test_byMSLevel_differs_from_global():
    """Per-level ranges should differ from global ranges when levels have different extents."""
    exp = _make_experiment()
    sr = exp.spectrumRanges()

    # Global MZ includes MS2's wider range [50, 500]
    assert sr.getMinMZ() == pytest.approx(50.0)
    assert sr.getMaxMZ() == pytest.approx(500.0)

    # MS1-only MZ is narrower [100, 250]
    ms1 = sr.byMSLevel(1)
    assert ms1.getMinMZ() == pytest.approx(100.0)
    assert ms1.getMaxMZ() == pytest.approx(250.0)


def test_byMSLevel_invalid_level_raises():
    exp = _make_experiment()
    sr = exp.spectrumRanges()

    with pytest.raises(RuntimeError):
        sr.byMSLevel(3)


def test_SpectrumRanges_importable():
    assert hasattr(pyopenms, "SpectrumRanges")
