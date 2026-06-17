"""
Tests for creating and reading ion mobility data via pyOpenMS.

Refactored from tools/scripts/create_im.py and
tools/scripts/create_im_with_swathIm.py — these were standalone scripts
that generated simulated ion mobility mzML files. The logic is preserved
here as round-trip tests: build an MSExperiment in memory, write it to
mzML, read it back, and verify the contents.
"""

import os

import numpy as np
import pytest
import pyopenms as poms


# ---------------------------------------------------------------------------
# Helpers (from create_im_with_swathIm.py)
# ---------------------------------------------------------------------------

def _add_ms1_precursor(sp, fda, mz, im_idx, base_intensity,
                       im_offset_idx=2, nr_peaks=10):
    """Add an MS1 precursor across ion-mobility bins to *sp* and *fda*."""
    for i in range(nr_peaks):
        p = poms.Peak1D()
        p.setMZ(mz)
        p.setIntensity(base_intensity + (i - 5))
        sp.push_back(p)
        fda.push_back((im_idx - im_offset_idx + i * im_offset_idx) / 500)


def _create_ms2_precursor(allmz, allint, allim, imCenterIdx, imWidth,
                          base_int, baseFragmentMz=100, nr_peaks=5,
                          im_bins=(300, 700), imWinLim=None):
    """Accumulate fragment peaks for one precursor across IM bins."""
    for im_idx in range(*im_bins):
        for i in range(nr_peaks):
            in_im = (imCenterIdx - imWidth) < im_idx < (imCenterIdx + imWidth)
            if imWinLim is not None:
                in_im = in_im and (imWinLim[0] < im_idx < imWinLim[1])
            if in_im:
                apex_dist = abs(imCenterIdx - im_idx)
                intensity = base_int * (i + 1) - base_int * (i + 1) * apex_dist / imWidth
                allmz.append(baseFragmentMz + i)
                allint.append(intensity)
                allim.append(im_idx / 500.0)


# ---------------------------------------------------------------------------
# Experiment builders
# ---------------------------------------------------------------------------

def _build_im_experiment():
    """Build a simple ion-mobility MSExperiment (from create_im.py)."""
    exp = poms.MSExperiment()

    NR_RT_SAMPLES = 50
    NR_IM_BINS = 300
    NR_PEAKS = 5

    # --- MS1 spectra ---
    for rt_idx in range(NR_RT_SAMPLES):
        sp = poms.MSSpectrum()
        sp.setRT(rt_idx)
        sp.setMSLevel(1)

        base_int = max(5, 100 - abs(20 - rt_idx) * 100 / 25.0)

        im_values = []

        for i in range(100):
            p = poms.Peak1D()
            p.setMZ(100 + i)
            p.setIntensity(base_int + i)
            sp.push_back(p)
            im_values.append(10.0)

        for i in range(10):
            p = poms.Peak1D()
            p.setMZ(412.502)
            p.setIntensity(base_int + (i - 5))
            sp.push_back(p)
            im_values.append(99.0 + (i - 5))

        for i in range(10):
            p = poms.Peak1D()
            p.setMZ(417.502)
            p.setIntensity(base_int + (i - 5))
            sp.push_back(p)
            im_values.append(152.0 + (i - 5))

        fda = poms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.array(im_values, dtype=np.float32))

        sp.setFloatDataArrays([fda])
        exp.addSpectrum(sp)

    # --- MS2 spectra (single SWATH window) ---
    for rt_idx in range(NR_RT_SAMPLES):
        base_int = 100 - abs(25 - rt_idx) * 100 / 25.0

        allmz, allint, allim = [], [], []

        for im_idx in range(NR_IM_BINS):
            for i in range(NR_PEAKS):
                if 90 < im_idx < 110:
                    apex_dist = abs(100 - im_idx)
                    allmz.append(100 + i)
                    allint.append(base_int * (i + 1) - base_int * (i + 1) * apex_dist / 10.0)
                    allim.append(im_idx / 500.0)

            for i in range(NR_PEAKS):
                if 130 < im_idx < 170:
                    apex_dist = abs(150 - im_idx)
                    allmz.append(100 + i)
                    allint.append(base_int * (i + 1) - base_int * (i + 1) * apex_dist / 20.0)
                    allim.append(im_idx / 500.0)

        fda = poms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.array(allim, dtype=np.float32))

        sframe = poms.MSSpectrum()
        sframe.setMSLevel(2)
        sframe.setRT(rt_idx + 0.2)
        sframe.setFloatDataArrays([fda])

        prec = poms.Precursor()
        prec.setMZ(412.5)
        prec.setIsolationWindowUpperOffset(12.5)
        prec.setIsolationWindowLowerOffset(12.5)
        sframe.setPrecursors([prec])
        sframe.set_peaks((allmz, allint))
        sframe.sortByPosition()
        exp.addSpectrum(sframe)

    exp.sortSpectra()
    return exp


def _build_im_swath_experiment():
    """Build an IM + dual-SWATH MSExperiment (from create_im_with_swathIm.py)."""
    exp = poms.MSExperiment()

    NR_RT_SAMPLES = 50
    NR_PEAKS = 5

    # --- MS1 spectra ---
    for rt_idx in range(NR_RT_SAMPLES):
        base_int = max(5, 100 - abs(20 - rt_idx) * 100 / 25.0)

        sp = poms.MSSpectrum()
        sp.setRT(rt_idx)
        sp.setMSLevel(1)

        im_values = []
        for i in range(100):
            p = poms.Peak1D()
            p.setMZ(100 + i)
            p.setIntensity(base_int + i)
            sp.push_back(p)
            im_values.append(300 / 500)

        fda = poms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.array(im_values, dtype=np.float32))

        _add_ms1_precursor(sp, fda, 412.502, 350, base_int)
        _add_ms1_precursor(sp, fda, 417.502, 450, base_int, im_offset_idx=5)
        _add_ms1_precursor(sp, fda, 422.502, 550, base_int)

        sp.setFloatDataArrays([fda])
        exp.addSpectrum(sp)

    # --- MS2 spectra, SWATH window 1 (IM 300–450) ---
    for rt_idx in range(NR_RT_SAMPLES):
        base_int = 100 - abs(25 - rt_idx) * 100 / 25.0
        allmz, allint, allim = [], [], []

        _create_ms2_precursor(allmz, allint, allim,
                              imCenterIdx=350, imWidth=10, base_int=base_int,
                              baseFragmentMz=100)
        _create_ms2_precursor(allmz, allint, allim,
                              imCenterIdx=450, imWidth=20, base_int=base_int,
                              imWinLim=(300, 450), baseFragmentMz=100.01)

        fda = poms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.array(allim, dtype=np.float32))

        sframe = poms.MSSpectrum()
        sframe.setMSLevel(2)
        sframe.setRT(rt_idx + 0.2)
        sframe.setMetaValue('ion mobility lower limit', 300 / 500)
        sframe.setMetaValue('ion mobility upper limit', 450 / 500)

        prec = poms.Precursor()
        prec.setMZ(412.5)
        prec.setIsolationWindowUpperOffset(12.5)
        prec.setIsolationWindowLowerOffset(12.5)
        sframe.setPrecursors([prec])
        sframe.set_peaks((allmz, allint))
        sframe.setFloatDataArrays([fda])
        sframe.sortByPosition()
        exp.addSpectrum(sframe)

    # --- MS2 spectra, SWATH window 2 (IM 350–600) ---
    for rt_idx in range(NR_RT_SAMPLES):
        base_int = 100 - abs(25 - rt_idx) * 100 / 25.0
        allmz, allint, allim = [], [], []

        _create_ms2_precursor(allmz, allint, allim,
                              imCenterIdx=350, imWidth=10, base_int=base_int,
                              imWinLim=(350, 600))
        _create_ms2_precursor(allmz, allint, allim,
                              imCenterIdx=450, imWidth=20, base_int=base_int,
                              imWinLim=(350, 600), baseFragmentMz=100.01)
        _create_ms2_precursor(allmz, allint, allim,
                              imCenterIdx=550, imWidth=10, base_int=base_int,
                              imWinLim=(350, 600), baseFragmentMz=100.02)

        fda = poms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.set_data(np.array(allim, dtype=np.float32))

        sframe = poms.MSSpectrum()
        sframe.setMSLevel(2)
        sframe.setRT(rt_idx + 0.3)
        sframe.setMetaValue('ion mobility lower limit', 350 / 500)
        sframe.setMetaValue('ion mobility upper limit', 600 / 500)

        prec = poms.Precursor()
        prec.setMZ(412.5)
        prec.setIsolationWindowUpperOffset(12.5)
        prec.setIsolationWindowLowerOffset(12.5)
        sframe.setPrecursors([prec])
        sframe.set_peaks((allmz, allint))
        sframe.setFloatDataArrays([fda])
        sframe.sortByPosition()
        exp.addSpectrum(sframe)

    exp.sortSpectra()
    return exp


# ---------------------------------------------------------------------------
# Tests — simple ion mobility experiment (create_im.py)
# ---------------------------------------------------------------------------

class TestIonMobilityExperiment:
    """Round-trip tests for a simple ion-mobility MSExperiment."""

    @pytest.fixture(scope="class")
    def im_exp(self):
        return _build_im_experiment()

    def test_spectrum_count(self, im_exp):
        """50 MS1 + 50 MS2 = 100 spectra."""
        assert im_exp.getNrSpectra() == 100

    def test_ms_levels(self, im_exp):
        ms1 = [s for s in im_exp if s.getMSLevel() == 1]
        ms2 = [s for s in im_exp if s.getMSLevel() == 2]
        assert len(ms1) == 50
        assert len(ms2) == 50

    def test_ms1_has_ion_mobility_fda(self, im_exp):
        ms1 = [s for s in im_exp if s.getMSLevel() == 1]
        for sp in ms1:
            fdas = sp.getFloatDataArrays()
            assert len(fdas) == 1
            assert fdas[0].getName() == "Ion Mobility"

    def test_ms2_has_ion_mobility_fda(self, im_exp):
        ms2 = [s for s in im_exp if s.getMSLevel() == 2]
        for sp in ms2:
            fdas = sp.getFloatDataArrays()
            assert len(fdas) == 1
            assert fdas[0].getName() == "Ion Mobility"

    def test_ms2_precursor(self, im_exp):
        ms2 = [s for s in im_exp if s.getMSLevel() == 2]
        for sp in ms2:
            precs = sp.getPrecursors()
            assert len(precs) == 1
            assert abs(precs[0].getMZ() - 412.5) < 0.1

    def test_ms1_peak_count(self, im_exp):
        """Each MS1 has 100 background + 10 + 10 = 120 peaks."""
        ms1 = [s for s in im_exp if s.getMSLevel() == 1]
        for sp in ms1:
            assert sp.size() == 120

    def test_mzml_roundtrip(self, im_exp, tmp_path):
        path = str(tmp_path / "im_test.mzML")

        f = poms.MzMLFile()
        opts = f.getOptions()
        opts.setCompression(True)
        f.setOptions(opts)
        f.store(path, im_exp)

        assert os.path.isfile(path)
        assert os.path.getsize(path) > 0

        loaded = poms.MSExperiment()
        f.load(path, loaded)
        assert loaded.getNrSpectra() == im_exp.getNrSpectra()


# ---------------------------------------------------------------------------
# Tests — ion mobility + SWATH experiment (create_im_with_swathIm.py)
# ---------------------------------------------------------------------------

class TestIonMobilitySwathExperiment:
    """Round-trip tests for an IM + dual-SWATH MSExperiment."""

    @pytest.fixture(scope="class")
    def swath_exp(self):
        return _build_im_swath_experiment()

    def test_spectrum_count(self, swath_exp):
        """50 MS1 + 50 SWATH-1 MS2 + 50 SWATH-2 MS2 = 150 spectra."""
        assert swath_exp.getNrSpectra() == 150

    def test_ms_levels(self, swath_exp):
        ms1 = [s for s in swath_exp if s.getMSLevel() == 1]
        ms2 = [s for s in swath_exp if s.getMSLevel() == 2]
        assert len(ms1) == 50
        assert len(ms2) == 100

    def test_ms2_im_metadata(self, swath_exp):
        """All MS2 spectra should have IM window limit meta values."""
        ms2 = [s for s in swath_exp if s.getMSLevel() == 2]
        for sp in ms2:
            assert sp.metaValueExists('ion mobility lower limit')
            assert sp.metaValueExists('ion mobility upper limit')
            lower = float(sp.getMetaValue('ion mobility lower limit'))
            upper = float(sp.getMetaValue('ion mobility upper limit'))
            assert upper > lower

    def test_swath_window_1_im_limits(self, swath_exp):
        """First 50 MS2 spectra (SWATH-1) have IM window 0.6–0.9."""
        ms2 = [s for s in swath_exp if s.getMSLevel() == 2]
        # SWATH-1 has RT offset 0.2, SWATH-2 has 0.3
        swath1 = [s for s in ms2
                  if abs(float(s.getMetaValue('ion mobility lower limit')) - 300 / 500) < 1e-6]
        assert len(swath1) == 50
        for sp in swath1:
            assert abs(float(sp.getMetaValue('ion mobility upper limit')) - 450 / 500) < 1e-6

    def test_swath_window_2_im_limits(self, swath_exp):
        """Second 50 MS2 spectra (SWATH-2) have IM window 0.7–1.2."""
        ms2 = [s for s in swath_exp if s.getMSLevel() == 2]
        swath2 = [s for s in ms2
                  if abs(float(s.getMetaValue('ion mobility lower limit')) - 350 / 500) < 1e-6]
        assert len(swath2) == 50
        for sp in swath2:
            assert abs(float(sp.getMetaValue('ion mobility upper limit')) - 600 / 500) < 1e-6

    def test_ms1_precursors_present(self, swath_exp):
        """MS1 spectra should have peaks around 412.5, 417.5, 422.5 m/z."""
        ms1 = [s for s in swath_exp if s.getMSLevel() == 1]
        sp = ms1[25]  # mid-experiment, good intensity
        mzs, _ = sp.get_peaks()
        mz_set = set(np.round(mzs, 1))
        assert 412.5 in mz_set
        assert 417.5 in mz_set
        assert 422.5 in mz_set

    def test_ms2_has_fragment_peaks(self, swath_exp):
        """MS2 spectra near the RT apex should have fragment peaks."""
        ms2 = [s for s in swath_exp if s.getMSLevel() == 2]
        # Pick one near RT=25 (apex for base_int)
        apex_spectra = [s for s in ms2 if 25.0 < s.getRT() < 26.0]
        assert len(apex_spectra) > 0
        for sp in apex_spectra:
            mzs, intens = sp.get_peaks()
            assert len(mzs) > 0
            assert np.any(np.array(intens) > 0)

    def test_mzml_roundtrip(self, swath_exp, tmp_path):
        path = str(tmp_path / "im_swath_test.mzML")

        f = poms.MzMLFile()
        opts = f.getOptions()
        opts.setCompression(True)
        f.setOptions(opts)
        f.store(path, swath_exp)

        assert os.path.isfile(path)
        assert os.path.getsize(path) > 0

        loaded = poms.MSExperiment()
        f.load(path, loaded)
        assert loaded.getNrSpectra() == swath_exp.getNrSpectra()

        # Verify MS levels survive round-trip
        ms1_loaded = [s for s in loaded if s.getMSLevel() == 1]
        ms2_loaded = [s for s in loaded if s.getMSLevel() == 2]
        assert len(ms1_loaded) == 50
        assert len(ms2_loaded) == 100
