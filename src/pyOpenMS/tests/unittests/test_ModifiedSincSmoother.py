import pytest
import pyopenms


def test_modified_sinc_smoother_init_default():
    smoother = pyopenms.ModifiedSincSmoother()
    params = smoother.getParameters().asDict()
    assert "degree" in params
    assert "m" in params
    assert "is_ms1" in params
    assert params["degree"] == 6
    assert params["m"] == 7
    assert params["is_ms1"] is False


def test_modified_sinc_smoother_init_parameterized():
    smoother = pyopenms.ModifiedSincSmoother(True, 4, 15)
    params = smoother.getParameters().asDict()
    assert params["degree"] == 4
    assert params["m"] == 15
    assert params["is_ms1"] is True


def test_modified_sinc_smoother_parameters():
    smoother = pyopenms.ModifiedSincSmoother()
    params = smoother.getParameters()
    params.setValue("degree", 4)
    params.setValue("m", 10)
    params.setValue("is_ms1", "true")
    smoother.setParameters(params)

    updated_params = smoother.getParameters().asDict()
    assert updated_params["degree"] == 4
    assert updated_params["m"] == 10
    assert updated_params["is_ms1"] is True


def test_modified_sinc_smoother_smooth_vector():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 5)
    # Triangular peak with flat baseline
    raw = [0.0, 0.0, 2.0, 10.0, 50.0, 10.0, 2.0, 0.0, 0.0]
    smoothed = smoother.smooth(raw)

    assert len(smoothed) == len(raw)
    # Peak apex should still be at the center (index 4)
    assert max(smoothed) == smoothed[4]
    # Smoothed peak apex intensity should be attenuated compared to sharp peak
    assert smoothed[4] < raw[4]
    assert smoothed[4] > raw[3]


def test_modified_sinc_smoother_filter_spectrum():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 5)
    spec = pyopenms.MSSpectrum()
    mzs = [100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0]
    ints = [1.0, 3.0, 15.0, 50.0, 15.0, 3.0, 1.0]
    spec.set_peaks((mzs, ints))

    smoother.filter(spec)

    assert spec.size() == len(mzs)
    filtered_mzs, filtered_ints = spec.get_peaks()
    # M/Z coordinates must not be altered by smoothing
    for orig_mz, filtered_mz in zip(mzs, filtered_mzs):
        assert abs(orig_mz - filtered_mz) < 1e-9
    # Apex intensity smoothed
    assert filtered_ints[3] < ints[3]


def test_modified_sinc_smoother_filter_chromatogram():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 5)
    chrom = pyopenms.MSChromatogram()
    rts = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
    ints = [0.5, 2.0, 10.0, 40.0, 10.0, 2.0, 0.5]
    chrom.set_peaks((rts, ints))

    smoother.filter(chrom)

    assert chrom.size() == len(rts)
    filtered_rts, filtered_ints = chrom.get_peaks()
    for orig_rt, filtered_rt in zip(rts, filtered_rts):
        assert abs(orig_rt - filtered_rt) < 1e-9
    assert filtered_ints[3] < ints[3]


def test_modified_sinc_smoother_filter_mobilogram():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 5)
    mob = pyopenms.Mobilogram()
    dts = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    ints = [1.0, 4.0, 12.0, 35.0, 12.0, 4.0, 1.0]
    mob.set_peaks((dts, ints))

    smoother.filter(mob)

    assert mob.size() == len(dts)
    filtered_dts, filtered_ints = mob.get_peaks()
    for orig_dt, filtered_dt in zip(dts, filtered_dts):
        assert abs(orig_dt - filtered_dt) < 1e-9
    assert filtered_ints[3] < ints[3]


def test_modified_sinc_smoother_filter_experiment():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 5)
    exp = pyopenms.MSExperiment()

    spec = pyopenms.MSSpectrum()
    spec.set_peaks(([100.0, 101.0, 102.0, 103.0, 104.0], [1.0, 5.0, 20.0, 5.0, 1.0]))
    exp.addSpectrum(spec)

    chrom = pyopenms.MSChromatogram()
    chrom.set_peaks(([10.0, 11.0, 12.0, 13.0, 14.0], [2.0, 8.0, 30.0, 8.0, 2.0]))
    exp.addChromatogram(chrom)

    smoother.filterExperiment(exp)

    assert exp.getNrSpectra() == 1
    assert exp.getNrChromatograms() == 1
    filtered_spec = exp.getSpectrum(0)
    filtered_chrom = exp.getChromatogram(0)
    assert filtered_spec.size() == 5
    assert filtered_chrom.size() == 5


def test_modified_sinc_smoother_static_methods():
    m = pyopenms.ModifiedSincSmoother.bandwidthToM(False, 6, 0.1)
    assert isinstance(m, int)
    assert m > 0

    noise_gain_m = pyopenms.ModifiedSincSmoother.noiseGainToM(False, 6, 0.2)
    assert isinstance(noise_gain_m, int)
    assert noise_gain_m > 0

    bw = pyopenms.ModifiedSincSmoother.savitzkyGolayBandwidth(4, 5)
    assert isinstance(bw, float)
    assert bw > 0.0


def test_modified_sinc_smoother_empty_containers():
    smoother = pyopenms.ModifiedSincSmoother()

    assert smoother.smooth([]) == []

    empty_spec = pyopenms.MSSpectrum()
    smoother.filter(empty_spec)
    assert empty_spec.size() == 0

    empty_chrom = pyopenms.MSChromatogram()
    smoother.filter(empty_chrom)
    assert empty_chrom.size() == 0

    empty_mob = pyopenms.Mobilogram()
    smoother.filter(empty_mob)
    assert empty_mob.size() == 0

    empty_exp = pyopenms.MSExperiment()
    smoother.filterExperiment(empty_exp)
    assert empty_exp.size() == 0


def test_modified_sinc_smoother_constant_signal():
    smoother = pyopenms.ModifiedSincSmoother(False, 4, 7)
    constant_level = 42.0
    raw = [constant_level] * 25
    smoothed = smoother.smooth(raw)

    assert len(smoothed) == len(raw)
    for val in smoothed:
        assert abs(val - constant_level) < 1e-4
