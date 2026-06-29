"""Tests for the BrukerTimsFile loader bindings.

The class is only exposed when OpenMS was built with ``WITH_OPENTIMS=ON``,
so every test feature-detects via ``hasattr(pyopenms, "BrukerTimsFile")``.
Real-data tests additionally require ``OPENTIMS_DDA_TEST_DATA`` /
``OPENTIMS_DIA_TEST_DATA`` env vars (set by CMake when those configure-time
paths were provided), mirroring the C++ ``BrukerTimsFile_test`` pattern.
"""
import os
import pytest
import pyopenms

pytestmark = pytest.mark.skipif(
    not hasattr(pyopenms, "BrukerTimsFile"),
    reason="pyOpenMS built without WITH_OPENTIMS",
)


def test_config_defaults():
    cfg = pyopenms.BrukerTimsFile.Config()
    assert cfg.calibrate is False
    assert cfg.load_ms1 is True
    assert cfg.ms1_centroid_mz_ppm == pytest.approx(10.0)
    assert cfg.ms1_centroid_im_pct == pytest.approx(0.0)
    assert cfg.dia_ms2_n_neighbors == 0
    assert cfg.dia_ms2_min_support == 1
    assert cfg.dia_ms2_centroid is False
    assert cfg.ms1_n_neighbors == 0
    assert cfg.ms1_min_support == 0
    assert cfg.export_mode == pyopenms.BrukerTimsFile.Config.ExportMode.AUTO
    assert cfg.tims_calibration_strategy == pyopenms.BrukerTimsFile.Config.TimsCalibrationStrategy.AUTO
    assert cfg.pressure_compensation == pyopenms.BrukerTimsFile.Config.PressureCompensation.NONE
    assert cfg.bruker_sdk_path == ""


def test_config_field_assignment():
    cfg = pyopenms.BrukerTimsFile.Config()
    cfg.load_ms1 = False
    cfg.ms1_centroid_mz_ppm = 5.0
    cfg.export_mode = pyopenms.BrukerTimsFile.Config.ExportMode.FRAME
    cfg.bruker_sdk_path = "/opt/bruker/libtimsdata.so"
    assert cfg.load_ms1 is False
    assert cfg.ms1_centroid_mz_ppm == pytest.approx(5.0)
    assert cfg.export_mode == pyopenms.BrukerTimsFile.Config.ExportMode.FRAME
    assert cfg.bruker_sdk_path == "/opt/bruker/libtimsdata.so"


def test_export_mode_enum_values():
    Em = pyopenms.BrukerTimsFile.Config.ExportMode
    assert {Em.AUTO, Em.SPECTRUM, Em.FRAME} == {Em.AUTO, Em.SPECTRUM, Em.FRAME}


def test_calibration_strategy_enum_values():
    Cs = pyopenms.BrukerTimsFile.Config.TimsCalibrationStrategy
    assert {Cs.AUTO, Cs.BRUKER_SDK, Cs.RATIONAL, Cs.LINEAR} == {
        Cs.AUTO, Cs.BRUKER_SDK, Cs.RATIONAL, Cs.LINEAR
    }


def test_pressure_compensation_enum_values():
    Pc = pyopenms.BrukerTimsFile.Config.PressureCompensation
    assert {Pc.NONE, Pc.GLOBAL, Pc.PER_FRAME} == {Pc.NONE, Pc.GLOBAL, Pc.PER_FRAME}


def test_centroid_algo_enum_values():
    Ca = pyopenms.BrukerTimsFile.Config.CentroidAlgo
    assert hasattr(Ca, "OFF")
    assert hasattr(Ca, "GREEDY2D")
    assert hasattr(Ca, "HILL_BASED")


def test_config_centroiding_defaults():
    cfg = pyopenms.BrukerTimsFile.Config()
    Ca = pyopenms.BrukerTimsFile.Config.CentroidAlgo
    assert cfg.ms1_centroid_algo == Ca.OFF
    assert cfg.ms2_centroid_algo == Ca.OFF
    assert cfg.ms2_centroid_mz_ppm == pytest.approx(20.0)
    assert cfg.centroid_valley_factor == pytest.approx(1.3)
    assert cfg.ms1_centroid_min_hill_length == 1
    assert cfg.ms2_centroid_min_hill_length == 2
    assert cfg.centroid_max_scan_gap == 0
    assert cfg.expose_hill_bounds is False
    assert cfg.isotopic_prefilter is False
    assert cfg.isotopic_prefilter_tol_ppm == pytest.approx(50.0)


def test_config_range_filter_defaults():
    cfg = pyopenms.BrukerTimsFile.Config()
    assert cfg.frame_id_min == 0
    assert cfg.frame_id_max > 0  # UINT32_MAX
    assert cfg.rt_min_sec == pytest.approx(0.0)
    import math
    assert math.isinf(cfg.rt_max_sec)


def test_config_centroid_algo_assignment():
    cfg = pyopenms.BrukerTimsFile.Config()
    Ca = pyopenms.BrukerTimsFile.Config.CentroidAlgo
    cfg.ms1_centroid_algo = Ca.HILL_BASED
    cfg.ms2_centroid_algo = Ca.GREEDY2D
    assert cfg.ms1_centroid_algo == Ca.HILL_BASED
    assert cfg.ms2_centroid_algo == Ca.GREEDY2D


def test_dia_streaming_metadata_defaults():
    meta = pyopenms.BrukerTimsFile.DIAStreamingMetadata()
    assert meta.nr_ms1_spectra == 0
    assert list(meta.nr_ms2_spectra) == []
    assert list(meta.boundaries) == []


def test_load_nonexistent_path_raises():
    f = pyopenms.BrukerTimsFile()
    with pytest.raises(Exception):
        f.load("/nonexistent/path.d")


@pytest.mark.skipif(
    not os.environ.get("OPENTIMS_DDA_TEST_DATA"),
    reason="OPENTIMS_DDA_TEST_DATA not set",
)
def test_load_dda_dataset():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    f = pyopenms.BrukerTimsFile()
    exp = f.load(path)
    assert exp.getNrSpectra() > 0


@pytest.mark.skipif(
    not os.environ.get("OPENTIMS_DDA_TEST_DATA"),
    reason="OPENTIMS_DDA_TEST_DATA not set",
)
def test_load_dda_ms2_only():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    f = pyopenms.BrukerTimsFile()
    cfg = pyopenms.BrukerTimsFile.Config()
    cfg.load_ms1 = False
    exp = f.load(path, cfg)
    # MS2-only: no MS1 spectra
    ms1_count = sum(1 for s in exp if s.getMSLevel() == 1)
    assert ms1_count == 0


@pytest.mark.skipif(
    not os.environ.get("OPENTIMS_DIA_TEST_DATA"),
    reason="OPENTIMS_DIA_TEST_DATA not set",
)
def test_read_dia_metadata():
    path = os.environ["OPENTIMS_DIA_TEST_DATA"]
    f = pyopenms.BrukerTimsFile()
    meta, settings = f.readDIAMetadata(path)
    assert meta.nr_ms1_spectra > 0
    assert len(meta.boundaries) > 0
    assert len(meta.nr_ms2_spectra) == len(meta.boundaries)
