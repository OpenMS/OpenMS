"""Tests for the BrukerTimsFile loader bindings.

The class is only exposed when OpenMS was built with ``WITH_OPENTIMS=ON``,
so every test feature-detects via ``hasattr(pyopenms, "BrukerTimsFile")``.
Real-data tests additionally require ``OPENTIMS_DDA_TEST_DATA`` /
``OPENTIMS_DIA_TEST_DATA`` env vars (set by CMake when those configure-time
paths were provided), mirroring the C++ ``BrukerTimsFile_test`` pattern.
"""
import math
import os
import pytest
import pyopenms

pytestmark = pytest.mark.skipif(
    not hasattr(pyopenms, "BrukerTimsFile"),
    reason="pyOpenMS built without WITH_OPENTIMS",
)

# skip markers of the real-data tests: the same conditions and reasons as the inline markers further down
needs_dda_data = pytest.mark.skipif(
    not os.environ.get("OPENTIMS_DDA_TEST_DATA"),
    reason="OPENTIMS_DDA_TEST_DATA not set",
)
needs_dia_data = pytest.mark.skipif(
    not os.environ.get("OPENTIMS_DIA_TEST_DATA"),
    reason="OPENTIMS_DIA_TEST_DATA not set",
)


def test_config_defaults():
    cfg = pyopenms.BrukerTimsFile.Config()
    assert cfg.calibrate is False
    assert cfg.load_ms1 is True
    assert cfg.ms1_centroid_mz_ppm == pytest.approx(10.0)
    assert cfg.ms1_centroid_im_pct == pytest.approx(0.0)
    assert cfg.ms2_centroid_algo == pyopenms.BrukerTimsFile.Config.CentroidAlgo.OFF
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


def test_hillbased_dia_no_rt_aggregation_config():
    cfg = pyopenms.BrukerTimsFile.Config()
    cfg.ms2_centroid_algo = pyopenms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED
    cfg.dia_ms2_n_neighbors = 0
    assert cfg.ms2_centroid_algo == pyopenms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED
    assert cfg.dia_ms2_n_neighbors == 0


def _enum_members(enum_type):
    """(name, value) of every member, in declaration order."""
    return [(member.name, member.value) for member in enum_type]


def _assert_members_distinct(enum_type, expected_count):
    members = list(enum_type)
    assert len(members) == expected_count
    assert len(set(members)) == expected_count
    for i, a in enumerate(members):
        for j, b in enumerate(members):
            assert (a == b) == (i == j)


def test_export_mode_enum_values():
    Em = pyopenms.BrukerTimsFile.Config.ExportMode
    # names and values of the C++ enumerators: enum ExportMode { AUTO, SPECTRUM, FRAME }
    assert _enum_members(Em) == [("AUTO", 0), ("SPECTRUM", 1), ("FRAME", 2)]
    _assert_members_distinct(Em, 3)
    assert Em.FRAME is Em(2)


def test_calibration_strategy_enum_values():
    Cs = pyopenms.BrukerTimsFile.Config.TimsCalibrationStrategy
    # enum class TimsCalibrationStrategy { AUTO, BRUKER_SDK, RATIONAL, LINEAR }
    assert _enum_members(Cs) == [("AUTO", 0), ("BRUKER_SDK", 1), ("RATIONAL", 2), ("LINEAR", 3)]
    _assert_members_distinct(Cs, 4)
    assert Cs.LINEAR is Cs(3)


def test_pressure_compensation_enum_values():
    Pc = pyopenms.BrukerTimsFile.Config.PressureCompensation
    # enum class PressureCompensation { NONE, GLOBAL, PER_FRAME }
    assert _enum_members(Pc) == [("NONE", 0), ("GLOBAL", 1), ("PER_FRAME", 2)]
    _assert_members_distinct(Pc, 3)
    assert Pc.PER_FRAME is Pc(2)


def test_centroid_algo_enum_values():
    Ca = pyopenms.BrukerTimsFile.Config.CentroidAlgo
    # enum class CentroidAlgo { OFF, GREEDY2D, HILL_BASED }
    assert _enum_members(Ca) == [("OFF", 0), ("GREEDY2D", 1), ("HILL_BASED", 2)]
    _assert_members_distinct(Ca, 3)
    assert Ca.HILL_BASED is Ca(2)


def _header_defaults():
    """Every field of BrukerTimsFile::Config with its default, in the order of BrukerTimsFile.h."""
    Config = pyopenms.BrukerTimsFile.Config
    return {
        "calibration_tolerance": 0.0,
        "calibrate": False,
        "load_ms1": True,
        "load_ms2": True,
        "frame_id_min": 0,
        "frame_id_max": 2**32 - 1,  # std::numeric_limits<uint32_t>::max()
        "rt_min_sec": 0.0,
        "rt_max_sec": math.inf,  # std::numeric_limits<double>::infinity()
        "ms1_centroid_algo": Config.CentroidAlgo.OFF,
        "ms2_centroid_algo": Config.CentroidAlgo.OFF,
        "ms1_centroid_mz_ppm": 10.0,
        "ms1_centroid_im_pct": 0.0,
        "ms1_centroid_max_peaks": 100000,
        "ms2_centroid_mz_ppm": 20.0,
        "centroid_valley_factor": 1.3,
        "ms1_centroid_min_hill_length": 1,
        "ms2_centroid_min_hill_length": 2,
        "centroid_max_scan_gap": 0,
        "expose_hill_bounds": False,
        "isotopic_prefilter": False,
        "isotopic_prefilter_tol_ppm": 50.0,
        "dia_ms2_n_neighbors": 0,
        "dia_ms2_min_support": 1,
        "dia_ms2_centroid": False,
        "ms1_n_neighbors": 0,
        "ms1_min_support": 0,
        "ms1_max_rt_distance_sec": 0.0,
        "export_mode": Config.ExportMode.AUTO,
        "tims_calibration_strategy": Config.TimsCalibrationStrategy.AUTO,
        "pressure_compensation": Config.PressureCompensation.NONE,
        "bruker_sdk_path": "",
    }


def test_config_binds_every_field():
    Config = pyopenms.BrukerTimsFile.Config
    bound = {
        name for name in dir(Config)
        if not name.startswith("_") and isinstance(getattr(Config, name), property)
    }
    assert bound == set(_header_defaults())


def test_config_defaults_match_header():
    cfg = pyopenms.BrukerTimsFile.Config()
    for name, expected in _header_defaults().items():
        value = getattr(cfg, name)
        # the type as well: False == 0 and 1 == 1.0 in Python
        assert type(value) is type(expected), name
        assert value == expected, name


def test_config_range_and_centroiding_field_assignment():
    Config = pyopenms.BrukerTimsFile.Config
    # the fields that were not bound before; every value differs from the default and is
    # exactly representable in the C++ type (float for ms2_centroid_mz_ppm)
    assigned = {
        "frame_id_min": 17,
        "frame_id_max": 4242,
        "rt_min_sec": 12.5,
        "rt_max_sec": 99.25,
        "ms1_centroid_algo": Config.CentroidAlgo.HILL_BASED,
        "ms2_centroid_mz_ppm": 7.5,
        "centroid_valley_factor": 2.75,
        "ms1_centroid_min_hill_length": 3,
        "ms2_centroid_min_hill_length": 5,
        "centroid_max_scan_gap": 4,
        "expose_hill_bounds": True,
        "isotopic_prefilter": True,
        "isotopic_prefilter_tol_ppm": 8.5,
    }
    defaults = _header_defaults()
    assert all(defaults[name] != value for name, value in assigned.items())

    # assign everything first and read back afterwards: two names bound to the same
    # C++ member would show up as a value overwritten by a later assignment
    cfg = Config()
    for name, value in assigned.items():
        setattr(cfg, name, value)
    for name, value in assigned.items():
        assert type(getattr(cfg, name)) is type(value), name
        assert getattr(cfg, name) == value, name
    # ... and no other field changed
    for name, value in defaults.items():
        if name not in assigned:
            assert getattr(cfg, name) == value, name

    # back to the defaults, including the two extreme ones
    for name in assigned:
        setattr(cfg, name, defaults[name])
    assert cfg.frame_id_max == 2**32 - 1
    assert cfg.rt_max_sec == math.inf
    for name, value in defaults.items():
        assert getattr(cfg, name) == value, name


def test_config_frame_id_is_uint32():
    cfg = pyopenms.BrukerTimsFile.Config()
    # frame IDs are uint32_t in C++: a value outside that range is refused, not wrapped
    for name in ("frame_id_min", "frame_id_max"):
        default = getattr(cfg, name)
        for out_of_range in (-1, 2**32):
            with pytest.raises(TypeError):
                setattr(cfg, name, out_of_range)
        assert getattr(cfg, name) == default
        setattr(cfg, name, 2**32 - 1)
        assert getattr(cfg, name) == 2**32 - 1


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


@needs_dda_data
def test_load_ms1_only():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    f = pyopenms.BrukerTimsFile()
    cfg = pyopenms.BrukerTimsFile.Config()
    cfg.frame_id_min, cfg.frame_id_max = 1, 40
    both = f.load(path, cfg)
    cfg.load_ms2 = False
    ms1_only = f.load(path, cfg)
    # MS1-only: no MS2 spectra, and the MS1 spectra are those of the load with both levels
    assert {s.getMSLevel() for s in both} == {1, 2}
    assert {s.getMSLevel() for s in ms1_only} == {1}
    expected = [s for s in both if s.getMSLevel() == 1]
    assert ms1_only.getNrSpectra() == len(expected)
    assert all(a == b for a, b in zip(ms1_only, expected))


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


def _config(**fields):
    cfg = pyopenms.BrukerTimsFile.Config()
    for name, value in fields.items():
        setattr(cfg, name, value)
    return cfg


def _frame_id(native_id):
    """Frame ID of a spectrum native ID ('frame=<id>' or 'frame=<id> <more>')."""
    first = native_id.split()[0]
    assert first.startswith("frame=")
    return int(first[len("frame="):])


def _load_summary(path, **fields):
    """(frame ID, MS level, RT, number of peaks) of every spectrum of a load, in order."""
    exp = pyopenms.BrukerTimsFile().load(path, _config(**fields))
    return [(_frame_id(s.getNativeID()), s.getMSLevel(), s.getRT(), s.size()) for s in exp]


def _ms_level_peaks(exp, ms_level):
    return sum(s.size() for s in exp if s.getMSLevel() == ms_level)


@needs_dda_data
def test_frame_id_range_restricts_load():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    FRAME = pyopenms.BrukerTimsFile.Config.ExportMode.FRAME

    # FRAME mode returns one spectrum per frame. The full load is restricted to the MS2
    # frames, which hold about 2% of the peaks of the run.
    full = _load_summary(path, export_mode=FRAME, load_ms1=False)
    assert len(full) > 1000
    assert {ms_level for _, ms_level, _, _ in full} == {2}

    lo, hi = full[100][0], full[400][0]
    sub = _load_summary(path, export_mode=FRAME, load_ms1=False, frame_id_min=lo, frame_id_max=hi)
    assert 0 < len(sub) < len(full)
    assert all(lo <= frame_id <= hi for frame_id, _, _, _ in sub)
    # the range selects spectra and leaves them as they are; both ends are inclusive
    assert sub == [s for s in full if lo <= s[0] <= hi]
    assert sub[0][0] == lo and sub[-1][0] == hi

    # with the MS1 frames: every frame of the range, in order (frame IDs are 1-based and dense)
    both = _load_summary(path, export_mode=FRAME, frame_id_min=1, frame_id_max=40)
    assert [frame_id for frame_id, _, _, _ in both] == list(range(1, 41))
    assert {ms_level for _, ms_level, _, _ in both} == {1, 2}
    assert [s for s in both if s[1] == 2] == [s for s in full if s[0] <= 40]

    # a single frame, and a range behind the last frame
    assert [s[0] for s in _load_summary(path, export_mode=FRAME, frame_id_min=7, frame_id_max=7)] == [7]
    behind = full[-1][0] + 1000
    assert _load_summary(path, export_mode=FRAME, frame_id_min=behind, frame_id_max=behind + 10) == []


@needs_dda_data
def test_rt_range_restricts_load():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    FRAME = pyopenms.BrukerTimsFile.Config.ExportMode.FRAME

    full = _load_summary(path, export_mode=FRAME, load_ms1=False)
    rts = [rt for _, _, rt, _ in full]
    assert rts == sorted(rts)

    # the retention times of two frames as limits: both ends are inclusive
    rt_lo, rt_hi = rts[100], rts[400]
    sub = _load_summary(path, export_mode=FRAME, load_ms1=False, rt_min_sec=rt_lo, rt_max_sec=rt_hi)
    assert 0 < len(sub) < len(full)
    assert all(rt_lo <= rt <= rt_hi for _, _, rt, _ in sub)
    assert sub == [s for s in full if rt_lo <= s[2] <= rt_hi]
    assert sub[0][2] == rt_lo and sub[-1][2] == rt_hi

    # limits between two frames
    between_lo = (rts[100] + rts[101]) / 2
    between_hi = (rts[399] + rts[400]) / 2
    sub = _load_summary(path, export_mode=FRAME, load_ms1=False, rt_min_sec=between_lo, rt_max_sec=between_hi)
    assert sub == full[101:400]

    # with the MS1 frames (FRAME mode: the retention time of a spectrum is that of its frame)
    both = _load_summary(path, export_mode=FRAME, rt_min_sec=rt_lo, rt_max_sec=rt_lo + 5.0)
    assert {ms_level for _, ms_level, _, _ in both} == {1, 2}
    assert all(rt_lo <= rt <= rt_lo + 5.0 for _, _, rt, _ in both)
    assert [s for s in both if s[1] == 2] == [s for s in full if rt_lo <= s[2] <= rt_lo + 5.0]
    # frame IDs are dense, and the range holds all frames between its first and its last one
    frame_ids = [frame_id for frame_id, _, _, _ in both]
    assert frame_ids == list(range(frame_ids[0], frame_ids[-1] + 1))

    # the retention time range is intersected with the frame ID range
    mid = full[250][0]
    sub = _load_summary(path, export_mode=FRAME, load_ms1=False,
                        rt_min_sec=rt_lo, rt_max_sec=rt_hi, frame_id_min=mid)
    assert sub == [s for s in full if rt_lo <= s[2] <= rt_hi and s[0] >= mid]
    assert 0 < len(sub) < len([s for s in full if rt_lo <= s[2] <= rt_hi])

    # nothing behind the end of the run
    assert _load_summary(path, export_mode=FRAME, rt_min_sec=rts[-1] + 1000.0) == []


@needs_dia_data
def test_frame_id_range_restricts_dia_load_and_metadata():
    path = os.environ["OPENTIMS_DIA_TEST_DATA"]
    f = pyopenms.BrukerTimsFile()
    cfg = _config(frame_id_min=1, frame_id_max=40)

    exp = f.load(path, cfg)
    assert exp.getNrSpectra() > 0
    assert all(1 <= _frame_id(s.getNativeID()) <= 40 for s in exp)
    n_ms1 = sum(1 for s in exp if s.getMSLevel() == 1)
    n_ms2 = sum(1 for s in exp if s.getMSLevel() == 2)
    assert n_ms1 > 0 and n_ms2 > 0

    # the metadata of the whole run needs no peak data, so it serves as the unrestricted reference
    meta_full, _ = f.readDIAMetadata(path)
    meta_sub, _ = f.readDIAMetadata(path, cfg)
    assert len(meta_sub.boundaries) == len(meta_full.boundaries)
    assert 0 < meta_sub.nr_ms1_spectra < meta_full.nr_ms1_spectra
    assert all(0 < sub < full for sub, full in zip(meta_sub.nr_ms2_spectra, meta_full.nr_ms2_spectra))

    # one MS1 spectrum per MS1 frame; the counts per window are the frames of the window, an
    # upper limit for the spectra (a frame without peaks in a window gives no spectrum)
    assert n_ms1 == meta_sub.nr_ms1_spectra
    assert n_ms2 <= sum(meta_sub.nr_ms2_spectra)
    # fewer spectra than the whole run has MS1 frames alone
    assert exp.getNrSpectra() < meta_full.nr_ms1_spectra


@needs_dda_data
def test_ms1_hill_centroiding_fields_reach_the_reader():
    path = os.environ["OPENTIMS_DDA_TEST_DATA"]
    Config = pyopenms.BrukerTimsFile.Config
    HILL = Config.CentroidAlgo.HILL_BASED
    f = pyopenms.BrukerTimsFile()

    def load(**fields):
        # frame 1 is an MS1 frame
        return f.load(path, _config(export_mode=Config.ExportMode.FRAME,
                                    frame_id_min=1, frame_id_max=1, **fields))

    def ms1_peaks(**fields):
        return _ms_level_peaks(load(**fields), 1)

    raw = ms1_peaks()
    hill = ms1_peaks(ms1_centroid_algo=HILL)
    assert 0 < hill < raw

    # a hill has to span more IM scans
    assert 0 < ms1_peaks(ms1_centroid_algo=HILL, ms1_centroid_min_hill_length=2) < hill
    # hills are joined across one empty IM scan
    assert 0 < ms1_peaks(ms1_centroid_algo=HILL, centroid_max_scan_gap=1) < hill
    # a smaller valley factor splits more hills
    assert (ms1_peaks(ms1_centroid_algo=HILL, centroid_valley_factor=1.0) > hill
            > ms1_peaks(ms1_centroid_algo=HILL, centroid_valley_factor=100.0) > 0)

    # the isotopic prefilter drops peaks without a partner, more of them with a narrow tolerance
    prefiltered = ms1_peaks(ms1_centroid_algo=HILL, isotopic_prefilter=True)
    assert 0 < prefiltered < hill
    assert 0 < ms1_peaks(ms1_centroid_algo=HILL, isotopic_prefilter=True,
                         isotopic_prefilter_tol_ppm=5.0) < prefiltered

    def ms1_array_names(exp):
        return [[a.getName() for a in s.getFloatDataArrays()] for s in exp if s.getMSLevel() == 1]

    im_array = "raw inverse reduced ion mobility array"
    assert ms1_array_names(load(ms1_centroid_algo=HILL)) == [[im_array]]
    assert ms1_array_names(load(ms1_centroid_algo=HILL, expose_hill_bounds=True)) == [
        [im_array, "im lower bound", "im upper bound", "m/z lower bound", "m/z upper bound"]
    ]


@needs_dia_data
def test_ms2_hill_centroiding_fields_reach_the_reader():
    path = os.environ["OPENTIMS_DIA_TEST_DATA"]
    HILL = pyopenms.BrukerTimsFile.Config.CentroidAlgo.HILL_BASED
    f = pyopenms.BrukerTimsFile()

    def ms2_peaks(**fields):
        # frames 2 to 9 are the eight MS2 frames of the first cycle
        exp = f.load(path, _config(load_ms1=False, frame_id_min=1, frame_id_max=10, **fields))
        return _ms_level_peaks(exp, 2)

    raw = ms2_peaks()
    hill = ms2_peaks(ms2_centroid_algo=HILL)
    assert 0 < hill < raw

    # the default minimum hill length for MS2 is 2
    assert hill < ms2_peaks(ms2_centroid_algo=HILL, ms2_centroid_min_hill_length=1) < raw
    # a narrower m/z tolerance links fewer peaks into hills of at least two scans
    assert 0 < ms2_peaks(ms2_centroid_algo=HILL, ms2_centroid_mz_ppm=5.0) < hill
    # no tolerance: hill centroiding is refused, the spectra stay as they are
    assert ms2_peaks(ms2_centroid_algo=HILL, ms2_centroid_mz_ppm=0.0) == raw
    # the isotopic prefilter applies to DIA MS2 with and without centroiding
    assert 0 < ms2_peaks(isotopic_prefilter=True) < raw
    assert 0 < ms2_peaks(ms2_centroid_algo=HILL, isotopic_prefilter=True) < hill
