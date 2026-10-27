"""Tests for the BrukerTimsImagingFile loader bindings.

The class is only exposed when OpenMS was built with ``WITH_OPENTIMS=ON``,
so every test feature-detects via ``hasattr(pyopenms, "BrukerTimsImagingFile")``.
The metadata readers and the error paths of load() run on a fake .d folder
(an ``analysis.tdf`` SQLite file without frame data), mirroring the C++
``BrukerTimsImagingFile_test``.
"""
import sqlite3
from pathlib import Path
import pytest
import pyopenms

pytestmark = pytest.mark.skipif(
    not hasattr(pyopenms, "BrukerTimsImagingFile"),
    reason="pyOpenMS built without WITH_OPENTIMS",
)


def _make_fake_d(parent: Path, name: str, maldi_application_type: str,
                 rows: list[tuple[int, int, int, int]],
                 with_z_column: bool = False, spot_size_um: str = "") -> str:
    """Writes ``<parent>/<name>.d/analysis.tdf`` with GlobalMetadata and MaldiFrameInfo.

    Each row in ``rows`` is (frame, x, y, z); z is only written with ``with_z_column``.
    """
    d = parent / f"{name}.d"
    d.mkdir(parents=True)
    conn = sqlite3.connect(str(d / "analysis.tdf"))
    conn.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    if maldi_application_type:
        conn.execute("INSERT INTO GlobalMetadata VALUES ('MaldiApplicationType', ?)", (maldi_application_type,))
    if spot_size_um:
        conn.execute("INSERT INTO GlobalMetadata VALUES ('MaldiSpotSize_um', ?)", (spot_size_um,))
    if with_z_column:
        conn.execute("CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER, ZIndexPos INTEGER)")
        conn.executemany("INSERT INTO MaldiFrameInfo VALUES (?, ?, ?, ?)", rows)
    else:
        conn.execute("CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER)")
        conn.executemany("INSERT INTO MaldiFrameInfo VALUES (?, ?, ?)", [r[:3] for r in rows])
    conn.commit()
    conn.close()
    return str(d)


def _load_both_ways(path: str, match: str,
                    config: "pyopenms.BrukerTimsImagingFile.Config | None" = None) -> None:
    """Calls both load styles; each must fail with the check named by ``match``.

    The fake .d has no frame data, so any load that gets past the check under test
    fails later in the frame reader; matching the message pins the check itself.
    """
    f = pyopenms.BrukerTimsImagingFile()
    args = (config,) if config is not None else ()
    with pytest.raises(Exception, match=match) as returned:
        f.load(path, *args)
    with pytest.raises(Exception, match=match) as filled:
        f.load(path, pyopenms.MSImagingExperiment(), *args)
    assert type(returned.value) is type(filled.value)


def test_config_defaults():
    cfg = pyopenms.BrukerTimsImagingFile.Config()
    assert cfg.strict_imaging_only is True
    assert isinstance(cfg.inner_config, pyopenms.BrukerTimsFile.Config)
    assert cfg.inner_config.export_mode == pyopenms.BrukerTimsFile.Config.ExportMode.AUTO


def test_config_field_assignment():
    cfg = pyopenms.BrukerTimsImagingFile.Config()
    cfg.strict_imaging_only = False
    assert cfg.strict_imaging_only is False

    # the nested config is edited in place
    cfg.inner_config.calibrate = True
    assert cfg.inner_config.calibrate is True

    inner = pyopenms.BrukerTimsFile.Config()
    inner.ms1_centroid_mz_ppm = 5.0
    cfg.inner_config = inner
    assert cfg.inner_config.ms1_centroid_mz_ppm == pytest.approx(5.0)


def test_maldi_pixel_fields():
    px = pyopenms.BrukerTimsImagingFile.MaldiPixel()
    assert (px.frame_id, px.x, px.y) == (0, 0, 0)
    px.frame_id, px.x, px.y = 3, 4, 5
    assert (px.frame_id, px.x, px.y) == (3, 4, 5)
    assert repr(px) == "MaldiPixel(frame_id=3, x=4, y=5)"


def test_is_imaging_dataset(tmp_path):
    d_img = _make_fake_d(tmp_path, "imaging", "Imaging",
                         [(1, 0, 0, 0), (2, 1, 0, 0), (3, 0, 1, 0), (4, 1, 1, 0)])
    assert pyopenms.BrukerTimsImagingFile.isImagingDataset(d_img) is True

    d_dd = _make_fake_d(tmp_path, "droplet", "SingleSpectra", [(1, 0, 0, 0)])
    assert pyopenms.BrukerTimsImagingFile.isImagingDataset(d_dd) is False

    # missing folder: False, not an exception
    assert pyopenms.BrukerTimsImagingFile.isImagingDataset(str(tmp_path / "does_not_exist.d")) is False


def test_read_global_metadata_value(tmp_path):
    d = _make_fake_d(tmp_path, "gm", "Imaging", [(1, 0, 0, 0)], spot_size_um="50")
    read = pyopenms.BrukerTimsImagingFile.readGlobalMetadataValue
    assert read(d, "MaldiApplicationType") == "Imaging"
    assert read(d, "MaldiSpotSize_um") == "50"
    assert read(d, "UnknownKey") == ""
    with pytest.raises(Exception):
        read(str(tmp_path / "nope.d"), "MaldiApplicationType")


def test_read_maldi_frame_info(tmp_path):
    # raw coordinates are returned, not normalized
    d = _make_fake_d(tmp_path, "grid2x2", "Imaging",
                     [(13, 11, 6, 0), (10, 10, 5, 0), (12, 10, 6, 0), (11, 11, 5, 0)])
    rows = pyopenms.BrukerTimsImagingFile.readMaldiFrameInfo(d)
    assert [(r.frame_id, r.x, r.y) for r in rows] == [(10, 10, 5), (11, 11, 5), (12, 10, 6), (13, 11, 6)]

    # a ZIndexPos column is ignored by the 2D reader
    d_z = _make_fake_d(tmp_path, "withz", "Imaging", [(1, 3, 4, 0), (2, 3, 5, 0)], with_z_column=True)
    assert [(r.x, r.y) for r in pyopenms.BrukerTimsImagingFile.readMaldiFrameInfo(d_z)] == [(3, 4), (3, 5)]

    d_empty = _make_fake_d(tmp_path, "empty", "Imaging", [])
    with pytest.raises(Exception):
        pyopenms.BrukerTimsImagingFile.readMaldiFrameInfo(d_empty)


def test_load_nonexistent_path_raises(tmp_path):
    _load_both_ways(str(tmp_path / "nope.d"), "is not readable")
    _load_both_ways(str(tmp_path / "nope.d"), "is not readable", pyopenms.BrukerTimsImagingFile.Config())


def test_load_non_imaging_strict_raises(tmp_path):
    d_dd = _make_fake_d(tmp_path, "rej", "SingleSpectra", [(1, 0, 0, 0), (2, 1, 0, 0)])
    _load_both_ways(d_dd, "Not a MALDI imaging dataset")
    _load_both_ways(d_dd, "Not a MALDI imaging dataset", pyopenms.BrukerTimsImagingFile.Config())


def test_load_multi_z_raises(tmp_path):
    d_multi = _make_fake_d(tmp_path, "multiz", "Imaging",
                           [(1, 0, 0, 0), (2, 1, 0, 0), (3, 0, 0, 1), (4, 1, 0, 1)],
                           with_z_column=True)
    _load_both_ways(d_multi, "Multi-section")
