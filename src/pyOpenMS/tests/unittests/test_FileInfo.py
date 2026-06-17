"""Tests for the library-level FileInfo helper (OpenMS#9596)."""

import os
import shutil
import tempfile

import pyopenms as oms

_HERE = os.path.dirname(os.path.abspath(__file__))


def _path(name):
    return os.path.join(_HERE, name)


def test_fileinfo_featurexml():
    r = oms.FileInfo().run_all(_path("BSA1_F1_idmapped.featureXML"))
    assert r.meta.file_type_name == "featureXML"
    assert r.feature is not None
    assert r.feature.is_consensus is False
    assert r.feature.num_features > 0
    # ranges are populated for feature maps
    assert r.ranges.is_experiment is False
    # the CLI text / tsv are reproduced verbatim and available from Python
    assert "-- General information --" in r.text
    assert "general: number of features" in r.tsv


def test_fileinfo_consensusxml():
    r = oms.FileInfo().run_all(_path("ProteinQuantifier_input.consensusXML"))
    assert r.meta.file_type_name == "consensusXML"
    assert r.feature is not None
    assert r.feature.is_consensus is True
    assert r.feature.num_features > 0
    # consensus map column headers are exposed as a typed list
    assert len(r.feature.map_columns) >= 1
    assert isinstance(r.feature.map_columns[0].filename, str)


def test_fileinfo_mzml_peaks():
    r = oms.FileInfo().run_all(_path("test.mzML"))
    assert r.meta.file_type_name == "mzML"
    assert r.peak is not None
    assert r.peak.num_spectra > 0
    assert len(r.peak.ms_levels) >= 1
    # spectra-per-MS-level is a typed dict
    assert sum(r.peak.spectra_per_ms_level.values()) == r.peak.num_spectra
    # MSExperiment ranges expose the combined / per-level structure
    assert r.ranges.is_experiment is True
    # activation methods come back as (ms_level, name, count) tuples
    for ms_level, name, count in r.peak.activation_methods_flat():
        assert isinstance(ms_level, int)
        assert isinstance(name, str)
        assert count >= 1


def test_fileinfo_idxml():
    r = oms.FileInfo().run_all(_path("test.idXML"))
    assert r.meta.file_type_name == "idXML"
    assert r.ident is not None
    assert r.ident.peptide_hits >= 0
    assert r.ident.num_runs >= 1


def test_fileinfo_options_and_renderers():
    fi = oms.FileInfo()

    # default options: meta/processing/statistics off
    r_default = fi.run(_path("test.featureXML"))
    assert "-- Statistics --" not in r_default.text

    # opt into statistics like the CLI -s flag
    opt = oms.FileInfo.Options()
    opt.statistics = True
    r_stats = fi.run(_path("test.featureXML"), opt)
    assert "-- Statistics --" in r_stats.text

    # static renderers return the cached CLI text / tsv
    assert oms.FileInfo.to_text(r_stats) == r_stats.text
    assert oms.FileInfo.to_tsv(r_stats) == r_stats.tsv


def test_fileinfo_forced_type():
    # Test that forced_type bypasses auto-detection by copying a featureXML
    # file to a .tmp extension and verifying it works with forced_type
    # but fails without it
    with tempfile.TemporaryDirectory() as tmpdir:
        # Copy the featureXML file to a neutral extension
        src = _path("test.featureXML")
        tmp_file = os.path.join(tmpdir, "test_data.tmp")
        shutil.copy(src, tmp_file)

        # With forced_type, it should succeed
        opt = oms.FileInfo.Options()
        opt.forced_type = oms.FileTypes.FileType.FEATUREXML
        r = oms.FileInfo().run(tmp_file, opt)
        assert r.meta.file_type == oms.FileTypes.FileType.FEATUREXML
        assert r.feature is not None

        # Without forced_type, it should fail to recognize the type
        # (the file type will be UNKNOWN)
        r_no_force = oms.FileInfo().run(tmp_file)
        assert r_no_force.meta.file_type == oms.FileTypes.FileType.UNKNOWN
