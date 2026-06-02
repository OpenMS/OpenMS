# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

import os

import pytest


def _repo_root():
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
    )


def _xipm_path():
    return os.path.join(
        _repo_root(),
        "src",
        "tests",
        "class_tests",
        "openms",
        "data",
        "XIPMParquetFile_1_input.xipm",
    )


def _get_xipm():
    import pyopenms as poms

    assert hasattr(poms, "XIPMParquetFile"), \
        "XIPMParquetFile missing from pyopenms — Arrow/Parquet bindings are required"

    return poms.XIPMParquetFile(_xipm_path())


def test_xipm_df_columns():
    xipm = _get_xipm()
    cols = xipm.df_columns()
    assert "RUN_ID" in cols
    assert "PEAKMAP_TYPE" in cols
    assert "MZ_DATA" in cols


def test_xipm_run_dict():
    xipm = _get_xipm()
    runs = xipm.get_run_dict()
    assert "run_id" in runs
    assert "source_file" in runs
    assert len(runs["run_id"]) > 0


def test_xipm_analyte_dict_default_nested():
    xipm = _get_xipm()
    nested = xipm.get_analyte_dict()
    exploded = xipm.get_analyte_dict(nest_transitions=False)
    assert len(nested["precursor_id"]) > 0
    assert len(exploded["precursor_id"]) >= len(nested["precursor_id"])
    assert hasattr(nested["transition_id"][0], "__iter__")


def test_xipm_analyte_columns_subset():
    xipm = _get_xipm()
    subset = xipm.get_analyte_dict(columns=["PRECURSOR_ID", "MODIFIED_SEQUENCE"])
    assert set(subset.keys()) == {"precursor_id", "modified_sequence"}


def test_xipm_analyte_columns_invalid():
    xipm = _get_xipm()
    with pytest.raises(RuntimeError):
        xipm.get_analyte_dict(columns=["PRECURSOR_CHARGE5"])


def test_xipm_to_df_summary_and_exploded():
    xipm = _get_xipm()
    df_summary = xipm.to_df()
    df_exploded = xipm.to_df(explode=True)
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert "target_rt" in df_summary.columns
    assert "target_rt" in df_exploded.columns
    assert hasattr(df_summary.iloc[0]["mz"], "__iter__")
    assert hasattr(df_summary.iloc[0]["rt"], "__iter__")
    assert hasattr(df_summary.iloc[0]["ion_mobility"], "__iter__")
    assert hasattr(df_summary.iloc[0]["intensity"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["mz"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["rt"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["ion_mobility"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["intensity"], "__iter__")


def test_xipm_query_builder_summary_and_exploded():
    xipm = _get_xipm()
    analytes = xipm.get_analyte_dict(columns=["PRECURSOR_ID"], nest_transitions=False)
    precursor_id = next((pid for pid in analytes["precursor_id"] if pid is not None), None)
    assert precursor_id is not None
    df_summary = (
        xipm.query_peak_maps()
        .filter_precursor_id(precursor_id)
        .filter_peakmap_type("transition")
        .to_df(explode=False)
    )
    df_exploded = (
        xipm.query_peak_maps()
        .filter_precursor_id(precursor_id)
        .filter_peakmap_type("transition")
        .to_df(explode=True)
    )
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert hasattr(df_summary.iloc[0]["mz"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["mz"], "__iter__")


def test_xipm_query_builder():
    xipm = _get_xipm()
    df = xipm.query_peak_maps().filter_peakmap_type("precursor").to_df(explode=False)
    assert len(df) > 0
