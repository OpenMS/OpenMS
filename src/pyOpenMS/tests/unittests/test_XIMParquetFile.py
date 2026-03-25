# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

import os

import pytest


def _repo_root():
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
    )


def _xim_path():
    return os.path.join(
        _repo_root(),
        "src",
        "tests",
        "class_tests",
        "openms",
        "data",
        "XIMParquetFile_23_input.xim",
    )


def _get_xim():
    import pyopenms as poms

    assert hasattr(poms, 'XIMParquetFile'), \
        "XIMParquetFile missing from pyopenms — Arrow/Parquet bindings are required"

    return poms.XIMParquetFile(_xim_path())


def test_xim_df_columns():
    xim = _get_xim()
    cols = xim.df_columns()
    assert "RUN_ID" in cols
    assert "MOBILITY_DATA" in cols


def test_xim_run_dict():
    xim = _get_xim()
    runs = xim.get_run_dict()
    assert "run_id" in runs
    assert "source_file" in runs
    assert len(runs["run_id"]) > 0


def test_xim_analyte_dict_default_nested():
    xim = _get_xim()
    nested = xim.get_analyte_dict()
    exploded = xim.get_analyte_dict(nest_transitions=False)
    assert len(nested["precursor_id"]) > 0
    assert len(exploded["precursor_id"]) >= len(nested["precursor_id"])
    # nested transition fields should be list-like (list or array)
    assert hasattr(nested["transition_id"][0], "__iter__")


def test_xim_analyte_columns_subset():
    xim = _get_xim()
    subset = xim.get_analyte_dict(columns=["PRECURSOR_ID", "MODIFIED_SEQUENCE"])
    assert set(subset.keys()) == {"precursor_id", "modified_sequence"}


def test_xim_analyte_columns_invalid():
    xim = _get_xim()
    with pytest.raises(RuntimeError):
        xim.get_analyte_dict(columns=["PRECURSOR_CHARGE5"])


def test_xim_to_df_summary_and_exploded():
    xim = _get_xim()
    df_summary = xim.to_df()
    df_exploded = xim.to_df(explode=True)
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert hasattr(df_summary.iloc[0]["mobility"], "__iter__")
    assert hasattr(df_summary.iloc[0]["intensity"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["mobility"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["intensity"], "__iter__")


def test_xim_query_builder_summary_and_exploded():
    xim = _get_xim()
    analytes = xim.get_analyte_dict(columns=["PRECURSOR_ID"], nest_transitions=False)
    precursor_id = next(
        (pid for pid in analytes["precursor_id"] if pid is not None), None
    )
    assert precursor_id is not None
    df_summary = (
        xim.query_mobilograms().filter_precursor_id(precursor_id).to_df(explode=False)
    )
    df_exploded = (
        xim.query_mobilograms().filter_precursor_id(precursor_id).to_df(explode=True)
    )
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert hasattr(df_summary.iloc[0]["mobility"], "__iter__")
    assert hasattr(df_summary.iloc[0]["intensity"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["mobility"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["intensity"], "__iter__")


def test_xim_query_builder():
    xim = _get_xim()
    analytes = xim.get_analyte_dict(columns=["PRECURSOR_ID"], nest_transitions=False)
    precursor_id = next(
        (pid for pid in analytes["precursor_id"] if pid is not None), None
    )
    assert precursor_id is not None
    df = xim.query_mobilograms().filter_precursor_id(precursor_id).to_df(explode=False)
    assert len(df) > 0
