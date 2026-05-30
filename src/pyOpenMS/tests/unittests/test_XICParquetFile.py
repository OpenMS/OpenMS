# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause

import os

import pytest


def _repo_root():
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
    )


def _xic_path():
    return os.path.join(
        _repo_root(),
        "src",
        "tests",
        "class_tests",
        "openms",
        "data",
        "XICParquetFile_1_input.xic",
    )


def _get_xic():
    import pyopenms as poms

    assert hasattr(poms, 'XICParquetFile'), \
        "XICParquetFile missing from pyopenms — Arrow/Parquet bindings are required"

    return poms.XICParquetFile(_xic_path())


def test_xic_df_columns():
    xic = _get_xic()
    cols = xic.df_columns()
    assert "RUN_ID" in cols
    assert "PRECURSOR_ID" in cols
    assert "RT_DATA" in cols


def test_xic_run_dict():
    xic = _get_xic()
    runs = xic.get_run_dict()
    assert "run_id" in runs
    assert "source_file" in runs
    assert len(runs["run_id"]) > 0


def test_xic_analyte_dict_default_nested():
    xic = _get_xic()
    nested = xic.get_analyte_dict()
    exploded = xic.get_analyte_dict(nest_transitions=False)
    assert len(nested["precursor_id"]) > 0
    assert len(exploded["precursor_id"]) >= len(nested["precursor_id"])
    # nested transition fields should be list-like (list or array)
    assert hasattr(nested["transition_id"][0], "__iter__")


def test_xic_analyte_columns_subset():
    xic = _get_xic()
    subset = xic.get_analyte_dict(columns=["PRECURSOR_ID", "MODIFIED_SEQUENCE"])
    assert set(subset.keys()) == {"precursor_id", "modified_sequence"}


def test_xic_analyte_columns_invalid():
    xic = _get_xic()
    with pytest.raises(RuntimeError):
        xic.get_analyte_dict(columns=["PRECURSOR_CHARGE5"])


def test_xic_to_df_summary_and_exploded():
    xic = _get_xic()
    df_summary = xic.to_df()
    df_exploded = xic.to_df(explode=True)
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert hasattr(df_summary.iloc[0]["rt"], "__iter__")
    assert hasattr(df_summary.iloc[0]["intensity"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["rt"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["intensity"], "__iter__")


def test_xic_query_builder_summary_and_exploded():
    xic = _get_xic()
    analytes = xic.get_analyte_dict(columns=["PRECURSOR_ID"], nest_transitions=False)
    precursor_id = analytes["precursor_id"][0]
    df_summary = xic.query_chromatograms().filter_precursor_id(precursor_id).to_df(explode=False)
    df_exploded = xic.query_chromatograms().filter_precursor_id(precursor_id).to_df(explode=True)
    assert len(df_summary) > 0
    assert len(df_exploded) >= len(df_summary)
    assert hasattr(df_summary.iloc[0]["rt"], "__iter__")
    assert hasattr(df_summary.iloc[0]["intensity"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["rt"], "__iter__")
    assert not hasattr(df_exploded.iloc[0]["intensity"], "__iter__")


def test_xic_query_builder():
    xic = _get_xic()
    analytes = xic.get_analyte_dict(columns=["PRECURSOR_ID"], nest_transitions=False)
    precursor_id = analytes["precursor_id"][0]
    df = xic.query_chromatograms().filter_precursor_id(precursor_id).to_df(explode=False)
    assert len(df) > 0
