"""DataFrame export utilities for pyOpenMS.

This module provides utility functions for converting OpenMS data structures to pandas DataFrames.
The get_df() methods are now implemented directly in the Cython classes (MSSpectrum, MSChromatogram,
Mobilogram, MSExperiment, ConsensusMap, FeatureMap, MRMTransitionGroupCP, PeptideIdentificationList).

This module provides backwards-compatible function aliases:
- peptide_identifications_to_df: Calls PeptideIdentificationList.get_df()
- update_scores_from_df: Calls PeptideIdentificationList.update_scores_from_df()

And utility functions:
- common_meta_value_types: Dictionary for numpy type mapping of common meta values
"""
from __future__ import annotations

__all__ = [
    'peptide_identifications_to_df',
    'update_scores_from_df',
    'common_meta_value_types',
    '_add_meta_values',
]

from typing import Any, TYPE_CHECKING

from . import PeptideIdentificationList as _PeptideIdentificationList
from . import DataValue as _DataValue

import numpy as _np

if TYPE_CHECKING:
    import pandas as _pd
else:
    class _PandasStub:
        DataFrame = Any

    _pd = _PandasStub()


# Common meta value types for numpy type mapping
# String types use 'object' dtype - memory efficient and safe (no truncation risk)
common_meta_value_types = {
    b'label': 'object',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'object',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f',
    b'Group': 'object',
    b'is_ungrouped_monoisotopic': 'i',
    b'leftWidth': 'f',
    b'rightWidth': 'f',
    b'total_xic': 'f',
    b'PeptideRef': 'object',
    b'peak_apices_sum': 'f'
}
"""Global dict to define which autoconversion to numpy types is tried for certain metavalues.

This can be changed to your liking but only affects future exports of any OpenMS datastructure to dataframes.
String types use 'object' dtype for memory efficiency (no fixed-width waste) and safety (no truncation).
"""


def peptide_identifications_to_df(peps: _PeptideIdentificationList, decode_ontology: bool = True,
                                  default_missing_values: dict = None,
                                  export_unidentified: bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    This is a backwards-compatible wrapper that calls PeptideIdentificationList.get_df().
    For new code, prefer calling peps.get_df() directly.

    :param peps: list of PeptideIdentification objects
    :type peps: PeptideIdentificationList
    :param decode_ontology: decode meta value names
    :type decode_ontology: bool
    :param default_missing_values: default value for missing values for each data type
    :type default_missing_values: dict
    :param export_unidentified: export PeptideIdentifications without PeptideHit
    :type export_unidentified: bool
    :return: peptide identifications in a DataFrame
    :rtype: pandas.DataFrame
    """
    return peps.get_df(decode_ontology=decode_ontology,
                       default_missing_values=default_missing_values,
                       export_unidentified=export_unidentified)


def update_scores_from_df(peps: _PeptideIdentificationList, df: _pd.DataFrame, main_score_name: str):
    """
    Updates the scores in PeptideIdentification objects using a pandas dataframe.

    This is a backwards-compatible wrapper that calls PeptideIdentificationList.update_scores_from_df().
    For new code, prefer calling peps.update_scores_from_df(df, main_score_name) directly.

    :param peps: list of PeptideIdentification objects
    :param df: pandas dataframe obtained by converting peps to a dataframe. Minimum required: P_ID column and column with name passed by main_score_name
    :param main_score_name: name of the score column
    :return: the updated list of peptide identifications
    """
    return peps.update_scores_from_df(df, main_score_name)


def _add_meta_values(df: _pd.DataFrame, object: Any) -> _pd.DataFrame:
    """
    Adds metavalues from given object to given DataFrame.

    :param df: DataFrame to which metavalues will be added.
    :type df: pandas.DataFrame
    :param object: Object from which metavalues will be extracted.
    :type object: Any
    :return: DataFrame with added meta values.
    :rtype: pandas.DataFrame
    """
    mvs = []
    object.getKeys(mvs)

    for k in mvs:
        dv = object.getMetaValue(k)
        col_name = k.decode()

        try:
            # Handle native Python types (returned by autowrap)
            if isinstance(dv, float):
                value = dv
                dtype = "float64"
            elif isinstance(dv, int):
                value = dv
                dtype = "int64"
            elif isinstance(dv, str):
                value = dv
                dtype = "object"
            elif isinstance(dv, bytes):
                value = dv.decode()
                dtype = "object"
            # Handle DataValue objects (if ever returned)
            elif hasattr(dv, 'valueType'):
                if dv.valueType() == _DataValue.STRING_VALUE:
                    value = dv.toString().decode()
                    dtype = "object"
                elif dv.valueType() == _DataValue.INT_VALUE:
                    value = dv.toInt()
                    dtype = "int32"
                elif dv.valueType() == _DataValue.DOUBLE_VALUE:
                    value = dv.toDouble()
                    dtype = "float64"
                elif dv.valueType() == _DataValue.EMPTY_VALUE:
                    continue
                else:
                    value = str(dv)
                    dtype = "object"
            else:
                value = str(dv)
                dtype = "object"

            df[col_name] = _np.full(df.shape[0], value, dtype=dtype)

        except Exception:
            df[col_name] = _np.full(df.shape[0], str(dv), dtype='object')

    return df
