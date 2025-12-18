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

__all__ = [
    'peptide_identifications_to_df',
    'update_scores_from_df',
    'common_meta_value_types',
    '_add_meta_values',
]

from typing import List, Union

from . import PeptideIdentificationList as _PeptideIdentificationList
from . import DataValue as _DataValue

import pandas as _pd
import numpy as _np


# Common meta value types for numpy type mapping
common_meta_value_types = {
    b'label': 'U50',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U100',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f',
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i',
    b'leftWidth': 'f',
    b'rightWidth': 'f',
    b'total_xic': 'f',
    b'PeptideRef': 'U100',
    b'peak_apices_sum': 'f'
}
"""Global dict to define which autoconversion to numpy types is tried for certain metavalues.

This can be changed to your liking but only affects future exports of any OpenMS datastructure to dataframes.
Especially string lengths (i.e., U types) benefit from adaption to save memory. The default type is currently
hardcoded to U50 (i.e., 50 unicode characters)
"""


def peptide_identifications_to_df(peps: _PeptideIdentificationList, decode_ontology: bool = True,
                                  default_missing_values: dict = None,
                                  export_unidentified: bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    This is a backwards-compatible wrapper that calls PeptideIdentificationList.get_df().
    For new code, prefer calling peps.get_df() directly.

    Parameters:
    peps (PeptideIdentificationList): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names
    default_missing_values: default value for missing values for each data type
    export_unidentified: export PeptideIdentifications without PeptideHit

    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
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


def _add_meta_values(df: _pd.DataFrame, object: any) -> _pd.DataFrame:
    """
    Adds metavalues from given object to given DataFrame.

    Args:
        df (pd.DataFrame): DataFrame to which metavalues will be added.
        object (any): Object from which metavalues will be extracted.

    Returns:
        pd.DataFrame: DataFrame with added meta values.
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
                dtype = f"U{max(1, len(value))}"
            elif isinstance(dv, bytes):
                value = dv.decode()
                dtype = f"U{max(1, len(value))}"
            # Handle DataValue objects (if ever returned)
            elif hasattr(dv, 'valueType'):
                if dv.valueType() == _DataValue.STRING_VALUE:
                    value = dv.toString().decode()
                    dtype = f"U{max(1, len(value))}"
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
