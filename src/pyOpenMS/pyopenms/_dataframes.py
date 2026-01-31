"""DataFrame export utilities for pyOpenMS.

This module provides utility functions for converting OpenMS data structures to pandas DataFrames.
The to_df() methods are now implemented directly in the Cython classes (MSSpectrum, MSChromatogram,
Mobilogram, MSExperiment, ConsensusMap, FeatureMap, MRMTransitionGroupCP, PeptideIdentificationList).

This module provides backwards-compatible function aliases:
- peptide_identifications_to_df: Calls PeptideIdentificationList.to_df()
- update_scores_from_df: Calls PeptideIdentificationList.update_scores_from_df()
"""
from __future__ import annotations

__all__ = [
    'peptide_identifications_to_df',
    'update_scores_from_df',
]

from typing import Any, TYPE_CHECKING

from . import PeptideIdentificationList as _PeptideIdentificationList

if TYPE_CHECKING:
    import pandas as _pd
else:
    class _PandasStub:
        DataFrame = Any

    _pd = _PandasStub()


def peptide_identifications_to_df(peps: _PeptideIdentificationList, decode_ontology: bool = True,
                                  default_missing_values: dict = None,
                                  export_unidentified: bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    .. deprecated::
        Use ``peps.to_df()`` instead.

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
    import warnings
    warnings.warn(
        "peptide_identifications_to_df() is deprecated and will be removed in a future version. "
        "Use peps.to_df() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return peps.to_df(decode_ontology=decode_ontology,
                      default_missing_values=default_missing_values,
                      export_unidentified=export_unidentified)


def update_scores_from_df(peps: _PeptideIdentificationList, df: _pd.DataFrame, main_score_name: str):
    """
    Updates the scores in PeptideIdentification objects using a pandas dataframe.

    .. deprecated::
        Use ``peps.update_scores_from_df(df, main_score_name)`` instead.

    :param peps: list of PeptideIdentification objects
    :param df: pandas dataframe obtained by converting peps to a dataframe. Minimum required: P_ID column and column with name passed by main_score_name
    :param main_score_name: name of the score column
    :return: the updated list of peptide identifications
    """
    import warnings
    warnings.warn(
        "update_scores_from_df() is deprecated and will be removed in a future version. "
        "Use peps.update_scores_from_df(df, main_score_name) instead.",
        DeprecationWarning,
        stacklevel=2
    )
    return peps.update_scores_from_df(df, main_score_name)
