"""Backward-compatible DataFrame utility functions."""

from __future__ import annotations


def peptide_identifications_to_df(peps, decode_ontology=True,
                                   default_missing_values=None,
                                   export_unidentified=True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    Deprecated: Use ``peps.to_df()`` instead.
    """
    import warnings
    warnings.warn(
        "peptide_identifications_to_df() is deprecated. Use peps.to_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return peps.to_df(decode_ontology=decode_ontology,
                      default_missing_values=default_missing_values,
                      export_unidentified=export_unidentified)


def update_scores_from_df(peps, df, main_score_name):
    """Updates scores in PeptideIdentification objects using a pandas DataFrame.

    Deprecated: Use ``peps.update_scores_from_df(df, main_score_name)`` instead.
    """
    import warnings
    warnings.warn(
        "update_scores_from_df() is deprecated. Use peps.update_scores_from_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return peps.update_scores_from_df(df, main_score_name)
