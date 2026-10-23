"""IdentificationData addon: the matches of a dataset as an Arrow table, and edits of them as patches."""
from __future__ import annotations

from . import addon

# The key columns of a match table, kept by every column selection.
_KEYS = ("match_id", "query_id", "run_uuid")


def _zerocopy():
    try:
        import pyopenms._arrow_zerocopy as zerocopy
    except ImportError as error:  # pragma: no cover - part of every build with Arrow
        raise ImportError("pyopenms._arrow_zerocopy is not available; Arrow support is required") from error
    return zerocopy


@addon("IdentificationData")
def to_arrow(self, columns=None):
    """
    One row per match, as a ``pyarrow.Table`` with the columns of ``matches.parquet`` in a native bundle.

    The columns are ``match_id``, ``query_id``, the molecule (``representation``, ``encoding``, ``charge``,
    ``calculated_mz``, ``target_decoy``, ``name``, ``formula``, ``identifiers``, ``adduct``),
    ``sequence_evidence``, ``peak_annotations``, ``metadata``, one ``score_<name>`` column per score
    definition (missing values are null), and ``run_uuid`` (dictionary-encoded). Runs follow each other,
    and the matches of a run follow their queries.

    The schema metadata describes the table as JSON: ``openms:score_definitions``,
    ``openms:metadata_descriptors`` (which the ``metadata`` column refers to by index),
    ``openms:primary_score`` (index or null) and ``openms:revisions`` (run UUID -> ``Run.getRevision()``).

    The table is a copy: editing it does not edit the dataset. Use :meth:`apply_patch` to write scores
    and target/decoy states back.

    Parameters
    ----------
    columns : list of str, optional
        Columns to keep besides the keys ``match_id``, ``query_id`` and ``run_uuid``; all if not given.

    Returns
    -------
    pyarrow.Table
    """
    table = _zerocopy().identification_data_matches_to_arrow(self)
    if columns is None:
        return table
    unknown = [name for name in columns if name not in table.column_names]
    if unknown:
        raise KeyError(f"Unknown match columns {unknown}; the table has {table.column_names}")
    wanted = set(columns) | set(_KEYS)
    return table.select([name for name in table.column_names if name in wanted])


@addon("IdentificationData")
def __arrow_c_stream__(self, requested_schema=None):
    """
    The Arrow stream of :meth:`to_arrow` (Arrow PyCapsule interface), so that Arrow consumers read the
    matches directly, e.g. ``pyarrow.table(data)``, ``polars.from_arrow(data)`` or DuckDB.
    """
    return self.to_arrow().__arrow_c_stream__(requested_schema)


@addon("IdentificationData")
def apply_patch(self, patch, add_scores=None, expected_revisions=None):
    """
    Edit matches from a table keyed by ``run_uuid`` and ``match_id``: all edits are applied, or none.

    The other columns of the patch are score columns, named as in :meth:`to_arrow` (null or NaN removes a
    supplementary value; the primary score cannot be removed), and ``target_decoy`` (0 unknown, 1 target,
    2 decoy, 3 both). Every key and value is checked before anything changes. Patched runs start a new
    state: match and score views bound before must be taken again.

    Parameters
    ----------
    patch : pyarrow.Table, pandas.DataFrame or any object with ``__arrow_c_stream__``
        The rows to edit. A pandas index is not part of the patch.
    add_scores : list of IdentificationData.ScoreDefinition, optional
        Score definitions added to every run with scores (not to sequence catalogs) before the patch is
        applied, e.g. after rescoring; existing definitions are reused. Their columns are named as in
        :meth:`to_arrow` (``score_<name>``).
    expected_revisions : dict of str to int, optional
        Revision that a run must still have, by run UUID (see :meth:`revisions`).

    Returns
    -------
    dict
        ``rows`` (matches edited), ``columns`` (patched columns) and ``added`` (score columns added).

    Raises
    ------
    ValueError
        For an invalid patch: an unknown column, run or match, a repeated key, a value that is not finite,
        a missing primary score, a target_decoy out of range or a changed revision. Nothing changes then.
    """
    module = type(patch).__module__ or ""
    if module.split(".")[0] == "pandas":
        import pyarrow as pa

        # The index of a frame is not a column of the patch (pandas would export it once it is not a range).
        patch = pa.Table.from_pandas(patch, preserve_index=False)
    return _zerocopy().identification_data_apply_patch(self, patch, list(add_scores or []), dict(expected_revisions or {}))
