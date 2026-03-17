"""Shared helpers for Arrow-based feature/PSM export."""

from __future__ import annotations


def join_best_psm_columns(features_table, psm_table, join_key, output_columns):
    """Join best-hit PSM columns onto a features Arrow table.

    Filters PSMs to rank-0 (best hit), deduplicates by feature ID, and appends
    the requested columns from the PSM table onto the features table.

    Parameters
    ----------
    features_table : pyarrow.Table
        Features table with a ``feature_id`` column.
    psm_table : pyarrow.Table
        PSM table with a join key column and a ``rank`` column.
    join_key : str
        Column name in ``psm_table`` that maps to ``feature_id`` in
        ``features_table`` (e.g. ``'feature_unique_id'``).
    output_columns : dict
        Mapping of ``{psm_column_name: output_column_name}``. Each PSM column
        is looked up and appended to the features table under the output name.
        Values of ``None`` use the PSM column name as-is.

    Returns
    -------
    pyarrow.Table
        Features table with additional columns from best-hit PSMs.
    """
    import pyarrow as pa
    import pyarrow.compute as pc

    if psm_table.num_rows == 0:
        # Append null columns
        n = features_table.num_rows
        for psm_col, out_name in output_columns.items():
            col_name = out_name or psm_col
            features_table = features_table.append_column(
                col_name, pa.nulls(n, type=pa.utf8()))
        return features_table

    # Keep only rank-0 (best) hits
    mask = pc.equal(psm_table.column('rank'), 0)
    best = psm_table.filter(mask)

    # Deduplicate: keep first PSM per feature
    seen = set()
    keep = []
    fids = best.column(join_key).to_pylist()
    for i, fid in enumerate(fids):
        if fid is not None and fid not in seen:
            seen.add(fid)
            keep.append(i)
    best = best.take(keep)

    # Build lookup: feature_id -> {col: value}
    lookup = {}
    psm_col_names = list(output_columns.keys())
    col_lists = {col: best.column(col).to_pylist() for col in psm_col_names}
    fid_list = best.column(join_key).to_pylist()
    for i, fid in enumerate(fid_list):
        lookup[fid] = {col: col_lists[col][i] for col in psm_col_names}

    # Map onto features table
    feature_ids = features_table.column('feature_id').to_pylist()
    for psm_col, out_name in output_columns.items():
        col_name = out_name or psm_col
        # Infer Arrow type from the PSM column
        psm_arrow_type = best.schema.field(psm_col).type
        values = [lookup.get(fid, {}).get(psm_col) for fid in feature_ids]
        features_table = features_table.append_column(
            col_name, pa.array(values, type=psm_arrow_type))

    return features_table
