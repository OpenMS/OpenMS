"""IdentificationDataFile addon: a flat PSM table from a native identification bundle."""
from __future__ import annotations

import json
import os

from . import addon

# Text of IdentificationData.TargetDecoy, by value (UNKNOWN, TARGET, DECOY, BOTH), as idXML writes it.
_TARGET_DECOY = ["", "target", "decoy", "target+decoy"]


def _string(column):
    """Decode a (dictionary-encoded) run UUID column to plain strings, for joins."""
    import pyarrow as pa
    import pyarrow.compute as pc

    return pc.cast(column, pa.string())


def _accessions(evidence):
    """The accessions of a sequence_evidence column (list of structs), as a list<string> column."""
    import pyarrow as pa
    import pyarrow.compute as pc

    chunks = [pa.ListArray.from_arrays(chunk.offsets, pc.struct_field(chunk.values, "accession")) for chunk in evidence.chunks]
    return pa.chunked_array(chunks, type=pa.list_(pa.string()))


@addon("IdentificationDataFile", "psm_table")
@staticmethod
def psm_table(path, runs=None):
    """
    One row per candidate of a native identification bundle (``.idparquet``), with the values of its query.

    This is the flat PSM table that ``psms.parquet`` of the four-table ``.idparquet`` layout of OpenMS 3.6
    provided. It reads the bundle's tables with pyarrow and does not load the identification data into
    pyOpenMS. Rows follow ``matches.parquet``; queries without candidates have no row.

    Columns:

    - ``run_uuid``, ``run_identifier``: the run of the candidate (UUID and display identifier)
    - ``reference_file_name``: path of the query's source file (empty if the file is not known)
    - ``query_id``, ``spectrum_reference`` (the query's data ID), ``rt``, ``observed_mz`` (the query's m/z)
    - ``match_id``, ``peptidoform`` (the molecular representation), ``encoding``, ``precursor_charge``,
      ``calculated_mz``
    - ``target_decoy`` (``"target"``, ``"decoy"``, ``"target+decoy"`` or ``""``) and ``is_decoy``
    - ``protein_accessions``: the accessions of the candidate's sequence evidence
    - ``score``, ``score_type``, ``higher_score_better``: the primary score of the dataset
    - every score column of ``matches.parquet`` (``score_<name>``), unchanged

    Parameters
    ----------
    path : str
        Directory of the bundle.
    runs : list of str, optional
        Runs to read, by UUID or display identifier; all runs if not given.

    Returns
    -------
    pyarrow.Table
        Use ``.to_pandas()`` for a DataFrame.
    """
    try:
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError("pyarrow is required for psm_table(). Install with: pip install pyarrow")

    manifest_path = os.path.join(path, "manifest.json")
    if not os.path.isfile(manifest_path):
        raise ValueError(
            f"{path} is not a native identification bundle (no manifest.json). The four-table .idparquet layout "
            "of OpenMS 3.6 is not supported; convert it to idXML with IDFileConverter of OpenMS 3.6."
        )
    with open(manifest_path, encoding="utf-8") as stream:
        manifest = json.load(stream)

    descriptors = manifest["runs"]
    if runs is not None:
        wanted = set(runs)
        descriptors = [run for run in descriptors if run["uuid"] in wanted or run["identifier"] in wanted]
        missing = wanted - {run["uuid"] for run in descriptors} - {run["identifier"] for run in descriptors}
        if missing:
            raise KeyError(f"Unknown runs: {sorted(missing)}")
    uuids = [run["uuid"] for run in descriptors]
    score_columns = list(manifest.get("score_columns", []))

    # The primary score is the same column in every run with scores (one score schema per dataset).
    primary = next((run for run in descriptors if run.get("primary_score") is not None), None)
    primary_column = score_columns[primary["primary_score"]] if primary else None
    primary_definition = primary["scores"][primary["primary_score"]] if primary else None

    match_columns = ["run_uuid", "match_id", "query_id", "representation", "encoding", "charge", "calculated_mz",
                     "target_decoy", "sequence_evidence"] + score_columns
    selection = [("run_uuid", "in", uuids)] if runs is not None else None
    matches = pq.read_table(os.path.join(path, "matches.parquet"), columns=match_columns, filters=selection)
    queries = pq.read_table(os.path.join(path, "queries.parquet"),
                            columns=["run_uuid", "query_id", "source_id", "data_id", "rt", "mz"], filters=selection)

    # Join flat keys only (joins do not take nested columns and do not keep the row order), then take rows.
    rows = pa.array(range(matches.num_rows), pa.int64())
    run_uuid = _string(matches["run_uuid"])
    keys = pa.table({"run_uuid": run_uuid, "query_id": matches["query_id"], "__match": rows})
    query_keys = pa.table({"run_uuid": _string(queries["run_uuid"]), "query_id": queries["query_id"],
                           "__query": pa.array(range(queries.num_rows), pa.int64())})
    found = keys.join(query_keys, keys=["run_uuid", "query_id"], join_type="left outer").sort_by("__match")
    query = queries.take(found["__query"])

    runs_table = pa.table({
        "run_uuid": pa.array(uuids, pa.string()),
        "run_identifier": pa.array([run["identifier"] for run in descriptors], pa.string()),
    })
    sources = [(run["uuid"], index, source.get("path", "")) for run in descriptors for index, source in enumerate(run["sources"])]
    sources_table = pa.table({
        "run_uuid": pa.array([s[0] for s in sources], pa.string()),
        "source_id": pa.array([s[1] for s in sources], pa.uint32()),
        "reference_file_name": pa.array([s[2] for s in sources], pa.string()),
    })
    names = pa.table({"run_uuid": run_uuid, "source_id": query["source_id"], "__match": rows})
    names = names.join(runs_table, keys="run_uuid", join_type="left outer")
    names = names.join(sources_table, keys=["run_uuid", "source_id"], join_type="left outer").sort_by("__match")

    target_decoy = pa.array(_TARGET_DECOY).take(pc.cast(matches["target_decoy"], pa.int64()))
    if primary_column is not None:
        score = matches[primary_column]
        score_type = pa.array([primary_definition["name"]] * matches.num_rows, pa.string())
        higher_better = pa.array([primary_definition["higher_better"]] * matches.num_rows, pa.bool_())
    else:
        score = pa.nulls(matches.num_rows, pa.float64())
        score_type = pa.nulls(matches.num_rows, pa.string())
        higher_better = pa.nulls(matches.num_rows, pa.bool_())

    columns = {
        "run_uuid": run_uuid,
        "run_identifier": names["run_identifier"],
        "reference_file_name": names["reference_file_name"],
        "query_id": matches["query_id"],
        "spectrum_reference": query["data_id"],
        "rt": query["rt"],
        "observed_mz": query["mz"],
        "match_id": matches["match_id"],
        "peptidoform": matches["representation"],
        "encoding": matches["encoding"],
        "precursor_charge": matches["charge"],
        "calculated_mz": matches["calculated_mz"],
        "target_decoy": target_decoy,
        "is_decoy": pc.equal(matches["target_decoy"], 2),
        "protein_accessions": _accessions(matches["sequence_evidence"]),
        "score": score,
        "score_type": score_type,
        "higher_score_better": higher_better,
    }
    for name in score_columns:
        columns[name] = matches[name]
    return pa.table(columns)
