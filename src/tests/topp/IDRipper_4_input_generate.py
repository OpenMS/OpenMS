#!/usr/bin/env python3
"""Regenerate IDRipper_4_input.idparquet from scratch.

Produces a minimal merged .idparquet directory whose 2 ProteinIdentifications
share stored run_identifier "TopPerc_test" — the post-Percolator-merged shape
described in the validation memory entry [[idparquet-quant-chain]]. The PSMs
split into 2 PeptideIdentifications, each carrying a `file_origin`
spectrum_metavalue, so IDRipper can rip them into per-file outputs.

Requires only pyarrow (no OpenMS build).
"""
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq

OUT_DIR = Path(__file__).parent / "IDRipper_4_input.idparquet"
RUN_ID = "TopPerc_test"

# === Schemas (mirror the OpenMS ArrowSchemaRegistry on develop) ===

MV_STRUCT = pa.struct([("name", pa.string()), ("value", pa.string()), ("value_type", pa.string())])
MV_LIST = pa.list_(MV_STRUCT)

# Non-null modifier: matches OpenMS schema registry.
def nn(field_name, dtype):
    return pa.field(field_name, dtype, nullable=False)

SEARCH_PARAMS_SCHEMA = pa.schema([
    nn("run_identifier", pa.string()),
    nn("search_engine", pa.string()),
    ("search_engine_version", pa.string()),
    ("inference_engine", pa.string()),
    ("inference_engine_version", pa.string()),
    ("date", pa.timestamp("ms")),
    nn("score_type", pa.string()),
    nn("higher_score_better", pa.bool_()),
    ("significance_threshold", pa.float64()),
    ("db", pa.string()),
    ("db_version", pa.string()),
    ("taxonomy", pa.string()),
    ("charges", pa.string()),
    nn("mass_type", pa.string()),
    nn("precursor_mass_tolerance", pa.float64()),
    nn("precursor_mass_tolerance_ppm", pa.bool_()),
    nn("fragment_mass_tolerance", pa.float64()),
    nn("fragment_mass_tolerance_ppm", pa.bool_()),
    ("digestion_enzyme", pa.string()),
    ("enzyme_term_specificity", pa.string()),
    nn("missed_cleavages", pa.int32()),
    nn("fixed_modifications", pa.list_(pa.string())),
    nn("variable_modifications", pa.list_(pa.string())),
    nn("primary_ms_run_paths", pa.list_(pa.string())),
    nn("metavalues", MV_LIST),
    nn("sp_metavalues", MV_LIST),
])

PROTEINS_SCHEMA = pa.schema([
    nn("accession", pa.string()),
    nn("score", pa.float64()),
    nn("rank", pa.int32()),
    ("coverage", pa.float64()),
    ("sequence", pa.string()),
    ("description", pa.string()),
    ("is_decoy", pa.bool_()),
    nn("run_identifier", pa.string()),
    ("modifications", pa.list_(pa.struct([("position", pa.int32()), ("modification", pa.string())]))),
    nn("metavalues", MV_LIST),
])

PROTEIN_GROUPS_SCHEMA = pa.schema([
    nn("group_type", pa.string()),
    nn("probability", pa.float64()),
    nn("accessions", pa.list_(pa.string())),
    nn("run_identifier", pa.string()),
    nn("group_index", pa.int32()),
    ("float_data", pa.list_(pa.struct([("name", pa.string()), ("values", pa.list_(pa.float64()))]))),
    ("string_data", pa.list_(pa.struct([("name", pa.string()), ("values", pa.list_(pa.string()))]))),
    ("integer_data", pa.list_(pa.struct([("name", pa.string()), ("values", pa.list_(pa.int64()))]))),
])

PSMS_SCHEMA = pa.schema([
    ("sequence", pa.string()),
    ("peptidoform", pa.string()),
    ("modifications", pa.list_(pa.struct([
        ("name", pa.string()), ("accession", pa.string()),
        ("positions", pa.list_(pa.struct([("position", pa.string()), ("scores", pa.float64())])))]))),
    ("precursor_charge", pa.int32()),
    ("posterior_error_probability", pa.float64()),
    ("is_decoy", pa.bool_()),
    ("calculated_mz", pa.float64()),
    ("observed_mz", pa.float64()),
    ("additional_scores", pa.list_(pa.struct([
        ("score_name", pa.string()), ("score_value", pa.float64()), ("higher_better", pa.bool_())]))),
    ("protein_accessions", pa.list_(pa.struct([
        pa.field("accession", pa.string(), nullable=False),
        ("aa_before", pa.string()), ("aa_after", pa.string()),
        ("start", pa.int32()), ("end", pa.int32())]))),
    ("predicted_rt", pa.float64()),
    ("reference_file_name", pa.string()),
    ("cv_params", pa.string()),
    ("scan", pa.int32()),
    ("rt", pa.float64()),
    ("ion_mobility", pa.float64()),
    ("spectrum_reference", pa.string()),
    ("score", pa.float64()),
    ("score_type", pa.string()),
    ("higher_score_better", pa.bool_()),
    ("hit_index", pa.int32()),
    ("peptide_identification_index", pa.int32()),
    ("psm_metavalues", MV_LIST),
    ("spectrum_metavalues", MV_LIST),
    ("run_identifier", pa.string()),
    ("mz_array", pa.list_(pa.float32())),
    ("intensity_array", pa.list_(pa.float32())),
    ("charge_array", pa.list_(pa.int32())),
    ("ion_type_array", pa.list_(pa.string())),
])

# === Build tables ===

def write(table, name):
    OUT_DIR.mkdir(exist_ok=True)
    # Minimal QPX metadata so OpenMS doesn't reject the file.
    md = {b"qpx_version": b"1.0", b"creator": b"OpenMS",
          b"file_type": name.encode(), b"software_provider": b"OpenMS",
          b"creation_date": b"2026-05-23T00:00:00Z",
          b"uuid": b"00000000-0000-4000-8000-000000000000"}
    table = table.replace_schema_metadata({**(table.schema.metadata or {}), **md})
    pq.write_table(table, OUT_DIR / f"{name}.parquet")

# search_params: 1 row
sp = pa.Table.from_pylist([{
    "run_identifier": RUN_ID, "search_engine": "Comet", "search_engine_version": None,
    "inference_engine": None, "inference_engine_version": None,
    "date": None, "score_type": "q-value", "higher_score_better": False,
    "significance_threshold": 0.0, "db": None, "db_version": None, "taxonomy": None,
    "charges": None, "mass_type": "MONOISOTOPIC",
    "precursor_mass_tolerance": 0.0, "precursor_mass_tolerance_ppm": False,
    "fragment_mass_tolerance": 0.0, "fragment_mass_tolerance_ppm": False,
    "digestion_enzyme": None, "enzyme_term_specificity": None, "missed_cleavages": 0,
    "fixed_modifications": [], "variable_modifications": [],
    "primary_ms_run_paths": [], "metavalues": [], "sp_metavalues": [],
}], schema=SEARCH_PARAMS_SCHEMA)
write(sp, "search_params")

# proteins: 2 rows (same run_identifier as the only IdentificationRun)
prot = pa.Table.from_pylist([
    {"accession": "PROT_A", "score": 0.01, "rank": 0, "coverage": None, "sequence": None,
     "description": None, "is_decoy": False, "run_identifier": RUN_ID,
     "modifications": None, "metavalues": []},
    {"accession": "PROT_B", "score": 0.02, "rank": 0, "coverage": None, "sequence": None,
     "description": None, "is_decoy": False, "run_identifier": RUN_ID,
     "modifications": None, "metavalues": []},
], schema=PROTEINS_SCHEMA)
write(prot, "proteins")

# protein_groups: empty
pg = pa.Table.from_pylist([], schema=PROTEIN_GROUPS_SCHEMA)
write(pg, "protein_groups")

# psms: 4 PSMs in 2 PeptideIdentifications, split by file_origin spectrum-metavalue.
def file_origin_mv(name):
    return [{"name": "file_origin", "value": name, "value_type": "string"}]

psm_rows = []
for pid_idx, (origin, base) in enumerate([("fileA.idXML", "AAAA"), ("fileB.idXML", "BBBB")]):
    for hit in range(2):
        # Suffix with all-letter tag (R/K) to keep sequences valid AAs.
        tag = ("R" * (hit + 1))  # "R" / "RR" instead of "0"/"1"
        psm_rows.append({
            "sequence": f"{base}PEPTIDE{tag}", "peptidoform": f"{base}PEPTIDE{tag}",
            "modifications": [], "precursor_charge": 2,
            "posterior_error_probability": None, "is_decoy": False,
            "calculated_mz": 500.0 + pid_idx + 0.1 * hit, "observed_mz": 500.0,
            "additional_scores": [], "protein_accessions": [
                {"accession": "PROT_A" if pid_idx == 0 else "PROT_B",
                 "aa_before": None, "aa_after": None, "start": None, "end": None}],
            "predicted_rt": None, "reference_file_name": None, "cv_params": None,
            "scan": None, "rt": 100.0 + pid_idx * 10, "ion_mobility": None,
            "spectrum_reference": None, "score": 0.001 + hit * 0.001, "score_type": "q-value",
            "higher_score_better": False, "hit_index": hit,
            "peptide_identification_index": pid_idx,
            "psm_metavalues": [], "spectrum_metavalues": file_origin_mv(origin),
            "run_identifier": RUN_ID,
            "mz_array": None, "intensity_array": None, "charge_array": None,
            "ion_type_array": None,
        })

psms = pa.Table.from_pylist(psm_rows, schema=PSMS_SCHEMA)
write(psms, "psms")

print(f"Generated {OUT_DIR}")
