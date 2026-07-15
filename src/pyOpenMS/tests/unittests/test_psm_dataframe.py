"""
Unit tests for PSM DataFrame export from PeptideIdentificationList.

Tests the PSM export methods:
- to_psm_df(): Export all PSMs as DataFrame
- to_psm_arrow(): Export all PSMs as Arrow Table
- to_parquet(): Export all PSMs to Parquet file
- to_psm_qpx(): Export as QPX format (dict with file_metadata and psms)
- psm_columns(): List available columns

All PSM export methods return results sorted by:
1. rt (retention time)
2. observed_mz (precursor m/z)
3. precursor_charge
4. hit_index (positional, hit index within identification)
"""

import pytest
import numpy as np
import pandas as pd

# Skip entire module if pyarrow is not installed (required for PSM export methods)
pytest.importorskip("pyarrow")


def to_list(val):
    """Convert numpy array to list for comparison. Handles Arrow->pandas conversion."""
    if hasattr(val, 'tolist'):
        return val.tolist()
    return val


def is_null(val):
    """Check if value is null (None or NaN). Handles Arrow->pandas conversion."""
    if val is None:
        return True
    try:
        return pd.isna(val)
    except (TypeError, ValueError):
        return False


def create_test_data():
    """Create test PeptideIdentifications with multiple hits."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Create a PeptideIdentification with multiple hits
    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.5)
    pep_id.setMZ(500.25)
    pep_id.setScoreType("Posterior Error Probability")
    pep_id.setMetaValue("spectrum_reference", "controllerType=0 controllerNumber=1 scan=1234")

    hits = []
    for i, (seq_str, score) in enumerate([
        ("PEPTIDE", 0.01),
        ("PEPTIDER", 0.05),
    ]):
        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq_str))
        hit.setCharge(2)
        hit.setScore(score)
        hit.setMetaValue("target_decoy", "target")
        hit.setMetaValue("some_score", 42.5)

        # Add protein evidence
        ev = oms.PeptideEvidence()
        ev.setProteinAccession(f"PROT{i}")
        hit.setPeptideEvidences([ev])

        hits.append(hit)

    pep_id.setHits(hits)
    pep_ids.append(pep_id)

    # Add second PeptideIdentification with single hit
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    pep_id2.setScoreType("PEP")
    pep_id2.setMetaValue("spectrum_reference", "scan=5678")

    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("TESTPEPTIDE"))
    hit2.setCharge(3)
    hit2.setScore(0.001)
    hit2.setMetaValue("target_decoy", "decoy")
    pep_id2.setHits([hit2])
    pep_ids.append(pep_id2)

    return pep_ids


def test_psm_df_basic():
    """Test basic PSM DataFrame export."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # Should have 3 rows (2 hits from first ID + 1 from second)
    assert len(df) == 3

    # Check required columns exist
    required_cols = [
        "sequence", "peptidoform", "precursor_charge",
        "observed_mz", "rt", "hit_index", "score", "score_type",
        "protein_accessions", "P_ID"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing column: {col}"

    # Check ranks are correct
    assert list(df["hit_index"]) == [0, 1, 0]

    # Check P_ID correctly tracks identification index
    assert list(df["P_ID"]) == [0, 0, 1]


def test_psm_df_top_hit_only():
    """Test exporting only top hit."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df(export_all_hits=False)

    # Should have 2 rows (1 per identification)
    assert len(df) == 2
    assert list(df["hit_index"]) == [0, 0]


def test_psm_df_columns_filter():
    """Test filtering to specific columns."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df(columns=['sequence', 'precursor_charge', 'score'])

    # Should only have requested columns
    assert list(df.columns) == ['sequence', 'precursor_charge', 'score']
    assert len(df) == 3

    # Verify data is correct
    assert df.iloc[0]["sequence"] == "PEPTIDE"
    assert df.iloc[0]["precursor_charge"] == 2


def test_psm_df_scan_parsing():
    """Test scan number extraction from spectrum reference."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # First two PSMs from first identification with scan=1234 (QPX uses string)
    assert df.iloc[0]["scan"] == "1234"
    assert df.iloc[1]["scan"] == "1234"

    # Third PSM from second identification with scan=5678
    assert df.iloc[2]["scan"] == "5678"


def test_psm_df_scan_parsing_native_id_formats():
    """Test scan extraction from various native ID formats using SpectrumLookup."""
    import pyopenms as oms

    test_cases = [
        # (spectrum_reference, expected_scan)
        ("controllerType=0 controllerNumber=1 scan=1234", "1234"),  # Thermo
        ("scan=5678", "5678"),  # Simple scan format
        ("function=1 process=0 scan=999", "999"),  # Waters
        ("index=42", "43"),  # index format (0-based, so scan = index+1)
        ("scanId=100", "100"),  # Agilent MassHunter
        ("spectrum=200", "200"),  # spectrum format
    ]

    for spec_ref, expected_scan in test_cases:
        pep_ids = oms.PeptideIdentificationList()
        pep_id = oms.PeptideIdentification()
        pep_id.setMetaValue("spectrum_reference", spec_ref)

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
        hit.setCharge(2)
        pep_id.setHits([hit])
        pep_ids.append(pep_id)

        df = pep_ids.to_psm_df()
        assert df.iloc[0]["scan"] == expected_scan, f"Failed for '{spec_ref}': expected {expected_scan}, got {df.iloc[0]['scan']}"


def test_psm_df_decoy_detection():
    """Test is_decoy field extraction (QPX uses int: 0=target, 1=decoy)."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # First two are targets (0)
    assert df.iloc[0]["is_decoy"] == 0
    assert df.iloc[1]["is_decoy"] == 0

    # Third is decoy (1)
    assert df.iloc[2]["is_decoy"] == 1


def test_psm_df_pep_score():
    """Test posterior_error_probability extraction."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # All should have PEP since score types contain "PEP"
    assert "posterior_error_probability" in df.columns

    # Check values match scores
    assert df.iloc[0]["posterior_error_probability"] == pytest.approx(0.01)
    assert df.iloc[2]["posterior_error_probability"] == pytest.approx(0.001)


def test_psm_df_protein_accessions():
    """Test protein accession extraction."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # Check protein accessions are lists
    assert df.iloc[0]["protein_accessions"] == ["PROT0"]
    assert df.iloc[1]["protein_accessions"] == ["PROT1"]


def test_psm_df_additional_scores():
    """Test additional scores extraction from metavalues (QPX array format).

    Only known score types (from Scores.getAllIDScoreNames()) go to additional_scores.
    Unknown numeric metavalues go to psm_metavalues instead.
    """
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.0)
    pep_id.setMZ(500.0)

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    # Known score type - should go to additional_scores
    hit.setMetaValue("q-value", 0.05)
    # Unknown score type - should go to psm_metavalues
    hit.setMetaValue("some_score", 42.5)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    # Check additional_scores has known score types
    scores = to_list(df.iloc[0]["additional_scores"])
    assert isinstance(scores, list)

    # Find q-value in the list (known score type)
    qvalue_entry = None
    for entry in scores:
        if entry["score_name"] == "q-value":
            qvalue_entry = entry
            break

    assert qvalue_entry is not None, "q-value should be in additional_scores"
    assert qvalue_entry["score_value"] == pytest.approx(0.05)
    # q-value is a known score type, so higher_better should be set (False for q-value)
    assert "higher_better" in qvalue_entry

    # Verify unknown score is NOT in additional_scores
    some_score_in_scores = any(e["score_name"] == "some_score" for e in scores)
    assert not some_score_in_scores, "Unknown score should NOT be in additional_scores"

    # Check psm_metavalues has unknown scores
    psm_mvs = to_list(df.iloc[0]["psm_metavalues"])
    assert isinstance(psm_mvs, list)

    # Find some_score in psm_metavalues
    some_score_entry = None
    for entry in psm_mvs:
        if entry["name"] == "some_score":
            some_score_entry = entry
            break

    assert some_score_entry is not None, "Unknown score should be in psm_metavalues"
    assert some_score_entry["value"] == "42.5"  # stringified in Arrow export


def test_psm_df_modifications():
    """Test modification extraction (QPX schema format)."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.0)
    pep_id.setMZ(600.0)
    pep_id.setScoreType("score")

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTM(Oxidation)IDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df(include_modifications=True)

    assert len(df) == 1
    assert "modifications" in df.columns

    mods = df.iloc[0]["modifications"]
    assert len(mods) == 1

    # QPX format: {"name", "accession", "positions": [{"position": "{AA}.{pos}", "scores"}]}
    mod = mods[0]
    assert mod["name"] == "Oxidation"
    assert "accession" in mod  # Should be UNIMOD:35 or None
    assert "positions" in mod
    assert len(mod["positions"]) == 1

    # Position format: "{AA}.{position}" where position is 1-indexed
    pos_entry = mod["positions"][0]
    assert pos_entry["position"] == "M.5"  # M is at position 5 (1-indexed)
    assert pos_entry["scores"] is None  # No localization scores


def test_psm_df_no_modifications():
    """Test with modifications disabled."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df(include_modifications=False)

    # When disabled, modifications column is still present but contains None
    assert "modifications" in df.columns
    assert df.iloc[0]["modifications"] is None


def test_psm_df_calculated_mz():
    """Test calculated m/z computation."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    # Check calculated_mz is present and reasonable
    assert "calculated_mz" in df.columns

    # For PEPTIDE with charge 2, calculated m/z should be around 400
    calc_mz = df.iloc[0]["calculated_mz"]
    assert calc_mz is not None
    assert 350 < calc_mz < 450


def test_psm_df_empty_identifications():
    """Test handling of empty PeptideIdentifications."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Add empty identification
    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.0)
    pep_id.setMZ(500.0)
    pep_id.setHits([])  # No hits
    pep_ids.append(pep_id)

    # Add identification with hit
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    pep_id2.setHits([hit])
    pep_ids.append(pep_id2)

    df = pep_ids.to_psm_df()

    # Should only have 1 row (empty identification skipped)
    assert len(df) == 1
    assert df.iloc[0]["P_ID"] == 1  # Second identification


def test_to_psm_qpx():
    """Test to_psm_qpx dict export with file_metadata and psms."""
    pep_ids = create_test_data()
    qpx_data = pep_ids.to_psm_qpx()

    # Check structure
    assert isinstance(qpx_data, dict)
    assert "file_metadata" in qpx_data
    assert "psms" in qpx_data

    # Check file_metadata
    metadata = qpx_data["file_metadata"]
    assert "qpx_version" in metadata
    assert "creator" in metadata
    assert "file_type" in metadata
    assert metadata["file_type"] == "psm"
    assert "creation_date" in metadata
    assert "uuid" in metadata
    assert "scan_format" in metadata
    assert "software_provider" in metadata

    # Check psms
    psms = qpx_data["psms"]
    assert len(psms) == 3
    # QPX schema fields
    assert "sequence" in psms[0]
    assert "peptidoform" in psms[0]
    assert "precursor_charge" in psms[0]
    assert "reference_file_name" in psms[0]
    assert "cv_params" in psms[0]


def test_to_psm_qpx_with_params():
    """Test to_psm_qpx with parameters."""
    pep_ids = create_test_data()
    qpx_data = pep_ids.to_psm_qpx(export_all_hits=False)

    assert len(qpx_data["psms"]) == 2


def test_to_psm_qpx_with_reference_file():
    """Test to_psm_qpx with reference_file_name parameter."""
    pep_ids = create_test_data()
    qpx_data = pep_ids.to_psm_qpx(reference_file_name="test.mzML")

    # All PSMs should have the reference_file_name
    for psm in qpx_data["psms"]:
        assert psm["reference_file_name"] == "test.mzML"


def test_to_psm_arrow():
    """Test to_psm_arrow Arrow Table export."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_psm_arrow()

    assert isinstance(table, pa.Table)
    assert table.num_rows == 3
    # Check expected columns
    assert "sequence" in table.schema.names
    assert "peptidoform" in table.schema.names
    assert "precursor_charge" in table.schema.names
    assert "reference_file_name" in table.schema.names
    assert "cv_params" in table.schema.names


def test_to_parquet(tmp_path):
    """Test to_parquet Parquet file export."""
    pytest.importorskip("pyarrow")
    import pandas as pd

    pep_ids = create_test_data()
    parquet_path = tmp_path / "psms.parquet"

    pep_ids.to_parquet(str(parquet_path))

    # Verify file exists and can be read back
    assert parquet_path.exists()
    df = pd.read_parquet(parquet_path)
    assert len(df) == 3
    assert "sequence" in df.columns
    assert "peptidoform" in df.columns
    assert "precursor_charge" in df.columns


def test_to_psm_qpx_file_metadata_params():
    """Test to_psm_qpx file_metadata customization."""
    pep_ids = create_test_data()
    qpx_data = pep_ids.to_psm_qpx(
        qpx_version="2.0",
        creator="TestCreator",
        software_provider="TestProvider",
        scan_format="nativeId"
    )

    metadata = qpx_data["file_metadata"]
    assert metadata["qpx_version"] == "2.0"
    assert metadata["creator"] == "TestCreator"
    assert metadata["software_provider"] == "TestProvider"
    assert metadata["scan_format"] == "nativeId"


def test_psm_columns():
    """Test psm_columns instance method."""
    pep_ids = create_test_data()
    columns = pep_ids.psm_columns()

    assert isinstance(columns, list)
    # QPX schema fields (in schema order)
    assert "sequence" in columns
    assert "peptidoform" in columns
    assert "modifications" in columns
    assert "precursor_charge" in columns
    assert "posterior_error_probability" in columns
    assert "is_decoy" in columns
    assert "calculated_mz" in columns
    assert "observed_mz" in columns
    assert "additional_scores" in columns
    assert "protein_accessions" in columns
    assert "predicted_rt" in columns
    assert "reference_file_name" in columns  # QPX required
    assert "cv_params" in columns  # QPX optional
    assert "scan" in columns
    assert "rt" in columns
    assert "ion_mobility" in columns
    assert "number_peaks" in columns
    assert "mz_array" in columns
    assert "intensity_array" in columns
    assert "charge_array" in columns
    assert "ion_type_array" in columns  # QPX naming
    assert "ion_mobility_array" in columns  # QPX field
    # OpenMS-specific fields
    assert "P_ID" in columns
    assert "hit_index" in columns
    assert "spectrum_reference" in columns
    assert "score" in columns
    assert "score_type" in columns


def test_psm_df_qpx_columns():
    """Test QPX schema columns are present."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df(reference_file_name="test.mzML")

    # Check QPX required columns
    assert "reference_file_name" in df.columns
    assert df.iloc[0]["reference_file_name"] == "test.mzML"

    # Check cv_params column (always None from OpenMS)
    assert "cv_params" in df.columns
    assert pd.isna(df.iloc[0]["cv_params"])


def test_psm_df_spectrum_reference():
    """Test spectrum reference preservation."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    assert "spectrum_reference" in df.columns
    assert df.iloc[0]["spectrum_reference"] == "controllerType=0 controllerNumber=1 scan=1234"
    assert df.iloc[2]["spectrum_reference"] == "scan=5678"


def test_psm_df_sequence_formats():
    """Test that both unmodified and modified sequence formats are present."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTM(Oxidation)IDE"))
    hit.setCharge(2)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    # sequence should be unmodified
    assert df.iloc[0]["sequence"] == "PEPTMIDE"

    # peptidoform should contain modification (ProForma uses UNIMOD accession)
    assert "UNIMOD:35" in df.iloc[0]["peptidoform"] or "Oxidation" in df.iloc[0]["peptidoform"]


def test_psm_df_terminal_modifications():
    """Test N-terminal and C-terminal modification extraction."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    # Create sequence with N-terminal acetylation
    hit.setSequence(oms.AASequence.fromString(".(Acetyl)PEPTIDE"))
    hit.setCharge(2)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df(include_modifications=True)

    mods = df.iloc[0]["modifications"]
    assert len(mods) >= 1
    # N-terminal mod should have position 0
    n_term_mods = [m for m in mods if m["positions"][0]["position"].endswith(".0")]
    assert len(n_term_mods) == 1
    assert "Acetyl" in n_term_mods[0]["name"]


def test_psm_df_ion_mobility():
    """Test ion mobility extraction from metavalues."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setMetaValue("ion_mobility", 0.95)

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert "ion_mobility" in df.columns
    assert df.iloc[0]["ion_mobility"] == pytest.approx(0.95)


def test_psm_df_ion_mobility_im_key():
    """Test ion mobility extraction using IM metavalue key."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setMetaValue("IM", 1.25)  # Alternative key

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert df.iloc[0]["ion_mobility"] == pytest.approx(1.25)


def test_psm_df_predicted_rt():
    """Test predicted_rt extraction from metavalues."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setMetaValue("predicted_RT", 150.5)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert "predicted_rt" in df.columns
    assert df.iloc[0]["predicted_rt"] == pytest.approx(150.5)


def test_psm_df_missing_target_decoy():
    """Test is_decoy is None when target_decoy metavalue is missing."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    # Note: NOT setting target_decoy metavalue
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert is_null(df.iloc[0]["is_decoy"])


def test_psm_df_charge_zero():
    """Test calculated_mz is None when charge is 0."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(0)  # Zero charge
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert is_null(df.iloc[0]["calculated_mz"])
    assert df.iloc[0]["precursor_charge"] == 0


def test_psm_df_empty_list():
    """Test handling of completely empty PeptideIdentificationList."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    # No identifications added

    df = pep_ids.to_psm_df()

    assert len(df) == 0
    # Columns should still be defined even with no data
    assert "sequence" in df.columns
    assert "peptidoform" in df.columns
    assert "precursor_charge" in df.columns
    assert "P_ID" in df.columns


def test_psm_df_non_pep_score_type():
    """Test posterior_error_probability is None when score type is not PEP."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setScoreType("Hyperscore")  # Not a PEP score type

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(100.0)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert is_null(df.iloc[0]["posterior_error_probability"])
    assert df.iloc[0]["score"] == pytest.approx(100.0)
    assert df.iloc[0]["score_type"] == "Hyperscore"


def test_psm_df_multiple_protein_accessions():
    """Test multiple protein accessions for a single hit."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)

    ev1 = oms.PeptideEvidence()
    ev1.setProteinAccession("PROT1")
    ev2 = oms.PeptideEvidence()
    ev2.setProteinAccession("PROT2")
    ev3 = oms.PeptideEvidence()
    ev3.setProteinAccession("PROT3")
    hit.setPeptideEvidences([ev1, ev2, ev3])
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert to_list(df.iloc[0]["protein_accessions"]) == ["PROT1", "PROT2", "PROT3"]


def test_psm_df_peak_annotations():
    """Test include_peak_annotations parameter."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)

    # Add peak annotations
    pa1 = oms.PeptideHit_PeakAnnotation()
    pa1.mz = 100.5
    pa1.intensity = 1000.0
    pa1.charge = 1
    pa1.annotation = "y1"

    pa2 = oms.PeptideHit_PeakAnnotation()
    pa2.mz = 200.5
    pa2.intensity = 2000.0
    pa2.charge = 1
    pa2.annotation = "y2"

    hit.setPeakAnnotations([pa1, pa2])
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df(include_peak_annotations=True)

    assert "number_peaks" in df.columns
    assert "mz_array" in df.columns
    assert "intensity_array" in df.columns
    assert "charge_array" in df.columns
    assert "ion_type_array" in df.columns

    assert df.iloc[0]["number_peaks"] == 2
    assert to_list(df.iloc[0]["mz_array"]) == [100.5, 200.5]
    assert to_list(df.iloc[0]["intensity_array"]) == [1000.0, 2000.0]
    assert to_list(df.iloc[0]["charge_array"]) == [1, 1]
    assert to_list(df.iloc[0]["ion_type_array"]) == ["y1", "y2"]


def test_psm_df_peak_annotations_empty():
    """Test include_peak_annotations with no annotations."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    # No peak annotations set
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df(include_peak_annotations=True)

    assert df.iloc[0]["number_peaks"] == 0
    assert len(to_list(df.iloc[0]["mz_array"])) == 0
    assert len(to_list(df.iloc[0]["intensity_array"])) == 0


def test_psm_df_peak_annotations_schema_consistency():
    """Test that peak annotation columns are included even when include_peak_annotations=False.

    This ensures schema consistency between psm_columns() and actual DataFrame columns.
    """
    pep_ids = create_test_data()

    # Get DataFrame without peak annotations (default)
    df = pep_ids.to_psm_df(include_peak_annotations=False)

    # Peak annotation columns should still be present for schema consistency
    peak_cols = ["number_peaks", "mz_array", "intensity_array",
                 "charge_array", "ion_type_array", "ion_mobility_array"]
    for col in peak_cols:
        assert col in df.columns, f"Column {col} should be present even with include_peak_annotations=False"
        # Values should be None/NaN when not requested
        assert is_null(df.iloc[0][col]), f"Column {col} should be null when include_peak_annotations=False"

    # Verify psm_columns() includes these columns
    expected_cols = pep_ids.psm_columns()
    for col in peak_cols:
        assert col in expected_cols, f"Column {col} should be in psm_columns()"


def test_psm_df_sorting():
    """Test that PSM DataFrame is sorted by rt, observed_mz, precursor_charge, rank."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Add identifications in reverse RT order to verify sorting
    # Third ID: RT=300, MZ=700
    pep_id3 = oms.PeptideIdentification()
    pep_id3.setRT(300.0)
    pep_id3.setMZ(700.0)
    hit3 = oms.PeptideHit()
    hit3.setSequence(oms.AASequence.fromString("THIRD"))
    hit3.setCharge(2)
    hit3.setScore(0.1)
    pep_id3.setHits([hit3])
    pep_ids.append(pep_id3)

    # First ID: RT=100, MZ=500
    pep_id1 = oms.PeptideIdentification()
    pep_id1.setRT(100.0)
    pep_id1.setMZ(500.0)
    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("FIRST"))
    hit1.setCharge(2)
    hit1.setScore(0.1)
    pep_id1.setHits([hit1])
    pep_ids.append(pep_id1)

    # Second ID: RT=200, MZ=600, with two hits (different ranks)
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    hit2a = oms.PeptideHit()
    hit2a.setSequence(oms.AASequence.fromString("SECNDA"))
    hit2a.setCharge(2)
    hit2a.setScore(0.1)  # rank 0
    hit2b = oms.PeptideHit()
    hit2b.setSequence(oms.AASequence.fromString("SECNDC"))
    hit2b.setCharge(2)
    hit2b.setScore(0.2)  # rank 1
    pep_id2.setHits([hit2a, hit2b])
    pep_ids.append(pep_id2)

    df = pep_ids.to_psm_df()

    # Should be sorted by RT first: 100, 200, 200, 300
    assert list(df["rt"]) == [100.0, 200.0, 200.0, 300.0]

    # For same RT (200), should be sorted by rank: 0, 1
    assert list(df["sequence"]) == ["FIRST", "SECNDA", "SECNDC", "THIRD"]
    assert list(df["hit_index"]) == [0, 0, 1, 0]


def test_psm_df_sorting_by_mz():
    """Test sorting by observed_mz when rt is the same."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Same RT but different MZ
    pep_id1 = oms.PeptideIdentification()
    pep_id1.setRT(100.0)
    pep_id1.setMZ(700.0)  # Higher MZ
    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("GAMMAAA"))
    hit1.setCharge(2)
    pep_id1.setHits([hit1])
    pep_ids.append(pep_id1)

    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(100.0)
    pep_id2.setMZ(500.0)  # Lower MZ
    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("ALPHAAA"))
    hit2.setCharge(2)
    pep_id2.setHits([hit2])
    pep_ids.append(pep_id2)

    df = pep_ids.to_psm_df()

    # Should be sorted by MZ when RT is same: 500, 700
    assert list(df["observed_mz"]) == [500.0, 700.0]
    assert list(df["sequence"]) == ["ALPHAAA", "GAMMAAA"]


def test_psm_df_sorting_by_charge():
    """Test sorting by precursor_charge when rt and mz are the same."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Same RT and MZ but different charge
    pep_id1 = oms.PeptideIdentification()
    pep_id1.setRT(100.0)
    pep_id1.setMZ(500.0)
    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("CHARGETHREE"))
    hit1.setCharge(3)
    pep_id1.setHits([hit1])
    pep_ids.append(pep_id1)

    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(100.0)
    pep_id2.setMZ(500.0)
    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("CHARGEAA"))
    hit2.setCharge(2)
    pep_id2.setHits([hit2])
    pep_ids.append(pep_id2)

    df = pep_ids.to_psm_df()

    # Should be sorted by charge when RT and MZ are same: 2, 3
    assert list(df["precursor_charge"]) == [2, 3]
    assert list(df["sequence"]) == ["CHARGEAA", "CHARGETHREE"]


# === Native Arrow Export Tests ===

def test_to_arrow_basic():
    """Test basic to_arrow() export (simpler format matching to_df)."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_arrow()

    assert isinstance(table, pa.Table)
    # to_arrow only exports first hit per identification
    assert table.num_rows == 2

    # Check expected columns
    assert "id" in table.schema.names
    assert "rt" in table.schema.names
    assert "mz" in table.schema.names
    assert "charge" in table.schema.names
    assert "protein_accession" in table.schema.names
    assert "P_ID" in table.schema.names
    assert "PSM_ID" in table.schema.names


def test_to_arrow_column_types():
    """Test Arrow column types are correct."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_arrow()

    schema = table.schema
    assert schema.field("rt").type == pa.float32()
    assert schema.field("mz").type == pa.float32()
    assert schema.field("charge").type == pa.int32()
    assert schema.field("P_ID").type == pa.int32()
    assert schema.field("id").type in (pa.utf8(), pa.large_utf8())


def test_to_arrow_column_filter():
    """Test to_arrow column filtering."""
    pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_arrow(columns=['id', 'rt', 'mz'])

    assert table.num_columns == 3
    assert set(table.schema.names) == {'id', 'rt', 'mz'}


def test_to_arrow_empty_list():
    """Test to_arrow with empty PeptideIdentificationList preserves schema."""
    pa = pytest.importorskip("pyarrow")
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    table = pep_ids.to_arrow()

    assert isinstance(table, pa.Table)
    assert table.num_rows == 0

    # Empty tables should still have full schema for consistency
    # Core columns must be present even with no data
    expected_columns = ["id", "rt", "mz", "charge", "protein_accession",
                        "start", "end", "P_ID", "PSM_ID"]
    for col in expected_columns:
        assert col in table.schema.names, f"Missing column {col} in empty table schema"

    # Verify types are correct
    assert table.schema.field("id").type == pa.utf8()
    assert table.schema.field("rt").type == pa.float32()
    assert table.schema.field("charge").type == pa.int32()


def test_to_psm_arrow_schema_types():
    """Test to_psm_arrow nested Arrow schema types."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_psm_arrow()

    schema = table.schema

    # Check scalar types
    assert schema.field("sequence").type == pa.utf8()
    assert schema.field("precursor_charge").type == pa.int32()
    assert schema.field("observed_mz").type == pa.float64()
    assert schema.field("rt").type == pa.float64()

    # Check list types
    assert pa.types.is_list(schema.field("protein_accessions").type)

    # Check nested struct types
    mods_type = schema.field("modifications").type
    assert pa.types.is_list(mods_type)
    # The inner type should be a struct
    inner_type = mods_type.value_type
    assert pa.types.is_struct(inner_type)


def test_to_psm_arrow_column_filter():
    """Test to_psm_arrow column filtering."""
    pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_psm_arrow(columns=['sequence', 'precursor_charge', 'score'])

    assert table.num_columns == 3
    assert set(table.schema.names) == {'sequence', 'precursor_charge', 'score'}


def test_to_psm_arrow_empty_list():
    """Test to_psm_arrow with empty PeptideIdentificationList preserves schema."""
    pa = pytest.importorskip("pyarrow")
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    table = pep_ids.to_psm_arrow()

    assert isinstance(table, pa.Table)
    assert table.num_rows == 0

    # Empty tables should have full schema for Parquet compatibility
    # and consumer code that relies on consistent column structure
    expected_columns = ["sequence", "peptidoform", "modifications",
                        "precursor_charge", "posterior_error_probability",
                        "is_decoy", "observed_mz", "protein_accessions",
                        "rt", "score", "hit_index", "P_ID"]
    for col in expected_columns:
        assert col in table.schema.names, f"Missing column {col} in empty table schema"

    # Verify nested types are correct even when empty
    mods_type = table.schema.field("modifications").type
    assert pa.types.is_list(mods_type)
    assert pa.types.is_struct(mods_type.value_type)

    scores_type = table.schema.field("additional_scores").type
    assert pa.types.is_list(scores_type)
    assert pa.types.is_struct(scores_type.value_type)


def test_to_psm_arrow_vs_psm_df_data_equivalence():
    """Test that to_psm_arrow produces equivalent data to to_psm_df."""
    pytest.importorskip("pyarrow")

    pep_ids = create_test_data()

    # Get data via both methods
    df = pep_ids.to_psm_df()
    table = pep_ids.to_psm_arrow()
    arrow_df = table.to_pandas()

    # Check row count matches
    assert len(df) == len(arrow_df)

    # Check scalar column values match
    for col in ['sequence', 'precursor_charge', 'observed_mz', 'rt', 'hit_index']:
        if col in df.columns and col in arrow_df.columns:
            assert list(df[col]) == list(arrow_df[col]), f"Mismatch in column {col}"


def test_to_psm_arrow_sorting():
    """Test that to_psm_arrow sorts by rt, observed_mz, precursor_charge, rank."""
    pytest.importorskip("pyarrow")
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()

    # Add identifications in reverse RT order
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("SECNDA"))
    hit2.setCharge(2)
    pep_id2.setHits([hit2])
    pep_ids.append(pep_id2)

    pep_id1 = oms.PeptideIdentification()
    pep_id1.setRT(100.0)
    pep_id1.setMZ(500.0)
    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("FIRST"))
    hit1.setCharge(2)
    pep_id1.setHits([hit1])
    pep_ids.append(pep_id1)

    table = pep_ids.to_psm_arrow()
    df = table.to_pandas()

    # Should be sorted by RT: 100, 200
    assert list(df["rt"]) == [100.0, 200.0]
    assert list(df["sequence"]) == ["FIRST", "SECNDA"]


def test_to_parquet_compression_options(tmp_path):
    """Test to_parquet with different compression options."""
    pytest.importorskip("pyarrow")
    import os

    pep_ids = create_test_data()

    # Test different compression algorithms
    for compression in ['none', 'snappy', 'gzip', 'zstd']:
        parquet_path = tmp_path / f"psms_{compression}.parquet"
        pep_ids.to_parquet(str(parquet_path), compression=compression)
        assert parquet_path.exists()
        # File should be readable
        import pandas as pd
        df = pd.read_parquet(parquet_path)
        assert len(df) == 3


def test_to_parquet_invalid_compression(tmp_path):
    """Test to_parquet with invalid compression raises error."""
    pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    parquet_path = tmp_path / "psms.parquet"

    with pytest.raises(ValueError, match="Unsupported compression"):
        pep_ids.to_parquet(str(parquet_path), compression='invalid')


def test_to_parquet_row_group_size(tmp_path):
    """Test to_parquet with custom row_group_size."""
    pq = pytest.importorskip("pyarrow.parquet")

    pep_ids = create_test_data()
    parquet_path = tmp_path / "psms.parquet"

    # Set row group size to 1 to get multiple row groups
    pep_ids.to_parquet(str(parquet_path), row_group_size=1)

    # Read metadata to verify row groups
    metadata = pq.read_metadata(parquet_path)
    assert metadata.num_row_groups == 3  # One per row


def test_to_parquet_with_peak_annotations(tmp_path):
    """Test to_parquet with peak annotations included."""
    pytest.importorskip("pyarrow")
    import pandas as pd
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)

    pa1 = oms.PeptideHit_PeakAnnotation()
    pa1.mz = 100.5
    pa1.intensity = 1000.0
    pa1.charge = 1
    pa1.annotation = "y1"
    hit.setPeakAnnotations([pa1])

    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    parquet_path = tmp_path / "psms_with_peaks.parquet"
    pep_ids.to_parquet(str(parquet_path), include_peak_annotations=True)

    df = pd.read_parquet(parquet_path)
    assert "mz_array" in df.columns
    assert "intensity_array" in df.columns
    assert df.iloc[0]["number_peaks"] == 1


def test_to_parquet_reference_file_name(tmp_path):
    """Test to_parquet with reference_file_name."""
    pytest.importorskip("pyarrow")
    import pandas as pd

    pep_ids = create_test_data()
    parquet_path = tmp_path / "psms.parquet"

    pep_ids.to_parquet(str(parquet_path), reference_file_name="test.mzML")

    df = pd.read_parquet(parquet_path)
    assert all(df["reference_file_name"] == "test.mzML")


# === Metavalue Export Tests ===

def test_psm_df_metavalue_columns_present():
    """Test that psm_metavalues and spectrum_metavalues columns are present."""
    pep_ids = create_test_data()
    df = pep_ids.to_psm_df()

    assert "psm_metavalues" in df.columns
    assert "spectrum_metavalues" in df.columns


def test_psm_df_spectrum_metavalues():
    """Test spectrum_metavalues (PeptideIdentification-level metavalues)."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.0)
    pep_id.setMZ(500.0)
    # Set some spectrum-level metavalues
    pep_id.setMetaValue("custom_spectrum_value", 42.5)
    pep_id.setMetaValue("spectrum_label", "test_spectrum")
    # spectrum_reference is excluded (has dedicated column)
    pep_id.setMetaValue("spectrum_reference", "scan=1234")

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    spectrum_mvs = to_list(df.iloc[0]["spectrum_metavalues"])
    assert isinstance(spectrum_mvs, list)

    # Should have custom_spectrum_value and spectrum_label
    mv_names = [mv["name"] for mv in spectrum_mvs]
    assert "custom_spectrum_value" in mv_names
    assert "spectrum_label" in mv_names
    # spectrum_reference should be excluded (has dedicated column)
    assert "spectrum_reference" not in mv_names

    # Check values and value_type (values are stringified in Arrow export)
    for mv in spectrum_mvs:
        if mv["name"] == "custom_spectrum_value":
            assert mv["value"] == "42.5"  # stringified
            assert mv["value_type"] == "float"
        elif mv["name"] == "spectrum_label":
            assert mv["value"] == "test_spectrum"
            assert mv["value_type"] == "str"


def test_psm_df_psm_metavalues():
    """Test psm_metavalues (PeptideHit-level non-score metavalues)."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    # Set various metavalues - some are known scores, some are not
    hit.setMetaValue("custom_property", "my_value")
    hit.setMetaValue("numeric_property", 123)
    hit.setMetaValue("target_decoy", "target")  # excluded - has dedicated column
    hit.setMetaValue("predicted_RT", 1000.0)  # excluded - has dedicated column

    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    psm_mvs = to_list(df.iloc[0]["psm_metavalues"])
    assert isinstance(psm_mvs, list)

    # Should have custom_property and numeric_property
    mv_names = [mv["name"] for mv in psm_mvs]
    assert "custom_property" in mv_names
    assert "numeric_property" in mv_names
    # Excluded metavalues should not be present
    assert "target_decoy" not in mv_names
    assert "predicted_RT" not in mv_names

    # Check values and value_type (values are stringified in Arrow export)
    for mv in psm_mvs:
        if mv["name"] == "custom_property":
            assert mv["value"] == "my_value"
            assert mv["value_type"] == "str"
        elif mv["name"] == "numeric_property":
            assert mv["value"] == "123"  # stringified
            assert mv["value_type"] == "int"


def test_psm_df_score_vs_metavalue_distinction():
    """Test that known scores go to additional_scores, others to psm_metavalues."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)

    # Known score types (should go to additional_scores)
    hit.setMetaValue("q-value", 0.05)
    hit.setMetaValue("E-Value", 0.001)

    # Non-score metavalues (should go to psm_metavalues)
    hit.setMetaValue("custom_metric", 42.0)
    hit.setMetaValue("string_property", "test")

    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    additional_scores = to_list(df.iloc[0]["additional_scores"])
    psm_metavalues = to_list(df.iloc[0]["psm_metavalues"])

    # Known scores should be in additional_scores
    score_names = [s["score_name"] for s in additional_scores]
    assert "q-value" in score_names
    assert "E-Value" in score_names

    # Non-scores should be in psm_metavalues
    mv_names = [mv["name"] for mv in psm_metavalues]
    assert "custom_metric" in mv_names
    assert "string_property" in mv_names

    # Cross-check: known scores should NOT be in psm_metavalues
    assert "q-value" not in mv_names
    assert "E-Value" not in mv_names


def test_psm_df_metavalues_empty_when_none():
    """Test metavalue columns are empty lists when no metavalues exist."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    # No metavalues set
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    df = pep_ids.to_psm_df()

    assert len(to_list(df.iloc[0]["psm_metavalues"])) == 0
    assert len(to_list(df.iloc[0]["spectrum_metavalues"])) == 0


def test_psm_columns_includes_metavalue_columns():
    """Test psm_columns() includes the new metavalue columns."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    columns = pep_ids.psm_columns()

    assert "psm_metavalues" in columns
    assert "spectrum_metavalues" in columns


def test_to_psm_arrow_metavalue_columns():
    """Test to_psm_arrow includes metavalue columns with correct types."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_psm_arrow()

    assert "psm_metavalues" in table.schema.names
    assert "spectrum_metavalues" in table.schema.names

    # Check types are list of structs
    psm_mv_type = table.schema.field("psm_metavalues").type
    assert pa.types.is_list(psm_mv_type)
    assert pa.types.is_struct(psm_mv_type.value_type)

    # Check struct has name, value, and value_type fields
    struct_type = psm_mv_type.value_type
    field_names = [struct_type.field(i).name for i in range(struct_type.num_fields)]
    assert "name" in field_names
    assert "value" in field_names
    assert "value_type" in field_names

    spectrum_mv_type = table.schema.field("spectrum_metavalues").type
    assert pa.types.is_list(spectrum_mv_type)
    assert pa.types.is_struct(spectrum_mv_type.value_type)


def test_to_psm_arrow_metavalue_data():
    """Test to_psm_arrow correctly exports metavalue data including value_type."""
    pytest.importorskip("pyarrow")
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setMetaValue("spectrum_string", "test_value")
    pep_id.setMetaValue("spectrum_float", 3.14)
    pep_id.setMetaValue("spectrum_int", 42)

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setMetaValue("hit_string", "hit_value")
    hit.setMetaValue("hit_float", 2.71)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    table = pep_ids.to_psm_arrow()
    df = table.to_pandas()

    # Check spectrum_metavalues has correct values and types
    spectrum_mvs = df.iloc[0]["spectrum_metavalues"]
    spectrum_mv_dict = {mv["name"]: mv for mv in spectrum_mvs}

    assert spectrum_mv_dict["spectrum_string"]["value"] == "test_value"
    assert spectrum_mv_dict["spectrum_string"]["value_type"] == "str"

    assert spectrum_mv_dict["spectrum_float"]["value"] == "3.14"
    assert spectrum_mv_dict["spectrum_float"]["value_type"] == "float"

    assert spectrum_mv_dict["spectrum_int"]["value"] == "42"
    assert spectrum_mv_dict["spectrum_int"]["value_type"] == "int"

    # Check psm_metavalues has correct values and types
    psm_mvs = df.iloc[0]["psm_metavalues"]
    psm_mv_dict = {mv["name"]: mv for mv in psm_mvs}

    assert psm_mv_dict["hit_string"]["value"] == "hit_value"
    assert psm_mv_dict["hit_string"]["value_type"] == "str"

    assert psm_mv_dict["hit_float"]["value"] == "2.71"
    assert psm_mv_dict["hit_float"]["value_type"] == "float"


def test_to_parquet_metavalue_columns(tmp_path):
    """Test to_parquet correctly exports metavalue columns."""
    pytest.importorskip("pyarrow")
    import pandas as pd
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()
    pep_id.setMetaValue("spectrum_tag", "sample_1")

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setMetaValue("hit_tag", "annotated")
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    parquet_path = tmp_path / "psms_with_metavalues.parquet"
    pep_ids.to_parquet(str(parquet_path))

    df = pd.read_parquet(parquet_path)
    assert "psm_metavalues" in df.columns
    assert "spectrum_metavalues" in df.columns

    # Verify data survived the roundtrip
    assert len(df.iloc[0]["spectrum_metavalues"]) == 1
    assert len(df.iloc[0]["psm_metavalues"]) == 1


def test_psm_df_additional_score_names_parameter():
    """Test additional_score_names parameter moves custom scores to additional_scores."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    # Custom score - not in known score names
    hit.setMetaValue("my_custom_score", 0.95)
    hit.setMetaValue("another_metric", 42.0)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    # Without additional_score_names: custom score goes to psm_metavalues
    df_without = pep_ids.to_psm_df()
    scores_without = df_without.iloc[0]["additional_scores"]
    psm_mvs_without = df_without.iloc[0]["psm_metavalues"]

    assert not any(s["score_name"] == "my_custom_score" for s in scores_without)
    assert any(mv["name"] == "my_custom_score" for mv in psm_mvs_without)

    # With additional_score_names: custom score goes to additional_scores
    df_with = pep_ids.to_psm_df(additional_score_names=["my_custom_score"])
    scores_with = df_with.iloc[0]["additional_scores"]
    psm_mvs_with = df_with.iloc[0]["psm_metavalues"]

    # my_custom_score should now be in additional_scores
    custom_score_entry = None
    for entry in scores_with:
        if entry["score_name"] == "my_custom_score":
            custom_score_entry = entry
            break

    assert custom_score_entry is not None, "Custom score should be in additional_scores"
    assert custom_score_entry["score_value"] == pytest.approx(0.95)

    # my_custom_score should NOT be in psm_metavalues
    assert not any(mv["name"] == "my_custom_score" for mv in psm_mvs_with)

    # another_metric should still be in psm_metavalues (not specified as score)
    assert any(mv["name"] == "another_metric" for mv in psm_mvs_with)


def test_to_psm_arrow_additional_score_names_parameter():
    """Test additional_score_names parameter works with to_psm_arrow."""
    pytest.importorskip("pyarrow")
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    pep_id = oms.PeptideIdentification()

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setMetaValue("custom_rescore", 0.99)
    pep_id.setHits([hit])
    pep_ids.append(pep_id)

    # With additional_score_names
    table = pep_ids.to_psm_arrow(additional_score_names=["custom_rescore"])
    df = table.to_pandas()

    scores = df.iloc[0]["additional_scores"]
    assert any(s["score_name"] == "custom_rescore" for s in scores)


def test_scan_format_native_id():
    """Test scan_format='nativeId' puts raw native ID in scan column."""
    pep_ids = create_test_data()

    # Default scan_format="scan" extracts scan numbers
    df_scan = pep_ids.to_psm_df()
    assert df_scan.iloc[0]["scan"] == "1234"

    # scan_format="nativeId" uses raw native ID string
    df_native = pep_ids.to_psm_df(scan_format="nativeId")
    assert df_native.iloc[0]["scan"] == "controllerType=0 controllerNumber=1 scan=1234"

    # Also works via to_psm_arrow
    table = pep_ids.to_psm_arrow(scan_format="nativeId")
    assert table.column("scan")[0].as_py() == "controllerType=0 controllerNumber=1 scan=1234"


def test_scan_format_invalid():
    """Test that invalid scan_format raises ValueError."""
    pep_ids = create_test_data()

    with pytest.raises(ValueError, match="scan_format must be"):
        pep_ids.to_psm_arrow(scan_format="invalid")

    with pytest.raises(ValueError, match="scan_format must be"):
        pep_ids.to_psm_df(scan_format="index")
