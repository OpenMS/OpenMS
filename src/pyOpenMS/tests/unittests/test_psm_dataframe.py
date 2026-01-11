"""
Unit tests for PSM DataFrame export from PeptideIdentificationList.

Tests the get_psm_df(), to_qpx(), and get_psm_columns() methods
that follow the QPX PSM schema.
"""

import pytest


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
    df = pep_ids.get_psm_df()

    # Should have 3 rows (2 hits from first ID + 1 from second)
    assert len(df) == 3

    # Check required columns exist
    required_cols = [
        "sequence", "peptidoform", "precursor_charge",
        "observed_mz", "rt", "rank", "score", "score_type",
        "protein_accessions", "P_ID"
    ]
    for col in required_cols:
        assert col in df.columns, f"Missing column: {col}"

    # Check ranks are correct
    assert list(df["rank"]) == [0, 1, 0]

    # Check P_ID correctly tracks identification index
    assert list(df["P_ID"]) == [0, 0, 1]


def test_psm_df_top_hit_only():
    """Test exporting only top hit."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df(export_all_hits=False)

    # Should have 2 rows (1 per identification)
    assert len(df) == 2
    assert list(df["rank"]) == [0, 0]
    assert list(df["P_ID"]) == [0, 1]


def test_psm_df_scan_parsing():
    """Test scan number extraction from spectrum reference."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

    # First two PSMs from first identification with scan=1234 (QPX uses string)
    assert df.iloc[0]["scan"] == "1234"
    assert df.iloc[1]["scan"] == "1234"

    # Third PSM from second identification with scan=5678
    assert df.iloc[2]["scan"] == "5678"


def test_psm_df_decoy_detection():
    """Test is_decoy field extraction (QPX uses int: 0=target, 1=decoy)."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

    # First two are targets (0)
    assert df.iloc[0]["is_decoy"] == 0
    assert df.iloc[1]["is_decoy"] == 0

    # Third is decoy (1)
    assert df.iloc[2]["is_decoy"] == 1


def test_psm_df_pep_score():
    """Test posterior_error_probability extraction."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

    # All should have PEP since score types contain "PEP"
    assert "posterior_error_probability" in df.columns

    # Check values match scores
    assert df.iloc[0]["posterior_error_probability"] == pytest.approx(0.01)
    assert df.iloc[2]["posterior_error_probability"] == pytest.approx(0.001)


def test_psm_df_protein_accessions():
    """Test protein accession extraction."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

    # Check protein accessions are lists
    assert df.iloc[0]["protein_accessions"] == ["PROT0"]
    assert df.iloc[1]["protein_accessions"] == ["PROT1"]


def test_psm_df_additional_scores():
    """Test additional scores extraction from metavalues."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

    # First hit has some_score metavalue
    scores = df.iloc[0]["additional_scores"]
    assert "some_score" in scores
    assert scores["some_score"] == pytest.approx(42.5)


def test_psm_df_modifications():
    """Test modification extraction."""
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

    df = pep_ids.get_psm_df(include_modifications=True)

    assert len(df) == 1
    assert "modifications" in df.columns

    mods = df.iloc[0]["modifications"]
    assert len(mods) == 1
    assert mods[0]["name"] == "Oxidation"
    assert mods[0]["position"] == 5  # M is at position 5 (1-indexed)
    assert mods[0]["mass"] == pytest.approx(15.9949, abs=0.01)


def test_psm_df_no_modifications():
    """Test with modifications disabled."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df(include_modifications=False)

    # When disabled, modifications column is still present but contains None
    assert "modifications" in df.columns
    assert df.iloc[0]["modifications"] is None


def test_psm_df_calculated_mz():
    """Test calculated m/z computation."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

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

    df = pep_ids.get_psm_df()

    # Should only have 1 row (empty identification skipped)
    assert len(df) == 1
    assert df.iloc[0]["P_ID"] == 1  # Second identification


def test_to_qpx():
    """Test to_qpx Arrow export."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_qpx()

    assert isinstance(table, pa.Table)
    assert table.num_rows == 3
    # QPX schema fields
    assert "sequence" in table.schema.names
    assert "peptidoform" in table.schema.names
    assert "precursor_charge" in table.schema.names
    assert "reference_file_name" in table.schema.names
    assert "cv_params" in table.schema.names


def test_to_qpx_with_params():
    """Test to_qpx Arrow export with parameters."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_qpx(export_all_hits=False)

    assert table.num_rows == 2


def test_to_qpx_with_reference_file():
    """Test to_qpx with reference_file_name parameter."""
    pa = pytest.importorskip("pyarrow")

    pep_ids = create_test_data()
    table = pep_ids.to_qpx(reference_file_name="test.mzML")

    # All rows should have the reference_file_name
    df = table.to_pandas()
    assert all(df["reference_file_name"] == "test.mzML")


def test_get_psm_columns():
    """Test get_psm_columns static method."""
    import pyopenms as oms

    columns = oms.PeptideIdentificationList.get_psm_columns()

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
    assert "rank" in columns
    assert "spectrum_reference" in columns
    assert "score" in columns
    assert "score_type" in columns


def test_psm_df_qpx_columns():
    """Test QPX schema columns are present."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df(reference_file_name="test.mzML")

    # Check QPX required columns
    assert "reference_file_name" in df.columns
    assert df.iloc[0]["reference_file_name"] == "test.mzML"

    # Check cv_params column (always None from OpenMS)
    assert "cv_params" in df.columns
    assert df.iloc[0]["cv_params"] is None


def test_psm_df_spectrum_reference():
    """Test spectrum reference preservation."""
    pep_ids = create_test_data()
    df = pep_ids.get_psm_df()

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

    df = pep_ids.get_psm_df()

    # sequence should be unmodified
    assert df.iloc[0]["sequence"] == "PEPTMIDE"

    # peptidoform should contain modification
    assert "Oxidation" in df.iloc[0]["peptidoform"]
