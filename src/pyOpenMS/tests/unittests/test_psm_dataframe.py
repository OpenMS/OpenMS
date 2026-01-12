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

        df = pep_ids.get_psm_df()
        assert df.iloc[0]["scan"] == expected_scan, f"Failed for '{spec_ref}': expected {expected_scan}, got {df.iloc[0]['scan']}"


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
    pytest.importorskip("pyarrow")

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

    df = pep_ids.get_psm_df(include_modifications=True)

    mods = df.iloc[0]["modifications"]
    assert len(mods) >= 1
    # N-terminal mod should have position 0
    n_term_mods = [m for m in mods if m["position"] == 0]
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

    df = pep_ids.get_psm_df()

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

    df = pep_ids.get_psm_df()

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

    df = pep_ids.get_psm_df()

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

    df = pep_ids.get_psm_df()

    assert df.iloc[0]["is_decoy"] is None


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

    df = pep_ids.get_psm_df()

    assert df.iloc[0]["calculated_mz"] is None
    assert df.iloc[0]["precursor_charge"] == 0


def test_psm_df_empty_list():
    """Test handling of completely empty PeptideIdentificationList."""
    import pyopenms as oms

    pep_ids = oms.PeptideIdentificationList()
    # No identifications added

    df = pep_ids.get_psm_df()

    assert len(df) == 0
    # Columns should still be defined even with no data
    assert "sequence" in df.columns or len(df.columns) == 0


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

    df = pep_ids.get_psm_df()

    assert df.iloc[0]["posterior_error_probability"] is None
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

    df = pep_ids.get_psm_df()

    assert df.iloc[0]["protein_accessions"] == ["PROT1", "PROT2", "PROT3"]


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

    df = pep_ids.get_psm_df(include_peak_annotations=True)

    assert "number_peaks" in df.columns
    assert "mz_array" in df.columns
    assert "intensity_array" in df.columns
    assert "charge_array" in df.columns
    assert "ion_type_array" in df.columns

    assert df.iloc[0]["number_peaks"] == 2
    assert df.iloc[0]["mz_array"] == [100.5, 200.5]
    assert df.iloc[0]["intensity_array"] == [1000.0, 2000.0]
    assert df.iloc[0]["charge_array"] == [1, 1]
    assert df.iloc[0]["ion_type_array"] == ["y1", "y2"]


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

    df = pep_ids.get_psm_df(include_peak_annotations=True)

    assert df.iloc[0]["number_peaks"] == 0
    assert df.iloc[0]["mz_array"] == []
    assert df.iloc[0]["intensity_array"] == []
