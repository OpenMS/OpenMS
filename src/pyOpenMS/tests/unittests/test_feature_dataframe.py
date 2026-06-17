"""
Unit tests for QPX feature DataFrame export from ConsensusMap.

Tests the feature export methods:
- to_feature_arrow(): Export features as Arrow Table
- to_feature_df(): Export features as DataFrame
- to_feature_parquet(): Export features to Parquet file
- to_feature_qpx(): Export as QPX format (dict with file_metadata and features)
- feature_columns(): List available columns
"""

import pytest

# Skip entire module if pyarrow is not installed
pytest.importorskip("pyarrow")


def to_list(val):
    """Convert numpy array to list for comparison."""
    if hasattr(val, 'tolist'):
        return val.tolist()
    return val


def is_null(val):
    """Check if value is null (None or NaN)."""
    import pandas as pd
    if val is None:
        return True
    try:
        return pd.isna(val)
    except (TypeError, ValueError):
        return False


def create_test_consensus_map(labelfree=True):
    """Create a ConsensusMap with identified features, column headers, and protein info."""
    import pyopenms as oms

    cmap = oms.ConsensusMap()

    if labelfree:
        cmap.setExperimentType("label-free")
    else:
        cmap.setExperimentType("labeled_MS2")

    # Set up column headers
    ch0 = oms.ColumnHeader()
    ch0.filename = "sample1.mzML"
    ch0.label = "" if labelfree else "126"

    ch1 = oms.ColumnHeader()
    ch1.filename = "sample2.mzML" if labelfree else "sample1.mzML"
    ch1.label = "" if labelfree else "127"

    cmap.setColumnHeaders({0: ch0, 1: ch1})

    # Set up protein identification
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("PI_0")

    ph1 = oms.ProteinHit()
    ph1.setAccession("PROT_A")
    ph1.setScore(0.99)
    ph2 = oms.ProteinHit()
    ph2.setAccession("PROT_B")
    ph2.setScore(0.95)
    prot_id.setHits([ph1, ph2])

    # Add protein group
    pg = oms.ProteinGroup()
    pg.probability = 0.01  # q-value
    pg.accessions = [b"PROT_A", b"PROT_B"]
    prot_id.insertProteinGroup(pg)

    cmap.setProteinIdentifications([prot_id])

    # Create ConsensusFeature 1 - identified
    cf1 = oms.ConsensusFeature()
    cf1.setMZ(500.25)
    cf1.setRT(100.5)
    cf1.setCharge(2)
    cf1.setQuality(0.8)

    # Add feature handles (intensities) using BaseFeature + insert(map_index, base_feature)
    bf1 = oms.BaseFeature()
    bf1.setIntensity(1000.0)
    bf1.setMZ(500.25)
    bf1.setRT(100.5)
    cf1.insert(0, bf1)

    bf2 = oms.BaseFeature()
    bf2.setIntensity(2000.0)
    bf2.setMZ(500.25)
    bf2.setRT(101.0)
    cf1.insert(1, bf2)

    # Add peptide identification
    pep_id1 = oms.PeptideIdentification()
    pep_id1.setRT(100.5)
    pep_id1.setMZ(500.25)
    pep_id1.setScoreType("Posterior Error Probability")
    pep_id1.setMetaValue("spectrum_reference",
                          "controllerType=0 controllerNumber=1 scan=1234")

    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit1.setCharge(2)
    hit1.setScore(0.01)
    hit1.setMetaValue("target_decoy", "target")

    ev1 = oms.PeptideEvidence()
    ev1.setProteinAccession("PROT_A")
    hit1.setPeptideEvidences([ev1])

    pep_id1.setHits([hit1])
    pep_id_list1 = oms.PeptideIdentificationList()
    pep_id_list1.append(pep_id1)
    cf1.setPeptideIdentifications(pep_id_list1)

    cmap.push_back(cf1)

    # Create ConsensusFeature 2 - identified, multiple proteins
    cf2 = oms.ConsensusFeature()
    cf2.setMZ(600.30)
    cf2.setRT(200.0)
    cf2.setCharge(3)
    cf2.setQuality(0.6)

    bf3 = oms.BaseFeature()
    bf3.setIntensity(500.0)
    bf3.setMZ(600.30)
    bf3.setRT(200.0)
    cf2.insert(0, bf3)

    pep_id2 = oms.PeptideIdentification()
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.30)
    pep_id2.setScoreType("Posterior Error Probability")

    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("PEPTIDER"))
    hit2.setCharge(3)
    hit2.setScore(0.05)
    hit2.setMetaValue("target_decoy", "decoy")

    ev2a = oms.PeptideEvidence()
    ev2a.setProteinAccession("PROT_A")
    ev2b = oms.PeptideEvidence()
    ev2b.setProteinAccession("PROT_B")
    hit2.setPeptideEvidences([ev2a, ev2b])

    pep_id2.setHits([hit2])
    pep_id_list2 = oms.PeptideIdentificationList()
    pep_id_list2.append(pep_id2)
    cf2.setPeptideIdentifications(pep_id_list2)

    cmap.push_back(cf2)

    # Create ConsensusFeature 3 - unidentified
    cf3 = oms.ConsensusFeature()
    cf3.setMZ(400.0)
    cf3.setRT(50.0)
    cf3.setCharge(1)
    cf3.setQuality(0.3)

    bf4 = oms.BaseFeature()
    bf4.setIntensity(100.0)
    bf4.setMZ(400.0)
    bf4.setRT(50.0)
    cf3.insert(0, bf4)

    cmap.push_back(cf3)

    return cmap


class TestFeatureArrow:
    """Tests for to_feature_arrow()."""

    def test_basic_export(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        assert table.num_rows == 3
        assert "sequence" in table.column_names
        assert "observed_mz" in table.column_names
        assert "intensities" in table.column_names

    def test_schema_columns(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        expected_cols = cmap.feature_columns()
        for col in expected_cols:
            assert col in table.column_names, f"Missing column: {col}"

    def test_sorted_by_rt(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        rt_col = table.column("rt").to_pylist()
        assert rt_col == sorted(rt_col)

    def test_sequence_values(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        seqs = table.column("sequence").to_pylist()
        # Sorted by RT: cf3(rt=50)=None, cf1(rt=100.5)=PEPTIDE, cf2(rt=200)=PEPTIDER
        assert seqs[0] is None  # unidentified
        assert seqs[1] == "PEPTIDE"
        assert seqs[2] == "PEPTIDER"

    def test_is_decoy(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        decoy = table.column("is_decoy").to_pylist()
        # Sorted by RT: unidentified(None), target(0), decoy(1)
        assert decoy[0] is None
        assert decoy[1] == 0
        assert decoy[2] == 1

    def test_unique(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        unique = table.column("unique").to_pylist()
        # Sorted by RT: unidentified(0), single protein(1), two proteins(0)
        assert unique[0] == 0
        assert unique[1] == 1
        assert unique[2] == 0

    def test_intensities_label_free(self):
        cmap = create_test_consensus_map(labelfree=True)
        table = cmap.to_feature_arrow()
        ints = table.column("intensities").to_pylist()
        # First row (rt=50) has one handle
        assert len(ints[0]) == 1
        assert ints[0][0]["sample_accession"] == "sample1.mzML"
        # Second row (rt=100.5) has two handles
        assert len(ints[1]) == 2

    def test_intensities_labeled(self):
        cmap = create_test_consensus_map(labelfree=False)
        table = cmap.to_feature_arrow()
        ints = table.column("intensities").to_pylist()
        # Second row (rt=100.5) should have channel info
        for entry in ints[1]:
            assert entry["channel"] in ("126", "127")

    def test_pg_accessions(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        pg = table.column("pg_accessions").to_pylist()
        # cf1 (second row after sort) has PROT_A
        assert pg[1] == ["PROT_A"]
        # cf2 (third row) has PROT_A and PROT_B
        assert set(pg[2]) == {"PROT_A", "PROT_B"}

    def test_pg_global_qvalue(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow()
        qvals = table.column("pg_global_qvalue").to_pylist()
        # cf1 -> PROT_A is in group with q=0.01
        assert qvals[1] == pytest.approx(0.01)

    def test_column_filtering(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow(columns=["sequence", "rt", "observed_mz"])
        assert set(table.column_names) == {"sequence", "rt", "observed_mz"}
        assert table.num_rows == 3

    def test_empty_consensus_map(self):
        import pyopenms as oms
        cmap = oms.ConsensusMap()
        table = cmap.to_feature_arrow()
        assert table.num_rows == 0
        # Should still have schema
        assert "sequence" in table.column_names

    def test_reference_file_name(self):
        cmap = create_test_consensus_map()
        table = cmap.to_feature_arrow(reference_file_name="test.mzML")
        refs = table.column("reference_file_name").to_pylist()
        assert all(r == "test.mzML" for r in refs)


class TestFeatureDf:
    """Tests for to_feature_df()."""

    def test_basic_export(self):
        import pandas as pd
        cmap = create_test_consensus_map()
        df = cmap.to_feature_df()
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3

    def test_round_trip(self):
        cmap = create_test_consensus_map()
        df = cmap.to_feature_df()
        table = cmap.to_feature_arrow()
        # Same number of rows and columns
        assert len(df) == table.num_rows
        assert len(df.columns) == len(table.column_names)


class TestFeatureParquet:
    """Tests for to_feature_parquet()."""

    def test_write_read(self, tmp_path):
        import pyarrow.parquet as pq
        cmap = create_test_consensus_map()
        path = str(tmp_path / "test.parquet")
        cmap.to_feature_parquet(path)
        table = pq.read_table(path)
        assert table.num_rows == 3
        assert "sequence" in table.column_names

    def test_compression_snappy(self, tmp_path):
        import pyarrow.parquet as pq
        cmap = create_test_consensus_map()
        path = str(tmp_path / "test.parquet")
        cmap.to_feature_parquet(path, compression='snappy')
        table = pq.read_table(path)
        assert table.num_rows == 3

    def test_compression_none(self, tmp_path):
        import pyarrow.parquet as pq
        cmap = create_test_consensus_map()
        path = str(tmp_path / "test.parquet")
        cmap.to_feature_parquet(path, compression='none')
        table = pq.read_table(path)
        assert table.num_rows == 3


class TestFeatureQPX:
    """Tests for to_feature_qpx()."""

    def test_structure(self):
        cmap = create_test_consensus_map()
        result = cmap.to_feature_qpx()
        assert "file_metadata" in result
        assert "features" in result
        assert result["file_metadata"]["file_type"] == "feature"
        assert result["file_metadata"]["qpx_version"] == "1.0"
        assert len(result["features"]) == 3

    def test_custom_metadata(self):
        cmap = create_test_consensus_map()
        result = cmap.to_feature_qpx(creator="test_tool", software_provider="TestSoft")
        assert result["file_metadata"]["creator"] == "test_tool"
        assert result["file_metadata"]["software_provider"] == "TestSoft"


class TestFeatureColumns:
    """Tests for feature_columns()."""

    def test_returns_list(self):
        import pyopenms as oms
        cmap = oms.ConsensusMap()
        cols = cmap.feature_columns()
        assert isinstance(cols, list)
        assert len(cols) > 0

    def test_contains_required_columns(self):
        import pyopenms as oms
        cmap = oms.ConsensusMap()
        cols = cmap.feature_columns()
        required = ["sequence", "peptidoform", "rt", "observed_mz",
                     "precursor_charge", "intensities", "pg_accessions"]
        for r in required:
            assert r in cols, f"Missing required column: {r}"

    def test_openms_specific_columns(self):
        import pyopenms as oms
        cmap = oms.ConsensusMap()
        cols = cmap.feature_columns()
        assert "quality" in cols
        assert "score" in cols
        assert "feature_metavalues" in cols
