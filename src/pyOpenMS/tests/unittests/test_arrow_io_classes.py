"""
Regression tests for Arrow IO nanobind bindings.

Tests the Parquet and zero-copy Arrow export/import for:
- FeatureMapArrowIO (bind_format.cpp + arrow_zerocopy.cpp)
- ConsensusMapArrowIO (bind_format.cpp + arrow_zerocopy.cpp)
- ProteinIdentificationArrowIO (bind_format.cpp + arrow_zerocopy.cpp)

Also verifies that CrossLinksDB.readFromOBOFile was removed.
"""

import pytest
import tempfile
import os

pa = pytest.importorskip("pyarrow")
import pyarrow.parquet as pq

import pyopenms as oms

# Skip if the Arrow IO classes were not built
if not hasattr(oms, "FeatureMapArrowIO"):
    pytest.skip("Arrow IO classes not available",
                allow_module_level=True)

# Import zero-copy functions
from pyopenms._arrow_zerocopy import (
    featuremap_features_to_arrow,
    featuremap_psms_to_arrow,
    featuremap_import_features_from_arrow,
    featuremap_import_psms_from_arrow,
    consensusmap_features_to_arrow,
    consensusmap_psms_to_arrow,
    consensusmap_import_features_from_arrow,
    consensusmap_import_psms_from_arrow,
    protein_ids_proteins_to_arrow,
    protein_ids_groups_to_arrow,
    protein_ids_search_params_to_arrow,
    protein_ids_import_search_params_from_arrow,
    protein_ids_import_proteins_from_arrow,
    protein_ids_import_groups_from_arrow,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def feature_map():
    """FeatureMap with 2 features and peptide identifications."""
    fmap = oms.FeatureMap()

    # Protein identification
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("PI_0")
    ph = oms.ProteinHit()
    ph.setAccession("PROT_A")
    ph.setScore(0.99)
    prot_id.setHits([ph])
    fmap.setProteinIdentifications([prot_id])

    # Feature 1 with PSM
    f1 = oms.Feature()
    f1.setRT(100.0)
    f1.setMZ(500.25)
    f1.setIntensity(1000.0)
    f1.setCharge(2)

    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.0)
    pep_id.setMZ(500.25)
    pep_id.setIdentifier("PI_0")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    ev = oms.PeptideEvidence()
    ev.setProteinAccession("PROT_A")
    hit.setPeptideEvidences([ev])
    pep_id.setHits([hit])
    pid_list = oms.PeptideIdentificationList()
    pid_list.append(pep_id)
    f1.setPeptideIdentifications(pid_list)
    fmap.push_back(f1)

    # Feature 2 without PSM
    f2 = oms.Feature()
    f2.setRT(200.0)
    f2.setMZ(600.30)
    f2.setIntensity(2000.0)
    f2.setCharge(3)
    fmap.push_back(f2)

    return fmap


@pytest.fixture
def consensus_map():
    """ConsensusMap with 2 consensus features and column headers."""
    cmap = oms.ConsensusMap()
    cmap.setExperimentType("label-free")

    ch0 = oms.ColumnHeader()
    ch0.filename = "sample1.mzML"
    ch0.label = ""
    cmap.setColumnHeaders({0: ch0})

    # Protein identification
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("PI_0")
    ph = oms.ProteinHit()
    ph.setAccession("PROT_A")
    ph.setScore(0.99)
    prot_id.setHits([ph])
    cmap.setProteinIdentifications([prot_id])

    # ConsensusFeature 1
    cf1 = oms.ConsensusFeature()
    cf1.setMZ(500.25)
    cf1.setRT(100.5)
    cf1.setCharge(2)
    cf1.setIntensity(1000.0)

    bf1 = oms.BaseFeature()
    bf1.setIntensity(1000.0)
    bf1.setMZ(500.25)
    bf1.setRT(100.5)
    cf1.insert(0, bf1)

    pep_id = oms.PeptideIdentification()
    pep_id.setRT(100.5)
    pep_id.setMZ(500.25)
    pep_id.setIdentifier("PI_0")
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
    hit.setCharge(2)
    hit.setScore(0.01)
    pep_id.setHits([hit])
    pid_list = oms.PeptideIdentificationList()
    pid_list.append(pep_id)
    cf1.setPeptideIdentifications(pid_list)
    cmap.push_back(cf1)

    # ConsensusFeature 2
    cf2 = oms.ConsensusFeature()
    cf2.setMZ(600.30)
    cf2.setRT(200.0)
    cf2.setCharge(3)
    cf2.setIntensity(2000.0)

    bf2 = oms.BaseFeature()
    bf2.setIntensity(2000.0)
    bf2.setMZ(600.30)
    bf2.setRT(200.0)
    cf2.insert(0, bf2)
    cmap.push_back(cf2)

    return cmap


@pytest.fixture
def protein_identifications():
    """List of ProteinIdentification with hits, groups, and search parameters."""
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("PI_0")
    prot_id.setSearchEngine("Comet")
    prot_id.setScoreType("q-value")

    # Search parameters
    sp = oms.SearchParameters()
    sp.db = "uniprot_human.fasta"
    sp.missed_cleavages = 2
    prot_id.setSearchParameters(sp)

    # Protein hits
    ph1 = oms.ProteinHit()
    ph1.setAccession("PROT_A")
    ph1.setScore(0.99)
    ph2 = oms.ProteinHit()
    ph2.setAccession("PROT_B")
    ph2.setScore(0.95)
    prot_id.setHits([ph1, ph2])

    # Protein group
    pg = oms.ProteinGroup()
    pg.probability = 0.01
    pg.accessions = [b"PROT_A", b"PROT_B"]
    prot_id.insertProteinGroup(pg)

    return [prot_id]


# ---------------------------------------------------------------------------
# CrossLinksDB — removed method regression
# ---------------------------------------------------------------------------

class TestCrossLinksDBRemoval:
    def test_readFromOBOFile_removed(self):
        """readFromOBOFile was removed from C++ in the ModificationsDB refactoring."""
        db = oms.CrossLinksDB.getInstance()
        assert not hasattr(db, "readFromOBOFile"), \
            "CrossLinksDB.readFromOBOFile should have been removed"


# ---------------------------------------------------------------------------
# FeatureMapArrowIO — Parquet round-trip
# ---------------------------------------------------------------------------

class TestFeatureMapArrowIOParquet:

    def test_parquet_round_trip(self, feature_map):
        """Export and re-import FeatureMap via Parquet directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            ok = oms.FeatureMapArrowIO.exportToParquet(feature_map, tmpdir)
            assert ok

            fmap2 = oms.FeatureMap()
            ok = oms.FeatureMapArrowIO.importFromParquet(tmpdir, fmap2)
            assert ok
            assert fmap2.size() == 2

    def test_parquet_preserves_rt_mz(self, feature_map):
        with tempfile.TemporaryDirectory() as tmpdir:
            oms.FeatureMapArrowIO.exportToParquet(feature_map, tmpdir)
            fmap2 = oms.FeatureMap()
            oms.FeatureMapArrowIO.importFromParquet(tmpdir, fmap2)

            assert fmap2[0].getRT() == pytest.approx(100.0)
            assert fmap2[0].getMZ() == pytest.approx(500.25)
            assert fmap2[1].getRT() == pytest.approx(200.0)
            assert fmap2[1].getMZ() == pytest.approx(600.30)


# ---------------------------------------------------------------------------
# FeatureMapArrowIO — zero-copy Arrow export/import
# ---------------------------------------------------------------------------

class TestFeatureMapArrowZerocopy:

    def test_features_to_arrow_returns_table(self, feature_map):
        table = featuremap_features_to_arrow(feature_map)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 2

    def test_features_to_arrow_columns(self, feature_map):
        table = featuremap_features_to_arrow(feature_map)
        names = table.column_names
        for col in ["rt", "mz", "intensity", "charge"]:
            assert col in names, f"Missing column: {col}"

    def test_features_to_arrow_values(self, feature_map):
        table = featuremap_features_to_arrow(feature_map)
        rt_col = table.column("rt").to_pylist()
        assert rt_col[0] == pytest.approx(100.0)
        assert rt_col[1] == pytest.approx(200.0)

    def test_psms_to_arrow_returns_table(self, feature_map):
        table = featuremap_psms_to_arrow(feature_map)
        assert isinstance(table, pa.Table)
        # Feature 1 has 1 PSM, Feature 2 has 0
        assert table.num_rows >= 1

    def test_features_arrow_round_trip(self, feature_map):
        """Export features to Arrow, import back, verify data preserved."""
        table = featuremap_features_to_arrow(feature_map)

        fmap2 = oms.FeatureMap()
        ok = featuremap_import_features_from_arrow(table, fmap2)
        assert ok
        assert fmap2.size() == 2
        assert fmap2[0].getRT() == pytest.approx(100.0)
        assert fmap2[0].getMZ() == pytest.approx(500.25)

    def test_empty_featuremap(self):
        fmap = oms.FeatureMap()
        table = featuremap_features_to_arrow(fmap)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 0


# ---------------------------------------------------------------------------
# ConsensusMapArrowIO — Parquet round-trip
# ---------------------------------------------------------------------------

class TestConsensusMapArrowIOParquet:

    def test_parquet_round_trip(self, consensus_map):
        with tempfile.TemporaryDirectory() as tmpdir:
            ok = oms.ConsensusMapArrowIO.exportToParquet(consensus_map, tmpdir)
            assert ok

            cmap2 = oms.ConsensusMap()
            ok = oms.ConsensusMapArrowIO.importFromParquet(tmpdir, cmap2)
            assert ok
            assert cmap2.size() == 2

    def test_parquet_preserves_rt_mz(self, consensus_map):
        with tempfile.TemporaryDirectory() as tmpdir:
            oms.ConsensusMapArrowIO.exportToParquet(consensus_map, tmpdir)
            cmap2 = oms.ConsensusMap()
            oms.ConsensusMapArrowIO.importFromParquet(tmpdir, cmap2)

            assert cmap2[0].getRT() == pytest.approx(100.5)
            assert cmap2[0].getMZ() == pytest.approx(500.25)


# ---------------------------------------------------------------------------
# ConsensusMapArrowIO — zero-copy Arrow export/import
# ---------------------------------------------------------------------------

class TestConsensusMapArrowZerocopy:

    def test_features_to_arrow_returns_table(self, consensus_map):
        table = consensusmap_features_to_arrow(consensus_map)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 2

    def test_features_to_arrow_columns(self, consensus_map):
        table = consensusmap_features_to_arrow(consensus_map)
        names = table.column_names
        for col in ["rt", "mz", "intensity", "charge"]:
            assert col in names, f"Missing column: {col}"

    def test_psms_to_arrow_returns_table(self, consensus_map):
        table = consensusmap_psms_to_arrow(consensus_map)
        assert isinstance(table, pa.Table)
        assert table.num_rows >= 1

    def test_features_arrow_round_trip(self, consensus_map):
        table = consensusmap_features_to_arrow(consensus_map)

        cmap2 = oms.ConsensusMap()
        ok = consensusmap_import_features_from_arrow(table, cmap2)
        assert ok
        assert cmap2.size() == 2
        assert cmap2[0].getRT() == pytest.approx(100.5)

    def test_empty_consensusmap(self):
        cmap = oms.ConsensusMap()
        table = consensusmap_features_to_arrow(cmap)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 0


# ---------------------------------------------------------------------------
# ProteinIdentificationArrowIO — Parquet round-trip
# ---------------------------------------------------------------------------

class TestProteinIdentificationArrowIOParquet:

    def test_export_proteins_to_parquet(self, protein_identifications, tmp_path):
        path = str(tmp_path / "proteins.parquet")
        ok = oms.ProteinIdentificationArrowIO.exportProteinsToParquet(
            protein_identifications, path)
        assert ok
        assert os.path.exists(path)
        table = pq.read_table(path)
        assert table.num_rows == 2  # 2 protein hits

    def test_export_groups_to_parquet(self, protein_identifications, tmp_path):
        path = str(tmp_path / "groups.parquet")
        ok = oms.ProteinIdentificationArrowIO.exportProteinGroupsToParquet(
            protein_identifications, path)
        assert ok
        table = pq.read_table(path)
        assert table.num_rows >= 1

    def test_export_search_params_to_parquet(self, protein_identifications, tmp_path):
        path = str(tmp_path / "params.parquet")
        ok = oms.ProteinIdentificationArrowIO.exportSearchParamsToParquet(
            protein_identifications, path)
        assert ok
        table = pq.read_table(path)
        assert table.num_rows == 1  # 1 search run

    def test_parquet_full_round_trip(self, protein_identifications):
        """Export all three files, reimport via importFromParquet."""
        with tempfile.TemporaryDirectory() as tmpdir:
            prot_f = os.path.join(tmpdir, "proteins.parquet")
            groups_f = os.path.join(tmpdir, "groups.parquet")
            params_f = os.path.join(tmpdir, "params.parquet")

            oms.ProteinIdentificationArrowIO.exportProteinsToParquet(
                protein_identifications, prot_f)
            oms.ProteinIdentificationArrowIO.exportProteinGroupsToParquet(
                protein_identifications, groups_f)
            oms.ProteinIdentificationArrowIO.exportSearchParamsToParquet(
                protein_identifications, params_f)

            result = oms.ProteinIdentificationArrowIO.importFromParquet(
                prot_f, groups_f, params_f)
            assert len(result) == 1
            assert len(result[0].getHits()) == 2
            assert result[0].getHits()[0].getAccession() == "PROT_A"

    def test_import_returns_list(self, protein_identifications):
        """Verify import methods return lists (not silently discard output)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            params_f = os.path.join(tmpdir, "params.parquet")
            prot_f = os.path.join(tmpdir, "proteins.parquet")

            oms.ProteinIdentificationArrowIO.exportSearchParamsToParquet(
                protein_identifications, params_f)
            oms.ProteinIdentificationArrowIO.exportProteinsToParquet(
                protein_identifications, prot_f)

            # Step 1: import search params (creates ProteinIdentification shells)
            prot_ids = oms.ProteinIdentificationArrowIO.importSearchParamsFromParquet(params_f)
            assert isinstance(prot_ids, list)
            assert len(prot_ids) == 1

            # Step 2: import proteins into existing identifications
            prot_ids = oms.ProteinIdentificationArrowIO.importProteinsFromParquet(prot_f, prot_ids)
            assert isinstance(prot_ids, list)
            assert len(prot_ids[0].getHits()) == 2


# ---------------------------------------------------------------------------
# ProteinIdentificationArrowIO — zero-copy Arrow export/import
# ---------------------------------------------------------------------------

class TestProteinIdentificationArrowZerocopy:

    def test_proteins_to_arrow(self, protein_identifications):
        table = protein_ids_proteins_to_arrow(protein_identifications)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 2

    def test_groups_to_arrow(self, protein_identifications):
        table = protein_ids_groups_to_arrow(protein_identifications)
        assert isinstance(table, pa.Table)
        assert table.num_rows >= 1

    def test_search_params_to_arrow(self, protein_identifications):
        table = protein_ids_search_params_to_arrow(protein_identifications)
        assert isinstance(table, pa.Table)
        assert table.num_rows == 1

    def test_search_params_round_trip(self, protein_identifications):
        """Export search params to Arrow, reimport, verify."""
        table = protein_ids_search_params_to_arrow(protein_identifications)
        result = protein_ids_import_search_params_from_arrow(table)
        assert isinstance(result, list)
        assert len(result) == 1

    def test_proteins_round_trip(self, protein_identifications):
        """Export proteins to Arrow, reimport into existing identifications."""
        params_table = protein_ids_search_params_to_arrow(protein_identifications)
        prot_ids = protein_ids_import_search_params_from_arrow(params_table)

        proteins_table = protein_ids_proteins_to_arrow(protein_identifications)
        prot_ids = protein_ids_import_proteins_from_arrow(proteins_table, prot_ids)
        assert len(prot_ids) == 1
        assert len(prot_ids[0].getHits()) == 2
        assert prot_ids[0].getHits()[0].getAccession() == "PROT_A"

    def test_groups_round_trip(self, protein_identifications):
        """Full three-step Arrow round-trip: params → proteins → groups."""
        params_table = protein_ids_search_params_to_arrow(protein_identifications)
        prot_ids = protein_ids_import_search_params_from_arrow(params_table)

        proteins_table = protein_ids_proteins_to_arrow(protein_identifications)
        prot_ids = protein_ids_import_proteins_from_arrow(proteins_table, prot_ids)

        groups_table = protein_ids_groups_to_arrow(protein_identifications)
        prot_ids = protein_ids_import_groups_from_arrow(groups_table, prot_ids)
        assert len(prot_ids) == 1
        assert len(prot_ids[0].getProteinGroups()) >= 1
