# The identifications of a map as peptide and protein identifications (MapPeptideIdentifications), apart from the map.
from pathlib import Path

import pytest

import pyopenms as oms

TOPP_DATA = Path(__file__).resolve().parents[3] / "tests" / "topp"


def _peptide(identifier, sequence):
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(1.0)
    peptide = oms.PeptideIdentification()
    peptide.setIdentifier(identifier)
    peptide.setScoreType("score")
    peptide.setHits([hit])
    return peptide


def test_view_import_export():
    run = oms.ProteinIdentification()
    run.setIdentifier("search")
    run.setScoreType("score")
    ids = oms.MapPeptideIdentifications()
    assert ids.empty()
    ids.proteins = [run]
    ids.setFeature([0], oms.PeptideIdentificationList([_peptide("search", "PEPTIDE")]))
    ids.setFeature([1, 0], oms.PeptideIdentificationList([_peptide("search", "PEPTIDEK")]))
    ids.unassigned = oms.PeptideIdentificationList([_peptide("search", "PEPTIDER")])
    assert ids.paths() == [[0], [1, 0]]
    assert ids.at([2]).size() == 0

    features = oms.FeatureMap()
    features.push_back(oms.Feature())
    parent = oms.Feature()
    parent.setSubordinates([oms.Feature()])
    features.push_back(parent)
    oms.IdentificationDataConverter.importMapIdentifications(features, ids)
    # the identifications are identification data: three identifications, two of them linked
    runs = features.getIdentificationData().getRuns()
    assert len(runs) == 1 and runs[0].getNumberOfIdentifications() == 3
    assert len(features[0].getIDMatches()) == 1

    exported = oms.IdentificationDataConverter.exportMapIdentifications(features)
    assert exported.paths() == [[0], [1, 0]]
    assert exported.at([1, 0])[0].getHits()[0].getSequence().toString() == "PEPTIDEK"
    assert exported.unassigned[0].getHits()[0].getSequence().toString() == "PEPTIDER"
    assert [p.getIdentifier() for p in exported.proteins] == ["search"]

    exported.eraseFeature([1])
    assert exported.paths() == [[0]]

    features.clearIdentifications()
    assert features.getIdentificationData().empty()
    assert len(features[0].getIDMatches()) == 0


@pytest.mark.skipif(not (TOPP_DATA / "FeatureFinderIdentification_1_output.featureXML").exists(), reason="needs the TOPP test data")
def test_loaded_map_view():
    features = oms.FeatureMap()
    oms.FeatureXMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_output.featureXML"), features)
    ids = oms.IdentificationDataConverter.exportMapIdentifications(features)
    assert len(ids.paths()) > 0
    assert all(ids.at(path).size() > 0 for path in ids.paths())
