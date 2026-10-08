# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Every PeptideIdentification needs the ProteinIdentification run (search run) its identifier names."""

import pytest

import pyopenms as oms

MESSAGE = "Peptide identification has no matching protein run: 'other'"


def _run(identifier):
    run = oms.ProteinIdentification()  # no protein hits, e.g. a de novo or peptidomics search
    run.setIdentifier(identifier)
    return run


def _peptide(identifier):
    pid = oms.PeptideIdentification()
    pid.setIdentifier(identifier)
    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString("PEPTIDEK"))
    pid.setHits([hit])
    return pid


def _peptides(*identifiers):
    peptides = oms.PeptideIdentificationList()
    for identifier in identifiers:
        peptides.push_back(_peptide(identifier))
    return peptides


def test_check():
    runs = [_run("search")]
    assert oms.ProteinRunReferences.missingRunMessage("other").startswith(MESSAGE)
    oms.ProteinRunReferences.check(runs, "search")
    oms.ProteinRunReferences.check(runs, _peptides("search"))
    with pytest.raises(Exception, match=MESSAGE):
        oms.ProteinRunReferences.check(runs, "other")
    with pytest.raises(Exception, match=MESSAGE):
        oms.ProteinRunReferences.check(runs, _peptides("search", "other"))

    fmap = oms.FeatureMap()
    fmap.setProteinIdentifications(runs)
    fmap.setUnassignedPeptideIdentifications(_peptides("search"))
    oms.ProteinRunReferences.check(fmap)
    fmap.setUnassignedPeptideIdentifications(_peptides("other"))
    with pytest.raises(Exception, match=MESSAGE):
        oms.ProteinRunReferences.check(fmap)


def test_idxml_store_does_not_omit(tmp_path):
    runs = [_run("search")]
    good = str(tmp_path / "good.idXML")
    oms.IdXMLFile().store(good, runs, _peptides("search"))
    runs_in, peptides_in = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(good, runs_in, peptides_in)
    assert len(runs_in) == 1 and peptides_in.size() == 1
    assert peptides_in[0].getIdentifier() == runs_in[0].getIdentifier()

    # a peptide identification without its run is not silently dropped: nothing is written
    bad = tmp_path / "bad.idXML"
    with pytest.raises(Exception, match=MESSAGE):
        oms.IdXMLFile().store(str(bad), runs, _peptides("search", "other"))
    assert not bad.exists()


def test_featurexml_store_does_not_omit(tmp_path):
    fmap = oms.FeatureMap()
    fmap.setProteinIdentifications([_run("search")])
    fmap.setUnassignedPeptideIdentifications(_peptides("other"))
    bad = tmp_path / "bad.featureXML"
    with pytest.raises(Exception, match=MESSAGE):
        oms.FeatureXMLFile().store(str(bad), fmap)
    assert not bad.exists()
