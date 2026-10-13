# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------
"""Every PeptideIdentification needs the ProteinIdentification run (search run) its identifier names."""

from pathlib import Path

import pytest

import pyopenms as oms

TOPP_DATA = Path(__file__).resolve().parents[3] / "tests" / "topp"

MESSAGE = "Peptide identification has no matching protein run: 'other'"


def _run(identifier):
    run = oms.ProteinIdentification()  # no protein hits, e.g. a de novo or peptidomics search
    run.setIdentifier(identifier)
    return run


def _peptide(identifier):
    pid = oms.PeptideIdentification()
    pid.setIdentifier(identifier)
    pid.setScoreType("score")  # identification data needs the name of a score
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

    # the identifications of a map, too
    fmap = oms.FeatureMap()
    ids = oms.MapPeptideIdentifications()
    ids.proteins = runs
    ids.unassigned = _peptides("search")
    oms.IdentificationDataConverter.importMapIdentifications(fmap, ids)
    ids.unassigned = _peptides("other")
    with pytest.raises(Exception, match=MESSAGE):
        oms.IdentificationDataConverter.importMapIdentifications(oms.FeatureMap(), ids)


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
    # a peptide identification without its run cannot be an identification of a map (so of a featureXML file)
    fmap = oms.FeatureMap()
    ids = oms.MapPeptideIdentifications()
    ids.proteins = [_run("search")]
    ids.unassigned = _peptides("other")
    with pytest.raises(Exception, match=MESSAGE):
        oms.IdentificationDataConverter.importMapIdentifications(fmap, ids)
    assert fmap.getIdentificationData().empty()


def test_featurexml_store_needs_the_proteins_in_the_run(tmp_path):
    run = _run("search")
    fmap = oms.FeatureMap()
    peptides = _peptides("search")
    hit = peptides[0].getHits()[0]
    hit.setPeptideEvidences([oms.PeptideEvidence("PROT_X", 0, 6, "-", "-")])
    pid = peptides[0]
    pid.setHits([hit])
    ids = oms.MapPeptideIdentifications()
    ids.proteins = [run]
    ids.unassigned = _list(pid)
    oms.IdentificationDataConverter.importMapIdentifications(fmap, ids)
    bad = tmp_path / "bad.featureXML"
    with pytest.raises(Exception, match="No accession PROT_X"):
        oms.FeatureXMLFile().store(str(bad), fmap)
    assert not bad.exists()


def _list(*pids):
    peptides = oms.PeptideIdentificationList()
    for pid in pids:
        peptides.push_back(pid)
    return peptides


@pytest.mark.skipif(not (TOPP_DATA / "FeatureFinderIdentification_1_input.mzML").exists(), reason="needs the TOPP test data")
def test_feature_finder_identification_seeds_without_runs(tmp_path):
    # The pseudo IDs of seeds need a protein identification run; without any (no peptide identifications), one is
    # added for them and removed again with them: the features neither keep pseudo IDs nor an empty placeholder run.
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_input.mzML"), exp)
    seeds = oms.FeatureMap()
    oms.FeatureXMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_output.featureXML"), seeds)
    algo = oms.FeatureFinderIdentificationAlgorithm()
    params = algo.getParameters()
    params.setValue("extract:mz_window", 0.1)
    params.setValue("extract:batch_size", 10)
    params.setValue("detect:peak_width", 60.0)
    params.setValue("model:type", "none")
    algo.setParameters(params)
    algo.setMSData(exp)
    features = oms.FeatureMap()
    algo.run(oms.PeptideIdentificationList(), [], features, seeds, "")
    assert features.size() > 0
    assert oms.IdentificationDataConverter.exportMapIdentifications(features).empty()
    oms.FeatureXMLFile().store(str(tmp_path / "seeds.featureXML"), features)


@pytest.mark.skipif(not (TOPP_DATA / "FeatureFinderIdentification_1_input.mzML").exists(), reason="needs the TOPP test data")
def test_feature_finder_identification_identification_data():
    # FeatureFinderIdentification takes and writes identification data: the features link their identifications, the
    # others stay unassigned in the data. Peptide and protein identifications are converted first, with the same result.
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_input.mzML"), exp)
    proteins, peptides = oms.FileHandler().loadIdentifications(str(TOPP_DATA / "FeatureFinderIdentification_1_input.idXML"))

    def run(*inputs):
        algo = oms.FeatureFinderIdentificationAlgorithm()
        params = algo.getParameters()
        params.setValue("extract:mz_window", 0.1)
        params.setValue("detect:peak_width", 60.0)
        params.setValue("model:type", "none")
        algo.setParameters(params)
        algo.setMSData(exp)
        features = oms.FeatureMap()
        algo.run(*inputs, features)
        return features

    native = run(oms.IdentificationDataAdapter.fromLegacy(proteins, peptides))
    converted = run(peptides, proteins)
    for features in (native, converted):
        assert features.size() > 0
        assert all(len(f.getIDMatches()) > 0 for f in features)
    assert native.size() == converted.size()

    def counts(features):
        return [(run.getNumberOfIdentifications(), run.getNumberOfMatches()) for run in features.getIdentificationData().getRuns()]

    assert counts(native) == counts(converted)
    # each identification keeps its best match only
    assert all(queries == matches for queries, matches in counts(native))
