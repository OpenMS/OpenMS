from pathlib import Path

import pytest

import pyopenms as oms

TOPP_DATA = Path(__file__).resolve().parents[3] / "tests" / "topp"
FFID_1 = {"extract:mz_window": 0.1, "extract:batch_size": 10, "detect:peak_width": 60.0, "model:type": "none"}

needs_data = pytest.mark.skipif(not (TOPP_DATA / "FeatureFinderIdentification_1_input.mzML").exists(),
                                reason="needs the TOPP test data")


def _run(mzml, idxml, keep, values=None):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(TOPP_DATA / mzml), exp)
    proteins, peptides = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(TOPP_DATA / idxml), proteins, peptides)
    algo = oms.FeatureFinderIdentificationAlgorithm()
    params = algo.getParameters()
    for key, value in (values or {}).items():
        params.setValue(key, value)
    params.setValue("keep_chromatograms", "true" if keep else "false")
    algo.setParameters(params)
    algo.setMSData(exp)
    features = oms.FeatureMap()
    algo.run(peptides, proteins, features, oms.FeatureMap(), "")
    return algo, features


@needs_data
def test_keep_chromatograms_returns_every_batch():
    algo, features = _run("FeatureFinderIdentification_1_input.mzML", "FeatureFinderIdentification_1_input.idXML",
                          True, FFID_1)
    chroms = algo.getChromatograms().getChromatograms()
    transition_ids = {t.getNativeID() for t in algo.getLibrary().getTransitions()}
    assert len(transition_ids) > 10
    assert sorted(c.getNativeID() for c in chroms) == sorted(transition_ids)
    assert all(c.size() > 0 for c in chroms)

    algo_default, features_default = _run("FeatureFinderIdentification_1_input.mzML",
                                          "FeatureFinderIdentification_1_input.idXML", False, FFID_1)
    assert algo_default.getChromatograms().getNrChromatograms() == 0
    assert features.size() > 0
    assert features.size() == features_default.size()


@pytest.mark.skipif(not (TOPP_DATA / "FAIMS_3CV_synthetic.mzML").exists(), reason="needs the TOPP test data")
def test_keep_chromatograms_faims_ids_are_unique_per_cv():
    algo, _ = _run("FAIMS_3CV_synthetic.mzML", "FAIMS_3CV_synthetic.idXML", True)
    chroms = algo.getChromatograms().getChromatograms()
    ids = [c.getNativeID() for c in chroms]
    assert len(ids) > 0 and len(set(ids)) == len(ids)
    for chrom in chroms:
        cv = chrom.getMetaValue("FAIMS_CV")
        assert chrom.getNativeID().endswith(f"_CV{cv}")
    assert {chrom.getMetaValue("FAIMS_CV") for chrom in chroms} == {-45.0, -55.0, -65.0}


@needs_data
def test_getChromatograms_returns_a_copy():
    algo, _ = _run("FeatureFinderIdentification_1_input.mzML", "FeatureFinderIdentification_1_input.idXML",
                   True, FFID_1)
    n = algo.getChromatograms().getNrChromatograms()
    algo.getChromatograms().clear(True)
    assert algo.getChromatograms().getNrChromatograms() == n
