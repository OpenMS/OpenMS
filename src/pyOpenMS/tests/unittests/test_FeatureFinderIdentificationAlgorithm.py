from pathlib import Path

import pytest

import pyopenms as oms

TOPP_DATA = Path(__file__).resolve().parents[3] / "tests" / "topp"


def _run(keep):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_input.mzML"), exp)
    proteins, peptides = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(TOPP_DATA / "FeatureFinderIdentification_1_input.idXML"), proteins, peptides)
    algo = oms.FeatureFinderIdentificationAlgorithm()
    params = algo.getParameters()
    params.setValue("extract:mz_window", 0.1)
    params.setValue("extract:batch_size", 10)
    params.setValue("detect:peak_width", 60.0)
    params.setValue("model:type", "none")
    if keep:
        params.setValue("keep_chromatograms", "true")
    algo.setParameters(params)
    algo.setMSData(exp)
    features = oms.FeatureMap()
    algo.run(peptides, proteins, features, oms.FeatureMap(), "")
    return algo, features


@pytest.mark.skipif(not (TOPP_DATA / "FeatureFinderIdentification_1_input.mzML").exists(), reason="needs the TOPP test data")
def test_keep_chromatograms():
    algo, features = _run(keep=True)
    assert features.size() > 0
    chroms = algo.getChromatograms()
    assert chroms.getNrChromatograms() == len(algo.getLibrary().getTransitions())
    assert chroms.getNrChromatograms() > 0

    algo_default, features_default = _run(keep=False)
    assert algo_default.getChromatograms().getNrChromatograms() == 0
    assert features_default.size() == features.size()
