"""Smoke tests for OpenMS Quality Control (QC) bindings.

Run directly:
    PYTHONPATH=OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/unittests/test_QC.py
"""

import pyopenms as oms


def test_import_qc():
    # Enums & Status
    assert hasattr(oms, "QCBaseRequires")
    assert hasattr(oms, "QCBaseToleranceUnit")
    assert hasattr(oms, "QCBaseStatus")
    assert hasattr(oms, "QCSpectraMap")
    assert hasattr(oms, "QCBase")

    # Metrics
    assert hasattr(oms, "FWHM")
    assert hasattr(oms, "TIC")
    assert hasattr(oms, "TICResult")
    assert hasattr(oms, "MissedCleavages")
    assert hasattr(oms, "FragmentMassError")
    assert hasattr(oms, "FragmentMassErrorStatistics")
    assert hasattr(oms, "SpectrumCount")
    assert hasattr(oms, "PeptideMass")
    assert hasattr(oms, "RTAlignment")
    assert hasattr(oms, "MzCalibration")
    assert hasattr(oms, "Contaminants")
    assert hasattr(oms, "ContaminantsSummary")
    assert hasattr(oms, "DBSuitability")
    assert hasattr(oms, "DBSuitabilityData")
    assert hasattr(oms, "FeatureSummary")
    assert hasattr(oms, "FeatureSummaryResult")
    assert hasattr(oms, "IdentificationSummary")
    assert hasattr(oms, "IdentificationSummaryResult")
    assert hasattr(oms, "UniqueID")
    assert hasattr(oms, "Ms2IdentificationRate")
    assert hasattr(oms, "IdentificationRateData")
    assert hasattr(oms, "Ms2SpectrumStats")
    assert hasattr(oms, "Ms2SpectrumStatsScanEvent")
    assert hasattr(oms, "PSMExplainedIonCurrent")
    assert hasattr(oms, "PSMExplainedIonCurrentStatistics")


def test_qcbase_status():
    # Default CTor
    status = oms.QCBaseStatus()
    assert status.empty()

    # Enum CTor
    status_raw = oms.QCBaseStatus(oms.QCBaseRequires.RAWMZML)
    assert not status_raw.empty()
    assert status_raw.value() > 0

    # Bitwise operations
    status_both = status_raw | oms.QCBaseRequires.ID
    assert status_both.isSuperSetOf(oms.QCBaseRequires.RAWMZML)
    assert status_both.isSuperSetOf(oms.QCBaseRequires.ID)

    status_both &= oms.QCBaseRequires.ID
    assert not status_both.isSuperSetOf(oms.QCBaseRequires.RAWMZML)
    assert status_both.isSuperSetOf(oms.QCBaseRequires.ID)


def test_spectra_map():
    map = oms.QCSpectraMap()
    assert map.empty()
    assert map.size() == 0


def test_tic():
    tic = oms.TIC()
    assert tic.getName() == "TIC"
    req = tic.requirements()
    assert req.isSuperSetOf(oms.QCBaseRequires.RAWMZML)

    exp = oms.MSExperiment()
    result = tic.compute(exp, 0.0, 1)
    assert result.area == 0
    assert len(result.intensities) == 0


def test_spectrum_count():
    sc = oms.SpectrumCount()
    assert sc.getName() == "SpectrumCount"
    exp = oms.MSExperiment()
    counts = sc.compute(exp)
    assert len(counts) == 0


def test_missed_cleavages():
    mc = oms.MissedCleavages()
    assert mc.getName() == "MissedCleavages"
    res = mc.getResults()
    assert len(res) == 0


def test_contaminants():
    cont = oms.Contaminants()
    assert cont.getName() == "Contaminants"
    res = cont.getResults()
    assert len(res) == 0


def test_dbsuitability():
    db = oms.DBSuitability()
    assert db.getName() == "DBSuitability"
    res = db.getResults()
    assert len(res) == 0


def test_feature_summary():
    fs = oms.FeatureSummary()
    assert fs.getName() == "Summary of features from featureXML file"
    req = fs.requirements()
    assert req.isSuperSetOf(oms.QCBaseRequires.POSTFDRFEAT)


if __name__ == "__main__":
    test_import_qc()
    print("test_import_qc: passed")
    test_qcbase_status()
    print("test_qcbase_status: passed")
    test_spectra_map()
    print("test_spectra_map: passed")
    test_tic()
    print("test_tic: passed")
    test_spectrum_count()
    print("test_spectrum_count: passed")
    test_missed_cleavages()
    print("test_missed_cleavages: passed")
    test_contaminants()
    print("test_contaminants: passed")
    test_dbsuitability()
    print("test_dbsuitability: passed")
    test_feature_summary()
    print("test_feature_summary: passed")
    print("\nAll QC pyOpenMS tests passed.")
