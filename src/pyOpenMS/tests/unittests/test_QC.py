"""Smoke tests for OpenMS Quality Control (QC) Python bindings.

Run directly (from OUTSIDE the source tree to avoid shadowing):
    PYTHONPATH=/path/to/OpenMS-build/pyOpenMS python3 src/pyOpenMS/tests/unittests/test_QC.py
"""

import copy
import pyopenms as oms


def test_import_qc():
    """Verify all expected symbols are exported from the qc module."""
    # Enums
    assert hasattr(oms, "QCBaseRequires")
    assert hasattr(oms, "QCBaseToleranceUnit")
    # Flag set / status
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


def test_qcbase_status_default():
    """QCBaseStatus should start empty."""
    s = oms.QCBaseStatus()
    assert s.empty()
    assert s.value() == 0


def test_qcbase_status_from_enum():
    """Constructing from enum should set the corresponding bit."""
    s = oms.QCBaseStatus(oms.QCBaseRequires.RAWMZML)
    assert not s.empty()
    assert s.value() > 0
    assert s.isSuperSetOf(oms.QCBaseRequires.RAWMZML)


def test_qcbase_status_bitwise_or():
    """Bitwise OR should combine two flags."""
    s_raw = oms.QCBaseStatus(oms.QCBaseRequires.RAWMZML)
    s_id  = oms.QCBaseStatus(oms.QCBaseRequires.ID)
    s_both = s_raw | s_id
    assert s_both.isSuperSetOf(oms.QCBaseRequires.RAWMZML)
    assert s_both.isSuperSetOf(oms.QCBaseRequires.ID)


def test_qcbase_status_bitwise_and():
    """Bitwise AND should mask flags."""
    s_both = (oms.QCBaseStatus(oms.QCBaseRequires.RAWMZML)
              | oms.QCBaseRequires.ID)
    s_masked = s_both & oms.QCBaseRequires.ID
    assert s_masked.isSuperSetOf(oms.QCBaseRequires.ID)
    assert not s_masked.isSuperSetOf(oms.QCBaseRequires.RAWMZML)


def test_qcbase_status_inplace_or():
    """In-place |= should mutate the object."""
    s = oms.QCBaseStatus()
    s |= oms.QCBaseRequires.RAWMZML
    assert not s.empty()
    assert s.isSuperSetOf(oms.QCBaseRequires.RAWMZML)


def test_qcbase_status_copy():
    """QCBaseStatus should support copy.copy and copy.deepcopy."""
    s = oms.QCBaseStatus(oms.QCBaseRequires.RAWMZML)
    s2 = copy.copy(s)
    s3 = copy.deepcopy(s)
    assert s2.isSuperSetOf(oms.QCBaseRequires.RAWMZML)
    assert s3.isSuperSetOf(oms.QCBaseRequires.RAWMZML)


def test_spectra_map_empty():
    """Freshly created SpectraMap should be empty."""
    m = oms.QCSpectraMap()
    assert m.empty()
    assert m.size() == 0


def test_spectra_map_copy():
    """QCSpectraMap should support copy.copy."""
    m = oms.QCSpectraMap()
    m2 = copy.copy(m)
    assert m2.empty()


def test_tic_name_and_requirements():
    """TIC should advertise its name and require a raw mzML."""
    tic = oms.TIC()
    assert tic.getName() == "TIC"
    req = tic.requirements()
    assert req.isSuperSetOf(oms.QCBaseRequires.RAWMZML)


def test_tic_compute_empty():
    """TIC.compute on an empty MSExperiment should return an empty result."""
    tic = oms.TIC()
    result = tic.compute(oms.MSExperiment(), 0.0, 1)
    assert result.area == 0.0
    assert len(result.intensities) == 0
    assert len(result.retention_times) == 0


def test_tic_copy():
    """TICResult and TIC should support copy."""
    result = oms.TICResult()
    result2 = copy.copy(result)
    assert result == result2

    tic = oms.TIC()
    tic2 = copy.copy(tic)
    assert tic2.getName() == "TIC"


def test_spectrum_count_name():
    """SpectrumCount should report the correct name."""
    sc = oms.SpectrumCount()
    assert sc.getName() == "SpectrumCount"


def test_spectrum_count_empty():
    """SpectrumCount.compute on empty experiment should return empty map."""
    sc = oms.SpectrumCount()
    result = sc.compute(oms.MSExperiment())
    assert len(result) == 0


def test_missed_cleavages_empty_results():
    """MissedCleavages.getResults() should be empty initially."""
    mc = oms.MissedCleavages()
    assert mc.getName() == "MissedCleavages"
    assert len(mc.getResults()) == 0


def test_missed_cleavages_copy():
    """MissedCleavages should support copy."""
    mc = oms.MissedCleavages()
    mc2 = copy.copy(mc)
    assert mc2.getName() == "MissedCleavages"


def test_fragment_mass_error_statistics_defaults():
    """FragmentMassErrorStatistics should have zero defaults."""
    stats = oms.FragmentMassErrorStatistics()
    assert stats.average_ppm == 0.0
    assert stats.variance_ppm == 0.0


def test_contaminants_empty_results():
    """Contaminants.getResults() should be empty initially."""
    cont = oms.Contaminants()
    assert cont.getName() == "Contaminants"
    assert len(cont.getResults()) == 0


def test_db_suitability_empty_results():
    """DBSuitability.getResults() should be empty initially."""
    db = oms.DBSuitability()
    # DBSuitability inherits DefaultParamHandler, not QCBase
    assert hasattr(db, "getResults")
    assert len(db.getResults()) == 0


def test_db_suitability_data_defaults():
    """DBSuitabilityData fields should have correct defaults."""
    data = oms.DBSuitabilityData()
    assert data.num_top_novo == 0
    assert data.num_top_db == 0
    assert data.suitability == 0.0


def test_feature_summary_requirements():
    """FeatureSummary should require POSTFDRFEAT."""
    fs = oms.FeatureSummary()
    req = fs.requirements()
    assert req.isSuperSetOf(oms.QCBaseRequires.POSTFDRFEAT)


def test_feature_summary_result_defaults():
    """FeatureSummaryResult should have zero defaults and support equality."""
    r = oms.FeatureSummaryResult()
    assert r.feature_count == 0
    r2 = copy.copy(r)
    assert r == r2


def test_identification_summary_unique_id_defaults():
    """UniqueID should have zero defaults."""
    uid = oms.UniqueID()
    assert uid.count == 0
    assert uid.fdr_threshold == -1.0


def test_ms2_identification_rate_empty_results():
    """Ms2IdentificationRate.getResults() should be empty initially."""
    rate = oms.Ms2IdentificationRate()
    assert rate.getName() == "Ms2IdentificationRate"
    assert len(rate.getResults()) == 0


def test_identification_rate_data_defaults():
    """IdentificationRateData should default to zeros."""
    d = oms.IdentificationRateData()
    assert d.num_peptide_identification == 0
    assert d.num_ms2_spectra == 0
    assert d.identification_rate == 0.0


def test_ms2_spectrum_stats_scan_event():
    """Ms2SpectrumStatsScanEvent should store fields."""
    se = oms.Ms2SpectrumStatsScanEvent(3, True)
    assert se.scan_event_number == 3
    assert se.ms2_presence is True
    se2 = copy.copy(se)
    assert se2.scan_event_number == 3


def test_psm_explained_ion_current_statistics_defaults():
    """PSMExplainedIonCurrentStatistics should default to zeros."""
    stats = oms.PSMExplainedIonCurrentStatistics()
    assert stats.average_correctness == 0.0
    assert stats.variance_correctness == 0.0
    stats2 = copy.copy(stats)
    assert stats2.average_correctness == 0.0


def test_psm_explained_ion_current_empty_results():
    """PSMExplainedIonCurrent.getResults() should be empty initially."""
    psm = oms.PSMExplainedIonCurrent()
    assert psm.getName() == "PSMExplainedIonCurrent"
    assert len(psm.getResults()) == 0


if __name__ == "__main__":
    tests = [
        test_import_qc,
        test_qcbase_status_default,
        test_qcbase_status_from_enum,
        test_qcbase_status_bitwise_or,
        test_qcbase_status_bitwise_and,
        test_qcbase_status_inplace_or,
        test_qcbase_status_copy,
        test_spectra_map_empty,
        test_spectra_map_copy,
        test_tic_name_and_requirements,
        test_tic_compute_empty,
        test_tic_copy,
        test_spectrum_count_name,
        test_spectrum_count_empty,
        test_missed_cleavages_empty_results,
        test_missed_cleavages_copy,
        test_fragment_mass_error_statistics_defaults,
        test_contaminants_empty_results,
        test_db_suitability_empty_results,
        test_db_suitability_data_defaults,
        test_feature_summary_requirements,
        test_feature_summary_result_defaults,
        test_identification_summary_unique_id_defaults,
        test_ms2_identification_rate_empty_results,
        test_identification_rate_data_defaults,
        test_ms2_spectrum_stats_scan_event,
        test_psm_explained_ion_current_statistics_defaults,
        test_psm_explained_ion_current_empty_results,
    ]
    for t in tests:
        t()
        print(f"{t.__name__}: passed")
    print("\nAll QC pyOpenMS tests passed.")
