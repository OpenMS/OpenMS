"""Tests for the FAIMSHelper static helpers exposed via pyOpenMS."""

import pyopenms as poms


def _make_spec(cv: float) -> poms.MSSpectrum:
    s = poms.MSSpectrum()
    s.setDriftTime(cv)
    s.setDriftTimeUnit(poms.DriftTimeUnit.FAIMS_COMPENSATION_VOLTAGE)
    return s


def test_class_exists_and_constructible():
    h = poms.FAIMSHelper()
    assert h is not None


def test_get_compensation_voltages_empty():
    exp = poms.MSExperiment()
    assert poms.FAIMSHelper.getCompensationVoltages(exp) == set()


def test_get_compensation_voltages_unique():
    exp = poms.MSExperiment()
    exp.addSpectrum(_make_spec(-45.0))
    exp.addSpectrum(_make_spec(-65.0))
    exp.addSpectrum(_make_spec(-45.0))
    cvs = poms.FAIMSHelper.getCompensationVoltages(exp)
    assert sorted(cvs) == [-65.0, -45.0]


def test_get_compensation_voltages_ignores_non_faims():
    exp = poms.MSExperiment()
    s = poms.MSSpectrum()
    s.setDriftTime(1.23)
    s.setDriftTimeUnit(poms.DriftTimeUnit.MILLISECOND)
    exp.addSpectrum(s)
    exp.addSpectrum(_make_spec(-45.0))
    assert poms.FAIMSHelper.getCompensationVoltages(exp) == {-45.0}


def _make_pid_list():
    pids = poms.PeptideIdentificationList()
    a = poms.PeptideIdentification()
    a.setMetaValue("FAIMS_CV", -45.0)
    b = poms.PeptideIdentification()
    b.setMetaValue("FAIMS_CV", -65.0)
    c = poms.PeptideIdentification()  # no FAIMS_CV — kept for backward compat
    pids.append(a)
    pids.append(b)
    pids.append(c)
    return pids


def test_filter_peptides_by_faims_cv_match():
    pids = _make_pid_list()
    out = poms.FAIMSHelper.filterPeptidesByFAIMSCV(pids, -45.0)
    cvs = [p.getMetaValue("FAIMS_CV") for p in out if p.metaValueExists("FAIMS_CV")]
    assert cvs == [-45.0]
    # untagged ID is kept
    assert sum(1 for p in out if not p.metaValueExists("FAIMS_CV")) == 1


def test_filter_peptides_by_faims_cv_tolerance():
    pids = _make_pid_list()
    # within tolerance
    out = poms.FAIMSHelper.filterPeptidesByFAIMSCV(pids, -45.005, 0.01)
    cvs = [p.getMetaValue("FAIMS_CV") for p in out if p.metaValueExists("FAIMS_CV")]
    assert cvs == [-45.0]
    # outside tolerance — only untagged ID survives
    out2 = poms.FAIMSHelper.filterPeptidesByFAIMSCV(pids, -50.0, 0.01)
    assert all(not p.metaValueExists("FAIMS_CV") for p in out2)
    assert len(out2) == 1
