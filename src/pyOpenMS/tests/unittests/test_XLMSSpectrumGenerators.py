"""SimpleTSGXLMS / TheoreticalSpectrumGeneratorXLMS fragment the peptide passed in.

The bindings used to construct an empty AASequence (or cross-link) internally and
return it alongside an empty spectrum. Expected values are those of the C++ class tests.
"""
import pytest

import pyopenms as oms

PEPTIDE = "IFSQVGK"
LINEAR_MZ = [43.55185, 57.54930, 74.06004, 86.09642, 102.57077, 114.09134, 117.08605,
             131.08351, 147.11280, 152.10497, 160.60207, 174.59953, 204.13426, 233.16484,
             261.15975, 303.20268, 320.19686, 348.19178]
XLINK_MZ = [442.55421, 551.94577, 566.94214, 580.95645, 599.96494, 618.97210, 629.97925,
            661.67042, 661.99842, 663.32768, 667.67394, 827.41502, 849.90957, 870.93103,
            899.44378, 927.95451, 944.46524]


def _by_ions_only(tsg, **extra):
    p = tsg.getParameters()
    for key, value in {"add_isotopes": "false", "max_isotope": 2, "add_a_ions": "false",
                       "add_b_ions": "true", "add_c_ions": "false", "add_x_ions": "false",
                       "add_y_ions": "true", "add_z_ions": "false", "add_losses": "false",
                       **extra}.items():
        p.setValue(key, value)
    tsg.setParameters(p)


class TestSimpleTSGXLMS:
    def test_linear_ion_spectrum(self):
        tsg = oms.SimpleTSGXLMS()
        spec = tsg.getLinearIonSpectrum([], oms.AASequence.fromString(PEPTIDE), 3, 2)
        assert [p.mz for p in spec] == pytest.approx(LINEAR_MZ, abs=1e-3)
        assert len(tsg.getLinearIonSpectrum([], oms.AASequence.fromString(PEPTIDE), 3, 3)) == 27

    def test_linear_ion_spectrum_updates_passed_list(self):
        spec = []
        oms.SimpleTSGXLMS().getLinearIonSpectrum(spec, oms.AASequence.fromString(PEPTIDE), 3, 2)
        assert len(spec) == len(LINEAR_MZ)

    def test_xlink_ion_spectrum_for_peptide(self):
        tsg = oms.SimpleTSGXLMS()
        _by_ions_only(tsg)
        spec = []
        result = tsg.getXLinkIonSpectrum(spec, oms.AASequence.fromString(PEPTIDE), 3, 2000.0, 2, 3)
        assert [p.mz for p in result] == pytest.approx(XLINK_MZ, abs=1e-3)
        assert len(spec) == len(XLINK_MZ)

    def test_xlink_ion_spectrum_without_alpha_peptide_is_empty(self):
        tsg = oms.SimpleTSGXLMS()
        assert tsg.getXLinkIonSpectrum([], oms.ProteinProteinCrossLink(), True, 2, 3) == []


class TestTheoreticalSpectrumGeneratorXLMS:
    def test_linear_ion_spectrum(self):
        spec = oms.MSSpectrum()
        oms.TheoreticalSpectrumGeneratorXLMS().getLinearIonSpectrum(
            spec, oms.AASequence.fromString(PEPTIDE), 3, True, 2)
        mz, _ = spec.get_peaks()
        assert list(mz) == pytest.approx(LINEAR_MZ, abs=1e-3)

    def test_xlink_ion_spectrum_for_peptide(self):
        tsg = oms.TheoreticalSpectrumGeneratorXLMS()
        _by_ions_only(tsg, add_metainfo="false")
        spec = oms.MSSpectrum()
        tsg.getXLinkIonSpectrum(spec, oms.AASequence.fromString(PEPTIDE), 3, 2000.0, True, 2, 3)
        assert spec.size() == len(XLINK_MZ)

    def test_xlink_ion_spectrum_without_alpha_peptide_is_empty(self):
        spec = oms.MSSpectrum()
        oms.TheoreticalSpectrumGeneratorXLMS().getXLinkIonSpectrum(
            spec, oms.ProteinProteinCrossLink(), True, 2, 3)
        assert spec.size() == 0
