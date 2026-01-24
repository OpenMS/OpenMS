"""
Tests for Universal Spectrum Identifier (USI) bindings.
"""

import unittest

import pyopenms as p


class TestUSIBindings(unittest.TestCase):
    """Test USI and PeptideIdentification USI helpers."""

    def test_build_usi_string_with_interpretation(self):
        pep_id = p.PeptideIdentification()
        pep_id.setSpectrumReference("scan=12345")

        hit = p.PeptideHit()
        hit.setSequence(p.AASequence.fromString("EM(Oxidation)K"))
        hit.setCharge(2)
        hit.setScore(100.0)
        pep_id.insertHit(hit)
        pep_id.sort()

        # Test with explicit dataset_id
        usi = pep_id.buildUSI("sample.mzML", "PXD000561", True)
        self.assertEqual(str(usi.toString()), "mzspec:PXD000561:sample.mzML:scan:12345:EM[UNIMOD:35]K/2")

        # Test with default dataset_id ("local")
        usi_local = pep_id.buildUSI("sample.mzML", "local", False)
        self.assertEqual(str(usi_local.toString()), "mzspec:local:sample.mzML:scan:12345")

