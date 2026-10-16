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

        # Test with explicit dataset_id (ms_run without extension)
        usi = pep_id.buildUSI("sample", "PXD000561", True)
        self.assertEqual(str(usi.toString()), "mzspec:PXD000561:sample:scan:12345:EM[UNIMOD:35]K/2")

        # Test with default dataset_id ("local")
        usi_local = pep_id.buildUSI("sample", "local", False)
        self.assertEqual(str(usi_local.toString()), "mzspec:local:sample:scan:12345")

    def test_extract_basename_strips_extension(self):
        self.assertEqual(str(p.USI.extractBasename("/path/to/sample.mzML")), "sample")
        self.assertEqual(str(p.USI.extractBasename("sample.raw")), "sample")
        self.assertEqual(str(p.USI.extractBasename("sample")), "sample")

    def test_native_id_to_usi_index(self):
        # Thermo standard -> scan
        result = p.USI.nativeIdToUSIIndex("controllerType=0 controllerNumber=1 scan=43920")
        self.assertIsNotNone(result)
        self.assertEqual(result[0], p.USI.IndexType.SCAN)
        self.assertEqual(str(result[1]), "43920")

        # Bruker TDF -> nativeId with comma-separated values
        result = p.USI.nativeIdToUSIIndex("frame=120 scan=475")
        self.assertIsNotNone(result)
        self.assertEqual(result[0], p.USI.IndexType.NATIVEID)
        self.assertEqual(str(result[1]), "120,475")

        # index= -> index
        result = p.USI.nativeIdToUSIIndex("index=2074")
        self.assertIsNotNone(result)
        self.assertEqual(result[0], p.USI.IndexType.INDEX)
        self.assertEqual(str(result[1]), "2074")

        # Empty -> None
        self.assertIsNone(p.USI.nativeIdToUSIIndex(""))

    def test_build_usi_bruker_tdf(self):
        """Reproduce the issue from #10591: Bruker TDF should use nativeId, not scan."""
        pep_id = p.PeptideIdentification()
        pep_id.setSpectrumReference("frame=120 scan=475")
        hit = p.PeptideHit()
        hit.setSequence(p.AASequence.fromString("EM(Oxidation)K"))
        hit.setCharge(2)
        pep_id.insertHit(hit)
        pep_id.sort()
        usi_str = str(pep_id.buildUSI("sample", "PXD000561", True).toString())
        self.assertEqual(usi_str, "mzspec:PXD000561:sample:nativeId:120,475:EM[UNIMOD:35]K/2")
