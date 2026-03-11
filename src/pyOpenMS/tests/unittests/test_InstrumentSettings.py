import unittest

import pyopenms


class TestInstrumentSettings(unittest.TestCase):

    def test_scanModeToString(self):
        """Test converting ScanMode enum to string."""
        SM = pyopenms.InstrumentSettings.ScanMode
        # Test known scan modes
        self.assertEqual(
            pyopenms.InstrumentSettings.scanModeToString(SM.MS1SPECTRUM),
            "MS1Spectrum"
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.scanModeToString(SM.MSNSPECTRUM),
            "MSnSpectrum"
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.scanModeToString(SM.SIM),
            "SelectedIonMonitoring"
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.scanModeToString(SM.SRM),
            "SelectedReactionMonitoring"
        )

    def test_toScanMode(self):
        """Test converting string to ScanMode enum."""
        SM = pyopenms.InstrumentSettings.ScanMode
        # Test known scan mode strings
        self.assertEqual(
            pyopenms.InstrumentSettings.toScanMode("MS1Spectrum"),
            SM.MS1SPECTRUM
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.toScanMode("MSnSpectrum"),
            SM.MSNSPECTRUM
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.toScanMode("SelectedIonMonitoring"),
            SM.SIM
        )
        self.assertEqual(
            pyopenms.InstrumentSettings.toScanMode("SelectedReactionMonitoring"),
            SM.SRM
        )

    def test_roundtrip_scanMode(self):
        """Test that enum -> string -> enum roundtrip works."""
        SM = pyopenms.InstrumentSettings.ScanMode
        modes = [
            SM.UNKNOWN,
            SM.MASSSPECTRUM,
            SM.MS1SPECTRUM,
            SM.MSNSPECTRUM,
            SM.SIM,
            SM.SRM,
            SM.CRM,
            SM.PRECURSOR,
        ]
        for mode in modes:
            mode_str = pyopenms.InstrumentSettings.scanModeToString(mode)
            result = pyopenms.InstrumentSettings.toScanMode(mode_str)
            self.assertEqual(result, mode)

    def test_getAllNamesOfScanMode(self):
        """Test getting all scan mode names."""
        names = pyopenms.InstrumentSettings.getAllNamesOfScanMode()
        self.assertGreater(len(names), 0)
        # getAllNames returns str
        self.assertIn("MS1Spectrum", names)
        self.assertIn("MSnSpectrum", names)

    def test_toScanMode_invalid(self):
        """Test that invalid string raises exception."""
        with self.assertRaises(Exception):
            pyopenms.InstrumentSettings.toScanMode("InvalidScanMode")


if __name__ == '__main__':
    unittest.main()
