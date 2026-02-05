import unittest

import pyopenms


class TestDataProcessing(unittest.TestCase):

    def test_processingActionToString(self):
        """Test converting ProcessingAction enum to string."""
        # Test known processing actions
        self.assertEqual(
            pyopenms.DataProcessing.processingActionToString(
                pyopenms.DataProcessing.ProcessingAction.PEAK_PICKING),
            "Peak picking"
        )
        self.assertEqual(
            pyopenms.DataProcessing.processingActionToString(
                pyopenms.DataProcessing.ProcessingAction.SMOOTHING),
            "Smoothing"
        )
        self.assertEqual(
            pyopenms.DataProcessing.processingActionToString(
                pyopenms.DataProcessing.ProcessingAction.CALIBRATION),
            "Calibration of m/z positions"
        )
        self.assertEqual(
            pyopenms.DataProcessing.processingActionToString(
                pyopenms.DataProcessing.ProcessingAction.DEISOTOPING),
            "Deisotoping"
        )

    def test_toProcessingAction(self):
        """Test converting string to ProcessingAction enum."""
        # Test known processing action strings
        self.assertEqual(
            pyopenms.DataProcessing.toProcessingAction("Peak picking"),
            pyopenms.DataProcessing.ProcessingAction.PEAK_PICKING
        )
        self.assertEqual(
            pyopenms.DataProcessing.toProcessingAction("Smoothing"),
            pyopenms.DataProcessing.ProcessingAction.SMOOTHING
        )
        self.assertEqual(
            pyopenms.DataProcessing.toProcessingAction("Calibration of m/z positions"),
            pyopenms.DataProcessing.ProcessingAction.CALIBRATION
        )
        self.assertEqual(
            pyopenms.DataProcessing.toProcessingAction("Deisotoping"),
            pyopenms.DataProcessing.ProcessingAction.DEISOTOPING
        )

    def test_roundtrip_processingAction(self):
        """Test that enum -> string -> enum roundtrip works."""
        PA = pyopenms.DataProcessing.ProcessingAction
        actions = [
            PA.DATA_PROCESSING,
            PA.CHARGE_DECONVOLUTION,
            PA.DEISOTOPING,
            PA.SMOOTHING,
            PA.PEAK_PICKING,
            PA.ALIGNMENT,
            PA.CALIBRATION,
            PA.NORMALIZATION,
            PA.FILTERING,
        ]
        for action in actions:
            action_str = pyopenms.DataProcessing.processingActionToString(action)
            result = pyopenms.DataProcessing.toProcessingAction(action_str)
            self.assertEqual(result, action)

    def test_getAllNamesOfProcessingAction(self):
        """Test getting all processing action names."""
        names = pyopenms.DataProcessing.getAllNamesOfProcessingAction()
        self.assertGreater(len(names), 0)
        # getAllNames returns str
        self.assertIn("Peak picking", names)
        self.assertIn("Smoothing", names)
        self.assertIn("Calibration of m/z positions", names)

    def test_toProcessingAction_invalid(self):
        """Test that invalid string raises exception."""
        with self.assertRaises(Exception):
            pyopenms.DataProcessing.toProcessingAction("InvalidProcessingAction")


if __name__ == '__main__':
    unittest.main()
