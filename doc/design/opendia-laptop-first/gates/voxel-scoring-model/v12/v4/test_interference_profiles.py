import unittest
import numpy as np

from interference_profiles import correct_batch


class InterferenceProfilesTest(unittest.TestCase):
    def setUp(self):
        self.t = np.arange(32)
        self.peak = np.exp(-0.5 * ((self.t - 16) / 2.2) ** 2).astype(np.float32)
        self.lib = np.ones((3, 6), dtype=np.float32)

    def test_caps_off_apex_interference_bump(self):
        x = np.tile(self.peak, (3, 6, 1))
        x[0, 3] += 0.8 * np.exp(-0.5 * ((self.t - 24) / 1.2) ** 2)
        corrected, features = correct_batch(x, self.lib)
        self.assertLess(corrected[0, 3, 24], x[0, 3, 24])
        self.assertTrue(np.isfinite(features).all())

    def test_single_active_trace_is_preserved(self):
        x = np.zeros((1, 6, 32), dtype=np.float32)
        x[0, 2] = self.peak
        corrected, features = correct_batch(x, self.lib[:1])
        self.assertGreater(float(corrected[0, 2].sum()), 0.0)
        self.assertTrue(np.isfinite(features).all())

    def test_empty_window_stays_empty(self):
        x = np.zeros((1, 6, 32), dtype=np.float32)
        corrected, features = correct_batch(x, self.lib[:1])
        self.assertFalse(np.any(corrected))
        self.assertTrue(np.isfinite(features).all())


if __name__ == '__main__':
    unittest.main()
