import unittest, numpy as np
from peak_core import local_peak, support

class T(unittest.TestCase):
    def test_gaussian_fwhm(self):
        rt = np.arange(0, 30, 1.65); s = 3.3 / 2.3548
        y = np.exp(-0.5 * ((rt - 15) / s) ** 2)
        a, w, p = local_peak(rt, y, 15.2)
        self.assertEqual(rt[a], rt[np.argmax(y)]); self.assertAlmostEqual(w, 3.3, delta=0.6); self.assertGreaterEqual(p, 1)
    def test_apex_search_limited(self):
        rt = np.arange(10.); y = np.array([0, 0, 1, 2, 1, 0, 0, 0, 9, 0.])
        a, w, p = local_peak(rt, y, 3.0)
        self.assertEqual(a, 3)               # the far stronger peak at 8 is outside +-2 scans
        self.assertAlmostEqual(w, 2.0)       # half max 1: crossings at 2 and 4 exactly
        self.assertEqual(p, 1)
    def test_edge_and_zero(self):
        rt = np.arange(5.); self.assertEqual(local_peak(rt, np.zeros(5), 2)[1:], (0.0, 0))
        a, w, p = local_peak(rt, np.array([5, 4, 3, 1, 0.]), 0)
        self.assertEqual(a, 0); self.assertAlmostEqual(w, 2.0 + (3 - 2.5) / 2)  # left edge at rt 0, right crossing at 2.25
    def test_support(self):
        M = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1], [0, 0, 0]], bool)
        self.assertEqual(support(M, 1), 3); self.assertEqual(support(M, 3), 0)

if __name__ == '__main__':
    unittest.main(verbosity=2)
