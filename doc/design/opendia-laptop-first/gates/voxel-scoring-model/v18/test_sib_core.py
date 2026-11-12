import unittest, numpy as np
from sib_core import sib_features, fold_z

class T(unittest.TestCase):
    def test_siblings(self):
        prec = np.array([1, 1, 2, 2, 3, 4])
        pep = np.array([10, 10, 10, 10, 20, 30])
        rt = np.array([100., 200., 103., 150., 100., 100.])
        z = np.array([5., 9., 2., 30., 7., 1.])
        f = sib_features(prec, pep, rt, z)
        # row0 (prec1 @100): sibling prec2 @103 within 3 s (boundary inclusive) -> z 2; own other window @200 excluded
        np.testing.assert_allclose(f[0], [1, 1, 2])
        # row1 (prec1 @200): no sibling within 3 s
        np.testing.assert_allclose(f[1], [1, 0, 0])
        # row2 (prec2 @103): sibling prec1 @100 -> 5; row3 @150 nothing near
        np.testing.assert_allclose(f[2], [1, 1, 5]); np.testing.assert_allclose(f[3], [1, 0, 0])
        # single-precursor peptides
        np.testing.assert_allclose(f[4], [0, 0, 0]); np.testing.assert_allclose(f[5], [0, 0, 0])
    def test_clip_and_max(self):
        f = sib_features(np.array([1, 2, 2]), np.array([5, 5, 5]), np.array([0., 1., 2.]), np.array([0., 50., 3.]))
        self.assertEqual(f[0, 2], 20.0)   # max of 50 and 3, clipped at 20
    def test_label_symmetry(self):
        # decoy flag is not an input: identical geometry gives identical features
        a = sib_features(np.array([1, 2]), np.array([7, 7]), np.array([0., 1.]), np.array([1., 2.]))
        b = sib_features(np.array([3, 4]), np.array([8, 8]), np.array([0., 1.]), np.array([1., 2.]))
        np.testing.assert_allclose(a, b)
    def test_fold_z(self):
        prec = np.array([1, 1, 2, 3, 4, 5]); sid = np.arange(6); dec = np.array([1, 1, 1, 1, 0, 1]); fold = np.array([0, 0, 0, 0, 0, 1])
        full = np.array([1., 3., 2., 4., 10., 0.]); w = np.array([1., -1e30, 2., 4., 10., 0.])
        with self.assertRaises(ValueError): fold_z(full, w, prec, sid, dec, fold)  # fold 1 has a single decoy -> MAD 0
        fold = np.zeros(6, int); w2 = w.copy(); w2[5] = 6.
        z = fold_z(full, w2, prec, sid, dec, fold)
        # decoy winners (W choice): prec1 row0 (1, not the higher full 3), prec2 2, prec3 4, prec5 0 -> median 1.5
        d = np.array([1., 2., 4., 0.]); med = np.median(d); mad = 1.4826 * np.median(np.abs(d - med))
        np.testing.assert_allclose(z, (full - med) / mad)

if __name__ == '__main__':
    unittest.main(verbosity=2)
