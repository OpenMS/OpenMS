"""Integrity tests for the V16 descriptive threshold audit."""
import ast
import unittest

import numpy as np
import pandas as pd

from threshold_audit import V, n_at, tail_percentile, tdc


class AuditTests(unittest.TestCase):
    def population(self):
        return pd.DataFrame({'sample_id': np.arange(10), 'precursor_id': [0,0,1,2,3,4,5,6,7,8],
                             'pair_key': ['a','a','a','b','b','c','c','d','d','e'],
                             'is_decoy': [0,0,1,0,1,0,1,0,1,1], 'fold': [0]*10})

    def test_original_evaluator_equivalence(self):
        tree = ast.parse((V/'evaluate_v2.py').read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'tdc')
        namespace = {'np': np, 'pd': pd}
        exec(compile(ast.Module(body=[node], type_ignores=[]), 'evaluate_v2.py', 'exec'), namespace)
        pop = self.population()
        scores = np.array([3,8,8,7,6,4,5,9,1,0], dtype=float)
        w, b = tdc(pop, scores)
        ref_w, ref_b = namespace['tdc'](pop, scores)
        pd.testing.assert_frame_equal(w.drop(columns='row'), ref_w)
        pd.testing.assert_frame_equal(b.drop(columns='row'), ref_b)
        self.assertTrue(b.loc[b.pair_key == 'a', 'is_decoy'].item() == 1)
        self.assertTrue(np.all(np.diff(b.q) >= 0))

    def test_sample_id_breaks_window_ties(self):
        pop = self.population()
        w, _ = tdc(pop, np.array([8,8,7,6,5,4,3,2,1,0], dtype=float))
        self.assertEqual(w.loc[w.precursor_id == 0, 'sample_id'].item(), 0)

    def test_tied_null_percentile(self):
        np.testing.assert_allclose(tail_percentile([0,1,2], [1,1]), [0,.5,1])

    def test_fdp_no_valid_prefix(self):
        self.assertEqual(n_at(np.array([3.,2.,1.]), np.ones(3), .01), 0)

    def test_fdp_all_human(self):
        self.assertEqual(n_at(np.array([3.,2.,1.]), np.zeros(3), .01), 3)

    def test_reject_degenerate_normalization(self):
        with self.assertRaises(ValueError):
            tdc(self.population(), np.ones(10))


if __name__ == '__main__':
    unittest.main()
