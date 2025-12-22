import unittest
import numpy as np
import pyopenms
import math


class TestMultipleTestingComputeModelFDR(unittest.TestCase):
    """Tests for compute_model_fdr function."""

    def test_basic_float64(self):
        """Test basic FDR computation with float64."""
        values = [0.1, 0.2, 0.3]
        result = pyopenms.compute_model_fdr(values, dtype="float64")
        self.assertEqual(len(result), 3)
        self.assertTrue(all(isinstance(x, (float, np.floating)) for x in result))
        # FDR should be monotone non-decreasing
        for i in range(1, len(result)):
            self.assertGreaterEqual(result[i], result[i-1])

    def test_dtype_float32(self):
        """Test with float32 dtype."""
        values = [0.05, 0.10, 0.15]
        result = pyopenms.compute_model_fdr(values, dtype="float32")
        self.assertEqual(len(result), 3)
        self.assertTrue(all(x >= 0 for x in result))

    def test_dtype_int32(self):
        """Test with int32 dtype."""
        values = [1, 2, 3, 4, 5]
        result = pyopenms.compute_model_fdr(values, dtype="int32")
        self.assertEqual(len(result), 5)
        self.assertTrue(all(x >= 0 for x in result))

    def test_invalid_dtype(self):
        """Test that invalid dtype raises ValueError."""
        with self.assertRaises(ValueError):
            pyopenms.compute_model_fdr([0.1], dtype="invalid")

    def test_single_value(self):
        """Test with single value."""
        result = pyopenms.compute_model_fdr([0.5], dtype="float64")
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0] >= 0)


class TestMultipleTesting_pemp(unittest.TestCase):
    """Tests for pemp (empirical p-value) function."""

    def test_basic_empirical_pvalue(self):
        """Test basic empirical p-value computation."""
        stat = [0.5, 1.0, 1.5, 2.0]
        stat0 = [0.0, 0.5, 1.0]
        result = pyopenms.pemp(stat, stat0, dtype="float64")
        
        self.assertEqual(len(result), len(stat))
        self.assertTrue(all(0 <= p <= 1 for p in result))

    def test_pemp_float32(self):
        """Test pemp with float32 dtype."""
        stat = [1.0, 2.0, 3.0]
        stat0 = [0.5, 1.5, 2.5]
        result = pyopenms.pemp(stat, stat0, dtype="float32")
        self.assertEqual(len(result), len(stat))

    def test_pemp_int32(self):
        """Test pemp with int32 dtype."""
        stat = [2, 3, 4]
        stat0 = [1, 2, 3]
        result = pyopenms.pemp(stat, stat0, dtype="int32")
        self.assertEqual(len(result), len(stat))

    def test_pemp_invalid_dtype(self):
        """Test that invalid dtype raises ValueError."""
        with self.assertRaises(ValueError):
            pyopenms.pemp([1.0], [0.5], dtype="bad_dtype")

    def test_pemp_single_points(self):
        """Test with single points."""
        stat = [2.0]
        stat0 = [1.0, 1.5, 2.0, 2.5, 3.0]
        result = pyopenms.pemp(stat, stat0)
        self.assertEqual(len(result), 1)
        self.assertTrue(0 <= result[0] <= 1)


class TestMultipleTesting_qvalue(unittest.TestCase):
    """Tests for qvalue function."""

    def test_basic_qvalue(self):
        """Test basic q-value computation."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        result = pyopenms.qvalue(p_values, pi0=1.0, pfdr=False)
        
        self.assertEqual(len(result), len(p_values))
        # q-values should be in [0,1] and monotone
        for i, q in enumerate(result):
            self.assertTrue(0 <= q <= 1 or math.isnan(q))
        # Monotone non-decreasing in p-value order
        for i in range(1, len(result)):
            self.assertGreaterEqual(result[i], result[i-1])

    def test_qvalue_with_pi0(self):
        """Test q-value with explicit pi0 estimate."""
        p_values = [0.001, 0.01, 0.05, 0.2]
        result = pyopenms.qvalue(p_values, pi0=0.8, pfdr=False)
        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= q <= 1 or math.isnan(q) for q in result))

    def test_qvalue_pfdr_mode(self):
        """Test q-value in positive FDR mode."""
        p_values = [0.001, 0.01, 0.05, 0.1]
        result = pyopenms.qvalue(p_values, pi0=1.0, pfdr=True)
        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= q <= 1 or math.isnan(q) for q in result))



    def test_qvalue_with_nan(self):
        """Test handling of NaN values."""
        p_values = [0.01, np.nan, 0.05]
        result = pyopenms.qvalue(p_values, pi0=1.0)
        self.assertEqual(len(result), len(p_values))
        # NaN in input typically results in NaN in output
        self.assertTrue(math.isnan(result[1]))


class TestMultipleTesting_pnorm(unittest.TestCase):
    """Tests for pnorm (parametric p-value) function."""

    def test_basic_pnorm(self):
        """Test basic parametric p-value computation."""
        stat = [0.0, 1.0, 2.0, -1.0]
        stat0 = [-1.0, -0.5, 0.0, 0.5, 1.0]  # null distribution
        
        result = pyopenms.pnorm(stat, stat0)
        
        self.assertEqual(len(result), len(stat))
        self.assertTrue(all(0 <= p <= 1 or math.isnan(p) for p in result))

    def test_pnorm_standard_normal_null(self):
        """Test with approximately standard normal null."""
        np.random.seed(42)
        stat0 = np.random.randn(100)  # N(0,1)
        stat = [0.0, 1.0, 2.0]  # should have tail probs roughly 0.5, 0.16, 0.023
        
        result = pyopenms.pnorm(stat, stat0)
        
        self.assertEqual(len(result), len(stat))
        # Check rough ordering (p-values should be decreasing for increasing stat)
        self.assertGreater(result[0], result[1])
        self.assertGreater(result[1], result[2])

    def test_pnorm_with_nan(self):
        """Test handling of NaN values."""
        stat = [0.0, np.nan, 1.0]
        stat0 = [0.0, 1.0]
        result = pyopenms.pnorm(stat, stat0)
        self.assertTrue(math.isnan(result[1]))


class TestMultipleTesting_pi0est(unittest.TestCase):
    """Tests for pi0est (pi0 estimation) function."""

    def test_basic_pi0est(self):
        """Test basic pi0 estimation."""
        np.random.seed(42)
        # Mix of true positives (uniform near 0) and nulls (uniform)
        true_pos = np.random.uniform(0, 0.01, 50)
        nulls = np.random.uniform(0, 1, 450)
        p_values = np.concatenate([true_pos, nulls])
        
        result = pyopenms.pi0est(p_values)
        
        self.assertIn('pi0', result)
        self.assertIn('pi0_lambda', result)
        self.assertIn('lambda', result)
        self.assertIn('pi0_smooth', result)
        
        # pi0 should be positive and <= 1
        self.assertTrue(0 < result['pi0'] <= 1.0)
        # Should be roughly 0.9 (450/500)
        self.assertTrue(0.75 < result['pi0'] <= 1.0)

    def test_pi0est_with_lambda(self):
        """Test pi0 estimation with explicit lambda values."""
        p_values = np.concatenate([
            np.random.uniform(0, 0.01, 30),
            np.random.uniform(0, 1, 70)
        ])
        lambda_vals = [0.3, 0.5, 0.7]
        
        result = pyopenms.pi0est(p_values, lambda_=lambda_vals)
        
        self.assertTrue(0 < result['pi0'] <= 1.0)
        self.assertEqual(len(result['pi0_lambda']), len(lambda_vals))

    def test_pi0est_default_lambda(self):
        """Test that default lambda is generated correctly."""
        p_values = np.random.uniform(0, 1, 100)
        
        result = pyopenms.pi0est(p_values)
        
        # Default lambda should be ~0.05, 0.10, ..., 0.95
        self.assertTrue(len(result['lambda']) > 0)
        self.assertTrue(result['lambda'][0] > 0)
        self.assertTrue(result['lambda'][-1] < 1.0)

    def test_pi0est_return_dict_structure(self):
        """Test that return dict has correct structure."""
        p_values = [0.01, 0.02, 0.05, 0.1, 0.5, 0.9]
        result = pyopenms.pi0est(p_values)
        
        required_keys = {'pi0', 'pi0_lambda', 'lambda', 'pi0_smooth'}
        self.assertEqual(set(result.keys()), required_keys)

    def test_pi0est_all_nulls(self):
        """Test with all uniform p-values (all nulls)."""
        p_values = np.random.uniform(0, 1, 100)
        result = pyopenms.pi0est(p_values)
        
        # Should estimate pi0 close to 1.0
        self.assertTrue(0.8 < result['pi0'] <= 1.0)


class TestMultipleTesting_lfdr(unittest.TestCase):
    """Tests for lfdr (local false discovery rate) function."""

    def test_basic_lfdr_probit(self):
        """Test basic local FDR with probit transformation."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9]
        pi0 = 0.7
        
        result = pyopenms.lfdr(p_values, pi0, trunc=True, monotone=True, transf="probit")
        
        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= lfdr <= 1 for lfdr in result))

    def test_lfdr_logit_transformation(self):
        """Test local FDR with logit transformation."""
        p_values = [0.01, 0.05, 0.1, 0.5]
        pi0 = 0.8
        
        result = pyopenms.lfdr(
            p_values, pi0, 
            trunc=True, monotone=True, 
            transf="logit"
        )
        
        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= x <= 1 for x in result))

    def test_lfdr_without_monotone(self):
        """Test local FDR without monotonicity constraint."""
        p_values = [0.01, 0.05, 0.1, 0.5, 0.9]
        pi0 = 0.8
        
        result = pyopenms.lfdr(
            p_values, pi0,
            trunc=True, monotone=False,
            transf="probit"
        )
        
        self.assertEqual(len(result), len(p_values))
        # Values should still be in [0,1] due to truncation
        self.assertTrue(all(0 <= x <= 1 for x in result))

    def test_lfdr_without_trunc(self):
        """Test local FDR without truncation."""
        p_values = [0.001, 0.01, 0.05]
        pi0 = 0.5
        
        result = pyopenms.lfdr(
            p_values, pi0,
            trunc=False, monotone=True,
            transf="probit"
        )
        
        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(math.isfinite(x) for x in result))

    def test_lfdr_custom_bandwidth_adj(self):
        """Test with custom bandwidth adjustment."""
        p_values = [0.01, 0.05, 0.1, 0.5]
        pi0 = 0.7
        
        result = pyopenms.lfdr(
            p_values, pi0,
            adj=2.0,  # wider bandwidth
            transf="probit"
        )
        
        self.assertEqual(len(result), len(p_values))

    def test_lfdr_custom_eps(self):
        """Test with custom eps parameter."""
        p_values = [0.001, 0.01, 0.1]
        pi0 = 0.6
        
        result = pyopenms.lfdr(
            p_values, pi0,
            eps=1e-6,
            transf="probit"
        )
        
        self.assertEqual(len(result), len(p_values))

    def test_lfdr_monotonicity(self):
        """Test that lfdr is monotone in p-value when monotone=True."""
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        pi0 = 0.7
        
        result = pyopenms.lfdr(p_values, pi0, trunc=True, monotone=True, transf="probit")
        
        # With monotone=True, lfdr should be non-decreasing in p-value
        for i in range(1, len(result)):
            self.assertGreaterEqual(result[i], result[i-1] - 1e-6)  # small tolerance for rounding


class TestMultipleTestingIntegration(unittest.TestCase):
    """Integration tests combining multiple functions."""

    def test_workflow_pi0_then_qvalue(self):
        """Test typical workflow: estimate pi0, then compute q-values."""
        np.random.seed(42)
        # Generate mixed p-values
        true_pos = np.random.uniform(0, 0.01, 50)
        nulls = np.random.uniform(0, 1, 450)
        p_values = np.concatenate([true_pos, nulls])
        
        # Step 1: Estimate pi0
        pi0_result = pyopenms.pi0est(p_values)
        pi0 = pi0_result['pi0']
        
        # Step 2: Compute q-values
        q_values = pyopenms.qvalue(p_values, pi0=pi0, pfdr=False)
        
        self.assertEqual(len(q_values), len(p_values))
        self.assertTrue(all(0 <= q <= 1 or math.isnan(q) for q in q_values))

    def test_workflow_pi0_then_lfdr(self):
        """Test workflow: estimate pi0, then compute local FDR."""
        np.random.seed(42)
        p_values = np.concatenate([
            np.random.uniform(0, 0.05, 30),
            np.random.uniform(0, 1, 170)
        ])
        
        # Estimate pi0
        pi0_result = pyopenms.pi0est(p_values)
        
        # Compute local FDR
        lfdr_result = pyopenms.lfdr(
            p_values, 
            pi0=pi0_result['pi0'],
            transf="probit"
        )
        
        self.assertEqual(len(lfdr_result), len(p_values))

    def test_array_conversion(self):
        """Test that functions handle different array types."""
        # Test with list
        result1 = pyopenms.qvalue([0.01, 0.05, 0.1])
        # Test with numpy array
        result2 = pyopenms.qvalue(np.array([0.01, 0.05, 0.1]))
        
        np.testing.assert_allclose(result1, result2)

    def test_large_input(self):
        """Test with larger input sizes."""
        np.random.seed(42)
        p_values = np.random.uniform(0, 1, 1000)
        
        # Should handle without issues
        pi0_result = pyopenms.pi0est(p_values)
        self.assertTrue(0 < pi0_result['pi0'] <= 1.0)
        
        q_values = pyopenms.qvalue(p_values)
        self.assertEqual(len(q_values), 1000)


if __name__ == '__main__':
    unittest.main()
