import unittest
import numpy as np
import pyopenms
import math


class TestPi0Result(unittest.TestCase):
    """Tests for Pi0Result struct."""

    def test_pi0result_exists(self):
        """Test that Pi0Result is accessible."""
        self.assertTrue(hasattr(pyopenms, 'Pi0Result'))

    def test_pi0result_default_values(self):
        """Test Pi0Result default values."""
        r = pyopenms.Pi0Result()
        self.assertEqual(r.pi0, 1.0)
        self.assertEqual(r.pi0_smooth, False)
        self.assertEqual(len(r.pi0_lambda), 0)
        self.assertEqual(len(r.lambda_), 0)

    def test_pi0result_fields_writable(self):
        """Test that Pi0Result fields can be modified."""
        r = pyopenms.Pi0Result()
        r.pi0 = 0.8
        r.pi0_smooth = True
        self.assertEqual(r.pi0, 0.8)
        self.assertEqual(r.pi0_smooth, True)


class TestMultipleTestingEnums(unittest.TestCase):
    """Tests for MultipleTesting enum classes."""

    def test_pi0method_enum_exists(self):
        """Test Pi0Method enum is attached to MultipleTesting."""
        self.assertTrue(hasattr(pyopenms.MultipleTesting, 'Pi0Method'))

    def test_pi0method_smoother(self):
        """Test Pi0Method.Smoother enum value."""
        method = pyopenms.MultipleTesting.Pi0Method.Smoother
        self.assertIsNotNone(method)

    def test_pi0method_bootstrap(self):
        """Test Pi0Method.Bootstrap enum value."""
        method = pyopenms.MultipleTesting.Pi0Method.Bootstrap
        self.assertIsNotNone(method)

    def test_pi0method_values_distinct(self):
        """Test that Pi0Method enum values are distinct."""
        smoother = pyopenms.MultipleTesting.Pi0Method.Smoother
        bootstrap = pyopenms.MultipleTesting.Pi0Method.Bootstrap
        self.assertNotEqual(smoother, bootstrap)

    def test_lfdr_transform_enum_exists(self):
        """Test LfdrTransform enum is attached to MultipleTesting."""
        self.assertTrue(hasattr(pyopenms.MultipleTesting, 'LfdrTransform'))

    def test_lfdr_transform_probit(self):
        """Test LfdrTransform.Probit enum value."""
        transf = pyopenms.MultipleTesting.LfdrTransform.Probit
        self.assertIsNotNone(transf)

    def test_lfdr_transform_logit(self):
        """Test LfdrTransform.Logit enum value."""
        transf = pyopenms.MultipleTesting.LfdrTransform.Logit
        self.assertIsNotNone(transf)

    def test_lfdr_transform_values_distinct(self):
        """Test that LfdrTransform enum values are distinct."""
        probit = pyopenms.MultipleTesting.LfdrTransform.Probit
        logit = pyopenms.MultipleTesting.LfdrTransform.Logit
        self.assertNotEqual(probit, logit)


class TestMultipleTestingEnumConversion(unittest.TestCase):
    """Tests for enum conversion static methods."""

    def test_pi0method_to_string_smoother(self):
        """Test pi0MethodToString for Smoother."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.pi0MethodToString(MT.Pi0Method.Smoother), "smoother")

    def test_pi0method_to_string_bootstrap(self):
        """Test pi0MethodToString for Bootstrap."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.pi0MethodToString(MT.Pi0Method.Bootstrap), "bootstrap")

    def test_string_to_pi0method_smoother(self):
        """Test toPi0Method for 'smoother'."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.toPi0Method("smoother"), MT.Pi0Method.Smoother)

    def test_string_to_pi0method_bootstrap(self):
        """Test toPi0Method for 'bootstrap'."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.toPi0Method("bootstrap"), MT.Pi0Method.Bootstrap)

    def test_string_to_pi0method_case_insensitive(self):
        """Test toPi0Method is case insensitive."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.toPi0Method("SMOOTHER"), MT.Pi0Method.Smoother)
        self.assertEqual(MT.toPi0Method("Bootstrap"), MT.Pi0Method.Bootstrap)

    def test_lfdr_transform_to_string_probit(self):
        """Test lfdrTransformToString for Probit."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.lfdrTransformToString(MT.LfdrTransform.Probit), "probit")

    def test_lfdr_transform_to_string_logit(self):
        """Test lfdrTransformToString for Logit."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.lfdrTransformToString(MT.LfdrTransform.Logit), "logit")

    def test_string_to_lfdr_transform_probit(self):
        """Test toLfdrTransform for 'probit'."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.toLfdrTransform("probit"), MT.LfdrTransform.Probit)

    def test_string_to_lfdr_transform_logit(self):
        """Test toLfdrTransform for 'logit'."""
        MT = pyopenms.MultipleTesting
        self.assertEqual(MT.toLfdrTransform("logit"), MT.LfdrTransform.Logit)


class TestMultipleTestingQValue(unittest.TestCase):
    """Tests for qValue static method."""

    def test_basic_qvalue(self):
        """Test basic q-value computation."""
        MT = pyopenms.MultipleTesting
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        result = MT.qValue(p_values, 1.0, False)

        self.assertEqual(len(result), len(p_values))
        # q-values should be in [0,1] and monotone non-decreasing
        for q in result:
            self.assertTrue(0 <= q <= 1 or math.isnan(q))
        for i in range(1, len(result)):
            self.assertGreaterEqual(result[i], result[i-1])

    def test_qvalue_with_pi0(self):
        """Test q-value with explicit pi0."""
        MT = pyopenms.MultipleTesting
        p_values = [0.001, 0.01, 0.05, 0.2]
        result = MT.qValue(p_values, 0.8, False)
        self.assertEqual(len(result), len(p_values))

    def test_qvalue_pfdr_mode(self):
        """Test q-value in positive FDR mode."""
        MT = pyopenms.MultipleTesting
        p_values = [0.001, 0.01, 0.05, 0.1]
        result = MT.qValue(p_values, 1.0, True)
        self.assertEqual(len(result), len(p_values))


class TestMultipleTestingPNorm(unittest.TestCase):
    """Tests for pNorm static method."""

    def test_basic_pnorm(self):
        """Test basic parametric p-value computation."""
        MT = pyopenms.MultipleTesting
        stat = [0.0, 1.0, 2.0, -1.0]
        stat0 = [-1.0, -0.5, 0.0, 0.5, 1.0]

        result = MT.pNorm(stat, stat0)

        self.assertEqual(len(result), len(stat))
        self.assertTrue(all(0 <= p <= 1 or math.isnan(p) for p in result))

    def test_pnorm_ordering(self):
        """Test that larger stats have smaller p-values."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        stat0 = list(np.random.randn(100))
        stat = [0.0, 1.0, 2.0]

        result = MT.pNorm(stat, stat0)

        # p-values should decrease as stat increases
        self.assertGreater(result[0], result[1])
        self.assertGreater(result[1], result[2])


class TestMultipleTestingPi0Est(unittest.TestCase):
    """Tests for pi0Est static method."""

    def test_basic_pi0est(self):
        """Test basic pi0 estimation."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        # Mix of true positives and nulls
        true_pos = list(np.random.uniform(0, 0.01, 50))
        nulls = list(np.random.uniform(0, 1, 450))
        p_values = true_pos + nulls

        result = MT.pi0Est(p_values, [], MT.Pi0Method.Smoother, 3, False)

        self.assertTrue(0 < result.pi0 <= 1.0)
        self.assertTrue(len(result.lambda_) > 0)

    def test_pi0est_with_lambda(self):
        """Test pi0 estimation with explicit lambda values."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        p_values = list(np.concatenate([
            np.random.uniform(0, 0.01, 30),
            np.random.uniform(0, 1, 70)
        ]))
        lambda_vals = [0.3, 0.5, 0.7]

        result = MT.pi0Est(p_values, lambda_vals, MT.Pi0Method.Smoother, 3, False)

        self.assertTrue(0 < result.pi0 <= 1.0)
        self.assertEqual(len(result.pi0_lambda), len(lambda_vals))

    def test_pi0est_bootstrap(self):
        """Test pi0 estimation with bootstrap method."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        p_values = list(np.random.uniform(0, 1, 100))

        result = MT.pi0Est(p_values, [], MT.Pi0Method.Bootstrap, 3, False)

        self.assertTrue(0 < result.pi0 <= 1.0)

    def test_pi0est_all_nulls(self):
        """Test with all uniform p-values (all nulls)."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        p_values = list(np.random.uniform(0, 1, 100))

        result = MT.pi0Est(p_values, [], MT.Pi0Method.Smoother, 3, False)

        # Should estimate pi0 close to 1.0
        self.assertTrue(0.8 < result.pi0 <= 1.0)


class TestMultipleTestingLfdr(unittest.TestCase):
    """Tests for lfdr static method."""

    def test_basic_lfdr_probit(self):
        """Test basic local FDR with probit transformation."""
        MT = pyopenms.MultipleTesting
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5, 0.9]
        pi0 = 0.7

        result = MT.lfdr(p_values, pi0, True, True,
                         MT.LfdrTransform.Probit, 1.5, 1e-8, 512, 3.0)

        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= lfdr <= 1 for lfdr in result))

    def test_lfdr_logit_transformation(self):
        """Test local FDR with logit transformation."""
        MT = pyopenms.MultipleTesting
        p_values = [0.01, 0.05, 0.1, 0.5]
        pi0 = 0.8

        result = MT.lfdr(p_values, pi0, True, True,
                         MT.LfdrTransform.Logit, 1.5, 1e-8, 512, 3.0)

        self.assertEqual(len(result), len(p_values))
        self.assertTrue(all(0 <= x <= 1 for x in result))

    def test_lfdr_monotonicity(self):
        """Test that lfdr is monotone when monotone=True."""
        MT = pyopenms.MultipleTesting
        p_values = [0.001, 0.01, 0.05, 0.1, 0.5]
        pi0 = 0.7

        result = MT.lfdr(p_values, pi0, True, True,
                         MT.LfdrTransform.Probit, 1.5, 1e-8, 512, 3.0)

        # With monotone=True, lfdr should be non-decreasing in p-value
        for i in range(1, len(result)):
            self.assertGreaterEqual(result[i], result[i-1] - 1e-6)


class TestMultipleTestingIntegration(unittest.TestCase):
    """Integration tests combining multiple functions."""

    def test_workflow_pi0_then_qvalue(self):
        """Test typical workflow: estimate pi0, then compute q-values."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        true_pos = list(np.random.uniform(0, 0.01, 50))
        nulls = list(np.random.uniform(0, 1, 450))
        p_values = true_pos + nulls

        # Step 1: Estimate pi0
        pi0_result = MT.pi0Est(p_values, [], MT.Pi0Method.Smoother, 3, False)
        pi0 = pi0_result.pi0

        # Step 2: Compute q-values
        q_values = MT.qValue(p_values, pi0, False)

        self.assertEqual(len(q_values), len(p_values))
        self.assertTrue(all(0 <= q <= 1 or math.isnan(q) for q in q_values))

    def test_workflow_pi0_then_lfdr(self):
        """Test workflow: estimate pi0, then compute local FDR."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        p_values = list(np.concatenate([
            np.random.uniform(0, 0.05, 30),
            np.random.uniform(0, 1, 170)
        ]))

        # Estimate pi0
        pi0_result = MT.pi0Est(p_values, [], MT.Pi0Method.Smoother, 3, False)

        # Compute local FDR
        lfdr_result = MT.lfdr(
            p_values,
            pi0_result.pi0,
            True, True,
            MT.LfdrTransform.Probit,
            1.5, 1e-8, 512, 3.0
        )

        self.assertEqual(len(lfdr_result), len(p_values))

    def test_large_input(self):
        """Test with larger input sizes."""
        MT = pyopenms.MultipleTesting
        np.random.seed(42)
        p_values = list(np.random.uniform(0, 1, 1000))

        # Should handle without issues
        pi0_result = MT.pi0Est(p_values, [], MT.Pi0Method.Smoother, 3, False)
        self.assertTrue(0 < pi0_result.pi0 <= 1.0)

        q_values = MT.qValue(p_values, pi0_result.pi0, False)
        self.assertEqual(len(q_values), 1000)


if __name__ == '__main__':
    unittest.main()
