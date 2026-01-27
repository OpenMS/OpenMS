"""Tests for FeatureFinderAlgorithmMetaboIdent DataFrame helper methods."""

import unittest

import numpy as np
import pandas as pd
import pyopenms


class TestCompoundsFromDF(unittest.TestCase):
    """Test cases for FeatureFinderAlgorithmMetaboIdent.compounds_from_df()."""

    def test_basic_conversion(self):
        """Test basic DataFrame to compound list conversion."""
        df = pd.DataFrame({
            'CompoundName': ['glucose', 'fructose'],
            'SumFormula': ['C6H12O6', 'C6H12O6'],
            'Mass': [0.0, 180.063],
            'Charge': [-1, 1],
            'RetentionTime': [123.4, 200.0],
            'RetentionTimeRange': [0.0, 10.0],
            'IsoDistribution': [0.0, 0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 2)

        # Check first compound
        self.assertEqual(compounds[0].getName(), 'glucose')
        self.assertEqual(compounds[0].getFormula(), 'C6H12O6')
        self.assertAlmostEqual(compounds[0].getMass(), 0.0)
        self.assertEqual(list(compounds[0].getCharges()), [-1])
        self.assertEqual(list(compounds[0].getRTs()), [123.4])
        self.assertEqual(list(compounds[0].getRTRanges()), [0.0])
        self.assertEqual(list(compounds[0].getIsotopeDistribution()), [0.0])

        # Check second compound
        self.assertEqual(compounds[1].getName(), 'fructose')
        self.assertAlmostEqual(compounds[1].getMass(), 180.063)
        self.assertEqual(list(compounds[1].getCharges()), [1])
        self.assertEqual(list(compounds[1].getRTs()), [200.0])
        self.assertEqual(list(compounds[1].getRTRanges()), [10.0])

    def test_multivalue_as_python_lists(self):
        """Test multi-value fields passed as Python lists."""
        df = pd.DataFrame({
            'CompoundName': ['compound1'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [[-1, 1, 2]],
            'RetentionTime': [[100.0, 200.0, 300.0]],
            'RetentionTimeRange': [[5.0, 10.0, 15.0]],
            'IsoDistribution': [[0.5, 0.3, 0.2]]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(list(compounds[0].getCharges()), [-1, 1, 2])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0, 300.0])
        self.assertEqual(list(compounds[0].getRTRanges()), [5.0, 10.0, 15.0])
        self.assertAlmostEqual(compounds[0].getIsotopeDistribution()[0], 0.5)
        self.assertAlmostEqual(compounds[0].getIsotopeDistribution()[1], 0.3)
        self.assertAlmostEqual(compounds[0].getIsotopeDistribution()[2], 0.2)

    def test_multivalue_as_comma_separated_strings(self):
        """Test multi-value fields passed as comma-separated strings."""
        df = pd.DataFrame({
            'CompoundName': ['compound1'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': ['-1, 1, 2'],
            'RetentionTime': ['100.0, 200.0, 300.0'],
            'RetentionTimeRange': ['5.0, 10.0'],
            'IsoDistribution': ['0.5, 0.3, 0.2']
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(list(compounds[0].getCharges()), [-1, 1, 2])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0, 300.0])
        self.assertEqual(list(compounds[0].getRTRanges()), [5.0, 10.0])

    def test_case_insensitive_columns(self):
        """Test that column names are case-insensitive."""
        df = pd.DataFrame({
            'compoundname': ['glucose'],
            'sumformula': ['C6H12O6'],
            'mass': [0.0],
            'charge': [-1],
            'retentiontime': [123.4],
            'retentiontimerange': [0.0],
            'isodistribution': [0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(compounds[0].getName(), 'glucose')

    def test_snake_case_columns(self):
        """Test that snake_case column names are accepted."""
        df = pd.DataFrame({
            'compound_name': ['glucose'],
            'sum_formula': ['C6H12O6'],
            'mass': [0.0],
            'charge': [-1],
            'retention_time': [123.4],
            'retention_time_range': [0.0],
            'iso_distribution': [0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(compounds[0].getName(), 'glucose')

    def test_mixed_case_columns(self):
        """Test mixed case column names."""
        df = pd.DataFrame({
            'COMPOUNDNAME': ['glucose'],
            'SumFormula': ['C6H12O6'],
            'MASS': [0.0],
            'Charge': [-1],
            'RetentionTime': [123.4],
            'retentiontimerange': [0.0],
            'ISO_DISTRIBUTION': [0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(compounds[0].getName(), 'glucose')

    def test_missing_required_column(self):
        """Test that missing required columns raise ValueError."""
        df = pd.DataFrame({
            'CompoundName': ['glucose'],
            'SumFormula': ['C6H12O6'],
            # Missing Mass, Charge, RetentionTime, etc.
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('Missing required columns', str(ctx.exception))

    def test_empty_compound_name(self):
        """Test that empty compound names raise ValueError."""
        df = pd.DataFrame({
            'CompoundName': [''],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [-1],
            'RetentionTime': [123.4],
            'RetentionTimeRange': [0.0],
            'IsoDistribution': [0.0]
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('CompoundName cannot be empty', str(ctx.exception))

    def test_nan_compound_name(self):
        """Test that NaN compound names raise ValueError."""
        df = pd.DataFrame({
            'CompoundName': [np.nan],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [-1],
            'RetentionTime': [123.4],
            'RetentionTimeRange': [0.0],
            'IsoDistribution': [0.0]
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('CompoundName cannot be empty', str(ctx.exception))

    def test_duplicate_compound_names(self):
        """Test that duplicate compound names raise ValueError."""
        df = pd.DataFrame({
            'CompoundName': ['glucose', 'glucose'],
            'SumFormula': ['C6H12O6', 'C6H12O6'],
            'Mass': [0.0, 0.0],
            'Charge': [-1, 1],
            'RetentionTime': [123.4, 200.0],
            'RetentionTimeRange': [0.0, 0.0],
            'IsoDistribution': [0.0, 0.0]
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('Duplicate CompoundName', str(ctx.exception))

    def test_empty_charge(self):
        """Test that empty charge raises ValueError."""
        df = pd.DataFrame({
            'CompoundName': ['glucose'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [''],
            'RetentionTime': [123.4],
            'RetentionTimeRange': [0.0],
            'IsoDistribution': [0.0]
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('Charge cannot be empty', str(ctx.exception))

    def test_empty_retention_time(self):
        """Test that empty retention time raises ValueError."""
        df = pd.DataFrame({
            'CompoundName': ['glucose'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [-1],
            'RetentionTime': [''],
            'RetentionTimeRange': [0.0],
            'IsoDistribution': [0.0]
        })

        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertIn('RetentionTime cannot be empty', str(ctx.exception))

    def test_nan_in_optional_fields(self):
        """Test that NaN values in optional fields default to [0.0]."""
        df = pd.DataFrame({
            'CompoundName': ['glucose'],
            'SumFormula': [np.nan],  # Optional - should become empty string
            'Mass': [np.nan],  # Should become 0.0
            'Charge': [-1],
            'RetentionTime': [123.4],
            'RetentionTimeRange': [np.nan],  # Should become [0.0]
            'IsoDistribution': [np.nan]  # Should become [0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(compounds[0].getFormula(), '')
        self.assertAlmostEqual(compounds[0].getMass(), 0.0)
        self.assertEqual(list(compounds[0].getRTRanges()), [0.0])
        self.assertEqual(list(compounds[0].getIsotopeDistribution()), [0.0])

    def test_numpy_arrays_in_multivalue_fields(self):
        """Test that numpy arrays work for multi-value fields."""
        df = pd.DataFrame({
            'CompoundName': ['compound1'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [np.array([-1, 1])],
            'RetentionTime': [np.array([100.0, 200.0])],
            'RetentionTimeRange': [np.array([5.0])],
            'IsoDistribution': [np.array([0.5, 0.5])]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(list(compounds[0].getCharges()), [-1, 1])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0])

    def test_rt_shorthand_column(self):
        """Test 'rt' as shorthand for RetentionTime."""
        df = pd.DataFrame({
            'CompoundName': ['glucose'],
            'SumFormula': ['C6H12O6'],
            'Mass': [0.0],
            'Charge': [-1],
            'rt': [123.4],
            'rt_range': [0.0],
            'IsoDistribution': [0.0]
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 1)
        self.assertEqual(list(compounds[0].getRTs()), [123.4])

    def test_empty_dataframe(self):
        """Test that empty DataFrame returns empty list."""
        df = pd.DataFrame({
            'CompoundName': [],
            'SumFormula': [],
            'Mass': [],
            'Charge': [],
            'RetentionTime': [],
            'RetentionTimeRange': [],
            'IsoDistribution': []
        })

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(len(compounds), 0)


class TestRunFromDF(unittest.TestCase):
    """Test cases for FeatureFinderAlgorithmMetaboIdent.run_from_df()."""

    def test_run_from_df_exists(self):
        """Test that run_from_df method exists and is callable."""
        ff = pyopenms.FeatureFinderAlgorithmMetaboIdent()
        self.assertTrue(hasattr(ff, 'run_from_df'))
        self.assertTrue(callable(getattr(ff, 'run_from_df')))


if __name__ == '__main__':
    unittest.main()
