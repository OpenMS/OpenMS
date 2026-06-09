"""Tests for FeatureFinderAlgorithmMetaboIdent DataFrame helper methods.

Covers the nanobind addon `compounds_from_df` / `run_from_df`, including the
optional IonMobility and Adduct columns.
"""

import unittest

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
import pyopenms


def _base_df(**overrides):
    """A minimal single-compound DataFrame; override individual columns."""
    data = {
        'CompoundName': ['glucose'],
        'SumFormula': ['C6H12O6'],
        'Mass': [0.0],
        'Charge': [-1],
        'RetentionTime': [123.4],
        'RetentionTimeRange': [0.0],
        'IsoDistribution': [0.0],
    }
    data.update(overrides)
    return pd.DataFrame(data)


class TestCompoundsFromDF(unittest.TestCase):
    """Test cases for FeatureFinderAlgorithmMetaboIdent.compounds_from_df()."""

    def test_basic_conversion(self):
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
        self.assertEqual(compounds[0].getName(), 'glucose')
        self.assertEqual(compounds[0].getFormula(), 'C6H12O6')
        self.assertAlmostEqual(compounds[0].getMass(), 0.0)
        self.assertEqual(list(compounds[0].getCharges()), [-1])
        self.assertEqual(list(compounds[0].getRTs()), [123.4])
        self.assertEqual(list(compounds[0].getRTRanges()), [0.0])
        self.assertEqual(list(compounds[0].getIsotopeDistribution()), [0.0])
        # No IM / adduct columns -> defaults
        self.assertEqual(list(compounds[0].getIonMobilities()), [])
        self.assertEqual(compounds[0].getAdduct(), '')

        self.assertEqual(compounds[1].getName(), 'fructose')
        self.assertAlmostEqual(compounds[1].getMass(), 180.063)
        self.assertEqual(list(compounds[1].getCharges()), [1])
        self.assertEqual(list(compounds[1].getRTs()), [200.0])
        self.assertEqual(list(compounds[1].getRTRanges()), [10.0])

    def test_multivalue_as_python_lists(self):
        df = _base_df(
            CompoundName=['compound1'],
            Charge=[[-1, 1, 2]],
            RetentionTime=[[100.0, 200.0, 300.0]],
            RetentionTimeRange=[[5.0, 10.0, 15.0]],
            IsoDistribution=[[0.5, 0.3, 0.2]],
        )

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(list(compounds[0].getCharges()), [-1, 1, 2])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0, 300.0])
        self.assertEqual(list(compounds[0].getRTRanges()), [5.0, 10.0, 15.0])
        self.assertAlmostEqual(compounds[0].getIsotopeDistribution()[0], 0.5)

    def test_multivalue_as_comma_separated_strings(self):
        df = _base_df(
            CompoundName=['compound1'],
            Charge=['-1, 1, 2'],
            RetentionTime=['100.0, 200.0, 300.0'],
            RetentionTimeRange=['5.0, 10.0'],
            IsoDistribution=['0.5, 0.3, 0.2'],
        )

        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)

        self.assertEqual(list(compounds[0].getCharges()), [-1, 1, 2])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0, 300.0])
        self.assertEqual(list(compounds[0].getRTRanges()), [5.0, 10.0])

    def test_ion_mobility_column(self):
        """Optional IonMobility column is exposed on the compound."""
        df = _base_df(IonMobility=[0.95])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getIonMobilities()), [0.95])

    def test_ion_mobility_multivalue(self):
        """IonMobility accepts one value per RT entry (list or comma string)."""
        df = _base_df(
            RetentionTime=[[100.0, 200.0]],
            IonMobility=['0.8, 0.9'],
        )
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getIonMobilities()), [0.8, 0.9])

    def test_ion_mobility_alias_and_zero(self):
        """The 'im' alias works and a 0/NaN value disables IM filtering."""
        df = _base_df(im=[np.nan])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getIonMobilities()), [])

    def test_adduct_column(self):
        """Optional Adduct column is exposed on the compound."""
        df = _base_df(Adduct=['[M-H]-'])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getAdduct(), '[M-H]-')

    def test_adduct_nan_becomes_empty(self):
        df = _base_df(Adduct=[np.nan])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getAdduct(), '')

    def test_case_insensitive_columns(self):
        df = pd.DataFrame({
            'compoundname': ['glucose'], 'sumformula': ['C6H12O6'], 'mass': [0.0],
            'charge': [-1], 'retentiontime': [123.4], 'retentiontimerange': [0.0],
            'isodistribution': [0.0]
        })
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getName(), 'glucose')

    def test_snake_case_columns(self):
        df = pd.DataFrame({
            'compound_name': ['glucose'], 'sum_formula': ['C6H12O6'], 'mass': [0.0],
            'charge': [-1], 'retention_time': [123.4], 'retention_time_range': [0.0],
            'iso_distribution': [0.0]
        })
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getName(), 'glucose')

    def test_rt_shorthand_column(self):
        df = pd.DataFrame({
            'CompoundName': ['glucose'], 'SumFormula': ['C6H12O6'], 'Mass': [0.0],
            'Charge': [-1], 'rt': [123.4], 'rt_range': [0.0], 'IsoDistribution': [0.0]
        })
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getRTs()), [123.4])

    def test_duplicate_column_mapping_raises(self):
        """Two columns mapping to the same canonical field is an error."""
        df = _base_df()
        df['rt'] = [999.0]  # both 'RetentionTime' and 'rt' present
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('map', str(ctx.exception))

    def test_missing_required_column(self):
        df = pd.DataFrame({'CompoundName': ['glucose'], 'SumFormula': ['C6H12O6']})
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('Missing required columns', str(ctx.exception))

    def test_empty_compound_name(self):
        df = _base_df(CompoundName=[''])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('CompoundName cannot be empty', str(ctx.exception))

    def test_nan_compound_name(self):
        df = _base_df(CompoundName=[np.nan])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('CompoundName cannot be empty', str(ctx.exception))

    def test_duplicate_compound_names(self):
        df = pd.DataFrame({
            'CompoundName': ['glucose', 'glucose'], 'SumFormula': ['C6H12O6', 'C6H12O6'],
            'Mass': [0.0, 0.0], 'Charge': [-1, 1], 'RetentionTime': [123.4, 200.0],
            'RetentionTimeRange': [0.0, 0.0], 'IsoDistribution': [0.0, 0.0]
        })
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('Duplicate CompoundName', str(ctx.exception))

    def test_empty_charge(self):
        df = _base_df(Charge=[''])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('Charge cannot be empty', str(ctx.exception))

    def test_empty_retention_time(self):
        df = _base_df(RetentionTime=[''])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('RetentionTime cannot be empty', str(ctx.exception))

    def test_unparseable_charge_reports_row(self):
        df = _base_df(Charge=['abc'])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('Charge', str(ctx.exception))

    def test_nan_in_optional_fields(self):
        # Keep a valid SumFormula so the m/z can still be derived (Mass NaN -> 0.0).
        df = _base_df(Mass=[np.nan],
                      RetentionTimeRange=[np.nan], IsoDistribution=[np.nan])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getFormula(), 'C6H12O6')
        self.assertAlmostEqual(compounds[0].getMass(), 0.0)
        self.assertEqual(list(compounds[0].getRTRanges()), [0.0])
        self.assertEqual(list(compounds[0].getIsotopeDistribution()), [0.0])

    def test_formula_nan_with_mass_ok(self):
        # No formula but a positive mass is still valid (m/z derives from mass).
        df = _base_df(SumFormula=[np.nan], Mass=[180.063])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(compounds[0].getFormula(), '')
        self.assertAlmostEqual(compounds[0].getMass(), 180.063)

    def test_no_mass_no_formula_raises(self):
        df = _base_df(SumFormula=[''], Mass=[0.0])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('Mass', str(ctx.exception))

    def test_zero_charge_raises(self):
        df = _base_df(Charge=[0])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('non-zero', str(ctx.exception))

    def test_rt_range_length_mismatch_raises(self):
        # 2 RT ranges for 3 RTs is neither 1 nor len(rts) -> error.
        df = _base_df(RetentionTime=[[100.0, 200.0, 300.0]],
                      RetentionTimeRange=[[5.0, 10.0]])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('RetentionTimeRange', str(ctx.exception))

    def test_ion_mobility_length_mismatch_raises(self):
        df = _base_df(RetentionTime=[[100.0, 200.0, 300.0]],
                      IonMobility=[[0.8, 0.9]])
        with self.assertRaises(ValueError) as ctx:
            pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertIn('IonMobility', str(ctx.exception))

    def test_rt_range_broadcast_ok(self):
        # A single RT range broadcasts to all RTs.
        df = _base_df(RetentionTime=[[100.0, 200.0, 300.0]],
                      RetentionTimeRange=[5.0])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getRTRanges()), [5.0])

    def test_numpy_arrays_in_multivalue_fields(self):
        df = _base_df(
            CompoundName=['compound1'],
            Charge=[np.array([-1, 1])],
            RetentionTime=[np.array([100.0, 200.0])],
            RetentionTimeRange=[np.array([5.0])],
            IsoDistribution=[np.array([0.5, 0.5])],
        )
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(list(compounds[0].getCharges()), [-1, 1])
        self.assertEqual(list(compounds[0].getRTs()), [100.0, 200.0])

    def test_empty_dataframe(self):
        df = pd.DataFrame({c: [] for c in (
            'CompoundName', 'SumFormula', 'Mass', 'Charge',
            'RetentionTime', 'RetentionTimeRange', 'IsoDistribution')})
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        self.assertEqual(len(compounds), 0)


class TestCompoundsToDF(unittest.TestCase):
    """Test cases for FeatureFinderAlgorithmMetaboIdent.compounds_to_df()."""

    def test_columns_and_values(self):
        df = _base_df(IonMobility=[0.95], Adduct=['[M-H]-'])
        compounds = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        out = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_to_df(compounds)
        self.assertEqual(list(out.columns), [
            'CompoundName', 'SumFormula', 'Mass', 'Charge', 'RetentionTime',
            'RetentionTimeRange', 'IsoDistribution', 'IonMobility', 'Adduct'])
        self.assertEqual(out.iloc[0]['CompoundName'], 'glucose')
        self.assertEqual(out.iloc[0]['Charge'], [-1])
        self.assertEqual(out.iloc[0]['IonMobility'], [0.95])
        self.assertEqual(out.iloc[0]['Adduct'], '[M-H]-')

    def test_empty_list(self):
        out = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_to_df([])
        self.assertEqual(len(out), 0)
        self.assertIn('IonMobility', out.columns)

    def test_round_trip(self):
        df = pd.DataFrame({
            'CompoundName': ['glucose', 'fructose'],
            'SumFormula': ['C6H12O6', ''],
            'Mass': [0.0, 180.063],
            'Charge': [[-1], [1, 2]],
            'RetentionTime': [[123.4], [200.0, 250.0]],
            'RetentionTimeRange': [[0.0], [5.0]],
            'IsoDistribution': [[0.0], [0.0]],
            'IonMobility': [[], [0.9]],
            'Adduct': ['[M-H]-', ''],
        })
        c1 = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
        df2 = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_to_df(c1)
        c2 = pyopenms.FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df2)

        self.assertEqual(len(c1), len(c2))
        for a, b in zip(c1, c2):
            self.assertEqual(a.getName(), b.getName())
            self.assertEqual(a.getFormula(), b.getFormula())
            self.assertAlmostEqual(a.getMass(), b.getMass())
            self.assertEqual(list(a.getCharges()), list(b.getCharges()))
            self.assertEqual(list(a.getRTs()), list(b.getRTs()))
            self.assertEqual(list(a.getRTRanges()), list(b.getRTRanges()))
            self.assertEqual(list(a.getIonMobilities()), list(b.getIonMobilities()))
            self.assertEqual(a.getAdduct(), b.getAdduct())


class TestRunFromDF(unittest.TestCase):
    """Test cases for FeatureFinderAlgorithmMetaboIdent.run_from_df()."""

    def test_run_from_df_exists(self):
        ff = pyopenms.FeatureFinderAlgorithmMetaboIdent()
        self.assertTrue(callable(getattr(ff, 'run_from_df', None)))


if __name__ == '__main__':
    unittest.main()
