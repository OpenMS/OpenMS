"""
Tests for DataFrame column selection feature in MSSpectrum and MSChromatogram.

Tests cover:
- get_df_columns() discovery method
- Column selection via get_df(columns=[...])
- Default behavior (backward compatibility)
- Non-default columns (ion_mobility_unit, chromatogram_type, comment)
- Meta value handling with column selection
- Edge cases (empty spectra, missing data)
"""

import pytest
import numpy as np

import pyopenms


class TestMSSpectrumColumnSelection:
    """Tests for MSSpectrum.get_df() column selection."""

    @pytest.fixture
    def spectrum_with_data(self):
        """Create a spectrum with peaks, precursor, IM data, and meta values."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(123.45)
        spec.setNativeID('scan=100')
        spec.setMetaValue('total_ion_current', 1000.0)
        spec.setMetaValue('base_peak_mz', 500.0)

        # Set peaks
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([10.0, 20.0, 30.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        # Set precursor
        precursor = pyopenms.Precursor()
        precursor.setMZ(500.0)
        precursor.setCharge(2)
        spec.setPrecursors([precursor])

        # Set ion mobility data using FloatDataArray
        fda = pyopenms.FloatDataArray()
        fda.setName('Ion Mobility')
        for val in [1.5, 2.0, 2.5]:
            fda.push_back(val)
        spec.setFloatDataArrays([fda])
        spec.setDriftTime(2.0)  # Set drift time for containsIMData()

        # Set ion annotations
        sda = pyopenms.StringDataArray()
        for ann in ['b2+', 'y3+', 'b4+']:
            sda.push_back(ann)
        sda.setName('IonNames')
        spec.setStringDataArrays([sda])

        return spec

    @pytest.fixture
    def simple_spectrum(self):
        """Create a simple MS1 spectrum without precursor or IM data."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(50.0)
        spec.setNativeID('scan=50')

        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([10.0, 20.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        return spec

    def test_get_df_columns_full_spectrum(self, spectrum_with_data):
        """Test get_df_columns() returns all expected columns for full spectrum."""
        cols = spectrum_with_data.get_df_columns()

        # Should have core columns
        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'rt' in cols
        assert 'ms_level' in cols
        assert 'native_id' in cols

        # Should have IM column (data present)
        assert 'ion_mobility' in cols

        # Should have precursor columns (MS2)
        assert 'precursor_mz' in cols
        assert 'precursor_charge' in cols

        # Should have annotation column
        assert 'ion_annotation' in cols

        # Should have meta values
        assert 'total_ion_current' in cols
        assert 'base_peak_mz' in cols

        # Should NOT have non-default columns
        assert 'ion_mobility_unit' not in cols

    def test_get_df_columns_simple_spectrum(self, simple_spectrum):
        """Test get_df_columns() for simple MS1 spectrum."""
        cols = simple_spectrum.get_df_columns()

        # Should have core columns
        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'rt' in cols

        # Should NOT have precursor columns (MS1)
        assert 'precursor_mz' not in cols
        assert 'precursor_charge' not in cols

        # Should NOT have IM columns (no IM data)
        assert 'ion_mobility' not in cols

    def test_get_df_columns_no_meta_values(self, spectrum_with_data):
        """Test get_df_columns() with export_meta_values=False."""
        cols = spectrum_with_data.get_df_columns(export_meta_values=False)

        assert 'mz' in cols
        assert 'total_ion_current' not in cols
        assert 'base_peak_mz' not in cols

    def test_get_df_default(self, spectrum_with_data):
        """Test get_df() default behavior returns all expected columns."""
        df = spectrum_with_data.get_df()

        # Check core columns present
        assert 'mz' in df.columns
        assert 'intensity' in df.columns
        assert 'rt' in df.columns
        assert 'ms_level' in df.columns
        assert 'native_id' in df.columns

        # Check data values
        assert len(df) == 3  # 3 peaks
        assert df.loc[0, 'mz'] == 100.0
        assert df.loc[0, 'rt'] == 123.45
        assert df.loc[0, 'ms_level'] == 2

        # Check meta values present
        assert 'total_ion_current' in df.columns

    def test_get_df_minimal_columns(self, spectrum_with_data):
        """Test get_df() with minimal column selection."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity'])

        assert list(df.columns) == ['mz', 'intensity']
        assert len(df) == 3
        assert df.loc[0, 'mz'] == 100.0
        assert df.loc[0, 'intensity'] == 10.0

    def test_get_df_custom_columns(self, spectrum_with_data):
        """Test get_df() with custom column selection."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity', 'rt', 'precursor_mz'])

        assert set(df.columns) == {'mz', 'intensity', 'rt', 'precursor_mz'}
        assert df.loc[0, 'precursor_mz'] == 500.0

    def test_get_df_with_ion_mobility(self, spectrum_with_data):
        """Test get_df() including ion mobility column."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity', 'ion_mobility'])

        assert 'ion_mobility' in df.columns
        assert df.loc[0, 'ion_mobility'] == 1.5
        assert df.loc[1, 'ion_mobility'] == 2.0

    def test_get_df_with_ion_mobility_unit(self, spectrum_with_data):
        """Test get_df() with non-default ion_mobility_unit column."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity', 'ion_mobility_unit'])

        assert 'ion_mobility_unit' in df.columns

    def test_get_df_with_ion_annotation(self, spectrum_with_data):
        """Test get_df() including ion annotation column."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity', 'ion_annotation'])

        assert 'ion_annotation' in df.columns
        assert df.loc[0, 'ion_annotation'] == 'b2+'
        assert df.loc[1, 'ion_annotation'] == 'y3+'

    def test_get_df_with_meta_value(self, spectrum_with_data):
        """Test get_df() with specific meta value column."""
        df = spectrum_with_data.get_df(columns=['mz', 'intensity', 'total_ion_current'])

        assert 'total_ion_current' in df.columns
        assert df.loc[0, 'total_ion_current'] == 1000.0

    def test_get_df_all_columns(self, spectrum_with_data):
        """Test get_df() requesting all available columns including non-defaults."""
        # Get all default columns
        cols = spectrum_with_data.get_df_columns()
        # Add non-default columns
        cols.append('ion_mobility_unit')

        df = spectrum_with_data.get_df(columns=cols)

        # Should have all columns
        assert 'mz' in df.columns
        assert 'ion_mobility_unit' in df.columns
        assert 'total_ion_current' in df.columns

    def test_get_df_missing_precursor(self, simple_spectrum):
        """Test get_df() requesting precursor columns when no precursor present."""
        df = simple_spectrum.get_df(columns=['mz', 'intensity', 'precursor_mz'])

        assert 'precursor_mz' in df.columns
        assert np.isnan(df.loc[0, 'precursor_mz'])

    def test_get_df_empty_spectrum(self):
        """Test get_df() with empty spectrum."""
        spec = pyopenms.MSSpectrum()
        df = spec.get_df()

        assert len(df) == 0

    def test_get_df_columns_empty_spectrum(self):
        """Test get_df_columns() with empty spectrum."""
        spec = pyopenms.MSSpectrum()
        cols = spec.get_df_columns()

        # Should still have core columns
        assert 'mz' in cols
        assert 'intensity' in cols


class TestMSChromatogramColumnSelection:
    """Tests for MSChromatogram.get_df() column selection."""

    @pytest.fixture
    def chromatogram_with_data(self):
        """Create a chromatogram with peaks and meta values."""
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID('chrom_1')
        chrom.setMetaValue('FWHM', 5.0)
        chrom.setMetaValue('peak_apex', 100.5)

        # Set precursor
        precursor = pyopenms.Precursor()
        precursor.setMZ(500.0)
        precursor.setCharge(2)
        chrom.setPrecursor(precursor)

        # Set product
        product = pyopenms.Product()
        product.setMZ(300.0)
        chrom.setProduct(product)

        # Set peaks
        rts = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        chrom.set_peaks([rts, ints])

        return chrom

    def test_get_df_columns(self, chromatogram_with_data):
        """Test get_df_columns() returns expected columns."""
        cols = chromatogram_with_data.get_df_columns()

        # Default columns
        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'precursor_mz' in cols
        assert 'precursor_charge' in cols
        assert 'product_mz' in cols
        assert 'native_id' in cols

        # Meta values
        assert 'FWHM' in cols
        assert 'peak_apex' in cols

        # Non-default columns should NOT be present
        assert 'chromatogram_type' not in cols
        assert 'comment' not in cols

    def test_get_df_columns_all(self, chromatogram_with_data):
        """Test get_df_columns('all') returns all columns including non-defaults."""
        cols = chromatogram_with_data.get_df_columns('all')

        # Default columns should be present
        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'precursor_mz' in cols

        # Non-default columns SHOULD be present with 'all'
        assert 'chromatogram_type' in cols
        assert 'comment' in cols

        # Meta values still present
        assert 'FWHM' in cols

    def test_get_df_default(self, chromatogram_with_data):
        """Test get_df() default behavior."""
        df = chromatogram_with_data.get_df()

        assert 'rt' in df.columns
        assert 'intensity' in df.columns
        assert 'precursor_mz' in df.columns
        assert 'native_id' in df.columns
        assert 'FWHM' in df.columns

        # Non-default should NOT be present
        assert 'chromatogram_type' not in df.columns
        assert 'comment' not in df.columns

        assert len(df) == 3
        assert df.loc[0, 'rt'] == 10.0

    def test_get_df_minimal_columns(self, chromatogram_with_data):
        """Test get_df() with minimal columns."""
        df = chromatogram_with_data.get_df(columns=['rt', 'intensity'])

        assert list(df.columns) == ['rt', 'intensity']
        assert len(df) == 3

    def test_get_df_with_non_default_columns(self, chromatogram_with_data):
        """Test get_df() with non-default columns."""
        df = chromatogram_with_data.get_df(columns=['rt', 'intensity', 'chromatogram_type', 'comment'])

        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns

    def test_get_df_with_meta_value(self, chromatogram_with_data):
        """Test get_df() with specific meta value."""
        df = chromatogram_with_data.get_df(columns=['rt', 'intensity', 'FWHM'])

        assert 'FWHM' in df.columns
        assert df.loc[0, 'FWHM'] == 5.0

    def test_get_df_all_columns(self, chromatogram_with_data):
        """Test get_df() with all columns including non-defaults using 'all' parameter."""
        # Use the cleaner API with 'all' parameter
        cols = chromatogram_with_data.get_df_columns('all')
        df = chromatogram_with_data.get_df(columns=cols)

        # Should have all columns
        assert 'rt' in df.columns
        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns
        assert 'FWHM' in df.columns


class TestMobilogramColumnSelection:
    """Tests for Mobilogram.get_df() column selection."""

    @pytest.fixture
    def mobilogram_with_data(self):
        """Create a mobilogram with peaks."""
        mob = pyopenms.Mobilogram()
        mob.setRT(100.0)

        # Set peaks
        mobilities = np.array([0.8, 0.9, 1.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        mob.set_peaks([mobilities, ints])

        return mob

    def test_get_df_columns(self, mobilogram_with_data):
        """Test get_df_columns() returns expected columns."""
        cols = mobilogram_with_data.get_df_columns()

        assert 'mobility' in cols
        assert 'intensity' in cols
        assert 'rt' in cols
        assert 'drift_time_unit' in cols

    def test_get_df_default(self, mobilogram_with_data):
        """Test get_df() default behavior."""
        df = mobilogram_with_data.get_df()

        assert 'mobility' in df.columns
        assert 'intensity' in df.columns
        assert 'rt' in df.columns
        assert len(df) == 3
        assert df.loc[0, 'mobility'] == 0.8
        assert df.loc[0, 'rt'] == 100.0

    def test_get_df_minimal_columns(self, mobilogram_with_data):
        """Test get_df() with minimal columns."""
        df = mobilogram_with_data.get_df(columns=['mobility', 'intensity'])

        assert list(df.columns) == ['mobility', 'intensity']
        assert len(df) == 3


class TestPeptideIdentificationListGetDF:
    """Tests for PeptideIdentificationList.get_df() method."""

    @pytest.fixture
    def peptide_id_list(self):
        """Create a PeptideIdentificationList with test data."""
        pep_list = pyopenms.PeptideIdentificationList()

        # Create first PeptideIdentification with a hit
        pep1 = pyopenms.PeptideIdentification()
        pep1.setRT(100.0)
        pep1.setMZ(500.25)
        pep1.setScoreType('Mascot')
        pep1.setIdentifier('test_id_1')

        hit1 = pyopenms.PeptideHit()
        hit1.setScore(50.0)
        hit1.setCharge(2)
        hit1.setSequence(pyopenms.AASequence.fromString('PEPTIDE'))
        hit1.setMetaValue('target_decoy', 'target')
        pep1.setHits([hit1])

        pep_list.append(pep1)

        # Create second PeptideIdentification with a hit
        pep2 = pyopenms.PeptideIdentification()
        pep2.setRT(200.0)
        pep2.setMZ(600.30)
        pep2.setScoreType('Mascot')
        pep2.setIdentifier('test_id_2')

        hit2 = pyopenms.PeptideHit()
        hit2.setScore(75.0)
        hit2.setCharge(3)
        hit2.setSequence(pyopenms.AASequence.fromString('ANOTHER'))
        hit2.setMetaValue('target_decoy', 'decoy')
        pep2.setHits([hit2])

        pep_list.append(pep2)

        return pep_list

    @pytest.fixture
    def empty_peptide_id_list(self):
        """Create an empty PeptideIdentificationList."""
        return pyopenms.PeptideIdentificationList()

    @pytest.fixture
    def peptide_id_list_with_unidentified(self):
        """Create a PeptideIdentificationList with unidentified entries."""
        pep_list = pyopenms.PeptideIdentificationList()

        # Create identified entry
        pep1 = pyopenms.PeptideIdentification()
        pep1.setRT(100.0)
        pep1.setMZ(500.25)
        pep1.setScoreType('Mascot')
        pep1.setIdentifier('test_id_1')

        hit1 = pyopenms.PeptideHit()
        hit1.setScore(50.0)
        hit1.setCharge(2)
        hit1.setSequence(pyopenms.AASequence.fromString('PEPTIDE'))
        pep1.setHits([hit1])
        pep_list.append(pep1)

        # Create unidentified entry (no hits)
        pep2 = pyopenms.PeptideIdentification()
        pep2.setRT(200.0)
        pep2.setMZ(600.30)
        pep2.setIdentifier('test_id_2')
        # No hits set
        pep_list.append(pep2)

        return pep_list

    def test_get_df_basic(self, peptide_id_list):
        """Test get_df() returns expected DataFrame structure."""
        df = peptide_id_list.get_df()

        assert len(df) == 2
        assert 'rt' in df.columns
        assert 'mz' in df.columns
        assert 'charge' in df.columns
        assert 'P_ID' in df.columns
        assert 'PSM_ID' in df.columns

    def test_get_df_values(self, peptide_id_list):
        """Test get_df() returns correct values."""
        df = peptide_id_list.get_df()

        assert df.loc[0, 'rt'] == pytest.approx(100.0, rel=1e-3)
        assert df.loc[0, 'mz'] == pytest.approx(500.25, rel=1e-3)
        assert df.loc[0, 'charge'] == 2

        assert df.loc[1, 'rt'] == pytest.approx(200.0, rel=1e-3)
        assert df.loc[1, 'mz'] == pytest.approx(600.30, rel=1e-3)
        assert df.loc[1, 'charge'] == 3

    def test_get_df_empty_list(self, empty_peptide_id_list):
        """Test get_df() with empty list."""
        df = empty_peptide_id_list.get_df()
        assert len(df) == 0

    def test_get_df_export_unidentified_true(self, peptide_id_list_with_unidentified):
        """Test get_df() exports unidentified entries when export_unidentified=True."""
        df = peptide_id_list_with_unidentified.get_df(export_unidentified=True)
        assert len(df) == 2

    def test_get_df_export_unidentified_false(self, peptide_id_list_with_unidentified):
        """Test get_df() skips unidentified entries when export_unidentified=False."""
        df = peptide_id_list_with_unidentified.get_df(export_unidentified=False)
        assert len(df) == 1

    def test_get_df_decode_ontology_false(self, peptide_id_list):
        """Test get_df() with decode_ontology=False."""
        df = peptide_id_list.get_df(decode_ontology=False)
        assert len(df) == 2

    def test_get_df_custom_missing_values(self, peptide_id_list_with_unidentified):
        """Test get_df() with custom default_missing_values."""
        custom_missing = {bool: False, int: 0, float: 0.0, str: 'N/A'}
        df = peptide_id_list_with_unidentified.get_df(
            default_missing_values=custom_missing,
            export_unidentified=True
        )
        assert len(df) == 2

    def test_update_scores_from_df(self, peptide_id_list):
        """Test update_scores_from_df() method."""
        df = peptide_id_list.get_df()

        # Modify scores in DataFrame
        df['new_score'] = [100.0, 200.0]

        # Update scores
        peptide_id_list.update_scores_from_df(df, 'new_score')

        # Verify scores were updated
        assert peptide_id_list[0].getHits()[0].getScore() == 100.0
        assert peptide_id_list[1].getHits()[0].getScore() == 200.0

        # Verify score type was set
        assert peptide_id_list[0].getScoreType() == 'new_score'


class TestBackwardsCompatibleFunctions:
    """Tests for backwards-compatible wrapper functions in _dataframes module."""

    @pytest.fixture
    def peptide_id_list(self):
        """Create a PeptideIdentificationList with test data."""
        pep_list = pyopenms.PeptideIdentificationList()

        pep1 = pyopenms.PeptideIdentification()
        pep1.setRT(100.0)
        pep1.setMZ(500.25)
        pep1.setScoreType('Mascot')
        pep1.setIdentifier('test_id_1')

        hit1 = pyopenms.PeptideHit()
        hit1.setScore(50.0)
        hit1.setCharge(2)
        hit1.setSequence(pyopenms.AASequence.fromString('PEPTIDE'))
        pep1.setHits([hit1])

        pep_list.append(pep1)
        return pep_list

    def test_peptide_identifications_to_df_wrapper(self, peptide_id_list):
        """Test that the wrapper function calls get_df() correctly."""
        # Call the wrapper function
        df = pyopenms.peptide_identifications_to_df(peptide_id_list)

        assert len(df) == 1
        assert 'rt' in df.columns
        assert 'mz' in df.columns

    def test_update_scores_from_df_wrapper(self, peptide_id_list):
        """Test that the wrapper function calls update_scores_from_df() correctly."""
        df = pyopenms.peptide_identifications_to_df(peptide_id_list)
        df['new_score'] = [100.0]

        # Call the wrapper function
        result = pyopenms.update_scores_from_df(peptide_id_list, df, 'new_score')

        assert result[0].getHits()[0].getScore() == 100.0


class TestBackwardCompatibility:
    """Tests to ensure backward compatibility with existing code."""

    def test_spectrum_get_df_no_args(self):
        """Test MSSpectrum.get_df() works without arguments."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([10.0, 20.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df = spec.get_df()  # No arguments - should work

        assert 'mz' in df.columns
        assert 'intensity' in df.columns

    def test_spectrum_get_df_export_meta_values_only(self):
        """Test MSSpectrum.get_df(export_meta_values=...) works."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setMetaValue('test', 123)
        mzs = np.array([100.0], dtype=np.float64)
        ints = np.array([10.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df_with = spec.get_df(export_meta_values=True)
        df_without = spec.get_df(export_meta_values=False)

        assert 'test' in df_with.columns
        assert 'test' not in df_without.columns

    def test_chromatogram_get_df_no_args(self):
        """Test MSChromatogram.get_df() works without arguments."""
        chrom = pyopenms.MSChromatogram()
        rts = np.array([10.0, 20.0], dtype=np.float64)
        ints = np.array([100.0, 200.0], dtype=np.float32)
        chrom.set_peaks([rts, ints])

        df = chrom.get_df()  # No arguments - should work

        assert 'rt' in df.columns
        assert 'intensity' in df.columns


class TestToArrowMethods:
    """Tests for to_arrow() methods that export to Apache Arrow Tables."""

    def test_spectrum_to_arrow(self):
        """Test MSSpectrum.to_arrow() returns valid Arrow Table."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(123.45)
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([10.0, 20.0, 30.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        table = spec.to_arrow()

        assert isinstance(table, pa.Table)
        assert table.num_rows == 3
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names

    def test_spectrum_to_arrow_column_selection(self):
        """Test MSSpectrum.to_arrow() with column selection."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        spec = pyopenms.MSSpectrum()
        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([10.0, 20.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        table = spec.to_arrow(columns=['mz', 'intensity'])

        assert table.num_columns == 2
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names
        assert 'rt' not in table.column_names

    def test_chromatogram_to_arrow(self):
        """Test MSChromatogram.to_arrow() returns valid Arrow Table."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        chrom = pyopenms.MSChromatogram()
        rts = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        chrom.set_peaks([rts, ints])

        table = chrom.to_arrow()

        assert isinstance(table, pa.Table)
        assert table.num_rows == 3
        assert 'rt' in table.column_names
        assert 'intensity' in table.column_names

    def test_mobilogram_to_arrow(self):
        """Test Mobilogram.to_arrow() returns valid Arrow Table."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        mob = pyopenms.Mobilogram()
        mobilities = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        mob.set_peaks([mobilities, ints])

        table = mob.to_arrow()

        assert isinstance(table, pa.Table)
        assert table.num_rows == 3
        assert 'mobility' in table.column_names
        assert 'intensity' in table.column_names

    def test_experiment_to_arrow_long_format(self):
        """Test MSExperiment.to_arrow() with long_format=True."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        exp = pyopenms.MSExperiment()
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([10.0, 20.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])
        exp.addSpectrum(spec)

        table = exp.to_arrow(long_format=True)

        assert isinstance(table, pa.Table)
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names
        assert 'ms_level' in table.column_names

    def test_feature_map_to_arrow(self):
        """Test FeatureMap.to_arrow() returns valid Arrow Table."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        fmap = pyopenms.FeatureMap()
        f = pyopenms.Feature()
        f.setRT(100.0)
        f.setMZ(500.0)
        f.setIntensity(1000.0)
        fmap.push_back(f)

        table = fmap.to_arrow()

        assert isinstance(table, pa.Table)
        assert table.num_rows == 1

    def test_to_arrow_to_pandas_roundtrip(self):
        """Test that to_arrow().to_pandas() produces same result as get_df()."""
        pytest.importorskip('pyarrow')

        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([10.0, 20.0, 30.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df_direct = spec.get_df(columns=['mz', 'intensity'])
        df_via_arrow = spec.to_arrow(columns=['mz', 'intensity']).to_pandas()

        assert list(df_direct['mz']) == list(df_via_arrow['mz'])
        assert list(df_direct['intensity']) == list(df_via_arrow['intensity'])


class TestBugFixes:
    """Tests for specific bug fixes."""

    def test_peptide_id_missing_metavalue_no_error(self):
        """Test that missing metavalues don't cause UnboundLocalError.

        Regression test for bug where accessing default_missing_values[type(val)]
        would fail if val was never assigned (metavalue didn't exist).
        """
        pep_list = pyopenms.PeptideIdentificationList()

        # Create two peptide IDs with different metavalues
        pep1 = pyopenms.PeptideIdentification()
        pep1.setRT(100.0)
        pep1.setMZ(500.0)
        hit1 = pyopenms.PeptideHit()
        hit1.setScore(10.0)
        hit1.setCharge(2)
        hit1.setSequence(pyopenms.AASequence.fromString('PEPTIDE'))
        hit1.setMetaValue('custom_score', 0.95)  # Only pep1 has this metavalue
        pep1.setHits([hit1])
        pep_list.append(pep1)

        pep2 = pyopenms.PeptideIdentification()
        pep2.setRT(200.0)
        pep2.setMZ(600.0)
        hit2 = pyopenms.PeptideHit()
        hit2.setScore(20.0)
        hit2.setCharge(3)
        hit2.setSequence(pyopenms.AASequence.fromString('ANOTHER'))
        # hit2 does NOT have 'custom_score' metavalue
        pep2.setHits([hit2])
        pep_list.append(pep2)

        # This should not raise UnboundLocalError
        df = pep_list.get_df()
        assert len(df) == 2

    def test_massql_df_empty_spectrum_no_division_by_zero(self):
        """Test that get_massql_df() handles empty/zero-intensity spectra.

        Regression test for division by zero when normalizing intensities.
        """
        exp = pyopenms.MSExperiment()

        # Add a normal spectrum
        spec1 = pyopenms.MSSpectrum()
        spec1.setMSLevel(1)
        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([1000.0, 2000.0], dtype=np.float32)
        spec1.set_peaks([mzs, ints])
        exp.addSpectrum(spec1)

        # Add an empty spectrum (edge case)
        spec2 = pyopenms.MSSpectrum()
        spec2.setMSLevel(1)
        spec2.set_peaks([np.array([], dtype=np.float64), np.array([], dtype=np.float32)])
        exp.addSpectrum(spec2)

        # This should not raise division by zero warning/error
        ms1_df, ms2_df = exp.get_massql_df()

        # Verify normalized columns are finite (no inf or nan from division)
        if len(ms1_df) > 0:
            assert ms1_df['i_norm'].notna().all() or (ms1_df['i_norm'] == 0).all()
            assert ms1_df['i_tic_norm'].notna().all() or (ms1_df['i_tic_norm'] == 0).all()

    def test_massql_df_zero_intensity_spectrum(self):
        """Test get_massql_df() with spectrum where all intensities are zero."""
        exp = pyopenms.MSExperiment()

        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([0.0, 0.0, 0.0], dtype=np.float32)  # All zeros
        spec.set_peaks([mzs, ints])
        exp.addSpectrum(spec)

        # This should handle division by zero gracefully
        ms1_df, ms2_df = exp.get_massql_df()

        # Normalized values should be 0 (not inf or nan)
        assert len(ms1_df) == 3
        assert (ms1_df['i_norm'] == 0).all()
        assert (ms1_df['i_tic_norm'] == 0).all()
