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
        assert 'time' in cols
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
        assert 'time' in cols
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

        assert 'time' in df.columns
        assert 'intensity' in df.columns
        assert 'precursor_mz' in df.columns
        assert 'native_id' in df.columns
        assert 'FWHM' in df.columns

        # Non-default should NOT be present
        assert 'chromatogram_type' not in df.columns
        assert 'comment' not in df.columns

        assert len(df) == 3
        assert df.loc[0, 'time'] == 10.0

    def test_get_df_minimal_columns(self, chromatogram_with_data):
        """Test get_df() with minimal columns."""
        df = chromatogram_with_data.get_df(columns=['time', 'intensity'])

        assert list(df.columns) == ['time', 'intensity']
        assert len(df) == 3

    def test_get_df_with_non_default_columns(self, chromatogram_with_data):
        """Test get_df() with non-default columns."""
        df = chromatogram_with_data.get_df(columns=['time', 'intensity', 'chromatogram_type', 'comment'])

        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns

    def test_get_df_with_meta_value(self, chromatogram_with_data):
        """Test get_df() with specific meta value."""
        df = chromatogram_with_data.get_df(columns=['time', 'intensity', 'FWHM'])

        assert 'FWHM' in df.columns
        assert df.loc[0, 'FWHM'] == 5.0

    def test_get_df_all_columns(self, chromatogram_with_data):
        """Test get_df() with all columns including non-defaults using 'all' parameter."""
        # Use the cleaner API with 'all' parameter
        cols = chromatogram_with_data.get_df_columns('all')
        df = chromatogram_with_data.get_df(columns=cols)

        # Should have all columns
        assert 'time' in df.columns
        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns
        assert 'FWHM' in df.columns


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

        assert 'time' in df.columns
        assert 'intensity' in df.columns
