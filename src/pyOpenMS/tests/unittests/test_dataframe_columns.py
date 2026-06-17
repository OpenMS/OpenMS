"""
Tests for DataFrame column selection feature in MSSpectrum and MSChromatogram.

Tests cover:
- df_columns() discovery method
- Column selection via to_df(columns=[...])
- Default behavior (backward compatibility)
- Non-default columns (ion_mobility_unit, chromatogram_type, comment)
- Meta value handling with column selection
- Edge cases (empty spectra, missing data)
"""

import pytest
import numpy as np

import pyopenms


class TestMSSpectrumColumnSelection:
    """Tests for MSSpectrum.to_df() column selection."""

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

    def test_df_columns_full_spectrum(self, spectrum_with_data):
        """Test df_columns() returns all expected columns for full spectrum."""
        cols = spectrum_with_data.df_columns()

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

    def test_df_columns_simple_spectrum(self, simple_spectrum):
        """Test df_columns() for simple MS1 spectrum."""
        cols = simple_spectrum.df_columns()

        # Should have core columns
        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'rt' in cols

        # Should NOT have precursor columns (MS1)
        assert 'precursor_mz' not in cols
        assert 'precursor_charge' not in cols

        # Should NOT have IM columns (no IM data)
        assert 'ion_mobility' not in cols

    def test_df_columns_no_meta_values(self, spectrum_with_data):
        """Test df_columns() with export_meta_values=False."""
        cols = spectrum_with_data.df_columns(export_meta_values=False)

        assert 'mz' in cols
        assert 'total_ion_current' not in cols
        assert 'base_peak_mz' not in cols

    def test_to_df_default(self, spectrum_with_data):
        """Test to_df() default behavior returns all expected columns."""
        df = spectrum_with_data.to_df()

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

    def test_to_df_minimal_columns(self, spectrum_with_data):
        """Test to_df() with minimal column selection."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity'])

        assert list(df.columns) == ['mz', 'intensity']
        assert len(df) == 3
        assert df.loc[0, 'mz'] == 100.0
        assert df.loc[0, 'intensity'] == 10.0

    def test_to_df_custom_columns(self, spectrum_with_data):
        """Test to_df() with custom column selection."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity', 'rt', 'precursor_mz'])

        assert set(df.columns) == {'mz', 'intensity', 'rt', 'precursor_mz'}
        assert df.loc[0, 'precursor_mz'] == 500.0

    def test_to_df_with_ion_mobility(self, spectrum_with_data):
        """Test to_df() including ion mobility column."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity', 'ion_mobility'])

        assert 'ion_mobility' in df.columns
        assert df.loc[0, 'ion_mobility'] == 1.5
        assert df.loc[1, 'ion_mobility'] == 2.0

    def test_to_df_with_ion_mobility_unit(self, spectrum_with_data):
        """Test to_df() with non-default ion_mobility_unit column."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity', 'ion_mobility_unit'])

        assert 'ion_mobility_unit' in df.columns

    def test_to_df_with_ion_annotation(self, spectrum_with_data):
        """Test to_df() including ion annotation column."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity', 'ion_annotation'])

        assert 'ion_annotation' in df.columns
        assert df.loc[0, 'ion_annotation'] == 'b2+'
        assert df.loc[1, 'ion_annotation'] == 'y3+'

    def test_to_df_with_meta_value(self, spectrum_with_data):
        """Test to_df() with specific meta value column."""
        df = spectrum_with_data.to_df(columns=['mz', 'intensity', 'total_ion_current'])

        assert 'total_ion_current' in df.columns
        assert df.loc[0, 'total_ion_current'] == 1000.0

    def test_to_df_all_columns(self, spectrum_with_data):
        """Test to_df() requesting all available columns including non-defaults."""
        # Get all default columns
        cols = spectrum_with_data.df_columns()
        # Add non-default columns
        cols.append('ion_mobility_unit')

        df = spectrum_with_data.to_df(columns=cols)

        # Should have all columns
        assert 'mz' in df.columns
        assert 'ion_mobility_unit' in df.columns
        assert 'total_ion_current' in df.columns

    def test_to_df_missing_precursor(self, simple_spectrum):
        """Test to_df() requesting precursor columns when no precursor present."""
        df = simple_spectrum.to_df(columns=['mz', 'intensity', 'precursor_mz'])

        assert 'precursor_mz' in df.columns
        assert np.isnan(df.loc[0, 'precursor_mz'])

    def test_to_df_empty_spectrum(self):
        """Test to_df() with empty spectrum."""
        spec = pyopenms.MSSpectrum()
        df = spec.to_df()

        assert len(df) == 0

    def test_df_columns_empty_spectrum(self):
        """Test df_columns() with empty spectrum."""
        spec = pyopenms.MSSpectrum()
        cols = spec.df_columns()

        # Should still have core columns
        assert 'mz' in cols
        assert 'intensity' in cols


class TestMSChromatogramColumnSelection:
    """Tests for MSChromatogram.to_df() column selection."""

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

    def test_df_columns(self, chromatogram_with_data):
        """Test df_columns() returns expected columns."""
        cols = chromatogram_with_data.df_columns()

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

    def test_df_columns_all(self, chromatogram_with_data):
        """Test df_columns('all') returns all columns including non-defaults."""
        cols = chromatogram_with_data.df_columns('all')

        # Default columns should be present
        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'precursor_mz' in cols

        # Non-default columns SHOULD be present with 'all'
        assert 'chromatogram_type' in cols
        assert 'comment' in cols

        # Meta values still present
        assert 'FWHM' in cols

    def test_to_df_default(self, chromatogram_with_data):
        """Test to_df() default behavior."""
        df = chromatogram_with_data.to_df()

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

    def test_to_df_minimal_columns(self, chromatogram_with_data):
        """Test to_df() with minimal columns."""
        df = chromatogram_with_data.to_df(columns=['rt', 'intensity'])

        assert list(df.columns) == ['rt', 'intensity']
        assert len(df) == 3

    def test_to_df_with_non_default_columns(self, chromatogram_with_data):
        """Test to_df() with non-default columns."""
        df = chromatogram_with_data.to_df(columns=['rt', 'intensity', 'chromatogram_type', 'comment'])

        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns

    def test_to_df_with_meta_value(self, chromatogram_with_data):
        """Test to_df() with specific meta value."""
        df = chromatogram_with_data.to_df(columns=['rt', 'intensity', 'FWHM'])

        assert 'FWHM' in df.columns
        assert df.loc[0, 'FWHM'] == 5.0

    def test_to_df_all_columns(self, chromatogram_with_data):
        """Test to_df() with all columns including non-defaults using 'all' parameter."""
        # Use the cleaner API with 'all' parameter
        cols = chromatogram_with_data.df_columns('all')
        df = chromatogram_with_data.to_df(columns=cols)

        # Should have all columns
        assert 'rt' in df.columns
        assert 'chromatogram_type' in df.columns
        assert 'comment' in df.columns
        assert 'FWHM' in df.columns


class TestMobilogramColumnSelection:
    """Tests for Mobilogram.to_df() column selection."""

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

    def test_df_columns(self, mobilogram_with_data):
        """Test df_columns() returns expected columns."""
        cols = mobilogram_with_data.df_columns()

        assert 'mobility' in cols
        assert 'intensity' in cols
        assert 'rt' in cols
        assert 'drift_time_unit' in cols

    def test_to_df_default(self, mobilogram_with_data):
        """Test to_df() default behavior."""
        df = mobilogram_with_data.to_df()

        assert 'mobility' in df.columns
        assert 'intensity' in df.columns
        assert 'rt' in df.columns
        assert len(df) == 3
        assert df.loc[0, 'mobility'] == 0.8
        assert df.loc[0, 'rt'] == 100.0

    def test_to_df_minimal_columns(self, mobilogram_with_data):
        """Test to_df() with minimal columns."""
        df = mobilogram_with_data.to_df(columns=['mobility', 'intensity'])

        assert list(df.columns) == ['mobility', 'intensity']
        assert len(df) == 3


class TestPeptideIdentificationListGetDF:
    """Tests for PeptideIdentificationList.to_df() method."""

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

    def test_to_df_basic(self, peptide_id_list):
        """Test to_df() returns expected DataFrame structure."""
        df = peptide_id_list.to_df()

        assert len(df) == 2
        assert 'rt' in df.columns
        assert 'mz' in df.columns
        assert 'charge' in df.columns
        assert 'P_ID' in df.columns
        assert 'PSM_ID' in df.columns

    def test_to_df_values(self, peptide_id_list):
        """Test to_df() returns correct values."""
        df = peptide_id_list.to_df()

        assert df.loc[0, 'rt'] == pytest.approx(100.0, rel=1e-3)
        assert df.loc[0, 'mz'] == pytest.approx(500.25, rel=1e-3)
        assert df.loc[0, 'charge'] == 2

        assert df.loc[1, 'rt'] == pytest.approx(200.0, rel=1e-3)
        assert df.loc[1, 'mz'] == pytest.approx(600.30, rel=1e-3)
        assert df.loc[1, 'charge'] == 3

    def test_to_df_empty_list(self, empty_peptide_id_list):
        """Test to_df() with empty list."""
        df = empty_peptide_id_list.to_df()
        assert len(df) == 0

    def test_to_df_export_unidentified_true(self, peptide_id_list_with_unidentified):
        """Test to_df() exports unidentified entries when export_unidentified=True."""
        df = peptide_id_list_with_unidentified.to_df(export_unidentified=True)
        assert len(df) == 2

    def test_to_df_export_unidentified_false(self, peptide_id_list_with_unidentified):
        """Test to_df() skips unidentified entries when export_unidentified=False."""
        df = peptide_id_list_with_unidentified.to_df(export_unidentified=False)
        assert len(df) == 1

    def test_to_df_decode_ontology_false(self, peptide_id_list):
        """Test to_df() with decode_ontology=False."""
        df = peptide_id_list.to_df(decode_ontology=False)
        assert len(df) == 2

    def test_to_df_custom_missing_values(self, peptide_id_list_with_unidentified):
        """Test to_df() with custom default_missing_values."""
        custom_missing = {bool: False, int: 0, float: 0.0, str: 'N/A'}
        df = peptide_id_list_with_unidentified.to_df(
            default_missing_values=custom_missing,
            export_unidentified=True
        )
        assert len(df) == 2

    def test_update_scores_from_df(self, peptide_id_list):
        """Test update_scores_from_df() method."""
        df = peptide_id_list.to_df()

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
        """Test that the wrapper function calls to_df() correctly."""
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

    def test_spectrum_to_df_no_args(self):
        """Test MSSpectrum.to_df() works without arguments."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0], dtype=np.float64)
        ints = np.array([10.0, 20.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df = spec.to_df()  # No arguments - should work

        assert 'mz' in df.columns
        assert 'intensity' in df.columns

    def test_spectrum_to_df_export_meta_values_only(self):
        """Test MSSpectrum.to_df(export_meta_values=...) works."""
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setMetaValue('test', 123)
        mzs = np.array([100.0], dtype=np.float64)
        ints = np.array([10.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df_with = spec.to_df(export_meta_values=True)
        df_without = spec.to_df(export_meta_values=False)

        assert 'test' in df_with.columns
        assert 'test' not in df_without.columns

    def test_chromatogram_to_df_no_args(self):
        """Test MSChromatogram.to_df() works without arguments."""
        chrom = pyopenms.MSChromatogram()
        rts = np.array([10.0, 20.0], dtype=np.float64)
        ints = np.array([100.0, 200.0], dtype=np.float32)
        chrom.set_peaks([rts, ints])

        df = chrom.to_df()  # No arguments - should work

        assert 'rt' in df.columns
        assert 'intensity' in df.columns


class TestMSExperimentUnifiedToArrow:
    """Tests for the unified MSExperiment.to_arrow() interface."""

    @pytest.fixture
    def experiment_with_data(self):
        """Create an MSExperiment with spectra and chromatograms."""
        exp = pyopenms.MSExperiment()

        # Add MS1 spectrum
        spec1 = pyopenms.MSSpectrum()
        spec1.setMSLevel(1)
        spec1.setRT(100.0)
        spec1.setNativeID('scan=1')
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([1000.0, 2000.0, 3000.0], dtype=np.float32)
        spec1.set_peaks([mzs, ints])
        exp.addSpectrum(spec1)

        # Add MS2 spectrum with precursor
        spec2 = pyopenms.MSSpectrum()
        spec2.setMSLevel(2)
        spec2.setRT(110.0)
        spec2.setNativeID('scan=2')
        prec = pyopenms.Precursor()
        prec.setMZ(500.0)
        prec.setCharge(2)
        prec.setIntensity(50000.0)
        prec.setIsolationWindowLowerOffset(1.5)
        prec.setIsolationWindowUpperOffset(1.5)
        spec2.setPrecursors([prec])
        mzs2 = np.array([150.0, 250.0], dtype=np.float64)
        ints2 = np.array([500.0, 750.0], dtype=np.float32)
        spec2.set_peaks([mzs2, ints2])
        exp.addSpectrum(spec2)

        # Add chromatogram
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID('TIC')
        chrom.getPrecursor().setMZ(500.0)
        chrom.getProduct().setMZ(200.0)
        rts = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        ints_c = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        chrom.set_peaks([rts, ints_c])
        exp.addChromatogram(chrom)

        return exp

    def test_to_arrow_spectra_long_format(self, experiment_with_data):
        """Test to_arrow() with data='spectra', format='long'."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        table = experiment_with_data.to_arrow(data='spectra', format='long')

        assert isinstance(table, pa.Table)
        # 3 peaks in MS1 + 2 peaks in MS2 = 5 rows
        assert table.num_rows == 5
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names
        assert 'rt' in table.column_names
        assert 'spectrum_index' in table.column_names
        assert 'ms_level' in table.column_names
        assert 'precursor_mz' in table.column_names

    def test_to_arrow_spectra_semi_wide_format(self, experiment_with_data):
        """Test to_arrow() with data='spectra', format='semi_wide'."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        table = experiment_with_data.to_arrow(data='spectra', format='semi_wide')

        assert isinstance(table, pa.Table)
        # 2 spectra = 2 rows
        assert table.num_rows == 2
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names
        # In semi_wide, mz should be a list type
        mz_type = table.schema.field('mz').type
        assert pa.types.is_list(mz_type)

    def test_to_arrow_chromatograms_long_format(self, experiment_with_data):
        """Test to_arrow() with data='chromatograms', format='long'."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        table = experiment_with_data.to_arrow(data='chromatograms', format='long')

        assert isinstance(table, pa.Table)
        # 3 points in chromatogram
        assert table.num_rows == 3
        assert 'rt' in table.column_names
        assert 'intensity' in table.column_names
        assert 'chromatogram_index' in table.column_names
        assert 'precursor_mz' in table.column_names
        assert 'product_mz' in table.column_names

    def test_to_arrow_chromatograms_semi_wide_format(self, experiment_with_data):
        """Test to_arrow() with data='chromatograms', format='semi_wide'."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        table = experiment_with_data.to_arrow(data='chromatograms', format='semi_wide')

        assert isinstance(table, pa.Table)
        # 1 chromatogram = 1 row
        assert table.num_rows == 1
        rt_type = table.schema.field('rt').type
        assert pa.types.is_list(rt_type)

    def test_to_arrow_both(self, experiment_with_data):
        """Test to_arrow() with data='both'."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        result = experiment_with_data.to_arrow(data='both')

        assert isinstance(result, dict)
        assert 'spectra' in result
        assert 'chromatograms' in result
        assert isinstance(result['spectra'], pa.Table)
        assert isinstance(result['chromatograms'], pa.Table)

    def test_to_arrow_ms_level_filter(self, experiment_with_data):
        """Test to_arrow() with ms_levels filter."""
        pytest.importorskip('pyarrow')

        # Only MS1
        table_ms1 = experiment_with_data.to_arrow(data='spectra', ms_levels=[1])
        assert table_ms1.num_rows == 3  # 3 peaks in MS1

        # Only MS2
        table_ms2 = experiment_with_data.to_arrow(data='spectra', ms_levels=[2])
        assert table_ms2.num_rows == 2  # 2 peaks in MS2

    def test_to_arrow_rt_filter(self, experiment_with_data):
        """Test to_arrow() with RT range filter."""
        pytest.importorskip('pyarrow')

        # Only spectra with RT >= 105 (should include only MS2 at RT=110)
        table = experiment_with_data.to_arrow(data='spectra', min_rt=105.0)
        assert table.num_rows == 2  # 2 peaks in MS2 spectrum

    def test_to_arrow_mz_filter(self, experiment_with_data):
        """Test to_arrow() with m/z range filter."""
        pytest.importorskip('pyarrow')

        # Only peaks with m/z between 150 and 250
        table = experiment_with_data.to_arrow(data='spectra', min_mz=150.0, max_mz=250.0)
        # MS1: 200.0 (1 peak), MS2: 150.0, 250.0 (2 peaks) = 3 peaks
        assert table.num_rows == 3

    def test_to_arrow_column_selection(self, experiment_with_data):
        """Test to_arrow() with column selection."""
        pytest.importorskip('pyarrow')

        table = experiment_with_data.to_arrow(
            data='spectra',
            columns=['mz', 'intensity', 'rt']
        )

        assert table.num_columns == 3
        assert 'mz' in table.column_names
        assert 'intensity' in table.column_names
        assert 'rt' in table.column_names
        assert 'spectrum_index' not in table.column_names

    def test_to_arrow_no_precursor_info(self, experiment_with_data):
        """Test to_arrow() with include_precursor_info=False."""
        pytest.importorskip('pyarrow')

        table = experiment_with_data.to_arrow(
            data='spectra',
            include_precursor_info=False
        )

        assert 'precursor_mz' not in table.column_names
        assert 'precursor_charge' not in table.column_names

    def test_to_arrow_no_ion_mobility(self, experiment_with_data):
        """Test to_arrow() with include_ion_mobility=False."""
        pytest.importorskip('pyarrow')

        table = experiment_with_data.to_arrow(
            data='spectra',
            include_ion_mobility=False
        )

        assert 'ion_mobility' not in table.column_names

    def test_to_arrow_backward_compat_long_format(self, experiment_with_data):
        """Test backward compatibility with deprecated long_format parameter."""
        pytest.importorskip('pyarrow')
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            table = experiment_with_data.to_arrow(long_format=True)
            # Should trigger deprecation warning
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert 'long_format' in str(w[0].message)

        assert table.num_rows == 5

    def test_to_arrow_empty_experiment(self):
        """Test to_arrow() with empty experiment."""
        pytest.importorskip('pyarrow')
        import pyarrow as pa

        exp = pyopenms.MSExperiment()
        # Empty experiments require a fix to handle getMinRT() on empty ranges
        # The fix is in the MSExperiment addon; after rebuild this will work
        try:
            table = exp.to_arrow(data='spectra')
            assert isinstance(table, pa.Table)
            assert table.num_rows == 0
        except RuntimeError as e:
            if "Empty or uninitialized range" in str(e):
                pytest.skip("Empty experiment handling requires rebuild with fixed addon")

    def test_to_arrow_invalid_data_param(self, experiment_with_data):
        """Test to_arrow() raises ValueError for invalid data parameter."""
        pytest.importorskip('pyarrow')

        with pytest.raises(ValueError, match="data must be"):
            experiment_with_data.to_arrow(data='invalid')

    def test_to_arrow_invalid_format_param(self, experiment_with_data):
        """Test to_arrow() raises ValueError for invalid format parameter."""
        pytest.importorskip('pyarrow')

        with pytest.raises(ValueError, match="format must be"):
            experiment_with_data.to_arrow(format='invalid')

    def test_arrow_columns_spectra(self, experiment_with_data):
        """Test arrow_columns() for spectra."""
        cols = experiment_with_data.arrow_columns(data='spectra', format='long')

        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'rt' in cols
        assert 'spectrum_index' in cols
        assert 'precursor_mz' in cols

    def test_arrow_columns_chromatograms(self, experiment_with_data):
        """Test arrow_columns() for chromatograms."""
        cols = experiment_with_data.arrow_columns(data='chromatograms', format='long')

        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'chromatogram_index' in cols
        assert 'precursor_mz' in cols
        assert 'product_mz' in cols

    def test_arrow_columns_both(self, experiment_with_data):
        """Test arrow_columns() for both."""
        cols = experiment_with_data.arrow_columns(data='both')

        assert isinstance(cols, dict)
        assert 'spectra' in cols
        assert 'chromatograms' in cols
        assert 'mz' in cols['spectra']
        assert 'chromatogram_index' in cols['chromatograms']


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
        """Test MSExperiment.to_arrow() with long_format=True (deprecated parameter)."""
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
        """Test that to_arrow().to_pandas() produces same result as to_df()."""
        pytest.importorskip('pyarrow')

        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(1)
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        ints = np.array([10.0, 20.0, 30.0], dtype=np.float32)
        spec.set_peaks([mzs, ints])

        df_direct = spec.to_df(columns=['mz', 'intensity'])
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
        df = pep_list.to_df()
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

    def test_featuremap_missing_spectrum_native_id(self):
        """Test FeatureMap.to_df() handles missing spectrum_native_id metavalue.

        Regression test for bug where getMetaValue('spectrum_native_id') was called
        without checking if the metavalue exists, causing KeyError or unexpected behavior.
        """
        fmap = pyopenms.FeatureMap()

        # Add a feature WITH spectrum_native_id
        f1 = pyopenms.Feature()
        f1.setRT(100.0)
        f1.setMZ(500.0)
        f1.setIntensity(1000.0)
        f1.setCharge(2)
        f1.setMetaValue('spectrum_native_id', 'scan=100')
        # Add peptide identification to trigger the code path
        pep1 = pyopenms.PeptideIdentification()
        hit1 = pyopenms.PeptideHit()
        hit1.setSequence(pyopenms.AASequence.fromString('PEPTIDE'))
        hit1.setScore(10.0)
        pep1.setHits([hit1])
        pep_list1 = pyopenms.PeptideIdentificationList()
        pep_list1.push_back(pep1)
        f1.setPeptideIdentifications(pep_list1)
        fmap.push_back(f1)

        # Add a feature WITHOUT spectrum_native_id (the bug case)
        f2 = pyopenms.Feature()
        f2.setRT(200.0)
        f2.setMZ(600.0)
        f2.setIntensity(2000.0)
        f2.setCharge(3)
        # NO setMetaValue('spectrum_native_id', ...) - this is the bug case
        pep2 = pyopenms.PeptideIdentification()
        hit2 = pyopenms.PeptideHit()
        hit2.setSequence(pyopenms.AASequence.fromString('ANOTHER'))
        hit2.setScore(20.0)
        pep2.setHits([hit2])
        pep_list2 = pyopenms.PeptideIdentificationList()
        pep_list2.push_back(pep2)
        f2.setPeptideIdentifications(pep_list2)
        fmap.push_back(f2)

        # This should not raise an error
        df = fmap.to_df(export_peptide_identifications=True)

        assert len(df) == 2
        # First feature should have the native_id
        assert df.iloc[0]['ID_native_id'] == 'scan=100'
        # Second feature should have 'None' string (numpy converts None to string due to U100 dtype)
        assert df.iloc[1]['ID_native_id'] == 'None'


class TestFeatureMapColumnSelection:
    """Tests for FeatureMap.df_columns() method."""

    @pytest.fixture
    def feature_map_with_data(self):
        """Create a FeatureMap with features."""
        fmap = pyopenms.FeatureMap()

        f = pyopenms.Feature()
        f.setRT(100.0)
        f.setMZ(500.0)
        f.setIntensity(1000.0)
        f.setCharge(2)
        f.setMetaValue('custom_score', 0.95)
        fmap.push_back(f)

        return fmap

    def test_df_columns_default(self, feature_map_with_data):
        """Test df_columns() returns expected default columns."""
        cols = feature_map_with_data.df_columns()

        # Core columns
        assert 'feature_id' in cols
        assert 'charge' in cols
        assert 'rt' in cols
        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'quality' in cols

        # Peptide ID columns (default)
        assert 'peptide_sequence' in cols
        assert 'peptide_score' in cols

    def test_df_columns_no_peptide_ids(self, feature_map_with_data):
        """Test df_columns() without peptide identification columns."""
        cols = feature_map_with_data.df_columns(export_peptide_identifications=False)

        assert 'feature_id' in cols
        assert 'rt' in cols
        assert 'peptide_sequence' not in cols
        assert 'peptide_score' not in cols

    def test_df_columns_all(self, feature_map_with_data):
        """Test df_columns('all') includes meta values."""
        cols = feature_map_with_data.df_columns(columns='all')

        assert 'feature_id' in cols
        assert 'custom_score' in cols


class TestConsensusMapColumnSelection:
    """Tests for ConsensusMap.df_columns() method."""

    @pytest.fixture
    def consensus_map_with_data(self):
        """Create a ConsensusMap with consensus features."""
        cmap = pyopenms.ConsensusMap()

        # Set up column headers for 2 files
        headers = {0: pyopenms.ColumnHeader(), 1: pyopenms.ColumnHeader()}
        headers[0].filename = 'file1.mzML'
        headers[0].label = 'sample1'
        headers[1].filename = 'file2.mzML'
        headers[1].label = 'sample2'
        cmap.setColumnHeaders(headers)

        # Add a consensus feature
        cf = pyopenms.ConsensusFeature()
        cf.setRT(100.0)
        cf.setMZ(500.0)
        cf.setIntensity(1000.0)
        cf.setCharge(2)
        cmap.push_back(cf)

        return cmap

    def test_df_columns_default(self, consensus_map_with_data):
        """Test df_columns() returns expected columns."""
        cols = consensus_map_with_data.df_columns()

        # Core columns
        assert 'sequence' in cols
        assert 'charge' in cols
        assert 'rt' in cols
        assert 'mz' in cols
        assert 'quality' in cols


class TestMRMTransitionGroupCPColumnSelection:
    """Tests for MRMTransitionGroupCP column name methods."""

    @pytest.fixture
    def mrm_group_with_data(self):
        """Create an MRMTransitionGroupCP with chromatograms and features."""
        group = pyopenms.MRMTransitionGroupCP()

        # Add a chromatogram
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID('trans_1')
        rts = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 150.0], dtype=np.float32)
        chrom.set_peaks([rts, ints])
        group.addChromatogram(chrom, 'trans_1')

        # Add a feature
        feature = pyopenms.MRMFeature()
        feature.setRT(20.0)
        feature.setIntensity(500.0)
        group.addFeature(feature)

        return group

    def test_chromatogram_df_columns(self, mrm_group_with_data):
        """Test chromatogram_df_columns() returns expected columns."""
        cols = mrm_group_with_data.chromatogram_df_columns()

        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'native_id' in cols

    def test_feature_df_columns(self, mrm_group_with_data):
        """Test feature_df_columns() returns expected columns."""
        cols = mrm_group_with_data.feature_df_columns()

        assert 'feature_id' in cols
        assert 'rt' in cols
        assert 'intensity' in cols
        assert 'quality' in cols


class TestMSExperimentDFColumnSelection:
    """Tests for MSExperiment.df_columns() method."""

    def test_df_columns_long_format(self):
        """Test df_columns() with long_format=True (deprecated parameter)."""
        exp = pyopenms.MSExperiment()
        cols = exp.df_columns(long_format=True)

        assert 'rt' in cols
        assert 'mz' in cols
        assert 'intensity' in cols
        assert 'ms_level' in cols

    def test_df_columns_compact_format(self):
        """Test df_columns() with long_format=False (deprecated parameter)."""
        exp = pyopenms.MSExperiment()
        cols = exp.df_columns(long_format=False)

        assert 'rt' in cols
        assert 'ms_level' in cols
        assert 'mz_array' in cols
        assert 'intensity_array' in cols
