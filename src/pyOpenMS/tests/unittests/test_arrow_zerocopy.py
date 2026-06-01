"""
Tests for the _arrow_zerocopy nanobind module.

These tests verify zero-copy Arrow export of spectra and chromatograms
from MSExperiment using the Arrow C Data Interface.
"""

import pytest
import numpy as np

pa = pytest.importorskip("pyarrow")

from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow

from pyopenms import MSExperiment, MSSpectrum, MSChromatogram, Precursor


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_experiment():
    """MSExperiment with 3 MS1 spectra, each with a few peaks."""
    exp = MSExperiment()
    for i in range(3):
        spec = MSSpectrum()
        spec.setRT(10.0 + i * 5.0)
        spec.setMSLevel(1)
        spec.setNativeID(f"scan={i}")
        mz = np.array([100.0 + j for j in range(4)], dtype=np.float64)
        intensity = np.array([1000.0 * (j + 1) for j in range(4)], dtype=np.float32)
        spec.set_peaks(mz, intensity)
        exp.addSpectrum(spec)
    return exp


@pytest.fixture
def ms2_experiment():
    """MSExperiment with MS1 and MS2 spectra including precursor info."""
    exp = MSExperiment()

    # MS1
    s1 = MSSpectrum()
    s1.setRT(5.0)
    s1.setMSLevel(1)
    s1.setNativeID("scan=0")
    s1.set_peaks(np.array([100.0, 200.0]), np.array([500.0, 600.0]))
    exp.addSpectrum(s1)

    # MS2 with precursor
    s2 = MSSpectrum()
    s2.setRT(6.0)
    s2.setMSLevel(2)
    s2.setNativeID("scan=1")
    s2.set_peaks(np.array([50.0, 75.0, 90.0]), np.array([100.0, 200.0, 300.0]))
    p = Precursor()
    p.setMZ(200.0)
    p.setCharge(2)
    p.setIntensity(600.0)
    s2.setPrecursors([p])
    exp.addSpectrum(s2)

    return exp


@pytest.fixture
def chrom_experiment():
    """MSExperiment with 2 chromatograms."""
    exp = MSExperiment()

    for i in range(2):
        chrom = MSChromatogram()
        chrom.setNativeID(f"chrom_{i}")
        prec = Precursor()
        prec.setMZ(500.0 + i * 100)
        chrom.setPrecursor(prec)
        rts = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        ints = np.array([10.0, 20.0, 30.0], dtype=np.float32)
        chrom.set_peaks(rts, ints)
        exp.addChromatogram(chrom)

    return exp


# ---------------------------------------------------------------------------
# spectra_to_arrow tests
# ---------------------------------------------------------------------------

class TestSpectraToArrow:

    def test_returns_pyarrow_table(self, simple_experiment):
        result = spectra_to_arrow(simple_experiment)
        assert isinstance(result, pa.Table)

    def test_long_format_row_count(self, simple_experiment):
        """Long format: one row per peak. 3 spectra * 4 peaks = 12 rows."""
        table = spectra_to_arrow(simple_experiment, format='long')
        assert table.num_rows == 12

    def test_long_format_has_expected_columns(self, simple_experiment):
        table = spectra_to_arrow(simple_experiment, format='long')
        names = table.column_names
        for col in ['mz', 'intensity', 'rt']:
            assert col in names, f"Missing column: {col}"

    def test_semi_wide_format(self, simple_experiment):
        """Semi-wide format: one row per spectrum."""
        table = spectra_to_arrow(simple_experiment, format='semi_wide')
        assert table.num_rows == 3

    def test_invalid_format_raises(self, simple_experiment):
        with pytest.raises((ValueError, RuntimeError)):
            spectra_to_arrow(simple_experiment, format='invalid')

    def test_ms_level_filter(self, ms2_experiment):
        """Only export MS2 spectra."""
        table = spectra_to_arrow(ms2_experiment, format='long', ms_levels=[2])
        # MS2 spectrum has 3 peaks
        assert table.num_rows == 3

    def test_ms_level_filter_ms1(self, ms2_experiment):
        table = spectra_to_arrow(ms2_experiment, format='long', ms_levels=[1])
        # MS1 spectrum has 2 peaks
        assert table.num_rows == 2

    def test_rt_filter(self, simple_experiment):
        """Filter by RT range — should only include spectra within range."""
        # Spectra at RT 10, 15, 20. Filter to [14, 16] should get 1 spectrum (4 peaks).
        table = spectra_to_arrow(simple_experiment, format='long',
                                 min_rt=14.0, max_rt=16.0)
        assert table.num_rows == 4

    def test_precursor_info(self, ms2_experiment):
        """MS2 export should include precursor columns."""
        table = spectra_to_arrow(ms2_experiment, format='long',
                                 ms_levels=[2], include_precursor_info=True)
        names = table.column_names
        assert 'precursor_mz' in names

    def test_column_filter(self, simple_experiment):
        """Only request specific columns."""
        table = spectra_to_arrow(simple_experiment, format='long',
                                 columns=['mz', 'intensity'])
        assert set(table.column_names) == {'mz', 'intensity'}

    def test_empty_experiment(self):
        exp = MSExperiment()
        table = spectra_to_arrow(exp, format='long')
        assert table.num_rows == 0

    def test_mz_values_match(self, simple_experiment):
        """Verify exported m/z values match what we put in."""
        table = spectra_to_arrow(simple_experiment, format='long')
        mz_col = table.column('mz').to_pylist()
        # First spectrum peaks: 100, 101, 102, 103
        assert mz_col[:4] == pytest.approx([100.0, 101.0, 102.0, 103.0])


# ---------------------------------------------------------------------------
# chromatograms_to_arrow tests
# ---------------------------------------------------------------------------

class TestChromatogramsToArrow:

    def test_returns_pyarrow_table(self, chrom_experiment):
        result = chromatograms_to_arrow(chrom_experiment)
        assert isinstance(result, pa.Table)

    def test_long_format_row_count(self, chrom_experiment):
        """2 chromatograms * 3 points = 6 rows."""
        table = chromatograms_to_arrow(chrom_experiment, format='long')
        assert table.num_rows == 6

    def test_long_format_columns(self, chrom_experiment):
        table = chromatograms_to_arrow(chrom_experiment, format='long')
        names = table.column_names
        for col in ['rt', 'intensity']:
            assert col in names

    def test_semi_wide_format(self, chrom_experiment):
        table = chromatograms_to_arrow(chrom_experiment, format='semi_wide')
        assert table.num_rows == 2

    def test_rt_filter(self, chrom_experiment):
        """Filter chromatogram points by RT."""
        table = chromatograms_to_arrow(chrom_experiment, format='long',
                                       min_rt=1.5, max_rt=2.5)
        # Each chromatogram has 1 point in [1.5, 2.5] (rt=2.0), so 2 total
        assert table.num_rows == 2

    def test_empty_experiment(self):
        exp = MSExperiment()
        table = chromatograms_to_arrow(exp, format='long')
        assert table.num_rows == 0

    def test_column_filter(self, chrom_experiment):
        table = chromatograms_to_arrow(chrom_experiment, format='long',
                                       columns=['rt', 'intensity'])
        assert set(table.column_names) == {'rt', 'intensity'}


# ---------------------------------------------------------------------------
# Integration: to_arrow() addon uses zero-copy path
# ---------------------------------------------------------------------------

class TestToArrowAddonIntegration:
    """Verify MSExperiment.to_arrow() dispatches to the zero-copy module."""

    def test_to_arrow_spectra(self, simple_experiment):
        simple_experiment.updateRanges()
        table = simple_experiment.to_arrow(data='spectra', format='long')
        assert isinstance(table, pa.Table)
        assert table.num_rows == 12

    def test_to_arrow_chromatograms(self, chrom_experiment):
        chrom_experiment.updateRanges()
        table = chrom_experiment.to_arrow(data='chromatograms', format='long')
        assert isinstance(table, pa.Table)
        assert table.num_rows == 6

    def test_to_arrow_both(self, ms2_experiment):
        # Add a chromatogram so 'both' has data
        chrom = MSChromatogram()
        chrom.setNativeID("chrom_0")
        chrom.set_peaks(np.array([1.0, 2.0]), np.array([10.0, 20.0]))
        ms2_experiment.addChromatogram(chrom)
        ms2_experiment.updateRanges()

        result = ms2_experiment.to_arrow(data='both', format='long')
        assert isinstance(result, dict)
        assert 'spectra' in result
        assert 'chromatograms' in result
        assert isinstance(result['spectra'], pa.Table)
        assert isinstance(result['chromatograms'], pa.Table)
