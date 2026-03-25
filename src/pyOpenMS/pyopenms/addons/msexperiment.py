"""Addon methods for MSExperiment class."""

from __future__ import annotations
import warnings
import numpy as np
from . import addon


@addon("MSExperiment")
def df_columns(self, long_format=False):
    """Returns a list of column names that to_df() would produce."""
    if long_format:
        return ['rt', 'mz', 'intensity', 'ms_level']
    else:
        return ['rt', 'ms_level', 'mz_array', 'intensity_array']


@addon("MSExperiment")
def to_df(self, columns=None, ms_levels=None, long_format=False):
    """Returns a pandas DataFrame of all peaks in the MSExperiment."""
    import pandas as pd

    if ms_levels is None:
        ms_levels = []
    self.updateRanges()
    if not ms_levels:
        ms_levels = self.getMSLevels()

    if long_format:
        # Fast C++ path: reuse get2DPeakDataSemiWide + np.repeat to expand
        # per-spectrum scalars to per-peak arrays
        rt_arr, ms_arr, mz_flat, int_flat, offsets = self.get2DPeakDataSemiWide(
            list(ms_levels)
        )
        counts = np.diff(offsets)
        df = pd.DataFrame({
            'rt': np.repeat(rt_arr, counts),
            'mz': mz_flat,
            'intensity': int_flat,
            'ms_level': np.repeat(ms_arr, counts),
        })
    else:
        # Fast C++ path: flat buffers + offsets, then slice for zero-copy views
        rt_arr, ms_arr, mz_flat, int_flat, offsets = self.get2DPeakDataSemiWide(
            list(ms_levels)
        )
        n = len(rt_arr)
        mz_arrays = [mz_flat[offsets[i]:offsets[i+1]] for i in range(n)]
        int_arrays = [int_flat[offsets[i]:offsets[i+1]] for i in range(n)]
        df = pd.DataFrame({
            'rt': rt_arr,
            'ms_level': ms_arr,
            'mz_array': mz_arrays,
            'intensity_array': int_arrays,
        })

    if columns is not None:
        available_cols = [c for c in columns if c in df.columns]
        df = df[available_cols]
    return df


@addon("MSExperiment")
def arrow_columns(self, data='spectra', format='long',
                  include_precursor_info=True, include_ion_mobility=True):
    """Returns column names that to_arrow() would produce."""
    spectra_long_cols = ['mz', 'intensity', 'rt']
    if include_ion_mobility:
        spectra_long_cols.append('ion_mobility')
    spectra_long_cols.extend(['spectrum_index', 'ms_level', 'native_id'])
    if include_precursor_info:
        spectra_long_cols.extend(['precursor_mz', 'precursor_charge', 'precursor_intensity',
                                  'isolation_lower', 'isolation_upper'])

    spectra_semi_wide_cols = ['spectrum_index', 'rt', 'ms_level', 'native_id', 'mz', 'intensity']
    if include_ion_mobility:
        spectra_semi_wide_cols.append('ion_mobility')
    if include_precursor_info:
        spectra_semi_wide_cols.extend(['precursor_mz', 'precursor_charge', 'precursor_intensity',
                                       'isolation_lower', 'isolation_upper'])

    chrom_long_cols = ['rt', 'intensity', 'chromatogram_index', 'native_id', 'precursor_mz', 'product_mz']
    chrom_semi_wide_cols = ['chromatogram_index', 'native_id', 'rt', 'intensity', 'precursor_mz', 'product_mz']

    if data == 'spectra':
        return spectra_long_cols if format == 'long' else spectra_semi_wide_cols
    elif data == 'chromatograms':
        return chrom_long_cols if format == 'long' else chrom_semi_wide_cols
    elif data == 'both':
        return {
            'spectra': spectra_long_cols if format == 'long' else spectra_semi_wide_cols,
            'chromatograms': chrom_long_cols if format == 'long' else chrom_semi_wide_cols
        }
    else:
        raise ValueError(f"data must be 'spectra', 'chromatograms', or 'both', got '{data}'")


@addon("MSExperiment")
def to_arrow(self, data='spectra', format='long', columns=None,
             ms_levels=None, min_rt=None, max_rt=None, min_mz=None,
             max_mz=None, include_precursor_info=True,
             include_ion_mobility=True, long_format=None):
    """Returns an Apache Arrow Table with spectra and/or chromatogram data."""
    import pyarrow as pa

    if long_format is not None:
        import warnings
        warnings.warn(
            "long_format parameter is deprecated. Use format='long' or format='semi_wide' instead.",
            DeprecationWarning, stacklevel=2
        )
        format = 'long' if long_format else 'semi_wide'

    if format not in ('long', 'semi_wide'):
        raise ValueError(f"format must be 'long' or 'semi_wide', got '{format}'")

    if data not in ('spectra', 'chromatograms', 'both'):
        raise ValueError(f"data must be 'spectra', 'chromatograms', or 'both', got '{data}'")

    if ms_levels is None:
        ms_levels = []
    else:
        ms_levels = list(ms_levels)

    if min_rt is None:
        min_rt = 0.0
    if max_rt is None:
        max_rt = float('inf')
    if min_mz is None:
        min_mz = 0.0
    if max_mz is None:
        max_mz = float('inf')

    if self.getNrSpectra() == 0 and self.getNrChromatograms() == 0:
        result = {}
        if data in ('spectra', 'both'):
            result['spectra'] = _build_spectra_arrow(
                self, format, columns, ms_levels,
                min_rt, max_rt, min_mz, max_mz,
                include_precursor_info, include_ion_mobility, pa
            )
        if data in ('chromatograms', 'both'):
            result['chromatograms'] = _build_chrom_arrow(
                self, format, columns, min_rt, max_rt, pa
            )
        if data == 'both':
            return result
        if data == 'spectra':
            return result['spectra']
        return result['chromatograms']

    self.updateRanges()

    if not ms_levels:
        ms_levels = self.getMSLevels() if self.getNrSpectra() > 0 else []

    # Try zero-copy C++ path first (requires Arrow/Parquet)
    try:
        from pyopenms._arrow_zerocopy import spectra_to_arrow, chromatograms_to_arrow
        _use_zerocopy = True
    except ImportError:
        _use_zerocopy = False
        warnings.warn(
            "pyopenms._arrow_zerocopy not found — Arrow/Parquet are required; "
            "this indicates a broken or incomplete installation, "
            "falling back to slower Python export.",
            stacklevel=2,
        )

    if _use_zerocopy:
        result = {}
        if data in ('spectra', 'both'):
            result['spectra'] = spectra_to_arrow(
                self, format=format, ms_levels=ms_levels,
                min_rt=min_rt, max_rt=max_rt, min_mz=min_mz, max_mz=max_mz,
                columns=columns, include_precursor_info=include_precursor_info,
                include_ion_mobility=include_ion_mobility)
        if data in ('chromatograms', 'both'):
            result['chromatograms'] = chromatograms_to_arrow(
                self, format=format, min_rt=min_rt, max_rt=max_rt, columns=columns)
        if data == 'both':
            return result
        elif data == 'spectra':
            return result['spectra']
        else:
            return result['chromatograms']

    # Fall back to pure Python implementation
    result = {}

    if data in ('spectra', 'both'):
        result['spectra'] = _build_spectra_arrow(
            self, format, columns, ms_levels,
            min_rt, max_rt, min_mz, max_mz,
            include_precursor_info, include_ion_mobility, pa
        )

    if data in ('chromatograms', 'both'):
        result['chromatograms'] = _build_chrom_arrow(
            self, format, columns, min_rt, max_rt, pa
        )

    if data == 'both':
        return result
    elif data == 'spectra':
        return result['spectra']
    else:
        return result['chromatograms']


def _build_spectra_arrow(exp, format, columns, ms_levels, min_rt, max_rt,
                         min_mz, max_mz, include_precursor_info, include_ion_mobility, pa):
    """Build Arrow table for spectra."""
    if format == 'long':
        all_mz, all_int, all_rt = [], [], []
        all_idx, all_ms, all_nid = [], [], []
        all_pmz, all_pch, all_pint = [], [], []
        all_ilo, all_iup, all_im = [], [], []

        for spec_idx, spec in enumerate(exp):
            ms_level = spec.getMSLevel()
            if ms_level not in ms_levels:
                continue
            rt = spec.getRT()
            if rt < min_rt or rt > max_rt:
                continue

            mzs, intensities = spec.get_peaks()
            if len(mzs) > 0:
                mask = (mzs >= min_mz) & (mzs <= max_mz)
                mzs = mzs[mask]
                intensities = intensities[mask]

            n = len(mzs)
            if n == 0:
                continue

            all_mz.append(mzs)
            all_int.append(intensities)
            all_rt.append(np.full(n, rt, dtype=np.float32))
            all_idx.append(np.full(n, spec_idx, dtype=np.uint32))
            all_ms.append(np.full(n, ms_level, dtype=np.uint8))
            all_nid.append(np.full(n, spec.getNativeID(), dtype=object))

            if include_precursor_info:
                precs = spec.getPrecursors()
                if precs:
                    p = precs[0]
                    all_pmz.extend([p.getMZ()] * n)
                    all_pch.extend([p.getCharge()] * n)
                    all_pint.extend([p.getIntensity()] * n)
                    all_ilo.extend([p.getIsolationWindowLowerOffset()] * n)
                    all_iup.extend([p.getIsolationWindowUpperOffset()] * n)
                else:
                    all_pmz.extend([None] * n)
                    all_pch.extend([None] * n)
                    all_pint.extend([None] * n)
                    all_ilo.extend([None] * n)
                    all_iup.extend([None] * n)

            if include_ion_mobility:
                if spec.containsIMData():
                    im_idx, _ = spec.getIMData()
                    fda = spec.getFloatDataArrays()[im_idx]
                    im_arr = np.asarray(fda.get_data(), dtype=np.float32)
                    if mask is not None and len(im_arr) == len(mask):
                        im_arr = im_arr[mask]
                    all_im.append(im_arr[:n])
                else:
                    all_im.append(np.full(n, np.nan, dtype=np.float32))

        d = {}
        d['mz'] = np.concatenate(all_mz) if all_mz else np.array([], dtype=np.float64)
        d['intensity'] = np.concatenate(all_int) if all_int else np.array([], dtype=np.float32)
        d['rt'] = np.concatenate(all_rt) if all_rt else np.array([], dtype=np.float32)
        if include_ion_mobility:
            d['ion_mobility'] = np.concatenate(all_im) if all_im else np.array([], dtype=np.float32)
        d['spectrum_index'] = np.concatenate(all_idx) if all_idx else np.array([], dtype=np.uint32)
        d['ms_level'] = np.concatenate(all_ms) if all_ms else np.array([], dtype=np.uint8)
        d['native_id'] = np.concatenate(all_nid) if all_nid else np.array([], dtype=object)
        if include_precursor_info:
            d['precursor_mz'] = pa.array(all_pmz, type=pa.float64())
            d['precursor_charge'] = pa.array(all_pch, type=pa.int16())
            d['precursor_intensity'] = pa.array(all_pint, type=pa.float32())
            d['isolation_lower'] = pa.array(all_ilo, type=pa.float64())
            d['isolation_upper'] = pa.array(all_iup, type=pa.float64())
    else:
        # Semi-wide
        all_mz, all_int = [], []
        all_rt, all_idx, all_ms, all_nid = [], [], [], []
        all_pmz, all_pch, all_pint = [], [], []
        all_ilo, all_iup, all_im = [], [], []

        for spec_idx, spec in enumerate(exp):
            ms_level = spec.getMSLevel()
            if ms_level not in ms_levels:
                continue
            rt = spec.getRT()
            if rt < min_rt or rt > max_rt:
                continue

            mzs, intensities = spec.get_peaks()
            if len(mzs) > 0:
                mask = (mzs >= min_mz) & (mzs <= max_mz)
                mzs = mzs[mask]
                intensities = intensities[mask]

            all_mz.append(mzs.tolist())
            all_int.append(intensities.tolist())
            all_rt.append(rt)
            all_idx.append(spec_idx)
            all_ms.append(ms_level)
            all_nid.append(spec.getNativeID())

            if include_precursor_info:
                precs = spec.getPrecursors()
                if precs:
                    p = precs[0]
                    all_pmz.append(p.getMZ())
                    all_pch.append(p.getCharge())
                    all_pint.append(p.getIntensity())
                    all_ilo.append(p.getIsolationWindowLowerOffset())
                    all_iup.append(p.getIsolationWindowUpperOffset())
                else:
                    all_pmz.append(None)
                    all_pch.append(None)
                    all_pint.append(None)
                    all_ilo.append(None)
                    all_iup.append(None)

            if include_ion_mobility:
                if spec.containsIMData():
                    im_idx, _ = spec.getIMData()
                    fda = spec.getFloatDataArrays()[im_idx]
                    all_im.append(fda.get_data().tolist())
                else:
                    all_im.append(None)

        d = {}
        d['spectrum_index'] = np.array(all_idx, dtype=np.uint32)
        d['rt'] = np.array(all_rt, dtype=np.float32)
        d['ms_level'] = np.array(all_ms, dtype=np.uint8)
        d['native_id'] = all_nid
        d['mz'] = all_mz
        d['intensity'] = all_int
        if include_ion_mobility:
            d['ion_mobility'] = all_im
        if include_precursor_info:
            d['precursor_mz'] = pa.array(all_pmz, type=pa.float64())
            d['precursor_charge'] = pa.array(all_pch, type=pa.int16())
            d['precursor_intensity'] = pa.array(all_pint, type=pa.float32())
            d['isolation_lower'] = pa.array(all_ilo, type=pa.float64())
            d['isolation_upper'] = pa.array(all_iup, type=pa.float64())

    if columns is not None:
        d = {k: v for k, v in d.items() if k in columns}
    return pa.Table.from_pydict(d)


def _build_chrom_arrow(exp, format, columns, min_rt, max_rt, pa):
    """Build Arrow table for chromatograms."""
    if format == 'long':
        all_rt, all_int, all_idx, all_nid = [], [], [], []
        all_pmz, all_prodmz = [], []

        for ci in range(exp.getNrChromatograms()):
            chrom = exp.getChromatogram(ci)
            rts, intensities = chrom.get_peaks()
            rts = np.asarray(rts, dtype=np.float64)
            intensities = np.asarray(intensities, dtype=np.float32)
            if len(rts) > 0:
                mask = (rts >= min_rt) & (rts <= max_rt)
                rts = rts[mask]
                intensities = intensities[mask]
            n = len(rts)
            if n == 0:
                continue
            all_rt.append(rts)
            all_int.append(intensities)
            all_idx.append(np.full(n, ci, dtype=np.uint32))
            all_nid.append(np.full(n, chrom.getNativeID(), dtype=object))
            all_pmz.append(np.full(n, chrom.getPrecursor().getMZ(), dtype=np.float64))
            all_prodmz.append(np.full(n, chrom.getProduct().getMZ(), dtype=np.float64))

        d = {
            'rt': np.concatenate(all_rt) if all_rt else np.array([], dtype=np.float64),
            'intensity': np.concatenate(all_int) if all_int else np.array([], dtype=np.float32),
            'chromatogram_index': np.concatenate(all_idx) if all_idx else np.array([], dtype=np.uint32),
            'native_id': np.concatenate(all_nid) if all_nid else np.array([], dtype=object),
            'precursor_mz': np.concatenate(all_pmz) if all_pmz else np.array([], dtype=np.float64),
            'product_mz': np.concatenate(all_prodmz) if all_prodmz else np.array([], dtype=np.float64),
        }
    else:
        all_rt, all_int, all_idx, all_nid = [], [], [], []
        all_pmz, all_prodmz = [], []

        for ci in range(exp.getNrChromatograms()):
            chrom = exp.getChromatogram(ci)
            rts, intensities = chrom.get_peaks()
            rts = np.asarray(rts, dtype=np.float64)
            intensities = np.asarray(intensities, dtype=np.float32)
            if len(rts) > 0:
                mask = (rts >= min_rt) & (rts <= max_rt)
                rts = rts[mask]
                intensities = intensities[mask]
            all_rt.append(rts.tolist())
            all_int.append(intensities.tolist())
            all_idx.append(ci)
            all_nid.append(chrom.getNativeID())
            all_pmz.append(chrom.getPrecursor().getMZ())
            all_prodmz.append(chrom.getProduct().getMZ())

        d = {
            'chromatogram_index': np.array(all_idx, dtype=np.uint32),
            'native_id': all_nid,
            'rt': all_rt,
            'intensity': all_int,
            'precursor_mz': np.array(all_pmz, dtype=np.float64),
            'product_mz': np.array(all_prodmz, dtype=np.float64),
        }

    if columns is not None:
        d = {k: v for k, v in d.items() if k in columns}
    return pa.Table.from_pydict(d)


@addon("MSExperiment")
def get_massql_df(self, ion_mobility=False):
    """Exports MS1 and MS2 data to pandas DataFrames for MassQL queries.

    Returns (ms1_df, ms2_df) tuple of DataFrames.
    """
    import pandas as pd
    import pyopenms

    self.updateRanges()

    def _get_polarity(spec):
        polarity = spec.getInstrumentSettings().getPolarity()
        if polarity == pyopenms.IonSource.Polarity.POLNULL:
            return 0
        elif polarity == pyopenms.IonSource.Polarity.POSITIVE:
            return 1
        elif polarity == pyopenms.IonSource.Polarity.NEGATIVE:
            return 2
        return 0

    def _get_spec_arrays(mslevel):
        for scan_num, spec in enumerate(self):
            if spec.getMSLevel() == mslevel:
                mz, inty = spec.get_peaks()
                mz = np.asarray(mz, dtype=np.float64)
                inty = np.asarray(inty, dtype=np.float32)
                max_inty = np.amax(inty, initial=0)
                sum_inty = np.sum(inty)
                i_norm = np.zeros_like(inty) if max_inty == 0 else inty / max_inty
                i_tic_norm = np.zeros_like(inty) if sum_inty == 0 else inty / sum_inty
                data = (inty, i_norm, i_tic_norm, mz, scan_num + 1, spec.getRT() / 60, _get_polarity(spec))
                cols = 7
                if mslevel == 2:
                    cols = 10
                    if spec.getPrecursors():
                        data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num) + 1, spec.getPrecursors()[0].getCharge())
                    else:
                        data += (-1, -1, -1)
                ndarr = np.empty(shape=(spec.size(), cols))
                for i in range(cols):
                    ndarr[:, i] = data[i]
                yield ndarr

    def _get_ion_spec_arrays(mslevel):
        for scan_num, spec in enumerate(self):
            if spec.getMSLevel() == mslevel:
                mz, inty = spec.get_peaks()
                mz = np.asarray(mz, dtype=np.float64)
                inty = np.asarray(inty, dtype=np.float32)
                ion_array_idx, ion_unit = spec.getIMData()
                ion_data_arr = spec.getFloatDataArrays()[ion_array_idx]
                ion_data = np.asarray(ion_data_arr.get_data(), dtype=np.float32)
                max_inty = np.amax(inty, initial=0)
                sum_inty = np.sum(inty)
                i_norm = np.zeros_like(inty) if max_inty == 0 else inty / max_inty
                i_tic_norm = np.zeros_like(inty) if sum_inty == 0 else inty / sum_inty
                data = (inty, i_norm, i_tic_norm, mz, scan_num + 1, spec.getRT() / 60, _get_polarity(spec), ion_data)
                cols = 8
                if mslevel == 2:
                    cols = 11
                    if spec.getPrecursors():
                        data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num) + 1, spec.getPrecursors()[0].getCharge())
                    else:
                        data += (-1, -1, -1)
                ndarr = np.empty(shape=(spec.size(), cols))
                for i in range(cols):
                    ndarr[:, i] = data[i]
                yield ndarr

    dtypes = {'i': 'float32', 'i_norm': 'float32', 'i_tic_norm': 'float32', 'mz': 'float64', 'scan': 'int32', 'rt': 'float32', 'polarity': 'int32'}
    if ion_mobility:
        dtypes = dict(dtypes, **{"ion_mobility": "float32"})

    if 1 in self.getMSLevels():
        spec_arrays = _get_spec_arrays(1) if not ion_mobility else _get_ion_spec_arrays(1)
        arrs = list(spec_arrays)
        if arrs:
            ms1_df = pd.DataFrame(np.concatenate(arrs, axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms1_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
    else:
        ms1_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

    dtypes = dict(dtypes, **{'precmz': 'float64', 'ms1scan': 'int32', 'charge': 'int32'})
    if 2 in self.getMSLevels():
        spec_arrays = _get_spec_arrays(2) if not ion_mobility else _get_ion_spec_arrays(2)
        arrs = list(spec_arrays)
        if arrs:
            ms2_df = pd.DataFrame(np.concatenate(arrs, axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms2_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
    else:
        ms2_df = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

    return ms1_df, ms2_df


@addon("MSExperiment")
def get_df(self, *args, **kwargs):
    """Deprecated: use to_df() instead."""
    warnings.warn("get_df() is deprecated. Use to_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_df(*args, **kwargs)


@addon("MSExperiment")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)


@addon("MSExperiment")
def get_ion_df(self):
    """Returns a DataFrame with RT, mz, intensity, and ion mobility columns for MS1 spectra."""
    import pandas as pd
    all_data = []
    for spec in self:
        if spec.getMSLevel() != 1:
            continue
        mz, inty = spec.get_peaks()
        n = len(mz)
        if n == 0:
            continue
        try:
            im_idx, im_unit = spec.getIMData()
        except RuntimeError:
            continue
        fdas = spec.getFloatDataArrays()
        if im_idx >= len(fdas):
            continue
        im_data = np.asarray(fdas[im_idx].get_data())
        if len(im_data) != n:
            continue
        rt_arr = np.full(n, spec.getRT())
        data = np.column_stack([rt_arr, mz, inty, im_data])
        all_data.append(data)
    if not all_data:
        return pd.DataFrame(columns=['RT', 'mz', 'inty', 'IM'])
    result = np.vstack(all_data)
    return pd.DataFrame(result, columns=['RT', 'mz', 'inty', 'IM'])
