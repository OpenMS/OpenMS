"""
Addon methods for MSSpectrum class.

Performance-critical methods (get_peaks, set_peaks) are implemented in C++.
This file contains pure Python helper methods including DataFrame export.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Tuple

import numpy as np

from . import addon

if TYPE_CHECKING:
    import pandas as pd


@addon("MSSpectrum")
def df_columns(self, columns='default', export_meta_values=True):
    """Returns a list of column names that to_df() would produce."""
    cols = ['mz', 'intensity', 'rt', 'ms_level', 'native_id']

    if self.containsIMData():
        cols.append('ion_mobility')
        if columns == 'all':
            cols.append('ion_mobility_unit')

    if len(self.getPrecursors()) > 0:
        cols.extend(['precursor_mz', 'precursor_charge'])

    for sda in self.getStringDataArrays():
        name = sda.getName()
        if name == 'IonNames':
            cols.append('ion_annotation')
            if columns == 'all':
                cols.append(f'string_array:{name}')
        elif columns == 'all':
            cols.append(f'string_array:{name}')

    if columns == 'all':
        for fda in self.getFloatDataArrays():
            cols.append(f'float_array:{fda.getName()}')
        for ida in self.getIntegerDataArrays():
            cols.append(f'int_array:{ida.getName()}')

    if export_meta_values:
        mvs = []
        self.getKeys(mvs)
        for k in mvs:
            k_str = k.decode() if isinstance(k, bytes) else k
            cols.append(k_str)

    return cols


@addon("MSSpectrum")
def get_data_dict(self, columns=None, export_meta_values=True):
    """Returns a dictionary of NumPy arrays with m/z, intensities, and metadata."""
    mzs, intensities = self.get_peaks()
    cnt = len(mzs)

    if columns is not None:
        requested = set(columns)
    else:
        requested = None

    def want(col):
        return requested is None or col in requested

    data_dict = {}

    if want('mz'):
        data_dict['mz'] = mzs
    if want('intensity'):
        data_dict['intensity'] = intensities

    if want('rt'):
        data_dict['rt'] = np.full(cnt, self.getRT(), dtype=np.float64)
    if want('ms_level'):
        data_dict['ms_level'] = np.full(cnt, self.getMSLevel(), dtype=np.uint16)
    if want('native_id'):
        data_dict['native_id'] = np.full(cnt, self.getNativeID(), dtype='U100')

    if want('ion_mobility') or want('ion_mobility_unit'):
        if self.containsIMData():
            im_index, drift_time_unit = self.getIMData()
            im_arrays = self.getFloatDataArrays()
            if want('ion_mobility'):
                if 0 <= im_index < len(im_arrays):
                    data_dict['ion_mobility'] = np.asarray(
                        im_arrays[im_index].get_data(), dtype=np.float64
                    )
                else:
                    data_dict['ion_mobility'] = np.full(cnt, np.nan, dtype=np.float64)
            if want('ion_mobility_unit'):
                data_dict['ion_mobility_unit'] = np.full(
                    cnt, self.getDriftTimeUnitAsString(), dtype='U50'
                )
        else:
            if requested is not None:
                if want('ion_mobility'):
                    data_dict['ion_mobility'] = np.full(cnt, np.nan, dtype=np.float64)
                if want('ion_mobility_unit'):
                    data_dict['ion_mobility_unit'] = np.full(cnt, '', dtype='U1')

    if want('precursor_mz') or want('precursor_charge'):
        precursors = self.getPrecursors()
        if len(precursors) > 0:
            precursor = precursors[0]
            if want('precursor_mz'):
                data_dict['precursor_mz'] = np.full(cnt, precursor.getMZ(), dtype=np.float64)
            if want('precursor_charge'):
                data_dict['precursor_charge'] = np.full(cnt, precursor.getCharge(), dtype=np.int16)
        else:
            if requested is not None:
                if want('precursor_mz'):
                    data_dict['precursor_mz'] = np.full(cnt, np.nan, dtype=np.float64)
                if want('precursor_charge'):
                    data_dict['precursor_charge'] = np.full(cnt, 0, dtype=np.int16)

    if want('ion_annotation'):
        ion_annotations = np.full(cnt, '', dtype='U1')
        for sda in self.getStringDataArrays():
            if sda.getName() == 'IonNames':
                if len(sda) == cnt:
                    annotations = sda.get_data()
                    max_len = max((len(s) for s in annotations), default=1)
                    ion_annotations = np.array(annotations, dtype=f'U{max_len}')
                break
        if requested is not None or any(ion_annotations != ''):
            data_dict['ion_annotation'] = ion_annotations

    # Meta values
    if requested is None and export_meta_values:
        mvs = []
        self.getKeys(mvs)
        for k in mvs:
            if not self.metaValueExists(k):
                continue
            v = self.getMetaValue(k)
            k_str = k.decode() if isinstance(k, bytes) else k
            try:
                if type(v) is type(True):
                    data_dict[k_str] = np.full(cnt, v, dtype=np.bool_)
                elif isinstance(v, int):
                    data_dict[k_str] = np.full(cnt, v, dtype=np.int64)
                elif isinstance(v, float):
                    data_dict[k_str] = np.full(cnt, v, dtype=np.float64)
                elif isinstance(v, str):
                    data_dict[k_str] = np.full(cnt, v, dtype=f"U{max(len(v), 1)}")
                else:
                    data_dict[k_str] = np.full(cnt, str(v), dtype='object')
            except Exception:
                data_dict[k_str] = np.full(cnt, str(v), dtype='object')
    elif requested is not None:
        mvs = []
        self.getKeys(mvs)
        mv_names = {(k.decode() if isinstance(k, bytes) else k): k for k in mvs}
        for col in requested:
            if col in mv_names:
                k = mv_names[col]
                if self.metaValueExists(k):
                    v = self.getMetaValue(k)
                    try:
                        if type(v) is type(True):
                            data_dict[col] = np.full(cnt, v, dtype=np.bool_)
                        elif isinstance(v, int):
                            data_dict[col] = np.full(cnt, v, dtype=np.int64)
                        elif isinstance(v, float):
                            data_dict[col] = np.full(cnt, v, dtype=np.float64)
                        elif isinstance(v, str):
                            data_dict[col] = np.full(cnt, v, dtype=f"U{max(len(v), 1)}")
                        else:
                            data_dict[col] = np.full(cnt, str(v), dtype='object')
                    except Exception:
                        data_dict[col] = np.full(cnt, str(v), dtype='object')

    # Custom data arrays - only when explicitly requested
    if requested is not None:
        for fda in self.getFloatDataArrays():
            col_name = f'float_array:{fda.getName()}'
            if col_name in requested:
                if len(fda) == cnt:
                    data_dict[col_name] = np.asarray(fda.get_data(), dtype=np.float32)
                else:
                    data_dict[col_name] = np.full(cnt, np.nan, dtype=np.float32)

        for ida in self.getIntegerDataArrays():
            col_name = f'int_array:{ida.getName()}'
            if col_name in requested:
                if len(ida) == cnt:
                    data_dict[col_name] = np.asarray(ida.get_data(), dtype=np.int64)
                else:
                    data_dict[col_name] = np.full(cnt, 0, dtype=np.int64)

        for sda in self.getStringDataArrays():
            col_name = f'string_array:{sda.getName()}'
            if col_name in requested:
                if len(sda) == cnt:
                    strings = sda.get_data()
                    max_len = max((len(s) for s in strings), default=1)
                    data_dict[col_name] = np.array(strings, dtype=f'U{max_len}')
                else:
                    data_dict[col_name] = np.full(cnt, '', dtype='U1')

    return data_dict


@addon("MSSpectrum")
def to_df(self, columns=None, export_meta_values=True):
    """Returns a pandas DataFrame representation of the MSSpectrum."""
    import pandas as pd
    data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
    return pd.DataFrame(data_dict)


@addon("MSSpectrum")
def to_arrow(self, columns=None, export_meta_values=True):
    """Returns an Apache Arrow Table representation of the MSSpectrum."""
    import pyarrow as pa
    data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
    return pa.Table.from_pydict(data_dict)


@addon("MSSpectrum")
def get_tic(self) -> float:
    """Get total ion current (sum of all intensities)."""
    return self.calculateTIC()


@addon("MSSpectrum")
def get_df(self, *args, **kwargs):
    """Deprecated: use to_df() instead."""
    warnings.warn("get_df() is deprecated. Use to_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_df(*args, **kwargs)


@addon("MSSpectrum")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)


@addon("MSSpectrum")
def intensityInRange(self, mz_min: float, mz_max: float) -> float:
    """Returns the sum of intensities of all peaks within the given m/z range.

    Parameters
    ----------
    mz_min : float
        Minimum m/z value (inclusive)
    mz_max : float
        Maximum m/z value (inclusive)

    Returns
    -------
    float
        Sum of intensities of peaks with mz_min <= m/z <= mz_max
    """
    if len(self) == 0:
        return 0.0
    mz, intensity = self.get_peaks()
    mask = (mz >= mz_min) & (mz <= mz_max)
    return float(np.sum(intensity[mask]))


@addon("MSSpectrum")
def get_base_peak(self) -> Tuple[float, float]:
    """Get the base peak (highest intensity peak)."""
    if len(self) == 0:
        return (0.0, 0.0)
    mz, intensity = self.get_peaks()
    idx = np.argmax(intensity)
    return (mz[idx], intensity[idx])

