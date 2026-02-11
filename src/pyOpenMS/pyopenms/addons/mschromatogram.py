"""Addon methods for MSChromatogram class."""

from __future__ import annotations
import warnings
import numpy as np
from . import addon


@addon("MSChromatogram")
def df_columns(self, columns='default', export_meta_values=True):
    """Returns a list of column names that to_df() would produce."""
    cols = ['rt', 'intensity', 'precursor_mz', 'precursor_charge',
            'product_mz', 'native_id']
    if columns == 'all':
        cols.extend(['chromatogram_type', 'comment'])
    if export_meta_values:
        mvs = []
        self.getKeys(mvs)
        for k in mvs:
            k_str = k.decode() if isinstance(k, bytes) else k
            cols.append(k_str)
    return cols


@addon("MSChromatogram")
def get_data_dict(self, columns=None, export_meta_values=True):
    """Returns a dictionary of NumPy arrays with RT, intensities, and metadata."""
    rts, intensities = self.get_peaks()
    cnt = len(rts)

    if columns is not None:
        requested = set(columns)
    else:
        requested = None

    def want(col):
        return requested is None or col in requested

    def want_explicit(col):
        return requested is not None and col in requested

    data_dict = {}

    if want('rt'):
        data_dict['rt'] = rts
    if want('intensity'):
        data_dict['intensity'] = intensities
    if want('precursor_mz'):
        data_dict['precursor_mz'] = np.full(cnt, self.getPrecursor().getMZ(), dtype=np.float64)
    if want('precursor_charge'):
        data_dict['precursor_charge'] = np.full(cnt, self.getPrecursor().getCharge(), dtype=np.uint16)
    if want('product_mz'):
        data_dict['product_mz'] = np.full(cnt, self.getProduct().getMZ(), dtype=np.float64)
    if want('native_id'):
        data_dict['native_id'] = np.full(cnt, self.getNativeID(), dtype='U100')

    if want_explicit('chromatogram_type'):
        chrom_type = self.getChromatogramType()
        type_names = {
            0: 'MASS_CHROMATOGRAM', 1: 'TOTAL_ION_CURRENT_CHROMATOGRAM',
            2: 'SELECTED_ION_CURRENT_CHROMATOGRAM', 3: 'BASEPEAK_CHROMATOGRAM',
            4: 'SELECTED_ION_MONITORING_CHROMATOGRAM',
            5: 'SELECTED_REACTION_MONITORING_CHROMATOGRAM',
            6: 'ELECTROMAGNETIC_RADIATION_CHROMATOGRAM',
            7: 'ABSORPTION_CHROMATOGRAM', 8: 'EMISSION_CHROMATOGRAM'
        }
        type_name = type_names.get(int(chrom_type), f'UNKNOWN_{chrom_type}')
        data_dict['chromatogram_type'] = np.full(cnt, type_name, dtype='U100')

    if want_explicit('comment'):
        data_dict['comment'] = np.full(cnt, self.getComment(), dtype='U100')

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

    return data_dict


@addon("MSChromatogram")
def to_df(self, columns=None, export_meta_values=True):
    """Returns a pandas DataFrame representation."""
    import pandas as pd
    return pd.DataFrame(self.get_data_dict(columns=columns, export_meta_values=export_meta_values))


@addon("MSChromatogram")
def get_df(self, *args, **kwargs):
    """Deprecated: use to_df() instead."""
    warnings.warn("get_df() is deprecated. Use to_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_df(*args, **kwargs)


@addon("MSChromatogram")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)


@addon("MSChromatogram")
def to_arrow(self, columns=None, export_meta_values=True):
    """Returns an Apache Arrow Table representation."""
    import pyarrow as pa
    return pa.Table.from_pydict(self.get_data_dict(columns=columns, export_meta_values=export_meta_values))
