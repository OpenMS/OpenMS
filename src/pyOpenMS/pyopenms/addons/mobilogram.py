"""Mobilogram addon methods for DataFrame support."""
import numpy as np
from . import addon


@addon("Mobilogram")
def df_columns(self, columns='default'):
    """Returns a list of column names that to_df() would produce."""
    return ['mobility', 'intensity', 'rt', 'drift_time_unit']


@addon("Mobilogram")
def get_data_dict(self, columns=None):
    """Returns a dictionary of NumPy arrays with mobility, intensities, and metadata."""
    mobilities, intensities = self.get_peaks()
    cnt = len(mobilities)

    if columns is not None:
        requested = set(columns)
    else:
        requested = None

    def want(col):
        return requested is None or col in requested

    data_dict = {}

    if want('mobility'):
        data_dict['mobility'] = mobilities
    if want('intensity'):
        data_dict['intensity'] = intensities
    if want('rt'):
        data_dict['rt'] = np.full(cnt, self.getRT(), dtype=np.float64)
    if want('drift_time_unit'):
        unit_str = self.getDriftTimeUnitAsString()
        unit_decoded = unit_str.decode('utf-8') if isinstance(unit_str, bytes) else str(unit_str)
        data_dict['drift_time_unit'] = np.full(cnt, unit_decoded, dtype='U50')

    return data_dict


@addon("Mobilogram")
def to_df(self, columns=None):
    """Returns a pandas DataFrame representation of the Mobilogram."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for to_df(). Install with: pip install pandas")
    data_dict = self.get_data_dict(columns=columns)
    return pd.DataFrame(data_dict)


@addon("Mobilogram")
def to_arrow(self, columns=None):
    """Returns an Apache Arrow Table representation of the Mobilogram."""
    try:
        import pyarrow as pa
    except ImportError:
        raise ImportError("pyarrow is required for to_arrow(). Install with: pip install pyarrow")
    data_dict = self.get_data_dict(columns=columns)
    return pa.Table.from_pydict(data_dict)


@addon("Mobilogram")
def get_df(self, *args, **kwargs):
    """Deprecated: Use to_df() instead."""
    import warnings
    warnings.warn(
        "get_df() is deprecated. Use to_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.to_df(*args, **kwargs)


@addon("Mobilogram")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    import warnings
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)
    
@addon("Mobilogram")
def get_peaks_struct(self):
    """
    Returns a zero-copy numpy structured array of the mobilogram's peaks (AoS layout).
    """
    # Get the raw 1D byte array from C++
    raw_view = self.get_peaks_struct_mv()

    peak_dtype = np.dtype({
        "names": ["mobility", "intensity"],
        "formats": [np.float64, np.float32],
        "offsets": [0, 8],
        "itemsize": 16
    })

    # Cast the byte array directly into our structured array
    return raw_view.view(peak_dtype)
