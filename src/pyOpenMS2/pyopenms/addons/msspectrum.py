"""
Addon methods for MSSpectrum class.

Performance-critical methods (get_peaks, set_peaks) are implemented in C++.
This file contains pure Python helper methods.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from . import addon

if TYPE_CHECKING:
    import pandas as pd


@addon("MSSpectrum")
def __repr__(self) -> str:
    """Return string representation of MSSpectrum."""
    return (
        f"MSSpectrum(n_peaks={len(self)}, "
        f"rt={self.getRT():.2f}s, "
        f"ms_level={self.getMSLevel()})"
    )


@addon("MSSpectrum")
def __str__(self) -> str:
    """Return user-friendly string of MSSpectrum."""
    return f"MSSpectrum with {len(self)} peaks at RT={self.getRT():.2f}s"


@addon("MSSpectrum")
def get_data_dict(
    self,
    columns: Optional[Sequence[str]] = None,
    export_meta_values: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Get spectrum data as a dictionary of numpy arrays.

    This is the foundation for DataFrame export. Use `to_df()` for
    pandas DataFrame output.

    Parameters
    ----------
    columns : sequence of str, optional
        Columns to include. If None, includes mz, intensity.
    export_meta_values : bool, default True
        Whether to include meta values from float/int/string data arrays.

    Returns
    -------
    dict
        Dictionary mapping column names to numpy arrays.
    """
    # Get peaks (implemented in C++ for performance)
    mz, intensity = self.get_peaks()

    data = {
        "mz": mz,
        "intensity": intensity,
    }

    # Filter columns if specified
    if columns is not None:
        columns_set = set(columns)
        data = {k: v for k, v in data.items() if k in columns_set}

    # Add meta values from data arrays
    if export_meta_values:
        # Float data arrays
        for fda in self.getFloatDataArrays():
            name = fda.getName()
            if columns is None or name in columns:
                arr = np.array([fda[i] for i in range(len(fda))], dtype=np.float32)
                data[name] = arr

        # Integer data arrays
        for ida in self.getIntegerDataArrays():
            name = ida.getName()
            if columns is None or name in columns:
                arr = np.array([ida[i] for i in range(len(ida))], dtype=np.int32)
                data[name] = arr

        # String data arrays
        for sda in self.getStringDataArrays():
            name = sda.getName()
            if columns is None or name in columns:
                arr = np.array([sda[i] for i in range(len(sda))], dtype=object)
                data[name] = arr

    return data


@addon("MSSpectrum")
def to_df(
    self,
    columns: Optional[Sequence[str]] = None,
    export_meta_values: bool = True,
) -> "pd.DataFrame":
    """
    Export spectrum to a pandas DataFrame.

    Parameters
    ----------
    columns : sequence of str, optional
        Columns to include.
    export_meta_values : bool, default True
        Whether to include meta values.

    Returns
    -------
    pd.DataFrame
        DataFrame with spectrum data.
    """
    import pandas as pd

    data = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
    return pd.DataFrame(data)


@addon("MSSpectrum")
def get_tic(self) -> float:
    """Get total ion current (sum of all intensities)."""
    return self.calculateTIC()


@addon("MSSpectrum")
def get_base_peak(self) -> Tuple[float, float]:
    """
    Get the base peak (highest intensity peak).

    Returns
    -------
    tuple
        (mz, intensity) of the base peak, or (0, 0) if spectrum is empty.
    """
    if len(self) == 0:
        return (0.0, 0.0)

    mz, intensity = self.get_peaks()
    idx = np.argmax(intensity)
    return (mz[idx], intensity[idx])
