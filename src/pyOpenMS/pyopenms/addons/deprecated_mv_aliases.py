"""Deprecated _mv aliases for zero-copy methods.

The _mv suffix (memory view) has been renamed to _view for clarity.
These aliases will be removed in a future release.
"""
import warnings
from . import addon


def _deprecated_mv(cls_name, old_name, new_name):
    """Create a deprecated _mv alias that delegates to the _view method."""
    def method(self, *args, **kwargs):
        warnings.warn(
            f"{old_name}() is deprecated. Use {new_name}() instead.",
            DeprecationWarning, stacklevel=2
        )
        return getattr(self, new_name)(*args, **kwargs)
    method.__name__ = old_name
    method.__qualname__ = f"{cls_name}.{old_name}"
    addon(cls_name)(method)


# FloatDataArray
_deprecated_mv("FloatDataArray", "get_data_mv", "data_view")

# IntegerDataArray
_deprecated_mv("IntegerDataArray", "get_data_mv", "data_view")

# MatrixDouble
_deprecated_mv("MatrixDouble", "get_matrix_mv", "matrix_view")
_deprecated_mv("MatrixDouble", "get_matrix_as_view", "matrix_view")

# MSSpectrum
_deprecated_mv("MSSpectrum", "get_drift_time_array_mv", "drift_time_array_view")

# OSBinaryDataArray
_deprecated_mv("OSBinaryDataArray", "get_data_mv", "data_view")

# OSSpectrum
_deprecated_mv("OSSpectrum", "get_mz_array_mv", "mz_array_view")
_deprecated_mv("OSSpectrum", "get_intensity_array_mv", "intensity_array_view")
_deprecated_mv("OSSpectrum", "get_drift_time_array_mv", "drift_time_array_view")

# OSChromatogram
_deprecated_mv("OSChromatogram", "get_time_array_mv", "time_array_view")
_deprecated_mv("OSChromatogram", "get_intensity_array_mv", "intensity_array_view")
