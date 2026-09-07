"""Deprecated method aliases kept for one release after a rename.

The ``_mv`` -> ``_view`` renames have their own module, ``deprecated_mv_aliases``; this one is
for individual renames that do not follow that pattern.
"""
import warnings
from . import addon


def _deprecated_alias(cls_name, old_name, new_name, reason=""):
    """Create a deprecated alias that delegates to the current method."""
    def method(self, *args, **kwargs):
        warnings.warn(
            f"{old_name}() is deprecated. Use {new_name}() instead.{reason}",
            DeprecationWarning, stacklevel=2
        )
        return getattr(self, new_name)(*args, **kwargs)
    method.__name__ = old_name
    method.__qualname__ = f"{cls_name}.{old_name}"
    addon(cls_name)(method)


# MSSpectrum: get_drift_time_unit() was ambiguous against the drift_time_unit scalar property,
# and used to report that scalar while testing for the per-peak array -- so it answered NONE for
# a frame whose array had a perfectly good unit. The replacement is array-derived throughout.
_deprecated_alias(
    "MSSpectrum", "get_drift_time_unit", "get_drift_time_array_unit",
    " Note the new method reports the unit of the per-peak ion mobility array;"
    " for the spectrum-wide scalar use getDriftTimeUnit() or the drift_time_unit property.",
)
