"""
Param dict-like API addon.

Adds keys(), values(), items(), asDict(), to_dict(), get(), __getitem__,
__setitem__, __contains__, __iter__, update(), and from_dict() to the
Param class for backward compatibility with the Cython bindings.
"""

from . import addon


def _normalize_key(key):
    """Normalize a bytes or str key to str."""
    return key.decode() if isinstance(key, bytes) else key


def _get_key_value_pairs(param):
    """Return list of (str_key, value) pairs using full hierarchical paths."""
    all_keys = param._get_all_keys()
    return [(k, param.getValue(k)) for k in all_keys]


@addon("Param")
def keys(self):
    """Return list of parameter keys as str."""
    return list(self._get_all_keys())


@addon("Param")
def values(self):
    """Return list of parameter values."""
    return [self.getValue(k) for k in self._get_all_keys()]


@addon("Param")
def items(self):
    """Return list of (str_key, value) tuples."""
    return [(k, self.getValue(k)) for k in self._get_all_keys()]


@addon("Param")
def asDict(self):
    """Return dict with str keys."""
    return {k: self.getValue(k) for k in self._get_all_keys()}


@addon("Param")
def to_dict(self):
    """Return dict with string keys."""
    return {k: self.getValue(k) for k in self._get_all_keys()}


@addon("Param")
def get(self, key, default=None):
    """Get a parameter value by key, returning default if not found."""
    key = _normalize_key(key)
    if self.exists(key):
        return self.getValue(key)
    return default


@addon("Param", "__getitem__")
def _getitem(self, key):
    """Get a parameter value by key, raising KeyError if not found."""
    key = _normalize_key(key)
    if not self.exists(key):
        raise KeyError(key)
    return self.getValue(key)


@addon("Param", "__setitem__")
def _setitem(self, key, value):
    """Set a parameter value by key."""
    key = _normalize_key(key)
    self.setValue(key, value)


@addon("Param", "__contains__")
def _contains(self, key):
    """Check if a parameter key exists."""
    key = _normalize_key(key)
    return self.exists(key)


@addon("Param", "__iter__")
def _iter(self):
    """Iterate over parameter keys as str."""
    yield from self._get_all_keys()


@addon("Param")
def update(self, source, flag=None):
    """Update parameters from a Param or dict.

    If source is a Param and flag is truthy (e.g., 1 or True),
    only update keys that already exist in self.
    If source is a dict, set all key-value pairs.
    """
    # Import here to avoid circular imports at module level
    import pyopenms

    if isinstance(source, pyopenms.Param):
        for k in source._get_all_keys():
            if flag and not self.exists(k):
                continue
            self.setValue(k, source.getValue(k))
    elif isinstance(source, dict):
        for k, v in source.items():
            k = _normalize_key(k)
            self.setValue(k, v)
    else:
        raise TypeError(f"expected Param or dict, got {type(source).__name__}")


def _from_dict(d):
    """Create a Param from a dict."""
    import pyopenms
    p = pyopenms.Param()
    for k, v in d.items():
        k = _normalize_key(k)
        p.setValue(k, v)
    return p


# Register from_dict as a staticmethod via the addon registry directly
# (staticmethod works reliably on nanobind types; classmethod may not)
from . import _addon_registry
if "Param" not in _addon_registry:
    _addon_registry["Param"] = {}
_addon_registry["Param"]["from_dict"] = staticmethod(_from_dict)
