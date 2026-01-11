# Add Pythonic dict-style parameter access for DefaultParamHandler subclasses

## Summary

Add a Pythonic interface for accessing and modifying parameters on classes that inherit from `DefaultParamHandler`. This would enable cleaner, more intuitive parameter management.

## Current API

```python
# Current workflow - verbose and C++-like
algorithm = PeakPickerHiRes()
defaults = algorithm.getDefaults()
defaults.setValue(b"threshold", 0.5, b"description")
algorithm.setParameters(defaults)

# Or with recent Param improvements
defaults = algorithm.getDefaults()
defaults["threshold"] = 0.5
algorithm.setParameters(defaults)
```

## Proposed API

```python
# Option 1: Direct dict-style access via property
algorithm = PeakPickerHiRes()
algorithm.parameters["threshold"] = 0.5

# Option 2: Assign dict directly
algorithm.parameters = {"threshold": 0.5, "iterations": 10}

# Option 3: Round-trip workflow
d = algorithm.parameters.to_dict()
d["threshold"] = 0.75
algorithm.parameters = d

# Option 4: Read-only defaults access
print(algorithm.defaults["threshold"])  # Get default value
```

## Implementation Plan

### Architecture Challenge

`DefaultParamHandler` is a C++ base class with `#wrap-ignore` in its `.pxd` file, meaning it's not wrapped as a standalone Python class. It's inherited by **180+ algorithm classes** via the `wrap-inherits` directive. This makes it impractical to add individual addon files for each class.

### Recommended Approach: Post-import Enhancement

Add Pythonic methods at Python import time via `_python_extras.py`. This is clean, requires no autowrap changes, and works immediately for all subclasses.

**File: `src/pyOpenMS/pyopenms/_python_extras.py`**

```python
class _ParametersProxy:
    """Proxy object for dict-style parameter access on DefaultParamHandler subclasses."""

    def __init__(self, obj):
        self._obj = obj

    def __getitem__(self, key):
        return self._obj.getParameters()[key]

    def __setitem__(self, key, value):
        p = self._obj.getParameters()
        p[key] = value
        self._obj.setParameters(p)

    def __contains__(self, key):
        return key in self._obj.getParameters()

    def __iter__(self):
        return iter(self._obj.getParameters())

    def __len__(self):
        return len(self._obj.getParameters())

    def to_dict(self):
        return self._obj.getParameters().to_dict()

    def update(self, d):
        p = self._obj.getParameters()
        p.update(d)
        self._obj.setParameters(p)


def _add_param_handler_interface(cls):
    """Add Pythonic parameter interface to a DefaultParamHandler subclass."""

    @property
    def parameters(self):
        return _ParametersProxy(self)

    @parameters.setter
    def parameters(self, value):
        if isinstance(value, dict):
            p = self.getDefaults()
            p.update(value)
            self.setParameters(p)
        else:
            self.setParameters(value)

    @property
    def defaults(self):
        return self.getDefaults()

    cls.parameters = parameters
    cls.defaults = defaults
    return cls


# Apply to all DefaultParamHandler subclasses
def _enhance_param_handlers():
    import pyopenms
    for name in dir(pyopenms):
        obj = getattr(pyopenms, name)
        if isinstance(obj, type) and hasattr(obj, 'getParameters') and hasattr(obj, 'setParameters'):
            _add_param_handler_interface(obj)

_enhance_param_handlers()
```

### How It Works

1. **`_ParametersProxy`**: A proxy object that wraps parameter access, providing dict-like interface while auto-syncing changes back to the algorithm via `setParameters()`.

2. **`parameters` property**:
   - **Getter**: Returns a `_ParametersProxy` for dict-style access
   - **Setter**: Accepts a `dict` or `Param` object, applies defaults, and sets parameters

3. **`defaults` property**: Read-only access to default parameters

4. **`_enhance_param_handlers()`**: Iterates over all pyopenms classes at import time and adds the interface to any class with `getParameters`/`setParameters` methods.

## Tasks

- [ ] Implement `_ParametersProxy` class in `_python_extras.py`
- [ ] Implement `_add_param_handler_interface` function
- [ ] Add `_enhance_param_handlers()` call at module load
- [ ] Add unit tests for new interface
- [ ] Update documentation

## Testing

```python
def test_param_handler_pythonic_interface():
    pp = pyopenms.PeakPickerHiRes()

    # Test property access
    assert hasattr(pp, 'parameters')
    assert hasattr(pp, 'defaults')

    # Test dict-style read
    assert "signal_to_noise" in pp.parameters

    # Test dict-style write
    pp.parameters["signal_to_noise"] = 2.0
    assert pp.getParameters()["signal_to_noise"] == 2.0

    # Test dict assignment
    pp.parameters = {"signal_to_noise": 3.0}
    assert pp.parameters["signal_to_noise"] == 3.0

    # Test to_dict round-trip
    d = pp.parameters.to_dict()
    d["signal_to_noise"] = 4.0
    pp.parameters = d
    assert pp.parameters["signal_to_noise"] == 4.0

    # Test defaults (read-only)
    default_val = pp.defaults["signal_to_noise"]
    assert default_val is not None
```

## Benefits

| Feature | Before | After |
|---------|--------|-------|
| Set parameter | `p = algo.getDefaults(); p["k"] = v; algo.setParameters(p)` | `algo.parameters["k"] = v` |
| Set multiple | Same multi-step | `algo.parameters = {"k1": v1, "k2": v2}` |
| Read parameter | `algo.getParameters()["key"]` | `algo.parameters["key"]` |
| Check existence | `"key" in algo.getParameters()` | `"key" in algo.parameters` |
| Round-trip | Manual bytes/string handling | `d = algo.parameters.to_dict(); d["k"] = v; algo.parameters = d` |

## Related

- Recent improvements to `Param` class: `to_dict()`, `from_dict()`, string key support
- `src/pyOpenMS/addons/Param.pyx` - Pythonic Param interface

## Future Enhancements

- Modify autowrap to support addon inheritance for `wrap-inherits` classes
- Add `configure(**kwargs)` method for keyword-based configuration
- Add context manager for temporary parameter changes
