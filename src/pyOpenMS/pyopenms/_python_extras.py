class ParamValue:
    """Python wrapper for ParamValue providing backward compatibility.

    In nanobind, ParamValue uses a type caster for transparent C++ <-> Python
    conversion, so we cannot use nb::class_<ParamValue>. This pure-Python
    wrapper provides the same API as the old Cython ParamValue.
    """
    # Value type constants matching C++ ParamValue::ValueType
    STRING_VALUE = 0
    INT_VALUE = 1
    DOUBLE_VALUE = 2
    STRING_LIST = 3
    INT_LIST = 4
    DOUBLE_LIST = 5
    EMPTY_VALUE = 6

    def __init__(self, value=None):
        self._value = value

    def isEmpty(self):
        return self._value is None

    def valueType(self):
        if self._value is None:
            return self.EMPTY_VALUE
        elif isinstance(self._value, str):
            return self.STRING_VALUE
        elif isinstance(self._value, bool):
            return self.INT_VALUE
        elif isinstance(self._value, int):
            return self.INT_VALUE
        elif isinstance(self._value, float):
            return self.DOUBLE_VALUE
        elif isinstance(self._value, list):
            if len(self._value) == 0:
                return self.STRING_LIST
            elif isinstance(self._value[0], str):
                return self.STRING_LIST
            elif isinstance(self._value[0], int):
                return self.INT_LIST
            elif isinstance(self._value[0], float):
                return self.DOUBLE_LIST
        return self.EMPTY_VALUE

    def toString(self):
        if self._value is None:
            return ""
        return str(self._value)

    def toBool(self):
        if isinstance(self._value, bool):
            return self._value
        if isinstance(self._value, str):
            return self._value.lower() in ("true", "1", "yes")
        return bool(self._value)

    def toInt(self):
        return int(self._value)

    def toDouble(self):
        return float(self._value)

    def toStringVector(self):
        return [str(x) for x in self._value]

    def toIntVector(self):
        return [int(x) for x in self._value]

    def toDoubleVector(self):
        return [float(x) for x in self._value]

    def __repr__(self):
        return f"ParamValue({self._value!r})"

    def __eq__(self, other):
        if isinstance(other, ParamValue):
            return self._value == other._value
        return self._value == other


class SimpleOpenMSSpectraFactory:
    """A factory that returns ISpectrumAccess implementations for an MSExperiment."""

    @staticmethod
    def getSpectrumAccessOpenMSPtr(exp):
        # Deferred import: this module is loaded from pyopenms/__init__.py, so the
        # binding names are not importable at module level. (Both branches below
        # raised NameError before this import existed -- nothing ever called them.)
        from . import SpectrumAccessOpenMS, SpectrumAccessOpenMSCached

        # Scans C++-side: exp[i] and getChromatograms() hand out owned copies
        # (see OWNERSHIP.md), so probing the marker from Python would copy
        # every spectrum in the run just to read DataProcessing metadata.
        if exp._contains_cached_data_marker():
            # Same contract as the C++ factory: the metadata file the
            # experiment was loaded from must have its ".cached" side-car
            # next to it, or the constructor raises.
            return SpectrumAccessOpenMSCached(exp.getLoadedFilePath())
        return SpectrumAccessOpenMS(exp)


