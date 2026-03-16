"""
Pure Python DataValue wrapper for pyOpenMS compatibility.

OpenMS::DataValue has a nanobind type caster that auto-converts to Python
primitives, so we cannot also register nb::class_<OpenMS::DataValue>.
This pure Python class provides the pyopenms.DataValue(...) API.
"""


class DataType:
    """Enum-like class for DataValue types (module-level access)."""
    STRING_VALUE = 0
    INT_VALUE = 1
    DOUBLE_VALUE = 2
    STRING_LIST = 3
    INT_LIST = 4
    DOUBLE_LIST = 5
    EMPTY_VALUE = 6


class DataValue:
    """Wrapper around OpenMS::DataValue, stored as a Python object internally."""

    # Enum constants as class attributes
    STRING_VALUE = DataType.STRING_VALUE
    INT_VALUE = DataType.INT_VALUE
    DOUBLE_VALUE = DataType.DOUBLE_VALUE
    STRING_LIST = DataType.STRING_LIST
    INT_LIST = DataType.INT_LIST
    DOUBLE_LIST = DataType.DOUBLE_LIST
    EMPTY_VALUE = DataType.EMPTY_VALUE

    __slots__ = ("_value", "_type")

    def __init__(self, value=None):
        if isinstance(value, DataValue):
            self._value = value._value
            self._type = value._type
        elif value is None:
            self._value = None
            self._type = DataType.EMPTY_VALUE
        elif isinstance(value, bool):
            # bool before int since bool is a subclass of int
            self._value = int(value)
            self._type = DataType.INT_VALUE
        elif isinstance(value, int):
            self._value = value
            self._type = DataType.INT_VALUE
        elif isinstance(value, float):
            self._value = value
            self._type = DataType.DOUBLE_VALUE
        elif isinstance(value, (str, bytes)):
            self._value = str(value) if isinstance(value, str) else value.decode("utf-8")
            self._type = DataType.STRING_VALUE
        elif isinstance(value, (list, tuple)):
            if len(value) == 0:
                self._value = []
                self._type = DataType.STRING_LIST
            else:
                first = value[0]
                if isinstance(first, str):
                    self._value = [str(v) for v in value]
                    self._type = DataType.STRING_LIST
                elif isinstance(first, float):
                    self._value = [float(v) for v in value]
                    self._type = DataType.DOUBLE_LIST
                elif isinstance(first, int):
                    self._value = [int(v) for v in value]
                    self._type = DataType.INT_LIST
                else:
                    self._value = [str(v) for v in value]
                    self._type = DataType.STRING_LIST
        else:
            self._value = str(value)
            self._type = DataType.STRING_VALUE

    def valueType(self):
        """Return the type of the stored value."""
        return self._type

    def isEmpty(self):
        """Return True if the value is empty."""
        return self._type == DataType.EMPTY_VALUE

    def toInt(self):
        """Return value as int."""
        return int(self._value)

    def toDouble(self):
        """Return value as float."""
        return float(self._value)

    def toString(self):
        """Return value as string."""
        if self._value is None:
            return ""
        return str(self._value)

    def toBool(self):
        """Return value as bool."""
        if isinstance(self._value, str):
            return self._value.lower() in ("true", "1")
        return bool(self._value)

    def toStringList(self):
        """Return value as list of strings."""
        if isinstance(self._value, list):
            return [str(v) for v in self._value]
        return [str(self._value)] if self._value is not None else []

    def toIntList(self):
        """Return value as list of ints."""
        if isinstance(self._value, list):
            return [int(v) for v in self._value]
        return [int(self._value)] if self._value is not None else []

    def toDoubleList(self):
        """Return value as list of floats."""
        if isinstance(self._value, list):
            return [float(v) for v in self._value]
        return [float(self._value)] if self._value is not None else []

    _EPSILON = 1e-6  # Match OpenMS DataValue::operator== epsilon

    def __eq__(self, other):
        if isinstance(other, DataValue):
            if self._type != other._type:
                return False
            if self._type == DataValue.DOUBLE_VALUE:
                return abs(self._value - other._value) < DataValue._EPSILON
            return self._value == other._value
        return NotImplemented

    def __hash__(self):
        if isinstance(self._value, list):
            return hash((self._type, tuple(self._value)))
        if self._type == DataValue.DOUBLE_VALUE:
            # Round to epsilon precision so equal values (within epsilon) hash equally
            return hash((self._type, round(self._value / DataValue._EPSILON)))
        return hash((self._type, self._value))

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __repr__(self):
        return f"DataValue({self._value!r})"

    def __str__(self):
        return self.toString()

    def __bool__(self):
        return not self.isEmpty()
