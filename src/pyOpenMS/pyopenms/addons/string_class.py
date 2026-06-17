"""
Pure Python String wrapper for pyOpenMS compatibility.

OpenMS::String has a nanobind type caster that auto-converts to Python str,
so we cannot also register nb::class_<OpenMS::String>. This pure Python class
provides the pyopenms.String(...) API that existing user code expects.
"""


class String:
    """Wrapper around OpenMS::String, stored as a Python str internally."""

    __slots__ = ("_value",)

    def __init__(self, value=""):
        if isinstance(value, String):
            self._value = value._value
        elif isinstance(value, bytes):
            self._value = value.decode("utf-8")
        else:
            self._value = str(value)

    def toString(self):
        """Return the string value."""
        return self._value

    def c_str(self):
        """Return the string as bytes (UTF-8 if stored as str)."""
        if isinstance(self._value, bytes):
            return self._value
        return self._value.encode("utf-8")

    def __str__(self):
        if isinstance(self._value, bytes):
            return self._value.decode("utf-8", errors="replace")
        return self._value

    def __repr__(self):
        return f"String({self._value!r})"

    def __len__(self):
        return len(self._value)

    def __eq__(self, other):
        if isinstance(other, String):
            return self._value == other._value
        if isinstance(other, str):
            return self._value == other
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        return hash(self._value)

    def __bytes__(self):
        return self.c_str()

    def __bool__(self):
        return len(self._value) > 0
