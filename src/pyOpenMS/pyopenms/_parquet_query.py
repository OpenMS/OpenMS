"""Typed query builders for Parquet-backed data access.

This module defines a small abstract query interface and a concrete
ChromatogramQuery for XICParquetFile. The builders are pure Python to
keep query builders separate from the C++ bindings while still using the C++
ParquetFilter under the hood for predicate pushdown.

Why this exists:
  - Keeps a stable, Python-friendly chaining API (filter_* methods).
  - Delegates predicate construction to the C++ ParquetFilterBuilder (bound in
    `_pyopenms_format`), enabling pushdown without string-encoded filters.
"""

from abc import ABC, abstractmethod


class ParquetQuery(ABC):
    """Abstract base class for Parquet query builders.

    Shared query methods:
      - and_(), or_()
      - filter(column, value, op="=")
      - to_dict(explode=False)
      - to_df(explode=True)

    Subclasses are expected to add domain-specific filter helpers and
    implement the data access methods.
    """

    def __init__(self, typed_builder=None):
        self._pending_connector = None
        self._condition_count = 0
        self._typed = typed_builder
        self._typed_has_condition = False
        if self._typed is None:
            try:
                from . import ParquetFilterBuilder
                self._typed = ParquetFilterBuilder()
            except Exception:
                self._typed = None

    def __repr__(self):
        return "%s(conditions=%d)" % (self.__class__.__name__, self._condition_count)

    def and_(self):
        """Explicitly set the next connector to AND."""
        self._pending_connector = "AND"
        return self

    def or_(self):
        """Explicitly set the next connector to OR."""
        self._pending_connector = "OR"
        return self

    def filter(self, column, value, op="="):
        """Generic filter using a column name.

        :param column: Column name (string)
        :param value: Value to compare (int or str, or list for IN)
        :param op: Operator, e.g. '=', '!=', '<', '<=', '>', '>=', 'IN'
        """
        if not isinstance(column, str):
            column = str(column)
        is_string = isinstance(value, (str, bytes)) or (
            isinstance(value, (list, tuple)) and value and isinstance(value[0], (str, bytes))
        )
        return self._add_condition(column, op, value, is_string)

    @abstractmethod
    def to_dict(self, explode=False):
        """Return results as a dict of numpy arrays."""
        raise NotImplementedError

    @abstractmethod
    def to_df(self, explode=True):
        """Return results as a pandas.DataFrame."""
        raise NotImplementedError

    def _consume_connector(self):
        if self._condition_count == 0:
            return ""
        if self._pending_connector is None:
            return "AND"
        connector = self._pending_connector
        self._pending_connector = None
        return connector

    def _quote_string(self, value):
        if value is None:
            raise ValueError("String filter value cannot be None")
        s = str(value)
        has_single = "'" in s
        has_double = '"' in s
        if has_single and has_double:
            raise ValueError("String filter value contains both quote types: %r" % s)
        if has_double:
            return "'%s'" % s
        return '"%s"' % s

    def _format_value(self, value, is_string):
        if is_string:
            return self._quote_string(value)
        try:
            return str(int(value))
        except Exception as exc:
            raise ValueError("Expected integer value, got %r" % (value,)) from exc

    def _add_condition(self, column, op, value, is_string):
        connector = self._consume_connector()
        op_norm = op.upper()
        is_list = isinstance(value, (list, tuple))
        if op_norm != "IN" and is_list:
            op_norm = "IN"
        if op_norm == "IN":
            if not is_list:
                raise ValueError("IN operator expects a list/tuple of values")
            if not value:
                raise ValueError("IN operator expects a non-empty list of values")
            if self._typed is not None:
                if self._typed_has_condition:
                    if connector == "OR":
                        self._typed.orNext()
                    else:
                        self._typed.andNext()
                if is_string:
                    self._typed.in_(column, [str(v) for v in value])
                else:
                    self._typed.in_(column, [int(v) for v in value])
                self._typed_has_condition = True
        else:
            if self._typed is not None:
                if self._typed_has_condition:
                    if connector == "OR":
                        self._typed.orNext()
                    else:
                        self._typed.andNext()
                if is_string:
                    sval = str(value)
                    if op_norm == "=":
                        self._typed.eq(column, sval)
                    elif op_norm == "!=":
                        self._typed.ne(column, sval)
                    elif op_norm == "<":
                        self._typed.lt(column, sval)
                    elif op_norm == "<=":
                        self._typed.le(column, sval)
                    elif op_norm == ">":
                        self._typed.gt(column, sval)
                    elif op_norm == ">=":
                        self._typed.ge(column, sval)
                    else:
                        raise ValueError("Unsupported operator: %s" % op)
                else:
                    ival = int(value)
                    if op_norm == "=":
                        self._typed.eq(column, ival)
                    elif op_norm == "!=":
                        self._typed.ne(column, ival)
                    elif op_norm == "<":
                        self._typed.lt(column, ival)
                    elif op_norm == "<=":
                        self._typed.le(column, ival)
                    elif op_norm == ">":
                        self._typed.gt(column, ival)
                    elif op_norm == ">=":
                        self._typed.ge(column, ival)
                    else:
                        raise ValueError("Unsupported operator: %s" % op)
                self._typed_has_condition = True
        if self._typed is None:
            raise RuntimeError("ParquetFilter is not available in this environment")
        self._condition_count += 1
        return self


class ChromatogramQuery(ParquetQuery):
    """Chainable query builder for XICParquetFile chromatograms.

    Inherited methods:
      - filter(column, value, op="=")
      - filter_precursor_id(value, op="=")
      - filter_transition_id(value, op="=")
      - filter_ms_level(value, op="=")
      - filter_run_id(value, op="=")
      - filter_precursor_charge(value, op="=")
      - filter_product_charge(value, op="=")
      - filter_annotation(value, op="=")
      - filter_modified_sequence(value, op="=")
      - filter_transition_type(value, op="=")
      - to_dict(explode=False)
      - to_df(explode=True)

    Example::

        q = xic.query_chromatograms()
        df = q.filter_precursor_id(1381).filter_annotation("y3^1").to_df()
        df = q.filter_precursor_id([1381,2025]).filter_annotation(["y3^1", "y4^1", "b9^1"]).to_df()
    """

    def __init__(self, xic):
        self._xic = xic
        super().__init__()

    def filter_precursor_id(self, value, op="="):
        return self._add_condition("precursor_id", op, value, False)

    def filter_transition_id(self, value, op="="):
        return self._add_condition("transition_id", op, value, False)

    def filter_ms_level(self, value, op="="):
        return self._add_condition("ms_level", op, value, False)

    def filter_run_id(self, value, op="="):
        return self._add_condition("run_id", op, value, False)

    def filter_precursor_charge(self, value, op="="):
        return self._add_condition("precursor_charge", op, value, False)

    def filter_product_charge(self, value, op="="):
        return self._add_condition("product_charge", op, value, False)

    def filter_annotation(self, value, op="="):
        return self._add_condition("annotation", op, value, True)

    def filter_modified_sequence(self, value, op="="):
        return self._add_condition("modified_sequence", op, value, True)

    def filter_transition_type(self, value, op="="):
        return self._add_condition("transition_type", op, value, True)

    def to_dict(self, explode=False):
        """Return chromatogram data as dict of numpy arrays."""
        if self._typed is None:
            raise RuntimeError("ParquetFilter is not available in this environment")
        if not self._typed_has_condition:
            return self._xic.get_data_dict(explode=explode)
        return self._xic._get_data_dict_from_parquet_filter(self._typed, explode=explode)

    def to_df(self, explode=True):
        """Return chromatogram data as pandas.DataFrame."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError("pandas is required for to_df(). Install with `pip install pandas`.") from exc
        data = self.to_dict(explode=explode)
        return pd.DataFrame({k: list(v) for k, v in data.items()})
