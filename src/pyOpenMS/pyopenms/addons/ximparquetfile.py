"""XIMParquetFile addon methods for DataFrame support."""

from . import addon


class _MobilogramQuery:
    """Chainable query builder for XIMParquetFile mobilograms.

    Example::
        q = xim.query_mobilograms()
        df = q.filter_precursor_id(123).to_df()
    """

    def __init__(self, xim):
        self._xim = xim
        self._filters = []

    def __repr__(self):
        return f"MobilogramQuery(conditions={len(self._filters)})"

    @staticmethod
    def _format_filter_value(value):
        if isinstance(value, str):
            escaped = value.replace("\\", "\\\\").replace('"', '\\"')
            return f'"{escaped}"'
        return str(value)

    def _add_filter(self, column, value, op="="):
        if isinstance(value, (list, tuple)):
            values_str = ",".join(self._format_filter_value(v) for v in value)
            self._filters.append(f"{column} IN ({values_str})")
        else:
            self._filters.append(f"{column}{op}{self._format_filter_value(value)}")
        return self

    def filter_precursor_id(self, value, op="="):
        return self._add_filter("precursor_id", value, op)

    def filter_transition_id(self, value, op="="):
        return self._add_filter("transition_id", value, op)

    def filter_ms_level(self, value, op="="):
        return self._add_filter("ms_level", value, op)

    def filter_run_id(self, value, op="="):
        return self._add_filter("run_id", value, op)

    def filter_precursor_charge(self, value, op="="):
        return self._add_filter("precursor_charge", value, op)

    def filter_product_charge(self, value, op="="):
        return self._add_filter("product_charge", value, op)

    def filter_annotation(self, value, op="="):
        return self._add_filter("annotation", value, op)

    def filter_modified_sequence(self, value, op="="):
        return self._add_filter("modified_sequence", value, op)

    def filter_transition_type(self, value, op="="):
        return self._add_filter("transition_type", value, op)

    def to_dict(self, explode=False):
        """Return mobilogram data as dict."""
        filter_str = " AND ".join(self._filters) if self._filters else ""
        return self._xim.get_data_dict(explode=explode, filter_expr=filter_str)

    def to_df(self, explode=True):
        """Return mobilogram data as a pandas DataFrame."""
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError(
                "pandas is required for to_df(). Install with `pip install pandas`."
            ) from e
        data = self.to_dict(explode=explode)
        return pd.DataFrame({k: list(v) for k, v in data.items()})


@addon("XIMParquetFile")
def __len__(self):
    """Return the number of files associated with this instance."""
    try:
        filenames = self.getFilenames()
        return len(filenames) if filenames is not None else 0
    except AttributeError:
        return 0


@addon("XIMParquetFile")
def df_columns(self):
    """
    Return the parquet schema column names.

    Returns
    -------
    list
        List of column names.
    """
    return self.getColumns()


@addon("XIMParquetFile")
def get_data_dict(
    self,
    explode=False,
    precursor_id=-1,
    transition_id=-1,
    modified_sequence="",
    precursor_charge=-1,
    product_charge=-1,
    ms_level=-1,
    run_id=-1,
    mobilogram_type="",
    feature_id=-1,
    feature_rt=-1.0,
    filter_expr="",
    **kwargs,
):
    """
    Return mobilogram data as a dict.

    If explode=True, returns long format with mobility/intensity rows.
    Otherwise, mobility and intensity are stored as lists.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per mobility/intensity.
    precursor_id : int
        Optional precursor id (-1 to ignore).
    transition_id : int
        Optional transition id (-1 to ignore).
    modified_sequence : str
        Optional modified sequence filter (empty to ignore).
    precursor_charge : int
        Optional precursor charge filter (-1 to ignore).
    product_charge : int
        Optional product charge filter (-1 to ignore).
    ms_level : int
        Optional MS level filter (-1 to ignore).
    run_id : int
        Optional run id filter (-1 to ignore).
    mobilogram_type : str
        Optional mobilogram type filter (empty to ignore).
    feature_id : int
        Optional feature id filter (-1 to ignore).
    feature_rt : float
        Optional feature RT filter (<0 to ignore).
    filter_expr : str
        Optional filter expression string.

    Returns
    -------
    dict
        Dict of lists keyed by column name.
    """
    return self.getMobilograms(
        precursor_id=precursor_id,
        transition_id=transition_id,
        modified_sequence=modified_sequence,
        precursor_charge=precursor_charge,
        product_charge=product_charge,
        ms_level=ms_level,
        run_id=run_id,
        mobilogram_type=mobilogram_type,
        feature_id=feature_id,
        feature_rt=feature_rt,
        filter=filter_expr,
        explode=explode,
    )


@addon("XIMParquetFile")
def get_run_dict(self):
    """
    Return unique run metadata as a dict.

    Returns
    -------
    dict
        Dict with run_id and source_file lists.
    """
    return self.getRuns()


@addon("XIMParquetFile")
def get_run_df(self):
    """
    Return unique run metadata as a pandas DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame with run_id and source_file.

    Raises
    ------
    ImportError
        If pandas is not installed.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for this method. Install with `pip install pandas`."
        ) from e
    data = self.get_run_dict()
    return pd.DataFrame(data)


@addon("XIMParquetFile")
def get_analyte_dict(self, nest_transitions=True, columns=None):
    """
    Return unique analyte metadata as a dict.

    If nest_transitions=False, each row represents a unique precursor-transition
    pair with scalar transition fields. If nest_transitions=True, each row
    represents a unique precursor and transition fields are lists.

    Parameters
    ----------
    nest_transitions : bool
        Aggregate transition fields per precursor.
    columns : list, optional
        List of column names to return. If None, uses all columns.

    Returns
    -------
    dict
        Dict of lists keyed by column name.
    """
    return self.getAnalytes(nest_transitions=nest_transitions, columns=columns)


@addon("XIMParquetFile")
def get_analyte_df(self, nest_transitions=True, columns=None):
    """
    Return unique analyte metadata as a pandas DataFrame.

    Parameters
    ----------
    nest_transitions : bool
        Aggregate transition fields per precursor.
    columns : list, optional
        List of column names to return.

    Returns
    -------
    pandas.DataFrame
        DataFrame with analyte metadata.

    Raises
    ------
    ImportError
        If pandas is not installed.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for this method. Install with `pip install pandas`."
        ) from e
    data = self.get_analyte_dict(nest_transitions=nest_transitions, columns=columns)
    return pd.DataFrame({k: list(v) for k, v in data.items()})


@addon("XIMParquetFile")
def to_df(self, explode=False):
    """
    Return mobilogram data as a pandas DataFrame.

    If explode=True, returns long format with mobility/intensity rows.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per mobility/intensity.

    Returns
    -------
    pandas.DataFrame
        DataFrame with mobilogram data.

    Raises
    ------
    ImportError
        If pandas is not installed.
    """
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for to_df(). Install with `pip install pandas`."
        ) from e
    data = self.get_data_dict(explode=explode)
    if explode:
        return pd.DataFrame(data)
    # Convert to plain Python lists to avoid pandas inferring 2D arrays.
    clean = {}
    for k, v in data.items():
        if hasattr(v, "ndim") and v.ndim > 1:
            clean[k] = [row for row in v]
        else:
            clean[k] = list(v)
    return pd.DataFrame(clean)


@addon("XIMParquetFile")
def to_arrow(self, explode=False):
    """
    Returns an Apache Arrow Table representation of the mobilograms.

    If explode=True, returns long format with mobility/intensity rows.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per mobility/intensity.

    Returns
    -------
    pyarrow.Table
        Arrow Table with mobilogram data.

    Raises
    ------
    ImportError
        If pyarrow is not installed.
    """
    try:
        import pyarrow as pa
    except ImportError as e:
        raise ImportError(
            "pyarrow is required for to_arrow(). Install with `pip install pyarrow`."
        ) from e
    return pa.Table.from_pydict(self.get_data_dict(explode=explode))


@addon("XIMParquetFile")
def query_mobilograms(self):
    """
    Return a chainable query builder for mobilograms.

    Example::
        df = xim.query_mobilograms().filter_precursor_id(123).to_df()

    Returns
    -------
    _MobilogramQuery
        Query builder with filter methods and to_df()/to_dict().
    """
    return _MobilogramQuery(self)
