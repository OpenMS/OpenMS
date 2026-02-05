"""XICParquetFile addon methods for DataFrame support."""
from . import addon


@addon("XICParquetFile")
def __len__(self):
    """Return the number of files associated with this instance."""
    try:
        return len(self.getFilenames())
    except Exception:
        return 0


@addon("XICParquetFile")
def __repr__(self):
    """Return a detailed string representation."""
    try:
        files = self.getFilenames()
        files = [f if isinstance(f, str) else str(f) for f in files]
        return f"XICParquetFile(n_files={len(files)}, files={files!r})"
    except Exception:
        return "XICParquetFile(<unavailable>)"


@addon("XICParquetFile")
def __str__(self):
    """Return a concise string representation."""
    try:
        return f"XICParquetFile(n_files={len(self.getFilenames())})"
    except Exception:
        return "XICParquetFile(n_files=?)"


@addon("XICParquetFile")
def df_columns(self):
    """
    Return the parquet schema column names.

    Returns
    -------
    list
        List of column names.
    """
    return self.getColumns()


@addon("XICParquetFile")
def get_data_dict(self, explode=False, precursor_id=-1, transition_id=-1,
                  modified_sequence="", precursor_charge=-1, product_charge=-1,
                  ms_level=-1, run_id=-1, filter=""):
    """
    Return chromatogram data as a dict.

    If explode=True, returns long format with rt/intensity rows.
    Otherwise, rt and intensity are stored as lists.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per RT/intensity.
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
    filter : str
        Optional filter expression string.

    Returns
    -------
    dict
        Dict of lists keyed by column name.
    """
    return self.getChromatograms(
        precursor_id=precursor_id,
        transition_id=transition_id,
        modified_sequence=modified_sequence,
        precursor_charge=precursor_charge,
        product_charge=product_charge,
        ms_level=ms_level,
        run_id=run_id,
        filter=filter,
        explode=explode
    )


@addon("XICParquetFile")
def get_run_dict(self):
    """
    Return unique run metadata as a dict.

    Returns
    -------
    dict
        Dict with run_id and source_file lists.
    """
    return self.getRuns()


@addon("XICParquetFile")
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


@addon("XICParquetFile")
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


@addon("XICParquetFile")
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


@addon("XICParquetFile")
def to_df(self, explode=False):
    """
    Return chromatogram data as a pandas DataFrame.

    If explode=True, returns long format with rt/intensity rows.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per RT/intensity.

    Returns
    -------
    pandas.DataFrame
        DataFrame with chromatogram data.

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


@addon("XICParquetFile")
def to_arrow(self, explode=False):
    """
    Returns an Apache Arrow Table representation of the chromatograms.

    If explode=True, returns long format with rt/intensity rows.

    Parameters
    ----------
    explode : bool
        If True, return long format with one row per RT/intensity.

    Returns
    -------
    pyarrow.Table
        Arrow Table with chromatogram data.

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
