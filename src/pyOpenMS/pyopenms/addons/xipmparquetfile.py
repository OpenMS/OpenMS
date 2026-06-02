"""XIPMParquetFile addon methods for DataFrame support."""

from . import addon


_ANALYTE_COLUMN_MAP = {
    "PRECURSOR_ID": "precursor_id",
    "MODIFIED_SEQUENCE": "modified_sequence",
    "PRECURSOR_CHARGE": "precursor_charge",
    "PRECURSOR_DECOY": "precursor_decoy",
    "TRANSITION_ID": "transition_id",
    "PRODUCT_CHARGE": "product_charge",
    "DETECTING_TRANSITION": "detecting_transition",
    "PRODUCT_DECOY": "product_decoy",
    "TRANSITION_ORDINAL": "transition_ordinal",
    "TRANSITION_TYPE": "transition_type",
    "ANNOTATION": "annotation",
}


def _compare(lhs, rhs, op):
    if op == "=":
        return lhs == rhs
    if op == "!=":
        return lhs != rhs
    if lhs is None:
        return False
    if op == "<":
        return lhs < rhs
    if op == "<=":
        return lhs <= rhs
    if op == ">":
        return lhs > rhs
    if op == ">=":
        return lhs >= rhs
    raise ValueError(f"Unsupported operator: {op}")


def _value_matches(lhs, rhs, op):
    if isinstance(rhs, (list, tuple)):
        if op not in ("=", "IN"):
            raise ValueError("List filters only support '=' or 'IN'")
        return lhs in rhs
    return _compare(lhs, rhs, op)


def _filter_peak_map_dict(data, filters):
    if not data:
        return {}
    n_rows = len(next(iter(data.values())))
    keep = []
    for i in range(n_rows):
        matched = True
        for column, value, op in filters:
            if column not in data:
                raise KeyError(f"Unknown XIPMParquetFile column: {column}")
            if not _value_matches(data[column][i], value, op):
                matched = False
                break
        if matched:
            keep.append(i)
    return {k: [v[i] for i in keep] for k, v in data.items()}


def _explode_peak_map_dict(data):
    if not data:
        return {}
    exploded = {k: [] for k in data}
    n_rows = len(data["run_id"])
    for i in range(n_rows):
        mzs = list(data["mz"][i])
        rts = list(data["rt"][i])
        ims = list(data["ion_mobility"][i])
        intensities = list(data["intensity"][i])
        if not (len(mzs) == len(rts) == len(ims) == len(intensities)):
            raise RuntimeError("XIPMParquetFile: mz/rt/ion_mobility/intensity length mismatch")
        if not mzs:
            continue
        for j in range(len(mzs)):
            for key, values in data.items():
                if key == "mz":
                    exploded[key].append(mzs[j])
                elif key == "rt":
                    exploded[key].append(rts[j])
                elif key == "ion_mobility":
                    exploded[key].append(ims[j])
                elif key == "intensity":
                    exploded[key].append(intensities[j])
                else:
                    exploded[key].append(values[i])
    return exploded


def _normalize_analyte_columns(columns):
    if columns is None:
        return list(_ANALYTE_COLUMN_MAP.values())

    normalized = []
    for column in columns:
        upper = str(column).upper()
        if upper not in _ANALYTE_COLUMN_MAP:
            raise RuntimeError(f"Unsupported XIPM analyte column: {column}")
        normalized.append(_ANALYTE_COLUMN_MAP[upper])
    return normalized


class _PeakMapQuery:
    """Chainable query builder for XIPMParquetFile peak maps."""

    def __init__(self, xipm):
        self._xipm = xipm
        self._filters = []

    def __repr__(self):
        return f"PeakMapQuery(conditions={len(self._filters)})"

    def _add_filter(self, column, value, op="="):
        self._filters.append((column, value, op))
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

    def filter_peakmap_type(self, value, op="="):
        return self._add_filter("peakmap_type", value, op)

    def to_dict(self, explode=False):
        """Return peak-map data as dict."""
        data = self._xipm.get_data_dict(explode=False)
        filtered = _filter_peak_map_dict(data, self._filters)
        if explode:
            return _explode_peak_map_dict(filtered)
        return filtered

    def to_df(self, explode=True):
        """Return peak-map data as a pandas DataFrame."""
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError(
                "pandas is required for to_df(). Install with `pip install pandas`."
            ) from e
        data = self.to_dict(explode=explode)
        return pd.DataFrame({k: list(v) for k, v in data.items()})


@addon("XIPMParquetFile")
def __len__(self):
    """Return the number of files associated with this instance."""
    try:
        return len(self.getFilenames())
    except Exception:
        return 0


@addon("XIPMParquetFile")
def df_columns(self):
    """Return the parquet schema column names."""
    return self.getColumns()


@addon("XIPMParquetFile")
def get_data_dict(self, explode=False, precursor_id=-1, transition_id=-1,
                  modified_sequence="", precursor_charge=-1, product_charge=-1,
                  ms_level=-1, run_id=-1, peakmap_type=""):
    """
    Return peak-map data as a dict.

    If explode=True, returns long format with one row per mz/rt/ion_mobility/intensity point.
    Otherwise, mz/rt/ion_mobility/intensity are stored as lists per extracted peak map.
    """
    if explode:
        return self.getPeakMaps(
            precursor_id=precursor_id,
            transition_id=transition_id,
            modified_sequence=modified_sequence,
            precursor_charge=precursor_charge,
            product_charge=product_charge,
            ms_level=ms_level,
            run_id=run_id,
            peakmap_type=peakmap_type,
            explode=True,
        )
    return self.getPeakMaps(
        precursor_id=precursor_id,
        transition_id=transition_id,
        modified_sequence=modified_sequence,
        precursor_charge=precursor_charge,
        product_charge=product_charge,
        ms_level=ms_level,
        run_id=run_id,
        peakmap_type=peakmap_type,
        explode=False,
    )


@addon("XIPMParquetFile")
def get_run_dict(self):
    """Return unique run metadata as a dict."""
    return self.getRuns()


@addon("XIPMParquetFile")
def get_run_df(self):
    """Return unique run metadata as a pandas DataFrame."""
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for this method. Install with `pip install pandas`."
        ) from e
    return pd.DataFrame(self.get_run_dict())


@addon("XIPMParquetFile")
def get_analyte_dict(self, nest_transitions=True, columns=None):
    """
    Return unique analyte metadata as a dict.

    If nest_transitions=False, each row represents a unique precursor-transition
    pair with scalar transition fields. If nest_transitions=True, each row
    represents a unique precursor and transition fields are lists.
    """
    requested = _normalize_analyte_columns(columns)
    data = self.get_data_dict(explode=False)

    if not data or not data["run_id"]:
        return {column: [] for column in requested}

    precursor_fields = [
        "precursor_id",
        "modified_sequence",
        "precursor_charge",
        "precursor_decoy",
    ]
    transition_fields = [
        "transition_id",
        "product_charge",
        "detecting_transition",
        "product_decoy",
        "transition_ordinal",
        "transition_type",
        "annotation",
    ]

    if not nest_transitions:
        seen = set()
        rows = []
        for i in range(len(data["run_id"])):
            record = {field: data[field][i] for field in precursor_fields + transition_fields}
            key = tuple(record[field] for field in precursor_fields + transition_fields)
            if key in seen:
                continue
            seen.add(key)
            rows.append(record)
        return {column: [row[column] for row in rows] for column in requested}

    grouped = {}
    order = []
    for i in range(len(data["run_id"])):
        key = tuple(data[field][i] for field in precursor_fields)
        if key not in grouped:
            grouped[key] = {
                "precursor_id": data["precursor_id"][i],
                "modified_sequence": data["modified_sequence"][i],
                "precursor_charge": data["precursor_charge"][i],
                "precursor_decoy": data["precursor_decoy"][i],
                "transition_id": [],
                "product_charge": [],
                "detecting_transition": [],
                "product_decoy": [],
                "transition_ordinal": [],
                "transition_type": [],
                "annotation": [],
            }
            order.append(key)

        if data["transition_id"][i] is not None:
            for field in transition_fields:
                grouped[key][field].append(data[field][i])

    return {column: [grouped[key][column] for key in order] for column in requested}


@addon("XIPMParquetFile")
def get_analyte_df(self, nest_transitions=True, columns=None):
    """Return unique analyte metadata as a pandas DataFrame."""
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for this method. Install with `pip install pandas`."
        ) from e
    data = self.get_analyte_dict(nest_transitions=nest_transitions, columns=columns)
    return pd.DataFrame({k: list(v) for k, v in data.items()})


@addon("XIPMParquetFile")
def to_df(self, explode=False):
    """Return peak-map data as a pandas DataFrame."""
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas is required for to_df(). Install with `pip install pandas`."
        ) from e
    data = self.get_data_dict(explode=explode)
    if explode:
        return pd.DataFrame(data)
    clean = {}
    for k, v in data.items():
        if hasattr(v, "ndim") and v.ndim > 1:
            clean[k] = [row for row in v]
        else:
            clean[k] = list(v)
    return pd.DataFrame(clean)


@addon("XIPMParquetFile")
def to_arrow(self, explode=False):
    """Return peak-map data as an Apache Arrow Table."""
    try:
        import pyarrow as pa
    except ImportError as e:
        raise ImportError(
            "pyarrow is required for to_arrow(). Install with `pip install pyarrow`."
        ) from e
    return pa.Table.from_pydict(self.get_data_dict(explode=explode))


@addon("XIPMParquetFile")
def query_peak_maps(self):
    """Return a chainable query builder for peak maps."""
    return _PeakMapQuery(self)
