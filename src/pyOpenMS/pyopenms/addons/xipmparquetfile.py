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

_ANALYTE_REVERSE_COLUMN_MAP = {v: k for k, v in _ANALYTE_COLUMN_MAP.items()}
_PRECURSOR_ANALYTE_FIELDS = (
    "precursor_id",
    "modified_sequence",
    "precursor_charge",
    "precursor_decoy",
)
_PRECURSOR_GROUPING_FIELDS = (
    "precursor_id",
    "modified_sequence",
    "precursor_charge",
)
_TRANSITION_ANALYTE_FIELDS = (
    "transition_id",
    "product_charge",
    "detecting_transition",
    "product_decoy",
    "transition_ordinal",
    "transition_type",
    "annotation",
)
_NUMERIC_PUSHDOWN_COLUMNS = {
    "precursor_id",
    "transition_id",
    "precursor_charge",
    "product_charge",
    "ms_level",
    "run_id",
}
_STRING_PUSHDOWN_COLUMNS = {
    "modified_sequence",
    "peakmap_type",
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


def _read_metadata_columns(self, columns):
    try:
        import pyarrow.dataset as ds
    except ImportError:
        return None

    filenames = [str(filename) for filename in self.getFilenames()]
    if not filenames:
        return {}

    table = ds.dataset(filenames, format="parquet").to_table(columns=list(columns))
    return {
        _ANALYTE_COLUMN_MAP.get(column, column.lower()): table[column].to_pylist()
        for column in columns
    }


def _deduplicate_transition_lists(group, requested_transition_fields):
    if not requested_transition_fields:
        return

    seen = set()
    keep = []
    length = len(group[requested_transition_fields[0]])
    for i in range(length):
        key = tuple(group[field][i] for field in requested_transition_fields)
        if key in seen:
            continue
        seen.add(key)
        keep.append(i)

    for field in requested_transition_fields:
        group[field] = [group[field][i] for i in keep]


def _build_analyte_dict(data, requested, nest_transitions):
    if not data:
        return {column: [] for column in requested}

    n_rows = len(next(iter(data.values()), []))
    if n_rows == 0:
        return {column: [] for column in requested}

    requested_precursor_fields = [field for field in _PRECURSOR_ANALYTE_FIELDS if field in requested]
    grouping_precursor_fields = [field for field in _PRECURSOR_GROUPING_FIELDS if field in requested]
    requested_transition_fields = [field for field in _TRANSITION_ANALYTE_FIELDS if field in requested]

    if not nest_transitions:
        seen = set()
        rows = []
        for i in range(n_rows):
            record = {field: data[field][i] for field in requested}
            key = tuple(record[field] for field in requested)
            if key in seen:
                continue
            seen.add(key)
            rows.append(record)
        return {column: [row[column] for row in rows] for column in requested}

    if not grouping_precursor_fields:
        raise RuntimeError(
            "nest_transitions=True requires at least one precursor discriminator column: "
            "PRECURSOR_ID, MODIFIED_SEQUENCE, or PRECURSOR_CHARGE"
        )

    grouped = {}
    order = []
    transition_ids = data.get("transition_id")
    for i in range(n_rows):
        key = tuple(data[field][i] for field in grouping_precursor_fields)
        if key not in grouped:
            grouped[key] = {
                field: data[field][i] for field in requested_precursor_fields
            }
            for field in requested_transition_fields:
                grouped[key][field] = []
            order.append(key)
        elif "precursor_decoy" in requested_precursor_fields:
            current_decoy = data["precursor_decoy"][i]
            if current_decoy is not None:
                grouped_decoy = grouped[key]["precursor_decoy"]
                if grouped_decoy is None or current_decoy == 0:
                    grouped[key]["precursor_decoy"] = current_decoy

        if requested_transition_fields and transition_ids is not None and transition_ids[i] is not None:
            for field in requested_transition_fields:
                grouped[key][field].append(data[field][i])

    for key in order:
        _deduplicate_transition_lists(grouped[key], requested_transition_fields)

    result = {}
    for column in requested:
        if column in requested_precursor_fields:
            result[column] = [grouped[key][column] for key in order]
        else:
            result[column] = [grouped[key].get(column, []) for key in order]
    return result


def _normalize_pushdown_value(column, value):
    if column in _NUMERIC_PUSHDOWN_COLUMNS:
        return int(value)
    if column in _STRING_PUSHDOWN_COLUMNS:
        return str(value)
    raise KeyError(f"Unsupported XIPM pushdown column: {column}")


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

    def _split_filters(self):
        pushdown = {
            "precursor_id": -1,
            "transition_id": -1,
            "modified_sequence": "",
            "precursor_charge": -1,
            "product_charge": -1,
            "ms_level": -1,
            "run_id": -1,
            "peakmap_type": "",
        }
        residual = []

        for column, value, op in self._filters:
            if op not in ("=", "=="):
                residual.append((column, value, op))
                continue
            if isinstance(value, (list, tuple)):
                residual.append((column, value, op))
                continue
            if column not in pushdown:
                residual.append((column, value, op))
                continue

            normalized = _normalize_pushdown_value(column, value)
            default_value = pushdown[column]
            if default_value in (-1, "") or default_value == normalized:
                pushdown[column] = normalized
            else:
                residual.append((column, value, op))

        return pushdown, residual

    def to_dict(self, explode=False):
        """Return peak-map data as dict."""
        pushdown, residual = self._split_filters()
        data = self._xipm.get_data_dict(explode=False, **pushdown)
        filtered = _filter_peak_map_dict(data, residual) if residual else data
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
    metadata = _read_metadata_columns(self, ["RUN_ID", "SOURCE_FILE"])
    if metadata is None:
        return self.getRuns()

    seen = set()
    run_ids = []
    source_files = []
    for run_id, source_file in zip(metadata["run_id"], metadata["source_file"]):
        key = (run_id, source_file)
        if key in seen:
            continue
        seen.add(key)
        run_ids.append(run_id)
        source_files.append(source_file)
    return {"run_id": run_ids, "source_file": source_files}


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
    parquet_columns = [_ANALYTE_REVERSE_COLUMN_MAP[column] for column in requested]
    if nest_transitions and requested and "TRANSITION_ID" not in parquet_columns:
        if any(field in requested for field in _TRANSITION_ANALYTE_FIELDS):
            parquet_columns.append("TRANSITION_ID")
    data = _read_metadata_columns(self, parquet_columns)
    if data is None:
        data = self.get_data_dict(explode=False)
    return _build_analyte_dict(data, requested, nest_transitions)


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
