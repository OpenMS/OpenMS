"""Pure Python addon for FeatureFinderAlgorithmMetaboIdent.

Adds a Pythonic, pandas-friendly interface for building the target compound
library, so users no longer have to construct ``FeatureFinderMetaboIdentCompound``
objects by hand.

The accepted DataFrame columns mirror the input TSV of the ``FeatureFinderMetaboIdent``
TOPP tool (see its documentation), including the optional ``IonMobility`` and
``Adduct`` columns:

    CompoundName  SumFormula  Mass  Charge  RetentionTime  RetentionTimeRange  IsoDistribution  [IonMobility]  [Adduct]
"""
from __future__ import annotations

from . import addon


# Canonical column name -> set of accepted aliases (all matched case-insensitively,
# ignoring surrounding whitespace). The canonical name itself is always accepted.
_COLUMN_ALIASES = {
    "CompoundName": ("compoundname", "compound_name", "name"),
    "SumFormula": ("sumformula", "sum_formula", "formula"),
    "Mass": ("mass",),
    "Charge": ("charge", "charges"),
    "RetentionTime": ("retentiontime", "retention_time", "rt"),
    "RetentionTimeRange": ("retentiontimerange", "retention_time_range", "rt_range", "rtrange"),
    "IsoDistribution": ("isodistribution", "iso_distribution", "isotope_distribution"),
    "IonMobility": ("ionmobility", "ion_mobility", "im"),
    "Adduct": ("adduct", "adducts"),
}

# Columns that must be present for a valid compound library.
_REQUIRED_COLUMNS = (
    "CompoundName", "SumFormula", "Mass", "Charge",
    "RetentionTime", "RetentionTimeRange", "IsoDistribution",
)

# Build a flat alias -> canonical lookup once at import time.
_ALIAS_TO_CANONICAL = {}
for _canonical, _aliases in _COLUMN_ALIASES.items():
    _ALIAS_TO_CANONICAL[_canonical.lower()] = _canonical
    for _alias in _aliases:
        _ALIAS_TO_CANONICAL[_alias] = _canonical


def _normalize_columns(df):
    """Return ``df`` with its columns renamed to canonical names.

    Raises ``ValueError`` if two source columns map to the same canonical name
    (e.g. both ``rt`` and ``RetentionTime``), which would otherwise silently
    create ambiguous duplicate columns.
    """
    rename = {}
    canonical_to_source = {}
    for col in df.columns:
        canonical = _ALIAS_TO_CANONICAL.get(str(col).lower().strip())
        if canonical is None:
            continue  # leave unrecognized columns untouched
        if canonical in canonical_to_source:
            raise ValueError(
                f"Columns '{canonical_to_source[canonical]}' and '{col}' both map "
                f"to '{canonical}'. Please provide only one."
            )
        canonical_to_source[canonical] = col
        rename[col] = canonical
    return df.rename(columns=rename)


def _to_list(value, dtype, np):
    """Convert a scalar / list / numpy array / comma-separated string to ``list[dtype]``.

    Empty strings, ``None`` and NaN scalars yield an empty list. Integer parsing
    accepts float-like tokens (``"1.0"`` -> ``1``) so spreadsheet exports work.
    """
    if value is None:
        return []
    # Scalar NaN (float or numpy float) -> treat as "not provided".
    if isinstance(value, float) and np.isnan(value):
        return []
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return []
        return [_coerce(tok.strip(), dtype) for tok in stripped.split(",") if tok.strip()]
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_coerce(x, dtype) for x in value]
    return [_coerce(value, dtype)]


def _coerce(x, dtype):
    """Coerce a single token to ``dtype`` (int parsing tolerates float-like tokens)."""
    if dtype is int:
        return int(float(x))
    return dtype(x)


@addon("FeatureFinderAlgorithmMetaboIdent", "compounds_from_df")
@staticmethod
def compounds_from_df(df):
    """
    Convert a pandas DataFrame to a list of ``FeatureFinderMetaboIdentCompound`` objects.

    Column names are case-insensitive and accept ``snake_case`` and common
    aliases (e.g. ``rt`` for ``RetentionTime``, ``im`` for ``IonMobility``).
    The columns mirror the input TSV of the ``FeatureFinderMetaboIdent`` TOPP tool.

    Parameters
    ----------
    df : pandas.DataFrame
        Table with the following columns (required unless noted):

        - ``CompoundName`` (``compound_name``, ``name``): str, unique identifier.
        - ``SumFormula`` (``sum_formula``, ``formula``): str, chemical formula
          (may be empty if ``Mass`` is given).
        - ``Mass`` (``mass``): float, neutral mass; ``0`` = derive from formula.
        - ``Charge`` (``charges``): int, list[int], or comma-separated str.
        - ``RetentionTime`` (``retention_time``, ``rt``): float, list[float], or
          comma-separated str.
        - ``RetentionTimeRange`` (``rt_range``): float/list/str; ``0`` = use the
          ``extract:rt_window`` parameter.
        - ``IsoDistribution`` (``iso_distribution``): float/list/str; ``0`` =
          calculate from formula.
        - ``IonMobility`` (``ion_mobility``, ``im``) *(optional)*: float/list/str;
          one value, or one per RT entry; absent or ``0`` disables IM filtering.
        - ``Adduct`` (``adducts``) *(optional)*: str, e.g. ``"[M+H]+"``,
          ``"[M+Na]+"``, ``"[M-H]-"``; if omitted, ``[M+H]+`` / ``[M-H]-`` is
          assumed from charge polarity.

    Returns
    -------
    list of FeatureFinderMetaboIdentCompound

    Raises
    ------
    ValueError
        If required columns are missing, two columns map to the same field,
        compound names are empty or duplicated, or ``Charge``/``RetentionTime``
        are empty for any compound.
    ImportError
        If pandas or numpy are not installed.

    Examples
    --------
    >>> import pandas as pd
    >>> from pyopenms import FeatureFinderAlgorithmMetaboIdent
    >>> df = pd.DataFrame({
    ...     'CompoundName': ['glucose', 'fructose'],
    ...     'SumFormula': ['C6H12O6', 'C6H12O6'],
    ...     'Mass': [0.0, 0.0],
    ...     'Charge': [-1, '-1,1'],
    ...     'RetentionTime': [123.4, '200.0,250.0'],
    ...     'RetentionTimeRange': [0.0, 0.0],
    ...     'IsoDistribution': [0.0, 0.0],
    ...     'IonMobility': [0.0, 0.95],   # optional
    ...     'Adduct': ['[M-H]-', ''],     # optional
    ... })
    >>> compounds = FeatureFinderAlgorithmMetaboIdent.compounds_from_df(df)
    """
    try:
        import numpy as np
        import pandas as pd
    except ImportError as e:
        raise ImportError(
            "pandas and numpy are required for compounds_from_df. "
            "Install with: pip install pandas numpy"
        ) from e

    from pyopenms._pyopenms_featurefinder import FeatureFinderMetaboIdentCompound

    df = _normalize_columns(df)

    missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns: {missing}. "
            f"Required: {list(_REQUIRED_COLUMNS)} (optional: IonMobility, Adduct)."
        )

    has_im = "IonMobility" in df.columns
    has_adduct = "Adduct" in df.columns

    compounds = []
    seen_names = set()

    for idx, row in df.iterrows():
        name = row["CompoundName"]
        if pd.isna(name) or (isinstance(name, str) and not name.strip()):
            raise ValueError(f"Row {idx}: CompoundName cannot be empty")
        name = str(name).strip()
        if name in seen_names:
            raise ValueError(f"Row {idx}: Duplicate CompoundName '{name}'")
        seen_names.add(name)

        formula = row["SumFormula"]
        formula = "" if pd.isna(formula) else str(formula).strip()

        mass = row["Mass"]
        mass = 0.0 if pd.isna(mass) else float(mass)

        try:
            charges = _to_list(row["Charge"], int, np)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Row {idx} ('{name}'): could not parse Charge "
                             f"value {row['Charge']!r}: {e}") from e
        if not charges:
            raise ValueError(f"Row {idx} ('{name}'): Charge cannot be empty")

        try:
            rts = _to_list(row["RetentionTime"], float, np)
        except (ValueError, TypeError) as e:
            raise ValueError(f"Row {idx} ('{name}'): could not parse RetentionTime "
                             f"value {row['RetentionTime']!r}: {e}") from e
        if not rts:
            raise ValueError(f"Row {idx} ('{name}'): RetentionTime cannot be empty")

        rt_ranges = _to_list(row["RetentionTimeRange"], float, np) or [0.0]
        iso_dist = _to_list(row["IsoDistribution"], float, np) or [0.0]

        ion_mobilities = _to_list(row["IonMobility"], float, np) if has_im else []

        adduct = ""
        if has_adduct:
            a = row["Adduct"]
            adduct = "" if pd.isna(a) else str(a).strip()

        compounds.append(
            FeatureFinderMetaboIdentCompound(
                name, formula, mass, charges, rts, rt_ranges, iso_dist,
                ion_mobilities, adduct,
            )
        )

    return compounds


@addon("FeatureFinderAlgorithmMetaboIdent")
def run_from_df(self, df, features, spectra_path=""):
    """
    Run targeted feature extraction using compounds from a pandas DataFrame.

    Convenience wrapper that combines :meth:`compounds_from_df` and :meth:`run`.
    The spectra must be set beforehand via :meth:`setMSData`.

    Parameters
    ----------
    df : pandas.DataFrame
        Compound table (see :meth:`compounds_from_df` for the accepted columns,
        including the optional ``IonMobility`` and ``Adduct`` columns).
    features : FeatureMap
        Output map; detected features are stored here (modified in-place).
    spectra_path : str, optional
        Path to the spectra file, used as a fallback for the primary MS run path
        annotated in the feature map.

    Examples
    --------
    >>> exp = MSExperiment()
    >>> MzMLFile().load("spectra.mzML", exp)
    >>> ff = FeatureFinderAlgorithmMetaboIdent()
    >>> ff.setMSData(exp)
    >>> features = FeatureMap()
    >>> ff.run_from_df(df, features, "spectra.mzML")
    """
    compounds = type(self).compounds_from_df(df)
    self.run(compounds, features, spectra_path)
