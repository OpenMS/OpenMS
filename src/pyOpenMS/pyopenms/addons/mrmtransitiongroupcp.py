"""MRMTransitionGroupCP addon methods for DataFrame support."""
import numpy as np
from . import addon, pin_string_dtype, register_element_views, string_dtype


@addon("MRMTransitionGroupCP")
def chromatogram_df_columns(self, columns='default', export_meta_values=True):
    """Returns a list of column names that to_chromatogram_df() would produce."""
    if self._chromatogram_count():
        return self.chromatogram_view(0).df_columns(columns=columns, export_meta_values=export_meta_values)
    return ['rt', 'intensity', 'precursor_mz', 'precursor_charge', 'product_mz', 'native_id']


@addon("MRMTransitionGroupCP")
def feature_df_columns(self, columns='default'):
    """Returns a list of column names that to_feature_df() would produce."""
    cols = ['feature_id', 'rt', 'intensity', 'quality']

    if columns == 'all':
        meta_values = set()
        for f in self.iter_feature_views():
            mvs = []
            f.getKeys(mvs)
            for m in mvs:
                meta_values.add(m.decode() if isinstance(m, bytes) else m)
        cols.extend(sorted(meta_values))

    return cols


@addon("MRMTransitionGroupCP")
def to_chromatogram_df(self, columns=None, export_meta_values=True):
    """Returns a DataFrame representation of the Chromatograms stored in MRMTransitionGroupCP."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for to_chromatogram_df(). Install with: pip install pandas")
    chroms = self.chromatogram_views()  # zero-copy read; to_df() only reads
    out = [c.to_df(columns=columns, export_meta_values=export_meta_values) for c in chroms]
    if out:
        return pd.concat(out, ignore_index=True)
    return pd.DataFrame()


@addon("MRMTransitionGroupCP")
def to_feature_df(self, columns=None, meta_values=None):
    """Returns a DataFrame representation of the Features stored in MRMTransitionGroupCP."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for to_feature_df(). Install with: pip install pandas")

    features = self.feature_views()  # zero-copy read; the generators below only read
    str_dtype = string_dtype(len(features))
    # keyed by the decoded (str) name: getKeys() and user-supplied meta_values may be str or bytes
    common_meta_value_types = {
        'label': str_dtype, 'spectrum_index': 'i', 'score_fit': 'f',
        'score_correlation': 'f', 'FWHM': 'f', 'spectrum_native_id': str_dtype,
        'max_height': 'f', 'num_of_masstraces': 'i', 'masstrace_intensity': 'f',
        'Group': str_dtype, 'is_ungrouped_monoisotopic': 'i', 'leftWidth': 'f',
        'rightWidth': 'f', 'total_xic': 'f', 'PeptideRef': str_dtype,
        'peak_apices_sum': 'f'
    }

    def gen(features, fun):
        for f in features:
            yield from fun(f)

    def extract_meta_data(f):
        vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values_list]
        yield tuple((f.getUniqueId(), f.getRT(), f.getIntensity(), f.getOverallQuality(), *vals))

    mddtypes = [('feature_id', np.dtype('uint64')), ('rt', 'f'), ('intensity', 'f'), ('quality', 'f')]

    if meta_values is not None:
        if meta_values == 'all':
            meta_values_list = set()
            for f in features:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values_list.add(m)
            meta_values_list = list(meta_values_list)
        else:
            meta_values_list = list(meta_values)

        for meta_value in meta_values_list:
            name = meta_value.decode() if isinstance(meta_value, bytes) else meta_value
            dtype = common_meta_value_types.get(name, str_dtype)
            if dtype == 'i' and not all(f.metaValueExists(meta_value) for f in features):
                # an integer field cannot hold NaN for the missing values; promote to float64
                # (what pandas does for an integer column with missing entries)
                dtype = 'd'
            mddtypes.append((name, dtype))
    else:
        meta_values_list = []

    mdarr = np.fromiter(iter=gen(features, extract_meta_data), dtype=mddtypes, count=len(features))

    df = pin_string_dtype(pd.DataFrame(mdarr), mdarr).set_index('feature_id')

    if columns is not None:
        available_cols = [c for c in columns if c in df.columns or c == 'feature_id']
        if 'feature_id' not in available_cols:
            available_cols = [c for c in columns if c in df.columns]
        df = df[[c for c in available_cols if c in df.columns]]

    return df


@addon("MRMTransitionGroupCP")
def to_arrow(self, columns=None, export_meta_values=True):
    """Returns an Apache Arrow Table representation of the Chromatograms."""
    try:
        import pyarrow as pa
    except ImportError:
        raise ImportError("pyarrow is required for to_arrow(). Install with: pip install pyarrow")
    df = self.to_chromatogram_df(columns=columns, export_meta_values=export_meta_values)
    return pa.Table.from_pandas(df)


@addon("MRMTransitionGroupCP")
def get_chromatogram_df(self, *args, **kwargs):
    """Deprecated: Use to_chromatogram_df() instead."""
    import warnings
    warnings.warn("get_chromatogram_df() is deprecated. Use to_chromatogram_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_chromatogram_df(*args, **kwargs)


@addon("MRMTransitionGroupCP")
def get_chromatogram_df_columns(self, *args, **kwargs):
    """Deprecated: Use chromatogram_df_columns() instead."""
    import warnings
    warnings.warn("get_chromatogram_df_columns() is deprecated. Use chromatogram_df_columns() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.chromatogram_df_columns(*args, **kwargs)


@addon("MRMTransitionGroupCP")
def get_feature_df(self, *args, **kwargs):
    """Deprecated: Use to_feature_df() instead."""
    import warnings
    warnings.warn("get_feature_df() is deprecated. Use to_feature_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_feature_df(*args, **kwargs)


@addon("MRMTransitionGroupCP")
def get_feature_df_columns(self, *args, **kwargs):
    """Deprecated: Use feature_df_columns() instead."""
    import warnings
    warnings.warn("get_feature_df_columns() is deprecated. Use feature_df_columns() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.feature_df_columns(*args, **kwargs)


# The plural/iterator view families are generated from one template so the
# naming and contract wording cannot drift between them.
register_element_views("MRMTransitionGroupCP", "feature", "_feature_count", "features")
register_element_views("MRMTransitionGroupCP", "chromatogram", "_chromatogram_count", "chromatograms")
register_element_views("LightMRMTransitionGroupCP", "feature", "_feature_count", "features")
register_element_views("LightMRMTransitionGroupCP", "chromatogram", "_chromatogram_count", "chromatograms")


