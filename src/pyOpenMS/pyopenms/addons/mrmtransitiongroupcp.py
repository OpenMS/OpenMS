"""MRMTransitionGroupCP addon methods for DataFrame support."""
import numpy as np
from . import addon


@addon("MRMTransitionGroupCP")
def chromatogram_df_columns(self, columns='default', export_meta_values=True):
    """Returns a list of column names that to_chromatogram_df() would produce."""
    chroms = self.getChromatograms()
    if chroms:
        return chroms[0].df_columns(columns=columns, export_meta_values=export_meta_values)
    return ['rt', 'intensity', 'precursor_mz', 'precursor_charge', 'product_mz', 'native_id']


@addon("MRMTransitionGroupCP")
def feature_df_columns(self, columns='default'):
    """Returns a list of column names that to_feature_df() would produce."""
    cols = ['feature_id', 'rt', 'intensity', 'quality']

    if columns == 'all':
        meta_values = set()
        for f in self.getFeatures():
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
    chroms = self.getChromatograms()
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

    common_meta_value_types = {
        b'label': 'U50', b'spectrum_index': 'i', b'score_fit': 'f',
        b'score_correlation': 'f', b'FWHM': 'f', b'spectrum_native_id': 'U100',
        b'max_height': 'f', b'num_of_masstraces': 'i', b'masstrace_intensity': 'f',
        b'Group': 'U50', b'is_ungrouped_monoisotopic': 'i', b'leftWidth': 'f',
        b'rightWidth': 'f', b'total_xic': 'f', b'PeptideRef': 'U100',
        b'peak_apices_sum': 'f'
    }

    def gen(features, fun):
        for f in features:
            yield from fun(f)

    def extract_meta_data(f):
        vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values_list]
        yield tuple((f.getUniqueId(), f.getRT(), f.getIntensity(), f.getOverallQuality(), *vals))

    features = self.getFeatures()
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
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode() if isinstance(meta_value, bytes) else meta_value,
                                common_meta_value_types.get(meta_value, 'U50')))
            else:
                mddtypes.append((meta_value.decode() if isinstance(meta_value, bytes) else meta_value, 'U50'))
    else:
        meta_values_list = []

    mdarr = np.fromiter(iter=gen(features, extract_meta_data), dtype=mddtypes, count=len(features))

    df = pd.DataFrame(mdarr).set_index('feature_id')

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
