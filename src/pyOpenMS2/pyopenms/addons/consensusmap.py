"""ConsensusMap addon methods for DataFrame support."""
import numpy as np
from collections import defaultdict as _defaultdict
from . import addon


@addon("ConsensusMap")
def __repr__(self) -> str:
    return f"ConsensusMap(num_consensus_features={len(self)})"


@addon("ConsensusMap")
def df_columns(self, columns='default'):
    """Returns a list of column names that to_df() would produce."""
    cols = ['sequence', 'charge', 'rt', 'mz', 'quality']

    labelfree = self.getExperimentType() == "label-free"
    filemeta = self.getColumnHeaders()

    if labelfree:
        files = list(set([header.filename for header in filemeta.values()]))
        cols.extend(files)
    else:
        labels = list(set([header.label for header in filemeta.values()]))
        if len(labels) == 1:
            labels[0] = "intensity"
        cols.extend(labels)
        cols.append('file')

    return cols


@addon("ConsensusMap")
def get_intensity_df(self):
    """Generates a pandas DataFrame with feature intensities from each sample in long format."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for get_intensity_df(). Install with: pip install pandas")

    labelfree = self.getExperimentType() == "label-free"
    filemeta = self.getColumnHeaders()

    labels = list(set([header.label for header in filemeta.values()]))
    files = list(set([header.filename for header in filemeta.values()]))
    label_to_idx = {k: v for v, k in enumerate(labels)}
    file_to_idx = {k: v for v, k in enumerate(files)}

    def gen(cmap, fun):
        for f in cmap:
            yield from fun(f)

    if not labelfree:
        def extract_row_blocks_channel_wide_file_long(f):
            subfeatures = f.getFeatureList()
            filerows = _defaultdict(lambda: [0] * len(labels))
            for fh in subfeatures:
                header = filemeta[fh.getMapIndex()]
                row = filerows[header.filename]
                row[label_to_idx[header.label]] = fh.getIntensity()
            return (f.getUniqueId(), filerows)

        def extract_rows_channel_wide_file_long(f):
            uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
            for file, row in rowdict.items():
                row.append(file)
                yield tuple([uniqueid] + row)

        if len(labels) == 1:
            labels[0] = "intensity"

        dtypes = [('id', np.dtype('uint64'))] + list(zip(labels, ['f'] * len(labels)))
        dtypes.append(('file', 'U300'))

        total_rows = sum(len(extract_row_blocks_channel_wide_file_long(f)[1]) for f in self)
        intyarr = np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=total_rows)

        return pd.DataFrame(intyarr).set_index('id')
    else:
        def extract_row_blocks_channel_long_file_wide_LF(f):
            subfeatures = f.getFeatureList()
            row = [0.] * len(files)
            for fh in subfeatures:
                header = filemeta[fh.getMapIndex()]
                row[file_to_idx[header.filename]] = fh.getIntensity()
            yield tuple([f.getUniqueId()] + row)

        dtypes = [('id', np.dtype('uint64'))] + list(zip(files, ['f'] * len(files)))
        intyarr = np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

        return pd.DataFrame(intyarr).set_index('id')


@addon("ConsensusMap")
def get_metadata_df(self):
    """Generates a pandas DataFrame with feature meta data."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for get_metadata_df(). Install with: pip install pandas")

    def gen(cmap, fun):
        for f in cmap:
            yield from fun(f)

    def extract_meta_data(f):
        pep = f.getPeptideIdentifications()
        if len(pep) != 0:
            hits = pep[0].getHits()
            if len(hits) != 0:
                besthit = hits[0]
                yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
        else:
            yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

    cnt = self.size()
    mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'),
                ('rt', np.dtype('double')), ('mz', np.dtype('double')), ('quality', 'f')]
    mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

    return pd.DataFrame(mdarr).set_index('id')


@addon("ConsensusMap")
def to_df(self, columns=None):
    """Generates a pandas DataFrame with both consensus feature meta data and intensities."""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for to_df(). Install with: pip install pandas")

    if columns is None:
        df = pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)
        return df

    requested = set(columns)
    metadata_cols = {'sequence', 'charge', 'rt', 'mz', 'quality'}

    labelfree = self.getExperimentType() == "label-free"
    filemeta = self.getColumnHeaders()
    if labelfree:
        intensity_cols = set([header.filename for header in filemeta.values()])
    else:
        labels_list = list(set([header.label for header in filemeta.values()]))
        if len(labels_list) == 1:
            labels_list[0] = "intensity"
        intensity_cols = set(labels_list)
        intensity_cols.add('file')

    need_metadata = len(requested & metadata_cols) > 0
    need_intensity = len(requested & intensity_cols) > 0

    dfs = []
    if need_metadata:
        dfs.append(self.get_metadata_df())
    if need_intensity:
        dfs.append(self.get_intensity_df())

    if not dfs:
        return pd.DataFrame(index=pd.Index([], name='id'))

    if len(dfs) == 1:
        df = dfs[0]
    else:
        df = pd.concat(dfs, axis=1)

    available_cols = [c for c in columns if c in df.columns]
    return df[available_cols]


@addon("ConsensusMap")
def to_arrow(self, columns=None):
    """Returns an Apache Arrow Table with consensus feature meta data and intensities."""
    try:
        import pyarrow as pa
    except ImportError:
        raise ImportError("pyarrow is required for to_arrow(). Install with: pip install pyarrow")
    df = self.to_df(columns=columns)
    return pa.Table.from_pandas(df)


@addon("ConsensusMap")
def get_df(self, *args, **kwargs):
    """Deprecated: Use to_df() instead."""
    import warnings
    warnings.warn(
        "get_df() is deprecated. Use to_df() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.to_df(*args, **kwargs)


@addon("ConsensusMap")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    import warnings
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)
