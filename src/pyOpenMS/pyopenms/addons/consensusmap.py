"""ConsensusMap addon methods for DataFrame support."""
import numpy as np
from collections import defaultdict as _defaultdict
from . import addon

try:
    from pyopenms._arrow_zerocopy import ConsensusFeatureExportSchema as _CFES
    _has_schema_constants = True
except ImportError:
    _has_schema_constants = False


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


@addon("ConsensusMap")
def feature_columns(self):
    """
    Return list of column names that to_feature_arrow() and to_feature_qpx() produce.

    Returns
    -------
    list
        List of column name strings.
    """
    if _has_schema_constants:
        return [
            # QPX feature schema fields
            _CFES.SEQUENCE,
            _CFES.PEPTIDOFORM,
            _CFES.MODIFICATIONS,
            _CFES.PRECURSOR_CHARGE,
            _CFES.CALCULATED_MZ,
            _CFES.OBSERVED_MZ,
            _CFES.RT,
            _CFES.POSTERIOR_ERROR_PROBABILITY,
            _CFES.IS_DECOY,
            _CFES.ADDITIONAL_SCORES,
            _CFES.PREDICTED_RT,
            _CFES.REFERENCE_FILE_NAME,
            _CFES.CV_PARAMS,
            _CFES.SCAN,
            _CFES.ION_MOBILITY,
            _CFES.START_ION_MOBILITY,
            _CFES.STOP_ION_MOBILITY,
            _CFES.INTENSITIES,
            _CFES.ADDITIONAL_INTENSITIES,
            _CFES.PG_ACCESSIONS,
            _CFES.ANCHOR_PROTEIN,
            _CFES.UNIQUE,
            _CFES.PG_GLOBAL_QVALUE,
            _CFES.GG_ACCESSIONS,
            _CFES.GG_NAMES,
            _CFES.SCAN_REFERENCE_FILE_NAME,
            _CFES.RT_START,
            _CFES.RT_STOP,
            # OpenMS-specific fields
            _CFES.QUALITY,
            _CFES.SCORE,
            _CFES.SCORE_TYPE,
            _CFES.SPECTRUM_REFERENCE,
            _CFES.FEATURE_METAVALUES,
        ]
    return [
        # QPX feature schema fields
        "sequence",
        "peptidoform",
        "modifications",
        "precursor_charge",
        "calculated_mz",
        "observed_mz",
        "rt",
        "posterior_error_probability",
        "is_decoy",
        "additional_scores",
        "predicted_rt",
        "reference_file_name",
        "cv_params",
        "scan",
        "ion_mobility",
        "start_ion_mobility",
        "stop_ion_mobility",
        "intensities",
        "additional_intensities",
        "pg_accessions",
        "anchor_protein",
        "unique",
        "pg_global_qvalue",
        "gg_accessions",
        "gg_names",
        "scan_reference_file_name",
        "rt_start",
        "rt_stop",
        # OpenMS-specific fields
        "quality",
        "score",
        "score_type",
        "spectrum_reference",
        "feature_metavalues",
    ]


@addon("ConsensusMap")
def to_feature_arrow(self, reference_file_name=None, columns=None,
                     include_modifications=True, scan_format='native_id'):
    """
    Export consensus features as Apache Arrow Table following QPX feature schema.

    Parameters
    ----------
    reference_file_name : str, optional
        Reference file name to use for unmatched features.
    columns : list, optional
        Columns to include. If None, includes all.
    include_modifications : bool
        Include modification information (default: True).
    scan_format : str
        Format for scan numbers (default: 'native_id').

    Returns
    -------
    pyarrow.Table
        Arrow table with feature data.
    """
    try:
        import pyarrow as pa
    except ImportError:
        raise ImportError("pyarrow is required for to_feature_arrow(). Install with: pip install pyarrow")

    import pyopenms

    # Column name constants (use schema constants when available, else strings)
    if _has_schema_constants:
        _c = _CFES
        _SEQUENCE = _c.SEQUENCE
        _PEPTIDOFORM = _c.PEPTIDOFORM
        _MODIFICATIONS = _c.MODIFICATIONS
        _PRECURSOR_CHARGE = _c.PRECURSOR_CHARGE
        _CALCULATED_MZ = _c.CALCULATED_MZ
        _OBSERVED_MZ = _c.OBSERVED_MZ
        _RT = _c.RT
        _PEP = _c.POSTERIOR_ERROR_PROBABILITY
        _IS_DECOY = _c.IS_DECOY
        _ADDITIONAL_SCORES = _c.ADDITIONAL_SCORES
        _PREDICTED_RT = _c.PREDICTED_RT
        _REFERENCE_FILE_NAME = _c.REFERENCE_FILE_NAME
        _CV_PARAMS = _c.CV_PARAMS
        _SCAN = _c.SCAN
        _ION_MOBILITY = _c.ION_MOBILITY
        _START_ION_MOBILITY = _c.START_ION_MOBILITY
        _STOP_ION_MOBILITY = _c.STOP_ION_MOBILITY
        _INTENSITIES = _c.INTENSITIES
        _ADDITIONAL_INTENSITIES = _c.ADDITIONAL_INTENSITIES
        _PG_ACCESSIONS = _c.PG_ACCESSIONS
        _ANCHOR_PROTEIN = _c.ANCHOR_PROTEIN
        _UNIQUE = _c.UNIQUE
        _PG_GLOBAL_QVALUE = _c.PG_GLOBAL_QVALUE
        _GG_ACCESSIONS = _c.GG_ACCESSIONS
        _GG_NAMES = _c.GG_NAMES
        _SCAN_REFERENCE_FILE_NAME = _c.SCAN_REFERENCE_FILE_NAME
        _RT_START = _c.RT_START
        _RT_STOP = _c.RT_STOP
        _QUALITY = _c.QUALITY
        _SCORE = _c.SCORE
        _SCORE_TYPE = _c.SCORE_TYPE
        _SPECTRUM_REFERENCE = _c.SPECTRUM_REFERENCE
        _FEATURE_METAVALUES = _c.FEATURE_METAVALUES
    else:
        _SEQUENCE = 'sequence'
        _PEPTIDOFORM = 'peptidoform'
        _MODIFICATIONS = 'modifications'
        _PRECURSOR_CHARGE = 'precursor_charge'
        _CALCULATED_MZ = 'calculated_mz'
        _OBSERVED_MZ = 'observed_mz'
        _RT = 'rt'
        _PEP = 'posterior_error_probability'
        _IS_DECOY = 'is_decoy'
        _ADDITIONAL_SCORES = 'additional_scores'
        _PREDICTED_RT = 'predicted_rt'
        _REFERENCE_FILE_NAME = 'reference_file_name'
        _CV_PARAMS = 'cv_params'
        _SCAN = 'scan'
        _ION_MOBILITY = 'ion_mobility'
        _START_ION_MOBILITY = 'start_ion_mobility'
        _STOP_ION_MOBILITY = 'stop_ion_mobility'
        _INTENSITIES = 'intensities'
        _ADDITIONAL_INTENSITIES = 'additional_intensities'
        _PG_ACCESSIONS = 'pg_accessions'
        _ANCHOR_PROTEIN = 'anchor_protein'
        _UNIQUE = 'unique'
        _PG_GLOBAL_QVALUE = 'pg_global_qvalue'
        _GG_ACCESSIONS = 'gg_accessions'
        _GG_NAMES = 'gg_names'
        _SCAN_REFERENCE_FILE_NAME = 'scan_reference_file_name'
        _RT_START = 'rt_start'
        _RT_STOP = 'rt_stop'
        _QUALITY = 'quality'
        _SCORE = 'score'
        _SCORE_TYPE = 'score_type'
        _SPECTRUM_REFERENCE = 'spectrum_reference'
        _FEATURE_METAVALUES = 'feature_metavalues'

    # Build data arrays
    data = {
        _SEQUENCE: [],
        _PEPTIDOFORM: [],
        _MODIFICATIONS: [],
        _PRECURSOR_CHARGE: [],
        _CALCULATED_MZ: [],
        _OBSERVED_MZ: [],
        _RT: [],
        _PEP: [],
        _IS_DECOY: [],
        _ADDITIONAL_SCORES: [],
        _PREDICTED_RT: [],
        _REFERENCE_FILE_NAME: [],
        _CV_PARAMS: [],
        _SCAN: [],
        _ION_MOBILITY: [],
        _START_ION_MOBILITY: [],
        _STOP_ION_MOBILITY: [],
        _INTENSITIES: [],
        _ADDITIONAL_INTENSITIES: [],
        _PG_ACCESSIONS: [],
        _ANCHOR_PROTEIN: [],
        _UNIQUE: [],
        _PG_GLOBAL_QVALUE: [],
        _GG_ACCESSIONS: [],
        _GG_NAMES: [],
        _SCAN_REFERENCE_FILE_NAME: [],
        _RT_START: [],
        _RT_STOP: [],
        _QUALITY: [],
        _SCORE: [],
        _SCORE_TYPE: [],
        _SPECTRUM_REFERENCE: [],
        _FEATURE_METAVALUES: [],
    }

    # Get column headers for intensity mapping
    col_headers = self.getColumnHeaders()
    prot_ids = self.getProteinIdentifications()

    # Build protein group lookup
    pg_lookup = {}
    pg_qvalue_lookup = {}
    if prot_ids:
        for prot_id in prot_ids:
            for pg in prot_id.getProteinGroups():
                for acc in pg.accessions:
                    acc_str = acc.decode() if isinstance(acc, bytes) else str(acc)
                    if acc_str not in pg_lookup:
                        pg_lookup[acc_str] = []
                    pg_lookup[acc_str].extend([a.decode() if isinstance(a, bytes) else str(a) for a in pg.accessions])
                    pg_qvalue_lookup[acc_str] = pg.probability

    for cf in self:
        pep_ids = cf.getPeptideIdentifications()
        best_hit = None
        if pep_ids:
            hits = pep_ids[0].getHits()
            if hits:
                best_hit = hits[0]

        # Sequence info
        if best_hit:
            seq = best_hit.getSequence()
            data[_SEQUENCE].append(seq.toUnmodifiedString())
            data[_PEPTIDOFORM].append(seq.toString())
            data[_MODIFICATIONS].append([])  # Simplified
            data[_CALCULATED_MZ].append(seq.getMZ(cf.getCharge()) if cf.getCharge() > 0 else None)
            data[_SCORE].append(best_hit.getScore())
            data[_SCORE_TYPE].append(pep_ids[0].getScoreType() if pep_ids else '')
            # Get protein accessions
            prot_accs = [ev.getProteinAccession() for ev in best_hit.getPeptideEvidences()]
            data[_PG_ACCESSIONS].append(prot_accs)
            data[_ANCHOR_PROTEIN].append(prot_accs[0] if prot_accs else None)
            # unique: 1 if single protein, 0 otherwise (QPX schema uses int)
            data[_UNIQUE].append(1 if len(set(prot_accs)) == 1 else 0)
            # Get q-value from protein groups
            qval = None
            for acc in prot_accs:
                if acc in pg_qvalue_lookup:
                    qval = pg_qvalue_lookup[acc]
                    break
            data[_PG_GLOBAL_QVALUE].append(qval)
            # Decoy status from metavalue (QPX schema uses int: 0=target, 1=decoy)
            if best_hit.metaValueExists('target_decoy'):
                td_val = best_hit.getMetaValue('target_decoy')
                is_decoy = 1 if td_val == 'decoy' else 0
            else:
                is_decoy = None
            data[_IS_DECOY].append(is_decoy)
        else:
            data[_SEQUENCE].append(None)
            data[_PEPTIDOFORM].append(None)
            data[_MODIFICATIONS].append(None)
            data[_CALCULATED_MZ].append(None)
            data[_SCORE].append(None)
            data[_SCORE_TYPE].append(None)
            data[_PG_ACCESSIONS].append(None)
            data[_ANCHOR_PROTEIN].append(None)
            data[_UNIQUE].append(0)  # QPX: 0 for unidentified
            data[_PG_GLOBAL_QVALUE].append(None)
            data[_IS_DECOY].append(None)

        data[_PRECURSOR_CHARGE].append(cf.getCharge())
        data[_OBSERVED_MZ].append(cf.getMZ())
        data[_RT].append(cf.getRT())
        data[_PEP].append(None)
        data[_ADDITIONAL_SCORES].append(None)
        data[_PREDICTED_RT].append(None)
        data[_REFERENCE_FILE_NAME].append(reference_file_name)
        data[_CV_PARAMS].append(None)
        data[_SCAN].append(None)
        data[_ION_MOBILITY].append(None)
        data[_START_ION_MOBILITY].append(None)
        data[_STOP_ION_MOBILITY].append(None)

        # Build intensities as list of dicts (QPX schema)
        intensities = []
        labelfree = self.getExperimentType() == "label-free"
        for fh in cf.getFeatureList():
            map_idx = fh.getMapIndex()
            if map_idx in col_headers:
                header = col_headers[map_idx]
                entry = {
                    "sample_accession": header.filename if header.filename else f"sample_{map_idx}",
                    "intensity": fh.getIntensity()
                }
                if not labelfree and header.label:
                    entry["channel"] = header.label
                intensities.append(entry)
        data[_INTENSITIES].append(intensities)
        data[_ADDITIONAL_INTENSITIES].append(None)

        data[_GG_ACCESSIONS].append(None)
        data[_GG_NAMES].append(None)
        data[_SCAN_REFERENCE_FILE_NAME].append(None)
        data[_RT_START].append(None)
        data[_RT_STOP].append(None)
        data[_QUALITY].append(cf.getQuality())
        data[_SPECTRUM_REFERENCE].append(None)
        data[_FEATURE_METAVALUES].append(None)

    # Filter columns if specified
    if columns:
        data = {k: v for k, v in data.items() if k in columns}

    table = pa.table(data)

    # Sort by RT for consistent ordering
    if _RT in table.column_names and table.num_rows > 0:
        import pyarrow.compute as pc
        indices = pc.sort_indices(table, sort_keys=[(_RT, 'ascending')])
        table = table.take(indices)

    return table


@addon("ConsensusMap")
def to_feature_df(self, **kwargs):
    """
    Export consensus features as pandas DataFrame following QPX feature schema.

    Parameters
    ----------
    **kwargs
        Passed to to_feature_arrow().

    Returns
    -------
    pandas.DataFrame
        DataFrame with feature data.
    """
    return self.to_feature_arrow(**kwargs).to_pandas()


@addon("ConsensusMap")
def to_feature_parquet(self, filename, compression='snappy', **kwargs):
    """
    Export consensus features to a Parquet file.

    Parameters
    ----------
    filename : str
        Output file path.
    compression : str
        Compression codec (default: 'snappy').
    **kwargs
        Passed to to_feature_arrow().
    """
    try:
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError("pyarrow is required for to_feature_parquet(). Install with: pip install pyarrow")

    table = self.to_feature_arrow(**kwargs)
    pq.write_table(table, filename, compression=compression)


@addon("ConsensusMap")
def to_feature_qpx(self, qpx_version="1.0", creator="pyopenms", software_provider="OpenMS",
                   **kwargs):
    """
    Export consensus features as QPX format (dict with file_metadata and features).

    Parameters
    ----------
    qpx_version : str
        QPX schema version.
    creator : str
        Creator identifier.
    software_provider : str
        Software provider name.
    **kwargs
        Passed to to_feature_arrow().

    Returns
    -------
    dict
        Dictionary with 'file_metadata' and 'features' keys.
    """
    import uuid
    from datetime import datetime

    table = self.to_feature_arrow(**kwargs)

    file_metadata = {
        "qpx_version": qpx_version,
        "creator": creator,
        "file_type": "feature",
        "creation_date": datetime.now().isoformat(),
        "uuid": str(uuid.uuid4()),
        "software_provider": software_provider,
    }

    # Convert to list of records (row-oriented) like legacy pandas to_dict(orient="records")
    return {
        "file_metadata": file_metadata,
        "features": table.to_pylist()
    }
