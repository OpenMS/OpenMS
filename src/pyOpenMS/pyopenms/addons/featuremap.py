"""Addon methods for FeatureMap class."""

from __future__ import annotations
import warnings
import numpy as np
from . import addon, register_element_views


@addon("FeatureMap")
def df_columns(self, columns='default', export_peptide_identifications=True):
    """Returns a list of column names that to_df() would produce."""
    cols = ['feature_id']
    if export_peptide_identifications:
        cols.extend(['peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'])
    cols.extend(['charge', 'rt', 'mz', 'rt_start', 'rt_end', 'mz_start', 'mz_end', 'quality', 'intensity'])

    if columns == 'all':
        meta_values = set()
        for f in self.iter_feature_views():
            mvs = []
            f.getKeys(mvs)
            meta_values.update(mvs)
        cols.extend(sorted(meta_values))
    return cols


@addon("FeatureMap")
def _get_prot_id_filename_from_pep_id(self, pep_id):
    """Gets the primary MS run path of the ProteinIdentification linked to the PeptideIdentification."""
    for prot in self.getProteinIdentifications():
        if prot.getIdentifier() == pep_id.getIdentifier():
            filenames = prot.getPrimaryMSRunPath()
            if filenames and filenames[0] != '':
                return filenames[0]
    return None


@addon("FeatureMap")
def to_df(self, columns=None, meta_values=None, export_peptide_identifications=True):
    """Returns a pandas DataFrame with feature information."""
    import pandas as pd

    pep_id_cols = {'peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'}
    if columns is not None:
        requested = set(columns)
        need_pep_ids = export_peptide_identifications and len(requested & pep_id_cols) > 0
    else:
        need_pep_ids = export_peptide_identifications

    if meta_values == 'all':
        meta_values_set = set()
        for f in self.iter_feature_views():
            mvs = []
            f.getKeys(mvs)
            for m in mvs:
                meta_values_set.add(m)
        meta_values = list(meta_values_set)
    elif not meta_values:
        meta_values = []

    rows = []
    for f in self.iter_feature_views():
        # Compute bounding box from hull points
        hull_pts = f.getConvexHull().getHullPoints()
        if len(hull_pts) > 0:
            pts = np.array(hull_pts)
            rt_start, mz_start = pts.min(axis=0)
            rt_end, mz_end = pts.max(axis=0)
        else:
            rt_start = rt_end = f.getRT()
            mz_start = mz_end = f.getMZ()

        vals = []
        for m in meta_values:
            if f.metaValueExists(m):
                vals.append(f.getMetaValue(m))
            else:
                vals.append(np.nan)

        if need_pep_ids:
            pep = f.getPeptideIdentifications()
            if len(pep) > 0:
                ID_filename = self._get_prot_id_filename_from_pep_id(pep[0])
                spec_id = None
                if f.metaValueExists('spectrum_native_id'):
                    spec_id = str(f.getMetaValue('spectrum_native_id'))
                hits = pep[0].getHits()
                if len(hits) > 0:
                    besthit = hits[0]
                    pep_values = [besthit.getSequence().toString(), besthit.getScore(), ID_filename, spec_id]
                else:
                    pep_values = [None, None, ID_filename, spec_id]
            else:
                pep_values = [None, None, None, None]
        else:
            pep_values = []

        row = [f.getUniqueId()] + pep_values + [
            f.getCharge(), f.getRT(), f.getMZ(),
            rt_start, rt_end, mz_start, mz_end,
            f.getOverallQuality(), f.getIntensity()
        ] + vals
        rows.append(row)

    col_names = ['feature_id']
    if need_pep_ids:
        col_names += ['peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id']
    col_names += ['charge', 'rt', 'mz', 'rt_start', 'rt_end', 'mz_start', 'mz_end', 'quality', 'intensity']
    for m in meta_values:  # the caller's meta_values may be bytes
        col_names.append(m.decode() if isinstance(m, bytes) else m)

    df = pd.DataFrame(rows, columns=col_names)
    # uint64 whatever the values: pandas would infer int64 when every ID happened to fit
    df['feature_id'] = df['feature_id'].astype(np.uint64)
    df = df.set_index('feature_id')

    if columns is not None:
        available_cols = [c for c in columns if c in df.columns]
        df = df[available_cols]
    return df


@addon("FeatureMap")
def get_df(self, *args, **kwargs):
    """Deprecated: use to_df() instead."""
    warnings.warn("get_df() is deprecated. Use to_df() instead.",
                  DeprecationWarning, stacklevel=2)
    return self.to_df(*args, **kwargs).reset_index()


@addon("FeatureMap")
def get_df_columns(self, *args, **kwargs):
    """Deprecated: Use df_columns() instead."""
    warnings.warn(
        "get_df_columns() is deprecated. Use df_columns() instead.",
        DeprecationWarning, stacklevel=2
    )
    return self.df_columns(*args, **kwargs)


@addon("FeatureMap")
def get_assigned_peptide_identifications(self):
    """Returns all PeptideIdentifications assigned to features in this map.

    The identifications are returned as the features store them, feature by feature.
    The list holds copies, so neither the map nor its features change. To relate them
    to their features, use to_peptide_df(), which adds each one's feature_id.
    """
    from pyopenms._pyopenms_metadata import PeptideIdentificationList
    result = PeptideIdentificationList()
    for f in self.iter_feature_views():
        for pid in f.getPeptideIdentifications():
            result.push_back(pid)
    return result


@addon("FeatureMap")
def peptide_df_columns(self, decode_ontology=True):
    """Returns a list of column names that to_peptide_df() would produce."""
    peps = self.get_assigned_peptide_identifications()
    return ['feature_id'] + [c for c in peps.df_columns(decode_ontology=decode_ontology) if c != 'feature_id']


@addon("FeatureMap")
def to_peptide_df(self, decode_ontology=True, default_missing_values=None, export_unidentified=True,
                  columns=None):
    """Returns the PeptideIdentifications assigned to features as a pandas DataFrame.

    One row per identification, as PeptideIdentificationList.to_df() writes the list
    that get_assigned_peptide_identifications() returns, preceded by a 'feature_id'
    column: the unique ID of the identification's feature, as the unsigned 64-bit
    integer that also indexes to_df(). 'P_ID' is the identification's position in
    get_assigned_peptide_identifications(). Merge the two frames with::

        merged = pd.merge(fmap.to_df().reset_index(),
                          fmap.to_peptide_df(export_unidentified=False),
                          on='feature_id', suffixes=('', '_psm'))

    The merge needs unique feature IDs. Features without one, such as features created
    in Python, have the ID 0 until FeatureMap.setUniqueIds() assigns them one. A
    'feature_id' meta value on the hits, as 3.5.0's get_assigned_peptide_identifications()
    added it, gives way to the 'feature_id' column.

    :param decode_ontology: Decode meta value names using the PSI-MS ontology.
    :param default_missing_values: Default values for missing data by type.
    :param export_unidentified: Export PeptideIdentifications without PeptideHit.
    :param columns: Columns to include after 'feature_id', which is always the first.
        If None, includes all.
    :return: DataFrame with one row per assigned peptide identification.
    """
    from pyopenms._pyopenms_metadata import PeptideIdentificationList
    peps = PeptideIdentificationList()
    feature_ids = []
    for f in self.iter_feature_views():
        feature_id = f.getUniqueId()
        for pid in f.getPeptideIdentifications():
            peps.push_back(pid)
            # to_df() writes a row for pid only if this holds
            if export_unidentified or pid.getHits():
                feature_ids.append(feature_id)
    if columns is not None:
        columns = [c for c in columns if c != 'feature_id']
    df = peps.to_df(decode_ontology=decode_ontology, default_missing_values=default_missing_values,
                    export_unidentified=export_unidentified, columns=columns)
    # a 'feature_id' meta value of the hits gives way to the feature's ID
    df = df.drop(columns='feature_id', errors='ignore')
    df.insert(0, 'feature_id', np.array(feature_ids, dtype=np.uint64))
    return df


@addon("FeatureMap")
def to_arrow(self, columns=None, meta_values=None, export_peptide_identifications=True):
    """Returns an Apache Arrow Table with feature information."""
    import pyarrow as pa
    df = self.to_df(columns=columns, meta_values=meta_values,
                    export_peptide_identifications=export_peptide_identifications)
    return pa.Table.from_pandas(df)


# The plural/iterator view families are generated from one template so the
# naming and contract wording cannot drift between them.
register_element_views("FeatureMap", "feature", "size", "features")


