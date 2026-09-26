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

    df = pd.DataFrame(rows, columns=col_names).set_index('feature_id')

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

    Every peptide hit carries the keys for merging with to_df() as meta values:
    'feature_id' (the unique ID of its feature), 'ID_native_id' (the feature's
    spectrum_native_id) and 'ID_filename' (the primary MS run path of the matching
    ProteinIdentification). A value that is not known, null in to_df(), is not set.
    'feature_id' is text, because a meta value cannot hold an unsigned 64-bit integer,
    so convert it before merging::

        peps = fmap.get_assigned_peptide_identifications()
        pep_df = peps.to_df(export_unidentified=False)
        pep_df['feature_id'] = pep_df['feature_id'].astype('uint64')
        merged = pd.merge(fmap.to_df().reset_index(), pep_df, on='feature_id')

    A PeptideIdentification without hits has no hit to carry the keys, hence
    export_unidentified=False. The feature map itself is not modified.
    """
    from pyopenms._pyopenms_metadata import PeptideIdentificationList
    result = PeptideIdentificationList()
    for f in self.iter_feature_views():
        feature_id = str(f.getUniqueId())
        native_id = None
        if f.metaValueExists('spectrum_native_id'):
            native_id = str(f.getMetaValue('spectrum_native_id'))
        for pep in f.getPeptideIdentifications():
            filename = self._get_prot_id_filename_from_pep_id(pep)
            hits = pep.getHits()
            for hit in hits:
                hit.setMetaValue('feature_id', feature_id)
                if native_id is not None:
                    hit.setMetaValue('ID_native_id', native_id)
                if filename is not None:
                    hit.setMetaValue('ID_filename', filename)
            pep.setHits(hits)
            result.push_back(pep)
    return result


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


