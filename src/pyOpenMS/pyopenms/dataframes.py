from collections import defaultdict
from typing import List

from . import ConsensusMap, ConsensusFeature, FeatureMap, Feature, MSExperiment, PeakMap, PeptideIdentification, ControlledVocabulary, File

import pandas as pd
import numpy as np

class ConsensusMapDF(ConsensusMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_intensity_df(self):
        """Generates a pandas DataFrame with feature intensities from each sample in long format (over files).
        For labelled analyses channel intensities will be in one row, therefore resulting in a semi-long/block format.
        Resulting DataFrame can be joined with result from get_metadata_df by their index 'id'.

        Returns:
        pandas.DataFrame: intensity DataFrame
        """
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()  # type: dict[int, ColumnHeader]

        labels = list(set([header.label for header in filemeta.values()]))
        files = list(set([header.filename for header in filemeta.values()]))
        label_to_idx = {k: v for v, k in enumerate(labels)}
        file_to_idx = {k: v for v, k in enumerate(files)}

        def gen(cmap: ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        if not labelfree:

            def extract_row_blocks_channel_wide_file_long(f: ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                filerows = defaultdict(lambda: [0] * len(labels))
                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row = filerows[header.filename]
                    row[label_to_idx[header.label]] = fh.getIntensity()
                return (f.getUniqueId(), filerows)

            def extract_rows_channel_wide_file_long(f: ConsensusFeature):
                uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
                for file, row in rowdict.items():
                    row.append(file)
                    yield tuple([uniqueid] + row)

            if len(labels) == 1:
                labels[0] = "intensity"

            dtypes = [('id', np.dtype('uint64'))] + list(zip(labels, ['f'] * len(labels)))
            dtypes.append(('file', 'U300'))

            intyarr = np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=self.size())

            return pd.DataFrame(intyarr).set_index('id')

        else:
            # Specialized for LabelFree which has to have only one channel
            def extract_row_blocks_channel_long_file_wide_LF(f: ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                row = [0.] * len(files)

                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row[file_to_idx[header.filename]] = fh.getIntensity()

                yield tuple([f.getUniqueId()] + row)

            dtypes = [('id', np.dtype('uint64'))] + list(zip(files, ['f'] * len(files)))

            intyarr = np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

            return pd.DataFrame(intyarr).set_index('id')

    def get_metadata_df(self):
        """Generates a pandas DataFrame with feature meta data (sequence, charge, mz, RT, quality).
        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        Returns:
        pandas.DataFrame: DataFrame with metadata for each feature (such as: best identified sequence, charge, centroid RT/mz, fitting quality)
        """

        def gen(cmap: ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        def extract_meta_data(f: ConsensusFeature):
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]

            if len(pep) != 0:
                hits = pep[0].getHits()

                if len(hits) != 0:
                    besthit = hits[0]  # type: PeptideHit
                    yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
                
                else:
                    yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()
            
            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), f.getQuality()

        cnt = self.size()

        mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'),
                    ('RT', np.dtype('double')), ('mz', np.dtype('double')), ('quality', 'f')]

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return pd.DataFrame(mdarr).set_index('id')
    
    def get_df(self):
        """Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        Returns:
        pandas.DataFrame: meta data and intensity DataFrame
        """
        return pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)

ConsensusMap = ConsensusMapDF

# TODO tell the advanced user that they could change this, in case they have different needs.
# TODO check if type could be inferred in the first pass
# TODO check if max. string lengths could be inferred in the first pass and how this affects runtime
# TODO check how runtime is affected if we use np.append instead of np.fromiter and use np.dyte = object for strings
common_meta_value_types = {
    b'label': 'U50',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U100',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f', # TODO this is actually a DoubleList. Think about what to do here. For np.fromiter we would need to set the length of the array.
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i' # TODO this sounds very boolean to me
}

class FeatureMapDF(FeatureMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    # meta_values = None (default), 'all' or list of meta value names
    def get_df(self, meta_values = None):
        """Generates a pandas DataFrame with information contained in the FeatureMap.

        Parameters:
        meta_values: meta values to include (None, [custom list of meta value names] or 'all')
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        # get all possible meta value keys in a set
        if meta_values == 'all':
            meta_values = set()
            for f in self:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m)

        elif not meta_values: # if None, set to empty list
            meta_values = []
        
        def gen(fmap: FeatureMap, fun):
            for f in fmap:
                yield from fun(f, meta_values)

        def extract_meta_data(f: Feature, meta_values):
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]
            bb = f.getConvexHull().getBoundingBox2D()
                
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.nan for m in meta_values]
            
            if len(pep) != 0:
                hits = pep[0].getHits()

                if len(hits) != 0:
                    besthit = hits[0]  # type: PeptideHit
                    yield (f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)
                else:
                    yield (f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)
            else:
                yield (f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)

        cnt = self.size()

        mddtypes = [('id', np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'), ('RT', np.dtype('double')), 
                    ('mz', np.dtype('double')), ('RTstart', np.dtype('double')), ('RTend', np.dtype('double')),
                    ('mzstart', np.dtype('double')), ('mzend', np.dtype('double')), ('quality', 'f'), ('intensity', 'f')]
        
        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return pd.DataFrame(mdarr).set_index('id')

FeatureMap = FeatureMapDF


class MSExperimentDF(MSExperiment):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, long : bool = False):
        """Generates a pandas DataFrame with all peaks in the MSExperiment

        Parameters:
        long: set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
            replicated RT information. If False, returns rows in the style: rt, np.array(mz), np.array(int)
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        if long:
            cols = ["RT", "mz", "inty"]
            self.updateRanges()
            spectraarrs2d = self.get2DPeakDataLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ())
            return pd.DataFrame(dict(zip(cols, spectraarrs2d)))

        cols = ["RT", "mzarray", "intarray"]

        return pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in self), columns=cols)

MSExperiment = MSExperimentDF
PeakMap = MSExperimentDF

# TODO think about the best way for such top-level function. IMHO in python, encapsulation in a stateless class in unnecessary.
#   We should probably not just import this whole submodule without prefix.
def peptide_identifications_to_df(peps: List[PeptideIdentification], decode_ontology : bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.

    Parameters:
    peps (List[PeptideIdentification]): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names

    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
    """
    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}
    metavals = []
    types = []
    mainscorename = "score"
    for pep in peps:
        hits = pep.getHits()
        if not len(hits) == 0:
            hits[0].getKeys(metavals)
            mainscorename = pep.getScoreType()
            for k in metavals:
                if k == b"target_decoy":
                    types.append('?')
                else:
                    mv = hits[0].getMetaValue(k)
                    types.append(switchDict[type(mv)])
            break

    decodedMVs = [m.decode("utf-8") for m in metavals] if decode_ontology else metavals
    cv = ControlledVocabulary()
    cv.loadFromOBO("psims", File.getOpenMSDataPath() + "/CV/psi-ms.obo")
    clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
    #cols = ["id", "RT", "mz", "score", "charge"] + decodedMVs
    clearcols = ["id", "RT", "mz", mainscorename, "charge"] + clearMVs
    coltypes = ['U100', 'f', 'f', 'f', 'i'] + types
    dt = list(zip(clearcols, coltypes))
    def extract(pep):
        hits = pep.getHits()
        if not hits:
            # default missing int value: -9999, replace with pandas.NA later
            return tuple([pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), np.nan, -9999] + [np.nan]*len(metavals))

        besthit = hits[0]
        ret = [pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()]
        for k in metavals:
            if besthit.metaValueExists(k):
                val = besthit.getMetaValue(k)
                if k == b"target_decoy":
                    if val[0] == 't':
                        ret.append(True)
                    else:
                        ret.append(False)
                else:
                    ret.append(val)
            else:
                ret.append(np.nan)
        return tuple(ret)

    psmarr = np.fromiter((extract(pep) for pep in peps), dtype=dt, count=len(peps))

    return pd.DataFrame(psmarr).replace(-9999, pd.NA)
