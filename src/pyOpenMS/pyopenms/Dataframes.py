from collections import defaultdict
from typing import List

from pyopenms import ConsensusMap, ConsensusFeature, FeatureMap, Feature, MSExperiment, PeptideIdentification, ControlledVocabulary, File

import pandas as pd
import numpy as np

class ConsensusMapDF(ConsensusMap):
    def __init__(self):
        super().__init__()

    def get_intensity_df(self):
        '''Generates a pandas DataFrame with feature intensities from each sample.

        Returns:
        pandas.DataFrame: intensity DataFrame
        '''
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

            # TODO write two functions for LF and labelled. One has only one channel, the other has only one file per CF
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
        '''Generates a pandas DataFrame with feature meta data (sequence, charge, mz, RT, quality).

        Returns:
        pandas.DataFrame: meta data DataFrame
        '''

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
        '''Generates a pandas DataFrame with feature meta data and intensities from each sample.

        Returns:
        pandas.DataFrame: meta data and intensity DataFrame
        '''
        return pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)

ConsensusMap = ConsensusMapDF

common_meta_value_types = {
    b'label': 'U30',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U30',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'i',
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i'
}

class FeatureMapDF(FeatureMap):
    def __init__(self):
        super().__init__()
    
    # meta_values = None (default), 'all' or list of meta value names
    def get_df(self, meta_values = None):
        '''Generates a pandas DataFrame with information contained in the FeatureMap.

        Parameters:
        meta_values: meta values to include (None, [custom list of meta value names] or 'all')
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        '''
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
                
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else np.NA for m in meta_values]   # find some NA or None value for numpy
            
            if len(pep) != 0:
                hits = pep[0].getHits()

                if len(hits) != 0:
                    besthit = hits[0]  # type: PeptideHit
                    yield f.getUniqueId(), besthit.getSequence().toString(), f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals
                else:
                    yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals
            else:
                yield f.getUniqueId(), None, f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals

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
    def __init__(self):
        super().__init__()

    ## TODO add chromatogram version
    ## TODO metadata df?

    def get_df(self, melt : bool = False):
        if melt:
            spectraarr = np.fromiter(((spec.getRT(), point[0], point[1]) for spec in self for point in zip(*spec.get_peaks())), dtype=[('RT', 'f'), ('mz', 'f'), ('inty', 'f')])

            return pd.DataFrame(data=spectraarr)

        else:
            cols = ["RT", "mzarray", "intarray"]

            return pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in self), columns=cols)

MSExperiment = MSExperimentDF

class DFConverter:
    '''Contains functions for the conversion of pyOpenMS data structures to a pandas DataFrame.'''
    
    def peptide_identifications_to_df(self, peps: List[PeptideIdentification], decode_ontology : bool = True):
        '''Converts a list of peptide identifications to a pandas DataFrame.

        Parameters:
        peps (List[PeptideIdentification]): list of PeptideIdentification objects
        decode_ontology (bool): decode meta value names
        
        Returns:
        pandas.DataFrame: peptide identifications in a DataFrame
        '''
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

        # TODO get score type name
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
                return tuple([pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), np.NA, np.NA] + [np.NA]*len(metavals))
            else:
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
                        ret.append(np.NA)
                return tuple(ret)

        #TODO implement hasHits function in C++
        psmarr = np.fromiter((extract(pep) for pep in peps), dtype=dt, count=len(peps))
        #TODO make spectrum_ref the index, if available?
        return pd.DataFrame(psmarr)
