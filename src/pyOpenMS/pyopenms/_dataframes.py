from collections import defaultdict as _defaultdict
from typing import List, Union

from . import ConsensusMap as _ConsensusMap
from . import ConsensusFeature as _ConsensusFeature
from . import FeatureMap as _FeatureMap
from . import Feature as _Feature
from . import MRMFeature as _MRMFeature
from . import MSExperiment as _MSExperiment
from . import PeakMap as _PeakMap
from . import PeptideIdentification as _PeptideIdentification
from . import ControlledVocabulary as _ControlledVocabulary
from . import File as _File
from . import IonSource as _IonSource
from . import MSSpectrum as _MSSpectrum
from . import MSChromatogram as _MSChromatogram
from . import MRMTransitionGroupCP as _MRMTransitionGroupCP

import pandas as _pd
import numpy as _np
from enum import Enum as _Enum

class _ConsensusMapDF(_ConsensusMap):
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

            def extract_row_blocks_channel_wide_file_long(f: _ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                filerows = _defaultdict(lambda: [0] * len(labels))
                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row = filerows[header.filename]
                    row[label_to_idx[header.label]] = fh.getIntensity()
                return (f.getUniqueId(), filerows)

            def extract_rows_channel_wide_file_long(f: _ConsensusFeature):
                uniqueid, rowdict = extract_row_blocks_channel_wide_file_long(f)
                for file, row in rowdict.items():
                    row.append(file)
                    yield tuple([uniqueid] + row)

            if len(labels) == 1:
                labels[0] = "intensity"

            dtypes = [('id', _np.dtype('uint64'))] + list(zip(labels, ['f'] * len(labels)))
            dtypes.append(('file', 'U300'))

            intyarr = _np.fromiter(iter=gen(self, extract_rows_channel_wide_file_long), dtype=dtypes, count=self.size())

            return _pd.DataFrame(intyarr).set_index('id')

        else:
            # Specialized for LabelFree which has to have only one channel
            def extract_row_blocks_channel_long_file_wide_LF(f: _ConsensusFeature):
                subfeatures = f.getFeatureList()  # type: list[FeatureHandle]
                row = [0.] * len(files)

                for fh in subfeatures:
                    header = filemeta[fh.getMapIndex()]
                    row[file_to_idx[header.filename]] = fh.getIntensity()

                yield tuple([f.getUniqueId()] + row)

            dtypes = [('id', _np.dtype('uint64'))] + list(zip(files, ['f'] * len(files)))

            intyarr = _np.fromiter(iter=gen(self, extract_row_blocks_channel_long_file_wide_LF), dtype=dtypes, count=self.size())

            return _pd.DataFrame(intyarr).set_index('id')

    def get_metadata_df(self):
        """Generates a pandas DataFrame with feature meta data (sequence, charge, mz, RT, quality).

        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        Returns:
        pandas.DataFrame: DataFrame with metadata for each feature (such as: best identified sequence, charge, centroid RT/mz, fitting quality)
        """

        def gen(cmap: _ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        def extract_meta_data(f: _ConsensusFeature):
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

        mddtypes = [('id', _np.dtype('uint64')), ('sequence', 'U200'), ('charge', 'i4'),
                    ('RT', _np.dtype('double')), ('mz', _np.dtype('double')), ('quality', 'f')]

        mdarr = _np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return _pd.DataFrame(mdarr).set_index('id')
    
    def get_df(self):
        """Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        Returns:
        pandas.DataFrame: meta data and intensity DataFrame
        """
        return _pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)

# fix class module and name to show up correctly in readthedocs page generated with sphinx autodoc
# needs to link back to rst page of original class, which is pyopenms.ConsensusMap, NOT pyopenms._dataframes._ConsensusMapDF (wh)
ConsensusMap = _ConsensusMapDF
ConsensusMap.__module__ = _ConsensusMap.__module__
ConsensusMap.__name__ = 'ConsensusMap'

# TODO tell the advanced user that they could change this, in case they have different needs.
# TODO check if type could be inferred in the first pass
# TODO check if max. string lengths could be inferred in the first pass and how this affects runtime
# TODO check how runtime is affected if we use _np.append instead of _np.fromiter and use _np.dyte = object for strings
common_meta_value_types = {
    b'label': 'U50',
    b'spectrum_index': 'i',
    b'score_fit': 'f',
    b'score_correlation': 'f',
    b'FWHM': 'f',
    b'spectrum_native_id': 'U100',
    b'max_height': 'f',
    b'num_of_masstraces': 'i',
    b'masstrace_intensity': 'f', # TODO this is actually a DoubleList. Think about what to do here. For _np.fromiter we would need to set the length of the array.
    b'Group': 'U50',
    b'is_ungrouped_monoisotopic': 'i', # TODO this sounds very boolean to me
    b'left_width': 'f',
    b'right_width': 'f',
    b'total_xic': 'f',
    b'PeptideRef': 'U100',
    b'peak_apices_sum': 'f'
}
"""Global dict to define which autoconversion to numpy types is tried for certain metavalues.

This can be changed to your liking but only affects future exports of any OpenMS datastructure to dataframes.
Especially string lengths (i.e., U types) benefit from adaption to save memory. The default type is currently
hardcoded to U50 (i.e., 50 unicode characters)
"""

class _FeatureMapDF(_FeatureMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __get_prot_id_filename_from_pep_id(self, pep_id: _PeptideIdentification) -> str:
        """Gets the primary MS run path of the ProteinIdentification linked with the given PeptideIdentification.

        Parameters:
        pep_id: PeptideIdentification

        Returns:
        str: primary MS run path (filename) of the ProteinIdentification with the same identifier as the given PeptideIdentification
        """
        for prot in self.getProteinIdentifications():
            if prot.getIdentifier() == pep_id.getIdentifier():
                filenames = []
                prot.getPrimaryMSRunPath(filenames)
                if filenames and filenames[0] != '':
                    return filenames[0]
        return 'unknown'
    
    # meta_values = None (default), 'all' or list of meta value names
    def get_df(self, meta_values: Union[None, List[str], str] = None, export_peptide_identifications: bool = True):
        """Generates a pandas DataFrame with information contained in the FeatureMap.

        Optionally the feature meta values and information for the assigned PeptideHit can be exported.

        Parameters:
        meta_values: meta values to include (None, [custom list of meta value names] or 'all')

        export_peptide_identifications (bool): export sequence and score for best PeptideHit assigned to a feature.
        Additionally the ID_filename (file name of the corresponding ProteinIdentification) and the ID_native_id 
        (spectrum ID of the corresponding Feature) are exported. They are also annotated as meta values when 
        collecting all assigned PeptideIdentifications from a FeatureMap with FeatureMap.get_assigned_peptide_identifications().
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])
        
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
                yield from fun(f)

        def extract_meta_data(f: _Feature):
            """Extracts feature meta data.
            
            Extracts information from a given feature with the requested meta values and, if requested,
            the sequence, score and ID_filename (primary MS run path of the linked ProteinIdentification)
            of the best PeptideHit (first) assigned to that feature.

            Parameters:
            f (Feature): feature from which to extract the meta data

            Yields:
            tuple: tuple containing feature information, peptide information (optional) and meta values (optional)
            """
            pep = f.getPeptideIdentifications()  # type: list[PeptideIdentification]
            bb = f.getConvexHull().getBoundingBox2D()
                
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else _np.nan for m in meta_values]
            
            if export_peptide_identifications:
                if len(pep) > 0:
                    ID_filename = self.__get_prot_id_filename_from_pep_id(pep[0])
                    hits = pep[0].getHits()
                    if len(hits) > 0:
                        besthit = hits[0]
                        pep_values = (besthit.getSequence().toString(), besthit.getScore(), ID_filename, f.getMetaValue('spectrum_native_id'))
                    else:
                        pep_values = (None, None, ID_filename, f.getMetaValue('spectrum_native_id'))
                else:
                    pep_values = (None, None, None, None)
            else:
                pep_values = ()

            yield tuple([f.getUniqueId()]) + pep_values + (f.getCharge(), f.getRT(), f.getMZ(), bb[0][0], bb[1][0], bb[0][1], bb[1][1], f.getOverallQuality(), f.getIntensity(), *vals)

        cnt = self.size()

        mddtypes = [('feature_id', 'U100')]
        if export_peptide_identifications:
            mddtypes += [('peptide_sequence', 'U200'), ('peptide_score', 'f'), ('ID_filename', 'U100'), ('ID_native_id', 'U100')]
        mddtypes += [('charge', 'i4'), ('RT', _np.dtype('double')), ('mz', _np.dtype('double')), ('RTstart', _np.dtype('double')), ('RTend', _np.dtype('double')),
                    ('MZstart', _np.dtype('double')), ('MZend', _np.dtype('double')), ('quality', 'f'), ('intensity', 'f')]
        
        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = _np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return _pd.DataFrame(mdarr).set_index('feature_id')

    def get_assigned_peptide_identifications(self):
        """Generates a list with peptide identifications assigned to a feature.

        Adds 'ID_native_id' (feature spectrum id), 'ID_filename' (primary MS run path of corresponding ProteinIdentification)
        and 'feature_id' (unique ID of corresponding Feature) as meta values to the peptide hits.
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = _pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])

        Returns:
        [PeptideIdentification]: list of PeptideIdentification objects
        """
        result = []
        for f in self:
            for pep in f.getPeptideIdentifications():
                hits = []
                for hit in pep.getHits():
                    hit.setMetaValue('feature_id', str(f.getUniqueId()))
                    hit.setMetaValue('ID_filename', self.__get_prot_id_filename_from_pep_id(pep))
                    if f.metaValueExists('spectrum_native_id'):
                        hit.setMetaValue('ID_native_id', f.getMetaValue('spectrum_native_id'))
                    else:
                        hit.setMetaValue('ID_native_id', 'unknown')
                    hits.append(hit)
                pep.setHits(hits)
                result.append(pep)
        return result

FeatureMap = _FeatureMapDF
FeatureMap.__module__ = _FeatureMap.__module__
FeatureMap.__name__ = 'FeatureMap'


class _MSExperimentDF(_MSExperiment):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, long : bool = False):
        """Generates a pandas DataFrame with all peaks in the MSExperiment

        Parameters:
        long: set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
            replicated RT information. If False, returns rows in the style: rt, _np.array(mz), _np.array(int)
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        if long:
            cols = ["RT", "mz", "inty"]
            self.updateRanges()
            spectraarrs2d = self.get2DPeakDataLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ())
            return _pd.DataFrame(dict(zip(cols, spectraarrs2d)))

        cols = ["RT", "mzarray", "intarray"]

        return _pd.DataFrame(data=((spec.getRT(), *spec.get_peaks()) for spec in self), columns=cols)

    def get_ion_df(self):
        """Generates a pandas DataFrame with all peaks and the ionic mobility in the MSExperiment
        
        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """
        
        cols = ["RT", "mz", "inty", "IM"]
        self.updateRanges()
        spectraarrs2d = self.get2DPeakDataLongIon(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ())
        return _pd.DataFrame(dict(zip(cols, spectraarrs2d)))

    def get_massql_df(self, ion_mobility=False):
        """Exports data from MSExperiment to pandas DataFrames to be used with MassQL.

        The Python module massql allows queries in mass spectrometry data (MS1 and MS2
        data frames) in a SQL like fashion (https://github.com/mwang87/MassQueryLanguage).
        
        Both dataframes contain the columns:
        'i': intensity of a peak
        'i_norm': intensity normalized by the maximun intensity in the spectrum
        'i_tic_norm': intensity normalized by the sum of intensities (TIC) in the spectrum
        'mz': mass to charge of a peak
        'scan': number of the spectrum
        'rt': retention time of the spectrum
        'polarity': ion mode of the spectrum as integer value (positive: 1, negative: 2)
        'ion': the ionic mobility of a peak if ion parameter is True
        
        The MS2 dataframe contains additional columns:
        'precmz': mass to charge of the precursor ion
        'ms1scan': number of the corresponding MS1 spectrum
        'charge': charge of the precursor ion
        
        Parameters:
            ion (bool): if True, returns the ion mobility of the peaks.

        Returns:
        ms1_df (pandas.DataFrame): peak data of MS1 spectra
        ms2_df (pandas.DataFrame): peak data of MS2 spectra with precursor information
        """
        self.updateRanges()

        def _get_polarity(spec):
            '''Returns polarity as an integer value for the massql dataframe.
            
            According to massql positive polarity is represented by 1 and negative by 2.

            Parameters:
            spec (MSSpectrum): the spectrum to extract polarity

            Returns:
            int: polarity as int value according to massql specification
            '''
            polarity = spec.getInstrumentSettings().getPolarity()
            if polarity == _IonSource.Polarity.POLNULL:
                return 0
            elif polarity == _IonSource.Polarity.POSITIVE:
                return 1
            elif polarity == _IonSource.Polarity.NEGATIVE:
                return 2

        def _get_spec_arrays(mslevel):
            '''Get spectrum data as a matrix.

            Generator yields peak data from each spectrum (with specified MS level) as a numpy.ndarray.
            Normalized intensity values are calculated and the placeholder values replaced. For 'i_norm' and
            'i_tic_norm' the intensity values are divided by the maximum intensity value in the spectrum and 
            the sum of intensity values, respectively.

            Parameters:
            mslevel (int): only spectra with the given MS level will be considered

            Yields:
            _np.ndarray: 2D array with peak data (rows) from each spectrum
            '''
            for scan_num, spec in enumerate(self):
                if spec.getMSLevel() == mslevel:
                    mz, inty = spec.get_peaks()
                    # data for both DataFrames: i, i_norm, i_tic_norm, mz, scan, rt, polarity
                    data = (inty, inty/_np.amax(inty, initial=0), inty/_np.sum(inty), mz, scan_num + 1, spec.getRT()/60, _get_polarity(spec))
                    cols = 7
                    if mslevel == 2:
                        cols = 10
                        # data for MS2 only: precmz, ms1scan, charge
                        # set fallback values if no precursor is annotated (-1)
                        if spec.getPrecursors():
                            data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num)+1, spec.getPrecursors()[0].getCharge())
                        else:
                            data += (-1, -1, -1)
                    # create empty ndarr with shape according to MS level
                    ndarr = _np.empty(shape=(spec.size(), cols))
                    # set column values
                    for i in range(cols):
                        ndarr[:,i] = data[i]
                    yield ndarr

        def _get_ion_spec_arrays(mslevel):
            '''Get spectrum data as a matrix.

            Generator yields peak data from each spectrum (with specified MS level) as a numpy.ndarray.
            Normalized intensity values are calculated and the placeholder values replaced. For 'i_norm' and
            'i_tic_norm' the intensity values are divided by the maximum intensity value in the spectrum and 
            the sum of intensity values, respectively.

            Parameters:
            mslevel (int): only spectra with the given MS level will be considered

            Yields:
            _np.ndarray: 2D array with peak data (rows) from each spectrum
            '''
            for scan_num, spec in enumerate(self):
                if spec.getMSLevel() == mslevel:
                    mz, inty = spec.get_peaks() 
                    ion_array_idx, ion_unit = spec.getIMData()
                    ion_data_arr = spec.getFloatDataArrays()[ion_array_idx]
                    ion_data = ion_data_arr.get_data()

                    # data for both DataFrames: i, i_norm, i_tic_norm, mz, scan, rt, polarity
                    data = (inty, inty/_np.amax(inty, initial=0), inty/_np.sum(inty), mz, scan_num + 1, spec.getRT()/60, _get_polarity(spec), ion_data)
                    cols = 8
                    if mslevel == 2:
                        cols = 11
                        # data for MS2 only: precmz, ms1scan, charge
                        # set fallback values if no precursor is annotated (-1)
                        if spec.getPrecursors():
                            data += (spec.getPrecursors()[0].getMZ(), self.getPrecursorSpectrum(scan_num)+1, spec.getPrecursors()[0].getCharge())
                        else:
                            data += (-1, -1, -1)
                    # create empty ndarr with shape according to MS level
                    ndarr = _np.empty(shape=(spec.size(), cols))
                    # set column values
                    for i in range(cols):
                        ndarr[:,i] = data[i]
                    yield ndarr

        # create DataFrame for MS1 and MS2 with according column names and data types
        # if there are no spectra of given MS level return an empty DataFrame
        dtypes = {'i': 'float32', 'i_norm': 'float32', 'i_tic_norm': 'float32', 'mz': 'float64', 'scan': 'int32', 'rt': 'float32', 'polarity': 'int32'}
        if ion_mobility:
            dtypes = dict(dtypes, **{"IM": "float32"})
    
        if 1 in self.getMSLevels():
            spec_arrays = _get_spec_arrays(1) if not ion_mobility else _get_ion_spec_arrays(1)
            ms1_df = _pd.DataFrame(_np.concatenate(list(spec_arrays), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms1_df = _pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        dtypes = dict(dtypes, **{'precmz': 'float64', 'ms1scan': 'int32', 'charge': 'int32'})
        if 2 in self.getMSLevels():
            spec_arrays = _get_spec_arrays(2) if not ion_mobility else _get_ion_spec_arrays(2)
            ms2_df = _pd.DataFrame(_np.concatenate(list(spec_arrays), axis=0), columns=dtypes.keys()).astype(dtypes)
        else:
            ms2_df = _pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        return ms1_df, ms2_df
    
PeakMap = _MSExperimentDF
PeakMap.__module__ = _PeakMap.__module__
PeakMap.__name__ = 'PeakMap'

MSExperiment = _MSExperimentDF
MSExperiment.__module__ = _MSExperiment.__module__
MSExperiment.__name__ = 'MSExperiment'


# TODO think about the best way for such top-level function. IMHO in python, encapsulation in a stateless class in unnecessary.
#   We should probably not just import this whole submodule without prefix.
def peptide_identifications_to_df(peps: List[_PeptideIdentification], decode_ontology : bool = True,
                                  default_missing_values: dict = {bool: False, int: -9999, float: _np.nan, str: ''},
                                  export_unidentified : bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.
    Parameters:
    peps (List[PeptideIdentification]): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names
    default_missing_values: default value for missing values for each data type
    export_unidentified: export PeptideIdentifications without PeptideHit
    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
    """
    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}

    # filter out PeptideIdentifications without PeptideHits if export_unidentified == False
    count = len(peps)
    if not export_unidentified:
        count = sum(len(pep.getHits()) > 0 for pep in peps)

    # get all possible metavalues
    metavals = []
    types = []
    mainscorename = "score"
    for pep in peps:
        hits = pep.getHits()
        if not len(hits) == 0:
            mvs = []
            hits[0].getKeys(mvs)
            metavals += mvs
            mainscorename = pep.getScoreType()

    metavals = list(set(metavals))

    # get type of all metavalues
    for k in metavals:
        if k == b"target_decoy":
            types.append('?')
        else:
            for p in peps:
                hits = p.getHits()
                if not len(hits) == 0:
                    mv = hits[0].getMetaValue(k)
                    types.append(switchDict[type(mv)])
                    break

    # get default value for each type in types to append if there are no hits in a PeptideIdentification
    def get_key(val):
        for key, value in switchDict.items():
            if val == value:
                return key
    dmv = [default_missing_values[get_key(t)] for t in types]

    decodedMVs = [m.decode("utf-8") for m in metavals]
    if decode_ontology:
        cv = _ControlledVocabulary()
        cv.loadFromOBO("psims", _File.getOpenMSDataPath() + "/CV/psi-ms.obo")
        clearMVs = [cv.getTerm(m).name if m.startswith("MS:") else m for m in decodedMVs]
    else:
        clearMVs = decodedMVs
        
    clearcols = ["id", "RT", "mz", mainscorename, "charge", "protein_accession", "start", "end", "P_ID", "PSM_ID"] + clearMVs
    coltypes = ['U100', 'f', 'f', 'f', 'i','U1000', 'U1000', 'U1000', 'i', 'i'] + types
    dt = list(zip(clearcols, coltypes))

    def extract(pep, pep_idx):
        hits = pep.getHits()
        if not hits:
            if export_unidentified:
                return (pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), default_missing_values[float], default_missing_values[int],
                        default_missing_values[str], default_missing_values[str], default_missing_values[str], pep_idx, default_missing_values[int], *dmv)
            else:
                return

        besthit = hits[0]
        ret = [pep.getIdentifier().encode('utf-8'), pep.getRT(), pep.getMZ(), besthit.getScore(), besthit.getCharge()] 
        # add accession, start and end positions of peptide evidences as comma separated str (like in mzTab)
        evs = besthit.getPeptideEvidences()
        ret += [','.join(v) if v else default_missing_values[str] for v in ([e.getProteinAccession() for e in evs],
                                                                            [str(e.getStart()) for e in evs],
                                                                            [str(e.getEnd()) for e in evs])]

        ret += [str(pep_idx), 0] # we currently only export the first hit

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
                ret.append(default_missing_values[type(val)])
        return tuple(ret)

    return _pd.DataFrame(_np.fromiter((extract(pep, pep_idx) for pep_idx, pep in enumerate(peps)), dtype=dt, count=count))


def update_scores_from_df(peps: List[_PeptideIdentification], df : _pd.DataFrame, main_score_name : str):
    """
    Updates the scores in PeptideIdentification objects using a pandas dataframe.
                
    :param peps: list of PeptideIdentification objects
    :param df: pandas dataframe obtained by converting peps to a dataframe. Minimum required: P_ID column and column with name passed by main_score_name
    :return: the updated list of peptide identifications
    """

    rets = peps

    for index, row in df.iterrows():
        pid_index = int(row["P_ID"])
        pi = _PeptideIdentification(peps[pid_index])
        pi.setScoreType(main_score_name)
        hits = pi.getHits() # type: list[PeptideHit]
        if len(hits) > 0:
            best_hit = hits[0]
            best_hit.setScore(float(row[main_score_name]))
            hits[0] = best_hit
            pi.setHits(hits)

        rets[pid_index] = pi

    return rets

def _add_meta_values(df: _pd.DataFrame, object: any) -> _pd.DataFrame:
    """
    Adds metavalues from given object to given DataFrame.
    
    Args:
        df (pd.DataFrame): DataFrame to which metavalues will be added.
        object (any): Object from which metavalues will be extracted.
    
    Returns:
        pd.DataFrame: DataFrame with added meta values.
    """
    mvs = []
    object.getKeys(mvs)
    for k in mvs:
        v = object.getMetaValue(k)
        dtype = 'U100'
        try:
            v = int(v)
            dtype = int
        except ValueError:
            try:
                v = float(v)
                dtype = 'double'
            except ValueError:
                dtype = f'U{len(v)}'
        
        df[k.decode()] = _np.full(df.shape[0], v, dtype=_np.dtype(dtype))

    return df

class _MSSpectrumDF(_MSSpectrum):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, export_meta_values: bool = True) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the MSSpectrum.

        Args:
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pd.DataFrame: DataFrame representation of the MSSpectrum.
        """
        mzs, intensities = self.get_peaks()

        df = _pd.DataFrame({'mz': mzs, 'intensity': intensities})

        cnt = df.shape[0]
        
        # ion mobility
        df['ion_mobility'] = _np.array([i for i in self.getFloatDataArrays()[0]]) if self.containsIMData() else _np.nan
        df['ion_mobility_unit'] = _np.full(cnt, self.getDriftTimeUnitAsString(), dtype=_np.dtype(f'U{len(self.getDriftTimeUnitAsString())}'))

        df['ms_level'] = _np.full(cnt, self.getMSLevel(), dtype=_np.dtype('uint16'))

        precs = self.getPrecursors()
        df['precursor_mz'] = _np.full(cnt, (precs[0].getMZ() if precs else 0.0), dtype=_np.dtype('double'))
        df['precursor_charge'] = _np.full(cnt, (precs[0].getCharge() if precs else 0), dtype=_np.dtype('uint16'))
        
        df['native_id'] = _np.full(cnt, self.getNativeID(), dtype=_np.dtype('U100'))

        # peptide sequence
        peps = self.getPeptideIdentifications()  # type: list[PeptideIdentification]
        seq = ''
        if peps:
            hits = peps[0].getHits()
            if hits:
                seq = hits[0].getSequence().toString()
        df['sequence'] = _np.full(cnt, seq, dtype=_np.dtype(f'U{len(seq)}'))

        # ion annotations in string data array with names IonName or IonNames
        ion_annotations = _np.full(cnt, '', dtype=_np.dtype('U1'))
        for sda in self.getStringDataArrays():
            if sda.getName() == 'IonNames':
                decoded = [ion.decode() for ion in sda]
                if len(decoded) == df.shape[0]:
                    ion_annotations = _np.array(decoded, dtype=_np.dtype(f'U{len(max(decoded))}'))
                    break
        df['ion_annotation'] = ion_annotations

        if export_meta_values:
            df = _add_meta_values(df, self)

        return df

MSSpectrum = _MSSpectrumDF
MSSpectrum.__module__ = _MSSpectrum.__module__
MSSpectrum.__name__ = 'MSSpectrum'

class _ChromatogramType(_Enum):
    MASS_CHROMATOGRAM = 0
    TOTAL_ION_CURRENT_CHROMATOGRAM = 1
    SELECTED_ION_CURRENT_CHROMATOGRAM = 2
    BASEPEAK_CHROMATOGRAM = 3
    SELECTED_ION_MONITORING_CHROMATOGRAM = 4
    SELECTED_REACTION_MONITORING_CHROMATOGRAM = 5
    ELECTROMAGNETIC_RADIATION_CHROMATOGRAM = 6
    ABSORPTION_CHROMATOGRAM = 7
    EMISSION_CHROMATOGRAM = 8

class _MSChromatogramDF(_MSChromatogram):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, export_meta_values: bool = True) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the MSChromatogram.

        time: The retention time (in seconds) of the chromatographic peaks.
        intensity: The intensity (abundance) of the signal at each time point.
        chromatogram_type: The type of chromatogram.
        precursor_mz: The mass-to-charge of the precursor ion.
        precursor_charge: The charge of the precursor ion.
        comment: A comment assigned to the chromatogram.
        native_id: The chromatogram native identifier.

        Args:
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pd.DataFrame: DataFrame representation of the MSChromatogram.
        """
        def extract_data(c: _MSChromatogram):
            rts, intys = c.get_peaks()
            for rt, inty in zip(rts, intys):
                yield rt, inty

        cnt = len(self.get_peaks()[0])

        dtypes = [('time', _np.dtype('double')), ('intensity', _np.dtype('uint64'))]

        arr = _np.fromiter(iter=extract_data(self), dtype=dtypes, count=cnt)

        df = _pd.DataFrame(arr)

        df['chromatogram_type'] = _np.full(cnt, _ChromatogramType(self.getChromatogramType()).name, dtype=_np.dtype('U100'))

        df['precursor_mz'] = _np.full(cnt, self.getPrecursor().getMZ(), dtype=_np.dtype('double'))
        df['precursor_charge'] = _np.full(cnt, self.getPrecursor().getCharge(), dtype=_np.dtype('uint16'))

        df['product_mz'] = _np.full(cnt, self.getProduct().getMZ(), dtype=_np.dtype('double'))

        df['comment'] = _np.full(cnt, self.getComment(), dtype=_np.dtype('U100'))

        df['native_id'] = _np.full(cnt, self.getNativeID(), dtype=_np.dtype('U100'))

        if export_meta_values:
            df = _add_meta_values(df, self)

        return df

MSChromatogram = _MSChromatogramDF
MSChromatogram.__module__ = _MSChromatogram.__module__
MSChromatogram.__name__ = 'MSChromatogram'

class _MRMTransitionGroupCPDF(_MRMTransitionGroupCP):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_chromatogram_df(self, export_meta_values: bool = True) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the Chromatograms stored in MRMTransitionGroupCP.

        rt: The retention time of the transition group.
        intensity: The intensity of the transition group.
        precursor_mz: The mass-to-charge ratio of the precursor ion.
        precursor_charge: The charge of the precursor ion.
        product_mz: The mass-to-charge ratio of the product ion.
        product_charge: The charge of the product ion.
        native_id: The native identifier of the transition group.

        Args:
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pd.DataFrame: DataFrame representation of the chromatograms stored in MRMTransitionGroupCP.
        """
        chroms = self.getChromatograms()
        out = [ _MSChromatogramDF(c).get_df(export_meta_values=export_meta_values) for c in chroms ]
        return _pd.concat(out)
    
    def get_feature_df(self, meta_values: Union[None, List[str], str] = None) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the Features stored in MRMTransitionGroupCP.

        rt: The retention time of the transition group.
        intensity: The intensity of the transition group.
        precursor_mz: The mass-to-charge ratio of the precursor ion.
        precursor_charge: The charge of the precursor ion.
        product_mz: The mass-to-charge ratio of the product ion.
        product_charge: The charge of the product ion.
        native_id: The native identifier of the transition group.

        Args:
            export_meta_values (bool): Whether to export meta values.

        Returns:
            pd.DataFrame: DataFrame representation of the Features stored in MRMTransitionGroupCP.
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

        features = self.getFeatures()
        
        def gen(features: List[_MRMFeature], fun):
            for f in features:
                yield from fun(f)

        def extract_meta_data(f: _MRMFeature):
            """Extracts feature meta data.
            
            Extracts information from a given feature with the requested meta values and, if requested,
            the sequence, score and ID_filename (primary MS run path of the linked ProteinIdentification)
            of the best PeptideHit (first) assigned to that feature.

            Parameters:
            f (Feature): feature from which to extract the meta data

            Yields:
            tuple: tuple containing feature information, and meta values (optional)
            """
            vals = [f.getMetaValue(m) if f.metaValueExists(m) else _np.nan for m in meta_values]
            
            yield tuple((f.getUniqueId(), f.getRT(), f.getIntensity(), f.getOverallQuality(), *vals))

        features = self.getFeatures()


        mddtypes = [('feature_id', _np.dtype('uint64')), ('RT', 'f'), ('intensity', 'f'), ('quality', 'f')]

        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = _np.fromiter(iter=gen(features, extract_meta_data), dtype=mddtypes, count=len(features))

        return _pd.DataFrame(mdarr).set_index('feature_id')

# fix class module and name to show up correctly in readthedocs page generated with sphinx autodoc
# needs to link back to rst page of original class, which is pyopenms.MRMTransitionGroupCP, NOT pyopenms._dataframes._MRMTransitionGroupCPDF (wh)
MRMTransitionGroupCP = _MRMTransitionGroupCPDF
MRMTransitionGroupCP.__module__ = _MRMTransitionGroupCP.__module__
MRMTransitionGroupCP.__name__ = 'MRMTransitionGroupCP'