from collections import defaultdict as _defaultdict
from typing import List, Union

from . import ConsensusMap as _ConsensusMap
from . import ConsensusFeature as _ConsensusFeature
from . import FeatureMap as _FeatureMap
from . import Feature as _Feature
from . import MRMFeature as _MRMFeature
from . import MSExperiment as _MSExperiment
from . import PeakMap as _PeakMap
from . import PeptideIdentificationList as _PeptideIdentificationList
from . import PeptideIdentification as _PeptideIdentification
from . import ControlledVocabulary as _ControlledVocabulary
from . import File as _File
from . import IonSource as _IonSource
from . import MSSpectrum as _MSSpectrum
from . import PeakSpectrum as _PeakSpectrum
from . import MSChromatogram as _MSChromatogram
from . import MRMTransitionGroupCP as _MRMTransitionGroupCP
from . import Mobilogram as _Mobilogram
from . import DataValue as _DataValue

import pandas as _pd
import numpy as _np
from enum import Enum as _Enum


class _MSSpectrumDF(_MSSpectrum):
    """MSSpectrum with DataFrame export capabilities.

    This class extends MSSpectrum with a get_df() method that converts
    spectrum data to a pandas DataFrame.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, columns: Union[None, List[str]] = None, export_meta_values: bool = True) -> _pd.DataFrame:
        """Returns a pandas DataFrame representation of the MSSpectrum.

        This method converts the spectrum data (peaks, metadata, precursor info,
        ion mobility) into a pandas DataFrame format.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns. Use get_df_columns()
                                   to discover available columns.
            export_meta_values (bool): Whether to include meta values. Only applies
                                       when columns=None. Defaults to True.

        Returns:
            pd.DataFrame: DataFrame with requested columns. Default columns include:
                - mz: m/z values of peaks
                - intensity: intensity values of peaks
                - rt: retention time (replicated for each peak)
                - ms_level: MS level (replicated for each peak)
                - native_id: native spectrum identifier
                - ion_mobility: ion mobility values (if IM data present)
                - precursor_mz: precursor m/z (if precursor present)
                - precursor_charge: precursor charge (if precursor present)
                - ion_annotation: ion annotation strings (if IonNames present)
                - Additional meta value columns (if export_meta_values=True)

        Example:
            >>> # Get all default columns
            >>> df = spectrum.get_df()

            >>> # Discover available columns
            >>> print(spectrum.get_df_columns())

            >>> # Get only specific columns (faster)
            >>> df = spectrum.get_df(columns=['mz', 'intensity'])

            >>> # Get all columns including non-defaults like ion_mobility_unit
            >>> cols = spectrum.get_df_columns()
            >>> cols.append('ion_mobility_unit')
            >>> df = spectrum.get_df(columns=cols)
        """
        # Try calling get_data_dict with columns parameter
        # inspect.signature doesn't work for Cython builtins
        try:
            data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
        except TypeError:
            # Fallback for older Cython builds without column selection
            data_dict = self.get_data_dict(export_meta_values=export_meta_values)
        return _pd.DataFrame(data_dict)


# Fix class module and name to show up correctly in documentation
MSSpectrum = _MSSpectrumDF
MSSpectrum.__module__ = _MSSpectrum.__module__
MSSpectrum.__name__ = 'MSSpectrum'

PeakSpectrum = _MSSpectrumDF
PeakSpectrum.__module__ = _PeakSpectrum.__module__
PeakSpectrum.__name__ = 'PeakSpectrum'


class _ConsensusMapDF(_ConsensusMap):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df_columns(self, columns: str = 'default') -> List[str]:
        """Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        Args:
            columns (str): 'default' for standard columns, 'all' for all available columns.

        Returns:
            list: List of column name strings.

        Example:
            >>> cmap.get_df_columns()
            ['sequence', 'charge', 'rt', 'mz', 'quality', 'intensity_file1', ...]
        """
        # Metadata columns
        cols = ['sequence', 'charge', 'rt', 'mz', 'quality']

        # Intensity columns depend on experiment type
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()

        if labelfree:
            # File-wide columns for label-free
            files = list(set([header.filename for header in filemeta.values()]))
            cols.extend(files)
        else:
            # Label columns for labelled experiments
            labels = list(set([header.label for header in filemeta.values()]))
            if len(labels) == 1:
                labels[0] = "intensity"
            cols.extend(labels)
            cols.append('file')

        return cols

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
        """Generates a pandas DataFrame with feature meta data.

        Columns: sequence, charge, rt, mz, quality (indexed by 'id').

        Resulting DataFrame can be joined with result from get_intensity_df by their index 'id'.

        Returns:
            pandas.DataFrame: DataFrame with metadata for each feature (sequence, charge,
                             rt, mz, quality). All column names are lowercase snake_case.
        """

        def gen(cmap: _ConsensusMap, fun):
            for f in cmap:
                yield from fun(f)

        def extract_meta_data(f: _ConsensusFeature):
            pep = f.getPeptideIdentifications()  # type: _PeptideIdentificationList

            if pep.size() != 0:
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
                    ('rt', _np.dtype('double')), ('mz', _np.dtype('double')), ('quality', 'f')]

        mdarr = _np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        return _pd.DataFrame(mdarr).set_index('id')

    def get_df(self, columns: Union[None, List[str]] = None):
        """Generates a pandas DataFrame with both consensus feature meta data and intensities from each sample.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_df_columns()
                                   to discover available columns.

        Returns:
            pandas.DataFrame: meta data and intensity DataFrame

        Example:
            >>> # Get all columns
            >>> df = cmap.get_df()

            >>> # Discover available columns
            >>> print(cmap.get_df_columns())

            >>> # Get only specific columns
            >>> df = cmap.get_df(columns=['sequence', 'mz', 'intensity'])
        """
        if columns is None:
            # No column selection - get everything
            df = _pd.concat([self.get_metadata_df(), self.get_intensity_df()], axis=1)
            return df

        # Efficient column selection: only compute what's needed
        requested = set(columns)
        metadata_cols = {'sequence', 'charge', 'rt', 'mz', 'quality'}

        # Get intensity column names
        labelfree = self.getExperimentType() == "label-free"
        filemeta = self.getColumnHeaders()
        if labelfree:
            intensity_cols = set([header.filename for header in filemeta.values()])
        else:
            labels = list(set([header.label for header in filemeta.values()]))
            if len(labels) == 1:
                labels[0] = "intensity"
            intensity_cols = set(labels)
            intensity_cols.add('file')

        need_metadata = bool(requested & metadata_cols)
        need_intensity = bool(requested & intensity_cols)

        dfs = []
        if need_metadata:
            dfs.append(self.get_metadata_df())
        if need_intensity:
            dfs.append(self.get_intensity_df())

        if not dfs:
            # No columns match - return empty DataFrame with index
            return _pd.DataFrame(index=_pd.Index([], name='id'))

        if len(dfs) == 1:
            df = dfs[0]
        else:
            df = _pd.concat(dfs, axis=1)

        # Filter to requested columns
        available_cols = [c for c in columns if c in df.columns]
        return df[available_cols]

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
    b'leftWidth': 'f',
    b'rightWidth': 'f',
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

    def get_df_columns(self, columns: str = 'default', export_peptide_identifications: bool = True) -> List[str]:
        """Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        Args:
            columns (str): 'default' for standard columns, 'all' to include all meta values.
            export_peptide_identifications (bool): Whether to include peptide ID columns.

        Returns:
            list: List of column name strings.

        Example:
            >>> fmap.get_df_columns()
            ['feature_id', 'peptide_sequence', 'charge', 'rt', 'mz', ...]
        """
        cols = ['feature_id']

        if export_peptide_identifications:
            cols.extend(['peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'])

        cols.extend(['charge', 'rt', 'mz', 'rt_start', 'rt_end', 'mz_start', 'mz_end', 'quality', 'intensity'])

        # Add meta values if 'all' requested
        if columns == 'all':
            meta_values = set()
            for f in self:
                mvs = []
                f.getKeys(mvs)
                for m in mvs:
                    meta_values.add(m.decode() if isinstance(m, bytes) else m)
            cols.extend(sorted(meta_values))

        return cols

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
    def get_df(self, columns: Union[None, List[str]] = None, meta_values: Union[None, List[str], str] = None, export_peptide_identifications: bool = True):
        """Generates a pandas DataFrame with information contained in the FeatureMap.

        Optionally the feature meta values and information for the assigned PeptideHit can be exported.

        Parameters:
        columns (list or None): List of column names to include. If None,
                               includes all columns. Use get_df_columns() to discover available columns.

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

        Example:
            >>> # Get all columns
            >>> df = fmap.get_df()

            >>> # Discover available columns
            >>> print(fmap.get_df_columns())

            >>> # Get only specific columns
            >>> df = fmap.get_df(columns=['feature_id', 'mz', 'rt', 'intensity'])
        """
        # Determine if peptide IDs are actually needed based on column selection
        pep_id_cols = {'peptide_sequence', 'peptide_score', 'ID_filename', 'ID_native_id'}
        if columns is not None:
            requested = set(columns)
            # Only export peptide IDs if explicitly requested or export_peptide_identifications is True
            # and at least one peptide column is requested
            need_pep_ids = export_peptide_identifications and bool(requested & pep_id_cols)
        else:
            need_pep_ids = export_peptide_identifications

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
            pep = f.getPeptideIdentifications()  # type: _PeptideIdentificationList
            bb = f.getConvexHull().getBoundingBox2D()

            vals = [f.getMetaValue(m) if f.metaValueExists(m) else _np.nan for m in meta_values]

            if need_pep_ids:
                if pep.size() > 0:
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
        if need_pep_ids:
            mddtypes += [('peptide_sequence', 'U200'), ('peptide_score', 'f'), ('ID_filename', 'U100'), ('ID_native_id', 'U100')]
        mddtypes += [('charge', 'i4'), ('rt', _np.dtype('double')), ('mz', _np.dtype('double')), ('rt_start', _np.dtype('double')), ('rt_end', _np.dtype('double')),
                    ('mz_start', _np.dtype('double')), ('mz_end', _np.dtype('double')), ('quality', 'f'), ('intensity', 'f')]

        for meta_value in meta_values:
            if meta_value in common_meta_value_types:
                mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
            else:
                mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = _np.fromiter(iter=gen(self, extract_meta_data), dtype=mddtypes, count=cnt)

        df = _pd.DataFrame(mdarr).set_index('feature_id')

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns or c == 'feature_id']
            # feature_id is the index, handle it specially
            if 'feature_id' not in available_cols:
                available_cols = [c for c in columns if c in df.columns]
            df = df[[c for c in available_cols if c in df.columns]]

        return df

    def get_assigned_peptide_identifications(self):
        """Generates a list with peptide identifications assigned to a feature.

        Adds 'ID_native_id' (feature spectrum id), 'ID_filename' (primary MS run path of corresponding ProteinIdentification)
        and 'feature_id' (unique ID of corresponding Feature) as meta values to the peptide hits.
        A DataFrame from the assigned peptides generated with peptide_identifications_to_df(assigned_peptides) can be
        merged with the FeatureMap DataFrame with:
        merged_df = _pd.merge(feature_df, assigned_peptide_df, on=['feature_id', 'ID_native_id', 'ID_filename'])

        Returns:
        _PeptideIdentificationList: list of PeptideIdentification objects
        """
        result = _PeptideIdentificationList()
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
                result.push_back(pep)
        return result

FeatureMap = _FeatureMapDF
FeatureMap.__module__ = _FeatureMap.__module__
FeatureMap.__name__ = 'FeatureMap'


class _MSExperimentDF(_MSExperiment):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df_columns(self, long: bool = False) -> List[str]:
        """Returns a list of column names that get_df() would produce.

        Useful for discovering available columns before export.

        Args:
            long (bool): If True, returns columns for long format.
                        If False, returns columns for compact format.

        Returns:
            list: List of column name strings.

        Example:
            >>> exp.get_df_columns(long=True)
            ['rt', 'mz', 'intensity', 'ms_level']

            >>> exp.get_df_columns(long=False)
            ['rt', 'ms_level', 'mz_array', 'intensity_array']
        """
        if long:
            return ['rt', 'mz', 'intensity', 'ms_level']
        else:
            return ['rt', 'ms_level', 'mz_array', 'intensity_array']

    def get_df(self, columns: Union[None, List[str]] = None, ms_levels: List[int] = [], long : bool = False):
        """Generates a pandas DataFrame with all peaks in the MSExperiment

        Parameters:
        columns (list or None): List of column names to include. If None,
                               includes all columns. Use get_df_columns() to discover available columns.
        ms_levels (List[int]): Get only spectra with the given MS levels. Default is an empty list, which means all MS levels will be included.
        long (bool): set to True if you want to have a long/expanded/melted dataframe with one row per peak. Faster but
            replicated RT information. If False, returns rows in the style: rt, _np.array(mz), _np.array(int)

        Returns:
        pandas.DataFrame: feature information stored in a DataFrame

        Example:
            >>> # Get all columns
            >>> df = exp.get_df()

            >>> # Discover available columns
            >>> print(exp.get_df_columns())

            >>> # Get only specific columns
            >>> df = exp.get_df(columns=['rt', 'mz', 'intensity'], long=True)
        """
        self.updateRanges()
        if not ms_levels:
            ms_levels = self.getMSLevels()
        if long:
            cols = ["rt", "mz", "intensity"]
            dfs = []
            for ms_level in ms_levels:
                spectraarrs2d = self.get2DPeakDataLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ(), ms_level)
                df = _pd.DataFrame(dict(zip(cols, spectraarrs2d)))
                df["ms_level"] = ms_level
                dfs.append(df)
            df = _pd.concat(dfs, ignore_index=True)
        else:
            cols = ["rt", "ms_level", "mz_array", "intensity_array"]
            df = _pd.DataFrame(data=((spec.getRT(), spec.getMSLevel(), *spec.get_peaks()) for spec in self if spec.getMSLevel() in ms_levels), columns=cols)

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns]
            df = df[available_cols]

        return df

    def get_ion_df(self):
        """Generates a pandas DataFrame with all peaks and the ion mobility in the MSExperiment

        Returns:
        pandas.DataFrame: feature information stored in a DataFrame
        """

        cols = ["rt", "mz", "intensity", "ion_mobility"]
        self.updateRanges()
        spectraarrs2d = self.get2DPeakDataIMLong(self.getMinRT(), self.getMaxRT(), self.getMinMZ(), self.getMaxMZ(), 1)
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
            dtypes = dict(dtypes, **{"ion_mobility": "float32"})
    
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
def peptide_identifications_to_df(peps: _PeptideIdentificationList, decode_ontology : bool = True,
                                  default_missing_values: dict = None,
                                  export_unidentified : bool = True):
    """Converts a list of peptide identifications to a pandas DataFrame.
    Parameters:
    peps (PeptideIdentificationList): list of PeptideIdentification objects
    decode_ontology (bool): decode meta value names
    default_missing_values: default value for missing values for each data type
    export_unidentified: export PeptideIdentifications without PeptideHit
    Returns:
    pandas.DataFrame: peptide identifications in a DataFrame
    """

    if default_missing_values is None:
        default_missing_values = {bool: False, int: -9999, float: _np.nan, str: ''}
        
    switchDict = {bool: '?', int: 'i', float: 'f', str: 'U100'}

    # filter out PeptideIdentifications without PeptideHits if export_unidentified == False
    count = peps.size()
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
        
    clearcols = ["id", "rt", "mz", mainscorename, "charge", "protein_accession", "start", "end", "P_ID", "PSM_ID"] + clearMVs
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


def update_scores_from_df(peps: _PeptideIdentificationList, df : _pd.DataFrame, main_score_name : str):
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
        dv = object.getMetaValue(k)
        col_name = k.decode()

        try:
            # Handle native Python types (returned by autowrap)
            if isinstance(dv, float):
                value = dv
                dtype = "float64"
            elif isinstance(dv, int):
                value = dv
                dtype = "int64"
            elif isinstance(dv, str):
                value = dv
                dtype = f"U{max(1, len(value))}"
            elif isinstance(dv, bytes):
                value = dv.decode()
                dtype = f"U{max(1, len(value))}"
            # Handle DataValue objects (if ever returned)
            elif hasattr(dv, 'valueType'):
                if dv.valueType() == _DataValue.STRING_VALUE:
                    value = dv.toString().decode()
                    dtype = f"U{max(1, len(value))}"
                elif dv.valueType() == _DataValue.INT_VALUE:
                    value = dv.toInt()
                    dtype = "int32"
                elif dv.valueType() == _DataValue.DOUBLE_VALUE:
                    value = dv.toDouble()
                    dtype = "float64"
                elif dv.valueType() == _DataValue.EMPTY_VALUE:
                    continue
                else:
                    value = str(dv)
                    dtype = "object"
            else:
                value = str(dv)
                dtype = "object"

            df[col_name] = _np.full(df.shape[0], value, dtype=dtype)

        except Exception:
            df[col_name] = _np.full(df.shape[0], str(dv), dtype='object')

    return df

class _MSChromatogramDF(_MSChromatogram):
    """MSChromatogram with DataFrame export capabilities.

    This class extends MSChromatogram with a get_df() method that converts
    chromatogram data to a pandas DataFrame.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, columns: Union[None, List[str]] = None, export_meta_values: bool = True) -> _pd.DataFrame:
        """Returns a pandas DataFrame representation of the MSChromatogram.

        This method converts the chromatogram data (peaks, metadata, precursor/product info)
        into a pandas DataFrame format.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns. Use get_df_columns()
                                   to discover available columns.
            export_meta_values (bool): Whether to include meta values. Only applies
                                       when columns=None. Defaults to True.

        Returns:
            pd.DataFrame: DataFrame with requested columns. Default columns include:
                - rt: retention time (in seconds)
                - intensity: signal intensity at each time point
                - precursor_mz: precursor m/z
                - precursor_charge: precursor charge
                - product_mz: product m/z
                - native_id: chromatogram native identifier
                - Additional meta value columns (if export_meta_values=True)

            Non-default columns (must be explicitly requested):
                - chromatogram_type: type of chromatogram
                - comment: chromatogram comment

        Example:
            >>> # Get all default columns
            >>> df = chrom.get_df()

            >>> # Discover available columns
            >>> print(chrom.get_df_columns())

            >>> # Get only specific columns (faster)
            >>> df = chrom.get_df(columns=['rt', 'intensity'])

            >>> # Get all columns including non-defaults
            >>> cols = chrom.get_df_columns('all')
            >>> df = chrom.get_df(columns=cols)
        """
        # Use get_data_dict from Cython addon (similar to MSSpectrum pattern)
        try:
            data_dict = self.get_data_dict(columns=columns, export_meta_values=export_meta_values)
        except TypeError:
            # Fallback for older Cython builds without column selection
            data_dict = self.get_data_dict(export_meta_values=export_meta_values)
        return _pd.DataFrame(data_dict)

MSChromatogram = _MSChromatogramDF
MSChromatogram.__module__ = _MSChromatogram.__module__
MSChromatogram.__name__ = 'MSChromatogram'


class _MobilogramDF(_Mobilogram):
    """Mobilogram with DataFrame export capabilities.

    This class extends Mobilogram with a get_df() method that converts
    mobilogram data to a pandas DataFrame.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_df(self, columns: Union[None, List[str]] = None) -> _pd.DataFrame:
        """Returns a pandas DataFrame representation of the Mobilogram.

        This method converts the mobilogram data (peaks, metadata)
        into a pandas DataFrame format.

        Note: Mobilogram does not support meta values (no MetaInfoInterface).

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns. Use get_df_columns()
                                   to discover available columns.

        Returns:
            pd.DataFrame: DataFrame with requested columns. Default columns include:
                - mobility: mobility values of peaks
                - intensity: intensity values of peaks
                - rt: retention time (replicated for each peak)
                - drift_time_unit: drift time unit string

        Example:
            >>> # Get all default columns
            >>> df = mobilogram.get_df()

            >>> # Discover available columns
            >>> print(mobilogram.get_df_columns())

            >>> # Get only specific columns (faster)
            >>> df = mobilogram.get_df(columns=['mobility', 'intensity'])
        """
        try:
            data_dict = self.get_data_dict(columns=columns)
        except TypeError:
            # Fallback for older Cython builds without column selection
            data_dict = self.get_data_dict()
        return _pd.DataFrame(data_dict)


Mobilogram = _MobilogramDF
Mobilogram.__module__ = _Mobilogram.__module__
Mobilogram.__name__ = 'Mobilogram'


class _MRMTransitionGroupCPDF(_MRMTransitionGroupCP):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_chromatogram_df_columns(self, columns: str = 'default', export_meta_values: bool = True) -> List[str]:
        """Returns a list of column names that get_chromatogram_df() would produce.

        Args:
            columns (str): 'default' for standard columns, 'all' for all available columns.
            export_meta_values (bool): Whether to include meta value columns.

        Returns:
            list: List of column name strings.
        """
        # Use the first chromatogram to get columns (all should have same structure)
        chroms = self.getChromatograms()
        if chroms:
            return _MSChromatogramDF(chroms[0]).get_df_columns(columns=columns, export_meta_values=export_meta_values)
        return ['rt', 'intensity', 'precursor_mz', 'precursor_charge', 'product_mz', 'native_id']

    def get_feature_df_columns(self, columns: str = 'default') -> List[str]:
        """Returns a list of column names that get_feature_df() would produce.

        Args:
            columns (str): 'default' for core columns, 'all' to include all meta values.

        Returns:
            list: List of column name strings.
        """
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

    def get_chromatogram_df(self, columns: Union[None, List[str]] = None, export_meta_values: bool = True) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the Chromatograms stored in MRMTransitionGroupCP.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all default columns.
            export_meta_values (bool): Whether to export meta values. Only applies
                                       when columns=None.

        Returns:
            pd.DataFrame: DataFrame representation of the chromatograms.

        Example:
            >>> # Get all default columns
            >>> df = mrm.get_chromatogram_df()

            >>> # Discover available columns
            >>> print(mrm.get_chromatogram_df_columns())

            >>> # Get only specific columns
            >>> df = mrm.get_chromatogram_df(columns=['rt', 'intensity'])
        """
        chroms = self.getChromatograms()
        out = [_MSChromatogramDF(c).get_df(columns=columns, export_meta_values=export_meta_values) for c in chroms]
        if out:
            return _pd.concat(out, ignore_index=True)
        return _pd.DataFrame()
    
    def get_feature_df(self, columns: Union[None, List[str]] = None, meta_values: Union[None, List[str], str] = None) -> _pd.DataFrame:
        """
        Returns a DataFrame representation of the Features stored in MRMTransitionGroupCP.

        Args:
            columns (list or None): List of column names to include. If None,
                                   includes all columns. Use get_feature_df_columns()
                                   to discover available columns.
            meta_values: meta values to include (None, [custom list of meta value names] or 'all')

        Returns:
            pd.DataFrame: DataFrame representation of the Features.

        Example:
            >>> # Get all columns
            >>> df = mrm.get_feature_df()

            >>> # Discover available columns
            >>> print(mrm.get_feature_df_columns())

            >>> # Get only specific columns
            >>> df = mrm.get_feature_df(columns=['feature_id', 'RT', 'intensity'])
        """
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

        # get all possible meta value keys in a set
        features = self.getFeatures()
        mddtypes = [('feature_id', _np.dtype('uint64')), ('rt', 'f'), ('intensity', 'f'), ('quality', 'f')]
        if meta_values is not None:
            # Add all possible meta values to meta_value array if 'all' is passed
            if meta_values == 'all':
                meta_values = set()
                for f in features:
                    mvs = []
                    f.getKeys(mvs)
                    for m in mvs:
                        meta_values.add(m)

            # Add meta_values to mddtypes
            for meta_value in meta_values:
                if meta_value in common_meta_value_types:
                    mddtypes.append((meta_value.decode(), common_meta_value_types[meta_value]))
                else:
                    mddtypes.append((meta_value.decode(), 'U50'))

        mdarr = _np.fromiter(iter=gen(features, extract_meta_data), dtype=mddtypes, count=len(features))

        df = _pd.DataFrame(mdarr).set_index('feature_id')

        # Filter columns if requested
        if columns is not None:
            available_cols = [c for c in columns if c in df.columns or c == 'feature_id']
            if 'feature_id' not in available_cols:
                available_cols = [c for c in columns if c in df.columns]
            df = df[[c for c in available_cols if c in df.columns]]

        return df

# fix class module and name to show up correctly in readthedocs page generated with sphinx autodoc
# needs to link back to rst page of original class, which is pyopenms.MRMTransitionGroupCP, NOT pyopenms._dataframes._MRMTransitionGroupCPDF (wh)
MRMTransitionGroupCP = _MRMTransitionGroupCPDF
MRMTransitionGroupCP.__module__ = _MRMTransitionGroupCP.__module__
MRMTransitionGroupCP.__name__ = 'MRMTransitionGroupCP'