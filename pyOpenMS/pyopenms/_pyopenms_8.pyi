from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_IonIdentityMolecularNetworking_annotateConsensusMap(consensus_map: ConsensusMap ) -> None:
    """
    Cython signature: void annotateConsensusMap(ConsensusMap & consensus_map)
        Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
        
        Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
        and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
        This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
        If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
        
        
        :param consensus_map: Input ConsensusMap without IIMN annotations.
    """
    ...

def __static_InternalCalibration_applyTransformation(pcs: List[Precursor] , trafo: MZTrafoModel ) -> None:
    """
    Cython signature: void applyTransformation(libcpp_vector[Precursor] & pcs, MZTrafoModel & trafo)
    """
    ...

def __static_InternalCalibration_applyTransformation(spec: MSSpectrum , target_mslvl: List[int] , trafo: MZTrafoModel ) -> None:
    """
    Cython signature: void applyTransformation(MSSpectrum & spec, IntList & target_mslvl, MZTrafoModel & trafo)
    """
    ...

def __static_InternalCalibration_applyTransformation(exp: MSExperiment , target_mslvl: List[int] , trafo: MZTrafoModel ) -> None:
    """
    Cython signature: void applyTransformation(MSExperiment & exp, IntList & target_mslvl, MZTrafoModel & trafo)
    """
    ...

def __static_FileHandler_computeFileHash(filename: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String computeFileHash(const String & filename)
    """
    ...

def __static_FileHandler_getType(filename: Union[bytes, str, String] ) -> int:
    """
    Cython signature: int getType(const String & filename)
    """
    ...

def __static_FileHandler_getTypeByContent(filename: Union[bytes, str, String] ) -> int:
    """
    Cython signature: FileType getTypeByContent(const String & filename)
    """
    ...

def __static_FileHandler_getTypeByFileName(filename: Union[bytes, str, String] ) -> int:
    """
    Cython signature: FileType getTypeByFileName(const String & filename)
    """
    ...

def __static_FileHandler_hasValidExtension(filename: Union[bytes, str, String] , type_: int ) -> bool:
    """
    Cython signature: bool hasValidExtension(const String & filename, FileType type_)
    """
    ...

def __static_FileHandler_isSupported(type_: int ) -> bool:
    """
    Cython signature: bool isSupported(FileType type_)
    """
    ...

def __static_FileHandler_stripExtension(file: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String stripExtension(String file)
    """
    ...

def __static_FileHandler_swapExtension(filename: Union[bytes, str, String] , new_type: int ) -> Union[bytes, str, String]:
    """
    Cython signature: String swapExtension(String filename, FileType new_type)
    """
    ...

def __static_IonIdentityMolecularNetworking_writeSupplementaryPairTable(consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None:
    """
    Cython signature: void writeSupplementaryPairTable(const ConsensusMap & consensus_map, const String & output_file)
        Write supplementary pair table (csv file) from a ConsensusMap with edge annotations for connected features. Required for GNPS IIMN.
        
        The table contains the columns "ID 1" (row ID of first feature), "ID 2" (row ID of second feature), "EdgeType" (MS1/2 annotation),
        "Score" (the number of direct partners from both connected features) and "Annotation" (adducts and delta m/z between two connected features).
        
        
        :param consensus_map: Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        :param output_file: Output file path for the supplementary pair table.
    """
    ...


class Base64:
    """
    Cython implementation of _Base64

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Base64.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Base64()
        Class to encode and decode Base64, it supports two precisions 32 bit (float) and 64 bit (double).
        """
        ...
    
    @overload
    def __init__(self, in_0: Base64 ) -> None:
        """
        Cython signature: void Base64(Base64 &)
        """
        ...
    
    def encodeIntegers(self, in_: List[int] , to_byte_order: int , out: String , zlib_compression: bool ) -> None:
        """
        Cython signature: void encodeIntegers(libcpp_vector[int] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)
        Encodes a vector of integer point numbers to a Base64 string
        """
        ...
    
    def decodeIntegers(self, in_: Union[bytes, str, String] , from_byte_order: int , out: List[int] , zlib_compression: bool ) -> None:
        """
        Cython signature: void decodeIntegers(const String & in_, ByteOrder from_byte_order, libcpp_vector[int] & out, bool zlib_compression)
        Decodes a Base64 string to a vector of integer numbers
        """
        ...
    
    def encodeStrings(self, in_: List[bytes] , out: String , zlib_compression: bool ) -> None:
        """
        Cython signature: void encodeStrings(libcpp_vector[String] & in_, String & out, bool zlib_compression)
        Encodes a vector of strings to a Base64 string
        """
        ...
    
    def decodeStrings(self, in_: Union[bytes, str, String] , out: List[bytes] , zlib_compression: bool ) -> None:
        """
        Cython signature: void decodeStrings(const String & in_, libcpp_vector[String] & out, bool zlib_compression)
        Decodes a Base64 string to a vector of (null-terminated) strings
        """
        ...
    ByteOrder : __ByteOrder 


class CVTermList:
    """
    Cython implementation of _CVTermList

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTermList.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVTermList()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVTermList ) -> None:
        """
        Cython signature: void CVTermList(CVTermList &)
        """
        ...
    
    def setCVTerms(self, terms: List[CVTerm] ) -> None:
        """
        Cython signature: void setCVTerms(libcpp_vector[CVTerm] & terms)
        Sets the CV terms
        """
        ...
    
    def replaceCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void replaceCVTerm(CVTerm & term)
        Replaces the specified CV term
        """
        ...
    
    def replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession)
        """
        ...
    
    def consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void consumeCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] cv_term_map)
        Merges the given map into the member map, no duplicate checking
        """
        ...
    
    def getCVTerms(self) -> Dict[bytes,List[CVTerm]]:
        """
        Cython signature: libcpp_map[String,libcpp_vector[CVTerm]] getCVTerms()
        Returns the accession string of the term
        """
        ...
    
    def addCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void addCVTerm(CVTerm & term)
        Adds a CV term
        """
        ...
    
    def hasCVTerm(self, accession: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasCVTerm(String accession)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def isMetaEmpty(self) -> bool:
        """
        Cython signature: bool isMetaEmpty()
        Returns if the MetaInfo is empty
        """
        ...
    
    def clearMetaInfo(self) -> None:
        """
        Cython signature: void clearMetaInfo()
        Removes all meta values
        """
        ...
    
    def metaRegistry(self) -> MetaInfoRegistry:
        """
        Cython signature: MetaInfoRegistry metaRegistry()
        Returns a reference to the MetaInfoRegistry
        """
        ...
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        Fills the given vector with a list of all keys for which a value is set
        """
        ...
    
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        Returns the value corresponding to a string, or
        """
        ...
    
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        Sets the DataValue corresponding to a name
        """
        ...
    
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        Returns whether an entry with the given name exists
        """
        ...
    
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        Removes the DataValue corresponding to `name` if it exists
        """
        ...
    
    def __richcmp__(self, other: CVTermList, op: int) -> Any:
        ... 


class DataFilter:
    """
    Cython implementation of _DataFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataFilter.html>`_
    """
    
    field: int
    
    op: int
    
    value: float
    
    value_string: Union[bytes, str, String]
    
    meta_name: Union[bytes, str, String]
    
    value_is_numerical: bool
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DataFilter()
        """
        ...
    
    @overload
    def __init__(self, in_0: DataFilter ) -> None:
        """
        Cython signature: void DataFilter(DataFilter &)
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    def fromString(self, filter_: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void fromString(const String & filter_)
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    def __richcmp__(self, other: DataFilter, op: int) -> Any:
        ... 


class DataFilters:
    """
    Cython implementation of _DataFilters

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataFilters.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DataFilters()
        """
        ...
    
    @overload
    def __init__(self, in_0: DataFilters ) -> None:
        """
        Cython signature: void DataFilters(DataFilters &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def __getitem__(self, in_0: int ) -> DataFilter:
        """
        Cython signature: DataFilter operator[](size_t)
        """
        ...
    
    def add(self, filter_: DataFilter ) -> None:
        """
        Cython signature: void add(DataFilter & filter_)
        """
        ...
    
    def remove(self, index: int ) -> None:
        """
        Cython signature: void remove(size_t index)
        """
        ...
    
    def replace(self, index: int , filter_: DataFilter ) -> None:
        """
        Cython signature: void replace(size_t index, DataFilter & filter_)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def setActive(self, is_active: bool ) -> None:
        """
        Cython signature: void setActive(bool is_active)
        """
        ...
    
    def isActive(self) -> bool:
        """
        Cython signature: bool isActive()
        """
        ...
    
    @overload
    def passes(self, feature: Feature ) -> bool:
        """
        Cython signature: bool passes(Feature & feature)
        """
        ...
    
    @overload
    def passes(self, consensus_feature: ConsensusFeature ) -> bool:
        """
        Cython signature: bool passes(ConsensusFeature & consensus_feature)
        """
        ...
    
    @overload
    def passes(self, spectrum: MSSpectrum , peak_index: int ) -> bool:
        """
        Cython signature: bool passes(MSSpectrum & spectrum, size_t peak_index)
        """
        ...
    FilterOperation : __FilterOperation
    FilterType : __FilterType 


class FeatureFinderAlgorithmPicked:
    """
    Cython implementation of _FeatureFinderAlgorithmPicked

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderAlgorithmPicked.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureFinderAlgorithmPicked()
        """
        ...
    
    def run(self, input_map: MSExperiment , output: FeatureMap , param: Param , seeds: FeatureMap ) -> None:
        """
        Cython signature: void run(MSExperiment & input_map, FeatureMap & output, Param & param, FeatureMap & seeds)
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class FileHandler:
    """
    Cython implementation of _FileHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FileHandler.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FileHandler()
        """
        ...
    
    def loadExperiment(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void loadExperiment(String, MSExperiment &)
        Loads a file into an MSExperiment
        
        
        :param filename: The file name of the file to load
        :param exp: The experiment to load the data into
        :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
        :param log: Progress logging mode
        :param rewrite_source_file: Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from
        :param compute_hash: If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile
        :return: true if the file could be loaded, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def storeExperiment(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void storeExperiment(String, MSExperiment)
        Stores an MSExperiment to a file\n
        
        The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used
        
        
        :param filename: The name of the file to store the data in
        :param exp: The experiment to store
        :param log: Progress logging mode
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be written
        """
        ...
    
    def loadFeatures(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None:
        """
        Cython signature: void loadFeatures(String, FeatureMap &)
        Loads a file into a FeatureMap
        
        
        :param filename: The file name of the file to load
        :param map: The FeatureMap to load the data into
        :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
        :return: true if the file could be loaded, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
        Access to the options for loading/storing
        """
        ...
    
    def setOptions(self, in_0: PeakFileOptions ) -> None:
        """
        Cython signature: void setOptions(PeakFileOptions)
        Sets options for loading/storing
        """
        ...
    
    computeFileHash: __static_FileHandler_computeFileHash
    
    getType: __static_FileHandler_getType
    
    getTypeByContent: __static_FileHandler_getTypeByContent
    
    getTypeByFileName: __static_FileHandler_getTypeByFileName
    
    hasValidExtension: __static_FileHandler_hasValidExtension
    
    isSupported: __static_FileHandler_isSupported
    
    stripExtension: __static_FileHandler_stripExtension
    
    swapExtension: __static_FileHandler_swapExtension 


class InternalCalibration:
    """
    Cython implementation of _InternalCalibration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InternalCalibration.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void InternalCalibration()
        A mass recalibration method using linear/quadratic interpolation (robust/weighted) of given reference masses
        """
        ...
    
    @overload
    def __init__(self, in_0: InternalCalibration ) -> None:
        """
        Cython signature: void InternalCalibration(InternalCalibration &)
        """
        ...
    
    @overload
    def fillCalibrants(self, in_0: MSExperiment , in_1: List[InternalCalibration_LockMass] , tol_ppm: float , lock_require_mono: bool , lock_require_iso: bool , failed_lock_masses: CalibrationData , verbose: bool ) -> int:
        """
        Cython signature: size_t fillCalibrants(MSExperiment, libcpp_vector[InternalCalibration_LockMass], double tol_ppm, bool lock_require_mono, bool lock_require_iso, CalibrationData & failed_lock_masses, bool verbose)
        Extract calibrants from Raw data (mzML)\n
        
        Lock masses are searched in each spectrum and added to the internal calibrant database\n
        
        Filters can be used to exclude spurious peaks, i.e. require the calibrant peak to be monoisotopic or
        to have a +1 isotope (should not be used for very low abundant calibrants)
        If a calibrant is not found, it is added to a 'failed_lock_masses' database which is returned and not stored internally.
        The intensity of the peaks describe the reason for failed detection: 0.0 - peak not found with the given ppm tolerance;
        1.0 - peak is not monoisotopic (can only occur if 'lock_require_mono' is true)
        2.0 - peak has no +1 isotope (can only occur if 'lock_require_iso' is true)
        
        
        :param exp: Peak map containing the lock masses
        :param ref_masses: List of lock masses
        :param tol_ppm: Search window for lock masses in 'exp'
        :param lock_require_mono: Require that a lock mass is the monoisotopic peak (i.e. not an isotope peak) -- lock mass is rejected otherwise
        :param lock_require_iso: Require that a lock mass has isotope peaks to its right -- lock mass is rejected otherwise
        :param failed_lock_masses: Set of calibration masses which were not found, i.e. their expected m/z and RT positions
        :param verbose: Print information on 'lock_require_XXX' matches during search
        :return: Number of calibration masses found
        """
        ...
    
    @overload
    def fillCalibrants(self, in_0: FeatureMap , in_1: float ) -> int:
        """
        Cython signature: size_t fillCalibrants(FeatureMap, double)
        Extract calibrants from identifications\n
        
        Extracts only the first hit from the first peptide identification of each feature
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter
        RT and m/z are taken from the features, not from the identifications (for an exception see below)!\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param fm: FeatureMap with peptide identifications
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
        """
        ...
    
    @overload
    def fillCalibrants(self, in_0: PeptideIdentificationList , in_1: float ) -> int:
        """
        Cython signature: size_t fillCalibrants(PeptideIdentificationList, double)
        Extract calibrants from identifications\n
        
        Extracts only the first hit from each peptide identification
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param pep_ids: Peptide ids (e.g. from an idXML file)
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
        """
        ...
    
    def getCalibrationPoints(self) -> CalibrationData:
        """
        Cython signature: CalibrationData getCalibrationPoints()
        Get container of calibration points\n
        
        Filled using fillCalibrants() methods
        
        
        :return: Container of calibration points
        """
        ...
    
    def calibrate(self, in_0: MSExperiment , in_1: List[int] , in_2: int , rt_chunk: float , use_RANSAC: bool , post_ppm_median: float , post_ppm_MAD: float , file_models: Union[bytes, str, String] , file_models_plot: Union[bytes, str, String] , file_residuals: Union[bytes, str, String] , file_residuals_plot: Union[bytes, str, String] , rscript_executable: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool calibrate(MSExperiment, libcpp_vector[int], MZTrafoModel_MODELTYPE, double rt_chunk, bool use_RANSAC, double post_ppm_median, double post_ppm_MAD, String file_models, String file_models_plot, String file_residuals, String file_residuals_plot, String rscript_executable)
        Apply calibration to data\n
        
        For each spectrum, a calibration model will be computed and applied.
        Make sure to call fillCalibrants() before, so a model can be created.\n
        
        The MSExperiment will be sorted by RT and m/z if unsorted.
        
        
        :param exp: MSExperiment holding the Raw data to calibrate
        :param target_mslvl: MS-levels where calibration should be applied to
        :param model_type: Linear or quadratic model; select based on your instrument
        :param rt_chunk: RT-window size (one-sided) of calibration points to collect around each spectrum. Set to negative values, to build one global model instead.
        :param use_RANSAC: Remove outliers before fitting a model?!
        :param post_ppm_median: The median ppm error of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
        :param post_ppm_MAD: The median absolute deviation of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
        :param file_models: Output CSV filename, where model parameters are written to (pass empty string to skip)
        :param file_models_plot: Output PNG image model parameters (pass empty string to skip)
        :param file_residuals: Output CSV filename, where ppm errors of calibrants before and after model fitting parameters are written to (pass empty string to skip)
        :param file_residuals_plot: Output PNG image of the ppm errors of calibrants (pass empty string to skip)
        :param rscript_executable: Full path to the Rscript executable
        :return: true upon successful calibration
        """
        ...
    
    def setLogType(self, in_0: int ) -> None:
        """
        Cython signature: void setLogType(LogType)
        Sets the progress log that should be used. The default type is NONE!
        """
        ...
    
    def getLogType(self) -> int:
        """
        Cython signature: LogType getLogType()
        Returns the type of progress log being used
        """
        ...
    
    def startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void startProgress(ptrdiff_t begin, ptrdiff_t end, String label)
        """
        ...
    
    def setProgress(self, value: int ) -> None:
        """
        Cython signature: void setProgress(ptrdiff_t value)
        Sets the current progress
        """
        ...
    
    def endProgress(self) -> None:
        """
        Cython signature: void endProgress()
        Ends the progress display
        """
        ...
    
    def nextProgress(self) -> None:
        """
        Cython signature: void nextProgress()
        Increment progress by 1 (according to range begin-end)
        """
        ...
    
    applyTransformation: __static_InternalCalibration_applyTransformation
    
    applyTransformation: __static_InternalCalibration_applyTransformation
    
    applyTransformation: __static_InternalCalibration_applyTransformation 


class InternalCalibration_LockMass:
    """
    Cython implementation of _InternalCalibration_LockMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InternalCalibration_LockMass.html>`_
    """
    
    mz: float
    
    ms_level: int
    
    charge: int
    
    @overload
    def __init__(self, mz_: float , lvl_: int , charge_: int ) -> None:
        """
        Cython signature: void InternalCalibration_LockMass(double mz_, int lvl_, int charge_)
        """
        ...
    
    @overload
    def __init__(self, in_0: InternalCalibration_LockMass ) -> None:
        """
        Cython signature: void InternalCalibration_LockMass(InternalCalibration_LockMass &)
        """
        ... 


class IonIdentityMolecularNetworking:
    """
    Cython implementation of _IonIdentityMolecularNetworking

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonIdentityMolecularNetworking.html>`_

    Includes the necessary functions to generate filed required for GNPS ion identity molecular networking (IIMN).
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void IonIdentityMolecularNetworking()
        """
        ...
    
    annotateConsensusMap: __static_IonIdentityMolecularNetworking_annotateConsensusMap
    
    writeSupplementaryPairTable: __static_IonIdentityMolecularNetworking_writeSupplementaryPairTable 


class MRMAssay:
    """
    Cython implementation of _MRMAssay

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMAssay.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMAssay()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMAssay ) -> None:
        """
        Cython signature: void MRMAssay(MRMAssay &)
        """
        ...
    
    def reannotateTransitions(self, exp: TargetedExperiment , precursor_mz_threshold: float , product_mz_threshold: float , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None:
        """
        Cython signature: void reannotateTransitions(TargetedExperiment & exp, double precursor_mz_threshold, double product_mz_threshold, libcpp_vector[String] fragment_types, libcpp_vector[size_t] fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow)
        Annotates and filters transitions in a TargetedExperiment
        
        
        :param exp: The input, unfiltered transitions
        :param precursor_mz_threshold: The precursor m/z threshold in Th for annotation
        :param product_mz_threshold: The product m/z threshold in Th for annotation
        :param fragment_types: The fragment types to consider for annotation
        :param fragment_charges: The fragment charges to consider for annotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param round_decPow: Round product m/z values to decimal power (default: -4)
        """
        ...
    
    def restrictTransitions(self, exp: TargetedExperiment , lower_mz_limit: float , upper_mz_limit: float , swathes: List[List[float, float]] ) -> None:
        """
        Cython signature: void restrictTransitions(TargetedExperiment & exp, double lower_mz_limit, double upper_mz_limit, libcpp_vector[libcpp_pair[double,double]] swathes)
        Restrict and filter transitions in a TargetedExperiment
        
        
        :param exp: The input, unfiltered transitions
        :param lower_mz_limit: The lower product m/z limit in Th
        :param upper_mz_limit: The upper product m/z limit in Th
        :param swathes: The swath window settings (to exclude fragment ions falling into the precursor isolation window)
        """
        ...
    
    def detectingTransitions(self, exp: TargetedExperiment , min_transitions: int , max_transitions: int ) -> None:
        """
        Cython signature: void detectingTransitions(TargetedExperiment & exp, int min_transitions, int max_transitions)
        Select detecting fragment ions
        
        
        :param exp: The input, unfiltered transitions
        :param min_transitions: The minimum number of transitions required per assay
        :param max_transitions: The maximum number of transitions required per assay
        """
        ...
    
    def filterMinMaxTransitionsCompound(self, exp: TargetedExperiment , min_transitions: int , max_transitions: int ) -> None:
        """
        Cython signature: void filterMinMaxTransitionsCompound(TargetedExperiment & exp, int min_transitions, int max_transitions)
        Filters target and decoy transitions by intensity, only keeping the top N transitions
        
        
        :param exp: The transition list which will be filtered
        :param min_transitions: The minimum number of transitions required per assay (targets only)
        :param max_transitions: The maximum number of transitions allowed per assay
        """
        ...
    
    def filterUnreferencedDecoysCompound(self, exp: TargetedExperiment ) -> None:
        """
        Cython signature: void filterUnreferencedDecoysCompound(TargetedExperiment & exp)
        Filters decoy transitions, which do not have respective target transition
        based on the transitionID.
        
        References between targets and decoys will be constructed based on the transitionsID
        and the "_decoy_" string. For example:
        
        target: 84_CompoundName_[M+H]+_88_22
        decoy: 84_CompoundName_decoy_[M+H]+_88_22
        
        
        :param exp: The transition list which will be filtered
        """
        ...
    
    def uisTransitions(self, exp: TargetedExperiment , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , enable_ms2_precursors: bool , mz_threshold: float , swathes: List[List[float, float]] , round_decPow: int , max_num_alternative_localizations: int , shuffle_seed: int ) -> None:
        """
        Cython signature: void uisTransitions(TargetedExperiment & exp, libcpp_vector[String] fragment_types, libcpp_vector[size_t] fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, bool enable_ms2_precursors, double mz_threshold, libcpp_vector[libcpp_pair[double,double]] swathes, int round_decPow, size_t max_num_alternative_localizations, int shuffle_seed)
        Annotate UIS / site-specific transitions
        
        Performs the following actions:
        
        - Step 1: For each peptide, compute all theoretical alternative peptidoforms; see transitions generateTargetInSilicoMap_()
        - Step 2: Generate target identification transitions; see generateTargetAssays_()
        
        - Step 3a: Generate decoy sequences that share peptidoform properties with targets; see generateDecoySequences_()
        - Step 3b: Generate decoy in silico peptide map containing theoretical transition; see generateDecoyInSilicoMap_()
        - Step 4: Generate decoy identification transitions; see generateDecoyAssays_()
        
        The IPF algorithm uses the concept of "identification transitions" that
        are used to discriminate different peptidoforms, these are generated in
        this function.  In brief, the algorithm takes the existing set of
        peptides and transitions and then appends these "identification
        transitions" for targets and decoys. The novel transitions are set to be
        non-detecting and non-quantifying and are annotated with the set of
        peptidoforms to which they map.
        
        
        :param exp: The input, unfiltered transitions
        :param fragment_types: The fragment types to consider for annotation
        :param fragment_charges: The fragment charges to consider for annotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param enable_ms2_precursors: Whether MS2 precursors should be considered
        :param mz_threshold: The product m/z threshold in Th for annotation
        :param swathes: The swath window settings (to exclude fragment ions falling
        :param round_decPow: Round product m/z values to decimal power (default: -4)
        :param max_num_alternative_localizations: Maximum number of allowed peptide sequence permutations
        :param shuffle_seed: Set seed for shuffle (-1: select seed based on time)
        :param disable_decoy_transitions: Whether to disable generation of decoy UIS transitions
        """
        ...
    
    def setLogType(self, in_0: int ) -> None:
        """
        Cython signature: void setLogType(LogType)
        Sets the progress log that should be used. The default type is NONE!
        """
        ...
    
    def getLogType(self) -> int:
        """
        Cython signature: LogType getLogType()
        Returns the type of progress log being used
        """
        ...
    
    def startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void startProgress(ptrdiff_t begin, ptrdiff_t end, String label)
        """
        ...
    
    def setProgress(self, value: int ) -> None:
        """
        Cython signature: void setProgress(ptrdiff_t value)
        Sets the current progress
        """
        ...
    
    def endProgress(self) -> None:
        """
        Cython signature: void endProgress()
        Ends the progress display
        """
        ...
    
    def nextProgress(self) -> None:
        """
        Cython signature: void nextProgress()
        Increment progress by 1 (according to range begin-end)
        """
        ... 


class MRMFeaturePickerFile:
    """
    Cython implementation of _MRMFeaturePickerFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeaturePickerFile.html>`_

    _MRMFeaturePickerFile_ loads components and components groups parameters from a .csv file
    
    The structures defined in [MRMFeaturePicker](@ref MRMFeaturePicker) are used
    
    It is required that columns `component_name` and `component_group_name` are present.
    Lines whose `component_name`'s or `component_group_name`'s value is an empty string, will be skipped.
    The class supports the absence of information within other columns.
    
    A reduced example of the expected format (fewer columns are shown here):
    > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
    > arg-L.arg-L_1.Heavy,arg-L,2,15
    > arg-L.arg-L_1.Light,arg-L,2,17
    > orn.orn_1.Heavy,orn,3,21
    > orn.orn_1.Light,orn,3,13
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeaturePickerFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeaturePickerFile ) -> None:
        """
        Cython signature: void MRMFeaturePickerFile(MRMFeaturePickerFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , cp_list: List[MRMFP_ComponentParams] , cgp_list: List[MRMFP_ComponentGroupParams] ) -> None:
        """
        Cython signature: void load(const String & filename, libcpp_vector[MRMFP_ComponentParams] & cp_list, libcpp_vector[MRMFP_ComponentGroupParams] & cgp_list)
        Loads the file's data and saves it into vectors of `ComponentParams` and `ComponentGroupParams`
        
        The file is expected to contain at least two columns: `component_name` and `component_group_name`. Otherwise,
        an exception is thrown
        
        If a component group (identified by its name) is found multiple times, only the first one is saved
        
        
        :param filename: Path to the .csv input file
        :param cp_list: Component params are saved in this list
        :param cgp_list: Component Group params are saved in this list
        :raises:
          Exception: MissingInformation If the required columns are not found
        :raises:
          Exception: FileNotFound If input file is not found
        """
        ... 


class MSDataAggregatingConsumer:
    """
    Cython implementation of _MSDataAggregatingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataAggregatingConsumer.html>`_
    """
    
    def __init__(self, in_0: MSDataAggregatingConsumer ) -> None:
        """
        Cython signature: void MSDataAggregatingConsumer(MSDataAggregatingConsumer &)
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        """
        ...
    
    def consumeChromatogram(self, in_0: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram &)
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        """
        ... 


class MultiplexDeltaMasses:
    """
    Cython implementation of _MultiplexDeltaMasses

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMasses.html>`_

    Data structure for mass shift pattern
    
    Groups of labelled peptides appear with characteristic mass shifts
    
    For example, for an Arg6 labeled SILAC peptide pair we expect to see
    mass shifts of 0 and 6 Da. Or as second example, for a
    peptide pair of a dimethyl labelled sample with a single lysine
    we will see mass shifts of 56 Da and 64 Da.
    28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K)
    for light and heavy partners respectively
    
    The data structure stores the mass shifts and corresponding labels
    for a group of matching peptide features
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MultiplexDeltaMasses()
        """
        ...
    
    @overload
    def __init__(self, in_0: MultiplexDeltaMasses ) -> None:
        """
        Cython signature: void MultiplexDeltaMasses(MultiplexDeltaMasses &)
        """
        ...
    
    @overload
    def __init__(self, dm: List[MultiplexDeltaMasses_DeltaMass] ) -> None:
        """
        Cython signature: void MultiplexDeltaMasses(libcpp_vector[MultiplexDeltaMasses_DeltaMass] & dm)
        """
        ...
    
    def getDeltaMasses(self) -> List[MultiplexDeltaMasses_DeltaMass]:
        """
        Cython signature: libcpp_vector[MultiplexDeltaMasses_DeltaMass] getDeltaMasses()
        """
        ... 


class MultiplexDeltaMasses_DeltaMass:
    """
    Cython implementation of _MultiplexDeltaMasses_DeltaMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMasses_DeltaMass.html>`_
    """
    
    delta_mass: float
    
    @overload
    def __init__(self, dm: float , l: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void MultiplexDeltaMasses_DeltaMass(double dm, String l)
        """
        ...
    
    @overload
    def __init__(self, in_0: MultiplexDeltaMasses_DeltaMass ) -> None:
        """
        Cython signature: void MultiplexDeltaMasses_DeltaMass(MultiplexDeltaMasses_DeltaMass &)
        """
        ... 


class Precursor:
    """
    Cython implementation of _Precursor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Precursor.html>`_
      -- Inherits from ['Peak1D', 'CVTermList']

    Precursor meta information
    
    This class contains precursor information:
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Precursor()
        """
        ...
    
    @overload
    def __init__(self, in_0: Precursor ) -> None:
        """
        Cython signature: void Precursor(Precursor &)
        """
        ...
    
    def getActivationMethods(self) -> Set[int]:
        """
        Cython signature: libcpp_set[ActivationMethod] getActivationMethods()
        Returns the activation methods
        """
        ...
    
    def getActivationMethodsAsString(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getActivationMethodsAsString()
        Returns the full names (e.g., "Collision-induced dissociation") of the activation methods set on this instance
        """
        ...
    
    def getActivationMethodsAsShortString(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getActivationMethodsAsShortString()
        Returns the abbreviations (e.g., "CID") of the activation methods set on this instance
        """
        ...
    
    def setActivationMethods(self, activation_methods: Set[int] ) -> None:
        """
        Cython signature: void setActivationMethods(libcpp_set[ActivationMethod] activation_methods)
        Sets the activation methods
        """
        ...
    
    @staticmethod
    def getAllNamesOfActivationMethods() -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getAllNamesOfActivationMethods()
        Returns the full names (e.g., "Collision-induced dissociation") of ALL possible activation methods, not just those set on this instance
        """
        ...
    
    @staticmethod
    def getAllShortNamesOfActivationMethods() -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getAllShortNamesOfActivationMethods()
        Returns the abbreviations (e.g., "CID") of ALL possible activation methods, not just those set on this instance
        """
        ...
    
    def getActivationEnergy(self) -> float:
        """
        Cython signature: double getActivationEnergy()
        Returns the activation energy (in electronvolt)
        """
        ...
    
    def setActivationEnergy(self, activation_energy: float ) -> None:
        """
        Cython signature: void setActivationEnergy(double activation_energy)
        Sets the activation energy (in electronvolt)
        """
        ...
    
    def getIsolationWindowLowerOffset(self) -> float:
        """
        Cython signature: double getIsolationWindowLowerOffset()
        Returns the lower offset from the target m/z
        """
        ...
    
    def setIsolationWindowLowerOffset(self, bound: float ) -> None:
        """
        Cython signature: void setIsolationWindowLowerOffset(double bound)
        Sets the lower offset from the target m/z
        """
        ...
    
    def getDriftTime(self) -> float:
        """
        Cython signature: double getDriftTime()
        Returns the ion mobility drift time in milliseconds (-1 means it is not set)
        """
        ...
    
    def setDriftTime(self, drift_time: float ) -> None:
        """
        Cython signature: void setDriftTime(double drift_time)
        Sets the ion mobility drift time in milliseconds
        """
        ...
    
    def getIsolationWindowUpperOffset(self) -> float:
        """
        Cython signature: double getIsolationWindowUpperOffset()
        Returns the upper offset from the target m/z
        """
        ...
    
    def setIsolationWindowUpperOffset(self, bound: float ) -> None:
        """
        Cython signature: void setIsolationWindowUpperOffset(double bound)
        Sets the upper offset from the target m/z
        """
        ...
    
    def getDriftTimeWindowLowerOffset(self) -> float:
        """
        Cython signature: double getDriftTimeWindowLowerOffset()
        Returns the lower offset from the target ion mobility in milliseconds
        """
        ...
    
    def setDriftTimeWindowLowerOffset(self, drift_time: float ) -> None:
        """
        Cython signature: void setDriftTimeWindowLowerOffset(double drift_time)
        Sets the lower offset from the target ion mobility
        """
        ...
    
    def getDriftTimeWindowUpperOffset(self) -> float:
        """
        Cython signature: double getDriftTimeWindowUpperOffset()
        Returns the upper offset from the target ion mobility in milliseconds
        """
        ...
    
    def setDriftTimeWindowUpperOffset(self, drift_time: float ) -> None:
        """
        Cython signature: void setDriftTimeWindowUpperOffset(double drift_time)
        Sets the upper offset from the target ion mobility
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        Returns the charge
        """
        ...
    
    def setCharge(self, charge: int ) -> None:
        """
        Cython signature: void setCharge(int charge)
        Sets the charge
        """
        ...
    
    def getPossibleChargeStates(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getPossibleChargeStates()
        Returns the possible charge states
        """
        ...
    
    def setPossibleChargeStates(self, possible_charge_states: List[int] ) -> None:
        """
        Cython signature: void setPossibleChargeStates(libcpp_vector[int] possible_charge_states)
        Sets the possible charge states
        """
        ...
    
    def getUnchargedMass(self) -> float:
        """
        Cython signature: double getUnchargedMass()
        Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
        """
        ...
    
    def getIntensity(self) -> float:
        """
        Cython signature: float getIntensity()
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        """
        ...
    
    def setMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setMZ(double)
        """
        ...
    
    def setIntensity(self, in_0: float ) -> None:
        """
        Cython signature: void setIntensity(float)
        """
        ...
    
    def getPos(self) -> float:
        """
        Cython signature: double getPos()
        """
        ...
    
    def setPos(self, pos: float ) -> None:
        """
        Cython signature: void setPos(double pos)
        """
        ...
    
    def setCVTerms(self, terms: List[CVTerm] ) -> None:
        """
        Cython signature: void setCVTerms(libcpp_vector[CVTerm] & terms)
        Sets the CV terms
        """
        ...
    
    def replaceCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void replaceCVTerm(CVTerm & term)
        Replaces the specified CV term
        """
        ...
    
    def replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession)
        """
        ...
    
    def consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void consumeCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] cv_term_map)
        Merges the given map into the member map, no duplicate checking
        """
        ...
    
    def getCVTerms(self) -> Dict[bytes,List[CVTerm]]:
        """
        Cython signature: libcpp_map[String,libcpp_vector[CVTerm]] getCVTerms()
        Returns the accession string of the term
        """
        ...
    
    def addCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void addCVTerm(CVTerm & term)
        Adds a CV term
        """
        ...
    
    def hasCVTerm(self, accession: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasCVTerm(String accession)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def isMetaEmpty(self) -> bool:
        """
        Cython signature: bool isMetaEmpty()
        Returns if the MetaInfo is empty
        """
        ...
    
    def clearMetaInfo(self) -> None:
        """
        Cython signature: void clearMetaInfo()
        Removes all meta values
        """
        ...
    
    def metaRegistry(self) -> MetaInfoRegistry:
        """
        Cython signature: MetaInfoRegistry metaRegistry()
        Returns a reference to the MetaInfoRegistry
        """
        ...
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        Fills the given vector with a list of all keys for which a value is set
        """
        ...
    
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        Returns the value corresponding to a string, or
        """
        ...
    
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        Sets the DataValue corresponding to a name
        """
        ...
    
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        Returns whether an entry with the given name exists
        """
        ...
    
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        Removes the DataValue corresponding to `name` if it exists
        """
        ...
    
    def __richcmp__(self, other: Precursor, op: int) -> Any:
        ...
    ActivationMethod : __ActivationMethod 


class TargetedExperiment:
    """
    Cython implementation of _TargetedExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TargetedExperiment.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TargetedExperiment()
        """
        ...
    
    @overload
    def __init__(self, in_0: TargetedExperiment ) -> None:
        """
        Cython signature: void TargetedExperiment(TargetedExperiment &)
        """
        ...
    
    def __add__(self: TargetedExperiment, other: TargetedExperiment) -> TargetedExperiment:
        ...
    
    def __iadd__(self: TargetedExperiment, other: TargetedExperiment) -> TargetedExperiment:
        ...
    
    def clear(self, clear_meta_data: bool ) -> None:
        """
        Cython signature: void clear(bool clear_meta_data)
        """
        ...
    
    def sortTransitionsByProductMZ(self) -> None:
        """
        Cython signature: void sortTransitionsByProductMZ()
        """
        ...
    
    def setCVs(self, cvs: List[CV] ) -> None:
        """
        Cython signature: void setCVs(libcpp_vector[CV] cvs)
        """
        ...
    
    def getCVs(self) -> List[CV]:
        """
        Cython signature: libcpp_vector[CV] getCVs()
        """
        ...
    
    def addCV(self, cv: CV ) -> None:
        """
        Cython signature: void addCV(CV cv)
        """
        ...
    
    def setContacts(self, contacts: List[Contact] ) -> None:
        """
        Cython signature: void setContacts(libcpp_vector[Contact] contacts)
        """
        ...
    
    def getContacts(self) -> List[Contact]:
        """
        Cython signature: libcpp_vector[Contact] getContacts()
        """
        ...
    
    def addContact(self, contact: Contact ) -> None:
        """
        Cython signature: void addContact(Contact contact)
        """
        ...
    
    def setPublications(self, publications: List[Publication] ) -> None:
        """
        Cython signature: void setPublications(libcpp_vector[Publication] publications)
        """
        ...
    
    def getPublications(self) -> List[Publication]:
        """
        Cython signature: libcpp_vector[Publication] getPublications()
        """
        ...
    
    def addPublication(self, publication: Publication ) -> None:
        """
        Cython signature: void addPublication(Publication publication)
        """
        ...
    
    def setTargetCVTerms(self, cv_terms: CVTermList ) -> None:
        """
        Cython signature: void setTargetCVTerms(CVTermList cv_terms)
        """
        ...
    
    def getTargetCVTerms(self) -> CVTermList:
        """
        Cython signature: CVTermList getTargetCVTerms()
        """
        ...
    
    def addTargetCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void addTargetCVTerm(CVTerm cv_term)
        """
        ...
    
    def setTargetMetaValue(self, name: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setTargetMetaValue(String name, DataValue value)
        """
        ...
    
    def setInstruments(self, instruments: List[TargetedExperiment_Instrument] ) -> None:
        """
        Cython signature: void setInstruments(libcpp_vector[TargetedExperiment_Instrument] instruments)
        """
        ...
    
    def getInstruments(self) -> List[TargetedExperiment_Instrument]:
        """
        Cython signature: libcpp_vector[TargetedExperiment_Instrument] getInstruments()
        """
        ...
    
    def addInstrument(self, instrument: TargetedExperiment_Instrument ) -> None:
        """
        Cython signature: void addInstrument(TargetedExperiment_Instrument instrument)
        """
        ...
    
    def setSoftware(self, software: List[Software] ) -> None:
        """
        Cython signature: void setSoftware(libcpp_vector[Software] software)
        """
        ...
    
    def getSoftware(self) -> List[Software]:
        """
        Cython signature: libcpp_vector[Software] getSoftware()
        """
        ...
    
    def addSoftware(self, software: Software ) -> None:
        """
        Cython signature: void addSoftware(Software software)
        """
        ...
    
    def setProteins(self, proteins: List[Protein] ) -> None:
        """
        Cython signature: void setProteins(libcpp_vector[Protein] proteins)
        """
        ...
    
    def getProteins(self) -> List[Protein]:
        """
        Cython signature: libcpp_vector[Protein] getProteins()
        """
        ...
    
    def getProteinByRef(self, ref: Union[bytes, str, String] ) -> Protein:
        """
        Cython signature: Protein getProteinByRef(String ref)
        """
        ...
    
    def hasProtein(self, ref: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasProtein(String ref)
        """
        ...
    
    def addProtein(self, protein: Protein ) -> None:
        """
        Cython signature: void addProtein(Protein protein)
        """
        ...
    
    def setCompounds(self, rhs: List[Compound] ) -> None:
        """
        Cython signature: void setCompounds(libcpp_vector[Compound] rhs)
        """
        ...
    
    def getCompounds(self) -> List[Compound]:
        """
        Cython signature: libcpp_vector[Compound] getCompounds()
        """
        ...
    
    def addCompound(self, rhs: Compound ) -> None:
        """
        Cython signature: void addCompound(Compound rhs)
        """
        ...
    
    def hasCompound(self, ref: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasCompound(String ref)
        """
        ...
    
    def getCompoundByRef(self, ref: Union[bytes, str, String] ) -> Compound:
        """
        Cython signature: Compound getCompoundByRef(String ref)
        """
        ...
    
    def setPeptides(self, rhs: List[Peptide] ) -> None:
        """
        Cython signature: void setPeptides(libcpp_vector[Peptide] rhs)
        """
        ...
    
    def getPeptides(self) -> List[Peptide]:
        """
        Cython signature: libcpp_vector[Peptide] getPeptides()
        """
        ...
    
    def hasPeptide(self, ref: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasPeptide(String ref)
        """
        ...
    
    def getPeptideByRef(self, ref: Union[bytes, str, String] ) -> Peptide:
        """
        Cython signature: Peptide getPeptideByRef(String ref)
        """
        ...
    
    def addPeptide(self, rhs: Peptide ) -> None:
        """
        Cython signature: void addPeptide(Peptide rhs)
        """
        ...
    
    def setTransitions(self, transitions: List[ReactionMonitoringTransition] ) -> None:
        """
        Cython signature: void setTransitions(libcpp_vector[ReactionMonitoringTransition] transitions)
        """
        ...
    
    def getTransitions(self) -> List[ReactionMonitoringTransition]:
        """
        Cython signature: libcpp_vector[ReactionMonitoringTransition] getTransitions()
        """
        ...
    
    def addTransition(self, transition: ReactionMonitoringTransition ) -> None:
        """
        Cython signature: void addTransition(ReactionMonitoringTransition transition)
        """
        ...
    
    def setIncludeTargets(self, targets: List[IncludeExcludeTarget] ) -> None:
        """
        Cython signature: void setIncludeTargets(libcpp_vector[IncludeExcludeTarget] targets)
        """
        ...
    
    def getIncludeTargets(self) -> List[IncludeExcludeTarget]:
        """
        Cython signature: libcpp_vector[IncludeExcludeTarget] getIncludeTargets()
        """
        ...
    
    def addIncludeTarget(self, target: IncludeExcludeTarget ) -> None:
        """
        Cython signature: void addIncludeTarget(IncludeExcludeTarget target)
        """
        ...
    
    def setExcludeTargets(self, targets: List[IncludeExcludeTarget] ) -> None:
        """
        Cython signature: void setExcludeTargets(libcpp_vector[IncludeExcludeTarget] targets)
        """
        ...
    
    def getExcludeTargets(self) -> List[IncludeExcludeTarget]:
        """
        Cython signature: libcpp_vector[IncludeExcludeTarget] getExcludeTargets()
        """
        ...
    
    def addExcludeTarget(self, target: IncludeExcludeTarget ) -> None:
        """
        Cython signature: void addExcludeTarget(IncludeExcludeTarget target)
        """
        ...
    
    def setSourceFiles(self, source_files: List[SourceFile] ) -> None:
        """
        Cython signature: void setSourceFiles(libcpp_vector[SourceFile] source_files)
        """
        ...
    
    def getSourceFiles(self) -> List[SourceFile]:
        """
        Cython signature: libcpp_vector[SourceFile] getSourceFiles()
        """
        ...
    
    def addSourceFile(self, source_file: SourceFile ) -> None:
        """
        Cython signature: void addSourceFile(SourceFile source_file)
        """
        ...
    
    def containsInvalidReferences(self) -> bool:
        """
        Cython signature: bool containsInvalidReferences()
        """
        ...
    
    def __richcmp__(self, other: TargetedExperiment, op: int) -> Any:
        ... 


class WindowMower:
    """
    Cython implementation of _WindowMower

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1WindowMower.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void WindowMower()
        """
        ...
    
    @overload
    def __init__(self, in_0: WindowMower ) -> None:
        """
        Cython signature: void WindowMower(WindowMower &)
        """
        ...
    
    def filterPeakSpectrumForTopNInSlidingWindow(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrumForTopNInSlidingWindow(MSSpectrum & spectrum)
        Sliding window version (slower)
        """
        ...
    
    def filterPeakSpectrumForTopNInJumpingWindow(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrumForTopNInJumpingWindow(MSSpectrum & spectrum)
        Jumping window version (faster)
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class __ActivationMethod:
    """
        Enum for activation/fragmentation methods of mass spectra
    
    - CID: Collision-induced dissociation (MS:1000133) (also CAD; parent term, but unless otherwise stated often used as synonym for trap-type CID)
    - PSD: Post-source decay
    - PD: Plasma desorption
    - SID: Surface-induced dissociation
    - BIRD: Blackbody infrared radiative dissociation
    - ECD: Electron capture dissociation (MS:1000250)
    - IMD: Infrared multiphoton dissociation
    - SORI: Sustained off-resonance irradiation
    - HCID: High-energy collision-induced dissociation
    - LCID: Low-energy collision-induced dissociation
    - PHD: Photodissociation
    - ETD: Electron transfer dissociation
    - ETciD: Electron transfer and collision-induced dissociation (MS:1003182)
    - EThcD: Electron transfer and higher-energy collision dissociation (MS:1002631)
    - PQD: Pulsed q dissociation (MS:1000599)
    - TRAP: trap-type collision-induced dissociation (MS:1002472)
    - HCD: beam-type collision-induced dissociation (MS:1000422)
    - INSOURCE: in-source collision-induced dissociation (MS:1001880)
    - LIFT: Bruker proprietary method (MS:1002000)
    """
    CID : int
    PSD : int
    PD : int
    SID : int
    BIRD : int
    ECD : int
    IMD : int
    SORI : int
    HCID : int
    LCID : int
    PHD : int
    ETD : int
    ETciD : int
    EThcD : int
    PQD : int
    TRAP : int
    HCD : int
    INSOURCE : int
    LIFT : int
    SIZE_OF_ACTIVATIONMETHOD : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ByteOrder:
    None
    BYTEORDER_BIGENDIAN : int
    BYTEORDER_LITTLEENDIAN : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __FilterOperation:
    None
    GREATER_EQUAL : int
    EQUAL : int
    LESS_EQUAL : int
    EXISTS : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __FilterType:
    None
    INTENSITY : int
    QUALITY : int
    CHARGE : int
    SIZE : int
    META_DATA : int

    def getMapping(self) -> Dict[int, str]:
       ... 

