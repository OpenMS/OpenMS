from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_PrecursorCorrection_correctToHighestIntensityMS1Peak(exp: MSExperiment , mz_tolerance: float , ppm: bool , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> Set[int]:
    """
    Cython signature: libcpp_set[size_t] correctToHighestIntensityMS1Peak(MSExperiment & exp, double mz_tolerance, bool ppm, libcpp_vector[double] & delta_mzs, libcpp_vector[double] & mzs, libcpp_vector[double] & rts)
    """
    ...

def __static_PrecursorCorrection_correctToNearestFeature(features: FeatureMap , exp: MSExperiment , rt_tolerance_s: float , mz_tolerance: float , ppm: bool , believe_charge: bool , keep_original: bool , all_matching_features: bool , max_trace: int , debug_level: int ) -> Set[int]:
    """
    Cython signature: libcpp_set[size_t] correctToNearestFeature(FeatureMap & features, MSExperiment & exp, double rt_tolerance_s, double mz_tolerance, bool ppm, bool believe_charge, bool keep_original, bool all_matching_features, int max_trace, int debug_level)
    """
    ...

def __static_PrecursorCorrection_correctToNearestMS1Peak(exp: MSExperiment , mz_tolerance: float , ppm: bool , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> Set[int]:
    """
    Cython signature: libcpp_set[size_t] correctToNearestMS1Peak(MSExperiment & exp, double mz_tolerance, bool ppm, libcpp_vector[double] & delta_mzs, libcpp_vector[double] & mzs, libcpp_vector[double] & rts)
    """
    ...

def __static_PrecursorCorrection_getPrecursors(exp: MSExperiment , precursors: List[Precursor] , precursors_rt: List[float] , precursor_scan_index: List[int] ) -> None:
    """
    Cython signature: void getPrecursors(MSExperiment & exp, libcpp_vector[Precursor] & precursors, libcpp_vector[double] & precursors_rt, libcpp_vector[size_t] & precursor_scan_index)
    """
    ...

def __static_PrecursorCorrection_writeHist(out_csv: String , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> None:
    """
    Cython signature: void writeHist(String & out_csv, libcpp_vector[double] & delta_mzs, libcpp_vector[double] & mzs, libcpp_vector[double] & rts)
    """
    ...


class DTA2DFile:
    """
    Cython implementation of _DTA2DFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DTA2DFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DTA2DFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: DTA2DFile ) -> None:
        """
        Cython signature: void DTA2DFile(DTA2DFile &)
        """
        ...
    
    def storeTIC(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None:
        """
        Cython signature: void storeTIC(String filename, MSExperiment & peakmap)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None:
        """
        Cython signature: void store(String filename, MSExperiment & peakmap)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None:
        """
        Cython signature: void load(String filename, MSExperiment & peakmap)
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
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


class DTAFile:
    """
    Cython implementation of _DTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DTAFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DTAFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: DTAFile ) -> None:
        """
        Cython signature: void DTAFile(DTAFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void load(String filename, MSSpectrum & spectrum)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void store(String filename, MSSpectrum & spectrum)
        """
        ... 


class FeatureHandle:
    """
    Cython implementation of _FeatureHandle

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureHandle.html>`_
      -- Inherits from ['Peak2D', 'UniqueIdInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureHandle()
        Representation of a Peak2D, RichPeak2D or Feature
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureHandle ) -> None:
        """
        Cython signature: void FeatureHandle(FeatureHandle &)
        """
        ...
    
    @overload
    def __init__(self, map_index: int , point: Peak2D , element_index: int ) -> None:
        """
        Cython signature: void FeatureHandle(uint64_t map_index, Peak2D & point, uint64_t element_index)
        """
        ...
    
    def getMapIndex(self) -> int:
        """
        Cython signature: uint64_t getMapIndex()
        Returns the map index
        """
        ...
    
    def setMapIndex(self, i: int ) -> None:
        """
        Cython signature: void setMapIndex(uint64_t i)
        Sets the map index
        """
        ...
    
    def setCharge(self, charge: int ) -> None:
        """
        Cython signature: void setCharge(int charge)
        Sets the charge
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        Returns the charge
        """
        ...
    
    def setWidth(self, width: float ) -> None:
        """
        Cython signature: void setWidth(float width)
        Sets the width (FWHM)
        """
        ...
    
    def getWidth(self) -> float:
        """
        Cython signature: float getWidth()
        Returns the width (FWHM)
        """
        ...
    
    def getIntensity(self) -> float:
        """
        Cython signature: float getIntensity()
        Returns the data point intensity (height)
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        Returns the m/z coordinate (index 1)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        Returns the RT coordinate (index 0)
        """
        ...
    
    def setMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setMZ(double)
        Returns the m/z coordinate (index 1)
        """
        ...
    
    def setRT(self, in_0: float ) -> None:
        """
        Cython signature: void setRT(double)
        Returns the RT coordinate (index 0)
        """
        ...
    
    def setIntensity(self, in_0: float ) -> None:
        """
        Cython signature: void setIntensity(float)
        Returns the data point intensity (height)
        """
        ...
    
    def getUniqueId(self) -> int:
        """
        Cython signature: size_t getUniqueId()
        Returns the unique id
        """
        ...
    
    def clearUniqueId(self) -> int:
        """
        Cython signature: size_t clearUniqueId()
        Clear the unique id. The new unique id will be invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        ...
    
    def hasValidUniqueId(self) -> int:
        """
        Cython signature: size_t hasValidUniqueId()
        Returns whether the unique id is valid. Returns 1 if the unique id is valid, 0 otherwise
        """
        ...
    
    def hasInvalidUniqueId(self) -> int:
        """
        Cython signature: size_t hasInvalidUniqueId()
        Returns whether the unique id is invalid. Returns 1 if the unique id is invalid, 0 otherwise
        """
        ...
    
    def setUniqueId(self, rhs: int ) -> None:
        """
        Cython signature: void setUniqueId(uint64_t rhs)
        Assigns a new, valid unique id. Always returns 1
        """
        ...
    
    def ensureUniqueId(self) -> int:
        """
        Cython signature: size_t ensureUniqueId()
        Assigns a valid unique id, but only if the present one is invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        ...
    
    def isValid(self, unique_id: int ) -> bool:
        """
        Cython signature: bool isValid(uint64_t unique_id)
        Returns true if the unique_id is valid, false otherwise
        """
        ...
    
    def __richcmp__(self, other: FeatureHandle, op: int) -> Any:
        ... 


class IncludeExcludeTarget:
    """
    Cython implementation of _IncludeExcludeTarget

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IncludeExcludeTarget.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IncludeExcludeTarget()
        This class stores a SRM/MRM transition
        """
        ...
    
    @overload
    def __init__(self, in_0: IncludeExcludeTarget ) -> None:
        """
        Cython signature: void IncludeExcludeTarget(IncludeExcludeTarget &)
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String & name)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def setPeptideRef(self, peptide_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPeptideRef(const String & peptide_ref)
        """
        ...
    
    def getPeptideRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPeptideRef()
        """
        ...
    
    def setNuctideRef(self, nuctide_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNuctideRef(const String & nuctide_ref)
        """
        ...
    
    def getNuctideRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNuctideRef()
        """
        ...
    
    def setCompoundRef(self, compound_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCompoundRef(const String & compound_ref)
        """
        ...
    
    def getCompoundRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCompoundRef()
        """
        ...
    
    def setPrecursorMZ(self, mz: float ) -> None:
        """
        Cython signature: void setPrecursorMZ(double mz)
        """
        ...
    
    def getPrecursorMZ(self) -> float:
        """
        Cython signature: double getPrecursorMZ()
        """
        ...
    
    def setPrecursorCVTermList(self, list_: CVTermList ) -> None:
        """
        Cython signature: void setPrecursorCVTermList(CVTermList & list_)
        """
        ...
    
    def addPrecursorCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void addPrecursorCVTerm(CVTerm & cv_term)
        """
        ...
    
    def getPrecursorCVTermList(self) -> CVTermList:
        """
        Cython signature: CVTermList getPrecursorCVTermList()
        """
        ...
    
    def setProductMZ(self, mz: float ) -> None:
        """
        Cython signature: void setProductMZ(double mz)
        """
        ...
    
    def getProductMZ(self) -> float:
        """
        Cython signature: double getProductMZ()
        """
        ...
    
    def setProductCVTermList(self, list_: CVTermList ) -> None:
        """
        Cython signature: void setProductCVTermList(CVTermList & list_)
        """
        ...
    
    def addProductCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void addProductCVTerm(CVTerm & cv_term)
        """
        ...
    
    def getProductCVTermList(self) -> CVTermList:
        """
        Cython signature: CVTermList getProductCVTermList()
        """
        ...
    
    def setInterpretations(self, interpretations: List[CVTermList] ) -> None:
        """
        Cython signature: void setInterpretations(libcpp_vector[CVTermList] & interpretations)
        """
        ...
    
    def getInterpretations(self) -> List[CVTermList]:
        """
        Cython signature: libcpp_vector[CVTermList] getInterpretations()
        """
        ...
    
    def addInterpretation(self, interpretation: CVTermList ) -> None:
        """
        Cython signature: void addInterpretation(CVTermList & interpretation)
        """
        ...
    
    def setConfigurations(self, configuration: List[Configuration] ) -> None:
        """
        Cython signature: void setConfigurations(libcpp_vector[Configuration] & configuration)
        """
        ...
    
    def getConfigurations(self) -> List[Configuration]:
        """
        Cython signature: libcpp_vector[Configuration] getConfigurations()
        """
        ...
    
    def addConfiguration(self, configuration: Configuration ) -> None:
        """
        Cython signature: void addConfiguration(Configuration & configuration)
        """
        ...
    
    def setPrediction(self, prediction: CVTermList ) -> None:
        """
        Cython signature: void setPrediction(CVTermList & prediction)
        """
        ...
    
    def addPredictionTerm(self, prediction: CVTerm ) -> None:
        """
        Cython signature: void addPredictionTerm(CVTerm & prediction)
        """
        ...
    
    def getPrediction(self) -> CVTermList:
        """
        Cython signature: CVTermList getPrediction()
        """
        ...
    
    def setRetentionTime(self, rt: RetentionTime ) -> None:
        """
        Cython signature: void setRetentionTime(RetentionTime rt)
        """
        ...
    
    def getRetentionTime(self) -> RetentionTime:
        """
        Cython signature: RetentionTime getRetentionTime()
        """
        ...
    
    def setCVTerms(self, terms: List[CVTerm] ) -> None:
        """
        Cython signature: void setCVTerms(libcpp_vector[CVTerm] & terms)
        """
        ...
    
    def replaceCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void replaceCVTerm(CVTerm & term)
        """
        ...
    
    @overload
    def replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_vector[CVTerm] cv_terms, String accession)
        """
        ...
    
    @overload
    def replaceCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] cv_term_map)
        """
        ...
    
    def getCVTerms(self) -> Dict[bytes,List[CVTerm]]:
        """
        Cython signature: libcpp_map[String,libcpp_vector[CVTerm]] getCVTerms()
        """
        ...
    
    def addCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void addCVTerm(CVTerm & term)
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
    
    def __richcmp__(self, other: IncludeExcludeTarget, op: int) -> Any:
        ... 


class IsobaricChannelExtractor:
    """
    Cython implementation of _IsobaricChannelExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricChannelExtractor.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, in_0: IsobaricChannelExtractor ) -> None:
        """
        Cython signature: void IsobaricChannelExtractor(IsobaricChannelExtractor &)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricChannelExtractor(ItraqEightPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricChannelExtractor(ItraqFourPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTSixPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricChannelExtractor(TMTSixPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTTenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricChannelExtractor(TMTTenPlexQuantitationMethod * quant_method)
        """
        ...
    
    def extractChannels(self, ms_exp_data: MSExperiment , consensus_map: ConsensusMap ) -> None:
        """
        Cython signature: void extractChannels(MSExperiment & ms_exp_data, ConsensusMap & consensus_map)
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


class KDTreeFeatureMaps:
    """
    Cython implementation of _KDTreeFeatureMaps

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KDTreeFeatureMaps.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void KDTreeFeatureMaps()
        Stores a set of features, together with a 2D tree for fast search
        """
        ...
    
    @overload
    def __init__(self, maps: List[FeatureMap] , param: Param ) -> None:
        """
        Cython signature: void KDTreeFeatureMaps(libcpp_vector[FeatureMap] & maps, Param & param)
        """
        ...
    
    @overload
    def __init__(self, maps: List[ConsensusMap] , param: Param ) -> None:
        """
        Cython signature: void KDTreeFeatureMaps(libcpp_vector[ConsensusMap] & maps, Param & param)
        """
        ...
    
    @overload
    def addMaps(self, maps: List[FeatureMap] ) -> None:
        """
        Cython signature: void addMaps(libcpp_vector[FeatureMap] & maps)
        Add `maps` and balance kd-tree
        """
        ...
    
    @overload
    def addMaps(self, maps: List[ConsensusMap] ) -> None:
        """
        Cython signature: void addMaps(libcpp_vector[ConsensusMap] & maps)
        """
        ...
    
    def rt(self, i: int ) -> float:
        """
        Cython signature: double rt(size_t i)
        """
        ...
    
    def mz(self, i: int ) -> float:
        """
        Cython signature: double mz(size_t i)
        """
        ...
    
    def intensity(self, i: int ) -> float:
        """
        Cython signature: float intensity(size_t i)
        """
        ...
    
    def charge(self, i: int ) -> int:
        """
        Cython signature: int charge(size_t i)
        """
        ...
    
    def mapIndex(self, i: int ) -> int:
        """
        Cython signature: size_t mapIndex(size_t i)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def treeSize(self) -> int:
        """
        Cython signature: size_t treeSize()
        """
        ...
    
    def numMaps(self) -> int:
        """
        Cython signature: size_t numMaps()
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def optimizeTree(self) -> None:
        """
        Cython signature: void optimizeTree()
        """
        ...
    
    def getNeighborhood(self, index: int , result_indices: List[int] , rt_tol: float , mz_tol: float , mz_ppm: bool , include_features_from_same_map: bool , max_pairwise_log_fc: float ) -> None:
        """
        Cython signature: void getNeighborhood(size_t index, libcpp_vector[size_t] & result_indices, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map, double max_pairwise_log_fc)
        Fill `result` with indices of all features compatible (wrt. RT, m/z, map index) to the feature with `index`
        """
        ...
    
    def queryRegion(self, rt_low: float , rt_high: float , mz_low: float , mz_high: float , result_indices: List[int] , ignored_map_index: int ) -> None:
        """
        Cython signature: void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, libcpp_vector[size_t] & result_indices, size_t ignored_map_index)
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


class KDTreeFeatureNode:
    """
    Cython implementation of _KDTreeFeatureNode

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KDTreeFeatureNode.html>`_

    A node of the kD-tree with pointer to corresponding data and index
    """
    
    @overload
    def __init__(self, in_0: KDTreeFeatureNode ) -> None:
        """
        Cython signature: void KDTreeFeatureNode(KDTreeFeatureNode &)
        """
        ...
    
    @overload
    def __init__(self, data: KDTreeFeatureMaps , idx: int ) -> None:
        """
        Cython signature: void KDTreeFeatureNode(KDTreeFeatureMaps * data, size_t idx)
        """
        ...
    
    def __getitem__(self, i: int ) -> float:
        """
        Cython signature: double operator[](size_t i)
        """
        ...
    
    def getIndex(self) -> int:
        """
        Cython signature: size_t getIndex()
        Returns index of corresponding feature in data_
        """
        ... 


class LogConfigHandler:
    """
    Cython implementation of _LogConfigHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LogConfigHandler.html>`_
    """
    
    def parse(self, setting: List[bytes] ) -> Param:
        """
        Cython signature: Param parse(const StringList & setting)
        Translates the given list of parameter settings into a LogStream configuration
        
        Translates the given list of parameter settings into a LogStream configuration.
        Usually this list stems from a command line call.
        
        Each element in the stringlist should follow this naming convention
        
        <LOG_NAME> <ACTION> <PARAMETER>
        
        with
        - LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
        - ACTION: add,remove,clear
        - PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter
        
        Example:
        `DEBUG add debug.log`
        
        This function will **not** apply to settings to the log handlers. Use configure() for that.
        
        :param setting: StringList containing the configuration options
        :raises ParseError: In case of an invalid configuration.
        :return: Param object containing all settings, that can be applied using the LogConfigHandler.configure() method
        """
        ...
    
    def configure(self, param: Param ) -> None:
        """
        Cython signature: void configure(const Param & param)
        Applies the given parameters (@p param) to the current configuration
        
        <LOG_NAME> <ACTION> <PARAMETER> <STREAMTYPE>
        
        LOG_NAME: DEBUG, INFO, WARNING, ERROR, FATAL_ERROR
        ACTION: add, remove, clear
        PARAMETER: for 'add'/'remove' it is the stream name ('cout', 'cerr' or a filename), 'clear' does not require any further parameter
        STREAMTYPE: FILE, STRING (for a StringStream, which you can grab by this name using getStream() )
        
        You cannot specify a file named "cout" or "cerr" even if you specify streamtype 'FILE' - the handler will mistake this for the
        internal streams, but you can use "./cout" to print to a file named cout.
        
        A classical configuration would contain a list of settings e.g.
        
        `DEBUG add debug.log FILE`
        `INFO remove cout FILE` (FILE will be ignored)
        `INFO add string_stream1 STRING`
        
        :raises ElementNotFound: If the LogStream (first argument) does not exist.
        :raises FileNotWritable: If a file (or stream) should be opened as log file (or stream) that is not accessible.
        :raises IllegalArgument: If a stream should be registered, that was already registered with a different type.
        """
        ...
    
    def setLogLevel(self, log_level: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLogLevel(const String & log_level)
        Sets a minimum log_level by removing all streams from loggers lower than that level.
        Valid levels are from low to high: "DEBUG", "INFO", "WARNING", "ERROR", "FATAL_ERROR"
        """
        ... 


class ModificationDefinition:
    """
    Cython implementation of _ModificationDefinition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationDefinition.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ModificationDefinition()
        """
        ...
    
    @overload
    def __init__(self, in_0: ModificationDefinition ) -> None:
        """
        Cython signature: void ModificationDefinition(ModificationDefinition &)
        """
        ...
    
    @overload
    def __init__(self, mod: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void ModificationDefinition(const String & mod)
        """
        ...
    
    @overload
    def __init__(self, mod: Union[bytes, str, String] , fixed: bool ) -> None:
        """
        Cython signature: void ModificationDefinition(const String & mod, bool fixed)
        """
        ...
    
    @overload
    def __init__(self, mod: Union[bytes, str, String] , fixed: bool , max_occur: int ) -> None:
        """
        Cython signature: void ModificationDefinition(const String & mod, bool fixed, unsigned int max_occur)
        """
        ...
    
    @overload
    def __init__(self, mod: ResidueModification ) -> None:
        """
        Cython signature: void ModificationDefinition(ResidueModification & mod)
        """
        ...
    
    @overload
    def __init__(self, mod: ResidueModification , fixed: bool ) -> None:
        """
        Cython signature: void ModificationDefinition(ResidueModification & mod, bool fixed)
        """
        ...
    
    @overload
    def __init__(self, mod: ResidueModification , fixed: bool , max_occur: int ) -> None:
        """
        Cython signature: void ModificationDefinition(ResidueModification & mod, bool fixed, unsigned int max_occur)
        """
        ...
    
    def setFixedModification(self, fixed: bool ) -> None:
        """
        Cython signature: void setFixedModification(bool fixed)
        Sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
        """
        ...
    
    def isFixedModification(self) -> bool:
        """
        Cython signature: bool isFixedModification()
        Returns if the modification if fixed true, else false
        """
        ...
    
    def setMaxOccurrences(self, num: int ) -> None:
        """
        Cython signature: void setMaxOccurrences(unsigned int num)
        Sets the maximal number of occurrences per peptide (unbounded if 0)
        """
        ...
    
    def getMaxOccurrences(self) -> int:
        """
        Cython signature: unsigned int getMaxOccurrences()
        Returns the maximal number of occurrences per peptide
        """
        ...
    
    def getModificationName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getModificationName()
        Returns the name of the modification
        """
        ...
    
    def setModification(self, modification: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setModification(const String & modification)
        Sets the modification, allowed are unique names provided by ModificationsDB
        """
        ...
    
    def getModification(self) -> ResidueModification:
        """
        Cython signature: ResidueModification getModification()
        """
        ...
    
    def __richcmp__(self, other: ModificationDefinition, op: int) -> Any:
        ... 


class MultiplexDeltaMassesGenerator:
    """
    Cython implementation of _MultiplexDeltaMassesGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMassesGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MultiplexDeltaMassesGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: MultiplexDeltaMassesGenerator ) -> None:
        """
        Cython signature: void MultiplexDeltaMassesGenerator(MultiplexDeltaMassesGenerator &)
        """
        ...
    
    @overload
    def __init__(self, labels: Union[bytes, str, String] , missed_cleavages: int , label_mass_shift: Dict[Union[bytes, str, String], float] ) -> None:
        """
        Cython signature: void MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, libcpp_map[String,double] label_mass_shift)
        """
        ...
    
    def generateKnockoutDeltaMasses(self) -> None:
        """
        Cython signature: void generateKnockoutDeltaMasses()
        """
        ...
    
    def getDeltaMassesList(self) -> List[MultiplexDeltaMasses]:
        """
        Cython signature: libcpp_vector[MultiplexDeltaMasses] getDeltaMassesList()
        """
        ...
    
    def getLabelShort(self, label: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String getLabelShort(String label)
        """
        ...
    
    def getLabelLong(self, label: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String getLabelLong(String label)
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


class MultiplexDeltaMassesGenerator_Label:
    """
    Cython implementation of _MultiplexDeltaMassesGenerator_Label

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMassesGenerator_Label.html>`_
    """
    
    short_name: Union[bytes, str, String]
    
    long_name: Union[bytes, str, String]
    
    description: Union[bytes, str, String]
    
    delta_mass: float
    
    def __init__(self, sn: Union[bytes, str, String] , ln: Union[bytes, str, String] , d: Union[bytes, str, String] , dm: float ) -> None:
        """
        Cython signature: void MultiplexDeltaMassesGenerator_Label(String sn, String ln, String d, double dm)
        """
        ... 


class MzQCFile:
    """
    Cython implementation of _MzQCFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzQCFile.html>`_

    File adapter for mzQC files used to load and store mzQC files
    
    This class collects the data for the mzQC File
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MzQCFile()
        """
        ...
    
    def store(self, input_file: Union[bytes, str, String] , output_file: Union[bytes, str, String] , exp: MSExperiment , contact_name: Union[bytes, str, String] , contact_address: Union[bytes, str, String] , description: Union[bytes, str, String] , label: Union[bytes, str, String] , feature_map: FeatureMap , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void store(String input_file, String output_file, MSExperiment & exp, String contact_name, String contact_address, String description, String label, FeatureMap & feature_map, libcpp_vector[ProteinIdentification] & prot_ids, PeptideIdentificationList & pep_ids)
        Stores QC data in mzQC file with JSON format
        
        
        :param input_file: MzML input file name
        :param output_file: MzQC output file name
        :param exp: MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
        :param contact_name: Name of the person creating the mzQC file
        :param contact_address: Contact address (mail/e-mail or phone) of the person creating the mzQC file
        :param description: Description and comments about the mzQC file contents
        :param label: Qnique and informative label for the run
        :param feature_map: FeatureMap from feature file (featureXML)
        :param prot_ids: Protein identifications from ID file (idXML)
        :param pep_ids: Protein identifications from ID file (idXML)
        """
        ... 


class PeptideEvidence:
    """
    Cython implementation of _PeptideEvidence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideEvidence.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideEvidence()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideEvidence ) -> None:
        """
        Cython signature: void PeptideEvidence(PeptideEvidence &)
        """
        ...
    
    def setStart(self, start: int ) -> None:
        """
        Cython signature: void setStart(int start)
        Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with `N_TERMINAL_AA`
        """
        ...
    
    def getStart(self) -> int:
        """
        Cython signature: int getStart()
        Returns the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        """
        ...
    
    def setEnd(self, end: int ) -> None:
        """
        Cython signature: void setEnd(int end)
        Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA
        """
        ...
    
    def getEnd(self) -> int:
        """
        Cython signature: int getEnd()
        Returns the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        """
        ...
    
    def setAABefore(self, rhs: bytes ) -> None:
        """
        Cython signature: void setAABefore(char rhs)
        Sets the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA
        """
        ...
    
    def getAABefore(self) -> bytes:
        """
        Cython signature: char getAABefore()
        Returns the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned
        """
        ...
    
    def setAAAfter(self, rhs: bytes ) -> None:
        """
        Cython signature: void setAAAfter(char rhs)
        Sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA
        """
        ...
    
    def getAAAfter(self) -> bytes:
        """
        Cython signature: char getAAAfter()
        Returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned
        """
        ...
    
    def setProteinAccession(self, s: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setProteinAccession(String s)
        Sets the protein accession the peptide matches to. If not available set to empty string
        """
        ...
    
    def getProteinAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getProteinAccession()
        Returns the protein accession the peptide matches to. If not available the empty string is returned
        """
        ...
    
    def hasValidLimits(self) -> bool:
        """
        Cython signature: bool hasValidLimits()
        Start and end numbers in evidence represent actual numeric indices
        """
        ...
    
    def __richcmp__(self, other: PeptideEvidence, op: int) -> Any:
        ... 


class PrecursorCorrection:
    """
    Cython implementation of _PrecursorCorrection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PrecursorCorrection.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PrecursorCorrection()
        """
        ...
    
    @overload
    def __init__(self, in_0: PrecursorCorrection ) -> None:
        """
        Cython signature: void PrecursorCorrection(PrecursorCorrection &)
        """
        ...
    
    correctToHighestIntensityMS1Peak: __static_PrecursorCorrection_correctToHighestIntensityMS1Peak
    
    correctToNearestFeature: __static_PrecursorCorrection_correctToNearestFeature
    
    correctToNearestMS1Peak: __static_PrecursorCorrection_correctToNearestMS1Peak
    
    getPrecursors: __static_PrecursorCorrection_getPrecursors
    
    writeHist: __static_PrecursorCorrection_writeHist 


class RibonucleotideDB:
    """
    Cython implementation of _RibonucleotideDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RibonucleotideDB.html>`_
    """
    
    def getRibonucleotide(self, code: bytes ) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * getRibonucleotide(const libcpp_string & code)
        """
        ...
    
    def getRibonucleotidePrefix(self, code: bytes ) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * getRibonucleotidePrefix(const libcpp_string & code)
        """
        ... 

