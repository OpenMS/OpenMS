from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_MetaboliteSpectralMatching_computeHyperScore(fragment_mass_error: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , db_spectrum: MSSpectrum , annotations: List[PeptideHit_PeakAnnotation] , mz_lower_bound: float ) -> float:
    """
    Cython signature: double computeHyperScore(double fragment_mass_error, bool fragment_mass_tolerance_unit_ppm, MSSpectrum exp_spectrum, MSSpectrum db_spectrum, libcpp_vector[PeptideHit_PeakAnnotation] & annotations, double mz_lower_bound)
    """
    ...


class CVMappingTerm:
    """
    Cython implementation of _CVMappingTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingTerm.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVMappingTerm()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVMappingTerm ) -> None:
        """
        Cython signature: void CVMappingTerm(CVMappingTerm &)
        """
        ...
    
    def setAccession(self, accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setAccession(String accession)
        Sets the accession string of the term
        """
        ...
    
    def getAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getAccession()
        Returns the accession string of the term
        """
        ...
    
    def setUseTermName(self, use_term_name: bool ) -> None:
        """
        Cython signature: void setUseTermName(bool use_term_name)
        Sets whether the term name should be used, instead of the accession
        """
        ...
    
    def getUseTermName(self) -> bool:
        """
        Cython signature: bool getUseTermName()
        Returns whether the term name should be used, instead of the accession
        """
        ...
    
    def setUseTerm(self, use_term: bool ) -> None:
        """
        Cython signature: void setUseTerm(bool use_term)
        Sets whether the term itself can be used (or only its children)
        """
        ...
    
    def getUseTerm(self) -> bool:
        """
        Cython signature: bool getUseTerm()
        Returns true if the term can be used, false if only children are allowed
        """
        ...
    
    def setTermName(self, term_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTermName(String term_name)
        Sets the name of the term
        """
        ...
    
    def getTermName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTermName()
        Returns the name of the term
        """
        ...
    
    def setIsRepeatable(self, is_repeatable: bool ) -> None:
        """
        Cython signature: void setIsRepeatable(bool is_repeatable)
        Sets whether this term can be repeated
        """
        ...
    
    def getIsRepeatable(self) -> bool:
        """
        Cython signature: bool getIsRepeatable()
        Returns true if this term can be repeated, false otherwise
        """
        ...
    
    def setAllowChildren(self, allow_children: bool ) -> None:
        """
        Cython signature: void setAllowChildren(bool allow_children)
        Sets whether children of this term are allowed
        """
        ...
    
    def getAllowChildren(self) -> bool:
        """
        Cython signature: bool getAllowChildren()
        Returns true if the children of this term are allowed to be used
        """
        ...
    
    def setCVIdentifierRef(self, cv_identifier_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCVIdentifierRef(String cv_identifier_ref)
        Sets the CV identifier reference string, e.g. UO for unit obo
        """
        ...
    
    def getCVIdentifierRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCVIdentifierRef()
        Returns the CV identifier reference string
        """
        ...
    
    def __richcmp__(self, other: CVMappingTerm, op: int) -> Any:
        ... 


class ConsensusIDAlgorithmAverage:
    """
    Cython implementation of _ConsensusIDAlgorithmAverage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmAverage.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmAverage()
        """
        ...
    
    def apply(self, ids: PeptideIdentificationList , number_of_runs: int ) -> None:
        """
        Cython signature: void apply(PeptideIdentificationList & ids, size_t number_of_runs)
        Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature
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


class ConsensusMapNormalizerAlgorithmMedian:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmMedian

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMapNormalizerAlgorithmMedian_1_1ConsensusMapNormalizerAlgorithmMedian.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusMapNormalizerAlgorithmMedian()
        """
        ...
    
    def computeMedians(self, input_map: ConsensusMap , medians: List[float] , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t computeMedians(ConsensusMap & input_map, libcpp_vector[double] & medians, const String & acc_filter, const String & desc_filter)
        Computes medians of all maps and returns index of map with most features
        """
        ...
    
    def normalizeMaps(self, input_map: ConsensusMap , method: int , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void normalizeMaps(ConsensusMap & input_map, NormalizationMethod method, const String & acc_filter, const String & desc_filter)
        Normalizes the maps of the consensusMap
        """
        ... 


class FeatureFinderAlgorithmMetaboIdent:
    """
    Cython implementation of _FeatureFinderAlgorithmMetaboIdent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderAlgorithmMetaboIdent.html>`_
      -- Inherits from ['DefaultParamHandler']

    Perform targeted feature extraction of compounds provided as table and stores them in features
    
    The algorithms detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications
    Internally, it uses algorithms for targeted data analysis from the OpenSWATH pipeline
    In the simplest case, only CompoundName, SumFormula, Charge and RetentionTime need to be given, all other values may be zero
    Every combination of compound (mass), RT and charge defines one target for feature detection
    Output:
    The main output is a feature map of detected features, with annotations in meta data entries
    Additional outputs are the extracted chromatograms/peak groups, the assay in TraML compatible format, and transformations
    that contain the error between provided and observed peaks
    
    Usage:
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureFinderAlgorithmMetaboIdent()
        """
        ...
    
    def setMSData(self, input: MSExperiment ) -> None:
        """
        Cython signature: void setMSData(MSExperiment & input)
        Sets spectra
        """
        ...
    
    def getMSData(self) -> MSExperiment:
        """
        Cython signature: const MSExperiment & getMSData()
        Returns spectra
        """
        ...
    
    def run(self, metaboIdentTable: List[FeatureFinderMetaboIdentCompound] , features: FeatureMap , spectra_path: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void run(const libcpp_vector[FeatureFinderMetaboIdentCompound] metaboIdentTable, FeatureMap & features, String spectra_path)
         Run feature extraction. spectra_path get's annotated as primaryMSRunPath in the resulting feature map.
        """
        ...
    
    def getChromatograms(self) -> MSExperiment:
        """
        Cython signature: MSExperiment & getChromatograms()
        Retrieves chromatograms (empty if run was not executed)
        """
        ...
    
    def getLibrary(self) -> TargetedExperiment:
        """
        Cython signature: const TargetedExperiment & getLibrary()
        Retrieves the assay library (e.g., to store as TraML, empty if run was not executed)
        """
        ...
    
    def getTransformations(self) -> TransformationDescription:
        """
        Cython signature: const TransformationDescription & getTransformations()
        Retrieves deviations between provided coordinates and extacted ones (e.g., to store as TrafoXML or for plotting)
        """
        ...
    
    def getNShared(self) -> int:
        """
        Cython signature: size_t getNShared()
        Retrieves number of features with shared identifications
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


class FeatureFinderMetaboIdentCompound:
    """
    Cython implementation of _FeatureFinderMetaboIdentCompound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderMetaboIdentCompound.html>`_
    """
    
    def __init__(self, name: Union[bytes, str, String] , formula: Union[bytes, str, String] , mass: float , charges: List[int] , rts: List[float] , rt_ranges: List[float] , iso_distrib: List[float] ) -> None:
        """
        Cython signature: void FeatureFinderMetaboIdentCompound(String name, String formula, double mass, libcpp_vector[int] charges, libcpp_vector[double] rts, libcpp_vector[double] rt_ranges, libcpp_vector[double] iso_distrib)
          Represents a compound in the in the FeatureFinderMetaboIdent library table.
        
        
          :param name: Unique name for the target compound.
          :param formula: Chemical sum formula.
          :param mass: Neutral mass; if zero calculated from formula.
          :param charges: List of possible charge states.
          :param rts: List of possible retention times.
          :param rt_ranges: List of possible retention time ranges (window around RT), either one value or one per RT entry.
          :param iso_distrib: List of relative abundances of isotopologues; if zero calculated from formula.
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
          Gets the compound name.
        
        
          :rtype: str
        """
        ...
    
    def getFormula(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFormula()
          Gets the compound chemical formula.
        
        
          :rtype: str
        """
        ...
    
    def getMass(self) -> float:
        """
        Cython signature: double getMass()
          Gets the compound mass.
        
        
          :rtype: float
        """
        ...
    
    def getCharges(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getCharges()
          Gets the compound charge states.
        
        
          :rtype: list of int
        """
        ...
    
    def getRTs(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getRTs()
          Gets the compound retention times.
        
        
          :rtype: list of float
        """
        ...
    
    def getRTRanges(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getRTRanges()
          Gets the compound retention time ranges.
        
        
          :rtype: list of float
        """
        ...
    
    def getIsotopeDistribution(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getIsotopeDistribution()
          Gets the compound isotopic distributions.
        
        
          :rtype: list of float
        """
        ... 


class FeatureFinderMultiplexAlgorithm:
    """
    Cython implementation of _FeatureFinderMultiplexAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderMultiplexAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureFinderMultiplexAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureFinderMultiplexAlgorithm ) -> None:
        """
        Cython signature: void FeatureFinderMultiplexAlgorithm(FeatureFinderMultiplexAlgorithm &)
        """
        ...
    
    def run(self, exp: MSExperiment , progress: bool ) -> None:
        """
        Cython signature: void run(MSExperiment & exp, bool progress)
        Main method for feature detection
        """
        ...
    
    def getFeatureMap(self) -> FeatureMap:
        """
        Cython signature: FeatureMap getFeatureMap()
        """
        ...
    
    def getConsensusMap(self) -> ConsensusMap:
        """
        Cython signature: ConsensusMap getConsensusMap()
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


class FeatureMap:
    """
    Cython implementation of _FeatureMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMap.html>`_
      -- Inherits from ['UniqueIdInterface', 'DocumentIdentifier', 'RangeManagerRtMzInt', 'MetaInfoInterface']

    A container for features.
    
    A feature map is a container holding features, which represent
    chemical entities (peptides, proteins, small molecules etc.) found
    in an LC-MS/MS experiment.
    
    This class supports direct iteration in Python.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureMap()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureMap ) -> None:
        """
        Cython signature: void FeatureMap(FeatureMap &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def __getitem__(self, in_0: int ) -> Feature:
        """
        Cython signature: Feature & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: Feature ) -> None:
        """Cython signature: Feature & operator[](size_t)"""
        ...
    
    @overload
    def push_back(self, spec: Feature ) -> None:
        """
        Cython signature: void push_back(Feature spec)
        """
        ...
    
    @overload
    def push_back(self, spec: MRMFeature ) -> None:
        """
        Cython signature: void push_back(MRMFeature spec)
        """
        ...
    
    @overload
    def sortByIntensity(self, ) -> None:
        """
        Cython signature: void sortByIntensity()
        Sorts the peaks according to ascending intensity
        """
        ...
    
    @overload
    def sortByIntensity(self, reverse: bool ) -> None:
        """
        Cython signature: void sortByIntensity(bool reverse)
        Sorts the peaks according to ascending intensity. Order is reversed if argument is `true` ( reverse = true )
        """
        ...
    
    def sortByPosition(self) -> None:
        """
        Cython signature: void sortByPosition()
        Sorts features by position. Lexicographical comparison (first RT then m/z) is done
        """
        ...
    
    def sortByRT(self) -> None:
        """
        Cython signature: void sortByRT()
        Sorts features by RT position
        """
        ...
    
    def sortByMZ(self) -> None:
        """
        Cython signature: void sortByMZ()
        Sorts features by m/z position
        """
        ...
    
    def sortByOverallQuality(self) -> None:
        """
        Cython signature: void sortByOverallQuality()
        Sorts features by ascending overall quality. Order is reversed if argument is `true` ( reverse = true )
        """
        ...
    
    def swap(self, in_0: FeatureMap ) -> None:
        """
        Cython signature: void swap(FeatureMap &)
        """
        ...
    
    def swapFeaturesOnly(self, swapfrom: FeatureMap ) -> None:
        """
        Cython signature: void swapFeaturesOnly(FeatureMap swapfrom)
        Swaps the feature content (plus its range information) of this map
        """
        ...
    
    @overload
    def clear(self, ) -> None:
        """
        Cython signature: void clear()
        Clears all data and meta data
        """
        ...
    
    @overload
    def clear(self, clear_meta_data: bool ) -> None:
        """
        Cython signature: void clear(bool clear_meta_data)
        Clears all data and meta data. If 'true' is passed as an argument, all meta data is cleared in addition to the data
        """
        ...
    
    def __add__(self: FeatureMap, other: FeatureMap) -> FeatureMap:
        ...
    
    def __iadd__(self: FeatureMap, other: FeatureMap) -> FeatureMap:
        ...
    
    def updateRanges(self) -> None:
        """
        Cython signature: void updateRanges()
        """
        ...
    
    def getProteinIdentifications(self) -> List[ProteinIdentification]:
        """
        Cython signature: libcpp_vector[ProteinIdentification] getProteinIdentifications()
        """
        ...
    
    def setProteinIdentifications(self, in_0: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void setProteinIdentifications(libcpp_vector[ProteinIdentification])
        Sets the protein identifications
        """
        ...
    
    def getUnassignedPeptideIdentifications(self) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList getUnassignedPeptideIdentifications()
        """
        ...
    
    def setUnassignedPeptideIdentifications(self, in_0: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setUnassignedPeptideIdentifications(PeptideIdentificationList)
        Sets the unassigned peptide identifications
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[DataProcessing] getDataProcessing()
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[DataProcessing])
        Sets the description of the applied data processing
        """
        ...
    
    @overload
    def setPrimaryMSRunPath(self, s: List[bytes] ) -> None:
        """
        Cython signature: void setPrimaryMSRunPath(StringList & s)
        Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
        """
        ...
    
    @overload
    def setPrimaryMSRunPath(self, s: List[bytes] , e: MSExperiment ) -> None:
        """
        Cython signature: void setPrimaryMSRunPath(StringList & s, MSExperiment & e)
        Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
        """
        ...
    
    def getPrimaryMSRunPath(self, toFill: List[bytes] ) -> None:
        """
        Cython signature: void getPrimaryMSRunPath(StringList & toFill)
        Returns the file path to the first MS run
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
    
    def setIdentifier(self, id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(String id)
        Sets document identifier (e.g. an LSID)
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        Retrieve document identifier (e.g. an LSID)
        """
        ...
    
    def setLoadedFileType(self, file_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLoadedFileType(String file_name)
        Sets the file_type according to the type of the file loaded from, preferably done whilst loading
        """
        ...
    
    def getLoadedFileType(self) -> int:
        """
        Cython signature: int getLoadedFileType()
        Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded
        """
        ...
    
    def setLoadedFilePath(self, file_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLoadedFilePath(String file_name)
        Sets the file_name according to absolute path of the file loaded, preferably done whilst loading
        """
        ...
    
    def getLoadedFilePath(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getLoadedFilePath()
        Returns the file_name which is the absolute path to the file loaded
        """
        ...
    
    def getMinRT(self) -> float:
        """
        Cython signature: double getMinRT()
        Returns the minimum RT
        """
        ...
    
    def getMaxRT(self) -> float:
        """
        Cython signature: double getMaxRT()
        Returns the maximum RT
        """
        ...
    
    def getMinMZ(self) -> float:
        """
        Cython signature: double getMinMZ()
        Returns the minimum m/z
        """
        ...
    
    def getMaxMZ(self) -> float:
        """
        Cython signature: double getMaxMZ()
        Returns the maximum m/z
        """
        ...
    
    def getMinIntensity(self) -> float:
        """
        Cython signature: double getMinIntensity()
        Returns the minimum intensity
        """
        ...
    
    def getMaxIntensity(self) -> float:
        """
        Cython signature: double getMaxIntensity()
        Returns the maximum intensity
        """
        ...
    
    def clearRanges(self) -> None:
        """
        Cython signature: void clearRanges()
        Resets all range dimensions as empty
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
    
    def __richcmp__(self, other: FeatureMap, op: int) -> Any:
        ...
    
    def __iter__(self) -> Feature:
       ... 


class GaussFilter:
    """
    Cython implementation of _GaussFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GaussFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GaussFilter()
        This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data
        """
        ...
    
    @overload
    def __init__(self, in_0: GaussFilter ) -> None:
        """
        Cython signature: void GaussFilter(GaussFilter &)
        """
        ...
    
    @overload
    def filter(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void filter(MSSpectrum & spectrum)
        Smoothes an MSSpectrum containing profile data
        """
        ...
    
    @overload
    def filter(self, chromatogram: MSChromatogram ) -> None:
        """
        Cython signature: void filter(MSChromatogram & chromatogram)
        """
        ...
    
    def filterExperiment(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterExperiment(MSExperiment & exp)
        Smoothes an MSExperiment containing profile data
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


class IndexedMzMLHandler:
    """
    Cython implementation of _IndexedMzMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLHandler.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IndexedMzMLHandler()
        """
        ...
    
    @overload
    def __init__(self, in_0: IndexedMzMLHandler ) -> None:
        """
        Cython signature: void IndexedMzMLHandler(IndexedMzMLHandler &)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void IndexedMzMLHandler(String filename)
        """
        ...
    
    def openFile(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void openFile(String filename)
        """
        ...
    
    def getParsingSuccess(self) -> bool:
        """
        Cython signature: bool getParsingSuccess()
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        """
        ...
    
    def getSpectrumById(self, id_: int ) -> _Interfaces_Spectrum:
        """
        Cython signature: shared_ptr[_Interfaces_Spectrum] getSpectrumById(int id_)
        """
        ...
    
    def getChromatogramById(self, id_: int ) -> _Interfaces_Chromatogram:
        """
        Cython signature: shared_ptr[_Interfaces_Chromatogram] getChromatogramById(int id_)
        """
        ...
    
    def getMSSpectrumById(self, id_: int ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum getMSSpectrumById(int id_)
        """
        ...
    
    def getMSSpectrumByNativeId(self, id_: bytes , spec: MSSpectrum ) -> None:
        """
        Cython signature: void getMSSpectrumByNativeId(libcpp_string id_, MSSpectrum & spec)
        """
        ...
    
    def getMSChromatogramById(self, id_: int ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getMSChromatogramById(int id_)
        """
        ...
    
    def getMSChromatogramByNativeId(self, id_: bytes , chrom: MSChromatogram ) -> None:
        """
        Cython signature: void getMSChromatogramByNativeId(libcpp_string id_, MSChromatogram & chrom)
        """
        ...
    
    def setSkipXMLChecks(self, skip: bool ) -> None:
        """
        Cython signature: void setSkipXMLChecks(bool skip)
        """
        ... 


class MapConversion:
    """
    Cython implementation of _MapConversion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapConversion.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MapConversion()
        """
        ...
    
    @overload
    def __init__(self, in_0: MapConversion ) -> None:
        """
        Cython signature: void MapConversion(MapConversion &)
        """
        ...
    
    @overload
    def convert(self, input_map_index: int , input_map: FeatureMap , output_map: ConsensusMap , n: int ) -> None:
        """
        Cython signature: void convert(uint64_t input_map_index, FeatureMap input_map, ConsensusMap & output_map, size_t n)
        """
        ...
    
    @overload
    def convert(self, input_map_index: int , input_map: MSExperiment , output_map: ConsensusMap , n: int ) -> None:
        """
        Cython signature: void convert(uint64_t input_map_index, MSExperiment & input_map, ConsensusMap & output_map, size_t n)
        """
        ...
    
    @overload
    def convert(self, input_map: ConsensusMap , keep_uids: bool , output_map: FeatureMap ) -> None:
        """
        Cython signature: void convert(ConsensusMap input_map, bool keep_uids, FeatureMap & output_map)
        """
        ... 


class MetaboliteSpectralMatching:
    """
    Cython implementation of _MetaboliteSpectralMatching

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboliteSpectralMatching.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboliteSpectralMatching()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboliteSpectralMatching ) -> None:
        """
        Cython signature: void MetaboliteSpectralMatching(MetaboliteSpectralMatching &)
        """
        ...
    
    def run(self, exp: MSExperiment , speclib: MSExperiment , mz_tab: MzTab , out_spectra: String ) -> None:
        """
        Cython signature: void run(MSExperiment & exp, MSExperiment & speclib, MzTab & mz_tab, String & out_spectra)
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
    
    computeHyperScore: __static_MetaboliteSpectralMatching_computeHyperScore 


class ModificationDefinitionsSet:
    """
    Cython implementation of _ModificationDefinitionsSet

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationDefinitionsSet.html>`_

    Representation of a set of modification definitions
    
    This class enhances the modification definitions as defined in the
    class ModificationDefinition into a set of definitions. This is also
    e.g. used as input parameters in search engines.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ModificationDefinitionsSet()
        """
        ...
    
    @overload
    def __init__(self, in_0: ModificationDefinitionsSet ) -> None:
        """
        Cython signature: void ModificationDefinitionsSet(ModificationDefinitionsSet &)
        """
        ...
    
    @overload
    def __init__(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None:
        """
        Cython signature: void ModificationDefinitionsSet(StringList fixed_modifications, StringList variable_modifications)
        """
        ...
    
    def setMaxModifications(self, max_mod: int ) -> None:
        """
        Cython signature: void setMaxModifications(size_t max_mod)
        Sets the maximal number of modifications allowed per peptide
        """
        ...
    
    def getMaxModifications(self) -> int:
        """
        Cython signature: size_t getMaxModifications()
        Return the maximal number of modifications allowed per peptide
        """
        ...
    
    def getNumberOfModifications(self) -> int:
        """
        Cython signature: size_t getNumberOfModifications()
        Returns the number of modifications stored in this set
        """
        ...
    
    def getNumberOfFixedModifications(self) -> int:
        """
        Cython signature: size_t getNumberOfFixedModifications()
        Returns the number of fixed modifications stored in this set
        """
        ...
    
    def getNumberOfVariableModifications(self) -> int:
        """
        Cython signature: size_t getNumberOfVariableModifications()
        Returns the number of variable modifications stored in this set
        """
        ...
    
    def addModification(self, mod_def: ModificationDefinition ) -> None:
        """
        Cython signature: void addModification(ModificationDefinition & mod_def)
        Adds a modification definition to the set
        """
        ...
    
    @overload
    def setModifications(self, mod_defs: Set[ModificationDefinition] ) -> None:
        """
        Cython signature: void setModifications(libcpp_set[ModificationDefinition] & mod_defs)
        Sets the modification definitions
        """
        ...
    
    @overload
    def setModifications(self, fixed_modifications: Union[bytes, str, String] , variable_modifications: String ) -> None:
        """
        Cython signature: void setModifications(const String & fixed_modifications, String & variable_modifications)
        Set the modification definitions from a string
        
        The strings should contain a comma separated list of modifications. The names
        can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
        definitions and other specific definitions are given by the modifications themselves.
        """
        ...
    
    @overload
    def setModifications(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None:
        """
        Cython signature: void setModifications(StringList & fixed_modifications, StringList & variable_modifications)
        Same as above, but using StringList instead of comma separated strings
        """
        ...
    
    def getModifications(self) -> Set[ModificationDefinition]:
        """
        Cython signature: libcpp_set[ModificationDefinition] getModifications()
        Returns the stored modification definitions
        """
        ...
    
    def getFixedModifications(self) -> Set[ModificationDefinition]:
        """
        Cython signature: libcpp_set[ModificationDefinition] getFixedModifications()
        Returns the stored fixed modification definitions
        """
        ...
    
    def getVariableModifications(self) -> Set[ModificationDefinition]:
        """
        Cython signature: libcpp_set[ModificationDefinition] getVariableModifications()
        Returns the stored variable modification definitions
        """
        ...
    
    @overload
    def getModificationNames(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None:
        """
        Cython signature: void getModificationNames(StringList & fixed_modifications, StringList & variable_modifications)
        Populates the output lists with the modification names (use e.g. for
        """
        ...
    
    @overload
    def getModificationNames(self, ) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getModificationNames()
        Returns only the names of the modifications stored in the set
        """
        ...
    
    def getFixedModificationNames(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getFixedModificationNames()
        Returns only the names of the fixed modifications
        """
        ...
    
    def getVariableModificationNames(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getVariableModificationNames()
        Returns only the names of the variable modifications
        """
        ...
    
    def isCompatible(self, peptide: AASequence ) -> bool:
        """
        Cython signature: bool isCompatible(AASequence & peptide)
        Returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications
        """
        ...
    
    def inferFromPeptides(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void inferFromPeptides(PeptideIdentificationList peptides)
        Infers the sets of defined modifications from the modifications present on peptide identifications
        """
        ... 


class NoiseEstimator:
    """
    Cython implementation of _NoiseEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NoiseEstimator.html>`_
    """
    
    nr_windows: int
    
    mz_start: float
    
    window_length: float
    
    result_windows_even: List[float]
    
    result_windows_odd: List[float]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NoiseEstimator()
        """
        ...
    
    @overload
    def __init__(self, in_0: NoiseEstimator ) -> None:
        """
        Cython signature: void NoiseEstimator(NoiseEstimator &)
        """
        ...
    
    @overload
    def __init__(self, nr_windows_: float , mz_start_: float , win_len_: float ) -> None:
        """
        Cython signature: void NoiseEstimator(double nr_windows_, double mz_start_, double win_len_)
        """
        ...
    
    def get_noise_value(self, mz: float ) -> float:
        """
        Cython signature: double get_noise_value(double mz)
        """
        ...
    
    def get_noise_even(self, mz: float ) -> float:
        """
        Cython signature: double get_noise_even(double mz)
        """
        ...
    
    def get_noise_odd(self, mz: float ) -> float:
        """
        Cython signature: double get_noise_odd(double mz)
        """
        ... 


class SignalToNoiseEstimatorMedianRapid:
    """
    Cython implementation of _SignalToNoiseEstimatorMedianRapid

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMedianRapid.html>`_
    """
    
    @overload
    def __init__(self, in_0: SignalToNoiseEstimatorMedianRapid ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMedianRapid(SignalToNoiseEstimatorMedianRapid &)
        """
        ...
    
    @overload
    def __init__(self, window_length: float ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMedianRapid(double window_length)
        """
        ...
    
    @overload
    def estimateNoise(self, in_0: _Interfaces_Spectrum ) -> NoiseEstimator:
        """
        Cython signature: NoiseEstimator estimateNoise(shared_ptr[_Interfaces_Spectrum])
        """
        ...
    
    @overload
    def estimateNoise(self, in_0: _Interfaces_Chromatogram ) -> NoiseEstimator:
        """
        Cython signature: NoiseEstimator estimateNoise(shared_ptr[_Interfaces_Chromatogram])
        """
        ...
    
    @overload
    def estimateNoise(self, mz_array: List[float] , int_array: List[float] ) -> NoiseEstimator:
        """
        Cython signature: NoiseEstimator estimateNoise(libcpp_vector[double] mz_array, libcpp_vector[double] int_array)
        """
        ... 


class SpectraSTSimilarityScore:
    """
    Cython implementation of _SpectraSTSimilarityScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectraSTSimilarityScore.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectraSTSimilarityScore()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectraSTSimilarityScore ) -> None:
        """
        Cython signature: void SpectraSTSimilarityScore(SpectraSTSimilarityScore &)
        """
        ...
    
    def preprocess(self, spec: MSSpectrum , remove_peak_intensity_threshold: float , cut_peaks_below: int , min_peak_number: int , max_peak_number: int ) -> bool:
        """
        Cython signature: bool preprocess(MSSpectrum & spec, float remove_peak_intensity_threshold, unsigned int cut_peaks_below, size_t min_peak_number, size_t max_peak_number)
        Preprocesses the spectrum
        
        The preprocessing removes peak below a intensity threshold, reject spectra that does
        not have enough peaks, and cuts peaks exceeding the max_peak_number most intense peaks
        
        :returns: true if spectrum passes filtering
        """
        ...
    
    def transform(self, spec: MSSpectrum ) -> BinnedSpectrum:
        """
        Cython signature: BinnedSpectrum transform(MSSpectrum & spec)
        Spectrum is transformed into a binned spectrum with bin size 1 and spread 1 and the intensities are normalized
        """
        ...
    
    def dot_bias(self, bin1: BinnedSpectrum , bin2: BinnedSpectrum , dot_product: float ) -> float:
        """
        Cython signature: double dot_bias(BinnedSpectrum & bin1, BinnedSpectrum & bin2, double dot_product)
        Calculates how much of the dot product is dominated by a few peaks
        
        :param dot_product: If -1 this value will be calculated as well.
        :param bin1: First spectrum in binned representation
        :param bin2: Second spectrum in binned representation
        """
        ...
    
    def delta_D(self, top_hit: float , runner_up: float ) -> float:
        """
        Cython signature: double delta_D(double top_hit, double runner_up)
        Calculates the normalized distance between top_hit and runner_up
        
        :param top_hit: Is the best score for a given match
        :param runner_up: A match with a worse score than top_hit, e.g. the second best score
        :returns: normalized distance
        """
        ...
    
    def compute_F(self, dot_product: float , delta_D: float , dot_bias: float ) -> float:
        """
        Cython signature: double compute_F(double dot_product, double delta_D, double dot_bias)
        Computes the overall all score
        
        :param dot_product: dot_product of a match
        :param delta_D: delta_D should be calculated after all dot products for a unidentified spectrum are computed
        :param dot_bias: the bias
        :returns: The SpectraST similarity score
        """
        ... 


class SpectralMatch:
    """
    Cython implementation of _SpectralMatch

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectralMatch.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectralMatch()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectralMatch ) -> None:
        """
        Cython signature: void SpectralMatch(SpectralMatch &)
        """
        ...
    
    def getObservedPrecursorMass(self) -> float:
        """
        Cython signature: double getObservedPrecursorMass()
        """
        ...
    
    def setObservedPrecursorMass(self, in_0: float ) -> None:
        """
        Cython signature: void setObservedPrecursorMass(double)
        """
        ...
    
    def getObservedPrecursorRT(self) -> float:
        """
        Cython signature: double getObservedPrecursorRT()
        """
        ...
    
    def setObservedPrecursorRT(self, in_0: float ) -> None:
        """
        Cython signature: void setObservedPrecursorRT(double)
        """
        ...
    
    def getFoundPrecursorMass(self) -> float:
        """
        Cython signature: double getFoundPrecursorMass()
        """
        ...
    
    def setFoundPrecursorMass(self, in_0: float ) -> None:
        """
        Cython signature: void setFoundPrecursorMass(double)
        """
        ...
    
    def getFoundPrecursorCharge(self) -> int:
        """
        Cython signature: int getFoundPrecursorCharge()
        """
        ...
    
    def setFoundPrecursorCharge(self, in_0: int ) -> None:
        """
        Cython signature: void setFoundPrecursorCharge(int)
        """
        ...
    
    def getMatchingScore(self) -> float:
        """
        Cython signature: double getMatchingScore()
        """
        ...
    
    def setMatchingScore(self, in_0: float ) -> None:
        """
        Cython signature: void setMatchingScore(double)
        """
        ...
    
    def getObservedSpectrumIndex(self) -> int:
        """
        Cython signature: size_t getObservedSpectrumIndex()
        """
        ...
    
    def setObservedSpectrumIndex(self, in_0: int ) -> None:
        """
        Cython signature: void setObservedSpectrumIndex(size_t)
        """
        ...
    
    def getMatchingSpectrumIndex(self) -> int:
        """
        Cython signature: size_t getMatchingSpectrumIndex()
        """
        ...
    
    def setMatchingSpectrumIndex(self, in_0: int ) -> None:
        """
        Cython signature: void setMatchingSpectrumIndex(size_t)
        """
        ...
    
    def getPrimaryIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPrimaryIdentifier()
        """
        ...
    
    def setPrimaryIdentifier(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPrimaryIdentifier(String)
        """
        ...
    
    def getSecondaryIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSecondaryIdentifier()
        """
        ...
    
    def setSecondaryIdentifier(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSecondaryIdentifier(String)
        """
        ...
    
    def getCommonName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCommonName()
        """
        ...
    
    def setCommonName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCommonName(String)
        """
        ...
    
    def getSumFormula(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSumFormula()
        """
        ...
    
    def setSumFormula(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSumFormula(String)
        """
        ...
    
    def getInchiString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getInchiString()
        """
        ...
    
    def setInchiString(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setInchiString(String)
        """
        ...
    
    def getSMILESString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSMILESString()
        """
        ...
    
    def setSMILESString(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSMILESString(String)
        """
        ...
    
    def getPrecursorAdduct(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPrecursorAdduct()
        """
        ...
    
    def setPrecursorAdduct(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPrecursorAdduct(String)
        """
        ... 


class SqMassConfig:
    """
    Cython implementation of _SqMassConfig

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqMassConfig.html>`_
    """
    
    write_full_meta: bool
    
    use_lossy_numpress: bool
    
    linear_fp_mass_acc: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SqMassConfig()
        """
        ...
    
    @overload
    def __init__(self, in_0: SqMassConfig ) -> None:
        """
        Cython signature: void SqMassConfig(SqMassConfig &)
        """
        ... 


class SqMassFile:
    """
    Cython implementation of _SqMassFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqMassFile.html>`_

    An class that uses on-disk SQLite database to read and write spectra and chromatograms
    
    This class provides functions to read and write spectra and chromatograms
    to disk using a SQLite database and store them in sqMass format. This
    allows users to access, select and filter spectra and chromatograms
    on-demand even in a large collection of data
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SqMassFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: SqMassFile ) -> None:
        """
        Cython signature: void SqMassFile(SqMassFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , map_: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & map_)
        Read / Write a complete mass spectrometric experiment
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , map_: MSExperiment ) -> None:
        """
        Cython signature: void store(const String & filename, MSExperiment & map_)
        Store an MSExperiment in sqMass format
        """
        ...
    
    def setConfig(self, config: SqMassConfig ) -> None:
        """
        Cython signature: void setConfig(SqMassConfig config)
        """
        ... 


class NormalizationMethod:
    None
    NM_SCALE : int
    NM_SHIFT : int

    def getMapping(self) -> Dict[int, str]:
       ... 

