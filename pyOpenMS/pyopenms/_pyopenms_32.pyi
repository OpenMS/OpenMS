from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AQS_featureConcentration:
    """
    Cython implementation of _AQS_featureConcentration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AQS_featureConcentration.html>`_
    """
    
    feature: Feature
    
    IS_feature: Feature
    
    actual_concentration: float
    
    IS_actual_concentration: float
    
    concentration_units: Union[bytes, str, String]
    
    dilution_factor: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AQS_featureConcentration()
        """
        ...
    
    @overload
    def __init__(self, in_0: AQS_featureConcentration ) -> None:
        """
        Cython signature: void AQS_featureConcentration(AQS_featureConcentration &)
        """
        ... 


class AQS_runConcentration:
    """
    Cython implementation of _AQS_runConcentration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AQS_runConcentration.html>`_
    """
    
    sample_name: Union[bytes, str, String]
    
    component_name: Union[bytes, str, String]
    
    IS_component_name: Union[bytes, str, String]
    
    actual_concentration: float
    
    IS_actual_concentration: float
    
    concentration_units: Union[bytes, str, String]
    
    dilution_factor: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AQS_runConcentration()
        """
        ...
    
    @overload
    def __init__(self, in_0: AQS_runConcentration ) -> None:
        """
        Cython signature: void AQS_runConcentration(AQS_runConcentration &)
        """
        ... 


class AbsoluteQuantitationStandards:
    """
    Cython implementation of _AbsoluteQuantitationStandards

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationStandards.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AbsoluteQuantitationStandards()
        """
        ...
    
    @overload
    def __init__(self, in_0: AbsoluteQuantitationStandards ) -> None:
        """
        Cython signature: void AbsoluteQuantitationStandards(AbsoluteQuantitationStandards &)
        """
        ...
    
    def getComponentFeatureConcentrations(self, run_concentrations: List[AQS_runConcentration] , feature_maps: List[FeatureMap] , component_name: Union[bytes, str, String] , feature_concentrations: List[AQS_featureConcentration] ) -> None:
        """
        Cython signature: void getComponentFeatureConcentrations(libcpp_vector[AQS_runConcentration] & run_concentrations, libcpp_vector[FeatureMap] & feature_maps, const String & component_name, libcpp_vector[AQS_featureConcentration] & feature_concentrations)
        """
        ... 


class ChromatogramRangeManager:
    """
    Cython implementation of _ChromatogramRangeManager

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramRangeManager.html>`_

    Range manager for chromatograms
    
    This class manages retention time, m/z, and intensity ranges for multiple chromatograms.
    It extends the basic RangeManager to provide specialized functionality for chromatogram data.
    
    The template parameters for the base RangeManager are ordered differently than in SpectrumRangeManager:
    - RangeRT (retention time) is the first parameter, as it's the primary dimension for chromatograms
    - RangeIntensity is the second parameter
    - RangeMZ is the third parameter
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramRangeManager()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramRangeManager ) -> None:
        """
        Cython signature: void ChromatogramRangeManager(ChromatogramRangeManager &)
        """
        ...
    
    def clearRanges(self) -> None:
        """
        Cython signature: void clearRanges()
        """
        ...
    
    def getMinRT(self) -> float:
        """
        Cython signature: double getMinRT()
        """
        ...
    
    def getMaxRT(self) -> float:
        """
        Cython signature: double getMaxRT()
        """
        ...
    
    def getMinMZ(self) -> float:
        """
        Cython signature: double getMinMZ()
        """
        ...
    
    def getMaxMZ(self) -> float:
        """
        Cython signature: double getMaxMZ()
        """
        ...
    
    def getMinIntensity(self) -> float:
        """
        Cython signature: double getMinIntensity()
        """
        ...
    
    def getMaxIntensity(self) -> float:
        """
        Cython signature: double getMaxIntensity()
        """
        ... 


class ColumnHeader:
    """
    Cython implementation of _ColumnHeader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMap_1_1ColumnHeader.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    filename: Union[bytes, str, String]
    
    label: Union[bytes, str, String]
    
    size: int
    
    unique_id: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ColumnHeader()
        """
        ...
    
    @overload
    def __init__(self, in_0: ColumnHeader ) -> None:
        """
        Cython signature: void ColumnHeader(ColumnHeader &)
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
    
    def __richcmp__(self, other: ColumnHeader, op: int) -> Any:
        ... 


class ConsensusMap:
    """
    Cython implementation of _ConsensusMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMap_1_1ConsensusMap.html>`_
      -- Inherits from ['UniqueIdInterface', 'DocumentIdentifier', 'RangeManagerRtMzInt', 'MetaInfoInterface']

    A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements
    (ConsensusFeature) which in turn represent analytes that have been
    quantified across multiple LC-MS/MS experiments. Each analyte in a
    ConsensusFeature is linked to its original LC-MS/MS run, the links are
    maintained by the ConsensusMap class.
    The map is implemented as a vector of elements of type ConsensusFeature.
    
    To be consistent, all maps who are referenced by ConsensusFeature objects
    (through a unique id) need to be registered in this class.
    
    This class supports direct iteration in Python.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ConsensusMap()
        """
        ...
    
    @overload
    def __init__(self, in_0: ConsensusMap ) -> None:
        """
        Cython signature: void ConsensusMap(ConsensusMap &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def reserve(self, s: int ) -> None:
        """
        Cython signature: void reserve(size_t s)
        """
        ...
    
    def __getitem__(self, in_0: int ) -> ConsensusFeature:
        """
        Cython signature: ConsensusFeature & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: ConsensusFeature ) -> None:
        """Cython signature: ConsensusFeature & operator[](size_t)"""
        ...
    
    def push_back(self, spec: ConsensusFeature ) -> None:
        """
        Cython signature: void push_back(ConsensusFeature spec)
        """
        ...
    
    def appendRows(self, in_0: ConsensusMap ) -> ConsensusMap:
        """
        Cython signature: ConsensusMap appendRows(ConsensusMap)
        Add consensus map entries as new rows
        """
        ...
    
    def appendColumns(self, in_0: ConsensusMap ) -> ConsensusMap:
        """
        Cython signature: ConsensusMap appendColumns(ConsensusMap)
        Add consensus map entries as new columns
        """
        ...
    
    @overload
    def clear(self, clear_meta_data: bool ) -> None:
        """
        Cython signature: void clear(bool clear_meta_data)
        Clears all data and meta data
        """
        ...
    
    @overload
    def clear(self, ) -> None:
        """
        Cython signature: void clear()
        """
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
    
    def setUnassignedPeptideIdentifications(self, unassigned_peptide_identifications: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setUnassignedPeptideIdentifications(PeptideIdentificationList unassigned_peptide_identifications)
        Sets the unassigned PeptideIdentificationList
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[DataProcessing] getDataProcessing()
        Returns a const reference to the description of the applied data processing
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
        Sets the file paths to the primary MS run (stored in ColumnHeaders)
        """
        ...
    
    @overload
    def setPrimaryMSRunPath(self, s: List[bytes] , e: MSExperiment ) -> None:
        """
        Cython signature: void setPrimaryMSRunPath(StringList & s, MSExperiment & e)
        """
        ...
    
    def getPrimaryMSRunPath(self, toFill: List[bytes] ) -> None:
        """
        Cython signature: void getPrimaryMSRunPath(StringList & toFill)
        Returns the MS run path (stored in ColumnHeaders)
        """
        ...
    
    @overload
    def sortByIntensity(self, reverse: bool ) -> None:
        """
        Cython signature: void sortByIntensity(bool reverse)
        Sorts the peaks according to ascending intensity.
        """
        ...
    
    @overload
    def sortByIntensity(self, ) -> None:
        """
        Cython signature: void sortByIntensity()
        """
        ...
    
    def sortByRT(self) -> None:
        """
        Cython signature: void sortByRT()
        Sorts the peaks according to RT position
        """
        ...
    
    def sortByMZ(self) -> None:
        """
        Cython signature: void sortByMZ()
        Sorts the peaks according to m/z position
        """
        ...
    
    def sortByPosition(self) -> None:
        """
        Cython signature: void sortByPosition()
        Lexicographically sorts the peaks by their position (First RT then m/z)
        """
        ...
    
    @overload
    def sortByQuality(self, reverse: bool ) -> None:
        """
        Cython signature: void sortByQuality(bool reverse)
        Sorts the peaks according to ascending quality.
        """
        ...
    
    @overload
    def sortByQuality(self, ) -> None:
        """
        Cython signature: void sortByQuality()
        """
        ...
    
    def sortBySize(self) -> None:
        """
        Cython signature: void sortBySize()
        Sorts with respect to the size (number of elements)
        """
        ...
    
    def sortByMaps(self) -> None:
        """
        Cython signature: void sortByMaps()
        Sorts with respect to the sets of maps covered by the consensus features (lexicographically)
        """
        ...
    
    def getExperimentType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getExperimentType()
        Non-mutable access to the experiment type
        """
        ...
    
    def setExperimentType(self, experiment_type: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setExperimentType(String experiment_type)
        Mutable access to the experiment type
        """
        ...
    
    def sortPeptideIdentificationsByMapIndex(self) -> None:
        """
        Cython signature: void sortPeptideIdentificationsByMapIndex()
        Sorts PeptideIdentifications of consensus features with respect to their map index.
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
    
    def __richcmp__(self, other: ConsensusMap, op: int) -> Any:
        ...
    
    def __iter__(self) -> ConsensusFeature:
       ... 


class DigestionEnzymeRNA:
    """
    Cython implementation of _DigestionEnzymeRNA

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzymeRNA.html>`_
      -- Inherits from ['DigestionEnzyme']

    Representation of a digestion enzyme for RNA (RNase)
    
    The cutting sites of these enzymes are defined using two different mechanisms:
    First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
    This is the same mechanism that is used for proteases (see ProteaseDigestion).
    However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
    Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
    If both expressions match, then there is a cutting site between the two ribonucleotides.
    
    There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
    A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DigestionEnzymeRNA()
        """
        ...
    
    @overload
    def __init__(self, in_0: DigestionEnzymeRNA ) -> None:
        """
        Cython signature: void DigestionEnzymeRNA(DigestionEnzymeRNA &)
        """
        ...
    
    def setCutsAfterRegEx(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCutsAfterRegEx(String value)
        Sets the "cuts after ..." regular expression
        """
        ...
    
    def getCutsAfterRegEx(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCutsAfterRegEx()
        Returns the "cuts after ..." regular expression
        """
        ...
    
    def setCutsBeforeRegEx(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCutsBeforeRegEx(String value)
        Sets the "cuts before ..." regular expression
        """
        ...
    
    def getCutsBeforeRegEx(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCutsBeforeRegEx()
        Returns the "cuts before ..." regular expression
        """
        ...
    
    def setThreePrimeGain(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setThreePrimeGain(String value)
        Sets the 3' gain
        """
        ...
    
    def setFivePrimeGain(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFivePrimeGain(String value)
        Sets the 5' gain
        """
        ...
    
    def getThreePrimeGain(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getThreePrimeGain()
        Returns the 3' gain
        """
        ...
    
    def getFivePrimeGain(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFivePrimeGain()
        Returns the 5' gain
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String & name)
        Sets the name of the enzyme
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the enzyme
        """
        ...
    
    def setSynonyms(self, synonyms: Set[bytes] ) -> None:
        """
        Cython signature: void setSynonyms(libcpp_set[String] & synonyms)
        Sets the synonyms
        """
        ...
    
    def addSynonym(self, synonym: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addSynonym(const String & synonym)
        Adds a synonym
        """
        ...
    
    def getSynonyms(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getSynonyms()
        Returns the synonyms
        """
        ...
    
    def setRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setRegEx(const String & cleavage_regex)
        Sets the cleavage regex
        """
        ...
    
    def getRegEx(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getRegEx()
        Returns the cleavage regex
        """
        ...
    
    def setRegExDescription(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setRegExDescription(const String & value)
        Sets the regex description
        """
        ...
    
    def getRegExDescription(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getRegExDescription()
        Returns the regex description
        """
        ...
    
    def setValueFromFile(self, key: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool setValueFromFile(String key, String value)
        Sets the value of a member variable based on an entry from an input file
        """
        ...
    
    def __richcmp__(self, other: DigestionEnzymeRNA, op: int) -> Any:
        ... 


class FeatureGroupingAlgorithmUnlabeled:
    """
    Cython implementation of _FeatureGroupingAlgorithmUnlabeled

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmUnlabeled.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureGroupingAlgorithmUnlabeled()
        """
        ...
    
    def group(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None:
        """
        Cython signature: void group(libcpp_vector[FeatureMap] & maps, ConsensusMap & out)
        """
        ...
    
    def addToGroup(self, map_id: int , feature_map: FeatureMap ) -> None:
        """
        Cython signature: void addToGroup(int map_id, FeatureMap feature_map)
        """
        ...
    
    def setReference(self, map_id: int , map: FeatureMap ) -> None:
        """
        Cython signature: void setReference(int map_id, FeatureMap map)
        """
        ...
    
    def getResultMap(self) -> ConsensusMap:
        """
        Cython signature: ConsensusMap getResultMap()
        """
        ...
    
    def transferSubelements(self, maps: List[ConsensusMap] , out: ConsensusMap ) -> None:
        """
        Cython signature: void transferSubelements(libcpp_vector[ConsensusMap] maps, ConsensusMap & out)
        Transfers subelements (grouped features) from input consensus maps to the result consensus map
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


class IdXMLFile:
    """
    Cython implementation of _IdXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IdXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void IdXMLFile()
        Used to load and store idXML files
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids)
        Loads the identifications of an idXML file without identifier
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids)
        Loads the identifications of an idXML file without identifier using PeptideIdentificationList
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: String ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String & document_id)
        Loads the identifications of an idXML file with identifier using PeptideIdentificationList
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String document_id)
        Stores the data in an idXML file
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids)
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String document_id)
        Stores the data in an idXML file using PeptideIdentificationList
        """
        ... 


class IndexedMzMLFileLoader:
    """
    Cython implementation of _IndexedMzMLFileLoader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLFileLoader.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void IndexedMzMLFileLoader()
        A class to load an indexedmzML file
        """
        ...
    
    def load(self, in_0: Union[bytes, str, String] , in_1: OnDiscMSExperiment ) -> bool:
        """
        Cython signature: bool load(String, OnDiscMSExperiment &)
        Load a file\n
        
        Tries to parse the file, success needs to be checked with the return value
        """
        ...
    
    @overload
    def store(self, in_0: Union[bytes, str, String] , in_1: OnDiscMSExperiment ) -> None:
        """
        Cython signature: void store(String, OnDiscMSExperiment &)
        Store a file from an on-disc data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
        """
        ...
    
    @overload
    def store(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void store(String, MSExperiment &)
        Store a file from an in-memory data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
        Returns the options for loading/storing
        """
        ...
    
    def setOptions(self, in_0: PeakFileOptions ) -> None:
        """
        Cython signature: void setOptions(PeakFileOptions)
        Returns the options for loading/storing
        """
        ... 


class IsobaricChannelInformation:
    """
    Cython implementation of _IsobaricChannelInformation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IsobaricQuantitationMethod_1_1IsobaricChannelInformation.html>`_
    """
    
    name: Union[bytes, str, String]
    
    id: int
    
    description: Union[bytes, str, String]
    
    center: float
    
    affected_channels: List[int]
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , id_: int , description: Union[bytes, str, String] , center: float , affected_channels: List[int] ) -> None:
        """
        Cython signature: void IsobaricChannelInformation(String name, int id_, String description, double center, libcpp_vector[int] affected_channels)
        """
        ...
    
    @overload
    def __init__(self, in_0: IsobaricChannelInformation ) -> None:
        """
        Cython signature: void IsobaricChannelInformation(IsobaricChannelInformation &)
        """
        ... 


class ItraqFourPlexQuantitationMethod:
    """
    Cython implementation of _ItraqFourPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqFourPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ItraqFourPlexQuantitationMethod()
        iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation
        """
        ...
    
    @overload
    def __init__(self, in_0: ItraqFourPlexQuantitationMethod ) -> None:
        """
        Cython signature: void ItraqFourPlexQuantitationMethod(ItraqFourPlexQuantitationMethod &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def getChannelInformation(self) -> List[IsobaricChannelInformation]:
        """
        Cython signature: libcpp_vector[IsobaricChannelInformation] getChannelInformation()
        """
        ...
    
    def getNumberOfChannels(self) -> int:
        """
        Cython signature: size_t getNumberOfChannels()
        """
        ...
    
    def getIsotopeCorrectionMatrix(self) -> MatrixDouble:
        """
        Cython signature: MatrixDouble getIsotopeCorrectionMatrix()
        """
        ...
    
    def getReferenceChannel(self) -> int:
        """
        Cython signature: size_t getReferenceChannel()
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
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class MSDataCachedConsumer:
    """
    Cython implementation of _MSDataCachedConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataCachedConsumer.html>`_

    Transforming and cached writing consumer of MS data
    
    Is able to transform a spectrum on the fly while it is read using a
    function pointer that can be set on the object. The spectra is then
    cached to disk using the functions provided in CachedMzMLHandler.
    """
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void MSDataCachedConsumer(String filename)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , clear: bool ) -> None:
        """
        Cython signature: void MSDataCachedConsumer(String filename, bool clear)
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        Write a spectrum to the output file
        
        May delete data from spectrum (if clearData is set)
        """
        ...
    
    def consumeChromatogram(self, c: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram & c)
        Write a chromatogram to the output file
        
        May delete data from chromatogram (if clearData is set)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        """
        ... 


class NonNegativeLeastSquaresSolver:
    """
    Cython implementation of _NonNegativeLeastSquaresSolver

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NonNegativeLeastSquaresSolver.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NonNegativeLeastSquaresSolver()
        """
        ...
    
    @overload
    def __init__(self, in_0: NonNegativeLeastSquaresSolver ) -> None:
        """
        Cython signature: void NonNegativeLeastSquaresSolver(NonNegativeLeastSquaresSolver &)
        """
        ...
    
    def solve(self, A: MatrixDouble , b: MatrixDouble , x: MatrixDouble ) -> int:
        """
        Cython signature: int solve(MatrixDouble & A, MatrixDouble & b, MatrixDouble & x)
        """
        ...
    RETURN_STATUS : __RETURN_STATUS 


class SplineInterpolatedPeaks:
    """
    Cython implementation of _SplineInterpolatedPeaks

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplineInterpolatedPeaks.html>`_
    """
    
    @overload
    def __init__(self, mz: List[float] , intensity: List[float] ) -> None:
        """
        Cython signature: void SplineInterpolatedPeaks(libcpp_vector[double] mz, libcpp_vector[double] intensity)
        """
        ...
    
    @overload
    def __init__(self, raw_spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void SplineInterpolatedPeaks(MSSpectrum raw_spectrum)
        """
        ...
    
    @overload
    def __init__(self, raw_chromatogram: MSChromatogram ) -> None:
        """
        Cython signature: void SplineInterpolatedPeaks(MSChromatogram raw_chromatogram)
        """
        ...
    
    @overload
    def __init__(self, in_0: SplineInterpolatedPeaks ) -> None:
        """
        Cython signature: void SplineInterpolatedPeaks(SplineInterpolatedPeaks &)
        """
        ...
    
    def getPosMin(self) -> float:
        """
        Cython signature: double getPosMin()
        """
        ...
    
    def getPosMax(self) -> float:
        """
        Cython signature: double getPosMax()
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def getNavigator(self, scaling: float ) -> SplineSpectrum_Navigator:
        """
        Cython signature: SplineSpectrum_Navigator getNavigator(double scaling)
        """
        ... 


class SplineSpectrum_Navigator:
    """
    Cython implementation of _SplineSpectrum_Navigator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplineSpectrum_Navigator.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SplineSpectrum_Navigator()
        """
        ...
    
    @overload
    def __init__(self, in_0: SplineSpectrum_Navigator ) -> None:
        """
        Cython signature: void SplineSpectrum_Navigator(SplineSpectrum_Navigator)
        """
        ...
    
    @overload
    def __init__(self, packages: List[SplinePackage] , posMax: float , scaling: float ) -> None:
        """
        Cython signature: void SplineSpectrum_Navigator(libcpp_vector[SplinePackage] * packages, double posMax, double scaling)
        """
        ...
    
    def eval(self, pos: float ) -> float:
        """
        Cython signature: double eval(double pos)
        """
        ...
    
    def getNextPos(self, pos: float ) -> float:
        """
        Cython signature: double getNextPos(double pos)
        """
        ... 


class TraMLFile:
    """
    Cython implementation of _TraMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TraMLFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TraMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: TraMLFile ) -> None:
        """
        Cython signature: void TraMLFile(TraMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , id: TargetedExperiment ) -> None:
        """
        Cython signature: void load(String filename, TargetedExperiment & id)
        Loads a map from a TraML file
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , id: TargetedExperiment ) -> None:
        """
        Cython signature: void store(String filename, TargetedExperiment & id)
        Stores a map in a TraML file
        """
        ...
    
    def isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool:
        """
        Cython signature: bool isSemanticallyValid(String filename, StringList & errors, StringList & warnings)
        Checks if a file is valid with respect to the mapping file and the controlled vocabulary
        
        :param filename: File name of the file to be checked
        :param errors: Errors during the validation are returned in this output parameter
        :param warnings: Warnings during the validation are returned in this output parameter
        """
        ... 


class __RETURN_STATUS:
    None
    SOLVED : int
    ITERATION_EXCEEDED : int

    def getMapping(self) -> Dict[int, str]:
       ... 

