from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_VersionDetails_create(in_0: Union[bytes, str, String] ) -> VersionDetails:
    """
    Cython signature: VersionDetails create(String)
    """
    ...

def __static_VersionInfo_getBranch() -> Union[bytes, str, String]:
    """
    Cython signature: String getBranch()
    """
    ...

def __static_VersionInfo_getRevision() -> Union[bytes, str, String]:
    """
    Cython signature: String getRevision()
    """
    ...

def __static_VersionInfo_getTime() -> Union[bytes, str, String]:
    """
    Cython signature: String getTime()
    """
    ...

def __static_VersionInfo_getVersion() -> Union[bytes, str, String]:
    """
    Cython signature: String getVersion()
    """
    ...

def __static_VersionInfo_getVersionStruct() -> VersionDetails:
    """
    Cython signature: VersionDetails getVersionStruct()
    """
    ...


class ConsensusIDAlgorithmPEPMatrix:
    """
    Cython implementation of _ConsensusIDAlgorithmPEPMatrix

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmPEPMatrix.html>`_
      -- Inherits from ['ConsensusIDAlgorithmSimilarity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmPEPMatrix()
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


class DRange1:
    """
    Cython implementation of _DRange1

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DRange1.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DRange1()
        """
        ...
    
    @overload
    def __init__(self, in_0: DRange1 ) -> None:
        """
        Cython signature: void DRange1(DRange1 &)
        """
        ...
    
    @overload
    def __init__(self, lower: DPosition1 , upper: DPosition1 ) -> None:
        """
        Cython signature: void DRange1(DPosition1 lower, DPosition1 upper)
        """
        ...
    
    def encloses(self, position: DPosition1 ) -> bool:
        """
        Cython signature: bool encloses(DPosition1 & position)
        """
        ...
    
    def united(self, other_range: DRange1 ) -> DRange1:
        """
        Cython signature: DRange1 united(DRange1 other_range)
        """
        ...
    
    def isIntersected(self, range_: DRange1 ) -> bool:
        """
        Cython signature: bool isIntersected(DRange1 & range_)
        """
        ...
    
    def isEmpty(self) -> bool:
        """
        Cython signature: bool isEmpty()
        """
        ...
    
    def __richcmp__(self, other: DRange1, op: int) -> Any:
        ... 


class DRange2:
    """
    Cython implementation of _DRange2

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DRange2.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DRange2()
        """
        ...
    
    @overload
    def __init__(self, in_0: DRange2 ) -> None:
        """
        Cython signature: void DRange2(DRange2 &)
        """
        ...
    
    @overload
    def __init__(self, lower: Union[Sequence[int], Sequence[float]] , upper: Union[Sequence[int], Sequence[float]] ) -> None:
        """
        Cython signature: void DRange2(DPosition2 lower, DPosition2 upper)
        """
        ...
    
    @overload
    def __init__(self, minx: float , miny: float , maxx: float , maxy: float ) -> None:
        """
        Cython signature: void DRange2(double minx, double miny, double maxx, double maxy)
        """
        ...
    
    def united(self, other_range: DRange2 ) -> DRange2:
        """
        Cython signature: DRange2 united(DRange2 other_range)
        """
        ...
    
    def isIntersected(self, range_: DRange2 ) -> bool:
        """
        Cython signature: bool isIntersected(DRange2 & range_)
        """
        ...
    
    def isEmpty(self) -> bool:
        """
        Cython signature: bool isEmpty()
        """
        ...
    
    def __richcmp__(self, other: DRange2, op: int) -> Any:
        ... 


class EnzymaticDigestion:
    """
    Cython implementation of _EnzymaticDigestion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EnzymaticDigestion.html>`_

      Class for the enzymatic digestion of proteins
    
      Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
      Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
      due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
      peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
      Thus no random selection of just n specific missed cleavage sites is performed.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EnzymaticDigestion()
        """
        ...
    
    @overload
    def __init__(self, in_0: EnzymaticDigestion ) -> None:
        """
        Cython signature: void EnzymaticDigestion(EnzymaticDigestion &)
        """
        ...
    
    def getMissedCleavages(self) -> int:
        """
        Cython signature: size_t getMissedCleavages()
        Returns the max. number of allowed missed cleavages for the digestion
        """
        ...
    
    def setMissedCleavages(self, missed_cleavages: int ) -> None:
        """
        Cython signature: void setMissedCleavages(size_t missed_cleavages)
        Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used
        """
        ...
    
    def countInternalCleavageSites(self, sequence: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t countInternalCleavageSites(String sequence)
        Returns the number of internal cleavage sites for this sequence.
        """
        ...
    
    def getEnzymeName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEnzymeName()
        Returns the enzyme for the digestion
        """
        ...
    
    def setEnzyme(self, enzyme: DigestionEnzyme ) -> None:
        """
        Cython signature: void setEnzyme(DigestionEnzyme * enzyme)
        Sets the enzyme for the digestion
        """
        ...
    
    def getSpecificity(self) -> int:
        """
        Cython signature: Specificity getSpecificity()
        Returns the specificity for the digestion
        """
        ...
    
    def setSpecificity(self, spec: int ) -> None:
        """
        Cython signature: void setSpecificity(Specificity spec)
        Sets the specificity for the digestion (default is SPEC_FULL)
        """
        ...
    
    def getSpecificityByName(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: Specificity getSpecificityByName(String name)
        Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid
        """
        ...
    
    def digestUnmodified(self, sequence: StringView , output: List[StringView] , min_length: int , max_length: int ) -> int:
        """
        Cython signature: size_t digestUnmodified(StringView sequence, libcpp_vector[StringView] & output, size_t min_length, size_t max_length)
        Performs the enzymatic digestion of an unmodified sequence\n
        By returning only references into the original string this is very fast
        
        
        :param sequence: Sequence to digest
        :param output: Digestion products
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        ...
    
    def isValidProduct(self, sequence: Union[bytes, str, String] , pos: int , length: int , ignore_missed_cleavages: bool ) -> bool:
        """
        Cython signature: bool isValidProduct(String sequence, int pos, int length, bool ignore_missed_cleavages)
        Boolean operator returns true if the peptide fragment starting at position `pos` with length `length` within the sequence `sequence` generated by the current enzyme\n
        Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here
        
        
        :param protein: Protein sequence
        :param pep_pos: Starting index of potential peptide
        :param pep_length: Length of potential peptide
        :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
        :return: True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
        """
        ...
    Specificity : __Specificity 


class FeatureGroupingAlgorithmQT:
    """
    Cython implementation of _FeatureGroupingAlgorithmQT

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmQT.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureGroupingAlgorithmQT()
        """
        ...
    
    @overload
    def group(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None:
        """
        Cython signature: void group(libcpp_vector[FeatureMap] & maps, ConsensusMap & out)
        """
        ...
    
    @overload
    def group(self, maps: List[ConsensusMap] , out: ConsensusMap ) -> None:
        """
        Cython signature: void group(libcpp_vector[ConsensusMap] & maps, ConsensusMap & out)
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


class LightCompound:
    """
    Cython implementation of _LightCompound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightCompound.html>`_
    """
    
    rt: float
    
    drift_time: float
    
    charge: int
    
    sequence: bytes
    
    protein_refs: List[bytes]
    
    peptide_group_label: bytes
    
    id: bytes
    
    sum_formula: bytes
    
    compound_name: bytes
    
    modifications: List[LightModification]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightCompound()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightCompound ) -> None:
        """
        Cython signature: void LightCompound(LightCompound &)
        """
        ...
    
    def setDriftTime(self, d: float ) -> None:
        """
        Cython signature: void setDriftTime(double d)
        """
        ...
    
    def getDriftTime(self) -> float:
        """
        Cython signature: double getDriftTime()
        """
        ...
    
    def getChargeState(self) -> int:
        """
        Cython signature: int getChargeState()
        """
        ...
    
    def isPeptide(self) -> bool:
        """
        Cython signature: bool isPeptide()
        """
        ...
    
    def setChargeState(self, ch: int ) -> None:
        """
        Cython signature: void setChargeState(int ch)
        """
        ... 


class LightMRMTransitionGroupCP:
    """
    Cython implementation of _MRMTransitionGroup[_MSChromatogram,_LightTransition]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroup[_MSChromatogram,_LightTransition].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightMRMTransitionGroupCP()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightMRMTransitionGroupCP ) -> None:
        """
        Cython signature: void LightMRMTransitionGroupCP(LightMRMTransitionGroupCP &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def getTransitionGroupID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTransitionGroupID()
        """
        ...
    
    def setTransitionGroupID(self, tr_gr_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTransitionGroupID(String tr_gr_id)
        """
        ...
    
    def getTransitions(self) -> List[LightTransition]:
        """
        Cython signature: libcpp_vector[LightTransition] getTransitions()
        """
        ...
    
    def getTransitionsMuteable(self) -> List[LightTransition]:
        """
        Cython signature: libcpp_vector[LightTransition] getTransitionsMuteable()
        """
        ...
    
    def addTransition(self, transition: LightTransition , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addTransition(LightTransition transition, String key)
        """
        ...
    
    def getTransition(self, key: Union[bytes, str, String] ) -> LightTransition:
        """
        Cython signature: LightTransition getTransition(String key)
        """
        ...
    
    def hasTransition(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasTransition(String key)
        """
        ...
    
    def getChromatograms(self) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] getChromatograms()
        """
        ...
    
    def addChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addChromatogram(MSChromatogram chromatogram, String key)
        """
        ...
    
    def getChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getChromatogram(String key)
        """
        ...
    
    def hasChromatogram(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasChromatogram(String key)
        """
        ...
    
    def getPrecursorChromatograms(self) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] getPrecursorChromatograms()
        """
        ...
    
    def addPrecursorChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addPrecursorChromatogram(MSChromatogram chromatogram, String key)
        """
        ...
    
    def getPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getPrecursorChromatogram(String key)
        """
        ...
    
    def hasPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasPrecursorChromatogram(String key)
        """
        ...
    
    def getFeatures(self) -> List[MRMFeature]:
        """
        Cython signature: libcpp_vector[MRMFeature] getFeatures()
        """
        ...
    
    def getFeaturesMuteable(self) -> List[MRMFeature]:
        """
        Cython signature: libcpp_vector[MRMFeature] getFeaturesMuteable()
        """
        ...
    
    def addFeature(self, feature: MRMFeature ) -> None:
        """
        Cython signature: void addFeature(MRMFeature feature)
        """
        ...
    
    def getBestFeature(self) -> MRMFeature:
        """
        Cython signature: MRMFeature getBestFeature()
        """
        ...
    
    def getLibraryIntensity(self, result: List[float] ) -> None:
        """
        Cython signature: void getLibraryIntensity(libcpp_vector[double] result)
        """
        ...
    
    def subset(self, tr_ids: List[Union[bytes, str]] ) -> LightMRMTransitionGroupCP:
        """
        Cython signature: LightMRMTransitionGroupCP subset(libcpp_vector[libcpp_utf8_string] tr_ids)
        """
        ...
    
    def isInternallyConsistent(self) -> bool:
        """
        Cython signature: bool isInternallyConsistent()
        """
        ...
    
    def chromatogramIdsMatch(self) -> bool:
        """
        Cython signature: bool chromatogramIdsMatch()
        """
        ... 


class LightModification:
    """
    Cython implementation of _LightModification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightModification.html>`_
    """
    
    location: int
    
    unimod_id: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightModification()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightModification ) -> None:
        """
        Cython signature: void LightModification(LightModification &)
        """
        ... 


class LightProtein:
    """
    Cython implementation of _LightProtein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightProtein.html>`_
    """
    
    id: bytes
    
    sequence: bytes
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightProtein()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightProtein ) -> None:
        """
        Cython signature: void LightProtein(LightProtein &)
        """
        ... 


class LightTargetedExperiment:
    """
    Cython implementation of _LightTargetedExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightTargetedExperiment.html>`_
    """
    
    transitions: List[LightTransition]
    
    compounds: List[LightCompound]
    
    proteins: List[LightProtein]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightTargetedExperiment()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightTargetedExperiment ) -> None:
        """
        Cython signature: void LightTargetedExperiment(LightTargetedExperiment &)
        """
        ...
    
    def getTransitions(self) -> List[LightTransition]:
        """
        Cython signature: libcpp_vector[LightTransition] getTransitions()
        """
        ...
    
    def getCompounds(self) -> List[LightCompound]:
        """
        Cython signature: libcpp_vector[LightCompound] getCompounds()
        """
        ...
    
    def getProteins(self) -> List[LightProtein]:
        """
        Cython signature: libcpp_vector[LightProtein] getProteins()
        """
        ...
    
    def getCompoundByRef(self, ref: bytes ) -> LightCompound:
        """
        Cython signature: LightCompound getCompoundByRef(libcpp_string & ref)
        """
        ...
    
    def getPeptideByRef(self, ref: bytes ) -> LightCompound:
        """
        Cython signature: LightCompound getPeptideByRef(libcpp_string & ref)
        """
        ... 


class LightTransition:
    """
    Cython implementation of _LightTransition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightTransition.html>`_
    """
    
    transition_name: bytes
    
    transition_ref: bytes
    
    library_intensity: float
    
    product_mz: float
    
    precursor_mz: float
    
    fragment_charge: int
    
    decoy: bool
    
    detecting_transition: bool
    
    quantifying_transition: bool
    
    identifying_transition: bool
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LightTransition()
        """
        ...
    
    @overload
    def __init__(self, in_0: LightTransition ) -> None:
        """
        Cython signature: void LightTransition(LightTransition &)
        """
        ...
    
    def getProductChargeState(self) -> int:
        """
        Cython signature: int getProductChargeState()
        """
        ...
    
    def isProductChargeStateSet(self) -> bool:
        """
        Cython signature: bool isProductChargeStateSet()
        """
        ...
    
    def getNativeID(self) -> bytes:
        """
        Cython signature: libcpp_string getNativeID()
        """
        ...
    
    def getTransRef(self) -> bytes:
        """
        Cython signature: libcpp_string getTransRef()
        """
        ...
    
    def getLibraryIntensity(self) -> float:
        """
        Cython signature: double getLibraryIntensity()
        """
        ...
    
    def setLibraryIntensity(self, l: float ) -> None:
        """
        Cython signature: void setLibraryIntensity(double l)
        """
        ...
    
    def getProductMZ(self) -> float:
        """
        Cython signature: double getProductMZ()
        """
        ...
    
    def getPrecursorMZ(self) -> float:
        """
        Cython signature: double getPrecursorMZ()
        """
        ...
    
    def setDetectingTransition(self, d: bool ) -> None:
        """
        Cython signature: void setDetectingTransition(bool d)
        """
        ...
    
    def isDetectingTransition(self) -> bool:
        """
        Cython signature: bool isDetectingTransition()
        """
        ...
    
    def setQuantifyingTransition(self, q: bool ) -> None:
        """
        Cython signature: void setQuantifyingTransition(bool q)
        """
        ...
    
    def isQuantifyingTransition(self) -> bool:
        """
        Cython signature: bool isQuantifyingTransition()
        """
        ...
    
    def setIdentifyingTransition(self, i: bool ) -> None:
        """
        Cython signature: void setIdentifyingTransition(bool i)
        """
        ...
    
    def isIdentifyingTransition(self) -> bool:
        """
        Cython signature: bool isIdentifyingTransition()
        """
        ... 


class MRMTransitionGroupCP:
    """
    Cython implementation of _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMTransitionGroupCP()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMTransitionGroupCP ) -> None:
        """
        Cython signature: void MRMTransitionGroupCP(MRMTransitionGroupCP &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def getTransitionGroupID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTransitionGroupID()
        """
        ...
    
    def setTransitionGroupID(self, tr_gr_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTransitionGroupID(String tr_gr_id)
        """
        ...
    
    def getTransitions(self) -> List[ReactionMonitoringTransition]:
        """
        Cython signature: libcpp_vector[ReactionMonitoringTransition] getTransitions()
        """
        ...
    
    def getTransitionsMuteable(self) -> List[ReactionMonitoringTransition]:
        """
        Cython signature: libcpp_vector[ReactionMonitoringTransition] getTransitionsMuteable()
        """
        ...
    
    def addTransition(self, transition: ReactionMonitoringTransition , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addTransition(ReactionMonitoringTransition transition, String key)
        """
        ...
    
    def getTransition(self, key: Union[bytes, str, String] ) -> ReactionMonitoringTransition:
        """
        Cython signature: ReactionMonitoringTransition getTransition(String key)
        """
        ...
    
    def hasTransition(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasTransition(String key)
        """
        ...
    
    def getChromatograms(self) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] getChromatograms()
        """
        ...
    
    def addChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addChromatogram(MSChromatogram chromatogram, String key)
        """
        ...
    
    def getChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getChromatogram(String key)
        """
        ...
    
    def hasChromatogram(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasChromatogram(String key)
        """
        ...
    
    def getPrecursorChromatograms(self) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] getPrecursorChromatograms()
        """
        ...
    
    def addPrecursorChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addPrecursorChromatogram(MSChromatogram chromatogram, String key)
        """
        ...
    
    def getPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getPrecursorChromatogram(String key)
        """
        ...
    
    def hasPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasPrecursorChromatogram(String key)
        """
        ...
    
    def getFeatures(self) -> List[MRMFeature]:
        """
        Cython signature: libcpp_vector[MRMFeature] getFeatures()
        """
        ...
    
    def getFeaturesMuteable(self) -> List[MRMFeature]:
        """
        Cython signature: libcpp_vector[MRMFeature] getFeaturesMuteable()
        """
        ...
    
    def addFeature(self, feature: MRMFeature ) -> None:
        """
        Cython signature: void addFeature(MRMFeature feature)
        """
        ...
    
    def getBestFeature(self) -> MRMFeature:
        """
        Cython signature: MRMFeature getBestFeature()
        """
        ...
    
    def getLibraryIntensity(self, result: List[float] ) -> None:
        """
        Cython signature: void getLibraryIntensity(libcpp_vector[double] result)
        """
        ...
    
    def subset(self, tr_ids: List[Union[bytes, str]] ) -> MRMTransitionGroupCP:
        """
        Cython signature: MRMTransitionGroupCP subset(libcpp_vector[libcpp_utf8_string] tr_ids)
        """
        ...
    
    def isInternallyConsistent(self) -> bool:
        """
        Cython signature: bool isInternallyConsistent()
        """
        ...
    
    def chromatogramIdsMatch(self) -> bool:
        """
        Cython signature: bool chromatogramIdsMatch()
        """
        ... 


class MzTabM:
    """
    Cython implementation of _MzTabM

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabM.html>`_

    Data model of MzTabM files
    
    Please see the official MzTabM specification at https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzTabM()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzTabM ) -> None:
        """
        Cython signature: void MzTabM(MzTabM &)
        """
        ...
    
    def exportFeatureMapToMzTabM(self, feature_map: FeatureMap ) -> MzTabM:
        """
        Cython signature: MzTabM exportFeatureMapToMzTabM(FeatureMap feature_map)
        Export FeatureMap with Identifications to MzTabM
        """
        ... 


class PeptideIndexing:
    """
    Cython implementation of _PeptideIndexing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIndexing.html>`_
      -- Inherits from ['DefaultParamHandler']

    Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information
    
    All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy",
    depending on whether the protein accession contains the decoy pattern (parameter `decoy_string`) as a suffix or prefix, respectively (see parameter `prefix`).
    For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
    only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
    (For FDR calculations, "target+decoy" peptide hits count as target hits.)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideIndexing()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideIndexing ) -> None:
        """
        Cython signature: void PeptideIndexing(PeptideIndexing &)
        """
        ...
    
    def run(self, proteins: List[FASTAEntry] , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> int:
        """
        Cython signature: PeptideIndexing_ExitCodes run(libcpp_vector[FASTAEntry] & proteins, libcpp_vector[ProteinIdentification] & prot_ids, PeptideIdentificationList & pep_ids)
        """
        ...
    
    def getDecoyString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDecoyString()
        """
        ...
    
    def isPrefix(self) -> bool:
        """
        Cython signature: bool isPrefix()
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
    PeptideIndexing_ExitCodes : __PeptideIndexing_ExitCodes 


class QTClusterFinder:
    """
    Cython implementation of _QTClusterFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QTClusterFinder.html>`_
      -- Inherits from ['BaseGroupFinder']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void QTClusterFinder()
        """
        ...
    
    @overload
    def run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None:
        """
        Cython signature: void run(libcpp_vector[ConsensusMap] & input_maps, ConsensusMap & result_map)
        """
        ...
    
    @overload
    def run(self, input_maps: List[FeatureMap] , result_map: ConsensusMap ) -> None:
        """
        Cython signature: void run(libcpp_vector[FeatureMap] & input_maps, ConsensusMap & result_map)
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


class RNaseDB:
    """
    Cython implementation of _RNaseDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RNaseDB.html>`_
    """
    
    def getEnzyme(self, name: Union[bytes, str, String] ) -> DigestionEnzymeRNA:
        """
        Cython signature: const DigestionEnzymeRNA * getEnzyme(const String & name)
        """
        ...
    
    def getEnzymeByRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> DigestionEnzymeRNA:
        """
        Cython signature: const DigestionEnzymeRNA * getEnzymeByRegEx(const String & cleavage_regex)
        """
        ...
    
    def getAllNames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllNames(libcpp_vector[String] & all_names)
        """
        ...
    
    def hasEnzyme(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasEnzyme(const String & name)
        """
        ... 


class Ribonucleotide:
    """
    Cython implementation of _Ribonucleotide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Ribonucleotide_1_1Ribonucleotide.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Ribonucleotide()
        """
        ...
    
    @overload
    def __init__(self, in_0: Ribonucleotide ) -> None:
        """
        Cython signature: void Ribonucleotide(Ribonucleotide &)
        """
        ...
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , code: Union[bytes, str, String] , new_code: Union[bytes, str, String] , html_code: Union[bytes, str, String] , formula: EmpiricalFormula , origin: bytes , mono_mass: float , avg_mass: float , term_spec: int , baseloss_formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void Ribonucleotide(String name, String code, String new_code, String html_code, EmpiricalFormula formula, char origin, double mono_mass, double avg_mass, TermSpecificityNuc term_spec, EmpiricalFormula baseloss_formula)
        """
        ...
    
    def getCode(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCode()
        Returns the short name
        """
        ...
    
    def setCode(self, code: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCode(String code)
        Sets the short name
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the ribonucleotide
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the ribonucleotide
        """
        ...
    
    def setFormula(self, formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void setFormula(EmpiricalFormula formula)
        Sets empirical formula of the ribonucleotide (must be full, with N and C-terminus)
        """
        ...
    
    def getFormula(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula()
        Returns the empirical formula of the residue
        """
        ...
    
    def setAvgMass(self, avg_mass: float ) -> None:
        """
        Cython signature: void setAvgMass(double avg_mass)
        Sets average mass of the ribonucleotide
        """
        ...
    
    def getAvgMass(self) -> float:
        """
        Cython signature: double getAvgMass()
        Returns average mass of the ribonucleotide
        """
        ...
    
    def setMonoMass(self, mono_mass: float ) -> None:
        """
        Cython signature: void setMonoMass(double mono_mass)
        Sets monoisotopic mass of the ribonucleotide
        """
        ...
    
    def getMonoMass(self) -> float:
        """
        Cython signature: double getMonoMass()
        Returns monoisotopic mass of the ribonucleotide
        """
        ...
    
    def getNewCode(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNewCode()
        Returns the new code
        """
        ...
    
    def setNewCode(self, code: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNewCode(String code)
        Sets the new code
        """
        ...
    
    def getOrigin(self) -> bytes:
        """
        Cython signature: char getOrigin()
        Returns the code of the unmodified base (e.g., "A", "C", ...)
        """
        ...
    
    def setOrigin(self, origin: bytes ) -> None:
        """
        Cython signature: void setOrigin(char origin)
        Sets the code of the unmodified base (e.g., "A", "C", ...)
        """
        ...
    
    def setHTMLCode(self, html_code: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setHTMLCode(String html_code)
        Sets the HTML (RNAMods) code
        """
        ...
    
    def getHTMLCode(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getHTMLCode()
        Returns the HTML (RNAMods) code
        """
        ...
    
    def setTermSpecificity(self, term_spec: int ) -> None:
        """
        Cython signature: void setTermSpecificity(TermSpecificityNuc term_spec)
        Sets the terminal specificity
        """
        ...
    
    def getTermSpecificity(self) -> int:
        """
        Cython signature: TermSpecificityNuc getTermSpecificity()
        Returns the terminal specificity
        """
        ...
    
    def getBaselossFormula(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getBaselossFormula()
        Returns sum formula after loss of the nucleobase
        """
        ...
    
    def setBaselossFormula(self, formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void setBaselossFormula(EmpiricalFormula formula)
        Sets sum formula after loss of the nucleobase
        """
        ...
    
    def isModified(self) -> bool:
        """
        Cython signature: bool isModified()
        True if the ribonucleotide is a modified one
        """
        ...
    
    def __richcmp__(self, other: Ribonucleotide, op: int) -> Any:
        ...
    TermSpecificityNuc : __TermSpecificityNuc 


class TMTEighteenPlexQuantitationMethod:
    """
    Cython implementation of _TMTEighteenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTEighteenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TMTEighteenPlexQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: TMTEighteenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void TMTEighteenPlexQuantitationMethod(TMTEighteenPlexQuantitationMethod &)
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


class UnimodXMLFile:
    """
    Cython implementation of _UnimodXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1UnimodXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void UnimodXMLFile()
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
        """
        ... 


class VersionDetails:
    """
    Cython implementation of _VersionDetails

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1VersionDetails.html>`_
    """
    
    version_major: int
    
    version_minor: int
    
    version_patch: int
    
    pre_release_identifier: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void VersionDetails()
        """
        ...
    
    @overload
    def __init__(self, in_0: VersionDetails ) -> None:
        """
        Cython signature: void VersionDetails(VersionDetails &)
        """
        ...
    
    def __richcmp__(self, other: VersionDetails, op: int) -> Any:
        ...
    
    create: __static_VersionDetails_create 


class VersionInfo:
    """
    Cython implementation of _VersionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1VersionInfo.html>`_
    """
    
    getBranch: __static_VersionInfo_getBranch
    
    getRevision: __static_VersionInfo_getRevision
    
    getTime: __static_VersionInfo_getTime
    
    getVersion: __static_VersionInfo_getVersion
    
    getVersionStruct: __static_VersionInfo_getVersionStruct 


class DRangeIntersection:
    None
    Disjoint : int
    Intersects : int
    Inside : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __PeptideIndexing_ExitCodes:
    None
    EXECUTION_OK : int
    DATABASE_EMPTY : int
    PEPTIDE_IDS_EMPTY : int
    ILLEGAL_PARAMETERS : int
    UNEXPECTED_RESULT : int
    DECOYSTRING_EMPTY : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __Specificity:
    None
    SPEC_NONE : int
    SPEC_SEMI : int
    SPEC_FULL : int
    SPEC_UNKNOWN : int
    SPEC_NOCTERM : int
    SPEC_NONTERM : int
    SIZE_OF_SPECIFICITY : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __TermSpecificityNuc:
    None
    ANYWHERE : int
    FIVE_PRIME : int
    THREE_PRIME : int
    NUMBER_OF_TERM_SPECIFICITY : int

    def getMapping(self) -> Dict[int, str]:
       ... 

