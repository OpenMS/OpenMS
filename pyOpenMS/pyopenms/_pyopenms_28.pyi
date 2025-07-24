from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_TransformationModelLinear_getDefaultParameters(in_0: Param ) -> None:
    """
    Cython signature: void getDefaultParameters(Param &)
    """
    ...


class AbsoluteQuantitationStandardsFile:
    """
    Cython implementation of _AbsoluteQuantitationStandardsFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationStandardsFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AbsoluteQuantitationStandardsFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: AbsoluteQuantitationStandardsFile ) -> None:
        """
        Cython signature: void AbsoluteQuantitationStandardsFile(AbsoluteQuantitationStandardsFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , run_concentrations: List[AQS_runConcentration] ) -> None:
        """
        Cython signature: void load(const String & filename, libcpp_vector[AQS_runConcentration] & run_concentrations)
        """
        ... 


class ChannelInfo:
    """
    Cython implementation of _ChannelInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChannelInfo.html>`_
    """
    
    description: bytes
    
    name: int
    
    id: int
    
    center: float
    
    active: bool 


class DBoundingBox2:
    """
    Cython implementation of _DBoundingBox2

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DBoundingBox2.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DBoundingBox2()
        """
        ...
    
    @overload
    def __init__(self, in_0: DBoundingBox2 ) -> None:
        """
        Cython signature: void DBoundingBox2(DBoundingBox2 &)
        """
        ...
    
    def minPosition(self) -> Union[Sequence[int], Sequence[float]]:
        """
        Cython signature: DPosition2 minPosition()
        """
        ...
    
    def maxPosition(self) -> Union[Sequence[int], Sequence[float]]:
        """
        Cython signature: DPosition2 maxPosition()
        """
        ... 


class DigestionEnzymeProtein:
    """
    Cython implementation of _DigestionEnzymeProtein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzymeProtein.html>`_
      -- Inherits from ['DigestionEnzyme']

    Representation of a digestion enzyme for proteins (protease)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DigestionEnzymeProtein()
        """
        ...
    
    @overload
    def __init__(self, in_0: DigestionEnzymeProtein ) -> None:
        """
        Cython signature: void DigestionEnzymeProtein(DigestionEnzymeProtein &)
        """
        ...
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] , n_term_gain: EmpiricalFormula , c_term_gain: EmpiricalFormula , psi_id: Union[bytes, str, String] , xtandem_id: Union[bytes, str, String] , comet_id: int , omssa_id: int ) -> None:
        """
        Cython signature: void DigestionEnzymeProtein(String name, String cleavage_regex, libcpp_set[String] synonyms, String regex_description, EmpiricalFormula n_term_gain, EmpiricalFormula c_term_gain, String psi_id, String xtandem_id, unsigned int comet_id, unsigned int omssa_id)
        """
        ...
    
    def setNTermGain(self, value: EmpiricalFormula ) -> None:
        """
        Cython signature: void setNTermGain(EmpiricalFormula value)
        Sets the N-term gain
        """
        ...
    
    def setCTermGain(self, value: EmpiricalFormula ) -> None:
        """
        Cython signature: void setCTermGain(EmpiricalFormula value)
        Sets the C-term gain
        """
        ...
    
    def getNTermGain(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getNTermGain()
        Returns the N-term gain
        """
        ...
    
    def getCTermGain(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getCTermGain()
        Returns the C-term gain
        """
        ...
    
    def setPSIID(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPSIID(String value)
        Sets the PSI ID
        """
        ...
    
    def getPSIID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPSIID()
        Returns the PSI ID
        """
        ...
    
    def setXTandemID(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setXTandemID(String value)
        Sets the X! Tandem enzyme ID
        """
        ...
    
    def getXTandemID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getXTandemID()
        Returns the X! Tandem enzyme ID
        """
        ...
    
    def setCometID(self, value: int ) -> None:
        """
        Cython signature: void setCometID(int value)
        Sets the Comet enzyme ID
        """
        ...
    
    def getCometID(self) -> int:
        """
        Cython signature: int getCometID()
        Returns the Comet enzyme ID
        """
        ...
    
    def setOMSSAID(self, value: int ) -> None:
        """
        Cython signature: void setOMSSAID(int value)
        Sets the OMSSA enzyme ID
        """
        ...
    
    def getOMSSAID(self) -> int:
        """
        Cython signature: int getOMSSAID()
        Returns the OMSSA enzyme ID
        """
        ...
    
    def setMSGFID(self, value: int ) -> None:
        """
        Cython signature: void setMSGFID(int value)
        Sets the MSGFPlus enzyme id
        """
        ...
    
    def getMSGFID(self) -> int:
        """
        Cython signature: int getMSGFID()
        Returns the MSGFPlus enzyme id
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
    
    def __richcmp__(self, other: DigestionEnzymeProtein, op: int) -> Any:
        ... 


class FIAMSScheduler:
    """
    Cython implementation of _FIAMSScheduler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FIAMSScheduler.html>`_

      ADD PYTHON DOCUMENTATION HERE
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FIAMSScheduler()
        Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor
        """
        ...
    
    @overload
    def __init__(self, in_0: FIAMSScheduler ) -> None:
        """
        Cython signature: void FIAMSScheduler(FIAMSScheduler &)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , base_dir: Union[bytes, str, String] , load_cached_: bool ) -> None:
        """
        Cython signature: void FIAMSScheduler(String filename, String base_dir, bool load_cached_)
        """
        ...
    
    def run(self) -> None:
        """
        Cython signature: void run()
        Run the FIA-MS data analysis for the batch defined in the @filename_
        """
        ...
    
    def getBaseDir(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getBaseDir()
        Returns the base directory for the relevant paths from the csv file
        """
        ... 


class FeatureGroupingAlgorithmLabeled:
    """
    Cython implementation of _FeatureGroupingAlgorithmLabeled

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmLabeled.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureGroupingAlgorithmLabeled()
        """
        ...
    
    def group(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None:
        """
        Cython signature: void group(libcpp_vector[FeatureMap] & maps, ConsensusMap & out)
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


class IsobaricQuantifier:
    """
    Cython implementation of _IsobaricQuantifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricQuantifier.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, in_0: IsobaricQuantifier ) -> None:
        """
        Cython signature: void IsobaricQuantifier(IsobaricQuantifier &)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricQuantifier(ItraqFourPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricQuantifier(ItraqEightPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTSixPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricQuantifier(TMTSixPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTTenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricQuantifier(TMTTenPlexQuantitationMethod * quant_method)
        """
        ...
    
    def quantify(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap ) -> None:
        """
        Cython signature: void quantify(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out)
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


class ItraqConstants:
    """
    Cython implementation of _ItraqConstants

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqConstants.html>`_

    Some constants used throughout iTRAQ classes
    
    Constants for iTRAQ experiments and a ChannelInfo structure to store information about a single channel
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ItraqConstants()
        """
        ...
    
    @overload
    def __init__(self, in_0: ItraqConstants ) -> None:
        """
        Cython signature: void ItraqConstants(ItraqConstants &)
        """
        ...
    
    def getIsotopeMatrixAsStringList(self, itraq_type: int , isotope_corrections: List[MatrixDouble] ) -> List[bytes]:
        """
        Cython signature: StringList getIsotopeMatrixAsStringList(int itraq_type, libcpp_vector[MatrixDouble] & isotope_corrections)
        Convert isotope correction matrix to stringlist\n
        
        Each line is converted into a string of the format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
        Useful for creating parameters or debug output
        
        
        :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
        :param isotope_corrections: Vector of the two matrices (4plex, 8plex)
        """
        ...
    
    def updateIsotopeMatrixFromStringList(self, itraq_type: int , channels: List[bytes] , isotope_corrections: List[MatrixDouble] ) -> None:
        """
        Cython signature: void updateIsotopeMatrixFromStringList(int itraq_type, StringList & channels, libcpp_vector[MatrixDouble] & isotope_corrections)
        Convert strings to isotope correction matrix rows\n
        
        Each string of format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
        is parsed and the corresponding channel(row) in the matrix is updated
        Not all channels need to be present, missing channels will be left untouched
        Useful to update the matrix with user isotope correction values
        
        
        :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
        :param channels: New channel isotope values as strings
        :param isotope_corrections: Vector of the two matrices (4plex, 8plex)
        """
        ...
    
    def translateIsotopeMatrix(self, itraq_type: int , isotope_corrections: List[MatrixDouble] ) -> MatrixDouble:
        """
        Cython signature: MatrixDouble translateIsotopeMatrix(int & itraq_type, libcpp_vector[MatrixDouble] & isotope_corrections)
        """
        ... 


class MRMFeature:
    """
    Cython implementation of _MRMFeature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeature.html>`_
      -- Inherits from ['Feature']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeature()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeature ) -> None:
        """
        Cython signature: void MRMFeature(MRMFeature &)
        """
        ...
    
    def getScores(self) -> OpenSwath_Scores:
        """
        Cython signature: OpenSwath_Scores getScores()
        Returns all peakgroup scores
        """
        ...
    
    def setScores(self, s: OpenSwath_Scores ) -> None:
        """
        Cython signature: void setScores(OpenSwath_Scores s)
        Sets all peakgroup scores
        """
        ...
    
    def getFeature(self, key: Union[bytes, str, String] ) -> Feature:
        """
        Cython signature: Feature getFeature(String key)
        Returns a specified feature
        """
        ...
    
    def addFeature(self, f: Feature , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addFeature(Feature & f, String key)
        Adds an feature from a single chromatogram into the feature
        """
        ...
    
    def getFeatures(self) -> List[Feature]:
        """
        Cython signature: libcpp_vector[Feature] getFeatures()
        Returns all the features
        """
        ...
    
    def getFeatureIDs(self, result: List[bytes] ) -> None:
        """
        Cython signature: void getFeatureIDs(libcpp_vector[String] & result)
        Returns a list of IDs of available features
        """
        ...
    
    def getPrecursorFeature(self, key: Union[bytes, str, String] ) -> Feature:
        """
        Cython signature: Feature getPrecursorFeature(String key)
        Returns a specified precursor feature
        """
        ...
    
    def addPrecursorFeature(self, f: Feature , key: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addPrecursorFeature(Feature & f, String key)
        Adds a precursor feature from a single chromatogram into the feature
        """
        ...
    
    def getPrecursorFeatureIDs(self, result: List[bytes] ) -> None:
        """
        Cython signature: void getPrecursorFeatureIDs(libcpp_vector[String] & result)
        Returns a list of IDs of available precursor features
        """
        ...
    
    def getQuality(self, index: int ) -> float:
        """
        Cython signature: float getQuality(size_t index)
        Returns the quality in dimension c
        """
        ...
    
    def setQuality(self, index: int , q: float ) -> None:
        """
        Cython signature: void setQuality(size_t index, float q)
        Sets the quality in dimension c
        """
        ...
    
    def getOverallQuality(self) -> float:
        """
        Cython signature: float getOverallQuality()
        Model and quality methods
        """
        ...
    
    def setOverallQuality(self, q: float ) -> None:
        """
        Cython signature: void setOverallQuality(float q)
        Sets the overall quality
        """
        ...
    
    def getSubordinates(self) -> List[Feature]:
        """
        Cython signature: libcpp_vector[Feature] getSubordinates()
        Returns the subordinate features
        """
        ...
    
    def setSubordinates(self, in_0: List[Feature] ) -> None:
        """
        Cython signature: void setSubordinates(libcpp_vector[Feature])
        Returns the subordinate features
        """
        ...
    
    def encloses(self, rt: float , mz: float ) -> bool:
        """
        Cython signature: bool encloses(double rt, double mz)
        Returns if the mass trace convex hulls of the feature enclose the position specified by `rt` and `mz`
        
        
        :param rt: Sequence to digest
        :param mz: Digestion products
        """
        ...
    
    def getConvexHull(self) -> ConvexHull2D:
        """
        Cython signature: ConvexHull2D getConvexHull()
        """
        ...
    
    def getConvexHulls(self) -> List[ConvexHull2D]:
        """
        Cython signature: libcpp_vector[ConvexHull2D] getConvexHulls()
        """
        ...
    
    def setConvexHulls(self, in_0: List[ConvexHull2D] ) -> None:
        """
        Cython signature: void setConvexHulls(libcpp_vector[ConvexHull2D])
        """
        ...
    
    def getWidth(self) -> float:
        """
        Cython signature: float getWidth()
        """
        ...
    
    def setWidth(self, q: float ) -> None:
        """
        Cython signature: void setWidth(float q)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        """
        ...
    
    def setCharge(self, q: int ) -> None:
        """
        Cython signature: void setCharge(int q)
        """
        ...
    
    def getAnnotationState(self) -> int:
        """
        Cython signature: AnnotationState getAnnotationState()
        """
        ...
    
    def getPeptideIdentifications(self) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList getPeptideIdentifications()
        Returns a reference to the PeptideIdentification vector
        """
        ...
    
    def setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setPeptideIdentifications(PeptideIdentificationList & peptides)
        Sets the PeptideIdentification vector
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
    
    def __richcmp__(self, other: MRMFeature, op: int) -> Any:
        ... 


class MascotGenericFile:
    """
    Cython implementation of _MascotGenericFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MascotGenericFile.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MascotGenericFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MascotGenericFile ) -> None:
        """
        Cython signature: void MascotGenericFile(MascotGenericFile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , experiment: MSExperiment ) -> None:
        """
        Cython signature: void store(const String & filename, MSExperiment & experiment)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & exp)
        Loads a Mascot Generic File into a PeakMap
        
        
        :param filename: File name which the map should be read from
        :param exp: The map which is filled with the data from the given file
        :raises:
          Exception: FileNotFound is thrown if the given file could not be found
        """
        ...
    
    def getHTTPPeakListEnclosure(self, filename: Union[bytes, str, String] ) -> List[Union[bytes, str, String], Union[bytes, str, String]]:
        """
        Cython signature: libcpp_pair[String,String] getHTTPPeakListEnclosure(const String & filename)
        Enclosing Strings of the peak list body for HTTP submission\n
        
        Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
        adding the peaks (in whatever format, e.g. mzXML) enclosed in this body
        The `filename` can later be found in the Mascot response
        """
        ...
    
    def updateMembers_(self) -> None:
        """
        Cython signature: void updateMembers_()
        Docu in base class
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


class SeedListGenerator:
    """
    Cython implementation of _SeedListGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SeedListGenerator.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SeedListGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: SeedListGenerator ) -> None:
        """
        Cython signature: void SeedListGenerator(SeedListGenerator &)
        """
        ...
    
    @overload
    def generateSeedList(self, exp: MSExperiment , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None:
        """
        Cython signature: void generateSeedList(MSExperiment exp, libcpp_vector[DPosition2] & seeds)
        Generate a seed list based on an MS experiment
        """
        ...
    
    @overload
    def generateSeedList(self, peptides: PeptideIdentificationList , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , use_peptide_mass: bool ) -> None:
        """
        Cython signature: void generateSeedList(PeptideIdentificationList & peptides, libcpp_vector[DPosition2] & seeds, bool use_peptide_mass)
        Generates a seed list based on a list of peptide identifications
        """
        ...
    
    @overload
    def convertSeedList(self, seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , features: FeatureMap ) -> None:
        """
        Cython signature: void convertSeedList(libcpp_vector[DPosition2] & seeds, FeatureMap & features)
        Converts a list of seed positions to a feature map (expected format for FeatureFinder)
        """
        ...
    
    @overload
    def convertSeedList(self, features: FeatureMap , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None:
        """
        Cython signature: void convertSeedList(FeatureMap & features, libcpp_vector[DPosition2] & seeds)
        Converts a feature map with seed positions back to a simple list
        """
        ... 


class SignalToNoiseEstimatorMedian:
    """
    Cython implementation of _SignalToNoiseEstimatorMedian[_MSSpectrum]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMedian[_MSSpectrum].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMedian()
        """
        ...
    
    @overload
    def __init__(self, in_0: SignalToNoiseEstimatorMedian ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMedian(SignalToNoiseEstimatorMedian &)
        """
        ...
    
    def init(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void init(MSSpectrum & spectrum)
        """
        ...
    
    def getSignalToNoise(self, index: int ) -> float:
        """
        Cython signature: double getSignalToNoise(size_t index)
        """
        ...
    
    def getSparseWindowPercent(self) -> float:
        """
        Cython signature: double getSparseWindowPercent()
        """
        ...
    
    def getHistogramRightmostPercent(self) -> float:
        """
        Cython signature: double getHistogramRightmostPercent()
        """
        ...
    IntensityThresholdCalculation : __IntensityThresholdCalculation 


class SiriusMSFile:
    """
    Cython implementation of _SiriusMSFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusMSFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusMSFile ) -> None:
        """
        Cython signature: void SiriusMSFile(SiriusMSFile &)
        """
        ... 


class SiriusMSFile_AccessionInfo:
    """
    Cython implementation of _SiriusMSFile_AccessionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile_AccessionInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusMSFile_AccessionInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusMSFile_AccessionInfo ) -> None:
        """
        Cython signature: void SiriusMSFile_AccessionInfo(SiriusMSFile_AccessionInfo &)
        """
        ... 


class SiriusMSFile_CompoundInfo:
    """
    Cython implementation of _SiriusMSFile_CompoundInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile_CompoundInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusMSFile_CompoundInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusMSFile_CompoundInfo ) -> None:
        """
        Cython signature: void SiriusMSFile_CompoundInfo(SiriusMSFile_CompoundInfo &)
        """
        ... 


class SpectrumAccessOpenMSCached:
    """
    Cython implementation of _SpectrumAccessOpenMSCached

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMSCached.html>`_
      -- Inherits from ['ISpectrumAccess']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSCached()
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSCached(String filename)
        An implementation of the Spectrum Access interface using on-disk caching
        
        This class implements the OpenSWATH Spectrum Access interface
        (ISpectrumAccess) using the CachedmzML class which is able to read and
        write a cached mzML file
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMSCached ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSCached(SpectrumAccessOpenMSCached &)
        """
        ...
    
    def getSpectrumById(self, id_: int ) -> OSSpectrum:
        """
        Cython signature: shared_ptr[OSSpectrum] getSpectrumById(int id_)
        Returns a pointer to a spectrum at the given string id
        """
        ...
    
    def getSpectraByRT(self, RT: float , deltaRT: float ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT)
        Returns a vector of ids of spectra that are within RT +/- deltaRT
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        Returns the number of spectra available
        """
        ...
    
    def getChromatogramById(self, id_: int ) -> OSChromatogram:
        """
        Cython signature: shared_ptr[OSChromatogram] getChromatogramById(int id_)
        Returns a pointer to a chromatogram at the given id
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        Returns the number of chromatograms available
        """
        ...
    
    def getChromatogramNativeID(self, id_: int ) -> str:
        """
        Cython signature: libcpp_utf8_output_string getChromatogramNativeID(int id_)
        """
        ... 


class TransformationModelLinear:
    """
    Cython implementation of _TransformationModelLinear

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelLinear.html>`_
      -- Inherits from ['TransformationModel']
    """
    
    def __init__(self, data: List[TM_DataPoint] , params: Param ) -> None:
        """
        Cython signature: void TransformationModelLinear(libcpp_vector[TM_DataPoint] & data, Param & params)
        """
        ...
    
    def evaluate(self, value: float ) -> float:
        """
        Cython signature: double evaluate(double value)
        """
        ...
    
    def invert(self) -> None:
        """
        Cython signature: void invert()
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        """
        ...
    
    def weightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void weightData(libcpp_vector[TM_DataPoint] & data)
        Weight the data by the given weight function
        """
        ...
    
    def checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool:
        """
        Cython signature: bool checkValidWeight(const String & weight, libcpp_vector[String] & valid_weights)
        Check for a valid weighting function string
        """
        ...
    
    def weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double weightDatum(double & datum, const String & weight)
        Weight the data according to the weighting function
        """
        ...
    
    def unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double unWeightDatum(double & datum, const String & weight)
        Apply the reverse of the weighting function to the data
        """
        ...
    
    def getValidXWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidXWeights()
        Returns a list of valid x weight function stringss
        """
        ...
    
    def getValidYWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidYWeights()
        Returns a list of valid y weight function strings
        """
        ...
    
    def unWeightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void unWeightData(libcpp_vector[TM_DataPoint] & data)
        Unweight the data by the given weight function
        """
        ...
    
    def checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float:
        """
        Cython signature: double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max)
        Check that the datum is within the valid min and max bounds
        """
        ...
    
    getDefaultParameters: __static_TransformationModelLinear_getDefaultParameters 


class ITRAQ_TYPES:
    None
    FOURPLEX : int
    EIGHTPLEX : int
    TMT_SIXPLEX : int
    SIZE_OF_ITRAQ_TYPES : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __IntensityThresholdCalculation:
    None
    MANUAL : int
    AUTOMAXBYSTDEV : int
    AUTOMAXBYPERCENT : int

    def getMapping(self) -> Dict[int, str]:
       ... 

