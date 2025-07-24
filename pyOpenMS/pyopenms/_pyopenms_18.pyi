from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AbsoluteQuantitationMethod:
    """
    Cython implementation of _AbsoluteQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationMethod.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AbsoluteQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: AbsoluteQuantitationMethod ) -> None:
        """
        Cython signature: void AbsoluteQuantitationMethod(AbsoluteQuantitationMethod &)
        """
        ...
    
    def setLLOD(self, llod: float ) -> None:
        """
        Cython signature: void setLLOD(double llod)
        """
        ...
    
    def setULOD(self, ulod: float ) -> None:
        """
        Cython signature: void setULOD(double ulod)
        """
        ...
    
    def getLLOD(self) -> float:
        """
        Cython signature: double getLLOD()
        """
        ...
    
    def getULOD(self) -> float:
        """
        Cython signature: double getULOD()
        """
        ...
    
    def setLLOQ(self, lloq: float ) -> None:
        """
        Cython signature: void setLLOQ(double lloq)
        """
        ...
    
    def setULOQ(self, uloq: float ) -> None:
        """
        Cython signature: void setULOQ(double uloq)
        """
        ...
    
    def getLLOQ(self) -> float:
        """
        Cython signature: double getLLOQ()
        """
        ...
    
    def getULOQ(self) -> float:
        """
        Cython signature: double getULOQ()
        """
        ...
    
    def checkLOD(self, value: float ) -> bool:
        """
        Cython signature: bool checkLOD(double value)
        """
        ...
    
    def checkLOQ(self, value: float ) -> bool:
        """
        Cython signature: bool checkLOQ(double value)
        """
        ...
    
    def setComponentName(self, component_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComponentName(const String & component_name)
        """
        ...
    
    def setISName(self, IS_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setISName(const String & IS_name)
        """
        ...
    
    def setFeatureName(self, feature_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFeatureName(const String & feature_name)
        """
        ...
    
    def getComponentName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComponentName()
        """
        ...
    
    def getISName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getISName()
        """
        ...
    
    def getFeatureName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFeatureName()
        """
        ...
    
    def setConcentrationUnits(self, concentration_units: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setConcentrationUnits(const String & concentration_units)
        """
        ...
    
    def getConcentrationUnits(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getConcentrationUnits()
        """
        ...
    
    def setTransformationModel(self, transformation_model: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTransformationModel(const String & transformation_model)
        """
        ...
    
    def setTransformationModelParams(self, transformation_model_param: Param ) -> None:
        """
        Cython signature: void setTransformationModelParams(Param transformation_model_param)
        """
        ...
    
    def getTransformationModel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTransformationModel()
        """
        ...
    
    def getTransformationModelParams(self) -> Param:
        """
        Cython signature: Param getTransformationModelParams()
        """
        ...
    
    def setNPoints(self, n_points: int ) -> None:
        """
        Cython signature: void setNPoints(int n_points)
        """
        ...
    
    def setCorrelationCoefficient(self, correlation_coefficient: float ) -> None:
        """
        Cython signature: void setCorrelationCoefficient(double correlation_coefficient)
        """
        ...
    
    def getNPoints(self) -> int:
        """
        Cython signature: int getNPoints()
        """
        ...
    
    def getCorrelationCoefficient(self) -> float:
        """
        Cython signature: double getCorrelationCoefficient()
        """
        ... 


class ClusterProxyKD:
    """
    Cython implementation of _ClusterProxyKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ClusterProxyKD.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ClusterProxyKD()
        """
        ...
    
    @overload
    def __init__(self, in_0: ClusterProxyKD ) -> None:
        """
        Cython signature: void ClusterProxyKD(ClusterProxyKD &)
        """
        ...
    
    @overload
    def __init__(self, size: int , avg_distance: float , center_index: int ) -> None:
        """
        Cython signature: void ClusterProxyKD(size_t size, double avg_distance, size_t center_index)
        """
        ...
    
    def getSize(self) -> int:
        """
        Cython signature: size_t getSize()
        """
        ...
    
    def isValid(self) -> bool:
        """
        Cython signature: bool isValid()
        """
        ...
    
    def getAvgDistance(self) -> float:
        """
        Cython signature: double getAvgDistance()
        """
        ...
    
    def getCenterIndex(self) -> int:
        """
        Cython signature: size_t getCenterIndex()
        """
        ...
    
    def __richcmp__(self, other: ClusterProxyKD, op: int) -> Any:
        ... 


class ConsensusIDAlgorithmWorst:
    """
    Cython implementation of _ConsensusIDAlgorithmWorst

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmWorst.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmWorst()
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


class ContactPerson:
    """
    Cython implementation of _ContactPerson

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ContactPerson.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ContactPerson()
        """
        ...
    
    @overload
    def __init__(self, in_0: ContactPerson ) -> None:
        """
        Cython signature: void ContactPerson(ContactPerson &)
        """
        ...
    
    def getFirstName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFirstName()
        Returns the first name of the person
        """
        ...
    
    def setFirstName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFirstName(String name)
        Sets the first name of the person
        """
        ...
    
    def getLastName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getLastName()
        Returns the last name of the person
        """
        ...
    
    def setLastName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLastName(String name)
        Sets the last name of the person
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the full name of the person (gets split into first and last name internally)
        """
        ...
    
    def getInstitution(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getInstitution()
        Returns the affiliation
        """
        ...
    
    def setInstitution(self, institution: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setInstitution(String institution)
        Sets the affiliation
        """
        ...
    
    def getEmail(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEmail()
        Returns the email address
        """
        ...
    
    def setEmail(self, email: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setEmail(String email)
        Sets the email address
        """
        ...
    
    def getURL(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getURL()
        Returns the URL associated with the contact person (e.g., the institute webpage
        """
        ...
    
    def setURL(self, email: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setURL(String email)
        Sets the URL associated with the contact person (e.g., the institute webpage
        """
        ...
    
    def getAddress(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getAddress()
        Returns the address
        """
        ...
    
    def setAddress(self, email: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setAddress(String email)
        Sets the address
        """
        ...
    
    def getContactInfo(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getContactInfo()
        Returns miscellaneous info about the contact person
        """
        ...
    
    def setContactInfo(self, contact_info: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setContactInfo(String contact_info)
        Sets miscellaneous info about the contact person
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
    
    def __richcmp__(self, other: ContactPerson, op: int) -> Any:
        ... 


class Element:
    """
    Cython implementation of _Element

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Element.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Element()
        """
        ...
    
    @overload
    def __init__(self, in_0: Element ) -> None:
        """
        Cython signature: void Element(Element &)
        """
        ...
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , symbol: Union[bytes, str, String] , atomic_number: int , average_weight: float , mono_weight: float , isotopes: IsotopeDistribution ) -> None:
        """
        Cython signature: void Element(String name, String symbol, unsigned int atomic_number, double average_weight, double mono_weight, IsotopeDistribution isotopes)
        """
        ...
    
    def setAtomicNumber(self, atomic_number: int ) -> None:
        """
        Cython signature: void setAtomicNumber(unsigned int atomic_number)
        Sets unique atomic number
        """
        ...
    
    def getAtomicNumber(self) -> int:
        """
        Cython signature: unsigned int getAtomicNumber()
        Returns the unique atomic number
        """
        ...
    
    def setAverageWeight(self, weight: float ) -> None:
        """
        Cython signature: void setAverageWeight(double weight)
        Sets the average weight of the element
        """
        ...
    
    def getAverageWeight(self) -> float:
        """
        Cython signature: double getAverageWeight()
        Returns the average weight of the element
        """
        ...
    
    def setMonoWeight(self, weight: float ) -> None:
        """
        Cython signature: void setMonoWeight(double weight)
        Sets the mono isotopic weight of the element
        """
        ...
    
    def getMonoWeight(self) -> float:
        """
        Cython signature: double getMonoWeight()
        Returns the mono isotopic weight of the element
        """
        ...
    
    def setIsotopeDistribution(self, isotopes: IsotopeDistribution ) -> None:
        """
        Cython signature: void setIsotopeDistribution(IsotopeDistribution isotopes)
        Sets the isotope distribution of the element
        """
        ...
    
    def getIsotopeDistribution(self) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution getIsotopeDistribution()
        Returns the isotope distribution of the element
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the element
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the element
        """
        ...
    
    def setSymbol(self, symbol: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSymbol(String symbol)
        Sets symbol of the element
        """
        ...
    
    def getSymbol(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSymbol()
        Returns symbol of the element
        """
        ... 


class IonDetector:
    """
    Cython implementation of _IonDetector

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonDetector.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IonDetector()
        Description of a ion detector (part of a MS Instrument)
        """
        ...
    
    @overload
    def __init__(self, in_0: IonDetector ) -> None:
        """
        Cython signature: void IonDetector(IonDetector &)
        """
        ...
    
    def getType(self) -> int:
        """
        Cython signature: Type_IonDetector getType()
        Returns the detector type
        """
        ...
    
    def setType(self, type_: int ) -> None:
        """
        Cython signature: void setType(Type_IonDetector type_)
        Sets the detector type
        """
        ...
    
    def getAcquisitionMode(self) -> int:
        """
        Cython signature: AcquisitionMode getAcquisitionMode()
        Returns the acquisition mode
        """
        ...
    
    def setAcquisitionMode(self, acquisition_mode: int ) -> None:
        """
        Cython signature: void setAcquisitionMode(AcquisitionMode acquisition_mode)
        Sets the acquisition mode
        """
        ...
    
    def getResolution(self) -> float:
        """
        Cython signature: double getResolution()
        Returns the resolution (in ns)
        """
        ...
    
    def setResolution(self, resolution: float ) -> None:
        """
        Cython signature: void setResolution(double resolution)
        Sets the resolution (in ns)
        """
        ...
    
    def getADCSamplingFrequency(self) -> float:
        """
        Cython signature: double getADCSamplingFrequency()
        Returns the analog-to-digital converter sampling frequency (in Hz)
        """
        ...
    
    def setADCSamplingFrequency(self, ADC_sampling_frequency: float ) -> None:
        """
        Cython signature: void setADCSamplingFrequency(double ADC_sampling_frequency)
        Sets the analog-to-digital converter sampling frequency (in Hz)
        """
        ...
    
    def getOrder(self) -> int:
        """
        Cython signature: int getOrder()
        Returns the order
        """
        ...
    
    def setOrder(self, order: int ) -> None:
        """
        Cython signature: void setOrder(int order)
        Sets the order
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
    
    def __richcmp__(self, other: IonDetector, op: int) -> Any:
        ...
    AcquisitionMode : __AcquisitionMode
    Type_IonDetector : __Type_IonDetector 


class MassDecomposition:
    """
    Cython implementation of _MassDecomposition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassDecomposition.html>`_

    Class represents a decomposition of a mass into amino acids
    
    This class represents a mass decomposition into amino acids. A
    decomposition are amino acids given with frequencies which add
    up to a specific mass.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassDecomposition()
        """
        ...
    
    @overload
    def __init__(self, in_0: MassDecomposition ) -> None:
        """
        Cython signature: void MassDecomposition(MassDecomposition &)
        """
        ...
    
    @overload
    def __init__(self, deco: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void MassDecomposition(const String & deco)
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the decomposition as a string
        """
        ...
    
    def toExpandedString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toExpandedString()
        Returns the decomposition as a string; instead of frequencies the amino acids are repeated
        """
        ...
    
    def getNumberOfMaxAA(self) -> int:
        """
        Cython signature: size_t getNumberOfMaxAA()
        Returns the max frequency of this composition
        """
        ...
    
    def containsTag(self, tag: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool containsTag(const String & tag)
        Returns true if tag is contained in the mass decomposition
        """
        ...
    
    def compatible(self, deco: MassDecomposition ) -> bool:
        """
        Cython signature: bool compatible(MassDecomposition & deco)
        Returns true if the mass decomposition if contained in this instance
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the decomposition as a string
        """
        ... 


class ParamNode:
    """
    Cython implementation of _ParamNode

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Param_1_1ParamNode.html>`_
    """
    
    name: Union[bytes, str, String]
    
    description: Union[bytes, str, String]
    
    entries: List[ParamEntry]
    
    nodes: List[ParamNode]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ParamNode()
        """
        ...
    
    @overload
    def __init__(self, in_0: ParamNode ) -> None:
        """
        Cython signature: void ParamNode(ParamNode &)
        """
        ...
    
    @overload
    def __init__(self, n: Union[bytes, str, String] , d: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void ParamNode(const String & n, const String & d)
        """
        ...
    
    def findParentOf(self, name: Union[bytes, str, String] ) -> ParamNode:
        """
        Cython signature: ParamNode * findParentOf(const String & name)
        """
        ...
    
    def findEntryRecursive(self, name: Union[bytes, str, String] ) -> ParamEntry:
        """
        Cython signature: ParamEntry * findEntryRecursive(const String & name)
        """
        ...
    
    @overload
    def insert(self, node: ParamNode , prefix: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void insert(ParamNode & node, const String & prefix)
        """
        ...
    
    @overload
    def insert(self, entry: ParamEntry , prefix: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void insert(ParamEntry & entry, const String & prefix)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def suffix(self, key: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String suffix(const String & key)
        """
        ...
    
    def __richcmp__(self, other: ParamNode, op: int) -> Any:
        ... 


class Peak1D:
    """
    Cython implementation of _Peak1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Peak1D.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Peak1D()
        """
        ...
    
    @overload
    def __init__(self, in_0: Peak1D ) -> None:
        """
        Cython signature: void Peak1D(Peak1D &)
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
    
    def __richcmp__(self, other: Peak1D, op: int) -> Any:
        ... 


class RangeIntensity:
    """
    Cython implementation of _RangeIntensity

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeIntensity.html>`_
      -- Inherits from ['RangeBase']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RangeIntensity()
        """
        ...
    
    @overload
    def __init__(self, in_0: RangeIntensity ) -> None:
        """
        Cython signature: void RangeIntensity(RangeIntensity &)
        """
        ...
    
    def setMinIntensity(self, min: float ) -> None:
        """
        Cython signature: void setMinIntensity(double min)
        """
        ...
    
    def setMaxIntensity(self, max: float ) -> None:
        """
        Cython signature: void setMaxIntensity(double max)
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
    
    def extendIntensity(self, value: float ) -> None:
        """
        Cython signature: void extendIntensity(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def containsIntensity(self, value: float ) -> bool:
        """
        Cython signature: bool containsIntensity(double value)
        Is value within [min, max]?
        """
        ...
    
    def setMin(self, min: float ) -> None:
        """
        Cython signature: void setMin(double min)
        """
        ...
    
    def setMax(self, max: float ) -> None:
        """
        Cython signature: void setMax(double max)
        """
        ...
    
    def getMin(self) -> float:
        """
        Cython signature: double getMin()
        """
        ...
    
    def getMax(self) -> float:
        """
        Cython signature: double getMax()
        """
        ...
    
    def extend(self, value: float ) -> None:
        """
        Cython signature: void extend(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def contains(self, value: float ) -> bool:
        """
        Cython signature: bool contains(double value)
        Is value within [min, max]?
        """
        ... 


class RangeMZ:
    """
    Cython implementation of _RangeMZ

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeMZ.html>`_
      -- Inherits from ['RangeBase']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RangeMZ()
        """
        ...
    
    @overload
    def __init__(self, in_0: RangeMZ ) -> None:
        """
        Cython signature: void RangeMZ(RangeMZ &)
        """
        ...
    
    def setMinMZ(self, min: float ) -> None:
        """
        Cython signature: void setMinMZ(double min)
        """
        ...
    
    def setMaxMZ(self, max: float ) -> None:
        """
        Cython signature: void setMaxMZ(double max)
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
    
    def extendMZ(self, value: float ) -> None:
        """
        Cython signature: void extendMZ(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def containsMZ(self, value: float ) -> bool:
        """
        Cython signature: bool containsMZ(double value)
        Is value within [min, max]?
        """
        ...
    
    def setMin(self, min: float ) -> None:
        """
        Cython signature: void setMin(double min)
        """
        ...
    
    def setMax(self, max: float ) -> None:
        """
        Cython signature: void setMax(double max)
        """
        ...
    
    def getMin(self) -> float:
        """
        Cython signature: double getMin()
        """
        ...
    
    def getMax(self) -> float:
        """
        Cython signature: double getMax()
        """
        ...
    
    def extend(self, value: float ) -> None:
        """
        Cython signature: void extend(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def contains(self, value: float ) -> bool:
        """
        Cython signature: bool contains(double value)
        Is value within [min, max]?
        """
        ... 


class RangeMobility:
    """
    Cython implementation of _RangeMobility

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeMobility.html>`_
      -- Inherits from ['RangeBase']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RangeMobility()
        """
        ...
    
    @overload
    def __init__(self, in_0: RangeMobility ) -> None:
        """
        Cython signature: void RangeMobility(RangeMobility &)
        """
        ...
    
    def setMinMobility(self, min: float ) -> None:
        """
        Cython signature: void setMinMobility(double min)
        """
        ...
    
    def setMaxMobility(self, max: float ) -> None:
        """
        Cython signature: void setMaxMobility(double max)
        """
        ...
    
    def getMinMobility(self) -> float:
        """
        Cython signature: double getMinMobility()
        Returns the minimum Mobility
        """
        ...
    
    def getMaxMobility(self) -> float:
        """
        Cython signature: double getMaxMobility()
        Returns the maximum Mobility
        """
        ...
    
    def extendMobility(self, value: float ) -> None:
        """
        Cython signature: void extendMobility(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def containsMobility(self, value: float ) -> bool:
        """
        Cython signature: bool containsMobility(double value)
        Is value within [min, max]?
        """
        ...
    
    def setMin(self, min: float ) -> None:
        """
        Cython signature: void setMin(double min)
        """
        ...
    
    def setMax(self, max: float ) -> None:
        """
        Cython signature: void setMax(double max)
        """
        ...
    
    def getMin(self) -> float:
        """
        Cython signature: double getMin()
        """
        ...
    
    def getMax(self) -> float:
        """
        Cython signature: double getMax()
        """
        ...
    
    def extend(self, value: float ) -> None:
        """
        Cython signature: void extend(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def contains(self, value: float ) -> bool:
        """
        Cython signature: bool contains(double value)
        Is value within [min, max]?
        """
        ... 


class RangeRT:
    """
    Cython implementation of _RangeRT

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeRT.html>`_
      -- Inherits from ['RangeBase']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RangeRT()
        """
        ...
    
    @overload
    def __init__(self, in_0: RangeRT ) -> None:
        """
        Cython signature: void RangeRT(RangeRT &)
        """
        ...
    
    def setMinRT(self, min: float ) -> None:
        """
        Cython signature: void setMinRT(double min)
        """
        ...
    
    def setMaxRT(self, max: float ) -> None:
        """
        Cython signature: void setMaxRT(double max)
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
    
    def extendRT(self, value: float ) -> None:
        """
        Cython signature: void extendRT(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def containsRT(self, value: float ) -> bool:
        """
        Cython signature: bool containsRT(double value)
        Is value within [min, max]?
        """
        ...
    
    def setMin(self, min: float ) -> None:
        """
        Cython signature: void setMin(double min)
        """
        ...
    
    def setMax(self, max: float ) -> None:
        """
        Cython signature: void setMax(double max)
        """
        ...
    
    def getMin(self) -> float:
        """
        Cython signature: double getMin()
        """
        ...
    
    def getMax(self) -> float:
        """
        Cython signature: double getMax()
        """
        ...
    
    def extend(self, value: float ) -> None:
        """
        Cython signature: void extend(double value)
        Extend the range such that it includes the given @p value
        """
        ...
    
    def contains(self, value: float ) -> bool:
        """
        Cython signature: bool contains(double value)
        Is value within [min, max]?
        """
        ... 


class SpectraMerger:
    """
    Cython implementation of _SpectraMerger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectraMerger.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectraMerger()
        Merges blocks of MS or MS2 spectra
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectraMerger ) -> None:
        """
        Cython signature: void SpectraMerger(SpectraMerger &)
        """
        ...
    
    def mergeSpectraBlockWise(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void mergeSpectraBlockWise(MSExperiment & exp)
        """
        ...
    
    def mergeSpectraPrecursors(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void mergeSpectraPrecursors(MSExperiment & exp)
        Merges spectra with similar precursors (must have MS2 level)
        """
        ...
    
    def average(self, exp: MSExperiment , average_type: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void average(MSExperiment & exp, String average_type)
        Average over neighbouring spectra
        
        :param exp: Experimental data to be averaged
        :param average_type: Averaging type to be used ("gaussian" or "tophat")
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


class SwathWindowLoader:
    """
    Cython implementation of _SwathWindowLoader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathWindowLoader.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SwathWindowLoader()
        """
        ...
    
    @overload
    def __init__(self, in_0: SwathWindowLoader ) -> None:
        """
        Cython signature: void SwathWindowLoader(SwathWindowLoader &)
        """
        ...
    
    def annotateSwathMapsFromFile(self, filename: Union[bytes, str, String] , swath_maps: List[SwathMap] , do_sort: bool , force: bool ) -> None:
        """
        Cython signature: void annotateSwathMapsFromFile(String filename, libcpp_vector[SwathMap] & swath_maps, bool do_sort, bool force)
        """
        ...
    
    def readSwathWindows(self, filename: Union[bytes, str, String] , swath_prec_lower: List[float] , swath_prec_upper: List[float] ) -> None:
        """
        Cython signature: void readSwathWindows(String filename, libcpp_vector[double] & swath_prec_lower, libcpp_vector[double] & swath_prec_upper)
        """
        ... 


class TransitionPQPFile:
    """
    Cython implementation of _TransitionPQPFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransitionPQPFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TransitionPQPFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: TransitionPQPFile ) -> None:
        """
        Cython signature: void TransitionPQPFile(TransitionPQPFile &)
        """
        ...
    
    def convertTargetedExperimentToPQP(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void convertTargetedExperimentToPQP(char * filename, TargetedExperiment & targeted_exp)
        Write out a targeted experiment (TraML structure) into a PQP file
        
        :param filename: The output file
        :param targeted_exp: The targeted experiment
        """
        ...
    
    @overload
    def convertPQPToTargetedExperiment(self, filename: bytes , targeted_exp: TargetedExperiment , legacy_traml_id: bool ) -> None:
        """
        Cython signature: void convertPQPToTargetedExperiment(char * filename, TargetedExperiment & targeted_exp, bool legacy_traml_id)
        Read in a PQP file and construct a targeted experiment (TraML structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
        """
        ...
    
    @overload
    def convertPQPToTargetedExperiment(self, filename: bytes , targeted_exp: LightTargetedExperiment , legacy_traml_id: bool ) -> None:
        """
        Cython signature: void convertPQPToTargetedExperiment(char * filename, LightTargetedExperiment & targeted_exp, bool legacy_traml_id)
        Read in a PQP file and construct a targeted experiment (Light transition structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
        """
        ...
    
    def convertTargetedExperimentToTSV(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void convertTargetedExperimentToTSV(char * filename, TargetedExperiment & targeted_exp)
        """
        ...
    
    @overload
    def convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment & targeted_exp)
        """
        ...
    
    @overload
    def convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None:
        """
        Cython signature: void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment & targeted_exp)
        """
        ...
    
    def validateTargetedExperiment(self, targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void validateTargetedExperiment(TargetedExperiment targeted_exp)
        """
        ... 


class __AcquisitionMode:
    None
    ACQMODENULL : int
    PULSECOUNTING : int
    ADC : int
    TDC : int
    TRANSIENTRECORDER : int
    SIZE_OF_ACQUISITIONMODE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __Type_IonDetector:
    None
    TYPENULL : int
    ELECTRONMULTIPLIER : int
    PHOTOMULTIPLIER : int
    FOCALPLANEARRAY : int
    FARADAYCUP : int
    CONVERSIONDYNODEELECTRONMULTIPLIER : int
    CONVERSIONDYNODEPHOTOMULTIPLIER : int
    MULTICOLLECTOR : int
    CHANNELELECTRONMULTIPLIER : int
    CHANNELTRON : int
    DALYDETECTOR : int
    MICROCHANNELPLATEDETECTOR : int
    ARRAYDETECTOR : int
    CONVERSIONDYNODE : int
    DYNODE : int
    FOCALPLANECOLLECTOR : int
    IONTOPHOTONDETECTOR : int
    POINTCOLLECTOR : int
    POSTACCELERATIONDETECTOR : int
    PHOTODIODEARRAYDETECTOR : int
    INDUCTIVEDETECTOR : int
    ELECTRONMULTIPLIERTUBE : int
    SIZE_OF_TYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

