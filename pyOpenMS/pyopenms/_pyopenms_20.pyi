from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_ExperimentalDesignFile_load(tsv_file: Union[bytes, str, String] , in_1: bool ) -> ExperimentalDesign:
    """
    Cython signature: ExperimentalDesign load(const String & tsv_file, bool)
    """
    ...


class AASeqWithMass:
    """
    Cython implementation of _AASeqWithMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AASeqWithMass.html>`_
    """
    
    peptide_mass: float
    
    peptide_seq: AASequence
    
    position: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AASeqWithMass()
        """
        ...
    
    @overload
    def __init__(self, in_0: AASeqWithMass ) -> None:
        """
        Cython signature: void AASeqWithMass(AASeqWithMass &)
        """
        ... 


class Acquisition:
    """
    Cython implementation of _Acquisition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Acquisition.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Acquisition()
        """
        ...
    
    @overload
    def __init__(self, in_0: Acquisition ) -> None:
        """
        Cython signature: void Acquisition(Acquisition &)
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        """
        ...
    
    def setIdentifier(self, identifier: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(const String & identifier)
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
    
    def __richcmp__(self, other: Acquisition, op: int) -> Any:
        ... 


class CV:
    """
    Cython implementation of _CV

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1CV.html>`_
    """
    
    id: Union[bytes, str, String]
    
    fullname: Union[bytes, str, String]
    
    version: Union[bytes, str, String]
    
    URI: Union[bytes, str, String]
    
    @overload
    def __init__(self, in_0: CV ) -> None:
        """
        Cython signature: void CV(CV &)
        """
        ...
    
    @overload
    def __init__(self, new_id: Union[bytes, str, String] , new_fullname: Union[bytes, str, String] , new_version: Union[bytes, str, String] , new_URI: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void CV(String new_id, String new_fullname, String new_version, String new_URI)
        """
        ...
    
    def __richcmp__(self, other: CV, op: int) -> Any:
        ... 


class Compound:
    """
    Cython implementation of _Compound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Compound.html>`_
      -- Inherits from ['CVTermList']
    """
    
    id: Union[bytes, str, String]
    
    molecular_formula: Union[bytes, str, String]
    
    smiles_string: Union[bytes, str, String]
    
    theoretical_mass: float
    
    rts: List[RetentionTime]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Compound()
        """
        ...
    
    @overload
    def __init__(self, in_0: Compound ) -> None:
        """
        Cython signature: void Compound(Compound &)
        """
        ...
    
    def setChargeState(self, charge: int ) -> None:
        """
        Cython signature: void setChargeState(int charge)
        Sets the peptide or compound charge state
        """
        ...
    
    def getChargeState(self) -> int:
        """
        Cython signature: int getChargeState()
        Returns the peptide or compound charge state
        """
        ...
    
    def hasCharge(self) -> bool:
        """
        Cython signature: bool hasCharge()
        Whether peptide or compound has set charge state
        """
        ...
    
    def getRetentionTime(self) -> float:
        """
        Cython signature: double getRetentionTime()
        Gets compound or peptide retention time
        """
        ...
    
    def hasRetentionTime(self) -> bool:
        """
        Cython signature: bool hasRetentionTime()
        Check whether compound or peptide has an annotated retention time
        """
        ...
    
    def getRetentionTimeType(self) -> int:
        """
        Cython signature: RTType getRetentionTimeType()
        Get compound or peptide retentiontime type
        """
        ...
    
    def getRetentionTimeUnit(self) -> int:
        """
        Cython signature: RTUnit getRetentionTimeUnit()
        Get compound or peptide retentiontime type
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
    
    def __richcmp__(self, other: Compound, op: int) -> Any:
        ... 


class Configuration:
    """
    Cython implementation of _Configuration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Configuration.html>`_
      -- Inherits from ['CVTermList']
    """
    
    contact_ref: Union[bytes, str, String]
    
    instrument_ref: Union[bytes, str, String]
    
    validations: List[CVTermList]
    
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
    
    def __richcmp__(self, other: Configuration, op: int) -> Any:
        ... 


class Contact:
    """
    Cython implementation of _Contact

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Contact.html>`_
      -- Inherits from ['CVTermList']
    """
    
    id: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Contact()
        """
        ...
    
    @overload
    def __init__(self, in_0: Contact ) -> None:
        """
        Cython signature: void Contact(Contact &)
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
    
    def __richcmp__(self, other: Contact, op: int) -> Any:
        ... 


class CrossLinksDB:
    """
    Cython implementation of _CrossLinksDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CrossLinksDB.html>`_
    """
    
    def getNumberOfModifications(self) -> int:
        """
        Cython signature: size_t getNumberOfModifications()
        """
        ...
    
    def searchModifications(self, mods: Set[ResidueModification] , mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> None:
        """
        Cython signature: void searchModifications(libcpp_set[const ResidueModification *] & mods, const String & mod_name, const String & residue, TermSpecificity term_spec)
        """
        ...
    
    @overload
    def getModification(self, index: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(size_t index)
        """
        ...
    
    @overload
    def getModification(self, mod_name: Union[bytes, str, String] ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(const String & mod_name)
        """
        ...
    
    @overload
    def getModification(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(const String & mod_name, const String & residue, TermSpecificity term_spec)
        """
        ...
    
    def has(self, modification: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool has(String modification)
        """
        ...
    
    def findModificationIndex(self, mod_name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t findModificationIndex(const String & mod_name)
        """
        ...
    
    def searchModificationsByDiffMonoMass(self, mods: List[bytes] , mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> None:
        """
        Cython signature: void searchModificationsByDiffMonoMass(libcpp_vector[String] & mods, double mass, double max_error, const String & residue, TermSpecificity term_spec)
        """
        ...
    
    def getBestModificationByDiffMonoMass(self, mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getBestModificationByDiffMonoMass(double mass, double max_error, const String residue, TermSpecificity term_spec)
        """
        ...
    
    def getAllSearchModifications(self, modifications: List[bytes] ) -> None:
        """
        Cython signature: void getAllSearchModifications(libcpp_vector[String] & modifications)
        """
        ...
    
    def readFromOBOFile(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readFromOBOFile(const String & filename)
        """
        ...
    
    def isInstantiated(self) -> bool:
        """
        Cython signature: bool isInstantiated()
        """
        ... 


class ExperimentalDesignFile:
    """
    Cython implementation of _ExperimentalDesignFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesignFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ExperimentalDesignFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: ExperimentalDesignFile ) -> None:
        """
        Cython signature: void ExperimentalDesignFile(ExperimentalDesignFile &)
        """
        ...
    
    load: __static_ExperimentalDesignFile_load 


class MRMTransitionGroupPicker:
    """
    Cython implementation of _MRMTransitionGroupPicker

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroupPicker.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMTransitionGroupPicker()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMTransitionGroupPicker ) -> None:
        """
        Cython signature: void MRMTransitionGroupPicker(MRMTransitionGroupPicker &)
        """
        ...
    
    @overload
    def pickTransitionGroup(self, transition_group: LightMRMTransitionGroupCP ) -> None:
        """
        Cython signature: void pickTransitionGroup(LightMRMTransitionGroupCP transition_group)
        """
        ...
    
    @overload
    def pickTransitionGroup(self, transition_group: MRMTransitionGroupCP ) -> None:
        """
        Cython signature: void pickTransitionGroup(MRMTransitionGroupCP transition_group)
        """
        ...
    
    def createMRMFeature(self, transition_group: LightMRMTransitionGroupCP , picked_chroms: List[MSChromatogram] , smoothed_chroms: List[MSChromatogram] , chr_idx: int , peak_idx: int ) -> MRMFeature:
        """
        Cython signature: MRMFeature createMRMFeature(LightMRMTransitionGroupCP transition_group, libcpp_vector[MSChromatogram] & picked_chroms, libcpp_vector[MSChromatogram] & smoothed_chroms, const int chr_idx, const int peak_idx)
        """
        ...
    
    def remove_overlapping_features(self, picked_chroms: List[MSChromatogram] , best_left: float , best_right: float ) -> None:
        """
        Cython signature: void remove_overlapping_features(libcpp_vector[MSChromatogram] & picked_chroms, double best_left, double best_right)
        """
        ...
    
    def findLargestPeak(self, picked_chroms: List[MSChromatogram] , chr_idx: int , peak_idx: int ) -> None:
        """
        Cython signature: void findLargestPeak(libcpp_vector[MSChromatogram] & picked_chroms, int & chr_idx, int & peak_idx)
        """
        ...
    
    def findWidestPeakIndices(self, picked_chroms: List[MSChromatogram] , chrom_idx: int , point_idx: int ) -> None:
        """
        Cython signature: void findWidestPeakIndices(libcpp_vector[MSChromatogram] & picked_chroms, int & chrom_idx, int & point_idx)
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


class MSSpectrum:
    """
    Cython implementation of _MSSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSSpectrum.html>`_
      -- Inherits from ['SpectrumSettings', 'RangeManagerMzInt']

    The representation of a 1D spectrum.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    Iterations yields access to underlying peak objects but is slower
    Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
    See help(SpectrumSettings) for information about meta-information
    
    Usage:
    
    .. code-block:: python
    
      ms_level = spectrum.getMSLevel()
      rt = spectrum.getRT()
      mz, intensities = spectrum.get_peaks()
    
    
    Usage:
    
    .. code-block:: python
    
      from pyopenms import *
    
      spectrum = MSSpectrum()
      spectrum.setDriftTime(25) # 25 ms
      spectrum.setRT(205.2) # 205.2 s
      spectrum.setMSLevel(3) # MS3
      p = Precursor()
      p.setIsolationWindowLowerOffset(1.5)
      p.setIsolationWindowUpperOffset(1.5)
      p.setMZ(600) # isolation at 600 +/- 1.5 Th
      p.setActivationEnergy(40) # 40 eV
      p.setCharge(4) # 4+ ion
      spectrum.setPrecursors( [p] )
    
      # Add raw data to spectrum
      spectrum.set_peaks( ([401.5], [900]) )
    
      # Additional data arrays / peak annotations
      fda = FloatDataArray()
      fda.setName("Signal to Noise Array")
      fda.push_back(15)
      sda = StringDataArray()
      sda.setName("Peak annotation")
      sda.push_back("y15++")
      spectrum.setFloatDataArrays( [fda] )
      spectrum.setStringDataArrays( [sda] )
    
      # Add spectrum to MSExperiment
      exp = MSExperiment()
      exp.addSpectrum(spectrum)
    
      # Add second spectrum and store as mzML file
      spectrum2 = MSSpectrum()
      spectrum2.set_peaks( ([1, 2], [1, 2]) )
      exp.addSpectrum(spectrum2)
    
      MzMLFile().store("testfile.mzML", exp)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSSpectrum()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSSpectrum ) -> None:
        """
        Cython signature: void MSSpectrum(MSSpectrum &)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        Returns the absolute retention time (in seconds)
        """
        ...
    
    def setRT(self, in_0: float ) -> None:
        """
        Cython signature: void setRT(double)
        Sets the absolute retention time (in seconds)
        """
        ...
    
    def getDriftTime(self) -> float:
        """
        Cython signature: double getDriftTime()
        Returns the drift time (-1 if not set)
        """
        ...
    
    def setDriftTime(self, in_0: float ) -> None:
        """
        Cython signature: void setDriftTime(double)
        Sets the drift time (-1 if not set)
        """
        ...
    
    def getDriftTimeUnit(self) -> int:
        """
        Cython signature: DriftTimeUnit getDriftTimeUnit()
        """
        ...
    
    def getDriftTimeUnitAsString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDriftTimeUnitAsString()
        """
        ...
    
    def setDriftTimeUnit(self, dt: int ) -> None:
        """
        Cython signature: void setDriftTimeUnit(DriftTimeUnit dt)
        """
        ...
    
    def containsIMData(self) -> bool:
        """
        Cython signature: bool containsIMData()
        """
        ...
    
    def getMSLevel(self) -> int:
        """
        Cython signature: unsigned int getMSLevel()
        Returns the MS level
        """
        ...
    
    def setMSLevel(self, in_0: int ) -> None:
        """
        Cython signature: void setMSLevel(unsigned int)
        Sets the MS level
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Returns the number of peaks in the spectrum
        """
        ...
    
    def reserve(self, n: int ) -> None:
        """
        Cython signature: void reserve(size_t n)
        """
        ...
    
    def resize(self, n: int ) -> None:
        """
        Cython signature: void resize(size_t n)
        Resize the peak array
        """
        ...
    
    def __getitem__(self, in_0: int ) -> Peak1D:
        """
        Cython signature: Peak1D & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: Peak1D ) -> None:
        """Cython signature: Peak1D & operator[](size_t)"""
        ...
    
    def updateRanges(self) -> None:
        """
        Cython signature: void updateRanges()
        """
        ...
    
    def clear(self, clear_meta_data: bool ) -> None:
        """
        Cython signature: void clear(bool clear_meta_data)
        Clears all data (and meta data if clear_meta_data is true)
        """
        ...
    
    def push_back(self, in_0: Peak1D ) -> None:
        """
        Cython signature: void push_back(Peak1D)
        Append a peak
        """
        ...
    
    def isSorted(self) -> bool:
        """
        Cython signature: bool isSorted()
        Returns true if the spectrum is sorte by m/z
        """
        ...
    
    @overload
    def findNearest(self, mz: float ) -> int:
        """
        Cython signature: int findNearest(double mz)
        Returns the index of the closest peak in m/z
        """
        ...
    
    @overload
    def findNearest(self, mz: float , tolerance: float ) -> int:
        """
        Cython signature: int findNearest(double mz, double tolerance)
        Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)
        """
        ...
    
    @overload
    def findNearest(self, mz: float , tolerance_left: float , tolerance_right: float ) -> int:
        """
        Cython signature: int findNearest(double mz, double tolerance_left, double tolerance_right)
        Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        """
        ...
    
    def findHighestInWindow(self, mz: float , tolerance_left: float , tolerance_right: float ) -> int:
        """
        Cython signature: int findHighestInWindow(double mz, double tolerance_left, double tolerance_right)
        Returns the index of the highest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        """
        ...
    
    def select(self, indices: List[int] ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum select(libcpp_vector[size_t] & indices)
        Subset the spectrum by indices. Also applies to associated data arrays if present.
        """
        ...
    
    def calculateTIC(self) -> float:
        """
        Cython signature: double calculateTIC()
        Returns the total ion current (=sum) of peak intensities in the spectrum
        """
        ...
    
    def sortByIntensity(self, reverse: bool ) -> None:
        """
        Cython signature: void sortByIntensity(bool reverse)
        """
        ...
    
    def sortByPosition(self) -> None:
        """
        Cython signature: void sortByPosition()
        """
        ...
    
    def getFloatDataArrays(self) -> List[FloatDataArray]:
        """
        Cython signature: libcpp_vector[FloatDataArray] getFloatDataArrays()
        Returns the additional float data arrays to store e.g. meta data
        """
        ...
    
    def getIntegerDataArrays(self) -> List[IntegerDataArray]:
        """
        Cython signature: libcpp_vector[IntegerDataArray] getIntegerDataArrays()
        Returns the additional int data arrays to store e.g. meta data
        """
        ...
    
    def getStringDataArrays(self) -> List[StringDataArray]:
        """
        Cython signature: libcpp_vector[StringDataArray] getStringDataArrays()
        Returns the additional string data arrays to store e.g. meta data
        """
        ...
    
    def setFloatDataArrays(self, fda: List[FloatDataArray] ) -> None:
        """
        Cython signature: void setFloatDataArrays(libcpp_vector[FloatDataArray] fda)
        Sets the additional float data arrays to store e.g. meta data
        """
        ...
    
    def setIntegerDataArrays(self, ida: List[IntegerDataArray] ) -> None:
        """
        Cython signature: void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida)
        Sets the additional int data arrays to store e.g. meta data
        """
        ...
    
    def setStringDataArrays(self, sda: List[StringDataArray] ) -> None:
        """
        Cython signature: void setStringDataArrays(libcpp_vector[StringDataArray] sda)
        Sets the additional string data arrays to store e.g. meta data
        """
        ...
    
    def unify(self, in_0: SpectrumSettings ) -> None:
        """
        Cython signature: void unify(SpectrumSettings)
        """
        ...
    
    def getType(self) -> int:
        """
        Cython signature: int getType()
        Returns the spectrum type (centroided (PEAKS) or profile data (RAW))
        """
        ...
    
    def setType(self, in_0: int ) -> None:
        """
        Cython signature: void setType(SpectrumType)
        Sets the spectrum type
        """
        ...
    
    def getNativeID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeID()
        Returns the native identifier for the spectrum, used by the acquisition software
        """
        ...
    
    def setNativeID(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeID(String)
        Sets the native identifier for the spectrum, used by the acquisition software
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the free-text comment
        """
        ...
    
    def setComment(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String)
        Sets the free-text comment
        """
        ...
    
    def getInstrumentSettings(self) -> InstrumentSettings:
        """
        Cython signature: InstrumentSettings getInstrumentSettings()
        Returns a const reference to the instrument settings of the current spectrum
        """
        ...
    
    def setInstrumentSettings(self, in_0: InstrumentSettings ) -> None:
        """
        Cython signature: void setInstrumentSettings(InstrumentSettings)
        Sets the instrument settings of the current spectrum
        """
        ...
    
    def getAcquisitionInfo(self) -> AcquisitionInfo:
        """
        Cython signature: AcquisitionInfo getAcquisitionInfo()
        Returns a const reference to the acquisition info
        """
        ...
    
    def setAcquisitionInfo(self, in_0: AcquisitionInfo ) -> None:
        """
        Cython signature: void setAcquisitionInfo(AcquisitionInfo)
        Sets the acquisition info
        """
        ...
    
    def getSourceFile(self) -> SourceFile:
        """
        Cython signature: SourceFile getSourceFile()
        Returns a const reference to the source file
        """
        ...
    
    def setSourceFile(self, in_0: SourceFile ) -> None:
        """
        Cython signature: void setSourceFile(SourceFile)
        Sets the source file
        """
        ...
    
    def getPrecursors(self) -> List[Precursor]:
        """
        Cython signature: libcpp_vector[Precursor] getPrecursors()
        Returns a const reference to the precursors
        """
        ...
    
    def setPrecursors(self, in_0: List[Precursor] ) -> None:
        """
        Cython signature: void setPrecursors(libcpp_vector[Precursor])
        Sets the precursors
        """
        ...
    
    def getProducts(self) -> List[Product]:
        """
        Cython signature: libcpp_vector[Product] getProducts()
        Returns a const reference to the products
        """
        ...
    
    def setProducts(self, in_0: List[Product] ) -> None:
        """
        Cython signature: void setProducts(libcpp_vector[Product])
        Sets the products
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[shared_ptr[DataProcessing]] getDataProcessing()
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[shared_ptr[DataProcessing]])
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
    
    def __richcmp__(self, other: MSSpectrum, op: int) -> Any:
        ...
    
    def __iter__(self) -> Peak1D:
       ... 


class MetaInfo:
    """
    Cython implementation of _MetaInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfo.html>`_

    A Type-Name-Value tuple class
    
    MetaInfo maps an index (an integer corresponding to a string) to
    DataValue objects.  The mapping of strings to the index is performed by
    the MetaInfoRegistry, which can be accessed by the method registry()
    
    There are two versions of nearly all members. One which operates with a
    string name and another one which operates on an index. The index version
    is always faster, as it does not need to look up the index corresponding
    to the string in the MetaInfoRegistry
    
    If you wish to add a MetaInfo member to a class, consider deriving that
    class from MetaInfoInterface, instead of simply adding MetaInfo as
    member. MetaInfoInterface implements a full interface to a MetaInfo
    member and is more memory efficient if no meta info gets added
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaInfo ) -> None:
        """
        Cython signature: void MetaInfo(MetaInfo &)
        """
        ...
    
    @overload
    def getValue(self, name: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getValue(String name)
        Returns the value corresponding to a string
        """
        ...
    
    @overload
    def getValue(self, index: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getValue(unsigned int index)
        Returns the value corresponding to an index
        """
        ...
    
    @overload
    def getValue(self, name: Union[bytes, str, String] , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getValue(String name, DataValue default_value)
        Returns the value corresponding to a string
        """
        ...
    
    @overload
    def getValue(self, index: int , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getValue(unsigned int index, DataValue default_value)
        Returns the value corresponding to an index
        """
        ...
    
    @overload
    def exists(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool exists(String name)
        Returns if this MetaInfo is set
        """
        ...
    
    @overload
    def exists(self, index: int ) -> bool:
        """
        Cython signature: bool exists(unsigned int index)
        Returns if this MetaInfo is set
        """
        ...
    
    @overload
    def setValue(self, name: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setValue(String name, DataValue value)
        Sets the DataValue corresponding to a name
        """
        ...
    
    @overload
    def setValue(self, index: int , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setValue(unsigned int index, DataValue value)
        Sets the DataValue corresponding to an index
        """
        ...
    
    @overload
    def removeValue(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeValue(String name)
        Removes the DataValue corresponding to `name` if it exists
        """
        ...
    
    @overload
    def removeValue(self, index: int ) -> None:
        """
        Cython signature: void removeValue(unsigned int index)
        Removes the DataValue corresponding to `index` if it exists
        """
        ...
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        Fills the given vector with a list of all keys for which a value is set
        """
        ...
    
    def getKeysAsIntegers(self, keys: List[int] ) -> None:
        """
        Cython signature: void getKeysAsIntegers(libcpp_vector[unsigned int] & keys)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Returns if the MetaInfo is empty
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Removes all meta values
        """
        ...
    
    def registry(self) -> MetaInfoRegistry:
        """
        Cython signature: MetaInfoRegistry registry()
        """
        ... 


class NLargest:
    """
    Cython implementation of _NLargest

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NLargest.html>`_
      -- Inherits from ['DefaultParamHandler']

    NLargest removes all but the n largest peaks
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NLargest()
        """
        ...
    
    @overload
    def __init__(self, in_0: NLargest ) -> None:
        """
        Cython signature: void NLargest(NLargest &)
        """
        ...
    
    def filterSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterSpectrum(MSSpectrum & spec)
        Keep only n-largest peaks in spectrum
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        Keep only n-largest peaks in spectrum
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
        Keep only n-largest peaks in each spectrum of a peak map
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


class Nuctide:
    """
    Cython implementation of _Nuctide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Nuctide.html>`_
      -- Inherits from ['CVTermList']
    """
    
    rts: List[RetentionTime]
    
    id: Union[bytes, str, String]
    
    oligo_refs: List[bytes]
    
    evidence: CVTermList
    
    sequence: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Nuctide()
        """
        ...
    
    @overload
    def __init__(self, in_0: Nuctide ) -> None:
        """
        Cython signature: void Nuctide(Nuctide &)
        """
        ...
    
    def setNuctideGroupLabel(self, label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNuctideGroupLabel(String label)
        Sets the nuctide group label
        """
        ...
    
    def getNuctideGroupLabel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNuctideGroupLabel()
        Get the nuctide group label
        """
        ...
    
    def setChargeState(self, charge: int ) -> None:
        """
        Cython signature: void setChargeState(int charge)
        Sets the nuctide or compound charge states
        """
        ...
    
    def getChargeState(self) -> int:
        """
        Cython signature: int getChargeState()
        Returns the nuctide charge state
        """
        ...
    
    def hasCharge(self) -> bool:
        """
        Cython signature: bool hasCharge()
        Whether product has set charge state
        """
        ...
    
    def getRetentionTime(self) -> float:
        """
        Cython signature: double getRetentionTime()
        Gets nuctide retention time
        """
        ...
    
    def hasRetentionTime(self) -> bool:
        """
        Cython signature: bool hasRetentionTime()
        Gets nuctide retention time
        """
        ...
    
    def getRetentionTimeType(self) -> int:
        """
        Cython signature: RTType getRetentionTimeType()
        Get nuctide retentiontime type
        """
        ...
    
    def getRetentionTimeUnit(self) -> int:
        """
        Cython signature: RTUnit getRetentionTimeUnit()
        Get nuctide retentiontime unit (minute/seconds)
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
    
    def __richcmp__(self, other: Nuctide, op: int) -> Any:
        ... 


class Peptide:
    """
    Cython implementation of _Peptide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Peptide.html>`_
      -- Inherits from ['CVTermList']
    """
    
    rts: List[RetentionTime]
    
    id: Union[bytes, str, String]
    
    protein_refs: List[bytes]
    
    evidence: CVTermList
    
    sequence: Union[bytes, str, String]
    
    mods: List[TargetedExperiment_Modification]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Peptide()
        """
        ...
    
    @overload
    def __init__(self, in_0: Peptide ) -> None:
        """
        Cython signature: void Peptide(Peptide &)
        """
        ...
    
    def setPeptideGroupLabel(self, label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPeptideGroupLabel(String label)
        Sets the peptide group label
        """
        ...
    
    def getPeptideGroupLabel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPeptideGroupLabel()
        Get the peptide group label
        """
        ...
    
    def setChargeState(self, charge: int ) -> None:
        """
        Cython signature: void setChargeState(int charge)
        Sets the peptide or compound charge states
        """
        ...
    
    def getChargeState(self) -> int:
        """
        Cython signature: int getChargeState()
        Returns the peptide or compound charge state
        """
        ...
    
    def hasCharge(self) -> bool:
        """
        Cython signature: bool hasCharge()
        Whether product has set charge state
        """
        ...
    
    def getRetentionTime(self) -> float:
        """
        Cython signature: double getRetentionTime()
        Gets compound or peptide retention time
        """
        ...
    
    def hasRetentionTime(self) -> bool:
        """
        Cython signature: bool hasRetentionTime()
        Gets compound or peptide retention time
        """
        ...
    
    def getRetentionTimeType(self) -> int:
        """
        Cython signature: RTType getRetentionTimeType()
        Get compound or peptide retentiontime type
        """
        ...
    
    def getRetentionTimeUnit(self) -> int:
        """
        Cython signature: RTUnit getRetentionTimeUnit()
        Get compound or peptide retentiontime unit (minute/seconds)
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
    
    def __richcmp__(self, other: Peptide, op: int) -> Any:
        ... 


class Prediction:
    """
    Cython implementation of _Prediction

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Prediction.html>`_
      -- Inherits from ['CVTermList']
    """
    
    software_ref: Union[bytes, str, String]
    
    contact_ref: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Prediction()
        """
        ...
    
    @overload
    def __init__(self, in_0: Prediction ) -> None:
        """
        Cython signature: void Prediction(Prediction &)
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
    
    def __richcmp__(self, other: Prediction, op: int) -> Any:
        ... 


class Protein:
    """
    Cython implementation of _Protein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Protein.html>`_
      -- Inherits from ['CVTermList']
    """
    
    id: Union[bytes, str, String]
    
    sequence: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Protein()
        """
        ...
    
    @overload
    def __init__(self, in_0: Protein ) -> None:
        """
        Cython signature: void Protein(Protein &)
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
    
    def __richcmp__(self, other: Protein, op: int) -> Any:
        ... 


class ProteinHit:
    """
    Cython implementation of _ProteinHit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinHit.html>`_
      -- Inherits from ['MetaInfoInterface']

    Representation of a protein hit
    
    It contains the fields score, score_type, rank, accession,
    sequence and coverage
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteinHit()
        """
        ...
    
    @overload
    def __init__(self, score: float , rank: int , accession: Union[bytes, str, String] , sequence: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void ProteinHit(double score, unsigned int rank, String accession, String sequence)
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteinHit ) -> None:
        """
        Cython signature: void ProteinHit(ProteinHit &)
        """
        ...
    
    def getScore(self) -> float:
        """
        Cython signature: float getScore()
        Returns the score of the protein hit
        """
        ...
    
    def getRank(self) -> int:
        """
        Cython signature: unsigned int getRank()
        Returns the rank of the protein hit
        """
        ...
    
    def getSequence(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSequence()
        Returns the protein sequence
        """
        ...
    
    def getAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getAccession()
        Returns the accession of the protein
        """
        ...
    
    def getDescription(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDescription()
        Returns the description of the protein
        """
        ...
    
    def getCoverage(self) -> float:
        """
        Cython signature: double getCoverage()
        Returns the coverage (in percent) of the protein hit based upon matched peptides
        """
        ...
    
    def setScore(self, in_0: float ) -> None:
        """
        Cython signature: void setScore(float)
        Sets the score of the protein hit
        """
        ...
    
    def setRank(self, in_0: int ) -> None:
        """
        Cython signature: void setRank(unsigned int)
        Sets the rank
        """
        ...
    
    def setSequence(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSequence(String)
        Sets the protein sequence
        """
        ...
    
    def setAccession(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setAccession(String)
        Sets the accession of the protein
        """
        ...
    
    def setDescription(self, description: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDescription(String description)
        Sets the description of the protein
        """
        ...
    
    def setCoverage(self, in_0: float ) -> None:
        """
        Cython signature: void setCoverage(double)
        Sets the coverage (in percent) of the protein hit based upon matched peptides
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
    
    def __richcmp__(self, other: ProteinHit, op: int) -> Any:
        ... 


class Publication:
    """
    Cython implementation of _Publication

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Publication.html>`_
      -- Inherits from ['CVTermList']
    """
    
    id: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Publication()
        """
        ...
    
    @overload
    def __init__(self, in_0: Publication ) -> None:
        """
        Cython signature: void Publication(Publication &)
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
    
    def __richcmp__(self, other: Publication, op: int) -> Any:
        ... 


class ReactionMonitoringTransition:
    """
    Cython implementation of _ReactionMonitoringTransition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ReactionMonitoringTransition.html>`_
      -- Inherits from ['CVTermList']

    This class stores a SRM/MRM transition
    
    This class is capable of representing a <Transition> tag in a TraML
    document completely and contains all associated information
    
    The default values for precursor m/z is 0.0 which indicates that it is
    uninitialized
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ReactionMonitoringTransition()
        """
        ...
    
    @overload
    def __init__(self, in_0: ReactionMonitoringTransition ) -> None:
        """
        Cython signature: void ReactionMonitoringTransition(ReactionMonitoringTransition &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def getNativeID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeID()
        """
        ...
    
    def getTransRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTransRef()
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        """
        ...
    
    def setNativeID(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeID(String name)
        """
        ...
    
    def setTransRef(self, trans_ref: Union[bytes, str, String] , type: int ) -> None:
        """
        Cython signature: void setTransRef(String trans_ref, TransType type)
        """
        ...
    
    def getTransType(self) -> int:
        """
        Cython signature: TransType getTransType()
        """
        ...
    
    def getProductMZ(self) -> float:
        """
        Cython signature: double getProductMZ()
        """
        ...
    
    def setProductMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setProductMZ(double)
        """
        ...
    
    def getPrecursorMZ(self) -> float:
        """
        Cython signature: double getPrecursorMZ()
        Returns the precursor mz (Q1 value)
        """
        ...
    
    def setPrecursorMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setPrecursorMZ(double)
        Sets the precursor mz (Q1 value)
        """
        ...
    
    def getDecoyTransitionType(self) -> int:
        """
        Cython signature: DecoyTransitionType getDecoyTransitionType()
        Returns the type of transition (target or decoy)
        """
        ...
    
    def hasPrecursorCVTerms(self) -> bool:
        """
        Cython signature: bool hasPrecursorCVTerms()
        Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
        """
        ...
    
    def setPrecursorCVTermList(self, list_: CVTermList ) -> None:
        """
        Cython signature: void setPrecursorCVTermList(CVTermList & list_)
        Sets a list of precursor CV Terms
        """
        ...
    
    def addPrecursorCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void addPrecursorCVTerm(CVTerm & cv_term)
        Adds precursor CV Term
        """
        ...
    
    def getPrecursorCVTermList(self) -> CVTermList:
        """
        Cython signature: CVTermList getPrecursorCVTermList()
        Obtains the list of CV Terms for the precursor
        """
        ...
    
    def addProductCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void addProductCVTerm(CVTerm & cv_term)
        """
        ...
    
    def getIntermediateProducts(self) -> List[TraMLProduct]:
        """
        Cython signature: libcpp_vector[TraMLProduct] getIntermediateProducts()
        """
        ...
    
    def addIntermediateProduct(self, product: TraMLProduct ) -> None:
        """
        Cython signature: void addIntermediateProduct(TraMLProduct product)
        """
        ...
    
    def setIntermediateProducts(self, products: List[TraMLProduct] ) -> None:
        """
        Cython signature: void setIntermediateProducts(libcpp_vector[TraMLProduct] & products)
        """
        ...
    
    def setProduct(self, product: TraMLProduct ) -> None:
        """
        Cython signature: void setProduct(TraMLProduct product)
        """
        ...
    
    def getProduct(self) -> TraMLProduct:
        """
        Cython signature: TraMLProduct getProduct()
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
    
    def setPrediction(self, prediction: Prediction ) -> None:
        """
        Cython signature: void setPrediction(Prediction & prediction)
        Sets prediction
        """
        ...
    
    def addPredictionTerm(self, prediction: CVTerm ) -> None:
        """
        Cython signature: void addPredictionTerm(CVTerm & prediction)
        Adds prediction term
        """
        ...
    
    def hasPrediction(self) -> bool:
        """
        Cython signature: bool hasPrediction()
        Returns true if a Prediction object exists (means it is safe to call getPrediction)
        """
        ...
    
    def getPrediction(self) -> Prediction:
        """
        Cython signature: Prediction getPrediction()
        Obtains the Prediction object
        """
        ...
    
    def setDecoyTransitionType(self, d: int ) -> None:
        """
        Cython signature: void setDecoyTransitionType(DecoyTransitionType & d)
        Sets the type of transition (target or decoy)
        """
        ...
    
    def getLibraryIntensity(self) -> float:
        """
        Cython signature: double getLibraryIntensity()
        Returns the library intensity (ion count or normalized ion count from a spectral library)
        """
        ...
    
    def setLibraryIntensity(self, intensity: float ) -> None:
        """
        Cython signature: void setLibraryIntensity(double intensity)
        Sets the library intensity (ion count or normalized ion count from a spectral library)
        """
        ...
    
    def getProductChargeState(self) -> int:
        """
        Cython signature: int getProductChargeState()
        Returns the charge state of the product
        """
        ...
    
    def isProductChargeStateSet(self) -> bool:
        """
        Cython signature: bool isProductChargeStateSet()
        Returns true if charge state of product is already set
        """
        ...
    
    def isDetectingTransition(self) -> bool:
        """
        Cython signature: bool isDetectingTransition()
        """
        ...
    
    def setDetectingTransition(self, val: bool ) -> None:
        """
        Cython signature: void setDetectingTransition(bool val)
        """
        ...
    
    def isIdentifyingTransition(self) -> bool:
        """
        Cython signature: bool isIdentifyingTransition()
        """
        ...
    
    def setIdentifyingTransition(self, val: bool ) -> None:
        """
        Cython signature: void setIdentifyingTransition(bool val)
        """
        ...
    
    def isQuantifyingTransition(self) -> bool:
        """
        Cython signature: bool isQuantifyingTransition()
        """
        ...
    
    def setQuantifyingTransition(self, val: bool ) -> None:
        """
        Cython signature: void setQuantifyingTransition(bool val)
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
    
    def __richcmp__(self, other: ReactionMonitoringTransition, op: int) -> Any:
        ... 


class RetentionTime:
    """
    Cython implementation of _RetentionTime

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1RetentionTime.html>`_
      -- Inherits from ['CVTermList']
    """
    
    software_ref: Union[bytes, str, String]
    
    retention_time_unit: int
    
    retention_time_type: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RetentionTime()
        """
        ...
    
    @overload
    def __init__(self, in_0: RetentionTime ) -> None:
        """
        Cython signature: void RetentionTime(RetentionTime &)
        """
        ...
    
    def isRTset(self) -> bool:
        """
        Cython signature: bool isRTset()
        """
        ...
    
    def setRT(self, rt: float ) -> None:
        """
        Cython signature: void setRT(double rt)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
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
    
    def __richcmp__(self, other: RetentionTime, op: int) -> Any:
        ...
    RTType : __RTType
    RTUnit : __RTUnit 


class SpectrumAlignmentScore:
    """
    Cython implementation of _SpectrumAlignmentScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAlignmentScore.html>`_
      -- Inherits from ['DefaultParamHandler']

    Similarity score via spectra alignment
    
    This class implements a simple scoring based on the alignment of spectra. This alignment
    is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
    of the peaks, minimizing the distances between the aligned peaks and maximizing the number
    of peak pairs
    
    The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the
    product of the intensities of the aligned peaks, with the given exponent (default is 2)
    sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAlignmentScore()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAlignmentScore ) -> None:
        """
        Cython signature: void SpectrumAlignmentScore(SpectrumAlignmentScore &)
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


class TargetedExperiment_Instrument:
    """
    Cython implementation of _TargetedExperiment_Instrument

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Instrument.html>`_
    """
    
    id: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TargetedExperiment_Instrument()
        """
        ...
    
    @overload
    def __init__(self, in_0: TargetedExperiment_Instrument ) -> None:
        """
        Cython signature: void TargetedExperiment_Instrument(TargetedExperiment_Instrument &)
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
    
    def consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void consumeCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] cv_term_map)
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
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        """
        ...
    
    def getKeysAsIntegers(self, keys: List[int] ) -> None:
        """
        Cython signature: void getKeysAsIntegers(libcpp_vector[unsigned int] & keys)
        """
        ...
    
    @overload
    def getMetaValue(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(unsigned int)
        """
        ...
    
    @overload
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        """
        ...
    
    @overload
    def setMetaValue(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(unsigned int, DataValue)
        """
        ...
    
    @overload
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        """
        ...
    
    @overload
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        """
        ...
    
    @overload
    def metaValueExists(self, in_0: int ) -> bool:
        """
        Cython signature: bool metaValueExists(unsigned int)
        """
        ...
    
    @overload
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        """
        ...
    
    @overload
    def removeMetaValue(self, in_0: int ) -> None:
        """
        Cython signature: void removeMetaValue(unsigned int)
        """
        ...
    
    def __richcmp__(self, other: TargetedExperiment_Instrument, op: int) -> Any:
        ... 


class TargetedExperiment_Interpretation:
    """
    Cython implementation of _TargetedExperiment_Interpretation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Interpretation.html>`_
    """
    
    ordinal: bytes
    
    rank: bytes
    
    iontype: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TargetedExperiment_Interpretation()
        """
        ...
    
    @overload
    def __init__(self, in_0: TargetedExperiment_Interpretation ) -> None:
        """
        Cython signature: void TargetedExperiment_Interpretation(TargetedExperiment_Interpretation &)
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
    
    def consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void consumeCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] cv_term_map)
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
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        """
        ...
    
    def getKeysAsIntegers(self, keys: List[int] ) -> None:
        """
        Cython signature: void getKeysAsIntegers(libcpp_vector[unsigned int] & keys)
        """
        ...
    
    @overload
    def getMetaValue(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(unsigned int)
        """
        ...
    
    @overload
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        """
        ...
    
    @overload
    def setMetaValue(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(unsigned int, DataValue)
        """
        ...
    
    @overload
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        """
        ...
    
    @overload
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        """
        ...
    
    @overload
    def metaValueExists(self, in_0: int ) -> bool:
        """
        Cython signature: bool metaValueExists(unsigned int)
        """
        ...
    
    @overload
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        """
        ...
    
    @overload
    def removeMetaValue(self, in_0: int ) -> None:
        """
        Cython signature: void removeMetaValue(unsigned int)
        """
        ...
    
    def __richcmp__(self, other: TargetedExperiment_Interpretation, op: int) -> Any:
        ... 


class TargetedExperiment_Modification:
    """
    Cython implementation of _TargetedExperiment_Modification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Modification.html>`_
    """
    
    avg_mass_delta: float
    
    mono_mass_delta: float
    
    location: int
    
    unimod_id: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TargetedExperiment_Modification()
        """
        ...
    
    @overload
    def __init__(self, in_0: TargetedExperiment_Modification ) -> None:
        """
        Cython signature: void TargetedExperiment_Modification(TargetedExperiment_Modification &)
        """
        ...
    
    def __richcmp__(self, other: TargetedExperiment_Modification, op: int) -> Any:
        ... 


class TheoreticalSpectrumGenerator:
    """
    Cython implementation of _TheoreticalSpectrumGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TheoreticalSpectrumGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TheoreticalSpectrumGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: TheoreticalSpectrumGenerator ) -> None:
        """
        Cython signature: void TheoreticalSpectrumGenerator(TheoreticalSpectrumGenerator &)
        """
        ...
    
    def getSpectrum(self, spec: MSSpectrum , peptide: AASequence , min_charge: int , max_charge: int ) -> None:
        """
        Cython signature: void getSpectrum(MSSpectrum & spec, AASequence & peptide, int min_charge, int max_charge)
        Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
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


class TraMLProduct:
    """
    Cython implementation of _TraMLProduct

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TraMLProduct.html>`_
      -- Inherits from ['CVTermList']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TraMLProduct()
        """
        ...
    
    @overload
    def __init__(self, in_0: TraMLProduct ) -> None:
        """
        Cython signature: void TraMLProduct(TraMLProduct &)
        """
        ...
    
    def setMZ(self, mz: float ) -> None:
        """
        Cython signature: void setMZ(double mz)
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        """
        ...
    
    def setChargeState(self, charge: int ) -> None:
        """
        Cython signature: void setChargeState(int charge)
        """
        ...
    
    def getChargeState(self) -> int:
        """
        Cython signature: int getChargeState()
        """
        ...
    
    def hasCharge(self) -> bool:
        """
        Cython signature: bool hasCharge()
        """
        ...
    
    def getConfigurationList(self) -> List[Configuration]:
        """
        Cython signature: libcpp_vector[Configuration] getConfigurationList()
        """
        ...
    
    def addConfiguration(self, configuration: Configuration ) -> None:
        """
        Cython signature: void addConfiguration(Configuration configuration)
        """
        ...
    
    def getInterpretationList(self) -> List[TargetedExperiment_Interpretation]:
        """
        Cython signature: libcpp_vector[TargetedExperiment_Interpretation] getInterpretationList()
        """
        ...
    
    def addInterpretation(self, interpretation: TargetedExperiment_Interpretation ) -> None:
        """
        Cython signature: void addInterpretation(TargetedExperiment_Interpretation interpretation)
        """
        ...
    
    def resetInterpretations(self) -> None:
        """
        Cython signature: void resetInterpretations()
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
    
    def __richcmp__(self, other: TraMLProduct, op: int) -> Any:
        ... 


class DecoyTransitionType:
    None
    UNKNOWN : int
    TARGET : int
    DECOY : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __RTType:
    None
    LOCAL : int
    NORMALIZED : int
    PREDICTED : int
    HPINS : int
    IRT : int
    UNKNOWN : int
    SIZE_OF_RTTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __RTUnit:
    None
    SECOND : int
    MINUTE : int
    UNKNOWN : int
    SIZE_OF_RTUNIT : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class TransType:
    None
    PEPTIDE : int
    NUCTIDE : int
    COMPOUND : int
    UNKNOWN : int

    def getMapping(self) -> Dict[int, str]:
       ... 

