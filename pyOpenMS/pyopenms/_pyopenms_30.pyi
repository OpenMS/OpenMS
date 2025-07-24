from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_AASequence_fromString(s: Union[bytes, str, String] ) -> AASequence:
    """
    Cython signature: AASequence fromString(String s)
        deprecated. Use AASequence(String) instead.
    """
    ...

def __static_AASequence_fromStringPermissive(s: Union[bytes, str, String] , permissive: bool ) -> AASequence:
    """
    Cython signature: AASequence fromStringPermissive(String s, bool permissive)
        deprecated. Use AASequence(String, bool) instead.
    """
    ...

def __static_DateTime_now() -> DateTime:
    """
    Cython signature: DateTime now()
    """
    ...


class AAIndex:
    """
    Cython implementation of _AAIndex

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AAIndex.html>`_
    """
    
    def aliphatic(self, aa: bytes ) -> float:
        """
        Cython signature: double aliphatic(char aa)
        """
        ...
    
    def acidic(self, aa: bytes ) -> float:
        """
        Cython signature: double acidic(char aa)
        """
        ...
    
    def basic(self, aa: bytes ) -> float:
        """
        Cython signature: double basic(char aa)
        """
        ...
    
    def polar(self, aa: bytes ) -> float:
        """
        Cython signature: double polar(char aa)
        """
        ...
    
    def getKHAG800101(self, aa: bytes ) -> float:
        """
        Cython signature: double getKHAG800101(char aa)
        """
        ...
    
    def getVASM830103(self, aa: bytes ) -> float:
        """
        Cython signature: double getVASM830103(char aa)
        """
        ...
    
    def getNADH010106(self, aa: bytes ) -> float:
        """
        Cython signature: double getNADH010106(char aa)
        """
        ...
    
    def getNADH010107(self, aa: bytes ) -> float:
        """
        Cython signature: double getNADH010107(char aa)
        """
        ...
    
    def getWILM950102(self, aa: bytes ) -> float:
        """
        Cython signature: double getWILM950102(char aa)
        """
        ...
    
    def getROBB760107(self, aa: bytes ) -> float:
        """
        Cython signature: double getROBB760107(char aa)
        """
        ...
    
    def getOOBM850104(self, aa: bytes ) -> float:
        """
        Cython signature: double getOOBM850104(char aa)
        """
        ...
    
    def getFAUJ880111(self, aa: bytes ) -> float:
        """
        Cython signature: double getFAUJ880111(char aa)
        """
        ...
    
    def getFINA770101(self, aa: bytes ) -> float:
        """
        Cython signature: double getFINA770101(char aa)
        """
        ...
    
    def getARGP820102(self, aa: bytes ) -> float:
        """
        Cython signature: double getARGP820102(char aa)
        """
        ...
    
    def calculateGB(self, seq: AASequence , T: float ) -> float:
        """
        Cython signature: double calculateGB(AASequence & seq, double T)
        """
        ... 


class AASequence:
    """
    Cython implementation of _AASequence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AASequence.html>`_

    Representation of a peptide/protein sequence
    This class represents amino acid sequences in OpenMS. An AASequence
    instance primarily contains a sequence of residues.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AASequence()
        """
        ...
    
    @overload
    def __init__(self, in_0: AASequence ) -> None:
        """
        Cython signature: void AASequence(AASequence &)
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void AASequence(const String &)
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE")
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[bytes, str, String] , permissive: bool ) -> None:
        """
        Cython signature: void AASequence(const String &, bool permissive)
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE"), permissive allows for '+', '*', and '#' in the sequence
        """
        ...
    
    def __add__(self: AASequence, other: AASequence) -> AASequence:
        ...
    
    def __iadd__(self: AASequence, other: AASequence) -> AASequence:
        ...
    
    def __getitem__(self, in_0: int ) -> Residue:
        """
        Cython signature: Residue operator[](size_t)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Check if sequence is empty
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the peptide as string with modifications embedded in brackets
        """
        ...
    
    def toUnmodifiedString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toUnmodifiedString()
        Returns the peptide as string without any modifications
        """
        ...
    
    def toUniModString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toUniModString()
        Returns the peptide as string with UniMod-style modifications embedded in brackets
        """
        ...
    
    @overload
    def toBracketString(self, ) -> Union[bytes, str, String]:
        """
        Cython signature: String toBracketString()
        Create a TPP compatible string of the modified sequence using bracket notation. Uses integer mass by default
        """
        ...
    
    @overload
    def toBracketString(self, integer_mass: bool ) -> Union[bytes, str, String]:
        """
        Cython signature: String toBracketString(bool integer_mass)
        Create a TPP compatible string of the modified sequence using bracket notation
        """
        ...
    
    @overload
    def toBracketString(self, integer_mass: bool , mass_delta: bool ) -> Union[bytes, str, String]:
        """
        Cython signature: String toBracketString(bool integer_mass, bool mass_delta)
        Create a TPP compatible string of the modified sequence using bracket notation.
        """
        ...
    
    @overload
    def toBracketString(self, integer_mass: bool , mass_delta: bool , fixed_modifications: List[bytes] ) -> Union[bytes, str, String]:
        """
        Cython signature: String toBracketString(bool integer_mass, bool mass_delta, libcpp_vector[String] fixed_modifications)
        Create a TPP compatible string of the modified sequence using bracket notation
        """
        ...
    
    @overload
    def setModification(self, index: int , modification: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setModification(size_t index, const String & modification)
        Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version
        """
        ...
    
    @overload
    def setModification(self, index: int , modification: ResidueModification ) -> None:
        """
        Cython signature: void setModification(size_t index, const ResidueModification & modification)
        Sets the modification of AA at index by providing a ResidueModification object. Stricter than just looking for the name and adds the Modification to the DB if not present
        """
        ...
    
    def setModificationByDiffMonoMass(self, index: int , diffMonoMass: float ) -> None:
        """
        Cython signature: void setModificationByDiffMonoMass(size_t index, double diffMonoMass)
        Modifies the residue at index in the sequence and potentially in the ResidueDB
        """
        ...
    
    @overload
    def setNTerminalModification(self, modification: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNTerminalModification(String modification)
        Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)
        """
        ...
    
    @overload
    def setNTerminalModification(self, mod: ResidueModification ) -> None:
        """
        Cython signature: void setNTerminalModification(const ResidueModification & mod)
        Sets the N-terminal modification (copies and adds to database if not present)
        """
        ...
    
    def setNTerminalModificationByDiffMonoMass(self, diffMonoMass: float , protein_term: bool ) -> None:
        """
        Cython signature: void setNTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term)
        Sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
        """
        ...
    
    @overload
    def setCTerminalModification(self, modification: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCTerminalModification(String modification)
        Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)
        """
        ...
    
    @overload
    def setCTerminalModification(self, mod: ResidueModification ) -> None:
        """
        Cython signature: void setCTerminalModification(const ResidueModification & mod)
        Sets the C-terminal modification (copies and adds to database if not present)
        """
        ...
    
    def setCTerminalModificationByDiffMonoMass(self, diffMonoMass: float , protein_term: bool ) -> None:
        """
        Cython signature: void setCTerminalModificationByDiffMonoMass(double diffMonoMass, bool protein_term)
        Sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
        """
        ...
    
    def getNTerminalModificationName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNTerminalModificationName()
        Returns the name (ID) of the N-terminal modification, or an empty string if none is set
        """
        ...
    
    def getNTerminalModification(self) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getNTerminalModification()
        Returns a copy of the name N-terminal modification object, or None
        """
        ...
    
    def getCTerminalModificationName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCTerminalModificationName()
        Returns the name (ID) of the C-terminal modification, or an empty string if none is set
        """
        ...
    
    def getCTerminalModification(self) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getCTerminalModification()
        Returns a copy of the name C-terminal modification object, or None
        """
        ...
    
    def getResidue(self, index: int ) -> Residue:
        """
        Cython signature: Residue getResidue(size_t index)
        Returns the residue at position index
        """
        ...
    
    @overload
    def getFormula(self, ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula()
        Convenience function with ResidueType=Full and charge = 0 by default
        """
        ...
    
    @overload
    def getFormula(self, type_: int , charge: int ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula(ResidueType type_, int charge)
        """
        ...
    
    @overload
    def getAverageWeight(self, ) -> float:
        """
        Cython signature: double getAverageWeight()
        Returns the average weight of the peptide
        """
        ...
    
    @overload
    def getAverageWeight(self, type_: int , charge: int ) -> float:
        """
        Cython signature: double getAverageWeight(ResidueType type_, int charge)
        """
        ...
    
    @overload
    def getMonoWeight(self, ) -> float:
        """
        Cython signature: double getMonoWeight()
        Returns the mono isotopic weight of the peptide
        """
        ...
    
    @overload
    def getMonoWeight(self, type_: int , charge: int ) -> float:
        """
        Cython signature: double getMonoWeight(ResidueType type_, int charge)
        """
        ...
    
    @overload
    def getMZ(self, charge: int ) -> float:
        """
        Cython signature: double getMZ(int charge)
        Returns the mass-to-charge ratio of the peptide
        """
        ...
    
    @overload
    def getMZ(self, charge: int , type_: int ) -> float:
        """
        Cython signature: double getMZ(int charge, ResidueType type_)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Returns the number of residues
        """
        ...
    
    def getPrefix(self, index: int ) -> AASequence:
        """
        Cython signature: AASequence getPrefix(size_t index)
        Returns a peptide sequence of the first index residues
        """
        ...
    
    def getSuffix(self, index: int ) -> AASequence:
        """
        Cython signature: AASequence getSuffix(size_t index)
        Returns a peptide sequence of the last index residues
        """
        ...
    
    def getSubsequence(self, index: int , number: int ) -> AASequence:
        """
        Cython signature: AASequence getSubsequence(size_t index, unsigned int number)
        Returns a peptide sequence of number residues, beginning at position index
        """
        ...
    
    def has(self, residue: Residue ) -> bool:
        """
        Cython signature: bool has(Residue residue)
        Returns true if the peptide contains the given residue
        """
        ...
    
    def hasSubsequence(self, peptide: AASequence ) -> bool:
        """
        Cython signature: bool hasSubsequence(AASequence peptide)
        Returns true if the peptide contains the given peptide
        """
        ...
    
    def hasPrefix(self, peptide: AASequence ) -> bool:
        """
        Cython signature: bool hasPrefix(AASequence peptide)
        Returns true if the peptide has the given prefix
        """
        ...
    
    def hasSuffix(self, peptide: AASequence ) -> bool:
        """
        Cython signature: bool hasSuffix(AASequence peptide)
        Returns true if the peptide has the given suffix
        """
        ...
    
    def hasNTerminalModification(self) -> bool:
        """
        Cython signature: bool hasNTerminalModification()
        Predicate which is true if the peptide is N-term modified
        """
        ...
    
    def hasCTerminalModification(self) -> bool:
        """
        Cython signature: bool hasCTerminalModification()
        Predicate which is true if the peptide is C-term modified
        """
        ...
    
    def isModified(self) -> bool:
        """
        Cython signature: bool isModified()
        Returns true if any of the residues or termini are modified
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the peptide as string with modifications embedded in brackets
        """
        ...
    
    def __richcmp__(self, other: AASequence, op: int) -> Any:
        ...
    
    fromString: __static_AASequence_fromString
    
    fromStringPermissive: __static_AASequence_fromStringPermissive 


class Compomer:
    """
    Cython implementation of _Compomer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Compomer.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Compomer()
        """
        ...
    
    @overload
    def __init__(self, in_0: Compomer ) -> None:
        """
        Cython signature: void Compomer(Compomer &)
        """
        ...
    
    def add(self, a: Adduct , side: int ) -> None:
        """
        Cython signature: void add(Adduct & a, unsigned int side)
        """
        ...
    
    def isConflicting(self, cmp: Compomer , side_this: int , side_other: int ) -> bool:
        """
        Cython signature: bool isConflicting(Compomer & cmp, unsigned int side_this, unsigned int side_other)
        """
        ...
    
    def setID(self, id: int ) -> None:
        """
        Cython signature: void setID(size_t id)
        Sets an Id which allows unique identification of a compomer
        """
        ...
    
    def getID(self) -> int:
        """
        Cython signature: size_t getID()
        Returns Id which allows unique identification of this compomer
        """
        ...
    
    def getNetCharge(self) -> int:
        """
        Cython signature: int getNetCharge()
        Net charge of compomer (i.e. difference between left and right side of compomer)
        """
        ...
    
    def getMass(self) -> float:
        """
        Cython signature: double getMass()
        Mass of all contained adducts
        """
        ...
    
    def getPositiveCharges(self) -> int:
        """
        Cython signature: int getPositiveCharges()
        Summed positive charges of contained adducts
        """
        ...
    
    def getNegativeCharges(self) -> int:
        """
        Cython signature: int getNegativeCharges()
        Summed negative charges of contained adducts
        """
        ...
    
    def getLogP(self) -> float:
        """
        Cython signature: double getLogP()
        Returns the log probability
        """
        ...
    
    def getRTShift(self) -> float:
        """
        Cython signature: double getRTShift()
        Returns the log probability
        """
        ...
    
    @overload
    def getAdductsAsString(self, ) -> Union[bytes, str, String]:
        """
        Cython signature: String getAdductsAsString()
        Get adducts with their abundance as compact string for both sides
        """
        ...
    
    @overload
    def getAdductsAsString(self, side: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getAdductsAsString(unsigned int side)
        Get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
        """
        ...
    
    def isSingleAdduct(self, a: Adduct , side: int ) -> bool:
        """
        Cython signature: bool isSingleAdduct(Adduct & a, unsigned int side)
        Check if Compomer only contains a single adduct on side @p side
        """
        ...
    
    @overload
    def removeAdduct(self, a: Adduct ) -> Compomer:
        """
        Cython signature: Compomer removeAdduct(Adduct & a)
        Remove ALL instances of the given adduct
        """
        ...
    
    @overload
    def removeAdduct(self, a: Adduct , side: int ) -> Compomer:
        """
        Cython signature: Compomer removeAdduct(Adduct & a, unsigned int side)
        """
        ...
    
    def getLabels(self, side: int ) -> List[bytes]:
        """
        Cython signature: StringList getLabels(unsigned int side)
        Returns the adduct labels from parameter(side) given. (LEFT or RIGHT)
        """
        ... 


class ConsensusIDAlgorithmBest:
    """
    Cython implementation of _ConsensusIDAlgorithmBest

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmBest.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmBest()
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


class DateTime:
    """
    Cython implementation of _DateTime

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DateTime.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DateTime()
        """
        ...
    
    @overload
    def __init__(self, in_0: DateTime ) -> None:
        """
        Cython signature: void DateTime(DateTime &)
        """
        ...
    
    def setDate(self, date: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDate(String date)
        """
        ...
    
    def setTime(self, date: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTime(String date)
        """
        ...
    
    def getDate(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDate()
        """
        ...
    
    def getTime(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTime()
        """
        ...
    
    def now(self) -> DateTime:
        """
        Cython signature: DateTime now()
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def get(self) -> Union[bytes, str, String]:
        """
        Cython signature: String get()
        """
        ...
    
    def set(self, date: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void set(String date)
        """
        ...
    
    now: __static_DateTime_now 


class ElementDB:
    """
    Cython implementation of _ElementDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElementDB.html>`_
    """
    
    @overload
    def getElement(self, name: Union[bytes, str, String] ) -> Element:
        """
        Cython signature: const Element * getElement(const String & name)
        """
        ...
    
    @overload
    def getElement(self, atomic_number: int ) -> Element:
        """
        Cython signature: const Element * getElement(unsigned int atomic_number)
        """
        ...
    
    def addElement(self, name: bytes , symbol: bytes , an: int , abundance: Dict[int, float] , mass: Dict[int, float] , replace_existing: bool ) -> None:
        """
        Cython signature: void addElement(libcpp_string name, libcpp_string symbol, unsigned int an, libcpp_map[unsigned int,double] abundance, libcpp_map[unsigned int,double] mass, bool replace_existing)
        """
        ...
    
    @overload
    def hasElement(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasElement(const String & name)
        Returns true if the db contains an element with the given name, else false
        """
        ...
    
    @overload
    def hasElement(self, atomic_number: int ) -> bool:
        """
        Cython signature: bool hasElement(unsigned int atomic_number)
        Returns true if the db contains an element with the given atomic_number, else false
        """
        ... 


class IsotopeFitter1D:
    """
    Cython implementation of _IsotopeFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeFitter1D.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsotopeFitter1D()
        Isotope distribution fitter (1-dim.) approximated using linear interpolation
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopeFitter1D ) -> None:
        """
        Cython signature: void IsotopeFitter1D(IsotopeFitter1D &)
        """
        ... 


class MetaInfoRegistry:
    """
    Cython implementation of _MetaInfoRegistry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoRegistry.html>`_

    Registry which assigns unique integer indices to strings
    
    When registering a new name an index >= 1024 is assigned.
    Indices from 1 to 1023 are reserved for fast access and will never change:
    1 - isotopic_range
    2 - cluster_id
    3 - label
    4 - icon
    5 - color
    6 - RT
    7 - MZ
    8 - predicted_RT
    9 - predicted_RT_p_value
    10 - spectrum_reference
    11 - ID
    12 - low_quality
    13 - charge
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaInfoRegistry()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaInfoRegistry ) -> None:
        """
        Cython signature: void MetaInfoRegistry(MetaInfoRegistry &)
        """
        ...
    
    def registerName(self, name: Union[bytes, str, String] , description: Union[bytes, str, String] , unit: Union[bytes, str, String] ) -> int:
        """
        Cython signature: unsigned int registerName(const String & name, const String & description, const String & unit)
        Registers a string, stores its description and unit, and returns the corresponding index. If the string is already registered, it returns the index of the string
        """
        ...
    
    @overload
    def setDescription(self, index: int , description: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDescription(unsigned int index, const String & description)
        Sets the description (String), corresponding to an index
        """
        ...
    
    @overload
    def setDescription(self, name: Union[bytes, str, String] , description: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDescription(const String & name, const String & description)
        Sets the description (String), corresponding to a name
        """
        ...
    
    @overload
    def setUnit(self, index: int , unit: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setUnit(unsigned int index, const String & unit)
        Sets the unit (String), corresponding to an index
        """
        ...
    
    @overload
    def setUnit(self, name: Union[bytes, str, String] , unit: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setUnit(const String & name, const String & unit)
        Sets the unit (String), corresponding to a name
        """
        ...
    
    def getIndex(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: unsigned int getIndex(const String & name)
        Returns the integer index corresponding to a string. If the string is not registered, returns UInt(-1) (= UINT_MAX)
        """
        ...
    
    def getName(self, index: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getName(unsigned int index)
        Returns the corresponding name to an index
        """
        ...
    
    @overload
    def getDescription(self, index: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getDescription(unsigned int index)
        Returns the description of an index
        """
        ...
    
    @overload
    def getDescription(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String getDescription(const String & name)
        Returns the description of a name
        """
        ...
    
    @overload
    def getUnit(self, index: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getUnit(unsigned int index)
        Returns the unit of an index
        """
        ...
    
    @overload
    def getUnit(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String getUnit(const String & name)
        Returns the unit of a name
        """
        ... 


class OPXLHelper:
    """
    Cython implementation of _OPXLHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLHelper.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OPXLHelper()
        """
        ...
    
    @overload
    def __init__(self, in_0: OPXLHelper ) -> None:
        """
        Cython signature: void OPXLHelper(OPXLHelper &)
        """
        ...
    
    def enumerateCrossLinksAndMasses(self, peptides: List[AASeqWithMass] , cross_link_mass_light: float , cross_link_mass_mono_link: List[float] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , spectrum_precursors: List[float] , precursor_correction_positions: List[int] , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool ) -> List[XLPrecursor]:
        """
        Cython signature: libcpp_vector[XLPrecursor] enumerateCrossLinksAndMasses(libcpp_vector[AASeqWithMass] peptides, double cross_link_mass_light, DoubleList cross_link_mass_mono_link, StringList cross_link_residue1, StringList cross_link_residue2, libcpp_vector[double] & spectrum_precursors, libcpp_vector[int] & precursor_correction_positions, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
        """
        ...
    
    def digestDatabase(self, fasta_db: List[FASTAEntry] , digestor: EnzymaticDigestion , min_peptide_length: int , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , fixed_modifications: ModifiedPeptideGenerator_MapToResidueType , variable_modifications: ModifiedPeptideGenerator_MapToResidueType , max_variable_mods_per_peptide: int ) -> List[AASeqWithMass]:
        """
        Cython signature: libcpp_vector[AASeqWithMass] digestDatabase(libcpp_vector[FASTAEntry] fasta_db, EnzymaticDigestion digestor, size_t min_peptide_length, StringList cross_link_residue1, StringList cross_link_residue2, ModifiedPeptideGenerator_MapToResidueType & fixed_modifications, ModifiedPeptideGenerator_MapToResidueType & variable_modifications, size_t max_variable_mods_per_peptide)
        """
        ...
    
    def buildCandidates(self, candidates: List[XLPrecursor] , precursor_corrections: List[int] , precursor_correction_positions: List[int] , peptide_masses: List[AASeqWithMass] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , cross_link_mass: float , cross_link_mass_mono_link: List[float] , spectrum_precursor_vector: List[float] , allowed_error_vector: List[float] , cross_link_name: Union[bytes, str, String] ) -> List[ProteinProteinCrossLink]:
        """
        Cython signature: libcpp_vector[ProteinProteinCrossLink] buildCandidates(libcpp_vector[XLPrecursor] & candidates, libcpp_vector[int] & precursor_corrections, libcpp_vector[int] & precursor_correction_positions, libcpp_vector[AASeqWithMass] & peptide_masses, const StringList & cross_link_residue1, const StringList & cross_link_residue2, double cross_link_mass, DoubleList cross_link_mass_mono_link, libcpp_vector[double] & spectrum_precursor_vector, libcpp_vector[double] & allowed_error_vector, String cross_link_name)
        """
        ...
    
    def buildFragmentAnnotations(self, frag_annotations: List[PeptideHit_PeakAnnotation] , matching: List[List[int, int]] , theoretical_spectrum: MSSpectrum , experiment_spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void buildFragmentAnnotations(libcpp_vector[PeptideHit_PeakAnnotation] & frag_annotations, libcpp_vector[libcpp_pair[size_t,size_t]] matching, MSSpectrum theoretical_spectrum, MSSpectrum experiment_spectrum)
        """
        ...
    
    def buildPeptideIDs(self, peptide_ids: PeptideIdentificationList , top_csms_spectrum: List[CrossLinkSpectrumMatch] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , all_top_csms_current_index: int , spectra: MSExperiment , scan_index: int , scan_index_heavy: int ) -> None:
        """
        Cython signature: void buildPeptideIDs(PeptideIdentificationList & peptide_ids, libcpp_vector[CrossLinkSpectrumMatch] top_csms_spectrum, libcpp_vector[libcpp_vector[CrossLinkSpectrumMatch]] & all_top_csms, size_t all_top_csms_current_index, MSExperiment spectra, size_t scan_index, size_t scan_index_heavy)
        """
        ...
    
    def addProteinPositionMetaValues(self, peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void addProteinPositionMetaValues(PeptideIdentificationList & peptide_ids)
        """
        ...
    
    def addXLTargetDecoyMV(self, peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void addXLTargetDecoyMV(PeptideIdentificationList & peptide_ids)
        """
        ...
    
    def addBetaAccessions(self, peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void addBetaAccessions(PeptideIdentificationList & peptide_ids)
        """
        ...
    
    def removeBetaPeptideHits(self, peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void removeBetaPeptideHits(PeptideIdentificationList & peptide_ids)
        """
        ...
    
    def addPercolatorFeatureList(self, prot_id: ProteinIdentification ) -> None:
        """
        Cython signature: void addPercolatorFeatureList(ProteinIdentification & prot_id)
        """
        ...
    
    def computeDeltaScores(self, peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void computeDeltaScores(PeptideIdentificationList & peptide_ids)
        """
        ...
    
    def combineTopRanksFromPairs(self, peptide_ids: PeptideIdentificationList , number_top_hits: int ) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList combineTopRanksFromPairs(PeptideIdentificationList & peptide_ids, size_t number_top_hits)
        """
        ...
    
    def collectPrecursorCandidates(self, precursor_correction_steps: List[int] , precursor_mass: float , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool , filtered_peptide_masses: List[AASeqWithMass] , cross_link_mass: float , cross_link_mass_mono_link: List[float] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , cross_link_name: Union[bytes, str, String] , use_sequence_tags: bool , tags: List[Union[bytes, str]] ) -> List[ProteinProteinCrossLink]:
        """
        Cython signature: libcpp_vector[ProteinProteinCrossLink] collectPrecursorCandidates(IntList precursor_correction_steps, double precursor_mass, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, libcpp_vector[AASeqWithMass] filtered_peptide_masses, double cross_link_mass, DoubleList cross_link_mass_mono_link, StringList cross_link_residue1, StringList cross_link_residue2, String cross_link_name, bool use_sequence_tags, const libcpp_vector[libcpp_utf8_string] & tags)
        """
        ...
    
    def computePrecursorError(self, csm: CrossLinkSpectrumMatch , precursor_mz: float , precursor_charge: int ) -> float:
        """
        Cython signature: double computePrecursorError(CrossLinkSpectrumMatch csm, double precursor_mz, int precursor_charge)
        """
        ...
    
    def isoPeakMeans(self, csm: CrossLinkSpectrumMatch , num_iso_peaks_array: IntegerDataArray , matched_spec_linear_alpha: List[List[int, int]] , matched_spec_linear_beta: List[List[int, int]] , matched_spec_xlinks_alpha: List[List[int, int]] , matched_spec_xlinks_beta: List[List[int, int]] ) -> None:
        """
        Cython signature: void isoPeakMeans(CrossLinkSpectrumMatch & csm, IntegerDataArray & num_iso_peaks_array, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_linear_alpha, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_linear_beta, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_xlinks_alpha, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_xlinks_beta)
        """
        ... 


class PepXMLFileMascot:
    """
    Cython implementation of _PepXMLFileMascot

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PepXMLFileMascot.html>`_

    Used to load Mascot PepXML files
    
    A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void PepXMLFileMascot()
        """
        ... 


class PosteriorErrorProbabilityModel:
    """
    Cython implementation of _PosteriorErrorProbabilityModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1PosteriorErrorProbabilityModel.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void PosteriorErrorProbabilityModel()
        """
        ...
    
    @overload
    def fit(self, search_engine_scores: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool fit(libcpp_vector[double] & search_engine_scores, String outlier_handling)
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
        """
        ...
    
    @overload
    def fit(self, search_engine_scores: List[float] , probabilities: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool fit(libcpp_vector[double] & search_engine_scores, libcpp_vector[double] & probabilities, String outlier_handling)
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :param probabilities: A vector which holds the probability for each data point after running this function. If it has some content it will be overwritten
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
        """
        ...
    
    def fillDensities(self, x_scores: List[float] , incorrect_density: List[float] , correct_density: List[float] ) -> None:
        """
        Cython signature: void fillDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density)
        Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences
        """
        ...
    
    def fillLogDensities(self, x_scores: List[float] , incorrect_density: List[float] , correct_density: List[float] ) -> None:
        """
        Cython signature: void fillLogDensities(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density)
        Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences
        """
        ...
    
    def computeLogLikelihood(self, incorrect_density: List[float] , correct_density: List[float] ) -> float:
        """
        Cython signature: double computeLogLikelihood(libcpp_vector[double] & incorrect_density, libcpp_vector[double] & correct_density)
        Computes the Maximum Likelihood with a log-likelihood function
        """
        ...
    
    def pos_neg_mean_weighted_posteriors(self, x_scores: List[float] , incorrect_posteriors: List[float] ) -> List[float, float]:
        """
        Cython signature: libcpp_pair[double,double] pos_neg_mean_weighted_posteriors(libcpp_vector[double] & x_scores, libcpp_vector[double] & incorrect_posteriors)
        """
        ...
    
    def getCorrectlyAssignedFitResult(self) -> GaussFitResult:
        """
        Cython signature: GaussFitResult getCorrectlyAssignedFitResult()
        Returns estimated parameters for correctly assigned sequences. Fit should be used before
        """
        ...
    
    def getIncorrectlyAssignedFitResult(self) -> GaussFitResult:
        """
        Cython signature: GaussFitResult getIncorrectlyAssignedFitResult()
        Returns estimated parameters for correctly assigned sequences. Fit should be used before
        """
        ...
    
    def getNegativePrior(self) -> float:
        """
        Cython signature: double getNegativePrior()
        Returns the estimated negative prior probability
        """
        ...
    
    def computeProbability(self, score: float ) -> float:
        """
        Cython signature: double computeProbability(double score)
        Returns the computed posterior error probability for a given score
        """
        ...
    
    def initPlots(self, x_scores: List[float] ) -> TextFile:
        """
        Cython signature: TextFile initPlots(libcpp_vector[double] & x_scores)
        Initializes the plots
        """
        ...
    
    def getGumbelGnuplotFormula(self, params: GaussFitResult ) -> Union[bytes, str, String]:
        """
        Cython signature: String getGumbelGnuplotFormula(GaussFitResult & params)
        Returns the gnuplot formula of the fitted gumbel distribution
        """
        ...
    
    def getGaussGnuplotFormula(self, params: GaussFitResult ) -> Union[bytes, str, String]:
        """
        Cython signature: String getGaussGnuplotFormula(GaussFitResult & params)
        Returns the gnuplot formula of the fitted gauss distribution
        """
        ...
    
    def getBothGnuplotFormula(self, incorrect: GaussFitResult , correct: GaussFitResult ) -> Union[bytes, str, String]:
        """
        Cython signature: String getBothGnuplotFormula(GaussFitResult & incorrect, GaussFitResult & correct)
        Returns the gnuplot formula of the fitted mixture distribution
        """
        ...
    
    def plotTargetDecoyEstimation(self, target: List[float] , decoy: List[float] ) -> None:
        """
        Cython signature: void plotTargetDecoyEstimation(libcpp_vector[double] & target, libcpp_vector[double] & decoy)
        Plots the estimated distribution against target and decoy hits
        """
        ...
    
    def getSmallestScore(self) -> float:
        """
        Cython signature: double getSmallestScore()
        Returns the smallest score used in the last fit
        """
        ...
    
    def tryGnuplot(self, gp_file: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void tryGnuplot(const String & gp_file)
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


class RNaseDigestion:
    """
    Cython implementation of _RNaseDigestion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RNaseDigestion.html>`_
      -- Inherits from ['EnzymaticDigestion']

    Class for the enzymatic digestion of RNA
    
    Usage:
    
    .. code-block:: python
    
          from pyopenms import *
          oligo = NASequence.fromString("pAUGUCGCAG");
    
          dig = RNaseDigestion()
          dig.setEnzyme("RNase_T1")
    
          result = []
          dig.digest(oligo, result)
          for fragment in result:
            print (fragment)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RNaseDigestion()
        """
        ...
    
    @overload
    def __init__(self, in_0: RNaseDigestion ) -> None:
        """
        Cython signature: void RNaseDigestion(RNaseDigestion &)
        """
        ...
    
    @overload
    def setEnzyme(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setEnzyme(String name)
        Sets the enzyme for the digestion (by name)
        """
        ...
    
    @overload
    def setEnzyme(self, enzyme: DigestionEnzyme ) -> None:
        """
        Cython signature: void setEnzyme(DigestionEnzyme * enzyme)
        Sets the enzyme for the digestion
        """
        ...
    
    @overload
    def digest(self, rna: NASequence , output: List[NASequence] ) -> None:
        """
        Cython signature: void digest(NASequence & rna, libcpp_vector[NASequence] & output)
        """
        ...
    
    @overload
    def digest(self, rna: NASequence , output: List[NASequence] , min_length: int , max_length: int ) -> None:
        """
        Cython signature: void digest(NASequence & rna, libcpp_vector[NASequence] & output, size_t min_length, size_t max_length)
        Performs the enzymatic digestion of a (potentially modified) RNA
        
        :param rna: Sequence to digest
        :param output: Digestion productsq
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
        Performs the enzymatic digestion of all RNA parent molecules in IdentificationData (id_data)
        
        :param id_data: IdentificationData object which includes sequences to digest
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
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


class TextFile:
    """
    Cython implementation of _TextFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TextFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TextFile()
        This class provides some basic file handling methods for text files
        """
        ...
    
    @overload
    def __init__(self, in_0: TextFile ) -> None:
        """
        Cython signature: void TextFile(TextFile &)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , trim_linesalse: bool , first_n1: int ) -> None:
        """
        Cython signature: void TextFile(const String & filename, bool trim_linesalse, int first_n1)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , trim_linesalse: bool , first_n1: int ) -> None:
        """
        Cython signature: void load(const String & filename, bool trim_linesalse, int first_n1)
        Loads data from a text file
        
        :param filename: The input file name
        :param trim_lines: Whether or not the lines are trimmed when reading them from file
        :param first_n: If set, only `first_n` lines the lines from the beginning of the file are read
        :param skip_empty_lines: Should empty lines be skipped? If used in conjunction with `trim_lines`, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & filename)
        Writes the data to a file
        """
        ...
    
    def addLine(self, line: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addLine(const String line)
        """
        ... 


class XMLFile:
    """
    Cython implementation of _XMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1XMLFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void XMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: XMLFile ) -> None:
        """
        Cython signature: void XMLFile(XMLFile &)
        """
        ...
    
    @overload
    def __init__(self, schema_location: Union[bytes, str, String] , version: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void XMLFile(const String & schema_location, const String & version)
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
        """
        ... 


class SIDE:
    None
    LEFT : int
    RIGHT : int
    BOTH : int

    def getMapping(self) -> Dict[int, str]:
       ... 

