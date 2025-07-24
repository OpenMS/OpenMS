from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_NASequence_fromString(s: Union[bytes, str, String] ) -> NASequence:
    """
    Cython signature: NASequence fromString(const String & s)
    """
    ...

def __static_PeptideProteinResolution_run(proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList ) -> None:
    """
    Cython signature: void run(libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList & peptides)
    """
    ...


class CVMappings:
    """
    Cython implementation of _CVMappings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappings.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVMappings()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVMappings ) -> None:
        """
        Cython signature: void CVMappings(CVMappings &)
        """
        ...
    
    def setMappingRules(self, cv_mapping_rules: List[CVMappingRule] ) -> None:
        """
        Cython signature: void setMappingRules(libcpp_vector[CVMappingRule] & cv_mapping_rules)
        Sets the mapping rules of the mapping file
        """
        ...
    
    def getMappingRules(self) -> List[CVMappingRule]:
        """
        Cython signature: libcpp_vector[CVMappingRule] getMappingRules()
        Returns the mapping rules
        """
        ...
    
    def addMappingRule(self, cv_mapping_rule: CVMappingRule ) -> None:
        """
        Cython signature: void addMappingRule(CVMappingRule & cv_mapping_rule)
        Adds a mapping rule
        """
        ...
    
    def setCVReferences(self, cv_references: List[CVReference] ) -> None:
        """
        Cython signature: void setCVReferences(libcpp_vector[CVReference] & cv_references)
        Sets the CV references
        """
        ...
    
    def getCVReferences(self) -> List[CVReference]:
        """
        Cython signature: libcpp_vector[CVReference] getCVReferences()
        Returns the CV references
        """
        ...
    
    def addCVReference(self, cv_reference: CVReference ) -> None:
        """
        Cython signature: void addCVReference(CVReference & cv_reference)
        Adds a CV reference
        """
        ...
    
    def hasCVReference(self, identifier: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasCVReference(const String & identifier)
        Returns true if a CV reference is given
        """
        ... 


class ChargePair:
    """
    Cython implementation of _ChargePair

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChargePair.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChargePair()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChargePair ) -> None:
        """
        Cython signature: void ChargePair(ChargePair &)
        """
        ...
    
    @overload
    def __init__(self, index0: int , index1: int , charge0: int , charge1: int , compomer: Compomer , mass_diff: float , active: bool ) -> None:
        """
        Cython signature: void ChargePair(size_t index0, size_t index1, int charge0, int charge1, Compomer compomer, double mass_diff, bool active)
        """
        ...
    
    def getCharge(self, pairID: int ) -> int:
        """
        Cython signature: int getCharge(unsigned int pairID)
        Returns the charge (for element 0 or 1)
        """
        ...
    
    def setCharge(self, pairID: int , e: int ) -> None:
        """
        Cython signature: void setCharge(unsigned int pairID, int e)
        Sets the charge (for element 0 or 1)
        """
        ...
    
    def getElementIndex(self, pairID: int ) -> int:
        """
        Cython signature: size_t getElementIndex(unsigned int pairID)
        Returns the element index (for element 0 or 1)
        """
        ...
    
    def setElementIndex(self, pairID: int , e: int ) -> None:
        """
        Cython signature: void setElementIndex(unsigned int pairID, size_t e)
        Sets the element index (for element 0 or 1)
        """
        ...
    
    def getCompomer(self) -> Compomer:
        """
        Cython signature: Compomer getCompomer()
        Returns the Id of the compomer that explains the mass difference
        """
        ...
    
    def setCompomer(self, compomer: Compomer ) -> None:
        """
        Cython signature: void setCompomer(Compomer & compomer)
        Sets the compomer id
        """
        ...
    
    def getMassDiff(self) -> float:
        """
        Cython signature: double getMassDiff()
        Returns the mass difference
        """
        ...
    
    def setMassDiff(self, mass_diff: float ) -> None:
        """
        Cython signature: void setMassDiff(double mass_diff)
        Sets the mass difference
        """
        ...
    
    def getEdgeScore(self) -> float:
        """
        Cython signature: double getEdgeScore()
        Returns the ILP edge score
        """
        ...
    
    def setEdgeScore(self, score: float ) -> None:
        """
        Cython signature: void setEdgeScore(double score)
        Sets the ILP edge score
        """
        ...
    
    def isActive(self) -> bool:
        """
        Cython signature: bool isActive()
        Is this pair realized?
        """
        ...
    
    def setActive(self, active: bool ) -> None:
        """
        Cython signature: void setActive(bool active)
        """
        ... 


class Date:
    """
    Cython implementation of _Date

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Date.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Date()
        """
        ...
    
    @overload
    def __init__(self, in_0: Date ) -> None:
        """
        Cython signature: void Date(Date &)
        """
        ...
    
    def set(self, date: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void set(const String & date)
        """
        ...
    
    def today(self) -> Date:
        """
        Cython signature: Date today()
        """
        ...
    
    def get(self) -> Union[bytes, str, String]:
        """
        Cython signature: String get()
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ... 


class Feature:
    """
    Cython implementation of _Feature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Feature.html>`_
      -- Inherits from ['UniqueIdInterface', 'RichPeak2D']

    An LC-MS feature
    
    The Feature class is used to describe the two-dimensional signal caused by an
    analyte. It can store a charge state and a list of peptide identifications
    (for peptides). The area occupied by the Feature in the LC-MS data set is
    represented by a list of convex hulls (one for each isotopic peak). There is
    also a convex hull for the entire Feature. The model description can store
    the parameters of a two-dimensional theoretical model of the underlying
    signal in LC-MS. Currently, non-peptide compounds are also represented as
    features
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Feature()
        """
        ...
    
    @overload
    def __init__(self, in_0: Feature ) -> None:
        """
        Cython signature: void Feature(Feature &)
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
    
    def __richcmp__(self, other: Feature, op: int) -> Any:
        ... 


class MRMFeatureQCFile:
    """
    Cython implementation of _MRMFeatureQCFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureQCFile.html>`_

    File adapter for MRMFeatureQC files
    
    Loads and stores .csv or .tsv files describing an MRMFeatureQC
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeatureQCFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeatureQCFile ) -> None:
        """
        Cython signature: void MRMFeatureQCFile(MRMFeatureQCFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , mrmfqc: MRMFeatureQC , is_component_group: bool ) -> None:
        """
        Cython signature: void load(const String & filename, MRMFeatureQC & mrmfqc, const bool is_component_group)
        Loads an MRMFeatureQC file
        
        
        :param filename: The path to the input file
        :param mrmfqc: The output class which will contain the criteria
        :param is_component_group: True if the user intends to load ComponentGroupQCs data, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ... 


class MultiplexIsotopicPeakPattern:
    """
    Cython implementation of _MultiplexIsotopicPeakPattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexIsotopicPeakPattern.html>`_
    """
    
    @overload
    def __init__(self, c: int , ppp: int , ms: MultiplexDeltaMasses , msi: int ) -> None:
        """
        Cython signature: void MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi)
        """
        ...
    
    @overload
    def __init__(self, in_0: MultiplexIsotopicPeakPattern ) -> None:
        """
        Cython signature: void MultiplexIsotopicPeakPattern(MultiplexIsotopicPeakPattern &)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        Returns charge
        """
        ...
    
    def getPeaksPerPeptide(self) -> int:
        """
        Cython signature: int getPeaksPerPeptide()
        Returns peaks per peptide
        """
        ...
    
    def getMassShifts(self) -> MultiplexDeltaMasses:
        """
        Cython signature: MultiplexDeltaMasses getMassShifts()
        Returns mass shifts
        """
        ...
    
    def getMassShiftIndex(self) -> int:
        """
        Cython signature: int getMassShiftIndex()
        Returns mass shift index
        """
        ...
    
    def getMassShiftCount(self) -> int:
        """
        Cython signature: unsigned int getMassShiftCount()
        Returns number of mass shifts i.e. the number of peptides in the multiplet
        """
        ...
    
    def getMassShiftAt(self, i: int ) -> float:
        """
        Cython signature: double getMassShiftAt(int i)
        Returns mass shift at position i
        """
        ...
    
    def getMZShiftAt(self, i: int ) -> float:
        """
        Cython signature: double getMZShiftAt(int i)
        Returns m/z shift at position i
        """
        ...
    
    def getMZShiftCount(self) -> int:
        """
        Cython signature: unsigned int getMZShiftCount()
        Returns number of m/z shifts
        """
        ... 


class NASequence:
    """
    Cython implementation of _NASequence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NASequence.html>`_

    Representation of an RNA sequence
    This class represents nucleic acid sequences in OpenMS. An NASequence
    instance primarily contains a sequence of ribonucleotides.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NASequence()
        """
        ...
    
    @overload
    def __init__(self, in_0: NASequence ) -> None:
        """
        Cython signature: void NASequence(NASequence &)
        """
        ...
    
    def getSequence(self) -> List[Ribonucleotide]:
        """
        Cython signature: libcpp_vector[const Ribonucleotide *] getSequence()
        """
        ...
    
    def __getitem__(self, index: int ) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * operator[](size_t index)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Check if sequence is empty
        """
        ...
    
    def setSequence(self, seq: List[Ribonucleotide] ) -> None:
        """
        Cython signature: void setSequence(const libcpp_vector[const Ribonucleotide *] & seq)
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the peptide as string with modifications embedded in brackets
        """
        ...
    
    def setFivePrimeMod(self, modification: Ribonucleotide ) -> None:
        """
        Cython signature: void setFivePrimeMod(const Ribonucleotide * modification)
        Sets the 5' modification
        """
        ...
    
    def getFivePrimeMod(self) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * getFivePrimeMod()
        Returns the name (ID) of the N-terminal modification, or an empty string if none is set
        """
        ...
    
    def setThreePrimeMod(self, modification: Ribonucleotide ) -> None:
        """
        Cython signature: void setThreePrimeMod(const Ribonucleotide * modification)
        Sets the 3' modification
        """
        ...
    
    def getThreePrimeMod(self) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * getThreePrimeMod()
        """
        ...
    
    def get(self, index: int ) -> Ribonucleotide:
        """
        Cython signature: const Ribonucleotide * get(size_t index)
        Returns the residue at position index
        """
        ...
    
    def set(self, index: int , r: Ribonucleotide ) -> None:
        """
        Cython signature: void set(size_t index, const Ribonucleotide * r)
        Sets the residue at position index
        """
        ...
    
    @overload
    def getFormula(self, ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula()
        Returns the formula of the peptide
        """
        ...
    
    @overload
    def getFormula(self, type_: int , charge: int ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula(NASFragmentType type_, int charge)
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
        Cython signature: double getAverageWeight(NASFragmentType type_, int charge)
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
        Cython signature: double getMonoWeight(NASFragmentType type_, int charge)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Returns the number of residues
        """
        ...
    
    def getPrefix(self, length: int ) -> NASequence:
        """
        Cython signature: NASequence getPrefix(size_t length)
        Returns a peptide sequence of the first index residues
        """
        ...
    
    def getSuffix(self, length: int ) -> NASequence:
        """
        Cython signature: NASequence getSuffix(size_t length)
        Returns a peptide sequence of the last index residues
        """
        ...
    
    def getSubsequence(self, start: int , length: int ) -> NASequence:
        """
        Cython signature: NASequence getSubsequence(size_t start, size_t length)
        Returns a peptide sequence of number residues, beginning at position index
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the peptide as string with modifications embedded in brackets
        """
        ...
    
    def __richcmp__(self, other: NASequence, op: int) -> Any:
        ...
    NASFragmentType : __NASFragmentType
    
    fromString: __static_NASequence_fromString 


class PeptideProteinResolution:
    """
    Cython implementation of _PeptideProteinResolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideProteinResolution.html>`_

    Resolves shared peptides based on protein scores
    
    Resolves connected components of the bipartite protein-peptide graph based
    on protein probabilities/scores and adds them as additional protein_groups
    to the protein identification run processed.
    Thereby greedily assigns shared peptides in this component uniquely to the
    proteins of the current @em best @em indistinguishable protein group, until
    every peptide is uniquely assigned. This effectively allows more peptides to
    be used in ProteinQuantifier at the cost of potentially additional noise in
    the peptides quantities.
    In accordance with most state-of-the-art protein inference tools, only the
    best hit (PSM) for a peptide ID is considered.  Probability ties are
    currently resolved by taking the protein with larger number of peptides
    
    The class could provide iterator for ConnectedComponents in the
    future. One could extend the graph to include all PeptideHits (not only the
    best). It becomes a tripartite graph with larger connected components then.
    Maybe extend it to work with MS1 features. Separate resolution and adding
    groups to output
    """
    
    @overload
    def __init__(self, statistics: bool ) -> None:
        """
        Cython signature: void PeptideProteinResolution(bool statistics)
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideProteinResolution ) -> None:
        """
        Cython signature: void PeptideProteinResolution(PeptideProteinResolution &)
        """
        ...
    
    def buildGraph(self, protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void buildGraph(ProteinIdentification & protein, PeptideIdentificationList & peptides)
        Initialize and store the graph (= maps), needs sorted groups for
        correct functionality. Therefore sorts the indist. protein groups
        if not skipped
        
        
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: Vector of ProteinIdentifications with links to the proteins
        :param skip_sort: Skips sorting of groups, nothing is modified then
        """
        ...
    
    def resolveGraph(self, protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void resolveGraph(ProteinIdentification & protein, PeptideIdentificationList & peptides)
        Applies resolveConnectedComponent to every component of the graph and
        is able to write statistics when specified. Parameters will
        both be mutated in this method
        
        
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: vector of ProteinIdentifications with links to the proteins
        """
        ...
    
    def findConnectedComponent(self, root_prot_grp: int ) -> PeptideProteinResolution_ConnectedComponent:
        """
        Cython signature: PeptideProteinResolution_ConnectedComponent findConnectedComponent(size_t & root_prot_grp)
        Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
        and peptides), switching from one to the other in each step
        
        
        :param root_prot_grp: Starts the BFS at this protein group index
        :return: Returns a Connected Component as set of group and peptide indices
        """
        ...
    
    def resolveConnectedComponent(self, conn_comp: PeptideProteinResolution_ConnectedComponent , protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void resolveConnectedComponent(PeptideProteinResolution_ConnectedComponent & conn_comp, ProteinIdentification & protein, PeptideIdentificationList & peptides)
        Resolves connected components based on posterior probabilities and adds them
        as additional protein_groups to the output idXML.
        Thereby greedily assigns shared peptides in this component uniquely to
        the proteins of the current BEST INDISTINGUISHABLE protein group,
        ready to be used in ProteinQuantifier then.
        This is achieved by removing all other evidence from the input
        PeptideIDs and iterating until each peptide is uniquely assigned.
        In accordance with Fido only the best hit (PSM) for an ID is considered.
        Probability ties resolved by taking protein with largest number of peptides
        
        
        :param conn_comp: The component to be resolved
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: Vector of ProteinIdentifications with links to the proteins
        """
        ...
    
    run: __static_PeptideProteinResolution_run 


class PeptideProteinResolution_ConnectedComponent:
    """
    Cython implementation of _PeptideProteinResolution_ConnectedComponent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideProteinResolution_ConnectedComponent.html>`_
    """
    
    prot_grp_indices: Set[int]
    
    pep_indices: Set[int]
    
    def __init__(self) -> None:
        """
        Cython signature: void PeptideProteinResolution_ConnectedComponent()
        """
        ... 


class PercolatorOutfile:
    """
    Cython implementation of _PercolatorOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorOutfile.html>`_

    Class for reading Percolator tab-delimited output files
    
    For PSM-level output, the file extension should be ".psms"
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PercolatorOutfile()
        """
        ...
    
    @overload
    def __init__(self, in_0: PercolatorOutfile ) -> None:
        """
        Cython signature: void PercolatorOutfile(PercolatorOutfile &)
        """
        ...
    
    def getScoreType(self, score_type_name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: PercolatorOutfile_ScoreType getScoreType(String score_type_name)
        Returns a score type given its name
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , proteins: ProteinIdentification , peptides: PeptideIdentificationList , lookup: SpectrumMetaDataLookup , output_score: int ) -> None:
        """
        Cython signature: void load(const String & filename, ProteinIdentification & proteins, PeptideIdentificationList & peptides, SpectrumMetaDataLookup & lookup, PercolatorOutfile_ScoreType output_score)
        Loads a Percolator output file
        """
        ...
    PercolatorOutfile_ScoreType : __PercolatorOutfile_ScoreType 


class QcMLFile:
    """
    Cython implementation of _QcMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QcMLFile.html>`_
      -- Inherits from ['XMLHandler', 'XMLFile', 'ProgressLogger']

    File adapter for QcML files used to load and store QcML files
    
    This Class is supposed to internally collect the data for the qcML File
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void QcMLFile()
        """
        ...
    
    def exportIDstats(self, filename: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String exportIDstats(const String & filename)
        """
        ...
    
    def addRunQualityParameter(self, r: Union[bytes, str, String] , qp: QualityParameter ) -> None:
        """
        Cython signature: void addRunQualityParameter(String r, QualityParameter qp)
        Adds a QualityParameter to run by the name r
        """
        ...
    
    def addRunAttachment(self, r: Union[bytes, str, String] , at: Attachment ) -> None:
        """
        Cython signature: void addRunAttachment(String r, Attachment at)
        Adds a attachment to run by the name r
        """
        ...
    
    def addSetQualityParameter(self, r: Union[bytes, str, String] , qp: QualityParameter ) -> None:
        """
        Cython signature: void addSetQualityParameter(String r, QualityParameter qp)
        Adds a QualityParameter to set by the name r
        """
        ...
    
    def addSetAttachment(self, r: Union[bytes, str, String] , at: Attachment ) -> None:
        """
        Cython signature: void addSetAttachment(String r, Attachment at)
        Adds a attachment to set by the name r
        """
        ...
    
    @overload
    def removeAttachment(self, r: Union[bytes, str, String] , ids: List[bytes] , at: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeAttachment(String r, libcpp_vector[String] & ids, String at)
        Removes attachments referencing an id given in ids, from run/set r. All attachments if no attachment name is given with at
        """
        ...
    
    @overload
    def removeAttachment(self, r: Union[bytes, str, String] , at: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeAttachment(String r, String at)
        Removes attachment with cv accession at from run/set r
        """
        ...
    
    def removeAllAttachments(self, at: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeAllAttachments(String at)
        Removes attachment with cv accession at from all runs/sets
        """
        ...
    
    def removeQualityParameter(self, r: Union[bytes, str, String] , ids: List[bytes] ) -> None:
        """
        Cython signature: void removeQualityParameter(String r, libcpp_vector[String] & ids)
        Removes QualityParameter going by one of the ID attributes given in ids
        """
        ...
    
    def merge(self, addendum: QcMLFile , setname: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void merge(QcMLFile & addendum, String setname)
        Merges the given QCFile into this one
        """
        ...
    
    def collectSetParameter(self, setname: Union[bytes, str, String] , qp: Union[bytes, str, String] , ret: List[bytes] ) -> None:
        """
        Cython signature: void collectSetParameter(String setname, String qp, libcpp_vector[String] & ret)
        Collects the values of given QPs (as CVid) of the given set
        """
        ...
    
    def exportAttachment(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String exportAttachment(String filename, String qpname)
        Returns a String of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname
        """
        ...
    
    def getRunNames(self, ids: List[bytes] ) -> None:
        """
        Cython signature: void getRunNames(libcpp_vector[String] & ids)
        Gives the names of the registered runs in the vector ids
        """
        ...
    
    def existsRun(self, filename: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool existsRun(String filename)
        Returns true if the given run id is present in this file, if checkname is true it also checks the names
        """
        ...
    
    def existsSet(self, filename: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool existsSet(String filename)
        Returns true if the given set id is present in this file, if checkname is true it also checks the names
        """
        ...
    
    def existsRunQualityParameter(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] , ids: List[bytes] ) -> None:
        """
        Cython signature: void existsRunQualityParameter(String filename, String qpname, libcpp_vector[String] & ids)
        Returns the ids of the parameter name given if found in given run empty else
        """
        ...
    
    def existsSetQualityParameter(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] , ids: List[bytes] ) -> None:
        """
        Cython signature: void existsSetQualityParameter(String filename, String qpname, libcpp_vector[String] & ids)
        Returns the ids of the parameter name given if found in given set, empty else
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & filename)
        Store the qcML file
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void load(const String & filename)
        Load a QCFile
        """
        ...
    
    def registerRun(self, id_: Union[bytes, str, String] , name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void registerRun(String id_, String name)
        Registers a run in the qcml file with the respective mappings
        """
        ...
    
    def registerSet(self, id_: Union[bytes, str, String] , name: Union[bytes, str, String] , names: Set[bytes] ) -> None:
        """
        Cython signature: void registerSet(String id_, String name, libcpp_set[String] & names)
        Registers a set in the qcml file with the respective mappings
        """
        ...
    
    def exportQP(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String exportQP(String filename, String qpname)
        Returns a String value in quotation of a QualityParameter by the name qpname in run/set by the name filename
        """
        ...
    
    def exportQPs(self, filename: Union[bytes, str, String] , qpnames: List[bytes] ) -> Union[bytes, str, String]:
        """
        Cython signature: String exportQPs(String filename, StringList qpnames)
        Returns a String of a tab separated QualityParameter by the name qpname in run/set by the name filename
        """
        ...
    
    def getRunIDs(self, ids: List[bytes] ) -> None:
        """
        Cython signature: void getRunIDs(libcpp_vector[String] & ids)
        Gives the ids of the registered runs in the vector ids
        """
        ...
    
    def reset(self) -> None:
        """
        Cython signature: void reset()
        """
        ...
    
    def error(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None:
        """
        Cython signature: void error(ActionMode mode, const String & msg, unsigned int line, unsigned int column)
        """
        ...
    
    def warning(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None:
        """
        Cython signature: void warning(ActionMode mode, const String & msg, unsigned int line, unsigned int column)
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
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


class QualityParameter:
    """
    Cython implementation of _QualityParameter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QualityParameter.html>`_
    """
    
    name: Union[bytes, str, String]
    
    id: Union[bytes, str, String]
    
    value: Union[bytes, str, String]
    
    cvRef: Union[bytes, str, String]
    
    cvAcc: Union[bytes, str, String]
    
    unitRef: Union[bytes, str, String]
    
    unitAcc: Union[bytes, str, String]
    
    flag: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void QualityParameter()
        """
        ...
    
    @overload
    def __init__(self, in_0: QualityParameter ) -> None:
        """
        Cython signature: void QualityParameter(QualityParameter &)
        """
        ...
    
    def toXMLString(self, indentation_level: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String toXMLString(unsigned int indentation_level)
        """
        ...
    
    def __richcmp__(self, other: QualityParameter, op: int) -> Any:
        ... 


class SiriusFragmentAnnotation:
    """
    Cython implementation of _SiriusFragmentAnnotation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusFragmentAnnotation.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusFragmentAnnotation()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusFragmentAnnotation ) -> None:
        """
        Cython signature: void SiriusFragmentAnnotation(SiriusFragmentAnnotation &)
        """
        ...
    
    def extractAnnotationsFromSiriusFile(self, path_to_sirius_workspace: String , max_rank: int , decoy: bool , use_exact_mass: bool ) -> List[MSSpectrum]:
        """
        Cython signature: libcpp_vector[MSSpectrum] extractAnnotationsFromSiriusFile(String & path_to_sirius_workspace, size_t max_rank, bool decoy, bool use_exact_mass)
        """
        ...
    
    def extractAndResolveSiriusAnnotations(self, sirius_workspace_subdirs: List[bytes] , score_threshold: float , use_exact_mass: bool , decoy_generation: bool ) -> List[SiriusFragmentAnnotation_SiriusTargetDecoySpectra]:
        """
        Cython signature: libcpp_vector[SiriusFragmentAnnotation_SiriusTargetDecoySpectra] extractAndResolveSiriusAnnotations(libcpp_vector[String] & sirius_workspace_subdirs, double score_threshold, bool use_exact_mass, bool decoy_generation)
        """
        ...
    
    def extract_columnname_to_columnindex(self, csvfile: CsvFile ) -> Dict[bytes, int]:
        """
        Cython signature: libcpp_map[libcpp_string,size_t] extract_columnname_to_columnindex(CsvFile & csvfile)
        """
        ... 


class SiriusFragmentAnnotation_SiriusTargetDecoySpectra:
    """
    Cython implementation of _SiriusFragmentAnnotation_SiriusTargetDecoySpectra

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusFragmentAnnotation_SiriusTargetDecoySpectra.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusFragmentAnnotation_SiriusTargetDecoySpectra()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusFragmentAnnotation_SiriusTargetDecoySpectra ) -> None:
        """
        Cython signature: void SiriusFragmentAnnotation_SiriusTargetDecoySpectra(SiriusFragmentAnnotation_SiriusTargetDecoySpectra &)
        """
        ... 


class SplinePackage:
    """
    Cython implementation of _SplinePackage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplinePackage.html>`_
    """
    
    @overload
    def __init__(self, pos: List[float] , intensity: List[float] ) -> None:
        """
        Cython signature: void SplinePackage(libcpp_vector[double] pos, libcpp_vector[double] intensity)
        """
        ...
    
    @overload
    def __init__(self, in_0: SplinePackage ) -> None:
        """
        Cython signature: void SplinePackage(SplinePackage &)
        """
        ...
    
    def getPosMin(self) -> float:
        """
        Cython signature: double getPosMin()
        Returns the minimum position for which the spline fit is valid
        """
        ...
    
    def getPosMax(self) -> float:
        """
        Cython signature: double getPosMax()
        Returns the maximum position for which the spline fit is valid
        """
        ...
    
    def getPosStepWidth(self) -> float:
        """
        Cython signature: double getPosStepWidth()
        Returns a sensible position step width for the package
        """
        ...
    
    def isInPackage(self, pos: float ) -> bool:
        """
        Cython signature: bool isInPackage(double pos)
        Returns true if position in
        """
        ...
    
    def eval(self, pos: float ) -> float:
        """
        Cython signature: double eval(double pos)
        Returns interpolated intensity position `pos`
        """
        ... 


class _Interfaces_BinaryDataArray:
    """
    Cython implementation of _BinaryDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1BinaryDataArray.html>`_
    """
    
    data: List[float]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void _Interfaces_BinaryDataArray()
        """
        ...
    
    @overload
    def __init__(self, in_0: _Interfaces_BinaryDataArray ) -> None:
        """
        Cython signature: void _Interfaces_BinaryDataArray(_Interfaces_BinaryDataArray &)
        """
        ... 


class _Interfaces_Chromatogram:
    """
    Cython implementation of _Chromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1Chromatogram.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void _Interfaces_Chromatogram()
        """
        ...
    
    @overload
    def __init__(self, in_0: _Interfaces_Chromatogram ) -> None:
        """
        Cython signature: void _Interfaces_Chromatogram(_Interfaces_Chromatogram &)
        """
        ... 


class _Interfaces_Spectrum:
    """
    Cython implementation of _Spectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1Spectrum.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void _Interfaces_Spectrum()
        """
        ...
    
    @overload
    def __init__(self, in_0: _Interfaces_Spectrum ) -> None:
        """
        Cython signature: void _Interfaces_Spectrum(_Interfaces_Spectrum &)
        """
        ... 


class __NASFragmentType:
    None
    Full : int
    Internal : int
    FivePrime : int
    ThreePrime : int
    AIon : int
    BIon : int
    CIon : int
    XIon : int
    YIon : int
    ZIon : int
    Precursor : int
    BIonMinusH20 : int
    YIonMinusH20 : int
    BIonMinusNH3 : int
    YIonMinusNH3 : int
    NonIdentified : int
    Unannotated : int
    WIon : int
    AminusB : int
    DIon : int
    SizeOfNASFragmentType : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __PercolatorOutfile_ScoreType:
    None
    QVALUE : int
    POSTERRPROB : int
    SCORE : int
    SIZE_OF_SCORETYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

