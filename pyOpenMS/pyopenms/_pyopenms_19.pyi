from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_MRMRTNormalizer_chauvenet(residuals: List[float] , pos: int ) -> bool:
    """
    Cython signature: bool chauvenet(libcpp_vector[double] residuals, int pos)
    """
    ...

def __static_MRMRTNormalizer_chauvenet_probability(residuals: List[float] , pos: int ) -> float:
    """
    Cython signature: double chauvenet_probability(libcpp_vector[double] residuals, int pos)
    """
    ...

def __static_MRMRTNormalizer_computeBinnedCoverage(rtRange: List[float, float] , pairs: List[List[float, float]] , nrBins: int , minPeptidesPerBin: int , minBinsFilled: int ) -> bool:
    """
    Cython signature: bool computeBinnedCoverage(libcpp_pair[double,double] rtRange, libcpp_vector[libcpp_pair[double,double]] & pairs, int nrBins, int minPeptidesPerBin, int minBinsFilled)
    """
    ...

def __static_MRMRTNormalizer_removeOutliersIterative(pairs: List[List[float, float]] , rsq_limit: float , coverage_limit: float , use_chauvenet: bool , outlier_detection_method: bytes ) -> List[List[float, float]]:
    """
    Cython signature: libcpp_vector[libcpp_pair[double,double]] removeOutliersIterative(libcpp_vector[libcpp_pair[double,double]] & pairs, double rsq_limit, double coverage_limit, bool use_chauvenet, libcpp_string outlier_detection_method)
    """
    ...

def __static_MRMRTNormalizer_removeOutliersRANSAC(pairs: List[List[float, float]] , rsq_limit: float , coverage_limit: float , max_iterations: int , max_rt_threshold: float , sampling_size: int ) -> List[List[float, float]]:
    """
    Cython signature: libcpp_vector[libcpp_pair[double,double]] removeOutliersRANSAC(libcpp_vector[libcpp_pair[double,double]] & pairs, double rsq_limit, double coverage_limit, size_t max_iterations, double max_rt_threshold, size_t sampling_size)
    """
    ...


class AnnotatedMSRun:
    """
    Cython implementation of _AnnotatedMSRun

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AnnotatedMSRun.html>`_

    Class for storing MS run data with peptide and protein identifications
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AnnotatedMSRun()
        """
        ...
    
    @overload
    def __init__(self, in_0: AnnotatedMSRun ) -> None:
        """
        Cython signature: void AnnotatedMSRun(AnnotatedMSRun)
        """
        ...
    
    def getProteinIdentifications(self) -> List[ProteinIdentification]:
        """
        Cython signature: libcpp_vector[ProteinIdentification] getProteinIdentifications()
        """
        ...
    
    def setProteinIdentifications(self, ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void setProteinIdentifications(libcpp_vector[ProteinIdentification] & ids)
        """
        ...
    
    def getPeptideIdentifications(self) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList getPeptideIdentifications()
        """
        ...
    
    def setPeptideIdentifications(self, ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setPeptideIdentifications(PeptideIdentificationList ids)
        """
        ...
    
    def getMSExperiment(self) -> MSExperiment:
        """
        Cython signature: MSExperiment getMSExperiment()
        """
        ...
    
    def setMSExperiment(self, experiment: MSExperiment ) -> None:
        """
        Cython signature: void setMSExperiment(MSExperiment & experiment)
        """
        ... 


class BinnedSpectrum:
    """
    Cython implementation of _BinnedSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BinnedSpectrum.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BinnedSpectrum()
        """
        ...
    
    @overload
    def __init__(self, in_0: BinnedSpectrum ) -> None:
        """
        Cython signature: void BinnedSpectrum(BinnedSpectrum &)
        """
        ...
    
    @overload
    def __init__(self, in_0: MSSpectrum , size: float , unit_ppm: bool , spread: int , offset: float ) -> None:
        """
        Cython signature: void BinnedSpectrum(MSSpectrum, float size, bool unit_ppm, unsigned int spread, float offset)
        """
        ...
    
    def getBinSize(self) -> float:
        """
        Cython signature: float getBinSize()
        Returns the bin size
        """
        ...
    
    def getBinSpread(self) -> int:
        """
        Cython signature: unsigned int getBinSpread()
        Returns the bin spread
        """
        ...
    
    def getBinIndex(self, mz: float ) -> int:
        """
        Cython signature: unsigned int getBinIndex(float mz)
        Returns the bin index of a given m/z position
        """
        ...
    
    def getBinLowerMZ(self, i: int ) -> float:
        """
        Cython signature: float getBinLowerMZ(size_t i)
        Returns the lower m/z of a bin given its index
        """
        ...
    
    def getBinIntensity(self, mz: float ) -> float:
        """
        Cython signature: float getBinIntensity(double mz)
        Returns the bin intensity at a given m/z position
        """
        ...
    
    def getPrecursors(self) -> List[Precursor]:
        """
        Cython signature: libcpp_vector[Precursor] getPrecursors()
        """
        ...
    
    def isCompatible(self, a: BinnedSpectrum , b: BinnedSpectrum ) -> bool:
        """
        Cython signature: bool isCompatible(BinnedSpectrum & a, BinnedSpectrum & b)
        """
        ...
    
    def getOffset(self) -> float:
        """
        Cython signature: float getOffset()
        Returns offset
        """
        ...
    
    def __richcmp__(self, other: BinnedSpectrum, op: int) -> Any:
        ... 


class ChromatogramPeak:
    """
    Cython implementation of _ChromatogramPeak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ChromatogramPeak_1_1ChromatogramPeak.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramPeak()
        A 1-dimensional raw data point or peak for chromatograms
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramPeak ) -> None:
        """
        Cython signature: void ChromatogramPeak(ChromatogramPeak &)
        """
        ...
    
    @overload
    def __init__(self, retention_time: DPosition1 , intensity: float ) -> None:
        """
        Cython signature: void ChromatogramPeak(DPosition1 retention_time, double intensity)
        """
        ...
    
    def getIntensity(self) -> float:
        """
        Cython signature: double getIntensity()
        Returns the intensity
        """
        ...
    
    def setIntensity(self, in_0: float ) -> None:
        """
        Cython signature: void setIntensity(double)
        Sets the intensity
        """
        ...
    
    def getPosition(self) -> DPosition1:
        """
        Cython signature: DPosition1 getPosition()
        """
        ...
    
    def setPosition(self, in_0: DPosition1 ) -> None:
        """
        Cython signature: void setPosition(DPosition1)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        Returns the retention time
        """
        ...
    
    def setRT(self, in_0: float ) -> None:
        """
        Cython signature: void setRT(double)
        Sets retention time
        """
        ...
    
    def getPos(self) -> float:
        """
        Cython signature: double getPos()
        Alias for getRT()
        """
        ...
    
    def setPos(self, in_0: float ) -> None:
        """
        Cython signature: void setPos(double)
        Alias for setRT()
        """
        ...
    
    def __richcmp__(self, other: ChromatogramPeak, op: int) -> Any:
        ... 


class ChromeleonFile:
    """
    Cython implementation of _ChromeleonFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromeleonFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromeleonFile()
        Load Chromeleon HPLC text file and save it into a `MSExperiment`.
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromeleonFile ) -> None:
        """
        Cython signature: void ChromeleonFile(ChromeleonFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , experiment: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & experiment)
        Load the file's data and metadata, and save it into a `MSExperiment`
        """
        ... 


class FalseDiscoveryRate:
    """
    Cython implementation of _FalseDiscoveryRate

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FalseDiscoveryRate.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FalseDiscoveryRate()
        """
        ...
    
    @overload
    def apply(self, forward_ids: PeptideIdentificationList , reverse_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void apply(PeptideIdentificationList & forward_ids, PeptideIdentificationList & reverse_ids)
        """
        ...
    
    @overload
    def apply(self, id: PeptideIdentificationList ) -> None:
        """
        Cython signature: void apply(PeptideIdentificationList & id)
        """
        ...
    
    @overload
    def apply(self, forward_ids: List[ProteinIdentification] , reverse_ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void apply(libcpp_vector[ProteinIdentification] & forward_ids, libcpp_vector[ProteinIdentification] & reverse_ids)
        """
        ...
    
    @overload
    def apply(self, id: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void apply(libcpp_vector[ProteinIdentification] & id)
        """
        ...
    
    def applyEstimated(self, ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void applyEstimated(libcpp_vector[ProteinIdentification] & ids)
        """
        ...
    
    @overload
    def applyEvaluateProteinIDs(self, ids: List[ProteinIdentification] , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float:
        """
        Cython signature: double applyEvaluateProteinIDs(libcpp_vector[ProteinIdentification] & ids, double pepCutoff, unsigned int fpCutoff, double diffWeight)
        """
        ...
    
    @overload
    def applyEvaluateProteinIDs(self, ids: ProteinIdentification , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float:
        """
        Cython signature: double applyEvaluateProteinIDs(ProteinIdentification & ids, double pepCutoff, unsigned int fpCutoff, double diffWeight)
        """
        ...
    
    @overload
    def applyBasic(self, run_info: List[ProteinIdentification] , ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void applyBasic(libcpp_vector[ProteinIdentification] & run_info, PeptideIdentificationList & ids)
        """
        ...
    
    @overload
    def applyBasic(self, ids: PeptideIdentificationList , higher_score_better: bool , charge: int , identifier: Union[bytes, str, String] , only_best_per_pep: bool ) -> None:
        """
        Cython signature: void applyBasic(PeptideIdentificationList & ids, bool higher_score_better, int charge, String identifier, bool only_best_per_pep)
        """
        ...
    
    @overload
    def applyBasic(self, cmap: ConsensusMap , use_unassigned_peptides: bool ) -> None:
        """
        Cython signature: void applyBasic(ConsensusMap & cmap, bool use_unassigned_peptides)
        """
        ...
    
    @overload
    def applyBasic(self, id: ProteinIdentification , groups_too: bool ) -> None:
        """
        Cython signature: void applyBasic(ProteinIdentification & id, bool groups_too)
        """
        ...
    
    def applyPickedProteinFDR(self, id: ProteinIdentification , decoy_string: String , decoy_prefix: bool , groups_too: bool ) -> None:
        """
        Cython signature: void applyPickedProteinFDR(ProteinIdentification & id, String & decoy_string, bool decoy_prefix, bool groups_too)
        """
        ...
    
    @overload
    def rocN(self, ids: PeptideIdentificationList , fp_cutoff: int ) -> float:
        """
        Cython signature: double rocN(PeptideIdentificationList & ids, size_t fp_cutoff)
        """
        ...
    
    @overload
    def rocN(self, ids: ConsensusMap , fp_cutoff: int , include_unassigned_peptides: bool ) -> float:
        """
        Cython signature: double rocN(ConsensusMap & ids, size_t fp_cutoff, bool include_unassigned_peptides)
        """
        ...
    
    @overload
    def rocN(self, ids: ConsensusMap , fp_cutoff: int , identifier: Union[bytes, str, String] , include_unassigned_peptides: bool ) -> float:
        """
        Cython signature: double rocN(ConsensusMap & ids, size_t fp_cutoff, const String & identifier, bool include_unassigned_peptides)
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


class ILPDCWrapper:
    """
    Cython implementation of _ILPDCWrapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ILPDCWrapper.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ILPDCWrapper()
        """
        ...
    
    @overload
    def __init__(self, in_0: ILPDCWrapper ) -> None:
        """
        Cython signature: void ILPDCWrapper(ILPDCWrapper &)
        """
        ...
    
    def compute(self, fm: FeatureMap , pairs: List[ChargePair] , verbose_level: int ) -> float:
        """
        Cython signature: double compute(FeatureMap fm, libcpp_vector[ChargePair] & pairs, size_t verbose_level)
        Compute optimal solution and return value of objective function. If the input feature map is empty, a warning is issued and -1 is returned
        """
        ... 


class ItraqEightPlexQuantitationMethod:
    """
    Cython implementation of _ItraqEightPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqEightPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ItraqEightPlexQuantitationMethod()
        iTRAQ 8 plex quantitation to be used with the IsobaricQuantitation
        """
        ...
    
    @overload
    def __init__(self, in_0: ItraqEightPlexQuantitationMethod ) -> None:
        """
        Cython signature: void ItraqEightPlexQuantitationMethod(ItraqEightPlexQuantitationMethod &)
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


class MRMRTNormalizer:
    """
    Cython implementation of _MRMRTNormalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMRTNormalizer.html>`_
    """
    
    chauvenet: __static_MRMRTNormalizer_chauvenet
    
    chauvenet_probability: __static_MRMRTNormalizer_chauvenet_probability
    
    computeBinnedCoverage: __static_MRMRTNormalizer_computeBinnedCoverage
    
    removeOutliersIterative: __static_MRMRTNormalizer_removeOutliersIterative
    
    removeOutliersRANSAC: __static_MRMRTNormalizer_removeOutliersRANSAC 


class ModificationsDB:
    """
    Cython implementation of _ModificationsDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationsDB.html>`_
    """
    
    def getNumberOfModifications(self) -> int:
        """
        Cython signature: size_t getNumberOfModifications()
        Returns the number of modifications read from the unimod.xml file
        """
        ...
    
    def searchModifications(self, mods: Set[ResidueModification] , mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> None:
        """
        Cython signature: void searchModifications(libcpp_set[const ResidueModification *] & mods, const String & mod_name, const String & residue, TermSpecificity term_spec)
        Collects all modifications which have the given name as synonym
        
        If `residue` is set, only modifications with matching residue of origin are considered
        If `term_spec` is set, only modifications with matching term specificity are considered
        The resulting set of modifications will be empty if no modification exists that fulfills the criteria
        """
        ...
    
    @overload
    def getModification(self, index: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(size_t index)
        Returns the modification with the given index
        """
        ...
    
    @overload
    def getModification(self, mod_name: Union[bytes, str, String] ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(const String & mod_name)
        Returns the modification with the given name
        """
        ...
    
    @overload
    def getModification(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification(const String & mod_name, const String & residue, TermSpecificity term_spec)
        Returns the modification with the given arguments
        """
        ...
    
    def has(self, modification: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool has(String modification)
        Returns true if the modification exists
        """
        ...
    
    def addModification(self, new_mod: ResidueModification ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * addModification(const ResidueModification & new_mod)
        Add a new modification to ModificationsDB. If the modification already exists (based on its fullID) it is not added. Returns the modification in the ModificationDB (which can differ from input if mod was already present).
        """
        ...
    
    def findModificationIndex(self, mod_name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t findModificationIndex(const String & mod_name)
        Returns the index of the modification in the mods_ vector; a unique name must be given
        """
        ...
    
    def searchModificationsByDiffMonoMass(self, mods: List[bytes] , mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> None:
        """
        Cython signature: void searchModificationsByDiffMonoMass(libcpp_vector[String] & mods, double mass, double max_error, const String & residue, TermSpecificity term_spec)
        Collects all modifications with delta mass inside a tolerance window
        """
        ...
    
    def getBestModificationByDiffMonoMass(self, mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getBestModificationByDiffMonoMass(double mass, double max_error, const String & residue, TermSpecificity term_spec)
        Returns the best matching modification for the given delta mass and residue
        
        Query the modifications DB to get the best matching modification with
        the given delta mass at the given residue (NULL pointer means no result,
        maybe the maximal error tolerance needs to be increased). Possible
        input for CAM modification would be a delta mass of 57 and a residue
        of "C".
        
        Note: If there are multiple possible matches with equal masses, it
        will choose the _first_ match which defaults to the first matching
        UniMod entry.
        
        
        :param residue: The residue at which the modifications occurs
        :param mass: The monoisotopic mass of the residue including the mass of the modification
        :param max_error: The maximal mass error in the modification search
        :return: A pointer to the best matching modification (or NULL if none was found)
        """
        ...
    
    def getAllSearchModifications(self, modifications: List[bytes] ) -> None:
        """
        Cython signature: void getAllSearchModifications(libcpp_vector[String] & modifications)
        Collects all modifications that can be used for identification searches
        """
        ...
    
    def isInstantiated(self) -> bool:
        """
        Cython signature: bool isInstantiated()
        Check whether ModificationsDB was instantiated before
        """
        ... 


class ProteinProteinCrossLink:
    """
    Cython implementation of _ProteinProteinCrossLink

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1ProteinProteinCrossLink.html>`_
    """
    
    alpha: AASequence
    
    beta: AASequence
    
    cross_link_position: List[int, int]
    
    cross_linker_mass: float
    
    cross_linker_name: Union[bytes, str, String]
    
    term_spec_alpha: int
    
    term_spec_beta: int
    
    precursor_correction: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteinProteinCrossLink()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteinProteinCrossLink ) -> None:
        """
        Cython signature: void ProteinProteinCrossLink(ProteinProteinCrossLink &)
        """
        ...
    
    def getType(self) -> int:
        """
        Cython signature: ProteinProteinCrossLinkType getType()
        """
        ...
    
    def __richcmp__(self, other: ProteinProteinCrossLink, op: int) -> Any:
        ... 


class RansacModelQuadratic:
    """
    Cython implementation of _RansacModelQuadratic

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RansacModelQuadratic.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RansacModelQuadratic()
        """
        ...
    
    @overload
    def __init__(self, in_0: RansacModelQuadratic ) -> None:
        """
        Cython signature: void RansacModelQuadratic(RansacModelQuadratic &)
        """
        ... 


class Sample:
    """
    Cython implementation of _Sample

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Sample.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Sample()
        """
        ...
    
    @overload
    def __init__(self, in_0: Sample ) -> None:
        """
        Cython signature: void Sample(Sample &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        """
        ...
    
    def getOrganism(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOrganism()
        """
        ...
    
    def setOrganism(self, organism: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setOrganism(String organism)
        """
        ...
    
    def getNumber(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNumber()
        Returns the sample number
        """
        ...
    
    def setNumber(self, number: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNumber(String number)
        Sets the sample number (e.g. sample ID)
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the comment (default "")
        """
        ...
    
    def setComment(self, comment: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String comment)
        Sets the comment (may contain newline characters)
        """
        ...
    
    def getState(self) -> int:
        """
        Cython signature: SampleState getState()
        Returns the state of aggregation (default SAMPLENULL)
        """
        ...
    
    def setState(self, state: int ) -> None:
        """
        Cython signature: void setState(SampleState state)
        Sets the state of aggregation
        """
        ...
    
    def getMass(self) -> float:
        """
        Cython signature: double getMass()
        Returns the mass (in gram) (default 0.0)
        """
        ...
    
    def setMass(self, mass: float ) -> None:
        """
        Cython signature: void setMass(double mass)
        Sets the mass (in gram)
        """
        ...
    
    def getVolume(self) -> float:
        """
        Cython signature: double getVolume()
        Returns the volume (in ml) (default 0.0)
        """
        ...
    
    def setVolume(self, volume: float ) -> None:
        """
        Cython signature: void setVolume(double volume)
        Sets the volume (in ml)
        """
        ...
    
    def getConcentration(self) -> float:
        """
        Cython signature: double getConcentration()
        Returns the concentration (in g/l) (default 0.0)
        """
        ...
    
    def setConcentration(self, concentration: float ) -> None:
        """
        Cython signature: void setConcentration(double concentration)
        Sets the concentration (in g/l)
        """
        ...
    
    def getSubsamples(self) -> List[Sample]:
        """
        Cython signature: libcpp_vector[Sample] getSubsamples()
        Returns a reference to the vector of subsamples that were combined to create this sample
        """
        ...
    
    def setSubsamples(self, subsamples: List[Sample] ) -> None:
        """
        Cython signature: void setSubsamples(libcpp_vector[Sample] subsamples)
        Sets the vector of subsamples that were combined to create this sample
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
    
    def __richcmp__(self, other: Sample, op: int) -> Any:
        ...
    SampleState : __SampleState 


class ScanWindow:
    """
    Cython implementation of _ScanWindow

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ScanWindow.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    begin: float
    
    end: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ScanWindow()
        """
        ...
    
    @overload
    def __init__(self, in_0: ScanWindow ) -> None:
        """
        Cython signature: void ScanWindow(ScanWindow &)
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
    
    def __richcmp__(self, other: ScanWindow, op: int) -> Any:
        ... 


class TMTSixteenPlexQuantitationMethod:
    """
    Cython implementation of _TMTSixteenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTSixteenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TMTSixteenPlexQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: TMTSixteenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void TMTSixteenPlexQuantitationMethod(TMTSixteenPlexQuantitationMethod &)
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


class __SampleState:
    None
    SAMPLENULL : int
    SOLID : int
    LIQUID : int
    GAS : int
    SOLUTION : int
    EMULSION : int
    SUSPENSION : int
    SIZE_OF_SAMPLESTATE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

