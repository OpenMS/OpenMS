from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AccurateMassSearchResult:
    """
    Cython implementation of _AccurateMassSearchResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AccurateMassSearchResult.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void AccurateMassSearchResult()
        """
        ...
    
    def getObservedMZ(self) -> float:
        """
        Cython signature: double getObservedMZ()
        """
        ...
    
    def setObservedMZ(self, m: float ) -> None:
        """
        Cython signature: void setObservedMZ(double & m)
        """
        ...
    
    def getCalculatedMZ(self) -> float:
        """
        Cython signature: double getCalculatedMZ()
        """
        ...
    
    def setCalculatedMZ(self, m: float ) -> None:
        """
        Cython signature: void setCalculatedMZ(double & m)
        """
        ...
    
    def getQueryMass(self) -> float:
        """
        Cython signature: double getQueryMass()
        """
        ...
    
    def setQueryMass(self, m: float ) -> None:
        """
        Cython signature: void setQueryMass(double & m)
        """
        ...
    
    def getFoundMass(self) -> float:
        """
        Cython signature: double getFoundMass()
        """
        ...
    
    def setFoundMass(self, m: float ) -> None:
        """
        Cython signature: void setFoundMass(double & m)
        """
        ...
    
    def getCharge(self) -> float:
        """
        Cython signature: double getCharge()
        """
        ...
    
    def setCharge(self, ch: float ) -> None:
        """
        Cython signature: void setCharge(double & ch)
        """
        ...
    
    def getMZErrorPPM(self) -> float:
        """
        Cython signature: double getMZErrorPPM()
        """
        ...
    
    def setMZErrorPPM(self, ppm: float ) -> None:
        """
        Cython signature: void setMZErrorPPM(double & ppm)
        """
        ...
    
    def getObservedRT(self) -> float:
        """
        Cython signature: double getObservedRT()
        """
        ...
    
    def setObservedRT(self, rt: float ) -> None:
        """
        Cython signature: void setObservedRT(double & rt)
        """
        ...
    
    def getObservedIntensity(self) -> float:
        """
        Cython signature: double getObservedIntensity()
        """
        ...
    
    def setObservedIntensity(self, intensity: float ) -> None:
        """
        Cython signature: void setObservedIntensity(double & intensity)
        """
        ...
    
    def getMatchingIndex(self) -> float:
        """
        Cython signature: double getMatchingIndex()
        """
        ...
    
    def setMatchingIndex(self, idx: float ) -> None:
        """
        Cython signature: void setMatchingIndex(double & idx)
        """
        ...
    
    def getFoundAdduct(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFoundAdduct()
        """
        ...
    
    def setFoundAdduct(self, add: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFoundAdduct(const String & add)
        """
        ...
    
    def getFormulaString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFormulaString()
        """
        ...
    
    def setEmpiricalFormula(self, ep: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setEmpiricalFormula(const String & ep)
        """
        ...
    
    def getMatchingHMDBids(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getMatchingHMDBids()
        """
        ...
    
    def setMatchingHMDBids(self, match_ids: List[bytes] ) -> None:
        """
        Cython signature: void setMatchingHMDBids(libcpp_vector[String] & match_ids)
        """
        ...
    
    def getIsotopesSimScore(self) -> float:
        """
        Cython signature: double getIsotopesSimScore()
        """
        ...
    
    def setIsotopesSimScore(self, sim_score: float ) -> None:
        """
        Cython signature: void setIsotopesSimScore(double & sim_score)
        """
        ...
    
    def getIndividualIntensities(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getIndividualIntensities()
        """
        ...
    
    def setIndividualIntensities(self, in_0: List[float] ) -> None:
        """
        Cython signature: void setIndividualIntensities(libcpp_vector[double])
        """
        ...
    
    def getSourceFeatureIndex(self) -> int:
        """
        Cython signature: size_t getSourceFeatureIndex()
        """
        ...
    
    def setSourceFeatureIndex(self, in_0: int ) -> None:
        """
        Cython signature: void setSourceFeatureIndex(size_t)
        """
        ...
    
    def getMasstraceIntensities(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getMasstraceIntensities()
        """
        ...
    
    def setMasstraceIntensities(self, in_0: List[float] ) -> None:
        """
        Cython signature: void setMasstraceIntensities(libcpp_vector[double] &)
        """
        ... 


class BasicProteinInferenceAlgorithm:
    """
    Cython implementation of _BasicProteinInferenceAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BasicProteinInferenceAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Algorithm class that implements simple protein inference by aggregation of peptide scores.
    
    It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
    and if you want to use shared peptides ("use_shared_peptides").
    First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
    Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
    - Modifications (modified sequence string)
    - Charge states
    The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
    Possible aggregation methods that can be set via the parameter "aggregation_method" are:
    - "best" (default)
    - "sum"
    - "product" (ignoring zeroes)
    Annotation of the number of peptides used for aggregation can be disabled (see parameters).
    Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
    Warning: Does not "link" the peptides to the resulting protein run. If you wish to do that you have to do
    it manually.
    
    Usage:
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void BasicProteinInferenceAlgorithm()
        """
        ...
    
    @overload
    def run(self, pep_ids: PeptideIdentificationList , prot_ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void run(PeptideIdentificationList & pep_ids, libcpp_vector[ProteinIdentification] & prot_ids)
        Performs basic aggregation-based inference per ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_ids: Vector of protein identification runs. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        """
        ...
    
    @overload
    def run(self, pep_ids: PeptideIdentificationList , prot_id: ProteinIdentification ) -> None:
        """
        Cython signature: void run(PeptideIdentificationList & pep_ids, ProteinIdentification & prot_id)
        Performs basic aggregation-based inference on single ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        """
        ...
    
    @overload
    def run(self, cmap: ConsensusMap , prot_id: ProteinIdentification , include_unassigned: bool ) -> None:
        """
        Cython signature: void run(ConsensusMap & cmap, ProteinIdentification & prot_id, bool include_unassigned)
        Performs basic aggregation-based inference on identifications in a ConsensusMap. See class help.\n
        `prot_id` should contain the union of all proteins in the map. E.g. use ConsensusMapMergerAlgorithm and
        then pass the first=merged run.
        
        
        :param cmap: ConsensusMap = Consensus features with metadata and peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
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
    AggregationMethod : __AggregationMethod 


class GNPSMGFFile:
    """
    Cython implementation of _GNPSMGFFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSMGFFile.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GNPSMGFFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: GNPSMGFFile ) -> None:
        """
        Cython signature: void GNPSMGFFile(GNPSMGFFile &)
        """
        ...
    
    def store(self, consensus_file_path: Union[bytes, str, String] , mzml_file_paths: List[bytes] , out: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & consensus_file_path, const StringList & mzml_file_paths, const String & out)
        Export consensus file from default workflow to GNPS MGF format
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


class InspectInfile:
    """
    Cython implementation of _InspectInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InspectInfile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void InspectInfile()
        Inspect input file adapter
        """
        ...
    
    @overload
    def __init__(self, in_0: InspectInfile ) -> None:
        """
        Cython signature: void InspectInfile(InspectInfile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & filename)
        Stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
        """
        ...
    
    def handlePTMs(self, modification_line: Union[bytes, str, String] , modifications_filename: Union[bytes, str, String] , monoisotopic: bool ) -> None:
        """
        Cython signature: void handlePTMs(const String & modification_line, const String & modifications_filename, bool monoisotopic)
        Retrieves the name, mass change, affected residues, type and position for all modifications from a string
        
        
        :param modification_line:
        :param modifications_filename:
        :param monoisotopic: if true, masses are considered to be monoisotopic
        :raises:
          Exception: FileNotReadable if the modifications_filename could not be read
        :raises:
          Exception: FileNotFound if modifications_filename could not be found
        :raises:
          Exception: ParseError if modifications_filename could not be parsed
        """
        ...
    
    def getSpectra(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSpectra()
        Specifies a spectrum file to search
        """
        ...
    
    def setSpectra(self, spectra: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSpectra(const String & spectra)
        Specifies a spectrum file to search
        """
        ...
    
    def getDb(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDb()
        Specifies the name of a database (.trie file) to search
        """
        ...
    
    def setDb(self, db: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDb(const String & db)
        Specifies the name of a database (.trie file) to search
        """
        ...
    
    def getEnzyme(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEnzyme()
        Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        """
        ...
    
    def setEnzyme(self, enzyme: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setEnzyme(const String & enzyme)
        Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        """
        ...
    
    def getModificationsPerPeptide(self) -> int:
        """
        Cython signature: int getModificationsPerPeptide()
        Number of PTMs permitted in a single peptide
        """
        ...
    
    def setModificationsPerPeptide(self, modifications_per_peptide: int ) -> None:
        """
        Cython signature: void setModificationsPerPeptide(int modifications_per_peptide)
        Number of PTMs permitted in a single peptide
        """
        ...
    
    def getBlind(self) -> int:
        """
        Cython signature: unsigned int getBlind()
        Run inspect in a blind mode
        """
        ...
    
    def setBlind(self, blind: int ) -> None:
        """
        Cython signature: void setBlind(unsigned int blind)
        Run inspect in a blind mode
        """
        ...
    
    def getMaxPTMsize(self) -> float:
        """
        Cython signature: float getMaxPTMsize()
        The maximum modification size (in Da) to consider in a blind search
        """
        ...
    
    def setMaxPTMsize(self, maxptmsize: float ) -> None:
        """
        Cython signature: void setMaxPTMsize(float maxptmsize)
        The maximum modification size (in Da) to consider in a blind search
        """
        ...
    
    def getPrecursorMassTolerance(self) -> float:
        """
        Cython signature: float getPrecursorMassTolerance()
        Specifies the parent mass tolerance, in Daltons
        """
        ...
    
    def setPrecursorMassTolerance(self, precursor_mass_tolerance: float ) -> None:
        """
        Cython signature: void setPrecursorMassTolerance(float precursor_mass_tolerance)
        Specifies the parent mass tolerance, in Daltons
        """
        ...
    
    def getPeakMassTolerance(self) -> float:
        """
        Cython signature: float getPeakMassTolerance()
        How far b and y peaks can be shifted from their expected masses.
        """
        ...
    
    def setPeakMassTolerance(self, peak_mass_tolerance: float ) -> None:
        """
        Cython signature: void setPeakMassTolerance(float peak_mass_tolerance)
        How far b and y peaks can be shifted from their expected masses
        """
        ...
    
    def getMulticharge(self) -> int:
        """
        Cython signature: unsigned int getMulticharge()
        If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        """
        ...
    
    def setMulticharge(self, multicharge: int ) -> None:
        """
        Cython signature: void setMulticharge(unsigned int multicharge)
        If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        """
        ...
    
    def getInstrument(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getInstrument()
        If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        """
        ...
    
    def setInstrument(self, instrument: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setInstrument(const String & instrument)
        If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        """
        ...
    
    def getTagCount(self) -> int:
        """
        Cython signature: int getTagCount()
        Number of tags to generate
        """
        ...
    
    def setTagCount(self, TagCount: int ) -> None:
        """
        Cython signature: void setTagCount(int TagCount)
        Number of tags to generate
        """
        ...
    
    def __richcmp__(self, other: InspectInfile, op: int) -> Any:
        ... 


class IsotopePattern:
    """
    Cython implementation of _IsotopePattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1IsotopePattern.html>`_
    """
    
    spectrum: List[int]
    
    intensity: List[float]
    
    mz_score: List[float]
    
    theoretical_mz: List[float]
    
    theoretical_pattern: TheoreticalIsotopePattern
    
    @overload
    def __init__(self, size: int ) -> None:
        """
        Cython signature: void IsotopePattern(size_t size)
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopePattern ) -> None:
        """
        Cython signature: void IsotopePattern(IsotopePattern &)
        """
        ... 


class MapAlignmentEvaluationAlgorithmRecall:
    """
    Cython implementation of _MapAlignmentEvaluationAlgorithmRecall

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentEvaluationAlgorithmRecall.html>`_
      -- Inherits from ['MapAlignmentEvaluationAlgorithm']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MapAlignmentEvaluationAlgorithmRecall()
        """
        ... 


class MassTrace:
    """
    Cython implementation of _MassTrace

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1MassTrace.html>`_
    """
    
    max_rt: float
    
    theoretical_int: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassTrace()
        """
        ...
    
    @overload
    def __init__(self, in_0: MassTrace ) -> None:
        """
        Cython signature: void MassTrace(MassTrace &)
        """
        ...
    
    def getConvexhull(self) -> ConvexHull2D:
        """
        Cython signature: ConvexHull2D getConvexhull()
        """
        ...
    
    def updateMaximum(self) -> None:
        """
        Cython signature: void updateMaximum()
        """
        ...
    
    def getAvgMZ(self) -> float:
        """
        Cython signature: double getAvgMZ()
        """
        ...
    
    def isValid(self) -> bool:
        """
        Cython signature: bool isValid()
        """
        ... 


class MassTraces:
    """
    Cython implementation of _MassTraces

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1MassTraces.html>`_
    """
    
    max_trace: int
    
    baseline: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassTraces()
        """
        ...
    
    @overload
    def __init__(self, in_0: MassTraces ) -> None:
        """
        Cython signature: void MassTraces(MassTraces &)
        """
        ...
    
    def getPeakCount(self) -> int:
        """
        Cython signature: size_t getPeakCount()
        """
        ...
    
    def isValid(self, seed_mz: float , trace_tolerance: float ) -> bool:
        """
        Cython signature: bool isValid(double seed_mz, double trace_tolerance)
        """
        ...
    
    def getTheoreticalmaxPosition(self) -> int:
        """
        Cython signature: size_t getTheoreticalmaxPosition()
        """
        ...
    
    def updateBaseline(self) -> None:
        """
        Cython signature: void updateBaseline()
        """
        ...
    
    def getRTBounds(self) -> List[float, float]:
        """
        Cython signature: libcpp_pair[double,double] getRTBounds()
        """
        ... 


class MetaboTargetedTargetDecoy:
    """
    Cython implementation of _MetaboTargetedTargetDecoy

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedTargetDecoy.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboTargetedTargetDecoy()
        Resolve overlapping fragments and missing decoys for experimental specific decoy generation in targeted/pseudo targeted metabolomics
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboTargetedTargetDecoy ) -> None:
        """
        Cython signature: void MetaboTargetedTargetDecoy(MetaboTargetedTargetDecoy &)
        """
        ...
    
    def constructTargetDecoyMassMapping(self, t_exp: TargetedExperiment ) -> List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping]:
        """
        Cython signature: libcpp_vector[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] constructTargetDecoyMassMapping(TargetedExperiment & t_exp)
        Constructs a mass mapping of targets and decoys using the unique m_id identifier
        
        
        :param t_exp: TransitionExperiment holds compound and transition information used for the mapping
        """
        ...
    
    def resolveOverlappingTargetDecoyMassesByDecoyMassShift(self, t_exp: TargetedExperiment , mappings: List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] , mass_to_add: float , mz_tol: float , mz_tol_unit: String ) -> None:
        """
        Cython signature: void resolveOverlappingTargetDecoyMassesByDecoyMassShift(TargetedExperiment & t_exp, libcpp_vector[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] & mappings, double & mass_to_add, double & mz_tol, String & mz_tol_unit)
        Resolves overlapping target and decoy transition masses by adding a specifiable mass (e.g. CH2) to the overlapping decoy fragment
        
        
        :param t_exp: TransitionExperiment holds compound and transition information
        :param mappings: Map of identifier to target and decoy masses
        :param mass_to_add: (e.g. CH2)
        :param mz_tol: m/z tolerarance for target and decoy transition masses to be considered overlapping
        :param mz_tol_unit: m/z tolerance unit
        """
        ...
    
    def generateMissingDecoysByMassShift(self, t_exp: TargetedExperiment , mappings: List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] , mass_to_add: float ) -> None:
        """
        Cython signature: void generateMissingDecoysByMassShift(TargetedExperiment & t_exp, libcpp_vector[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] & mappings, double & mass_to_add)
        Generate a decoy for targets where fragmentation tree re-rooting was not possible, by adding a specifiable mass to the target fragments
        
        
        :param t_exp: TransitionExperiment holds compound and transition information
        :param mappings: Map of identifier to target and decoy masses
        :param mass_to_add: The maximum number of transitions required per assay
        """
        ... 


class MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping:
    """
    Cython implementation of _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ) -> None:
        """
        Cython signature: void MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping &)
        """
        ... 


class MzMLSqliteHandler:
    """
    Cython implementation of _MzMLSqliteHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1MzMLSqliteHandler.html>`_
    """
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , run_id: int ) -> None:
        """
        Cython signature: void MzMLSqliteHandler(String filename, uint64_t run_id)
        """
        ...
    
    @overload
    def __init__(self, in_0: MzMLSqliteHandler ) -> None:
        """
        Cython signature: void MzMLSqliteHandler(MzMLSqliteHandler &)
        """
        ...
    
    def readExperiment(self, exp: MSExperiment , meta_only: bool ) -> None:
        """
        Cython signature: void readExperiment(MSExperiment & exp, bool meta_only)
        Read an experiment into an MSExperiment structure
        
        
        :param exp: The result data structure
        :param meta_only: Only read the meta data
        """
        ...
    
    def readSpectra(self, exp: List[MSSpectrum] , indices: List[int] , meta_only: bool ) -> None:
        """
        Cython signature: void readSpectra(libcpp_vector[MSSpectrum] & exp, libcpp_vector[int] indices, bool meta_only)
        Read a set of spectra (potentially restricted to a subset)
        
        
        :param exp: The result data structure
        :param indices: A list of indices restricting the resulting spectra only to those specified here
        :param meta_only: Only read the meta data
        """
        ...
    
    def readChromatograms(self, exp: List[MSChromatogram] , indices: List[int] , meta_only: bool ) -> None:
        """
        Cython signature: void readChromatograms(libcpp_vector[MSChromatogram] & exp, libcpp_vector[int] indices, bool meta_only)
        Read a set of chromatograms (potentially restricted to a subset)
        
        
        :param exp: The result data structure
        :param indices: A list of indices restricting the resulting spectra only to those specified here
        :param meta_only: Only read the meta data
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        Returns number of spectra in the file, reutrns the number of spectra
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        Returns the number of chromatograms in the file
        """
        ...
    
    def setConfig(self, write_full_meta: bool , use_lossy_compression: bool , linear_abs_mass_acc: float ) -> None:
        """
        Cython signature: void setConfig(bool write_full_meta, bool use_lossy_compression, double linear_abs_mass_acc)
        Sets file configuration
        
        
        :param write_full_meta: Whether to write a complete mzML meta data structure into the RUN_EXTRA field (allows complete recovery of the input file)
        :param use_lossy_compression: Whether to use lossy compression (ms numpress)
        :param linear_abs_mass_acc: Accepted loss in mass accuracy (absolute m/z, in Th)
        """
        ...
    
    def getSpectraIndicesbyRT(self, RT: float , deltaRT: float , indices: List[int] ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] getSpectraIndicesbyRT(double RT, double deltaRT, libcpp_vector[int] indices)
        Returns spectral indices around a specific retention time
        
        :param RT: The retention time
        :param deltaRT: Tolerance window around RT (if less or equal than zero, only the first spectrum *after* RT is returned)
        :param indices: Spectra to consider (if empty, all spectra are considered)
        :return: The indices of the spectra within RT +/- deltaRT
        """
        ...
    
    def writeExperiment(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void writeExperiment(MSExperiment exp)
        Write an MSExperiment to disk
        """
        ...
    
    def createTables(self) -> None:
        """
        Cython signature: void createTables()
        Create data tables for a new file
        """
        ...
    
    def writeSpectra(self, spectra: List[MSSpectrum] ) -> None:
        """
        Cython signature: void writeSpectra(libcpp_vector[MSSpectrum] spectra)
        Writes a set of spectra to disk
        """
        ...
    
    def writeChromatograms(self, chroms: List[MSChromatogram] ) -> None:
        """
        Cython signature: void writeChromatograms(libcpp_vector[MSChromatogram] chroms)
        Writes a set of chromatograms to disk
        """
        ...
    
    def writeRunLevelInformation(self, exp: MSExperiment , write_full_meta: bool ) -> None:
        """
        Cython signature: void writeRunLevelInformation(MSExperiment exp, bool write_full_meta)
        Write the run-level information for an experiment into tables
        
        This is a low level function, do not call this function unless you know what you are doing
        
        
        :param exp: The result data structure
        :param meta_only: Only read the meta data
        """
        ...
    
    def getRunID(self) -> int:
        """
        Cython signature: uint64_t getRunID()
        Extract the `RUN` ID from the sqMass file
        """
        ... 


class OpenPepXLAlgorithm:
    """
    Cython implementation of _OpenPepXLAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenPepXLAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenPepXLAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenPepXLAlgorithm ) -> None:
        """
        Cython signature: void OpenPepXLAlgorithm(OpenPepXLAlgorithm &)
        """
        ...
    
    def run(self, unprocessed_spectra: MSExperiment , cfeatures: ConsensusMap , fasta_db: List[FASTAEntry] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , preprocessed_pair_spectra: OPXL_PreprocessedPairSpectra , spectrum_pairs: List[List[int, int]] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment ) -> int:
        """
        Cython signature: OpenPepXLAlgorithm_ExitCodes run(MSExperiment & unprocessed_spectra, ConsensusMap & cfeatures, libcpp_vector[FASTAEntry] & fasta_db, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, OPXL_PreprocessedPairSpectra & preprocessed_pair_spectra, libcpp_vector[libcpp_pair[size_t,size_t]] & spectrum_pairs, libcpp_vector[libcpp_vector[CrossLinkSpectrumMatch]] & all_top_csms, MSExperiment & spectra)
        Performs the main function of this class, the search for cross-linked peptides
        
        
        :param unprocessed_spectra: The input PeakMap of experimental spectra
        :param cfeatures: The input cfeatures
        :param fasta_db: The protein database containing targets and decoys
        :param protein_ids: A result vector containing search settings. Should contain one PeptideIdentification
        :param peptide_ids: A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty
        :param preprocessed_pair_spectra: A result structure containing linear and cross-linked ion spectra. Will be overwritten. This is only necessary for writing out xQuest type spectrum files
        :param spectrum_pairs: A result vector containing paired spectra indices. Should be empty. This is only necessary for writing out xQuest type spectrum files
        :param all_top_csms: A result vector containing cross-link spectrum matches as CrossLinkSpectrumMatches. Should be empty. This is only necessary for writing out xQuest type spectrum files
        :param spectra: A result vector containing the input spectra after preprocessing and filtering. Should be empty. This is only necessary for writing out xQuest type spectrum files
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
    OpenPepXLAlgorithm_ExitCodes : __OpenPepXLAlgorithm_ExitCodes 


class PeakWidthEstimator:
    """
    Cython implementation of _PeakWidthEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakWidthEstimator.html>`_

    Rough estimation of the peak width at m/z
    
    Based on the peaks of the dataset (peak position & width) and the peak
    boundaries as reported by the PeakPickerHiRes, the typical peak width is
    estimated for arbitrary m/z using a spline interpolationThis struct can be used to store both peak or feature indices`
    """
    
    @overload
    def __init__(self, in_0: PeakWidthEstimator ) -> None:
        """
        Cython signature: void PeakWidthEstimator(PeakWidthEstimator &)
        """
        ...
    
    @overload
    def __init__(self, exp_picked: MSExperiment , boundaries: List[List[PeakBoundary]] ) -> None:
        """
        Cython signature: void PeakWidthEstimator(MSExperiment exp_picked, libcpp_vector[libcpp_vector[PeakBoundary]] & boundaries)
        """
        ...
    
    def getPeakWidth(self, mz: float ) -> float:
        """
        Cython signature: double getPeakWidth(double mz)
        Returns the estimated peak width at m/z
        """
        ... 


class Product:
    """
    Cython implementation of _Product

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Product.html>`_

    This class describes the product isolation window for special scan types, such as MRM
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Product()
        """
        ...
    
    @overload
    def __init__(self, in_0: Product ) -> None:
        """
        Cython signature: void Product(Product &)
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        Returns the target m/z
        """
        ...
    
    def setMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setMZ(double)
        Sets the target m/z
        """
        ...
    
    def getIsolationWindowLowerOffset(self) -> float:
        """
        Cython signature: double getIsolationWindowLowerOffset()
        Returns the lower offset from the target m/z
        """
        ...
    
    def setIsolationWindowLowerOffset(self, bound: float ) -> None:
        """
        Cython signature: void setIsolationWindowLowerOffset(double bound)
        Sets the lower offset from the target m/z
        """
        ...
    
    def getIsolationWindowUpperOffset(self) -> float:
        """
        Cython signature: double getIsolationWindowUpperOffset()
        Returns the upper offset from the target m/z
        """
        ...
    
    def setIsolationWindowUpperOffset(self, bound: float ) -> None:
        """
        Cython signature: void setIsolationWindowUpperOffset(double bound)
        Sets the upper offset from the target m/z
        """
        ...
    
    def __richcmp__(self, other: Product, op: int) -> Any:
        ... 


class RankScaler:
    """
    Cython implementation of _RankScaler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RankScaler.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void RankScaler()
        """
        ...
    
    def filterSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
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


class Seed:
    """
    Cython implementation of _Seed

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1Seed.html>`_
    """
    
    spectrum: int
    
    peak: int
    
    intensity: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Seed()
        """
        ...
    
    @overload
    def __init__(self, in_0: Seed ) -> None:
        """
        Cython signature: void Seed(Seed &)
        """
        ...
    
    def __richcmp__(self, other: Seed, op: int) -> Any:
        ... 


class SourceFile:
    """
    Cython implementation of _SourceFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SourceFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SourceFile()
        Description of a file location, used to store the origin of (meta) data
        """
        ...
    
    @overload
    def __init__(self, in_0: SourceFile ) -> None:
        """
        Cython signature: void SourceFile(SourceFile &)
        """
        ...
    
    def getNameOfFile(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNameOfFile()
        Returns the file name
        """
        ...
    
    def setNameOfFile(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNameOfFile(String)
        Sets the file name
        """
        ...
    
    def getPathToFile(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPathToFile()
        Returns the file path
        """
        ...
    
    def setPathToFile(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPathToFile(String)
        Sets the file path
        """
        ...
    
    def getFileSize(self) -> float:
        """
        Cython signature: float getFileSize()
        Returns the file size in MB
        """
        ...
    
    def setFileSize(self, in_0: float ) -> None:
        """
        Cython signature: void setFileSize(float)
        Sets the file size in MB
        """
        ...
    
    def getFileType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFileType()
        Returns the file type
        """
        ...
    
    def setFileType(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFileType(String)
        Sets the file type
        """
        ...
    
    def getChecksum(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getChecksum()
        Returns the file's checksum
        """
        ...
    
    def setChecksum(self, in_0: Union[bytes, str, String] , in_1: int ) -> None:
        """
        Cython signature: void setChecksum(String, ChecksumType)
        Sets the file's checksum
        """
        ...
    
    def getChecksumType(self) -> int:
        """
        Cython signature: ChecksumType getChecksumType()
        Returns the checksum type
        """
        ...
    
    def getNativeIDType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeIDType()
        Returns the native ID type of the spectra
        """
        ...
    
    def setNativeIDType(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeIDType(String)
        Sets the native ID type of the spectra
        """
        ...
    
    def getNativeIDTypeAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeIDTypeAccession()
        Returns the nativeID of the spectra
        """
        ...
    
    def setNativeIDTypeAccession(self, accesssion: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeIDTypeAccession(const String & accesssion)
        Sets the native ID of the spectra
        """
        ... 


class Tagger:
    """
    Cython implementation of _Tagger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Tagger.html>`_

    Constructor for Tagger
    
    The parameter `max_charge_` should be >= `min_charge_`
    Also `max_tag_length` should be >= `min_tag_length`
    
    :param min_tag_length: The minimal sequence tag length
    :param ppm: The tolerance for matching residue masses to peak delta masses
    :param max_tag_length: The maximal sequence tag length
    :param min_charge: Minimal fragment charge considered for each sequence tag
    :param max_charge: Maximal fragment charge considered for each sequence tag
    :param fixed_mods: A list of modification names. The modified residues replace the unmodified versions
    :param var_mods: A list of modification names. The modified residues are added as additional entries to the list of residues
    """
    
    @overload
    def __init__(self, in_0: Tagger ) -> None:
        """
        Cython signature: void Tagger(Tagger &)
        """
        ...
    
    @overload
    def __init__(self, min_tag_length: int , ppm: float , max_tag_length: int , min_charge: int , max_charge: int , fixed_mods: List[bytes] , var_mods: List[bytes] ) -> None:
        """
        Cython signature: void Tagger(size_t min_tag_length, double ppm, size_t max_tag_length, size_t min_charge, size_t max_charge, const StringList & fixed_mods, const StringList & var_mods)
        """
        ...
    
    @overload
    def getTag(self, mzs: List[float] , tags: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void getTag(const libcpp_vector[double] & mzs, libcpp_vector[libcpp_utf8_string] & tags)
        Generate tags from mass vector `mzs`
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param mzs: A vector of mz values, containing the mz values from a centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
        """
        ...
    
    @overload
    def getTag(self, spec: MSSpectrum , tags: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void getTag(const MSSpectrum & spec, libcpp_vector[libcpp_utf8_string] & tags)
        Generate tags from an MSSpectrum
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param spec: A centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
        """
        ...
    
    def setMaxCharge(self, max_charge: int ) -> None:
        """
        Cython signature: void setMaxCharge(size_t max_charge)
        Change the maximal charge considered by the tagger
        
        Allows to change the maximal considered charge e.g. based on a spectra
        precursor charge without calling the constructor multiple times
        
        :param max_charge: The new maximal charge
        """
        ... 


class TheoreticalIsotopePattern:
    """
    Cython implementation of _TheoreticalIsotopePattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1TheoreticalIsotopePattern.html>`_
    """
    
    intensity: List[float]
    
    optional_begin: int
    
    optional_end: int
    
    max: float
    
    trimmed_left: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TheoreticalIsotopePattern()
        """
        ...
    
    @overload
    def __init__(self, in_0: TheoreticalIsotopePattern ) -> None:
        """
        Cython signature: void TheoreticalIsotopePattern(TheoreticalIsotopePattern &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ... 


class __AggregationMethod:
    """
          Aggregation method
    """
    PROD : int
    SUM : int
    BEST : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class ChecksumType:
    None
    UNKNOWN_CHECKSUM : int
    SHA1 : int
    MD5 : int
    SIZE_OF_CHECKSUMTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __OpenPepXLAlgorithm_ExitCodes:
    None
    EXECUTION_OK : int
    ILLEGAL_PARAMETERS : int
    UNEXPECTED_RESULT : int
    INCOMPATIBLE_INPUT_DATA : int

    def getMapping(self) -> Dict[int, str]:
       ... 

