from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_CachedmzML_load(filename: Union[bytes, str, String] , exp: CachedmzML ) -> None:
    """
    Cython signature: void load(const String & filename, CachedmzML & exp)
    """
    ...

def __static_CachedmzML_store(filename: Union[bytes, str, String] , exp: MSExperiment ) -> None:
    """
    Cython signature: void store(const String & filename, MSExperiment exp)
    """
    ...


class AnnotationStatistics:
    """
    Cython implementation of _AnnotationStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AnnotationStatistics.html>`_
    """
    
    states: List[int]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AnnotationStatistics()
        """
        ...
    
    @overload
    def __init__(self, in_0: AnnotationStatistics ) -> None:
        """
        Cython signature: void AnnotationStatistics(AnnotationStatistics &)
        """
        ...
    
    def __richcmp__(self, other: AnnotationStatistics, op: int) -> Any:
        ... 


class CachedmzML:
    """
    Cython implementation of _CachedmzML

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CachedmzML.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CachedmzML()
        A class that uses on-disk caching to read and write spectra and chromatograms
        """
        ...
    
    @overload
    def __init__(self, in_0: CachedmzML ) -> None:
        """
        Cython signature: void CachedmzML(CachedmzML &)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void CachedmzML(String filename)
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
    
    def getSpectrum(self, idx: int ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum getSpectrum(size_t idx)
        """
        ...
    
    def getChromatogram(self, idx: int ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getChromatogram(size_t idx)
        """
        ...
    
    def getMetaData(self) -> MSExperiment:
        """
        Cython signature: MSExperiment getMetaData()
        """
        ...
    
    load: __static_CachedmzML_load
    
    store: __static_CachedmzML_store 


class DecoyGenerator:
    """
    Cython implementation of _DecoyGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DecoyGenerator.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DecoyGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: DecoyGenerator ) -> None:
        """
        Cython signature: void DecoyGenerator(DecoyGenerator &)
        """
        ...
    
    def setSeed(self, in_0: int ) -> None:
        """
        Cython signature: void setSeed(uint64_t)
        """
        ...
    
    def reverseProtein(self, protein: AASequence ) -> AASequence:
        """
        Cython signature: AASequence reverseProtein(const AASequence & protein)
        Reverses the protein sequence
        """
        ...
    
    def reversePeptides(self, protein: AASequence , protease: Union[bytes, str, String] ) -> AASequence:
        """
        Cython signature: AASequence reversePeptides(const AASequence & protein, const String & protease)
        Reverses the protein's peptide sequences between enzymatic cutting positions
        """
        ...
    
    def shufflePeptides(self, aas: AASequence , protease: Union[bytes, str, String] , max_attempts: int ) -> AASequence:
        """
        Cython signature: AASequence shufflePeptides(const AASequence & aas, const String & protease, const int max_attempts)
        Shuffle the protein's peptide sequences between enzymatic cutting positions, each peptide is shuffled @param max_attempts times to minimize sequence identity
        """
        ... 


class EmpiricalFormula:
    """
    Cython implementation of _EmpiricalFormula

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmpiricalFormula.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EmpiricalFormula()
        Representation of an empirical formula
        """
        ...
    
    @overload
    def __init__(self, in_0: EmpiricalFormula ) -> None:
        """
        Cython signature: void EmpiricalFormula(EmpiricalFormula &)
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void EmpiricalFormula(String)
        EmpiricalFormula Constructor from string
        """
        ...
    
    @overload
    def __init__(self, number: int , element: Element , charge: int ) -> None:
        """
        Cython signature: void EmpiricalFormula(ptrdiff_t number, Element * element, ptrdiff_t charge)
        EmpiricalFormula Constructor with element pointer and number
        """
        ...
    
    def getMonoWeight(self) -> float:
        """
        Cython signature: double getMonoWeight()
        Returns the mono isotopic weight of the formula (includes proton charges)
        """
        ...
    
    def getAverageWeight(self) -> float:
        """
        Cython signature: double getAverageWeight()
        Returns the average weight of the formula (includes proton charges)
        """
        ...
    
    def estimateFromWeightAndComp(self, average_weight: float , C: float , H: float , N: float , O: float , S: float , P: float ) -> bool:
        """
        Cython signature: bool estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P)
        Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry
        """
        ...
    
    def estimateFromWeightAndCompAndS(self, average_weight: float , S: int , C: float , H: float , N: float , O: float , P: float ) -> bool:
        """
        Cython signature: bool estimateFromWeightAndCompAndS(double average_weight, unsigned int S, double C, double H, double N, double O, double P)
        Fills this EmpiricalFormula with an approximate elemental composition for a given average weight, exact number of sulfurs, and approximate elemental stoichiometry
        """
        ...
    
    @overload
    def getIsotopeDistribution(self, in_0: CoarseIsotopePatternGenerator ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution getIsotopeDistribution(CoarseIsotopePatternGenerator)
        Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator or the FineIsotopePatternGenerator method
        """
        ...
    
    @overload
    def getIsotopeDistribution(self, in_0: FineIsotopePatternGenerator ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution getIsotopeDistribution(FineIsotopePatternGenerator)
        """
        ...
    
    def getConditionalFragmentIsotopeDist(self, precursor: EmpiricalFormula , precursor_isotopes: Set[int] , method: CoarseIsotopePatternGenerator ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution getConditionalFragmentIsotopeDist(EmpiricalFormula & precursor, libcpp_set[unsigned int] & precursor_isotopes, CoarseIsotopePatternGenerator method)
        """
        ...
    
    def getNumberOfAtoms(self) -> int:
        """
        Cython signature: size_t getNumberOfAtoms()
        Returns the total number of atoms
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: ptrdiff_t getCharge()
        Returns the total charge
        """
        ...
    
    def setCharge(self, charge: int ) -> None:
        """
        Cython signature: void setCharge(ptrdiff_t charge)
        Sets the charge
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the formula as a string (charges are not included)
        """
        ...
    
    def getElementalComposition(self) -> Dict[bytes, int]:
        """
        Cython signature: libcpp_map[libcpp_string,int] getElementalComposition()
        Get elemental composition as a hash {'Symbol' -> NrAtoms}
        """
        ...
    
    def isEmpty(self) -> bool:
        """
        Cython signature: bool isEmpty()
        Returns true if the formula does not contain a element
        """
        ...
    
    def isCharged(self) -> bool:
        """
        Cython signature: bool isCharged()
        Returns true if charge is not equal to zero
        """
        ...
    
    def hasElement(self, element: Element ) -> bool:
        """
        Cython signature: bool hasElement(Element * element)
        Returns true if the formula contains the element
        """
        ...
    
    def contains(self, ef: EmpiricalFormula ) -> bool:
        """
        Cython signature: bool contains(EmpiricalFormula ef)
        Returns true if all elements from `ef` ( empirical formula ) are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula
        """
        ...
    
    def __add__(self: EmpiricalFormula, other: EmpiricalFormula) -> EmpiricalFormula:
        ...
    
    def __sub__(self: EmpiricalFormula, other: EmpiricalFormula) -> EmpiricalFormula:
        ...
    
    def __iadd__(self: EmpiricalFormula, other: EmpiricalFormula) -> EmpiricalFormula:
        ...
    
    def __isub__(self: EmpiricalFormula, other: EmpiricalFormula) -> EmpiricalFormula:
        ...
    
    def calculateTheoreticalIsotopesNumber(self) -> float:
        """
        Cython signature: double calculateTheoreticalIsotopesNumber()
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        Returns the formula as a string (charges are not included)
        """
        ...
    
    def __richcmp__(self, other: EmpiricalFormula, op: int) -> Any:
        ... 


class FeatureDistance:
    """
    Cython implementation of _FeatureDistance

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureDistance.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, max_intensity: float , force_constraints: bool ) -> None:
        """
        Cython signature: void FeatureDistance(double max_intensity, bool force_constraints)
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureDistance ) -> None:
        """
        Cython signature: void FeatureDistance(FeatureDistance &)
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


class FeatureGroupingAlgorithmKD:
    """
    Cython implementation of _FeatureGroupingAlgorithmKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmKD.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm', 'ProgressLogger']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureGroupingAlgorithmKD()
        A feature grouping algorithm for unlabeled data
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


class IsotopeLabelingMDVs:
    """
    Cython implementation of _IsotopeLabelingMDVs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeLabelingMDVs.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsotopeLabelingMDVs()
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopeLabelingMDVs ) -> None:
        """
        Cython signature: void IsotopeLabelingMDVs(IsotopeLabelingMDVs &)
        """
        ...
    
    def isotopicCorrection(self, normalized_feature: Feature , corrected_feature: Feature , correction_matrix: MatrixDouble , correction_matrix_agent: int ) -> None:
        """
        Cython signature: void isotopicCorrection(const Feature & normalized_feature, Feature & corrected_feature, MatrixDouble & correction_matrix, const DerivatizationAgent & correction_matrix_agent)
        This function performs an isotopic correction to account for unlabeled abundances coming from
        the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
        
        
        :param normalized_feature: Feature with normalized values for each component and unlabeled chemical formula for each component group
        :param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
        :param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
        :return: corrected_feature: Feature with corrected values for each component
        """
        ...
    
    def isotopicCorrections(self, normalized_featureMap: FeatureMap , corrected_featureMap: FeatureMap , correction_matrix: MatrixDouble , correction_matrix_agent: int ) -> None:
        """
        Cython signature: void isotopicCorrections(const FeatureMap & normalized_featureMap, FeatureMap & corrected_featureMap, MatrixDouble & correction_matrix, const DerivatizationAgent & correction_matrix_agent)
        This function performs an isotopic correction to account for unlabeled abundances coming from
        the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
        
        
        :param normalized_featuremap: FeatureMap with normalized values for each component and unlabeled chemical formula for each component group
        :param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
        :param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
        :return corrected_featuremap: FeatureMap with corrected values for each component
        """
        ...
    
    def calculateIsotopicPurity(self, normalized_feature: Feature , experiment_data: List[float] , isotopic_purity_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void calculateIsotopicPurity(const Feature & normalized_feature, const libcpp_vector[double] & experiment_data, const String & isotopic_purity_name)
        This function calculates the isotopic purity of the MDV using the following formula:
        isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
        where n in M+n is represented as the index of the result
        The formula is extracted from "High-resolution 13C metabolic flux analysis",
        Long et al, doi:10.1038/s41596-019-0204-0
        
        
        :param normalized_feature: Feature with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a Feature with the calculated isotopic purity for the component group
        :param experiment_data: Vector of experiment data in percent
        :param isotopic_purity_name: Name of the isotopic purity tracer to be saved as a meta value
        """
        ...
    
    def calculateMDVAccuracy(self, normalized_feature: Feature , feature_name: Union[bytes, str, String] , fragment_isotopomer_theoretical_formula: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void calculateMDVAccuracy(const Feature & normalized_feature, const String & feature_name, const String & fragment_isotopomer_theoretical_formula)
        This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
        using average deviation to the mean. The result is mapped to the meta value "average_accuracy" in the updated feature
        
        
        :param normalized_feature: Feature with normalized values for each component and the chemical formula of the component group. Out is a Feature with the component group accuracy and accuracy for the error for each component
        :param fragment_isotopomer_measured: Measured scan values
        :param fragment_isotopomer_theoretical_formula: Empirical formula from which the theoretical values will be generated
        """
        ...
    
    def calculateMDVAccuracies(self, normalized_featureMap: FeatureMap , feature_name: Union[bytes, str, String] , fragment_isotopomer_theoretical_formulas: Dict[Union[bytes, str], Union[bytes, str]] ) -> None:
        """
        Cython signature: void calculateMDVAccuracies(const FeatureMap & normalized_featureMap, const String & feature_name, const libcpp_map[libcpp_utf8_string,libcpp_utf8_string] & fragment_isotopomer_theoretical_formulas)
        This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
        using average deviation to the mean
        
        
        param normalized_featuremap: FeatureMap with normalized values for each component and the chemical formula of the component group. Out is a FeatureMap with the component group accuracy and accuracy for the error for each component
        param fragment_isotopomer_measured: Measured scan values
        param fragment_isotopomer_theoretical_formula: A map of ProteinName/peptideRef to Empirical formula from which the theoretical values will be generated
        """
        ...
    
    def calculateMDV(self, measured_feature: Feature , normalized_feature: Feature , mass_intensity_type: int , feature_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void calculateMDV(const Feature & measured_feature, Feature & normalized_feature, const MassIntensityType & mass_intensity_type, const String & feature_name)
        """
        ...
    
    def calculateMDVs(self, measured_featureMap: FeatureMap , normalized_featureMap: FeatureMap , mass_intensity_type: int , feature_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void calculateMDVs(const FeatureMap & measured_featureMap, FeatureMap & normalized_featureMap, const MassIntensityType & mass_intensity_type, const String & feature_name)
        """
        ... 


class Kernel_MassTrace:
    """
    Cython implementation of _Kernel_MassTrace

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Kernel_MassTrace.html>`_
    """
    
    fwhm_mz_avg: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Kernel_MassTrace()
        """
        ...
    
    @overload
    def __init__(self, in_0: Kernel_MassTrace ) -> None:
        """
        Cython signature: void Kernel_MassTrace(Kernel_MassTrace &)
        """
        ...
    
    @overload
    def __init__(self, trace_peaks: List[Peak2D] ) -> None:
        """
        Cython signature: void Kernel_MassTrace(const libcpp_vector[Peak2D] & trace_peaks)
        """
        ...
    
    def getSize(self) -> int:
        """
        Cython signature: size_t getSize()
        Returns the number of peaks contained in the mass trace
        """
        ...
    
    def getLabel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getLabel()
        Returns label of mass trace
        """
        ...
    
    def setLabel(self, label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLabel(String label)
        Sets label of mass trace
        """
        ...
    
    def getCentroidMZ(self) -> float:
        """
        Cython signature: double getCentroidMZ()
        Returns the centroid m/z
        """
        ...
    
    def getCentroidRT(self) -> float:
        """
        Cython signature: double getCentroidRT()
        Returns the centroid RT
        """
        ...
    
    def getCentroidSD(self) -> float:
        """
        Cython signature: double getCentroidSD()
        Returns the centroid SD
        """
        ...
    
    def getFWHM(self) -> float:
        """
        Cython signature: double getFWHM()
        Returns FWHM
        """
        ...
    
    def getTraceLength(self) -> float:
        """
        Cython signature: double getTraceLength()
        Returns the length of the trace (as difference in RT)
        """
        ...
    
    def getFWHMborders(self) -> List[int, int]:
        """
        Cython signature: libcpp_pair[size_t,size_t] getFWHMborders()
        Returns FWHM boarders
        """
        ...
    
    def getSmoothedIntensities(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getSmoothedIntensities()
        Returns smoothed intensities (empty if no smoothing was explicitly done beforehand!)
        """
        ...
    
    def getAverageMS1CycleTime(self) -> float:
        """
        Cython signature: double getAverageMS1CycleTime()
        Returns average scan time of mass trace
        """
        ...
    
    def computeSmoothedPeakArea(self) -> float:
        """
        Cython signature: double computeSmoothedPeakArea()
        Sums all non-negative (smoothed!) intensities in the mass trace
        """
        ...
    
    def computePeakArea(self) -> float:
        """
        Cython signature: double computePeakArea()
        Sums intensities of all peaks in the mass trace
        """
        ...
    
    def computeIntensitySum(self) -> float:
        """
        Cython signature: double computeIntensitySum()
        Sum all peak intensities in the mass trace
        """
        ...
    
    def findMaxByIntPeak(self, in_0: bool ) -> int:
        """
        Cython signature: size_t findMaxByIntPeak(bool)
        Returns the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities)
        """
        ...
    
    def estimateFWHM(self, in_0: bool ) -> int:
        """
        Cython signature: size_t estimateFWHM(bool)
        Estimates FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities)
        """
        ...
    
    def computeFwhmArea(self) -> float:
        """
        Cython signature: double computeFwhmArea()
        """
        ...
    
    def computeFwhmAreaSmooth(self) -> float:
        """
        Cython signature: double computeFwhmAreaSmooth()
        Computes chromatographic peak area within the FWHM range.
        """
        ...
    
    def getIntensity(self, in_0: bool ) -> float:
        """
        Cython signature: double getIntensity(bool)
        Returns the intensity
        """
        ...
    
    def getMaxIntensity(self, in_0: bool ) -> float:
        """
        Cython signature: double getMaxIntensity(bool)
        Returns the max intensity
        """
        ...
    
    def getConvexhull(self) -> ConvexHull2D:
        """
        Cython signature: ConvexHull2D getConvexhull()
        Returns the mass trace's convex hull
        """
        ...
    
    def setCentroidSD(self, tmp_sd: float ) -> None:
        """
        Cython signature: void setCentroidSD(double & tmp_sd)
        """
        ...
    
    def setSmoothedIntensities(self, db_vec: List[float] ) -> None:
        """
        Cython signature: void setSmoothedIntensities(libcpp_vector[double] & db_vec)
        Sets smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing)
        """
        ...
    
    def updateSmoothedMaxRT(self) -> None:
        """
        Cython signature: void updateSmoothedMaxRT()
        """
        ...
    
    def updateWeightedMeanRT(self) -> None:
        """
        Cython signature: void updateWeightedMeanRT()
        Compute & update centroid RT as a intensity-weighted mean of RTs
        """
        ...
    
    def updateSmoothedWeightedMeanRT(self) -> None:
        """
        Cython signature: void updateSmoothedWeightedMeanRT()
        """
        ...
    
    def updateMedianRT(self) -> None:
        """
        Cython signature: void updateMedianRT()
        Compute & update centroid RT as median position of intensities
        """
        ...
    
    def updateMedianMZ(self) -> None:
        """
        Cython signature: void updateMedianMZ()
        Compute & update centroid m/z as median of m/z values
        """
        ...
    
    def updateMeanMZ(self) -> None:
        """
        Cython signature: void updateMeanMZ()
        Compute & update centroid m/z as mean of m/z values
        """
        ...
    
    def updateWeightedMeanMZ(self) -> None:
        """
        Cython signature: void updateWeightedMeanMZ()
        Compute & update centroid m/z as weighted mean of m/z values
        """
        ...
    
    def updateWeightedMZsd(self) -> None:
        """
        Cython signature: void updateWeightedMZsd()
        Compute & update m/z standard deviation of mass trace as weighted mean of m/z values
        
        Make sure to call update(Weighted)(Mean|Median)MZ() first! <br>
        use getCentroidSD() to get result
        """
        ...
    
    def setQuantMethod(self, method: int ) -> None:
        """
        Cython signature: void setQuantMethod(MT_QUANTMETHOD method)
        Determine if area or median is used for quantification
        """
        ...
    
    def getQuantMethod(self) -> int:
        """
        Cython signature: MT_QUANTMETHOD getQuantMethod()
        Check if area or median is used for quantification
        """
        ... 


class MSChromatogram:
    """
    Cython implementation of _MSChromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSChromatogram.html>`_
      -- Inherits from ['ChromatogramSettings', 'RangeManagerRtInt']

    The representation of a chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    Iterations yields access to underlying peak objects but is slower
    Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
    See help(ChromatogramSettings) for information about meta-information
    
    Usage:
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSChromatogram()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSChromatogram ) -> None:
        """
        Cython signature: void MSChromatogram(MSChromatogram &)
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        Returns the mz of the product entry, makes sense especially for MRM scans
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
        Cython signature: void setName(String)
        Sets the name
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
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
    
    def __getitem__(self, in_0: int ) -> ChromatogramPeak:
        """
        Cython signature: ChromatogramPeak & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: ChromatogramPeak ) -> None:
        """Cython signature: ChromatogramPeak & operator[](size_t)"""
        ...
    
    def updateRanges(self) -> None:
        """
        Cython signature: void updateRanges()
        """
        ...
    
    def clear(self, in_0: int ) -> None:
        """
        Cython signature: void clear(int)
        Clears all data and meta data
        
        
        :param clear_meta_data: If true, all meta data is cleared in addition to the data
        """
        ...
    
    def push_back(self, in_0: ChromatogramPeak ) -> None:
        """
        Cython signature: void push_back(ChromatogramPeak)
        Append a peak
        """
        ...
    
    def isSorted(self) -> bool:
        """
        Cython signature: bool isSorted()
        Checks if all peaks are sorted with respect to ascending RT
        """
        ...
    
    def sortByIntensity(self, reverse: bool ) -> None:
        """
        Cython signature: void sortByIntensity(bool reverse)
        Lexicographically sorts the peaks by their intensity
        
        
        Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly
        """
        ...
    
    def sortByPosition(self) -> None:
        """
        Cython signature: void sortByPosition()
        Lexicographically sorts the peaks by their position
        
        
        The chromatogram is sorted with respect to position. Meta data arrays will be sorted accordingly
        """
        ...
    
    def findNearest(self, in_0: float ) -> int:
        """
        Cython signature: int findNearest(double)
        Binary search for the peak nearest to a specific RT
        :note: Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined
        
        
        :param rt: The searched for mass-to-charge ratio searched
        :return: Returns the index of the peak.
        :raises:
          Exception: Precondition is thrown if the chromatogram is empty (not only in debug mode)
        """
        ...
    
    def getFloatDataArrays(self) -> List[FloatDataArray]:
        """
        Cython signature: libcpp_vector[FloatDataArray] getFloatDataArrays()
        Returns a reference to the float meta data arrays
        """
        ...
    
    def getIntegerDataArrays(self) -> List[IntegerDataArray]:
        """
        Cython signature: libcpp_vector[IntegerDataArray] getIntegerDataArrays()
        Returns a reference to the integer meta data arrays
        """
        ...
    
    def getStringDataArrays(self) -> List[StringDataArray]:
        """
        Cython signature: libcpp_vector[StringDataArray] getStringDataArrays()
        Returns a reference to the string meta data arrays
        """
        ...
    
    def setFloatDataArrays(self, fda: List[FloatDataArray] ) -> None:
        """
        Cython signature: void setFloatDataArrays(libcpp_vector[FloatDataArray] fda)
        Sets the float meta data arrays
        """
        ...
    
    def setIntegerDataArrays(self, ida: List[IntegerDataArray] ) -> None:
        """
        Cython signature: void setIntegerDataArrays(libcpp_vector[IntegerDataArray] ida)
        Sets the integer meta data arrays
        """
        ...
    
    def setStringDataArrays(self, sda: List[StringDataArray] ) -> None:
        """
        Cython signature: void setStringDataArrays(libcpp_vector[StringDataArray] sda)
        Sets the string meta data arrays
        """
        ...
    
    def getProduct(self) -> Product:
        """
        Cython signature: Product getProduct()
        Returns the product ion
        """
        ...
    
    def setProduct(self, p: Product ) -> None:
        """
        Cython signature: void setProduct(Product p)
        Sets the product ion
        """
        ...
    
    def getNativeID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeID()
        Returns the native identifier for the spectrum, used by the acquisition software.
        """
        ...
    
    def setNativeID(self, native_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeID(String native_id)
        Sets the native identifier for the spectrum, used by the acquisition software.
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the free-text comment
        """
        ...
    
    def setComment(self, comment: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String comment)
        Sets the free-text comment
        """
        ...
    
    def getInstrumentSettings(self) -> InstrumentSettings:
        """
        Cython signature: InstrumentSettings getInstrumentSettings()
        Returns the instrument settings of the current spectrum
        """
        ...
    
    def setInstrumentSettings(self, instrument_settings: InstrumentSettings ) -> None:
        """
        Cython signature: void setInstrumentSettings(InstrumentSettings instrument_settings)
        Sets the instrument settings of the current spectrum
        """
        ...
    
    def getAcquisitionInfo(self) -> AcquisitionInfo:
        """
        Cython signature: AcquisitionInfo getAcquisitionInfo()
        Returns the acquisition info
        """
        ...
    
    def setAcquisitionInfo(self, acquisition_info: AcquisitionInfo ) -> None:
        """
        Cython signature: void setAcquisitionInfo(AcquisitionInfo acquisition_info)
        Sets the acquisition info
        """
        ...
    
    def getSourceFile(self) -> SourceFile:
        """
        Cython signature: SourceFile getSourceFile()
        Returns the source file
        """
        ...
    
    def setSourceFile(self, source_file: SourceFile ) -> None:
        """
        Cython signature: void setSourceFile(SourceFile source_file)
        Sets the source file
        """
        ...
    
    def getPrecursor(self) -> Precursor:
        """
        Cython signature: Precursor getPrecursor()
        Returns the precursors
        """
        ...
    
    def setPrecursor(self, precursor: Precursor ) -> None:
        """
        Cython signature: void setPrecursor(Precursor precursor)
        Sets the precursors
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[shared_ptr[DataProcessing]] getDataProcessing()
        Returns the description of the applied processing
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[shared_ptr[DataProcessing]])
        Sets the description of the applied processing
        """
        ...
    
    def setChromatogramType(self, type: int ) -> None:
        """
        Cython signature: void setChromatogramType(ChromatogramType type)
        Sets the chromatogram type
        """
        ...
    
    def getChromatogramType(self) -> int:
        """
        Cython signature: ChromatogramType getChromatogramType()
        Get the chromatogram type
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
    
    def __richcmp__(self, other: MSChromatogram, op: int) -> Any:
        ...
    
    def __iter__(self) -> ChromatogramPeak:
       ... 


class MzXMLFile:
    """
    Cython implementation of _MzXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzXMLFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzXMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzXMLFile ) -> None:
        """
        Cython signature: void MzXMLFile(MzXMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None:
        """
        Cython signature: void load(String filename, MSExperiment & exp)
        Loads a MSExperiment from a MzXML file
        
        
        :param exp: MSExperiment
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None:
        """
        Cython signature: void store(String filename, MSExperiment & exp)
        Stores a MSExperiment in a MzXML file
        
        
        :param exp: MSExperiment
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
        Sets options for loading/storing
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


class OPXL_PreprocessedPairSpectra:
    """
    Cython implementation of _OPXL_PreprocessedPairSpectra

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1OPXL_PreprocessedPairSpectra.html>`_
    """
    
    spectra_linear_peaks: MSExperiment
    
    spectra_xlink_peaks: MSExperiment
    
    spectra_all_peaks: MSExperiment
    
    @overload
    def __init__(self, size: int ) -> None:
        """
        Cython signature: void OPXL_PreprocessedPairSpectra(size_t size)
        """
        ...
    
    @overload
    def __init__(self, in_0: OPXL_PreprocessedPairSpectra ) -> None:
        """
        Cython signature: void OPXL_PreprocessedPairSpectra(OPXL_PreprocessedPairSpectra &)
        """
        ... 


class OpenSwathHelper:
    """
    Cython implementation of _OpenSwathHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathHelper.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenSwathHelper()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenSwathHelper ) -> None:
        """
        Cython signature: void OpenSwathHelper(OpenSwathHelper &)
        """
        ...
    
    def checkSwathMapAndSelectTransitions(self, exp: MSExperiment , targeted_exp: TargetedExperiment , transition_exp_used: TargetedExperiment , min_upper_edge_dist: float ) -> bool:
        """
        Cython signature: bool checkSwathMapAndSelectTransitions(MSExperiment & exp, TargetedExperiment & targeted_exp, TargetedExperiment & transition_exp_used, double min_upper_edge_dist)
        """
        ...
    
    def estimateRTRange(self, exp: LightTargetedExperiment ) -> List[float, float]:
        """
        Cython signature: libcpp_pair[double,double] estimateRTRange(LightTargetedExperiment exp)
        Computes the min and max retention time value
        
        Estimate the retention time span of a targeted experiment by returning the min/max values in retention time as a pair
        
        
        :return: A std `pair` that contains (min,max)
        """
        ...
    
    def computePrecursorId(self, transition_group_id: Union[bytes, str, String] , isotope: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String computePrecursorId(const String & transition_group_id, int isotope)
        Computes unique precursor identifier
        
        Uses transition_group_id and isotope number to compute a unique precursor
        id of the form "groupID_Precursor_ix" where x is the isotope number, e.g.
        the monoisotopic precursor would become "groupID_Precursor_i0"
        
        
        :param transition_group_id: Unique id of the transition group (peptide/compound)
        :param isotope: Precursor isotope number
        :return: Unique precursor identifier
        """
        ... 


class PI_PeakArea:
    """
    Cython implementation of _PI_PeakArea

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakArea.html>`_
    """
    
    area: float
    
    height: float
    
    apex_pos: float
    
    hull_points: '_np.ndarray[Any, _np.dtype[_np.float32]]'
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PI_PeakArea()
        """
        ...
    
    @overload
    def __init__(self, in_0: PI_PeakArea ) -> None:
        """
        Cython signature: void PI_PeakArea(PI_PeakArea &)
        """
        ... 


class PI_PeakBackground:
    """
    Cython implementation of _PI_PeakBackground

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakBackground.html>`_
    """
    
    area: float
    
    height: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PI_PeakBackground()
        """
        ...
    
    @overload
    def __init__(self, in_0: PI_PeakBackground ) -> None:
        """
        Cython signature: void PI_PeakBackground(PI_PeakBackground &)
        """
        ... 


class PI_PeakShapeMetrics:
    """
    Cython implementation of _PI_PeakShapeMetrics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakShapeMetrics.html>`_
    """
    
    width_at_5: float
    
    width_at_10: float
    
    width_at_50: float
    
    start_position_at_5: float
    
    start_position_at_10: float
    
    start_position_at_50: float
    
    end_position_at_5: float
    
    end_position_at_10: float
    
    end_position_at_50: float
    
    total_width: float
    
    tailing_factor: float
    
    asymmetry_factor: float
    
    slope_of_baseline: float
    
    baseline_delta_2_height: float
    
    points_across_baseline: int
    
    points_across_half_height: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PI_PeakShapeMetrics()
        """
        ...
    
    @overload
    def __init__(self, in_0: PI_PeakShapeMetrics ) -> None:
        """
        Cython signature: void PI_PeakShapeMetrics(PI_PeakShapeMetrics &)
        """
        ... 


class ParamEntry:
    """
    Cython implementation of _ParamEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Param_1_1ParamEntry.html>`_
    """
    
    name: bytes
    
    description: bytes
    
    value: Union[int, float, bytes, str, List[int], List[float], List[bytes]]
    
    tags: Set[bytes]
    
    valid_strings: List[bytes]
    
    max_float: float
    
    min_float: float
    
    max_int: int
    
    min_int: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ParamEntry()
        """
        ...
    
    @overload
    def __init__(self, in_0: ParamEntry ) -> None:
        """
        Cython signature: void ParamEntry(ParamEntry &)
        """
        ...
    
    @overload
    def __init__(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes , t: List[bytes] ) -> None:
        """
        Cython signature: void ParamEntry(libcpp_string n, ParamValue v, libcpp_string d, libcpp_vector[libcpp_string] t)
        """
        ...
    
    @overload
    def __init__(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes ) -> None:
        """
        Cython signature: void ParamEntry(libcpp_string n, ParamValue v, libcpp_string d)
        """
        ...
    
    def __richcmp__(self, other: ParamEntry, op: int) -> Any:
        ... 


class PeakIntegrator:
    """
    Cython implementation of _PeakIntegrator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakIntegrator.html>`_
      -- Inherits from ['DefaultParamHandler']

    Compute the area, background and shape metrics of a peak
    
    The area computation is performed in integratePeak() and it supports
    integration by simple sum of the intensity, integration by Simpson's rule
    implementations for an odd number of unequally spaced points or integration
    by the trapezoid rule
    
    The background computation is performed in estimateBackground() and it
    supports three different approaches to baseline correction, namely
    computing a rectangular shape under the peak based on the minimum value of
    the peak borders (vertical_division_min), a rectangular shape based on the
    maximum value of the beak borders (vertical_division_max) or a trapezoidal
    shape based on a straight line between the peak borders (base_to_base)
    
    Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
    metrics are supported
    
    The containers supported by the methods are MSChromatogram and MSSpectrum
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakIntegrator()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakIntegrator ) -> None:
        """
        Cython signature: void PeakIntegrator(PeakIntegrator &)
        """
        ...
    
    def getDefaultParameters(self, in_0: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param)
        """
        ...
    
    @overload
    def integratePeak(self, chromatogram: MSChromatogram , left: float , right: float ) -> PI_PeakArea:
        """
        Cython signature: PI_PeakArea integratePeak(MSChromatogram & chromatogram, double left, double right)
        """
        ...
    
    @overload
    def integratePeak(self, spectrum: MSSpectrum , left: float , right: float ) -> PI_PeakArea:
        """
        Cython signature: PI_PeakArea integratePeak(MSSpectrum & spectrum, double left, double right)
        """
        ...
    
    @overload
    def estimateBackground(self, chromatogram: MSChromatogram , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground:
        """
        Cython signature: PI_PeakBackground estimateBackground(MSChromatogram & chromatogram, double left, double right, double peak_apex_pos)
        """
        ...
    
    @overload
    def estimateBackground(self, spectrum: MSSpectrum , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground:
        """
        Cython signature: PI_PeakBackground estimateBackground(MSSpectrum & spectrum, double left, double right, double peak_apex_pos)
        """
        ...
    
    @overload
    def calculatePeakShapeMetrics(self, chromatogram: MSChromatogram , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics:
        """
        Cython signature: PI_PeakShapeMetrics calculatePeakShapeMetrics(MSChromatogram & chromatogram, double left, double right, double peak_height, double peak_apex_pos)
        """
        ...
    
    @overload
    def calculatePeakShapeMetrics(self, spectrum: MSSpectrum , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics:
        """
        Cython signature: PI_PeakShapeMetrics calculatePeakShapeMetrics(MSSpectrum & spectrum, double left, double right, double peak_height, double peak_apex_pos)
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


class __DerivatizationAgent:
    None
    NOT_SELECTED : int
    TBDMS : int
    SIZE_OF_DERIVATIZATIONAGENT : int

    def getMapping(self) -> Dict[int, str]:
       ...
    DerivatizationAgent : __DerivatizationAgent 


class MT_QUANTMETHOD:
    None
    MT_QUANT_AREA : int
    MT_QUANT_MEDIAN : int
    SIZE_OF_MT_QUANTMETHOD : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __MassIntensityType:
    None
    NORM_MAX : int
    NORM_SUM : int
    SIZE_OF_MASSINTENSITYTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ...
    MassIntensityType : __MassIntensityType 

