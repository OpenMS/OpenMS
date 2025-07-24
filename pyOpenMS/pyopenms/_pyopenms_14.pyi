from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class BiGaussFitter1D:
    """
    Cython implementation of _BiGaussFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BiGaussFitter1D.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BiGaussFitter1D()
        """
        ...
    
    @overload
    def __init__(self, in_0: BiGaussFitter1D ) -> None:
        """
        Cython signature: void BiGaussFitter1D(BiGaussFitter1D &)
        """
        ... 


class EmgScoring:
    """
    Cython implementation of _EmgScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgScoring.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EmgScoring()
        Helps in scoring of an elution peak using an exponentially modified gaussian distribution model
        """
        ...
    
    @overload
    def __init__(self, in_0: EmgScoring ) -> None:
        """
        Cython signature: void EmgScoring(EmgScoring &)
        """
        ...
    
    def setFitterParam(self, param: Param ) -> None:
        """
        Cython signature: void setFitterParam(Param param)
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        """
        ...
    
    def elutionModelFit(self, current_section: '_np.ndarray[Any, _np.dtype[_np.float32]]' , smooth_data: bool ) -> float:
        """
        Cython signature: double elutionModelFit(libcpp_vector[DPosition2] current_section, bool smooth_data)
        """
        ... 


class GaussFitResult:
    """
    Cython implementation of _GaussFitResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1GaussFitResult.html>`_
    """
    
    A: float
    
    x0: float
    
    sigma: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GaussFitResult()
        """
        ...
    
    @overload
    def __init__(self, in_0: float , in_1: float , in_2: float ) -> None:
        """
        Cython signature: void GaussFitResult(double, double, double)
        """
        ...
    
    @overload
    def __init__(self, in_0: GaussFitResult ) -> None:
        """
        Cython signature: void GaussFitResult(GaussFitResult &)
        """
        ...
    
    def eval(self, in_0: float ) -> float:
        """
        Cython signature: double eval(double)
        """
        ... 


class GaussFitter:
    """
    Cython implementation of _GaussFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1GaussFitter.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void GaussFitter()
        Implements a fitter for Gaussian functions
        """
        ...
    
    def setInitialParameters(self, result: GaussFitResult ) -> None:
        """
        Cython signature: void setInitialParameters(GaussFitResult & result)
        Sets the initial parameters used by the fit method as initial guess for the Gaussian
        """
        ...
    
    def fit(self, points: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> GaussFitResult:
        """
        Cython signature: GaussFitResult fit(libcpp_vector[DPosition2] points)
        Fits a Gaussian distribution to the given data points
        """
        ... 


class LowessSmoothing:
    """
    Cython implementation of _LowessSmoothing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LowessSmoothing.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void LowessSmoothing()
        """
        ...
    
    def smoothData(self, x: List[float] , y: List[float] , y_smoothed: List[float] ) -> None:
        """
        Cython signature: void smoothData(libcpp_vector[double] x, libcpp_vector[double] y, libcpp_vector[double] & y_smoothed)
        Smoothing method that receives x and y coordinates (e.g., RT and intensities) and computes smoothed intensities
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


class MRMFQC_ComponentGroupPairQCs:
    """
    Cython implementation of _MRMFQC_ComponentGroupPairQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentGroupPairQCs.html>`_
    """
    
    component_group_name: Union[bytes, str, String]
    
    resolution_pair_name: Union[bytes, str, String]
    
    resolution_l: float
    
    resolution_u: float
    
    rt_diff_l: float
    
    rt_diff_u: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFQC_ComponentGroupPairQCs()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFQC_ComponentGroupPairQCs ) -> None:
        """
        Cython signature: void MRMFQC_ComponentGroupPairQCs(MRMFQC_ComponentGroupPairQCs &)
        """
        ... 


class MRMFQC_ComponentGroupQCs:
    """
    Cython implementation of _MRMFQC_ComponentGroupQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentGroupQCs.html>`_
    """
    
    component_group_name: Union[bytes, str, String]
    
    retention_time_l: float
    
    retention_time_u: float
    
    intensity_l: float
    
    intensity_u: float
    
    overall_quality_l: float
    
    overall_quality_u: float
    
    n_heavy_l: int
    
    n_heavy_u: int
    
    n_light_l: int
    
    n_light_u: int
    
    n_detecting_l: int
    
    n_detecting_u: int
    
    n_quantifying_l: int
    
    n_quantifying_u: int
    
    n_identifying_l: int
    
    n_identifying_u: int
    
    n_transitions_l: int
    
    n_transitions_u: int
    
    ion_ratio_pair_name_1: Union[bytes, str, String]
    
    ion_ratio_pair_name_2: Union[bytes, str, String]
    
    ion_ratio_l: float
    
    ion_ratio_u: float
    
    ion_ratio_feature_name: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFQC_ComponentGroupQCs()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFQC_ComponentGroupQCs ) -> None:
        """
        Cython signature: void MRMFQC_ComponentGroupQCs(MRMFQC_ComponentGroupQCs &)
        """
        ... 


class MRMFQC_ComponentQCs:
    """
    Cython implementation of _MRMFQC_ComponentQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentQCs.html>`_
    """
    
    component_name: Union[bytes, str, String]
    
    retention_time_l: float
    
    retention_time_u: float
    
    intensity_l: float
    
    intensity_u: float
    
    overall_quality_l: float
    
    overall_quality_u: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFQC_ComponentQCs()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFQC_ComponentQCs ) -> None:
        """
        Cython signature: void MRMFQC_ComponentQCs(MRMFQC_ComponentQCs &)
        """
        ... 


class MRMFeatureQC:
    """
    Cython implementation of _MRMFeatureQC

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureQC.html>`_
    """
    
    component_qcs: List[MRMFQC_ComponentQCs]
    
    component_group_qcs: List[MRMFQC_ComponentGroupQCs]
    
    component_group_pair_qcs: List[MRMFQC_ComponentGroupPairQCs]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeatureQC()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeatureQC ) -> None:
        """
        Cython signature: void MRMFeatureQC(MRMFeatureQC &)
        """
        ... 


class MRMMapping:
    """
    Cython implementation of _MRMMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMMapping.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MRMMapping()
        """
        ...
    
    def mapExperiment(self, input_chromatograms: MSExperiment , targeted_exp: TargetedExperiment , output: MSExperiment ) -> None:
        """
        Cython signature: void mapExperiment(MSExperiment input_chromatograms, TargetedExperiment targeted_exp, MSExperiment & output)
        Maps input chromatograms to assays in a targeted experiment
        
        The output chromatograms are an annotated copy of the input chromatograms
        with native id, precursor information and peptide sequence (if available)
        annotated in the chromatogram files
        
        The algorithm tries to match a given set of chromatograms and targeted
        assays. It iterates through all the chromatograms retrieves one or more
        matching targeted assay for the chromatogram. By default, the algorithm
        assumes that a 1:1 mapping exists. If a chromatogram cannot be mapped
        (does not have a corresponding assay) the algorithm issues a warning, the
        user can specify that the program should abort in such a case (see
        error_on_unmapped)
        
        :note If multiple mapping is enabled (see map_multiple_assays parameter)
        then each mapped assay will get its own chromatogram that contains the
        same raw data but different meta-annotation. This *can* be useful if the
        same transition is used to monitor multiple analytes but may also
        indicate a problem with too wide mapping tolerances
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


class MapAlignmentAlgorithmPoseClustering:
    """
    Cython implementation of _MapAlignmentAlgorithmPoseClustering

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmPoseClustering.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MapAlignmentAlgorithmPoseClustering()
        """
        ...
    
    @overload
    def align(self, in_0: FeatureMap , in_1: TransformationDescription ) -> None:
        """
        Cython signature: void align(FeatureMap, TransformationDescription &)
        """
        ...
    
    @overload
    def align(self, in_0: MSExperiment , in_1: TransformationDescription ) -> None:
        """
        Cython signature: void align(MSExperiment, TransformationDescription &)
        """
        ...
    
    @overload
    def setReference(self, in_0: FeatureMap ) -> None:
        """
        Cython signature: void setReference(FeatureMap)
        Sets the reference for the alignment
        """
        ...
    
    @overload
    def setReference(self, in_0: MSExperiment ) -> None:
        """
        Cython signature: void setReference(MSExperiment)
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


class MassExplainer:
    """
    Cython implementation of _MassExplainer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassExplainer.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassExplainer()
        Computes empirical formulas for given mass differences using a set of allowed elements
        """
        ...
    
    @overload
    def __init__(self, in_0: MassExplainer ) -> None:
        """
        Cython signature: void MassExplainer(MassExplainer &)
        """
        ...
    
    @overload
    def __init__(self, adduct_base: List[Adduct] ) -> None:
        """
        Cython signature: void MassExplainer(libcpp_vector[Adduct] adduct_base)
        """
        ...
    
    @overload
    def __init__(self, q_min: int , q_max: int , max_span: int , thresh_logp: float ) -> None:
        """
        Cython signature: void MassExplainer(int q_min, int q_max, int max_span, double thresh_logp)
        """
        ...
    
    def setAdductBase(self, adduct_base: List[Adduct] ) -> None:
        """
        Cython signature: void setAdductBase(libcpp_vector[Adduct] adduct_base)
        Sets the set of possible adducts
        """
        ...
    
    def getAdductBase(self) -> List[Adduct]:
        """
        Cython signature: libcpp_vector[Adduct] getAdductBase()
        Returns the set of adducts
        """
        ...
    
    def getCompomerById(self, id: int ) -> Compomer:
        """
        Cython signature: Compomer getCompomerById(size_t id)
        Returns a compomer by its Id (useful after a query() )
        """
        ...
    
    def compute(self) -> None:
        """
        Cython signature: void compute()
        Fill map with possible mass-differences along with their explanation
        """
        ... 


class PrecursorPurity:
    """
    Cython implementation of _PrecursorPurity

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PrecursorPurity.html>`_

    Precursor purity or noise estimation
    
    This class computes metrics for precursor isolation window purity (or noise)
    The function extracts the peaks from an isolation window targeted for fragmentation
    and determines which peaks are isotopes of the target and which come from other sources
    The intensities of the assumed target peaks are summed up as the target intensity
    Using this information it calculates an intensity ratio for the relative intensity of the target
    compared to other sources
    These metrics are combined over the previous and the next MS1 spectrum
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PrecursorPurity()
        """
        ...
    
    @overload
    def __init__(self, in_0: PrecursorPurity ) -> None:
        """
        Cython signature: void PrecursorPurity(PrecursorPurity &)
        """
        ...
    
    def computePrecursorPurity(self, ms1: MSSpectrum , pre: Precursor , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool ) -> PurityScores:
        """
        Cython signature: PurityScores computePrecursorPurity(MSSpectrum ms1, Precursor pre, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm)
        Compute precursor purity metrics for one MS2 precursor
        
        Note: This function is implemented in a general way and can also be used for e.g. MS3 precursor isolation windows in MS2 spectra
        Spectra annotated with charge 0 will be treated as charge 1.
        
        
        :param ms1: The Spectrum containing the isolation window
        :param pre: The precursor containing the definition the isolation window
        :param precursor_mass_tolerance: The precursor tolerance. Is used for determining the targeted peak and deisotoping
        :param precursor_mass_tolerance_unit_ppm: The unit of the precursor tolerance
        """
        ... 


class ProteinGroup:
    """
    Cython implementation of _ProteinGroup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinGroup.html>`_
    """
    
    probability: float
    
    accessions: List[bytes]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteinGroup()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteinGroup ) -> None:
        """
        Cython signature: void ProteinGroup(ProteinGroup &)
        """
        ... 


class ProteinIdentification:
    """
    Cython implementation of _ProteinIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinIdentification.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteinIdentification()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteinIdentification ) -> None:
        """
        Cython signature: void ProteinIdentification(ProteinIdentification &)
        """
        ...
    
    def getHits(self) -> List[ProteinHit]:
        """
        Cython signature: libcpp_vector[ProteinHit] getHits()
        Returns the protein hits
        """
        ...
    
    def insertHit(self, input: ProteinHit ) -> None:
        """
        Cython signature: void insertHit(ProteinHit input)
        Appends a protein hit
        """
        ...
    
    def setHits(self, hits: List[ProteinHit] ) -> None:
        """
        Cython signature: void setHits(libcpp_vector[ProteinHit] hits)
        Sets the protein hits
        """
        ...
    
    def getProteinGroups(self) -> List[ProteinGroup]:
        """
        Cython signature: libcpp_vector[ProteinGroup] getProteinGroups()
        Returns the protein groups
        """
        ...
    
    def insertProteinGroup(self, group: ProteinGroup ) -> None:
        """
        Cython signature: void insertProteinGroup(ProteinGroup group)
        Appends a new protein group
        """
        ...
    
    def getIndistinguishableProteins(self) -> List[ProteinGroup]:
        """
        Cython signature: libcpp_vector[ProteinGroup] getIndistinguishableProteins()
        Returns the indistinguishable proteins
        """
        ...
    
    def insertIndistinguishableProteins(self, group: ProteinGroup ) -> None:
        """
        Cython signature: void insertIndistinguishableProteins(ProteinGroup group)
        Appends new indistinguishable proteins
        """
        ...
    
    def getSignificanceThreshold(self) -> float:
        """
        Cython signature: double getSignificanceThreshold()
        Returns the protein significance threshold value
        """
        ...
    
    def setSignificanceThreshold(self, value: float ) -> None:
        """
        Cython signature: void setSignificanceThreshold(double value)
        Sets the protein significance threshold value
        """
        ...
    
    def getScoreType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getScoreType()
        Returns the protein score type
        """
        ...
    
    def setScoreType(self, type: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setScoreType(String type)
        Sets the protein score type
        """
        ...
    
    def isHigherScoreBetter(self) -> bool:
        """
        Cython signature: bool isHigherScoreBetter()
        Returns true if a higher score represents a better score
        """
        ...
    
    def setHigherScoreBetter(self, higher_is_better: bool ) -> None:
        """
        Cython signature: void setHigherScoreBetter(bool higher_is_better)
        Sets the orientation of the score (is higher better?)
        """
        ...
    
    def sort(self) -> None:
        """
        Cython signature: void sort()
        Sorts the protein hits according to their score
        """
        ...
    
    def computeCoverage(self, pep_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void computeCoverage(PeptideIdentificationList pep_ids)
        Compute the coverage (in percent) of all ProteinHits given PeptideHits
        """
        ...
    
    def getDateTime(self) -> DateTime:
        """
        Cython signature: DateTime getDateTime()
        Returns the date of the protein identification run
        """
        ...
    
    def setDateTime(self, date: DateTime ) -> None:
        """
        Cython signature: void setDateTime(DateTime date)
        Sets the date of the protein identification run
        """
        ...
    
    def setSearchEngine(self, search_engine: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSearchEngine(String search_engine)
        Sets the search engine type
        """
        ...
    
    def getSearchEngine(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSearchEngine()
        Returns the type of search engine used
        """
        ...
    
    def setSearchEngineVersion(self, search_engine_version: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSearchEngineVersion(String search_engine_version)
        Sets the search engine version
        """
        ...
    
    def getSearchEngineVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSearchEngineVersion()
        Returns the search engine version
        """
        ...
    
    def setSearchParameters(self, search_parameters: SearchParameters ) -> None:
        """
        Cython signature: void setSearchParameters(SearchParameters search_parameters)
        Sets the search parameters
        """
        ...
    
    def getSearchParameters(self) -> SearchParameters:
        """
        Cython signature: SearchParameters getSearchParameters()
        Returns the search parameters
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        Returns the identifier
        """
        ...
    
    def setIdentifier(self, id_: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(String id_)
        Sets the identifier
        """
        ...
    
    @overload
    def setPrimaryMSRunPath(self, s: List[bytes] ) -> None:
        """
        Cython signature: void setPrimaryMSRunPath(StringList & s)
        Set the file paths to the primary MS runs (usually the mzML files obtained after data conversion from raw files)
        
        
        :param raw: Store paths to the raw files (or equivalent) rather than mzMLs
        """
        ...
    
    @overload
    def setPrimaryMSRunPath(self, s: List[bytes] , raw: bool ) -> None:
        """
        Cython signature: void setPrimaryMSRunPath(StringList & s, bool raw)
        """
        ...
    
    @overload
    def addPrimaryMSRunPath(self, s: List[bytes] ) -> None:
        """
        Cython signature: void addPrimaryMSRunPath(StringList & s)
        """
        ...
    
    @overload
    def addPrimaryMSRunPath(self, s: List[bytes] , raw: bool ) -> None:
        """
        Cython signature: void addPrimaryMSRunPath(StringList & s, bool raw)
        """
        ...
    
    @overload
    def getPrimaryMSRunPath(self, output: List[bytes] ) -> None:
        """
        Cython signature: void getPrimaryMSRunPath(StringList & output)
        """
        ...
    
    @overload
    def getPrimaryMSRunPath(self, output: List[bytes] , raw: bool ) -> None:
        """
        Cython signature: void getPrimaryMSRunPath(StringList & output, bool raw)
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
    
    def __richcmp__(self, other: ProteinIdentification, op: int) -> Any:
        ...
    PeakMassType : __PeakMassType 


class PurityScores:
    """
    Cython implementation of _PurityScores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PurityScores.html>`_
    """
    
    total_intensity: float
    
    target_intensity: float
    
    signal_proportion: float
    
    target_peak_count: int
    
    interfering_peak_count: int
    
    interfering_peaks: MSSpectrum
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PurityScores()
        """
        ...
    
    @overload
    def __init__(self, in_0: PurityScores ) -> None:
        """
        Cython signature: void PurityScores(PurityScores &)
        """
        ... 


class SearchParameters:
    """
    Cython implementation of _SearchParameters

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SearchParameters.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    db: Union[bytes, str, String]
    
    db_version: Union[bytes, str, String]
    
    taxonomy: Union[bytes, str, String]
    
    charges: Union[bytes, str, String]
    
    mass_type: int
    
    fixed_modifications: List[bytes]
    
    variable_modifications: List[bytes]
    
    missed_cleavages: int
    
    fragment_mass_tolerance: float
    
    fragment_mass_tolerance_ppm: bool
    
    precursor_mass_tolerance: float
    
    precursor_mass_tolerance_ppm: bool
    
    digestion_enzyme: DigestionEnzymeProtein
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SearchParameters()
        """
        ...
    
    @overload
    def __init__(self, in_0: SearchParameters ) -> None:
        """
        Cython signature: void SearchParameters(SearchParameters &)
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
    
    def __richcmp__(self, other: SearchParameters, op: int) -> Any:
        ... 


class SpectrumSettings:
    """
    Cython implementation of _SpectrumSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumSettings.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumSettings()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumSettings ) -> None:
        """
        Cython signature: void SpectrumSettings(SpectrumSettings &)
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
    
    def __richcmp__(self, other: SpectrumSettings, op: int) -> Any:
        ...
    SpectrumType : __SpectrumType 


class SwathMap:
    """
    Cython implementation of _SwathMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1SwathMap.html>`_
    """
    
    lower: float
    
    upper: float
    
    center: float
    
    ms1: bool
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SwathMap()
        Data structure to hold one SWATH map with information about upper / lower isolation window and whether the map is MS1 or MS2
        """
        ...
    
    @overload
    def __init__(self, in_0: SwathMap ) -> None:
        """
        Cython signature: void SwathMap(SwathMap &)
        """
        ...
    
    @overload
    def __init__(self, mz_start: float , mz_end: float , mz_center: float , is_ms1: bool ) -> None:
        """
        Cython signature: void SwathMap(double mz_start, double mz_end, double mz_center, bool is_ms1)
        """
        ... 


class streampos:
    """
    Cython implementation of _streampos

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classstd_1_1streampos.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void streampos()
        """
        ...
    
    @overload
    def __init__(self, in_0: streampos ) -> None:
        """
        Cython signature: void streampos(streampos &)
        """
        ... 


class __PeakMassType:
    None
    MONOISOTOPIC : int
    AVERAGE : int
    SIZE_OF_PEAKMASSTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __SpectrumType:
    None
    UNKNOWN : int
    CENTROID : int
    PROFILE : int
    SIZE_OF_SPECTRUMTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

