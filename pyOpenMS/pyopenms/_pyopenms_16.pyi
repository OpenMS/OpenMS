from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_Deisotoper_deisotopeAndSingleCharge(spectra: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool , min_charge: int , max_charge: int , keep_only_deisotoped: bool , min_isopeaks: int , max_isopeaks: int , make_single_charged: bool , annotate_charge: bool , annotate_iso_peak_count: bool , use_decreasing_model: bool , start_intensity_check: int , add_up_intensity: bool , annotate_features: bool ) -> None:
    """
    Cython signature: void deisotopeAndSingleCharge(MSSpectrum & spectra, double fragment_tolerance, bool fragment_unit_ppm, int min_charge, int max_charge, bool keep_only_deisotoped, unsigned int min_isopeaks, unsigned int max_isopeaks, bool make_single_charged, bool annotate_charge, bool annotate_iso_peak_count, bool use_decreasing_model, unsigned int start_intensity_check, bool add_up_intensity, bool annotate_features)
    """
    ...

def __static_Deisotoper_deisotopeAndSingleChargeDefault(spectra: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool ) -> None:
    """
    Cython signature: void deisotopeAndSingleChargeDefault(MSSpectrum & spectra, double fragment_tolerance, bool fragment_unit_ppm)
    """
    ...

def __static_Deisotoper_deisotopeWithAveragineModel(spectrum: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool , number_of_final_peaks: int , min_charge: int , max_charge: int , keep_only_deisotoped: bool , min_isopeaks: int , max_isopeaks: int , make_single_charged: bool , annotate_charge: bool , annotate_iso_peak_count: bool , add_up_intensity: bool ) -> None:
    """
    Cython signature: void deisotopeWithAveragineModel(MSSpectrum & spectrum, double fragment_tolerance, bool fragment_unit_ppm, int number_of_final_peaks, int min_charge, int max_charge, bool keep_only_deisotoped, unsigned int min_isopeaks, unsigned int max_isopeaks, bool make_single_charged, bool annotate_charge, bool annotate_iso_peak_count, bool add_up_intensity)
    """
    ...


class AcquisitionInfo:
    """
    Cython implementation of _AcquisitionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AcquisitionInfo.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AcquisitionInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: AcquisitionInfo ) -> None:
        """
        Cython signature: void AcquisitionInfo(AcquisitionInfo &)
        """
        ...
    
    def getMethodOfCombination(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getMethodOfCombination()
        Returns the method of combination
        """
        ...
    
    def setMethodOfCombination(self, method: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setMethodOfCombination(String method)
        Sets the method of combination
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Number a Acquisition objects
        """
        ...
    
    def __getitem__(self, in_0: int ) -> Acquisition:
        """
        Cython signature: Acquisition & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: Acquisition ) -> None:
        """Cython signature: Acquisition & operator[](size_t)"""
        ...
    
    def push_back(self, in_0: Acquisition ) -> None:
        """
        Cython signature: void push_back(Acquisition)
        Append a Acquisition object
        """
        ...
    
    def resize(self, n: int ) -> None:
        """
        Cython signature: void resize(size_t n)
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
    
    def __richcmp__(self, other: AcquisitionInfo, op: int) -> Any:
        ... 


class BayesianProteinInferenceAlgorithm:
    """
    Cython implementation of _BayesianProteinInferenceAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BayesianProteinInferenceAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Performs a Bayesian protein inference on Protein/Peptide identifications or ConsensusMap.
    
    - Filters for best n PSMs per spectrum.
    - Calculates and filters for best peptide per spectrum.
    - Builds a k-partite graph from the structures.
    - Finds and splits into connected components by DFS
    - Extends the graph by adding layers from indist. protein groups, peptides with the same parents and optionally
      some additional layers (peptide sequence, charge, replicate -> extended model = experimental)
    - Builds a factor graph representation of a Bayesian network using the Evergreen library
      See model param section. It is based on the Fido noisy-OR model with an option for
      regularizing the number of proteins per peptide.
    - Performs loopy belief propagation on the graph and queries protein, protein group and/or peptide posteriors
      See loopy_belief_propagation param section.
    - Learns best parameters via grid search if the parameters were not given in the param section.
    - Writes posteriors to peptides and/or proteins and adds indistinguishable protein groups to the underlying
      data structures.
    - Can make use of OpenMP to parallelize over connected components.
    
    Usage:
    
    .. code-block:: python
    
      from pyopenms import *
      from urllib.request import urlretrieve
      urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/BayesianProteinInference_test.idXML", "BayesianProteinInference_test.idXML")
      proteins = []
      peptides = []
      idf = IdXMLFile()
      idf.load("BayesianProteinInference_test.idXML", proteins, peptides)
      bpia = BayesianProteinInferenceAlgorithm()
      p = bpia.getParameters()
      p.setValue("update_PSM_probabilities", "false")
      bpia.setParameters(p)
      bpia.inferPosteriorProbabilities(proteins, peptides)
      #
      print(len(peptides)) # 9
      print(peptides[0].getHits()[0].getScore()) # 0.6
      print(proteins[0].getHits()[0].getScore()) # 0.624641
      print(proteins[0].getHits()[1].getScore()) # 0.648346
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BayesianProteinInferenceAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, debug_lvl: int ) -> None:
        """
        Cython signature: void BayesianProteinInferenceAlgorithm(unsigned int debug_lvl)
        """
        ...
    
    @overload
    def inferPosteriorProbabilities(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool ) -> None:
        """
        Cython signature: void inferPosteriorProbabilities(libcpp_vector[ProteinIdentification] & proteinIDs, PeptideIdentificationList & peptideIDs, bool greedy_group_resolution)
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        """
        ...
    
    @overload
    def inferPosteriorProbabilities(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None:
        """
        Cython signature: void inferPosteriorProbabilities(libcpp_vector[ProteinIdentification] & proteinIDs, PeptideIdentificationList & peptideIDs, bool greedy_group_resolution, ExperimentalDesign exp_des)
        Writes its results into protein and (optionally also) peptide hits (as new score).
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        Experimental design can be used to create an extended graph with replicate information. (experimental)
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :param exp_des: Experimental Design
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        """
        ...
    
    @overload
    def inferPosteriorProbabilities(self, cmap: ConsensusMap , greedy_group_resolution: bool ) -> None:
        """
        Cython signature: void inferPosteriorProbabilities(ConsensusMap & cmap, bool greedy_group_resolution)
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :return: Writes its protein ID results into the ConsensusMap
        """
        ...
    
    @overload
    def inferPosteriorProbabilities(self, cmap: ConsensusMap , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None:
        """
        Cython signature: void inferPosteriorProbabilities(ConsensusMap & cmap, bool greedy_group_resolution, ExperimentalDesign exp_des)
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs.
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :param exp_des: Experimental Design
        :return: Writes its protein ID results into the ConsensusMap
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


class Deisotoper:
    """
    Cython implementation of _Deisotoper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Deisotoper.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Deisotoper()
        """
        ...
    
    @overload
    def __init__(self, in_0: Deisotoper ) -> None:
        """
        Cython signature: void Deisotoper(Deisotoper &)
        """
        ...
    
    deisotopeAndSingleCharge: __static_Deisotoper_deisotopeAndSingleCharge
    
    deisotopeAndSingleChargeDefault: __static_Deisotoper_deisotopeAndSingleChargeDefault
    
    deisotopeWithAveragineModel: __static_Deisotoper_deisotopeWithAveragineModel 


class ElutionModelFitter:
    """
    Cython implementation of _ElutionModelFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElutionModelFitter.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ElutionModelFitter()
        Helper class for fitting elution models to features
        """
        ...
    
    @overload
    def __init__(self, in_0: ElutionModelFitter ) -> None:
        """
        Cython signature: void ElutionModelFitter(ElutionModelFitter &)
        """
        ...
    
    def fitElutionModels(self, features: FeatureMap ) -> None:
        """
        Cython signature: void fitElutionModels(FeatureMap & features)
        Fit models of elution profiles to all features (and validate them)
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


class FIAMSDataProcessor:
    """
    Cython implementation of _FIAMSDataProcessor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FIAMSDataProcessor.html>`_
      -- Inherits from ['DefaultParamHandler']

      ADD PYTHON DOCUMENTATION HERE
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FIAMSDataProcessor()
        Data processing for FIA-MS data
        """
        ...
    
    @overload
    def __init__(self, in_0: FIAMSDataProcessor ) -> None:
        """
        Cython signature: void FIAMSDataProcessor(FIAMSDataProcessor &)
        """
        ...
    
    def run(self, experiment: MSExperiment , n_seconds: float , output: MzTab , load_cached_spectrum: bool ) -> bool:
        """
        Cython signature: bool run(MSExperiment & experiment, float & n_seconds, MzTab & output, bool load_cached_spectrum)
        Run the full analysis for the experiment for the given time interval\n
        
        The workflow steps are:
        - the time axis of the experiment is cut to the interval from 0 to n_seconds
        - the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
        - data is smoothed by applying the Savitzky-Golay filter
        - peaks are picked
        - the accurate mass search for all the picked peaks is performed
        
        The intermediate summed spectra and picked peaks can be saved to the filesystem.
        Also, the results of the accurate mass search and the signal-to-noise information
        of the resulting spectrum is saved.
        
        
        :param experiment: Input MSExperiment
        :param n_seconds: Input number of seconds
        :param load_cached_spectrum: Load the cached picked spectrum if exists
        :param output: Output of the accurate mass search results
        :return: A boolean indicating if the picked spectrum was loaded from the cached file
        """
        ...
    
    def extractPeaks(self, input_: MSSpectrum ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum extractPeaks(MSSpectrum & input_)
        Pick peaks from the summed spectrum
        
        
        :param input: Input vector of spectra
        :return: A spectrum with picked peaks
        """
        ...
    
    def convertToFeatureMap(self, input_: MSSpectrum ) -> FeatureMap:
        """
        Cython signature: FeatureMap convertToFeatureMap(MSSpectrum & input_)
        Convert a spectrum to a feature map with the corresponding polarity\n
        
        Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`
        
        
        :param input: Input a picked spectrum
        :return: A feature map with the peaks converted to features and polarity from the parameters
        """
        ...
    
    def trackNoise(self, input_: MSSpectrum ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum trackNoise(MSSpectrum & input_)
        Estimate noise for each peak\n
        
        Uses `SignalToNoiseEstimatorMedianRapid`
        
        
        :param input: Input a picked spectrum
        :return: A spectrum object storing logSN information
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


class FeatureDeconvolution:
    """
    Cython implementation of _FeatureDeconvolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureDeconvolution.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureDeconvolution()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureDeconvolution ) -> None:
        """
        Cython signature: void FeatureDeconvolution(FeatureDeconvolution &)
        """
        ...
    
    def compute(self, input: FeatureMap , output: FeatureMap , cmap1: ConsensusMap , cmap2: ConsensusMap ) -> None:
        """
        Cython signature: void compute(FeatureMap & input, FeatureMap & output, ConsensusMap & cmap1, ConsensusMap & cmap2)
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
    CHARGEMODE_FD : __CHARGEMODE_FD 


class ModifiedPeptideGenerator:
    """
    Cython implementation of _ModifiedPeptideGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModifiedPeptideGenerator.html>`_

    Generates modified peptides/proteins.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ModifiedPeptideGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: ModifiedPeptideGenerator ) -> None:
        """
        Cython signature: void ModifiedPeptideGenerator(ModifiedPeptideGenerator &)
        """
        ...
    
    @staticmethod
    def getModifications(modNames: List[bytes] ) -> ModifiedPeptideGenerator_MapToResidueType:
        """
        Cython signature: ModifiedPeptideGenerator_MapToResidueType getModifications(const StringList & modNames)
        """
        ...
    
    @staticmethod
    def applyFixedModifications(fixed_mods: ModifiedPeptideGenerator_MapToResidueType , peptide: AASequence ) -> None:
        """
        Cython signature: void applyFixedModifications(const ModifiedPeptideGenerator_MapToResidueType & fixed_mods, AASequence & peptide)
        """
        ...
    
    @staticmethod
    def applyVariableModifications(var_mods: ModifiedPeptideGenerator_MapToResidueType , peptide: AASequence , max_variable_mods_per_peptide: int , all_modified_peptides: List[AASequence] , keep_original: bool ) -> None:
        """
        Cython signature: void applyVariableModifications(const ModifiedPeptideGenerator_MapToResidueType & var_mods, const AASequence & peptide, size_t max_variable_mods_per_peptide, libcpp_vector[AASequence] & all_modified_peptides, bool keep_original)
        """
        ... 


class ModifiedPeptideGenerator_MapToResidueType:
    """
    Cython implementation of _ModifiedPeptideGenerator_MapToResidueType

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModifiedPeptideGenerator_MapToResidueType.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ModifiedPeptideGenerator_MapToResidueType()
        """
        ...
    
    @overload
    def __init__(self, in_0: ModifiedPeptideGenerator_MapToResidueType ) -> None:
        """
        Cython signature: void ModifiedPeptideGenerator_MapToResidueType(ModifiedPeptideGenerator_MapToResidueType &)
        """
        ... 


class OpenSwathScoring:
    """
    Cython implementation of _OpenSwathScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathScoring.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenSwathScoring()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenSwathScoring ) -> None:
        """
        Cython signature: void OpenSwathScoring(OpenSwathScoring &)
        """
        ...
    
    def initialize(self, rt_normalization_factor: float , add_up_spectra: int , spacing_for_spectra_resampling: float , merge_spectra_by_peak_width_fraction: float , drift_extra: float , su: OpenSwath_Scores_Usage , spectrum_addition_method: bytes , spectrum_merge_method_type: bytes , use_ms1_ion_mobility: bool , apply_im_peak_picking: bool ) -> None:
        """
        Cython signature: void initialize(double rt_normalization_factor, int add_up_spectra, double spacing_for_spectra_resampling, double merge_spectra_by_peak_width_fraction, double drift_extra, OpenSwath_Scores_Usage su, libcpp_string spectrum_addition_method, libcpp_string spectrum_merge_method_type, bool use_ms1_ion_mobility, bool apply_im_peak_picking)
        Initialize the scoring object\n
        Sets the parameters for the scoring
        
        
        :param rt_normalization_factor: Specifies the range of the normalized retention time space
        :param add_up_spectra: How many spectra to add up (default 1)
        :param spacing_for_spectra_resampling: Spacing factor for spectra addition
        :param merge_spectra_by_peak_width_fraction: Fraction of peak width to construct the number of spectra to add
        :param drift_extra: Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)
        :param su: Which scores to actually compute
        :param spectrum_addition_method: Method to use for spectrum addition (valid: "simple", "resample")
        :param spectrum_merge_method_type: Type of method to use for spectrum addition. (valid: "fixed", "dynamic")
        :param use_ms1_ion_mobility: Use MS1 ion mobility extraction in DIA scores
        :param apply_im_peak_picking: Apply peak picking to the  extracted ion mobilograms
        """
        ...
    
    def getNormalized_library_intensities_(self, transitions: List[LightTransition] , normalized_library_intensity: List[float] ) -> None:
        """
        Cython signature: void getNormalized_library_intensities_(libcpp_vector[LightTransition] transitions, libcpp_vector[double] normalized_library_intensity)
        """
        ... 


class OpenSwath_Scores:
    """
    Cython implementation of _OpenSwath_Scores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwath_Scores.html>`_
    """
    
    elution_model_fit_score: float
    
    library_corr: float
    
    library_norm_manhattan: float
    
    library_rootmeansquare: float
    
    library_sangle: float
    
    norm_rt_score: float
    
    isotope_correlation: float
    
    isotope_overlap: float
    
    massdev_score: float
    
    xcorr_coelution_score: float
    
    xcorr_shape_score: float
    
    yseries_score: float
    
    bseries_score: float
    
    log_sn_score: float
    
    weighted_coelution_score: float
    
    weighted_xcorr_shape: float
    
    weighted_massdev_score: float
    
    ms1_xcorr_coelution_score: float
    
    ms1_xcorr_coelution_contrast_score: float
    
    ms1_xcorr_coelution_combined_score: float
    
    ms1_xcorr_shape_score: float
    
    ms1_xcorr_shape_contrast_score: float
    
    ms1_xcorr_shape_combined_score: float
    
    ms1_ppm_score: float
    
    ms1_isotope_correlation: float
    
    ms1_isotope_overlap: float
    
    ms1_mi_score: float
    
    ms1_mi_contrast_score: float
    
    ms1_mi_combined_score: float
    
    library_manhattan: float
    
    library_dotprod: float
    
    intensity: float
    
    total_xic: float
    
    nr_peaks: float
    
    sn_ratio: float
    
    mi_score: float
    
    weighted_mi_score: float
    
    rt_difference: float
    
    normalized_experimental_rt: float
    
    raw_rt_score: float
    
    dotprod_score_dia: float
    
    manhatt_score_dia: float
    
    def __init__(self) -> None:
        """
        Cython signature: void OpenSwath_Scores()
        """
        ...
    
    def get_quick_lda_score(self, library_corr_: float , library_norm_manhattan_: float , norm_rt_score_: float , xcorr_coelution_score_: float , xcorr_shape_score_: float , log_sn_score_: float ) -> float:
        """
        Cython signature: double get_quick_lda_score(double library_corr_, double library_norm_manhattan_, double norm_rt_score_, double xcorr_coelution_score_, double xcorr_shape_score_, double log_sn_score_)
        """
        ...
    
    def calculate_lda_prescore(self, scores: OpenSwath_Scores ) -> float:
        """
        Cython signature: double calculate_lda_prescore(OpenSwath_Scores scores)
        """
        ...
    
    def calculate_swath_lda_prescore(self, scores: OpenSwath_Scores ) -> float:
        """
        Cython signature: double calculate_swath_lda_prescore(OpenSwath_Scores scores)
        """
        ... 


class OpenSwath_Scores_Usage:
    """
    Cython implementation of _OpenSwath_Scores_Usage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwath_Scores_Usage.html>`_
    """
    
    use_coelution_score_: bool
    
    use_shape_score_: bool
    
    use_rt_score_: bool
    
    use_library_score_: bool
    
    use_elution_model_score_: bool
    
    use_intensity_score_: bool
    
    use_total_xic_score_: bool
    
    use_total_mi_score_: bool
    
    use_nr_peaks_score_: bool
    
    use_sn_score_: bool
    
    use_mi_score_: bool
    
    use_dia_scores_: bool
    
    use_ms1_correlation: bool
    
    use_ms1_fullscan: bool
    
    use_ms1_mi: bool
    
    use_uis_scores: bool
    
    def __init__(self) -> None:
        """
        Cython signature: void OpenSwath_Scores_Usage()
        """
        ... 


class Peak2D:
    """
    Cython implementation of _Peak2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Peak2D.html>`_

    A 2-dimensional raw data point or peak.
    
    This data structure is intended for continuous data or peak data.
    If you want to annotated single peaks with meta data, use RichPeak2D instead
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Peak2D()
        """
        ...
    
    @overload
    def __init__(self, in_0: Peak2D ) -> None:
        """
        Cython signature: void Peak2D(Peak2D &)
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
    
    def __richcmp__(self, other: Peak2D, op: int) -> Any:
        ... 


class PeptideIdentification:
    """
    Cython implementation of _PeptideIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIdentification.html>`_
      -- Inherits from ['MetaInfoInterface']

    Represents the peptide hits for a spectrum
    
    This class is closely related to ProteinIdentification, which stores the protein hits
    and the general information about the identification run. More than one PeptideIdentification
    can belong to one ProteinIdentification. The general information about a
    PeptideIdentification has to be looked up in the corresponding ProteinIndentification, using
    the unique `identifier` that links the two.
    When loading PeptideHit instances from a File, the retention time and mass-to-charge ratio
    of the precursor spectrum can be accessed using getRT() and getMZ().
    This information can be used to map the peptide hits to an MSExperiment, a FeatureMap
    or a ConsensusMap using the IDMapper class
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideIdentification()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideIdentification ) -> None:
        """
        Cython signature: void PeptideIdentification(PeptideIdentification &)
        """
        ...
    
    def getHits(self) -> List[PeptideHit]:
        """
        Cython signature: libcpp_vector[PeptideHit] getHits()
        Returns the peptide hits as const
        """
        ...
    
    def insertHit(self, in_0: PeptideHit ) -> None:
        """
        Cython signature: void insertHit(PeptideHit)
        Appends a peptide hit
        """
        ...
    
    def setHits(self, in_0: List[PeptideHit] ) -> None:
        """
        Cython signature: void setHits(libcpp_vector[PeptideHit])
        Sets the peptide hits
        """
        ...
    
    def getSignificanceThreshold(self) -> float:
        """
        Cython signature: double getSignificanceThreshold()
        Returns the peptide significance threshold value
        """
        ...
    
    def setSignificanceThreshold(self, value: float ) -> None:
        """
        Cython signature: void setSignificanceThreshold(double value)
        Setting of the peptide significance threshold value
        """
        ...
    
    def getScoreType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getScoreType()
        """
        ...
    
    def setScoreType(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setScoreType(String)
        """
        ...
    
    def isHigherScoreBetter(self) -> bool:
        """
        Cython signature: bool isHigherScoreBetter()
        """
        ...
    
    def setHigherScoreBetter(self, in_0: bool ) -> None:
        """
        Cython signature: void setHigherScoreBetter(bool)
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        """
        ...
    
    def setIdentifier(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(String)
        """
        ...
    
    def hasMZ(self) -> bool:
        """
        Cython signature: bool hasMZ()
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
    
    def hasRT(self) -> bool:
        """
        Cython signature: bool hasRT()
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        """
        ...
    
    def setRT(self, in_0: float ) -> None:
        """
        Cython signature: void setRT(double)
        """
        ...
    
    def getBaseName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getBaseName()
        """
        ...
    
    def setBaseName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setBaseName(String)
        """
        ...
    
    def getExperimentLabel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getExperimentLabel()
        """
        ...
    
    def setExperimentLabel(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setExperimentLabel(String)
        """
        ...
    
    def sort(self) -> None:
        """
        Cython signature: void sort()
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def getReferencingHits(self, in_0: List[PeptideHit] , in_1: Set[bytes] ) -> List[PeptideHit]:
        """
        Cython signature: libcpp_vector[PeptideHit] getReferencingHits(libcpp_vector[PeptideHit], libcpp_set[String] &)
        Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)
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
    
    def __richcmp__(self, other: PeptideIdentification, op: int) -> Any:
        ... 


class PercolatorFeatureSetHelper:
    """
    Cython implementation of _PercolatorFeatureSetHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorFeatureSetHelper.html>`_

    Percolator feature set and integration helper
    
    This class contains functions to handle (compute, aggregate, integrate)
    Percolator features. This includes the calculation or extraction of
    Percolator features depending on the search engine(s) for later use with
    PercolatorAdapter. It also includes handling the reintegration of the
    percolator result into the set of Identifications
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PercolatorFeatureSetHelper()
        """
        ...
    
    @overload
    def __init__(self, in_0: PercolatorFeatureSetHelper ) -> None:
        """
        Cython signature: void PercolatorFeatureSetHelper(PercolatorFeatureSetHelper &)
        """
        ...
    
    def concatMULTISEPeptideIds(self, all_peptide_ids: PeptideIdentificationList , new_peptide_ids: PeptideIdentificationList , search_engine: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void concatMULTISEPeptideIds(PeptideIdentificationList & all_peptide_ids, PeptideIdentificationList & new_peptide_ids, String search_engine)
        Appends a vector of PeptideIdentification to another and prepares Percolator features in MetaInfo (With the respective key "CONCAT:" + search_engine)
        
        
        :param all_peptide_ids: PeptideIdentification vector to append to
        :param new_peptide_ids: PeptideIdentification vector to be appended
        :param search_engine: Search engine to depend on for feature creation
        """
        ...
    
    def mergeMULTISEPeptideIds(self, all_peptide_ids: PeptideIdentificationList , new_peptide_ids: PeptideIdentificationList , search_engine: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void mergeMULTISEPeptideIds(PeptideIdentificationList & all_peptide_ids, PeptideIdentificationList & new_peptide_ids, String search_engine)
        Merges a vector of PeptideIdentification into another and prepares the merged MetaInfo and scores for collection in addMULTISEFeatures for feature registration
        
        
        :param all_peptide_idsL: PeptideIdentification vector to be merged into
        :param new_peptide_idsL: PeptideIdentification vector to merge
        :param search_engineL: Search engine to create features from their scores
        """
        ...
    
    def mergeMULTISEProteinIds(self, all_protein_ids: List[ProteinIdentification] , new_protein_ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void mergeMULTISEProteinIds(libcpp_vector[ProteinIdentification] & all_protein_ids, libcpp_vector[ProteinIdentification] & new_protein_ids)
        Concatenates SearchParameter of multiple search engine runs and merges PeptideEvidences, collects used search engines in MetaInfo for collection in addMULTISEFeatures for feature registration
        
        
        :param all_protein_ids: ProteinIdentification vector to be merged into
        :param new_protein_ids: ProteinIdentification vector to merge
        """
        ...
    
    def addMSGFFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None:
        """
        Cython signature: void addMSGFFeatures(PeptideIdentificationList & peptide_ids, StringList & feature_set)
        Creates and adds MSGF+ specific Percolator features and registers them in feature_set. MSGF+ should be run with the addFeatures flag enabled
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        ...
    
    def addXTANDEMFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None:
        """
        Cython signature: void addXTANDEMFeatures(PeptideIdentificationList & peptide_ids, StringList & feature_set)
        Creates and adds X!Tandem specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        ...
    
    def addCOMETFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None:
        """
        Cython signature: void addCOMETFeatures(PeptideIdentificationList & peptide_ids, StringList & feature_set)
        Creates and adds Comet specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        ...
    
    def addMASCOTFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None:
        """
        Cython signature: void addMASCOTFeatures(PeptideIdentificationList & peptide_ids, StringList & feature_set)
        Creates and adds Mascot specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        ...
    
    def addMULTISEFeatures(self, peptide_ids: PeptideIdentificationList , search_engines_used: List[bytes] , feature_set: List[bytes] , complete_only: bool , limits_imputation: bool ) -> None:
        """
        Cython signature: void addMULTISEFeatures(PeptideIdentificationList & peptide_ids, StringList & search_engines_used, StringList & feature_set, bool complete_only, bool limits_imputation)
        Adds multiple search engine specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param search_engines_used: The list of search engines to be considered
        :param feature_set: Register of added features
        :param complete_only: Will only add features for PeptideIdentifications where all given search engines identified something
        :param limits_imputation: Uses C++ numeric limits as imputed values instead of min/max of that feature
        """
        ...
    
    def addCONCATSEFeatures(self, peptide_id_list: PeptideIdentificationList , search_engines_used: List[bytes] , feature_set: List[bytes] ) -> None:
        """
        Cython signature: void addCONCATSEFeatures(PeptideIdentificationList & peptide_id_list, StringList & search_engines_used, StringList & feature_set)
        Adds multiple search engine specific Percolator features and registers them in feature_set
        
        This struct can be used to store both peak or feature indices
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param search_engines_used: The list of search engines to be considered
        :param feature_set: Register of added features
        """
        ...
    
    def checkExtraFeatures(self, psms: List[PeptideHit] , extra_features: List[bytes] ) -> None:
        """
        Cython signature: void checkExtraFeatures(libcpp_vector[PeptideHit] & psms, StringList & extra_features)
        Checks and removes requested extra Percolator features that are actually unavailable (to compute)
        
        
        :param psms: The vector of PeptideHit to be checked
        :param extra_features: The list of requested extra features
        """
        ... 


class ResidueModification:
    """
    Cython implementation of _ResidueModification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ResidueModification.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ResidueModification()
        """
        ...
    
    @overload
    def __init__(self, in_0: ResidueModification ) -> None:
        """
        Cython signature: void ResidueModification(ResidueModification &)
        """
        ...
    
    def setId(self, id_: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setId(const String & id_)
        Sets the identifier of the modification
        """
        ...
    
    def getId(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getId()
        Returns the identifier of the modification
        """
        ...
    
    def setFullId(self, full_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFullId(const String & full_id)
        Sets the full identifier (Unimod Accession + origin, if available)
        """
        ...
    
    def getFullId(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFullId()
        """
        ...
    
    def getUniModRecordId(self) -> int:
        """
        Cython signature: int getUniModRecordId()
        Gets the unimod record id
        """
        ...
    
    def setUniModRecordId(self, id_: int ) -> None:
        """
        Cython signature: void setUniModRecordId(int id_)
        Sets the unimod record id
        """
        ...
    
    def getUniModAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getUniModAccession()
        Returns the unimod accession if available
        """
        ...
    
    def setPSIMODAccession(self, id_: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPSIMODAccession(const String & id_)
        Sets the MOD-XXXXX accession of PSI-MOD
        """
        ...
    
    def getPSIMODAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPSIMODAccession()
        Returns the PSI-MOD accession if available
        """
        ...
    
    def setFullName(self, full_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFullName(const String & full_name)
        Sets the full name of the modification; must NOT contain the origin (or . for terminals!)
        """
        ...
    
    def getFullName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFullName()
        Returns the full name of the modification
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String & name)
        Sets the name of modification
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the PSI-MS-label if available; e.g. Mascot uses this name
        """
        ...
    
    @overload
    def setTermSpecificity(self, term_spec: int ) -> None:
        """
        Cython signature: void setTermSpecificity(TermSpecificity term_spec)
        Sets the term specificity
        """
        ...
    
    @overload
    def setTermSpecificity(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTermSpecificity(const String & name)
        Sets the terminal specificity using a name
        """
        ...
    
    def getTermSpecificity(self) -> int:
        """
        Cython signature: TermSpecificity getTermSpecificity()
        Returns terminal specificity
        """
        ...
    
    def getTermSpecificityName(self, in_0: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getTermSpecificityName(TermSpecificity)
        Returns the name of the terminal specificity
        """
        ...
    
    def setOrigin(self, origin: bytes ) -> None:
        """
        Cython signature: void setOrigin(char origin)
        Sets the origin (i.e. modified amino acid)
        """
        ...
    
    def getOrigin(self) -> bytes:
        """
        Cython signature: char getOrigin()
        Returns the origin (i.e. modified amino acid)
        """
        ...
    
    @overload
    def setSourceClassification(self, classification: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSourceClassification(const String & classification)
        Classification as defined by the PSI-MOD
        """
        ...
    
    @overload
    def setSourceClassification(self, classification: int ) -> None:
        """
        Cython signature: void setSourceClassification(SourceClassification classification)
        Sets the source classification
        """
        ...
    
    def getSourceClassification(self) -> int:
        """
        Cython signature: SourceClassification getSourceClassification()
        Returns the source classification, if none was set, it is unspecific
        """
        ...
    
    def getSourceClassificationName(self, classification: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getSourceClassificationName(SourceClassification classification)
        Returns the classification
        """
        ...
    
    def setAverageMass(self, mass: float ) -> None:
        """
        Cython signature: void setAverageMass(double mass)
        Sets the average mass
        """
        ...
    
    def getAverageMass(self) -> float:
        """
        Cython signature: double getAverageMass()
        Returns the average mass if set
        """
        ...
    
    def setMonoMass(self, mass: float ) -> None:
        """
        Cython signature: void setMonoMass(double mass)
        Sets the monoisotopic mass (this must include the weight of the residue itself!)
        """
        ...
    
    def getMonoMass(self) -> float:
        """
        Cython signature: double getMonoMass()
        Return the monoisotopic mass, or 0.0 if not set
        """
        ...
    
    def setDiffAverageMass(self, mass: float ) -> None:
        """
        Cython signature: void setDiffAverageMass(double mass)
        Sets the difference average mass
        """
        ...
    
    def getDiffAverageMass(self) -> float:
        """
        Cython signature: double getDiffAverageMass()
        Returns the difference average mass, or 0.0 if not set
        """
        ...
    
    def setDiffMonoMass(self, mass: float ) -> None:
        """
        Cython signature: void setDiffMonoMass(double mass)
        Sets the difference monoisotopic mass
        """
        ...
    
    def getDiffMonoMass(self) -> float:
        """
        Cython signature: double getDiffMonoMass()
        Returns the diff monoisotopic mass, or 0.0 if not set
        """
        ...
    
    def setFormula(self, composition: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFormula(const String & composition)
        Sets the formula (no masses will be changed)
        """
        ...
    
    def getFormula(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFormula()
        Returns the chemical formula if set
        """
        ...
    
    def setDiffFormula(self, diff_formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void setDiffFormula(EmpiricalFormula & diff_formula)
        Sets diff formula (no masses will be changed)
        """
        ...
    
    def getDiffFormula(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getDiffFormula()
        Returns the diff formula if one was set
        """
        ...
    
    def setSynonyms(self, synonyms: Set[bytes] ) -> None:
        """
        Cython signature: void setSynonyms(libcpp_set[String] & synonyms)
        Sets the synonyms of that modification
        """
        ...
    
    def addSynonym(self, synonym: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addSynonym(const String & synonym)
        Adds a synonym to the unique list
        """
        ...
    
    def getSynonyms(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getSynonyms()
        Returns the set of synonyms
        """
        ...
    
    def setNeutralLossDiffFormulas(self, diff_formulas: List[EmpiricalFormula] ) -> None:
        """
        Cython signature: void setNeutralLossDiffFormulas(libcpp_vector[EmpiricalFormula] & diff_formulas)
        Sets the neutral loss formula
        """
        ...
    
    def getNeutralLossDiffFormulas(self) -> List[EmpiricalFormula]:
        """
        Cython signature: libcpp_vector[EmpiricalFormula] getNeutralLossDiffFormulas()
        Returns the neutral loss diff formula (if available)
        """
        ...
    
    def setNeutralLossMonoMasses(self, mono_masses: List[float] ) -> None:
        """
        Cython signature: void setNeutralLossMonoMasses(libcpp_vector[double] mono_masses)
        Sets the neutral loss mono weight
        """
        ...
    
    def getNeutralLossMonoMasses(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getNeutralLossMonoMasses()
        Returns the neutral loss mono weight
        """
        ...
    
    def setNeutralLossAverageMasses(self, average_masses: List[float] ) -> None:
        """
        Cython signature: void setNeutralLossAverageMasses(libcpp_vector[double] average_masses)
        Sets the neutral loss average weight
        """
        ...
    
    def getNeutralLossAverageMasses(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getNeutralLossAverageMasses()
        Returns the neutral loss average weight
        """
        ...
    
    def hasNeutralLoss(self) -> bool:
        """
        Cython signature: bool hasNeutralLoss()
        Returns true if a neutral loss formula is set
        """
        ...
    
    def isUserDefined(self) -> bool:
        """
        Cython signature: bool isUserDefined()
        Returns true if it is a user-defined modification (empty id)
        """
        ...
    
    def __richcmp__(self, other: ResidueModification, op: int) -> Any:
        ...
    SourceClassification : __SourceClassification
    TermSpecificity : __TermSpecificity 


class SimpleSearchEngineAlgorithm:
    """
    Cython implementation of _SimpleSearchEngineAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimpleSearchEngineAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SimpleSearchEngineAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: SimpleSearchEngineAlgorithm ) -> None:
        """
        Cython signature: void SimpleSearchEngineAlgorithm(SimpleSearchEngineAlgorithm &)
        """
        ...
    
    def search(self, in_mzML: Union[bytes, str, String] , in_db: Union[bytes, str, String] , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void search(const String & in_mzML, const String & in_db, libcpp_vector[ProteinIdentification] & prot_ids, PeptideIdentificationList & pep_ids)
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


class SpectrumAccessOpenMS:
    """
    Cython implementation of _SpectrumAccessOpenMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMS.html>`_
      -- Inherits from ['ISpectrumAccess']

    An implementation of the OpenSWATH Spectrum Access interface using OpenMS
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMS()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMS ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMS(SpectrumAccessOpenMS &)
        """
        ...
    
    @overload
    def __init__(self, ms_experiment: MSExperiment ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMS(shared_ptr[MSExperiment] & ms_experiment)
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


class TM_DataPoint:
    """
    Cython implementation of _TM_DataPoint

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TM_DataPoint.html>`_
    """
    
    first: float
    
    second: float
    
    note: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TM_DataPoint()
        """
        ...
    
    @overload
    def __init__(self, in_0: float , in_1: float ) -> None:
        """
        Cython signature: void TM_DataPoint(double, double)
        """
        ...
    
    @overload
    def __init__(self, in_0: float , in_1: float , in_2: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void TM_DataPoint(double, double, const String &)
        """
        ...
    
    def __richcmp__(self, other: TM_DataPoint, op: int) -> Any:
        ... 


class __CHARGEMODE_FD:
    None
    QFROMFEATURE : int
    QHEURISTIC : int
    QALL : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class DimensionDescription:
    None
    RT : int
    MZ : int
    DIMENSION : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __SourceClassification:
    None
    ARTIFACT : int
    HYPOTHETICAL : int
    NATURAL : int
    POSTTRANSLATIONAL : int
    MULTIPLE : int
    CHEMICAL_DERIVATIVE : int
    ISOTOPIC_LABEL : int
    PRETRANSLATIONAL : int
    OTHER_GLYCOSYLATION : int
    NLINKED_GLYCOSYLATION : int
    AA_SUBSTITUTION : int
    OTHER : int
    NONSTANDARD_RESIDUE : int
    COTRANSLATIONAL : int
    OLINKED_GLYCOSYLATION : int
    UNKNOWN : int
    NUMBER_OF_SOURCE_CLASSIFICATIONS : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __TermSpecificity:
    None
    ANYWHERE : int
    C_TERM : int
    N_TERM : int
    PROTEIN_C_TERM : int
    PROTEIN_N_TERM : int
    NUMBER_OF_TERM_SPECIFICITY : int

    def getMapping(self) -> Dict[int, str]:
       ... 

