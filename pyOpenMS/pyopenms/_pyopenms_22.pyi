from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_PercolatorInfile_store(pin_file: Union[bytes, str, String] , peptide_ids: PeptideIdentificationList , feature_set: List[bytes] , in_3: bytes , min_charge: int , max_charge: int ) -> None:
    """
    Cython signature: void store(String pin_file, PeptideIdentificationList peptide_ids, StringList feature_set, libcpp_string, int min_charge, int max_charge)
    """
    ...


class ChargedIndexSet:
    """
    Cython implementation of _ChargedIndexSet

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChargedIndexSet.html>`_
    """
    
    charge: int
    
    def __init__(self) -> None:
        """
        Cython signature: void ChargedIndexSet()
        Index set with associated charge estimate
        """
        ... 


class EmgGradientDescent:
    """
    Cython implementation of _EmgGradientDescent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgGradientDescent.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EmgGradientDescent()
        Compute the area, background and shape metrics of a peak
        """
        ...
    
    @overload
    def __init__(self, in_0: EmgGradientDescent ) -> None:
        """
        Cython signature: void EmgGradientDescent(EmgGradientDescent &)
        """
        ...
    
    def getDefaultParameters(self, in_0: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param &)
        """
        ...
    
    @overload
    def fitEMGPeakModel(self, input_peak: MSChromatogram , output_peak: MSChromatogram ) -> None:
        """
        Cython signature: void fitEMGPeakModel(MSChromatogram & input_peak, MSChromatogram & output_peak)
        """
        ...
    
    @overload
    def fitEMGPeakModel(self, input_peak: MSSpectrum , output_peak: MSSpectrum ) -> None:
        """
        Cython signature: void fitEMGPeakModel(MSSpectrum & input_peak, MSSpectrum & output_peak)
        """
        ...
    
    @overload
    def fitEMGPeakModel(self, input_peak: MSChromatogram , output_peak: MSChromatogram , left_pos: float , right_pos: float ) -> None:
        """
        Cython signature: void fitEMGPeakModel(MSChromatogram & input_peak, MSChromatogram & output_peak, double left_pos, double right_pos)
        """
        ...
    
    @overload
    def fitEMGPeakModel(self, input_peak: MSSpectrum , output_peak: MSSpectrum , left_pos: float , right_pos: float ) -> None:
        """
        Cython signature: void fitEMGPeakModel(MSSpectrum & input_peak, MSSpectrum & output_peak, double left_pos, double right_pos)
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


class FeatureFileOptions:
    """
    Cython implementation of _FeatureFileOptions

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFileOptions.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureFileOptions()
        Options for loading files containing features.
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureFileOptions ) -> None:
        """
        Cython signature: void FeatureFileOptions(FeatureFileOptions &)
        """
        ...
    
    def setMetadataOnly(self, in_0: bool ) -> None:
        """
        Cython signature: void setMetadataOnly(bool)
        Sets whether or not to load only meta data
        """
        ...
    
    def getMetadataOnly(self) -> bool:
        """
        Cython signature: bool getMetadataOnly()
        Returns whether or not to load only meta data
        """
        ...
    
    def setSizeOnly(self, in_0: bool ) -> None:
        """
        Cython signature: void setSizeOnly(bool)
        Sets whether or not to load only feature count
        """
        ...
    
    def getSizeOnly(self) -> bool:
        """
        Cython signature: bool getSizeOnly()
        Returns whether or not to load only meta data
        """
        ...
    
    def setLoadConvexHull(self, in_0: bool ) -> None:
        """
        Cython signature: void setLoadConvexHull(bool)
        Sets whether or not to load convex hull
        """
        ...
    
    def getLoadConvexHull(self) -> bool:
        """
        Cython signature: bool getLoadConvexHull()
        Returns whether or not to load convex hull
        """
        ...
    
    def setLoadSubordinates(self, in_0: bool ) -> None:
        """
        Cython signature: void setLoadSubordinates(bool)
        Sets whether or not load subordinates
        """
        ...
    
    def getLoadSubordinates(self) -> bool:
        """
        Cython signature: bool getLoadSubordinates()
        Returns whether or not to load subordinates
        """
        ...
    
    def setRTRange(self, range_: DRange1 ) -> None:
        """
        Cython signature: void setRTRange(DRange1 & range_)
        Restricts the range of RT values for peaks to load
        """
        ...
    
    def hasRTRange(self) -> bool:
        """
        Cython signature: bool hasRTRange()
        Returns true if an RT range has been set
        """
        ...
    
    def getRTRange(self) -> DRange1:
        """
        Cython signature: DRange1 getRTRange()
        Returns the RT range
        """
        ...
    
    def setMZRange(self, range_: DRange1 ) -> None:
        """
        Cython signature: void setMZRange(DRange1 & range_)
        Restricts the range of MZ values for peaks to load
        """
        ...
    
    def hasMZRange(self) -> bool:
        """
        Cython signature: bool hasMZRange()
        Returns true if an MZ range has been set
        """
        ...
    
    def getMZRange(self) -> DRange1:
        """
        Cython signature: DRange1 getMZRange()
        Returns the MZ range
        """
        ...
    
    def setIntensityRange(self, range_: DRange1 ) -> None:
        """
        Cython signature: void setIntensityRange(DRange1 & range_)
        Restricts the range of intensity values for peaks to load
        """
        ...
    
    def hasIntensityRange(self) -> bool:
        """
        Cython signature: bool hasIntensityRange()
        Returns true if an intensity range has been set
        """
        ...
    
    def getIntensityRange(self) -> DRange1:
        """
        Cython signature: DRange1 getIntensityRange()
        Returns the intensity range
        """
        ... 


class FeatureFinderIdentificationAlgorithm:
    """
    Cython implementation of _FeatureFinderIdentificationAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderIdentificationAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']

    Algorithm class for FeatureFinderIdentification
    
    External IDs (peptides_ext, proteins_ext) may be empty,
    in which case no machine learning or FDR estimation will be performed.
    Optional seeds from e.g. untargeted FeatureFinders can be added with
    seeds.
    Results will be written to features .
    Caution: peptide IDs will be shrunk to best hit, FFid metavalues added
    and potential seed IDs added.
    
    Usage:
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureFinderIdentificationAlgorithm()
        """
        ...
    
    @overload
    def run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap ) -> None:
        """
        Cython signature: void run(PeptideIdentificationList peptides, libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList peptides_ext, libcpp_vector[ProteinIdentification] proteins_ext, FeatureMap & features)
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        """
        ...
    
    @overload
    def run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap ) -> None:
        """
        Cython signature: void run(PeptideIdentificationList peptides, libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList peptides_ext, libcpp_vector[ProteinIdentification] proteins_ext, FeatureMap & features, FeatureMap & seeds)
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        """
        ...
    
    @overload
    def run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap , spectra_file: String ) -> None:
        """
        Cython signature: void run(PeptideIdentificationList peptides, libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList peptides_ext, libcpp_vector[ProteinIdentification] proteins_ext, FeatureMap & features, FeatureMap & seeds, String & spectra_file)
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        :param spectra_file: Path will be stored in features in case the MSExperiment has no proper primaryMSRunPath
        """
        ...
    
    def runOnCandidates(self, features: FeatureMap ) -> None:
        """
        Cython signature: void runOnCandidates(FeatureMap & features)
        Run feature detection on identified features (e.g. loaded from an IdXML file)
        """
        ...
    
    def setMSData(self, in_0: MSExperiment ) -> None:
        """
        Cython signature: void setMSData(const MSExperiment &)
        Sets ms data
        """
        ...
    
    def getMSData(self) -> MSExperiment:
        """
        Cython signature: MSExperiment getMSData()
        Returns ms data as MSExperiment
        """
        ...
    
    def getChromatograms(self) -> MSExperiment:
        """
        Cython signature: MSExperiment getChromatograms()
        Returns chromatogram data as MSExperiment
        """
        ...
    
    def getLibrary(self) -> TargetedExperiment:
        """
        Cython signature: TargetedExperiment getLibrary()
        Returns constructed assay library
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


class IBSpectraFile:
    """
    Cython implementation of _IBSpectraFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IBSpectraFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IBSpectraFile()
        Implements the export of consensusmaps into the IBSpectra format used by isobar to load quantification results
        """
        ...
    
    @overload
    def __init__(self, in_0: IBSpectraFile ) -> None:
        """
        Cython signature: void IBSpectraFile(IBSpectraFile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , cm: ConsensusMap ) -> None:
        """
        Cython signature: void store(const String & filename, ConsensusMap & cm)
        Writes the contents of the ConsensusMap cm into the file named by filename
        
        
        :param filename: The name of the file where the contents of cm should be stored
        :param cm: The ConsensusMap that should be exported to filename
        :raises:
          Exception: InvalidParameter if the ConsensusMap does not hold the result of an isobaric quantification experiment (e.g., itraq)
        """
        ... 


class IsotopeCluster:
    """
    Cython implementation of _IsotopeCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeCluster.html>`_
    """
    
    peaks: ChargedIndexSet
    
    scans: List[int]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsotopeCluster()
        Stores information about an isotopic cluster (i.e. potential peptide charge variants)
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopeCluster ) -> None:
        """
        Cython signature: void IsotopeCluster(IsotopeCluster &)
        """
        ... 


class MapAlignmentAlgorithmKD:
    """
    Cython implementation of _MapAlignmentAlgorithmKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmKD.html>`_

    An efficient reference-free feature map alignment algorithm for unlabeled data
    
    This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
    in a compatibility graph on feature data. This graph is comprised of nodes corresponding
    to features and edges connecting features f and f' iff both are within each other's tolerance
    windows (wrt. RT and m/z difference). CCCs are those CCs that do not contain multiple features
    from the same input map, and whose features all have the same charge state
    
    All CCCs above a user-specified minimum size are considered true sets of corresponding features
    and based on these, LOWESS transformations are computed for each input map such that the average
    deviation from the mean retention time within all CCCs is minimized
    """
    
    @overload
    def __init__(self, num_maps: int , param: Param ) -> None:
        """
        Cython signature: void MapAlignmentAlgorithmKD(size_t num_maps, Param & param)
        """
        ...
    
    @overload
    def __init__(self, in_0: MapAlignmentAlgorithmKD ) -> None:
        """
        Cython signature: void MapAlignmentAlgorithmKD(MapAlignmentAlgorithmKD &)
        """
        ...
    
    def addRTFitData(self, kd_data: KDTreeFeatureMaps ) -> None:
        """
        Cython signature: void addRTFitData(KDTreeFeatureMaps & kd_data)
        Compute data points needed for RT transformation in the current `kd_data`, add to `fit_data_`
        """
        ...
    
    def fitLOWESS(self) -> None:
        """
        Cython signature: void fitLOWESS()
        Fit LOWESS to fit_data_, store final models in `transformations_`
        """
        ...
    
    def transform(self, kd_data: KDTreeFeatureMaps ) -> None:
        """
        Cython signature: void transform(KDTreeFeatureMaps & kd_data)
        Transform RTs for `kd_data`
        """
        ... 


class MetaInfoDescription:
    """
    Cython implementation of _MetaInfoDescription

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoDescription.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaInfoDescription()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaInfoDescription ) -> None:
        """
        Cython signature: void MetaInfoDescription(MetaInfoDescription &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the peak annotations
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the peak annotations
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[shared_ptr[DataProcessing]] getDataProcessing()
        Returns a reference to the description of the applied processing
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[shared_ptr[DataProcessing]])
        Sets the description of the applied processing
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
    
    def __richcmp__(self, other: MetaInfoDescription, op: int) -> Any:
        ... 


class MetaboTargetedAssay:
    """
    Cython implementation of _MetaboTargetedAssay

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedAssay.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboTargetedAssay()
        This class provides methods for the extraction of targeted assays for metabolomics
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboTargetedAssay ) -> None:
        """
        Cython signature: void MetaboTargetedAssay(MetaboTargetedAssay &)
        """
        ...
    
    def extractMetaboTargetedAssay(self, spectra: MSExperiment , feature_ms2_index: FeatureMapping_FeatureToMs2Indices , precursor_rt_tol: float , precursor_mz_distance: float , cosine_sim_threshold: float , transition_threshold: float , min_fragment_mz: float , max_fragment_mz: float , method_consensus_spectrum: bool , exclude_ms2_precursor: bool , file_counter: int ) -> List[MetaboTargetedAssay]:
        """
        Cython signature: libcpp_vector[MetaboTargetedAssay] extractMetaboTargetedAssay(MSExperiment & spectra, FeatureMapping_FeatureToMs2Indices & feature_ms2_index, double & precursor_rt_tol, double & precursor_mz_distance, double & cosine_sim_threshold, double & transition_threshold, double & min_fragment_mz, double & max_fragment_mz, bool & method_consensus_spectrum, bool & exclude_ms2_precursor, unsigned int & file_counter)
        Extract a vector of MetaboTargetedAssays without using fragment annotation
        
        
        :param spectra: Input of MSExperiment with spectra information
        :param feature_ms2_spectra_map: FeatureMapping class with associated MS2 spectra
        :param precursor_rt_tol: Retention time tolerance of the precursor
        :param precursor_mz_distance: Max m/z distance of the precursor entries of two spectra to be merged
        :param cosine_sim_threshold: Cosine similarty threshold for the usage of SpectraMerger
        :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
        :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
        :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
        :param method_consensus_spectrum: Boolean to use consensus spectrum method
        :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
        :return: Vector of MetaboTargetedAssay
        """
        ...
    
    def extractMetaboTargetedAssayFragmentAnnotation(self, v_cmp_spec: List[MetaboTargetedAssay_CompoundTargetDecoyPair] , transition_threshold: float , min_fragment_mz: float , max_fragment_mz: float , use_exact_mass: bool , exclude_ms2_precursor: bool ) -> List[MetaboTargetedAssay]:
        """
        Cython signature: libcpp_vector[MetaboTargetedAssay] extractMetaboTargetedAssayFragmentAnnotation(libcpp_vector[MetaboTargetedAssay_CompoundTargetDecoyPair] & v_cmp_spec, double & transition_threshold, double & min_fragment_mz, double & max_fragment_mz, bool & use_exact_mass, bool & exclude_ms2_precursor)
        Extract a vector of MetaboTargetedAssays using fragment
        
        
        :param v_cmp_spec: Vector of CompoundInfo with associated fragment annotated MSspectrum
        :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
        :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
        :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
        :param use_exact_mass: Boolean if exact mass should be used as peak mass for annotated fragments
        :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
        :param file_counter: Count if multiple files are used.
        :return: Vector of MetaboTargetedAssay
        """
        ...
    
    def pairCompoundWithAnnotatedTDSpectraPairs(self, v_cmpinfo: List[SiriusMSFile_CompoundInfo] , annotated_spectra: List[SiriusFragmentAnnotation_SiriusTargetDecoySpectra] ) -> List[MetaboTargetedAssay_CompoundTargetDecoyPair]:
        """
        Cython signature: libcpp_vector[MetaboTargetedAssay_CompoundTargetDecoyPair] pairCompoundWithAnnotatedTDSpectraPairs(libcpp_vector[SiriusMSFile_CompoundInfo] & v_cmpinfo, libcpp_vector[SiriusFragmentAnnotation_SiriusTargetDecoySpectra] & annotated_spectra)
        Pair compound information (SiriusMSFile) with the annotated target and decoy spectrum from SIRIUS/Passatutto based on the m_id (unique identifier composed of description_filepath_native_id_k introduced in the SiriusMSConverter)
        
        
        :param v_cmpinfo: Vector of SiriusMSFile::CompoundInfo
        :param annotated_spectra: Vector of SiriusTargetDecoySpectra
        :return: Vector of MetaboTargetedAssay::CompoundTargetDecoyPair
        """
        ... 


class MetaboTargetedAssay_CompoundTargetDecoyPair:
    """
    Cython implementation of _MetaboTargetedAssay_CompoundTargetDecoyPair

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedAssay_CompoundTargetDecoyPair.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboTargetedAssay_CompoundTargetDecoyPair()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboTargetedAssay_CompoundTargetDecoyPair ) -> None:
        """
        Cython signature: void MetaboTargetedAssay_CompoundTargetDecoyPair(MetaboTargetedAssay_CompoundTargetDecoyPair &)
        """
        ... 


class MorpheusScore:
    """
    Cython implementation of _MorpheusScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorpheusScore.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MorpheusScore()
        """
        ...
    
    @overload
    def __init__(self, in_0: MorpheusScore ) -> None:
        """
        Cython signature: void MorpheusScore(MorpheusScore &)
        """
        ...
    
    def compute(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , theo_spectrum: MSSpectrum ) -> MorpheusScore_Result:
        """
        Cython signature: MorpheusScore_Result compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const MSSpectrum & exp_spectrum, const MSSpectrum & theo_spectrum)
        Returns Morpheus Score
        """
        ... 


class MorpheusScore_Result:
    """
    Cython implementation of _MorpheusScore_Result

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorpheusScore_Result.html>`_
    """
    
    matches: int
    
    n_peaks: int
    
    score: float
    
    MIC: float
    
    TIC: float
    
    err: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MorpheusScore_Result()
        """
        ...
    
    @overload
    def __init__(self, in_0: MorpheusScore_Result ) -> None:
        """
        Cython signature: void MorpheusScore_Result(MorpheusScore_Result &)
        """
        ... 


class MzTabMFile:
    """
    Cython implementation of _MzTabMFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabMFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzTabMFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzTabMFile ) -> None:
        """
        Cython signature: void MzTabMFile(MzTabMFile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , mztab_m: MzTabM ) -> None:
        """
        Cython signature: void store(String filename, MzTabM & mztab_m)
        Store MzTabM file
        """
        ... 


class ParamCTDFile:
    """
    Cython implementation of _ParamCTDFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamCTDFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ParamCTDFile()
        """
        ...
    
    def store(self, filename: Union[bytes, str] , param: Param , tool_info: ToolInfo ) -> None:
        """
        Cython signature: void store(libcpp_utf8_string filename, Param param, ToolInfo tool_info)
        """
        ... 


class PercolatorInfile:
    """
    Cython implementation of _PercolatorInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorInfile.html>`_

    Class for storing Percolator tab-delimited input files
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PercolatorInfile()
        """
        ...
    
    @overload
    def __init__(self, in_0: PercolatorInfile ) -> None:
        """
        Cython signature: void PercolatorInfile(PercolatorInfile &)
        """
        ...
    
    store: __static_PercolatorInfile_store 


class SiriusExportAlgorithm:
    """
    Cython implementation of _SiriusExportAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusExportAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SiriusExportAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: SiriusExportAlgorithm ) -> None:
        """
        Cython signature: void SiriusExportAlgorithm(SiriusExportAlgorithm &)
        """
        ...
    
    def isFeatureOnly(self) -> bool:
        """
        Cython signature: bool isFeatureOnly()
        """
        ...
    
    def getFilterByNumMassTraces(self) -> int:
        """
        Cython signature: unsigned int getFilterByNumMassTraces()
        """
        ...
    
    def getPrecursorMzTolerance(self) -> float:
        """
        Cython signature: double getPrecursorMzTolerance()
        """
        ...
    
    def getPrecursorRtTolerance(self) -> float:
        """
        Cython signature: double getPrecursorRtTolerance()
        """
        ...
    
    def precursorMzToleranceUnitIsPPM(self) -> bool:
        """
        Cython signature: bool precursorMzToleranceUnitIsPPM()
        """
        ...
    
    def isNoMasstraceInfoIsotopePattern(self) -> bool:
        """
        Cython signature: bool isNoMasstraceInfoIsotopePattern()
        """
        ...
    
    def getIsotopePatternIterations(self) -> int:
        """
        Cython signature: int getIsotopePatternIterations()
        """
        ...
    
    def preprocessing(self, featureXML_path: Union[bytes, str, String] , spectra: MSExperiment , feature_mapping_info: FeatureMapping_FeatureMappingInfo , feature_ms2_indices: FeatureMapping_FeatureToMs2Indices ) -> None:
        """
        Cython signature: void preprocessing(const String & featureXML_path, MSExperiment & spectra, FeatureMapping_FeatureMappingInfo & feature_mapping_info, FeatureMapping_FeatureToMs2Indices & feature_ms2_indices)
        Preprocessing needed for SIRIUS
        
        Filter number of masstraces and perform feature mapping
        
        :param featureXML_path: Path to featureXML
        :param spectra: Input of MSExperiment with spectra information
        :param feature_mapping_info: Emtpy - stores FeatureMaps and KDTreeMaps internally
        :param feature_ms2_indices: Empty FeatureToMs2Indices
        """
        ...
    
    def logFeatureSpectraNumber(self, featureXML_path: Union[bytes, str, String] , feature_ms2_indices: FeatureMapping_FeatureToMs2Indices , spectra: MSExperiment ) -> None:
        """
        Cython signature: void logFeatureSpectraNumber(const String & featureXML_path, FeatureMapping_FeatureToMs2Indices & feature_ms2_indices, MSExperiment & spectra)
        Logs number of features and spectra used
        
        Prints the number of features and spectra used (OPENMS_LOG_INFO)
        
        :param featureXML_path: Path to featureXML
        :param feature_ms2_indices: FeatureToMs2Indices with feature mapping
        :param spectra: Input of MSExperiment with spectra information
        """
        ...
    
    def run(self, mzML_files: List[bytes] , featureXML_files: List[bytes] , out_ms: Union[bytes, str, String] , out_compoundinfo: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void run(const StringList & mzML_files, const StringList & featureXML_files, const String & out_ms, const String & out_compoundinfo)
        Runs SiriusExport with mzML and featureXML (optional) files as input.
        
        Generates a SIRIUS .ms file and compound info table (optional).
        
        :param mzML_files: List with paths to mzML files
        :param featureXML_files: List with paths to featureXML files
        :param out_ms: Output file name for SIRIUS .ms file
        :param out_compoundinfo: Output file name for tsv file with compound info
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


class ToolInfo:
    """
    Cython implementation of _ToolInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ToolInfo.html>`_
    """
    
    def __init__(self, in_0: ToolInfo ) -> None:
        """
        Cython signature: void ToolInfo(ToolInfo)
        """
        ... 

