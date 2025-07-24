from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_FeatureMapping_assignMS2IndexToFeature(spectra: MSExperiment , fm_info: FeatureMapping_FeatureMappingInfo , precursor_mz_tolerance: float , precursor_rt_tolerance: float , ppm: bool ) -> FeatureMapping_FeatureToMs2Indices:
    """
    Cython signature: FeatureMapping_FeatureToMs2Indices assignMS2IndexToFeature(MSExperiment & spectra, FeatureMapping_FeatureMappingInfo & fm_info, double precursor_mz_tolerance, double precursor_rt_tolerance, bool ppm)
    """
    ...


class Adduct:
    """
    Cython implementation of _Adduct

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Adduct.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Adduct()
        """
        ...
    
    @overload
    def __init__(self, in_0: Adduct ) -> None:
        """
        Cython signature: void Adduct(Adduct &)
        """
        ...
    
    @overload
    def __init__(self, charge: int ) -> None:
        """
        Cython signature: void Adduct(int charge)
        """
        ...
    
    @overload
    def __init__(self, charge: int , amount: int , singleMass: float , formula: Union[bytes, str, String] , log_prob: float , rt_shift: float , label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void Adduct(int charge, int amount, double singleMass, String formula, double log_prob, double rt_shift, String label)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        """
        ...
    
    def setCharge(self, charge: int ) -> None:
        """
        Cython signature: void setCharge(int charge)
        """
        ...
    
    def getAmount(self) -> int:
        """
        Cython signature: int getAmount()
        """
        ...
    
    def setAmount(self, amount: int ) -> None:
        """
        Cython signature: void setAmount(int amount)
        """
        ...
    
    def getSingleMass(self) -> float:
        """
        Cython signature: double getSingleMass()
        """
        ...
    
    def setSingleMass(self, singleMass: float ) -> None:
        """
        Cython signature: void setSingleMass(double singleMass)
        """
        ...
    
    def getLogProb(self) -> float:
        """
        Cython signature: double getLogProb()
        """
        ...
    
    def setLogProb(self, log_prob: float ) -> None:
        """
        Cython signature: void setLogProb(double log_prob)
        """
        ...
    
    def getFormula(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFormula()
        """
        ...
    
    def setFormula(self, formula: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFormula(String formula)
        """
        ...
    
    def getRTShift(self) -> float:
        """
        Cython signature: double getRTShift()
        """
        ...
    
    def getLabel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getLabel()
        """
        ... 


class AverageLinkage:
    """
    Cython implementation of _AverageLinkage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AverageLinkage.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AverageLinkage()
        """
        ...
    
    @overload
    def __init__(self, in_0: AverageLinkage ) -> None:
        """
        Cython signature: void AverageLinkage(AverageLinkage &)
        """
        ... 


class ChromatogramExtractorAlgorithm:
    """
    Cython implementation of _ChromatogramExtractorAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramExtractorAlgorithm.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramExtractorAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramExtractorAlgorithm ) -> None:
        """
        Cython signature: void ChromatogramExtractorAlgorithm(ChromatogramExtractorAlgorithm &)
        """
        ...
    
    def extractChromatograms(self, input: SpectrumAccessOpenMS , output: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , mz_extraction_window: float , ppm: bool , im_extraction_window: float , filter: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void extractChromatograms(shared_ptr[SpectrumAccessOpenMS] input, libcpp_vector[shared_ptr[OSChromatogram]] & output, libcpp_vector[ExtractionCoordinates] extraction_coordinates, double mz_extraction_window, bool ppm, double im_extraction_window, String filter)
          Extract chromatograms at the m/z and RT defined by the ExtractionCoordinates
        
        
        :param input: Input spectral map
        :param output: Output chromatograms (XICs)
        :param extraction_coordinates: Extracts around these coordinates (from
         rt_start to rt_end in seconds - extracts the whole chromatogram if
         rt_end - rt_start < 0).
        :param mz_extraction_window: Extracts a window of this size in m/z
          dimension in Th or ppm (e.g. a window of 50 ppm means an extraction of
          25 ppm on either side)
        :param ppm: Whether mz_extraction_window is in ppm or in Th
        :param filter: Which function to apply in m/z space (currently "tophat" only)
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


class ConsensusFeature:
    """
    Cython implementation of _ConsensusFeature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusFeature.html>`_
      -- Inherits from ['UniqueIdInterface', 'BaseFeature']

    A consensus feature spanning multiple LC-MS/MS experiments.
    
    A ConsensusFeature represents analytes that have been
    quantified across multiple LC-MS/MS experiments. Each analyte in a
    ConsensusFeature is linked to its original LC-MS/MS run through a
    unique identifier.
    
    Get access to the underlying features through getFeatureList()
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ConsensusFeature()
        """
        ...
    
    @overload
    def __init__(self, in_0: ConsensusFeature ) -> None:
        """
        Cython signature: void ConsensusFeature(ConsensusFeature &)
        """
        ...
    
    @overload
    def __init__(self, in_0: int , in_1: Peak2D , in_2: int ) -> None:
        """
        Cython signature: void ConsensusFeature(uint64_t, Peak2D, uint64_t)
        """
        ...
    
    @overload
    def __init__(self, in_0: int , in_1: BaseFeature ) -> None:
        """
        Cython signature: void ConsensusFeature(uint64_t, BaseFeature)
        """
        ...
    
    @overload
    def __init__(self, in_0: int , in_1: ConsensusFeature ) -> None:
        """
        Cython signature: void ConsensusFeature(uint64_t, ConsensusFeature)
        """
        ...
    
    def computeConsensus(self) -> None:
        """
        Cython signature: void computeConsensus()
        Computes and updates the consensus position, intensity, and charge
        """
        ...
    
    def computeMonoisotopicConsensus(self) -> None:
        """
        Cython signature: void computeMonoisotopicConsensus()
        Computes and updates the consensus position, intensity, and charge
        """
        ...
    
    def computeDechargeConsensus(self, in_0: FeatureMap , in_1: bool ) -> None:
        """
        Cython signature: void computeDechargeConsensus(FeatureMap, bool)
        Computes the uncharged parent RT & mass, assuming the handles are charge variants
        """
        ...
    
    @overload
    def insert(self, map_idx: int , in_1: Peak2D , element_idx: int ) -> None:
        """
        Cython signature: void insert(uint64_t map_idx, Peak2D, uint64_t element_idx)
        """
        ...
    
    @overload
    def insert(self, map_idx: int , in_1: BaseFeature ) -> None:
        """
        Cython signature: void insert(uint64_t map_idx, BaseFeature)
        """
        ...
    
    @overload
    def insert(self, map_idx: int , in_1: ConsensusFeature ) -> None:
        """
        Cython signature: void insert(uint64_t map_idx, ConsensusFeature)
        """
        ...
    
    def getFeatureList(self) -> List[FeatureHandle]:
        """
        Cython signature: libcpp_vector[FeatureHandle] getFeatureList()
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def addRatio(self, r: Ratio ) -> None:
        """
        Cython signature: void addRatio(Ratio r)
        Connects a ratio to the ConsensusFeature.
        """
        ...
    
    def setRatios(self, rs: List[Ratio] ) -> None:
        """
        Cython signature: void setRatios(libcpp_vector[Ratio] rs)
        Connects the ratios to the ConsensusFeature.
        """
        ...
    
    def getRatios(self) -> List[Ratio]:
        """
        Cython signature: libcpp_vector[Ratio] getRatios()
        Get the ratio vector.
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
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
    
    def getQuality(self) -> float:
        """
        Cython signature: float getQuality()
        Returns the overall quality
        """
        ...
    
    def setQuality(self, q: float ) -> None:
        """
        Cython signature: void setQuality(float q)
        Sets the overall quality
        """
        ...
    
    def getWidth(self) -> float:
        """
        Cython signature: float getWidth()
        Returns the features width (full width at half max, FWHM)
        """
        ...
    
    def setWidth(self, q: float ) -> None:
        """
        Cython signature: void setWidth(float q)
        Sets the width of the feature (FWHM)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        Returns the charge state
        """
        ...
    
    def setCharge(self, q: int ) -> None:
        """
        Cython signature: void setCharge(int q)
        Sets the charge state
        """
        ...
    
    def getAnnotationState(self) -> int:
        """
        Cython signature: AnnotationState getAnnotationState()
        State of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only
        """
        ...
    
    def getPeptideIdentifications(self) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList getPeptideIdentifications()
        Returns the PeptideIdentification vector
        """
        ...
    
    def setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setPeptideIdentifications(PeptideIdentificationList & peptides)
        Sets the PeptideIdentification vector
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
    
    def __richcmp__(self, other: ConsensusFeature, op: int) -> Any:
        ... 


class ConsensusXMLFile:
    """
    Cython implementation of _ConsensusXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusXMLFile()
        """
        ...
    
    def load(self, in_0: Union[bytes, str, String] , in_1: ConsensusMap ) -> None:
        """
        Cython signature: void load(String, ConsensusMap &)
        Loads a consensus map from file and calls updateRanges
        """
        ...
    
    def store(self, in_0: Union[bytes, str, String] , in_1: ConsensusMap ) -> None:
        """
        Cython signature: void store(String, ConsensusMap &)
        Stores a consensus map to file
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
        Mutable access to the options for loading/storing
        """
        ... 


class ExtractionCoordinates:
    """
    Cython implementation of _ExtractionCoordinates

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExtractionCoordinates.html>`_
    """
    
    mz: float
    
    mz_precursor: float
    
    rt_start: float
    
    rt_end: float
    
    ion_mobility: float
    
    id: bytes
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ExtractionCoordinates()
        """
        ...
    
    @overload
    def __init__(self, in_0: ExtractionCoordinates ) -> None:
        """
        Cython signature: void ExtractionCoordinates(ExtractionCoordinates)
        """
        ... 


class FeatureMapping:
    """
    Cython implementation of _FeatureMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureMapping()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureMapping ) -> None:
        """
        Cython signature: void FeatureMapping(FeatureMapping &)
        """
        ...
    
    assignMS2IndexToFeature: __static_FeatureMapping_assignMS2IndexToFeature 


class FeatureMapping_FeatureMappingInfo:
    """
    Cython implementation of _FeatureMapping_FeatureMappingInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping_FeatureMappingInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureMapping_FeatureMappingInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureMapping_FeatureMappingInfo ) -> None:
        """
        Cython signature: void FeatureMapping_FeatureMappingInfo(FeatureMapping_FeatureMappingInfo &)
        """
        ... 


class FeatureMapping_FeatureToMs2Indices:
    """
    Cython implementation of _FeatureMapping_FeatureToMs2Indices

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping_FeatureToMs2Indices.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureMapping_FeatureToMs2Indices()
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureMapping_FeatureToMs2Indices ) -> None:
        """
        Cython signature: void FeatureMapping_FeatureToMs2Indices(FeatureMapping_FeatureToMs2Indices &)
        """
        ... 


class GNPSQuantificationFile:
    """
    Cython implementation of _GNPSQuantificationFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSQuantificationFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GNPSQuantificationFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: GNPSQuantificationFile ) -> None:
        """
        Cython signature: void GNPSQuantificationFile(GNPSQuantificationFile &)
        """
        ...
    
    def store(self, consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const ConsensusMap & consensus_map, const String & output_file)
        Write feature quantification table (txt file) from a ConsensusMap. Required for GNPS FBMN.
        
        The table contains map information on the featureXML files from which the ConsensusMap was generated as well as
        a row for every consensus feature with information on rt, mz, intensity, width and quality. The same information is
        added for each original feature in the consensus feature.
        
        :param consensus_map: Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        :param output_file: Output file path for the feature quantification table.
        """
        ... 


class IDMapper:
    """
    Cython implementation of _IDMapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDMapper.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IDMapper()
        Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications
        """
        ...
    
    @overload
    def __init__(self, in_0: IDMapper ) -> None:
        """
        Cython signature: void IDMapper(IDMapper &)
        """
        ...
    
    @overload
    def annotate(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None:
        """
        Cython signature: void annotate(AnnotatedMSRun & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool clear_ids, bool mapMS1)
        Mapping method for peak maps\n
        
        The identifications stored in a PeptideIdentification instance can be added to the
        corresponding spectrum
        Note that a PeptideIdentication is added to ALL spectra which are within the allowed RT and MZ boundaries
        
        
        :param map: AnnotatedMSRun to receive the identifications
        :param peptide_ids: PeptideIdentification for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param clear_ids: Reset peptide and protein identifications of each scan before annotating
        :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
        :raises:
          Exception: MissingInformation is thrown if entries of 'peptide_ids' do not contain 'MZ' and 'RT' information
        """
        ...
    
    @overload
    def annotate(self, map_: AnnotatedMSRun , fmap: FeatureMap , clear_ids: bool , mapMS1: bool ) -> None:
        """
        Cython signature: void annotate(AnnotatedMSRun & map_, FeatureMap & fmap, bool clear_ids, bool mapMS1)
        Mapping method for peak maps\n
        
        Add peptide identifications stored in a feature map to their
        corresponding spectrum
        This function converts the feature map to a vector of peptide identifications (all peptide IDs from each feature are taken)
        and calls the respective annotate() function
        RT and m/z are taken from the peptides, or (if missing) from the feature itself
        
        
        :param map: AnnotatedMSRun to receive the identifications
        :param fmap: FeatureMap with PeptideIdentifications for the MSExperiment
        :param clear_ids: Reset peptide and protein identifications of each scan before annotating
        :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
        """
        ...
    
    @overload
    def annotate(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None:
        """
        Cython signature: void annotate(FeatureMap & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool use_centroid_rt, bool use_centroid_mz, MSExperiment & spectra)
        Mapping method for peak maps\n
        
        If all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if 'use_centroid_rt' or 'use_centroid_mz' are true\n
        
        In any case, tolerance in RT and m/z dimension is applied according to the global parameters 'rt_tolerance' and 'mz_tolerance'. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value\n
        
        If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them
        
        
        :param map: MSExperiment to receive the identifications
        :param ids: PeptideIdentification for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param use_centroid_rt: Whether to use the RT value of feature centroids even if convex hulls are present
        :param use_centroid_mz: Whether to use the m/z value of feature centroids even if convex hulls are present
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        :raises:
          Exception: MissingInformation is thrown if entries of 'ids' do not contain 'MZ' and 'RT' information
        """
        ...
    
    @overload
    def annotate(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None:
        """
        Cython signature: void annotate(ConsensusMap & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool measure_from_subelements, bool annotate_ids_with_subelements, MSExperiment & spectra)
        Mapping method for peak maps\n
        
        If all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if 'use_centroid_rt' or 'use_centroid_mz' are true\n
        
        In any case, tolerance in RT and m/z dimension is applied according to the global parameters 'rt_tolerance' and 'mz_tolerance'. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value\n
        
        If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them
        
        
        :param map: MSExperiment to receive the identifications
        :param ids: PeptideIdentification for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param measure_from_subelements: Boolean operator set to true if distance estimate from FeatureHandles instead of Centroid
        :param annotate_ids_with_subelements: Boolean operator set to true if store map index of FeatureHandle in peptide identification
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        :raises:
          Exception: MissingInformation is thrown if entries of 'ids' do not contain 'MZ' and 'RT' information
        """
        ...
    
    @overload
    def annotate(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None:
        """
        Cython signature: void annotate(AnnotatedMSRun & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool clear_ids, bool mapMS1)
        Mapping method using PeptideIdentificationList\n
        
        The identifications stored in a PeptideIdentificationList instance can be added to the
        corresponding spectrum
        
        
        :param map: AnnotatedMSRun to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param clear_ids: Reset peptide and protein identifications of each scan before annotating
        :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
        """
        ...
    
    @overload
    def annotate(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None:
        """
        Cython signature: void annotate(FeatureMap & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool use_centroid_rt, bool use_centroid_mz, MSExperiment & spectra)
        Mapping method using PeptideIdentificationList\n
        
        :param map: FeatureMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param use_centroid_rt: Whether to use the RT value of feature centroids even if convex hulls are present
        :param use_centroid_mz: Whether to use the m/z value of feature centroids even if convex hulls are present
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        """
        ...
    
    @overload
    def annotate(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None:
        """
        Cython signature: void annotate(ConsensusMap & map_, PeptideIdentificationList & ids, libcpp_vector[ProteinIdentification] & protein_ids, bool measure_from_subelements, bool annotate_ids_with_subelements, MSExperiment & spectra)
        Mapping method using PeptideIdentificationList\n
        
        :param map: ConsensusMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param measure_from_subelements: Boolean operator set to true if distance estimate from FeatureHandles instead of Centroid
        :param annotate_ids_with_subelements: Boolean operator set to true if store map index of FeatureHandle in peptide identification
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        """
        ...
    
    @overload
    def mapPrecursorsToIdentifications(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState:
        """
        Cython signature: IDMapper_PeptideIdentificationListState mapPrecursorsToIdentifications(MSExperiment spectra, PeptideIdentificationList & ids, double mz_tol, double rt_tol)
        Mapping of peptide identifications to spectra\n
        This helper function partitions all spectra into those that had:
        - no precursor (e.g. MS1 spectra),
        - at least one identified precursor,
        - or only unidentified precursor
        
        
        :param spectra: The mass spectra
        :param ids: The peptide identifications
        :param mz_tol: Tolerance used to map to precursor m/z
        :param rt_tol: Tolerance used to map to spectrum retention time
        :return: A struct of vectors holding spectra indices of the partitioning
        """
        ...
    
    @overload
    def mapPrecursorsToIdentifications(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState:
        """
        Cython signature: IDMapper_PeptideIdentificationListState mapPrecursorsToIdentifications(MSExperiment spectra, PeptideIdentificationList & ids, double mz_tol, double rt_tol)
        Mapping of PeptideIdentificationList to spectra\n
        This helper function partitions all spectra into those that had:
        - no precursor (e.g. MS1 spectra),
        - at least one identified precursor,
        - or only unidentified precursor
        
        
        :param spectra: The mass spectra
        :param ids: The PeptideIdentificationList
        :param mz_tol: Tolerance used to map to precursor m/z
        :param rt_tol: Tolerance used to map to spectrum retention time
        :return: A struct of vectors holding spectra indices of the partitioning
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


class IDMapper_PeptideIdentificationListState:
    """
    Cython implementation of _IDMapper_PeptideIdentificationListState

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDMapper_PeptideIdentificationListState.html>`_
    """
    
    no_precursors: List[int]
    
    identified: List[int]
    
    unidentified: List[int]
    
    def __init__(self) -> None:
        """
        Cython signature: void IDMapper_PeptideIdentificationListState()
        """
        ... 


class IMSElement:
    """
    Cython implementation of _IMSElement

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSElement_1_1IMSElement.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IMSElement()
        Represents a chemical atom with name and isotope distribution
        """
        ...
    
    @overload
    def __init__(self, in_0: IMSElement ) -> None:
        """
        Cython signature: void IMSElement(IMSElement &)
        """
        ...
    
    @overload
    def __init__(self, name: bytes , isotopes: IMSIsotopeDistribution ) -> None:
        """
        Cython signature: void IMSElement(libcpp_string & name, IMSIsotopeDistribution & isotopes)
        """
        ...
    
    @overload
    def __init__(self, name: bytes , mass: float ) -> None:
        """
        Cython signature: void IMSElement(libcpp_string & name, double mass)
        """
        ...
    
    @overload
    def __init__(self, name: bytes , nominal_mass: int ) -> None:
        """
        Cython signature: void IMSElement(libcpp_string & name, unsigned int nominal_mass)
        """
        ...
    
    def getName(self) -> bytes:
        """
        Cython signature: libcpp_string getName()
        Gets element's name
        """
        ...
    
    def setName(self, name: bytes ) -> None:
        """
        Cython signature: void setName(libcpp_string & name)
        Sets element's name
        """
        ...
    
    def getSequence(self) -> bytes:
        """
        Cython signature: libcpp_string getSequence()
        Gets element's sequence
        """
        ...
    
    def setSequence(self, sequence: bytes ) -> None:
        """
        Cython signature: void setSequence(libcpp_string & sequence)
        Sets element's sequence
        """
        ...
    
    def getNominalMass(self) -> int:
        """
        Cython signature: unsigned int getNominalMass()
        Gets element's nominal mass
        """
        ...
    
    def getMass(self, index: int ) -> float:
        """
        Cython signature: double getMass(int index)
        Gets mass of element's isotope 'index'
        """
        ...
    
    def getAverageMass(self) -> float:
        """
        Cython signature: double getAverageMass()
        Gets element's average mass
        """
        ...
    
    def getIonMass(self, electrons_number: int ) -> float:
        """
        Cython signature: double getIonMass(int electrons_number)
        Gets ion mass of element. By default ion lacks 1 electron, but this can be changed by setting other 'electrons_number'
        """
        ...
    
    def getIsotopeDistribution(self) -> IMSIsotopeDistribution:
        """
        Cython signature: IMSIsotopeDistribution getIsotopeDistribution()
        Gets element's isotope distribution
        """
        ...
    
    def setIsotopeDistribution(self, isotopes: IMSIsotopeDistribution ) -> None:
        """
        Cython signature: void setIsotopeDistribution(IMSIsotopeDistribution & isotopes)
        Sets element's isotope distribution
        """
        ...
    
    def __richcmp__(self, other: IMSElement, op: int) -> Any:
        ... 


class Instrument:
    """
    Cython implementation of _Instrument

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Instrument.html>`_
      -- Inherits from ['MetaInfoInterface']

    Description of a MS instrument
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Instrument()
        Description of a MS instrument
        """
        ...
    
    @overload
    def __init__(self, in_0: Instrument ) -> None:
        """
        Cython signature: void Instrument(Instrument &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the instrument
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the instrument
        """
        ...
    
    def getVendor(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVendor()
        Returns the instrument vendor
        """
        ...
    
    def setVendor(self, vendor: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setVendor(String vendor)
        Sets the instrument vendor
        """
        ...
    
    def getModel(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getModel()
        Returns the instrument model
        """
        ...
    
    def setModel(self, model: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setModel(String model)
        Sets the instrument model
        """
        ...
    
    def getCustomizations(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCustomizations()
        Returns a description of customizations
        """
        ...
    
    def setCustomizations(self, customizations: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCustomizations(String customizations)
        Sets the a description of customizations
        """
        ...
    
    def getIonSources(self) -> List[IonSource]:
        """
        Cython signature: libcpp_vector[IonSource] getIonSources()
        Returns the ion source list
        """
        ...
    
    def setIonSources(self, ion_sources: List[IonSource] ) -> None:
        """
        Cython signature: void setIonSources(libcpp_vector[IonSource] ion_sources)
        Sets the ion source list
        """
        ...
    
    def getMassAnalyzers(self) -> List[MassAnalyzer]:
        """
        Cython signature: libcpp_vector[MassAnalyzer] getMassAnalyzers()
        Returns the mass analyzer list
        """
        ...
    
    def setMassAnalyzers(self, mass_analyzers: List[MassAnalyzer] ) -> None:
        """
        Cython signature: void setMassAnalyzers(libcpp_vector[MassAnalyzer] mass_analyzers)
        Sets the mass analyzer list
        """
        ...
    
    def getIonDetectors(self) -> List[IonDetector]:
        """
        Cython signature: libcpp_vector[IonDetector] getIonDetectors()
        Returns the ion detector list
        """
        ...
    
    def setIonDetectors(self, ion_detectors: List[IonDetector] ) -> None:
        """
        Cython signature: void setIonDetectors(libcpp_vector[IonDetector] ion_detectors)
        Sets the ion detector list
        """
        ...
    
    def getSoftware(self) -> Software:
        """
        Cython signature: Software getSoftware()
        Returns the instrument software
        """
        ...
    
    def setSoftware(self, software: Software ) -> None:
        """
        Cython signature: void setSoftware(Software software)
        Sets the instrument software
        """
        ...
    
    def getIonOptics(self) -> int:
        """
        Cython signature: IonOpticsType getIonOptics()
        Returns the ion optics type
        """
        ...
    
    def setIonOptics(self, ion_optics: int ) -> None:
        """
        Cython signature: void setIonOptics(IonOpticsType ion_optics)
        Sets the ion optics type
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
    
    def __richcmp__(self, other: Instrument, op: int) -> Any:
        ... 


class IonSource:
    """
    Cython implementation of _IonSource

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonSource.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IonSource()
        Description of an ion source (part of a MS Instrument)
        """
        ...
    
    @overload
    def __init__(self, in_0: IonSource ) -> None:
        """
        Cython signature: void IonSource(IonSource &)
        """
        ...
    
    def getPolarity(self) -> int:
        """
        Cython signature: Polarity getPolarity()
        Returns the ionization mode
        """
        ...
    
    def setPolarity(self, polarity: int ) -> None:
        """
        Cython signature: void setPolarity(Polarity polarity)
        Sets the ionization mode
        """
        ...
    
    def getInletType(self) -> int:
        """
        Cython signature: InletType getInletType()
        Returns the inlet type
        """
        ...
    
    def setInletType(self, inlet_type: int ) -> None:
        """
        Cython signature: void setInletType(InletType inlet_type)
        Sets the inlet type
        """
        ...
    
    def getIonizationMethod(self) -> int:
        """
        Cython signature: IonizationMethod getIonizationMethod()
        Returns the ionization method
        """
        ...
    
    def setIonizationMethod(self, ionization_type: int ) -> None:
        """
        Cython signature: void setIonizationMethod(IonizationMethod ionization_type)
        Sets the ionization method
        """
        ...
    
    def getOrder(self) -> int:
        """
        Cython signature: int getOrder()
        Returns the position of this part in the whole Instrument
        
        Order can be ignored, as long the instrument has this default setup:
          - one ion source
          - one or many mass analyzers
          - one ion detector
        
        For more complex instruments, the order should be defined.
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
    
    def __richcmp__(self, other: IonSource, op: int) -> Any:
        ...
    InletType : __InletType
    IonizationMethod : __IonizationMethod
    Polarity : __Polarity 


class OpenPepXLLFAlgorithm:
    """
    Cython implementation of _OpenPepXLLFAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenPepXLLFAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']

    Search for cross-linked peptide pairs in tandem MS spectra
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenPepXLLFAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenPepXLLFAlgorithm ) -> None:
        """
        Cython signature: void OpenPepXLLFAlgorithm(OpenPepXLLFAlgorithm &)
        """
        ...
    
    def run(self, unprocessed_spectra: MSExperiment , fasta_db: List[FASTAEntry] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment ) -> int:
        """
        Cython signature: OpenPepXLLFAlgorithm_ExitCodes run(MSExperiment & unprocessed_spectra, libcpp_vector[FASTAEntry] & fasta_db, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, libcpp_vector[libcpp_vector[CrossLinkSpectrumMatch]] & all_top_csms, MSExperiment & spectra)
        Performs the main function of this class, the search for cross-linked peptides
        
        
        :param unprocessed_spectra: The input PeakMap of experimental spectra
        :param fasta_db: The protein database containing targets and decoys
        :param protein_ids: A result vector containing search settings. Should contain one PeptideIdentification
        :param peptide_ids: A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty
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
    OpenPepXLLFAlgorithm_ExitCodes : __OpenPepXLLFAlgorithm_ExitCodes 


class Param:
    """
    Cython implementation of _Param

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Param.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Param()
        """
        ...
    
    @overload
    def __init__(self, in_0: Param ) -> None:
        """
        Cython signature: void Param(Param &)
        """
        ...
    
    @overload
    def setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] , tags: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc, libcpp_vector[libcpp_utf8_string] tags)
        """
        ...
    
    @overload
    def setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] ) -> None:
        """
        Cython signature: void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc)
        """
        ...
    
    @overload
    def setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setValue(libcpp_utf8_string key, ParamValue val)
        """
        ...
    
    def getValue(self, key: Union[bytes, str] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: ParamValue getValue(libcpp_utf8_string key)
        """
        ...
    
    def getValueType(self, key: Union[bytes, str] ) -> int:
        """
        Cython signature: ValueType getValueType(libcpp_utf8_string key)
        """
        ...
    
    def getEntry(self, in_0: Union[bytes, str] ) -> ParamEntry:
        """
        Cython signature: ParamEntry getEntry(libcpp_utf8_string)
        """
        ...
    
    def exists(self, key: Union[bytes, str] ) -> bool:
        """
        Cython signature: bool exists(libcpp_utf8_string key)
        """
        ...
    
    def addTag(self, key: Union[bytes, str] , tag: Union[bytes, str] ) -> None:
        """
        Cython signature: void addTag(libcpp_utf8_string key, libcpp_utf8_string tag)
        """
        ...
    
    def addTags(self, key: Union[bytes, str] , tags: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void addTags(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] tags)
        """
        ...
    
    def hasTag(self, key: Union[bytes, str] , tag: Union[bytes, str] ) -> int:
        """
        Cython signature: int hasTag(libcpp_utf8_string key, libcpp_utf8_string tag)
        """
        ...
    
    def getTags(self, key: Union[bytes, str] ) -> List[bytes]:
        """
        Cython signature: libcpp_vector[libcpp_string] getTags(libcpp_utf8_string key)
        """
        ...
    
    def clearTags(self, key: Union[bytes, str] ) -> None:
        """
        Cython signature: void clearTags(libcpp_utf8_string key)
        """
        ...
    
    def getDescription(self, key: Union[bytes, str] ) -> str:
        """
        Cython signature: libcpp_utf8_output_string getDescription(libcpp_utf8_string key)
        """
        ...
    
    def setSectionDescription(self, key: Union[bytes, str] , desc: Union[bytes, str] ) -> None:
        """
        Cython signature: void setSectionDescription(libcpp_utf8_string key, libcpp_utf8_string desc)
        """
        ...
    
    def getSectionDescription(self, key: Union[bytes, str] ) -> str:
        """
        Cython signature: libcpp_utf8_output_string getSectionDescription(libcpp_utf8_string key)
        """
        ...
    
    def addSection(self, key: Union[bytes, str] , desc: Union[bytes, str] ) -> None:
        """
        Cython signature: void addSection(libcpp_utf8_string key, libcpp_utf8_string desc)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def insert(self, prefix: Union[bytes, str] , param: Param ) -> None:
        """
        Cython signature: void insert(libcpp_utf8_string prefix, Param param)
        """
        ...
    
    def remove(self, key: Union[bytes, str] ) -> None:
        """
        Cython signature: void remove(libcpp_utf8_string key)
        """
        ...
    
    def removeAll(self, prefix: Union[bytes, str] ) -> None:
        """
        Cython signature: void removeAll(libcpp_utf8_string prefix)
        """
        ...
    
    @overload
    def copy(self, prefix: Union[bytes, str] , in_1: bool ) -> Param:
        """
        Cython signature: Param copy(libcpp_utf8_string prefix, bool)
        """
        ...
    
    @overload
    def copy(self, prefix: Union[bytes, str] ) -> Param:
        """
        Cython signature: Param copy(libcpp_utf8_string prefix)
        """
        ...
    
    def merge(self, toMerge: Param ) -> None:
        """
        Cython signature: void merge(Param toMerge)
        """
        ...
    
    @overload
    def setDefaults(self, defaults: Param , prefix: Union[bytes, str] , showMessage: bool ) -> None:
        """
        Cython signature: void setDefaults(Param defaults, libcpp_utf8_string prefix, bool showMessage)
        """
        ...
    
    @overload
    def setDefaults(self, defaults: Param , prefix: Union[bytes, str] ) -> None:
        """
        Cython signature: void setDefaults(Param defaults, libcpp_utf8_string prefix)
        """
        ...
    
    @overload
    def setDefaults(self, defaults: Param ) -> None:
        """
        Cython signature: void setDefaults(Param defaults)
        """
        ...
    
    @overload
    def checkDefaults(self, name: Union[bytes, str] , defaults: Param , prefix: Union[bytes, str] ) -> None:
        """
        Cython signature: void checkDefaults(libcpp_utf8_string name, Param defaults, libcpp_utf8_string prefix)
        """
        ...
    
    @overload
    def checkDefaults(self, name: Union[bytes, str] , defaults: Param ) -> None:
        """
        Cython signature: void checkDefaults(libcpp_utf8_string name, Param defaults)
        """
        ...
    
    def getValidStrings(self, key: Union[bytes, str] ) -> List[Union[bytes, str]]:
        """
        Cython signature: libcpp_vector[libcpp_utf8_string] getValidStrings(libcpp_utf8_string key)
        """
        ...
    
    def setValidStrings(self, key: Union[bytes, str] , strings: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void setValidStrings(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] strings)
        """
        ...
    
    def setMinInt(self, key: Union[bytes, str] , min: int ) -> None:
        """
        Cython signature: void setMinInt(libcpp_utf8_string key, int min)
        """
        ...
    
    def setMaxInt(self, key: Union[bytes, str] , max: int ) -> None:
        """
        Cython signature: void setMaxInt(libcpp_utf8_string key, int max)
        """
        ...
    
    def setMinFloat(self, key: Union[bytes, str] , min: float ) -> None:
        """
        Cython signature: void setMinFloat(libcpp_utf8_string key, double min)
        """
        ...
    
    def setMaxFloat(self, key: Union[bytes, str] , max: float ) -> None:
        """
        Cython signature: void setMaxFloat(libcpp_utf8_string key, double max)
        """
        ...
    
    def __richcmp__(self, other: Param, op: int) -> Any:
        ... 


class Ratio:
    """
    Cython implementation of _Ratio

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Ratio.html>`_
    """
    
    ratio_value_: float
    
    denominator_ref_: Union[bytes, str, String]
    
    numerator_ref_: Union[bytes, str, String]
    
    description_: List[bytes]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Ratio()
        """
        ...
    
    @overload
    def __init__(self, rhs: Ratio ) -> None:
        """
        Cython signature: void Ratio(Ratio rhs)
        """
        ... 


class TraceInfo:
    """
    Cython implementation of _TraceInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TraceInfo.html>`_
    """
    
    name: bytes
    
    description: bytes
    
    opened: bool
    
    @overload
    def __init__(self, n: Union[bytes, str] , d: Union[bytes, str] , o: bool ) -> None:
        """
        Cython signature: void TraceInfo(libcpp_utf8_string n, libcpp_utf8_string d, bool o)
        """
        ...
    
    @overload
    def __init__(self, in_0: TraceInfo ) -> None:
        """
        Cython signature: void TraceInfo(TraceInfo)
        """
        ... 


class XTandemXMLFile:
    """
    Cython implementation of _XTandemXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XTandemXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void XTandemXMLFile()
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , mod_def_set: ModificationDefinitionsSet ) -> None:
        """
        Cython signature: void load(String filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data, ModificationDefinitionsSet & mod_def_set)
        """
        ... 


class __InletType:
    None
    INLETNULL : int
    DIRECT : int
    BATCH : int
    CHROMATOGRAPHY : int
    PARTICLEBEAM : int
    MEMBRANESEPARATOR : int
    OPENSPLIT : int
    JETSEPARATOR : int
    SEPTUM : int
    RESERVOIR : int
    MOVINGBELT : int
    MOVINGWIRE : int
    FLOWINJECTIONANALYSIS : int
    ELECTROSPRAYINLET : int
    THERMOSPRAYINLET : int
    INFUSION : int
    CONTINUOUSFLOWFASTATOMBOMBARDMENT : int
    INDUCTIVELYCOUPLEDPLASMA : int
    MEMBRANE : int
    NANOSPRAY : int
    SIZE_OF_INLETTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class IonOpticsType:
    None
    UNKNOWN : int
    MAGNETIC_DEFLECTION : int
    DELAYED_EXTRACTION : int
    COLLISION_QUADRUPOLE : int
    SELECTED_ION_FLOW_TUBE : int
    TIME_LAG_FOCUSING : int
    REFLECTRON : int
    EINZEL_LENS : int
    FIRST_STABILITY_REGION : int
    FRINGING_FIELD : int
    KINETIC_ENERGY_ANALYZER : int
    STATIC_FIELD : int
    SIZE_OF_IONOPTICSTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __IonizationMethod:
    None
    IONMETHODNULL : int
    ESI : int
    EI : int
    CI : int
    FAB : int
    TSP : int
    LD : int
    FD : int
    FI : int
    PD : int
    SI : int
    TI : int
    API : int
    ISI : int
    CID : int
    CAD : int
    HN : int
    APCI : int
    APPI : int
    ICP : int
    NESI : int
    MESI : int
    SELDI : int
    SEND : int
    FIB : int
    MALDI : int
    MPI : int
    DI : int
    FA : int
    FII : int
    GD_MS : int
    NICI : int
    NRMS : int
    PI : int
    PYMS : int
    REMPI : int
    AI : int
    ASI : int
    AD : int
    AUI : int
    CEI : int
    CHEMI : int
    DISSI : int
    LSI : int
    PEI : int
    SOI : int
    SPI : int
    SUI : int
    VI : int
    AP_MALDI : int
    SILI : int
    SALDI : int
    SIZE_OF_IONIZATIONMETHOD : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class Measure:
    None
    MEASURE_PPM : int
    MEASURE_DA : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __OpenPepXLLFAlgorithm_ExitCodes:
    None
    EXECUTION_OK : int
    ILLEGAL_PARAMETERS : int
    UNEXPECTED_RESULT : int
    INCOMPATIBLE_INPUT_DATA : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __Polarity:
    None
    POLNULL : int
    POSITIVE : int
    NEGATIVE : int
    SIZE_OF_POLARITY : int

    def getMapping(self) -> Dict[int, str]:
       ... 

