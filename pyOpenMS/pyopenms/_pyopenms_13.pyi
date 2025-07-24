from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class ConfidenceScoring:
    """
    Cython implementation of _ConfidenceScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConfidenceScoring.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ConfidenceScoring()
        """
        ...
    
    @overload
    def __init__(self, in_0: ConfidenceScoring ) -> None:
        """
        Cython signature: void ConfidenceScoring(ConfidenceScoring &)
        """
        ...
    
    def initialize(self, targeted: TargetedExperiment , n_decoys: int , n_transitions: int , trafo: TransformationDescription ) -> None:
        """
        Cython signature: void initialize(TargetedExperiment & targeted, size_t n_decoys, size_t n_transitions, TransformationDescription trafo)
        """
        ...
    
    def initializeGlm(self, intercept: float , rt_coef: float , int_coef: float ) -> None:
        """
        Cython signature: void initializeGlm(double intercept, double rt_coef, double int_coef)
        """
        ...
    
    def scoreMap(self, map: FeatureMap ) -> None:
        """
        Cython signature: void scoreMap(FeatureMap & map)
        Score a feature map -> make sure the class is properly initialized
        """
        ... 


class ConsensusIDAlgorithmRanks:
    """
    Cython implementation of _ConsensusIDAlgorithmRanks

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmRanks.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmRanks()
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


class ConsensusMapNormalizerAlgorithmQuantile:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmQuantile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusMapNormalizerAlgorithmQuantile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusMapNormalizerAlgorithmQuantile()
        """
        ...
    
    def normalizeMaps(self, input_map: ConsensusMap ) -> None:
        """
        Cython signature: void normalizeMaps(ConsensusMap & input_map)
        """
        ...
    
    def resample(self, data_in: List[float] , data_out: List[float] , n_resampling_points: int ) -> None:
        """
        Cython signature: void resample(libcpp_vector[double] & data_in, libcpp_vector[double] & data_out, unsigned int n_resampling_points)
        Resamples data_in and writes the results to data_out
        """
        ...
    
    def extractIntensityVectors(self, map_: ConsensusMap , out_intensities: List[List[float]] ) -> None:
        """
        Cython signature: void extractIntensityVectors(ConsensusMap & map_, libcpp_vector[libcpp_vector[double]] & out_intensities)
        Extracts the intensities of the features of the different maps
        """
        ...
    
    def setNormalizedIntensityValues(self, feature_ints: List[List[float]] , map_: ConsensusMap ) -> None:
        """
        Cython signature: void setNormalizedIntensityValues(libcpp_vector[libcpp_vector[double]] & feature_ints, ConsensusMap & map_)
        Writes the intensity values in feature_ints to the corresponding features in map
        """
        ... 


class DataValue:
    """
    Cython implementation of _DataValue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataValue.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DataValue()
        """
        ...
    
    @overload
    def __init__(self, in_0: DataValue ) -> None:
        """
        Cython signature: void DataValue(DataValue &)
        """
        ...
    
    @overload
    def __init__(self, in_0: bytes ) -> None:
        """
        Cython signature: void DataValue(char *)
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void DataValue(const String &)
        """
        ...
    
    @overload
    def __init__(self, in_0: int ) -> None:
        """
        Cython signature: void DataValue(int)
        """
        ...
    
    @overload
    def __init__(self, in_0: float ) -> None:
        """
        Cython signature: void DataValue(double)
        """
        ...
    
    @overload
    def __init__(self, in_0: List[bytes] ) -> None:
        """
        Cython signature: void DataValue(StringList)
        """
        ...
    
    @overload
    def __init__(self, in_0: List[int] ) -> None:
        """
        Cython signature: void DataValue(IntList)
        """
        ...
    
    @overload
    def __init__(self, in_0: List[float] ) -> None:
        """
        Cython signature: void DataValue(DoubleList)
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void DataValue(ParamValue)
        """
        ...
    
    def toStringList(self) -> List[bytes]:
        """
        Cython signature: StringList toStringList()
        """
        ...
    
    def toDoubleList(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] toDoubleList()
        """
        ...
    
    def toIntList(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] toIntList()
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    def toBool(self) -> bool:
        """
        Cython signature: bool toBool()
        """
        ...
    
    def valueType(self) -> int:
        """
        Cython signature: DataType valueType()
        """
        ...
    
    def isEmpty(self) -> int:
        """
        Cython signature: int isEmpty()
        """
        ...
    
    def getUnitType(self) -> int:
        """
        Cython signature: UnitType getUnitType()
        """
        ...
    
    def setUnitType(self, u: int ) -> None:
        """
        Cython signature: void setUnitType(UnitType u)
        """
        ...
    
    def hasUnit(self) -> bool:
        """
        Cython signature: bool hasUnit()
        """
        ...
    
    def getUnit(self) -> int:
        """
        Cython signature: int getUnit()
        """
        ...
    
    def setUnit(self, unit_id: int ) -> None:
        """
        Cython signature: void setUnit(int unit_id)
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ... 


class LabeledPairFinder:
    """
    Cython implementation of _LabeledPairFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LabeledPairFinder.html>`_
      -- Inherits from ['BaseGroupFinder']

    The LabeledPairFinder allows the matching of labeled features (features with a fixed distance)
    
    Finds feature pairs that have a defined distance in RT and m/z in the same map
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void LabeledPairFinder()
        """
        ...
    
    def run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None:
        """
        Cython signature: void run(libcpp_vector[ConsensusMap] & input_maps, ConsensusMap & result_map)
        Runs the LabeledPairFinder algorithm
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


class MRMDecoy:
    """
    Cython implementation of _MRMDecoy

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMDecoy.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMDecoy()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMDecoy ) -> None:
        """
        Cython signature: void MRMDecoy(MRMDecoy &)
        """
        ...
    
    def generateDecoys(self, exp: TargetedExperiment , dec: TargetedExperiment , method: Union[bytes, str, String] , aim_decoy_fraction: float , switchKR: bool , decoy_tag: Union[bytes, str, String] , max_attempts: int , identity_threshold: float , precursor_mz_shift: float , product_mz_shift: float , product_mz_threshold: float , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None:
        """
        Cython signature: void generateDecoys(TargetedExperiment & exp, TargetedExperiment & dec, String method, double aim_decoy_fraction, bool switchKR, String decoy_tag, int max_attempts, double identity_threshold, double precursor_mz_shift, double product_mz_shift, double product_mz_threshold, libcpp_vector[String] fragment_types, libcpp_vector[size_t] fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow)
        Generate decoys from a TargetedExperiment
        
        Will generate decoy peptides for each target peptide provided in exp and
        write them into the decoy experiment
        
        Valid methods: shuffle, reverse, pseudo-reverse
        
        If theoretical is true, the target transitions will be returned but their
        masses will be adjusted to match the theoretical value of the fragment ion
        that is the most likely explanation for the product
        
        `mz_threshold` is used for the matching of theoretical ion series to the observed one
        
        To generate decoys with different precursor mass, use the "switchKR" flag
        which switches terminal K/R (switches K to R and R to K). This generates
        different precursor m/z and ensures that the y ion series has a different
        mass. For a description of the procedure, see (supplemental material)
        
        Bruderer et al. Mol Cell Proteomics. 2017. 10.1074/mcp.RA117.000314.
        """
        ...
    
    def findFixedResidues(self, sequence: Union[bytes, str, String] , keepN: bool , keepC: bool , keep_const_pattern: Union[bytes, str, String] ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] findFixedResidues(const String & sequence, bool keepN, bool keepC, const String & keep_const_pattern)
        Find all residues in a sequence that should not be reversed / shuffled
        
        
        :param sequence: The amino acid sequence
        :param keepN: Whether to keep N terminus constant
        :param keepC: Whether to keep C terminus constant
        :param keep_const_pattern: A string containing the AA to not change (e.g. 'KRP')
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


class MapAlignmentEvaluationAlgorithmPrecision:
    """
    Cython implementation of _MapAlignmentEvaluationAlgorithmPrecision

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentEvaluationAlgorithmPrecision.html>`_
      -- Inherits from ['MapAlignmentEvaluationAlgorithm']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MapAlignmentEvaluationAlgorithmPrecision()
        """
        ... 


class MzMLFile:
    """
    Cython implementation of _MzMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLFile.html>`_
      -- Inherits from ['ProgressLogger']

    File adapter for MzML files
    
    Provides methods to load and store MzML files.
    PeakFileOptions allow to load a reduced subset of the data into an MSExperiment.
    
    See help(MSExperiment) how data is stored after loading.
    See help(PeakFileOptions) for available options.
    
    Usage:
    
    .. code-block:: python
    
      exp = MSExperiment()
      MzMLFile().load("test.mzML", exp)
      spec = []
      for s in exp.getSpectra():
        if s.getMSLevel() != 1:
          spec.append(s)
      exp.setSpectra(spec)
      MzMLFile().store("filtered.mzML", exp)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzMLFile ) -> None:
        """
        Cython signature: void MzMLFile(MzMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment &)
        Loads from an MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void store(const String & filename, MSExperiment &)
        Stores a MSExperiment in an MzML file
        """
        ...
    
    def storeBuffer(self, output: String , exp: MSExperiment ) -> None:
        """
        Cython signature: void storeBuffer(String & output, MSExperiment exp)
        Stores a map in an output string
        
        
        :param output: An empty string to store the result
        :param exp: Has to be an MSExperiment
        """
        ...
    
    def loadBuffer(self, input: Union[bytes, str, String] , exp: MSExperiment ) -> None:
        """
        Cython signature: void loadBuffer(const String & input, MSExperiment & exp)
        Loads a map from a MzML file stored in a buffer (in memory)
        
        
        :param buffer: The buffer with the data (i.e. string with content of an mzML file)
        :param exp: Is an MSExperiment
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
        """
        ...
    
    def setOptions(self, in_0: PeakFileOptions ) -> None:
        """
        Cython signature: void setOptions(PeakFileOptions)
        Set PeakFileOptions to perform filtering during loading. E.g., to load only MS1 spectra or meta data only
        """
        ...
    
    def isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool:
        """
        Cython signature: bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
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


class MzMLSpectrumDecoder:
    """
    Cython implementation of _MzMLSpectrumDecoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLSpectrumDecoder.html>`_

    A class to decode input strings that contain an mzML chromatogram or spectrum tag
    
    It uses xercesc to parse a string containing either a exactly one mzML
    spectrum or chromatogram (from <chromatogram> to </chromatogram> or
    <spectrum> to </spectrum> tag). It returns the data contained in the
    binaryDataArray for Intensity / mass-to-charge or Intensity / time
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzMLSpectrumDecoder()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzMLSpectrumDecoder ) -> None:
        """
        Cython signature: void MzMLSpectrumDecoder(MzMLSpectrumDecoder &)
        """
        ...
    
    def domParseChromatogram(self, in_: Union[bytes, str, String] , cptr: _Interfaces_Chromatogram ) -> None:
        """
        Cython signature: void domParseChromatogram(String in_, shared_ptr[_Interfaces_Chromatogram] & cptr)
        Extract data from a string which contains a full mzML chromatogram
        
        Extracts data from the input string which is expected to contain exactly
        one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
        function will extract the contained binaryDataArray and provide the
        result as Chromatogram
        
        
        :param in: Input string containing the raw XML
        :param cptr: Resulting chromatogram
        """
        ...
    
    def domParseSpectrum(self, in_: Union[bytes, str, String] , cptr: _Interfaces_Spectrum ) -> None:
        """
        Cython signature: void domParseSpectrum(String in_, shared_ptr[_Interfaces_Spectrum] & cptr)
        Extract data from a string which contains a full mzML spectrum
        
        Extracts data from the input string which is expected to contain exactly
        one <spectrum> tag (from <spectrum> to </spectrum>). This function will
        extract the contained binaryDataArray and provide the result as Spectrum
        
        
        :param in: Input string containing the raw XML
        :param cptr: Resulting spectrum
        """
        ...
    
    def setSkipXMLChecks(self, only: bool ) -> None:
        """
        Cython signature: void setSkipXMLChecks(bool only)
        Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead
        """
        ... 


class ProteinInference:
    """
    Cython implementation of _ProteinInference

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinInference.html>`_

    [experimental class] given a peptide quantitation, infer corresponding protein quantities
    
    Infers protein ratios from peptide ratios (currently using unique peptides only).
    Use the IDMapper class to add protein and peptide information to a
    quantitative ConsensusMap prior to this step
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteinInference()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteinInference ) -> None:
        """
        Cython signature: void ProteinInference(ProteinInference &)
        """
        ...
    
    def infer(self, consensus_map: ConsensusMap , reference_map: int ) -> None:
        """
        Cython signature: void infer(ConsensusMap & consensus_map, unsigned int reference_map)
        Given a peptide quantitation, infer corresponding protein quantities
        
        Infers protein ratios from peptide ratios (currently using unique peptides only).
        Use the IDMapper class to add protein and peptide information to a
        quantitative ConsensusMap prior to this step
        
        
        :param consensus_map: Peptide quantitation with ProteinIdentifications attached, where protein quantitation will be attached
        :param reference_map: Index of (iTRAQ) reference channel within the consensus map
        """
        ... 


class SpectrumAccessOpenMSInMemory:
    """
    Cython implementation of _SpectrumAccessOpenMSInMemory

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMSInMemory.html>`_
      -- Inherits from ['ISpectrumAccess']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSInMemory()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMS ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMS &)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMSCached ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSCached &)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMSInMemory ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSInMemory(SpectrumAccessOpenMSInMemory &)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessQuadMZTransforming ) -> None:
        """
        Cython signature: void SpectrumAccessOpenMSInMemory(SpectrumAccessQuadMZTransforming &)
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


class SpectrumLookup:
    """
    Cython implementation of _SpectrumLookup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumLookup.html>`_
    """
    
    rt_tolerance: float
    
    def __init__(self) -> None:
        """
        Cython signature: void SpectrumLookup()
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Check if any spectra were set
        """
        ...
    
    def readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readSpectra(MSExperiment spectra, String scan_regexp)
        Read and index spectra for later look-up
        
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
        """
        ...
    
    def findByRT(self, rt: float ) -> int:
        """
        Cython signature: size_t findByRT(double rt)
        Look up spectrum by retention time (RT)
        
        :param rt: Retention time to look up
        :returns: Index of the spectrum that matched
        """
        ...
    
    def findByNativeID(self, native_id: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t findByNativeID(String native_id)
        Look up spectrum by native ID
        
        :param native_id: Native ID to look up
        :returns: Index of the spectrum that matched
        """
        ...
    
    def findByIndex(self, index: int , count_from_one: bool ) -> int:
        """
        Cython signature: size_t findByIndex(size_t index, bool count_from_one)
        Look up spectrum by index (position in the vector of spectra)
        
        :param index: Index to look up
        :param count_from_one: Do indexes start counting at one (default zero)?
        :returns: Index of the spectrum that matched
        """
        ...
    
    def findByScanNumber(self, scan_number: int ) -> int:
        """
        Cython signature: size_t findByScanNumber(size_t scan_number)
        Look up spectrum by scan number (extracted from the native ID)
        
        :param scan_number: Scan number to look up
        :returns: Index of the spectrum that matched
        """
        ...
    
    def findByReference(self, spectrum_ref: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t findByReference(String spectrum_ref)
        Look up spectrum by reference
        
        :param spectrum_ref: Spectrum reference to parse
        :returns: Index of the spectrum that matched
        """
        ...
    
    def addReferenceFormat(self, regexp: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addReferenceFormat(String regexp)
        Register a possible format for a spectrum reference
        
        :param regexp: Regular expression defining the format
        """
        ...
    
    def extractScanNumber(self, native_id: Union[bytes, str, String] , native_id_type_accession: Union[bytes, str, String] ) -> int:
        """
        Cython signature: int extractScanNumber(const String & native_id, const String & native_id_type_accession)
        """
        ... 


class SpectrumRangeManager:
    """
    Cython implementation of _SpectrumRangeManager

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumRangeManager.html>`_

    Advanced range manager for MS spectra with separate ranges for each MS level
    
    This class extends the basic RangeManager to provide separate range tracking for different MS levels
    (MS1, MS2, etc.). It manages four types of ranges:
    - m/z (mass-to-charge ratio)
    - intensity
    - retention time (RT)
    - ion mobility
    
    A global range is tracked for all MS levels, and additional ranges are maintained for each specific MS level.
    This allows for efficient querying of ranges for specific MS levels, which is useful for visualization,
    filtering, and processing operations that need to work with specific MS levels.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumRangeManager()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumRangeManager ) -> None:
        """
        Cython signature: void SpectrumRangeManager(SpectrumRangeManager &)
        """
        ...
    
    def clearRanges(self) -> None:
        """
        Cython signature: void clearRanges()
        """
        ...
    
    def getMSLevels(self) -> Set[int]:
        """
        Cython signature: libcpp_set[unsigned int] getMSLevels()
        """
        ...
    
    def extendRT(self, rt: float , ms_level: int ) -> None:
        """
        Cython signature: void extendRT(double rt, unsigned int ms_level)
        """
        ...
    
    def extendMZ(self, mz: float , ms_level: int ) -> None:
        """
        Cython signature: void extendMZ(double mz, unsigned int ms_level)
        """
        ...
    
    def extendUnsafe(self, spectrum: MSSpectrum , ms_level: int ) -> None:
        """
        Cython signature: void extendUnsafe(const MSSpectrum & spectrum, unsigned int ms_level)
        """
        ...
    
    def getMinRT(self) -> float:
        """
        Cython signature: double getMinRT()
        """
        ...
    
    def getMaxRT(self) -> float:
        """
        Cython signature: double getMaxRT()
        """
        ...
    
    def getMinMZ(self) -> float:
        """
        Cython signature: double getMinMZ()
        """
        ...
    
    def getMaxMZ(self) -> float:
        """
        Cython signature: double getMaxMZ()
        """
        ...
    
    def getMinIntensity(self) -> float:
        """
        Cython signature: double getMinIntensity()
        """
        ...
    
    def getMaxIntensity(self) -> float:
        """
        Cython signature: double getMaxIntensity()
        """
        ...
    
    def getMinMobility(self) -> float:
        """
        Cython signature: double getMinMobility()
        """
        ...
    
    def getMaxMobility(self) -> float:
        """
        Cython signature: double getMaxMobility()
        """
        ...
    
    def __richcmp__(self, other: SpectrumRangeManager, op: int) -> Any:
        ... 


class DataType:
    None
    STRING_VALUE : int
    INT_VALUE : int
    DOUBLE_VALUE : int
    STRING_LIST : int
    INT_LIST : int
    DOUBLE_LIST : int
    EMPTY_VALUE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class UnitType:
    None
    UNIT_ONTOLOGY : int
    MS_ONTOLOGY : int
    OTHER : int

    def getMapping(self) -> Dict[int, str]:
       ... 

