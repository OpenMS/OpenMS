from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_OpenMSOSInfo_getBinaryArchitecture() -> Union[bytes, str, String]:
    """
    Cython signature: String getBinaryArchitecture()
    """
    ...

def __static_OpenMSBuildInfo_getBuildType() -> Union[bytes, str, String]:
    """
    Cython signature: String getBuildType()
    """
    ...

def __static_TransformationModelLowess_getDefaultParameters(params: Param ) -> None:
    """
    Cython signature: void getDefaultParameters(Param & params)
    """
    ...

def __static_OpenMSOSInfo_getOSInfo() -> OpenMSOSInfo:
    """
    Cython signature: OpenMSOSInfo getOSInfo()
    """
    ...

def __static_OpenMSBuildInfo_getOpenMPMaxNumThreads() -> int:
    """
    Cython signature: size_t getOpenMPMaxNumThreads()
    """
    ...

def __static_OpenMSBuildInfo_isOpenMPEnabled() -> bool:
    """
    Cython signature: bool isOpenMPEnabled()
    """
    ...

def __static_OpenMSBuildInfo_setOpenMPNumThreads(num_threads: int ) -> None:
    """
    Cython signature: void setOpenMPNumThreads(int num_threads)
    """
    ...


class ConsensusIDAlgorithmPEPIons:
    """
    Cython implementation of _ConsensusIDAlgorithmPEPIons

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmPEPIons.html>`_
      -- Inherits from ['ConsensusIDAlgorithmSimilarity']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusIDAlgorithmPEPIons()
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


class MRMFeatureFilter:
    """
    Cython implementation of _MRMFeatureFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureFilter.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeatureFilter()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeatureFilter ) -> None:
        """
        Cython signature: void MRMFeatureFilter(MRMFeatureFilter &)
        """
        ...
    
    def FilterFeatureMap(self, features: FeatureMap , filter_criteria: MRMFeatureQC , transitions: TargetedExperiment ) -> None:
        """
        Cython signature: void FilterFeatureMap(FeatureMap features, MRMFeatureQC filter_criteria, TargetedExperiment transitions)
        Flags or filters features and subordinates in a FeatureMap
        
        
        :param features: FeatureMap to flag or filter
        :param filter_criteria: MRMFeatureQC class defining QC parameters
        :param transitions: Transitions from a TargetedExperiment
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


class MRMIonSeries:
    """
    Cython implementation of _MRMIonSeries

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMIonSeries.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMIonSeries()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMIonSeries ) -> None:
        """
        Cython signature: void MRMIonSeries(MRMIonSeries &)
        """
        ...
    
    def annotateTransitionCV(self, tr: ReactionMonitoringTransition , annotation: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void annotateTransitionCV(ReactionMonitoringTransition & tr, String annotation)
        Annotates transition with CV terms
        
        
        :param tr: The transition to annotate
        :param annotation: The fragment ion annotation
        """
        ...
    
    def annotateTransition(self, tr: ReactionMonitoringTransition , peptide: Peptide , precursor_mz_threshold: float , product_mz_threshold: float , enable_reannotation: bool , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None:
        """
        Cython signature: void annotateTransition(ReactionMonitoringTransition & tr, Peptide peptide, double precursor_mz_threshold, double product_mz_threshold, bool enable_reannotation, libcpp_vector[String] fragment_types, libcpp_vector[size_t] fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow)
        Annotates transition
        
        
        :param tr: The transition to annotate
        :param peptide: The corresponding peptide
        :param precursor_mz_threshold: The m/z threshold for annotation of the precursor ion
        :param product_mz_threshold: The m/z threshold for annotation of the fragment ion
        :param enable_reannotation: Whether the original (e.g. SpectraST) annotation should be used or reannotation should be conducted
        :param fragment_types: The fragment ion types for reannotation
        :param fragment_charges: The fragment ion charges for reannotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param round_decPow: Round precursor and product m/z values to decimal power (default: -4)
        """
        ... 


class MS2File:
    """
    Cython implementation of _MS2File

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MS2File.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MS2File()
        """
        ...
    
    @overload
    def __init__(self, in_0: MS2File ) -> None:
        """
        Cython signature: void MS2File(MS2File &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & exp)
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


class OMSSACSVFile:
    """
    Cython implementation of _OMSSACSVFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OMSSACSVFile.html>`_

    File adapter for OMSSACSV files
    
    The files contain the results of the OMSSA algorithm in a comma separated manner. This file adapter is able to
    load the data from such a file into the structures of OpenMS
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OMSSACSVFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: OMSSACSVFile ) -> None:
        """
        Cython signature: void OMSSACSVFile(OMSSACSVFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList ) -> None:
        """
        Cython signature: void load(const String & filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data)
        Loads a OMSSA file
        
        The content of the file is stored in `features`
        
        
        :param filename: The name of the file to read from
        :param protein_identification: The protein ProteinIdentification data
        :param id_data: The peptide ids of the file
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ... 


class OpenMSBuildInfo:
    """
    Cython implementation of _OpenMSBuildInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1OpenMSBuildInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenMSBuildInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenMSBuildInfo ) -> None:
        """
        Cython signature: void OpenMSBuildInfo(OpenMSBuildInfo &)
        """
        ...
    
    getBuildType: __static_OpenMSBuildInfo_getBuildType
    
    getOpenMPMaxNumThreads: __static_OpenMSBuildInfo_getOpenMPMaxNumThreads
    
    isOpenMPEnabled: __static_OpenMSBuildInfo_isOpenMPEnabled
    
    setOpenMPNumThreads: __static_OpenMSBuildInfo_setOpenMPNumThreads 


class OpenMSOSInfo:
    """
    Cython implementation of _OpenMSOSInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1OpenMSOSInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenMSOSInfo()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenMSOSInfo ) -> None:
        """
        Cython signature: void OpenMSOSInfo(OpenMSOSInfo &)
        """
        ...
    
    def getOSAsString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOSAsString()
        """
        ...
    
    def getArchAsString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getArchAsString()
        """
        ...
    
    def getOSVersionAsString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOSVersionAsString()
        """
        ...
    
    getBinaryArchitecture: __static_OpenMSOSInfo_getBinaryArchitecture
    
    getOSInfo: __static_OpenMSOSInfo_getOSInfo 


class PeakIndex:
    """
    Cython implementation of _PeakIndex

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakIndex.html>`_

    Index of a peak or feature
    
    This struct can be used to store both peak or feature indices
    """
    
    peak: int
    
    spectrum: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakIndex()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakIndex ) -> None:
        """
        Cython signature: void PeakIndex(PeakIndex &)
        """
        ...
    
    @overload
    def __init__(self, peak: int ) -> None:
        """
        Cython signature: void PeakIndex(size_t peak)
        """
        ...
    
    @overload
    def __init__(self, spectrum: int , peak: int ) -> None:
        """
        Cython signature: void PeakIndex(size_t spectrum, size_t peak)
        """
        ...
    
    def isValid(self) -> bool:
        """
        Cython signature: bool isValid()
        Returns if the current peak ref is valid
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Invalidates the current index
        """
        ...
    
    def getFeature(self, map_: FeatureMap ) -> Feature:
        """
        Cython signature: Feature getFeature(FeatureMap & map_)
        Returns the feature (or consensus feature) corresponding to this index
        
        This method is intended for arrays of features e.g. FeatureMap
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        ...
    
    def getPeak(self, map_: MSExperiment ) -> Peak1D:
        """
        Cython signature: Peak1D getPeak(MSExperiment & map_)
        Returns a peak corresponding to this index
        
        This method is intended for arrays of DSpectra e.g. MSExperiment
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        ...
    
    def getSpectrum(self, map_: MSExperiment ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum getSpectrum(MSExperiment & map_)
        Returns a spectrum corresponding to this index
        
        This method is intended for arrays of DSpectra e.g. MSExperiment
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        ...
    
    def __richcmp__(self, other: PeakIndex, op: int) -> Any:
        ... 


class ProtXMLFile:
    """
    Cython implementation of _ProtXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProtXMLFile.html>`_

    Used to load (storing not supported, yet) ProtXML files
    
    This class is used to load (storing not supported, yet) documents that implement
    the schema of ProtXML files
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ProtXMLFile()
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , protein_ids: ProteinIdentification , peptide_ids: PeptideIdentification ) -> None:
        """
        Cython signature: void load(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids)
        Loads the identifications of an ProtXML file without identifier
        
        The information is read in and the information is stored in the
        corresponding variables
        
        :raises:
          Exception: FileNotFound is thrown if the file could not be found
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , protein_ids: ProteinIdentification , peptide_ids: PeptideIdentification , document_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(String filename, ProteinIdentification & protein_ids, PeptideIdentification & peptide_ids, String document_id)
        """
        ... 


class ProteaseDB:
    """
    Cython implementation of _ProteaseDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteaseDB.html>`_
    """
    
    def getEnzyme(self, name: Union[bytes, str, String] ) -> DigestionEnzymeProtein:
        """
        Cython signature: const DigestionEnzymeProtein * getEnzyme(const String & name)
        """
        ...
    
    def getEnzymeByRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> DigestionEnzymeProtein:
        """
        Cython signature: const DigestionEnzymeProtein * getEnzymeByRegEx(const String & cleavage_regex)
        """
        ...
    
    def getAllNames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllNames(libcpp_vector[String] & all_names)
        """
        ...
    
    def getAllXTandemNames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllXTandemNames(libcpp_vector[String] & all_names)
        Returns all the enzyme names available for XTandem
        """
        ...
    
    def getAllOMSSANames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllOMSSANames(libcpp_vector[String] & all_names)
        Returns all the enzyme names available for OMSSA
        """
        ...
    
    def getAllCometNames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllCometNames(libcpp_vector[String] & all_names)
        Returns all the enzyme names available for Comet
        """
        ...
    
    def getAllMSGFNames(self, all_names: List[bytes] ) -> None:
        """
        Cython signature: void getAllMSGFNames(libcpp_vector[String] & all_names)
        Returns all the enzyme names available for MSGFPlus
        """
        ...
    
    def hasEnzyme(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasEnzyme(const String & name)
        """
        ...
    
    def hasRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasRegEx(const String & cleavage_regex)
        """
        ... 


class RansacModelLinear:
    """
    Cython implementation of _RansacModelLinear

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RansacModelLinear.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RansacModelLinear()
        """
        ...
    
    @overload
    def __init__(self, in_0: RansacModelLinear ) -> None:
        """
        Cython signature: void RansacModelLinear(RansacModelLinear &)
        """
        ... 


class TSE_Match:
    """
    Cython implementation of _TSE_Match

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TSE_Match.html>`_
    """
    
    spectrum: MSSpectrum
    
    score: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TSE_Match()
        """
        ...
    
    @overload
    def __init__(self, in_0: TSE_Match ) -> None:
        """
        Cython signature: void TSE_Match(TSE_Match &)
        """
        ...
    
    @overload
    def __init__(self, spectrum: MSSpectrum , score: float ) -> None:
        """
        Cython signature: void TSE_Match(MSSpectrum & spectrum, double score)
        """
        ... 


class TargetedSpectraExtractor:
    """
    Cython implementation of _TargetedSpectraExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TargetedSpectraExtractor.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TargetedSpectraExtractor()
        """
        ...
    
    @overload
    def __init__(self, in_0: TargetedSpectraExtractor ) -> None:
        """
        Cython signature: void TargetedSpectraExtractor(TargetedSpectraExtractor &)
        """
        ...
    
    def getDefaultParameters(self, in_0: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param &)
        """
        ...
    
    @overload
    def annotateSpectra(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None:
        """
        Cython signature: void annotateSpectra(libcpp_vector[MSSpectrum] &, TargetedExperiment &, libcpp_vector[MSSpectrum] &, FeatureMap &)
        """
        ...
    
    @overload
    def annotateSpectra(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None:
        """
        Cython signature: void annotateSpectra(libcpp_vector[MSSpectrum] &, TargetedExperiment &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    @overload
    def annotateSpectra(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None:
        """
        Cython signature: void annotateSpectra(libcpp_vector[MSSpectrum] &, FeatureMap &, FeatureMap &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    def searchSpectrum(self, in_0: FeatureMap , in_1: FeatureMap , in_2: bool ) -> None:
        """
        Cython signature: void searchSpectrum(FeatureMap &, FeatureMap &, bool)
        """
        ...
    
    def pickSpectrum(self, in_0: MSSpectrum , in_1: MSSpectrum ) -> None:
        """
        Cython signature: void pickSpectrum(MSSpectrum &, MSSpectrum &)
        """
        ...
    
    @overload
    def scoreSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None:
        """
        Cython signature: void scoreSpectra(libcpp_vector[MSSpectrum] &, libcpp_vector[MSSpectrum] &, FeatureMap &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    @overload
    def scoreSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: List[MSSpectrum] ) -> None:
        """
        Cython signature: void scoreSpectra(libcpp_vector[MSSpectrum] &, libcpp_vector[MSSpectrum] &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    @overload
    def selectSpectra(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None:
        """
        Cython signature: void selectSpectra(libcpp_vector[MSSpectrum] &, FeatureMap &, libcpp_vector[MSSpectrum] &, FeatureMap &)
        """
        ...
    
    @overload
    def selectSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] ) -> None:
        """
        Cython signature: void selectSpectra(libcpp_vector[MSSpectrum] &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    @overload
    def extractSpectra(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap , in_4: bool ) -> None:
        """
        Cython signature: void extractSpectra(MSExperiment &, TargetedExperiment &, libcpp_vector[MSSpectrum] &, FeatureMap &, bool)
        """
        ...
    
    @overload
    def extractSpectra(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None:
        """
        Cython signature: void extractSpectra(MSExperiment &, TargetedExperiment &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    @overload
    def extractSpectra(self, in_0: MSExperiment , in_1: FeatureMap , in_2: List[MSSpectrum] ) -> None:
        """
        Cython signature: void extractSpectra(MSExperiment &, FeatureMap &, libcpp_vector[MSSpectrum] &)
        """
        ...
    
    def constructTransitionsList(self, in_0: FeatureMap , in_1: FeatureMap , in_2: TargetedExperiment ) -> None:
        """
        Cython signature: void constructTransitionsList(FeatureMap &, FeatureMap &, TargetedExperiment &)
        """
        ...
    
    def storeSpectraMSP(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None:
        """
        Cython signature: void storeSpectraMSP(const String &, MSExperiment &)
        """
        ...
    
    def mergeFeatures(self, in_0: FeatureMap , in_1: FeatureMap ) -> None:
        """
        Cython signature: void mergeFeatures(FeatureMap &, FeatureMap &)
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


class TransformationModelInterpolated:
    """
    Cython implementation of _TransformationModelInterpolated

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelInterpolated.html>`_
      -- Inherits from ['TransformationModel']
    """
    
    def __init__(self, data: List[TM_DataPoint] , params: Param ) -> None:
        """
        Cython signature: void TransformationModelInterpolated(libcpp_vector[TM_DataPoint] & data, Param & params)
        """
        ...
    
    def getDefaultParameters(self, in_0: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param &)
        Gets the default parameters
        """
        ...
    
    def evaluate(self, value: float ) -> float:
        """
        Cython signature: double evaluate(double value)
        Evaluate the interpolation model at the given value
        
        :param value: The position where the interpolation should be evaluated
        :returns: The interpolated value
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        """
        ...
    
    def weightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void weightData(libcpp_vector[TM_DataPoint] & data)
        Weight the data by the given weight function
        """
        ...
    
    def checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool:
        """
        Cython signature: bool checkValidWeight(const String & weight, libcpp_vector[String] & valid_weights)
        Check for a valid weighting function string
        """
        ...
    
    def weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double weightDatum(double & datum, const String & weight)
        Weight the data according to the weighting function
        """
        ...
    
    def unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double unWeightDatum(double & datum, const String & weight)
        Apply the reverse of the weighting function to the data
        """
        ...
    
    def getValidXWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidXWeights()
        Returns a list of valid x weight function stringss
        """
        ...
    
    def getValidYWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidYWeights()
        Returns a list of valid y weight function strings
        """
        ...
    
    def unWeightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void unWeightData(libcpp_vector[TM_DataPoint] & data)
        Unweight the data by the given weight function
        """
        ...
    
    def checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float:
        """
        Cython signature: double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max)
        Check that the datum is within the valid min and max bounds
        """
        ... 


class TransformationModelLowess:
    """
    Cython implementation of _TransformationModelLowess

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelLowess.html>`_
      -- Inherits from ['TransformationModel']
    """
    
    def __init__(self, data: List[TM_DataPoint] , params: Param ) -> None:
        """
        Cython signature: void TransformationModelLowess(libcpp_vector[TM_DataPoint] & data, Param & params)
        """
        ...
    
    def getDefaultParameters(self, in_0: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param &)
        """
        ...
    
    def evaluate(self, value: float ) -> float:
        """
        Cython signature: double evaluate(double value)
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        """
        ...
    
    def weightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void weightData(libcpp_vector[TM_DataPoint] & data)
        Weight the data by the given weight function
        """
        ...
    
    def checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool:
        """
        Cython signature: bool checkValidWeight(const String & weight, libcpp_vector[String] & valid_weights)
        Check for a valid weighting function string
        """
        ...
    
    def weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double weightDatum(double & datum, const String & weight)
        Weight the data according to the weighting function
        """
        ...
    
    def unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double unWeightDatum(double & datum, const String & weight)
        Apply the reverse of the weighting function to the data
        """
        ...
    
    def getValidXWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidXWeights()
        Returns a list of valid x weight function stringss
        """
        ...
    
    def getValidYWeights(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getValidYWeights()
        Returns a list of valid y weight function strings
        """
        ...
    
    def unWeightData(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void unWeightData(libcpp_vector[TM_DataPoint] & data)
        Unweight the data by the given weight function
        """
        ...
    
    def checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float:
        """
        Cython signature: double checkDatumRange(const double & datum, const double & datum_min, const double & datum_max)
        Check that the datum is within the valid min and max bounds
        """
        ...
    
    getDefaultParameters: __static_TransformationModelLowess_getDefaultParameters 

