from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_SpectrumMetaDataLookup_addMissingIMToPeptideIDs(in_0: PeptideIdentificationList , exp: MSExperiment ) -> bool:
    """
    Cython signature: bool addMissingIMToPeptideIDs(PeptideIdentificationList, MSExperiment exp)
    """
    ...

def __static_SpectrumMetaDataLookup_addMissingRTsToPeptideIDs(in_0: PeptideIdentificationList , exp: MSExperiment ) -> bool:
    """
    Cython signature: bool addMissingRTsToPeptideIDs(PeptideIdentificationList, MSExperiment exp)
    """
    ...

def __static_SpectrumMetaDataLookup_addMissingSpectrumReferences(in_0: PeptideIdentificationList , filename: Union[bytes, str, String] , stop_on_error: bool , override_spectra_data: bool , override_spectra_references: bool , proteins: List[ProteinIdentification] ) -> bool:
    """
    Cython signature: bool addMissingSpectrumReferences(PeptideIdentificationList, String filename, bool stop_on_error, bool override_spectra_data, bool override_spectra_references, libcpp_vector[ProteinIdentification] proteins)
    """
    ...

def __static_SpectrumMetaDataLookup_getSpectrumMetaData(spectrum: MSSpectrum , meta: SpectrumMetaData ) -> None:
    """
    Cython signature: void getSpectrumMetaData(MSSpectrum spectrum, SpectrumMetaData & meta)
    """
    ...


class BilinearInterpolation:
    """
    Cython implementation of _BilinearInterpolation[double,double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1BilinearInterpolation[double,double].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BilinearInterpolation()
        """
        ...
    
    @overload
    def __init__(self, in_0: BilinearInterpolation ) -> None:
        """
        Cython signature: void BilinearInterpolation(BilinearInterpolation &)
        """
        ...
    
    def value(self, arg_pos_0: float , arg_pos_1: float ) -> float:
        """
        Cython signature: double value(double arg_pos_0, double arg_pos_1)
        """
        ...
    
    def addValue(self, arg_pos_0: float , arg_pos_1: float , arg_value: float ) -> None:
        """
        Cython signature: void addValue(double arg_pos_0, double arg_pos_1, double arg_value)
        Performs bilinear resampling. The arg_value is split up and added to the data points around arg_pos. ("forward resampling")
        """
        ...
    
    def getData(self) -> MatrixDouble:
        """
        Cython signature: MatrixDouble getData()
        """
        ...
    
    def setData(self, data: MatrixDouble ) -> None:
        """
        Cython signature: void setData(MatrixDouble & data)
        Assigns data to the internal random access container storing the data. SourceContainer must be assignable to ContainerType
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def key2index_0(self, pos: float ) -> float:
        """
        Cython signature: double key2index_0(double pos)
        The transformation from "outside" to "inside" coordinates
        """
        ...
    
    def index2key_0(self, pos: float ) -> float:
        """
        Cython signature: double index2key_0(double pos)
        The transformation from "inside" to "outside" coordinates
        """
        ...
    
    def key2index_1(self, pos: float ) -> float:
        """
        Cython signature: double key2index_1(double pos)
        The transformation from "outside" to "inside" coordinates
        """
        ...
    
    def index2key_1(self, pos: float ) -> float:
        """
        Cython signature: double index2key_1(double pos)
        The transformation from "inside" to "outside" coordinates
        """
        ...
    
    def getScale_0(self) -> float:
        """
        Cython signature: double getScale_0()
        """
        ...
    
    def setScale_0(self, scale: float ) -> None:
        """
        Cython signature: void setScale_0(double & scale)
        """
        ...
    
    def getScale_1(self) -> float:
        """
        Cython signature: double getScale_1()
        """
        ...
    
    def setScale_1(self, scale: float ) -> None:
        """
        Cython signature: void setScale_1(double & scale)
        """
        ...
    
    def getOffset_0(self) -> float:
        """
        Cython signature: double getOffset_0()
        Accessor. "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        """
        ...
    
    def setOffset_0(self, offset: float ) -> None:
        """
        Cython signature: void setOffset_0(double & offset)
        """
        ...
    
    def getOffset_1(self) -> float:
        """
        Cython signature: double getOffset_1()
        Accessor. "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        """
        ...
    
    def setOffset_1(self, offset: float ) -> None:
        """
        Cython signature: void setOffset_1(double & offset)
        """
        ...
    
    @overload
    def setMapping_0(self, scale: float , inside: float , outside: float ) -> None:
        """
        Cython signature: void setMapping_0(double & scale, double & inside, double & outside)
        """
        ...
    
    @overload
    def setMapping_0(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None:
        """
        Cython signature: void setMapping_0(double & inside_low, double & outside_low, double & inside_high, double & outside_high)
        """
        ...
    
    @overload
    def setMapping_1(self, scale: float , inside: float , outside: float ) -> None:
        """
        Cython signature: void setMapping_1(double & scale, double & inside, double & outside)
        """
        ...
    
    @overload
    def setMapping_1(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None:
        """
        Cython signature: void setMapping_1(double & inside_low, double & outside_low, double & inside_high, double & outside_high)
        """
        ...
    
    def getInsideReferencePoint_0(self) -> float:
        """
        Cython signature: double getInsideReferencePoint_0()
        """
        ...
    
    def getInsideReferencePoint_1(self) -> float:
        """
        Cython signature: double getInsideReferencePoint_1()
        """
        ...
    
    def getOutsideReferencePoint_0(self) -> float:
        """
        Cython signature: double getOutsideReferencePoint_0()
        """
        ...
    
    def getOutsideReferencePoint_1(self) -> float:
        """
        Cython signature: double getOutsideReferencePoint_1()
        """
        ...
    
    def supportMin_0(self) -> float:
        """
        Cython signature: double supportMin_0()
        Lower boundary of the support, in "outside" coordinates
        """
        ...
    
    def supportMin_1(self) -> float:
        """
        Cython signature: double supportMin_1()
        Lower boundary of the support, in "outside" coordinates
        """
        ...
    
    def supportMax_0(self) -> float:
        """
        Cython signature: double supportMax_0()
        Upper boundary of the support, in "outside" coordinates
        """
        ...
    
    def supportMax_1(self) -> float:
        """
        Cython signature: double supportMax_1()
        Upper boundary of the support, in "outside" coordinates
        """
        ... 


class CachedMzMLHandler:
    """
    Cython implementation of _CachedMzMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1CachedMzMLHandler.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CachedMzMLHandler()
        An internal class that handles single spectra and chromatograms
        """
        ...
    
    @overload
    def __init__(self, in_0: CachedMzMLHandler ) -> None:
        """
        Cython signature: void CachedMzMLHandler(CachedMzMLHandler &)
        """
        ...
    
    def writeMemdump(self, exp: MSExperiment , out: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void writeMemdump(MSExperiment exp, String out)
        Write complete spectra as a dump to the disk
        """
        ...
    
    def writeMetadata(self, exp: MSExperiment , out_meta: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void writeMetadata(MSExperiment exp, String out_meta)
        Write only the meta data of an MSExperiment
        """
        ...
    
    def readMemdump(self, exp: MSExperiment , filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readMemdump(MSExperiment exp, String filename)
        Read all spectra from a dump from the disk
        """
        ...
    
    def getSpectraIndex(self) -> List[streampos]:
        """
        Cython signature: libcpp_vector[streampos] getSpectraIndex()
        """
        ...
    
    def getChromatogramIndex(self) -> List[streampos]:
        """
        Cython signature: libcpp_vector[streampos] getChromatogramIndex()
        """
        ...
    
    def createMemdumpIndex(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void createMemdumpIndex(String filename)
        Create an index on the location of all the spectra and chromatograms
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


class IndexedMzMLDecoder:
    """
    Cython implementation of _IndexedMzMLDecoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLDecoder.html>`_

    A class to analyze indexedmzML files and extract the offsets of individual tags
    
    Specifically, this class allows one to extract the offsets of the <indexList>
    tag and of all <spectrum> and <chromatogram> tag using the indices found at
    the end of the indexedmzML XML structure
    
    While findIndexListOffset tries extracts the offset of the indexList tag from
    the last 1024 bytes of the file, this offset allows the function parseOffsets
    to extract all elements contained in the <indexList> tag and thus get access
    to all spectra and chromatogram offsets
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IndexedMzMLDecoder()
        """
        ...
    
    @overload
    def __init__(self, in_0: IndexedMzMLDecoder ) -> None:
        """
        Cython signature: void IndexedMzMLDecoder(IndexedMzMLDecoder &)
        """
        ...
    
    def findIndexListOffset(self, in_: Union[bytes, str, String] , buffersize: int ) -> streampos:
        """
        Cython signature: streampos findIndexListOffset(String in_, int buffersize)
        Tries to extract the indexList offset from an indexedmzML\n
        
        This function reads by default the last few (1024) bytes of the given
        input file and tries to read the content of the <indexListOffset> tag
        The idea is that somewhere in the last parts of the file specified by the
        input string, the string <indexListOffset>xxx</indexListOffset> occurs
        This function returns the xxx part converted to an integer\n
        
        Since this function cannot determine where it will start reading
        the XML, no regular XML parser can be used for this. Therefore it uses
        regex to do its job. It matches the <indexListOffset> part and any
        numerical characters that follow
        
        
        :param in: Filename of the input indexedmzML file
        :param buffersize: How many bytes of the input file should be searched for the tag
        :return: A positive integer containing the content of the indexListOffset tag, returns -1 in case of failure no tag was found (you can re-try with a larger buffersize but most likely its not an indexed mzML). Using -1 is what the reference docu recommends: http://en.cppreference.com/w/cpp/io/streamoff
        :raises:
          Exception: FileNotFound is thrown if file cannot be found
        :raises:
          Exception: ParseError if offset cannot be parsed
        """
        ... 


class IsobaricIsotopeCorrector:
    """
    Cython implementation of _IsobaricIsotopeCorrector

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricIsotopeCorrector.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsobaricIsotopeCorrector()
        """
        ...
    
    @overload
    def __init__(self, in_0: IsobaricIsotopeCorrector ) -> None:
        """
        Cython signature: void IsobaricIsotopeCorrector(IsobaricIsotopeCorrector &)
        """
        ...
    
    @overload
    def correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqEightPlexQuantitationMethod ) -> IsobaricQuantifierStatistics:
        """
        Cython signature: IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, ItraqEightPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqFourPlexQuantitationMethod ) -> IsobaricQuantifierStatistics:
        """
        Cython signature: IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, ItraqFourPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTSixPlexQuantitationMethod ) -> IsobaricQuantifierStatistics:
        """
        Cython signature: IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, TMTSixPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTTenPlexQuantitationMethod ) -> IsobaricQuantifierStatistics:
        """
        Cython signature: IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, TMTTenPlexQuantitationMethod * quant_method)
        """
        ... 


class IsobaricQuantifierStatistics:
    """
    Cython implementation of _IsobaricQuantifierStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricQuantifierStatistics.html>`_
    """
    
    channel_count: int
    
    iso_number_ms2_negative: int
    
    iso_number_reporter_negative: int
    
    iso_number_reporter_different: int
    
    iso_solution_different_intensity: float
    
    iso_total_intensity_negative: float
    
    number_ms2_total: int
    
    number_ms2_empty: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsobaricQuantifierStatistics()
        """
        ...
    
    @overload
    def __init__(self, in_0: IsobaricQuantifierStatistics ) -> None:
        """
        Cython signature: void IsobaricQuantifierStatistics(IsobaricQuantifierStatistics &)
        """
        ...
    
    def reset(self) -> None:
        """
        Cython signature: void reset()
        """
        ... 


class MetaboliteFeatureDeconvolution:
    """
    Cython implementation of _MetaboliteFeatureDeconvolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboliteFeatureDeconvolution.html>`_
      -- Inherits from ['DefaultParamHandler']

    An algorithm to decharge small molecule features (i.e. as found by FeatureFinder)
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaboliteFeatureDeconvolution()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaboliteFeatureDeconvolution ) -> None:
        """
        Cython signature: void MetaboliteFeatureDeconvolution(MetaboliteFeatureDeconvolution &)
        """
        ...
    
    def compute(self, fm_in: FeatureMap , fm_out: FeatureMap , cons_map: ConsensusMap , cons_map_p: ConsensusMap ) -> None:
        """
        Cython signature: void compute(FeatureMap & fm_in, FeatureMap & fm_out, ConsensusMap & cons_map, ConsensusMap & cons_map_p)
        Compute a zero-charge feature map from a set of charged features
        
        Find putative ChargePairs, then score them and hand over to ILP
        
        
        :param fm_in: Input feature-map
        :param fm_out: Output feature-map (sorted by position and augmented with user params)
        :param cons_map: Output of grouped features belonging to a charge group
        :param cons_map_p: Output of paired features connected by an edge
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
    CHARGEMODE_MFD : __CHARGEMODE_MFD 


class MzIdentMLFile:
    """
    Cython implementation of _MzIdentMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzIdentMLFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzIdentMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzIdentMLFile ) -> None:
        """
        Cython signature: void MzIdentMLFile(MzIdentMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & poid, PeptideIdentificationList & peid)
        Loads the identifications from a MzIdentML file
        
        
        :param filename: File name of the file to be checked
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsin
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & poid, PeptideIdentificationList & peid)
        Stores the identifications in a MzIdentML file
        
        
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be created
        """
        ...
    
    def isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool:
        """
        Cython signature: bool isSemanticallyValid(String filename, StringList errors, StringList warnings)
        Checks if a file is valid with respect to the mapping file and the controlled vocabulary
        
        
        :param filename: File name of the file to be checked
        :param errors: Errors during the validation are returned in this output parameter
        :param warnings: Warnings during the validation are returned in this output parameter
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
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


class Normalizer:
    """
    Cython implementation of _Normalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Normalizer.html>`_
      -- Inherits from ['DefaultParamHandler']

    Normalizes the peak intensities spectrum-wise
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Normalizer()
        """
        ...
    
    @overload
    def __init__(self, in_0: Normalizer ) -> None:
        """
        Cython signature: void Normalizer(Normalizer)
        """
        ...
    
    def filterSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterSpectrum(MSSpectrum & spec)
        Normalizes the spectrum
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        Normalizes the peak spectrum
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
        Normalizes the peak map
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


class OPXLDataStructs:
    """
    Cython implementation of _OPXLDataStructs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLDataStructs.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OPXLDataStructs()
        """
        ...
    
    @overload
    def __init__(self, in_0: OPXLDataStructs ) -> None:
        """
        Cython signature: void OPXLDataStructs(OPXLDataStructs &)
        """
        ...
    PeptidePosition : __PeptidePosition
    ProteinProteinCrossLinkType : __ProteinProteinCrossLinkType 


class PeptideAndProteinQuant:
    """
    Cython implementation of _PeptideAndProteinQuant

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant()
        Helper class for peptide and protein quantification based on feature data annotated with IDs
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideAndProteinQuant ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant(PeptideAndProteinQuant &)
        """
        ...
    
    @overload
    def readQuantData(self, map_in: FeatureMap , ed: ExperimentalDesign ) -> None:
        """
        Cython signature: void readQuantData(FeatureMap & map_in, ExperimentalDesign & ed)
        Read quantitative data from a feature map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        ...
    
    @overload
    def readQuantData(self, map_in: ConsensusMap , ed: ExperimentalDesign ) -> None:
        """
        Cython signature: void readQuantData(ConsensusMap & map_in, ExperimentalDesign & ed)
        Read quantitative data from a consensus map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        ...
    
    @overload
    def readQuantData(self, proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList , ed: ExperimentalDesign ) -> None:
        """
        Cython signature: void readQuantData(libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList & peptides, ExperimentalDesign & ed)
        Read quantitative data from identification results (for quantification via spectral counting)
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        ...
    
    def quantifyPeptides(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void quantifyPeptides(PeptideIdentificationList & peptides)
        Compute peptide abundances
        
        Based on quantitative data for individual charge states (in member `pep_quant_`), overall abundances for peptides are computed (and stored again in `pep_quant_`)
        Quantitative data must first be read via readQuantData()
        Optional (peptide-level) protein inference information (e.g. from Fido or ProteinProphet) can be supplied via `peptides`. In that case, peptide-to-protein associations - the basis for protein-level quantification - will also be read from `peptides`!
        """
        ...
    
    def quantifyProteins(self, proteins: ProteinIdentification ) -> None:
        """
        Cython signature: void quantifyProteins(ProteinIdentification & proteins)
        Compute protein abundances
        
        Peptide abundances must be computed first with quantifyPeptides(). Optional protein inference information (e.g. from Fido or ProteinProphet) can be supplied via `proteins`
        """
        ...
    
    def getStatistics(self) -> PeptideAndProteinQuant_Statistics:
        """
        Cython signature: PeptideAndProteinQuant_Statistics getStatistics()
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


class PeptideAndProteinQuant_PeptideData:
    """
    Cython implementation of _PeptideAndProteinQuant_PeptideData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_PeptideData.html>`_
    """
    
    accessions: Set[bytes]
    
    psm_count: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_PeptideData()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideAndProteinQuant_PeptideData ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_PeptideData(PeptideAndProteinQuant_PeptideData &)
        """
        ... 


class PeptideAndProteinQuant_ProteinData:
    """
    Cython implementation of _PeptideAndProteinQuant_ProteinData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_ProteinData.html>`_
    """
    
    psm_count: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_ProteinData()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideAndProteinQuant_ProteinData ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_ProteinData(PeptideAndProteinQuant_ProteinData &)
        """
        ... 


class PeptideAndProteinQuant_Statistics:
    """
    Cython implementation of _PeptideAndProteinQuant_Statistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_Statistics.html>`_
    """
    
    n_samples: int
    
    quant_proteins: int
    
    too_few_peptides: int
    
    quant_peptides: int
    
    total_peptides: int
    
    quant_features: int
    
    total_features: int
    
    blank_features: int
    
    ambig_features: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_Statistics()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideAndProteinQuant_Statistics ) -> None:
        """
        Cython signature: void PeptideAndProteinQuant_Statistics(PeptideAndProteinQuant_Statistics &)
        """
        ... 


class ProgressLogger:
    """
    Cython implementation of _ProgressLogger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProgressLogger.html>`_

    Base class for all classes that want to report their progress
    
    Per default the progress log is disabled. Use setLogType to enable it
    
    Use startProgress, setProgress and endProgress for the actual logging
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProgressLogger()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProgressLogger ) -> None:
        """
        Cython signature: void ProgressLogger(ProgressLogger &)
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


class SpectrumMetaData:
    """
    Cython implementation of _SpectrumMetaData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumMetaData.html>`_
    """
    
    rt: float
    
    precursor_rt: float
    
    precursor_mz: float
    
    precursor_charge: int
    
    ms_level: int
    
    scan_number: int
    
    native_id: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumMetaData()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumMetaData ) -> None:
        """
        Cython signature: void SpectrumMetaData(SpectrumMetaData &)
        """
        ... 


class SpectrumMetaDataLookup:
    """
    Cython implementation of _SpectrumMetaDataLookup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumMetaDataLookup.html>`_
      -- Inherits from ['SpectrumLookup']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void SpectrumMetaDataLookup()
        """
        ...
    
    @overload
    def readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] , get_precursor_rt: bool ) -> None:
        """
        Cython signature: void readSpectra(MSExperiment spectra, String scan_regexp, bool get_precursor_rt)
        Read spectra and store their meta data
        
        :param SpectrumContainer: Spectrum container class, must support `size` and `operator[]`
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")
        :param get_precursor_rt: Assign precursor retention times? (This relies on all precursor spectra being present and in the right order.)
        """
        ...
    
    @overload
    def readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readSpectra(MSExperiment spectra, String scan_regexp)
        Read and index spectra for later look-up
        
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
        """
        ...
    
    @overload
    def getSpectrumMetaData(self, index: int , meta: SpectrumMetaData ) -> None:
        """
        Cython signature: void getSpectrumMetaData(size_t index, SpectrumMetaData & meta)
        Look up meta data of a spectrum
        
        :param index: Index of the spectrum
        :param meta: Meta data output
        """
        ...
    
    @overload
    def getSpectrumMetaData(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData ) -> None:
        """
        Cython signature: void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData & meta)
        Extract meta data from a spectrum
        
        :param spectrum: Spectrum input
        :param meta: Meta data output
        :param scan_regexp: Regular expression for extracting scan number from spectrum native ID
        :param precursor_rts: RTs of potential precursor spectra of different MS levels
        """
        ...
    
    @overload
    def getSpectrumMetaData(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData , flags: bytes ) -> None:
        """
        Cython signature: void getSpectrumMetaData(String spectrum_ref, SpectrumMetaData & meta, unsigned char flags)
        Extract meta data via a spectrum reference
        
        :param spectrum_ref: Spectrum reference to parse
        :param metadata: Meta data output
        :param flags: What meta data to extract
        """
        ...
    
    def setSpectraDataRef(self, spectra_data: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSpectraDataRef(const String & spectra_data)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Check if any spectra were set
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
    
    addMissingIMToPeptideIDs: __static_SpectrumMetaDataLookup_addMissingIMToPeptideIDs
    
    addMissingRTsToPeptideIDs: __static_SpectrumMetaDataLookup_addMissingRTsToPeptideIDs
    
    addMissingSpectrumReferences: __static_SpectrumMetaDataLookup_addMissingSpectrumReferences
    
    getSpectrumMetaData: __static_SpectrumMetaDataLookup_getSpectrumMetaData 


class StringDataArray:
    """
    Cython implementation of _StringDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1StringDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra string data attached to a spectrum or chromatogram.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void StringDataArray()
        """
        ...
    
    @overload
    def __init__(self, in_0: StringDataArray ) -> None:
        """
        Cython signature: void StringDataArray(StringDataArray &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def resize(self, n: int ) -> None:
        """
        Cython signature: void resize(size_t n)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def push_back(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void push_back(String)
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
    
    def __richcmp__(self, other: StringDataArray, op: int) -> Any:
        ... 


class __CHARGEMODE_MFD:
    None
    QFROMFEATURE : int
    QHEURISTIC : int
    QALL : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class LogType:
    None
    CMD : int
    GUI : int
    NONE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __PeptidePosition:
    None
    INTERNAL : int
    C_TERM : int
    N_TERM : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ProteinProteinCrossLinkType:
    None
    CROSS : int
    MONO : int
    LOOP : int
    NUMBER_OF_CROSS_LINK_TYPES : int

    def getMapping(self) -> Dict[int, str]:
       ... 

