from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_TransformationDescription_getModelTypes(result: List[bytes] ) -> None:
    """
    Cython signature: void getModelTypes(StringList result)
    """
    ...

def __static_SpectrumHelper_removePeaks(p: MSChromatogram , pos_start: float , pos_end: float ) -> None:
    """
    Cython signature: void removePeaks(MSChromatogram & p, double pos_start, double pos_end)
    """
    ...

def __static_SpectrumHelper_removePeaks(p: MSSpectrum , pos_start: float , pos_end: float ) -> None:
    """
    Cython signature: void removePeaks(MSSpectrum & p, double pos_start, double pos_end)
    """
    ...

def __static_SpectrumHelper_subtractMinimumIntensity(p: MSChromatogram ) -> None:
    """
    Cython signature: void subtractMinimumIntensity(MSChromatogram & p)
    """
    ...

def __static_SpectrumHelper_subtractMinimumIntensity(p: MSSpectrum ) -> None:
    """
    Cython signature: void subtractMinimumIntensity(MSSpectrum & p)
    """
    ...

def __static_IMTypes_toDriftTimeUnit(dtu_string: bytes ) -> int:
    """
    Cython signature: DriftTimeUnit toDriftTimeUnit(const libcpp_string & dtu_string)
    """
    ...

def __static_IMTypes_toIMFormat(IM_format: bytes ) -> int:
    """
    Cython signature: IMFormat toIMFormat(const libcpp_string & IM_format)
    """
    ...

def __static_IMTypes_toString(value: int ) -> bytes:
    """
    Cython signature: libcpp_string toString(const DriftTimeUnit value)
    """
    ...

def __static_IMTypes_toString(value: int ) -> bytes:
    """
    Cython signature: libcpp_string toString(const IMFormat value)
    """
    ...


class CVTerm_ControlledVocabulary:
    """
    Cython implementation of _CVTerm_ControlledVocabulary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTerm_ControlledVocabulary.html>`_
    """
    
    name: Union[bytes, str, String]
    
    id: Union[bytes, str, String]
    
    parents: Set[bytes]
    
    children: Set[bytes]
    
    obsolete: bool
    
    description: Union[bytes, str, String]
    
    synonyms: List[bytes]
    
    unparsed: List[bytes]
    
    xref_type: int
    
    xref_binary: List[bytes]
    
    units: Set[bytes]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVTerm_ControlledVocabulary()
        """
        ...
    
    @overload
    def __init__(self, rhs: CVTerm_ControlledVocabulary ) -> None:
        """
        Cython signature: void CVTerm_ControlledVocabulary(CVTerm_ControlledVocabulary rhs)
        """
        ...
    
    @overload
    def toXMLString(self, ref: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String toXMLString(String ref, String value)
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
        """
        ...
    
    @overload
    def toXMLString(self, ref: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[bytes, str, String]:
        """
        Cython signature: String toXMLString(String ref, DataValue value)
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
        """
        ...
    
    def getXRefTypeName(self, type: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getXRefTypeName(XRefType_CVTerm_ControlledVocabulary type)
        """
        ...
    
    def isHigherBetterScore(self, term: CVTerm_ControlledVocabulary ) -> bool:
        """
        Cython signature: bool isHigherBetterScore(CVTerm_ControlledVocabulary term)
        """
        ... 


class ControlledVocabulary:
    """
    Cython implementation of _ControlledVocabulary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ControlledVocabulary.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ControlledVocabulary()
        """
        ...
    
    @overload
    def __init__(self, in_0: ControlledVocabulary ) -> None:
        """
        Cython signature: void ControlledVocabulary(ControlledVocabulary &)
        """
        ...
    
    def name(self) -> Union[bytes, str, String]:
        """
        Cython signature: String name()
        Returns the CV name (set in the load method)
        """
        ...
    
    def loadFromOBO(self, name: Union[bytes, str, String] , filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void loadFromOBO(String name, String filename)
        Loads the CV from an OBO file
        """
        ...
    
    def exists(self, id: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool exists(String id)
        Returns true if the term is in the CV. Returns false otherwise.
        """
        ...
    
    def hasTermWithName(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasTermWithName(String name)
        Returns true if a term with the given name is in the CV. Returns false otherwise
        """
        ...
    
    def getTerm(self, id: Union[bytes, str, String] ) -> CVTerm_ControlledVocabulary:
        """
        Cython signature: CVTerm_ControlledVocabulary getTerm(String id)
        Returns a term specified by ID
        """
        ...
    
    def getTermByName(self, name: Union[bytes, str, String] , desc: Union[bytes, str, String] ) -> CVTerm_ControlledVocabulary:
        """
        Cython signature: CVTerm_ControlledVocabulary getTermByName(String name, String desc)
        Returns a term specified by name
        """
        ...
    
    def getAllChildTerms(self, terms: Set[bytes] , parent: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void getAllChildTerms(libcpp_set[String] terms, String parent)
        Writes all child terms recursively into terms
        """
        ...
    
    def isChildOf(self, child: Union[bytes, str, String] , parent: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool isChildOf(String child, String parent)
        Returns True if `child` is a child of `parent`
        """
        ... 


class FloatDataArray:
    """
    Cython implementation of _FloatDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1FloatDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra float data attached to a spectrum or chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FloatDataArray()
        """
        ...
    
    @overload
    def __init__(self, in_0: FloatDataArray ) -> None:
        """
        Cython signature: void FloatDataArray(FloatDataArray &)
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
    
    def reserve(self, n: int ) -> None:
        """
        Cython signature: void reserve(size_t n)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def push_back(self, in_0: float ) -> None:
        """
        Cython signature: void push_back(float)
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
    
    def __richcmp__(self, other: FloatDataArray, op: int) -> Any:
        ... 


class IDRipper:
    """
    Cython implementation of _IDRipper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1IDRipper.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void IDRipper()
        Ripping protein/peptide identification according their file origin
        """
        ...
    
    def rip(self, rfis: List[RipFileIdentifier] , rfcs: List[RipFileContent] , proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList , full_split: bool , split_ident_runs: bool ) -> None:
        """
        Cython signature: void rip(libcpp_vector[RipFileIdentifier] & rfis, libcpp_vector[RipFileContent] & rfcs, libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList & peptides, bool full_split, bool split_ident_runs)
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


class IMTypes:
    """
    Cython implementation of _IMTypes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IMTypes.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IMTypes()
        """
        ...
    
    @overload
    def __init__(self, in_0: IMTypes ) -> None:
        """
        Cython signature: void IMTypes(IMTypes &)
        """
        ...
    
    @overload
    def determineIMFormat(self, exp: MSExperiment ) -> int:
        """
        Cython signature: IMFormat determineIMFormat(const MSExperiment & exp)
        """
        ...
    
    @overload
    def determineIMFormat(self, spec: MSSpectrum ) -> int:
        """
        Cython signature: IMFormat determineIMFormat(const MSSpectrum & spec)
        """
        ...
    
    toDriftTimeUnit: __static_IMTypes_toDriftTimeUnit
    
    toIMFormat: __static_IMTypes_toIMFormat
    
    toString: __static_IMTypes_toString
    
    toString: __static_IMTypes_toString 


class IdentificationRuns:
    """
    Cython implementation of _IdentificationRuns

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1IdentificationRuns.html>`_
    """
    
    def __init__(self, prot_ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void IdentificationRuns(libcpp_vector[ProteinIdentification] & prot_ids)
        """
        ... 


class InstrumentSettings:
    """
    Cython implementation of _InstrumentSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InstrumentSettings.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void InstrumentSettings()
        Description of the settings a MS Instrument was run with
        """
        ...
    
    @overload
    def __init__(self, in_0: InstrumentSettings ) -> None:
        """
        Cython signature: void InstrumentSettings(InstrumentSettings &)
        """
        ...
    
    def getPolarity(self) -> int:
        """
        Cython signature: Polarity getPolarity()
        Returns the polarity
        """
        ...
    
    def setPolarity(self, in_0: int ) -> None:
        """
        Cython signature: void setPolarity(Polarity)
        Sets the polarity
        """
        ...
    
    def getScanMode(self) -> int:
        """
        Cython signature: ScanMode getScanMode()
        Returns the scan mode
        """
        ...
    
    def setScanMode(self, scan_mode: int ) -> None:
        """
        Cython signature: void setScanMode(ScanMode scan_mode)
        Sets the scan mode
        """
        ...
    
    def getZoomScan(self) -> bool:
        """
        Cython signature: bool getZoomScan()
        Returns if this scan is a zoom (enhanced resolution) scan
        """
        ...
    
    def setZoomScan(self, zoom_scan: bool ) -> None:
        """
        Cython signature: void setZoomScan(bool zoom_scan)
        Sets if this scan is a zoom (enhanced resolution) scan
        """
        ...
    
    def getScanWindows(self) -> List[ScanWindow]:
        """
        Cython signature: libcpp_vector[ScanWindow] getScanWindows()
        Returns the m/z scan windows
        """
        ...
    
    def setScanWindows(self, scan_windows: List[ScanWindow] ) -> None:
        """
        Cython signature: void setScanWindows(libcpp_vector[ScanWindow] scan_windows)
        Sets the m/z scan windows
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
    
    def __richcmp__(self, other: InstrumentSettings, op: int) -> Any:
        ... 


class MRMFP_ComponentGroupParams:
    """
    Cython implementation of _MRMFP_ComponentGroupParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFP_ComponentGroupParams.html>`_
    """
    
    component_group_name: Union[bytes, str, String]
    
    params: Param
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFP_ComponentGroupParams()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFP_ComponentGroupParams ) -> None:
        """
        Cython signature: void MRMFP_ComponentGroupParams(MRMFP_ComponentGroupParams &)
        """
        ... 


class MRMFP_ComponentParams:
    """
    Cython implementation of _MRMFP_ComponentParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFP_ComponentParams.html>`_
    """
    
    component_name: Union[bytes, str, String]
    
    component_group_name: Union[bytes, str, String]
    
    params: Param
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFP_ComponentParams()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFP_ComponentParams ) -> None:
        """
        Cython signature: void MRMFP_ComponentParams(MRMFP_ComponentParams &)
        """
        ... 


class MRMFeaturePicker:
    """
    Cython implementation of _MRMFeaturePicker

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeaturePicker.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMFeaturePicker()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMFeaturePicker ) -> None:
        """
        Cython signature: void MRMFeaturePicker(MRMFeaturePicker &)
        """
        ... 


class MSNumpressCoder:
    """
    Cython implementation of _MSNumpressCoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSNumpressCoder.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSNumpressCoder()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSNumpressCoder ) -> None:
        """
        Cython signature: void MSNumpressCoder(MSNumpressCoder &)
        """
        ...
    
    def encodeNP(self, in_: List[float] , result: String , zlib_compression: bool , config: NumpressConfig ) -> None:
        """
        Cython signature: void encodeNP(libcpp_vector[double] in_, String & result, bool zlib_compression, NumpressConfig config)
        Encodes a vector of floating point numbers into a Base64 string using numpress
        
        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
        
        This function will first apply the numpress encoding to the data, then
        encode the result in base64 (with optional zlib compression before
        base64 encoding)
        
        :note In case of error, result string is empty
        
        
        :param in: The vector of floating point numbers to be encoded
        :param result: The resulting string
        :param zlib_compression: Whether to apply zlib compression after numpress compression
        :param config: The numpress configuration defining the compression strategy
        """
        ...
    
    def decodeNP(self, in_: Union[bytes, str, String] , out: List[float] , zlib_compression: bool , config: NumpressConfig ) -> None:
        """
        Cython signature: void decodeNP(const String & in_, libcpp_vector[double] & out, bool zlib_compression, NumpressConfig config)
        Decodes a Base64 string to a vector of floating point numbers using numpress
        
        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
        
        This function will first decode the input base64 string (with optional
        zlib decompression after decoding) and then apply numpress decoding to
        the data
        
        
        :param in: The base64 encoded string
        :param out: The resulting vector of doubles
        :param zlib_compression: Whether to apply zlib de-compression before numpress de-compression
        :param config: The numpress configuration defining the compression strategy
        :raises:
          Exception: ConversionError if the string cannot be converted
        """
        ...
    
    def encodeNPRaw(self, in_: List[float] , result: String , config: NumpressConfig ) -> None:
        """
        Cython signature: void encodeNPRaw(libcpp_vector[double] in_, String & result, NumpressConfig config)
        Encode the data vector "in" to a raw byte array
        
        :note In case of error, "result" is given back unmodified
        :note The result is not a string but a raw byte array and may contain zero bytes
        
        This performs the raw numpress encoding on a set of data and does no
        Base64 encoding on the result. Therefore the result string is likely
        *unsafe* to handle and is a raw byte array.
        
        Please use the safe versions above unless you need access to the raw
        byte arrays
        
        
        :param in: The vector of floating point numbers to be encoded
        :param result: The resulting string
        :param config: The numpress configuration defining the compression strategy
        """
        ...
    
    def decodeNPRaw(self, in_: Union[bytes, str, String] , out: List[float] , config: NumpressConfig ) -> None:
        """
        Cython signature: void decodeNPRaw(const String & in_, libcpp_vector[double] & out, NumpressConfig config)
        Decode the raw byte array "in" to the result vector "out"
        
        :note The string in should *only* contain the data and _no_ extra
        null terminating byte
        
        This performs the raw numpress decoding on a raw byte array (not Base64
        encoded). Therefore the input string is likely *unsafe* to handle and is
        basically a byte container
        
        Please use the safe versions above unless you need access to the raw
        byte arrays
        
        
        :param in: The base64 encoded string
        :param out: The resulting vector of doubles
        :param config: The numpress configuration defining the compression strategy
        """
        ...
    NumpressCompression : __NumpressCompression 


class MSstatsFile:
    """
    Cython implementation of _MSstatsFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSstatsFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSstatsFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSstatsFile ) -> None:
        """
        Cython signature: void MSstatsFile(MSstatsFile &)
        """
        ...
    
    def storeLFQ(self, filename: String , consensus_map: ConsensusMap , design: ExperimentalDesign , reannotate_filenames: List[bytes] , is_isotope_label_type: bool , bioreplicate: String , condition: String , retention_time_summarization_method: String ) -> None:
        """
        Cython signature: void storeLFQ(String & filename, ConsensusMap & consensus_map, ExperimentalDesign & design, StringList & reannotate_filenames, bool is_isotope_label_type, String & bioreplicate, String & condition, String & retention_time_summarization_method)
        Store label free experiment (MSstats)
        """
        ...
    
    def storeISO(self, filename: String , consensus_map: ConsensusMap , design: ExperimentalDesign , reannotate_filenames: List[bytes] , bioreplicate: String , condition: String , mixture: String , retention_time_summarization_method: String ) -> None:
        """
        Cython signature: void storeISO(String & filename, ConsensusMap & consensus_map, ExperimentalDesign & design, StringList & reannotate_filenames, String & bioreplicate, String & condition, String & mixture, String & retention_time_summarization_method)
        Store isobaric experiment (MSstatsTMT)
        """
        ... 


class NumpressConfig:
    """
    Cython implementation of _NumpressConfig

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NumpressConfig.html>`_
    """
    
    numpressFixedPoint: float
    
    numpressErrorTolerance: float
    
    np_compression: int
    
    estimate_fixed_point: bool
    
    linear_fp_mass_acc: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NumpressConfig()
        """
        ...
    
    @overload
    def __init__(self, in_0: NumpressConfig ) -> None:
        """
        Cython signature: void NumpressConfig(NumpressConfig &)
        """
        ...
    
    def setCompression(self, compression: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCompression(const String & compression)
        """
        ... 


class OPXLSpectrumProcessingAlgorithms:
    """
    Cython implementation of _OPXLSpectrumProcessingAlgorithms

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLSpectrumProcessingAlgorithms.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OPXLSpectrumProcessingAlgorithms()
        """
        ...
    
    @overload
    def __init__(self, in_0: OPXLSpectrumProcessingAlgorithms ) -> None:
        """
        Cython signature: void OPXLSpectrumProcessingAlgorithms(OPXLSpectrumProcessingAlgorithms &)
        """
        ...
    
    def mergeAnnotatedSpectra(self, first_spectrum: MSSpectrum , second_spectrum: MSSpectrum ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum mergeAnnotatedSpectra(MSSpectrum & first_spectrum, MSSpectrum & second_spectrum)
        """
        ...
    
    def preprocessSpectra(self, exp: MSExperiment , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peptide_min_size: int , min_precursor_charge: int , max_precursor_charge: int , deisotope: bool , labeled: bool ) -> MSExperiment:
        """
        Cython signature: MSExperiment preprocessSpectra(MSExperiment & exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, size_t peptide_min_size, int min_precursor_charge, int max_precursor_charge, bool deisotope, bool labeled)
        """
        ...
    
    def getSpectrumAlignmentFastCharge(self, alignment: List[List[int, int]] , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , theo_spectrum: MSSpectrum , exp_spectrum: MSSpectrum , theo_charges: IntegerDataArray , exp_charges: IntegerDataArray , ppm_error_array: FloatDataArray , intensity_cutoff: float ) -> None:
        """
        Cython signature: void getSpectrumAlignmentFastCharge(libcpp_vector[libcpp_pair[size_t,size_t]] & alignment, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const MSSpectrum & theo_spectrum, const MSSpectrum & exp_spectrum, const IntegerDataArray & theo_charges, const IntegerDataArray & exp_charges, FloatDataArray & ppm_error_array, double intensity_cutoff)
        """
        ...
    
    def getSpectrumAlignmentSimple(self, alignment: List[List[int, int]] , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , theo_spectrum: List[SimplePeak] , exp_spectrum: MSSpectrum , exp_charges: IntegerDataArray ) -> None:
        """
        Cython signature: void getSpectrumAlignmentSimple(libcpp_vector[libcpp_pair[size_t,size_t]] & alignment, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const libcpp_vector[SimplePeak] & theo_spectrum, const MSSpectrum & exp_spectrum, const IntegerDataArray & exp_charges)
        """
        ... 


class OSBinaryDataArray:
    """
    Cython implementation of _OSBinaryDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSBinaryDataArray.html>`_
    """
    
    data: List[float]
    
    description: bytes
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OSBinaryDataArray()
        """
        ...
    
    @overload
    def __init__(self, in_0: OSBinaryDataArray ) -> None:
        """
        Cython signature: void OSBinaryDataArray(OSBinaryDataArray &)
        """
        ... 


class OSChromatogram:
    """
    Cython implementation of _OSChromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSChromatogram.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OSChromatogram()
        """
        ...
    
    @overload
    def __init__(self, in_0: OSChromatogram ) -> None:
        """
        Cython signature: void OSChromatogram(OSChromatogram &)
        """
        ... 


class OSSpectrum:
    """
    Cython implementation of _OSSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSSpectrum.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OSSpectrum()
        """
        ...
    
    @overload
    def __init__(self, in_0: OSSpectrum ) -> None:
        """
        Cython signature: void OSSpectrum(OSSpectrum &)
        """
        ... 


class RipFileContent:
    """
    Cython implementation of _RipFileContent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1RipFileContent.html>`_
    """
    
    def __init__(self, prot_idents: List[ProteinIdentification] , pep_idents: PeptideIdentificationList ) -> None:
        """
        Cython signature: void RipFileContent(libcpp_vector[ProteinIdentification] & prot_idents, PeptideIdentificationList & pep_idents)
        """
        ...
    
    def getProteinIdentifications(self) -> List[ProteinIdentification]:
        """
        Cython signature: libcpp_vector[ProteinIdentification] getProteinIdentifications()
        """
        ...
    
    def getPeptideIdentifications(self) -> PeptideIdentificationList:
        """
        Cython signature: PeptideIdentificationList getPeptideIdentifications()
        """
        ... 


class RipFileIdentifier:
    """
    Cython implementation of _RipFileIdentifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1RipFileIdentifier.html>`_
    """
    
    def __init__(self, id_runs: IdentificationRuns , pep_id: PeptideIdentification , file_origin_map: Dict[Union[bytes, str, String], int] , origin_annotation_fmt: int , split_ident_runs: bool ) -> None:
        """
        Cython signature: void RipFileIdentifier(IdentificationRuns & id_runs, PeptideIdentification & pep_id, libcpp_map[String,unsigned int] & file_origin_map, OriginAnnotationFormat origin_annotation_fmt, bool split_ident_runs)
        """
        ...
    
    def getIdentRunIdx(self) -> int:
        """
        Cython signature: unsigned int getIdentRunIdx()
        """
        ...
    
    def getFileOriginIdx(self) -> int:
        """
        Cython signature: unsigned int getFileOriginIdx()
        """
        ...
    
    def getOriginFullname(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOriginFullname()
        """
        ...
    
    def getOutputBasename(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOutputBasename()
        """
        ... 


class SpectrumAccessSqMass:
    """
    Cython implementation of _SpectrumAccessSqMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessSqMass.html>`_
      -- Inherits from ['ISpectrumAccess']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAccessSqMass()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessSqMass ) -> None:
        """
        Cython signature: void SpectrumAccessSqMass(SpectrumAccessSqMass &)
        """
        ...
    
    @overload
    def __init__(self, in_0: MzMLSqliteHandler , indices: List[int] ) -> None:
        """
        Cython signature: void SpectrumAccessSqMass(MzMLSqliteHandler, libcpp_vector[int] indices)
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


class SpectrumHelper:
    """
    Cython implementation of _SpectrumHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/class_1_1SpectrumHelper.html>`_
    """
    
    removePeaks: __static_SpectrumHelper_removePeaks
    
    removePeaks: __static_SpectrumHelper_removePeaks
    
    subtractMinimumIntensity: __static_SpectrumHelper_subtractMinimumIntensity
    
    subtractMinimumIntensity: __static_SpectrumHelper_subtractMinimumIntensity 


class SwathMapMassCorrection:
    """
    Cython implementation of _SwathMapMassCorrection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathMapMassCorrection.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SwathMapMassCorrection()
        """
        ...
    
    @overload
    def __init__(self, in_0: SwathMapMassCorrection ) -> None:
        """
        Cython signature: void SwathMapMassCorrection(SwathMapMassCorrection)
        """
        ... 


class TransformationDescription:
    """
    Cython implementation of _TransformationDescription

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TransformationDescription_1_1TransformationDescription.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TransformationDescription()
        """
        ...
    
    @overload
    def __init__(self, in_0: TransformationDescription ) -> None:
        """
        Cython signature: void TransformationDescription(TransformationDescription &)
        """
        ...
    
    def getDataPoints(self) -> List[TM_DataPoint]:
        """
        Cython signature: libcpp_vector[TM_DataPoint] getDataPoints()
        Returns the data points
        """
        ...
    
    @overload
    def setDataPoints(self, data: List[TM_DataPoint] ) -> None:
        """
        Cython signature: void setDataPoints(libcpp_vector[TM_DataPoint] & data)
        Sets the data points. Removes the model that was previously fitted to the data (if any)
        """
        ...
    
    @overload
    def setDataPoints(self, data: List[List[float, float]] ) -> None:
        """
        Cython signature: void setDataPoints(libcpp_vector[libcpp_pair[double,double]] & data)
        Sets the data points (backwards-compatible overload). Removes the model that was previously fitted to the data (if any)
        """
        ...
    
    def apply(self, in_0: float ) -> float:
        """
        Cython signature: double apply(double)
        Applies the transformation to `value`
        """
        ...
    
    @overload
    def fitModel(self, model_type: Union[bytes, str, String] , params: Param ) -> None:
        """
        Cython signature: void fitModel(String model_type, Param params)
        Fits a model to the data
        """
        ...
    
    @overload
    def fitModel(self, model_type: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void fitModel(String model_type)
        Fits a model to the data
        """
        ...
    
    def getModelType(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getModelType()
        Gets the type of the fitted model
        """
        ...
    
    def getModelParameters(self) -> Param:
        """
        Cython signature: Param getModelParameters()
        Returns the model parameters
        """
        ...
    
    def invert(self) -> None:
        """
        Cython signature: void invert()
        Computes an (approximate) inverse of the transformation
        """
        ...
    
    def getDeviations(self, diffs: List[float] , do_apply: bool , do_sort: bool ) -> None:
        """
        Cython signature: void getDeviations(libcpp_vector[double] & diffs, bool do_apply, bool do_sort)
        Get the deviations between the data pairs
        
        :param diffs: Output
        :param do_apply: Get deviations after applying the model?
        :param do_sort: Sort `diffs` before returning?
        """
        ...
    
    def getStatistics(self) -> TransformationStatistics:
        """
        Cython signature: TransformationStatistics getStatistics()
        """
        ...
    
    getModelTypes: __static_TransformationDescription_getModelTypes 


class TransformationStatistics:
    """
    Cython implementation of _TransformationStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TransformationDescription_1_1TransformationStatistics.html>`_
    """
    
    xmin: float
    
    xmax: float
    
    ymin: float
    
    ymax: float
    
    percentiles_before: Dict[int, float]
    
    percentiles_after: Dict[int, float]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TransformationStatistics()
        """
        ...
    
    @overload
    def __init__(self, in_0: TransformationStatistics ) -> None:
        """
        Cython signature: void TransformationStatistics(TransformationStatistics &)
        """
        ... 


class DriftTimeUnit:
    None
    NONE : int
    MILLISECOND : int
    VSSC : int
    FAIMS_COMPENSATION_VOLTAGE : int
    SIZE_OF_DRIFTTIMEUNIT : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class IMFormat:
    None
    NONE : int
    CONCATENATED : int
    MULTIPLE_SPECTRA : int
    MIXED : int
    SIZE_OF_IMFORMAT : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __NumpressCompression:
    None
    NONE : int
    LINEAR : int
    PIC : int
    SLOF : int
    SIZE_OF_NUMPRESSCOMPRESSION : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class OriginAnnotationFormat:
    None
    FILE_ORIGIN : int
    MAP_INDEX : int
    ID_MERGE_INDEX : int
    UNKNOWN_OAF : int
    SIZE_OF_ORIGIN_ANNOTATION_FORMAT : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class ScanMode:
    None
    UNKNOWN : int
    MASSSPECTRUM : int
    MS1SPECTRUM : int
    MSNSPECTRUM : int
    SIM : int
    SRM : int
    CRM : int
    CNG : int
    CNL : int
    PRECURSOR : int
    EMC : int
    TDF : int
    EMR : int
    EMISSION : int
    ABSORPTION : int
    SIZE_OF_SCANMODE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class XRefType_CVTerm_ControlledVocabulary:
    None
    XSD_STRING : int
    XSD_INTEGER : int
    XSD_DECIMAL : int
    XSD_NEGATIVE_INTEGER : int
    XSD_POSITIVE_INTEGER : int
    XSD_NON_NEGATIVE_INTEGER : int
    XSD_NON_POSITIVE_INTEGER : int
    XSD_BOOLEAN : int
    XSD_DATE : int
    XSD_ANYURI : int
    NONE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

