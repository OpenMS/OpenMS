from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_CalibrationData_getMetaValues() -> List[bytes]:
    """
    Cython signature: StringList getMetaValues()
    """
    ...


class AbsoluteQuantitationMethodFile:
    """
    Cython implementation of _AbsoluteQuantitationMethodFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationMethodFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AbsoluteQuantitationMethodFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: AbsoluteQuantitationMethodFile ) -> None:
        """
        Cython signature: void AbsoluteQuantitationMethodFile(AbsoluteQuantitationMethodFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , aqm_list: List[AbsoluteQuantitationMethod] ) -> None:
        """
        Cython signature: void load(const String & filename, libcpp_vector[AbsoluteQuantitationMethod] & aqm_list)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , aqm_list: List[AbsoluteQuantitationMethod] ) -> None:
        """
        Cython signature: void store(const String & filename, libcpp_vector[AbsoluteQuantitationMethod] & aqm_list)
        """
        ... 


class CalibrationData:
    """
    Cython implementation of _CalibrationData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CalibrationData.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CalibrationData()
        """
        ...
    
    @overload
    def __init__(self, in_0: CalibrationData ) -> None:
        """
        Cython signature: void CalibrationData(CalibrationData &)
        """
        ...
    
    def getMZ(self, in_0: int ) -> float:
        """
        Cython signature: double getMZ(size_t)
        Retrieve the observed m/z of the i'th calibration point
        """
        ...
    
    def getRT(self, in_0: int ) -> float:
        """
        Cython signature: double getRT(size_t)
        Retrieve the observed RT of the i'th calibration point
        """
        ...
    
    def getIntensity(self, in_0: int ) -> float:
        """
        Cython signature: double getIntensity(size_t)
        Retrieve the intensity of the i'th calibration point
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Number of calibration points
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Returns `True` if there are no peaks
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Remove all calibration points
        """
        ...
    
    def setUsePPM(self, in_0: bool ) -> None:
        """
        Cython signature: void setUsePPM(bool)
        """
        ...
    
    def usePPM(self) -> bool:
        """
        Cython signature: bool usePPM()
        Current error unit (ppm or Th)
        """
        ...
    
    def insertCalibrationPoint(self, rt: float , mz_obs: float , intensity: float , mz_ref: float , weight: float , group: int ) -> None:
        """
        Cython signature: void insertCalibrationPoint(double rt, double mz_obs, float intensity, double mz_ref, double weight, int group)
        """
        ...
    
    def getNrOfGroups(self) -> int:
        """
        Cython signature: size_t getNrOfGroups()
        Number of peak groups (can be 0)
        """
        ...
    
    def getError(self, in_0: int ) -> float:
        """
        Cython signature: double getError(size_t)
        Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())
        """
        ...
    
    def getRefMZ(self, in_0: int ) -> float:
        """
        Cython signature: double getRefMZ(size_t)
        Retrieve the theoretical m/z of the i'th calibration point
        """
        ...
    
    def getWeight(self, in_0: int ) -> float:
        """
        Cython signature: double getWeight(size_t)
        Retrieve the weight of the i'th calibration point
        """
        ...
    
    def getGroup(self, i: int ) -> int:
        """
        Cython signature: int getGroup(size_t i)
        Retrieve the group of the i'th calibration point
        """
        ...
    
    def median(self, in_0: float , in_1: float ) -> CalibrationData:
        """
        Cython signature: CalibrationData median(double, double)
        Compute the median in the given RT range for every peak group
        """
        ...
    
    def sortByRT(self) -> None:
        """
        Cython signature: void sortByRT()
        Sort calibration points by RT, to allow for valid RT chunking
        """
        ...
    
    getMetaValues: __static_CalibrationData_getMetaValues 


class DigestionEnzyme:
    """
    Cython implementation of _DigestionEnzyme

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzyme.html>`_

      Base class for digestion enzymes
    """
    
    @overload
    def __init__(self, in_0: DigestionEnzyme ) -> None:
        """
        Cython signature: void DigestionEnzyme(DigestionEnzyme &)
        """
        ...
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void DigestionEnzyme(const String & name, const String & cleavage_regex, libcpp_set[String] & synonyms, String regex_description)
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String & name)
        Sets the name of the enzyme
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the enzyme
        """
        ...
    
    def setSynonyms(self, synonyms: Set[bytes] ) -> None:
        """
        Cython signature: void setSynonyms(libcpp_set[String] & synonyms)
        Sets the synonyms
        """
        ...
    
    def addSynonym(self, synonym: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addSynonym(const String & synonym)
        Adds a synonym
        """
        ...
    
    def getSynonyms(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getSynonyms()
        Returns the synonyms
        """
        ...
    
    def setRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setRegEx(const String & cleavage_regex)
        Sets the cleavage regex
        """
        ...
    
    def getRegEx(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getRegEx()
        Returns the cleavage regex
        """
        ...
    
    def setRegExDescription(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setRegExDescription(const String & value)
        Sets the regex description
        """
        ...
    
    def getRegExDescription(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getRegExDescription()
        Returns the regex description
        """
        ...
    
    def setValueFromFile(self, key: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool setValueFromFile(String key, String value)
        Sets the value of a member variable based on an entry from an input file
        """
        ...
    
    def __richcmp__(self, other: DigestionEnzyme, op: int) -> Any:
        ... 


class ElutionPeakDetection:
    """
    Cython implementation of _ElutionPeakDetection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElutionPeakDetection.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ElutionPeakDetection()
        """
        ...
    
    @overload
    def __init__(self, in_0: ElutionPeakDetection ) -> None:
        """
        Cython signature: void ElutionPeakDetection(ElutionPeakDetection &)
        """
        ...
    
    @overload
    def detectPeaks(self, in_: Kernel_MassTrace , out: List[Kernel_MassTrace] ) -> None:
        """
        Cython signature: void detectPeaks(Kernel_MassTrace & in_, libcpp_vector[Kernel_MassTrace] & out)
        """
        ...
    
    @overload
    def detectPeaks(self, in_: List[Kernel_MassTrace] , out: List[Kernel_MassTrace] ) -> None:
        """
        Cython signature: void detectPeaks(libcpp_vector[Kernel_MassTrace] & in_, libcpp_vector[Kernel_MassTrace] & out)
        """
        ...
    
    def filterByPeakWidth(self, in_: List[Kernel_MassTrace] , out: List[Kernel_MassTrace] ) -> None:
        """
        Cython signature: void filterByPeakWidth(libcpp_vector[Kernel_MassTrace] & in_, libcpp_vector[Kernel_MassTrace] & out)
        """
        ...
    
    def computeMassTraceNoise(self, in_0: Kernel_MassTrace ) -> float:
        """
        Cython signature: double computeMassTraceNoise(Kernel_MassTrace &)
        Compute noise level (as RMSE of the actual signal and the smoothed signal)
        """
        ...
    
    def computeMassTraceSNR(self, in_0: Kernel_MassTrace ) -> float:
        """
        Cython signature: double computeMassTraceSNR(Kernel_MassTrace &)
        Compute the signal to noise ratio (estimated by computeMassTraceNoise)
        """
        ...
    
    def computeApexSNR(self, in_0: Kernel_MassTrace ) -> float:
        """
        Cython signature: double computeApexSNR(Kernel_MassTrace &)
        Compute the signal to noise ratio at the apex (estimated by computeMassTraceNoise)
        """
        ...
    
    def findLocalExtrema(self, in_0: Kernel_MassTrace , in_1: int , in_2: List[int] , in_3: List[int] ) -> None:
        """
        Cython signature: void findLocalExtrema(Kernel_MassTrace &, size_t &, libcpp_vector[size_t] &, libcpp_vector[size_t] &)
        """
        ...
    
    def smoothData(self, mt: Kernel_MassTrace , win_size: int ) -> None:
        """
        Cython signature: void smoothData(Kernel_MassTrace & mt, int win_size)
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


class FeatureXMLFile:
    """
    Cython implementation of _FeatureXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FeatureXMLFile()
        This class provides Input/Output functionality for feature maps
        """
        ...
    
    def load(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None:
        """
        Cython signature: void load(String, FeatureMap &)
        Loads the file with name `filename` into `map` and calls updateRanges()
        """
        ...
    
    def store(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None:
        """
        Cython signature: void store(String, FeatureMap &)
        Stores the map `feature_map` in file with name `filename`
        """
        ...
    
    def getOptions(self) -> FeatureFileOptions:
        """
        Cython signature: FeatureFileOptions getOptions()
        Access to the options for loading/storing
        """
        ...
    
    def setOptions(self, in_0: FeatureFileOptions ) -> None:
        """
        Cython signature: void setOptions(FeatureFileOptions)
        Setter for options for loading/storing
        """
        ...
    
    def loadSize(self, path: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t loadSize(String path)
        """
        ... 


class GNPSMetaValueFile:
    """
    Cython implementation of _GNPSMetaValueFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSMetaValueFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GNPSMetaValueFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: GNPSMetaValueFile ) -> None:
        """
        Cython signature: void GNPSMetaValueFile(GNPSMetaValueFile &)
        """
        ...
    
    def store(self, consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const ConsensusMap & consensus_map, const String & output_file)
        Write meta value table (tsv file) from a list of mzML files. Required for GNPS FBMN.
        
        This will produce the minimal required meta values and can be extended manually.
        
        :param consensus_map: Input ConsensusMap from which the input mzML files will be determined.
        :param output_file: Output file path for the meta value table.
        """
        ... 


class HyperScore:
    """
    Cython implementation of _HyperScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1HyperScore.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void HyperScore()
        An implementation of the X!Tandem HyperScore PSM scoring function
        """
        ...
    
    @overload
    def __init__(self, in_0: HyperScore ) -> None:
        """
        Cython signature: void HyperScore(HyperScore &)
        """
        ...
    
    def compute(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , theo_spectrum: MSSpectrum ) -> float:
        """
        Cython signature: double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, MSSpectrum & exp_spectrum, MSSpectrum & theo_spectrum)
        Compute the (ln transformed) X!Tandem HyperScore\n
        
        1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
        2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
        
        
        :note: Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
        :param fragment_mass_tolerance: Mass tolerance applied left and right of the theoretical spectrum peak position
        :param fragment_mass_tolerance_unit_ppm: Unit of the mass tolerance is: Thomson if false, ppm if true
        :param exp_spectrum: Measured spectrum
        :param theo_spectrum: Theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator
        """
        ... 


class IMSWeights:
    """
    Cython implementation of _IMSWeights

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::Weights_1_1IMSWeights.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IMSWeights()
        """
        ...
    
    @overload
    def __init__(self, in_0: IMSWeights ) -> None:
        """
        Cython signature: void IMSWeights(IMSWeights)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        Gets size of a set of weights
        """
        ...
    
    def getWeight(self, i: int ) -> int:
        """
        Cython signature: unsigned int getWeight(int i)
        Gets a scaled integer weight by index
        """
        ...
    
    def setPrecision(self, precision: float ) -> None:
        """
        Cython signature: void setPrecision(double precision)
        Sets a new precision to scale double values to integer
        """
        ...
    
    def getPrecision(self) -> float:
        """
        Cython signature: double getPrecision()
        Gets precision.
        """
        ...
    
    def back(self) -> int:
        """
        Cython signature: unsigned int back()
        Gets a last weight
        """
        ...
    
    def getAlphabetMass(self, i: int ) -> float:
        """
        Cython signature: double getAlphabetMass(int i)
        Gets an original (double) alphabet mass by index
        """
        ...
    
    def getParentMass(self, decomposition: List[int] ) -> float:
        """
        Cython signature: double getParentMass(libcpp_vector[unsigned int] & decomposition)
        Returns a parent mass for a given `decomposition`
        """
        ...
    
    def swap(self, index1: int , index2: int ) -> None:
        """
        Cython signature: void swap(int index1, int index2)
        Exchanges weight and mass at index1 with weight and mass at index2
        """
        ...
    
    def divideByGCD(self) -> bool:
        """
        Cython signature: bool divideByGCD()
        Divides the integer weights by their gcd. The precision is also adjusted
        """
        ...
    
    def getMinRoundingError(self) -> float:
        """
        Cython signature: double getMinRoundingError()
        """
        ...
    
    def getMaxRoundingError(self) -> float:
        """
        Cython signature: double getMaxRoundingError()
        """
        ... 


class IsobaricNormalizer:
    """
    Cython implementation of _IsobaricNormalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricNormalizer.html>`_
    """
    
    @overload
    def __init__(self, in_0: IsobaricNormalizer ) -> None:
        """
        Cython signature: void IsobaricNormalizer(IsobaricNormalizer &)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricNormalizer(ItraqFourPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricNormalizer(ItraqEightPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTSixPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricNormalizer(TMTSixPlexQuantitationMethod * quant_method)
        """
        ...
    
    @overload
    def __init__(self, quant_method: TMTTenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void IsobaricNormalizer(TMTTenPlexQuantitationMethod * quant_method)
        """
        ...
    
    def normalize(self, consensus_map: ConsensusMap ) -> None:
        """
        Cython signature: void normalize(ConsensusMap & consensus_map)
        """
        ... 


class JavaInfo:
    """
    Cython implementation of _JavaInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1JavaInfo.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void JavaInfo()
        Detect Java and retrieve information
        """
        ...
    
    @overload
    def __init__(self, in_0: JavaInfo ) -> None:
        """
        Cython signature: void JavaInfo(JavaInfo &)
        """
        ...
    
    def canRun(self, java_executable: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool canRun(String java_executable)
        Determine if Java is installed and reachable\n
        
        The call fails if either Java is not installed or if a relative location is given and Java is not on the search PATH
        
        
        :param java_executable: Path to Java executable. Can be absolute, relative or just a filename
        :return: Returns false if Java executable can not be called; true if Java executable can be executed
        """
        ... 


class LPWrapper:
    """
    Cython implementation of _LPWrapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LPWrapper.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void LPWrapper()
        """
        ...
    
    @overload
    def addRow(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: int addRow(libcpp_vector[int] row_indices, libcpp_vector[double] row_values, const String & name)
        Adds a row to the LP matrix, returns index
        """
        ...
    
    @overload
    def addRow(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int:
        """
        Cython signature: int addRow(libcpp_vector[int] & row_indices, libcpp_vector[double] & row_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_)
        Adds a row with boundaries to the LP matrix, returns index
        """
        ...
    
    @overload
    def addColumn(self, ) -> int:
        """
        Cython signature: int addColumn()
        Adds an empty column to the LP matrix, returns index
        """
        ...
    
    @overload
    def addColumn(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: int addColumn(libcpp_vector[int] column_indices, libcpp_vector[double] column_values, const String & name)
        Adds a column to the LP matrix, returns index
        """
        ...
    
    @overload
    def addColumn(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int:
        """
        Cython signature: int addColumn(libcpp_vector[int] & column_indices, libcpp_vector[double] & column_values, const String & name, double lower_bound, double upper_bound, LPWrapper_Type type_)
        Adds a column with boundaries to the LP matrix, returns index
        """
        ...
    
    def deleteRow(self, index: int ) -> None:
        """
        Cython signature: void deleteRow(int index)
        Delete index-th row
        """
        ...
    
    def setColumnName(self, index: int , name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setColumnName(int index, const String & name)
        Sets name of the index-th column
        """
        ...
    
    def getColumnName(self, index: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getColumnName(int index)
        Returns name of the index-th column
        """
        ...
    
    def getRowName(self, index: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getRowName(int index)
        Sets name of the index-th row
        """
        ...
    
    def getRowIndex(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: int getRowIndex(const String & name)
        Returns index of the row with name
        """
        ...
    
    def getColumnIndex(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: int getColumnIndex(const String & name)
        Returns index of the column with name
        """
        ...
    
    def getColumnUpperBound(self, index: int ) -> float:
        """
        Cython signature: double getColumnUpperBound(int index)
        Returns column's upper bound
        """
        ...
    
    def getColumnLowerBound(self, index: int ) -> float:
        """
        Cython signature: double getColumnLowerBound(int index)
        Returns column's lower bound
        """
        ...
    
    def getRowUpperBound(self, index: int ) -> float:
        """
        Cython signature: double getRowUpperBound(int index)
        Returns row's upper bound
        """
        ...
    
    def getRowLowerBound(self, index: int ) -> float:
        """
        Cython signature: double getRowLowerBound(int index)
        Returns row's lower bound
        """
        ...
    
    def setRowName(self, index: int , name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setRowName(int index, const String & name)
        Sets name of the index-th row
        """
        ...
    
    def setColumnBounds(self, index: int , lower_bound: float , upper_bound: float , type_: int ) -> None:
        """
        Cython signature: void setColumnBounds(int index, double lower_bound, double upper_bound, LPWrapper_Type type_)
        Sets column bounds
        """
        ...
    
    def setRowBounds(self, index: int , lower_bound: float , upper_bound: float , type_: int ) -> None:
        """
        Cython signature: void setRowBounds(int index, double lower_bound, double upper_bound, LPWrapper_Type type_)
        Sets row bounds
        """
        ...
    
    def setColumnType(self, index: int , type_: int ) -> None:
        """
        Cython signature: void setColumnType(int index, VariableType type_)
        Sets column/variable type.
        """
        ...
    
    def getColumnType(self, index: int ) -> int:
        """
        Cython signature: VariableType getColumnType(int index)
        Returns column/variable type.
        """
        ...
    
    def setObjective(self, index: int , obj_value: float ) -> None:
        """
        Cython signature: void setObjective(int index, double obj_value)
        Sets objective value for column with index
        """
        ...
    
    def getObjective(self, index: int ) -> float:
        """
        Cython signature: double getObjective(int index)
        Returns objective value for column with index
        """
        ...
    
    def setObjectiveSense(self, sense: int ) -> None:
        """
        Cython signature: void setObjectiveSense(Sense sense)
        Sets objective direction
        """
        ...
    
    def getObjectiveSense(self) -> int:
        """
        Cython signature: Sense getObjectiveSense()
        Returns objective sense
        """
        ...
    
    def getNumberOfColumns(self) -> int:
        """
        Cython signature: int getNumberOfColumns()
        Returns number of columns
        """
        ...
    
    def getNumberOfRows(self) -> int:
        """
        Cython signature: int getNumberOfRows()
        Returns number of rows
        """
        ...
    
    def setElement(self, row_index: int , column_index: int , value: float ) -> None:
        """
        Cython signature: void setElement(int row_index, int column_index, double value)
        Sets the element
        """
        ...
    
    def getElement(self, row_index: int , column_index: int ) -> float:
        """
        Cython signature: double getElement(int row_index, int column_index)
        Returns the element
        """
        ...
    
    def readProblem(self, filename: Union[bytes, str, String] , format_: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readProblem(String filename, String format_)
        Read LP from file
        
        
        :param filename: Filename where to store the LP problem
        :param format: LP, MPS or GLPK
        """
        ...
    
    def writeProblem(self, filename: Union[bytes, str, String] , format_: int ) -> None:
        """
        Cython signature: void writeProblem(const String & filename, WriteFormat format_)
        Write LP formulation to a file
        
        
        :param filename: Output filename, if the filename ends with '.gz' it will be compressed
        :param format: MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK
        """
        ...
    
    def solve(self, solver_param: SolverParam , verbose_level: int ) -> int:
        """
        Cython signature: int solve(SolverParam & solver_param, size_t verbose_level)
        Solve problems, parameters like enabled heuristics can be given via solver_param\n
        
        The verbose level (0,1,2) determines if the solver prints status messages and internals
        
        
        :param solver_param: Parameters of the solver introduced by SolverParam
        :param verbose_level: Sets verbose level
        :return: solver dependent
        """
        ...
    
    def getStatus(self) -> int:
        """
        Cython signature: SolverStatus getStatus()
        Returns solution status
        
        
        :return: status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution
        """
        ...
    
    def getObjectiveValue(self) -> float:
        """
        Cython signature: double getObjectiveValue()
        """
        ...
    
    def getColumnValue(self, index: int ) -> float:
        """
        Cython signature: double getColumnValue(int index)
        """
        ...
    
    def getNumberOfNonZeroEntriesInRow(self, idx: int ) -> int:
        """
        Cython signature: int getNumberOfNonZeroEntriesInRow(int idx)
        """
        ...
    
    def getMatrixRow(self, idx: int , indexes: List[int] ) -> None:
        """
        Cython signature: void getMatrixRow(int idx, libcpp_vector[int] & indexes)
        """
        ...
    
    def getSolver(self) -> int:
        """
        Cython signature: SOLVER getSolver()
        Returns currently active solver
        """
        ...
    LPWrapper_Type : __LPWrapper_Type
    SOLVER : __SOLVER
    Sense : __Sense
    SolverStatus : __SolverStatus
    VariableType : __VariableType
    WriteFormat : __WriteFormat 


class OpenSwathDataAccessHelper:
    """
    Cython implementation of _OpenSwathDataAccessHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathDataAccessHelper.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OpenSwathDataAccessHelper()
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenSwathDataAccessHelper ) -> None:
        """
        Cython signature: void OpenSwathDataAccessHelper(OpenSwathDataAccessHelper &)
        """
        ...
    
    def convertToOpenMSSpectrum(self, sptr: OSSpectrum , spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void convertToOpenMSSpectrum(shared_ptr[OSSpectrum] sptr, MSSpectrum & spectrum)
        Converts a SpectrumPtr to an OpenMS Spectrum
        """
        ...
    
    def convertToOpenMSChromatogram(self, cptr: OSChromatogram , chromatogram: MSChromatogram ) -> None:
        """
        Cython signature: void convertToOpenMSChromatogram(shared_ptr[OSChromatogram] cptr, MSChromatogram & chromatogram)
        Converts a ChromatogramPtr to an OpenMS Chromatogram
        """
        ...
    
    def convertToOpenMSChromatogramFilter(self, chromatogram: MSChromatogram , cptr: OSChromatogram , rt_min: float , rt_max: float ) -> None:
        """
        Cython signature: void convertToOpenMSChromatogramFilter(MSChromatogram & chromatogram, shared_ptr[OSChromatogram] cptr, double rt_min, double rt_max)
        """
        ...
    
    def convertTargetedExp(self, transition_exp_: TargetedExperiment , transition_exp: LightTargetedExperiment ) -> None:
        """
        Cython signature: void convertTargetedExp(TargetedExperiment & transition_exp_, LightTargetedExperiment & transition_exp)
        Converts from the OpenMS TargetedExperiment to the OpenMs LightTargetedExperiment
        """
        ...
    
    def convertPeptideToAASequence(self, peptide: LightCompound , aa_sequence: AASequence ) -> None:
        """
        Cython signature: void convertPeptideToAASequence(LightCompound & peptide, AASequence & aa_sequence)
        Converts from the LightCompound to an OpenMS AASequence (with correct modifications)
        """
        ...
    
    def convertTargetedPeptide(self, pep: Peptide , p: LightCompound ) -> None:
        """
        Cython signature: void convertTargetedPeptide(Peptide pep, LightCompound & p)
        Converts from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide
        """
        ...
    
    def convertTargetedNuctide(self, pep: Nuctide , p: LightCompound ) -> None:
        """
        Cython signature: void convertTargetedNuctide(Nuctide pep, LightCompound & p)
        Converts from the OpenMS TargetedExperiment Nuctide to the LightTargetedExperiment Nuctide
        """
        ... 


class ParamValue:
    """
    Cython implementation of _ParamValue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamValue.html>`_

    Class to hold strings, numeric values, vectors of strings and vectors of numeric values using the stl types
    
    - To choose one of these types, just use the appropriate constructor
    - Automatic conversion is supported and throws Exceptions in case of invalid conversions
    - An empty object is created with the default constructor
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ParamValue()
        """
        ...
    
    @overload
    def __init__(self, in_0: ParamValue ) -> None:
        """
        Cython signature: void ParamValue(ParamValue &)
        """
        ...
    
    @overload
    def __init__(self, in_0: bytes ) -> None:
        """
        Cython signature: void ParamValue(char *)
        """
        ...
    
    @overload
    def __init__(self, in_0: Union[bytes, str] ) -> None:
        """
        Cython signature: void ParamValue(const libcpp_utf8_string &)
        """
        ...
    
    @overload
    def __init__(self, in_0: int ) -> None:
        """
        Cython signature: void ParamValue(int)
        """
        ...
    
    @overload
    def __init__(self, in_0: float ) -> None:
        """
        Cython signature: void ParamValue(double)
        """
        ...
    
    @overload
    def __init__(self, in_0: List[Union[bytes, str]] ) -> None:
        """
        Cython signature: void ParamValue(libcpp_vector[libcpp_utf8_string])
        """
        ...
    
    @overload
    def __init__(self, in_0: List[int] ) -> None:
        """
        Cython signature: void ParamValue(libcpp_vector[int])
        """
        ...
    
    @overload
    def __init__(self, in_0: List[float] ) -> None:
        """
        Cython signature: void ParamValue(libcpp_vector[double])
        """
        ...
    
    def toStringVector(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[libcpp_string] toStringVector()
        Explicitly convert ParamValue to string vector
        """
        ...
    
    def toDoubleVector(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] toDoubleVector()
        Explicitly convert ParamValue to DoubleList
        """
        ...
    
    def toIntVector(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] toIntVector()
        Explicitly convert ParamValue to IntList
        """
        ...
    
    def toBool(self) -> bool:
        """
        Cython signature: bool toBool()
        Converts the strings 'true' and 'false' to a bool
        """
        ...
    
    def valueType(self) -> int:
        """
        Cython signature: ValueType valueType()
        """
        ...
    
    def isEmpty(self) -> int:
        """
        Cython signature: int isEmpty()
        Test if the value is empty
        """
        ... 


class SolverParam:
    """
    Cython implementation of _SolverParam

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SolverParam.html>`_
    """
    
    message_level: int
    
    branching_tech: int
    
    backtrack_tech: int
    
    preprocessing_tech: int
    
    enable_feas_pump_heuristic: bool
    
    enable_gmi_cuts: bool
    
    enable_mir_cuts: bool
    
    enable_cov_cuts: bool
    
    enable_clq_cuts: bool
    
    mip_gap: float
    
    time_limit: int
    
    output_freq: int
    
    output_delay: int
    
    enable_presolve: bool
    
    enable_binarization: bool
    
    def __init__(self) -> None:
        """
        Cython signature: void SolverParam()
        Hold the parameters of the LP solver
        """
        ... 


class __LPWrapper_Type:
    None
    UNBOUNDED : int
    LOWER_BOUND_ONLY : int
    UPPER_BOUND_ONLY : int
    DOUBLE_BOUNDED : int
    FIXED : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __SOLVER:
    None
    SOLVER_GLPK : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __Sense:
    None
    MIN : int
    MAX : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __SolverStatus:
    None
    UNDEFINED : int
    OPTIMAL : int
    FEASIBLE : int
    NO_FEASIBLE_SOL : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class ValueType:
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


class __VariableType:
    None
    CONTINUOUS : int
    INTEGER : int
    BINARY : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __WriteFormat:
    None
    FORMAT_LP : int
    FORMAT_MPS : int
    FORMAT_GLPK : int

    def getMapping(self) -> Dict[int, str]:
       ... 

