from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class CVTerm:
    """
    Cython implementation of _CVTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTerm.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVTerm()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVTerm ) -> None:
        """
        Cython signature: void CVTerm(CVTerm &)
        """
        ...
    
    def setAccession(self, accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setAccession(String accession)
        Sets the accession string of the term
        """
        ...
    
    def getAccession(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getAccession()
        Returns the accession string of the term
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the term
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the term
        """
        ...
    
    def setCVIdentifierRef(self, cv_id_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCVIdentifierRef(String cv_id_ref)
        Sets the CV identifier reference string, e.g. UO for unit obo
        """
        ...
    
    def getCVIdentifierRef(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCVIdentifierRef()
        Returns the CV identifier reference string
        """
        ...
    
    def getValue(self) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getValue()
        Returns the value of the term
        """
        ...
    
    def setValue(self, value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setValue(DataValue value)
        Sets the value of the term
        """
        ...
    
    def setUnit(self, unit: Unit ) -> None:
        """
        Cython signature: void setUnit(Unit & unit)
        Sets the unit of the term
        """
        ...
    
    def getUnit(self) -> Unit:
        """
        Cython signature: Unit getUnit()
        Returns the unit
        """
        ...
    
    def hasValue(self) -> bool:
        """
        Cython signature: bool hasValue()
        Checks whether the term has a value
        """
        ...
    
    def hasUnit(self) -> bool:
        """
        Cython signature: bool hasUnit()
        Checks whether the term has a unit
        """
        ...
    
    def __richcmp__(self, other: CVTerm, op: int) -> Any:
        ... 


class ConvexHull2D:
    """
    Cython implementation of _ConvexHull2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConvexHull2D.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ConvexHull2D()
        """
        ...
    
    @overload
    def __init__(self, in_0: ConvexHull2D ) -> None:
        """
        Cython signature: void ConvexHull2D(ConvexHull2D &)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Removes all points
        """
        ...
    
    def compress(self) -> int:
        """
        Cython signature: size_t compress()
        Allows to reduce the disk/memory footprint of a hull
        """
        ...
    
    def expandToBoundingBox(self) -> None:
        """
        Cython signature: void expandToBoundingBox()
        Expand a convex hull to its bounding box.
        """
        ...
    
    def addPoint(self, point: Union[Sequence[int], Sequence[float]] ) -> bool:
        """
        Cython signature: bool addPoint(DPosition2 point)
        Adds a point to the hull if it is not already contained. Returns if the point was added. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
        """
        ...
    
    def addPoints(self, points: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None:
        """
        Cython signature: void addPoints(libcpp_vector[DPosition2] points)
        Adds points to the hull if it is not already contained. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
        """
        ...
    
    def encloses(self, in_0: Union[Sequence[int], Sequence[float]] ) -> bool:
        """
        Cython signature: bool encloses(DPosition2)
        Returns if the `point` lies in the feature hull
        """
        ...
    
    def getHullPoints(self) -> '_np.ndarray[Any, _np.dtype[_np.float32]]':
        """
        Cython signature: libcpp_vector[DPosition2] getHullPoints()
        Accessor for the outer points
        """
        ...
    
    def setHullPoints(self, in_0: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None:
        """
        Cython signature: void setHullPoints(libcpp_vector[DPosition2])
        Accessor for the outer(!) points (no checking is performed if this is actually a convex hull)
        """
        ...
    
    def getBoundingBox(self) -> DBoundingBox2:
        """
        Cython signature: DBoundingBox2 getBoundingBox()
        Returns the bounding box of the feature hull points
        """
        ...
    
    def __richcmp__(self, other: ConvexHull2D, op: int) -> Any:
        ... 


class CrossLinkSpectrumMatch:
    """
    Cython implementation of _CrossLinkSpectrumMatch

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1CrossLinkSpectrumMatch.html>`_
    """
    
    cross_link: ProteinProteinCrossLink
    
    scan_index_light: int
    
    scan_index_heavy: int
    
    score: float
    
    rank: int
    
    xquest_score: float
    
    pre_score: float
    
    percTIC: float
    
    wTIC: float
    
    wTICold: float
    
    int_sum: float
    
    intsum_alpha: float
    
    intsum_beta: float
    
    total_current: float
    
    precursor_error_ppm: float
    
    match_odds: float
    
    match_odds_alpha: float
    
    match_odds_beta: float
    
    log_occupancy: float
    
    log_occupancy_alpha: float
    
    log_occupancy_beta: float
    
    xcorrx_max: float
    
    xcorrc_max: float
    
    matched_linear_alpha: int
    
    matched_linear_beta: int
    
    matched_xlink_alpha: int
    
    matched_xlink_beta: int
    
    num_iso_peaks_mean: float
    
    num_iso_peaks_mean_linear_alpha: float
    
    num_iso_peaks_mean_linear_beta: float
    
    num_iso_peaks_mean_xlinks_alpha: float
    
    num_iso_peaks_mean_xlinks_beta: float
    
    ppm_error_abs_sum_linear_alpha: float
    
    ppm_error_abs_sum_linear_beta: float
    
    ppm_error_abs_sum_xlinks_alpha: float
    
    ppm_error_abs_sum_xlinks_beta: float
    
    ppm_error_abs_sum_linear: float
    
    ppm_error_abs_sum_xlinks: float
    
    ppm_error_abs_sum_alpha: float
    
    ppm_error_abs_sum_beta: float
    
    ppm_error_abs_sum: float
    
    precursor_correction: int
    
    precursor_total_intensity: float
    
    precursor_target_intensity: float
    
    precursor_signal_proportion: float
    
    precursor_target_peak_count: int
    
    precursor_residual_peak_count: int
    
    frag_annotations: List[PeptideHit_PeakAnnotation]
    
    peptide_id_index: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CrossLinkSpectrumMatch()
        """
        ...
    
    @overload
    def __init__(self, in_0: CrossLinkSpectrumMatch ) -> None:
        """
        Cython signature: void CrossLinkSpectrumMatch(CrossLinkSpectrumMatch &)
        """
        ... 


class CsvFile:
    """
    Cython implementation of _CsvFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CsvFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CsvFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: CsvFile ) -> None:
        """
        Cython signature: void CsvFile(CsvFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , is_: bytes , ie_: bool , first_n: int ) -> None:
        """
        Cython signature: void load(const String & filename, char is_, bool ie_, int first_n)
        Loads data from a text file
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & filename)
        Stores the buffer's content into a file
        """
        ...
    
    def addRow(self, list: List[bytes] ) -> None:
        """
        Cython signature: void addRow(const StringList & list)
        Add a row to the buffer
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Clears the buffer
        """
        ...
    
    def getRow(self, row: int , list: List[bytes] ) -> bool:
        """
        Cython signature: bool getRow(int row, StringList & list)
        Writes all items from a row to list
        """
        ... 


class FeatureFindingMetabo:
    """
    Cython implementation of _FeatureFindingMetabo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFindingMetabo.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FeatureFindingMetabo()
        Method for the assembly of mass traces belonging to the same isotope
        pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
        and isotope abundances
        """
        ...
    
    @overload
    def __init__(self, in_0: FeatureFindingMetabo ) -> None:
        """
        Cython signature: void FeatureFindingMetabo(FeatureFindingMetabo &)
        """
        ...
    
    def run(self, input_mtraces: List[Kernel_MassTrace] , output_featmap: FeatureMap , output_chromatograms: List[List[MSChromatogram]] ) -> None:
        """
        Cython signature: void run(libcpp_vector[Kernel_MassTrace] input_mtraces, FeatureMap & output_featmap, libcpp_vector[libcpp_vector[MSChromatogram]] & output_chromatograms)
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


class IsotopeModel:
    """
    Cython implementation of _IsotopeModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeModel.html>`_

    Isotope distribution approximated using linear interpolation
    
    This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query)
    If you only want the distribution (no widening), use either
    EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
    or
    IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine
    
    Peak widening is achieved by either a Gaussian or Lorentzian shape
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsotopeModel()
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopeModel ) -> None:
        """
        Cython signature: void IsotopeModel(IsotopeModel &)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: unsigned int getCharge()
        """
        ...
    
    def setOffset(self, offset: float ) -> None:
        """
        Cython signature: void setOffset(double offset)
        Set the offset of the model
        
        The whole model will be shifted to the new offset without being computing all over
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model
        """
        ...
    
    def getOffset(self) -> float:
        """
        Cython signature: double getOffset()
        Get the offset of the model
        """
        ...
    
    def getFormula(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula()
        Return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
        """
        ...
    
    def setSamples(self, formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void setSamples(EmpiricalFormula & formula)
        Set sample/supporting points of interpolation
        """
        ...
    
    def getCenter(self) -> float:
        """
        Cython signature: double getCenter()
        Get the center of the Isotope model
        
        This is a m/z-value not necessarily the monoisotopic mass
        """
        ...
    
    def getIsotopeDistribution(self) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution getIsotopeDistribution()
        Get the Isotope distribution (without widening) from the last setSamples() call
        
        Useful to determine the number of isotopes that the model contains and their position
        """
        ...
    Averagines : __Averagines 


class MassAnalyzer:
    """
    Cython implementation of _MassAnalyzer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassAnalyzer.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassAnalyzer()
        """
        ...
    
    @overload
    def __init__(self, in_0: MassAnalyzer ) -> None:
        """
        Cython signature: void MassAnalyzer(MassAnalyzer &)
        """
        ...
    
    def getType(self) -> int:
        """
        Cython signature: AnalyzerType getType()
        Returns the analyzer type
        """
        ...
    
    def setType(self, type: int ) -> None:
        """
        Cython signature: void setType(AnalyzerType type)
        Sets the analyzer type
        """
        ...
    
    def getResolutionMethod(self) -> int:
        """
        Cython signature: ResolutionMethod getResolutionMethod()
        Returns the method used for determination of the resolution
        """
        ...
    
    def setResolutionMethod(self, resolution_method: int ) -> None:
        """
        Cython signature: void setResolutionMethod(ResolutionMethod resolution_method)
        Sets the method used for determination of the resolution
        """
        ...
    
    def getResolutionType(self) -> int:
        """
        Cython signature: ResolutionType getResolutionType()
        Returns the resolution type
        """
        ...
    
    def setResolutionType(self, resolution_type: int ) -> None:
        """
        Cython signature: void setResolutionType(ResolutionType resolution_type)
        Sets the resolution type
        """
        ...
    
    def getScanDirection(self) -> int:
        """
        Cython signature: ScanDirection getScanDirection()
        Returns the direction of scanning
        """
        ...
    
    def setScanDirection(self, scan_direction: int ) -> None:
        """
        Cython signature: void setScanDirection(ScanDirection scan_direction)
        Sets the direction of scanning
        """
        ...
    
    def getScanLaw(self) -> int:
        """
        Cython signature: ScanLaw getScanLaw()
        Returns the scan law
        """
        ...
    
    def setScanLaw(self, scan_law: int ) -> None:
        """
        Cython signature: void setScanLaw(ScanLaw scan_law)
        Sets the scan law
        """
        ...
    
    def getReflectronState(self) -> int:
        """
        Cython signature: ReflectronState getReflectronState()
        Returns the reflectron state (for TOF)
        """
        ...
    
    def setReflectronState(self, reflecton_state: int ) -> None:
        """
        Cython signature: void setReflectronState(ReflectronState reflecton_state)
        Sets the reflectron state (for TOF)
        """
        ...
    
    def getResolution(self) -> float:
        """
        Cython signature: double getResolution()
        Returns the resolution. The maximum m/z value at which two peaks can be resolved, according to one of the standard measures
        """
        ...
    
    def setResolution(self, resolution: float ) -> None:
        """
        Cython signature: void setResolution(double resolution)
        Sets the resolution
        """
        ...
    
    def getAccuracy(self) -> float:
        """
        Cython signature: double getAccuracy()
        Returns the mass accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)
        """
        ...
    
    def setAccuracy(self, accuracy: float ) -> None:
        """
        Cython signature: void setAccuracy(double accuracy)
        Sets the accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)
        """
        ...
    
    def getScanRate(self) -> float:
        """
        Cython signature: double getScanRate()
        Returns the scan rate (in s)
        """
        ...
    
    def setScanRate(self, scan_rate: float ) -> None:
        """
        Cython signature: void setScanRate(double scan_rate)
        Sets the scan rate (in s)
        """
        ...
    
    def getScanTime(self) -> float:
        """
        Cython signature: double getScanTime()
        Returns the scan time for a single scan (in s)
        """
        ...
    
    def setScanTime(self, scan_time: float ) -> None:
        """
        Cython signature: void setScanTime(double scan_time)
        Sets the scan time for a single scan (in s)
        """
        ...
    
    def getTOFTotalPathLength(self) -> float:
        """
        Cython signature: double getTOFTotalPathLength()
        Returns the path length for a TOF mass analyzer (in meter)
        """
        ...
    
    def setTOFTotalPathLength(self, TOF_total_path_length: float ) -> None:
        """
        Cython signature: void setTOFTotalPathLength(double TOF_total_path_length)
        Sets the path length for a TOF mass analyzer (in meter)
        """
        ...
    
    def getIsolationWidth(self) -> float:
        """
        Cython signature: double getIsolationWidth()
        Returns the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
        """
        ...
    
    def setIsolationWidth(self, isolation_width: float ) -> None:
        """
        Cython signature: void setIsolationWidth(double isolation_width)
        Sets the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
        """
        ...
    
    def getFinalMSExponent(self) -> int:
        """
        Cython signature: int getFinalMSExponent()
        Returns the final MS exponent
        """
        ...
    
    def setFinalMSExponent(self, final_MS_exponent: int ) -> None:
        """
        Cython signature: void setFinalMSExponent(int final_MS_exponent)
        Sets the final MS exponent
        """
        ...
    
    def getMagneticFieldStrength(self) -> float:
        """
        Cython signature: double getMagneticFieldStrength()
        Returns the strength of the magnetic field (in T)
        """
        ...
    
    def setMagneticFieldStrength(self, magnetic_field_strength: float ) -> None:
        """
        Cython signature: void setMagneticFieldStrength(double magnetic_field_strength)
        Sets the strength of the magnetic field (in T)
        """
        ...
    
    def getOrder(self) -> int:
        """
        Cython signature: int getOrder()
        Returns the position of this part in the whole Instrument
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
    
    def __richcmp__(self, other: MassAnalyzer, op: int) -> Any:
        ...
    AnalyzerType : __AnalyzerType
    ReflectronState : __ReflectronState
    ResolutionMethod : __ResolutionMethod
    ResolutionType : __ResolutionType
    ScanDirection : __ScanDirection
    ScanLaw : __ScanLaw 


class MzDataFile:
    """
    Cython implementation of _MzDataFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzDataFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzDataFile()
        File adapter for MzData files
        """
        ...
    
    @overload
    def __init__(self, in_0: MzDataFile ) -> None:
        """
        Cython signature: void MzDataFile(MzDataFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , map: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & map)
        Loads a map from a MzData file
        
        
        :param filename: Directory of the file with the file name
        :param map: It has to be a MSExperiment or have the same interface
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , map: MSExperiment ) -> None:
        """
        Cython signature: void store(const String & filename, MSExperiment & map)
        Stores a map in a MzData file
        
        
        :param filename: Directory of the file with the file name
        :param map: It has to be a MSExperiment or have the same interface
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be created
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
    
    def isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool:
        """
        Cython signature: bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
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


class NucleicAcidSpectrumGenerator:
    """
    Cython implementation of _NucleicAcidSpectrumGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NucleicAcidSpectrumGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void NucleicAcidSpectrumGenerator()
        """
        ...
    
    @overload
    def __init__(self, in_0: NucleicAcidSpectrumGenerator ) -> None:
        """
        Cython signature: void NucleicAcidSpectrumGenerator(NucleicAcidSpectrumGenerator &)
        """
        ...
    
    def getSpectrum(self, spec: MSSpectrum , oligo: NASequence , min_charge: int , max_charge: int ) -> None:
        """
        Cython signature: void getSpectrum(MSSpectrum & spec, NASequence & oligo, int min_charge, int max_charge)
        Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
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


class OnDiscMSExperiment:
    """
    Cython implementation of _OnDiscMSExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OnDiscMSExperiment.html>`_

    Representation of a mass spectrometry experiment on disk.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OnDiscMSExperiment()
        """
        ...
    
    @overload
    def __init__(self, in_0: OnDiscMSExperiment ) -> None:
        """
        Cython signature: void OnDiscMSExperiment(OnDiscMSExperiment &)
        """
        ...
    
    @overload
    def openFile(self, filename: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool openFile(String filename)
        """
        ...
    
    @overload
    def openFile(self, filename: Union[bytes, str, String] , skipLoadingMetaData: bool ) -> bool:
        """
        Cython signature: bool openFile(String filename, bool skipLoadingMetaData)
        Open a specific file on disk
        
        This tries to read the indexed mzML by parsing the index and then reading the meta information into memory
        
        returns: Whether the parsing of the file was successful (if false, the file most likely was not an indexed mzML file)
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        Returns the total number of spectra available
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        Returns the total number of chromatograms available
        """
        ...
    
    def getExperimentalSettings(self) -> ExperimentalSettings:
        """
        Cython signature: shared_ptr[const ExperimentalSettings] getExperimentalSettings()
        Returns the meta information of this experiment (const access)
        """
        ...
    
    def getMetaData(self) -> MSExperiment:
        """
        Cython signature: shared_ptr[MSExperiment] getMetaData()
        Returns the meta information of this experiment
        """
        ...
    
    def getSpectrum(self, id: int ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum getSpectrum(size_t id)
        Returns a single spectrum
        
        
        :param id: The index of the spectrum
        """
        ...
    
    def getSpectrumByNativeId(self, id: Union[bytes, str, String] ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum getSpectrumByNativeId(String id)
        Returns a single spectrum
        
        
        :param id: The native identifier of the spectrum
        """
        ...
    
    def getChromatogram(self, id: int ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getChromatogram(size_t id)
        Returns a single chromatogram
        
        
        :param id: The index of the chromatogram
        """
        ...
    
    def getChromatogramByNativeId(self, id: Union[bytes, str, String] ) -> MSChromatogram:
        """
        Cython signature: MSChromatogram getChromatogramByNativeId(String id)
        Returns a single chromatogram
        
        
        :param id: The native identifier of the chromatogram
        """
        ...
    
    def getSpectrumById(self, id_: int ) -> _Interfaces_Spectrum:
        """
        Cython signature: shared_ptr[_Interfaces_Spectrum] getSpectrumById(int id_)
        Returns a single spectrum
        """
        ...
    
    def getChromatogramById(self, id_: int ) -> _Interfaces_Chromatogram:
        """
        Cython signature: shared_ptr[_Interfaces_Chromatogram] getChromatogramById(int id_)
        Returns a single chromatogram
        """
        ...
    
    def setSkipXMLChecks(self, skip: bool ) -> None:
        """
        Cython signature: void setSkipXMLChecks(bool skip)
        Sets whether to skip some XML checks and be fast instead
        """
        ... 


class PeakFileOptions:
    """
    Cython implementation of _PeakFileOptions

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakFileOptions.html>`_

    Options for loading files containing peak data
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakFileOptions()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakFileOptions ) -> None:
        """
        Cython signature: void PeakFileOptions(PeakFileOptions &)
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
    
    def setWriteSupplementalData(self, in_0: bool ) -> None:
        """
        Cython signature: void setWriteSupplementalData(bool)
        Sets whether or not to write supplemental peak data in MzData files
        """
        ...
    
    def getWriteSupplementalData(self) -> bool:
        """
        Cython signature: bool getWriteSupplementalData()
        Returns whether or not to write supplemental peak data in MzData files
        """
        ...
    
    def setMSLevels(self, levels: List[int] ) -> None:
        """
        Cython signature: void setMSLevels(libcpp_vector[int] levels)
        Sets the desired MS levels for peaks to load
        """
        ...
    
    def addMSLevel(self, level: int ) -> None:
        """
        Cython signature: void addMSLevel(int level)
        Adds a desired MS level for peaks to load
        """
        ...
    
    def clearMSLevels(self) -> None:
        """
        Cython signature: void clearMSLevels()
        Clears the MS levels
        """
        ...
    
    def hasMSLevels(self) -> bool:
        """
        Cython signature: bool hasMSLevels()
        Returns true, if MS levels have been set
        """
        ...
    
    def containsMSLevel(self, level: int ) -> bool:
        """
        Cython signature: bool containsMSLevel(int level)
        Returns true, if MS level `level` has been set
        """
        ...
    
    def getMSLevels(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getMSLevels()
        Returns the set MS levels
        """
        ...
    
    def setCompression(self, in_0: bool ) -> None:
        """
        Cython signature: void setCompression(bool)
        Sets if data should be compressed when writing
        """
        ...
    
    def getCompression(self) -> bool:
        """
        Cython signature: bool getCompression()
        Returns true, if data should be compressed when writing
        """
        ...
    
    def setMz32Bit(self, mz_32_bit: bool ) -> None:
        """
        Cython signature: void setMz32Bit(bool mz_32_bit)
        Sets if mz-data and rt-data should be stored with 32bit or 64bit precision
        """
        ...
    
    def getMz32Bit(self) -> bool:
        """
        Cython signature: bool getMz32Bit()
        Returns true, if mz-data and rt-data should be stored with 32bit precision
        """
        ...
    
    def setIntensity32Bit(self, int_32_bit: bool ) -> None:
        """
        Cython signature: void setIntensity32Bit(bool int_32_bit)
        Sets if intensity data should be stored with 32bit or 64bit precision
        """
        ...
    
    def getIntensity32Bit(self) -> bool:
        """
        Cython signature: bool getIntensity32Bit()
        Returns true, if intensity data should be stored with 32bit precision
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
    
    def getMaxDataPoolSize(self) -> int:
        """
        Cython signature: size_t getMaxDataPoolSize()
        Returns maximal size of the data pool
        """
        ...
    
    def setMaxDataPoolSize(self, s: int ) -> None:
        """
        Cython signature: void setMaxDataPoolSize(size_t s)
        Sets maximal size of the data pool
        """
        ...
    
    def setSortSpectraByMZ(self, doSort: bool ) -> None:
        """
        Cython signature: void setSortSpectraByMZ(bool doSort)
        Sets whether or not to sort peaks in spectra
        """
        ...
    
    def getSortSpectraByMZ(self) -> bool:
        """
        Cython signature: bool getSortSpectraByMZ()
        Returns whether or not peaks in spectra should be sorted
        """
        ...
    
    def setSortChromatogramsByRT(self, doSort: bool ) -> None:
        """
        Cython signature: void setSortChromatogramsByRT(bool doSort)
        Sets whether or not to sort peaks in chromatograms
        """
        ...
    
    def getSortChromatogramsByRT(self) -> bool:
        """
        Cython signature: bool getSortChromatogramsByRT()
        Returns whether or not peaks in chromatograms should be sorted
        """
        ...
    
    def hasFilters(self) -> bool:
        """
        Cython signature: bool hasFilters()
        """
        ...
    
    def setFillData(self, only: bool ) -> None:
        """
        Cython signature: void setFillData(bool only)
        Sets whether to fill the actual data into the container (spectrum/chromatogram)
        """
        ...
    
    def getFillData(self) -> bool:
        """
        Cython signature: bool getFillData()
        Returns whether to fill the actual data into the container (spectrum/chromatogram)
        """
        ...
    
    def setSkipXMLChecks(self, only: bool ) -> None:
        """
        Cython signature: void setSkipXMLChecks(bool only)
        Sets whether to skip some XML checks and be fast instead
        """
        ...
    
    def getSkipXMLChecks(self) -> bool:
        """
        Cython signature: bool getSkipXMLChecks()
        Returns whether to skip some XML checks and be fast instead
        """
        ...
    
    def getWriteIndex(self) -> bool:
        """
        Cython signature: bool getWriteIndex()
        Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        """
        ...
    
    def setWriteIndex(self, write_index: bool ) -> None:
        """
        Cython signature: void setWriteIndex(bool write_index)
        Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        """
        ...
    
    def getNumpressConfigurationMassTime(self) -> NumpressConfig:
        """
        Cython signature: NumpressConfig getNumpressConfigurationMassTime()
        Sets numpress configuration options for m/z or rt dimension
        """
        ...
    
    def setNumpressConfigurationMassTime(self, config: NumpressConfig ) -> None:
        """
        Cython signature: void setNumpressConfigurationMassTime(NumpressConfig config)
        Returns numpress configuration options for m/z or rt dimension
        """
        ...
    
    def getNumpressConfigurationIntensity(self) -> NumpressConfig:
        """
        Cython signature: NumpressConfig getNumpressConfigurationIntensity()
        Sets numpress configuration options for intensity dimension
        """
        ...
    
    def setNumpressConfigurationIntensity(self, config: NumpressConfig ) -> None:
        """
        Cython signature: void setNumpressConfigurationIntensity(NumpressConfig config)
        Returns numpress configuration options for intensity dimension
        """
        ...
    
    def getNumpressConfigurationFloatDataArray(self) -> NumpressConfig:
        """
        Cython signature: NumpressConfig getNumpressConfigurationFloatDataArray()
        Sets numpress configuration options for float data arrays
        """
        ...
    
    def setNumpressConfigurationFloatDataArray(self, config: NumpressConfig ) -> None:
        """
        Cython signature: void setNumpressConfigurationFloatDataArray(NumpressConfig config)
        Returns numpress configuration options for float data arrays
        """
        ...
    
    def setForceMQCompatability(self, forceMQ: bool ) -> None:
        """
        Cython signature: void setForceMQCompatability(bool forceMQ)
        [mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        """
        ...
    
    def getForceMQCompatability(self) -> bool:
        """
        Cython signature: bool getForceMQCompatability()
        [mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        """
        ...
    
    def setForceTPPCompatability(self, forceTPP: bool ) -> None:
        """
        Cython signature: void setForceTPPCompatability(bool forceTPP)
        [ mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        """
        ...
    
    def getForceTPPCompatability(self) -> bool:
        """
        Cython signature: bool getForceTPPCompatability()
        [mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        """
        ... 


class SequestOutfile:
    """
    Cython implementation of _SequestOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SequestOutfile.html>`_

    Representation of a Sequest output file
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SequestOutfile()
        Representation of a Sequest output file
        """
        ...
    
    @overload
    def __init__(self, in_0: SequestOutfile ) -> None:
        """
        Cython signature: void SequestOutfile(SequestOutfile &)
        """
        ...
    
    def load(self, result_filename: Union[bytes, str, String] , peptide_identifications: PeptideIdentificationList , protein_identification: ProteinIdentification , p_value_threshold: float , pvalues: List[float] , database: Union[bytes, str, String] , ignore_proteins_per_peptide: bool ) -> None:
        """
        Cython signature: void load(const String & result_filename, PeptideIdentificationList & peptide_identifications, ProteinIdentification & protein_identification, double p_value_threshold, libcpp_vector[double] & pvalues, const String & database, bool ignore_proteins_per_peptide)
        Loads data from a Sequest outfile
        
        :param result_filename: The file to be loaded
        :param peptide_identifications: The identifications
        :param protein_identification: The protein identifications
        :param p_value_threshold: The significance level (for the peptide hit scores)
        :param pvalues: A list with the pvalues of the peptides (pvalues computed with peptide prophet)
        :param database: The database used for the search
        :param ignore_proteins_per_peptide: This is a hack to deal with files that use a suffix like "+1" in column "Reference", but do not actually list extra protein references in subsequent lines
        """
        ...
    
    def getColumns(self, line: Union[bytes, str, String] , substrings: List[bytes] , number_of_columns: int , reference_column: int ) -> bool:
        """
        Cython signature: bool getColumns(const String & line, libcpp_vector[String] & substrings, size_t number_of_columns, size_t reference_column)
        Retrieves columns from a Sequest outfile line
        """
        ...
    
    def getACAndACType(self, line: Union[bytes, str, String] , accession: String , accession_type: String ) -> None:
        """
        Cython signature: void getACAndACType(String line, String & accession, String & accession_type)
        Retrieves the accession type and accession number from a protein description line
        """
        ...
    
    def __richcmp__(self, other: SequestOutfile, op: int) -> Any:
        ... 


class StringView:
    """
    Cython implementation of _StringView

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1StringView.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void StringView()
        """
        ...
    
    @overload
    def __init__(self, in_0: bytes ) -> None:
        """
        Cython signature: void StringView(const libcpp_string &)
        """
        ...
    
    @overload
    def __init__(self, in_0: StringView ) -> None:
        """
        Cython signature: void StringView(StringView &)
        """
        ...
    
    def substr(self, start: int , end: int ) -> StringView:
        """
        Cython signature: StringView substr(size_t start, size_t end)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def getString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getString()
        """
        ...
    
    def __richcmp__(self, other: StringView, op: int) -> Any:
        ... 


class Unit:
    """
    Cython implementation of _Unit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Unit.html>`_
    """
    
    accession: Union[bytes, str, String]
    
    name: Union[bytes, str, String]
    
    cv_ref: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Unit()
        """
        ...
    
    @overload
    def __init__(self, in_0: Unit ) -> None:
        """
        Cython signature: void Unit(Unit)
        """
        ...
    
    @overload
    def __init__(self, p_accession: Union[bytes, str, String] , p_name: Union[bytes, str, String] , p_cv_ref: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void Unit(const String & p_accession, const String & p_name, const String & p_cv_ref)
        """
        ...
    
    def __richcmp__(self, other: Unit, op: int) -> Any:
        ... 


class XMLHandler:
    """
    Cython implementation of _XMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1XMLHandler.html>`_
    """
    
    def __init__(self, filename: Union[bytes, str, String] , version: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void XMLHandler(const String & filename, const String & version)
        """
        ...
    
    def reset(self) -> None:
        """
        Cython signature: void reset()
        """
        ...
    
    def error(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None:
        """
        Cython signature: void error(ActionMode mode, const String & msg, unsigned int line, unsigned int column)
        """
        ...
    
    def warning(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None:
        """
        Cython signature: void warning(ActionMode mode, const String & msg, unsigned int line, unsigned int column)
        """
        ...
    ActionMode : __ActionMode 


class __ActionMode:
    None
    LOAD : int
    STORE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __AnalyzerType:
    None
    ANALYZERNULL : int
    QUADRUPOLE : int
    PAULIONTRAP : int
    RADIALEJECTIONLINEARIONTRAP : int
    AXIALEJECTIONLINEARIONTRAP : int
    TOF : int
    SECTOR : int
    FOURIERTRANSFORM : int
    IONSTORAGE : int
    ESA : int
    IT : int
    SWIFT : int
    CYCLOTRON : int
    ORBITRAP : int
    LIT : int
    SIZE_OF_ANALYZERTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __Averagines:
    None
    C : int
    H : int
    N : int
    O : int
    S : int
    AVERAGINE_NUM : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ReflectronState:
    None
    REFLSTATENULL : int
    ON : int
    OFF : int
    NONE : int
    SIZE_OF_REFLECTRONSTATE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ResolutionMethod:
    None
    RESMETHNULL : int
    FWHM : int
    TENPERCENTVALLEY : int
    BASELINE : int
    SIZE_OF_RESOLUTIONMETHOD : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ResolutionType:
    None
    RESTYPENULL : int
    CONSTANT : int
    PROPORTIONAL : int
    SIZE_OF_RESOLUTIONTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ScanDirection:
    None
    SCANDIRNULL : int
    UP : int
    DOWN : int
    SIZE_OF_SCANDIRECTION : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ScanLaw:
    None
    SCANLAWNULL : int
    EXPONENTIAL : int
    LINEAR : int
    QUADRATIC : int
    SIZE_OF_SCANLAW : int

    def getMapping(self) -> Dict[int, str]:
       ... 

