from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class Attachment:
    """
    Cython implementation of _Attachment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::QcMLFile_1_1Attachment.html>`_
    """
    
    name: Union[bytes, str, String]
    
    id: Union[bytes, str, String]
    
    value: Union[bytes, str, String]
    
    cvRef: Union[bytes, str, String]
    
    cvAcc: Union[bytes, str, String]
    
    unitRef: Union[bytes, str, String]
    
    unitAcc: Union[bytes, str, String]
    
    binary: Union[bytes, str, String]
    
    qualityRef: Union[bytes, str, String]
    
    colTypes: List[bytes]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Attachment()
        """
        ...
    
    @overload
    def __init__(self, in_0: Attachment ) -> None:
        """
        Cython signature: void Attachment(Attachment &)
        """
        ...
    
    def toXMLString(self, indentation_level: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String toXMLString(unsigned int indentation_level)
        """
        ...
    
    def toCSVString(self, separator: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String toCSVString(String separator)
        """
        ...
    
    def __richcmp__(self, other: Attachment, op: int) -> Any:
        ... 


class BaseFeature:
    """
    Cython implementation of _BaseFeature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BaseFeature.html>`_
      -- Inherits from ['UniqueIdInterface', 'RichPeak2D']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BaseFeature()
        """
        ...
    
    @overload
    def __init__(self, in_0: BaseFeature ) -> None:
        """
        Cython signature: void BaseFeature(BaseFeature &)
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
    
    @overload
    def setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setPeptideIdentifications(PeptideIdentificationList & peptides)
        Sets the PeptideIdentification vector
        """
        ...
    
    @overload
    def setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setPeptideIdentifications(PeptideIdentificationList & peptides)
        Sets the PeptideIdentificationList
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
    
    def __richcmp__(self, other: BaseFeature, op: int) -> Any:
        ... 


class CachedSwathFileConsumer:
    """
    Cython implementation of _CachedSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CachedSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """
    
    @overload
    def __init__(self, in_0: CachedSwathFileConsumer ) -> None:
        """
        Cython signature: void CachedSwathFileConsumer(CachedSwathFileConsumer &)
        """
        ...
    
    @overload
    def __init__(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None:
        """
        Cython signature: void CachedSwathFileConsumer(String cachedir, String basename, size_t nr_ms1_spectra, libcpp_vector[int] nr_ms2_spectra)
        """
        ...
    
    def setExpectedSize(self, s: int , c: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t s, size_t c)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings exp)
        """
        ...
    
    def retrieveSwathMaps(self, maps: List[SwathMap] ) -> None:
        """
        Cython signature: void retrieveSwathMaps(libcpp_vector[SwathMap] & maps)
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        """
        ...
    
    def consumeChromatogram(self, c: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram & c)
        """
        ... 


class EmgModel:
    """
    Cython implementation of _EmgModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgModel.html>`_
      -- Inherits from ['InterpolationModel']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EmgModel()
        Exponentially modified gaussian distribution model for elution profiles
        """
        ...
    
    @overload
    def __init__(self, in_0: EmgModel ) -> None:
        """
        Cython signature: void EmgModel(EmgModel &)
        """
        ...
    
    def getIntensity(self, coord: float ) -> float:
        """
        Cython signature: double getIntensity(double coord)
        Access model predicted intensity at position 'pos'
        """
        ...
    
    def getScalingFactor(self) -> float:
        """
        Cython signature: double getScalingFactor()
        Returns the interpolation class
        """
        ...
    
    def setOffset(self, offset: float ) -> None:
        """
        Cython signature: void setOffset(double offset)
        Sets the offset of the model
        """
        ...
    
    def getCenter(self) -> float:
        """
        Cython signature: double getCenter()
        Returns the "center" of the model, particular definition (depends on the derived model)
        """
        ...
    
    def setSamples(self) -> None:
        """
        Cython signature: void setSamples()
        Sets sample/supporting points of interpolation wrt params
        """
        ...
    
    def setInterpolationStep(self, interpolation_step: float ) -> None:
        """
        Cython signature: void setInterpolationStep(double interpolation_step)
        Sets the interpolation step for the linear interpolation of the model
        """
        ...
    
    def setScalingFactor(self, scaling: float ) -> None:
        """
        Cython signature: void setScalingFactor(double scaling)
        Sets the scaling factor of the model
        """
        ...
    
    def getInterpolation(self) -> LinearInterpolation:
        """
        Cython signature: LinearInterpolation getInterpolation()
        Returns the interpolation class
        """
        ... 


class GaussTraceFitter:
    """
    Cython implementation of _GaussTraceFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GaussTraceFitter.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void GaussTraceFitter()
        Fitter for RT profiles using a Gaussian background model
        """
        ...
    
    @overload
    def __init__(self, in_0: GaussTraceFitter ) -> None:
        """
        Cython signature: void GaussTraceFitter(GaussTraceFitter &)
        """
        ...
    
    def fit(self, traces: MassTraces ) -> None:
        """
        Cython signature: void fit(MassTraces & traces)
        Override important methods
        """
        ...
    
    def getLowerRTBound(self) -> float:
        """
        Cython signature: double getLowerRTBound()
        Returns the lower RT bound
        """
        ...
    
    def getUpperRTBound(self) -> float:
        """
        Cython signature: double getUpperRTBound()
        Returns the upper RT bound
        """
        ...
    
    def getHeight(self) -> float:
        """
        Cython signature: double getHeight()
        Returns height of the fitted gaussian model
        """
        ...
    
    def getCenter(self) -> float:
        """
        Cython signature: double getCenter()
        Returns center of the fitted gaussian model
        """
        ...
    
    def getFWHM(self) -> float:
        """
        Cython signature: double getFWHM()
        Returns FWHM of the fitted gaussian model
        """
        ...
    
    def getSigma(self) -> float:
        """
        Cython signature: double getSigma()
        Returns Sigma of the fitted gaussian model
        """
        ...
    
    def checkMaximalRTSpan(self, max_rt_span: float ) -> bool:
        """
        Cython signature: bool checkMaximalRTSpan(double max_rt_span)
        """
        ...
    
    def checkMinimalRTSpan(self, rt_bounds: List[float, float] , min_rt_span: float ) -> bool:
        """
        Cython signature: bool checkMinimalRTSpan(libcpp_pair[double,double] & rt_bounds, double min_rt_span)
        """
        ...
    
    def computeTheoretical(self, trace: MassTrace , k: int ) -> float:
        """
        Cython signature: double computeTheoretical(MassTrace & trace, size_t k)
        """
        ...
    
    def getArea(self) -> float:
        """
        Cython signature: double getArea()
        Returns area of the fitted gaussian model
        """
        ...
    
    def getGnuplotFormula(self, trace: MassTrace , function_name: bytes , baseline: float , rt_shift: float ) -> Union[bytes, str, String]:
        """
        Cython signature: String getGnuplotFormula(MassTrace & trace, char function_name, double baseline, double rt_shift)
        """
        ...
    
    def getValue(self, rt: float ) -> float:
        """
        Cython signature: double getValue(double rt)
        Returns value of the fitted gaussian model
        """
        ... 


class InterpolationModel:
    """
    Cython implementation of _InterpolationModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InterpolationModel.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void InterpolationModel()
        Abstract class for 1D-models that are approximated using linear interpolation
        """
        ...
    
    @overload
    def __init__(self, in_0: InterpolationModel ) -> None:
        """
        Cython signature: void InterpolationModel(InterpolationModel &)
        """
        ...
    
    def getIntensity(self, coord: float ) -> float:
        """
        Cython signature: double getIntensity(double coord)
        Access model predicted intensity at position 'pos'
        """
        ...
    
    def getScalingFactor(self) -> float:
        """
        Cython signature: double getScalingFactor()
        Returns the interpolation class
        """
        ...
    
    def setOffset(self, offset: float ) -> None:
        """
        Cython signature: void setOffset(double offset)
        Sets the offset of the model
        """
        ...
    
    def getCenter(self) -> float:
        """
        Cython signature: double getCenter()
        Returns the "center" of the model, particular definition (depends on the derived model)
        """
        ...
    
    def setSamples(self) -> None:
        """
        Cython signature: void setSamples()
        Sets sample/supporting points of interpolation wrt params
        """
        ...
    
    def setInterpolationStep(self, interpolation_step: float ) -> None:
        """
        Cython signature: void setInterpolationStep(double interpolation_step)
        Sets the interpolation step for the linear interpolation of the model
        """
        ...
    
    def setScalingFactor(self, scaling: float ) -> None:
        """
        Cython signature: void setScalingFactor(double scaling)
        Sets the scaling factor of the model
        """
        ...
    
    def getInterpolation(self) -> LinearInterpolation:
        """
        Cython signature: LinearInterpolation getInterpolation()
        Returns the interpolation class
        """
        ... 


class MassDecompositionAlgorithm:
    """
    Cython implementation of _MassDecompositionAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassDecompositionAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MassDecompositionAlgorithm()
        """
        ...
    
    def getDecompositions(self, decomps: List[MassDecomposition] , weight: float ) -> None:
        """
        Cython signature: void getDecompositions(libcpp_vector[MassDecomposition] & decomps, double weight)
        Returns the possible decompositions given the weight
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


class MassTraceDetection:
    """
    Cython implementation of _MassTraceDetection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassTraceDetection.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MassTraceDetection()
        """
        ...
    
    @overload
    def __init__(self, in_0: MassTraceDetection ) -> None:
        """
        Cython signature: void MassTraceDetection(MassTraceDetection &)
        """
        ...
    
    def run(self, input_map: MSExperiment , traces: List[Kernel_MassTrace] , max_traces: int ) -> None:
        """
        Cython signature: void run(MSExperiment & input_map, libcpp_vector[Kernel_MassTrace] & traces, size_t max_traces)
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


class MetaInfoInterface:
    """
    Cython implementation of _MetaInfoInterface

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoInterface.html>`_

    Interface for classes that can store arbitrary meta information
    (Type-Name-Value tuples).
    
    MetaInfoInterface is a base class for all classes that use one MetaInfo
    object as member.  If you want to add meta information to a class, let it
    publicly inherit the MetaInfoInterface.  Meta information is an array of
    Type-Name-Value tuples.
    
    Usage:
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MetaInfoInterface()
        """
        ...
    
    @overload
    def __init__(self, in_0: MetaInfoInterface ) -> None:
        """
        Cython signature: void MetaInfoInterface(MetaInfoInterface &)
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
    
    def __richcmp__(self, other: MetaInfoInterface, op: int) -> Any:
        ... 


class MorphologicalFilter:
    """
    Cython implementation of _MorphologicalFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorphologicalFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MorphologicalFilter()
        """
        ...
    
    def filter(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void filter(MSSpectrum & spectrum)
        Applies the morphological filtering operation to an MSSpectrum
        
        If the size of the structuring element is given in 'Thomson', the number of data points for
        the structuring element is computed as follows:
        """
        ...
    
    def filterExperiment(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterExperiment(MSExperiment & exp)
        Applies the morphological filtering operation to an MSExperiment
        
        The size of the structuring element is computed for each spectrum individually, if it is given in 'Thomson'
        See the filtering method for MSSpectrum for details
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


class MzMLSwathFileConsumer:
    """
    Cython implementation of _MzMLSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """
    
    @overload
    def __init__(self, in_0: MzMLSwathFileConsumer ) -> None:
        """
        Cython signature: void MzMLSwathFileConsumer(MzMLSwathFileConsumer)
        """
        ...
    
    @overload
    def __init__(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None:
        """
        Cython signature: void MzMLSwathFileConsumer(String cachedir, String basename, size_t nr_ms1_spectra, libcpp_vector[int] nr_ms2_spectra)
        """
        ...
    
    @overload
    def __init__(self, known_window_boundaries: List[SwathMap] , cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None:
        """
        Cython signature: void MzMLSwathFileConsumer(libcpp_vector[SwathMap] known_window_boundaries, String cachedir, String basename, size_t nr_ms1_spectra, libcpp_vector[int] nr_ms2_spectra)
        """
        ...
    
    def setExpectedSize(self, s: int , c: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t s, size_t c)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings exp)
        """
        ...
    
    def retrieveSwathMaps(self, maps: List[SwathMap] ) -> None:
        """
        Cython signature: void retrieveSwathMaps(libcpp_vector[SwathMap] & maps)
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        """
        ...
    
    def consumeChromatogram(self, c: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram & c)
        """
        ... 


class RegularSwathFileConsumer:
    """
    Cython implementation of _RegularSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RegularSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RegularSwathFileConsumer()
        """
        ...
    
    @overload
    def __init__(self, in_0: RegularSwathFileConsumer ) -> None:
        """
        Cython signature: void RegularSwathFileConsumer(RegularSwathFileConsumer &)
        """
        ...
    
    def setExpectedSize(self, s: int , c: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t s, size_t c)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings exp)
        """
        ...
    
    def retrieveSwathMaps(self, maps: List[SwathMap] ) -> None:
        """
        Cython signature: void retrieveSwathMaps(libcpp_vector[SwathMap] & maps)
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        """
        ...
    
    def consumeChromatogram(self, c: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram & c)
        """
        ... 


class SimplePeak:
    """
    Cython implementation of _SimplePeak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimplePeak.html>`_
    """
    
    mz: float
    
    charge: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SimplePeak()
        A simple struct to represent peaks with mz and charge and sort them easily
        """
        ...
    
    @overload
    def __init__(self, mz: float , charge: int ) -> None:
        """
        Cython signature: void SimplePeak(double mz, int charge)
        """
        ...
    
    @overload
    def __init__(self, in_0: SimplePeak ) -> None:
        """
        Cython signature: void SimplePeak(SimplePeak &)
        """
        ... 


class SimpleTSGXLMS:
    """
    Cython implementation of _SimpleTSGXLMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimpleTSGXLMS.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SimpleTSGXLMS()
        Generates theoretical spectra for cross-linked peptides
        
        The spectra this class generates are vectors of SimplePeaks
        This class generates the same peak types as TheoreticalSpectrumGeneratorXLMS
        and the interface is very similar, but it is simpler and faster
        SimplePeak only contains an mz value and a charge. No intensity values
        or String annotations or other additional DataArrays are generated
        """
        ...
    
    @overload
    def __init__(self, in_0: SimpleTSGXLMS ) -> None:
        """
        Cython signature: void SimpleTSGXLMS(SimpleTSGXLMS &)
        """
        ...
    
    def getLinearIonSpectrum(self, spectrum: List[SimplePeak] , peptide: AASequence , link_pos: int , charge: int , link_pos_2: int ) -> None:
        """
        Cython signature: void getLinearIonSpectrum(libcpp_vector[SimplePeak] & spectrum, AASequence peptide, size_t link_pos, int charge, size_t link_pos_2)
        Generates fragment ions not containing the cross-linker for one peptide
        
        B-ions are generated from the beginning of the peptide up to the first linked position,
        y-ions are generated from the second linked position up the end of the peptide
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
        The generated ion types and other additional settings are determined by the tool parameters
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param peptide: The peptide to fragment
        :param link_pos: The position of the cross-linker on the given peptide
        :param charge: The maximal charge of the ions
        :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        ...
    
    @overload
    def getXLinkIonSpectrum(self, spectrum: List[SimplePeak] , peptide: AASequence , link_pos: int , precursor_mass: float , mincharge: int , maxcharge: int , link_pos_2: int ) -> None:
        """
        Cython signature: void getXLinkIonSpectrum(libcpp_vector[SimplePeak] & spectrum, AASequence peptide, size_t link_pos, double precursor_mass, int mincharge, int maxcharge, size_t link_pos_2)
        Generates fragment ions containing the cross-linker for one peptide
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
        Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
        a precursor mass for the two peptides and the linker is needed
        In the case of a loop link the precursor mass is the mass of the only peptide and the linker
        Although this function is more general, currently it is mainly used for loop-links and mono-links,
        because residues in the second, unknown peptide cannot be considered for possible neutral losses
        The generated ion types and other additional settings are determined by the tool parameters
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param peptide: The peptide to fragment
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        ...
    
    @overload
    def getXLinkIonSpectrum(self, spectrum: List[SimplePeak] , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None:
        """
        Cython signature: void getXLinkIonSpectrum(libcpp_vector[SimplePeak] & spectrum, ProteinProteinCrossLink crosslink, bool frag_alpha, int mincharge, int maxcharge)
        Generates fragment ions containing the cross-linker for a pair of peptides
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        This function generates neutral loss ions by considering both linked peptides
        Only one of the peptides, decided by @frag_alpha, is fragmented
        This simplifies the function, but it has to be called twice to get all fragments of a peptide pair
        The generated ion types and other additional settings are determined by the tool parameters
        This function is not suitable to generate fragments for mono-links or loop-links
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param crosslink: ProteinProteinCrossLink to be fragmented
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
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


class XLPrecursor:
    """
    Cython implementation of _XLPrecursor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XLPrecursor.html>`_
    """
    
    precursor_mass: float
    
    alpha_index: int
    
    beta_index: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void XLPrecursor()
        """
        ...
    
    @overload
    def __init__(self, in_0: XLPrecursor ) -> None:
        """
        Cython signature: void XLPrecursor(XLPrecursor &)
        """
        ... 


class AnnotationState:
    None
    FEATURE_ID_NONE : int
    FEATURE_ID_SINGLE : int
    FEATURE_ID_MULTIPLE_SAME : int
    FEATURE_ID_MULTIPLE_DIVERGENT : int
    SIZE_OF_ANNOTATIONSTATE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

