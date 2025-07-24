from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_TransformationModelBSpline_getDefaultParameters(params: Param ) -> None:
    """
    Cython signature: void getDefaultParameters(Param & params)
    """
    ...

def __static_ChromatogramExtractor_prepare_coordinates(output_chromatograms: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , targeted: TargetedExperiment , rt_extraction_window: float , ms1: bool , ms1_isotopes: int ) -> None:
    """
    Cython signature: void prepare_coordinates(libcpp_vector[shared_ptr[OSChromatogram]] & output_chromatograms, libcpp_vector[ExtractionCoordinates] & extraction_coordinates, TargetedExperiment & targeted, double rt_extraction_window, bool ms1, int ms1_isotopes)
    """
    ...


class ChromatogramExtractor:
    """
    Cython implementation of _ChromatogramExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramExtractor.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramExtractor()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramExtractor ) -> None:
        """
        Cython signature: void ChromatogramExtractor(ChromatogramExtractor &)
        """
        ...
    
    def extractChromatograms(self, input: SpectrumAccessOpenMS , output: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , mz_extraction_window: float , ppm: bool , im_extraction_window: float , filter: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void extractChromatograms(shared_ptr[SpectrumAccessOpenMS] input, libcpp_vector[shared_ptr[OSChromatogram]] & output, libcpp_vector[ExtractionCoordinates] extraction_coordinates, double mz_extraction_window, bool ppm, double im_extraction_window, String filter)
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
    
    prepare_coordinates: __static_ChromatogramExtractor_prepare_coordinates 


class DIAScoring:
    """
    Cython implementation of _DIAScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DIAScoring.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void DIAScoring()
        """
        ...
    
    def dia_ms1_massdiff_score(self, precursor_mz: float , spectrum: List[OSSpectrum] , im_range: RangeMobility , ppm_score: float ) -> bool:
        """
        Cython signature: bool dia_ms1_massdiff_score(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum]] spectrum, RangeMobility & im_range, double & ppm_score)
        """
        ...
    
    def dia_ms1_isotope_scores_averagine(self, precursor_mz: float , spectrum: List[OSSpectrum] , charge_state: int , im_range: RangeMobility , isotope_corr: float , isotope_overlap: float ) -> None:
        """
        Cython signature: void dia_ms1_isotope_scores_averagine(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum]] spectrum, int charge_state, RangeMobility & im_range, double & isotope_corr, double & isotope_overlap)
        """
        ...
    
    def dia_ms1_isotope_scores(self, precursor_mz: float , spectrum: List[OSSpectrum] , im_range: RangeMobility , isotope_corr: float , isotope_overlap: float , sum_formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void dia_ms1_isotope_scores(double precursor_mz, libcpp_vector[shared_ptr[OSSpectrum]] spectrum, RangeMobility & im_range, double & isotope_corr, double & isotope_overlap, EmpiricalFormula & sum_formula)
        """
        ...
    
    def score_with_isotopes(self, spectrum: List[OSSpectrum] , transitions: List[LightTransition] , im_range: RangeMobility , dotprod: float , manhattan: float ) -> None:
        """
        Cython signature: void score_with_isotopes(libcpp_vector[shared_ptr[OSSpectrum]] spectrum, libcpp_vector[LightTransition] transitions, RangeMobility & im_range, double & dotprod, double & manhattan)
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


class MSDataStoringConsumer:
    """
    Cython implementation of _MSDataStoringConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataStoringConsumer.html>`_

    Consumer class that simply stores the data
    
    This class is able to keep spectra and chromatograms passed to it in memory
    and the data can be accessed through getData()
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSDataStoringConsumer()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSDataStoringConsumer ) -> None:
        """
        Cython signature: void MSDataStoringConsumer(MSDataStoringConsumer &)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        Sets experimental settings
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        Sets expected size
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        """
        ...
    
    def consumeChromatogram(self, in_0: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram &)
        """
        ...
    
    def getData(self) -> MSExperiment:
        """
        Cython signature: MSExperiment getData()
        """
        ... 


class MasstraceCorrelator:
    """
    Cython implementation of _MasstraceCorrelator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MasstraceCorrelator.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MasstraceCorrelator()
        """
        ...
    
    @overload
    def __init__(self, in_0: MasstraceCorrelator ) -> None:
        """
        Cython signature: void MasstraceCorrelator(MasstraceCorrelator &)
        """
        ...
    
    def createPseudoSpectra(self, map_: ConsensusMap , pseudo_spectra: MSExperiment , min_peak_nr: int , min_correlation: float , max_lag: int , max_rt_apex_difference: float ) -> None:
        """
        Cython signature: void createPseudoSpectra(const ConsensusMap & map_, MSExperiment & pseudo_spectra, size_t min_peak_nr, double min_correlation, int max_lag, double max_rt_apex_difference)
        Compute pseudo-spectra from a set of (MS2) masstraces
        
        This function will take a set of masstraces (consensus map) as input and
        produce a vector of pseudo spectra as output (pseudo_spectra result
        vector).
        
        It basically makes an all-vs-all comparison of all masstraces against
        each other and scores them on how similar they are in their mass traces.
        
        This assumes that the consensus feature is only from one (SWATH) map
        This assumes that the consensus map is sorted by intensity
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


class MatrixDouble:
    """
    Cython implementation of _Matrix[double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Matrix[double].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MatrixDouble()
        """
        ...
    
    @overload
    def __init__(self, in_0: MatrixDouble ) -> None:
        """
        Cython signature: void MatrixDouble(MatrixDouble)
        """
        ...
    
    @overload
    def __init__(self, rows: int , cols: int , value: float ) -> None:
        """
        Cython signature: void MatrixDouble(size_t rows, size_t cols, double value)
        """
        ...
    
    def getValue(self, i: int , j: int ) -> float:
        """
        Cython signature: double getValue(size_t i, size_t j)
        """
        ...
    
    def setValue(self, i: int , j: int , value: float ) -> None:
        """
        Cython signature: void setValue(size_t i, size_t j, double value)
        """
        ...
    
    def rows(self) -> int:
        """
        Cython signature: size_t rows()
        """
        ...
    
    def cols(self) -> int:
        """
        Cython signature: size_t cols()
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def resize(self, rows: int , cols: int ) -> None:
        """
        Cython signature: void resize(size_t rows, size_t cols)
        """
        ... 


class NoopMSDataWritingConsumer:
    """
    Cython implementation of _NoopMSDataWritingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NoopMSDataWritingConsumer.html>`_

    Consumer class that perform no operation
    
    This is sometimes necessary to fulfill the requirement of passing an
    valid MSDataWritingConsumer object or pointer but no operation is
    required
    """
    
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void NoopMSDataWritingConsumer(String filename)
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
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        """
        ...
    
    def addDataProcessing(self, d: DataProcessing ) -> None:
        """
        Cython signature: void addDataProcessing(DataProcessing d)
        """
        ...
    
    def getNrSpectraWritten(self) -> int:
        """
        Cython signature: size_t getNrSpectraWritten()
        """
        ...
    
    def getNrChromatogramsWritten(self) -> int:
        """
        Cython signature: size_t getNrChromatogramsWritten()
        """
        ... 


class PlainMSDataWritingConsumer:
    """
    Cython implementation of _PlainMSDataWritingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PlainMSDataWritingConsumer.html>`_
    """
    
    def __init__(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void PlainMSDataWritingConsumer(String filename)
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
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        Set experimental settings for the whole file
        
        
        :param exp: Experimental settings to be used for this file (from this and the first spectrum/chromatogram, the class will deduce most of the header of the mzML file)
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        Set expected size of spectra and chromatograms to be written
        
        These numbers will be written in the spectrumList and chromatogramList
        tag in the mzML file. Therefore, these will contain wrong numbers if
        the expected size is not set correctly
        
        
        :param expectedSpectra: Number of spectra expected
        :param expectedChromatograms: Number of chromatograms expected
        """
        ...
    
    def addDataProcessing(self, d: DataProcessing ) -> None:
        """
        Cython signature: void addDataProcessing(DataProcessing d)
        Optionally add a data processing method to each chromatogram and spectrum
        
        The provided DataProcessing object will be added to each chromatogram
        and spectrum written to to the mzML file
        
        
        :param d: The DataProcessing object to be added
        """
        ...
    
    def getNrSpectraWritten(self) -> int:
        """
        Cython signature: size_t getNrSpectraWritten()
        Returns the number of spectra written
        """
        ...
    
    def getNrChromatogramsWritten(self) -> int:
        """
        Cython signature: size_t getNrChromatogramsWritten()
        Returns the number of chromatograms written
        """
        ...
    
    def setOptions(self, opt: PeakFileOptions ) -> None:
        """
        Cython signature: void setOptions(PeakFileOptions opt)
        """
        ...
    
    def getOptions(self) -> PeakFileOptions:
        """
        Cython signature: PeakFileOptions getOptions()
        """
        ... 


class QTCluster:
    """
    Cython implementation of _QTCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QTCluster.html>`_
    """
    
    def __init__(self, in_0: QTCluster ) -> None:
        """
        Cython signature: void QTCluster(QTCluster &)
        """
        ...
    
    def getCenterRT(self) -> float:
        """
        Cython signature: double getCenterRT()
        Returns the RT value of the cluster
        """
        ...
    
    def getCenterMZ(self) -> float:
        """
        Cython signature: double getCenterMZ()
        Returns the m/z value of the cluster center
        """
        ...
    
    def getXCoord(self) -> int:
        """
        Cython signature: int getXCoord()
        Returns the x coordinate in the grid
        """
        ...
    
    def getYCoord(self) -> int:
        """
        Cython signature: int getYCoord()
        Returns the y coordinate in the grid
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Returns the size of the cluster (number of elements, incl. center)
        """
        ...
    
    def getQuality(self) -> float:
        """
        Cython signature: double getQuality()
        Returns the cluster quality and recomputes if necessary
        """
        ...
    
    def getAnnotations(self) -> Set[AASequence]:
        """
        Cython signature: libcpp_set[AASequence] getAnnotations()
        Returns the set of peptide sequences annotated to the cluster center
        """
        ...
    
    def setInvalid(self) -> None:
        """
        Cython signature: void setInvalid()
        Sets current cluster as invalid (also frees some memory)
        """
        ...
    
    def isInvalid(self) -> bool:
        """
        Cython signature: bool isInvalid()
        Whether current cluster is invalid
        """
        ...
    
    def initializeCluster(self) -> None:
        """
        Cython signature: void initializeCluster()
        Has to be called before adding elements (calling
        """
        ...
    
    def finalizeCluster(self) -> None:
        """
        Cython signature: void finalizeCluster()
        Has to be called after adding elements (after calling
        """
        ...
    
    def __richcmp__(self, other: QTCluster, op: int) -> Any:
        ... 


class SwathFile:
    """
    Cython implementation of _SwathFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SwathFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: SwathFile ) -> None:
        """
        Cython signature: void SwathFile(SwathFile &)
        """
        ...
    
    def loadSplit(self, file_list: List[bytes] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]:
        """
        Cython signature: libcpp_vector[SwathMap] loadSplit(StringList file_list, String tmp, shared_ptr[ExperimentalSettings] exp_meta, String readoptions)
        Loads a Swath run from a list of split mzML files
        """
        ...
    
    def loadMzML(self, file_: Union[bytes, str, String] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]:
        """
        Cython signature: libcpp_vector[SwathMap] loadMzML(String file_, String tmp, shared_ptr[ExperimentalSettings] exp_meta, String readoptions)
        Loads a Swath run from a single mzML file
        
        Using the `plugin_consumer`, you can provide a custom consumer which will be chained
        into the process of loading the data and making it available (depending on `readoptions`).
        This is useful if you want to modify the data a priori or extract some other information using
        MSDataTransformingConsumer (for example). Make sure it leaves the data intact, such that the
        returned SwathMaps are actually useful
        
        :param file: Input filename
        :param tmp: Temporary directory (for cached data)
        :param exp_meta: Experimental metadata from mzML file
        :param readoptions: How are spectra accessed after reading - tradeoff between memory usage and time (disk caching)
        :param plugin_consumer: An intermediate custom consumer
        :returns: Swath maps for MS2 and MS1 (unless readoptions == split, which returns no data)
        """
        ...
    
    def loadMzXML(self, file_: Union[bytes, str, String] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]:
        """
        Cython signature: libcpp_vector[SwathMap] loadMzXML(String file_, String tmp, shared_ptr[ExperimentalSettings] exp_meta, String readoptions)
        Loads a Swath run from a single mzXML file
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


class TMTElevenPlexQuantitationMethod:
    """
    Cython implementation of _TMTElevenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTElevenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TMTElevenPlexQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: TMTElevenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void TMTElevenPlexQuantitationMethod(TMTElevenPlexQuantitationMethod &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def getChannelInformation(self) -> List[IsobaricChannelInformation]:
        """
        Cython signature: libcpp_vector[IsobaricChannelInformation] getChannelInformation()
        """
        ...
    
    def getNumberOfChannels(self) -> int:
        """
        Cython signature: size_t getNumberOfChannels()
        """
        ...
    
    def getIsotopeCorrectionMatrix(self) -> MatrixDouble:
        """
        Cython signature: MatrixDouble getIsotopeCorrectionMatrix()
        """
        ...
    
    def getReferenceChannel(self) -> int:
        """
        Cython signature: size_t getReferenceChannel()
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
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class TransformationModelBSpline:
    """
    Cython implementation of _TransformationModelBSpline

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelBSpline.html>`_
      -- Inherits from ['TransformationModel']
    """
    
    def __init__(self, data: List[TM_DataPoint] , params: Param ) -> None:
        """
        Cython signature: void TransformationModelBSpline(libcpp_vector[TM_DataPoint] & data, Param & params)
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
        Evaluates the model at the given values
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
    
    getDefaultParameters: __static_TransformationModelBSpline_getDefaultParameters 


class TransitionTSVFile:
    """
    Cython implementation of _TransitionTSVFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransitionTSVFile.html>`_
      -- Inherits from ['ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TransitionTSVFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: TransitionTSVFile ) -> None:
        """
        Cython signature: void TransitionTSVFile(TransitionTSVFile &)
        """
        ...
    
    def convertTargetedExperimentToTSV(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void convertTargetedExperimentToTSV(char * filename, TargetedExperiment & targeted_exp)
        Write out a targeted experiment (TraML structure) into a tsv file
        """
        ...
    
    @overload
    def convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void convertTSVToTargetedExperiment(char * filename, FileType filetype, TargetedExperiment & targeted_exp)
        Read in a tsv/mrm file and construct a targeted experiment (TraML structure)
        """
        ...
    
    @overload
    def convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None:
        """
        Cython signature: void convertTSVToTargetedExperiment(char * filename, FileType filetype, LightTargetedExperiment & targeted_exp)
        Read in a tsv file and construct a targeted experiment (Light transition structure)
        """
        ...
    
    def validateTargetedExperiment(self, targeted_exp: TargetedExperiment ) -> None:
        """
        Cython signature: void validateTargetedExperiment(TargetedExperiment targeted_exp)
        Validate a TargetedExperiment (check that all ids are unique)
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


class UniqueIdGenerator:
    """
    Cython implementation of _UniqueIdGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1UniqueIdGenerator.html>`_
    """
    
    def getUniqueId(self) -> int:
        """
        Cython signature: uint64_t getUniqueId()
        """
        ...
    
    def setSeed(self, in_0: int ) -> None:
        """
        Cython signature: void setSeed(uint64_t)
        """
        ...
    
    def getSeed(self) -> int:
        """
        Cython signature: uint64_t getSeed()
        """
        ... 


class XTandemInfile:
    """
    Cython implementation of _XTandemInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XTandemInfile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void XTandemInfile()
        """
        ...
    
    def setFragmentMassTolerance(self, tolerance: float ) -> None:
        """
        Cython signature: void setFragmentMassTolerance(double tolerance)
        """
        ...
    
    def getFragmentMassTolerance(self) -> float:
        """
        Cython signature: double getFragmentMassTolerance()
        """
        ...
    
    def setPrecursorMassTolerancePlus(self, tol: float ) -> None:
        """
        Cython signature: void setPrecursorMassTolerancePlus(double tol)
        """
        ...
    
    def getPrecursorMassTolerancePlus(self) -> float:
        """
        Cython signature: double getPrecursorMassTolerancePlus()
        """
        ...
    
    def setPrecursorMassToleranceMinus(self, tol: float ) -> None:
        """
        Cython signature: void setPrecursorMassToleranceMinus(double tol)
        """
        ...
    
    def getPrecursorMassToleranceMinus(self) -> float:
        """
        Cython signature: double getPrecursorMassToleranceMinus()
        """
        ...
    
    def setPrecursorErrorType(self, mono_isotopic: int ) -> None:
        """
        Cython signature: void setPrecursorErrorType(MassType mono_isotopic)
        """
        ...
    
    def getPrecursorErrorType(self) -> int:
        """
        Cython signature: MassType getPrecursorErrorType()
        """
        ...
    
    def setFragmentMassErrorUnit(self, unit: int ) -> None:
        """
        Cython signature: void setFragmentMassErrorUnit(ErrorUnit unit)
        """
        ...
    
    def getFragmentMassErrorUnit(self) -> int:
        """
        Cython signature: ErrorUnit getFragmentMassErrorUnit()
        """
        ...
    
    def setPrecursorMassErrorUnit(self, unit: int ) -> None:
        """
        Cython signature: void setPrecursorMassErrorUnit(ErrorUnit unit)
        """
        ...
    
    def getPrecursorMassErrorUnit(self) -> int:
        """
        Cython signature: ErrorUnit getPrecursorMassErrorUnit()
        """
        ...
    
    def setNumberOfThreads(self, threads: int ) -> None:
        """
        Cython signature: void setNumberOfThreads(unsigned int threads)
        """
        ...
    
    def getNumberOfThreads(self) -> int:
        """
        Cython signature: unsigned int getNumberOfThreads()
        """
        ...
    
    def setModifications(self, mods: ModificationDefinitionsSet ) -> None:
        """
        Cython signature: void setModifications(ModificationDefinitionsSet & mods)
        """
        ...
    
    def getModifications(self) -> ModificationDefinitionsSet:
        """
        Cython signature: ModificationDefinitionsSet getModifications()
        """
        ...
    
    def setOutputFilename(self, output: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setOutputFilename(const String & output)
        """
        ...
    
    def getOutputFilename(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOutputFilename()
        """
        ...
    
    def setInputFilename(self, input_file: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setInputFilename(const String & input_file)
        """
        ...
    
    def getInputFilename(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getInputFilename()
        """
        ...
    
    def setTaxonomyFilename(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTaxonomyFilename(const String & filename)
        """
        ...
    
    def getTaxonomyFilename(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTaxonomyFilename()
        """
        ...
    
    def setDefaultParametersFilename(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDefaultParametersFilename(const String & filename)
        """
        ...
    
    def getDefaultParametersFilename(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDefaultParametersFilename()
        """
        ...
    
    def setTaxon(self, taxon: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTaxon(const String & taxon)
        """
        ...
    
    def getTaxon(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getTaxon()
        """
        ...
    
    def setMaxPrecursorCharge(self, max_charge: int ) -> None:
        """
        Cython signature: void setMaxPrecursorCharge(int max_charge)
        """
        ...
    
    def getMaxPrecursorCharge(self) -> int:
        """
        Cython signature: int getMaxPrecursorCharge()
        """
        ...
    
    def setNumberOfMissedCleavages(self, missed_cleavages: int ) -> None:
        """
        Cython signature: void setNumberOfMissedCleavages(unsigned int missed_cleavages)
        """
        ...
    
    def getNumberOfMissedCleavages(self) -> int:
        """
        Cython signature: unsigned int getNumberOfMissedCleavages()
        """
        ...
    
    def setOutputResults(self, result: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setOutputResults(String result)
        """
        ...
    
    def getOutputResults(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOutputResults()
        """
        ...
    
    def setMaxValidEValue(self, value: float ) -> None:
        """
        Cython signature: void setMaxValidEValue(double value)
        """
        ...
    
    def getMaxValidEValue(self) -> float:
        """
        Cython signature: double getMaxValidEValue()
        """
        ...
    
    def setSemiCleavage(self, semi_cleavage: bool ) -> None:
        """
        Cython signature: void setSemiCleavage(bool semi_cleavage)
        """
        ...
    
    def setAllowIsotopeError(self, allow_isotope_error: bool ) -> None:
        """
        Cython signature: void setAllowIsotopeError(bool allow_isotope_error)
        """
        ...
    
    def write(self, filename: Union[bytes, str, String] , ignore_member_parameters: bool , force_default_mods: bool ) -> None:
        """
        Cython signature: void write(String filename, bool ignore_member_parameters, bool force_default_mods)
        """
        ...
    
    def setCleavageSite(self, cleavage_site: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setCleavageSite(String cleavage_site)
        """
        ...
    
    def getCleavageSite(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getCleavageSite()
        """
        ...
    ErrorUnit : __ErrorUnit
    MassType : __MassType 


class __ErrorUnit:
    None
    DALTONS : int
    PPM : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __MassType:
    None
    MONOISOTOPIC : int
    AVERAGE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

