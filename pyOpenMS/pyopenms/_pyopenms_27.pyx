#Generated with autowrap 0.23.0 and Cython (Parser) 3.1.2
#cython: c_string_encoding=ascii
#cython: embedsignature=False
from  enum            import Enum as _PyEnum
from  cpython         cimport Py_buffer
from  cpython         cimport bool as pybool_t
from  libcpp.string   cimport string as libcpp_string
from  libcpp.string   cimport string as libcpp_utf8_string
from  libcpp.string   cimport string as libcpp_utf8_output_string
from  libcpp.set      cimport set as libcpp_set
from  libcpp.vector   cimport vector as libcpp_vector
from  libcpp.pair     cimport pair as libcpp_pair
from  libcpp.map      cimport map  as libcpp_map
from  libcpp          cimport bool
from  libc.string     cimport const_char
from  cython.operator cimport dereference as deref, preincrement as inc, address as address
from  AutowrapRefHolder      cimport AutowrapRefHolder
from  AutowrapPtrHolder      cimport AutowrapPtrHolder
from  AutowrapConstPtrHolder cimport AutowrapConstPtrHolder
from  smart_ptr       cimport shared_ptr
cimport numpy as np
import numpy as np
cimport numpy as numpy
import numpy as numpy
def __static_TransformationModelBSpline_getDefaultParameters(Param params ):
    """
    __static_TransformationModelBSpline_getDefaultParameters(params: Param ) -> None
    """
    assert isinstance(params, Param), 'arg params wrong type'

    _getDefaultParameters_TransformationModelBSpline((deref(params.inst.get())))

def __static_ChromatogramExtractor_prepare_coordinates(list output_chromatograms , list extraction_coordinates , TargetedExperiment targeted , double rt_extraction_window , bool ms1 ,  ms1_isotopes ):
    """
    __static_ChromatogramExtractor_prepare_coordinates(output_chromatograms: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , targeted: TargetedExperiment , rt_extraction_window: float , ms1: bool , ms1_isotopes: int ) -> None
    """
    assert isinstance(output_chromatograms, list) and all(isinstance(elemt_rec, OSChromatogram) for elemt_rec in output_chromatograms), 'arg output_chromatograms wrong type'
    assert isinstance(extraction_coordinates, list) and all(isinstance(elemt_rec, ExtractionCoordinates) for elemt_rec in extraction_coordinates), 'arg extraction_coordinates wrong type'
    assert isinstance(targeted, TargetedExperiment), 'arg targeted wrong type'
    assert isinstance(rt_extraction_window, float), 'arg rt_extraction_window wrong type'
    assert isinstance(ms1, pybool_t), 'arg ms1 wrong type'
    assert isinstance(ms1_isotopes, int), 'arg ms1_isotopes wrong type'
    cdef libcpp_vector[shared_ptr[_OSChromatogram]] v0
    cdef OSChromatogram output_chromatograms_rec
    for output_chromatograms_rec in output_chromatograms:
        v0.push_back(output_chromatograms_rec.inst)
    # call
    cdef libcpp_vector[_ExtractionCoordinates] * v1 = new libcpp_vector[_ExtractionCoordinates]()
    cdef ExtractionCoordinates item1
    for item1 in extraction_coordinates:
        v1.push_back(deref(item1.inst.get()))




    _prepare_coordinates_ChromatogramExtractor(v0, deref(v1), (deref(targeted.inst.get())), (<double>rt_extraction_window), (<bool>ms1), (<int>ms1_isotopes))
    cdef libcpp_vector[_ExtractionCoordinates].iterator it_extraction_coordinates = v1.begin()
    replace_0 = []
    while it_extraction_coordinates != v1.end():
        item1 = ExtractionCoordinates.__new__(ExtractionCoordinates)
        item1.inst = shared_ptr[_ExtractionCoordinates](new _ExtractionCoordinates(deref(it_extraction_coordinates)))
        replace_0.append(item1)
        inc(it_extraction_coordinates)
    extraction_coordinates[:] = replace_0
    del v1
    # gather results
    replace = list()
    cdef libcpp_vector[shared_ptr[_OSChromatogram]].iterator it_output_chromatograms = v0.begin()
    while it_output_chromatograms != v0.end():
       output_chromatograms_rec = OSChromatogram.__new__(OSChromatogram)
       output_chromatograms_rec.inst = deref(it_output_chromatograms)
       replace.append(output_chromatograms_rec)
       inc(it_output_chromatograms)
    # replace old vector with new contents
    output_chromatograms[:] = []
    for output_chromatograms_rec_b in replace:
        output_chromatograms.append(output_chromatograms_rec_b) 

cdef class __ErrorUnit:
    None
    DALTONS = 0
    PPM = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __MassType:
    None
    MONOISOTOPIC = 0
    AVERAGE = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ChromatogramExtractor:
    """
    Cython implementation of _ChromatogramExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramExtractor.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramExtractor rv = ChromatogramExtractor.__new__(ChromatogramExtractor)
       rv.inst = shared_ptr[_ChromatogramExtractor](new _ChromatogramExtractor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramExtractor rv = ChromatogramExtractor.__new__(ChromatogramExtractor)
       rv.inst = shared_ptr[_ChromatogramExtractor](new _ChromatogramExtractor(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChromatogramExtractor](new _ChromatogramExtractor())
    
    def _init_1(self, ChromatogramExtractor in_0 ):
        """
        _init_1(self, in_0: ChromatogramExtractor ) -> None
        """
        assert isinstance(in_0, ChromatogramExtractor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramExtractor](new _ChromatogramExtractor((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramExtractor ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramExtractor)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractChromatograms(self, SpectrumAccessOpenMS input , list output , list extraction_coordinates , double mz_extraction_window , bool ppm , double im_extraction_window ,  filter ):
        """
        extractChromatograms(self, input: SpectrumAccessOpenMS , output: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , mz_extraction_window: float , ppm: bool , im_extraction_window: float , filter: Union[bytes, str, String] ) -> None
        """
        assert isinstance(input, SpectrumAccessOpenMS), 'arg input wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, OSChromatogram) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(extraction_coordinates, list) and all(isinstance(elemt_rec, ExtractionCoordinates) for elemt_rec in extraction_coordinates), 'arg extraction_coordinates wrong type'
        assert isinstance(mz_extraction_window, float), 'arg mz_extraction_window wrong type'
        assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
        assert isinstance(im_extraction_window, float), 'arg im_extraction_window wrong type'
        assert (isinstance(filter, str) or isinstance(filter, bytes) or isinstance(filter, String)), 'arg filter wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMS] input_input = input.inst
        cdef libcpp_vector[shared_ptr[_OSChromatogram]] v1
        cdef OSChromatogram output_rec
        for output_rec in output:
            v1.push_back(output_rec.inst)
        # call
        cdef libcpp_vector[_ExtractionCoordinates] * v2 = new libcpp_vector[_ExtractionCoordinates]()
        cdef ExtractionCoordinates item2
        for item2 in extraction_coordinates:
            v2.push_back(deref(item2.inst.get()))
    
    
    
    
        self.inst.get().extractChromatograms(input_input, v1, deref(v2), (<double>mz_extraction_window), (<bool>ppm), (<double>im_extraction_window), deref((convString(filter)).get()))
        del v2
        # gather results
        replace = list()
        cdef libcpp_vector[shared_ptr[_OSChromatogram]].iterator it_output = v1.begin()
        while it_output != v1.end():
           output_rec = OSChromatogram.__new__(OSChromatogram)
           output_rec.inst = deref(it_output)
           replace.append(output_rec)
           inc(it_output)
        # replace old vector with new contents
        output[:] = []
        for output_rec_b in replace:
            output.append(output_rec_b)
    
    def setLogType(self, int in_0 ):
        """
        setLogType(self, in_0: int ) -> None
        Sets the progress log that should be used. The default type is NONE!
        """
        assert in_0 in [0, 1, 2], 'arg in_0 wrong type'
    
        self.inst.get().setLogType((<_LogType>in_0))
    
    def getLogType(self):
        """
        getLogType(self) -> int
        Returns the type of progress log being used
        """
        cdef _LogType _r = self.inst.get().getLogType()
        py_result = <int>_r
        return py_result
    
    def startProgress(self,  begin ,  end ,  label ):
        """
        startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None
        """
        assert isinstance(begin, int) and begin >= 0, 'arg begin wrong type'
        assert isinstance(end, int) and end >= 0, 'arg end wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
    
    
        self.inst.get().startProgress((<ptrdiff_t>begin), (<ptrdiff_t>end), deref((convString(label)).get()))
    
    def setProgress(self,  value ):
        """
        setProgress(self, value: int ) -> None
        Sets the current progress
        """
        assert isinstance(value, int) and value >= 0, 'arg value wrong type'
    
        self.inst.get().setProgress((<ptrdiff_t>value))
    
    def endProgress(self):
        """
        endProgress(self) -> None
        Ends the progress display
        """
        self.inst.get().endProgress()
    
    def nextProgress(self):
        """
        nextProgress(self) -> None
        Increment progress by 1 (according to range begin-end)
        """
        self.inst.get().nextProgress()
    prepare_coordinates = __static_ChromatogramExtractor_prepare_coordinates 

cdef class DIAScoring:
    """
    Cython implementation of _DIAScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DIAScoring.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_DIAScoring](new _DIAScoring())
    
    def dia_ms1_massdiff_score(self, double precursor_mz , list spectrum , RangeMobility im_range , double ppm_score ):
        """
        dia_ms1_massdiff_score(self, precursor_mz: float , spectrum: List[OSSpectrum] , im_range: RangeMobility , ppm_score: float ) -> bool
        """
        assert isinstance(precursor_mz, float), 'arg precursor_mz wrong type'
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(im_range, RangeMobility), 'arg im_range wrong type'
        assert isinstance(ppm_score, float), 'arg ppm_score wrong type'
    
        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v1
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v1.push_back(spectrum_rec.inst)
        # call
    
    
        cdef bool _r = self.inst.get().dia_ms1_massdiff_score((<double>precursor_mz), v1, (deref(im_range.inst.get())), (<double &>ppm_score))
        
        py_result = <bool>_r
        return py_result
    
    def dia_ms1_isotope_scores_averagine(self, double precursor_mz , list spectrum ,  charge_state , RangeMobility im_range , double isotope_corr , double isotope_overlap ):
        """
        dia_ms1_isotope_scores_averagine(self, precursor_mz: float , spectrum: List[OSSpectrum] , charge_state: int , im_range: RangeMobility , isotope_corr: float , isotope_overlap: float ) -> None
        """
        assert isinstance(precursor_mz, float), 'arg precursor_mz wrong type'
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(charge_state, int), 'arg charge_state wrong type'
        assert isinstance(im_range, RangeMobility), 'arg im_range wrong type'
        assert isinstance(isotope_corr, float), 'arg isotope_corr wrong type'
        assert isinstance(isotope_overlap, float), 'arg isotope_overlap wrong type'
    
        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v1
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v1.push_back(spectrum_rec.inst)
        # call
    
    
    
    
        self.inst.get().dia_ms1_isotope_scores_averagine((<double>precursor_mz), v1, (<int>charge_state), (deref(im_range.inst.get())), (<double &>isotope_corr), (<double &>isotope_overlap))
        
    
    def dia_ms1_isotope_scores(self, double precursor_mz , list spectrum , RangeMobility im_range , double isotope_corr , double isotope_overlap , EmpiricalFormula sum_formula ):
        """
        dia_ms1_isotope_scores(self, precursor_mz: float , spectrum: List[OSSpectrum] , im_range: RangeMobility , isotope_corr: float , isotope_overlap: float , sum_formula: EmpiricalFormula ) -> None
        """
        assert isinstance(precursor_mz, float), 'arg precursor_mz wrong type'
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(im_range, RangeMobility), 'arg im_range wrong type'
        assert isinstance(isotope_corr, float), 'arg isotope_corr wrong type'
        assert isinstance(isotope_overlap, float), 'arg isotope_overlap wrong type'
        assert isinstance(sum_formula, EmpiricalFormula), 'arg sum_formula wrong type'
    
        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v1
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v1.push_back(spectrum_rec.inst)
        # call
    
    
    
    
        self.inst.get().dia_ms1_isotope_scores((<double>precursor_mz), v1, (deref(im_range.inst.get())), (<double &>isotope_corr), (<double &>isotope_overlap), (deref(sum_formula.inst.get())))
        
    
    def score_with_isotopes(self, list spectrum , list transitions , RangeMobility im_range , double dotprod , double manhattan ):
        """
        score_with_isotopes(self, spectrum: List[OSSpectrum] , transitions: List[LightTransition] , im_range: RangeMobility , dotprod: float , manhattan: float ) -> None
        """
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(transitions, list) and all(isinstance(elemt_rec, LightTransition) for elemt_rec in transitions), 'arg transitions wrong type'
        assert isinstance(im_range, RangeMobility), 'arg im_range wrong type'
        assert isinstance(dotprod, float), 'arg dotprod wrong type'
        assert isinstance(manhattan, float), 'arg manhattan wrong type'
        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v0
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v0.push_back(spectrum_rec.inst)
        # call
        cdef libcpp_vector[_LightTransition] * v1 = new libcpp_vector[_LightTransition]()
        cdef LightTransition item1
        for item1 in transitions:
            v1.push_back(deref(item1.inst.get()))
    
    
    
        self.inst.get().score_with_isotopes(v0, deref(v1), (deref(im_range.inst.get())), (<double &>dotprod), (<double &>manhattan))
        del v1
        
    
    def getSubsections(self):
        """
        getSubsections(self) -> List[bytes]
        """
        _r = self.inst.get().getSubsections()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setParameters(self, Param param ):
        """
        setParameters(self, param: Param ) -> None
        Sets the parameters
        """
        assert isinstance(param, Param), 'arg param wrong type'
    
        self.inst.get().setParameters((deref(param.inst.get())))
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        Returns the parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def getDefaults(self):
        """
        getDefaults(self) -> Param
        Returns the default parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getDefaults())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        Sets the name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get()))
    
    def dia_by_ion_score(self, list spectrum , AASequence sequence, charge, RangeMobility im_range, float bseries_score , float yseries_score):
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, OSSpectrum) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(sequence, AASequence), 'arg sequence wrong type'
        assert isinstance(im_range, RangeMobility), 'arg sequence wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(bseries_score, float), 'arg bseries_score wrong type'
        assert isinstance(yseries_score, float), 'arg yseries_score wrong type'

        cdef libcpp_vector[shared_ptr[_OSSpectrum]] v1
        cdef OSSpectrum spectrum_rec
        for spectrum_rec in spectrum:
            v1.push_back(spectrum_rec.inst)

        cdef double input_bseries_score = (<double>bseries_score)
        cdef double input_yseries_score = (<double>yseries_score)
        self.inst.get().dia_by_ion_score(v1, (deref(sequence.inst.get())), (<int>charge), (deref(im_range.inst.get())), input_bseries_score, input_yseries_score)
        yseries_score = input_yseries_score
        bseries_score = input_bseries_score
        return (bseries_score,yseries_score) 

cdef class MSDataStoringConsumer:
    """
    Cython implementation of _MSDataStoringConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataStoringConsumer.html>`_

    Consumer class that simply stores the data
    
    This class is able to keep spectra and chromatograms passed to it in memory
    and the data can be accessed through getData()
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSDataStoringConsumer rv = MSDataStoringConsumer.__new__(MSDataStoringConsumer)
       rv.inst = shared_ptr[_MSDataStoringConsumer](new _MSDataStoringConsumer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSDataStoringConsumer rv = MSDataStoringConsumer.__new__(MSDataStoringConsumer)
       rv.inst = shared_ptr[_MSDataStoringConsumer](new _MSDataStoringConsumer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSDataStoringConsumer](new _MSDataStoringConsumer())
    
    def _init_1(self, MSDataStoringConsumer in_0 ):
        """
        _init_1(self, in_0: MSDataStoringConsumer ) -> None
        """
        assert isinstance(in_0, MSDataStoringConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSDataStoringConsumer](new _MSDataStoringConsumer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSDataStoringConsumer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSDataStoringConsumer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        Sets experimental settings
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        Sets expected size
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms))
    
    def consumeSpectrum(self, MSSpectrum s ):
        """
        consumeSpectrum(self, s: MSSpectrum ) -> None
        """
        assert isinstance(s, MSSpectrum), 'arg s wrong type'
    
        self.inst.get().consumeSpectrum((deref(s.inst.get())))
    
    def consumeChromatogram(self, MSChromatogram in_0 ):
        """
        consumeChromatogram(self, in_0: MSChromatogram ) -> None
        """
        assert isinstance(in_0, MSChromatogram), 'arg in_0 wrong type'
    
        self.inst.get().consumeChromatogram((deref(in_0.inst.get())))
    
    def getData(self):
        """
        getData(self) -> MSExperiment
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getData())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result 

cdef class MasstraceCorrelator:
    """
    Cython implementation of _MasstraceCorrelator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MasstraceCorrelator.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MasstraceCorrelator rv = MasstraceCorrelator.__new__(MasstraceCorrelator)
       rv.inst = shared_ptr[_MasstraceCorrelator](new _MasstraceCorrelator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MasstraceCorrelator rv = MasstraceCorrelator.__new__(MasstraceCorrelator)
       rv.inst = shared_ptr[_MasstraceCorrelator](new _MasstraceCorrelator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MasstraceCorrelator](new _MasstraceCorrelator())
    
    def _init_1(self, MasstraceCorrelator in_0 ):
        """
        _init_1(self, in_0: MasstraceCorrelator ) -> None
        """
        assert isinstance(in_0, MasstraceCorrelator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MasstraceCorrelator](new _MasstraceCorrelator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MasstraceCorrelator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MasstraceCorrelator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def createPseudoSpectra(self, ConsensusMap map_ , MSExperiment pseudo_spectra ,  min_peak_nr , double min_correlation ,  max_lag , double max_rt_apex_difference ):
        """
        createPseudoSpectra(self, map_: ConsensusMap , pseudo_spectra: MSExperiment , min_peak_nr: int , min_correlation: float , max_lag: int , max_rt_apex_difference: float ) -> None
        Compute pseudo-spectra from a set of (MS2) masstraces
        
        This function will take a set of masstraces (consensus map) as input and
        produce a vector of pseudo spectra as output (pseudo_spectra result
        vector).
        
        It basically makes an all-vs-all comparison of all masstraces against
        each other and scores them on how similar they are in their mass traces.
        
        This assumes that the consensus feature is only from one (SWATH) map
        This assumes that the consensus map is sorted by intensity
        """
        assert isinstance(map_, ConsensusMap), 'arg map_ wrong type'
        assert isinstance(pseudo_spectra, MSExperiment), 'arg pseudo_spectra wrong type'
        assert isinstance(min_peak_nr, int) and min_peak_nr >= 0, 'arg min_peak_nr wrong type'
        assert isinstance(min_correlation, float), 'arg min_correlation wrong type'
        assert isinstance(max_lag, int), 'arg max_lag wrong type'
        assert isinstance(max_rt_apex_difference, float), 'arg max_rt_apex_difference wrong type'
    
    
    
    
    
    
        self.inst.get().createPseudoSpectra((deref(map_.inst.get())), (deref(pseudo_spectra.inst.get())), (<size_t>min_peak_nr), (<double>min_correlation), (<int>max_lag), (<double>max_rt_apex_difference))
    
    def getSubsections(self):
        """
        getSubsections(self) -> List[bytes]
        """
        _r = self.inst.get().getSubsections()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setParameters(self, Param param ):
        """
        setParameters(self, param: Param ) -> None
        Sets the parameters
        """
        assert isinstance(param, Param), 'arg param wrong type'
    
        self.inst.get().setParameters((deref(param.inst.get())))
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        Returns the parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def getDefaults(self):
        """
        getDefaults(self) -> Param
        Returns the default parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getDefaults())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        Sets the name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get()))
    
    def setLogType(self, int in_0 ):
        """
        setLogType(self, in_0: int ) -> None
        Sets the progress log that should be used. The default type is NONE!
        """
        assert in_0 in [0, 1, 2], 'arg in_0 wrong type'
    
        self.inst.get().setLogType((<_LogType>in_0))
    
    def getLogType(self):
        """
        getLogType(self) -> int
        Returns the type of progress log being used
        """
        cdef _LogType _r = self.inst.get().getLogType()
        py_result = <int>_r
        return py_result
    
    def startProgress(self,  begin ,  end ,  label ):
        """
        startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None
        """
        assert isinstance(begin, int) and begin >= 0, 'arg begin wrong type'
        assert isinstance(end, int) and end >= 0, 'arg end wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
    
    
        self.inst.get().startProgress((<ptrdiff_t>begin), (<ptrdiff_t>end), deref((convString(label)).get()))
    
    def setProgress(self,  value ):
        """
        setProgress(self, value: int ) -> None
        Sets the current progress
        """
        assert isinstance(value, int) and value >= 0, 'arg value wrong type'
    
        self.inst.get().setProgress((<ptrdiff_t>value))
    
    def endProgress(self):
        """
        endProgress(self) -> None
        Ends the progress display
        """
        self.inst.get().endProgress()
    
    def nextProgress(self):
        """
        nextProgress(self) -> None
        Increment progress by 1 (according to range begin-end)
        """
        self.inst.get().nextProgress() 

cdef class MatrixDouble:
    """
    Cython implementation of _Matrix[double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Matrix[double].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Matrix[double]](new _Matrix[double]())
    
    def _init_1(self, MatrixDouble in_0 ):
        """
        _init_1(self, in_0: MatrixDouble ) -> None
        """
        assert isinstance(in_0, MatrixDouble), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Matrix[double]](new _Matrix[double]((deref(in_0.inst.get()))))
    
    def _init_2(self,  rows ,  cols , double value ):
        """
        _init_2(self, rows: int , cols: int , value: float ) -> None
        """
        assert isinstance(rows, int) and rows >= 0, 'arg rows wrong type'
        assert isinstance(cols, int) and cols >= 0, 'arg cols wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
    
        self.inst = shared_ptr[_Matrix[double]](new _Matrix[double]((<size_t>rows), (<size_t>cols), (<double>value)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MatrixDouble ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, rows: int , cols: int , value: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MatrixDouble)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getValue(self,  i ,  j ):
        """
        getValue(self, i: int , j: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
    
    
        cdef double _r = self.inst.get().getValue((<size_t>i), (<size_t>j))
        py_result = <double>_r
        return py_result
    
    def setValue(self,  i ,  j , double value ):
        """
        setValue(self, i: int , j: int , value: float ) -> None
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
    
        self.inst.get().setValue((<size_t>i), (<size_t>j), (<double>value))
    
    def rows(self):
        """
        rows(self) -> int
        """
        cdef size_t _r = self.inst.get().rows()
        py_result = <size_t>_r
        return py_result
    
    def cols(self):
        """
        cols(self) -> int
        """
        cdef size_t _r = self.inst.get().cols()
        py_result = <size_t>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def resize(self,  rows ,  cols ):
        """
        resize(self, rows: int , cols: int ) -> None
        """
        assert isinstance(rows, int) and rows >= 0, 'arg rows wrong type'
        assert isinstance(cols, int) and cols >= 0, 'arg cols wrong type'
    
    
        self.inst.get().resize((<size_t>rows), (<size_t>cols))
    
# continue with extra code if needed
            

    
    @staticmethod
    def fromNdArray(np.ndarray[double, ndim=2] data not None):
        """Creates a new Matrix from a numpy ndarray."""
        cdef MatrixDouble mat = MatrixDouble()
        mat.set_matrix(data)
        return mat

    def get_matrix_as_view(self):
        """get_matrix(self) -> np.ndarray[double, ndim=2]

        Returns a view on the underlying Matrix as a 2D numpy ndarray.
        .. caution::
           Future changes to the Matrix will affect the ndarray and vice versa.
           Make sure that the Matrix does not go out of scope before the last use
           of your ndarray.
        """

        cdef _Matrix[double] * mat_ = self.inst.get()
        cdef unsigned int rows = mat_.rows()
        cdef unsigned int cols = mat_.cols()
        cdef double* data = mat_.data()
        cdef double[:,:] mem_view = <double[:rows,:cols]>data
        dtype = 'double'
        cdef int itemsize = np.dtype(dtype).itemsize
        cdef unsigned int row_stride, col_stride
        o = 'F'
        if mat_.rowMajor():
            row_stride = mat_.outerStride() if mat_.outerStride() > 0 else cols
            col_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            o = 'F'
        else:
            row_stride = mat_.innerStride() if mat_.innerStride() > 0 else 1
            col_stride = mat_.outerStride() if mat_.outerStride() > 0 else rows
            o = 'C'

        return np.lib.stride_tricks.as_strided(np.asarray(mem_view, dtype=dtype, order=o), strides=[row_stride*itemsize, col_stride*itemsize])


    def get_matrix(self):
        """get_matrix(self) -> np.ndarray[double, ndim=2]

        Returns a copy of the underlying Matrix as a 2D numpy ndarray.
        """

        return np.copy(self.get_matrix_as_view())


    def set_matrix(self, np.ndarray[double, ndim=2] data not None):
        """set_matrix(self, data: np.ndarray[double, ndim=2]) -> None

        Copies the values from the numpy ndarray into the Matrix.
        """

        cdef _Matrix[double] * mat_ = self.inst.get()

        cdef unsigned int rows = data.shape[0]
        cdef unsigned int cols = data.shape[1]
        mat_.resize(rows, cols)

        cdef int i = 0
        cdef int j = 0
        for i in range(int(rows)):
            for j in range(int(cols)):
                mat_.setValue(i, j, data[i][j]) 

cdef class NoopMSDataWritingConsumer:
    """
    Cython implementation of _NoopMSDataWritingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NoopMSDataWritingConsumer.html>`_

    Consumer class that perform no operation
    
    This is sometimes necessary to fulfill the requirement of passing an
    valid MSDataWritingConsumer object or pointer but no operation is
    required
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self,  filename ):
        """
        __init__(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_NoopMSDataWritingConsumer](new _NoopMSDataWritingConsumer(deref((convString(filename)).get())))
    
    def consumeSpectrum(self, MSSpectrum s ):
        """
        consumeSpectrum(self, s: MSSpectrum ) -> None
        """
        assert isinstance(s, MSSpectrum), 'arg s wrong type'
    
        self.inst.get().consumeSpectrum((deref(s.inst.get())))
    
    def consumeChromatogram(self, MSChromatogram c ):
        """
        consumeChromatogram(self, c: MSChromatogram ) -> None
        """
        assert isinstance(c, MSChromatogram), 'arg c wrong type'
    
        self.inst.get().consumeChromatogram((deref(c.inst.get())))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms))
    
    def addDataProcessing(self, DataProcessing d ):
        """
        addDataProcessing(self, d: DataProcessing ) -> None
        """
        assert isinstance(d, DataProcessing), 'arg d wrong type'
    
        self.inst.get().addDataProcessing((deref(d.inst.get())))
    
    def getNrSpectraWritten(self):
        """
        getNrSpectraWritten(self) -> int
        """
        cdef size_t _r = self.inst.get().getNrSpectraWritten()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatogramsWritten(self):
        """
        getNrChromatogramsWritten(self) -> int
        """
        cdef size_t _r = self.inst.get().getNrChromatogramsWritten()
        py_result = <size_t>_r
        return py_result 

cdef class PlainMSDataWritingConsumer:
    """
    Cython implementation of _PlainMSDataWritingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PlainMSDataWritingConsumer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self,  filename ):
        """
        __init__(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_PlainMSDataWritingConsumer](new _PlainMSDataWritingConsumer(deref((convString(filename)).get())))
    
    def consumeSpectrum(self, MSSpectrum s ):
        """
        consumeSpectrum(self, s: MSSpectrum ) -> None
        """
        assert isinstance(s, MSSpectrum), 'arg s wrong type'
    
        self.inst.get().consumeSpectrum((deref(s.inst.get())))
    
    def consumeChromatogram(self, MSChromatogram c ):
        """
        consumeChromatogram(self, c: MSChromatogram ) -> None
        """
        assert isinstance(c, MSChromatogram), 'arg c wrong type'
    
        self.inst.get().consumeChromatogram((deref(c.inst.get())))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        Set experimental settings for the whole file
        
        
        :param exp: Experimental settings to be used for this file (from this and the first spectrum/chromatogram, the class will deduce most of the header of the mzML file)
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        Set expected size of spectra and chromatograms to be written
        
        These numbers will be written in the spectrumList and chromatogramList
        tag in the mzML file. Therefore, these will contain wrong numbers if
        the expected size is not set correctly
        
        
        :param expectedSpectra: Number of spectra expected
        :param expectedChromatograms: Number of chromatograms expected
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms))
    
    def addDataProcessing(self, DataProcessing d ):
        """
        addDataProcessing(self, d: DataProcessing ) -> None
        Optionally add a data processing method to each chromatogram and spectrum
        
        The provided DataProcessing object will be added to each chromatogram
        and spectrum written to to the mzML file
        
        
        :param d: The DataProcessing object to be added
        """
        assert isinstance(d, DataProcessing), 'arg d wrong type'
    
        self.inst.get().addDataProcessing((deref(d.inst.get())))
    
    def getNrSpectraWritten(self):
        """
        getNrSpectraWritten(self) -> int
        Returns the number of spectra written
        """
        cdef size_t _r = self.inst.get().getNrSpectraWritten()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatogramsWritten(self):
        """
        getNrChromatogramsWritten(self) -> int
        Returns the number of chromatograms written
        """
        cdef size_t _r = self.inst.get().getNrChromatogramsWritten()
        py_result = <size_t>_r
        return py_result
    
    def setOptions(self, PeakFileOptions opt ):
        """
        setOptions(self, opt: PeakFileOptions ) -> None
        """
        assert isinstance(opt, PeakFileOptions), 'arg opt wrong type'
    
        self.inst.get().setOptions((deref(opt.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
        return py_result 

cdef class QTCluster:
    """
    Cython implementation of _QTCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QTCluster.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef QTCluster rv = QTCluster.__new__(QTCluster)
       rv.inst = shared_ptr[_QTCluster](new _QTCluster(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef QTCluster rv = QTCluster.__new__(QTCluster)
       rv.inst = shared_ptr[_QTCluster](new _QTCluster(deref(self.inst.get())))
       return rv
    
    def __init__(self, QTCluster in_0 ):
        """
        __init__(self, in_0: QTCluster ) -> None
        """
        assert isinstance(in_0, QTCluster), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_QTCluster](new _QTCluster((deref(in_0.inst.get()))))
    
    def getCenterRT(self):
        """
        getCenterRT(self) -> float
        Returns the RT value of the cluster
        """
        cdef double _r = self.inst.get().getCenterRT()
        py_result = <double>_r
        return py_result
    
    def getCenterMZ(self):
        """
        getCenterMZ(self) -> float
        Returns the m/z value of the cluster center
        """
        cdef double _r = self.inst.get().getCenterMZ()
        py_result = <double>_r
        return py_result
    
    def getXCoord(self):
        """
        getXCoord(self) -> int
        Returns the x coordinate in the grid
        """
        cdef int _r = self.inst.get().getXCoord()
        py_result = <int>_r
        return py_result
    
    def getYCoord(self):
        """
        getYCoord(self) -> int
        Returns the y coordinate in the grid
        """
        cdef int _r = self.inst.get().getYCoord()
        py_result = <int>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        Returns the size of the cluster (number of elements, incl. center)
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def getQuality(self):
        """
        getQuality(self) -> float
        Returns the cluster quality and recomputes if necessary
        """
        cdef double _r = self.inst.get().getQuality()
        py_result = <double>_r
        return py_result
    
    def getAnnotations(self):
        """
        getAnnotations(self) -> Set[AASequence]
        Returns the set of peptide sequences annotated to the cluster center
        """
        _r = self.inst.get().getAnnotations()
        py_result = set()
        cdef libcpp_set[_AASequence].iterator it__r = _r.begin()
        cdef AASequence item_py_result
        while it__r != _r.end():
           item_py_result = AASequence.__new__(AASequence)
           item_py_result.inst = shared_ptr[_AASequence](new _AASequence(deref(it__r)))
           py_result.add(item_py_result)
           inc(it__r)
        return py_result
    
    def setInvalid(self):
        """
        setInvalid(self) -> None
        Sets current cluster as invalid (also frees some memory)
        """
        self.inst.get().setInvalid()
    
    def isInvalid(self):
        """
        isInvalid(self) -> bool
        Whether current cluster is invalid
        """
        cdef bool _r = self.inst.get().isInvalid()
        py_result = <bool>_r
        return py_result
    
    def initializeCluster(self):
        """
        initializeCluster(self) -> None
        Has to be called before adding elements (calling
        """
        self.inst.get().initializeCluster()
    
    def finalizeCluster(self):
        """
        finalizeCluster(self) -> None
        Has to be called after adding elements (after calling
        """
        self.inst.get().finalizeCluster()
    
    def __richcmp__(self, other, op):
        if op not in (0,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, QTCluster):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef QTCluster other_casted = other
        cdef QTCluster self_casted = self
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class SwathFile:
    """
    Cython implementation of _SwathFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SwathFile rv = SwathFile.__new__(SwathFile)
       rv.inst = shared_ptr[_SwathFile](new _SwathFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SwathFile rv = SwathFile.__new__(SwathFile)
       rv.inst = shared_ptr[_SwathFile](new _SwathFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SwathFile](new _SwathFile())
    
    def _init_1(self, SwathFile in_0 ):
        """
        _init_1(self, in_0: SwathFile ) -> None
        """
        assert isinstance(in_0, SwathFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SwathFile](new _SwathFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SwathFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SwathFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def loadSplit(self, list file_list ,  tmp , ExperimentalSettings exp_meta ,  readoptions ):
        """
        loadSplit(self, file_list: List[bytes] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]
        Loads a Swath run from a list of split mzML files
        """
        assert isinstance(file_list, list) and all(isinstance(li, bytes) for li in file_list), 'arg file_list wrong type'
        assert (isinstance(tmp, str) or isinstance(tmp, bytes) or isinstance(tmp, String)), 'arg tmp wrong type'
        assert isinstance(exp_meta, ExperimentalSettings), 'arg exp_meta wrong type'
        assert (isinstance(readoptions, str) or isinstance(readoptions, bytes) or isinstance(readoptions, String)), 'arg readoptions wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in file_list:
           v0.push_back(_String(<char *>item0))
    
        cdef shared_ptr[_ExperimentalSettings] input_exp_meta = exp_meta.inst
    
        _r = self.inst.get().loadSplit(deref(v0), deref((convString(tmp)).get()), input_exp_meta, deref((convString(readoptions)).get()))
        del v0
        py_result = []
        cdef libcpp_vector[_SwathMap].iterator it__r = _r.begin()
        cdef SwathMap item_py_result
        while it__r != _r.end():
           item_py_result = SwathMap.__new__(SwathMap)
           item_py_result.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def loadMzML(self,  file_ ,  tmp , ExperimentalSettings exp_meta ,  readoptions ):
        """
        loadMzML(self, file_: Union[bytes, str, String] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]
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
        assert (isinstance(file_, str) or isinstance(file_, bytes) or isinstance(file_, String)), 'arg file_ wrong type'
        assert (isinstance(tmp, str) or isinstance(tmp, bytes) or isinstance(tmp, String)), 'arg tmp wrong type'
        assert isinstance(exp_meta, ExperimentalSettings), 'arg exp_meta wrong type'
        assert (isinstance(readoptions, str) or isinstance(readoptions, bytes) or isinstance(readoptions, String)), 'arg readoptions wrong type'
    
    
        cdef shared_ptr[_ExperimentalSettings] input_exp_meta = exp_meta.inst
    
        _r = self.inst.get().loadMzML(deref((convString(file_)).get()), deref((convString(tmp)).get()), input_exp_meta, deref((convString(readoptions)).get()))
        py_result = []
        cdef libcpp_vector[_SwathMap].iterator it__r = _r.begin()
        cdef SwathMap item_py_result
        while it__r != _r.end():
           item_py_result = SwathMap.__new__(SwathMap)
           item_py_result.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def loadMzXML(self,  file_ ,  tmp , ExperimentalSettings exp_meta ,  readoptions ):
        """
        loadMzXML(self, file_: Union[bytes, str, String] , tmp: Union[bytes, str, String] , exp_meta: ExperimentalSettings , readoptions: Union[bytes, str, String] ) -> List[SwathMap]
        Loads a Swath run from a single mzXML file
        """
        assert (isinstance(file_, str) or isinstance(file_, bytes) or isinstance(file_, String)), 'arg file_ wrong type'
        assert (isinstance(tmp, str) or isinstance(tmp, bytes) or isinstance(tmp, String)), 'arg tmp wrong type'
        assert isinstance(exp_meta, ExperimentalSettings), 'arg exp_meta wrong type'
        assert (isinstance(readoptions, str) or isinstance(readoptions, bytes) or isinstance(readoptions, String)), 'arg readoptions wrong type'
    
    
        cdef shared_ptr[_ExperimentalSettings] input_exp_meta = exp_meta.inst
    
        _r = self.inst.get().loadMzXML(deref((convString(file_)).get()), deref((convString(tmp)).get()), input_exp_meta, deref((convString(readoptions)).get()))
        py_result = []
        cdef libcpp_vector[_SwathMap].iterator it__r = _r.begin()
        cdef SwathMap item_py_result
        while it__r != _r.end():
           item_py_result = SwathMap.__new__(SwathMap)
           item_py_result.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setLogType(self, int in_0 ):
        """
        setLogType(self, in_0: int ) -> None
        Sets the progress log that should be used. The default type is NONE!
        """
        assert in_0 in [0, 1, 2], 'arg in_0 wrong type'
    
        self.inst.get().setLogType((<_LogType>in_0))
    
    def getLogType(self):
        """
        getLogType(self) -> int
        Returns the type of progress log being used
        """
        cdef _LogType _r = self.inst.get().getLogType()
        py_result = <int>_r
        return py_result
    
    def startProgress(self,  begin ,  end ,  label ):
        """
        startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None
        """
        assert isinstance(begin, int) and begin >= 0, 'arg begin wrong type'
        assert isinstance(end, int) and end >= 0, 'arg end wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
    
    
        self.inst.get().startProgress((<ptrdiff_t>begin), (<ptrdiff_t>end), deref((convString(label)).get()))
    
    def setProgress(self,  value ):
        """
        setProgress(self, value: int ) -> None
        Sets the current progress
        """
        assert isinstance(value, int) and value >= 0, 'arg value wrong type'
    
        self.inst.get().setProgress((<ptrdiff_t>value))
    
    def endProgress(self):
        """
        endProgress(self) -> None
        Ends the progress display
        """
        self.inst.get().endProgress()
    
    def nextProgress(self):
        """
        nextProgress(self) -> None
        Increment progress by 1 (according to range begin-end)
        """
        self.inst.get().nextProgress() 

cdef class TMTElevenPlexQuantitationMethod:
    """
    Cython implementation of _TMTElevenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTElevenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TMTElevenPlexQuantitationMethod rv = TMTElevenPlexQuantitationMethod.__new__(TMTElevenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTElevenPlexQuantitationMethod](new _TMTElevenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TMTElevenPlexQuantitationMethod rv = TMTElevenPlexQuantitationMethod.__new__(TMTElevenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTElevenPlexQuantitationMethod](new _TMTElevenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TMTElevenPlexQuantitationMethod](new _TMTElevenPlexQuantitationMethod())
    
    def _init_1(self, TMTElevenPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: TMTElevenPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, TMTElevenPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TMTElevenPlexQuantitationMethod](new _TMTElevenPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TMTElevenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTElevenPlexQuantitationMethod)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def getChannelInformation(self):
        """
        getChannelInformation(self) -> List[IsobaricChannelInformation]
        """
        _r = self.inst.get().getChannelInformation()
        py_result = []
        cdef libcpp_vector[_IsobaricChannelInformation].iterator it__r = _r.begin()
        cdef IsobaricChannelInformation item_py_result
        while it__r != _r.end():
           item_py_result = IsobaricChannelInformation.__new__(IsobaricChannelInformation)
           item_py_result.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getNumberOfChannels(self):
        """
        getNumberOfChannels(self) -> int
        """
        cdef size_t _r = self.inst.get().getNumberOfChannels()
        py_result = <size_t>_r
        return py_result
    
    def getIsotopeCorrectionMatrix(self):
        """
        getIsotopeCorrectionMatrix(self) -> MatrixDouble
        """
        cdef _Matrix[double] * _r = new _Matrix[double](self.inst.get().getIsotopeCorrectionMatrix())
        cdef MatrixDouble py_result = MatrixDouble.__new__(MatrixDouble)
        py_result.inst = shared_ptr[_Matrix[double]](_r)
        return py_result
    
    def getReferenceChannel(self):
        """
        getReferenceChannel(self) -> int
        """
        cdef size_t _r = self.inst.get().getReferenceChannel()
        py_result = <size_t>_r
        return py_result
    
    def getSubsections(self):
        """
        getSubsections(self) -> List[bytes]
        """
        _r = self.inst.get().getSubsections()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setParameters(self, Param param ):
        """
        setParameters(self, param: Param ) -> None
        Sets the parameters
        """
        assert isinstance(param, Param), 'arg param wrong type'
    
        self.inst.get().setParameters((deref(param.inst.get())))
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        Returns the parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def getDefaults(self):
        """
        getDefaults(self) -> Param
        Returns the default parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getDefaults())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        Sets the name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get())) 

cdef class TransformationModelBSpline:
    """
    Cython implementation of _TransformationModelBSpline

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelBSpline.html>`_
      -- Inherits from ['TransformationModel']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list data , Param params ):
        """
        __init__(self, data: List[TM_DataPoint] , params: Param ) -> None
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        assert isinstance(params, Param), 'arg params wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_TransformationModelBSpline](new _TransformationModelBSpline(deref(v0), (deref(params.inst.get()))))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        Gets the default parameters
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def evaluate(self, double value ):
        """
        evaluate(self, value: float ) -> float
        Evaluates the model at the given values
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef double _r = self.inst.get().evaluate((<double>value))
        py_result = <double>_r
        return py_result
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def weightData(self, list data ):
        """
        weightData(self, data: List[TM_DataPoint] ) -> None
        Weight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().weightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkValidWeight(self,  weight , list valid_weights ):
        """
        checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool
        Check for a valid weighting function string
        """
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
        assert isinstance(valid_weights, list) and all(isinstance(i, bytes) for i in valid_weights), 'arg valid_weights wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in valid_weights:
           v1.push_back(_String(<char *>item1))
        cdef bool _r = self.inst.get().checkValidWeight(deref((convString(weight)).get()), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        valid_weights[:] = replace
        del v1
        py_result = <bool>_r
        return py_result
    
    def weightDatum(self, double datum ,  weight ):
        """
        weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Weight the data according to the weighting function
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().weightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def unWeightDatum(self, double datum ,  weight ):
        """
        unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Apply the reverse of the weighting function to the data
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().unWeightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def getValidXWeights(self):
        """
        getValidXWeights(self) -> List[bytes]
        Returns a list of valid x weight function stringss
        """
        _r = self.inst.get().getValidXWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getValidYWeights(self):
        """
        getValidYWeights(self) -> List[bytes]
        Returns a list of valid y weight function strings
        """
        _r = self.inst.get().getValidYWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def unWeightData(self, list data ):
        """
        unWeightData(self, data: List[TM_DataPoint] ) -> None
        Unweight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().unWeightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkDatumRange(self, double datum , double datum_min , double datum_max ):
        """
        checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float
        Check that the datum is within the valid min and max bounds
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert isinstance(datum_min, float), 'arg datum_min wrong type'
        assert isinstance(datum_max, float), 'arg datum_max wrong type'
    
    
    
        cdef double _r = self.inst.get().checkDatumRange((<const double &>datum), (<const double &>datum_min), (<const double &>datum_max))
        py_result = <double>_r
        return py_result
    getDefaultParameters = __static_TransformationModelBSpline_getDefaultParameters 

cdef class TransitionTSVFile:
    """
    Cython implementation of _TransitionTSVFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransitionTSVFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TransitionTSVFile rv = TransitionTSVFile.__new__(TransitionTSVFile)
       rv.inst = shared_ptr[_TransitionTSVFile](new _TransitionTSVFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TransitionTSVFile rv = TransitionTSVFile.__new__(TransitionTSVFile)
       rv.inst = shared_ptr[_TransitionTSVFile](new _TransitionTSVFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TransitionTSVFile](new _TransitionTSVFile())
    
    def _init_1(self, TransitionTSVFile in_0 ):
        """
        _init_1(self, in_0: TransitionTSVFile ) -> None
        """
        assert isinstance(in_0, TransitionTSVFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TransitionTSVFile](new _TransitionTSVFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TransitionTSVFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TransitionTSVFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def convertTargetedExperimentToTSV(self, bytes filename , TargetedExperiment targeted_exp ):
        """
        convertTargetedExperimentToTSV(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None
        Write out a targeted experiment (TraML structure) into a tsv file
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
    
        self.inst.get().convertTargetedExperimentToTSV((<char *>filename), (deref(targeted_exp.inst.get())))
    
    def _convertTSVToTargetedExperiment_0(self, bytes filename , int filetype , TargetedExperiment targeted_exp ):
        """
        _convertTSVToTargetedExperiment_0(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None
        Read in a tsv/mrm file and construct a targeted experiment (TraML structure)
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert filetype in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg filetype wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
    
    
        self.inst.get().convertTSVToTargetedExperiment((<char *>filename), (<_FileType>filetype), (deref(targeted_exp.inst.get())))
    
    def _convertTSVToTargetedExperiment_1(self, bytes filename , int filetype , LightTargetedExperiment targeted_exp ):
        """
        _convertTSVToTargetedExperiment_1(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None
        Read in a tsv file and construct a targeted experiment (Light transition structure)
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert filetype in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg filetype wrong type'
        assert isinstance(targeted_exp, LightTargetedExperiment), 'arg targeted_exp wrong type'
    
    
    
        self.inst.get().convertTSVToTargetedExperiment((<char *>filename), (<_FileType>filetype), (deref(targeted_exp.inst.get())))
    
    def convertTSVToTargetedExperiment(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None
          :noindex:
        
        Read in a tsv/mrm file and construct a targeted experiment (TraML structure)

        
        .. rubric:: Overload:
        .. py:function:: convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None
          :noindex:
        
        Read in a tsv file and construct a targeted experiment (Light transition structure)
    
        """
        if (len(args)==3) and (isinstance(args[0], bytes)) and (args[1] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]) and (isinstance(args[2], TargetedExperiment)):
            return self._convertTSVToTargetedExperiment_0(*args)
        elif (len(args)==3) and (isinstance(args[0], bytes)) and (args[1] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]) and (isinstance(args[2], LightTargetedExperiment)):
            return self._convertTSVToTargetedExperiment_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def validateTargetedExperiment(self, TargetedExperiment targeted_exp ):
        """
        validateTargetedExperiment(self, targeted_exp: TargetedExperiment ) -> None
        Validate a TargetedExperiment (check that all ids are unique)
        """
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
        self.inst.get().validateTargetedExperiment((deref(targeted_exp.inst.get())))
    
    def setLogType(self, int in_0 ):
        """
        setLogType(self, in_0: int ) -> None
        Sets the progress log that should be used. The default type is NONE!
        """
        assert in_0 in [0, 1, 2], 'arg in_0 wrong type'
    
        self.inst.get().setLogType((<_LogType>in_0))
    
    def getLogType(self):
        """
        getLogType(self) -> int
        Returns the type of progress log being used
        """
        cdef _LogType _r = self.inst.get().getLogType()
        py_result = <int>_r
        return py_result
    
    def startProgress(self,  begin ,  end ,  label ):
        """
        startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None
        """
        assert isinstance(begin, int) and begin >= 0, 'arg begin wrong type'
        assert isinstance(end, int) and end >= 0, 'arg end wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
    
    
        self.inst.get().startProgress((<ptrdiff_t>begin), (<ptrdiff_t>end), deref((convString(label)).get()))
    
    def setProgress(self,  value ):
        """
        setProgress(self, value: int ) -> None
        Sets the current progress
        """
        assert isinstance(value, int) and value >= 0, 'arg value wrong type'
    
        self.inst.get().setProgress((<ptrdiff_t>value))
    
    def endProgress(self):
        """
        endProgress(self) -> None
        Ends the progress display
        """
        self.inst.get().endProgress()
    
    def nextProgress(self):
        """
        nextProgress(self) -> None
        Increment progress by 1 (according to range begin-end)
        """
        self.inst.get().nextProgress() 

cdef class UniqueIdGenerator:
    """
    Cython implementation of _UniqueIdGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1UniqueIdGenerator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def getUniqueId(self):
        """
        getUniqueId(self) -> int
        """
        cdef uint64_t _r = self.inst.get().getUniqueId()
        py_result = <uint64_t>_r
        return py_result
    
    def setSeed(self,  in_0 ):
        """
        setSeed(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        self.inst.get().setSeed((<uint64_t>in_0))
    
    def getSeed(self):
        """
        getSeed(self) -> int
        """
        cdef uint64_t _r = self.inst.get().getSeed()
        py_result = <uint64_t>_r
        return py_result 

cdef class XTandemInfile:
    """
    Cython implementation of _XTandemInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XTandemInfile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_XTandemInfile](new _XTandemInfile())
    
    def setFragmentMassTolerance(self, double tolerance ):
        """
        setFragmentMassTolerance(self, tolerance: float ) -> None
        """
        assert isinstance(tolerance, float), 'arg tolerance wrong type'
    
        self.inst.get().setFragmentMassTolerance((<double>tolerance))
    
    def getFragmentMassTolerance(self):
        """
        getFragmentMassTolerance(self) -> float
        """
        cdef double _r = self.inst.get().getFragmentMassTolerance()
        py_result = <double>_r
        return py_result
    
    def setPrecursorMassTolerancePlus(self, double tol ):
        """
        setPrecursorMassTolerancePlus(self, tol: float ) -> None
        """
        assert isinstance(tol, float), 'arg tol wrong type'
    
        self.inst.get().setPrecursorMassTolerancePlus((<double>tol))
    
    def getPrecursorMassTolerancePlus(self):
        """
        getPrecursorMassTolerancePlus(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorMassTolerancePlus()
        py_result = <double>_r
        return py_result
    
    def setPrecursorMassToleranceMinus(self, double tol ):
        """
        setPrecursorMassToleranceMinus(self, tol: float ) -> None
        """
        assert isinstance(tol, float), 'arg tol wrong type'
    
        self.inst.get().setPrecursorMassToleranceMinus((<double>tol))
    
    def getPrecursorMassToleranceMinus(self):
        """
        getPrecursorMassToleranceMinus(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorMassToleranceMinus()
        py_result = <double>_r
        return py_result
    
    def setPrecursorErrorType(self, int mono_isotopic ):
        """
        setPrecursorErrorType(self, mono_isotopic: int ) -> None
        """
        assert mono_isotopic in [0, 1], 'arg mono_isotopic wrong type'
    
        self.inst.get().setPrecursorErrorType((<_MassType>mono_isotopic))
    
    def getPrecursorErrorType(self):
        """
        getPrecursorErrorType(self) -> int
        """
        cdef _MassType _r = self.inst.get().getPrecursorErrorType()
        py_result = <int>_r
        return py_result
    
    def setFragmentMassErrorUnit(self, int unit ):
        """
        setFragmentMassErrorUnit(self, unit: int ) -> None
        """
        assert unit in [0, 1], 'arg unit wrong type'
    
        self.inst.get().setFragmentMassErrorUnit((<_ErrorUnit>unit))
    
    def getFragmentMassErrorUnit(self):
        """
        getFragmentMassErrorUnit(self) -> int
        """
        cdef _ErrorUnit _r = self.inst.get().getFragmentMassErrorUnit()
        py_result = <int>_r
        return py_result
    
    def setPrecursorMassErrorUnit(self, int unit ):
        """
        setPrecursorMassErrorUnit(self, unit: int ) -> None
        """
        assert unit in [0, 1], 'arg unit wrong type'
    
        self.inst.get().setPrecursorMassErrorUnit((<_ErrorUnit>unit))
    
    def getPrecursorMassErrorUnit(self):
        """
        getPrecursorMassErrorUnit(self) -> int
        """
        cdef _ErrorUnit _r = self.inst.get().getPrecursorMassErrorUnit()
        py_result = <int>_r
        return py_result
    
    def setNumberOfThreads(self,  threads ):
        """
        setNumberOfThreads(self, threads: int ) -> None
        """
        assert isinstance(threads, int), 'arg threads wrong type'
    
        self.inst.get().setNumberOfThreads((<unsigned int>threads))
    
    def getNumberOfThreads(self):
        """
        getNumberOfThreads(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getNumberOfThreads()
        py_result = <unsigned int>_r
        return py_result
    
    def setModifications(self, ModificationDefinitionsSet mods ):
        """
        setModifications(self, mods: ModificationDefinitionsSet ) -> None
        """
        assert isinstance(mods, ModificationDefinitionsSet), 'arg mods wrong type'
    
        self.inst.get().setModifications((deref(mods.inst.get())))
    
    def getModifications(self):
        """
        getModifications(self) -> ModificationDefinitionsSet
        """
        cdef _ModificationDefinitionsSet * _r = new _ModificationDefinitionsSet(self.inst.get().getModifications())
        cdef ModificationDefinitionsSet py_result = ModificationDefinitionsSet.__new__(ModificationDefinitionsSet)
        py_result.inst = shared_ptr[_ModificationDefinitionsSet](_r)
        return py_result
    
    def setOutputFilename(self,  output ):
        """
        setOutputFilename(self, output: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(output, str) or isinstance(output, bytes) or isinstance(output, String)), 'arg output wrong type'
    
        self.inst.get().setOutputFilename(deref((convString(output)).get()))
    
    def getOutputFilename(self):
        """
        getOutputFilename(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOutputFilename()
        py_result = convOutputString(_r)
        return py_result
    
    def setInputFilename(self,  input_file ):
        """
        setInputFilename(self, input_file: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(input_file, str) or isinstance(input_file, bytes) or isinstance(input_file, String)), 'arg input_file wrong type'
    
        self.inst.get().setInputFilename(deref((convString(input_file)).get()))
    
    def getInputFilename(self):
        """
        getInputFilename(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getInputFilename()
        py_result = convOutputString(_r)
        return py_result
    
    def setTaxonomyFilename(self,  filename ):
        """
        setTaxonomyFilename(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().setTaxonomyFilename(deref((convString(filename)).get()))
    
    def getTaxonomyFilename(self):
        """
        getTaxonomyFilename(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTaxonomyFilename()
        py_result = convOutputString(_r)
        return py_result
    
    def setDefaultParametersFilename(self,  filename ):
        """
        setDefaultParametersFilename(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().setDefaultParametersFilename(deref((convString(filename)).get()))
    
    def getDefaultParametersFilename(self):
        """
        getDefaultParametersFilename(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getDefaultParametersFilename()
        py_result = convOutputString(_r)
        return py_result
    
    def setTaxon(self,  taxon ):
        """
        setTaxon(self, taxon: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(taxon, str) or isinstance(taxon, bytes) or isinstance(taxon, String)), 'arg taxon wrong type'
    
        self.inst.get().setTaxon(deref((convString(taxon)).get()))
    
    def getTaxon(self):
        """
        getTaxon(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTaxon()
        py_result = convOutputString(_r)
        return py_result
    
    def setMaxPrecursorCharge(self,  max_charge ):
        """
        setMaxPrecursorCharge(self, max_charge: int ) -> None
        """
        assert isinstance(max_charge, int), 'arg max_charge wrong type'
    
        self.inst.get().setMaxPrecursorCharge((<int>max_charge))
    
    def getMaxPrecursorCharge(self):
        """
        getMaxPrecursorCharge(self) -> int
        """
        cdef int _r = self.inst.get().getMaxPrecursorCharge()
        py_result = <int>_r
        return py_result
    
    def setNumberOfMissedCleavages(self,  missed_cleavages ):
        """
        setNumberOfMissedCleavages(self, missed_cleavages: int ) -> None
        """
        assert isinstance(missed_cleavages, int), 'arg missed_cleavages wrong type'
    
        self.inst.get().setNumberOfMissedCleavages((<unsigned int>missed_cleavages))
    
    def getNumberOfMissedCleavages(self):
        """
        getNumberOfMissedCleavages(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getNumberOfMissedCleavages()
        py_result = <unsigned int>_r
        return py_result
    
    def setOutputResults(self,  result ):
        """
        setOutputResults(self, result: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(result, str) or isinstance(result, bytes) or isinstance(result, String)), 'arg result wrong type'
    
        self.inst.get().setOutputResults(deref((convString(result)).get()))
    
    def getOutputResults(self):
        """
        getOutputResults(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOutputResults()
        py_result = convOutputString(_r)
        return py_result
    
    def setMaxValidEValue(self, double value ):
        """
        setMaxValidEValue(self, value: float ) -> None
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setMaxValidEValue((<double>value))
    
    def getMaxValidEValue(self):
        """
        getMaxValidEValue(self) -> float
        """
        cdef double _r = self.inst.get().getMaxValidEValue()
        py_result = <double>_r
        return py_result
    
    def setSemiCleavage(self, bool semi_cleavage ):
        """
        setSemiCleavage(self, semi_cleavage: bool ) -> None
        """
        assert isinstance(semi_cleavage, pybool_t), 'arg semi_cleavage wrong type'
    
        self.inst.get().setSemiCleavage((<bool>semi_cleavage))
    
    def setAllowIsotopeError(self, bool allow_isotope_error ):
        """
        setAllowIsotopeError(self, allow_isotope_error: bool ) -> None
        """
        assert isinstance(allow_isotope_error, pybool_t), 'arg allow_isotope_error wrong type'
    
        self.inst.get().setAllowIsotopeError((<bool>allow_isotope_error))
    
    def write(self,  filename , bool ignore_member_parameters , bool force_default_mods ):
        """
        write(self, filename: Union[bytes, str, String] , ignore_member_parameters: bool , force_default_mods: bool ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(ignore_member_parameters, pybool_t), 'arg ignore_member_parameters wrong type'
        assert isinstance(force_default_mods, pybool_t), 'arg force_default_mods wrong type'
    
    
    
        self.inst.get().write(deref((convString(filename)).get()), (<bool>ignore_member_parameters), (<bool>force_default_mods))
    
    def setCleavageSite(self,  cleavage_site ):
        """
        setCleavageSite(self, cleavage_site: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(cleavage_site, str) or isinstance(cleavage_site, bytes) or isinstance(cleavage_site, String)), 'arg cleavage_site wrong type'
    
        self.inst.get().setCleavageSite(deref((convString(cleavage_site)).get()))
    
    def getCleavageSite(self):
        """
        getCleavageSite(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getCleavageSite()
        py_result = convOutputString(_r)
        return py_result
    ErrorUnit = __ErrorUnit
    MassType = __MassType 
from libc.stddef cimport ptrdiff_t
from libc.stdint cimport *
from libcpp.map cimport map as libcpp_map
cimport numpy as np
import numpy as np
ctypedef libcpp_vector[ double ] _DoubleList
ctypedef libcpp_vector[ int ] _IntList


# The format field of the Python buffer interface is a string that specifies the type of the elements in the buffer. It uses the same syntax as the struct module in Python. Here are some possible values:
# 'b': signed byte
# 'B': unsigned byte
# 'h': signed short
# 'H': unsigned short
# 'i': signed int
# 'I': unsigned int
# 'l': signed long
# 'L': unsigned long
# 'q': signed long long
# 'Q': unsigned long long
# 'f': float
# 'd': double
# '?': bool
# 'c': char
# 's': char[]
# 'p': char[] (Pascal string)
# 'P': void *

# define ArrayWrapperFloat as holding in a vector
cdef class ArrayWrapperFloat:
    cdef libcpp_vector[float] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]

    # constructor and destructor are fairly unimportant now since
    # vec will be destroyed automatically.

    cdef set_data(self, libcpp_vector[float]& data):
       self.vec.swap(data)

    # now implement the buffer protocol for the class
    # which makes it generally useful to anything that expects an array
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # relevant documentation http://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
        cdef Py_ssize_t itemsize = sizeof(self.vec[0])

        self.shape[0] = self.vec.size()
        self.strides[0] = itemsize
        buffer.buf = <char *>&(self.vec[0])
        buffer.format = 'f'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.vec.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

# define ArrayWrapperDouble as holding in a vector
cdef class ArrayWrapperDouble:
    cdef libcpp_vector[double] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]

    # constructor and destructor are fairly unimportant now since
    # vec will be destroyed automatically.

    cdef set_data(self, libcpp_vector[double]& data):
       self.vec.swap(data)

    # now implement the buffer protocol for the class
    # which makes it generally useful to anything that expects an array
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # relevant documentation http://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
        cdef Py_ssize_t itemsize = sizeof(self.vec[0])

        self.shape[0] = self.vec.size()
        self.strides[0] = itemsize
        buffer.buf = <char *>&(self.vec[0])
        buffer.format = 'd'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.vec.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL 
 
 
 

## Add this to all modules
class Interfaces:

    BinaryDataArray = _Interfaces_BinaryDataArray
    Spectrum = _Interfaces_Spectrum
    Chromatogram = _Interfaces_Chromatogram 

# Add this to all modules
PeakSpectrum = MSSpectrum

PeakMap = MSExperiment 
