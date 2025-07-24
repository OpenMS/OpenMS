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
def __static_FeatureMapping_assignMS2IndexToFeature(MSExperiment spectra , FeatureMapping_FeatureMappingInfo fm_info , double precursor_mz_tolerance , double precursor_rt_tolerance , bool ppm ):
    """
    __static_FeatureMapping_assignMS2IndexToFeature(spectra: MSExperiment , fm_info: FeatureMapping_FeatureMappingInfo , precursor_mz_tolerance: float , precursor_rt_tolerance: float , ppm: bool ) -> FeatureMapping_FeatureToMs2Indices
    """
    assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    assert isinstance(fm_info, FeatureMapping_FeatureMappingInfo), 'arg fm_info wrong type'
    assert isinstance(precursor_mz_tolerance, float), 'arg precursor_mz_tolerance wrong type'
    assert isinstance(precursor_rt_tolerance, float), 'arg precursor_rt_tolerance wrong type'
    assert isinstance(ppm, pybool_t), 'arg ppm wrong type'





    cdef _FeatureMapping_FeatureToMs2Indices * _r = new _FeatureMapping_FeatureToMs2Indices(_assignMS2IndexToFeature_FeatureMapping((deref(spectra.inst.get())), (deref(fm_info.inst.get())), (<double>precursor_mz_tolerance), (<double>precursor_rt_tolerance), (<bool>ppm)))
    cdef FeatureMapping_FeatureToMs2Indices py_result = FeatureMapping_FeatureToMs2Indices.__new__(FeatureMapping_FeatureToMs2Indices)
    py_result.inst = shared_ptr[_FeatureMapping_FeatureToMs2Indices](_r)
    return py_result 

cdef class __InletType:
    None
    INLETNULL = 0
    DIRECT = 1
    BATCH = 2
    CHROMATOGRAPHY = 3
    PARTICLEBEAM = 4
    MEMBRANESEPARATOR = 5
    OPENSPLIT = 6
    JETSEPARATOR = 7
    SEPTUM = 8
    RESERVOIR = 9
    MOVINGBELT = 10
    MOVINGWIRE = 11
    FLOWINJECTIONANALYSIS = 12
    ELECTROSPRAYINLET = 13
    THERMOSPRAYINLET = 14
    INFUSION = 15
    CONTINUOUSFLOWFASTATOMBOMBARDMENT = 16
    INDUCTIVELYCOUPLEDPLASMA = 17
    MEMBRANE = 18
    NANOSPRAY = 19
    SIZE_OF_INLETTYPE = 20

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class IonOpticsType:
    None
    UNKNOWN = 0
    MAGNETIC_DEFLECTION = 1
    DELAYED_EXTRACTION = 2
    COLLISION_QUADRUPOLE = 3
    SELECTED_ION_FLOW_TUBE = 4
    TIME_LAG_FOCUSING = 5
    REFLECTRON = 6
    EINZEL_LENS = 7
    FIRST_STABILITY_REGION = 8
    FRINGING_FIELD = 9
    KINETIC_ENERGY_ANALYZER = 10
    STATIC_FIELD = 11
    SIZE_OF_IONOPTICSTYPE = 12

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __IonizationMethod:
    None
    IONMETHODNULL = 0
    ESI = 1
    EI = 2
    CI = 3
    FAB = 4
    TSP = 5
    LD = 6
    FD = 7
    FI = 8
    PD = 9
    SI = 10
    TI = 11
    API = 12
    ISI = 13
    CID = 14
    CAD = 15
    HN = 16
    APCI = 17
    APPI = 18
    ICP = 19
    NESI = 20
    MESI = 21
    SELDI = 22
    SEND = 23
    FIB = 24
    MALDI = 25
    MPI = 26
    DI = 27
    FA = 28
    FII = 29
    GD_MS = 30
    NICI = 31
    NRMS = 32
    PI = 33
    PYMS = 34
    REMPI = 35
    AI = 36
    ASI = 37
    AD = 38
    AUI = 39
    CEI = 40
    CHEMI = 41
    DISSI = 42
    LSI = 43
    PEI = 44
    SOI = 45
    SPI = 46
    SUI = 47
    VI = 48
    AP_MALDI = 49
    SILI = 50
    SALDI = 51
    SIZE_OF_IONIZATIONMETHOD = 52

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class Measure:
    None
    MEASURE_PPM = 0
    MEASURE_DA = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __OpenPepXLLFAlgorithm_ExitCodes:
    None
    EXECUTION_OK = 0
    ILLEGAL_PARAMETERS = 1
    UNEXPECTED_RESULT = 2
    INCOMPATIBLE_INPUT_DATA = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __Polarity:
    None
    POLNULL = 0
    POSITIVE = 1
    NEGATIVE = 2
    SIZE_OF_POLARITY = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class Adduct:
    """
    Cython implementation of _Adduct

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Adduct.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Adduct rv = Adduct.__new__(Adduct)
       rv.inst = shared_ptr[_Adduct](new _Adduct(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Adduct rv = Adduct.__new__(Adduct)
       rv.inst = shared_ptr[_Adduct](new _Adduct(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Adduct](new _Adduct())
    
    def _init_1(self, Adduct in_0 ):
        """
        _init_1(self, in_0: Adduct ) -> None
        """
        assert isinstance(in_0, Adduct), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Adduct](new _Adduct((deref(in_0.inst.get()))))
    
    def _init_2(self,  charge ):
        """
        _init_2(self, charge: int ) -> None
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst = shared_ptr[_Adduct](new _Adduct((<int>charge)))
    
    def _init_3(self,  charge ,  amount , double singleMass ,  formula , double log_prob , double rt_shift ,  label ):
        """
        _init_3(self, charge: int , amount: int , singleMass: float , formula: Union[bytes, str, String] , log_prob: float , rt_shift: float , label: Union[bytes, str, String] ) -> None
        """
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(amount, int), 'arg amount wrong type'
        assert isinstance(singleMass, float), 'arg singleMass wrong type'
        assert (isinstance(formula, str) or isinstance(formula, bytes) or isinstance(formula, String)), 'arg formula wrong type'
        assert isinstance(log_prob, float), 'arg log_prob wrong type'
        assert isinstance(rt_shift, float), 'arg rt_shift wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
    
    
    
    
    
    
        self.inst = shared_ptr[_Adduct](new _Adduct((<int>charge), (<int>amount), (<double>singleMass), deref((convString(formula)).get()), (<double>log_prob), (<double>rt_shift), deref((convString(label)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Adduct ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, charge: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, charge: int , amount: int , singleMass: float , formula: Union[bytes, str, String] , log_prob: float , rt_shift: float , label: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Adduct)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
             self._init_2(*args)
        elif (len(args)==7) and (isinstance(args[0], int)) and (isinstance(args[1], int)) and (isinstance(args[2], float)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and (isinstance(args[4], float)) and (isinstance(args[5], float)) and ((isinstance(args[6], str) or isinstance(args[6], bytes) or isinstance(args[6], String))):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def setCharge(self,  charge ):
        """
        setCharge(self, charge: int ) -> None
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setCharge((<int>charge))
    
    def getAmount(self):
        """
        getAmount(self) -> int
        """
        cdef int _r = self.inst.get().getAmount()
        py_result = <int>_r
        return py_result
    
    def setAmount(self,  amount ):
        """
        setAmount(self, amount: int ) -> None
        """
        assert isinstance(amount, int), 'arg amount wrong type'
    
        self.inst.get().setAmount((<int>amount))
    
    def getSingleMass(self):
        """
        getSingleMass(self) -> float
        """
        cdef double _r = self.inst.get().getSingleMass()
        py_result = <double>_r
        return py_result
    
    def setSingleMass(self, double singleMass ):
        """
        setSingleMass(self, singleMass: float ) -> None
        """
        assert isinstance(singleMass, float), 'arg singleMass wrong type'
    
        self.inst.get().setSingleMass((<double>singleMass))
    
    def getLogProb(self):
        """
        getLogProb(self) -> float
        """
        cdef double _r = self.inst.get().getLogProb()
        py_result = <double>_r
        return py_result
    
    def setLogProb(self, double log_prob ):
        """
        setLogProb(self, log_prob: float ) -> None
        """
        assert isinstance(log_prob, float), 'arg log_prob wrong type'
    
        self.inst.get().setLogProb((<double>log_prob))
    
    def getFormula(self):
        """
        getFormula(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getFormula()
        py_result = convOutputString(_r)
        return py_result
    
    def setFormula(self,  formula ):
        """
        setFormula(self, formula: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(formula, str) or isinstance(formula, bytes) or isinstance(formula, String)), 'arg formula wrong type'
    
        self.inst.get().setFormula(deref((convString(formula)).get()))
    
    def getRTShift(self):
        """
        getRTShift(self) -> float
        """
        cdef double _r = self.inst.get().getRTShift()
        py_result = <double>_r
        return py_result
    
    def getLabel(self):
        """
        getLabel(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getLabel()
        py_result = convOutputString(_r)
        return py_result 

cdef class AverageLinkage:
    """
    Cython implementation of _AverageLinkage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AverageLinkage.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AverageLinkage rv = AverageLinkage.__new__(AverageLinkage)
       rv.inst = shared_ptr[_AverageLinkage](new _AverageLinkage(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AverageLinkage rv = AverageLinkage.__new__(AverageLinkage)
       rv.inst = shared_ptr[_AverageLinkage](new _AverageLinkage(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AverageLinkage](new _AverageLinkage())
    
    def _init_1(self, AverageLinkage in_0 ):
        """
        _init_1(self, in_0: AverageLinkage ) -> None
        """
        assert isinstance(in_0, AverageLinkage), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AverageLinkage](new _AverageLinkage((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AverageLinkage ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AverageLinkage)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ChromatogramExtractorAlgorithm:
    """
    Cython implementation of _ChromatogramExtractorAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramExtractorAlgorithm.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramExtractorAlgorithm rv = ChromatogramExtractorAlgorithm.__new__(ChromatogramExtractorAlgorithm)
       rv.inst = shared_ptr[_ChromatogramExtractorAlgorithm](new _ChromatogramExtractorAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramExtractorAlgorithm rv = ChromatogramExtractorAlgorithm.__new__(ChromatogramExtractorAlgorithm)
       rv.inst = shared_ptr[_ChromatogramExtractorAlgorithm](new _ChromatogramExtractorAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChromatogramExtractorAlgorithm](new _ChromatogramExtractorAlgorithm())
    
    def _init_1(self, ChromatogramExtractorAlgorithm in_0 ):
        """
        _init_1(self, in_0: ChromatogramExtractorAlgorithm ) -> None
        """
        assert isinstance(in_0, ChromatogramExtractorAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramExtractorAlgorithm](new _ChromatogramExtractorAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramExtractorAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramExtractorAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractChromatograms(self, SpectrumAccessOpenMS input , list output , list extraction_coordinates , double mz_extraction_window , bool ppm , double im_extraction_window ,  filter ):
        """
        extractChromatograms(self, input: SpectrumAccessOpenMS , output: List[OSChromatogram] , extraction_coordinates: List[ExtractionCoordinates] , mz_extraction_window: float , ppm: bool , im_extraction_window: float , filter: Union[bytes, str, String] ) -> None
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

cdef class ConsensusFeature:
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

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ConsensusFeature rv = ConsensusFeature.__new__(ConsensusFeature)
       rv.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ConsensusFeature rv = ConsensusFeature.__new__(ConsensusFeature)
       rv.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature())
    
    def _init_1(self, ConsensusFeature in_0 ):
        """
        _init_1(self, in_0: ConsensusFeature ) -> None
        """
        assert isinstance(in_0, ConsensusFeature), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature((deref(in_0.inst.get()))))
    
    def _init_2(self,  in_0 , Peak2D in_1 ,  in_2 ):
        """
        _init_2(self, in_0: int , in_1: Peak2D , in_2: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
        assert isinstance(in_1, Peak2D), 'arg in_1 wrong type'
        assert isinstance(in_2, int) and in_2 >= 0, 'arg in_2 wrong type'
    
    
    
        self.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature((<uint64_t>in_0), (deref(in_1.inst.get())), (<uint64_t>in_2)))
    
    def _init_3(self,  in_0 , BaseFeature in_1 ):
        """
        _init_3(self, in_0: int , in_1: BaseFeature ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
        assert isinstance(in_1, BaseFeature), 'arg in_1 wrong type'
    
    
        self.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature((<uint64_t>in_0), (deref(in_1.inst.get()))))
    
    def _init_4(self,  in_0 , ConsensusFeature in_1 ):
        """
        _init_4(self, in_0: int , in_1: ConsensusFeature ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
        assert isinstance(in_1, ConsensusFeature), 'arg in_1 wrong type'
    
    
        self.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature((<uint64_t>in_0), (deref(in_1.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ConsensusFeature ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: int , in_1: Peak2D , in_2: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: int , in_1: BaseFeature ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: int , in_1: ConsensusFeature ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusFeature)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], Peak2D)) and (isinstance(args[2], int) and args[2] >= 0):
             self._init_2(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], BaseFeature)):
             self._init_3(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], ConsensusFeature)):
             self._init_4(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def computeConsensus(self):
        """
        computeConsensus(self) -> None
        Computes and updates the consensus position, intensity, and charge
        """
        self.inst.get().computeConsensus()
    
    def computeMonoisotopicConsensus(self):
        """
        computeMonoisotopicConsensus(self) -> None
        Computes and updates the consensus position, intensity, and charge
        """
        self.inst.get().computeMonoisotopicConsensus()
    
    def computeDechargeConsensus(self, FeatureMap in_0 , bool in_1 ):
        """
        computeDechargeConsensus(self, in_0: FeatureMap , in_1: bool ) -> None
        Computes the uncharged parent RT & mass, assuming the handles are charge variants
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, pybool_t), 'arg in_1 wrong type'
    
    
        self.inst.get().computeDechargeConsensus((deref(in_0.inst.get())), (<bool>in_1))
    
    def _insert_0(self,  map_idx , Peak2D in_1 ,  element_idx ):
        """
        _insert_0(self, map_idx: int , in_1: Peak2D , element_idx: int ) -> None
        """
        assert isinstance(map_idx, int) and map_idx >= 0, 'arg map_idx wrong type'
        assert isinstance(in_1, Peak2D), 'arg in_1 wrong type'
        assert isinstance(element_idx, int) and element_idx >= 0, 'arg element_idx wrong type'
    
    
    
        self.inst.get().insert((<uint64_t>map_idx), (deref(in_1.inst.get())), (<uint64_t>element_idx))
    
    def _insert_1(self,  map_idx , BaseFeature in_1 ):
        """
        _insert_1(self, map_idx: int , in_1: BaseFeature ) -> None
        """
        assert isinstance(map_idx, int) and map_idx >= 0, 'arg map_idx wrong type'
        assert isinstance(in_1, BaseFeature), 'arg in_1 wrong type'
    
    
        self.inst.get().insert((<uint64_t>map_idx), (deref(in_1.inst.get())))
    
    def _insert_2(self,  map_idx , ConsensusFeature in_1 ):
        """
        _insert_2(self, map_idx: int , in_1: ConsensusFeature ) -> None
        """
        assert isinstance(map_idx, int) and map_idx >= 0, 'arg map_idx wrong type'
        assert isinstance(in_1, ConsensusFeature), 'arg in_1 wrong type'
    
    
        self.inst.get().insert((<uint64_t>map_idx), (deref(in_1.inst.get())))
    
    def insert(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: insert(self, map_idx: int , in_1: Peak2D , element_idx: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: insert(self, map_idx: int , in_1: BaseFeature ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: insert(self, map_idx: int , in_1: ConsensusFeature ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], Peak2D)) and (isinstance(args[2], int) and args[2] >= 0):
            return self._insert_0(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], BaseFeature)):
            return self._insert_1(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], ConsensusFeature)):
            return self._insert_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getFeatureList(self):
        """
        getFeatureList(self) -> List[FeatureHandle]
        """
        _r = self.inst.get().getFeatureList()
        py_result = []
        cdef libcpp_vector[_FeatureHandle].iterator it__r = _r.begin()
        cdef FeatureHandle item_py_result
        while it__r != _r.end():
           item_py_result = FeatureHandle.__new__(FeatureHandle)
           item_py_result.inst = shared_ptr[_FeatureHandle](new _FeatureHandle(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def addRatio(self, Ratio r ):
        """
        addRatio(self, r: Ratio ) -> None
        Connects a ratio to the ConsensusFeature.
        """
        assert isinstance(r, Ratio), 'arg r wrong type'
    
        self.inst.get().addRatio((deref(r.inst.get())))
    
    def setRatios(self, list rs ):
        """
        setRatios(self, rs: List[Ratio] ) -> None
        Connects the ratios to the ConsensusFeature.
        """
        assert isinstance(rs, list) and all(isinstance(elemt_rec, Ratio) for elemt_rec in rs), 'arg rs wrong type'
        cdef libcpp_vector[_Ratio] * v0 = new libcpp_vector[_Ratio]()
        cdef Ratio item0
        for item0 in rs:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setRatios(deref(v0))
        del v0
    
    def getRatios(self):
        """
        getRatios(self) -> List[Ratio]
        Get the ratio vector.
        """
        _r = self.inst.get().getRatios()
        py_result = []
        cdef libcpp_vector[_Ratio].iterator it__r = _r.begin()
        cdef Ratio item_py_result
        while it__r != _r.end():
           item_py_result = Ratio.__new__(Ratio)
           item_py_result.inst = shared_ptr[_Ratio](new _Ratio(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def getUniqueId(self):
        """
        getUniqueId(self) -> int
        Returns the unique id
        """
        cdef size_t _r = self.inst.get().getUniqueId()
        py_result = <size_t>_r
        return py_result
    
    def clearUniqueId(self):
        """
        clearUniqueId(self) -> int
        Clear the unique id. The new unique id will be invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        cdef size_t _r = self.inst.get().clearUniqueId()
        py_result = <size_t>_r
        return py_result
    
    def hasValidUniqueId(self):
        """
        hasValidUniqueId(self) -> int
        Returns whether the unique id is valid. Returns 1 if the unique id is valid, 0 otherwise
        """
        cdef size_t _r = self.inst.get().hasValidUniqueId()
        py_result = <size_t>_r
        return py_result
    
    def hasInvalidUniqueId(self):
        """
        hasInvalidUniqueId(self) -> int
        Returns whether the unique id is invalid. Returns 1 if the unique id is invalid, 0 otherwise
        """
        cdef size_t _r = self.inst.get().hasInvalidUniqueId()
        py_result = <size_t>_r
        return py_result
    
    def setUniqueId(self,  rhs ):
        """
        setUniqueId(self, rhs: int ) -> None
        Assigns a new, valid unique id. Always returns 1
        """
        assert isinstance(rhs, int) and rhs >= 0, 'arg rhs wrong type'
    
        self.inst.get().setUniqueId((<uint64_t>rhs))
    
    def ensureUniqueId(self):
        """
        ensureUniqueId(self) -> int
        Assigns a valid unique id, but only if the present one is invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        cdef size_t _r = self.inst.get().ensureUniqueId()
        py_result = <size_t>_r
        return py_result
    
    def isValid(self,  unique_id ):
        """
        isValid(self, unique_id: int ) -> bool
        Returns true if the unique_id is valid, false otherwise
        """
        assert isinstance(unique_id, int) and unique_id >= 0, 'arg unique_id wrong type'
    
        cdef bool _r = self.inst.get().isValid((<uint64_t>unique_id))
        py_result = <bool>_r
        return py_result
    
    def getQuality(self):
        """
        getQuality(self) -> float
        Returns the overall quality
        """
        cdef float _r = self.inst.get().getQuality()
        py_result = <float>_r
        return py_result
    
    def setQuality(self, float q ):
        """
        setQuality(self, q: float ) -> None
        Sets the overall quality
        """
        assert isinstance(q, float), 'arg q wrong type'
    
        self.inst.get().setQuality((<float>q))
    
    def getWidth(self):
        """
        getWidth(self) -> float
        Returns the features width (full width at half max, FWHM)
        """
        cdef float _r = self.inst.get().getWidth()
        py_result = <float>_r
        return py_result
    
    def setWidth(self, float q ):
        """
        setWidth(self, q: float ) -> None
        Sets the width of the feature (FWHM)
        """
        assert isinstance(q, float), 'arg q wrong type'
    
        self.inst.get().setWidth((<float>q))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns the charge state
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def setCharge(self,  q ):
        """
        setCharge(self, q: int ) -> None
        Sets the charge state
        """
        assert isinstance(q, int), 'arg q wrong type'
    
        self.inst.get().setCharge((<int>q))
    
    def getAnnotationState(self):
        """
        getAnnotationState(self) -> int
        State of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only
        """
        cdef _AnnotationState _r = self.inst.get().getAnnotationState()
        py_result = <int>_r
        return py_result
    
    def getPeptideIdentifications(self):
        """
        getPeptideIdentifications(self) -> PeptideIdentificationList
        Returns the PeptideIdentification vector
        """
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().getPeptideIdentifications())
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result
    
    def setPeptideIdentifications(self, PeptideIdentificationList peptides ):
        """
        setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None
        Sets the PeptideIdentification vector
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().setPeptideIdentifications((deref(peptides.inst.get())))
    
    def getIntensity(self):
        """
        getIntensity(self) -> float
        Returns the data point intensity (height)
        """
        cdef float _r = self.inst.get().getIntensity()
        py_result = <float>_r
        return py_result
    
    def getMZ(self):
        """
        getMZ(self) -> float
        Returns the m/z coordinate (index 1)
        """
        cdef double _r = self.inst.get().getMZ()
        py_result = <double>_r
        return py_result
    
    def getRT(self):
        """
        getRT(self) -> float
        Returns the RT coordinate (index 0)
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def setMZ(self, double in_0 ):
        """
        setMZ(self, in_0: float ) -> None
        Returns the m/z coordinate (index 1)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setMZ((<double>in_0))
    
    def setRT(self, double in_0 ):
        """
        setRT(self, in_0: float ) -> None
        Returns the RT coordinate (index 0)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setRT((<double>in_0))
    
    def setIntensity(self, float in_0 ):
        """
        setIntensity(self, in_0: float ) -> None
        Returns the data point intensity (height)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setIntensity((<float>in_0))
    
    def isMetaEmpty(self):
        """
        isMetaEmpty(self) -> bool
        Returns if the MetaInfo is empty
        """
        cdef bool _r = self.inst.get().isMetaEmpty()
        py_result = <bool>_r
        return py_result
    
    def clearMetaInfo(self):
        """
        clearMetaInfo(self) -> None
        Removes all meta values
        """
        self.inst.get().clearMetaInfo()
    
    def metaRegistry(self):
        """
        metaRegistry(self) -> MetaInfoRegistry
        Returns a reference to the MetaInfoRegistry
        """
        cdef _MetaInfoRegistry * _r = new _MetaInfoRegistry(self.inst.get().metaRegistry())
        cdef MetaInfoRegistry py_result = MetaInfoRegistry.__new__(MetaInfoRegistry)
        py_result.inst = shared_ptr[_MetaInfoRegistry](_r)
        return py_result
    
    def getKeys(self, list keys ):
        """
        getKeys(self, keys: List[bytes] ) -> None
        Fills the given vector with a list of all keys for which a value is set
        """
        assert isinstance(keys, list) and all(isinstance(i, bytes) for i in keys), 'arg keys wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in keys:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getKeys(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        keys[:] = replace
        del v0
    
    def getMetaValue(self,  in_0 ):
        """
        getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to a string, or
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef _DataValue _r = self.inst.get().getMetaValue(deref((convString(in_0)).get()))
        cdef DataValue _value = DataValue.__new__(DataValue)
        _value.inst = shared_ptr[_DataValue](new _DataValue(_r))
        cdef int _type = _r.valueType()
        cdef object py_result
        if _type == DataType.STRING_VALUE:
            py_result = _value.toString()
        elif _type == DataType.INT_VALUE:
            py_result = _value.toInt()
        elif _type == DataType.DOUBLE_VALUE:
            py_result = _value.toDouble()
        elif _type == DataType.INT_LIST:
            py_result = _value.toIntList()
        elif _type == DataType.DOUBLE_LIST:
            py_result = _value.toDoubleList()
        elif _type == DataType.STRING_LIST:
            py_result = _value.toStringList()
        elif _type == DataType.EMPTY_VALUE:
            py_result = None
        else:
            raise Exception("DataValue instance has invalid value type %d" % _type)
        return py_result
    
    def setMetaValue(self,  in_0 ,  in_1 ):
        """
        setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the DataValue corresponding to a name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue(deref((convString(in_0)).get()), deref(DataValue(in_1).inst.get()))
    
    def metaValueExists(self,  in_0 ):
        """
        metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool
        Returns whether an entry with the given name exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists(deref((convString(in_0)).get()))
        py_result = <bool>_r
        return py_result
    
    def removeMetaValue(self,  in_0 ):
        """
        removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None
        Removes the DataValue corresponding to `name` if it exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue(deref((convString(in_0)).get()))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ConsensusFeature):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ConsensusFeature other_casted = other
        cdef ConsensusFeature self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ConsensusXMLFile:
    """
    Cython implementation of _ConsensusXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusXMLFile](new _ConsensusXMLFile())
    
    def load(self,  in_0 , ConsensusMap in_1 ):
        """
        load(self, in_0: Union[bytes, str, String] , in_1: ConsensusMap ) -> None
        Loads a consensus map from file and calls updateRanges
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, ConsensusMap), 'arg in_1 wrong type'
    
    
        self.inst.get().load(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def store(self,  in_0 , ConsensusMap in_1 ):
        """
        store(self, in_0: Union[bytes, str, String] , in_1: ConsensusMap ) -> None
        Stores a consensus map to file
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, ConsensusMap), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        Mutable access to the options for loading/storing
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
        return py_result 

cdef class ExtractionCoordinates:
    """
    Cython implementation of _ExtractionCoordinates

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExtractionCoordinates.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property mz:
        def __set__(self, double mz):
        
            self.inst.get().mz = (<double>mz)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz
            py_result = <double>_r
            return py_result
    
    property mz_precursor:
        def __set__(self, double mz_precursor):
        
            self.inst.get().mz_precursor = (<double>mz_precursor)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz_precursor
            py_result = <double>_r
            return py_result
    
    property rt_start:
        def __set__(self, double rt_start):
        
            self.inst.get().rt_start = (<double>rt_start)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_start
            py_result = <double>_r
            return py_result
    
    property rt_end:
        def __set__(self, double rt_end):
        
            self.inst.get().rt_end = (<double>rt_end)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_end
            py_result = <double>_r
            return py_result
    
    property ion_mobility:
        def __set__(self, double ion_mobility):
        
            self.inst.get().ion_mobility = (<double>ion_mobility)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ion_mobility
            py_result = <double>_r
            return py_result
    
    property id:
        def __set__(self, bytes id):
        
            self.inst.get().id = (<libcpp_string>id)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().id
            py_result = <libcpp_string>_r
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ExtractionCoordinates](new _ExtractionCoordinates())
    
    def _init_1(self, ExtractionCoordinates in_0 ):
        """
        _init_1(self, in_0: ExtractionCoordinates ) -> None
        """
        assert isinstance(in_0, ExtractionCoordinates), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ExtractionCoordinates](new _ExtractionCoordinates((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ExtractionCoordinates ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ExtractionCoordinates)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class FeatureMapping:
    """
    Cython implementation of _FeatureMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureMapping rv = FeatureMapping.__new__(FeatureMapping)
       rv.inst = shared_ptr[_FeatureMapping](new _FeatureMapping(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureMapping rv = FeatureMapping.__new__(FeatureMapping)
       rv.inst = shared_ptr[_FeatureMapping](new _FeatureMapping(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureMapping](new _FeatureMapping())
    
    def _init_1(self, FeatureMapping in_0 ):
        """
        _init_1(self, in_0: FeatureMapping ) -> None
        """
        assert isinstance(in_0, FeatureMapping), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureMapping](new _FeatureMapping((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureMapping ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureMapping)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    assignMS2IndexToFeature = __static_FeatureMapping_assignMS2IndexToFeature 

cdef class FeatureMapping_FeatureMappingInfo:
    """
    Cython implementation of _FeatureMapping_FeatureMappingInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping_FeatureMappingInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureMapping_FeatureMappingInfo rv = FeatureMapping_FeatureMappingInfo.__new__(FeatureMapping_FeatureMappingInfo)
       rv.inst = shared_ptr[_FeatureMapping_FeatureMappingInfo](new _FeatureMapping_FeatureMappingInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureMapping_FeatureMappingInfo rv = FeatureMapping_FeatureMappingInfo.__new__(FeatureMapping_FeatureMappingInfo)
       rv.inst = shared_ptr[_FeatureMapping_FeatureMappingInfo](new _FeatureMapping_FeatureMappingInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureMapping_FeatureMappingInfo](new _FeatureMapping_FeatureMappingInfo())
    
    def _init_1(self, FeatureMapping_FeatureMappingInfo in_0 ):
        """
        _init_1(self, in_0: FeatureMapping_FeatureMappingInfo ) -> None
        """
        assert isinstance(in_0, FeatureMapping_FeatureMappingInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureMapping_FeatureMappingInfo](new _FeatureMapping_FeatureMappingInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureMapping_FeatureMappingInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureMapping_FeatureMappingInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class FeatureMapping_FeatureToMs2Indices:
    """
    Cython implementation of _FeatureMapping_FeatureToMs2Indices

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMapping_FeatureToMs2Indices.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureMapping_FeatureToMs2Indices rv = FeatureMapping_FeatureToMs2Indices.__new__(FeatureMapping_FeatureToMs2Indices)
       rv.inst = shared_ptr[_FeatureMapping_FeatureToMs2Indices](new _FeatureMapping_FeatureToMs2Indices(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureMapping_FeatureToMs2Indices rv = FeatureMapping_FeatureToMs2Indices.__new__(FeatureMapping_FeatureToMs2Indices)
       rv.inst = shared_ptr[_FeatureMapping_FeatureToMs2Indices](new _FeatureMapping_FeatureToMs2Indices(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureMapping_FeatureToMs2Indices](new _FeatureMapping_FeatureToMs2Indices())
    
    def _init_1(self, FeatureMapping_FeatureToMs2Indices in_0 ):
        """
        _init_1(self, in_0: FeatureMapping_FeatureToMs2Indices ) -> None
        """
        assert isinstance(in_0, FeatureMapping_FeatureToMs2Indices), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureMapping_FeatureToMs2Indices](new _FeatureMapping_FeatureToMs2Indices((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureMapping_FeatureToMs2Indices ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureMapping_FeatureToMs2Indices)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class GNPSQuantificationFile:
    """
    Cython implementation of _GNPSQuantificationFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSQuantificationFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GNPSQuantificationFile rv = GNPSQuantificationFile.__new__(GNPSQuantificationFile)
       rv.inst = shared_ptr[_GNPSQuantificationFile](new _GNPSQuantificationFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GNPSQuantificationFile rv = GNPSQuantificationFile.__new__(GNPSQuantificationFile)
       rv.inst = shared_ptr[_GNPSQuantificationFile](new _GNPSQuantificationFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_GNPSQuantificationFile](new _GNPSQuantificationFile())
    
    def _init_1(self, GNPSQuantificationFile in_0 ):
        """
        _init_1(self, in_0: GNPSQuantificationFile ) -> None
        """
        assert isinstance(in_0, GNPSQuantificationFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GNPSQuantificationFile](new _GNPSQuantificationFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GNPSQuantificationFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], GNPSQuantificationFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self, ConsensusMap consensus_map ,  output_file ):
        """
        store(self, consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None
        Write feature quantification table (txt file) from a ConsensusMap. Required for GNPS FBMN.
        
        The table contains map information on the featureXML files from which the ConsensusMap was generated as well as
        a row for every consensus feature with information on rt, mz, intensity, width and quality. The same information is
        added for each original feature in the consensus feature.
        
        :param consensus_map: Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        :param output_file: Output file path for the feature quantification table.
        """
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
        assert (isinstance(output_file, str) or isinstance(output_file, bytes) or isinstance(output_file, String)), 'arg output_file wrong type'
    
    
        self.inst.get().store((deref(consensus_map.inst.get())), deref((convString(output_file)).get())) 

cdef class IDMapper:
    """
    Cython implementation of _IDMapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDMapper.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IDMapper rv = IDMapper.__new__(IDMapper)
       rv.inst = shared_ptr[_IDMapper](new _IDMapper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IDMapper rv = IDMapper.__new__(IDMapper)
       rv.inst = shared_ptr[_IDMapper](new _IDMapper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications
        """
        self.inst = shared_ptr[_IDMapper](new _IDMapper())
    
    def _init_1(self, IDMapper in_0 ):
        """
        _init_1(self, in_0: IDMapper ) -> None
        """
        assert isinstance(in_0, IDMapper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IDMapper](new _IDMapper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IDMapper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IDMapper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _annotate_0(self, AnnotatedMSRun map_ , PeptideIdentificationList ids , list protein_ids , bool clear_ids , bool mapMS1 ):
        """
        _annotate_0(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None
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
        assert isinstance(map_, AnnotatedMSRun), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(clear_ids, pybool_t), 'arg clear_ids wrong type'
        assert isinstance(mapMS1, pybool_t), 'arg mapMS1 wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>clear_ids), (<bool>mapMS1))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def _annotate_1(self, AnnotatedMSRun map_ , FeatureMap fmap , bool clear_ids , bool mapMS1 ):
        """
        _annotate_1(self, map_: AnnotatedMSRun , fmap: FeatureMap , clear_ids: bool , mapMS1: bool ) -> None
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
        assert isinstance(map_, AnnotatedMSRun), 'arg map_ wrong type'
        assert isinstance(fmap, FeatureMap), 'arg fmap wrong type'
        assert isinstance(clear_ids, pybool_t), 'arg clear_ids wrong type'
        assert isinstance(mapMS1, pybool_t), 'arg mapMS1 wrong type'
    
    
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(fmap.inst.get())), (<bool>clear_ids), (<bool>mapMS1))
    
    def _annotate_2(self, FeatureMap map_ , PeptideIdentificationList ids , list protein_ids , bool use_centroid_rt , bool use_centroid_mz , MSExperiment spectra ):
        """
        _annotate_2(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None
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
        assert isinstance(map_, FeatureMap), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(use_centroid_rt, pybool_t), 'arg use_centroid_rt wrong type'
        assert isinstance(use_centroid_mz, pybool_t), 'arg use_centroid_mz wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>use_centroid_rt), (<bool>use_centroid_mz), (deref(spectra.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def _annotate_3(self, ConsensusMap map_ , PeptideIdentificationList ids , list protein_ids , bool measure_from_subelements , bool annotate_ids_with_subelements , MSExperiment spectra ):
        """
        _annotate_3(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None
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
        assert isinstance(map_, ConsensusMap), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(measure_from_subelements, pybool_t), 'arg measure_from_subelements wrong type'
        assert isinstance(annotate_ids_with_subelements, pybool_t), 'arg annotate_ids_with_subelements wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>measure_from_subelements), (<bool>annotate_ids_with_subelements), (deref(spectra.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def _annotate_4(self, AnnotatedMSRun map_ , PeptideIdentificationList ids , list protein_ids , bool clear_ids , bool mapMS1 ):
        """
        _annotate_4(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None
        Mapping method using PeptideIdentificationList\n
        
        The identifications stored in a PeptideIdentificationList instance can be added to the
        corresponding spectrum
        
        
        :param map: AnnotatedMSRun to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param clear_ids: Reset peptide and protein identifications of each scan before annotating
        :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
        """
        assert isinstance(map_, AnnotatedMSRun), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(clear_ids, pybool_t), 'arg clear_ids wrong type'
        assert isinstance(mapMS1, pybool_t), 'arg mapMS1 wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>clear_ids), (<bool>mapMS1))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def _annotate_5(self, FeatureMap map_ , PeptideIdentificationList ids , list protein_ids , bool use_centroid_rt , bool use_centroid_mz , MSExperiment spectra ):
        """
        _annotate_5(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None
        Mapping method using PeptideIdentificationList\n
        
        :param map: FeatureMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param use_centroid_rt: Whether to use the RT value of feature centroids even if convex hulls are present
        :param use_centroid_mz: Whether to use the m/z value of feature centroids even if convex hulls are present
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        """
        assert isinstance(map_, FeatureMap), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(use_centroid_rt, pybool_t), 'arg use_centroid_rt wrong type'
        assert isinstance(use_centroid_mz, pybool_t), 'arg use_centroid_mz wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>use_centroid_rt), (<bool>use_centroid_mz), (deref(spectra.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def _annotate_6(self, ConsensusMap map_ , PeptideIdentificationList ids , list protein_ids , bool measure_from_subelements , bool annotate_ids_with_subelements , MSExperiment spectra ):
        """
        _annotate_6(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None
        Mapping method using PeptideIdentificationList\n
        
        :param map: ConsensusMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param measure_from_subelements: Boolean operator set to true if distance estimate from FeatureHandles instead of Centroid
        :param annotate_ids_with_subelements: Boolean operator set to true if store map index of FeatureHandle in peptide identification
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        """
        assert isinstance(map_, ConsensusMap), 'arg map_ wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(measure_from_subelements, pybool_t), 'arg measure_from_subelements wrong type'
        assert isinstance(annotate_ids_with_subelements, pybool_t), 'arg annotate_ids_with_subelements wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
    
    
        self.inst.get().annotate((deref(map_.inst.get())), (deref(ids.inst.get())), deref(v2), (<bool>measure_from_subelements), (<bool>annotate_ids_with_subelements), (deref(spectra.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
    
    def annotate(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None
          :noindex:
        
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
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: AnnotatedMSRun , fmap: FeatureMap , clear_ids: bool , mapMS1: bool ) -> None
          :noindex:
        
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
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None
          :noindex:
        
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
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None
          :noindex:
        
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
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: AnnotatedMSRun , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , clear_ids: bool , mapMS1: bool ) -> None
          :noindex:
        
        Mapping method using PeptideIdentificationList\n
        
        The identifications stored in a PeptideIdentificationList instance can be added to the
        corresponding spectrum
        
        
        :param map: AnnotatedMSRun to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param clear_ids: Reset peptide and protein identifications of each scan before annotating
        :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: FeatureMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , use_centroid_rt: bool , use_centroid_mz: bool , spectra: MSExperiment ) -> None
          :noindex:
        
        Mapping method using PeptideIdentificationList\n
        
        :param map: FeatureMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param use_centroid_rt: Whether to use the RT value of feature centroids even if convex hulls are present
        :param use_centroid_mz: Whether to use the m/z value of feature centroids even if convex hulls are present
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
        
        .. rubric:: Overload:
        .. py:function:: annotate(self, map_: ConsensusMap , ids: PeptideIdentificationList , protein_ids: List[ProteinIdentification] , measure_from_subelements: bool , annotate_ids_with_subelements: bool , spectra: MSExperiment ) -> None
          :noindex:
        
        Mapping method using PeptideIdentificationList\n
        
        :param map: ConsensusMap to receive the identifications
        :param ids: PeptideIdentificationList for the MSExperiment
        :param protein_ids: ProteinIdentification for the MSExperiment
        :param measure_from_subelements: Boolean operator set to true if distance estimate from FeatureHandles instead of Centroid
        :param annotate_ids_with_subelements: Boolean operator set to true if store map index of FeatureHandle in peptide identification
        :param spectra: Whether precursors not contained in the identifications are annotated with an empty PeptideIdentification object containing the scan index
    
        """
        if (len(args)==5) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)):
            return self._annotate_0(*args)
        elif (len(args)==4) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], FeatureMap)) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], pybool_t)):
            return self._annotate_1(*args)
        elif (len(args)==6) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], MSExperiment)):
            return self._annotate_2(*args)
        elif (len(args)==6) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], MSExperiment)):
            return self._annotate_3(*args)
        elif (len(args)==5) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)):
            return self._annotate_4(*args)
        elif (len(args)==6) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], MSExperiment)):
            return self._annotate_5(*args)
        elif (len(args)==6) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[2])) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], MSExperiment)):
            return self._annotate_6(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _mapPrecursorsToIdentifications_0(self, MSExperiment spectra , PeptideIdentificationList ids , double mz_tol , double rt_tol ):
        """
        _mapPrecursorsToIdentifications_0(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState
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
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(mz_tol, float), 'arg mz_tol wrong type'
        assert isinstance(rt_tol, float), 'arg rt_tol wrong type'
    
    
    
    
        cdef _IDMapper_PeptideIdentificationListState * _r = new _IDMapper_PeptideIdentificationListState(self.inst.get().mapPrecursorsToIdentifications((deref(spectra.inst.get())), (deref(ids.inst.get())), (<double>mz_tol), (<double>rt_tol)))
        cdef IDMapper_PeptideIdentificationListState py_result = IDMapper_PeptideIdentificationListState.__new__(IDMapper_PeptideIdentificationListState)
        py_result.inst = shared_ptr[_IDMapper_PeptideIdentificationListState](_r)
        return py_result
    
    def _mapPrecursorsToIdentifications_1(self, MSExperiment spectra , PeptideIdentificationList ids , double mz_tol , double rt_tol ):
        """
        _mapPrecursorsToIdentifications_1(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState
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
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(mz_tol, float), 'arg mz_tol wrong type'
        assert isinstance(rt_tol, float), 'arg rt_tol wrong type'
    
    
    
    
        cdef _IDMapper_PeptideIdentificationListState * _r = new _IDMapper_PeptideIdentificationListState(self.inst.get().mapPrecursorsToIdentifications((deref(spectra.inst.get())), (deref(ids.inst.get())), (<double>mz_tol), (<double>rt_tol)))
        cdef IDMapper_PeptideIdentificationListState py_result = IDMapper_PeptideIdentificationListState.__new__(IDMapper_PeptideIdentificationListState)
        py_result.inst = shared_ptr[_IDMapper_PeptideIdentificationListState](_r)
        return py_result
    
    def mapPrecursorsToIdentifications(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: mapPrecursorsToIdentifications(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState
          :noindex:
        
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
        
        .. rubric:: Overload:
        .. py:function:: mapPrecursorsToIdentifications(self, spectra: MSExperiment , ids: PeptideIdentificationList , mz_tol: float , rt_tol: float ) -> IDMapper_PeptideIdentificationListState
          :noindex:
        
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
        if (len(args)==4) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._mapPrecursorsToIdentifications_0(*args)
        elif (len(args)==4) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._mapPrecursorsToIdentifications_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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

cdef class IDMapper_PeptideIdentificationListState:
    """
    Cython implementation of _IDMapper_PeptideIdentificationListState

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDMapper_PeptideIdentificationListState.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property no_precursors:
        def __set__(self, list no_precursors):
            cdef libcpp_vector[size_t] v0 = no_precursors
            self.inst.get().no_precursors = v0
            
    
        def __get__(self):
            _r = self.inst.get().no_precursors
            cdef list py_result = _r
            return py_result
    
    property identified:
        def __set__(self, list identified):
            cdef libcpp_vector[size_t] v0 = identified
            self.inst.get().identified = v0
            
    
        def __get__(self):
            _r = self.inst.get().identified
            cdef list py_result = _r
            return py_result
    
    property unidentified:
        def __set__(self, list unidentified):
            cdef libcpp_vector[size_t] v0 = unidentified
            self.inst.get().unidentified = v0
            
    
        def __get__(self):
            _r = self.inst.get().unidentified
            cdef list py_result = _r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_IDMapper_PeptideIdentificationListState](new _IDMapper_PeptideIdentificationListState()) 

cdef class IMSElement:
    """
    Cython implementation of _IMSElement

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSElement_1_1IMSElement.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IMSElement rv = IMSElement.__new__(IMSElement)
       rv.inst = shared_ptr[_IMSElement](new _IMSElement(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IMSElement rv = IMSElement.__new__(IMSElement)
       rv.inst = shared_ptr[_IMSElement](new _IMSElement(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Represents a chemical atom with name and isotope distribution
        """
        self.inst = shared_ptr[_IMSElement](new _IMSElement())
    
    def _init_1(self, IMSElement in_0 ):
        """
        _init_1(self, in_0: IMSElement ) -> None
        """
        assert isinstance(in_0, IMSElement), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IMSElement](new _IMSElement((deref(in_0.inst.get()))))
    
    def _init_2(self, bytes name , IMSIsotopeDistribution isotopes ):
        """
        _init_2(self, name: bytes , isotopes: IMSIsotopeDistribution ) -> None
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(isotopes, IMSIsotopeDistribution), 'arg isotopes wrong type'
    
    
        self.inst = shared_ptr[_IMSElement](new _IMSElement((<libcpp_string>name), (deref(isotopes.inst.get()))))
    
    def _init_3(self, bytes name , double mass ):
        """
        _init_3(self, name: bytes , mass: float ) -> None
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(mass, float), 'arg mass wrong type'
    
    
        self.inst = shared_ptr[_IMSElement](new _IMSElement((<libcpp_string>name), (<double>mass)))
    
    def _init_4(self, bytes name ,  nominal_mass ):
        """
        _init_4(self, name: bytes , nominal_mass: int ) -> None
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(nominal_mass, int), 'arg nominal_mass wrong type'
    
    
        self.inst = shared_ptr[_IMSElement](new _IMSElement((<libcpp_string>name), (<unsigned int>nominal_mass)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Represents a chemical atom with name and isotope distribution

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IMSElement ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: bytes , isotopes: IMSIsotopeDistribution ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: bytes , mass: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: bytes , nominal_mass: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSElement)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], bytes)) and (isinstance(args[1], IMSIsotopeDistribution)):
             self._init_2(*args)
        elif (len(args)==2) and (isinstance(args[0], bytes)) and (isinstance(args[1], float)):
             self._init_3(*args)
        elif (len(args)==2) and (isinstance(args[0], bytes)) and (isinstance(args[1], int)):
             self._init_4(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getName(self):
        """
        getName(self) -> bytes
        Gets element's name
        """
        cdef libcpp_string _r = self.inst.get().getName()
        py_result = <libcpp_string>_r
        return py_result
    
    def setName(self, bytes name ):
        """
        setName(self, name: bytes ) -> None
        Sets element's name
        """
        assert isinstance(name, bytes), 'arg name wrong type'
    
        self.inst.get().setName((<libcpp_string>name))
    
    def getSequence(self):
        """
        getSequence(self) -> bytes
        Gets element's sequence
        """
        cdef libcpp_string _r = self.inst.get().getSequence()
        py_result = <libcpp_string>_r
        return py_result
    
    def setSequence(self, bytes sequence ):
        """
        setSequence(self, sequence: bytes ) -> None
        Sets element's sequence
        """
        assert isinstance(sequence, bytes), 'arg sequence wrong type'
    
        self.inst.get().setSequence((<libcpp_string>sequence))
    
    def getNominalMass(self):
        """
        getNominalMass(self) -> int
        Gets element's nominal mass
        """
        cdef unsigned int _r = self.inst.get().getNominalMass()
        py_result = <unsigned int>_r
        return py_result
    
    def getMass(self,  index ):
        """
        getMass(self, index: int ) -> float
        Gets mass of element's isotope 'index'
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getMass((<int>index))
        py_result = <double>_r
        return py_result
    
    def getAverageMass(self):
        """
        getAverageMass(self) -> float
        Gets element's average mass
        """
        cdef double _r = self.inst.get().getAverageMass()
        py_result = <double>_r
        return py_result
    
    def getIonMass(self,  electrons_number ):
        """
        getIonMass(self, electrons_number: int ) -> float
        Gets ion mass of element. By default ion lacks 1 electron, but this can be changed by setting other 'electrons_number'
        """
        assert isinstance(electrons_number, int), 'arg electrons_number wrong type'
    
        cdef double _r = self.inst.get().getIonMass((<int>electrons_number))
        py_result = <double>_r
        return py_result
    
    def getIsotopeDistribution(self):
        """
        getIsotopeDistribution(self) -> IMSIsotopeDistribution
        Gets element's isotope distribution
        """
        cdef _IMSIsotopeDistribution * _r = new _IMSIsotopeDistribution(self.inst.get().getIsotopeDistribution())
        cdef IMSIsotopeDistribution py_result = IMSIsotopeDistribution.__new__(IMSIsotopeDistribution)
        py_result.inst = shared_ptr[_IMSIsotopeDistribution](_r)
        return py_result
    
    def setIsotopeDistribution(self, IMSIsotopeDistribution isotopes ):
        """
        setIsotopeDistribution(self, isotopes: IMSIsotopeDistribution ) -> None
        Sets element's isotope distribution
        """
        assert isinstance(isotopes, IMSIsotopeDistribution), 'arg isotopes wrong type'
    
        self.inst.get().setIsotopeDistribution((deref(isotopes.inst.get())))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, IMSElement):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IMSElement other_casted = other
        cdef IMSElement self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Instrument:
    """
    Cython implementation of _Instrument

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Instrument.html>`_
      -- Inherits from ['MetaInfoInterface']

    Description of a MS instrument
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Instrument rv = Instrument.__new__(Instrument)
       rv.inst = shared_ptr[_Instrument](new _Instrument(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Instrument rv = Instrument.__new__(Instrument)
       rv.inst = shared_ptr[_Instrument](new _Instrument(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Description of a MS instrument
        """
        self.inst = shared_ptr[_Instrument](new _Instrument())
    
    def _init_1(self, Instrument in_0 ):
        """
        _init_1(self, in_0: Instrument ) -> None
        """
        assert isinstance(in_0, Instrument), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Instrument](new _Instrument((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Description of a MS instrument

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Instrument ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Instrument)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the instrument
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the instrument
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getVendor(self):
        """
        getVendor(self) -> Union[bytes, str, String]
        Returns the instrument vendor
        """
        cdef _String _r = self.inst.get().getVendor()
        py_result = convOutputString(_r)
        return py_result
    
    def setVendor(self,  vendor ):
        """
        setVendor(self, vendor: Union[bytes, str, String] ) -> None
        Sets the instrument vendor
        """
        assert (isinstance(vendor, str) or isinstance(vendor, bytes) or isinstance(vendor, String)), 'arg vendor wrong type'
    
        self.inst.get().setVendor(deref((convString(vendor)).get()))
    
    def getModel(self):
        """
        getModel(self) -> Union[bytes, str, String]
        Returns the instrument model
        """
        cdef _String _r = self.inst.get().getModel()
        py_result = convOutputString(_r)
        return py_result
    
    def setModel(self,  model ):
        """
        setModel(self, model: Union[bytes, str, String] ) -> None
        Sets the instrument model
        """
        assert (isinstance(model, str) or isinstance(model, bytes) or isinstance(model, String)), 'arg model wrong type'
    
        self.inst.get().setModel(deref((convString(model)).get()))
    
    def getCustomizations(self):
        """
        getCustomizations(self) -> Union[bytes, str, String]
        Returns a description of customizations
        """
        cdef _String _r = self.inst.get().getCustomizations()
        py_result = convOutputString(_r)
        return py_result
    
    def setCustomizations(self,  customizations ):
        """
        setCustomizations(self, customizations: Union[bytes, str, String] ) -> None
        Sets the a description of customizations
        """
        assert (isinstance(customizations, str) or isinstance(customizations, bytes) or isinstance(customizations, String)), 'arg customizations wrong type'
    
        self.inst.get().setCustomizations(deref((convString(customizations)).get()))
    
    def getIonSources(self):
        """
        getIonSources(self) -> List[IonSource]
        Returns the ion source list
        """
        _r = self.inst.get().getIonSources()
        py_result = []
        cdef libcpp_vector[_IonSource].iterator it__r = _r.begin()
        cdef IonSource item_py_result
        while it__r != _r.end():
           item_py_result = IonSource.__new__(IonSource)
           item_py_result.inst = shared_ptr[_IonSource](new _IonSource(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setIonSources(self, list ion_sources ):
        """
        setIonSources(self, ion_sources: List[IonSource] ) -> None
        Sets the ion source list
        """
        assert isinstance(ion_sources, list) and all(isinstance(elemt_rec, IonSource) for elemt_rec in ion_sources), 'arg ion_sources wrong type'
        cdef libcpp_vector[_IonSource] * v0 = new libcpp_vector[_IonSource]()
        cdef IonSource item0
        for item0 in ion_sources:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIonSources(deref(v0))
        del v0
    
    def getMassAnalyzers(self):
        """
        getMassAnalyzers(self) -> List[MassAnalyzer]
        Returns the mass analyzer list
        """
        _r = self.inst.get().getMassAnalyzers()
        py_result = []
        cdef libcpp_vector[_MassAnalyzer].iterator it__r = _r.begin()
        cdef MassAnalyzer item_py_result
        while it__r != _r.end():
           item_py_result = MassAnalyzer.__new__(MassAnalyzer)
           item_py_result.inst = shared_ptr[_MassAnalyzer](new _MassAnalyzer(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setMassAnalyzers(self, list mass_analyzers ):
        """
        setMassAnalyzers(self, mass_analyzers: List[MassAnalyzer] ) -> None
        Sets the mass analyzer list
        """
        assert isinstance(mass_analyzers, list) and all(isinstance(elemt_rec, MassAnalyzer) for elemt_rec in mass_analyzers), 'arg mass_analyzers wrong type'
        cdef libcpp_vector[_MassAnalyzer] * v0 = new libcpp_vector[_MassAnalyzer]()
        cdef MassAnalyzer item0
        for item0 in mass_analyzers:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setMassAnalyzers(deref(v0))
        del v0
    
    def getIonDetectors(self):
        """
        getIonDetectors(self) -> List[IonDetector]
        Returns the ion detector list
        """
        _r = self.inst.get().getIonDetectors()
        py_result = []
        cdef libcpp_vector[_IonDetector].iterator it__r = _r.begin()
        cdef IonDetector item_py_result
        while it__r != _r.end():
           item_py_result = IonDetector.__new__(IonDetector)
           item_py_result.inst = shared_ptr[_IonDetector](new _IonDetector(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setIonDetectors(self, list ion_detectors ):
        """
        setIonDetectors(self, ion_detectors: List[IonDetector] ) -> None
        Sets the ion detector list
        """
        assert isinstance(ion_detectors, list) and all(isinstance(elemt_rec, IonDetector) for elemt_rec in ion_detectors), 'arg ion_detectors wrong type'
        cdef libcpp_vector[_IonDetector] * v0 = new libcpp_vector[_IonDetector]()
        cdef IonDetector item0
        for item0 in ion_detectors:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIonDetectors(deref(v0))
        del v0
    
    def getSoftware(self):
        """
        getSoftware(self) -> Software
        Returns the instrument software
        """
        cdef _Software * _r = new _Software(self.inst.get().getSoftware())
        cdef Software py_result = Software.__new__(Software)
        py_result.inst = shared_ptr[_Software](_r)
        return py_result
    
    def setSoftware(self, Software software ):
        """
        setSoftware(self, software: Software ) -> None
        Sets the instrument software
        """
        assert isinstance(software, Software), 'arg software wrong type'
    
        self.inst.get().setSoftware((deref(software.inst.get())))
    
    def getIonOptics(self):
        """
        getIonOptics(self) -> int
        Returns the ion optics type
        """
        cdef _IonOpticsType _r = self.inst.get().getIonOptics()
        py_result = <int>_r
        return py_result
    
    def setIonOptics(self, int ion_optics ):
        """
        setIonOptics(self, ion_optics: int ) -> None
        Sets the ion optics type
        """
        assert ion_optics in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'arg ion_optics wrong type'
    
        self.inst.get().setIonOptics((<_IonOpticsType>ion_optics))
    
    def isMetaEmpty(self):
        """
        isMetaEmpty(self) -> bool
        Returns if the MetaInfo is empty
        """
        cdef bool _r = self.inst.get().isMetaEmpty()
        py_result = <bool>_r
        return py_result
    
    def clearMetaInfo(self):
        """
        clearMetaInfo(self) -> None
        Removes all meta values
        """
        self.inst.get().clearMetaInfo()
    
    def metaRegistry(self):
        """
        metaRegistry(self) -> MetaInfoRegistry
        Returns a reference to the MetaInfoRegistry
        """
        cdef _MetaInfoRegistry * _r = new _MetaInfoRegistry(self.inst.get().metaRegistry())
        cdef MetaInfoRegistry py_result = MetaInfoRegistry.__new__(MetaInfoRegistry)
        py_result.inst = shared_ptr[_MetaInfoRegistry](_r)
        return py_result
    
    def getKeys(self, list keys ):
        """
        getKeys(self, keys: List[bytes] ) -> None
        Fills the given vector with a list of all keys for which a value is set
        """
        assert isinstance(keys, list) and all(isinstance(i, bytes) for i in keys), 'arg keys wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in keys:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getKeys(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        keys[:] = replace
        del v0
    
    def getMetaValue(self,  in_0 ):
        """
        getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to a string, or
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef _DataValue _r = self.inst.get().getMetaValue(deref((convString(in_0)).get()))
        cdef DataValue _value = DataValue.__new__(DataValue)
        _value.inst = shared_ptr[_DataValue](new _DataValue(_r))
        cdef int _type = _r.valueType()
        cdef object py_result
        if _type == DataType.STRING_VALUE:
            py_result = _value.toString()
        elif _type == DataType.INT_VALUE:
            py_result = _value.toInt()
        elif _type == DataType.DOUBLE_VALUE:
            py_result = _value.toDouble()
        elif _type == DataType.INT_LIST:
            py_result = _value.toIntList()
        elif _type == DataType.DOUBLE_LIST:
            py_result = _value.toDoubleList()
        elif _type == DataType.STRING_LIST:
            py_result = _value.toStringList()
        elif _type == DataType.EMPTY_VALUE:
            py_result = None
        else:
            raise Exception("DataValue instance has invalid value type %d" % _type)
        return py_result
    
    def setMetaValue(self,  in_0 ,  in_1 ):
        """
        setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the DataValue corresponding to a name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue(deref((convString(in_0)).get()), deref(DataValue(in_1).inst.get()))
    
    def metaValueExists(self,  in_0 ):
        """
        metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool
        Returns whether an entry with the given name exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists(deref((convString(in_0)).get()))
        py_result = <bool>_r
        return py_result
    
    def removeMetaValue(self,  in_0 ):
        """
        removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None
        Removes the DataValue corresponding to `name` if it exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue(deref((convString(in_0)).get()))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Instrument):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Instrument other_casted = other
        cdef Instrument self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class IonSource:
    """
    Cython implementation of _IonSource

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonSource.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IonSource rv = IonSource.__new__(IonSource)
       rv.inst = shared_ptr[_IonSource](new _IonSource(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IonSource rv = IonSource.__new__(IonSource)
       rv.inst = shared_ptr[_IonSource](new _IonSource(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Description of an ion source (part of a MS Instrument)
        """
        self.inst = shared_ptr[_IonSource](new _IonSource())
    
    def _init_1(self, IonSource in_0 ):
        """
        _init_1(self, in_0: IonSource ) -> None
        """
        assert isinstance(in_0, IonSource), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IonSource](new _IonSource((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Description of an ion source (part of a MS Instrument)

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IonSource ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IonSource)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPolarity(self):
        """
        getPolarity(self) -> int
        Returns the ionization mode
        """
        cdef _Polarity _r = self.inst.get().getPolarity()
        py_result = <int>_r
        return py_result
    
    def setPolarity(self, int polarity ):
        """
        setPolarity(self, polarity: int ) -> None
        Sets the ionization mode
        """
        assert polarity in [0, 1, 2, 3], 'arg polarity wrong type'
    
        self.inst.get().setPolarity((<_Polarity>polarity))
    
    def getInletType(self):
        """
        getInletType(self) -> int
        Returns the inlet type
        """
        cdef _InletType _r = self.inst.get().getInletType()
        py_result = <int>_r
        return py_result
    
    def setInletType(self, int inlet_type ):
        """
        setInletType(self, inlet_type: int ) -> None
        Sets the inlet type
        """
        assert inlet_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'arg inlet_type wrong type'
    
        self.inst.get().setInletType((<_InletType>inlet_type))
    
    def getIonizationMethod(self):
        """
        getIonizationMethod(self) -> int
        Returns the ionization method
        """
        cdef _IonizationMethod _r = self.inst.get().getIonizationMethod()
        py_result = <int>_r
        return py_result
    
    def setIonizationMethod(self, int ionization_type ):
        """
        setIonizationMethod(self, ionization_type: int ) -> None
        Sets the ionization method
        """
        assert ionization_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52], 'arg ionization_type wrong type'
    
        self.inst.get().setIonizationMethod((<_IonizationMethod>ionization_type))
    
    def getOrder(self):
        """
        getOrder(self) -> int
        Returns the position of this part in the whole Instrument
        
        Order can be ignored, as long the instrument has this default setup:
          - one ion source
          - one or many mass analyzers
          - one ion detector
        
        For more complex instruments, the order should be defined.
        """
        cdef int _r = self.inst.get().getOrder()
        py_result = <int>_r
        return py_result
    
    def setOrder(self,  order ):
        """
        setOrder(self, order: int ) -> None
        Sets the order
        """
        assert isinstance(order, int), 'arg order wrong type'
    
        self.inst.get().setOrder((<int>order))
    
    def isMetaEmpty(self):
        """
        isMetaEmpty(self) -> bool
        Returns if the MetaInfo is empty
        """
        cdef bool _r = self.inst.get().isMetaEmpty()
        py_result = <bool>_r
        return py_result
    
    def clearMetaInfo(self):
        """
        clearMetaInfo(self) -> None
        Removes all meta values
        """
        self.inst.get().clearMetaInfo()
    
    def metaRegistry(self):
        """
        metaRegistry(self) -> MetaInfoRegistry
        Returns a reference to the MetaInfoRegistry
        """
        cdef _MetaInfoRegistry * _r = new _MetaInfoRegistry(self.inst.get().metaRegistry())
        cdef MetaInfoRegistry py_result = MetaInfoRegistry.__new__(MetaInfoRegistry)
        py_result.inst = shared_ptr[_MetaInfoRegistry](_r)
        return py_result
    
    def getKeys(self, list keys ):
        """
        getKeys(self, keys: List[bytes] ) -> None
        Fills the given vector with a list of all keys for which a value is set
        """
        assert isinstance(keys, list) and all(isinstance(i, bytes) for i in keys), 'arg keys wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in keys:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getKeys(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        keys[:] = replace
        del v0
    
    def getMetaValue(self,  in_0 ):
        """
        getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to a string, or
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef _DataValue _r = self.inst.get().getMetaValue(deref((convString(in_0)).get()))
        cdef DataValue _value = DataValue.__new__(DataValue)
        _value.inst = shared_ptr[_DataValue](new _DataValue(_r))
        cdef int _type = _r.valueType()
        cdef object py_result
        if _type == DataType.STRING_VALUE:
            py_result = _value.toString()
        elif _type == DataType.INT_VALUE:
            py_result = _value.toInt()
        elif _type == DataType.DOUBLE_VALUE:
            py_result = _value.toDouble()
        elif _type == DataType.INT_LIST:
            py_result = _value.toIntList()
        elif _type == DataType.DOUBLE_LIST:
            py_result = _value.toDoubleList()
        elif _type == DataType.STRING_LIST:
            py_result = _value.toStringList()
        elif _type == DataType.EMPTY_VALUE:
            py_result = None
        else:
            raise Exception("DataValue instance has invalid value type %d" % _type)
        return py_result
    
    def setMetaValue(self,  in_0 ,  in_1 ):
        """
        setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the DataValue corresponding to a name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue(deref((convString(in_0)).get()), deref(DataValue(in_1).inst.get()))
    
    def metaValueExists(self,  in_0 ):
        """
        metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool
        Returns whether an entry with the given name exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists(deref((convString(in_0)).get()))
        py_result = <bool>_r
        return py_result
    
    def removeMetaValue(self,  in_0 ):
        """
        removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None
        Removes the DataValue corresponding to `name` if it exists
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue(deref((convString(in_0)).get()))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, IonSource):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IonSource other_casted = other
        cdef IonSource self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    InletType = __InletType
    IonizationMethod = __IonizationMethod
    Polarity = __Polarity 

cdef class OpenPepXLLFAlgorithm:
    """
    Cython implementation of _OpenPepXLLFAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenPepXLLFAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']

    Search for cross-linked peptide pairs in tandem MS spectra
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenPepXLLFAlgorithm rv = OpenPepXLLFAlgorithm.__new__(OpenPepXLLFAlgorithm)
       rv.inst = shared_ptr[_OpenPepXLLFAlgorithm](new _OpenPepXLLFAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenPepXLLFAlgorithm rv = OpenPepXLLFAlgorithm.__new__(OpenPepXLLFAlgorithm)
       rv.inst = shared_ptr[_OpenPepXLLFAlgorithm](new _OpenPepXLLFAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenPepXLLFAlgorithm](new _OpenPepXLLFAlgorithm())
    
    def _init_1(self, OpenPepXLLFAlgorithm in_0 ):
        """
        _init_1(self, in_0: OpenPepXLLFAlgorithm ) -> None
        """
        assert isinstance(in_0, OpenPepXLLFAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenPepXLLFAlgorithm](new _OpenPepXLLFAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenPepXLLFAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenPepXLLFAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment unprocessed_spectra , list fasta_db , list protein_ids , PeptideIdentificationList peptide_ids , list all_top_csms , MSExperiment spectra ):
        """
        run(self, unprocessed_spectra: MSExperiment , fasta_db: List[FASTAEntry] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment ) -> int
        Performs the main function of this class, the search for cross-linked peptides
        
        
        :param unprocessed_spectra: The input PeakMap of experimental spectra
        :param fasta_db: The protein database containing targets and decoys
        :param protein_ids: A result vector containing search settings. Should contain one PeptideIdentification
        :param peptide_ids: A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty
        :param all_top_csms: A result vector containing cross-link spectrum matches as CrossLinkSpectrumMatches. Should be empty. This is only necessary for writing out xQuest type spectrum files
        :param spectra: A result vector containing the input spectra after preprocessing and filtering. Should be empty. This is only necessary for writing out xQuest type spectrum files
        """
        assert isinstance(unprocessed_spectra, MSExperiment), 'arg unprocessed_spectra wrong type'
        assert isinstance(fasta_db, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in fasta_db), 'arg fasta_db wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(all_top_csms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in all_top_csms), 'arg all_top_csms wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
        cdef libcpp_vector[_FASTAEntry] * v1 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item1
        for item1 in fasta_db:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in protein_ids:
            v2.push_back(deref(item2.inst.get()))
    
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]] * v4 = new libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]]()
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v4_rec = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item4_rec
        cdef libcpp_vector[_CrossLinkSpectrumMatch].iterator it_all_top_csms_rec
        for all_top_csms_rec in all_top_csms:
            v4_rec.clear()
            for item4_rec in all_top_csms_rec:
                v4_rec.push_back(deref(item4_rec.inst.get()))
            v4.push_back(deref(v4_rec))
    
        cdef _OpenPepXLLFAlgorithm_ExitCodes _r = self.inst.get().run((deref(unprocessed_spectra.inst.get())), deref(v1), deref(v2), (deref(peptide_ids.inst.get())), deref(v4), (deref(spectra.inst.get())))
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]].iterator it_all_top_csms = v4.begin()
        replace_0 = []
        while it_all_top_csms != v4.end():
            it_all_top_csms_rec = deref(it_all_top_csms).begin()
            replace_1 = []
            while it_all_top_csms_rec != deref(it_all_top_csms).end():
                item4_rec = CrossLinkSpectrumMatch.__new__(CrossLinkSpectrumMatch)
                item4_rec.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch(deref(it_all_top_csms_rec)))
                replace_1.append(item4_rec)
                inc(it_all_top_csms_rec)
            replace_0.append(replace_1)
            inc(it_all_top_csms)
        all_top_csms[:] = replace_0
        del v4
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v2.begin()
        replace_0 = []
        while it_protein_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item2)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v2
        cdef libcpp_vector[_FASTAEntry].iterator it_fasta_db = v1.begin()
        replace_0 = []
        while it_fasta_db != v1.end():
            item1 = FASTAEntry.__new__(FASTAEntry)
            item1.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_fasta_db)))
            replace_0.append(item1)
            inc(it_fasta_db)
        fasta_db[:] = replace_0
        del v1
        py_result = <int>_r
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
    OpenPepXLLFAlgorithm_ExitCodes = __OpenPepXLLFAlgorithm_ExitCodes 

cdef class Param:
    """
    Cython implementation of _Param

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Param.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Param rv = Param.__new__(Param)
       rv.inst = shared_ptr[_Param](new _Param(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Param rv = Param.__new__(Param)
       rv.inst = shared_ptr[_Param](new _Param(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Param](new _Param())
    
    def _init_1(self, Param in_0 ):
        """
        _init_1(self, in_0: Param ) -> None
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Param](new _Param((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Param ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Param)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setValue_0(self,  key ,  val ,  desc , list tags ):
        """
        _setValue_0(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] , tags: List[Union[bytes, str]] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(val, (int, float, list, bytes, str)), 'arg val wrong type'
        assert isinstance(desc, (bytes, str)), 'arg desc wrong type'
        assert isinstance(tags, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tags), 'arg tags wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        if isinstance(desc, str):
            desc = desc.encode('utf-8')
        cdef libcpp_vector[libcpp_utf8_string] v3 = tags
        self.inst.get().setValue((<libcpp_string>key), deref(ParamValue(val).inst.get()), (<libcpp_string>desc), v3)
        
    
    def _setValue_1(self,  key ,  val ,  desc ):
        """
        _setValue_1(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(val, (int, float, list, bytes, str)), 'arg val wrong type'
        assert isinstance(desc, (bytes, str)), 'arg desc wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        if isinstance(desc, str):
            desc = desc.encode('utf-8')
        self.inst.get().setValue((<libcpp_string>key), deref(ParamValue(val).inst.get()), (<libcpp_string>desc))
    
    def _setValue_2(self,  key ,  val ):
        """
        _setValue_2(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(val, (int, float, list, bytes, str)), 'arg val wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        self.inst.get().setValue((<libcpp_string>key), deref(ParamValue(val).inst.get()))
    
    def setValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] , tags: List[Union[bytes, str]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , desc: Union[bytes, str] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setValue(self, key: Union[bytes, str] , val: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], (int, float, list, bytes, str))) and (isinstance(args[2], (bytes, str))) and (isinstance(args[3], list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in args[3])):
            return self._setValue_0(*args)
        elif (len(args)==3) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], (int, float, list, bytes, str))) and (isinstance(args[2], (bytes, str))):
            return self._setValue_1(*args)
        elif (len(args)==2) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setValue_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getValue(self,  key ):
        """
        getValue(self, key: Union[bytes, str] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef _ParamValue _r = self.inst.get().getValue((<libcpp_string>key))
        cdef ParamValue _value = ParamValue.__new__(ParamValue)
        _value.inst = shared_ptr[_ParamValue](new _ParamValue(_r))
        cdef int _type = _r.valueType()
        cdef object py_result
        if _type == ValueType.STRING_VALUE:
            py_result = _value.toString()
        elif _type == ValueType.INT_VALUE:
            py_result = _value.toInt()
        elif _type == ValueType.DOUBLE_VALUE:
            py_result = _value.toDouble()
        elif _type == ValueType.INT_LIST:
            py_result = _value.toIntVector()
        elif _type == ValueType.DOUBLE_LIST:
            py_result = _value.toDoubleVector()
        elif _type == ValueType.STRING_LIST:
            py_result = _value.toStringVector()
        elif _type == ValueType.EMPTY_VALUE:
            py_result = None
        else:
            raise Exception("ParamValue instance has invalid value type %d" % _type)
        return py_result
    
    def getValueType(self,  key ):
        """
        getValueType(self, key: Union[bytes, str] ) -> int
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef _ValueType _r = self.inst.get().getValueType((<libcpp_string>key))
        py_result = <int>_r
        return py_result
    
    def getEntry(self,  in_0 ):
        """
        getEntry(self, in_0: Union[bytes, str] ) -> ParamEntry
        """
        assert isinstance(in_0, (bytes, str)), 'arg in_0 wrong type'
        if isinstance(in_0, str):
            in_0 = in_0.encode('utf-8')
        cdef _ParamEntry * _r = new _ParamEntry(self.inst.get().getEntry((<libcpp_string>in_0)))
        cdef ParamEntry py_result = ParamEntry.__new__(ParamEntry)
        py_result.inst = shared_ptr[_ParamEntry](_r)
        return py_result
    
    def exists(self,  key ):
        """
        exists(self, key: Union[bytes, str] ) -> bool
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef bool _r = self.inst.get().exists((<libcpp_string>key))
        py_result = <bool>_r
        return py_result
    
    def addTag(self,  key ,  tag ):
        """
        addTag(self, key: Union[bytes, str] , tag: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(tag, (bytes, str)), 'arg tag wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(tag, str):
            tag = tag.encode('utf-8')
        self.inst.get().addTag((<libcpp_string>key), (<libcpp_string>tag))
    
    def addTags(self,  key , list tags ):
        """
        addTags(self, key: Union[bytes, str] , tags: List[Union[bytes, str]] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(tags, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tags), 'arg tags wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef libcpp_vector[libcpp_utf8_string] v1 = tags
        self.inst.get().addTags((<libcpp_string>key), v1)
        
    
    def hasTag(self,  key ,  tag ):
        """
        hasTag(self, key: Union[bytes, str] , tag: Union[bytes, str] ) -> int
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(tag, (bytes, str)), 'arg tag wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(tag, str):
            tag = tag.encode('utf-8')
        cdef int _r = self.inst.get().hasTag((<libcpp_string>key), (<libcpp_string>tag))
        py_result = <int>_r
        return py_result
    
    def getTags(self,  key ):
        """
        getTags(self, key: Union[bytes, str] ) -> List[bytes]
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        _r = self.inst.get().getTags((<libcpp_string>key))
        cdef list py_result = _r
        return py_result
    
    def clearTags(self,  key ):
        """
        clearTags(self, key: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        self.inst.get().clearTags((<libcpp_string>key))
    
    def getDescription(self,  key ):
        """
        getDescription(self, key: Union[bytes, str] ) -> str
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef libcpp_utf8_output_string _r = self.inst.get().getDescription((<libcpp_string>key))
        py_result = _r.decode('utf-8')
        return py_result
    
    def setSectionDescription(self,  key ,  desc ):
        """
        setSectionDescription(self, key: Union[bytes, str] , desc: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(desc, (bytes, str)), 'arg desc wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(desc, str):
            desc = desc.encode('utf-8')
        self.inst.get().setSectionDescription((<libcpp_string>key), (<libcpp_string>desc))
    
    def getSectionDescription(self,  key ):
        """
        getSectionDescription(self, key: Union[bytes, str] ) -> str
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef libcpp_utf8_output_string _r = self.inst.get().getSectionDescription((<libcpp_string>key))
        py_result = _r.decode('utf-8')
        return py_result
    
    def addSection(self,  key ,  desc ):
        """
        addSection(self, key: Union[bytes, str] , desc: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(desc, (bytes, str)), 'arg desc wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        if isinstance(desc, str):
            desc = desc.encode('utf-8')
        self.inst.get().addSection((<libcpp_string>key), (<libcpp_string>desc))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def insert(self,  prefix , Param param ):
        """
        insert(self, prefix: Union[bytes, str] , param: Param ) -> None
        """
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
    
        self.inst.get().insert((<libcpp_string>prefix), (deref(param.inst.get())))
    
    def remove(self,  key ):
        """
        remove(self, key: Union[bytes, str] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        self.inst.get().remove((<libcpp_string>key))
    
    def removeAll(self,  prefix ):
        """
        removeAll(self, prefix: Union[bytes, str] ) -> None
        """
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
        self.inst.get().removeAll((<libcpp_string>prefix))
    
    def _copy_0(self,  prefix , bool in_1 ):
        """
        _copy_0(self, prefix: Union[bytes, str] , in_1: bool ) -> Param
        """
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        assert isinstance(in_1, pybool_t), 'arg in_1 wrong type'
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
    
        cdef _Param * _r = new _Param(self.inst.get().copy((<libcpp_string>prefix), (<bool>in_1)))
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def _copy_1(self,  prefix ):
        """
        _copy_1(self, prefix: Union[bytes, str] ) -> Param
        """
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
        cdef _Param * _r = new _Param(self.inst.get().copy((<libcpp_string>prefix)))
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def copy(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: copy(self, prefix: Union[bytes, str] , in_1: bool ) -> Param
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: copy(self, prefix: Union[bytes, str] ) -> Param
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], pybool_t)):
            return self._copy_0(*args)
        elif (len(args)==1) and (isinstance(args[0], (bytes, str))):
            return self._copy_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def merge(self, Param toMerge ):
        """
        merge(self, toMerge: Param ) -> None
        """
        assert isinstance(toMerge, Param), 'arg toMerge wrong type'
    
        self.inst.get().merge((deref(toMerge.inst.get())))
    
    def _setDefaults_0(self, Param defaults ,  prefix , bool showMessage ):
        """
        _setDefaults_0(self, defaults: Param , prefix: Union[bytes, str] , showMessage: bool ) -> None
        """
        assert isinstance(defaults, Param), 'arg defaults wrong type'
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        assert isinstance(showMessage, pybool_t), 'arg showMessage wrong type'
    
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
    
        self.inst.get().setDefaults((deref(defaults.inst.get())), (<libcpp_string>prefix), (<bool>showMessage))
    
    def _setDefaults_1(self, Param defaults ,  prefix ):
        """
        _setDefaults_1(self, defaults: Param , prefix: Union[bytes, str] ) -> None
        """
        assert isinstance(defaults, Param), 'arg defaults wrong type'
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
    
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
        self.inst.get().setDefaults((deref(defaults.inst.get())), (<libcpp_string>prefix))
    
    def _setDefaults_2(self, Param defaults ):
        """
        _setDefaults_2(self, defaults: Param ) -> None
        """
        assert isinstance(defaults, Param), 'arg defaults wrong type'
    
        self.inst.get().setDefaults((deref(defaults.inst.get())))
    
    def setDefaults(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setDefaults(self, defaults: Param , prefix: Union[bytes, str] , showMessage: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setDefaults(self, defaults: Param , prefix: Union[bytes, str] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setDefaults(self, defaults: Param ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], Param)) and (isinstance(args[1], (bytes, str))) and (isinstance(args[2], pybool_t)):
            return self._setDefaults_0(*args)
        elif (len(args)==2) and (isinstance(args[0], Param)) and (isinstance(args[1], (bytes, str))):
            return self._setDefaults_1(*args)
        elif (len(args)==1) and (isinstance(args[0], Param)):
            return self._setDefaults_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _checkDefaults_0(self,  name , Param defaults ,  prefix ):
        """
        _checkDefaults_0(self, name: Union[bytes, str] , defaults: Param , prefix: Union[bytes, str] ) -> None
        """
        assert isinstance(name, (bytes, str)), 'arg name wrong type'
        assert isinstance(defaults, Param), 'arg defaults wrong type'
        assert isinstance(prefix, (bytes, str)), 'arg prefix wrong type'
        if isinstance(name, str):
            name = name.encode('utf-8')
    
        if isinstance(prefix, str):
            prefix = prefix.encode('utf-8')
        self.inst.get().checkDefaults((<libcpp_string>name), (deref(defaults.inst.get())), (<libcpp_string>prefix))
    
    def _checkDefaults_1(self,  name , Param defaults ):
        """
        _checkDefaults_1(self, name: Union[bytes, str] , defaults: Param ) -> None
        """
        assert isinstance(name, (bytes, str)), 'arg name wrong type'
        assert isinstance(defaults, Param), 'arg defaults wrong type'
        if isinstance(name, str):
            name = name.encode('utf-8')
    
        self.inst.get().checkDefaults((<libcpp_string>name), (deref(defaults.inst.get())))
    
    def checkDefaults(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: checkDefaults(self, name: Union[bytes, str] , defaults: Param , prefix: Union[bytes, str] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: checkDefaults(self, name: Union[bytes, str] , defaults: Param ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], Param)) and (isinstance(args[2], (bytes, str))):
            return self._checkDefaults_0(*args)
        elif (len(args)==2) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], Param)):
            return self._checkDefaults_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getValidStrings(self,  key ):
        """
        getValidStrings(self, key: Union[bytes, str] ) -> List[Union[bytes, str]]
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        _r = self.inst.get().getValidStrings((<libcpp_string>key))
        cdef list py_result = _r
        return py_result
    
    def setValidStrings(self,  key , list strings ):
        """
        setValidStrings(self, key: Union[bytes, str] , strings: List[Union[bytes, str]] ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(strings, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in strings), 'arg strings wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
        cdef libcpp_vector[libcpp_utf8_string] v1 = strings
        self.inst.get().setValidStrings((<libcpp_string>key), v1)
        
    
    def setMinInt(self,  key ,  min ):
        """
        setMinInt(self, key: Union[bytes, str] , min: int ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(min, int), 'arg min wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        self.inst.get().setMinInt((<libcpp_string>key), (<int>min))
    
    def setMaxInt(self,  key ,  max ):
        """
        setMaxInt(self, key: Union[bytes, str] , max: int ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(max, int), 'arg max wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        self.inst.get().setMaxInt((<libcpp_string>key), (<int>max))
    
    def setMinFloat(self,  key , double min ):
        """
        setMinFloat(self, key: Union[bytes, str] , min: float ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(min, float), 'arg min wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        self.inst.get().setMinFloat((<libcpp_string>key), (<double>min))
    
    def setMaxFloat(self,  key , double max ):
        """
        setMaxFloat(self, key: Union[bytes, str] , max: float ) -> None
        """
        assert isinstance(key, (bytes, str)), 'arg key wrong type'
        assert isinstance(max, float), 'arg max wrong type'
        if isinstance(key, str):
            key = key.encode('utf-8')
    
        self.inst.get().setMaxFloat((<libcpp_string>key), (<double>max))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Param):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Param other_casted = other
        cdef Param self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    
    def initPluginParam(self, name, version):
        self.addSection(name, "Main section for Plugin '" + name + "'")
        self.setValue(name + ":version", version, "Version of Plugin '" + name + "'")
        self.addSection(name + ":1", "Instance '1' section for Plugin '" + name + "'")

        self.setValue(name + ":1:in", "", "The input file", [b"required", b"input file"])
        self.setValue(name + ":1:out", "", "The output file", [b"required", b"output file"])


    def asDict(self):
        return dict(self.items())

    def keys(self):
        keys = list()
        cdef _ParamIterator param_it = self.inst.get().begin()
        while param_it != self.inst.get().end():
            key = param_it.getName().c_str()
            keys.append(key)
            inc(param_it)
        return keys

    def items(self):
        return [(k, self[k]) for k in self.keys()]

    def values(self):
        return [self[k] for k in self.keys()]
    
    def descriptions(self):
        return [self.getDescription(k) for k in self.keys()]

    def update(self, *a):
        """
        use cases:

           p.update(dict d)
           p.update(Param p)
           p.update(Param p, int flag)
        """

        cdef Param p
        cdef int flag
        cdef bool r

        if len(a) == 1 and isinstance(a[0], dict):
            dd, = a
            for key, v in dd.items():
                self[key] = v
            return True
        elif len(a) == 1 and isinstance(a[0], Param):
            p, = a
            r = self.inst.get().update(<_Param>deref(p.inst.get()))
            return r
        elif len(a) == 2 and isinstance(a[0], Param) and isinstance(a[1], int):
            p, flag = a
            r = self.inst.get().update(<_Param>deref(p.inst.get()), <bool> flag)
            return r
        else:
            raise Exception("can not handle parameters of type %s" % (map(type, a)))

    def get(self, bytes key, default=None):
        if self.exists(key):
            return self.getValue(key)
        return default

    def __getitem__(self, bytes key):
        return self.getValue(key)

    def __setitem__(self, bytes key, value):
        tags = self.getTags(key)
        desc = self.getDescription(key)
        self.setValue(key, value, desc, tags)
        
    def __str__(self):
        return str(list(zip([k.decode() for k in self.keys()], self.values(), self.descriptions())))

    def __repr__(self):
        return self.__str__() 

cdef class Ratio:
    """
    Cython implementation of _Ratio

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Ratio.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property ratio_value_:
        def __set__(self, double ratio_value_):
        
            self.inst.get().ratio_value_ = (<double>ratio_value_)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ratio_value_
            py_result = <double>_r
            return py_result
    
    property denominator_ref_:
        def __set__(self,  denominator_ref_):
        
            self.inst.get().denominator_ref_ = deref((convString(denominator_ref_)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().denominator_ref_
            py_result = convOutputString(_r)
            return py_result
    
    property numerator_ref_:
        def __set__(self,  numerator_ref_):
        
            self.inst.get().numerator_ref_ = deref((convString(numerator_ref_)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().numerator_ref_
            py_result = convOutputString(_r)
            return py_result
    
    property description_:
        def __set__(self, list description_):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in description_:
               v0.push_back(_String(<char *>item0))
            self.inst.get().description_ = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().description_
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Ratio](new _Ratio())
    
    def _init_1(self, Ratio rhs ):
        """
        _init_1(self, rhs: Ratio ) -> None
        """
        assert isinstance(rhs, Ratio), 'arg rhs wrong type'
    
        self.inst = shared_ptr[_Ratio](new _Ratio((deref(rhs.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, rhs: Ratio ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Ratio)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class TraceInfo:
    """
    Cython implementation of _TraceInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TraceInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property name:
        def __set__(self, bytes name):
        
            self.inst.get().name = (<libcpp_string>name)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().name
            py_result = <libcpp_string>_r
            return py_result
    
    property description:
        def __set__(self, bytes description):
        
            self.inst.get().description = (<libcpp_string>description)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().description
            py_result = <libcpp_string>_r
            return py_result
    
    property opened:
        def __set__(self, bool opened):
        
            self.inst.get().opened = (<bool>opened)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().opened
            py_result = <bool>_r
            return py_result
    
    def _init_0(self,  n ,  d , bool o ):
        """
        _init_0(self, n: Union[bytes, str] , d: Union[bytes, str] , o: bool ) -> None
        """
        assert isinstance(n, (bytes, str)), 'arg n wrong type'
        assert isinstance(d, (bytes, str)), 'arg d wrong type'
        assert isinstance(o, pybool_t), 'arg o wrong type'
        if isinstance(n, str):
            n = n.encode('utf-8')
        if isinstance(d, str):
            d = d.encode('utf-8')
    
        self.inst = shared_ptr[_TraceInfo](new _TraceInfo((<libcpp_string>n), (<libcpp_string>d), (<bool>o)))
    
    def _init_1(self, TraceInfo in_0 ):
        """
        _init_1(self, in_0: TraceInfo ) -> None
        """
        assert isinstance(in_0, TraceInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TraceInfo](new _TraceInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, n: Union[bytes, str] , d: Union[bytes, str] , o: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TraceInfo ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], (bytes, str))) and (isinstance(args[1], (bytes, str))) and (isinstance(args[2], pybool_t)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TraceInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class XTandemXMLFile:
    """
    Cython implementation of _XTandemXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XTandemXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_XTandemXMLFile](new _XTandemXMLFile())
    
    def load(self,  filename , ProteinIdentification protein_identification , PeptideIdentificationList id_data , ModificationDefinitionsSet mod_def_set ):
        """
        load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , mod_def_set: ModificationDefinitionsSet ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(id_data, PeptideIdentificationList), 'arg id_data wrong type'
        assert isinstance(mod_def_set, ModificationDefinitionsSet), 'arg mod_def_set wrong type'
    
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(protein_identification.inst.get())), (deref(id_data.inst.get())), (deref(mod_def_set.inst.get()))) 
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
