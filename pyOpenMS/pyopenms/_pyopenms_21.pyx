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
def __static_PrecursorCorrection_correctToHighestIntensityMS1Peak(MSExperiment exp , double mz_tolerance , bool ppm , list delta_mzs , list mzs , list rts ):
    """
    __static_PrecursorCorrection_correctToHighestIntensityMS1Peak(exp: MSExperiment , mz_tolerance: float , ppm: bool , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> Set[int]
    """
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    assert isinstance(mz_tolerance, float), 'arg mz_tolerance wrong type'
    assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
    assert isinstance(delta_mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in delta_mzs), 'arg delta_mzs wrong type'
    assert isinstance(mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in mzs), 'arg mzs wrong type'
    assert isinstance(rts, list) and all(isinstance(elemt_rec, float) for elemt_rec in rts), 'arg rts wrong type'



    cdef libcpp_vector[double] v3 = delta_mzs
    cdef libcpp_vector[double] v4 = mzs
    cdef libcpp_vector[double] v5 = rts
    _r = _correctToHighestIntensityMS1Peak_PrecursorCorrection((deref(exp.inst.get())), (<double>mz_tolerance), (<bool>ppm), v3, v4, v5)
    rts[:] = v5
    mzs[:] = v4
    delta_mzs[:] = v3
    cdef set py_result = _r
    return py_result

def __static_PrecursorCorrection_correctToNearestFeature(FeatureMap features , MSExperiment exp , double rt_tolerance_s , double mz_tolerance , bool ppm , bool believe_charge , bool keep_original , bool all_matching_features ,  max_trace ,  debug_level ):
    """
    __static_PrecursorCorrection_correctToNearestFeature(features: FeatureMap , exp: MSExperiment , rt_tolerance_s: float , mz_tolerance: float , ppm: bool , believe_charge: bool , keep_original: bool , all_matching_features: bool , max_trace: int , debug_level: int ) -> Set[int]
    """
    assert isinstance(features, FeatureMap), 'arg features wrong type'
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    assert isinstance(rt_tolerance_s, float), 'arg rt_tolerance_s wrong type'
    assert isinstance(mz_tolerance, float), 'arg mz_tolerance wrong type'
    assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
    assert isinstance(believe_charge, pybool_t), 'arg believe_charge wrong type'
    assert isinstance(keep_original, pybool_t), 'arg keep_original wrong type'
    assert isinstance(all_matching_features, pybool_t), 'arg all_matching_features wrong type'
    assert isinstance(max_trace, int), 'arg max_trace wrong type'
    assert isinstance(debug_level, int), 'arg debug_level wrong type'










    _r = _correctToNearestFeature_PrecursorCorrection((deref(features.inst.get())), (deref(exp.inst.get())), (<double>rt_tolerance_s), (<double>mz_tolerance), (<bool>ppm), (<bool>believe_charge), (<bool>keep_original), (<bool>all_matching_features), (<int>max_trace), (<int>debug_level))
    cdef set py_result = _r
    return py_result

def __static_PrecursorCorrection_correctToNearestMS1Peak(MSExperiment exp , double mz_tolerance , bool ppm , list delta_mzs , list mzs , list rts ):
    """
    __static_PrecursorCorrection_correctToNearestMS1Peak(exp: MSExperiment , mz_tolerance: float , ppm: bool , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> Set[int]
    """
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    assert isinstance(mz_tolerance, float), 'arg mz_tolerance wrong type'
    assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
    assert isinstance(delta_mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in delta_mzs), 'arg delta_mzs wrong type'
    assert isinstance(mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in mzs), 'arg mzs wrong type'
    assert isinstance(rts, list) and all(isinstance(elemt_rec, float) for elemt_rec in rts), 'arg rts wrong type'



    cdef libcpp_vector[double] v3 = delta_mzs
    cdef libcpp_vector[double] v4 = mzs
    cdef libcpp_vector[double] v5 = rts
    _r = _correctToNearestMS1Peak_PrecursorCorrection((deref(exp.inst.get())), (<double>mz_tolerance), (<bool>ppm), v3, v4, v5)
    rts[:] = v5
    mzs[:] = v4
    delta_mzs[:] = v3
    cdef set py_result = _r
    return py_result

def __static_PrecursorCorrection_getPrecursors(MSExperiment exp , list precursors , list precursors_rt , list precursor_scan_index ):
    """
    __static_PrecursorCorrection_getPrecursors(exp: MSExperiment , precursors: List[Precursor] , precursors_rt: List[float] , precursor_scan_index: List[int] ) -> None
    """
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    assert isinstance(precursors, list) and all(isinstance(elemt_rec, Precursor) for elemt_rec in precursors), 'arg precursors wrong type'
    assert isinstance(precursors_rt, list) and all(isinstance(elemt_rec, float) for elemt_rec in precursors_rt), 'arg precursors_rt wrong type'
    assert isinstance(precursor_scan_index, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in precursor_scan_index), 'arg precursor_scan_index wrong type'

    cdef libcpp_vector[_Precursor] * v1 = new libcpp_vector[_Precursor]()
    cdef Precursor item1
    for item1 in precursors:
        v1.push_back(deref(item1.inst.get()))
    cdef libcpp_vector[double] v2 = precursors_rt
    cdef libcpp_vector[size_t] v3 = precursor_scan_index
    _getPrecursors_PrecursorCorrection((deref(exp.inst.get())), deref(v1), v2, v3)
    precursor_scan_index[:] = v3
    precursors_rt[:] = v2
    cdef libcpp_vector[_Precursor].iterator it_precursors = v1.begin()
    replace_0 = []
    while it_precursors != v1.end():
        item1 = Precursor.__new__(Precursor)
        item1.inst = shared_ptr[_Precursor](new _Precursor(deref(it_precursors)))
        replace_0.append(item1)
        inc(it_precursors)
    precursors[:] = replace_0
    del v1

def __static_PrecursorCorrection_writeHist( out_csv , list delta_mzs , list mzs , list rts ):
    """
    __static_PrecursorCorrection_writeHist(out_csv: String , delta_mzs: List[float] , mzs: List[float] , rts: List[float] ) -> None
    """
    assert isinstance(out_csv, String), 'arg out_csv wrong type'
    assert isinstance(delta_mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in delta_mzs), 'arg delta_mzs wrong type'
    assert isinstance(mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in mzs), 'arg mzs wrong type'
    assert isinstance(rts, list) and all(isinstance(elemt_rec, float) for elemt_rec in rts), 'arg rts wrong type'

    cdef libcpp_vector[double] v1 = delta_mzs
    cdef libcpp_vector[double] v2 = mzs
    cdef libcpp_vector[double] v3 = rts
    _writeHist_PrecursorCorrection(deref((<String>out_csv).inst.get()), v1, v2, v3)
    rts[:] = v3
    mzs[:] = v2
    delta_mzs[:] = v1 

cdef class DTA2DFile:
    """
    Cython implementation of _DTA2DFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DTA2DFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DTA2DFile rv = DTA2DFile.__new__(DTA2DFile)
       rv.inst = shared_ptr[_DTA2DFile](new _DTA2DFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DTA2DFile rv = DTA2DFile.__new__(DTA2DFile)
       rv.inst = shared_ptr[_DTA2DFile](new _DTA2DFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DTA2DFile](new _DTA2DFile())
    
    def _init_1(self, DTA2DFile in_0 ):
        """
        _init_1(self, in_0: DTA2DFile ) -> None
        """
        assert isinstance(in_0, DTA2DFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DTA2DFile](new _DTA2DFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DTA2DFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DTA2DFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def storeTIC(self,  filename , MSExperiment peakmap ):
        """
        storeTIC(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(peakmap, MSExperiment), 'arg peakmap wrong type'
    
    
        self.inst.get().storeTIC(deref((convString(filename)).get()), (deref(peakmap.inst.get())))
    
    def store(self,  filename , MSExperiment peakmap ):
        """
        store(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(peakmap, MSExperiment), 'arg peakmap wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(peakmap.inst.get())))
    
    def load(self,  filename , MSExperiment peakmap ):
        """
        load(self, filename: Union[bytes, str, String] , peakmap: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(peakmap, MSExperiment), 'arg peakmap wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(peakmap.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
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

cdef class DTAFile:
    """
    Cython implementation of _DTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DTAFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DTAFile rv = DTAFile.__new__(DTAFile)
       rv.inst = shared_ptr[_DTAFile](new _DTAFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DTAFile rv = DTAFile.__new__(DTAFile)
       rv.inst = shared_ptr[_DTAFile](new _DTAFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DTAFile](new _DTAFile())
    
    def _init_1(self, DTAFile in_0 ):
        """
        _init_1(self, in_0: DTAFile ) -> None
        """
        assert isinstance(in_0, DTAFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DTAFile](new _DTAFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DTAFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DTAFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSSpectrum spectrum ):
        """
        load(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(spectrum.inst.get())))
    
    def store(self,  filename , MSSpectrum spectrum ):
        """
        store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(spectrum.inst.get()))) 

cdef class FeatureHandle:
    """
    Cython implementation of _FeatureHandle

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureHandle.html>`_
      -- Inherits from ['Peak2D', 'UniqueIdInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureHandle rv = FeatureHandle.__new__(FeatureHandle)
       rv.inst = shared_ptr[_FeatureHandle](new _FeatureHandle(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureHandle rv = FeatureHandle.__new__(FeatureHandle)
       rv.inst = shared_ptr[_FeatureHandle](new _FeatureHandle(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Representation of a Peak2D, RichPeak2D or Feature
        """
        self.inst = shared_ptr[_FeatureHandle](new _FeatureHandle())
    
    def _init_1(self, FeatureHandle in_0 ):
        """
        _init_1(self, in_0: FeatureHandle ) -> None
        """
        assert isinstance(in_0, FeatureHandle), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureHandle](new _FeatureHandle((deref(in_0.inst.get()))))
    
    def _init_2(self,  map_index , Peak2D point ,  element_index ):
        """
        _init_2(self, map_index: int , point: Peak2D , element_index: int ) -> None
        """
        assert isinstance(map_index, int) and map_index >= 0, 'arg map_index wrong type'
        assert isinstance(point, Peak2D), 'arg point wrong type'
        assert isinstance(element_index, int) and element_index >= 0, 'arg element_index wrong type'
    
    
    
        self.inst = shared_ptr[_FeatureHandle](new _FeatureHandle((<uint64_t>map_index), (deref(point.inst.get())), (<uint64_t>element_index)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Representation of a Peak2D, RichPeak2D or Feature

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureHandle ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, map_index: int , point: Peak2D , element_index: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureHandle)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], Peak2D)) and (isinstance(args[2], int) and args[2] >= 0):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMapIndex(self):
        """
        getMapIndex(self) -> int
        Returns the map index
        """
        cdef uint64_t _r = self.inst.get().getMapIndex()
        py_result = <uint64_t>_r
        return py_result
    
    def setMapIndex(self,  i ):
        """
        setMapIndex(self, i: int ) -> None
        Sets the map index
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        self.inst.get().setMapIndex((<uint64_t>i))
    
    def setCharge(self,  charge ):
        """
        setCharge(self, charge: int ) -> None
        Sets the charge
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setCharge((<int>charge))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns the charge
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def setWidth(self, float width ):
        """
        setWidth(self, width: float ) -> None
        Sets the width (FWHM)
        """
        assert isinstance(width, float), 'arg width wrong type'
    
        self.inst.get().setWidth((<float>width))
    
    def getWidth(self):
        """
        getWidth(self) -> float
        Returns the width (FWHM)
        """
        cdef float _r = self.inst.get().getWidth()
        py_result = <float>_r
        return py_result
    
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
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, FeatureHandle):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef FeatureHandle other_casted = other
        cdef FeatureHandle self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class IncludeExcludeTarget:
    """
    Cython implementation of _IncludeExcludeTarget

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IncludeExcludeTarget.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IncludeExcludeTarget rv = IncludeExcludeTarget.__new__(IncludeExcludeTarget)
       rv.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IncludeExcludeTarget rv = IncludeExcludeTarget.__new__(IncludeExcludeTarget)
       rv.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        This class stores a SRM/MRM transition
        """
        self.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget())
    
    def _init_1(self, IncludeExcludeTarget in_0 ):
        """
        _init_1(self, in_0: IncludeExcludeTarget ) -> None
        """
        assert isinstance(in_0, IncludeExcludeTarget), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        This class stores a SRM/MRM transition

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IncludeExcludeTarget ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IncludeExcludeTarget)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setPeptideRef(self,  peptide_ref ):
        """
        setPeptideRef(self, peptide_ref: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(peptide_ref, str) or isinstance(peptide_ref, bytes) or isinstance(peptide_ref, String)), 'arg peptide_ref wrong type'
    
        self.inst.get().setPeptideRef(deref((convString(peptide_ref)).get()))
    
    def getPeptideRef(self):
        """
        getPeptideRef(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getPeptideRef()
        py_result = convOutputString(_r)
        return py_result
    
    def setNuctideRef(self,  nuctide_ref ):
        """
        setNuctideRef(self, nuctide_ref: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(nuctide_ref, str) or isinstance(nuctide_ref, bytes) or isinstance(nuctide_ref, String)), 'arg nuctide_ref wrong type'
    
        self.inst.get().setNuctideRef(deref((convString(nuctide_ref)).get()))
    
    def getNuctideRef(self):
        """
        getNuctideRef(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getNuctideRef()
        py_result = convOutputString(_r)
        return py_result
    
    def setCompoundRef(self,  compound_ref ):
        """
        setCompoundRef(self, compound_ref: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(compound_ref, str) or isinstance(compound_ref, bytes) or isinstance(compound_ref, String)), 'arg compound_ref wrong type'
    
        self.inst.get().setCompoundRef(deref((convString(compound_ref)).get()))
    
    def getCompoundRef(self):
        """
        getCompoundRef(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getCompoundRef()
        py_result = convOutputString(_r)
        return py_result
    
    def setPrecursorMZ(self, double mz ):
        """
        setPrecursorMZ(self, mz: float ) -> None
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        self.inst.get().setPrecursorMZ((<double>mz))
    
    def getPrecursorMZ(self):
        """
        getPrecursorMZ(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorMZ()
        py_result = <double>_r
        return py_result
    
    def setPrecursorCVTermList(self, CVTermList list_ ):
        """
        setPrecursorCVTermList(self, list_: CVTermList ) -> None
        """
        assert isinstance(list_, CVTermList), 'arg list_ wrong type'
    
        self.inst.get().setPrecursorCVTermList((deref(list_.inst.get())))
    
    def addPrecursorCVTerm(self, CVTerm cv_term ):
        """
        addPrecursorCVTerm(self, cv_term: CVTerm ) -> None
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().addPrecursorCVTerm((deref(cv_term.inst.get())))
    
    def getPrecursorCVTermList(self):
        """
        getPrecursorCVTermList(self) -> CVTermList
        """
        cdef _CVTermList * _r = new _CVTermList(self.inst.get().getPrecursorCVTermList())
        cdef CVTermList py_result = CVTermList.__new__(CVTermList)
        py_result.inst = shared_ptr[_CVTermList](_r)
        return py_result
    
    def setProductMZ(self, double mz ):
        """
        setProductMZ(self, mz: float ) -> None
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        self.inst.get().setProductMZ((<double>mz))
    
    def getProductMZ(self):
        """
        getProductMZ(self) -> float
        """
        cdef double _r = self.inst.get().getProductMZ()
        py_result = <double>_r
        return py_result
    
    def setProductCVTermList(self, CVTermList list_ ):
        """
        setProductCVTermList(self, list_: CVTermList ) -> None
        """
        assert isinstance(list_, CVTermList), 'arg list_ wrong type'
    
        self.inst.get().setProductCVTermList((deref(list_.inst.get())))
    
    def addProductCVTerm(self, CVTerm cv_term ):
        """
        addProductCVTerm(self, cv_term: CVTerm ) -> None
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().addProductCVTerm((deref(cv_term.inst.get())))
    
    def getProductCVTermList(self):
        """
        getProductCVTermList(self) -> CVTermList
        """
        cdef _CVTermList * _r = new _CVTermList(self.inst.get().getProductCVTermList())
        cdef CVTermList py_result = CVTermList.__new__(CVTermList)
        py_result.inst = shared_ptr[_CVTermList](_r)
        return py_result
    
    def setInterpretations(self, list interpretations ):
        """
        setInterpretations(self, interpretations: List[CVTermList] ) -> None
        """
        assert isinstance(interpretations, list) and all(isinstance(elemt_rec, CVTermList) for elemt_rec in interpretations), 'arg interpretations wrong type'
        cdef libcpp_vector[_CVTermList] * v0 = new libcpp_vector[_CVTermList]()
        cdef CVTermList item0
        for item0 in interpretations:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setInterpretations(deref(v0))
        cdef libcpp_vector[_CVTermList].iterator it_interpretations = v0.begin()
        replace_0 = []
        while it_interpretations != v0.end():
            item0 = CVTermList.__new__(CVTermList)
            item0.inst = shared_ptr[_CVTermList](new _CVTermList(deref(it_interpretations)))
            replace_0.append(item0)
            inc(it_interpretations)
        interpretations[:] = replace_0
        del v0
    
    def getInterpretations(self):
        """
        getInterpretations(self) -> List[CVTermList]
        """
        _r = self.inst.get().getInterpretations()
        py_result = []
        cdef libcpp_vector[_CVTermList].iterator it__r = _r.begin()
        cdef CVTermList item_py_result
        while it__r != _r.end():
           item_py_result = CVTermList.__new__(CVTermList)
           item_py_result.inst = shared_ptr[_CVTermList](new _CVTermList(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addInterpretation(self, CVTermList interpretation ):
        """
        addInterpretation(self, interpretation: CVTermList ) -> None
        """
        assert isinstance(interpretation, CVTermList), 'arg interpretation wrong type'
    
        self.inst.get().addInterpretation((deref(interpretation.inst.get())))
    
    def setConfigurations(self, list configuration ):
        """
        setConfigurations(self, configuration: List[Configuration] ) -> None
        """
        assert isinstance(configuration, list) and all(isinstance(elemt_rec, Configuration) for elemt_rec in configuration), 'arg configuration wrong type'
        cdef libcpp_vector[_Configuration] * v0 = new libcpp_vector[_Configuration]()
        cdef Configuration item0
        for item0 in configuration:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setConfigurations(deref(v0))
        cdef libcpp_vector[_Configuration].iterator it_configuration = v0.begin()
        replace_0 = []
        while it_configuration != v0.end():
            item0 = Configuration.__new__(Configuration)
            item0.inst = shared_ptr[_Configuration](new _Configuration(deref(it_configuration)))
            replace_0.append(item0)
            inc(it_configuration)
        configuration[:] = replace_0
        del v0
    
    def getConfigurations(self):
        """
        getConfigurations(self) -> List[Configuration]
        """
        _r = self.inst.get().getConfigurations()
        py_result = []
        cdef libcpp_vector[_Configuration].iterator it__r = _r.begin()
        cdef Configuration item_py_result
        while it__r != _r.end():
           item_py_result = Configuration.__new__(Configuration)
           item_py_result.inst = shared_ptr[_Configuration](new _Configuration(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addConfiguration(self, Configuration configuration ):
        """
        addConfiguration(self, configuration: Configuration ) -> None
        """
        assert isinstance(configuration, Configuration), 'arg configuration wrong type'
    
        self.inst.get().addConfiguration((deref(configuration.inst.get())))
    
    def setPrediction(self, CVTermList prediction ):
        """
        setPrediction(self, prediction: CVTermList ) -> None
        """
        assert isinstance(prediction, CVTermList), 'arg prediction wrong type'
    
        self.inst.get().setPrediction((deref(prediction.inst.get())))
    
    def addPredictionTerm(self, CVTerm prediction ):
        """
        addPredictionTerm(self, prediction: CVTerm ) -> None
        """
        assert isinstance(prediction, CVTerm), 'arg prediction wrong type'
    
        self.inst.get().addPredictionTerm((deref(prediction.inst.get())))
    
    def getPrediction(self):
        """
        getPrediction(self) -> CVTermList
        """
        cdef _CVTermList * _r = new _CVTermList(self.inst.get().getPrediction())
        cdef CVTermList py_result = CVTermList.__new__(CVTermList)
        py_result.inst = shared_ptr[_CVTermList](_r)
        return py_result
    
    def setRetentionTime(self, RetentionTime rt ):
        """
        setRetentionTime(self, rt: RetentionTime ) -> None
        """
        assert isinstance(rt, RetentionTime), 'arg rt wrong type'
    
        self.inst.get().setRetentionTime((deref(rt.inst.get())))
    
    def getRetentionTime(self):
        """
        getRetentionTime(self) -> RetentionTime
        """
        cdef _RetentionTime * _r = new _RetentionTime(self.inst.get().getRetentionTime())
        cdef RetentionTime py_result = RetentionTime.__new__(RetentionTime)
        py_result.inst = shared_ptr[_RetentionTime](_r)
        return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        """
        assert isinstance(terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in terms), 'arg terms wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in terms:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setCVTerms(deref(v0))
        cdef libcpp_vector[_CVTerm].iterator it_terms = v0.begin()
        replace_0 = []
        while it_terms != v0.end():
            item0 = CVTerm.__new__(CVTerm)
            item0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(it_terms)))
            replace_0.append(item0)
            inc(it_terms)
        terms[:] = replace_0
        del v0
    
    def replaceCVTerm(self, CVTerm term ):
        """
        replaceCVTerm(self, term: CVTerm ) -> None
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def _replaceCVTerms_0(self, list cv_terms ,  accession ):
        """
        _replaceCVTerms_0(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
    def _replaceCVTerms_1(self, dict cv_term_map ):
        """
        _replaceCVTerms_1(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None
        """
        assert isinstance(cv_term_map, dict) and all(isinstance(k, bytes) for k in cv_term_map.keys()) and all(isinstance(v, list) for v in cv_term_map.values()) and all(isinstance(vi, CVTerm) for v in cv_term_map.values() for vi in
          v), 'arg cv_term_map wrong type'
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]] _map_0
        cdef libcpp_vector[_CVTerm] _v_vec_0
        cdef _String _k_str_0
        cdef CVTerm _v_i_0
        for k, v in cv_term_map.items():
            _v_vec_0.clear()
            for _v_i_0 in v:
                _v_vec_0.push_back(deref(_v_i_0.inst.get()))
            _map_0[_String(<char *>k)] = _v_vec_0
        self.inst.get().replaceCVTerms(_map_0)
    
    def replaceCVTerms(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: replaceCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in args[0])) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._replaceCVTerms_0(*args)
        elif (len(args)==1) and (isinstance(args[0], dict) and all(isinstance(k, bytes) for k in args[0].keys()) and all(isinstance(v, list) for v in args[0].values()) and all(isinstance(vi, CVTerm) for v in args[0].values() for vi in
          v)):
            return self._replaceCVTerms_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773790660001753366646 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773790660001753366646
        cdef CVTerm item_1268773790660001753366646
        cdef bytes inner_key_1268773790660001753366646
        cdef list inner_values_1268773790660001753366646
        while outer_it_1268773790660001753366646 != _r.end():
           inner_key_1268773790660001753366646 = deref(outer_it_1268773790660001753366646).first.c_str()
           inner_values_1268773790660001753366646 = []
           inner_it_1268773790660001753366646 = deref(outer_it_1268773790660001753366646).second.begin()
           while inner_it_1268773790660001753366646 != deref(outer_it_1268773790660001753366646).second.end():
               item_1268773790660001753366646 = CVTerm.__new__(CVTerm)
               item_1268773790660001753366646.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773790660001753366646)))
               inner_values_1268773790660001753366646.append(item_1268773790660001753366646)
               inc(inner_it_1268773790660001753366646)
           py_result[inner_key_1268773790660001753366646] = inner_values_1268773790660001753366646
           inc(outer_it_1268773790660001753366646)
        return py_result
    
    def addCVTerm(self, CVTerm term ):
        """
        addCVTerm(self, term: CVTerm ) -> None
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().addCVTerm((deref(term.inst.get())))
    
    def hasCVTerm(self,  accession ):
        """
        hasCVTerm(self, accession: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
    
        cdef bool _r = self.inst.get().hasCVTerm(deref((convString(accession)).get()))
        py_result = <bool>_r
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, IncludeExcludeTarget):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IncludeExcludeTarget other_casted = other
        cdef IncludeExcludeTarget self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class IsobaricChannelExtractor:
    """
    Cython implementation of _IsobaricChannelExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricChannelExtractor.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsobaricChannelExtractor rv = IsobaricChannelExtractor.__new__(IsobaricChannelExtractor)
       rv.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricChannelExtractor rv = IsobaricChannelExtractor.__new__(IsobaricChannelExtractor)
       rv.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor(deref(self.inst.get())))
       return rv
    
    def _init_0(self, IsobaricChannelExtractor in_0 ):
        """
        _init_0(self, in_0: IsobaricChannelExtractor ) -> None
        """
        assert isinstance(in_0, IsobaricChannelExtractor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor((deref(in_0.inst.get()))))
    
    def _init_1(self, ItraqEightPlexQuantitationMethod quant_method ):
        """
        _init_1(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqEightPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor((quant_method.inst.get())))
    
    def _init_2(self, ItraqFourPlexQuantitationMethod quant_method ):
        """
        _init_2(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqFourPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor((quant_method.inst.get())))
    
    def _init_3(self, TMTSixPlexQuantitationMethod quant_method ):
        """
        _init_3(self, quant_method: TMTSixPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTSixPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor((quant_method.inst.get())))
    
    def _init_4(self, TMTTenPlexQuantitationMethod quant_method ):
        """
        _init_4(self, quant_method: TMTTenPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTTenPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelExtractor](new _IsobaricChannelExtractor((quant_method.inst.get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricChannelExtractor ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: TMTSixPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: TMTTenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], IsobaricChannelExtractor)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqEightPlexQuantitationMethod)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqFourPlexQuantitationMethod)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTSixPlexQuantitationMethod)):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTTenPlexQuantitationMethod)):
             self._init_4(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractChannels(self, MSExperiment ms_exp_data , ConsensusMap consensus_map ):
        """
        extractChannels(self, ms_exp_data: MSExperiment , consensus_map: ConsensusMap ) -> None
        """
        assert isinstance(ms_exp_data, MSExperiment), 'arg ms_exp_data wrong type'
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
    
    
        self.inst.get().extractChannels((deref(ms_exp_data.inst.get())), (deref(consensus_map.inst.get())))
    
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

cdef class KDTreeFeatureMaps:
    """
    Cython implementation of _KDTreeFeatureMaps

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KDTreeFeatureMaps.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        Stores a set of features, together with a 2D tree for fast search
        """
        self.inst = shared_ptr[_KDTreeFeatureMaps](new _KDTreeFeatureMaps())
    
    def _init_1(self, list maps , Param param ):
        """
        _init_1(self, maps: List[FeatureMap] , param: Param ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in maps), 'arg maps wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
        cdef libcpp_vector[_FeatureMap] * v0 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_KDTreeFeatureMaps](new _KDTreeFeatureMaps(deref(v0), (deref(param.inst.get()))))
        cdef libcpp_vector[_FeatureMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = FeatureMap.__new__(FeatureMap)
            item0.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def _init_2(self, list maps , Param param ):
        """
        _init_2(self, maps: List[ConsensusMap] , param: Param ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in maps), 'arg maps wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_KDTreeFeatureMaps](new _KDTreeFeatureMaps(deref(v0), (deref(param.inst.get()))))
        cdef libcpp_vector[_ConsensusMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Stores a set of features, together with a 2D tree for fast search

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, maps: List[FeatureMap] , param: Param ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, maps: List[ConsensusMap] , param: Param ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in args[0])) and (isinstance(args[1], Param)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in args[0])) and (isinstance(args[1], Param)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _addMaps_0(self, list maps ):
        """
        _addMaps_0(self, maps: List[FeatureMap] ) -> None
        Add `maps` and balance kd-tree
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in maps), 'arg maps wrong type'
        cdef libcpp_vector[_FeatureMap] * v0 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().addMaps(deref(v0))
        cdef libcpp_vector[_FeatureMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = FeatureMap.__new__(FeatureMap)
            item0.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def _addMaps_1(self, list maps ):
        """
        _addMaps_1(self, maps: List[ConsensusMap] ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in maps), 'arg maps wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().addMaps(deref(v0))
        cdef libcpp_vector[_ConsensusMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def addMaps(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: addMaps(self, maps: List[FeatureMap] ) -> None
          :noindex:
        
        Add `maps` and balance kd-tree

        
        .. rubric:: Overload:
        .. py:function:: addMaps(self, maps: List[ConsensusMap] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in args[0])):
            return self._addMaps_0(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in args[0])):
            return self._addMaps_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def rt(self,  i ):
        """
        rt(self, i: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef double _r = self.inst.get().rt((<size_t>i))
        py_result = <double>_r
        return py_result
    
    def mz(self,  i ):
        """
        mz(self, i: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef double _r = self.inst.get().mz((<size_t>i))
        py_result = <double>_r
        return py_result
    
    def intensity(self,  i ):
        """
        intensity(self, i: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef float _r = self.inst.get().intensity((<size_t>i))
        py_result = <float>_r
        return py_result
    
    def charge(self,  i ):
        """
        charge(self, i: int ) -> int
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef int _r = self.inst.get().charge((<size_t>i))
        py_result = <int>_r
        return py_result
    
    def mapIndex(self,  i ):
        """
        mapIndex(self, i: int ) -> int
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef size_t _r = self.inst.get().mapIndex((<size_t>i))
        py_result = <size_t>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def treeSize(self):
        """
        treeSize(self) -> int
        """
        cdef size_t _r = self.inst.get().treeSize()
        py_result = <size_t>_r
        return py_result
    
    def numMaps(self):
        """
        numMaps(self) -> int
        """
        cdef size_t _r = self.inst.get().numMaps()
        py_result = <size_t>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def optimizeTree(self):
        """
        optimizeTree(self) -> None
        """
        self.inst.get().optimizeTree()
    
    def getNeighborhood(self,  index , list result_indices , double rt_tol , double mz_tol , bool mz_ppm , bool include_features_from_same_map , double max_pairwise_log_fc ):
        """
        getNeighborhood(self, index: int , result_indices: List[int] , rt_tol: float , mz_tol: float , mz_ppm: bool , include_features_from_same_map: bool , max_pairwise_log_fc: float ) -> None
        Fill `result` with indices of all features compatible (wrt. RT, m/z, map index) to the feature with `index`
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(result_indices, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in result_indices), 'arg result_indices wrong type'
        assert isinstance(rt_tol, float), 'arg rt_tol wrong type'
        assert isinstance(mz_tol, float), 'arg mz_tol wrong type'
        assert isinstance(mz_ppm, pybool_t), 'arg mz_ppm wrong type'
        assert isinstance(include_features_from_same_map, pybool_t), 'arg include_features_from_same_map wrong type'
        assert isinstance(max_pairwise_log_fc, float), 'arg max_pairwise_log_fc wrong type'
    
        cdef libcpp_vector[size_t] v1 = result_indices
    
    
    
    
    
        self.inst.get().getNeighborhood((<size_t>index), v1, (<double>rt_tol), (<double>mz_tol), (<bool>mz_ppm), (<bool>include_features_from_same_map), (<double>max_pairwise_log_fc))
        result_indices[:] = v1
    
    def queryRegion(self, double rt_low , double rt_high , double mz_low , double mz_high , list result_indices ,  ignored_map_index ):
        """
        queryRegion(self, rt_low: float , rt_high: float , mz_low: float , mz_high: float , result_indices: List[int] , ignored_map_index: int ) -> None
        """
        assert isinstance(rt_low, float), 'arg rt_low wrong type'
        assert isinstance(rt_high, float), 'arg rt_high wrong type'
        assert isinstance(mz_low, float), 'arg mz_low wrong type'
        assert isinstance(mz_high, float), 'arg mz_high wrong type'
        assert isinstance(result_indices, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in result_indices), 'arg result_indices wrong type'
        assert isinstance(ignored_map_index, int) and ignored_map_index >= 0, 'arg ignored_map_index wrong type'
    
    
    
    
        cdef libcpp_vector[size_t] v4 = result_indices
    
        self.inst.get().queryRegion((<double>rt_low), (<double>rt_high), (<double>mz_low), (<double>mz_high), v4, (<size_t>ignored_map_index))
        result_indices[:] = v4
    
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

cdef class KDTreeFeatureNode:
    """
    Cython implementation of _KDTreeFeatureNode

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KDTreeFeatureNode.html>`_

    A node of the kD-tree with pointer to corresponding data and index
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef KDTreeFeatureNode rv = KDTreeFeatureNode.__new__(KDTreeFeatureNode)
       rv.inst = shared_ptr[_KDTreeFeatureNode](new _KDTreeFeatureNode(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef KDTreeFeatureNode rv = KDTreeFeatureNode.__new__(KDTreeFeatureNode)
       rv.inst = shared_ptr[_KDTreeFeatureNode](new _KDTreeFeatureNode(deref(self.inst.get())))
       return rv
    
    def _init_0(self, KDTreeFeatureNode in_0 ):
        """
        _init_0(self, in_0: KDTreeFeatureNode ) -> None
        """
        assert isinstance(in_0, KDTreeFeatureNode), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_KDTreeFeatureNode](new _KDTreeFeatureNode((deref(in_0.inst.get()))))
    
    def _init_1(self, KDTreeFeatureMaps data ,  idx ):
        """
        _init_1(self, data: KDTreeFeatureMaps , idx: int ) -> None
        """
        assert isinstance(data, KDTreeFeatureMaps), 'arg data wrong type'
        assert isinstance(idx, int) and idx >= 0, 'arg idx wrong type'
    
    
        self.inst = shared_ptr[_KDTreeFeatureNode](new _KDTreeFeatureNode((data.inst.get()), (<size_t>idx)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: KDTreeFeatureNode ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, data: KDTreeFeatureMaps , idx: int ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], KDTreeFeatureNode)):
             self._init_0(*args)
        elif (len(args)==2) and (isinstance(args[0], KDTreeFeatureMaps)) and (isinstance(args[1], int) and args[1] >= 0):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __getitem__(self,  i ):
        """
        __getitem__(self, i: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef long _idx = (<size_t>i)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        cdef double _r = deref(self.inst.get())[(<size_t>i)]
        py_result = <double>_r
        return py_result
    
    def getIndex(self):
        """
        getIndex(self) -> int
        Returns index of corresponding feature in data_
        """
        cdef size_t _r = self.inst.get().getIndex()
        py_result = <size_t>_r
        return py_result 

cdef class LogConfigHandler:
    """
    Cython implementation of _LogConfigHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LogConfigHandler.html>`_
    """

    
    def parse(self, list setting ):
        """
        parse(self, setting: List[bytes] ) -> Param
        Translates the given list of parameter settings into a LogStream configuration
        
        Translates the given list of parameter settings into a LogStream configuration.
        Usually this list stems from a command line call.
        
        Each element in the stringlist should follow this naming convention
        
        <LOG_NAME> <ACTION> <PARAMETER>
        
        with
        - LOG_NAME: DEBUG,INFO,WARNING,ERROR,FATAL_ERROR
        - ACTION: add,remove,clear
        - PARAMETER: for 'add'/'remove' it is the stream name (cout, cerr or a filename), 'clear' does not require any further parameter
        
        Example:
        `DEBUG add debug.log`
        
        This function will **not** apply to settings to the log handlers. Use configure() for that.
        
        :param setting: StringList containing the configuration options
        :raises ParseError: In case of an invalid configuration.
        :return: Param object containing all settings, that can be applied using the LogConfigHandler.configure() method
        """
        assert isinstance(setting, list) and all(isinstance(li, bytes) for li in setting), 'arg setting wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in setting:
           v0.push_back(_String(<char *>item0))
        cdef _Param * _r = new _Param(self.inst.get().parse(deref(v0)))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        setting[:] = replace
        del v0
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def configure(self, Param param ):
        """
        configure(self, param: Param ) -> None
        Applies the given parameters (@p param) to the current configuration
        
        <LOG_NAME> <ACTION> <PARAMETER> <STREAMTYPE>
        
        LOG_NAME: DEBUG, INFO, WARNING, ERROR, FATAL_ERROR
        ACTION: add, remove, clear
        PARAMETER: for 'add'/'remove' it is the stream name ('cout', 'cerr' or a filename), 'clear' does not require any further parameter
        STREAMTYPE: FILE, STRING (for a StringStream, which you can grab by this name using getStream() )
        
        You cannot specify a file named "cout" or "cerr" even if you specify streamtype 'FILE' - the handler will mistake this for the
        internal streams, but you can use "./cout" to print to a file named cout.
        
        A classical configuration would contain a list of settings e.g.
        
        `DEBUG add debug.log FILE`
        `INFO remove cout FILE` (FILE will be ignored)
        `INFO add string_stream1 STRING`
        
        :raises ElementNotFound: If the LogStream (first argument) does not exist.
        :raises FileNotWritable: If a file (or stream) should be opened as log file (or stream) that is not accessible.
        :raises IllegalArgument: If a stream should be registered, that was already registered with a different type.
        """
        assert isinstance(param, Param), 'arg param wrong type'
    
        self.inst.get().configure((deref(param.inst.get())))
    
    def setLogLevel(self,  log_level ):
        """
        setLogLevel(self, log_level: Union[bytes, str, String] ) -> None
        Sets a minimum log_level by removing all streams from loggers lower than that level.
        Valid levels are from low to high: "DEBUG", "INFO", "WARNING", "ERROR", "FATAL_ERROR"
        """
        assert (isinstance(log_level, str) or isinstance(log_level, bytes) or isinstance(log_level, String)), 'arg log_level wrong type'
    
        self.inst.get().setLogLevel(deref((convString(log_level)).get()))
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_LogConfigHandler] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_LogConfigHandler](_getInstance_LogConfigHandler())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class ModificationDefinition:
    """
    Cython implementation of _ModificationDefinition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationDefinition.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).getModificationName().c_str() )

    
    def __copy__(self):
       cdef ModificationDefinition rv = ModificationDefinition.__new__(ModificationDefinition)
       rv.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ModificationDefinition rv = ModificationDefinition.__new__(ModificationDefinition)
       rv.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition())
    
    def _init_1(self, ModificationDefinition in_0 ):
        """
        _init_1(self, in_0: ModificationDefinition ) -> None
        """
        assert isinstance(in_0, ModificationDefinition), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition((deref(in_0.inst.get()))))
    
    def _init_2(self,  mod ):
        """
        _init_2(self, mod: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(mod, str) or isinstance(mod, bytes) or isinstance(mod, String)), 'arg mod wrong type'
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref((convString(mod)).get())))
    
    def _init_3(self,  mod , bool fixed ):
        """
        _init_3(self, mod: Union[bytes, str, String] , fixed: bool ) -> None
        """
        assert (isinstance(mod, str) or isinstance(mod, bytes) or isinstance(mod, String)), 'arg mod wrong type'
        assert isinstance(fixed, pybool_t), 'arg fixed wrong type'
    
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref((convString(mod)).get()), (<bool>fixed)))
    
    def _init_4(self,  mod , bool fixed ,  max_occur ):
        """
        _init_4(self, mod: Union[bytes, str, String] , fixed: bool , max_occur: int ) -> None
        """
        assert (isinstance(mod, str) or isinstance(mod, bytes) or isinstance(mod, String)), 'arg mod wrong type'
        assert isinstance(fixed, pybool_t), 'arg fixed wrong type'
        assert isinstance(max_occur, int), 'arg max_occur wrong type'
    
    
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref((convString(mod)).get()), (<bool>fixed), (<unsigned int>max_occur)))
    
    def _init_5(self, ResidueModification mod ):
        """
        _init_5(self, mod: ResidueModification ) -> None
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition((deref(mod.inst.get()))))
    
    def _init_6(self, ResidueModification mod , bool fixed ):
        """
        _init_6(self, mod: ResidueModification , fixed: bool ) -> None
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
        assert isinstance(fixed, pybool_t), 'arg fixed wrong type'
    
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition((deref(mod.inst.get())), (<bool>fixed)))
    
    def _init_7(self, ResidueModification mod , bool fixed ,  max_occur ):
        """
        _init_7(self, mod: ResidueModification , fixed: bool , max_occur: int ) -> None
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
        assert isinstance(fixed, pybool_t), 'arg fixed wrong type'
        assert isinstance(max_occur, int), 'arg max_occur wrong type'
    
    
    
        self.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition((deref(mod.inst.get())), (<bool>fixed), (<unsigned int>max_occur)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ModificationDefinition ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: Union[bytes, str, String] , fixed: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: Union[bytes, str, String] , fixed: bool , max_occur: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: ResidueModification ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: ResidueModification , fixed: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mod: ResidueModification , fixed: bool , max_occur: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ModificationDefinition)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)):
             self._init_3(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], int)):
             self._init_4(*args)
        elif (len(args)==1) and (isinstance(args[0], ResidueModification)):
             self._init_5(*args)
        elif (len(args)==2) and (isinstance(args[0], ResidueModification)) and (isinstance(args[1], pybool_t)):
             self._init_6(*args)
        elif (len(args)==3) and (isinstance(args[0], ResidueModification)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], int)):
             self._init_7(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setFixedModification(self, bool fixed ):
        """
        setFixedModification(self, fixed: bool ) -> None
        Sets whether this modification definition is fixed or variable (modification must occur vs. can occur)
        """
        assert isinstance(fixed, pybool_t), 'arg fixed wrong type'
    
        self.inst.get().setFixedModification((<bool>fixed))
    
    def isFixedModification(self):
        """
        isFixedModification(self) -> bool
        Returns if the modification if fixed true, else false
        """
        cdef bool _r = self.inst.get().isFixedModification()
        py_result = <bool>_r
        return py_result
    
    def setMaxOccurrences(self,  num ):
        """
        setMaxOccurrences(self, num: int ) -> None
        Sets the maximal number of occurrences per peptide (unbounded if 0)
        """
        assert isinstance(num, int), 'arg num wrong type'
    
        self.inst.get().setMaxOccurrences((<unsigned int>num))
    
    def getMaxOccurrences(self):
        """
        getMaxOccurrences(self) -> int
        Returns the maximal number of occurrences per peptide
        """
        cdef unsigned int _r = self.inst.get().getMaxOccurrences()
        py_result = <unsigned int>_r
        return py_result
    
    def getModificationName(self):
        """
        getModificationName(self) -> Union[bytes, str, String]
        Returns the name of the modification
        """
        cdef _String _r = self.inst.get().getModificationName()
        py_result = convOutputString(_r)
        return py_result
    
    def setModification(self,  modification ):
        """
        setModification(self, modification: Union[bytes, str, String] ) -> None
        Sets the modification, allowed are unique names provided by ModificationsDB
        """
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
        self.inst.get().setModification(deref((convString(modification)).get()))
    
    def getModification(self):
        """
        getModification(self) -> ResidueModification
        """
        cdef _ResidueModification * _r = new _ResidueModification(self.inst.get().getModification())
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ModificationDefinition):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ModificationDefinition other_casted = other
        cdef ModificationDefinition self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class MultiplexDeltaMassesGenerator:
    """
    Cython implementation of _MultiplexDeltaMassesGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMassesGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MultiplexDeltaMassesGenerator rv = MultiplexDeltaMassesGenerator.__new__(MultiplexDeltaMassesGenerator)
       rv.inst = shared_ptr[_MultiplexDeltaMassesGenerator](new _MultiplexDeltaMassesGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MultiplexDeltaMassesGenerator rv = MultiplexDeltaMassesGenerator.__new__(MultiplexDeltaMassesGenerator)
       rv.inst = shared_ptr[_MultiplexDeltaMassesGenerator](new _MultiplexDeltaMassesGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MultiplexDeltaMassesGenerator](new _MultiplexDeltaMassesGenerator())
    
    def _init_1(self, MultiplexDeltaMassesGenerator in_0 ):
        """
        _init_1(self, in_0: MultiplexDeltaMassesGenerator ) -> None
        """
        assert isinstance(in_0, MultiplexDeltaMassesGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MultiplexDeltaMassesGenerator](new _MultiplexDeltaMassesGenerator((deref(in_0.inst.get()))))
    
    def _init_2(self,  labels ,  missed_cleavages , dict label_mass_shift ):
        """
        _init_2(self, labels: Union[bytes, str, String] , missed_cleavages: int , label_mass_shift: Dict[Union[bytes, str, String], float] ) -> None
        """
        assert (isinstance(labels, str) or isinstance(labels, bytes) or isinstance(labels, String)), 'arg labels wrong type'
        assert isinstance(missed_cleavages, int), 'arg missed_cleavages wrong type'
        assert isinstance(label_mass_shift, dict) and all((isinstance(k, str) or isinstance(k, bytes) or isinstance(k, String)) for k in label_mass_shift.keys()) and all(isinstance(v, float) for v in label_mass_shift.values()), 'arg label_mass_shift wrong type'
    
    
        cdef libcpp_map[_String, double] * v2 = new libcpp_map[_String, double]()
        for key, value in label_mass_shift.items():
        
        
            deref(v2)[ deref(<_String *> (<String> key).inst.get()) ] = (<double>value)
        
        
        self.inst = shared_ptr[_MultiplexDeltaMassesGenerator](new _MultiplexDeltaMassesGenerator(deref((convString(labels)).get()), (<int>missed_cleavages), deref(v2)))
        del v2
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MultiplexDeltaMassesGenerator ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, labels: Union[bytes, str, String] , missed_cleavages: int , label_mass_shift: Dict[Union[bytes, str, String], float] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MultiplexDeltaMassesGenerator)):
             self._init_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int)) and (isinstance(args[2], dict) and all((isinstance(k, str) or isinstance(k, bytes) or isinstance(k, String)) for k in args[2].keys()) and all(isinstance(v, float) for v in args[2].values())):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def generateKnockoutDeltaMasses(self):
        """
        generateKnockoutDeltaMasses(self) -> None
        """
        self.inst.get().generateKnockoutDeltaMasses()
    
    def getDeltaMassesList(self):
        """
        getDeltaMassesList(self) -> List[MultiplexDeltaMasses]
        """
        _r = self.inst.get().getDeltaMassesList()
        py_result = []
        cdef libcpp_vector[_MultiplexDeltaMasses].iterator it__r = _r.begin()
        cdef MultiplexDeltaMasses item_py_result
        while it__r != _r.end():
           item_py_result = MultiplexDeltaMasses.__new__(MultiplexDeltaMasses)
           item_py_result.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getLabelShort(self,  label ):
        """
        getLabelShort(self, label: Union[bytes, str, String] ) -> Union[bytes, str, String]
        """
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
        cdef _String _r = self.inst.get().getLabelShort(deref((convString(label)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def getLabelLong(self,  label ):
        """
        getLabelLong(self, label: Union[bytes, str, String] ) -> Union[bytes, str, String]
        """
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
        cdef _String _r = self.inst.get().getLabelLong(deref((convString(label)).get()))
        py_result = convOutputString(_r)
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

cdef class MultiplexDeltaMassesGenerator_Label:
    """
    Cython implementation of _MultiplexDeltaMassesGenerator_Label

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMassesGenerator_Label.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property short_name:
        def __set__(self,  short_name):
        
            self.inst.get().short_name = deref((convString(short_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().short_name
            py_result = convOutputString(_r)
            return py_result
    
    property long_name:
        def __set__(self,  long_name):
        
            self.inst.get().long_name = deref((convString(long_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().long_name
            py_result = convOutputString(_r)
            return py_result
    
    property description:
        def __set__(self,  description):
        
            self.inst.get().description = deref((convString(description)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().description
            py_result = convOutputString(_r)
            return py_result
    
    property delta_mass:
        def __set__(self, double delta_mass):
        
            self.inst.get().delta_mass = (<double>delta_mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().delta_mass
            py_result = <double>_r
            return py_result
    
    def __init__(self,  sn ,  ln ,  d , double dm ):
        """
        __init__(self, sn: Union[bytes, str, String] , ln: Union[bytes, str, String] , d: Union[bytes, str, String] , dm: float ) -> None
        """
        assert (isinstance(sn, str) or isinstance(sn, bytes) or isinstance(sn, String)), 'arg sn wrong type'
        assert (isinstance(ln, str) or isinstance(ln, bytes) or isinstance(ln, String)), 'arg ln wrong type'
        assert (isinstance(d, str) or isinstance(d, bytes) or isinstance(d, String)), 'arg d wrong type'
        assert isinstance(dm, float), 'arg dm wrong type'
    
    
    
    
        self.inst = shared_ptr[_MultiplexDeltaMassesGenerator_Label](new _MultiplexDeltaMassesGenerator_Label(deref((convString(sn)).get()), deref((convString(ln)).get()), deref((convString(d)).get()), (<double>dm))) 

cdef class MzQCFile:
    """
    Cython implementation of _MzQCFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzQCFile.html>`_

    File adapter for mzQC files used to load and store mzQC files
    
    This class collects the data for the mzQC File
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MzQCFile](new _MzQCFile())
    
    def store(self,  input_file ,  output_file , MSExperiment exp ,  contact_name ,  contact_address ,  description ,  label , FeatureMap feature_map , list prot_ids , PeptideIdentificationList pep_ids ):
        """
        store(self, input_file: Union[bytes, str, String] , output_file: Union[bytes, str, String] , exp: MSExperiment , contact_name: Union[bytes, str, String] , contact_address: Union[bytes, str, String] , description: Union[bytes, str, String] , label: Union[bytes, str, String] , feature_map: FeatureMap , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> None
        Stores QC data in mzQC file with JSON format
        
        
        :param input_file: MzML input file name
        :param output_file: MzQC output file name
        :param exp: MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
        :param contact_name: Name of the person creating the mzQC file
        :param contact_address: Contact address (mail/e-mail or phone) of the person creating the mzQC file
        :param description: Description and comments about the mzQC file contents
        :param label: Qnique and informative label for the run
        :param feature_map: FeatureMap from feature file (featureXML)
        :param prot_ids: Protein identifications from ID file (idXML)
        :param pep_ids: Protein identifications from ID file (idXML)
        """
        assert (isinstance(input_file, str) or isinstance(input_file, bytes) or isinstance(input_file, String)), 'arg input_file wrong type'
        assert (isinstance(output_file, str) or isinstance(output_file, bytes) or isinstance(output_file, String)), 'arg output_file wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(contact_name, str) or isinstance(contact_name, bytes) or isinstance(contact_name, String)), 'arg contact_name wrong type'
        assert (isinstance(contact_address, str) or isinstance(contact_address, bytes) or isinstance(contact_address, String)), 'arg contact_address wrong type'
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
        assert isinstance(feature_map, FeatureMap), 'arg feature_map wrong type'
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
    
    
    
    
    
    
    
    
        cdef libcpp_vector[_ProteinIdentification] * v8 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item8
        for item8 in prot_ids:
            v8.push_back(deref(item8.inst.get()))
    
        self.inst.get().store(deref((convString(input_file)).get()), deref((convString(output_file)).get()), (deref(exp.inst.get())), deref((convString(contact_name)).get()), deref((convString(contact_address)).get()), deref((convString(description)).get()), deref((convString(label)).get()), (deref(feature_map.inst.get())), deref(v8), (deref(pep_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v8.begin()
        replace_0 = []
        while it_prot_ids != v8.end():
            item8 = ProteinIdentification.__new__(ProteinIdentification)
            item8.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item8)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v8 

cdef class PeptideEvidence:
    """
    Cython implementation of _PeptideEvidence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideEvidence.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideEvidence rv = PeptideEvidence.__new__(PeptideEvidence)
       rv.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideEvidence rv = PeptideEvidence.__new__(PeptideEvidence)
       rv.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence())
    
    def _init_1(self, PeptideEvidence in_0 ):
        """
        _init_1(self, in_0: PeptideEvidence ) -> None
        """
        assert isinstance(in_0, PeptideEvidence), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideEvidence ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideEvidence)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setStart(self,  start ):
        """
        setStart(self, start: int ) -> None
        Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with `N_TERMINAL_AA`
        """
        assert isinstance(start, int), 'arg start wrong type'
    
        self.inst.get().setStart((<int>start))
    
    def getStart(self):
        """
        getStart(self) -> int
        Returns the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        """
        cdef int _r = self.inst.get().getStart()
        py_result = <int>_r
        return py_result
    
    def setEnd(self,  end ):
        """
        setEnd(self, end: int ) -> None
        Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA
        """
        assert isinstance(end, int), 'arg end wrong type'
    
        self.inst.get().setEnd((<int>end))
    
    def getEnd(self):
        """
        getEnd(self) -> int
        Returns the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned
        """
        cdef int _r = self.inst.get().getEnd()
        py_result = <int>_r
        return py_result
    
    def setAABefore(self, bytes rhs ):
        """
        setAABefore(self, rhs: bytes ) -> None
        Sets the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA
        """
        assert isinstance(rhs, bytes) and len(rhs) == 1, 'arg rhs wrong type'
    
        self.inst.get().setAABefore((<char>((rhs)[0])))
    
    def getAABefore(self):
        """
        getAABefore(self) -> bytes
        Returns the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned
        """
        cdef char  _r = self.inst.get().getAABefore()
        py_result = chr(<char>(_r))
        return py_result
    
    def setAAAfter(self, bytes rhs ):
        """
        setAAAfter(self, rhs: bytes ) -> None
        Sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA
        """
        assert isinstance(rhs, bytes) and len(rhs) == 1, 'arg rhs wrong type'
    
        self.inst.get().setAAAfter((<char>((rhs)[0])))
    
    def getAAAfter(self):
        """
        getAAAfter(self) -> bytes
        Returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned
        """
        cdef char  _r = self.inst.get().getAAAfter()
        py_result = chr(<char>(_r))
        return py_result
    
    def setProteinAccession(self,  s ):
        """
        setProteinAccession(self, s: Union[bytes, str, String] ) -> None
        Sets the protein accession the peptide matches to. If not available set to empty string
        """
        assert (isinstance(s, str) or isinstance(s, bytes) or isinstance(s, String)), 'arg s wrong type'
    
        self.inst.get().setProteinAccession(deref((convString(s)).get()))
    
    def getProteinAccession(self):
        """
        getProteinAccession(self) -> Union[bytes, str, String]
        Returns the protein accession the peptide matches to. If not available the empty string is returned
        """
        cdef _String _r = self.inst.get().getProteinAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def hasValidLimits(self):
        """
        hasValidLimits(self) -> bool
        Start and end numbers in evidence represent actual numeric indices
        """
        cdef bool _r = self.inst.get().hasValidLimits()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, PeptideEvidence):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef PeptideEvidence other_casted = other
        cdef PeptideEvidence self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class PrecursorCorrection:
    """
    Cython implementation of _PrecursorCorrection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PrecursorCorrection.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PrecursorCorrection rv = PrecursorCorrection.__new__(PrecursorCorrection)
       rv.inst = shared_ptr[_PrecursorCorrection](new _PrecursorCorrection(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PrecursorCorrection rv = PrecursorCorrection.__new__(PrecursorCorrection)
       rv.inst = shared_ptr[_PrecursorCorrection](new _PrecursorCorrection(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PrecursorCorrection](new _PrecursorCorrection())
    
    def _init_1(self, PrecursorCorrection in_0 ):
        """
        _init_1(self, in_0: PrecursorCorrection ) -> None
        """
        assert isinstance(in_0, PrecursorCorrection), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PrecursorCorrection](new _PrecursorCorrection((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PrecursorCorrection ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PrecursorCorrection)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    correctToHighestIntensityMS1Peak = __static_PrecursorCorrection_correctToHighestIntensityMS1Peak
    correctToNearestFeature = __static_PrecursorCorrection_correctToNearestFeature
    correctToNearestMS1Peak = __static_PrecursorCorrection_correctToNearestMS1Peak
    getPrecursors = __static_PrecursorCorrection_getPrecursors
    writeHist = __static_PrecursorCorrection_writeHist 

cdef class RibonucleotideDB:
    """
    Cython implementation of _RibonucleotideDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RibonucleotideDB.html>`_
    """

    
    def getRibonucleotide(self, bytes code ):
        """
        getRibonucleotide(self, code: bytes ) -> Ribonucleotide
        """
        assert isinstance(code, bytes), 'arg code wrong type'
    
        cdef const _Ribonucleotide * __r = (self.inst.get().getRibonucleotide((<libcpp_string>code)))
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    
    def getRibonucleotidePrefix(self, bytes code ):
        """
        getRibonucleotidePrefix(self, code: bytes ) -> Ribonucleotide
        """
        assert isinstance(code, bytes), 'arg code wrong type'
    
        cdef const _Ribonucleotide * __r = (self.inst.get().getRibonucleotidePrefix((<libcpp_string>code)))
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    


    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_RibonucleotideDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_RibonucleotideDB](_getInstance_RibonucleotideDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass


    def getRibonucleotideAlternatives(self, bytes code ):
        """Cython signature: libcpp_pair[const Ribonucleotide *,const Ribonucleotide *] getRibonucleotideAlternatives(const libcpp_string & code)"""
        assert isinstance(code, bytes), 'arg code wrong type'
    
        _r = self.inst.get().getRibonucleotideAlternatives((<libcpp_string>code))
        cdef const _Ribonucleotide * out_ptr1 = _r.first
        cdef const _Ribonucleotide * out_ptr2 = _r.second

        cdef Ribonucleotide out1 = Ribonucleotide.__new__(Ribonucleotide)
        cdef Ribonucleotide out2 = Ribonucleotide.__new__(Ribonucleotide)

        out1.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(_r.first)))
        out2.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(_r.second)))

        cdef list py_result = [out1, out2]
        return py_result 
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
