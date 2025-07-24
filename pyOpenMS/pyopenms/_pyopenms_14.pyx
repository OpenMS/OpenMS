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

cdef class __PeakMassType:
    None
    MONOISOTOPIC = 0
    AVERAGE = 1
    SIZE_OF_PEAKMASSTYPE = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __SpectrumType:
    None
    UNKNOWN = 0
    CENTROID = 1
    PROFILE = 2
    SIZE_OF_SPECTRUMTYPE = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class BiGaussFitter1D:
    """
    Cython implementation of _BiGaussFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BiGaussFitter1D.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef BiGaussFitter1D rv = BiGaussFitter1D.__new__(BiGaussFitter1D)
       rv.inst = shared_ptr[_BiGaussFitter1D](new _BiGaussFitter1D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef BiGaussFitter1D rv = BiGaussFitter1D.__new__(BiGaussFitter1D)
       rv.inst = shared_ptr[_BiGaussFitter1D](new _BiGaussFitter1D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BiGaussFitter1D](new _BiGaussFitter1D())
    
    def _init_1(self, BiGaussFitter1D in_0 ):
        """
        _init_1(self, in_0: BiGaussFitter1D ) -> None
        """
        assert isinstance(in_0, BiGaussFitter1D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BiGaussFitter1D](new _BiGaussFitter1D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: BiGaussFitter1D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], BiGaussFitter1D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class EmgScoring:
    """
    Cython implementation of _EmgScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgScoring.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EmgScoring rv = EmgScoring.__new__(EmgScoring)
       rv.inst = shared_ptr[_EmgScoring](new _EmgScoring(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EmgScoring rv = EmgScoring.__new__(EmgScoring)
       rv.inst = shared_ptr[_EmgScoring](new _EmgScoring(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Helps in scoring of an elution peak using an exponentially modified gaussian distribution model
        """
        self.inst = shared_ptr[_EmgScoring](new _EmgScoring())
    
    def _init_1(self, EmgScoring in_0 ):
        """
        _init_1(self, in_0: EmgScoring ) -> None
        """
        assert isinstance(in_0, EmgScoring), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmgScoring](new _EmgScoring((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Helps in scoring of an elution peak using an exponentially modified gaussian distribution model

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EmgScoring ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EmgScoring)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setFitterParam(self, Param param ):
        """
        setFitterParam(self, param: Param ) -> None
        """
        assert isinstance(param, Param), 'arg param wrong type'
    
        self.inst.get().setFitterParam((deref(param.inst.get())))
    
    def getDefaults(self):
        """
        getDefaults(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getDefaults())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def elutionModelFit(self, np.ndarray[np.float32_t,ndim=2] current_section , bool smooth_data ):
        """
        elutionModelFit(self, current_section: '_np.ndarray[Any, _np.dtype[_np.float32]]' , smooth_data: bool ) -> float
        """
        assert current_section.shape[1] == 2, 'arg current_section wrong type'
        assert isinstance(smooth_data, pybool_t), 'arg smooth_data wrong type'
        cdef libcpp_vector[_DPosition2] _dp_vec_0
        cdef _DPosition2 _dp_0
        cdef int _dp_ii_0
        cdef int _dp_N_0 = current_section.shape[0]
        for _dp_ii_0 in range(_dp_N_0):
            _dp_0[0] = current_section[_dp_ii_0,0]
            _dp_0[1] = current_section[_dp_ii_0,1]
            _dp_vec_0.push_back(_dp_0)
    
        cdef double _r = self.inst.get().elutionModelFit(_dp_vec_0, (<bool>smooth_data))
        py_result = <double>_r
        return py_result 

cdef class GaussFitResult:
    """
    Cython implementation of _GaussFitResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1GaussFitResult.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property A:
        def __set__(self, double A):
        
            self.inst.get().A = (<double>A)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().A
            py_result = <double>_r
            return py_result
    
    property x0:
        def __set__(self, double x0):
        
            self.inst.get().x0 = (<double>x0)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().x0
            py_result = <double>_r
            return py_result
    
    property sigma:
        def __set__(self, double sigma):
        
            self.inst.get().sigma = (<double>sigma)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().sigma
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef GaussFitResult rv = GaussFitResult.__new__(GaussFitResult)
       rv.inst = shared_ptr[_GaussFitResult](new _GaussFitResult(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GaussFitResult rv = GaussFitResult.__new__(GaussFitResult)
       rv.inst = shared_ptr[_GaussFitResult](new _GaussFitResult(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_GaussFitResult](new _GaussFitResult())
    
    def _init_1(self, double in_0 , double in_1 , double in_2 ):
        """
        _init_1(self, in_0: float , in_1: float , in_2: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
        assert isinstance(in_2, float), 'arg in_2 wrong type'
    
    
    
        self.inst = shared_ptr[_GaussFitResult](new _GaussFitResult((<double>in_0), (<double>in_1), (<double>in_2)))
    
    def _init_2(self, GaussFitResult in_0 ):
        """
        _init_2(self, in_0: GaussFitResult ) -> None
        """
        assert isinstance(in_0, GaussFitResult), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GaussFitResult](new _GaussFitResult((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: float , in_1: float , in_2: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GaussFitResult ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], GaussFitResult)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def eval(self, double in_0 ):
        """
        eval(self, in_0: float ) -> float
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().eval((<double>in_0))
        py_result = <double>_r
        return py_result 

cdef class GaussFitter:
    """
    Cython implementation of _GaussFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1GaussFitter.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        Implements a fitter for Gaussian functions
        """
        self.inst = shared_ptr[_GaussFitter](new _GaussFitter())
    
    def setInitialParameters(self, GaussFitResult result ):
        """
        setInitialParameters(self, result: GaussFitResult ) -> None
        Sets the initial parameters used by the fit method as initial guess for the Gaussian
        """
        assert isinstance(result, GaussFitResult), 'arg result wrong type'
    
        self.inst.get().setInitialParameters((deref(result.inst.get())))
    
    def fit(self, np.ndarray[np.float32_t,ndim=2] points ):
        """
        fit(self, points: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> GaussFitResult
        Fits a Gaussian distribution to the given data points
        """
        assert points.shape[1] == 2, 'arg points wrong type'
        cdef libcpp_vector[_DPosition2] _dp_vec_0
        cdef _DPosition2 _dp_0
        cdef int _dp_ii_0
        cdef int _dp_N_0 = points.shape[0]
        for _dp_ii_0 in range(_dp_N_0):
            _dp_0[0] = points[_dp_ii_0,0]
            _dp_0[1] = points[_dp_ii_0,1]
            _dp_vec_0.push_back(_dp_0)
        cdef _GaussFitResult * _r = new _GaussFitResult(self.inst.get().fit(_dp_vec_0))
        cdef GaussFitResult py_result = GaussFitResult.__new__(GaussFitResult)
        py_result.inst = shared_ptr[_GaussFitResult](_r)
        return py_result 

cdef class LowessSmoothing:
    """
    Cython implementation of _LowessSmoothing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LowessSmoothing.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_LowessSmoothing](new _LowessSmoothing())
    
    def smoothData(self, list x , list y , list y_smoothed ):
        """
        smoothData(self, x: List[float] , y: List[float] , y_smoothed: List[float] ) -> None
        Smoothing method that receives x and y coordinates (e.g., RT and intensities) and computes smoothed intensities
        """
        assert isinstance(x, list) and all(isinstance(elemt_rec, float) for elemt_rec in x), 'arg x wrong type'
        assert isinstance(y, list) and all(isinstance(elemt_rec, float) for elemt_rec in y), 'arg y wrong type'
        assert isinstance(y_smoothed, list) and all(isinstance(elemt_rec, float) for elemt_rec in y_smoothed), 'arg y_smoothed wrong type'
        cdef libcpp_vector[double] v0 = x
        cdef libcpp_vector[double] v1 = y
        cdef libcpp_vector[double] v2 = y_smoothed
        self.inst.get().smoothData(v0, v1, v2)
        y_smoothed[:] = v2
        
        
    
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

cdef class MRMFQC_ComponentGroupPairQCs:
    """
    Cython implementation of _MRMFQC_ComponentGroupPairQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentGroupPairQCs.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property component_group_name:
        def __set__(self,  component_group_name):
        
            self.inst.get().component_group_name = deref((convString(component_group_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().component_group_name
            py_result = convOutputString(_r)
            return py_result
    
    property resolution_pair_name:
        def __set__(self,  resolution_pair_name):
        
            self.inst.get().resolution_pair_name = deref((convString(resolution_pair_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().resolution_pair_name
            py_result = convOutputString(_r)
            return py_result
    
    property resolution_l:
        def __set__(self, double resolution_l):
        
            self.inst.get().resolution_l = (<double>resolution_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().resolution_l
            py_result = <double>_r
            return py_result
    
    property resolution_u:
        def __set__(self, double resolution_u):
        
            self.inst.get().resolution_u = (<double>resolution_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().resolution_u
            py_result = <double>_r
            return py_result
    
    property rt_diff_l:
        def __set__(self, double rt_diff_l):
        
            self.inst.get().rt_diff_l = (<double>rt_diff_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_diff_l
            py_result = <double>_r
            return py_result
    
    property rt_diff_u:
        def __set__(self, double rt_diff_u):
        
            self.inst.get().rt_diff_u = (<double>rt_diff_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_diff_u
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef MRMFQC_ComponentGroupPairQCs rv = MRMFQC_ComponentGroupPairQCs.__new__(MRMFQC_ComponentGroupPairQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentGroupPairQCs](new _MRMFQC_ComponentGroupPairQCs(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFQC_ComponentGroupPairQCs rv = MRMFQC_ComponentGroupPairQCs.__new__(MRMFQC_ComponentGroupPairQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentGroupPairQCs](new _MRMFQC_ComponentGroupPairQCs(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFQC_ComponentGroupPairQCs](new _MRMFQC_ComponentGroupPairQCs())
    
    def _init_1(self, MRMFQC_ComponentGroupPairQCs in_0 ):
        """
        _init_1(self, in_0: MRMFQC_ComponentGroupPairQCs ) -> None
        """
        assert isinstance(in_0, MRMFQC_ComponentGroupPairQCs), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFQC_ComponentGroupPairQCs](new _MRMFQC_ComponentGroupPairQCs((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFQC_ComponentGroupPairQCs ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFQC_ComponentGroupPairQCs)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMFQC_ComponentGroupQCs:
    """
    Cython implementation of _MRMFQC_ComponentGroupQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentGroupQCs.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property component_group_name:
        def __set__(self,  component_group_name):
        
            self.inst.get().component_group_name = deref((convString(component_group_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().component_group_name
            py_result = convOutputString(_r)
            return py_result
    
    property retention_time_l:
        def __set__(self, double retention_time_l):
        
            self.inst.get().retention_time_l = (<double>retention_time_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().retention_time_l
            py_result = <double>_r
            return py_result
    
    property retention_time_u:
        def __set__(self, double retention_time_u):
        
            self.inst.get().retention_time_u = (<double>retention_time_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().retention_time_u
            py_result = <double>_r
            return py_result
    
    property intensity_l:
        def __set__(self, double intensity_l):
        
            self.inst.get().intensity_l = (<double>intensity_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity_l
            py_result = <double>_r
            return py_result
    
    property intensity_u:
        def __set__(self, double intensity_u):
        
            self.inst.get().intensity_u = (<double>intensity_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity_u
            py_result = <double>_r
            return py_result
    
    property overall_quality_l:
        def __set__(self, double overall_quality_l):
        
            self.inst.get().overall_quality_l = (<double>overall_quality_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().overall_quality_l
            py_result = <double>_r
            return py_result
    
    property overall_quality_u:
        def __set__(self, double overall_quality_u):
        
            self.inst.get().overall_quality_u = (<double>overall_quality_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().overall_quality_u
            py_result = <double>_r
            return py_result
    
    property n_heavy_l:
        def __set__(self,  n_heavy_l):
        
            self.inst.get().n_heavy_l = (<int>n_heavy_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_heavy_l
            py_result = <int>_r
            return py_result
    
    property n_heavy_u:
        def __set__(self,  n_heavy_u):
        
            self.inst.get().n_heavy_u = (<int>n_heavy_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_heavy_u
            py_result = <int>_r
            return py_result
    
    property n_light_l:
        def __set__(self,  n_light_l):
        
            self.inst.get().n_light_l = (<int>n_light_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_light_l
            py_result = <int>_r
            return py_result
    
    property n_light_u:
        def __set__(self,  n_light_u):
        
            self.inst.get().n_light_u = (<int>n_light_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_light_u
            py_result = <int>_r
            return py_result
    
    property n_detecting_l:
        def __set__(self,  n_detecting_l):
        
            self.inst.get().n_detecting_l = (<int>n_detecting_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_detecting_l
            py_result = <int>_r
            return py_result
    
    property n_detecting_u:
        def __set__(self,  n_detecting_u):
        
            self.inst.get().n_detecting_u = (<int>n_detecting_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_detecting_u
            py_result = <int>_r
            return py_result
    
    property n_quantifying_l:
        def __set__(self,  n_quantifying_l):
        
            self.inst.get().n_quantifying_l = (<int>n_quantifying_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_quantifying_l
            py_result = <int>_r
            return py_result
    
    property n_quantifying_u:
        def __set__(self,  n_quantifying_u):
        
            self.inst.get().n_quantifying_u = (<int>n_quantifying_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_quantifying_u
            py_result = <int>_r
            return py_result
    
    property n_identifying_l:
        def __set__(self,  n_identifying_l):
        
            self.inst.get().n_identifying_l = (<int>n_identifying_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_identifying_l
            py_result = <int>_r
            return py_result
    
    property n_identifying_u:
        def __set__(self,  n_identifying_u):
        
            self.inst.get().n_identifying_u = (<int>n_identifying_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_identifying_u
            py_result = <int>_r
            return py_result
    
    property n_transitions_l:
        def __set__(self,  n_transitions_l):
        
            self.inst.get().n_transitions_l = (<int>n_transitions_l)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_transitions_l
            py_result = <int>_r
            return py_result
    
    property n_transitions_u:
        def __set__(self,  n_transitions_u):
        
            self.inst.get().n_transitions_u = (<int>n_transitions_u)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().n_transitions_u
            py_result = <int>_r
            return py_result
    
    property ion_ratio_pair_name_1:
        def __set__(self,  ion_ratio_pair_name_1):
        
            self.inst.get().ion_ratio_pair_name_1 = deref((convString(ion_ratio_pair_name_1)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().ion_ratio_pair_name_1
            py_result = convOutputString(_r)
            return py_result
    
    property ion_ratio_pair_name_2:
        def __set__(self,  ion_ratio_pair_name_2):
        
            self.inst.get().ion_ratio_pair_name_2 = deref((convString(ion_ratio_pair_name_2)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().ion_ratio_pair_name_2
            py_result = convOutputString(_r)
            return py_result
    
    property ion_ratio_l:
        def __set__(self, double ion_ratio_l):
        
            self.inst.get().ion_ratio_l = (<double>ion_ratio_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ion_ratio_l
            py_result = <double>_r
            return py_result
    
    property ion_ratio_u:
        def __set__(self, double ion_ratio_u):
        
            self.inst.get().ion_ratio_u = (<double>ion_ratio_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ion_ratio_u
            py_result = <double>_r
            return py_result
    
    property ion_ratio_feature_name:
        def __set__(self,  ion_ratio_feature_name):
        
            self.inst.get().ion_ratio_feature_name = deref((convString(ion_ratio_feature_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().ion_ratio_feature_name
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef MRMFQC_ComponentGroupQCs rv = MRMFQC_ComponentGroupQCs.__new__(MRMFQC_ComponentGroupQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentGroupQCs](new _MRMFQC_ComponentGroupQCs(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFQC_ComponentGroupQCs rv = MRMFQC_ComponentGroupQCs.__new__(MRMFQC_ComponentGroupQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentGroupQCs](new _MRMFQC_ComponentGroupQCs(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFQC_ComponentGroupQCs](new _MRMFQC_ComponentGroupQCs())
    
    def _init_1(self, MRMFQC_ComponentGroupQCs in_0 ):
        """
        _init_1(self, in_0: MRMFQC_ComponentGroupQCs ) -> None
        """
        assert isinstance(in_0, MRMFQC_ComponentGroupQCs), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFQC_ComponentGroupQCs](new _MRMFQC_ComponentGroupQCs((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFQC_ComponentGroupQCs ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFQC_ComponentGroupQCs)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMFQC_ComponentQCs:
    """
    Cython implementation of _MRMFQC_ComponentQCs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFQC_ComponentQCs.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property component_name:
        def __set__(self,  component_name):
        
            self.inst.get().component_name = deref((convString(component_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().component_name
            py_result = convOutputString(_r)
            return py_result
    
    property retention_time_l:
        def __set__(self, double retention_time_l):
        
            self.inst.get().retention_time_l = (<double>retention_time_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().retention_time_l
            py_result = <double>_r
            return py_result
    
    property retention_time_u:
        def __set__(self, double retention_time_u):
        
            self.inst.get().retention_time_u = (<double>retention_time_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().retention_time_u
            py_result = <double>_r
            return py_result
    
    property intensity_l:
        def __set__(self, double intensity_l):
        
            self.inst.get().intensity_l = (<double>intensity_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity_l
            py_result = <double>_r
            return py_result
    
    property intensity_u:
        def __set__(self, double intensity_u):
        
            self.inst.get().intensity_u = (<double>intensity_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity_u
            py_result = <double>_r
            return py_result
    
    property overall_quality_l:
        def __set__(self, double overall_quality_l):
        
            self.inst.get().overall_quality_l = (<double>overall_quality_l)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().overall_quality_l
            py_result = <double>_r
            return py_result
    
    property overall_quality_u:
        def __set__(self, double overall_quality_u):
        
            self.inst.get().overall_quality_u = (<double>overall_quality_u)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().overall_quality_u
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef MRMFQC_ComponentQCs rv = MRMFQC_ComponentQCs.__new__(MRMFQC_ComponentQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentQCs](new _MRMFQC_ComponentQCs(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFQC_ComponentQCs rv = MRMFQC_ComponentQCs.__new__(MRMFQC_ComponentQCs)
       rv.inst = shared_ptr[_MRMFQC_ComponentQCs](new _MRMFQC_ComponentQCs(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFQC_ComponentQCs](new _MRMFQC_ComponentQCs())
    
    def _init_1(self, MRMFQC_ComponentQCs in_0 ):
        """
        _init_1(self, in_0: MRMFQC_ComponentQCs ) -> None
        """
        assert isinstance(in_0, MRMFQC_ComponentQCs), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFQC_ComponentQCs](new _MRMFQC_ComponentQCs((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFQC_ComponentQCs ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFQC_ComponentQCs)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMFeatureQC:
    """
    Cython implementation of _MRMFeatureQC

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureQC.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property component_qcs:
        def __set__(self, list component_qcs):
            cdef libcpp_vector[_MRMFQC_ComponentQCs] * v0 = new libcpp_vector[_MRMFQC_ComponentQCs]()
            cdef MRMFQC_ComponentQCs item0
            for item0 in component_qcs:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().component_qcs = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().component_qcs
            py_result = []
            cdef libcpp_vector[_MRMFQC_ComponentQCs].iterator it__r = _r.begin()
            cdef MRMFQC_ComponentQCs item_py_result
            while it__r != _r.end():
               item_py_result = MRMFQC_ComponentQCs.__new__(MRMFQC_ComponentQCs)
               item_py_result.inst = shared_ptr[_MRMFQC_ComponentQCs](new _MRMFQC_ComponentQCs(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property component_group_qcs:
        def __set__(self, list component_group_qcs):
            cdef libcpp_vector[_MRMFQC_ComponentGroupQCs] * v0 = new libcpp_vector[_MRMFQC_ComponentGroupQCs]()
            cdef MRMFQC_ComponentGroupQCs item0
            for item0 in component_group_qcs:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().component_group_qcs = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().component_group_qcs
            py_result = []
            cdef libcpp_vector[_MRMFQC_ComponentGroupQCs].iterator it__r = _r.begin()
            cdef MRMFQC_ComponentGroupQCs item_py_result
            while it__r != _r.end():
               item_py_result = MRMFQC_ComponentGroupQCs.__new__(MRMFQC_ComponentGroupQCs)
               item_py_result.inst = shared_ptr[_MRMFQC_ComponentGroupQCs](new _MRMFQC_ComponentGroupQCs(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property component_group_pair_qcs:
        def __set__(self, list component_group_pair_qcs):
            cdef libcpp_vector[_MRMFQC_ComponentGroupPairQCs] * v0 = new libcpp_vector[_MRMFQC_ComponentGroupPairQCs]()
            cdef MRMFQC_ComponentGroupPairQCs item0
            for item0 in component_group_pair_qcs:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().component_group_pair_qcs = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().component_group_pair_qcs
            py_result = []
            cdef libcpp_vector[_MRMFQC_ComponentGroupPairQCs].iterator it__r = _r.begin()
            cdef MRMFQC_ComponentGroupPairQCs item_py_result
            while it__r != _r.end():
               item_py_result = MRMFQC_ComponentGroupPairQCs.__new__(MRMFQC_ComponentGroupPairQCs)
               item_py_result.inst = shared_ptr[_MRMFQC_ComponentGroupPairQCs](new _MRMFQC_ComponentGroupPairQCs(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef MRMFeatureQC rv = MRMFeatureQC.__new__(MRMFeatureQC)
       rv.inst = shared_ptr[_MRMFeatureQC](new _MRMFeatureQC(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeatureQC rv = MRMFeatureQC.__new__(MRMFeatureQC)
       rv.inst = shared_ptr[_MRMFeatureQC](new _MRMFeatureQC(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeatureQC](new _MRMFeatureQC())
    
    def _init_1(self, MRMFeatureQC in_0 ):
        """
        _init_1(self, in_0: MRMFeatureQC ) -> None
        """
        assert isinstance(in_0, MRMFeatureQC), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeatureQC](new _MRMFeatureQC((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeatureQC ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeatureQC)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMMapping:
    """
    Cython implementation of _MRMMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMMapping.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MRMMapping](new _MRMMapping())
    
    def mapExperiment(self, MSExperiment input_chromatograms , TargetedExperiment targeted_exp , MSExperiment output ):
        """
        mapExperiment(self, input_chromatograms: MSExperiment , targeted_exp: TargetedExperiment , output: MSExperiment ) -> None
        Maps input chromatograms to assays in a targeted experiment
        
        The output chromatograms are an annotated copy of the input chromatograms
        with native id, precursor information and peptide sequence (if available)
        annotated in the chromatogram files
        
        The algorithm tries to match a given set of chromatograms and targeted
        assays. It iterates through all the chromatograms retrieves one or more
        matching targeted assay for the chromatogram. By default, the algorithm
        assumes that a 1:1 mapping exists. If a chromatogram cannot be mapped
        (does not have a corresponding assay) the algorithm issues a warning, the
        user can specify that the program should abort in such a case (see
        error_on_unmapped)
        
        :note If multiple mapping is enabled (see map_multiple_assays parameter)
        then each mapped assay will get its own chromatogram that contains the
        same raw data but different meta-annotation. This *can* be useful if the
        same transition is used to monitor multiple analytes but may also
        indicate a problem with too wide mapping tolerances
        """
        assert isinstance(input_chromatograms, MSExperiment), 'arg input_chromatograms wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
        assert isinstance(output, MSExperiment), 'arg output wrong type'
    
    
    
        self.inst.get().mapExperiment((deref(input_chromatograms.inst.get())), (deref(targeted_exp.inst.get())), (deref(output.inst.get())))
    
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

cdef class MapAlignmentAlgorithmPoseClustering:
    """
    Cython implementation of _MapAlignmentAlgorithmPoseClustering

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmPoseClustering.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MapAlignmentAlgorithmPoseClustering](new _MapAlignmentAlgorithmPoseClustering())
    
    def _align_0(self, FeatureMap in_0 , TransformationDescription in_1 ):
        """
        _align_0(self, in_0: FeatureMap , in_1: TransformationDescription ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
    
    
        self.inst.get().align((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def _align_1(self, MSExperiment in_0 , TransformationDescription in_1 ):
        """
        _align_1(self, in_0: MSExperiment , in_1: TransformationDescription ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
    
    
        self.inst.get().align((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def align(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: align(self, in_0: FeatureMap , in_1: TransformationDescription ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: align(self, in_0: MSExperiment , in_1: TransformationDescription ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], TransformationDescription)):
            return self._align_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], TransformationDescription)):
            return self._align_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setReference_0(self, FeatureMap in_0 ):
        """
        _setReference_0(self, in_0: FeatureMap ) -> None
        Sets the reference for the alignment
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
    
        self.inst.get().setReference((deref(in_0.inst.get())))
    
    def _setReference_1(self, MSExperiment in_0 ):
        """
        _setReference_1(self, in_0: MSExperiment ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
    
        self.inst.get().setReference((deref(in_0.inst.get())))
    
    def setReference(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setReference(self, in_0: FeatureMap ) -> None
          :noindex:
        
        Sets the reference for the alignment

        
        .. rubric:: Overload:
        .. py:function:: setReference(self, in_0: MSExperiment ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], FeatureMap)):
            return self._setReference_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSExperiment)):
            return self._setReference_1(*args)
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

cdef class MassExplainer:
    """
    Cython implementation of _MassExplainer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassExplainer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MassExplainer rv = MassExplainer.__new__(MassExplainer)
       rv.inst = shared_ptr[_MassExplainer](new _MassExplainer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassExplainer rv = MassExplainer.__new__(MassExplainer)
       rv.inst = shared_ptr[_MassExplainer](new _MassExplainer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Computes empirical formulas for given mass differences using a set of allowed elements
        """
        self.inst = shared_ptr[_MassExplainer](new _MassExplainer())
    
    def _init_1(self, MassExplainer in_0 ):
        """
        _init_1(self, in_0: MassExplainer ) -> None
        """
        assert isinstance(in_0, MassExplainer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassExplainer](new _MassExplainer((deref(in_0.inst.get()))))
    
    def _init_2(self, list adduct_base ):
        """
        _init_2(self, adduct_base: List[Adduct] ) -> None
        """
        assert isinstance(adduct_base, list) and all(isinstance(elemt_rec, Adduct) for elemt_rec in adduct_base), 'arg adduct_base wrong type'
        cdef libcpp_vector[_Adduct] * v0 = new libcpp_vector[_Adduct]()
        cdef Adduct item0
        for item0 in adduct_base:
            v0.push_back(deref(item0.inst.get()))
        self.inst = shared_ptr[_MassExplainer](new _MassExplainer(deref(v0)))
        del v0
    
    def _init_3(self,  q_min ,  q_max ,  max_span , double thresh_logp ):
        """
        _init_3(self, q_min: int , q_max: int , max_span: int , thresh_logp: float ) -> None
        """
        assert isinstance(q_min, int), 'arg q_min wrong type'
        assert isinstance(q_max, int), 'arg q_max wrong type'
        assert isinstance(max_span, int), 'arg max_span wrong type'
        assert isinstance(thresh_logp, float), 'arg thresh_logp wrong type'
    
    
    
    
        self.inst = shared_ptr[_MassExplainer](new _MassExplainer((<int>q_min), (<int>q_max), (<int>max_span), (<double>thresh_logp)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Computes empirical formulas for given mass differences using a set of allowed elements

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassExplainer ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, adduct_base: List[Adduct] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, q_min: int , q_max: int , max_span: int , thresh_logp: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassExplainer)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, Adduct) for elemt_rec in args[0])):
             self._init_2(*args)
        elif (len(args)==4) and (isinstance(args[0], int)) and (isinstance(args[1], int)) and (isinstance(args[2], int)) and (isinstance(args[3], float)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setAdductBase(self, list adduct_base ):
        """
        setAdductBase(self, adduct_base: List[Adduct] ) -> None
        Sets the set of possible adducts
        """
        assert isinstance(adduct_base, list) and all(isinstance(elemt_rec, Adduct) for elemt_rec in adduct_base), 'arg adduct_base wrong type'
        cdef libcpp_vector[_Adduct] * v0 = new libcpp_vector[_Adduct]()
        cdef Adduct item0
        for item0 in adduct_base:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setAdductBase(deref(v0))
        del v0
    
    def getAdductBase(self):
        """
        getAdductBase(self) -> List[Adduct]
        Returns the set of adducts
        """
        _r = self.inst.get().getAdductBase()
        py_result = []
        cdef libcpp_vector[_Adduct].iterator it__r = _r.begin()
        cdef Adduct item_py_result
        while it__r != _r.end():
           item_py_result = Adduct.__new__(Adduct)
           item_py_result.inst = shared_ptr[_Adduct](new _Adduct(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getCompomerById(self,  id ):
        """
        getCompomerById(self, id: int ) -> Compomer
        Returns a compomer by its Id (useful after a query() )
        """
        assert isinstance(id, int) and id >= 0, 'arg id wrong type'
    
        cdef _Compomer * _r = new _Compomer(self.inst.get().getCompomerById((<size_t>id)))
        cdef Compomer py_result = Compomer.__new__(Compomer)
        py_result.inst = shared_ptr[_Compomer](_r)
        return py_result
    
    def compute(self):
        """
        compute(self) -> None
        Fill map with possible mass-differences along with their explanation
        """
        self.inst.get().compute() 

cdef class PrecursorPurity:
    """
    Cython implementation of _PrecursorPurity

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PrecursorPurity.html>`_

    Precursor purity or noise estimation
    
    This class computes metrics for precursor isolation window purity (or noise)
    The function extracts the peaks from an isolation window targeted for fragmentation
    and determines which peaks are isotopes of the target and which come from other sources
    The intensities of the assumed target peaks are summed up as the target intensity
    Using this information it calculates an intensity ratio for the relative intensity of the target
    compared to other sources
    These metrics are combined over the previous and the next MS1 spectrum
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PrecursorPurity rv = PrecursorPurity.__new__(PrecursorPurity)
       rv.inst = shared_ptr[_PrecursorPurity](new _PrecursorPurity(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PrecursorPurity rv = PrecursorPurity.__new__(PrecursorPurity)
       rv.inst = shared_ptr[_PrecursorPurity](new _PrecursorPurity(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PrecursorPurity](new _PrecursorPurity())
    
    def _init_1(self, PrecursorPurity in_0 ):
        """
        _init_1(self, in_0: PrecursorPurity ) -> None
        """
        assert isinstance(in_0, PrecursorPurity), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PrecursorPurity](new _PrecursorPurity((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PrecursorPurity ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PrecursorPurity)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def computePrecursorPurity(self, MSSpectrum ms1 , Precursor pre , double precursor_mass_tolerance , bool precursor_mass_tolerance_unit_ppm ):
        """
        computePrecursorPurity(self, ms1: MSSpectrum , pre: Precursor , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool ) -> PurityScores
        Compute precursor purity metrics for one MS2 precursor
        
        Note: This function is implemented in a general way and can also be used for e.g. MS3 precursor isolation windows in MS2 spectra
        Spectra annotated with charge 0 will be treated as charge 1.
        
        
        :param ms1: The Spectrum containing the isolation window
        :param pre: The precursor containing the definition the isolation window
        :param precursor_mass_tolerance: The precursor tolerance. Is used for determining the targeted peak and deisotoping
        :param precursor_mass_tolerance_unit_ppm: The unit of the precursor tolerance
        """
        assert isinstance(ms1, MSSpectrum), 'arg ms1 wrong type'
        assert isinstance(pre, Precursor), 'arg pre wrong type'
        assert isinstance(precursor_mass_tolerance, float), 'arg precursor_mass_tolerance wrong type'
        assert isinstance(precursor_mass_tolerance_unit_ppm, pybool_t), 'arg precursor_mass_tolerance_unit_ppm wrong type'
    
    
    
    
        cdef _PurityScores * _r = new _PurityScores(self.inst.get().computePrecursorPurity((deref(ms1.inst.get())), (deref(pre.inst.get())), (<double>precursor_mass_tolerance), (<bool>precursor_mass_tolerance_unit_ppm)))
        cdef PurityScores py_result = PurityScores.__new__(PurityScores)
        py_result.inst = shared_ptr[_PurityScores](_r)
        return py_result 

cdef class ProteinGroup:
    """
    Cython implementation of _ProteinGroup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinGroup.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property probability:
        def __set__(self, double probability):
        
            self.inst.get().probability = (<double>probability)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().probability
            py_result = <double>_r
            return py_result
    
    property accessions:
        def __set__(self, list accessions):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in accessions:
               v0.push_back(_String(<char *>item0))
            self.inst.get().accessions = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().accessions
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef ProteinGroup rv = ProteinGroup.__new__(ProteinGroup)
       rv.inst = shared_ptr[_ProteinGroup](new _ProteinGroup(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteinGroup rv = ProteinGroup.__new__(ProteinGroup)
       rv.inst = shared_ptr[_ProteinGroup](new _ProteinGroup(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteinGroup](new _ProteinGroup())
    
    def _init_1(self, ProteinGroup in_0 ):
        """
        _init_1(self, in_0: ProteinGroup ) -> None
        """
        assert isinstance(in_0, ProteinGroup), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteinGroup](new _ProteinGroup((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteinGroup ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteinGroup)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ProteinIdentification:
    """
    Cython implementation of _ProteinIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinIdentification.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ProteinIdentification rv = ProteinIdentification.__new__(ProteinIdentification)
       rv.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteinIdentification rv = ProteinIdentification.__new__(ProteinIdentification)
       rv.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification())
    
    def _init_1(self, ProteinIdentification in_0 ):
        """
        _init_1(self, in_0: ProteinIdentification ) -> None
        """
        assert isinstance(in_0, ProteinIdentification), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteinIdentification ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteinIdentification)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getHits(self):
        """
        getHits(self) -> List[ProteinHit]
        Returns the protein hits
        """
        _r = self.inst.get().getHits()
        py_result = []
        cdef libcpp_vector[_ProteinHit].iterator it__r = _r.begin()
        cdef ProteinHit item_py_result
        while it__r != _r.end():
           item_py_result = ProteinHit.__new__(ProteinHit)
           item_py_result.inst = shared_ptr[_ProteinHit](new _ProteinHit(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def insertHit(self, ProteinHit input ):
        """
        insertHit(self, input: ProteinHit ) -> None
        Appends a protein hit
        """
        assert isinstance(input, ProteinHit), 'arg input wrong type'
    
        self.inst.get().insertHit((deref(input.inst.get())))
    
    def setHits(self, list hits ):
        """
        setHits(self, hits: List[ProteinHit] ) -> None
        Sets the protein hits
        """
        assert isinstance(hits, list) and all(isinstance(elemt_rec, ProteinHit) for elemt_rec in hits), 'arg hits wrong type'
        cdef libcpp_vector[_ProteinHit] * v0 = new libcpp_vector[_ProteinHit]()
        cdef ProteinHit item0
        for item0 in hits:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setHits(deref(v0))
        del v0
    
    def getProteinGroups(self):
        """
        getProteinGroups(self) -> List[ProteinGroup]
        Returns the protein groups
        """
        _r = self.inst.get().getProteinGroups()
        py_result = []
        cdef libcpp_vector[_ProteinGroup].iterator it__r = _r.begin()
        cdef ProteinGroup item_py_result
        while it__r != _r.end():
           item_py_result = ProteinGroup.__new__(ProteinGroup)
           item_py_result.inst = shared_ptr[_ProteinGroup](new _ProteinGroup(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def insertProteinGroup(self, ProteinGroup group ):
        """
        insertProteinGroup(self, group: ProteinGroup ) -> None
        Appends a new protein group
        """
        assert isinstance(group, ProteinGroup), 'arg group wrong type'
    
        self.inst.get().insertProteinGroup((deref(group.inst.get())))
    
    def getIndistinguishableProteins(self):
        """
        getIndistinguishableProteins(self) -> List[ProteinGroup]
        Returns the indistinguishable proteins
        """
        _r = self.inst.get().getIndistinguishableProteins()
        py_result = []
        cdef libcpp_vector[_ProteinGroup].iterator it__r = _r.begin()
        cdef ProteinGroup item_py_result
        while it__r != _r.end():
           item_py_result = ProteinGroup.__new__(ProteinGroup)
           item_py_result.inst = shared_ptr[_ProteinGroup](new _ProteinGroup(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def insertIndistinguishableProteins(self, ProteinGroup group ):
        """
        insertIndistinguishableProteins(self, group: ProteinGroup ) -> None
        Appends new indistinguishable proteins
        """
        assert isinstance(group, ProteinGroup), 'arg group wrong type'
    
        self.inst.get().insertIndistinguishableProteins((deref(group.inst.get())))
    
    def getSignificanceThreshold(self):
        """
        getSignificanceThreshold(self) -> float
        Returns the protein significance threshold value
        """
        cdef double _r = self.inst.get().getSignificanceThreshold()
        py_result = <double>_r
        return py_result
    
    def setSignificanceThreshold(self, double value ):
        """
        setSignificanceThreshold(self, value: float ) -> None
        Sets the protein significance threshold value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setSignificanceThreshold((<double>value))
    
    def getScoreType(self):
        """
        getScoreType(self) -> Union[bytes, str, String]
        Returns the protein score type
        """
        cdef _String _r = self.inst.get().getScoreType()
        py_result = convOutputString(_r)
        return py_result
    
    def setScoreType(self,  type ):
        """
        setScoreType(self, type: Union[bytes, str, String] ) -> None
        Sets the protein score type
        """
        assert (isinstance(type, str) or isinstance(type, bytes) or isinstance(type, String)), 'arg type wrong type'
    
        self.inst.get().setScoreType(deref((convString(type)).get()))
    
    def isHigherScoreBetter(self):
        """
        isHigherScoreBetter(self) -> bool
        Returns true if a higher score represents a better score
        """
        cdef bool _r = self.inst.get().isHigherScoreBetter()
        py_result = <bool>_r
        return py_result
    
    def setHigherScoreBetter(self, bool higher_is_better ):
        """
        setHigherScoreBetter(self, higher_is_better: bool ) -> None
        Sets the orientation of the score (is higher better?)
        """
        assert isinstance(higher_is_better, pybool_t), 'arg higher_is_better wrong type'
    
        self.inst.get().setHigherScoreBetter((<bool>higher_is_better))
    
    def sort(self):
        """
        sort(self) -> None
        Sorts the protein hits according to their score
        """
        self.inst.get().sort()
    
    def computeCoverage(self, PeptideIdentificationList pep_ids ):
        """
        computeCoverage(self, pep_ids: PeptideIdentificationList ) -> None
        Compute the coverage (in percent) of all ProteinHits given PeptideHits
        """
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
    
        self.inst.get().computeCoverage((deref(pep_ids.inst.get())))
    
    def getDateTime(self):
        """
        getDateTime(self) -> DateTime
        Returns the date of the protein identification run
        """
        cdef _DateTime * _r = new _DateTime(self.inst.get().getDateTime())
        cdef DateTime py_result = DateTime.__new__(DateTime)
        py_result.inst = shared_ptr[_DateTime](_r)
        return py_result
    
    def setDateTime(self, DateTime date ):
        """
        setDateTime(self, date: DateTime ) -> None
        Sets the date of the protein identification run
        """
        assert isinstance(date, DateTime), 'arg date wrong type'
    
        self.inst.get().setDateTime((deref(date.inst.get())))
    
    def setSearchEngine(self,  search_engine ):
        """
        setSearchEngine(self, search_engine: Union[bytes, str, String] ) -> None
        Sets the search engine type
        """
        assert (isinstance(search_engine, str) or isinstance(search_engine, bytes) or isinstance(search_engine, String)), 'arg search_engine wrong type'
    
        self.inst.get().setSearchEngine(deref((convString(search_engine)).get()))
    
    def getSearchEngine(self):
        """
        getSearchEngine(self) -> Union[bytes, str, String]
        Returns the type of search engine used
        """
        cdef _String _r = self.inst.get().getSearchEngine()
        py_result = convOutputString(_r)
        return py_result
    
    def setSearchEngineVersion(self,  search_engine_version ):
        """
        setSearchEngineVersion(self, search_engine_version: Union[bytes, str, String] ) -> None
        Sets the search engine version
        """
        assert (isinstance(search_engine_version, str) or isinstance(search_engine_version, bytes) or isinstance(search_engine_version, String)), 'arg search_engine_version wrong type'
    
        self.inst.get().setSearchEngineVersion(deref((convString(search_engine_version)).get()))
    
    def getSearchEngineVersion(self):
        """
        getSearchEngineVersion(self) -> Union[bytes, str, String]
        Returns the search engine version
        """
        cdef _String _r = self.inst.get().getSearchEngineVersion()
        py_result = convOutputString(_r)
        return py_result
    
    def setSearchParameters(self, SearchParameters search_parameters ):
        """
        setSearchParameters(self, search_parameters: SearchParameters ) -> None
        Sets the search parameters
        """
        assert isinstance(search_parameters, SearchParameters), 'arg search_parameters wrong type'
    
        self.inst.get().setSearchParameters((deref(search_parameters.inst.get())))
    
    def getSearchParameters(self):
        """
        getSearchParameters(self) -> SearchParameters
        Returns the search parameters
        """
        cdef _SearchParameters * _r = new _SearchParameters(self.inst.get().getSearchParameters())
        cdef SearchParameters py_result = SearchParameters.__new__(SearchParameters)
        py_result.inst = shared_ptr[_SearchParameters](_r)
        return py_result
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        Returns the identifier
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setIdentifier(self,  id_ ):
        """
        setIdentifier(self, id_: Union[bytes, str, String] ) -> None
        Sets the identifier
        """
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
    
        self.inst.get().setIdentifier(deref((convString(id_)).get()))
    
    def _setPrimaryMSRunPath_0(self, list s ):
        """
        _setPrimaryMSRunPath_0(self, s: List[bytes] ) -> None
        Set the file paths to the primary MS runs (usually the mzML files obtained after data conversion from raw files)
        
        
        :param raw: Store paths to the raw files (or equivalent) rather than mzMLs
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
        self.inst.get().setPrimaryMSRunPath(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def _setPrimaryMSRunPath_1(self, list s , bool raw ):
        """
        _setPrimaryMSRunPath_1(self, s: List[bytes] , raw: bool ) -> None
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        assert isinstance(raw, pybool_t), 'arg raw wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
    
        self.inst.get().setPrimaryMSRunPath(deref(v0), (<bool>raw))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def setPrimaryMSRunPath(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setPrimaryMSRunPath(self, s: List[bytes] ) -> None
          :noindex:
        
        Set the file paths to the primary MS runs (usually the mzML files obtained after data conversion from raw files)
        
        
        :param raw: Store paths to the raw files (or equivalent) rather than mzMLs
        
        .. rubric:: Overload:
        .. py:function:: setPrimaryMSRunPath(self, s: List[bytes] , raw: bool ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])):
            return self._setPrimaryMSRunPath_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], pybool_t)):
            return self._setPrimaryMSRunPath_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _addPrimaryMSRunPath_0(self, list s ):
        """
        _addPrimaryMSRunPath_0(self, s: List[bytes] ) -> None
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
        self.inst.get().addPrimaryMSRunPath(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def _addPrimaryMSRunPath_1(self, list s , bool raw ):
        """
        _addPrimaryMSRunPath_1(self, s: List[bytes] , raw: bool ) -> None
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        assert isinstance(raw, pybool_t), 'arg raw wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
    
        self.inst.get().addPrimaryMSRunPath(deref(v0), (<bool>raw))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def addPrimaryMSRunPath(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: addPrimaryMSRunPath(self, s: List[bytes] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: addPrimaryMSRunPath(self, s: List[bytes] , raw: bool ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])):
            return self._addPrimaryMSRunPath_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], pybool_t)):
            return self._addPrimaryMSRunPath_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getPrimaryMSRunPath_0(self, list output ):
        """
        _getPrimaryMSRunPath_0(self, output: List[bytes] ) -> None
        """
        assert isinstance(output, list) and all(isinstance(li, bytes) for li in output), 'arg output wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in output:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getPrimaryMSRunPath(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        output[:] = replace
        del v0
    
    def _getPrimaryMSRunPath_1(self, list output , bool raw ):
        """
        _getPrimaryMSRunPath_1(self, output: List[bytes] , raw: bool ) -> None
        """
        assert isinstance(output, list) and all(isinstance(li, bytes) for li in output), 'arg output wrong type'
        assert isinstance(raw, pybool_t), 'arg raw wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in output:
           v0.push_back(_String(<char *>item0))
    
        self.inst.get().getPrimaryMSRunPath(deref(v0), (<bool>raw))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        output[:] = replace
        del v0
    
    def getPrimaryMSRunPath(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getPrimaryMSRunPath(self, output: List[bytes] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getPrimaryMSRunPath(self, output: List[bytes] , raw: bool ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])):
            return self._getPrimaryMSRunPath_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], pybool_t)):
            return self._getPrimaryMSRunPath_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
        if not isinstance(other, ProteinIdentification):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ProteinIdentification other_casted = other
        cdef ProteinIdentification self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    PeakMassType = __PeakMassType 

cdef class PurityScores:
    """
    Cython implementation of _PurityScores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PurityScores.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property total_intensity:
        def __set__(self, double total_intensity):
        
            self.inst.get().total_intensity = (<double>total_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().total_intensity
            py_result = <double>_r
            return py_result
    
    property target_intensity:
        def __set__(self, double target_intensity):
        
            self.inst.get().target_intensity = (<double>target_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().target_intensity
            py_result = <double>_r
            return py_result
    
    property signal_proportion:
        def __set__(self, double signal_proportion):
        
            self.inst.get().signal_proportion = (<double>signal_proportion)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().signal_proportion
            py_result = <double>_r
            return py_result
    
    property target_peak_count:
        def __set__(self,  target_peak_count):
        
            self.inst.get().target_peak_count = (<size_t>target_peak_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().target_peak_count
            py_result = <size_t>_r
            return py_result
    
    property interfering_peak_count:
        def __set__(self,  interfering_peak_count):
        
            self.inst.get().interfering_peak_count = (<size_t>interfering_peak_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().interfering_peak_count
            py_result = <size_t>_r
            return py_result
    
    property interfering_peaks:
        def __set__(self, MSSpectrum interfering_peaks):
        
            self.inst.get().interfering_peaks = (deref(interfering_peaks.inst.get()))
        
    
        def __get__(self):
            cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().interfering_peaks)
            cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
            py_result.inst = shared_ptr[_MSSpectrum](_r)
            return py_result
    
    def __copy__(self):
       cdef PurityScores rv = PurityScores.__new__(PurityScores)
       rv.inst = shared_ptr[_PurityScores](new _PurityScores(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PurityScores rv = PurityScores.__new__(PurityScores)
       rv.inst = shared_ptr[_PurityScores](new _PurityScores(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PurityScores](new _PurityScores())
    
    def _init_1(self, PurityScores in_0 ):
        """
        _init_1(self, in_0: PurityScores ) -> None
        """
        assert isinstance(in_0, PurityScores), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PurityScores](new _PurityScores((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PurityScores ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PurityScores)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SearchParameters:
    """
    Cython implementation of _SearchParameters

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SearchParameters.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property db:
        def __set__(self,  db):
        
            self.inst.get().db = deref((convString(db)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().db
            py_result = convOutputString(_r)
            return py_result
    
    property db_version:
        def __set__(self,  db_version):
        
            self.inst.get().db_version = deref((convString(db_version)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().db_version
            py_result = convOutputString(_r)
            return py_result
    
    property taxonomy:
        def __set__(self,  taxonomy):
        
            self.inst.get().taxonomy = deref((convString(taxonomy)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().taxonomy
            py_result = convOutputString(_r)
            return py_result
    
    property charges:
        def __set__(self,  charges):
        
            self.inst.get().charges = deref((convString(charges)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().charges
            py_result = convOutputString(_r)
            return py_result
    
    property mass_type:
        def __set__(self, int mass_type):
        
            self.inst.get().mass_type = (<_PeakMassType>mass_type)
        
    
        def __get__(self):
            cdef _PeakMassType _r = self.inst.get().mass_type
            py_result = <int>_r
            return py_result
    
    property fixed_modifications:
        def __set__(self, list fixed_modifications):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in fixed_modifications:
               v0.push_back(_String(<char *>item0))
            self.inst.get().fixed_modifications = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().fixed_modifications
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property variable_modifications:
        def __set__(self, list variable_modifications):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in variable_modifications:
               v0.push_back(_String(<char *>item0))
            self.inst.get().variable_modifications = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().variable_modifications
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property missed_cleavages:
        def __set__(self,  missed_cleavages):
        
            self.inst.get().missed_cleavages = (<unsigned int>missed_cleavages)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().missed_cleavages
            py_result = <unsigned int>_r
            return py_result
    
    property fragment_mass_tolerance:
        def __set__(self, double fragment_mass_tolerance):
        
            self.inst.get().fragment_mass_tolerance = (<double>fragment_mass_tolerance)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().fragment_mass_tolerance
            py_result = <double>_r
            return py_result
    
    property fragment_mass_tolerance_ppm:
        def __set__(self, bool fragment_mass_tolerance_ppm):
        
            self.inst.get().fragment_mass_tolerance_ppm = (<bool>fragment_mass_tolerance_ppm)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().fragment_mass_tolerance_ppm
            py_result = <bool>_r
            return py_result
    
    property precursor_mass_tolerance:
        def __set__(self, double precursor_mass_tolerance):
        
            self.inst.get().precursor_mass_tolerance = (<double>precursor_mass_tolerance)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_mass_tolerance
            py_result = <double>_r
            return py_result
    
    property precursor_mass_tolerance_ppm:
        def __set__(self, bool precursor_mass_tolerance_ppm):
        
            self.inst.get().precursor_mass_tolerance_ppm = (<bool>precursor_mass_tolerance_ppm)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().precursor_mass_tolerance_ppm
            py_result = <bool>_r
            return py_result
    
    property digestion_enzyme:
        def __set__(self, DigestionEnzymeProtein digestion_enzyme):
        
            self.inst.get().digestion_enzyme = (deref(digestion_enzyme.inst.get()))
        
    
        def __get__(self):
            cdef _DigestionEnzymeProtein * _r = new _DigestionEnzymeProtein(self.inst.get().digestion_enzyme)
            cdef DigestionEnzymeProtein py_result = DigestionEnzymeProtein.__new__(DigestionEnzymeProtein)
            py_result.inst = shared_ptr[_DigestionEnzymeProtein](_r)
            return py_result
    
    def __copy__(self):
       cdef SearchParameters rv = SearchParameters.__new__(SearchParameters)
       rv.inst = shared_ptr[_SearchParameters](new _SearchParameters(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SearchParameters rv = SearchParameters.__new__(SearchParameters)
       rv.inst = shared_ptr[_SearchParameters](new _SearchParameters(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SearchParameters](new _SearchParameters())
    
    def _init_1(self, SearchParameters in_0 ):
        """
        _init_1(self, in_0: SearchParameters ) -> None
        """
        assert isinstance(in_0, SearchParameters), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SearchParameters](new _SearchParameters((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SearchParameters ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SearchParameters)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
        if not isinstance(other, SearchParameters):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef SearchParameters other_casted = other
        cdef SearchParameters self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class SpectrumSettings:
    """
    Cython implementation of _SpectrumSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumSettings.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumSettings rv = SpectrumSettings.__new__(SpectrumSettings)
       rv.inst = shared_ptr[_SpectrumSettings](new _SpectrumSettings(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumSettings rv = SpectrumSettings.__new__(SpectrumSettings)
       rv.inst = shared_ptr[_SpectrumSettings](new _SpectrumSettings(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectrumSettings](new _SpectrumSettings())
    
    def _init_1(self, SpectrumSettings in_0 ):
        """
        _init_1(self, in_0: SpectrumSettings ) -> None
        """
        assert isinstance(in_0, SpectrumSettings), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumSettings](new _SpectrumSettings((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumSettings ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumSettings)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def unify(self, SpectrumSettings in_0 ):
        """
        unify(self, in_0: SpectrumSettings ) -> None
        """
        assert isinstance(in_0, SpectrumSettings), 'arg in_0 wrong type'
    
        self.inst.get().unify((deref(in_0.inst.get())))
    
    def getType(self):
        """
        getType(self) -> int
        Returns the spectrum type (centroided (PEAKS) or profile data (RAW))
        """
        cdef int _r = self.inst.get().getType()
        py_result = <int>_r
        return py_result
    
    def setType(self, int in_0 ):
        """
        setType(self, in_0: int ) -> None
        Sets the spectrum type
        """
        assert in_0 in [0, 1, 2, 3], 'arg in_0 wrong type'
    
        self.inst.get().setType((<_SpectrumType>in_0))
    
    def getNativeID(self):
        """
        getNativeID(self) -> Union[bytes, str, String]
        Returns the native identifier for the spectrum, used by the acquisition software
        """
        cdef _String _r = self.inst.get().getNativeID()
        py_result = convOutputString(_r)
        return py_result
    
    def setNativeID(self,  in_0 ):
        """
        setNativeID(self, in_0: Union[bytes, str, String] ) -> None
        Sets the native identifier for the spectrum, used by the acquisition software
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setNativeID(deref((convString(in_0)).get()))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the free-text comment
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  in_0 ):
        """
        setComment(self, in_0: Union[bytes, str, String] ) -> None
        Sets the free-text comment
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setComment(deref((convString(in_0)).get()))
    
    def getInstrumentSettings(self):
        """
        getInstrumentSettings(self) -> InstrumentSettings
        Returns a const reference to the instrument settings of the current spectrum
        """
        cdef _InstrumentSettings * _r = new _InstrumentSettings(self.inst.get().getInstrumentSettings())
        cdef InstrumentSettings py_result = InstrumentSettings.__new__(InstrumentSettings)
        py_result.inst = shared_ptr[_InstrumentSettings](_r)
        return py_result
    
    def setInstrumentSettings(self, InstrumentSettings in_0 ):
        """
        setInstrumentSettings(self, in_0: InstrumentSettings ) -> None
        Sets the instrument settings of the current spectrum
        """
        assert isinstance(in_0, InstrumentSettings), 'arg in_0 wrong type'
    
        self.inst.get().setInstrumentSettings((deref(in_0.inst.get())))
    
    def getAcquisitionInfo(self):
        """
        getAcquisitionInfo(self) -> AcquisitionInfo
        Returns a const reference to the acquisition info
        """
        cdef _AcquisitionInfo * _r = new _AcquisitionInfo(self.inst.get().getAcquisitionInfo())
        cdef AcquisitionInfo py_result = AcquisitionInfo.__new__(AcquisitionInfo)
        py_result.inst = shared_ptr[_AcquisitionInfo](_r)
        return py_result
    
    def setAcquisitionInfo(self, AcquisitionInfo in_0 ):
        """
        setAcquisitionInfo(self, in_0: AcquisitionInfo ) -> None
        Sets the acquisition info
        """
        assert isinstance(in_0, AcquisitionInfo), 'arg in_0 wrong type'
    
        self.inst.get().setAcquisitionInfo((deref(in_0.inst.get())))
    
    def getSourceFile(self):
        """
        getSourceFile(self) -> SourceFile
        Returns a const reference to the source file
        """
        cdef _SourceFile * _r = new _SourceFile(self.inst.get().getSourceFile())
        cdef SourceFile py_result = SourceFile.__new__(SourceFile)
        py_result.inst = shared_ptr[_SourceFile](_r)
        return py_result
    
    def setSourceFile(self, SourceFile in_0 ):
        """
        setSourceFile(self, in_0: SourceFile ) -> None
        Sets the source file
        """
        assert isinstance(in_0, SourceFile), 'arg in_0 wrong type'
    
        self.inst.get().setSourceFile((deref(in_0.inst.get())))
    
    def getPrecursors(self):
        """
        getPrecursors(self) -> List[Precursor]
        Returns a const reference to the precursors
        """
        _r = self.inst.get().getPrecursors()
        py_result = []
        cdef libcpp_vector[_Precursor].iterator it__r = _r.begin()
        cdef Precursor item_py_result
        while it__r != _r.end():
           item_py_result = Precursor.__new__(Precursor)
           item_py_result.inst = shared_ptr[_Precursor](new _Precursor(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setPrecursors(self, list in_0 ):
        """
        setPrecursors(self, in_0: List[Precursor] ) -> None
        Sets the precursors
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, Precursor) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_Precursor] * v0 = new libcpp_vector[_Precursor]()
        cdef Precursor item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setPrecursors(deref(v0))
        del v0
    
    def getProducts(self):
        """
        getProducts(self) -> List[Product]
        Returns a const reference to the products
        """
        _r = self.inst.get().getProducts()
        py_result = []
        cdef libcpp_vector[_Product].iterator it__r = _r.begin()
        cdef Product item_py_result
        while it__r != _r.end():
           item_py_result = Product.__new__(Product)
           item_py_result.inst = shared_ptr[_Product](new _Product(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setProducts(self, list in_0 ):
        """
        setProducts(self, in_0: List[Product] ) -> None
        Sets the products
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, Product) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_Product] * v0 = new libcpp_vector[_Product]()
        cdef Product item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setProducts(deref(v0))
        del v0
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        """
        _r = self.inst.get().getDataProcessing()
        py_result = []
        cdef libcpp_vector[shared_ptr[_DataProcessing]].iterator it__r = _r.begin()
        cdef DataProcessing item_py_result
        while it__r != _r.end():
           item_py_result = DataProcessing.__new__(DataProcessing)
           item_py_result.inst = deref(it__r)
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setDataProcessing(self, list in_0 ):
        """
        setDataProcessing(self, in_0: List[DataProcessing] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, DataProcessing) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[shared_ptr[_DataProcessing]] v0
        cdef DataProcessing in_0_rec
        for in_0_rec in in_0:
            v0.push_back(in_0_rec.inst)
        # call
        self.inst.get().setDataProcessing(v0)
        
    
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
        if not isinstance(other, SpectrumSettings):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef SpectrumSettings other_casted = other
        cdef SpectrumSettings self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    SpectrumType = __SpectrumType 

cdef class SwathMap:
    """
    Cython implementation of _SwathMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1SwathMap.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property lower:
        def __set__(self, double lower):
        
            self.inst.get().lower = (<double>lower)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().lower
            py_result = <double>_r
            return py_result
    
    property upper:
        def __set__(self, double upper):
        
            self.inst.get().upper = (<double>upper)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().upper
            py_result = <double>_r
            return py_result
    
    property center:
        def __set__(self, double center):
        
            self.inst.get().center = (<double>center)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().center
            py_result = <double>_r
            return py_result
    
    property ms1:
        def __set__(self, bool ms1):
        
            self.inst.get().ms1 = (<bool>ms1)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().ms1
            py_result = <bool>_r
            return py_result
    
    def __copy__(self):
       cdef SwathMap rv = SwathMap.__new__(SwathMap)
       rv.inst = shared_ptr[_SwathMap](new _SwathMap(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SwathMap rv = SwathMap.__new__(SwathMap)
       rv.inst = shared_ptr[_SwathMap](new _SwathMap(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Data structure to hold one SWATH map with information about upper / lower isolation window and whether the map is MS1 or MS2
        """
        self.inst = shared_ptr[_SwathMap](new _SwathMap())
    
    def _init_1(self, SwathMap in_0 ):
        """
        _init_1(self, in_0: SwathMap ) -> None
        """
        assert isinstance(in_0, SwathMap), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SwathMap](new _SwathMap((deref(in_0.inst.get()))))
    
    def _init_2(self, double mz_start , double mz_end , double mz_center , bool is_ms1 ):
        """
        _init_2(self, mz_start: float , mz_end: float , mz_center: float , is_ms1: bool ) -> None
        """
        assert isinstance(mz_start, float), 'arg mz_start wrong type'
        assert isinstance(mz_end, float), 'arg mz_end wrong type'
        assert isinstance(mz_center, float), 'arg mz_center wrong type'
        assert isinstance(is_ms1, pybool_t), 'arg is_ms1 wrong type'
    
    
    
    
        self.inst = shared_ptr[_SwathMap](new _SwathMap((<double>mz_start), (<double>mz_end), (<double>mz_center), (<bool>is_ms1)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Data structure to hold one SWATH map with information about upper / lower isolation window and whether the map is MS1 or MS2

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SwathMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mz_start: float , mz_end: float , mz_center: float , is_ms1: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SwathMap)):
             self._init_1(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], pybool_t)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSpectrumPtr(self):
        _r = self.inst.get().sptr

        if (_r.get() == NULL):
          return None

        cdef _SpectrumAccessOpenMS * ptr_sa = dynamic_cast[ _SpectrumAccessOpenMSPtr ](_r.get() )
        cdef _SpectrumAccessOpenMSInMemory * ptr_inmem = dynamic_cast[ _SpectrumAccessOpenMSInMemoryPtr ](_r.get() )
        cdef _SpectrumAccessOpenMSCached * ptr_cached = dynamic_cast[ _SpectrumAccessOpenMSCachedPtr ](_r.get() )
        cdef _SpectrumAccessQuadMZTransforming * ptr_quad = dynamic_cast[ _SpectrumAccessQuadMZTransformingPtr ](_r.get() )

        if (ptr_sa != NULL):
          res_sa = SpectrumAccessOpenMS(__createUnsafeObject__=True)
          res_sa.inst = dynamic_pointer_cast[_SpectrumAccessOpenMS, _ISpectrumAccess](_r)
          return res_sa
        elif (ptr_inmem != NULL):
          res_inmem = SpectrumAccessOpenMSInMemory(__createUnsafeObject__=True)
          res_inmem.inst = dynamic_pointer_cast[_SpectrumAccessOpenMSInMemory, _ISpectrumAccess](_r)
          return res_inmem
        elif (ptr_cached != NULL):
          res_cached = SpectrumAccessOpenMSCached(__createUnsafeObject__=True)
          res_cached.inst = dynamic_pointer_cast[_SpectrumAccessOpenMSCached, _ISpectrumAccess](_r)
          return res_cached
        elif (ptr_quad != NULL):
          res_quad = SpectrumAccessQuadMZTransforming(__createUnsafeObject__=True)
          res_quad.inst = dynamic_pointer_cast[_SpectrumAccessQuadMZTransforming, _ISpectrumAccess](_r)
          return res_quad
        else:
          raise Exception("Did not find suitable conversion to Python object")

    def setSpectrumPtr(self, arg):
        cdef SpectrumAccessOpenMS arg_sa 
        cdef SpectrumAccessOpenMSCached arg_cached 
        cdef SpectrumAccessOpenMSInMemory arg_inmem 
        cdef SpectrumAccessQuadMZTransforming arg_quad 
        if isinstance(arg, SpectrumAccessOpenMS):
            arg_sa = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMS](arg_sa.inst)
        elif isinstance(arg, SpectrumAccessOpenMSCached):
            arg_cached = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMSCached](arg_cached.inst)
        elif isinstance(arg, SpectrumAccessOpenMSInMemory):
            arg_inmem = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessOpenMSInMemory](arg_inmem.inst)
        elif isinstance(arg, SpectrumAccessQuadMZTransforming):
            arg_quad = arg
            self.inst.get().sptr = dynamic_pointer_cast[_ISpectrumAccess,_SpectrumAccessQuadMZTransforming](arg_quad.inst)
        else:
          raise Exception("Need to provide suitable ISpectrumAccess-derived child class") 

cdef class streampos:
    """
    Cython implementation of _streampos

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classstd_1_1streampos.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef streampos rv = streampos.__new__(streampos)
       rv.inst = shared_ptr[_streampos](new _streampos(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef streampos rv = streampos.__new__(streampos)
       rv.inst = shared_ptr[_streampos](new _streampos(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_streampos](new _streampos())
    
    def _init_1(self, streampos in_0 ):
        """
        _init_1(self, in_0: streampos ) -> None
        """
        assert isinstance(in_0, streampos), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_streampos](new _streampos((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: streampos ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], streampos)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 
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
