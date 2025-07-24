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

cdef class __ActionMode:
    None
    LOAD = 0
    STORE = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __AnalyzerType:
    None
    ANALYZERNULL = 0
    QUADRUPOLE = 1
    PAULIONTRAP = 2
    RADIALEJECTIONLINEARIONTRAP = 3
    AXIALEJECTIONLINEARIONTRAP = 4
    TOF = 5
    SECTOR = 6
    FOURIERTRANSFORM = 7
    IONSTORAGE = 8
    ESA = 9
    IT = 10
    SWIFT = 11
    CYCLOTRON = 12
    ORBITRAP = 13
    LIT = 14
    SIZE_OF_ANALYZERTYPE = 15

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __Averagines:
    None
    C = 0
    H = 1
    N = 2
    O = 3
    S = 4
    AVERAGINE_NUM = 5

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ReflectronState:
    None
    REFLSTATENULL = 0
    ON = 1
    OFF = 2
    NONE = 3
    SIZE_OF_REFLECTRONSTATE = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ResolutionMethod:
    None
    RESMETHNULL = 0
    FWHM = 1
    TENPERCENTVALLEY = 2
    BASELINE = 3
    SIZE_OF_RESOLUTIONMETHOD = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ResolutionType:
    None
    RESTYPENULL = 0
    CONSTANT = 1
    PROPORTIONAL = 2
    SIZE_OF_RESOLUTIONTYPE = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ScanDirection:
    None
    SCANDIRNULL = 0
    UP = 1
    DOWN = 2
    SIZE_OF_SCANDIRECTION = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ScanLaw:
    None
    SCANLAWNULL = 0
    EXPONENTIAL = 1
    LINEAR = 2
    QUADRATIC = 3
    SIZE_OF_SCANLAW = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class CVTerm:
    """
    Cython implementation of _CVTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTerm.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVTerm rv = CVTerm.__new__(CVTerm)
       rv.inst = shared_ptr[_CVTerm](new _CVTerm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVTerm rv = CVTerm.__new__(CVTerm)
       rv.inst = shared_ptr[_CVTerm](new _CVTerm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVTerm](new _CVTerm())
    
    def _init_1(self, CVTerm in_0 ):
        """
        _init_1(self, in_0: CVTerm ) -> None
        """
        assert isinstance(in_0, CVTerm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVTerm](new _CVTerm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVTerm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVTerm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setAccession(self,  accession ):
        """
        setAccession(self, accession: Union[bytes, str, String] ) -> None
        Sets the accession string of the term
        """
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
    
        self.inst.get().setAccession(deref((convString(accession)).get()))
    
    def getAccession(self):
        """
        getAccession(self) -> Union[bytes, str, String]
        Returns the accession string of the term
        """
        cdef _String _r = self.inst.get().getAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the term
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the term
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setCVIdentifierRef(self,  cv_id_ref ):
        """
        setCVIdentifierRef(self, cv_id_ref: Union[bytes, str, String] ) -> None
        Sets the CV identifier reference string, e.g. UO for unit obo
        """
        assert (isinstance(cv_id_ref, str) or isinstance(cv_id_ref, bytes) or isinstance(cv_id_ref, String)), 'arg cv_id_ref wrong type'
    
        self.inst.get().setCVIdentifierRef(deref((convString(cv_id_ref)).get()))
    
    def getCVIdentifierRef(self):
        """
        getCVIdentifierRef(self) -> Union[bytes, str, String]
        Returns the CV identifier reference string
        """
        cdef _String _r = self.inst.get().getCVIdentifierRef()
        py_result = convOutputString(_r)
        return py_result
    
    def getValue(self):
        """
        getValue(self) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value of the term
        """
        cdef _DataValue _r = self.inst.get().getValue()
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
    
    def setValue(self,  value ):
        """
        setValue(self, value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the value of the term
        """
        assert isinstance(value, (int, float, list, bytes, str)), 'arg value wrong type'
    
        self.inst.get().setValue(deref(DataValue(value).inst.get()))
    
    def setUnit(self, Unit unit ):
        """
        setUnit(self, unit: Unit ) -> None
        Sets the unit of the term
        """
        assert isinstance(unit, Unit), 'arg unit wrong type'
    
        self.inst.get().setUnit((deref(unit.inst.get())))
    
    def getUnit(self):
        """
        getUnit(self) -> Unit
        Returns the unit
        """
        cdef _Unit * _r = new _Unit(self.inst.get().getUnit())
        cdef Unit py_result = Unit.__new__(Unit)
        py_result.inst = shared_ptr[_Unit](_r)
        return py_result
    
    def hasValue(self):
        """
        hasValue(self) -> bool
        Checks whether the term has a value
        """
        cdef bool _r = self.inst.get().hasValue()
        py_result = <bool>_r
        return py_result
    
    def hasUnit(self):
        """
        hasUnit(self) -> bool
        Checks whether the term has a unit
        """
        cdef bool _r = self.inst.get().hasUnit()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, CVTerm):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVTerm other_casted = other
        cdef CVTerm self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class ConvexHull2D:
    """
    Cython implementation of _ConvexHull2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConvexHull2D.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ConvexHull2D rv = ConvexHull2D.__new__(ConvexHull2D)
       rv.inst = shared_ptr[_ConvexHull2D](new _ConvexHull2D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ConvexHull2D rv = ConvexHull2D.__new__(ConvexHull2D)
       rv.inst = shared_ptr[_ConvexHull2D](new _ConvexHull2D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ConvexHull2D](new _ConvexHull2D())
    
    def _init_1(self, ConvexHull2D in_0 ):
        """
        _init_1(self, in_0: ConvexHull2D ) -> None
        """
        assert isinstance(in_0, ConvexHull2D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ConvexHull2D](new _ConvexHull2D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ConvexHull2D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConvexHull2D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def clear(self):
        """
        clear(self) -> None
        Removes all points
        """
        self.inst.get().clear()
    
    def compress(self):
        """
        compress(self) -> int
        Allows to reduce the disk/memory footprint of a hull
        """
        cdef size_t _r = self.inst.get().compress()
        py_result = <size_t>_r
        return py_result
    
    def expandToBoundingBox(self):
        """
        expandToBoundingBox(self) -> None
        Expand a convex hull to its bounding box.
        """
        self.inst.get().expandToBoundingBox()
    
    def addPoint(self,  point ):
        """
        addPoint(self, point: Union[Sequence[int], Sequence[float]] ) -> bool
        Adds a point to the hull if it is not already contained. Returns if the point was added. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
        """
        assert len(point) == 2 and isinstance(point[0], (int, float)) and isinstance(point[1], (int, float)), 'arg point wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>point[0]
        _dp_0[1] = <float>point[1]
        cdef bool _r = self.inst.get().addPoint(_dp_0)
        py_result = <bool>_r
        return py_result
    
    def addPoints(self, np.ndarray[np.float32_t,ndim=2] points ):
        """
        addPoints(self, points: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
        Adds points to the hull if it is not already contained. This will trigger recomputation of the outer hull points (thus points set with setHullPoints() will be lost)
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
        self.inst.get().addPoints(_dp_vec_0)
    
    def encloses(self,  in_0 ):
        """
        encloses(self, in_0: Union[Sequence[int], Sequence[float]] ) -> bool
        Returns if the `point` lies in the feature hull
        """
        assert len(in_0) == 2 and isinstance(in_0[0], (int, float)) and isinstance(in_0[1], (int, float)), 'arg in_0 wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>in_0[0]
        _dp_0[1] = <float>in_0[1]
        cdef bool _r = self.inst.get().encloses(_dp_0)
        py_result = <bool>_r
        return py_result
    
    def getHullPoints(self):
        """
        getHullPoints(self) -> '_np.ndarray[Any, _np.dtype[_np.float32]]'
        Accessor for the outer points
        """
        cdef libcpp_vector[_DPosition2] _r = self.inst.get().getHullPoints()
        cdef int _out_n_dpos_vec = _r.size()
        cdef py_result = np.zeros([_out_n_dpos_vec,2], dtype=np.float32)
        cdef libcpp_vector[_DPosition2].iterator _out_it_dpos_vec = _r.begin()
        cdef int _out_ii_dpos_vec = 0
        while _out_it_dpos_vec != _r.end():
             py_result[_out_ii_dpos_vec, 0] = deref(_out_it_dpos_vec)[0]
             py_result[_out_ii_dpos_vec, 1] = deref(_out_it_dpos_vec)[1]
             inc(_out_it_dpos_vec)
             _out_ii_dpos_vec += 1
        return py_result
    
    def setHullPoints(self, np.ndarray[np.float32_t,ndim=2] in_0 ):
        """
        setHullPoints(self, in_0: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
        Accessor for the outer(!) points (no checking is performed if this is actually a convex hull)
        """
        assert in_0.shape[1] == 2, 'arg in_0 wrong type'
        cdef libcpp_vector[_DPosition2] _dp_vec_0
        cdef _DPosition2 _dp_0
        cdef int _dp_ii_0
        cdef int _dp_N_0 = in_0.shape[0]
        for _dp_ii_0 in range(_dp_N_0):
            _dp_0[0] = in_0[_dp_ii_0,0]
            _dp_0[1] = in_0[_dp_ii_0,1]
            _dp_vec_0.push_back(_dp_0)
        self.inst.get().setHullPoints(_dp_vec_0)
    
    def getBoundingBox(self):
        """
        getBoundingBox(self) -> DBoundingBox2
        Returns the bounding box of the feature hull points
        """
        cdef _DBoundingBox2 * _r = new _DBoundingBox2(self.inst.get().getBoundingBox())
        cdef DBoundingBox2 py_result = DBoundingBox2.__new__(DBoundingBox2)
        py_result.inst = shared_ptr[_DBoundingBox2](_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ConvexHull2D):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ConvexHull2D other_casted = other
        cdef ConvexHull2D self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    
    def enclosesXY(self, float x, float y):
        """
        Parameters:
        x (float)
        y (float)
        
        Returns:
        int
        """
        cdef _DPosition2 pos
        pos[0] = x
        pos[1] = y
        return self.inst.get().encloses(pos)


    def getHullPointsNPY(self):
        """
        Returns:
        result (np.ndarray[np.float32_t, ndim=2])
        """
        cdef libcpp_vector[_DPosition2] points = self.inst.get().getHullPoints()
        cdef np.ndarray[np.float32_t, ndim=2] result
        cdef n = points.size()
        result = np.zeros( [n,2], dtype=np.float32)
        cdef libcpp_vector[_DPosition2].iterator it = points.begin()
        cdef int i = 0
        while it != points.end():
            result[i,0] = deref(it)[0]
            result[i,1] = deref(it)[1]
            inc(it)
            i += 1
        return result

    def setHullPointsNPY(self, np.ndarray[np.float32_t, ndim=2] points):
        """
        Parameters:
        points (np.ndarray[np.float32_t, ndim=2])
        """
        cdef _ConvexHull2D * hull = self.inst.get()
        cdef int N = points.shape[0]
        cdef int i
        cdef libcpp_vector[_DPosition2] vec
        cdef _DPosition2 p
        for i in range(N):
            p[0] = points[i,0]
            p[1] = points[i,1]
            vec.push_back(p)
        self.inst.get().setHullPoints(vec)

    def getBoundingBox2D(self):
        """
        Returns:
        ((double,double),(double,double))
        """
        cdef _DBoundingBox2 box = self.inst.get().getBoundingBox()
        cdef _DPosition2 minp = box.minPosition()
        cdef _DPosition2 maxp = box.maxPosition()
        return (minp[0], minp[1]), (maxp[0], maxp[1])

    def addPointXY(self, x, y):
        """
        Parameters:
        x (double)
        y (double)
        """
        cdef _DPosition2 p
        p[0] = x
        p[1] = y
        self.inst.get().addPoint(p)

    def addPointsNPY(self, np.ndarray[np.float32_t, ndim=2] points):
        """
        Parameters:
        points (np.ndarray[np.float32_t, ndim=2])
        """
        cdef _ConvexHull2D * hull = self.inst.get()
        cdef int N = points.shape[0]
        cdef int i
        cdef libcpp_vector[_DPosition2] vec
        cdef _DPosition2 p
        for i in range(N):
            p[0] = points[i,0]
            p[1] = points[i,1]
            vec.push_back(p)
        self.inst.get().addPoints(vec) 

cdef class CrossLinkSpectrumMatch:
    """
    Cython implementation of _CrossLinkSpectrumMatch

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1CrossLinkSpectrumMatch.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property cross_link:
        def __set__(self, ProteinProteinCrossLink cross_link):
        
            self.inst.get().cross_link = (deref(cross_link.inst.get()))
        
    
        def __get__(self):
            cdef _ProteinProteinCrossLink * _r = new _ProteinProteinCrossLink(self.inst.get().cross_link)
            cdef ProteinProteinCrossLink py_result = ProteinProteinCrossLink.__new__(ProteinProteinCrossLink)
            py_result.inst = shared_ptr[_ProteinProteinCrossLink](_r)
            return py_result
    
    property scan_index_light:
        def __set__(self,  scan_index_light):
        
            self.inst.get().scan_index_light = (<size_t>scan_index_light)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().scan_index_light
            py_result = <size_t>_r
            return py_result
    
    property scan_index_heavy:
        def __set__(self,  scan_index_heavy):
        
            self.inst.get().scan_index_heavy = (<size_t>scan_index_heavy)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().scan_index_heavy
            py_result = <size_t>_r
            return py_result
    
    property score:
        def __set__(self, double score):
        
            self.inst.get().score = (<double>score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().score
            py_result = <double>_r
            return py_result
    
    property rank:
        def __set__(self,  rank):
        
            self.inst.get().rank = (<size_t>rank)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().rank
            py_result = <size_t>_r
            return py_result
    
    property xquest_score:
        def __set__(self, double xquest_score):
        
            self.inst.get().xquest_score = (<double>xquest_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xquest_score
            py_result = <double>_r
            return py_result
    
    property pre_score:
        def __set__(self, double pre_score):
        
            self.inst.get().pre_score = (<double>pre_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().pre_score
            py_result = <double>_r
            return py_result
    
    property percTIC:
        def __set__(self, double percTIC):
        
            self.inst.get().percTIC = (<double>percTIC)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().percTIC
            py_result = <double>_r
            return py_result
    
    property wTIC:
        def __set__(self, double wTIC):
        
            self.inst.get().wTIC = (<double>wTIC)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().wTIC
            py_result = <double>_r
            return py_result
    
    property wTICold:
        def __set__(self, double wTICold):
        
            self.inst.get().wTICold = (<double>wTICold)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().wTICold
            py_result = <double>_r
            return py_result
    
    property int_sum:
        def __set__(self, double int_sum):
        
            self.inst.get().int_sum = (<double>int_sum)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().int_sum
            py_result = <double>_r
            return py_result
    
    property intsum_alpha:
        def __set__(self, double intsum_alpha):
        
            self.inst.get().intsum_alpha = (<double>intsum_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intsum_alpha
            py_result = <double>_r
            return py_result
    
    property intsum_beta:
        def __set__(self, double intsum_beta):
        
            self.inst.get().intsum_beta = (<double>intsum_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intsum_beta
            py_result = <double>_r
            return py_result
    
    property total_current:
        def __set__(self, double total_current):
        
            self.inst.get().total_current = (<double>total_current)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().total_current
            py_result = <double>_r
            return py_result
    
    property precursor_error_ppm:
        def __set__(self, double precursor_error_ppm):
        
            self.inst.get().precursor_error_ppm = (<double>precursor_error_ppm)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_error_ppm
            py_result = <double>_r
            return py_result
    
    property match_odds:
        def __set__(self, double match_odds):
        
            self.inst.get().match_odds = (<double>match_odds)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().match_odds
            py_result = <double>_r
            return py_result
    
    property match_odds_alpha:
        def __set__(self, double match_odds_alpha):
        
            self.inst.get().match_odds_alpha = (<double>match_odds_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().match_odds_alpha
            py_result = <double>_r
            return py_result
    
    property match_odds_beta:
        def __set__(self, double match_odds_beta):
        
            self.inst.get().match_odds_beta = (<double>match_odds_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().match_odds_beta
            py_result = <double>_r
            return py_result
    
    property log_occupancy:
        def __set__(self, double log_occupancy):
        
            self.inst.get().log_occupancy = (<double>log_occupancy)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().log_occupancy
            py_result = <double>_r
            return py_result
    
    property log_occupancy_alpha:
        def __set__(self, double log_occupancy_alpha):
        
            self.inst.get().log_occupancy_alpha = (<double>log_occupancy_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().log_occupancy_alpha
            py_result = <double>_r
            return py_result
    
    property log_occupancy_beta:
        def __set__(self, double log_occupancy_beta):
        
            self.inst.get().log_occupancy_beta = (<double>log_occupancy_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().log_occupancy_beta
            py_result = <double>_r
            return py_result
    
    property xcorrx_max:
        def __set__(self, double xcorrx_max):
        
            self.inst.get().xcorrx_max = (<double>xcorrx_max)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xcorrx_max
            py_result = <double>_r
            return py_result
    
    property xcorrc_max:
        def __set__(self, double xcorrc_max):
        
            self.inst.get().xcorrc_max = (<double>xcorrc_max)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xcorrc_max
            py_result = <double>_r
            return py_result
    
    property matched_linear_alpha:
        def __set__(self,  matched_linear_alpha):
        
            self.inst.get().matched_linear_alpha = (<size_t>matched_linear_alpha)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().matched_linear_alpha
            py_result = <size_t>_r
            return py_result
    
    property matched_linear_beta:
        def __set__(self,  matched_linear_beta):
        
            self.inst.get().matched_linear_beta = (<size_t>matched_linear_beta)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().matched_linear_beta
            py_result = <size_t>_r
            return py_result
    
    property matched_xlink_alpha:
        def __set__(self,  matched_xlink_alpha):
        
            self.inst.get().matched_xlink_alpha = (<size_t>matched_xlink_alpha)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().matched_xlink_alpha
            py_result = <size_t>_r
            return py_result
    
    property matched_xlink_beta:
        def __set__(self,  matched_xlink_beta):
        
            self.inst.get().matched_xlink_beta = (<size_t>matched_xlink_beta)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().matched_xlink_beta
            py_result = <size_t>_r
            return py_result
    
    property num_iso_peaks_mean:
        def __set__(self, double num_iso_peaks_mean):
        
            self.inst.get().num_iso_peaks_mean = (<double>num_iso_peaks_mean)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().num_iso_peaks_mean
            py_result = <double>_r
            return py_result
    
    property num_iso_peaks_mean_linear_alpha:
        def __set__(self, double num_iso_peaks_mean_linear_alpha):
        
            self.inst.get().num_iso_peaks_mean_linear_alpha = (<double>num_iso_peaks_mean_linear_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().num_iso_peaks_mean_linear_alpha
            py_result = <double>_r
            return py_result
    
    property num_iso_peaks_mean_linear_beta:
        def __set__(self, double num_iso_peaks_mean_linear_beta):
        
            self.inst.get().num_iso_peaks_mean_linear_beta = (<double>num_iso_peaks_mean_linear_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().num_iso_peaks_mean_linear_beta
            py_result = <double>_r
            return py_result
    
    property num_iso_peaks_mean_xlinks_alpha:
        def __set__(self, double num_iso_peaks_mean_xlinks_alpha):
        
            self.inst.get().num_iso_peaks_mean_xlinks_alpha = (<double>num_iso_peaks_mean_xlinks_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().num_iso_peaks_mean_xlinks_alpha
            py_result = <double>_r
            return py_result
    
    property num_iso_peaks_mean_xlinks_beta:
        def __set__(self, double num_iso_peaks_mean_xlinks_beta):
        
            self.inst.get().num_iso_peaks_mean_xlinks_beta = (<double>num_iso_peaks_mean_xlinks_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().num_iso_peaks_mean_xlinks_beta
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_linear_alpha:
        def __set__(self, double ppm_error_abs_sum_linear_alpha):
        
            self.inst.get().ppm_error_abs_sum_linear_alpha = (<double>ppm_error_abs_sum_linear_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_linear_alpha
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_linear_beta:
        def __set__(self, double ppm_error_abs_sum_linear_beta):
        
            self.inst.get().ppm_error_abs_sum_linear_beta = (<double>ppm_error_abs_sum_linear_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_linear_beta
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_xlinks_alpha:
        def __set__(self, double ppm_error_abs_sum_xlinks_alpha):
        
            self.inst.get().ppm_error_abs_sum_xlinks_alpha = (<double>ppm_error_abs_sum_xlinks_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_xlinks_alpha
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_xlinks_beta:
        def __set__(self, double ppm_error_abs_sum_xlinks_beta):
        
            self.inst.get().ppm_error_abs_sum_xlinks_beta = (<double>ppm_error_abs_sum_xlinks_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_xlinks_beta
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_linear:
        def __set__(self, double ppm_error_abs_sum_linear):
        
            self.inst.get().ppm_error_abs_sum_linear = (<double>ppm_error_abs_sum_linear)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_linear
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_xlinks:
        def __set__(self, double ppm_error_abs_sum_xlinks):
        
            self.inst.get().ppm_error_abs_sum_xlinks = (<double>ppm_error_abs_sum_xlinks)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_xlinks
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_alpha:
        def __set__(self, double ppm_error_abs_sum_alpha):
        
            self.inst.get().ppm_error_abs_sum_alpha = (<double>ppm_error_abs_sum_alpha)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_alpha
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum_beta:
        def __set__(self, double ppm_error_abs_sum_beta):
        
            self.inst.get().ppm_error_abs_sum_beta = (<double>ppm_error_abs_sum_beta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum_beta
            py_result = <double>_r
            return py_result
    
    property ppm_error_abs_sum:
        def __set__(self, double ppm_error_abs_sum):
        
            self.inst.get().ppm_error_abs_sum = (<double>ppm_error_abs_sum)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ppm_error_abs_sum
            py_result = <double>_r
            return py_result
    
    property precursor_correction:
        def __set__(self,  precursor_correction):
        
            self.inst.get().precursor_correction = (<int>precursor_correction)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().precursor_correction
            py_result = <int>_r
            return py_result
    
    property precursor_total_intensity:
        def __set__(self, double precursor_total_intensity):
        
            self.inst.get().precursor_total_intensity = (<double>precursor_total_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_total_intensity
            py_result = <double>_r
            return py_result
    
    property precursor_target_intensity:
        def __set__(self, double precursor_target_intensity):
        
            self.inst.get().precursor_target_intensity = (<double>precursor_target_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_target_intensity
            py_result = <double>_r
            return py_result
    
    property precursor_signal_proportion:
        def __set__(self, double precursor_signal_proportion):
        
            self.inst.get().precursor_signal_proportion = (<double>precursor_signal_proportion)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_signal_proportion
            py_result = <double>_r
            return py_result
    
    property precursor_target_peak_count:
        def __set__(self,  precursor_target_peak_count):
        
            self.inst.get().precursor_target_peak_count = (<size_t>precursor_target_peak_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().precursor_target_peak_count
            py_result = <size_t>_r
            return py_result
    
    property precursor_residual_peak_count:
        def __set__(self,  precursor_residual_peak_count):
        
            self.inst.get().precursor_residual_peak_count = (<size_t>precursor_residual_peak_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().precursor_residual_peak_count
            py_result = <size_t>_r
            return py_result
    
    property frag_annotations:
        def __set__(self, list frag_annotations):
            cdef libcpp_vector[_PeptideHit_PeakAnnotation] * v0 = new libcpp_vector[_PeptideHit_PeakAnnotation]()
            cdef PeptideHit_PeakAnnotation item0
            for item0 in frag_annotations:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().frag_annotations = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().frag_annotations
            py_result = []
            cdef libcpp_vector[_PeptideHit_PeakAnnotation].iterator it__r = _r.begin()
            cdef PeptideHit_PeakAnnotation item_py_result
            while it__r != _r.end():
               item_py_result = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
               item_py_result.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property peptide_id_index:
        def __set__(self,  peptide_id_index):
        
            self.inst.get().peptide_id_index = (<size_t>peptide_id_index)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().peptide_id_index
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef CrossLinkSpectrumMatch rv = CrossLinkSpectrumMatch.__new__(CrossLinkSpectrumMatch)
       rv.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CrossLinkSpectrumMatch rv = CrossLinkSpectrumMatch.__new__(CrossLinkSpectrumMatch)
       rv.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch())
    
    def _init_1(self, CrossLinkSpectrumMatch in_0 ):
        """
        _init_1(self, in_0: CrossLinkSpectrumMatch ) -> None
        """
        assert isinstance(in_0, CrossLinkSpectrumMatch), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CrossLinkSpectrumMatch ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CrossLinkSpectrumMatch)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class CsvFile:
    """
    Cython implementation of _CsvFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CsvFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CsvFile rv = CsvFile.__new__(CsvFile)
       rv.inst = shared_ptr[_CsvFile](new _CsvFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CsvFile rv = CsvFile.__new__(CsvFile)
       rv.inst = shared_ptr[_CsvFile](new _CsvFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CsvFile](new _CsvFile())
    
    def _init_1(self, CsvFile in_0 ):
        """
        _init_1(self, in_0: CsvFile ) -> None
        """
        assert isinstance(in_0, CsvFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CsvFile](new _CsvFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CsvFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CsvFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , bytes is_ , bool ie_ ,  first_n ):
        """
        load(self, filename: Union[bytes, str, String] , is_: bytes , ie_: bool , first_n: int ) -> None
        Loads data from a text file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(is_, bytes) and len(is_) == 1, 'arg is_ wrong type'
        assert isinstance(ie_, pybool_t), 'arg ie_ wrong type'
        assert isinstance(first_n, int), 'arg first_n wrong type'
    
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (<char>((is_)[0])), (<bool>ie_), (<int>first_n))
    
    def store(self,  filename ):
        """
        store(self, filename: Union[bytes, str, String] ) -> None
        Stores the buffer's content into a file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().store(deref((convString(filename)).get()))
    
    def addRow(self, list list ):
        """
        addRow(self, list: List[bytes] ) -> None
        Add a row to the buffer
        """
        assert isinstance(list, list) and all(isinstance(li, bytes) for li in list), 'arg list wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in list:
           v0.push_back(_String(<char *>item0))
        self.inst.get().addRow(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        list[:] = replace
        del v0
    
    def clear(self):
        """
        clear(self) -> None
        Clears the buffer
        """
        self.inst.get().clear()
    
    def getRow(self,  row , list list ):
        """
        getRow(self, row: int , list: List[bytes] ) -> bool
        Writes all items from a row to list
        """
        assert isinstance(row, int), 'arg row wrong type'
        assert isinstance(list, list) and all(isinstance(li, bytes) for li in list), 'arg list wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in list:
           v1.push_back(_String(<char *>item1))
        cdef bool _r = self.inst.get().getRow((<int>row), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        list[:] = replace
        del v1
        py_result = <bool>_r
        return py_result 

cdef class FeatureFindingMetabo:
    """
    Cython implementation of _FeatureFindingMetabo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFindingMetabo.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureFindingMetabo rv = FeatureFindingMetabo.__new__(FeatureFindingMetabo)
       rv.inst = shared_ptr[_FeatureFindingMetabo](new _FeatureFindingMetabo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureFindingMetabo rv = FeatureFindingMetabo.__new__(FeatureFindingMetabo)
       rv.inst = shared_ptr[_FeatureFindingMetabo](new _FeatureFindingMetabo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Method for the assembly of mass traces belonging to the same isotope
        pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
        and isotope abundances
        """
        self.inst = shared_ptr[_FeatureFindingMetabo](new _FeatureFindingMetabo())
    
    def _init_1(self, FeatureFindingMetabo in_0 ):
        """
        _init_1(self, in_0: FeatureFindingMetabo ) -> None
        """
        assert isinstance(in_0, FeatureFindingMetabo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureFindingMetabo](new _FeatureFindingMetabo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Method for the assembly of mass traces belonging to the same isotope
        pattern, i.e., that are compatible in retention times, mass-to-charge ratios,
        and isotope abundances
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureFindingMetabo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureFindingMetabo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, list input_mtraces , FeatureMap output_featmap , list output_chromatograms ):
        """
        run(self, input_mtraces: List[Kernel_MassTrace] , output_featmap: FeatureMap , output_chromatograms: List[List[MSChromatogram]] ) -> None
        """
        assert isinstance(input_mtraces, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in input_mtraces), 'arg input_mtraces wrong type'
        assert isinstance(output_featmap, FeatureMap), 'arg output_featmap wrong type'
        assert isinstance(output_chromatograms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, MSChromatogram) for elemt_rec_rec in elemt_rec) for elemt_rec in output_chromatograms), 'arg output_chromatograms wrong type'
        cdef libcpp_vector[_Kernel_MassTrace] * v0 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item0
        for item0 in input_mtraces:
            v0.push_back(deref(item0.inst.get()))
    
        cdef libcpp_vector[libcpp_vector[_MSChromatogram]] * v2 = new libcpp_vector[libcpp_vector[_MSChromatogram]]()
        cdef libcpp_vector[_MSChromatogram] * v2_rec = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item2_rec
        cdef libcpp_vector[_MSChromatogram].iterator it_output_chromatograms_rec
        for output_chromatograms_rec in output_chromatograms:
            v2_rec.clear()
            for item2_rec in output_chromatograms_rec:
                v2_rec.push_back(deref(item2_rec.inst.get()))
            v2.push_back(deref(v2_rec))
        self.inst.get().run(deref(v0), (deref(output_featmap.inst.get())), deref(v2))
        cdef libcpp_vector[libcpp_vector[_MSChromatogram]].iterator it_output_chromatograms = v2.begin()
        replace_0 = []
        while it_output_chromatograms != v2.end():
            it_output_chromatograms_rec = deref(it_output_chromatograms).begin()
            replace_1 = []
            while it_output_chromatograms_rec != deref(it_output_chromatograms).end():
                item2_rec = MSChromatogram.__new__(MSChromatogram)
                item2_rec.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_output_chromatograms_rec)))
                replace_1.append(item2_rec)
                inc(it_output_chromatograms_rec)
            replace_0.append(replace_1)
            inc(it_output_chromatograms)
        output_chromatograms[:] = replace_0
        del v2
        del v0
    
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

cdef class IsotopeModel:
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

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsotopeModel rv = IsotopeModel.__new__(IsotopeModel)
       rv.inst = shared_ptr[_IsotopeModel](new _IsotopeModel(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopeModel rv = IsotopeModel.__new__(IsotopeModel)
       rv.inst = shared_ptr[_IsotopeModel](new _IsotopeModel(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IsotopeModel](new _IsotopeModel())
    
    def _init_1(self, IsotopeModel in_0 ):
        """
        _init_1(self, in_0: IsotopeModel ) -> None
        """
        assert isinstance(in_0, IsotopeModel), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopeModel](new _IsotopeModel((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopeModel ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopeModel)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getCharge()
        py_result = <unsigned int>_r
        return py_result
    
    def setOffset(self, double offset ):
        """
        setOffset(self, offset: float ) -> None
        Set the offset of the model
        
        The whole model will be shifted to the new offset without being computing all over
        This leaves a discrepancy which is minor in small shifts (i.e. shifting by one or two
        standard deviations) but can get significant otherwise. In that case use setParameters()
        which enforces a recomputation of the model
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset((<double>offset))
    
    def getOffset(self):
        """
        getOffset(self) -> float
        Get the offset of the model
        """
        cdef double _r = self.inst.get().getOffset()
        py_result = <double>_r
        return py_result
    
    def getFormula(self):
        """
        getFormula(self) -> EmpiricalFormula
        Return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def setSamples(self, EmpiricalFormula formula ):
        """
        setSamples(self, formula: EmpiricalFormula ) -> None
        Set sample/supporting points of interpolation
        """
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
    
        self.inst.get().setSamples((deref(formula.inst.get())))
    
    def getCenter(self):
        """
        getCenter(self) -> float
        Get the center of the Isotope model
        
        This is a m/z-value not necessarily the monoisotopic mass
        """
        cdef double _r = self.inst.get().getCenter()
        py_result = <double>_r
        return py_result
    
    def getIsotopeDistribution(self):
        """
        getIsotopeDistribution(self) -> IsotopeDistribution
        Get the Isotope distribution (without widening) from the last setSamples() call
        
        Useful to determine the number of isotopes that the model contains and their position
        """
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().getIsotopeDistribution())
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    Averagines = __Averagines 

cdef class MassAnalyzer:
    """
    Cython implementation of _MassAnalyzer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassAnalyzer.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MassAnalyzer rv = MassAnalyzer.__new__(MassAnalyzer)
       rv.inst = shared_ptr[_MassAnalyzer](new _MassAnalyzer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassAnalyzer rv = MassAnalyzer.__new__(MassAnalyzer)
       rv.inst = shared_ptr[_MassAnalyzer](new _MassAnalyzer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MassAnalyzer](new _MassAnalyzer())
    
    def _init_1(self, MassAnalyzer in_0 ):
        """
        _init_1(self, in_0: MassAnalyzer ) -> None
        """
        assert isinstance(in_0, MassAnalyzer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassAnalyzer](new _MassAnalyzer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassAnalyzer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassAnalyzer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getType(self):
        """
        getType(self) -> int
        Returns the analyzer type
        """
        cdef _AnalyzerType _r = self.inst.get().getType()
        py_result = <int>_r
        return py_result
    
    def setType(self, int type ):
        """
        setType(self, type: int ) -> None
        Sets the analyzer type
        """
        assert type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'arg type wrong type'
    
        self.inst.get().setType((<_AnalyzerType>type))
    
    def getResolutionMethod(self):
        """
        getResolutionMethod(self) -> int
        Returns the method used for determination of the resolution
        """
        cdef _ResolutionMethod _r = self.inst.get().getResolutionMethod()
        py_result = <int>_r
        return py_result
    
    def setResolutionMethod(self, int resolution_method ):
        """
        setResolutionMethod(self, resolution_method: int ) -> None
        Sets the method used for determination of the resolution
        """
        assert resolution_method in [0, 1, 2, 3, 4], 'arg resolution_method wrong type'
    
        self.inst.get().setResolutionMethod((<_ResolutionMethod>resolution_method))
    
    def getResolutionType(self):
        """
        getResolutionType(self) -> int
        Returns the resolution type
        """
        cdef _ResolutionType _r = self.inst.get().getResolutionType()
        py_result = <int>_r
        return py_result
    
    def setResolutionType(self, int resolution_type ):
        """
        setResolutionType(self, resolution_type: int ) -> None
        Sets the resolution type
        """
        assert resolution_type in [0, 1, 2, 3], 'arg resolution_type wrong type'
    
        self.inst.get().setResolutionType((<_ResolutionType>resolution_type))
    
    def getScanDirection(self):
        """
        getScanDirection(self) -> int
        Returns the direction of scanning
        """
        cdef _ScanDirection _r = self.inst.get().getScanDirection()
        py_result = <int>_r
        return py_result
    
    def setScanDirection(self, int scan_direction ):
        """
        setScanDirection(self, scan_direction: int ) -> None
        Sets the direction of scanning
        """
        assert scan_direction in [0, 1, 2, 3], 'arg scan_direction wrong type'
    
        self.inst.get().setScanDirection((<_ScanDirection>scan_direction))
    
    def getScanLaw(self):
        """
        getScanLaw(self) -> int
        Returns the scan law
        """
        cdef _ScanLaw _r = self.inst.get().getScanLaw()
        py_result = <int>_r
        return py_result
    
    def setScanLaw(self, int scan_law ):
        """
        setScanLaw(self, scan_law: int ) -> None
        Sets the scan law
        """
        assert scan_law in [0, 1, 2, 3, 4], 'arg scan_law wrong type'
    
        self.inst.get().setScanLaw((<_ScanLaw>scan_law))
    
    def getReflectronState(self):
        """
        getReflectronState(self) -> int
        Returns the reflectron state (for TOF)
        """
        cdef _ReflectronState _r = self.inst.get().getReflectronState()
        py_result = <int>_r
        return py_result
    
    def setReflectronState(self, int reflecton_state ):
        """
        setReflectronState(self, reflecton_state: int ) -> None
        Sets the reflectron state (for TOF)
        """
        assert reflecton_state in [0, 1, 2, 3, 4], 'arg reflecton_state wrong type'
    
        self.inst.get().setReflectronState((<_ReflectronState>reflecton_state))
    
    def getResolution(self):
        """
        getResolution(self) -> float
        Returns the resolution. The maximum m/z value at which two peaks can be resolved, according to one of the standard measures
        """
        cdef double _r = self.inst.get().getResolution()
        py_result = <double>_r
        return py_result
    
    def setResolution(self, double resolution ):
        """
        setResolution(self, resolution: float ) -> None
        Sets the resolution
        """
        assert isinstance(resolution, float), 'arg resolution wrong type'
    
        self.inst.get().setResolution((<double>resolution))
    
    def getAccuracy(self):
        """
        getAccuracy(self) -> float
        Returns the mass accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)
        """
        cdef double _r = self.inst.get().getAccuracy()
        py_result = <double>_r
        return py_result
    
    def setAccuracy(self, double accuracy ):
        """
        setAccuracy(self, accuracy: float ) -> None
        Sets the accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)
        """
        assert isinstance(accuracy, float), 'arg accuracy wrong type'
    
        self.inst.get().setAccuracy((<double>accuracy))
    
    def getScanRate(self):
        """
        getScanRate(self) -> float
        Returns the scan rate (in s)
        """
        cdef double _r = self.inst.get().getScanRate()
        py_result = <double>_r
        return py_result
    
    def setScanRate(self, double scan_rate ):
        """
        setScanRate(self, scan_rate: float ) -> None
        Sets the scan rate (in s)
        """
        assert isinstance(scan_rate, float), 'arg scan_rate wrong type'
    
        self.inst.get().setScanRate((<double>scan_rate))
    
    def getScanTime(self):
        """
        getScanTime(self) -> float
        Returns the scan time for a single scan (in s)
        """
        cdef double _r = self.inst.get().getScanTime()
        py_result = <double>_r
        return py_result
    
    def setScanTime(self, double scan_time ):
        """
        setScanTime(self, scan_time: float ) -> None
        Sets the scan time for a single scan (in s)
        """
        assert isinstance(scan_time, float), 'arg scan_time wrong type'
    
        self.inst.get().setScanTime((<double>scan_time))
    
    def getTOFTotalPathLength(self):
        """
        getTOFTotalPathLength(self) -> float
        Returns the path length for a TOF mass analyzer (in meter)
        """
        cdef double _r = self.inst.get().getTOFTotalPathLength()
        py_result = <double>_r
        return py_result
    
    def setTOFTotalPathLength(self, double TOF_total_path_length ):
        """
        setTOFTotalPathLength(self, TOF_total_path_length: float ) -> None
        Sets the path length for a TOF mass analyzer (in meter)
        """
        assert isinstance(TOF_total_path_length, float), 'arg TOF_total_path_length wrong type'
    
        self.inst.get().setTOFTotalPathLength((<double>TOF_total_path_length))
    
    def getIsolationWidth(self):
        """
        getIsolationWidth(self) -> float
        Returns the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
        """
        cdef double _r = self.inst.get().getIsolationWidth()
        py_result = <double>_r
        return py_result
    
    def setIsolationWidth(self, double isolation_width ):
        """
        setIsolationWidth(self, isolation_width: float ) -> None
        Sets the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
        """
        assert isinstance(isolation_width, float), 'arg isolation_width wrong type'
    
        self.inst.get().setIsolationWidth((<double>isolation_width))
    
    def getFinalMSExponent(self):
        """
        getFinalMSExponent(self) -> int
        Returns the final MS exponent
        """
        cdef int _r = self.inst.get().getFinalMSExponent()
        py_result = <int>_r
        return py_result
    
    def setFinalMSExponent(self,  final_MS_exponent ):
        """
        setFinalMSExponent(self, final_MS_exponent: int ) -> None
        Sets the final MS exponent
        """
        assert isinstance(final_MS_exponent, int), 'arg final_MS_exponent wrong type'
    
        self.inst.get().setFinalMSExponent((<int>final_MS_exponent))
    
    def getMagneticFieldStrength(self):
        """
        getMagneticFieldStrength(self) -> float
        Returns the strength of the magnetic field (in T)
        """
        cdef double _r = self.inst.get().getMagneticFieldStrength()
        py_result = <double>_r
        return py_result
    
    def setMagneticFieldStrength(self, double magnetic_field_strength ):
        """
        setMagneticFieldStrength(self, magnetic_field_strength: float ) -> None
        Sets the strength of the magnetic field (in T)
        """
        assert isinstance(magnetic_field_strength, float), 'arg magnetic_field_strength wrong type'
    
        self.inst.get().setMagneticFieldStrength((<double>magnetic_field_strength))
    
    def getOrder(self):
        """
        getOrder(self) -> int
        Returns the position of this part in the whole Instrument
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
        if not isinstance(other, MassAnalyzer):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MassAnalyzer other_casted = other
        cdef MassAnalyzer self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    AnalyzerType = __AnalyzerType
    ReflectronState = __ReflectronState
    ResolutionMethod = __ResolutionMethod
    ResolutionType = __ResolutionType
    ScanDirection = __ScanDirection
    ScanLaw = __ScanLaw 

cdef class MzDataFile:
    """
    Cython implementation of _MzDataFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzDataFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzDataFile rv = MzDataFile.__new__(MzDataFile)
       rv.inst = shared_ptr[_MzDataFile](new _MzDataFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzDataFile rv = MzDataFile.__new__(MzDataFile)
       rv.inst = shared_ptr[_MzDataFile](new _MzDataFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        File adapter for MzData files
        """
        self.inst = shared_ptr[_MzDataFile](new _MzDataFile())
    
    def _init_1(self, MzDataFile in_0 ):
        """
        _init_1(self, in_0: MzDataFile ) -> None
        """
        assert isinstance(in_0, MzDataFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzDataFile](new _MzDataFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        File adapter for MzData files

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzDataFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzDataFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment map ):
        """
        load(self, filename: Union[bytes, str, String] , map: MSExperiment ) -> None
        Loads a map from a MzData file
        
        
        :param filename: Directory of the file with the file name
        :param map: It has to be a MSExperiment or have the same interface
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map, MSExperiment), 'arg map wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(map.inst.get())))
    
    def store(self,  filename , MSExperiment map ):
        """
        store(self, filename: Union[bytes, str, String] , map: MSExperiment ) -> None
        Stores a map in a MzData file
        
        
        :param filename: Directory of the file with the file name
        :param map: It has to be a MSExperiment or have the same interface
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be created
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map, MSExperiment), 'arg map wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(map.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        Returns the options for loading/storing
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
        return py_result
    
    def setOptions(self, PeakFileOptions in_0 ):
        """
        setOptions(self, in_0: PeakFileOptions ) -> None
        Sets options for loading/storing
        """
        assert isinstance(in_0, PeakFileOptions), 'arg in_0 wrong type'
    
        self.inst.get().setOptions((deref(in_0.inst.get())))
    
    def isSemanticallyValid(self,  filename , list errors , list warnings ):
        """
        isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool
        Checks if a file is valid with respect to the mapping file and the controlled vocabulary
        
        
        :param filename: File name of the file to be checked
        :param errors: Errors during the validation are returned in this output parameter
        :param warnings: Warnings during the validation are returned in this output parameter
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(errors, list) and all(isinstance(li, bytes) for li in errors), 'arg errors wrong type'
        assert isinstance(warnings, list) and all(isinstance(li, bytes) for li in warnings), 'arg warnings wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in errors:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in warnings:
           v2.push_back(_String(<char *>item2))
        cdef bool _r = self.inst.get().isSemanticallyValid(deref((convString(filename)).get()), deref(v1), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        warnings[:] = replace
        del v2
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        errors[:] = replace
        del v1
        py_result = <bool>_r
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

cdef class NucleicAcidSpectrumGenerator:
    """
    Cython implementation of _NucleicAcidSpectrumGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NucleicAcidSpectrumGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef NucleicAcidSpectrumGenerator rv = NucleicAcidSpectrumGenerator.__new__(NucleicAcidSpectrumGenerator)
       rv.inst = shared_ptr[_NucleicAcidSpectrumGenerator](new _NucleicAcidSpectrumGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NucleicAcidSpectrumGenerator rv = NucleicAcidSpectrumGenerator.__new__(NucleicAcidSpectrumGenerator)
       rv.inst = shared_ptr[_NucleicAcidSpectrumGenerator](new _NucleicAcidSpectrumGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NucleicAcidSpectrumGenerator](new _NucleicAcidSpectrumGenerator())
    
    def _init_1(self, NucleicAcidSpectrumGenerator in_0 ):
        """
        _init_1(self, in_0: NucleicAcidSpectrumGenerator ) -> None
        """
        assert isinstance(in_0, NucleicAcidSpectrumGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NucleicAcidSpectrumGenerator](new _NucleicAcidSpectrumGenerator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NucleicAcidSpectrumGenerator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NucleicAcidSpectrumGenerator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSpectrum(self, MSSpectrum spec , NASequence oligo ,  min_charge ,  max_charge ):
        """
        getSpectrum(self, spec: MSSpectrum , oligo: NASequence , min_charge: int , max_charge: int ) -> None
        Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(oligo, NASequence), 'arg oligo wrong type'
        assert isinstance(min_charge, int), 'arg min_charge wrong type'
        assert isinstance(max_charge, int), 'arg max_charge wrong type'
    
    
    
    
        self.inst.get().getSpectrum((deref(spec.inst.get())), (deref(oligo.inst.get())), (<int>min_charge), (<int>max_charge))
    
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

cdef class OnDiscMSExperiment:
    """
    Cython implementation of _OnDiscMSExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OnDiscMSExperiment.html>`_

    Representation of a mass spectrometry experiment on disk.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OnDiscMSExperiment rv = OnDiscMSExperiment.__new__(OnDiscMSExperiment)
       rv.inst = shared_ptr[_OnDiscMSExperiment](new _OnDiscMSExperiment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OnDiscMSExperiment rv = OnDiscMSExperiment.__new__(OnDiscMSExperiment)
       rv.inst = shared_ptr[_OnDiscMSExperiment](new _OnDiscMSExperiment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OnDiscMSExperiment](new _OnDiscMSExperiment())
    
    def _init_1(self, OnDiscMSExperiment in_0 ):
        """
        _init_1(self, in_0: OnDiscMSExperiment ) -> None
        """
        assert isinstance(in_0, OnDiscMSExperiment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OnDiscMSExperiment](new _OnDiscMSExperiment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OnDiscMSExperiment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OnDiscMSExperiment)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _openFile_0(self,  filename ):
        """
        _openFile_0(self, filename: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        cdef bool _r = self.inst.get().openFile(deref((convString(filename)).get()))
        py_result = <bool>_r
        return py_result
    
    def _openFile_1(self,  filename , bool skipLoadingMetaData ):
        """
        _openFile_1(self, filename: Union[bytes, str, String] , skipLoadingMetaData: bool ) -> bool
        Open a specific file on disk
        
        This tries to read the indexed mzML by parsing the index and then reading the meta information into memory
        
        returns: Whether the parsing of the file was successful (if false, the file most likely was not an indexed mzML file)
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(skipLoadingMetaData, pybool_t), 'arg skipLoadingMetaData wrong type'
    
    
        cdef bool _r = self.inst.get().openFile(deref((convString(filename)).get()), (<bool>skipLoadingMetaData))
        py_result = <bool>_r
        return py_result
    
    def openFile(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: openFile(self, filename: Union[bytes, str, String] ) -> bool
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: openFile(self, filename: Union[bytes, str, String] , skipLoadingMetaData: bool ) -> bool
          :noindex:
        
        Open a specific file on disk
        
        This tries to read the indexed mzML by parsing the index and then reading the meta information into memory
        
        returns: Whether the parsing of the file was successful (if false, the file most likely was not an indexed mzML file)
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._openFile_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)):
            return self._openFile_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getNrSpectra(self):
        """
        getNrSpectra(self) -> int
        Returns the total number of spectra available
        """
        cdef size_t _r = self.inst.get().getNrSpectra()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatograms(self):
        """
        getNrChromatograms(self) -> int
        Returns the total number of chromatograms available
        """
        cdef size_t _r = self.inst.get().getNrChromatograms()
        py_result = <size_t>_r
        return py_result
    
    def getExperimentalSettings(self):
        """
        getExperimentalSettings(self) -> ExperimentalSettings
        Returns the meta information of this experiment (const access)
        """
        cdef shared_ptr[const _ExperimentalSettings] _r = self.inst.get().getExperimentalSettings()
        # Const shared_ptr detected, we need to produce a non-const copy to stick into Python object
        cdef _ExperimentalSettings * raw__r = new _ExperimentalSettings((deref(<_ExperimentalSettings * const>_r.get())))
        cdef ExperimentalSettings py_result
        py_result = ExperimentalSettings.__new__(ExperimentalSettings)
        py_result.inst = shared_ptr[_ExperimentalSettings](raw__r)
        return py_result
    
    def getMetaData(self):
        """
        getMetaData(self) -> MSExperiment
        Returns the meta information of this experiment
        """
        cdef shared_ptr[_MSExperiment] _r = self.inst.get().getMetaData()
        cdef MSExperiment py_result
        py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = _r
        return py_result
    
    def getSpectrum(self,  id ):
        """
        getSpectrum(self, id: int ) -> MSSpectrum
        Returns a single spectrum
        
        
        :param id: The index of the spectrum
        """
        assert isinstance(id, int) and id >= 0, 'arg id wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>id)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def getSpectrumByNativeId(self,  id ):
        """
        getSpectrumByNativeId(self, id: Union[bytes, str, String] ) -> MSSpectrum
        Returns a single spectrum
        
        
        :param id: The native identifier of the spectrum
        """
        assert (isinstance(id, str) or isinstance(id, bytes) or isinstance(id, String)), 'arg id wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrumByNativeId(deref((convString(id)).get())))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def getChromatogram(self,  id ):
        """
        getChromatogram(self, id: int ) -> MSChromatogram
        Returns a single chromatogram
        
        
        :param id: The index of the chromatogram
        """
        assert isinstance(id, int) and id >= 0, 'arg id wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>id)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def getChromatogramByNativeId(self,  id ):
        """
        getChromatogramByNativeId(self, id: Union[bytes, str, String] ) -> MSChromatogram
        Returns a single chromatogram
        
        
        :param id: The native identifier of the chromatogram
        """
        assert (isinstance(id, str) or isinstance(id, bytes) or isinstance(id, String)), 'arg id wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogramByNativeId(deref((convString(id)).get())))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def getSpectrumById(self,  id_ ):
        """
        getSpectrumById(self, id_: int ) -> _Interfaces_Spectrum
        Returns a single spectrum
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_Spectrum] _r = self.inst.get().getSpectrumById((<int>id_))
        cdef _Interfaces_Spectrum py_result
        py_result = _Interfaces_Spectrum.__new__(_Interfaces_Spectrum)
        py_result.inst = _r
        return py_result
    
    def getChromatogramById(self,  id_ ):
        """
        getChromatogramById(self, id_: int ) -> _Interfaces_Chromatogram
        Returns a single chromatogram
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_Chromatogram] _r = self.inst.get().getChromatogramById((<int>id_))
        cdef _Interfaces_Chromatogram py_result
        py_result = _Interfaces_Chromatogram.__new__(_Interfaces_Chromatogram)
        py_result.inst = _r
        return py_result
    
    def setSkipXMLChecks(self, bool skip ):
        """
        setSkipXMLChecks(self, skip: bool ) -> None
        Sets whether to skip some XML checks and be fast instead
        """
        assert isinstance(skip, pybool_t), 'arg skip wrong type'
    
        self.inst.get().setSkipXMLChecks((<bool>skip)) 

cdef class PeakFileOptions:
    """
    Cython implementation of _PeakFileOptions

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakFileOptions.html>`_

    Options for loading files containing peak data
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakFileOptions rv = PeakFileOptions.__new__(PeakFileOptions)
       rv.inst = shared_ptr[_PeakFileOptions](new _PeakFileOptions(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakFileOptions rv = PeakFileOptions.__new__(PeakFileOptions)
       rv.inst = shared_ptr[_PeakFileOptions](new _PeakFileOptions(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakFileOptions](new _PeakFileOptions())
    
    def _init_1(self, PeakFileOptions in_0 ):
        """
        _init_1(self, in_0: PeakFileOptions ) -> None
        """
        assert isinstance(in_0, PeakFileOptions), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakFileOptions](new _PeakFileOptions((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakFileOptions ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakFileOptions)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMetadataOnly(self, bool in_0 ):
        """
        setMetadataOnly(self, in_0: bool ) -> None
        Sets whether or not to load only meta data
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setMetadataOnly((<bool>in_0))
    
    def getMetadataOnly(self):
        """
        getMetadataOnly(self) -> bool
        Returns whether or not to load only meta data
        """
        cdef bool _r = self.inst.get().getMetadataOnly()
        py_result = <bool>_r
        return py_result
    
    def setWriteSupplementalData(self, bool in_0 ):
        """
        setWriteSupplementalData(self, in_0: bool ) -> None
        Sets whether or not to write supplemental peak data in MzData files
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setWriteSupplementalData((<bool>in_0))
    
    def getWriteSupplementalData(self):
        """
        getWriteSupplementalData(self) -> bool
        Returns whether or not to write supplemental peak data in MzData files
        """
        cdef bool _r = self.inst.get().getWriteSupplementalData()
        py_result = <bool>_r
        return py_result
    
    def setMSLevels(self, list levels ):
        """
        setMSLevels(self, levels: List[int] ) -> None
        Sets the desired MS levels for peaks to load
        """
        assert isinstance(levels, list) and all(isinstance(elemt_rec, int) for elemt_rec in levels), 'arg levels wrong type'
        cdef libcpp_vector[int] v0 = levels
        self.inst.get().setMSLevels(v0)
        
    
    def addMSLevel(self,  level ):
        """
        addMSLevel(self, level: int ) -> None
        Adds a desired MS level for peaks to load
        """
        assert isinstance(level, int), 'arg level wrong type'
    
        self.inst.get().addMSLevel((<int>level))
    
    def clearMSLevels(self):
        """
        clearMSLevels(self) -> None
        Clears the MS levels
        """
        self.inst.get().clearMSLevels()
    
    def hasMSLevels(self):
        """
        hasMSLevels(self) -> bool
        Returns true, if MS levels have been set
        """
        cdef bool _r = self.inst.get().hasMSLevels()
        py_result = <bool>_r
        return py_result
    
    def containsMSLevel(self,  level ):
        """
        containsMSLevel(self, level: int ) -> bool
        Returns true, if MS level `level` has been set
        """
        assert isinstance(level, int), 'arg level wrong type'
    
        cdef bool _r = self.inst.get().containsMSLevel((<int>level))
        py_result = <bool>_r
        return py_result
    
    def getMSLevels(self):
        """
        getMSLevels(self) -> List[int]
        Returns the set MS levels
        """
        _r = self.inst.get().getMSLevels()
        cdef list py_result = _r
        return py_result
    
    def setCompression(self, bool in_0 ):
        """
        setCompression(self, in_0: bool ) -> None
        Sets if data should be compressed when writing
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setCompression((<bool>in_0))
    
    def getCompression(self):
        """
        getCompression(self) -> bool
        Returns true, if data should be compressed when writing
        """
        cdef bool _r = self.inst.get().getCompression()
        py_result = <bool>_r
        return py_result
    
    def setMz32Bit(self, bool mz_32_bit ):
        """
        setMz32Bit(self, mz_32_bit: bool ) -> None
        Sets if mz-data and rt-data should be stored with 32bit or 64bit precision
        """
        assert isinstance(mz_32_bit, pybool_t), 'arg mz_32_bit wrong type'
    
        self.inst.get().setMz32Bit((<bool>mz_32_bit))
    
    def getMz32Bit(self):
        """
        getMz32Bit(self) -> bool
        Returns true, if mz-data and rt-data should be stored with 32bit precision
        """
        cdef bool _r = self.inst.get().getMz32Bit()
        py_result = <bool>_r
        return py_result
    
    def setIntensity32Bit(self, bool int_32_bit ):
        """
        setIntensity32Bit(self, int_32_bit: bool ) -> None
        Sets if intensity data should be stored with 32bit or 64bit precision
        """
        assert isinstance(int_32_bit, pybool_t), 'arg int_32_bit wrong type'
    
        self.inst.get().setIntensity32Bit((<bool>int_32_bit))
    
    def getIntensity32Bit(self):
        """
        getIntensity32Bit(self) -> bool
        Returns true, if intensity data should be stored with 32bit precision
        """
        cdef bool _r = self.inst.get().getIntensity32Bit()
        py_result = <bool>_r
        return py_result
    
    def setRTRange(self, DRange1 range_ ):
        """
        setRTRange(self, range_: DRange1 ) -> None
        Restricts the range of RT values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setRTRange((deref(range_.inst.get())))
    
    def hasRTRange(self):
        """
        hasRTRange(self) -> bool
        Returns true if an RT range has been set
        """
        cdef bool _r = self.inst.get().hasRTRange()
        py_result = <bool>_r
        return py_result
    
    def getRTRange(self):
        """
        getRTRange(self) -> DRange1
        Returns the RT range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getRTRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def setMZRange(self, DRange1 range_ ):
        """
        setMZRange(self, range_: DRange1 ) -> None
        Restricts the range of MZ values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setMZRange((deref(range_.inst.get())))
    
    def hasMZRange(self):
        """
        hasMZRange(self) -> bool
        Returns true if an MZ range has been set
        """
        cdef bool _r = self.inst.get().hasMZRange()
        py_result = <bool>_r
        return py_result
    
    def getMZRange(self):
        """
        getMZRange(self) -> DRange1
        Returns the MZ range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getMZRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def setIntensityRange(self, DRange1 range_ ):
        """
        setIntensityRange(self, range_: DRange1 ) -> None
        Restricts the range of intensity values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setIntensityRange((deref(range_.inst.get())))
    
    def hasIntensityRange(self):
        """
        hasIntensityRange(self) -> bool
        Returns true if an intensity range has been set
        """
        cdef bool _r = self.inst.get().hasIntensityRange()
        py_result = <bool>_r
        return py_result
    
    def getIntensityRange(self):
        """
        getIntensityRange(self) -> DRange1
        Returns the intensity range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getIntensityRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def getMaxDataPoolSize(self):
        """
        getMaxDataPoolSize(self) -> int
        Returns maximal size of the data pool
        """
        cdef size_t _r = self.inst.get().getMaxDataPoolSize()
        py_result = <size_t>_r
        return py_result
    
    def setMaxDataPoolSize(self,  s ):
        """
        setMaxDataPoolSize(self, s: int ) -> None
        Sets maximal size of the data pool
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
    
        self.inst.get().setMaxDataPoolSize((<size_t>s))
    
    def setSortSpectraByMZ(self, bool doSort ):
        """
        setSortSpectraByMZ(self, doSort: bool ) -> None
        Sets whether or not to sort peaks in spectra
        """
        assert isinstance(doSort, pybool_t), 'arg doSort wrong type'
    
        self.inst.get().setSortSpectraByMZ((<bool>doSort))
    
    def getSortSpectraByMZ(self):
        """
        getSortSpectraByMZ(self) -> bool
        Returns whether or not peaks in spectra should be sorted
        """
        cdef bool _r = self.inst.get().getSortSpectraByMZ()
        py_result = <bool>_r
        return py_result
    
    def setSortChromatogramsByRT(self, bool doSort ):
        """
        setSortChromatogramsByRT(self, doSort: bool ) -> None
        Sets whether or not to sort peaks in chromatograms
        """
        assert isinstance(doSort, pybool_t), 'arg doSort wrong type'
    
        self.inst.get().setSortChromatogramsByRT((<bool>doSort))
    
    def getSortChromatogramsByRT(self):
        """
        getSortChromatogramsByRT(self) -> bool
        Returns whether or not peaks in chromatograms should be sorted
        """
        cdef bool _r = self.inst.get().getSortChromatogramsByRT()
        py_result = <bool>_r
        return py_result
    
    def hasFilters(self):
        """
        hasFilters(self) -> bool
        """
        cdef bool _r = self.inst.get().hasFilters()
        py_result = <bool>_r
        return py_result
    
    def setFillData(self, bool only ):
        """
        setFillData(self, only: bool ) -> None
        Sets whether to fill the actual data into the container (spectrum/chromatogram)
        """
        assert isinstance(only, pybool_t), 'arg only wrong type'
    
        self.inst.get().setFillData((<bool>only))
    
    def getFillData(self):
        """
        getFillData(self) -> bool
        Returns whether to fill the actual data into the container (spectrum/chromatogram)
        """
        cdef bool _r = self.inst.get().getFillData()
        py_result = <bool>_r
        return py_result
    
    def setSkipXMLChecks(self, bool only ):
        """
        setSkipXMLChecks(self, only: bool ) -> None
        Sets whether to skip some XML checks and be fast instead
        """
        assert isinstance(only, pybool_t), 'arg only wrong type'
    
        self.inst.get().setSkipXMLChecks((<bool>only))
    
    def getSkipXMLChecks(self):
        """
        getSkipXMLChecks(self) -> bool
        Returns whether to skip some XML checks and be fast instead
        """
        cdef bool _r = self.inst.get().getSkipXMLChecks()
        py_result = <bool>_r
        return py_result
    
    def getWriteIndex(self):
        """
        getWriteIndex(self) -> bool
        Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        """
        cdef bool _r = self.inst.get().getWriteIndex()
        py_result = <bool>_r
        return py_result
    
    def setWriteIndex(self, bool write_index ):
        """
        setWriteIndex(self, write_index: bool ) -> None
        Returns whether to write an index at the end of the file (e.g. indexedmzML file format)
        """
        assert isinstance(write_index, pybool_t), 'arg write_index wrong type'
    
        self.inst.get().setWriteIndex((<bool>write_index))
    
    def getNumpressConfigurationMassTime(self):
        """
        getNumpressConfigurationMassTime(self) -> NumpressConfig
        Sets numpress configuration options for m/z or rt dimension
        """
        cdef _NumpressConfig * _r = new _NumpressConfig(self.inst.get().getNumpressConfigurationMassTime())
        cdef NumpressConfig py_result = NumpressConfig.__new__(NumpressConfig)
        py_result.inst = shared_ptr[_NumpressConfig](_r)
        return py_result
    
    def setNumpressConfigurationMassTime(self, NumpressConfig config ):
        """
        setNumpressConfigurationMassTime(self, config: NumpressConfig ) -> None
        Returns numpress configuration options for m/z or rt dimension
        """
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
    
        self.inst.get().setNumpressConfigurationMassTime((deref(config.inst.get())))
    
    def getNumpressConfigurationIntensity(self):
        """
        getNumpressConfigurationIntensity(self) -> NumpressConfig
        Sets numpress configuration options for intensity dimension
        """
        cdef _NumpressConfig * _r = new _NumpressConfig(self.inst.get().getNumpressConfigurationIntensity())
        cdef NumpressConfig py_result = NumpressConfig.__new__(NumpressConfig)
        py_result.inst = shared_ptr[_NumpressConfig](_r)
        return py_result
    
    def setNumpressConfigurationIntensity(self, NumpressConfig config ):
        """
        setNumpressConfigurationIntensity(self, config: NumpressConfig ) -> None
        Returns numpress configuration options for intensity dimension
        """
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
    
        self.inst.get().setNumpressConfigurationIntensity((deref(config.inst.get())))
    
    def getNumpressConfigurationFloatDataArray(self):
        """
        getNumpressConfigurationFloatDataArray(self) -> NumpressConfig
        Sets numpress configuration options for float data arrays
        """
        cdef _NumpressConfig * _r = new _NumpressConfig(self.inst.get().getNumpressConfigurationFloatDataArray())
        cdef NumpressConfig py_result = NumpressConfig.__new__(NumpressConfig)
        py_result.inst = shared_ptr[_NumpressConfig](_r)
        return py_result
    
    def setNumpressConfigurationFloatDataArray(self, NumpressConfig config ):
        """
        setNumpressConfigurationFloatDataArray(self, config: NumpressConfig ) -> None
        Returns numpress configuration options for float data arrays
        """
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
    
        self.inst.get().setNumpressConfigurationFloatDataArray((deref(config.inst.get())))
    
    def setForceMQCompatability(self, bool forceMQ ):
        """
        setForceMQCompatability(self, forceMQ: bool ) -> None
        [mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        """
        assert isinstance(forceMQ, pybool_t), 'arg forceMQ wrong type'
    
        self.inst.get().setForceMQCompatability((<bool>forceMQ))
    
    def getForceMQCompatability(self):
        """
        getForceMQCompatability(self) -> bool
        [mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
        """
        cdef bool _r = self.inst.get().getForceMQCompatability()
        py_result = <bool>_r
        return py_result
    
    def setForceTPPCompatability(self, bool forceTPP ):
        """
        setForceTPPCompatability(self, forceTPP: bool ) -> None
        [ mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        """
        assert isinstance(forceTPP, pybool_t), 'arg forceTPP wrong type'
    
        self.inst.get().setForceTPPCompatability((<bool>forceTPP))
    
    def getForceTPPCompatability(self):
        """
        getForceTPPCompatability(self) -> bool
        [mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
        """
        cdef bool _r = self.inst.get().getForceTPPCompatability()
        py_result = <bool>_r
        return py_result 

cdef class SequestOutfile:
    """
    Cython implementation of _SequestOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SequestOutfile.html>`_

    Representation of a Sequest output file
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SequestOutfile rv = SequestOutfile.__new__(SequestOutfile)
       rv.inst = shared_ptr[_SequestOutfile](new _SequestOutfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SequestOutfile rv = SequestOutfile.__new__(SequestOutfile)
       rv.inst = shared_ptr[_SequestOutfile](new _SequestOutfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Representation of a Sequest output file
        """
        self.inst = shared_ptr[_SequestOutfile](new _SequestOutfile())
    
    def _init_1(self, SequestOutfile in_0 ):
        """
        _init_1(self, in_0: SequestOutfile ) -> None
        """
        assert isinstance(in_0, SequestOutfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SequestOutfile](new _SequestOutfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Representation of a Sequest output file

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SequestOutfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SequestOutfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  result_filename , PeptideIdentificationList peptide_identifications , ProteinIdentification protein_identification , double p_value_threshold , list pvalues ,  database , bool ignore_proteins_per_peptide ):
        """
        load(self, result_filename: Union[bytes, str, String] , peptide_identifications: PeptideIdentificationList , protein_identification: ProteinIdentification , p_value_threshold: float , pvalues: List[float] , database: Union[bytes, str, String] , ignore_proteins_per_peptide: bool ) -> None
        Loads data from a Sequest outfile
        
        :param result_filename: The file to be loaded
        :param peptide_identifications: The identifications
        :param protein_identification: The protein identifications
        :param p_value_threshold: The significance level (for the peptide hit scores)
        :param pvalues: A list with the pvalues of the peptides (pvalues computed with peptide prophet)
        :param database: The database used for the search
        :param ignore_proteins_per_peptide: This is a hack to deal with files that use a suffix like "+1" in column "Reference", but do not actually list extra protein references in subsequent lines
        """
        assert (isinstance(result_filename, str) or isinstance(result_filename, bytes) or isinstance(result_filename, String)), 'arg result_filename wrong type'
        assert isinstance(peptide_identifications, PeptideIdentificationList), 'arg peptide_identifications wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(p_value_threshold, float), 'arg p_value_threshold wrong type'
        assert isinstance(pvalues, list) and all(isinstance(elemt_rec, float) for elemt_rec in pvalues), 'arg pvalues wrong type'
        assert (isinstance(database, str) or isinstance(database, bytes) or isinstance(database, String)), 'arg database wrong type'
        assert isinstance(ignore_proteins_per_peptide, pybool_t), 'arg ignore_proteins_per_peptide wrong type'
    
    
    
    
        cdef libcpp_vector[double] v4 = pvalues
    
    
        self.inst.get().load(deref((convString(result_filename)).get()), (deref(peptide_identifications.inst.get())), (deref(protein_identification.inst.get())), (<double>p_value_threshold), v4, deref((convString(database)).get()), (<bool>ignore_proteins_per_peptide))
        pvalues[:] = v4
    
    def getColumns(self,  line , list substrings ,  number_of_columns ,  reference_column ):
        """
        getColumns(self, line: Union[bytes, str, String] , substrings: List[bytes] , number_of_columns: int , reference_column: int ) -> bool
        Retrieves columns from a Sequest outfile line
        """
        assert (isinstance(line, str) or isinstance(line, bytes) or isinstance(line, String)), 'arg line wrong type'
        assert isinstance(substrings, list) and all(isinstance(i, bytes) for i in substrings), 'arg substrings wrong type'
        assert isinstance(number_of_columns, int) and number_of_columns >= 0, 'arg number_of_columns wrong type'
        assert isinstance(reference_column, int) and reference_column >= 0, 'arg reference_column wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in substrings:
           v1.push_back(_String(<char *>item1))
    
    
        cdef bool _r = self.inst.get().getColumns(deref((convString(line)).get()), deref(v1), (<size_t>number_of_columns), (<size_t>reference_column))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        substrings[:] = replace
        del v1
        py_result = <bool>_r
        return py_result
    
    def getACAndACType(self,  line ,  accession ,  accession_type ):
        """
        getACAndACType(self, line: Union[bytes, str, String] , accession: String , accession_type: String ) -> None
        Retrieves the accession type and accession number from a protein description line
        """
        assert (isinstance(line, str) or isinstance(line, bytes) or isinstance(line, String)), 'arg line wrong type'
        assert isinstance(accession, String), 'arg accession wrong type'
        assert isinstance(accession_type, String), 'arg accession_type wrong type'
    
    
    
        self.inst.get().getACAndACType(deref((convString(line)).get()), deref((<String>accession).inst.get()), deref((<String>accession_type).inst.get()))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, SequestOutfile):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef SequestOutfile other_casted = other
        cdef SequestOutfile self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    

    def getSequences(self, String database_filename , dict ac_position_map, list sequences, list found, dict not_found):
        # void getSequences(String &database_filename,
        # libcpp_map[ String, Size ] &ac_position_map,
        # libcpp_vector[ String ] &sequences,
        # libcpp_vector[ libcpp_pair[ String, Size ] ] &found,
        # libcpp_map[ String, Size ] &not_found) except + nogil 
        assert isinstance(database_filename, String), 'arg database_filename wrong type'
        assert isinstance(sequences, list) and all(isinstance(i, bytes) for i in sequences), 'arg sequences wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in sequences:
           v1.push_back(_String(<char *>item1))

        cdef libcpp_map[_String, size_t ] c_dict_ac_position
        for k,v in ac_position_map.iteritems():
            c_dict_ac_position[ _String(<char *>k) ] = v

        cdef libcpp_map[_String, size_t ] c_dict_not_found
        for k,v in not_found.iteritems():
            c_dict_not_found[ _String(<char *>k) ] = v

        cdef libcpp_vector[ libcpp_pair[_String, size_t] ] * v4 = new libcpp_vector[ libcpp_pair[_String, size_t] ]()
        cdef list item4
        cdef libcpp_pair[_String, size_t] * aPair
        for item4 in found:

           aPair = new libcpp_pair[ _String, size_t] ( _String(<char *>item4[0]), item4[1])
           v4.push_back( deref(aPair) )
           del aPair
 
        #
        # Call C++
        #
        self.inst.get().getSequences(deref(database_filename.inst.get()), c_dict_ac_position, deref(v1), deref(v4), c_dict_not_found)
        #
        # Get arguments back from C++
        #

        # libcpp_map[ String, Size ] &ac_position_map,
        replace_2 = dict()
        cdef libcpp_map[_String, size_t].iterator it_dict2 = c_dict_ac_position.begin()
        while it_dict2 != c_dict_ac_position.end():
            replace_2[ <libcpp_string>deref(it_dict2).first ] = deref(it_dict2).second
            inc(it_dict2)
        ac_position_map.update(replace_2)

        # libcpp_vector[ String ] &sequences,
        cdef replace = []
        cdef libcpp_vector[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.append(<char*>deref(it).c_str())
           inc(it)
        sequences[:] = replace

        # libcpp_vector[ libcpp_pair[ String, Size ] ] &found,
        replace_4 = list()
        cdef libcpp_vector[ libcpp_pair[_String, size_t] ].iterator it_list4 = deref(v4).begin()
        cdef libcpp_pair[_String, size_t] anotherPair
        while it_list4 != deref(v4).end():
            # get the pair back
            anotherPair = deref(it_list4)
            replace_4.append( [<libcpp_string>anotherPair.first, anotherPair.second] )
            inc(it_list4)
        found[:] = replace_4

        # libcpp_map[ String, Size ] &not_found
        replace_5 = dict()
        cdef libcpp_map[_String, size_t].iterator it_dict5 = c_dict_not_found.begin()
        while it_dict5 != c_dict_not_found.end():
            replace_2[ <libcpp_string>deref(it_dict5).first ] = deref(it_dict5).second
            inc(it_dict5)
        not_found.update(replace_2)


        del v1
        del v4 

cdef class StringView:
    """
    Cython implementation of _StringView

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1StringView.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef StringView rv = StringView.__new__(StringView)
       rv.inst = shared_ptr[_StringView](new _StringView(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef StringView rv = StringView.__new__(StringView)
       rv.inst = shared_ptr[_StringView](new _StringView(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_StringView](new _StringView())
    
    def _init_1(self, bytes in_0 ):
        """
        _init_1(self, in_0: bytes ) -> None
        """
        assert isinstance(in_0, bytes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_StringView](new _StringView((<libcpp_string>in_0)))
    
    def _init_2(self, StringView in_0 ):
        """
        _init_2(self, in_0: StringView ) -> None
        """
        assert isinstance(in_0, StringView), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_StringView](new _StringView((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: bytes ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: StringView ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], bytes)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], StringView)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def substr(self,  start ,  end ):
        """
        substr(self, start: int , end: int ) -> StringView
        """
        assert isinstance(start, int) and start >= 0, 'arg start wrong type'
        assert isinstance(end, int) and end >= 0, 'arg end wrong type'
    
    
        cdef _StringView * _r = new _StringView(self.inst.get().substr((<size_t>start), (<size_t>end)))
        cdef StringView py_result = StringView.__new__(StringView)
        py_result.inst = shared_ptr[_StringView](_r)
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def getString(self):
        """
        getString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getString()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (0,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, StringView):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef StringView other_casted = other
        cdef StringView self_casted = self
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class Unit:
    """
    Cython implementation of _Unit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Unit.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property accession:
        def __set__(self,  accession):
        
            self.inst.get().accession = deref((convString(accession)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().accession
            py_result = convOutputString(_r)
            return py_result
    
    property name:
        def __set__(self,  name):
        
            self.inst.get().name = deref((convString(name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().name
            py_result = convOutputString(_r)
            return py_result
    
    property cv_ref:
        def __set__(self,  cv_ref):
        
            self.inst.get().cv_ref = deref((convString(cv_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().cv_ref
            py_result = convOutputString(_r)
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Unit](new _Unit())
    
    def _init_1(self, Unit in_0 ):
        """
        _init_1(self, in_0: Unit ) -> None
        """
        assert isinstance(in_0, Unit), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Unit](new _Unit((deref(in_0.inst.get()))))
    
    def _init_2(self,  p_accession ,  p_name ,  p_cv_ref ):
        """
        _init_2(self, p_accession: Union[bytes, str, String] , p_name: Union[bytes, str, String] , p_cv_ref: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(p_accession, str) or isinstance(p_accession, bytes) or isinstance(p_accession, String)), 'arg p_accession wrong type'
        assert (isinstance(p_name, str) or isinstance(p_name, bytes) or isinstance(p_name, String)), 'arg p_name wrong type'
        assert (isinstance(p_cv_ref, str) or isinstance(p_cv_ref, bytes) or isinstance(p_cv_ref, String)), 'arg p_cv_ref wrong type'
    
    
    
        self.inst = shared_ptr[_Unit](new _Unit(deref((convString(p_accession)).get()), deref((convString(p_name)).get()), deref((convString(p_cv_ref)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Unit ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, p_accession: Union[bytes, str, String] , p_name: Union[bytes, str, String] , p_cv_ref: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Unit)):
             self._init_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Unit):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Unit other_casted = other
        cdef Unit self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class XMLHandler:
    """
    Cython implementation of _XMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1XMLHandler.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self,  filename ,  version ):
        """
        __init__(self, filename: Union[bytes, str, String] , version: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(version, str) or isinstance(version, bytes) or isinstance(version, String)), 'arg version wrong type'
    
    
        self.inst = shared_ptr[_XMLHandler](new _XMLHandler(deref((convString(filename)).get()), deref((convString(version)).get())))
    
    def reset(self):
        """
        reset(self) -> None
        """
        self.inst.get().reset()
    
    def error(self, int mode ,  msg ,  line ,  column ):
        """
        error(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None
        """
        assert mode in [0, 1], 'arg mode wrong type'
        assert (isinstance(msg, str) or isinstance(msg, bytes) or isinstance(msg, String)), 'arg msg wrong type'
        assert isinstance(line, int), 'arg line wrong type'
        assert isinstance(column, int), 'arg column wrong type'
    
    
    
    
        self.inst.get().error((<_ActionMode>mode), deref((convString(msg)).get()), (<unsigned int>line), (<unsigned int>column))
    
    def warning(self, int mode ,  msg ,  line ,  column ):
        """
        warning(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None
        """
        assert mode in [0, 1], 'arg mode wrong type'
        assert (isinstance(msg, str) or isinstance(msg, bytes) or isinstance(msg, String)), 'arg msg wrong type'
        assert isinstance(line, int), 'arg line wrong type'
        assert isinstance(column, int), 'arg column wrong type'
    
    
    
    
        self.inst.get().warning((<_ActionMode>mode), deref((convString(msg)).get()), (<unsigned int>line), (<unsigned int>column))
    ActionMode = __ActionMode 
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
