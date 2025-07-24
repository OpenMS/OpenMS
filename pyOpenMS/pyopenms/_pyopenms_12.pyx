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
def __static_CalibrationData_getMetaValues():
    """
    __static_CalibrationData_getMetaValues() -> List[bytes]
    """
    _r = _getMetaValues_CalibrationData()
    py_result = []
    cdef libcpp_vector[_String].iterator it__r = _r.begin()
    while it__r != _r.end():
       py_result.append(<char*>deref(it__r).c_str())
       inc(it__r)
    return py_result 

cdef class __LPWrapper_Type:
    None
    UNBOUNDED = 0
    LOWER_BOUND_ONLY = 1
    UPPER_BOUND_ONLY = 2
    DOUBLE_BOUNDED = 3
    FIXED = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __SOLVER:
    None
    SOLVER_GLPK = 0

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __Sense:
    None
    MIN = 0
    MAX = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __SolverStatus:
    None
    UNDEFINED = 0
    OPTIMAL = 1
    FEASIBLE = 2
    NO_FEASIBLE_SOL = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ValueType:
    None
    STRING_VALUE = 0
    INT_VALUE = 1
    DOUBLE_VALUE = 2
    STRING_LIST = 3
    INT_LIST = 4
    DOUBLE_LIST = 5
    EMPTY_VALUE = 6

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __VariableType:
    None
    CONTINUOUS = 0
    INTEGER = 1
    BINARY = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __WriteFormat:
    None
    FORMAT_LP = 0
    FORMAT_MPS = 1
    FORMAT_GLPK = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AbsoluteQuantitationMethodFile:
    """
    Cython implementation of _AbsoluteQuantitationMethodFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationMethodFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AbsoluteQuantitationMethodFile rv = AbsoluteQuantitationMethodFile.__new__(AbsoluteQuantitationMethodFile)
       rv.inst = shared_ptr[_AbsoluteQuantitationMethodFile](new _AbsoluteQuantitationMethodFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AbsoluteQuantitationMethodFile rv = AbsoluteQuantitationMethodFile.__new__(AbsoluteQuantitationMethodFile)
       rv.inst = shared_ptr[_AbsoluteQuantitationMethodFile](new _AbsoluteQuantitationMethodFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AbsoluteQuantitationMethodFile](new _AbsoluteQuantitationMethodFile())
    
    def _init_1(self, AbsoluteQuantitationMethodFile in_0 ):
        """
        _init_1(self, in_0: AbsoluteQuantitationMethodFile ) -> None
        """
        assert isinstance(in_0, AbsoluteQuantitationMethodFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AbsoluteQuantitationMethodFile](new _AbsoluteQuantitationMethodFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AbsoluteQuantitationMethodFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AbsoluteQuantitationMethodFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , list aqm_list ):
        """
        load(self, filename: Union[bytes, str, String] , aqm_list: List[AbsoluteQuantitationMethod] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(aqm_list, list) and all(isinstance(elemt_rec, AbsoluteQuantitationMethod) for elemt_rec in aqm_list), 'arg aqm_list wrong type'
    
        cdef libcpp_vector[_AbsoluteQuantitationMethod] * v1 = new libcpp_vector[_AbsoluteQuantitationMethod]()
        cdef AbsoluteQuantitationMethod item1
        for item1 in aqm_list:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().load(deref((convString(filename)).get()), deref(v1))
        cdef libcpp_vector[_AbsoluteQuantitationMethod].iterator it_aqm_list = v1.begin()
        replace_0 = []
        while it_aqm_list != v1.end():
            item1 = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
            item1.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(it_aqm_list)))
            replace_0.append(item1)
            inc(it_aqm_list)
        aqm_list[:] = replace_0
        del v1
    
    def store(self,  filename , list aqm_list ):
        """
        store(self, filename: Union[bytes, str, String] , aqm_list: List[AbsoluteQuantitationMethod] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(aqm_list, list) and all(isinstance(elemt_rec, AbsoluteQuantitationMethod) for elemt_rec in aqm_list), 'arg aqm_list wrong type'
    
        cdef libcpp_vector[_AbsoluteQuantitationMethod] * v1 = new libcpp_vector[_AbsoluteQuantitationMethod]()
        cdef AbsoluteQuantitationMethod item1
        for item1 in aqm_list:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().store(deref((convString(filename)).get()), deref(v1))
        cdef libcpp_vector[_AbsoluteQuantitationMethod].iterator it_aqm_list = v1.begin()
        replace_0 = []
        while it_aqm_list != v1.end():
            item1 = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
            item1.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(it_aqm_list)))
            replace_0.append(item1)
            inc(it_aqm_list)
        aqm_list[:] = replace_0
        del v1 

cdef class CalibrationData:
    """
    Cython implementation of _CalibrationData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CalibrationData.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CalibrationData rv = CalibrationData.__new__(CalibrationData)
       rv.inst = shared_ptr[_CalibrationData](new _CalibrationData(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CalibrationData rv = CalibrationData.__new__(CalibrationData)
       rv.inst = shared_ptr[_CalibrationData](new _CalibrationData(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CalibrationData](new _CalibrationData())
    
    def _init_1(self, CalibrationData in_0 ):
        """
        _init_1(self, in_0: CalibrationData ) -> None
        """
        assert isinstance(in_0, CalibrationData), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CalibrationData](new _CalibrationData((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CalibrationData ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CalibrationData)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMZ(self,  in_0 ):
        """
        getMZ(self, in_0: int ) -> float
        Retrieve the observed m/z of the i'th calibration point
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getMZ((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def getRT(self,  in_0 ):
        """
        getRT(self, in_0: int ) -> float
        Retrieve the observed RT of the i'th calibration point
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getRT((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def getIntensity(self,  in_0 ):
        """
        getIntensity(self, in_0: int ) -> float
        Retrieve the intensity of the i'th calibration point
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getIntensity((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        Number of calibration points
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        Returns `True` if there are no peaks
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        Remove all calibration points
        """
        self.inst.get().clear()
    
    def setUsePPM(self, bool in_0 ):
        """
        setUsePPM(self, in_0: bool ) -> None
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setUsePPM((<bool>in_0))
    
    def usePPM(self):
        """
        usePPM(self) -> bool
        Current error unit (ppm or Th)
        """
        cdef bool _r = self.inst.get().usePPM()
        py_result = <bool>_r
        return py_result
    
    def insertCalibrationPoint(self, double rt , double mz_obs , float intensity , double mz_ref , double weight ,  group ):
        """
        insertCalibrationPoint(self, rt: float , mz_obs: float , intensity: float , mz_ref: float , weight: float , group: int ) -> None
        """
        assert isinstance(rt, float), 'arg rt wrong type'
        assert isinstance(mz_obs, float), 'arg mz_obs wrong type'
        assert isinstance(intensity, float), 'arg intensity wrong type'
        assert isinstance(mz_ref, float), 'arg mz_ref wrong type'
        assert isinstance(weight, float), 'arg weight wrong type'
        assert isinstance(group, int), 'arg group wrong type'
    
    
    
    
    
    
        self.inst.get().insertCalibrationPoint((<double>rt), (<double>mz_obs), (<float>intensity), (<double>mz_ref), (<double>weight), (<int>group))
    
    def getNrOfGroups(self):
        """
        getNrOfGroups(self) -> int
        Number of peak groups (can be 0)
        """
        cdef size_t _r = self.inst.get().getNrOfGroups()
        py_result = <size_t>_r
        return py_result
    
    def getError(self,  in_0 ):
        """
        getError(self, in_0: int ) -> float
        Retrieve the error for i'th calibrant in either ppm or Th (depending on usePPM())
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getError((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def getRefMZ(self,  in_0 ):
        """
        getRefMZ(self, in_0: int ) -> float
        Retrieve the theoretical m/z of the i'th calibration point
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getRefMZ((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def getWeight(self,  in_0 ):
        """
        getWeight(self, in_0: int ) -> float
        Retrieve the weight of the i'th calibration point
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getWeight((<size_t>in_0))
        py_result = <double>_r
        return py_result
    
    def getGroup(self,  i ):
        """
        getGroup(self, i: int ) -> int
        Retrieve the group of the i'th calibration point
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef int _r = self.inst.get().getGroup((<size_t>i))
        py_result = <int>_r
        return py_result
    
    def median(self, double in_0 , double in_1 ):
        """
        median(self, in_0: float , in_1: float ) -> CalibrationData
        Compute the median in the given RT range for every peak group
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        cdef _CalibrationData * _r = new _CalibrationData(self.inst.get().median((<double>in_0), (<double>in_1)))
        cdef CalibrationData py_result = CalibrationData.__new__(CalibrationData)
        py_result.inst = shared_ptr[_CalibrationData](_r)
        return py_result
    
    def sortByRT(self):
        """
        sortByRT(self) -> None
        Sort calibration points by RT, to allow for valid RT chunking
        """
        self.inst.get().sortByRT()
    getMetaValues = __static_CalibrationData_getMetaValues 

cdef class DigestionEnzyme:
    """
    Cython implementation of _DigestionEnzyme

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzyme.html>`_

      Base class for digestion enzymes
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DigestionEnzyme rv = DigestionEnzyme.__new__(DigestionEnzyme)
       rv.inst = shared_ptr[_DigestionEnzyme](new _DigestionEnzyme(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DigestionEnzyme rv = DigestionEnzyme.__new__(DigestionEnzyme)
       rv.inst = shared_ptr[_DigestionEnzyme](new _DigestionEnzyme(deref(self.inst.get())))
       return rv
    
    def _init_0(self, DigestionEnzyme in_0 ):
        """
        _init_0(self, in_0: DigestionEnzyme ) -> None
        """
        assert isinstance(in_0, DigestionEnzyme), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DigestionEnzyme](new _DigestionEnzyme((deref(in_0.inst.get()))))
    
    def _init_1(self,  name ,  cleavage_regex , set synonyms ,  regex_description ):
        """
        _init_1(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
        assert isinstance(synonyms, set) and all(isinstance(i, bytes) for i in synonyms), 'arg synonyms wrong type'
        assert (isinstance(regex_description, str) or isinstance(regex_description, bytes) or isinstance(regex_description, String)), 'arg regex_description wrong type'
    
    
        cdef libcpp_set[_String] * v2 = new libcpp_set[_String]()
        cdef bytes item2
        for item2 in synonyms:
           v2.insert(_String(<char *>item2))
    
        self.inst = shared_ptr[_DigestionEnzyme](new _DigestionEnzyme(deref((convString(name)).get()), deref((convString(cleavage_regex)).get()), deref(v2), deref((convString(regex_description)).get())))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v2.begin()
        while it != v2.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        synonyms.clear()
        synonyms.update(replace)
        del v2
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DigestionEnzyme ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], DigestionEnzyme)):
             self._init_0(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], set) and all(isinstance(i, bytes) for i in args[2])) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the enzyme
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the enzyme
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setSynonyms(self, set synonyms ):
        """
        setSynonyms(self, synonyms: Set[bytes] ) -> None
        Sets the synonyms
        """
        assert isinstance(synonyms, set) and all(isinstance(i, bytes) for i in synonyms), 'arg synonyms wrong type'
        cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
        cdef bytes item0
        for item0 in synonyms:
           v0.insert(_String(<char *>item0))
        self.inst.get().setSynonyms(deref(v0))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v0.begin()
        while it != v0.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        synonyms.clear()
        synonyms.update(replace)
        del v0
    
    def addSynonym(self,  synonym ):
        """
        addSynonym(self, synonym: Union[bytes, str, String] ) -> None
        Adds a synonym
        """
        assert (isinstance(synonym, str) or isinstance(synonym, bytes) or isinstance(synonym, String)), 'arg synonym wrong type'
    
        self.inst.get().addSynonym(deref((convString(synonym)).get()))
    
    def getSynonyms(self):
        """
        getSynonyms(self) -> Set[bytes]
        Returns the synonyms
        """
        _r = self.inst.get().getSynonyms()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setRegEx(self,  cleavage_regex ):
        """
        setRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> None
        Sets the cleavage regex
        """
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
    
        self.inst.get().setRegEx(deref((convString(cleavage_regex)).get()))
    
    def getRegEx(self):
        """
        getRegEx(self) -> Union[bytes, str, String]
        Returns the cleavage regex
        """
        cdef _String _r = self.inst.get().getRegEx()
        py_result = convOutputString(_r)
        return py_result
    
    def setRegExDescription(self,  value ):
        """
        setRegExDescription(self, value: Union[bytes, str, String] ) -> None
        Sets the regex description
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setRegExDescription(deref((convString(value)).get()))
    
    def getRegExDescription(self):
        """
        getRegExDescription(self) -> Union[bytes, str, String]
        Returns the regex description
        """
        cdef _String _r = self.inst.get().getRegExDescription()
        py_result = convOutputString(_r)
        return py_result
    
    def setValueFromFile(self,  key ,  value ):
        """
        setValueFromFile(self, key: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> bool
        Sets the value of a member variable based on an entry from an input file
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
    
        cdef bool _r = self.inst.get().setValueFromFile(deref((convString(key)).get()), deref((convString(value)).get()))
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, DigestionEnzyme):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DigestionEnzyme other_casted = other
        cdef DigestionEnzyme self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class ElutionPeakDetection:
    """
    Cython implementation of _ElutionPeakDetection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElutionPeakDetection.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ElutionPeakDetection rv = ElutionPeakDetection.__new__(ElutionPeakDetection)
       rv.inst = shared_ptr[_ElutionPeakDetection](new _ElutionPeakDetection(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ElutionPeakDetection rv = ElutionPeakDetection.__new__(ElutionPeakDetection)
       rv.inst = shared_ptr[_ElutionPeakDetection](new _ElutionPeakDetection(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ElutionPeakDetection](new _ElutionPeakDetection())
    
    def _init_1(self, ElutionPeakDetection in_0 ):
        """
        _init_1(self, in_0: ElutionPeakDetection ) -> None
        """
        assert isinstance(in_0, ElutionPeakDetection), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ElutionPeakDetection](new _ElutionPeakDetection((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ElutionPeakDetection ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ElutionPeakDetection)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _detectPeaks_0(self, Kernel_MassTrace in_ , list out ):
        """
        _detectPeaks_0(self, in_: Kernel_MassTrace , out: List[Kernel_MassTrace] ) -> None
        """
        assert isinstance(in_, Kernel_MassTrace), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in out), 'arg out wrong type'
    
        cdef libcpp_vector[_Kernel_MassTrace] * v1 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item1
        for item1 in out:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().detectPeaks((deref(in_.inst.get())), deref(v1))
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_out = v1.begin()
        replace_0 = []
        while it_out != v1.end():
            item1 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item1.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_out)))
            replace_0.append(item1)
            inc(it_out)
        out[:] = replace_0
        del v1
    
    def _detectPeaks_1(self, list in_ , list out ):
        """
        _detectPeaks_1(self, in_: List[Kernel_MassTrace] , out: List[Kernel_MassTrace] ) -> None
        """
        assert isinstance(in_, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in in_), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in out), 'arg out wrong type'
        cdef libcpp_vector[_Kernel_MassTrace] * v0 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item0
        for item0 in in_:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_Kernel_MassTrace] * v1 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item1
        for item1 in out:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().detectPeaks(deref(v0), deref(v1))
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_out = v1.begin()
        replace_0 = []
        while it_out != v1.end():
            item1 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item1.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_out)))
            replace_0.append(item1)
            inc(it_out)
        out[:] = replace_0
        del v1
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_in_ = v0.begin()
        replace_0 = []
        while it_in_ != v0.end():
            item0 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item0.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_in_)))
            replace_0.append(item0)
            inc(it_in_)
        in_[:] = replace_0
        del v0
    
    def detectPeaks(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: detectPeaks(self, in_: Kernel_MassTrace , out: List[Kernel_MassTrace] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: detectPeaks(self, in_: List[Kernel_MassTrace] , out: List[Kernel_MassTrace] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], Kernel_MassTrace)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in args[1])):
            return self._detectPeaks_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in args[1])):
            return self._detectPeaks_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filterByPeakWidth(self, list in_ , list out ):
        """
        filterByPeakWidth(self, in_: List[Kernel_MassTrace] , out: List[Kernel_MassTrace] ) -> None
        """
        assert isinstance(in_, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in in_), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in out), 'arg out wrong type'
        cdef libcpp_vector[_Kernel_MassTrace] * v0 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item0
        for item0 in in_:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_Kernel_MassTrace] * v1 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item1
        for item1 in out:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().filterByPeakWidth(deref(v0), deref(v1))
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_out = v1.begin()
        replace_0 = []
        while it_out != v1.end():
            item1 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item1.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_out)))
            replace_0.append(item1)
            inc(it_out)
        out[:] = replace_0
        del v1
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_in_ = v0.begin()
        replace_0 = []
        while it_in_ != v0.end():
            item0 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item0.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_in_)))
            replace_0.append(item0)
            inc(it_in_)
        in_[:] = replace_0
        del v0
    
    def computeMassTraceNoise(self, Kernel_MassTrace in_0 ):
        """
        computeMassTraceNoise(self, in_0: Kernel_MassTrace ) -> float
        Compute noise level (as RMSE of the actual signal and the smoothed signal)
        """
        assert isinstance(in_0, Kernel_MassTrace), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().computeMassTraceNoise((deref(in_0.inst.get())))
        py_result = <double>_r
        return py_result
    
    def computeMassTraceSNR(self, Kernel_MassTrace in_0 ):
        """
        computeMassTraceSNR(self, in_0: Kernel_MassTrace ) -> float
        Compute the signal to noise ratio (estimated by computeMassTraceNoise)
        """
        assert isinstance(in_0, Kernel_MassTrace), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().computeMassTraceSNR((deref(in_0.inst.get())))
        py_result = <double>_r
        return py_result
    
    def computeApexSNR(self, Kernel_MassTrace in_0 ):
        """
        computeApexSNR(self, in_0: Kernel_MassTrace ) -> float
        Compute the signal to noise ratio at the apex (estimated by computeMassTraceNoise)
        """
        assert isinstance(in_0, Kernel_MassTrace), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().computeApexSNR((deref(in_0.inst.get())))
        py_result = <double>_r
        return py_result
    
    def findLocalExtrema(self, Kernel_MassTrace in_0 ,  in_1 , list in_2 , list in_3 ):
        """
        findLocalExtrema(self, in_0: Kernel_MassTrace , in_1: int , in_2: List[int] , in_3: List[int] ) -> None
        """
        assert isinstance(in_0, Kernel_MassTrace), 'arg in_0 wrong type'
        assert isinstance(in_1, int) and in_1 >= 0, 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in in_2), 'arg in_2 wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in in_3), 'arg in_3 wrong type'
    
    
        cdef libcpp_vector[size_t] v2 = in_2
        cdef libcpp_vector[size_t] v3 = in_3
        self.inst.get().findLocalExtrema((deref(in_0.inst.get())), (<size_t &>in_1), v2, v3)
        in_3[:] = v3
        in_2[:] = v2
    
    def smoothData(self, Kernel_MassTrace mt ,  win_size ):
        """
        smoothData(self, mt: Kernel_MassTrace , win_size: int ) -> None
        """
        assert isinstance(mt, Kernel_MassTrace), 'arg mt wrong type'
        assert isinstance(win_size, int), 'arg win_size wrong type'
    
    
        self.inst.get().smoothData((deref(mt.inst.get())), (<int>win_size))
    
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

cdef class FeatureXMLFile:
    """
    Cython implementation of _FeatureXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        This class provides Input/Output functionality for feature maps
        """
        self.inst = shared_ptr[_FeatureXMLFile](new _FeatureXMLFile())
    
    def load(self,  in_0 , FeatureMap in_1 ):
        """
        load(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None
        Loads the file with name `filename` into `map` and calls updateRanges()
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
    
    
        self.inst.get().load(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def store(self,  in_0 , FeatureMap in_1 ):
        """
        store(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None
        Stores the map `feature_map` in file with name `filename`
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> FeatureFileOptions
        Access to the options for loading/storing
        """
        cdef _FeatureFileOptions * _r = new _FeatureFileOptions(self.inst.get().getOptions())
        cdef FeatureFileOptions py_result = FeatureFileOptions.__new__(FeatureFileOptions)
        py_result.inst = shared_ptr[_FeatureFileOptions](_r)
        return py_result
    
    def setOptions(self, FeatureFileOptions in_0 ):
        """
        setOptions(self, in_0: FeatureFileOptions ) -> None
        Setter for options for loading/storing
        """
        assert isinstance(in_0, FeatureFileOptions), 'arg in_0 wrong type'
    
        self.inst.get().setOptions((deref(in_0.inst.get())))
    
    def loadSize(self,  path ):
        """
        loadSize(self, path: Union[bytes, str, String] ) -> int
        """
        assert (isinstance(path, str) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'
    
        cdef size_t _r = self.inst.get().loadSize(deref((convString(path)).get()))
        py_result = <size_t>_r
        return py_result 

cdef class GNPSMetaValueFile:
    """
    Cython implementation of _GNPSMetaValueFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSMetaValueFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GNPSMetaValueFile rv = GNPSMetaValueFile.__new__(GNPSMetaValueFile)
       rv.inst = shared_ptr[_GNPSMetaValueFile](new _GNPSMetaValueFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GNPSMetaValueFile rv = GNPSMetaValueFile.__new__(GNPSMetaValueFile)
       rv.inst = shared_ptr[_GNPSMetaValueFile](new _GNPSMetaValueFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_GNPSMetaValueFile](new _GNPSMetaValueFile())
    
    def _init_1(self, GNPSMetaValueFile in_0 ):
        """
        _init_1(self, in_0: GNPSMetaValueFile ) -> None
        """
        assert isinstance(in_0, GNPSMetaValueFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GNPSMetaValueFile](new _GNPSMetaValueFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GNPSMetaValueFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], GNPSMetaValueFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self, ConsensusMap consensus_map ,  output_file ):
        """
        store(self, consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None
        Write meta value table (tsv file) from a list of mzML files. Required for GNPS FBMN.
        
        This will produce the minimal required meta values and can be extended manually.
        
        :param consensus_map: Input ConsensusMap from which the input mzML files will be determined.
        :param output_file: Output file path for the meta value table.
        """
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
        assert (isinstance(output_file, str) or isinstance(output_file, bytes) or isinstance(output_file, String)), 'arg output_file wrong type'
    
    
        self.inst.get().store((deref(consensus_map.inst.get())), deref((convString(output_file)).get())) 

cdef class HyperScore:
    """
    Cython implementation of _HyperScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1HyperScore.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef HyperScore rv = HyperScore.__new__(HyperScore)
       rv.inst = shared_ptr[_HyperScore](new _HyperScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef HyperScore rv = HyperScore.__new__(HyperScore)
       rv.inst = shared_ptr[_HyperScore](new _HyperScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        An implementation of the X!Tandem HyperScore PSM scoring function
        """
        self.inst = shared_ptr[_HyperScore](new _HyperScore())
    
    def _init_1(self, HyperScore in_0 ):
        """
        _init_1(self, in_0: HyperScore ) -> None
        """
        assert isinstance(in_0, HyperScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_HyperScore](new _HyperScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        An implementation of the X!Tandem HyperScore PSM scoring function

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: HyperScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], HyperScore)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , MSSpectrum exp_spectrum , MSSpectrum theo_spectrum ):
        """
        compute(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , theo_spectrum: MSSpectrum ) -> float
        Compute the (ln transformed) X!Tandem HyperScore\n
        
        1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
        2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
        
        
        :note: Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
        :param fragment_mass_tolerance: Mass tolerance applied left and right of the theoretical spectrum peak position
        :param fragment_mass_tolerance_unit_ppm: Unit of the mass tolerance is: Thomson if false, ppm if true
        :param exp_spectrum: Measured spectrum
        :param theo_spectrum: Theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator
        """
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(exp_spectrum, MSSpectrum), 'arg exp_spectrum wrong type'
        assert isinstance(theo_spectrum, MSSpectrum), 'arg theo_spectrum wrong type'
    
    
    
    
        cdef double _r = self.inst.get().compute((<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), (deref(exp_spectrum.inst.get())), (deref(theo_spectrum.inst.get())))
        py_result = <double>_r
        return py_result 

cdef class IMSWeights:
    """
    Cython implementation of _IMSWeights

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::Weights_1_1IMSWeights.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IMSWeights](new _IMSWeights())
    
    def _init_1(self, IMSWeights in_0 ):
        """
        _init_1(self, in_0: IMSWeights ) -> None
        """
        assert isinstance(in_0, IMSWeights), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IMSWeights](new _IMSWeights((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IMSWeights ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSWeights)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        Gets size of a set of weights
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def getWeight(self,  i ):
        """
        getWeight(self, i: int ) -> int
        Gets a scaled integer weight by index
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef unsigned int _r = self.inst.get().getWeight((<int>i))
        py_result = <unsigned int>_r
        return py_result
    
    def setPrecision(self, double precision ):
        """
        setPrecision(self, precision: float ) -> None
        Sets a new precision to scale double values to integer
        """
        assert isinstance(precision, float), 'arg precision wrong type'
    
        self.inst.get().setPrecision((<double>precision))
    
    def getPrecision(self):
        """
        getPrecision(self) -> float
        Gets precision.
        """
        cdef double _r = self.inst.get().getPrecision()
        py_result = <double>_r
        return py_result
    
    def back(self):
        """
        back(self) -> int
        Gets a last weight
        """
        cdef unsigned int _r = self.inst.get().back()
        py_result = <unsigned int>_r
        return py_result
    
    def getAlphabetMass(self,  i ):
        """
        getAlphabetMass(self, i: int ) -> float
        Gets an original (double) alphabet mass by index
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef double _r = self.inst.get().getAlphabetMass((<int>i))
        py_result = <double>_r
        return py_result
    
    def getParentMass(self, list decomposition ):
        """
        getParentMass(self, decomposition: List[int] ) -> float
        Returns a parent mass for a given `decomposition`
        """
        assert isinstance(decomposition, list) and all(isinstance(elemt_rec, int) for elemt_rec in decomposition), 'arg decomposition wrong type'
        cdef libcpp_vector[unsigned int] v0 = decomposition
        cdef double _r = self.inst.get().getParentMass(v0)
        decomposition[:] = v0
        py_result = <double>_r
        return py_result
    
    def swap(self,  index1 ,  index2 ):
        """
        swap(self, index1: int , index2: int ) -> None
        Exchanges weight and mass at index1 with weight and mass at index2
        """
        assert isinstance(index1, int), 'arg index1 wrong type'
        assert isinstance(index2, int), 'arg index2 wrong type'
    
    
        self.inst.get().swap((<int>index1), (<int>index2))
    
    def divideByGCD(self):
        """
        divideByGCD(self) -> bool
        Divides the integer weights by their gcd. The precision is also adjusted
        """
        cdef bool _r = self.inst.get().divideByGCD()
        py_result = <bool>_r
        return py_result
    
    def getMinRoundingError(self):
        """
        getMinRoundingError(self) -> float
        """
        cdef double _r = self.inst.get().getMinRoundingError()
        py_result = <double>_r
        return py_result
    
    def getMaxRoundingError(self):
        """
        getMaxRoundingError(self) -> float
        """
        cdef double _r = self.inst.get().getMaxRoundingError()
        py_result = <double>_r
        return py_result 

cdef class IsobaricNormalizer:
    """
    Cython implementation of _IsobaricNormalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricNormalizer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsobaricNormalizer rv = IsobaricNormalizer.__new__(IsobaricNormalizer)
       rv.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricNormalizer rv = IsobaricNormalizer.__new__(IsobaricNormalizer)
       rv.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer(deref(self.inst.get())))
       return rv
    
    def _init_0(self, IsobaricNormalizer in_0 ):
        """
        _init_0(self, in_0: IsobaricNormalizer ) -> None
        """
        assert isinstance(in_0, IsobaricNormalizer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer((deref(in_0.inst.get()))))
    
    def _init_1(self, ItraqFourPlexQuantitationMethod quant_method ):
        """
        _init_1(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqFourPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer((quant_method.inst.get())))
    
    def _init_2(self, ItraqEightPlexQuantitationMethod quant_method ):
        """
        _init_2(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqEightPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer((quant_method.inst.get())))
    
    def _init_3(self, TMTSixPlexQuantitationMethod quant_method ):
        """
        _init_3(self, quant_method: TMTSixPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTSixPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer((quant_method.inst.get())))
    
    def _init_4(self, TMTTenPlexQuantitationMethod quant_method ):
        """
        _init_4(self, quant_method: TMTTenPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTTenPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricNormalizer](new _IsobaricNormalizer((quant_method.inst.get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricNormalizer ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: TMTSixPlexQuantitationMethod ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, quant_method: TMTTenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], IsobaricNormalizer)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqFourPlexQuantitationMethod)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqEightPlexQuantitationMethod)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTSixPlexQuantitationMethod)):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTTenPlexQuantitationMethod)):
             self._init_4(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def normalize(self, ConsensusMap consensus_map ):
        """
        normalize(self, consensus_map: ConsensusMap ) -> None
        """
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
    
        self.inst.get().normalize((deref(consensus_map.inst.get()))) 

cdef class JavaInfo:
    """
    Cython implementation of _JavaInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1JavaInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef JavaInfo rv = JavaInfo.__new__(JavaInfo)
       rv.inst = shared_ptr[_JavaInfo](new _JavaInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef JavaInfo rv = JavaInfo.__new__(JavaInfo)
       rv.inst = shared_ptr[_JavaInfo](new _JavaInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Detect Java and retrieve information
        """
        self.inst = shared_ptr[_JavaInfo](new _JavaInfo())
    
    def _init_1(self, JavaInfo in_0 ):
        """
        _init_1(self, in_0: JavaInfo ) -> None
        """
        assert isinstance(in_0, JavaInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_JavaInfo](new _JavaInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Detect Java and retrieve information

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: JavaInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], JavaInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def canRun(self,  java_executable ):
        """
        canRun(self, java_executable: Union[bytes, str, String] ) -> bool
        Determine if Java is installed and reachable\n
        
        The call fails if either Java is not installed or if a relative location is given and Java is not on the search PATH
        
        
        :param java_executable: Path to Java executable. Can be absolute, relative or just a filename
        :return: Returns false if Java executable can not be called; true if Java executable can be executed
        """
        assert (isinstance(java_executable, str) or isinstance(java_executable, bytes) or isinstance(java_executable, String)), 'arg java_executable wrong type'
    
        cdef bool _r = self.inst.get().canRun(deref((convString(java_executable)).get()))
        py_result = <bool>_r
        return py_result 

cdef class LPWrapper:
    """
    Cython implementation of _LPWrapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LPWrapper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_LPWrapper](new _LPWrapper())
    
    def _addRow_0(self, list row_indices , list row_values ,  name ):
        """
        _addRow_0(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] ) -> int
        Adds a row to the LP matrix, returns index
        """
        assert isinstance(row_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in row_indices), 'arg row_indices wrong type'
        assert isinstance(row_values, list) and all(isinstance(elemt_rec, float) for elemt_rec in row_values), 'arg row_values wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        cdef libcpp_vector[int] v0 = row_indices
        cdef libcpp_vector[double] v1 = row_values
    
        cdef int _r = self.inst.get().addRow(v0, v1, deref((convString(name)).get()))
        
        
        py_result = <int>_r
        return py_result
    
    def _addRow_1(self, list row_indices , list row_values ,  name , double lower_bound , double upper_bound , int type_ ):
        """
        _addRow_1(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int
        Adds a row with boundaries to the LP matrix, returns index
        """
        assert isinstance(row_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in row_indices), 'arg row_indices wrong type'
        assert isinstance(row_values, list) and all(isinstance(elemt_rec, float) for elemt_rec in row_values), 'arg row_values wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(lower_bound, float), 'arg lower_bound wrong type'
        assert isinstance(upper_bound, float), 'arg upper_bound wrong type'
        assert type_ in [0, 1, 2, 3, 4], 'arg type_ wrong type'
        cdef libcpp_vector[int] v0 = row_indices
        cdef libcpp_vector[double] v1 = row_values
    
    
    
    
        cdef int _r = self.inst.get().addRow(v0, v1, deref((convString(name)).get()), (<double>lower_bound), (<double>upper_bound), (<_LPWrapper_Type>type_))
        row_values[:] = v1
        row_indices[:] = v0
        py_result = <int>_r
        return py_result
    
    def addRow(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: addRow(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] ) -> int
          :noindex:
        
        Adds a row to the LP matrix, returns index

        
        .. rubric:: Overload:
        .. py:function:: addRow(self, row_indices: List[int] , row_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int
          :noindex:
        
        Adds a row with boundaries to the LP matrix, returns index
    
        """
        if (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
            return self._addRow_0(*args)
        elif (len(args)==6) and (isinstance(args[0], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], float)) and (isinstance(args[4], float)) and (args[5] in [0, 1, 2, 3, 4]):
            return self._addRow_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _addColumn_0(self):
        """
        _addColumn_0(self) -> int
        Adds an empty column to the LP matrix, returns index
        """
        cdef int _r = self.inst.get().addColumn()
        py_result = <int>_r
        return py_result
    
    def _addColumn_1(self, list column_indices , list column_values ,  name ):
        """
        _addColumn_1(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] ) -> int
        Adds a column to the LP matrix, returns index
        """
        assert isinstance(column_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in column_indices), 'arg column_indices wrong type'
        assert isinstance(column_values, list) and all(isinstance(elemt_rec, float) for elemt_rec in column_values), 'arg column_values wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        cdef libcpp_vector[int] v0 = column_indices
        cdef libcpp_vector[double] v1 = column_values
    
        cdef int _r = self.inst.get().addColumn(v0, v1, deref((convString(name)).get()))
        
        
        py_result = <int>_r
        return py_result
    
    def _addColumn_2(self, list column_indices , list column_values ,  name , double lower_bound , double upper_bound , int type_ ):
        """
        _addColumn_2(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int
        Adds a column with boundaries to the LP matrix, returns index
        """
        assert isinstance(column_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in column_indices), 'arg column_indices wrong type'
        assert isinstance(column_values, list) and all(isinstance(elemt_rec, float) for elemt_rec in column_values), 'arg column_values wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(lower_bound, float), 'arg lower_bound wrong type'
        assert isinstance(upper_bound, float), 'arg upper_bound wrong type'
        assert type_ in [0, 1, 2, 3, 4], 'arg type_ wrong type'
        cdef libcpp_vector[int] v0 = column_indices
        cdef libcpp_vector[double] v1 = column_values
    
    
    
    
        cdef int _r = self.inst.get().addColumn(v0, v1, deref((convString(name)).get()), (<double>lower_bound), (<double>upper_bound), (<_LPWrapper_Type>type_))
        column_values[:] = v1
        column_indices[:] = v0
        py_result = <int>_r
        return py_result
    
    def addColumn(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: addColumn(self, ) -> int
          :noindex:
        
        Adds an empty column to the LP matrix, returns index

        
        .. rubric:: Overload:
        .. py:function:: addColumn(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] ) -> int
          :noindex:
        
        Adds a column to the LP matrix, returns index

        
        .. rubric:: Overload:
        .. py:function:: addColumn(self, column_indices: List[int] , column_values: List[float] , name: Union[bytes, str, String] , lower_bound: float , upper_bound: float , type_: int ) -> int
          :noindex:
        
        Adds a column with boundaries to the LP matrix, returns index
    
        """
        if not args:
            return self._addColumn_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
            return self._addColumn_1(*args)
        elif (len(args)==6) and (isinstance(args[0], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], float)) and (isinstance(args[4], float)) and (args[5] in [0, 1, 2, 3, 4]):
            return self._addColumn_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def deleteRow(self,  index ):
        """
        deleteRow(self, index: int ) -> None
        Delete index-th row
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        self.inst.get().deleteRow((<int>index))
    
    def setColumnName(self,  index ,  name ):
        """
        setColumnName(self, index: int , name: Union[bytes, str, String] ) -> None
        Sets name of the index-th column
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
    
        self.inst.get().setColumnName((<int>index), deref((convString(name)).get()))
    
    def getColumnName(self,  index ):
        """
        getColumnName(self, index: int ) -> Union[bytes, str, String]
        Returns name of the index-th column
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _String _r = self.inst.get().getColumnName((<int>index))
        py_result = convOutputString(_r)
        return py_result
    
    def getRowName(self,  index ):
        """
        getRowName(self, index: int ) -> Union[bytes, str, String]
        Sets name of the index-th row
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _String _r = self.inst.get().getRowName((<int>index))
        py_result = convOutputString(_r)
        return py_result
    
    def getRowIndex(self,  name ):
        """
        getRowIndex(self, name: Union[bytes, str, String] ) -> int
        Returns index of the row with name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef int _r = self.inst.get().getRowIndex(deref((convString(name)).get()))
        py_result = <int>_r
        return py_result
    
    def getColumnIndex(self,  name ):
        """
        getColumnIndex(self, name: Union[bytes, str, String] ) -> int
        Returns index of the column with name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef int _r = self.inst.get().getColumnIndex(deref((convString(name)).get()))
        py_result = <int>_r
        return py_result
    
    def getColumnUpperBound(self,  index ):
        """
        getColumnUpperBound(self, index: int ) -> float
        Returns column's upper bound
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getColumnUpperBound((<int>index))
        py_result = <double>_r
        return py_result
    
    def getColumnLowerBound(self,  index ):
        """
        getColumnLowerBound(self, index: int ) -> float
        Returns column's lower bound
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getColumnLowerBound((<int>index))
        py_result = <double>_r
        return py_result
    
    def getRowUpperBound(self,  index ):
        """
        getRowUpperBound(self, index: int ) -> float
        Returns row's upper bound
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getRowUpperBound((<int>index))
        py_result = <double>_r
        return py_result
    
    def getRowLowerBound(self,  index ):
        """
        getRowLowerBound(self, index: int ) -> float
        Returns row's lower bound
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getRowLowerBound((<int>index))
        py_result = <double>_r
        return py_result
    
    def setRowName(self,  index ,  name ):
        """
        setRowName(self, index: int , name: Union[bytes, str, String] ) -> None
        Sets name of the index-th row
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
    
        self.inst.get().setRowName((<int>index), deref((convString(name)).get()))
    
    def setColumnBounds(self,  index , double lower_bound , double upper_bound , int type_ ):
        """
        setColumnBounds(self, index: int , lower_bound: float , upper_bound: float , type_: int ) -> None
        Sets column bounds
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert isinstance(lower_bound, float), 'arg lower_bound wrong type'
        assert isinstance(upper_bound, float), 'arg upper_bound wrong type'
        assert type_ in [0, 1, 2, 3, 4], 'arg type_ wrong type'
    
    
    
    
        self.inst.get().setColumnBounds((<int>index), (<double>lower_bound), (<double>upper_bound), (<_LPWrapper_Type>type_))
    
    def setRowBounds(self,  index , double lower_bound , double upper_bound , int type_ ):
        """
        setRowBounds(self, index: int , lower_bound: float , upper_bound: float , type_: int ) -> None
        Sets row bounds
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert isinstance(lower_bound, float), 'arg lower_bound wrong type'
        assert isinstance(upper_bound, float), 'arg upper_bound wrong type'
        assert type_ in [0, 1, 2, 3, 4], 'arg type_ wrong type'
    
    
    
    
        self.inst.get().setRowBounds((<int>index), (<double>lower_bound), (<double>upper_bound), (<_LPWrapper_Type>type_))
    
    def setColumnType(self,  index , int type_ ):
        """
        setColumnType(self, index: int , type_: int ) -> None
        Sets column/variable type.
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert type_ in [0, 1, 2], 'arg type_ wrong type'
    
    
        self.inst.get().setColumnType((<int>index), (<_VariableType>type_))
    
    def getColumnType(self,  index ):
        """
        getColumnType(self, index: int ) -> int
        Returns column/variable type.
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _VariableType _r = self.inst.get().getColumnType((<int>index))
        py_result = <int>_r
        return py_result
    
    def setObjective(self,  index , double obj_value ):
        """
        setObjective(self, index: int , obj_value: float ) -> None
        Sets objective value for column with index
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert isinstance(obj_value, float), 'arg obj_value wrong type'
    
    
        self.inst.get().setObjective((<int>index), (<double>obj_value))
    
    def getObjective(self,  index ):
        """
        getObjective(self, index: int ) -> float
        Returns objective value for column with index
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getObjective((<int>index))
        py_result = <double>_r
        return py_result
    
    def setObjectiveSense(self, int sense ):
        """
        setObjectiveSense(self, sense: int ) -> None
        Sets objective direction
        """
        assert sense in [0, 1], 'arg sense wrong type'
    
        self.inst.get().setObjectiveSense((<_Sense>sense))
    
    def getObjectiveSense(self):
        """
        getObjectiveSense(self) -> int
        Returns objective sense
        """
        cdef _Sense _r = self.inst.get().getObjectiveSense()
        py_result = <int>_r
        return py_result
    
    def getNumberOfColumns(self):
        """
        getNumberOfColumns(self) -> int
        Returns number of columns
        """
        cdef int _r = self.inst.get().getNumberOfColumns()
        py_result = <int>_r
        return py_result
    
    def getNumberOfRows(self):
        """
        getNumberOfRows(self) -> int
        Returns number of rows
        """
        cdef int _r = self.inst.get().getNumberOfRows()
        py_result = <int>_r
        return py_result
    
    def setElement(self,  row_index ,  column_index , double value ):
        """
        setElement(self, row_index: int , column_index: int , value: float ) -> None
        Sets the element
        """
        assert isinstance(row_index, int), 'arg row_index wrong type'
        assert isinstance(column_index, int), 'arg column_index wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
    
        self.inst.get().setElement((<int>row_index), (<int>column_index), (<double>value))
    
    def getElement(self,  row_index ,  column_index ):
        """
        getElement(self, row_index: int , column_index: int ) -> float
        Returns the element
        """
        assert isinstance(row_index, int), 'arg row_index wrong type'
        assert isinstance(column_index, int), 'arg column_index wrong type'
    
    
        cdef double _r = self.inst.get().getElement((<int>row_index), (<int>column_index))
        py_result = <double>_r
        return py_result
    
    def readProblem(self,  filename ,  format_ ):
        """
        readProblem(self, filename: Union[bytes, str, String] , format_: Union[bytes, str, String] ) -> None
        Read LP from file
        
        
        :param filename: Filename where to store the LP problem
        :param format: LP, MPS or GLPK
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(format_, str) or isinstance(format_, bytes) or isinstance(format_, String)), 'arg format_ wrong type'
    
    
        self.inst.get().readProblem(deref((convString(filename)).get()), deref((convString(format_)).get()))
    
    def writeProblem(self,  filename , int format_ ):
        """
        writeProblem(self, filename: Union[bytes, str, String] , format_: int ) -> None
        Write LP formulation to a file
        
        
        :param filename: Output filename, if the filename ends with '.gz' it will be compressed
        :param format: MPS-format is supported by GLPK and COIN-OR; LP and GLPK-formats only by GLPK
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert format_ in [0, 1, 2], 'arg format_ wrong type'
    
    
        self.inst.get().writeProblem(deref((convString(filename)).get()), (<_WriteFormat>format_))
    
    def solve(self, SolverParam solver_param ,  verbose_level ):
        """
        solve(self, solver_param: SolverParam , verbose_level: int ) -> int
        Solve problems, parameters like enabled heuristics can be given via solver_param\n
        
        The verbose level (0,1,2) determines if the solver prints status messages and internals
        
        
        :param solver_param: Parameters of the solver introduced by SolverParam
        :param verbose_level: Sets verbose level
        :return: solver dependent
        """
        assert isinstance(solver_param, SolverParam), 'arg solver_param wrong type'
        assert isinstance(verbose_level, int) and verbose_level >= 0, 'arg verbose_level wrong type'
    
    
        cdef int _r = self.inst.get().solve((deref(solver_param.inst.get())), (<size_t>verbose_level))
        py_result = <int>_r
        return py_result
    
    def getStatus(self):
        """
        getStatus(self) -> int
        Returns solution status
        
        
        :return: status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution
        """
        cdef _SolverStatus _r = self.inst.get().getStatus()
        py_result = <int>_r
        return py_result
    
    def getObjectiveValue(self):
        """
        getObjectiveValue(self) -> float
        """
        cdef double _r = self.inst.get().getObjectiveValue()
        py_result = <double>_r
        return py_result
    
    def getColumnValue(self,  index ):
        """
        getColumnValue(self, index: int ) -> float
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getColumnValue((<int>index))
        py_result = <double>_r
        return py_result
    
    def getNumberOfNonZeroEntriesInRow(self,  idx ):
        """
        getNumberOfNonZeroEntriesInRow(self, idx: int ) -> int
        """
        assert isinstance(idx, int), 'arg idx wrong type'
    
        cdef int _r = self.inst.get().getNumberOfNonZeroEntriesInRow((<int>idx))
        py_result = <int>_r
        return py_result
    
    def getMatrixRow(self,  idx , list indexes ):
        """
        getMatrixRow(self, idx: int , indexes: List[int] ) -> None
        """
        assert isinstance(idx, int), 'arg idx wrong type'
        assert isinstance(indexes, list) and all(isinstance(elemt_rec, int) for elemt_rec in indexes), 'arg indexes wrong type'
    
        cdef libcpp_vector[int] v1 = indexes
        self.inst.get().getMatrixRow((<int>idx), v1)
        indexes[:] = v1
    
    def getSolver(self):
        """
        getSolver(self) -> int
        Returns currently active solver
        """
        cdef _SOLVER _r = self.inst.get().getSolver()
        py_result = <int>_r
        return py_result
    LPWrapper_Type = __LPWrapper_Type
    SOLVER = __SOLVER
    Sense = __Sense
    SolverStatus = __SolverStatus
    VariableType = __VariableType
    WriteFormat = __WriteFormat 

cdef class OpenSwathDataAccessHelper:
    """
    Cython implementation of _OpenSwathDataAccessHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathDataAccessHelper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenSwathDataAccessHelper rv = OpenSwathDataAccessHelper.__new__(OpenSwathDataAccessHelper)
       rv.inst = shared_ptr[_OpenSwathDataAccessHelper](new _OpenSwathDataAccessHelper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenSwathDataAccessHelper rv = OpenSwathDataAccessHelper.__new__(OpenSwathDataAccessHelper)
       rv.inst = shared_ptr[_OpenSwathDataAccessHelper](new _OpenSwathDataAccessHelper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenSwathDataAccessHelper](new _OpenSwathDataAccessHelper())
    
    def _init_1(self, OpenSwathDataAccessHelper in_0 ):
        """
        _init_1(self, in_0: OpenSwathDataAccessHelper ) -> None
        """
        assert isinstance(in_0, OpenSwathDataAccessHelper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenSwathDataAccessHelper](new _OpenSwathDataAccessHelper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenSwathDataAccessHelper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenSwathDataAccessHelper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def convertToOpenMSSpectrum(self, OSSpectrum sptr , MSSpectrum spectrum ):
        """
        convertToOpenMSSpectrum(self, sptr: OSSpectrum , spectrum: MSSpectrum ) -> None
        Converts a SpectrumPtr to an OpenMS Spectrum
        """
        assert isinstance(sptr, OSSpectrum), 'arg sptr wrong type'
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        cdef shared_ptr[_OSSpectrum] input_sptr = sptr.inst
    
        self.inst.get().convertToOpenMSSpectrum(input_sptr, (deref(spectrum.inst.get())))
    
    def convertToOpenMSChromatogram(self, OSChromatogram cptr , MSChromatogram chromatogram ):
        """
        convertToOpenMSChromatogram(self, cptr: OSChromatogram , chromatogram: MSChromatogram ) -> None
        Converts a ChromatogramPtr to an OpenMS Chromatogram
        """
        assert isinstance(cptr, OSChromatogram), 'arg cptr wrong type'
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        cdef shared_ptr[_OSChromatogram] input_cptr = cptr.inst
    
        self.inst.get().convertToOpenMSChromatogram(input_cptr, (deref(chromatogram.inst.get())))
    
    def convertToOpenMSChromatogramFilter(self, MSChromatogram chromatogram , OSChromatogram cptr , double rt_min , double rt_max ):
        """
        convertToOpenMSChromatogramFilter(self, chromatogram: MSChromatogram , cptr: OSChromatogram , rt_min: float , rt_max: float ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert isinstance(cptr, OSChromatogram), 'arg cptr wrong type'
        assert isinstance(rt_min, float), 'arg rt_min wrong type'
        assert isinstance(rt_max, float), 'arg rt_max wrong type'
    
        cdef shared_ptr[_OSChromatogram] input_cptr = cptr.inst
    
    
        self.inst.get().convertToOpenMSChromatogramFilter((deref(chromatogram.inst.get())), input_cptr, (<double>rt_min), (<double>rt_max))
    
    def convertTargetedExp(self, TargetedExperiment transition_exp_ , LightTargetedExperiment transition_exp ):
        """
        convertTargetedExp(self, transition_exp_: TargetedExperiment , transition_exp: LightTargetedExperiment ) -> None
        Converts from the OpenMS TargetedExperiment to the OpenMs LightTargetedExperiment
        """
        assert isinstance(transition_exp_, TargetedExperiment), 'arg transition_exp_ wrong type'
        assert isinstance(transition_exp, LightTargetedExperiment), 'arg transition_exp wrong type'
    
    
        self.inst.get().convertTargetedExp((deref(transition_exp_.inst.get())), (deref(transition_exp.inst.get())))
    
    def convertPeptideToAASequence(self, LightCompound peptide , AASequence aa_sequence ):
        """
        convertPeptideToAASequence(self, peptide: LightCompound , aa_sequence: AASequence ) -> None
        Converts from the LightCompound to an OpenMS AASequence (with correct modifications)
        """
        assert isinstance(peptide, LightCompound), 'arg peptide wrong type'
        assert isinstance(aa_sequence, AASequence), 'arg aa_sequence wrong type'
    
    
        self.inst.get().convertPeptideToAASequence((deref(peptide.inst.get())), (deref(aa_sequence.inst.get())))
    
    def convertTargetedPeptide(self, Peptide pep , LightCompound p ):
        """
        convertTargetedPeptide(self, pep: Peptide , p: LightCompound ) -> None
        Converts from the OpenMS TargetedExperiment Peptide to the LightTargetedExperiment Peptide
        """
        assert isinstance(pep, Peptide), 'arg pep wrong type'
        assert isinstance(p, LightCompound), 'arg p wrong type'
    
    
        self.inst.get().convertTargetedPeptide((deref(pep.inst.get())), (deref(p.inst.get())))
    
    def convertTargetedNuctide(self, Nuctide pep , LightCompound p ):
        """
        convertTargetedNuctide(self, pep: Nuctide , p: LightCompound ) -> None
        Converts from the OpenMS TargetedExperiment Nuctide to the LightTargetedExperiment Nuctide
        """
        assert isinstance(pep, Nuctide), 'arg pep wrong type'
        assert isinstance(p, LightCompound), 'arg p wrong type'
    
    
        self.inst.get().convertTargetedNuctide((deref(pep.inst.get())), (deref(p.inst.get())))
    
    def convertToSpectrumPtr(self, MSSpectrum spectrum):
        assert isinstance(spectrum, MSSpectrum), 'arg spec wrong type'

        _r = self.inst.get().convertToSpectrumPtr((deref(spectrum.inst.get())))
        # cdef shared_ptr[_OSSpectrum] py_result = _r
        spec = OSSpectrum()
        spec.inst = _r
        return spec

    def convertToChromatogramPtr(self, MSChromatogram chrom):
        assert isinstance(chrom, MSChromatogram), 'arg spec wrong type'

        _r = self.inst.get().convertToChromatogramPtr((deref(chrom.inst.get())))
        # cdef shared_ptr[_OSChromatogram] py_result = _r
        ch = OSChromatogram()
        ch.inst = _r
        return ch 

cdef class ParamValue:
    """
    Cython implementation of _ParamValue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamValue.html>`_

    Class to hold strings, numeric values, vectors of strings and vectors of numeric values using the stl types
    
    - To choose one of these types, just use the appropriate constructor
    - Automatic conversion is supported and throws Exceptions in case of invalid conversions
    - An empty object is created with the default constructor
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ParamValue rv = ParamValue.__new__(ParamValue)
       rv.inst = shared_ptr[_ParamValue](new _ParamValue(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ParamValue rv = ParamValue.__new__(ParamValue)
       rv.inst = shared_ptr[_ParamValue](new _ParamValue(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ParamValue](new _ParamValue())
    
    def _init_1(self, ParamValue in_0 ):
        """
        _init_1(self, in_0: ParamValue ) -> None
        """
        assert isinstance(in_0, ParamValue), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamValue](new _ParamValue((deref(in_0.inst.get()))))
    
    def _init_2(self, bytes in_0 ):
        """
        _init_2(self, in_0: bytes ) -> None
        """
        assert isinstance(in_0, bytes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamValue](new _ParamValue((<char *>in_0)))
    
    def _init_3(self,  in_0 ):
        """
        _init_3(self, in_0: Union[bytes, str] ) -> None
        """
        assert isinstance(in_0, (bytes, str)), 'arg in_0 wrong type'
        if isinstance(in_0, str):
            in_0 = in_0.encode('utf-8')
        self.inst = shared_ptr[_ParamValue](new _ParamValue((<libcpp_string>in_0)))
    
    def _init_4(self,  in_0 ):
        """
        _init_4(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamValue](new _ParamValue((<int>in_0)))
    
    def _init_5(self, double in_0 ):
        """
        _init_5(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamValue](new _ParamValue((<double>in_0)))
    
    def _init_6(self, list in_0 ):
        """
        _init_6(self, in_0: List[Union[bytes, str]] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[libcpp_utf8_string] v0 = in_0
        self.inst = shared_ptr[_ParamValue](new _ParamValue(v0))
        
    
    def _init_7(self, list in_0 ):
        """
        _init_7(self, in_0: List[int] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, int) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[int] v0 = in_0
        self.inst = shared_ptr[_ParamValue](new _ParamValue(v0))
        
    
    def _init_8(self, list in_0 ):
        """
        _init_8(self, in_0: List[float] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] v0 = in_0
        self.inst = shared_ptr[_ParamValue](new _ParamValue(v0))
        
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ParamValue ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: bytes ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[bytes, str] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[Union[bytes, str]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[float] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ParamValue)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], bytes)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], (bytes, str))):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
             self._init_4(*args)
        elif (len(args)==1) and (isinstance(args[0], float)):
             self._init_5(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in args[0])):
             self._init_6(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[0])):
             self._init_7(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])):
             self._init_8(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toInt(self):
        cdef int _r = <int>(deref(self.inst.get()))
        py_res = <int>_r
        return py_res
    
    def toString(self):
        cdef libcpp_utf8_output_string _r = <libcpp_utf8_output_string>(deref(self.inst.get()))
        py_res = _r.decode('utf-8')
        return py_res
    
    def toDouble(self):
        cdef double _r = <double>(deref(self.inst.get()))
        py_res = <double>_r
        return py_res
    
    def toStringVector(self):
        """
        toStringVector(self) -> List[bytes]
        Explicitly convert ParamValue to string vector
        """
        _r = self.inst.get().toStringVector()
        cdef list py_result = _r
        return py_result
    
    def toDoubleVector(self):
        """
        toDoubleVector(self) -> List[float]
        Explicitly convert ParamValue to DoubleList
        """
        _r = self.inst.get().toDoubleVector()
        cdef list py_result = _r
        return py_result
    
    def toIntVector(self):
        """
        toIntVector(self) -> List[int]
        Explicitly convert ParamValue to IntList
        """
        _r = self.inst.get().toIntVector()
        cdef list py_result = _r
        return py_result
    
    def toBool(self):
        """
        toBool(self) -> bool
        Converts the strings 'true' and 'false' to a bool
        """
        cdef bool _r = self.inst.get().toBool()
        py_result = <bool>_r
        return py_result
    
    def valueType(self):
        """
        valueType(self) -> int
        """
        cdef _ValueType _r = self.inst.get().valueType()
        py_result = <int>_r
        return py_result
    
    def isEmpty(self):
        """
        isEmpty(self) -> int
        Test if the value is empty
        """
        cdef int _r = self.inst.get().isEmpty()
        py_result = <int>_r
        return py_result 

cdef class SolverParam:
    """
    Cython implementation of _SolverParam

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SolverParam.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property message_level:
        def __set__(self,  message_level):
        
            self.inst.get().message_level = (<int>message_level)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().message_level
            py_result = <int>_r
            return py_result
    
    property branching_tech:
        def __set__(self,  branching_tech):
        
            self.inst.get().branching_tech = (<int>branching_tech)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().branching_tech
            py_result = <int>_r
            return py_result
    
    property backtrack_tech:
        def __set__(self,  backtrack_tech):
        
            self.inst.get().backtrack_tech = (<int>backtrack_tech)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().backtrack_tech
            py_result = <int>_r
            return py_result
    
    property preprocessing_tech:
        def __set__(self,  preprocessing_tech):
        
            self.inst.get().preprocessing_tech = (<int>preprocessing_tech)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().preprocessing_tech
            py_result = <int>_r
            return py_result
    
    property enable_feas_pump_heuristic:
        def __set__(self, bool enable_feas_pump_heuristic):
        
            self.inst.get().enable_feas_pump_heuristic = (<bool>enable_feas_pump_heuristic)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_feas_pump_heuristic
            py_result = <bool>_r
            return py_result
    
    property enable_gmi_cuts:
        def __set__(self, bool enable_gmi_cuts):
        
            self.inst.get().enable_gmi_cuts = (<bool>enable_gmi_cuts)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_gmi_cuts
            py_result = <bool>_r
            return py_result
    
    property enable_mir_cuts:
        def __set__(self, bool enable_mir_cuts):
        
            self.inst.get().enable_mir_cuts = (<bool>enable_mir_cuts)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_mir_cuts
            py_result = <bool>_r
            return py_result
    
    property enable_cov_cuts:
        def __set__(self, bool enable_cov_cuts):
        
            self.inst.get().enable_cov_cuts = (<bool>enable_cov_cuts)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_cov_cuts
            py_result = <bool>_r
            return py_result
    
    property enable_clq_cuts:
        def __set__(self, bool enable_clq_cuts):
        
            self.inst.get().enable_clq_cuts = (<bool>enable_clq_cuts)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_clq_cuts
            py_result = <bool>_r
            return py_result
    
    property mip_gap:
        def __set__(self, double mip_gap):
        
            self.inst.get().mip_gap = (<double>mip_gap)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mip_gap
            py_result = <double>_r
            return py_result
    
    property time_limit:
        def __set__(self,  time_limit):
        
            self.inst.get().time_limit = (<int>time_limit)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().time_limit
            py_result = <int>_r
            return py_result
    
    property output_freq:
        def __set__(self,  output_freq):
        
            self.inst.get().output_freq = (<int>output_freq)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().output_freq
            py_result = <int>_r
            return py_result
    
    property output_delay:
        def __set__(self,  output_delay):
        
            self.inst.get().output_delay = (<int>output_delay)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().output_delay
            py_result = <int>_r
            return py_result
    
    property enable_presolve:
        def __set__(self, bool enable_presolve):
        
            self.inst.get().enable_presolve = (<bool>enable_presolve)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_presolve
            py_result = <bool>_r
            return py_result
    
    property enable_binarization:
        def __set__(self, bool enable_binarization):
        
            self.inst.get().enable_binarization = (<bool>enable_binarization)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().enable_binarization
            py_result = <bool>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        Hold the parameters of the LP solver
        """
        self.inst = shared_ptr[_SolverParam](new _SolverParam()) 
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
