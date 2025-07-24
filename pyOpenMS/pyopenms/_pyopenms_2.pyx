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
def __static_TransformationDescription_getModelTypes(list result ):
    """
    __static_TransformationDescription_getModelTypes(result: List[bytes] ) -> None
    """
    assert isinstance(result, list) and all(isinstance(li, bytes) for li in result), 'arg result wrong type'
    cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
    cdef bytes item0
    for item0 in result:
       v0.push_back(_String(<char *>item0))
    _getModelTypes_TransformationDescription(deref(v0))
    del v0

def __static_SpectrumHelper_removePeaks(MSChromatogram p , double pos_start , double pos_end ):
    """
    __static_SpectrumHelper_removePeaks(p: MSChromatogram , pos_start: float , pos_end: float ) -> None
    """
    assert isinstance(p, MSChromatogram), 'arg p wrong type'
    assert isinstance(pos_start, float), 'arg pos_start wrong type'
    assert isinstance(pos_end, float), 'arg pos_end wrong type'



    _removePeaks_SpectrumHelper((deref(p.inst.get())), (<double>pos_start), (<double>pos_end))

def __static_SpectrumHelper_removePeaks(MSSpectrum p , double pos_start , double pos_end ):
    """
    __static_SpectrumHelper_removePeaks(p: MSSpectrum , pos_start: float , pos_end: float ) -> None
    """
    assert isinstance(p, MSSpectrum), 'arg p wrong type'
    assert isinstance(pos_start, float), 'arg pos_start wrong type'
    assert isinstance(pos_end, float), 'arg pos_end wrong type'



    _removePeaks_SpectrumHelper((deref(p.inst.get())), (<double>pos_start), (<double>pos_end))

def __static_SpectrumHelper_subtractMinimumIntensity(MSChromatogram p ):
    """
    __static_SpectrumHelper_subtractMinimumIntensity(p: MSChromatogram ) -> None
    """
    assert isinstance(p, MSChromatogram), 'arg p wrong type'

    _subtractMinimumIntensity_SpectrumHelper((deref(p.inst.get())))

def __static_SpectrumHelper_subtractMinimumIntensity(MSSpectrum p ):
    """
    __static_SpectrumHelper_subtractMinimumIntensity(p: MSSpectrum ) -> None
    """
    assert isinstance(p, MSSpectrum), 'arg p wrong type'

    _subtractMinimumIntensity_SpectrumHelper((deref(p.inst.get())))

def __static_IMTypes_toDriftTimeUnit(bytes dtu_string ):
    """
    __static_IMTypes_toDriftTimeUnit(dtu_string: bytes ) -> int
    """
    assert isinstance(dtu_string, bytes), 'arg dtu_string wrong type'

    cdef _DriftTimeUnit _r = _toDriftTimeUnit_IMTypes((<libcpp_string>dtu_string))
    py_result = <int>_r
    return py_result

def __static_IMTypes_toIMFormat(bytes IM_format ):
    """
    __static_IMTypes_toIMFormat(IM_format: bytes ) -> int
    """
    assert isinstance(IM_format, bytes), 'arg IM_format wrong type'

    cdef _IMFormat _r = _toIMFormat_IMTypes((<libcpp_string>IM_format))
    py_result = <int>_r
    return py_result

def __static_IMTypes_toString(int value ):
    """
    __static_IMTypes_toString(value: int ) -> bytes
    """
    assert value in [0, 1, 2, 3, 4], 'arg value wrong type'

    cdef libcpp_string _r = _toString_IMTypes((<_DriftTimeUnit>value))
    py_result = <libcpp_string>_r
    return py_result

def __static_IMTypes_toString(int value ):
    """
    __static_IMTypes_toString(value: int ) -> bytes
    """
    assert value in [0, 1, 2, 3, 4], 'arg value wrong type'

    cdef libcpp_string _r = _toString_IMTypes((<_IMFormat>value))
    py_result = <libcpp_string>_r
    return py_result 

cdef class DriftTimeUnit:
    None
    NONE = 0
    MILLISECOND = 1
    VSSC = 2
    FAIMS_COMPENSATION_VOLTAGE = 3
    SIZE_OF_DRIFTTIMEUNIT = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class IMFormat:
    None
    NONE = 0
    CONCATENATED = 1
    MULTIPLE_SPECTRA = 2
    MIXED = 3
    SIZE_OF_IMFORMAT = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __NumpressCompression:
    None
    NONE = 0
    LINEAR = 1
    PIC = 2
    SLOF = 3
    SIZE_OF_NUMPRESSCOMPRESSION = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class OriginAnnotationFormat:
    None
    FILE_ORIGIN = 0
    MAP_INDEX = 1
    ID_MERGE_INDEX = 2
    UNKNOWN_OAF = 3
    SIZE_OF_ORIGIN_ANNOTATION_FORMAT = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ScanMode:
    None
    UNKNOWN = 0
    MASSSPECTRUM = 1
    MS1SPECTRUM = 2
    MSNSPECTRUM = 3
    SIM = 4
    SRM = 5
    CRM = 6
    CNG = 7
    CNL = 8
    PRECURSOR = 9
    EMC = 10
    TDF = 11
    EMR = 12
    EMISSION = 13
    ABSORPTION = 14
    SIZE_OF_SCANMODE = 15

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class XRefType_CVTerm_ControlledVocabulary:
    None
    XSD_STRING = 0
    XSD_INTEGER = 1
    XSD_DECIMAL = 2
    XSD_NEGATIVE_INTEGER = 3
    XSD_POSITIVE_INTEGER = 4
    XSD_NON_NEGATIVE_INTEGER = 5
    XSD_NON_POSITIVE_INTEGER = 6
    XSD_BOOLEAN = 7
    XSD_DATE = 8
    XSD_ANYURI = 9
    NONE = 10

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class CVTerm_ControlledVocabulary:
    """
    Cython implementation of _CVTerm_ControlledVocabulary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTerm_ControlledVocabulary.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property name:
        def __set__(self,  name):
        
            self.inst.get().name = deref((convString(name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().name
            py_result = convOutputString(_r)
            return py_result
    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property parents:
        def __set__(self, set parents):
            cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
            cdef bytes item0
            for item0 in parents:
               v0.insert(_String(<char *>item0))
            self.inst.get().parents = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().parents
            py_result = set()
            cdef libcpp_set[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.add(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property children:
        def __set__(self, set children):
            cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
            cdef bytes item0
            for item0 in children:
               v0.insert(_String(<char *>item0))
            self.inst.get().children = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().children
            py_result = set()
            cdef libcpp_set[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.add(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property obsolete:
        def __set__(self, bool obsolete):
        
            self.inst.get().obsolete = (<bool>obsolete)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().obsolete
            py_result = <bool>_r
            return py_result
    
    property description:
        def __set__(self,  description):
        
            self.inst.get().description = deref((convString(description)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().description
            py_result = convOutputString(_r)
            return py_result
    
    property synonyms:
        def __set__(self, list synonyms):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in synonyms:
               v0.push_back(_String(<char *>item0))
            self.inst.get().synonyms = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().synonyms
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property unparsed:
        def __set__(self, list unparsed):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in unparsed:
               v0.push_back(_String(<char *>item0))
            self.inst.get().unparsed = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().unparsed
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property xref_type:
        def __set__(self, int xref_type):
        
            self.inst.get().xref_type = (<_XRefType_CVTerm_ControlledVocabulary>xref_type)
        
    
        def __get__(self):
            cdef _XRefType_CVTerm_ControlledVocabulary _r = self.inst.get().xref_type
            py_result = <int>_r
            return py_result
    
    property xref_binary:
        def __set__(self, list xref_binary):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in xref_binary:
               v0.push_back(_String(<char *>item0))
            self.inst.get().xref_binary = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().xref_binary
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property units:
        def __set__(self, set units):
            cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
            cdef bytes item0
            for item0 in units:
               v0.insert(_String(<char *>item0))
            self.inst.get().units = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().units
            py_result = set()
            cdef libcpp_set[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.add(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVTerm_ControlledVocabulary](new _CVTerm_ControlledVocabulary())
    
    def _init_1(self, CVTerm_ControlledVocabulary rhs ):
        """
        _init_1(self, rhs: CVTerm_ControlledVocabulary ) -> None
        """
        assert isinstance(rhs, CVTerm_ControlledVocabulary), 'arg rhs wrong type'
    
        self.inst = shared_ptr[_CVTerm_ControlledVocabulary](new _CVTerm_ControlledVocabulary((deref(rhs.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, rhs: CVTerm_ControlledVocabulary ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVTerm_ControlledVocabulary)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _toXMLString_0(self,  ref ,  value ):
        """
        _toXMLString_0(self, ref: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
    
        cdef _String _r = self.inst.get().toXMLString(deref((convString(ref)).get()), deref((convString(value)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def _toXMLString_1(self,  ref ,  value ):
        """
        _toXMLString_1(self, ref: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[bytes, str, String]
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
        assert isinstance(value, (int, float, list, bytes, str)), 'arg value wrong type'
    
    
        cdef _String _r = self.inst.get().toXMLString(deref((convString(ref)).get()), deref(DataValue(value).inst.get()))
        py_result = convOutputString(_r)
        return py_result
    
    def toXMLString(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: toXMLString(self, ref: Union[bytes, str, String] , value: Union[bytes, str, String] ) -> Union[bytes, str, String]
          :noindex:
        
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available

        
        .. rubric:: Overload:
        .. py:function:: toXMLString(self, ref: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[bytes, str, String]
          :noindex:
        
        Get mzidentml formatted string. i.e. a cvparam xml element, ref should be the name of the ControlledVocabulary (i.e. cv.name()) containing the CVTerm (e.g. PSI-MS for the psi-ms.obo - gets loaded in all cases like that??), value can be empty if not available
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._toXMLString_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._toXMLString_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getXRefTypeName(self, int type ):
        """
        getXRefTypeName(self, type: int ) -> Union[bytes, str, String]
        """
        assert type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 'arg type wrong type'
    
        cdef _String _r = self.inst.get().getXRefTypeName((<_XRefType_CVTerm_ControlledVocabulary>type))
        py_result = convOutputString(_r)
        return py_result
    
    def isHigherBetterScore(self, CVTerm_ControlledVocabulary term ):
        """
        isHigherBetterScore(self, term: CVTerm_ControlledVocabulary ) -> bool
        """
        assert isinstance(term, CVTerm_ControlledVocabulary), 'arg term wrong type'
    
        cdef bool _r = self.inst.get().isHigherBetterScore((deref(term.inst.get())))
        py_result = <bool>_r
        return py_result 

cdef class ControlledVocabulary:
    """
    Cython implementation of _ControlledVocabulary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ControlledVocabulary.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ControlledVocabulary rv = ControlledVocabulary.__new__(ControlledVocabulary)
       rv.inst = shared_ptr[_ControlledVocabulary](new _ControlledVocabulary(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ControlledVocabulary rv = ControlledVocabulary.__new__(ControlledVocabulary)
       rv.inst = shared_ptr[_ControlledVocabulary](new _ControlledVocabulary(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ControlledVocabulary](new _ControlledVocabulary())
    
    def _init_1(self, ControlledVocabulary in_0 ):
        """
        _init_1(self, in_0: ControlledVocabulary ) -> None
        """
        assert isinstance(in_0, ControlledVocabulary), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ControlledVocabulary](new _ControlledVocabulary((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ControlledVocabulary ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ControlledVocabulary)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def name(self):
        """
        name(self) -> Union[bytes, str, String]
        Returns the CV name (set in the load method)
        """
        cdef _String _r = self.inst.get().name()
        py_result = convOutputString(_r)
        return py_result
    
    def loadFromOBO(self,  name ,  filename ):
        """
        loadFromOBO(self, name: Union[bytes, str, String] , filename: Union[bytes, str, String] ) -> None
        Loads the CV from an OBO file
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
    
        self.inst.get().loadFromOBO(deref((convString(name)).get()), deref((convString(filename)).get()))
    
    def exists(self,  id ):
        """
        exists(self, id: Union[bytes, str, String] ) -> bool
        Returns true if the term is in the CV. Returns false otherwise.
        """
        assert (isinstance(id, str) or isinstance(id, bytes) or isinstance(id, String)), 'arg id wrong type'
    
        cdef bool _r = self.inst.get().exists(deref((convString(id)).get()))
        py_result = <bool>_r
        return py_result
    
    def hasTermWithName(self,  name ):
        """
        hasTermWithName(self, name: Union[bytes, str, String] ) -> bool
        Returns true if a term with the given name is in the CV. Returns false otherwise
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasTermWithName(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    def getTerm(self,  id ):
        """
        getTerm(self, id: Union[bytes, str, String] ) -> CVTerm_ControlledVocabulary
        Returns a term specified by ID
        """
        assert (isinstance(id, str) or isinstance(id, bytes) or isinstance(id, String)), 'arg id wrong type'
    
        cdef _CVTerm_ControlledVocabulary * _r = new _CVTerm_ControlledVocabulary(self.inst.get().getTerm(deref((convString(id)).get())))
        cdef CVTerm_ControlledVocabulary py_result = CVTerm_ControlledVocabulary.__new__(CVTerm_ControlledVocabulary)
        py_result.inst = shared_ptr[_CVTerm_ControlledVocabulary](_r)
        return py_result
    
    def getTermByName(self,  name ,  desc ):
        """
        getTermByName(self, name: Union[bytes, str, String] , desc: Union[bytes, str, String] ) -> CVTerm_ControlledVocabulary
        Returns a term specified by name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(desc, str) or isinstance(desc, bytes) or isinstance(desc, String)), 'arg desc wrong type'
    
    
        cdef _CVTerm_ControlledVocabulary * _r = new _CVTerm_ControlledVocabulary(self.inst.get().getTermByName(deref((convString(name)).get()), deref((convString(desc)).get())))
        cdef CVTerm_ControlledVocabulary py_result = CVTerm_ControlledVocabulary.__new__(CVTerm_ControlledVocabulary)
        py_result.inst = shared_ptr[_CVTerm_ControlledVocabulary](_r)
        return py_result
    
    def getAllChildTerms(self, set terms ,  parent ):
        """
        getAllChildTerms(self, terms: Set[bytes] , parent: Union[bytes, str, String] ) -> None
        Writes all child terms recursively into terms
        """
        assert isinstance(terms, set) and all(isinstance(i, bytes) for i in terms), 'arg terms wrong type'
        assert (isinstance(parent, str) or isinstance(parent, bytes) or isinstance(parent, String)), 'arg parent wrong type'
        cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
        cdef bytes item0
        for item0 in terms:
           v0.insert(_String(<char *>item0))
    
        self.inst.get().getAllChildTerms(deref(v0), deref((convString(parent)).get()))
        del v0
    
    def isChildOf(self,  child ,  parent ):
        """
        isChildOf(self, child: Union[bytes, str, String] , parent: Union[bytes, str, String] ) -> bool
        Returns True if `child` is a child of `parent`
        """
        assert (isinstance(child, str) or isinstance(child, bytes) or isinstance(child, String)), 'arg child wrong type'
        assert (isinstance(parent, str) or isinstance(parent, bytes) or isinstance(parent, String)), 'arg parent wrong type'
    
    
        cdef bool _r = self.inst.get().isChildOf(deref((convString(child)).get()), deref((convString(parent)).get()))
        py_result = <bool>_r
        return py_result 

cdef class FloatDataArray:
    """
    Cython implementation of _FloatDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1FloatDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra float data attached to a spectrum or chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FloatDataArray rv = FloatDataArray.__new__(FloatDataArray)
       rv.inst = shared_ptr[_FloatDataArray](new _FloatDataArray(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FloatDataArray rv = FloatDataArray.__new__(FloatDataArray)
       rv.inst = shared_ptr[_FloatDataArray](new _FloatDataArray(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FloatDataArray](new _FloatDataArray())
    
    def _init_1(self, FloatDataArray in_0 ):
        """
        _init_1(self, in_0: FloatDataArray ) -> None
        """
        assert isinstance(in_0, FloatDataArray), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FloatDataArray](new _FloatDataArray((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FloatDataArray ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FloatDataArray)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def resize(self,  n ):
        """
        resize(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().resize((<size_t>n))
    
    def reserve(self,  n ):
        """
        reserve(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().reserve((<size_t>n))
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def push_back(self, float in_0 ):
        """
        push_back(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().push_back((<float>in_0))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the peak annotations
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the peak annotations
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        Returns a reference to the description of the applied processing
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
        Sets the description of the applied processing
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
        if not isinstance(other, FloatDataArray):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef FloatDataArray other_casted = other
        cdef FloatDataArray self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    


    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert in_0 >= 0, 'arg in_0 cannot be negative'

        cdef unsigned int _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <float>_r
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, int), 'arg key wrong type'
        assert isinstance(value, (int, float)), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef float _v = (<float>value)
        deref( self.inst.get() )[(<int>key)] = _v

    def get_data(self):
        """
        Gets the raw data for the float data array.
        .. warning::
           This returns a fast but unsafe view on the underlying vector data. Make sure that the object you are getting the data from survives the usage of this view! Specifically, do NOT use it like this: :code:`data = spectrum.getFloatDataArrays()[0].get_data()` since the underlying FloatDataArray is temporary for this line and will most likely be garbage collected right after that line.

        Example usage: 
        
        .. code-block:: python
        
            fd = pyopenms.FloatDataArray()
            data = fd.get_data()

        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        cdef unsigned int n = fda_.size()
         
        # Obtain a raw ptr to the beginning of the C++ array
        cdef libcpp_vector[float] * vec_ptr = <libcpp_vector[float]*> fda_
        cdef float * raw_ptr =  address(deref(vec_ptr)[0]) 

        # We use a memory view to get the data from the raw data
        # See https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html 
        # See https://stackoverflow.com/questions/43021574/cast-c-array-into-numpy-array-cython-typed-memoryview-in-cython-code
        cdef float[:] fda_view = <float[:n]>raw_ptr # cast to memoryview, refer to the underlying buffer without copy
        xarr = np.asarray(fda_view) # numpy array refer to the underlying buffer without copy
        return xarr

        ## # Alternatively, we could also use the assignment method:
        ## cdef _FloatDataArray * fda_ = self.inst.get()
        ## cdef libcpp_vector[float].iterator it = fda_.begin()
        ## cdef unsigned int n = fda_.size()
        ## cdef np.ndarray[np.float32_t, ndim=1] data
        ## data = np.zeros( (n,), dtype=np.float32)
        ## 
        ## cdef int i = 0
        ## while it != fda_.end():
        ##     data[i] = deref(it)
        ##     inc(it)
        ##     i += 1
        ## 
        ## return data
        ##

        # We see an improvement of ca 50x using the memory view method (and a
        # 2600x fold improvement compared to pure Python).

    def set_data(self, np.ndarray[float, ndim=1, mode="c"] data not None):
        """
        Sets the raw data for the float data array

        Example usage: 

          fd = pyopenms.FloatDataArray()
          data = numpy.array( [1, 2, 3, 5 ,6] ).astype(numpy.float32)
          fd.set_data(data)

        """
        cdef _FloatDataArray * fda_ = self.inst.get()
        fda_.clear()
        # Note: The np.ndarray[float, ndim=1, mode="c"] assures that we get a C-contiguous numpy array
        # See: https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC

        # We use "assign" to directly to copy the numpy array (stored in
        # data.data) into the std::vector<float> in our FloatDataArray object
        cdef float * array_start = <float*>data.data
        cdef int N
        N = data.size
        fda_.assign(array_start, array_start + N)


        ## # Alternatively, we could also "push_back" each element:
        ## fda_.reserve(N)
        ## for i in range(N):
        ##     fda_.push_back(data[i])
        
        # We see an improvement of ca 10x using the assign method (and an 600x
        # improvement compared to pure Python). 

cdef class IDRipper:
    """
    Cython implementation of _IDRipper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1IDRipper.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        Ripping protein/peptide identification according their file origin
        """
        self.inst = shared_ptr[_IDRipper](new _IDRipper())
    
    def rip(self, list rfis , list rfcs , list proteins , PeptideIdentificationList peptides , bool full_split , bool split_ident_runs ):
        """
        rip(self, rfis: List[RipFileIdentifier] , rfcs: List[RipFileContent] , proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList , full_split: bool , split_ident_runs: bool ) -> None
        """
        assert isinstance(rfis, list) and all(isinstance(elemt_rec, RipFileIdentifier) for elemt_rec in rfis), 'arg rfis wrong type'
        assert isinstance(rfcs, list) and all(isinstance(elemt_rec, RipFileContent) for elemt_rec in rfcs), 'arg rfcs wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(full_split, pybool_t), 'arg full_split wrong type'
        assert isinstance(split_ident_runs, pybool_t), 'arg split_ident_runs wrong type'
        cdef libcpp_vector[_RipFileIdentifier] * v0 = new libcpp_vector[_RipFileIdentifier]()
        cdef RipFileIdentifier item0
        for item0 in rfis:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_RipFileContent] * v1 = new libcpp_vector[_RipFileContent]()
        cdef RipFileContent item1
        for item1 in rfcs:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in proteins:
            v2.push_back(deref(item2.inst.get()))
    
    
    
        self.inst.get().rip(deref(v0), deref(v1), deref(v2), (deref(peptides.inst.get())), (<bool>full_split), (<bool>split_ident_runs))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v2.begin()
        replace_0 = []
        while it_proteins != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item2)
            inc(it_proteins)
        proteins[:] = replace_0
        del v2
        cdef libcpp_vector[_RipFileContent].iterator it_rfcs = v1.begin()
        replace_0 = []
        while it_rfcs != v1.end():
            item1 = RipFileContent.__new__(RipFileContent)
            item1.inst = shared_ptr[_RipFileContent](new _RipFileContent(deref(it_rfcs)))
            replace_0.append(item1)
            inc(it_rfcs)
        rfcs[:] = replace_0
        del v1
        cdef libcpp_vector[_RipFileIdentifier].iterator it_rfis = v0.begin()
        replace_0 = []
        while it_rfis != v0.end():
            item0 = RipFileIdentifier.__new__(RipFileIdentifier)
            item0.inst = shared_ptr[_RipFileIdentifier](new _RipFileIdentifier(deref(it_rfis)))
            replace_0.append(item0)
            inc(it_rfis)
        rfis[:] = replace_0
        del v0
    
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

cdef class IMTypes:
    """
    Cython implementation of _IMTypes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IMTypes.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IMTypes rv = IMTypes.__new__(IMTypes)
       rv.inst = shared_ptr[_IMTypes](new _IMTypes(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IMTypes rv = IMTypes.__new__(IMTypes)
       rv.inst = shared_ptr[_IMTypes](new _IMTypes(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IMTypes](new _IMTypes())
    
    def _init_1(self, IMTypes in_0 ):
        """
        _init_1(self, in_0: IMTypes ) -> None
        """
        assert isinstance(in_0, IMTypes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IMTypes](new _IMTypes((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IMTypes ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMTypes)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _determineIMFormat_0(self, MSExperiment exp ):
        """
        _determineIMFormat_0(self, exp: MSExperiment ) -> int
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        cdef _IMFormat _r = self.inst.get().determineIMFormat((deref(exp.inst.get())))
        py_result = <int>_r
        return py_result
    
    def _determineIMFormat_1(self, MSSpectrum spec ):
        """
        _determineIMFormat_1(self, spec: MSSpectrum ) -> int
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        cdef _IMFormat _r = self.inst.get().determineIMFormat((deref(spec.inst.get())))
        py_result = <int>_r
        return py_result
    
    def determineIMFormat(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: determineIMFormat(self, exp: MSExperiment ) -> int
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: determineIMFormat(self, spec: MSSpectrum ) -> int
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], MSExperiment)):
            return self._determineIMFormat_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSSpectrum)):
            return self._determineIMFormat_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    toDriftTimeUnit = __static_IMTypes_toDriftTimeUnit
    toIMFormat = __static_IMTypes_toIMFormat
    toString = __static_IMTypes_toString
    toString = __static_IMTypes_toString 

cdef class IdentificationRuns:
    """
    Cython implementation of _IdentificationRuns

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1IdentificationRuns.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list prot_ids ):
        """
        __init__(self, prot_ids: List[ProteinIdentification] ) -> None
        """
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in prot_ids:
            v0.push_back(deref(item0.inst.get()))
        self.inst = shared_ptr[_IdentificationRuns](new _IdentificationRuns(deref(v0)))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v0.begin()
        replace_0 = []
        while it_prot_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item0)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v0 

cdef class InstrumentSettings:
    """
    Cython implementation of _InstrumentSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InstrumentSettings.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef InstrumentSettings rv = InstrumentSettings.__new__(InstrumentSettings)
       rv.inst = shared_ptr[_InstrumentSettings](new _InstrumentSettings(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InstrumentSettings rv = InstrumentSettings.__new__(InstrumentSettings)
       rv.inst = shared_ptr[_InstrumentSettings](new _InstrumentSettings(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Description of the settings a MS Instrument was run with
        """
        self.inst = shared_ptr[_InstrumentSettings](new _InstrumentSettings())
    
    def _init_1(self, InstrumentSettings in_0 ):
        """
        _init_1(self, in_0: InstrumentSettings ) -> None
        """
        assert isinstance(in_0, InstrumentSettings), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InstrumentSettings](new _InstrumentSettings((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Description of the settings a MS Instrument was run with

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InstrumentSettings ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InstrumentSettings)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPolarity(self):
        """
        getPolarity(self) -> int
        Returns the polarity
        """
        cdef _Polarity _r = self.inst.get().getPolarity()
        py_result = <int>_r
        return py_result
    
    def setPolarity(self, int in_0 ):
        """
        setPolarity(self, in_0: int ) -> None
        Sets the polarity
        """
        assert in_0 in [0, 1, 2, 3], 'arg in_0 wrong type'
    
        self.inst.get().setPolarity((<_Polarity>in_0))
    
    def getScanMode(self):
        """
        getScanMode(self) -> int
        Returns the scan mode
        """
        cdef _ScanMode _r = self.inst.get().getScanMode()
        py_result = <int>_r
        return py_result
    
    def setScanMode(self, int scan_mode ):
        """
        setScanMode(self, scan_mode: int ) -> None
        Sets the scan mode
        """
        assert scan_mode in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], 'arg scan_mode wrong type'
    
        self.inst.get().setScanMode((<_ScanMode>scan_mode))
    
    def getZoomScan(self):
        """
        getZoomScan(self) -> bool
        Returns if this scan is a zoom (enhanced resolution) scan
        """
        cdef bool _r = self.inst.get().getZoomScan()
        py_result = <bool>_r
        return py_result
    
    def setZoomScan(self, bool zoom_scan ):
        """
        setZoomScan(self, zoom_scan: bool ) -> None
        Sets if this scan is a zoom (enhanced resolution) scan
        """
        assert isinstance(zoom_scan, pybool_t), 'arg zoom_scan wrong type'
    
        self.inst.get().setZoomScan((<bool>zoom_scan))
    
    def getScanWindows(self):
        """
        getScanWindows(self) -> List[ScanWindow]
        Returns the m/z scan windows
        """
        _r = self.inst.get().getScanWindows()
        py_result = []
        cdef libcpp_vector[_ScanWindow].iterator it__r = _r.begin()
        cdef ScanWindow item_py_result
        while it__r != _r.end():
           item_py_result = ScanWindow.__new__(ScanWindow)
           item_py_result.inst = shared_ptr[_ScanWindow](new _ScanWindow(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setScanWindows(self, list scan_windows ):
        """
        setScanWindows(self, scan_windows: List[ScanWindow] ) -> None
        Sets the m/z scan windows
        """
        assert isinstance(scan_windows, list) and all(isinstance(elemt_rec, ScanWindow) for elemt_rec in scan_windows), 'arg scan_windows wrong type'
        cdef libcpp_vector[_ScanWindow] * v0 = new libcpp_vector[_ScanWindow]()
        cdef ScanWindow item0
        for item0 in scan_windows:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setScanWindows(deref(v0))
        del v0
    
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
        if not isinstance(other, InstrumentSettings):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef InstrumentSettings other_casted = other
        cdef InstrumentSettings self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class MRMFP_ComponentGroupParams:
    """
    Cython implementation of _MRMFP_ComponentGroupParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFP_ComponentGroupParams.html>`_
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
    
    property params:
        def __set__(self, Param params):
        
            self.inst.get().params = (deref(params.inst.get()))
        
    
        def __get__(self):
            cdef _Param * _r = new _Param(self.inst.get().params)
            cdef Param py_result = Param.__new__(Param)
            py_result.inst = shared_ptr[_Param](_r)
            return py_result
    
    def __copy__(self):
       cdef MRMFP_ComponentGroupParams rv = MRMFP_ComponentGroupParams.__new__(MRMFP_ComponentGroupParams)
       rv.inst = shared_ptr[_MRMFP_ComponentGroupParams](new _MRMFP_ComponentGroupParams(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFP_ComponentGroupParams rv = MRMFP_ComponentGroupParams.__new__(MRMFP_ComponentGroupParams)
       rv.inst = shared_ptr[_MRMFP_ComponentGroupParams](new _MRMFP_ComponentGroupParams(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFP_ComponentGroupParams](new _MRMFP_ComponentGroupParams())
    
    def _init_1(self, MRMFP_ComponentGroupParams in_0 ):
        """
        _init_1(self, in_0: MRMFP_ComponentGroupParams ) -> None
        """
        assert isinstance(in_0, MRMFP_ComponentGroupParams), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFP_ComponentGroupParams](new _MRMFP_ComponentGroupParams((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFP_ComponentGroupParams ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFP_ComponentGroupParams)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMFP_ComponentParams:
    """
    Cython implementation of _MRMFP_ComponentParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFP_ComponentParams.html>`_
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
    
    property component_group_name:
        def __set__(self,  component_group_name):
        
            self.inst.get().component_group_name = deref((convString(component_group_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().component_group_name
            py_result = convOutputString(_r)
            return py_result
    
    property params:
        def __set__(self, Param params):
        
            self.inst.get().params = (deref(params.inst.get()))
        
    
        def __get__(self):
            cdef _Param * _r = new _Param(self.inst.get().params)
            cdef Param py_result = Param.__new__(Param)
            py_result.inst = shared_ptr[_Param](_r)
            return py_result
    
    def __copy__(self):
       cdef MRMFP_ComponentParams rv = MRMFP_ComponentParams.__new__(MRMFP_ComponentParams)
       rv.inst = shared_ptr[_MRMFP_ComponentParams](new _MRMFP_ComponentParams(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFP_ComponentParams rv = MRMFP_ComponentParams.__new__(MRMFP_ComponentParams)
       rv.inst = shared_ptr[_MRMFP_ComponentParams](new _MRMFP_ComponentParams(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFP_ComponentParams](new _MRMFP_ComponentParams())
    
    def _init_1(self, MRMFP_ComponentParams in_0 ):
        """
        _init_1(self, in_0: MRMFP_ComponentParams ) -> None
        """
        assert isinstance(in_0, MRMFP_ComponentParams), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFP_ComponentParams](new _MRMFP_ComponentParams((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFP_ComponentParams ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFP_ComponentParams)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MRMFeaturePicker:
    """
    Cython implementation of _MRMFeaturePicker

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeaturePicker.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMFeaturePicker rv = MRMFeaturePicker.__new__(MRMFeaturePicker)
       rv.inst = shared_ptr[_MRMFeaturePicker](new _MRMFeaturePicker(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeaturePicker rv = MRMFeaturePicker.__new__(MRMFeaturePicker)
       rv.inst = shared_ptr[_MRMFeaturePicker](new _MRMFeaturePicker(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeaturePicker](new _MRMFeaturePicker())
    
    def _init_1(self, MRMFeaturePicker in_0 ):
        """
        _init_1(self, in_0: MRMFeaturePicker ) -> None
        """
        assert isinstance(in_0, MRMFeaturePicker), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeaturePicker](new _MRMFeaturePicker((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeaturePicker ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeaturePicker)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MSNumpressCoder:
    """
    Cython implementation of _MSNumpressCoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSNumpressCoder.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSNumpressCoder rv = MSNumpressCoder.__new__(MSNumpressCoder)
       rv.inst = shared_ptr[_MSNumpressCoder](new _MSNumpressCoder(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSNumpressCoder rv = MSNumpressCoder.__new__(MSNumpressCoder)
       rv.inst = shared_ptr[_MSNumpressCoder](new _MSNumpressCoder(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSNumpressCoder](new _MSNumpressCoder())
    
    def _init_1(self, MSNumpressCoder in_0 ):
        """
        _init_1(self, in_0: MSNumpressCoder ) -> None
        """
        assert isinstance(in_0, MSNumpressCoder), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSNumpressCoder](new _MSNumpressCoder((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSNumpressCoder ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSNumpressCoder)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def encodeNP(self, list in_ ,  result , bool zlib_compression , NumpressConfig config ):
        """
        encodeNP(self, in_: List[float] , result: String , zlib_compression: bool , config: NumpressConfig ) -> None
        Encodes a vector of floating point numbers into a Base64 string using numpress
        
        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
        
        This function will first apply the numpress encoding to the data, then
        encode the result in base64 (with optional zlib compression before
        base64 encoding)
        
        :note In case of error, result string is empty
        
        
        :param in: The vector of floating point numbers to be encoded
        :param result: The resulting string
        :param zlib_compression: Whether to apply zlib compression after numpress compression
        :param config: The numpress configuration defining the compression strategy
        """
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert isinstance(result, String), 'arg result wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
        cdef libcpp_vector[double] v0 = in_
    
    
    
        self.inst.get().encodeNP(v0, deref((<String>result).inst.get()), (<bool>zlib_compression), (deref(config.inst.get())))
        
    
    def decodeNP(self,  in_ , list out , bool zlib_compression , NumpressConfig config ):
        """
        decodeNP(self, in_: Union[bytes, str, String] , out: List[float] , zlib_compression: bool , config: NumpressConfig ) -> None
        Decodes a Base64 string to a vector of floating point numbers using numpress
        
        This code is obtained from the proteowizard implementation
        ./pwiz/pwiz/data/msdata/BinaryDataEncoder.cpp (adapted by Hannes Roest)
        
        This function will first decode the input base64 string (with optional
        zlib decompression after decoding) and then apply numpress decoding to
        the data
        
        
        :param in: The base64 encoded string
        :param out: The resulting vector of doubles
        :param zlib_compression: Whether to apply zlib de-compression before numpress de-compression
        :param config: The numpress configuration defining the compression strategy
        :raises:
          Exception: ConversionError if the string cannot be converted
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
    
        cdef libcpp_vector[double] v1 = out
    
    
        self.inst.get().decodeNP(deref((convString(in_)).get()), v1, (<bool>zlib_compression), (deref(config.inst.get())))
        out[:] = v1
    
    def encodeNPRaw(self, list in_ ,  result , NumpressConfig config ):
        """
        encodeNPRaw(self, in_: List[float] , result: String , config: NumpressConfig ) -> None
        Encode the data vector "in" to a raw byte array
        
        :note In case of error, "result" is given back unmodified
        :note The result is not a string but a raw byte array and may contain zero bytes
        
        This performs the raw numpress encoding on a set of data and does no
        Base64 encoding on the result. Therefore the result string is likely
        *unsafe* to handle and is a raw byte array.
        
        Please use the safe versions above unless you need access to the raw
        byte arrays
        
        
        :param in: The vector of floating point numbers to be encoded
        :param result: The resulting string
        :param config: The numpress configuration defining the compression strategy
        """
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert isinstance(result, String), 'arg result wrong type'
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
        cdef libcpp_vector[double] v0 = in_
    
    
        self.inst.get().encodeNPRaw(v0, deref((<String>result).inst.get()), (deref(config.inst.get())))
        
    
    def decodeNPRaw(self,  in_ , list out , NumpressConfig config ):
        """
        decodeNPRaw(self, in_: Union[bytes, str, String] , out: List[float] , config: NumpressConfig ) -> None
        Decode the raw byte array "in" to the result vector "out"
        
        :note The string in should *only* contain the data and _no_ extra
        null terminating byte
        
        This performs the raw numpress decoding on a raw byte array (not Base64
        encoded). Therefore the input string is likely *unsafe* to handle and is
        basically a byte container
        
        Please use the safe versions above unless you need access to the raw
        byte arrays
        
        
        :param in: The base64 encoded string
        :param out: The resulting vector of doubles
        :param config: The numpress configuration defining the compression strategy
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(config, NumpressConfig), 'arg config wrong type'
    
        cdef libcpp_vector[double] v1 = out
    
        self.inst.get().decodeNPRaw(deref((convString(in_)).get()), v1, (deref(config.inst.get())))
        out[:] = v1
    NumpressCompression = __NumpressCompression 

cdef class MSstatsFile:
    """
    Cython implementation of _MSstatsFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSstatsFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSstatsFile rv = MSstatsFile.__new__(MSstatsFile)
       rv.inst = shared_ptr[_MSstatsFile](new _MSstatsFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSstatsFile rv = MSstatsFile.__new__(MSstatsFile)
       rv.inst = shared_ptr[_MSstatsFile](new _MSstatsFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSstatsFile](new _MSstatsFile())
    
    def _init_1(self, MSstatsFile in_0 ):
        """
        _init_1(self, in_0: MSstatsFile ) -> None
        """
        assert isinstance(in_0, MSstatsFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSstatsFile](new _MSstatsFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSstatsFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSstatsFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def storeLFQ(self,  filename , ConsensusMap consensus_map , ExperimentalDesign design , list reannotate_filenames , bool is_isotope_label_type ,  bioreplicate ,  condition ,  retention_time_summarization_method ):
        """
        storeLFQ(self, filename: String , consensus_map: ConsensusMap , design: ExperimentalDesign , reannotate_filenames: List[bytes] , is_isotope_label_type: bool , bioreplicate: String , condition: String , retention_time_summarization_method: String ) -> None
        Store label free experiment (MSstats)
        """
        assert isinstance(filename, String), 'arg filename wrong type'
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
        assert isinstance(design, ExperimentalDesign), 'arg design wrong type'
        assert isinstance(reannotate_filenames, list) and all(isinstance(li, bytes) for li in reannotate_filenames), 'arg reannotate_filenames wrong type'
        assert isinstance(is_isotope_label_type, pybool_t), 'arg is_isotope_label_type wrong type'
        assert isinstance(bioreplicate, String), 'arg bioreplicate wrong type'
        assert isinstance(condition, String), 'arg condition wrong type'
        assert isinstance(retention_time_summarization_method, String), 'arg retention_time_summarization_method wrong type'
    
    
    
        cdef libcpp_vector[_String] * v3 = new libcpp_vector[_String]()
        cdef bytes item3
        for item3 in reannotate_filenames:
           v3.push_back(_String(<char *>item3))
    
    
    
    
        self.inst.get().storeLFQ(deref((<String>filename).inst.get()), (deref(consensus_map.inst.get())), (deref(design.inst.get())), deref(v3), (<bool>is_isotope_label_type), deref((<String>bioreplicate).inst.get()), deref((<String>condition).inst.get()), deref((<String>retention_time_summarization_method).inst.get()))
        replace = []
        cdef libcpp_vector[_String].iterator it_3 = v3.begin()
        while it_3 != v3.end():
           replace.append(<char*>deref(it_3).c_str())
           inc(it_3)
        reannotate_filenames[:] = replace
        del v3
    
    def storeISO(self,  filename , ConsensusMap consensus_map , ExperimentalDesign design , list reannotate_filenames ,  bioreplicate ,  condition ,  mixture ,  retention_time_summarization_method ):
        """
        storeISO(self, filename: String , consensus_map: ConsensusMap , design: ExperimentalDesign , reannotate_filenames: List[bytes] , bioreplicate: String , condition: String , mixture: String , retention_time_summarization_method: String ) -> None
        Store isobaric experiment (MSstatsTMT)
        """
        assert isinstance(filename, String), 'arg filename wrong type'
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
        assert isinstance(design, ExperimentalDesign), 'arg design wrong type'
        assert isinstance(reannotate_filenames, list) and all(isinstance(li, bytes) for li in reannotate_filenames), 'arg reannotate_filenames wrong type'
        assert isinstance(bioreplicate, String), 'arg bioreplicate wrong type'
        assert isinstance(condition, String), 'arg condition wrong type'
        assert isinstance(mixture, String), 'arg mixture wrong type'
        assert isinstance(retention_time_summarization_method, String), 'arg retention_time_summarization_method wrong type'
    
    
    
        cdef libcpp_vector[_String] * v3 = new libcpp_vector[_String]()
        cdef bytes item3
        for item3 in reannotate_filenames:
           v3.push_back(_String(<char *>item3))
    
    
    
    
        self.inst.get().storeISO(deref((<String>filename).inst.get()), (deref(consensus_map.inst.get())), (deref(design.inst.get())), deref(v3), deref((<String>bioreplicate).inst.get()), deref((<String>condition).inst.get()), deref((<String>mixture).inst.get()), deref((<String>retention_time_summarization_method).inst.get()))
        replace = []
        cdef libcpp_vector[_String].iterator it_3 = v3.begin()
        while it_3 != v3.end():
           replace.append(<char*>deref(it_3).c_str())
           inc(it_3)
        reannotate_filenames[:] = replace
        del v3 

cdef class NumpressConfig:
    """
    Cython implementation of _NumpressConfig

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NumpressConfig.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property numpressFixedPoint:
        def __set__(self, double numpressFixedPoint):
        
            self.inst.get().numpressFixedPoint = (<double>numpressFixedPoint)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().numpressFixedPoint
            py_result = <double>_r
            return py_result
    
    property numpressErrorTolerance:
        def __set__(self, double numpressErrorTolerance):
        
            self.inst.get().numpressErrorTolerance = (<double>numpressErrorTolerance)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().numpressErrorTolerance
            py_result = <double>_r
            return py_result
    
    property np_compression:
        def __set__(self, int np_compression):
        
            self.inst.get().np_compression = (<_NumpressCompression>np_compression)
        
    
        def __get__(self):
            cdef _NumpressCompression _r = self.inst.get().np_compression
            py_result = <int>_r
            return py_result
    
    property estimate_fixed_point:
        def __set__(self, bool estimate_fixed_point):
        
            self.inst.get().estimate_fixed_point = (<bool>estimate_fixed_point)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().estimate_fixed_point
            py_result = <bool>_r
            return py_result
    
    property linear_fp_mass_acc:
        def __set__(self, double linear_fp_mass_acc):
        
            self.inst.get().linear_fp_mass_acc = (<double>linear_fp_mass_acc)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().linear_fp_mass_acc
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef NumpressConfig rv = NumpressConfig.__new__(NumpressConfig)
       rv.inst = shared_ptr[_NumpressConfig](new _NumpressConfig(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NumpressConfig rv = NumpressConfig.__new__(NumpressConfig)
       rv.inst = shared_ptr[_NumpressConfig](new _NumpressConfig(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NumpressConfig](new _NumpressConfig())
    
    def _init_1(self, NumpressConfig in_0 ):
        """
        _init_1(self, in_0: NumpressConfig ) -> None
        """
        assert isinstance(in_0, NumpressConfig), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NumpressConfig](new _NumpressConfig((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NumpressConfig ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NumpressConfig)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCompression(self,  compression ):
        """
        setCompression(self, compression: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(compression, str) or isinstance(compression, bytes) or isinstance(compression, String)), 'arg compression wrong type'
    
        self.inst.get().setCompression(deref((convString(compression)).get())) 

cdef class OPXLSpectrumProcessingAlgorithms:
    """
    Cython implementation of _OPXLSpectrumProcessingAlgorithms

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLSpectrumProcessingAlgorithms.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OPXLSpectrumProcessingAlgorithms rv = OPXLSpectrumProcessingAlgorithms.__new__(OPXLSpectrumProcessingAlgorithms)
       rv.inst = shared_ptr[_OPXLSpectrumProcessingAlgorithms](new _OPXLSpectrumProcessingAlgorithms(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OPXLSpectrumProcessingAlgorithms rv = OPXLSpectrumProcessingAlgorithms.__new__(OPXLSpectrumProcessingAlgorithms)
       rv.inst = shared_ptr[_OPXLSpectrumProcessingAlgorithms](new _OPXLSpectrumProcessingAlgorithms(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OPXLSpectrumProcessingAlgorithms](new _OPXLSpectrumProcessingAlgorithms())
    
    def _init_1(self, OPXLSpectrumProcessingAlgorithms in_0 ):
        """
        _init_1(self, in_0: OPXLSpectrumProcessingAlgorithms ) -> None
        """
        assert isinstance(in_0, OPXLSpectrumProcessingAlgorithms), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OPXLSpectrumProcessingAlgorithms](new _OPXLSpectrumProcessingAlgorithms((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OPXLSpectrumProcessingAlgorithms ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OPXLSpectrumProcessingAlgorithms)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def mergeAnnotatedSpectra(self, MSSpectrum first_spectrum , MSSpectrum second_spectrum ):
        """
        mergeAnnotatedSpectra(self, first_spectrum: MSSpectrum , second_spectrum: MSSpectrum ) -> MSSpectrum
        """
        assert isinstance(first_spectrum, MSSpectrum), 'arg first_spectrum wrong type'
        assert isinstance(second_spectrum, MSSpectrum), 'arg second_spectrum wrong type'
    
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().mergeAnnotatedSpectra((deref(first_spectrum.inst.get())), (deref(second_spectrum.inst.get()))))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def preprocessSpectra(self, MSExperiment exp , double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm ,  peptide_min_size ,  min_precursor_charge ,  max_precursor_charge , bool deisotope , bool labeled ):
        """
        preprocessSpectra(self, exp: MSExperiment , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peptide_min_size: int , min_precursor_charge: int , max_precursor_charge: int , deisotope: bool , labeled: bool ) -> MSExperiment
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(peptide_min_size, int) and peptide_min_size >= 0, 'arg peptide_min_size wrong type'
        assert isinstance(min_precursor_charge, int), 'arg min_precursor_charge wrong type'
        assert isinstance(max_precursor_charge, int), 'arg max_precursor_charge wrong type'
        assert isinstance(deisotope, pybool_t), 'arg deisotope wrong type'
        assert isinstance(labeled, pybool_t), 'arg labeled wrong type'
    
    
    
    
    
    
    
    
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().preprocessSpectra((deref(exp.inst.get())), (<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), (<size_t>peptide_min_size), (<int>min_precursor_charge), (<int>max_precursor_charge), (<bool>deisotope), (<bool>labeled)))
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def getSpectrumAlignmentFastCharge(self, list alignment , double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , MSSpectrum theo_spectrum , MSSpectrum exp_spectrum , IntegerDataArray theo_charges , IntegerDataArray exp_charges , FloatDataArray ppm_error_array , double intensity_cutoff ):
        """
        getSpectrumAlignmentFastCharge(self, alignment: List[List[int, int]] , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , theo_spectrum: MSSpectrum , exp_spectrum: MSSpectrum , theo_charges: IntegerDataArray , exp_charges: IntegerDataArray , ppm_error_array: FloatDataArray , intensity_cutoff: float ) -> None
        """
        assert isinstance(alignment, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in alignment), 'arg alignment wrong type'
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(theo_spectrum, MSSpectrum), 'arg theo_spectrum wrong type'
        assert isinstance(exp_spectrum, MSSpectrum), 'arg exp_spectrum wrong type'
        assert isinstance(theo_charges, IntegerDataArray), 'arg theo_charges wrong type'
        assert isinstance(exp_charges, IntegerDataArray), 'arg exp_charges wrong type'
        assert isinstance(ppm_error_array, FloatDataArray), 'arg ppm_error_array wrong type'
        assert isinstance(intensity_cutoff, float), 'arg intensity_cutoff wrong type'
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v0 = alignment
    
    
    
    
    
    
    
    
        self.inst.get().getSpectrumAlignmentFastCharge(v0, (<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), (deref(theo_spectrum.inst.get())), (deref(exp_spectrum.inst.get())), (deref(theo_charges.inst.get())), (deref(exp_charges.inst.get())), (deref(ppm_error_array.inst.get())), (<double>intensity_cutoff))
        alignment[:] = v0
    
    def getSpectrumAlignmentSimple(self, list alignment , double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , list theo_spectrum , MSSpectrum exp_spectrum , IntegerDataArray exp_charges ):
        """
        getSpectrumAlignmentSimple(self, alignment: List[List[int, int]] , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , theo_spectrum: List[SimplePeak] , exp_spectrum: MSSpectrum , exp_charges: IntegerDataArray ) -> None
        """
        assert isinstance(alignment, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in alignment), 'arg alignment wrong type'
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(theo_spectrum, list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in theo_spectrum), 'arg theo_spectrum wrong type'
        assert isinstance(exp_spectrum, MSSpectrum), 'arg exp_spectrum wrong type'
        assert isinstance(exp_charges, IntegerDataArray), 'arg exp_charges wrong type'
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v0 = alignment
    
    
        cdef libcpp_vector[_SimplePeak] * v3 = new libcpp_vector[_SimplePeak]()
        cdef SimplePeak item3
        for item3 in theo_spectrum:
            v3.push_back(deref(item3.inst.get()))
    
    
        self.inst.get().getSpectrumAlignmentSimple(v0, (<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), deref(v3), (deref(exp_spectrum.inst.get())), (deref(exp_charges.inst.get())))
        del v3
        alignment[:] = v0 

cdef class OSBinaryDataArray:
    """
    Cython implementation of _OSBinaryDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSBinaryDataArray.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property data:
        def __set__(self, list data):
            cdef libcpp_vector[double] v0 = data
            self.inst.get().data = v0
            
    
        def __get__(self):
            _r = self.inst.get().data
            cdef list py_result = _r
            return py_result
    
    property description:
        def __set__(self, bytes description):
        
            self.inst.get().description = (<libcpp_string>description)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().description
            py_result = <libcpp_string>_r
            return py_result
    
    def __copy__(self):
       cdef OSBinaryDataArray rv = OSBinaryDataArray.__new__(OSBinaryDataArray)
       rv.inst = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSBinaryDataArray rv = OSBinaryDataArray.__new__(OSBinaryDataArray)
       rv.inst = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray())
    
    def _init_1(self, OSBinaryDataArray in_0 ):
        """
        _init_1(self, in_0: OSBinaryDataArray ) -> None
        """
        assert isinstance(in_0, OSBinaryDataArray), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSBinaryDataArray ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSBinaryDataArray)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OSChromatogram:
    """
    Cython implementation of _OSChromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSChromatogram.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OSChromatogram rv = OSChromatogram.__new__(OSChromatogram)
       rv.inst = shared_ptr[_OSChromatogram](new _OSChromatogram(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSChromatogram rv = OSChromatogram.__new__(OSChromatogram)
       rv.inst = shared_ptr[_OSChromatogram](new _OSChromatogram(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OSChromatogram](new _OSChromatogram())
    
    def _init_1(self, OSChromatogram in_0 ):
        """
        _init_1(self, in_0: OSChromatogram ) -> None
        """
        assert isinstance(in_0, OSChromatogram), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSChromatogram](new _OSChromatogram((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSChromatogram ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSChromatogram)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getTimeArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def setTimeArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setTimeArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0) 

cdef class OSSpectrum:
    """
    Cython implementation of _OSSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSSpectrum.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OSSpectrum rv = OSSpectrum.__new__(OSSpectrum)
       rv.inst = shared_ptr[_OSSpectrum](new _OSSpectrum(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSSpectrum rv = OSSpectrum.__new__(OSSpectrum)
       rv.inst = shared_ptr[_OSSpectrum](new _OSSpectrum(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OSSpectrum](new _OSSpectrum())
    
    def _init_1(self, OSSpectrum in_0 ):
        """
        _init_1(self, in_0: OSSpectrum ) -> None
        """
        assert isinstance(in_0, OSSpectrum), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSSpectrum](new _OSSpectrum((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSSpectrum ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSSpectrum)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMZArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getMZArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        cdef shared_ptr[_OSBinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def setMZArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setMZArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_OSBinaryDataArray] v0 = shared_ptr[_OSBinaryDataArray](new _OSBinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0) 

cdef class RipFileContent:
    """
    Cython implementation of _RipFileContent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1RipFileContent.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list prot_idents , PeptideIdentificationList pep_idents ):
        """
        __init__(self, prot_idents: List[ProteinIdentification] , pep_idents: PeptideIdentificationList ) -> None
        """
        assert isinstance(prot_idents, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_idents), 'arg prot_idents wrong type'
        assert isinstance(pep_idents, PeptideIdentificationList), 'arg pep_idents wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in prot_idents:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_RipFileContent](new _RipFileContent(deref(v0), (deref(pep_idents.inst.get()))))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_idents = v0.begin()
        replace_0 = []
        while it_prot_idents != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_idents)))
            replace_0.append(item0)
            inc(it_prot_idents)
        prot_idents[:] = replace_0
        del v0
    
    def getProteinIdentifications(self):
        """
        getProteinIdentifications(self) -> List[ProteinIdentification]
        """
        _r = self.inst.get().getProteinIdentifications()
        py_result = []
        cdef libcpp_vector[_ProteinIdentification].iterator it__r = _r.begin()
        cdef ProteinIdentification item_py_result
        while it__r != _r.end():
           item_py_result = ProteinIdentification.__new__(ProteinIdentification)
           item_py_result.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getPeptideIdentifications(self):
        """
        getPeptideIdentifications(self) -> PeptideIdentificationList
        """
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().getPeptideIdentifications())
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result 

cdef class RipFileIdentifier:
    """
    Cython implementation of _RipFileIdentifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IDRipper_1_1RipFileIdentifier.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, IdentificationRuns id_runs , PeptideIdentification pep_id , dict file_origin_map , int origin_annotation_fmt , bool split_ident_runs ):
        """
        __init__(self, id_runs: IdentificationRuns , pep_id: PeptideIdentification , file_origin_map: Dict[Union[bytes, str, String], int] , origin_annotation_fmt: int , split_ident_runs: bool ) -> None
        """
        assert isinstance(id_runs, IdentificationRuns), 'arg id_runs wrong type'
        assert isinstance(pep_id, PeptideIdentification), 'arg pep_id wrong type'
        assert isinstance(file_origin_map, dict) and all((isinstance(k, str) or isinstance(k, bytes) or isinstance(k, String)) for k in file_origin_map.keys()) and all(isinstance(v, int) for v in file_origin_map.values()), 'arg file_origin_map wrong type'
        assert origin_annotation_fmt in [0, 1, 2, 3, 4], 'arg origin_annotation_fmt wrong type'
        assert isinstance(split_ident_runs, pybool_t), 'arg split_ident_runs wrong type'
    
    
        cdef libcpp_map[_String, unsigned int] * v2 = new libcpp_map[_String, unsigned int]()
        for key, value in file_origin_map.items():
        
        
            deref(v2)[ deref(<_String *> (<String> key).inst.get()) ] = (<unsigned int>value)
        
        
    
    
        self.inst = shared_ptr[_RipFileIdentifier](new _RipFileIdentifier((deref(id_runs.inst.get())), (deref(pep_id.inst.get())), deref(v2), (<_OriginAnnotationFormat>origin_annotation_fmt), (<bool>split_ident_runs)))
        replace = dict()
        cdef libcpp_map[_String, unsigned int].iterator it_file_origin_map = v2.begin()
        cdef String itemk_file_origin_map
        while it_file_origin_map != v2.end():
           itemk_file_origin_map = String.__new__(String)
           itemk_file_origin_map.inst = shared_ptr[_String](new _String((deref(it_file_origin_map)).first))
           replace[itemk_file_origin_map] = <unsigned int> deref(it_file_origin_map).second
           inc(it_file_origin_map)
        file_origin_map.clear()
        file_origin_map.update(replace)
        del v2
    
    def getIdentRunIdx(self):
        """
        getIdentRunIdx(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getIdentRunIdx()
        py_result = <unsigned int>_r
        return py_result
    
    def getFileOriginIdx(self):
        """
        getFileOriginIdx(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getFileOriginIdx()
        py_result = <unsigned int>_r
        return py_result
    
    def getOriginFullname(self):
        """
        getOriginFullname(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOriginFullname()
        py_result = convOutputString(_r)
        return py_result
    
    def getOutputBasename(self):
        """
        getOutputBasename(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOutputBasename()
        py_result = convOutputString(_r)
        return py_result 

cdef class SpectrumAccessSqMass:
    """
    Cython implementation of _SpectrumAccessSqMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessSqMass.html>`_
      -- Inherits from ['ISpectrumAccess']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAccessSqMass rv = SpectrumAccessSqMass.__new__(SpectrumAccessSqMass)
       rv.inst = shared_ptr[_SpectrumAccessSqMass](new _SpectrumAccessSqMass(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAccessSqMass rv = SpectrumAccessSqMass.__new__(SpectrumAccessSqMass)
       rv.inst = shared_ptr[_SpectrumAccessSqMass](new _SpectrumAccessSqMass(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        pass
    
    def _init_1(self, SpectrumAccessSqMass in_0 ):
        """
        _init_1(self, in_0: SpectrumAccessSqMass ) -> None
        """
        assert isinstance(in_0, SpectrumAccessSqMass), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessSqMass](new _SpectrumAccessSqMass((deref(in_0.inst.get()))))
    
    def _init_2(self, MzMLSqliteHandler in_0 , list indices ):
        """
        _init_2(self, in_0: MzMLSqliteHandler , indices: List[int] ) -> None
        """
        assert isinstance(in_0, MzMLSqliteHandler), 'arg in_0 wrong type'
        assert isinstance(indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in indices), 'arg indices wrong type'
    
        cdef libcpp_vector[int] v1 = indices
        self.inst = shared_ptr[_SpectrumAccessSqMass](new _SpectrumAccessSqMass((deref(in_0.inst.get())), v1))
        
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessSqMass ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzMLSqliteHandler , indices: List[int] ) -> None
          :noindex:
    
        """
        if kwargs.get("__createUnsafeObject__") is True:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessSqMass)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], MzMLSqliteHandler)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[1])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSpectrumById(self,  id_ ):
        """
        getSpectrumById(self, id_: int ) -> OSSpectrum
        Returns a pointer to a spectrum at the given string id
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_OSSpectrum] _r = self.inst.get().getSpectrumById((<int>id_))
        cdef OSSpectrum py_result
        py_result = OSSpectrum.__new__(OSSpectrum)
        py_result.inst = _r
        return py_result
    
    def getSpectraByRT(self, double RT , double deltaRT ):
        """
        getSpectraByRT(self, RT: float , deltaRT: float ) -> List[int]
        Returns a vector of ids of spectra that are within RT +/- deltaRT
        """
        assert isinstance(RT, float), 'arg RT wrong type'
        assert isinstance(deltaRT, float), 'arg deltaRT wrong type'
    
    
        _r = self.inst.get().getSpectraByRT((<double>RT), (<double>deltaRT))
        cdef list py_result = _r
        return py_result
    
    def getNrSpectra(self):
        """
        getNrSpectra(self) -> int
        Returns the number of spectra available
        """
        cdef size_t _r = self.inst.get().getNrSpectra()
        py_result = <size_t>_r
        return py_result
    
    def getChromatogramById(self,  id_ ):
        """
        getChromatogramById(self, id_: int ) -> OSChromatogram
        Returns a pointer to a chromatogram at the given id
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_OSChromatogram] _r = self.inst.get().getChromatogramById((<int>id_))
        cdef OSChromatogram py_result
        py_result = OSChromatogram.__new__(OSChromatogram)
        py_result.inst = _r
        return py_result
    
    def getNrChromatograms(self):
        """
        getNrChromatograms(self) -> int
        Returns the number of chromatograms available
        """
        cdef size_t _r = self.inst.get().getNrChromatograms()
        py_result = <size_t>_r
        return py_result
    
    def getChromatogramNativeID(self,  id_ ):
        """
        getChromatogramNativeID(self, id_: int ) -> str
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef libcpp_utf8_output_string _r = self.inst.get().getChromatogramNativeID((<int>id_))
        py_result = _r.decode('utf-8')
        return py_result 

cdef class SpectrumHelper:
    """
    Cython implementation of _SpectrumHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/class_1_1SpectrumHelper.html>`_
    """

    removePeaks = __static_SpectrumHelper_removePeaks
    removePeaks = __static_SpectrumHelper_removePeaks
    subtractMinimumIntensity = __static_SpectrumHelper_subtractMinimumIntensity
    subtractMinimumIntensity = __static_SpectrumHelper_subtractMinimumIntensity 

cdef class SwathMapMassCorrection:
    """
    Cython implementation of _SwathMapMassCorrection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathMapMassCorrection.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SwathMapMassCorrection](new _SwathMapMassCorrection())
    
    def _init_1(self, SwathMapMassCorrection in_0 ):
        """
        _init_1(self, in_0: SwathMapMassCorrection ) -> None
        """
        assert isinstance(in_0, SwathMapMassCorrection), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SwathMapMassCorrection](new _SwathMapMassCorrection((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SwathMapMassCorrection ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SwathMapMassCorrection)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class TransformationDescription:
    """
    Cython implementation of _TransformationDescription

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TransformationDescription_1_1TransformationDescription.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TransformationDescription rv = TransformationDescription.__new__(TransformationDescription)
       rv.inst = shared_ptr[_TransformationDescription](new _TransformationDescription(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TransformationDescription rv = TransformationDescription.__new__(TransformationDescription)
       rv.inst = shared_ptr[_TransformationDescription](new _TransformationDescription(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TransformationDescription](new _TransformationDescription())
    
    def _init_1(self, TransformationDescription in_0 ):
        """
        _init_1(self, in_0: TransformationDescription ) -> None
        """
        assert isinstance(in_0, TransformationDescription), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TransformationDescription](new _TransformationDescription((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TransformationDescription ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TransformationDescription)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getDataPoints(self):
        """
        getDataPoints(self) -> List[TM_DataPoint]
        Returns the data points
        """
        _r = self.inst.get().getDataPoints()
        py_result = []
        cdef libcpp_vector[_TM_DataPoint].iterator it__r = _r.begin()
        cdef TM_DataPoint item_py_result
        while it__r != _r.end():
           item_py_result = TM_DataPoint.__new__(TM_DataPoint)
           item_py_result.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def _setDataPoints_0(self, list data ):
        """
        _setDataPoints_0(self, data: List[TM_DataPoint] ) -> None
        Sets the data points. Removes the model that was previously fitted to the data (if any)
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setDataPoints(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def _setDataPoints_1(self, list data ):
        """
        _setDataPoints_1(self, data: List[List[float, float]] ) -> None
        Sets the data points (backwards-compatible overload). Removes the model that was previously fitted to the data (if any)
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[libcpp_pair[double,double]] v0 = data
        self.inst.get().setDataPoints(v0)
        data[:] = v0
    
    def setDataPoints(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setDataPoints(self, data: List[TM_DataPoint] ) -> None
          :noindex:
        
        Sets the data points. Removes the model that was previously fitted to the data (if any)

        
        .. rubric:: Overload:
        .. py:function:: setDataPoints(self, data: List[List[float, float]] ) -> None
          :noindex:
        
        Sets the data points (backwards-compatible overload). Removes the model that was previously fitted to the data (if any)
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in args[0])):
            return self._setDataPoints_0(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in args[0])):
            return self._setDataPoints_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def apply(self, double in_0 ):
        """
        apply(self, in_0: float ) -> float
        Applies the transformation to `value`
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().apply((<double>in_0))
        py_result = <double>_r
        return py_result
    
    def _fitModel_0(self,  model_type , Param params ):
        """
        _fitModel_0(self, model_type: Union[bytes, str, String] , params: Param ) -> None
        Fits a model to the data
        """
        assert (isinstance(model_type, str) or isinstance(model_type, bytes) or isinstance(model_type, String)), 'arg model_type wrong type'
        assert isinstance(params, Param), 'arg params wrong type'
    
    
        self.inst.get().fitModel(deref((convString(model_type)).get()), (deref(params.inst.get())))
    
    def _fitModel_1(self,  model_type ):
        """
        _fitModel_1(self, model_type: Union[bytes, str, String] ) -> None
        Fits a model to the data
        """
        assert (isinstance(model_type, str) or isinstance(model_type, bytes) or isinstance(model_type, String)), 'arg model_type wrong type'
    
        self.inst.get().fitModel(deref((convString(model_type)).get()))
    
    def fitModel(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: fitModel(self, model_type: Union[bytes, str, String] , params: Param ) -> None
          :noindex:
        
        Fits a model to the data

        
        .. rubric:: Overload:
        .. py:function:: fitModel(self, model_type: Union[bytes, str, String] ) -> None
          :noindex:
        
        Fits a model to the data
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], Param)):
            return self._fitModel_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._fitModel_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getModelType(self):
        """
        getModelType(self) -> Union[bytes, str, String]
        Gets the type of the fitted model
        """
        cdef _String _r = self.inst.get().getModelType()
        py_result = convOutputString(_r)
        return py_result
    
    def getModelParameters(self):
        """
        getModelParameters(self) -> Param
        Returns the model parameters
        """
        cdef _Param * _r = new _Param(self.inst.get().getModelParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def invert(self):
        """
        invert(self) -> None
        Computes an (approximate) inverse of the transformation
        """
        self.inst.get().invert()
    
    def getDeviations(self, list diffs , bool do_apply , bool do_sort ):
        """
        getDeviations(self, diffs: List[float] , do_apply: bool , do_sort: bool ) -> None
        Get the deviations between the data pairs
        
        :param diffs: Output
        :param do_apply: Get deviations after applying the model?
        :param do_sort: Sort `diffs` before returning?
        """
        assert isinstance(diffs, list) and all(isinstance(elemt_rec, float) for elemt_rec in diffs), 'arg diffs wrong type'
        assert isinstance(do_apply, pybool_t), 'arg do_apply wrong type'
        assert isinstance(do_sort, pybool_t), 'arg do_sort wrong type'
        cdef libcpp_vector[double] v0 = diffs
    
    
        self.inst.get().getDeviations(v0, (<bool>do_apply), (<bool>do_sort))
        diffs[:] = v0
    
    def getStatistics(self):
        """
        getStatistics(self) -> TransformationStatistics
        """
        cdef _TransformationStatistics * _r = new _TransformationStatistics(self.inst.get().getStatistics())
        cdef TransformationStatistics py_result = TransformationStatistics.__new__(TransformationStatistics)
        py_result.inst = shared_ptr[_TransformationStatistics](_r)
        return py_result
    getModelTypes = __static_TransformationDescription_getModelTypes 

cdef class TransformationStatistics:
    """
    Cython implementation of _TransformationStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TransformationDescription_1_1TransformationStatistics.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property xmin:
        def __set__(self, double xmin):
        
            self.inst.get().xmin = (<double>xmin)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xmin
            py_result = <double>_r
            return py_result
    
    property xmax:
        def __set__(self, double xmax):
        
            self.inst.get().xmax = (<double>xmax)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xmax
            py_result = <double>_r
            return py_result
    
    property ymin:
        def __set__(self, double ymin):
        
            self.inst.get().ymin = (<double>ymin)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ymin
            py_result = <double>_r
            return py_result
    
    property ymax:
        def __set__(self, double ymax):
        
            self.inst.get().ymax = (<double>ymax)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ymax
            py_result = <double>_r
            return py_result
    
    property percentiles_before:
        def __set__(self, dict percentiles_before):
            cdef libcpp_map[size_t, double] * v0 = new libcpp_map[size_t, double]()
            for key, value in percentiles_before.items():
            
            
                deref(v0)[ (<size_t>key) ] = (<double>value)
            
            
            self.inst.get().percentiles_before = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().percentiles_before
            py_result = dict()
            cdef libcpp_map[size_t, double].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result[<size_t>(deref(it__r).first)] = <double>(deref(it__r).second)
               inc(it__r)
            return py_result
    
    property percentiles_after:
        def __set__(self, dict percentiles_after):
            cdef libcpp_map[size_t, double] * v0 = new libcpp_map[size_t, double]()
            for key, value in percentiles_after.items():
            
            
                deref(v0)[ (<size_t>key) ] = (<double>value)
            
            
            self.inst.get().percentiles_after = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().percentiles_after
            py_result = dict()
            cdef libcpp_map[size_t, double].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result[<size_t>(deref(it__r).first)] = <double>(deref(it__r).second)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef TransformationStatistics rv = TransformationStatistics.__new__(TransformationStatistics)
       rv.inst = shared_ptr[_TransformationStatistics](new _TransformationStatistics(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TransformationStatistics rv = TransformationStatistics.__new__(TransformationStatistics)
       rv.inst = shared_ptr[_TransformationStatistics](new _TransformationStatistics(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TransformationStatistics](new _TransformationStatistics())
    
    def _init_1(self, TransformationStatistics in_0 ):
        """
        _init_1(self, in_0: TransformationStatistics ) -> None
        """
        assert isinstance(in_0, TransformationStatistics), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TransformationStatistics](new _TransformationStatistics((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TransformationStatistics ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TransformationStatistics)):
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
