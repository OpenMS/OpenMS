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

cdef class FileType:
    None
    UNKNOWN = 0
    DTA = 1
    DTA2D = 2
    MZDATA = 3
    MZXML = 4
    FEATUREXML = 5
    IDXML = 6
    CONSENSUSXML = 7
    MGF = 8
    INI = 9
    TOPPAS = 10
    TRANSFORMATIONXML = 11
    MZML = 12
    CACHEDMZML = 13
    MS2 = 14
    PEPXML = 15
    PROTXML = 16
    MZIDENTML = 17
    QCML = 18
    GELML = 19
    TRAML = 20
    MSP = 21
    OMSSAXML = 22
    MASCOTXML = 23
    PNG = 24
    XMASS = 25
    TSV = 26
    PEPLIST = 27
    HARDKLOER = 28
    KROENIK = 29
    FASTA = 30
    EDTA = 31
    CSV = 32
    TXT = 33
    OBO = 34
    HTML = 35
    XML = 36
    ANALYSISXML = 37
    XSD = 38
    PSQ = 39
    MRM = 40
    SQMASS = 41
    PQP = 42
    OSW = 43
    PSMS = 44
    PARAMXML = 45
    SIZE_OF_TYPE = 46

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AScore:
    """
    Cython implementation of _AScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AScore.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AScore rv = AScore.__new__(AScore)
       rv.inst = shared_ptr[_AScore](new _AScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AScore rv = AScore.__new__(AScore)
       rv.inst = shared_ptr[_AScore](new _AScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AScore](new _AScore())
    
    def _init_1(self, AScore in_0 ):
        """
        _init_1(self, in_0: AScore ) -> None
        """
        assert isinstance(in_0, AScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AScore](new _AScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AScore)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, PeptideHit hit , MSSpectrum real_spectrum ):
        """
        compute(self, hit: PeptideHit , real_spectrum: MSSpectrum ) -> PeptideHit
        """
        assert isinstance(hit, PeptideHit), 'arg hit wrong type'
        assert isinstance(real_spectrum, MSSpectrum), 'arg real_spectrum wrong type'
    
    
        cdef _PeptideHit * _r = new _PeptideHit(self.inst.get().compute((deref(hit.inst.get())), (deref(real_spectrum.inst.get()))))
        cdef PeptideHit py_result = PeptideHit.__new__(PeptideHit)
        py_result.inst = shared_ptr[_PeptideHit](_r)
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

cdef class CVTermListInterface:
    """
    Cython implementation of _CVTermListInterface

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTermListInterface.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVTermListInterface rv = CVTermListInterface.__new__(CVTermListInterface)
       rv.inst = shared_ptr[_CVTermListInterface](new _CVTermListInterface(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVTermListInterface rv = CVTermListInterface.__new__(CVTermListInterface)
       rv.inst = shared_ptr[_CVTermListInterface](new _CVTermListInterface(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVTermListInterface](new _CVTermListInterface())
    
    def _init_1(self, CVTermListInterface in_0 ):
        """
        _init_1(self, in_0: CVTermListInterface ) -> None
        """
        assert isinstance(in_0, CVTermListInterface), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVTermListInterface](new _CVTermListInterface((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVTermListInterface ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVTermListInterface)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _replaceCVTerms_0(self, dict cv_terms ):
        """
        _replaceCVTerms_0(self, cv_terms: Dict[bytes,List[CVTerm]] ) -> None
        """
        assert isinstance(cv_terms, dict) and all(isinstance(k, bytes) for k in cv_terms.keys()) and all(isinstance(v, list) for v in cv_terms.values()) and all(isinstance(vi, CVTerm) for v in cv_terms.values() for vi in
          v), 'arg cv_terms wrong type'
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]] _map_0
        cdef libcpp_vector[_CVTerm] _v_vec_0
        cdef _String _k_str_0
        cdef CVTerm _v_i_0
        for k, v in cv_terms.items():
            _v_vec_0.clear()
            for _v_i_0 in v:
                _v_vec_0.push_back(deref(_v_i_0.inst.get()))
            _map_0[_String(<char *>k)] = _v_vec_0
        self.inst.get().replaceCVTerms(_map_0)
        cdef _replace_0 = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_0 = _map_0.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_0
        cdef CVTerm item_0
        cdef bytes inner_key_0
        cdef list inner_values_0
        while outer_it_0 != _map_0.end():
           inner_key_0 = deref(outer_it_0).first.c_str()
           inner_values_0 = []
           inner_it_0 = deref(outer_it_0).second.begin()
           while inner_it_0 != deref(outer_it_0).second.end():
               item_0 = CVTerm.__new__(CVTerm)
               item_0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_0)))
               inner_values_0.append(item_0)
               inc(inner_it_0)
           _replace_0[inner_key_0] = inner_values_0
           inc(outer_it_0)
        cv_terms.clear()
        cv_terms.update(_replace_0)
    
    def _replaceCVTerms_1(self, list cv_terms ,  accession ):
        """
        _replaceCVTerms_1(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        cdef libcpp_vector[_CVTerm].iterator it_cv_terms = v0.begin()
        replace_0 = []
        while it_cv_terms != v0.end():
            item0 = CVTerm.__new__(CVTerm)
            item0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(it_cv_terms)))
            replace_0.append(item0)
            inc(it_cv_terms)
        cv_terms[:] = replace_0
        del v0
    
    def replaceCVTerms(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: replaceCVTerms(self, cv_terms: Dict[bytes,List[CVTerm]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], dict) and all(isinstance(k, bytes) for k in args[0].keys()) and all(isinstance(v, list) for v in args[0].values()) and all(isinstance(vi, CVTerm) for v in args[0].values() for vi in
          v)):
            return self._replaceCVTerms_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in args[0])) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._replaceCVTerms_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def replaceCVTerm(self, CVTerm cv_term ):
        """
        replaceCVTerm(self, cv_term: CVTerm ) -> None
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().replaceCVTerm((deref(cv_term.inst.get())))
    
    def consumeCVTerms(self, dict cv_term_map ):
        """
        consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None
        Merges the given map into the member map, no duplicate checking
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
        self.inst.get().consumeCVTerms(_map_0)
        cdef _replace_0 = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_0 = _map_0.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_0
        cdef CVTerm item_0
        cdef bytes inner_key_0
        cdef list inner_values_0
        while outer_it_0 != _map_0.end():
           inner_key_0 = deref(outer_it_0).first.c_str()
           inner_values_0 = []
           inner_it_0 = deref(outer_it_0).second.begin()
           while inner_it_0 != deref(outer_it_0).second.end():
               item_0 = CVTerm.__new__(CVTerm)
               item_0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_0)))
               inner_values_0.append(item_0)
               inc(inner_it_0)
           _replace_0[inner_key_0] = inner_values_0
           inc(outer_it_0)
        cv_term_map.clear()
        cv_term_map.update(_replace_0)
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268774070928321753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268774070928321753366645
        cdef CVTerm item_1268774070928321753366645
        cdef bytes inner_key_1268774070928321753366645
        cdef list inner_values_1268774070928321753366645
        while outer_it_1268774070928321753366645 != _r.end():
           inner_key_1268774070928321753366645 = deref(outer_it_1268774070928321753366645).first.c_str()
           inner_values_1268774070928321753366645 = []
           inner_it_1268774070928321753366645 = deref(outer_it_1268774070928321753366645).second.begin()
           while inner_it_1268774070928321753366645 != deref(outer_it_1268774070928321753366645).second.end():
               item_1268774070928321753366645 = CVTerm.__new__(CVTerm)
               item_1268774070928321753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268774070928321753366645)))
               inner_values_1268774070928321753366645.append(item_1268774070928321753366645)
               inc(inner_it_1268774070928321753366645)
           py_result[inner_key_1268774070928321753366645] = inner_values_1268774070928321753366645
           inc(outer_it_1268774070928321753366645)
        return py_result
    
    def addCVTerm(self, CVTerm term ):
        """
        addCVTerm(self, term: CVTerm ) -> None
        Adds a CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().addCVTerm((deref(term.inst.get())))
    
    def hasCVTerm(self,  accession ):
        """
        hasCVTerm(self, accession: Union[bytes, str, String] ) -> bool
        Checks whether the term has a value
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
        if not isinstance(other, CVTermListInterface):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVTermListInterface other_casted = other
        cdef CVTermListInterface self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ChromatogramTools:
    """
    Cython implementation of _ChromatogramTools

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramTools.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramTools rv = ChromatogramTools.__new__(ChromatogramTools)
       rv.inst = shared_ptr[_ChromatogramTools](new _ChromatogramTools(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramTools rv = ChromatogramTools.__new__(ChromatogramTools)
       rv.inst = shared_ptr[_ChromatogramTools](new _ChromatogramTools(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChromatogramTools](new _ChromatogramTools())
    
    def _init_1(self, ChromatogramTools in_0 ):
        """
        _init_1(self, in_0: ChromatogramTools ) -> None
        """
        assert isinstance(in_0, ChromatogramTools), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramTools](new _ChromatogramTools((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramTools ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramTools)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def convertChromatogramsToSpectra(self, MSExperiment epx ):
        """
        convertChromatogramsToSpectra(self, epx: MSExperiment ) -> None
        Converts the chromatogram to a list of spectra with instrument settings
        """
        assert isinstance(epx, MSExperiment), 'arg epx wrong type'
    
        self.inst.get().convertChromatogramsToSpectra((deref(epx.inst.get())))
    
    def convertSpectraToChromatograms(self, MSExperiment epx , bool remove_spectra , bool force_conversion ):
        """
        convertSpectraToChromatograms(self, epx: MSExperiment , remove_spectra: bool , force_conversion: bool ) -> None
        Converts e.g. SRM spectra to chromatograms
        """
        assert isinstance(epx, MSExperiment), 'arg epx wrong type'
        assert isinstance(remove_spectra, pybool_t), 'arg remove_spectra wrong type'
        assert isinstance(force_conversion, pybool_t), 'arg force_conversion wrong type'
    
    
    
        self.inst.get().convertSpectraToChromatograms((deref(epx.inst.get())), (<bool>remove_spectra), (<bool>force_conversion)) 

cdef class EDTAFile:
    """
    Cython implementation of _EDTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EDTAFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EDTAFile rv = EDTAFile.__new__(EDTAFile)
       rv.inst = shared_ptr[_EDTAFile](new _EDTAFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EDTAFile rv = EDTAFile.__new__(EDTAFile)
       rv.inst = shared_ptr[_EDTAFile](new _EDTAFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_EDTAFile](new _EDTAFile())
    
    def _init_1(self, EDTAFile in_0 ):
        """
        _init_1(self, in_0: EDTAFile ) -> None
        """
        assert isinstance(in_0, EDTAFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EDTAFile](new _EDTAFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EDTAFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EDTAFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _store_0(self,  filename , FeatureMap map ):
        """
        _store_0(self, filename: Union[bytes, str, String] , map: FeatureMap ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map, FeatureMap), 'arg map wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(map.inst.get())))
    
    def _store_1(self,  filename , ConsensusMap map ):
        """
        _store_1(self, filename: Union[bytes, str, String] , map: ConsensusMap ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map, ConsensusMap), 'arg map wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(map.inst.get())))
    
    def store(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , map: FeatureMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , map: ConsensusMap ) -> None
          :noindex:
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], FeatureMap)):
            return self._store_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], ConsensusMap)):
            return self._store_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , ConsensusMap consensus_map ):
        """
        load(self, filename: Union[bytes, str, String] , consensus_map: ConsensusMap ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(consensus_map.inst.get()))) 

cdef class FileTypes:
    """
    Cython implementation of _FileTypes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FileTypes.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FileTypes rv = FileTypes.__new__(FileTypes)
       rv.inst = shared_ptr[_FileTypes](new _FileTypes(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FileTypes rv = FileTypes.__new__(FileTypes)
       rv.inst = shared_ptr[_FileTypes](new _FileTypes(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Centralizes the file types recognized by FileHandler
        """
        self.inst = shared_ptr[_FileTypes](new _FileTypes())
    
    def _init_1(self, FileTypes in_0 ):
        """
        _init_1(self, in_0: FileTypes ) -> None
        """
        assert isinstance(in_0, FileTypes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FileTypes](new _FileTypes((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Centralizes the file types recognized by FileHandler

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FileTypes ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FileTypes)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def typeToName(self, int t ):
        """
        typeToName(self, t: int ) -> Union[bytes, str, String]
        Returns the name/extension of the type
        """
        assert t in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg t wrong type'
    
        cdef _String _r = self.inst.get().typeToName((<_FileType>t))
        py_result = convOutputString(_r)
        return py_result
    
    def typeToMZML(self, int t ):
        """
        typeToMZML(self, t: int ) -> Union[bytes, str, String]
        Returns the mzML name
        """
        assert t in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg t wrong type'
    
        cdef _String _r = self.inst.get().typeToMZML((<_FileType>t))
        py_result = convOutputString(_r)
        return py_result
    
    def nameToType(self,  name ):
        """
        nameToType(self, name: Union[bytes, str, String] ) -> int
        Converts a file type name into a Type
        
        
        :param name: A case-insensitive name (e.g. FASTA or Fasta, etc.)
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _FileType _r = self.inst.get().nameToType(deref((convString(name)).get()))
        py_result = <int>_r
        return py_result 

cdef class IDDecoyProbability:
    """
    Cython implementation of _IDDecoyProbability

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDDecoyProbability.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        IDDecoyProbability calculates probabilities using decoy approach
        """
        self.inst = shared_ptr[_IDDecoyProbability](new _IDDecoyProbability())
    
    def _init_1(self, IDDecoyProbability in_0 ):
        """
        _init_1(self, in_0: IDDecoyProbability ) -> None
        """
        assert isinstance(in_0, IDDecoyProbability), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IDDecoyProbability](new _IDDecoyProbability((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        IDDecoyProbability calculates probabilities using decoy approach

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IDDecoyProbability ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IDDecoyProbability)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _apply_0(self, PeptideIdentificationList prob_ids , PeptideIdentificationList fwd_ids , PeptideIdentificationList rev_ids ):
        """
        _apply_0(self, prob_ids: PeptideIdentificationList , fwd_ids: PeptideIdentificationList , rev_ids: PeptideIdentificationList ) -> None
        Converts the forward and reverse identification into probabilities
        
        
        :param prob_ids: Output of the algorithm which includes identifications with probability based scores
        :param fwd_ids: Input parameter which represents the identifications of the forward search
        :param rev_ids: Input parameter which represents the identifications of the reversed search
        """
        assert isinstance(prob_ids, PeptideIdentificationList), 'arg prob_ids wrong type'
        assert isinstance(fwd_ids, PeptideIdentificationList), 'arg fwd_ids wrong type'
        assert isinstance(rev_ids, PeptideIdentificationList), 'arg rev_ids wrong type'
    
    
    
        self.inst.get().apply((deref(prob_ids.inst.get())), (deref(fwd_ids.inst.get())), (deref(rev_ids.inst.get())))
    
    def _apply_1(self, PeptideIdentificationList ids ):
        """
        _apply_1(self, ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
    
        self.inst.get().apply((deref(ids.inst.get())))
    
    def apply(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: apply(self, prob_ids: PeptideIdentificationList , fwd_ids: PeptideIdentificationList , rev_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        Converts the forward and reverse identification into probabilities
        
        
        :param prob_ids: Output of the algorithm which includes identifications with probability based scores
        :param fwd_ids: Input parameter which represents the identifications of the forward search
        :param rev_ids: Input parameter which represents the identifications of the reversed search
        
        .. rubric:: Overload:
        .. py:function:: apply(self, ids: PeptideIdentificationList ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], PeptideIdentificationList)):
            return self._apply_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._apply_1(*args)
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

cdef class KroenikFile:
    """
    Cython implementation of _KroenikFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KroenikFile.html>`_

    File adapter for Kroenik (HardKloer sibling) files
    
    The first line is the header and contains the column names:
    File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
    
    Every subsequent line is a feature
    
    All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
    metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications"
    
    The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_KroenikFile](new _KroenikFile())
    
    def store(self,  filename , MSSpectrum spectrum ):
        """
        store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None
        Stores a MSExperiment into a Kroenik file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(spectrum.inst.get())))
    
    def load(self,  filename , FeatureMap feature_map ):
        """
        load(self, filename: Union[bytes, str, String] , feature_map: FeatureMap ) -> None
        Loads a Kroenik file into a featureXML
        
        The content of the file is stored in `features`
        
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(feature_map, FeatureMap), 'arg feature_map wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(feature_map.inst.get()))) 

cdef class MSPGenericFile:
    """
    Cython implementation of _MSPGenericFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSPGenericFile.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSPGenericFile rv = MSPGenericFile.__new__(MSPGenericFile)
       rv.inst = shared_ptr[_MSPGenericFile](new _MSPGenericFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSPGenericFile rv = MSPGenericFile.__new__(MSPGenericFile)
       rv.inst = shared_ptr[_MSPGenericFile](new _MSPGenericFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSPGenericFile](new _MSPGenericFile())
    
    def _init_1(self, MSPGenericFile in_0 ):
        """
        _init_1(self, in_0: MSPGenericFile ) -> None
        """
        assert isinstance(in_0, MSPGenericFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSPGenericFile](new _MSPGenericFile((deref(in_0.inst.get()))))
    
    def _init_2(self,  filename , MSExperiment library ):
        """
        _init_2(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(library, MSExperiment), 'arg library wrong type'
    
    
        self.inst = shared_ptr[_MSPGenericFile](new _MSPGenericFile(deref((convString(filename)).get()), (deref(library.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSPGenericFile ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSPGenericFile)):
             self._init_1(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], MSExperiment)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment library ):
        """
        load(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None
        Load the file's data and metadata, and save it into an `MSExperiment`
        
        
        :param filename: Path to the MSP input file
        :param library: The variable into which the extracted information will be saved
        :raises:
          Exception: FileNotFound If the file could not be found
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(library, MSExperiment), 'arg library wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(library.inst.get())))
    
    def store(self,  filename , MSExperiment library ):
        """
        store(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None
        Save data and metadata into a file
        
        
        :param filename: Path to the MSP input file
        :param library: The variable from which extracted information will be saved
        :raises:
          Exception: FileNotWritable If the file is not writable
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(library, MSExperiment), 'arg library wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(library.inst.get())))
    
    def getDefaultParameters(self, Param params ):
        """
        getDefaultParameters(self, params: Param ) -> None
        Returns the class' default parameters
        """
        assert isinstance(params, Param), 'arg params wrong type'
    
        self.inst.get().getDefaultParameters((deref(params.inst.get())))
    
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

cdef class PeakPickerIterative:
    """
    Cython implementation of _PeakPickerIterative

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerIterative.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakPickerIterative rv = PeakPickerIterative.__new__(PeakPickerIterative)
       rv.inst = shared_ptr[_PeakPickerIterative](new _PeakPickerIterative(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakPickerIterative rv = PeakPickerIterative.__new__(PeakPickerIterative)
       rv.inst = shared_ptr[_PeakPickerIterative](new _PeakPickerIterative(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakPickerIterative](new _PeakPickerIterative())
    
    def _init_1(self, PeakPickerIterative in_0 ):
        """
        _init_1(self, in_0: PeakPickerIterative ) -> None
        """
        assert isinstance(in_0, PeakPickerIterative), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakPickerIterative](new _PeakPickerIterative((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakPickerIterative ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakPickerIterative)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def pick(self, MSSpectrum input , MSSpectrum output ):
        """
        pick(self, input: MSSpectrum , output: MSSpectrum ) -> None
        This will pick one single spectrum. The PeakPickerHiRes is used to
        generate seeds, these seeds are then used to re-center the mass and
        compute peak width and integrated intensity of the peak
        
        Finally, other peaks that would fall within the primary peak are
        discarded
        
        The output are the remaining peaks
        """
        assert isinstance(input, MSSpectrum), 'arg input wrong type'
        assert isinstance(output, MSSpectrum), 'arg output wrong type'
    
    
        self.inst.get().pick((deref(input.inst.get())), (deref(output.inst.get())))
    
    def pickExperiment(self, MSExperiment input , MSExperiment output ):
        """
        pickExperiment(self, input: MSExperiment , output: MSExperiment ) -> None
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
        assert isinstance(output, MSExperiment), 'arg output wrong type'
    
    
        self.inst.get().pickExperiment((deref(input.inst.get())), (deref(output.inst.get())))
    
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

cdef class PeptideIdentificationList:
    """
    Cython implementation of _PeptideIdentificationList

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIdentificationList.html>`_

    A container for peptide identifications from multiple spectra.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideIdentificationList](new _PeptideIdentificationList())
    
    def _init_1(self, PeptideIdentificationList in_0 ):
        """
        _init_1(self, in_0: PeptideIdentificationList ) -> None
        """
        assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideIdentificationList](new _PeptideIdentificationList((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideIdentificationList ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        Returns the number of peptide identifications
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        Returns true if the container is empty
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        Removes all peptide identifications from the container
        """
        self.inst.get().clear()
    
    def push_back(self, PeptideIdentification in_0 ):
        """
        push_back(self, in_0: PeptideIdentification ) -> None
        Adds a peptide identification to the end of the container
        """
        assert isinstance(in_0, PeptideIdentification), 'arg in_0 wrong type'
    
        self.inst.get().push_back((deref(in_0.inst.get())))
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> PeptideIdentification
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _PeptideIdentification * _r = new _PeptideIdentification(deref(self.inst.get())[(<size_t>in_0)])
        cdef PeptideIdentification py_result = PeptideIdentification.__new__(PeptideIdentification)
        py_result.inst = shared_ptr[_PeptideIdentification](_r)
        return py_result
    def __setitem__(self, key, PeptideIdentification value):
        """Cython signature: PeptideIdentification & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def at(self,  in_0 ):
        """
        at(self, in_0: int ) -> PeptideIdentification
        Returns the peptide identification at the given index with bounds checking
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef _PeptideIdentification * _r = new _PeptideIdentification(self.inst.get().at((<size_t>in_0)))
        cdef PeptideIdentification py_result = PeptideIdentification.__new__(PeptideIdentification)
        py_result.inst = shared_ptr[_PeptideIdentification](_r)
        return py_result
    
    def back(self):
        """
        back(self) -> PeptideIdentification
        Returns the last peptide identification in the container
        """
        cdef _PeptideIdentification * _r = new _PeptideIdentification(self.inst.get().back())
        cdef PeptideIdentification py_result = PeptideIdentification.__new__(PeptideIdentification)
        py_result.inst = shared_ptr[_PeptideIdentification](_r)
        return py_result
    
    def front(self):
        """
        front(self) -> PeptideIdentification
        Returns the first peptide identification in the container
        """
        cdef _PeptideIdentification * _r = new _PeptideIdentification(self.inst.get().front())
        cdef PeptideIdentification py_result = PeptideIdentification.__new__(PeptideIdentification)
        py_result.inst = shared_ptr[_PeptideIdentification](_r)
        return py_result
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef PeptideIdentification out
        while it != self.inst.get().end():
            out = PeptideIdentification.__new__(PeptideIdentification)
            out.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification(deref(it)))
            yield out
            inc(it) 

cdef class ProbablePhosphoSites:
    """
    Cython implementation of _ProbablePhosphoSites

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProbablePhosphoSites.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property first:
        def __set__(self,  first):
        
            self.inst.get().first = (<size_t>first)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().first
            py_result = <size_t>_r
            return py_result
    
    property second:
        def __set__(self,  second):
        
            self.inst.get().second = (<size_t>second)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().second
            py_result = <size_t>_r
            return py_result
    
    property seq_1:
        def __set__(self,  seq_1):
        
            self.inst.get().seq_1 = (<size_t>seq_1)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().seq_1
            py_result = <size_t>_r
            return py_result
    
    property seq_2:
        def __set__(self,  seq_2):
        
            self.inst.get().seq_2 = (<size_t>seq_2)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().seq_2
            py_result = <size_t>_r
            return py_result
    
    property peak_depth:
        def __set__(self,  peak_depth):
        
            self.inst.get().peak_depth = (<size_t>peak_depth)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().peak_depth
            py_result = <size_t>_r
            return py_result
    
    property AScore:
        def __set__(self,  AScore):
        
            self.inst.get().AScore = (<size_t>AScore)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().AScore
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef ProbablePhosphoSites rv = ProbablePhosphoSites.__new__(ProbablePhosphoSites)
       rv.inst = shared_ptr[_ProbablePhosphoSites](new _ProbablePhosphoSites(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProbablePhosphoSites rv = ProbablePhosphoSites.__new__(ProbablePhosphoSites)
       rv.inst = shared_ptr[_ProbablePhosphoSites](new _ProbablePhosphoSites(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProbablePhosphoSites](new _ProbablePhosphoSites())
    
    def _init_1(self, ProbablePhosphoSites in_0 ):
        """
        _init_1(self, in_0: ProbablePhosphoSites ) -> None
        """
        assert isinstance(in_0, ProbablePhosphoSites), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProbablePhosphoSites](new _ProbablePhosphoSites((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProbablePhosphoSites ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProbablePhosphoSites)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class RealMassDecomposer:
    """
    Cython implementation of _RealMassDecomposer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims_1_1RealMassDecomposer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self, RealMassDecomposer in_0 ):
        """
        _init_0(self, in_0: RealMassDecomposer ) -> None
        """
        assert isinstance(in_0, RealMassDecomposer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RealMassDecomposer](new _RealMassDecomposer((deref(in_0.inst.get()))))
    
    def _init_1(self, IMSWeights weights ):
        """
        _init_1(self, weights: IMSWeights ) -> None
        """
        assert isinstance(weights, IMSWeights), 'arg weights wrong type'
    
        self.inst = shared_ptr[_RealMassDecomposer](new _RealMassDecomposer((deref(weights.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RealMassDecomposer ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, weights: IMSWeights ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], RealMassDecomposer)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSWeights)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getNumberOfDecompositions(self, double mass , double error ):
        """
        getNumberOfDecompositions(self, mass: float , error: float ) -> int
        Gets a number of all decompositions for amass with an error
        allowed. It's similar to thegetDecompositions(double,double) function
        but less space consuming, since doesn't use container to store decompositions
        
        
        :param mass: Mass to be decomposed
        :param error: Error allowed between given and result decomposition
        :return: Number of all decompositions for a given mass and error
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(error, float), 'arg error wrong type'
    
    
        cdef uint64_t _r = self.inst.get().getNumberOfDecompositions((<double>mass), (<double>error))
        py_result = <uint64_t>_r
        return py_result 

cdef class SemanticValidator:
    """
    Cython implementation of _SemanticValidator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1SemanticValidator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, CVMappings mapping , ControlledVocabulary cv ):
        """
        __init__(self, mapping: CVMappings , cv: ControlledVocabulary ) -> None
        """
        assert isinstance(mapping, CVMappings), 'arg mapping wrong type'
        assert isinstance(cv, ControlledVocabulary), 'arg cv wrong type'
    
    
        self.inst = shared_ptr[_SemanticValidator](new _SemanticValidator((deref(mapping.inst.get())), (deref(cv.inst.get()))))
    
    def validate(self,  filename , list errors , list warnings ):
        """
        validate(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool
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
        cdef bool _r = self.inst.get().validate(deref((convString(filename)).get()), deref(v1), deref(v2))
        del v2
        del v1
        py_result = <bool>_r
        return py_result
    
    def locateTerm(self,  path , SemanticValidator_CVTerm parsed_term ):
        """
        locateTerm(self, path: Union[bytes, str, String] , parsed_term: SemanticValidator_CVTerm ) -> bool
        Checks if a CVTerm is allowed in a given path
        """
        assert (isinstance(path, str) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'
        assert isinstance(parsed_term, SemanticValidator_CVTerm), 'arg parsed_term wrong type'
    
    
        cdef bool _r = self.inst.get().locateTerm(deref((convString(path)).get()), (deref(parsed_term.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def setTag(self,  tag ):
        """
        setTag(self, tag: Union[bytes, str, String] ) -> None
        Sets the CV parameter tag name (default 'cvParam')
        """
        assert (isinstance(tag, str) or isinstance(tag, bytes) or isinstance(tag, String)), 'arg tag wrong type'
    
        self.inst.get().setTag(deref((convString(tag)).get()))
    
    def setAccessionAttribute(self,  accession ):
        """
        setAccessionAttribute(self, accession: Union[bytes, str, String] ) -> None
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'accession')
        """
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
    
        self.inst.get().setAccessionAttribute(deref((convString(accession)).get()))
    
    def setNameAttribute(self,  name ):
        """
        setNameAttribute(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'name')
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setNameAttribute(deref((convString(name)).get()))
    
    def setValueAttribute(self,  value ):
        """
        setValueAttribute(self, value: Union[bytes, str, String] ) -> None
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'value')
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setValueAttribute(deref((convString(value)).get()))
    
    def setCheckTermValueTypes(self, bool check ):
        """
        setCheckTermValueTypes(self, check: bool ) -> None
        Sets if CV term value types should be check (enabled by default)
        """
        assert isinstance(check, pybool_t), 'arg check wrong type'
    
        self.inst.get().setCheckTermValueTypes((<bool>check))
    
    def setCheckUnits(self, bool check ):
        """
        setCheckUnits(self, check: bool ) -> None
        Sets if CV term units should be check (disabled by default)
        """
        assert isinstance(check, pybool_t), 'arg check wrong type'
    
        self.inst.get().setCheckUnits((<bool>check))
    
    def setUnitAccessionAttribute(self,  accession ):
        """
        setUnitAccessionAttribute(self, accession: Union[bytes, str, String] ) -> None
        Sets the name of the unit accession attribute (default 'unitAccession')
        """
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
    
        self.inst.get().setUnitAccessionAttribute(deref((convString(accession)).get()))
    
    def setUnitNameAttribute(self,  name ):
        """
        setUnitNameAttribute(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the unit name attribute (default 'unitName')
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setUnitNameAttribute(deref((convString(name)).get())) 

cdef class SemanticValidator_CVTerm:
    """
    Cython implementation of _SemanticValidator_CVTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1SemanticValidator_CVTerm.html>`_
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
    
    property value:
        def __set__(self,  value):
        
            self.inst.get().value = deref((convString(value)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().value
            py_result = convOutputString(_r)
            return py_result
    
    property has_value:
        def __set__(self, bool has_value):
        
            self.inst.get().has_value = (<bool>has_value)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().has_value
            py_result = <bool>_r
            return py_result
    
    property unit_accession:
        def __set__(self,  unit_accession):
        
            self.inst.get().unit_accession = deref((convString(unit_accession)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().unit_accession
            py_result = convOutputString(_r)
            return py_result
    
    property has_unit_accession:
        def __set__(self, bool has_unit_accession):
        
            self.inst.get().has_unit_accession = (<bool>has_unit_accession)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().has_unit_accession
            py_result = <bool>_r
            return py_result
    
    property unit_name:
        def __set__(self,  unit_name):
        
            self.inst.get().unit_name = deref((convString(unit_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().unit_name
            py_result = convOutputString(_r)
            return py_result
    
    property has_unit_name:
        def __set__(self, bool has_unit_name):
        
            self.inst.get().has_unit_name = (<bool>has_unit_name)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().has_unit_name
            py_result = <bool>_r
            return py_result
    
    def __copy__(self):
       cdef SemanticValidator_CVTerm rv = SemanticValidator_CVTerm.__new__(SemanticValidator_CVTerm)
       rv.inst = shared_ptr[_SemanticValidator_CVTerm](new _SemanticValidator_CVTerm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SemanticValidator_CVTerm rv = SemanticValidator_CVTerm.__new__(SemanticValidator_CVTerm)
       rv.inst = shared_ptr[_SemanticValidator_CVTerm](new _SemanticValidator_CVTerm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SemanticValidator_CVTerm](new _SemanticValidator_CVTerm())
    
    def _init_1(self, SemanticValidator_CVTerm in_0 ):
        """
        _init_1(self, in_0: SemanticValidator_CVTerm ) -> None
        """
        assert isinstance(in_0, SemanticValidator_CVTerm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SemanticValidator_CVTerm](new _SemanticValidator_CVTerm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SemanticValidator_CVTerm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SemanticValidator_CVTerm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SpectrumAlignment:
    """
    Cython implementation of _SpectrumAlignment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAlignment.html>`_
      -- Inherits from ['DefaultParamHandler']

    Aligns the peaks of two sorted spectra
    Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given
        Scoring function is the m/z distance between peaks. Intensity does not play a role!
    Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
    Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window
    
    note: A peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2
    note: Intensity is ignored
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAlignment rv = SpectrumAlignment.__new__(SpectrumAlignment)
       rv.inst = shared_ptr[_SpectrumAlignment](new _SpectrumAlignment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAlignment rv = SpectrumAlignment.__new__(SpectrumAlignment)
       rv.inst = shared_ptr[_SpectrumAlignment](new _SpectrumAlignment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectrumAlignment](new _SpectrumAlignment())
    
    def _init_1(self, SpectrumAlignment in_0 ):
        """
        _init_1(self, in_0: SpectrumAlignment ) -> None
        """
        assert isinstance(in_0, SpectrumAlignment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAlignment](new _SpectrumAlignment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAlignment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAlignment)):
             self._init_1(*args)
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
    
    def getSpectrumAlignment(self, list result, spec1, spec2):
        # TODO: use autowrap for this directly
        assert isinstance(spec1, (MSSpectrum))
        assert isinstance(spec2, (MSSpectrum))
        cdef _MSSpectrum * _new_spec1 = new _MSSpectrum()
        cdef _MSSpectrum * _new_spec2 = new _MSSpectrum()
        cdef _Peak1D _peak

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec1).inst.get()):
                _new_spec1.push_back(<_Peak1D>_peak)

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec2).inst.get()):
                _new_spec2.push_back(<_Peak1D>_peak)

        cdef libcpp_vector[libcpp_pair[Size, Size]] _result

        self.inst.get().getSpectrumAlignment(_result, deref(_new_spec1), deref(_new_spec2))

        result[:] = []

        cdef libcpp_vector[libcpp_pair[Size, Size]].iterator _it =  _result.begin()
        while _it != _result.end():
            result.append((<int>(deref(_it).first), <int>(deref(_it).second)))
            inc(_it)

        del _new_spec1
        del _new_spec2 

cdef class SpectrumAnnotator:
    """
    Cython implementation of _SpectrumAnnotator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAnnotator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAnnotator rv = SpectrumAnnotator.__new__(SpectrumAnnotator)
       rv.inst = shared_ptr[_SpectrumAnnotator](new _SpectrumAnnotator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAnnotator rv = SpectrumAnnotator.__new__(SpectrumAnnotator)
       rv.inst = shared_ptr[_SpectrumAnnotator](new _SpectrumAnnotator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Annotates spectra from identifications and theoretical spectra or
        identifications from spectra and theoretical spectra matching
        with various options
        """
        self.inst = shared_ptr[_SpectrumAnnotator](new _SpectrumAnnotator())
    
    def _init_1(self, SpectrumAnnotator in_0 ):
        """
        _init_1(self, in_0: SpectrumAnnotator ) -> None
        """
        assert isinstance(in_0, SpectrumAnnotator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAnnotator](new _SpectrumAnnotator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Annotates spectra from identifications and theoretical spectra or
        identifications from spectra and theoretical spectra matching
        with various options
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAnnotator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAnnotator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def annotateMatches(self, MSSpectrum spec , PeptideHit ph , TheoreticalSpectrumGenerator tg , SpectrumAlignment sa ):
        """
        annotateMatches(self, spec: MSSpectrum , ph: PeptideHit , tg: TheoreticalSpectrumGenerator , sa: SpectrumAlignment ) -> None
        Adds ion match annotation to the `spec` input spectrum
        
        :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
        :param ph: A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
        :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
        :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(ph, PeptideHit), 'arg ph wrong type'
        assert isinstance(tg, TheoreticalSpectrumGenerator), 'arg tg wrong type'
        assert isinstance(sa, SpectrumAlignment), 'arg sa wrong type'
    
    
    
    
        self.inst.get().annotateMatches((deref(spec.inst.get())), (deref(ph.inst.get())), (deref(tg.inst.get())), (deref(sa.inst.get())))
    
    def addIonMatchStatistics(self, PeptideIdentification pi , MSSpectrum spec , TheoreticalSpectrumGenerator tg , SpectrumAlignment sa ):
        """
        addIonMatchStatistics(self, pi: PeptideIdentification , spec: MSSpectrum , tg: TheoreticalSpectrumGenerator , sa: SpectrumAlignment ) -> None
        Adds ion match statistics to `pi` PeptideIdentifcation
        
        :param pi: A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
        :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
        :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
        :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
        """
        assert isinstance(pi, PeptideIdentification), 'arg pi wrong type'
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(tg, TheoreticalSpectrumGenerator), 'arg tg wrong type'
        assert isinstance(sa, SpectrumAlignment), 'arg sa wrong type'
    
    
    
    
        self.inst.get().addIonMatchStatistics((deref(pi.inst.get())), (deref(spec.inst.get())), (deref(tg.inst.get())), (deref(sa.inst.get())))
    
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


cdef inline convOutputString(_String _r):
    # Generic function to convert OpenMS::String to Python unicode strings
    #
    # This assumes that all data encoded in OpenMS::String is encoded using
    # UTF8 (in practice OpenMS does not make any assumption about encoding and
    # OpenMS::String simply contains a byte sequence, that is implicitely
    # assumed to be ASCII encoded).
    #
    cdef char* c_string = _cast_const_away(<char*> _r.c_str())
    cdef Py_ssize_t length = _r.length()

    try:
        py_result = c_string[:length].decode('UTF-8')
    except UnicodeDecodeError:
        py_result = c_string[:length]

    return py_result

cdef inline shared_ptr[_String] convString(argument_var):
    # Generic function to convert Python strings to OpenMS::String
    # 
    # This allows us to either extract the already existing shared_ptr from a
    # Python-type OpenMS::String holder or create a new shared_ptr from a
    # Python string.
    # In case the user only provides a str, unicode or bytes argument, we will
    # have to create a new OpenMS::String object, the shared_ptr will
    # automatically delete it after it has been used for a function call
    # (presumably the user does not want to have the value returned by
    # reference).
    #
    # Note: even though Python 3.x does not know the unicode keyword,
    # Cython does and uses the PyUnicode_Check call. 
    #
    cdef char* c_string_argument_var
    cdef shared_ptr[_String] res
    if isinstance(argument_var, String):
        res = (<String>argument_var).inst
        return res
    elif isinstance(argument_var, bytes):
        # Simple, convert directly to char* - however there may be zero bytes
        # in the array when using certain encodings, we therefore need to make
        # sure we capture the full length.
        res = shared_ptr[_String](new _String(<char*>argument_var, len(argument_var)))
        return res
    elif isinstance(argument_var, str) or isinstance(argument_var, unicode):
        # First encode with UTF8 (note that output encoding as well as
        # String.toString both decode with UTF8), then convert the result to
        # char*
        py_byte_string = argument_var.encode('UTF-8')
        c_string_argument_var = py_byte_string
        res = shared_ptr[_String](new _String(<char*>c_string_argument_var))
        return res
    else:
        raise Exception("Can only convert the following types to String: pyopenms.String, bytes, str, unicode") 
 

## Add this to all modules
class Interfaces:

    BinaryDataArray = _Interfaces_BinaryDataArray
    Spectrum = _Interfaces_Spectrum
    Chromatogram = _Interfaces_Chromatogram 

# Add this to all modules
PeakSpectrum = MSSpectrum

PeakMap = MSExperiment 
