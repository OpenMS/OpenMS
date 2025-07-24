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
def __static_ExperimentalDesignFile_load( tsv_file , bool in_1 ):
    """
    __static_ExperimentalDesignFile_load(tsv_file: Union[bytes, str, String] , in_1: bool ) -> ExperimentalDesign
    """
    assert (isinstance(tsv_file, str) or isinstance(tsv_file, bytes) or isinstance(tsv_file, String)), 'arg tsv_file wrong type'
    assert isinstance(in_1, pybool_t), 'arg in_1 wrong type'


    cdef _ExperimentalDesign * _r = new _ExperimentalDesign(_load_ExperimentalDesignFile(deref((convString(tsv_file)).get()), (<bool>in_1)))
    cdef ExperimentalDesign py_result = ExperimentalDesign.__new__(ExperimentalDesign)
    py_result.inst = shared_ptr[_ExperimentalDesign](_r)
    return py_result 

cdef class DecoyTransitionType:
    None
    UNKNOWN = 0
    TARGET = 1
    DECOY = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __RTType:
    None
    LOCAL = 0
    NORMALIZED = 1
    PREDICTED = 2
    HPINS = 3
    IRT = 4
    UNKNOWN = 5
    SIZE_OF_RTTYPE = 6

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __RTUnit:
    None
    SECOND = 0
    MINUTE = 1
    UNKNOWN = 2
    SIZE_OF_RTUNIT = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class TransType:
    None
    PEPTIDE = 0
    NUCTIDE = 1
    COMPOUND = 2
    UNKNOWN = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AASeqWithMass:
    """
    Cython implementation of _AASeqWithMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AASeqWithMass.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property peptide_mass:
        def __set__(self, double peptide_mass):
        
            self.inst.get().peptide_mass = (<double>peptide_mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().peptide_mass
            py_result = <double>_r
            return py_result
    
    property peptide_seq:
        def __set__(self, AASequence peptide_seq):
        
            self.inst.get().peptide_seq = (deref(peptide_seq.inst.get()))
        
    
        def __get__(self):
            cdef _AASequence * _r = new _AASequence(self.inst.get().peptide_seq)
            cdef AASequence py_result = AASequence.__new__(AASequence)
            py_result.inst = shared_ptr[_AASequence](_r)
            return py_result
    
    property position:
        def __set__(self, int position):
        
            self.inst.get().position = (<_PeptidePosition>position)
        
    
        def __get__(self):
            cdef _PeptidePosition _r = self.inst.get().position
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef AASeqWithMass rv = AASeqWithMass.__new__(AASeqWithMass)
       rv.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AASeqWithMass rv = AASeqWithMass.__new__(AASeqWithMass)
       rv.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass())
    
    def _init_1(self, AASeqWithMass in_0 ):
        """
        _init_1(self, in_0: AASeqWithMass ) -> None
        """
        assert isinstance(in_0, AASeqWithMass), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AASeqWithMass ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AASeqWithMass)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class Acquisition:
    """
    Cython implementation of _Acquisition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Acquisition.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Acquisition rv = Acquisition.__new__(Acquisition)
       rv.inst = shared_ptr[_Acquisition](new _Acquisition(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Acquisition rv = Acquisition.__new__(Acquisition)
       rv.inst = shared_ptr[_Acquisition](new _Acquisition(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Acquisition](new _Acquisition())
    
    def _init_1(self, Acquisition in_0 ):
        """
        _init_1(self, in_0: Acquisition ) -> None
        """
        assert isinstance(in_0, Acquisition), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Acquisition](new _Acquisition((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Acquisition ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Acquisition)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setIdentifier(self,  identifier ):
        """
        setIdentifier(self, identifier: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
    
        self.inst.get().setIdentifier(deref((convString(identifier)).get()))
    
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
        if not isinstance(other, Acquisition):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Acquisition other_casted = other
        cdef Acquisition self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class CV:
    """
    Cython implementation of _CV

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1CV.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property fullname:
        def __set__(self,  fullname):
        
            self.inst.get().fullname = deref((convString(fullname)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().fullname
            py_result = convOutputString(_r)
            return py_result
    
    property version:
        def __set__(self,  version):
        
            self.inst.get().version = deref((convString(version)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().version
            py_result = convOutputString(_r)
            return py_result
    
    property URI:
        def __set__(self,  URI):
        
            self.inst.get().URI = deref((convString(URI)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().URI
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef CV rv = CV.__new__(CV)
       rv.inst = shared_ptr[_CV](new _CV(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CV rv = CV.__new__(CV)
       rv.inst = shared_ptr[_CV](new _CV(deref(self.inst.get())))
       return rv
    
    def _init_0(self, CV in_0 ):
        """
        _init_0(self, in_0: CV ) -> None
        """
        assert isinstance(in_0, CV), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CV](new _CV((deref(in_0.inst.get()))))
    
    def _init_1(self,  new_id ,  new_fullname ,  new_version ,  new_URI ):
        """
        _init_1(self, new_id: Union[bytes, str, String] , new_fullname: Union[bytes, str, String] , new_version: Union[bytes, str, String] , new_URI: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(new_id, str) or isinstance(new_id, bytes) or isinstance(new_id, String)), 'arg new_id wrong type'
        assert (isinstance(new_fullname, str) or isinstance(new_fullname, bytes) or isinstance(new_fullname, String)), 'arg new_fullname wrong type'
        assert (isinstance(new_version, str) or isinstance(new_version, bytes) or isinstance(new_version, String)), 'arg new_version wrong type'
        assert (isinstance(new_URI, str) or isinstance(new_URI, bytes) or isinstance(new_URI, String)), 'arg new_URI wrong type'
    
    
    
    
        self.inst = shared_ptr[_CV](new _CV(deref((convString(new_id)).get()), deref((convString(new_fullname)).get()), deref((convString(new_version)).get()), deref((convString(new_URI)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CV ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, new_id: Union[bytes, str, String] , new_fullname: Union[bytes, str, String] , new_version: Union[bytes, str, String] , new_URI: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], CV)):
             self._init_0(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, CV):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CV other_casted = other
        cdef CV self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Compound:
    """
    Cython implementation of _Compound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Compound.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property molecular_formula:
        def __set__(self,  molecular_formula):
        
            self.inst.get().molecular_formula = deref((convString(molecular_formula)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().molecular_formula
            py_result = convOutputString(_r)
            return py_result
    
    property smiles_string:
        def __set__(self,  smiles_string):
        
            self.inst.get().smiles_string = deref((convString(smiles_string)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().smiles_string
            py_result = convOutputString(_r)
            return py_result
    
    property theoretical_mass:
        def __set__(self, double theoretical_mass):
        
            self.inst.get().theoretical_mass = (<double>theoretical_mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().theoretical_mass
            py_result = <double>_r
            return py_result
    
    property rts:
        def __set__(self, list rts):
            cdef libcpp_vector[_RetentionTime] * v0 = new libcpp_vector[_RetentionTime]()
            cdef RetentionTime item0
            for item0 in rts:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().rts = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().rts
            py_result = []
            cdef libcpp_vector[_RetentionTime].iterator it__r = _r.begin()
            cdef RetentionTime item_py_result
            while it__r != _r.end():
               item_py_result = RetentionTime.__new__(RetentionTime)
               item_py_result.inst = shared_ptr[_RetentionTime](new _RetentionTime(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef Compound rv = Compound.__new__(Compound)
       rv.inst = shared_ptr[_Compound](new _Compound(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Compound rv = Compound.__new__(Compound)
       rv.inst = shared_ptr[_Compound](new _Compound(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Compound](new _Compound())
    
    def _init_1(self, Compound in_0 ):
        """
        _init_1(self, in_0: Compound ) -> None
        """
        assert isinstance(in_0, Compound), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Compound](new _Compound((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Compound ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Compound)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setChargeState(self,  charge ):
        """
        setChargeState(self, charge: int ) -> None
        Sets the peptide or compound charge state
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setChargeState((<int>charge))
    
    def getChargeState(self):
        """
        getChargeState(self) -> int
        Returns the peptide or compound charge state
        """
        cdef int _r = self.inst.get().getChargeState()
        py_result = <int>_r
        return py_result
    
    def hasCharge(self):
        """
        hasCharge(self) -> bool
        Whether peptide or compound has set charge state
        """
        cdef bool _r = self.inst.get().hasCharge()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTime(self):
        """
        getRetentionTime(self) -> float
        Gets compound or peptide retention time
        """
        cdef double _r = self.inst.get().getRetentionTime()
        py_result = <double>_r
        return py_result
    
    def hasRetentionTime(self):
        """
        hasRetentionTime(self) -> bool
        Check whether compound or peptide has an annotated retention time
        """
        cdef bool _r = self.inst.get().hasRetentionTime()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTimeType(self):
        """
        getRetentionTimeType(self) -> int
        Get compound or peptide retentiontime type
        """
        cdef _RTType _r = self.inst.get().getRetentionTimeType()
        py_result = <int>_r
        return py_result
    
    def getRetentionTimeUnit(self):
        """
        getRetentionTimeUnit(self) -> int
        Get compound or peptide retentiontime type
        """
        cdef _RTUnit _r = self.inst.get().getRetentionTimeUnit()
        py_result = <int>_r
        return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Compound):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Compound other_casted = other
        cdef Compound self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Configuration:
    """
    Cython implementation of _Configuration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Configuration.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property contact_ref:
        def __set__(self,  contact_ref):
        
            self.inst.get().contact_ref = deref((convString(contact_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().contact_ref
            py_result = convOutputString(_r)
            return py_result
    
    property instrument_ref:
        def __set__(self,  instrument_ref):
        
            self.inst.get().instrument_ref = deref((convString(instrument_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().instrument_ref
            py_result = convOutputString(_r)
            return py_result
    
    property validations:
        def __set__(self, list validations):
            cdef libcpp_vector[_CVTermList] * v0 = new libcpp_vector[_CVTermList]()
            cdef CVTermList item0
            for item0 in validations:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().validations = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().validations
            py_result = []
            cdef libcpp_vector[_CVTermList].iterator it__r = _r.begin()
            cdef CVTermList item_py_result
            while it__r != _r.end():
               item_py_result = CVTermList.__new__(CVTermList)
               item_py_result.inst = shared_ptr[_CVTermList](new _CVTermList(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Configuration):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Configuration other_casted = other
        cdef Configuration self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Contact:
    """
    Cython implementation of _Contact

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Contact.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef Contact rv = Contact.__new__(Contact)
       rv.inst = shared_ptr[_Contact](new _Contact(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Contact rv = Contact.__new__(Contact)
       rv.inst = shared_ptr[_Contact](new _Contact(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Contact](new _Contact())
    
    def _init_1(self, Contact in_0 ):
        """
        _init_1(self, in_0: Contact ) -> None
        """
        assert isinstance(in_0, Contact), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Contact](new _Contact((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Contact ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Contact)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Contact):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Contact other_casted = other
        cdef Contact self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class CrossLinksDB:
    """
    Cython implementation of _CrossLinksDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CrossLinksDB.html>`_
    """

    
    def getNumberOfModifications(self):
        """
        getNumberOfModifications(self) -> int
        """
        cdef size_t _r = self.inst.get().getNumberOfModifications()
        py_result = <size_t>_r
        return py_result
    
    def searchModifications(self, set mods ,  mod_name ,  residue , int term_spec ):
        """
        searchModifications(self, mods: Set[ResidueModification] , mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> None
        """
        assert isinstance(mods, set) and all(isinstance(li, ResidueModification) for li in mods), 'arg mods wrong type'
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
        cdef libcpp_set[const _ResidueModification *] * v0 = new libcpp_set[const _ResidueModification *]()
        cdef ResidueModification item0
        for item0 in mods:
           v0.insert((item0.inst.get()))
    
    
    
        self.inst.get().searchModifications(deref(v0), deref((convString(mod_name)).get()), deref((convString(residue)).get()), (<_TermSpecificity>term_spec))
        replace = set()
        cdef libcpp_set[const _ResidueModification *].iterator it_mods = v0.begin()
        while it_mods != v0.end():
           item0 = ResidueModification.__new__(ResidueModification)
           item0.inst = shared_ptr[_ResidueModification](new _ResidueModification(deref(deref(it_mods))))
           replace.add(item0)
           inc(it_mods)
        mods.clear()
        mods.update(replace)
        del v0
    
    def _getModification_0(self,  index ):
        """
        _getModification_0(self, index: int ) -> ResidueModification
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification((<size_t>index)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def _getModification_1(self,  mod_name ):
        """
        _getModification_1(self, mod_name: Union[bytes, str, String] ) -> ResidueModification
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification(deref((convString(mod_name)).get())))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def _getModification_2(self,  mod_name ,  residue , int term_spec ):
        """
        _getModification_2(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
    
    
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification(deref((convString(mod_name)).get()), deref((convString(residue)).get()), (<_TermSpecificity>term_spec)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getModification(self, index: int ) -> ResidueModification
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getModification(self, mod_name: Union[bytes, str, String] ) -> ResidueModification
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getModification(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
            return self._getModification_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getModification_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (args[2] in [0, 1, 2, 3, 4, 5]):
            return self._getModification_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def has(self,  modification ):
        """
        has(self, modification: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
        cdef bool _r = self.inst.get().has(deref((convString(modification)).get()))
        py_result = <bool>_r
        return py_result
    
    def findModificationIndex(self,  mod_name ):
        """
        findModificationIndex(self, mod_name: Union[bytes, str, String] ) -> int
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
    
        cdef size_t _r = self.inst.get().findModificationIndex(deref((convString(mod_name)).get()))
        py_result = <size_t>_r
        return py_result
    
    def searchModificationsByDiffMonoMass(self, list mods , double mass , double max_error ,  residue , int term_spec ):
        """
        searchModificationsByDiffMonoMass(self, mods: List[bytes] , mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> None
        """
        assert isinstance(mods, list) and all(isinstance(i, bytes) for i in mods), 'arg mods wrong type'
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(max_error, float), 'arg max_error wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in mods:
           v0.push_back(_String(<char *>item0))
    
    
    
    
        self.inst.get().searchModificationsByDiffMonoMass(deref(v0), (<double>mass), (<double>max_error), deref((convString(residue)).get()), (<_TermSpecificity>term_spec))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        mods[:] = replace
        del v0
    
    def getBestModificationByDiffMonoMass(self, double mass , double max_error ,  residue , int term_spec ):
        """
        getBestModificationByDiffMonoMass(self, mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(max_error, float), 'arg max_error wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
    
    
    
    
        cdef const _ResidueModification * __r = (self.inst.get().getBestModificationByDiffMonoMass((<double>mass), (<double>max_error), deref((convString(residue)).get()), (<_TermSpecificity>term_spec)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getAllSearchModifications(self, list modifications ):
        """
        getAllSearchModifications(self, modifications: List[bytes] ) -> None
        """
        assert isinstance(modifications, list) and all(isinstance(i, bytes) for i in modifications), 'arg modifications wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in modifications:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllSearchModifications(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        modifications[:] = replace
        del v0
    
    def readFromOBOFile(self,  filename ):
        """
        readFromOBOFile(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().readFromOBOFile(deref((convString(filename)).get()))
    
    def isInstantiated(self):
        """
        isInstantiated(self) -> bool
        """
        cdef bool _r = self.inst.get().isInstantiated()
        py_result = <bool>_r
        return py_result
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # Provides the same interface as expected by autowrap (get/assign) but
    # simply holds the ptr and does not do any memory management.
    # see autowrap/data_files/autowrap/README.md for implementation
    # cdef AutowrapPtrHolder[_CrossLinksDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_CrossLinksDB](_getInstance_CrossLinksDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class ExperimentalDesignFile:
    """
    Cython implementation of _ExperimentalDesignFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesignFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ExperimentalDesignFile rv = ExperimentalDesignFile.__new__(ExperimentalDesignFile)
       rv.inst = shared_ptr[_ExperimentalDesignFile](new _ExperimentalDesignFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ExperimentalDesignFile rv = ExperimentalDesignFile.__new__(ExperimentalDesignFile)
       rv.inst = shared_ptr[_ExperimentalDesignFile](new _ExperimentalDesignFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ExperimentalDesignFile](new _ExperimentalDesignFile())
    
    def _init_1(self, ExperimentalDesignFile in_0 ):
        """
        _init_1(self, in_0: ExperimentalDesignFile ) -> None
        """
        assert isinstance(in_0, ExperimentalDesignFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ExperimentalDesignFile](new _ExperimentalDesignFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ExperimentalDesignFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ExperimentalDesignFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    load = __static_ExperimentalDesignFile_load 

cdef class MRMTransitionGroupPicker:
    """
    Cython implementation of _MRMTransitionGroupPicker

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroupPicker.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMTransitionGroupPicker rv = MRMTransitionGroupPicker.__new__(MRMTransitionGroupPicker)
       rv.inst = shared_ptr[_MRMTransitionGroupPicker](new _MRMTransitionGroupPicker(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMTransitionGroupPicker rv = MRMTransitionGroupPicker.__new__(MRMTransitionGroupPicker)
       rv.inst = shared_ptr[_MRMTransitionGroupPicker](new _MRMTransitionGroupPicker(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMTransitionGroupPicker](new _MRMTransitionGroupPicker())
    
    def _init_1(self, MRMTransitionGroupPicker in_0 ):
        """
        _init_1(self, in_0: MRMTransitionGroupPicker ) -> None
        """
        assert isinstance(in_0, MRMTransitionGroupPicker), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMTransitionGroupPicker](new _MRMTransitionGroupPicker((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMTransitionGroupPicker ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMTransitionGroupPicker)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _pickTransitionGroup_0(self, LightMRMTransitionGroupCP transition_group ):
        """
        _pickTransitionGroup_0(self, transition_group: LightMRMTransitionGroupCP ) -> None
        """
        assert isinstance(transition_group, LightMRMTransitionGroupCP), 'arg transition_group wrong type'
    
        self.inst.get().pickTransitionGroup((deref(transition_group.inst.get())))
    
    def _pickTransitionGroup_1(self, MRMTransitionGroupCP transition_group ):
        """
        _pickTransitionGroup_1(self, transition_group: MRMTransitionGroupCP ) -> None
        """
        assert isinstance(transition_group, MRMTransitionGroupCP), 'arg transition_group wrong type'
    
        self.inst.get().pickTransitionGroup((deref(transition_group.inst.get())))
    
    def pickTransitionGroup(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: pickTransitionGroup(self, transition_group: LightMRMTransitionGroupCP ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: pickTransitionGroup(self, transition_group: MRMTransitionGroupCP ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], LightMRMTransitionGroupCP)):
            return self._pickTransitionGroup_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMTransitionGroupCP)):
            return self._pickTransitionGroup_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def createMRMFeature(self, LightMRMTransitionGroupCP transition_group , list picked_chroms , list smoothed_chroms ,  chr_idx ,  peak_idx ):
        """
        createMRMFeature(self, transition_group: LightMRMTransitionGroupCP , picked_chroms: List[MSChromatogram] , smoothed_chroms: List[MSChromatogram] , chr_idx: int , peak_idx: int ) -> MRMFeature
        """
        assert isinstance(transition_group, LightMRMTransitionGroupCP), 'arg transition_group wrong type'
        assert isinstance(picked_chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in picked_chroms), 'arg picked_chroms wrong type'
        assert isinstance(smoothed_chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in smoothed_chroms), 'arg smoothed_chroms wrong type'
        assert isinstance(chr_idx, int), 'arg chr_idx wrong type'
        assert isinstance(peak_idx, int), 'arg peak_idx wrong type'
    
        cdef libcpp_vector[_MSChromatogram] * v1 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item1
        for item1 in picked_chroms:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[_MSChromatogram] * v2 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item2
        for item2 in smoothed_chroms:
            v2.push_back(deref(item2.inst.get()))
    
    
        cdef _MRMFeature * _r = new _MRMFeature(self.inst.get().createMRMFeature((deref(transition_group.inst.get())), deref(v1), deref(v2), (<const int>chr_idx), (<const int>peak_idx)))
        cdef libcpp_vector[_MSChromatogram].iterator it_smoothed_chroms = v2.begin()
        replace_0 = []
        while it_smoothed_chroms != v2.end():
            item2 = MSChromatogram.__new__(MSChromatogram)
            item2.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_smoothed_chroms)))
            replace_0.append(item2)
            inc(it_smoothed_chroms)
        smoothed_chroms[:] = replace_0
        del v2
        cdef libcpp_vector[_MSChromatogram].iterator it_picked_chroms = v1.begin()
        replace_0 = []
        while it_picked_chroms != v1.end():
            item1 = MSChromatogram.__new__(MSChromatogram)
            item1.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_picked_chroms)))
            replace_0.append(item1)
            inc(it_picked_chroms)
        picked_chroms[:] = replace_0
        del v1
        cdef MRMFeature py_result = MRMFeature.__new__(MRMFeature)
        py_result.inst = shared_ptr[_MRMFeature](_r)
        return py_result
    
    def remove_overlapping_features(self, list picked_chroms , double best_left , double best_right ):
        """
        remove_overlapping_features(self, picked_chroms: List[MSChromatogram] , best_left: float , best_right: float ) -> None
        """
        assert isinstance(picked_chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in picked_chroms), 'arg picked_chroms wrong type'
        assert isinstance(best_left, float), 'arg best_left wrong type'
        assert isinstance(best_right, float), 'arg best_right wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in picked_chroms:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().remove_overlapping_features(deref(v0), (<double>best_left), (<double>best_right))
        cdef libcpp_vector[_MSChromatogram].iterator it_picked_chroms = v0.begin()
        replace_0 = []
        while it_picked_chroms != v0.end():
            item0 = MSChromatogram.__new__(MSChromatogram)
            item0.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_picked_chroms)))
            replace_0.append(item0)
            inc(it_picked_chroms)
        picked_chroms[:] = replace_0
        del v0
    
    def findLargestPeak(self, list picked_chroms ,  chr_idx ,  peak_idx ):
        """
        findLargestPeak(self, picked_chroms: List[MSChromatogram] , chr_idx: int , peak_idx: int ) -> None
        """
        assert isinstance(picked_chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in picked_chroms), 'arg picked_chroms wrong type'
        assert isinstance(chr_idx, int), 'arg chr_idx wrong type'
        assert isinstance(peak_idx, int), 'arg peak_idx wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in picked_chroms:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().findLargestPeak(deref(v0), (<int &>chr_idx), (<int &>peak_idx))
        cdef libcpp_vector[_MSChromatogram].iterator it_picked_chroms = v0.begin()
        replace_0 = []
        while it_picked_chroms != v0.end():
            item0 = MSChromatogram.__new__(MSChromatogram)
            item0.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_picked_chroms)))
            replace_0.append(item0)
            inc(it_picked_chroms)
        picked_chroms[:] = replace_0
        del v0
    
    def findWidestPeakIndices(self, list picked_chroms ,  chrom_idx ,  point_idx ):
        """
        findWidestPeakIndices(self, picked_chroms: List[MSChromatogram] , chrom_idx: int , point_idx: int ) -> None
        """
        assert isinstance(picked_chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in picked_chroms), 'arg picked_chroms wrong type'
        assert isinstance(chrom_idx, int), 'arg chrom_idx wrong type'
        assert isinstance(point_idx, int), 'arg point_idx wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in picked_chroms:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().findWidestPeakIndices(deref(v0), (<int &>chrom_idx), (<int &>point_idx))
        cdef libcpp_vector[_MSChromatogram].iterator it_picked_chroms = v0.begin()
        replace_0 = []
        while it_picked_chroms != v0.end():
            item0 = MSChromatogram.__new__(MSChromatogram)
            item0.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_picked_chroms)))
            replace_0.append(item0)
            inc(it_picked_chroms)
        picked_chroms[:] = replace_0
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

cdef class MSSpectrum:
    """
    Cython implementation of _MSSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSSpectrum.html>`_
      -- Inherits from ['SpectrumSettings', 'RangeManagerMzInt']

    The representation of a 1D spectrum.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    Iterations yields access to underlying peak objects but is slower
    Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
    See help(SpectrumSettings) for information about meta-information
    
    Usage:
    
    .. code-block:: python
    
      ms_level = spectrum.getMSLevel()
      rt = spectrum.getRT()
      mz, intensities = spectrum.get_peaks()
    
    
    Usage:
    
    .. code-block:: python
    
      from pyopenms import *
    
      spectrum = MSSpectrum()
      spectrum.setDriftTime(25) # 25 ms
      spectrum.setRT(205.2) # 205.2 s
      spectrum.setMSLevel(3) # MS3
      p = Precursor()
      p.setIsolationWindowLowerOffset(1.5)
      p.setIsolationWindowUpperOffset(1.5)
      p.setMZ(600) # isolation at 600 +/- 1.5 Th
      p.setActivationEnergy(40) # 40 eV
      p.setCharge(4) # 4+ ion
      spectrum.setPrecursors( [p] )
    
      # Add raw data to spectrum
      spectrum.set_peaks( ([401.5], [900]) )
    
      # Additional data arrays / peak annotations
      fda = FloatDataArray()
      fda.setName("Signal to Noise Array")
      fda.push_back(15)
      sda = StringDataArray()
      sda.setName("Peak annotation")
      sda.push_back("y15++")
      spectrum.setFloatDataArrays( [fda] )
      spectrum.setStringDataArrays( [sda] )
    
      # Add spectrum to MSExperiment
      exp = MSExperiment()
      exp.addSpectrum(spectrum)
    
      # Add second spectrum and store as mzML file
      spectrum2 = MSSpectrum()
      spectrum2.set_peaks( ([1, 2], [1, 2]) )
      exp.addSpectrum(spectrum2)
    
      MzMLFile().store("testfile.mzML", exp)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSSpectrum rv = MSSpectrum.__new__(MSSpectrum)
       rv.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSSpectrum rv = MSSpectrum.__new__(MSSpectrum)
       rv.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSSpectrum](new _MSSpectrum())
    
    def _init_1(self, MSSpectrum in_0 ):
        """
        _init_1(self, in_0: MSSpectrum ) -> None
        """
        assert isinstance(in_0, MSSpectrum), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSSpectrum](new _MSSpectrum((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSSpectrum ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSSpectrum)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getRT(self):
        """
        getRT(self) -> float
        Returns the absolute retention time (in seconds)
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def setRT(self, double in_0 ):
        """
        setRT(self, in_0: float ) -> None
        Sets the absolute retention time (in seconds)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setRT((<double>in_0))
    
    def getDriftTime(self):
        """
        getDriftTime(self) -> float
        Returns the drift time (-1 if not set)
        """
        cdef double _r = self.inst.get().getDriftTime()
        py_result = <double>_r
        return py_result
    
    def setDriftTime(self, double in_0 ):
        """
        setDriftTime(self, in_0: float ) -> None
        Sets the drift time (-1 if not set)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setDriftTime((<double>in_0))
    
    def getDriftTimeUnit(self):
        """
        getDriftTimeUnit(self) -> int
        """
        cdef _DriftTimeUnit _r = self.inst.get().getDriftTimeUnit()
        py_result = <int>_r
        return py_result
    
    def getDriftTimeUnitAsString(self):
        """
        getDriftTimeUnitAsString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getDriftTimeUnitAsString()
        py_result = convOutputString(_r)
        return py_result
    
    def setDriftTimeUnit(self, int dt ):
        """
        setDriftTimeUnit(self, dt: int ) -> None
        """
        assert dt in [0, 1, 2, 3, 4], 'arg dt wrong type'
    
        self.inst.get().setDriftTimeUnit((<_DriftTimeUnit>dt))
    
    def containsIMData(self):
        """
        containsIMData(self) -> bool
        """
        cdef bool _r = self.inst.get().containsIMData()
        py_result = <bool>_r
        return py_result
    
    def getMSLevel(self):
        """
        getMSLevel(self) -> int
        Returns the MS level
        """
        cdef unsigned int _r = self.inst.get().getMSLevel()
        py_result = <unsigned int>_r
        return py_result
    
    def setMSLevel(self,  in_0 ):
        """
        setMSLevel(self, in_0: int ) -> None
        Sets the MS level
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().setMSLevel((<unsigned int>in_0))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get()))
    
    def size(self):
        """
        size(self) -> int
        Returns the number of peaks in the spectrum
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def reserve(self,  n ):
        """
        reserve(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().reserve((<size_t>n))
    
    def resize(self,  n ):
        """
        resize(self, n: int ) -> None
        Resize the peak array
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().resize((<size_t>n))
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> Peak1D
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _Peak1D * _r = new _Peak1D(deref(self.inst.get())[(<size_t>in_0)])
        cdef Peak1D py_result = Peak1D.__new__(Peak1D)
        py_result.inst = shared_ptr[_Peak1D](_r)
        return py_result
    def __setitem__(self, key, Peak1D value):
        """Cython signature: Peak1D & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def updateRanges(self):
        """
        updateRanges(self) -> None
        """
        self.inst.get().updateRanges()
    
    def clear(self, bool clear_meta_data ):
        """
        clear(self, clear_meta_data: bool ) -> None
        Clears all data (and meta data if clear_meta_data is true)
        """
        assert isinstance(clear_meta_data, pybool_t), 'arg clear_meta_data wrong type'
    
        self.inst.get().clear((<bool>clear_meta_data))
    
    def push_back(self, Peak1D in_0 ):
        """
        push_back(self, in_0: Peak1D ) -> None
        Append a peak
        """
        assert isinstance(in_0, Peak1D), 'arg in_0 wrong type'
    
        self.inst.get().push_back((deref(in_0.inst.get())))
    
    def isSorted(self):
        """
        isSorted(self) -> bool
        Returns true if the spectrum is sorte by m/z
        """
        cdef bool _r = self.inst.get().isSorted()
        py_result = <bool>_r
        return py_result
    
    def _findNearest_0(self, double mz ):
        """
        _findNearest_0(self, mz: float ) -> int
        Returns the index of the closest peak in m/z
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef int _r = self.inst.get().findNearest((<double>mz))
        py_result = <int>_r
        return py_result
    
    def _findNearest_1(self, double mz , double tolerance ):
        """
        _findNearest_1(self, mz: float , tolerance: float ) -> int
        Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)
        """
        assert isinstance(mz, float), 'arg mz wrong type'
        assert isinstance(tolerance, float), 'arg tolerance wrong type'
    
    
        cdef int _r = self.inst.get().findNearest((<double>mz), (<double>tolerance))
        py_result = <int>_r
        return py_result
    
    def _findNearest_2(self, double mz , double tolerance_left , double tolerance_right ):
        """
        _findNearest_2(self, mz: float , tolerance_left: float , tolerance_right: float ) -> int
        Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        """
        assert isinstance(mz, float), 'arg mz wrong type'
        assert isinstance(tolerance_left, float), 'arg tolerance_left wrong type'
        assert isinstance(tolerance_right, float), 'arg tolerance_right wrong type'
    
    
    
        cdef int _r = self.inst.get().findNearest((<double>mz), (<double>tolerance_left), (<double>tolerance_right))
        py_result = <int>_r
        return py_result
    
    def findNearest(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: findNearest(self, mz: float ) -> int
          :noindex:
        
        Returns the index of the closest peak in m/z

        
        .. rubric:: Overload:
        .. py:function:: findNearest(self, mz: float , tolerance: float ) -> int
          :noindex:
        
        Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)

        
        .. rubric:: Overload:
        .. py:function:: findNearest(self, mz: float , tolerance_left: float , tolerance_right: float ) -> int
          :noindex:
        
        Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
    
        """
        if (len(args)==1) and (isinstance(args[0], float)):
            return self._findNearest_0(*args)
        elif (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], float)):
            return self._findNearest_1(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._findNearest_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def findHighestInWindow(self, double mz , double tolerance_left , double tolerance_right ):
        """
        findHighestInWindow(self, mz: float , tolerance_left: float , tolerance_right: float ) -> int
        Returns the index of the highest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)
        """
        assert isinstance(mz, float), 'arg mz wrong type'
        assert isinstance(tolerance_left, float), 'arg tolerance_left wrong type'
        assert isinstance(tolerance_right, float), 'arg tolerance_right wrong type'
    
    
    
        cdef int _r = self.inst.get().findHighestInWindow((<double>mz), (<double>tolerance_left), (<double>tolerance_right))
        py_result = <int>_r
        return py_result
    
    def select(self, list indices ):
        """
        select(self, indices: List[int] ) -> MSSpectrum
        Subset the spectrum by indices. Also applies to associated data arrays if present.
        """
        assert isinstance(indices, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in indices), 'arg indices wrong type'
        cdef libcpp_vector[size_t] v0 = indices
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().select(v0))
        indices[:] = v0
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def calculateTIC(self):
        """
        calculateTIC(self) -> float
        Returns the total ion current (=sum) of peak intensities in the spectrum
        """
        cdef double _r = self.inst.get().calculateTIC()
        py_result = <double>_r
        return py_result
    
    def sortByIntensity(self, bool reverse ):
        """
        sortByIntensity(self, reverse: bool ) -> None
        """
        assert isinstance(reverse, pybool_t), 'arg reverse wrong type'
    
        self.inst.get().sortByIntensity((<bool>reverse))
    
    def sortByPosition(self):
        """
        sortByPosition(self) -> None
        """
        self.inst.get().sortByPosition()
    
    def getFloatDataArrays(self):
        """
        getFloatDataArrays(self) -> List[FloatDataArray]
        Returns the additional float data arrays to store e.g. meta data
        """
        _r = self.inst.get().getFloatDataArrays()
        py_result = []
        cdef libcpp_vector[_FloatDataArray].iterator it__r = _r.begin()
        cdef FloatDataArray item_py_result
        while it__r != _r.end():
           item_py_result = FloatDataArray.__new__(FloatDataArray)
           item_py_result.inst = shared_ptr[_FloatDataArray](new _FloatDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getIntegerDataArrays(self):
        """
        getIntegerDataArrays(self) -> List[IntegerDataArray]
        Returns the additional int data arrays to store e.g. meta data
        """
        _r = self.inst.get().getIntegerDataArrays()
        py_result = []
        cdef libcpp_vector[_IntegerDataArray].iterator it__r = _r.begin()
        cdef IntegerDataArray item_py_result
        while it__r != _r.end():
           item_py_result = IntegerDataArray.__new__(IntegerDataArray)
           item_py_result.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getStringDataArrays(self):
        """
        getStringDataArrays(self) -> List[StringDataArray]
        Returns the additional string data arrays to store e.g. meta data
        """
        _r = self.inst.get().getStringDataArrays()
        py_result = []
        cdef libcpp_vector[_StringDataArray].iterator it__r = _r.begin()
        cdef StringDataArray item_py_result
        while it__r != _r.end():
           item_py_result = StringDataArray.__new__(StringDataArray)
           item_py_result.inst = shared_ptr[_StringDataArray](new _StringDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setFloatDataArrays(self, list fda ):
        """
        setFloatDataArrays(self, fda: List[FloatDataArray] ) -> None
        Sets the additional float data arrays to store e.g. meta data
        """
        assert isinstance(fda, list) and all(isinstance(elemt_rec, FloatDataArray) for elemt_rec in fda), 'arg fda wrong type'
        cdef libcpp_vector[_FloatDataArray] * v0 = new libcpp_vector[_FloatDataArray]()
        cdef FloatDataArray item0
        for item0 in fda:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setFloatDataArrays(deref(v0))
        del v0
    
    def setIntegerDataArrays(self, list ida ):
        """
        setIntegerDataArrays(self, ida: List[IntegerDataArray] ) -> None
        Sets the additional int data arrays to store e.g. meta data
        """
        assert isinstance(ida, list) and all(isinstance(elemt_rec, IntegerDataArray) for elemt_rec in ida), 'arg ida wrong type'
        cdef libcpp_vector[_IntegerDataArray] * v0 = new libcpp_vector[_IntegerDataArray]()
        cdef IntegerDataArray item0
        for item0 in ida:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIntegerDataArrays(deref(v0))
        del v0
    
    def setStringDataArrays(self, list sda ):
        """
        setStringDataArrays(self, sda: List[StringDataArray] ) -> None
        Sets the additional string data arrays to store e.g. meta data
        """
        assert isinstance(sda, list) and all(isinstance(elemt_rec, StringDataArray) for elemt_rec in sda), 'arg sda wrong type'
        cdef libcpp_vector[_StringDataArray] * v0 = new libcpp_vector[_StringDataArray]()
        cdef StringDataArray item0
        for item0 in sda:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setStringDataArrays(deref(v0))
        del v0
    
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
    
    def getMinMZ(self):
        """
        getMinMZ(self) -> float
        Returns the minimum m/z
        """
        cdef double _r = self.inst.get().getMinMZ()
        py_result = <double>_r
        return py_result
    
    def getMaxMZ(self):
        """
        getMaxMZ(self) -> float
        Returns the maximum m/z
        """
        cdef double _r = self.inst.get().getMaxMZ()
        py_result = <double>_r
        return py_result
    
    def getMinIntensity(self):
        """
        getMinIntensity(self) -> float
        Returns the minimum intensity
        """
        cdef double _r = self.inst.get().getMinIntensity()
        py_result = <double>_r
        return py_result
    
    def getMaxIntensity(self):
        """
        getMaxIntensity(self) -> float
        Returns the maximum intensity
        """
        cdef double _r = self.inst.get().getMaxIntensity()
        py_result = <double>_r
        return py_result
    
    def clearRanges(self):
        """
        clearRanges(self) -> None
        Resets all range dimensions as empty
        """
        self.inst.get().clearRanges()
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, MSSpectrum):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MSSpectrum other_casted = other
        cdef MSSpectrum self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef Peak1D out
        while it != self.inst.get().end():
            out = Peak1D.__new__(Peak1D)
            out.inst = shared_ptr[_Peak1D](new _Peak1D(deref(it)))
            yield out
            inc(it)
    
    def get_mz_array(MSSpectrum self):
        """Cython signature: numpy_vector get_mz_array()
        
        Will return a numpy array corresponding
        to the mz values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mzs
        mzs = np.empty( (n,), dtype=np.float64)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            inc(it)
            i += 1

        return mzs


    def get_intensity_array(MSSpectrum self):
        """Cython signature: numpy_vector get_intensity_array()
        
        Will return a numpy array corresponding
        to the intensity values in the MSSpectrum.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.empty( (n,), dtype=np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return intensities

    def get_peaks(self):
        """Cython signature: numpy_vector, numpy_vector get_peaks()
        
        Will return a tuple of two numpy arrays (m/z, intensity) corresponding
        to the peaks in the MSSpectrum. Provides fast access to peaks.
        """

        cdef _MSSpectrum * spec_ = self.inst.get()

        cdef unsigned int n = spec_.size()
        cdef np.ndarray[np.float64_t, ndim=1] mzs
        mzs = np.empty( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.empty( (n,), dtype=np.float32)
        cdef _Peak1D p

        cdef libcpp_vector[_Peak1D].iterator it = spec_.begin()
        cdef int i = 0
        while it != spec_.end():
            mzs[i] = deref(it).getMZ()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return mzs, intensities

    def set_peaks(self, peaks):
        """Cython signature: set_peaks((numpy_vector, numpy_vector))
        
        Takes a tuple or list of two arrays (m/z, intensity) and populates the
        MSSpectrum. The arrays can be numpy arrays (faster).
        """

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (mz and intensity vector)"

        mzs, intensities = peaks
        assert len(mzs) == len(intensities), "Input vectors for set_peaks need to have the same length (mz and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float32 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(mzs, intensities)
        elif isinstance(mzs, np.ndarray) and isinstance(intensities, np.ndarray) and \
          mzs.dtype == np.float64 and intensities.dtype == np.float64 and \
          mzs.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(mzs, intensities)
        else:
            self._set_peaks_orig(mzs, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef double intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_mz not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(data_mz)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(data_mz)

        for i in range(N):
            mz = data_mz[i]
            intensity = data_i[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()


    def _set_peaks_orig(self, mzs, intensities):


        cdef _MSSpectrum * spec_ = self.inst.get()

        spec_.resize(0) # empty vector, keep meta data and data arrays
        spec_.reserve(<int>len(mzs)) # allocate space for incoming data
        cdef _Peak1D p = _Peak1D()
        cdef double mz
        cdef float intensity
        cdef int N
        N = len(mzs)

        for i in range(N):
            mz = mzs[i]
            intensity = intensities[i]
            p.setMZ(<double>mz)
            p.setIntensity(<float>intensity)
            spec_.push_back(p)

        spec_.updateRanges()

    def intensityInRange(self, float mzmin, float mzmax):

        cdef double I

        cdef _MSSpectrum * spec_ = self.inst.get()
        cdef int N = spec_.size()

        I = 0.0
        for i in range(N):
                if deref(spec_)[i].getMZ() >= mzmin:
                    break

        cdef _Peak1D * p
        for j in range(i, N):
                p = address(deref(spec_)[j])
                if p.getMZ() > mzmax:
                    break
                I += p.getIntensity()

        return I

    def getIMData(self):

        cdef libcpp_pair[Size, _DriftTimeUnit] r = self.inst.get().getIMData()

        pos = r.first
        unit = <int>r.second

        return (pos, unit) 

cdef class MetaInfo:
    """
    Cython implementation of _MetaInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfo.html>`_

    A Type-Name-Value tuple class
    
    MetaInfo maps an index (an integer corresponding to a string) to
    DataValue objects.  The mapping of strings to the index is performed by
    the MetaInfoRegistry, which can be accessed by the method registry()
    
    There are two versions of nearly all members. One which operates with a
    string name and another one which operates on an index. The index version
    is always faster, as it does not need to look up the index corresponding
    to the string in the MetaInfoRegistry
    
    If you wish to add a MetaInfo member to a class, consider deriving that
    class from MetaInfoInterface, instead of simply adding MetaInfo as
    member. MetaInfoInterface implements a full interface to a MetaInfo
    member and is more memory efficient if no meta info gets added
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaInfo rv = MetaInfo.__new__(MetaInfo)
       rv.inst = shared_ptr[_MetaInfo](new _MetaInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaInfo rv = MetaInfo.__new__(MetaInfo)
       rv.inst = shared_ptr[_MetaInfo](new _MetaInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaInfo](new _MetaInfo())
    
    def _init_1(self, MetaInfo in_0 ):
        """
        _init_1(self, in_0: MetaInfo ) -> None
        """
        assert isinstance(in_0, MetaInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaInfo](new _MetaInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getValue_0(self,  name ):
        """
        _getValue_0(self, name: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to a string
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _DataValue _r = self.inst.get().getValue(deref((convString(name)).get()))
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
    
    def _getValue_1(self,  index ):
        """
        _getValue_1(self, index: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _DataValue _r = self.inst.get().getValue((<unsigned int>index))
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
    
    def _getValue_2(self,  name ,  default_value ):
        """
        _getValue_2(self, name: Union[bytes, str, String] , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to a string
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(default_value, (int, float, list, bytes, str)), 'arg default_value wrong type'
    
    
        cdef _DataValue _r = self.inst.get().getValue(deref((convString(name)).get()), deref(DataValue(default_value).inst.get()))
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
    
    def _getValue_3(self,  index ,  default_value ):
        """
        _getValue_3(self, index: int , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        Returns the value corresponding to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert isinstance(default_value, (int, float, list, bytes, str)), 'arg default_value wrong type'
    
    
        cdef _DataValue _r = self.inst.get().getValue((<unsigned int>index), deref(DataValue(default_value).inst.get()))
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
    
    def getValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getValue(self, name: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        Returns the value corresponding to a string

        
        .. rubric:: Overload:
        .. py:function:: getValue(self, index: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        Returns the value corresponding to an index

        
        .. rubric:: Overload:
        .. py:function:: getValue(self, name: Union[bytes, str, String] , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        Returns the value corresponding to a string

        
        .. rubric:: Overload:
        .. py:function:: getValue(self, index: int , default_value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        Returns the value corresponding to an index
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getValue_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._getValue_1(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._getValue_2(*args)
        elif (len(args)==2) and (isinstance(args[0], int)) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._getValue_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _exists_0(self,  name ):
        """
        _exists_0(self, name: Union[bytes, str, String] ) -> bool
        Returns if this MetaInfo is set
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().exists(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    def _exists_1(self,  index ):
        """
        _exists_1(self, index: int ) -> bool
        Returns if this MetaInfo is set
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef bool _r = self.inst.get().exists((<unsigned int>index))
        py_result = <bool>_r
        return py_result
    
    def exists(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: exists(self, name: Union[bytes, str, String] ) -> bool
          :noindex:
        
        Returns if this MetaInfo is set

        
        .. rubric:: Overload:
        .. py:function:: exists(self, index: int ) -> bool
          :noindex:
        
        Returns if this MetaInfo is set
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._exists_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._exists_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setValue_0(self,  name ,  value ):
        """
        _setValue_0(self, name: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the DataValue corresponding to a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(value, (int, float, list, bytes, str)), 'arg value wrong type'
    
    
        self.inst.get().setValue(deref((convString(name)).get()), deref(DataValue(value).inst.get()))
    
    def _setValue_1(self,  index ,  value ):
        """
        _setValue_1(self, index: int , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        Sets the DataValue corresponding to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert isinstance(value, (int, float, list, bytes, str)), 'arg value wrong type'
    
    
        self.inst.get().setValue((<unsigned int>index), deref(DataValue(value).inst.get()))
    
    def setValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setValue(self, name: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
        
        Sets the DataValue corresponding to a name

        
        .. rubric:: Overload:
        .. py:function:: setValue(self, index: int , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
        
        Sets the DataValue corresponding to an index
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setValue_0(*args)
        elif (len(args)==2) and (isinstance(args[0], int)) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _removeValue_0(self,  name ):
        """
        _removeValue_0(self, name: Union[bytes, str, String] ) -> None
        Removes the DataValue corresponding to `name` if it exists
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().removeValue(deref((convString(name)).get()))
    
    def _removeValue_1(self,  index ):
        """
        _removeValue_1(self, index: int ) -> None
        Removes the DataValue corresponding to `index` if it exists
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        self.inst.get().removeValue((<unsigned int>index))
    
    def removeValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeValue(self, name: Union[bytes, str, String] ) -> None
          :noindex:
        
        Removes the DataValue corresponding to `name` if it exists

        
        .. rubric:: Overload:
        .. py:function:: removeValue(self, index: int ) -> None
          :noindex:
        
        Removes the DataValue corresponding to `index` if it exists
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._removeValue_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._removeValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def getKeysAsIntegers(self, list keys ):
        """
        getKeysAsIntegers(self, keys: List[int] ) -> None
        """
        assert isinstance(keys, list) and all(isinstance(elemt_rec, int) for elemt_rec in keys), 'arg keys wrong type'
        cdef libcpp_vector[unsigned int] v0 = keys
        self.inst.get().getKeys(v0)
        keys[:] = v0
    
    def empty(self):
        """
        empty(self) -> bool
        Returns if the MetaInfo is empty
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        Removes all meta values
        """
        self.inst.get().clear()
    
    def registry(self):
        """
        registry(self) -> MetaInfoRegistry
        """
        cdef _MetaInfoRegistry * _r = new _MetaInfoRegistry(self.inst.get().registry())
        cdef MetaInfoRegistry py_result = MetaInfoRegistry.__new__(MetaInfoRegistry)
        py_result.inst = shared_ptr[_MetaInfoRegistry](_r)
        return py_result 

cdef class NLargest:
    """
    Cython implementation of _NLargest

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NLargest.html>`_
      -- Inherits from ['DefaultParamHandler']

    NLargest removes all but the n largest peaks
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef NLargest rv = NLargest.__new__(NLargest)
       rv.inst = shared_ptr[_NLargest](new _NLargest(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NLargest rv = NLargest.__new__(NLargest)
       rv.inst = shared_ptr[_NLargest](new _NLargest(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NLargest](new _NLargest())
    
    def _init_1(self, NLargest in_0 ):
        """
        _init_1(self, in_0: NLargest ) -> None
        """
        assert isinstance(in_0, NLargest), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NLargest](new _NLargest((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NLargest ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NLargest)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filterSpectrum(self, MSSpectrum spec ):
        """
        filterSpectrum(self, spec: MSSpectrum ) -> None
        Keep only n-largest peaks in spectrum
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterSpectrum((deref(spec.inst.get())))
    
    def filterPeakSpectrum(self, MSSpectrum spec ):
        """
        filterPeakSpectrum(self, spec: MSSpectrum ) -> None
        Keep only n-largest peaks in spectrum
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterPeakSpectrum((deref(spec.inst.get())))
    
    def filterPeakMap(self, MSExperiment exp ):
        """
        filterPeakMap(self, exp: MSExperiment ) -> None
        Keep only n-largest peaks in each spectrum of a peak map
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().filterPeakMap((deref(exp.inst.get())))
    
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

cdef class Nuctide:
    """
    Cython implementation of _Nuctide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Nuctide.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property rts:
        def __set__(self, list rts):
            cdef libcpp_vector[_RetentionTime] * v0 = new libcpp_vector[_RetentionTime]()
            cdef RetentionTime item0
            for item0 in rts:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().rts = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().rts
            py_result = []
            cdef libcpp_vector[_RetentionTime].iterator it__r = _r.begin()
            cdef RetentionTime item_py_result
            while it__r != _r.end():
               item_py_result = RetentionTime.__new__(RetentionTime)
               item_py_result.inst = shared_ptr[_RetentionTime](new _RetentionTime(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property oligo_refs:
        def __set__(self, list oligo_refs):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in oligo_refs:
               v0.push_back(_String(<char *>item0))
            self.inst.get().oligo_refs = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().oligo_refs
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property evidence:
        def __set__(self, CVTermList evidence):
        
            self.inst.get().evidence = (deref(evidence.inst.get()))
        
    
        def __get__(self):
            cdef _CVTermList * _r = new _CVTermList(self.inst.get().evidence)
            cdef CVTermList py_result = CVTermList.__new__(CVTermList)
            py_result.inst = shared_ptr[_CVTermList](_r)
            return py_result
    
    property sequence:
        def __set__(self,  sequence):
        
            self.inst.get().sequence = deref((convString(sequence)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().sequence
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef Nuctide rv = Nuctide.__new__(Nuctide)
       rv.inst = shared_ptr[_Nuctide](new _Nuctide(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Nuctide rv = Nuctide.__new__(Nuctide)
       rv.inst = shared_ptr[_Nuctide](new _Nuctide(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Nuctide](new _Nuctide())
    
    def _init_1(self, Nuctide in_0 ):
        """
        _init_1(self, in_0: Nuctide ) -> None
        """
        assert isinstance(in_0, Nuctide), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Nuctide](new _Nuctide((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Nuctide ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Nuctide)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setNuctideGroupLabel(self,  label ):
        """
        setNuctideGroupLabel(self, label: Union[bytes, str, String] ) -> None
        Sets the nuctide group label
        """
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
        self.inst.get().setNuctideGroupLabel(deref((convString(label)).get()))
    
    def getNuctideGroupLabel(self):
        """
        getNuctideGroupLabel(self) -> Union[bytes, str, String]
        Get the nuctide group label
        """
        cdef _String _r = self.inst.get().getNuctideGroupLabel()
        py_result = convOutputString(_r)
        return py_result
    
    def setChargeState(self,  charge ):
        """
        setChargeState(self, charge: int ) -> None
        Sets the nuctide or compound charge states
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setChargeState((<int>charge))
    
    def getChargeState(self):
        """
        getChargeState(self) -> int
        Returns the nuctide charge state
        """
        cdef int _r = self.inst.get().getChargeState()
        py_result = <int>_r
        return py_result
    
    def hasCharge(self):
        """
        hasCharge(self) -> bool
        Whether product has set charge state
        """
        cdef bool _r = self.inst.get().hasCharge()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTime(self):
        """
        getRetentionTime(self) -> float
        Gets nuctide retention time
        """
        cdef double _r = self.inst.get().getRetentionTime()
        py_result = <double>_r
        return py_result
    
    def hasRetentionTime(self):
        """
        hasRetentionTime(self) -> bool
        Gets nuctide retention time
        """
        cdef bool _r = self.inst.get().hasRetentionTime()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTimeType(self):
        """
        getRetentionTimeType(self) -> int
        Get nuctide retentiontime type
        """
        cdef _RTType _r = self.inst.get().getRetentionTimeType()
        py_result = <int>_r
        return py_result
    
    def getRetentionTimeUnit(self):
        """
        getRetentionTimeUnit(self) -> int
        Get nuctide retentiontime unit (minute/seconds)
        """
        cdef _RTUnit _r = self.inst.get().getRetentionTimeUnit()
        py_result = <int>_r
        return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Nuctide):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Nuctide other_casted = other
        cdef Nuctide self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Peptide:
    """
    Cython implementation of _Peptide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Peptide.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property rts:
        def __set__(self, list rts):
            cdef libcpp_vector[_RetentionTime] * v0 = new libcpp_vector[_RetentionTime]()
            cdef RetentionTime item0
            for item0 in rts:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().rts = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().rts
            py_result = []
            cdef libcpp_vector[_RetentionTime].iterator it__r = _r.begin()
            cdef RetentionTime item_py_result
            while it__r != _r.end():
               item_py_result = RetentionTime.__new__(RetentionTime)
               item_py_result.inst = shared_ptr[_RetentionTime](new _RetentionTime(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property protein_refs:
        def __set__(self, list protein_refs):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in protein_refs:
               v0.push_back(_String(<char *>item0))
            self.inst.get().protein_refs = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().protein_refs
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property evidence:
        def __set__(self, CVTermList evidence):
        
            self.inst.get().evidence = (deref(evidence.inst.get()))
        
    
        def __get__(self):
            cdef _CVTermList * _r = new _CVTermList(self.inst.get().evidence)
            cdef CVTermList py_result = CVTermList.__new__(CVTermList)
            py_result.inst = shared_ptr[_CVTermList](_r)
            return py_result
    
    property sequence:
        def __set__(self,  sequence):
        
            self.inst.get().sequence = deref((convString(sequence)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().sequence
            py_result = convOutputString(_r)
            return py_result
    
    property mods:
        def __set__(self, list mods):
            cdef libcpp_vector[_TargetedExperiment_Modification] * v0 = new libcpp_vector[_TargetedExperiment_Modification]()
            cdef TargetedExperiment_Modification item0
            for item0 in mods:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().mods = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().mods
            py_result = []
            cdef libcpp_vector[_TargetedExperiment_Modification].iterator it__r = _r.begin()
            cdef TargetedExperiment_Modification item_py_result
            while it__r != _r.end():
               item_py_result = TargetedExperiment_Modification.__new__(TargetedExperiment_Modification)
               item_py_result.inst = shared_ptr[_TargetedExperiment_Modification](new _TargetedExperiment_Modification(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef Peptide rv = Peptide.__new__(Peptide)
       rv.inst = shared_ptr[_Peptide](new _Peptide(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Peptide rv = Peptide.__new__(Peptide)
       rv.inst = shared_ptr[_Peptide](new _Peptide(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Peptide](new _Peptide())
    
    def _init_1(self, Peptide in_0 ):
        """
        _init_1(self, in_0: Peptide ) -> None
        """
        assert isinstance(in_0, Peptide), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Peptide](new _Peptide((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Peptide ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Peptide)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setPeptideGroupLabel(self,  label ):
        """
        setPeptideGroupLabel(self, label: Union[bytes, str, String] ) -> None
        Sets the peptide group label
        """
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
        self.inst.get().setPeptideGroupLabel(deref((convString(label)).get()))
    
    def getPeptideGroupLabel(self):
        """
        getPeptideGroupLabel(self) -> Union[bytes, str, String]
        Get the peptide group label
        """
        cdef _String _r = self.inst.get().getPeptideGroupLabel()
        py_result = convOutputString(_r)
        return py_result
    
    def setChargeState(self,  charge ):
        """
        setChargeState(self, charge: int ) -> None
        Sets the peptide or compound charge states
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setChargeState((<int>charge))
    
    def getChargeState(self):
        """
        getChargeState(self) -> int
        Returns the peptide or compound charge state
        """
        cdef int _r = self.inst.get().getChargeState()
        py_result = <int>_r
        return py_result
    
    def hasCharge(self):
        """
        hasCharge(self) -> bool
        Whether product has set charge state
        """
        cdef bool _r = self.inst.get().hasCharge()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTime(self):
        """
        getRetentionTime(self) -> float
        Gets compound or peptide retention time
        """
        cdef double _r = self.inst.get().getRetentionTime()
        py_result = <double>_r
        return py_result
    
    def hasRetentionTime(self):
        """
        hasRetentionTime(self) -> bool
        Gets compound or peptide retention time
        """
        cdef bool _r = self.inst.get().hasRetentionTime()
        py_result = <bool>_r
        return py_result
    
    def getRetentionTimeType(self):
        """
        getRetentionTimeType(self) -> int
        Get compound or peptide retentiontime type
        """
        cdef _RTType _r = self.inst.get().getRetentionTimeType()
        py_result = <int>_r
        return py_result
    
    def getRetentionTimeUnit(self):
        """
        getRetentionTimeUnit(self) -> int
        Get compound or peptide retentiontime unit (minute/seconds)
        """
        cdef _RTUnit _r = self.inst.get().getRetentionTimeUnit()
        py_result = <int>_r
        return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Peptide):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Peptide other_casted = other
        cdef Peptide self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Prediction:
    """
    Cython implementation of _Prediction

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Prediction.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property software_ref:
        def __set__(self,  software_ref):
        
            self.inst.get().software_ref = deref((convString(software_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().software_ref
            py_result = convOutputString(_r)
            return py_result
    
    property contact_ref:
        def __set__(self,  contact_ref):
        
            self.inst.get().contact_ref = deref((convString(contact_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().contact_ref
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef Prediction rv = Prediction.__new__(Prediction)
       rv.inst = shared_ptr[_Prediction](new _Prediction(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Prediction rv = Prediction.__new__(Prediction)
       rv.inst = shared_ptr[_Prediction](new _Prediction(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Prediction](new _Prediction())
    
    def _init_1(self, Prediction in_0 ):
        """
        _init_1(self, in_0: Prediction ) -> None
        """
        assert isinstance(in_0, Prediction), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Prediction](new _Prediction((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Prediction ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Prediction)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Prediction):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Prediction other_casted = other
        cdef Prediction self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Protein:
    """
    Cython implementation of _Protein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Protein.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    property sequence:
        def __set__(self,  sequence):
        
            self.inst.get().sequence = deref((convString(sequence)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().sequence
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef Protein rv = Protein.__new__(Protein)
       rv.inst = shared_ptr[_Protein](new _Protein(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Protein rv = Protein.__new__(Protein)
       rv.inst = shared_ptr[_Protein](new _Protein(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Protein](new _Protein())
    
    def _init_1(self, Protein in_0 ):
        """
        _init_1(self, in_0: Protein ) -> None
        """
        assert isinstance(in_0, Protein), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Protein](new _Protein((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Protein ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Protein)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Protein):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Protein other_casted = other
        cdef Protein self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ProteinHit:
    """
    Cython implementation of _ProteinHit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinHit.html>`_
      -- Inherits from ['MetaInfoInterface']

    Representation of a protein hit
    
    It contains the fields score, score_type, rank, accession,
    sequence and coverage
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ProteinHit rv = ProteinHit.__new__(ProteinHit)
       rv.inst = shared_ptr[_ProteinHit](new _ProteinHit(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteinHit rv = ProteinHit.__new__(ProteinHit)
       rv.inst = shared_ptr[_ProteinHit](new _ProteinHit(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteinHit](new _ProteinHit())
    
    def _init_1(self, double score ,  rank ,  accession ,  sequence ):
        """
        _init_1(self, score: float , rank: int , accession: Union[bytes, str, String] , sequence: Union[bytes, str, String] ) -> None
        """
        assert isinstance(score, float), 'arg score wrong type'
        assert isinstance(rank, int), 'arg rank wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
    
    
    
    
        self.inst = shared_ptr[_ProteinHit](new _ProteinHit((<double>score), (<unsigned int>rank), deref((convString(accession)).get()), deref((convString(sequence)).get())))
    
    def _init_2(self, ProteinHit in_0 ):
        """
        _init_2(self, in_0: ProteinHit ) -> None
        """
        assert isinstance(in_0, ProteinHit), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteinHit](new _ProteinHit((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, score: float , rank: int , accession: Union[bytes, str, String] , sequence: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteinHit ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], int)) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteinHit)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getScore(self):
        """
        getScore(self) -> float
        Returns the score of the protein hit
        """
        cdef float _r = self.inst.get().getScore()
        py_result = <float>_r
        return py_result
    
    def getRank(self):
        """
        getRank(self) -> int
        Returns the rank of the protein hit
        """
        cdef unsigned int _r = self.inst.get().getRank()
        py_result = <unsigned int>_r
        return py_result
    
    def getSequence(self):
        """
        getSequence(self) -> Union[bytes, str, String]
        Returns the protein sequence
        """
        cdef _String _r = self.inst.get().getSequence()
        py_result = convOutputString(_r)
        return py_result
    
    def getAccession(self):
        """
        getAccession(self) -> Union[bytes, str, String]
        Returns the accession of the protein
        """
        cdef _String _r = self.inst.get().getAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def getDescription(self):
        """
        getDescription(self) -> Union[bytes, str, String]
        Returns the description of the protein
        """
        cdef _String _r = self.inst.get().getDescription()
        py_result = convOutputString(_r)
        return py_result
    
    def getCoverage(self):
        """
        getCoverage(self) -> float
        Returns the coverage (in percent) of the protein hit based upon matched peptides
        """
        cdef double _r = self.inst.get().getCoverage()
        py_result = <double>_r
        return py_result
    
    def setScore(self, float in_0 ):
        """
        setScore(self, in_0: float ) -> None
        Sets the score of the protein hit
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setScore((<float>in_0))
    
    def setRank(self,  in_0 ):
        """
        setRank(self, in_0: int ) -> None
        Sets the rank
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().setRank((<unsigned int>in_0))
    
    def setSequence(self,  in_0 ):
        """
        setSequence(self, in_0: Union[bytes, str, String] ) -> None
        Sets the protein sequence
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setSequence(deref((convString(in_0)).get()))
    
    def setAccession(self,  in_0 ):
        """
        setAccession(self, in_0: Union[bytes, str, String] ) -> None
        Sets the accession of the protein
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setAccession(deref((convString(in_0)).get()))
    
    def setDescription(self,  description ):
        """
        setDescription(self, description: Union[bytes, str, String] ) -> None
        Sets the description of the protein
        """
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
    
        self.inst.get().setDescription(deref((convString(description)).get()))
    
    def setCoverage(self, double in_0 ):
        """
        setCoverage(self, in_0: float ) -> None
        Sets the coverage (in percent) of the protein hit based upon matched peptides
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setCoverage((<double>in_0))
    
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
        if not isinstance(other, ProteinHit):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ProteinHit other_casted = other
        cdef ProteinHit self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Publication:
    """
    Cython implementation of _Publication

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1Publication.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef Publication rv = Publication.__new__(Publication)
       rv.inst = shared_ptr[_Publication](new _Publication(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Publication rv = Publication.__new__(Publication)
       rv.inst = shared_ptr[_Publication](new _Publication(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Publication](new _Publication())
    
    def _init_1(self, Publication in_0 ):
        """
        _init_1(self, in_0: Publication ) -> None
        """
        assert isinstance(in_0, Publication), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Publication](new _Publication((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Publication ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Publication)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, Publication):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Publication other_casted = other
        cdef Publication self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ReactionMonitoringTransition:
    """
    Cython implementation of _ReactionMonitoringTransition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ReactionMonitoringTransition.html>`_
      -- Inherits from ['CVTermList']

    This class stores a SRM/MRM transition
    
    This class is capable of representing a <Transition> tag in a TraML
    document completely and contains all associated information
    
    The default values for precursor m/z is 0.0 which indicates that it is
    uninitialized
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ReactionMonitoringTransition rv = ReactionMonitoringTransition.__new__(ReactionMonitoringTransition)
       rv.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ReactionMonitoringTransition rv = ReactionMonitoringTransition.__new__(ReactionMonitoringTransition)
       rv.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition())
    
    def _init_1(self, ReactionMonitoringTransition in_0 ):
        """
        _init_1(self, in_0: ReactionMonitoringTransition ) -> None
        """
        assert isinstance(in_0, ReactionMonitoringTransition), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ReactionMonitoringTransition ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ReactionMonitoringTransition)):
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
    
    def getNativeID(self):
        """
        getNativeID(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getNativeID()
        py_result = convOutputString(_r)
        return py_result
    
    def getTransRef(self):
        """
        getTransRef(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTransRef()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def setNativeID(self,  name ):
        """
        setNativeID(self, name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setNativeID(deref((convString(name)).get()))
    
    def setTransRef(self,  trans_ref , int type ):
        """
        setTransRef(self, trans_ref: Union[bytes, str, String] , type: int ) -> None
        """
        assert (isinstance(trans_ref, str) or isinstance(trans_ref, bytes) or isinstance(trans_ref, String)), 'arg trans_ref wrong type'
        assert type in [0, 1, 2, 3], 'arg type wrong type'
    
    
        self.inst.get().setTransRef(deref((convString(trans_ref)).get()), (<_TransType>type))
    
    def getTransType(self):
        """
        getTransType(self) -> int
        """
        cdef _TransType _r = self.inst.get().getTransType()
        py_result = <int>_r
        return py_result
    
    def getProductMZ(self):
        """
        getProductMZ(self) -> float
        """
        cdef double _r = self.inst.get().getProductMZ()
        py_result = <double>_r
        return py_result
    
    def setProductMZ(self, double in_0 ):
        """
        setProductMZ(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setProductMZ((<double>in_0))
    
    def getPrecursorMZ(self):
        """
        getPrecursorMZ(self) -> float
        Returns the precursor mz (Q1 value)
        """
        cdef double _r = self.inst.get().getPrecursorMZ()
        py_result = <double>_r
        return py_result
    
    def setPrecursorMZ(self, double in_0 ):
        """
        setPrecursorMZ(self, in_0: float ) -> None
        Sets the precursor mz (Q1 value)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setPrecursorMZ((<double>in_0))
    
    def getDecoyTransitionType(self):
        """
        getDecoyTransitionType(self) -> int
        Returns the type of transition (target or decoy)
        """
        cdef _DecoyTransitionType _r = self.inst.get().getDecoyTransitionType()
        py_result = <int>_r
        return py_result
    
    def hasPrecursorCVTerms(self):
        """
        hasPrecursorCVTerms(self) -> bool
        Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)
        """
        cdef bool _r = self.inst.get().hasPrecursorCVTerms()
        py_result = <bool>_r
        return py_result
    
    def setPrecursorCVTermList(self, CVTermList list_ ):
        """
        setPrecursorCVTermList(self, list_: CVTermList ) -> None
        Sets a list of precursor CV Terms
        """
        assert isinstance(list_, CVTermList), 'arg list_ wrong type'
    
        self.inst.get().setPrecursorCVTermList((deref(list_.inst.get())))
    
    def addPrecursorCVTerm(self, CVTerm cv_term ):
        """
        addPrecursorCVTerm(self, cv_term: CVTerm ) -> None
        Adds precursor CV Term
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().addPrecursorCVTerm((deref(cv_term.inst.get())))
    
    def getPrecursorCVTermList(self):
        """
        getPrecursorCVTermList(self) -> CVTermList
        Obtains the list of CV Terms for the precursor
        """
        cdef _CVTermList * _r = new _CVTermList(self.inst.get().getPrecursorCVTermList())
        cdef CVTermList py_result = CVTermList.__new__(CVTermList)
        py_result.inst = shared_ptr[_CVTermList](_r)
        return py_result
    
    def addProductCVTerm(self, CVTerm cv_term ):
        """
        addProductCVTerm(self, cv_term: CVTerm ) -> None
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().addProductCVTerm((deref(cv_term.inst.get())))
    
    def getIntermediateProducts(self):
        """
        getIntermediateProducts(self) -> List[TraMLProduct]
        """
        _r = self.inst.get().getIntermediateProducts()
        py_result = []
        cdef libcpp_vector[_TraMLProduct].iterator it__r = _r.begin()
        cdef TraMLProduct item_py_result
        while it__r != _r.end():
           item_py_result = TraMLProduct.__new__(TraMLProduct)
           item_py_result.inst = shared_ptr[_TraMLProduct](new _TraMLProduct(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addIntermediateProduct(self, TraMLProduct product ):
        """
        addIntermediateProduct(self, product: TraMLProduct ) -> None
        """
        assert isinstance(product, TraMLProduct), 'arg product wrong type'
    
        self.inst.get().addIntermediateProduct((deref(product.inst.get())))
    
    def setIntermediateProducts(self, list products ):
        """
        setIntermediateProducts(self, products: List[TraMLProduct] ) -> None
        """
        assert isinstance(products, list) and all(isinstance(elemt_rec, TraMLProduct) for elemt_rec in products), 'arg products wrong type'
        cdef libcpp_vector[_TraMLProduct] * v0 = new libcpp_vector[_TraMLProduct]()
        cdef TraMLProduct item0
        for item0 in products:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIntermediateProducts(deref(v0))
        cdef libcpp_vector[_TraMLProduct].iterator it_products = v0.begin()
        replace_0 = []
        while it_products != v0.end():
            item0 = TraMLProduct.__new__(TraMLProduct)
            item0.inst = shared_ptr[_TraMLProduct](new _TraMLProduct(deref(it_products)))
            replace_0.append(item0)
            inc(it_products)
        products[:] = replace_0
        del v0
    
    def setProduct(self, TraMLProduct product ):
        """
        setProduct(self, product: TraMLProduct ) -> None
        """
        assert isinstance(product, TraMLProduct), 'arg product wrong type'
    
        self.inst.get().setProduct((deref(product.inst.get())))
    
    def getProduct(self):
        """
        getProduct(self) -> TraMLProduct
        """
        cdef _TraMLProduct * _r = new _TraMLProduct(self.inst.get().getProduct())
        cdef TraMLProduct py_result = TraMLProduct.__new__(TraMLProduct)
        py_result.inst = shared_ptr[_TraMLProduct](_r)
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
    
    def setPrediction(self, Prediction prediction ):
        """
        setPrediction(self, prediction: Prediction ) -> None
        Sets prediction
        """
        assert isinstance(prediction, Prediction), 'arg prediction wrong type'
    
        self.inst.get().setPrediction((deref(prediction.inst.get())))
    
    def addPredictionTerm(self, CVTerm prediction ):
        """
        addPredictionTerm(self, prediction: CVTerm ) -> None
        Adds prediction term
        """
        assert isinstance(prediction, CVTerm), 'arg prediction wrong type'
    
        self.inst.get().addPredictionTerm((deref(prediction.inst.get())))
    
    def hasPrediction(self):
        """
        hasPrediction(self) -> bool
        Returns true if a Prediction object exists (means it is safe to call getPrediction)
        """
        cdef bool _r = self.inst.get().hasPrediction()
        py_result = <bool>_r
        return py_result
    
    def getPrediction(self):
        """
        getPrediction(self) -> Prediction
        Obtains the Prediction object
        """
        cdef _Prediction * _r = new _Prediction(self.inst.get().getPrediction())
        cdef Prediction py_result = Prediction.__new__(Prediction)
        py_result.inst = shared_ptr[_Prediction](_r)
        return py_result
    
    def setDecoyTransitionType(self, int d ):
        """
        setDecoyTransitionType(self, d: int ) -> None
        Sets the type of transition (target or decoy)
        """
        assert d in [0, 1, 2], 'arg d wrong type'
    
        self.inst.get().setDecoyTransitionType((<_DecoyTransitionType>d))
    
    def getLibraryIntensity(self):
        """
        getLibraryIntensity(self) -> float
        Returns the library intensity (ion count or normalized ion count from a spectral library)
        """
        cdef double _r = self.inst.get().getLibraryIntensity()
        py_result = <double>_r
        return py_result
    
    def setLibraryIntensity(self, double intensity ):
        """
        setLibraryIntensity(self, intensity: float ) -> None
        Sets the library intensity (ion count or normalized ion count from a spectral library)
        """
        assert isinstance(intensity, float), 'arg intensity wrong type'
    
        self.inst.get().setLibraryIntensity((<double>intensity))
    
    def getProductChargeState(self):
        """
        getProductChargeState(self) -> int
        Returns the charge state of the product
        """
        cdef int _r = self.inst.get().getProductChargeState()
        py_result = <int>_r
        return py_result
    
    def isProductChargeStateSet(self):
        """
        isProductChargeStateSet(self) -> bool
        Returns true if charge state of product is already set
        """
        cdef bool _r = self.inst.get().isProductChargeStateSet()
        py_result = <bool>_r
        return py_result
    
    def isDetectingTransition(self):
        """
        isDetectingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isDetectingTransition()
        py_result = <bool>_r
        return py_result
    
    def setDetectingTransition(self, bool val ):
        """
        setDetectingTransition(self, val: bool ) -> None
        """
        assert isinstance(val, pybool_t), 'arg val wrong type'
    
        self.inst.get().setDetectingTransition((<bool>val))
    
    def isIdentifyingTransition(self):
        """
        isIdentifyingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isIdentifyingTransition()
        py_result = <bool>_r
        return py_result
    
    def setIdentifyingTransition(self, bool val ):
        """
        setIdentifyingTransition(self, val: bool ) -> None
        """
        assert isinstance(val, pybool_t), 'arg val wrong type'
    
        self.inst.get().setIdentifyingTransition((<bool>val))
    
    def isQuantifyingTransition(self):
        """
        isQuantifyingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isQuantifyingTransition()
        py_result = <bool>_r
        return py_result
    
    def setQuantifyingTransition(self, bool val ):
        """
        setQuantifyingTransition(self, val: bool ) -> None
        """
        assert isinstance(val, pybool_t), 'arg val wrong type'
    
        self.inst.get().setQuantifyingTransition((<bool>val))
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, ReactionMonitoringTransition):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ReactionMonitoringTransition other_casted = other
        cdef ReactionMonitoringTransition self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class RetentionTime:
    """
    Cython implementation of _RetentionTime

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1RetentionTime.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property software_ref:
        def __set__(self,  software_ref):
        
            self.inst.get().software_ref = deref((convString(software_ref)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().software_ref
            py_result = convOutputString(_r)
            return py_result
    
    property retention_time_unit:
        def __set__(self, int retention_time_unit):
        
            self.inst.get().retention_time_unit = (<_RTUnit>retention_time_unit)
        
    
        def __get__(self):
            cdef _RTUnit _r = self.inst.get().retention_time_unit
            py_result = <int>_r
            return py_result
    
    property retention_time_type:
        def __set__(self, int retention_time_type):
        
            self.inst.get().retention_time_type = (<_RTType>retention_time_type)
        
    
        def __get__(self):
            cdef _RTType _r = self.inst.get().retention_time_type
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef RetentionTime rv = RetentionTime.__new__(RetentionTime)
       rv.inst = shared_ptr[_RetentionTime](new _RetentionTime(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RetentionTime rv = RetentionTime.__new__(RetentionTime)
       rv.inst = shared_ptr[_RetentionTime](new _RetentionTime(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RetentionTime](new _RetentionTime())
    
    def _init_1(self, RetentionTime in_0 ):
        """
        _init_1(self, in_0: RetentionTime ) -> None
        """
        assert isinstance(in_0, RetentionTime), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RetentionTime](new _RetentionTime((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RetentionTime ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RetentionTime)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isRTset(self):
        """
        isRTset(self) -> bool
        """
        cdef bool _r = self.inst.get().isRTset()
        py_result = <bool>_r
        return py_result
    
    def setRT(self, double rt ):
        """
        setRT(self, rt: float ) -> None
        """
        assert isinstance(rt, float), 'arg rt wrong type'
    
        self.inst.get().setRT((<double>rt))
    
    def getRT(self):
        """
        getRT(self) -> float
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
        if not isinstance(other, RetentionTime):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef RetentionTime other_casted = other
        cdef RetentionTime self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    RTType = __RTType
    RTUnit = __RTUnit 

cdef class SpectrumAlignmentScore:
    """
    Cython implementation of _SpectrumAlignmentScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAlignmentScore.html>`_
      -- Inherits from ['DefaultParamHandler']

    Similarity score via spectra alignment
    
    This class implements a simple scoring based on the alignment of spectra. This alignment
    is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
    of the peaks, minimizing the distances between the aligned peaks and maximizing the number
    of peak pairs
    
    The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the
    product of the intensities of the aligned peaks, with the given exponent (default is 2)
    sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAlignmentScore rv = SpectrumAlignmentScore.__new__(SpectrumAlignmentScore)
       rv.inst = shared_ptr[_SpectrumAlignmentScore](new _SpectrumAlignmentScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAlignmentScore rv = SpectrumAlignmentScore.__new__(SpectrumAlignmentScore)
       rv.inst = shared_ptr[_SpectrumAlignmentScore](new _SpectrumAlignmentScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectrumAlignmentScore](new _SpectrumAlignmentScore())
    
    def _init_1(self, SpectrumAlignmentScore in_0 ):
        """
        _init_1(self, in_0: SpectrumAlignmentScore ) -> None
        """
        assert isinstance(in_0, SpectrumAlignmentScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAlignmentScore](new _SpectrumAlignmentScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAlignmentScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAlignmentScore)):
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
    
    def __call__(self, *args):
        assert len(args) in (1, 2), len(args)
        if len(args) == 1:
            return self._score_1(args[0])
        return self._score_2(args[0], args[1])

    def _score_1(self, spec1):
        assert isinstance(spec1, (MSSpectrum,))
        cdef _MSSpectrum * _new_spec1 = new _MSSpectrum()
        cdef _Peak1D _peak

        if True:
            for _peak in deref(<_MSSpectrum *>(<MSSpectrum>spec1).inst.get()):
                _new_spec1.push_back(<_Peak1D>_peak)

        cdef _SpectrumAlignmentScore * scorer = self.inst.get()
        cdef double score = <double>deref(scorer)(deref(_new_spec1))
        del _new_spec1
        return score

    def _score_2(self, spec1, spec2):
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

        cdef _SpectrumAlignmentScore * scorer = self.inst.get()
        cdef double score = <double>deref(scorer)(deref(_new_spec1), deref(_new_spec2))
        del _new_spec1
        del _new_spec2
        return score 

cdef class TargetedExperiment_Instrument:
    """
    Cython implementation of _TargetedExperiment_Instrument

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Instrument.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = deref((convString(id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().id
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef TargetedExperiment_Instrument rv = TargetedExperiment_Instrument.__new__(TargetedExperiment_Instrument)
       rv.inst = shared_ptr[_TargetedExperiment_Instrument](new _TargetedExperiment_Instrument(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TargetedExperiment_Instrument rv = TargetedExperiment_Instrument.__new__(TargetedExperiment_Instrument)
       rv.inst = shared_ptr[_TargetedExperiment_Instrument](new _TargetedExperiment_Instrument(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TargetedExperiment_Instrument](new _TargetedExperiment_Instrument())
    
    def _init_1(self, TargetedExperiment_Instrument in_0 ):
        """
        _init_1(self, in_0: TargetedExperiment_Instrument ) -> None
        """
        assert isinstance(in_0, TargetedExperiment_Instrument), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TargetedExperiment_Instrument](new _TargetedExperiment_Instrument((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TargetedExperiment_Instrument ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TargetedExperiment_Instrument)):
             self._init_1(*args)
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
    
    def consumeCVTerms(self, dict cv_term_map ):
        """
        consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
    
    def getKeys(self, list keys ):
        """
        getKeys(self, keys: List[bytes] ) -> None
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
    
    def getKeysAsIntegers(self, list keys ):
        """
        getKeysAsIntegers(self, keys: List[int] ) -> None
        """
        assert isinstance(keys, list) and all(isinstance(elemt_rec, int) for elemt_rec in keys), 'arg keys wrong type'
        cdef libcpp_vector[unsigned int] v0 = keys
        self.inst.get().getKeys(v0)
        keys[:] = v0
    
    def _getMetaValue_0(self,  in_0 ):
        """
        _getMetaValue_0(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        cdef _DataValue _r = self.inst.get().getMetaValue((<unsigned int>in_0))
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
    
    def _getMetaValue_1(self,  in_0 ):
        """
        _getMetaValue_1(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
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
    
    def getMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMetaValue(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int)):
            return self._getMetaValue_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setMetaValue_0(self,  in_0 ,  in_1 ):
        """
        _setMetaValue_0(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue((<unsigned int>in_0), deref(DataValue(in_1).inst.get()))
    
    def _setMetaValue_1(self,  in_0 ,  in_1 ):
        """
        _setMetaValue_1(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue(deref((convString(in_0)).get()), deref(DataValue(in_1).inst.get()))
    
    def setMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMetaValue(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], int)) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setMetaValue_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _metaValueExists_0(self,  in_0 ):
        """
        _metaValueExists_0(self, in_0: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists(deref((convString(in_0)).get()))
        py_result = <bool>_r
        return py_result
    
    def _metaValueExists_1(self,  in_0 ):
        """
        _metaValueExists_1(self, in_0: int ) -> bool
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists((<unsigned int>in_0))
        py_result = <bool>_r
        return py_result
    
    def metaValueExists(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: metaValueExists(self, in_0: int ) -> bool
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._metaValueExists_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._metaValueExists_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _removeMetaValue_0(self,  in_0 ):
        """
        _removeMetaValue_0(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue(deref((convString(in_0)).get()))
    
    def _removeMetaValue_1(self,  in_0 ):
        """
        _removeMetaValue_1(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue((<unsigned int>in_0))
    
    def removeMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: removeMetaValue(self, in_0: int ) -> None
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._removeMetaValue_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._removeMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, TargetedExperiment_Instrument):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TargetedExperiment_Instrument other_casted = other
        cdef TargetedExperiment_Instrument self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class TargetedExperiment_Interpretation:
    """
    Cython implementation of _TargetedExperiment_Interpretation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Interpretation.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property ordinal:
        def __set__(self, bytes ordinal):
        
            self.inst.get().ordinal = (<char>((ordinal)[0]))
        
    
        def __get__(self):
            cdef char  _r = self.inst.get().ordinal
            py_result = chr(<char>(_r))
            return py_result
    
    property rank:
        def __set__(self, bytes rank):
        
            self.inst.get().rank = (<char>((rank)[0]))
        
    
        def __get__(self):
            cdef char  _r = self.inst.get().rank
            py_result = chr(<char>(_r))
            return py_result
    
    property iontype:
        def __set__(self, int iontype):
        
            self.inst.get().iontype = (<_ResidueType>iontype)
        
    
        def __get__(self):
            cdef _ResidueType _r = self.inst.get().iontype
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef TargetedExperiment_Interpretation rv = TargetedExperiment_Interpretation.__new__(TargetedExperiment_Interpretation)
       rv.inst = shared_ptr[_TargetedExperiment_Interpretation](new _TargetedExperiment_Interpretation(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TargetedExperiment_Interpretation rv = TargetedExperiment_Interpretation.__new__(TargetedExperiment_Interpretation)
       rv.inst = shared_ptr[_TargetedExperiment_Interpretation](new _TargetedExperiment_Interpretation(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TargetedExperiment_Interpretation](new _TargetedExperiment_Interpretation())
    
    def _init_1(self, TargetedExperiment_Interpretation in_0 ):
        """
        _init_1(self, in_0: TargetedExperiment_Interpretation ) -> None
        """
        assert isinstance(in_0, TargetedExperiment_Interpretation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TargetedExperiment_Interpretation](new _TargetedExperiment_Interpretation((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TargetedExperiment_Interpretation ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TargetedExperiment_Interpretation)):
             self._init_1(*args)
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
    
    def consumeCVTerms(self, dict cv_term_map ):
        """
        consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366645
        cdef CVTerm item_1268773805887521753366645
        cdef bytes inner_key_1268773805887521753366645
        cdef list inner_values_1268773805887521753366645
        while outer_it_1268773805887521753366645 != _r.end():
           inner_key_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).first.c_str()
           inner_values_1268773805887521753366645 = []
           inner_it_1268773805887521753366645 = deref(outer_it_1268773805887521753366645).second.begin()
           while inner_it_1268773805887521753366645 != deref(outer_it_1268773805887521753366645).second.end():
               item_1268773805887521753366645 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366645)))
               inner_values_1268773805887521753366645.append(item_1268773805887521753366645)
               inc(inner_it_1268773805887521753366645)
           py_result[inner_key_1268773805887521753366645] = inner_values_1268773805887521753366645
           inc(outer_it_1268773805887521753366645)
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
    
    def getKeys(self, list keys ):
        """
        getKeys(self, keys: List[bytes] ) -> None
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
    
    def getKeysAsIntegers(self, list keys ):
        """
        getKeysAsIntegers(self, keys: List[int] ) -> None
        """
        assert isinstance(keys, list) and all(isinstance(elemt_rec, int) for elemt_rec in keys), 'arg keys wrong type'
        cdef libcpp_vector[unsigned int] v0 = keys
        self.inst.get().getKeys(v0)
        keys[:] = v0
    
    def _getMetaValue_0(self,  in_0 ):
        """
        _getMetaValue_0(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        cdef _DataValue _r = self.inst.get().getMetaValue((<unsigned int>in_0))
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
    
    def _getMetaValue_1(self,  in_0 ):
        """
        _getMetaValue_1(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
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
    
    def getMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMetaValue(self, in_0: int ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int)):
            return self._getMetaValue_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setMetaValue_0(self,  in_0 ,  in_1 ):
        """
        _setMetaValue_0(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue((<unsigned int>in_0), deref(DataValue(in_1).inst.get()))
    
    def _setMetaValue_1(self,  in_0 ,  in_1 ):
        """
        _setMetaValue_1(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, (int, float, list, bytes, str)), 'arg in_1 wrong type'
    
    
        self.inst.get().setMetaValue(deref((convString(in_0)).get()), deref(DataValue(in_1).inst.get()))
    
    def setMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMetaValue(self, in_0: int , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], int)) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setMetaValue_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], (int, float, list, bytes, str))):
            return self._setMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _metaValueExists_0(self,  in_0 ):
        """
        _metaValueExists_0(self, in_0: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists(deref((convString(in_0)).get()))
        py_result = <bool>_r
        return py_result
    
    def _metaValueExists_1(self,  in_0 ):
        """
        _metaValueExists_1(self, in_0: int ) -> bool
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        cdef bool _r = self.inst.get().metaValueExists((<unsigned int>in_0))
        py_result = <bool>_r
        return py_result
    
    def metaValueExists(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: metaValueExists(self, in_0: int ) -> bool
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._metaValueExists_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._metaValueExists_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _removeMetaValue_0(self,  in_0 ):
        """
        _removeMetaValue_0(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue(deref((convString(in_0)).get()))
    
    def _removeMetaValue_1(self,  in_0 ):
        """
        _removeMetaValue_1(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().removeMetaValue((<unsigned int>in_0))
    
    def removeMetaValue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: removeMetaValue(self, in_0: int ) -> None
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._removeMetaValue_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._removeMetaValue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, TargetedExperiment_Interpretation):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TargetedExperiment_Interpretation other_casted = other
        cdef TargetedExperiment_Interpretation self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class TargetedExperiment_Modification:
    """
    Cython implementation of _TargetedExperiment_Modification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TargetedExperiment_Modification.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property avg_mass_delta:
        def __set__(self, double avg_mass_delta):
        
            self.inst.get().avg_mass_delta = (<double>avg_mass_delta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().avg_mass_delta
            py_result = <double>_r
            return py_result
    
    property mono_mass_delta:
        def __set__(self, double mono_mass_delta):
        
            self.inst.get().mono_mass_delta = (<double>mono_mass_delta)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mono_mass_delta
            py_result = <double>_r
            return py_result
    
    property location:
        def __set__(self,  location):
        
            self.inst.get().location = (<int>location)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().location
            py_result = <int>_r
            return py_result
    
    property unimod_id:
        def __set__(self,  unimod_id):
        
            self.inst.get().unimod_id = (<int>unimod_id)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().unimod_id
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef TargetedExperiment_Modification rv = TargetedExperiment_Modification.__new__(TargetedExperiment_Modification)
       rv.inst = shared_ptr[_TargetedExperiment_Modification](new _TargetedExperiment_Modification(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TargetedExperiment_Modification rv = TargetedExperiment_Modification.__new__(TargetedExperiment_Modification)
       rv.inst = shared_ptr[_TargetedExperiment_Modification](new _TargetedExperiment_Modification(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TargetedExperiment_Modification](new _TargetedExperiment_Modification())
    
    def _init_1(self, TargetedExperiment_Modification in_0 ):
        """
        _init_1(self, in_0: TargetedExperiment_Modification ) -> None
        """
        assert isinstance(in_0, TargetedExperiment_Modification), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TargetedExperiment_Modification](new _TargetedExperiment_Modification((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TargetedExperiment_Modification ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TargetedExperiment_Modification)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, TargetedExperiment_Modification):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TargetedExperiment_Modification other_casted = other
        cdef TargetedExperiment_Modification self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class TheoreticalSpectrumGenerator:
    """
    Cython implementation of _TheoreticalSpectrumGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TheoreticalSpectrumGenerator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TheoreticalSpectrumGenerator rv = TheoreticalSpectrumGenerator.__new__(TheoreticalSpectrumGenerator)
       rv.inst = shared_ptr[_TheoreticalSpectrumGenerator](new _TheoreticalSpectrumGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TheoreticalSpectrumGenerator rv = TheoreticalSpectrumGenerator.__new__(TheoreticalSpectrumGenerator)
       rv.inst = shared_ptr[_TheoreticalSpectrumGenerator](new _TheoreticalSpectrumGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TheoreticalSpectrumGenerator](new _TheoreticalSpectrumGenerator())
    
    def _init_1(self, TheoreticalSpectrumGenerator in_0 ):
        """
        _init_1(self, in_0: TheoreticalSpectrumGenerator ) -> None
        """
        assert isinstance(in_0, TheoreticalSpectrumGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TheoreticalSpectrumGenerator](new _TheoreticalSpectrumGenerator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TheoreticalSpectrumGenerator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TheoreticalSpectrumGenerator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSpectrum(self, MSSpectrum spec , AASequence peptide ,  min_charge ,  max_charge ):
        """
        getSpectrum(self, spec: MSSpectrum , peptide: AASequence , min_charge: int , max_charge: int ) -> None
        Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters. If precursor_charge is set to 0 max_charge + 1 will be used
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(min_charge, int), 'arg min_charge wrong type'
        assert isinstance(max_charge, int), 'arg max_charge wrong type'
    
    
    
    
        self.inst.get().getSpectrum((deref(spec.inst.get())), (deref(peptide.inst.get())), (<int>min_charge), (<int>max_charge))
    
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

cdef class TraMLProduct:
    """
    Cython implementation of _TraMLProduct

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::TargetedExperimentHelper::RetentionTime::RTUnit_1_1TraMLProduct.html>`_
      -- Inherits from ['CVTermList']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TraMLProduct rv = TraMLProduct.__new__(TraMLProduct)
       rv.inst = shared_ptr[_TraMLProduct](new _TraMLProduct(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TraMLProduct rv = TraMLProduct.__new__(TraMLProduct)
       rv.inst = shared_ptr[_TraMLProduct](new _TraMLProduct(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TraMLProduct](new _TraMLProduct())
    
    def _init_1(self, TraMLProduct in_0 ):
        """
        _init_1(self, in_0: TraMLProduct ) -> None
        """
        assert isinstance(in_0, TraMLProduct), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TraMLProduct](new _TraMLProduct((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TraMLProduct ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TraMLProduct)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMZ(self, double mz ):
        """
        setMZ(self, mz: float ) -> None
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        self.inst.get().setMZ((<double>mz))
    
    def getMZ(self):
        """
        getMZ(self) -> float
        """
        cdef double _r = self.inst.get().getMZ()
        py_result = <double>_r
        return py_result
    
    def setChargeState(self,  charge ):
        """
        setChargeState(self, charge: int ) -> None
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setChargeState((<int>charge))
    
    def getChargeState(self):
        """
        getChargeState(self) -> int
        """
        cdef int _r = self.inst.get().getChargeState()
        py_result = <int>_r
        return py_result
    
    def hasCharge(self):
        """
        hasCharge(self) -> bool
        """
        cdef bool _r = self.inst.get().hasCharge()
        py_result = <bool>_r
        return py_result
    
    def getConfigurationList(self):
        """
        getConfigurationList(self) -> List[Configuration]
        """
        _r = self.inst.get().getConfigurationList()
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
    
    def getInterpretationList(self):
        """
        getInterpretationList(self) -> List[TargetedExperiment_Interpretation]
        """
        _r = self.inst.get().getInterpretationList()
        py_result = []
        cdef libcpp_vector[_TargetedExperiment_Interpretation].iterator it__r = _r.begin()
        cdef TargetedExperiment_Interpretation item_py_result
        while it__r != _r.end():
           item_py_result = TargetedExperiment_Interpretation.__new__(TargetedExperiment_Interpretation)
           item_py_result.inst = shared_ptr[_TargetedExperiment_Interpretation](new _TargetedExperiment_Interpretation(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addInterpretation(self, TargetedExperiment_Interpretation interpretation ):
        """
        addInterpretation(self, interpretation: TargetedExperiment_Interpretation ) -> None
        """
        assert isinstance(interpretation, TargetedExperiment_Interpretation), 'arg interpretation wrong type'
    
        self.inst.get().addInterpretation((deref(interpretation.inst.get())))
    
    def resetInterpretations(self):
        """
        resetInterpretations(self) -> None
        """
        self.inst.get().resetInterpretations()
    
    def setCVTerms(self, list terms ):
        """
        setCVTerms(self, terms: List[CVTerm] ) -> None
        Sets the CV terms
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
        Replaces the specified CV term
        """
        assert isinstance(term, CVTerm), 'arg term wrong type'
    
        self.inst.get().replaceCVTerm((deref(term.inst.get())))
    
    def replaceCVTerms(self, list cv_terms ,  accession ):
        """
        replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
        cdef libcpp_vector[_CVTerm] * v0 = new libcpp_vector[_CVTerm]()
        cdef CVTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().replaceCVTerms(deref(v0), deref((convString(accession)).get()))
        del v0
    
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
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> Dict[bytes,List[CVTerm]]
        Returns the accession string of the term
        """
        _r = self.inst.get().getCVTerms()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773805887521753366646 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773805887521753366646
        cdef CVTerm item_1268773805887521753366646
        cdef bytes inner_key_1268773805887521753366646
        cdef list inner_values_1268773805887521753366646
        while outer_it_1268773805887521753366646 != _r.end():
           inner_key_1268773805887521753366646 = deref(outer_it_1268773805887521753366646).first.c_str()
           inner_values_1268773805887521753366646 = []
           inner_it_1268773805887521753366646 = deref(outer_it_1268773805887521753366646).second.begin()
           while inner_it_1268773805887521753366646 != deref(outer_it_1268773805887521753366646).second.end():
               item_1268773805887521753366646 = CVTerm.__new__(CVTerm)
               item_1268773805887521753366646.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773805887521753366646)))
               inner_values_1268773805887521753366646.append(item_1268773805887521753366646)
               inc(inner_it_1268773805887521753366646)
           py_result[inner_key_1268773805887521753366646] = inner_values_1268773805887521753366646
           inc(outer_it_1268773805887521753366646)
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
        if not isinstance(other, TraMLProduct):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TraMLProduct other_casted = other
        cdef TraMLProduct self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 
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
