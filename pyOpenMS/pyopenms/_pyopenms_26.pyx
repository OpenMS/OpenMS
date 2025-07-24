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

cdef class AnnotationState:
    None
    FEATURE_ID_NONE = 0
    FEATURE_ID_SINGLE = 1
    FEATURE_ID_MULTIPLE_SAME = 2
    FEATURE_ID_MULTIPLE_DIVERGENT = 3
    SIZE_OF_ANNOTATIONSTATE = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class Attachment:
    """
    Cython implementation of _Attachment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::QcMLFile_1_1Attachment.html>`_
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
    
    property value:
        def __set__(self,  value):
        
            self.inst.get().value = deref((convString(value)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().value
            py_result = convOutputString(_r)
            return py_result
    
    property cvRef:
        def __set__(self,  cvRef):
        
            self.inst.get().cvRef = deref((convString(cvRef)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().cvRef
            py_result = convOutputString(_r)
            return py_result
    
    property cvAcc:
        def __set__(self,  cvAcc):
        
            self.inst.get().cvAcc = deref((convString(cvAcc)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().cvAcc
            py_result = convOutputString(_r)
            return py_result
    
    property unitRef:
        def __set__(self,  unitRef):
        
            self.inst.get().unitRef = deref((convString(unitRef)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().unitRef
            py_result = convOutputString(_r)
            return py_result
    
    property unitAcc:
        def __set__(self,  unitAcc):
        
            self.inst.get().unitAcc = deref((convString(unitAcc)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().unitAcc
            py_result = convOutputString(_r)
            return py_result
    
    property binary:
        def __set__(self,  binary):
        
            self.inst.get().binary = deref((convString(binary)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().binary
            py_result = convOutputString(_r)
            return py_result
    
    property qualityRef:
        def __set__(self,  qualityRef):
        
            self.inst.get().qualityRef = deref((convString(qualityRef)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().qualityRef
            py_result = convOutputString(_r)
            return py_result
    
    property colTypes:
        def __set__(self, list colTypes):
            cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
            cdef bytes item0
            for item0 in colTypes:
               v0.push_back(_String(<char *>item0))
            self.inst.get().colTypes = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().colTypes
            py_result = []
            cdef libcpp_vector[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.append(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef Attachment rv = Attachment.__new__(Attachment)
       rv.inst = shared_ptr[_Attachment](new _Attachment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Attachment rv = Attachment.__new__(Attachment)
       rv.inst = shared_ptr[_Attachment](new _Attachment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Attachment](new _Attachment())
    
    def _init_1(self, Attachment in_0 ):
        """
        _init_1(self, in_0: Attachment ) -> None
        """
        assert isinstance(in_0, Attachment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Attachment](new _Attachment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Attachment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Attachment)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toXMLString(self,  indentation_level ):
        """
        toXMLString(self, indentation_level: int ) -> Union[bytes, str, String]
        """
        assert isinstance(indentation_level, int), 'arg indentation_level wrong type'
    
        cdef _String _r = self.inst.get().toXMLString((<unsigned int>indentation_level))
        py_result = convOutputString(_r)
        return py_result
    
    def toCSVString(self,  separator ):
        """
        toCSVString(self, separator: Union[bytes, str, String] ) -> Union[bytes, str, String]
        """
        assert (isinstance(separator, str) or isinstance(separator, bytes) or isinstance(separator, String)), 'arg separator wrong type'
    
        cdef _String _r = self.inst.get().toCSVString(deref((convString(separator)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 0, 4):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Attachment):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Attachment other_casted = other
        cdef Attachment self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get())
        if op==4:
            return deref(self_casted.inst.get()) > deref(other_casted.inst.get())
    

    property tableRows:
    
        def __set__(self, list tableRows):
            # libcpp_vector[ libcpp_vector[String] ]
            cdef libcpp_vector[libcpp_vector[_String]] * v0 = new libcpp_vector[libcpp_vector[_String]]()
            cdef libcpp_vector[_String] * v0_rec = new libcpp_vector[_String]()
            cdef bytes item0
            for tableRows_rec in tableRows:
                v0_rec.clear()
                for item0 in tableRows_rec:
                    v0_rec.push_back(_String(<char *>item0))
                v0.push_back(deref(v0_rec))
            self.inst.get().tableRows = deref(v0)
            del v0
            del v0_rec
    

        def __get__(self):
            c_res = self.inst.get().tableRows
            # libcpp_vector[ libcpp_vector[String] ]
            cdef replace = []
            cdef libcpp_vector[libcpp_vector[_String] ].iterator it = c_res.begin()
            cdef libcpp_vector[_String ].iterator it_inner
            while it != c_res.end():
                replace_inner = []
                it_inner = deref(it).begin()
                while it_inner != deref(it).end():
                    replace_inner.append(<char*>deref(it_inner).c_str())
                    inc(it_inner)
                replace.append(replace_inner)
                inc(it)
            return replace 

cdef class BaseFeature:
    """
    Cython implementation of _BaseFeature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BaseFeature.html>`_
      -- Inherits from ['UniqueIdInterface', 'RichPeak2D']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef BaseFeature rv = BaseFeature.__new__(BaseFeature)
       rv.inst = shared_ptr[_BaseFeature](new _BaseFeature(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef BaseFeature rv = BaseFeature.__new__(BaseFeature)
       rv.inst = shared_ptr[_BaseFeature](new _BaseFeature(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BaseFeature](new _BaseFeature())
    
    def _init_1(self, BaseFeature in_0 ):
        """
        _init_1(self, in_0: BaseFeature ) -> None
        """
        assert isinstance(in_0, BaseFeature), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BaseFeature](new _BaseFeature((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: BaseFeature ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], BaseFeature)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def _setPeptideIdentifications_0(self, PeptideIdentificationList peptides ):
        """
        _setPeptideIdentifications_0(self, peptides: PeptideIdentificationList ) -> None
        Sets the PeptideIdentification vector
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().setPeptideIdentifications((deref(peptides.inst.get())))
    
    def _setPeptideIdentifications_1(self, PeptideIdentificationList peptides ):
        """
        _setPeptideIdentifications_1(self, peptides: PeptideIdentificationList ) -> None
        Sets the PeptideIdentificationList
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().setPeptideIdentifications((deref(peptides.inst.get())))
    
    def setPeptideIdentifications(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None
          :noindex:
        
        Sets the PeptideIdentification vector

        
        .. rubric:: Overload:
        .. py:function:: setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None
          :noindex:
        
        Sets the PeptideIdentificationList
    
        """
        if (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._setPeptideIdentifications_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._setPeptideIdentifications_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
        if not isinstance(other, BaseFeature):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef BaseFeature other_casted = other
        cdef BaseFeature self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class CachedSwathFileConsumer:
    """
    Cython implementation of _CachedSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CachedSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CachedSwathFileConsumer rv = CachedSwathFileConsumer.__new__(CachedSwathFileConsumer)
       rv.inst = shared_ptr[_CachedSwathFileConsumer](new _CachedSwathFileConsumer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CachedSwathFileConsumer rv = CachedSwathFileConsumer.__new__(CachedSwathFileConsumer)
       rv.inst = shared_ptr[_CachedSwathFileConsumer](new _CachedSwathFileConsumer(deref(self.inst.get())))
       return rv
    
    def _init_0(self, CachedSwathFileConsumer in_0 ):
        """
        _init_0(self, in_0: CachedSwathFileConsumer ) -> None
        """
        assert isinstance(in_0, CachedSwathFileConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CachedSwathFileConsumer](new _CachedSwathFileConsumer((deref(in_0.inst.get()))))
    
    def _init_1(self,  cachedir ,  basename ,  nr_ms1_spectra , list nr_ms2_spectra ):
        """
        _init_1(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
        """
        assert (isinstance(cachedir, str) or isinstance(cachedir, bytes) or isinstance(cachedir, String)), 'arg cachedir wrong type'
        assert (isinstance(basename, str) or isinstance(basename, bytes) or isinstance(basename, String)), 'arg basename wrong type'
        assert isinstance(nr_ms1_spectra, int) and nr_ms1_spectra >= 0, 'arg nr_ms1_spectra wrong type'
        assert isinstance(nr_ms2_spectra, list) and all(isinstance(elemt_rec, int) for elemt_rec in nr_ms2_spectra), 'arg nr_ms2_spectra wrong type'
    
    
    
        cdef libcpp_vector[int] v3 = nr_ms2_spectra
        self.inst = shared_ptr[_CachedSwathFileConsumer](new _CachedSwathFileConsumer(deref((convString(cachedir)).get()), deref((convString(basename)).get()), (<size_t>nr_ms1_spectra), v3))
        
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CachedSwathFileConsumer ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], CachedSwathFileConsumer)):
             self._init_0(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[3])):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setExpectedSize(self,  s ,  c ):
        """
        setExpectedSize(self, s: int , c: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
        assert isinstance(c, int) and c >= 0, 'arg c wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>s), (<size_t>c))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def retrieveSwathMaps(self, list maps ):
        """
        retrieveSwathMaps(self, maps: List[SwathMap] ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in maps), 'arg maps wrong type'
        cdef libcpp_vector[_SwathMap] * v0 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().retrieveSwathMaps(deref(v0))
        cdef libcpp_vector[_SwathMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = SwathMap.__new__(SwathMap)
            item0.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
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

cdef class EmgModel:
    """
    Cython implementation of _EmgModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgModel.html>`_
      -- Inherits from ['InterpolationModel']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EmgModel rv = EmgModel.__new__(EmgModel)
       rv.inst = shared_ptr[_EmgModel](new _EmgModel(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EmgModel rv = EmgModel.__new__(EmgModel)
       rv.inst = shared_ptr[_EmgModel](new _EmgModel(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Exponentially modified gaussian distribution model for elution profiles
        """
        self.inst = shared_ptr[_EmgModel](new _EmgModel())
    
    def _init_1(self, EmgModel in_0 ):
        """
        _init_1(self, in_0: EmgModel ) -> None
        """
        assert isinstance(in_0, EmgModel), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmgModel](new _EmgModel((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Exponentially modified gaussian distribution model for elution profiles

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EmgModel ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EmgModel)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIntensity(self, double coord ):
        """
        getIntensity(self, coord: float ) -> float
        Access model predicted intensity at position 'pos'
        """
        assert isinstance(coord, float), 'arg coord wrong type'
    
        cdef double _r = self.inst.get().getIntensity((<double>coord))
        py_result = <double>_r
        return py_result
    
    def getScalingFactor(self):
        """
        getScalingFactor(self) -> float
        Returns the interpolation class
        """
        cdef double _r = self.inst.get().getScalingFactor()
        py_result = <double>_r
        return py_result
    
    def setOffset(self, double offset ):
        """
        setOffset(self, offset: float ) -> None
        Sets the offset of the model
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset((<double>offset))
    
    def getCenter(self):
        """
        getCenter(self) -> float
        Returns the "center" of the model, particular definition (depends on the derived model)
        """
        cdef double _r = self.inst.get().getCenter()
        py_result = <double>_r
        return py_result
    
    def setSamples(self):
        """
        setSamples(self) -> None
        Sets sample/supporting points of interpolation wrt params
        """
        self.inst.get().setSamples()
    
    def setInterpolationStep(self, double interpolation_step ):
        """
        setInterpolationStep(self, interpolation_step: float ) -> None
        Sets the interpolation step for the linear interpolation of the model
        """
        assert isinstance(interpolation_step, float), 'arg interpolation_step wrong type'
    
        self.inst.get().setInterpolationStep((<double>interpolation_step))
    
    def setScalingFactor(self, double scaling ):
        """
        setScalingFactor(self, scaling: float ) -> None
        Sets the scaling factor of the model
        """
        assert isinstance(scaling, float), 'arg scaling wrong type'
    
        self.inst.get().setScalingFactor((<double>scaling))
    
    def getInterpolation(self):
        """
        getInterpolation(self) -> LinearInterpolation
        Returns the interpolation class
        """
        cdef _LinearInterpolation[double,double] * _r = new _LinearInterpolation[double,double](self.inst.get().getInterpolation())
        cdef LinearInterpolation py_result = LinearInterpolation.__new__(LinearInterpolation)
        py_result.inst = shared_ptr[_LinearInterpolation[double,double]](_r)
        return py_result 

cdef class GaussTraceFitter:
    """
    Cython implementation of _GaussTraceFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GaussTraceFitter.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GaussTraceFitter rv = GaussTraceFitter.__new__(GaussTraceFitter)
       rv.inst = shared_ptr[_GaussTraceFitter](new _GaussTraceFitter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GaussTraceFitter rv = GaussTraceFitter.__new__(GaussTraceFitter)
       rv.inst = shared_ptr[_GaussTraceFitter](new _GaussTraceFitter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Fitter for RT profiles using a Gaussian background model
        """
        self.inst = shared_ptr[_GaussTraceFitter](new _GaussTraceFitter())
    
    def _init_1(self, GaussTraceFitter in_0 ):
        """
        _init_1(self, in_0: GaussTraceFitter ) -> None
        """
        assert isinstance(in_0, GaussTraceFitter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GaussTraceFitter](new _GaussTraceFitter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Fitter for RT profiles using a Gaussian background model

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GaussTraceFitter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], GaussTraceFitter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def fit(self, MassTraces traces ):
        """
        fit(self, traces: MassTraces ) -> None
        Override important methods
        """
        assert isinstance(traces, MassTraces), 'arg traces wrong type'
    
        self.inst.get().fit((deref(traces.inst.get())))
    
    def getLowerRTBound(self):
        """
        getLowerRTBound(self) -> float
        Returns the lower RT bound
        """
        cdef double _r = self.inst.get().getLowerRTBound()
        py_result = <double>_r
        return py_result
    
    def getUpperRTBound(self):
        """
        getUpperRTBound(self) -> float
        Returns the upper RT bound
        """
        cdef double _r = self.inst.get().getUpperRTBound()
        py_result = <double>_r
        return py_result
    
    def getHeight(self):
        """
        getHeight(self) -> float
        Returns height of the fitted gaussian model
        """
        cdef double _r = self.inst.get().getHeight()
        py_result = <double>_r
        return py_result
    
    def getCenter(self):
        """
        getCenter(self) -> float
        Returns center of the fitted gaussian model
        """
        cdef double _r = self.inst.get().getCenter()
        py_result = <double>_r
        return py_result
    
    def getFWHM(self):
        """
        getFWHM(self) -> float
        Returns FWHM of the fitted gaussian model
        """
        cdef double _r = self.inst.get().getFWHM()
        py_result = <double>_r
        return py_result
    
    def getSigma(self):
        """
        getSigma(self) -> float
        Returns Sigma of the fitted gaussian model
        """
        cdef double _r = self.inst.get().getSigma()
        py_result = <double>_r
        return py_result
    
    def checkMaximalRTSpan(self, double max_rt_span ):
        """
        checkMaximalRTSpan(self, max_rt_span: float ) -> bool
        """
        assert isinstance(max_rt_span, float), 'arg max_rt_span wrong type'
    
        cdef bool _r = self.inst.get().checkMaximalRTSpan((<double>max_rt_span))
        py_result = <bool>_r
        return py_result
    
    def checkMinimalRTSpan(self, list rt_bounds , double min_rt_span ):
        """
        checkMinimalRTSpan(self, rt_bounds: List[float, float] , min_rt_span: float ) -> bool
        """
        assert isinstance(rt_bounds, list) and len(rt_bounds) == 2 and isinstance(rt_bounds[0], float) and isinstance(rt_bounds[1], float), 'arg rt_bounds wrong type'
        assert isinstance(min_rt_span, float), 'arg min_rt_span wrong type'
        cdef libcpp_pair[double, double] v0
        v0.first = rt_bounds[0]
        v0.second = rt_bounds[1]
    
        cdef bool _r = self.inst.get().checkMinimalRTSpan(v0, (<double>min_rt_span))
        rt_bounds[:] = [v0.first, v0.second]
        py_result = <bool>_r
        return py_result
    
    def computeTheoretical(self, MassTrace trace ,  k ):
        """
        computeTheoretical(self, trace: MassTrace , k: int ) -> float
        """
        assert isinstance(trace, MassTrace), 'arg trace wrong type'
        assert isinstance(k, int) and k >= 0, 'arg k wrong type'
    
    
        cdef double _r = self.inst.get().computeTheoretical((deref(trace.inst.get())), (<size_t>k))
        py_result = <double>_r
        return py_result
    
    def getArea(self):
        """
        getArea(self) -> float
        Returns area of the fitted gaussian model
        """
        cdef double _r = self.inst.get().getArea()
        py_result = <double>_r
        return py_result
    
    def getGnuplotFormula(self, MassTrace trace , bytes function_name , double baseline , double rt_shift ):
        """
        getGnuplotFormula(self, trace: MassTrace , function_name: bytes , baseline: float , rt_shift: float ) -> Union[bytes, str, String]
        """
        assert isinstance(trace, MassTrace), 'arg trace wrong type'
        assert isinstance(function_name, bytes) and len(function_name) == 1, 'arg function_name wrong type'
        assert isinstance(baseline, float), 'arg baseline wrong type'
        assert isinstance(rt_shift, float), 'arg rt_shift wrong type'
    
    
    
    
        cdef _String _r = self.inst.get().getGnuplotFormula((deref(trace.inst.get())), (<char>((function_name)[0])), (<double>baseline), (<double>rt_shift))
        py_result = convOutputString(_r)
        return py_result
    
    def getValue(self, double rt ):
        """
        getValue(self, rt: float ) -> float
        Returns value of the fitted gaussian model
        """
        assert isinstance(rt, float), 'arg rt wrong type'
    
        cdef double _r = self.inst.get().getValue((<double>rt))
        py_result = <double>_r
        return py_result 

cdef class InterpolationModel:
    """
    Cython implementation of _InterpolationModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InterpolationModel.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef InterpolationModel rv = InterpolationModel.__new__(InterpolationModel)
       rv.inst = shared_ptr[_InterpolationModel](new _InterpolationModel(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InterpolationModel rv = InterpolationModel.__new__(InterpolationModel)
       rv.inst = shared_ptr[_InterpolationModel](new _InterpolationModel(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Abstract class for 1D-models that are approximated using linear interpolation
        """
        self.inst = shared_ptr[_InterpolationModel](new _InterpolationModel())
    
    def _init_1(self, InterpolationModel in_0 ):
        """
        _init_1(self, in_0: InterpolationModel ) -> None
        """
        assert isinstance(in_0, InterpolationModel), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InterpolationModel](new _InterpolationModel((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Abstract class for 1D-models that are approximated using linear interpolation

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InterpolationModel ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InterpolationModel)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIntensity(self, double coord ):
        """
        getIntensity(self, coord: float ) -> float
        Access model predicted intensity at position 'pos'
        """
        assert isinstance(coord, float), 'arg coord wrong type'
    
        cdef double _r = self.inst.get().getIntensity((<double>coord))
        py_result = <double>_r
        return py_result
    
    def getScalingFactor(self):
        """
        getScalingFactor(self) -> float
        Returns the interpolation class
        """
        cdef double _r = self.inst.get().getScalingFactor()
        py_result = <double>_r
        return py_result
    
    def setOffset(self, double offset ):
        """
        setOffset(self, offset: float ) -> None
        Sets the offset of the model
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset((<double>offset))
    
    def getCenter(self):
        """
        getCenter(self) -> float
        Returns the "center" of the model, particular definition (depends on the derived model)
        """
        cdef double _r = self.inst.get().getCenter()
        py_result = <double>_r
        return py_result
    
    def setSamples(self):
        """
        setSamples(self) -> None
        Sets sample/supporting points of interpolation wrt params
        """
        self.inst.get().setSamples()
    
    def setInterpolationStep(self, double interpolation_step ):
        """
        setInterpolationStep(self, interpolation_step: float ) -> None
        Sets the interpolation step for the linear interpolation of the model
        """
        assert isinstance(interpolation_step, float), 'arg interpolation_step wrong type'
    
        self.inst.get().setInterpolationStep((<double>interpolation_step))
    
    def setScalingFactor(self, double scaling ):
        """
        setScalingFactor(self, scaling: float ) -> None
        Sets the scaling factor of the model
        """
        assert isinstance(scaling, float), 'arg scaling wrong type'
    
        self.inst.get().setScalingFactor((<double>scaling))
    
    def getInterpolation(self):
        """
        getInterpolation(self) -> LinearInterpolation
        Returns the interpolation class
        """
        cdef _LinearInterpolation[double,double] * _r = new _LinearInterpolation[double,double](self.inst.get().getInterpolation())
        cdef LinearInterpolation py_result = LinearInterpolation.__new__(LinearInterpolation)
        py_result.inst = shared_ptr[_LinearInterpolation[double,double]](_r)
        return py_result 

cdef class MassDecompositionAlgorithm:
    """
    Cython implementation of _MassDecompositionAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassDecompositionAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MassDecompositionAlgorithm](new _MassDecompositionAlgorithm())
    
    def getDecompositions(self, list decomps , double weight ):
        """
        getDecompositions(self, decomps: List[MassDecomposition] , weight: float ) -> None
        Returns the possible decompositions given the weight
        """
        assert isinstance(decomps, list) and all(isinstance(elemt_rec, MassDecomposition) for elemt_rec in decomps), 'arg decomps wrong type'
        assert isinstance(weight, float), 'arg weight wrong type'
        cdef libcpp_vector[_MassDecomposition] * v0 = new libcpp_vector[_MassDecomposition]()
        cdef MassDecomposition item0
        for item0 in decomps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().getDecompositions(deref(v0), (<double>weight))
        cdef libcpp_vector[_MassDecomposition].iterator it_decomps = v0.begin()
        replace_0 = []
        while it_decomps != v0.end():
            item0 = MassDecomposition.__new__(MassDecomposition)
            item0.inst = shared_ptr[_MassDecomposition](new _MassDecomposition(deref(it_decomps)))
            replace_0.append(item0)
            inc(it_decomps)
        decomps[:] = replace_0
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

cdef class MassTraceDetection:
    """
    Cython implementation of _MassTraceDetection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassTraceDetection.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MassTraceDetection rv = MassTraceDetection.__new__(MassTraceDetection)
       rv.inst = shared_ptr[_MassTraceDetection](new _MassTraceDetection(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassTraceDetection rv = MassTraceDetection.__new__(MassTraceDetection)
       rv.inst = shared_ptr[_MassTraceDetection](new _MassTraceDetection(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MassTraceDetection](new _MassTraceDetection())
    
    def _init_1(self, MassTraceDetection in_0 ):
        """
        _init_1(self, in_0: MassTraceDetection ) -> None
        """
        assert isinstance(in_0, MassTraceDetection), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassTraceDetection](new _MassTraceDetection((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassTraceDetection ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassTraceDetection)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment input_map , list traces ,  max_traces ):
        """
        run(self, input_map: MSExperiment , traces: List[Kernel_MassTrace] , max_traces: int ) -> None
        """
        assert isinstance(input_map, MSExperiment), 'arg input_map wrong type'
        assert isinstance(traces, list) and all(isinstance(elemt_rec, Kernel_MassTrace) for elemt_rec in traces), 'arg traces wrong type'
        assert isinstance(max_traces, int) and max_traces >= 0, 'arg max_traces wrong type'
    
        cdef libcpp_vector[_Kernel_MassTrace] * v1 = new libcpp_vector[_Kernel_MassTrace]()
        cdef Kernel_MassTrace item1
        for item1 in traces:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().run((deref(input_map.inst.get())), deref(v1), (<size_t>max_traces))
        cdef libcpp_vector[_Kernel_MassTrace].iterator it_traces = v1.begin()
        replace_0 = []
        while it_traces != v1.end():
            item1 = Kernel_MassTrace.__new__(Kernel_MassTrace)
            item1.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(it_traces)))
            replace_0.append(item1)
            inc(it_traces)
        traces[:] = replace_0
        del v1
    
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

cdef class MetaInfoInterface:
    """
    Cython implementation of _MetaInfoInterface

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoInterface.html>`_

    Interface for classes that can store arbitrary meta information
    (Type-Name-Value tuples).
    
    MetaInfoInterface is a base class for all classes that use one MetaInfo
    object as member.  If you want to add meta information to a class, let it
    publicly inherit the MetaInfoInterface.  Meta information is an array of
    Type-Name-Value tuples.
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaInfoInterface rv = MetaInfoInterface.__new__(MetaInfoInterface)
       rv.inst = shared_ptr[_MetaInfoInterface](new _MetaInfoInterface(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaInfoInterface rv = MetaInfoInterface.__new__(MetaInfoInterface)
       rv.inst = shared_ptr[_MetaInfoInterface](new _MetaInfoInterface(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaInfoInterface](new _MetaInfoInterface())
    
    def _init_1(self, MetaInfoInterface in_0 ):
        """
        _init_1(self, in_0: MetaInfoInterface ) -> None
        """
        assert isinstance(in_0, MetaInfoInterface), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaInfoInterface](new _MetaInfoInterface((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaInfoInterface ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaInfoInterface)):
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
        if not isinstance(other, MetaInfoInterface):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MetaInfoInterface other_casted = other
        cdef MetaInfoInterface self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    


    def getMetaValues(self):
        """Cython signature: dict getMetaValues()

        Returns all meta values in a Python dictionary.
        """
        # cdef _DataValue _c_value
        # cdef DataValue _value 
        # cdef int _type

        mmap = {}
        cdef libcpp_vector[_String] keys 
        cdef object py_result

        # get keys for iteration
        self.inst.get().getKeys(keys)

        cdef libcpp_vector[_String].iterator k_it = keys.begin()
        while k_it != keys.end():
            # easy approach: call Python fxn
            py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            py_result = self.getMetaValue(py_str)
            # # hard approach: do it ourselves
            # _c_value = self.inst.get().getMetaValue(deref(k_it))
            # py_str = _cast_const_away(<char*>(deref(k_it)).c_str())
            # _type = _c_value.valueType()
            # _value = DataValue.__new__(DataValue)
            # _value.inst = shared_ptr[_DataValue](new _DataValue(_c_value))
            # if _type == DataType.STRING_VALUE:
            #     py_result = _value.toString()
            # elif _type == DataType.INT_VALUE:
            #     py_result = _value.toInt()
            # elif _type == DataType.DOUBLE_VALUE:
            #     py_result = _value.toDouble()
            # elif _type == DataType.INT_LIST:
            #     py_result = _value.toIntList()
            # elif _type == DataType.DOUBLE_LIST:
            #     py_result = _value.toDoubleList()
            # elif _type == DataType.STRING_LIST:
            #     py_result = _value.toStringList()
            # elif _type == DataType.EMPTY_VALUE:
            #     py_result = None
            # else:
            #     raise Exception("DataValue instance has invalid value type %d" % _type)

            mmap[ py_str ] = py_result
            inc(k_it)

        return mmap

    def setMetaValues(self, dict mmap):
        """Cython signature: setMetaValues(dict values)

        Sets the meta values given in the Python dictionary.
        """
        self.inst.get().clearMetaInfo() # ensure its empty first
        for k, v in mmap.iteritems():
            self.inst.get().setMetaValue(deref((convString(k)).get()), deref(DataValue(v).inst.get()))
            # self.setMetaValue(k, v) 

cdef class MorphologicalFilter:
    """
    Cython implementation of _MorphologicalFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorphologicalFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MorphologicalFilter](new _MorphologicalFilter())
    
    def filter(self, MSSpectrum spectrum ):
        """
        filter(self, spectrum: MSSpectrum ) -> None
        Applies the morphological filtering operation to an MSSpectrum
        
        If the size of the structuring element is given in 'Thomson', the number of data points for
        the structuring element is computed as follows:
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().filter((deref(spectrum.inst.get())))
    
    def filterExperiment(self, MSExperiment exp ):
        """
        filterExperiment(self, exp: MSExperiment ) -> None
        Applies the morphological filtering operation to an MSExperiment
        
        The size of the structuring element is computed for each spectrum individually, if it is given in 'Thomson'
        See the filtering method for MSSpectrum for details
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().filterExperiment((deref(exp.inst.get())))
    
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

cdef class MzMLSwathFileConsumer:
    """
    Cython implementation of _MzMLSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self, MzMLSwathFileConsumer in_0 ):
        """
        _init_0(self, in_0: MzMLSwathFileConsumer ) -> None
        """
        assert isinstance(in_0, MzMLSwathFileConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzMLSwathFileConsumer](new _MzMLSwathFileConsumer((deref(in_0.inst.get()))))
    
    def _init_1(self,  cachedir ,  basename ,  nr_ms1_spectra , list nr_ms2_spectra ):
        """
        _init_1(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
        """
        assert (isinstance(cachedir, str) or isinstance(cachedir, bytes) or isinstance(cachedir, String)), 'arg cachedir wrong type'
        assert (isinstance(basename, str) or isinstance(basename, bytes) or isinstance(basename, String)), 'arg basename wrong type'
        assert isinstance(nr_ms1_spectra, int) and nr_ms1_spectra >= 0, 'arg nr_ms1_spectra wrong type'
        assert isinstance(nr_ms2_spectra, list) and all(isinstance(elemt_rec, int) for elemt_rec in nr_ms2_spectra), 'arg nr_ms2_spectra wrong type'
    
    
    
        cdef libcpp_vector[int] v3 = nr_ms2_spectra
        self.inst = shared_ptr[_MzMLSwathFileConsumer](new _MzMLSwathFileConsumer(deref((convString(cachedir)).get()), deref((convString(basename)).get()), (<size_t>nr_ms1_spectra), v3))
        
    
    def _init_2(self, list known_window_boundaries ,  cachedir ,  basename ,  nr_ms1_spectra , list nr_ms2_spectra ):
        """
        _init_2(self, known_window_boundaries: List[SwathMap] , cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
        """
        assert isinstance(known_window_boundaries, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in known_window_boundaries), 'arg known_window_boundaries wrong type'
        assert (isinstance(cachedir, str) or isinstance(cachedir, bytes) or isinstance(cachedir, String)), 'arg cachedir wrong type'
        assert (isinstance(basename, str) or isinstance(basename, bytes) or isinstance(basename, String)), 'arg basename wrong type'
        assert isinstance(nr_ms1_spectra, int) and nr_ms1_spectra >= 0, 'arg nr_ms1_spectra wrong type'
        assert isinstance(nr_ms2_spectra, list) and all(isinstance(elemt_rec, int) for elemt_rec in nr_ms2_spectra), 'arg nr_ms2_spectra wrong type'
        cdef libcpp_vector[_SwathMap] * v0 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item0
        for item0 in known_window_boundaries:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        cdef libcpp_vector[int] v4 = nr_ms2_spectra
        self.inst = shared_ptr[_MzMLSwathFileConsumer](new _MzMLSwathFileConsumer(deref(v0), deref((convString(cachedir)).get()), deref((convString(basename)).get()), (<size_t>nr_ms1_spectra), v4))
        
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzMLSwathFileConsumer ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, known_window_boundaries: List[SwathMap] , cachedir: Union[bytes, str, String] , basename: Union[bytes, str, String] , nr_ms1_spectra: int , nr_ms2_spectra: List[int] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], MzMLSwathFileConsumer)):
             self._init_0(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[3])):
             self._init_1(*args)
        elif (len(args)==5) and (isinstance(args[0], list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in args[0])) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], int) and args[3] >= 0) and (isinstance(args[4], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[4])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setExpectedSize(self,  s ,  c ):
        """
        setExpectedSize(self, s: int , c: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
        assert isinstance(c, int) and c >= 0, 'arg c wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>s), (<size_t>c))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def retrieveSwathMaps(self, list maps ):
        """
        retrieveSwathMaps(self, maps: List[SwathMap] ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in maps), 'arg maps wrong type'
        cdef libcpp_vector[_SwathMap] * v0 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().retrieveSwathMaps(deref(v0))
        cdef libcpp_vector[_SwathMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = SwathMap.__new__(SwathMap)
            item0.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
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

cdef class RegularSwathFileConsumer:
    """
    Cython implementation of _RegularSwathFileConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RegularSwathFileConsumer.html>`_
      -- Inherits from ['FullSwathFileConsumer']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RegularSwathFileConsumer rv = RegularSwathFileConsumer.__new__(RegularSwathFileConsumer)
       rv.inst = shared_ptr[_RegularSwathFileConsumer](new _RegularSwathFileConsumer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RegularSwathFileConsumer rv = RegularSwathFileConsumer.__new__(RegularSwathFileConsumer)
       rv.inst = shared_ptr[_RegularSwathFileConsumer](new _RegularSwathFileConsumer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RegularSwathFileConsumer](new _RegularSwathFileConsumer())
    
    def _init_1(self, RegularSwathFileConsumer in_0 ):
        """
        _init_1(self, in_0: RegularSwathFileConsumer ) -> None
        """
        assert isinstance(in_0, RegularSwathFileConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RegularSwathFileConsumer](new _RegularSwathFileConsumer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RegularSwathFileConsumer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RegularSwathFileConsumer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setExpectedSize(self,  s ,  c ):
        """
        setExpectedSize(self, s: int , c: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
        assert isinstance(c, int) and c >= 0, 'arg c wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>s), (<size_t>c))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def retrieveSwathMaps(self, list maps ):
        """
        retrieveSwathMaps(self, maps: List[SwathMap] ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in maps), 'arg maps wrong type'
        cdef libcpp_vector[_SwathMap] * v0 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().retrieveSwathMaps(deref(v0))
        cdef libcpp_vector[_SwathMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = SwathMap.__new__(SwathMap)
            item0.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
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

cdef class SimplePeak:
    """
    Cython implementation of _SimplePeak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimplePeak.html>`_
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
    
    property charge:
        def __set__(self,  charge):
        
            self.inst.get().charge = (<int>charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().charge
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef SimplePeak rv = SimplePeak.__new__(SimplePeak)
       rv.inst = shared_ptr[_SimplePeak](new _SimplePeak(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SimplePeak rv = SimplePeak.__new__(SimplePeak)
       rv.inst = shared_ptr[_SimplePeak](new _SimplePeak(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A simple struct to represent peaks with mz and charge and sort them easily
        """
        self.inst = shared_ptr[_SimplePeak](new _SimplePeak())
    
    def _init_1(self, double mz ,  charge ):
        """
        _init_1(self, mz: float , charge: int ) -> None
        """
        assert isinstance(mz, float), 'arg mz wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        self.inst = shared_ptr[_SimplePeak](new _SimplePeak((<double>mz), (<int>charge)))
    
    def _init_2(self, SimplePeak in_0 ):
        """
        _init_2(self, in_0: SimplePeak ) -> None
        """
        assert isinstance(in_0, SimplePeak), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SimplePeak](new _SimplePeak((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A simple struct to represent peaks with mz and charge and sort them easily

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mz: float , charge: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SimplePeak ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], int)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], SimplePeak)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SimpleTSGXLMS:
    """
    Cython implementation of _SimpleTSGXLMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimpleTSGXLMS.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SimpleTSGXLMS rv = SimpleTSGXLMS.__new__(SimpleTSGXLMS)
       rv.inst = shared_ptr[_SimpleTSGXLMS](new _SimpleTSGXLMS(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SimpleTSGXLMS rv = SimpleTSGXLMS.__new__(SimpleTSGXLMS)
       rv.inst = shared_ptr[_SimpleTSGXLMS](new _SimpleTSGXLMS(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Generates theoretical spectra for cross-linked peptides
        
        The spectra this class generates are vectors of SimplePeaks
        This class generates the same peak types as TheoreticalSpectrumGeneratorXLMS
        and the interface is very similar, but it is simpler and faster
        SimplePeak only contains an mz value and a charge. No intensity values
        or String annotations or other additional DataArrays are generated
        """
        self.inst = shared_ptr[_SimpleTSGXLMS](new _SimpleTSGXLMS())
    
    def _init_1(self, SimpleTSGXLMS in_0 ):
        """
        _init_1(self, in_0: SimpleTSGXLMS ) -> None
        """
        assert isinstance(in_0, SimpleTSGXLMS), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SimpleTSGXLMS](new _SimpleTSGXLMS((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Generates theoretical spectra for cross-linked peptides
        
        The spectra this class generates are vectors of SimplePeaks
        This class generates the same peak types as TheoreticalSpectrumGeneratorXLMS
        and the interface is very similar, but it is simpler and faster
        SimplePeak only contains an mz value and a charge. No intensity values
        or String annotations or other additional DataArrays are generated
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SimpleTSGXLMS ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SimpleTSGXLMS)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getLinearIonSpectrum(self, list spectrum , AASequence peptide ,  link_pos ,  charge ,  link_pos_2 ):
        """
        getLinearIonSpectrum(self, spectrum: List[SimplePeak] , peptide: AASequence , link_pos: int , charge: int , link_pos_2: int ) -> None
        Generates fragment ions not containing the cross-linker for one peptide
        
        B-ions are generated from the beginning of the peptide up to the first linked position,
        y-ions are generated from the second linked position up the end of the peptide
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
        The generated ion types and other additional settings are determined by the tool parameters
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param peptide: The peptide to fragment
        :param link_pos: The position of the cross-linker on the given peptide
        :param charge: The maximal charge of the ions
        :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(link_pos, int) and link_pos >= 0, 'arg link_pos wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(link_pos_2, int) and link_pos_2 >= 0, 'arg link_pos_2 wrong type'
        cdef libcpp_vector[_SimplePeak] * v0 = new libcpp_vector[_SimplePeak]()
        cdef SimplePeak item0
        for item0 in spectrum:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
        self.inst.get().getLinearIonSpectrum(deref(v0), (deref(peptide.inst.get())), (<size_t>link_pos), (<int>charge), (<size_t>link_pos_2))
        cdef libcpp_vector[_SimplePeak].iterator it_spectrum = v0.begin()
        replace_0 = []
        while it_spectrum != v0.end():
            item0 = SimplePeak.__new__(SimplePeak)
            item0.inst = shared_ptr[_SimplePeak](new _SimplePeak(deref(it_spectrum)))
            replace_0.append(item0)
            inc(it_spectrum)
        spectrum[:] = replace_0
        del v0
    
    def _getXLinkIonSpectrum_0(self, list spectrum , AASequence peptide ,  link_pos , double precursor_mass ,  mincharge ,  maxcharge ,  link_pos_2 ):
        """
        _getXLinkIonSpectrum_0(self, spectrum: List[SimplePeak] , peptide: AASequence , link_pos: int , precursor_mass: float , mincharge: int , maxcharge: int , link_pos_2: int ) -> None
        Generates fragment ions containing the cross-linker for one peptide
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
        Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
        a precursor mass for the two peptides and the linker is needed
        In the case of a loop link the precursor mass is the mass of the only peptide and the linker
        Although this function is more general, currently it is mainly used for loop-links and mono-links,
        because residues in the second, unknown peptide cannot be considered for possible neutral losses
        The generated ion types and other additional settings are determined by the tool parameters
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param peptide: The peptide to fragment
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(link_pos, int) and link_pos >= 0, 'arg link_pos wrong type'
        assert isinstance(precursor_mass, float), 'arg precursor_mass wrong type'
        assert isinstance(mincharge, int), 'arg mincharge wrong type'
        assert isinstance(maxcharge, int), 'arg maxcharge wrong type'
        assert isinstance(link_pos_2, int) and link_pos_2 >= 0, 'arg link_pos_2 wrong type'
        cdef libcpp_vector[_SimplePeak] * v0 = new libcpp_vector[_SimplePeak]()
        cdef SimplePeak item0
        for item0 in spectrum:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
    
    
        self.inst.get().getXLinkIonSpectrum(deref(v0), (deref(peptide.inst.get())), (<size_t>link_pos), (<double>precursor_mass), (<int>mincharge), (<int>maxcharge), (<size_t>link_pos_2))
        cdef libcpp_vector[_SimplePeak].iterator it_spectrum = v0.begin()
        replace_0 = []
        while it_spectrum != v0.end():
            item0 = SimplePeak.__new__(SimplePeak)
            item0.inst = shared_ptr[_SimplePeak](new _SimplePeak(deref(it_spectrum)))
            replace_0.append(item0)
            inc(it_spectrum)
        spectrum[:] = replace_0
        del v0
    
    def _getXLinkIonSpectrum_1(self, list spectrum , ProteinProteinCrossLink crosslink , bool frag_alpha ,  mincharge ,  maxcharge ):
        """
        _getXLinkIonSpectrum_1(self, spectrum: List[SimplePeak] , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None
        Generates fragment ions containing the cross-linker for a pair of peptides
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        This function generates neutral loss ions by considering both linked peptides
        Only one of the peptides, decided by @frag_alpha, is fragmented
        This simplifies the function, but it has to be called twice to get all fragments of a peptide pair
        The generated ion types and other additional settings are determined by the tool parameters
        This function is not suitable to generate fragments for mono-links or loop-links
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param crosslink: ProteinProteinCrossLink to be fragmented
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        """
        assert isinstance(spectrum, list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in spectrum), 'arg spectrum wrong type'
        assert isinstance(crosslink, ProteinProteinCrossLink), 'arg crosslink wrong type'
        assert isinstance(frag_alpha, pybool_t), 'arg frag_alpha wrong type'
        assert isinstance(mincharge, int), 'arg mincharge wrong type'
        assert isinstance(maxcharge, int), 'arg maxcharge wrong type'
        cdef libcpp_vector[_SimplePeak] * v0 = new libcpp_vector[_SimplePeak]()
        cdef SimplePeak item0
        for item0 in spectrum:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
        self.inst.get().getXLinkIonSpectrum(deref(v0), (deref(crosslink.inst.get())), (<bool>frag_alpha), (<int>mincharge), (<int>maxcharge))
        cdef libcpp_vector[_SimplePeak].iterator it_spectrum = v0.begin()
        replace_0 = []
        while it_spectrum != v0.end():
            item0 = SimplePeak.__new__(SimplePeak)
            item0.inst = shared_ptr[_SimplePeak](new _SimplePeak(deref(it_spectrum)))
            replace_0.append(item0)
            inc(it_spectrum)
        spectrum[:] = replace_0
        del v0
    
    def getXLinkIonSpectrum(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getXLinkIonSpectrum(self, spectrum: List[SimplePeak] , peptide: AASequence , link_pos: int , precursor_mass: float , mincharge: int , maxcharge: int , link_pos_2: int ) -> None
          :noindex:
        
        Generates fragment ions containing the cross-linker for one peptide
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
        Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
        a precursor mass for the two peptides and the linker is needed
        In the case of a loop link the precursor mass is the mass of the only peptide and the linker
        Although this function is more general, currently it is mainly used for loop-links and mono-links,
        because residues in the second, unknown peptide cannot be considered for possible neutral losses
        The generated ion types and other additional settings are determined by the tool parameters
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param peptide: The peptide to fragment
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        :param link_pos_2: A second position for the linker, in case it is a loop link
        
        .. rubric:: Overload:
        .. py:function:: getXLinkIonSpectrum(self, spectrum: List[SimplePeak] , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None
          :noindex:
        
        Generates fragment ions containing the cross-linker for a pair of peptides
        
        B-ions are generated from the first linked position up to the end of the peptide,
        y-ions are generated from the beginning of the peptide up to the second linked position
        This function generates neutral loss ions by considering both linked peptides
        Only one of the peptides, decided by @frag_alpha, is fragmented
        This simplifies the function, but it has to be called twice to get all fragments of a peptide pair
        The generated ion types and other additional settings are determined by the tool parameters
        This function is not suitable to generate fragments for mono-links or loop-links
        
        :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
        :param crosslink: ProteinProteinCrossLink to be fragmented
        :param link_pos: The position of the cross-linker on the given peptide
        :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
        :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
        :param mincharge: The minimal charge of the ions
        :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
    
        """
        if (len(args)==7) and (isinstance(args[0], list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in args[0])) and (isinstance(args[1], AASequence)) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], float)) and (isinstance(args[4], int)) and (isinstance(args[5], int)) and (isinstance(args[6], int) and args[6] >= 0):
            return self._getXLinkIonSpectrum_0(*args)
        elif (len(args)==5) and (isinstance(args[0], list) and all(isinstance(elemt_rec, SimplePeak) for elemt_rec in args[0])) and (isinstance(args[1], ProteinProteinCrossLink)) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], int)) and (isinstance(args[4], int)):
            return self._getXLinkIonSpectrum_1(*args)
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

cdef class XLPrecursor:
    """
    Cython implementation of _XLPrecursor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XLPrecursor.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property precursor_mass:
        def __set__(self, float precursor_mass):
        
            self.inst.get().precursor_mass = (<float>precursor_mass)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().precursor_mass
            py_result = <float>_r
            return py_result
    
    property alpha_index:
        def __set__(self,  alpha_index):
        
            self.inst.get().alpha_index = (<unsigned int>alpha_index)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().alpha_index
            py_result = <unsigned int>_r
            return py_result
    
    property beta_index:
        def __set__(self,  beta_index):
        
            self.inst.get().beta_index = (<unsigned int>beta_index)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().beta_index
            py_result = <unsigned int>_r
            return py_result
    
    def __copy__(self):
       cdef XLPrecursor rv = XLPrecursor.__new__(XLPrecursor)
       rv.inst = shared_ptr[_XLPrecursor](new _XLPrecursor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef XLPrecursor rv = XLPrecursor.__new__(XLPrecursor)
       rv.inst = shared_ptr[_XLPrecursor](new _XLPrecursor(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_XLPrecursor](new _XLPrecursor())
    
    def _init_1(self, XLPrecursor in_0 ):
        """
        _init_1(self, in_0: XLPrecursor ) -> None
        """
        assert isinstance(in_0, XLPrecursor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_XLPrecursor](new _XLPrecursor((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: XLPrecursor ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], XLPrecursor)):
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
