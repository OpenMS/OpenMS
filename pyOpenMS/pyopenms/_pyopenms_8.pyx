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
def __static_IonIdentityMolecularNetworking_annotateConsensusMap(ConsensusMap consensus_map ):
    """
    __static_IonIdentityMolecularNetworking_annotateConsensusMap(consensus_map: ConsensusMap ) -> None
        Annotate ConsensusMap for ion identity molecular networking (IIMN) workflow by GNPS.
        
        Adds meta values Constants::UserParams::IIMN_ROW_ID (unique index for each feature), Constants::UserParams::IIMN_ADDUCT_PARTNERS (related features row IDs)
        and Constants::UserParams::IIMN_ANNOTATION_NETWORK_NUMBER (all related features with different adduct states) get the same network number).
        This method requires the features annotated with the Constants::UserParams::IIMN_LINKED_GROUPS meta value.
        If at least one of the features has an annotation for Constants::UserParam::IIMN_LINKED_GROUPS, annotate ConsensusMap for IIMN.
        
        
        :param consensus_map: Input ConsensusMap without IIMN annotations.
    """
    assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'

    _annotateConsensusMap_IonIdentityMolecularNetworking((deref(consensus_map.inst.get())))

def __static_InternalCalibration_applyTransformation(list pcs , MZTrafoModel trafo ):
    """
    __static_InternalCalibration_applyTransformation(pcs: List[Precursor] , trafo: MZTrafoModel ) -> None
    """
    assert isinstance(pcs, list) and all(isinstance(elemt_rec, Precursor) for elemt_rec in pcs), 'arg pcs wrong type'
    assert isinstance(trafo, MZTrafoModel), 'arg trafo wrong type'
    cdef libcpp_vector[_Precursor] * v0 = new libcpp_vector[_Precursor]()
    cdef Precursor item0
    for item0 in pcs:
        v0.push_back(deref(item0.inst.get()))

    _applyTransformation_InternalCalibration(deref(v0), (deref(trafo.inst.get())))
    cdef libcpp_vector[_Precursor].iterator it_pcs = v0.begin()
    replace_0 = []
    while it_pcs != v0.end():
        item0 = Precursor.__new__(Precursor)
        item0.inst = shared_ptr[_Precursor](new _Precursor(deref(it_pcs)))
        replace_0.append(item0)
        inc(it_pcs)
    pcs[:] = replace_0
    del v0

def __static_InternalCalibration_applyTransformation(MSSpectrum spec , list target_mslvl , MZTrafoModel trafo ):
    """
    __static_InternalCalibration_applyTransformation(spec: MSSpectrum , target_mslvl: List[int] , trafo: MZTrafoModel ) -> None
    """
    assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    assert isinstance(target_mslvl, list) and all(isinstance(li, int) for li in target_mslvl), 'arg target_mslvl wrong type'
    assert isinstance(trafo, MZTrafoModel), 'arg trafo wrong type'

    cdef libcpp_vector[int] _v1 = target_mslvl
    cdef _IntList v1 = _IntList(_v1)

    _applyTransformation_InternalCalibration((deref(spec.inst.get())), (v1), (deref(trafo.inst.get())))

def __static_InternalCalibration_applyTransformation(MSExperiment exp , list target_mslvl , MZTrafoModel trafo ):
    """
    __static_InternalCalibration_applyTransformation(exp: MSExperiment , target_mslvl: List[int] , trafo: MZTrafoModel ) -> None
    """
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    assert isinstance(target_mslvl, list) and all(isinstance(li, int) for li in target_mslvl), 'arg target_mslvl wrong type'
    assert isinstance(trafo, MZTrafoModel), 'arg trafo wrong type'

    cdef libcpp_vector[int] _v1 = target_mslvl
    cdef _IntList v1 = _IntList(_v1)

    _applyTransformation_InternalCalibration((deref(exp.inst.get())), (v1), (deref(trafo.inst.get())))

def __static_FileHandler_computeFileHash( filename ):
    """
    __static_FileHandler_computeFileHash(filename: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'

    cdef _String _r = _computeFileHash_FileHandler(deref((convString(filename)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_FileHandler_getType( filename ):
    """
    __static_FileHandler_getType(filename: Union[bytes, str, String] ) -> int
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'

    cdef int _r = _getType_FileHandler(deref((convString(filename)).get()))
    py_result = <int>_r
    return py_result

def __static_FileHandler_getTypeByContent( filename ):
    """
    __static_FileHandler_getTypeByContent(filename: Union[bytes, str, String] ) -> int
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'

    cdef _FileType _r = _getTypeByContent_FileHandler(deref((convString(filename)).get()))
    py_result = <int>_r
    return py_result

def __static_FileHandler_getTypeByFileName( filename ):
    """
    __static_FileHandler_getTypeByFileName(filename: Union[bytes, str, String] ) -> int
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'

    cdef _FileType _r = _getTypeByFileName_FileHandler(deref((convString(filename)).get()))
    py_result = <int>_r
    return py_result

def __static_FileHandler_hasValidExtension( filename , int type_ ):
    """
    __static_FileHandler_hasValidExtension(filename: Union[bytes, str, String] , type_: int ) -> bool
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg type_ wrong type'


    cdef bool _r = _hasValidExtension_FileHandler(deref((convString(filename)).get()), (<_FileType>type_))
    py_result = <bool>_r
    return py_result

def __static_FileHandler_isSupported(int type_ ):
    """
    __static_FileHandler_isSupported(type_: int ) -> bool
    """
    assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg type_ wrong type'

    cdef bool _r = _isSupported_FileHandler((<_FileType>type_))
    py_result = <bool>_r
    return py_result

def __static_FileHandler_stripExtension( file ):
    """
    __static_FileHandler_stripExtension(file: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef _String _r = _stripExtension_FileHandler(deref((convString(file)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_FileHandler_swapExtension( filename , int new_type ):
    """
    __static_FileHandler_swapExtension(filename: Union[bytes, str, String] , new_type: int ) -> Union[bytes, str, String]
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert new_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg new_type wrong type'


    cdef _String _r = _swapExtension_FileHandler(deref((convString(filename)).get()), (<_FileType>new_type))
    py_result = convOutputString(_r)
    return py_result

def __static_IonIdentityMolecularNetworking_writeSupplementaryPairTable(ConsensusMap consensus_map ,  output_file ):
    """
    __static_IonIdentityMolecularNetworking_writeSupplementaryPairTable(consensus_map: ConsensusMap , output_file: Union[bytes, str, String] ) -> None
        Write supplementary pair table (csv file) from a ConsensusMap with edge annotations for connected features. Required for GNPS IIMN.
        
        The table contains the columns "ID 1" (row ID of first feature), "ID 2" (row ID of second feature), "EdgeType" (MS1/2 annotation),
        "Score" (the number of direct partners from both connected features) and "Annotation" (adducts and delta m/z between two connected features).
        
        
        :param consensus_map: Input ConsensusMap annotated with IonIdentityMolecularNetworking.annotateConsensusMap.
        :param output_file: Output file path for the supplementary pair table.
    """
    assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
    assert (isinstance(output_file, str) or isinstance(output_file, bytes) or isinstance(output_file, String)), 'arg output_file wrong type'


    _writeSupplementaryPairTable_IonIdentityMolecularNetworking((deref(consensus_map.inst.get())), deref((convString(output_file)).get())) 

cdef class __ActivationMethod:
    """
        Enum for activation/fragmentation methods of mass spectra
    
    - CID: Collision-induced dissociation (MS:1000133) (also CAD; parent term, but unless otherwise stated often used as synonym for trap-type CID)
    - PSD: Post-source decay
    - PD: Plasma desorption
    - SID: Surface-induced dissociation
    - BIRD: Blackbody infrared radiative dissociation
    - ECD: Electron capture dissociation (MS:1000250)
    - IMD: Infrared multiphoton dissociation
    - SORI: Sustained off-resonance irradiation
    - HCID: High-energy collision-induced dissociation
    - LCID: Low-energy collision-induced dissociation
    - PHD: Photodissociation
    - ETD: Electron transfer dissociation
    - ETciD: Electron transfer and collision-induced dissociation (MS:1003182)
    - EThcD: Electron transfer and higher-energy collision dissociation (MS:1002631)
    - PQD: Pulsed q dissociation (MS:1000599)
    - TRAP: trap-type collision-induced dissociation (MS:1002472)
    - HCD: beam-type collision-induced dissociation (MS:1000422)
    - INSOURCE: in-source collision-induced dissociation (MS:1001880)
    - LIFT: Bruker proprietary method (MS:1002000)
    """
    CID = 0
    PSD = 1
    PD = 2
    SID = 3
    BIRD = 4
    ECD = 5
    IMD = 6
    SORI = 7
    HCID = 8
    LCID = 9
    PHD = 10
    ETD = 11
    ETciD = 12
    EThcD = 13
    PQD = 14
    TRAP = 15
    HCD = 16
    INSOURCE = 17
    LIFT = 18
    SIZE_OF_ACTIVATIONMETHOD = 19

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ByteOrder:
    None
    BYTEORDER_BIGENDIAN = 0
    BYTEORDER_LITTLEENDIAN = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __FilterOperation:
    None
    GREATER_EQUAL = 0
    EQUAL = 1
    LESS_EQUAL = 2
    EXISTS = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __FilterType:
    None
    INTENSITY = 0
    QUALITY = 1
    CHARGE = 2
    SIZE = 3
    META_DATA = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class Base64:
    """
    Cython implementation of _Base64

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Base64.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Base64 rv = Base64.__new__(Base64)
       rv.inst = shared_ptr[_Base64](new _Base64(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Base64 rv = Base64.__new__(Base64)
       rv.inst = shared_ptr[_Base64](new _Base64(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Class to encode and decode Base64, it supports two precisions 32 bit (float) and 64 bit (double).
        """
        self.inst = shared_ptr[_Base64](new _Base64())
    
    def _init_1(self, Base64 in_0 ):
        """
        _init_1(self, in_0: Base64 ) -> None
        """
        assert isinstance(in_0, Base64), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Base64](new _Base64((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Class to encode and decode Base64, it supports two precisions 32 bit (float) and 64 bit (double).

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Base64 ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Base64)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def encodeIntegers(self, list in_ , int to_byte_order ,  out , bool zlib_compression ):
        """
        encodeIntegers(self, in_: List[int] , to_byte_order: int , out: String , zlib_compression: bool ) -> None
        Encodes a vector of integer point numbers to a Base64 string
        """
        assert isinstance(in_, list) and all(isinstance(elemt_rec, int) for elemt_rec in in_), 'arg in_ wrong type'
        assert to_byte_order in [0, 1], 'arg to_byte_order wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
        cdef libcpp_vector[int] v0 = in_
    
    
    
        self.inst.get().encodeIntegers(v0, (<_ByteOrder>to_byte_order), deref((<String>out).inst.get()), (<bool>zlib_compression))
        in_[:] = v0
    
    def decodeIntegers(self,  in_ , int from_byte_order , list out , bool zlib_compression ):
        """
        decodeIntegers(self, in_: Union[bytes, str, String] , from_byte_order: int , out: List[int] , zlib_compression: bool ) -> None
        Decodes a Base64 string to a vector of integer numbers
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert from_byte_order in [0, 1], 'arg from_byte_order wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, int) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
    
    
        cdef libcpp_vector[int] v2 = out
    
        self.inst.get().decodeIntegers(deref((convString(in_)).get()), (<_ByteOrder>from_byte_order), v2, (<bool>zlib_compression))
        out[:] = v2
    
    def encodeStrings(self, list in_ ,  out , bool zlib_compression ):
        """
        encodeStrings(self, in_: List[bytes] , out: String , zlib_compression: bool ) -> None
        Encodes a vector of strings to a Base64 string
        """
        assert isinstance(in_, list) and all(isinstance(i, bytes) for i in in_), 'arg in_ wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in in_:
           v0.push_back(_String(<char *>item0))
    
    
        self.inst.get().encodeStrings(deref(v0), deref((<String>out).inst.get()), (<bool>zlib_compression))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        in_[:] = replace
        del v0
    
    def decodeStrings(self,  in_ , list out , bool zlib_compression ):
        """
        decodeStrings(self, in_: Union[bytes, str, String] , out: List[bytes] , zlib_compression: bool ) -> None
        Decodes a Base64 string to a vector of (null-terminated) strings
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(out, list) and all(isinstance(i, bytes) for i in out), 'arg out wrong type'
        assert isinstance(zlib_compression, pybool_t), 'arg zlib_compression wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in out:
           v1.push_back(_String(<char *>item1))
    
        self.inst.get().decodeStrings(deref((convString(in_)).get()), deref(v1), (<bool>zlib_compression))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        out[:] = replace
        del v1
    


    def encode64(self, list in_ , int to_byte_order ,  out ,  zlib_compression ):
        """Cython signature: `void encode64(libcpp_vector[double] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)`"""
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert to_byte_order in [0, 1], 'arg to_byte_order wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, int), 'arg zlib_compression wrong type'
        cdef libcpp_vector[double] v0 = in_
    
    
    
        self.inst.get().encode(v0, (<_ByteOrder>to_byte_order), deref((<String>out).inst.get()), (<bool>zlib_compression))
    
    def decode64(self,  in_ , int from_byte_order , list out ,  zlib_compression ):
        """Cython signature: `void decode64(const String & in_, ByteOrder from_byte_order, libcpp_vector[double] & out, bool zlib_compression)`"""
        assert (isinstance(in_, str) or isinstance(in_, unicode) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert from_byte_order in [0, 1], 'arg from_byte_order wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, int), 'arg zlib_compression wrong type'
    
    
        cdef libcpp_vector[double] v2 = out
    
        self.inst.get().decode(deref((convString(in_)).get()), (<_ByteOrder>from_byte_order), v2, (<bool>zlib_compression))
        out[:] = v2

    def encode32(self, list in_ , int to_byte_order ,  out ,  zlib_compression ):
        """Cython signature: `void encode32(libcpp_vector[float] & in_, ByteOrder to_byte_order, String & out, bool zlib_compression)`"""
        assert isinstance(in_, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_), 'arg in_ wrong type'
        assert to_byte_order in [0, 1], 'arg to_byte_order wrong type'
        assert isinstance(out, String), 'arg out wrong type'
        assert isinstance(zlib_compression, int), 'arg zlib_compression wrong type'
        cdef libcpp_vector[float] v0 = in_
    
    
    
        self.inst.get().encode(v0, (<_ByteOrder>to_byte_order), deref((<String>out).inst.get()), (<bool>zlib_compression))
    
    def decode32(self,  in_ , int from_byte_order , list out ,  zlib_compression ):
        """Cython signature: `void decode32(const String & in_, ByteOrder from_byte_order, libcpp_vector[float] & out, bool zlib_compression)`"""
        assert (isinstance(in_, str) or isinstance(in_, unicode) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert from_byte_order in [0, 1], 'arg from_byte_order wrong type'
        assert isinstance(out, list) and all(isinstance(elemt_rec, float) for elemt_rec in out), 'arg out wrong type'
        assert isinstance(zlib_compression, int), 'arg zlib_compression wrong type'
    
    
        cdef libcpp_vector[float] v2 = out
    
        self.inst.get().decode(deref((convString(in_)).get()), (<_ByteOrder>from_byte_order), v2, (<bool>zlib_compression))
        out[:] = v2
    ByteOrder = __ByteOrder 

cdef class CVTermList:
    """
    Cython implementation of _CVTermList

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTermList.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVTermList rv = CVTermList.__new__(CVTermList)
       rv.inst = shared_ptr[_CVTermList](new _CVTermList(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVTermList rv = CVTermList.__new__(CVTermList)
       rv.inst = shared_ptr[_CVTermList](new _CVTermList(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVTermList](new _CVTermList())
    
    def _init_1(self, CVTermList in_0 ):
        """
        _init_1(self, in_0: CVTermList ) -> None
        """
        assert isinstance(in_0, CVTermList), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVTermList](new _CVTermList((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVTermList ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVTermList)):
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
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773863406401753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773863406401753366645
        cdef CVTerm item_1268773863406401753366645
        cdef bytes inner_key_1268773863406401753366645
        cdef list inner_values_1268773863406401753366645
        while outer_it_1268773863406401753366645 != _r.end():
           inner_key_1268773863406401753366645 = deref(outer_it_1268773863406401753366645).first.c_str()
           inner_values_1268773863406401753366645 = []
           inner_it_1268773863406401753366645 = deref(outer_it_1268773863406401753366645).second.begin()
           while inner_it_1268773863406401753366645 != deref(outer_it_1268773863406401753366645).second.end():
               item_1268773863406401753366645 = CVTerm.__new__(CVTerm)
               item_1268773863406401753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773863406401753366645)))
               inner_values_1268773863406401753366645.append(item_1268773863406401753366645)
               inc(inner_it_1268773863406401753366645)
           py_result[inner_key_1268773863406401753366645] = inner_values_1268773863406401753366645
           inc(outer_it_1268773863406401753366645)
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
        if not isinstance(other, CVTermList):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVTermList other_casted = other
        cdef CVTermList self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class DataFilter:
    """
    Cython implementation of _DataFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataFilter.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property field:
        def __set__(self, int field):
        
            self.inst.get().field = (<_FilterType>field)
        
    
        def __get__(self):
            cdef _FilterType _r = self.inst.get().field
            py_result = <int>_r
            return py_result
    
    property op:
        def __set__(self, int op):
        
            self.inst.get().op = (<_FilterOperation>op)
        
    
        def __get__(self):
            cdef _FilterOperation _r = self.inst.get().op
            py_result = <int>_r
            return py_result
    
    property value:
        def __set__(self, double value):
        
            self.inst.get().value = (<double>value)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().value
            py_result = <double>_r
            return py_result
    
    property value_string:
        def __set__(self,  value_string):
        
            self.inst.get().value_string = deref((convString(value_string)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().value_string
            py_result = convOutputString(_r)
            return py_result
    
    property meta_name:
        def __set__(self,  meta_name):
        
            self.inst.get().meta_name = deref((convString(meta_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().meta_name
            py_result = convOutputString(_r)
            return py_result
    
    property value_is_numerical:
        def __set__(self, bool value_is_numerical):
        
            self.inst.get().value_is_numerical = (<bool>value_is_numerical)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().value_is_numerical
            py_result = <bool>_r
            return py_result
    
    def __copy__(self):
       cdef DataFilter rv = DataFilter.__new__(DataFilter)
       rv.inst = shared_ptr[_DataFilter](new _DataFilter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DataFilter rv = DataFilter.__new__(DataFilter)
       rv.inst = shared_ptr[_DataFilter](new _DataFilter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DataFilter](new _DataFilter())
    
    def _init_1(self, DataFilter in_0 ):
        """
        _init_1(self, in_0: DataFilter ) -> None
        """
        assert isinstance(in_0, DataFilter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataFilter](new _DataFilter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DataFilter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DataFilter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def fromString(self,  filter_ ):
        """
        fromString(self, filter_: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filter_, str) or isinstance(filter_, bytes) or isinstance(filter_, String)), 'arg filter_ wrong type'
    
        self.inst.get().fromString(deref((convString(filter_)).get()))
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, DataFilter):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DataFilter other_casted = other
        cdef DataFilter self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class DataFilters:
    """
    Cython implementation of _DataFilters

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataFilters.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DataFilters rv = DataFilters.__new__(DataFilters)
       rv.inst = shared_ptr[_DataFilters](new _DataFilters(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DataFilters rv = DataFilters.__new__(DataFilters)
       rv.inst = shared_ptr[_DataFilters](new _DataFilters(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DataFilters](new _DataFilters())
    
    def _init_1(self, DataFilters in_0 ):
        """
        _init_1(self, in_0: DataFilters ) -> None
        """
        assert isinstance(in_0, DataFilters), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataFilters](new _DataFilters((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DataFilters ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DataFilters)):
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
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> DataFilter
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _DataFilter * _r = new _DataFilter(deref(self.inst.get())[(<size_t>in_0)])
        cdef DataFilter py_result = DataFilter.__new__(DataFilter)
        py_result.inst = shared_ptr[_DataFilter](_r)
        return py_result
    
    def add(self, DataFilter filter_ ):
        """
        add(self, filter_: DataFilter ) -> None
        """
        assert isinstance(filter_, DataFilter), 'arg filter_ wrong type'
    
        self.inst.get().add((deref(filter_.inst.get())))
    
    def remove(self,  index ):
        """
        remove(self, index: int ) -> None
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        self.inst.get().remove((<size_t>index))
    
    def replace(self,  index , DataFilter filter_ ):
        """
        replace(self, index: int , filter_: DataFilter ) -> None
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(filter_, DataFilter), 'arg filter_ wrong type'
    
    
        self.inst.get().replace((<size_t>index), (deref(filter_.inst.get())))
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def setActive(self, bool is_active ):
        """
        setActive(self, is_active: bool ) -> None
        """
        assert isinstance(is_active, pybool_t), 'arg is_active wrong type'
    
        self.inst.get().setActive((<bool>is_active))
    
    def isActive(self):
        """
        isActive(self) -> bool
        """
        cdef bool _r = self.inst.get().isActive()
        py_result = <bool>_r
        return py_result
    
    def _passes_0(self, Feature feature ):
        """
        _passes_0(self, feature: Feature ) -> bool
        """
        assert isinstance(feature, Feature), 'arg feature wrong type'
    
        cdef bool _r = self.inst.get().passes((deref(feature.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def _passes_1(self, ConsensusFeature consensus_feature ):
        """
        _passes_1(self, consensus_feature: ConsensusFeature ) -> bool
        """
        assert isinstance(consensus_feature, ConsensusFeature), 'arg consensus_feature wrong type'
    
        cdef bool _r = self.inst.get().passes((deref(consensus_feature.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def _passes_2(self, MSSpectrum spectrum ,  peak_index ):
        """
        _passes_2(self, spectrum: MSSpectrum , peak_index: int ) -> bool
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(peak_index, int) and peak_index >= 0, 'arg peak_index wrong type'
    
    
        cdef bool _r = self.inst.get().passes((deref(spectrum.inst.get())), (<size_t>peak_index))
        py_result = <bool>_r
        return py_result
    
    def passes(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: passes(self, feature: Feature ) -> bool
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: passes(self, consensus_feature: ConsensusFeature ) -> bool
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: passes(self, spectrum: MSSpectrum , peak_index: int ) -> bool
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], Feature)):
            return self._passes_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusFeature)):
            return self._passes_1(*args)
        elif (len(args)==2) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], int) and args[1] >= 0):
            return self._passes_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    FilterOperation = __FilterOperation
    FilterType = __FilterType 

cdef class FeatureFinderAlgorithmPicked:
    """
    Cython implementation of _FeatureFinderAlgorithmPicked

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderAlgorithmPicked.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureFinderAlgorithmPicked](new _FeatureFinderAlgorithmPicked())
    
    def run(self, MSExperiment input_map , FeatureMap output , Param param , FeatureMap seeds ):
        """
        run(self, input_map: MSExperiment , output: FeatureMap , param: Param , seeds: FeatureMap ) -> None
        """
        assert isinstance(input_map, MSExperiment), 'arg input_map wrong type'
        assert isinstance(output, FeatureMap), 'arg output wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
        assert isinstance(seeds, FeatureMap), 'arg seeds wrong type'
    
    
    
    
        self.inst.get().run((deref(input_map.inst.get())), (deref(output.inst.get())), (deref(param.inst.get())), (deref(seeds.inst.get())))
    
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

cdef class FileHandler:
    """
    Cython implementation of _FileHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FileHandler.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FileHandler](new _FileHandler())
    
    def loadExperiment(self,  in_0 , MSExperiment in_1 ):
        """
        loadExperiment(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        Loads a file into an MSExperiment
        
        
        :param filename: The file name of the file to load
        :param exp: The experiment to load the data into
        :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
        :param log: Progress logging mode
        :param rewrite_source_file: Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from
        :param compute_hash: If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile
        :return: true if the file could be loaded, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().loadExperiment(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def storeExperiment(self,  in_0 , MSExperiment in_1 ):
        """
        storeExperiment(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        Stores an MSExperiment to a file\n
        
        The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used
        
        
        :param filename: The name of the file to store the data in
        :param exp: The experiment to store
        :param log: Progress logging mode
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be written
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().storeExperiment(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def loadFeatures(self,  in_0 , FeatureMap in_1 ):
        """
        loadFeatures(self, in_0: Union[bytes, str, String] , in_1: FeatureMap ) -> None
        Loads a file into a FeatureMap
        
        
        :param filename: The file name of the file to load
        :param map: The FeatureMap to load the data into
        :param force_type: Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails)
        :return: true if the file could be loaded, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
    
    
        self.inst.get().loadFeatures(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        Access to the options for loading/storing
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
    computeFileHash = __static_FileHandler_computeFileHash
    getType = __static_FileHandler_getType
    getTypeByContent = __static_FileHandler_getTypeByContent
    getTypeByFileName = __static_FileHandler_getTypeByFileName
    hasValidExtension = __static_FileHandler_hasValidExtension
    isSupported = __static_FileHandler_isSupported
    stripExtension = __static_FileHandler_stripExtension
    swapExtension = __static_FileHandler_swapExtension 

cdef class InternalCalibration:
    """
    Cython implementation of _InternalCalibration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InternalCalibration.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef InternalCalibration rv = InternalCalibration.__new__(InternalCalibration)
       rv.inst = shared_ptr[_InternalCalibration](new _InternalCalibration(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InternalCalibration rv = InternalCalibration.__new__(InternalCalibration)
       rv.inst = shared_ptr[_InternalCalibration](new _InternalCalibration(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A mass recalibration method using linear/quadratic interpolation (robust/weighted) of given reference masses
        """
        self.inst = shared_ptr[_InternalCalibration](new _InternalCalibration())
    
    def _init_1(self, InternalCalibration in_0 ):
        """
        _init_1(self, in_0: InternalCalibration ) -> None
        """
        assert isinstance(in_0, InternalCalibration), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InternalCalibration](new _InternalCalibration((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A mass recalibration method using linear/quadratic interpolation (robust/weighted) of given reference masses

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InternalCalibration ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InternalCalibration)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _fillCalibrants_0(self, MSExperiment in_0 , list in_1 , double tol_ppm , bool lock_require_mono , bool lock_require_iso , CalibrationData failed_lock_masses , bool verbose ):
        """
        _fillCalibrants_0(self, in_0: MSExperiment , in_1: List[InternalCalibration_LockMass] , tol_ppm: float , lock_require_mono: bool , lock_require_iso: bool , failed_lock_masses: CalibrationData , verbose: bool ) -> int
        Extract calibrants from Raw data (mzML)\n
        
        Lock masses are searched in each spectrum and added to the internal calibrant database\n
        
        Filters can be used to exclude spurious peaks, i.e. require the calibrant peak to be monoisotopic or
        to have a +1 isotope (should not be used for very low abundant calibrants)
        If a calibrant is not found, it is added to a 'failed_lock_masses' database which is returned and not stored internally.
        The intensity of the peaks describe the reason for failed detection: 0.0 - peak not found with the given ppm tolerance;
        1.0 - peak is not monoisotopic (can only occur if 'lock_require_mono' is true)
        2.0 - peak has no +1 isotope (can only occur if 'lock_require_iso' is true)
        
        
        :param exp: Peak map containing the lock masses
        :param ref_masses: List of lock masses
        :param tol_ppm: Search window for lock masses in 'exp'
        :param lock_require_mono: Require that a lock mass is the monoisotopic peak (i.e. not an isotope peak) -- lock mass is rejected otherwise
        :param lock_require_iso: Require that a lock mass has isotope peaks to its right -- lock mass is rejected otherwise
        :param failed_lock_masses: Set of calibration masses which were not found, i.e. their expected m/z and RT positions
        :param verbose: Print information on 'lock_require_XXX' matches during search
        :return: Number of calibration masses found
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, InternalCalibration_LockMass) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(tol_ppm, float), 'arg tol_ppm wrong type'
        assert isinstance(lock_require_mono, pybool_t), 'arg lock_require_mono wrong type'
        assert isinstance(lock_require_iso, pybool_t), 'arg lock_require_iso wrong type'
        assert isinstance(failed_lock_masses, CalibrationData), 'arg failed_lock_masses wrong type'
        assert isinstance(verbose, pybool_t), 'arg verbose wrong type'
    
        cdef libcpp_vector[_InternalCalibration_LockMass] * v1 = new libcpp_vector[_InternalCalibration_LockMass]()
        cdef InternalCalibration_LockMass item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
    
    
    
    
    
        cdef size_t _r = self.inst.get().fillCalibrants((deref(in_0.inst.get())), deref(v1), (<double>tol_ppm), (<bool>lock_require_mono), (<bool>lock_require_iso), (deref(failed_lock_masses.inst.get())), (<bool>verbose))
        del v1
        py_result = <size_t>_r
        return py_result
    
    def _fillCalibrants_1(self, FeatureMap in_0 , double in_1 ):
        """
        _fillCalibrants_1(self, in_0: FeatureMap , in_1: float ) -> int
        Extract calibrants from identifications\n
        
        Extracts only the first hit from the first peptide identification of each feature
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter
        RT and m/z are taken from the features, not from the identifications (for an exception see below)!\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param fm: FeatureMap with peptide identifications
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        cdef size_t _r = self.inst.get().fillCalibrants((deref(in_0.inst.get())), (<double>in_1))
        py_result = <size_t>_r
        return py_result
    
    def _fillCalibrants_2(self, PeptideIdentificationList in_0 , double in_1 ):
        """
        _fillCalibrants_2(self, in_0: PeptideIdentificationList , in_1: float ) -> int
        Extract calibrants from identifications\n
        
        Extracts only the first hit from each peptide identification
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param pep_ids: Peptide ids (e.g. from an idXML file)
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
        """
        assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        cdef size_t _r = self.inst.get().fillCalibrants((deref(in_0.inst.get())), (<double>in_1))
        py_result = <size_t>_r
        return py_result
    
    def fillCalibrants(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: fillCalibrants(self, in_0: MSExperiment , in_1: List[InternalCalibration_LockMass] , tol_ppm: float , lock_require_mono: bool , lock_require_iso: bool , failed_lock_masses: CalibrationData , verbose: bool ) -> int
          :noindex:
        
        Extract calibrants from Raw data (mzML)\n
        
        Lock masses are searched in each spectrum and added to the internal calibrant database\n
        
        Filters can be used to exclude spurious peaks, i.e. require the calibrant peak to be monoisotopic or
        to have a +1 isotope (should not be used for very low abundant calibrants)
        If a calibrant is not found, it is added to a 'failed_lock_masses' database which is returned and not stored internally.
        The intensity of the peaks describe the reason for failed detection: 0.0 - peak not found with the given ppm tolerance;
        1.0 - peak is not monoisotopic (can only occur if 'lock_require_mono' is true)
        2.0 - peak has no +1 isotope (can only occur if 'lock_require_iso' is true)
        
        
        :param exp: Peak map containing the lock masses
        :param ref_masses: List of lock masses
        :param tol_ppm: Search window for lock masses in 'exp'
        :param lock_require_mono: Require that a lock mass is the monoisotopic peak (i.e. not an isotope peak) -- lock mass is rejected otherwise
        :param lock_require_iso: Require that a lock mass has isotope peaks to its right -- lock mass is rejected otherwise
        :param failed_lock_masses: Set of calibration masses which were not found, i.e. their expected m/z and RT positions
        :param verbose: Print information on 'lock_require_XXX' matches during search
        :return: Number of calibration masses found
        
        .. rubric:: Overload:
        .. py:function:: fillCalibrants(self, in_0: FeatureMap , in_1: float ) -> int
          :noindex:
        
        Extract calibrants from identifications\n
        
        Extracts only the first hit from the first peptide identification of each feature
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter
        RT and m/z are taken from the features, not from the identifications (for an exception see below)!\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param fm: FeatureMap with peptide identifications
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
        
        .. rubric:: Overload:
        .. py:function:: fillCalibrants(self, in_0: PeptideIdentificationList , in_1: float ) -> int
          :noindex:
        
        Extract calibrants from identifications\n
        
        Extracts only the first hit from each peptide identification
        Hits are sorted beforehand
        Ambiguities should be resolved before, e.g. using IDFilter\n
        
        Unassigned peptide identifications are also taken into account!
        RT and m/z are naturally taken from the IDs, since to feature is assigned
        If you do not want these IDs, remove them from the feature map before calling this function\n
        
        A filtering step is done in the m/z dimension using 'tol_ppm'
        Since precursor masses could be annotated wrongly (e.g. isotope peak instead of mono),
        larger outliers are removed before accepting an ID as calibrant
        
        
        :param pep_ids: Peptide ids (e.g. from an idXML file)
        :param tol_ppm: Only accept ID's whose theoretical mass deviates at most this much from annotated
        :return: Number of calibration masses found
    
        """
        if (len(args)==7) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, InternalCalibration_LockMass) for elemt_rec in args[1])) and (isinstance(args[2], float)) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], CalibrationData)) and (isinstance(args[6], pybool_t)):
            return self._fillCalibrants_0(*args)
        elif (len(args)==2) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], float)):
            return self._fillCalibrants_1(*args)
        elif (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], float)):
            return self._fillCalibrants_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCalibrationPoints(self):
        """
        getCalibrationPoints(self) -> CalibrationData
        Get container of calibration points\n
        
        Filled using fillCalibrants() methods
        
        
        :return: Container of calibration points
        """
        cdef _CalibrationData * _r = new _CalibrationData(self.inst.get().getCalibrationPoints())
        cdef CalibrationData py_result = CalibrationData.__new__(CalibrationData)
        py_result.inst = shared_ptr[_CalibrationData](_r)
        return py_result
    
    def calibrate(self, MSExperiment in_0 , list in_1 , int in_2 , double rt_chunk , bool use_RANSAC , double post_ppm_median , double post_ppm_MAD ,  file_models ,  file_models_plot ,  file_residuals ,  file_residuals_plot ,  rscript_executable ):
        """
        calibrate(self, in_0: MSExperiment , in_1: List[int] , in_2: int , rt_chunk: float , use_RANSAC: bool , post_ppm_median: float , post_ppm_MAD: float , file_models: Union[bytes, str, String] , file_models_plot: Union[bytes, str, String] , file_residuals: Union[bytes, str, String] , file_residuals_plot: Union[bytes, str, String] , rscript_executable: Union[bytes, str, String] ) -> bool
        Apply calibration to data\n
        
        For each spectrum, a calibration model will be computed and applied.
        Make sure to call fillCalibrants() before, so a model can be created.\n
        
        The MSExperiment will be sorted by RT and m/z if unsorted.
        
        
        :param exp: MSExperiment holding the Raw data to calibrate
        :param target_mslvl: MS-levels where calibration should be applied to
        :param model_type: Linear or quadratic model; select based on your instrument
        :param rt_chunk: RT-window size (one-sided) of calibration points to collect around each spectrum. Set to negative values, to build one global model instead.
        :param use_RANSAC: Remove outliers before fitting a model?!
        :param post_ppm_median: The median ppm error of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
        :param post_ppm_MAD: The median absolute deviation of the calibrants must be at least this good after calibration; otherwise this method returns false(fail)
        :param file_models: Output CSV filename, where model parameters are written to (pass empty string to skip)
        :param file_models_plot: Output PNG image model parameters (pass empty string to skip)
        :param file_residuals: Output CSV filename, where ppm errors of calibrants before and after model fitting parameters are written to (pass empty string to skip)
        :param file_residuals_plot: Output PNG image of the ppm errors of calibrants (pass empty string to skip)
        :param rscript_executable: Full path to the Rscript executable
        :return: true upon successful calibration
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, int) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert in_2 in [0, 1, 2, 3, 4], 'arg in_2 wrong type'
        assert isinstance(rt_chunk, float), 'arg rt_chunk wrong type'
        assert isinstance(use_RANSAC, pybool_t), 'arg use_RANSAC wrong type'
        assert isinstance(post_ppm_median, float), 'arg post_ppm_median wrong type'
        assert isinstance(post_ppm_MAD, float), 'arg post_ppm_MAD wrong type'
        assert (isinstance(file_models, str) or isinstance(file_models, bytes) or isinstance(file_models, String)), 'arg file_models wrong type'
        assert (isinstance(file_models_plot, str) or isinstance(file_models_plot, bytes) or isinstance(file_models_plot, String)), 'arg file_models_plot wrong type'
        assert (isinstance(file_residuals, str) or isinstance(file_residuals, bytes) or isinstance(file_residuals, String)), 'arg file_residuals wrong type'
        assert (isinstance(file_residuals_plot, str) or isinstance(file_residuals_plot, bytes) or isinstance(file_residuals_plot, String)), 'arg file_residuals_plot wrong type'
        assert (isinstance(rscript_executable, str) or isinstance(rscript_executable, bytes) or isinstance(rscript_executable, String)), 'arg rscript_executable wrong type'
    
        cdef libcpp_vector[int] v1 = in_1
    
    
    
    
    
    
    
    
    
    
        cdef bool _r = self.inst.get().calibrate((deref(in_0.inst.get())), v1, (<_MZTrafoModel_MODELTYPE>in_2), (<double>rt_chunk), (<bool>use_RANSAC), (<double>post_ppm_median), (<double>post_ppm_MAD), deref((convString(file_models)).get()), deref((convString(file_models_plot)).get()), deref((convString(file_residuals)).get()), deref((convString(file_residuals_plot)).get()), deref((convString(rscript_executable)).get()))
        
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
    applyTransformation = __static_InternalCalibration_applyTransformation
    applyTransformation = __static_InternalCalibration_applyTransformation
    applyTransformation = __static_InternalCalibration_applyTransformation 

cdef class InternalCalibration_LockMass:
    """
    Cython implementation of _InternalCalibration_LockMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InternalCalibration_LockMass.html>`_
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
    
    property ms_level:
        def __set__(self,  ms_level):
        
            self.inst.get().ms_level = (<unsigned int>ms_level)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().ms_level
            py_result = <unsigned int>_r
            return py_result
    
    property charge:
        def __set__(self,  charge):
        
            self.inst.get().charge = (<int>charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().charge
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef InternalCalibration_LockMass rv = InternalCalibration_LockMass.__new__(InternalCalibration_LockMass)
       rv.inst = shared_ptr[_InternalCalibration_LockMass](new _InternalCalibration_LockMass(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InternalCalibration_LockMass rv = InternalCalibration_LockMass.__new__(InternalCalibration_LockMass)
       rv.inst = shared_ptr[_InternalCalibration_LockMass](new _InternalCalibration_LockMass(deref(self.inst.get())))
       return rv
    
    def _init_0(self, double mz_ ,  lvl_ ,  charge_ ):
        """
        _init_0(self, mz_: float , lvl_: int , charge_: int ) -> None
        """
        assert isinstance(mz_, float), 'arg mz_ wrong type'
        assert isinstance(lvl_, int), 'arg lvl_ wrong type'
        assert isinstance(charge_, int), 'arg charge_ wrong type'
    
    
    
        self.inst = shared_ptr[_InternalCalibration_LockMass](new _InternalCalibration_LockMass((<double>mz_), (<int>lvl_), (<int>charge_)))
    
    def _init_1(self, InternalCalibration_LockMass in_0 ):
        """
        _init_1(self, in_0: InternalCalibration_LockMass ) -> None
        """
        assert isinstance(in_0, InternalCalibration_LockMass), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InternalCalibration_LockMass](new _InternalCalibration_LockMass((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, mz_: float , lvl_: int , charge_: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InternalCalibration_LockMass ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], int)) and (isinstance(args[2], int)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InternalCalibration_LockMass)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class IonIdentityMolecularNetworking:
    """
    Cython implementation of _IonIdentityMolecularNetworking

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonIdentityMolecularNetworking.html>`_

    Includes the necessary functions to generate filed required for GNPS ion identity molecular networking (IIMN).
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_IonIdentityMolecularNetworking](new _IonIdentityMolecularNetworking())
    annotateConsensusMap = __static_IonIdentityMolecularNetworking_annotateConsensusMap
    writeSupplementaryPairTable = __static_IonIdentityMolecularNetworking_writeSupplementaryPairTable 

cdef class MRMAssay:
    """
    Cython implementation of _MRMAssay

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMAssay.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMAssay rv = MRMAssay.__new__(MRMAssay)
       rv.inst = shared_ptr[_MRMAssay](new _MRMAssay(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMAssay rv = MRMAssay.__new__(MRMAssay)
       rv.inst = shared_ptr[_MRMAssay](new _MRMAssay(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMAssay](new _MRMAssay())
    
    def _init_1(self, MRMAssay in_0 ):
        """
        _init_1(self, in_0: MRMAssay ) -> None
        """
        assert isinstance(in_0, MRMAssay), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMAssay](new _MRMAssay((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMAssay ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMAssay)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def reannotateTransitions(self, TargetedExperiment exp , double precursor_mz_threshold , double product_mz_threshold , list fragment_types , list fragment_charges , bool enable_specific_losses , bool enable_unspecific_losses ,  round_decPow ):
        """
        reannotateTransitions(self, exp: TargetedExperiment , precursor_mz_threshold: float , product_mz_threshold: float , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None
        Annotates and filters transitions in a TargetedExperiment
        
        
        :param exp: The input, unfiltered transitions
        :param precursor_mz_threshold: The precursor m/z threshold in Th for annotation
        :param product_mz_threshold: The product m/z threshold in Th for annotation
        :param fragment_types: The fragment types to consider for annotation
        :param fragment_charges: The fragment charges to consider for annotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param round_decPow: Round product m/z values to decimal power (default: -4)
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(precursor_mz_threshold, float), 'arg precursor_mz_threshold wrong type'
        assert isinstance(product_mz_threshold, float), 'arg product_mz_threshold wrong type'
        assert isinstance(fragment_types, list) and all(isinstance(i, bytes) for i in fragment_types), 'arg fragment_types wrong type'
        assert isinstance(fragment_charges, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in fragment_charges), 'arg fragment_charges wrong type'
        assert isinstance(enable_specific_losses, pybool_t), 'arg enable_specific_losses wrong type'
        assert isinstance(enable_unspecific_losses, pybool_t), 'arg enable_unspecific_losses wrong type'
        assert isinstance(round_decPow, int), 'arg round_decPow wrong type'
    
    
    
        cdef libcpp_vector[_String] * v3 = new libcpp_vector[_String]()
        cdef bytes item3
        for item3 in fragment_types:
           v3.push_back(_String(<char *>item3))
        cdef libcpp_vector[size_t] v4 = fragment_charges
    
    
    
        self.inst.get().reannotateTransitions((deref(exp.inst.get())), (<double>precursor_mz_threshold), (<double>product_mz_threshold), deref(v3), v4, (<bool>enable_specific_losses), (<bool>enable_unspecific_losses), (<int>round_decPow))
        
        del v3
    
    def restrictTransitions(self, TargetedExperiment exp , double lower_mz_limit , double upper_mz_limit , list swathes ):
        """
        restrictTransitions(self, exp: TargetedExperiment , lower_mz_limit: float , upper_mz_limit: float , swathes: List[List[float, float]] ) -> None
        Restrict and filter transitions in a TargetedExperiment
        
        
        :param exp: The input, unfiltered transitions
        :param lower_mz_limit: The lower product m/z limit in Th
        :param upper_mz_limit: The upper product m/z limit in Th
        :param swathes: The swath window settings (to exclude fragment ions falling into the precursor isolation window)
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(lower_mz_limit, float), 'arg lower_mz_limit wrong type'
        assert isinstance(upper_mz_limit, float), 'arg upper_mz_limit wrong type'
        assert isinstance(swathes, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in swathes), 'arg swathes wrong type'
    
    
    
        cdef libcpp_vector[libcpp_pair[double,double]] v3 = swathes
        self.inst.get().restrictTransitions((deref(exp.inst.get())), (<double>lower_mz_limit), (<double>upper_mz_limit), v3)
        
    
    def detectingTransitions(self, TargetedExperiment exp ,  min_transitions ,  max_transitions ):
        """
        detectingTransitions(self, exp: TargetedExperiment , min_transitions: int , max_transitions: int ) -> None
        Select detecting fragment ions
        
        
        :param exp: The input, unfiltered transitions
        :param min_transitions: The minimum number of transitions required per assay
        :param max_transitions: The maximum number of transitions required per assay
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(min_transitions, int), 'arg min_transitions wrong type'
        assert isinstance(max_transitions, int), 'arg max_transitions wrong type'
    
    
    
        self.inst.get().detectingTransitions((deref(exp.inst.get())), (<int>min_transitions), (<int>max_transitions))
    
    def filterMinMaxTransitionsCompound(self, TargetedExperiment exp ,  min_transitions ,  max_transitions ):
        """
        filterMinMaxTransitionsCompound(self, exp: TargetedExperiment , min_transitions: int , max_transitions: int ) -> None
        Filters target and decoy transitions by intensity, only keeping the top N transitions
        
        
        :param exp: The transition list which will be filtered
        :param min_transitions: The minimum number of transitions required per assay (targets only)
        :param max_transitions: The maximum number of transitions allowed per assay
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(min_transitions, int), 'arg min_transitions wrong type'
        assert isinstance(max_transitions, int), 'arg max_transitions wrong type'
    
    
    
        self.inst.get().filterMinMaxTransitionsCompound((deref(exp.inst.get())), (<int>min_transitions), (<int>max_transitions))
    
    def filterUnreferencedDecoysCompound(self, TargetedExperiment exp ):
        """
        filterUnreferencedDecoysCompound(self, exp: TargetedExperiment ) -> None
        Filters decoy transitions, which do not have respective target transition
        based on the transitionID.
        
        References between targets and decoys will be constructed based on the transitionsID
        and the "_decoy_" string. For example:
        
        target: 84_CompoundName_[M+H]+_88_22
        decoy: 84_CompoundName_decoy_[M+H]+_88_22
        
        
        :param exp: The transition list which will be filtered
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
    
        self.inst.get().filterUnreferencedDecoysCompound((deref(exp.inst.get())))
    
    def uisTransitions(self, TargetedExperiment exp , list fragment_types , list fragment_charges , bool enable_specific_losses , bool enable_unspecific_losses , bool enable_ms2_precursors , double mz_threshold , list swathes ,  round_decPow ,  max_num_alternative_localizations ,  shuffle_seed ):
        """
        uisTransitions(self, exp: TargetedExperiment , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , enable_ms2_precursors: bool , mz_threshold: float , swathes: List[List[float, float]] , round_decPow: int , max_num_alternative_localizations: int , shuffle_seed: int ) -> None
        Annotate UIS / site-specific transitions
        
        Performs the following actions:
        
        - Step 1: For each peptide, compute all theoretical alternative peptidoforms; see transitions generateTargetInSilicoMap_()
        - Step 2: Generate target identification transitions; see generateTargetAssays_()
        
        - Step 3a: Generate decoy sequences that share peptidoform properties with targets; see generateDecoySequences_()
        - Step 3b: Generate decoy in silico peptide map containing theoretical transition; see generateDecoyInSilicoMap_()
        - Step 4: Generate decoy identification transitions; see generateDecoyAssays_()
        
        The IPF algorithm uses the concept of "identification transitions" that
        are used to discriminate different peptidoforms, these are generated in
        this function.  In brief, the algorithm takes the existing set of
        peptides and transitions and then appends these "identification
        transitions" for targets and decoys. The novel transitions are set to be
        non-detecting and non-quantifying and are annotated with the set of
        peptidoforms to which they map.
        
        
        :param exp: The input, unfiltered transitions
        :param fragment_types: The fragment types to consider for annotation
        :param fragment_charges: The fragment charges to consider for annotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param enable_ms2_precursors: Whether MS2 precursors should be considered
        :param mz_threshold: The product m/z threshold in Th for annotation
        :param swathes: The swath window settings (to exclude fragment ions falling
        :param round_decPow: Round product m/z values to decimal power (default: -4)
        :param max_num_alternative_localizations: Maximum number of allowed peptide sequence permutations
        :param shuffle_seed: Set seed for shuffle (-1: select seed based on time)
        :param disable_decoy_transitions: Whether to disable generation of decoy UIS transitions
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(fragment_types, list) and all(isinstance(i, bytes) for i in fragment_types), 'arg fragment_types wrong type'
        assert isinstance(fragment_charges, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in fragment_charges), 'arg fragment_charges wrong type'
        assert isinstance(enable_specific_losses, pybool_t), 'arg enable_specific_losses wrong type'
        assert isinstance(enable_unspecific_losses, pybool_t), 'arg enable_unspecific_losses wrong type'
        assert isinstance(enable_ms2_precursors, pybool_t), 'arg enable_ms2_precursors wrong type'
        assert isinstance(mz_threshold, float), 'arg mz_threshold wrong type'
        assert isinstance(swathes, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in swathes), 'arg swathes wrong type'
        assert isinstance(round_decPow, int), 'arg round_decPow wrong type'
        assert isinstance(max_num_alternative_localizations, int) and max_num_alternative_localizations >= 0, 'arg max_num_alternative_localizations wrong type'
        assert isinstance(shuffle_seed, int), 'arg shuffle_seed wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in fragment_types:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[size_t] v2 = fragment_charges
    
    
    
    
        cdef libcpp_vector[libcpp_pair[double,double]] v7 = swathes
    
    
    
        self.inst.get().uisTransitions((deref(exp.inst.get())), deref(v1), v2, (<bool>enable_specific_losses), (<bool>enable_unspecific_losses), (<bool>enable_ms2_precursors), (<double>mz_threshold), v7, (<int>round_decPow), (<size_t>max_num_alternative_localizations), (<int>shuffle_seed))
        
        
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

cdef class MRMFeaturePickerFile:
    """
    Cython implementation of _MRMFeaturePickerFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeaturePickerFile.html>`_

    _MRMFeaturePickerFile_ loads components and components groups parameters from a .csv file
    
    The structures defined in [MRMFeaturePicker](@ref MRMFeaturePicker) are used
    
    It is required that columns `component_name` and `component_group_name` are present.
    Lines whose `component_name`'s or `component_group_name`'s value is an empty string, will be skipped.
    The class supports the absence of information within other columns.
    
    A reduced example of the expected format (fewer columns are shown here):
    > component_name,component_group_name,TransitionGroupPicker:stop_after_feature,TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length
    > arg-L.arg-L_1.Heavy,arg-L,2,15
    > arg-L.arg-L_1.Light,arg-L,2,17
    > orn.orn_1.Heavy,orn,3,21
    > orn.orn_1.Light,orn,3,13
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMFeaturePickerFile rv = MRMFeaturePickerFile.__new__(MRMFeaturePickerFile)
       rv.inst = shared_ptr[_MRMFeaturePickerFile](new _MRMFeaturePickerFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeaturePickerFile rv = MRMFeaturePickerFile.__new__(MRMFeaturePickerFile)
       rv.inst = shared_ptr[_MRMFeaturePickerFile](new _MRMFeaturePickerFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeaturePickerFile](new _MRMFeaturePickerFile())
    
    def _init_1(self, MRMFeaturePickerFile in_0 ):
        """
        _init_1(self, in_0: MRMFeaturePickerFile ) -> None
        """
        assert isinstance(in_0, MRMFeaturePickerFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeaturePickerFile](new _MRMFeaturePickerFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeaturePickerFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeaturePickerFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , list cp_list , list cgp_list ):
        """
        load(self, filename: Union[bytes, str, String] , cp_list: List[MRMFP_ComponentParams] , cgp_list: List[MRMFP_ComponentGroupParams] ) -> None
        Loads the file's data and saves it into vectors of `ComponentParams` and `ComponentGroupParams`
        
        The file is expected to contain at least two columns: `component_name` and `component_group_name`. Otherwise,
        an exception is thrown
        
        If a component group (identified by its name) is found multiple times, only the first one is saved
        
        
        :param filename: Path to the .csv input file
        :param cp_list: Component params are saved in this list
        :param cgp_list: Component Group params are saved in this list
        :raises:
          Exception: MissingInformation If the required columns are not found
        :raises:
          Exception: FileNotFound If input file is not found
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(cp_list, list) and all(isinstance(elemt_rec, MRMFP_ComponentParams) for elemt_rec in cp_list), 'arg cp_list wrong type'
        assert isinstance(cgp_list, list) and all(isinstance(elemt_rec, MRMFP_ComponentGroupParams) for elemt_rec in cgp_list), 'arg cgp_list wrong type'
    
        cdef libcpp_vector[_MRMFP_ComponentParams] * v1 = new libcpp_vector[_MRMFP_ComponentParams]()
        cdef MRMFP_ComponentParams item1
        for item1 in cp_list:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[_MRMFP_ComponentGroupParams] * v2 = new libcpp_vector[_MRMFP_ComponentGroupParams]()
        cdef MRMFP_ComponentGroupParams item2
        for item2 in cgp_list:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), deref(v2))
        cdef libcpp_vector[_MRMFP_ComponentGroupParams].iterator it_cgp_list = v2.begin()
        replace_0 = []
        while it_cgp_list != v2.end():
            item2 = MRMFP_ComponentGroupParams.__new__(MRMFP_ComponentGroupParams)
            item2.inst = shared_ptr[_MRMFP_ComponentGroupParams](new _MRMFP_ComponentGroupParams(deref(it_cgp_list)))
            replace_0.append(item2)
            inc(it_cgp_list)
        cgp_list[:] = replace_0
        del v2
        cdef libcpp_vector[_MRMFP_ComponentParams].iterator it_cp_list = v1.begin()
        replace_0 = []
        while it_cp_list != v1.end():
            item1 = MRMFP_ComponentParams.__new__(MRMFP_ComponentParams)
            item1.inst = shared_ptr[_MRMFP_ComponentParams](new _MRMFP_ComponentParams(deref(it_cp_list)))
            replace_0.append(item1)
            inc(it_cp_list)
        cp_list[:] = replace_0
        del v1 

cdef class MSDataAggregatingConsumer:
    """
    Cython implementation of _MSDataAggregatingConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataAggregatingConsumer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSDataAggregatingConsumer rv = MSDataAggregatingConsumer.__new__(MSDataAggregatingConsumer)
       rv.inst = shared_ptr[_MSDataAggregatingConsumer](new _MSDataAggregatingConsumer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSDataAggregatingConsumer rv = MSDataAggregatingConsumer.__new__(MSDataAggregatingConsumer)
       rv.inst = shared_ptr[_MSDataAggregatingConsumer](new _MSDataAggregatingConsumer(deref(self.inst.get())))
       return rv
    
    def __init__(self, MSDataAggregatingConsumer in_0 ):
        """
        __init__(self, in_0: MSDataAggregatingConsumer ) -> None
        """
        assert isinstance(in_0, MSDataAggregatingConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSDataAggregatingConsumer](new _MSDataAggregatingConsumer((deref(in_0.inst.get()))))
    
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
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get()))) 

cdef class MultiplexDeltaMasses:
    """
    Cython implementation of _MultiplexDeltaMasses

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMasses.html>`_

    Data structure for mass shift pattern
    
    Groups of labelled peptides appear with characteristic mass shifts
    
    For example, for an Arg6 labeled SILAC peptide pair we expect to see
    mass shifts of 0 and 6 Da. Or as second example, for a
    peptide pair of a dimethyl labelled sample with a single lysine
    we will see mass shifts of 56 Da and 64 Da.
    28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K)
    for light and heavy partners respectively
    
    The data structure stores the mass shifts and corresponding labels
    for a group of matching peptide features
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MultiplexDeltaMasses rv = MultiplexDeltaMasses.__new__(MultiplexDeltaMasses)
       rv.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MultiplexDeltaMasses rv = MultiplexDeltaMasses.__new__(MultiplexDeltaMasses)
       rv.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses())
    
    def _init_1(self, MultiplexDeltaMasses in_0 ):
        """
        _init_1(self, in_0: MultiplexDeltaMasses ) -> None
        """
        assert isinstance(in_0, MultiplexDeltaMasses), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses((deref(in_0.inst.get()))))
    
    def _init_2(self, list dm ):
        """
        _init_2(self, dm: List[MultiplexDeltaMasses_DeltaMass] ) -> None
        """
        assert isinstance(dm, list) and all(isinstance(elemt_rec, MultiplexDeltaMasses_DeltaMass) for elemt_rec in dm), 'arg dm wrong type'
        cdef libcpp_vector[_MultiplexDeltaMasses_DeltaMass] * v0 = new libcpp_vector[_MultiplexDeltaMasses_DeltaMass]()
        cdef MultiplexDeltaMasses_DeltaMass item0
        for item0 in dm:
            v0.push_back(deref(item0.inst.get()))
        self.inst = shared_ptr[_MultiplexDeltaMasses](new _MultiplexDeltaMasses(deref(v0)))
        cdef libcpp_vector[_MultiplexDeltaMasses_DeltaMass].iterator it_dm = v0.begin()
        replace_0 = []
        while it_dm != v0.end():
            item0 = MultiplexDeltaMasses_DeltaMass.__new__(MultiplexDeltaMasses_DeltaMass)
            item0.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass(deref(it_dm)))
            replace_0.append(item0)
            inc(it_dm)
        dm[:] = replace_0
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MultiplexDeltaMasses ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, dm: List[MultiplexDeltaMasses_DeltaMass] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MultiplexDeltaMasses)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MultiplexDeltaMasses_DeltaMass) for elemt_rec in args[0])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getDeltaMasses(self):
        """
        getDeltaMasses(self) -> List[MultiplexDeltaMasses_DeltaMass]
        """
        _r = self.inst.get().getDeltaMasses()
        py_result = []
        cdef libcpp_vector[_MultiplexDeltaMasses_DeltaMass].iterator it__r = _r.begin()
        cdef MultiplexDeltaMasses_DeltaMass item_py_result
        while it__r != _r.end():
           item_py_result = MultiplexDeltaMasses_DeltaMass.__new__(MultiplexDeltaMasses_DeltaMass)
           item_py_result.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result 

cdef class MultiplexDeltaMasses_DeltaMass:
    """
    Cython implementation of _MultiplexDeltaMasses_DeltaMass

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexDeltaMasses_DeltaMass.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property delta_mass:
        def __set__(self, double delta_mass):
        
            self.inst.get().delta_mass = (<double>delta_mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().delta_mass
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef MultiplexDeltaMasses_DeltaMass rv = MultiplexDeltaMasses_DeltaMass.__new__(MultiplexDeltaMasses_DeltaMass)
       rv.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MultiplexDeltaMasses_DeltaMass rv = MultiplexDeltaMasses_DeltaMass.__new__(MultiplexDeltaMasses_DeltaMass)
       rv.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass(deref(self.inst.get())))
       return rv
    
    def _init_0(self, double dm ,  l ):
        """
        _init_0(self, dm: float , l: Union[bytes, str, String] ) -> None
        """
        assert isinstance(dm, float), 'arg dm wrong type'
        assert (isinstance(l, str) or isinstance(l, bytes) or isinstance(l, String)), 'arg l wrong type'
    
    
        self.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass((<double>dm), deref((convString(l)).get())))
    
    def _init_1(self, MultiplexDeltaMasses_DeltaMass in_0 ):
        """
        _init_1(self, in_0: MultiplexDeltaMasses_DeltaMass ) -> None
        """
        assert isinstance(in_0, MultiplexDeltaMasses_DeltaMass), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MultiplexDeltaMasses_DeltaMass](new _MultiplexDeltaMasses_DeltaMass((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, dm: float , l: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MultiplexDeltaMasses_DeltaMass ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], float)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MultiplexDeltaMasses_DeltaMass)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class Precursor:
    """
    Cython implementation of _Precursor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Precursor.html>`_
      -- Inherits from ['Peak1D', 'CVTermList']

    Precursor meta information
    
    This class contains precursor information:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Precursor rv = Precursor.__new__(Precursor)
       rv.inst = shared_ptr[_Precursor](new _Precursor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Precursor rv = Precursor.__new__(Precursor)
       rv.inst = shared_ptr[_Precursor](new _Precursor(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Precursor](new _Precursor())
    
    def _init_1(self, Precursor in_0 ):
        """
        _init_1(self, in_0: Precursor ) -> None
        """
        assert isinstance(in_0, Precursor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Precursor](new _Precursor((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Precursor ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Precursor)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getActivationMethods(self):
        """
        getActivationMethods(self) -> Set[int]
        Returns the activation methods
        """
        _r = self.inst.get().getActivationMethods()
        py_result = set()
        cdef libcpp_set[_ActivationMethod].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<int>deref(it__r))
           inc(it__r)
        return py_result
    
    def getActivationMethodsAsString(self):
        """
        getActivationMethodsAsString(self) -> List[bytes]
        Returns the full names (e.g., "Collision-induced dissociation") of the activation methods set on this instance
        """
        _r = self.inst.get().getActivationMethodsAsString()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getActivationMethodsAsShortString(self):
        """
        getActivationMethodsAsShortString(self) -> List[bytes]
        Returns the abbreviations (e.g., "CID") of the activation methods set on this instance
        """
        _r = self.inst.get().getActivationMethodsAsShortString()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setActivationMethods(self, set activation_methods ):
        """
        setActivationMethods(self, activation_methods: Set[int] ) -> None
        Sets the activation methods
        """
        assert isinstance(activation_methods, set) and all(li in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19] for li in activation_methods), 'arg activation_methods wrong type'
        cdef libcpp_set[_ActivationMethod] * v0 = new libcpp_set[_ActivationMethod]()
        cdef int item0
        for item0 in activation_methods:
           v0.insert(<_ActivationMethod> item0)
        self.inst.get().setActivationMethods(deref(v0))
        del v0
    
    @staticmethod
    def getAllNamesOfActivationMethods():
        """
        getAllNamesOfActivationMethods() -> List[bytes]
        Returns the full names (e.g., "Collision-induced dissociation") of ALL possible activation methods, not just those set on this instance
        """
        _r = _Precursor.getAllNamesOfActivationMethods()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    @staticmethod
    def getAllShortNamesOfActivationMethods():
        """
        getAllShortNamesOfActivationMethods() -> List[bytes]
        Returns the abbreviations (e.g., "CID") of ALL possible activation methods, not just those set on this instance
        """
        _r = _Precursor.getAllShortNamesOfActivationMethods()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getActivationEnergy(self):
        """
        getActivationEnergy(self) -> float
        Returns the activation energy (in electronvolt)
        """
        cdef double _r = self.inst.get().getActivationEnergy()
        py_result = <double>_r
        return py_result
    
    def setActivationEnergy(self, double activation_energy ):
        """
        setActivationEnergy(self, activation_energy: float ) -> None
        Sets the activation energy (in electronvolt)
        """
        assert isinstance(activation_energy, float), 'arg activation_energy wrong type'
    
        self.inst.get().setActivationEnergy((<double>activation_energy))
    
    def getIsolationWindowLowerOffset(self):
        """
        getIsolationWindowLowerOffset(self) -> float
        Returns the lower offset from the target m/z
        """
        cdef double _r = self.inst.get().getIsolationWindowLowerOffset()
        py_result = <double>_r
        return py_result
    
    def setIsolationWindowLowerOffset(self, double bound ):
        """
        setIsolationWindowLowerOffset(self, bound: float ) -> None
        Sets the lower offset from the target m/z
        """
        assert isinstance(bound, float), 'arg bound wrong type'
    
        self.inst.get().setIsolationWindowLowerOffset((<double>bound))
    
    def getDriftTime(self):
        """
        getDriftTime(self) -> float
        Returns the ion mobility drift time in milliseconds (-1 means it is not set)
        """
        cdef double _r = self.inst.get().getDriftTime()
        py_result = <double>_r
        return py_result
    
    def setDriftTime(self, double drift_time ):
        """
        setDriftTime(self, drift_time: float ) -> None
        Sets the ion mobility drift time in milliseconds
        """
        assert isinstance(drift_time, float), 'arg drift_time wrong type'
    
        self.inst.get().setDriftTime((<double>drift_time))
    
    def getIsolationWindowUpperOffset(self):
        """
        getIsolationWindowUpperOffset(self) -> float
        Returns the upper offset from the target m/z
        """
        cdef double _r = self.inst.get().getIsolationWindowUpperOffset()
        py_result = <double>_r
        return py_result
    
    def setIsolationWindowUpperOffset(self, double bound ):
        """
        setIsolationWindowUpperOffset(self, bound: float ) -> None
        Sets the upper offset from the target m/z
        """
        assert isinstance(bound, float), 'arg bound wrong type'
    
        self.inst.get().setIsolationWindowUpperOffset((<double>bound))
    
    def getDriftTimeWindowLowerOffset(self):
        """
        getDriftTimeWindowLowerOffset(self) -> float
        Returns the lower offset from the target ion mobility in milliseconds
        """
        cdef double _r = self.inst.get().getDriftTimeWindowLowerOffset()
        py_result = <double>_r
        return py_result
    
    def setDriftTimeWindowLowerOffset(self, double drift_time ):
        """
        setDriftTimeWindowLowerOffset(self, drift_time: float ) -> None
        Sets the lower offset from the target ion mobility
        """
        assert isinstance(drift_time, float), 'arg drift_time wrong type'
    
        self.inst.get().setDriftTimeWindowLowerOffset((<double>drift_time))
    
    def getDriftTimeWindowUpperOffset(self):
        """
        getDriftTimeWindowUpperOffset(self) -> float
        Returns the upper offset from the target ion mobility in milliseconds
        """
        cdef double _r = self.inst.get().getDriftTimeWindowUpperOffset()
        py_result = <double>_r
        return py_result
    
    def setDriftTimeWindowUpperOffset(self, double drift_time ):
        """
        setDriftTimeWindowUpperOffset(self, drift_time: float ) -> None
        Sets the upper offset from the target ion mobility
        """
        assert isinstance(drift_time, float), 'arg drift_time wrong type'
    
        self.inst.get().setDriftTimeWindowUpperOffset((<double>drift_time))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns the charge
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def setCharge(self,  charge ):
        """
        setCharge(self, charge: int ) -> None
        Sets the charge
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        self.inst.get().setCharge((<int>charge))
    
    def getPossibleChargeStates(self):
        """
        getPossibleChargeStates(self) -> List[int]
        Returns the possible charge states
        """
        _r = self.inst.get().getPossibleChargeStates()
        cdef list py_result = _r
        return py_result
    
    def setPossibleChargeStates(self, list possible_charge_states ):
        """
        setPossibleChargeStates(self, possible_charge_states: List[int] ) -> None
        Sets the possible charge states
        """
        assert isinstance(possible_charge_states, list) and all(isinstance(elemt_rec, int) for elemt_rec in possible_charge_states), 'arg possible_charge_states wrong type'
        cdef libcpp_vector[int] v0 = possible_charge_states
        self.inst.get().setPossibleChargeStates(v0)
        
    
    def getUnchargedMass(self):
        """
        getUnchargedMass(self) -> float
        Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged
        """
        cdef double _r = self.inst.get().getUnchargedMass()
        py_result = <double>_r
        return py_result
    
    def getIntensity(self):
        """
        getIntensity(self) -> float
        """
        cdef float _r = self.inst.get().getIntensity()
        py_result = <float>_r
        return py_result
    
    def getMZ(self):
        """
        getMZ(self) -> float
        """
        cdef double _r = self.inst.get().getMZ()
        py_result = <double>_r
        return py_result
    
    def setMZ(self, double in_0 ):
        """
        setMZ(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setMZ((<double>in_0))
    
    def setIntensity(self, float in_0 ):
        """
        setIntensity(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setIntensity((<float>in_0))
    
    def getPos(self):
        """
        getPos(self) -> float
        """
        cdef double _r = self.inst.get().getPos()
        py_result = <double>_r
        return py_result
    
    def setPos(self, double pos ):
        """
        setPos(self, pos: float ) -> None
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        self.inst.get().setPos((<double>pos))
    
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
        cdef libcpp_map[_String, libcpp_vector[_CVTerm]].iterator outer_it_1268773863406401753366645 = _r.begin()
        cdef libcpp_vector[_CVTerm].iterator inner_it_1268773863406401753366645
        cdef CVTerm item_1268773863406401753366645
        cdef bytes inner_key_1268773863406401753366645
        cdef list inner_values_1268773863406401753366645
        while outer_it_1268773863406401753366645 != _r.end():
           inner_key_1268773863406401753366645 = deref(outer_it_1268773863406401753366645).first.c_str()
           inner_values_1268773863406401753366645 = []
           inner_it_1268773863406401753366645 = deref(outer_it_1268773863406401753366645).second.begin()
           while inner_it_1268773863406401753366645 != deref(outer_it_1268773863406401753366645).second.end():
               item_1268773863406401753366645 = CVTerm.__new__(CVTerm)
               item_1268773863406401753366645.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it_1268773863406401753366645)))
               inner_values_1268773863406401753366645.append(item_1268773863406401753366645)
               inc(inner_it_1268773863406401753366645)
           py_result[inner_key_1268773863406401753366645] = inner_values_1268773863406401753366645
           inc(outer_it_1268773863406401753366645)
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
        if not isinstance(other, Precursor):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Precursor other_casted = other
        cdef Precursor self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    ActivationMethod = __ActivationMethod 

cdef class TargetedExperiment:
    """
    Cython implementation of _TargetedExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TargetedExperiment.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TargetedExperiment rv = TargetedExperiment.__new__(TargetedExperiment)
       rv.inst = shared_ptr[_TargetedExperiment](new _TargetedExperiment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TargetedExperiment rv = TargetedExperiment.__new__(TargetedExperiment)
       rv.inst = shared_ptr[_TargetedExperiment](new _TargetedExperiment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TargetedExperiment](new _TargetedExperiment())
    
    def _init_1(self, TargetedExperiment in_0 ):
        """
        _init_1(self, in_0: TargetedExperiment ) -> None
        """
        assert isinstance(in_0, TargetedExperiment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TargetedExperiment](new _TargetedExperiment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TargetedExperiment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TargetedExperiment)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __add__(TargetedExperiment self, TargetedExperiment other not None):
        cdef _TargetedExperiment * this = self.inst.get()
        cdef _TargetedExperiment * that = other.inst.get()
        cdef _TargetedExperiment applied = deref(this) + deref(that)
        cdef TargetedExperiment result = TargetedExperiment.__new__(TargetedExperiment)
        result.inst = shared_ptr[_TargetedExperiment](new _TargetedExperiment(applied))
        return result
    
    def __iadd__(TargetedExperiment self, TargetedExperiment other not None):
        cdef _TargetedExperiment * this = self.inst.get()
        cdef _TargetedExperiment * that = other.inst.get()
        _iadd(this, that)
        return self
    
    def clear(self, bool clear_meta_data ):
        """
        clear(self, clear_meta_data: bool ) -> None
        """
        assert isinstance(clear_meta_data, pybool_t), 'arg clear_meta_data wrong type'
    
        self.inst.get().clear((<bool>clear_meta_data))
    
    def sortTransitionsByProductMZ(self):
        """
        sortTransitionsByProductMZ(self) -> None
        """
        self.inst.get().sortTransitionsByProductMZ()
    
    def setCVs(self, list cvs ):
        """
        setCVs(self, cvs: List[CV] ) -> None
        """
        assert isinstance(cvs, list) and all(isinstance(elemt_rec, CV) for elemt_rec in cvs), 'arg cvs wrong type'
        cdef libcpp_vector[_CV] * v0 = new libcpp_vector[_CV]()
        cdef CV item0
        for item0 in cvs:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setCVs(deref(v0))
        del v0
    
    def getCVs(self):
        """
        getCVs(self) -> List[CV]
        """
        _r = self.inst.get().getCVs()
        py_result = []
        cdef libcpp_vector[_CV].iterator it__r = _r.begin()
        cdef CV item_py_result
        while it__r != _r.end():
           item_py_result = CV.__new__(CV)
           item_py_result.inst = shared_ptr[_CV](new _CV(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addCV(self, CV cv ):
        """
        addCV(self, cv: CV ) -> None
        """
        assert isinstance(cv, CV), 'arg cv wrong type'
    
        self.inst.get().addCV((deref(cv.inst.get())))
    
    def setContacts(self, list contacts ):
        """
        setContacts(self, contacts: List[Contact] ) -> None
        """
        assert isinstance(contacts, list) and all(isinstance(elemt_rec, Contact) for elemt_rec in contacts), 'arg contacts wrong type'
        cdef libcpp_vector[_Contact] * v0 = new libcpp_vector[_Contact]()
        cdef Contact item0
        for item0 in contacts:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setContacts(deref(v0))
        del v0
    
    def getContacts(self):
        """
        getContacts(self) -> List[Contact]
        """
        _r = self.inst.get().getContacts()
        py_result = []
        cdef libcpp_vector[_Contact].iterator it__r = _r.begin()
        cdef Contact item_py_result
        while it__r != _r.end():
           item_py_result = Contact.__new__(Contact)
           item_py_result.inst = shared_ptr[_Contact](new _Contact(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addContact(self, Contact contact ):
        """
        addContact(self, contact: Contact ) -> None
        """
        assert isinstance(contact, Contact), 'arg contact wrong type'
    
        self.inst.get().addContact((deref(contact.inst.get())))
    
    def setPublications(self, list publications ):
        """
        setPublications(self, publications: List[Publication] ) -> None
        """
        assert isinstance(publications, list) and all(isinstance(elemt_rec, Publication) for elemt_rec in publications), 'arg publications wrong type'
        cdef libcpp_vector[_Publication] * v0 = new libcpp_vector[_Publication]()
        cdef Publication item0
        for item0 in publications:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setPublications(deref(v0))
        del v0
    
    def getPublications(self):
        """
        getPublications(self) -> List[Publication]
        """
        _r = self.inst.get().getPublications()
        py_result = []
        cdef libcpp_vector[_Publication].iterator it__r = _r.begin()
        cdef Publication item_py_result
        while it__r != _r.end():
           item_py_result = Publication.__new__(Publication)
           item_py_result.inst = shared_ptr[_Publication](new _Publication(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addPublication(self, Publication publication ):
        """
        addPublication(self, publication: Publication ) -> None
        """
        assert isinstance(publication, Publication), 'arg publication wrong type'
    
        self.inst.get().addPublication((deref(publication.inst.get())))
    
    def setTargetCVTerms(self, CVTermList cv_terms ):
        """
        setTargetCVTerms(self, cv_terms: CVTermList ) -> None
        """
        assert isinstance(cv_terms, CVTermList), 'arg cv_terms wrong type'
    
        self.inst.get().setTargetCVTerms((deref(cv_terms.inst.get())))
    
    def getTargetCVTerms(self):
        """
        getTargetCVTerms(self) -> CVTermList
        """
        cdef _CVTermList * _r = new _CVTermList(self.inst.get().getTargetCVTerms())
        cdef CVTermList py_result = CVTermList.__new__(CVTermList)
        py_result.inst = shared_ptr[_CVTermList](_r)
        return py_result
    
    def addTargetCVTerm(self, CVTerm cv_term ):
        """
        addTargetCVTerm(self, cv_term: CVTerm ) -> None
        """
        assert isinstance(cv_term, CVTerm), 'arg cv_term wrong type'
    
        self.inst.get().addTargetCVTerm((deref(cv_term.inst.get())))
    
    def setTargetMetaValue(self,  name ,  value ):
        """
        setTargetMetaValue(self, name: Union[bytes, str, String] , value: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(value, (int, float, list, bytes, str)), 'arg value wrong type'
    
    
        self.inst.get().setTargetMetaValue(deref((convString(name)).get()), deref(DataValue(value).inst.get()))
    
    def setInstruments(self, list instruments ):
        """
        setInstruments(self, instruments: List[TargetedExperiment_Instrument] ) -> None
        """
        assert isinstance(instruments, list) and all(isinstance(elemt_rec, TargetedExperiment_Instrument) for elemt_rec in instruments), 'arg instruments wrong type'
        cdef libcpp_vector[_TargetedExperiment_Instrument] * v0 = new libcpp_vector[_TargetedExperiment_Instrument]()
        cdef TargetedExperiment_Instrument item0
        for item0 in instruments:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setInstruments(deref(v0))
        del v0
    
    def getInstruments(self):
        """
        getInstruments(self) -> List[TargetedExperiment_Instrument]
        """
        _r = self.inst.get().getInstruments()
        py_result = []
        cdef libcpp_vector[_TargetedExperiment_Instrument].iterator it__r = _r.begin()
        cdef TargetedExperiment_Instrument item_py_result
        while it__r != _r.end():
           item_py_result = TargetedExperiment_Instrument.__new__(TargetedExperiment_Instrument)
           item_py_result.inst = shared_ptr[_TargetedExperiment_Instrument](new _TargetedExperiment_Instrument(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addInstrument(self, TargetedExperiment_Instrument instrument ):
        """
        addInstrument(self, instrument: TargetedExperiment_Instrument ) -> None
        """
        assert isinstance(instrument, TargetedExperiment_Instrument), 'arg instrument wrong type'
    
        self.inst.get().addInstrument((deref(instrument.inst.get())))
    
    def setSoftware(self, list software ):
        """
        setSoftware(self, software: List[Software] ) -> None
        """
        assert isinstance(software, list) and all(isinstance(elemt_rec, Software) for elemt_rec in software), 'arg software wrong type'
        cdef libcpp_vector[_Software] * v0 = new libcpp_vector[_Software]()
        cdef Software item0
        for item0 in software:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSoftware(deref(v0))
        del v0
    
    def getSoftware(self):
        """
        getSoftware(self) -> List[Software]
        """
        _r = self.inst.get().getSoftware()
        py_result = []
        cdef libcpp_vector[_Software].iterator it__r = _r.begin()
        cdef Software item_py_result
        while it__r != _r.end():
           item_py_result = Software.__new__(Software)
           item_py_result.inst = shared_ptr[_Software](new _Software(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addSoftware(self, Software software ):
        """
        addSoftware(self, software: Software ) -> None
        """
        assert isinstance(software, Software), 'arg software wrong type'
    
        self.inst.get().addSoftware((deref(software.inst.get())))
    
    def setProteins(self, list proteins ):
        """
        setProteins(self, proteins: List[Protein] ) -> None
        """
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, Protein) for elemt_rec in proteins), 'arg proteins wrong type'
        cdef libcpp_vector[_Protein] * v0 = new libcpp_vector[_Protein]()
        cdef Protein item0
        for item0 in proteins:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setProteins(deref(v0))
        del v0
    
    def getProteins(self):
        """
        getProteins(self) -> List[Protein]
        """
        _r = self.inst.get().getProteins()
        py_result = []
        cdef libcpp_vector[_Protein].iterator it__r = _r.begin()
        cdef Protein item_py_result
        while it__r != _r.end():
           item_py_result = Protein.__new__(Protein)
           item_py_result.inst = shared_ptr[_Protein](new _Protein(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getProteinByRef(self,  ref ):
        """
        getProteinByRef(self, ref: Union[bytes, str, String] ) -> Protein
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef _Protein * _r = new _Protein(self.inst.get().getProteinByRef(deref((convString(ref)).get())))
        cdef Protein py_result = Protein.__new__(Protein)
        py_result.inst = shared_ptr[_Protein](_r)
        return py_result
    
    def hasProtein(self,  ref ):
        """
        hasProtein(self, ref: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef bool _r = self.inst.get().hasProtein(deref((convString(ref)).get()))
        py_result = <bool>_r
        return py_result
    
    def addProtein(self, Protein protein ):
        """
        addProtein(self, protein: Protein ) -> None
        """
        assert isinstance(protein, Protein), 'arg protein wrong type'
    
        self.inst.get().addProtein((deref(protein.inst.get())))
    
    def setCompounds(self, list rhs ):
        """
        setCompounds(self, rhs: List[Compound] ) -> None
        """
        assert isinstance(rhs, list) and all(isinstance(elemt_rec, Compound) for elemt_rec in rhs), 'arg rhs wrong type'
        cdef libcpp_vector[_Compound] * v0 = new libcpp_vector[_Compound]()
        cdef Compound item0
        for item0 in rhs:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setCompounds(deref(v0))
        del v0
    
    def getCompounds(self):
        """
        getCompounds(self) -> List[Compound]
        """
        _r = self.inst.get().getCompounds()
        py_result = []
        cdef libcpp_vector[_Compound].iterator it__r = _r.begin()
        cdef Compound item_py_result
        while it__r != _r.end():
           item_py_result = Compound.__new__(Compound)
           item_py_result.inst = shared_ptr[_Compound](new _Compound(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addCompound(self, Compound rhs ):
        """
        addCompound(self, rhs: Compound ) -> None
        """
        assert isinstance(rhs, Compound), 'arg rhs wrong type'
    
        self.inst.get().addCompound((deref(rhs.inst.get())))
    
    def hasCompound(self,  ref ):
        """
        hasCompound(self, ref: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef bool _r = self.inst.get().hasCompound(deref((convString(ref)).get()))
        py_result = <bool>_r
        return py_result
    
    def getCompoundByRef(self,  ref ):
        """
        getCompoundByRef(self, ref: Union[bytes, str, String] ) -> Compound
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef _Compound * _r = new _Compound(self.inst.get().getCompoundByRef(deref((convString(ref)).get())))
        cdef Compound py_result = Compound.__new__(Compound)
        py_result.inst = shared_ptr[_Compound](_r)
        return py_result
    
    def setPeptides(self, list rhs ):
        """
        setPeptides(self, rhs: List[Peptide] ) -> None
        """
        assert isinstance(rhs, list) and all(isinstance(elemt_rec, Peptide) for elemt_rec in rhs), 'arg rhs wrong type'
        cdef libcpp_vector[_Peptide] * v0 = new libcpp_vector[_Peptide]()
        cdef Peptide item0
        for item0 in rhs:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setPeptides(deref(v0))
        del v0
    
    def getPeptides(self):
        """
        getPeptides(self) -> List[Peptide]
        """
        _r = self.inst.get().getPeptides()
        py_result = []
        cdef libcpp_vector[_Peptide].iterator it__r = _r.begin()
        cdef Peptide item_py_result
        while it__r != _r.end():
           item_py_result = Peptide.__new__(Peptide)
           item_py_result.inst = shared_ptr[_Peptide](new _Peptide(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def hasPeptide(self,  ref ):
        """
        hasPeptide(self, ref: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef bool _r = self.inst.get().hasPeptide(deref((convString(ref)).get()))
        py_result = <bool>_r
        return py_result
    
    def getPeptideByRef(self,  ref ):
        """
        getPeptideByRef(self, ref: Union[bytes, str, String] ) -> Peptide
        """
        assert (isinstance(ref, str) or isinstance(ref, bytes) or isinstance(ref, String)), 'arg ref wrong type'
    
        cdef _Peptide * _r = new _Peptide(self.inst.get().getPeptideByRef(deref((convString(ref)).get())))
        cdef Peptide py_result = Peptide.__new__(Peptide)
        py_result.inst = shared_ptr[_Peptide](_r)
        return py_result
    
    def addPeptide(self, Peptide rhs ):
        """
        addPeptide(self, rhs: Peptide ) -> None
        """
        assert isinstance(rhs, Peptide), 'arg rhs wrong type'
    
        self.inst.get().addPeptide((deref(rhs.inst.get())))
    
    def setTransitions(self, list transitions ):
        """
        setTransitions(self, transitions: List[ReactionMonitoringTransition] ) -> None
        """
        assert isinstance(transitions, list) and all(isinstance(elemt_rec, ReactionMonitoringTransition) for elemt_rec in transitions), 'arg transitions wrong type'
        cdef libcpp_vector[_ReactionMonitoringTransition] * v0 = new libcpp_vector[_ReactionMonitoringTransition]()
        cdef ReactionMonitoringTransition item0
        for item0 in transitions:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setTransitions(deref(v0))
        del v0
    
    def getTransitions(self):
        """
        getTransitions(self) -> List[ReactionMonitoringTransition]
        """
        _r = self.inst.get().getTransitions()
        py_result = []
        cdef libcpp_vector[_ReactionMonitoringTransition].iterator it__r = _r.begin()
        cdef ReactionMonitoringTransition item_py_result
        while it__r != _r.end():
           item_py_result = ReactionMonitoringTransition.__new__(ReactionMonitoringTransition)
           item_py_result.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addTransition(self, ReactionMonitoringTransition transition ):
        """
        addTransition(self, transition: ReactionMonitoringTransition ) -> None
        """
        assert isinstance(transition, ReactionMonitoringTransition), 'arg transition wrong type'
    
        self.inst.get().addTransition((deref(transition.inst.get())))
    
    def setIncludeTargets(self, list targets ):
        """
        setIncludeTargets(self, targets: List[IncludeExcludeTarget] ) -> None
        """
        assert isinstance(targets, list) and all(isinstance(elemt_rec, IncludeExcludeTarget) for elemt_rec in targets), 'arg targets wrong type'
        cdef libcpp_vector[_IncludeExcludeTarget] * v0 = new libcpp_vector[_IncludeExcludeTarget]()
        cdef IncludeExcludeTarget item0
        for item0 in targets:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIncludeTargets(deref(v0))
        del v0
    
    def getIncludeTargets(self):
        """
        getIncludeTargets(self) -> List[IncludeExcludeTarget]
        """
        _r = self.inst.get().getIncludeTargets()
        py_result = []
        cdef libcpp_vector[_IncludeExcludeTarget].iterator it__r = _r.begin()
        cdef IncludeExcludeTarget item_py_result
        while it__r != _r.end():
           item_py_result = IncludeExcludeTarget.__new__(IncludeExcludeTarget)
           item_py_result.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addIncludeTarget(self, IncludeExcludeTarget target ):
        """
        addIncludeTarget(self, target: IncludeExcludeTarget ) -> None
        """
        assert isinstance(target, IncludeExcludeTarget), 'arg target wrong type'
    
        self.inst.get().addIncludeTarget((deref(target.inst.get())))
    
    def setExcludeTargets(self, list targets ):
        """
        setExcludeTargets(self, targets: List[IncludeExcludeTarget] ) -> None
        """
        assert isinstance(targets, list) and all(isinstance(elemt_rec, IncludeExcludeTarget) for elemt_rec in targets), 'arg targets wrong type'
        cdef libcpp_vector[_IncludeExcludeTarget] * v0 = new libcpp_vector[_IncludeExcludeTarget]()
        cdef IncludeExcludeTarget item0
        for item0 in targets:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setExcludeTargets(deref(v0))
        del v0
    
    def getExcludeTargets(self):
        """
        getExcludeTargets(self) -> List[IncludeExcludeTarget]
        """
        _r = self.inst.get().getExcludeTargets()
        py_result = []
        cdef libcpp_vector[_IncludeExcludeTarget].iterator it__r = _r.begin()
        cdef IncludeExcludeTarget item_py_result
        while it__r != _r.end():
           item_py_result = IncludeExcludeTarget.__new__(IncludeExcludeTarget)
           item_py_result.inst = shared_ptr[_IncludeExcludeTarget](new _IncludeExcludeTarget(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addExcludeTarget(self, IncludeExcludeTarget target ):
        """
        addExcludeTarget(self, target: IncludeExcludeTarget ) -> None
        """
        assert isinstance(target, IncludeExcludeTarget), 'arg target wrong type'
    
        self.inst.get().addExcludeTarget((deref(target.inst.get())))
    
    def setSourceFiles(self, list source_files ):
        """
        setSourceFiles(self, source_files: List[SourceFile] ) -> None
        """
        assert isinstance(source_files, list) and all(isinstance(elemt_rec, SourceFile) for elemt_rec in source_files), 'arg source_files wrong type'
        cdef libcpp_vector[_SourceFile] * v0 = new libcpp_vector[_SourceFile]()
        cdef SourceFile item0
        for item0 in source_files:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSourceFiles(deref(v0))
        del v0
    
    def getSourceFiles(self):
        """
        getSourceFiles(self) -> List[SourceFile]
        """
        _r = self.inst.get().getSourceFiles()
        py_result = []
        cdef libcpp_vector[_SourceFile].iterator it__r = _r.begin()
        cdef SourceFile item_py_result
        while it__r != _r.end():
           item_py_result = SourceFile.__new__(SourceFile)
           item_py_result.inst = shared_ptr[_SourceFile](new _SourceFile(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addSourceFile(self, SourceFile source_file ):
        """
        addSourceFile(self, source_file: SourceFile ) -> None
        """
        assert isinstance(source_file, SourceFile), 'arg source_file wrong type'
    
        self.inst.get().addSourceFile((deref(source_file.inst.get())))
    
    def containsInvalidReferences(self):
        """
        containsInvalidReferences(self) -> bool
        """
        cdef bool _r = self.inst.get().containsInvalidReferences()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, TargetedExperiment):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TargetedExperiment other_casted = other
        cdef TargetedExperiment self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class WindowMower:
    """
    Cython implementation of _WindowMower

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1WindowMower.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef WindowMower rv = WindowMower.__new__(WindowMower)
       rv.inst = shared_ptr[_WindowMower](new _WindowMower(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef WindowMower rv = WindowMower.__new__(WindowMower)
       rv.inst = shared_ptr[_WindowMower](new _WindowMower(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_WindowMower](new _WindowMower())
    
    def _init_1(self, WindowMower in_0 ):
        """
        _init_1(self, in_0: WindowMower ) -> None
        """
        assert isinstance(in_0, WindowMower), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_WindowMower](new _WindowMower((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: WindowMower ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], WindowMower)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filterPeakSpectrumForTopNInSlidingWindow(self, MSSpectrum spectrum ):
        """
        filterPeakSpectrumForTopNInSlidingWindow(self, spectrum: MSSpectrum ) -> None
        Sliding window version (slower)
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().filterPeakSpectrumForTopNInSlidingWindow((deref(spectrum.inst.get())))
    
    def filterPeakSpectrumForTopNInJumpingWindow(self, MSSpectrum spectrum ):
        """
        filterPeakSpectrumForTopNInJumpingWindow(self, spectrum: MSSpectrum ) -> None
        Jumping window version (faster)
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().filterPeakSpectrumForTopNInJumpingWindow((deref(spectrum.inst.get())))
    
    def filterPeakSpectrum(self, MSSpectrum spec ):
        """
        filterPeakSpectrum(self, spec: MSSpectrum ) -> None
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterPeakSpectrum((deref(spec.inst.get())))
    
    def filterPeakMap(self, MSExperiment exp ):
        """
        filterPeakMap(self, exp: MSExperiment ) -> None
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
