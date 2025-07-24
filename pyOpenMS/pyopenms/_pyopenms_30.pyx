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
def __static_AASequence_fromString( s ):
    """
    __static_AASequence_fromString(s: Union[bytes, str, String] ) -> AASequence
        deprecated. Use AASequence(String) instead.
    """
    assert (isinstance(s, str) or isinstance(s, bytes) or isinstance(s, String)), 'arg s wrong type'

    cdef _AASequence * _r = new _AASequence(_fromString_AASequence(deref((convString(s)).get())))
    cdef AASequence py_result = AASequence.__new__(AASequence)
    py_result.inst = shared_ptr[_AASequence](_r)
    return py_result

def __static_AASequence_fromStringPermissive( s , bool permissive ):
    """
    __static_AASequence_fromStringPermissive(s: Union[bytes, str, String] , permissive: bool ) -> AASequence
        deprecated. Use AASequence(String, bool) instead.
    """
    assert (isinstance(s, str) or isinstance(s, bytes) or isinstance(s, String)), 'arg s wrong type'
    assert isinstance(permissive, pybool_t), 'arg permissive wrong type'


    cdef _AASequence * _r = new _AASequence(_fromString_AASequence(deref((convString(s)).get()), (<bool>permissive)))
    cdef AASequence py_result = AASequence.__new__(AASequence)
    py_result.inst = shared_ptr[_AASequence](_r)
    return py_result

def __static_DateTime_now():
    """
    __static_DateTime_now() -> DateTime
    """
    cdef _DateTime * _r = new _DateTime(_now_DateTime())
    cdef DateTime py_result = DateTime.__new__(DateTime)
    py_result.inst = shared_ptr[_DateTime](_r)
    return py_result 

cdef class SIDE:
    None
    LEFT = 0
    RIGHT = 1
    BOTH = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AAIndex:
    """
    Cython implementation of _AAIndex

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AAIndex.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def aliphatic(self, bytes aa ):
        """
        aliphatic(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().aliphatic((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def acidic(self, bytes aa ):
        """
        acidic(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().acidic((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def basic(self, bytes aa ):
        """
        basic(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().basic((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def polar(self, bytes aa ):
        """
        polar(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().polar((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getKHAG800101(self, bytes aa ):
        """
        getKHAG800101(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getKHAG800101((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getVASM830103(self, bytes aa ):
        """
        getVASM830103(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getVASM830103((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getNADH010106(self, bytes aa ):
        """
        getNADH010106(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getNADH010106((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getNADH010107(self, bytes aa ):
        """
        getNADH010107(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getNADH010107((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getWILM950102(self, bytes aa ):
        """
        getWILM950102(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getWILM950102((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getROBB760107(self, bytes aa ):
        """
        getROBB760107(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getROBB760107((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getOOBM850104(self, bytes aa ):
        """
        getOOBM850104(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getOOBM850104((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getFAUJ880111(self, bytes aa ):
        """
        getFAUJ880111(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getFAUJ880111((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getFINA770101(self, bytes aa ):
        """
        getFINA770101(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getFINA770101((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def getARGP820102(self, bytes aa ):
        """
        getARGP820102(self, aa: bytes ) -> float
        """
        assert isinstance(aa, bytes) and len(aa) == 1, 'arg aa wrong type'
    
        cdef double _r = self.inst.get().getARGP820102((<char>((aa)[0])))
        py_result = <double>_r
        return py_result
    
    def calculateGB(self, AASequence seq , double T ):
        """
        calculateGB(self, seq: AASequence , T: float ) -> float
        """
        assert isinstance(seq, AASequence), 'arg seq wrong type'
        assert isinstance(T, float), 'arg T wrong type'
    
    
        cdef double _r = self.inst.get().calculateGB((deref(seq.inst.get())), (<double>T))
        py_result = <double>_r
        return py_result 

cdef class AASequence:
    """
    Cython implementation of _AASequence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AASequence.html>`_

    Representation of a peptide/protein sequence
    This class represents amino acid sequences in OpenMS. An AASequence
    instance primarily contains a sequence of residues.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).toString().c_str() )

    
    def __copy__(self):
       cdef AASequence rv = AASequence.__new__(AASequence)
       rv.inst = shared_ptr[_AASequence](new _AASequence(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AASequence rv = AASequence.__new__(AASequence)
       rv.inst = shared_ptr[_AASequence](new _AASequence(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AASequence](new _AASequence())
    
    def _init_1(self, AASequence in_0 ):
        """
        _init_1(self, in_0: AASequence ) -> None
        """
        assert isinstance(in_0, AASequence), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AASequence](new _AASequence((deref(in_0.inst.get()))))
    
    def _init_2(self,  in_0 ):
        """
        _init_2(self, in_0: Union[bytes, str, String] ) -> None
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE")
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AASequence](new _AASequence(deref((convString(in_0)).get())))
    
    def _init_3(self,  in_0 , bool permissive ):
        """
        _init_3(self, in_0: Union[bytes, str, String] , permissive: bool ) -> None
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE"), permissive allows for '+', '*', and '#' in the sequence
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(permissive, pybool_t), 'arg permissive wrong type'
    
    
        self.inst = shared_ptr[_AASequence](new _AASequence(deref((convString(in_0)).get()), (<bool>permissive)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AASequence ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[bytes, str, String] ) -> None
          :noindex:
        
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE")

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[bytes, str, String] , permissive: bool ) -> None
          :noindex:
        
        Constructor from amino acid sequence (e.g. "PEPTM(Oxidatio)IDE"), permissive allows for '+', '*', and '#' in the sequence
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AASequence)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __add__(AASequence self, AASequence other not None):
        cdef _AASequence * this = self.inst.get()
        cdef _AASequence * that = other.inst.get()
        cdef _AASequence applied = deref(this) + deref(that)
        cdef AASequence result = AASequence.__new__(AASequence)
        result.inst = shared_ptr[_AASequence](new _AASequence(applied))
        return result
    
    def __iadd__(AASequence self, AASequence other not None):
        cdef _AASequence * this = self.inst.get()
        cdef _AASequence * that = other.inst.get()
        _iadd(this, that)
        return self
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> Residue
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _Residue * _r = new _Residue(deref(self.inst.get())[(<size_t>in_0)])
        cdef Residue py_result = Residue.__new__(Residue)
        py_result.inst = shared_ptr[_Residue](_r)
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        Check if sequence is empty
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        Returns the peptide as string with modifications embedded in brackets
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def toUnmodifiedString(self):
        """
        toUnmodifiedString(self) -> Union[bytes, str, String]
        Returns the peptide as string without any modifications
        """
        cdef _String _r = self.inst.get().toUnmodifiedString()
        py_result = convOutputString(_r)
        return py_result
    
    def toUniModString(self):
        """
        toUniModString(self) -> Union[bytes, str, String]
        Returns the peptide as string with UniMod-style modifications embedded in brackets
        """
        cdef _String _r = self.inst.get().toUniModString()
        py_result = convOutputString(_r)
        return py_result
    
    def _toBracketString_0(self):
        """
        _toBracketString_0(self) -> Union[bytes, str, String]
        Create a TPP compatible string of the modified sequence using bracket notation. Uses integer mass by default
        """
        cdef _String _r = self.inst.get().toBracketString()
        py_result = convOutputString(_r)
        return py_result
    
    def _toBracketString_1(self, bool integer_mass ):
        """
        _toBracketString_1(self, integer_mass: bool ) -> Union[bytes, str, String]
        Create a TPP compatible string of the modified sequence using bracket notation
        """
        assert isinstance(integer_mass, pybool_t), 'arg integer_mass wrong type'
    
        cdef _String _r = self.inst.get().toBracketString((<bool>integer_mass))
        py_result = convOutputString(_r)
        return py_result
    
    def _toBracketString_2(self, bool integer_mass , bool mass_delta ):
        """
        _toBracketString_2(self, integer_mass: bool , mass_delta: bool ) -> Union[bytes, str, String]
        Create a TPP compatible string of the modified sequence using bracket notation.
        """
        assert isinstance(integer_mass, pybool_t), 'arg integer_mass wrong type'
        assert isinstance(mass_delta, pybool_t), 'arg mass_delta wrong type'
    
    
        cdef _String _r = self.inst.get().toBracketString((<bool>integer_mass), (<bool>mass_delta))
        py_result = convOutputString(_r)
        return py_result
    
    def _toBracketString_3(self, bool integer_mass , bool mass_delta , list fixed_modifications ):
        """
        _toBracketString_3(self, integer_mass: bool , mass_delta: bool , fixed_modifications: List[bytes] ) -> Union[bytes, str, String]
        Create a TPP compatible string of the modified sequence using bracket notation
        """
        assert isinstance(integer_mass, pybool_t), 'arg integer_mass wrong type'
        assert isinstance(mass_delta, pybool_t), 'arg mass_delta wrong type'
        assert isinstance(fixed_modifications, list) and all(isinstance(i, bytes) for i in fixed_modifications), 'arg fixed_modifications wrong type'
    
    
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in fixed_modifications:
           v2.push_back(_String(<char *>item2))
        cdef _String _r = self.inst.get().toBracketString((<bool>integer_mass), (<bool>mass_delta), deref(v2))
        del v2
        py_result = convOutputString(_r)
        return py_result
    
    def toBracketString(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: toBracketString(self, ) -> Union[bytes, str, String]
          :noindex:
        
        Create a TPP compatible string of the modified sequence using bracket notation. Uses integer mass by default

        
        .. rubric:: Overload:
        .. py:function:: toBracketString(self, integer_mass: bool ) -> Union[bytes, str, String]
          :noindex:
        
        Create a TPP compatible string of the modified sequence using bracket notation

        
        .. rubric:: Overload:
        .. py:function:: toBracketString(self, integer_mass: bool , mass_delta: bool ) -> Union[bytes, str, String]
          :noindex:
        
        Create a TPP compatible string of the modified sequence using bracket notation.

        
        .. rubric:: Overload:
        .. py:function:: toBracketString(self, integer_mass: bool , mass_delta: bool , fixed_modifications: List[bytes] ) -> Union[bytes, str, String]
          :noindex:
        
        Create a TPP compatible string of the modified sequence using bracket notation
    
        """
        if not args:
            return self._toBracketString_0(*args)
        elif (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._toBracketString_1(*args)
        elif (len(args)==2) and (isinstance(args[0], pybool_t)) and (isinstance(args[1], pybool_t)):
            return self._toBracketString_2(*args)
        elif (len(args)==3) and (isinstance(args[0], pybool_t)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], list) and all(isinstance(i, bytes) for i in args[2])):
            return self._toBracketString_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setModification_0(self,  index ,  modification ):
        """
        _setModification_0(self, index: int , modification: Union[bytes, str, String] ) -> None
        Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
    
        self.inst.get().setModification((<size_t>index), deref((convString(modification)).get()))
    
    def _setModification_1(self,  index , ResidueModification modification ):
        """
        _setModification_1(self, index: int , modification: ResidueModification ) -> None
        Sets the modification of AA at index by providing a ResidueModification object. Stricter than just looking for the name and adds the Modification to the DB if not present
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(modification, ResidueModification), 'arg modification wrong type'
    
    
        self.inst.get().setModification((<size_t>index), (deref(modification.inst.get())))
    
    def setModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setModification(self, index: int , modification: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version

        
        .. rubric:: Overload:
        .. py:function:: setModification(self, index: int , modification: ResidueModification ) -> None
          :noindex:
        
        Sets the modification of AA at index by providing a ResidueModification object. Stricter than just looking for the name and adds the Modification to the DB if not present
    
        """
        if (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._setModification_0(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], ResidueModification)):
            return self._setModification_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setModificationByDiffMonoMass(self,  index , double diffMonoMass ):
        """
        setModificationByDiffMonoMass(self, index: int , diffMonoMass: float ) -> None
        Modifies the residue at index in the sequence and potentially in the ResidueDB
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(diffMonoMass, float), 'arg diffMonoMass wrong type'
    
    
        self.inst.get().setModificationByDiffMonoMass((<size_t>index), (<double>diffMonoMass))
    
    def _setNTerminalModification_0(self,  modification ):
        """
        _setNTerminalModification_0(self, modification: Union[bytes, str, String] ) -> None
        Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)
        """
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
        self.inst.get().setNTerminalModification(deref((convString(modification)).get()))
    
    def _setNTerminalModification_1(self, ResidueModification mod ):
        """
        _setNTerminalModification_1(self, mod: ResidueModification ) -> None
        Sets the N-terminal modification (copies and adds to database if not present)
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
    
        self.inst.get().setNTerminalModification((deref(mod.inst.get())))
    
    def setNTerminalModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setNTerminalModification(self, modification: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        
        .. rubric:: Overload:
        .. py:function:: setNTerminalModification(self, mod: ResidueModification ) -> None
          :noindex:
        
        Sets the N-terminal modification (copies and adds to database if not present)
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setNTerminalModification_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ResidueModification)):
            return self._setNTerminalModification_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setNTerminalModificationByDiffMonoMass(self, double diffMonoMass , bool protein_term ):
        """
        setNTerminalModificationByDiffMonoMass(self, diffMonoMass: float , protein_term: bool ) -> None
        Sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
        """
        assert isinstance(diffMonoMass, float), 'arg diffMonoMass wrong type'
        assert isinstance(protein_term, pybool_t), 'arg protein_term wrong type'
    
    
        self.inst.get().setNTerminalModificationByDiffMonoMass((<double>diffMonoMass), (<bool>protein_term))
    
    def _setCTerminalModification_0(self,  modification ):
        """
        _setCTerminalModification_0(self, modification: Union[bytes, str, String] ) -> None
        Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)
        """
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
        self.inst.get().setCTerminalModification(deref((convString(modification)).get()))
    
    def _setCTerminalModification_1(self, ResidueModification mod ):
        """
        _setCTerminalModification_1(self, mod: ResidueModification ) -> None
        Sets the C-terminal modification (copies and adds to database if not present)
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
    
        self.inst.get().setCTerminalModification((deref(mod.inst.get())))
    
    def setCTerminalModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setCTerminalModification(self, modification: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)

        
        .. rubric:: Overload:
        .. py:function:: setCTerminalModification(self, mod: ResidueModification ) -> None
          :noindex:
        
        Sets the C-terminal modification (copies and adds to database if not present)
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setCTerminalModification_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ResidueModification)):
            return self._setCTerminalModification_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCTerminalModificationByDiffMonoMass(self, double diffMonoMass , bool protein_term ):
        """
        setCTerminalModificationByDiffMonoMass(self, diffMonoMass: float , protein_term: bool ) -> None
        Sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
        """
        assert isinstance(diffMonoMass, float), 'arg diffMonoMass wrong type'
        assert isinstance(protein_term, pybool_t), 'arg protein_term wrong type'
    
    
        self.inst.get().setCTerminalModificationByDiffMonoMass((<double>diffMonoMass), (<bool>protein_term))
    
    def getNTerminalModificationName(self):
        """
        getNTerminalModificationName(self) -> Union[bytes, str, String]
        Returns the name (ID) of the N-terminal modification, or an empty string if none is set
        """
        cdef _String _r = self.inst.get().getNTerminalModificationName()
        py_result = convOutputString(_r)
        return py_result
    
    def getNTerminalModification(self):
        """
        getNTerminalModification(self) -> ResidueModification
        Returns a copy of the name N-terminal modification object, or None
        """
        cdef const _ResidueModification * __r = (self.inst.get().getNTerminalModification())
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getCTerminalModificationName(self):
        """
        getCTerminalModificationName(self) -> Union[bytes, str, String]
        Returns the name (ID) of the C-terminal modification, or an empty string if none is set
        """
        cdef _String _r = self.inst.get().getCTerminalModificationName()
        py_result = convOutputString(_r)
        return py_result
    
    def getCTerminalModification(self):
        """
        getCTerminalModification(self) -> ResidueModification
        Returns a copy of the name C-terminal modification object, or None
        """
        cdef const _ResidueModification * __r = (self.inst.get().getCTerminalModification())
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getResidue(self,  index ):
        """
        getResidue(self, index: int ) -> Residue
        Returns the residue at position index
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef _Residue * _r = new _Residue(self.inst.get().getResidue((<size_t>index)))
        cdef Residue py_result = Residue.__new__(Residue)
        py_result.inst = shared_ptr[_Residue](_r)
        return py_result
    
    def _getFormula_0(self):
        """
        _getFormula_0(self) -> EmpiricalFormula
        Convenience function with ResidueType=Full and charge = 0 by default
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def _getFormula_1(self, int type_ ,  charge ):
        """
        _getFormula_1(self, type_: int , charge: int ) -> EmpiricalFormula
        """
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula((<_ResidueType>type_), (<int>charge)))
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getFormula(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getFormula(self, ) -> EmpiricalFormula
          :noindex:
        
        Convenience function with ResidueType=Full and charge = 0 by default

        
        .. rubric:: Overload:
        .. py:function:: getFormula(self, type_: int , charge: int ) -> EmpiricalFormula
          :noindex:
    
        """
        if not args:
            return self._getFormula_0(*args)
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]) and (isinstance(args[1], int)):
            return self._getFormula_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getAverageWeight_0(self):
        """
        _getAverageWeight_0(self) -> float
        Returns the average weight of the peptide
        """
        cdef double _r = self.inst.get().getAverageWeight()
        py_result = <double>_r
        return py_result
    
    def _getAverageWeight_1(self, int type_ ,  charge ):
        """
        _getAverageWeight_1(self, type_: int , charge: int ) -> float
        """
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef double _r = self.inst.get().getAverageWeight((<_ResidueType>type_), (<int>charge))
        py_result = <double>_r
        return py_result
    
    def getAverageWeight(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getAverageWeight(self, ) -> float
          :noindex:
        
        Returns the average weight of the peptide

        
        .. rubric:: Overload:
        .. py:function:: getAverageWeight(self, type_: int , charge: int ) -> float
          :noindex:
    
        """
        if not args:
            return self._getAverageWeight_0(*args)
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]) and (isinstance(args[1], int)):
            return self._getAverageWeight_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getMonoWeight_0(self):
        """
        _getMonoWeight_0(self) -> float
        Returns the mono isotopic weight of the peptide
        """
        cdef double _r = self.inst.get().getMonoWeight()
        py_result = <double>_r
        return py_result
    
    def _getMonoWeight_1(self, int type_ ,  charge ):
        """
        _getMonoWeight_1(self, type_: int , charge: int ) -> float
        """
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef double _r = self.inst.get().getMonoWeight((<_ResidueType>type_), (<int>charge))
        py_result = <double>_r
        return py_result
    
    def getMonoWeight(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMonoWeight(self, ) -> float
          :noindex:
        
        Returns the mono isotopic weight of the peptide

        
        .. rubric:: Overload:
        .. py:function:: getMonoWeight(self, type_: int , charge: int ) -> float
          :noindex:
    
        """
        if not args:
            return self._getMonoWeight_0(*args)
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]) and (isinstance(args[1], int)):
            return self._getMonoWeight_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getMZ_0(self,  charge ):
        """
        _getMZ_0(self, charge: int ) -> float
        Returns the mass-to-charge ratio of the peptide
        """
        assert isinstance(charge, int), 'arg charge wrong type'
    
        cdef double _r = self.inst.get().getMZ((<int>charge))
        py_result = <double>_r
        return py_result
    
    def _getMZ_1(self,  charge , int type_ ):
        """
        _getMZ_1(self, charge: int , type_: int ) -> float
        """
        assert isinstance(charge, int), 'arg charge wrong type'
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg type_ wrong type'
    
    
        cdef double _r = self.inst.get().getMZ((<int>charge), (<_ResidueType>type_))
        py_result = <double>_r
        return py_result
    
    def getMZ(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMZ(self, charge: int ) -> float
          :noindex:
        
        Returns the mass-to-charge ratio of the peptide

        
        .. rubric:: Overload:
        .. py:function:: getMZ(self, charge: int , type_: int ) -> float
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int)):
            return self._getMZ_0(*args)
        elif (len(args)==2) and (isinstance(args[0], int)) and (args[1] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]):
            return self._getMZ_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        Returns the number of residues
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def getPrefix(self,  index ):
        """
        getPrefix(self, index: int ) -> AASequence
        Returns a peptide sequence of the first index residues
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().getPrefix((<size_t>index)))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def getSuffix(self,  index ):
        """
        getSuffix(self, index: int ) -> AASequence
        Returns a peptide sequence of the last index residues
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().getSuffix((<size_t>index)))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def getSubsequence(self,  index ,  number ):
        """
        getSubsequence(self, index: int , number: int ) -> AASequence
        Returns a peptide sequence of number residues, beginning at position index
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(number, int), 'arg number wrong type'
    
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().getSubsequence((<size_t>index), (<unsigned int>number)))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def has(self, Residue residue ):
        """
        has(self, residue: Residue ) -> bool
        Returns true if the peptide contains the given residue
        """
        assert isinstance(residue, Residue), 'arg residue wrong type'
    
        cdef bool _r = self.inst.get().has((deref(residue.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def hasSubsequence(self, AASequence peptide ):
        """
        hasSubsequence(self, peptide: AASequence ) -> bool
        Returns true if the peptide contains the given peptide
        """
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
    
        cdef bool _r = self.inst.get().hasSubsequence((deref(peptide.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def hasPrefix(self, AASequence peptide ):
        """
        hasPrefix(self, peptide: AASequence ) -> bool
        Returns true if the peptide has the given prefix
        """
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
    
        cdef bool _r = self.inst.get().hasPrefix((deref(peptide.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def hasSuffix(self, AASequence peptide ):
        """
        hasSuffix(self, peptide: AASequence ) -> bool
        Returns true if the peptide has the given suffix
        """
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
    
        cdef bool _r = self.inst.get().hasSuffix((deref(peptide.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def hasNTerminalModification(self):
        """
        hasNTerminalModification(self) -> bool
        Predicate which is true if the peptide is N-term modified
        """
        cdef bool _r = self.inst.get().hasNTerminalModification()
        py_result = <bool>_r
        return py_result
    
    def hasCTerminalModification(self):
        """
        hasCTerminalModification(self) -> bool
        Predicate which is true if the peptide is C-term modified
        """
        cdef bool _r = self.inst.get().hasCTerminalModification()
        py_result = <bool>_r
        return py_result
    
    def isModified(self):
        """
        isModified(self) -> bool
        Returns true if any of the residues or termini are modified
        """
        cdef bool _r = self.inst.get().isModified()
        py_result = <bool>_r
        return py_result
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        Returns the peptide as string with modifications embedded in brackets
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, AASequence):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef AASequence other_casted = other
        cdef AASequence self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def getAAFrequencies(self, dict mmap):
        cdef libcpp_map[_String, size_t] c_mmap
        self.inst.get().getAAFrequencies(c_mmap)
        for k,v in mmap.iteritems():
            v = c_mmap[ _String(<char *>k) ]

    def __iter__(self):
        cdef unsigned int n = self.inst.get().size()
        cdef unsigned int i = 0

        cdef Residue py_result 
        while i < n:

            py_result = Residue.__new__(Residue)
            py_result.inst = shared_ptr[_Residue](new _Residue( deref(self.inst.get())[i] ))

            yield py_result

            i += 1
    fromString = __static_AASequence_fromString
    fromStringPermissive = __static_AASequence_fromStringPermissive 

cdef class Compomer:
    """
    Cython implementation of _Compomer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Compomer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Compomer rv = Compomer.__new__(Compomer)
       rv.inst = shared_ptr[_Compomer](new _Compomer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Compomer rv = Compomer.__new__(Compomer)
       rv.inst = shared_ptr[_Compomer](new _Compomer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Compomer](new _Compomer())
    
    def _init_1(self, Compomer in_0 ):
        """
        _init_1(self, in_0: Compomer ) -> None
        """
        assert isinstance(in_0, Compomer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Compomer](new _Compomer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Compomer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Compomer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def add(self, Adduct a ,  side ):
        """
        add(self, a: Adduct , side: int ) -> None
        """
        assert isinstance(a, Adduct), 'arg a wrong type'
        assert isinstance(side, int), 'arg side wrong type'
    
    
        self.inst.get().add((deref(a.inst.get())), (<unsigned int>side))
    
    def isConflicting(self, Compomer cmp ,  side_this ,  side_other ):
        """
        isConflicting(self, cmp: Compomer , side_this: int , side_other: int ) -> bool
        """
        assert isinstance(cmp, Compomer), 'arg cmp wrong type'
        assert isinstance(side_this, int), 'arg side_this wrong type'
        assert isinstance(side_other, int), 'arg side_other wrong type'
    
    
    
        cdef bool _r = self.inst.get().isConflicting((deref(cmp.inst.get())), (<unsigned int>side_this), (<unsigned int>side_other))
        py_result = <bool>_r
        return py_result
    
    def setID(self,  id ):
        """
        setID(self, id: int ) -> None
        Sets an Id which allows unique identification of a compomer
        """
        assert isinstance(id, int) and id >= 0, 'arg id wrong type'
    
        self.inst.get().setID((<size_t>id))
    
    def getID(self):
        """
        getID(self) -> int
        Returns Id which allows unique identification of this compomer
        """
        cdef size_t _r = self.inst.get().getID()
        py_result = <size_t>_r
        return py_result
    
    def getNetCharge(self):
        """
        getNetCharge(self) -> int
        Net charge of compomer (i.e. difference between left and right side of compomer)
        """
        cdef int _r = self.inst.get().getNetCharge()
        py_result = <int>_r
        return py_result
    
    def getMass(self):
        """
        getMass(self) -> float
        Mass of all contained adducts
        """
        cdef double _r = self.inst.get().getMass()
        py_result = <double>_r
        return py_result
    
    def getPositiveCharges(self):
        """
        getPositiveCharges(self) -> int
        Summed positive charges of contained adducts
        """
        cdef int _r = self.inst.get().getPositiveCharges()
        py_result = <int>_r
        return py_result
    
    def getNegativeCharges(self):
        """
        getNegativeCharges(self) -> int
        Summed negative charges of contained adducts
        """
        cdef int _r = self.inst.get().getNegativeCharges()
        py_result = <int>_r
        return py_result
    
    def getLogP(self):
        """
        getLogP(self) -> float
        Returns the log probability
        """
        cdef double _r = self.inst.get().getLogP()
        py_result = <double>_r
        return py_result
    
    def getRTShift(self):
        """
        getRTShift(self) -> float
        Returns the log probability
        """
        cdef double _r = self.inst.get().getRTShift()
        py_result = <double>_r
        return py_result
    
    def _getAdductsAsString_0(self):
        """
        _getAdductsAsString_0(self) -> Union[bytes, str, String]
        Get adducts with their abundance as compact string for both sides
        """
        cdef _String _r = self.inst.get().getAdductsAsString()
        py_result = convOutputString(_r)
        return py_result
    
    def _getAdductsAsString_1(self,  side ):
        """
        _getAdductsAsString_1(self, side: int ) -> Union[bytes, str, String]
        Get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
        """
        assert isinstance(side, int), 'arg side wrong type'
    
        cdef _String _r = self.inst.get().getAdductsAsString((<unsigned int>side))
        py_result = convOutputString(_r)
        return py_result
    
    def getAdductsAsString(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getAdductsAsString(self, ) -> Union[bytes, str, String]
          :noindex:
        
        Get adducts with their abundance as compact string for both sides

        
        .. rubric:: Overload:
        .. py:function:: getAdductsAsString(self, side: int ) -> Union[bytes, str, String]
          :noindex:
        
        Get adducts with their abundance as compact string (amounts are absolute unless side=BOTH)
    
        """
        if not args:
            return self._getAdductsAsString_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._getAdductsAsString_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isSingleAdduct(self, Adduct a ,  side ):
        """
        isSingleAdduct(self, a: Adduct , side: int ) -> bool
        Check if Compomer only contains a single adduct on side @p side
        """
        assert isinstance(a, Adduct), 'arg a wrong type'
        assert isinstance(side, int), 'arg side wrong type'
    
    
        cdef bool _r = self.inst.get().isSingleAdduct((deref(a.inst.get())), (<unsigned int>side))
        py_result = <bool>_r
        return py_result
    
    def _removeAdduct_0(self, Adduct a ):
        """
        _removeAdduct_0(self, a: Adduct ) -> Compomer
        Remove ALL instances of the given adduct
        """
        assert isinstance(a, Adduct), 'arg a wrong type'
    
        cdef _Compomer * _r = new _Compomer(self.inst.get().removeAdduct((deref(a.inst.get()))))
        cdef Compomer py_result = Compomer.__new__(Compomer)
        py_result.inst = shared_ptr[_Compomer](_r)
        return py_result
    
    def _removeAdduct_1(self, Adduct a ,  side ):
        """
        _removeAdduct_1(self, a: Adduct , side: int ) -> Compomer
        """
        assert isinstance(a, Adduct), 'arg a wrong type'
        assert isinstance(side, int), 'arg side wrong type'
    
    
        cdef _Compomer * _r = new _Compomer(self.inst.get().removeAdduct((deref(a.inst.get())), (<unsigned int>side)))
        cdef Compomer py_result = Compomer.__new__(Compomer)
        py_result.inst = shared_ptr[_Compomer](_r)
        return py_result
    
    def removeAdduct(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeAdduct(self, a: Adduct ) -> Compomer
          :noindex:
        
        Remove ALL instances of the given adduct

        
        .. rubric:: Overload:
        .. py:function:: removeAdduct(self, a: Adduct , side: int ) -> Compomer
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], Adduct)):
            return self._removeAdduct_0(*args)
        elif (len(args)==2) and (isinstance(args[0], Adduct)) and (isinstance(args[1], int)):
            return self._removeAdduct_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getLabels(self,  side ):
        """
        getLabels(self, side: int ) -> List[bytes]
        Returns the adduct labels from parameter(side) given. (LEFT or RIGHT)
        """
        assert isinstance(side, int), 'arg side wrong type'
    
        _r = self.inst.get().getLabels((<unsigned int>side))
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result 

cdef class ConsensusIDAlgorithmBest:
    """
    Cython implementation of _ConsensusIDAlgorithmBest

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmBest.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmBest](new _ConsensusIDAlgorithmBest())
    
    def apply(self, PeptideIdentificationList ids ,  number_of_runs ):
        """
        apply(self, ids: PeptideIdentificationList , number_of_runs: int ) -> None
        Calculates the consensus ID for a set of peptide identifications of one spectrum or (consensus) feature
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(number_of_runs, int) and number_of_runs >= 0, 'arg number_of_runs wrong type'
    
    
        self.inst.get().apply((deref(ids.inst.get())), (<size_t>number_of_runs))
    
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

cdef class DateTime:
    """
    Cython implementation of _DateTime

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DateTime.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DateTime rv = DateTime.__new__(DateTime)
       rv.inst = shared_ptr[_DateTime](new _DateTime(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DateTime rv = DateTime.__new__(DateTime)
       rv.inst = shared_ptr[_DateTime](new _DateTime(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DateTime](new _DateTime())
    
    def _init_1(self, DateTime in_0 ):
        """
        _init_1(self, in_0: DateTime ) -> None
        """
        assert isinstance(in_0, DateTime), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DateTime](new _DateTime((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DateTime ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DateTime)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setDate(self,  date ):
        """
        setDate(self, date: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(date, str) or isinstance(date, bytes) or isinstance(date, String)), 'arg date wrong type'
    
        self.inst.get().setDate(deref((convString(date)).get()))
    
    def setTime(self,  date ):
        """
        setTime(self, date: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(date, str) or isinstance(date, bytes) or isinstance(date, String)), 'arg date wrong type'
    
        self.inst.get().setTime(deref((convString(date)).get()))
    
    def getDate(self):
        """
        getDate(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getDate()
        py_result = convOutputString(_r)
        return py_result
    
    def getTime(self):
        """
        getTime(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTime()
        py_result = convOutputString(_r)
        return py_result
    
    def now(self):
        """
        now(self) -> DateTime
        """
        cdef _DateTime * _r = new _DateTime(self.inst.get().now())
        cdef DateTime py_result = DateTime.__new__(DateTime)
        py_result.inst = shared_ptr[_DateTime](_r)
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def get(self):
        """
        get(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().get()
        py_result = convOutputString(_r)
        return py_result
    
    def set(self,  date ):
        """
        set(self, date: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(date, str) or isinstance(date, bytes) or isinstance(date, String)), 'arg date wrong type'
    
        self.inst.get().set(deref((convString(date)).get()))
    now = __static_DateTime_now 

cdef class ElementDB:
    """
    Cython implementation of _ElementDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElementDB.html>`_
    """

    
    def _getElement_0(self,  name ):
        """
        _getElement_0(self, name: Union[bytes, str, String] ) -> Element
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef const _Element * __r = (self.inst.get().getElement(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _Element * _r = new _Element(deref(__r))
        cdef Element py_result = Element.__new__(Element)
        py_result.inst = shared_ptr[_Element](_r)
        return py_result
    
    def _getElement_1(self,  atomic_number ):
        """
        _getElement_1(self, atomic_number: int ) -> Element
        """
        assert isinstance(atomic_number, int), 'arg atomic_number wrong type'
    
        cdef const _Element * __r = (self.inst.get().getElement((<unsigned int>atomic_number)))
        if __r == NULL:
            return None
        cdef _Element * _r = new _Element(deref(__r))
        cdef Element py_result = Element.__new__(Element)
        py_result.inst = shared_ptr[_Element](_r)
        return py_result
    
    def getElement(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getElement(self, name: Union[bytes, str, String] ) -> Element
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: getElement(self, atomic_number: int ) -> Element
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getElement_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._getElement_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def addElement(self, bytes name , bytes symbol ,  an , dict abundance , dict mass , bool replace_existing ):
        """
        addElement(self, name: bytes , symbol: bytes , an: int , abundance: Dict[int, float] , mass: Dict[int, float] , replace_existing: bool ) -> None
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(symbol, bytes), 'arg symbol wrong type'
        assert isinstance(an, int), 'arg an wrong type'
        assert isinstance(abundance, dict) and all(isinstance(k, int) for k in abundance.keys()) and all(isinstance(v, float) for v in abundance.values()), 'arg abundance wrong type'
        assert isinstance(mass, dict) and all(isinstance(k, int) for k in mass.keys()) and all(isinstance(v, float) for v in mass.values()), 'arg mass wrong type'
        assert isinstance(replace_existing, pybool_t), 'arg replace_existing wrong type'
    
    
    
        cdef libcpp_map[unsigned int, double] * v3 = new libcpp_map[unsigned int, double]()
        for key, value in abundance.items():
        
        
            deref(v3)[ (<unsigned int>key) ] = (<double>value)
        
        
        cdef libcpp_map[unsigned int, double] * v4 = new libcpp_map[unsigned int, double]()
        for key, value in mass.items():
        
        
            deref(v4)[ (<unsigned int>key) ] = (<double>value)
        
        
    
        self.inst.get().addElement((<libcpp_string>name), (<libcpp_string>symbol), (<unsigned int>an), deref(v3), deref(v4), (<bool>replace_existing))
        del v4
        del v3
    
    def _hasElement_0(self,  name ):
        """
        _hasElement_0(self, name: Union[bytes, str, String] ) -> bool
        Returns true if the db contains an element with the given name, else false
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasElement(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    def _hasElement_1(self,  atomic_number ):
        """
        _hasElement_1(self, atomic_number: int ) -> bool
        Returns true if the db contains an element with the given atomic_number, else false
        """
        assert isinstance(atomic_number, int), 'arg atomic_number wrong type'
    
        cdef bool _r = self.inst.get().hasElement((<unsigned int>atomic_number))
        py_result = <bool>_r
        return py_result
    
    def hasElement(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: hasElement(self, name: Union[bytes, str, String] ) -> bool
          :noindex:
        
        Returns true if the db contains an element with the given name, else false

        
        .. rubric:: Overload:
        .. py:function:: hasElement(self, atomic_number: int ) -> bool
          :noindex:
        
        Returns true if the db contains an element with the given atomic_number, else false
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._hasElement_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._hasElement_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_ElementDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_ElementDB](_getInstance_ElementDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class IsotopeFitter1D:
    """
    Cython implementation of _IsotopeFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeFitter1D.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsotopeFitter1D rv = IsotopeFitter1D.__new__(IsotopeFitter1D)
       rv.inst = shared_ptr[_IsotopeFitter1D](new _IsotopeFitter1D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopeFitter1D rv = IsotopeFitter1D.__new__(IsotopeFitter1D)
       rv.inst = shared_ptr[_IsotopeFitter1D](new _IsotopeFitter1D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Isotope distribution fitter (1-dim.) approximated using linear interpolation
        """
        self.inst = shared_ptr[_IsotopeFitter1D](new _IsotopeFitter1D())
    
    def _init_1(self, IsotopeFitter1D in_0 ):
        """
        _init_1(self, in_0: IsotopeFitter1D ) -> None
        """
        assert isinstance(in_0, IsotopeFitter1D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopeFitter1D](new _IsotopeFitter1D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Isotope distribution fitter (1-dim.) approximated using linear interpolation

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopeFitter1D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopeFitter1D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MetaInfoRegistry:
    """
    Cython implementation of _MetaInfoRegistry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoRegistry.html>`_

    Registry which assigns unique integer indices to strings
    
    When registering a new name an index >= 1024 is assigned.
    Indices from 1 to 1023 are reserved for fast access and will never change:
    1 - isotopic_range
    2 - cluster_id
    3 - label
    4 - icon
    5 - color
    6 - RT
    7 - MZ
    8 - predicted_RT
    9 - predicted_RT_p_value
    10 - spectrum_reference
    11 - ID
    12 - low_quality
    13 - charge
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaInfoRegistry rv = MetaInfoRegistry.__new__(MetaInfoRegistry)
       rv.inst = shared_ptr[_MetaInfoRegistry](new _MetaInfoRegistry(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaInfoRegistry rv = MetaInfoRegistry.__new__(MetaInfoRegistry)
       rv.inst = shared_ptr[_MetaInfoRegistry](new _MetaInfoRegistry(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaInfoRegistry](new _MetaInfoRegistry())
    
    def _init_1(self, MetaInfoRegistry in_0 ):
        """
        _init_1(self, in_0: MetaInfoRegistry ) -> None
        """
        assert isinstance(in_0, MetaInfoRegistry), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaInfoRegistry](new _MetaInfoRegistry((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaInfoRegistry ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaInfoRegistry)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def registerName(self,  name ,  description ,  unit ):
        """
        registerName(self, name: Union[bytes, str, String] , description: Union[bytes, str, String] , unit: Union[bytes, str, String] ) -> int
        Registers a string, stores its description and unit, and returns the corresponding index. If the string is already registered, it returns the index of the string
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
        assert (isinstance(unit, str) or isinstance(unit, bytes) or isinstance(unit, String)), 'arg unit wrong type'
    
    
    
        cdef unsigned int _r = self.inst.get().registerName(deref((convString(name)).get()), deref((convString(description)).get()), deref((convString(unit)).get()))
        py_result = <unsigned int>_r
        return py_result
    
    def _setDescription_0(self,  index ,  description ):
        """
        _setDescription_0(self, index: int , description: Union[bytes, str, String] ) -> None
        Sets the description (String), corresponding to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
    
    
        self.inst.get().setDescription((<unsigned int>index), deref((convString(description)).get()))
    
    def _setDescription_1(self,  name ,  description ):
        """
        _setDescription_1(self, name: Union[bytes, str, String] , description: Union[bytes, str, String] ) -> None
        Sets the description (String), corresponding to a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
    
    
        self.inst.get().setDescription(deref((convString(name)).get()), deref((convString(description)).get()))
    
    def setDescription(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setDescription(self, index: int , description: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the description (String), corresponding to an index

        
        .. rubric:: Overload:
        .. py:function:: setDescription(self, name: Union[bytes, str, String] , description: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the description (String), corresponding to a name
    
        """
        if (len(args)==2) and (isinstance(args[0], int)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._setDescription_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._setDescription_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setUnit_0(self,  index ,  unit ):
        """
        _setUnit_0(self, index: int , unit: Union[bytes, str, String] ) -> None
        Sets the unit (String), corresponding to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
        assert (isinstance(unit, str) or isinstance(unit, bytes) or isinstance(unit, String)), 'arg unit wrong type'
    
    
        self.inst.get().setUnit((<unsigned int>index), deref((convString(unit)).get()))
    
    def _setUnit_1(self,  name ,  unit ):
        """
        _setUnit_1(self, name: Union[bytes, str, String] , unit: Union[bytes, str, String] ) -> None
        Sets the unit (String), corresponding to a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(unit, str) or isinstance(unit, bytes) or isinstance(unit, String)), 'arg unit wrong type'
    
    
        self.inst.get().setUnit(deref((convString(name)).get()), deref((convString(unit)).get()))
    
    def setUnit(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setUnit(self, index: int , unit: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the unit (String), corresponding to an index

        
        .. rubric:: Overload:
        .. py:function:: setUnit(self, name: Union[bytes, str, String] , unit: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the unit (String), corresponding to a name
    
        """
        if (len(args)==2) and (isinstance(args[0], int)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._setUnit_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._setUnit_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIndex(self,  name ):
        """
        getIndex(self, name: Union[bytes, str, String] ) -> int
        Returns the integer index corresponding to a string. If the string is not registered, returns UInt(-1) (= UINT_MAX)
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef unsigned int _r = self.inst.get().getIndex(deref((convString(name)).get()))
        py_result = <unsigned int>_r
        return py_result
    
    def getName(self,  index ):
        """
        getName(self, index: int ) -> Union[bytes, str, String]
        Returns the corresponding name to an index
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _String _r = self.inst.get().getName((<unsigned int>index))
        py_result = convOutputString(_r)
        return py_result
    
    def _getDescription_0(self,  index ):
        """
        _getDescription_0(self, index: int ) -> Union[bytes, str, String]
        Returns the description of an index
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _String _r = self.inst.get().getDescription((<unsigned int>index))
        py_result = convOutputString(_r)
        return py_result
    
    def _getDescription_1(self,  name ):
        """
        _getDescription_1(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Returns the description of a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _String _r = self.inst.get().getDescription(deref((convString(name)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def getDescription(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getDescription(self, index: int ) -> Union[bytes, str, String]
          :noindex:
        
        Returns the description of an index

        
        .. rubric:: Overload:
        .. py:function:: getDescription(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]
          :noindex:
        
        Returns the description of a name
    
        """
        if (len(args)==1) and (isinstance(args[0], int)):
            return self._getDescription_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getDescription_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getUnit_0(self,  index ):
        """
        _getUnit_0(self, index: int ) -> Union[bytes, str, String]
        Returns the unit of an index
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _String _r = self.inst.get().getUnit((<unsigned int>index))
        py_result = convOutputString(_r)
        return py_result
    
    def _getUnit_1(self,  name ):
        """
        _getUnit_1(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Returns the unit of a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _String _r = self.inst.get().getUnit(deref((convString(name)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def getUnit(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getUnit(self, index: int ) -> Union[bytes, str, String]
          :noindex:
        
        Returns the unit of an index

        
        .. rubric:: Overload:
        .. py:function:: getUnit(self, name: Union[bytes, str, String] ) -> Union[bytes, str, String]
          :noindex:
        
        Returns the unit of a name
    
        """
        if (len(args)==1) and (isinstance(args[0], int)):
            return self._getUnit_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getUnit_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OPXLHelper:
    """
    Cython implementation of _OPXLHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLHelper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OPXLHelper rv = OPXLHelper.__new__(OPXLHelper)
       rv.inst = shared_ptr[_OPXLHelper](new _OPXLHelper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OPXLHelper rv = OPXLHelper.__new__(OPXLHelper)
       rv.inst = shared_ptr[_OPXLHelper](new _OPXLHelper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OPXLHelper](new _OPXLHelper())
    
    def _init_1(self, OPXLHelper in_0 ):
        """
        _init_1(self, in_0: OPXLHelper ) -> None
        """
        assert isinstance(in_0, OPXLHelper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OPXLHelper](new _OPXLHelper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OPXLHelper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OPXLHelper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def enumerateCrossLinksAndMasses(self, list peptides , double cross_link_mass_light , list cross_link_mass_mono_link , list cross_link_residue1 , list cross_link_residue2 , list spectrum_precursors , list precursor_correction_positions , double precursor_mass_tolerance , bool precursor_mass_tolerance_unit_ppm ):
        """
        enumerateCrossLinksAndMasses(self, peptides: List[AASeqWithMass] , cross_link_mass_light: float , cross_link_mass_mono_link: List[float] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , spectrum_precursors: List[float] , precursor_correction_positions: List[int] , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool ) -> List[XLPrecursor]
        """
        assert isinstance(peptides, list) and all(isinstance(elemt_rec, AASeqWithMass) for elemt_rec in peptides), 'arg peptides wrong type'
        assert isinstance(cross_link_mass_light, float), 'arg cross_link_mass_light wrong type'
        assert isinstance(cross_link_mass_mono_link, list) and all(isinstance(li, float) for li in cross_link_mass_mono_link), 'arg cross_link_mass_mono_link wrong type'
        assert isinstance(cross_link_residue1, list) and all(isinstance(li, bytes) for li in cross_link_residue1), 'arg cross_link_residue1 wrong type'
        assert isinstance(cross_link_residue2, list) and all(isinstance(li, bytes) for li in cross_link_residue2), 'arg cross_link_residue2 wrong type'
        assert isinstance(spectrum_precursors, list) and all(isinstance(elemt_rec, float) for elemt_rec in spectrum_precursors), 'arg spectrum_precursors wrong type'
        assert isinstance(precursor_correction_positions, list) and all(isinstance(elemt_rec, int) for elemt_rec in precursor_correction_positions), 'arg precursor_correction_positions wrong type'
        assert isinstance(precursor_mass_tolerance, float), 'arg precursor_mass_tolerance wrong type'
        assert isinstance(precursor_mass_tolerance_unit_ppm, pybool_t), 'arg precursor_mass_tolerance_unit_ppm wrong type'
        cdef libcpp_vector[_AASeqWithMass] * v0 = new libcpp_vector[_AASeqWithMass]()
        cdef AASeqWithMass item0
        for item0 in peptides:
            v0.push_back(deref(item0.inst.get()))
    
        cdef libcpp_vector[double] _v2 = cross_link_mass_mono_link
        cdef _DoubleList v2 = _DoubleList(_v2)
        cdef libcpp_vector[_String] * v3 = new libcpp_vector[_String]()
        cdef bytes item3
        for item3 in cross_link_residue1:
           v3.push_back(_String(<char *>item3))
        cdef libcpp_vector[_String] * v4 = new libcpp_vector[_String]()
        cdef bytes item4
        for item4 in cross_link_residue2:
           v4.push_back(_String(<char *>item4))
        cdef libcpp_vector[double] v5 = spectrum_precursors
        cdef libcpp_vector[int] v6 = precursor_correction_positions
    
    
        _r = self.inst.get().enumerateCrossLinksAndMasses(deref(v0), (<double>cross_link_mass_light), (v2), deref(v3), deref(v4), v5, v6, (<double>precursor_mass_tolerance), (<bool>precursor_mass_tolerance_unit_ppm))
        precursor_correction_positions[:] = v6
        spectrum_precursors[:] = v5
        del v4
        del v3
        del v0
        py_result = []
        cdef libcpp_vector[_XLPrecursor].iterator it__r = _r.begin()
        cdef XLPrecursor item_py_result
        while it__r != _r.end():
           item_py_result = XLPrecursor.__new__(XLPrecursor)
           item_py_result.inst = shared_ptr[_XLPrecursor](new _XLPrecursor(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def digestDatabase(self, list fasta_db , EnzymaticDigestion digestor ,  min_peptide_length , list cross_link_residue1 , list cross_link_residue2 , ModifiedPeptideGenerator_MapToResidueType fixed_modifications , ModifiedPeptideGenerator_MapToResidueType variable_modifications ,  max_variable_mods_per_peptide ):
        """
        digestDatabase(self, fasta_db: List[FASTAEntry] , digestor: EnzymaticDigestion , min_peptide_length: int , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , fixed_modifications: ModifiedPeptideGenerator_MapToResidueType , variable_modifications: ModifiedPeptideGenerator_MapToResidueType , max_variable_mods_per_peptide: int ) -> List[AASeqWithMass]
        """
        assert isinstance(fasta_db, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in fasta_db), 'arg fasta_db wrong type'
        assert isinstance(digestor, EnzymaticDigestion), 'arg digestor wrong type'
        assert isinstance(min_peptide_length, int) and min_peptide_length >= 0, 'arg min_peptide_length wrong type'
        assert isinstance(cross_link_residue1, list) and all(isinstance(li, bytes) for li in cross_link_residue1), 'arg cross_link_residue1 wrong type'
        assert isinstance(cross_link_residue2, list) and all(isinstance(li, bytes) for li in cross_link_residue2), 'arg cross_link_residue2 wrong type'
        assert isinstance(fixed_modifications, ModifiedPeptideGenerator_MapToResidueType), 'arg fixed_modifications wrong type'
        assert isinstance(variable_modifications, ModifiedPeptideGenerator_MapToResidueType), 'arg variable_modifications wrong type'
        assert isinstance(max_variable_mods_per_peptide, int) and max_variable_mods_per_peptide >= 0, 'arg max_variable_mods_per_peptide wrong type'
        cdef libcpp_vector[_FASTAEntry] * v0 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item0
        for item0 in fasta_db:
            v0.push_back(deref(item0.inst.get()))
    
    
        cdef libcpp_vector[_String] * v3 = new libcpp_vector[_String]()
        cdef bytes item3
        for item3 in cross_link_residue1:
           v3.push_back(_String(<char *>item3))
        cdef libcpp_vector[_String] * v4 = new libcpp_vector[_String]()
        cdef bytes item4
        for item4 in cross_link_residue2:
           v4.push_back(_String(<char *>item4))
    
    
    
        _r = self.inst.get().digestDatabase(deref(v0), (deref(digestor.inst.get())), (<size_t>min_peptide_length), deref(v3), deref(v4), (deref(fixed_modifications.inst.get())), (deref(variable_modifications.inst.get())), (<size_t>max_variable_mods_per_peptide))
        del v4
        del v3
        del v0
        py_result = []
        cdef libcpp_vector[_AASeqWithMass].iterator it__r = _r.begin()
        cdef AASeqWithMass item_py_result
        while it__r != _r.end():
           item_py_result = AASeqWithMass.__new__(AASeqWithMass)
           item_py_result.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def buildCandidates(self, list candidates , list precursor_corrections , list precursor_correction_positions , list peptide_masses , list cross_link_residue1 , list cross_link_residue2 , double cross_link_mass , list cross_link_mass_mono_link , list spectrum_precursor_vector , list allowed_error_vector ,  cross_link_name ):
        """
        buildCandidates(self, candidates: List[XLPrecursor] , precursor_corrections: List[int] , precursor_correction_positions: List[int] , peptide_masses: List[AASeqWithMass] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , cross_link_mass: float , cross_link_mass_mono_link: List[float] , spectrum_precursor_vector: List[float] , allowed_error_vector: List[float] , cross_link_name: Union[bytes, str, String] ) -> List[ProteinProteinCrossLink]
        """
        assert isinstance(candidates, list) and all(isinstance(elemt_rec, XLPrecursor) for elemt_rec in candidates), 'arg candidates wrong type'
        assert isinstance(precursor_corrections, list) and all(isinstance(elemt_rec, int) for elemt_rec in precursor_corrections), 'arg precursor_corrections wrong type'
        assert isinstance(precursor_correction_positions, list) and all(isinstance(elemt_rec, int) for elemt_rec in precursor_correction_positions), 'arg precursor_correction_positions wrong type'
        assert isinstance(peptide_masses, list) and all(isinstance(elemt_rec, AASeqWithMass) for elemt_rec in peptide_masses), 'arg peptide_masses wrong type'
        assert isinstance(cross_link_residue1, list) and all(isinstance(li, bytes) for li in cross_link_residue1), 'arg cross_link_residue1 wrong type'
        assert isinstance(cross_link_residue2, list) and all(isinstance(li, bytes) for li in cross_link_residue2), 'arg cross_link_residue2 wrong type'
        assert isinstance(cross_link_mass, float), 'arg cross_link_mass wrong type'
        assert isinstance(cross_link_mass_mono_link, list) and all(isinstance(li, float) for li in cross_link_mass_mono_link), 'arg cross_link_mass_mono_link wrong type'
        assert isinstance(spectrum_precursor_vector, list) and all(isinstance(elemt_rec, float) for elemt_rec in spectrum_precursor_vector), 'arg spectrum_precursor_vector wrong type'
        assert isinstance(allowed_error_vector, list) and all(isinstance(elemt_rec, float) for elemt_rec in allowed_error_vector), 'arg allowed_error_vector wrong type'
        assert (isinstance(cross_link_name, str) or isinstance(cross_link_name, bytes) or isinstance(cross_link_name, String)), 'arg cross_link_name wrong type'
        cdef libcpp_vector[_XLPrecursor] * v0 = new libcpp_vector[_XLPrecursor]()
        cdef XLPrecursor item0
        for item0 in candidates:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[int] v1 = precursor_corrections
        cdef libcpp_vector[int] v2 = precursor_correction_positions
        cdef libcpp_vector[_AASeqWithMass] * v3 = new libcpp_vector[_AASeqWithMass]()
        cdef AASeqWithMass item3
        for item3 in peptide_masses:
            v3.push_back(deref(item3.inst.get()))
        cdef libcpp_vector[_String] * v4 = new libcpp_vector[_String]()
        cdef bytes item4
        for item4 in cross_link_residue1:
           v4.push_back(_String(<char *>item4))
        cdef libcpp_vector[_String] * v5 = new libcpp_vector[_String]()
        cdef bytes item5
        for item5 in cross_link_residue2:
           v5.push_back(_String(<char *>item5))
    
        cdef libcpp_vector[double] _v7 = cross_link_mass_mono_link
        cdef _DoubleList v7 = _DoubleList(_v7)
        cdef libcpp_vector[double] v8 = spectrum_precursor_vector
        cdef libcpp_vector[double] v9 = allowed_error_vector
    
        _r = self.inst.get().buildCandidates(deref(v0), v1, v2, deref(v3), deref(v4), deref(v5), (<double>cross_link_mass), (v7), v8, v9, deref((convString(cross_link_name)).get()))
        allowed_error_vector[:] = v9
        spectrum_precursor_vector[:] = v8
        replace = []
        cdef libcpp_vector[_String].iterator it_5 = v5.begin()
        while it_5 != v5.end():
           replace.append(<char*>deref(it_5).c_str())
           inc(it_5)
        cross_link_residue2[:] = replace
        del v5
        replace = []
        cdef libcpp_vector[_String].iterator it_4 = v4.begin()
        while it_4 != v4.end():
           replace.append(<char*>deref(it_4).c_str())
           inc(it_4)
        cross_link_residue1[:] = replace
        del v4
        cdef libcpp_vector[_AASeqWithMass].iterator it_peptide_masses = v3.begin()
        replace_0 = []
        while it_peptide_masses != v3.end():
            item3 = AASeqWithMass.__new__(AASeqWithMass)
            item3.inst = shared_ptr[_AASeqWithMass](new _AASeqWithMass(deref(it_peptide_masses)))
            replace_0.append(item3)
            inc(it_peptide_masses)
        peptide_masses[:] = replace_0
        del v3
        precursor_correction_positions[:] = v2
        precursor_corrections[:] = v1
        cdef libcpp_vector[_XLPrecursor].iterator it_candidates = v0.begin()
        replace_0 = []
        while it_candidates != v0.end():
            item0 = XLPrecursor.__new__(XLPrecursor)
            item0.inst = shared_ptr[_XLPrecursor](new _XLPrecursor(deref(it_candidates)))
            replace_0.append(item0)
            inc(it_candidates)
        candidates[:] = replace_0
        del v0
        py_result = []
        cdef libcpp_vector[_ProteinProteinCrossLink].iterator it__r = _r.begin()
        cdef ProteinProteinCrossLink item_py_result
        while it__r != _r.end():
           item_py_result = ProteinProteinCrossLink.__new__(ProteinProteinCrossLink)
           item_py_result.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def buildFragmentAnnotations(self, list frag_annotations , list matching , MSSpectrum theoretical_spectrum , MSSpectrum experiment_spectrum ):
        """
        buildFragmentAnnotations(self, frag_annotations: List[PeptideHit_PeakAnnotation] , matching: List[List[int, int]] , theoretical_spectrum: MSSpectrum , experiment_spectrum: MSSpectrum ) -> None
        """
        assert isinstance(frag_annotations, list) and all(isinstance(elemt_rec, PeptideHit_PeakAnnotation) for elemt_rec in frag_annotations), 'arg frag_annotations wrong type'
        assert isinstance(matching, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matching), 'arg matching wrong type'
        assert isinstance(theoretical_spectrum, MSSpectrum), 'arg theoretical_spectrum wrong type'
        assert isinstance(experiment_spectrum, MSSpectrum), 'arg experiment_spectrum wrong type'
        cdef libcpp_vector[_PeptideHit_PeakAnnotation] * v0 = new libcpp_vector[_PeptideHit_PeakAnnotation]()
        cdef PeptideHit_PeakAnnotation item0
        for item0 in frag_annotations:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v1 = matching
    
    
        self.inst.get().buildFragmentAnnotations(deref(v0), v1, (deref(theoretical_spectrum.inst.get())), (deref(experiment_spectrum.inst.get())))
        
        cdef libcpp_vector[_PeptideHit_PeakAnnotation].iterator it_frag_annotations = v0.begin()
        replace_0 = []
        while it_frag_annotations != v0.end():
            item0 = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
            item0.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(it_frag_annotations)))
            replace_0.append(item0)
            inc(it_frag_annotations)
        frag_annotations[:] = replace_0
        del v0
    
    def buildPeptideIDs(self, PeptideIdentificationList peptide_ids , list top_csms_spectrum , list all_top_csms ,  all_top_csms_current_index , MSExperiment spectra ,  scan_index ,  scan_index_heavy ):
        """
        buildPeptideIDs(self, peptide_ids: PeptideIdentificationList , top_csms_spectrum: List[CrossLinkSpectrumMatch] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , all_top_csms_current_index: int , spectra: MSExperiment , scan_index: int , scan_index_heavy: int ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(top_csms_spectrum, list) and all(isinstance(elemt_rec, CrossLinkSpectrumMatch) for elemt_rec in top_csms_spectrum), 'arg top_csms_spectrum wrong type'
        assert isinstance(all_top_csms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in all_top_csms), 'arg all_top_csms wrong type'
        assert isinstance(all_top_csms_current_index, int) and all_top_csms_current_index >= 0, 'arg all_top_csms_current_index wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(scan_index, int) and scan_index >= 0, 'arg scan_index wrong type'
        assert isinstance(scan_index_heavy, int) and scan_index_heavy >= 0, 'arg scan_index_heavy wrong type'
    
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v1 = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item1
        for item1 in top_csms_spectrum:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]] * v2 = new libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]]()
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v2_rec = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item2_rec
        cdef libcpp_vector[_CrossLinkSpectrumMatch].iterator it_all_top_csms_rec
        for all_top_csms_rec in all_top_csms:
            v2_rec.clear()
            for item2_rec in all_top_csms_rec:
                v2_rec.push_back(deref(item2_rec.inst.get()))
            v2.push_back(deref(v2_rec))
    
    
    
    
        self.inst.get().buildPeptideIDs((deref(peptide_ids.inst.get())), deref(v1), deref(v2), (<size_t>all_top_csms_current_index), (deref(spectra.inst.get())), (<size_t>scan_index), (<size_t>scan_index_heavy))
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]].iterator it_all_top_csms = v2.begin()
        replace_0 = []
        while it_all_top_csms != v2.end():
            it_all_top_csms_rec = deref(it_all_top_csms).begin()
            replace_1 = []
            while it_all_top_csms_rec != deref(it_all_top_csms).end():
                item2_rec = CrossLinkSpectrumMatch.__new__(CrossLinkSpectrumMatch)
                item2_rec.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch(deref(it_all_top_csms_rec)))
                replace_1.append(item2_rec)
                inc(it_all_top_csms_rec)
            replace_0.append(replace_1)
            inc(it_all_top_csms)
        all_top_csms[:] = replace_0
        del v2
        del v1
    
    def addProteinPositionMetaValues(self, PeptideIdentificationList peptide_ids ):
        """
        addProteinPositionMetaValues(self, peptide_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        self.inst.get().addProteinPositionMetaValues((deref(peptide_ids.inst.get())))
    
    def addXLTargetDecoyMV(self, PeptideIdentificationList peptide_ids ):
        """
        addXLTargetDecoyMV(self, peptide_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        self.inst.get().addXLTargetDecoyMV((deref(peptide_ids.inst.get())))
    
    def addBetaAccessions(self, PeptideIdentificationList peptide_ids ):
        """
        addBetaAccessions(self, peptide_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        self.inst.get().addBetaAccessions((deref(peptide_ids.inst.get())))
    
    def removeBetaPeptideHits(self, PeptideIdentificationList peptide_ids ):
        """
        removeBetaPeptideHits(self, peptide_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        self.inst.get().removeBetaPeptideHits((deref(peptide_ids.inst.get())))
    
    def addPercolatorFeatureList(self, ProteinIdentification prot_id ):
        """
        addPercolatorFeatureList(self, prot_id: ProteinIdentification ) -> None
        """
        assert isinstance(prot_id, ProteinIdentification), 'arg prot_id wrong type'
    
        self.inst.get().addPercolatorFeatureList((deref(prot_id.inst.get())))
    
    def computeDeltaScores(self, PeptideIdentificationList peptide_ids ):
        """
        computeDeltaScores(self, peptide_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        self.inst.get().computeDeltaScores((deref(peptide_ids.inst.get())))
    
    def combineTopRanksFromPairs(self, PeptideIdentificationList peptide_ids ,  number_top_hits ):
        """
        combineTopRanksFromPairs(self, peptide_ids: PeptideIdentificationList , number_top_hits: int ) -> PeptideIdentificationList
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(number_top_hits, int) and number_top_hits >= 0, 'arg number_top_hits wrong type'
    
    
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().combineTopRanksFromPairs((deref(peptide_ids.inst.get())), (<size_t>number_top_hits)))
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result
    
    def collectPrecursorCandidates(self, list precursor_correction_steps , double precursor_mass , double precursor_mass_tolerance , bool precursor_mass_tolerance_unit_ppm , list filtered_peptide_masses , double cross_link_mass , list cross_link_mass_mono_link , list cross_link_residue1 , list cross_link_residue2 ,  cross_link_name , bool use_sequence_tags , list tags ):
        """
        collectPrecursorCandidates(self, precursor_correction_steps: List[int] , precursor_mass: float , precursor_mass_tolerance: float , precursor_mass_tolerance_unit_ppm: bool , filtered_peptide_masses: List[AASeqWithMass] , cross_link_mass: float , cross_link_mass_mono_link: List[float] , cross_link_residue1: List[bytes] , cross_link_residue2: List[bytes] , cross_link_name: Union[bytes, str, String] , use_sequence_tags: bool , tags: List[Union[bytes, str]] ) -> List[ProteinProteinCrossLink]
        """
        assert isinstance(precursor_correction_steps, list) and all(isinstance(li, int) for li in precursor_correction_steps), 'arg precursor_correction_steps wrong type'
        assert isinstance(precursor_mass, float), 'arg precursor_mass wrong type'
        assert isinstance(precursor_mass_tolerance, float), 'arg precursor_mass_tolerance wrong type'
        assert isinstance(precursor_mass_tolerance_unit_ppm, pybool_t), 'arg precursor_mass_tolerance_unit_ppm wrong type'
        assert isinstance(filtered_peptide_masses, list) and all(isinstance(elemt_rec, AASeqWithMass) for elemt_rec in filtered_peptide_masses), 'arg filtered_peptide_masses wrong type'
        assert isinstance(cross_link_mass, float), 'arg cross_link_mass wrong type'
        assert isinstance(cross_link_mass_mono_link, list) and all(isinstance(li, float) for li in cross_link_mass_mono_link), 'arg cross_link_mass_mono_link wrong type'
        assert isinstance(cross_link_residue1, list) and all(isinstance(li, bytes) for li in cross_link_residue1), 'arg cross_link_residue1 wrong type'
        assert isinstance(cross_link_residue2, list) and all(isinstance(li, bytes) for li in cross_link_residue2), 'arg cross_link_residue2 wrong type'
        assert (isinstance(cross_link_name, str) or isinstance(cross_link_name, bytes) or isinstance(cross_link_name, String)), 'arg cross_link_name wrong type'
        assert isinstance(use_sequence_tags, pybool_t), 'arg use_sequence_tags wrong type'
        assert isinstance(tags, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tags), 'arg tags wrong type'
        cdef libcpp_vector[int] _v0 = precursor_correction_steps
        cdef _IntList v0 = _IntList(_v0)
    
    
    
        cdef libcpp_vector[_AASeqWithMass] * v4 = new libcpp_vector[_AASeqWithMass]()
        cdef AASeqWithMass item4
        for item4 in filtered_peptide_masses:
            v4.push_back(deref(item4.inst.get()))
    
        cdef libcpp_vector[double] _v6 = cross_link_mass_mono_link
        cdef _DoubleList v6 = _DoubleList(_v6)
        cdef libcpp_vector[_String] * v7 = new libcpp_vector[_String]()
        cdef bytes item7
        for item7 in cross_link_residue1:
           v7.push_back(_String(<char *>item7))
        cdef libcpp_vector[_String] * v8 = new libcpp_vector[_String]()
        cdef bytes item8
        for item8 in cross_link_residue2:
           v8.push_back(_String(<char *>item8))
    
    
        cdef libcpp_vector[libcpp_utf8_string] v11 = tags
        _r = self.inst.get().collectPrecursorCandidates((v0), (<double>precursor_mass), (<double>precursor_mass_tolerance), (<bool>precursor_mass_tolerance_unit_ppm), deref(v4), (<double>cross_link_mass), (v6), deref(v7), deref(v8), deref((convString(cross_link_name)).get()), (<bool>use_sequence_tags), v11)
        
        del v8
        del v7
        del v4
        py_result = []
        cdef libcpp_vector[_ProteinProteinCrossLink].iterator it__r = _r.begin()
        cdef ProteinProteinCrossLink item_py_result
        while it__r != _r.end():
           item_py_result = ProteinProteinCrossLink.__new__(ProteinProteinCrossLink)
           item_py_result.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def computePrecursorError(self, CrossLinkSpectrumMatch csm , double precursor_mz ,  precursor_charge ):
        """
        computePrecursorError(self, csm: CrossLinkSpectrumMatch , precursor_mz: float , precursor_charge: int ) -> float
        """
        assert isinstance(csm, CrossLinkSpectrumMatch), 'arg csm wrong type'
        assert isinstance(precursor_mz, float), 'arg precursor_mz wrong type'
        assert isinstance(precursor_charge, int), 'arg precursor_charge wrong type'
    
    
    
        cdef double _r = self.inst.get().computePrecursorError((deref(csm.inst.get())), (<double>precursor_mz), (<int>precursor_charge))
        py_result = <double>_r
        return py_result
    
    def isoPeakMeans(self, CrossLinkSpectrumMatch csm , IntegerDataArray num_iso_peaks_array , list matched_spec_linear_alpha , list matched_spec_linear_beta , list matched_spec_xlinks_alpha , list matched_spec_xlinks_beta ):
        """
        isoPeakMeans(self, csm: CrossLinkSpectrumMatch , num_iso_peaks_array: IntegerDataArray , matched_spec_linear_alpha: List[List[int, int]] , matched_spec_linear_beta: List[List[int, int]] , matched_spec_xlinks_alpha: List[List[int, int]] , matched_spec_xlinks_beta: List[List[int, int]] ) -> None
        """
        assert isinstance(csm, CrossLinkSpectrumMatch), 'arg csm wrong type'
        assert isinstance(num_iso_peaks_array, IntegerDataArray), 'arg num_iso_peaks_array wrong type'
        assert isinstance(matched_spec_linear_alpha, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_linear_alpha), 'arg matched_spec_linear_alpha wrong type'
        assert isinstance(matched_spec_linear_beta, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_linear_beta), 'arg matched_spec_linear_beta wrong type'
        assert isinstance(matched_spec_xlinks_alpha, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_xlinks_alpha), 'arg matched_spec_xlinks_alpha wrong type'
        assert isinstance(matched_spec_xlinks_beta, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_xlinks_beta), 'arg matched_spec_xlinks_beta wrong type'
    
    
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v2 = matched_spec_linear_alpha
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v3 = matched_spec_linear_beta
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v4 = matched_spec_xlinks_alpha
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v5 = matched_spec_xlinks_beta
        self.inst.get().isoPeakMeans((deref(csm.inst.get())), (deref(num_iso_peaks_array.inst.get())), v2, v3, v4, v5)
        matched_spec_xlinks_beta[:] = v5
        matched_spec_xlinks_alpha[:] = v4
        matched_spec_linear_beta[:] = v3
        matched_spec_linear_alpha[:] = v2 

cdef class PepXMLFileMascot:
    """
    Cython implementation of _PepXMLFileMascot

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PepXMLFileMascot.html>`_

    Used to load Mascot PepXML files
    
    A schema for this format can be found at http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_PepXMLFileMascot](new _PepXMLFileMascot()) 

cdef class PosteriorErrorProbabilityModel:
    """
    Cython implementation of _PosteriorErrorProbabilityModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1PosteriorErrorProbabilityModel.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_PosteriorErrorProbabilityModel](new _PosteriorErrorProbabilityModel())
    
    def _fit_0(self, list search_engine_scores ,  outlier_handling ):
        """
        _fit_0(self, search_engine_scores: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
        """
        assert isinstance(search_engine_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in search_engine_scores), 'arg search_engine_scores wrong type'
        assert (isinstance(outlier_handling, str) or isinstance(outlier_handling, bytes) or isinstance(outlier_handling, String)), 'arg outlier_handling wrong type'
        cdef libcpp_vector[double] v0 = search_engine_scores
    
        cdef bool _r = self.inst.get().fit(v0, deref((convString(outlier_handling)).get()))
        search_engine_scores[:] = v0
        py_result = <bool>_r
        return py_result
    
    def _fit_1(self, list search_engine_scores , list probabilities ,  outlier_handling ):
        """
        _fit_1(self, search_engine_scores: List[float] , probabilities: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :param probabilities: A vector which holds the probability for each data point after running this function. If it has some content it will be overwritten
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
        """
        assert isinstance(search_engine_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in search_engine_scores), 'arg search_engine_scores wrong type'
        assert isinstance(probabilities, list) and all(isinstance(elemt_rec, float) for elemt_rec in probabilities), 'arg probabilities wrong type'
        assert (isinstance(outlier_handling, str) or isinstance(outlier_handling, bytes) or isinstance(outlier_handling, String)), 'arg outlier_handling wrong type'
        cdef libcpp_vector[double] v0 = search_engine_scores
        cdef libcpp_vector[double] v1 = probabilities
    
        cdef bool _r = self.inst.get().fit(v0, v1, deref((convString(outlier_handling)).get()))
        probabilities[:] = v1
        search_engine_scores[:] = v0
        py_result = <bool>_r
        return py_result
    
    def fit(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: fit(self, search_engine_scores: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool
          :noindex:
        
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
        
        .. rubric:: Overload:
        .. py:function:: fit(self, search_engine_scores: List[float] , probabilities: List[float] , outlier_handling: Union[bytes, str, String] ) -> bool
          :noindex:
        
        Fits the distributions to the data points(search_engine_scores). Estimated parameters for the distributions are saved in member variables
        computeProbability can be used afterwards
        Uses two Gaussians to fit. And Gauss+Gauss or Gumbel+Gauss to plot and calculate final probabilities
        
        
        :param search_engine_scores: A vector which holds the data points
        :param probabilities: A vector which holds the probability for each data point after running this function. If it has some content it will be overwritten
        :return: `true` if algorithm has run through. Else false will be returned. In that case no plot and no probabilities are calculated
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._fit_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
            return self._fit_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def fillDensities(self, list x_scores , list incorrect_density , list correct_density ):
        """
        fillDensities(self, x_scores: List[float] , incorrect_density: List[float] , correct_density: List[float] ) -> None
        Writes the distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences
        """
        assert isinstance(x_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in x_scores), 'arg x_scores wrong type'
        assert isinstance(incorrect_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in incorrect_density), 'arg incorrect_density wrong type'
        assert isinstance(correct_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in correct_density), 'arg correct_density wrong type'
        cdef libcpp_vector[double] v0 = x_scores
        cdef libcpp_vector[double] v1 = incorrect_density
        cdef libcpp_vector[double] v2 = correct_density
        self.inst.get().fillDensities(v0, v1, v2)
        correct_density[:] = v2
        incorrect_density[:] = v1
        x_scores[:] = v0
    
    def fillLogDensities(self, list x_scores , list incorrect_density , list correct_density ):
        """
        fillLogDensities(self, x_scores: List[float] , incorrect_density: List[float] , correct_density: List[float] ) -> None
        Writes the log distributions densities into the two vectors for a set of scores. Incorrect_densities represent the incorrectly assigned sequences
        """
        assert isinstance(x_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in x_scores), 'arg x_scores wrong type'
        assert isinstance(incorrect_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in incorrect_density), 'arg incorrect_density wrong type'
        assert isinstance(correct_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in correct_density), 'arg correct_density wrong type'
        cdef libcpp_vector[double] v0 = x_scores
        cdef libcpp_vector[double] v1 = incorrect_density
        cdef libcpp_vector[double] v2 = correct_density
        self.inst.get().fillLogDensities(v0, v1, v2)
        correct_density[:] = v2
        incorrect_density[:] = v1
        x_scores[:] = v0
    
    def computeLogLikelihood(self, list incorrect_density , list correct_density ):
        """
        computeLogLikelihood(self, incorrect_density: List[float] , correct_density: List[float] ) -> float
        Computes the Maximum Likelihood with a log-likelihood function
        """
        assert isinstance(incorrect_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in incorrect_density), 'arg incorrect_density wrong type'
        assert isinstance(correct_density, list) and all(isinstance(elemt_rec, float) for elemt_rec in correct_density), 'arg correct_density wrong type'
        cdef libcpp_vector[double] v0 = incorrect_density
        cdef libcpp_vector[double] v1 = correct_density
        cdef double _r = self.inst.get().computeLogLikelihood(v0, v1)
        correct_density[:] = v1
        incorrect_density[:] = v0
        py_result = <double>_r
        return py_result
    
    def pos_neg_mean_weighted_posteriors(self, list x_scores , list incorrect_posteriors ):
        """
        pos_neg_mean_weighted_posteriors(self, x_scores: List[float] , incorrect_posteriors: List[float] ) -> List[float, float]
        """
        assert isinstance(x_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in x_scores), 'arg x_scores wrong type'
        assert isinstance(incorrect_posteriors, list) and all(isinstance(elemt_rec, float) for elemt_rec in incorrect_posteriors), 'arg incorrect_posteriors wrong type'
        cdef libcpp_vector[double] v0 = x_scores
        cdef libcpp_vector[double] v1 = incorrect_posteriors
        _r = self.inst.get().pos_neg_mean_weighted_posteriors(v0, v1)
        incorrect_posteriors[:] = v1
        x_scores[:] = v0
        cdef list py_result = [_r.first, _r.second]
        return py_result
    
    def getCorrectlyAssignedFitResult(self):
        """
        getCorrectlyAssignedFitResult(self) -> GaussFitResult
        Returns estimated parameters for correctly assigned sequences. Fit should be used before
        """
        cdef _GaussFitResult * _r = new _GaussFitResult(self.inst.get().getCorrectlyAssignedFitResult())
        cdef GaussFitResult py_result = GaussFitResult.__new__(GaussFitResult)
        py_result.inst = shared_ptr[_GaussFitResult](_r)
        return py_result
    
    def getIncorrectlyAssignedFitResult(self):
        """
        getIncorrectlyAssignedFitResult(self) -> GaussFitResult
        Returns estimated parameters for correctly assigned sequences. Fit should be used before
        """
        cdef _GaussFitResult * _r = new _GaussFitResult(self.inst.get().getIncorrectlyAssignedFitResult())
        cdef GaussFitResult py_result = GaussFitResult.__new__(GaussFitResult)
        py_result.inst = shared_ptr[_GaussFitResult](_r)
        return py_result
    
    def getNegativePrior(self):
        """
        getNegativePrior(self) -> float
        Returns the estimated negative prior probability
        """
        cdef double _r = self.inst.get().getNegativePrior()
        py_result = <double>_r
        return py_result
    
    def computeProbability(self, double score ):
        """
        computeProbability(self, score: float ) -> float
        Returns the computed posterior error probability for a given score
        """
        assert isinstance(score, float), 'arg score wrong type'
    
        cdef double _r = self.inst.get().computeProbability((<double>score))
        py_result = <double>_r
        return py_result
    
    def initPlots(self, list x_scores ):
        """
        initPlots(self, x_scores: List[float] ) -> TextFile
        Initializes the plots
        """
        assert isinstance(x_scores, list) and all(isinstance(elemt_rec, float) for elemt_rec in x_scores), 'arg x_scores wrong type'
        cdef libcpp_vector[double] v0 = x_scores
        cdef _TextFile * _r = new _TextFile(self.inst.get().initPlots(v0))
        x_scores[:] = v0
        cdef TextFile py_result = TextFile.__new__(TextFile)
        py_result.inst = shared_ptr[_TextFile](_r)
        return py_result
    
    def getGumbelGnuplotFormula(self, GaussFitResult params ):
        """
        getGumbelGnuplotFormula(self, params: GaussFitResult ) -> Union[bytes, str, String]
        Returns the gnuplot formula of the fitted gumbel distribution
        """
        assert isinstance(params, GaussFitResult), 'arg params wrong type'
    
        cdef _String _r = self.inst.get().getGumbelGnuplotFormula((deref(params.inst.get())))
        py_result = convOutputString(_r)
        return py_result
    
    def getGaussGnuplotFormula(self, GaussFitResult params ):
        """
        getGaussGnuplotFormula(self, params: GaussFitResult ) -> Union[bytes, str, String]
        Returns the gnuplot formula of the fitted gauss distribution
        """
        assert isinstance(params, GaussFitResult), 'arg params wrong type'
    
        cdef _String _r = self.inst.get().getGaussGnuplotFormula((deref(params.inst.get())))
        py_result = convOutputString(_r)
        return py_result
    
    def getBothGnuplotFormula(self, GaussFitResult incorrect , GaussFitResult correct ):
        """
        getBothGnuplotFormula(self, incorrect: GaussFitResult , correct: GaussFitResult ) -> Union[bytes, str, String]
        Returns the gnuplot formula of the fitted mixture distribution
        """
        assert isinstance(incorrect, GaussFitResult), 'arg incorrect wrong type'
        assert isinstance(correct, GaussFitResult), 'arg correct wrong type'
    
    
        cdef _String _r = self.inst.get().getBothGnuplotFormula((deref(incorrect.inst.get())), (deref(correct.inst.get())))
        py_result = convOutputString(_r)
        return py_result
    
    def plotTargetDecoyEstimation(self, list target , list decoy ):
        """
        plotTargetDecoyEstimation(self, target: List[float] , decoy: List[float] ) -> None
        Plots the estimated distribution against target and decoy hits
        """
        assert isinstance(target, list) and all(isinstance(elemt_rec, float) for elemt_rec in target), 'arg target wrong type'
        assert isinstance(decoy, list) and all(isinstance(elemt_rec, float) for elemt_rec in decoy), 'arg decoy wrong type'
        cdef libcpp_vector[double] v0 = target
        cdef libcpp_vector[double] v1 = decoy
        self.inst.get().plotTargetDecoyEstimation(v0, v1)
        decoy[:] = v1
        target[:] = v0
    
    def getSmallestScore(self):
        """
        getSmallestScore(self) -> float
        Returns the smallest score used in the last fit
        """
        cdef double _r = self.inst.get().getSmallestScore()
        py_result = <double>_r
        return py_result
    
    def tryGnuplot(self,  gp_file ):
        """
        tryGnuplot(self, gp_file: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(gp_file, str) or isinstance(gp_file, bytes) or isinstance(gp_file, String)), 'arg gp_file wrong type'
    
        self.inst.get().tryGnuplot(deref((convString(gp_file)).get()))
    
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

cdef class RNaseDigestion:
    """
    Cython implementation of _RNaseDigestion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RNaseDigestion.html>`_
      -- Inherits from ['EnzymaticDigestion']

    Class for the enzymatic digestion of RNA
    
    Usage:
    
    .. code-block:: python
    
          from pyopenms import *
          oligo = NASequence.fromString("pAUGUCGCAG");
    
          dig = RNaseDigestion()
          dig.setEnzyme("RNase_T1")
    
          result = []
          dig.digest(oligo, result)
          for fragment in result:
            print (fragment)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RNaseDigestion rv = RNaseDigestion.__new__(RNaseDigestion)
       rv.inst = shared_ptr[_RNaseDigestion](new _RNaseDigestion(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RNaseDigestion rv = RNaseDigestion.__new__(RNaseDigestion)
       rv.inst = shared_ptr[_RNaseDigestion](new _RNaseDigestion(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RNaseDigestion](new _RNaseDigestion())
    
    def _init_1(self, RNaseDigestion in_0 ):
        """
        _init_1(self, in_0: RNaseDigestion ) -> None
        """
        assert isinstance(in_0, RNaseDigestion), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RNaseDigestion](new _RNaseDigestion((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RNaseDigestion ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RNaseDigestion)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setEnzyme_0(self,  name ):
        """
        _setEnzyme_0(self, name: Union[bytes, str, String] ) -> None
        Sets the enzyme for the digestion (by name)
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setEnzyme(deref((convString(name)).get()))
    
    def _setEnzyme_1(self, DigestionEnzyme enzyme ):
        """
        _setEnzyme_1(self, enzyme: DigestionEnzyme ) -> None
        Sets the enzyme for the digestion
        """
        assert isinstance(enzyme, DigestionEnzyme), 'arg enzyme wrong type'
    
        self.inst.get().setEnzyme((enzyme.inst.get()))
    
    def setEnzyme(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setEnzyme(self, name: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the enzyme for the digestion (by name)

        
        .. rubric:: Overload:
        .. py:function:: setEnzyme(self, enzyme: DigestionEnzyme ) -> None
          :noindex:
        
        Sets the enzyme for the digestion
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setEnzyme_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DigestionEnzyme)):
            return self._setEnzyme_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _digest_0(self, NASequence rna , list output ):
        """
        _digest_0(self, rna: NASequence , output: List[NASequence] ) -> None
        """
        assert isinstance(rna, NASequence), 'arg rna wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, NASequence) for elemt_rec in output), 'arg output wrong type'
    
        cdef libcpp_vector[_NASequence] * v1 = new libcpp_vector[_NASequence]()
        cdef NASequence item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().digest((deref(rna.inst.get())), deref(v1))
        cdef libcpp_vector[_NASequence].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = NASequence.__new__(NASequence)
            item1.inst = shared_ptr[_NASequence](new _NASequence(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
    
    def _digest_1(self, NASequence rna , list output ,  min_length ,  max_length ):
        """
        _digest_1(self, rna: NASequence , output: List[NASequence] , min_length: int , max_length: int ) -> None
        Performs the enzymatic digestion of a (potentially modified) RNA
        
        :param rna: Sequence to digest
        :param output: Digestion productsq
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
        Performs the enzymatic digestion of all RNA parent molecules in IdentificationData (id_data)
        
        :param id_data: IdentificationData object which includes sequences to digest
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
        """
        assert isinstance(rna, NASequence), 'arg rna wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, NASequence) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(min_length, int) and min_length >= 0, 'arg min_length wrong type'
        assert isinstance(max_length, int) and max_length >= 0, 'arg max_length wrong type'
    
        cdef libcpp_vector[_NASequence] * v1 = new libcpp_vector[_NASequence]()
        cdef NASequence item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().digest((deref(rna.inst.get())), deref(v1), (<size_t>min_length), (<size_t>max_length))
        cdef libcpp_vector[_NASequence].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = NASequence.__new__(NASequence)
            item1.inst = shared_ptr[_NASequence](new _NASequence(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
    
    def digest(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: digest(self, rna: NASequence , output: List[NASequence] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: digest(self, rna: NASequence , output: List[NASequence] , min_length: int , max_length: int ) -> None
          :noindex:
        
        Performs the enzymatic digestion of a (potentially modified) RNA
        
        :param rna: Sequence to digest
        :param output: Digestion productsq
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
        Performs the enzymatic digestion of all RNA parent molecules in IdentificationData (id_data)
        
        :param id_data: IdentificationData object which includes sequences to digest
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :returns: Number of discarded digestion products (which are not matching length restrictions)
    
        """
        if (len(args)==2) and (isinstance(args[0], NASequence)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, NASequence) for elemt_rec in args[1])):
            return self._digest_0(*args)
        elif (len(args)==4) and (isinstance(args[0], NASequence)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, NASequence) for elemt_rec in args[1])) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], int) and args[3] >= 0):
            return self._digest_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMissedCleavages(self):
        """
        getMissedCleavages(self) -> int
        Returns the max. number of allowed missed cleavages for the digestion
        """
        cdef size_t _r = self.inst.get().getMissedCleavages()
        py_result = <size_t>_r
        return py_result
    
    def setMissedCleavages(self,  missed_cleavages ):
        """
        setMissedCleavages(self, missed_cleavages: int ) -> None
        Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used
        """
        assert isinstance(missed_cleavages, int) and missed_cleavages >= 0, 'arg missed_cleavages wrong type'
    
        self.inst.get().setMissedCleavages((<size_t>missed_cleavages))
    
    def countInternalCleavageSites(self,  sequence ):
        """
        countInternalCleavageSites(self, sequence: Union[bytes, str, String] ) -> int
        Returns the number of internal cleavage sites for this sequence.
        """
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
    
        cdef size_t _r = self.inst.get().countInternalCleavageSites(deref((convString(sequence)).get()))
        py_result = <size_t>_r
        return py_result
    
    def getEnzymeName(self):
        """
        getEnzymeName(self) -> Union[bytes, str, String]
        Returns the enzyme for the digestion
        """
        cdef _String _r = self.inst.get().getEnzymeName()
        py_result = convOutputString(_r)
        return py_result
    
    def getSpecificity(self):
        """
        getSpecificity(self) -> int
        Returns the specificity for the digestion
        """
        cdef _Specificity _r = self.inst.get().getSpecificity()
        py_result = <int>_r
        return py_result
    
    def setSpecificity(self, int spec ):
        """
        setSpecificity(self, spec: int ) -> None
        Sets the specificity for the digestion (default is SPEC_FULL)
        """
        assert spec in [0, 1, 2, 3, 8, 9, 10], 'arg spec wrong type'
    
        self.inst.get().setSpecificity((<_Specificity>spec))
    
    def getSpecificityByName(self,  name ):
        """
        getSpecificityByName(self, name: Union[bytes, str, String] ) -> int
        Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _Specificity _r = self.inst.get().getSpecificityByName(deref((convString(name)).get()))
        py_result = <int>_r
        return py_result
    
    def digestUnmodified(self, StringView sequence , list output ,  min_length ,  max_length ):
        """
        digestUnmodified(self, sequence: StringView , output: List[StringView] , min_length: int , max_length: int ) -> int
        Performs the enzymatic digestion of an unmodified sequence\n
        By returning only references into the original string this is very fast
        
        
        :param sequence: Sequence to digest
        :param output: Digestion products
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        assert isinstance(sequence, StringView), 'arg sequence wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, StringView) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(min_length, int) and min_length >= 0, 'arg min_length wrong type'
        assert isinstance(max_length, int) and max_length >= 0, 'arg max_length wrong type'
    
        cdef libcpp_vector[_StringView] * v1 = new libcpp_vector[_StringView]()
        cdef StringView item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
    
    
        cdef size_t _r = self.inst.get().digestUnmodified((deref(sequence.inst.get())), deref(v1), (<size_t>min_length), (<size_t>max_length))
        cdef libcpp_vector[_StringView].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = StringView.__new__(StringView)
            item1.inst = shared_ptr[_StringView](new _StringView(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
        py_result = <size_t>_r
        return py_result
    
    def isValidProduct(self,  sequence ,  pos ,  length , bool ignore_missed_cleavages ):
        """
        isValidProduct(self, sequence: Union[bytes, str, String] , pos: int , length: int , ignore_missed_cleavages: bool ) -> bool
        Boolean operator returns true if the peptide fragment starting at position `pos` with length `length` within the sequence `sequence` generated by the current enzyme\n
        Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here
        
        
        :param protein: Protein sequence
        :param pep_pos: Starting index of potential peptide
        :param pep_length: Length of potential peptide
        :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
        :return: True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
        """
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
        assert isinstance(pos, int), 'arg pos wrong type'
        assert isinstance(length, int), 'arg length wrong type'
        assert isinstance(ignore_missed_cleavages, pybool_t), 'arg ignore_missed_cleavages wrong type'
    
    
    
    
        cdef bool _r = self.inst.get().isValidProduct(deref((convString(sequence)).get()), (<int>pos), (<int>length), (<bool>ignore_missed_cleavages))
        py_result = <bool>_r
        return py_result 

cdef class TextFile:
    """
    Cython implementation of _TextFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TextFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TextFile rv = TextFile.__new__(TextFile)
       rv.inst = shared_ptr[_TextFile](new _TextFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TextFile rv = TextFile.__new__(TextFile)
       rv.inst = shared_ptr[_TextFile](new _TextFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        This class provides some basic file handling methods for text files
        """
        self.inst = shared_ptr[_TextFile](new _TextFile())
    
    def _init_1(self, TextFile in_0 ):
        """
        _init_1(self, in_0: TextFile ) -> None
        """
        assert isinstance(in_0, TextFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TextFile](new _TextFile((deref(in_0.inst.get()))))
    
    def _init_2(self,  filename , bool trim_linesalse ,  first_n1 ):
        """
        _init_2(self, filename: Union[bytes, str, String] , trim_linesalse: bool , first_n1: int ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(trim_linesalse, pybool_t), 'arg trim_linesalse wrong type'
        assert isinstance(first_n1, int), 'arg first_n1 wrong type'
    
    
    
        self.inst = shared_ptr[_TextFile](new _TextFile(deref((convString(filename)).get()), (<bool>trim_linesalse), (<int>first_n1)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        This class provides some basic file handling methods for text files

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TextFile ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , trim_linesalse: bool , first_n1: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TextFile)):
             self._init_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], int)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , bool trim_linesalse ,  first_n1 ):
        """
        load(self, filename: Union[bytes, str, String] , trim_linesalse: bool , first_n1: int ) -> None
        Loads data from a text file
        
        :param filename: The input file name
        :param trim_lines: Whether or not the lines are trimmed when reading them from file
        :param first_n: If set, only `first_n` lines the lines from the beginning of the file are read
        :param skip_empty_lines: Should empty lines be skipped? If used in conjunction with `trim_lines`, also lines with only whitespace will be skipped. Skipped lines do not count towards the total number of read lines
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(trim_linesalse, pybool_t), 'arg trim_linesalse wrong type'
        assert isinstance(first_n1, int), 'arg first_n1 wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (<bool>trim_linesalse), (<int>first_n1))
    
    def store(self,  filename ):
        """
        store(self, filename: Union[bytes, str, String] ) -> None
        Writes the data to a file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().store(deref((convString(filename)).get()))
    
    def addLine(self,  line ):
        """
        addLine(self, line: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(line, str) or isinstance(line, bytes) or isinstance(line, String)), 'arg line wrong type'
    
        self.inst.get().addLine(deref((convString(line)).get())) 

cdef class XMLFile:
    """
    Cython implementation of _XMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1XMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef XMLFile rv = XMLFile.__new__(XMLFile)
       rv.inst = shared_ptr[_XMLFile](new _XMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef XMLFile rv = XMLFile.__new__(XMLFile)
       rv.inst = shared_ptr[_XMLFile](new _XMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_XMLFile](new _XMLFile())
    
    def _init_1(self, XMLFile in_0 ):
        """
        _init_1(self, in_0: XMLFile ) -> None
        """
        assert isinstance(in_0, XMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_XMLFile](new _XMLFile((deref(in_0.inst.get()))))
    
    def _init_2(self,  schema_location ,  version ):
        """
        _init_2(self, schema_location: Union[bytes, str, String] , version: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(schema_location, str) or isinstance(schema_location, bytes) or isinstance(schema_location, String)), 'arg schema_location wrong type'
        assert (isinstance(version, str) or isinstance(version, bytes) or isinstance(version, String)), 'arg version wrong type'
    
    
        self.inst = shared_ptr[_XMLFile](new _XMLFile(deref((convString(schema_location)).get()), deref((convString(version)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: XMLFile ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, schema_location: Union[bytes, str, String] , version: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], XMLFile)):
             self._init_1(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
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
