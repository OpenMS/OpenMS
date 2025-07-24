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
def __static_VersionDetails_create( in_0 ):
    """
    __static_VersionDetails_create(in_0: Union[bytes, str, String] ) -> VersionDetails
    """
    assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'

    cdef _VersionDetails * _r = new _VersionDetails(_create_VersionInfo(deref((convString(in_0)).get())))
    cdef VersionDetails py_result = VersionDetails.__new__(VersionDetails)
    py_result.inst = shared_ptr[_VersionDetails](_r)
    return py_result

def __static_VersionInfo_getBranch():
    """
    __static_VersionInfo_getBranch() -> Union[bytes, str, String]
    """
    cdef _String _r = _getBranch_VersionInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_VersionInfo_getRevision():
    """
    __static_VersionInfo_getRevision() -> Union[bytes, str, String]
    """
    cdef _String _r = _getRevision_VersionInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_VersionInfo_getTime():
    """
    __static_VersionInfo_getTime() -> Union[bytes, str, String]
    """
    cdef _String _r = _getTime_VersionInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_VersionInfo_getVersion():
    """
    __static_VersionInfo_getVersion() -> Union[bytes, str, String]
    """
    cdef _String _r = _getVersion_VersionInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_VersionInfo_getVersionStruct():
    """
    __static_VersionInfo_getVersionStruct() -> VersionDetails
    """
    cdef _VersionDetails * _r = new _VersionDetails(_getVersionStruct_VersionInfo())
    cdef VersionDetails py_result = VersionDetails.__new__(VersionDetails)
    py_result.inst = shared_ptr[_VersionDetails](_r)
    return py_result 

cdef class DRangeIntersection:
    None
    Disjoint = 0
    Intersects = 1
    Inside = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __PeptideIndexing_ExitCodes:
    None
    EXECUTION_OK = 0
    DATABASE_EMPTY = 1
    PEPTIDE_IDS_EMPTY = 2
    ILLEGAL_PARAMETERS = 3
    UNEXPECTED_RESULT = 4
    DECOYSTRING_EMPTY = 5

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __Specificity:
    None
    SPEC_NONE = 0
    SPEC_SEMI = 1
    SPEC_FULL = 2
    SPEC_UNKNOWN = 3
    SPEC_NOCTERM = 8
    SPEC_NONTERM = 9
    SIZE_OF_SPECIFICITY = 10

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __TermSpecificityNuc:
    None
    ANYWHERE = 0
    FIVE_PRIME = 1
    THREE_PRIME = 2
    NUMBER_OF_TERM_SPECIFICITY = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ConsensusIDAlgorithmPEPMatrix:
    """
    Cython implementation of _ConsensusIDAlgorithmPEPMatrix

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmPEPMatrix.html>`_
      -- Inherits from ['ConsensusIDAlgorithmSimilarity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmPEPMatrix](new _ConsensusIDAlgorithmPEPMatrix())
    
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

cdef class DRange1:
    """
    Cython implementation of _DRange1

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DRange1.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DRange1 rv = DRange1.__new__(DRange1)
       rv.inst = shared_ptr[_DRange1](new _DRange1(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DRange1 rv = DRange1.__new__(DRange1)
       rv.inst = shared_ptr[_DRange1](new _DRange1(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DRange1](new _DRange1())
    
    def _init_1(self, DRange1 in_0 ):
        """
        _init_1(self, in_0: DRange1 ) -> None
        """
        assert isinstance(in_0, DRange1), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DRange1](new _DRange1((deref(in_0.inst.get()))))
    
    def _init_2(self, DPosition1 lower , DPosition1 upper ):
        """
        _init_2(self, lower: DPosition1 , upper: DPosition1 ) -> None
        """
        assert isinstance(lower, DPosition1), 'arg lower wrong type'
        assert isinstance(upper, DPosition1), 'arg upper wrong type'
    
    
        self.inst = shared_ptr[_DRange1](new _DRange1((deref(lower.inst.get())), (deref(upper.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DRange1 ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, lower: DPosition1 , upper: DPosition1 ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DRange1)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], DPosition1)) and (isinstance(args[1], DPosition1)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def encloses(self, DPosition1 position ):
        """
        encloses(self, position: DPosition1 ) -> bool
        """
        assert isinstance(position, DPosition1), 'arg position wrong type'
    
        cdef bool _r = self.inst.get().encloses((deref(position.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def united(self, DRange1 other_range ):
        """
        united(self, other_range: DRange1 ) -> DRange1
        """
        assert isinstance(other_range, DRange1), 'arg other_range wrong type'
    
        cdef _DRange1 * _r = new _DRange1(self.inst.get().united((deref(other_range.inst.get()))))
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def isIntersected(self, DRange1 range_ ):
        """
        isIntersected(self, range_: DRange1 ) -> bool
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        cdef bool _r = self.inst.get().isIntersected((deref(range_.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def isEmpty(self):
        """
        isEmpty(self) -> bool
        """
        cdef bool _r = self.inst.get().isEmpty()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, DRange1):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DRange1 other_casted = other
        cdef DRange1 self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class DRange2:
    """
    Cython implementation of _DRange2

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DRange2.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DRange2 rv = DRange2.__new__(DRange2)
       rv.inst = shared_ptr[_DRange2](new _DRange2(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DRange2 rv = DRange2.__new__(DRange2)
       rv.inst = shared_ptr[_DRange2](new _DRange2(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DRange2](new _DRange2())
    
    def _init_1(self, DRange2 in_0 ):
        """
        _init_1(self, in_0: DRange2 ) -> None
        """
        assert isinstance(in_0, DRange2), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DRange2](new _DRange2((deref(in_0.inst.get()))))
    
    def _init_2(self,  lower ,  upper ):
        """
        _init_2(self, lower: Union[Sequence[int], Sequence[float]] , upper: Union[Sequence[int], Sequence[float]] ) -> None
        """
        assert len(lower) == 2 and isinstance(lower[0], (int, float)) and isinstance(lower[1], (int, float)), 'arg lower wrong type'
        assert len(upper) == 2 and isinstance(upper[0], (int, float)) and isinstance(upper[1], (int, float)), 'arg upper wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>lower[0]
        _dp_0[1] = <float>lower[1]
        cdef _DPosition2 _dp_1
        _dp_1[0] = <float>upper[0]
        _dp_1[1] = <float>upper[1]
        self.inst = shared_ptr[_DRange2](new _DRange2(_dp_0, _dp_1))
    
    def _init_3(self, double minx , double miny , double maxx , double maxy ):
        """
        _init_3(self, minx: float , miny: float , maxx: float , maxy: float ) -> None
        """
        assert isinstance(minx, float), 'arg minx wrong type'
        assert isinstance(miny, float), 'arg miny wrong type'
        assert isinstance(maxx, float), 'arg maxx wrong type'
        assert isinstance(maxy, float), 'arg maxy wrong type'
    
    
    
    
        self.inst = shared_ptr[_DRange2](new _DRange2((<double>minx), (<double>miny), (<double>maxx), (<double>maxy)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DRange2 ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, lower: Union[Sequence[int], Sequence[float]] , upper: Union[Sequence[int], Sequence[float]] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, minx: float , miny: float , maxx: float , maxy: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DRange2)):
             self._init_1(*args)
        elif (len(args)==2) and (len(args[0]) == 2 and isinstance(args[0][0], (int, float)) and isinstance(args[0][1], (int, float))) and (len(args[1]) == 2 and isinstance(args[1][0], (int, float)) and isinstance(args[1][1], (int, float))):
             self._init_2(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def united(self, DRange2 other_range ):
        """
        united(self, other_range: DRange2 ) -> DRange2
        """
        assert isinstance(other_range, DRange2), 'arg other_range wrong type'
    
        cdef _DRange2 * _r = new _DRange2(self.inst.get().united((deref(other_range.inst.get()))))
        cdef DRange2 py_result = DRange2.__new__(DRange2)
        py_result.inst = shared_ptr[_DRange2](_r)
        return py_result
    
    def isIntersected(self, DRange2 range_ ):
        """
        isIntersected(self, range_: DRange2 ) -> bool
        """
        assert isinstance(range_, DRange2), 'arg range_ wrong type'
    
        cdef bool _r = self.inst.get().isIntersected((deref(range_.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def isEmpty(self):
        """
        isEmpty(self) -> bool
        """
        cdef bool _r = self.inst.get().isEmpty()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, DRange2):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DRange2 other_casted = other
        cdef DRange2 self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class EnzymaticDigestion:
    """
    Cython implementation of _EnzymaticDigestion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EnzymaticDigestion.html>`_

      Class for the enzymatic digestion of proteins
    
      Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
      Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
      due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
      peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
      Thus no random selection of just n specific missed cleavage sites is performed.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EnzymaticDigestion rv = EnzymaticDigestion.__new__(EnzymaticDigestion)
       rv.inst = shared_ptr[_EnzymaticDigestion](new _EnzymaticDigestion(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EnzymaticDigestion rv = EnzymaticDigestion.__new__(EnzymaticDigestion)
       rv.inst = shared_ptr[_EnzymaticDigestion](new _EnzymaticDigestion(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_EnzymaticDigestion](new _EnzymaticDigestion())
    
    def _init_1(self, EnzymaticDigestion in_0 ):
        """
        _init_1(self, in_0: EnzymaticDigestion ) -> None
        """
        assert isinstance(in_0, EnzymaticDigestion), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EnzymaticDigestion](new _EnzymaticDigestion((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EnzymaticDigestion ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EnzymaticDigestion)):
             self._init_1(*args)
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
    
    def setEnzyme(self, DigestionEnzyme enzyme ):
        """
        setEnzyme(self, enzyme: DigestionEnzyme ) -> None
        Sets the enzyme for the digestion
        """
        assert isinstance(enzyme, DigestionEnzyme), 'arg enzyme wrong type'
    
        self.inst.get().setEnzyme((enzyme.inst.get()))
    
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
    Specificity = __Specificity 

cdef class FeatureGroupingAlgorithmQT:
    """
    Cython implementation of _FeatureGroupingAlgorithmQT

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmQT.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureGroupingAlgorithmQT](new _FeatureGroupingAlgorithmQT())
    
    def _group_0(self, list maps , ConsensusMap out ):
        """
        _group_0(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in maps), 'arg maps wrong type'
        assert isinstance(out, ConsensusMap), 'arg out wrong type'
        cdef libcpp_vector[_FeatureMap] * v0 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().group(deref(v0), (deref(out.inst.get())))
        cdef libcpp_vector[_FeatureMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = FeatureMap.__new__(FeatureMap)
            item0.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def _group_1(self, list maps , ConsensusMap out ):
        """
        _group_1(self, maps: List[ConsensusMap] , out: ConsensusMap ) -> None
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in maps), 'arg maps wrong type'
        assert isinstance(out, ConsensusMap), 'arg out wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().group(deref(v0), (deref(out.inst.get())))
        cdef libcpp_vector[_ConsensusMap].iterator it_maps = v0.begin()
        replace_0 = []
        while it_maps != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_maps)))
            replace_0.append(item0)
            inc(it_maps)
        maps[:] = replace_0
        del v0
    
    def group(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: group(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: group(self, maps: List[ConsensusMap] , out: ConsensusMap ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in args[0])) and (isinstance(args[1], ConsensusMap)):
            return self._group_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in args[0])) and (isinstance(args[1], ConsensusMap)):
            return self._group_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def transferSubelements(self, list maps , ConsensusMap out ):
        """
        transferSubelements(self, maps: List[ConsensusMap] , out: ConsensusMap ) -> None
        Transfers subelements (grouped features) from input consensus maps to the result consensus map
        """
        assert isinstance(maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in maps), 'arg maps wrong type'
        assert isinstance(out, ConsensusMap), 'arg out wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().transferSubelements(deref(v0), (deref(out.inst.get())))
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

cdef class LightCompound:
    """
    Cython implementation of _LightCompound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightCompound.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property rt:
        def __set__(self, double rt):
        
            self.inst.get().rt = (<double>rt)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt
            py_result = <double>_r
            return py_result
    
    property drift_time:
        def __set__(self, double drift_time):
        
            self.inst.get().drift_time = (<double>drift_time)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().drift_time
            py_result = <double>_r
            return py_result
    
    property charge:
        def __set__(self,  charge):
        
            self.inst.get().charge = (<int>charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().charge
            py_result = <int>_r
            return py_result
    
    property sequence:
        def __set__(self, bytes sequence):
        
            self.inst.get().sequence = (<libcpp_string>sequence)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().sequence
            py_result = <libcpp_string>_r
            return py_result
    
    property protein_refs:
        def __set__(self, list protein_refs):
            cdef libcpp_vector[libcpp_string] v0 = protein_refs
            self.inst.get().protein_refs = v0
            
    
        def __get__(self):
            _r = self.inst.get().protein_refs
            cdef list py_result = _r
            return py_result
    
    property peptide_group_label:
        def __set__(self, bytes peptide_group_label):
        
            self.inst.get().peptide_group_label = (<libcpp_string>peptide_group_label)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().peptide_group_label
            py_result = <libcpp_string>_r
            return py_result
    
    property id:
        def __set__(self, bytes id):
        
            self.inst.get().id = (<libcpp_string>id)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().id
            py_result = <libcpp_string>_r
            return py_result
    
    property sum_formula:
        def __set__(self, bytes sum_formula):
        
            self.inst.get().sum_formula = (<libcpp_string>sum_formula)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().sum_formula
            py_result = <libcpp_string>_r
            return py_result
    
    property compound_name:
        def __set__(self, bytes compound_name):
        
            self.inst.get().compound_name = (<libcpp_string>compound_name)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().compound_name
            py_result = <libcpp_string>_r
            return py_result
    
    property modifications:
        def __set__(self, list modifications):
            cdef libcpp_vector[_LightModification] * v0 = new libcpp_vector[_LightModification]()
            cdef LightModification item0
            for item0 in modifications:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().modifications = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().modifications
            py_result = []
            cdef libcpp_vector[_LightModification].iterator it__r = _r.begin()
            cdef LightModification item_py_result
            while it__r != _r.end():
               item_py_result = LightModification.__new__(LightModification)
               item_py_result.inst = shared_ptr[_LightModification](new _LightModification(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef LightCompound rv = LightCompound.__new__(LightCompound)
       rv.inst = shared_ptr[_LightCompound](new _LightCompound(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightCompound rv = LightCompound.__new__(LightCompound)
       rv.inst = shared_ptr[_LightCompound](new _LightCompound(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LightCompound](new _LightCompound())
    
    def _init_1(self, LightCompound in_0 ):
        """
        _init_1(self, in_0: LightCompound ) -> None
        """
        assert isinstance(in_0, LightCompound), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LightCompound](new _LightCompound((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightCompound ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightCompound)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setDriftTime(self, double d ):
        """
        setDriftTime(self, d: float ) -> None
        """
        assert isinstance(d, float), 'arg d wrong type'
    
        self.inst.get().setDriftTime((<double>d))
    
    def getDriftTime(self):
        """
        getDriftTime(self) -> float
        """
        cdef double _r = self.inst.get().getDriftTime()
        py_result = <double>_r
        return py_result
    
    def getChargeState(self):
        """
        getChargeState(self) -> int
        """
        cdef int _r = self.inst.get().getChargeState()
        py_result = <int>_r
        return py_result
    
    def isPeptide(self):
        """
        isPeptide(self) -> bool
        """
        cdef bool _r = self.inst.get().isPeptide()
        py_result = <bool>_r
        return py_result
    
    def setChargeState(self,  ch ):
        """
        setChargeState(self, ch: int ) -> None
        """
        assert isinstance(ch, int), 'arg ch wrong type'
    
        self.inst.get().setChargeState((<int>ch)) 

cdef class LightMRMTransitionGroupCP:
    """
    Cython implementation of _MRMTransitionGroup[_MSChromatogram,_LightTransition]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroup[_MSChromatogram,_LightTransition].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef LightMRMTransitionGroupCP rv = LightMRMTransitionGroupCP.__new__(LightMRMTransitionGroupCP)
       rv.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_LightTransition]](new _MRMTransitionGroup[_MSChromatogram,_LightTransition](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightMRMTransitionGroupCP rv = LightMRMTransitionGroupCP.__new__(LightMRMTransitionGroupCP)
       rv.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_LightTransition]](new _MRMTransitionGroup[_MSChromatogram,_LightTransition](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_LightTransition]](new _MRMTransitionGroup[_MSChromatogram,_LightTransition]())
    
    def _init_1(self, LightMRMTransitionGroupCP in_0 ):
        """
        _init_1(self, in_0: LightMRMTransitionGroupCP ) -> None
        """
        assert isinstance(in_0, LightMRMTransitionGroupCP), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_LightTransition]](new _MRMTransitionGroup[_MSChromatogram,_LightTransition]((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightMRMTransitionGroupCP ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightMRMTransitionGroupCP)):
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
    
    def getTransitionGroupID(self):
        """
        getTransitionGroupID(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTransitionGroupID()
        py_result = convOutputString(_r)
        return py_result
    
    def setTransitionGroupID(self,  tr_gr_id ):
        """
        setTransitionGroupID(self, tr_gr_id: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(tr_gr_id, str) or isinstance(tr_gr_id, bytes) or isinstance(tr_gr_id, String)), 'arg tr_gr_id wrong type'
    
        self.inst.get().setTransitionGroupID(deref((convString(tr_gr_id)).get()))
    
    def getTransitions(self):
        """
        getTransitions(self) -> List[LightTransition]
        """
        _r = self.inst.get().getTransitions()
        py_result = []
        cdef libcpp_vector[_LightTransition].iterator it__r = _r.begin()
        cdef LightTransition item_py_result
        while it__r != _r.end():
           item_py_result = LightTransition.__new__(LightTransition)
           item_py_result.inst = shared_ptr[_LightTransition](new _LightTransition(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getTransitionsMuteable(self):
        """
        getTransitionsMuteable(self) -> List[LightTransition]
        """
        _r = self.inst.get().getTransitionsMuteable()
        py_result = []
        cdef libcpp_vector[_LightTransition].iterator it__r = _r.begin()
        cdef LightTransition item_py_result
        while it__r != _r.end():
           item_py_result = LightTransition.__new__(LightTransition)
           item_py_result.inst = shared_ptr[_LightTransition](new _LightTransition(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addTransition(self, LightTransition transition ,  key ):
        """
        addTransition(self, transition: LightTransition , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(transition, LightTransition), 'arg transition wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addTransition((deref(transition.inst.get())), deref((convString(key)).get()))
    
    def getTransition(self,  key ):
        """
        getTransition(self, key: Union[bytes, str, String] ) -> LightTransition
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _LightTransition * _r = new _LightTransition(self.inst.get().getTransition(deref((convString(key)).get())))
        cdef LightTransition py_result = LightTransition.__new__(LightTransition)
        py_result.inst = shared_ptr[_LightTransition](_r)
        return py_result
    
    def hasTransition(self,  key ):
        """
        hasTransition(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasTransition(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getChromatograms(self):
        """
        getChromatograms(self) -> List[MSChromatogram]
        """
        _r = self.inst.get().getChromatograms()
        py_result = []
        cdef libcpp_vector[_MSChromatogram].iterator it__r = _r.begin()
        cdef MSChromatogram item_py_result
        while it__r != _r.end():
           item_py_result = MSChromatogram.__new__(MSChromatogram)
           item_py_result.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addChromatogram(self, MSChromatogram chromatogram ,  key ):
        """
        addChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addChromatogram((deref(chromatogram.inst.get())), deref((convString(key)).get()))
    
    def getChromatogram(self,  key ):
        """
        getChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram(deref((convString(key)).get())))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def hasChromatogram(self,  key ):
        """
        hasChromatogram(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasChromatogram(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getPrecursorChromatograms(self):
        """
        getPrecursorChromatograms(self) -> List[MSChromatogram]
        """
        _r = self.inst.get().getPrecursorChromatograms()
        py_result = []
        cdef libcpp_vector[_MSChromatogram].iterator it__r = _r.begin()
        cdef MSChromatogram item_py_result
        while it__r != _r.end():
           item_py_result = MSChromatogram.__new__(MSChromatogram)
           item_py_result.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addPrecursorChromatogram(self, MSChromatogram chromatogram ,  key ):
        """
        addPrecursorChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addPrecursorChromatogram((deref(chromatogram.inst.get())), deref((convString(key)).get()))
    
    def getPrecursorChromatogram(self,  key ):
        """
        getPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getPrecursorChromatogram(deref((convString(key)).get())))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def hasPrecursorChromatogram(self,  key ):
        """
        hasPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasPrecursorChromatogram(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getFeatures(self):
        """
        getFeatures(self) -> List[MRMFeature]
        """
        _r = self.inst.get().getFeatures()
        py_result = []
        cdef libcpp_vector[_MRMFeature].iterator it__r = _r.begin()
        cdef MRMFeature item_py_result
        while it__r != _r.end():
           item_py_result = MRMFeature.__new__(MRMFeature)
           item_py_result.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getFeaturesMuteable(self):
        """
        getFeaturesMuteable(self) -> List[MRMFeature]
        """
        _r = self.inst.get().getFeaturesMuteable()
        py_result = []
        cdef libcpp_vector[_MRMFeature].iterator it__r = _r.begin()
        cdef MRMFeature item_py_result
        while it__r != _r.end():
           item_py_result = MRMFeature.__new__(MRMFeature)
           item_py_result.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addFeature(self, MRMFeature feature ):
        """
        addFeature(self, feature: MRMFeature ) -> None
        """
        assert isinstance(feature, MRMFeature), 'arg feature wrong type'
    
        self.inst.get().addFeature((deref(feature.inst.get())))
    
    def getBestFeature(self):
        """
        getBestFeature(self) -> MRMFeature
        """
        cdef _MRMFeature * _r = new _MRMFeature(self.inst.get().getBestFeature())
        cdef MRMFeature py_result = MRMFeature.__new__(MRMFeature)
        py_result.inst = shared_ptr[_MRMFeature](_r)
        return py_result
    
    def getLibraryIntensity(self, list result ):
        """
        getLibraryIntensity(self, result: List[float] ) -> None
        """
        assert isinstance(result, list) and all(isinstance(elemt_rec, float) for elemt_rec in result), 'arg result wrong type'
        cdef libcpp_vector[double] v0 = result
        self.inst.get().getLibraryIntensity(v0)
        
    
    def subset(self, list tr_ids ):
        """
        subset(self, tr_ids: List[Union[bytes, str]] ) -> LightMRMTransitionGroupCP
        """
        assert isinstance(tr_ids, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tr_ids), 'arg tr_ids wrong type'
        cdef libcpp_vector[libcpp_utf8_string] v0 = tr_ids
        cdef _MRMTransitionGroup[_MSChromatogram,_LightTransition] * _r = new _MRMTransitionGroup[_MSChromatogram,_LightTransition](self.inst.get().subset(v0))
        
        cdef LightMRMTransitionGroupCP py_result = LightMRMTransitionGroupCP.__new__(LightMRMTransitionGroupCP)
        py_result.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_LightTransition]](_r)
        return py_result
    
    def isInternallyConsistent(self):
        """
        isInternallyConsistent(self) -> bool
        """
        cdef bool _r = self.inst.get().isInternallyConsistent()
        py_result = <bool>_r
        return py_result
    
    def chromatogramIdsMatch(self):
        """
        chromatogramIdsMatch(self) -> bool
        """
        cdef bool _r = self.inst.get().chromatogramIdsMatch()
        py_result = <bool>_r
        return py_result 

cdef class LightModification:
    """
    Cython implementation of _LightModification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightModification.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
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
       cdef LightModification rv = LightModification.__new__(LightModification)
       rv.inst = shared_ptr[_LightModification](new _LightModification(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightModification rv = LightModification.__new__(LightModification)
       rv.inst = shared_ptr[_LightModification](new _LightModification(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LightModification](new _LightModification())
    
    def _init_1(self, LightModification in_0 ):
        """
        _init_1(self, in_0: LightModification ) -> None
        """
        assert isinstance(in_0, LightModification), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LightModification](new _LightModification((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightModification ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightModification)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class LightProtein:
    """
    Cython implementation of _LightProtein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightProtein.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property id:
        def __set__(self, bytes id):
        
            self.inst.get().id = (<libcpp_string>id)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().id
            py_result = <libcpp_string>_r
            return py_result
    
    property sequence:
        def __set__(self, bytes sequence):
        
            self.inst.get().sequence = (<libcpp_string>sequence)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().sequence
            py_result = <libcpp_string>_r
            return py_result
    
    def __copy__(self):
       cdef LightProtein rv = LightProtein.__new__(LightProtein)
       rv.inst = shared_ptr[_LightProtein](new _LightProtein(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightProtein rv = LightProtein.__new__(LightProtein)
       rv.inst = shared_ptr[_LightProtein](new _LightProtein(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LightProtein](new _LightProtein())
    
    def _init_1(self, LightProtein in_0 ):
        """
        _init_1(self, in_0: LightProtein ) -> None
        """
        assert isinstance(in_0, LightProtein), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LightProtein](new _LightProtein((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightProtein ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightProtein)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class LightTargetedExperiment:
    """
    Cython implementation of _LightTargetedExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightTargetedExperiment.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property transitions:
        def __set__(self, list transitions):
            cdef libcpp_vector[_LightTransition] * v0 = new libcpp_vector[_LightTransition]()
            cdef LightTransition item0
            for item0 in transitions:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().transitions = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().transitions
            py_result = []
            cdef libcpp_vector[_LightTransition].iterator it__r = _r.begin()
            cdef LightTransition item_py_result
            while it__r != _r.end():
               item_py_result = LightTransition.__new__(LightTransition)
               item_py_result.inst = shared_ptr[_LightTransition](new _LightTransition(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property compounds:
        def __set__(self, list compounds):
            cdef libcpp_vector[_LightCompound] * v0 = new libcpp_vector[_LightCompound]()
            cdef LightCompound item0
            for item0 in compounds:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().compounds = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().compounds
            py_result = []
            cdef libcpp_vector[_LightCompound].iterator it__r = _r.begin()
            cdef LightCompound item_py_result
            while it__r != _r.end():
               item_py_result = LightCompound.__new__(LightCompound)
               item_py_result.inst = shared_ptr[_LightCompound](new _LightCompound(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property proteins:
        def __set__(self, list proteins):
            cdef libcpp_vector[_LightProtein] * v0 = new libcpp_vector[_LightProtein]()
            cdef LightProtein item0
            for item0 in proteins:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().proteins = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().proteins
            py_result = []
            cdef libcpp_vector[_LightProtein].iterator it__r = _r.begin()
            cdef LightProtein item_py_result
            while it__r != _r.end():
               item_py_result = LightProtein.__new__(LightProtein)
               item_py_result.inst = shared_ptr[_LightProtein](new _LightProtein(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef LightTargetedExperiment rv = LightTargetedExperiment.__new__(LightTargetedExperiment)
       rv.inst = shared_ptr[_LightTargetedExperiment](new _LightTargetedExperiment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightTargetedExperiment rv = LightTargetedExperiment.__new__(LightTargetedExperiment)
       rv.inst = shared_ptr[_LightTargetedExperiment](new _LightTargetedExperiment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LightTargetedExperiment](new _LightTargetedExperiment())
    
    def _init_1(self, LightTargetedExperiment in_0 ):
        """
        _init_1(self, in_0: LightTargetedExperiment ) -> None
        """
        assert isinstance(in_0, LightTargetedExperiment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LightTargetedExperiment](new _LightTargetedExperiment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightTargetedExperiment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightTargetedExperiment)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getTransitions(self):
        """
        getTransitions(self) -> List[LightTransition]
        """
        _r = self.inst.get().getTransitions()
        py_result = []
        cdef libcpp_vector[_LightTransition].iterator it__r = _r.begin()
        cdef LightTransition item_py_result
        while it__r != _r.end():
           item_py_result = LightTransition.__new__(LightTransition)
           item_py_result.inst = shared_ptr[_LightTransition](new _LightTransition(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getCompounds(self):
        """
        getCompounds(self) -> List[LightCompound]
        """
        _r = self.inst.get().getCompounds()
        py_result = []
        cdef libcpp_vector[_LightCompound].iterator it__r = _r.begin()
        cdef LightCompound item_py_result
        while it__r != _r.end():
           item_py_result = LightCompound.__new__(LightCompound)
           item_py_result.inst = shared_ptr[_LightCompound](new _LightCompound(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getProteins(self):
        """
        getProteins(self) -> List[LightProtein]
        """
        _r = self.inst.get().getProteins()
        py_result = []
        cdef libcpp_vector[_LightProtein].iterator it__r = _r.begin()
        cdef LightProtein item_py_result
        while it__r != _r.end():
           item_py_result = LightProtein.__new__(LightProtein)
           item_py_result.inst = shared_ptr[_LightProtein](new _LightProtein(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getCompoundByRef(self, bytes ref ):
        """
        getCompoundByRef(self, ref: bytes ) -> LightCompound
        """
        assert isinstance(ref, bytes), 'arg ref wrong type'
    
        cdef _LightCompound * _r = new _LightCompound(self.inst.get().getCompoundByRef((<libcpp_string>ref)))
        cdef LightCompound py_result = LightCompound.__new__(LightCompound)
        py_result.inst = shared_ptr[_LightCompound](_r)
        return py_result
    
    def getPeptideByRef(self, bytes ref ):
        """
        getPeptideByRef(self, ref: bytes ) -> LightCompound
        """
        assert isinstance(ref, bytes), 'arg ref wrong type'
    
        cdef _LightCompound * _r = new _LightCompound(self.inst.get().getPeptideByRef((<libcpp_string>ref)))
        cdef LightCompound py_result = LightCompound.__new__(LightCompound)
        py_result.inst = shared_ptr[_LightCompound](_r)
        return py_result 

cdef class LightTransition:
    """
    Cython implementation of _LightTransition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1LightTransition.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property transition_name:
        def __set__(self, bytes transition_name):
        
            self.inst.get().transition_name = (<libcpp_string>transition_name)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().transition_name
            py_result = <libcpp_string>_r
            return py_result
    
    property transition_ref:
        def __set__(self, bytes transition_ref):
        
            self.inst.get().transition_ref = (<libcpp_string>transition_ref)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().transition_ref
            py_result = <libcpp_string>_r
            return py_result
    
    property library_intensity:
        def __set__(self, double library_intensity):
        
            self.inst.get().library_intensity = (<double>library_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_intensity
            py_result = <double>_r
            return py_result
    
    property product_mz:
        def __set__(self, double product_mz):
        
            self.inst.get().product_mz = (<double>product_mz)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().product_mz
            py_result = <double>_r
            return py_result
    
    property precursor_mz:
        def __set__(self, double precursor_mz):
        
            self.inst.get().precursor_mz = (<double>precursor_mz)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_mz
            py_result = <double>_r
            return py_result
    
    property fragment_charge:
        def __set__(self,  fragment_charge):
        
            self.inst.get().fragment_charge = (<int>fragment_charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().fragment_charge
            py_result = <int>_r
            return py_result
    
    property decoy:
        def __set__(self, bool decoy):
        
            self.inst.get().decoy = (<bool>decoy)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().decoy
            py_result = <bool>_r
            return py_result
    
    property detecting_transition:
        def __set__(self, bool detecting_transition):
        
            self.inst.get().detecting_transition = (<bool>detecting_transition)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().detecting_transition
            py_result = <bool>_r
            return py_result
    
    property quantifying_transition:
        def __set__(self, bool quantifying_transition):
        
            self.inst.get().quantifying_transition = (<bool>quantifying_transition)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().quantifying_transition
            py_result = <bool>_r
            return py_result
    
    property identifying_transition:
        def __set__(self, bool identifying_transition):
        
            self.inst.get().identifying_transition = (<bool>identifying_transition)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().identifying_transition
            py_result = <bool>_r
            return py_result
    
    def __copy__(self):
       cdef LightTransition rv = LightTransition.__new__(LightTransition)
       rv.inst = shared_ptr[_LightTransition](new _LightTransition(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LightTransition rv = LightTransition.__new__(LightTransition)
       rv.inst = shared_ptr[_LightTransition](new _LightTransition(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LightTransition](new _LightTransition())
    
    def _init_1(self, LightTransition in_0 ):
        """
        _init_1(self, in_0: LightTransition ) -> None
        """
        assert isinstance(in_0, LightTransition), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LightTransition](new _LightTransition((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LightTransition ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LightTransition)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getProductChargeState(self):
        """
        getProductChargeState(self) -> int
        """
        cdef int _r = self.inst.get().getProductChargeState()
        py_result = <int>_r
        return py_result
    
    def isProductChargeStateSet(self):
        """
        isProductChargeStateSet(self) -> bool
        """
        cdef bool _r = self.inst.get().isProductChargeStateSet()
        py_result = <bool>_r
        return py_result
    
    def getNativeID(self):
        """
        getNativeID(self) -> bytes
        """
        cdef libcpp_string _r = self.inst.get().getNativeID()
        py_result = <libcpp_string>_r
        return py_result
    
    def getTransRef(self):
        """
        getTransRef(self) -> bytes
        """
        cdef libcpp_string _r = self.inst.get().getTransRef()
        py_result = <libcpp_string>_r
        return py_result
    
    def getLibraryIntensity(self):
        """
        getLibraryIntensity(self) -> float
        """
        cdef double _r = self.inst.get().getLibraryIntensity()
        py_result = <double>_r
        return py_result
    
    def setLibraryIntensity(self, double l ):
        """
        setLibraryIntensity(self, l: float ) -> None
        """
        assert isinstance(l, float), 'arg l wrong type'
    
        self.inst.get().setLibraryIntensity((<double>l))
    
    def getProductMZ(self):
        """
        getProductMZ(self) -> float
        """
        cdef double _r = self.inst.get().getProductMZ()
        py_result = <double>_r
        return py_result
    
    def getPrecursorMZ(self):
        """
        getPrecursorMZ(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorMZ()
        py_result = <double>_r
        return py_result
    
    def setDetectingTransition(self, bool d ):
        """
        setDetectingTransition(self, d: bool ) -> None
        """
        assert isinstance(d, pybool_t), 'arg d wrong type'
    
        self.inst.get().setDetectingTransition((<bool>d))
    
    def isDetectingTransition(self):
        """
        isDetectingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isDetectingTransition()
        py_result = <bool>_r
        return py_result
    
    def setQuantifyingTransition(self, bool q ):
        """
        setQuantifyingTransition(self, q: bool ) -> None
        """
        assert isinstance(q, pybool_t), 'arg q wrong type'
    
        self.inst.get().setQuantifyingTransition((<bool>q))
    
    def isQuantifyingTransition(self):
        """
        isQuantifyingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isQuantifyingTransition()
        py_result = <bool>_r
        return py_result
    
    def setIdentifyingTransition(self, bool i ):
        """
        setIdentifyingTransition(self, i: bool ) -> None
        """
        assert isinstance(i, pybool_t), 'arg i wrong type'
    
        self.inst.get().setIdentifyingTransition((<bool>i))
    
    def isIdentifyingTransition(self):
        """
        isIdentifyingTransition(self) -> bool
        """
        cdef bool _r = self.inst.get().isIdentifyingTransition()
        py_result = <bool>_r
        return py_result 

cdef class MRMTransitionGroupCP:
    """
    Cython implementation of _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMTransitionGroupCP rv = MRMTransitionGroupCP.__new__(MRMTransitionGroupCP)
       rv.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]](new _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMTransitionGroupCP rv = MRMTransitionGroupCP.__new__(MRMTransitionGroupCP)
       rv.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]](new _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]](new _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]())
    
    def _init_1(self, MRMTransitionGroupCP in_0 ):
        """
        _init_1(self, in_0: MRMTransitionGroupCP ) -> None
        """
        assert isinstance(in_0, MRMTransitionGroupCP), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]](new _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMTransitionGroupCP ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMTransitionGroupCP)):
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
    
    def getTransitionGroupID(self):
        """
        getTransitionGroupID(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTransitionGroupID()
        py_result = convOutputString(_r)
        return py_result
    
    def setTransitionGroupID(self,  tr_gr_id ):
        """
        setTransitionGroupID(self, tr_gr_id: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(tr_gr_id, str) or isinstance(tr_gr_id, bytes) or isinstance(tr_gr_id, String)), 'arg tr_gr_id wrong type'
    
        self.inst.get().setTransitionGroupID(deref((convString(tr_gr_id)).get()))
    
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
    
    def getTransitionsMuteable(self):
        """
        getTransitionsMuteable(self) -> List[ReactionMonitoringTransition]
        """
        _r = self.inst.get().getTransitionsMuteable()
        py_result = []
        cdef libcpp_vector[_ReactionMonitoringTransition].iterator it__r = _r.begin()
        cdef ReactionMonitoringTransition item_py_result
        while it__r != _r.end():
           item_py_result = ReactionMonitoringTransition.__new__(ReactionMonitoringTransition)
           item_py_result.inst = shared_ptr[_ReactionMonitoringTransition](new _ReactionMonitoringTransition(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addTransition(self, ReactionMonitoringTransition transition ,  key ):
        """
        addTransition(self, transition: ReactionMonitoringTransition , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(transition, ReactionMonitoringTransition), 'arg transition wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addTransition((deref(transition.inst.get())), deref((convString(key)).get()))
    
    def getTransition(self,  key ):
        """
        getTransition(self, key: Union[bytes, str, String] ) -> ReactionMonitoringTransition
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _ReactionMonitoringTransition * _r = new _ReactionMonitoringTransition(self.inst.get().getTransition(deref((convString(key)).get())))
        cdef ReactionMonitoringTransition py_result = ReactionMonitoringTransition.__new__(ReactionMonitoringTransition)
        py_result.inst = shared_ptr[_ReactionMonitoringTransition](_r)
        return py_result
    
    def hasTransition(self,  key ):
        """
        hasTransition(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasTransition(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getChromatograms(self):
        """
        getChromatograms(self) -> List[MSChromatogram]
        """
        _r = self.inst.get().getChromatograms()
        py_result = []
        cdef libcpp_vector[_MSChromatogram].iterator it__r = _r.begin()
        cdef MSChromatogram item_py_result
        while it__r != _r.end():
           item_py_result = MSChromatogram.__new__(MSChromatogram)
           item_py_result.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addChromatogram(self, MSChromatogram chromatogram ,  key ):
        """
        addChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addChromatogram((deref(chromatogram.inst.get())), deref((convString(key)).get()))
    
    def getChromatogram(self,  key ):
        """
        getChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram(deref((convString(key)).get())))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def hasChromatogram(self,  key ):
        """
        hasChromatogram(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasChromatogram(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getPrecursorChromatograms(self):
        """
        getPrecursorChromatograms(self) -> List[MSChromatogram]
        """
        _r = self.inst.get().getPrecursorChromatograms()
        py_result = []
        cdef libcpp_vector[_MSChromatogram].iterator it__r = _r.begin()
        cdef MSChromatogram item_py_result
        while it__r != _r.end():
           item_py_result = MSChromatogram.__new__(MSChromatogram)
           item_py_result.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addPrecursorChromatogram(self, MSChromatogram chromatogram ,  key ):
        """
        addPrecursorChromatogram(self, chromatogram: MSChromatogram , key: Union[bytes, str, String] ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addPrecursorChromatogram((deref(chromatogram.inst.get())), deref((convString(key)).get()))
    
    def getPrecursorChromatogram(self,  key ):
        """
        getPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> MSChromatogram
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getPrecursorChromatogram(deref((convString(key)).get())))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def hasPrecursorChromatogram(self,  key ):
        """
        hasPrecursorChromatogram(self, key: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef bool _r = self.inst.get().hasPrecursorChromatogram(deref((convString(key)).get()))
        py_result = <bool>_r
        return py_result
    
    def getFeatures(self):
        """
        getFeatures(self) -> List[MRMFeature]
        """
        _r = self.inst.get().getFeatures()
        py_result = []
        cdef libcpp_vector[_MRMFeature].iterator it__r = _r.begin()
        cdef MRMFeature item_py_result
        while it__r != _r.end():
           item_py_result = MRMFeature.__new__(MRMFeature)
           item_py_result.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getFeaturesMuteable(self):
        """
        getFeaturesMuteable(self) -> List[MRMFeature]
        """
        _r = self.inst.get().getFeaturesMuteable()
        py_result = []
        cdef libcpp_vector[_MRMFeature].iterator it__r = _r.begin()
        cdef MRMFeature item_py_result
        while it__r != _r.end():
           item_py_result = MRMFeature.__new__(MRMFeature)
           item_py_result.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addFeature(self, MRMFeature feature ):
        """
        addFeature(self, feature: MRMFeature ) -> None
        """
        assert isinstance(feature, MRMFeature), 'arg feature wrong type'
    
        self.inst.get().addFeature((deref(feature.inst.get())))
    
    def getBestFeature(self):
        """
        getBestFeature(self) -> MRMFeature
        """
        cdef _MRMFeature * _r = new _MRMFeature(self.inst.get().getBestFeature())
        cdef MRMFeature py_result = MRMFeature.__new__(MRMFeature)
        py_result.inst = shared_ptr[_MRMFeature](_r)
        return py_result
    
    def getLibraryIntensity(self, list result ):
        """
        getLibraryIntensity(self, result: List[float] ) -> None
        """
        assert isinstance(result, list) and all(isinstance(elemt_rec, float) for elemt_rec in result), 'arg result wrong type'
        cdef libcpp_vector[double] v0 = result
        self.inst.get().getLibraryIntensity(v0)
        
    
    def subset(self, list tr_ids ):
        """
        subset(self, tr_ids: List[Union[bytes, str]] ) -> MRMTransitionGroupCP
        """
        assert isinstance(tr_ids, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tr_ids), 'arg tr_ids wrong type'
        cdef libcpp_vector[libcpp_utf8_string] v0 = tr_ids
        cdef _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition] * _r = new _MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition](self.inst.get().subset(v0))
        
        cdef MRMTransitionGroupCP py_result = MRMTransitionGroupCP.__new__(MRMTransitionGroupCP)
        py_result.inst = shared_ptr[_MRMTransitionGroup[_MSChromatogram,_ReactionMonitoringTransition]](_r)
        return py_result
    
    def isInternallyConsistent(self):
        """
        isInternallyConsistent(self) -> bool
        """
        cdef bool _r = self.inst.get().isInternallyConsistent()
        py_result = <bool>_r
        return py_result
    
    def chromatogramIdsMatch(self):
        """
        chromatogramIdsMatch(self) -> bool
        """
        cdef bool _r = self.inst.get().chromatogramIdsMatch()
        py_result = <bool>_r
        return py_result 

cdef class MzTabM:
    """
    Cython implementation of _MzTabM

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabM.html>`_

    Data model of MzTabM files
    
    Please see the official MzTabM specification at https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzTabM rv = MzTabM.__new__(MzTabM)
       rv.inst = shared_ptr[_MzTabM](new _MzTabM(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzTabM rv = MzTabM.__new__(MzTabM)
       rv.inst = shared_ptr[_MzTabM](new _MzTabM(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzTabM](new _MzTabM())
    
    def _init_1(self, MzTabM in_0 ):
        """
        _init_1(self, in_0: MzTabM ) -> None
        """
        assert isinstance(in_0, MzTabM), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzTabM](new _MzTabM((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzTabM ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzTabM)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def exportFeatureMapToMzTabM(self, FeatureMap feature_map ):
        """
        exportFeatureMapToMzTabM(self, feature_map: FeatureMap ) -> MzTabM
        Export FeatureMap with Identifications to MzTabM
        """
        assert isinstance(feature_map, FeatureMap), 'arg feature_map wrong type'
    
        cdef _MzTabM * _r = new _MzTabM(self.inst.get().exportFeatureMapToMzTabM((deref(feature_map.inst.get()))))
        cdef MzTabM py_result = MzTabM.__new__(MzTabM)
        py_result.inst = shared_ptr[_MzTabM](_r)
        return py_result 

cdef class PeptideIndexing:
    """
    Cython implementation of _PeptideIndexing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIndexing.html>`_
      -- Inherits from ['DefaultParamHandler']

    Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information
    
    All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy",
    depending on whether the protein accession contains the decoy pattern (parameter `decoy_string`) as a suffix or prefix, respectively (see parameter `prefix`).
    For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
    only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
    (For FDR calculations, "target+decoy" peptide hits count as target hits.)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideIndexing rv = PeptideIndexing.__new__(PeptideIndexing)
       rv.inst = shared_ptr[_PeptideIndexing](new _PeptideIndexing(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideIndexing rv = PeptideIndexing.__new__(PeptideIndexing)
       rv.inst = shared_ptr[_PeptideIndexing](new _PeptideIndexing(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideIndexing](new _PeptideIndexing())
    
    def _init_1(self, PeptideIndexing in_0 ):
        """
        _init_1(self, in_0: PeptideIndexing ) -> None
        """
        assert isinstance(in_0, PeptideIndexing), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideIndexing](new _PeptideIndexing((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideIndexing ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIndexing)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, list proteins , list prot_ids , PeptideIdentificationList pep_ids ):
        """
        run(self, proteins: List[FASTAEntry] , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> int
        """
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
        cdef libcpp_vector[_FASTAEntry] * v0 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item0
        for item0 in proteins:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in prot_ids:
            v1.push_back(deref(item1.inst.get()))
    
        cdef _PeptideIndexing_ExitCodes _r = self.inst.get().run(deref(v0), deref(v1), (deref(pep_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v1.begin()
        replace_0 = []
        while it_prot_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item1)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v1
        cdef libcpp_vector[_FASTAEntry].iterator it_proteins = v0.begin()
        replace_0 = []
        while it_proteins != v0.end():
            item0 = FASTAEntry.__new__(FASTAEntry)
            item0.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_proteins)))
            replace_0.append(item0)
            inc(it_proteins)
        proteins[:] = replace_0
        del v0
        py_result = <int>_r
        return py_result
    
    def getDecoyString(self):
        """
        getDecoyString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getDecoyString()
        py_result = convOutputString(_r)
        return py_result
    
    def isPrefix(self):
        """
        isPrefix(self) -> bool
        """
        cdef bool _r = self.inst.get().isPrefix()
        py_result = <bool>_r
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
    PeptideIndexing_ExitCodes = __PeptideIndexing_ExitCodes 

cdef class QTClusterFinder:
    """
    Cython implementation of _QTClusterFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QTClusterFinder.html>`_
      -- Inherits from ['BaseGroupFinder']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_QTClusterFinder](new _QTClusterFinder())
    
    def _run_0(self, list input_maps , ConsensusMap result_map ):
        """
        _run_0(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None
        """
        assert isinstance(input_maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in input_maps), 'arg input_maps wrong type'
        assert isinstance(result_map, ConsensusMap), 'arg result_map wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in input_maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().run(deref(v0), (deref(result_map.inst.get())))
        cdef libcpp_vector[_ConsensusMap].iterator it_input_maps = v0.begin()
        replace_0 = []
        while it_input_maps != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_input_maps)))
            replace_0.append(item0)
            inc(it_input_maps)
        input_maps[:] = replace_0
        del v0
    
    def _run_1(self, list input_maps , ConsensusMap result_map ):
        """
        _run_1(self, input_maps: List[FeatureMap] , result_map: ConsensusMap ) -> None
        """
        assert isinstance(input_maps, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in input_maps), 'arg input_maps wrong type'
        assert isinstance(result_map, ConsensusMap), 'arg result_map wrong type'
        cdef libcpp_vector[_FeatureMap] * v0 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item0
        for item0 in input_maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().run(deref(v0), (deref(result_map.inst.get())))
        cdef libcpp_vector[_FeatureMap].iterator it_input_maps = v0.begin()
        replace_0 = []
        while it_input_maps != v0.end():
            item0 = FeatureMap.__new__(FeatureMap)
            item0.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_input_maps)))
            replace_0.append(item0)
            inc(it_input_maps)
        input_maps[:] = replace_0
        del v0
    
    def run(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: run(self, input_maps: List[FeatureMap] , result_map: ConsensusMap ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in args[0])) and (isinstance(args[1], ConsensusMap)):
            return self._run_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in args[0])) and (isinstance(args[1], ConsensusMap)):
            return self._run_1(*args)
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

cdef class RNaseDB:
    """
    Cython implementation of _RNaseDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RNaseDB.html>`_
    """

    
    def getEnzyme(self,  name ):
        """
        getEnzyme(self, name: Union[bytes, str, String] ) -> DigestionEnzymeRNA
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef const _DigestionEnzymeRNA * __r = (self.inst.get().getEnzyme(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _DigestionEnzymeRNA * _r = new _DigestionEnzymeRNA(deref(__r))
        cdef DigestionEnzymeRNA py_result = DigestionEnzymeRNA.__new__(DigestionEnzymeRNA)
        py_result.inst = shared_ptr[_DigestionEnzymeRNA](_r)
        return py_result
    
    def getEnzymeByRegEx(self,  cleavage_regex ):
        """
        getEnzymeByRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> DigestionEnzymeRNA
        """
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
    
        cdef const _DigestionEnzymeRNA * __r = (self.inst.get().getEnzymeByRegEx(deref((convString(cleavage_regex)).get())))
        if __r == NULL:
            return None
        cdef _DigestionEnzymeRNA * _r = new _DigestionEnzymeRNA(deref(__r))
        cdef DigestionEnzymeRNA py_result = DigestionEnzymeRNA.__new__(DigestionEnzymeRNA)
        py_result.inst = shared_ptr[_DigestionEnzymeRNA](_r)
        return py_result
    
    def getAllNames(self, list all_names ):
        """
        getAllNames(self, all_names: List[bytes] ) -> None
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def hasEnzyme(self,  name ):
        """
        hasEnzyme(self, name: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasEnzyme(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_RNaseDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_RNaseDB](_getInstance_RNaseDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class Ribonucleotide:
    """
    Cython implementation of _Ribonucleotide

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Ribonucleotide_1_1Ribonucleotide.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).getName().c_str() )

    
    def __copy__(self):
       cdef Ribonucleotide rv = Ribonucleotide.__new__(Ribonucleotide)
       rv.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Ribonucleotide rv = Ribonucleotide.__new__(Ribonucleotide)
       rv.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide())
    
    def _init_1(self, Ribonucleotide in_0 ):
        """
        _init_1(self, in_0: Ribonucleotide ) -> None
        """
        assert isinstance(in_0, Ribonucleotide), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide((deref(in_0.inst.get()))))
    
    def _init_2(self,  name ,  code ,  new_code ,  html_code , EmpiricalFormula formula , bytes origin , double mono_mass , double avg_mass , int term_spec , EmpiricalFormula baseloss_formula ):
        """
        _init_2(self, name: Union[bytes, str, String] , code: Union[bytes, str, String] , new_code: Union[bytes, str, String] , html_code: Union[bytes, str, String] , formula: EmpiricalFormula , origin: bytes , mono_mass: float , avg_mass: float , term_spec: int , baseloss_formula: EmpiricalFormula ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(code, str) or isinstance(code, bytes) or isinstance(code, String)), 'arg code wrong type'
        assert (isinstance(new_code, str) or isinstance(new_code, bytes) or isinstance(new_code, String)), 'arg new_code wrong type'
        assert (isinstance(html_code, str) or isinstance(html_code, bytes) or isinstance(html_code, String)), 'arg html_code wrong type'
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
        assert isinstance(origin, bytes) and len(origin) == 1, 'arg origin wrong type'
        assert isinstance(mono_mass, float), 'arg mono_mass wrong type'
        assert isinstance(avg_mass, float), 'arg avg_mass wrong type'
        assert term_spec in [0, 1, 2, 3], 'arg term_spec wrong type'
        assert isinstance(baseloss_formula, EmpiricalFormula), 'arg baseloss_formula wrong type'
    
    
    
    
    
    
    
    
    
    
        self.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref((convString(name)).get()), deref((convString(code)).get()), deref((convString(new_code)).get()), deref((convString(html_code)).get()), (deref(formula.inst.get())), (<char>((origin)[0])), (<double>mono_mass), (<double>avg_mass), (<_TermSpecificityNuc>term_spec), (deref(baseloss_formula.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Ribonucleotide ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , code: Union[bytes, str, String] , new_code: Union[bytes, str, String] , html_code: Union[bytes, str, String] , formula: EmpiricalFormula , origin: bytes , mono_mass: float , avg_mass: float , term_spec: int , baseloss_formula: EmpiricalFormula ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Ribonucleotide)):
             self._init_1(*args)
        elif (len(args)==10) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and (isinstance(args[4], EmpiricalFormula)) and (isinstance(args[5], bytes) and len(args[5]) == 1) and (isinstance(args[6], float)) and (isinstance(args[7], float)) and (args[8] in [0, 1, 2, 3]) and (isinstance(args[9], EmpiricalFormula)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCode(self):
        """
        getCode(self) -> Union[bytes, str, String]
        Returns the short name
        """
        cdef _String _r = self.inst.get().getCode()
        py_result = convOutputString(_r)
        return py_result
    
    def setCode(self,  code ):
        """
        setCode(self, code: Union[bytes, str, String] ) -> None
        Sets the short name
        """
        assert (isinstance(code, str) or isinstance(code, bytes) or isinstance(code, String)), 'arg code wrong type'
    
        self.inst.get().setCode(deref((convString(code)).get()))
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the ribonucleotide
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the ribonucleotide
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setFormula(self, EmpiricalFormula formula ):
        """
        setFormula(self, formula: EmpiricalFormula ) -> None
        Sets empirical formula of the ribonucleotide (must be full, with N and C-terminus)
        """
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
    
        self.inst.get().setFormula((deref(formula.inst.get())))
    
    def getFormula(self):
        """
        getFormula(self) -> EmpiricalFormula
        Returns the empirical formula of the residue
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def setAvgMass(self, double avg_mass ):
        """
        setAvgMass(self, avg_mass: float ) -> None
        Sets average mass of the ribonucleotide
        """
        assert isinstance(avg_mass, float), 'arg avg_mass wrong type'
    
        self.inst.get().setAvgMass((<double>avg_mass))
    
    def getAvgMass(self):
        """
        getAvgMass(self) -> float
        Returns average mass of the ribonucleotide
        """
        cdef double _r = self.inst.get().getAvgMass()
        py_result = <double>_r
        return py_result
    
    def setMonoMass(self, double mono_mass ):
        """
        setMonoMass(self, mono_mass: float ) -> None
        Sets monoisotopic mass of the ribonucleotide
        """
        assert isinstance(mono_mass, float), 'arg mono_mass wrong type'
    
        self.inst.get().setMonoMass((<double>mono_mass))
    
    def getMonoMass(self):
        """
        getMonoMass(self) -> float
        Returns monoisotopic mass of the ribonucleotide
        """
        cdef double _r = self.inst.get().getMonoMass()
        py_result = <double>_r
        return py_result
    
    def getNewCode(self):
        """
        getNewCode(self) -> Union[bytes, str, String]
        Returns the new code
        """
        cdef _String _r = self.inst.get().getNewCode()
        py_result = convOutputString(_r)
        return py_result
    
    def setNewCode(self,  code ):
        """
        setNewCode(self, code: Union[bytes, str, String] ) -> None
        Sets the new code
        """
        assert (isinstance(code, str) or isinstance(code, bytes) or isinstance(code, String)), 'arg code wrong type'
    
        self.inst.get().setNewCode(deref((convString(code)).get()))
    
    def getOrigin(self):
        """
        getOrigin(self) -> bytes
        Returns the code of the unmodified base (e.g., "A", "C", ...)
        """
        cdef char  _r = self.inst.get().getOrigin()
        py_result = chr(<char>(_r))
        return py_result
    
    def setOrigin(self, bytes origin ):
        """
        setOrigin(self, origin: bytes ) -> None
        Sets the code of the unmodified base (e.g., "A", "C", ...)
        """
        assert isinstance(origin, bytes) and len(origin) == 1, 'arg origin wrong type'
    
        self.inst.get().setOrigin((<char>((origin)[0])))
    
    def setHTMLCode(self,  html_code ):
        """
        setHTMLCode(self, html_code: Union[bytes, str, String] ) -> None
        Sets the HTML (RNAMods) code
        """
        assert (isinstance(html_code, str) or isinstance(html_code, bytes) or isinstance(html_code, String)), 'arg html_code wrong type'
    
        self.inst.get().setHTMLCode(deref((convString(html_code)).get()))
    
    def getHTMLCode(self):
        """
        getHTMLCode(self) -> Union[bytes, str, String]
        Returns the HTML (RNAMods) code
        """
        cdef _String _r = self.inst.get().getHTMLCode()
        py_result = convOutputString(_r)
        return py_result
    
    def setTermSpecificity(self, int term_spec ):
        """
        setTermSpecificity(self, term_spec: int ) -> None
        Sets the terminal specificity
        """
        assert term_spec in [0, 1, 2, 3], 'arg term_spec wrong type'
    
        self.inst.get().setTermSpecificity((<_TermSpecificityNuc>term_spec))
    
    def getTermSpecificity(self):
        """
        getTermSpecificity(self) -> int
        Returns the terminal specificity
        """
        cdef _TermSpecificityNuc _r = self.inst.get().getTermSpecificity()
        py_result = <int>_r
        return py_result
    
    def getBaselossFormula(self):
        """
        getBaselossFormula(self) -> EmpiricalFormula
        Returns sum formula after loss of the nucleobase
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getBaselossFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def setBaselossFormula(self, EmpiricalFormula formula ):
        """
        setBaselossFormula(self, formula: EmpiricalFormula ) -> None
        Sets sum formula after loss of the nucleobase
        """
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
    
        self.inst.get().setBaselossFormula((deref(formula.inst.get())))
    
    def isModified(self):
        """
        isModified(self) -> bool
        True if the ribonucleotide is a modified one
        """
        cdef bool _r = self.inst.get().isModified()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Ribonucleotide):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Ribonucleotide other_casted = other
        cdef Ribonucleotide self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    TermSpecificityNuc = __TermSpecificityNuc 

cdef class TMTEighteenPlexQuantitationMethod:
    """
    Cython implementation of _TMTEighteenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTEighteenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TMTEighteenPlexQuantitationMethod rv = TMTEighteenPlexQuantitationMethod.__new__(TMTEighteenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTEighteenPlexQuantitationMethod](new _TMTEighteenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TMTEighteenPlexQuantitationMethod rv = TMTEighteenPlexQuantitationMethod.__new__(TMTEighteenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTEighteenPlexQuantitationMethod](new _TMTEighteenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TMTEighteenPlexQuantitationMethod](new _TMTEighteenPlexQuantitationMethod())
    
    def _init_1(self, TMTEighteenPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: TMTEighteenPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, TMTEighteenPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TMTEighteenPlexQuantitationMethod](new _TMTEighteenPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TMTEighteenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTEighteenPlexQuantitationMethod)):
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
    
    def getChannelInformation(self):
        """
        getChannelInformation(self) -> List[IsobaricChannelInformation]
        """
        _r = self.inst.get().getChannelInformation()
        py_result = []
        cdef libcpp_vector[_IsobaricChannelInformation].iterator it__r = _r.begin()
        cdef IsobaricChannelInformation item_py_result
        while it__r != _r.end():
           item_py_result = IsobaricChannelInformation.__new__(IsobaricChannelInformation)
           item_py_result.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getNumberOfChannels(self):
        """
        getNumberOfChannels(self) -> int
        """
        cdef size_t _r = self.inst.get().getNumberOfChannels()
        py_result = <size_t>_r
        return py_result
    
    def getIsotopeCorrectionMatrix(self):
        """
        getIsotopeCorrectionMatrix(self) -> MatrixDouble
        """
        cdef _Matrix[double] * _r = new _Matrix[double](self.inst.get().getIsotopeCorrectionMatrix())
        cdef MatrixDouble py_result = MatrixDouble.__new__(MatrixDouble)
        py_result.inst = shared_ptr[_Matrix[double]](_r)
        return py_result
    
    def getReferenceChannel(self):
        """
        getReferenceChannel(self) -> int
        """
        cdef size_t _r = self.inst.get().getReferenceChannel()
        py_result = <size_t>_r
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
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        Sets the name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get())) 

cdef class UnimodXMLFile:
    """
    Cython implementation of _UnimodXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1UnimodXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_UnimodXMLFile](new _UnimodXMLFile())
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
        return py_result 

cdef class VersionDetails:
    """
    Cython implementation of _VersionDetails

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1VersionDetails.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property version_major:
        def __set__(self,  version_major):
        
            self.inst.get().version_major = (<int>version_major)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().version_major
            py_result = <int>_r
            return py_result
    
    property version_minor:
        def __set__(self,  version_minor):
        
            self.inst.get().version_minor = (<int>version_minor)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().version_minor
            py_result = <int>_r
            return py_result
    
    property version_patch:
        def __set__(self,  version_patch):
        
            self.inst.get().version_patch = (<int>version_patch)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().version_patch
            py_result = <int>_r
            return py_result
    
    property pre_release_identifier:
        def __set__(self,  pre_release_identifier):
        
            self.inst.get().pre_release_identifier = deref((convString(pre_release_identifier)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().pre_release_identifier
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef VersionDetails rv = VersionDetails.__new__(VersionDetails)
       rv.inst = shared_ptr[_VersionDetails](new _VersionDetails(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef VersionDetails rv = VersionDetails.__new__(VersionDetails)
       rv.inst = shared_ptr[_VersionDetails](new _VersionDetails(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_VersionDetails](new _VersionDetails())
    
    def _init_1(self, VersionDetails in_0 ):
        """
        _init_1(self, in_0: VersionDetails ) -> None
        """
        assert isinstance(in_0, VersionDetails), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_VersionDetails](new _VersionDetails((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: VersionDetails ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], VersionDetails)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 0, 4):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, VersionDetails):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef VersionDetails other_casted = other
        cdef VersionDetails self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get())
        if op==4:
            return deref(self_casted.inst.get()) > deref(other_casted.inst.get())
    create = __static_VersionDetails_create 

cdef class VersionInfo:
    """
    Cython implementation of _VersionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1VersionInfo.html>`_
    """

    getBranch = __static_VersionInfo_getBranch
    getRevision = __static_VersionInfo_getRevision
    getTime = __static_VersionInfo_getTime
    getVersion = __static_VersionInfo_getVersion
    getVersionStruct = __static_VersionInfo_getVersionStruct 
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
