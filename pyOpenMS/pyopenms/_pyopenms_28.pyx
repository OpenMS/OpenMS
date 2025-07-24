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
def __static_TransformationModelLinear_getDefaultParameters(Param in_0 ):
    """
    __static_TransformationModelLinear_getDefaultParameters(in_0: Param ) -> None
    """
    assert isinstance(in_0, Param), 'arg in_0 wrong type'

    _getDefaultParameters_TransformationModelLinear((deref(in_0.inst.get()))) 

cdef class ITRAQ_TYPES:
    None
    FOURPLEX = 0
    EIGHTPLEX = 1
    TMT_SIXPLEX = 2
    SIZE_OF_ITRAQ_TYPES = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __IntensityThresholdCalculation:
    None
    MANUAL = 0
    AUTOMAXBYSTDEV = 1
    AUTOMAXBYPERCENT = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AbsoluteQuantitationStandardsFile:
    """
    Cython implementation of _AbsoluteQuantitationStandardsFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationStandardsFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AbsoluteQuantitationStandardsFile rv = AbsoluteQuantitationStandardsFile.__new__(AbsoluteQuantitationStandardsFile)
       rv.inst = shared_ptr[_AbsoluteQuantitationStandardsFile](new _AbsoluteQuantitationStandardsFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AbsoluteQuantitationStandardsFile rv = AbsoluteQuantitationStandardsFile.__new__(AbsoluteQuantitationStandardsFile)
       rv.inst = shared_ptr[_AbsoluteQuantitationStandardsFile](new _AbsoluteQuantitationStandardsFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AbsoluteQuantitationStandardsFile](new _AbsoluteQuantitationStandardsFile())
    
    def _init_1(self, AbsoluteQuantitationStandardsFile in_0 ):
        """
        _init_1(self, in_0: AbsoluteQuantitationStandardsFile ) -> None
        """
        assert isinstance(in_0, AbsoluteQuantitationStandardsFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AbsoluteQuantitationStandardsFile](new _AbsoluteQuantitationStandardsFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AbsoluteQuantitationStandardsFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AbsoluteQuantitationStandardsFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , list run_concentrations ):
        """
        load(self, filename: Union[bytes, str, String] , run_concentrations: List[AQS_runConcentration] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(run_concentrations, list) and all(isinstance(elemt_rec, AQS_runConcentration) for elemt_rec in run_concentrations), 'arg run_concentrations wrong type'
    
        cdef libcpp_vector[_AQS_runConcentration] * v1 = new libcpp_vector[_AQS_runConcentration]()
        cdef AQS_runConcentration item1
        for item1 in run_concentrations:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().load(deref((convString(filename)).get()), deref(v1))
        cdef libcpp_vector[_AQS_runConcentration].iterator it_run_concentrations = v1.begin()
        replace_0 = []
        while it_run_concentrations != v1.end():
            item1 = AQS_runConcentration.__new__(AQS_runConcentration)
            item1.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration(deref(it_run_concentrations)))
            replace_0.append(item1)
            inc(it_run_concentrations)
        run_concentrations[:] = replace_0
        del v1 

cdef class ChannelInfo:
    """
    Cython implementation of _ChannelInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChannelInfo.html>`_
    """

    
    property description:
        def __set__(self, bytes description):
        
            self.inst.get().description = (<libcpp_string>description)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().description
            py_result = <libcpp_string>_r
            return py_result
    
    property name:
        def __set__(self,  name):
        
            self.inst.get().name = (<int>name)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().name
            py_result = <int>_r
            return py_result
    
    property id:
        def __set__(self,  id):
        
            self.inst.get().id = (<int>id)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().id
            py_result = <int>_r
            return py_result
    
    property center:
        def __set__(self, double center):
        
            self.inst.get().center = (<double>center)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().center
            py_result = <double>_r
            return py_result
    
    property active:
        def __set__(self, bool active):
        
            self.inst.get().active = (<bool>active)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().active
            py_result = <bool>_r
            return py_result 

cdef class DBoundingBox2:
    """
    Cython implementation of _DBoundingBox2

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DBoundingBox2.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DBoundingBox2 rv = DBoundingBox2.__new__(DBoundingBox2)
       rv.inst = shared_ptr[_DBoundingBox2](new _DBoundingBox2(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DBoundingBox2 rv = DBoundingBox2.__new__(DBoundingBox2)
       rv.inst = shared_ptr[_DBoundingBox2](new _DBoundingBox2(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DBoundingBox2](new _DBoundingBox2())
    
    def _init_1(self, DBoundingBox2 in_0 ):
        """
        _init_1(self, in_0: DBoundingBox2 ) -> None
        """
        assert isinstance(in_0, DBoundingBox2), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DBoundingBox2](new _DBoundingBox2((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DBoundingBox2 ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DBoundingBox2)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def minPosition(self):
        """
        minPosition(self) -> Union[Sequence[int], Sequence[float]]
        """
        cdef _DPosition2 _r = self.inst.get().minPosition()
        py_result = [_r[0], _r[1]]
        return py_result
    
    def maxPosition(self):
        """
        maxPosition(self) -> Union[Sequence[int], Sequence[float]]
        """
        cdef _DPosition2 _r = self.inst.get().maxPosition()
        py_result = [_r[0], _r[1]]
        return py_result 

cdef class DigestionEnzymeProtein:
    """
    Cython implementation of _DigestionEnzymeProtein

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzymeProtein.html>`_
      -- Inherits from ['DigestionEnzyme']

    Representation of a digestion enzyme for proteins (protease)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DigestionEnzymeProtein rv = DigestionEnzymeProtein.__new__(DigestionEnzymeProtein)
       rv.inst = shared_ptr[_DigestionEnzymeProtein](new _DigestionEnzymeProtein(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DigestionEnzymeProtein rv = DigestionEnzymeProtein.__new__(DigestionEnzymeProtein)
       rv.inst = shared_ptr[_DigestionEnzymeProtein](new _DigestionEnzymeProtein(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DigestionEnzymeProtein](new _DigestionEnzymeProtein())
    
    def _init_1(self, DigestionEnzymeProtein in_0 ):
        """
        _init_1(self, in_0: DigestionEnzymeProtein ) -> None
        """
        assert isinstance(in_0, DigestionEnzymeProtein), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DigestionEnzymeProtein](new _DigestionEnzymeProtein((deref(in_0.inst.get()))))
    
    def _init_2(self,  name ,  cleavage_regex , set synonyms ,  regex_description , EmpiricalFormula n_term_gain , EmpiricalFormula c_term_gain ,  psi_id ,  xtandem_id ,  comet_id ,  omssa_id ):
        """
        _init_2(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] , n_term_gain: EmpiricalFormula , c_term_gain: EmpiricalFormula , psi_id: Union[bytes, str, String] , xtandem_id: Union[bytes, str, String] , comet_id: int , omssa_id: int ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
        assert isinstance(synonyms, set) and all(isinstance(i, bytes) for i in synonyms), 'arg synonyms wrong type'
        assert (isinstance(regex_description, str) or isinstance(regex_description, bytes) or isinstance(regex_description, String)), 'arg regex_description wrong type'
        assert isinstance(n_term_gain, EmpiricalFormula), 'arg n_term_gain wrong type'
        assert isinstance(c_term_gain, EmpiricalFormula), 'arg c_term_gain wrong type'
        assert (isinstance(psi_id, str) or isinstance(psi_id, bytes) or isinstance(psi_id, String)), 'arg psi_id wrong type'
        assert (isinstance(xtandem_id, str) or isinstance(xtandem_id, bytes) or isinstance(xtandem_id, String)), 'arg xtandem_id wrong type'
        assert isinstance(comet_id, int), 'arg comet_id wrong type'
        assert isinstance(omssa_id, int), 'arg omssa_id wrong type'
    
    
        cdef libcpp_set[_String] * v2 = new libcpp_set[_String]()
        cdef bytes item2
        for item2 in synonyms:
           v2.insert(_String(<char *>item2))
    
    
    
    
    
    
    
        self.inst = shared_ptr[_DigestionEnzymeProtein](new _DigestionEnzymeProtein(deref((convString(name)).get()), deref((convString(cleavage_regex)).get()), deref(v2), deref((convString(regex_description)).get()), (deref(n_term_gain.inst.get())), (deref(c_term_gain.inst.get())), deref((convString(psi_id)).get()), deref((convString(xtandem_id)).get()), (<unsigned int>comet_id), (<unsigned int>omssa_id)))
        del v2
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DigestionEnzymeProtein ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , cleavage_regex: Union[bytes, str, String] , synonyms: Set[bytes] , regex_description: Union[bytes, str, String] , n_term_gain: EmpiricalFormula , c_term_gain: EmpiricalFormula , psi_id: Union[bytes, str, String] , xtandem_id: Union[bytes, str, String] , comet_id: int , omssa_id: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DigestionEnzymeProtein)):
             self._init_1(*args)
        elif (len(args)==10) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], set) and all(isinstance(i, bytes) for i in args[2])) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and (isinstance(args[4], EmpiricalFormula)) and (isinstance(args[5], EmpiricalFormula)) and ((isinstance(args[6], str) or isinstance(args[6], bytes) or isinstance(args[6], String))) and ((isinstance(args[7], str) or isinstance(args[7], bytes) or isinstance(args[7], String))) and (isinstance(args[8], int)) and (isinstance(args[9], int)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setNTermGain(self, EmpiricalFormula value ):
        """
        setNTermGain(self, value: EmpiricalFormula ) -> None
        Sets the N-term gain
        """
        assert isinstance(value, EmpiricalFormula), 'arg value wrong type'
    
        self.inst.get().setNTermGain((deref(value.inst.get())))
    
    def setCTermGain(self, EmpiricalFormula value ):
        """
        setCTermGain(self, value: EmpiricalFormula ) -> None
        Sets the C-term gain
        """
        assert isinstance(value, EmpiricalFormula), 'arg value wrong type'
    
        self.inst.get().setCTermGain((deref(value.inst.get())))
    
    def getNTermGain(self):
        """
        getNTermGain(self) -> EmpiricalFormula
        Returns the N-term gain
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getNTermGain())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getCTermGain(self):
        """
        getCTermGain(self) -> EmpiricalFormula
        Returns the C-term gain
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getCTermGain())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def setPSIID(self,  value ):
        """
        setPSIID(self, value: Union[bytes, str, String] ) -> None
        Sets the PSI ID
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setPSIID(deref((convString(value)).get()))
    
    def getPSIID(self):
        """
        getPSIID(self) -> Union[bytes, str, String]
        Returns the PSI ID
        """
        cdef _String _r = self.inst.get().getPSIID()
        py_result = convOutputString(_r)
        return py_result
    
    def setXTandemID(self,  value ):
        """
        setXTandemID(self, value: Union[bytes, str, String] ) -> None
        Sets the X! Tandem enzyme ID
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setXTandemID(deref((convString(value)).get()))
    
    def getXTandemID(self):
        """
        getXTandemID(self) -> Union[bytes, str, String]
        Returns the X! Tandem enzyme ID
        """
        cdef _String _r = self.inst.get().getXTandemID()
        py_result = convOutputString(_r)
        return py_result
    
    def setCometID(self,  value ):
        """
        setCometID(self, value: int ) -> None
        Sets the Comet enzyme ID
        """
        assert isinstance(value, int), 'arg value wrong type'
    
        self.inst.get().setCometID((<int>value))
    
    def getCometID(self):
        """
        getCometID(self) -> int
        Returns the Comet enzyme ID
        """
        cdef int _r = self.inst.get().getCometID()
        py_result = <int>_r
        return py_result
    
    def setOMSSAID(self,  value ):
        """
        setOMSSAID(self, value: int ) -> None
        Sets the OMSSA enzyme ID
        """
        assert isinstance(value, int), 'arg value wrong type'
    
        self.inst.get().setOMSSAID((<int>value))
    
    def getOMSSAID(self):
        """
        getOMSSAID(self) -> int
        Returns the OMSSA enzyme ID
        """
        cdef int _r = self.inst.get().getOMSSAID()
        py_result = <int>_r
        return py_result
    
    def setMSGFID(self,  value ):
        """
        setMSGFID(self, value: int ) -> None
        Sets the MSGFPlus enzyme id
        """
        assert isinstance(value, int), 'arg value wrong type'
    
        self.inst.get().setMSGFID((<int>value))
    
    def getMSGFID(self):
        """
        getMSGFID(self) -> int
        Returns the MSGFPlus enzyme id
        """
        cdef int _r = self.inst.get().getMSGFID()
        py_result = <int>_r
        return py_result
    
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
        if not isinstance(other, DigestionEnzymeProtein):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DigestionEnzymeProtein other_casted = other
        cdef DigestionEnzymeProtein self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class FIAMSScheduler:
    """
    Cython implementation of _FIAMSScheduler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FIAMSScheduler.html>`_

      ADD PYTHON DOCUMENTATION HERE
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FIAMSScheduler rv = FIAMSScheduler.__new__(FIAMSScheduler)
       rv.inst = shared_ptr[_FIAMSScheduler](new _FIAMSScheduler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FIAMSScheduler rv = FIAMSScheduler.__new__(FIAMSScheduler)
       rv.inst = shared_ptr[_FIAMSScheduler](new _FIAMSScheduler(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor
        """
        self.inst = shared_ptr[_FIAMSScheduler](new _FIAMSScheduler())
    
    def _init_1(self, FIAMSScheduler in_0 ):
        """
        _init_1(self, in_0: FIAMSScheduler ) -> None
        """
        assert isinstance(in_0, FIAMSScheduler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FIAMSScheduler](new _FIAMSScheduler((deref(in_0.inst.get()))))
    
    def _init_2(self,  filename ,  base_dir , bool load_cached_ ):
        """
        _init_2(self, filename: Union[bytes, str, String] , base_dir: Union[bytes, str, String] , load_cached_: bool ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(base_dir, str) or isinstance(base_dir, bytes) or isinstance(base_dir, String)), 'arg base_dir wrong type'
        assert isinstance(load_cached_, pybool_t), 'arg load_cached_ wrong type'
    
    
    
        self.inst = shared_ptr[_FIAMSScheduler](new _FIAMSScheduler(deref((convString(filename)).get()), deref((convString(base_dir)).get()), (<bool>load_cached_)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Scheduler for FIA-MS data batches. Works with FIAMSDataProcessor

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FIAMSScheduler ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , base_dir: Union[bytes, str, String] , load_cached_: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FIAMSScheduler)):
             self._init_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], pybool_t)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self):
        """
        run(self) -> None
        Run the FIA-MS data analysis for the batch defined in the @filename_
        """
        self.inst.get().run()
    
    def getBaseDir(self):
        """
        getBaseDir(self) -> Union[bytes, str, String]
        Returns the base directory for the relevant paths from the csv file
        """
        cdef _String _r = self.inst.get().getBaseDir()
        py_result = convOutputString(_r)
        return py_result 

cdef class FeatureGroupingAlgorithmLabeled:
    """
    Cython implementation of _FeatureGroupingAlgorithmLabeled

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmLabeled.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureGroupingAlgorithmLabeled](new _FeatureGroupingAlgorithmLabeled())
    
    def group(self, list maps , ConsensusMap out ):
        """
        group(self, maps: List[FeatureMap] , out: ConsensusMap ) -> None
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

cdef class IsobaricQuantifier:
    """
    Cython implementation of _IsobaricQuantifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricQuantifier.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsobaricQuantifier rv = IsobaricQuantifier.__new__(IsobaricQuantifier)
       rv.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricQuantifier rv = IsobaricQuantifier.__new__(IsobaricQuantifier)
       rv.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier(deref(self.inst.get())))
       return rv
    
    def _init_0(self, IsobaricQuantifier in_0 ):
        """
        _init_0(self, in_0: IsobaricQuantifier ) -> None
        """
        assert isinstance(in_0, IsobaricQuantifier), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier((deref(in_0.inst.get()))))
    
    def _init_1(self, ItraqFourPlexQuantitationMethod quant_method ):
        """
        _init_1(self, quant_method: ItraqFourPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqFourPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier((quant_method.inst.get())))
    
    def _init_2(self, ItraqEightPlexQuantitationMethod quant_method ):
        """
        _init_2(self, quant_method: ItraqEightPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, ItraqEightPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier((quant_method.inst.get())))
    
    def _init_3(self, TMTSixPlexQuantitationMethod quant_method ):
        """
        _init_3(self, quant_method: TMTSixPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTSixPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier((quant_method.inst.get())))
    
    def _init_4(self, TMTTenPlexQuantitationMethod quant_method ):
        """
        _init_4(self, quant_method: TMTTenPlexQuantitationMethod ) -> None
        """
        assert isinstance(quant_method, TMTTenPlexQuantitationMethod), 'arg quant_method wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifier](new _IsobaricQuantifier((quant_method.inst.get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricQuantifier ) -> None
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
        if (len(args)==1) and (isinstance(args[0], IsobaricQuantifier)):
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
    
    def quantify(self, ConsensusMap consensus_map_in , ConsensusMap consensus_map_out ):
        """
        quantify(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap ) -> None
        """
        assert isinstance(consensus_map_in, ConsensusMap), 'arg consensus_map_in wrong type'
        assert isinstance(consensus_map_out, ConsensusMap), 'arg consensus_map_out wrong type'
    
    
        self.inst.get().quantify((deref(consensus_map_in.inst.get())), (deref(consensus_map_out.inst.get())))
    
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

cdef class ItraqConstants:
    """
    Cython implementation of _ItraqConstants

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqConstants.html>`_

    Some constants used throughout iTRAQ classes
    
    Constants for iTRAQ experiments and a ChannelInfo structure to store information about a single channel
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ItraqConstants rv = ItraqConstants.__new__(ItraqConstants)
       rv.inst = shared_ptr[_ItraqConstants](new _ItraqConstants(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ItraqConstants rv = ItraqConstants.__new__(ItraqConstants)
       rv.inst = shared_ptr[_ItraqConstants](new _ItraqConstants(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ItraqConstants](new _ItraqConstants())
    
    def _init_1(self, ItraqConstants in_0 ):
        """
        _init_1(self, in_0: ItraqConstants ) -> None
        """
        assert isinstance(in_0, ItraqConstants), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ItraqConstants](new _ItraqConstants((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ItraqConstants ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqConstants)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIsotopeMatrixAsStringList(self,  itraq_type , list isotope_corrections ):
        """
        getIsotopeMatrixAsStringList(self, itraq_type: int , isotope_corrections: List[MatrixDouble] ) -> List[bytes]
        Convert isotope correction matrix to stringlist\n
        
        Each line is converted into a string of the format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
        Useful for creating parameters or debug output
        
        
        :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
        :param isotope_corrections: Vector of the two matrices (4plex, 8plex)
        """
        assert isinstance(itraq_type, int), 'arg itraq_type wrong type'
        assert isinstance(isotope_corrections, list) and all(isinstance(elemt_rec, MatrixDouble) for elemt_rec in isotope_corrections), 'arg isotope_corrections wrong type'
    
        cdef libcpp_vector[_Matrix[double]] * v1 = new libcpp_vector[_Matrix[double]]()
        cdef MatrixDouble item1
        for item1 in isotope_corrections:
            v1.push_back(deref(item1.inst.get()))
        _r = self.inst.get().getIsotopeMatrixAsStringList((<int>itraq_type), deref(v1))
        cdef libcpp_vector[_Matrix[double]].iterator it_isotope_corrections = v1.begin()
        replace_0 = []
        while it_isotope_corrections != v1.end():
            item1 = MatrixDouble.__new__(MatrixDouble)
            item1.inst = shared_ptr[_Matrix[double]](new _Matrix[double](deref(it_isotope_corrections)))
            replace_0.append(item1)
            inc(it_isotope_corrections)
        isotope_corrections[:] = replace_0
        del v1
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def updateIsotopeMatrixFromStringList(self,  itraq_type , list channels , list isotope_corrections ):
        """
        updateIsotopeMatrixFromStringList(self, itraq_type: int , channels: List[bytes] , isotope_corrections: List[MatrixDouble] ) -> None
        Convert strings to isotope correction matrix rows\n
        
        Each string of format channel:-2Da/-1Da/+1Da/+2Da ; e.g. '114:0/0.3/4/0'
        is parsed and the corresponding channel(row) in the matrix is updated
        Not all channels need to be present, missing channels will be left untouched
        Useful to update the matrix with user isotope correction values
        
        
        :param itraq_type: Which matrix to stringify. Should be of values from enum ITRAQ_TYPES
        :param channels: New channel isotope values as strings
        :param isotope_corrections: Vector of the two matrices (4plex, 8plex)
        """
        assert isinstance(itraq_type, int), 'arg itraq_type wrong type'
        assert isinstance(channels, list) and all(isinstance(li, bytes) for li in channels), 'arg channels wrong type'
        assert isinstance(isotope_corrections, list) and all(isinstance(elemt_rec, MatrixDouble) for elemt_rec in isotope_corrections), 'arg isotope_corrections wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in channels:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[_Matrix[double]] * v2 = new libcpp_vector[_Matrix[double]]()
        cdef MatrixDouble item2
        for item2 in isotope_corrections:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().updateIsotopeMatrixFromStringList((<int>itraq_type), deref(v1), deref(v2))
        cdef libcpp_vector[_Matrix[double]].iterator it_isotope_corrections = v2.begin()
        replace_0 = []
        while it_isotope_corrections != v2.end():
            item2 = MatrixDouble.__new__(MatrixDouble)
            item2.inst = shared_ptr[_Matrix[double]](new _Matrix[double](deref(it_isotope_corrections)))
            replace_0.append(item2)
            inc(it_isotope_corrections)
        isotope_corrections[:] = replace_0
        del v2
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        channels[:] = replace
        del v1
    
    def translateIsotopeMatrix(self,  itraq_type , list isotope_corrections ):
        """
        translateIsotopeMatrix(self, itraq_type: int , isotope_corrections: List[MatrixDouble] ) -> MatrixDouble
        """
        assert isinstance(itraq_type, int), 'arg itraq_type wrong type'
        assert isinstance(isotope_corrections, list) and all(isinstance(elemt_rec, MatrixDouble) for elemt_rec in isotope_corrections), 'arg isotope_corrections wrong type'
    
        cdef libcpp_vector[_Matrix[double]] * v1 = new libcpp_vector[_Matrix[double]]()
        cdef MatrixDouble item1
        for item1 in isotope_corrections:
            v1.push_back(deref(item1.inst.get()))
        cdef _Matrix[double] * _r = new _Matrix[double](self.inst.get().translateIsotopeMatrix((<int &>itraq_type), deref(v1)))
        cdef libcpp_vector[_Matrix[double]].iterator it_isotope_corrections = v1.begin()
        replace_0 = []
        while it_isotope_corrections != v1.end():
            item1 = MatrixDouble.__new__(MatrixDouble)
            item1.inst = shared_ptr[_Matrix[double]](new _Matrix[double](deref(it_isotope_corrections)))
            replace_0.append(item1)
            inc(it_isotope_corrections)
        isotope_corrections[:] = replace_0
        del v1
        cdef MatrixDouble py_result = MatrixDouble.__new__(MatrixDouble)
        py_result.inst = shared_ptr[_Matrix[double]](_r)
        return py_result 

cdef class MRMFeature:
    """
    Cython implementation of _MRMFeature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeature.html>`_
      -- Inherits from ['Feature']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMFeature rv = MRMFeature.__new__(MRMFeature)
       rv.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeature rv = MRMFeature.__new__(MRMFeature)
       rv.inst = shared_ptr[_MRMFeature](new _MRMFeature(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeature](new _MRMFeature())
    
    def _init_1(self, MRMFeature in_0 ):
        """
        _init_1(self, in_0: MRMFeature ) -> None
        """
        assert isinstance(in_0, MRMFeature), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeature](new _MRMFeature((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeature ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeature)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getScores(self):
        """
        getScores(self) -> OpenSwath_Scores
        Returns all peakgroup scores
        """
        cdef _OpenSwath_Scores * _r = new _OpenSwath_Scores(self.inst.get().getScores())
        cdef OpenSwath_Scores py_result = OpenSwath_Scores.__new__(OpenSwath_Scores)
        py_result.inst = shared_ptr[_OpenSwath_Scores](_r)
        return py_result
    
    def setScores(self, OpenSwath_Scores s ):
        """
        setScores(self, s: OpenSwath_Scores ) -> None
        Sets all peakgroup scores
        """
        assert isinstance(s, OpenSwath_Scores), 'arg s wrong type'
    
        self.inst.get().setScores((deref(s.inst.get())))
    
    def getFeature(self,  key ):
        """
        getFeature(self, key: Union[bytes, str, String] ) -> Feature
        Returns a specified feature
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _Feature * _r = new _Feature(self.inst.get().getFeature(deref((convString(key)).get())))
        cdef Feature py_result = Feature.__new__(Feature)
        py_result.inst = shared_ptr[_Feature](_r)
        return py_result
    
    def addFeature(self, Feature f ,  key ):
        """
        addFeature(self, f: Feature , key: Union[bytes, str, String] ) -> None
        Adds an feature from a single chromatogram into the feature
        """
        assert isinstance(f, Feature), 'arg f wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addFeature((deref(f.inst.get())), deref((convString(key)).get()))
    
    def getFeatures(self):
        """
        getFeatures(self) -> List[Feature]
        Returns all the features
        """
        _r = self.inst.get().getFeatures()
        py_result = []
        cdef libcpp_vector[_Feature].iterator it__r = _r.begin()
        cdef Feature item_py_result
        while it__r != _r.end():
           item_py_result = Feature.__new__(Feature)
           item_py_result.inst = shared_ptr[_Feature](new _Feature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getFeatureIDs(self, list result ):
        """
        getFeatureIDs(self, result: List[bytes] ) -> None
        Returns a list of IDs of available features
        """
        assert isinstance(result, list) and all(isinstance(i, bytes) for i in result), 'arg result wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in result:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getFeatureIDs(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        result[:] = replace
        del v0
    
    def getPrecursorFeature(self,  key ):
        """
        getPrecursorFeature(self, key: Union[bytes, str, String] ) -> Feature
        Returns a specified precursor feature
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _Feature * _r = new _Feature(self.inst.get().getPrecursorFeature(deref((convString(key)).get())))
        cdef Feature py_result = Feature.__new__(Feature)
        py_result.inst = shared_ptr[_Feature](_r)
        return py_result
    
    def addPrecursorFeature(self, Feature f ,  key ):
        """
        addPrecursorFeature(self, f: Feature , key: Union[bytes, str, String] ) -> None
        Adds a precursor feature from a single chromatogram into the feature
        """
        assert isinstance(f, Feature), 'arg f wrong type'
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
    
        self.inst.get().addPrecursorFeature((deref(f.inst.get())), deref((convString(key)).get()))
    
    def getPrecursorFeatureIDs(self, list result ):
        """
        getPrecursorFeatureIDs(self, result: List[bytes] ) -> None
        Returns a list of IDs of available precursor features
        """
        assert isinstance(result, list) and all(isinstance(i, bytes) for i in result), 'arg result wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in result:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getPrecursorFeatureIDs(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        result[:] = replace
        del v0
    
    def getQuality(self,  index ):
        """
        getQuality(self, index: int ) -> float
        Returns the quality in dimension c
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef float _r = self.inst.get().getQuality((<size_t>index))
        py_result = <float>_r
        return py_result
    
    def setQuality(self,  index , float q ):
        """
        setQuality(self, index: int , q: float ) -> None
        Sets the quality in dimension c
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(q, float), 'arg q wrong type'
    
    
        self.inst.get().setQuality((<size_t>index), (<float>q))
    
    def getOverallQuality(self):
        """
        getOverallQuality(self) -> float
        Model and quality methods
        """
        cdef float _r = self.inst.get().getOverallQuality()
        py_result = <float>_r
        return py_result
    
    def setOverallQuality(self, float q ):
        """
        setOverallQuality(self, q: float ) -> None
        Sets the overall quality
        """
        assert isinstance(q, float), 'arg q wrong type'
    
        self.inst.get().setOverallQuality((<float>q))
    
    def getSubordinates(self):
        """
        getSubordinates(self) -> List[Feature]
        Returns the subordinate features
        """
        _r = self.inst.get().getSubordinates()
        py_result = []
        cdef libcpp_vector[_Feature].iterator it__r = _r.begin()
        cdef Feature item_py_result
        while it__r != _r.end():
           item_py_result = Feature.__new__(Feature)
           item_py_result.inst = shared_ptr[_Feature](new _Feature(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setSubordinates(self, list in_0 ):
        """
        setSubordinates(self, in_0: List[Feature] ) -> None
        Returns the subordinate features
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, Feature) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_Feature] * v0 = new libcpp_vector[_Feature]()
        cdef Feature item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSubordinates(deref(v0))
        del v0
    
    def encloses(self, double rt , double mz ):
        """
        encloses(self, rt: float , mz: float ) -> bool
        Returns if the mass trace convex hulls of the feature enclose the position specified by `rt` and `mz`
        
        
        :param rt: Sequence to digest
        :param mz: Digestion products
        """
        assert isinstance(rt, float), 'arg rt wrong type'
        assert isinstance(mz, float), 'arg mz wrong type'
    
    
        cdef bool _r = self.inst.get().encloses((<double>rt), (<double>mz))
        py_result = <bool>_r
        return py_result
    
    def getConvexHull(self):
        """
        getConvexHull(self) -> ConvexHull2D
        """
        cdef _ConvexHull2D * _r = new _ConvexHull2D(self.inst.get().getConvexHull())
        cdef ConvexHull2D py_result = ConvexHull2D.__new__(ConvexHull2D)
        py_result.inst = shared_ptr[_ConvexHull2D](_r)
        return py_result
    
    def getConvexHulls(self):
        """
        getConvexHulls(self) -> List[ConvexHull2D]
        """
        _r = self.inst.get().getConvexHulls()
        py_result = []
        cdef libcpp_vector[_ConvexHull2D].iterator it__r = _r.begin()
        cdef ConvexHull2D item_py_result
        while it__r != _r.end():
           item_py_result = ConvexHull2D.__new__(ConvexHull2D)
           item_py_result.inst = shared_ptr[_ConvexHull2D](new _ConvexHull2D(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setConvexHulls(self, list in_0 ):
        """
        setConvexHulls(self, in_0: List[ConvexHull2D] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, ConvexHull2D) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_ConvexHull2D] * v0 = new libcpp_vector[_ConvexHull2D]()
        cdef ConvexHull2D item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setConvexHulls(deref(v0))
        del v0
    
    def getWidth(self):
        """
        getWidth(self) -> float
        """
        cdef float _r = self.inst.get().getWidth()
        py_result = <float>_r
        return py_result
    
    def setWidth(self, float q ):
        """
        setWidth(self, q: float ) -> None
        """
        assert isinstance(q, float), 'arg q wrong type'
    
        self.inst.get().setWidth((<float>q))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def setCharge(self,  q ):
        """
        setCharge(self, q: int ) -> None
        """
        assert isinstance(q, int), 'arg q wrong type'
    
        self.inst.get().setCharge((<int>q))
    
    def getAnnotationState(self):
        """
        getAnnotationState(self) -> int
        """
        cdef _AnnotationState _r = self.inst.get().getAnnotationState()
        py_result = <int>_r
        return py_result
    
    def getPeptideIdentifications(self):
        """
        getPeptideIdentifications(self) -> PeptideIdentificationList
        Returns a reference to the PeptideIdentification vector
        """
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().getPeptideIdentifications())
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result
    
    def setPeptideIdentifications(self, PeptideIdentificationList peptides ):
        """
        setPeptideIdentifications(self, peptides: PeptideIdentificationList ) -> None
        Sets the PeptideIdentification vector
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().setPeptideIdentifications((deref(peptides.inst.get())))
    
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
        if not isinstance(other, MRMFeature):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MRMFeature other_casted = other
        cdef MRMFeature self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class MascotGenericFile:
    """
    Cython implementation of _MascotGenericFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MascotGenericFile.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MascotGenericFile rv = MascotGenericFile.__new__(MascotGenericFile)
       rv.inst = shared_ptr[_MascotGenericFile](new _MascotGenericFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MascotGenericFile rv = MascotGenericFile.__new__(MascotGenericFile)
       rv.inst = shared_ptr[_MascotGenericFile](new _MascotGenericFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MascotGenericFile](new _MascotGenericFile())
    
    def _init_1(self, MascotGenericFile in_0 ):
        """
        _init_1(self, in_0: MascotGenericFile ) -> None
        """
        assert isinstance(in_0, MascotGenericFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MascotGenericFile](new _MascotGenericFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MascotGenericFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MascotGenericFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename , MSExperiment experiment ):
        """
        store(self, filename: Union[bytes, str, String] , experiment: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(experiment, MSExperiment), 'arg experiment wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(experiment.inst.get())))
    
    def load(self,  filename , MSExperiment exp ):
        """
        load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None
        Loads a Mascot Generic File into a PeakMap
        
        
        :param filename: File name which the map should be read from
        :param exp: The map which is filled with the data from the given file
        :raises:
          Exception: FileNotFound is thrown if the given file could not be found
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(exp.inst.get())))
    
    def getHTTPPeakListEnclosure(self,  filename ):
        """
        getHTTPPeakListEnclosure(self, filename: Union[bytes, str, String] ) -> List[Union[bytes, str, String], Union[bytes, str, String]]
        Enclosing Strings of the peak list body for HTTP submission\n
        
        Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
        adding the peaks (in whatever format, e.g. mzXML) enclosed in this body
        The `filename` can later be found in the Mascot response
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        _r = self.inst.get().getHTTPPeakListEnclosure(deref((convString(filename)).get()))
        cdef String out1 = String.__new__(String)
        out1.inst = shared_ptr[_String](new _String(_r.first))
        cdef String out2 = String.__new__(String)
        out2.inst = shared_ptr[_String](new _String(_r.second))
        cdef list py_result = [out1, out2]
        return py_result
    
    def updateMembers_(self):
        """
        updateMembers_(self) -> None
        Docu in base class
        """
        self.inst.get().updateMembers_()
    
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

cdef class SeedListGenerator:
    """
    Cython implementation of _SeedListGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SeedListGenerator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SeedListGenerator rv = SeedListGenerator.__new__(SeedListGenerator)
       rv.inst = shared_ptr[_SeedListGenerator](new _SeedListGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SeedListGenerator rv = SeedListGenerator.__new__(SeedListGenerator)
       rv.inst = shared_ptr[_SeedListGenerator](new _SeedListGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SeedListGenerator](new _SeedListGenerator())
    
    def _init_1(self, SeedListGenerator in_0 ):
        """
        _init_1(self, in_0: SeedListGenerator ) -> None
        """
        assert isinstance(in_0, SeedListGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SeedListGenerator](new _SeedListGenerator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SeedListGenerator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SeedListGenerator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _generateSeedList_0(self, MSExperiment exp , np.ndarray[np.float32_t,ndim=2] seeds ):
        """
        _generateSeedList_0(self, exp: MSExperiment , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
        Generate a seed list based on an MS experiment
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert seeds.shape[1] == 2, 'arg seeds wrong type'
    
        cdef libcpp_vector[_DPosition2] _dp_vec_1
        cdef _DPosition2 _dp_1
        cdef int _dp_ii_1
        cdef int _dp_N_1 = seeds.shape[0]
        for _dp_ii_1 in range(_dp_N_1):
            _dp_1[0] = seeds[_dp_ii_1,0]
            _dp_1[1] = seeds[_dp_ii_1,1]
            _dp_vec_1.push_back(_dp_1)
        self.inst.get().generateSeedList((deref(exp.inst.get())), _dp_vec_1)
        _dp_n_1 = _dp_vec_1.size()
        seeds.resize((_dp_n_1,2))
        _dp_n_1 = _dp_vec_1.size()
        cdef libcpp_vector[_DPosition2].iterator _dp_it_1 = _dp_vec_1.begin()
        _dp_ii_1 = 0
        while _dp_it_1 != _dp_vec_1.end():
             seeds[_dp_ii_1, 0] = deref(_dp_it_1)[0]
             seeds[_dp_ii_1, 1] = deref(_dp_it_1)[1]
             inc(_dp_it_1)
             _dp_ii_1 += 1
    
    def _generateSeedList_1(self, PeptideIdentificationList peptides , np.ndarray[np.float32_t,ndim=2] seeds , bool use_peptide_mass ):
        """
        _generateSeedList_1(self, peptides: PeptideIdentificationList , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , use_peptide_mass: bool ) -> None
        Generates a seed list based on a list of peptide identifications
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert seeds.shape[1] == 2, 'arg seeds wrong type'
        assert isinstance(use_peptide_mass, pybool_t), 'arg use_peptide_mass wrong type'
    
        cdef libcpp_vector[_DPosition2] _dp_vec_1
        cdef _DPosition2 _dp_1
        cdef int _dp_ii_1
        cdef int _dp_N_1 = seeds.shape[0]
        for _dp_ii_1 in range(_dp_N_1):
            _dp_1[0] = seeds[_dp_ii_1,0]
            _dp_1[1] = seeds[_dp_ii_1,1]
            _dp_vec_1.push_back(_dp_1)
    
        self.inst.get().generateSeedList((deref(peptides.inst.get())), _dp_vec_1, (<bool>use_peptide_mass))
        _dp_n_1 = _dp_vec_1.size()
        seeds.resize((_dp_n_1,2))
        _dp_n_1 = _dp_vec_1.size()
        cdef libcpp_vector[_DPosition2].iterator _dp_it_1 = _dp_vec_1.begin()
        _dp_ii_1 = 0
        while _dp_it_1 != _dp_vec_1.end():
             seeds[_dp_ii_1, 0] = deref(_dp_it_1)[0]
             seeds[_dp_ii_1, 1] = deref(_dp_it_1)[1]
             inc(_dp_it_1)
             _dp_ii_1 += 1
    
    def generateSeedList(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: generateSeedList(self, exp: MSExperiment , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
          :noindex:
        
        Generate a seed list based on an MS experiment

        
        .. rubric:: Overload:
        .. py:function:: generateSeedList(self, peptides: PeptideIdentificationList , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , use_peptide_mass: bool ) -> None
          :noindex:
        
        Generates a seed list based on a list of peptide identifications
    
        """
        if (len(args)==2) and (isinstance(args[0], MSExperiment)) and (args[1].shape[1] == 2):
            return self._generateSeedList_0(*args)
        elif (len(args)==3) and (isinstance(args[0], PeptideIdentificationList)) and (args[1].shape[1] == 2) and (isinstance(args[2], pybool_t)):
            return self._generateSeedList_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _convertSeedList_0(self, np.ndarray[np.float32_t,ndim=2] seeds , FeatureMap features ):
        """
        _convertSeedList_0(self, seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , features: FeatureMap ) -> None
        Converts a list of seed positions to a feature map (expected format for FeatureFinder)
        """
        assert seeds.shape[1] == 2, 'arg seeds wrong type'
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        cdef libcpp_vector[_DPosition2] _dp_vec_0
        cdef _DPosition2 _dp_0
        cdef int _dp_ii_0
        cdef int _dp_N_0 = seeds.shape[0]
        for _dp_ii_0 in range(_dp_N_0):
            _dp_0[0] = seeds[_dp_ii_0,0]
            _dp_0[1] = seeds[_dp_ii_0,1]
            _dp_vec_0.push_back(_dp_0)
    
        self.inst.get().convertSeedList(_dp_vec_0, (deref(features.inst.get())))
        _dp_n_0 = _dp_vec_0.size()
        seeds.resize((_dp_n_0,2))
        _dp_n_0 = _dp_vec_0.size()
        cdef libcpp_vector[_DPosition2].iterator _dp_it_0 = _dp_vec_0.begin()
        _dp_ii_0 = 0
        while _dp_it_0 != _dp_vec_0.end():
             seeds[_dp_ii_0, 0] = deref(_dp_it_0)[0]
             seeds[_dp_ii_0, 1] = deref(_dp_it_0)[1]
             inc(_dp_it_0)
             _dp_ii_0 += 1
    
    def _convertSeedList_1(self, FeatureMap features , np.ndarray[np.float32_t,ndim=2] seeds ):
        """
        _convertSeedList_1(self, features: FeatureMap , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
        Converts a feature map with seed positions back to a simple list
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        assert seeds.shape[1] == 2, 'arg seeds wrong type'
    
        cdef libcpp_vector[_DPosition2] _dp_vec_1
        cdef _DPosition2 _dp_1
        cdef int _dp_ii_1
        cdef int _dp_N_1 = seeds.shape[0]
        for _dp_ii_1 in range(_dp_N_1):
            _dp_1[0] = seeds[_dp_ii_1,0]
            _dp_1[1] = seeds[_dp_ii_1,1]
            _dp_vec_1.push_back(_dp_1)
        self.inst.get().convertSeedList((deref(features.inst.get())), _dp_vec_1)
        _dp_n_1 = _dp_vec_1.size()
        seeds.resize((_dp_n_1,2))
        _dp_n_1 = _dp_vec_1.size()
        cdef libcpp_vector[_DPosition2].iterator _dp_it_1 = _dp_vec_1.begin()
        _dp_ii_1 = 0
        while _dp_it_1 != _dp_vec_1.end():
             seeds[_dp_ii_1, 0] = deref(_dp_it_1)[0]
             seeds[_dp_ii_1, 1] = deref(_dp_it_1)[1]
             inc(_dp_it_1)
             _dp_ii_1 += 1
    
    def convertSeedList(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: convertSeedList(self, seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' , features: FeatureMap ) -> None
          :noindex:
        
        Converts a list of seed positions to a feature map (expected format for FeatureFinder)

        
        .. rubric:: Overload:
        .. py:function:: convertSeedList(self, features: FeatureMap , seeds: '_np.ndarray[Any, _np.dtype[_np.float32]]' ) -> None
          :noindex:
        
        Converts a feature map with seed positions back to a simple list
    
        """
        if (len(args)==2) and (args[0].shape[1] == 2) and (isinstance(args[1], FeatureMap)):
            return self._convertSeedList_0(*args)
        elif (len(args)==2) and (isinstance(args[0], FeatureMap)) and (args[1].shape[1] == 2):
            return self._convertSeedList_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SignalToNoiseEstimatorMedian:
    """
    Cython implementation of _SignalToNoiseEstimatorMedian[_MSSpectrum]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMedian[_MSSpectrum].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SignalToNoiseEstimatorMedian rv = SignalToNoiseEstimatorMedian.__new__(SignalToNoiseEstimatorMedian)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSSpectrum]](new _SignalToNoiseEstimatorMedian[_MSSpectrum](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SignalToNoiseEstimatorMedian rv = SignalToNoiseEstimatorMedian.__new__(SignalToNoiseEstimatorMedian)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSSpectrum]](new _SignalToNoiseEstimatorMedian[_MSSpectrum](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSSpectrum]](new _SignalToNoiseEstimatorMedian[_MSSpectrum]())
    
    def _init_1(self, SignalToNoiseEstimatorMedian in_0 ):
        """
        _init_1(self, in_0: SignalToNoiseEstimatorMedian ) -> None
        """
        assert isinstance(in_0, SignalToNoiseEstimatorMedian), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSSpectrum]](new _SignalToNoiseEstimatorMedian[_MSSpectrum]((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SignalToNoiseEstimatorMedian ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SignalToNoiseEstimatorMedian)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def init(self, MSSpectrum spectrum ):
        """
        init(self, spectrum: MSSpectrum ) -> None
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().init((deref(spectrum.inst.get())))
    
    def getSignalToNoise(self,  index ):
        """
        getSignalToNoise(self, index: int ) -> float
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef double _r = self.inst.get().getSignalToNoise((<size_t>index))
        py_result = <double>_r
        return py_result
    
    def getSparseWindowPercent(self):
        """
        getSparseWindowPercent(self) -> float
        """
        cdef double _r = self.inst.get().getSparseWindowPercent()
        py_result = <double>_r
        return py_result
    
    def getHistogramRightmostPercent(self):
        """
        getHistogramRightmostPercent(self) -> float
        """
        cdef double _r = self.inst.get().getHistogramRightmostPercent()
        py_result = <double>_r
        return py_result
    
cdef class SignalToNoiseEstimatorMedianChrom:

    cdef shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram]] inst

    def __dealloc__(self):
         self.inst.reset()


    def init(self, MSChromatogram chromatogram ):
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'

        self.inst.get().init((deref(chromatogram.inst.get())))

    def _init_0(self):
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram]](new _SignalToNoiseEstimatorMedian[_MSChromatogram]())

    def _init_1(self, SignalToNoiseEstimatorMedianChrom in_0 ):
        assert isinstance(in_0, SignalToNoiseEstimatorMedianChrom), 'arg in_0 wrong type'

        self.inst = shared_ptr[_SignalToNoiseEstimatorMedian[_MSChromatogram]](new _SignalToNoiseEstimatorMedian[_MSChromatogram]((deref(in_0.inst.get()))))

    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SignalToNoiseEstimatorMedian)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))

    def getSignalToNoise(self, index ):
        assert isinstance(index, int), 'arg data_point wrong type'

        cdef double r = self.inst.get().getSignalToNoise(index)
        return r
    IntensityThresholdCalculation = __IntensityThresholdCalculation 

cdef class SiriusMSFile:
    """
    Cython implementation of _SiriusMSFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusMSFile rv = SiriusMSFile.__new__(SiriusMSFile)
       rv.inst = shared_ptr[_SiriusMSFile](new _SiriusMSFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusMSFile rv = SiriusMSFile.__new__(SiriusMSFile)
       rv.inst = shared_ptr[_SiriusMSFile](new _SiriusMSFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusMSFile](new _SiriusMSFile())
    
    def _init_1(self, SiriusMSFile in_0 ):
        """
        _init_1(self, in_0: SiriusMSFile ) -> None
        """
        assert isinstance(in_0, SiriusMSFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusMSFile](new _SiriusMSFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusMSFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusMSFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SiriusMSFile_AccessionInfo:
    """
    Cython implementation of _SiriusMSFile_AccessionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile_AccessionInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusMSFile_AccessionInfo rv = SiriusMSFile_AccessionInfo.__new__(SiriusMSFile_AccessionInfo)
       rv.inst = shared_ptr[_SiriusMSFile_AccessionInfo](new _SiriusMSFile_AccessionInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusMSFile_AccessionInfo rv = SiriusMSFile_AccessionInfo.__new__(SiriusMSFile_AccessionInfo)
       rv.inst = shared_ptr[_SiriusMSFile_AccessionInfo](new _SiriusMSFile_AccessionInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusMSFile_AccessionInfo](new _SiriusMSFile_AccessionInfo())
    
    def _init_1(self, SiriusMSFile_AccessionInfo in_0 ):
        """
        _init_1(self, in_0: SiriusMSFile_AccessionInfo ) -> None
        """
        assert isinstance(in_0, SiriusMSFile_AccessionInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusMSFile_AccessionInfo](new _SiriusMSFile_AccessionInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusMSFile_AccessionInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusMSFile_AccessionInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SiriusMSFile_CompoundInfo:
    """
    Cython implementation of _SiriusMSFile_CompoundInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusMSFile_CompoundInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusMSFile_CompoundInfo rv = SiriusMSFile_CompoundInfo.__new__(SiriusMSFile_CompoundInfo)
       rv.inst = shared_ptr[_SiriusMSFile_CompoundInfo](new _SiriusMSFile_CompoundInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusMSFile_CompoundInfo rv = SiriusMSFile_CompoundInfo.__new__(SiriusMSFile_CompoundInfo)
       rv.inst = shared_ptr[_SiriusMSFile_CompoundInfo](new _SiriusMSFile_CompoundInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusMSFile_CompoundInfo](new _SiriusMSFile_CompoundInfo())
    
    def _init_1(self, SiriusMSFile_CompoundInfo in_0 ):
        """
        _init_1(self, in_0: SiriusMSFile_CompoundInfo ) -> None
        """
        assert isinstance(in_0, SiriusMSFile_CompoundInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusMSFile_CompoundInfo](new _SiriusMSFile_CompoundInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusMSFile_CompoundInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusMSFile_CompoundInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SpectrumAccessOpenMSCached:
    """
    Cython implementation of _SpectrumAccessOpenMSCached

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMSCached.html>`_
      -- Inherits from ['ISpectrumAccess']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAccessOpenMSCached rv = SpectrumAccessOpenMSCached.__new__(SpectrumAccessOpenMSCached)
       rv.inst = shared_ptr[_SpectrumAccessOpenMSCached](new _SpectrumAccessOpenMSCached(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAccessOpenMSCached rv = SpectrumAccessOpenMSCached.__new__(SpectrumAccessOpenMSCached)
       rv.inst = shared_ptr[_SpectrumAccessOpenMSCached](new _SpectrumAccessOpenMSCached(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        pass
    
    def _init_1(self,  filename ):
        """
        _init_1(self, filename: Union[bytes, str, String] ) -> None
        An implementation of the Spectrum Access interface using on-disk caching
        
        This class implements the OpenSWATH Spectrum Access interface
        (ISpectrumAccess) using the CachedmzML class which is able to read and
        write a cached mzML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSCached](new _SpectrumAccessOpenMSCached(deref((convString(filename)).get())))
    
    def _init_2(self, SpectrumAccessOpenMSCached in_0 ):
        """
        _init_2(self, in_0: SpectrumAccessOpenMSCached ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMSCached), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSCached](new _SpectrumAccessOpenMSCached((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] ) -> None
          :noindex:
        
        An implementation of the Spectrum Access interface using on-disk caching
        
        This class implements the OpenSWATH Spectrum Access interface
        (ISpectrumAccess) using the CachedmzML class which is able to read and
        write a cached mzML file
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMSCached ) -> None
          :noindex:
    
        """
        if kwargs.get("__createUnsafeObject__") is True:
             self._init_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMSCached)):
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

cdef class TransformationModelLinear:
    """
    Cython implementation of _TransformationModelLinear

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelLinear.html>`_
      -- Inherits from ['TransformationModel']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list data , Param params ):
        """
        __init__(self, data: List[TM_DataPoint] , params: Param ) -> None
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        assert isinstance(params, Param), 'arg params wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_TransformationModelLinear](new _TransformationModelLinear(deref(v0), (deref(params.inst.get()))))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def evaluate(self, double value ):
        """
        evaluate(self, value: float ) -> float
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef double _r = self.inst.get().evaluate((<double>value))
        py_result = <double>_r
        return py_result
    
    def invert(self):
        """
        invert(self) -> None
        """
        self.inst.get().invert()
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def weightData(self, list data ):
        """
        weightData(self, data: List[TM_DataPoint] ) -> None
        Weight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().weightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkValidWeight(self,  weight , list valid_weights ):
        """
        checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool
        Check for a valid weighting function string
        """
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
        assert isinstance(valid_weights, list) and all(isinstance(i, bytes) for i in valid_weights), 'arg valid_weights wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in valid_weights:
           v1.push_back(_String(<char *>item1))
        cdef bool _r = self.inst.get().checkValidWeight(deref((convString(weight)).get()), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        valid_weights[:] = replace
        del v1
        py_result = <bool>_r
        return py_result
    
    def weightDatum(self, double datum ,  weight ):
        """
        weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Weight the data according to the weighting function
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().weightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def unWeightDatum(self, double datum ,  weight ):
        """
        unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Apply the reverse of the weighting function to the data
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().unWeightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def getValidXWeights(self):
        """
        getValidXWeights(self) -> List[bytes]
        Returns a list of valid x weight function stringss
        """
        _r = self.inst.get().getValidXWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getValidYWeights(self):
        """
        getValidYWeights(self) -> List[bytes]
        Returns a list of valid y weight function strings
        """
        _r = self.inst.get().getValidYWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def unWeightData(self, list data ):
        """
        unWeightData(self, data: List[TM_DataPoint] ) -> None
        Unweight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().unWeightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkDatumRange(self, double datum , double datum_min , double datum_max ):
        """
        checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float
        Check that the datum is within the valid min and max bounds
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert isinstance(datum_min, float), 'arg datum_min wrong type'
        assert isinstance(datum_max, float), 'arg datum_max wrong type'
    
    
    
        cdef double _r = self.inst.get().checkDatumRange((<const double &>datum), (<const double &>datum_min), (<const double &>datum_max))
        py_result = <double>_r
        return py_result
    getDefaultParameters = __static_TransformationModelLinear_getDefaultParameters 
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
