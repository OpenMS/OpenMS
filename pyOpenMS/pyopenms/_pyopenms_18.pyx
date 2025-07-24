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

cdef class __AcquisitionMode:
    None
    ACQMODENULL = 0
    PULSECOUNTING = 1
    ADC = 2
    TDC = 3
    TRANSIENTRECORDER = 4
    SIZE_OF_ACQUISITIONMODE = 5

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __Type_IonDetector:
    None
    TYPENULL = 0
    ELECTRONMULTIPLIER = 1
    PHOTOMULTIPLIER = 2
    FOCALPLANEARRAY = 3
    FARADAYCUP = 4
    CONVERSIONDYNODEELECTRONMULTIPLIER = 5
    CONVERSIONDYNODEPHOTOMULTIPLIER = 6
    MULTICOLLECTOR = 7
    CHANNELELECTRONMULTIPLIER = 8
    CHANNELTRON = 9
    DALYDETECTOR = 10
    MICROCHANNELPLATEDETECTOR = 11
    ARRAYDETECTOR = 12
    CONVERSIONDYNODE = 13
    DYNODE = 14
    FOCALPLANECOLLECTOR = 15
    IONTOPHOTONDETECTOR = 16
    POINTCOLLECTOR = 17
    POSTACCELERATIONDETECTOR = 18
    PHOTODIODEARRAYDETECTOR = 19
    INDUCTIVEDETECTOR = 20
    ELECTRONMULTIPLIERTUBE = 21
    SIZE_OF_TYPE = 22

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AbsoluteQuantitationMethod:
    """
    Cython implementation of _AbsoluteQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationMethod.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AbsoluteQuantitationMethod rv = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
       rv.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AbsoluteQuantitationMethod rv = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
       rv.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod())
    
    def _init_1(self, AbsoluteQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: AbsoluteQuantitationMethod ) -> None
        """
        assert isinstance(in_0, AbsoluteQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AbsoluteQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AbsoluteQuantitationMethod)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setLLOD(self, double llod ):
        """
        setLLOD(self, llod: float ) -> None
        """
        assert isinstance(llod, float), 'arg llod wrong type'
    
        self.inst.get().setLLOD((<double>llod))
    
    def setULOD(self, double ulod ):
        """
        setULOD(self, ulod: float ) -> None
        """
        assert isinstance(ulod, float), 'arg ulod wrong type'
    
        self.inst.get().setULOD((<double>ulod))
    
    def getLLOD(self):
        """
        getLLOD(self) -> float
        """
        cdef double _r = self.inst.get().getLLOD()
        py_result = <double>_r
        return py_result
    
    def getULOD(self):
        """
        getULOD(self) -> float
        """
        cdef double _r = self.inst.get().getULOD()
        py_result = <double>_r
        return py_result
    
    def setLLOQ(self, double lloq ):
        """
        setLLOQ(self, lloq: float ) -> None
        """
        assert isinstance(lloq, float), 'arg lloq wrong type'
    
        self.inst.get().setLLOQ((<double>lloq))
    
    def setULOQ(self, double uloq ):
        """
        setULOQ(self, uloq: float ) -> None
        """
        assert isinstance(uloq, float), 'arg uloq wrong type'
    
        self.inst.get().setULOQ((<double>uloq))
    
    def getLLOQ(self):
        """
        getLLOQ(self) -> float
        """
        cdef double _r = self.inst.get().getLLOQ()
        py_result = <double>_r
        return py_result
    
    def getULOQ(self):
        """
        getULOQ(self) -> float
        """
        cdef double _r = self.inst.get().getULOQ()
        py_result = <double>_r
        return py_result
    
    def checkLOD(self, double value ):
        """
        checkLOD(self, value: float ) -> bool
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().checkLOD((<double>value))
        py_result = <bool>_r
        return py_result
    
    def checkLOQ(self, double value ):
        """
        checkLOQ(self, value: float ) -> bool
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().checkLOQ((<double>value))
        py_result = <bool>_r
        return py_result
    
    def setComponentName(self,  component_name ):
        """
        setComponentName(self, component_name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(component_name, str) or isinstance(component_name, bytes) or isinstance(component_name, String)), 'arg component_name wrong type'
    
        self.inst.get().setComponentName(deref((convString(component_name)).get()))
    
    def setISName(self,  IS_name ):
        """
        setISName(self, IS_name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(IS_name, str) or isinstance(IS_name, bytes) or isinstance(IS_name, String)), 'arg IS_name wrong type'
    
        self.inst.get().setISName(deref((convString(IS_name)).get()))
    
    def setFeatureName(self,  feature_name ):
        """
        setFeatureName(self, feature_name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
    
        self.inst.get().setFeatureName(deref((convString(feature_name)).get()))
    
    def getComponentName(self):
        """
        getComponentName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getComponentName()
        py_result = convOutputString(_r)
        return py_result
    
    def getISName(self):
        """
        getISName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getISName()
        py_result = convOutputString(_r)
        return py_result
    
    def getFeatureName(self):
        """
        getFeatureName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getFeatureName()
        py_result = convOutputString(_r)
        return py_result
    
    def setConcentrationUnits(self,  concentration_units ):
        """
        setConcentrationUnits(self, concentration_units: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(concentration_units, str) or isinstance(concentration_units, bytes) or isinstance(concentration_units, String)), 'arg concentration_units wrong type'
    
        self.inst.get().setConcentrationUnits(deref((convString(concentration_units)).get()))
    
    def getConcentrationUnits(self):
        """
        getConcentrationUnits(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getConcentrationUnits()
        py_result = convOutputString(_r)
        return py_result
    
    def setTransformationModel(self,  transformation_model ):
        """
        setTransformationModel(self, transformation_model: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(transformation_model, str) or isinstance(transformation_model, bytes) or isinstance(transformation_model, String)), 'arg transformation_model wrong type'
    
        self.inst.get().setTransformationModel(deref((convString(transformation_model)).get()))
    
    def setTransformationModelParams(self, Param transformation_model_param ):
        """
        setTransformationModelParams(self, transformation_model_param: Param ) -> None
        """
        assert isinstance(transformation_model_param, Param), 'arg transformation_model_param wrong type'
    
        self.inst.get().setTransformationModelParams((deref(transformation_model_param.inst.get())))
    
    def getTransformationModel(self):
        """
        getTransformationModel(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getTransformationModel()
        py_result = convOutputString(_r)
        return py_result
    
    def getTransformationModelParams(self):
        """
        getTransformationModelParams(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getTransformationModelParams())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def setNPoints(self,  n_points ):
        """
        setNPoints(self, n_points: int ) -> None
        """
        assert isinstance(n_points, int), 'arg n_points wrong type'
    
        self.inst.get().setNPoints((<int>n_points))
    
    def setCorrelationCoefficient(self, double correlation_coefficient ):
        """
        setCorrelationCoefficient(self, correlation_coefficient: float ) -> None
        """
        assert isinstance(correlation_coefficient, float), 'arg correlation_coefficient wrong type'
    
        self.inst.get().setCorrelationCoefficient((<double>correlation_coefficient))
    
    def getNPoints(self):
        """
        getNPoints(self) -> int
        """
        cdef int _r = self.inst.get().getNPoints()
        py_result = <int>_r
        return py_result
    
    def getCorrelationCoefficient(self):
        """
        getCorrelationCoefficient(self) -> float
        """
        cdef double _r = self.inst.get().getCorrelationCoefficient()
        py_result = <double>_r
        return py_result 

cdef class ClusterProxyKD:
    """
    Cython implementation of _ClusterProxyKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ClusterProxyKD.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ClusterProxyKD rv = ClusterProxyKD.__new__(ClusterProxyKD)
       rv.inst = shared_ptr[_ClusterProxyKD](new _ClusterProxyKD(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ClusterProxyKD rv = ClusterProxyKD.__new__(ClusterProxyKD)
       rv.inst = shared_ptr[_ClusterProxyKD](new _ClusterProxyKD(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ClusterProxyKD](new _ClusterProxyKD())
    
    def _init_1(self, ClusterProxyKD in_0 ):
        """
        _init_1(self, in_0: ClusterProxyKD ) -> None
        """
        assert isinstance(in_0, ClusterProxyKD), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ClusterProxyKD](new _ClusterProxyKD((deref(in_0.inst.get()))))
    
    def _init_2(self,  size , double avg_distance ,  center_index ):
        """
        _init_2(self, size: int , avg_distance: float , center_index: int ) -> None
        """
        assert isinstance(size, int) and size >= 0, 'arg size wrong type'
        assert isinstance(avg_distance, float), 'arg avg_distance wrong type'
        assert isinstance(center_index, int) and center_index >= 0, 'arg center_index wrong type'
    
    
    
        self.inst = shared_ptr[_ClusterProxyKD](new _ClusterProxyKD((<size_t>size), (<double>avg_distance), (<size_t>center_index)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ClusterProxyKD ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, size: int , avg_distance: float , center_index: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ClusterProxyKD)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], float)) and (isinstance(args[2], int) and args[2] >= 0):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSize(self):
        """
        getSize(self) -> int
        """
        cdef size_t _r = self.inst.get().getSize()
        py_result = <size_t>_r
        return py_result
    
    def isValid(self):
        """
        isValid(self) -> bool
        """
        cdef bool _r = self.inst.get().isValid()
        py_result = <bool>_r
        return py_result
    
    def getAvgDistance(self):
        """
        getAvgDistance(self) -> float
        """
        cdef double _r = self.inst.get().getAvgDistance()
        py_result = <double>_r
        return py_result
    
    def getCenterIndex(self):
        """
        getCenterIndex(self) -> int
        """
        cdef size_t _r = self.inst.get().getCenterIndex()
        py_result = <size_t>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ClusterProxyKD):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ClusterProxyKD other_casted = other
        cdef ClusterProxyKD self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class ConsensusIDAlgorithmWorst:
    """
    Cython implementation of _ConsensusIDAlgorithmWorst

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmWorst.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmWorst](new _ConsensusIDAlgorithmWorst())
    
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

cdef class ContactPerson:
    """
    Cython implementation of _ContactPerson

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ContactPerson.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ContactPerson rv = ContactPerson.__new__(ContactPerson)
       rv.inst = shared_ptr[_ContactPerson](new _ContactPerson(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ContactPerson rv = ContactPerson.__new__(ContactPerson)
       rv.inst = shared_ptr[_ContactPerson](new _ContactPerson(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ContactPerson](new _ContactPerson())
    
    def _init_1(self, ContactPerson in_0 ):
        """
        _init_1(self, in_0: ContactPerson ) -> None
        """
        assert isinstance(in_0, ContactPerson), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ContactPerson](new _ContactPerson((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ContactPerson ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ContactPerson)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getFirstName(self):
        """
        getFirstName(self) -> Union[bytes, str, String]
        Returns the first name of the person
        """
        cdef _String _r = self.inst.get().getFirstName()
        py_result = convOutputString(_r)
        return py_result
    
    def setFirstName(self,  name ):
        """
        setFirstName(self, name: Union[bytes, str, String] ) -> None
        Sets the first name of the person
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setFirstName(deref((convString(name)).get()))
    
    def getLastName(self):
        """
        getLastName(self) -> Union[bytes, str, String]
        Returns the last name of the person
        """
        cdef _String _r = self.inst.get().getLastName()
        py_result = convOutputString(_r)
        return py_result
    
    def setLastName(self,  name ):
        """
        setLastName(self, name: Union[bytes, str, String] ) -> None
        Sets the last name of the person
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setLastName(deref((convString(name)).get()))
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the full name of the person (gets split into first and last name internally)
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getInstitution(self):
        """
        getInstitution(self) -> Union[bytes, str, String]
        Returns the affiliation
        """
        cdef _String _r = self.inst.get().getInstitution()
        py_result = convOutputString(_r)
        return py_result
    
    def setInstitution(self,  institution ):
        """
        setInstitution(self, institution: Union[bytes, str, String] ) -> None
        Sets the affiliation
        """
        assert (isinstance(institution, str) or isinstance(institution, bytes) or isinstance(institution, String)), 'arg institution wrong type'
    
        self.inst.get().setInstitution(deref((convString(institution)).get()))
    
    def getEmail(self):
        """
        getEmail(self) -> Union[bytes, str, String]
        Returns the email address
        """
        cdef _String _r = self.inst.get().getEmail()
        py_result = convOutputString(_r)
        return py_result
    
    def setEmail(self,  email ):
        """
        setEmail(self, email: Union[bytes, str, String] ) -> None
        Sets the email address
        """
        assert (isinstance(email, str) or isinstance(email, bytes) or isinstance(email, String)), 'arg email wrong type'
    
        self.inst.get().setEmail(deref((convString(email)).get()))
    
    def getURL(self):
        """
        getURL(self) -> Union[bytes, str, String]
        Returns the URL associated with the contact person (e.g., the institute webpage
        """
        cdef _String _r = self.inst.get().getURL()
        py_result = convOutputString(_r)
        return py_result
    
    def setURL(self,  email ):
        """
        setURL(self, email: Union[bytes, str, String] ) -> None
        Sets the URL associated with the contact person (e.g., the institute webpage
        """
        assert (isinstance(email, str) or isinstance(email, bytes) or isinstance(email, String)), 'arg email wrong type'
    
        self.inst.get().setURL(deref((convString(email)).get()))
    
    def getAddress(self):
        """
        getAddress(self) -> Union[bytes, str, String]
        Returns the address
        """
        cdef _String _r = self.inst.get().getAddress()
        py_result = convOutputString(_r)
        return py_result
    
    def setAddress(self,  email ):
        """
        setAddress(self, email: Union[bytes, str, String] ) -> None
        Sets the address
        """
        assert (isinstance(email, str) or isinstance(email, bytes) or isinstance(email, String)), 'arg email wrong type'
    
        self.inst.get().setAddress(deref((convString(email)).get()))
    
    def getContactInfo(self):
        """
        getContactInfo(self) -> Union[bytes, str, String]
        Returns miscellaneous info about the contact person
        """
        cdef _String _r = self.inst.get().getContactInfo()
        py_result = convOutputString(_r)
        return py_result
    
    def setContactInfo(self,  contact_info ):
        """
        setContactInfo(self, contact_info: Union[bytes, str, String] ) -> None
        Sets miscellaneous info about the contact person
        """
        assert (isinstance(contact_info, str) or isinstance(contact_info, bytes) or isinstance(contact_info, String)), 'arg contact_info wrong type'
    
        self.inst.get().setContactInfo(deref((convString(contact_info)).get()))
    
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
        if not isinstance(other, ContactPerson):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ContactPerson other_casted = other
        cdef ContactPerson self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class Element:
    """
    Cython implementation of _Element

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Element.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Element rv = Element.__new__(Element)
       rv.inst = shared_ptr[_Element](new _Element(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Element rv = Element.__new__(Element)
       rv.inst = shared_ptr[_Element](new _Element(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Element](new _Element())
    
    def _init_1(self, Element in_0 ):
        """
        _init_1(self, in_0: Element ) -> None
        """
        assert isinstance(in_0, Element), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Element](new _Element((deref(in_0.inst.get()))))
    
    def _init_2(self,  name ,  symbol ,  atomic_number , double average_weight , double mono_weight , IsotopeDistribution isotopes ):
        """
        _init_2(self, name: Union[bytes, str, String] , symbol: Union[bytes, str, String] , atomic_number: int , average_weight: float , mono_weight: float , isotopes: IsotopeDistribution ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(symbol, str) or isinstance(symbol, bytes) or isinstance(symbol, String)), 'arg symbol wrong type'
        assert isinstance(atomic_number, int), 'arg atomic_number wrong type'
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(mono_weight, float), 'arg mono_weight wrong type'
        assert isinstance(isotopes, IsotopeDistribution), 'arg isotopes wrong type'
    
    
    
    
    
    
        self.inst = shared_ptr[_Element](new _Element(deref((convString(name)).get()), deref((convString(symbol)).get()), (<unsigned int>atomic_number), (<double>average_weight), (<double>mono_weight), (deref(isotopes.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Element ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , symbol: Union[bytes, str, String] , atomic_number: int , average_weight: float , mono_weight: float , isotopes: IsotopeDistribution ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Element)):
             self._init_1(*args)
        elif (len(args)==6) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], int)) and (isinstance(args[3], float)) and (isinstance(args[4], float)) and (isinstance(args[5], IsotopeDistribution)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setAtomicNumber(self,  atomic_number ):
        """
        setAtomicNumber(self, atomic_number: int ) -> None
        Sets unique atomic number
        """
        assert isinstance(atomic_number, int), 'arg atomic_number wrong type'
    
        self.inst.get().setAtomicNumber((<unsigned int>atomic_number))
    
    def getAtomicNumber(self):
        """
        getAtomicNumber(self) -> int
        Returns the unique atomic number
        """
        cdef unsigned int _r = self.inst.get().getAtomicNumber()
        py_result = <unsigned int>_r
        return py_result
    
    def setAverageWeight(self, double weight ):
        """
        setAverageWeight(self, weight: float ) -> None
        Sets the average weight of the element
        """
        assert isinstance(weight, float), 'arg weight wrong type'
    
        self.inst.get().setAverageWeight((<double>weight))
    
    def getAverageWeight(self):
        """
        getAverageWeight(self) -> float
        Returns the average weight of the element
        """
        cdef double _r = self.inst.get().getAverageWeight()
        py_result = <double>_r
        return py_result
    
    def setMonoWeight(self, double weight ):
        """
        setMonoWeight(self, weight: float ) -> None
        Sets the mono isotopic weight of the element
        """
        assert isinstance(weight, float), 'arg weight wrong type'
    
        self.inst.get().setMonoWeight((<double>weight))
    
    def getMonoWeight(self):
        """
        getMonoWeight(self) -> float
        Returns the mono isotopic weight of the element
        """
        cdef double _r = self.inst.get().getMonoWeight()
        py_result = <double>_r
        return py_result
    
    def setIsotopeDistribution(self, IsotopeDistribution isotopes ):
        """
        setIsotopeDistribution(self, isotopes: IsotopeDistribution ) -> None
        Sets the isotope distribution of the element
        """
        assert isinstance(isotopes, IsotopeDistribution), 'arg isotopes wrong type'
    
        self.inst.get().setIsotopeDistribution((deref(isotopes.inst.get())))
    
    def getIsotopeDistribution(self):
        """
        getIsotopeDistribution(self) -> IsotopeDistribution
        Returns the isotope distribution of the element
        """
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().getIsotopeDistribution())
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the element
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the element
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setSymbol(self,  symbol ):
        """
        setSymbol(self, symbol: Union[bytes, str, String] ) -> None
        Sets symbol of the element
        """
        assert (isinstance(symbol, str) or isinstance(symbol, bytes) or isinstance(symbol, String)), 'arg symbol wrong type'
    
        self.inst.get().setSymbol(deref((convString(symbol)).get()))
    
    def getSymbol(self):
        """
        getSymbol(self) -> Union[bytes, str, String]
        Returns symbol of the element
        """
        cdef _String _r = self.inst.get().getSymbol()
        py_result = convOutputString(_r)
        return py_result 

cdef class IonDetector:
    """
    Cython implementation of _IonDetector

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IonDetector.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IonDetector rv = IonDetector.__new__(IonDetector)
       rv.inst = shared_ptr[_IonDetector](new _IonDetector(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IonDetector rv = IonDetector.__new__(IonDetector)
       rv.inst = shared_ptr[_IonDetector](new _IonDetector(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Description of a ion detector (part of a MS Instrument)
        """
        self.inst = shared_ptr[_IonDetector](new _IonDetector())
    
    def _init_1(self, IonDetector in_0 ):
        """
        _init_1(self, in_0: IonDetector ) -> None
        """
        assert isinstance(in_0, IonDetector), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IonDetector](new _IonDetector((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Description of a ion detector (part of a MS Instrument)

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IonDetector ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IonDetector)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getType(self):
        """
        getType(self) -> int
        Returns the detector type
        """
        cdef _Type_IonDetector _r = self.inst.get().getType()
        py_result = <int>_r
        return py_result
    
    def setType(self, int type_ ):
        """
        setType(self, type_: int ) -> None
        Sets the detector type
        """
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], 'arg type_ wrong type'
    
        self.inst.get().setType((<_Type_IonDetector>type_))
    
    def getAcquisitionMode(self):
        """
        getAcquisitionMode(self) -> int
        Returns the acquisition mode
        """
        cdef _AcquisitionMode _r = self.inst.get().getAcquisitionMode()
        py_result = <int>_r
        return py_result
    
    def setAcquisitionMode(self, int acquisition_mode ):
        """
        setAcquisitionMode(self, acquisition_mode: int ) -> None
        Sets the acquisition mode
        """
        assert acquisition_mode in [0, 1, 2, 3, 4, 5], 'arg acquisition_mode wrong type'
    
        self.inst.get().setAcquisitionMode((<_AcquisitionMode>acquisition_mode))
    
    def getResolution(self):
        """
        getResolution(self) -> float
        Returns the resolution (in ns)
        """
        cdef double _r = self.inst.get().getResolution()
        py_result = <double>_r
        return py_result
    
    def setResolution(self, double resolution ):
        """
        setResolution(self, resolution: float ) -> None
        Sets the resolution (in ns)
        """
        assert isinstance(resolution, float), 'arg resolution wrong type'
    
        self.inst.get().setResolution((<double>resolution))
    
    def getADCSamplingFrequency(self):
        """
        getADCSamplingFrequency(self) -> float
        Returns the analog-to-digital converter sampling frequency (in Hz)
        """
        cdef double _r = self.inst.get().getADCSamplingFrequency()
        py_result = <double>_r
        return py_result
    
    def setADCSamplingFrequency(self, double ADC_sampling_frequency ):
        """
        setADCSamplingFrequency(self, ADC_sampling_frequency: float ) -> None
        Sets the analog-to-digital converter sampling frequency (in Hz)
        """
        assert isinstance(ADC_sampling_frequency, float), 'arg ADC_sampling_frequency wrong type'
    
        self.inst.get().setADCSamplingFrequency((<double>ADC_sampling_frequency))
    
    def getOrder(self):
        """
        getOrder(self) -> int
        Returns the order
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
        if not isinstance(other, IonDetector):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IonDetector other_casted = other
        cdef IonDetector self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    AcquisitionMode = __AcquisitionMode
    Type_IonDetector = __Type_IonDetector 

cdef class MassDecomposition:
    """
    Cython implementation of _MassDecomposition

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MassDecomposition.html>`_

    Class represents a decomposition of a mass into amino acids
    
    This class represents a mass decomposition into amino acids. A
    decomposition are amino acids given with frequencies which add
    up to a specific mass.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MassDecomposition rv = MassDecomposition.__new__(MassDecomposition)
       rv.inst = shared_ptr[_MassDecomposition](new _MassDecomposition(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassDecomposition rv = MassDecomposition.__new__(MassDecomposition)
       rv.inst = shared_ptr[_MassDecomposition](new _MassDecomposition(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MassDecomposition](new _MassDecomposition())
    
    def _init_1(self, MassDecomposition in_0 ):
        """
        _init_1(self, in_0: MassDecomposition ) -> None
        """
        assert isinstance(in_0, MassDecomposition), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassDecomposition](new _MassDecomposition((deref(in_0.inst.get()))))
    
    def _init_2(self,  deco ):
        """
        _init_2(self, deco: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(deco, str) or isinstance(deco, bytes) or isinstance(deco, String)), 'arg deco wrong type'
    
        self.inst = shared_ptr[_MassDecomposition](new _MassDecomposition(deref((convString(deco)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassDecomposition ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, deco: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassDecomposition)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        Returns the decomposition as a string
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def toExpandedString(self):
        """
        toExpandedString(self) -> Union[bytes, str, String]
        Returns the decomposition as a string; instead of frequencies the amino acids are repeated
        """
        cdef _String _r = self.inst.get().toExpandedString()
        py_result = convOutputString(_r)
        return py_result
    
    def getNumberOfMaxAA(self):
        """
        getNumberOfMaxAA(self) -> int
        Returns the max frequency of this composition
        """
        cdef size_t _r = self.inst.get().getNumberOfMaxAA()
        py_result = <size_t>_r
        return py_result
    
    def containsTag(self,  tag ):
        """
        containsTag(self, tag: Union[bytes, str, String] ) -> bool
        Returns true if tag is contained in the mass decomposition
        """
        assert (isinstance(tag, str) or isinstance(tag, bytes) or isinstance(tag, String)), 'arg tag wrong type'
    
        cdef bool _r = self.inst.get().containsTag(deref((convString(tag)).get()))
        py_result = <bool>_r
        return py_result
    
    def compatible(self, MassDecomposition deco ):
        """
        compatible(self, deco: MassDecomposition ) -> bool
        Returns true if the mass decomposition if contained in this instance
        """
        assert isinstance(deco, MassDecomposition), 'arg deco wrong type'
    
        cdef bool _r = self.inst.get().compatible((deref(deco.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        Returns the decomposition as a string
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result 

cdef class ParamNode:
    """
    Cython implementation of _ParamNode

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Param_1_1ParamNode.html>`_
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
    
    property description:
        def __set__(self,  description):
        
            self.inst.get().description = deref((convString(description)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().description
            py_result = convOutputString(_r)
            return py_result
    
    property entries:
        def __set__(self, list entries):
            cdef libcpp_vector[_ParamEntry] * v0 = new libcpp_vector[_ParamEntry]()
            cdef ParamEntry item0
            for item0 in entries:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().entries = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().entries
            py_result = []
            cdef libcpp_vector[_ParamEntry].iterator it__r = _r.begin()
            cdef ParamEntry item_py_result
            while it__r != _r.end():
               item_py_result = ParamEntry.__new__(ParamEntry)
               item_py_result.inst = shared_ptr[_ParamEntry](new _ParamEntry(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    property nodes:
        def __set__(self, list nodes):
            cdef libcpp_vector[_ParamNode] * v0 = new libcpp_vector[_ParamNode]()
            cdef ParamNode item0
            for item0 in nodes:
                v0.push_back(deref(item0.inst.get()))
            self.inst.get().nodes = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().nodes
            py_result = []
            cdef libcpp_vector[_ParamNode].iterator it__r = _r.begin()
            cdef ParamNode item_py_result
            while it__r != _r.end():
               item_py_result = ParamNode.__new__(ParamNode)
               item_py_result.inst = shared_ptr[_ParamNode](new _ParamNode(deref(it__r)))
               py_result.append(item_py_result)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef ParamNode rv = ParamNode.__new__(ParamNode)
       rv.inst = shared_ptr[_ParamNode](new _ParamNode(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ParamNode rv = ParamNode.__new__(ParamNode)
       rv.inst = shared_ptr[_ParamNode](new _ParamNode(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ParamNode](new _ParamNode())
    
    def _init_1(self, ParamNode in_0 ):
        """
        _init_1(self, in_0: ParamNode ) -> None
        """
        assert isinstance(in_0, ParamNode), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamNode](new _ParamNode((deref(in_0.inst.get()))))
    
    def _init_2(self,  n ,  d ):
        """
        _init_2(self, n: Union[bytes, str, String] , d: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(n, str) or isinstance(n, bytes) or isinstance(n, String)), 'arg n wrong type'
        assert (isinstance(d, str) or isinstance(d, bytes) or isinstance(d, String)), 'arg d wrong type'
    
    
        self.inst = shared_ptr[_ParamNode](new _ParamNode(deref((convString(n)).get()), deref((convString(d)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ParamNode ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, n: Union[bytes, str, String] , d: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ParamNode)):
             self._init_1(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def findParentOf(self,  name ):
        """
        findParentOf(self, name: Union[bytes, str, String] ) -> ParamNode
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef  _ParamNode * __r = (self.inst.get().findParentOf(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _ParamNode * _r = new _ParamNode(deref(__r))
        cdef ParamNode py_result = ParamNode.__new__(ParamNode)
        py_result.inst = shared_ptr[_ParamNode](_r)
        return py_result
    
    def findEntryRecursive(self,  name ):
        """
        findEntryRecursive(self, name: Union[bytes, str, String] ) -> ParamEntry
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef  _ParamEntry * __r = (self.inst.get().findEntryRecursive(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _ParamEntry * _r = new _ParamEntry(deref(__r))
        cdef ParamEntry py_result = ParamEntry.__new__(ParamEntry)
        py_result.inst = shared_ptr[_ParamEntry](_r)
        return py_result
    
    def _insert_0(self, ParamNode node ,  prefix ):
        """
        _insert_0(self, node: ParamNode , prefix: Union[bytes, str, String] ) -> None
        """
        assert isinstance(node, ParamNode), 'arg node wrong type'
        assert (isinstance(prefix, str) or isinstance(prefix, bytes) or isinstance(prefix, String)), 'arg prefix wrong type'
    
    
        self.inst.get().insert((deref(node.inst.get())), deref((convString(prefix)).get()))
    
    def _insert_1(self, ParamEntry entry ,  prefix ):
        """
        _insert_1(self, entry: ParamEntry , prefix: Union[bytes, str, String] ) -> None
        """
        assert isinstance(entry, ParamEntry), 'arg entry wrong type'
        assert (isinstance(prefix, str) or isinstance(prefix, bytes) or isinstance(prefix, String)), 'arg prefix wrong type'
    
    
        self.inst.get().insert((deref(entry.inst.get())), deref((convString(prefix)).get()))
    
    def insert(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: insert(self, node: ParamNode , prefix: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: insert(self, entry: ParamEntry , prefix: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], ParamNode)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._insert_0(*args)
        elif (len(args)==2) and (isinstance(args[0], ParamEntry)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._insert_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def suffix(self,  key ):
        """
        suffix(self, key: Union[bytes, str, String] ) -> Union[bytes, str, String]
        """
        assert (isinstance(key, str) or isinstance(key, bytes) or isinstance(key, String)), 'arg key wrong type'
    
        cdef _String _r = self.inst.get().suffix(deref((convString(key)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ParamNode):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ParamNode other_casted = other
        cdef ParamNode self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class Peak1D:
    """
    Cython implementation of _Peak1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Peak1D.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Peak1D rv = Peak1D.__new__(Peak1D)
       rv.inst = shared_ptr[_Peak1D](new _Peak1D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Peak1D rv = Peak1D.__new__(Peak1D)
       rv.inst = shared_ptr[_Peak1D](new _Peak1D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Peak1D](new _Peak1D())
    
    def _init_1(self, Peak1D in_0 ):
        """
        _init_1(self, in_0: Peak1D ) -> None
        """
        assert isinstance(in_0, Peak1D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Peak1D](new _Peak1D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Peak1D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Peak1D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Peak1D):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Peak1D other_casted = other
        cdef Peak1D self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class RangeIntensity:
    """
    Cython implementation of _RangeIntensity

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeIntensity.html>`_
      -- Inherits from ['RangeBase']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RangeIntensity rv = RangeIntensity.__new__(RangeIntensity)
       rv.inst = shared_ptr[_RangeIntensity](new _RangeIntensity(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RangeIntensity rv = RangeIntensity.__new__(RangeIntensity)
       rv.inst = shared_ptr[_RangeIntensity](new _RangeIntensity(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RangeIntensity](new _RangeIntensity())
    
    def _init_1(self, RangeIntensity in_0 ):
        """
        _init_1(self, in_0: RangeIntensity ) -> None
        """
        assert isinstance(in_0, RangeIntensity), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RangeIntensity](new _RangeIntensity((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RangeIntensity ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RangeIntensity)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMinIntensity(self, double min ):
        """
        setMinIntensity(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMinIntensity((<double>min))
    
    def setMaxIntensity(self, double max ):
        """
        setMaxIntensity(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMaxIntensity((<double>max))
    
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
    
    def extendIntensity(self, double value ):
        """
        extendIntensity(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extendIntensity((<double>value))
    
    def containsIntensity(self, double value ):
        """
        containsIntensity(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().containsIntensity((<double>value))
        py_result = <bool>_r
        return py_result
    
    def setMin(self, double min ):
        """
        setMin(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMin((<double>min))
    
    def setMax(self, double max ):
        """
        setMax(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMax((<double>max))
    
    def getMin(self):
        """
        getMin(self) -> float
        """
        cdef double _r = self.inst.get().getMin()
        py_result = <double>_r
        return py_result
    
    def getMax(self):
        """
        getMax(self) -> float
        """
        cdef double _r = self.inst.get().getMax()
        py_result = <double>_r
        return py_result
    
    def extend(self, double value ):
        """
        extend(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extend((<double>value))
    
    def contains(self, double value ):
        """
        contains(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().contains((<double>value))
        py_result = <bool>_r
        return py_result 

cdef class RangeMZ:
    """
    Cython implementation of _RangeMZ

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeMZ.html>`_
      -- Inherits from ['RangeBase']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RangeMZ rv = RangeMZ.__new__(RangeMZ)
       rv.inst = shared_ptr[_RangeMZ](new _RangeMZ(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RangeMZ rv = RangeMZ.__new__(RangeMZ)
       rv.inst = shared_ptr[_RangeMZ](new _RangeMZ(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RangeMZ](new _RangeMZ())
    
    def _init_1(self, RangeMZ in_0 ):
        """
        _init_1(self, in_0: RangeMZ ) -> None
        """
        assert isinstance(in_0, RangeMZ), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RangeMZ](new _RangeMZ((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RangeMZ ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RangeMZ)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMinMZ(self, double min ):
        """
        setMinMZ(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMinMZ((<double>min))
    
    def setMaxMZ(self, double max ):
        """
        setMaxMZ(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMaxMZ((<double>max))
    
    def getMinMZ(self):
        """
        getMinMZ(self) -> float
        """
        cdef double _r = self.inst.get().getMinMZ()
        py_result = <double>_r
        return py_result
    
    def getMaxMZ(self):
        """
        getMaxMZ(self) -> float
        """
        cdef double _r = self.inst.get().getMaxMZ()
        py_result = <double>_r
        return py_result
    
    def extendMZ(self, double value ):
        """
        extendMZ(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extendMZ((<double>value))
    
    def containsMZ(self, double value ):
        """
        containsMZ(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().containsMZ((<double>value))
        py_result = <bool>_r
        return py_result
    
    def setMin(self, double min ):
        """
        setMin(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMin((<double>min))
    
    def setMax(self, double max ):
        """
        setMax(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMax((<double>max))
    
    def getMin(self):
        """
        getMin(self) -> float
        """
        cdef double _r = self.inst.get().getMin()
        py_result = <double>_r
        return py_result
    
    def getMax(self):
        """
        getMax(self) -> float
        """
        cdef double _r = self.inst.get().getMax()
        py_result = <double>_r
        return py_result
    
    def extend(self, double value ):
        """
        extend(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extend((<double>value))
    
    def contains(self, double value ):
        """
        contains(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().contains((<double>value))
        py_result = <bool>_r
        return py_result 

cdef class RangeMobility:
    """
    Cython implementation of _RangeMobility

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeMobility.html>`_
      -- Inherits from ['RangeBase']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RangeMobility rv = RangeMobility.__new__(RangeMobility)
       rv.inst = shared_ptr[_RangeMobility](new _RangeMobility(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RangeMobility rv = RangeMobility.__new__(RangeMobility)
       rv.inst = shared_ptr[_RangeMobility](new _RangeMobility(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RangeMobility](new _RangeMobility())
    
    def _init_1(self, RangeMobility in_0 ):
        """
        _init_1(self, in_0: RangeMobility ) -> None
        """
        assert isinstance(in_0, RangeMobility), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RangeMobility](new _RangeMobility((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RangeMobility ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RangeMobility)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMinMobility(self, double min ):
        """
        setMinMobility(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMinMobility((<double>min))
    
    def setMaxMobility(self, double max ):
        """
        setMaxMobility(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMaxMobility((<double>max))
    
    def getMinMobility(self):
        """
        getMinMobility(self) -> float
        Returns the minimum Mobility
        """
        cdef double _r = self.inst.get().getMinMobility()
        py_result = <double>_r
        return py_result
    
    def getMaxMobility(self):
        """
        getMaxMobility(self) -> float
        Returns the maximum Mobility
        """
        cdef double _r = self.inst.get().getMaxMobility()
        py_result = <double>_r
        return py_result
    
    def extendMobility(self, double value ):
        """
        extendMobility(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extendMobility((<double>value))
    
    def containsMobility(self, double value ):
        """
        containsMobility(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().containsMobility((<double>value))
        py_result = <bool>_r
        return py_result
    
    def setMin(self, double min ):
        """
        setMin(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMin((<double>min))
    
    def setMax(self, double max ):
        """
        setMax(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMax((<double>max))
    
    def getMin(self):
        """
        getMin(self) -> float
        """
        cdef double _r = self.inst.get().getMin()
        py_result = <double>_r
        return py_result
    
    def getMax(self):
        """
        getMax(self) -> float
        """
        cdef double _r = self.inst.get().getMax()
        py_result = <double>_r
        return py_result
    
    def extend(self, double value ):
        """
        extend(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extend((<double>value))
    
    def contains(self, double value ):
        """
        contains(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().contains((<double>value))
        py_result = <bool>_r
        return py_result 

cdef class RangeRT:
    """
    Cython implementation of _RangeRT

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RangeRT.html>`_
      -- Inherits from ['RangeBase']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RangeRT rv = RangeRT.__new__(RangeRT)
       rv.inst = shared_ptr[_RangeRT](new _RangeRT(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RangeRT rv = RangeRT.__new__(RangeRT)
       rv.inst = shared_ptr[_RangeRT](new _RangeRT(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RangeRT](new _RangeRT())
    
    def _init_1(self, RangeRT in_0 ):
        """
        _init_1(self, in_0: RangeRT ) -> None
        """
        assert isinstance(in_0, RangeRT), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RangeRT](new _RangeRT((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RangeRT ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RangeRT)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMinRT(self, double min ):
        """
        setMinRT(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMinRT((<double>min))
    
    def setMaxRT(self, double max ):
        """
        setMaxRT(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMaxRT((<double>max))
    
    def getMinRT(self):
        """
        getMinRT(self) -> float
        """
        cdef double _r = self.inst.get().getMinRT()
        py_result = <double>_r
        return py_result
    
    def getMaxRT(self):
        """
        getMaxRT(self) -> float
        """
        cdef double _r = self.inst.get().getMaxRT()
        py_result = <double>_r
        return py_result
    
    def extendRT(self, double value ):
        """
        extendRT(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extendRT((<double>value))
    
    def containsRT(self, double value ):
        """
        containsRT(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().containsRT((<double>value))
        py_result = <bool>_r
        return py_result
    
    def setMin(self, double min ):
        """
        setMin(self, min: float ) -> None
        """
        assert isinstance(min, float), 'arg min wrong type'
    
        self.inst.get().setMin((<double>min))
    
    def setMax(self, double max ):
        """
        setMax(self, max: float ) -> None
        """
        assert isinstance(max, float), 'arg max wrong type'
    
        self.inst.get().setMax((<double>max))
    
    def getMin(self):
        """
        getMin(self) -> float
        """
        cdef double _r = self.inst.get().getMin()
        py_result = <double>_r
        return py_result
    
    def getMax(self):
        """
        getMax(self) -> float
        """
        cdef double _r = self.inst.get().getMax()
        py_result = <double>_r
        return py_result
    
    def extend(self, double value ):
        """
        extend(self, value: float ) -> None
        Extend the range such that it includes the given @p value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().extend((<double>value))
    
    def contains(self, double value ):
        """
        contains(self, value: float ) -> bool
        Is value within [min, max]?
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef bool _r = self.inst.get().contains((<double>value))
        py_result = <bool>_r
        return py_result 

cdef class SpectraMerger:
    """
    Cython implementation of _SpectraMerger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectraMerger.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectraMerger rv = SpectraMerger.__new__(SpectraMerger)
       rv.inst = shared_ptr[_SpectraMerger](new _SpectraMerger(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectraMerger rv = SpectraMerger.__new__(SpectraMerger)
       rv.inst = shared_ptr[_SpectraMerger](new _SpectraMerger(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Merges blocks of MS or MS2 spectra
        """
        self.inst = shared_ptr[_SpectraMerger](new _SpectraMerger())
    
    def _init_1(self, SpectraMerger in_0 ):
        """
        _init_1(self, in_0: SpectraMerger ) -> None
        """
        assert isinstance(in_0, SpectraMerger), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectraMerger](new _SpectraMerger((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Merges blocks of MS or MS2 spectra

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectraMerger ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectraMerger)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def mergeSpectraBlockWise(self, MSExperiment exp ):
        """
        mergeSpectraBlockWise(self, exp: MSExperiment ) -> None
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().mergeSpectraBlockWise((deref(exp.inst.get())))
    
    def mergeSpectraPrecursors(self, MSExperiment exp ):
        """
        mergeSpectraPrecursors(self, exp: MSExperiment ) -> None
        Merges spectra with similar precursors (must have MS2 level)
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().mergeSpectraPrecursors((deref(exp.inst.get())))
    
    def average(self, MSExperiment exp ,  average_type ):
        """
        average(self, exp: MSExperiment , average_type: Union[bytes, str, String] ) -> None
        Average over neighbouring spectra
        
        :param exp: Experimental data to be averaged
        :param average_type: Averaging type to be used ("gaussian" or "tophat")
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(average_type, str) or isinstance(average_type, bytes) or isinstance(average_type, String)), 'arg average_type wrong type'
    
    
        self.inst.get().average((deref(exp.inst.get())), deref((convString(average_type)).get()))
    
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

cdef class SwathWindowLoader:
    """
    Cython implementation of _SwathWindowLoader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SwathWindowLoader.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SwathWindowLoader rv = SwathWindowLoader.__new__(SwathWindowLoader)
       rv.inst = shared_ptr[_SwathWindowLoader](new _SwathWindowLoader(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SwathWindowLoader rv = SwathWindowLoader.__new__(SwathWindowLoader)
       rv.inst = shared_ptr[_SwathWindowLoader](new _SwathWindowLoader(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SwathWindowLoader](new _SwathWindowLoader())
    
    def _init_1(self, SwathWindowLoader in_0 ):
        """
        _init_1(self, in_0: SwathWindowLoader ) -> None
        """
        assert isinstance(in_0, SwathWindowLoader), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SwathWindowLoader](new _SwathWindowLoader((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SwathWindowLoader ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SwathWindowLoader)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def annotateSwathMapsFromFile(self,  filename , list swath_maps , bool do_sort , bool force ):
        """
        annotateSwathMapsFromFile(self, filename: Union[bytes, str, String] , swath_maps: List[SwathMap] , do_sort: bool , force: bool ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(swath_maps, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in swath_maps), 'arg swath_maps wrong type'
        assert isinstance(do_sort, pybool_t), 'arg do_sort wrong type'
        assert isinstance(force, pybool_t), 'arg force wrong type'
    
        cdef libcpp_vector[_SwathMap] * v1 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item1
        for item1 in swath_maps:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().annotateSwathMapsFromFile(deref((convString(filename)).get()), deref(v1), (<bool>do_sort), (<bool>force))
        cdef libcpp_vector[_SwathMap].iterator it_swath_maps = v1.begin()
        replace_0 = []
        while it_swath_maps != v1.end():
            item1 = SwathMap.__new__(SwathMap)
            item1.inst = shared_ptr[_SwathMap](new _SwathMap(deref(it_swath_maps)))
            replace_0.append(item1)
            inc(it_swath_maps)
        swath_maps[:] = replace_0
        del v1
    
    def readSwathWindows(self,  filename , list swath_prec_lower , list swath_prec_upper ):
        """
        readSwathWindows(self, filename: Union[bytes, str, String] , swath_prec_lower: List[float] , swath_prec_upper: List[float] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(swath_prec_lower, list) and all(isinstance(elemt_rec, float) for elemt_rec in swath_prec_lower), 'arg swath_prec_lower wrong type'
        assert isinstance(swath_prec_upper, list) and all(isinstance(elemt_rec, float) for elemt_rec in swath_prec_upper), 'arg swath_prec_upper wrong type'
    
        cdef libcpp_vector[double] v1 = swath_prec_lower
        cdef libcpp_vector[double] v2 = swath_prec_upper
        self.inst.get().readSwathWindows(deref((convString(filename)).get()), v1, v2)
        swath_prec_upper[:] = v2
        swath_prec_lower[:] = v1 

cdef class TransitionPQPFile:
    """
    Cython implementation of _TransitionPQPFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransitionPQPFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TransitionPQPFile rv = TransitionPQPFile.__new__(TransitionPQPFile)
       rv.inst = shared_ptr[_TransitionPQPFile](new _TransitionPQPFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TransitionPQPFile rv = TransitionPQPFile.__new__(TransitionPQPFile)
       rv.inst = shared_ptr[_TransitionPQPFile](new _TransitionPQPFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TransitionPQPFile](new _TransitionPQPFile())
    
    def _init_1(self, TransitionPQPFile in_0 ):
        """
        _init_1(self, in_0: TransitionPQPFile ) -> None
        """
        assert isinstance(in_0, TransitionPQPFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TransitionPQPFile](new _TransitionPQPFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TransitionPQPFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TransitionPQPFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def convertTargetedExperimentToPQP(self, bytes filename , TargetedExperiment targeted_exp ):
        """
        convertTargetedExperimentToPQP(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None
        Write out a targeted experiment (TraML structure) into a PQP file
        
        :param filename: The output file
        :param targeted_exp: The targeted experiment
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
    
        self.inst.get().convertTargetedExperimentToPQP((<char *>filename), (deref(targeted_exp.inst.get())))
    
    def _convertPQPToTargetedExperiment_0(self, bytes filename , TargetedExperiment targeted_exp , bool legacy_traml_id ):
        """
        _convertPQPToTargetedExperiment_0(self, filename: bytes , targeted_exp: TargetedExperiment , legacy_traml_id: bool ) -> None
        Read in a PQP file and construct a targeted experiment (TraML structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
        assert isinstance(legacy_traml_id, pybool_t), 'arg legacy_traml_id wrong type'
    
    
    
        self.inst.get().convertPQPToTargetedExperiment((<char *>filename), (deref(targeted_exp.inst.get())), (<bool>legacy_traml_id))
    
    def _convertPQPToTargetedExperiment_1(self, bytes filename , LightTargetedExperiment targeted_exp , bool legacy_traml_id ):
        """
        _convertPQPToTargetedExperiment_1(self, filename: bytes , targeted_exp: LightTargetedExperiment , legacy_traml_id: bool ) -> None
        Read in a PQP file and construct a targeted experiment (Light transition structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(targeted_exp, LightTargetedExperiment), 'arg targeted_exp wrong type'
        assert isinstance(legacy_traml_id, pybool_t), 'arg legacy_traml_id wrong type'
    
    
    
        self.inst.get().convertPQPToTargetedExperiment((<char *>filename), (deref(targeted_exp.inst.get())), (<bool>legacy_traml_id))
    
    def convertPQPToTargetedExperiment(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: convertPQPToTargetedExperiment(self, filename: bytes , targeted_exp: TargetedExperiment , legacy_traml_id: bool ) -> None
          :noindex:
        
        Read in a PQP file and construct a targeted experiment (TraML structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
        
        .. rubric:: Overload:
        .. py:function:: convertPQPToTargetedExperiment(self, filename: bytes , targeted_exp: LightTargetedExperiment , legacy_traml_id: bool ) -> None
          :noindex:
        
        Read in a PQP file and construct a targeted experiment (Light transition structure)
        
        :param filename: The input file
        :param targeted_exp: The output targeted experiment
        :param legacy_traml_id: Should legacy TraML IDs be used (boolean)?
    
        """
        if (len(args)==3) and (isinstance(args[0], bytes)) and (isinstance(args[1], TargetedExperiment)) and (isinstance(args[2], pybool_t)):
            return self._convertPQPToTargetedExperiment_0(*args)
        elif (len(args)==3) and (isinstance(args[0], bytes)) and (isinstance(args[1], LightTargetedExperiment)) and (isinstance(args[2], pybool_t)):
            return self._convertPQPToTargetedExperiment_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def convertTargetedExperimentToTSV(self, bytes filename , TargetedExperiment targeted_exp ):
        """
        convertTargetedExperimentToTSV(self, filename: bytes , targeted_exp: TargetedExperiment ) -> None
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
    
        self.inst.get().convertTargetedExperimentToTSV((<char *>filename), (deref(targeted_exp.inst.get())))
    
    def _convertTSVToTargetedExperiment_0(self, bytes filename , int filetype , TargetedExperiment targeted_exp ):
        """
        _convertTSVToTargetedExperiment_0(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert filetype in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg filetype wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
    
    
        self.inst.get().convertTSVToTargetedExperiment((<char *>filename), (<_FileType>filetype), (deref(targeted_exp.inst.get())))
    
    def _convertTSVToTargetedExperiment_1(self, bytes filename , int filetype , LightTargetedExperiment targeted_exp ):
        """
        _convertTSVToTargetedExperiment_1(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None
        """
        assert isinstance(filename, bytes), 'arg filename wrong type'
        assert filetype in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], 'arg filetype wrong type'
        assert isinstance(targeted_exp, LightTargetedExperiment), 'arg targeted_exp wrong type'
    
    
    
        self.inst.get().convertTSVToTargetedExperiment((<char *>filename), (<_FileType>filetype), (deref(targeted_exp.inst.get())))
    
    def convertTSVToTargetedExperiment(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: TargetedExperiment ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: convertTSVToTargetedExperiment(self, filename: bytes , filetype: int , targeted_exp: LightTargetedExperiment ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], bytes)) and (args[1] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]) and (isinstance(args[2], TargetedExperiment)):
            return self._convertTSVToTargetedExperiment_0(*args)
        elif (len(args)==3) and (isinstance(args[0], bytes)) and (args[1] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46]) and (isinstance(args[2], LightTargetedExperiment)):
            return self._convertTSVToTargetedExperiment_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def validateTargetedExperiment(self, TargetedExperiment targeted_exp ):
        """
        validateTargetedExperiment(self, targeted_exp: TargetedExperiment ) -> None
        """
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
    
        self.inst.get().validateTargetedExperiment((deref(targeted_exp.inst.get()))) 
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
