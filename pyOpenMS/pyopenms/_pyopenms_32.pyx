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

cdef class __RETURN_STATUS:
    None
    SOLVED = 0
    ITERATION_EXCEEDED = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AQS_featureConcentration:
    """
    Cython implementation of _AQS_featureConcentration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AQS_featureConcentration.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property feature:
        def __set__(self, Feature feature):
        
            self.inst.get().feature = (deref(feature.inst.get()))
        
    
        def __get__(self):
            cdef _Feature * _r = new _Feature(self.inst.get().feature)
            cdef Feature py_result = Feature.__new__(Feature)
            py_result.inst = shared_ptr[_Feature](_r)
            return py_result
    
    property IS_feature:
        def __set__(self, Feature IS_feature):
        
            self.inst.get().IS_feature = (deref(IS_feature.inst.get()))
        
    
        def __get__(self):
            cdef _Feature * _r = new _Feature(self.inst.get().IS_feature)
            cdef Feature py_result = Feature.__new__(Feature)
            py_result.inst = shared_ptr[_Feature](_r)
            return py_result
    
    property actual_concentration:
        def __set__(self, double actual_concentration):
        
            self.inst.get().actual_concentration = (<double>actual_concentration)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().actual_concentration
            py_result = <double>_r
            return py_result
    
    property IS_actual_concentration:
        def __set__(self, double IS_actual_concentration):
        
            self.inst.get().IS_actual_concentration = (<double>IS_actual_concentration)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().IS_actual_concentration
            py_result = <double>_r
            return py_result
    
    property concentration_units:
        def __set__(self,  concentration_units):
        
            self.inst.get().concentration_units = deref((convString(concentration_units)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().concentration_units
            py_result = convOutputString(_r)
            return py_result
    
    property dilution_factor:
        def __set__(self, double dilution_factor):
        
            self.inst.get().dilution_factor = (<double>dilution_factor)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().dilution_factor
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef AQS_featureConcentration rv = AQS_featureConcentration.__new__(AQS_featureConcentration)
       rv.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AQS_featureConcentration rv = AQS_featureConcentration.__new__(AQS_featureConcentration)
       rv.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration())
    
    def _init_1(self, AQS_featureConcentration in_0 ):
        """
        _init_1(self, in_0: AQS_featureConcentration ) -> None
        """
        assert isinstance(in_0, AQS_featureConcentration), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AQS_featureConcentration ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AQS_featureConcentration)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class AQS_runConcentration:
    """
    Cython implementation of _AQS_runConcentration

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AQS_runConcentration.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property sample_name:
        def __set__(self,  sample_name):
        
            self.inst.get().sample_name = deref((convString(sample_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().sample_name
            py_result = convOutputString(_r)
            return py_result
    
    property component_name:
        def __set__(self,  component_name):
        
            self.inst.get().component_name = deref((convString(component_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().component_name
            py_result = convOutputString(_r)
            return py_result
    
    property IS_component_name:
        def __set__(self,  IS_component_name):
        
            self.inst.get().IS_component_name = deref((convString(IS_component_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().IS_component_name
            py_result = convOutputString(_r)
            return py_result
    
    property actual_concentration:
        def __set__(self, double actual_concentration):
        
            self.inst.get().actual_concentration = (<double>actual_concentration)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().actual_concentration
            py_result = <double>_r
            return py_result
    
    property IS_actual_concentration:
        def __set__(self, double IS_actual_concentration):
        
            self.inst.get().IS_actual_concentration = (<double>IS_actual_concentration)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().IS_actual_concentration
            py_result = <double>_r
            return py_result
    
    property concentration_units:
        def __set__(self,  concentration_units):
        
            self.inst.get().concentration_units = deref((convString(concentration_units)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().concentration_units
            py_result = convOutputString(_r)
            return py_result
    
    property dilution_factor:
        def __set__(self, double dilution_factor):
        
            self.inst.get().dilution_factor = (<double>dilution_factor)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().dilution_factor
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef AQS_runConcentration rv = AQS_runConcentration.__new__(AQS_runConcentration)
       rv.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AQS_runConcentration rv = AQS_runConcentration.__new__(AQS_runConcentration)
       rv.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration())
    
    def _init_1(self, AQS_runConcentration in_0 ):
        """
        _init_1(self, in_0: AQS_runConcentration ) -> None
        """
        assert isinstance(in_0, AQS_runConcentration), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AQS_runConcentration ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AQS_runConcentration)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class AbsoluteQuantitationStandards:
    """
    Cython implementation of _AbsoluteQuantitationStandards

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitationStandards.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AbsoluteQuantitationStandards rv = AbsoluteQuantitationStandards.__new__(AbsoluteQuantitationStandards)
       rv.inst = shared_ptr[_AbsoluteQuantitationStandards](new _AbsoluteQuantitationStandards(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AbsoluteQuantitationStandards rv = AbsoluteQuantitationStandards.__new__(AbsoluteQuantitationStandards)
       rv.inst = shared_ptr[_AbsoluteQuantitationStandards](new _AbsoluteQuantitationStandards(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AbsoluteQuantitationStandards](new _AbsoluteQuantitationStandards())
    
    def _init_1(self, AbsoluteQuantitationStandards in_0 ):
        """
        _init_1(self, in_0: AbsoluteQuantitationStandards ) -> None
        """
        assert isinstance(in_0, AbsoluteQuantitationStandards), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AbsoluteQuantitationStandards](new _AbsoluteQuantitationStandards((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AbsoluteQuantitationStandards ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AbsoluteQuantitationStandards)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getComponentFeatureConcentrations(self, list run_concentrations , list feature_maps ,  component_name , list feature_concentrations ):
        """
        getComponentFeatureConcentrations(self, run_concentrations: List[AQS_runConcentration] , feature_maps: List[FeatureMap] , component_name: Union[bytes, str, String] , feature_concentrations: List[AQS_featureConcentration] ) -> None
        """
        assert isinstance(run_concentrations, list) and all(isinstance(elemt_rec, AQS_runConcentration) for elemt_rec in run_concentrations), 'arg run_concentrations wrong type'
        assert isinstance(feature_maps, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in feature_maps), 'arg feature_maps wrong type'
        assert (isinstance(component_name, str) or isinstance(component_name, bytes) or isinstance(component_name, String)), 'arg component_name wrong type'
        assert isinstance(feature_concentrations, list) and all(isinstance(elemt_rec, AQS_featureConcentration) for elemt_rec in feature_concentrations), 'arg feature_concentrations wrong type'
        cdef libcpp_vector[_AQS_runConcentration] * v0 = new libcpp_vector[_AQS_runConcentration]()
        cdef AQS_runConcentration item0
        for item0 in run_concentrations:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_FeatureMap] * v1 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item1
        for item1 in feature_maps:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_AQS_featureConcentration] * v3 = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration item3
        for item3 in feature_concentrations:
            v3.push_back(deref(item3.inst.get()))
        self.inst.get().getComponentFeatureConcentrations(deref(v0), deref(v1), deref((convString(component_name)).get()), deref(v3))
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_feature_concentrations = v3.begin()
        replace_0 = []
        while it_feature_concentrations != v3.end():
            item3 = AQS_featureConcentration.__new__(AQS_featureConcentration)
            item3.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_feature_concentrations)))
            replace_0.append(item3)
            inc(it_feature_concentrations)
        feature_concentrations[:] = replace_0
        del v3
        cdef libcpp_vector[_FeatureMap].iterator it_feature_maps = v1.begin()
        replace_0 = []
        while it_feature_maps != v1.end():
            item1 = FeatureMap.__new__(FeatureMap)
            item1.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_feature_maps)))
            replace_0.append(item1)
            inc(it_feature_maps)
        feature_maps[:] = replace_0
        del v1
        cdef libcpp_vector[_AQS_runConcentration].iterator it_run_concentrations = v0.begin()
        replace_0 = []
        while it_run_concentrations != v0.end():
            item0 = AQS_runConcentration.__new__(AQS_runConcentration)
            item0.inst = shared_ptr[_AQS_runConcentration](new _AQS_runConcentration(deref(it_run_concentrations)))
            replace_0.append(item0)
            inc(it_run_concentrations)
        run_concentrations[:] = replace_0
        del v0 

cdef class ChromatogramRangeManager:
    """
    Cython implementation of _ChromatogramRangeManager

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramRangeManager.html>`_

    Range manager for chromatograms
    
    This class manages retention time, m/z, and intensity ranges for multiple chromatograms.
    It extends the basic RangeManager to provide specialized functionality for chromatogram data.
    
    The template parameters for the base RangeManager are ordered differently than in SpectrumRangeManager:
    - RangeRT (retention time) is the first parameter, as it's the primary dimension for chromatograms
    - RangeIntensity is the second parameter
    - RangeMZ is the third parameter
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramRangeManager rv = ChromatogramRangeManager.__new__(ChromatogramRangeManager)
       rv.inst = shared_ptr[_ChromatogramRangeManager](new _ChromatogramRangeManager(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramRangeManager rv = ChromatogramRangeManager.__new__(ChromatogramRangeManager)
       rv.inst = shared_ptr[_ChromatogramRangeManager](new _ChromatogramRangeManager(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChromatogramRangeManager](new _ChromatogramRangeManager())
    
    def _init_1(self, ChromatogramRangeManager in_0 ):
        """
        _init_1(self, in_0: ChromatogramRangeManager ) -> None
        """
        assert isinstance(in_0, ChromatogramRangeManager), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramRangeManager](new _ChromatogramRangeManager((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramRangeManager ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramRangeManager)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def clearRanges(self):
        """
        clearRanges(self) -> None
        """
        self.inst.get().clearRanges()
    
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
    
    def getMinIntensity(self):
        """
        getMinIntensity(self) -> float
        """
        cdef double _r = self.inst.get().getMinIntensity()
        py_result = <double>_r
        return py_result
    
    def getMaxIntensity(self):
        """
        getMaxIntensity(self) -> float
        """
        cdef double _r = self.inst.get().getMaxIntensity()
        py_result = <double>_r
        return py_result 

cdef class ColumnHeader:
    """
    Cython implementation of _ColumnHeader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMap_1_1ColumnHeader.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property filename:
        def __set__(self,  filename):
        
            self.inst.get().filename = deref((convString(filename)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().filename
            py_result = convOutputString(_r)
            return py_result
    
    property label:
        def __set__(self,  label):
        
            self.inst.get().label = deref((convString(label)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().label
            py_result = convOutputString(_r)
            return py_result
    
    property size:
        def __set__(self,  size):
        
            self.inst.get().size = (<size_t>size)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().size
            py_result = <size_t>_r
            return py_result
    
    property unique_id:
        def __set__(self,  unique_id):
        
            self.inst.get().unique_id = (<uint64_t>unique_id)
        
    
        def __get__(self):
            cdef uint64_t _r = self.inst.get().unique_id
            py_result = <uint64_t>_r
            return py_result
    
    def __copy__(self):
       cdef ColumnHeader rv = ColumnHeader.__new__(ColumnHeader)
       rv.inst = shared_ptr[_ColumnHeader](new _ColumnHeader(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ColumnHeader rv = ColumnHeader.__new__(ColumnHeader)
       rv.inst = shared_ptr[_ColumnHeader](new _ColumnHeader(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ColumnHeader](new _ColumnHeader())
    
    def _init_1(self, ColumnHeader in_0 ):
        """
        _init_1(self, in_0: ColumnHeader ) -> None
        """
        assert isinstance(in_0, ColumnHeader), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ColumnHeader](new _ColumnHeader((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ColumnHeader ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ColumnHeader)):
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
        if not isinstance(other, ColumnHeader):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ColumnHeader other_casted = other
        cdef ColumnHeader self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ConsensusMap:
    """
    Cython implementation of _ConsensusMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMap_1_1ConsensusMap.html>`_
      -- Inherits from ['UniqueIdInterface', 'DocumentIdentifier', 'RangeManagerRtMzInt', 'MetaInfoInterface']

    A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements
    (ConsensusFeature) which in turn represent analytes that have been
    quantified across multiple LC-MS/MS experiments. Each analyte in a
    ConsensusFeature is linked to its original LC-MS/MS run, the links are
    maintained by the ConsensusMap class.
    The map is implemented as a vector of elements of type ConsensusFeature.
    
    To be consistent, all maps who are referenced by ConsensusFeature objects
    (through a unique id) need to be registered in this class.
    
    This class supports direct iteration in Python.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ConsensusMap rv = ConsensusMap.__new__(ConsensusMap)
       rv.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ConsensusMap rv = ConsensusMap.__new__(ConsensusMap)
       rv.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ConsensusMap](new _ConsensusMap())
    
    def _init_1(self, ConsensusMap in_0 ):
        """
        _init_1(self, in_0: ConsensusMap ) -> None
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ConsensusMap](new _ConsensusMap((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ConsensusMap ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusMap)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def reserve(self,  s ):
        """
        reserve(self, s: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
    
        self.inst.get().reserve((<size_t>s))
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> ConsensusFeature
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _ConsensusFeature * _r = new _ConsensusFeature(deref(self.inst.get())[(<size_t>in_0)])
        cdef ConsensusFeature py_result = ConsensusFeature.__new__(ConsensusFeature)
        py_result.inst = shared_ptr[_ConsensusFeature](_r)
        return py_result
    def __setitem__(self, key, ConsensusFeature value):
        """Cython signature: ConsensusFeature & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def push_back(self, ConsensusFeature spec ):
        """
        push_back(self, spec: ConsensusFeature ) -> None
        """
        assert isinstance(spec, ConsensusFeature), 'arg spec wrong type'
    
        self.inst.get().push_back((deref(spec.inst.get())))
    
    def appendRows(self, ConsensusMap in_0 ):
        """
        appendRows(self, in_0: ConsensusMap ) -> ConsensusMap
        Add consensus map entries as new rows
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
    
        cdef _ConsensusMap * _r = new _ConsensusMap(self.inst.get().appendRows((deref(in_0.inst.get()))))
        cdef ConsensusMap py_result = ConsensusMap.__new__(ConsensusMap)
        py_result.inst = shared_ptr[_ConsensusMap](_r)
        return py_result
    
    def appendColumns(self, ConsensusMap in_0 ):
        """
        appendColumns(self, in_0: ConsensusMap ) -> ConsensusMap
        Add consensus map entries as new columns
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
    
        cdef _ConsensusMap * _r = new _ConsensusMap(self.inst.get().appendColumns((deref(in_0.inst.get()))))
        cdef ConsensusMap py_result = ConsensusMap.__new__(ConsensusMap)
        py_result.inst = shared_ptr[_ConsensusMap](_r)
        return py_result
    
    def _clear_0(self, bool clear_meta_data ):
        """
        _clear_0(self, clear_meta_data: bool ) -> None
        Clears all data and meta data
        """
        assert isinstance(clear_meta_data, pybool_t), 'arg clear_meta_data wrong type'
    
        self.inst.get().clear((<bool>clear_meta_data))
    
    def _clear_1(self):
        """
        _clear_1(self) -> None
        """
        self.inst.get().clear()
    
    def clear(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: clear(self, clear_meta_data: bool ) -> None
          :noindex:
        
        Clears all data and meta data

        
        .. rubric:: Overload:
        .. py:function:: clear(self, ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._clear_0(*args)
        elif not args:
            return self._clear_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def updateRanges(self):
        """
        updateRanges(self) -> None
        """
        self.inst.get().updateRanges()
    
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
    
    def setProteinIdentifications(self, list in_0 ):
        """
        setProteinIdentifications(self, in_0: List[ProteinIdentification] ) -> None
        Sets the protein identifications
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setProteinIdentifications(deref(v0))
        del v0
    
    def getUnassignedPeptideIdentifications(self):
        """
        getUnassignedPeptideIdentifications(self) -> PeptideIdentificationList
        """
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().getUnassignedPeptideIdentifications())
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result
    
    def setUnassignedPeptideIdentifications(self, PeptideIdentificationList unassigned_peptide_identifications ):
        """
        setUnassignedPeptideIdentifications(self, unassigned_peptide_identifications: PeptideIdentificationList ) -> None
        Sets the unassigned PeptideIdentificationList
        """
        assert isinstance(unassigned_peptide_identifications, PeptideIdentificationList), 'arg unassigned_peptide_identifications wrong type'
    
        self.inst.get().setUnassignedPeptideIdentifications((deref(unassigned_peptide_identifications.inst.get())))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        Returns a const reference to the description of the applied data processing
        """
        _r = self.inst.get().getDataProcessing()
        py_result = []
        cdef libcpp_vector[_DataProcessing].iterator it__r = _r.begin()
        cdef DataProcessing item_py_result
        while it__r != _r.end():
           item_py_result = DataProcessing.__new__(DataProcessing)
           item_py_result.inst = shared_ptr[_DataProcessing](new _DataProcessing(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setDataProcessing(self, list in_0 ):
        """
        setDataProcessing(self, in_0: List[DataProcessing] ) -> None
        Sets the description of the applied data processing
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, DataProcessing) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_DataProcessing] * v0 = new libcpp_vector[_DataProcessing]()
        cdef DataProcessing item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setDataProcessing(deref(v0))
        del v0
    
    def _setPrimaryMSRunPath_0(self, list s ):
        """
        _setPrimaryMSRunPath_0(self, s: List[bytes] ) -> None
        Sets the file paths to the primary MS run (stored in ColumnHeaders)
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
        self.inst.get().setPrimaryMSRunPath(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def _setPrimaryMSRunPath_1(self, list s , MSExperiment e ):
        """
        _setPrimaryMSRunPath_1(self, s: List[bytes] , e: MSExperiment ) -> None
        """
        assert isinstance(s, list) and all(isinstance(li, bytes) for li in s), 'arg s wrong type'
        assert isinstance(e, MSExperiment), 'arg e wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in s:
           v0.push_back(_String(<char *>item0))
    
        self.inst.get().setPrimaryMSRunPath(deref(v0), (deref(e.inst.get())))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        s[:] = replace
        del v0
    
    def setPrimaryMSRunPath(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setPrimaryMSRunPath(self, s: List[bytes] ) -> None
          :noindex:
        
        Sets the file paths to the primary MS run (stored in ColumnHeaders)

        
        .. rubric:: Overload:
        .. py:function:: setPrimaryMSRunPath(self, s: List[bytes] , e: MSExperiment ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])):
            return self._setPrimaryMSRunPath_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], MSExperiment)):
            return self._setPrimaryMSRunPath_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPrimaryMSRunPath(self, list toFill ):
        """
        getPrimaryMSRunPath(self, toFill: List[bytes] ) -> None
        Returns the MS run path (stored in ColumnHeaders)
        """
        assert isinstance(toFill, list) and all(isinstance(li, bytes) for li in toFill), 'arg toFill wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in toFill:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getPrimaryMSRunPath(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        toFill[:] = replace
        del v0
    
    def _sortByIntensity_0(self, bool reverse ):
        """
        _sortByIntensity_0(self, reverse: bool ) -> None
        Sorts the peaks according to ascending intensity.
        """
        assert isinstance(reverse, pybool_t), 'arg reverse wrong type'
    
        self.inst.get().sortByIntensity((<bool>reverse))
    
    def _sortByIntensity_1(self):
        """
        _sortByIntensity_1(self) -> None
        """
        self.inst.get().sortByIntensity()
    
    def sortByIntensity(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: sortByIntensity(self, reverse: bool ) -> None
          :noindex:
        
        Sorts the peaks according to ascending intensity.

        
        .. rubric:: Overload:
        .. py:function:: sortByIntensity(self, ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._sortByIntensity_0(*args)
        elif not args:
            return self._sortByIntensity_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def sortByRT(self):
        """
        sortByRT(self) -> None
        Sorts the peaks according to RT position
        """
        self.inst.get().sortByRT()
    
    def sortByMZ(self):
        """
        sortByMZ(self) -> None
        Sorts the peaks according to m/z position
        """
        self.inst.get().sortByMZ()
    
    def sortByPosition(self):
        """
        sortByPosition(self) -> None
        Lexicographically sorts the peaks by their position (First RT then m/z)
        """
        self.inst.get().sortByPosition()
    
    def _sortByQuality_0(self, bool reverse ):
        """
        _sortByQuality_0(self, reverse: bool ) -> None
        Sorts the peaks according to ascending quality.
        """
        assert isinstance(reverse, pybool_t), 'arg reverse wrong type'
    
        self.inst.get().sortByQuality((<bool>reverse))
    
    def _sortByQuality_1(self):
        """
        _sortByQuality_1(self) -> None
        """
        self.inst.get().sortByQuality()
    
    def sortByQuality(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: sortByQuality(self, reverse: bool ) -> None
          :noindex:
        
        Sorts the peaks according to ascending quality.

        
        .. rubric:: Overload:
        .. py:function:: sortByQuality(self, ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._sortByQuality_0(*args)
        elif not args:
            return self._sortByQuality_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def sortBySize(self):
        """
        sortBySize(self) -> None
        Sorts with respect to the size (number of elements)
        """
        self.inst.get().sortBySize()
    
    def sortByMaps(self):
        """
        sortByMaps(self) -> None
        Sorts with respect to the sets of maps covered by the consensus features (lexicographically)
        """
        self.inst.get().sortByMaps()
    
    def getExperimentType(self):
        """
        getExperimentType(self) -> Union[bytes, str, String]
        Non-mutable access to the experiment type
        """
        cdef _String _r = self.inst.get().getExperimentType()
        py_result = convOutputString(_r)
        return py_result
    
    def setExperimentType(self,  experiment_type ):
        """
        setExperimentType(self, experiment_type: Union[bytes, str, String] ) -> None
        Mutable access to the experiment type
        """
        assert (isinstance(experiment_type, str) or isinstance(experiment_type, bytes) or isinstance(experiment_type, String)), 'arg experiment_type wrong type'
    
        self.inst.get().setExperimentType(deref((convString(experiment_type)).get()))
    
    def sortPeptideIdentificationsByMapIndex(self):
        """
        sortPeptideIdentificationsByMapIndex(self) -> None
        Sorts PeptideIdentifications of consensus features with respect to their map index.
        """
        self.inst.get().sortPeptideIdentificationsByMapIndex()
    
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
    
    def setIdentifier(self,  id ):
        """
        setIdentifier(self, id: Union[bytes, str, String] ) -> None
        Sets document identifier (e.g. an LSID)
        """
        assert (isinstance(id, str) or isinstance(id, bytes) or isinstance(id, String)), 'arg id wrong type'
    
        self.inst.get().setIdentifier(deref((convString(id)).get()))
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        Retrieve document identifier (e.g. an LSID)
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setLoadedFileType(self,  file_name ):
        """
        setLoadedFileType(self, file_name: Union[bytes, str, String] ) -> None
        Sets the file_type according to the type of the file loaded from, preferably done whilst loading
        """
        assert (isinstance(file_name, str) or isinstance(file_name, bytes) or isinstance(file_name, String)), 'arg file_name wrong type'
    
        self.inst.get().setLoadedFileType(deref((convString(file_name)).get()))
    
    def getLoadedFileType(self):
        """
        getLoadedFileType(self) -> int
        Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded
        """
        cdef int _r = self.inst.get().getLoadedFileType()
        py_result = <int>_r
        return py_result
    
    def setLoadedFilePath(self,  file_name ):
        """
        setLoadedFilePath(self, file_name: Union[bytes, str, String] ) -> None
        Sets the file_name according to absolute path of the file loaded, preferably done whilst loading
        """
        assert (isinstance(file_name, str) or isinstance(file_name, bytes) or isinstance(file_name, String)), 'arg file_name wrong type'
    
        self.inst.get().setLoadedFilePath(deref((convString(file_name)).get()))
    
    def getLoadedFilePath(self):
        """
        getLoadedFilePath(self) -> Union[bytes, str, String]
        Returns the file_name which is the absolute path to the file loaded
        """
        cdef _String _r = self.inst.get().getLoadedFilePath()
        py_result = convOutputString(_r)
        return py_result
    
    def getMinRT(self):
        """
        getMinRT(self) -> float
        Returns the minimum RT
        """
        cdef double _r = self.inst.get().getMinRT()
        py_result = <double>_r
        return py_result
    
    def getMaxRT(self):
        """
        getMaxRT(self) -> float
        Returns the maximum RT
        """
        cdef double _r = self.inst.get().getMaxRT()
        py_result = <double>_r
        return py_result
    
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
        if not isinstance(other, ConsensusMap):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ConsensusMap other_casted = other
        cdef ConsensusMap self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef ConsensusFeature out
        while it != self.inst.get().end():
            out = ConsensusFeature.__new__(ConsensusFeature)
            out.inst = shared_ptr[_ConsensusFeature](new _ConsensusFeature(deref(it)))
            yield out
            inc(it)
    
    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId))

    def getColumnHeaders(self):

        # ColumnHeaders is a type alias for Map<..> which can not be
        # handled by autowrap automatically. So we have to provide a manual
        # converter here:
        #
        # (the wrapper works on linux, but msvc complains
        # a lot about the generated code... uwe schmitt)

        cdef ColumnHeaders _r = self.inst.get().getColumnHeaders()
        py_result = dict()
        cdef ColumnHeaders_iterator it__r = _r.begin()
        cdef ColumnHeader item_py_result
        while it__r != _r.end():
           item_py_result = ColumnHeader.__new__(ColumnHeader)
           item_py_result.inst = shared_ptr[_ColumnHeader](new _ColumnHeader((deref(it__r)).second))
           py_result[<UInt64>(deref(it__r).first)] = item_py_result
           inc(it__r)
        return py_result

    def setColumnHeaders(self, dict in_0 ):
        assert isinstance(in_0, dict) and all(isinstance(k, int) for k in in_0.keys()) and all(isinstance(v, ColumnHeader) for v in in_0.values()), 'arg in_0 wrong type'
        cdef ColumnHeaders v0
        for key, value in in_0.items():
           v0[<UInt64> key] = deref((<ColumnHeader>value).inst.get())

        # we have to utilize AutowrapRefHolder here, because cython does not
        # like
        #
        #     self.inst.get().getColumnHeaders() = v0
        #
        # and if you choose
        #
        #     cdef ColumnHeaders & ref = self.inst.get().getColumnHeaders()
        #     ref = v0
        #
        # the c++ compiler will complain, as cython generates invalid code in
        # this case: Cython creates two statements from the first line:
        #
        #     ColumnHeader & ref;   // invalid: ref not initailized
        #     ref = self.inst.get().getColumnHeaders()
        #
        # wrapping this ref into a class holding the ref and using pointers
        # works. the following code results in c++ code simialar to
        #
        #    AutowrapRefHolder<ColumnHeader> * refholder;
        #    refholder = new AutowrapRefHolder<ColumnHeader>(self.inst...)
        #    refholder.assign(v0)

        cdef AutowrapRefHolder[ColumnHeaders] * refholder
        refholder = new AutowrapRefHolder[ColumnHeaders](self.inst.get().getColumnHeaders())
        refholder.assign(v0)

        cdef replace_in_0 = dict()
        cdef ColumnHeaders_iterator it_in_0 = v0.begin()
        cdef ColumnHeader item_in_0
        while it_in_0 != v0.end():
           item_in_0 = ColumnHeader.__new__(ColumnHeader)
           item_in_0.inst = shared_ptr[_ColumnHeader](new _ColumnHeader((deref(it_in_0)).second))
           replace_in_0[<UInt64> deref(it_in_0).first] = item_in_0
           inc(it_in_0)
        in_0.clear()
        in_0.update(replace_in_0) 

cdef class DigestionEnzymeRNA:
    """
    Cython implementation of _DigestionEnzymeRNA

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionEnzymeRNA.html>`_
      -- Inherits from ['DigestionEnzyme']

    Representation of a digestion enzyme for RNA (RNase)
    
    The cutting sites of these enzymes are defined using two different mechanisms:
    First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
    This is the same mechanism that is used for proteases (see ProteaseDigestion).
    However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
    Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
    If both expressions match, then there is a cutting site between the two ribonucleotides.
    
    There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
    A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DigestionEnzymeRNA rv = DigestionEnzymeRNA.__new__(DigestionEnzymeRNA)
       rv.inst = shared_ptr[_DigestionEnzymeRNA](new _DigestionEnzymeRNA(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DigestionEnzymeRNA rv = DigestionEnzymeRNA.__new__(DigestionEnzymeRNA)
       rv.inst = shared_ptr[_DigestionEnzymeRNA](new _DigestionEnzymeRNA(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DigestionEnzymeRNA](new _DigestionEnzymeRNA())
    
    def _init_1(self, DigestionEnzymeRNA in_0 ):
        """
        _init_1(self, in_0: DigestionEnzymeRNA ) -> None
        """
        assert isinstance(in_0, DigestionEnzymeRNA), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DigestionEnzymeRNA](new _DigestionEnzymeRNA((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DigestionEnzymeRNA ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DigestionEnzymeRNA)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setCutsAfterRegEx(self,  value ):
        """
        setCutsAfterRegEx(self, value: Union[bytes, str, String] ) -> None
        Sets the "cuts after ..." regular expression
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setCutsAfterRegEx(deref((convString(value)).get()))
    
    def getCutsAfterRegEx(self):
        """
        getCutsAfterRegEx(self) -> Union[bytes, str, String]
        Returns the "cuts after ..." regular expression
        """
        cdef _String _r = self.inst.get().getCutsAfterRegEx()
        py_result = convOutputString(_r)
        return py_result
    
    def setCutsBeforeRegEx(self,  value ):
        """
        setCutsBeforeRegEx(self, value: Union[bytes, str, String] ) -> None
        Sets the "cuts before ..." regular expression
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setCutsBeforeRegEx(deref((convString(value)).get()))
    
    def getCutsBeforeRegEx(self):
        """
        getCutsBeforeRegEx(self) -> Union[bytes, str, String]
        Returns the "cuts before ..." regular expression
        """
        cdef _String _r = self.inst.get().getCutsBeforeRegEx()
        py_result = convOutputString(_r)
        return py_result
    
    def setThreePrimeGain(self,  value ):
        """
        setThreePrimeGain(self, value: Union[bytes, str, String] ) -> None
        Sets the 3' gain
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setThreePrimeGain(deref((convString(value)).get()))
    
    def setFivePrimeGain(self,  value ):
        """
        setFivePrimeGain(self, value: Union[bytes, str, String] ) -> None
        Sets the 5' gain
        """
        assert (isinstance(value, str) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
    
        self.inst.get().setFivePrimeGain(deref((convString(value)).get()))
    
    def getThreePrimeGain(self):
        """
        getThreePrimeGain(self) -> Union[bytes, str, String]
        Returns the 3' gain
        """
        cdef _String _r = self.inst.get().getThreePrimeGain()
        py_result = convOutputString(_r)
        return py_result
    
    def getFivePrimeGain(self):
        """
        getFivePrimeGain(self) -> Union[bytes, str, String]
        Returns the 5' gain
        """
        cdef _String _r = self.inst.get().getFivePrimeGain()
        py_result = convOutputString(_r)
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
        if not isinstance(other, DigestionEnzymeRNA):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DigestionEnzymeRNA other_casted = other
        cdef DigestionEnzymeRNA self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class FeatureGroupingAlgorithmUnlabeled:
    """
    Cython implementation of _FeatureGroupingAlgorithmUnlabeled

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmUnlabeled.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureGroupingAlgorithmUnlabeled](new _FeatureGroupingAlgorithmUnlabeled())
    
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
    
    def addToGroup(self,  map_id , FeatureMap feature_map ):
        """
        addToGroup(self, map_id: int , feature_map: FeatureMap ) -> None
        """
        assert isinstance(map_id, int), 'arg map_id wrong type'
        assert isinstance(feature_map, FeatureMap), 'arg feature_map wrong type'
    
    
        self.inst.get().addToGroup((<int>map_id), (deref(feature_map.inst.get())))
    
    def setReference(self,  map_id , FeatureMap map ):
        """
        setReference(self, map_id: int , map: FeatureMap ) -> None
        """
        assert isinstance(map_id, int), 'arg map_id wrong type'
        assert isinstance(map, FeatureMap), 'arg map wrong type'
    
    
        self.inst.get().setReference((<int>map_id), (deref(map.inst.get())))
    
    def getResultMap(self):
        """
        getResultMap(self) -> ConsensusMap
        """
        cdef _ConsensusMap * _r = new _ConsensusMap(self.inst.get().getResultMap())
        cdef ConsensusMap py_result = ConsensusMap.__new__(ConsensusMap)
        py_result.inst = shared_ptr[_ConsensusMap](_r)
        return py_result
    
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

cdef class IdXMLFile:
    """
    Cython implementation of _IdXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IdXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        Used to load and store idXML files
        """
        self.inst = shared_ptr[_IdXMLFile](new _IdXMLFile())
    
    def _load_0(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ):
        """
        _load_0(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
        Loads the identifications of an idXML file without identifier
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _load_1(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ):
        """
        _load_1(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
        Loads the identifications of an idXML file without identifier using PeptideIdentificationList
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _load_2(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  document_id ):
        """
        _load_2(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: String ) -> None
        Loads the identifications of an idXML file with identifier using PeptideIdentificationList
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(document_id, String), 'arg document_id wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((<String>document_id).inst.get()))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def load(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        Loads the identifications of an idXML file without identifier

        
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        Loads the identifications of an idXML file without identifier using PeptideIdentificationList

        
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: String ) -> None
          :noindex:
        
        Loads the identifications of an idXML file with identifier using PeptideIdentificationList
    
        """
        if (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)):
            return self._load_0(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)):
            return self._load_1(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and (isinstance(args[3], String)):
            return self._load_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _store_0(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  document_id ):
        """
        _store_0(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None
        Stores the data in an idXML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert (isinstance(document_id, str) or isinstance(document_id, bytes) or isinstance(document_id, String)), 'arg document_id wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((convString(document_id)).get()))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _store_1(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ):
        """
        _store_1(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _store_2(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  document_id ):
        """
        _store_2(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None
        Stores the data in an idXML file using PeptideIdentificationList
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert (isinstance(document_id, str) or isinstance(document_id, bytes) or isinstance(document_id, String)), 'arg document_id wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((convString(document_id)).get()))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def store(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None
          :noindex:
        
        Stores the data in an idXML file

        
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , document_id: Union[bytes, str, String] ) -> None
          :noindex:
        
        Stores the data in an idXML file using PeptideIdentificationList
    
        """
        if (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
            return self._store_0(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)):
            return self._store_1(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
            return self._store_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class IndexedMzMLFileLoader:
    """
    Cython implementation of _IndexedMzMLFileLoader

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLFileLoader.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        A class to load an indexedmzML file
        """
        self.inst = shared_ptr[_IndexedMzMLFileLoader](new _IndexedMzMLFileLoader())
    
    def load(self,  in_0 , OnDiscMSExperiment in_1 ):
        """
        load(self, in_0: Union[bytes, str, String] , in_1: OnDiscMSExperiment ) -> bool
        Load a file\n
        
        Tries to parse the file, success needs to be checked with the return value
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, OnDiscMSExperiment), 'arg in_1 wrong type'
    
    
        cdef bool _r = self.inst.get().load(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def _store_0(self,  in_0 , OnDiscMSExperiment in_1 ):
        """
        _store_0(self, in_0: Union[bytes, str, String] , in_1: OnDiscMSExperiment ) -> None
        Store a file from an on-disc data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, OnDiscMSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def _store_1(self,  in_0 , MSExperiment in_1 ):
        """
        _store_1(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        Store a file from an in-memory data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def store(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: store(self, in_0: Union[bytes, str, String] , in_1: OnDiscMSExperiment ) -> None
          :noindex:
        
        Store a file from an on-disc data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
        
        .. rubric:: Overload:
        .. py:function:: store(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None
          :noindex:
        
        Store a file from an in-memory data-structure
        
        
        :param filename: Filename determines where the file will be stored
        :param exp: MS data to be stored
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], OnDiscMSExperiment)):
            return self._store_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], MSExperiment)):
            return self._store_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        Returns the options for loading/storing
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
        return py_result
    
    def setOptions(self, PeakFileOptions in_0 ):
        """
        setOptions(self, in_0: PeakFileOptions ) -> None
        Returns the options for loading/storing
        """
        assert isinstance(in_0, PeakFileOptions), 'arg in_0 wrong type'
    
        self.inst.get().setOptions((deref(in_0.inst.get()))) 

cdef class IsobaricChannelInformation:
    """
    Cython implementation of _IsobaricChannelInformation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::IsobaricQuantitationMethod_1_1IsobaricChannelInformation.html>`_
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
        
            self.inst.get().id = (<int>id)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().id
            py_result = <int>_r
            return py_result
    
    property description:
        def __set__(self,  description):
        
            self.inst.get().description = deref((convString(description)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().description
            py_result = convOutputString(_r)
            return py_result
    
    property center:
        def __set__(self, double center):
        
            self.inst.get().center = (<double>center)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().center
            py_result = <double>_r
            return py_result
    
    property affected_channels:
        def __set__(self, list affected_channels):
            cdef libcpp_vector[int] v0 = affected_channels
            self.inst.get().affected_channels = v0
            
    
        def __get__(self):
            _r = self.inst.get().affected_channels
            cdef list py_result = _r
            return py_result
    
    def __copy__(self):
       cdef IsobaricChannelInformation rv = IsobaricChannelInformation.__new__(IsobaricChannelInformation)
       rv.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricChannelInformation rv = IsobaricChannelInformation.__new__(IsobaricChannelInformation)
       rv.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  name ,  id_ ,  description , double center , list affected_channels ):
        """
        _init_0(self, name: Union[bytes, str, String] , id_: int , description: Union[bytes, str, String] , center: float , affected_channels: List[int] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert (isinstance(description, str) or isinstance(description, bytes) or isinstance(description, String)), 'arg description wrong type'
        assert isinstance(center, float), 'arg center wrong type'
        assert isinstance(affected_channels, list) and all(isinstance(elemt_rec, int) for elemt_rec in affected_channels), 'arg affected_channels wrong type'
    
    
    
    
        cdef libcpp_vector[int] v4 = affected_channels
        self.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation(deref((convString(name)).get()), (<int>id_), deref((convString(description)).get()), (<double>center), v4))
        
    
    def _init_1(self, IsobaricChannelInformation in_0 ):
        """
        _init_1(self, in_0: IsobaricChannelInformation ) -> None
        """
        assert isinstance(in_0, IsobaricChannelInformation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricChannelInformation](new _IsobaricChannelInformation((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , id_: int , description: Union[bytes, str, String] , center: float , affected_channels: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricChannelInformation ) -> None
          :noindex:
    
        """
        if (len(args)==5) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int)) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], float)) and (isinstance(args[4], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[4])):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsobaricChannelInformation)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ItraqFourPlexQuantitationMethod:
    """
    Cython implementation of _ItraqFourPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqFourPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ItraqFourPlexQuantitationMethod rv = ItraqFourPlexQuantitationMethod.__new__(ItraqFourPlexQuantitationMethod)
       rv.inst = shared_ptr[_ItraqFourPlexQuantitationMethod](new _ItraqFourPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ItraqFourPlexQuantitationMethod rv = ItraqFourPlexQuantitationMethod.__new__(ItraqFourPlexQuantitationMethod)
       rv.inst = shared_ptr[_ItraqFourPlexQuantitationMethod](new _ItraqFourPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation
        """
        self.inst = shared_ptr[_ItraqFourPlexQuantitationMethod](new _ItraqFourPlexQuantitationMethod())
    
    def _init_1(self, ItraqFourPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: ItraqFourPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, ItraqFourPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ItraqFourPlexQuantitationMethod](new _ItraqFourPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        iTRAQ 4 plex quantitation to be used with the IsobaricQuantitation

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ItraqFourPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqFourPlexQuantitationMethod)):
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

cdef class MSDataCachedConsumer:
    """
    Cython implementation of _MSDataCachedConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataCachedConsumer.html>`_

    Transforming and cached writing consumer of MS data
    
    Is able to transform a spectrum on the fly while it is read using a
    function pointer that can be set on the object. The spectra is then
    cached to disk using the functions provided in CachedMzMLHandler.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self,  filename ):
        """
        _init_0(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_MSDataCachedConsumer](new _MSDataCachedConsumer(deref((convString(filename)).get())))
    
    def _init_1(self,  filename , bool clear ):
        """
        _init_1(self, filename: Union[bytes, str, String] , clear: bool ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(clear, pybool_t), 'arg clear wrong type'
    
    
        self.inst = shared_ptr[_MSDataCachedConsumer](new _MSDataCachedConsumer(deref((convString(filename)).get()), (<bool>clear)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , clear: bool ) -> None
          :noindex:
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], pybool_t)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def consumeSpectrum(self, MSSpectrum s ):
        """
        consumeSpectrum(self, s: MSSpectrum ) -> None
        Write a spectrum to the output file
        
        May delete data from spectrum (if clearData is set)
        """
        assert isinstance(s, MSSpectrum), 'arg s wrong type'
    
        self.inst.get().consumeSpectrum((deref(s.inst.get())))
    
    def consumeChromatogram(self, MSChromatogram c ):
        """
        consumeChromatogram(self, c: MSChromatogram ) -> None
        Write a chromatogram to the output file
        
        May delete data from chromatogram (if clearData is set)
        """
        assert isinstance(c, MSChromatogram), 'arg c wrong type'
    
        self.inst.get().consumeChromatogram((deref(c.inst.get())))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get())))
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms)) 

cdef class NonNegativeLeastSquaresSolver:
    """
    Cython implementation of _NonNegativeLeastSquaresSolver

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NonNegativeLeastSquaresSolver.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef NonNegativeLeastSquaresSolver rv = NonNegativeLeastSquaresSolver.__new__(NonNegativeLeastSquaresSolver)
       rv.inst = shared_ptr[_NonNegativeLeastSquaresSolver](new _NonNegativeLeastSquaresSolver(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NonNegativeLeastSquaresSolver rv = NonNegativeLeastSquaresSolver.__new__(NonNegativeLeastSquaresSolver)
       rv.inst = shared_ptr[_NonNegativeLeastSquaresSolver](new _NonNegativeLeastSquaresSolver(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NonNegativeLeastSquaresSolver](new _NonNegativeLeastSquaresSolver())
    
    def _init_1(self, NonNegativeLeastSquaresSolver in_0 ):
        """
        _init_1(self, in_0: NonNegativeLeastSquaresSolver ) -> None
        """
        assert isinstance(in_0, NonNegativeLeastSquaresSolver), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NonNegativeLeastSquaresSolver](new _NonNegativeLeastSquaresSolver((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NonNegativeLeastSquaresSolver ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NonNegativeLeastSquaresSolver)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def solve(self, MatrixDouble A , MatrixDouble b , MatrixDouble x ):
        """
        solve(self, A: MatrixDouble , b: MatrixDouble , x: MatrixDouble ) -> int
        """
        assert isinstance(A, MatrixDouble), 'arg A wrong type'
        assert isinstance(b, MatrixDouble), 'arg b wrong type'
        assert isinstance(x, MatrixDouble), 'arg x wrong type'
    
    
    
        cdef int _r = self.inst.get().solve((deref(A.inst.get())), (deref(b.inst.get())), (deref(x.inst.get())))
        py_result = <int>_r
        return py_result
    RETURN_STATUS = __RETURN_STATUS 

cdef class SplineInterpolatedPeaks:
    """
    Cython implementation of _SplineInterpolatedPeaks

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplineInterpolatedPeaks.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SplineInterpolatedPeaks rv = SplineInterpolatedPeaks.__new__(SplineInterpolatedPeaks)
       rv.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SplineInterpolatedPeaks rv = SplineInterpolatedPeaks.__new__(SplineInterpolatedPeaks)
       rv.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks(deref(self.inst.get())))
       return rv
    
    def _init_0(self, list mz , list intensity ):
        """
        _init_0(self, mz: List[float] , intensity: List[float] ) -> None
        """
        assert isinstance(mz, list) and all(isinstance(elemt_rec, float) for elemt_rec in mz), 'arg mz wrong type'
        assert isinstance(intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in intensity), 'arg intensity wrong type'
        cdef libcpp_vector[double] v0 = mz
        cdef libcpp_vector[double] v1 = intensity
        self.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks(v0, v1))
        
        
    
    def _init_1(self, MSSpectrum raw_spectrum ):
        """
        _init_1(self, raw_spectrum: MSSpectrum ) -> None
        """
        assert isinstance(raw_spectrum, MSSpectrum), 'arg raw_spectrum wrong type'
    
        self.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks((deref(raw_spectrum.inst.get()))))
    
    def _init_2(self, MSChromatogram raw_chromatogram ):
        """
        _init_2(self, raw_chromatogram: MSChromatogram ) -> None
        """
        assert isinstance(raw_chromatogram, MSChromatogram), 'arg raw_chromatogram wrong type'
    
        self.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks((deref(raw_chromatogram.inst.get()))))
    
    def _init_3(self, SplineInterpolatedPeaks in_0 ):
        """
        _init_3(self, in_0: SplineInterpolatedPeaks ) -> None
        """
        assert isinstance(in_0, SplineInterpolatedPeaks), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SplineInterpolatedPeaks](new _SplineInterpolatedPeaks((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, mz: List[float] , intensity: List[float] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, raw_spectrum: MSSpectrum ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, raw_chromatogram: MSChromatogram ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SplineInterpolatedPeaks ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSSpectrum)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], MSChromatogram)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], SplineInterpolatedPeaks)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPosMin(self):
        """
        getPosMin(self) -> float
        """
        cdef double _r = self.inst.get().getPosMin()
        py_result = <double>_r
        return py_result
    
    def getPosMax(self):
        """
        getPosMax(self) -> float
        """
        cdef double _r = self.inst.get().getPosMax()
        py_result = <double>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def getNavigator(self, double scaling ):
        """
        getNavigator(self, scaling: float ) -> SplineSpectrum_Navigator
        """
        assert isinstance(scaling, float), 'arg scaling wrong type'
    
        cdef _SplineSpectrum_Navigator * _r = new _SplineSpectrum_Navigator(self.inst.get().getNavigator((<double>scaling)))
        cdef SplineSpectrum_Navigator py_result = SplineSpectrum_Navigator.__new__(SplineSpectrum_Navigator)
        py_result.inst = shared_ptr[_SplineSpectrum_Navigator](_r)
        return py_result 

cdef class SplineSpectrum_Navigator:
    """
    Cython implementation of _SplineSpectrum_Navigator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplineSpectrum_Navigator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SplineSpectrum_Navigator](new _SplineSpectrum_Navigator())
    
    def _init_1(self, SplineSpectrum_Navigator in_0 ):
        """
        _init_1(self, in_0: SplineSpectrum_Navigator ) -> None
        """
        assert isinstance(in_0, SplineSpectrum_Navigator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SplineSpectrum_Navigator](new _SplineSpectrum_Navigator((deref(in_0.inst.get()))))
    
    def _init_2(self, list packages , double posMax , double scaling ):
        """
        _init_2(self, packages: List[SplinePackage] , posMax: float , scaling: float ) -> None
        """
        assert isinstance(packages, list) and all(isinstance(elemt_rec, SplinePackage) for elemt_rec in packages), 'arg packages wrong type'
        assert isinstance(posMax, float), 'arg posMax wrong type'
        assert isinstance(scaling, float), 'arg scaling wrong type'
        cdef libcpp_vector[_SplinePackage] * v0 = new libcpp_vector[_SplinePackage]()
        cdef SplinePackage item0
        for item0 in packages:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst = shared_ptr[_SplineSpectrum_Navigator](new _SplineSpectrum_Navigator(v0, (<double>posMax), (<double>scaling)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SplineSpectrum_Navigator ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, packages: List[SplinePackage] , posMax: float , scaling: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SplineSpectrum_Navigator)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, SplinePackage) for elemt_rec in args[0])) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def eval(self, double pos ):
        """
        eval(self, pos: float ) -> float
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().eval((<double>pos))
        py_result = <double>_r
        return py_result
    
    def getNextPos(self, double pos ):
        """
        getNextPos(self, pos: float ) -> float
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().getNextPos((<double>pos))
        py_result = <double>_r
        return py_result 

cdef class TraMLFile:
    """
    Cython implementation of _TraMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TraMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TraMLFile rv = TraMLFile.__new__(TraMLFile)
       rv.inst = shared_ptr[_TraMLFile](new _TraMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TraMLFile rv = TraMLFile.__new__(TraMLFile)
       rv.inst = shared_ptr[_TraMLFile](new _TraMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TraMLFile](new _TraMLFile())
    
    def _init_1(self, TraMLFile in_0 ):
        """
        _init_1(self, in_0: TraMLFile ) -> None
        """
        assert isinstance(in_0, TraMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TraMLFile](new _TraMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TraMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TraMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , TargetedExperiment id ):
        """
        load(self, filename: Union[bytes, str, String] , id: TargetedExperiment ) -> None
        Loads a map from a TraML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(id, TargetedExperiment), 'arg id wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(id.inst.get())))
    
    def store(self,  filename , TargetedExperiment id ):
        """
        store(self, filename: Union[bytes, str, String] , id: TargetedExperiment ) -> None
        Stores a map in a TraML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(id, TargetedExperiment), 'arg id wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(id.inst.get())))
    
    def isSemanticallyValid(self,  filename , list errors , list warnings ):
        """
        isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool
        Checks if a file is valid with respect to the mapping file and the controlled vocabulary
        
        :param filename: File name of the file to be checked
        :param errors: Errors during the validation are returned in this output parameter
        :param warnings: Warnings during the validation are returned in this output parameter
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
        cdef bool _r = self.inst.get().isSemanticallyValid(deref((convString(filename)).get()), deref(v1), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        warnings[:] = replace
        del v2
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        errors[:] = replace
        del v1
        py_result = <bool>_r
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
 
 

cdef class DPosition1:

    def __dealloc__(self):
        self.inst.reset()

    def __init__(self, *args , **kwargs):
        if not args:
             self._init_0(*args)
        elif (len(args)==1):
             self._init_1(*args)
        else:
             raise Exception('can not handle type of %s' % (args,))
 
    def _init_0(self):
        self.inst = shared_ptr[_DPosition1](new _DPosition1())

    def _init_1(self, double a):
        self.inst = shared_ptr[_DPosition1](new _DPosition1(a))

    def __getitem__(self, ix):
        if ix != 0:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[0]


cdef class DPosition2:

    def __dealloc__(self):
        self.inst.reset()


    def __init__(self, *args , **kwargs):
        if not args:
             self._init_0(*args)
        elif (len(args)==2):
             self._init_1(*args)
        else:
             raise Exception('can not handle type of %s' % (args,))
 
    def _init_0(self):
        self.inst = shared_ptr[_DPosition2](new _DPosition2())

    def _init_1(self, double a, double b):
        self.inst = shared_ptr[_DPosition2](new _DPosition2(a, b))

    def __getitem__(self, ix):
        if ix != 0 and ix != 1:
            raise IndexError("invalid index %d" % ix)
        return deref(self.inst.get())[ix] 

## Add this to all modules
class Interfaces:

    BinaryDataArray = _Interfaces_BinaryDataArray
    Spectrum = _Interfaces_Spectrum
    Chromatogram = _Interfaces_Chromatogram 

# Add this to all modules
PeakSpectrum = MSSpectrum

PeakMap = MSExperiment 
