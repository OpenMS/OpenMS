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

cdef class DataType:
    None
    STRING_VALUE = 0
    INT_VALUE = 1
    DOUBLE_VALUE = 2
    STRING_LIST = 3
    INT_LIST = 4
    DOUBLE_LIST = 5
    EMPTY_VALUE = 6

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class UnitType:
    None
    UNIT_ONTOLOGY = 0
    MS_ONTOLOGY = 1
    OTHER = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ConfidenceScoring:
    """
    Cython implementation of _ConfidenceScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConfidenceScoring.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ConfidenceScoring rv = ConfidenceScoring.__new__(ConfidenceScoring)
       rv.inst = shared_ptr[_ConfidenceScoring](new _ConfidenceScoring(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ConfidenceScoring rv = ConfidenceScoring.__new__(ConfidenceScoring)
       rv.inst = shared_ptr[_ConfidenceScoring](new _ConfidenceScoring(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ConfidenceScoring](new _ConfidenceScoring())
    
    def _init_1(self, ConfidenceScoring in_0 ):
        """
        _init_1(self, in_0: ConfidenceScoring ) -> None
        """
        assert isinstance(in_0, ConfidenceScoring), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ConfidenceScoring](new _ConfidenceScoring((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ConfidenceScoring ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConfidenceScoring)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def initialize(self, TargetedExperiment targeted ,  n_decoys ,  n_transitions , TransformationDescription trafo ):
        """
        initialize(self, targeted: TargetedExperiment , n_decoys: int , n_transitions: int , trafo: TransformationDescription ) -> None
        """
        assert isinstance(targeted, TargetedExperiment), 'arg targeted wrong type'
        assert isinstance(n_decoys, int) and n_decoys >= 0, 'arg n_decoys wrong type'
        assert isinstance(n_transitions, int) and n_transitions >= 0, 'arg n_transitions wrong type'
        assert isinstance(trafo, TransformationDescription), 'arg trafo wrong type'
    
    
    
    
        self.inst.get().initialize((deref(targeted.inst.get())), (<size_t>n_decoys), (<size_t>n_transitions), (deref(trafo.inst.get())))
    
    def initializeGlm(self, double intercept , double rt_coef , double int_coef ):
        """
        initializeGlm(self, intercept: float , rt_coef: float , int_coef: float ) -> None
        """
        assert isinstance(intercept, float), 'arg intercept wrong type'
        assert isinstance(rt_coef, float), 'arg rt_coef wrong type'
        assert isinstance(int_coef, float), 'arg int_coef wrong type'
    
    
    
        self.inst.get().initializeGlm((<double>intercept), (<double>rt_coef), (<double>int_coef))
    
    def scoreMap(self, FeatureMap map ):
        """
        scoreMap(self, map: FeatureMap ) -> None
        Score a feature map -> make sure the class is properly initialized
        """
        assert isinstance(map, FeatureMap), 'arg map wrong type'
    
        self.inst.get().scoreMap((deref(map.inst.get()))) 

cdef class ConsensusIDAlgorithmRanks:
    """
    Cython implementation of _ConsensusIDAlgorithmRanks

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmRanks.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmRanks](new _ConsensusIDAlgorithmRanks())
    
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

cdef class ConsensusMapNormalizerAlgorithmQuantile:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmQuantile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusMapNormalizerAlgorithmQuantile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusMapNormalizerAlgorithmQuantile](new _ConsensusMapNormalizerAlgorithmQuantile())
    
    def normalizeMaps(self, ConsensusMap input_map ):
        """
        normalizeMaps(self, input_map: ConsensusMap ) -> None
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
    
        self.inst.get().normalizeMaps((deref(input_map.inst.get())))
    
    def resample(self, list data_in , list data_out ,  n_resampling_points ):
        """
        resample(self, data_in: List[float] , data_out: List[float] , n_resampling_points: int ) -> None
        Resamples data_in and writes the results to data_out
        """
        assert isinstance(data_in, list) and all(isinstance(elemt_rec, float) for elemt_rec in data_in), 'arg data_in wrong type'
        assert isinstance(data_out, list) and all(isinstance(elemt_rec, float) for elemt_rec in data_out), 'arg data_out wrong type'
        assert isinstance(n_resampling_points, int), 'arg n_resampling_points wrong type'
        cdef libcpp_vector[double] v0 = data_in
        cdef libcpp_vector[double] v1 = data_out
    
        self.inst.get().resample(v0, v1, (<unsigned int>n_resampling_points))
        data_out[:] = v1
        data_in[:] = v0
    
    def extractIntensityVectors(self, ConsensusMap map_ , list out_intensities ):
        """
        extractIntensityVectors(self, map_: ConsensusMap , out_intensities: List[List[float]] ) -> None
        Extracts the intensities of the features of the different maps
        """
        assert isinstance(map_, ConsensusMap), 'arg map_ wrong type'
        assert isinstance(out_intensities, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, float) for elemt_rec_rec in elemt_rec) for elemt_rec in out_intensities), 'arg out_intensities wrong type'
    
        cdef libcpp_vector[libcpp_vector[double]] v1 = out_intensities
        self.inst.get().extractIntensityVectors((deref(map_.inst.get())), v1)
        out_intensities[:] = v1
    
    def setNormalizedIntensityValues(self, list feature_ints , ConsensusMap map_ ):
        """
        setNormalizedIntensityValues(self, feature_ints: List[List[float]] , map_: ConsensusMap ) -> None
        Writes the intensity values in feature_ints to the corresponding features in map
        """
        assert isinstance(feature_ints, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, float) for elemt_rec_rec in elemt_rec) for elemt_rec in feature_ints), 'arg feature_ints wrong type'
        assert isinstance(map_, ConsensusMap), 'arg map_ wrong type'
        cdef libcpp_vector[libcpp_vector[double]] v0 = feature_ints
    
        self.inst.get().setNormalizedIntensityValues(v0, (deref(map_.inst.get())))
        feature_ints[:] = v0 

cdef class DataValue:
    """
    Cython implementation of _DataValue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataValue.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DataValue rv = DataValue.__new__(DataValue)
       rv.inst = shared_ptr[_DataValue](new _DataValue(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DataValue rv = DataValue.__new__(DataValue)
       rv.inst = shared_ptr[_DataValue](new _DataValue(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DataValue](new _DataValue())
    
    def _init_1(self, DataValue in_0 ):
        """
        _init_1(self, in_0: DataValue ) -> None
        """
        assert isinstance(in_0, DataValue), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue((deref(in_0.inst.get()))))
    
    def _init_2(self, bytes in_0 ):
        """
        _init_2(self, in_0: bytes ) -> None
        """
        assert isinstance(in_0, bytes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue((<char *>in_0)))
    
    def _init_3(self,  in_0 ):
        """
        _init_3(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue(deref((convString(in_0)).get())))
    
    def _init_4(self,  in_0 ):
        """
        _init_4(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue((<int>in_0)))
    
    def _init_5(self, double in_0 ):
        """
        _init_5(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue((<double>in_0)))
    
    def _init_6(self, list in_0 ):
        """
        _init_6(self, in_0: List[bytes] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(li, bytes) for li in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in in_0:
           v0.push_back(_String(<char *>item0))
        self.inst = shared_ptr[_DataValue](new _DataValue(deref(v0)))
        del v0
    
    def _init_7(self, list in_0 ):
        """
        _init_7(self, in_0: List[int] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(li, int) for li in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[int] _v0 = in_0
        cdef _IntList v0 = _IntList(_v0)
        self.inst = shared_ptr[_DataValue](new _DataValue((v0)))
    
    def _init_8(self, list in_0 ):
        """
        _init_8(self, in_0: List[float] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(li, float) for li in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] _v0 = in_0
        cdef _DoubleList v0 = _DoubleList(_v0)
        self.inst = shared_ptr[_DataValue](new _DataValue((v0)))
    
    def _init_9(self,  in_0 ):
        """
        _init_9(self, in_0: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
        """
        assert isinstance(in_0, (int, float, list, bytes, str)), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataValue](new _DataValue(deref(ParamValue(in_0).inst.get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DataValue ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: bytes ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[bytes] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: List[float] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DataValue)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], bytes)):
             self._init_2(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
             self._init_4(*args)
        elif (len(args)==1) and (isinstance(args[0], float)):
             self._init_5(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])):
             self._init_6(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, int) for li in args[0])):
             self._init_7(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(li, float) for li in args[0])):
             self._init_8(*args)
        elif (len(args)==1) and (isinstance(args[0], (int, float, list, bytes, str))):
             self._init_9(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toInt(self):
        cdef int _r = <int>(deref(self.inst.get()))
        py_res = <int>_r
        return py_res
    
    def toString(self):
        cdef _String _r = <_String>(deref(self.inst.get()))
        py_res = convOutputString(_r)
        return py_res
    
    def toDouble(self):
        cdef double _r = <double>(deref(self.inst.get()))
        py_res = <double>_r
        return py_res
    
    def toStringList(self):
        """
        toStringList(self) -> List[bytes]
        """
        _r = self.inst.get().toStringList()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def toDoubleList(self):
        """
        toDoubleList(self) -> List[float]
        """
        _r = self.inst.get().toDoubleList()
        cdef list py_result = _r
        return py_result
    
    def toIntList(self):
        """
        toIntList(self) -> List[int]
        """
        _r = self.inst.get().toIntList()
        cdef list py_result = _r
        return py_result
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def toBool(self):
        """
        toBool(self) -> bool
        """
        cdef bool _r = self.inst.get().toBool()
        py_result = <bool>_r
        return py_result
    
    def valueType(self):
        """
        valueType(self) -> int
        """
        cdef _DataType _r = self.inst.get().valueType()
        py_result = <int>_r
        return py_result
    
    def isEmpty(self):
        """
        isEmpty(self) -> int
        """
        cdef int _r = self.inst.get().isEmpty()
        py_result = <int>_r
        return py_result
    
    def getUnitType(self):
        """
        getUnitType(self) -> int
        """
        cdef _UnitType _r = self.inst.get().getUnitType()
        py_result = <int>_r
        return py_result
    
    def setUnitType(self, int u ):
        """
        setUnitType(self, u: int ) -> None
        """
        assert u in [0, 1, 2], 'arg u wrong type'
    
        self.inst.get().setUnitType((<_UnitType>u))
    
    def hasUnit(self):
        """
        hasUnit(self) -> bool
        """
        cdef bool _r = self.inst.get().hasUnit()
        py_result = <bool>_r
        return py_result
    
    def getUnit(self):
        """
        getUnit(self) -> int
        """
        cdef int _r = self.inst.get().getUnit()
        py_result = <int>_r
        return py_result
    
    def setUnit(self,  unit_id ):
        """
        setUnit(self, unit_id: int ) -> None
        """
        assert isinstance(unit_id, int), 'arg unit_id wrong type'
    
        self.inst.get().setUnit((<int>unit_id))
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result 

cdef class LabeledPairFinder:
    """
    Cython implementation of _LabeledPairFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LabeledPairFinder.html>`_
      -- Inherits from ['BaseGroupFinder']

    The LabeledPairFinder allows the matching of labeled features (features with a fixed distance)
    
    Finds feature pairs that have a defined distance in RT and m/z in the same map
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_LabeledPairFinder](new _LabeledPairFinder())
    
    def run(self, list input_maps , ConsensusMap result_map ):
        """
        run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None
        Runs the LabeledPairFinder algorithm
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

cdef class MRMDecoy:
    """
    Cython implementation of _MRMDecoy

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMDecoy.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMDecoy rv = MRMDecoy.__new__(MRMDecoy)
       rv.inst = shared_ptr[_MRMDecoy](new _MRMDecoy(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMDecoy rv = MRMDecoy.__new__(MRMDecoy)
       rv.inst = shared_ptr[_MRMDecoy](new _MRMDecoy(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMDecoy](new _MRMDecoy())
    
    def _init_1(self, MRMDecoy in_0 ):
        """
        _init_1(self, in_0: MRMDecoy ) -> None
        """
        assert isinstance(in_0, MRMDecoy), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMDecoy](new _MRMDecoy((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMDecoy ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMDecoy)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def generateDecoys(self, TargetedExperiment exp , TargetedExperiment dec ,  method , double aim_decoy_fraction , bool switchKR ,  decoy_tag ,  max_attempts , double identity_threshold , double precursor_mz_shift , double product_mz_shift , double product_mz_threshold , list fragment_types , list fragment_charges , bool enable_specific_losses , bool enable_unspecific_losses ,  round_decPow ):
        """
        generateDecoys(self, exp: TargetedExperiment , dec: TargetedExperiment , method: Union[bytes, str, String] , aim_decoy_fraction: float , switchKR: bool , decoy_tag: Union[bytes, str, String] , max_attempts: int , identity_threshold: float , precursor_mz_shift: float , product_mz_shift: float , product_mz_threshold: float , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None
        Generate decoys from a TargetedExperiment
        
        Will generate decoy peptides for each target peptide provided in exp and
        write them into the decoy experiment
        
        Valid methods: shuffle, reverse, pseudo-reverse
        
        If theoretical is true, the target transitions will be returned but their
        masses will be adjusted to match the theoretical value of the fragment ion
        that is the most likely explanation for the product
        
        `mz_threshold` is used for the matching of theoretical ion series to the observed one
        
        To generate decoys with different precursor mass, use the "switchKR" flag
        which switches terminal K/R (switches K to R and R to K). This generates
        different precursor m/z and ensures that the y ion series has a different
        mass. For a description of the procedure, see (supplemental material)
        
        Bruderer et al. Mol Cell Proteomics. 2017. 10.1074/mcp.RA117.000314.
        """
        assert isinstance(exp, TargetedExperiment), 'arg exp wrong type'
        assert isinstance(dec, TargetedExperiment), 'arg dec wrong type'
        assert (isinstance(method, str) or isinstance(method, bytes) or isinstance(method, String)), 'arg method wrong type'
        assert isinstance(aim_decoy_fraction, float), 'arg aim_decoy_fraction wrong type'
        assert isinstance(switchKR, pybool_t), 'arg switchKR wrong type'
        assert (isinstance(decoy_tag, str) or isinstance(decoy_tag, bytes) or isinstance(decoy_tag, String)), 'arg decoy_tag wrong type'
        assert isinstance(max_attempts, int), 'arg max_attempts wrong type'
        assert isinstance(identity_threshold, float), 'arg identity_threshold wrong type'
        assert isinstance(precursor_mz_shift, float), 'arg precursor_mz_shift wrong type'
        assert isinstance(product_mz_shift, float), 'arg product_mz_shift wrong type'
        assert isinstance(product_mz_threshold, float), 'arg product_mz_threshold wrong type'
        assert isinstance(fragment_types, list) and all(isinstance(i, bytes) for i in fragment_types), 'arg fragment_types wrong type'
        assert isinstance(fragment_charges, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in fragment_charges), 'arg fragment_charges wrong type'
        assert isinstance(enable_specific_losses, pybool_t), 'arg enable_specific_losses wrong type'
        assert isinstance(enable_unspecific_losses, pybool_t), 'arg enable_unspecific_losses wrong type'
        assert isinstance(round_decPow, int), 'arg round_decPow wrong type'
    
    
    
    
    
    
    
    
    
    
    
        cdef libcpp_vector[_String] * v11 = new libcpp_vector[_String]()
        cdef bytes item11
        for item11 in fragment_types:
           v11.push_back(_String(<char *>item11))
        cdef libcpp_vector[size_t] v12 = fragment_charges
    
    
    
        self.inst.get().generateDecoys((deref(exp.inst.get())), (deref(dec.inst.get())), deref((convString(method)).get()), (<double>aim_decoy_fraction), (<bool>switchKR), deref((convString(decoy_tag)).get()), (<int>max_attempts), (<double>identity_threshold), (<double>precursor_mz_shift), (<double>product_mz_shift), (<double>product_mz_threshold), deref(v11), v12, (<bool>enable_specific_losses), (<bool>enable_unspecific_losses), (<int>round_decPow))
        
        del v11
    
    def findFixedResidues(self,  sequence , bool keepN , bool keepC ,  keep_const_pattern ):
        """
        findFixedResidues(self, sequence: Union[bytes, str, String] , keepN: bool , keepC: bool , keep_const_pattern: Union[bytes, str, String] ) -> List[int]
        Find all residues in a sequence that should not be reversed / shuffled
        
        
        :param sequence: The amino acid sequence
        :param keepN: Whether to keep N terminus constant
        :param keepC: Whether to keep C terminus constant
        :param keep_const_pattern: A string containing the AA to not change (e.g. 'KRP')
        """
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
        assert isinstance(keepN, pybool_t), 'arg keepN wrong type'
        assert isinstance(keepC, pybool_t), 'arg keepC wrong type'
        assert (isinstance(keep_const_pattern, str) or isinstance(keep_const_pattern, bytes) or isinstance(keep_const_pattern, String)), 'arg keep_const_pattern wrong type'
    
    
    
    
        _r = self.inst.get().findFixedResidues(deref((convString(sequence)).get()), (<bool>keepN), (<bool>keepC), deref((convString(keep_const_pattern)).get()))
        cdef list py_result = _r
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

cdef class MapAlignmentEvaluationAlgorithmPrecision:
    """
    Cython implementation of _MapAlignmentEvaluationAlgorithmPrecision

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentEvaluationAlgorithmPrecision.html>`_
      -- Inherits from ['MapAlignmentEvaluationAlgorithm']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MapAlignmentEvaluationAlgorithmPrecision](new _MapAlignmentEvaluationAlgorithmPrecision()) 

cdef class MzMLFile:
    """
    Cython implementation of _MzMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLFile.html>`_
      -- Inherits from ['ProgressLogger']

    File adapter for MzML files
    
    Provides methods to load and store MzML files.
    PeakFileOptions allow to load a reduced subset of the data into an MSExperiment.
    
    See help(MSExperiment) how data is stored after loading.
    See help(PeakFileOptions) for available options.
    
    Usage:
    
    .. code-block:: python
    
      exp = MSExperiment()
      MzMLFile().load("test.mzML", exp)
      spec = []
      for s in exp.getSpectra():
        if s.getMSLevel() != 1:
          spec.append(s)
      exp.setSpectra(spec)
      MzMLFile().store("filtered.mzML", exp)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzMLFile rv = MzMLFile.__new__(MzMLFile)
       rv.inst = shared_ptr[_MzMLFile](new _MzMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzMLFile rv = MzMLFile.__new__(MzMLFile)
       rv.inst = shared_ptr[_MzMLFile](new _MzMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzMLFile](new _MzMLFile())
    
    def _init_1(self, MzMLFile in_0 ):
        """
        _init_1(self, in_0: MzMLFile ) -> None
        """
        assert isinstance(in_0, MzMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzMLFile](new _MzMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment in_1 ):
        """
        load(self, filename: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        Loads from an MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions)
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(in_1.inst.get())))
    
    def store(self,  filename , MSExperiment in_1 ):
        """
        store(self, filename: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        Stores a MSExperiment in an MzML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(in_1.inst.get())))
    
    def storeBuffer(self,  output , MSExperiment exp ):
        """
        storeBuffer(self, output: String , exp: MSExperiment ) -> None
        Stores a map in an output string
        
        
        :param output: An empty string to store the result
        :param exp: Has to be an MSExperiment
        """
        assert isinstance(output, String), 'arg output wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().storeBuffer(deref((<String>output).inst.get()), (deref(exp.inst.get())))
    
    def loadBuffer(self,  input , MSExperiment exp ):
        """
        loadBuffer(self, input: Union[bytes, str, String] , exp: MSExperiment ) -> None
        Loads a map from a MzML file stored in a buffer (in memory)
        
        
        :param buffer: The buffer with the data (i.e. string with content of an mzML file)
        :param exp: Is an MSExperiment
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(input, str) or isinstance(input, bytes) or isinstance(input, String)), 'arg input wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().loadBuffer(deref((convString(input)).get()), (deref(exp.inst.get())))
    
    def getOptions(self):
        """
        getOptions(self) -> PeakFileOptions
        """
        cdef _PeakFileOptions * _r = new _PeakFileOptions(self.inst.get().getOptions())
        cdef PeakFileOptions py_result = PeakFileOptions.__new__(PeakFileOptions)
        py_result.inst = shared_ptr[_PeakFileOptions](_r)
        return py_result
    
    def setOptions(self, PeakFileOptions in_0 ):
        """
        setOptions(self, in_0: PeakFileOptions ) -> None
        Set PeakFileOptions to perform filtering during loading. E.g., to load only MS1 spectra or meta data only
        """
        assert isinstance(in_0, PeakFileOptions), 'arg in_0 wrong type'
    
        self.inst.get().setOptions((deref(in_0.inst.get())))
    
    def isSemanticallyValid(self,  filename , list errors , list warnings ):
        """
        isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool
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
    

    def transform(self, *args):
        if (len(args)==2):
             self._transform_1(*args)
        elif (len(args)==3):
             self._transform_2(*args)
        elif (len(args)==4):
             self._transform_3(*args)
        elif (len(args)==5):
             self._transform_4(*args)
        else:
           raise Exception('can not handle type of %s' % (args,))

        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *) except + nogil  # wrap-ignore
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, bool skip_full_count, bool skip_first_pass) except + nogil 
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e) except + nogil  # wrap-ignore
        # void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e, bool skip_full_count, bool skip_first_pass) except + nogil  # wrap-ignore

    def _transform_4(self, path, transformer, MSExperiment exp, bool skip_full_count, bool skip_first_pass):
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzml,
                                             _wrap_MSChromatogram_mzml,
                                             _wrap_ExperimentalSettings_mzml)

        try:
            self.inst.get().transform(deref((convString(path)).get()), consumer, deref(exp.inst.get()), skip_full_count, skip_first_pass)
        finally:
            del consumer

    def _transform_2(self, path, transformer, MSExperiment exp):
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzml,
                                             _wrap_MSChromatogram_mzml,
                                             _wrap_ExperimentalSettings_mzml)

        try:
            self.inst.get().transform(deref((convString(path)).get()), consumer, deref(exp.inst.get()) )
        finally:
            del consumer

    def _transform_3(self, path, transformer, bool skip_full_count, bool skip_first_pass):
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzml,
                                             _wrap_MSChromatogram_mzml,
                                             _wrap_ExperimentalSettings_mzml)

        try:
            self.inst.get().transform(deref((convString(path)).get()), consumer, skip_full_count, skip_first_pass)
        finally:
            del consumer


    def _transform_1(self, path, transformer):
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

        assert hasattr(transformer, "consumeSpectrum"), "expected method consumeSpectrum"
        assert hasattr(transformer, "consumeChromatogram"), "expected method consumeChromatogram"
        assert hasattr(transformer, "setExpectedSize"), "expected method setExpectedSize"
        assert hasattr(transformer, "setExperimentalSettings"), "expected method setExperimentalSettings"
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzml,
                                             _wrap_MSChromatogram_mzml,
                                             _wrap_ExperimentalSettings_mzml)

        try:
            self.inst.get().transform(deref((convString(path)).get()), consumer)
        finally:
            del consumer


cdef _wrap_MSSpectrum_mzml(const _MSSpectrum & _spec):
    cdef MSSpectrum spec = MSSpectrum.__new__(MSSpectrum)
    spec.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(_spec))
    return spec


cdef _wrap_MSChromatogram_mzml(const _MSChromatogram & _chromo):
    cdef MSChromatogram chromo = MSChromatogram.__new__(MSChromatogram)
    chromo.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(_chromo))
    return chromo


cdef _wrap_ExperimentalSettings_mzml(const _ExperimentalSettings & _exp):
    cdef ExperimentalSettings exp = ExperimentalSettings.__new__(ExperimentalSettings)
    exp.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(_exp))
    return exp 

cdef class MzMLSpectrumDecoder:
    """
    Cython implementation of _MzMLSpectrumDecoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzMLSpectrumDecoder.html>`_

    A class to decode input strings that contain an mzML chromatogram or spectrum tag
    
    It uses xercesc to parse a string containing either a exactly one mzML
    spectrum or chromatogram (from <chromatogram> to </chromatogram> or
    <spectrum> to </spectrum> tag). It returns the data contained in the
    binaryDataArray for Intensity / mass-to-charge or Intensity / time
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzMLSpectrumDecoder rv = MzMLSpectrumDecoder.__new__(MzMLSpectrumDecoder)
       rv.inst = shared_ptr[_MzMLSpectrumDecoder](new _MzMLSpectrumDecoder(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzMLSpectrumDecoder rv = MzMLSpectrumDecoder.__new__(MzMLSpectrumDecoder)
       rv.inst = shared_ptr[_MzMLSpectrumDecoder](new _MzMLSpectrumDecoder(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzMLSpectrumDecoder](new _MzMLSpectrumDecoder())
    
    def _init_1(self, MzMLSpectrumDecoder in_0 ):
        """
        _init_1(self, in_0: MzMLSpectrumDecoder ) -> None
        """
        assert isinstance(in_0, MzMLSpectrumDecoder), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzMLSpectrumDecoder](new _MzMLSpectrumDecoder((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzMLSpectrumDecoder ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzMLSpectrumDecoder)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def domParseChromatogram(self,  in_ , _Interfaces_Chromatogram cptr ):
        """
        domParseChromatogram(self, in_: Union[bytes, str, String] , cptr: _Interfaces_Chromatogram ) -> None
        Extract data from a string which contains a full mzML chromatogram
        
        Extracts data from the input string which is expected to contain exactly
        one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
        function will extract the contained binaryDataArray and provide the
        result as Chromatogram
        
        
        :param in: Input string containing the raw XML
        :param cptr: Resulting chromatogram
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(cptr, _Interfaces_Chromatogram), 'arg cptr wrong type'
    
        cdef shared_ptr[_Chromatogram] input_cptr = cptr.inst
        self.inst.get().domParseChromatogram(deref((convString(in_)).get()), input_cptr)
        cptr.inst = input_cptr
    
    def domParseSpectrum(self,  in_ , _Interfaces_Spectrum cptr ):
        """
        domParseSpectrum(self, in_: Union[bytes, str, String] , cptr: _Interfaces_Spectrum ) -> None
        Extract data from a string which contains a full mzML spectrum
        
        Extracts data from the input string which is expected to contain exactly
        one <spectrum> tag (from <spectrum> to </spectrum>). This function will
        extract the contained binaryDataArray and provide the result as Spectrum
        
        
        :param in: Input string containing the raw XML
        :param cptr: Resulting spectrum
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(cptr, _Interfaces_Spectrum), 'arg cptr wrong type'
    
        cdef shared_ptr[_Spectrum] input_cptr = cptr.inst
        self.inst.get().domParseSpectrum(deref((convString(in_)).get()), input_cptr)
        cptr.inst = input_cptr
    
    def setSkipXMLChecks(self, bool only ):
        """
        setSkipXMLChecks(self, only: bool ) -> None
        Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead
        """
        assert isinstance(only, pybool_t), 'arg only wrong type'
    
        self.inst.get().setSkipXMLChecks((<bool>only)) 

cdef class ProteinInference:
    """
    Cython implementation of _ProteinInference

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteinInference.html>`_

    [experimental class] given a peptide quantitation, infer corresponding protein quantities
    
    Infers protein ratios from peptide ratios (currently using unique peptides only).
    Use the IDMapper class to add protein and peptide information to a
    quantitative ConsensusMap prior to this step
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ProteinInference rv = ProteinInference.__new__(ProteinInference)
       rv.inst = shared_ptr[_ProteinInference](new _ProteinInference(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteinInference rv = ProteinInference.__new__(ProteinInference)
       rv.inst = shared_ptr[_ProteinInference](new _ProteinInference(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteinInference](new _ProteinInference())
    
    def _init_1(self, ProteinInference in_0 ):
        """
        _init_1(self, in_0: ProteinInference ) -> None
        """
        assert isinstance(in_0, ProteinInference), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteinInference](new _ProteinInference((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteinInference ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteinInference)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def infer(self, ConsensusMap consensus_map ,  reference_map ):
        """
        infer(self, consensus_map: ConsensusMap , reference_map: int ) -> None
        Given a peptide quantitation, infer corresponding protein quantities
        
        Infers protein ratios from peptide ratios (currently using unique peptides only).
        Use the IDMapper class to add protein and peptide information to a
        quantitative ConsensusMap prior to this step
        
        
        :param consensus_map: Peptide quantitation with ProteinIdentifications attached, where protein quantitation will be attached
        :param reference_map: Index of (iTRAQ) reference channel within the consensus map
        """
        assert isinstance(consensus_map, ConsensusMap), 'arg consensus_map wrong type'
        assert isinstance(reference_map, int), 'arg reference_map wrong type'
    
    
        self.inst.get().infer((deref(consensus_map.inst.get())), (<unsigned int>reference_map)) 

cdef class SpectrumAccessOpenMSInMemory:
    """
    Cython implementation of _SpectrumAccessOpenMSInMemory

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMSInMemory.html>`_
      -- Inherits from ['ISpectrumAccess']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAccessOpenMSInMemory rv = SpectrumAccessOpenMSInMemory.__new__(SpectrumAccessOpenMSInMemory)
       rv.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAccessOpenMSInMemory rv = SpectrumAccessOpenMSInMemory.__new__(SpectrumAccessOpenMSInMemory)
       rv.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        pass
    
    def _init_1(self, SpectrumAccessOpenMS in_0 ):
        """
        _init_1(self, in_0: SpectrumAccessOpenMS ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMS), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory((deref(in_0.inst.get()))))
    
    def _init_2(self, SpectrumAccessOpenMSCached in_0 ):
        """
        _init_2(self, in_0: SpectrumAccessOpenMSCached ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMSCached), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory((deref(in_0.inst.get()))))
    
    def _init_3(self, SpectrumAccessOpenMSInMemory in_0 ):
        """
        _init_3(self, in_0: SpectrumAccessOpenMSInMemory ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMSInMemory), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory((deref(in_0.inst.get()))))
    
    def _init_4(self, SpectrumAccessQuadMZTransforming in_0 ):
        """
        _init_4(self, in_0: SpectrumAccessQuadMZTransforming ) -> None
        """
        assert isinstance(in_0, SpectrumAccessQuadMZTransforming), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessOpenMSInMemory](new _SpectrumAccessOpenMSInMemory((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMS ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMSCached ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMSInMemory ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessQuadMZTransforming ) -> None
          :noindex:
    
        """
        if kwargs.get("__createUnsafeObject__") is True:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMS)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMSCached)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMSInMemory)):
             self._init_3(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessQuadMZTransforming)):
             self._init_4(*args)
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

cdef class SpectrumLookup:
    """
    Cython implementation of _SpectrumLookup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumLookup.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property rt_tolerance:
        def __set__(self, double rt_tolerance):
        
            self.inst.get().rt_tolerance = (<double>rt_tolerance)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_tolerance
            py_result = <double>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_SpectrumLookup](new _SpectrumLookup())
    
    def empty(self):
        """
        empty(self) -> bool
        Check if any spectra were set
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def readSpectra(self, MSExperiment spectra ,  scan_regexp ):
        """
        readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] ) -> None
        Read and index spectra for later look-up
        
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
        """
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert (isinstance(scan_regexp, str) or isinstance(scan_regexp, bytes) or isinstance(scan_regexp, String)), 'arg scan_regexp wrong type'
    
    
        self.inst.get().readSpectra((deref(spectra.inst.get())), deref((convString(scan_regexp)).get()))
    
    def findByRT(self, double rt ):
        """
        findByRT(self, rt: float ) -> int
        Look up spectrum by retention time (RT)
        
        :param rt: Retention time to look up
        :returns: Index of the spectrum that matched
        """
        assert isinstance(rt, float), 'arg rt wrong type'
    
        cdef size_t _r = self.inst.get().findByRT((<double>rt))
        py_result = <size_t>_r
        return py_result
    
    def findByNativeID(self,  native_id ):
        """
        findByNativeID(self, native_id: Union[bytes, str, String] ) -> int
        Look up spectrum by native ID
        
        :param native_id: Native ID to look up
        :returns: Index of the spectrum that matched
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
    
        cdef size_t _r = self.inst.get().findByNativeID(deref((convString(native_id)).get()))
        py_result = <size_t>_r
        return py_result
    
    def findByIndex(self,  index , bool count_from_one ):
        """
        findByIndex(self, index: int , count_from_one: bool ) -> int
        Look up spectrum by index (position in the vector of spectra)
        
        :param index: Index to look up
        :param count_from_one: Do indexes start counting at one (default zero)?
        :returns: Index of the spectrum that matched
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(count_from_one, pybool_t), 'arg count_from_one wrong type'
    
    
        cdef size_t _r = self.inst.get().findByIndex((<size_t>index), (<bool>count_from_one))
        py_result = <size_t>_r
        return py_result
    
    def findByScanNumber(self,  scan_number ):
        """
        findByScanNumber(self, scan_number: int ) -> int
        Look up spectrum by scan number (extracted from the native ID)
        
        :param scan_number: Scan number to look up
        :returns: Index of the spectrum that matched
        """
        assert isinstance(scan_number, int) and scan_number >= 0, 'arg scan_number wrong type'
    
        cdef size_t _r = self.inst.get().findByScanNumber((<size_t>scan_number))
        py_result = <size_t>_r
        return py_result
    
    def findByReference(self,  spectrum_ref ):
        """
        findByReference(self, spectrum_ref: Union[bytes, str, String] ) -> int
        Look up spectrum by reference
        
        :param spectrum_ref: Spectrum reference to parse
        :returns: Index of the spectrum that matched
        """
        assert (isinstance(spectrum_ref, str) or isinstance(spectrum_ref, bytes) or isinstance(spectrum_ref, String)), 'arg spectrum_ref wrong type'
    
        cdef size_t _r = self.inst.get().findByReference(deref((convString(spectrum_ref)).get()))
        py_result = <size_t>_r
        return py_result
    
    def addReferenceFormat(self,  regexp ):
        """
        addReferenceFormat(self, regexp: Union[bytes, str, String] ) -> None
        Register a possible format for a spectrum reference
        
        :param regexp: Regular expression defining the format
        """
        assert (isinstance(regexp, str) or isinstance(regexp, bytes) or isinstance(regexp, String)), 'arg regexp wrong type'
    
        self.inst.get().addReferenceFormat(deref((convString(regexp)).get()))
    
    def extractScanNumber(self,  native_id ,  native_id_type_accession ):
        """
        extractScanNumber(self, native_id: Union[bytes, str, String] , native_id_type_accession: Union[bytes, str, String] ) -> int
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
        assert (isinstance(native_id_type_accession, str) or isinstance(native_id_type_accession, bytes) or isinstance(native_id_type_accession, String)), 'arg native_id_type_accession wrong type'
    
    
        cdef int _r = self.inst.get().extractScanNumber(deref((convString(native_id)).get()), deref((convString(native_id_type_accession)).get()))
        py_result = <int>_r
        return py_result 

cdef class SpectrumRangeManager:
    """
    Cython implementation of _SpectrumRangeManager

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumRangeManager.html>`_

    Advanced range manager for MS spectra with separate ranges for each MS level
    
    This class extends the basic RangeManager to provide separate range tracking for different MS levels
    (MS1, MS2, etc.). It manages four types of ranges:
    - m/z (mass-to-charge ratio)
    - intensity
    - retention time (RT)
    - ion mobility
    
    A global range is tracked for all MS levels, and additional ranges are maintained for each specific MS level.
    This allows for efficient querying of ranges for specific MS levels, which is useful for visualization,
    filtering, and processing operations that need to work with specific MS levels.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumRangeManager rv = SpectrumRangeManager.__new__(SpectrumRangeManager)
       rv.inst = shared_ptr[_SpectrumRangeManager](new _SpectrumRangeManager(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumRangeManager rv = SpectrumRangeManager.__new__(SpectrumRangeManager)
       rv.inst = shared_ptr[_SpectrumRangeManager](new _SpectrumRangeManager(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectrumRangeManager](new _SpectrumRangeManager())
    
    def _init_1(self, SpectrumRangeManager in_0 ):
        """
        _init_1(self, in_0: SpectrumRangeManager ) -> None
        """
        assert isinstance(in_0, SpectrumRangeManager), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumRangeManager](new _SpectrumRangeManager((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumRangeManager ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumRangeManager)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def clearRanges(self):
        """
        clearRanges(self) -> None
        """
        self.inst.get().clearRanges()
    
    def getMSLevels(self):
        """
        getMSLevels(self) -> Set[int]
        """
        _r = self.inst.get().getMSLevels()
        cdef set py_result = _r
        return py_result
    
    def extendRT(self, double rt ,  ms_level ):
        """
        extendRT(self, rt: float , ms_level: int ) -> None
        """
        assert isinstance(rt, float), 'arg rt wrong type'
        assert isinstance(ms_level, int), 'arg ms_level wrong type'
    
    
        self.inst.get().extendRT((<double>rt), (<unsigned int>ms_level))
    
    def extendMZ(self, double mz ,  ms_level ):
        """
        extendMZ(self, mz: float , ms_level: int ) -> None
        """
        assert isinstance(mz, float), 'arg mz wrong type'
        assert isinstance(ms_level, int), 'arg ms_level wrong type'
    
    
        self.inst.get().extendMZ((<double>mz), (<unsigned int>ms_level))
    
    def extendUnsafe(self, MSSpectrum spectrum ,  ms_level ):
        """
        extendUnsafe(self, spectrum: MSSpectrum , ms_level: int ) -> None
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(ms_level, int), 'arg ms_level wrong type'
    
    
        self.inst.get().extendUnsafe((deref(spectrum.inst.get())), (<unsigned int>ms_level))
    
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
    
    def getMinMobility(self):
        """
        getMinMobility(self) -> float
        """
        cdef double _r = self.inst.get().getMinMobility()
        py_result = <double>_r
        return py_result
    
    def getMaxMobility(self):
        """
        getMaxMobility(self) -> float
        """
        cdef double _r = self.inst.get().getMaxMobility()
        py_result = <double>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, SpectrumRangeManager):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef SpectrumRangeManager other_casted = other
        cdef SpectrumRangeManager self_casted = self
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
