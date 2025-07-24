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

cdef class __ProcessingAction:
    None
    DATA_PROCESSING = 0
    CHARGE_DECONVOLUTION = 1
    DEISOTOPING = 2
    SMOOTHING = 3
    CHARGE_CALCULATION = 4
    PRECURSOR_RECALCULATION = 5
    BASELINE_REDUCTION = 6
    PEAK_PICKING = 7
    ALIGNMENT = 8
    CALIBRATION = 9
    NORMALIZATION = 10
    FILTERING = 11
    QUANTITATION = 12
    FEATURE_GROUPING = 13
    IDENTIFICATION_MAPPING = 14
    FORMAT_CONVERSION = 15
    CONVERSION_MZDATA = 16
    CONVERSION_MZML = 17
    CONVERSION_MZXML = 18
    CONVERSION_DTA = 19
    IDENTIFICATION = 20
    SIZE_OF_PROCESSINGACTION = 21

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class BiGaussModel:
    """
    Cython implementation of _BiGaussModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BiGaussModel.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef BiGaussModel rv = BiGaussModel.__new__(BiGaussModel)
       rv.inst = shared_ptr[_BiGaussModel](new _BiGaussModel(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef BiGaussModel rv = BiGaussModel.__new__(BiGaussModel)
       rv.inst = shared_ptr[_BiGaussModel](new _BiGaussModel(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BiGaussModel](new _BiGaussModel())
    
    def _init_1(self, BiGaussModel in_0 ):
        """
        _init_1(self, in_0: BiGaussModel ) -> None
        """
        assert isinstance(in_0, BiGaussModel), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BiGaussModel](new _BiGaussModel((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: BiGaussModel ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], BiGaussModel)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setOffset(self, double offset ):
        """
        setOffset(self, offset: float ) -> None
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset((<double>offset))
    
    def setSamples(self):
        """
        setSamples(self) -> None
        """
        self.inst.get().setSamples()
    
    def getCenter(self):
        """
        getCenter(self) -> float
        """
        cdef double _r = self.inst.get().getCenter()
        py_result = <double>_r
        return py_result 

cdef class DataProcessing:
    """
    Cython implementation of _DataProcessing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataProcessing.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DataProcessing rv = DataProcessing.__new__(DataProcessing)
       rv.inst = shared_ptr[_DataProcessing](new _DataProcessing(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DataProcessing rv = DataProcessing.__new__(DataProcessing)
       rv.inst = shared_ptr[_DataProcessing](new _DataProcessing(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DataProcessing](new _DataProcessing())
    
    def _init_1(self, DataProcessing in_0 ):
        """
        _init_1(self, in_0: DataProcessing ) -> None
        """
        assert isinstance(in_0, DataProcessing), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DataProcessing](new _DataProcessing((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DataProcessing ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DataProcessing)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setProcessingActions(self, set in_0 ):
        """
        setProcessingActions(self, in_0: Set[int] ) -> None
        """
        assert isinstance(in_0, set) and all(li in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21] for li in in_0), 'arg in_0 wrong type'
        cdef libcpp_set[_ProcessingAction] * v0 = new libcpp_set[_ProcessingAction]()
        cdef int item0
        for item0 in in_0:
           v0.insert(<_ProcessingAction> item0)
        self.inst.get().setProcessingActions(deref(v0))
        del v0
    
    def getProcessingActions(self):
        """
        getProcessingActions(self) -> Set[int]
        """
        _r = self.inst.get().getProcessingActions()
        py_result = set()
        cdef libcpp_set[_ProcessingAction].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<int>deref(it__r))
           inc(it__r)
        return py_result
    
    def getSoftware(self):
        """
        getSoftware(self) -> Software
        """
        cdef _Software * _r = new _Software(self.inst.get().getSoftware())
        cdef Software py_result = Software.__new__(Software)
        py_result.inst = shared_ptr[_Software](_r)
        return py_result
    
    def setSoftware(self, Software s ):
        """
        setSoftware(self, s: Software ) -> None
        """
        assert isinstance(s, Software), 'arg s wrong type'
    
        self.inst.get().setSoftware((deref(s.inst.get())))
    
    def getCompletionTime(self):
        """
        getCompletionTime(self) -> DateTime
        """
        cdef _DateTime * _r = new _DateTime(self.inst.get().getCompletionTime())
        cdef DateTime py_result = DateTime.__new__(DateTime)
        py_result.inst = shared_ptr[_DateTime](_r)
        return py_result
    
    def setCompletionTime(self, DateTime t ):
        """
        setCompletionTime(self, t: DateTime ) -> None
        """
        assert isinstance(t, DateTime), 'arg t wrong type'
    
        self.inst.get().setCompletionTime((deref(t.inst.get())))
    
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
        if not isinstance(other, DataProcessing):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DataProcessing other_casted = other
        cdef DataProcessing self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    ProcessingAction = __ProcessingAction 

cdef class EmgFitter1D:
    """
    Cython implementation of _EmgFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgFitter1D.html>`_
      -- Inherits from ['LevMarqFitter1D']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EmgFitter1D rv = EmgFitter1D.__new__(EmgFitter1D)
       rv.inst = shared_ptr[_EmgFitter1D](new _EmgFitter1D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EmgFitter1D rv = EmgFitter1D.__new__(EmgFitter1D)
       rv.inst = shared_ptr[_EmgFitter1D](new _EmgFitter1D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization
        """
        self.inst = shared_ptr[_EmgFitter1D](new _EmgFitter1D())
    
    def _init_1(self, EmgFitter1D in_0 ):
        """
        _init_1(self, in_0: EmgFitter1D ) -> None
        """
        assert isinstance(in_0, EmgFitter1D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmgFitter1D](new _EmgFitter1D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EmgFitter1D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EmgFitter1D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class GridBasedCluster:
    """
    Cython implementation of _GridBasedCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GridBasedCluster.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GridBasedCluster rv = GridBasedCluster.__new__(GridBasedCluster)
       rv.inst = shared_ptr[_GridBasedCluster](new _GridBasedCluster(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GridBasedCluster rv = GridBasedCluster.__new__(GridBasedCluster)
       rv.inst = shared_ptr[_GridBasedCluster](new _GridBasedCluster(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  centre , DBoundingBox2 bounding_box , list point_indices ,  property_A , list properties_B ):
        """
        _init_0(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] , property_A: int , properties_B: List[int] ) -> None
        """
        assert len(centre) == 2 and isinstance(centre[0], (int, float)) and isinstance(centre[1], (int, float)), 'arg centre wrong type'
        assert isinstance(bounding_box, DBoundingBox2), 'arg bounding_box wrong type'
        assert isinstance(point_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in point_indices), 'arg point_indices wrong type'
        assert isinstance(property_A, int), 'arg property_A wrong type'
        assert isinstance(properties_B, list) and all(isinstance(elemt_rec, int) for elemt_rec in properties_B), 'arg properties_B wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>centre[0]
        _dp_0[1] = <float>centre[1]
    
        cdef libcpp_vector[int] v2 = point_indices
    
        cdef libcpp_vector[int] v4 = properties_B
        self.inst = shared_ptr[_GridBasedCluster](new _GridBasedCluster(_dp_0, (deref(bounding_box.inst.get())), v2, (<int>property_A), v4))
        
        
    
    def _init_1(self,  centre , DBoundingBox2 bounding_box , list point_indices ):
        """
        _init_1(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] ) -> None
        """
        assert len(centre) == 2 and isinstance(centre[0], (int, float)) and isinstance(centre[1], (int, float)), 'arg centre wrong type'
        assert isinstance(bounding_box, DBoundingBox2), 'arg bounding_box wrong type'
        assert isinstance(point_indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in point_indices), 'arg point_indices wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>centre[0]
        _dp_0[1] = <float>centre[1]
    
        cdef libcpp_vector[int] v2 = point_indices
        self.inst = shared_ptr[_GridBasedCluster](new _GridBasedCluster(_dp_0, (deref(bounding_box.inst.get())), v2))
        
    
    def _init_2(self, GridBasedCluster in_0 ):
        """
        _init_2(self, in_0: GridBasedCluster ) -> None
        """
        assert isinstance(in_0, GridBasedCluster), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GridBasedCluster](new _GridBasedCluster((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] , property_A: int , properties_B: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GridBasedCluster ) -> None
          :noindex:
    
        """
        if (len(args)==5) and (len(args[0]) == 2 and isinstance(args[0][0], (int, float)) and isinstance(args[0][1], (int, float))) and (isinstance(args[1], DBoundingBox2)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[2])) and (isinstance(args[3], int)) and (isinstance(args[4], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[4])):
             self._init_0(*args)
        elif (len(args)==3) and (len(args[0]) == 2 and isinstance(args[0][0], (int, float)) and isinstance(args[0][1], (int, float))) and (isinstance(args[1], DBoundingBox2)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, int) for elemt_rec in args[2])):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], GridBasedCluster)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCentre(self):
        """
        getCentre(self) -> Union[Sequence[int], Sequence[float]]
        Returns cluster centre
        """
        cdef _DPosition2 _r = self.inst.get().getCentre()
        py_result = [_r[0], _r[1]]
        return py_result
    
    def getBoundingBox(self):
        """
        getBoundingBox(self) -> DBoundingBox2
        Returns bounding box
        """
        cdef _DBoundingBox2 * _r = new _DBoundingBox2(self.inst.get().getBoundingBox())
        cdef DBoundingBox2 py_result = DBoundingBox2.__new__(DBoundingBox2)
        py_result.inst = shared_ptr[_DBoundingBox2](_r)
        return py_result
    
    def getPoints(self):
        """
        getPoints(self) -> List[int]
        Returns indices of points in cluster
        """
        _r = self.inst.get().getPoints()
        cdef list py_result = _r
        return py_result
    
    def getPropertyA(self):
        """
        getPropertyA(self) -> int
        Returns property A
        """
        cdef int _r = self.inst.get().getPropertyA()
        py_result = <int>_r
        return py_result
    
    def getPropertiesB(self):
        """
        getPropertiesB(self) -> List[int]
        Returns properties B of all points
        """
        _r = self.inst.get().getPropertiesB()
        cdef list py_result = _r
        return py_result 

cdef class IDConflictResolverAlgorithm:
    """
    Cython implementation of _IDConflictResolverAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDConflictResolverAlgorithm.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IDConflictResolverAlgorithm rv = IDConflictResolverAlgorithm.__new__(IDConflictResolverAlgorithm)
       rv.inst = shared_ptr[_IDConflictResolverAlgorithm](new _IDConflictResolverAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IDConflictResolverAlgorithm rv = IDConflictResolverAlgorithm.__new__(IDConflictResolverAlgorithm)
       rv.inst = shared_ptr[_IDConflictResolverAlgorithm](new _IDConflictResolverAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Resolves ambiguous annotations of features with peptide identifications
        """
        self.inst = shared_ptr[_IDConflictResolverAlgorithm](new _IDConflictResolverAlgorithm())
    
    def _init_1(self, IDConflictResolverAlgorithm in_0 ):
        """
        _init_1(self, in_0: IDConflictResolverAlgorithm ) -> None
        """
        assert isinstance(in_0, IDConflictResolverAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IDConflictResolverAlgorithm](new _IDConflictResolverAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Resolves ambiguous annotations of features with peptide identifications

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IDConflictResolverAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IDConflictResolverAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _resolve_0(self, FeatureMap features ):
        """
        _resolve_0(self, features: FeatureMap ) -> None
        Resolves ambiguous annotations of features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
    
        self.inst.get().resolve((deref(features.inst.get())))
    
    def _resolve_1(self, ConsensusMap features ):
        """
        _resolve_1(self, features: ConsensusMap ) -> None
        Resolves ambiguous annotations of consensus features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
        """
        assert isinstance(features, ConsensusMap), 'arg features wrong type'
    
        self.inst.get().resolve((deref(features.inst.get())))
    
    def resolve(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: resolve(self, features: FeatureMap ) -> None
          :noindex:
        
        Resolves ambiguous annotations of features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
        
        .. rubric:: Overload:
        .. py:function:: resolve(self, features: ConsensusMap ) -> None
          :noindex:
        
        Resolves ambiguous annotations of consensus features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
    
        """
        if (len(args)==1) and (isinstance(args[0], FeatureMap)):
            return self._resolve_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusMap)):
            return self._resolve_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _resolveBetweenFeatures_0(self, FeatureMap features ):
        """
        _resolveBetweenFeatures_0(self, features: FeatureMap ) -> None
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
    
        self.inst.get().resolveBetweenFeatures((deref(features.inst.get())))
    
    def _resolveBetweenFeatures_1(self, ConsensusMap features ):
        """
        _resolveBetweenFeatures_1(self, features: ConsensusMap ) -> None
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
        """
        assert isinstance(features, ConsensusMap), 'arg features wrong type'
    
        self.inst.get().resolveBetweenFeatures((deref(features.inst.get())))
    
    def resolveBetweenFeatures(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: resolveBetweenFeatures(self, features: FeatureMap ) -> None
          :noindex:
        
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
        
        .. rubric:: Overload:
        .. py:function:: resolveBetweenFeatures(self, features: ConsensusMap ) -> None
          :noindex:
        
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
    
        """
        if (len(args)==1) and (isinstance(args[0], FeatureMap)):
            return self._resolveBetweenFeatures_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusMap)):
            return self._resolveBetweenFeatures_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MSExperiment:
    """
    Cython implementation of _MSExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSExperiment.html>`_
      -- Inherits from ['ExperimentalSettings']

    In-Memory representation of a mass spectrometry experiment.
    
    Contains the data and metadata of an experiment performed with an MS (or
    HPLC and MS). This representation of an MS experiment is organized as list
    of spectra and chromatograms and provides an in-memory representation of
    popular mass-spectrometric file formats such as mzXML or mzML. The
    meta-data associated with an experiment is contained in
    ExperimentalSettings (by inheritance) while the raw data (as well as
    spectra and chromatogram level meta data) is stored in objects of type
    MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
    and getChromatogram functions.
    
    Spectra can be accessed by direct iteration or by getSpectrum(),
    while chromatograms are accessed through getChromatogram().
    See help(ExperimentalSettings) for information about meta-data.
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSExperiment rv = MSExperiment.__new__(MSExperiment)
       rv.inst = shared_ptr[_MSExperiment](new _MSExperiment(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSExperiment rv = MSExperiment.__new__(MSExperiment)
       rv.inst = shared_ptr[_MSExperiment](new _MSExperiment(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSExperiment](new _MSExperiment())
    
    def _init_1(self, MSExperiment in_0 ):
        """
        _init_1(self, in_0: MSExperiment ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSExperiment](new _MSExperiment((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSExperiment ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSExperiment)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getExperimentalSettings(self):
        """
        getExperimentalSettings(self) -> ExperimentalSettings
        """
        cdef _ExperimentalSettings * _r = new _ExperimentalSettings(self.inst.get().getExperimentalSettings())
        cdef ExperimentalSettings py_result = ExperimentalSettings.__new__(ExperimentalSettings)
        py_result.inst = shared_ptr[_ExperimentalSettings](_r)
        return py_result
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> MSSpectrum
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _MSSpectrum * _r = new _MSSpectrum(deref(self.inst.get())[(<size_t>in_0)])
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    def __setitem__(self, key, MSSpectrum value):
        """Cython signature: MSSpectrum & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def addSpectrum(self, MSSpectrum spec ):
        """
        addSpectrum(self, spec: MSSpectrum ) -> None
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().addSpectrum((deref(spec.inst.get())))
    
    def setSpectra(self, list spectra ):
        """
        setSpectra(self, spectra: List[MSSpectrum] ) -> None
        """
        assert isinstance(spectra, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in spectra), 'arg spectra wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in spectra:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSpectra(deref(v0))
        cdef libcpp_vector[_MSSpectrum].iterator it_spectra = v0.begin()
        replace_0 = []
        while it_spectra != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_spectra)))
            replace_0.append(item0)
            inc(it_spectra)
        spectra[:] = replace_0
        del v0
    
    def getSpectra(self):
        """
        getSpectra(self) -> List[MSSpectrum]
        """
        _r = self.inst.get().getSpectra()
        py_result = []
        cdef libcpp_vector[_MSSpectrum].iterator it__r = _r.begin()
        cdef MSSpectrum item_py_result
        while it__r != _r.end():
           item_py_result = MSSpectrum.__new__(MSSpectrum)
           item_py_result.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def aggregateFromMatrix(self, MatrixDouble ranges ,  ms_level , bytes mz_agg ):
        """
        aggregateFromMatrix(self, ranges: MatrixDouble , ms_level: int , mz_agg: bytes ) -> List[List[float]]
        Aggregates intensity values for multiple m/z and RT ranges specified in a matrix
        """
        assert isinstance(ranges, MatrixDouble), 'arg ranges wrong type'
        assert isinstance(ms_level, int), 'arg ms_level wrong type'
        assert isinstance(mz_agg, bytes), 'arg mz_agg wrong type'
    
    
    
        _r = self.inst.get().aggregateFromMatrix((deref(ranges.inst.get())), (<unsigned int>ms_level), (<libcpp_string>mz_agg))
        cdef list py_result = _r
        return py_result
    
    def extractXICsFromMatrix(self, MatrixDouble ranges ,  ms_level , bytes mz_agg ):
        """
        extractXICsFromMatrix(self, ranges: MatrixDouble , ms_level: int , mz_agg: bytes ) -> List[MSChromatogram]
        Extracts XIC chromatograms for multiple m/z and RT ranges specified in a matrix
        """
        assert isinstance(ranges, MatrixDouble), 'arg ranges wrong type'
        assert isinstance(ms_level, int), 'arg ms_level wrong type'
        assert isinstance(mz_agg, bytes), 'arg mz_agg wrong type'
    
    
    
        _r = self.inst.get().extractXICsFromMatrix((deref(ranges.inst.get())), (<unsigned int>ms_level), (<libcpp_string>mz_agg))
        py_result = []
        cdef libcpp_vector[_MSChromatogram].iterator it__r = _r.begin()
        cdef MSChromatogram item_py_result
        while it__r != _r.end():
           item_py_result = MSChromatogram.__new__(MSChromatogram)
           item_py_result.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addChromatogram(self, MSChromatogram chromatogram ):
        """
        addChromatogram(self, chromatogram: MSChromatogram ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
    
        self.inst.get().addChromatogram((deref(chromatogram.inst.get())))
    
    def setChromatograms(self, list chromatograms ):
        """
        setChromatograms(self, chromatograms: List[MSChromatogram] ) -> None
        """
        assert isinstance(chromatograms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in chromatograms), 'arg chromatograms wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in chromatograms:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setChromatograms(deref(v0))
        del v0
    
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
    
    def calculateTIC(self):
        """
        calculateTIC(self) -> MSChromatogram
        Returns the total ion chromatogram
        """
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().calculateTIC())
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def clear(self, bool clear_meta_data ):
        """
        clear(self, clear_meta_data: bool ) -> None
        Clear all spectra data and meta data (if called with True)
        """
        assert isinstance(clear_meta_data, pybool_t), 'arg clear_meta_data wrong type'
    
        self.inst.get().clear((<bool>clear_meta_data))
    
    def updateRanges(self):
        """
        updateRanges(self) -> None
        Recalculate global ranges for both spectra and chromatrograms after changes to the data has been made.
        """
        self.inst.get().updateRanges()
    
    def reserveSpaceSpectra(self,  s ):
        """
        reserveSpaceSpectra(self, s: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
    
        self.inst.get().reserveSpaceSpectra((<size_t>s))
    
    def reserveSpaceChromatograms(self,  s ):
        """
        reserveSpaceChromatograms(self, s: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
    
        self.inst.get().reserveSpaceChromatograms((<size_t>s))
    
    def getSize(self):
        """
        getSize(self) -> int
        Returns the total number of peaks
        """
        cdef uint64_t _r = self.inst.get().getSize()
        py_result = <uint64_t>_r
        return py_result
    
    def size(self):
        """
        size(self) -> int
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def resize(self,  s ):
        """
        resize(self, s: int ) -> None
        """
        assert isinstance(s, int) and s >= 0, 'arg s wrong type'
    
        self.inst.get().resize((<size_t>s))
    
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
    
    def getNrSpectra(self):
        """
        getNrSpectra(self) -> int
        Returns the number of MS spectra
        """
        cdef size_t _r = self.inst.get().getNrSpectra()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatograms(self):
        """
        getNrChromatograms(self) -> int
        Returns the number of chromatograms
        """
        cdef size_t _r = self.inst.get().getNrChromatograms()
        py_result = <size_t>_r
        return py_result
    
    def _sortSpectra_0(self, bool sort_mz ):
        """
        _sortSpectra_0(self, sort_mz: bool ) -> None
        Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
        """
        assert isinstance(sort_mz, pybool_t), 'arg sort_mz wrong type'
    
        self.inst.get().sortSpectra((<bool>sort_mz))
    
    def _sortSpectra_1(self):
        """
        _sortSpectra_1(self) -> None
        """
        self.inst.get().sortSpectra()
    
    def sortSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: sortSpectra(self, sort_mz: bool ) -> None
          :noindex:
        
        Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z

        
        .. rubric:: Overload:
        .. py:function:: sortSpectra(self, ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._sortSpectra_0(*args)
        elif not args:
            return self._sortSpectra_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _sortChromatograms_0(self, bool sort_rt ):
        """
        _sortChromatograms_0(self, sort_rt: bool ) -> None
        Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT
        """
        assert isinstance(sort_rt, pybool_t), 'arg sort_rt wrong type'
    
        self.inst.get().sortChromatograms((<bool>sort_rt))
    
    def _sortChromatograms_1(self):
        """
        _sortChromatograms_1(self) -> None
        """
        self.inst.get().sortChromatograms()
    
    def sortChromatograms(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: sortChromatograms(self, sort_rt: bool ) -> None
          :noindex:
        
        Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT

        
        .. rubric:: Overload:
        .. py:function:: sortChromatograms(self, ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._sortChromatograms_0(*args)
        elif not args:
            return self._sortChromatograms_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _isSorted_0(self, bool check_mz ):
        """
        _isSorted_0(self, check_mz: bool ) -> bool
        Checks if all spectra are sorted with respect to ascending RT
        """
        assert isinstance(check_mz, pybool_t), 'arg check_mz wrong type'
    
        cdef bool _r = self.inst.get().isSorted((<bool>check_mz))
        py_result = <bool>_r
        return py_result
    
    def _isSorted_1(self):
        """
        _isSorted_1(self) -> bool
        """
        cdef bool _r = self.inst.get().isSorted()
        py_result = <bool>_r
        return py_result
    
    def isSorted(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: isSorted(self, check_mz: bool ) -> bool
          :noindex:
        
        Checks if all spectra are sorted with respect to ascending RT

        
        .. rubric:: Overload:
        .. py:function:: isSorted(self, ) -> bool
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._isSorted_0(*args)
        elif not args:
            return self._isSorted_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPrimaryMSRunPath(self, list toFill ):
        """
        getPrimaryMSRunPath(self, toFill: List[bytes] ) -> None
        References to the first MS file(s) after conversions. Used to trace results back to original data.
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
    
    def swap(self, MSExperiment in_0 ):
        """
        swap(self, in_0: MSExperiment ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
    
        self.inst.get().swap((deref(in_0.inst.get())))
    
    def reset(self):
        """
        reset(self) -> None
        """
        self.inst.get().reset()
    
    def clearMetaDataArrays(self):
        """
        clearMetaDataArrays(self) -> bool
        """
        cdef bool _r = self.inst.get().clearMetaDataArrays()
        py_result = <bool>_r
        return py_result
    
    def getPrecursorSpectrum(self,  zero_based_index ):
        """
        getPrecursorSpectrum(self, zero_based_index: int ) -> int
        Returns the index of the precursor spectrum for spectrum at index @p zero_based_index
        """
        assert isinstance(zero_based_index, int), 'arg zero_based_index wrong type'
    
        cdef int _r = self.inst.get().getPrecursorSpectrum((<int>zero_based_index))
        py_result = <int>_r
        return py_result
    
    def spectrumRanges(self):
        """
        spectrumRanges(self) -> SpectrumRangeManager
        Returns a reference to the spectrum range manager
        """
        cdef _SpectrumRangeManager * _r = new _SpectrumRangeManager(self.inst.get().spectrumRanges())
        cdef SpectrumRangeManager py_result = SpectrumRangeManager.__new__(SpectrumRangeManager)
        py_result.inst = shared_ptr[_SpectrumRangeManager](_r)
        return py_result
    
    def chromatogramRanges(self):
        """
        chromatogramRanges(self) -> ChromatogramRangeManager
        Returns a reference to the chromatogram range manager
        """
        cdef _ChromatogramRangeManager * _r = new _ChromatogramRangeManager(self.inst.get().chromatogramRanges())
        cdef ChromatogramRangeManager py_result = ChromatogramRangeManager.__new__(ChromatogramRangeManager)
        py_result.inst = shared_ptr[_ChromatogramRangeManager](_r)
        return py_result
    
    def getMinRT(self):
        """
        getMinRT(self) -> float
        Get the minimum RT value from the combined ranges
        """
        cdef double _r = self.inst.get().getMinRT()
        py_result = <double>_r
        return py_result
    
    def getMaxRT(self):
        """
        getMaxRT(self) -> float
        Get the maximum RT value from the combined ranges
        """
        cdef double _r = self.inst.get().getMaxRT()
        py_result = <double>_r
        return py_result
    
    def getMinMZ(self):
        """
        getMinMZ(self) -> float
        Get the minimum m/z value from the combined ranges
        """
        cdef double _r = self.inst.get().getMinMZ()
        py_result = <double>_r
        return py_result
    
    def getMaxMZ(self):
        """
        getMaxMZ(self) -> float
        Get the maximum m/z value from the combined ranges
        """
        cdef double _r = self.inst.get().getMaxMZ()
        py_result = <double>_r
        return py_result
    
    def getMinIntensity(self):
        """
        getMinIntensity(self) -> float
        Get the minimum intensity value from the combined ranges
        """
        cdef double _r = self.inst.get().getMinIntensity()
        py_result = <double>_r
        return py_result
    
    def getMaxIntensity(self):
        """
        getMaxIntensity(self) -> float
        Get the maximum intensity value from the combined ranges
        """
        cdef double _r = self.inst.get().getMaxIntensity()
        py_result = <double>_r
        return py_result
    
    def getMinMobility(self):
        """
        getMinMobility(self) -> float
        Get the minimum mobility value from the combined ranges
        """
        cdef double _r = self.inst.get().getMinMobility()
        py_result = <double>_r
        return py_result
    
    def getMaxMobility(self):
        """
        getMaxMobility(self) -> float
        Get the maximum mobility value from the combined ranges
        """
        cdef double _r = self.inst.get().getMaxMobility()
        py_result = <double>_r
        return py_result
    
    def clearRanges(self):
        """
        clearRanges(self) -> None
        Clear all ranges in all range managers
        """
        self.inst.get().clearRanges()
    
    def getSourceFiles(self):
        """
        getSourceFiles(self) -> List[SourceFile]
        Returns a reference to the source data file
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
    
    def setSourceFiles(self, list source_files ):
        """
        setSourceFiles(self, source_files: List[SourceFile] ) -> None
        Sets the source data file
        """
        assert isinstance(source_files, list) and all(isinstance(elemt_rec, SourceFile) for elemt_rec in source_files), 'arg source_files wrong type'
        cdef libcpp_vector[_SourceFile] * v0 = new libcpp_vector[_SourceFile]()
        cdef SourceFile item0
        for item0 in source_files:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSourceFiles(deref(v0))
        del v0
    
    def getDateTime(self):
        """
        getDateTime(self) -> DateTime
        Returns the date the experiment was performed
        """
        cdef _DateTime * _r = new _DateTime(self.inst.get().getDateTime())
        cdef DateTime py_result = DateTime.__new__(DateTime)
        py_result.inst = shared_ptr[_DateTime](_r)
        return py_result
    
    def setDateTime(self, DateTime date_time ):
        """
        setDateTime(self, date_time: DateTime ) -> None
        Sets the date the experiment was performed
        """
        assert isinstance(date_time, DateTime), 'arg date_time wrong type'
    
        self.inst.get().setDateTime((deref(date_time.inst.get())))
    
    def getSample(self):
        """
        getSample(self) -> Sample
        Returns a reference to the sample description
        """
        cdef _Sample * _r = new _Sample(self.inst.get().getSample())
        cdef Sample py_result = Sample.__new__(Sample)
        py_result.inst = shared_ptr[_Sample](_r)
        return py_result
    
    def setSample(self, Sample sample ):
        """
        setSample(self, sample: Sample ) -> None
        Sets the sample description
        """
        assert isinstance(sample, Sample), 'arg sample wrong type'
    
        self.inst.get().setSample((deref(sample.inst.get())))
    
    def getContacts(self):
        """
        getContacts(self) -> List[ContactPerson]
        Returns a reference to the list of contact persons
        """
        _r = self.inst.get().getContacts()
        py_result = []
        cdef libcpp_vector[_ContactPerson].iterator it__r = _r.begin()
        cdef ContactPerson item_py_result
        while it__r != _r.end():
           item_py_result = ContactPerson.__new__(ContactPerson)
           item_py_result.inst = shared_ptr[_ContactPerson](new _ContactPerson(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setContacts(self, list contacts ):
        """
        setContacts(self, contacts: List[ContactPerson] ) -> None
        Sets the list of contact persons
        """
        assert isinstance(contacts, list) and all(isinstance(elemt_rec, ContactPerson) for elemt_rec in contacts), 'arg contacts wrong type'
        cdef libcpp_vector[_ContactPerson] * v0 = new libcpp_vector[_ContactPerson]()
        cdef ContactPerson item0
        for item0 in contacts:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setContacts(deref(v0))
        del v0
    
    def getInstrument(self):
        """
        getInstrument(self) -> Instrument
        Returns a reference to the MS instrument description
        """
        cdef _Instrument * _r = new _Instrument(self.inst.get().getInstrument())
        cdef Instrument py_result = Instrument.__new__(Instrument)
        py_result.inst = shared_ptr[_Instrument](_r)
        return py_result
    
    def setInstrument(self, Instrument instrument ):
        """
        setInstrument(self, instrument: Instrument ) -> None
        Sets the MS instrument description
        """
        assert isinstance(instrument, Instrument), 'arg instrument wrong type'
    
        self.inst.get().setInstrument((deref(instrument.inst.get())))
    
    def getHPLC(self):
        """
        getHPLC(self) -> HPLC
        Returns a reference to the description of the HPLC run
        """
        cdef _HPLC * _r = new _HPLC(self.inst.get().getHPLC())
        cdef HPLC py_result = HPLC.__new__(HPLC)
        py_result.inst = shared_ptr[_HPLC](_r)
        return py_result
    
    def setHPLC(self, HPLC hplc ):
        """
        setHPLC(self, hplc: HPLC ) -> None
        Sets the description of the HPLC run
        """
        assert isinstance(hplc, HPLC), 'arg hplc wrong type'
    
        self.inst.get().setHPLC((deref(hplc.inst.get())))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the free-text comment
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  comment ):
        """
        setComment(self, comment: Union[bytes, str, String] ) -> None
        Sets the free-text comment
        """
        assert (isinstance(comment, str) or isinstance(comment, bytes) or isinstance(comment, String)), 'arg comment wrong type'
    
        self.inst.get().setComment(deref((convString(comment)).get()))
    
    def getFractionIdentifier(self):
        """
        getFractionIdentifier(self) -> Union[bytes, str, String]
        Returns fraction identifier
        """
        cdef _String _r = self.inst.get().getFractionIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setFractionIdentifier(self,  fraction_identifier ):
        """
        setFractionIdentifier(self, fraction_identifier: Union[bytes, str, String] ) -> None
        Sets the fraction identifier
        """
        assert (isinstance(fraction_identifier, str) or isinstance(fraction_identifier, bytes) or isinstance(fraction_identifier, String)), 'arg fraction_identifier wrong type'
    
        self.inst.get().setFractionIdentifier(deref((convString(fraction_identifier)).get()))
    
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
        if not isinstance(other, MSExperiment):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MSExperiment other_casted = other
        cdef MSExperiment self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef MSSpectrum out
        while it != self.inst.get().end():
            out = MSSpectrum.__new__(MSSpectrum)
            out.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it)))
            yield out
            inc(it)
    
    def get2DPeakData(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakData(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""        

        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        exp_.get2DPeakDataPerSpectrum(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty)



    def get2DPeakDataLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty] get2DPeakDataLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        exp_.get2DPeakData(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap))
    
    def get2DPeakDataIM(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIM(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[libcpp_vector[float]] mz
        cdef libcpp_vector[libcpp_vector[float]] inty
        cdef libcpp_vector[libcpp_vector[float]] ion_mobility
        exp_.get2DPeakDataIMPerSpectrum(min_rt, max_rt, min_mz, max_mz,  ms_level, rt, mz, inty, ion_mobility)

        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)

        cdef np.ndarray all_mz = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_inty = np.empty(rt.size(), dtype=object)
        cdef np.ndarray all_ion = np.empty(rt.size(), dtype=object)
        cdef ArrayWrapperFloat mz_wrap
        cdef ArrayWrapperFloat inty_wrap
        cdef ArrayWrapperFloat ion_mobility_wrap

        cdef unsigned int i
        for i in range(0, mz.size()):
            mz_wrap = ArrayWrapperFloat()
            inty_wrap = ArrayWrapperFloat()
            ion_mobility_wrap = ArrayWrapperFloat()
            mz_wrap.set_data(mz[i])
            inty_wrap.set_data(inty[i])
            ion_mobility_wrap.set_data(ion_mobility[i])
            all_mz[i] = np.frombuffer(mz_wrap)
            all_inty[i] = np.frombuffer(inty_wrap)
            all_ion[i] = np.frombuffer(ion_mobility_wrap)

        return (np.frombuffer(rt_wrap), all_mz, all_inty, all_ion)

    def get2DPeakDataIMLong(MSExperiment self, float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level):
        """Cython signature: tuple[np.array[float] rt, np.array[float] mz, np.array[float] inty, np.array[float] ion_mobility] get2DPeakDataIMLong(float min_rt, float max_rt, float min_mz, float max_mz, unsigned int ms_level)"""
        cdef _MSExperiment * exp_ = self.inst.get()
        cdef libcpp_vector[float] rt
        cdef libcpp_vector[float] mz
        cdef libcpp_vector[float] inty
        cdef libcpp_vector[float] ion_mobility
        exp_.get2DPeakDataIM(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, inty, ion_mobility)
       
        cdef ArrayWrapperFloat rt_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat mz_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat inty_wrap = ArrayWrapperFloat()
        cdef ArrayWrapperFloat ion_mobility_wrap = ArrayWrapperFloat()
        rt_wrap.set_data(rt)
        mz_wrap.set_data(mz)
        inty_wrap.set_data(inty)
        ion_mobility_wrap.set_data(ion_mobility)

        return (np.asarray(rt_wrap), np.asarray(mz_wrap), np.asarray(inty_wrap), np.asarray(ion_mobility_wrap))

    def getMSLevels(self):
        """Cython signature: list[int] getMSLevels()"""
        cdef libcpp_vector[unsigned int] _r = self.inst.get().getMSLevels()
        cdef libcpp_vector[unsigned int].iterator it__r = _r.begin()
        cdef list result = []
        while it__r != _r.end():
            result.append(deref(it__r))
            inc(it__r)
        return result

    def getChromatogram(self, id_):
        """Cython signature: `MSChromatogram getChromatogram(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrChromatograms(), 'Requested chromatogram %s does not exist, there are only %s chromatograms' % (id_, self.getNrChromatograms() )
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>id_)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result

    def getSpectrum(self, id_):
        """Cython signature: `MSSpectrum getSpectrum(size_t id_)`"""
        assert isinstance(id_, int), 'arg id_ wrong type'
        assert id_ < self.getNrSpectra(), 'Requested spectrum %s does not exist, there are only %s spectra' % (id_, self.getNrSpectra() )
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>id_)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result 

cdef class MzTab:
    """
    Cython implementation of _MzTab

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTab.html>`_

    Data model of MzTab files
    
    Please see the official MzTab specification at https://code.google.com/p/mztab/
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzTab rv = MzTab.__new__(MzTab)
       rv.inst = shared_ptr[_MzTab](new _MzTab(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzTab rv = MzTab.__new__(MzTab)
       rv.inst = shared_ptr[_MzTab](new _MzTab(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzTab](new _MzTab())
    
    def _init_1(self, MzTab in_0 ):
        """
        _init_1(self, in_0: MzTab ) -> None
        """
        assert isinstance(in_0, MzTab), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzTab](new _MzTab((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzTab ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzTab)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OSW_ChromExtractParams:
    """
    Cython implementation of _OSW_ChromExtractParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OSW_ChromExtractParams.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property min_upper_edge_dist:
        def __set__(self, double min_upper_edge_dist):
        
            self.inst.get().min_upper_edge_dist = (<double>min_upper_edge_dist)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().min_upper_edge_dist
            py_result = <double>_r
            return py_result
    
    property mz_extraction_window:
        def __set__(self, double mz_extraction_window):
        
            self.inst.get().mz_extraction_window = (<double>mz_extraction_window)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz_extraction_window
            py_result = <double>_r
            return py_result
    
    property ppm:
        def __set__(self, bool ppm):
        
            self.inst.get().ppm = (<bool>ppm)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().ppm
            py_result = <bool>_r
            return py_result
    
    property extraction_function:
        def __set__(self, bytes extraction_function):
        
            self.inst.get().extraction_function = (<libcpp_string>extraction_function)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().extraction_function
            py_result = <libcpp_string>_r
            return py_result
    
    property rt_extraction_window:
        def __set__(self, double rt_extraction_window):
        
            self.inst.get().rt_extraction_window = (<double>rt_extraction_window)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_extraction_window
            py_result = <double>_r
            return py_result
    
    property extra_rt_extract:
        def __set__(self, double extra_rt_extract):
        
            self.inst.get().extra_rt_extract = (<double>extra_rt_extract)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().extra_rt_extract
            py_result = <double>_r
            return py_result
    
    property im_extraction_window:
        def __set__(self, double im_extraction_window):
        
            self.inst.get().im_extraction_window = (<double>im_extraction_window)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().im_extraction_window
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef OSW_ChromExtractParams rv = OSW_ChromExtractParams.__new__(OSW_ChromExtractParams)
       rv.inst = shared_ptr[_OSW_ChromExtractParams](new _OSW_ChromExtractParams(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSW_ChromExtractParams rv = OSW_ChromExtractParams.__new__(OSW_ChromExtractParams)
       rv.inst = shared_ptr[_OSW_ChromExtractParams](new _OSW_ChromExtractParams(deref(self.inst.get())))
       return rv
    
    def __init__(self, OSW_ChromExtractParams in_0 ):
        """
        __init__(self, in_0: OSW_ChromExtractParams ) -> None
        """
        assert isinstance(in_0, OSW_ChromExtractParams), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSW_ChromExtractParams](new _OSW_ChromExtractParams((deref(in_0.inst.get())))) 

cdef class OpenSwathOSWWriter:
    """
    Cython implementation of _OpenSwathOSWWriter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathOSWWriter.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenSwathOSWWriter rv = OpenSwathOSWWriter.__new__(OpenSwathOSWWriter)
       rv.inst = shared_ptr[_OpenSwathOSWWriter](new _OpenSwathOSWWriter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenSwathOSWWriter rv = OpenSwathOSWWriter.__new__(OpenSwathOSWWriter)
       rv.inst = shared_ptr[_OpenSwathOSWWriter](new _OpenSwathOSWWriter(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  output_filename ,  run_id ,  input_filename , bool uis_scores ):
        """
        _init_0(self, output_filename: Union[bytes, str, String] , run_id: int , input_filename: Union[bytes, str, String] , uis_scores: bool ) -> None
        """
        assert (isinstance(output_filename, str) or isinstance(output_filename, bytes) or isinstance(output_filename, String)), 'arg output_filename wrong type'
        assert isinstance(run_id, int) and run_id >= 0, 'arg run_id wrong type'
        assert (isinstance(input_filename, str) or isinstance(input_filename, bytes) or isinstance(input_filename, String)), 'arg input_filename wrong type'
        assert isinstance(uis_scores, pybool_t), 'arg uis_scores wrong type'
    
    
    
    
        self.inst = shared_ptr[_OpenSwathOSWWriter](new _OpenSwathOSWWriter(deref((convString(output_filename)).get()), (<uint64_t>run_id), deref((convString(input_filename)).get()), (<bool>uis_scores)))
    
    def _init_1(self, OpenSwathOSWWriter in_0 ):
        """
        _init_1(self, in_0: OpenSwathOSWWriter ) -> None
        """
        assert isinstance(in_0, OpenSwathOSWWriter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenSwathOSWWriter](new _OpenSwathOSWWriter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, output_filename: Union[bytes, str, String] , run_id: int , input_filename: Union[bytes, str, String] , uis_scores: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenSwathOSWWriter ) -> None
          :noindex:
    
        """
        if (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int) and args[1] >= 0) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], pybool_t)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenSwathOSWWriter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isActive(self):
        """
        isActive(self) -> bool
        """
        cdef bool _r = self.inst.get().isActive()
        py_result = <bool>_r
        return py_result
    
    def writeHeader(self):
        """
        writeHeader(self) -> None
        Initializes file by generating SQLite tables
        """
        self.inst.get().writeHeader()
    
    def prepareLine(self, LightCompound compound , LightTransition tr , FeatureMap output ,  id_ ):
        """
        prepareLine(self, compound: LightCompound , tr: LightTransition , output: FeatureMap , id_: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Prepare a single line (feature) for output
        
        The result can be flushed to disk using writeLines (either line by line or after collecting several lines)
        
        
        :param pep: The compound (peptide/metabolite) used for extraction
        :param transition: The transition used for extraction
        :param output: The feature map containing all features (each feature will generate one entry in the output)
        :param id: The transition group identifier (peptide/metabolite id)
        :return: A String to be written using writeLines
        """
        assert isinstance(compound, LightCompound), 'arg compound wrong type'
        assert isinstance(tr, LightTransition), 'arg tr wrong type'
        assert isinstance(output, FeatureMap), 'arg output wrong type'
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
    
    
    
    
        cdef _String _r = self.inst.get().prepareLine((deref(compound.inst.get())), (tr.inst.get()), (deref(output.inst.get())), deref((convString(id_)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def writeLines(self, list to_osw_output ):
        """
        writeLines(self, to_osw_output: List[bytes] ) -> None
        Write data to disk
        
        Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
        
        
        :param to_osw_output: Statements generated by prepareLine
        """
        assert isinstance(to_osw_output, list) and all(isinstance(i, bytes) for i in to_osw_output), 'arg to_osw_output wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in to_osw_output:
           v0.push_back(_String(<char *>item0))
        self.inst.get().writeLines(deref(v0))
        del v0 

cdef class PeakTypeEstimator:
    """
    Cython implementation of _PeakTypeEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakTypeEstimator.html>`_

    Estimates if the data of a spectrum is raw data or peak data
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakTypeEstimator rv = PeakTypeEstimator.__new__(PeakTypeEstimator)
       rv.inst = shared_ptr[_PeakTypeEstimator](new _PeakTypeEstimator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakTypeEstimator rv = PeakTypeEstimator.__new__(PeakTypeEstimator)
       rv.inst = shared_ptr[_PeakTypeEstimator](new _PeakTypeEstimator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakTypeEstimator](new _PeakTypeEstimator())
    
    def _init_1(self, PeakTypeEstimator in_0 ):
        """
        _init_1(self, in_0: PeakTypeEstimator ) -> None
        """
        assert isinstance(in_0, PeakTypeEstimator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakTypeEstimator](new _PeakTypeEstimator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakTypeEstimator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakTypeEstimator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def estimateType(self, MSSpectrum spec):
        return self.inst.get().estimateType(spec.inst.get().begin(), spec.inst.get().end()) 

cdef class SequestInfile:
    """
    Cython implementation of _SequestInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SequestInfile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SequestInfile rv = SequestInfile.__new__(SequestInfile)
       rv.inst = shared_ptr[_SequestInfile](new _SequestInfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SequestInfile rv = SequestInfile.__new__(SequestInfile)
       rv.inst = shared_ptr[_SequestInfile](new _SequestInfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Sequest input file adapter
        """
        self.inst = shared_ptr[_SequestInfile](new _SequestInfile())
    
    def _init_1(self, SequestInfile in_0 ):
        """
        _init_1(self, in_0: SequestInfile ) -> None
        """
        assert isinstance(in_0, SequestInfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SequestInfile](new _SequestInfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Sequest input file adapter

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SequestInfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SequestInfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename ):
        """
        store(self, filename: Union[bytes, str, String] ) -> None
        Stores the experiment data in a Sequest input file that can be used as input for Sequest shell execution
        
        :param filename: the name of the file in which the infile is stored into
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().store(deref((convString(filename)).get()))
    
    def getEnzymeInfoAsString(self):
        """
        getEnzymeInfoAsString(self) -> Union[bytes, str, String]
        Returns the enzyme list as a string
        """
        cdef _String _r = self.inst.get().getEnzymeInfoAsString()
        py_result = convOutputString(_r)
        return py_result
    
    def getDatabase(self):
        """
        getDatabase(self) -> Union[bytes, str, String]
        Returns the used database
        """
        cdef _String _r = self.inst.get().getDatabase()
        py_result = convOutputString(_r)
        return py_result
    
    def setDatabase(self,  database ):
        """
        setDatabase(self, database: Union[bytes, str, String] ) -> None
        Sets the used database
        """
        assert (isinstance(database, str) or isinstance(database, bytes) or isinstance(database, String)), 'arg database wrong type'
    
        self.inst.get().setDatabase(deref((convString(database)).get()))
    
    def getNeutralLossesForIons(self):
        """
        getNeutralLossesForIons(self) -> Union[bytes, str, String]
        Returns whether neutral losses are considered for the a-, b- and y-ions
        """
        cdef _String _r = self.inst.get().getNeutralLossesForIons()
        py_result = convOutputString(_r)
        return py_result
    
    def setNeutralLossesForIons(self,  neutral_losses_for_ions ):
        """
        setNeutralLossesForIons(self, neutral_losses_for_ions: Union[bytes, str, String] ) -> None
        Sets whether neutral losses are considered for the a-, b- and y-ions
        """
        assert (isinstance(neutral_losses_for_ions, str) or isinstance(neutral_losses_for_ions, bytes) or isinstance(neutral_losses_for_ions, String)), 'arg neutral_losses_for_ions wrong type'
    
        self.inst.get().setNeutralLossesForIons(deref((convString(neutral_losses_for_ions)).get()))
    
    def getIonSeriesWeights(self):
        """
        getIonSeriesWeights(self) -> Union[bytes, str, String]
        Returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        """
        cdef _String _r = self.inst.get().getIonSeriesWeights()
        py_result = convOutputString(_r)
        return py_result
    
    def setIonSeriesWeights(self,  ion_series_weights ):
        """
        setIonSeriesWeights(self, ion_series_weights: Union[bytes, str, String] ) -> None
        Sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        """
        assert (isinstance(ion_series_weights, str) or isinstance(ion_series_weights, bytes) or isinstance(ion_series_weights, String)), 'arg ion_series_weights wrong type'
    
        self.inst.get().setIonSeriesWeights(deref((convString(ion_series_weights)).get()))
    
    def getPartialSequence(self):
        """
        getPartialSequence(self) -> Union[bytes, str, String]
        Returns the partial sequences (space delimited) that have to occur in the theoretical spectra
        """
        cdef _String _r = self.inst.get().getPartialSequence()
        py_result = convOutputString(_r)
        return py_result
    
    def setPartialSequence(self,  partial_sequence ):
        """
        setPartialSequence(self, partial_sequence: Union[bytes, str, String] ) -> None
        Sets the partial sequences (space delimited) that have to occur in the theoretical spectra
        """
        assert (isinstance(partial_sequence, str) or isinstance(partial_sequence, bytes) or isinstance(partial_sequence, String)), 'arg partial_sequence wrong type'
    
        self.inst.get().setPartialSequence(deref((convString(partial_sequence)).get()))
    
    def getSequenceHeaderFilter(self):
        """
        getSequenceHeaderFilter(self) -> Union[bytes, str, String]
        Returns the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        """
        cdef _String _r = self.inst.get().getSequenceHeaderFilter()
        py_result = convOutputString(_r)
        return py_result
    
    def setSequenceHeaderFilter(self,  sequence_header_filter ):
        """
        setSequenceHeaderFilter(self, sequence_header_filter: Union[bytes, str, String] ) -> None
        Sets the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        """
        assert (isinstance(sequence_header_filter, str) or isinstance(sequence_header_filter, bytes) or isinstance(sequence_header_filter, String)), 'arg sequence_header_filter wrong type'
    
        self.inst.get().setSequenceHeaderFilter(deref((convString(sequence_header_filter)).get()))
    
    def getProteinMassFilter(self):
        """
        getProteinMassFilter(self) -> Union[bytes, str, String]
        Returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        """
        cdef _String _r = self.inst.get().getProteinMassFilter()
        py_result = convOutputString(_r)
        return py_result
    
    def setProteinMassFilter(self,  protein_mass_filter ):
        """
        setProteinMassFilter(self, protein_mass_filter: Union[bytes, str, String] ) -> None
        Sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        """
        assert (isinstance(protein_mass_filter, str) or isinstance(protein_mass_filter, bytes) or isinstance(protein_mass_filter, String)), 'arg protein_mass_filter wrong type'
    
        self.inst.get().setProteinMassFilter(deref((convString(protein_mass_filter)).get()))
    
    def getPeakMassTolerance(self):
        """
        getPeakMassTolerance(self) -> float
        Returns the peak mass tolerance
        """
        cdef float _r = self.inst.get().getPeakMassTolerance()
        py_result = <float>_r
        return py_result
    
    def setPeakMassTolerance(self, float peak_mass_tolerance ):
        """
        setPeakMassTolerance(self, peak_mass_tolerance: float ) -> None
        Sets the peak mass tolerance
        """
        assert isinstance(peak_mass_tolerance, float), 'arg peak_mass_tolerance wrong type'
    
        self.inst.get().setPeakMassTolerance((<float>peak_mass_tolerance))
    
    def getPrecursorMassTolerance(self):
        """
        getPrecursorMassTolerance(self) -> float
        Returns the precursor mass tolerance
        """
        cdef float _r = self.inst.get().getPrecursorMassTolerance()
        py_result = <float>_r
        return py_result
    
    def setPrecursorMassTolerance(self, float precursor_mass_tolerance ):
        """
        setPrecursorMassTolerance(self, precursor_mass_tolerance: float ) -> None
        Sets the precursor mass tolerance
        """
        assert isinstance(precursor_mass_tolerance, float), 'arg precursor_mass_tolerance wrong type'
    
        self.inst.get().setPrecursorMassTolerance((<float>precursor_mass_tolerance))
    
    def getMatchPeakTolerance(self):
        """
        getMatchPeakTolerance(self) -> float
        Returns the match peak tolerance
        """
        cdef float _r = self.inst.get().getMatchPeakTolerance()
        py_result = <float>_r
        return py_result
    
    def setMatchPeakTolerance(self, float match_peak_tolerance ):
        """
        setMatchPeakTolerance(self, match_peak_tolerance: float ) -> None
        Sets the match peak tolerance
        """
        assert isinstance(match_peak_tolerance, float), 'arg match_peak_tolerance wrong type'
    
        self.inst.get().setMatchPeakTolerance((<float>match_peak_tolerance))
    
    def getIonCutoffPercentage(self):
        """
        getIonCutoffPercentage(self) -> float
        Returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks
        """
        cdef float _r = self.inst.get().getIonCutoffPercentage()
        py_result = <float>_r
        return py_result
    
    def setIonCutoffPercentage(self, float ion_cutoff_percentage ):
        """
        setIonCutoffPercentage(self, ion_cutoff_percentage: float ) -> None
        Sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks
        """
        assert isinstance(ion_cutoff_percentage, float), 'arg ion_cutoff_percentage wrong type'
    
        self.inst.get().setIonCutoffPercentage((<float>ion_cutoff_percentage))
    
    def getPeptideMassUnit(self):
        """
        getPeptideMassUnit(self) -> int
        Returns the peptide mass unit
        """
        cdef size_t _r = self.inst.get().getPeptideMassUnit()
        py_result = <size_t>_r
        return py_result
    
    def setPeptideMassUnit(self,  peptide_mass_unit ):
        """
        setPeptideMassUnit(self, peptide_mass_unit: int ) -> None
        Sets the peptide mass unit
        """
        assert isinstance(peptide_mass_unit, int) and peptide_mass_unit >= 0, 'arg peptide_mass_unit wrong type'
    
        self.inst.get().setPeptideMassUnit((<size_t>peptide_mass_unit))
    
    def getOutputLines(self):
        """
        getOutputLines(self) -> int
        Returns the number of peptides to be displayed
        """
        cdef size_t _r = self.inst.get().getOutputLines()
        py_result = <size_t>_r
        return py_result
    
    def setOutputLines(self,  output_lines ):
        """
        setOutputLines(self, output_lines: int ) -> None
        Sets the number of peptides to be displayed
        """
        assert isinstance(output_lines, int) and output_lines >= 0, 'arg output_lines wrong type'
    
        self.inst.get().setOutputLines((<size_t>output_lines))
    
    def getEnzymeNumber(self):
        """
        getEnzymeNumber(self) -> int
        Returns the enzyme used for cleavage (by means of the number from a list of enzymes)
        """
        cdef size_t _r = self.inst.get().getEnzymeNumber()
        py_result = <size_t>_r
        return py_result
    
    def getEnzymeName(self):
        """
        getEnzymeName(self) -> Union[bytes, str, String]
        Returns the enzyme used for cleavage
        """
        cdef _String _r = self.inst.get().getEnzymeName()
        py_result = convOutputString(_r)
        return py_result
    
    def setEnzyme(self,  enzyme_name ):
        """
        setEnzyme(self, enzyme_name: Union[bytes, str, String] ) -> int
        Sets the enzyme used for cleavage (by means of the number from a list of enzymes)
        """
        assert (isinstance(enzyme_name, str) or isinstance(enzyme_name, bytes) or isinstance(enzyme_name, String)), 'arg enzyme_name wrong type'
    
        cdef size_t _r = self.inst.get().setEnzyme(deref((convString(enzyme_name)).get()))
        py_result = <size_t>_r
        return py_result
    
    def getMaxAAPerModPerPeptide(self):
        """
        getMaxAAPerModPerPeptide(self) -> int
        Returns the maximum number of amino acids containing the same modification in a peptide
        """
        cdef size_t _r = self.inst.get().getMaxAAPerModPerPeptide()
        py_result = <size_t>_r
        return py_result
    
    def setMaxAAPerModPerPeptide(self,  max_aa_per_mod_per_peptide ):
        """
        setMaxAAPerModPerPeptide(self, max_aa_per_mod_per_peptide: int ) -> None
        Sets the maximum number of amino acids containing the same modification in a peptide
        """
        assert isinstance(max_aa_per_mod_per_peptide, int) and max_aa_per_mod_per_peptide >= 0, 'arg max_aa_per_mod_per_peptide wrong type'
    
        self.inst.get().setMaxAAPerModPerPeptide((<size_t>max_aa_per_mod_per_peptide))
    
    def getMaxModsPerPeptide(self):
        """
        getMaxModsPerPeptide(self) -> int
        Returns the maximum number of modifications that are allowed in a peptide
        """
        cdef size_t _r = self.inst.get().getMaxModsPerPeptide()
        py_result = <size_t>_r
        return py_result
    
    def setMaxModsPerPeptide(self,  max_mods_per_peptide ):
        """
        setMaxModsPerPeptide(self, max_mods_per_peptide: int ) -> None
        Sets the maximum number of modifications that are allowed in a peptide
        """
        assert isinstance(max_mods_per_peptide, int) and max_mods_per_peptide >= 0, 'arg max_mods_per_peptide wrong type'
    
        self.inst.get().setMaxModsPerPeptide((<size_t>max_mods_per_peptide))
    
    def getNucleotideReadingFrame(self):
        """
        getNucleotideReadingFrame(self) -> int
        Returns the nucleotide reading frame
        """
        cdef size_t _r = self.inst.get().getNucleotideReadingFrame()
        py_result = <size_t>_r
        return py_result
    
    def setNucleotideReadingFrame(self,  nucleotide_reading_frame ):
        """
        setNucleotideReadingFrame(self, nucleotide_reading_frame: int ) -> None
        Sets the nucleotide reading frame
        """
        assert isinstance(nucleotide_reading_frame, int) and nucleotide_reading_frame >= 0, 'arg nucleotide_reading_frame wrong type'
    
        self.inst.get().setNucleotideReadingFrame((<size_t>nucleotide_reading_frame))
    
    def getMaxInternalCleavageSites(self):
        """
        getMaxInternalCleavageSites(self) -> int
        Returns the maximum number of internal cleavage sites
        """
        cdef size_t _r = self.inst.get().getMaxInternalCleavageSites()
        py_result = <size_t>_r
        return py_result
    
    def setMaxInternalCleavageSites(self,  max_internal_cleavage_sites ):
        """
        setMaxInternalCleavageSites(self, max_internal_cleavage_sites: int ) -> None
        Sets the maximum number of internal cleavage sites
        """
        assert isinstance(max_internal_cleavage_sites, int) and max_internal_cleavage_sites >= 0, 'arg max_internal_cleavage_sites wrong type'
    
        self.inst.get().setMaxInternalCleavageSites((<size_t>max_internal_cleavage_sites))
    
    def getMatchPeakCount(self):
        """
        getMatchPeakCount(self) -> int
        Returns the number of top abundant peaks to match with theoretical ones
        """
        cdef size_t _r = self.inst.get().getMatchPeakCount()
        py_result = <size_t>_r
        return py_result
    
    def setMatchPeakCount(self,  match_peak_count ):
        """
        setMatchPeakCount(self, match_peak_count: int ) -> None
        Sets the number of top abundant peaks to with theoretical ones
        """
        assert isinstance(match_peak_count, int) and match_peak_count >= 0, 'arg match_peak_count wrong type'
    
        self.inst.get().setMatchPeakCount((<size_t>match_peak_count))
    
    def getMatchPeakAllowedError(self):
        """
        getMatchPeakAllowedError(self) -> int
        Returns the number of top abundant peaks that are allowed not to match with a theoretical peak
        """
        cdef size_t _r = self.inst.get().getMatchPeakAllowedError()
        py_result = <size_t>_r
        return py_result
    
    def setMatchPeakAllowedError(self,  match_peak_allowed_error ):
        """
        setMatchPeakAllowedError(self, match_peak_allowed_error: int ) -> None
        Sets the number of top abundant peaks that are allowed not to match with a theoretical peak
        """
        assert isinstance(match_peak_allowed_error, int) and match_peak_allowed_error >= 0, 'arg match_peak_allowed_error wrong type'
    
        self.inst.get().setMatchPeakAllowedError((<size_t>match_peak_allowed_error))
    
    def getShowFragmentIons(self):
        """
        getShowFragmentIons(self) -> bool
        Returns whether fragment ions shall be displayed
        """
        cdef bool _r = self.inst.get().getShowFragmentIons()
        py_result = <bool>_r
        return py_result
    
    def setShowFragmentIons(self, bool show_fragments ):
        """
        setShowFragmentIons(self, show_fragments: bool ) -> None
        Sets whether fragment ions shall be displayed
        """
        assert isinstance(show_fragments, pybool_t), 'arg show_fragments wrong type'
    
        self.inst.get().setShowFragmentIons((<bool>show_fragments))
    
    def getPrintDuplicateReferences(self):
        """
        getPrintDuplicateReferences(self) -> bool
        Returns whether all proteins containing a found peptide should be displayed
        """
        cdef bool _r = self.inst.get().getPrintDuplicateReferences()
        py_result = <bool>_r
        return py_result
    
    def setPrintDuplicateReferences(self, bool print_duplicate_references ):
        """
        setPrintDuplicateReferences(self, print_duplicate_references: bool ) -> None
        Sets whether all proteins containing a found peptide should be displayed
        """
        assert isinstance(print_duplicate_references, pybool_t), 'arg print_duplicate_references wrong type'
    
        self.inst.get().setPrintDuplicateReferences((<bool>print_duplicate_references))
    
    def getRemovePrecursorNearPeaks(self):
        """
        getRemovePrecursorNearPeaks(self) -> bool
        Returns whether peaks near (15 amu) the precursor peak are removed
        """
        cdef bool _r = self.inst.get().getRemovePrecursorNearPeaks()
        py_result = <bool>_r
        return py_result
    
    def setRemovePrecursorNearPeaks(self, bool remove_precursor_near_peaks ):
        """
        setRemovePrecursorNearPeaks(self, remove_precursor_near_peaks: bool ) -> None
        Sets whether peaks near (15 amu) the precursor peak are removed
        """
        assert isinstance(remove_precursor_near_peaks, pybool_t), 'arg remove_precursor_near_peaks wrong type'
    
        self.inst.get().setRemovePrecursorNearPeaks((<bool>remove_precursor_near_peaks))
    
    def getMassTypeParent(self):
        """
        getMassTypeParent(self) -> bool
        Returns the mass type of the parent (0 - monoisotopic, 1 - average mass)
        """
        cdef bool _r = self.inst.get().getMassTypeParent()
        py_result = <bool>_r
        return py_result
    
    def setMassTypeParent(self, bool mass_type_parent ):
        """
        setMassTypeParent(self, mass_type_parent: bool ) -> None
        Sets the mass type of the parent (0 - monoisotopic, 1 - average mass)
        """
        assert isinstance(mass_type_parent, pybool_t), 'arg mass_type_parent wrong type'
    
        self.inst.get().setMassTypeParent((<bool>mass_type_parent))
    
    def getMassTypeFragment(self):
        """
        getMassTypeFragment(self) -> bool
        Returns the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        """
        cdef bool _r = self.inst.get().getMassTypeFragment()
        py_result = <bool>_r
        return py_result
    
    def setMassTypeFragment(self, bool mass_type_fragment ):
        """
        setMassTypeFragment(self, mass_type_fragment: bool ) -> None
        Sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        """
        assert isinstance(mass_type_fragment, pybool_t), 'arg mass_type_fragment wrong type'
    
        self.inst.get().setMassTypeFragment((<bool>mass_type_fragment))
    
    def getNormalizeXcorr(self):
        """
        getNormalizeXcorr(self) -> bool
        Returns whether normalized xcorr values are displayed
        """
        cdef bool _r = self.inst.get().getNormalizeXcorr()
        py_result = <bool>_r
        return py_result
    
    def setNormalizeXcorr(self, bool normalize_xcorr ):
        """
        setNormalizeXcorr(self, normalize_xcorr: bool ) -> None
        Sets whether normalized xcorr values are displayed
        """
        assert isinstance(normalize_xcorr, pybool_t), 'arg normalize_xcorr wrong type'
    
        self.inst.get().setNormalizeXcorr((<bool>normalize_xcorr))
    
    def getResiduesInUpperCase(self):
        """
        getResiduesInUpperCase(self) -> bool
        Returns whether residues are in upper case
        """
        cdef bool _r = self.inst.get().getResiduesInUpperCase()
        py_result = <bool>_r
        return py_result
    
    def setResiduesInUpperCase(self, bool residues_in_upper_case ):
        """
        setResiduesInUpperCase(self, residues_in_upper_case: bool ) -> None
        Sets whether residues are in upper case
        """
        assert isinstance(residues_in_upper_case, pybool_t), 'arg residues_in_upper_case wrong type'
    
        self.inst.get().setResiduesInUpperCase((<bool>residues_in_upper_case))
    
    def addEnzymeInfo(self, list enzyme_info ):
        """
        addEnzymeInfo(self, enzyme_info: List[bytes] ) -> None
        Adds an enzyme to the list and sets is as used
        """
        assert isinstance(enzyme_info, list) and all(isinstance(i, bytes) for i in enzyme_info), 'arg enzyme_info wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in enzyme_info:
           v0.push_back(_String(<char *>item0))
        self.inst.get().addEnzymeInfo(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        enzyme_info[:] = replace
        del v0
    
    def handlePTMs(self,  modification_line ,  modifications_filename , bool monoisotopic ):
        """
        handlePTMs(self, modification_line: Union[bytes, str, String] , modifications_filename: Union[bytes, str, String] , monoisotopic: bool ) -> None
        """
        assert (isinstance(modification_line, str) or isinstance(modification_line, bytes) or isinstance(modification_line, String)), 'arg modification_line wrong type'
        assert (isinstance(modifications_filename, str) or isinstance(modifications_filename, bytes) or isinstance(modifications_filename, String)), 'arg modifications_filename wrong type'
        assert isinstance(monoisotopic, pybool_t), 'arg monoisotopic wrong type'
    
    
    
        self.inst.get().handlePTMs(deref((convString(modification_line)).get()), deref((convString(modifications_filename)).get()), (<bool>monoisotopic))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, SequestInfile):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef SequestInfile other_casted = other
        cdef SequestInfile self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    
    def getModifications(self):
        cdef libcpp_map[_String, libcpp_vector[_String] ] c_res 
        c_res = self.inst.get().getModifications()
        pass 
        # 
        ## Get the data back from C++
        #
        replace = dict()
        # Make sure each C++ and each Cython objects you use are declared here
        # (only pure Python objects need no declaration).
        cdef libcpp_map[_String, libcpp_vector[_String] ].iterator it_outer = c_res.begin()
        cdef libcpp_vector[_String] c_inner_vec
        cdef libcpp_vector[_String].iterator it_inner
        cdef _String myString
        cdef String py_String
        while it_outer != c_res.end():
            
            inner_vec = []
            c_inner_vec = deref(it_outer).second
            it_inner = c_inner_vec.begin()
            while it_inner != c_inner_vec.end():

                py_String = String.__new__(String)   
                py_String.inst = shared_ptr[_String](new _String(deref(it_inner)))
                inner_vec.append(py_String)
                inc(it_inner)

            replace[ <libcpp_string>deref(it_outer).first ] = inner_vec
            inc(it_outer)

        return replace 

cdef class TMTTenPlexQuantitationMethod:
    """
    Cython implementation of _TMTTenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTTenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TMTTenPlexQuantitationMethod rv = TMTTenPlexQuantitationMethod.__new__(TMTTenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTTenPlexQuantitationMethod](new _TMTTenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TMTTenPlexQuantitationMethod rv = TMTTenPlexQuantitationMethod.__new__(TMTTenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTTenPlexQuantitationMethod](new _TMTTenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TMTTenPlexQuantitationMethod](new _TMTTenPlexQuantitationMethod())
    
    def _init_1(self, TMTTenPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: TMTTenPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, TMTTenPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TMTTenPlexQuantitationMethod](new _TMTTenPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TMTTenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTTenPlexQuantitationMethod)):
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

cdef class TheoreticalSpectrumGeneratorXLMS:
    """
    Cython implementation of _TheoreticalSpectrumGeneratorXLMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TheoreticalSpectrumGeneratorXLMS.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TheoreticalSpectrumGeneratorXLMS rv = TheoreticalSpectrumGeneratorXLMS.__new__(TheoreticalSpectrumGeneratorXLMS)
       rv.inst = shared_ptr[_TheoreticalSpectrumGeneratorXLMS](new _TheoreticalSpectrumGeneratorXLMS(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TheoreticalSpectrumGeneratorXLMS rv = TheoreticalSpectrumGeneratorXLMS.__new__(TheoreticalSpectrumGeneratorXLMS)
       rv.inst = shared_ptr[_TheoreticalSpectrumGeneratorXLMS](new _TheoreticalSpectrumGeneratorXLMS(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TheoreticalSpectrumGeneratorXLMS](new _TheoreticalSpectrumGeneratorXLMS())
    
    def _init_1(self, TheoreticalSpectrumGeneratorXLMS in_0 ):
        """
        _init_1(self, in_0: TheoreticalSpectrumGeneratorXLMS ) -> None
        """
        assert isinstance(in_0, TheoreticalSpectrumGeneratorXLMS), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TheoreticalSpectrumGeneratorXLMS](new _TheoreticalSpectrumGeneratorXLMS((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TheoreticalSpectrumGeneratorXLMS ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TheoreticalSpectrumGeneratorXLMS)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getLinearIonSpectrum(self, MSSpectrum spectrum , AASequence peptide ,  link_pos , bool frag_alpha ,  charge ,  link_pos_2 ):
        """
        getLinearIonSpectrum(self, spectrum: MSSpectrum , peptide: AASequence , link_pos: int , frag_alpha: bool , charge: int , link_pos_2: int ) -> None
            Generates fragment ions not containing the cross-linker for one peptide
        
            B-ions are generated from the beginning of the peptide up to the first linked position,
            y-ions are generated from the second linked position up the end of the peptide.
            If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
            For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param peptide: The peptide to fragment
            :param link_pos: The position of the cross-linker on the given peptide
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation
            :param charge: The maximal charge of the ions
            :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(link_pos, int) and link_pos >= 0, 'arg link_pos wrong type'
        assert isinstance(frag_alpha, pybool_t), 'arg frag_alpha wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(link_pos_2, int) and link_pos_2 >= 0, 'arg link_pos_2 wrong type'
    
    
    
    
    
    
        self.inst.get().getLinearIonSpectrum((deref(spectrum.inst.get())), (deref(peptide.inst.get())), (<size_t>link_pos), (<bool>frag_alpha), (<int>charge), (<size_t>link_pos_2))
    
    def _getXLinkIonSpectrum_0(self, MSSpectrum spectrum , AASequence peptide ,  link_pos , double precursor_mass , bool frag_alpha ,  mincharge ,  maxcharge ,  link_pos_2 ):
        """
        _getXLinkIonSpectrum_0(self, spectrum: MSSpectrum , peptide: AASequence , link_pos: int , precursor_mass: float , frag_alpha: bool , mincharge: int , maxcharge: int , link_pos_2: int ) -> None
            Generates fragment ions containing the cross-linker for one peptide
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
            For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
            Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
            a precursor mass for the two peptides and the linker is needed.
            In the case of a loop link the precursor mass is the mass of the only peptide and the linker.
            Although this function is more general, currently it is mainly used for loop-links and mono-links,
            because residues in the second, unknown peptide cannot be considered for possible neutral losses.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param peptide: The peptide to fragment
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
            :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(link_pos, int) and link_pos >= 0, 'arg link_pos wrong type'
        assert isinstance(precursor_mass, float), 'arg precursor_mass wrong type'
        assert isinstance(frag_alpha, pybool_t), 'arg frag_alpha wrong type'
        assert isinstance(mincharge, int), 'arg mincharge wrong type'
        assert isinstance(maxcharge, int), 'arg maxcharge wrong type'
        assert isinstance(link_pos_2, int) and link_pos_2 >= 0, 'arg link_pos_2 wrong type'
    
    
    
    
    
    
    
    
        self.inst.get().getXLinkIonSpectrum((deref(spectrum.inst.get())), (deref(peptide.inst.get())), (<size_t>link_pos), (<double>precursor_mass), (<bool>frag_alpha), (<int>mincharge), (<int>maxcharge), (<size_t>link_pos_2))
    
    def _getXLinkIonSpectrum_1(self, MSSpectrum spectrum , ProteinProteinCrossLink crosslink , bool frag_alpha ,  mincharge ,  maxcharge ):
        """
        _getXLinkIonSpectrum_1(self, spectrum: MSSpectrum , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None
            Generates fragment ions containing the cross-linker for a pair of peptides
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            This function generates neutral loss ions by considering both linked peptides.
            Only one of the peptides, decided by @frag_alpha, is fragmented.
            This function is not suitable to generate fragments for mono-links or loop-links.
            This simplifies the function, but it has to be called twice to get all fragments of a peptide pair.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param crosslink: ProteinProteinCrossLink to be fragmented
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(crosslink, ProteinProteinCrossLink), 'arg crosslink wrong type'
        assert isinstance(frag_alpha, pybool_t), 'arg frag_alpha wrong type'
        assert isinstance(mincharge, int), 'arg mincharge wrong type'
        assert isinstance(maxcharge, int), 'arg maxcharge wrong type'
    
    
    
    
    
        self.inst.get().getXLinkIonSpectrum((deref(spectrum.inst.get())), (deref(crosslink.inst.get())), (<bool>frag_alpha), (<int>mincharge), (<int>maxcharge))
    
    def getXLinkIonSpectrum(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getXLinkIonSpectrum(self, spectrum: MSSpectrum , peptide: AASequence , link_pos: int , precursor_mass: float , frag_alpha: bool , mincharge: int , maxcharge: int , link_pos_2: int ) -> None
          :noindex:
        
            Generates fragment ions containing the cross-linker for one peptide
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
            For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
            Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
            a precursor mass for the two peptides and the linker is needed.
            In the case of a loop link the precursor mass is the mass of the only peptide and the linker.
            Although this function is more general, currently it is mainly used for loop-links and mono-links,
            because residues in the second, unknown peptide cannot be considered for possible neutral losses.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param peptide: The peptide to fragment
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
            :param link_pos_2: A second position for the linker, in case it is a loop link
        
        .. rubric:: Overload:
        .. py:function:: getXLinkIonSpectrum(self, spectrum: MSSpectrum , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None
          :noindex:
        
            Generates fragment ions containing the cross-linker for a pair of peptides
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            This function generates neutral loss ions by considering both linked peptides.
            Only one of the peptides, decided by @frag_alpha, is fragmented.
            This function is not suitable to generate fragments for mono-links or loop-links.
            This simplifies the function, but it has to be called twice to get all fragments of a peptide pair.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param crosslink: ProteinProteinCrossLink to be fragmented
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
    
        """
        if (len(args)==8) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], AASequence)) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], float)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], int)) and (isinstance(args[6], int)) and (isinstance(args[7], int) and args[7] >= 0):
            return self._getXLinkIonSpectrum_0(*args)
        elif (len(args)==5) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], ProteinProteinCrossLink)) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], int)) and (isinstance(args[4], int)):
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
