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
def __static_MZTrafoModel_enumToName(int mt ):
    """
    __static_MZTrafoModel_enumToName(mt: int ) -> bytes
    """
    assert mt in [0, 1, 2, 3, 4], 'arg mt wrong type'

    cdef libcpp_string _r = _enumToName_MZTrafoModel((<_MZTrafoModel_MODELTYPE>mt))
    py_result = <libcpp_string>_r
    return py_result

def __static_MZTrafoModel_findNearest(list tms , double rt ):
    """
    __static_MZTrafoModel_findNearest(tms: List[MZTrafoModel] , rt: float ) -> int
    """
    assert isinstance(tms, list) and all(isinstance(elemt_rec, MZTrafoModel) for elemt_rec in tms), 'arg tms wrong type'
    assert isinstance(rt, float), 'arg rt wrong type'
    cdef libcpp_vector[_MZTrafoModel] * v0 = new libcpp_vector[_MZTrafoModel]()
    cdef MZTrafoModel item0
    for item0 in tms:
        v0.push_back(deref(item0.inst.get()))

    cdef size_t _r = _findNearest_MZTrafoModel(deref(v0), (<double>rt))
    cdef libcpp_vector[_MZTrafoModel].iterator it_tms = v0.begin()
    replace_0 = []
    while it_tms != v0.end():
        item0 = MZTrafoModel.__new__(MZTrafoModel)
        item0.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel(deref(it_tms)))
        replace_0.append(item0)
        inc(it_tms)
    tms[:] = replace_0
    del v0
    py_result = <size_t>_r
    return py_result

def __static_ExperimentalDesign_fromConsensusMap(ConsensusMap c ):
    """
    __static_ExperimentalDesign_fromConsensusMap(c: ConsensusMap ) -> ExperimentalDesign
    """
    assert isinstance(c, ConsensusMap), 'arg c wrong type'

    cdef _ExperimentalDesign * _r = new _ExperimentalDesign(_fromConsensusMap_ExperimentalDesign((deref(c.inst.get()))))
    cdef ExperimentalDesign py_result = ExperimentalDesign.__new__(ExperimentalDesign)
    py_result.inst = shared_ptr[_ExperimentalDesign](_r)
    return py_result

def __static_ExperimentalDesign_fromFeatureMap(FeatureMap f ):
    """
    __static_ExperimentalDesign_fromFeatureMap(f: FeatureMap ) -> ExperimentalDesign
    """
    assert isinstance(f, FeatureMap), 'arg f wrong type'

    cdef _ExperimentalDesign * _r = new _ExperimentalDesign(_fromFeatureMap_ExperimentalDesign((deref(f.inst.get()))))
    cdef ExperimentalDesign py_result = ExperimentalDesign.__new__(ExperimentalDesign)
    py_result.inst = shared_ptr[_ExperimentalDesign](_r)
    return py_result

def __static_ExperimentalDesign_fromIdentifications(list proteins ):
    """
    __static_ExperimentalDesign_fromIdentifications(proteins: List[ProteinIdentification] ) -> ExperimentalDesign
    """
    assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
    cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
    cdef ProteinIdentification item0
    for item0 in proteins:
        v0.push_back(deref(item0.inst.get()))
    cdef _ExperimentalDesign * _r = new _ExperimentalDesign(_fromIdentifications_ExperimentalDesign(deref(v0)))
    cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v0.begin()
    replace_0 = []
    while it_proteins != v0.end():
        item0 = ProteinIdentification.__new__(ProteinIdentification)
        item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
        replace_0.append(item0)
        inc(it_proteins)
    proteins[:] = replace_0
    del v0
    cdef ExperimentalDesign py_result = ExperimentalDesign.__new__(ExperimentalDesign)
    py_result.inst = shared_ptr[_ExperimentalDesign](_r)
    return py_result

def __static_MZTrafoModel_isValidModel(MZTrafoModel trafo ):
    """
    __static_MZTrafoModel_isValidModel(trafo: MZTrafoModel ) -> bool
    """
    assert isinstance(trafo, MZTrafoModel), 'arg trafo wrong type'

    cdef bool _r = _isValidModel_MZTrafoModel((deref(trafo.inst.get())))
    py_result = <bool>_r
    return py_result

def __static_MZTrafoModel_nameToEnum(bytes name ):
    """
    __static_MZTrafoModel_nameToEnum(name: bytes ) -> int
    """
    assert isinstance(name, bytes), 'arg name wrong type'

    cdef _MZTrafoModel_MODELTYPE _r = _nameToEnum_MZTrafoModel((<libcpp_string>name))
    py_result = <int>_r
    return py_result

def __static_MZTrafoModel_setCoefficientLimits(double offset , double scale , double power ):
    """
    __static_MZTrafoModel_setCoefficientLimits(offset: float , scale: float , power: float ) -> None
    """
    assert isinstance(offset, float), 'arg offset wrong type'
    assert isinstance(scale, float), 'arg scale wrong type'
    assert isinstance(power, float), 'arg power wrong type'



    _setCoefficientLimits_MZTrafoModel((<double>offset), (<double>scale), (<double>power))

def __static_MZTrafoModel_setRANSACParams(RANSACParam p ):
    """
    __static_MZTrafoModel_setRANSACParams(p: RANSACParam ) -> None
    """
    assert isinstance(p, RANSACParam), 'arg p wrong type'

    _setRANSACParams_MZTrafoModel((deref(p.inst.get()))) 

cdef class BoundaryCondition:
    None
    BC_ZERO_ENDPOINTS = 0
    BC_ZERO_FIRST = 1
    BC_ZERO_SECOND = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class MZTrafoModel_MODELTYPE:
    None
    LINEAR = 0
    LINEAR_WEIGHTED = 1
    QUADRATIC = 2
    QUADRATIC_WEIGHTED = 3
    SIZE_OF_MODELTYPE = 4

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class BSpline2d:
    """
    Cython implementation of _BSpline2d

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BSpline2d.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list x , list y , double wave_length , int boundary_condition ,  num_nodes ):
        """
        __init__(self, x: List[float] , y: List[float] , wave_length: float , boundary_condition: int , num_nodes: int ) -> None
        """
        assert isinstance(x, list) and all(isinstance(elemt_rec, float) for elemt_rec in x), 'arg x wrong type'
        assert isinstance(y, list) and all(isinstance(elemt_rec, float) for elemt_rec in y), 'arg y wrong type'
        assert isinstance(wave_length, float), 'arg wave_length wrong type'
        assert boundary_condition in [0, 1, 2], 'arg boundary_condition wrong type'
        assert isinstance(num_nodes, int) and num_nodes >= 0, 'arg num_nodes wrong type'
        cdef libcpp_vector[double] v0 = x
        cdef libcpp_vector[double] v1 = y
    
    
    
        self.inst = shared_ptr[_BSpline2d](new _BSpline2d(v0, v1, (<double>wave_length), (<_BoundaryCondition>boundary_condition), (<size_t>num_nodes)))
        
        
    
    def solve(self, list y ):
        """
        solve(self, y: List[float] ) -> bool
        Solve the spline curve for a new set of y values. Returns false if the solution fails
        """
        assert isinstance(y, list) and all(isinstance(elemt_rec, float) for elemt_rec in y), 'arg y wrong type'
        cdef libcpp_vector[double] v0 = y
        cdef bool _r = self.inst.get().solve(v0)
        
        py_result = <bool>_r
        return py_result
    
    def eval(self, double x ):
        """
        eval(self, x: float ) -> float
        Returns the evaluation of the smoothed curve at a particular x value. If current state is not ok(), returns zero
        """
        assert isinstance(x, float), 'arg x wrong type'
    
        cdef double _r = self.inst.get().eval((<double>x))
        py_result = <double>_r
        return py_result
    
    def derivative(self, double x ):
        """
        derivative(self, x: float ) -> float
        Returns the first derivative of the spline curve at the given position x. Returns zero if the current state is not ok()
        """
        assert isinstance(x, float), 'arg x wrong type'
    
        cdef double _r = self.inst.get().derivative((<double>x))
        py_result = <double>_r
        return py_result
    
    def ok(self):
        """
        ok(self) -> bool
        Returns whether the spline fit was successful
        """
        cdef bool _r = self.inst.get().ok()
        py_result = <bool>_r
        return py_result
    
    def debug(self, bool enable ):
        """
        debug(self, enable: bool ) -> None
        Enable or disable debug messages from the B-spline library
        """
        assert isinstance(enable, pybool_t), 'arg enable wrong type'
    
        self.inst.get().debug((<bool>enable)) 

cdef class CVReference:
    """
    Cython implementation of _CVReference

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVReference.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVReference rv = CVReference.__new__(CVReference)
       rv.inst = shared_ptr[_CVReference](new _CVReference(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVReference rv = CVReference.__new__(CVReference)
       rv.inst = shared_ptr[_CVReference](new _CVReference(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVReference](new _CVReference())
    
    def _init_1(self, CVReference in_0 ):
        """
        _init_1(self, in_0: CVReference ) -> None
        """
        assert isinstance(in_0, CVReference), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVReference](new _CVReference((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVReference ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVReference)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the CV reference
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the CV reference
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setIdentifier(self,  identifier ):
        """
        setIdentifier(self, identifier: Union[bytes, str, String] ) -> None
        Sets the CV identifier which is referenced
        """
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
    
        self.inst.get().setIdentifier(deref((convString(identifier)).get()))
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        Returns the CV identifier which is referenced
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, CVReference):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVReference other_casted = other
        cdef CVReference self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ExperimentalDesign:
    """
    Cython implementation of _ExperimentalDesign

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign.html>`_

    Representation of an experimental design in OpenMS. Instances can be loaded with the ExperimentalDesignFile class
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ExperimentalDesign rv = ExperimentalDesign.__new__(ExperimentalDesign)
       rv.inst = shared_ptr[_ExperimentalDesign](new _ExperimentalDesign(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ExperimentalDesign rv = ExperimentalDesign.__new__(ExperimentalDesign)
       rv.inst = shared_ptr[_ExperimentalDesign](new _ExperimentalDesign(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ExperimentalDesign](new _ExperimentalDesign())
    
    def _init_1(self, ExperimentalDesign in_0 ):
        """
        _init_1(self, in_0: ExperimentalDesign ) -> None
        """
        assert isinstance(in_0, ExperimentalDesign), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ExperimentalDesign](new _ExperimentalDesign((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ExperimentalDesign ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ExperimentalDesign)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMSFileSection(self):
        """
        getMSFileSection(self) -> List[ExperimentalDesign_MSFileSectionEntry]
        """
        _r = self.inst.get().getMSFileSection()
        py_result = []
        cdef libcpp_vector[_ExperimentalDesign_MSFileSectionEntry].iterator it__r = _r.begin()
        cdef ExperimentalDesign_MSFileSectionEntry item_py_result
        while it__r != _r.end():
           item_py_result = ExperimentalDesign_MSFileSectionEntry.__new__(ExperimentalDesign_MSFileSectionEntry)
           item_py_result.inst = shared_ptr[_ExperimentalDesign_MSFileSectionEntry](new _ExperimentalDesign_MSFileSectionEntry(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setMSFileSection(self, list msfile_section ):
        """
        setMSFileSection(self, msfile_section: List[ExperimentalDesign_MSFileSectionEntry] ) -> None
        """
        assert isinstance(msfile_section, list) and all(isinstance(elemt_rec, ExperimentalDesign_MSFileSectionEntry) for elemt_rec in msfile_section), 'arg msfile_section wrong type'
        cdef libcpp_vector[_ExperimentalDesign_MSFileSectionEntry] * v0 = new libcpp_vector[_ExperimentalDesign_MSFileSectionEntry]()
        cdef ExperimentalDesign_MSFileSectionEntry item0
        for item0 in msfile_section:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setMSFileSection(deref(v0))
        del v0
    
    def getSampleSection(self):
        """
        getSampleSection(self) -> ExperimentalDesign_SampleSection
        Returns the Sample Section of the experimental design file
        """
        cdef _ExperimentalDesign_SampleSection * _r = new _ExperimentalDesign_SampleSection(self.inst.get().getSampleSection())
        cdef ExperimentalDesign_SampleSection py_result = ExperimentalDesign_SampleSection.__new__(ExperimentalDesign_SampleSection)
        py_result.inst = shared_ptr[_ExperimentalDesign_SampleSection](_r)
        return py_result
    
    def setSampleSection(self, ExperimentalDesign_SampleSection sample_section ):
        """
        setSampleSection(self, sample_section: ExperimentalDesign_SampleSection ) -> None
        Sets the Sample Section of the experimental design file
        """
        assert isinstance(sample_section, ExperimentalDesign_SampleSection), 'arg sample_section wrong type'
    
        self.inst.get().setSampleSection((deref(sample_section.inst.get())))
    
    def getNumberOfSamples(self):
        """
        getNumberOfSamples(self) -> int
        Returns the number of samples measured (= highest sample index)
        """
        cdef unsigned int _r = self.inst.get().getNumberOfSamples()
        py_result = <unsigned int>_r
        return py_result
    
    def getNumberOfFractions(self):
        """
        getNumberOfFractions(self) -> int
        Returns the number of fractions (= highest fraction index)
        """
        cdef unsigned int _r = self.inst.get().getNumberOfFractions()
        py_result = <unsigned int>_r
        return py_result
    
    def getNumberOfLabels(self):
        """
        getNumberOfLabels(self) -> int
        Returns the number of labels per file
        """
        cdef unsigned int _r = self.inst.get().getNumberOfLabels()
        py_result = <unsigned int>_r
        return py_result
    
    def getNumberOfMSFiles(self):
        """
        getNumberOfMSFiles(self) -> int
        Returns the number of MS files (= fractions * fraction_groups)
        """
        cdef unsigned int _r = self.inst.get().getNumberOfMSFiles()
        py_result = <unsigned int>_r
        return py_result
    
    def getNumberOfFractionGroups(self):
        """
        getNumberOfFractionGroups(self) -> int
        Allows to group fraction ids and source files. Return the number of fraction_groups
        """
        cdef unsigned int _r = self.inst.get().getNumberOfFractionGroups()
        py_result = <unsigned int>_r
        return py_result
    
    def getSample(self,  fraction_group ,  label ):
        """
        getSample(self, fraction_group: int , label: int ) -> int
        Returns sample index (depends on fraction_group and label)
        """
        assert isinstance(fraction_group, int), 'arg fraction_group wrong type'
        assert isinstance(label, int), 'arg label wrong type'
    
    
        cdef unsigned int _r = self.inst.get().getSample((<unsigned int>fraction_group), (<unsigned int>label))
        py_result = <unsigned int>_r
        return py_result
    
    def isFractionated(self):
        """
        isFractionated(self) -> bool
        Returns whether at least one fraction_group in this experimental design is fractionated
        """
        cdef bool _r = self.inst.get().isFractionated()
        py_result = <bool>_r
        return py_result
    
    def sameNrOfMSFilesPerFraction(self):
        """
        sameNrOfMSFilesPerFraction(self) -> bool
        Returns if each fraction number is associated with the same number of fraction_group
        """
        cdef bool _r = self.inst.get().sameNrOfMSFilesPerFraction()
        py_result = <bool>_r
        return py_result
    fromConsensusMap = __static_ExperimentalDesign_fromConsensusMap
    fromFeatureMap = __static_ExperimentalDesign_fromFeatureMap
    fromIdentifications = __static_ExperimentalDesign_fromIdentifications 

cdef class ExperimentalDesign_MSFileSectionEntry:
    """
    Cython implementation of _ExperimentalDesign_MSFileSectionEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign_MSFileSectionEntry.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property path:
        def __set__(self, bytes path):
        
            self.inst.get().path = (<libcpp_string>path)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().path
            py_result = <libcpp_string>_r
            return py_result
    
    property fraction_group:
        def __set__(self,  fraction_group):
        
            self.inst.get().fraction_group = (<unsigned int>fraction_group)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().fraction_group
            py_result = <unsigned int>_r
            return py_result
    
    property fraction:
        def __set__(self,  fraction):
        
            self.inst.get().fraction = (<unsigned int>fraction)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().fraction
            py_result = <unsigned int>_r
            return py_result
    
    property label:
        def __set__(self,  label):
        
            self.inst.get().label = (<unsigned int>label)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().label
            py_result = <unsigned int>_r
            return py_result
    
    property sample:
        def __set__(self,  sample):
        
            self.inst.get().sample = (<unsigned int>sample)
        
    
        def __get__(self):
            cdef unsigned int _r = self.inst.get().sample
            py_result = <unsigned int>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ExperimentalDesign_MSFileSectionEntry](new _ExperimentalDesign_MSFileSectionEntry()) 

cdef class ExperimentalDesign_SampleSection:
    """
    Cython implementation of _ExperimentalDesign_SampleSection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign_SampleSection.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ExperimentalDesign_SampleSection](new _ExperimentalDesign_SampleSection())
    
    def _init_1(self, ExperimentalDesign_SampleSection in_0 ):
        """
        _init_1(self, in_0: ExperimentalDesign_SampleSection ) -> None
        """
        assert isinstance(in_0, ExperimentalDesign_SampleSection), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ExperimentalDesign_SampleSection](new _ExperimentalDesign_SampleSection((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ExperimentalDesign_SampleSection ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ExperimentalDesign_SampleSection)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSamples(self):
        """
        getSamples(self) -> Set[bytes]
        Returns a set of all samples that are present in the sample section
        """
        _r = self.inst.get().getSamples()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getFactors(self):
        """
        getFactors(self) -> Set[bytes]
        Returns a set of all factors (column names) that were defined for the sample section
        """
        _r = self.inst.get().getFactors()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def hasSample(self,  sample ):
        """
        hasSample(self, sample: int ) -> bool
        Checks whether sample section has row for a sample number
        """
        assert isinstance(sample, int), 'arg sample wrong type'
    
        cdef bool _r = self.inst.get().hasSample((<unsigned int>sample))
        py_result = <bool>_r
        return py_result
    
    def hasFactor(self,  factor ):
        """
        hasFactor(self, factor: String ) -> bool
        Checks whether Sample Section has a specific factor (i.e. column name)
        """
        assert isinstance(factor, String), 'arg factor wrong type'
    
        cdef bool _r = self.inst.get().hasFactor(deref((<String>factor).inst.get()))
        py_result = <bool>_r
        return py_result
    
    def getFactorValue(self,  sample ,  factor ):
        """
        getFactorValue(self, sample: int , factor: String ) -> Union[bytes, str, String]
        Returns value of factor for given sample and factor name
        """
        assert isinstance(sample, int), 'arg sample wrong type'
        assert isinstance(factor, String), 'arg factor wrong type'
    
    
        cdef _String _r = self.inst.get().getFactorValue((<unsigned int>sample), deref((<String>factor).inst.get()))
        py_result = convOutputString(_r)
        return py_result 

cdef class IMSIsotopeDistribution:
    """
    Cython implementation of _IMSIsotopeDistribution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSIsotopeDistribution_1_1IMSIsotopeDistribution.html>`_

    Represents a distribution of isotopes restricted to the first K elements
    
    Represents a distribution of isotopes of chemical elements as a list
    of peaks each as a pair of mass and abundance. 'IsotopeDistribution'
    unlike 'IsotopeSpecies' has one abundance per a nominal mass.
    Here is an example in the format (mass; abundance %)
    for molecule H2O (values are taken randomly):
    
    - IsotopeDistribution
        (18.00221; 99.03 %)
        (19.00334; 0.8 %)
        (20.00476; 0.17 %)
    
    - IsotopeSpecies
        (18.00197; 98.012 %)
        (18.00989; 1.018 %)
        (19.00312; 0.683 %)
        (19.00531; 0.117 %)
        (20.00413; 0.134 %)
        (20.00831; 0.036 %)
    
    To the sake of faster computations distribution is restricted
    to the first K elements, where K can be set by adjusting size
    'SIZE' of distribution. @note For the elements most abundant in
    living beings (CHNOPS) this restriction is negligible, since abundances
    decrease dramatically in isotopes order and are usually of no interest
    starting from +10 isotope.
    
    'IsotopeDistribution' implements folding with other distribution using an
    algorithm described in details in paper:
    Boecker et al. "Decomposing metabolic isotope patterns" WABI 2006. doi: 10.1007/11851561_2
    
    Folding with itself is done using Russian Multiplication Scheme
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property ABUNDANCES_SUM_ERROR:
        def __set__(self, double ABUNDANCES_SUM_ERROR):
        
            self.inst.get().ABUNDANCES_SUM_ERROR = (<double>ABUNDANCES_SUM_ERROR)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ABUNDANCES_SUM_ERROR
            py_result = <double>_r
            return py_result
    
    property SIZE:
        def __set__(self,  SIZE):
        
            self.inst.get().SIZE = (<int>SIZE)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().SIZE
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef IMSIsotopeDistribution rv = IMSIsotopeDistribution.__new__(IMSIsotopeDistribution)
       rv.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IMSIsotopeDistribution rv = IMSIsotopeDistribution.__new__(IMSIsotopeDistribution)
       rv.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution())
    
    def _init_1(self, IMSIsotopeDistribution in_0 ):
        """
        _init_1(self, in_0: IMSIsotopeDistribution ) -> None
        """
        assert isinstance(in_0, IMSIsotopeDistribution), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution((deref(in_0.inst.get()))))
    
    def _init_2(self,  nominalMass ):
        """
        _init_2(self, nominalMass: int ) -> None
        """
        assert isinstance(nominalMass, int), 'arg nominalMass wrong type'
    
        self.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution((<unsigned int>nominalMass)))
    
    def _init_3(self, double mass ):
        """
        _init_3(self, mass: float ) -> None
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution((<double>mass)))
    
    def _init_4(self, list peaks ,  nominalMass ):
        """
        _init_4(self, peaks: List[IMSIsotopeDistribution_Peak] , nominalMass: int ) -> None
        """
        assert isinstance(peaks, list) and all(isinstance(elemt_rec, IMSIsotopeDistribution_Peak) for elemt_rec in peaks), 'arg peaks wrong type'
        assert isinstance(nominalMass, int), 'arg nominalMass wrong type'
        cdef libcpp_vector[_IMSIsotopeDistribution_Peak] * v0 = new libcpp_vector[_IMSIsotopeDistribution_Peak]()
        cdef IMSIsotopeDistribution_Peak item0
        for item0 in peaks:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_IMSIsotopeDistribution](new _IMSIsotopeDistribution(deref(v0), (<unsigned int>nominalMass)))
        cdef libcpp_vector[_IMSIsotopeDistribution_Peak].iterator it_peaks = v0.begin()
        replace_0 = []
        while it_peaks != v0.end():
            item0 = IMSIsotopeDistribution_Peak.__new__(IMSIsotopeDistribution_Peak)
            item0.inst = shared_ptr[_IMSIsotopeDistribution_Peak](new _IMSIsotopeDistribution_Peak(deref(it_peaks)))
            replace_0.append(item0)
            inc(it_peaks)
        peaks[:] = replace_0
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IMSIsotopeDistribution ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, nominalMass: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, mass: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, peaks: List[IMSIsotopeDistribution_Peak] , nominalMass: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSIsotopeDistribution)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
             self._init_2(*args)
        elif (len(args)==1) and (isinstance(args[0], float)):
             self._init_3(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, IMSIsotopeDistribution_Peak) for elemt_rec in args[0])) and (isinstance(args[1], int)):
             self._init_4(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def getMass(self,  i ):
        """
        getMass(self, i: int ) -> float
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef double _r = self.inst.get().getMass((<int>i))
        py_result = <double>_r
        return py_result
    
    def getAbundance(self,  i ):
        """
        getAbundance(self, i: int ) -> float
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef double _r = self.inst.get().getAbundance((<int>i))
        py_result = <double>_r
        return py_result
    
    def getAverageMass(self):
        """
        getAverageMass(self) -> float
        """
        cdef double _r = self.inst.get().getAverageMass()
        py_result = <double>_r
        return py_result
    
    def getNominalMass(self):
        """
        getNominalMass(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getNominalMass()
        py_result = <unsigned int>_r
        return py_result
    
    def setNominalMass(self,  nominalMass ):
        """
        setNominalMass(self, nominalMass: int ) -> None
        """
        assert isinstance(nominalMass, int), 'arg nominalMass wrong type'
    
        self.inst.get().setNominalMass((<unsigned int>nominalMass))
    
    def getMasses(self):
        """
        getMasses(self) -> List[float]
        Gets a mass of isotope 'i'
        """
        _r = self.inst.get().getMasses()
        cdef list py_result = _r
        return py_result
    
    def getAbundances(self):
        """
        getAbundances(self) -> List[float]
        Gets an abundance of isotope 'i'
        """
        _r = self.inst.get().getAbundances()
        cdef list py_result = _r
        return py_result
    
    def normalize(self):
        """
        normalize(self) -> None
        Normalizes distribution, i.e. scaling abundances to be summed up to 1 with an error
        """
        self.inst.get().normalize()
    
    def empty(self):
        """
        empty(self) -> bool
        Returns true if the distribution has no peaks, false - otherwise
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, IMSIsotopeDistribution):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IMSIsotopeDistribution other_casted = other
        cdef IMSIsotopeDistribution self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class IMSIsotopeDistribution_Peak:
    """
    Cython implementation of _IMSIsotopeDistribution_Peak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSIsotopeDistribution_1_1IMSIsotopeDistribution_Peak.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property mass:
        def __set__(self, double mass):
        
            self.inst.get().mass = (<double>mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mass
            py_result = <double>_r
            return py_result
    
    property abundance:
        def __set__(self, double abundance):
        
            self.inst.get().abundance = (<double>abundance)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().abundance
            py_result = <double>_r
            return py_result
    
    def __init__(self, double mass , double abundance ):
        """
        __init__(self, mass: float , abundance: float ) -> None
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(abundance, float), 'arg abundance wrong type'
    
    
        self.inst = shared_ptr[_IMSIsotopeDistribution_Peak](new _IMSIsotopeDistribution_Peak((<double>mass), (<double>abundance)))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, IMSIsotopeDistribution_Peak):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IMSIsotopeDistribution_Peak other_casted = other
        cdef IMSIsotopeDistribution_Peak self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class LinearInterpolation:
    """
    Cython implementation of _LinearInterpolation[double,double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1LinearInterpolation[double,double].html>`_

    Provides access to linearly interpolated values (and
    derivatives) from discrete data points.  Values beyond the given range
    of data points are implicitly taken as zero.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef LinearInterpolation rv = LinearInterpolation.__new__(LinearInterpolation)
       rv.inst = shared_ptr[_LinearInterpolation[double,double]](new _LinearInterpolation[double,double](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LinearInterpolation rv = LinearInterpolation.__new__(LinearInterpolation)
       rv.inst = shared_ptr[_LinearInterpolation[double,double]](new _LinearInterpolation[double,double](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LinearInterpolation[double,double]](new _LinearInterpolation[double,double]())
    
    def _init_1(self, LinearInterpolation in_0 ):
        """
        _init_1(self, in_0: LinearInterpolation ) -> None
        """
        assert isinstance(in_0, LinearInterpolation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LinearInterpolation[double,double]](new _LinearInterpolation[double,double]((deref(in_0.inst.get()))))
    
    def _init_2(self, double scale , double offset ):
        """
        _init_2(self, scale: float , offset: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
        assert isinstance(offset, float), 'arg offset wrong type'
    
    
        self.inst = shared_ptr[_LinearInterpolation[double,double]](new _LinearInterpolation[double,double]((<double>scale), (<double>offset)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LinearInterpolation ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, scale: float , offset: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LinearInterpolation)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def value(self, double arg_pos ):
        """
        value(self, arg_pos: float ) -> float
        Returns the interpolated value
        """
        assert isinstance(arg_pos, float), 'arg arg_pos wrong type'
    
        cdef double _r = self.inst.get().value((<double>arg_pos))
        py_result = <double>_r
        return py_result
    
    def addValue(self, double arg_pos , double arg_value ):
        """
        addValue(self, arg_pos: float , arg_value: float ) -> None
        Performs linear resampling. The `arg_value` is split up and added to the data points around `arg_pos`
        """
        assert isinstance(arg_pos, float), 'arg arg_pos wrong type'
        assert isinstance(arg_value, float), 'arg arg_value wrong type'
    
    
        self.inst.get().addValue((<double>arg_pos), (<double>arg_value))
    
    def derivative(self, double arg_pos ):
        """
        derivative(self, arg_pos: float ) -> float
        Returns the interpolated derivative
        """
        assert isinstance(arg_pos, float), 'arg arg_pos wrong type'
    
        cdef double _r = self.inst.get().derivative((<double>arg_pos))
        py_result = <double>_r
        return py_result
    
    def getData(self):
        """
        getData(self) -> List[float]
        Returns the internal random access container from which interpolated values are being sampled
        """
        _r = self.inst.get().getData()
        cdef list py_result = _r
        return py_result
    
    def setData(self, list data ):
        """
        setData(self, data: List[float] ) -> None
        Assigns data to the internal random access container from which interpolated values are being sampled
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, float) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[double] v0 = data
        self.inst.get().setData(v0)
        data[:] = v0
    
    def empty(self):
        """
        empty(self) -> bool
        Returns `true` if getData() is empty
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def key2index(self, double pos ):
        """
        key2index(self, pos: float ) -> float
        The transformation from "outside" to "inside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().key2index((<double>pos))
        py_result = <double>_r
        return py_result
    
    def index2key(self, double pos ):
        """
        index2key(self, pos: float ) -> float
        The transformation from "inside" to "outside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().index2key((<double>pos))
        py_result = <double>_r
        return py_result
    
    def getScale(self):
        """
        getScale(self) -> float
        "Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        """
        cdef double _r = self.inst.get().getScale()
        py_result = <double>_r
        return py_result
    
    def setScale(self, double scale ):
        """
        setScale(self, scale: float ) -> None
        "Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        """
        assert isinstance(scale, float), 'arg scale wrong type'
    
        self.inst.get().setScale((<double &>scale))
    
    def getOffset(self):
        """
        getOffset(self) -> float
        "Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        """
        cdef double _r = self.inst.get().getOffset()
        py_result = <double>_r
        return py_result
    
    def setOffset(self, double offset ):
        """
        setOffset(self, offset: float ) -> None
        "Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset((<double &>offset))
    
    def _setMapping_0(self, double scale , double inside , double outside ):
        """
        _setMapping_0(self, scale: float , inside: float , outside: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
        assert isinstance(inside, float), 'arg inside wrong type'
        assert isinstance(outside, float), 'arg outside wrong type'
    
    
    
        self.inst.get().setMapping((<double &>scale), (<double &>inside), (<double &>outside))
    
    def _setMapping_1(self, double inside_low , double outside_low , double inside_high , double outside_high ):
        """
        _setMapping_1(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
        """
        assert isinstance(inside_low, float), 'arg inside_low wrong type'
        assert isinstance(outside_low, float), 'arg outside_low wrong type'
        assert isinstance(inside_high, float), 'arg inside_high wrong type'
        assert isinstance(outside_high, float), 'arg outside_high wrong type'
    
    
    
    
        self.inst.get().setMapping((<double &>inside_low), (<double &>outside_low), (<double &>inside_high), (<double &>outside_high))
    
    def setMapping(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMapping(self, scale: float , inside: float , outside: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMapping(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._setMapping_0(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._setMapping_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getInsideReferencePoint(self):
        """
        getInsideReferencePoint(self) -> float
        """
        cdef double _r = self.inst.get().getInsideReferencePoint()
        py_result = <double>_r
        return py_result
    
    def getOutsideReferencePoint(self):
        """
        getOutsideReferencePoint(self) -> float
        """
        cdef double _r = self.inst.get().getOutsideReferencePoint()
        py_result = <double>_r
        return py_result
    
    def supportMin(self):
        """
        supportMin(self) -> float
        """
        cdef double _r = self.inst.get().supportMin()
        py_result = <double>_r
        return py_result
    
    def supportMax(self):
        """
        supportMax(self) -> float
        """
        cdef double _r = self.inst.get().supportMax()
        py_result = <double>_r
        return py_result 

cdef class LinearResamplerAlign:
    """
    Cython implementation of _LinearResamplerAlign

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LinearResamplerAlign.html>`_
      -- Inherits from ['LinearResampler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef LinearResamplerAlign rv = LinearResamplerAlign.__new__(LinearResamplerAlign)
       rv.inst = shared_ptr[_LinearResamplerAlign](new _LinearResamplerAlign(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LinearResamplerAlign rv = LinearResamplerAlign.__new__(LinearResamplerAlign)
       rv.inst = shared_ptr[_LinearResamplerAlign](new _LinearResamplerAlign(deref(self.inst.get())))
       return rv
    
    def __init__(self, LinearResamplerAlign in_0 ):
        """
        __init__(self, in_0: LinearResamplerAlign ) -> None
        """
        assert isinstance(in_0, LinearResamplerAlign), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LinearResamplerAlign](new _LinearResamplerAlign((deref(in_0.inst.get()))))
    
    def raster(self, MSSpectrum input ):
        """
        raster(self, input: MSSpectrum ) -> None
        Applies the resampling algorithm to an MSSpectrum
        """
        assert isinstance(input, MSSpectrum), 'arg input wrong type'
    
        self.inst.get().raster((deref(input.inst.get())))
    
    def rasterExperiment(self, MSExperiment input ):
        """
        rasterExperiment(self, input: MSExperiment ) -> None
        Resamples the data in an MSExperiment
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
    
        self.inst.get().rasterExperiment((deref(input.inst.get())))
    
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

cdef class MRMFeatureFinderScoring:
    """
    Cython implementation of _MRMFeatureFinderScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureFinderScoring.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MRMFeatureFinderScoring](new _MRMFeatureFinderScoring())
    
    def pickExperiment(self, MSExperiment chromatograms , FeatureMap output , TargetedExperiment transition_exp_ , TransformationDescription trafo , MSExperiment swath_map ):
        """
        pickExperiment(self, chromatograms: MSExperiment , output: FeatureMap , transition_exp_: TargetedExperiment , trafo: TransformationDescription , swath_map: MSExperiment ) -> None
        Pick features in one experiment containing chromatogram
        
        Function for for wrapping in Python, only uses OpenMS datastructures and does not return the map
        
        
        :param chromatograms: The input chromatograms
        :param output: The output features with corresponding scores
        :param transition_exp: The transition list describing the experiment
        :param trafo: Optional transformation of the experimental retention time to the normalized retention time space used in the transition list
        :param swath_map: Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted
        """
        assert isinstance(chromatograms, MSExperiment), 'arg chromatograms wrong type'
        assert isinstance(output, FeatureMap), 'arg output wrong type'
        assert isinstance(transition_exp_, TargetedExperiment), 'arg transition_exp_ wrong type'
        assert isinstance(trafo, TransformationDescription), 'arg trafo wrong type'
        assert isinstance(swath_map, MSExperiment), 'arg swath_map wrong type'
    
    
    
    
    
        self.inst.get().pickExperiment((deref(chromatograms.inst.get())), (deref(output.inst.get())), (deref(transition_exp_.inst.get())), (deref(trafo.inst.get())), (deref(swath_map.inst.get())))
    
    def setStrictFlag(self, bool flag ):
        """
        setStrictFlag(self, flag: bool ) -> None
        """
        assert isinstance(flag, pybool_t), 'arg flag wrong type'
    
        self.inst.get().setStrictFlag((<bool>flag))
    
    def _setMS1Map_0(self, SpectrumAccessOpenMS ms1_map ):
        """
        _setMS1Map_0(self, ms1_map: SpectrumAccessOpenMS ) -> None
        """
        assert isinstance(ms1_map, SpectrumAccessOpenMS), 'arg ms1_map wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMS] input_ms1_map = ms1_map.inst
        self.inst.get().setMS1Map(input_ms1_map)
    
    def _setMS1Map_1(self, SpectrumAccessOpenMSCached ms1_map ):
        """
        _setMS1Map_1(self, ms1_map: SpectrumAccessOpenMSCached ) -> None
        """
        assert isinstance(ms1_map, SpectrumAccessOpenMSCached), 'arg ms1_map wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMSCached] input_ms1_map = ms1_map.inst
        self.inst.get().setMS1Map(input_ms1_map)
    
    def setMS1Map(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMS1Map(self, ms1_map: SpectrumAccessOpenMS ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMS1Map(self, ms1_map: SpectrumAccessOpenMSCached ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMS)):
            return self._setMS1Map_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMSCached)):
            return self._setMS1Map_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def scorePeakgroups(self, LightMRMTransitionGroupCP transition_group , TransformationDescription trafo , list swath_maps , FeatureMap output , bool ms1only ):
        """
        scorePeakgroups(self, transition_group: LightMRMTransitionGroupCP , trafo: TransformationDescription , swath_maps: List[SwathMap] , output: FeatureMap , ms1only: bool ) -> None
        Score all peak groups of a transition group
        
        Iterate through all features found along the chromatograms of the transition group and score each one individually
        
        
        :param transition_group: The MRMTransitionGroup to be scored (input)
        :param trafo: Optional transformation of the experimental retention time
            to the normalized retention time space used in thetransition list
        :param swath_maps: Optional SWATH-MS (DIA) map corresponding from which
            the chromatograms were extracted. Use empty map if no data is available
        :param output: The output features with corresponding scores (the found
            features will be added to this FeatureMap)
        :param ms1only: Whether to only do MS1 scoring and skip all MS2 scoring
        """
        assert isinstance(transition_group, LightMRMTransitionGroupCP), 'arg transition_group wrong type'
        assert isinstance(trafo, TransformationDescription), 'arg trafo wrong type'
        assert isinstance(swath_maps, list) and all(isinstance(elemt_rec, SwathMap) for elemt_rec in swath_maps), 'arg swath_maps wrong type'
        assert isinstance(output, FeatureMap), 'arg output wrong type'
        assert isinstance(ms1only, pybool_t), 'arg ms1only wrong type'
    
    
        cdef libcpp_vector[_SwathMap] * v2 = new libcpp_vector[_SwathMap]()
        cdef SwathMap item2
        for item2 in swath_maps:
            v2.push_back(deref(item2.inst.get()))
    
    
        self.inst.get().scorePeakgroups((deref(transition_group.inst.get())), (deref(trafo.inst.get())), deref(v2), (deref(output.inst.get())), (<bool>ms1only))
        del v2
    
    def prepareProteinPeptideMaps_(self, LightTargetedExperiment transition_exp ):
        """
        prepareProteinPeptideMaps_(self, transition_exp: LightTargetedExperiment ) -> None
        Prepares the internal mappings of peptides and proteins
        
        Calling this method _is_ required before calling scorePeakgroups
        
        
        :param transition_exp: The transition list describing the experiment
        """
        assert isinstance(transition_exp, LightTargetedExperiment), 'arg transition_exp wrong type'
    
        self.inst.get().prepareProteinPeptideMaps_((deref(transition_exp.inst.get())))
    
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

cdef class MZTrafoModel:
    """
    Cython implementation of _MZTrafoModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MZTrafoModel.html>`_

    Create and apply models of a mass recalibration function
    
    The input is a list of calibration points (ideally spanning a wide m/z range to prevent extrapolation when applying to model)
    
    Models (LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED) can be trained using CalData points (or a subset of them)
    Calibration points can have different retention time points, and a model should be build such that it captures
    the local (in time) decalibration of the instrument, i.e. choose appropriate time windows along RT to calibrate the
    spectra in this RT region
    From the available calibrant data, a model is build. Later, any uncalibrated m/z value can be fed to the model, to obtain
    a calibrated m/z
    
    The input domain can either be absolute mass differences in [Th], or relative differences in [ppm]
    The models are build based on this input
    
    Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MZTrafoModel rv = MZTrafoModel.__new__(MZTrafoModel)
       rv.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MZTrafoModel rv = MZTrafoModel.__new__(MZTrafoModel)
       rv.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel())
    
    def _init_1(self, MZTrafoModel in_0 ):
        """
        _init_1(self, in_0: MZTrafoModel ) -> None
        """
        assert isinstance(in_0, MZTrafoModel), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel((deref(in_0.inst.get()))))
    
    def _init_2(self, bool in_0 ):
        """
        _init_2(self, in_0: bool ) -> None
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MZTrafoModel](new _MZTrafoModel((<bool>in_0)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MZTrafoModel ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MZTrafoModel)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], pybool_t)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isTrained(self):
        """
        isTrained(self) -> bool
        Returns true if the model have coefficients (i.e. was trained successfully)
        """
        cdef bool _r = self.inst.get().isTrained()
        py_result = <bool>_r
        return py_result
    
    def getRT(self):
        """
        getRT(self) -> float
        Get RT associated with the model (training region)
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def predict(self, double mz ):
        """
        predict(self, mz: float ) -> float
        Apply the model to an uncalibrated m/z value
        
        Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!
        
        Applies the function y = intercept + slope*mz + power*mz^2
        and returns y
        
        
        :param mz: The uncalibrated m/z value
        :return: The calibrated m/z value
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef double _r = self.inst.get().predict((<double>mz))
        py_result = <double>_r
        return py_result
    
    def _train_0(self, CalibrationData cd , int md , bool use_RANSAC , double rt_left , double rt_right ):
        """
        _train_0(self, cd: CalibrationData , md: int , use_RANSAC: bool , rt_left: float , rt_right: float ) -> bool
        Train a model using calibrant data
        
        If the CalibrationData was created using peak groups (usually corresponding to mass traces),
        the median for each group is used as a group representative. This
        is more robust, and reduces the number of data points drastically, i.e. one value per group
        
        Internally, these steps take place:
        - apply RT filter
        - [compute median per group] (only if groups were given in 'cd')
        - set Model's rt position
        - call train() (see overloaded method)
        
        
        :param cd: List of calibrants
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :param rt_left: Filter 'cd' by RT; all calibrants with RT < 'rt_left' are removed
        :param rt_right: Filter 'cd' by RT; all calibrants with RT > 'rt_right' are removed
        :return: True if model was build, false otherwise
        """
        assert isinstance(cd, CalibrationData), 'arg cd wrong type'
        assert md in [0, 1, 2, 3, 4], 'arg md wrong type'
        assert isinstance(use_RANSAC, pybool_t), 'arg use_RANSAC wrong type'
        assert isinstance(rt_left, float), 'arg rt_left wrong type'
        assert isinstance(rt_right, float), 'arg rt_right wrong type'
    
    
    
    
    
        cdef bool _r = self.inst.get().train((deref(cd.inst.get())), (<_MZTrafoModel_MODELTYPE>md), (<bool>use_RANSAC), (<double>rt_left), (<double>rt_right))
        py_result = <bool>_r
        return py_result
    
    def _train_1(self, list error_mz , list theo_mz , list weights , int md , bool use_RANSAC ):
        """
        _train_1(self, error_mz: List[float] , theo_mz: List[float] , weights: List[float] , md: int , use_RANSAC: bool ) -> bool
        Train a model using calibrant data
        
        Given theoretical and observed mass values (and corresponding weights),
        a model (linear, quadratic, ...) is build
        Outlier removal is applied before
        The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm]
        The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
        rather on how 'predict()' works internally
        
        Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
        
        Internally, these steps take place:
        - [apply RANSAC] (depending on 'use_RANSAC')
        - build model and store its parameters internally
        
        
        :param error_mz: Observed Mass error (in ppm or Th)
        :param theo_mz: Theoretical m/z values, corresponding to 'error_mz'
        :param weights: For weighted models only: weight of calibrants; ignored otherwise
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :return: True if model was build, false otherwise
        """
        assert isinstance(error_mz, list) and all(isinstance(elemt_rec, float) for elemt_rec in error_mz), 'arg error_mz wrong type'
        assert isinstance(theo_mz, list) and all(isinstance(elemt_rec, float) for elemt_rec in theo_mz), 'arg theo_mz wrong type'
        assert isinstance(weights, list) and all(isinstance(elemt_rec, float) for elemt_rec in weights), 'arg weights wrong type'
        assert md in [0, 1, 2, 3, 4], 'arg md wrong type'
        assert isinstance(use_RANSAC, pybool_t), 'arg use_RANSAC wrong type'
        cdef libcpp_vector[double] v0 = error_mz
        cdef libcpp_vector[double] v1 = theo_mz
        cdef libcpp_vector[double] v2 = weights
    
    
        cdef bool _r = self.inst.get().train(v0, v1, v2, (<_MZTrafoModel_MODELTYPE>md), (<bool>use_RANSAC))
        
        
        
        py_result = <bool>_r
        return py_result
    
    def train(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: train(self, cd: CalibrationData , md: int , use_RANSAC: bool , rt_left: float , rt_right: float ) -> bool
          :noindex:
        
        Train a model using calibrant data
        
        If the CalibrationData was created using peak groups (usually corresponding to mass traces),
        the median for each group is used as a group representative. This
        is more robust, and reduces the number of data points drastically, i.e. one value per group
        
        Internally, these steps take place:
        - apply RT filter
        - [compute median per group] (only if groups were given in 'cd')
        - set Model's rt position
        - call train() (see overloaded method)
        
        
        :param cd: List of calibrants
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :param rt_left: Filter 'cd' by RT; all calibrants with RT < 'rt_left' are removed
        :param rt_right: Filter 'cd' by RT; all calibrants with RT > 'rt_right' are removed
        :return: True if model was build, false otherwise
        
        .. rubric:: Overload:
        .. py:function:: train(self, error_mz: List[float] , theo_mz: List[float] , weights: List[float] , md: int , use_RANSAC: bool ) -> bool
          :noindex:
        
        Train a model using calibrant data
        
        Given theoretical and observed mass values (and corresponding weights),
        a model (linear, quadratic, ...) is build
        Outlier removal is applied before
        The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm]
        The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
        rather on how 'predict()' works internally
        
        Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
        
        Internally, these steps take place:
        - [apply RANSAC] (depending on 'use_RANSAC')
        - build model and store its parameters internally
        
        
        :param error_mz: Observed Mass error (in ppm or Th)
        :param theo_mz: Theoretical m/z values, corresponding to 'error_mz'
        :param weights: For weighted models only: weight of calibrants; ignored otherwise
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :return: True if model was build, false otherwise
    
        """
        if (len(args)==5) and (isinstance(args[0], CalibrationData)) and (args[1] in [0, 1, 2, 3, 4]) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], float)) and (isinstance(args[4], float)):
            return self._train_0(*args)
        elif (len(args)==5) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])) and (isinstance(args[2], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[2])) and (args[3] in [0, 1, 2, 3, 4]) and (isinstance(args[4], pybool_t)):
            return self._train_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCoefficients(self, double intercept , double slope , double power ):
        """
        getCoefficients(self, intercept: float , slope: float , power: float ) -> None
        Get model coefficients
        
        Parameters will be filled with internal model parameters
        The model must be trained before; Exception is thrown otherwise!
        
        
        :param intercept: The intercept
        :param slope: The slope
        :param power: The coefficient for x*x (will be 0 for linear models)
        """
        assert isinstance(intercept, float), 'arg intercept wrong type'
        assert isinstance(slope, float), 'arg slope wrong type'
        assert isinstance(power, float), 'arg power wrong type'
    
    
    
        self.inst.get().getCoefficients((<double &>intercept), (<double &>slope), (<double &>power))
    
    def _setCoefficients_0(self, MZTrafoModel in_0 ):
        """
        _setCoefficients_0(self, in_0: MZTrafoModel ) -> None
        Copy model coefficients from another model
        """
        assert isinstance(in_0, MZTrafoModel), 'arg in_0 wrong type'
    
        self.inst.get().setCoefficients((deref(in_0.inst.get())))
    
    def _setCoefficients_1(self, double in_0 , double in_1 , double in_2 ):
        """
        _setCoefficients_1(self, in_0: float , in_1: float , in_2: float ) -> None
        Manually set model coefficients
        
        Can be used instead of train(), so manually set coefficients
        It must be exactly three values. If you want a linear model, set 'power' to zero
        If you want a constant model, set slope to zero in addition
        
        
        :param intercept: The offset
        :param slope: The slope
        :param power: The x*x coefficient (for quadratic models)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
        assert isinstance(in_2, float), 'arg in_2 wrong type'
    
    
    
        self.inst.get().setCoefficients((<double>in_0), (<double>in_1), (<double>in_2))
    
    def setCoefficients(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setCoefficients(self, in_0: MZTrafoModel ) -> None
          :noindex:
        
        Copy model coefficients from another model

        
        .. rubric:: Overload:
        .. py:function:: setCoefficients(self, in_0: float , in_1: float , in_2: float ) -> None
          :noindex:
        
        Manually set model coefficients
        
        Can be used instead of train(), so manually set coefficients
        It must be exactly three values. If you want a linear model, set 'power' to zero
        If you want a constant model, set slope to zero in addition
        
        
        :param intercept: The offset
        :param slope: The slope
        :param power: The x*x coefficient (for quadratic models)
    
        """
        if (len(args)==1) and (isinstance(args[0], MZTrafoModel)):
            return self._setCoefficients_0(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._setCoefficients_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    enumToName = __static_MZTrafoModel_enumToName
    findNearest = __static_MZTrafoModel_findNearest
    isValidModel = __static_MZTrafoModel_isValidModel
    nameToEnum = __static_MZTrafoModel_nameToEnum
    setCoefficientLimits = __static_MZTrafoModel_setCoefficientLimits
    setRANSACParams = __static_MZTrafoModel_setRANSACParams 

cdef class MascotXMLFile:
    """
    Cython implementation of _MascotXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MascotXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MascotXMLFile rv = MascotXMLFile.__new__(MascotXMLFile)
       rv.inst = shared_ptr[_MascotXMLFile](new _MascotXMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MascotXMLFile rv = MascotXMLFile.__new__(MascotXMLFile)
       rv.inst = shared_ptr[_MascotXMLFile](new _MascotXMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MascotXMLFile](new _MascotXMLFile())
    
    def _init_1(self, MascotXMLFile in_0 ):
        """
        _init_1(self, in_0: MascotXMLFile ) -> None
        """
        assert isinstance(in_0, MascotXMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MascotXMLFile](new _MascotXMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MascotXMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MascotXMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , ProteinIdentification protein_identification , PeptideIdentificationList id_data , SpectrumMetaDataLookup rt_mapping ):
        """
        load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , rt_mapping: SpectrumMetaDataLookup ) -> None
        Loads data from a Mascot XML file
        
        
        :param filename: The file to be loaded
        :param protein_identification: Protein identifications belonging to the whole experiment
        :param id_data: The identifications with m/z and RT
        :param lookup: Helper object for looking up spectrum meta data
        :raises:
          Exception: FileNotFound is thrown if the file does not exists
        :raises:
          Exception: ParseError is thrown if the file does not suit to the standard
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(id_data, PeptideIdentificationList), 'arg id_data wrong type'
        assert isinstance(rt_mapping, SpectrumMetaDataLookup), 'arg rt_mapping wrong type'
    
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(protein_identification.inst.get())), (deref(id_data.inst.get())), (deref(rt_mapping.inst.get())))
    
    def initializeLookup(self, SpectrumMetaDataLookup lookup , MSExperiment experiment ,  scan_regex ):
        """
        initializeLookup(self, lookup: SpectrumMetaDataLookup , experiment: MSExperiment , scan_regex: Union[bytes, str, String] ) -> None
        Initializes a helper object for looking up spectrum meta data (RT, m/z)
        
        
        :param lookup: Helper object to initialize
        :param experiment: Experiment containing the spectra
        :param scan_regex: Optional regular expression for extracting information from references to spectra
        """
        assert isinstance(lookup, SpectrumMetaDataLookup), 'arg lookup wrong type'
        assert isinstance(experiment, MSExperiment), 'arg experiment wrong type'
        assert (isinstance(scan_regex, str) or isinstance(scan_regex, bytes) or isinstance(scan_regex, String)), 'arg scan_regex wrong type'
    
    
    
        self.inst.get().initializeLookup((deref(lookup.inst.get())), (deref(experiment.inst.get())), deref((convString(scan_regex)).get()))
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
        return py_result 

cdef class OSWFile:
    """
    Cython implementation of _OSWFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OSWFile.html>`_

    This class serves for reading in and writing OpenSWATH OSW files
    
    See OpenSwathOSWWriter for more functionality
    
    The reader and writer returns data in a format suitable for PercolatorAdapter.
    OSW files have a flexible data structure. They contain all peptide query
    parameters of TraML/PQP files with the detected and quantified features of
    OpenSwathWorkflow (feature, feature_ms1, feature_ms2 & feature_transition)
    
    The OSWFile reader extracts the feature information from the OSW file for
    each level (MS1, MS2 & transition) separately and generates Percolator input
    files. For each of the three Percolator reports, OSWFile writer adds a table
    (score_ms1, score_ms2, score_transition) with the respective confidence metrics.
    These tables can be mapped to the corresponding feature tables, are very similar
    to PyProphet results and can thus be used interchangeably
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OSWFile rv = OSWFile.__new__(OSWFile)
       rv.inst = shared_ptr[_OSWFile](new _OSWFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSWFile rv = OSWFile.__new__(OSWFile)
       rv.inst = shared_ptr[_OSWFile](new _OSWFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  filename ):
        """
        _init_0(self, filename: Union[bytes, str] ) -> None
        """
        assert isinstance(filename, (bytes, str)), 'arg filename wrong type'
        if isinstance(filename, str):
            filename = filename.encode('utf-8')
        self.inst = shared_ptr[_OSWFile](new _OSWFile((<libcpp_string>filename)))
    
    def _init_1(self, OSWFile in_0 ):
        """
        _init_1(self, in_0: OSWFile ) -> None
        """
        assert isinstance(in_0, OSWFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSWFile](new _OSWFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSWFile ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], (bytes, str))):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSWFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeptideHit:
    """
    Cython implementation of _PeptideHit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit.html>`_
      -- Inherits from ['MetaInfoInterface']

    Representation of a peptide hit
    
    It contains the fields score, score_type, rank, and sequence
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideHit rv = PeptideHit.__new__(PeptideHit)
       rv.inst = shared_ptr[_PeptideHit](new _PeptideHit(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideHit rv = PeptideHit.__new__(PeptideHit)
       rv.inst = shared_ptr[_PeptideHit](new _PeptideHit(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideHit](new _PeptideHit())
    
    def _init_1(self, double score ,  rank ,  charge , AASequence sequence ):
        """
        _init_1(self, score: float , rank: int , charge: int , sequence: AASequence ) -> None
        """
        assert isinstance(score, float), 'arg score wrong type'
        assert isinstance(rank, int), 'arg rank wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(sequence, AASequence), 'arg sequence wrong type'
    
    
    
    
        self.inst = shared_ptr[_PeptideHit](new _PeptideHit((<double>score), (<unsigned int>rank), (<int>charge), (deref(sequence.inst.get()))))
    
    def _init_2(self, PeptideHit in_0 ):
        """
        _init_2(self, in_0: PeptideHit ) -> None
        """
        assert isinstance(in_0, PeptideHit), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideHit](new _PeptideHit((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, score: float , rank: int , charge: int , sequence: AASequence ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideHit ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], int)) and (isinstance(args[2], int)) and (isinstance(args[3], AASequence)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideHit)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getScore(self):
        """
        getScore(self) -> float
        Returns the PSM score
        """
        cdef float _r = self.inst.get().getScore()
        py_result = <float>_r
        return py_result
    
    def getRank(self):
        """
        getRank(self) -> int
        Returns the PSM rank
        """
        cdef unsigned int _r = self.inst.get().getRank()
        py_result = <unsigned int>_r
        return py_result
    
    def getSequence(self):
        """
        getSequence(self) -> AASequence
        Returns the peptide sequence without trailing or following spaces
        """
        cdef _AASequence * _r = new _AASequence(self.inst.get().getSequence())
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns the charge of the peptide
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def getPeptideEvidences(self):
        """
        getPeptideEvidences(self) -> List[PeptideEvidence]
        Returns information on peptides (potentially) identified by this PSM
        """
        _r = self.inst.get().getPeptideEvidences()
        py_result = []
        cdef libcpp_vector[_PeptideEvidence].iterator it__r = _r.begin()
        cdef PeptideEvidence item_py_result
        while it__r != _r.end():
           item_py_result = PeptideEvidence.__new__(PeptideEvidence)
           item_py_result.inst = shared_ptr[_PeptideEvidence](new _PeptideEvidence(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setPeptideEvidences(self, list in_0 ):
        """
        setPeptideEvidences(self, in_0: List[PeptideEvidence] ) -> None
        Sets information on peptides (potentially) identified by this PSM
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, PeptideEvidence) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_PeptideEvidence] * v0 = new libcpp_vector[_PeptideEvidence]()
        cdef PeptideEvidence item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setPeptideEvidences(deref(v0))
        del v0
    
    def addPeptideEvidence(self, PeptideEvidence in_0 ):
        """
        addPeptideEvidence(self, in_0: PeptideEvidence ) -> None
        Adds information on a peptide that is (potentially) identified by this PSM
        """
        assert isinstance(in_0, PeptideEvidence), 'arg in_0 wrong type'
    
        self.inst.get().addPeptideEvidence((deref(in_0.inst.get())))
    
    def extractProteinAccessionsSet(self):
        """
        extractProteinAccessionsSet(self) -> Set[bytes]
        Extracts the set of non-empty protein accessions from peptide evidences
        """
        _r = self.inst.get().extractProteinAccessionsSet()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setAnalysisResults(self, list aresult ):
        """
        setAnalysisResults(self, aresult: List[PeptideHit_AnalysisResult] ) -> None
        Sets information on (search engine) sub scores associated with this PSM
        """
        assert isinstance(aresult, list) and all(isinstance(elemt_rec, PeptideHit_AnalysisResult) for elemt_rec in aresult), 'arg aresult wrong type'
        cdef libcpp_vector[_PeptideHit_AnalysisResult] * v0 = new libcpp_vector[_PeptideHit_AnalysisResult]()
        cdef PeptideHit_AnalysisResult item0
        for item0 in aresult:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setAnalysisResults(deref(v0))
        del v0
    
    def addAnalysisResults(self, PeptideHit_AnalysisResult aresult ):
        """
        addAnalysisResults(self, aresult: PeptideHit_AnalysisResult ) -> None
        Add information on (search engine) sub scores associated with this PSM
        """
        assert isinstance(aresult, PeptideHit_AnalysisResult), 'arg aresult wrong type'
    
        self.inst.get().addAnalysisResults((deref(aresult.inst.get())))
    
    def getAnalysisResults(self):
        """
        getAnalysisResults(self) -> List[PeptideHit_AnalysisResult]
        Returns information on (search engine) sub scores associated with this PSM
        """
        _r = self.inst.get().getAnalysisResults()
        py_result = []
        cdef libcpp_vector[_PeptideHit_AnalysisResult].iterator it__r = _r.begin()
        cdef PeptideHit_AnalysisResult item_py_result
        while it__r != _r.end():
           item_py_result = PeptideHit_AnalysisResult.__new__(PeptideHit_AnalysisResult)
           item_py_result.inst = shared_ptr[_PeptideHit_AnalysisResult](new _PeptideHit_AnalysisResult(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setPeakAnnotations(self, list in_0 ):
        """
        setPeakAnnotations(self, in_0: List[PeptideHit_PeakAnnotation] ) -> None
        Sets the fragment annotations
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, PeptideHit_PeakAnnotation) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_PeptideHit_PeakAnnotation] * v0 = new libcpp_vector[_PeptideHit_PeakAnnotation]()
        cdef PeptideHit_PeakAnnotation item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setPeakAnnotations(deref(v0))
        del v0
    
    def getPeakAnnotations(self):
        """
        getPeakAnnotations(self) -> List[PeptideHit_PeakAnnotation]
        Returns the fragment annotations
        """
        _r = self.inst.get().getPeakAnnotations()
        py_result = []
        cdef libcpp_vector[_PeptideHit_PeakAnnotation].iterator it__r = _r.begin()
        cdef PeptideHit_PeakAnnotation item_py_result
        while it__r != _r.end():
           item_py_result = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
           item_py_result.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setScore(self, double in_0 ):
        """
        setScore(self, in_0: float ) -> None
        Sets the PSM score
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setScore((<double>in_0))
    
    def setRank(self,  in_0 ):
        """
        setRank(self, in_0: int ) -> None
        Sets the PSM rank
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().setRank((<unsigned int>in_0))
    
    def setSequence(self, AASequence in_0 ):
        """
        setSequence(self, in_0: AASequence ) -> None
        Sets the peptide sequence
        """
        assert isinstance(in_0, AASequence), 'arg in_0 wrong type'
    
        self.inst.get().setSequence((deref(in_0.inst.get())))
    
    def setCharge(self,  in_0 ):
        """
        setCharge(self, in_0: int ) -> None
        Sets the charge of the peptide
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().setCharge((<int>in_0))
    
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
        if not isinstance(other, PeptideHit):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef PeptideHit other_casted = other
        cdef PeptideHit self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class PeptideHit_AnalysisResult:
    """
    Cython implementation of _PeptideHit_AnalysisResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit_AnalysisResult.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property score_type:
        def __set__(self,  score_type):
        
            self.inst.get().score_type = deref((convString(score_type)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().score_type
            py_result = convOutputString(_r)
            return py_result
    
    property higher_is_better:
        def __set__(self, bool higher_is_better):
        
            self.inst.get().higher_is_better = (<bool>higher_is_better)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().higher_is_better
            py_result = <bool>_r
            return py_result
    
    property main_score:
        def __set__(self, double main_score):
        
            self.inst.get().main_score = (<double>main_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().main_score
            py_result = <double>_r
            return py_result
    
    property sub_scores:
        def __set__(self, dict sub_scores):
            cdef libcpp_map[_String, double] * v0 = new libcpp_map[_String, double]()
            for key, value in sub_scores.items():
            
            
                deref(v0)[ deref(<_String *> (<String> key).inst.get()) ] = (<double>value)
            
            
            self.inst.get().sub_scores = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().sub_scores
            py_result = dict()
            cdef libcpp_map[_String, double].iterator it__r = _r.begin()
            cdef String itemk_py_result
            while it__r != _r.end():
               #py_result[deref(<_String *> (<String> key).inst.get())] = <double>(deref(it__r).second)
               itemk_py_result = String.__new__(String)
               itemk_py_result.inst = shared_ptr[_String](new _String((deref(it__r)).first))
               # py_result[deref(<_String *> (<String> key).inst.get())] = <double>(deref(it__r).second)
               py_result[itemk_py_result] = <double>(deref(it__r).second)
               inc(it__r)
            return py_result
    
    def __copy__(self):
       cdef PeptideHit_AnalysisResult rv = PeptideHit_AnalysisResult.__new__(PeptideHit_AnalysisResult)
       rv.inst = shared_ptr[_PeptideHit_AnalysisResult](new _PeptideHit_AnalysisResult(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideHit_AnalysisResult rv = PeptideHit_AnalysisResult.__new__(PeptideHit_AnalysisResult)
       rv.inst = shared_ptr[_PeptideHit_AnalysisResult](new _PeptideHit_AnalysisResult(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideHit_AnalysisResult](new _PeptideHit_AnalysisResult())
    
    def _init_1(self, PeptideHit_AnalysisResult in_0 ):
        """
        _init_1(self, in_0: PeptideHit_AnalysisResult ) -> None
        """
        assert isinstance(in_0, PeptideHit_AnalysisResult), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideHit_AnalysisResult](new _PeptideHit_AnalysisResult((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideHit_AnalysisResult ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideHit_AnalysisResult)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeptideHit_PeakAnnotation:
    """
    Cython implementation of _PeptideHit_PeakAnnotation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit_PeakAnnotation.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property annotation:
        def __set__(self,  annotation):
        
            self.inst.get().annotation = deref((convString(annotation)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().annotation
            py_result = convOutputString(_r)
            return py_result
    
    property charge:
        def __set__(self,  charge):
        
            self.inst.get().charge = (<int>charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().charge
            py_result = <int>_r
            return py_result
    
    property mz:
        def __set__(self, double mz):
        
            self.inst.get().mz = (<double>mz)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz
            py_result = <double>_r
            return py_result
    
    property intensity:
        def __set__(self, double intensity):
        
            self.inst.get().intensity = (<double>intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef PeptideHit_PeakAnnotation rv = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
       rv.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideHit_PeakAnnotation rv = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
       rv.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation())
    
    def _init_1(self, PeptideHit_PeakAnnotation in_0 ):
        """
        _init_1(self, in_0: PeptideHit_PeakAnnotation ) -> None
        """
        assert isinstance(in_0, PeptideHit_PeakAnnotation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideHit_PeakAnnotation ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideHit_PeakAnnotation)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def writePeakAnnotationsString_(self,  annotation_string , list annotations ):
        """
        writePeakAnnotationsString_(self, annotation_string: String , annotations: List[PeptideHit_PeakAnnotation] ) -> None
        """
        assert isinstance(annotation_string, String), 'arg annotation_string wrong type'
        assert isinstance(annotations, list) and all(isinstance(elemt_rec, PeptideHit_PeakAnnotation) for elemt_rec in annotations), 'arg annotations wrong type'
    
        cdef libcpp_vector[_PeptideHit_PeakAnnotation] * v1 = new libcpp_vector[_PeptideHit_PeakAnnotation]()
        cdef PeptideHit_PeakAnnotation item1
        for item1 in annotations:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().writePeakAnnotationsString_(deref((<String>annotation_string).inst.get()), deref(v1))
        del v1
    
    def __richcmp__(self, other, op):
        if op not in (2, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, PeptideHit_PeakAnnotation):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef PeptideHit_PeakAnnotation other_casted = other
        cdef PeptideHit_PeakAnnotation self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class RANSAC:
    """
    Cython implementation of _RANSAC[_RansacModelLinear]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSAC[_RansacModelLinear].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RANSAC[_RansacModelLinear]](new _RANSAC[_RansacModelLinear]())
    
    def _init_1(self,  seed ):
        """
        _init_1(self, seed: int ) -> None
        """
        assert isinstance(seed, int) and seed >= 0, 'arg seed wrong type'
    
        self.inst = shared_ptr[_RANSAC[_RansacModelLinear]](new _RANSAC[_RansacModelLinear]((<uint64_t>seed)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, seed: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def ransac(self, list pairs ,  n ,  k , double t ,  d , bool relative_d ):
        """
        ransac(self, pairs: List[List[float, float]] , n: int , k: int , t: float , d: int , relative_d: bool ) -> List[List[float, float]]
        """
        assert isinstance(pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in pairs), 'arg pairs wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
        assert isinstance(k, int) and k >= 0, 'arg k wrong type'
        assert isinstance(t, float), 'arg t wrong type'
        assert isinstance(d, int) and d >= 0, 'arg d wrong type'
        assert isinstance(relative_d, pybool_t), 'arg relative_d wrong type'
        cdef libcpp_vector[libcpp_pair[double,double]] v0 = pairs
    
    
    
    
    
        _r = self.inst.get().ransac(v0, (<size_t>n), (<size_t>k), (<double>t), (<size_t>d), (<bool>relative_d))
        
        cdef list py_result = _r
        return py_result 

cdef class RANSACParam:
    """
    Cython implementation of _RANSACParam

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSACParam.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property n:
        def __set__(self,  n):
        
            self.inst.get().n = (<size_t>n)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().n
            py_result = <size_t>_r
            return py_result
    
    property k:
        def __set__(self,  k):
        
            self.inst.get().k = (<size_t>k)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().k
            py_result = <size_t>_r
            return py_result
    
    property t:
        def __set__(self, double t):
        
            self.inst.get().t = (<double>t)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().t
            py_result = <double>_r
            return py_result
    
    property d:
        def __set__(self,  d):
        
            self.inst.get().d = (<size_t>d)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().d
            py_result = <size_t>_r
            return py_result
    
    property relative_d:
        def __set__(self, bool relative_d):
        
            self.inst.get().relative_d = (<bool>relative_d)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().relative_d
            py_result = <bool>_r
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A simple struct to carry all the parameters required for a RANSAC run
        """
        self.inst = shared_ptr[_RANSACParam](new _RANSACParam())
    
    def _init_1(self,  p_n ,  p_k , double p_t ,  p_d , bool p_relative_d ):
        """
        _init_1(self, p_n: int , p_k: int , p_t: float , p_d: int , p_relative_d: bool ) -> None
        """
        assert isinstance(p_n, int) and p_n >= 0, 'arg p_n wrong type'
        assert isinstance(p_k, int) and p_k >= 0, 'arg p_k wrong type'
        assert isinstance(p_t, float), 'arg p_t wrong type'
        assert isinstance(p_d, int) and p_d >= 0, 'arg p_d wrong type'
        assert isinstance(p_relative_d, pybool_t), 'arg p_relative_d wrong type'
    
    
    
    
    
        self.inst = shared_ptr[_RANSACParam](new _RANSACParam((<size_t>p_n), (<size_t>p_k), (<double>p_t), (<size_t>p_d), (<bool>p_relative_d)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A simple struct to carry all the parameters required for a RANSAC run

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, p_n: int , p_k: int , p_t: float , p_d: int , p_relative_d: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==5) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], float)) and (isinstance(args[3], int) and args[3] >= 0) and (isinstance(args[4], pybool_t)):
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
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result 

cdef class RANSACQuadratic:
    """
    Cython implementation of _RANSAC[_RansacModelQuadratic]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSAC[_RansacModelQuadratic].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RANSAC[_RansacModelQuadratic]](new _RANSAC[_RansacModelQuadratic]())
    
    def _init_1(self,  seed ):
        """
        _init_1(self, seed: int ) -> None
        """
        assert isinstance(seed, int) and seed >= 0, 'arg seed wrong type'
    
        self.inst = shared_ptr[_RANSAC[_RansacModelQuadratic]](new _RANSAC[_RansacModelQuadratic]((<uint64_t>seed)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, seed: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def ransac(self, list pairs ,  n ,  k , double t ,  d , bool relative_d ):
        """
        ransac(self, pairs: List[List[float, float]] , n: int , k: int , t: float , d: int , relative_d: bool ) -> List[List[float, float]]
        """
        assert isinstance(pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in pairs), 'arg pairs wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
        assert isinstance(k, int) and k >= 0, 'arg k wrong type'
        assert isinstance(t, float), 'arg t wrong type'
        assert isinstance(d, int) and d >= 0, 'arg d wrong type'
        assert isinstance(relative_d, pybool_t), 'arg relative_d wrong type'
        cdef libcpp_vector[libcpp_pair[double,double]] v0 = pairs
    
    
    
    
    
        _r = self.inst.get().ransac(v0, (<size_t>n), (<size_t>k), (<double>t), (<size_t>d), (<bool>relative_d))
        
        cdef list py_result = _r
        return py_result 

cdef class Software:
    """
    Cython implementation of _Software

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Software.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Software rv = Software.__new__(Software)
       rv.inst = shared_ptr[_Software](new _Software(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Software rv = Software.__new__(Software)
       rv.inst = shared_ptr[_Software](new _Software(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Software](new _Software())
    
    def _init_1(self, Software in_0 ):
        """
        _init_1(self, in_0: Software ) -> None
        """
        assert isinstance(in_0, Software), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Software](new _Software((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Software ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Software)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the software
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Returns the software version
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  in_0 ):
        """
        setName(self, in_0: Union[bytes, str, String] ) -> None
        Sets the name of the software
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setName(deref((convString(in_0)).get()))
    
    def setVersion(self,  in_0 ):
        """
        setVersion(self, in_0: Union[bytes, str, String] ) -> None
        Sets the software version
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setVersion(deref((convString(in_0)).get())) 

cdef class TransformationXMLFile:
    """
    Cython implementation of _TransformationXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_TransformationXMLFile](new _TransformationXMLFile())
    
    def load(self,  in_0 , TransformationDescription in_1 , bool fit_model ):
        """
        load(self, in_0: Union[bytes, str, String] , in_1: TransformationDescription , fit_model: bool ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
        assert isinstance(fit_model, pybool_t), 'arg fit_model wrong type'
    
    
    
        self.inst.get().load(deref((convString(in_0)).get()), (deref(in_1.inst.get())), (<bool>fit_model))
    
    def store(self,  in_0 , TransformationDescription in_1 ):
        """
        store(self, in_0: Union[bytes, str, String] , in_1: TransformationDescription ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
    
    
        self.inst.get().store(deref((convString(in_0)).get()), (deref(in_1.inst.get()))) 
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
