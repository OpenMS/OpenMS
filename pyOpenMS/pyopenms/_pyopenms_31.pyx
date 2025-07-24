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
def __static_NASequence_fromString( s ):
    """
    __static_NASequence_fromString(s: Union[bytes, str, String] ) -> NASequence
    """
    assert (isinstance(s, str) or isinstance(s, bytes) or isinstance(s, String)), 'arg s wrong type'

    cdef _NASequence * _r = new _NASequence(_fromString_NASequence(deref((convString(s)).get())))
    cdef NASequence py_result = NASequence.__new__(NASequence)
    py_result.inst = shared_ptr[_NASequence](_r)
    return py_result

def __static_PeptideProteinResolution_run(list proteins , PeptideIdentificationList peptides ):
    """
    __static_PeptideProteinResolution_run(proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList ) -> None
    """
    assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
    assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
    cdef ProteinIdentification item0
    for item0 in proteins:
        v0.push_back(deref(item0.inst.get()))

    _run_PeptideProteinResolution(deref(v0), (deref(peptides.inst.get())))
    cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v0.begin()
    replace_0 = []
    while it_proteins != v0.end():
        item0 = ProteinIdentification.__new__(ProteinIdentification)
        item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
        replace_0.append(item0)
        inc(it_proteins)
    proteins[:] = replace_0
    del v0 

cdef class __NASFragmentType:
    None
    Full = 0
    Internal = 1
    FivePrime = 2
    ThreePrime = 3
    AIon = 4
    BIon = 5
    CIon = 6
    XIon = 7
    YIon = 8
    ZIon = 9
    Precursor = 10
    BIonMinusH20 = 11
    YIonMinusH20 = 12
    BIonMinusNH3 = 13
    YIonMinusNH3 = 14
    NonIdentified = 15
    Unannotated = 16
    WIon = 17
    AminusB = 18
    DIon = 19
    SizeOfNASFragmentType = 20

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __PercolatorOutfile_ScoreType:
    None
    QVALUE = 0
    POSTERRPROB = 1
    SCORE = 2
    SIZE_OF_SCORETYPE = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class CVMappings:
    """
    Cython implementation of _CVMappings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappings.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVMappings rv = CVMappings.__new__(CVMappings)
       rv.inst = shared_ptr[_CVMappings](new _CVMappings(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVMappings rv = CVMappings.__new__(CVMappings)
       rv.inst = shared_ptr[_CVMappings](new _CVMappings(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVMappings](new _CVMappings())
    
    def _init_1(self, CVMappings in_0 ):
        """
        _init_1(self, in_0: CVMappings ) -> None
        """
        assert isinstance(in_0, CVMappings), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVMappings](new _CVMappings((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVMappings ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVMappings)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMappingRules(self, list cv_mapping_rules ):
        """
        setMappingRules(self, cv_mapping_rules: List[CVMappingRule] ) -> None
        Sets the mapping rules of the mapping file
        """
        assert isinstance(cv_mapping_rules, list) and all(isinstance(elemt_rec, CVMappingRule) for elemt_rec in cv_mapping_rules), 'arg cv_mapping_rules wrong type'
        cdef libcpp_vector[_CVMappingRule] * v0 = new libcpp_vector[_CVMappingRule]()
        cdef CVMappingRule item0
        for item0 in cv_mapping_rules:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setMappingRules(deref(v0))
        cdef libcpp_vector[_CVMappingRule].iterator it_cv_mapping_rules = v0.begin()
        replace_0 = []
        while it_cv_mapping_rules != v0.end():
            item0 = CVMappingRule.__new__(CVMappingRule)
            item0.inst = shared_ptr[_CVMappingRule](new _CVMappingRule(deref(it_cv_mapping_rules)))
            replace_0.append(item0)
            inc(it_cv_mapping_rules)
        cv_mapping_rules[:] = replace_0
        del v0
    
    def getMappingRules(self):
        """
        getMappingRules(self) -> List[CVMappingRule]
        Returns the mapping rules
        """
        _r = self.inst.get().getMappingRules()
        py_result = []
        cdef libcpp_vector[_CVMappingRule].iterator it__r = _r.begin()
        cdef CVMappingRule item_py_result
        while it__r != _r.end():
           item_py_result = CVMappingRule.__new__(CVMappingRule)
           item_py_result.inst = shared_ptr[_CVMappingRule](new _CVMappingRule(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addMappingRule(self, CVMappingRule cv_mapping_rule ):
        """
        addMappingRule(self, cv_mapping_rule: CVMappingRule ) -> None
        Adds a mapping rule
        """
        assert isinstance(cv_mapping_rule, CVMappingRule), 'arg cv_mapping_rule wrong type'
    
        self.inst.get().addMappingRule((deref(cv_mapping_rule.inst.get())))
    
    def setCVReferences(self, list cv_references ):
        """
        setCVReferences(self, cv_references: List[CVReference] ) -> None
        Sets the CV references
        """
        assert isinstance(cv_references, list) and all(isinstance(elemt_rec, CVReference) for elemt_rec in cv_references), 'arg cv_references wrong type'
        cdef libcpp_vector[_CVReference] * v0 = new libcpp_vector[_CVReference]()
        cdef CVReference item0
        for item0 in cv_references:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setCVReferences(deref(v0))
        cdef libcpp_vector[_CVReference].iterator it_cv_references = v0.begin()
        replace_0 = []
        while it_cv_references != v0.end():
            item0 = CVReference.__new__(CVReference)
            item0.inst = shared_ptr[_CVReference](new _CVReference(deref(it_cv_references)))
            replace_0.append(item0)
            inc(it_cv_references)
        cv_references[:] = replace_0
        del v0
    
    def getCVReferences(self):
        """
        getCVReferences(self) -> List[CVReference]
        Returns the CV references
        """
        _r = self.inst.get().getCVReferences()
        py_result = []
        cdef libcpp_vector[_CVReference].iterator it__r = _r.begin()
        cdef CVReference item_py_result
        while it__r != _r.end():
           item_py_result = CVReference.__new__(CVReference)
           item_py_result.inst = shared_ptr[_CVReference](new _CVReference(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addCVReference(self, CVReference cv_reference ):
        """
        addCVReference(self, cv_reference: CVReference ) -> None
        Adds a CV reference
        """
        assert isinstance(cv_reference, CVReference), 'arg cv_reference wrong type'
    
        self.inst.get().addCVReference((deref(cv_reference.inst.get())))
    
    def hasCVReference(self,  identifier ):
        """
        hasCVReference(self, identifier: Union[bytes, str, String] ) -> bool
        Returns true if a CV reference is given
        """
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
    
        cdef bool _r = self.inst.get().hasCVReference(deref((convString(identifier)).get()))
        py_result = <bool>_r
        return py_result 

cdef class ChargePair:
    """
    Cython implementation of _ChargePair

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChargePair.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChargePair rv = ChargePair.__new__(ChargePair)
       rv.inst = shared_ptr[_ChargePair](new _ChargePair(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChargePair rv = ChargePair.__new__(ChargePair)
       rv.inst = shared_ptr[_ChargePair](new _ChargePair(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChargePair](new _ChargePair())
    
    def _init_1(self, ChargePair in_0 ):
        """
        _init_1(self, in_0: ChargePair ) -> None
        """
        assert isinstance(in_0, ChargePair), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChargePair](new _ChargePair((deref(in_0.inst.get()))))
    
    def _init_2(self,  index0 ,  index1 ,  charge0 ,  charge1 , Compomer compomer , double mass_diff , bool active ):
        """
        _init_2(self, index0: int , index1: int , charge0: int , charge1: int , compomer: Compomer , mass_diff: float , active: bool ) -> None
        """
        assert isinstance(index0, int) and index0 >= 0, 'arg index0 wrong type'
        assert isinstance(index1, int) and index1 >= 0, 'arg index1 wrong type'
        assert isinstance(charge0, int), 'arg charge0 wrong type'
        assert isinstance(charge1, int), 'arg charge1 wrong type'
        assert isinstance(compomer, Compomer), 'arg compomer wrong type'
        assert isinstance(mass_diff, float), 'arg mass_diff wrong type'
        assert isinstance(active, pybool_t), 'arg active wrong type'
    
    
    
    
    
    
    
        self.inst = shared_ptr[_ChargePair](new _ChargePair((<size_t>index0), (<size_t>index1), (<int>charge0), (<int>charge1), (deref(compomer.inst.get())), (<double>mass_diff), (<bool>active)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChargePair ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, index0: int , index1: int , charge0: int , charge1: int , compomer: Compomer , mass_diff: float , active: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChargePair)):
             self._init_1(*args)
        elif (len(args)==7) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int)) and (isinstance(args[3], int)) and (isinstance(args[4], Compomer)) and (isinstance(args[5], float)) and (isinstance(args[6], pybool_t)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCharge(self,  pairID ):
        """
        getCharge(self, pairID: int ) -> int
        Returns the charge (for element 0 or 1)
        """
        assert isinstance(pairID, int), 'arg pairID wrong type'
    
        cdef int _r = self.inst.get().getCharge((<unsigned int>pairID))
        py_result = <int>_r
        return py_result
    
    def setCharge(self,  pairID ,  e ):
        """
        setCharge(self, pairID: int , e: int ) -> None
        Sets the charge (for element 0 or 1)
        """
        assert isinstance(pairID, int), 'arg pairID wrong type'
        assert isinstance(e, int), 'arg e wrong type'
    
    
        self.inst.get().setCharge((<unsigned int>pairID), (<int>e))
    
    def getElementIndex(self,  pairID ):
        """
        getElementIndex(self, pairID: int ) -> int
        Returns the element index (for element 0 or 1)
        """
        assert isinstance(pairID, int), 'arg pairID wrong type'
    
        cdef size_t _r = self.inst.get().getElementIndex((<unsigned int>pairID))
        py_result = <size_t>_r
        return py_result
    
    def setElementIndex(self,  pairID ,  e ):
        """
        setElementIndex(self, pairID: int , e: int ) -> None
        Sets the element index (for element 0 or 1)
        """
        assert isinstance(pairID, int), 'arg pairID wrong type'
        assert isinstance(e, int) and e >= 0, 'arg e wrong type'
    
    
        self.inst.get().setElementIndex((<unsigned int>pairID), (<size_t>e))
    
    def getCompomer(self):
        """
        getCompomer(self) -> Compomer
        Returns the Id of the compomer that explains the mass difference
        """
        cdef _Compomer * _r = new _Compomer(self.inst.get().getCompomer())
        cdef Compomer py_result = Compomer.__new__(Compomer)
        py_result.inst = shared_ptr[_Compomer](_r)
        return py_result
    
    def setCompomer(self, Compomer compomer ):
        """
        setCompomer(self, compomer: Compomer ) -> None
        Sets the compomer id
        """
        assert isinstance(compomer, Compomer), 'arg compomer wrong type'
    
        self.inst.get().setCompomer((deref(compomer.inst.get())))
    
    def getMassDiff(self):
        """
        getMassDiff(self) -> float
        Returns the mass difference
        """
        cdef double _r = self.inst.get().getMassDiff()
        py_result = <double>_r
        return py_result
    
    def setMassDiff(self, double mass_diff ):
        """
        setMassDiff(self, mass_diff: float ) -> None
        Sets the mass difference
        """
        assert isinstance(mass_diff, float), 'arg mass_diff wrong type'
    
        self.inst.get().setMassDiff((<double>mass_diff))
    
    def getEdgeScore(self):
        """
        getEdgeScore(self) -> float
        Returns the ILP edge score
        """
        cdef double _r = self.inst.get().getEdgeScore()
        py_result = <double>_r
        return py_result
    
    def setEdgeScore(self, double score ):
        """
        setEdgeScore(self, score: float ) -> None
        Sets the ILP edge score
        """
        assert isinstance(score, float), 'arg score wrong type'
    
        self.inst.get().setEdgeScore((<double>score))
    
    def isActive(self):
        """
        isActive(self) -> bool
        Is this pair realized?
        """
        cdef bool _r = self.inst.get().isActive()
        py_result = <bool>_r
        return py_result
    
    def setActive(self, bool active ):
        """
        setActive(self, active: bool ) -> None
        """
        assert isinstance(active, pybool_t), 'arg active wrong type'
    
        self.inst.get().setActive((<bool>active)) 

cdef class Date:
    """
    Cython implementation of _Date

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Date.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Date rv = Date.__new__(Date)
       rv.inst = shared_ptr[_Date](new _Date(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Date rv = Date.__new__(Date)
       rv.inst = shared_ptr[_Date](new _Date(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Date](new _Date())
    
    def _init_1(self, Date in_0 ):
        """
        _init_1(self, in_0: Date ) -> None
        """
        assert isinstance(in_0, Date), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Date](new _Date((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Date ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Date)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def set(self,  date ):
        """
        set(self, date: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(date, str) or isinstance(date, bytes) or isinstance(date, String)), 'arg date wrong type'
    
        self.inst.get().set(deref((convString(date)).get()))
    
    def today(self):
        """
        today(self) -> Date
        """
        cdef _Date * _r = new _Date(self.inst.get().today())
        cdef Date py_result = Date.__new__(Date)
        py_result.inst = shared_ptr[_Date](_r)
        return py_result
    
    def get(self):
        """
        get(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().get()
        py_result = convOutputString(_r)
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear() 

cdef class Feature:
    """
    Cython implementation of _Feature

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Feature.html>`_
      -- Inherits from ['UniqueIdInterface', 'RichPeak2D']

    An LC-MS feature
    
    The Feature class is used to describe the two-dimensional signal caused by an
    analyte. It can store a charge state and a list of peptide identifications
    (for peptides). The area occupied by the Feature in the LC-MS data set is
    represented by a list of convex hulls (one for each isotopic peak). There is
    also a convex hull for the entire Feature. The model description can store
    the parameters of a two-dimensional theoretical model of the underlying
    signal in LC-MS. Currently, non-peptide compounds are also represented as
    features
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Feature rv = Feature.__new__(Feature)
       rv.inst = shared_ptr[_Feature](new _Feature(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Feature rv = Feature.__new__(Feature)
       rv.inst = shared_ptr[_Feature](new _Feature(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Feature](new _Feature())
    
    def _init_1(self, Feature in_0 ):
        """
        _init_1(self, in_0: Feature ) -> None
        """
        assert isinstance(in_0, Feature), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Feature](new _Feature((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Feature ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Feature)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
        if not isinstance(other, Feature):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Feature other_casted = other
        cdef Feature self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class MRMFeatureQCFile:
    """
    Cython implementation of _MRMFeatureQCFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureQCFile.html>`_

    File adapter for MRMFeatureQC files
    
    Loads and stores .csv or .tsv files describing an MRMFeatureQC
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMFeatureQCFile rv = MRMFeatureQCFile.__new__(MRMFeatureQCFile)
       rv.inst = shared_ptr[_MRMFeatureQCFile](new _MRMFeatureQCFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeatureQCFile rv = MRMFeatureQCFile.__new__(MRMFeatureQCFile)
       rv.inst = shared_ptr[_MRMFeatureQCFile](new _MRMFeatureQCFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeatureQCFile](new _MRMFeatureQCFile())
    
    def _init_1(self, MRMFeatureQCFile in_0 ):
        """
        _init_1(self, in_0: MRMFeatureQCFile ) -> None
        """
        assert isinstance(in_0, MRMFeatureQCFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeatureQCFile](new _MRMFeatureQCFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeatureQCFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeatureQCFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MRMFeatureQC mrmfqc , bool is_component_group ):
        """
        load(self, filename: Union[bytes, str, String] , mrmfqc: MRMFeatureQC , is_component_group: bool ) -> None
        Loads an MRMFeatureQC file
        
        
        :param filename: The path to the input file
        :param mrmfqc: The output class which will contain the criteria
        :param is_component_group: True if the user intends to load ComponentGroupQCs data, false otherwise
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(mrmfqc, MRMFeatureQC), 'arg mrmfqc wrong type'
        assert isinstance(is_component_group, pybool_t), 'arg is_component_group wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(mrmfqc.inst.get())), (<const bool>is_component_group)) 

cdef class MultiplexIsotopicPeakPattern:
    """
    Cython implementation of _MultiplexIsotopicPeakPattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MultiplexIsotopicPeakPattern.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MultiplexIsotopicPeakPattern rv = MultiplexIsotopicPeakPattern.__new__(MultiplexIsotopicPeakPattern)
       rv.inst = shared_ptr[_MultiplexIsotopicPeakPattern](new _MultiplexIsotopicPeakPattern(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MultiplexIsotopicPeakPattern rv = MultiplexIsotopicPeakPattern.__new__(MultiplexIsotopicPeakPattern)
       rv.inst = shared_ptr[_MultiplexIsotopicPeakPattern](new _MultiplexIsotopicPeakPattern(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  c ,  ppp , MultiplexDeltaMasses ms ,  msi ):
        """
        _init_0(self, c: int , ppp: int , ms: MultiplexDeltaMasses , msi: int ) -> None
        """
        assert isinstance(c, int), 'arg c wrong type'
        assert isinstance(ppp, int), 'arg ppp wrong type'
        assert isinstance(ms, MultiplexDeltaMasses), 'arg ms wrong type'
        assert isinstance(msi, int), 'arg msi wrong type'
    
    
    
    
        self.inst = shared_ptr[_MultiplexIsotopicPeakPattern](new _MultiplexIsotopicPeakPattern((<int>c), (<int>ppp), (deref(ms.inst.get())), (<int>msi)))
    
    def _init_1(self, MultiplexIsotopicPeakPattern in_0 ):
        """
        _init_1(self, in_0: MultiplexIsotopicPeakPattern ) -> None
        """
        assert isinstance(in_0, MultiplexIsotopicPeakPattern), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MultiplexIsotopicPeakPattern](new _MultiplexIsotopicPeakPattern((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, c: int , ppp: int , ms: MultiplexDeltaMasses , msi: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MultiplexIsotopicPeakPattern ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], int)) and (isinstance(args[1], int)) and (isinstance(args[2], MultiplexDeltaMasses)) and (isinstance(args[3], int)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MultiplexIsotopicPeakPattern)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns charge
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def getPeaksPerPeptide(self):
        """
        getPeaksPerPeptide(self) -> int
        Returns peaks per peptide
        """
        cdef int _r = self.inst.get().getPeaksPerPeptide()
        py_result = <int>_r
        return py_result
    
    def getMassShifts(self):
        """
        getMassShifts(self) -> MultiplexDeltaMasses
        Returns mass shifts
        """
        cdef _MultiplexDeltaMasses * _r = new _MultiplexDeltaMasses(self.inst.get().getMassShifts())
        cdef MultiplexDeltaMasses py_result = MultiplexDeltaMasses.__new__(MultiplexDeltaMasses)
        py_result.inst = shared_ptr[_MultiplexDeltaMasses](_r)
        return py_result
    
    def getMassShiftIndex(self):
        """
        getMassShiftIndex(self) -> int
        Returns mass shift index
        """
        cdef int _r = self.inst.get().getMassShiftIndex()
        py_result = <int>_r
        return py_result
    
    def getMassShiftCount(self):
        """
        getMassShiftCount(self) -> int
        Returns number of mass shifts i.e. the number of peptides in the multiplet
        """
        cdef unsigned int _r = self.inst.get().getMassShiftCount()
        py_result = <unsigned int>_r
        return py_result
    
    def getMassShiftAt(self,  i ):
        """
        getMassShiftAt(self, i: int ) -> float
        Returns mass shift at position i
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef double _r = self.inst.get().getMassShiftAt((<int>i))
        py_result = <double>_r
        return py_result
    
    def getMZShiftAt(self,  i ):
        """
        getMZShiftAt(self, i: int ) -> float
        Returns m/z shift at position i
        """
        assert isinstance(i, int), 'arg i wrong type'
    
        cdef double _r = self.inst.get().getMZShiftAt((<int>i))
        py_result = <double>_r
        return py_result
    
    def getMZShiftCount(self):
        """
        getMZShiftCount(self) -> int
        Returns number of m/z shifts
        """
        cdef unsigned int _r = self.inst.get().getMZShiftCount()
        py_result = <unsigned int>_r
        return py_result 

cdef class NASequence:
    """
    Cython implementation of _NASequence

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NASequence.html>`_

    Representation of an RNA sequence
    This class represents nucleic acid sequences in OpenMS. An NASequence
    instance primarily contains a sequence of ribonucleotides.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).toString().c_str() )

    
    def __copy__(self):
       cdef NASequence rv = NASequence.__new__(NASequence)
       rv.inst = shared_ptr[_NASequence](new _NASequence(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NASequence rv = NASequence.__new__(NASequence)
       rv.inst = shared_ptr[_NASequence](new _NASequence(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NASequence](new _NASequence())
    
    def _init_1(self, NASequence in_0 ):
        """
        _init_1(self, in_0: NASequence ) -> None
        """
        assert isinstance(in_0, NASequence), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NASequence](new _NASequence((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NASequence ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NASequence)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSequence(self):
        """
        getSequence(self) -> List[Ribonucleotide]
        """
        _r = self.inst.get().getSequence()
        py_result = []
        cdef libcpp_vector[const _Ribonucleotide *].iterator it__r = _r.begin()
        cdef Ribonucleotide item_py_result
        while it__r != _r.end():
           item_py_result = Ribonucleotide.__new__(Ribonucleotide)
           item_py_result.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(deref(it__r))))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def __getitem__(self,  index ):
        """
        __getitem__(self, index: int ) -> Ribonucleotide
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef long _idx = (<size_t>index)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        cdef const _Ribonucleotide * __r = (deref(self.inst.get())[(<size_t>index)])
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    
    def empty(self):
        """
        empty(self) -> bool
        Check if sequence is empty
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def setSequence(self, list seq ):
        """
        setSequence(self, seq: List[Ribonucleotide] ) -> None
        """
        assert isinstance(seq, list) and all(isinstance(elemt_rec, Ribonucleotide) for elemt_rec in seq), 'arg seq wrong type'
        cdef libcpp_vector[const _Ribonucleotide *] * v0 = new libcpp_vector[const _Ribonucleotide *]()
        cdef Ribonucleotide item0
        for item0 in seq:
            v0.push_back((item0.inst.get()))
        self.inst.get().setSequence(deref(v0))
        del v0
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        Returns the peptide as string with modifications embedded in brackets
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def setFivePrimeMod(self, Ribonucleotide modification ):
        """
        setFivePrimeMod(self, modification: Ribonucleotide ) -> None
        Sets the 5' modification
        """
        assert isinstance(modification, Ribonucleotide), 'arg modification wrong type'
    
        self.inst.get().setFivePrimeMod((modification.inst.get()))
    
    def getFivePrimeMod(self):
        """
        getFivePrimeMod(self) -> Ribonucleotide
        Returns the name (ID) of the N-terminal modification, or an empty string if none is set
        """
        cdef const _Ribonucleotide * __r = (self.inst.get().getFivePrimeMod())
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    
    def setThreePrimeMod(self, Ribonucleotide modification ):
        """
        setThreePrimeMod(self, modification: Ribonucleotide ) -> None
        Sets the 3' modification
        """
        assert isinstance(modification, Ribonucleotide), 'arg modification wrong type'
    
        self.inst.get().setThreePrimeMod((modification.inst.get()))
    
    def getThreePrimeMod(self):
        """
        getThreePrimeMod(self) -> Ribonucleotide
        """
        cdef const _Ribonucleotide * __r = (self.inst.get().getThreePrimeMod())
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    
    def get(self,  index ):
        """
        get(self, index: int ) -> Ribonucleotide
        Returns the residue at position index
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef const _Ribonucleotide * __r = (self.inst.get().get((<size_t>index)))
        if __r == NULL:
            return None
        cdef _Ribonucleotide * _r = new _Ribonucleotide(deref(__r))
        cdef Ribonucleotide py_result = Ribonucleotide.__new__(Ribonucleotide)
        py_result.inst = shared_ptr[_Ribonucleotide](_r)
        return py_result
    
    def set(self,  index , Ribonucleotide r ):
        """
        set(self, index: int , r: Ribonucleotide ) -> None
        Sets the residue at position index
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(r, Ribonucleotide), 'arg r wrong type'
    
    
        self.inst.get().set((<size_t>index), (r.inst.get()))
    
    def _getFormula_0(self):
        """
        _getFormula_0(self) -> EmpiricalFormula
        Returns the formula of the peptide
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def _getFormula_1(self, int type_ ,  charge ):
        """
        _getFormula_1(self, type_: int , charge: int ) -> EmpiricalFormula
        """
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula((<_NASFragmentType>type_), (<int>charge)))
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getFormula(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getFormula(self, ) -> EmpiricalFormula
          :noindex:
        
        Returns the formula of the peptide

        
        .. rubric:: Overload:
        .. py:function:: getFormula(self, type_: int , charge: int ) -> EmpiricalFormula
          :noindex:
    
        """
        if not args:
            return self._getFormula_0(*args)
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) and (isinstance(args[1], int)):
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
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef double _r = self.inst.get().getAverageWeight((<_NASFragmentType>type_), (<int>charge))
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
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) and (isinstance(args[1], int)):
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
        assert type_ in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20], 'arg type_ wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
        cdef double _r = self.inst.get().getMonoWeight((<_NASFragmentType>type_), (<int>charge))
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
        elif (len(args)==2) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) and (isinstance(args[1], int)):
            return self._getMonoWeight_1(*args)
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
    
    def getPrefix(self,  length ):
        """
        getPrefix(self, length: int ) -> NASequence
        Returns a peptide sequence of the first index residues
        """
        assert isinstance(length, int) and length >= 0, 'arg length wrong type'
    
        cdef _NASequence * _r = new _NASequence(self.inst.get().getPrefix((<size_t>length)))
        cdef NASequence py_result = NASequence.__new__(NASequence)
        py_result.inst = shared_ptr[_NASequence](_r)
        return py_result
    
    def getSuffix(self,  length ):
        """
        getSuffix(self, length: int ) -> NASequence
        Returns a peptide sequence of the last index residues
        """
        assert isinstance(length, int) and length >= 0, 'arg length wrong type'
    
        cdef _NASequence * _r = new _NASequence(self.inst.get().getSuffix((<size_t>length)))
        cdef NASequence py_result = NASequence.__new__(NASequence)
        py_result.inst = shared_ptr[_NASequence](_r)
        return py_result
    
    def getSubsequence(self,  start ,  length ):
        """
        getSubsequence(self, start: int , length: int ) -> NASequence
        Returns a peptide sequence of number residues, beginning at position index
        """
        assert isinstance(start, int) and start >= 0, 'arg start wrong type'
        assert isinstance(length, int) and length >= 0, 'arg length wrong type'
    
    
        cdef _NASequence * _r = new _NASequence(self.inst.get().getSubsequence((<size_t>start), (<size_t>length)))
        cdef NASequence py_result = NASequence.__new__(NASequence)
        py_result.inst = shared_ptr[_NASequence](_r)
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
        if op not in (2, 3, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, NASequence):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef NASequence other_casted = other
        cdef NASequence self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get())
    



    ## def __iter__(self):
    ##     _r = self.inst.get().getSequence()
    ##     cdef libcpp_vector[const _Ribonucleotide *].iterator it__r = _r.begin()
    ##     cdef Ribonucleotide item_py_result
    ##     while it__r != _r.end():
    ##        item_py_result = Ribonucleotide.__new__(Ribonucleotide)
    ##        item_py_result.inst = shared_ptr[_Ribonucleotide](new _Ribonucleotide(deref(deref(it__r))))
    ##        yield py_result
    ##        inc(it__r)

    def __iter__(self):

        for r in self.getSequence():
          yield r
    NASFragmentType = __NASFragmentType
    fromString = __static_NASequence_fromString 

cdef class PeptideProteinResolution:
    """
    Cython implementation of _PeptideProteinResolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideProteinResolution.html>`_

    Resolves shared peptides based on protein scores
    
    Resolves connected components of the bipartite protein-peptide graph based
    on protein probabilities/scores and adds them as additional protein_groups
    to the protein identification run processed.
    Thereby greedily assigns shared peptides in this component uniquely to the
    proteins of the current @em best @em indistinguishable protein group, until
    every peptide is uniquely assigned. This effectively allows more peptides to
    be used in ProteinQuantifier at the cost of potentially additional noise in
    the peptides quantities.
    In accordance with most state-of-the-art protein inference tools, only the
    best hit (PSM) for a peptide ID is considered.  Probability ties are
    currently resolved by taking the protein with larger number of peptides
    
    The class could provide iterator for ConnectedComponents in the
    future. One could extend the graph to include all PeptideHits (not only the
    best). It becomes a tripartite graph with larger connected components then.
    Maybe extend it to work with MS1 features. Separate resolution and adding
    groups to output
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideProteinResolution rv = PeptideProteinResolution.__new__(PeptideProteinResolution)
       rv.inst = shared_ptr[_PeptideProteinResolution](new _PeptideProteinResolution(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideProteinResolution rv = PeptideProteinResolution.__new__(PeptideProteinResolution)
       rv.inst = shared_ptr[_PeptideProteinResolution](new _PeptideProteinResolution(deref(self.inst.get())))
       return rv
    
    def _init_0(self, bool statistics ):
        """
        _init_0(self, statistics: bool ) -> None
        """
        assert isinstance(statistics, pybool_t), 'arg statistics wrong type'
    
        self.inst = shared_ptr[_PeptideProteinResolution](new _PeptideProteinResolution((<bool>statistics)))
    
    def _init_1(self, PeptideProteinResolution in_0 ):
        """
        _init_1(self, in_0: PeptideProteinResolution ) -> None
        """
        assert isinstance(in_0, PeptideProteinResolution), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideProteinResolution](new _PeptideProteinResolution((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, statistics: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideProteinResolution ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], pybool_t)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideProteinResolution)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def buildGraph(self, ProteinIdentification protein , PeptideIdentificationList peptides ):
        """
        buildGraph(self, protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None
        Initialize and store the graph (= maps), needs sorted groups for
        correct functionality. Therefore sorts the indist. protein groups
        if not skipped
        
        
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: Vector of ProteinIdentifications with links to the proteins
        :param skip_sort: Skips sorting of groups, nothing is modified then
        """
        assert isinstance(protein, ProteinIdentification), 'arg protein wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
    
        self.inst.get().buildGraph((deref(protein.inst.get())), (deref(peptides.inst.get())))
    
    def resolveGraph(self, ProteinIdentification protein , PeptideIdentificationList peptides ):
        """
        resolveGraph(self, protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None
        Applies resolveConnectedComponent to every component of the graph and
        is able to write statistics when specified. Parameters will
        both be mutated in this method
        
        
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: vector of ProteinIdentifications with links to the proteins
        """
        assert isinstance(protein, ProteinIdentification), 'arg protein wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
    
        self.inst.get().resolveGraph((deref(protein.inst.get())), (deref(peptides.inst.get())))
    
    def findConnectedComponent(self,  root_prot_grp ):
        """
        findConnectedComponent(self, root_prot_grp: int ) -> PeptideProteinResolution_ConnectedComponent
        Does a BFS on the two maps (= two parts of the graph; indist. prot. groups
        and peptides), switching from one to the other in each step
        
        
        :param root_prot_grp: Starts the BFS at this protein group index
        :return: Returns a Connected Component as set of group and peptide indices
        """
        assert isinstance(root_prot_grp, int) and root_prot_grp >= 0, 'arg root_prot_grp wrong type'
    
        cdef _PeptideProteinResolution_ConnectedComponent * _r = new _PeptideProteinResolution_ConnectedComponent(self.inst.get().findConnectedComponent((<size_t &>root_prot_grp)))
        cdef PeptideProteinResolution_ConnectedComponent py_result = PeptideProteinResolution_ConnectedComponent.__new__(PeptideProteinResolution_ConnectedComponent)
        py_result.inst = shared_ptr[_PeptideProteinResolution_ConnectedComponent](_r)
        return py_result
    
    def resolveConnectedComponent(self, PeptideProteinResolution_ConnectedComponent conn_comp , ProteinIdentification protein , PeptideIdentificationList peptides ):
        """
        resolveConnectedComponent(self, conn_comp: PeptideProteinResolution_ConnectedComponent , protein: ProteinIdentification , peptides: PeptideIdentificationList ) -> None
        Resolves connected components based on posterior probabilities and adds them
        as additional protein_groups to the output idXML.
        Thereby greedily assigns shared peptides in this component uniquely to
        the proteins of the current BEST INDISTINGUISHABLE protein group,
        ready to be used in ProteinQuantifier then.
        This is achieved by removing all other evidence from the input
        PeptideIDs and iterating until each peptide is uniquely assigned.
        In accordance with Fido only the best hit (PSM) for an ID is considered.
        Probability ties resolved by taking protein with largest number of peptides
        
        
        :param conn_comp: The component to be resolved
        :param protein: ProteinIdentification object storing IDs and groups
        :param peptides: Vector of ProteinIdentifications with links to the proteins
        """
        assert isinstance(conn_comp, PeptideProteinResolution_ConnectedComponent), 'arg conn_comp wrong type'
        assert isinstance(protein, ProteinIdentification), 'arg protein wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
    
    
        self.inst.get().resolveConnectedComponent((deref(conn_comp.inst.get())), (deref(protein.inst.get())), (deref(peptides.inst.get())))
    run = __static_PeptideProteinResolution_run 

cdef class PeptideProteinResolution_ConnectedComponent:
    """
    Cython implementation of _PeptideProteinResolution_ConnectedComponent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideProteinResolution_ConnectedComponent.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property prot_grp_indices:
        def __set__(self, set prot_grp_indices):
            cdef libcpp_set[size_t] v0 = prot_grp_indices
            self.inst.get().prot_grp_indices = v0
        
    
        def __get__(self):
            _r = self.inst.get().prot_grp_indices
            cdef set py_result = _r
            return py_result
    
    property pep_indices:
        def __set__(self, set pep_indices):
            cdef libcpp_set[size_t] v0 = pep_indices
            self.inst.get().pep_indices = v0
        
    
        def __get__(self):
            _r = self.inst.get().pep_indices
            cdef set py_result = _r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_PeptideProteinResolution_ConnectedComponent](new _PeptideProteinResolution_ConnectedComponent()) 

cdef class PercolatorOutfile:
    """
    Cython implementation of _PercolatorOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorOutfile.html>`_

    Class for reading Percolator tab-delimited output files
    
    For PSM-level output, the file extension should be ".psms"
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PercolatorOutfile rv = PercolatorOutfile.__new__(PercolatorOutfile)
       rv.inst = shared_ptr[_PercolatorOutfile](new _PercolatorOutfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PercolatorOutfile rv = PercolatorOutfile.__new__(PercolatorOutfile)
       rv.inst = shared_ptr[_PercolatorOutfile](new _PercolatorOutfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PercolatorOutfile](new _PercolatorOutfile())
    
    def _init_1(self, PercolatorOutfile in_0 ):
        """
        _init_1(self, in_0: PercolatorOutfile ) -> None
        """
        assert isinstance(in_0, PercolatorOutfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PercolatorOutfile](new _PercolatorOutfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PercolatorOutfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PercolatorOutfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getScoreType(self,  score_type_name ):
        """
        getScoreType(self, score_type_name: Union[bytes, str, String] ) -> int
        Returns a score type given its name
        """
        assert (isinstance(score_type_name, str) or isinstance(score_type_name, bytes) or isinstance(score_type_name, String)), 'arg score_type_name wrong type'
    
        cdef _PercolatorOutfile_ScoreType _r = self.inst.get().getScoreType(deref((convString(score_type_name)).get()))
        py_result = <int>_r
        return py_result
    
    def load(self,  filename , ProteinIdentification proteins , PeptideIdentificationList peptides , SpectrumMetaDataLookup lookup , int output_score ):
        """
        load(self, filename: Union[bytes, str, String] , proteins: ProteinIdentification , peptides: PeptideIdentificationList , lookup: SpectrumMetaDataLookup , output_score: int ) -> None
        Loads a Percolator output file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(proteins, ProteinIdentification), 'arg proteins wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(lookup, SpectrumMetaDataLookup), 'arg lookup wrong type'
        assert output_score in [0, 1, 2, 3], 'arg output_score wrong type'
    
    
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(proteins.inst.get())), (deref(peptides.inst.get())), (deref(lookup.inst.get())), (<_PercolatorOutfile_ScoreType>output_score))
    PercolatorOutfile_ScoreType = __PercolatorOutfile_ScoreType 

cdef class QcMLFile:
    """
    Cython implementation of _QcMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QcMLFile.html>`_
      -- Inherits from ['XMLHandler', 'XMLFile', 'ProgressLogger']

    File adapter for QcML files used to load and store QcML files
    
    This Class is supposed to internally collect the data for the qcML File
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_QcMLFile](new _QcMLFile())
    
    def exportIDstats(self,  filename ):
        """
        exportIDstats(self, filename: Union[bytes, str, String] ) -> Union[bytes, str, String]
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        cdef _String _r = self.inst.get().exportIDstats(deref((convString(filename)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def addRunQualityParameter(self,  r , QualityParameter qp ):
        """
        addRunQualityParameter(self, r: Union[bytes, str, String] , qp: QualityParameter ) -> None
        Adds a QualityParameter to run by the name r
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(qp, QualityParameter), 'arg qp wrong type'
    
    
        self.inst.get().addRunQualityParameter(deref((convString(r)).get()), (deref(qp.inst.get())))
    
    def addRunAttachment(self,  r , Attachment at ):
        """
        addRunAttachment(self, r: Union[bytes, str, String] , at: Attachment ) -> None
        Adds a attachment to run by the name r
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(at, Attachment), 'arg at wrong type'
    
    
        self.inst.get().addRunAttachment(deref((convString(r)).get()), (deref(at.inst.get())))
    
    def addSetQualityParameter(self,  r , QualityParameter qp ):
        """
        addSetQualityParameter(self, r: Union[bytes, str, String] , qp: QualityParameter ) -> None
        Adds a QualityParameter to set by the name r
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(qp, QualityParameter), 'arg qp wrong type'
    
    
        self.inst.get().addSetQualityParameter(deref((convString(r)).get()), (deref(qp.inst.get())))
    
    def addSetAttachment(self,  r , Attachment at ):
        """
        addSetAttachment(self, r: Union[bytes, str, String] , at: Attachment ) -> None
        Adds a attachment to set by the name r
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(at, Attachment), 'arg at wrong type'
    
    
        self.inst.get().addSetAttachment(deref((convString(r)).get()), (deref(at.inst.get())))
    
    def _removeAttachment_0(self,  r , list ids ,  at ):
        """
        _removeAttachment_0(self, r: Union[bytes, str, String] , ids: List[bytes] , at: Union[bytes, str, String] ) -> None
        Removes attachments referencing an id given in ids, from run/set r. All attachments if no attachment name is given with at
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
        assert (isinstance(at, str) or isinstance(at, bytes) or isinstance(at, String)), 'arg at wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in ids:
           v1.push_back(_String(<char *>item1))
    
        self.inst.get().removeAttachment(deref((convString(r)).get()), deref(v1), deref((convString(at)).get()))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        ids[:] = replace
        del v1
    
    def _removeAttachment_1(self,  r ,  at ):
        """
        _removeAttachment_1(self, r: Union[bytes, str, String] , at: Union[bytes, str, String] ) -> None
        Removes attachment with cv accession at from run/set r
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert (isinstance(at, str) or isinstance(at, bytes) or isinstance(at, String)), 'arg at wrong type'
    
    
        self.inst.get().removeAttachment(deref((convString(r)).get()), deref((convString(at)).get()))
    
    def removeAttachment(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeAttachment(self, r: Union[bytes, str, String] , ids: List[bytes] , at: Union[bytes, str, String] ) -> None
          :noindex:
        
        Removes attachments referencing an id given in ids, from run/set r. All attachments if no attachment name is given with at

        
        .. rubric:: Overload:
        .. py:function:: removeAttachment(self, r: Union[bytes, str, String] , at: Union[bytes, str, String] ) -> None
          :noindex:
        
        Removes attachment with cv accession at from run/set r
    
        """
        if (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(i, bytes) for i in args[1])) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
            return self._removeAttachment_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._removeAttachment_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def removeAllAttachments(self,  at ):
        """
        removeAllAttachments(self, at: Union[bytes, str, String] ) -> None
        Removes attachment with cv accession at from all runs/sets
        """
        assert (isinstance(at, str) or isinstance(at, bytes) or isinstance(at, String)), 'arg at wrong type'
    
        self.inst.get().removeAllAttachments(deref((convString(at)).get()))
    
    def removeQualityParameter(self,  r , list ids ):
        """
        removeQualityParameter(self, r: Union[bytes, str, String] , ids: List[bytes] ) -> None
        Removes QualityParameter going by one of the ID attributes given in ids
        """
        assert (isinstance(r, str) or isinstance(r, bytes) or isinstance(r, String)), 'arg r wrong type'
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in ids:
           v1.push_back(_String(<char *>item1))
        self.inst.get().removeQualityParameter(deref((convString(r)).get()), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        ids[:] = replace
        del v1
    
    def merge(self, QcMLFile addendum ,  setname ):
        """
        merge(self, addendum: QcMLFile , setname: Union[bytes, str, String] ) -> None
        Merges the given QCFile into this one
        """
        assert isinstance(addendum, QcMLFile), 'arg addendum wrong type'
        assert (isinstance(setname, str) or isinstance(setname, bytes) or isinstance(setname, String)), 'arg setname wrong type'
    
    
        self.inst.get().merge((deref(addendum.inst.get())), deref((convString(setname)).get()))
    
    def collectSetParameter(self,  setname ,  qp , list ret ):
        """
        collectSetParameter(self, setname: Union[bytes, str, String] , qp: Union[bytes, str, String] , ret: List[bytes] ) -> None
        Collects the values of given QPs (as CVid) of the given set
        """
        assert (isinstance(setname, str) or isinstance(setname, bytes) or isinstance(setname, String)), 'arg setname wrong type'
        assert (isinstance(qp, str) or isinstance(qp, bytes) or isinstance(qp, String)), 'arg qp wrong type'
        assert isinstance(ret, list) and all(isinstance(i, bytes) for i in ret), 'arg ret wrong type'
    
    
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in ret:
           v2.push_back(_String(<char *>item2))
        self.inst.get().collectSetParameter(deref((convString(setname)).get()), deref((convString(qp)).get()), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        ret[:] = replace
        del v2
    
    def exportAttachment(self,  filename ,  qpname ):
        """
        exportAttachment(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Returns a String of a tab separated rows if found empty string else from run/set by the name filename of the qualityparameter by the name qpname
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(qpname, str) or isinstance(qpname, bytes) or isinstance(qpname, String)), 'arg qpname wrong type'
    
    
        cdef _String _r = self.inst.get().exportAttachment(deref((convString(filename)).get()), deref((convString(qpname)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def getRunNames(self, list ids ):
        """
        getRunNames(self, ids: List[bytes] ) -> None
        Gives the names of the registered runs in the vector ids
        """
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in ids:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getRunNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        ids[:] = replace
        del v0
    
    def existsRun(self,  filename ):
        """
        existsRun(self, filename: Union[bytes, str, String] ) -> bool
        Returns true if the given run id is present in this file, if checkname is true it also checks the names
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        cdef bool _r = self.inst.get().existsRun(deref((convString(filename)).get()))
        py_result = <bool>_r
        return py_result
    
    def existsSet(self,  filename ):
        """
        existsSet(self, filename: Union[bytes, str, String] ) -> bool
        Returns true if the given set id is present in this file, if checkname is true it also checks the names
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        cdef bool _r = self.inst.get().existsSet(deref((convString(filename)).get()))
        py_result = <bool>_r
        return py_result
    
    def existsRunQualityParameter(self,  filename ,  qpname , list ids ):
        """
        existsRunQualityParameter(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] , ids: List[bytes] ) -> None
        Returns the ids of the parameter name given if found in given run empty else
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(qpname, str) or isinstance(qpname, bytes) or isinstance(qpname, String)), 'arg qpname wrong type'
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
    
    
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in ids:
           v2.push_back(_String(<char *>item2))
        self.inst.get().existsRunQualityParameter(deref((convString(filename)).get()), deref((convString(qpname)).get()), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        ids[:] = replace
        del v2
    
    def existsSetQualityParameter(self,  filename ,  qpname , list ids ):
        """
        existsSetQualityParameter(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] , ids: List[bytes] ) -> None
        Returns the ids of the parameter name given if found in given set, empty else
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(qpname, str) or isinstance(qpname, bytes) or isinstance(qpname, String)), 'arg qpname wrong type'
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
    
    
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in ids:
           v2.push_back(_String(<char *>item2))
        self.inst.get().existsSetQualityParameter(deref((convString(filename)).get()), deref((convString(qpname)).get()), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        ids[:] = replace
        del v2
    
    def store(self,  filename ):
        """
        store(self, filename: Union[bytes, str, String] ) -> None
        Store the qcML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().store(deref((convString(filename)).get()))
    
    def load(self,  filename ):
        """
        load(self, filename: Union[bytes, str, String] ) -> None
        Load a QCFile
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().load(deref((convString(filename)).get()))
    
    def registerRun(self,  id_ ,  name ):
        """
        registerRun(self, id_: Union[bytes, str, String] , name: Union[bytes, str, String] ) -> None
        Registers a run in the qcml file with the respective mappings
        """
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
    
        self.inst.get().registerRun(deref((convString(id_)).get()), deref((convString(name)).get()))
    
    def registerSet(self,  id_ ,  name , set names ):
        """
        registerSet(self, id_: Union[bytes, str, String] , name: Union[bytes, str, String] , names: Set[bytes] ) -> None
        Registers a set in the qcml file with the respective mappings
        """
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(names, set) and all(isinstance(i, bytes) for i in names), 'arg names wrong type'
    
    
        cdef libcpp_set[_String] * v2 = new libcpp_set[_String]()
        cdef bytes item2
        for item2 in names:
           v2.insert(_String(<char *>item2))
        self.inst.get().registerSet(deref((convString(id_)).get()), deref((convString(name)).get()), deref(v2))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v2.begin()
        while it != v2.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        names.clear()
        names.update(replace)
        del v2
    
    def exportQP(self,  filename ,  qpname ):
        """
        exportQP(self, filename: Union[bytes, str, String] , qpname: Union[bytes, str, String] ) -> Union[bytes, str, String]
        Returns a String value in quotation of a QualityParameter by the name qpname in run/set by the name filename
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(qpname, str) or isinstance(qpname, bytes) or isinstance(qpname, String)), 'arg qpname wrong type'
    
    
        cdef _String _r = self.inst.get().exportQP(deref((convString(filename)).get()), deref((convString(qpname)).get()))
        py_result = convOutputString(_r)
        return py_result
    
    def exportQPs(self,  filename , list qpnames ):
        """
        exportQPs(self, filename: Union[bytes, str, String] , qpnames: List[bytes] ) -> Union[bytes, str, String]
        Returns a String of a tab separated QualityParameter by the name qpname in run/set by the name filename
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(qpnames, list) and all(isinstance(li, bytes) for li in qpnames), 'arg qpnames wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in qpnames:
           v1.push_back(_String(<char *>item1))
        cdef _String _r = self.inst.get().exportQPs(deref((convString(filename)).get()), deref(v1))
        del v1
        py_result = convOutputString(_r)
        return py_result
    
    def getRunIDs(self, list ids ):
        """
        getRunIDs(self, ids: List[bytes] ) -> None
        Gives the ids of the registered runs in the vector ids
        """
        assert isinstance(ids, list) and all(isinstance(i, bytes) for i in ids), 'arg ids wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in ids:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getRunIDs(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        ids[:] = replace
        del v0
    
    def reset(self):
        """
        reset(self) -> None
        """
        self.inst.get().reset()
    
    def error(self, int mode ,  msg ,  line ,  column ):
        """
        error(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None
        """
        assert mode in [0, 1], 'arg mode wrong type'
        assert (isinstance(msg, str) or isinstance(msg, bytes) or isinstance(msg, String)), 'arg msg wrong type'
        assert isinstance(line, int), 'arg line wrong type'
        assert isinstance(column, int), 'arg column wrong type'
    
    
    
    
        self.inst.get().error((<_ActionMode>mode), deref((convString(msg)).get()), (<unsigned int>line), (<unsigned int>column))
    
    def warning(self, int mode ,  msg ,  line ,  column ):
        """
        warning(self, mode: int , msg: Union[bytes, str, String] , line: int , column: int ) -> None
        """
        assert mode in [0, 1], 'arg mode wrong type'
        assert (isinstance(msg, str) or isinstance(msg, bytes) or isinstance(msg, String)), 'arg msg wrong type'
        assert isinstance(line, int), 'arg line wrong type'
        assert isinstance(column, int), 'arg column wrong type'
    
    
    
    
        self.inst.get().warning((<_ActionMode>mode), deref((convString(msg)).get()), (<unsigned int>line), (<unsigned int>column))
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
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
    
    def map2csv(self, dict csv_table, String separator):
        cdef libcpp_map[_String, libcpp_map[_String, _String] ] c_dict_outer
        cdef libcpp_map[_String, _String] c_dict_inner
        for k,v in csv_table.iteritems():
            c_dict_inner.clear()
            for k_i,v_i in v.iteritems():
                c_dict_inner[ _String(<char *>k_i) ] = _String(<char *>v_i)
            c_dict_outer[ _String(<char *>k) ] = c_dict_inner

        
        cdef _String _r = self.inst.get().map2csv(c_dict_outer, deref(separator.inst.get()) )
        py_result = _cast_const_away(<char*>_r.c_str())
        return py_result 

cdef class QualityParameter:
    """
    Cython implementation of _QualityParameter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1QualityParameter.html>`_
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
    
    property flag:
        def __set__(self,  flag):
        
            self.inst.get().flag = deref((convString(flag)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().flag
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef QualityParameter rv = QualityParameter.__new__(QualityParameter)
       rv.inst = shared_ptr[_QualityParameter](new _QualityParameter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef QualityParameter rv = QualityParameter.__new__(QualityParameter)
       rv.inst = shared_ptr[_QualityParameter](new _QualityParameter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_QualityParameter](new _QualityParameter())
    
    def _init_1(self, QualityParameter in_0 ):
        """
        _init_1(self, in_0: QualityParameter ) -> None
        """
        assert isinstance(in_0, QualityParameter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_QualityParameter](new _QualityParameter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: QualityParameter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], QualityParameter)):
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
    
    def __richcmp__(self, other, op):
        if op not in (2, 0, 4):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, QualityParameter):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef QualityParameter other_casted = other
        cdef QualityParameter self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get())
        if op==4:
            return deref(self_casted.inst.get()) > deref(other_casted.inst.get()) 

cdef class SiriusFragmentAnnotation:
    """
    Cython implementation of _SiriusFragmentAnnotation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusFragmentAnnotation.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusFragmentAnnotation rv = SiriusFragmentAnnotation.__new__(SiriusFragmentAnnotation)
       rv.inst = shared_ptr[_SiriusFragmentAnnotation](new _SiriusFragmentAnnotation(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusFragmentAnnotation rv = SiriusFragmentAnnotation.__new__(SiriusFragmentAnnotation)
       rv.inst = shared_ptr[_SiriusFragmentAnnotation](new _SiriusFragmentAnnotation(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusFragmentAnnotation](new _SiriusFragmentAnnotation())
    
    def _init_1(self, SiriusFragmentAnnotation in_0 ):
        """
        _init_1(self, in_0: SiriusFragmentAnnotation ) -> None
        """
        assert isinstance(in_0, SiriusFragmentAnnotation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusFragmentAnnotation](new _SiriusFragmentAnnotation((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusFragmentAnnotation ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusFragmentAnnotation)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractAnnotationsFromSiriusFile(self,  path_to_sirius_workspace ,  max_rank , bool decoy , bool use_exact_mass ):
        """
        extractAnnotationsFromSiriusFile(self, path_to_sirius_workspace: String , max_rank: int , decoy: bool , use_exact_mass: bool ) -> List[MSSpectrum]
        """
        assert isinstance(path_to_sirius_workspace, String), 'arg path_to_sirius_workspace wrong type'
        assert isinstance(max_rank, int) and max_rank >= 0, 'arg max_rank wrong type'
        assert isinstance(decoy, pybool_t), 'arg decoy wrong type'
        assert isinstance(use_exact_mass, pybool_t), 'arg use_exact_mass wrong type'
    
    
    
    
        _r = self.inst.get().extractAnnotationsFromSiriusFile(deref((<String>path_to_sirius_workspace).inst.get()), (<size_t>max_rank), (<bool>decoy), (<bool>use_exact_mass))
        py_result = []
        cdef libcpp_vector[_MSSpectrum].iterator it__r = _r.begin()
        cdef MSSpectrum item_py_result
        while it__r != _r.end():
           item_py_result = MSSpectrum.__new__(MSSpectrum)
           item_py_result.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def extractAndResolveSiriusAnnotations(self, list sirius_workspace_subdirs , double score_threshold , bool use_exact_mass , bool decoy_generation ):
        """
        extractAndResolveSiriusAnnotations(self, sirius_workspace_subdirs: List[bytes] , score_threshold: float , use_exact_mass: bool , decoy_generation: bool ) -> List[SiriusFragmentAnnotation_SiriusTargetDecoySpectra]
        """
        assert isinstance(sirius_workspace_subdirs, list) and all(isinstance(i, bytes) for i in sirius_workspace_subdirs), 'arg sirius_workspace_subdirs wrong type'
        assert isinstance(score_threshold, float), 'arg score_threshold wrong type'
        assert isinstance(use_exact_mass, pybool_t), 'arg use_exact_mass wrong type'
        assert isinstance(decoy_generation, pybool_t), 'arg decoy_generation wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in sirius_workspace_subdirs:
           v0.push_back(_String(<char *>item0))
    
    
    
        _r = self.inst.get().extractAndResolveSiriusAnnotations(deref(v0), (<double>score_threshold), (<bool>use_exact_mass), (<bool>decoy_generation))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        sirius_workspace_subdirs[:] = replace
        del v0
        py_result = []
        cdef libcpp_vector[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra].iterator it__r = _r.begin()
        cdef SiriusFragmentAnnotation_SiriusTargetDecoySpectra item_py_result
        while it__r != _r.end():
           item_py_result = SiriusFragmentAnnotation_SiriusTargetDecoySpectra.__new__(SiriusFragmentAnnotation_SiriusTargetDecoySpectra)
           item_py_result.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def extract_columnname_to_columnindex(self, CsvFile csvfile ):
        """
        extract_columnname_to_columnindex(self, csvfile: CsvFile ) -> Dict[bytes, int]
        """
        assert isinstance(csvfile, CsvFile), 'arg csvfile wrong type'
    
        _r = self.inst.get().extract_columnname_to_columnindex((deref(csvfile.inst.get())))
        py_result = dict()
        cdef libcpp_map[libcpp_string, size_t].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result[<libcpp_string>(deref(it__r).first)] = <size_t>(deref(it__r).second)
           inc(it__r)
        return py_result 

cdef class SiriusFragmentAnnotation_SiriusTargetDecoySpectra:
    """
    Cython implementation of _SiriusFragmentAnnotation_SiriusTargetDecoySpectra

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusFragmentAnnotation_SiriusTargetDecoySpectra.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusFragmentAnnotation_SiriusTargetDecoySpectra rv = SiriusFragmentAnnotation_SiriusTargetDecoySpectra.__new__(SiriusFragmentAnnotation_SiriusTargetDecoySpectra)
       rv.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusFragmentAnnotation_SiriusTargetDecoySpectra rv = SiriusFragmentAnnotation_SiriusTargetDecoySpectra.__new__(SiriusFragmentAnnotation_SiriusTargetDecoySpectra)
       rv.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra())
    
    def _init_1(self, SiriusFragmentAnnotation_SiriusTargetDecoySpectra in_0 ):
        """
        _init_1(self, in_0: SiriusFragmentAnnotation_SiriusTargetDecoySpectra ) -> None
        """
        assert isinstance(in_0, SiriusFragmentAnnotation_SiriusTargetDecoySpectra), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusFragmentAnnotation_SiriusTargetDecoySpectra ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusFragmentAnnotation_SiriusTargetDecoySpectra)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SplinePackage:
    """
    Cython implementation of _SplinePackage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SplinePackage.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SplinePackage rv = SplinePackage.__new__(SplinePackage)
       rv.inst = shared_ptr[_SplinePackage](new _SplinePackage(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SplinePackage rv = SplinePackage.__new__(SplinePackage)
       rv.inst = shared_ptr[_SplinePackage](new _SplinePackage(deref(self.inst.get())))
       return rv
    
    def _init_0(self, list pos , list intensity ):
        """
        _init_0(self, pos: List[float] , intensity: List[float] ) -> None
        """
        assert isinstance(pos, list) and all(isinstance(elemt_rec, float) for elemt_rec in pos), 'arg pos wrong type'
        assert isinstance(intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in intensity), 'arg intensity wrong type'
        cdef libcpp_vector[double] v0 = pos
        cdef libcpp_vector[double] v1 = intensity
        self.inst = shared_ptr[_SplinePackage](new _SplinePackage(v0, v1))
        
        
    
    def _init_1(self, SplinePackage in_0 ):
        """
        _init_1(self, in_0: SplinePackage ) -> None
        """
        assert isinstance(in_0, SplinePackage), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SplinePackage](new _SplinePackage((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, pos: List[float] , intensity: List[float] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SplinePackage ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SplinePackage)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPosMin(self):
        """
        getPosMin(self) -> float
        Returns the minimum position for which the spline fit is valid
        """
        cdef double _r = self.inst.get().getPosMin()
        py_result = <double>_r
        return py_result
    
    def getPosMax(self):
        """
        getPosMax(self) -> float
        Returns the maximum position for which the spline fit is valid
        """
        cdef double _r = self.inst.get().getPosMax()
        py_result = <double>_r
        return py_result
    
    def getPosStepWidth(self):
        """
        getPosStepWidth(self) -> float
        Returns a sensible position step width for the package
        """
        cdef double _r = self.inst.get().getPosStepWidth()
        py_result = <double>_r
        return py_result
    
    def isInPackage(self, double pos ):
        """
        isInPackage(self, pos: float ) -> bool
        Returns true if position in
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef bool _r = self.inst.get().isInPackage((<double>pos))
        py_result = <bool>_r
        return py_result
    
    def eval(self, double pos ):
        """
        eval(self, pos: float ) -> float
        Returns interpolated intensity position `pos`
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().eval((<double>pos))
        py_result = <double>_r
        return py_result 

cdef class _Interfaces_BinaryDataArray:
    """
    Cython implementation of _BinaryDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1BinaryDataArray.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property data:
        def __set__(self, list data):
            cdef libcpp_vector[double] v0 = data
            self.inst.get().data = v0
            
    
        def __get__(self):
            _r = self.inst.get().data
            cdef list py_result = _r
            return py_result
    
    def __copy__(self):
       cdef _Interfaces_BinaryDataArray rv = _Interfaces_BinaryDataArray.__new__(_Interfaces_BinaryDataArray)
       rv.inst = shared_ptr[_BinaryDataArray](new _BinaryDataArray(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef _Interfaces_BinaryDataArray rv = _Interfaces_BinaryDataArray.__new__(_Interfaces_BinaryDataArray)
       rv.inst = shared_ptr[_BinaryDataArray](new _BinaryDataArray(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BinaryDataArray](new _BinaryDataArray())
    
    def _init_1(self, _Interfaces_BinaryDataArray in_0 ):
        """
        _init_1(self, in_0: _Interfaces_BinaryDataArray ) -> None
        """
        assert isinstance(in_0, _Interfaces_BinaryDataArray), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BinaryDataArray](new _BinaryDataArray((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: _Interfaces_BinaryDataArray ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], _Interfaces_BinaryDataArray)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class _Interfaces_Chromatogram:
    """
    Cython implementation of _Chromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1Chromatogram.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef _Interfaces_Chromatogram rv = _Interfaces_Chromatogram.__new__(_Interfaces_Chromatogram)
       rv.inst = shared_ptr[_Chromatogram](new _Chromatogram(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef _Interfaces_Chromatogram rv = _Interfaces_Chromatogram.__new__(_Interfaces_Chromatogram)
       rv.inst = shared_ptr[_Chromatogram](new _Chromatogram(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Chromatogram](new _Chromatogram())
    
    def _init_1(self, _Interfaces_Chromatogram in_0 ):
        """
        _init_1(self, in_0: _Interfaces_Chromatogram ) -> None
        """
        assert isinstance(in_0, _Interfaces_Chromatogram), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Chromatogram](new _Chromatogram((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: _Interfaces_Chromatogram ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], _Interfaces_Chromatogram)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getTimeArray(self):
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getTimeArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def setTimeArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setTimeArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0) 

cdef class _Interfaces_Spectrum:
    """
    Cython implementation of _Spectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Interfaces_1_1Spectrum.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef _Interfaces_Spectrum rv = _Interfaces_Spectrum.__new__(_Interfaces_Spectrum)
       rv.inst = shared_ptr[_Spectrum](new _Spectrum(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef _Interfaces_Spectrum rv = _Interfaces_Spectrum.__new__(_Interfaces_Spectrum)
       rv.inst = shared_ptr[_Spectrum](new _Spectrum(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Spectrum](new _Spectrum())
    
    def _init_1(self, _Interfaces_Spectrum in_0 ):
        """
        _init_1(self, in_0: _Interfaces_Spectrum ) -> None
        """
        assert isinstance(in_0, _Interfaces_Spectrum), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Spectrum](new _Spectrum((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: _Interfaces_Spectrum ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], _Interfaces_Spectrum)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMZArray(self):
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getMZArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def getIntensityArray(self):
        cdef shared_ptr[_BinaryDataArray] _r = self.inst.get().getIntensityArray()
        cdef libcpp_vector[double] _vec = _r.get().data
        cdef list py_result = _vec
        return py_result

    def setMZArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setMZArray(v0)

    def setIntensityArray(self, list data):
        assert isinstance(data, list), 'arg transitions wrong type'

        cdef shared_ptr[_BinaryDataArray] v0 = shared_ptr[_BinaryDataArray](new _BinaryDataArray() ) 
        cdef libcpp_vector[double] _vec = data
        v0.get().data = data
        self.inst.get().setIntensityArray(v0) 
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
