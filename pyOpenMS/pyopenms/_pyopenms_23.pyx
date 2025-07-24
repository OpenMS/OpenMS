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
def __static_CachedmzML_load( filename , CachedmzML exp ):
    """
    __static_CachedmzML_load(filename: Union[bytes, str, String] , exp: CachedmzML ) -> None
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert isinstance(exp, CachedmzML), 'arg exp wrong type'


    _load_CachedmzML(deref((convString(filename)).get()), (deref(exp.inst.get())))

def __static_CachedmzML_store( filename , MSExperiment exp ):
    """
    __static_CachedmzML_store(filename: Union[bytes, str, String] , exp: MSExperiment ) -> None
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'


    _store_CachedmzML(deref((convString(filename)).get()), (deref(exp.inst.get()))) 

cdef class __DerivatizationAgent:
    None
    NOT_SELECTED = 0
    TBDMS = 1
    SIZE_OF_DERIVATIZATIONAGENT = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class MT_QUANTMETHOD:
    None
    MT_QUANT_AREA = 0
    MT_QUANT_MEDIAN = 1
    SIZE_OF_MT_QUANTMETHOD = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __MassIntensityType:
    None
    NORM_MAX = 0
    NORM_SUM = 1
    SIZE_OF_MASSINTENSITYTYPE = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AnnotationStatistics:
    """
    Cython implementation of _AnnotationStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AnnotationStatistics.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property states:
        def __set__(self, list states):
            cdef libcpp_vector[size_t] v0 = states
            self.inst.get().states = v0
            
    
        def __get__(self):
            _r = self.inst.get().states
            cdef list py_result = _r
            return py_result
    
    def __copy__(self):
       cdef AnnotationStatistics rv = AnnotationStatistics.__new__(AnnotationStatistics)
       rv.inst = shared_ptr[_AnnotationStatistics](new _AnnotationStatistics(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AnnotationStatistics rv = AnnotationStatistics.__new__(AnnotationStatistics)
       rv.inst = shared_ptr[_AnnotationStatistics](new _AnnotationStatistics(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AnnotationStatistics](new _AnnotationStatistics())
    
    def _init_1(self, AnnotationStatistics in_0 ):
        """
        _init_1(self, in_0: AnnotationStatistics ) -> None
        """
        assert isinstance(in_0, AnnotationStatistics), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AnnotationStatistics](new _AnnotationStatistics((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AnnotationStatistics ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AnnotationStatistics)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, AnnotationStatistics):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef AnnotationStatistics other_casted = other
        cdef AnnotationStatistics self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class CachedmzML:
    """
    Cython implementation of _CachedmzML

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CachedmzML.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CachedmzML rv = CachedmzML.__new__(CachedmzML)
       rv.inst = shared_ptr[_CachedmzML](new _CachedmzML(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CachedmzML rv = CachedmzML.__new__(CachedmzML)
       rv.inst = shared_ptr[_CachedmzML](new _CachedmzML(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A class that uses on-disk caching to read and write spectra and chromatograms
        """
        self.inst = shared_ptr[_CachedmzML](new _CachedmzML())
    
    def _init_1(self, CachedmzML in_0 ):
        """
        _init_1(self, in_0: CachedmzML ) -> None
        """
        assert isinstance(in_0, CachedmzML), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CachedmzML](new _CachedmzML((deref(in_0.inst.get()))))
    
    def _init_2(self,  filename ):
        """
        _init_2(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_CachedmzML](new _CachedmzML(deref((convString(filename)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A class that uses on-disk caching to read and write spectra and chromatograms

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CachedmzML ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CachedmzML)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getNrSpectra(self):
        """
        getNrSpectra(self) -> int
        """
        cdef size_t _r = self.inst.get().getNrSpectra()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatograms(self):
        """
        getNrChromatograms(self) -> int
        """
        cdef size_t _r = self.inst.get().getNrChromatograms()
        py_result = <size_t>_r
        return py_result
    
    def getSpectrum(self,  idx ):
        """
        getSpectrum(self, idx: int ) -> MSSpectrum
        """
        assert isinstance(idx, int) and idx >= 0, 'arg idx wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((<size_t>idx)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def getChromatogram(self,  idx ):
        """
        getChromatogram(self, idx: int ) -> MSChromatogram
        """
        assert isinstance(idx, int) and idx >= 0, 'arg idx wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getChromatogram((<size_t>idx)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def getMetaData(self):
        """
        getMetaData(self) -> MSExperiment
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getMetaData())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    load = __static_CachedmzML_load
    store = __static_CachedmzML_store 

cdef class DecoyGenerator:
    """
    Cython implementation of _DecoyGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DecoyGenerator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DecoyGenerator rv = DecoyGenerator.__new__(DecoyGenerator)
       rv.inst = shared_ptr[_DecoyGenerator](new _DecoyGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DecoyGenerator rv = DecoyGenerator.__new__(DecoyGenerator)
       rv.inst = shared_ptr[_DecoyGenerator](new _DecoyGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DecoyGenerator](new _DecoyGenerator())
    
    def _init_1(self, DecoyGenerator in_0 ):
        """
        _init_1(self, in_0: DecoyGenerator ) -> None
        """
        assert isinstance(in_0, DecoyGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DecoyGenerator](new _DecoyGenerator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DecoyGenerator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DecoyGenerator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setSeed(self,  in_0 ):
        """
        setSeed(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        self.inst.get().setSeed((<uint64_t>in_0))
    
    def reverseProtein(self, AASequence protein ):
        """
        reverseProtein(self, protein: AASequence ) -> AASequence
        Reverses the protein sequence
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().reverseProtein((deref(protein.inst.get()))))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def reversePeptides(self, AASequence protein ,  protease ):
        """
        reversePeptides(self, protein: AASequence , protease: Union[bytes, str, String] ) -> AASequence
        Reverses the protein's peptide sequences between enzymatic cutting positions
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
        assert (isinstance(protease, str) or isinstance(protease, bytes) or isinstance(protease, String)), 'arg protease wrong type'
    
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().reversePeptides((deref(protein.inst.get())), deref((convString(protease)).get())))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result
    
    def shufflePeptides(self, AASequence aas ,  protease ,  max_attempts ):
        """
        shufflePeptides(self, aas: AASequence , protease: Union[bytes, str, String] , max_attempts: int ) -> AASequence
        Shuffle the protein's peptide sequences between enzymatic cutting positions, each peptide is shuffled @param max_attempts times to minimize sequence identity
        """
        assert isinstance(aas, AASequence), 'arg aas wrong type'
        assert (isinstance(protease, str) or isinstance(protease, bytes) or isinstance(protease, String)), 'arg protease wrong type'
        assert isinstance(max_attempts, int), 'arg max_attempts wrong type'
    
    
    
        cdef _AASequence * _r = new _AASequence(self.inst.get().shufflePeptides((deref(aas.inst.get())), deref((convString(protease)).get()), (<const int>max_attempts)))
        cdef AASequence py_result = AASequence.__new__(AASequence)
        py_result.inst = shared_ptr[_AASequence](_r)
        return py_result 

cdef class EmpiricalFormula:
    """
    Cython implementation of _EmpiricalFormula

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmpiricalFormula.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EmpiricalFormula rv = EmpiricalFormula.__new__(EmpiricalFormula)
       rv.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EmpiricalFormula rv = EmpiricalFormula.__new__(EmpiricalFormula)
       rv.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Representation of an empirical formula
        """
        self.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula())
    
    def _init_1(self, EmpiricalFormula in_0 ):
        """
        _init_1(self, in_0: EmpiricalFormula ) -> None
        """
        assert isinstance(in_0, EmpiricalFormula), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula((deref(in_0.inst.get()))))
    
    def _init_2(self,  in_0 ):
        """
        _init_2(self, in_0: Union[bytes, str, String] ) -> None
        EmpiricalFormula Constructor from string
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref((convString(in_0)).get())))
    
    def _init_3(self,  number , Element element ,  charge ):
        """
        _init_3(self, number: int , element: Element , charge: int ) -> None
        EmpiricalFormula Constructor with element pointer and number
        """
        assert isinstance(number, int) and number >= 0, 'arg number wrong type'
        assert isinstance(element, Element), 'arg element wrong type'
        assert isinstance(charge, int) and charge >= 0, 'arg charge wrong type'
    
    
    
        self.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula((<ptrdiff_t>number), (element.inst.get()), (<ptrdiff_t>charge)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Representation of an empirical formula

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EmpiricalFormula ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Union[bytes, str, String] ) -> None
          :noindex:
        
        EmpiricalFormula Constructor from string

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, number: int , element: Element , charge: int ) -> None
          :noindex:
        
        EmpiricalFormula Constructor with element pointer and number
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EmpiricalFormula)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        elif (len(args)==3) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], Element)) and (isinstance(args[2], int) and args[2] >= 0):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMonoWeight(self):
        """
        getMonoWeight(self) -> float
        Returns the mono isotopic weight of the formula (includes proton charges)
        """
        cdef double _r = self.inst.get().getMonoWeight()
        py_result = <double>_r
        return py_result
    
    def getAverageWeight(self):
        """
        getAverageWeight(self) -> float
        Returns the average weight of the formula (includes proton charges)
        """
        cdef double _r = self.inst.get().getAverageWeight()
        py_result = <double>_r
        return py_result
    
    def estimateFromWeightAndComp(self, double average_weight , double C , double H , double N , double O , double S , double P ):
        """
        estimateFromWeightAndComp(self, average_weight: float , C: float , H: float , N: float , O: float , S: float , P: float ) -> bool
        Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(C, float), 'arg C wrong type'
        assert isinstance(H, float), 'arg H wrong type'
        assert isinstance(N, float), 'arg N wrong type'
        assert isinstance(O, float), 'arg O wrong type'
        assert isinstance(S, float), 'arg S wrong type'
        assert isinstance(P, float), 'arg P wrong type'
    
    
    
    
    
    
    
        cdef bool _r = self.inst.get().estimateFromWeightAndComp((<double>average_weight), (<double>C), (<double>H), (<double>N), (<double>O), (<double>S), (<double>P))
        py_result = <bool>_r
        return py_result
    
    def estimateFromWeightAndCompAndS(self, double average_weight ,  S , double C , double H , double N , double O , double P ):
        """
        estimateFromWeightAndCompAndS(self, average_weight: float , S: int , C: float , H: float , N: float , O: float , P: float ) -> bool
        Fills this EmpiricalFormula with an approximate elemental composition for a given average weight, exact number of sulfurs, and approximate elemental stoichiometry
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(S, int), 'arg S wrong type'
        assert isinstance(C, float), 'arg C wrong type'
        assert isinstance(H, float), 'arg H wrong type'
        assert isinstance(N, float), 'arg N wrong type'
        assert isinstance(O, float), 'arg O wrong type'
        assert isinstance(P, float), 'arg P wrong type'
    
    
    
    
    
    
    
        cdef bool _r = self.inst.get().estimateFromWeightAndCompAndS((<double>average_weight), (<unsigned int>S), (<double>C), (<double>H), (<double>N), (<double>O), (<double>P))
        py_result = <bool>_r
        return py_result
    
    def _getIsotopeDistribution_0(self, CoarseIsotopePatternGenerator in_0 ):
        """
        _getIsotopeDistribution_0(self, in_0: CoarseIsotopePatternGenerator ) -> IsotopeDistribution
        Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator or the FineIsotopePatternGenerator method
        """
        assert isinstance(in_0, CoarseIsotopePatternGenerator), 'arg in_0 wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().getIsotopeDistribution((deref(in_0.inst.get()))))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def _getIsotopeDistribution_1(self, FineIsotopePatternGenerator in_0 ):
        """
        _getIsotopeDistribution_1(self, in_0: FineIsotopePatternGenerator ) -> IsotopeDistribution
        """
        assert isinstance(in_0, FineIsotopePatternGenerator), 'arg in_0 wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().getIsotopeDistribution((deref(in_0.inst.get()))))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def getIsotopeDistribution(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getIsotopeDistribution(self, in_0: CoarseIsotopePatternGenerator ) -> IsotopeDistribution
          :noindex:
        
        Computes the isotope distribution of an empirical formula using the CoarseIsotopePatternGenerator or the FineIsotopePatternGenerator method

        
        .. rubric:: Overload:
        .. py:function:: getIsotopeDistribution(self, in_0: FineIsotopePatternGenerator ) -> IsotopeDistribution
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], CoarseIsotopePatternGenerator)):
            return self._getIsotopeDistribution_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FineIsotopePatternGenerator)):
            return self._getIsotopeDistribution_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getConditionalFragmentIsotopeDist(self, EmpiricalFormula precursor , set precursor_isotopes , CoarseIsotopePatternGenerator method ):
        """
        getConditionalFragmentIsotopeDist(self, precursor: EmpiricalFormula , precursor_isotopes: Set[int] , method: CoarseIsotopePatternGenerator ) -> IsotopeDistribution
        """
        assert isinstance(precursor, EmpiricalFormula), 'arg precursor wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
        assert isinstance(method, CoarseIsotopePatternGenerator), 'arg method wrong type'
    
        cdef libcpp_set[unsigned int] v1 = precursor_isotopes
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().getConditionalFragmentIsotopeDist((deref(precursor.inst.get())), v1, (deref(method.inst.get()))))
        precursor_isotopes.clear()
        precursor_isotopes.update(v1)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def getNumberOfAtoms(self):
        """
        getNumberOfAtoms(self) -> int
        Returns the total number of atoms
        """
        cdef size_t _r = self.inst.get().getNumberOfAtoms()
        py_result = <size_t>_r
        return py_result
    
    def getCharge(self):
        """
        getCharge(self) -> int
        Returns the total charge
        """
        cdef ptrdiff_t _r = self.inst.get().getCharge()
        py_result = <ptrdiff_t>_r
        return py_result
    
    def setCharge(self,  charge ):
        """
        setCharge(self, charge: int ) -> None
        Sets the charge
        """
        assert isinstance(charge, int) and charge >= 0, 'arg charge wrong type'
    
        self.inst.get().setCharge((<ptrdiff_t>charge))
    
    def toString(self):
        """
        toString(self) -> Union[bytes, str, String]
        Returns the formula as a string (charges are not included)
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def getElementalComposition(self):
        """
        getElementalComposition(self) -> Dict[bytes, int]
        Get elemental composition as a hash {'Symbol' -> NrAtoms}
        """
        _r = self.inst.get().toMap()
        py_result = dict()
        cdef libcpp_map[libcpp_string, int].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result[<libcpp_string>(deref(it__r).first)] = <int>(deref(it__r).second)
           inc(it__r)
        return py_result
    
    def isEmpty(self):
        """
        isEmpty(self) -> bool
        Returns true if the formula does not contain a element
        """
        cdef bool _r = self.inst.get().isEmpty()
        py_result = <bool>_r
        return py_result
    
    def isCharged(self):
        """
        isCharged(self) -> bool
        Returns true if charge is not equal to zero
        """
        cdef bool _r = self.inst.get().isCharged()
        py_result = <bool>_r
        return py_result
    
    def hasElement(self, Element element ):
        """
        hasElement(self, element: Element ) -> bool
        Returns true if the formula contains the element
        """
        assert isinstance(element, Element), 'arg element wrong type'
    
        cdef bool _r = self.inst.get().hasElement((element.inst.get()))
        py_result = <bool>_r
        return py_result
    
    def contains(self, EmpiricalFormula ef ):
        """
        contains(self, ef: EmpiricalFormula ) -> bool
        Returns true if all elements from `ef` ( empirical formula ) are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula
        """
        assert isinstance(ef, EmpiricalFormula), 'arg ef wrong type'
    
        cdef bool _r = self.inst.get().contains((deref(ef.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def __add__(EmpiricalFormula self, EmpiricalFormula other not None):
        cdef _EmpiricalFormula * this = self.inst.get()
        cdef _EmpiricalFormula * that = other.inst.get()
        cdef _EmpiricalFormula applied = deref(this) + deref(that)
        cdef EmpiricalFormula result = EmpiricalFormula.__new__(EmpiricalFormula)
        result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(applied))
        return result
    
    def __sub__(EmpiricalFormula self, EmpiricalFormula other not None):
        cdef _EmpiricalFormula * this = self.inst.get()
        cdef _EmpiricalFormula * that = other.inst.get()
        cdef _EmpiricalFormula applied = deref(this) - deref(that)
        cdef EmpiricalFormula result = EmpiricalFormula.__new__(EmpiricalFormula)
        result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(applied))
        return result
    
    def __iadd__(EmpiricalFormula self, EmpiricalFormula other not None):
        cdef _EmpiricalFormula * this = self.inst.get()
        cdef _EmpiricalFormula * that = other.inst.get()
        _iadd(this, that)
        return self
    
    def __isub__(EmpiricalFormula self, EmpiricalFormula other not None):
        cdef _EmpiricalFormula * this = self.inst.get()
        cdef _EmpiricalFormula * that = other.inst.get()
        _isub(this, that)
        return self
    
    def calculateTheoreticalIsotopesNumber(self):
        """
        calculateTheoreticalIsotopesNumber(self) -> float
        """
        cdef double _r = self.inst.get().calculateTheoreticalIsotopesNumber()
        py_result = <double>_r
        return py_result
    
    def __str__(self):
        """
        __str__(self) -> Union[bytes, str, String]
        Returns the formula as a string (charges are not included)
        """
        cdef _String _r = self.inst.get().toString()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, EmpiricalFormula):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef EmpiricalFormula other_casted = other
        cdef EmpiricalFormula self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class FeatureDistance:
    """
    Cython implementation of _FeatureDistance

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureDistance.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureDistance rv = FeatureDistance.__new__(FeatureDistance)
       rv.inst = shared_ptr[_FeatureDistance](new _FeatureDistance(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureDistance rv = FeatureDistance.__new__(FeatureDistance)
       rv.inst = shared_ptr[_FeatureDistance](new _FeatureDistance(deref(self.inst.get())))
       return rv
    
    def _init_0(self, double max_intensity , bool force_constraints ):
        """
        _init_0(self, max_intensity: float , force_constraints: bool ) -> None
        """
        assert isinstance(max_intensity, float), 'arg max_intensity wrong type'
        assert isinstance(force_constraints, pybool_t), 'arg force_constraints wrong type'
    
    
        self.inst = shared_ptr[_FeatureDistance](new _FeatureDistance((<double>max_intensity), (<bool>force_constraints)))
    
    def _init_1(self, FeatureDistance in_0 ):
        """
        _init_1(self, in_0: FeatureDistance ) -> None
        """
        assert isinstance(in_0, FeatureDistance), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureDistance](new _FeatureDistance((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, max_intensity: float , force_constraints: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureDistance ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], pybool_t)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureDistance)):
             self._init_1(*args)
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

cdef class FeatureGroupingAlgorithmKD:
    """
    Cython implementation of _FeatureGroupingAlgorithmKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureGroupingAlgorithmKD.html>`_
      -- Inherits from ['FeatureGroupingAlgorithm', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        A feature grouping algorithm for unlabeled data
        """
        self.inst = shared_ptr[_FeatureGroupingAlgorithmKD](new _FeatureGroupingAlgorithmKD())
    
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

cdef class IsotopeLabelingMDVs:
    """
    Cython implementation of _IsotopeLabelingMDVs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeLabelingMDVs.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsotopeLabelingMDVs rv = IsotopeLabelingMDVs.__new__(IsotopeLabelingMDVs)
       rv.inst = shared_ptr[_IsotopeLabelingMDVs](new _IsotopeLabelingMDVs(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopeLabelingMDVs rv = IsotopeLabelingMDVs.__new__(IsotopeLabelingMDVs)
       rv.inst = shared_ptr[_IsotopeLabelingMDVs](new _IsotopeLabelingMDVs(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IsotopeLabelingMDVs](new _IsotopeLabelingMDVs())
    
    def _init_1(self, IsotopeLabelingMDVs in_0 ):
        """
        _init_1(self, in_0: IsotopeLabelingMDVs ) -> None
        """
        assert isinstance(in_0, IsotopeLabelingMDVs), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopeLabelingMDVs](new _IsotopeLabelingMDVs((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopeLabelingMDVs ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopeLabelingMDVs)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isotopicCorrection(self, Feature normalized_feature , Feature corrected_feature , MatrixDouble correction_matrix , int correction_matrix_agent ):
        """
        isotopicCorrection(self, normalized_feature: Feature , corrected_feature: Feature , correction_matrix: MatrixDouble , correction_matrix_agent: int ) -> None
        This function performs an isotopic correction to account for unlabeled abundances coming from
        the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
        
        
        :param normalized_feature: Feature with normalized values for each component and unlabeled chemical formula for each component group
        :param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
        :param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
        :return: corrected_feature: Feature with corrected values for each component
        """
        assert isinstance(normalized_feature, Feature), 'arg normalized_feature wrong type'
        assert isinstance(corrected_feature, Feature), 'arg corrected_feature wrong type'
        assert isinstance(correction_matrix, MatrixDouble), 'arg correction_matrix wrong type'
        assert correction_matrix_agent in [0, 1, 2], 'arg correction_matrix_agent wrong type'
    
    
    
    
        self.inst.get().isotopicCorrection((deref(normalized_feature.inst.get())), (deref(corrected_feature.inst.get())), (deref(correction_matrix.inst.get())), (<_DerivatizationAgent>correction_matrix_agent))
    
    def isotopicCorrections(self, FeatureMap normalized_featureMap , FeatureMap corrected_featureMap , MatrixDouble correction_matrix , int correction_matrix_agent ):
        """
        isotopicCorrections(self, normalized_featureMap: FeatureMap , corrected_featureMap: FeatureMap , correction_matrix: MatrixDouble , correction_matrix_agent: int ) -> None
        This function performs an isotopic correction to account for unlabeled abundances coming from
        the derivatization agent (e.g., tBDMS) using correction matrix method and is calculated as follows:
        
        
        :param normalized_featuremap: FeatureMap with normalized values for each component and unlabeled chemical formula for each component group
        :param correction_matrix: Square matrix holding correction factors derived either experimentally or theoretically which describe how spectral peaks of naturally abundant 13C contribute to spectral peaks that overlap (or convolve) the spectral peaks of the corrected MDV of the derivatization agent
        :param correction_matrix_agent: Name of the derivatization agent, the internally stored correction matrix if the name of the agent is supplied, only "TBDMS" is supported for now
        :return corrected_featuremap: FeatureMap with corrected values for each component
        """
        assert isinstance(normalized_featureMap, FeatureMap), 'arg normalized_featureMap wrong type'
        assert isinstance(corrected_featureMap, FeatureMap), 'arg corrected_featureMap wrong type'
        assert isinstance(correction_matrix, MatrixDouble), 'arg correction_matrix wrong type'
        assert correction_matrix_agent in [0, 1, 2], 'arg correction_matrix_agent wrong type'
    
    
    
    
        self.inst.get().isotopicCorrections((deref(normalized_featureMap.inst.get())), (deref(corrected_featureMap.inst.get())), (deref(correction_matrix.inst.get())), (<_DerivatizationAgent>correction_matrix_agent))
    
    def calculateIsotopicPurity(self, Feature normalized_feature , list experiment_data ,  isotopic_purity_name ):
        """
        calculateIsotopicPurity(self, normalized_feature: Feature , experiment_data: List[float] , isotopic_purity_name: Union[bytes, str, String] ) -> None
        This function calculates the isotopic purity of the MDV using the following formula:
        isotopic purity of tracer (atom % 13C) = n / [n + (M + n-1)/(M + n)],
        where n in M+n is represented as the index of the result
        The formula is extracted from "High-resolution 13C metabolic flux analysis",
        Long et al, doi:10.1038/s41596-019-0204-0
        
        
        :param normalized_feature: Feature with normalized values for each component and the number of heavy labeled e.g., carbons. Out is a Feature with the calculated isotopic purity for the component group
        :param experiment_data: Vector of experiment data in percent
        :param isotopic_purity_name: Name of the isotopic purity tracer to be saved as a meta value
        """
        assert isinstance(normalized_feature, Feature), 'arg normalized_feature wrong type'
        assert isinstance(experiment_data, list) and all(isinstance(elemt_rec, float) for elemt_rec in experiment_data), 'arg experiment_data wrong type'
        assert (isinstance(isotopic_purity_name, str) or isinstance(isotopic_purity_name, bytes) or isinstance(isotopic_purity_name, String)), 'arg isotopic_purity_name wrong type'
    
        cdef libcpp_vector[double] v1 = experiment_data
    
        self.inst.get().calculateIsotopicPurity((deref(normalized_feature.inst.get())), v1, deref((convString(isotopic_purity_name)).get()))
        
    
    def calculateMDVAccuracy(self, Feature normalized_feature ,  feature_name ,  fragment_isotopomer_theoretical_formula ):
        """
        calculateMDVAccuracy(self, normalized_feature: Feature , feature_name: Union[bytes, str, String] , fragment_isotopomer_theoretical_formula: Union[bytes, str, String] ) -> None
        This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
        using average deviation to the mean. The result is mapped to the meta value "average_accuracy" in the updated feature
        
        
        :param normalized_feature: Feature with normalized values for each component and the chemical formula of the component group. Out is a Feature with the component group accuracy and accuracy for the error for each component
        :param fragment_isotopomer_measured: Measured scan values
        :param fragment_isotopomer_theoretical_formula: Empirical formula from which the theoretical values will be generated
        """
        assert isinstance(normalized_feature, Feature), 'arg normalized_feature wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert (isinstance(fragment_isotopomer_theoretical_formula, str) or isinstance(fragment_isotopomer_theoretical_formula, bytes) or isinstance(fragment_isotopomer_theoretical_formula, String)), 'arg fragment_isotopomer_theoretical_formula wrong type'
    
    
    
        self.inst.get().calculateMDVAccuracy((deref(normalized_feature.inst.get())), deref((convString(feature_name)).get()), deref((convString(fragment_isotopomer_theoretical_formula)).get()))
    
    def calculateMDVAccuracies(self, FeatureMap normalized_featureMap ,  feature_name , dict fragment_isotopomer_theoretical_formulas ):
        """
        calculateMDVAccuracies(self, normalized_featureMap: FeatureMap , feature_name: Union[bytes, str, String] , fragment_isotopomer_theoretical_formulas: Dict[Union[bytes, str], Union[bytes, str]] ) -> None
        This function calculates the accuracy of the MDV as compared to the theoretical MDV (only for 12C quality control experiments)
        using average deviation to the mean
        
        
        param normalized_featuremap: FeatureMap with normalized values for each component and the chemical formula of the component group. Out is a FeatureMap with the component group accuracy and accuracy for the error for each component
        param fragment_isotopomer_measured: Measured scan values
        param fragment_isotopomer_theoretical_formula: A map of ProteinName/peptideRef to Empirical formula from which the theoretical values will be generated
        """
        assert isinstance(normalized_featureMap, FeatureMap), 'arg normalized_featureMap wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert isinstance(fragment_isotopomer_theoretical_formulas, dict) and all(isinstance(k, (bytes, str)) for k in fragment_isotopomer_theoretical_formulas.keys()) and all(isinstance(v, (bytes, str)) for v in fragment_isotopomer_theoretical_formulas.values()), 'arg fragment_isotopomer_theoretical_formulas wrong type'
    
    
        cdef libcpp_map[libcpp_utf8_string, libcpp_utf8_string] * v2 = new libcpp_map[libcpp_utf8_string, libcpp_utf8_string]()
        for key, value in fragment_isotopomer_theoretical_formulas.items():
            if isinstance(key, str):
                key = key.encode('utf-8')
            if isinstance(value, str):
                value = value.encode('utf-8')
            deref(v2)[ (<libcpp_string>key) ] = (<libcpp_string>value)
        
        
        self.inst.get().calculateMDVAccuracies((deref(normalized_featureMap.inst.get())), deref((convString(feature_name)).get()), deref(v2))
        del v2
    
    def calculateMDV(self, Feature measured_feature , Feature normalized_feature , int mass_intensity_type ,  feature_name ):
        """
        calculateMDV(self, measured_feature: Feature , normalized_feature: Feature , mass_intensity_type: int , feature_name: Union[bytes, str, String] ) -> None
        """
        assert isinstance(measured_feature, Feature), 'arg measured_feature wrong type'
        assert isinstance(normalized_feature, Feature), 'arg normalized_feature wrong type'
        assert mass_intensity_type in [0, 1, 2], 'arg mass_intensity_type wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
    
    
    
    
        self.inst.get().calculateMDV((deref(measured_feature.inst.get())), (deref(normalized_feature.inst.get())), (<_MassIntensityType>mass_intensity_type), deref((convString(feature_name)).get()))
    
    def calculateMDVs(self, FeatureMap measured_featureMap , FeatureMap normalized_featureMap , int mass_intensity_type ,  feature_name ):
        """
        calculateMDVs(self, measured_featureMap: FeatureMap , normalized_featureMap: FeatureMap , mass_intensity_type: int , feature_name: Union[bytes, str, String] ) -> None
        """
        assert isinstance(measured_featureMap, FeatureMap), 'arg measured_featureMap wrong type'
        assert isinstance(normalized_featureMap, FeatureMap), 'arg normalized_featureMap wrong type'
        assert mass_intensity_type in [0, 1, 2], 'arg mass_intensity_type wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
    
    
    
    
        self.inst.get().calculateMDVs((deref(measured_featureMap.inst.get())), (deref(normalized_featureMap.inst.get())), (<_MassIntensityType>mass_intensity_type), deref((convString(feature_name)).get())) 

cdef class Kernel_MassTrace:
    """
    Cython implementation of _Kernel_MassTrace

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Kernel_MassTrace.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property fwhm_mz_avg:
        def __set__(self, double fwhm_mz_avg):
        
            self.inst.get().fwhm_mz_avg = (<double>fwhm_mz_avg)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().fwhm_mz_avg
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef Kernel_MassTrace rv = Kernel_MassTrace.__new__(Kernel_MassTrace)
       rv.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Kernel_MassTrace rv = Kernel_MassTrace.__new__(Kernel_MassTrace)
       rv.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace())
    
    def _init_1(self, Kernel_MassTrace in_0 ):
        """
        _init_1(self, in_0: Kernel_MassTrace ) -> None
        """
        assert isinstance(in_0, Kernel_MassTrace), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace((deref(in_0.inst.get()))))
    
    def _init_2(self, list trace_peaks ):
        """
        _init_2(self, trace_peaks: List[Peak2D] ) -> None
        """
        assert isinstance(trace_peaks, list) and all(isinstance(elemt_rec, Peak2D) for elemt_rec in trace_peaks), 'arg trace_peaks wrong type'
        cdef libcpp_vector[_Peak2D] * v0 = new libcpp_vector[_Peak2D]()
        cdef Peak2D item0
        for item0 in trace_peaks:
            v0.push_back(deref(item0.inst.get()))
        self.inst = shared_ptr[_Kernel_MassTrace](new _Kernel_MassTrace(deref(v0)))
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Kernel_MassTrace ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, trace_peaks: List[Peak2D] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Kernel_MassTrace)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, Peak2D) for elemt_rec in args[0])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSize(self):
        """
        getSize(self) -> int
        Returns the number of peaks contained in the mass trace
        """
        cdef size_t _r = self.inst.get().getSize()
        py_result = <size_t>_r
        return py_result
    
    def getLabel(self):
        """
        getLabel(self) -> Union[bytes, str, String]
        Returns label of mass trace
        """
        cdef _String _r = self.inst.get().getLabel()
        py_result = convOutputString(_r)
        return py_result
    
    def setLabel(self,  label ):
        """
        setLabel(self, label: Union[bytes, str, String] ) -> None
        Sets label of mass trace
        """
        assert (isinstance(label, str) or isinstance(label, bytes) or isinstance(label, String)), 'arg label wrong type'
    
        self.inst.get().setLabel(deref((convString(label)).get()))
    
    def getCentroidMZ(self):
        """
        getCentroidMZ(self) -> float
        Returns the centroid m/z
        """
        cdef double _r = self.inst.get().getCentroidMZ()
        py_result = <double>_r
        return py_result
    
    def getCentroidRT(self):
        """
        getCentroidRT(self) -> float
        Returns the centroid RT
        """
        cdef double _r = self.inst.get().getCentroidRT()
        py_result = <double>_r
        return py_result
    
    def getCentroidSD(self):
        """
        getCentroidSD(self) -> float
        Returns the centroid SD
        """
        cdef double _r = self.inst.get().getCentroidSD()
        py_result = <double>_r
        return py_result
    
    def getFWHM(self):
        """
        getFWHM(self) -> float
        Returns FWHM
        """
        cdef double _r = self.inst.get().getFWHM()
        py_result = <double>_r
        return py_result
    
    def getTraceLength(self):
        """
        getTraceLength(self) -> float
        Returns the length of the trace (as difference in RT)
        """
        cdef double _r = self.inst.get().getTraceLength()
        py_result = <double>_r
        return py_result
    
    def getFWHMborders(self):
        """
        getFWHMborders(self) -> List[int, int]
        Returns FWHM boarders
        """
        _r = self.inst.get().getFWHMborders()
        cdef list py_result = [_r.first, _r.second]
        return py_result
    
    def getSmoothedIntensities(self):
        """
        getSmoothedIntensities(self) -> List[float]
        Returns smoothed intensities (empty if no smoothing was explicitly done beforehand!)
        """
        _r = self.inst.get().getSmoothedIntensities()
        cdef list py_result = _r
        return py_result
    
    def getAverageMS1CycleTime(self):
        """
        getAverageMS1CycleTime(self) -> float
        Returns average scan time of mass trace
        """
        cdef double _r = self.inst.get().getAverageMS1CycleTime()
        py_result = <double>_r
        return py_result
    
    def computeSmoothedPeakArea(self):
        """
        computeSmoothedPeakArea(self) -> float
        Sums all non-negative (smoothed!) intensities in the mass trace
        """
        cdef double _r = self.inst.get().computeSmoothedPeakArea()
        py_result = <double>_r
        return py_result
    
    def computePeakArea(self):
        """
        computePeakArea(self) -> float
        Sums intensities of all peaks in the mass trace
        """
        cdef double _r = self.inst.get().computePeakArea()
        py_result = <double>_r
        return py_result
    
    def computeIntensitySum(self):
        """
        computeIntensitySum(self) -> float
        Sum all peak intensities in the mass trace
        """
        cdef double _r = self.inst.get().computeIntensitySum()
        py_result = <double>_r
        return py_result
    
    def findMaxByIntPeak(self, bool in_0 ):
        """
        findMaxByIntPeak(self, in_0: bool ) -> int
        Returns the index of the mass trace's highest peak within the MassTrace container (based either on raw or smoothed intensities)
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        cdef size_t _r = self.inst.get().findMaxByIntPeak((<bool>in_0))
        py_result = <size_t>_r
        return py_result
    
    def estimateFWHM(self, bool in_0 ):
        """
        estimateFWHM(self, in_0: bool ) -> int
        Estimates FWHM of chromatographic peak in seconds (based on either raw or smoothed intensities)
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        cdef size_t _r = self.inst.get().estimateFWHM((<bool>in_0))
        py_result = <size_t>_r
        return py_result
    
    def computeFwhmArea(self):
        """
        computeFwhmArea(self) -> float
        """
        cdef double _r = self.inst.get().computeFwhmArea()
        py_result = <double>_r
        return py_result
    
    def computeFwhmAreaSmooth(self):
        """
        computeFwhmAreaSmooth(self) -> float
        Computes chromatographic peak area within the FWHM range.
        """
        cdef double _r = self.inst.get().computeFwhmAreaSmooth()
        py_result = <double>_r
        return py_result
    
    def getIntensity(self, bool in_0 ):
        """
        getIntensity(self, in_0: bool ) -> float
        Returns the intensity
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getIntensity((<bool>in_0))
        py_result = <double>_r
        return py_result
    
    def getMaxIntensity(self, bool in_0 ):
        """
        getMaxIntensity(self, in_0: bool ) -> float
        Returns the max intensity
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        cdef double _r = self.inst.get().getMaxIntensity((<bool>in_0))
        py_result = <double>_r
        return py_result
    
    def getConvexhull(self):
        """
        getConvexhull(self) -> ConvexHull2D
        Returns the mass trace's convex hull
        """
        cdef _ConvexHull2D * _r = new _ConvexHull2D(self.inst.get().getConvexhull())
        cdef ConvexHull2D py_result = ConvexHull2D.__new__(ConvexHull2D)
        py_result.inst = shared_ptr[_ConvexHull2D](_r)
        return py_result
    
    def setCentroidSD(self, double tmp_sd ):
        """
        setCentroidSD(self, tmp_sd: float ) -> None
        """
        assert isinstance(tmp_sd, float), 'arg tmp_sd wrong type'
    
        self.inst.get().setCentroidSD((<double &>tmp_sd))
    
    def setSmoothedIntensities(self, list db_vec ):
        """
        setSmoothedIntensities(self, db_vec: List[float] ) -> None
        Sets smoothed intensities (smoothing is done externally, e.g. by LowessSmoothing)
        """
        assert isinstance(db_vec, list) and all(isinstance(elemt_rec, float) for elemt_rec in db_vec), 'arg db_vec wrong type'
        cdef libcpp_vector[double] v0 = db_vec
        self.inst.get().setSmoothedIntensities(v0)
        db_vec[:] = v0
    
    def updateSmoothedMaxRT(self):
        """
        updateSmoothedMaxRT(self) -> None
        """
        self.inst.get().updateSmoothedMaxRT()
    
    def updateWeightedMeanRT(self):
        """
        updateWeightedMeanRT(self) -> None
        Compute & update centroid RT as a intensity-weighted mean of RTs
        """
        self.inst.get().updateWeightedMeanRT()
    
    def updateSmoothedWeightedMeanRT(self):
        """
        updateSmoothedWeightedMeanRT(self) -> None
        """
        self.inst.get().updateSmoothedWeightedMeanRT()
    
    def updateMedianRT(self):
        """
        updateMedianRT(self) -> None
        Compute & update centroid RT as median position of intensities
        """
        self.inst.get().updateMedianRT()
    
    def updateMedianMZ(self):
        """
        updateMedianMZ(self) -> None
        Compute & update centroid m/z as median of m/z values
        """
        self.inst.get().updateMedianMZ()
    
    def updateMeanMZ(self):
        """
        updateMeanMZ(self) -> None
        Compute & update centroid m/z as mean of m/z values
        """
        self.inst.get().updateMeanMZ()
    
    def updateWeightedMeanMZ(self):
        """
        updateWeightedMeanMZ(self) -> None
        Compute & update centroid m/z as weighted mean of m/z values
        """
        self.inst.get().updateWeightedMeanMZ()
    
    def updateWeightedMZsd(self):
        """
        updateWeightedMZsd(self) -> None
        Compute & update m/z standard deviation of mass trace as weighted mean of m/z values
        
        Make sure to call update(Weighted)(Mean|Median)MZ() first! <br>
        use getCentroidSD() to get result
        """
        self.inst.get().updateWeightedMZsd()
    
    def setQuantMethod(self, int method ):
        """
        setQuantMethod(self, method: int ) -> None
        Determine if area or median is used for quantification
        """
        assert method in [0, 1, 2], 'arg method wrong type'
    
        self.inst.get().setQuantMethod((<_MT_QUANTMETHOD>method))
    
    def getQuantMethod(self):
        """
        getQuantMethod(self) -> int
        Check if area or median is used for quantification
        """
        cdef _MT_QUANTMETHOD _r = self.inst.get().getQuantMethod()
        py_result = <int>_r
        return py_result 

cdef class MSChromatogram:
    """
    Cython implementation of _MSChromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSChromatogram.html>`_
      -- Inherits from ['ChromatogramSettings', 'RangeManagerRtInt']

    The representation of a chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    Iterations yields access to underlying peak objects but is slower
    Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
    See help(ChromatogramSettings) for information about meta-information
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSChromatogram rv = MSChromatogram.__new__(MSChromatogram)
       rv.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSChromatogram rv = MSChromatogram.__new__(MSChromatogram)
       rv.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MSChromatogram](new _MSChromatogram())
    
    def _init_1(self, MSChromatogram in_0 ):
        """
        _init_1(self, in_0: MSChromatogram ) -> None
        """
        assert isinstance(in_0, MSChromatogram), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSChromatogram](new _MSChromatogram((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSChromatogram ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSChromatogram)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMZ(self):
        """
        getMZ(self) -> float
        Returns the mz of the product entry, makes sense especially for MRM scans
        """
        cdef double _r = self.inst.get().getMZ()
        py_result = <double>_r
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
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def reserve(self,  n ):
        """
        reserve(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().reserve((<size_t>n))
    
    def resize(self,  n ):
        """
        resize(self, n: int ) -> None
        Resize the peak array
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().resize((<size_t>n))
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> ChromatogramPeak
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        cdef _ChromatogramPeak * _r = new _ChromatogramPeak(deref(self.inst.get())[(<size_t>in_0)])
        cdef ChromatogramPeak py_result = ChromatogramPeak.__new__(ChromatogramPeak)
        py_result.inst = shared_ptr[_ChromatogramPeak](_r)
        return py_result
    def __setitem__(self, key, ChromatogramPeak value):
        """Cython signature: ChromatogramPeak & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def updateRanges(self):
        """
        updateRanges(self) -> None
        """
        self.inst.get().updateRanges()
    
    def clear(self,  in_0 ):
        """
        clear(self, in_0: int ) -> None
        Clears all data and meta data
        
        
        :param clear_meta_data: If true, all meta data is cleared in addition to the data
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().clear((<int>in_0))
    
    def push_back(self, ChromatogramPeak in_0 ):
        """
        push_back(self, in_0: ChromatogramPeak ) -> None
        Append a peak
        """
        assert isinstance(in_0, ChromatogramPeak), 'arg in_0 wrong type'
    
        self.inst.get().push_back((deref(in_0.inst.get())))
    
    def isSorted(self):
        """
        isSorted(self) -> bool
        Checks if all peaks are sorted with respect to ascending RT
        """
        cdef bool _r = self.inst.get().isSorted()
        py_result = <bool>_r
        return py_result
    
    def sortByIntensity(self, bool reverse ):
        """
        sortByIntensity(self, reverse: bool ) -> None
        Lexicographically sorts the peaks by their intensity
        
        
        Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly
        """
        assert isinstance(reverse, pybool_t), 'arg reverse wrong type'
    
        self.inst.get().sortByIntensity((<bool>reverse))
    
    def sortByPosition(self):
        """
        sortByPosition(self) -> None
        Lexicographically sorts the peaks by their position
        
        
        The chromatogram is sorted with respect to position. Meta data arrays will be sorted accordingly
        """
        self.inst.get().sortByPosition()
    
    def findNearest(self, double in_0 ):
        """
        findNearest(self, in_0: float ) -> int
        Binary search for the peak nearest to a specific RT
        :note: Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined
        
        
        :param rt: The searched for mass-to-charge ratio searched
        :return: Returns the index of the peak.
        :raises:
          Exception: Precondition is thrown if the chromatogram is empty (not only in debug mode)
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        cdef int _r = self.inst.get().findNearest((<double>in_0))
        py_result = <int>_r
        return py_result
    
    def getFloatDataArrays(self):
        """
        getFloatDataArrays(self) -> List[FloatDataArray]
        Returns a reference to the float meta data arrays
        """
        _r = self.inst.get().getFloatDataArrays()
        py_result = []
        cdef libcpp_vector[_FloatDataArray].iterator it__r = _r.begin()
        cdef FloatDataArray item_py_result
        while it__r != _r.end():
           item_py_result = FloatDataArray.__new__(FloatDataArray)
           item_py_result.inst = shared_ptr[_FloatDataArray](new _FloatDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getIntegerDataArrays(self):
        """
        getIntegerDataArrays(self) -> List[IntegerDataArray]
        Returns a reference to the integer meta data arrays
        """
        _r = self.inst.get().getIntegerDataArrays()
        py_result = []
        cdef libcpp_vector[_IntegerDataArray].iterator it__r = _r.begin()
        cdef IntegerDataArray item_py_result
        while it__r != _r.end():
           item_py_result = IntegerDataArray.__new__(IntegerDataArray)
           item_py_result.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getStringDataArrays(self):
        """
        getStringDataArrays(self) -> List[StringDataArray]
        Returns a reference to the string meta data arrays
        """
        _r = self.inst.get().getStringDataArrays()
        py_result = []
        cdef libcpp_vector[_StringDataArray].iterator it__r = _r.begin()
        cdef StringDataArray item_py_result
        while it__r != _r.end():
           item_py_result = StringDataArray.__new__(StringDataArray)
           item_py_result.inst = shared_ptr[_StringDataArray](new _StringDataArray(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setFloatDataArrays(self, list fda ):
        """
        setFloatDataArrays(self, fda: List[FloatDataArray] ) -> None
        Sets the float meta data arrays
        """
        assert isinstance(fda, list) and all(isinstance(elemt_rec, FloatDataArray) for elemt_rec in fda), 'arg fda wrong type'
        cdef libcpp_vector[_FloatDataArray] * v0 = new libcpp_vector[_FloatDataArray]()
        cdef FloatDataArray item0
        for item0 in fda:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setFloatDataArrays(deref(v0))
        del v0
    
    def setIntegerDataArrays(self, list ida ):
        """
        setIntegerDataArrays(self, ida: List[IntegerDataArray] ) -> None
        Sets the integer meta data arrays
        """
        assert isinstance(ida, list) and all(isinstance(elemt_rec, IntegerDataArray) for elemt_rec in ida), 'arg ida wrong type'
        cdef libcpp_vector[_IntegerDataArray] * v0 = new libcpp_vector[_IntegerDataArray]()
        cdef IntegerDataArray item0
        for item0 in ida:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setIntegerDataArrays(deref(v0))
        del v0
    
    def setStringDataArrays(self, list sda ):
        """
        setStringDataArrays(self, sda: List[StringDataArray] ) -> None
        Sets the string meta data arrays
        """
        assert isinstance(sda, list) and all(isinstance(elemt_rec, StringDataArray) for elemt_rec in sda), 'arg sda wrong type'
        cdef libcpp_vector[_StringDataArray] * v0 = new libcpp_vector[_StringDataArray]()
        cdef StringDataArray item0
        for item0 in sda:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setStringDataArrays(deref(v0))
        del v0
    
    def getProduct(self):
        """
        getProduct(self) -> Product
        Returns the product ion
        """
        cdef _Product * _r = new _Product(self.inst.get().getProduct())
        cdef Product py_result = Product.__new__(Product)
        py_result.inst = shared_ptr[_Product](_r)
        return py_result
    
    def setProduct(self, Product p ):
        """
        setProduct(self, p: Product ) -> None
        Sets the product ion
        """
        assert isinstance(p, Product), 'arg p wrong type'
    
        self.inst.get().setProduct((deref(p.inst.get())))
    
    def getNativeID(self):
        """
        getNativeID(self) -> Union[bytes, str, String]
        Returns the native identifier for the spectrum, used by the acquisition software.
        """
        cdef _String _r = self.inst.get().getNativeID()
        py_result = convOutputString(_r)
        return py_result
    
    def setNativeID(self,  native_id ):
        """
        setNativeID(self, native_id: Union[bytes, str, String] ) -> None
        Sets the native identifier for the spectrum, used by the acquisition software.
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
    
        self.inst.get().setNativeID(deref((convString(native_id)).get()))
    
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
    
    def getInstrumentSettings(self):
        """
        getInstrumentSettings(self) -> InstrumentSettings
        Returns the instrument settings of the current spectrum
        """
        cdef _InstrumentSettings * _r = new _InstrumentSettings(self.inst.get().getInstrumentSettings())
        cdef InstrumentSettings py_result = InstrumentSettings.__new__(InstrumentSettings)
        py_result.inst = shared_ptr[_InstrumentSettings](_r)
        return py_result
    
    def setInstrumentSettings(self, InstrumentSettings instrument_settings ):
        """
        setInstrumentSettings(self, instrument_settings: InstrumentSettings ) -> None
        Sets the instrument settings of the current spectrum
        """
        assert isinstance(instrument_settings, InstrumentSettings), 'arg instrument_settings wrong type'
    
        self.inst.get().setInstrumentSettings((deref(instrument_settings.inst.get())))
    
    def getAcquisitionInfo(self):
        """
        getAcquisitionInfo(self) -> AcquisitionInfo
        Returns the acquisition info
        """
        cdef _AcquisitionInfo * _r = new _AcquisitionInfo(self.inst.get().getAcquisitionInfo())
        cdef AcquisitionInfo py_result = AcquisitionInfo.__new__(AcquisitionInfo)
        py_result.inst = shared_ptr[_AcquisitionInfo](_r)
        return py_result
    
    def setAcquisitionInfo(self, AcquisitionInfo acquisition_info ):
        """
        setAcquisitionInfo(self, acquisition_info: AcquisitionInfo ) -> None
        Sets the acquisition info
        """
        assert isinstance(acquisition_info, AcquisitionInfo), 'arg acquisition_info wrong type'
    
        self.inst.get().setAcquisitionInfo((deref(acquisition_info.inst.get())))
    
    def getSourceFile(self):
        """
        getSourceFile(self) -> SourceFile
        Returns the source file
        """
        cdef _SourceFile * _r = new _SourceFile(self.inst.get().getSourceFile())
        cdef SourceFile py_result = SourceFile.__new__(SourceFile)
        py_result.inst = shared_ptr[_SourceFile](_r)
        return py_result
    
    def setSourceFile(self, SourceFile source_file ):
        """
        setSourceFile(self, source_file: SourceFile ) -> None
        Sets the source file
        """
        assert isinstance(source_file, SourceFile), 'arg source_file wrong type'
    
        self.inst.get().setSourceFile((deref(source_file.inst.get())))
    
    def getPrecursor(self):
        """
        getPrecursor(self) -> Precursor
        Returns the precursors
        """
        cdef _Precursor * _r = new _Precursor(self.inst.get().getPrecursor())
        cdef Precursor py_result = Precursor.__new__(Precursor)
        py_result.inst = shared_ptr[_Precursor](_r)
        return py_result
    
    def setPrecursor(self, Precursor precursor ):
        """
        setPrecursor(self, precursor: Precursor ) -> None
        Sets the precursors
        """
        assert isinstance(precursor, Precursor), 'arg precursor wrong type'
    
        self.inst.get().setPrecursor((deref(precursor.inst.get())))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        Returns the description of the applied processing
        """
        _r = self.inst.get().getDataProcessing()
        py_result = []
        cdef libcpp_vector[shared_ptr[_DataProcessing]].iterator it__r = _r.begin()
        cdef DataProcessing item_py_result
        while it__r != _r.end():
           item_py_result = DataProcessing.__new__(DataProcessing)
           item_py_result.inst = deref(it__r)
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setDataProcessing(self, list in_0 ):
        """
        setDataProcessing(self, in_0: List[DataProcessing] ) -> None
        Sets the description of the applied processing
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, DataProcessing) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[shared_ptr[_DataProcessing]] v0
        cdef DataProcessing in_0_rec
        for in_0_rec in in_0:
            v0.push_back(in_0_rec.inst)
        # call
        self.inst.get().setDataProcessing(v0)
        
    
    def setChromatogramType(self, int type ):
        """
        setChromatogramType(self, type: int ) -> None
        Sets the chromatogram type
        """
        assert type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'arg type wrong type'
    
        self.inst.get().setChromatogramType((<_ChromatogramType>type))
    
    def getChromatogramType(self):
        """
        getChromatogramType(self) -> int
        Get the chromatogram type
        """
        cdef _ChromatogramType _r = self.inst.get().getChromatogramType()
        py_result = <int>_r
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
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, MSChromatogram):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MSChromatogram other_casted = other
        cdef MSChromatogram self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef ChromatogramPeak out
        while it != self.inst.get().end():
            out = ChromatogramPeak.__new__(ChromatogramPeak)
            out.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak(deref(it)))
            yield out
            inc(it)
    


    def get_peaks(self):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        cdef unsigned int n = chrom_.size()
        cdef np.ndarray[np.float64_t, ndim=1] rts
        rts = np.zeros( (n,), dtype=np.float64)
        cdef np.ndarray[np.float32_t, ndim=1] intensities
        intensities = np.zeros( (n,), dtype=np.float32)
        cdef _ChromatogramPeak p

        cdef libcpp_vector[_ChromatogramPeak].iterator it = chrom_.begin()
        cdef int i = 0
        while it != chrom_.end():
            rts[i] = deref(it).getRT()
            intensities[i] = deref(it).getIntensity()
            inc(it)
            i += 1

        return rts, intensities

    def set_peaks(self, peaks):

        assert isinstance(peaks, (tuple, list)), "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"
        assert len(peaks) == 2, "Input for set_peaks needs to be a tuple or a list of size 2 (rt and intensity vector)"

        rts, intensities = peaks
        assert len(rts) == len(intensities), "Input vectors for set_peaks need to have the same length (rt and intensity vector)"

        # Select which function to use for set_peaks:
        # If we have numpy arrays, it helps to use optimized functions
        if isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float32 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_df(rts, intensities)
        elif isinstance(rts, np.ndarray) and isinstance(intensities, np.ndarray) and \
          rts.dtype == np.float64 and intensities.dtype == np.float64 and \
          rts.flags["C_CONTIGUOUS"] and intensities.flags["C_CONTIGUOUS"]  :
            self._set_peaks_fast_dd(rts, intensities)
        else:
            self._set_peaks_orig(rts, intensities)



    def _set_peaks_fast_dd(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[double, ndim=1, mode="c"] data_i not None):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef double intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_fast_df(self, np.ndarray[double, ndim=1, mode="c"] data_rt not None, np.ndarray[float, ndim=1, mode="c"] data_i not None):

        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(data_rt)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(data_rt)

        for i in range(N):
            rt = data_rt[i]
            intensity = data_i[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges()


    def _set_peaks_orig(self, rts, intensities):


        cdef _MSChromatogram * chrom_ = self.inst.get()

        chrom_.resize(0) # empty vector, keep meta data and data arrays
        chrom_.reserve(<int>len(rts)) # allocate space for incoming data
        cdef _ChromatogramPeak p = _ChromatogramPeak()
        cdef double rt
        cdef float intensity
        cdef int N
        N = len(rts)

        for i in range(N):
            rt = rts[i]
            intensity = intensities[i]
            p.setRT(<double>rt)
            p.setIntensity(<float>intensity)
            chrom_.push_back(p)

        chrom_.updateRanges() 

cdef class MzXMLFile:
    """
    Cython implementation of _MzXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzXMLFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzXMLFile rv = MzXMLFile.__new__(MzXMLFile)
       rv.inst = shared_ptr[_MzXMLFile](new _MzXMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzXMLFile rv = MzXMLFile.__new__(MzXMLFile)
       rv.inst = shared_ptr[_MzXMLFile](new _MzXMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzXMLFile](new _MzXMLFile())
    
    def _init_1(self, MzXMLFile in_0 ):
        """
        _init_1(self, in_0: MzXMLFile ) -> None
        """
        assert isinstance(in_0, MzXMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzXMLFile](new _MzXMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzXMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzXMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment exp ):
        """
        load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None
        Loads a MSExperiment from a MzXML file
        
        
        :param exp: MSExperiment
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(exp.inst.get())))
    
    def store(self,  filename , MSExperiment exp ):
        """
        store(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None
        Stores a MSExperiment in a MzXML file
        
        
        :param exp: MSExperiment
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(exp.inst.get())))
    
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
        Sets options for loading/storing
        """
        assert isinstance(in_0, PeakFileOptions), 'arg in_0 wrong type'
    
        self.inst.get().setOptions((deref(in_0.inst.get())))
    
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
    
    def transform(self, path, transformer):
        assert (isinstance(path, str) or isinstance(path, unicode) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

        #
        # the referenced functions _wrap_MSSpectrum, _wrap_MSChromatogram and
        # _wrap_ExperimentalSettings are declared in the MzMLFile.pyx file!
        #

        assert hasattr(transformer, "consumeSpectrum")
        assert hasattr(transformer, "consumeChromatogram")
        assert hasattr(transformer, "setExpectedSize")
        assert hasattr(transformer, "setExperimentalSettings")
        cdef _PythonMSDataConsumer * consumer
        consumer = new _PythonMSDataConsumer(transformer,
                                             _wrap_MSSpectrum_mzxml,
                                             _wrap_MSChromatogram_mzxml,
                                             _wrap_ExperimentalSettings_mzxml)

        try:
            self.inst.get().transform(deref((convString(path)).get()), consumer)
        finally:
            del consumer

cdef _wrap_MSSpectrum_mzxml(const _MSSpectrum & _spec):
    cdef MSSpectrum spec = MSSpectrum.__new__(MSSpectrum)
    spec.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(_spec))
    return spec


cdef _wrap_MSChromatogram_mzxml(const _MSChromatogram & _chromo):
    cdef MSChromatogram chromo = MSChromatogram.__new__(MSChromatogram)
    chromo.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(_chromo))
    return chromo


cdef _wrap_ExperimentalSettings_mzxml(const _ExperimentalSettings & _exp):
    cdef ExperimentalSettings exp = ExperimentalSettings.__new__(ExperimentalSettings)
    exp.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(_exp))
    return exp 

cdef class OPXL_PreprocessedPairSpectra:
    """
    Cython implementation of _OPXL_PreprocessedPairSpectra

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1OPXL_PreprocessedPairSpectra.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property spectra_linear_peaks:
        def __set__(self, MSExperiment spectra_linear_peaks):
        
            self.inst.get().spectra_linear_peaks = (deref(spectra_linear_peaks.inst.get()))
        
    
        def __get__(self):
            cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().spectra_linear_peaks)
            cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
            py_result.inst = shared_ptr[_MSExperiment](_r)
            return py_result
    
    property spectra_xlink_peaks:
        def __set__(self, MSExperiment spectra_xlink_peaks):
        
            self.inst.get().spectra_xlink_peaks = (deref(spectra_xlink_peaks.inst.get()))
        
    
        def __get__(self):
            cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().spectra_xlink_peaks)
            cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
            py_result.inst = shared_ptr[_MSExperiment](_r)
            return py_result
    
    property spectra_all_peaks:
        def __set__(self, MSExperiment spectra_all_peaks):
        
            self.inst.get().spectra_all_peaks = (deref(spectra_all_peaks.inst.get()))
        
    
        def __get__(self):
            cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().spectra_all_peaks)
            cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
            py_result.inst = shared_ptr[_MSExperiment](_r)
            return py_result
    
    def __copy__(self):
       cdef OPXL_PreprocessedPairSpectra rv = OPXL_PreprocessedPairSpectra.__new__(OPXL_PreprocessedPairSpectra)
       rv.inst = shared_ptr[_OPXL_PreprocessedPairSpectra](new _OPXL_PreprocessedPairSpectra(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OPXL_PreprocessedPairSpectra rv = OPXL_PreprocessedPairSpectra.__new__(OPXL_PreprocessedPairSpectra)
       rv.inst = shared_ptr[_OPXL_PreprocessedPairSpectra](new _OPXL_PreprocessedPairSpectra(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  size ):
        """
        _init_0(self, size: int ) -> None
        """
        assert isinstance(size, int) and size >= 0, 'arg size wrong type'
    
        self.inst = shared_ptr[_OPXL_PreprocessedPairSpectra](new _OPXL_PreprocessedPairSpectra((<size_t>size)))
    
    def _init_1(self, OPXL_PreprocessedPairSpectra in_0 ):
        """
        _init_1(self, in_0: OPXL_PreprocessedPairSpectra ) -> None
        """
        assert isinstance(in_0, OPXL_PreprocessedPairSpectra), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OPXL_PreprocessedPairSpectra](new _OPXL_PreprocessedPairSpectra((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, size: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OPXL_PreprocessedPairSpectra ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OPXL_PreprocessedPairSpectra)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OpenSwathHelper:
    """
    Cython implementation of _OpenSwathHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathHelper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenSwathHelper rv = OpenSwathHelper.__new__(OpenSwathHelper)
       rv.inst = shared_ptr[_OpenSwathHelper](new _OpenSwathHelper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenSwathHelper rv = OpenSwathHelper.__new__(OpenSwathHelper)
       rv.inst = shared_ptr[_OpenSwathHelper](new _OpenSwathHelper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenSwathHelper](new _OpenSwathHelper())
    
    def _init_1(self, OpenSwathHelper in_0 ):
        """
        _init_1(self, in_0: OpenSwathHelper ) -> None
        """
        assert isinstance(in_0, OpenSwathHelper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenSwathHelper](new _OpenSwathHelper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenSwathHelper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenSwathHelper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def checkSwathMapAndSelectTransitions(self, MSExperiment exp , TargetedExperiment targeted_exp , TargetedExperiment transition_exp_used , double min_upper_edge_dist ):
        """
        checkSwathMapAndSelectTransitions(self, exp: MSExperiment , targeted_exp: TargetedExperiment , transition_exp_used: TargetedExperiment , min_upper_edge_dist: float ) -> bool
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(targeted_exp, TargetedExperiment), 'arg targeted_exp wrong type'
        assert isinstance(transition_exp_used, TargetedExperiment), 'arg transition_exp_used wrong type'
        assert isinstance(min_upper_edge_dist, float), 'arg min_upper_edge_dist wrong type'
    
    
    
    
        cdef bool _r = self.inst.get().checkSwathMapAndSelectTransitions((deref(exp.inst.get())), (deref(targeted_exp.inst.get())), (deref(transition_exp_used.inst.get())), (<double>min_upper_edge_dist))
        py_result = <bool>_r
        return py_result
    
    def estimateRTRange(self, LightTargetedExperiment exp ):
        """
        estimateRTRange(self, exp: LightTargetedExperiment ) -> List[float, float]
        Computes the min and max retention time value
        
        Estimate the retention time span of a targeted experiment by returning the min/max values in retention time as a pair
        
        
        :return: A std `pair` that contains (min,max)
        """
        assert isinstance(exp, LightTargetedExperiment), 'arg exp wrong type'
    
        _r = self.inst.get().estimateRTRange((deref(exp.inst.get())))
        cdef list py_result = [_r.first, _r.second]
        return py_result
    
    def computePrecursorId(self,  transition_group_id ,  isotope ):
        """
        computePrecursorId(self, transition_group_id: Union[bytes, str, String] , isotope: int ) -> Union[bytes, str, String]
        Computes unique precursor identifier
        
        Uses transition_group_id and isotope number to compute a unique precursor
        id of the form "groupID_Precursor_ix" where x is the isotope number, e.g.
        the monoisotopic precursor would become "groupID_Precursor_i0"
        
        
        :param transition_group_id: Unique id of the transition group (peptide/compound)
        :param isotope: Precursor isotope number
        :return: Unique precursor identifier
        """
        assert (isinstance(transition_group_id, str) or isinstance(transition_group_id, bytes) or isinstance(transition_group_id, String)), 'arg transition_group_id wrong type'
        assert isinstance(isotope, int), 'arg isotope wrong type'
    
    
        cdef _String _r = self.inst.get().computePrecursorId(deref((convString(transition_group_id)).get()), (<int>isotope))
        py_result = convOutputString(_r)
        return py_result 

cdef class PI_PeakArea:
    """
    Cython implementation of _PI_PeakArea

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakArea.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property area:
        def __set__(self, double area):
        
            self.inst.get().area = (<double>area)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().area
            py_result = <double>_r
            return py_result
    
    property height:
        def __set__(self, double height):
        
            self.inst.get().height = (<double>height)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().height
            py_result = <double>_r
            return py_result
    
    property apex_pos:
        def __set__(self, double apex_pos):
        
            self.inst.get().apex_pos = (<double>apex_pos)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().apex_pos
            py_result = <double>_r
            return py_result
    
    property hull_points:
        def __set__(self, np.ndarray[np.float32_t,ndim=2] hull_points):
            cdef libcpp_vector[_DPosition2] _dp_vec_0
            cdef _DPosition2 _dp_0
            cdef int _dp_ii_0
            cdef int _dp_N_0 = hull_points.shape[0]
            for _dp_ii_0 in range(_dp_N_0):
                _dp_0[0] = hull_points[_dp_ii_0,0]
                _dp_0[1] = hull_points[_dp_ii_0,1]
                _dp_vec_0.push_back(_dp_0)
            self.inst.get().hull_points = _dp_vec_0
        
    
        def __get__(self):
            cdef libcpp_vector[_DPosition2] _r = self.inst.get().hull_points
            cdef int _out_n_dpos_vec = _r.size()
            cdef py_result = np.zeros([_out_n_dpos_vec,2], dtype=np.float32)
            cdef libcpp_vector[_DPosition2].iterator _out_it_dpos_vec = _r.begin()
            cdef int _out_ii_dpos_vec = 0
            while _out_it_dpos_vec != _r.end():
                 py_result[_out_ii_dpos_vec, 0] = deref(_out_it_dpos_vec)[0]
                 py_result[_out_ii_dpos_vec, 1] = deref(_out_it_dpos_vec)[1]
                 inc(_out_it_dpos_vec)
                 _out_ii_dpos_vec += 1
            return py_result
    
    def __copy__(self):
       cdef PI_PeakArea rv = PI_PeakArea.__new__(PI_PeakArea)
       rv.inst = shared_ptr[_PI_PeakArea](new _PI_PeakArea(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PI_PeakArea rv = PI_PeakArea.__new__(PI_PeakArea)
       rv.inst = shared_ptr[_PI_PeakArea](new _PI_PeakArea(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PI_PeakArea](new _PI_PeakArea())
    
    def _init_1(self, PI_PeakArea in_0 ):
        """
        _init_1(self, in_0: PI_PeakArea ) -> None
        """
        assert isinstance(in_0, PI_PeakArea), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PI_PeakArea](new _PI_PeakArea((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PI_PeakArea ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PI_PeakArea)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PI_PeakBackground:
    """
    Cython implementation of _PI_PeakBackground

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakBackground.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property area:
        def __set__(self, double area):
        
            self.inst.get().area = (<double>area)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().area
            py_result = <double>_r
            return py_result
    
    property height:
        def __set__(self, double height):
        
            self.inst.get().height = (<double>height)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().height
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef PI_PeakBackground rv = PI_PeakBackground.__new__(PI_PeakBackground)
       rv.inst = shared_ptr[_PI_PeakBackground](new _PI_PeakBackground(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PI_PeakBackground rv = PI_PeakBackground.__new__(PI_PeakBackground)
       rv.inst = shared_ptr[_PI_PeakBackground](new _PI_PeakBackground(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PI_PeakBackground](new _PI_PeakBackground())
    
    def _init_1(self, PI_PeakBackground in_0 ):
        """
        _init_1(self, in_0: PI_PeakBackground ) -> None
        """
        assert isinstance(in_0, PI_PeakBackground), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PI_PeakBackground](new _PI_PeakBackground((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PI_PeakBackground ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PI_PeakBackground)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PI_PeakShapeMetrics:
    """
    Cython implementation of _PI_PeakShapeMetrics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PI_PeakShapeMetrics.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property width_at_5:
        def __set__(self, double width_at_5):
        
            self.inst.get().width_at_5 = (<double>width_at_5)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().width_at_5
            py_result = <double>_r
            return py_result
    
    property width_at_10:
        def __set__(self, double width_at_10):
        
            self.inst.get().width_at_10 = (<double>width_at_10)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().width_at_10
            py_result = <double>_r
            return py_result
    
    property width_at_50:
        def __set__(self, double width_at_50):
        
            self.inst.get().width_at_50 = (<double>width_at_50)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().width_at_50
            py_result = <double>_r
            return py_result
    
    property start_position_at_5:
        def __set__(self, double start_position_at_5):
        
            self.inst.get().start_position_at_5 = (<double>start_position_at_5)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().start_position_at_5
            py_result = <double>_r
            return py_result
    
    property start_position_at_10:
        def __set__(self, double start_position_at_10):
        
            self.inst.get().start_position_at_10 = (<double>start_position_at_10)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().start_position_at_10
            py_result = <double>_r
            return py_result
    
    property start_position_at_50:
        def __set__(self, double start_position_at_50):
        
            self.inst.get().start_position_at_50 = (<double>start_position_at_50)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().start_position_at_50
            py_result = <double>_r
            return py_result
    
    property end_position_at_5:
        def __set__(self, double end_position_at_5):
        
            self.inst.get().end_position_at_5 = (<double>end_position_at_5)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().end_position_at_5
            py_result = <double>_r
            return py_result
    
    property end_position_at_10:
        def __set__(self, double end_position_at_10):
        
            self.inst.get().end_position_at_10 = (<double>end_position_at_10)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().end_position_at_10
            py_result = <double>_r
            return py_result
    
    property end_position_at_50:
        def __set__(self, double end_position_at_50):
        
            self.inst.get().end_position_at_50 = (<double>end_position_at_50)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().end_position_at_50
            py_result = <double>_r
            return py_result
    
    property total_width:
        def __set__(self, double total_width):
        
            self.inst.get().total_width = (<double>total_width)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().total_width
            py_result = <double>_r
            return py_result
    
    property tailing_factor:
        def __set__(self, double tailing_factor):
        
            self.inst.get().tailing_factor = (<double>tailing_factor)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().tailing_factor
            py_result = <double>_r
            return py_result
    
    property asymmetry_factor:
        def __set__(self, double asymmetry_factor):
        
            self.inst.get().asymmetry_factor = (<double>asymmetry_factor)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().asymmetry_factor
            py_result = <double>_r
            return py_result
    
    property slope_of_baseline:
        def __set__(self, double slope_of_baseline):
        
            self.inst.get().slope_of_baseline = (<double>slope_of_baseline)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().slope_of_baseline
            py_result = <double>_r
            return py_result
    
    property baseline_delta_2_height:
        def __set__(self, double baseline_delta_2_height):
        
            self.inst.get().baseline_delta_2_height = (<double>baseline_delta_2_height)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().baseline_delta_2_height
            py_result = <double>_r
            return py_result
    
    property points_across_baseline:
        def __set__(self,  points_across_baseline):
        
            self.inst.get().points_across_baseline = (<int>points_across_baseline)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().points_across_baseline
            py_result = <int>_r
            return py_result
    
    property points_across_half_height:
        def __set__(self,  points_across_half_height):
        
            self.inst.get().points_across_half_height = (<int>points_across_half_height)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().points_across_half_height
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef PI_PeakShapeMetrics rv = PI_PeakShapeMetrics.__new__(PI_PeakShapeMetrics)
       rv.inst = shared_ptr[_PI_PeakShapeMetrics](new _PI_PeakShapeMetrics(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PI_PeakShapeMetrics rv = PI_PeakShapeMetrics.__new__(PI_PeakShapeMetrics)
       rv.inst = shared_ptr[_PI_PeakShapeMetrics](new _PI_PeakShapeMetrics(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PI_PeakShapeMetrics](new _PI_PeakShapeMetrics())
    
    def _init_1(self, PI_PeakShapeMetrics in_0 ):
        """
        _init_1(self, in_0: PI_PeakShapeMetrics ) -> None
        """
        assert isinstance(in_0, PI_PeakShapeMetrics), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PI_PeakShapeMetrics](new _PI_PeakShapeMetrics((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PI_PeakShapeMetrics ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PI_PeakShapeMetrics)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ParamEntry:
    """
    Cython implementation of _ParamEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Param_1_1ParamEntry.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property name:
        def __set__(self, bytes name):
        
            self.inst.get().name = (<libcpp_string>name)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().name
            py_result = <libcpp_string>_r
            return py_result
    
    property description:
        def __set__(self, bytes description):
        
            self.inst.get().description = (<libcpp_string>description)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().description
            py_result = <libcpp_string>_r
            return py_result
    
    property value:
        def __set__(self,  value):
        
            self.inst.get().value = deref(ParamValue(value).inst.get())
        
    
        def __get__(self):
            cdef _ParamValue _r = self.inst.get().value
            cdef ParamValue _value = ParamValue.__new__(ParamValue)
            _value.inst = shared_ptr[_ParamValue](new _ParamValue(_r))
            cdef int _type = _r.valueType()
            cdef object py_result
            if _type == ValueType.STRING_VALUE:
                py_result = _value.toString()
            elif _type == ValueType.INT_VALUE:
                py_result = _value.toInt()
            elif _type == ValueType.DOUBLE_VALUE:
                py_result = _value.toDouble()
            elif _type == ValueType.INT_LIST:
                py_result = _value.toIntVector()
            elif _type == ValueType.DOUBLE_LIST:
                py_result = _value.toDoubleVector()
            elif _type == ValueType.STRING_LIST:
                py_result = _value.toStringVector()
            elif _type == ValueType.EMPTY_VALUE:
                py_result = None
            else:
                raise Exception("ParamValue instance has invalid value type %d" % _type)
            return py_result
    
    property tags:
        def __set__(self, set tags):
            cdef libcpp_set[libcpp_string] v0 = tags
            self.inst.get().tags = v0
        
    
        def __get__(self):
            _r = self.inst.get().tags
            cdef set py_result = _r
            return py_result
    
    property valid_strings:
        def __set__(self, list valid_strings):
            cdef libcpp_vector[libcpp_string] v0 = valid_strings
            self.inst.get().valid_strings = v0
            
    
        def __get__(self):
            _r = self.inst.get().valid_strings
            cdef list py_result = _r
            return py_result
    
    property max_float:
        def __set__(self, double max_float):
        
            self.inst.get().max_float = (<double>max_float)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().max_float
            py_result = <double>_r
            return py_result
    
    property min_float:
        def __set__(self, double min_float):
        
            self.inst.get().min_float = (<double>min_float)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().min_float
            py_result = <double>_r
            return py_result
    
    property max_int:
        def __set__(self,  max_int):
        
            self.inst.get().max_int = (<int>max_int)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().max_int
            py_result = <int>_r
            return py_result
    
    property min_int:
        def __set__(self,  min_int):
        
            self.inst.get().min_int = (<int>min_int)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().min_int
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef ParamEntry rv = ParamEntry.__new__(ParamEntry)
       rv.inst = shared_ptr[_ParamEntry](new _ParamEntry(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ParamEntry rv = ParamEntry.__new__(ParamEntry)
       rv.inst = shared_ptr[_ParamEntry](new _ParamEntry(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ParamEntry](new _ParamEntry())
    
    def _init_1(self, ParamEntry in_0 ):
        """
        _init_1(self, in_0: ParamEntry ) -> None
        """
        assert isinstance(in_0, ParamEntry), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamEntry](new _ParamEntry((deref(in_0.inst.get()))))
    
    def _init_2(self, bytes n ,  v , bytes d , list t ):
        """
        _init_2(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes , t: List[bytes] ) -> None
        """
        assert isinstance(n, bytes), 'arg n wrong type'
        assert isinstance(v, (int, float, list, bytes, str)), 'arg v wrong type'
        assert isinstance(d, bytes), 'arg d wrong type'
        assert isinstance(t, list) and all(isinstance(elemt_rec, bytes) for elemt_rec in t), 'arg t wrong type'
    
    
    
        cdef libcpp_vector[libcpp_string] v3 = t
        self.inst = shared_ptr[_ParamEntry](new _ParamEntry((<libcpp_string>n), deref(ParamValue(v).inst.get()), (<libcpp_string>d), v3))
        
    
    def _init_3(self, bytes n ,  v , bytes d ):
        """
        _init_3(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes ) -> None
        """
        assert isinstance(n, bytes), 'arg n wrong type'
        assert isinstance(v, (int, float, list, bytes, str)), 'arg v wrong type'
        assert isinstance(d, bytes), 'arg d wrong type'
    
    
    
        self.inst = shared_ptr[_ParamEntry](new _ParamEntry((<libcpp_string>n), deref(ParamValue(v).inst.get()), (<libcpp_string>d)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ParamEntry ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes , t: List[bytes] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, n: bytes , v: Union[int, float, bytes, str, List[int], List[float], List[bytes]] , d: bytes ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ParamEntry)):
             self._init_1(*args)
        elif (len(args)==4) and (isinstance(args[0], bytes)) and (isinstance(args[1], (int, float, list, bytes, str))) and (isinstance(args[2], bytes)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, bytes) for elemt_rec in args[3])):
             self._init_2(*args)
        elif (len(args)==3) and (isinstance(args[0], bytes)) and (isinstance(args[1], (int, float, list, bytes, str))) and (isinstance(args[2], bytes)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ParamEntry):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ParamEntry other_casted = other
        cdef ParamEntry self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class PeakIntegrator:
    """
    Cython implementation of _PeakIntegrator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakIntegrator.html>`_
      -- Inherits from ['DefaultParamHandler']

    Compute the area, background and shape metrics of a peak
    
    The area computation is performed in integratePeak() and it supports
    integration by simple sum of the intensity, integration by Simpson's rule
    implementations for an odd number of unequally spaced points or integration
    by the trapezoid rule
    
    The background computation is performed in estimateBackground() and it
    supports three different approaches to baseline correction, namely
    computing a rectangular shape under the peak based on the minimum value of
    the peak borders (vertical_division_min), a rectangular shape based on the
    maximum value of the beak borders (vertical_division_max) or a trapezoidal
    shape based on a straight line between the peak borders (base_to_base)
    
    Peak shape metrics are computed in calculatePeakShapeMetrics() and multiple
    metrics are supported
    
    The containers supported by the methods are MSChromatogram and MSSpectrum
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakIntegrator rv = PeakIntegrator.__new__(PeakIntegrator)
       rv.inst = shared_ptr[_PeakIntegrator](new _PeakIntegrator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakIntegrator rv = PeakIntegrator.__new__(PeakIntegrator)
       rv.inst = shared_ptr[_PeakIntegrator](new _PeakIntegrator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakIntegrator](new _PeakIntegrator())
    
    def _init_1(self, PeakIntegrator in_0 ):
        """
        _init_1(self, in_0: PeakIntegrator ) -> None
        """
        assert isinstance(in_0, PeakIntegrator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakIntegrator](new _PeakIntegrator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakIntegrator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakIntegrator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def _integratePeak_0(self, MSChromatogram chromatogram , double left , double right ):
        """
        _integratePeak_0(self, chromatogram: MSChromatogram , left: float , right: float ) -> PI_PeakArea
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
    
    
    
        cdef _PI_PeakArea * _r = new _PI_PeakArea(self.inst.get().integratePeak((deref(chromatogram.inst.get())), (<double>left), (<double>right)))
        cdef PI_PeakArea py_result = PI_PeakArea.__new__(PI_PeakArea)
        py_result.inst = shared_ptr[_PI_PeakArea](_r)
        return py_result
    
    def _integratePeak_1(self, MSSpectrum spectrum , double left , double right ):
        """
        _integratePeak_1(self, spectrum: MSSpectrum , left: float , right: float ) -> PI_PeakArea
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
    
    
    
        cdef _PI_PeakArea * _r = new _PI_PeakArea(self.inst.get().integratePeak((deref(spectrum.inst.get())), (<double>left), (<double>right)))
        cdef PI_PeakArea py_result = PI_PeakArea.__new__(PI_PeakArea)
        py_result.inst = shared_ptr[_PI_PeakArea](_r)
        return py_result
    
    def integratePeak(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: integratePeak(self, chromatogram: MSChromatogram , left: float , right: float ) -> PI_PeakArea
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: integratePeak(self, spectrum: MSSpectrum , left: float , right: float ) -> PI_PeakArea
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._integratePeak_0(*args)
        elif (len(args)==3) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._integratePeak_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _estimateBackground_0(self, MSChromatogram chromatogram , double left , double right , double peak_apex_pos ):
        """
        _estimateBackground_0(self, chromatogram: MSChromatogram , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
        assert isinstance(peak_apex_pos, float), 'arg peak_apex_pos wrong type'
    
    
    
    
        cdef _PI_PeakBackground * _r = new _PI_PeakBackground(self.inst.get().estimateBackground((deref(chromatogram.inst.get())), (<double>left), (<double>right), (<double>peak_apex_pos)))
        cdef PI_PeakBackground py_result = PI_PeakBackground.__new__(PI_PeakBackground)
        py_result.inst = shared_ptr[_PI_PeakBackground](_r)
        return py_result
    
    def _estimateBackground_1(self, MSSpectrum spectrum , double left , double right , double peak_apex_pos ):
        """
        _estimateBackground_1(self, spectrum: MSSpectrum , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
        assert isinstance(peak_apex_pos, float), 'arg peak_apex_pos wrong type'
    
    
    
    
        cdef _PI_PeakBackground * _r = new _PI_PeakBackground(self.inst.get().estimateBackground((deref(spectrum.inst.get())), (<double>left), (<double>right), (<double>peak_apex_pos)))
        cdef PI_PeakBackground py_result = PI_PeakBackground.__new__(PI_PeakBackground)
        py_result.inst = shared_ptr[_PI_PeakBackground](_r)
        return py_result
    
    def estimateBackground(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: estimateBackground(self, chromatogram: MSChromatogram , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: estimateBackground(self, spectrum: MSSpectrum , left: float , right: float , peak_apex_pos: float ) -> PI_PeakBackground
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._estimateBackground_0(*args)
        elif (len(args)==4) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._estimateBackground_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _calculatePeakShapeMetrics_0(self, MSChromatogram chromatogram , double left , double right , double peak_height , double peak_apex_pos ):
        """
        _calculatePeakShapeMetrics_0(self, chromatogram: MSChromatogram , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
        assert isinstance(peak_height, float), 'arg peak_height wrong type'
        assert isinstance(peak_apex_pos, float), 'arg peak_apex_pos wrong type'
    
    
    
    
    
        cdef _PI_PeakShapeMetrics * _r = new _PI_PeakShapeMetrics(self.inst.get().calculatePeakShapeMetrics((deref(chromatogram.inst.get())), (<double>left), (<double>right), (<double>peak_height), (<double>peak_apex_pos)))
        cdef PI_PeakShapeMetrics py_result = PI_PeakShapeMetrics.__new__(PI_PeakShapeMetrics)
        py_result.inst = shared_ptr[_PI_PeakShapeMetrics](_r)
        return py_result
    
    def _calculatePeakShapeMetrics_1(self, MSSpectrum spectrum , double left , double right , double peak_height , double peak_apex_pos ):
        """
        _calculatePeakShapeMetrics_1(self, spectrum: MSSpectrum , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(left, float), 'arg left wrong type'
        assert isinstance(right, float), 'arg right wrong type'
        assert isinstance(peak_height, float), 'arg peak_height wrong type'
        assert isinstance(peak_apex_pos, float), 'arg peak_apex_pos wrong type'
    
    
    
    
    
        cdef _PI_PeakShapeMetrics * _r = new _PI_PeakShapeMetrics(self.inst.get().calculatePeakShapeMetrics((deref(spectrum.inst.get())), (<double>left), (<double>right), (<double>peak_height), (<double>peak_apex_pos)))
        cdef PI_PeakShapeMetrics py_result = PI_PeakShapeMetrics.__new__(PI_PeakShapeMetrics)
        py_result.inst = shared_ptr[_PI_PeakShapeMetrics](_r)
        return py_result
    
    def calculatePeakShapeMetrics(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: calculatePeakShapeMetrics(self, chromatogram: MSChromatogram , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: calculatePeakShapeMetrics(self, spectrum: MSSpectrum , left: float , right: float , peak_height: float , peak_apex_pos: float ) -> PI_PeakShapeMetrics
          :noindex:
    
        """
        if (len(args)==5) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)) and (isinstance(args[4], float)):
            return self._calculatePeakShapeMetrics_0(*args)
        elif (len(args)==5) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)) and (isinstance(args[4], float)):
            return self._calculatePeakShapeMetrics_1(*args)
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
    DerivatizationAgent = __DerivatizationAgent 
    MassIntensityType = __MassIntensityType 
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
