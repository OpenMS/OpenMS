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
def __static_MRMRTNormalizer_chauvenet(list residuals ,  pos ):
    """
    __static_MRMRTNormalizer_chauvenet(residuals: List[float] , pos: int ) -> bool
    """
    assert isinstance(residuals, list) and all(isinstance(elemt_rec, float) for elemt_rec in residuals), 'arg residuals wrong type'
    assert isinstance(pos, int), 'arg pos wrong type'
    cdef libcpp_vector[double] v0 = residuals

    cdef bool _r = _chauvenet_MRMRTNormalizer(v0, (<int>pos))
    
    py_result = <bool>_r
    return py_result

def __static_MRMRTNormalizer_chauvenet_probability(list residuals ,  pos ):
    """
    __static_MRMRTNormalizer_chauvenet_probability(residuals: List[float] , pos: int ) -> float
    """
    assert isinstance(residuals, list) and all(isinstance(elemt_rec, float) for elemt_rec in residuals), 'arg residuals wrong type'
    assert isinstance(pos, int), 'arg pos wrong type'
    cdef libcpp_vector[double] v0 = residuals

    cdef double _r = _chauvenet_probability_MRMRTNormalizer(v0, (<int>pos))
    
    py_result = <double>_r
    return py_result

def __static_MRMRTNormalizer_computeBinnedCoverage(list rtRange , list pairs ,  nrBins ,  minPeptidesPerBin ,  minBinsFilled ):
    """
    __static_MRMRTNormalizer_computeBinnedCoverage(rtRange: List[float, float] , pairs: List[List[float, float]] , nrBins: int , minPeptidesPerBin: int , minBinsFilled: int ) -> bool
    """
    assert isinstance(rtRange, list) and len(rtRange) == 2 and isinstance(rtRange[0], float) and isinstance(rtRange[1], float), 'arg rtRange wrong type'
    assert isinstance(pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in pairs), 'arg pairs wrong type'
    assert isinstance(nrBins, int), 'arg nrBins wrong type'
    assert isinstance(minPeptidesPerBin, int), 'arg minPeptidesPerBin wrong type'
    assert isinstance(minBinsFilled, int), 'arg minBinsFilled wrong type'
    cdef libcpp_pair[double, double] v0
    v0.first = rtRange[0]
    v0.second = rtRange[1]
    cdef libcpp_vector[libcpp_pair[double,double]] v1 = pairs



    cdef bool _r = _computeBinnedCoverage_MRMRTNormalizer(v0, v1, (<int>nrBins), (<int>minPeptidesPerBin), (<int>minBinsFilled))
    pairs[:] = v1
    py_result = <bool>_r
    return py_result

def __static_MRMRTNormalizer_removeOutliersIterative(list pairs , double rsq_limit , double coverage_limit , bool use_chauvenet , bytes outlier_detection_method ):
    """
    __static_MRMRTNormalizer_removeOutliersIterative(pairs: List[List[float, float]] , rsq_limit: float , coverage_limit: float , use_chauvenet: bool , outlier_detection_method: bytes ) -> List[List[float, float]]
    """
    assert isinstance(pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in pairs), 'arg pairs wrong type'
    assert isinstance(rsq_limit, float), 'arg rsq_limit wrong type'
    assert isinstance(coverage_limit, float), 'arg coverage_limit wrong type'
    assert isinstance(use_chauvenet, pybool_t), 'arg use_chauvenet wrong type'
    assert isinstance(outlier_detection_method, bytes), 'arg outlier_detection_method wrong type'
    cdef libcpp_vector[libcpp_pair[double,double]] v0 = pairs




    _r = _removeOutliersIterative_MRMRTNormalizer(v0, (<double>rsq_limit), (<double>coverage_limit), (<bool>use_chauvenet), (<libcpp_string>outlier_detection_method))
    pairs[:] = v0
    cdef list py_result = _r
    return py_result

def __static_MRMRTNormalizer_removeOutliersRANSAC(list pairs , double rsq_limit , double coverage_limit ,  max_iterations , double max_rt_threshold ,  sampling_size ):
    """
    __static_MRMRTNormalizer_removeOutliersRANSAC(pairs: List[List[float, float]] , rsq_limit: float , coverage_limit: float , max_iterations: int , max_rt_threshold: float , sampling_size: int ) -> List[List[float, float]]
    """
    assert isinstance(pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], float) and isinstance(elemt_rec[1], float) for elemt_rec in pairs), 'arg pairs wrong type'
    assert isinstance(rsq_limit, float), 'arg rsq_limit wrong type'
    assert isinstance(coverage_limit, float), 'arg coverage_limit wrong type'
    assert isinstance(max_iterations, int) and max_iterations >= 0, 'arg max_iterations wrong type'
    assert isinstance(max_rt_threshold, float), 'arg max_rt_threshold wrong type'
    assert isinstance(sampling_size, int) and sampling_size >= 0, 'arg sampling_size wrong type'
    cdef libcpp_vector[libcpp_pair[double,double]] v0 = pairs





    _r = _removeOutliersRANSAC_MRMRTNormalizer(v0, (<double>rsq_limit), (<double>coverage_limit), (<size_t>max_iterations), (<double>max_rt_threshold), (<size_t>sampling_size))
    pairs[:] = v0
    cdef list py_result = _r
    return py_result 

cdef class __SampleState:
    None
    SAMPLENULL = 0
    SOLID = 1
    LIQUID = 2
    GAS = 3
    SOLUTION = 4
    EMULSION = 5
    SUSPENSION = 6
    SIZE_OF_SAMPLESTATE = 7

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AnnotatedMSRun:
    """
    Cython implementation of _AnnotatedMSRun

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AnnotatedMSRun.html>`_

    Class for storing MS run data with peptide and protein identifications
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AnnotatedMSRun](new _AnnotatedMSRun())
    
    def _init_1(self, AnnotatedMSRun in_0 ):
        """
        _init_1(self, in_0: AnnotatedMSRun ) -> None
        """
        assert isinstance(in_0, AnnotatedMSRun), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AnnotatedMSRun](new _AnnotatedMSRun((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AnnotatedMSRun ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AnnotatedMSRun)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def setProteinIdentifications(self, list ids ):
        """
        setProteinIdentifications(self, ids: List[ProteinIdentification] ) -> None
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setProteinIdentifications(deref(v0))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def getPeptideIdentifications(self):
        """
        getPeptideIdentifications(self) -> PeptideIdentificationList
        """
        cdef _PeptideIdentificationList * _r = new _PeptideIdentificationList(self.inst.get().getPeptideIdentifications())
        cdef PeptideIdentificationList py_result = PeptideIdentificationList.__new__(PeptideIdentificationList)
        py_result.inst = shared_ptr[_PeptideIdentificationList](_r)
        return py_result
    
    def setPeptideIdentifications(self, PeptideIdentificationList ids ):
        """
        setPeptideIdentifications(self, ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
    
        self.inst.get().setPeptideIdentifications((deref(ids.inst.get())))
    
    def getMSExperiment(self):
        """
        getMSExperiment(self) -> MSExperiment
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getMSExperiment())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def setMSExperiment(self, MSExperiment experiment ):
        """
        setMSExperiment(self, experiment: MSExperiment ) -> None
        """
        assert isinstance(experiment, MSExperiment), 'arg experiment wrong type'
    
        self.inst.get().setMSExperiment((deref(experiment.inst.get()))) 

cdef class BinnedSpectrum:
    """
    Cython implementation of _BinnedSpectrum

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BinnedSpectrum.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef BinnedSpectrum rv = BinnedSpectrum.__new__(BinnedSpectrum)
       rv.inst = shared_ptr[_BinnedSpectrum](new _BinnedSpectrum(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef BinnedSpectrum rv = BinnedSpectrum.__new__(BinnedSpectrum)
       rv.inst = shared_ptr[_BinnedSpectrum](new _BinnedSpectrum(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BinnedSpectrum](new _BinnedSpectrum())
    
    def _init_1(self, BinnedSpectrum in_0 ):
        """
        _init_1(self, in_0: BinnedSpectrum ) -> None
        """
        assert isinstance(in_0, BinnedSpectrum), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BinnedSpectrum](new _BinnedSpectrum((deref(in_0.inst.get()))))
    
    def _init_2(self, MSSpectrum in_0 , float size , bool unit_ppm ,  spread , float offset ):
        """
        _init_2(self, in_0: MSSpectrum , size: float , unit_ppm: bool , spread: int , offset: float ) -> None
        """
        assert isinstance(in_0, MSSpectrum), 'arg in_0 wrong type'
        assert isinstance(size, float), 'arg size wrong type'
        assert isinstance(unit_ppm, pybool_t), 'arg unit_ppm wrong type'
        assert isinstance(spread, int), 'arg spread wrong type'
        assert isinstance(offset, float), 'arg offset wrong type'
    
    
    
    
    
        self.inst = shared_ptr[_BinnedSpectrum](new _BinnedSpectrum((deref(in_0.inst.get())), (<float>size), (<bool>unit_ppm), (<unsigned int>spread), (<float>offset)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: BinnedSpectrum ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSSpectrum , size: float , unit_ppm: bool , spread: int , offset: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], BinnedSpectrum)):
             self._init_1(*args)
        elif (len(args)==5) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], float)) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], int)) and (isinstance(args[4], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getBinSize(self):
        """
        getBinSize(self) -> float
        Returns the bin size
        """
        cdef float _r = self.inst.get().getBinSize()
        py_result = <float>_r
        return py_result
    
    def getBinSpread(self):
        """
        getBinSpread(self) -> int
        Returns the bin spread
        """
        cdef unsigned int _r = self.inst.get().getBinSpread()
        py_result = <unsigned int>_r
        return py_result
    
    def getBinIndex(self, float mz ):
        """
        getBinIndex(self, mz: float ) -> int
        Returns the bin index of a given m/z position
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef unsigned int _r = self.inst.get().getBinIndex((<float>mz))
        py_result = <unsigned int>_r
        return py_result
    
    def getBinLowerMZ(self,  i ):
        """
        getBinLowerMZ(self, i: int ) -> float
        Returns the lower m/z of a bin given its index
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
    
        cdef float _r = self.inst.get().getBinLowerMZ((<size_t>i))
        py_result = <float>_r
        return py_result
    
    def getBinIntensity(self, double mz ):
        """
        getBinIntensity(self, mz: float ) -> float
        Returns the bin intensity at a given m/z position
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef float _r = self.inst.get().getBinIntensity((<double>mz))
        py_result = <float>_r
        return py_result
    
    def getPrecursors(self):
        """
        getPrecursors(self) -> List[Precursor]
        """
        _r = self.inst.get().getPrecursors()
        py_result = []
        cdef libcpp_vector[_Precursor].iterator it__r = _r.begin()
        cdef Precursor item_py_result
        while it__r != _r.end():
           item_py_result = Precursor.__new__(Precursor)
           item_py_result.inst = shared_ptr[_Precursor](new _Precursor(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def isCompatible(self, BinnedSpectrum a , BinnedSpectrum b ):
        """
        isCompatible(self, a: BinnedSpectrum , b: BinnedSpectrum ) -> bool
        """
        assert isinstance(a, BinnedSpectrum), 'arg a wrong type'
        assert isinstance(b, BinnedSpectrum), 'arg b wrong type'
    
    
        cdef bool _r = self.inst.get().isCompatible((deref(a.inst.get())), (deref(b.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def getOffset(self):
        """
        getOffset(self) -> float
        Returns offset
        """
        cdef float _r = self.inst.get().getOffset()
        py_result = <float>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, BinnedSpectrum):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef BinnedSpectrum other_casted = other
        cdef BinnedSpectrum self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ChromatogramPeak:
    """
    Cython implementation of _ChromatogramPeak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ChromatogramPeak_1_1ChromatogramPeak.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramPeak rv = ChromatogramPeak.__new__(ChromatogramPeak)
       rv.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramPeak rv = ChromatogramPeak.__new__(ChromatogramPeak)
       rv.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A 1-dimensional raw data point or peak for chromatograms
        """
        self.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak())
    
    def _init_1(self, ChromatogramPeak in_0 ):
        """
        _init_1(self, in_0: ChromatogramPeak ) -> None
        """
        assert isinstance(in_0, ChromatogramPeak), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak((deref(in_0.inst.get()))))
    
    def _init_2(self, DPosition1 retention_time , double intensity ):
        """
        _init_2(self, retention_time: DPosition1 , intensity: float ) -> None
        """
        assert isinstance(retention_time, DPosition1), 'arg retention_time wrong type'
        assert isinstance(intensity, float), 'arg intensity wrong type'
    
    
        self.inst = shared_ptr[_ChromatogramPeak](new _ChromatogramPeak((deref(retention_time.inst.get())), (<double>intensity)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A 1-dimensional raw data point or peak for chromatograms

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramPeak ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, retention_time: DPosition1 , intensity: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramPeak)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], DPosition1)) and (isinstance(args[1], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getIntensity(self):
        """
        getIntensity(self) -> float
        Returns the intensity
        """
        cdef double _r = self.inst.get().getIntensity()
        py_result = <double>_r
        return py_result
    
    def setIntensity(self, double in_0 ):
        """
        setIntensity(self, in_0: float ) -> None
        Sets the intensity
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setIntensity((<double>in_0))
    
    def getPosition(self):
        """
        getPosition(self) -> DPosition1
        """
        cdef _DPosition1 * _r = new _DPosition1(self.inst.get().getPosition())
        cdef DPosition1 py_result = DPosition1.__new__(DPosition1)
        py_result.inst = shared_ptr[_DPosition1](_r)
        return py_result
    
    def setPosition(self, DPosition1 in_0 ):
        """
        setPosition(self, in_0: DPosition1 ) -> None
        """
        assert isinstance(in_0, DPosition1), 'arg in_0 wrong type'
    
        self.inst.get().setPosition((deref(in_0.inst.get())))
    
    def getRT(self):
        """
        getRT(self) -> float
        Returns the retention time
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def setRT(self, double in_0 ):
        """
        setRT(self, in_0: float ) -> None
        Sets retention time
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setRT((<double>in_0))
    
    def getPos(self):
        """
        getPos(self) -> float
        Alias for getRT()
        """
        cdef double _r = self.inst.get().getPos()
        py_result = <double>_r
        return py_result
    
    def setPos(self, double in_0 ):
        """
        setPos(self, in_0: float ) -> None
        Alias for setRT()
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setPos((<double>in_0))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ChromatogramPeak):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ChromatogramPeak other_casted = other
        cdef ChromatogramPeak self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ChromeleonFile:
    """
    Cython implementation of _ChromeleonFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromeleonFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromeleonFile rv = ChromeleonFile.__new__(ChromeleonFile)
       rv.inst = shared_ptr[_ChromeleonFile](new _ChromeleonFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromeleonFile rv = ChromeleonFile.__new__(ChromeleonFile)
       rv.inst = shared_ptr[_ChromeleonFile](new _ChromeleonFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Load Chromeleon HPLC text file and save it into a `MSExperiment`.
        """
        self.inst = shared_ptr[_ChromeleonFile](new _ChromeleonFile())
    
    def _init_1(self, ChromeleonFile in_0 ):
        """
        _init_1(self, in_0: ChromeleonFile ) -> None
        """
        assert isinstance(in_0, ChromeleonFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromeleonFile](new _ChromeleonFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Load Chromeleon HPLC text file and save it into a `MSExperiment`.

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromeleonFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromeleonFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment experiment ):
        """
        load(self, filename: Union[bytes, str, String] , experiment: MSExperiment ) -> None
        Load the file's data and metadata, and save it into a `MSExperiment`
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(experiment, MSExperiment), 'arg experiment wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(experiment.inst.get()))) 

cdef class FalseDiscoveryRate:
    """
    Cython implementation of _FalseDiscoveryRate

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FalseDiscoveryRate.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FalseDiscoveryRate](new _FalseDiscoveryRate())
    
    def _apply_0(self, PeptideIdentificationList forward_ids , PeptideIdentificationList reverse_ids ):
        """
        _apply_0(self, forward_ids: PeptideIdentificationList , reverse_ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(forward_ids, PeptideIdentificationList), 'arg forward_ids wrong type'
        assert isinstance(reverse_ids, PeptideIdentificationList), 'arg reverse_ids wrong type'
    
    
        self.inst.get().apply((deref(forward_ids.inst.get())), (deref(reverse_ids.inst.get())))
    
    def _apply_1(self, PeptideIdentificationList id ):
        """
        _apply_1(self, id: PeptideIdentificationList ) -> None
        """
        assert isinstance(id, PeptideIdentificationList), 'arg id wrong type'
    
        self.inst.get().apply((deref(id.inst.get())))
    
    def _apply_2(self, list forward_ids , list reverse_ids ):
        """
        _apply_2(self, forward_ids: List[ProteinIdentification] , reverse_ids: List[ProteinIdentification] ) -> None
        """
        assert isinstance(forward_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in forward_ids), 'arg forward_ids wrong type'
        assert isinstance(reverse_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in reverse_ids), 'arg reverse_ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in forward_ids:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in reverse_ids:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().apply(deref(v0), deref(v1))
        cdef libcpp_vector[_ProteinIdentification].iterator it_reverse_ids = v1.begin()
        replace_0 = []
        while it_reverse_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_reverse_ids)))
            replace_0.append(item1)
            inc(it_reverse_ids)
        reverse_ids[:] = replace_0
        del v1
        cdef libcpp_vector[_ProteinIdentification].iterator it_forward_ids = v0.begin()
        replace_0 = []
        while it_forward_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_forward_ids)))
            replace_0.append(item0)
            inc(it_forward_ids)
        forward_ids[:] = replace_0
        del v0
    
    def _apply_3(self, list id ):
        """
        _apply_3(self, id: List[ProteinIdentification] ) -> None
        """
        assert isinstance(id, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in id), 'arg id wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in id:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().apply(deref(v0))
        cdef libcpp_vector[_ProteinIdentification].iterator it_id = v0.begin()
        replace_0 = []
        while it_id != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_id)))
            replace_0.append(item0)
            inc(it_id)
        id[:] = replace_0
        del v0
    
    def apply(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: apply(self, forward_ids: PeptideIdentificationList , reverse_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: apply(self, id: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: apply(self, forward_ids: List[ProteinIdentification] , reverse_ids: List[ProteinIdentification] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: apply(self, id: List[ProteinIdentification] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], PeptideIdentificationList)):
            return self._apply_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._apply_1(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])):
            return self._apply_2(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])):
            return self._apply_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def applyEstimated(self, list ids ):
        """
        applyEstimated(self, ids: List[ProteinIdentification] ) -> None
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().applyEstimated(deref(v0))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def _applyEvaluateProteinIDs_0(self, list ids , double pepCutoff ,  fpCutoff , double diffWeight ):
        """
        _applyEvaluateProteinIDs_0(self, ids: List[ProteinIdentification] , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(pepCutoff, float), 'arg pepCutoff wrong type'
        assert isinstance(fpCutoff, int), 'arg fpCutoff wrong type'
        assert isinstance(diffWeight, float), 'arg diffWeight wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        cdef double _r = self.inst.get().applyEvaluateProteinIDs(deref(v0), (<double>pepCutoff), (<unsigned int>fpCutoff), (<double>diffWeight))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
        py_result = <double>_r
        return py_result
    
    def _applyEvaluateProteinIDs_1(self, ProteinIdentification ids , double pepCutoff ,  fpCutoff , double diffWeight ):
        """
        _applyEvaluateProteinIDs_1(self, ids: ProteinIdentification , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float
        """
        assert isinstance(ids, ProteinIdentification), 'arg ids wrong type'
        assert isinstance(pepCutoff, float), 'arg pepCutoff wrong type'
        assert isinstance(fpCutoff, int), 'arg fpCutoff wrong type'
        assert isinstance(diffWeight, float), 'arg diffWeight wrong type'
    
    
    
    
        cdef double _r = self.inst.get().applyEvaluateProteinIDs((deref(ids.inst.get())), (<double>pepCutoff), (<unsigned int>fpCutoff), (<double>diffWeight))
        py_result = <double>_r
        return py_result
    
    def applyEvaluateProteinIDs(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: applyEvaluateProteinIDs(self, ids: List[ProteinIdentification] , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: applyEvaluateProteinIDs(self, ids: ProteinIdentification , pepCutoff: float , fpCutoff: int , diffWeight: float ) -> float
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], float)) and (isinstance(args[2], int)) and (isinstance(args[3], float)):
            return self._applyEvaluateProteinIDs_0(*args)
        elif (len(args)==4) and (isinstance(args[0], ProteinIdentification)) and (isinstance(args[1], float)) and (isinstance(args[2], int)) and (isinstance(args[3], float)):
            return self._applyEvaluateProteinIDs_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _applyBasic_0(self, list run_info , PeptideIdentificationList ids ):
        """
        _applyBasic_0(self, run_info: List[ProteinIdentification] , ids: PeptideIdentificationList ) -> None
        """
        assert isinstance(run_info, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in run_info), 'arg run_info wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in run_info:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().applyBasic(deref(v0), (deref(ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_run_info = v0.begin()
        replace_0 = []
        while it_run_info != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_run_info)))
            replace_0.append(item0)
            inc(it_run_info)
        run_info[:] = replace_0
        del v0
    
    def _applyBasic_1(self, PeptideIdentificationList ids , bool higher_score_better ,  charge ,  identifier , bool only_best_per_pep ):
        """
        _applyBasic_1(self, ids: PeptideIdentificationList , higher_score_better: bool , charge: int , identifier: Union[bytes, str, String] , only_best_per_pep: bool ) -> None
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(higher_score_better, pybool_t), 'arg higher_score_better wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
        assert isinstance(only_best_per_pep, pybool_t), 'arg only_best_per_pep wrong type'
    
    
    
    
    
        self.inst.get().applyBasic((deref(ids.inst.get())), (<bool>higher_score_better), (<int>charge), deref((convString(identifier)).get()), (<bool>only_best_per_pep))
    
    def _applyBasic_2(self, ConsensusMap cmap , bool use_unassigned_peptides ):
        """
        _applyBasic_2(self, cmap: ConsensusMap , use_unassigned_peptides: bool ) -> None
        """
        assert isinstance(cmap, ConsensusMap), 'arg cmap wrong type'
        assert isinstance(use_unassigned_peptides, pybool_t), 'arg use_unassigned_peptides wrong type'
    
    
        self.inst.get().applyBasic((deref(cmap.inst.get())), (<bool>use_unassigned_peptides))
    
    def _applyBasic_3(self, ProteinIdentification id , bool groups_too ):
        """
        _applyBasic_3(self, id: ProteinIdentification , groups_too: bool ) -> None
        """
        assert isinstance(id, ProteinIdentification), 'arg id wrong type'
        assert isinstance(groups_too, pybool_t), 'arg groups_too wrong type'
    
    
        self.inst.get().applyBasic((deref(id.inst.get())), (<bool>groups_too))
    
    def applyBasic(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: applyBasic(self, run_info: List[ProteinIdentification] , ids: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: applyBasic(self, ids: PeptideIdentificationList , higher_score_better: bool , charge: int , identifier: Union[bytes, str, String] , only_best_per_pep: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: applyBasic(self, cmap: ConsensusMap , use_unassigned_peptides: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: applyBasic(self, id: ProteinIdentification , groups_too: bool ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], PeptideIdentificationList)):
            return self._applyBasic_0(*args)
        elif (len(args)==5) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], int)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and (isinstance(args[4], pybool_t)):
            return self._applyBasic_1(*args)
        elif (len(args)==2) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], pybool_t)):
            return self._applyBasic_2(*args)
        elif (len(args)==2) and (isinstance(args[0], ProteinIdentification)) and (isinstance(args[1], pybool_t)):
            return self._applyBasic_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def applyPickedProteinFDR(self, ProteinIdentification id ,  decoy_string , bool decoy_prefix , bool groups_too ):
        """
        applyPickedProteinFDR(self, id: ProteinIdentification , decoy_string: String , decoy_prefix: bool , groups_too: bool ) -> None
        """
        assert isinstance(id, ProteinIdentification), 'arg id wrong type'
        assert isinstance(decoy_string, String), 'arg decoy_string wrong type'
        assert isinstance(decoy_prefix, pybool_t), 'arg decoy_prefix wrong type'
        assert isinstance(groups_too, pybool_t), 'arg groups_too wrong type'
    
    
    
    
        self.inst.get().applyPickedProteinFDR((deref(id.inst.get())), deref((<String>decoy_string).inst.get()), (<bool>decoy_prefix), (<bool>groups_too))
    
    def _rocN_0(self, PeptideIdentificationList ids ,  fp_cutoff ):
        """
        _rocN_0(self, ids: PeptideIdentificationList , fp_cutoff: int ) -> float
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(fp_cutoff, int) and fp_cutoff >= 0, 'arg fp_cutoff wrong type'
    
    
        cdef double _r = self.inst.get().rocN((deref(ids.inst.get())), (<size_t>fp_cutoff))
        py_result = <double>_r
        return py_result
    
    def _rocN_1(self, ConsensusMap ids ,  fp_cutoff , bool include_unassigned_peptides ):
        """
        _rocN_1(self, ids: ConsensusMap , fp_cutoff: int , include_unassigned_peptides: bool ) -> float
        """
        assert isinstance(ids, ConsensusMap), 'arg ids wrong type'
        assert isinstance(fp_cutoff, int) and fp_cutoff >= 0, 'arg fp_cutoff wrong type'
        assert isinstance(include_unassigned_peptides, pybool_t), 'arg include_unassigned_peptides wrong type'
    
    
    
        cdef double _r = self.inst.get().rocN((deref(ids.inst.get())), (<size_t>fp_cutoff), (<bool>include_unassigned_peptides))
        py_result = <double>_r
        return py_result
    
    def _rocN_2(self, ConsensusMap ids ,  fp_cutoff ,  identifier , bool include_unassigned_peptides ):
        """
        _rocN_2(self, ids: ConsensusMap , fp_cutoff: int , identifier: Union[bytes, str, String] , include_unassigned_peptides: bool ) -> float
        """
        assert isinstance(ids, ConsensusMap), 'arg ids wrong type'
        assert isinstance(fp_cutoff, int) and fp_cutoff >= 0, 'arg fp_cutoff wrong type'
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
        assert isinstance(include_unassigned_peptides, pybool_t), 'arg include_unassigned_peptides wrong type'
    
    
    
    
        cdef double _r = self.inst.get().rocN((deref(ids.inst.get())), (<size_t>fp_cutoff), deref((convString(identifier)).get()), (<bool>include_unassigned_peptides))
        py_result = <double>_r
        return py_result
    
    def rocN(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: rocN(self, ids: PeptideIdentificationList , fp_cutoff: int ) -> float
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: rocN(self, ids: ConsensusMap , fp_cutoff: int , include_unassigned_peptides: bool ) -> float
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: rocN(self, ids: ConsensusMap , fp_cutoff: int , identifier: Union[bytes, str, String] , include_unassigned_peptides: bool ) -> float
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], int) and args[1] >= 0):
            return self._rocN_0(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], pybool_t)):
            return self._rocN_1(*args)
        elif (len(args)==4) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], int) and args[1] >= 0) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], pybool_t)):
            return self._rocN_2(*args)
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

cdef class ILPDCWrapper:
    """
    Cython implementation of _ILPDCWrapper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ILPDCWrapper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ILPDCWrapper rv = ILPDCWrapper.__new__(ILPDCWrapper)
       rv.inst = shared_ptr[_ILPDCWrapper](new _ILPDCWrapper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ILPDCWrapper rv = ILPDCWrapper.__new__(ILPDCWrapper)
       rv.inst = shared_ptr[_ILPDCWrapper](new _ILPDCWrapper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ILPDCWrapper](new _ILPDCWrapper())
    
    def _init_1(self, ILPDCWrapper in_0 ):
        """
        _init_1(self, in_0: ILPDCWrapper ) -> None
        """
        assert isinstance(in_0, ILPDCWrapper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ILPDCWrapper](new _ILPDCWrapper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ILPDCWrapper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ILPDCWrapper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, FeatureMap fm , list pairs ,  verbose_level ):
        """
        compute(self, fm: FeatureMap , pairs: List[ChargePair] , verbose_level: int ) -> float
        Compute optimal solution and return value of objective function. If the input feature map is empty, a warning is issued and -1 is returned
        """
        assert isinstance(fm, FeatureMap), 'arg fm wrong type'
        assert isinstance(pairs, list) and all(isinstance(elemt_rec, ChargePair) for elemt_rec in pairs), 'arg pairs wrong type'
        assert isinstance(verbose_level, int) and verbose_level >= 0, 'arg verbose_level wrong type'
    
        cdef libcpp_vector[_ChargePair] * v1 = new libcpp_vector[_ChargePair]()
        cdef ChargePair item1
        for item1 in pairs:
            v1.push_back(deref(item1.inst.get()))
    
        cdef double _r = self.inst.get().compute((deref(fm.inst.get())), deref(v1), (<size_t>verbose_level))
        cdef libcpp_vector[_ChargePair].iterator it_pairs = v1.begin()
        replace_0 = []
        while it_pairs != v1.end():
            item1 = ChargePair.__new__(ChargePair)
            item1.inst = shared_ptr[_ChargePair](new _ChargePair(deref(it_pairs)))
            replace_0.append(item1)
            inc(it_pairs)
        pairs[:] = replace_0
        del v1
        py_result = <double>_r
        return py_result 

cdef class ItraqEightPlexQuantitationMethod:
    """
    Cython implementation of _ItraqEightPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ItraqEightPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ItraqEightPlexQuantitationMethod rv = ItraqEightPlexQuantitationMethod.__new__(ItraqEightPlexQuantitationMethod)
       rv.inst = shared_ptr[_ItraqEightPlexQuantitationMethod](new _ItraqEightPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ItraqEightPlexQuantitationMethod rv = ItraqEightPlexQuantitationMethod.__new__(ItraqEightPlexQuantitationMethod)
       rv.inst = shared_ptr[_ItraqEightPlexQuantitationMethod](new _ItraqEightPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        iTRAQ 8 plex quantitation to be used with the IsobaricQuantitation
        """
        self.inst = shared_ptr[_ItraqEightPlexQuantitationMethod](new _ItraqEightPlexQuantitationMethod())
    
    def _init_1(self, ItraqEightPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: ItraqEightPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, ItraqEightPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ItraqEightPlexQuantitationMethod](new _ItraqEightPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        iTRAQ 8 plex quantitation to be used with the IsobaricQuantitation

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ItraqEightPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ItraqEightPlexQuantitationMethod)):
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

cdef class MRMRTNormalizer:
    """
    Cython implementation of _MRMRTNormalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMRTNormalizer.html>`_
    """

    chauvenet = __static_MRMRTNormalizer_chauvenet
    chauvenet_probability = __static_MRMRTNormalizer_chauvenet_probability
    computeBinnedCoverage = __static_MRMRTNormalizer_computeBinnedCoverage
    removeOutliersIterative = __static_MRMRTNormalizer_removeOutliersIterative
    removeOutliersRANSAC = __static_MRMRTNormalizer_removeOutliersRANSAC 

cdef class ModificationsDB:
    """
    Cython implementation of _ModificationsDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationsDB.html>`_
    """

    
    def getNumberOfModifications(self):
        """
        getNumberOfModifications(self) -> int
        Returns the number of modifications read from the unimod.xml file
        """
        cdef size_t _r = self.inst.get().getNumberOfModifications()
        py_result = <size_t>_r
        return py_result
    
    def searchModifications(self, set mods ,  mod_name ,  residue , int term_spec ):
        """
        searchModifications(self, mods: Set[ResidueModification] , mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> None
        Collects all modifications which have the given name as synonym
        
        If `residue` is set, only modifications with matching residue of origin are considered
        If `term_spec` is set, only modifications with matching term specificity are considered
        The resulting set of modifications will be empty if no modification exists that fulfills the criteria
        """
        assert isinstance(mods, set) and all(isinstance(li, ResidueModification) for li in mods), 'arg mods wrong type'
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
        cdef libcpp_set[const _ResidueModification *] * v0 = new libcpp_set[const _ResidueModification *]()
        cdef ResidueModification item0
        for item0 in mods:
           v0.insert((item0.inst.get()))
    
    
    
        self.inst.get().searchModifications(deref(v0), deref((convString(mod_name)).get()), deref((convString(residue)).get()), (<_TermSpecificity>term_spec))
        replace = set()
        cdef libcpp_set[const _ResidueModification *].iterator it_mods = v0.begin()
        while it_mods != v0.end():
           item0 = ResidueModification.__new__(ResidueModification)
           item0.inst = shared_ptr[_ResidueModification](new _ResidueModification(deref(deref(it_mods))))
           replace.add(item0)
           inc(it_mods)
        mods.clear()
        mods.update(replace)
        del v0
    
    def _getModification_0(self,  index ):
        """
        _getModification_0(self, index: int ) -> ResidueModification
        Returns the modification with the given index
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification((<size_t>index)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def _getModification_1(self,  mod_name ):
        """
        _getModification_1(self, mod_name: Union[bytes, str, String] ) -> ResidueModification
        Returns the modification with the given name
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification(deref((convString(mod_name)).get())))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def _getModification_2(self,  mod_name ,  residue , int term_spec ):
        """
        _getModification_2(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
        Returns the modification with the given arguments
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
    
    
    
        cdef const _ResidueModification * __r = (self.inst.get().getModification(deref((convString(mod_name)).get()), deref((convString(residue)).get()), (<_TermSpecificity>term_spec)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getModification(self, index: int ) -> ResidueModification
          :noindex:
        
        Returns the modification with the given index

        
        .. rubric:: Overload:
        .. py:function:: getModification(self, mod_name: Union[bytes, str, String] ) -> ResidueModification
          :noindex:
        
        Returns the modification with the given name

        
        .. rubric:: Overload:
        .. py:function:: getModification(self, mod_name: Union[bytes, str, String] , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
          :noindex:
        
        Returns the modification with the given arguments
    
        """
        if (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
            return self._getModification_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getModification_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (args[2] in [0, 1, 2, 3, 4, 5]):
            return self._getModification_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def has(self,  modification ):
        """
        has(self, modification: Union[bytes, str, String] ) -> bool
        Returns true if the modification exists
        """
        assert (isinstance(modification, str) or isinstance(modification, bytes) or isinstance(modification, String)), 'arg modification wrong type'
    
        cdef bool _r = self.inst.get().has(deref((convString(modification)).get()))
        py_result = <bool>_r
        return py_result
    
    def addModification(self, ResidueModification new_mod ):
        """
        addModification(self, new_mod: ResidueModification ) -> ResidueModification
        Add a new modification to ModificationsDB. If the modification already exists (based on its fullID) it is not added. Returns the modification in the ModificationDB (which can differ from input if mod was already present).
        """
        assert isinstance(new_mod, ResidueModification), 'arg new_mod wrong type'
    
        cdef const _ResidueModification * __r = (self.inst.get().addModification((deref(new_mod.inst.get()))))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def findModificationIndex(self,  mod_name ):
        """
        findModificationIndex(self, mod_name: Union[bytes, str, String] ) -> int
        Returns the index of the modification in the mods_ vector; a unique name must be given
        """
        assert (isinstance(mod_name, str) or isinstance(mod_name, bytes) or isinstance(mod_name, String)), 'arg mod_name wrong type'
    
        cdef size_t _r = self.inst.get().findModificationIndex(deref((convString(mod_name)).get()))
        py_result = <size_t>_r
        return py_result
    
    def searchModificationsByDiffMonoMass(self, list mods , double mass , double max_error ,  residue , int term_spec ):
        """
        searchModificationsByDiffMonoMass(self, mods: List[bytes] , mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> None
        Collects all modifications with delta mass inside a tolerance window
        """
        assert isinstance(mods, list) and all(isinstance(i, bytes) for i in mods), 'arg mods wrong type'
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(max_error, float), 'arg max_error wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in mods:
           v0.push_back(_String(<char *>item0))
    
    
    
    
        self.inst.get().searchModificationsByDiffMonoMass(deref(v0), (<double>mass), (<double>max_error), deref((convString(residue)).get()), (<_TermSpecificity>term_spec))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        mods[:] = replace
        del v0
    
    def getBestModificationByDiffMonoMass(self, double mass , double max_error ,  residue , int term_spec ):
        """
        getBestModificationByDiffMonoMass(self, mass: float , max_error: float , residue: Union[bytes, str, String] , term_spec: int ) -> ResidueModification
        Returns the best matching modification for the given delta mass and residue
        
        Query the modifications DB to get the best matching modification with
        the given delta mass at the given residue (NULL pointer means no result,
        maybe the maximal error tolerance needs to be increased). Possible
        input for CAM modification would be a delta mass of 57 and a residue
        of "C".
        
        Note: If there are multiple possible matches with equal masses, it
        will choose the _first_ match which defaults to the first matching
        UniMod entry.
        
        
        :param residue: The residue at which the modifications occurs
        :param mass: The monoisotopic mass of the residue including the mass of the modification
        :param max_error: The maximal mass error in the modification search
        :return: A pointer to the best matching modification (or NULL if none was found)
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(max_error, float), 'arg max_error wrong type'
        assert (isinstance(residue, str) or isinstance(residue, bytes) or isinstance(residue, String)), 'arg residue wrong type'
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
    
    
    
    
        cdef const _ResidueModification * __r = (self.inst.get().getBestModificationByDiffMonoMass((<double>mass), (<double>max_error), deref((convString(residue)).get()), (<_TermSpecificity>term_spec)))
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def getAllSearchModifications(self, list modifications ):
        """
        getAllSearchModifications(self, modifications: List[bytes] ) -> None
        Collects all modifications that can be used for identification searches
        """
        assert isinstance(modifications, list) and all(isinstance(i, bytes) for i in modifications), 'arg modifications wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in modifications:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllSearchModifications(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        modifications[:] = replace
        del v0
    
    def isInstantiated(self):
        """
        isInstantiated(self) -> bool
        Check whether ModificationsDB was instantiated before
        """
        cdef bool _r = self.inst.get().isInstantiated()
        py_result = <bool>_r
        return py_result
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_ModificationsDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_ModificationsDB](_getInstance_ModificationsDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class ProteinProteinCrossLink:
    """
    Cython implementation of _ProteinProteinCrossLink

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::OPXLDataStructs_1_1ProteinProteinCrossLink.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property alpha:
        def __set__(self, AASequence alpha):
           raise AttributeError("Cannot set constant")
    
        def __get__(self):
            if self.inst.get().alpha is NULL:
                 raise Exception("Cannot access pointer that is NULL")
            cdef const _AASequence * __r = (self.inst.get().alpha)
            if __r == NULL:
                return None
            cdef _AASequence * _r = new _AASequence(deref(__r))
            cdef AASequence py_result = AASequence.__new__(AASequence)
            py_result.inst = shared_ptr[_AASequence](_r)
            return py_result
    
    property beta:
        def __set__(self, AASequence beta):
           raise AttributeError("Cannot set constant")
    
        def __get__(self):
            if self.inst.get().beta is NULL:
                 raise Exception("Cannot access pointer that is NULL")
            cdef const _AASequence * __r = (self.inst.get().beta)
            if __r == NULL:
                return None
            cdef _AASequence * _r = new _AASequence(deref(__r))
            cdef AASequence py_result = AASequence.__new__(AASequence)
            py_result.inst = shared_ptr[_AASequence](_r)
            return py_result
    
    property cross_link_position:
        def __set__(self, list cross_link_position):
            cdef libcpp_pair[ptrdiff_t, ptrdiff_t] v0
            v0.first = cross_link_position[0]
            v0.second = cross_link_position[1]
            self.inst.get().cross_link_position = v0
    
        def __get__(self):
            _r = self.inst.get().cross_link_position
            cdef list py_result = [_r.first, _r.second]
            return py_result
    
    property cross_linker_mass:
        def __set__(self, double cross_linker_mass):
        
            self.inst.get().cross_linker_mass = (<double>cross_linker_mass)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().cross_linker_mass
            py_result = <double>_r
            return py_result
    
    property cross_linker_name:
        def __set__(self,  cross_linker_name):
        
            self.inst.get().cross_linker_name = deref((convString(cross_linker_name)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().cross_linker_name
            py_result = convOutputString(_r)
            return py_result
    
    property term_spec_alpha:
        def __set__(self, int term_spec_alpha):
        
            self.inst.get().term_spec_alpha = (<_TermSpecificity>term_spec_alpha)
        
    
        def __get__(self):
            cdef _TermSpecificity _r = self.inst.get().term_spec_alpha
            py_result = <int>_r
            return py_result
    
    property term_spec_beta:
        def __set__(self, int term_spec_beta):
        
            self.inst.get().term_spec_beta = (<_TermSpecificity>term_spec_beta)
        
    
        def __get__(self):
            cdef _TermSpecificity _r = self.inst.get().term_spec_beta
            py_result = <int>_r
            return py_result
    
    property precursor_correction:
        def __set__(self,  precursor_correction):
        
            self.inst.get().precursor_correction = (<int>precursor_correction)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().precursor_correction
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef ProteinProteinCrossLink rv = ProteinProteinCrossLink.__new__(ProteinProteinCrossLink)
       rv.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteinProteinCrossLink rv = ProteinProteinCrossLink.__new__(ProteinProteinCrossLink)
       rv.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink())
    
    def _init_1(self, ProteinProteinCrossLink in_0 ):
        """
        _init_1(self, in_0: ProteinProteinCrossLink ) -> None
        """
        assert isinstance(in_0, ProteinProteinCrossLink), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteinProteinCrossLink](new _ProteinProteinCrossLink((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteinProteinCrossLink ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteinProteinCrossLink)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getType(self):
        """
        getType(self) -> int
        """
        cdef _ProteinProteinCrossLinkType _r = self.inst.get().getType()
        py_result = <int>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ProteinProteinCrossLink):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ProteinProteinCrossLink other_casted = other
        cdef ProteinProteinCrossLink self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class RansacModelQuadratic:
    """
    Cython implementation of _RansacModelQuadratic

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RansacModelQuadratic.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RansacModelQuadratic rv = RansacModelQuadratic.__new__(RansacModelQuadratic)
       rv.inst = shared_ptr[_RansacModelQuadratic](new _RansacModelQuadratic(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RansacModelQuadratic rv = RansacModelQuadratic.__new__(RansacModelQuadratic)
       rv.inst = shared_ptr[_RansacModelQuadratic](new _RansacModelQuadratic(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RansacModelQuadratic](new _RansacModelQuadratic())
    
    def _init_1(self, RansacModelQuadratic in_0 ):
        """
        _init_1(self, in_0: RansacModelQuadratic ) -> None
        """
        assert isinstance(in_0, RansacModelQuadratic), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RansacModelQuadratic](new _RansacModelQuadratic((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RansacModelQuadratic ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RansacModelQuadratic)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class Sample:
    """
    Cython implementation of _Sample

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Sample.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Sample rv = Sample.__new__(Sample)
       rv.inst = shared_ptr[_Sample](new _Sample(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Sample rv = Sample.__new__(Sample)
       rv.inst = shared_ptr[_Sample](new _Sample(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Sample](new _Sample())
    
    def _init_1(self, Sample in_0 ):
        """
        _init_1(self, in_0: Sample ) -> None
        """
        assert isinstance(in_0, Sample), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Sample](new _Sample((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Sample ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Sample)):
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
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getOrganism(self):
        """
        getOrganism(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOrganism()
        py_result = convOutputString(_r)
        return py_result
    
    def setOrganism(self,  organism ):
        """
        setOrganism(self, organism: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(organism, str) or isinstance(organism, bytes) or isinstance(organism, String)), 'arg organism wrong type'
    
        self.inst.get().setOrganism(deref((convString(organism)).get()))
    
    def getNumber(self):
        """
        getNumber(self) -> Union[bytes, str, String]
        Returns the sample number
        """
        cdef _String _r = self.inst.get().getNumber()
        py_result = convOutputString(_r)
        return py_result
    
    def setNumber(self,  number ):
        """
        setNumber(self, number: Union[bytes, str, String] ) -> None
        Sets the sample number (e.g. sample ID)
        """
        assert (isinstance(number, str) or isinstance(number, bytes) or isinstance(number, String)), 'arg number wrong type'
    
        self.inst.get().setNumber(deref((convString(number)).get()))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the comment (default "")
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  comment ):
        """
        setComment(self, comment: Union[bytes, str, String] ) -> None
        Sets the comment (may contain newline characters)
        """
        assert (isinstance(comment, str) or isinstance(comment, bytes) or isinstance(comment, String)), 'arg comment wrong type'
    
        self.inst.get().setComment(deref((convString(comment)).get()))
    
    def getState(self):
        """
        getState(self) -> int
        Returns the state of aggregation (default SAMPLENULL)
        """
        cdef _SampleState _r = self.inst.get().getState()
        py_result = <int>_r
        return py_result
    
    def setState(self, int state ):
        """
        setState(self, state: int ) -> None
        Sets the state of aggregation
        """
        assert state in [0, 1, 2, 3, 4, 5, 6, 7], 'arg state wrong type'
    
        self.inst.get().setState((<_SampleState>state))
    
    def getMass(self):
        """
        getMass(self) -> float
        Returns the mass (in gram) (default 0.0)
        """
        cdef double _r = self.inst.get().getMass()
        py_result = <double>_r
        return py_result
    
    def setMass(self, double mass ):
        """
        setMass(self, mass: float ) -> None
        Sets the mass (in gram)
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst.get().setMass((<double>mass))
    
    def getVolume(self):
        """
        getVolume(self) -> float
        Returns the volume (in ml) (default 0.0)
        """
        cdef double _r = self.inst.get().getVolume()
        py_result = <double>_r
        return py_result
    
    def setVolume(self, double volume ):
        """
        setVolume(self, volume: float ) -> None
        Sets the volume (in ml)
        """
        assert isinstance(volume, float), 'arg volume wrong type'
    
        self.inst.get().setVolume((<double>volume))
    
    def getConcentration(self):
        """
        getConcentration(self) -> float
        Returns the concentration (in g/l) (default 0.0)
        """
        cdef double _r = self.inst.get().getConcentration()
        py_result = <double>_r
        return py_result
    
    def setConcentration(self, double concentration ):
        """
        setConcentration(self, concentration: float ) -> None
        Sets the concentration (in g/l)
        """
        assert isinstance(concentration, float), 'arg concentration wrong type'
    
        self.inst.get().setConcentration((<double>concentration))
    
    def getSubsamples(self):
        """
        getSubsamples(self) -> List[Sample]
        Returns a reference to the vector of subsamples that were combined to create this sample
        """
        _r = self.inst.get().getSubsamples()
        py_result = []
        cdef libcpp_vector[_Sample].iterator it__r = _r.begin()
        cdef Sample item_py_result
        while it__r != _r.end():
           item_py_result = Sample.__new__(Sample)
           item_py_result.inst = shared_ptr[_Sample](new _Sample(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setSubsamples(self, list subsamples ):
        """
        setSubsamples(self, subsamples: List[Sample] ) -> None
        Sets the vector of subsamples that were combined to create this sample
        """
        assert isinstance(subsamples, list) and all(isinstance(elemt_rec, Sample) for elemt_rec in subsamples), 'arg subsamples wrong type'
        cdef libcpp_vector[_Sample] * v0 = new libcpp_vector[_Sample]()
        cdef Sample item0
        for item0 in subsamples:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSubsamples(deref(v0))
        del v0
    
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
        if not isinstance(other, Sample):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Sample other_casted = other
        cdef Sample self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    SampleState = __SampleState 

cdef class ScanWindow:
    """
    Cython implementation of _ScanWindow

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ScanWindow.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property begin:
        def __set__(self, double begin):
        
            self.inst.get().begin = (<double>begin)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().begin
            py_result = <double>_r
            return py_result
    
    property end:
        def __set__(self, double end):
        
            self.inst.get().end = (<double>end)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().end
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef ScanWindow rv = ScanWindow.__new__(ScanWindow)
       rv.inst = shared_ptr[_ScanWindow](new _ScanWindow(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ScanWindow rv = ScanWindow.__new__(ScanWindow)
       rv.inst = shared_ptr[_ScanWindow](new _ScanWindow(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ScanWindow](new _ScanWindow())
    
    def _init_1(self, ScanWindow in_0 ):
        """
        _init_1(self, in_0: ScanWindow ) -> None
        """
        assert isinstance(in_0, ScanWindow), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ScanWindow](new _ScanWindow((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ScanWindow ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ScanWindow)):
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
        if not isinstance(other, ScanWindow):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ScanWindow other_casted = other
        cdef ScanWindow self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class TMTSixteenPlexQuantitationMethod:
    """
    Cython implementation of _TMTSixteenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTSixteenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TMTSixteenPlexQuantitationMethod rv = TMTSixteenPlexQuantitationMethod.__new__(TMTSixteenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTSixteenPlexQuantitationMethod](new _TMTSixteenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TMTSixteenPlexQuantitationMethod rv = TMTSixteenPlexQuantitationMethod.__new__(TMTSixteenPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTSixteenPlexQuantitationMethod](new _TMTSixteenPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TMTSixteenPlexQuantitationMethod](new _TMTSixteenPlexQuantitationMethod())
    
    def _init_1(self, TMTSixteenPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: TMTSixteenPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, TMTSixteenPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TMTSixteenPlexQuantitationMethod](new _TMTSixteenPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TMTSixteenPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTSixteenPlexQuantitationMethod)):
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
