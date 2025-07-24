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
def __static_MetaboliteSpectralMatching_computeHyperScore(double fragment_mass_error , bool fragment_mass_tolerance_unit_ppm , MSSpectrum exp_spectrum , MSSpectrum db_spectrum , list annotations , double mz_lower_bound ):
    """
    __static_MetaboliteSpectralMatching_computeHyperScore(fragment_mass_error: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , db_spectrum: MSSpectrum , annotations: List[PeptideHit_PeakAnnotation] , mz_lower_bound: float ) -> float
    """
    assert isinstance(fragment_mass_error, float), 'arg fragment_mass_error wrong type'
    assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
    assert isinstance(exp_spectrum, MSSpectrum), 'arg exp_spectrum wrong type'
    assert isinstance(db_spectrum, MSSpectrum), 'arg db_spectrum wrong type'
    assert isinstance(annotations, list) and all(isinstance(elemt_rec, PeptideHit_PeakAnnotation) for elemt_rec in annotations), 'arg annotations wrong type'
    assert isinstance(mz_lower_bound, float), 'arg mz_lower_bound wrong type'




    cdef libcpp_vector[_PeptideHit_PeakAnnotation] * v4 = new libcpp_vector[_PeptideHit_PeakAnnotation]()
    cdef PeptideHit_PeakAnnotation item4
    for item4 in annotations:
        v4.push_back(deref(item4.inst.get()))

    cdef double _r = _computeHyperScore_MetaboliteSpectralMatching((<double>fragment_mass_error), (<bool>fragment_mass_tolerance_unit_ppm), (deref(exp_spectrum.inst.get())), (deref(db_spectrum.inst.get())), deref(v4), (<double>mz_lower_bound))
    cdef libcpp_vector[_PeptideHit_PeakAnnotation].iterator it_annotations = v4.begin()
    replace_0 = []
    while it_annotations != v4.end():
        item4 = PeptideHit_PeakAnnotation.__new__(PeptideHit_PeakAnnotation)
        item4.inst = shared_ptr[_PeptideHit_PeakAnnotation](new _PeptideHit_PeakAnnotation(deref(it_annotations)))
        replace_0.append(item4)
        inc(it_annotations)
    annotations[:] = replace_0
    del v4
    py_result = <double>_r
    return py_result 

cdef class NormalizationMethod:
    None
    NM_SCALE = 0
    NM_SHIFT = 1

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class CVMappingTerm:
    """
    Cython implementation of _CVMappingTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingTerm.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVMappingTerm rv = CVMappingTerm.__new__(CVMappingTerm)
       rv.inst = shared_ptr[_CVMappingTerm](new _CVMappingTerm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVMappingTerm rv = CVMappingTerm.__new__(CVMappingTerm)
       rv.inst = shared_ptr[_CVMappingTerm](new _CVMappingTerm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVMappingTerm](new _CVMappingTerm())
    
    def _init_1(self, CVMappingTerm in_0 ):
        """
        _init_1(self, in_0: CVMappingTerm ) -> None
        """
        assert isinstance(in_0, CVMappingTerm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVMappingTerm](new _CVMappingTerm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVMappingTerm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVMappingTerm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setAccession(self,  accession ):
        """
        setAccession(self, accession: Union[bytes, str, String] ) -> None
        Sets the accession string of the term
        """
        assert (isinstance(accession, str) or isinstance(accession, bytes) or isinstance(accession, String)), 'arg accession wrong type'
    
        self.inst.get().setAccession(deref((convString(accession)).get()))
    
    def getAccession(self):
        """
        getAccession(self) -> Union[bytes, str, String]
        Returns the accession string of the term
        """
        cdef _String _r = self.inst.get().getAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def setUseTermName(self, bool use_term_name ):
        """
        setUseTermName(self, use_term_name: bool ) -> None
        Sets whether the term name should be used, instead of the accession
        """
        assert isinstance(use_term_name, pybool_t), 'arg use_term_name wrong type'
    
        self.inst.get().setUseTermName((<bool>use_term_name))
    
    def getUseTermName(self):
        """
        getUseTermName(self) -> bool
        Returns whether the term name should be used, instead of the accession
        """
        cdef bool _r = self.inst.get().getUseTermName()
        py_result = <bool>_r
        return py_result
    
    def setUseTerm(self, bool use_term ):
        """
        setUseTerm(self, use_term: bool ) -> None
        Sets whether the term itself can be used (or only its children)
        """
        assert isinstance(use_term, pybool_t), 'arg use_term wrong type'
    
        self.inst.get().setUseTerm((<bool>use_term))
    
    def getUseTerm(self):
        """
        getUseTerm(self) -> bool
        Returns true if the term can be used, false if only children are allowed
        """
        cdef bool _r = self.inst.get().getUseTerm()
        py_result = <bool>_r
        return py_result
    
    def setTermName(self,  term_name ):
        """
        setTermName(self, term_name: Union[bytes, str, String] ) -> None
        Sets the name of the term
        """
        assert (isinstance(term_name, str) or isinstance(term_name, bytes) or isinstance(term_name, String)), 'arg term_name wrong type'
    
        self.inst.get().setTermName(deref((convString(term_name)).get()))
    
    def getTermName(self):
        """
        getTermName(self) -> Union[bytes, str, String]
        Returns the name of the term
        """
        cdef _String _r = self.inst.get().getTermName()
        py_result = convOutputString(_r)
        return py_result
    
    def setIsRepeatable(self, bool is_repeatable ):
        """
        setIsRepeatable(self, is_repeatable: bool ) -> None
        Sets whether this term can be repeated
        """
        assert isinstance(is_repeatable, pybool_t), 'arg is_repeatable wrong type'
    
        self.inst.get().setIsRepeatable((<bool>is_repeatable))
    
    def getIsRepeatable(self):
        """
        getIsRepeatable(self) -> bool
        Returns true if this term can be repeated, false otherwise
        """
        cdef bool _r = self.inst.get().getIsRepeatable()
        py_result = <bool>_r
        return py_result
    
    def setAllowChildren(self, bool allow_children ):
        """
        setAllowChildren(self, allow_children: bool ) -> None
        Sets whether children of this term are allowed
        """
        assert isinstance(allow_children, pybool_t), 'arg allow_children wrong type'
    
        self.inst.get().setAllowChildren((<bool>allow_children))
    
    def getAllowChildren(self):
        """
        getAllowChildren(self) -> bool
        Returns true if the children of this term are allowed to be used
        """
        cdef bool _r = self.inst.get().getAllowChildren()
        py_result = <bool>_r
        return py_result
    
    def setCVIdentifierRef(self,  cv_identifier_ref ):
        """
        setCVIdentifierRef(self, cv_identifier_ref: Union[bytes, str, String] ) -> None
        Sets the CV identifier reference string, e.g. UO for unit obo
        """
        assert (isinstance(cv_identifier_ref, str) or isinstance(cv_identifier_ref, bytes) or isinstance(cv_identifier_ref, String)), 'arg cv_identifier_ref wrong type'
    
        self.inst.get().setCVIdentifierRef(deref((convString(cv_identifier_ref)).get()))
    
    def getCVIdentifierRef(self):
        """
        getCVIdentifierRef(self) -> Union[bytes, str, String]
        Returns the CV identifier reference string
        """
        cdef _String _r = self.inst.get().getCVIdentifierRef()
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, CVMappingTerm):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVMappingTerm other_casted = other
        cdef CVMappingTerm self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ConsensusIDAlgorithmAverage:
    """
    Cython implementation of _ConsensusIDAlgorithmAverage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmAverage.html>`_
      -- Inherits from ['ConsensusIDAlgorithmIdentity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmAverage](new _ConsensusIDAlgorithmAverage())
    
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

cdef class ConsensusMapNormalizerAlgorithmMedian:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmMedian

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ConsensusMapNormalizerAlgorithmMedian_1_1ConsensusMapNormalizerAlgorithmMedian.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusMapNormalizerAlgorithmMedian](new _ConsensusMapNormalizerAlgorithmMedian())
    
    def computeMedians(self, ConsensusMap input_map , list medians ,  acc_filter ,  desc_filter ):
        """
        computeMedians(self, input_map: ConsensusMap , medians: List[float] , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> int
        Computes medians of all maps and returns index of map with most features
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
        assert isinstance(medians, list) and all(isinstance(elemt_rec, float) for elemt_rec in medians), 'arg medians wrong type'
        assert (isinstance(acc_filter, str) or isinstance(acc_filter, bytes) or isinstance(acc_filter, String)), 'arg acc_filter wrong type'
        assert (isinstance(desc_filter, str) or isinstance(desc_filter, bytes) or isinstance(desc_filter, String)), 'arg desc_filter wrong type'
    
        cdef libcpp_vector[double] v1 = medians
    
    
        cdef size_t _r = self.inst.get().computeMedians((deref(input_map.inst.get())), v1, deref((convString(acc_filter)).get()), deref((convString(desc_filter)).get()))
        medians[:] = v1
        py_result = <size_t>_r
        return py_result
    
    def normalizeMaps(self, ConsensusMap input_map , int method ,  acc_filter ,  desc_filter ):
        """
        normalizeMaps(self, input_map: ConsensusMap , method: int , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> None
        Normalizes the maps of the consensusMap
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
        assert method in [0, 1], 'arg method wrong type'
        assert (isinstance(acc_filter, str) or isinstance(acc_filter, bytes) or isinstance(acc_filter, String)), 'arg acc_filter wrong type'
        assert (isinstance(desc_filter, str) or isinstance(desc_filter, bytes) or isinstance(desc_filter, String)), 'arg desc_filter wrong type'
    
    
    
    
        self.inst.get().normalizeMaps((deref(input_map.inst.get())), (<_NormalizationMethod>method), deref((convString(acc_filter)).get()), deref((convString(desc_filter)).get())) 

cdef class FeatureFinderAlgorithmMetaboIdent:
    """
    Cython implementation of _FeatureFinderAlgorithmMetaboIdent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderAlgorithmMetaboIdent.html>`_
      -- Inherits from ['DefaultParamHandler']

    Perform targeted feature extraction of compounds provided as table and stores them in features
    
    The algorithms detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications
    Internally, it uses algorithms for targeted data analysis from the OpenSWATH pipeline
    In the simplest case, only CompoundName, SumFormula, Charge and RetentionTime need to be given, all other values may be zero
    Every combination of compound (mass), RT and charge defines one target for feature detection
    Output:
    The main output is a feature map of detected features, with annotations in meta data entries
    Additional outputs are the extracted chromatograms/peak groups, the assay in TraML compatible format, and transformations
    that contain the error between provided and observed peaks
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureFinderAlgorithmMetaboIdent](new _FeatureFinderAlgorithmMetaboIdent())
    
    def setMSData(self, MSExperiment input ):
        """
        setMSData(self, input: MSExperiment ) -> None
        Sets spectra
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
    
        self.inst.get().setMSData((deref(input.inst.get())))
    
    def getMSData(self):
        """
        getMSData(self) -> MSExperiment
        Returns spectra
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getMSData())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def run(self, list metaboIdentTable , FeatureMap features ,  spectra_path ):
        """
        run(self, metaboIdentTable: List[FeatureFinderMetaboIdentCompound] , features: FeatureMap , spectra_path: Union[bytes, str, String] ) -> None
         Run feature extraction. spectra_path get's annotated as primaryMSRunPath in the resulting feature map.
        """
        assert isinstance(metaboIdentTable, list) and all(isinstance(elemt_rec, FeatureFinderMetaboIdentCompound) for elemt_rec in metaboIdentTable), 'arg metaboIdentTable wrong type'
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        assert (isinstance(spectra_path, str) or isinstance(spectra_path, bytes) or isinstance(spectra_path, String)), 'arg spectra_path wrong type'
        cdef libcpp_vector[_FeatureFinderMetaboIdentCompound] * v0 = new libcpp_vector[_FeatureFinderMetaboIdentCompound]()
        cdef FeatureFinderMetaboIdentCompound item0
        for item0 in metaboIdentTable:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().run(deref(v0), (deref(features.inst.get())), deref((convString(spectra_path)).get()))
        del v0
    
    def getChromatograms(self):
        """
        getChromatograms(self) -> MSExperiment
        Retrieves chromatograms (empty if run was not executed)
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getChromatograms())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def getLibrary(self):
        """
        getLibrary(self) -> TargetedExperiment
        Retrieves the assay library (e.g., to store as TraML, empty if run was not executed)
        """
        cdef _TargetedExperiment * _r = new _TargetedExperiment(self.inst.get().getLibrary())
        cdef TargetedExperiment py_result = TargetedExperiment.__new__(TargetedExperiment)
        py_result.inst = shared_ptr[_TargetedExperiment](_r)
        return py_result
    
    def getTransformations(self):
        """
        getTransformations(self) -> TransformationDescription
        Retrieves deviations between provided coordinates and extacted ones (e.g., to store as TrafoXML or for plotting)
        """
        cdef _TransformationDescription * _r = new _TransformationDescription(self.inst.get().getTransformations())
        cdef TransformationDescription py_result = TransformationDescription.__new__(TransformationDescription)
        py_result.inst = shared_ptr[_TransformationDescription](_r)
        return py_result
    
    def getNShared(self):
        """
        getNShared(self) -> int
        Retrieves number of features with shared identifications
        """
        cdef size_t _r = self.inst.get().getNShared()
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

cdef class FeatureFinderMetaboIdentCompound:
    """
    Cython implementation of _FeatureFinderMetaboIdentCompound

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderMetaboIdentCompound.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self,  name ,  formula , double mass , list charges , list rts , list rt_ranges , list iso_distrib ):
        """
        __init__(self, name: Union[bytes, str, String] , formula: Union[bytes, str, String] , mass: float , charges: List[int] , rts: List[float] , rt_ranges: List[float] , iso_distrib: List[float] ) -> None
          Represents a compound in the in the FeatureFinderMetaboIdent library table.
        
        
          :param name: Unique name for the target compound.
          :param formula: Chemical sum formula.
          :param mass: Neutral mass; if zero calculated from formula.
          :param charges: List of possible charge states.
          :param rts: List of possible retention times.
          :param rt_ranges: List of possible retention time ranges (window around RT), either one value or one per RT entry.
          :param iso_distrib: List of relative abundances of isotopologues; if zero calculated from formula.
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(formula, str) or isinstance(formula, bytes) or isinstance(formula, String)), 'arg formula wrong type'
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(charges, list) and all(isinstance(elemt_rec, int) for elemt_rec in charges), 'arg charges wrong type'
        assert isinstance(rts, list) and all(isinstance(elemt_rec, float) for elemt_rec in rts), 'arg rts wrong type'
        assert isinstance(rt_ranges, list) and all(isinstance(elemt_rec, float) for elemt_rec in rt_ranges), 'arg rt_ranges wrong type'
        assert isinstance(iso_distrib, list) and all(isinstance(elemt_rec, float) for elemt_rec in iso_distrib), 'arg iso_distrib wrong type'
    
    
    
        cdef libcpp_vector[int] v3 = charges
        cdef libcpp_vector[double] v4 = rts
        cdef libcpp_vector[double] v5 = rt_ranges
        cdef libcpp_vector[double] v6 = iso_distrib
        self.inst = shared_ptr[_FeatureFinderMetaboIdentCompound](new _FeatureFinderMetaboIdentCompound(deref((convString(name)).get()), deref((convString(formula)).get()), (<double>mass), v3, v4, v5, v6))
        
        
        
        
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
          Gets the compound name.
        
        
          :rtype: str
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def getFormula(self):
        """
        getFormula(self) -> Union[bytes, str, String]
          Gets the compound chemical formula.
        
        
          :rtype: str
        """
        cdef _String _r = self.inst.get().getFormula()
        py_result = convOutputString(_r)
        return py_result
    
    def getMass(self):
        """
        getMass(self) -> float
          Gets the compound mass.
        
        
          :rtype: float
        """
        cdef double _r = self.inst.get().getMass()
        py_result = <double>_r
        return py_result
    
    def getCharges(self):
        """
        getCharges(self) -> List[int]
          Gets the compound charge states.
        
        
          :rtype: list of int
        """
        _r = self.inst.get().getCharges()
        cdef list py_result = _r
        return py_result
    
    def getRTs(self):
        """
        getRTs(self) -> List[float]
          Gets the compound retention times.
        
        
          :rtype: list of float
        """
        _r = self.inst.get().getRTs()
        cdef list py_result = _r
        return py_result
    
    def getRTRanges(self):
        """
        getRTRanges(self) -> List[float]
          Gets the compound retention time ranges.
        
        
          :rtype: list of float
        """
        _r = self.inst.get().getRTRanges()
        cdef list py_result = _r
        return py_result
    
    def getIsotopeDistribution(self):
        """
        getIsotopeDistribution(self) -> List[float]
          Gets the compound isotopic distributions.
        
        
          :rtype: list of float
        """
        _r = self.inst.get().getIsotopeDistribution()
        cdef list py_result = _r
        return py_result 

cdef class FeatureFinderMultiplexAlgorithm:
    """
    Cython implementation of _FeatureFinderMultiplexAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderMultiplexAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureFinderMultiplexAlgorithm rv = FeatureFinderMultiplexAlgorithm.__new__(FeatureFinderMultiplexAlgorithm)
       rv.inst = shared_ptr[_FeatureFinderMultiplexAlgorithm](new _FeatureFinderMultiplexAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureFinderMultiplexAlgorithm rv = FeatureFinderMultiplexAlgorithm.__new__(FeatureFinderMultiplexAlgorithm)
       rv.inst = shared_ptr[_FeatureFinderMultiplexAlgorithm](new _FeatureFinderMultiplexAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureFinderMultiplexAlgorithm](new _FeatureFinderMultiplexAlgorithm())
    
    def _init_1(self, FeatureFinderMultiplexAlgorithm in_0 ):
        """
        _init_1(self, in_0: FeatureFinderMultiplexAlgorithm ) -> None
        """
        assert isinstance(in_0, FeatureFinderMultiplexAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureFinderMultiplexAlgorithm](new _FeatureFinderMultiplexAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureFinderMultiplexAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureFinderMultiplexAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment exp , bool progress ):
        """
        run(self, exp: MSExperiment , progress: bool ) -> None
        Main method for feature detection
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(progress, pybool_t), 'arg progress wrong type'
    
    
        self.inst.get().run((deref(exp.inst.get())), (<bool>progress))
    
    def getFeatureMap(self):
        """
        getFeatureMap(self) -> FeatureMap
        """
        cdef _FeatureMap * _r = new _FeatureMap(self.inst.get().getFeatureMap())
        cdef FeatureMap py_result = FeatureMap.__new__(FeatureMap)
        py_result.inst = shared_ptr[_FeatureMap](_r)
        return py_result
    
    def getConsensusMap(self):
        """
        getConsensusMap(self) -> ConsensusMap
        """
        cdef _ConsensusMap * _r = new _ConsensusMap(self.inst.get().getConsensusMap())
        cdef ConsensusMap py_result = ConsensusMap.__new__(ConsensusMap)
        py_result.inst = shared_ptr[_ConsensusMap](_r)
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

cdef class FeatureMap:
    """
    Cython implementation of _FeatureMap

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureMap.html>`_
      -- Inherits from ['UniqueIdInterface', 'DocumentIdentifier', 'RangeManagerRtMzInt', 'MetaInfoInterface']

    A container for features.
    
    A feature map is a container holding features, which represent
    chemical entities (peptides, proteins, small molecules etc.) found
    in an LC-MS/MS experiment.
    
    This class supports direct iteration in Python.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureMap rv = FeatureMap.__new__(FeatureMap)
       rv.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureMap rv = FeatureMap.__new__(FeatureMap)
       rv.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureMap](new _FeatureMap())
    
    def _init_1(self, FeatureMap in_0 ):
        """
        _init_1(self, in_0: FeatureMap ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureMap](new _FeatureMap((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureMap ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureMap)):
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
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> Feature
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _Feature * _r = new _Feature(deref(self.inst.get())[(<size_t>in_0)])
        cdef Feature py_result = Feature.__new__(Feature)
        py_result.inst = shared_ptr[_Feature](_r)
        return py_result
    def __setitem__(self, key, Feature value):
        """Cython signature: Feature & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def _push_back_0(self, Feature spec ):
        """
        _push_back_0(self, spec: Feature ) -> None
        """
        assert isinstance(spec, Feature), 'arg spec wrong type'
    
        self.inst.get().push_back((deref(spec.inst.get())))
    
    def _push_back_1(self, MRMFeature spec ):
        """
        _push_back_1(self, spec: MRMFeature ) -> None
        """
        assert isinstance(spec, MRMFeature), 'arg spec wrong type'
    
        self.inst.get().push_back((deref(spec.inst.get())))
    
    def push_back(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: push_back(self, spec: Feature ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: push_back(self, spec: MRMFeature ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], Feature)):
            return self._push_back_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeature)):
            return self._push_back_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _sortByIntensity_0(self):
        """
        _sortByIntensity_0(self) -> None
        Sorts the peaks according to ascending intensity
        """
        self.inst.get().sortByIntensity()
    
    def _sortByIntensity_1(self, bool reverse ):
        """
        _sortByIntensity_1(self, reverse: bool ) -> None
        Sorts the peaks according to ascending intensity. Order is reversed if argument is `true` ( reverse = true )
        """
        assert isinstance(reverse, pybool_t), 'arg reverse wrong type'
    
        self.inst.get().sortByIntensity((<bool>reverse))
    
    def sortByIntensity(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: sortByIntensity(self, ) -> None
          :noindex:
        
        Sorts the peaks according to ascending intensity

        
        .. rubric:: Overload:
        .. py:function:: sortByIntensity(self, reverse: bool ) -> None
          :noindex:
        
        Sorts the peaks according to ascending intensity. Order is reversed if argument is `true` ( reverse = true )
    
        """
        if not args:
            return self._sortByIntensity_0(*args)
        elif (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._sortByIntensity_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def sortByPosition(self):
        """
        sortByPosition(self) -> None
        Sorts features by position. Lexicographical comparison (first RT then m/z) is done
        """
        self.inst.get().sortByPosition()
    
    def sortByRT(self):
        """
        sortByRT(self) -> None
        Sorts features by RT position
        """
        self.inst.get().sortByRT()
    
    def sortByMZ(self):
        """
        sortByMZ(self) -> None
        Sorts features by m/z position
        """
        self.inst.get().sortByMZ()
    
    def sortByOverallQuality(self):
        """
        sortByOverallQuality(self) -> None
        Sorts features by ascending overall quality. Order is reversed if argument is `true` ( reverse = true )
        """
        self.inst.get().sortByOverallQuality()
    
    def swap(self, FeatureMap in_0 ):
        """
        swap(self, in_0: FeatureMap ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
    
        self.inst.get().swap((deref(in_0.inst.get())))
    
    def swapFeaturesOnly(self, FeatureMap swapfrom ):
        """
        swapFeaturesOnly(self, swapfrom: FeatureMap ) -> None
        Swaps the feature content (plus its range information) of this map
        """
        assert isinstance(swapfrom, FeatureMap), 'arg swapfrom wrong type'
    
        self.inst.get().swapFeaturesOnly((deref(swapfrom.inst.get())))
    
    def _clear_0(self):
        """
        _clear_0(self) -> None
        Clears all data and meta data
        """
        self.inst.get().clear()
    
    def _clear_1(self, bool clear_meta_data ):
        """
        _clear_1(self, clear_meta_data: bool ) -> None
        Clears all data and meta data. If 'true' is passed as an argument, all meta data is cleared in addition to the data
        """
        assert isinstance(clear_meta_data, pybool_t), 'arg clear_meta_data wrong type'
    
        self.inst.get().clear((<bool>clear_meta_data))
    
    def clear(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: clear(self, ) -> None
          :noindex:
        
        Clears all data and meta data

        
        .. rubric:: Overload:
        .. py:function:: clear(self, clear_meta_data: bool ) -> None
          :noindex:
        
        Clears all data and meta data. If 'true' is passed as an argument, all meta data is cleared in addition to the data
    
        """
        if not args:
            return self._clear_0(*args)
        elif (len(args)==1) and (isinstance(args[0], pybool_t)):
            return self._clear_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __add__(FeatureMap self, FeatureMap other not None):
        cdef _FeatureMap * this = self.inst.get()
        cdef _FeatureMap * that = other.inst.get()
        cdef _FeatureMap applied = deref(this) + deref(that)
        cdef FeatureMap result = FeatureMap.__new__(FeatureMap)
        result.inst = shared_ptr[_FeatureMap](new _FeatureMap(applied))
        return result
    
    def __iadd__(FeatureMap self, FeatureMap other not None):
        cdef _FeatureMap * this = self.inst.get()
        cdef _FeatureMap * that = other.inst.get()
        _iadd(this, that)
        return self
    
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
    
    def setUnassignedPeptideIdentifications(self, PeptideIdentificationList in_0 ):
        """
        setUnassignedPeptideIdentifications(self, in_0: PeptideIdentificationList ) -> None
        Sets the unassigned peptide identifications
        """
        assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    
        self.inst.get().setUnassignedPeptideIdentifications((deref(in_0.inst.get())))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
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
        Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
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
        Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
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
        
        Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)

        
        .. rubric:: Overload:
        .. py:function:: setPrimaryMSRunPath(self, s: List[bytes] , e: MSExperiment ) -> None
          :noindex:
        
        Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`
    
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
        Returns the file path to the first MS run
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
        if not isinstance(other, FeatureMap):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef FeatureMap other_casted = other
        cdef FeatureMap self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef Feature out
        while it != self.inst.get().end():
            out = Feature.__new__(Feature)
            out.inst = shared_ptr[_Feature](new _Feature(deref(it)))
            yield out
            inc(it)
    
    def setUniqueIds(self):
        self.inst.get().applyMemberFunction(address(_setUniqueId)) 

cdef class GaussFilter:
    """
    Cython implementation of _GaussFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GaussFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GaussFilter rv = GaussFilter.__new__(GaussFilter)
       rv.inst = shared_ptr[_GaussFilter](new _GaussFilter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GaussFilter rv = GaussFilter.__new__(GaussFilter)
       rv.inst = shared_ptr[_GaussFilter](new _GaussFilter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data
        """
        self.inst = shared_ptr[_GaussFilter](new _GaussFilter())
    
    def _init_1(self, GaussFilter in_0 ):
        """
        _init_1(self, in_0: GaussFilter ) -> None
        """
        assert isinstance(in_0, GaussFilter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GaussFilter](new _GaussFilter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        This class represents a Gaussian lowpass-filter which works on uniform as well as on non-uniform profile data

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GaussFilter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], GaussFilter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _filter_0(self, MSSpectrum spectrum ):
        """
        _filter_0(self, spectrum: MSSpectrum ) -> None
        Smoothes an MSSpectrum containing profile data
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().filter((deref(spectrum.inst.get())))
    
    def _filter_1(self, MSChromatogram chromatogram ):
        """
        _filter_1(self, chromatogram: MSChromatogram ) -> None
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
    
        self.inst.get().filter((deref(chromatogram.inst.get())))
    
    def filter(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: filter(self, spectrum: MSSpectrum ) -> None
          :noindex:
        
        Smoothes an MSSpectrum containing profile data

        
        .. rubric:: Overload:
        .. py:function:: filter(self, chromatogram: MSChromatogram ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], MSSpectrum)):
            return self._filter_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSChromatogram)):
            return self._filter_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filterExperiment(self, MSExperiment exp ):
        """
        filterExperiment(self, exp: MSExperiment ) -> None
        Smoothes an MSExperiment containing profile data
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().filterExperiment((deref(exp.inst.get())))
    
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

cdef class IndexedMzMLHandler:
    """
    Cython implementation of _IndexedMzMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLHandler.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IndexedMzMLHandler rv = IndexedMzMLHandler.__new__(IndexedMzMLHandler)
       rv.inst = shared_ptr[_IndexedMzMLHandler](new _IndexedMzMLHandler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IndexedMzMLHandler rv = IndexedMzMLHandler.__new__(IndexedMzMLHandler)
       rv.inst = shared_ptr[_IndexedMzMLHandler](new _IndexedMzMLHandler(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IndexedMzMLHandler](new _IndexedMzMLHandler())
    
    def _init_1(self, IndexedMzMLHandler in_0 ):
        """
        _init_1(self, in_0: IndexedMzMLHandler ) -> None
        """
        assert isinstance(in_0, IndexedMzMLHandler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IndexedMzMLHandler](new _IndexedMzMLHandler((deref(in_0.inst.get()))))
    
    def _init_2(self,  filename ):
        """
        _init_2(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst = shared_ptr[_IndexedMzMLHandler](new _IndexedMzMLHandler(deref((convString(filename)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IndexedMzMLHandler ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IndexedMzMLHandler)):
             self._init_1(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def openFile(self,  filename ):
        """
        openFile(self, filename: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().openFile(deref((convString(filename)).get()))
    
    def getParsingSuccess(self):
        """
        getParsingSuccess(self) -> bool
        """
        cdef bool _r = self.inst.get().getParsingSuccess()
        py_result = <bool>_r
        return py_result
    
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
    
    def getSpectrumById(self,  id_ ):
        """
        getSpectrumById(self, id_: int ) -> _Interfaces_Spectrum
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_Spectrum] _r = self.inst.get().getSpectrumById((<int>id_))
        cdef _Interfaces_Spectrum py_result
        py_result = _Interfaces_Spectrum.__new__(_Interfaces_Spectrum)
        py_result.inst = _r
        return py_result
    
    def getChromatogramById(self,  id_ ):
        """
        getChromatogramById(self, id_: int ) -> _Interfaces_Chromatogram
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef shared_ptr[_Chromatogram] _r = self.inst.get().getChromatogramById((<int>id_))
        cdef _Interfaces_Chromatogram py_result
        py_result = _Interfaces_Chromatogram.__new__(_Interfaces_Chromatogram)
        py_result.inst = _r
        return py_result
    
    def getMSSpectrumById(self,  id_ ):
        """
        getMSSpectrumById(self, id_: int ) -> MSSpectrum
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getMSSpectrumById((<int>id_)))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def getMSSpectrumByNativeId(self, bytes id_ , MSSpectrum spec ):
        """
        getMSSpectrumByNativeId(self, id_: bytes , spec: MSSpectrum ) -> None
        """
        assert isinstance(id_, bytes), 'arg id_ wrong type'
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
    
        self.inst.get().getMSSpectrumByNativeId((<libcpp_string>id_), (deref(spec.inst.get())))
    
    def getMSChromatogramById(self,  id_ ):
        """
        getMSChromatogramById(self, id_: int ) -> MSChromatogram
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        cdef _MSChromatogram * _r = new _MSChromatogram(self.inst.get().getMSChromatogramById((<int>id_)))
        cdef MSChromatogram py_result = MSChromatogram.__new__(MSChromatogram)
        py_result.inst = shared_ptr[_MSChromatogram](_r)
        return py_result
    
    def getMSChromatogramByNativeId(self, bytes id_ , MSChromatogram chrom ):
        """
        getMSChromatogramByNativeId(self, id_: bytes , chrom: MSChromatogram ) -> None
        """
        assert isinstance(id_, bytes), 'arg id_ wrong type'
        assert isinstance(chrom, MSChromatogram), 'arg chrom wrong type'
    
    
        self.inst.get().getMSChromatogramByNativeId((<libcpp_string>id_), (deref(chrom.inst.get())))
    
    def setSkipXMLChecks(self, bool skip ):
        """
        setSkipXMLChecks(self, skip: bool ) -> None
        """
        assert isinstance(skip, pybool_t), 'arg skip wrong type'
    
        self.inst.get().setSkipXMLChecks((<bool>skip)) 

cdef class MapConversion:
    """
    Cython implementation of _MapConversion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapConversion.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MapConversion rv = MapConversion.__new__(MapConversion)
       rv.inst = shared_ptr[_MapConversion](new _MapConversion(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MapConversion rv = MapConversion.__new__(MapConversion)
       rv.inst = shared_ptr[_MapConversion](new _MapConversion(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MapConversion](new _MapConversion())
    
    def _init_1(self, MapConversion in_0 ):
        """
        _init_1(self, in_0: MapConversion ) -> None
        """
        assert isinstance(in_0, MapConversion), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MapConversion](new _MapConversion((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MapConversion ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MapConversion)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _convert_0(self,  input_map_index , FeatureMap input_map , ConsensusMap output_map ,  n ):
        """
        _convert_0(self, input_map_index: int , input_map: FeatureMap , output_map: ConsensusMap , n: int ) -> None
        """
        assert isinstance(input_map_index, int) and input_map_index >= 0, 'arg input_map_index wrong type'
        assert isinstance(input_map, FeatureMap), 'arg input_map wrong type'
        assert isinstance(output_map, ConsensusMap), 'arg output_map wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
    
    
    
        self.inst.get().convert((<uint64_t>input_map_index), (deref(input_map.inst.get())), (deref(output_map.inst.get())), (<size_t>n))
    
    def _convert_1(self,  input_map_index , MSExperiment input_map , ConsensusMap output_map ,  n ):
        """
        _convert_1(self, input_map_index: int , input_map: MSExperiment , output_map: ConsensusMap , n: int ) -> None
        """
        assert isinstance(input_map_index, int) and input_map_index >= 0, 'arg input_map_index wrong type'
        assert isinstance(input_map, MSExperiment), 'arg input_map wrong type'
        assert isinstance(output_map, ConsensusMap), 'arg output_map wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
    
    
    
        self.inst.get().convert((<uint64_t>input_map_index), (deref(input_map.inst.get())), (deref(output_map.inst.get())), (<size_t>n))
    
    def _convert_2(self, ConsensusMap input_map , bool keep_uids , FeatureMap output_map ):
        """
        _convert_2(self, input_map: ConsensusMap , keep_uids: bool , output_map: FeatureMap ) -> None
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
        assert isinstance(keep_uids, pybool_t), 'arg keep_uids wrong type'
        assert isinstance(output_map, FeatureMap), 'arg output_map wrong type'
    
    
    
        self.inst.get().convert((deref(input_map.inst.get())), (<bool>keep_uids), (deref(output_map.inst.get())))
    
    def convert(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: convert(self, input_map_index: int , input_map: FeatureMap , output_map: ConsensusMap , n: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: convert(self, input_map_index: int , input_map: MSExperiment , output_map: ConsensusMap , n: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: convert(self, input_map: ConsensusMap , keep_uids: bool , output_map: FeatureMap ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], FeatureMap)) and (isinstance(args[2], ConsensusMap)) and (isinstance(args[3], int) and args[3] >= 0):
            return self._convert_0(*args)
        elif (len(args)==4) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], MSExperiment)) and (isinstance(args[2], ConsensusMap)) and (isinstance(args[3], int) and args[3] >= 0):
            return self._convert_1(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], FeatureMap)):
            return self._convert_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MetaboliteSpectralMatching:
    """
    Cython implementation of _MetaboliteSpectralMatching

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboliteSpectralMatching.html>`_
      -- Inherits from ['ProgressLogger', 'DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboliteSpectralMatching rv = MetaboliteSpectralMatching.__new__(MetaboliteSpectralMatching)
       rv.inst = shared_ptr[_MetaboliteSpectralMatching](new _MetaboliteSpectralMatching(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboliteSpectralMatching rv = MetaboliteSpectralMatching.__new__(MetaboliteSpectralMatching)
       rv.inst = shared_ptr[_MetaboliteSpectralMatching](new _MetaboliteSpectralMatching(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaboliteSpectralMatching](new _MetaboliteSpectralMatching())
    
    def _init_1(self, MetaboliteSpectralMatching in_0 ):
        """
        _init_1(self, in_0: MetaboliteSpectralMatching ) -> None
        """
        assert isinstance(in_0, MetaboliteSpectralMatching), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboliteSpectralMatching](new _MetaboliteSpectralMatching((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboliteSpectralMatching ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboliteSpectralMatching)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment exp , MSExperiment speclib , MzTab mz_tab ,  out_spectra ):
        """
        run(self, exp: MSExperiment , speclib: MSExperiment , mz_tab: MzTab , out_spectra: String ) -> None
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(speclib, MSExperiment), 'arg speclib wrong type'
        assert isinstance(mz_tab, MzTab), 'arg mz_tab wrong type'
        assert isinstance(out_spectra, String), 'arg out_spectra wrong type'
    
    
    
    
        self.inst.get().run((deref(exp.inst.get())), (deref(speclib.inst.get())), (deref(mz_tab.inst.get())), deref((<String>out_spectra).inst.get()))
    
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
    computeHyperScore = __static_MetaboliteSpectralMatching_computeHyperScore 

cdef class ModificationDefinitionsSet:
    """
    Cython implementation of _ModificationDefinitionsSet

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModificationDefinitionsSet.html>`_

    Representation of a set of modification definitions
    
    This class enhances the modification definitions as defined in the
    class ModificationDefinition into a set of definitions. This is also
    e.g. used as input parameters in search engines.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ModificationDefinitionsSet rv = ModificationDefinitionsSet.__new__(ModificationDefinitionsSet)
       rv.inst = shared_ptr[_ModificationDefinitionsSet](new _ModificationDefinitionsSet(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ModificationDefinitionsSet rv = ModificationDefinitionsSet.__new__(ModificationDefinitionsSet)
       rv.inst = shared_ptr[_ModificationDefinitionsSet](new _ModificationDefinitionsSet(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ModificationDefinitionsSet](new _ModificationDefinitionsSet())
    
    def _init_1(self, ModificationDefinitionsSet in_0 ):
        """
        _init_1(self, in_0: ModificationDefinitionsSet ) -> None
        """
        assert isinstance(in_0, ModificationDefinitionsSet), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ModificationDefinitionsSet](new _ModificationDefinitionsSet((deref(in_0.inst.get()))))
    
    def _init_2(self, list fixed_modifications , list variable_modifications ):
        """
        _init_2(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
        """
        assert isinstance(fixed_modifications, list) and all(isinstance(li, bytes) for li in fixed_modifications), 'arg fixed_modifications wrong type'
        assert isinstance(variable_modifications, list) and all(isinstance(li, bytes) for li in variable_modifications), 'arg variable_modifications wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in fixed_modifications:
           v0.push_back(_String(<char *>item0))
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in variable_modifications:
           v1.push_back(_String(<char *>item1))
        self.inst = shared_ptr[_ModificationDefinitionsSet](new _ModificationDefinitionsSet(deref(v0), deref(v1)))
        del v1
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ModificationDefinitionsSet ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ModificationDefinitionsSet)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], list) and all(isinstance(li, bytes) for li in args[1])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMaxModifications(self,  max_mod ):
        """
        setMaxModifications(self, max_mod: int ) -> None
        Sets the maximal number of modifications allowed per peptide
        """
        assert isinstance(max_mod, int) and max_mod >= 0, 'arg max_mod wrong type'
    
        self.inst.get().setMaxModifications((<size_t>max_mod))
    
    def getMaxModifications(self):
        """
        getMaxModifications(self) -> int
        Return the maximal number of modifications allowed per peptide
        """
        cdef size_t _r = self.inst.get().getMaxModifications()
        py_result = <size_t>_r
        return py_result
    
    def getNumberOfModifications(self):
        """
        getNumberOfModifications(self) -> int
        Returns the number of modifications stored in this set
        """
        cdef size_t _r = self.inst.get().getNumberOfModifications()
        py_result = <size_t>_r
        return py_result
    
    def getNumberOfFixedModifications(self):
        """
        getNumberOfFixedModifications(self) -> int
        Returns the number of fixed modifications stored in this set
        """
        cdef size_t _r = self.inst.get().getNumberOfFixedModifications()
        py_result = <size_t>_r
        return py_result
    
    def getNumberOfVariableModifications(self):
        """
        getNumberOfVariableModifications(self) -> int
        Returns the number of variable modifications stored in this set
        """
        cdef size_t _r = self.inst.get().getNumberOfVariableModifications()
        py_result = <size_t>_r
        return py_result
    
    def addModification(self, ModificationDefinition mod_def ):
        """
        addModification(self, mod_def: ModificationDefinition ) -> None
        Adds a modification definition to the set
        """
        assert isinstance(mod_def, ModificationDefinition), 'arg mod_def wrong type'
    
        self.inst.get().addModification((deref(mod_def.inst.get())))
    
    def _setModifications_0(self, set mod_defs ):
        """
        _setModifications_0(self, mod_defs: Set[ModificationDefinition] ) -> None
        Sets the modification definitions
        """
        assert isinstance(mod_defs, set) and all(isinstance(li, ModificationDefinition) for li in mod_defs), 'arg mod_defs wrong type'
        cdef libcpp_set[_ModificationDefinition] * v0 = new libcpp_set[_ModificationDefinition]()
        cdef ModificationDefinition item0
        for item0 in mod_defs:
           v0.insert(deref(item0.inst.get()))
        self.inst.get().setModifications(deref(v0))
        replace = set()
        cdef libcpp_set[_ModificationDefinition].iterator it_mod_defs = v0.begin()
        while it_mod_defs != v0.end():
           item0 = ModificationDefinition.__new__(ModificationDefinition)
           item0.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(it_mod_defs)))
           replace.add(item0)
           inc(it_mod_defs)
        mod_defs.clear()
        mod_defs.update(replace)
        del v0
    
    def _setModifications_1(self,  fixed_modifications ,  variable_modifications ):
        """
        _setModifications_1(self, fixed_modifications: Union[bytes, str, String] , variable_modifications: String ) -> None
        Set the modification definitions from a string
        
        The strings should contain a comma separated list of modifications. The names
        can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
        definitions and other specific definitions are given by the modifications themselves.
        """
        assert (isinstance(fixed_modifications, str) or isinstance(fixed_modifications, bytes) or isinstance(fixed_modifications, String)), 'arg fixed_modifications wrong type'
        assert isinstance(variable_modifications, String), 'arg variable_modifications wrong type'
    
    
        self.inst.get().setModifications(deref((convString(fixed_modifications)).get()), deref((<String>variable_modifications).inst.get()))
    
    def _setModifications_2(self, list fixed_modifications , list variable_modifications ):
        """
        _setModifications_2(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
        Same as above, but using StringList instead of comma separated strings
        """
        assert isinstance(fixed_modifications, list) and all(isinstance(li, bytes) for li in fixed_modifications), 'arg fixed_modifications wrong type'
        assert isinstance(variable_modifications, list) and all(isinstance(li, bytes) for li in variable_modifications), 'arg variable_modifications wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in fixed_modifications:
           v0.push_back(_String(<char *>item0))
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in variable_modifications:
           v1.push_back(_String(<char *>item1))
        self.inst.get().setModifications(deref(v0), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        variable_modifications[:] = replace
        del v1
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        fixed_modifications[:] = replace
        del v0
    
    def setModifications(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setModifications(self, mod_defs: Set[ModificationDefinition] ) -> None
          :noindex:
        
        Sets the modification definitions

        
        .. rubric:: Overload:
        .. py:function:: setModifications(self, fixed_modifications: Union[bytes, str, String] , variable_modifications: String ) -> None
          :noindex:
        
        Set the modification definitions from a string
        
        The strings should contain a comma separated list of modifications. The names
        can be PSI-MOD identifier or any other unique name supported by PSI-MOD. TermSpec
        definitions and other specific definitions are given by the modifications themselves.
        
        .. rubric:: Overload:
        .. py:function:: setModifications(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
          :noindex:
        
        Same as above, but using StringList instead of comma separated strings
    
        """
        if (len(args)==1) and (isinstance(args[0], set) and all(isinstance(li, ModificationDefinition) for li in args[0])):
            return self._setModifications_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], String)):
            return self._setModifications_1(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], list) and all(isinstance(li, bytes) for li in args[1])):
            return self._setModifications_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getModifications(self):
        """
        getModifications(self) -> Set[ModificationDefinition]
        Returns the stored modification definitions
        """
        _r = self.inst.get().getModifications()
        py_result = set()
        cdef libcpp_set[_ModificationDefinition].iterator it__r = _r.begin()
        cdef ModificationDefinition item_py_result
        while it__r != _r.end():
           item_py_result = ModificationDefinition.__new__(ModificationDefinition)
           item_py_result.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(it__r)))
           py_result.add(item_py_result)
           inc(it__r)
        return py_result
    
    def getFixedModifications(self):
        """
        getFixedModifications(self) -> Set[ModificationDefinition]
        Returns the stored fixed modification definitions
        """
        _r = self.inst.get().getFixedModifications()
        py_result = set()
        cdef libcpp_set[_ModificationDefinition].iterator it__r = _r.begin()
        cdef ModificationDefinition item_py_result
        while it__r != _r.end():
           item_py_result = ModificationDefinition.__new__(ModificationDefinition)
           item_py_result.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(it__r)))
           py_result.add(item_py_result)
           inc(it__r)
        return py_result
    
    def getVariableModifications(self):
        """
        getVariableModifications(self) -> Set[ModificationDefinition]
        Returns the stored variable modification definitions
        """
        _r = self.inst.get().getVariableModifications()
        py_result = set()
        cdef libcpp_set[_ModificationDefinition].iterator it__r = _r.begin()
        cdef ModificationDefinition item_py_result
        while it__r != _r.end():
           item_py_result = ModificationDefinition.__new__(ModificationDefinition)
           item_py_result.inst = shared_ptr[_ModificationDefinition](new _ModificationDefinition(deref(it__r)))
           py_result.add(item_py_result)
           inc(it__r)
        return py_result
    
    def _getModificationNames_0(self, list fixed_modifications , list variable_modifications ):
        """
        _getModificationNames_0(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
        Populates the output lists with the modification names (use e.g. for
        """
        assert isinstance(fixed_modifications, list) and all(isinstance(li, bytes) for li in fixed_modifications), 'arg fixed_modifications wrong type'
        assert isinstance(variable_modifications, list) and all(isinstance(li, bytes) for li in variable_modifications), 'arg variable_modifications wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in fixed_modifications:
           v0.push_back(_String(<char *>item0))
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in variable_modifications:
           v1.push_back(_String(<char *>item1))
        self.inst.get().getModificationNames(deref(v0), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        variable_modifications[:] = replace
        del v1
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        fixed_modifications[:] = replace
        del v0
    
    def _getModificationNames_1(self):
        """
        _getModificationNames_1(self) -> Set[bytes]
        Returns only the names of the modifications stored in the set
        """
        _r = self.inst.get().getModificationNames()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getModificationNames(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getModificationNames(self, fixed_modifications: List[bytes] , variable_modifications: List[bytes] ) -> None
          :noindex:
        
        Populates the output lists with the modification names (use e.g. for
        
        .. rubric:: Overload:
        .. py:function:: getModificationNames(self, ) -> Set[bytes]
          :noindex:
        
        Returns only the names of the modifications stored in the set
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(li, bytes) for li in args[0])) and (isinstance(args[1], list) and all(isinstance(li, bytes) for li in args[1])):
            return self._getModificationNames_0(*args)
        elif not args:
            return self._getModificationNames_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getFixedModificationNames(self):
        """
        getFixedModificationNames(self) -> Set[bytes]
        Returns only the names of the fixed modifications
        """
        _r = self.inst.get().getFixedModificationNames()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getVariableModificationNames(self):
        """
        getVariableModificationNames(self) -> Set[bytes]
        Returns only the names of the variable modifications
        """
        _r = self.inst.get().getVariableModificationNames()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def isCompatible(self, AASequence peptide ):
        """
        isCompatible(self, peptide: AASequence ) -> bool
        Returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications
        """
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
    
        cdef bool _r = self.inst.get().isCompatible((deref(peptide.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def inferFromPeptides(self, PeptideIdentificationList peptides ):
        """
        inferFromPeptides(self, peptides: PeptideIdentificationList ) -> None
        Infers the sets of defined modifications from the modifications present on peptide identifications
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().inferFromPeptides((deref(peptides.inst.get()))) 

cdef class NoiseEstimator:
    """
    Cython implementation of _NoiseEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1NoiseEstimator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property nr_windows:
        def __set__(self,  nr_windows):
        
            self.inst.get().nr_windows = (<int>nr_windows)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().nr_windows
            py_result = <int>_r
            return py_result
    
    property mz_start:
        def __set__(self, double mz_start):
        
            self.inst.get().mz_start = (<double>mz_start)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz_start
            py_result = <double>_r
            return py_result
    
    property window_length:
        def __set__(self, double window_length):
        
            self.inst.get().window_length = (<double>window_length)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().window_length
            py_result = <double>_r
            return py_result
    
    property result_windows_even:
        def __set__(self, list result_windows_even):
            cdef libcpp_vector[double] v0 = result_windows_even
            self.inst.get().result_windows_even = v0
            
    
        def __get__(self):
            _r = self.inst.get().result_windows_even
            cdef list py_result = _r
            return py_result
    
    property result_windows_odd:
        def __set__(self, list result_windows_odd):
            cdef libcpp_vector[double] v0 = result_windows_odd
            self.inst.get().result_windows_odd = v0
            
    
        def __get__(self):
            _r = self.inst.get().result_windows_odd
            cdef list py_result = _r
            return py_result
    
    def __copy__(self):
       cdef NoiseEstimator rv = NoiseEstimator.__new__(NoiseEstimator)
       rv.inst = shared_ptr[_NoiseEstimator](new _NoiseEstimator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef NoiseEstimator rv = NoiseEstimator.__new__(NoiseEstimator)
       rv.inst = shared_ptr[_NoiseEstimator](new _NoiseEstimator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_NoiseEstimator](new _NoiseEstimator())
    
    def _init_1(self, NoiseEstimator in_0 ):
        """
        _init_1(self, in_0: NoiseEstimator ) -> None
        """
        assert isinstance(in_0, NoiseEstimator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_NoiseEstimator](new _NoiseEstimator((deref(in_0.inst.get()))))
    
    def _init_2(self, double nr_windows_ , double mz_start_ , double win_len_ ):
        """
        _init_2(self, nr_windows_: float , mz_start_: float , win_len_: float ) -> None
        """
        assert isinstance(nr_windows_, float), 'arg nr_windows_ wrong type'
        assert isinstance(mz_start_, float), 'arg mz_start_ wrong type'
        assert isinstance(win_len_, float), 'arg win_len_ wrong type'
    
    
    
        self.inst = shared_ptr[_NoiseEstimator](new _NoiseEstimator((<double>nr_windows_), (<double>mz_start_), (<double>win_len_)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: NoiseEstimator ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, nr_windows_: float , mz_start_: float , win_len_: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], NoiseEstimator)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def get_noise_value(self, double mz ):
        """
        get_noise_value(self, mz: float ) -> float
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef double _r = self.inst.get().get_noise_value((<double>mz))
        py_result = <double>_r
        return py_result
    
    def get_noise_even(self, double mz ):
        """
        get_noise_even(self, mz: float ) -> float
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef double _r = self.inst.get().get_noise_even((<double>mz))
        py_result = <double>_r
        return py_result
    
    def get_noise_odd(self, double mz ):
        """
        get_noise_odd(self, mz: float ) -> float
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef double _r = self.inst.get().get_noise_odd((<double>mz))
        py_result = <double>_r
        return py_result 

cdef class SignalToNoiseEstimatorMedianRapid:
    """
    Cython implementation of _SignalToNoiseEstimatorMedianRapid

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMedianRapid.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SignalToNoiseEstimatorMedianRapid rv = SignalToNoiseEstimatorMedianRapid.__new__(SignalToNoiseEstimatorMedianRapid)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMedianRapid](new _SignalToNoiseEstimatorMedianRapid(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SignalToNoiseEstimatorMedianRapid rv = SignalToNoiseEstimatorMedianRapid.__new__(SignalToNoiseEstimatorMedianRapid)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMedianRapid](new _SignalToNoiseEstimatorMedianRapid(deref(self.inst.get())))
       return rv
    
    def _init_0(self, SignalToNoiseEstimatorMedianRapid in_0 ):
        """
        _init_0(self, in_0: SignalToNoiseEstimatorMedianRapid ) -> None
        """
        assert isinstance(in_0, SignalToNoiseEstimatorMedianRapid), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedianRapid](new _SignalToNoiseEstimatorMedianRapid((deref(in_0.inst.get()))))
    
    def _init_1(self, double window_length ):
        """
        _init_1(self, window_length: float ) -> None
        """
        assert isinstance(window_length, float), 'arg window_length wrong type'
    
        self.inst = shared_ptr[_SignalToNoiseEstimatorMedianRapid](new _SignalToNoiseEstimatorMedianRapid((<double>window_length)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SignalToNoiseEstimatorMedianRapid ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, window_length: float ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], SignalToNoiseEstimatorMedianRapid)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], float)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _estimateNoise_0(self, _Interfaces_Spectrum in_0 ):
        """
        _estimateNoise_0(self, in_0: _Interfaces_Spectrum ) -> NoiseEstimator
        """
        assert isinstance(in_0, _Interfaces_Spectrum), 'arg in_0 wrong type'
        cdef shared_ptr[_Spectrum] input_in_0 = in_0.inst
        cdef _NoiseEstimator * _r = new _NoiseEstimator(self.inst.get().estimateNoise(input_in_0))
        cdef NoiseEstimator py_result = NoiseEstimator.__new__(NoiseEstimator)
        py_result.inst = shared_ptr[_NoiseEstimator](_r)
        return py_result
    
    def _estimateNoise_1(self, _Interfaces_Chromatogram in_0 ):
        """
        _estimateNoise_1(self, in_0: _Interfaces_Chromatogram ) -> NoiseEstimator
        """
        assert isinstance(in_0, _Interfaces_Chromatogram), 'arg in_0 wrong type'
        cdef shared_ptr[_Chromatogram] input_in_0 = in_0.inst
        cdef _NoiseEstimator * _r = new _NoiseEstimator(self.inst.get().estimateNoise(input_in_0))
        cdef NoiseEstimator py_result = NoiseEstimator.__new__(NoiseEstimator)
        py_result.inst = shared_ptr[_NoiseEstimator](_r)
        return py_result
    
    def _estimateNoise_2(self, list mz_array , list int_array ):
        """
        _estimateNoise_2(self, mz_array: List[float] , int_array: List[float] ) -> NoiseEstimator
        """
        assert isinstance(mz_array, list) and all(isinstance(elemt_rec, float) for elemt_rec in mz_array), 'arg mz_array wrong type'
        assert isinstance(int_array, list) and all(isinstance(elemt_rec, float) for elemt_rec in int_array), 'arg int_array wrong type'
        cdef libcpp_vector[double] v0 = mz_array
        cdef libcpp_vector[double] v1 = int_array
        cdef _NoiseEstimator * _r = new _NoiseEstimator(self.inst.get().estimateNoise(v0, v1))
        
        
        cdef NoiseEstimator py_result = NoiseEstimator.__new__(NoiseEstimator)
        py_result.inst = shared_ptr[_NoiseEstimator](_r)
        return py_result
    
    def estimateNoise(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: estimateNoise(self, in_0: _Interfaces_Spectrum ) -> NoiseEstimator
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: estimateNoise(self, in_0: _Interfaces_Chromatogram ) -> NoiseEstimator
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: estimateNoise(self, mz_array: List[float] , int_array: List[float] ) -> NoiseEstimator
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], _Interfaces_Spectrum)):
            return self._estimateNoise_0(*args)
        elif (len(args)==1) and (isinstance(args[0], _Interfaces_Chromatogram)):
            return self._estimateNoise_1(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])):
            return self._estimateNoise_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SpectraSTSimilarityScore:
    """
    Cython implementation of _SpectraSTSimilarityScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectraSTSimilarityScore.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectraSTSimilarityScore rv = SpectraSTSimilarityScore.__new__(SpectraSTSimilarityScore)
       rv.inst = shared_ptr[_SpectraSTSimilarityScore](new _SpectraSTSimilarityScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectraSTSimilarityScore rv = SpectraSTSimilarityScore.__new__(SpectraSTSimilarityScore)
       rv.inst = shared_ptr[_SpectraSTSimilarityScore](new _SpectraSTSimilarityScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectraSTSimilarityScore](new _SpectraSTSimilarityScore())
    
    def _init_1(self, SpectraSTSimilarityScore in_0 ):
        """
        _init_1(self, in_0: SpectraSTSimilarityScore ) -> None
        """
        assert isinstance(in_0, SpectraSTSimilarityScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectraSTSimilarityScore](new _SpectraSTSimilarityScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectraSTSimilarityScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectraSTSimilarityScore)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def preprocess(self, MSSpectrum spec , float remove_peak_intensity_threshold ,  cut_peaks_below ,  min_peak_number ,  max_peak_number ):
        """
        preprocess(self, spec: MSSpectrum , remove_peak_intensity_threshold: float , cut_peaks_below: int , min_peak_number: int , max_peak_number: int ) -> bool
        Preprocesses the spectrum
        
        The preprocessing removes peak below a intensity threshold, reject spectra that does
        not have enough peaks, and cuts peaks exceeding the max_peak_number most intense peaks
        
        :returns: true if spectrum passes filtering
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(remove_peak_intensity_threshold, float), 'arg remove_peak_intensity_threshold wrong type'
        assert isinstance(cut_peaks_below, int), 'arg cut_peaks_below wrong type'
        assert isinstance(min_peak_number, int) and min_peak_number >= 0, 'arg min_peak_number wrong type'
        assert isinstance(max_peak_number, int) and max_peak_number >= 0, 'arg max_peak_number wrong type'
    
    
    
    
    
        cdef bool _r = self.inst.get().preprocess((deref(spec.inst.get())), (<float>remove_peak_intensity_threshold), (<unsigned int>cut_peaks_below), (<size_t>min_peak_number), (<size_t>max_peak_number))
        py_result = <bool>_r
        return py_result
    
    def transform(self, MSSpectrum spec ):
        """
        transform(self, spec: MSSpectrum ) -> BinnedSpectrum
        Spectrum is transformed into a binned spectrum with bin size 1 and spread 1 and the intensities are normalized
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        cdef _BinnedSpectrum * _r = new _BinnedSpectrum(self.inst.get().transform((deref(spec.inst.get()))))
        cdef BinnedSpectrum py_result = BinnedSpectrum.__new__(BinnedSpectrum)
        py_result.inst = shared_ptr[_BinnedSpectrum](_r)
        return py_result
    
    def dot_bias(self, BinnedSpectrum bin1 , BinnedSpectrum bin2 , double dot_product ):
        """
        dot_bias(self, bin1: BinnedSpectrum , bin2: BinnedSpectrum , dot_product: float ) -> float
        Calculates how much of the dot product is dominated by a few peaks
        
        :param dot_product: If -1 this value will be calculated as well.
        :param bin1: First spectrum in binned representation
        :param bin2: Second spectrum in binned representation
        """
        assert isinstance(bin1, BinnedSpectrum), 'arg bin1 wrong type'
        assert isinstance(bin2, BinnedSpectrum), 'arg bin2 wrong type'
        assert isinstance(dot_product, float), 'arg dot_product wrong type'
    
    
    
        cdef double _r = self.inst.get().dot_bias((deref(bin1.inst.get())), (deref(bin2.inst.get())), (<double>dot_product))
        py_result = <double>_r
        return py_result
    
    def delta_D(self, double top_hit , double runner_up ):
        """
        delta_D(self, top_hit: float , runner_up: float ) -> float
        Calculates the normalized distance between top_hit and runner_up
        
        :param top_hit: Is the best score for a given match
        :param runner_up: A match with a worse score than top_hit, e.g. the second best score
        :returns: normalized distance
        """
        assert isinstance(top_hit, float), 'arg top_hit wrong type'
        assert isinstance(runner_up, float), 'arg runner_up wrong type'
    
    
        cdef double _r = self.inst.get().delta_D((<double>top_hit), (<double>runner_up))
        py_result = <double>_r
        return py_result
    
    def compute_F(self, double dot_product , double delta_D , double dot_bias ):
        """
        compute_F(self, dot_product: float , delta_D: float , dot_bias: float ) -> float
        Computes the overall all score
        
        :param dot_product: dot_product of a match
        :param delta_D: delta_D should be calculated after all dot products for a unidentified spectrum are computed
        :param dot_bias: the bias
        :returns: The SpectraST similarity score
        """
        assert isinstance(dot_product, float), 'arg dot_product wrong type'
        assert isinstance(delta_D, float), 'arg delta_D wrong type'
        assert isinstance(dot_bias, float), 'arg dot_bias wrong type'
    
    
    
        cdef double _r = self.inst.get().compute_F((<double>dot_product), (<double>delta_D), (<double>dot_bias))
        py_result = <double>_r
        return py_result 

cdef class SpectralMatch:
    """
    Cython implementation of _SpectralMatch

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectralMatch.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectralMatch rv = SpectralMatch.__new__(SpectralMatch)
       rv.inst = shared_ptr[_SpectralMatch](new _SpectralMatch(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectralMatch rv = SpectralMatch.__new__(SpectralMatch)
       rv.inst = shared_ptr[_SpectralMatch](new _SpectralMatch(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectralMatch](new _SpectralMatch())
    
    def _init_1(self, SpectralMatch in_0 ):
        """
        _init_1(self, in_0: SpectralMatch ) -> None
        """
        assert isinstance(in_0, SpectralMatch), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectralMatch](new _SpectralMatch((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectralMatch ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectralMatch)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getObservedPrecursorMass(self):
        """
        getObservedPrecursorMass(self) -> float
        """
        cdef double _r = self.inst.get().getObservedPrecursorMass()
        py_result = <double>_r
        return py_result
    
    def setObservedPrecursorMass(self, double in_0 ):
        """
        setObservedPrecursorMass(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setObservedPrecursorMass((<double>in_0))
    
    def getObservedPrecursorRT(self):
        """
        getObservedPrecursorRT(self) -> float
        """
        cdef double _r = self.inst.get().getObservedPrecursorRT()
        py_result = <double>_r
        return py_result
    
    def setObservedPrecursorRT(self, double in_0 ):
        """
        setObservedPrecursorRT(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setObservedPrecursorRT((<double>in_0))
    
    def getFoundPrecursorMass(self):
        """
        getFoundPrecursorMass(self) -> float
        """
        cdef double _r = self.inst.get().getFoundPrecursorMass()
        py_result = <double>_r
        return py_result
    
    def setFoundPrecursorMass(self, double in_0 ):
        """
        setFoundPrecursorMass(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setFoundPrecursorMass((<double>in_0))
    
    def getFoundPrecursorCharge(self):
        """
        getFoundPrecursorCharge(self) -> int
        """
        cdef int _r = self.inst.get().getFoundPrecursorCharge()
        py_result = <int>_r
        return py_result
    
    def setFoundPrecursorCharge(self,  in_0 ):
        """
        setFoundPrecursorCharge(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().setFoundPrecursorCharge((<int>in_0))
    
    def getMatchingScore(self):
        """
        getMatchingScore(self) -> float
        """
        cdef double _r = self.inst.get().getMatchingScore()
        py_result = <double>_r
        return py_result
    
    def setMatchingScore(self, double in_0 ):
        """
        setMatchingScore(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setMatchingScore((<double>in_0))
    
    def getObservedSpectrumIndex(self):
        """
        getObservedSpectrumIndex(self) -> int
        """
        cdef size_t _r = self.inst.get().getObservedSpectrumIndex()
        py_result = <size_t>_r
        return py_result
    
    def setObservedSpectrumIndex(self,  in_0 ):
        """
        setObservedSpectrumIndex(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        self.inst.get().setObservedSpectrumIndex((<size_t>in_0))
    
    def getMatchingSpectrumIndex(self):
        """
        getMatchingSpectrumIndex(self) -> int
        """
        cdef size_t _r = self.inst.get().getMatchingSpectrumIndex()
        py_result = <size_t>_r
        return py_result
    
    def setMatchingSpectrumIndex(self,  in_0 ):
        """
        setMatchingSpectrumIndex(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        self.inst.get().setMatchingSpectrumIndex((<size_t>in_0))
    
    def getPrimaryIdentifier(self):
        """
        getPrimaryIdentifier(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getPrimaryIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setPrimaryIdentifier(self,  in_0 ):
        """
        setPrimaryIdentifier(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setPrimaryIdentifier(deref((convString(in_0)).get()))
    
    def getSecondaryIdentifier(self):
        """
        getSecondaryIdentifier(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getSecondaryIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setSecondaryIdentifier(self,  in_0 ):
        """
        setSecondaryIdentifier(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setSecondaryIdentifier(deref((convString(in_0)).get()))
    
    def getCommonName(self):
        """
        getCommonName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getCommonName()
        py_result = convOutputString(_r)
        return py_result
    
    def setCommonName(self,  in_0 ):
        """
        setCommonName(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setCommonName(deref((convString(in_0)).get()))
    
    def getSumFormula(self):
        """
        getSumFormula(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getSumFormula()
        py_result = convOutputString(_r)
        return py_result
    
    def setSumFormula(self,  in_0 ):
        """
        setSumFormula(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setSumFormula(deref((convString(in_0)).get()))
    
    def getInchiString(self):
        """
        getInchiString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getInchiString()
        py_result = convOutputString(_r)
        return py_result
    
    def setInchiString(self,  in_0 ):
        """
        setInchiString(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setInchiString(deref((convString(in_0)).get()))
    
    def getSMILESString(self):
        """
        getSMILESString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getSMILESString()
        py_result = convOutputString(_r)
        return py_result
    
    def setSMILESString(self,  in_0 ):
        """
        setSMILESString(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setSMILESString(deref((convString(in_0)).get()))
    
    def getPrecursorAdduct(self):
        """
        getPrecursorAdduct(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getPrecursorAdduct()
        py_result = convOutputString(_r)
        return py_result
    
    def setPrecursorAdduct(self,  in_0 ):
        """
        setPrecursorAdduct(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setPrecursorAdduct(deref((convString(in_0)).get())) 

cdef class SqMassConfig:
    """
    Cython implementation of _SqMassConfig

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqMassConfig.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property write_full_meta:
        def __set__(self, bool write_full_meta):
        
            self.inst.get().write_full_meta = (<bool>write_full_meta)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().write_full_meta
            py_result = <bool>_r
            return py_result
    
    property use_lossy_numpress:
        def __set__(self, bool use_lossy_numpress):
        
            self.inst.get().use_lossy_numpress = (<bool>use_lossy_numpress)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_lossy_numpress
            py_result = <bool>_r
            return py_result
    
    property linear_fp_mass_acc:
        def __set__(self, double linear_fp_mass_acc):
        
            self.inst.get().linear_fp_mass_acc = (<double>linear_fp_mass_acc)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().linear_fp_mass_acc
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef SqMassConfig rv = SqMassConfig.__new__(SqMassConfig)
       rv.inst = shared_ptr[_SqMassConfig](new _SqMassConfig(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SqMassConfig rv = SqMassConfig.__new__(SqMassConfig)
       rv.inst = shared_ptr[_SqMassConfig](new _SqMassConfig(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SqMassConfig](new _SqMassConfig())
    
    def _init_1(self, SqMassConfig in_0 ):
        """
        _init_1(self, in_0: SqMassConfig ) -> None
        """
        assert isinstance(in_0, SqMassConfig), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SqMassConfig](new _SqMassConfig((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SqMassConfig ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SqMassConfig)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SqMassFile:
    """
    Cython implementation of _SqMassFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqMassFile.html>`_

    An class that uses on-disk SQLite database to read and write spectra and chromatograms
    
    This class provides functions to read and write spectra and chromatograms
    to disk using a SQLite database and store them in sqMass format. This
    allows users to access, select and filter spectra and chromatograms
    on-demand even in a large collection of data
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SqMassFile rv = SqMassFile.__new__(SqMassFile)
       rv.inst = shared_ptr[_SqMassFile](new _SqMassFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SqMassFile rv = SqMassFile.__new__(SqMassFile)
       rv.inst = shared_ptr[_SqMassFile](new _SqMassFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SqMassFile](new _SqMassFile())
    
    def _init_1(self, SqMassFile in_0 ):
        """
        _init_1(self, in_0: SqMassFile ) -> None
        """
        assert isinstance(in_0, SqMassFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SqMassFile](new _SqMassFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SqMassFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SqMassFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment map_ ):
        """
        load(self, filename: Union[bytes, str, String] , map_: MSExperiment ) -> None
        Read / Write a complete mass spectrometric experiment
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map_, MSExperiment), 'arg map_ wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(map_.inst.get())))
    
    def store(self,  filename , MSExperiment map_ ):
        """
        store(self, filename: Union[bytes, str, String] , map_: MSExperiment ) -> None
        Store an MSExperiment in sqMass format
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(map_, MSExperiment), 'arg map_ wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(map_.inst.get())))
    
    def setConfig(self, SqMassConfig config ):
        """
        setConfig(self, config: SqMassConfig ) -> None
        """
        assert isinstance(config, SqMassConfig), 'arg config wrong type'
    
        self.inst.get().setConfig((deref(config.inst.get()))) 
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
