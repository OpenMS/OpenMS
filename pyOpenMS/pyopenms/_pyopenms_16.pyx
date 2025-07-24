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
def __static_Deisotoper_deisotopeAndSingleCharge(MSSpectrum spectra , double fragment_tolerance , bool fragment_unit_ppm ,  min_charge ,  max_charge , bool keep_only_deisotoped ,  min_isopeaks ,  max_isopeaks , bool make_single_charged , bool annotate_charge , bool annotate_iso_peak_count , bool use_decreasing_model ,  start_intensity_check , bool add_up_intensity , bool annotate_features ):
    """
    __static_Deisotoper_deisotopeAndSingleCharge(spectra: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool , min_charge: int , max_charge: int , keep_only_deisotoped: bool , min_isopeaks: int , max_isopeaks: int , make_single_charged: bool , annotate_charge: bool , annotate_iso_peak_count: bool , use_decreasing_model: bool , start_intensity_check: int , add_up_intensity: bool , annotate_features: bool ) -> None
    """
    assert isinstance(spectra, MSSpectrum), 'arg spectra wrong type'
    assert isinstance(fragment_tolerance, float), 'arg fragment_tolerance wrong type'
    assert isinstance(fragment_unit_ppm, pybool_t), 'arg fragment_unit_ppm wrong type'
    assert isinstance(min_charge, int), 'arg min_charge wrong type'
    assert isinstance(max_charge, int), 'arg max_charge wrong type'
    assert isinstance(keep_only_deisotoped, pybool_t), 'arg keep_only_deisotoped wrong type'
    assert isinstance(min_isopeaks, int), 'arg min_isopeaks wrong type'
    assert isinstance(max_isopeaks, int), 'arg max_isopeaks wrong type'
    assert isinstance(make_single_charged, pybool_t), 'arg make_single_charged wrong type'
    assert isinstance(annotate_charge, pybool_t), 'arg annotate_charge wrong type'
    assert isinstance(annotate_iso_peak_count, pybool_t), 'arg annotate_iso_peak_count wrong type'
    assert isinstance(use_decreasing_model, pybool_t), 'arg use_decreasing_model wrong type'
    assert isinstance(start_intensity_check, int), 'arg start_intensity_check wrong type'
    assert isinstance(add_up_intensity, pybool_t), 'arg add_up_intensity wrong type'
    assert isinstance(annotate_features, pybool_t), 'arg annotate_features wrong type'















    _deisotopeAndSingleCharge_Deisotoper((deref(spectra.inst.get())), (<double>fragment_tolerance), (<bool>fragment_unit_ppm), (<int>min_charge), (<int>max_charge), (<bool>keep_only_deisotoped), (<unsigned int>min_isopeaks), (<unsigned int>max_isopeaks), (<bool>make_single_charged), (<bool>annotate_charge), (<bool>annotate_iso_peak_count), (<bool>use_decreasing_model), (<unsigned int>start_intensity_check), (<bool>add_up_intensity), (<bool>annotate_features))

def __static_Deisotoper_deisotopeAndSingleChargeDefault(MSSpectrum spectra , double fragment_tolerance , bool fragment_unit_ppm ):
    """
    __static_Deisotoper_deisotopeAndSingleChargeDefault(spectra: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool ) -> None
    """
    assert isinstance(spectra, MSSpectrum), 'arg spectra wrong type'
    assert isinstance(fragment_tolerance, float), 'arg fragment_tolerance wrong type'
    assert isinstance(fragment_unit_ppm, pybool_t), 'arg fragment_unit_ppm wrong type'



    _deisotopeAndSingleCharge_Deisotoper((deref(spectra.inst.get())), (<double>fragment_tolerance), (<bool>fragment_unit_ppm))

def __static_Deisotoper_deisotopeWithAveragineModel(MSSpectrum spectrum , double fragment_tolerance , bool fragment_unit_ppm ,  number_of_final_peaks ,  min_charge ,  max_charge , bool keep_only_deisotoped ,  min_isopeaks ,  max_isopeaks , bool make_single_charged , bool annotate_charge , bool annotate_iso_peak_count , bool add_up_intensity ):
    """
    __static_Deisotoper_deisotopeWithAveragineModel(spectrum: MSSpectrum , fragment_tolerance: float , fragment_unit_ppm: bool , number_of_final_peaks: int , min_charge: int , max_charge: int , keep_only_deisotoped: bool , min_isopeaks: int , max_isopeaks: int , make_single_charged: bool , annotate_charge: bool , annotate_iso_peak_count: bool , add_up_intensity: bool ) -> None
    """
    assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    assert isinstance(fragment_tolerance, float), 'arg fragment_tolerance wrong type'
    assert isinstance(fragment_unit_ppm, pybool_t), 'arg fragment_unit_ppm wrong type'
    assert isinstance(number_of_final_peaks, int), 'arg number_of_final_peaks wrong type'
    assert isinstance(min_charge, int), 'arg min_charge wrong type'
    assert isinstance(max_charge, int), 'arg max_charge wrong type'
    assert isinstance(keep_only_deisotoped, pybool_t), 'arg keep_only_deisotoped wrong type'
    assert isinstance(min_isopeaks, int), 'arg min_isopeaks wrong type'
    assert isinstance(max_isopeaks, int), 'arg max_isopeaks wrong type'
    assert isinstance(make_single_charged, pybool_t), 'arg make_single_charged wrong type'
    assert isinstance(annotate_charge, pybool_t), 'arg annotate_charge wrong type'
    assert isinstance(annotate_iso_peak_count, pybool_t), 'arg annotate_iso_peak_count wrong type'
    assert isinstance(add_up_intensity, pybool_t), 'arg add_up_intensity wrong type'













    _deisotopeWithAveragineModel_Deisotoper((deref(spectrum.inst.get())), (<double>fragment_tolerance), (<bool>fragment_unit_ppm), (<int>number_of_final_peaks), (<int>min_charge), (<int>max_charge), (<bool>keep_only_deisotoped), (<unsigned int>min_isopeaks), (<unsigned int>max_isopeaks), (<bool>make_single_charged), (<bool>annotate_charge), (<bool>annotate_iso_peak_count), (<bool>add_up_intensity)) 

cdef class __CHARGEMODE_FD:
    None
    QFROMFEATURE = 0
    QHEURISTIC = 1
    QALL = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class DimensionDescription:
    None
    RT = 0
    MZ = 1
    DIMENSION = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __SourceClassification:
    None
    ARTIFACT = 0
    HYPOTHETICAL = 1
    NATURAL = 2
    POSTTRANSLATIONAL = 3
    MULTIPLE = 4
    CHEMICAL_DERIVATIVE = 5
    ISOTOPIC_LABEL = 6
    PRETRANSLATIONAL = 7
    OTHER_GLYCOSYLATION = 8
    NLINKED_GLYCOSYLATION = 9
    AA_SUBSTITUTION = 10
    OTHER = 11
    NONSTANDARD_RESIDUE = 12
    COTRANSLATIONAL = 13
    OLINKED_GLYCOSYLATION = 14
    UNKNOWN = 15
    NUMBER_OF_SOURCE_CLASSIFICATIONS = 16

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __TermSpecificity:
    None
    ANYWHERE = 0
    C_TERM = 1
    N_TERM = 2
    PROTEIN_C_TERM = 3
    PROTEIN_N_TERM = 4
    NUMBER_OF_TERM_SPECIFICITY = 5

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AcquisitionInfo:
    """
    Cython implementation of _AcquisitionInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AcquisitionInfo.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AcquisitionInfo rv = AcquisitionInfo.__new__(AcquisitionInfo)
       rv.inst = shared_ptr[_AcquisitionInfo](new _AcquisitionInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AcquisitionInfo rv = AcquisitionInfo.__new__(AcquisitionInfo)
       rv.inst = shared_ptr[_AcquisitionInfo](new _AcquisitionInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AcquisitionInfo](new _AcquisitionInfo())
    
    def _init_1(self, AcquisitionInfo in_0 ):
        """
        _init_1(self, in_0: AcquisitionInfo ) -> None
        """
        assert isinstance(in_0, AcquisitionInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AcquisitionInfo](new _AcquisitionInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AcquisitionInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AcquisitionInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMethodOfCombination(self):
        """
        getMethodOfCombination(self) -> Union[bytes, str, String]
        Returns the method of combination
        """
        cdef _String _r = self.inst.get().getMethodOfCombination()
        py_result = convOutputString(_r)
        return py_result
    
    def setMethodOfCombination(self,  method ):
        """
        setMethodOfCombination(self, method: Union[bytes, str, String] ) -> None
        Sets the method of combination
        """
        assert (isinstance(method, str) or isinstance(method, bytes) or isinstance(method, String)), 'arg method wrong type'
    
        self.inst.get().setMethodOfCombination(deref((convString(method)).get()))
    
    def size(self):
        """
        size(self) -> int
        Number a Acquisition objects
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def __getitem__(self,  in_0 ):
        """
        __getitem__(self, in_0: int ) -> Acquisition
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        cdef long _idx = (<size_t>in_0)
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        cdef _Acquisition * _r = new _Acquisition(deref(self.inst.get())[(<size_t>in_0)])
        cdef Acquisition py_result = Acquisition.__new__(Acquisition)
        py_result.inst = shared_ptr[_Acquisition](_r)
        return py_result
    def __setitem__(self, key, Acquisition value):
        """Cython signature: Acquisition & operator[](size_t)"""
        assert isinstance(key, int), 'arg index wrong type'
    
        cdef long _idx = key
        if _idx < 0:
            raise IndexError("invalid index %d" % _idx)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)
        deref(self.inst.get())[key] = (deref(value.inst.get()))
    
    def push_back(self, Acquisition in_0 ):
        """
        push_back(self, in_0: Acquisition ) -> None
        Append a Acquisition object
        """
        assert isinstance(in_0, Acquisition), 'arg in_0 wrong type'
    
        self.inst.get().push_back((deref(in_0.inst.get())))
    
    def resize(self,  n ):
        """
        resize(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().resize((<size_t>n))
    
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
        if not isinstance(other, AcquisitionInfo):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef AcquisitionInfo other_casted = other
        cdef AcquisitionInfo self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class BayesianProteinInferenceAlgorithm:
    """
    Cython implementation of _BayesianProteinInferenceAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BayesianProteinInferenceAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Performs a Bayesian protein inference on Protein/Peptide identifications or ConsensusMap.
    
    - Filters for best n PSMs per spectrum.
    - Calculates and filters for best peptide per spectrum.
    - Builds a k-partite graph from the structures.
    - Finds and splits into connected components by DFS
    - Extends the graph by adding layers from indist. protein groups, peptides with the same parents and optionally
      some additional layers (peptide sequence, charge, replicate -> extended model = experimental)
    - Builds a factor graph representation of a Bayesian network using the Evergreen library
      See model param section. It is based on the Fido noisy-OR model with an option for
      regularizing the number of proteins per peptide.
    - Performs loopy belief propagation on the graph and queries protein, protein group and/or peptide posteriors
      See loopy_belief_propagation param section.
    - Learns best parameters via grid search if the parameters were not given in the param section.
    - Writes posteriors to peptides and/or proteins and adds indistinguishable protein groups to the underlying
      data structures.
    - Can make use of OpenMP to parallelize over connected components.
    
    Usage:
    
    .. code-block:: python
    
      from pyopenms import *
      from urllib.request import urlretrieve
      urlretrieve("https://raw.githubusercontent.com/OpenMS/OpenMS/develop/src/tests/class_tests/openms/data/BayesianProteinInference_test.idXML", "BayesianProteinInference_test.idXML")
      proteins = []
      peptides = []
      idf = IdXMLFile()
      idf.load("BayesianProteinInference_test.idXML", proteins, peptides)
      bpia = BayesianProteinInferenceAlgorithm()
      p = bpia.getParameters()
      p.setValue("update_PSM_probabilities", "false")
      bpia.setParameters(p)
      bpia.inferPosteriorProbabilities(proteins, peptides)
      #
      print(len(peptides)) # 9
      print(peptides[0].getHits()[0].getScore()) # 0.6
      print(proteins[0].getHits()[0].getScore()) # 0.624641
      print(proteins[0].getHits()[1].getScore()) # 0.648346
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BayesianProteinInferenceAlgorithm](new _BayesianProteinInferenceAlgorithm())
    
    def _init_1(self,  debug_lvl ):
        """
        _init_1(self, debug_lvl: int ) -> None
        """
        assert isinstance(debug_lvl, int), 'arg debug_lvl wrong type'
    
        self.inst = shared_ptr[_BayesianProteinInferenceAlgorithm](new _BayesianProteinInferenceAlgorithm((<unsigned int>debug_lvl)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, debug_lvl: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _inferPosteriorProbabilities_0(self, list proteinIDs , PeptideIdentificationList peptideIDs , bool greedy_group_resolution ):
        """
        _inferPosteriorProbabilities_0(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool ) -> None
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        """
        assert isinstance(proteinIDs, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteinIDs), 'arg proteinIDs wrong type'
        assert isinstance(peptideIDs, PeptideIdentificationList), 'arg peptideIDs wrong type'
        assert isinstance(greedy_group_resolution, pybool_t), 'arg greedy_group_resolution wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in proteinIDs:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().inferPosteriorProbabilities(deref(v0), (deref(peptideIDs.inst.get())), (<bool>greedy_group_resolution))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteinIDs = v0.begin()
        replace_0 = []
        while it_proteinIDs != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteinIDs)))
            replace_0.append(item0)
            inc(it_proteinIDs)
        proteinIDs[:] = replace_0
        del v0
    
    def _inferPosteriorProbabilities_1(self, list proteinIDs , PeptideIdentificationList peptideIDs , bool greedy_group_resolution , ExperimentalDesign exp_des ):
        """
        _inferPosteriorProbabilities_1(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None
        Writes its results into protein and (optionally also) peptide hits (as new score).
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        Experimental design can be used to create an extended graph with replicate information. (experimental)
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :param exp_des: Experimental Design
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        """
        assert isinstance(proteinIDs, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteinIDs), 'arg proteinIDs wrong type'
        assert isinstance(peptideIDs, PeptideIdentificationList), 'arg peptideIDs wrong type'
        assert isinstance(greedy_group_resolution, pybool_t), 'arg greedy_group_resolution wrong type'
        assert isinstance(exp_des, ExperimentalDesign), 'arg exp_des wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in proteinIDs:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        self.inst.get().inferPosteriorProbabilities(deref(v0), (deref(peptideIDs.inst.get())), (<bool>greedy_group_resolution), (deref(exp_des.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteinIDs = v0.begin()
        replace_0 = []
        while it_proteinIDs != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteinIDs)))
            replace_0.append(item0)
            inc(it_proteinIDs)
        proteinIDs[:] = replace_0
        del v0
    
    def _inferPosteriorProbabilities_2(self, ConsensusMap cmap , bool greedy_group_resolution ):
        """
        _inferPosteriorProbabilities_2(self, cmap: ConsensusMap , greedy_group_resolution: bool ) -> None
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :return: Writes its protein ID results into the ConsensusMap
        """
        assert isinstance(cmap, ConsensusMap), 'arg cmap wrong type'
        assert isinstance(greedy_group_resolution, pybool_t), 'arg greedy_group_resolution wrong type'
    
    
        self.inst.get().inferPosteriorProbabilities((deref(cmap.inst.get())), (<bool>greedy_group_resolution))
    
    def _inferPosteriorProbabilities_3(self, ConsensusMap cmap , bool greedy_group_resolution , ExperimentalDesign exp_des ):
        """
        _inferPosteriorProbabilities_3(self, cmap: ConsensusMap , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs.
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :param exp_des: Experimental Design
        :return: Writes its protein ID results into the ConsensusMap
        """
        assert isinstance(cmap, ConsensusMap), 'arg cmap wrong type'
        assert isinstance(greedy_group_resolution, pybool_t), 'arg greedy_group_resolution wrong type'
        assert isinstance(exp_des, ExperimentalDesign), 'arg exp_des wrong type'
    
    
    
        self.inst.get().inferPosteriorProbabilities((deref(cmap.inst.get())), (<bool>greedy_group_resolution), (deref(exp_des.inst.get())))
    
    def inferPosteriorProbabilities(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: inferPosteriorProbabilities(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool ) -> None
          :noindex:
        
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        
        .. rubric:: Overload:
        .. py:function:: inferPosteriorProbabilities(self, proteinIDs: List[ProteinIdentification] , peptideIDs: PeptideIdentificationList , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None
          :noindex:
        
        Writes its results into protein and (optionally also) peptide hits (as new score).
        Optionally adds indistinguishable protein groups with separate scores, too
        Currently only takes first proteinID run and all peptides
        Experimental design can be used to create an extended graph with replicate information. (experimental)
        
        
        :param proteinIDs: Vector of protein identifications
        :param peptideIDs: Vector of peptide identifications
        :param exp_des: Experimental Design
        :return: Writes its results into protein and (optionally also) peptide hits (as new score)
        
        .. rubric:: Overload:
        .. py:function:: inferPosteriorProbabilities(self, cmap: ConsensusMap , greedy_group_resolution: bool ) -> None
          :noindex:
        
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :return: Writes its protein ID results into the ConsensusMap
        
        .. rubric:: Overload:
        .. py:function:: inferPosteriorProbabilities(self, cmap: ConsensusMap , greedy_group_resolution: bool , exp_des: ExperimentalDesign ) -> None
          :noindex:
        
        Writes its results into protein and (optionally also) peptide hits (as new score)
        Optionally adds indistinguishable protein groups with separate scores, too
        Loops over all runs in the ConsensusMaps' protein IDs (experimental)
        
        
        :param cmap: ConsensusMaps with protein IDs.
        :param greedy_group_resolution: Adds indistinguishable protein groups with separate scores
        :param exp_des: Experimental Design
        :return: Writes its protein ID results into the ConsensusMap
    
        """
        if (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], pybool_t)):
            return self._inferPosteriorProbabilities_0(*args)
        elif (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], pybool_t)) and (isinstance(args[3], ExperimentalDesign)):
            return self._inferPosteriorProbabilities_1(*args)
        elif (len(args)==2) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], pybool_t)):
            return self._inferPosteriorProbabilities_2(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], ExperimentalDesign)):
            return self._inferPosteriorProbabilities_3(*args)
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

cdef class Deisotoper:
    """
    Cython implementation of _Deisotoper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Deisotoper.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Deisotoper rv = Deisotoper.__new__(Deisotoper)
       rv.inst = shared_ptr[_Deisotoper](new _Deisotoper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Deisotoper rv = Deisotoper.__new__(Deisotoper)
       rv.inst = shared_ptr[_Deisotoper](new _Deisotoper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Deisotoper](new _Deisotoper())
    
    def _init_1(self, Deisotoper in_0 ):
        """
        _init_1(self, in_0: Deisotoper ) -> None
        """
        assert isinstance(in_0, Deisotoper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Deisotoper](new _Deisotoper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Deisotoper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Deisotoper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    deisotopeAndSingleCharge = __static_Deisotoper_deisotopeAndSingleCharge
    deisotopeAndSingleChargeDefault = __static_Deisotoper_deisotopeAndSingleChargeDefault
    deisotopeWithAveragineModel = __static_Deisotoper_deisotopeWithAveragineModel 

cdef class ElutionModelFitter:
    """
    Cython implementation of _ElutionModelFitter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ElutionModelFitter.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ElutionModelFitter rv = ElutionModelFitter.__new__(ElutionModelFitter)
       rv.inst = shared_ptr[_ElutionModelFitter](new _ElutionModelFitter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ElutionModelFitter rv = ElutionModelFitter.__new__(ElutionModelFitter)
       rv.inst = shared_ptr[_ElutionModelFitter](new _ElutionModelFitter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Helper class for fitting elution models to features
        """
        self.inst = shared_ptr[_ElutionModelFitter](new _ElutionModelFitter())
    
    def _init_1(self, ElutionModelFitter in_0 ):
        """
        _init_1(self, in_0: ElutionModelFitter ) -> None
        """
        assert isinstance(in_0, ElutionModelFitter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ElutionModelFitter](new _ElutionModelFitter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Helper class for fitting elution models to features

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ElutionModelFitter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ElutionModelFitter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def fitElutionModels(self, FeatureMap features ):
        """
        fitElutionModels(self, features: FeatureMap ) -> None
        Fit models of elution profiles to all features (and validate them)
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
    
        self.inst.get().fitElutionModels((deref(features.inst.get())))
    
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

cdef class FIAMSDataProcessor:
    """
    Cython implementation of _FIAMSDataProcessor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FIAMSDataProcessor.html>`_
      -- Inherits from ['DefaultParamHandler']

      ADD PYTHON DOCUMENTATION HERE
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FIAMSDataProcessor rv = FIAMSDataProcessor.__new__(FIAMSDataProcessor)
       rv.inst = shared_ptr[_FIAMSDataProcessor](new _FIAMSDataProcessor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FIAMSDataProcessor rv = FIAMSDataProcessor.__new__(FIAMSDataProcessor)
       rv.inst = shared_ptr[_FIAMSDataProcessor](new _FIAMSDataProcessor(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Data processing for FIA-MS data
        """
        self.inst = shared_ptr[_FIAMSDataProcessor](new _FIAMSDataProcessor())
    
    def _init_1(self, FIAMSDataProcessor in_0 ):
        """
        _init_1(self, in_0: FIAMSDataProcessor ) -> None
        """
        assert isinstance(in_0, FIAMSDataProcessor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FIAMSDataProcessor](new _FIAMSDataProcessor((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Data processing for FIA-MS data

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FIAMSDataProcessor ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FIAMSDataProcessor)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment experiment , float n_seconds , MzTab output , bool load_cached_spectrum ):
        """
        run(self, experiment: MSExperiment , n_seconds: float , output: MzTab , load_cached_spectrum: bool ) -> bool
        Run the full analysis for the experiment for the given time interval\n
        
        The workflow steps are:
        - the time axis of the experiment is cut to the interval from 0 to n_seconds
        - the spectra are summed into one along the time axis with the bin size determined by mz and instrument resolution
        - data is smoothed by applying the Savitzky-Golay filter
        - peaks are picked
        - the accurate mass search for all the picked peaks is performed
        
        The intermediate summed spectra and picked peaks can be saved to the filesystem.
        Also, the results of the accurate mass search and the signal-to-noise information
        of the resulting spectrum is saved.
        
        
        :param experiment: Input MSExperiment
        :param n_seconds: Input number of seconds
        :param load_cached_spectrum: Load the cached picked spectrum if exists
        :param output: Output of the accurate mass search results
        :return: A boolean indicating if the picked spectrum was loaded from the cached file
        """
        assert isinstance(experiment, MSExperiment), 'arg experiment wrong type'
        assert isinstance(n_seconds, float), 'arg n_seconds wrong type'
        assert isinstance(output, MzTab), 'arg output wrong type'
        assert isinstance(load_cached_spectrum, pybool_t), 'arg load_cached_spectrum wrong type'
    
    
    
    
        cdef bool _r = self.inst.get().run((deref(experiment.inst.get())), (<float &>n_seconds), (deref(output.inst.get())), (<bool>load_cached_spectrum))
        py_result = <bool>_r
        return py_result
    
    def extractPeaks(self, MSSpectrum input_ ):
        """
        extractPeaks(self, input_: MSSpectrum ) -> MSSpectrum
        Pick peaks from the summed spectrum
        
        
        :param input: Input vector of spectra
        :return: A spectrum with picked peaks
        """
        assert isinstance(input_, MSSpectrum), 'arg input_ wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().extractPeaks((deref(input_.inst.get()))))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def convertToFeatureMap(self, MSSpectrum input_ ):
        """
        convertToFeatureMap(self, input_: MSSpectrum ) -> FeatureMap
        Convert a spectrum to a feature map with the corresponding polarity\n
        
        Applies `SavitzkyGolayFilter` and `PeakPickerHiRes`
        
        
        :param input: Input a picked spectrum
        :return: A feature map with the peaks converted to features and polarity from the parameters
        """
        assert isinstance(input_, MSSpectrum), 'arg input_ wrong type'
    
        cdef _FeatureMap * _r = new _FeatureMap(self.inst.get().convertToFeatureMap((deref(input_.inst.get()))))
        cdef FeatureMap py_result = FeatureMap.__new__(FeatureMap)
        py_result.inst = shared_ptr[_FeatureMap](_r)
        return py_result
    
    def trackNoise(self, MSSpectrum input_ ):
        """
        trackNoise(self, input_: MSSpectrum ) -> MSSpectrum
        Estimate noise for each peak\n
        
        Uses `SignalToNoiseEstimatorMedianRapid`
        
        
        :param input: Input a picked spectrum
        :return: A spectrum object storing logSN information
        """
        assert isinstance(input_, MSSpectrum), 'arg input_ wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().trackNoise((deref(input_.inst.get()))))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
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

cdef class FeatureDeconvolution:
    """
    Cython implementation of _FeatureDeconvolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureDeconvolution.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureDeconvolution rv = FeatureDeconvolution.__new__(FeatureDeconvolution)
       rv.inst = shared_ptr[_FeatureDeconvolution](new _FeatureDeconvolution(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureDeconvolution rv = FeatureDeconvolution.__new__(FeatureDeconvolution)
       rv.inst = shared_ptr[_FeatureDeconvolution](new _FeatureDeconvolution(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FeatureDeconvolution](new _FeatureDeconvolution())
    
    def _init_1(self, FeatureDeconvolution in_0 ):
        """
        _init_1(self, in_0: FeatureDeconvolution ) -> None
        """
        assert isinstance(in_0, FeatureDeconvolution), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureDeconvolution](new _FeatureDeconvolution((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureDeconvolution ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureDeconvolution)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, FeatureMap input , FeatureMap output , ConsensusMap cmap1 , ConsensusMap cmap2 ):
        """
        compute(self, input: FeatureMap , output: FeatureMap , cmap1: ConsensusMap , cmap2: ConsensusMap ) -> None
        """
        assert isinstance(input, FeatureMap), 'arg input wrong type'
        assert isinstance(output, FeatureMap), 'arg output wrong type'
        assert isinstance(cmap1, ConsensusMap), 'arg cmap1 wrong type'
        assert isinstance(cmap2, ConsensusMap), 'arg cmap2 wrong type'
    
    
    
    
        self.inst.get().compute((deref(input.inst.get())), (deref(output.inst.get())), (deref(cmap1.inst.get())), (deref(cmap2.inst.get())))
    
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
    CHARGEMODE_FD = __CHARGEMODE_FD 

cdef class ModifiedPeptideGenerator:
    """
    Cython implementation of _ModifiedPeptideGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModifiedPeptideGenerator.html>`_

    Generates modified peptides/proteins.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ModifiedPeptideGenerator rv = ModifiedPeptideGenerator.__new__(ModifiedPeptideGenerator)
       rv.inst = shared_ptr[_ModifiedPeptideGenerator](new _ModifiedPeptideGenerator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ModifiedPeptideGenerator rv = ModifiedPeptideGenerator.__new__(ModifiedPeptideGenerator)
       rv.inst = shared_ptr[_ModifiedPeptideGenerator](new _ModifiedPeptideGenerator(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ModifiedPeptideGenerator](new _ModifiedPeptideGenerator())
    
    def _init_1(self, ModifiedPeptideGenerator in_0 ):
        """
        _init_1(self, in_0: ModifiedPeptideGenerator ) -> None
        """
        assert isinstance(in_0, ModifiedPeptideGenerator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ModifiedPeptideGenerator](new _ModifiedPeptideGenerator((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ModifiedPeptideGenerator ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ModifiedPeptideGenerator)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    @staticmethod
    def getModifications(list modNames ):
        """
        getModifications(modNames: List[bytes] ) -> ModifiedPeptideGenerator_MapToResidueType
        """
        assert isinstance(modNames, list) and all(isinstance(li, bytes) for li in modNames), 'arg modNames wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in modNames:
           v0.push_back(_String(<char *>item0))
        cdef _ModifiedPeptideGenerator_MapToResidueType * _r = new _ModifiedPeptideGenerator_MapToResidueType(_ModifiedPeptideGenerator.getModifications(deref(v0)))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        modNames[:] = replace
        del v0
        cdef ModifiedPeptideGenerator_MapToResidueType py_result = ModifiedPeptideGenerator_MapToResidueType.__new__(ModifiedPeptideGenerator_MapToResidueType)
        py_result.inst = shared_ptr[_ModifiedPeptideGenerator_MapToResidueType](_r)
        return py_result
    
    @staticmethod
    def applyFixedModifications(ModifiedPeptideGenerator_MapToResidueType fixed_mods , AASequence peptide ):
        """
        applyFixedModifications(fixed_mods: ModifiedPeptideGenerator_MapToResidueType , peptide: AASequence ) -> None
        """
        assert isinstance(fixed_mods, ModifiedPeptideGenerator_MapToResidueType), 'arg fixed_mods wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
    
    
        _ModifiedPeptideGenerator.applyFixedModifications((deref(fixed_mods.inst.get())), (deref(peptide.inst.get())))
    
    @staticmethod
    def applyVariableModifications(ModifiedPeptideGenerator_MapToResidueType var_mods , AASequence peptide ,  max_variable_mods_per_peptide , list all_modified_peptides , bool keep_original ):
        """
        applyVariableModifications(var_mods: ModifiedPeptideGenerator_MapToResidueType , peptide: AASequence , max_variable_mods_per_peptide: int , all_modified_peptides: List[AASequence] , keep_original: bool ) -> None
        """
        assert isinstance(var_mods, ModifiedPeptideGenerator_MapToResidueType), 'arg var_mods wrong type'
        assert isinstance(peptide, AASequence), 'arg peptide wrong type'
        assert isinstance(max_variable_mods_per_peptide, int) and max_variable_mods_per_peptide >= 0, 'arg max_variable_mods_per_peptide wrong type'
        assert isinstance(all_modified_peptides, list) and all(isinstance(elemt_rec, AASequence) for elemt_rec in all_modified_peptides), 'arg all_modified_peptides wrong type'
        assert isinstance(keep_original, pybool_t), 'arg keep_original wrong type'
    
    
    
        cdef libcpp_vector[_AASequence] * v3 = new libcpp_vector[_AASequence]()
        cdef AASequence item3
        for item3 in all_modified_peptides:
            v3.push_back(deref(item3.inst.get()))
    
        _ModifiedPeptideGenerator.applyVariableModifications((deref(var_mods.inst.get())), (deref(peptide.inst.get())), (<size_t>max_variable_mods_per_peptide), deref(v3), (<bool>keep_original))
        cdef libcpp_vector[_AASequence].iterator it_all_modified_peptides = v3.begin()
        replace_0 = []
        while it_all_modified_peptides != v3.end():
            item3 = AASequence.__new__(AASequence)
            item3.inst = shared_ptr[_AASequence](new _AASequence(deref(it_all_modified_peptides)))
            replace_0.append(item3)
            inc(it_all_modified_peptides)
        all_modified_peptides[:] = replace_0
        del v3 

cdef class ModifiedPeptideGenerator_MapToResidueType:
    """
    Cython implementation of _ModifiedPeptideGenerator_MapToResidueType

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ModifiedPeptideGenerator_MapToResidueType.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ModifiedPeptideGenerator_MapToResidueType rv = ModifiedPeptideGenerator_MapToResidueType.__new__(ModifiedPeptideGenerator_MapToResidueType)
       rv.inst = shared_ptr[_ModifiedPeptideGenerator_MapToResidueType](new _ModifiedPeptideGenerator_MapToResidueType(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ModifiedPeptideGenerator_MapToResidueType rv = ModifiedPeptideGenerator_MapToResidueType.__new__(ModifiedPeptideGenerator_MapToResidueType)
       rv.inst = shared_ptr[_ModifiedPeptideGenerator_MapToResidueType](new _ModifiedPeptideGenerator_MapToResidueType(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ModifiedPeptideGenerator_MapToResidueType](new _ModifiedPeptideGenerator_MapToResidueType())
    
    def _init_1(self, ModifiedPeptideGenerator_MapToResidueType in_0 ):
        """
        _init_1(self, in_0: ModifiedPeptideGenerator_MapToResidueType ) -> None
        """
        assert isinstance(in_0, ModifiedPeptideGenerator_MapToResidueType), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ModifiedPeptideGenerator_MapToResidueType](new _ModifiedPeptideGenerator_MapToResidueType((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ModifiedPeptideGenerator_MapToResidueType ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ModifiedPeptideGenerator_MapToResidueType)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OpenSwathScoring:
    """
    Cython implementation of _OpenSwathScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathScoring.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenSwathScoring rv = OpenSwathScoring.__new__(OpenSwathScoring)
       rv.inst = shared_ptr[_OpenSwathScoring](new _OpenSwathScoring(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenSwathScoring rv = OpenSwathScoring.__new__(OpenSwathScoring)
       rv.inst = shared_ptr[_OpenSwathScoring](new _OpenSwathScoring(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenSwathScoring](new _OpenSwathScoring())
    
    def _init_1(self, OpenSwathScoring in_0 ):
        """
        _init_1(self, in_0: OpenSwathScoring ) -> None
        """
        assert isinstance(in_0, OpenSwathScoring), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenSwathScoring](new _OpenSwathScoring((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenSwathScoring ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenSwathScoring)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def initialize(self, double rt_normalization_factor ,  add_up_spectra , double spacing_for_spectra_resampling , double merge_spectra_by_peak_width_fraction , double drift_extra , OpenSwath_Scores_Usage su , bytes spectrum_addition_method , bytes spectrum_merge_method_type , bool use_ms1_ion_mobility , bool apply_im_peak_picking ):
        """
        initialize(self, rt_normalization_factor: float , add_up_spectra: int , spacing_for_spectra_resampling: float , merge_spectra_by_peak_width_fraction: float , drift_extra: float , su: OpenSwath_Scores_Usage , spectrum_addition_method: bytes , spectrum_merge_method_type: bytes , use_ms1_ion_mobility: bool , apply_im_peak_picking: bool ) -> None
        Initialize the scoring object\n
        Sets the parameters for the scoring
        
        
        :param rt_normalization_factor: Specifies the range of the normalized retention time space
        :param add_up_spectra: How many spectra to add up (default 1)
        :param spacing_for_spectra_resampling: Spacing factor for spectra addition
        :param merge_spectra_by_peak_width_fraction: Fraction of peak width to construct the number of spectra to add
        :param drift_extra: Extend the extraction window to gain a larger field of view beyond drift_upper - drift_lower (in percent)
        :param su: Which scores to actually compute
        :param spectrum_addition_method: Method to use for spectrum addition (valid: "simple", "resample")
        :param spectrum_merge_method_type: Type of method to use for spectrum addition. (valid: "fixed", "dynamic")
        :param use_ms1_ion_mobility: Use MS1 ion mobility extraction in DIA scores
        :param apply_im_peak_picking: Apply peak picking to the  extracted ion mobilograms
        """
        assert isinstance(rt_normalization_factor, float), 'arg rt_normalization_factor wrong type'
        assert isinstance(add_up_spectra, int), 'arg add_up_spectra wrong type'
        assert isinstance(spacing_for_spectra_resampling, float), 'arg spacing_for_spectra_resampling wrong type'
        assert isinstance(merge_spectra_by_peak_width_fraction, float), 'arg merge_spectra_by_peak_width_fraction wrong type'
        assert isinstance(drift_extra, float), 'arg drift_extra wrong type'
        assert isinstance(su, OpenSwath_Scores_Usage), 'arg su wrong type'
        assert isinstance(spectrum_addition_method, bytes), 'arg spectrum_addition_method wrong type'
        assert isinstance(spectrum_merge_method_type, bytes), 'arg spectrum_merge_method_type wrong type'
        assert isinstance(use_ms1_ion_mobility, pybool_t), 'arg use_ms1_ion_mobility wrong type'
        assert isinstance(apply_im_peak_picking, pybool_t), 'arg apply_im_peak_picking wrong type'
    
    
    
    
    
    
    
    
    
    
        self.inst.get().initialize((<double>rt_normalization_factor), (<int>add_up_spectra), (<double>spacing_for_spectra_resampling), (<double>merge_spectra_by_peak_width_fraction), (<double>drift_extra), (deref(su.inst.get())), (<libcpp_string>spectrum_addition_method), (<libcpp_string>spectrum_merge_method_type), (<bool>use_ms1_ion_mobility), (<bool>apply_im_peak_picking))
    
    def getNormalized_library_intensities_(self, list transitions , list normalized_library_intensity ):
        """
        getNormalized_library_intensities_(self, transitions: List[LightTransition] , normalized_library_intensity: List[float] ) -> None
        """
        assert isinstance(transitions, list) and all(isinstance(elemt_rec, LightTransition) for elemt_rec in transitions), 'arg transitions wrong type'
        assert isinstance(normalized_library_intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in normalized_library_intensity), 'arg normalized_library_intensity wrong type'
        cdef libcpp_vector[_LightTransition] * v0 = new libcpp_vector[_LightTransition]()
        cdef LightTransition item0
        for item0 in transitions:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[double] v1 = normalized_library_intensity
        self.inst.get().getNormalized_library_intensities_(deref(v0), v1)
        
        del v0 

cdef class OpenSwath_Scores:
    """
    Cython implementation of _OpenSwath_Scores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwath_Scores.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property elution_model_fit_score:
        def __set__(self, double elution_model_fit_score):
        
            self.inst.get().elution_model_fit_score = (<double>elution_model_fit_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().elution_model_fit_score
            py_result = <double>_r
            return py_result
    
    property library_corr:
        def __set__(self, double library_corr):
        
            self.inst.get().library_corr = (<double>library_corr)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_corr
            py_result = <double>_r
            return py_result
    
    property library_norm_manhattan:
        def __set__(self, double library_norm_manhattan):
        
            self.inst.get().library_norm_manhattan = (<double>library_norm_manhattan)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_norm_manhattan
            py_result = <double>_r
            return py_result
    
    property library_rootmeansquare:
        def __set__(self, double library_rootmeansquare):
        
            self.inst.get().library_rootmeansquare = (<double>library_rootmeansquare)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_rootmeansquare
            py_result = <double>_r
            return py_result
    
    property library_sangle:
        def __set__(self, double library_sangle):
        
            self.inst.get().library_sangle = (<double>library_sangle)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_sangle
            py_result = <double>_r
            return py_result
    
    property norm_rt_score:
        def __set__(self, double norm_rt_score):
        
            self.inst.get().norm_rt_score = (<double>norm_rt_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().norm_rt_score
            py_result = <double>_r
            return py_result
    
    property isotope_correlation:
        def __set__(self, double isotope_correlation):
        
            self.inst.get().isotope_correlation = (<double>isotope_correlation)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().isotope_correlation
            py_result = <double>_r
            return py_result
    
    property isotope_overlap:
        def __set__(self, double isotope_overlap):
        
            self.inst.get().isotope_overlap = (<double>isotope_overlap)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().isotope_overlap
            py_result = <double>_r
            return py_result
    
    property massdev_score:
        def __set__(self, double massdev_score):
        
            self.inst.get().massdev_score = (<double>massdev_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().massdev_score
            py_result = <double>_r
            return py_result
    
    property xcorr_coelution_score:
        def __set__(self, double xcorr_coelution_score):
        
            self.inst.get().xcorr_coelution_score = (<double>xcorr_coelution_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xcorr_coelution_score
            py_result = <double>_r
            return py_result
    
    property xcorr_shape_score:
        def __set__(self, double xcorr_shape_score):
        
            self.inst.get().xcorr_shape_score = (<double>xcorr_shape_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().xcorr_shape_score
            py_result = <double>_r
            return py_result
    
    property yseries_score:
        def __set__(self, double yseries_score):
        
            self.inst.get().yseries_score = (<double>yseries_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().yseries_score
            py_result = <double>_r
            return py_result
    
    property bseries_score:
        def __set__(self, double bseries_score):
        
            self.inst.get().bseries_score = (<double>bseries_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().bseries_score
            py_result = <double>_r
            return py_result
    
    property log_sn_score:
        def __set__(self, double log_sn_score):
        
            self.inst.get().log_sn_score = (<double>log_sn_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().log_sn_score
            py_result = <double>_r
            return py_result
    
    property weighted_coelution_score:
        def __set__(self, double weighted_coelution_score):
        
            self.inst.get().weighted_coelution_score = (<double>weighted_coelution_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().weighted_coelution_score
            py_result = <double>_r
            return py_result
    
    property weighted_xcorr_shape:
        def __set__(self, double weighted_xcorr_shape):
        
            self.inst.get().weighted_xcorr_shape = (<double>weighted_xcorr_shape)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().weighted_xcorr_shape
            py_result = <double>_r
            return py_result
    
    property weighted_massdev_score:
        def __set__(self, double weighted_massdev_score):
        
            self.inst.get().weighted_massdev_score = (<double>weighted_massdev_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().weighted_massdev_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_coelution_score:
        def __set__(self, double ms1_xcorr_coelution_score):
        
            self.inst.get().ms1_xcorr_coelution_score = (<double>ms1_xcorr_coelution_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_coelution_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_coelution_contrast_score:
        def __set__(self, double ms1_xcorr_coelution_contrast_score):
        
            self.inst.get().ms1_xcorr_coelution_contrast_score = (<double>ms1_xcorr_coelution_contrast_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_coelution_contrast_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_coelution_combined_score:
        def __set__(self, double ms1_xcorr_coelution_combined_score):
        
            self.inst.get().ms1_xcorr_coelution_combined_score = (<double>ms1_xcorr_coelution_combined_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_coelution_combined_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_shape_score:
        def __set__(self, double ms1_xcorr_shape_score):
        
            self.inst.get().ms1_xcorr_shape_score = (<double>ms1_xcorr_shape_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_shape_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_shape_contrast_score:
        def __set__(self, double ms1_xcorr_shape_contrast_score):
        
            self.inst.get().ms1_xcorr_shape_contrast_score = (<double>ms1_xcorr_shape_contrast_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_shape_contrast_score
            py_result = <double>_r
            return py_result
    
    property ms1_xcorr_shape_combined_score:
        def __set__(self, double ms1_xcorr_shape_combined_score):
        
            self.inst.get().ms1_xcorr_shape_combined_score = (<double>ms1_xcorr_shape_combined_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_xcorr_shape_combined_score
            py_result = <double>_r
            return py_result
    
    property ms1_ppm_score:
        def __set__(self, double ms1_ppm_score):
        
            self.inst.get().ms1_ppm_score = (<double>ms1_ppm_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_ppm_score
            py_result = <double>_r
            return py_result
    
    property ms1_isotope_correlation:
        def __set__(self, double ms1_isotope_correlation):
        
            self.inst.get().ms1_isotope_correlation = (<double>ms1_isotope_correlation)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_isotope_correlation
            py_result = <double>_r
            return py_result
    
    property ms1_isotope_overlap:
        def __set__(self, double ms1_isotope_overlap):
        
            self.inst.get().ms1_isotope_overlap = (<double>ms1_isotope_overlap)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_isotope_overlap
            py_result = <double>_r
            return py_result
    
    property ms1_mi_score:
        def __set__(self, double ms1_mi_score):
        
            self.inst.get().ms1_mi_score = (<double>ms1_mi_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_mi_score
            py_result = <double>_r
            return py_result
    
    property ms1_mi_contrast_score:
        def __set__(self, double ms1_mi_contrast_score):
        
            self.inst.get().ms1_mi_contrast_score = (<double>ms1_mi_contrast_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_mi_contrast_score
            py_result = <double>_r
            return py_result
    
    property ms1_mi_combined_score:
        def __set__(self, double ms1_mi_combined_score):
        
            self.inst.get().ms1_mi_combined_score = (<double>ms1_mi_combined_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().ms1_mi_combined_score
            py_result = <double>_r
            return py_result
    
    property library_manhattan:
        def __set__(self, double library_manhattan):
        
            self.inst.get().library_manhattan = (<double>library_manhattan)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_manhattan
            py_result = <double>_r
            return py_result
    
    property library_dotprod:
        def __set__(self, double library_dotprod):
        
            self.inst.get().library_dotprod = (<double>library_dotprod)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().library_dotprod
            py_result = <double>_r
            return py_result
    
    property intensity:
        def __set__(self, double intensity):
        
            self.inst.get().intensity = (<double>intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().intensity
            py_result = <double>_r
            return py_result
    
    property total_xic:
        def __set__(self, double total_xic):
        
            self.inst.get().total_xic = (<double>total_xic)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().total_xic
            py_result = <double>_r
            return py_result
    
    property nr_peaks:
        def __set__(self, double nr_peaks):
        
            self.inst.get().nr_peaks = (<double>nr_peaks)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().nr_peaks
            py_result = <double>_r
            return py_result
    
    property sn_ratio:
        def __set__(self, double sn_ratio):
        
            self.inst.get().sn_ratio = (<double>sn_ratio)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().sn_ratio
            py_result = <double>_r
            return py_result
    
    property mi_score:
        def __set__(self, double mi_score):
        
            self.inst.get().mi_score = (<double>mi_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mi_score
            py_result = <double>_r
            return py_result
    
    property weighted_mi_score:
        def __set__(self, double weighted_mi_score):
        
            self.inst.get().weighted_mi_score = (<double>weighted_mi_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().weighted_mi_score
            py_result = <double>_r
            return py_result
    
    property rt_difference:
        def __set__(self, double rt_difference):
        
            self.inst.get().rt_difference = (<double>rt_difference)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt_difference
            py_result = <double>_r
            return py_result
    
    property normalized_experimental_rt:
        def __set__(self, double normalized_experimental_rt):
        
            self.inst.get().normalized_experimental_rt = (<double>normalized_experimental_rt)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().normalized_experimental_rt
            py_result = <double>_r
            return py_result
    
    property raw_rt_score:
        def __set__(self, double raw_rt_score):
        
            self.inst.get().raw_rt_score = (<double>raw_rt_score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().raw_rt_score
            py_result = <double>_r
            return py_result
    
    property dotprod_score_dia:
        def __set__(self, double dotprod_score_dia):
        
            self.inst.get().dotprod_score_dia = (<double>dotprod_score_dia)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().dotprod_score_dia
            py_result = <double>_r
            return py_result
    
    property manhatt_score_dia:
        def __set__(self, double manhatt_score_dia):
        
            self.inst.get().manhatt_score_dia = (<double>manhatt_score_dia)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().manhatt_score_dia
            py_result = <double>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_OpenSwath_Scores](new _OpenSwath_Scores())
    
    def get_quick_lda_score(self, double library_corr_ , double library_norm_manhattan_ , double norm_rt_score_ , double xcorr_coelution_score_ , double xcorr_shape_score_ , double log_sn_score_ ):
        """
        get_quick_lda_score(self, library_corr_: float , library_norm_manhattan_: float , norm_rt_score_: float , xcorr_coelution_score_: float , xcorr_shape_score_: float , log_sn_score_: float ) -> float
        """
        assert isinstance(library_corr_, float), 'arg library_corr_ wrong type'
        assert isinstance(library_norm_manhattan_, float), 'arg library_norm_manhattan_ wrong type'
        assert isinstance(norm_rt_score_, float), 'arg norm_rt_score_ wrong type'
        assert isinstance(xcorr_coelution_score_, float), 'arg xcorr_coelution_score_ wrong type'
        assert isinstance(xcorr_shape_score_, float), 'arg xcorr_shape_score_ wrong type'
        assert isinstance(log_sn_score_, float), 'arg log_sn_score_ wrong type'
    
    
    
    
    
    
        cdef double _r = self.inst.get().get_quick_lda_score((<double>library_corr_), (<double>library_norm_manhattan_), (<double>norm_rt_score_), (<double>xcorr_coelution_score_), (<double>xcorr_shape_score_), (<double>log_sn_score_))
        py_result = <double>_r
        return py_result
    
    def calculate_lda_prescore(self, OpenSwath_Scores scores ):
        """
        calculate_lda_prescore(self, scores: OpenSwath_Scores ) -> float
        """
        assert isinstance(scores, OpenSwath_Scores), 'arg scores wrong type'
    
        cdef double _r = self.inst.get().calculate_lda_prescore((deref(scores.inst.get())))
        py_result = <double>_r
        return py_result
    
    def calculate_swath_lda_prescore(self, OpenSwath_Scores scores ):
        """
        calculate_swath_lda_prescore(self, scores: OpenSwath_Scores ) -> float
        """
        assert isinstance(scores, OpenSwath_Scores), 'arg scores wrong type'
    
        cdef double _r = self.inst.get().calculate_swath_lda_prescore((deref(scores.inst.get())))
        py_result = <double>_r
        return py_result 

cdef class OpenSwath_Scores_Usage:
    """
    Cython implementation of _OpenSwath_Scores_Usage

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwath_Scores_Usage.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property use_coelution_score_:
        def __set__(self, bool use_coelution_score_):
        
            self.inst.get().use_coelution_score_ = (<bool>use_coelution_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_coelution_score_
            py_result = <bool>_r
            return py_result
    
    property use_shape_score_:
        def __set__(self, bool use_shape_score_):
        
            self.inst.get().use_shape_score_ = (<bool>use_shape_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_shape_score_
            py_result = <bool>_r
            return py_result
    
    property use_rt_score_:
        def __set__(self, bool use_rt_score_):
        
            self.inst.get().use_rt_score_ = (<bool>use_rt_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_rt_score_
            py_result = <bool>_r
            return py_result
    
    property use_library_score_:
        def __set__(self, bool use_library_score_):
        
            self.inst.get().use_library_score_ = (<bool>use_library_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_library_score_
            py_result = <bool>_r
            return py_result
    
    property use_elution_model_score_:
        def __set__(self, bool use_elution_model_score_):
        
            self.inst.get().use_elution_model_score_ = (<bool>use_elution_model_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_elution_model_score_
            py_result = <bool>_r
            return py_result
    
    property use_intensity_score_:
        def __set__(self, bool use_intensity_score_):
        
            self.inst.get().use_intensity_score_ = (<bool>use_intensity_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_intensity_score_
            py_result = <bool>_r
            return py_result
    
    property use_total_xic_score_:
        def __set__(self, bool use_total_xic_score_):
        
            self.inst.get().use_total_xic_score_ = (<bool>use_total_xic_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_total_xic_score_
            py_result = <bool>_r
            return py_result
    
    property use_total_mi_score_:
        def __set__(self, bool use_total_mi_score_):
        
            self.inst.get().use_total_mi_score_ = (<bool>use_total_mi_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_total_mi_score_
            py_result = <bool>_r
            return py_result
    
    property use_nr_peaks_score_:
        def __set__(self, bool use_nr_peaks_score_):
        
            self.inst.get().use_nr_peaks_score_ = (<bool>use_nr_peaks_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_nr_peaks_score_
            py_result = <bool>_r
            return py_result
    
    property use_sn_score_:
        def __set__(self, bool use_sn_score_):
        
            self.inst.get().use_sn_score_ = (<bool>use_sn_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_sn_score_
            py_result = <bool>_r
            return py_result
    
    property use_mi_score_:
        def __set__(self, bool use_mi_score_):
        
            self.inst.get().use_mi_score_ = (<bool>use_mi_score_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_mi_score_
            py_result = <bool>_r
            return py_result
    
    property use_dia_scores_:
        def __set__(self, bool use_dia_scores_):
        
            self.inst.get().use_dia_scores_ = (<bool>use_dia_scores_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_dia_scores_
            py_result = <bool>_r
            return py_result
    
    property use_ms1_correlation:
        def __set__(self, bool use_ms1_correlation):
        
            self.inst.get().use_ms1_correlation = (<bool>use_ms1_correlation)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_ms1_correlation
            py_result = <bool>_r
            return py_result
    
    property use_ms1_fullscan:
        def __set__(self, bool use_ms1_fullscan):
        
            self.inst.get().use_ms1_fullscan = (<bool>use_ms1_fullscan)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_ms1_fullscan
            py_result = <bool>_r
            return py_result
    
    property use_ms1_mi:
        def __set__(self, bool use_ms1_mi):
        
            self.inst.get().use_ms1_mi = (<bool>use_ms1_mi)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_ms1_mi
            py_result = <bool>_r
            return py_result
    
    property use_uis_scores:
        def __set__(self, bool use_uis_scores):
        
            self.inst.get().use_uis_scores = (<bool>use_uis_scores)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().use_uis_scores
            py_result = <bool>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_OpenSwath_Scores_Usage](new _OpenSwath_Scores_Usage()) 

cdef class Peak2D:
    """
    Cython implementation of _Peak2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Peak2D.html>`_

    A 2-dimensional raw data point or peak.
    
    This data structure is intended for continuous data or peak data.
    If you want to annotated single peaks with meta data, use RichPeak2D instead
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Peak2D rv = Peak2D.__new__(Peak2D)
       rv.inst = shared_ptr[_Peak2D](new _Peak2D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Peak2D rv = Peak2D.__new__(Peak2D)
       rv.inst = shared_ptr[_Peak2D](new _Peak2D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Peak2D](new _Peak2D())
    
    def _init_1(self, Peak2D in_0 ):
        """
        _init_1(self, in_0: Peak2D ) -> None
        """
        assert isinstance(in_0, Peak2D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Peak2D](new _Peak2D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Peak2D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Peak2D)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Peak2D):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Peak2D other_casted = other
        cdef Peak2D self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class PeptideIdentification:
    """
    Cython implementation of _PeptideIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIdentification.html>`_
      -- Inherits from ['MetaInfoInterface']

    Represents the peptide hits for a spectrum
    
    This class is closely related to ProteinIdentification, which stores the protein hits
    and the general information about the identification run. More than one PeptideIdentification
    can belong to one ProteinIdentification. The general information about a
    PeptideIdentification has to be looked up in the corresponding ProteinIndentification, using
    the unique `identifier` that links the two.
    When loading PeptideHit instances from a File, the retention time and mass-to-charge ratio
    of the precursor spectrum can be accessed using getRT() and getMZ().
    This information can be used to map the peptide hits to an MSExperiment, a FeatureMap
    or a ConsensusMap using the IDMapper class
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideIdentification rv = PeptideIdentification.__new__(PeptideIdentification)
       rv.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideIdentification rv = PeptideIdentification.__new__(PeptideIdentification)
       rv.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification())
    
    def _init_1(self, PeptideIdentification in_0 ):
        """
        _init_1(self, in_0: PeptideIdentification ) -> None
        """
        assert isinstance(in_0, PeptideIdentification), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideIdentification](new _PeptideIdentification((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideIdentification ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentification)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getHits(self):
        """
        getHits(self) -> List[PeptideHit]
        Returns the peptide hits as const
        """
        _r = self.inst.get().getHits()
        py_result = []
        cdef libcpp_vector[_PeptideHit].iterator it__r = _r.begin()
        cdef PeptideHit item_py_result
        while it__r != _r.end():
           item_py_result = PeptideHit.__new__(PeptideHit)
           item_py_result.inst = shared_ptr[_PeptideHit](new _PeptideHit(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def insertHit(self, PeptideHit in_0 ):
        """
        insertHit(self, in_0: PeptideHit ) -> None
        Appends a peptide hit
        """
        assert isinstance(in_0, PeptideHit), 'arg in_0 wrong type'
    
        self.inst.get().insertHit((deref(in_0.inst.get())))
    
    def setHits(self, list in_0 ):
        """
        setHits(self, in_0: List[PeptideHit] ) -> None
        Sets the peptide hits
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, PeptideHit) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_PeptideHit] * v0 = new libcpp_vector[_PeptideHit]()
        cdef PeptideHit item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setHits(deref(v0))
        del v0
    
    def getSignificanceThreshold(self):
        """
        getSignificanceThreshold(self) -> float
        Returns the peptide significance threshold value
        """
        cdef double _r = self.inst.get().getSignificanceThreshold()
        py_result = <double>_r
        return py_result
    
    def setSignificanceThreshold(self, double value ):
        """
        setSignificanceThreshold(self, value: float ) -> None
        Setting of the peptide significance threshold value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setSignificanceThreshold((<double>value))
    
    def getScoreType(self):
        """
        getScoreType(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getScoreType()
        py_result = convOutputString(_r)
        return py_result
    
    def setScoreType(self,  in_0 ):
        """
        setScoreType(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setScoreType(deref((convString(in_0)).get()))
    
    def isHigherScoreBetter(self):
        """
        isHigherScoreBetter(self) -> bool
        """
        cdef bool _r = self.inst.get().isHigherScoreBetter()
        py_result = <bool>_r
        return py_result
    
    def setHigherScoreBetter(self, bool in_0 ):
        """
        setHigherScoreBetter(self, in_0: bool ) -> None
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setHigherScoreBetter((<bool>in_0))
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setIdentifier(self,  in_0 ):
        """
        setIdentifier(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setIdentifier(deref((convString(in_0)).get()))
    
    def hasMZ(self):
        """
        hasMZ(self) -> bool
        """
        cdef bool _r = self.inst.get().hasMZ()
        py_result = <bool>_r
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
    
    def hasRT(self):
        """
        hasRT(self) -> bool
        """
        cdef bool _r = self.inst.get().hasRT()
        py_result = <bool>_r
        return py_result
    
    def getRT(self):
        """
        getRT(self) -> float
        """
        cdef double _r = self.inst.get().getRT()
        py_result = <double>_r
        return py_result
    
    def setRT(self, double in_0 ):
        """
        setRT(self, in_0: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setRT((<double>in_0))
    
    def getBaseName(self):
        """
        getBaseName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getBaseName()
        py_result = convOutputString(_r)
        return py_result
    
    def setBaseName(self,  in_0 ):
        """
        setBaseName(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setBaseName(deref((convString(in_0)).get()))
    
    def getExperimentLabel(self):
        """
        getExperimentLabel(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getExperimentLabel()
        py_result = convOutputString(_r)
        return py_result
    
    def setExperimentLabel(self,  in_0 ):
        """
        setExperimentLabel(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setExperimentLabel(deref((convString(in_0)).get()))
    
    def sort(self):
        """
        sort(self) -> None
        """
        self.inst.get().sort()
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def getReferencingHits(self, list in_0 , set in_1 ):
        """
        getReferencingHits(self, in_0: List[PeptideHit] , in_1: Set[bytes] ) -> List[PeptideHit]
        Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, PeptideHit) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, set) and all(isinstance(i, bytes) for i in in_1), 'arg in_1 wrong type'
        cdef libcpp_vector[_PeptideHit] * v0 = new libcpp_vector[_PeptideHit]()
        cdef PeptideHit item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in in_1:
           v1.insert(_String(<char *>item1))
        _r = self.inst.get().getReferencingHits(deref(v0), deref(v1))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        in_1.clear()
        in_1.update(replace)
        del v1
        del v0
        py_result = []
        cdef libcpp_vector[_PeptideHit].iterator it__r = _r.begin()
        cdef PeptideHit item_py_result
        while it__r != _r.end():
           item_py_result = PeptideHit.__new__(PeptideHit)
           item_py_result.inst = shared_ptr[_PeptideHit](new _PeptideHit(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
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
        if not isinstance(other, PeptideIdentification):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef PeptideIdentification other_casted = other
        cdef PeptideIdentification self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class PercolatorFeatureSetHelper:
    """
    Cython implementation of _PercolatorFeatureSetHelper

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorFeatureSetHelper.html>`_

    Percolator feature set and integration helper
    
    This class contains functions to handle (compute, aggregate, integrate)
    Percolator features. This includes the calculation or extraction of
    Percolator features depending on the search engine(s) for later use with
    PercolatorAdapter. It also includes handling the reintegration of the
    percolator result into the set of Identifications
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PercolatorFeatureSetHelper rv = PercolatorFeatureSetHelper.__new__(PercolatorFeatureSetHelper)
       rv.inst = shared_ptr[_PercolatorFeatureSetHelper](new _PercolatorFeatureSetHelper(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PercolatorFeatureSetHelper rv = PercolatorFeatureSetHelper.__new__(PercolatorFeatureSetHelper)
       rv.inst = shared_ptr[_PercolatorFeatureSetHelper](new _PercolatorFeatureSetHelper(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PercolatorFeatureSetHelper](new _PercolatorFeatureSetHelper())
    
    def _init_1(self, PercolatorFeatureSetHelper in_0 ):
        """
        _init_1(self, in_0: PercolatorFeatureSetHelper ) -> None
        """
        assert isinstance(in_0, PercolatorFeatureSetHelper), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PercolatorFeatureSetHelper](new _PercolatorFeatureSetHelper((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PercolatorFeatureSetHelper ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PercolatorFeatureSetHelper)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def concatMULTISEPeptideIds(self, PeptideIdentificationList all_peptide_ids , PeptideIdentificationList new_peptide_ids ,  search_engine ):
        """
        concatMULTISEPeptideIds(self, all_peptide_ids: PeptideIdentificationList , new_peptide_ids: PeptideIdentificationList , search_engine: Union[bytes, str, String] ) -> None
        Appends a vector of PeptideIdentification to another and prepares Percolator features in MetaInfo (With the respective key "CONCAT:" + search_engine)
        
        
        :param all_peptide_ids: PeptideIdentification vector to append to
        :param new_peptide_ids: PeptideIdentification vector to be appended
        :param search_engine: Search engine to depend on for feature creation
        """
        assert isinstance(all_peptide_ids, PeptideIdentificationList), 'arg all_peptide_ids wrong type'
        assert isinstance(new_peptide_ids, PeptideIdentificationList), 'arg new_peptide_ids wrong type'
        assert (isinstance(search_engine, str) or isinstance(search_engine, bytes) or isinstance(search_engine, String)), 'arg search_engine wrong type'
    
    
    
        self.inst.get().concatMULTISEPeptideIds((deref(all_peptide_ids.inst.get())), (deref(new_peptide_ids.inst.get())), deref((convString(search_engine)).get()))
    
    def mergeMULTISEPeptideIds(self, PeptideIdentificationList all_peptide_ids , PeptideIdentificationList new_peptide_ids ,  search_engine ):
        """
        mergeMULTISEPeptideIds(self, all_peptide_ids: PeptideIdentificationList , new_peptide_ids: PeptideIdentificationList , search_engine: Union[bytes, str, String] ) -> None
        Merges a vector of PeptideIdentification into another and prepares the merged MetaInfo and scores for collection in addMULTISEFeatures for feature registration
        
        
        :param all_peptide_idsL: PeptideIdentification vector to be merged into
        :param new_peptide_idsL: PeptideIdentification vector to merge
        :param search_engineL: Search engine to create features from their scores
        """
        assert isinstance(all_peptide_ids, PeptideIdentificationList), 'arg all_peptide_ids wrong type'
        assert isinstance(new_peptide_ids, PeptideIdentificationList), 'arg new_peptide_ids wrong type'
        assert (isinstance(search_engine, str) or isinstance(search_engine, bytes) or isinstance(search_engine, String)), 'arg search_engine wrong type'
    
    
    
        self.inst.get().mergeMULTISEPeptideIds((deref(all_peptide_ids.inst.get())), (deref(new_peptide_ids.inst.get())), deref((convString(search_engine)).get()))
    
    def mergeMULTISEProteinIds(self, list all_protein_ids , list new_protein_ids ):
        """
        mergeMULTISEProteinIds(self, all_protein_ids: List[ProteinIdentification] , new_protein_ids: List[ProteinIdentification] ) -> None
        Concatenates SearchParameter of multiple search engine runs and merges PeptideEvidences, collects used search engines in MetaInfo for collection in addMULTISEFeatures for feature registration
        
        
        :param all_protein_ids: ProteinIdentification vector to be merged into
        :param new_protein_ids: ProteinIdentification vector to merge
        """
        assert isinstance(all_protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in all_protein_ids), 'arg all_protein_ids wrong type'
        assert isinstance(new_protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in new_protein_ids), 'arg new_protein_ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in all_protein_ids:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in new_protein_ids:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().mergeMULTISEProteinIds(deref(v0), deref(v1))
        cdef libcpp_vector[_ProteinIdentification].iterator it_new_protein_ids = v1.begin()
        replace_0 = []
        while it_new_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_new_protein_ids)))
            replace_0.append(item1)
            inc(it_new_protein_ids)
        new_protein_ids[:] = replace_0
        del v1
        cdef libcpp_vector[_ProteinIdentification].iterator it_all_protein_ids = v0.begin()
        replace_0 = []
        while it_all_protein_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_all_protein_ids)))
            replace_0.append(item0)
            inc(it_all_protein_ids)
        all_protein_ids[:] = replace_0
        del v0
    
    def addMSGFFeatures(self, PeptideIdentificationList peptide_ids , list feature_set ):
        """
        addMSGFFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None
        Creates and adds MSGF+ specific Percolator features and registers them in feature_set. MSGF+ should be run with the addFeatures flag enabled
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in feature_set:
           v1.push_back(_String(<char *>item1))
        self.inst.get().addMSGFFeatures((deref(peptide_ids.inst.get())), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        feature_set[:] = replace
        del v1
    
    def addXTANDEMFeatures(self, PeptideIdentificationList peptide_ids , list feature_set ):
        """
        addXTANDEMFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None
        Creates and adds X!Tandem specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in feature_set:
           v1.push_back(_String(<char *>item1))
        self.inst.get().addXTANDEMFeatures((deref(peptide_ids.inst.get())), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        feature_set[:] = replace
        del v1
    
    def addCOMETFeatures(self, PeptideIdentificationList peptide_ids , list feature_set ):
        """
        addCOMETFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None
        Creates and adds Comet specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in feature_set:
           v1.push_back(_String(<char *>item1))
        self.inst.get().addCOMETFeatures((deref(peptide_ids.inst.get())), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        feature_set[:] = replace
        del v1
    
    def addMASCOTFeatures(self, PeptideIdentificationList peptide_ids , list feature_set ):
        """
        addMASCOTFeatures(self, peptide_ids: PeptideIdentificationList , feature_set: List[bytes] ) -> None
        Creates and adds Mascot specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param feature_set: Register of added features
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in feature_set:
           v1.push_back(_String(<char *>item1))
        self.inst.get().addMASCOTFeatures((deref(peptide_ids.inst.get())), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        feature_set[:] = replace
        del v1
    
    def addMULTISEFeatures(self, PeptideIdentificationList peptide_ids , list search_engines_used , list feature_set , bool complete_only , bool limits_imputation ):
        """
        addMULTISEFeatures(self, peptide_ids: PeptideIdentificationList , search_engines_used: List[bytes] , feature_set: List[bytes] , complete_only: bool , limits_imputation: bool ) -> None
        Adds multiple search engine specific Percolator features and registers them in feature_set
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param search_engines_used: The list of search engines to be considered
        :param feature_set: Register of added features
        :param complete_only: Will only add features for PeptideIdentifications where all given search engines identified something
        :param limits_imputation: Uses C++ numeric limits as imputed values instead of min/max of that feature
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(search_engines_used, list) and all(isinstance(li, bytes) for li in search_engines_used), 'arg search_engines_used wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
        assert isinstance(complete_only, pybool_t), 'arg complete_only wrong type'
        assert isinstance(limits_imputation, pybool_t), 'arg limits_imputation wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in search_engines_used:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in feature_set:
           v2.push_back(_String(<char *>item2))
    
    
        self.inst.get().addMULTISEFeatures((deref(peptide_ids.inst.get())), deref(v1), deref(v2), (<bool>complete_only), (<bool>limits_imputation))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        feature_set[:] = replace
        del v2
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        search_engines_used[:] = replace
        del v1
    
    def addCONCATSEFeatures(self, PeptideIdentificationList peptide_id_list , list search_engines_used , list feature_set ):
        """
        addCONCATSEFeatures(self, peptide_id_list: PeptideIdentificationList , search_engines_used: List[bytes] , feature_set: List[bytes] ) -> None
        Adds multiple search engine specific Percolator features and registers them in feature_set
        
        This struct can be used to store both peak or feature indices
        
        
        :param peptide_ids: PeptideIdentification vector to create Percolator features in
        :param search_engines_used: The list of search engines to be considered
        :param feature_set: Register of added features
        """
        assert isinstance(peptide_id_list, PeptideIdentificationList), 'arg peptide_id_list wrong type'
        assert isinstance(search_engines_used, list) and all(isinstance(li, bytes) for li in search_engines_used), 'arg search_engines_used wrong type'
        assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in search_engines_used:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in feature_set:
           v2.push_back(_String(<char *>item2))
        self.inst.get().addCONCATSEFeatures((deref(peptide_id_list.inst.get())), deref(v1), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        feature_set[:] = replace
        del v2
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        search_engines_used[:] = replace
        del v1
    
    def checkExtraFeatures(self, list psms , list extra_features ):
        """
        checkExtraFeatures(self, psms: List[PeptideHit] , extra_features: List[bytes] ) -> None
        Checks and removes requested extra Percolator features that are actually unavailable (to compute)
        
        
        :param psms: The vector of PeptideHit to be checked
        :param extra_features: The list of requested extra features
        """
        assert isinstance(psms, list) and all(isinstance(elemt_rec, PeptideHit) for elemt_rec in psms), 'arg psms wrong type'
        assert isinstance(extra_features, list) and all(isinstance(li, bytes) for li in extra_features), 'arg extra_features wrong type'
        cdef libcpp_vector[_PeptideHit] * v0 = new libcpp_vector[_PeptideHit]()
        cdef PeptideHit item0
        for item0 in psms:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in extra_features:
           v1.push_back(_String(<char *>item1))
        self.inst.get().checkExtraFeatures(deref(v0), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        extra_features[:] = replace
        del v1
        cdef libcpp_vector[_PeptideHit].iterator it_psms = v0.begin()
        replace_0 = []
        while it_psms != v0.end():
            item0 = PeptideHit.__new__(PeptideHit)
            item0.inst = shared_ptr[_PeptideHit](new _PeptideHit(deref(it_psms)))
            replace_0.append(item0)
            inc(it_psms)
        psms[:] = replace_0
        del v0 

cdef class ResidueModification:
    """
    Cython implementation of _ResidueModification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ResidueModification.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).getFullId().c_str() )

    
    def __copy__(self):
       cdef ResidueModification rv = ResidueModification.__new__(ResidueModification)
       rv.inst = shared_ptr[_ResidueModification](new _ResidueModification(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ResidueModification rv = ResidueModification.__new__(ResidueModification)
       rv.inst = shared_ptr[_ResidueModification](new _ResidueModification(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ResidueModification](new _ResidueModification())
    
    def _init_1(self, ResidueModification in_0 ):
        """
        _init_1(self, in_0: ResidueModification ) -> None
        """
        assert isinstance(in_0, ResidueModification), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ResidueModification](new _ResidueModification((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ResidueModification ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ResidueModification)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setId(self,  id_ ):
        """
        setId(self, id_: Union[bytes, str, String] ) -> None
        Sets the identifier of the modification
        """
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
    
        self.inst.get().setId(deref((convString(id_)).get()))
    
    def getId(self):
        """
        getId(self) -> Union[bytes, str, String]
        Returns the identifier of the modification
        """
        cdef _String _r = self.inst.get().getId()
        py_result = convOutputString(_r)
        return py_result
    
    def setFullId(self,  full_id ):
        """
        setFullId(self, full_id: Union[bytes, str, String] ) -> None
        Sets the full identifier (Unimod Accession + origin, if available)
        """
        assert (isinstance(full_id, str) or isinstance(full_id, bytes) or isinstance(full_id, String)), 'arg full_id wrong type'
    
        self.inst.get().setFullId(deref((convString(full_id)).get()))
    
    def getFullId(self):
        """
        getFullId(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getFullId()
        py_result = convOutputString(_r)
        return py_result
    
    def getUniModRecordId(self):
        """
        getUniModRecordId(self) -> int
        Gets the unimod record id
        """
        cdef int _r = self.inst.get().getUniModRecordId()
        py_result = <int>_r
        return py_result
    
    def setUniModRecordId(self,  id_ ):
        """
        setUniModRecordId(self, id_: int ) -> None
        Sets the unimod record id
        """
        assert isinstance(id_, int), 'arg id_ wrong type'
    
        self.inst.get().setUniModRecordId((<int>id_))
    
    def getUniModAccession(self):
        """
        getUniModAccession(self) -> Union[bytes, str, String]
        Returns the unimod accession if available
        """
        cdef _String _r = self.inst.get().getUniModAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def setPSIMODAccession(self,  id_ ):
        """
        setPSIMODAccession(self, id_: Union[bytes, str, String] ) -> None
        Sets the MOD-XXXXX accession of PSI-MOD
        """
        assert (isinstance(id_, str) or isinstance(id_, bytes) or isinstance(id_, String)), 'arg id_ wrong type'
    
        self.inst.get().setPSIMODAccession(deref((convString(id_)).get()))
    
    def getPSIMODAccession(self):
        """
        getPSIMODAccession(self) -> Union[bytes, str, String]
        Returns the PSI-MOD accession if available
        """
        cdef _String _r = self.inst.get().getPSIMODAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def setFullName(self,  full_name ):
        """
        setFullName(self, full_name: Union[bytes, str, String] ) -> None
        Sets the full name of the modification; must NOT contain the origin (or . for terminals!)
        """
        assert (isinstance(full_name, str) or isinstance(full_name, bytes) or isinstance(full_name, String)), 'arg full_name wrong type'
    
        self.inst.get().setFullName(deref((convString(full_name)).get()))
    
    def getFullName(self):
        """
        getFullName(self) -> Union[bytes, str, String]
        Returns the full name of the modification
        """
        cdef _String _r = self.inst.get().getFullName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of modification
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the PSI-MS-label if available; e.g. Mascot uses this name
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def _setTermSpecificity_0(self, int term_spec ):
        """
        _setTermSpecificity_0(self, term_spec: int ) -> None
        Sets the term specificity
        """
        assert term_spec in [0, 1, 2, 3, 4, 5], 'arg term_spec wrong type'
    
        self.inst.get().setTermSpecificity((<_TermSpecificity>term_spec))
    
    def _setTermSpecificity_1(self,  name ):
        """
        _setTermSpecificity_1(self, name: Union[bytes, str, String] ) -> None
        Sets the terminal specificity using a name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setTermSpecificity(deref((convString(name)).get()))
    
    def setTermSpecificity(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setTermSpecificity(self, term_spec: int ) -> None
          :noindex:
        
        Sets the term specificity

        
        .. rubric:: Overload:
        .. py:function:: setTermSpecificity(self, name: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the terminal specificity using a name
    
        """
        if (len(args)==1) and (args[0] in [0, 1, 2, 3, 4, 5]):
            return self._setTermSpecificity_0(*args)
        elif (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setTermSpecificity_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getTermSpecificity(self):
        """
        getTermSpecificity(self) -> int
        Returns terminal specificity
        """
        cdef _TermSpecificity _r = self.inst.get().getTermSpecificity()
        py_result = <int>_r
        return py_result
    
    def getTermSpecificityName(self, int in_0 ):
        """
        getTermSpecificityName(self, in_0: int ) -> Union[bytes, str, String]
        Returns the name of the terminal specificity
        """
        assert in_0 in [0, 1, 2, 3, 4, 5], 'arg in_0 wrong type'
    
        cdef _String _r = self.inst.get().getTermSpecificityName((<_TermSpecificity>in_0))
        py_result = convOutputString(_r)
        return py_result
    
    def setOrigin(self, bytes origin ):
        """
        setOrigin(self, origin: bytes ) -> None
        Sets the origin (i.e. modified amino acid)
        """
        assert isinstance(origin, bytes) and len(origin) == 1, 'arg origin wrong type'
    
        self.inst.get().setOrigin((<char>((origin)[0])))
    
    def getOrigin(self):
        """
        getOrigin(self) -> bytes
        Returns the origin (i.e. modified amino acid)
        """
        cdef char  _r = self.inst.get().getOrigin()
        py_result = chr(<char>(_r))
        return py_result
    
    def _setSourceClassification_0(self,  classification ):
        """
        _setSourceClassification_0(self, classification: Union[bytes, str, String] ) -> None
        Classification as defined by the PSI-MOD
        """
        assert (isinstance(classification, str) or isinstance(classification, bytes) or isinstance(classification, String)), 'arg classification wrong type'
    
        self.inst.get().setSourceClassification(deref((convString(classification)).get()))
    
    def _setSourceClassification_1(self, int classification ):
        """
        _setSourceClassification_1(self, classification: int ) -> None
        Sets the source classification
        """
        assert classification in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 'arg classification wrong type'
    
        self.inst.get().setSourceClassification((<_SourceClassification>classification))
    
    def setSourceClassification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setSourceClassification(self, classification: Union[bytes, str, String] ) -> None
          :noindex:
        
        Classification as defined by the PSI-MOD

        
        .. rubric:: Overload:
        .. py:function:: setSourceClassification(self, classification: int ) -> None
          :noindex:
        
        Sets the source classification
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setSourceClassification_0(*args)
        elif (len(args)==1) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]):
            return self._setSourceClassification_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSourceClassification(self):
        """
        getSourceClassification(self) -> int
        Returns the source classification, if none was set, it is unspecific
        """
        cdef _SourceClassification _r = self.inst.get().getSourceClassification()
        py_result = <int>_r
        return py_result
    
    def getSourceClassificationName(self, int classification ):
        """
        getSourceClassificationName(self, classification: int ) -> Union[bytes, str, String]
        Returns the classification
        """
        assert classification in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], 'arg classification wrong type'
    
        cdef _String _r = self.inst.get().getSourceClassificationName((<_SourceClassification>classification))
        py_result = convOutputString(_r)
        return py_result
    
    def setAverageMass(self, double mass ):
        """
        setAverageMass(self, mass: float ) -> None
        Sets the average mass
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst.get().setAverageMass((<double>mass))
    
    def getAverageMass(self):
        """
        getAverageMass(self) -> float
        Returns the average mass if set
        """
        cdef double _r = self.inst.get().getAverageMass()
        py_result = <double>_r
        return py_result
    
    def setMonoMass(self, double mass ):
        """
        setMonoMass(self, mass: float ) -> None
        Sets the monoisotopic mass (this must include the weight of the residue itself!)
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst.get().setMonoMass((<double>mass))
    
    def getMonoMass(self):
        """
        getMonoMass(self) -> float
        Return the monoisotopic mass, or 0.0 if not set
        """
        cdef double _r = self.inst.get().getMonoMass()
        py_result = <double>_r
        return py_result
    
    def setDiffAverageMass(self, double mass ):
        """
        setDiffAverageMass(self, mass: float ) -> None
        Sets the difference average mass
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst.get().setDiffAverageMass((<double>mass))
    
    def getDiffAverageMass(self):
        """
        getDiffAverageMass(self) -> float
        Returns the difference average mass, or 0.0 if not set
        """
        cdef double _r = self.inst.get().getDiffAverageMass()
        py_result = <double>_r
        return py_result
    
    def setDiffMonoMass(self, double mass ):
        """
        setDiffMonoMass(self, mass: float ) -> None
        Sets the difference monoisotopic mass
        """
        assert isinstance(mass, float), 'arg mass wrong type'
    
        self.inst.get().setDiffMonoMass((<double>mass))
    
    def getDiffMonoMass(self):
        """
        getDiffMonoMass(self) -> float
        Returns the diff monoisotopic mass, or 0.0 if not set
        """
        cdef double _r = self.inst.get().getDiffMonoMass()
        py_result = <double>_r
        return py_result
    
    def setFormula(self,  composition ):
        """
        setFormula(self, composition: Union[bytes, str, String] ) -> None
        Sets the formula (no masses will be changed)
        """
        assert (isinstance(composition, str) or isinstance(composition, bytes) or isinstance(composition, String)), 'arg composition wrong type'
    
        self.inst.get().setFormula(deref((convString(composition)).get()))
    
    def getFormula(self):
        """
        getFormula(self) -> Union[bytes, str, String]
        Returns the chemical formula if set
        """
        cdef _String _r = self.inst.get().getFormula()
        py_result = convOutputString(_r)
        return py_result
    
    def setDiffFormula(self, EmpiricalFormula diff_formula ):
        """
        setDiffFormula(self, diff_formula: EmpiricalFormula ) -> None
        Sets diff formula (no masses will be changed)
        """
        assert isinstance(diff_formula, EmpiricalFormula), 'arg diff_formula wrong type'
    
        self.inst.get().setDiffFormula((deref(diff_formula.inst.get())))
    
    def getDiffFormula(self):
        """
        getDiffFormula(self) -> EmpiricalFormula
        Returns the diff formula if one was set
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getDiffFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def setSynonyms(self, set synonyms ):
        """
        setSynonyms(self, synonyms: Set[bytes] ) -> None
        Sets the synonyms of that modification
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
        Adds a synonym to the unique list
        """
        assert (isinstance(synonym, str) or isinstance(synonym, bytes) or isinstance(synonym, String)), 'arg synonym wrong type'
    
        self.inst.get().addSynonym(deref((convString(synonym)).get()))
    
    def getSynonyms(self):
        """
        getSynonyms(self) -> Set[bytes]
        Returns the set of synonyms
        """
        _r = self.inst.get().getSynonyms()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setNeutralLossDiffFormulas(self, list diff_formulas ):
        """
        setNeutralLossDiffFormulas(self, diff_formulas: List[EmpiricalFormula] ) -> None
        Sets the neutral loss formula
        """
        assert isinstance(diff_formulas, list) and all(isinstance(elemt_rec, EmpiricalFormula) for elemt_rec in diff_formulas), 'arg diff_formulas wrong type'
        cdef libcpp_vector[_EmpiricalFormula] * v0 = new libcpp_vector[_EmpiricalFormula]()
        cdef EmpiricalFormula item0
        for item0 in diff_formulas:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setNeutralLossDiffFormulas(deref(v0))
        cdef libcpp_vector[_EmpiricalFormula].iterator it_diff_formulas = v0.begin()
        replace_0 = []
        while it_diff_formulas != v0.end():
            item0 = EmpiricalFormula.__new__(EmpiricalFormula)
            item0.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(it_diff_formulas)))
            replace_0.append(item0)
            inc(it_diff_formulas)
        diff_formulas[:] = replace_0
        del v0
    
    def getNeutralLossDiffFormulas(self):
        """
        getNeutralLossDiffFormulas(self) -> List[EmpiricalFormula]
        Returns the neutral loss diff formula (if available)
        """
        _r = self.inst.get().getNeutralLossDiffFormulas()
        py_result = []
        cdef libcpp_vector[_EmpiricalFormula].iterator it__r = _r.begin()
        cdef EmpiricalFormula item_py_result
        while it__r != _r.end():
           item_py_result = EmpiricalFormula.__new__(EmpiricalFormula)
           item_py_result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setNeutralLossMonoMasses(self, list mono_masses ):
        """
        setNeutralLossMonoMasses(self, mono_masses: List[float] ) -> None
        Sets the neutral loss mono weight
        """
        assert isinstance(mono_masses, list) and all(isinstance(elemt_rec, float) for elemt_rec in mono_masses), 'arg mono_masses wrong type'
        cdef libcpp_vector[double] v0 = mono_masses
        self.inst.get().setNeutralLossMonoMasses(v0)
        
    
    def getNeutralLossMonoMasses(self):
        """
        getNeutralLossMonoMasses(self) -> List[float]
        Returns the neutral loss mono weight
        """
        _r = self.inst.get().getNeutralLossMonoMasses()
        cdef list py_result = _r
        return py_result
    
    def setNeutralLossAverageMasses(self, list average_masses ):
        """
        setNeutralLossAverageMasses(self, average_masses: List[float] ) -> None
        Sets the neutral loss average weight
        """
        assert isinstance(average_masses, list) and all(isinstance(elemt_rec, float) for elemt_rec in average_masses), 'arg average_masses wrong type'
        cdef libcpp_vector[double] v0 = average_masses
        self.inst.get().setNeutralLossAverageMasses(v0)
        
    
    def getNeutralLossAverageMasses(self):
        """
        getNeutralLossAverageMasses(self) -> List[float]
        Returns the neutral loss average weight
        """
        _r = self.inst.get().getNeutralLossAverageMasses()
        cdef list py_result = _r
        return py_result
    
    def hasNeutralLoss(self):
        """
        hasNeutralLoss(self) -> bool
        Returns true if a neutral loss formula is set
        """
        cdef bool _r = self.inst.get().hasNeutralLoss()
        py_result = <bool>_r
        return py_result
    
    def isUserDefined(self):
        """
        isUserDefined(self) -> bool
        Returns true if it is a user-defined modification (empty id)
        """
        cdef bool _r = self.inst.get().isUserDefined()
        py_result = <bool>_r
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, ResidueModification):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ResidueModification other_casted = other
        cdef ResidueModification self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    SourceClassification = __SourceClassification
    TermSpecificity = __TermSpecificity 

cdef class SimpleSearchEngineAlgorithm:
    """
    Cython implementation of _SimpleSearchEngineAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SimpleSearchEngineAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SimpleSearchEngineAlgorithm rv = SimpleSearchEngineAlgorithm.__new__(SimpleSearchEngineAlgorithm)
       rv.inst = shared_ptr[_SimpleSearchEngineAlgorithm](new _SimpleSearchEngineAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SimpleSearchEngineAlgorithm rv = SimpleSearchEngineAlgorithm.__new__(SimpleSearchEngineAlgorithm)
       rv.inst = shared_ptr[_SimpleSearchEngineAlgorithm](new _SimpleSearchEngineAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SimpleSearchEngineAlgorithm](new _SimpleSearchEngineAlgorithm())
    
    def _init_1(self, SimpleSearchEngineAlgorithm in_0 ):
        """
        _init_1(self, in_0: SimpleSearchEngineAlgorithm ) -> None
        """
        assert isinstance(in_0, SimpleSearchEngineAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SimpleSearchEngineAlgorithm](new _SimpleSearchEngineAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SimpleSearchEngineAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SimpleSearchEngineAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def search(self,  in_mzML ,  in_db , list prot_ids , PeptideIdentificationList pep_ids ):
        """
        search(self, in_mzML: Union[bytes, str, String] , in_db: Union[bytes, str, String] , prot_ids: List[ProteinIdentification] , pep_ids: PeptideIdentificationList ) -> None
        """
        assert (isinstance(in_mzML, str) or isinstance(in_mzML, bytes) or isinstance(in_mzML, String)), 'arg in_mzML wrong type'
        assert (isinstance(in_db, str) or isinstance(in_db, bytes) or isinstance(in_db, String)), 'arg in_db wrong type'
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in prot_ids:
            v2.push_back(deref(item2.inst.get()))
    
        self.inst.get().search(deref((convString(in_mzML)).get()), deref((convString(in_db)).get()), deref(v2), (deref(pep_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v2.begin()
        replace_0 = []
        while it_prot_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item2)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v2
    
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

cdef class SpectrumAccessOpenMS:
    """
    Cython implementation of _SpectrumAccessOpenMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessOpenMS.html>`_
      -- Inherits from ['ISpectrumAccess']

    An implementation of the OpenSWATH Spectrum Access interface using OpenMS
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAccessOpenMS rv = SpectrumAccessOpenMS.__new__(SpectrumAccessOpenMS)
       rv.inst = shared_ptr[_SpectrumAccessOpenMS](new _SpectrumAccessOpenMS(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAccessOpenMS rv = SpectrumAccessOpenMS.__new__(SpectrumAccessOpenMS)
       rv.inst = shared_ptr[_SpectrumAccessOpenMS](new _SpectrumAccessOpenMS(deref(self.inst.get())))
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
    
        self.inst = shared_ptr[_SpectrumAccessOpenMS](new _SpectrumAccessOpenMS((deref(in_0.inst.get()))))
    
    def _init_2(self, MSExperiment ms_experiment ):
        """
        _init_2(self, ms_experiment: MSExperiment ) -> None
        """
        assert isinstance(ms_experiment, MSExperiment), 'arg ms_experiment wrong type'
        cdef shared_ptr[_MSExperiment] input_ms_experiment = ms_experiment.inst
        self.inst = shared_ptr[_SpectrumAccessOpenMS](new _SpectrumAccessOpenMS(input_ms_experiment))
        ms_experiment.inst = input_ms_experiment
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMS ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, ms_experiment: MSExperiment ) -> None
          :noindex:
    
        """
        if kwargs.get("__createUnsafeObject__") is True:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessOpenMS)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], MSExperiment)):
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

cdef class TM_DataPoint:
    """
    Cython implementation of _TM_DataPoint

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TM_DataPoint.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property first:
        def __set__(self, double first):
        
            self.inst.get().first = (<double>first)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().first
            py_result = <double>_r
            return py_result
    
    property second:
        def __set__(self, double second):
        
            self.inst.get().second = (<double>second)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().second
            py_result = <double>_r
            return py_result
    
    property note:
        def __set__(self,  note):
        
            self.inst.get().note = deref((convString(note)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().note
            py_result = convOutputString(_r)
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint())
    
    def _init_1(self, double in_0 , double in_1 ):
        """
        _init_1(self, in_0: float , in_1: float ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        self.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint((<double>in_0), (<double>in_1)))
    
    def _init_2(self, double in_0 , double in_1 ,  in_2 ):
        """
        _init_2(self, in_0: float , in_1: float , in_2: Union[bytes, str, String] ) -> None
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
        assert (isinstance(in_2, str) or isinstance(in_2, bytes) or isinstance(in_2, String)), 'arg in_2 wrong type'
    
    
    
        self.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint((<double>in_0), (<double>in_1), deref((convString(in_2)).get())))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: float , in_1: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: float , in_1: float , in_2: Union[bytes, str, String] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], float)):
             self._init_1(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (2, 0):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, TM_DataPoint):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef TM_DataPoint other_casted = other
        cdef TM_DataPoint self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 
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
