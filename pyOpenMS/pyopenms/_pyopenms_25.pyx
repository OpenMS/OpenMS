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

cdef class __AggregationMethod:
    """
          Aggregation method
    """
    PROD = 0
    SUM = 1
    BEST = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ChecksumType:
    None
    UNKNOWN_CHECKSUM = 0
    SHA1 = 1
    MD5 = 2
    SIZE_OF_CHECKSUMTYPE = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __OpenPepXLAlgorithm_ExitCodes:
    None
    EXECUTION_OK = 0
    ILLEGAL_PARAMETERS = 1
    UNEXPECTED_RESULT = 2
    INCOMPATIBLE_INPUT_DATA = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AccurateMassSearchResult:
    """
    Cython implementation of _AccurateMassSearchResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AccurateMassSearchResult.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_AccurateMassSearchResult](new _AccurateMassSearchResult())
    
    def getObservedMZ(self):
        """
        getObservedMZ(self) -> float
        """
        cdef double _r = self.inst.get().getObservedMZ()
        py_result = <double>_r
        return py_result
    
    def setObservedMZ(self, double m ):
        """
        setObservedMZ(self, m: float ) -> None
        """
        assert isinstance(m, float), 'arg m wrong type'
    
        self.inst.get().setObservedMZ((<double &>m))
    
    def getCalculatedMZ(self):
        """
        getCalculatedMZ(self) -> float
        """
        cdef double _r = self.inst.get().getCalculatedMZ()
        py_result = <double>_r
        return py_result
    
    def setCalculatedMZ(self, double m ):
        """
        setCalculatedMZ(self, m: float ) -> None
        """
        assert isinstance(m, float), 'arg m wrong type'
    
        self.inst.get().setCalculatedMZ((<double &>m))
    
    def getQueryMass(self):
        """
        getQueryMass(self) -> float
        """
        cdef double _r = self.inst.get().getQueryMass()
        py_result = <double>_r
        return py_result
    
    def setQueryMass(self, double m ):
        """
        setQueryMass(self, m: float ) -> None
        """
        assert isinstance(m, float), 'arg m wrong type'
    
        self.inst.get().setQueryMass((<double &>m))
    
    def getFoundMass(self):
        """
        getFoundMass(self) -> float
        """
        cdef double _r = self.inst.get().getFoundMass()
        py_result = <double>_r
        return py_result
    
    def setFoundMass(self, double m ):
        """
        setFoundMass(self, m: float ) -> None
        """
        assert isinstance(m, float), 'arg m wrong type'
    
        self.inst.get().setFoundMass((<double &>m))
    
    def getCharge(self):
        """
        getCharge(self) -> float
        """
        cdef double _r = self.inst.get().getCharge()
        py_result = <double>_r
        return py_result
    
    def setCharge(self, double ch ):
        """
        setCharge(self, ch: float ) -> None
        """
        assert isinstance(ch, float), 'arg ch wrong type'
    
        self.inst.get().setCharge((<double &>ch))
    
    def getMZErrorPPM(self):
        """
        getMZErrorPPM(self) -> float
        """
        cdef double _r = self.inst.get().getMZErrorPPM()
        py_result = <double>_r
        return py_result
    
    def setMZErrorPPM(self, double ppm ):
        """
        setMZErrorPPM(self, ppm: float ) -> None
        """
        assert isinstance(ppm, float), 'arg ppm wrong type'
    
        self.inst.get().setMZErrorPPM((<double &>ppm))
    
    def getObservedRT(self):
        """
        getObservedRT(self) -> float
        """
        cdef double _r = self.inst.get().getObservedRT()
        py_result = <double>_r
        return py_result
    
    def setObservedRT(self, double rt ):
        """
        setObservedRT(self, rt: float ) -> None
        """
        assert isinstance(rt, float), 'arg rt wrong type'
    
        self.inst.get().setObservedRT((<double &>rt))
    
    def getObservedIntensity(self):
        """
        getObservedIntensity(self) -> float
        """
        cdef double _r = self.inst.get().getObservedIntensity()
        py_result = <double>_r
        return py_result
    
    def setObservedIntensity(self, double intensity ):
        """
        setObservedIntensity(self, intensity: float ) -> None
        """
        assert isinstance(intensity, float), 'arg intensity wrong type'
    
        self.inst.get().setObservedIntensity((<double &>intensity))
    
    def getMatchingIndex(self):
        """
        getMatchingIndex(self) -> float
        """
        cdef double _r = self.inst.get().getMatchingIndex()
        py_result = <double>_r
        return py_result
    
    def setMatchingIndex(self, double idx ):
        """
        setMatchingIndex(self, idx: float ) -> None
        """
        assert isinstance(idx, float), 'arg idx wrong type'
    
        self.inst.get().setMatchingIndex((<double &>idx))
    
    def getFoundAdduct(self):
        """
        getFoundAdduct(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getFoundAdduct()
        py_result = convOutputString(_r)
        return py_result
    
    def setFoundAdduct(self,  add ):
        """
        setFoundAdduct(self, add: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(add, str) or isinstance(add, bytes) or isinstance(add, String)), 'arg add wrong type'
    
        self.inst.get().setFoundAdduct(deref((convString(add)).get()))
    
    def getFormulaString(self):
        """
        getFormulaString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getFormulaString()
        py_result = convOutputString(_r)
        return py_result
    
    def setEmpiricalFormula(self,  ep ):
        """
        setEmpiricalFormula(self, ep: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(ep, str) or isinstance(ep, bytes) or isinstance(ep, String)), 'arg ep wrong type'
    
        self.inst.get().setEmpiricalFormula(deref((convString(ep)).get()))
    
    def getMatchingHMDBids(self):
        """
        getMatchingHMDBids(self) -> List[bytes]
        """
        _r = self.inst.get().getMatchingHMDBids()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setMatchingHMDBids(self, list match_ids ):
        """
        setMatchingHMDBids(self, match_ids: List[bytes] ) -> None
        """
        assert isinstance(match_ids, list) and all(isinstance(i, bytes) for i in match_ids), 'arg match_ids wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in match_ids:
           v0.push_back(_String(<char *>item0))
        self.inst.get().setMatchingHMDBids(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        match_ids[:] = replace
        del v0
    
    def getIsotopesSimScore(self):
        """
        getIsotopesSimScore(self) -> float
        """
        cdef double _r = self.inst.get().getIsotopesSimScore()
        py_result = <double>_r
        return py_result
    
    def setIsotopesSimScore(self, double sim_score ):
        """
        setIsotopesSimScore(self, sim_score: float ) -> None
        """
        assert isinstance(sim_score, float), 'arg sim_score wrong type'
    
        self.inst.get().setIsotopesSimScore((<double &>sim_score))
    
    def getIndividualIntensities(self):
        """
        getIndividualIntensities(self) -> List[float]
        """
        _r = self.inst.get().getIndividualIntensities()
        cdef list py_result = _r
        return py_result
    
    def setIndividualIntensities(self, list in_0 ):
        """
        setIndividualIntensities(self, in_0: List[float] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] v0 = in_0
        self.inst.get().setIndividualIntensities(v0)
        
    
    def getSourceFeatureIndex(self):
        """
        getSourceFeatureIndex(self) -> int
        """
        cdef size_t _r = self.inst.get().getSourceFeatureIndex()
        py_result = <size_t>_r
        return py_result
    
    def setSourceFeatureIndex(self,  in_0 ):
        """
        setSourceFeatureIndex(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int) and in_0 >= 0, 'arg in_0 wrong type'
    
        self.inst.get().setSourceFeatureIndex((<size_t>in_0))
    
    def getMasstraceIntensities(self):
        """
        getMasstraceIntensities(self) -> List[float]
        """
        _r = self.inst.get().getMasstraceIntensities()
        cdef list py_result = _r
        return py_result
    
    def setMasstraceIntensities(self, list in_0 ):
        """
        setMasstraceIntensities(self, in_0: List[float] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, float) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[double] v0 = in_0
        self.inst.get().setMasstraceIntensities(v0)
        in_0[:] = v0 

cdef class BasicProteinInferenceAlgorithm:
    """
    Cython implementation of _BasicProteinInferenceAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BasicProteinInferenceAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Algorithm class that implements simple protein inference by aggregation of peptide scores.
    
    It has multiple parameter options like the aggregation method, when to distinguish peptidoforms,
    and if you want to use shared peptides ("use_shared_peptides").
    First, the best PSM per spectrum is used, then only the best PSM per peptidoform is aggregated.
    Peptidoforms can optionally be distinguished via the treat_X_separate parameters:
    - Modifications (modified sequence string)
    - Charge states
    The algorithm assumes posteriors or posterior error probabilities and converts to posteriors initially.
    Possible aggregation methods that can be set via the parameter "aggregation_method" are:
    - "best" (default)
    - "sum"
    - "product" (ignoring zeroes)
    Annotation of the number of peptides used for aggregation can be disabled (see parameters).
    Supports multiple runs but goes through them one by one iterating over the full PeptideIdentification vector.
    Warning: Does not "link" the peptides to the resulting protein run. If you wish to do that you have to do
    it manually.
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_BasicProteinInferenceAlgorithm](new _BasicProteinInferenceAlgorithm())
    
    def _run_0(self, PeptideIdentificationList pep_ids , list prot_ids ):
        """
        _run_0(self, pep_ids: PeptideIdentificationList , prot_ids: List[ProteinIdentification] ) -> None
        Performs basic aggregation-based inference per ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_ids: Vector of protein identification runs. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        """
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in prot_ids:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().run((deref(pep_ids.inst.get())), deref(v1))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v1.begin()
        replace_0 = []
        while it_prot_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item1)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v1
    
    def _run_1(self, PeptideIdentificationList pep_ids , ProteinIdentification prot_id ):
        """
        _run_1(self, pep_ids: PeptideIdentificationList , prot_id: ProteinIdentification ) -> None
        Performs basic aggregation-based inference on single ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        """
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
        assert isinstance(prot_id, ProteinIdentification), 'arg prot_id wrong type'
    
    
        self.inst.get().run((deref(pep_ids.inst.get())), (deref(prot_id.inst.get())))
    
    def _run_2(self, ConsensusMap cmap , ProteinIdentification prot_id , bool include_unassigned ):
        """
        _run_2(self, cmap: ConsensusMap , prot_id: ProteinIdentification , include_unassigned: bool ) -> None
        Performs basic aggregation-based inference on identifications in a ConsensusMap. See class help.\n
        `prot_id` should contain the union of all proteins in the map. E.g. use ConsensusMapMergerAlgorithm and
        then pass the first=merged run.
        
        
        :param cmap: ConsensusMap = Consensus features with metadata and peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        """
        assert isinstance(cmap, ConsensusMap), 'arg cmap wrong type'
        assert isinstance(prot_id, ProteinIdentification), 'arg prot_id wrong type'
        assert isinstance(include_unassigned, pybool_t), 'arg include_unassigned wrong type'
    
    
    
        self.inst.get().run((deref(cmap.inst.get())), (deref(prot_id.inst.get())), (<bool>include_unassigned))
    
    def run(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: run(self, pep_ids: PeptideIdentificationList , prot_ids: List[ProteinIdentification] ) -> None
          :noindex:
        
        Performs basic aggregation-based inference per ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_ids: Vector of protein identification runs. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        
        .. rubric:: Overload:
        .. py:function:: run(self, pep_ids: PeptideIdentificationList , prot_id: ProteinIdentification ) -> None
          :noindex:
        
        Performs basic aggregation-based inference on single ProteinIdentification run. See class help.
        
        
        :param pep_ids: Vector of peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
        
        .. rubric:: Overload:
        .. py:function:: run(self, cmap: ConsensusMap , prot_id: ProteinIdentification , include_unassigned: bool ) -> None
          :noindex:
        
        Performs basic aggregation-based inference on identifications in a ConsensusMap. See class help.\n
        `prot_id` should contain the union of all proteins in the map. E.g. use ConsensusMapMergerAlgorithm and
        then pass the first=merged run.
        
        
        :param cmap: ConsensusMap = Consensus features with metadata and peptide identifications
        :param prot_id: ProteinIdentification run with possible proteins. Scores will be overwritten and groups added.
        :return: Writes its results into prot_ids
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])):
            return self._run_0(*args)
        elif (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], ProteinIdentification)):
            return self._run_1(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ProteinIdentification)) and (isinstance(args[2], pybool_t)):
            return self._run_2(*args)
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
    AggregationMethod = __AggregationMethod 

cdef class GNPSMGFFile:
    """
    Cython implementation of _GNPSMGFFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GNPSMGFFile.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef GNPSMGFFile rv = GNPSMGFFile.__new__(GNPSMGFFile)
       rv.inst = shared_ptr[_GNPSMGFFile](new _GNPSMGFFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef GNPSMGFFile rv = GNPSMGFFile.__new__(GNPSMGFFile)
       rv.inst = shared_ptr[_GNPSMGFFile](new _GNPSMGFFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_GNPSMGFFile](new _GNPSMGFFile())
    
    def _init_1(self, GNPSMGFFile in_0 ):
        """
        _init_1(self, in_0: GNPSMGFFile ) -> None
        """
        assert isinstance(in_0, GNPSMGFFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_GNPSMGFFile](new _GNPSMGFFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: GNPSMGFFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], GNPSMGFFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  consensus_file_path , list mzml_file_paths ,  out ):
        """
        store(self, consensus_file_path: Union[bytes, str, String] , mzml_file_paths: List[bytes] , out: Union[bytes, str, String] ) -> None
        Export consensus file from default workflow to GNPS MGF format
        """
        assert (isinstance(consensus_file_path, str) or isinstance(consensus_file_path, bytes) or isinstance(consensus_file_path, String)), 'arg consensus_file_path wrong type'
        assert isinstance(mzml_file_paths, list) and all(isinstance(li, bytes) for li in mzml_file_paths), 'arg mzml_file_paths wrong type'
        assert (isinstance(out, str) or isinstance(out, bytes) or isinstance(out, String)), 'arg out wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in mzml_file_paths:
           v1.push_back(_String(<char *>item1))
    
        self.inst.get().store(deref((convString(consensus_file_path)).get()), deref(v1), deref((convString(out)).get()))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        mzml_file_paths[:] = replace
        del v1
    
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

cdef class InspectInfile:
    """
    Cython implementation of _InspectInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InspectInfile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef InspectInfile rv = InspectInfile.__new__(InspectInfile)
       rv.inst = shared_ptr[_InspectInfile](new _InspectInfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InspectInfile rv = InspectInfile.__new__(InspectInfile)
       rv.inst = shared_ptr[_InspectInfile](new _InspectInfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Inspect input file adapter
        """
        self.inst = shared_ptr[_InspectInfile](new _InspectInfile())
    
    def _init_1(self, InspectInfile in_0 ):
        """
        _init_1(self, in_0: InspectInfile ) -> None
        """
        assert isinstance(in_0, InspectInfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InspectInfile](new _InspectInfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Inspect input file adapter

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InspectInfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InspectInfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename ):
        """
        store(self, filename: Union[bytes, str, String] ) -> None
        Stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().store(deref((convString(filename)).get()))
    
    def handlePTMs(self,  modification_line ,  modifications_filename , bool monoisotopic ):
        """
        handlePTMs(self, modification_line: Union[bytes, str, String] , modifications_filename: Union[bytes, str, String] , monoisotopic: bool ) -> None
        Retrieves the name, mass change, affected residues, type and position for all modifications from a string
        
        
        :param modification_line:
        :param modifications_filename:
        :param monoisotopic: if true, masses are considered to be monoisotopic
        :raises:
          Exception: FileNotReadable if the modifications_filename could not be read
        :raises:
          Exception: FileNotFound if modifications_filename could not be found
        :raises:
          Exception: ParseError if modifications_filename could not be parsed
        """
        assert (isinstance(modification_line, str) or isinstance(modification_line, bytes) or isinstance(modification_line, String)), 'arg modification_line wrong type'
        assert (isinstance(modifications_filename, str) or isinstance(modifications_filename, bytes) or isinstance(modifications_filename, String)), 'arg modifications_filename wrong type'
        assert isinstance(monoisotopic, pybool_t), 'arg monoisotopic wrong type'
    
    
    
        self.inst.get().handlePTMs(deref((convString(modification_line)).get()), deref((convString(modifications_filename)).get()), (<bool>monoisotopic))
    
    def getSpectra(self):
        """
        getSpectra(self) -> Union[bytes, str, String]
        Specifies a spectrum file to search
        """
        cdef _String _r = self.inst.get().getSpectra()
        py_result = convOutputString(_r)
        return py_result
    
    def setSpectra(self,  spectra ):
        """
        setSpectra(self, spectra: Union[bytes, str, String] ) -> None
        Specifies a spectrum file to search
        """
        assert (isinstance(spectra, str) or isinstance(spectra, bytes) or isinstance(spectra, String)), 'arg spectra wrong type'
    
        self.inst.get().setSpectra(deref((convString(spectra)).get()))
    
    def getDb(self):
        """
        getDb(self) -> Union[bytes, str, String]
        Specifies the name of a database (.trie file) to search
        """
        cdef _String _r = self.inst.get().getDb()
        py_result = convOutputString(_r)
        return py_result
    
    def setDb(self,  db ):
        """
        setDb(self, db: Union[bytes, str, String] ) -> None
        Specifies the name of a database (.trie file) to search
        """
        assert (isinstance(db, str) or isinstance(db, bytes) or isinstance(db, String)), 'arg db wrong type'
    
        self.inst.get().setDb(deref((convString(db)).get()))
    
    def getEnzyme(self):
        """
        getEnzyme(self) -> Union[bytes, str, String]
        Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        """
        cdef _String _r = self.inst.get().getEnzyme()
        py_result = convOutputString(_r)
        return py_result
    
    def setEnzyme(self,  enzyme ):
        """
        setEnzyme(self, enzyme: Union[bytes, str, String] ) -> None
        Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
        """
        assert (isinstance(enzyme, str) or isinstance(enzyme, bytes) or isinstance(enzyme, String)), 'arg enzyme wrong type'
    
        self.inst.get().setEnzyme(deref((convString(enzyme)).get()))
    
    def getModificationsPerPeptide(self):
        """
        getModificationsPerPeptide(self) -> int
        Number of PTMs permitted in a single peptide
        """
        cdef int _r = self.inst.get().getModificationsPerPeptide()
        py_result = <int>_r
        return py_result
    
    def setModificationsPerPeptide(self,  modifications_per_peptide ):
        """
        setModificationsPerPeptide(self, modifications_per_peptide: int ) -> None
        Number of PTMs permitted in a single peptide
        """
        assert isinstance(modifications_per_peptide, int), 'arg modifications_per_peptide wrong type'
    
        self.inst.get().setModificationsPerPeptide((<int>modifications_per_peptide))
    
    def getBlind(self):
        """
        getBlind(self) -> int
        Run inspect in a blind mode
        """
        cdef unsigned int _r = self.inst.get().getBlind()
        py_result = <unsigned int>_r
        return py_result
    
    def setBlind(self,  blind ):
        """
        setBlind(self, blind: int ) -> None
        Run inspect in a blind mode
        """
        assert isinstance(blind, int), 'arg blind wrong type'
    
        self.inst.get().setBlind((<unsigned int>blind))
    
    def getMaxPTMsize(self):
        """
        getMaxPTMsize(self) -> float
        The maximum modification size (in Da) to consider in a blind search
        """
        cdef float _r = self.inst.get().getMaxPTMsize()
        py_result = <float>_r
        return py_result
    
    def setMaxPTMsize(self, float maxptmsize ):
        """
        setMaxPTMsize(self, maxptmsize: float ) -> None
        The maximum modification size (in Da) to consider in a blind search
        """
        assert isinstance(maxptmsize, float), 'arg maxptmsize wrong type'
    
        self.inst.get().setMaxPTMsize((<float>maxptmsize))
    
    def getPrecursorMassTolerance(self):
        """
        getPrecursorMassTolerance(self) -> float
        Specifies the parent mass tolerance, in Daltons
        """
        cdef float _r = self.inst.get().getPrecursorMassTolerance()
        py_result = <float>_r
        return py_result
    
    def setPrecursorMassTolerance(self, float precursor_mass_tolerance ):
        """
        setPrecursorMassTolerance(self, precursor_mass_tolerance: float ) -> None
        Specifies the parent mass tolerance, in Daltons
        """
        assert isinstance(precursor_mass_tolerance, float), 'arg precursor_mass_tolerance wrong type'
    
        self.inst.get().setPrecursorMassTolerance((<float>precursor_mass_tolerance))
    
    def getPeakMassTolerance(self):
        """
        getPeakMassTolerance(self) -> float
        How far b and y peaks can be shifted from their expected masses.
        """
        cdef float _r = self.inst.get().getPeakMassTolerance()
        py_result = <float>_r
        return py_result
    
    def setPeakMassTolerance(self, float peak_mass_tolerance ):
        """
        setPeakMassTolerance(self, peak_mass_tolerance: float ) -> None
        How far b and y peaks can be shifted from their expected masses
        """
        assert isinstance(peak_mass_tolerance, float), 'arg peak_mass_tolerance wrong type'
    
        self.inst.get().setPeakMassTolerance((<float>peak_mass_tolerance))
    
    def getMulticharge(self):
        """
        getMulticharge(self) -> int
        If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        """
        cdef unsigned int _r = self.inst.get().getMulticharge()
        py_result = <unsigned int>_r
        return py_result
    
    def setMulticharge(self,  multicharge ):
        """
        setMulticharge(self, multicharge: int ) -> None
        If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible
        """
        assert isinstance(multicharge, int), 'arg multicharge wrong type'
    
        self.inst.get().setMulticharge((<unsigned int>multicharge))
    
    def getInstrument(self):
        """
        getInstrument(self) -> Union[bytes, str, String]
        If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        """
        cdef _String _r = self.inst.get().getInstrument()
        py_result = convOutputString(_r)
        return py_result
    
    def setInstrument(self,  instrument ):
        """
        setInstrument(self, instrument: Union[bytes, str, String] ) -> None
        If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass
        """
        assert (isinstance(instrument, str) or isinstance(instrument, bytes) or isinstance(instrument, String)), 'arg instrument wrong type'
    
        self.inst.get().setInstrument(deref((convString(instrument)).get()))
    
    def getTagCount(self):
        """
        getTagCount(self) -> int
        Number of tags to generate
        """
        cdef int _r = self.inst.get().getTagCount()
        py_result = <int>_r
        return py_result
    
    def setTagCount(self,  TagCount ):
        """
        setTagCount(self, TagCount: int ) -> None
        Number of tags to generate
        """
        assert isinstance(TagCount, int), 'arg TagCount wrong type'
    
        self.inst.get().setTagCount((<int>TagCount))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, InspectInfile):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef InspectInfile other_casted = other
        cdef InspectInfile self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
    
    def getModifications(self):
        _r = self.inst.get().getModifications()
        py_result = dict()
        cdef libcpp_map[_String, libcpp_vector[_String]].iterator outer_it = _r.begin()
        cdef libcpp_vector[_String].iterator inner_it
        cdef String item_0
        cdef bytes inner_key
        cdef list inner_values
        while outer_it != _r.end():
           inner_key = deref(outer_it).first.c_str()
           inner_values = []
           inner_it = deref(outer_it).second.begin()
           while inner_it != deref(outer_it).second.end():
               # item_0 = CVTerm.__new__(CVTerm)
               # item_0.inst = shared_ptr[_CVTerm](new _CVTerm(deref(inner_it)))
               inner_values.append(  deref(inner_it).c_str() )
               inc(inner_it)
           py_result[inner_key] = inner_values
           inc(outer_it)
        return py_result 

cdef class IsotopePattern:
    """
    Cython implementation of _IsotopePattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1IsotopePattern.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property spectrum:
        def __set__(self, list spectrum):
            cdef libcpp_vector[size_t] v0 = spectrum
            self.inst.get().spectrum = v0
            
    
        def __get__(self):
            _r = self.inst.get().spectrum
            cdef list py_result = _r
            return py_result
    
    property intensity:
        def __set__(self, list intensity):
            cdef libcpp_vector[double] v0 = intensity
            self.inst.get().intensity = v0
            
    
        def __get__(self):
            _r = self.inst.get().intensity
            cdef list py_result = _r
            return py_result
    
    property mz_score:
        def __set__(self, list mz_score):
            cdef libcpp_vector[double] v0 = mz_score
            self.inst.get().mz_score = v0
            
    
        def __get__(self):
            _r = self.inst.get().mz_score
            cdef list py_result = _r
            return py_result
    
    property theoretical_mz:
        def __set__(self, list theoretical_mz):
            cdef libcpp_vector[double] v0 = theoretical_mz
            self.inst.get().theoretical_mz = v0
            
    
        def __get__(self):
            _r = self.inst.get().theoretical_mz
            cdef list py_result = _r
            return py_result
    
    property theoretical_pattern:
        def __set__(self, TheoreticalIsotopePattern theoretical_pattern):
        
            self.inst.get().theoretical_pattern = (deref(theoretical_pattern.inst.get()))
        
    
        def __get__(self):
            cdef _TheoreticalIsotopePattern * _r = new _TheoreticalIsotopePattern(self.inst.get().theoretical_pattern)
            cdef TheoreticalIsotopePattern py_result = TheoreticalIsotopePattern.__new__(TheoreticalIsotopePattern)
            py_result.inst = shared_ptr[_TheoreticalIsotopePattern](_r)
            return py_result
    
    def __copy__(self):
       cdef IsotopePattern rv = IsotopePattern.__new__(IsotopePattern)
       rv.inst = shared_ptr[_IsotopePattern](new _IsotopePattern(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopePattern rv = IsotopePattern.__new__(IsotopePattern)
       rv.inst = shared_ptr[_IsotopePattern](new _IsotopePattern(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  size ):
        """
        _init_0(self, size: int ) -> None
        """
        assert isinstance(size, int) and size >= 0, 'arg size wrong type'
    
        self.inst = shared_ptr[_IsotopePattern](new _IsotopePattern((<size_t>size)))
    
    def _init_1(self, IsotopePattern in_0 ):
        """
        _init_1(self, in_0: IsotopePattern ) -> None
        """
        assert isinstance(in_0, IsotopePattern), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopePattern](new _IsotopePattern((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, size: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopePattern ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopePattern)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MapAlignmentEvaluationAlgorithmRecall:
    """
    Cython implementation of _MapAlignmentEvaluationAlgorithmRecall

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentEvaluationAlgorithmRecall.html>`_
      -- Inherits from ['MapAlignmentEvaluationAlgorithm']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MapAlignmentEvaluationAlgorithmRecall](new _MapAlignmentEvaluationAlgorithmRecall()) 

cdef class MassTrace:
    """
    Cython implementation of _MassTrace

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1MassTrace.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property max_rt:
        def __set__(self, double max_rt):
        
            self.inst.get().max_rt = (<double>max_rt)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().max_rt
            py_result = <double>_r
            return py_result
    
    property theoretical_int:
        def __set__(self, double theoretical_int):
        
            self.inst.get().theoretical_int = (<double>theoretical_int)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().theoretical_int
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef MassTrace rv = MassTrace.__new__(MassTrace)
       rv.inst = shared_ptr[_MassTrace](new _MassTrace(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassTrace rv = MassTrace.__new__(MassTrace)
       rv.inst = shared_ptr[_MassTrace](new _MassTrace(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MassTrace](new _MassTrace())
    
    def _init_1(self, MassTrace in_0 ):
        """
        _init_1(self, in_0: MassTrace ) -> None
        """
        assert isinstance(in_0, MassTrace), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassTrace](new _MassTrace((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassTrace ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassTrace)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getConvexhull(self):
        """
        getConvexhull(self) -> ConvexHull2D
        """
        cdef _ConvexHull2D * _r = new _ConvexHull2D(self.inst.get().getConvexhull())
        cdef ConvexHull2D py_result = ConvexHull2D.__new__(ConvexHull2D)
        py_result.inst = shared_ptr[_ConvexHull2D](_r)
        return py_result
    
    def updateMaximum(self):
        """
        updateMaximum(self) -> None
        """
        self.inst.get().updateMaximum()
    
    def getAvgMZ(self):
        """
        getAvgMZ(self) -> float
        """
        cdef double _r = self.inst.get().getAvgMZ()
        py_result = <double>_r
        return py_result
    
    def isValid(self):
        """
        isValid(self) -> bool
        """
        cdef bool _r = self.inst.get().isValid()
        py_result = <bool>_r
        return py_result 

cdef class MassTraces:
    """
    Cython implementation of _MassTraces

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1MassTraces.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property max_trace:
        def __set__(self,  max_trace):
        
            self.inst.get().max_trace = (<size_t>max_trace)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().max_trace
            py_result = <size_t>_r
            return py_result
    
    property baseline:
        def __set__(self, double baseline):
        
            self.inst.get().baseline = (<double>baseline)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().baseline
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef MassTraces rv = MassTraces.__new__(MassTraces)
       rv.inst = shared_ptr[_MassTraces](new _MassTraces(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MassTraces rv = MassTraces.__new__(MassTraces)
       rv.inst = shared_ptr[_MassTraces](new _MassTraces(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MassTraces](new _MassTraces())
    
    def _init_1(self, MassTraces in_0 ):
        """
        _init_1(self, in_0: MassTraces ) -> None
        """
        assert isinstance(in_0, MassTraces), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MassTraces](new _MassTraces((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MassTraces ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MassTraces)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPeakCount(self):
        """
        getPeakCount(self) -> int
        """
        cdef size_t _r = self.inst.get().getPeakCount()
        py_result = <size_t>_r
        return py_result
    
    def isValid(self, double seed_mz , double trace_tolerance ):
        """
        isValid(self, seed_mz: float , trace_tolerance: float ) -> bool
        """
        assert isinstance(seed_mz, float), 'arg seed_mz wrong type'
        assert isinstance(trace_tolerance, float), 'arg trace_tolerance wrong type'
    
    
        cdef bool _r = self.inst.get().isValid((<double>seed_mz), (<double>trace_tolerance))
        py_result = <bool>_r
        return py_result
    
    def getTheoreticalmaxPosition(self):
        """
        getTheoreticalmaxPosition(self) -> int
        """
        cdef size_t _r = self.inst.get().getTheoreticalmaxPosition()
        py_result = <size_t>_r
        return py_result
    
    def updateBaseline(self):
        """
        updateBaseline(self) -> None
        """
        self.inst.get().updateBaseline()
    
    def getRTBounds(self):
        """
        getRTBounds(self) -> List[float, float]
        """
        _r = self.inst.get().getRTBounds()
        cdef list py_result = [_r.first, _r.second]
        return py_result 

cdef class MetaboTargetedTargetDecoy:
    """
    Cython implementation of _MetaboTargetedTargetDecoy

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedTargetDecoy.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboTargetedTargetDecoy rv = MetaboTargetedTargetDecoy.__new__(MetaboTargetedTargetDecoy)
       rv.inst = shared_ptr[_MetaboTargetedTargetDecoy](new _MetaboTargetedTargetDecoy(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboTargetedTargetDecoy rv = MetaboTargetedTargetDecoy.__new__(MetaboTargetedTargetDecoy)
       rv.inst = shared_ptr[_MetaboTargetedTargetDecoy](new _MetaboTargetedTargetDecoy(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Resolve overlapping fragments and missing decoys for experimental specific decoy generation in targeted/pseudo targeted metabolomics
        """
        self.inst = shared_ptr[_MetaboTargetedTargetDecoy](new _MetaboTargetedTargetDecoy())
    
    def _init_1(self, MetaboTargetedTargetDecoy in_0 ):
        """
        _init_1(self, in_0: MetaboTargetedTargetDecoy ) -> None
        """
        assert isinstance(in_0, MetaboTargetedTargetDecoy), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboTargetedTargetDecoy](new _MetaboTargetedTargetDecoy((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Resolve overlapping fragments and missing decoys for experimental specific decoy generation in targeted/pseudo targeted metabolomics

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboTargetedTargetDecoy ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboTargetedTargetDecoy)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def constructTargetDecoyMassMapping(self, TargetedExperiment t_exp ):
        """
        constructTargetDecoyMassMapping(self, t_exp: TargetedExperiment ) -> List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping]
        Constructs a mass mapping of targets and decoys using the unique m_id identifier
        
        
        :param t_exp: TransitionExperiment holds compound and transition information used for the mapping
        """
        assert isinstance(t_exp, TargetedExperiment), 'arg t_exp wrong type'
    
        _r = self.inst.get().constructTargetDecoyMassMapping((deref(t_exp.inst.get())))
        py_result = []
        cdef libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping].iterator it__r = _r.begin()
        cdef MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping item_py_result
        while it__r != _r.end():
           item_py_result = MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.__new__(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)
           item_py_result.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def resolveOverlappingTargetDecoyMassesByDecoyMassShift(self, TargetedExperiment t_exp , list mappings , double mass_to_add , double mz_tol ,  mz_tol_unit ):
        """
        resolveOverlappingTargetDecoyMassesByDecoyMassShift(self, t_exp: TargetedExperiment , mappings: List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] , mass_to_add: float , mz_tol: float , mz_tol_unit: String ) -> None
        Resolves overlapping target and decoy transition masses by adding a specifiable mass (e.g. CH2) to the overlapping decoy fragment
        
        
        :param t_exp: TransitionExperiment holds compound and transition information
        :param mappings: Map of identifier to target and decoy masses
        :param mass_to_add: (e.g. CH2)
        :param mz_tol: m/z tolerarance for target and decoy transition masses to be considered overlapping
        :param mz_tol_unit: m/z tolerance unit
        """
        assert isinstance(t_exp, TargetedExperiment), 'arg t_exp wrong type'
        assert isinstance(mappings, list) and all(isinstance(elemt_rec, MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping) for elemt_rec in mappings), 'arg mappings wrong type'
        assert isinstance(mass_to_add, float), 'arg mass_to_add wrong type'
        assert isinstance(mz_tol, float), 'arg mz_tol wrong type'
        assert isinstance(mz_tol_unit, String), 'arg mz_tol_unit wrong type'
    
        cdef libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] * v1 = new libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping]()
        cdef MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping item1
        for item1 in mappings:
            v1.push_back(deref(item1.inst.get()))
    
    
    
        self.inst.get().resolveOverlappingTargetDecoyMassesByDecoyMassShift((deref(t_exp.inst.get())), deref(v1), (<double &>mass_to_add), (<double &>mz_tol), deref((<String>mz_tol_unit).inst.get()))
        cdef libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping].iterator it_mappings = v1.begin()
        replace_0 = []
        while it_mappings != v1.end():
            item1 = MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.__new__(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)
            item1.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(deref(it_mappings)))
            replace_0.append(item1)
            inc(it_mappings)
        mappings[:] = replace_0
        del v1
    
    def generateMissingDecoysByMassShift(self, TargetedExperiment t_exp , list mappings , double mass_to_add ):
        """
        generateMissingDecoysByMassShift(self, t_exp: TargetedExperiment , mappings: List[MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] , mass_to_add: float ) -> None
        Generate a decoy for targets where fragmentation tree re-rooting was not possible, by adding a specifiable mass to the target fragments
        
        
        :param t_exp: TransitionExperiment holds compound and transition information
        :param mappings: Map of identifier to target and decoy masses
        :param mass_to_add: The maximum number of transitions required per assay
        """
        assert isinstance(t_exp, TargetedExperiment), 'arg t_exp wrong type'
        assert isinstance(mappings, list) and all(isinstance(elemt_rec, MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping) for elemt_rec in mappings), 'arg mappings wrong type'
        assert isinstance(mass_to_add, float), 'arg mass_to_add wrong type'
    
        cdef libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping] * v1 = new libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping]()
        cdef MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping item1
        for item1 in mappings:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().generateMissingDecoysByMassShift((deref(t_exp.inst.get())), deref(v1), (<double &>mass_to_add))
        cdef libcpp_vector[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping].iterator it_mappings = v1.begin()
        replace_0 = []
        while it_mappings != v1.end():
            item1 = MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.__new__(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)
            item1.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(deref(it_mappings)))
            replace_0.append(item1)
            inc(it_mappings)
        mappings[:] = replace_0
        del v1 

cdef class MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping:
    """
    Cython implementation of _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping rv = MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.__new__(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)
       rv.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping rv = MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping.__new__(MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)
       rv.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping())
    
    def _init_1(self, MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping in_0 ):
        """
        _init_1(self, in_0: MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ) -> None
        """
        assert isinstance(in_0, MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping](new _MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboTargetedTargetDecoy_MetaboTargetDecoyMassMapping)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MzMLSqliteHandler:
    """
    Cython implementation of _MzMLSqliteHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1MzMLSqliteHandler.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzMLSqliteHandler rv = MzMLSqliteHandler.__new__(MzMLSqliteHandler)
       rv.inst = shared_ptr[_MzMLSqliteHandler](new _MzMLSqliteHandler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzMLSqliteHandler rv = MzMLSqliteHandler.__new__(MzMLSqliteHandler)
       rv.inst = shared_ptr[_MzMLSqliteHandler](new _MzMLSqliteHandler(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  filename ,  run_id ):
        """
        _init_0(self, filename: Union[bytes, str, String] , run_id: int ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(run_id, int) and run_id >= 0, 'arg run_id wrong type'
    
    
        self.inst = shared_ptr[_MzMLSqliteHandler](new _MzMLSqliteHandler(deref((convString(filename)).get()), (<uint64_t>run_id)))
    
    def _init_1(self, MzMLSqliteHandler in_0 ):
        """
        _init_1(self, in_0: MzMLSqliteHandler ) -> None
        """
        assert isinstance(in_0, MzMLSqliteHandler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzMLSqliteHandler](new _MzMLSqliteHandler((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , run_id: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzMLSqliteHandler ) -> None
          :noindex:
    
        """
        if (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int) and args[1] >= 0):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzMLSqliteHandler)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def readExperiment(self, MSExperiment exp , bool meta_only ):
        """
        readExperiment(self, exp: MSExperiment , meta_only: bool ) -> None
        Read an experiment into an MSExperiment structure
        
        
        :param exp: The result data structure
        :param meta_only: Only read the meta data
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(meta_only, pybool_t), 'arg meta_only wrong type'
    
    
        self.inst.get().readExperiment((deref(exp.inst.get())), (<bool>meta_only))
    
    def readSpectra(self, list exp , list indices , bool meta_only ):
        """
        readSpectra(self, exp: List[MSSpectrum] , indices: List[int] , meta_only: bool ) -> None
        Read a set of spectra (potentially restricted to a subset)
        
        
        :param exp: The result data structure
        :param indices: A list of indices restricting the resulting spectra only to those specified here
        :param meta_only: Only read the meta data
        """
        assert isinstance(exp, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in exp), 'arg exp wrong type'
        assert isinstance(indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in indices), 'arg indices wrong type'
        assert isinstance(meta_only, pybool_t), 'arg meta_only wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in exp:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[int] v1 = indices
    
        self.inst.get().readSpectra(deref(v0), v1, (<bool>meta_only))
        
        cdef libcpp_vector[_MSSpectrum].iterator it_exp = v0.begin()
        replace_0 = []
        while it_exp != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_exp)))
            replace_0.append(item0)
            inc(it_exp)
        exp[:] = replace_0
        del v0
    
    def readChromatograms(self, list exp , list indices , bool meta_only ):
        """
        readChromatograms(self, exp: List[MSChromatogram] , indices: List[int] , meta_only: bool ) -> None
        Read a set of chromatograms (potentially restricted to a subset)
        
        
        :param exp: The result data structure
        :param indices: A list of indices restricting the resulting spectra only to those specified here
        :param meta_only: Only read the meta data
        """
        assert isinstance(exp, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in exp), 'arg exp wrong type'
        assert isinstance(indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in indices), 'arg indices wrong type'
        assert isinstance(meta_only, pybool_t), 'arg meta_only wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in exp:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[int] v1 = indices
    
        self.inst.get().readChromatograms(deref(v0), v1, (<bool>meta_only))
        
        cdef libcpp_vector[_MSChromatogram].iterator it_exp = v0.begin()
        replace_0 = []
        while it_exp != v0.end():
            item0 = MSChromatogram.__new__(MSChromatogram)
            item0.inst = shared_ptr[_MSChromatogram](new _MSChromatogram(deref(it_exp)))
            replace_0.append(item0)
            inc(it_exp)
        exp[:] = replace_0
        del v0
    
    def getNrSpectra(self):
        """
        getNrSpectra(self) -> int
        Returns number of spectra in the file, reutrns the number of spectra
        """
        cdef size_t _r = self.inst.get().getNrSpectra()
        py_result = <size_t>_r
        return py_result
    
    def getNrChromatograms(self):
        """
        getNrChromatograms(self) -> int
        Returns the number of chromatograms in the file
        """
        cdef size_t _r = self.inst.get().getNrChromatograms()
        py_result = <size_t>_r
        return py_result
    
    def setConfig(self, bool write_full_meta , bool use_lossy_compression , double linear_abs_mass_acc ):
        """
        setConfig(self, write_full_meta: bool , use_lossy_compression: bool , linear_abs_mass_acc: float ) -> None
        Sets file configuration
        
        
        :param write_full_meta: Whether to write a complete mzML meta data structure into the RUN_EXTRA field (allows complete recovery of the input file)
        :param use_lossy_compression: Whether to use lossy compression (ms numpress)
        :param linear_abs_mass_acc: Accepted loss in mass accuracy (absolute m/z, in Th)
        """
        assert isinstance(write_full_meta, pybool_t), 'arg write_full_meta wrong type'
        assert isinstance(use_lossy_compression, pybool_t), 'arg use_lossy_compression wrong type'
        assert isinstance(linear_abs_mass_acc, float), 'arg linear_abs_mass_acc wrong type'
    
    
    
        self.inst.get().setConfig((<bool>write_full_meta), (<bool>use_lossy_compression), (<double>linear_abs_mass_acc))
    
    def getSpectraIndicesbyRT(self, double RT , double deltaRT , list indices ):
        """
        getSpectraIndicesbyRT(self, RT: float , deltaRT: float , indices: List[int] ) -> List[int]
        Returns spectral indices around a specific retention time
        
        :param RT: The retention time
        :param deltaRT: Tolerance window around RT (if less or equal than zero, only the first spectrum *after* RT is returned)
        :param indices: Spectra to consider (if empty, all spectra are considered)
        :return: The indices of the spectra within RT +/- deltaRT
        """
        assert isinstance(RT, float), 'arg RT wrong type'
        assert isinstance(deltaRT, float), 'arg deltaRT wrong type'
        assert isinstance(indices, list) and all(isinstance(elemt_rec, int) for elemt_rec in indices), 'arg indices wrong type'
    
    
        cdef libcpp_vector[int] v2 = indices
        _r = self.inst.get().getSpectraIndicesbyRT((<double>RT), (<double>deltaRT), v2)
        
        cdef list py_result = _r
        return py_result
    
    def writeExperiment(self, MSExperiment exp ):
        """
        writeExperiment(self, exp: MSExperiment ) -> None
        Write an MSExperiment to disk
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().writeExperiment((deref(exp.inst.get())))
    
    def createTables(self):
        """
        createTables(self) -> None
        Create data tables for a new file
        """
        self.inst.get().createTables()
    
    def writeSpectra(self, list spectra ):
        """
        writeSpectra(self, spectra: List[MSSpectrum] ) -> None
        Writes a set of spectra to disk
        """
        assert isinstance(spectra, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in spectra), 'arg spectra wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in spectra:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().writeSpectra(deref(v0))
        del v0
    
    def writeChromatograms(self, list chroms ):
        """
        writeChromatograms(self, chroms: List[MSChromatogram] ) -> None
        Writes a set of chromatograms to disk
        """
        assert isinstance(chroms, list) and all(isinstance(elemt_rec, MSChromatogram) for elemt_rec in chroms), 'arg chroms wrong type'
        cdef libcpp_vector[_MSChromatogram] * v0 = new libcpp_vector[_MSChromatogram]()
        cdef MSChromatogram item0
        for item0 in chroms:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().writeChromatograms(deref(v0))
        del v0
    
    def writeRunLevelInformation(self, MSExperiment exp , bool write_full_meta ):
        """
        writeRunLevelInformation(self, exp: MSExperiment , write_full_meta: bool ) -> None
        Write the run-level information for an experiment into tables
        
        This is a low level function, do not call this function unless you know what you are doing
        
        
        :param exp: The result data structure
        :param meta_only: Only read the meta data
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(write_full_meta, pybool_t), 'arg write_full_meta wrong type'
    
    
        self.inst.get().writeRunLevelInformation((deref(exp.inst.get())), (<bool>write_full_meta))
    
    def getRunID(self):
        """
        getRunID(self) -> int
        Extract the `RUN` ID from the sqMass file
        """
        cdef uint64_t _r = self.inst.get().getRunID()
        py_result = <uint64_t>_r
        return py_result 

cdef class OpenPepXLAlgorithm:
    """
    Cython implementation of _OpenPepXLAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenPepXLAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenPepXLAlgorithm rv = OpenPepXLAlgorithm.__new__(OpenPepXLAlgorithm)
       rv.inst = shared_ptr[_OpenPepXLAlgorithm](new _OpenPepXLAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenPepXLAlgorithm rv = OpenPepXLAlgorithm.__new__(OpenPepXLAlgorithm)
       rv.inst = shared_ptr[_OpenPepXLAlgorithm](new _OpenPepXLAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenPepXLAlgorithm](new _OpenPepXLAlgorithm())
    
    def _init_1(self, OpenPepXLAlgorithm in_0 ):
        """
        _init_1(self, in_0: OpenPepXLAlgorithm ) -> None
        """
        assert isinstance(in_0, OpenPepXLAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenPepXLAlgorithm](new _OpenPepXLAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenPepXLAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenPepXLAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, MSExperiment unprocessed_spectra , ConsensusMap cfeatures , list fasta_db , list protein_ids , PeptideIdentificationList peptide_ids , OPXL_PreprocessedPairSpectra preprocessed_pair_spectra , list spectrum_pairs , list all_top_csms , MSExperiment spectra ):
        """
        run(self, unprocessed_spectra: MSExperiment , cfeatures: ConsensusMap , fasta_db: List[FASTAEntry] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , preprocessed_pair_spectra: OPXL_PreprocessedPairSpectra , spectrum_pairs: List[List[int, int]] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment ) -> int
        Performs the main function of this class, the search for cross-linked peptides
        
        
        :param unprocessed_spectra: The input PeakMap of experimental spectra
        :param cfeatures: The input cfeatures
        :param fasta_db: The protein database containing targets and decoys
        :param protein_ids: A result vector containing search settings. Should contain one PeptideIdentification
        :param peptide_ids: A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty
        :param preprocessed_pair_spectra: A result structure containing linear and cross-linked ion spectra. Will be overwritten. This is only necessary for writing out xQuest type spectrum files
        :param spectrum_pairs: A result vector containing paired spectra indices. Should be empty. This is only necessary for writing out xQuest type spectrum files
        :param all_top_csms: A result vector containing cross-link spectrum matches as CrossLinkSpectrumMatches. Should be empty. This is only necessary for writing out xQuest type spectrum files
        :param spectra: A result vector containing the input spectra after preprocessing and filtering. Should be empty. This is only necessary for writing out xQuest type spectrum files
        """
        assert isinstance(unprocessed_spectra, MSExperiment), 'arg unprocessed_spectra wrong type'
        assert isinstance(cfeatures, ConsensusMap), 'arg cfeatures wrong type'
        assert isinstance(fasta_db, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in fasta_db), 'arg fasta_db wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(preprocessed_pair_spectra, OPXL_PreprocessedPairSpectra), 'arg preprocessed_pair_spectra wrong type'
        assert isinstance(spectrum_pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in spectrum_pairs), 'arg spectrum_pairs wrong type'
        assert isinstance(all_top_csms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in all_top_csms), 'arg all_top_csms wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
        cdef libcpp_vector[_FASTAEntry] * v2 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item2
        for item2 in fasta_db:
            v2.push_back(deref(item2.inst.get()))
        cdef libcpp_vector[_ProteinIdentification] * v3 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item3
        for item3 in protein_ids:
            v3.push_back(deref(item3.inst.get()))
    
    
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v6 = spectrum_pairs
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]] * v7 = new libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]]()
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v7_rec = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item7_rec
        cdef libcpp_vector[_CrossLinkSpectrumMatch].iterator it_all_top_csms_rec
        for all_top_csms_rec in all_top_csms:
            v7_rec.clear()
            for item7_rec in all_top_csms_rec:
                v7_rec.push_back(deref(item7_rec.inst.get()))
            v7.push_back(deref(v7_rec))
    
        cdef _OpenPepXLAlgorithm_ExitCodes _r = self.inst.get().run((deref(unprocessed_spectra.inst.get())), (deref(cfeatures.inst.get())), deref(v2), deref(v3), (deref(peptide_ids.inst.get())), (deref(preprocessed_pair_spectra.inst.get())), v6, deref(v7), (deref(spectra.inst.get())))
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]].iterator it_all_top_csms = v7.begin()
        replace_0 = []
        while it_all_top_csms != v7.end():
            it_all_top_csms_rec = deref(it_all_top_csms).begin()
            replace_1 = []
            while it_all_top_csms_rec != deref(it_all_top_csms).end():
                item7_rec = CrossLinkSpectrumMatch.__new__(CrossLinkSpectrumMatch)
                item7_rec.inst = shared_ptr[_CrossLinkSpectrumMatch](new _CrossLinkSpectrumMatch(deref(it_all_top_csms_rec)))
                replace_1.append(item7_rec)
                inc(it_all_top_csms_rec)
            replace_0.append(replace_1)
            inc(it_all_top_csms)
        all_top_csms[:] = replace_0
        del v7
        spectrum_pairs[:] = v6
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v3.begin()
        replace_0 = []
        while it_protein_ids != v3.end():
            item3 = ProteinIdentification.__new__(ProteinIdentification)
            item3.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item3)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v3
        cdef libcpp_vector[_FASTAEntry].iterator it_fasta_db = v2.begin()
        replace_0 = []
        while it_fasta_db != v2.end():
            item2 = FASTAEntry.__new__(FASTAEntry)
            item2.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_fasta_db)))
            replace_0.append(item2)
            inc(it_fasta_db)
        fasta_db[:] = replace_0
        del v2
        py_result = <int>_r
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
    OpenPepXLAlgorithm_ExitCodes = __OpenPepXLAlgorithm_ExitCodes 

cdef class PeakWidthEstimator:
    """
    Cython implementation of _PeakWidthEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakWidthEstimator.html>`_

    Rough estimation of the peak width at m/z
    
    Based on the peaks of the dataset (peak position & width) and the peak
    boundaries as reported by the PeakPickerHiRes, the typical peak width is
    estimated for arbitrary m/z using a spline interpolationThis struct can be used to store both peak or feature indices`
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakWidthEstimator rv = PeakWidthEstimator.__new__(PeakWidthEstimator)
       rv.inst = shared_ptr[_PeakWidthEstimator](new _PeakWidthEstimator(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakWidthEstimator rv = PeakWidthEstimator.__new__(PeakWidthEstimator)
       rv.inst = shared_ptr[_PeakWidthEstimator](new _PeakWidthEstimator(deref(self.inst.get())))
       return rv
    
    def _init_0(self, PeakWidthEstimator in_0 ):
        """
        _init_0(self, in_0: PeakWidthEstimator ) -> None
        """
        assert isinstance(in_0, PeakWidthEstimator), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakWidthEstimator](new _PeakWidthEstimator((deref(in_0.inst.get()))))
    
    def _init_1(self, MSExperiment exp_picked , list boundaries ):
        """
        _init_1(self, exp_picked: MSExperiment , boundaries: List[List[PeakBoundary]] ) -> None
        """
        assert isinstance(exp_picked, MSExperiment), 'arg exp_picked wrong type'
        assert isinstance(boundaries, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in boundaries), 'arg boundaries wrong type'
    
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]] * v1 = new libcpp_vector[libcpp_vector[_PeakBoundary]]()
        cdef libcpp_vector[_PeakBoundary] * v1_rec = new libcpp_vector[_PeakBoundary]()
        cdef PeakBoundary item1_rec
        cdef libcpp_vector[_PeakBoundary].iterator it_boundaries_rec
        for boundaries_rec in boundaries:
            v1_rec.clear()
            for item1_rec in boundaries_rec:
                v1_rec.push_back(deref(item1_rec.inst.get()))
            v1.push_back(deref(v1_rec))
        self.inst = shared_ptr[_PeakWidthEstimator](new _PeakWidthEstimator((deref(exp_picked.inst.get())), deref(v1)))
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]].iterator it_boundaries = v1.begin()
        replace_0 = []
        while it_boundaries != v1.end():
            it_boundaries_rec = deref(it_boundaries).begin()
            replace_1 = []
            while it_boundaries_rec != deref(it_boundaries).end():
                item1_rec = PeakBoundary.__new__(PeakBoundary)
                item1_rec.inst = shared_ptr[_PeakBoundary](new _PeakBoundary(deref(it_boundaries_rec)))
                replace_1.append(item1_rec)
                inc(it_boundaries_rec)
            replace_0.append(replace_1)
            inc(it_boundaries)
        boundaries[:] = replace_0
        del v1
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakWidthEstimator ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, exp_picked: MSExperiment , boundaries: List[List[PeakBoundary]] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], PeakWidthEstimator)):
             self._init_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in args[1])):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getPeakWidth(self, double mz ):
        """
        getPeakWidth(self, mz: float ) -> float
        Returns the estimated peak width at m/z
        """
        assert isinstance(mz, float), 'arg mz wrong type'
    
        cdef double _r = self.inst.get().getPeakWidth((<double>mz))
        py_result = <double>_r
        return py_result 

cdef class Product:
    """
    Cython implementation of _Product

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Product.html>`_

    This class describes the product isolation window for special scan types, such as MRM
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Product rv = Product.__new__(Product)
       rv.inst = shared_ptr[_Product](new _Product(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Product rv = Product.__new__(Product)
       rv.inst = shared_ptr[_Product](new _Product(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Product](new _Product())
    
    def _init_1(self, Product in_0 ):
        """
        _init_1(self, in_0: Product ) -> None
        """
        assert isinstance(in_0, Product), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Product](new _Product((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Product ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Product)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMZ(self):
        """
        getMZ(self) -> float
        Returns the target m/z
        """
        cdef double _r = self.inst.get().getMZ()
        py_result = <double>_r
        return py_result
    
    def setMZ(self, double in_0 ):
        """
        setMZ(self, in_0: float ) -> None
        Sets the target m/z
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setMZ((<double>in_0))
    
    def getIsolationWindowLowerOffset(self):
        """
        getIsolationWindowLowerOffset(self) -> float
        Returns the lower offset from the target m/z
        """
        cdef double _r = self.inst.get().getIsolationWindowLowerOffset()
        py_result = <double>_r
        return py_result
    
    def setIsolationWindowLowerOffset(self, double bound ):
        """
        setIsolationWindowLowerOffset(self, bound: float ) -> None
        Sets the lower offset from the target m/z
        """
        assert isinstance(bound, float), 'arg bound wrong type'
    
        self.inst.get().setIsolationWindowLowerOffset((<double>bound))
    
    def getIsolationWindowUpperOffset(self):
        """
        getIsolationWindowUpperOffset(self) -> float
        Returns the upper offset from the target m/z
        """
        cdef double _r = self.inst.get().getIsolationWindowUpperOffset()
        py_result = <double>_r
        return py_result
    
    def setIsolationWindowUpperOffset(self, double bound ):
        """
        setIsolationWindowUpperOffset(self, bound: float ) -> None
        Sets the upper offset from the target m/z
        """
        assert isinstance(bound, float), 'arg bound wrong type'
    
        self.inst.get().setIsolationWindowUpperOffset((<double>bound))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Product):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Product other_casted = other
        cdef Product self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class RankScaler:
    """
    Cython implementation of _RankScaler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RankScaler.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_RankScaler](new _RankScaler())
    
    def filterSpectrum(self, MSSpectrum spec ):
        """
        filterSpectrum(self, spec: MSSpectrum ) -> None
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterSpectrum((deref(spec.inst.get())))
    
    def filterPeakSpectrum(self, MSSpectrum spec ):
        """
        filterPeakSpectrum(self, spec: MSSpectrum ) -> None
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterPeakSpectrum((deref(spec.inst.get())))
    
    def filterPeakMap(self, MSExperiment exp ):
        """
        filterPeakMap(self, exp: MSExperiment ) -> None
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
        self.inst.get().filterPeakMap((deref(exp.inst.get())))
    
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

cdef class Seed:
    """
    Cython implementation of _Seed

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1Seed.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property spectrum:
        def __set__(self,  spectrum):
        
            self.inst.get().spectrum = (<size_t>spectrum)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().spectrum
            py_result = <size_t>_r
            return py_result
    
    property peak:
        def __set__(self,  peak):
        
            self.inst.get().peak = (<size_t>peak)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().peak
            py_result = <size_t>_r
            return py_result
    
    property intensity:
        def __set__(self, float intensity):
        
            self.inst.get().intensity = (<float>intensity)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().intensity
            py_result = <float>_r
            return py_result
    
    def __copy__(self):
       cdef Seed rv = Seed.__new__(Seed)
       rv.inst = shared_ptr[_Seed](new _Seed(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Seed rv = Seed.__new__(Seed)
       rv.inst = shared_ptr[_Seed](new _Seed(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Seed](new _Seed())
    
    def _init_1(self, Seed in_0 ):
        """
        _init_1(self, in_0: Seed ) -> None
        """
        assert isinstance(in_0, Seed), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Seed](new _Seed((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Seed ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Seed)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def __richcmp__(self, other, op):
        if op not in (0,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Seed):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Seed other_casted = other
        cdef Seed self_casted = self
        if op==0:
            return deref(self_casted.inst.get()) < deref(other_casted.inst.get()) 

cdef class SourceFile:
    """
    Cython implementation of _SourceFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SourceFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SourceFile rv = SourceFile.__new__(SourceFile)
       rv.inst = shared_ptr[_SourceFile](new _SourceFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SourceFile rv = SourceFile.__new__(SourceFile)
       rv.inst = shared_ptr[_SourceFile](new _SourceFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Description of a file location, used to store the origin of (meta) data
        """
        self.inst = shared_ptr[_SourceFile](new _SourceFile())
    
    def _init_1(self, SourceFile in_0 ):
        """
        _init_1(self, in_0: SourceFile ) -> None
        """
        assert isinstance(in_0, SourceFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SourceFile](new _SourceFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Description of a file location, used to store the origin of (meta) data

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SourceFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SourceFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getNameOfFile(self):
        """
        getNameOfFile(self) -> Union[bytes, str, String]
        Returns the file name
        """
        cdef _String _r = self.inst.get().getNameOfFile()
        py_result = convOutputString(_r)
        return py_result
    
    def setNameOfFile(self,  in_0 ):
        """
        setNameOfFile(self, in_0: Union[bytes, str, String] ) -> None
        Sets the file name
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setNameOfFile(deref((convString(in_0)).get()))
    
    def getPathToFile(self):
        """
        getPathToFile(self) -> Union[bytes, str, String]
        Returns the file path
        """
        cdef _String _r = self.inst.get().getPathToFile()
        py_result = convOutputString(_r)
        return py_result
    
    def setPathToFile(self,  in_0 ):
        """
        setPathToFile(self, in_0: Union[bytes, str, String] ) -> None
        Sets the file path
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setPathToFile(deref((convString(in_0)).get()))
    
    def getFileSize(self):
        """
        getFileSize(self) -> float
        Returns the file size in MB
        """
        cdef float _r = self.inst.get().getFileSize()
        py_result = <float>_r
        return py_result
    
    def setFileSize(self, float in_0 ):
        """
        setFileSize(self, in_0: float ) -> None
        Sets the file size in MB
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
    
        self.inst.get().setFileSize((<float>in_0))
    
    def getFileType(self):
        """
        getFileType(self) -> Union[bytes, str, String]
        Returns the file type
        """
        cdef _String _r = self.inst.get().getFileType()
        py_result = convOutputString(_r)
        return py_result
    
    def setFileType(self,  in_0 ):
        """
        setFileType(self, in_0: Union[bytes, str, String] ) -> None
        Sets the file type
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setFileType(deref((convString(in_0)).get()))
    
    def getChecksum(self):
        """
        getChecksum(self) -> Union[bytes, str, String]
        Returns the file's checksum
        """
        cdef _String _r = self.inst.get().getChecksum()
        py_result = convOutputString(_r)
        return py_result
    
    def setChecksum(self,  in_0 , int in_1 ):
        """
        setChecksum(self, in_0: Union[bytes, str, String] , in_1: int ) -> None
        Sets the file's checksum
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert in_1 in [0, 1, 2, 3], 'arg in_1 wrong type'
    
    
        self.inst.get().setChecksum(deref((convString(in_0)).get()), (<_ChecksumType>in_1))
    
    def getChecksumType(self):
        """
        getChecksumType(self) -> int
        Returns the checksum type
        """
        cdef _ChecksumType _r = self.inst.get().getChecksumType()
        py_result = <int>_r
        return py_result
    
    def getNativeIDType(self):
        """
        getNativeIDType(self) -> Union[bytes, str, String]
        Returns the native ID type of the spectra
        """
        cdef _String _r = self.inst.get().getNativeIDType()
        py_result = convOutputString(_r)
        return py_result
    
    def setNativeIDType(self,  in_0 ):
        """
        setNativeIDType(self, in_0: Union[bytes, str, String] ) -> None
        Sets the native ID type of the spectra
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().setNativeIDType(deref((convString(in_0)).get()))
    
    def getNativeIDTypeAccession(self):
        """
        getNativeIDTypeAccession(self) -> Union[bytes, str, String]
        Returns the nativeID of the spectra
        """
        cdef _String _r = self.inst.get().getNativeIDTypeAccession()
        py_result = convOutputString(_r)
        return py_result
    
    def setNativeIDTypeAccession(self,  accesssion ):
        """
        setNativeIDTypeAccession(self, accesssion: Union[bytes, str, String] ) -> None
        Sets the native ID of the spectra
        """
        assert (isinstance(accesssion, str) or isinstance(accesssion, bytes) or isinstance(accesssion, String)), 'arg accesssion wrong type'
    
        self.inst.get().setNativeIDTypeAccession(deref((convString(accesssion)).get())) 

cdef class Tagger:
    """
    Cython implementation of _Tagger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Tagger.html>`_

    Constructor for Tagger
    
    The parameter `max_charge_` should be >= `min_charge_`
    Also `max_tag_length` should be >= `min_tag_length`
    
    :param min_tag_length: The minimal sequence tag length
    :param ppm: The tolerance for matching residue masses to peak delta masses
    :param max_tag_length: The maximal sequence tag length
    :param min_charge: Minimal fragment charge considered for each sequence tag
    :param max_charge: Maximal fragment charge considered for each sequence tag
    :param fixed_mods: A list of modification names. The modified residues replace the unmodified versions
    :param var_mods: A list of modification names. The modified residues are added as additional entries to the list of residues
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Tagger rv = Tagger.__new__(Tagger)
       rv.inst = shared_ptr[_Tagger](new _Tagger(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Tagger rv = Tagger.__new__(Tagger)
       rv.inst = shared_ptr[_Tagger](new _Tagger(deref(self.inst.get())))
       return rv
    
    def _init_0(self, Tagger in_0 ):
        """
        _init_0(self, in_0: Tagger ) -> None
        """
        assert isinstance(in_0, Tagger), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Tagger](new _Tagger((deref(in_0.inst.get()))))
    
    def _init_1(self,  min_tag_length , double ppm ,  max_tag_length ,  min_charge ,  max_charge , list fixed_mods , list var_mods ):
        """
        _init_1(self, min_tag_length: int , ppm: float , max_tag_length: int , min_charge: int , max_charge: int , fixed_mods: List[bytes] , var_mods: List[bytes] ) -> None
        """
        assert isinstance(min_tag_length, int) and min_tag_length >= 0, 'arg min_tag_length wrong type'
        assert isinstance(ppm, float), 'arg ppm wrong type'
        assert isinstance(max_tag_length, int) and max_tag_length >= 0, 'arg max_tag_length wrong type'
        assert isinstance(min_charge, int) and min_charge >= 0, 'arg min_charge wrong type'
        assert isinstance(max_charge, int) and max_charge >= 0, 'arg max_charge wrong type'
        assert isinstance(fixed_mods, list) and all(isinstance(li, bytes) for li in fixed_mods), 'arg fixed_mods wrong type'
        assert isinstance(var_mods, list) and all(isinstance(li, bytes) for li in var_mods), 'arg var_mods wrong type'
    
    
    
    
    
        cdef libcpp_vector[_String] * v5 = new libcpp_vector[_String]()
        cdef bytes item5
        for item5 in fixed_mods:
           v5.push_back(_String(<char *>item5))
        cdef libcpp_vector[_String] * v6 = new libcpp_vector[_String]()
        cdef bytes item6
        for item6 in var_mods:
           v6.push_back(_String(<char *>item6))
        self.inst = shared_ptr[_Tagger](new _Tagger((<size_t>min_tag_length), (<double>ppm), (<size_t>max_tag_length), (<size_t>min_charge), (<size_t>max_charge), deref(v5), deref(v6)))
        replace = []
        cdef libcpp_vector[_String].iterator it_6 = v6.begin()
        while it_6 != v6.end():
           replace.append(<char*>deref(it_6).c_str())
           inc(it_6)
        var_mods[:] = replace
        del v6
        replace = []
        cdef libcpp_vector[_String].iterator it_5 = v5.begin()
        while it_5 != v5.end():
           replace.append(<char*>deref(it_5).c_str())
           inc(it_5)
        fixed_mods[:] = replace
        del v5
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Tagger ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, min_tag_length: int , ppm: float , max_tag_length: int , min_charge: int , max_charge: int , fixed_mods: List[bytes] , var_mods: List[bytes] ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], Tagger)):
             self._init_0(*args)
        elif (len(args)==7) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], float)) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], int) and args[3] >= 0) and (isinstance(args[4], int) and args[4] >= 0) and (isinstance(args[5], list) and all(isinstance(li, bytes) for li in args[5])) and (isinstance(args[6], list) and all(isinstance(li, bytes) for li in args[6])):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getTag_0(self, list mzs , list tags ):
        """
        _getTag_0(self, mzs: List[float] , tags: List[Union[bytes, str]] ) -> None
        Generate tags from mass vector `mzs`
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param mzs: A vector of mz values, containing the mz values from a centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
        """
        assert isinstance(mzs, list) and all(isinstance(elemt_rec, float) for elemt_rec in mzs), 'arg mzs wrong type'
        assert isinstance(tags, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tags), 'arg tags wrong type'
        cdef libcpp_vector[double] v0 = mzs
        cdef libcpp_vector[libcpp_utf8_string] v1 = tags
        self.inst.get().getTag(v0, v1)
        tags[:] = v1
        
    
    def _getTag_1(self, MSSpectrum spec , list tags ):
        """
        _getTag_1(self, spec: MSSpectrum , tags: List[Union[bytes, str]] ) -> None
        Generate tags from an MSSpectrum
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param spec: A centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(tags, list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in tags), 'arg tags wrong type'
    
        cdef libcpp_vector[libcpp_utf8_string] v1 = tags
        self.inst.get().getTag((deref(spec.inst.get())), v1)
        tags[:] = v1
    
    def getTag(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getTag(self, mzs: List[float] , tags: List[Union[bytes, str]] ) -> None
          :noindex:
        
        Generate tags from mass vector `mzs`
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param mzs: A vector of mz values, containing the mz values from a centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
        
        .. rubric:: Overload:
        .. py:function:: getTag(self, spec: MSSpectrum , tags: List[Union[bytes, str]] ) -> None
          :noindex:
        
        Generate tags from an MSSpectrum
        
        The parameter `tags` is filled with one string per sequence tag
        It uses the standard residues from ResidueDB including
        the fixed and variable modifications given to the constructor
        
        :param spec: A centroided fragment spectrum
        :param tags: The vector of tags, that is filled with this function
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in args[1])):
            return self._getTag_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, (bytes, str)) for elemt_rec in args[1])):
            return self._getTag_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMaxCharge(self,  max_charge ):
        """
        setMaxCharge(self, max_charge: int ) -> None
        Change the maximal charge considered by the tagger
        
        Allows to change the maximal considered charge e.g. based on a spectra
        precursor charge without calling the constructor multiple times
        
        :param max_charge: The new maximal charge
        """
        assert isinstance(max_charge, int) and max_charge >= 0, 'arg max_charge wrong type'
    
        self.inst.get().setMaxCharge((<size_t>max_charge)) 

cdef class TheoreticalIsotopePattern:
    """
    Cython implementation of _TheoreticalIsotopePattern

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::FeatureFinderAlgorithmPickedHelperStructs_1_1TheoreticalIsotopePattern.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property intensity:
        def __set__(self, list intensity):
            cdef libcpp_vector[double] v0 = intensity
            self.inst.get().intensity = v0
            
    
        def __get__(self):
            _r = self.inst.get().intensity
            cdef list py_result = _r
            return py_result
    
    property optional_begin:
        def __set__(self,  optional_begin):
        
            self.inst.get().optional_begin = (<size_t>optional_begin)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().optional_begin
            py_result = <size_t>_r
            return py_result
    
    property optional_end:
        def __set__(self,  optional_end):
        
            self.inst.get().optional_end = (<size_t>optional_end)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().optional_end
            py_result = <size_t>_r
            return py_result
    
    property max:
        def __set__(self, double max):
        
            self.inst.get().max = (<double>max)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().max
            py_result = <double>_r
            return py_result
    
    property trimmed_left:
        def __set__(self,  trimmed_left):
        
            self.inst.get().trimmed_left = (<size_t>trimmed_left)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().trimmed_left
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef TheoreticalIsotopePattern rv = TheoreticalIsotopePattern.__new__(TheoreticalIsotopePattern)
       rv.inst = shared_ptr[_TheoreticalIsotopePattern](new _TheoreticalIsotopePattern(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TheoreticalIsotopePattern rv = TheoreticalIsotopePattern.__new__(TheoreticalIsotopePattern)
       rv.inst = shared_ptr[_TheoreticalIsotopePattern](new _TheoreticalIsotopePattern(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TheoreticalIsotopePattern](new _TheoreticalIsotopePattern())
    
    def _init_1(self, TheoreticalIsotopePattern in_0 ):
        """
        _init_1(self, in_0: TheoreticalIsotopePattern ) -> None
        """
        assert isinstance(in_0, TheoreticalIsotopePattern), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TheoreticalIsotopePattern](new _TheoreticalIsotopePattern((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TheoreticalIsotopePattern ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TheoreticalIsotopePattern)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
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
 

## Add this to all modules
class Interfaces:

    BinaryDataArray = _Interfaces_BinaryDataArray
    Spectrum = _Interfaces_Spectrum
    Chromatogram = _Interfaces_Chromatogram 

# Add this to all modules
PeakSpectrum = MSSpectrum

PeakMap = MSExperiment 
