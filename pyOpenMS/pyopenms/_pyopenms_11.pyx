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
def __static_OpenMSOSInfo_getBinaryArchitecture():
    """
    __static_OpenMSOSInfo_getBinaryArchitecture() -> Union[bytes, str, String]
    """
    cdef _String _r = _getBinaryArchitecture_BuildInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_OpenMSBuildInfo_getBuildType():
    """
    __static_OpenMSBuildInfo_getBuildType() -> Union[bytes, str, String]
    """
    cdef _String _r = _getBuildType_BuildInfo()
    py_result = convOutputString(_r)
    return py_result

def __static_TransformationModelLowess_getDefaultParameters(Param params ):
    """
    __static_TransformationModelLowess_getDefaultParameters(params: Param ) -> None
    """
    assert isinstance(params, Param), 'arg params wrong type'

    _getDefaultParameters_TransformationModelLowess((deref(params.inst.get())))

def __static_OpenMSOSInfo_getOSInfo():
    """
    __static_OpenMSOSInfo_getOSInfo() -> OpenMSOSInfo
    """
    cdef _OpenMSOSInfo * _r = new _OpenMSOSInfo(_getOSInfo_BuildInfo())
    cdef OpenMSOSInfo py_result = OpenMSOSInfo.__new__(OpenMSOSInfo)
    py_result.inst = shared_ptr[_OpenMSOSInfo](_r)
    return py_result

def __static_OpenMSBuildInfo_getOpenMPMaxNumThreads():
    """
    __static_OpenMSBuildInfo_getOpenMPMaxNumThreads() -> int
    """
    cdef size_t _r = _getOpenMPMaxNumThreads_BuildInfo()
    py_result = <size_t>_r
    return py_result

def __static_OpenMSBuildInfo_isOpenMPEnabled():
    """
    __static_OpenMSBuildInfo_isOpenMPEnabled() -> bool
    """
    cdef bool _r = _isOpenMPEnabled_BuildInfo()
    py_result = <bool>_r
    return py_result

def __static_OpenMSBuildInfo_setOpenMPNumThreads( num_threads ):
    """
    __static_OpenMSBuildInfo_setOpenMPNumThreads(num_threads: int ) -> None
    """
    assert isinstance(num_threads, int), 'arg num_threads wrong type'

    _setOpenMPNumThreads_BuildInfo((<int>num_threads)) 

cdef class ConsensusIDAlgorithmPEPIons:
    """
    Cython implementation of _ConsensusIDAlgorithmPEPIons

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusIDAlgorithmPEPIons.html>`_
      -- Inherits from ['ConsensusIDAlgorithmSimilarity']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusIDAlgorithmPEPIons](new _ConsensusIDAlgorithmPEPIons())
    
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

cdef class MRMFeatureFilter:
    """
    Cython implementation of _MRMFeatureFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureFilter.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMFeatureFilter rv = MRMFeatureFilter.__new__(MRMFeatureFilter)
       rv.inst = shared_ptr[_MRMFeatureFilter](new _MRMFeatureFilter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMFeatureFilter rv = MRMFeatureFilter.__new__(MRMFeatureFilter)
       rv.inst = shared_ptr[_MRMFeatureFilter](new _MRMFeatureFilter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMFeatureFilter](new _MRMFeatureFilter())
    
    def _init_1(self, MRMFeatureFilter in_0 ):
        """
        _init_1(self, in_0: MRMFeatureFilter ) -> None
        """
        assert isinstance(in_0, MRMFeatureFilter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMFeatureFilter](new _MRMFeatureFilter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMFeatureFilter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMFeatureFilter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def FilterFeatureMap(self, FeatureMap features , MRMFeatureQC filter_criteria , TargetedExperiment transitions ):
        """
        FilterFeatureMap(self, features: FeatureMap , filter_criteria: MRMFeatureQC , transitions: TargetedExperiment ) -> None
        Flags or filters features and subordinates in a FeatureMap
        
        
        :param features: FeatureMap to flag or filter
        :param filter_criteria: MRMFeatureQC class defining QC parameters
        :param transitions: Transitions from a TargetedExperiment
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        assert isinstance(filter_criteria, MRMFeatureQC), 'arg filter_criteria wrong type'
        assert isinstance(transitions, TargetedExperiment), 'arg transitions wrong type'
    
    
    
        self.inst.get().FilterFeatureMap((deref(features.inst.get())), (deref(filter_criteria.inst.get())), (deref(transitions.inst.get())))
    
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

cdef class MRMIonSeries:
    """
    Cython implementation of _MRMIonSeries

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMIonSeries.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMIonSeries rv = MRMIonSeries.__new__(MRMIonSeries)
       rv.inst = shared_ptr[_MRMIonSeries](new _MRMIonSeries(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMIonSeries rv = MRMIonSeries.__new__(MRMIonSeries)
       rv.inst = shared_ptr[_MRMIonSeries](new _MRMIonSeries(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMIonSeries](new _MRMIonSeries())
    
    def _init_1(self, MRMIonSeries in_0 ):
        """
        _init_1(self, in_0: MRMIonSeries ) -> None
        """
        assert isinstance(in_0, MRMIonSeries), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMIonSeries](new _MRMIonSeries((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMIonSeries ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMIonSeries)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def annotateTransitionCV(self, ReactionMonitoringTransition tr ,  annotation ):
        """
        annotateTransitionCV(self, tr: ReactionMonitoringTransition , annotation: Union[bytes, str, String] ) -> None
        Annotates transition with CV terms
        
        
        :param tr: The transition to annotate
        :param annotation: The fragment ion annotation
        """
        assert isinstance(tr, ReactionMonitoringTransition), 'arg tr wrong type'
        assert (isinstance(annotation, str) or isinstance(annotation, bytes) or isinstance(annotation, String)), 'arg annotation wrong type'
    
    
        self.inst.get().annotateTransitionCV((deref(tr.inst.get())), deref((convString(annotation)).get()))
    
    def annotateTransition(self, ReactionMonitoringTransition tr , Peptide peptide , double precursor_mz_threshold , double product_mz_threshold , bool enable_reannotation , list fragment_types , list fragment_charges , bool enable_specific_losses , bool enable_unspecific_losses ,  round_decPow ):
        """
        annotateTransition(self, tr: ReactionMonitoringTransition , peptide: Peptide , precursor_mz_threshold: float , product_mz_threshold: float , enable_reannotation: bool , fragment_types: List[bytes] , fragment_charges: List[int] , enable_specific_losses: bool , enable_unspecific_losses: bool , round_decPow: int ) -> None
        Annotates transition
        
        
        :param tr: The transition to annotate
        :param peptide: The corresponding peptide
        :param precursor_mz_threshold: The m/z threshold for annotation of the precursor ion
        :param product_mz_threshold: The m/z threshold for annotation of the fragment ion
        :param enable_reannotation: Whether the original (e.g. SpectraST) annotation should be used or reannotation should be conducted
        :param fragment_types: The fragment ion types for reannotation
        :param fragment_charges: The fragment ion charges for reannotation
        :param enable_specific_losses: Whether specific neutral losses should be considered
        :param enable_unspecific_losses: Whether unspecific neutral losses (H2O1, H3N1, C1H2N2, C1H2N1O1) should be considered
        :param round_decPow: Round precursor and product m/z values to decimal power (default: -4)
        """
        assert isinstance(tr, ReactionMonitoringTransition), 'arg tr wrong type'
        assert isinstance(peptide, Peptide), 'arg peptide wrong type'
        assert isinstance(precursor_mz_threshold, float), 'arg precursor_mz_threshold wrong type'
        assert isinstance(product_mz_threshold, float), 'arg product_mz_threshold wrong type'
        assert isinstance(enable_reannotation, pybool_t), 'arg enable_reannotation wrong type'
        assert isinstance(fragment_types, list) and all(isinstance(i, bytes) for i in fragment_types), 'arg fragment_types wrong type'
        assert isinstance(fragment_charges, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in fragment_charges), 'arg fragment_charges wrong type'
        assert isinstance(enable_specific_losses, pybool_t), 'arg enable_specific_losses wrong type'
        assert isinstance(enable_unspecific_losses, pybool_t), 'arg enable_unspecific_losses wrong type'
        assert isinstance(round_decPow, int), 'arg round_decPow wrong type'
    
    
    
    
    
        cdef libcpp_vector[_String] * v5 = new libcpp_vector[_String]()
        cdef bytes item5
        for item5 in fragment_types:
           v5.push_back(_String(<char *>item5))
        cdef libcpp_vector[size_t] v6 = fragment_charges
    
    
    
        self.inst.get().annotateTransition((deref(tr.inst.get())), (deref(peptide.inst.get())), (<double>precursor_mz_threshold), (<double>product_mz_threshold), (<bool>enable_reannotation), deref(v5), v6, (<bool>enable_specific_losses), (<bool>enable_unspecific_losses), (<int>round_decPow))
        
        del v5 

cdef class MS2File:
    """
    Cython implementation of _MS2File

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MS2File.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MS2File rv = MS2File.__new__(MS2File)
       rv.inst = shared_ptr[_MS2File](new _MS2File(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MS2File rv = MS2File.__new__(MS2File)
       rv.inst = shared_ptr[_MS2File](new _MS2File(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MS2File](new _MS2File())
    
    def _init_1(self, MS2File in_0 ):
        """
        _init_1(self, in_0: MS2File ) -> None
        """
        assert isinstance(in_0, MS2File), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MS2File](new _MS2File((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MS2File ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MS2File)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , MSExperiment exp ):
        """
        load(self, filename: Union[bytes, str, String] , exp: MSExperiment ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(exp.inst.get())))
    
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

cdef class OMSSACSVFile:
    """
    Cython implementation of _OMSSACSVFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OMSSACSVFile.html>`_

    File adapter for OMSSACSV files
    
    The files contain the results of the OMSSA algorithm in a comma separated manner. This file adapter is able to
    load the data from such a file into the structures of OpenMS
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OMSSACSVFile rv = OMSSACSVFile.__new__(OMSSACSVFile)
       rv.inst = shared_ptr[_OMSSACSVFile](new _OMSSACSVFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OMSSACSVFile rv = OMSSACSVFile.__new__(OMSSACSVFile)
       rv.inst = shared_ptr[_OMSSACSVFile](new _OMSSACSVFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OMSSACSVFile](new _OMSSACSVFile())
    
    def _init_1(self, OMSSACSVFile in_0 ):
        """
        _init_1(self, in_0: OMSSACSVFile ) -> None
        """
        assert isinstance(in_0, OMSSACSVFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OMSSACSVFile](new _OMSSACSVFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OMSSACSVFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OMSSACSVFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , ProteinIdentification protein_identification , PeptideIdentificationList id_data ):
        """
        load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList ) -> None
        Loads a OMSSA file
        
        The content of the file is stored in `features`
        
        
        :param filename: The name of the file to read from
        :param protein_identification: The protein ProteinIdentification data
        :param id_data: The peptide ids of the file
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(id_data, PeptideIdentificationList), 'arg id_data wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(protein_identification.inst.get())), (deref(id_data.inst.get()))) 

cdef class OpenMSBuildInfo:
    """
    Cython implementation of _OpenMSBuildInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1OpenMSBuildInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenMSBuildInfo rv = OpenMSBuildInfo.__new__(OpenMSBuildInfo)
       rv.inst = shared_ptr[_OpenMSBuildInfo](new _OpenMSBuildInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenMSBuildInfo rv = OpenMSBuildInfo.__new__(OpenMSBuildInfo)
       rv.inst = shared_ptr[_OpenMSBuildInfo](new _OpenMSBuildInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenMSBuildInfo](new _OpenMSBuildInfo())
    
    def _init_1(self, OpenMSBuildInfo in_0 ):
        """
        _init_1(self, in_0: OpenMSBuildInfo ) -> None
        """
        assert isinstance(in_0, OpenMSBuildInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenMSBuildInfo](new _OpenMSBuildInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenMSBuildInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenMSBuildInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    getBuildType = __static_OpenMSBuildInfo_getBuildType
    getOpenMPMaxNumThreads = __static_OpenMSBuildInfo_getOpenMPMaxNumThreads
    isOpenMPEnabled = __static_OpenMSBuildInfo_isOpenMPEnabled
    setOpenMPNumThreads = __static_OpenMSBuildInfo_setOpenMPNumThreads 

cdef class OpenMSOSInfo:
    """
    Cython implementation of _OpenMSOSInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1OpenMSOSInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OpenMSOSInfo rv = OpenMSOSInfo.__new__(OpenMSOSInfo)
       rv.inst = shared_ptr[_OpenMSOSInfo](new _OpenMSOSInfo(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OpenMSOSInfo rv = OpenMSOSInfo.__new__(OpenMSOSInfo)
       rv.inst = shared_ptr[_OpenMSOSInfo](new _OpenMSOSInfo(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OpenMSOSInfo](new _OpenMSOSInfo())
    
    def _init_1(self, OpenMSOSInfo in_0 ):
        """
        _init_1(self, in_0: OpenMSOSInfo ) -> None
        """
        assert isinstance(in_0, OpenMSOSInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OpenMSOSInfo](new _OpenMSOSInfo((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OpenMSOSInfo ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OpenMSOSInfo)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getOSAsString(self):
        """
        getOSAsString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOSAsString()
        py_result = convOutputString(_r)
        return py_result
    
    def getArchAsString(self):
        """
        getArchAsString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getArchAsString()
        py_result = convOutputString(_r)
        return py_result
    
    def getOSVersionAsString(self):
        """
        getOSVersionAsString(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getOSVersionAsString()
        py_result = convOutputString(_r)
        return py_result
    getBinaryArchitecture = __static_OpenMSOSInfo_getBinaryArchitecture
    getOSInfo = __static_OpenMSOSInfo_getOSInfo 

cdef class PeakIndex:
    """
    Cython implementation of _PeakIndex

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakIndex.html>`_

    Index of a peak or feature
    
    This struct can be used to store both peak or feature indices
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property peak:
        def __set__(self,  peak):
        
            self.inst.get().peak = (<size_t>peak)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().peak
            py_result = <size_t>_r
            return py_result
    
    property spectrum:
        def __set__(self,  spectrum):
        
            self.inst.get().spectrum = (<size_t>spectrum)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().spectrum
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef PeakIndex rv = PeakIndex.__new__(PeakIndex)
       rv.inst = shared_ptr[_PeakIndex](new _PeakIndex(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakIndex rv = PeakIndex.__new__(PeakIndex)
       rv.inst = shared_ptr[_PeakIndex](new _PeakIndex(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakIndex](new _PeakIndex())
    
    def _init_1(self, PeakIndex in_0 ):
        """
        _init_1(self, in_0: PeakIndex ) -> None
        """
        assert isinstance(in_0, PeakIndex), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakIndex](new _PeakIndex((deref(in_0.inst.get()))))
    
    def _init_2(self,  peak ):
        """
        _init_2(self, peak: int ) -> None
        """
        assert isinstance(peak, int) and peak >= 0, 'arg peak wrong type'
    
        self.inst = shared_ptr[_PeakIndex](new _PeakIndex((<size_t>peak)))
    
    def _init_3(self,  spectrum ,  peak ):
        """
        _init_3(self, spectrum: int , peak: int ) -> None
        """
        assert isinstance(spectrum, int) and spectrum >= 0, 'arg spectrum wrong type'
        assert isinstance(peak, int) and peak >= 0, 'arg peak wrong type'
    
    
        self.inst = shared_ptr[_PeakIndex](new _PeakIndex((<size_t>spectrum), (<size_t>peak)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakIndex ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, peak: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, spectrum: int , peak: int ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakIndex)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_2(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isValid(self):
        """
        isValid(self) -> bool
        Returns if the current peak ref is valid
        """
        cdef bool _r = self.inst.get().isValid()
        py_result = <bool>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        Invalidates the current index
        """
        self.inst.get().clear()
    
    def getFeature(self, FeatureMap map_ ):
        """
        getFeature(self, map_: FeatureMap ) -> Feature
        Returns the feature (or consensus feature) corresponding to this index
        
        This method is intended for arrays of features e.g. FeatureMap
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        assert isinstance(map_, FeatureMap), 'arg map_ wrong type'
    
        cdef _Feature * _r = new _Feature(self.inst.get().getFeature((deref(map_.inst.get()))))
        cdef Feature py_result = Feature.__new__(Feature)
        py_result.inst = shared_ptr[_Feature](_r)
        return py_result
    
    def getPeak(self, MSExperiment map_ ):
        """
        getPeak(self, map_: MSExperiment ) -> Peak1D
        Returns a peak corresponding to this index
        
        This method is intended for arrays of DSpectra e.g. MSExperiment
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        assert isinstance(map_, MSExperiment), 'arg map_ wrong type'
    
        cdef _Peak1D * _r = new _Peak1D(self.inst.get().getPeak((deref(map_.inst.get()))))
        cdef Peak1D py_result = Peak1D.__new__(Peak1D)
        py_result.inst = shared_ptr[_Peak1D](_r)
        return py_result
    
    def getSpectrum(self, MSExperiment map_ ):
        """
        getSpectrum(self, map_: MSExperiment ) -> MSSpectrum
        Returns a spectrum corresponding to this index
        
        This method is intended for arrays of DSpectra e.g. MSExperiment
        
        The main advantage of using this method instead accessing the data directly is that range
        check performed in debug mode
        
        :raises:
          Exception: Precondition is thrown if this index is invalid for the `map` (only in debug mode)
        """
        assert isinstance(map_, MSExperiment), 'arg map_ wrong type'
    
        cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().getSpectrum((deref(map_.inst.get()))))
        cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
        py_result.inst = shared_ptr[_MSSpectrum](_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, PeakIndex):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef PeakIndex other_casted = other
        cdef PeakIndex self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class ProtXMLFile:
    """
    Cython implementation of _ProtXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProtXMLFile.html>`_

    Used to load (storing not supported, yet) ProtXML files
    
    This class is used to load (storing not supported, yet) documents that implement
    the schema of ProtXML files
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ProtXMLFile](new _ProtXMLFile())
    
    def load(self,  filename , ProteinIdentification protein_ids , PeptideIdentification peptide_ids ):
        """
        load(self, filename: Union[bytes, str, String] , protein_ids: ProteinIdentification , peptide_ids: PeptideIdentification ) -> None
        Loads the identifications of an ProtXML file without identifier
        
        The information is read in and the information is stored in the
        corresponding variables
        
        :raises:
          Exception: FileNotFound is thrown if the file could not be found
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, ProteinIdentification), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentification), 'arg peptide_ids wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(protein_ids.inst.get())), (deref(peptide_ids.inst.get())))
    
    def store(self,  filename , ProteinIdentification protein_ids , PeptideIdentification peptide_ids ,  document_id ):
        """
        store(self, filename: Union[bytes, str, String] , protein_ids: ProteinIdentification , peptide_ids: PeptideIdentification , document_id: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, ProteinIdentification), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentification), 'arg peptide_ids wrong type'
        assert (isinstance(document_id, str) or isinstance(document_id, bytes) or isinstance(document_id, String)), 'arg document_id wrong type'
    
    
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(protein_ids.inst.get())), (deref(peptide_ids.inst.get())), deref((convString(document_id)).get())) 

cdef class ProteaseDB:
    """
    Cython implementation of _ProteaseDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteaseDB.html>`_
    """

    
    def getEnzyme(self,  name ):
        """
        getEnzyme(self, name: Union[bytes, str, String] ) -> DigestionEnzymeProtein
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef const _DigestionEnzymeProtein * __r = (self.inst.get().getEnzyme(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _DigestionEnzymeProtein * _r = new _DigestionEnzymeProtein(deref(__r))
        cdef DigestionEnzymeProtein py_result = DigestionEnzymeProtein.__new__(DigestionEnzymeProtein)
        py_result.inst = shared_ptr[_DigestionEnzymeProtein](_r)
        return py_result
    
    def getEnzymeByRegEx(self,  cleavage_regex ):
        """
        getEnzymeByRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> DigestionEnzymeProtein
        """
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
    
        cdef const _DigestionEnzymeProtein * __r = (self.inst.get().getEnzymeByRegEx(deref((convString(cleavage_regex)).get())))
        if __r == NULL:
            return None
        cdef _DigestionEnzymeProtein * _r = new _DigestionEnzymeProtein(deref(__r))
        cdef DigestionEnzymeProtein py_result = DigestionEnzymeProtein.__new__(DigestionEnzymeProtein)
        py_result.inst = shared_ptr[_DigestionEnzymeProtein](_r)
        return py_result
    
    def getAllNames(self, list all_names ):
        """
        getAllNames(self, all_names: List[bytes] ) -> None
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def getAllXTandemNames(self, list all_names ):
        """
        getAllXTandemNames(self, all_names: List[bytes] ) -> None
        Returns all the enzyme names available for XTandem
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllXTandemNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def getAllOMSSANames(self, list all_names ):
        """
        getAllOMSSANames(self, all_names: List[bytes] ) -> None
        Returns all the enzyme names available for OMSSA
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllOMSSANames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def getAllCometNames(self, list all_names ):
        """
        getAllCometNames(self, all_names: List[bytes] ) -> None
        Returns all the enzyme names available for Comet
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllCometNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def getAllMSGFNames(self, list all_names ):
        """
        getAllMSGFNames(self, all_names: List[bytes] ) -> None
        Returns all the enzyme names available for MSGFPlus
        """
        assert isinstance(all_names, list) and all(isinstance(i, bytes) for i in all_names), 'arg all_names wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in all_names:
           v0.push_back(_String(<char *>item0))
        self.inst.get().getAllMSGFNames(deref(v0))
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        all_names[:] = replace
        del v0
    
    def hasEnzyme(self,  name ):
        """
        hasEnzyme(self, name: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasEnzyme(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    def hasRegEx(self,  cleavage_regex ):
        """
        hasRegEx(self, cleavage_regex: Union[bytes, str, String] ) -> bool
        """
        assert (isinstance(cleavage_regex, str) or isinstance(cleavage_regex, bytes) or isinstance(cleavage_regex, String)), 'arg cleavage_regex wrong type'
    
        cdef bool _r = self.inst.get().hasRegEx(deref((convString(cleavage_regex)).get()))
        py_result = <bool>_r
        return py_result
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_ProteaseDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_ProteaseDB](_getInstance_ProteaseDB()) 

cdef class RansacModelLinear:
    """
    Cython implementation of _RansacModelLinear

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RansacModelLinear.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RansacModelLinear rv = RansacModelLinear.__new__(RansacModelLinear)
       rv.inst = shared_ptr[_RansacModelLinear](new _RansacModelLinear(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RansacModelLinear rv = RansacModelLinear.__new__(RansacModelLinear)
       rv.inst = shared_ptr[_RansacModelLinear](new _RansacModelLinear(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_RansacModelLinear](new _RansacModelLinear())
    
    def _init_1(self, RansacModelLinear in_0 ):
        """
        _init_1(self, in_0: RansacModelLinear ) -> None
        """
        assert isinstance(in_0, RansacModelLinear), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RansacModelLinear](new _RansacModelLinear((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RansacModelLinear ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RansacModelLinear)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class TSE_Match:
    """
    Cython implementation of _TSE_Match

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TSE_Match.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property spectrum:
        def __set__(self, MSSpectrum spectrum):
        
            self.inst.get().spectrum = (deref(spectrum.inst.get()))
        
    
        def __get__(self):
            cdef _MSSpectrum * _r = new _MSSpectrum(self.inst.get().spectrum)
            cdef MSSpectrum py_result = MSSpectrum.__new__(MSSpectrum)
            py_result.inst = shared_ptr[_MSSpectrum](_r)
            return py_result
    
    property score:
        def __set__(self, double score):
        
            self.inst.get().score = (<double>score)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().score
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef TSE_Match rv = TSE_Match.__new__(TSE_Match)
       rv.inst = shared_ptr[_TSE_Match](new _TSE_Match(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TSE_Match rv = TSE_Match.__new__(TSE_Match)
       rv.inst = shared_ptr[_TSE_Match](new _TSE_Match(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TSE_Match](new _TSE_Match())
    
    def _init_1(self, TSE_Match in_0 ):
        """
        _init_1(self, in_0: TSE_Match ) -> None
        """
        assert isinstance(in_0, TSE_Match), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TSE_Match](new _TSE_Match((deref(in_0.inst.get()))))
    
    def _init_2(self, MSSpectrum spectrum , double score ):
        """
        _init_2(self, spectrum: MSSpectrum , score: float ) -> None
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
        assert isinstance(score, float), 'arg score wrong type'
    
    
        self.inst = shared_ptr[_TSE_Match](new _TSE_Match((deref(spectrum.inst.get())), (<double>score)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TSE_Match ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, spectrum: MSSpectrum , score: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TSE_Match)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class TargetedSpectraExtractor:
    """
    Cython implementation of _TargetedSpectraExtractor

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TargetedSpectraExtractor.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TargetedSpectraExtractor rv = TargetedSpectraExtractor.__new__(TargetedSpectraExtractor)
       rv.inst = shared_ptr[_TargetedSpectraExtractor](new _TargetedSpectraExtractor(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TargetedSpectraExtractor rv = TargetedSpectraExtractor.__new__(TargetedSpectraExtractor)
       rv.inst = shared_ptr[_TargetedSpectraExtractor](new _TargetedSpectraExtractor(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TargetedSpectraExtractor](new _TargetedSpectraExtractor())
    
    def _init_1(self, TargetedSpectraExtractor in_0 ):
        """
        _init_1(self, in_0: TargetedSpectraExtractor ) -> None
        """
        assert isinstance(in_0, TargetedSpectraExtractor), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TargetedSpectraExtractor](new _TargetedSpectraExtractor((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TargetedSpectraExtractor ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TargetedSpectraExtractor)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def _annotateSpectra_0(self, list in_0 , TargetedExperiment in_1 , list in_2 , FeatureMap in_3 ):
        """
        _annotateSpectra_0(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, TargetedExperiment), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
        assert isinstance(in_3, FeatureMap), 'arg in_3 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
    
        self.inst.get().annotateSpectra(deref(v0), (deref(in_1.inst.get())), deref(v2), (deref(in_3.inst.get())))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def _annotateSpectra_1(self, list in_0 , TargetedExperiment in_1 , list in_2 ):
        """
        _annotateSpectra_1(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, TargetedExperiment), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().annotateSpectra(deref(v0), (deref(in_1.inst.get())), deref(v2))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def _annotateSpectra_2(self, list in_0 , FeatureMap in_1 , FeatureMap in_2 , list in_3 ):
        """
        _annotateSpectra_2(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
        assert isinstance(in_2, FeatureMap), 'arg in_2 wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_3), 'arg in_3 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
    
    
        cdef libcpp_vector[_MSSpectrum] * v3 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item3
        for item3 in in_3:
            v3.push_back(deref(item3.inst.get()))
        self.inst.get().annotateSpectra(deref(v0), (deref(in_1.inst.get())), (deref(in_2.inst.get())), deref(v3))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_3 = v3.begin()
        replace_0 = []
        while it_in_3 != v3.end():
            item3 = MSSpectrum.__new__(MSSpectrum)
            item3.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_3)))
            replace_0.append(item3)
            inc(it_in_3)
        in_3[:] = replace_0
        del v3
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def annotateSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: annotateSpectra(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: annotateSpectra(self, in_0: List[MSSpectrum] , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: annotateSpectra(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], TargetedExperiment)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])) and (isinstance(args[3], FeatureMap)):
            return self._annotateSpectra_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], TargetedExperiment)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])):
            return self._annotateSpectra_1(*args)
        elif (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], FeatureMap)) and (isinstance(args[2], FeatureMap)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[3])):
            return self._annotateSpectra_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def searchSpectrum(self, FeatureMap in_0 , FeatureMap in_1 , bool in_2 ):
        """
        searchSpectrum(self, in_0: FeatureMap , in_1: FeatureMap , in_2: bool ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
        assert isinstance(in_2, pybool_t), 'arg in_2 wrong type'
    
    
    
        self.inst.get().searchSpectrum((deref(in_0.inst.get())), (deref(in_1.inst.get())), (<bool>in_2))
    
    def pickSpectrum(self, MSSpectrum in_0 , MSSpectrum in_1 ):
        """
        pickSpectrum(self, in_0: MSSpectrum , in_1: MSSpectrum ) -> None
        """
        assert isinstance(in_0, MSSpectrum), 'arg in_0 wrong type'
        assert isinstance(in_1, MSSpectrum), 'arg in_1 wrong type'
    
    
        self.inst.get().pickSpectrum((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def _scoreSpectra_0(self, list in_0 , list in_1 , FeatureMap in_2 , list in_3 ):
        """
        _scoreSpectra_0(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(in_2, FeatureMap), 'arg in_2 wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_3), 'arg in_3 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_MSSpectrum] * v1 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_MSSpectrum] * v3 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item3
        for item3 in in_3:
            v3.push_back(deref(item3.inst.get()))
        self.inst.get().scoreSpectra(deref(v0), deref(v1), (deref(in_2.inst.get())), deref(v3))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_3 = v3.begin()
        replace_0 = []
        while it_in_3 != v3.end():
            item3 = MSSpectrum.__new__(MSSpectrum)
            item3.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_3)))
            replace_0.append(item3)
            inc(it_in_3)
        in_3[:] = replace_0
        del v3
        cdef libcpp_vector[_MSSpectrum].iterator it_in_1 = v1.begin()
        replace_0 = []
        while it_in_1 != v1.end():
            item1 = MSSpectrum.__new__(MSSpectrum)
            item1.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_1)))
            replace_0.append(item1)
            inc(it_in_1)
        in_1[:] = replace_0
        del v1
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def _scoreSpectra_1(self, list in_0 , list in_1 , list in_2 ):
        """
        _scoreSpectra_1(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_MSSpectrum] * v1 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().scoreSpectra(deref(v0), deref(v1), deref(v2))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
        cdef libcpp_vector[_MSSpectrum].iterator it_in_1 = v1.begin()
        replace_0 = []
        while it_in_1 != v1.end():
            item1 = MSSpectrum.__new__(MSSpectrum)
            item1.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_1)))
            replace_0.append(item1)
            inc(it_in_1)
        in_1[:] = replace_0
        del v1
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def scoreSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: scoreSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: FeatureMap , in_3: List[MSSpectrum] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: scoreSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] , in_2: List[MSSpectrum] ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[1])) and (isinstance(args[2], FeatureMap)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[3])):
            return self._scoreSpectra_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[1])) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])):
            return self._scoreSpectra_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _selectSpectra_0(self, list in_0 , FeatureMap in_1 , list in_2 , FeatureMap in_3 ):
        """
        _selectSpectra_0(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
        assert isinstance(in_3, FeatureMap), 'arg in_3 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
    
        self.inst.get().selectSpectra(deref(v0), (deref(in_1.inst.get())), deref(v2), (deref(in_3.inst.get())))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def _selectSpectra_1(self, list in_0 , list in_1 ):
        """
        _selectSpectra_1(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_1), 'arg in_1 wrong type'
        cdef libcpp_vector[_MSSpectrum] * v0 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_MSSpectrum] * v1 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().selectSpectra(deref(v0), deref(v1))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_1 = v1.begin()
        replace_0 = []
        while it_in_1 != v1.end():
            item1 = MSSpectrum.__new__(MSSpectrum)
            item1.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_1)))
            replace_0.append(item1)
            inc(it_in_1)
        in_1[:] = replace_0
        del v1
        cdef libcpp_vector[_MSSpectrum].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = MSSpectrum.__new__(MSSpectrum)
            item0.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def selectSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: selectSpectra(self, in_0: List[MSSpectrum] , in_1: FeatureMap , in_2: List[MSSpectrum] , in_3: FeatureMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: selectSpectra(self, in_0: List[MSSpectrum] , in_1: List[MSSpectrum] ) -> None
          :noindex:
    
        """
        if (len(args)==4) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], FeatureMap)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])) and (isinstance(args[3], FeatureMap)):
            return self._selectSpectra_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[1])):
            return self._selectSpectra_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _extractSpectra_0(self, MSExperiment in_0 , TargetedExperiment in_1 , list in_2 , FeatureMap in_3 , bool in_4 ):
        """
        _extractSpectra_0(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap , in_4: bool ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, TargetedExperiment), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
        assert isinstance(in_3, FeatureMap), 'arg in_3 wrong type'
        assert isinstance(in_4, pybool_t), 'arg in_4 wrong type'
    
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
    
    
        self.inst.get().extractSpectra((deref(in_0.inst.get())), (deref(in_1.inst.get())), deref(v2), (deref(in_3.inst.get())), (<bool>in_4))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
    
    def _extractSpectra_1(self, MSExperiment in_0 , TargetedExperiment in_1 , list in_2 ):
        """
        _extractSpectra_1(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, TargetedExperiment), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
    
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().extractSpectra((deref(in_0.inst.get())), (deref(in_1.inst.get())), deref(v2))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
    
    def _extractSpectra_2(self, MSExperiment in_0 , FeatureMap in_1 , list in_2 ):
        """
        _extractSpectra_2(self, in_0: MSExperiment , in_1: FeatureMap , in_2: List[MSSpectrum] ) -> None
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
        assert isinstance(in_2, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in in_2), 'arg in_2 wrong type'
    
    
        cdef libcpp_vector[_MSSpectrum] * v2 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item2
        for item2 in in_2:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().extractSpectra((deref(in_0.inst.get())), (deref(in_1.inst.get())), deref(v2))
        cdef libcpp_vector[_MSSpectrum].iterator it_in_2 = v2.begin()
        replace_0 = []
        while it_in_2 != v2.end():
            item2 = MSSpectrum.__new__(MSSpectrum)
            item2.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_in_2)))
            replace_0.append(item2)
            inc(it_in_2)
        in_2[:] = replace_0
        del v2
    
    def extractSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: extractSpectra(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] , in_3: FeatureMap , in_4: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: extractSpectra(self, in_0: MSExperiment , in_1: TargetedExperiment , in_2: List[MSSpectrum] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: extractSpectra(self, in_0: MSExperiment , in_1: FeatureMap , in_2: List[MSSpectrum] ) -> None
          :noindex:
    
        """
        if (len(args)==5) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], TargetedExperiment)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])) and (isinstance(args[3], FeatureMap)) and (isinstance(args[4], pybool_t)):
            return self._extractSpectra_0(*args)
        elif (len(args)==3) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], TargetedExperiment)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])):
            return self._extractSpectra_1(*args)
        elif (len(args)==3) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], FeatureMap)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[2])):
            return self._extractSpectra_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def constructTransitionsList(self, FeatureMap in_0 , FeatureMap in_1 , TargetedExperiment in_2 ):
        """
        constructTransitionsList(self, in_0: FeatureMap , in_1: FeatureMap , in_2: TargetedExperiment ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
        assert isinstance(in_2, TargetedExperiment), 'arg in_2 wrong type'
    
    
    
        self.inst.get().constructTransitionsList((deref(in_0.inst.get())), (deref(in_1.inst.get())), (deref(in_2.inst.get())))
    
    def storeSpectraMSP(self,  in_0 , MSExperiment in_1 ):
        """
        storeSpectraMSP(self, in_0: Union[bytes, str, String] , in_1: MSExperiment ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
        assert isinstance(in_1, MSExperiment), 'arg in_1 wrong type'
    
    
        self.inst.get().storeSpectraMSP(deref((convString(in_0)).get()), (deref(in_1.inst.get())))
    
    def mergeFeatures(self, FeatureMap in_0 , FeatureMap in_1 ):
        """
        mergeFeatures(self, in_0: FeatureMap , in_1: FeatureMap ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, FeatureMap), 'arg in_1 wrong type'
    
    
        self.inst.get().mergeFeatures((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
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

cdef class TransformationModelInterpolated:
    """
    Cython implementation of _TransformationModelInterpolated

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelInterpolated.html>`_
      -- Inherits from ['TransformationModel']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list data , Param params ):
        """
        __init__(self, data: List[TM_DataPoint] , params: Param ) -> None
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        assert isinstance(params, Param), 'arg params wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_TransformationModelInterpolated](new _TransformationModelInterpolated(deref(v0), (deref(params.inst.get()))))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        Gets the default parameters
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def evaluate(self, double value ):
        """
        evaluate(self, value: float ) -> float
        Evaluate the interpolation model at the given value
        
        :param value: The position where the interpolation should be evaluated
        :returns: The interpolated value
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef double _r = self.inst.get().evaluate((<double>value))
        py_result = <double>_r
        return py_result
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def weightData(self, list data ):
        """
        weightData(self, data: List[TM_DataPoint] ) -> None
        Weight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().weightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkValidWeight(self,  weight , list valid_weights ):
        """
        checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool
        Check for a valid weighting function string
        """
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
        assert isinstance(valid_weights, list) and all(isinstance(i, bytes) for i in valid_weights), 'arg valid_weights wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in valid_weights:
           v1.push_back(_String(<char *>item1))
        cdef bool _r = self.inst.get().checkValidWeight(deref((convString(weight)).get()), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        valid_weights[:] = replace
        del v1
        py_result = <bool>_r
        return py_result
    
    def weightDatum(self, double datum ,  weight ):
        """
        weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Weight the data according to the weighting function
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().weightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def unWeightDatum(self, double datum ,  weight ):
        """
        unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Apply the reverse of the weighting function to the data
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().unWeightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def getValidXWeights(self):
        """
        getValidXWeights(self) -> List[bytes]
        Returns a list of valid x weight function stringss
        """
        _r = self.inst.get().getValidXWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getValidYWeights(self):
        """
        getValidYWeights(self) -> List[bytes]
        Returns a list of valid y weight function strings
        """
        _r = self.inst.get().getValidYWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def unWeightData(self, list data ):
        """
        unWeightData(self, data: List[TM_DataPoint] ) -> None
        Unweight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().unWeightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkDatumRange(self, double datum , double datum_min , double datum_max ):
        """
        checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float
        Check that the datum is within the valid min and max bounds
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert isinstance(datum_min, float), 'arg datum_min wrong type'
        assert isinstance(datum_max, float), 'arg datum_max wrong type'
    
    
    
        cdef double _r = self.inst.get().checkDatumRange((<const double &>datum), (<const double &>datum_min), (<const double &>datum_max))
        py_result = <double>_r
        return py_result 

cdef class TransformationModelLowess:
    """
    Cython implementation of _TransformationModelLowess

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationModelLowess.html>`_
      -- Inherits from ['TransformationModel']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, list data , Param params ):
        """
        __init__(self, data: List[TM_DataPoint] , params: Param ) -> None
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        assert isinstance(params, Param), 'arg params wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst = shared_ptr[_TransformationModelLowess](new _TransformationModelLowess(deref(v0), (deref(params.inst.get()))))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def evaluate(self, double value ):
        """
        evaluate(self, value: float ) -> float
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        cdef double _r = self.inst.get().evaluate((<double>value))
        py_result = <double>_r
        return py_result
    
    def getParameters(self):
        """
        getParameters(self) -> Param
        """
        cdef _Param * _r = new _Param(self.inst.get().getParameters())
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def weightData(self, list data ):
        """
        weightData(self, data: List[TM_DataPoint] ) -> None
        Weight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().weightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkValidWeight(self,  weight , list valid_weights ):
        """
        checkValidWeight(self, weight: Union[bytes, str, String] , valid_weights: List[bytes] ) -> bool
        Check for a valid weighting function string
        """
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
        assert isinstance(valid_weights, list) and all(isinstance(i, bytes) for i in valid_weights), 'arg valid_weights wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in valid_weights:
           v1.push_back(_String(<char *>item1))
        cdef bool _r = self.inst.get().checkValidWeight(deref((convString(weight)).get()), deref(v1))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        valid_weights[:] = replace
        del v1
        py_result = <bool>_r
        return py_result
    
    def weightDatum(self, double datum ,  weight ):
        """
        weightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Weight the data according to the weighting function
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().weightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def unWeightDatum(self, double datum ,  weight ):
        """
        unWeightDatum(self, datum: float , weight: Union[bytes, str, String] ) -> float
        Apply the reverse of the weighting function to the data
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert (isinstance(weight, str) or isinstance(weight, bytes) or isinstance(weight, String)), 'arg weight wrong type'
    
    
        cdef double _r = self.inst.get().unWeightDatum((<double &>datum), deref((convString(weight)).get()))
        py_result = <double>_r
        return py_result
    
    def getValidXWeights(self):
        """
        getValidXWeights(self) -> List[bytes]
        Returns a list of valid x weight function stringss
        """
        _r = self.inst.get().getValidXWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getValidYWeights(self):
        """
        getValidYWeights(self) -> List[bytes]
        Returns a list of valid y weight function strings
        """
        _r = self.inst.get().getValidYWeights()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def unWeightData(self, list data ):
        """
        unWeightData(self, data: List[TM_DataPoint] ) -> None
        Unweight the data by the given weight function
        """
        assert isinstance(data, list) and all(isinstance(elemt_rec, TM_DataPoint) for elemt_rec in data), 'arg data wrong type'
        cdef libcpp_vector[_TM_DataPoint] * v0 = new libcpp_vector[_TM_DataPoint]()
        cdef TM_DataPoint item0
        for item0 in data:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().unWeightData(deref(v0))
        cdef libcpp_vector[_TM_DataPoint].iterator it_data = v0.begin()
        replace_0 = []
        while it_data != v0.end():
            item0 = TM_DataPoint.__new__(TM_DataPoint)
            item0.inst = shared_ptr[_TM_DataPoint](new _TM_DataPoint(deref(it_data)))
            replace_0.append(item0)
            inc(it_data)
        data[:] = replace_0
        del v0
    
    def checkDatumRange(self, double datum , double datum_min , double datum_max ):
        """
        checkDatumRange(self, datum: float , datum_min: float , datum_max: float ) -> float
        Check that the datum is within the valid min and max bounds
        """
        assert isinstance(datum, float), 'arg datum wrong type'
        assert isinstance(datum_min, float), 'arg datum_min wrong type'
        assert isinstance(datum_max, float), 'arg datum_max wrong type'
    
    
    
        cdef double _r = self.inst.get().checkDatumRange((<const double &>datum), (<const double &>datum_min), (<const double &>datum_max))
        py_result = <double>_r
        return py_result
    getDefaultParameters = __static_TransformationModelLowess_getDefaultParameters 
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
