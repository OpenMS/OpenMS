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

cdef class __IntensityThresholdCalculation:
    None
    MANUAL = 0
    AUTOMAXBYSTDEV = 1
    AUTOMAXBYPERCENT = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ResidueType:
    None
    Full = 0
    Internal = 1
    NTerminal = 2
    CTerminal = 3
    AIon = 4
    BIon = 5
    CIon = 6
    XIon = 7
    YIon = 8
    ZIon = 9
    Precursor_ion = 10
    BIonMinusH20 = 11
    YIonMinusH20 = 12
    BIonMinusNH3 = 13
    YIonMinusNH3 = 14
    NonIdentified = 15
    Unannotated = 16
    SizeOfResidueType = 17

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AccurateMassSearchEngine:
    """
    Cython implementation of _AccurateMassSearchEngine

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AccurateMassSearchEngine.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AccurateMassSearchEngine rv = AccurateMassSearchEngine.__new__(AccurateMassSearchEngine)
       rv.inst = shared_ptr[_AccurateMassSearchEngine](new _AccurateMassSearchEngine(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AccurateMassSearchEngine rv = AccurateMassSearchEngine.__new__(AccurateMassSearchEngine)
       rv.inst = shared_ptr[_AccurateMassSearchEngine](new _AccurateMassSearchEngine(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AccurateMassSearchEngine](new _AccurateMassSearchEngine())
    
    def _init_1(self, AccurateMassSearchEngine in_0 ):
        """
        _init_1(self, in_0: AccurateMassSearchEngine ) -> None
        """
        assert isinstance(in_0, AccurateMassSearchEngine), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AccurateMassSearchEngine](new _AccurateMassSearchEngine((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AccurateMassSearchEngine ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AccurateMassSearchEngine)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def queryByMZ(self, double observed_mz ,  observed_charge ,  ion_mode , list in_3 , EmpiricalFormula observed_adduct ):
        """
        queryByMZ(self, observed_mz: float , observed_charge: int , ion_mode: Union[bytes, str, String] , in_3: List[AccurateMassSearchResult] , observed_adduct: EmpiricalFormula ) -> None
        """
        assert isinstance(observed_mz, float), 'arg observed_mz wrong type'
        assert isinstance(observed_charge, int), 'arg observed_charge wrong type'
        assert (isinstance(ion_mode, str) or isinstance(ion_mode, bytes) or isinstance(ion_mode, String)), 'arg ion_mode wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, AccurateMassSearchResult) for elemt_rec in in_3), 'arg in_3 wrong type'
        assert isinstance(observed_adduct, EmpiricalFormula), 'arg observed_adduct wrong type'
    
    
    
        cdef libcpp_vector[_AccurateMassSearchResult] * v3 = new libcpp_vector[_AccurateMassSearchResult]()
        cdef AccurateMassSearchResult item3
        for item3 in in_3:
            v3.push_back(deref(item3.inst.get()))
    
        self.inst.get().queryByMZ((<double>observed_mz), (<int>observed_charge), deref((convString(ion_mode)).get()), deref(v3), (deref(observed_adduct.inst.get())))
        cdef libcpp_vector[_AccurateMassSearchResult].iterator it_in_3 = v3.begin()
        replace_0 = []
        while it_in_3 != v3.end():
            item3 = AccurateMassSearchResult.__new__(AccurateMassSearchResult)
            item3.inst = shared_ptr[_AccurateMassSearchResult](new _AccurateMassSearchResult(deref(it_in_3)))
            replace_0.append(item3)
            inc(it_in_3)
        in_3[:] = replace_0
        del v3
    
    def queryByFeature(self, Feature feature ,  feature_index ,  ion_mode , list in_3 ):
        """
        queryByFeature(self, feature: Feature , feature_index: int , ion_mode: Union[bytes, str, String] , in_3: List[AccurateMassSearchResult] ) -> None
        """
        assert isinstance(feature, Feature), 'arg feature wrong type'
        assert isinstance(feature_index, int) and feature_index >= 0, 'arg feature_index wrong type'
        assert (isinstance(ion_mode, str) or isinstance(ion_mode, bytes) or isinstance(ion_mode, String)), 'arg ion_mode wrong type'
        assert isinstance(in_3, list) and all(isinstance(elemt_rec, AccurateMassSearchResult) for elemt_rec in in_3), 'arg in_3 wrong type'
    
    
    
        cdef libcpp_vector[_AccurateMassSearchResult] * v3 = new libcpp_vector[_AccurateMassSearchResult]()
        cdef AccurateMassSearchResult item3
        for item3 in in_3:
            v3.push_back(deref(item3.inst.get()))
        self.inst.get().queryByFeature((deref(feature.inst.get())), (<size_t>feature_index), deref((convString(ion_mode)).get()), deref(v3))
        cdef libcpp_vector[_AccurateMassSearchResult].iterator it_in_3 = v3.begin()
        replace_0 = []
        while it_in_3 != v3.end():
            item3 = AccurateMassSearchResult.__new__(AccurateMassSearchResult)
            item3.inst = shared_ptr[_AccurateMassSearchResult](new _AccurateMassSearchResult(deref(it_in_3)))
            replace_0.append(item3)
            inc(it_in_3)
        in_3[:] = replace_0
        del v3
    
    def queryByConsensusFeature(self, ConsensusFeature cfeat ,  cf_index ,  number_of_maps ,  ion_mode , list results ):
        """
        queryByConsensusFeature(self, cfeat: ConsensusFeature , cf_index: int , number_of_maps: int , ion_mode: Union[bytes, str, String] , results: List[AccurateMassSearchResult] ) -> None
        """
        assert isinstance(cfeat, ConsensusFeature), 'arg cfeat wrong type'
        assert isinstance(cf_index, int) and cf_index >= 0, 'arg cf_index wrong type'
        assert isinstance(number_of_maps, int) and number_of_maps >= 0, 'arg number_of_maps wrong type'
        assert (isinstance(ion_mode, str) or isinstance(ion_mode, bytes) or isinstance(ion_mode, String)), 'arg ion_mode wrong type'
        assert isinstance(results, list) and all(isinstance(elemt_rec, AccurateMassSearchResult) for elemt_rec in results), 'arg results wrong type'
    
    
    
    
        cdef libcpp_vector[_AccurateMassSearchResult] * v4 = new libcpp_vector[_AccurateMassSearchResult]()
        cdef AccurateMassSearchResult item4
        for item4 in results:
            v4.push_back(deref(item4.inst.get()))
        self.inst.get().queryByConsensusFeature((deref(cfeat.inst.get())), (<size_t>cf_index), (<size_t>number_of_maps), deref((convString(ion_mode)).get()), deref(v4))
        cdef libcpp_vector[_AccurateMassSearchResult].iterator it_results = v4.begin()
        replace_0 = []
        while it_results != v4.end():
            item4 = AccurateMassSearchResult.__new__(AccurateMassSearchResult)
            item4.inst = shared_ptr[_AccurateMassSearchResult](new _AccurateMassSearchResult(deref(it_results)))
            replace_0.append(item4)
            inc(it_results)
        results[:] = replace_0
        del v4
    
    def _run_0(self, FeatureMap in_0 , MzTab in_1 ):
        """
        _run_0(self, in_0: FeatureMap , in_1: MzTab ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, MzTab), 'arg in_1 wrong type'
    
    
        self.inst.get().run((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def _run_1(self, FeatureMap in_0 , MzTabM in_1 ):
        """
        _run_1(self, in_0: FeatureMap , in_1: MzTabM ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, MzTabM), 'arg in_1 wrong type'
    
    
        self.inst.get().run((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def _run_2(self, ConsensusMap in_0 , MzTab in_1 ):
        """
        _run_2(self, in_0: ConsensusMap , in_1: MzTab ) -> None
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
        assert isinstance(in_1, MzTab), 'arg in_1 wrong type'
    
    
        self.inst.get().run((deref(in_0.inst.get())), (deref(in_1.inst.get())))
    
    def run(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: run(self, in_0: FeatureMap , in_1: MzTab ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: run(self, in_0: FeatureMap , in_1: MzTabM ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: run(self, in_0: ConsensusMap , in_1: MzTab ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], MzTab)):
            return self._run_0(*args)
        elif (len(args)==2) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], MzTabM)):
            return self._run_1(*args)
        elif (len(args)==2) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], MzTab)):
            return self._run_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def init(self):
        """
        init(self) -> None
        """
        self.inst.get().init()
    
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

cdef class CVMappingFile:
    """
    Cython implementation of _CVMappingFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_CVMappingFile](new _CVMappingFile())
    
    def load(self,  filename , CVMappings cv_mappings , bool strip_namespaces ):
        """
        load(self, filename: Union[bytes, str, String] , cv_mappings: CVMappings , strip_namespaces: bool ) -> None
        Loads CvMappings from the given file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(cv_mappings, CVMappings), 'arg cv_mappings wrong type'
        assert isinstance(strip_namespaces, pybool_t), 'arg strip_namespaces wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(cv_mappings.inst.get())), (<bool>strip_namespaces)) 

cdef class FASTAEntry:
    """
    Cython implementation of _FASTAEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FASTAEntry.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property identifier:
        def __set__(self,  identifier):
        
            self.inst.get().identifier = deref((convString(identifier)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().identifier
            py_result = convOutputString(_r)
            return py_result
    
    property description:
        def __set__(self,  description):
        
            self.inst.get().description = deref((convString(description)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().description
            py_result = convOutputString(_r)
            return py_result
    
    property sequence:
        def __set__(self,  sequence):
        
            self.inst.get().sequence = deref((convString(sequence)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().sequence
            py_result = convOutputString(_r)
            return py_result
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FASTAEntry](new _FASTAEntry())
    
    def _init_1(self, FASTAEntry in_0 ):
        """
        _init_1(self, in_0: FASTAEntry ) -> None
        """
        assert isinstance(in_0, FASTAEntry), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FASTAEntry](new _FASTAEntry((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FASTAEntry ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FASTAEntry)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def headerMatches(self, FASTAEntry rhs ):
        """
        headerMatches(self, rhs: FASTAEntry ) -> bool
        """
        assert isinstance(rhs, FASTAEntry), 'arg rhs wrong type'
    
        cdef bool _r = self.inst.get().headerMatches((deref(rhs.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def sequenceMatches(self, FASTAEntry rhs ):
        """
        sequenceMatches(self, rhs: FASTAEntry ) -> bool
        """
        assert isinstance(rhs, FASTAEntry), 'arg rhs wrong type'
    
        cdef bool _r = self.inst.get().sequenceMatches((deref(rhs.inst.get())))
        py_result = <bool>_r
        return py_result 

cdef class FASTAFile:
    """
    Cython implementation of _FASTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FASTAFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        This class serves for reading in and writing FASTA files
        """
        self.inst = shared_ptr[_FASTAFile](new _FASTAFile())
    
    def load(self,  filename , list data ):
        """
        load(self, filename: Union[bytes, str, String] , data: List[FASTAEntry] ) -> None
        Loads a FASTA file given by 'filename' and stores the information in 'data'
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(data, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in data), 'arg data wrong type'
    
        cdef libcpp_vector[_FASTAEntry] * v1 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item1
        for item1 in data:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().load(deref((convString(filename)).get()), deref(v1))
        cdef libcpp_vector[_FASTAEntry].iterator it_data = v1.begin()
        replace_0 = []
        while it_data != v1.end():
            item1 = FASTAEntry.__new__(FASTAEntry)
            item1.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_data)))
            replace_0.append(item1)
            inc(it_data)
        data[:] = replace_0
        del v1
    
    def store(self,  filename , list data ):
        """
        store(self, filename: Union[bytes, str, String] , data: List[FASTAEntry] ) -> None
        Stores the data given by 'data' at the file 'filename'
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(data, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in data), 'arg data wrong type'
    
        cdef libcpp_vector[_FASTAEntry] * v1 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item1
        for item1 in data:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().store(deref((convString(filename)).get()), deref(v1))
        cdef libcpp_vector[_FASTAEntry].iterator it_data = v1.begin()
        replace_0 = []
        while it_data != v1.end():
            item1 = FASTAEntry.__new__(FASTAEntry)
            item1.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_data)))
            replace_0.append(item1)
            inc(it_data)
        data[:] = replace_0
        del v1
    
    def readStart(self,  filename ):
        """
        readStart(self, filename: Union[bytes, str, String] ) -> None
        Prepares a FASTA file given by 'filename' for streamed reading using readNext()
        
        :raises:
            Exception:FileNotFound is thrown if the file does not exists
        :raises:
            Exception:ParseError is thrown if the file does not suit to the standard
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().readStart(deref((convString(filename)).get()))
    
    def readNext(self, FASTAEntry protein ):
        """
        readNext(self, protein: FASTAEntry ) -> bool
        Reads the next FASTA entry from file
        
        If you want to read all entries in one go, use load()
        
        :return: true if entry was read; false if eof was reached
        :raises:
            Exception:FileNotFound is thrown if the file does not exists
        :raises:
            Exception:ParseError is thrown if the file does not suit to the standard
        """
        assert isinstance(protein, FASTAEntry), 'arg protein wrong type'
    
        cdef bool _r = self.inst.get().readNext((deref(protein.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def atEnd(self):
        """
        atEnd(self) -> bool
        Boolean function to check if streams is at end of file
        """
        cdef bool _r = self.inst.get().atEnd()
        py_result = <bool>_r
        return py_result
    
    def writeStart(self,  filename ):
        """
        writeStart(self, filename: Union[bytes, str, String] ) -> None
        Prepares a FASTA file given by 'filename' for streamed writing using writeNext()
        
        :raises:
            Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().writeStart(deref((convString(filename)).get()))
    
    def writeNext(self, FASTAEntry protein ):
        """
        writeNext(self, protein: FASTAEntry ) -> None
        Stores the data given by `protein`. Call writeStart() once before calling writeNext()
        
        Call writeEnd() when done to close the file!
        
        :raises:
            Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        """
        assert isinstance(protein, FASTAEntry), 'arg protein wrong type'
    
        self.inst.get().writeNext((deref(protein.inst.get())))
    
    def writeEnd(self):
        """
        writeEnd(self) -> None
        Closes the file (flush). Called implicitly when FASTAFile object does out of scope
        """
        self.inst.get().writeEnd() 

cdef class IMSAlphabet:
    """
    Cython implementation of _IMSAlphabet

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSAlphabet_1_1IMSAlphabet.html>`_

    Holds an indexed list of bio-chemical elements.\n
    
    Presents an indexed list of bio-chemical elements of type (or derived from
    type) 'Element'. Due to indexed structure 'Alphabet' can be used similar
    to std::vector, for example to add a new element to 'Alphabet' function
    push_back(element_type) can be used. Elements or their properties (such
    as element's mass) can be accessed by index in a constant time. On the other
    hand accessing elements by their names takes linear time. Due to this and
    also the fact that 'Alphabet' is 'heavy-weighted' (consisting of
    'Element' -s or their derivatives where the depth of derivation as well is
    undefined resulting in possibly 'heavy' access operations) it is recommended
    not use 'Alphabet' directly in operations where fast access to
    'Element' 's properties is required. Instead consider to use
    'light-weighted' equivalents, such as 'Weights'
    
    
    :param map: MSExperiment to receive the identifications
    :param fmap: FeatureMap with PeptideIdentifications for the MSExperiment
    :param clear_ids: Reset peptide and protein identifications of each scan before annotating
    :param map_ms1: Attach Ids to MS1 spectra using RT mapping only (without precursor, without m/z)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IMSAlphabet rv = IMSAlphabet.__new__(IMSAlphabet)
       rv.inst = shared_ptr[_IMSAlphabet](new _IMSAlphabet(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IMSAlphabet rv = IMSAlphabet.__new__(IMSAlphabet)
       rv.inst = shared_ptr[_IMSAlphabet](new _IMSAlphabet(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IMSAlphabet](new _IMSAlphabet())
    
    def _init_1(self, IMSAlphabet in_0 ):
        """
        _init_1(self, in_0: IMSAlphabet ) -> None
        """
        assert isinstance(in_0, IMSAlphabet), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IMSAlphabet](new _IMSAlphabet((deref(in_0.inst.get()))))
    
    def _init_2(self, list elements ):
        """
        _init_2(self, elements: List[IMSElement] ) -> None
        """
        assert isinstance(elements, list) and all(isinstance(elemt_rec, IMSElement) for elemt_rec in elements), 'arg elements wrong type'
        cdef libcpp_vector[_IMSElement] * v0 = new libcpp_vector[_IMSElement]()
        cdef IMSElement item0
        for item0 in elements:
            v0.push_back(deref(item0.inst.get()))
        self.inst = shared_ptr[_IMSAlphabet](new _IMSAlphabet(deref(v0)))
        cdef libcpp_vector[_IMSElement].iterator it_elements = v0.begin()
        replace_0 = []
        while it_elements != v0.end():
            item0 = IMSElement.__new__(IMSElement)
            item0.inst = shared_ptr[_IMSElement](new _IMSElement(deref(it_elements)))
            replace_0.append(item0)
            inc(it_elements)
        elements[:] = replace_0
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IMSAlphabet ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, elements: List[IMSElement] ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSAlphabet)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, IMSElement) for elemt_rec in args[0])):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getElement_0(self, bytes name ):
        """
        _getElement_0(self, name: bytes ) -> IMSElement
        Gets the element with 'index' and returns element with the given index in alphabet
        """
        assert isinstance(name, bytes), 'arg name wrong type'
    
        cdef _IMSElement * _r = new _IMSElement(self.inst.get().getElement((<libcpp_string>name)))
        cdef IMSElement py_result = IMSElement.__new__(IMSElement)
        py_result.inst = shared_ptr[_IMSElement](_r)
        return py_result
    
    def _getElement_1(self,  index ):
        """
        _getElement_1(self, index: int ) -> IMSElement
        Gets the element with 'index'
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef _IMSElement * _r = new _IMSElement(self.inst.get().getElement((<int>index)))
        cdef IMSElement py_result = IMSElement.__new__(IMSElement)
        py_result.inst = shared_ptr[_IMSElement](_r)
        return py_result
    
    def getElement(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getElement(self, name: bytes ) -> IMSElement
          :noindex:
        
        Gets the element with 'index' and returns element with the given index in alphabet

        
        .. rubric:: Overload:
        .. py:function:: getElement(self, index: int ) -> IMSElement
          :noindex:
        
        Gets the element with 'index'
    
        """
        if (len(args)==1) and (isinstance(args[0], bytes)):
            return self._getElement_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._getElement_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getName(self,  index ):
        """
        getName(self, index: int ) -> bytes
        Gets the symbol of the element with an 'index' in alphabet
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef libcpp_string _r = self.inst.get().getName((<int>index))
        py_result = <libcpp_string>_r
        return py_result
    
    def _getMass_0(self, bytes name ):
        """
        _getMass_0(self, name: bytes ) -> float
        Gets mono isotopic mass of the element with the symbol 'name'
        """
        assert isinstance(name, bytes), 'arg name wrong type'
    
        cdef double _r = self.inst.get().getMass((<libcpp_string>name))
        py_result = <double>_r
        return py_result
    
    def _getMass_1(self,  index ):
        """
        _getMass_1(self, index: int ) -> float
        Gets mass of the element with an 'index' in alphabet
        """
        assert isinstance(index, int), 'arg index wrong type'
    
        cdef double _r = self.inst.get().getMass((<int>index))
        py_result = <double>_r
        return py_result
    
    def getMass(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMass(self, name: bytes ) -> float
          :noindex:
        
        Gets mono isotopic mass of the element with the symbol 'name'
        
        .. rubric:: Overload:
        .. py:function:: getMass(self, index: int ) -> float
          :noindex:
        
        Gets mass of the element with an 'index' in alphabet
    
        """
        if (len(args)==1) and (isinstance(args[0], bytes)):
            return self._getMass_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int)):
            return self._getMass_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMasses(self,  isotope_index ):
        """
        getMasses(self, isotope_index: int ) -> List[float]
        Gets masses of elements isotopes given by 'isotope_index'
        """
        assert isinstance(isotope_index, int), 'arg isotope_index wrong type'
    
        _r = self.inst.get().getMasses((<int>isotope_index))
        cdef list py_result = _r
        return py_result
    
    def getAverageMasses(self):
        """
        getAverageMasses(self) -> List[float]
        Gets average masses of elements
        """
        _r = self.inst.get().getAverageMasses()
        cdef list py_result = _r
        return py_result
    
    def hasName(self, bytes name ):
        """
        hasName(self, name: bytes ) -> bool
        Returns true if there is an element with symbol 'name' in the alphabet, false - otherwise
        """
        assert isinstance(name, bytes), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasName((<libcpp_string>name))
        py_result = <bool>_r
        return py_result
    
    def _push_back_0(self, bytes name , double value ):
        """
        _push_back_0(self, name: bytes , value: float ) -> None
        Adds a new element with 'name' and mass 'value'
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
        self.inst.get().push_back((<libcpp_string>name), (<double>value))
    
    def _push_back_1(self, IMSElement element ):
        """
        _push_back_1(self, element: IMSElement ) -> None
        Adds a new 'element' to the alphabet
        """
        assert isinstance(element, IMSElement), 'arg element wrong type'
    
        self.inst.get().push_back((deref(element.inst.get())))
    
    def push_back(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: push_back(self, name: bytes , value: float ) -> None
          :noindex:
        
        Adds a new element with 'name' and mass 'value'

        
        .. rubric:: Overload:
        .. py:function:: push_back(self, element: IMSElement ) -> None
          :noindex:
        
        Adds a new 'element' to the alphabet
    
        """
        if (len(args)==2) and (isinstance(args[0], bytes)) and (isinstance(args[1], float)):
            return self._push_back_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IMSElement)):
            return self._push_back_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def clear(self):
        """
        clear(self) -> None
        Clears the alphabet data
        """
        self.inst.get().clear()
    
    def sortByNames(self):
        """
        sortByNames(self) -> None
        Sorts the alphabet by names
        """
        self.inst.get().sortByNames()
    
    def sortByValues(self):
        """
        sortByValues(self) -> None
        Sorts the alphabet by mass values
        """
        self.inst.get().sortByValues()
    
    def load(self,  fname ):
        """
        load(self, fname: String ) -> None
        Loads the alphabet data from the file 'fname' using the default parser. If there is no file 'fname', throws an 'IOException'
        """
        assert isinstance(fname, String), 'arg fname wrong type'
    
        self.inst.get().load(deref((<String>fname).inst.get()))
    
    def size(self):
        """
        size(self) -> int
        """
        cdef int _r = self.inst.get().size()
        py_result = <int>_r
        return py_result
    
    def setElement(self, bytes name , double mass , bool forced ):
        """
        setElement(self, name: bytes , mass: float , forced: bool ) -> None
        Overwrites an element in the alphabet with the 'name' with a new element constructed from the given 'name' and 'mass'
        """
        assert isinstance(name, bytes), 'arg name wrong type'
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(forced, pybool_t), 'arg forced wrong type'
    
    
    
        self.inst.get().setElement((<libcpp_string>name), (<double>mass), (<bool>forced))
    
    def erase(self, bytes name ):
        """
        erase(self, name: bytes ) -> bool
        Removes the element with 'name' from the alphabet
        """
        assert isinstance(name, bytes), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().erase((<libcpp_string>name))
        py_result = <bool>_r
        return py_result 

cdef class MapAlignmentAlgorithmIdentification:
    """
    Cython implementation of _MapAlignmentAlgorithmIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmIdentification.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_MapAlignmentAlgorithmIdentification](new _MapAlignmentAlgorithmIdentification())
    
    def _align_0(self, list in_0 , list in_1 ,  in_2 ):
        """
        _align_0(self, in_0: List[FeatureMap] , in_1: List[TransformationDescription] , in_2: int ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, TransformationDescription) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(in_2, int), 'arg in_2 wrong type'
        cdef libcpp_vector[_FeatureMap] * v0 = new libcpp_vector[_FeatureMap]()
        cdef FeatureMap item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_TransformationDescription] * v1 = new libcpp_vector[_TransformationDescription]()
        cdef TransformationDescription item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().align(deref(v0), deref(v1), (<int>in_2))
        cdef libcpp_vector[_TransformationDescription].iterator it_in_1 = v1.begin()
        replace_0 = []
        while it_in_1 != v1.end():
            item1 = TransformationDescription.__new__(TransformationDescription)
            item1.inst = shared_ptr[_TransformationDescription](new _TransformationDescription(deref(it_in_1)))
            replace_0.append(item1)
            inc(it_in_1)
        in_1[:] = replace_0
        del v1
        cdef libcpp_vector[_FeatureMap].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = FeatureMap.__new__(FeatureMap)
            item0.inst = shared_ptr[_FeatureMap](new _FeatureMap(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def _align_1(self, list in_0 , list in_1 ,  in_2 ):
        """
        _align_1(self, in_0: List[ConsensusMap] , in_1: List[TransformationDescription] , in_2: int ) -> None
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in in_0), 'arg in_0 wrong type'
        assert isinstance(in_1, list) and all(isinstance(elemt_rec, TransformationDescription) for elemt_rec in in_1), 'arg in_1 wrong type'
        assert isinstance(in_2, int), 'arg in_2 wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_TransformationDescription] * v1 = new libcpp_vector[_TransformationDescription]()
        cdef TransformationDescription item1
        for item1 in in_1:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().align(deref(v0), deref(v1), (<int>in_2))
        cdef libcpp_vector[_TransformationDescription].iterator it_in_1 = v1.begin()
        replace_0 = []
        while it_in_1 != v1.end():
            item1 = TransformationDescription.__new__(TransformationDescription)
            item1.inst = shared_ptr[_TransformationDescription](new _TransformationDescription(deref(it_in_1)))
            replace_0.append(item1)
            inc(it_in_1)
        in_1[:] = replace_0
        del v1
        cdef libcpp_vector[_ConsensusMap].iterator it_in_0 = v0.begin()
        replace_0 = []
        while it_in_0 != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_in_0)))
            replace_0.append(item0)
            inc(it_in_0)
        in_0[:] = replace_0
        del v0
    
    def align(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: align(self, in_0: List[FeatureMap] , in_1: List[TransformationDescription] , in_2: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: align(self, in_0: List[ConsensusMap] , in_1: List[TransformationDescription] , in_2: int ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, FeatureMap) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, TransformationDescription) for elemt_rec in args[1])) and (isinstance(args[2], int)):
            return self._align_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, TransformationDescription) for elemt_rec in args[1])) and (isinstance(args[2], int)):
            return self._align_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setReference_0(self, FeatureMap in_0 ):
        """
        _setReference_0(self, in_0: FeatureMap ) -> None
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
    
        self.inst.get().setReference((deref(in_0.inst.get())))
    
    def _setReference_1(self, ConsensusMap in_0 ):
        """
        _setReference_1(self, in_0: ConsensusMap ) -> None
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
    
        self.inst.get().setReference((deref(in_0.inst.get())))
    
    def _setReference_2(self, PeptideIdentificationList in_0 ):
        """
        _setReference_2(self, in_0: PeptideIdentificationList ) -> None
        """
        assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    
        self.inst.get().setReference((deref(in_0.inst.get())))
    
    def setReference(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setReference(self, in_0: FeatureMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setReference(self, in_0: ConsensusMap ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setReference(self, in_0: PeptideIdentificationList ) -> None
          :noindex:
    
        """
        if (len(args)==1) and (isinstance(args[0], FeatureMap)):
            return self._setReference_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ConsensusMap)):
            return self._setReference_1(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._setReference_2(*args)
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

cdef class OSSpectrumMeta:
    """
    Cython implementation of _OSSpectrumMeta

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSSpectrumMeta.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property index:
        def __set__(self,  index):
        
            self.inst.get().index = (<size_t>index)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().index
            py_result = <size_t>_r
            return py_result
    
    property id:
        def __set__(self, bytes id):
        
            self.inst.get().id = (<libcpp_string>id)
        
    
        def __get__(self):
            cdef libcpp_string _r = self.inst.get().id
            py_result = <libcpp_string>_r
            return py_result
    
    property RT:
        def __set__(self, double RT):
        
            self.inst.get().RT = (<double>RT)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().RT
            py_result = <double>_r
            return py_result
    
    property ms_level:
        def __set__(self,  ms_level):
        
            self.inst.get().ms_level = (<int>ms_level)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().ms_level
            py_result = <int>_r
            return py_result
    
    def __copy__(self):
       cdef OSSpectrumMeta rv = OSSpectrumMeta.__new__(OSSpectrumMeta)
       rv.inst = shared_ptr[_OSSpectrumMeta](new _OSSpectrumMeta(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSSpectrumMeta rv = OSSpectrumMeta.__new__(OSSpectrumMeta)
       rv.inst = shared_ptr[_OSSpectrumMeta](new _OSSpectrumMeta(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OSSpectrumMeta](new _OSSpectrumMeta())
    
    def _init_1(self, OSSpectrumMeta in_0 ):
        """
        _init_1(self, in_0: OSSpectrumMeta ) -> None
        """
        assert isinstance(in_0, OSSpectrumMeta), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSSpectrumMeta](new _OSSpectrumMeta((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSSpectrumMeta ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSSpectrumMeta)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PScore:
    """
    Cython implementation of _PScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PScore.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PScore rv = PScore.__new__(PScore)
       rv.inst = shared_ptr[_PScore](new _PScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PScore rv = PScore.__new__(PScore)
       rv.inst = shared_ptr[_PScore](new _PScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PScore](new _PScore())
    
    def _init_1(self, PScore in_0 ):
        """
        _init_1(self, in_0: PScore ) -> None
        """
        assert isinstance(in_0, PScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PScore](new _PScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PScore)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def calculateIntensityRankInMZWindow(self, list mz , list intensities , double mz_window ):
        """
        calculateIntensityRankInMZWindow(self, mz: List[float] , intensities: List[float] , mz_window: float ) -> List[int]
        Calculate local (windowed) peak ranks
        
        The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
        The result can be used to efficiently filter spectra for top 1..n peaks in mass windows
        
        
        :param mz: The m/z positions of the peaks
        :param intensities: The intensities of the peaks
        :param mz_window: The window in Thomson centered at each peak
        """
        assert isinstance(mz, list) and all(isinstance(elemt_rec, float) for elemt_rec in mz), 'arg mz wrong type'
        assert isinstance(intensities, list) and all(isinstance(elemt_rec, float) for elemt_rec in intensities), 'arg intensities wrong type'
        assert isinstance(mz_window, float), 'arg mz_window wrong type'
        cdef libcpp_vector[double] v0 = mz
        cdef libcpp_vector[double] v1 = intensities
    
        _r = self.inst.get().calculateIntensityRankInMZWindow(v0, v1, (<double>mz_window))
        intensities[:] = v1
        mz[:] = v0
        cdef list py_result = _r
        return py_result
    
    def calculateRankMap(self, MSExperiment peak_map , double mz_window ):
        """
        calculateRankMap(self, peak_map: MSExperiment , mz_window: float ) -> List[List[int]]
        Precalculated, windowed peak ranks for a whole experiment
        
        The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
        
        
        :param peak_map: Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra
        :param mz_window: Window in Thomson centered at each peak
        """
        assert isinstance(peak_map, MSExperiment), 'arg peak_map wrong type'
        assert isinstance(mz_window, float), 'arg mz_window wrong type'
    
    
        _r = self.inst.get().calculateRankMap((deref(peak_map.inst.get())), (<double>mz_window))
        cdef list py_result = _r
        return py_result
    
    def calculatePeakLevelSpectra(self, MSSpectrum spec , list ranks ,  min_level ,  max_level ):
        """
        calculatePeakLevelSpectra(self, spec: MSSpectrum , ranks: List[int] , min_level: int , max_level: int ) -> Dict[int, MSSpectrum]
        Calculates spectra for peak level between min_level to max_level and stores them in the map
        
        A spectrum of peak level n retains the (n+1) top intensity peaks in a sliding mz_window centered at each peak
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
        assert isinstance(ranks, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in ranks), 'arg ranks wrong type'
        assert isinstance(min_level, int) and min_level >= 0, 'arg min_level wrong type'
        assert isinstance(max_level, int) and max_level >= 0, 'arg max_level wrong type'
    
        cdef libcpp_vector[size_t] v1 = ranks
    
    
        _r = self.inst.get().calculatePeakLevelSpectra((deref(spec.inst.get())), v1, (<size_t>min_level), (<size_t>max_level))
        ranks[:] = v1
        py_result = dict()
        cdef libcpp_map[size_t, _MSSpectrum].iterator it__r = _r.begin()
        cdef MSSpectrum item_py_result
        while it__r != _r.end():
           item_py_result = MSSpectrum.__new__(MSSpectrum)
           item_py_result.inst = shared_ptr[_MSSpectrum](new _MSSpectrum((deref(it__r)).second))
           py_result[<size_t>(deref(it__r).first)] = item_py_result
           inc(it__r)
        return py_result
    
    def _computePScore_0(self, double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , dict peak_level_spectra , list theo_spectra , double mz_window ):
        """
        _computePScore_0(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectra: List[MSSpectrum] , mz_window: float ) -> float
        Computes the PScore for a vector of theoretical spectra
        
        Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
        The best score obtained by scoring all those theoretical spectra against the experimental ones is returned
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank).
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
        """
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(peak_level_spectra, dict) and all(isinstance(k, int) and k >= 0 for k in peak_level_spectra.keys()) and all(isinstance(v, MSSpectrum) for v in peak_level_spectra.values()), 'arg peak_level_spectra wrong type'
        assert isinstance(theo_spectra, list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in theo_spectra), 'arg theo_spectra wrong type'
        assert isinstance(mz_window, float), 'arg mz_window wrong type'
    
    
        cdef libcpp_map[size_t, _MSSpectrum] * v2 = new libcpp_map[size_t, _MSSpectrum]()
        for key, value in peak_level_spectra.items():
        
        
            deref(v2)[ (<size_t>key) ] = deref((<MSSpectrum>value).inst.get())
        
        
        cdef libcpp_vector[_MSSpectrum] * v3 = new libcpp_vector[_MSSpectrum]()
        cdef MSSpectrum item3
        for item3 in theo_spectra:
            v3.push_back(deref(item3.inst.get()))
    
        cdef double _r = self.inst.get().computePScore((<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), deref(v2), deref(v3), (<double>mz_window))
        cdef libcpp_vector[_MSSpectrum].iterator it_theo_spectra = v3.begin()
        replace_0 = []
        while it_theo_spectra != v3.end():
            item3 = MSSpectrum.__new__(MSSpectrum)
            item3.inst = shared_ptr[_MSSpectrum](new _MSSpectrum(deref(it_theo_spectra)))
            replace_0.append(item3)
            inc(it_theo_spectra)
        theo_spectra[:] = replace_0
        del v3
        replace = dict()
        cdef libcpp_map[size_t, _MSSpectrum].iterator it_peak_level_spectra = v2.begin()
        cdef MSSpectrum item_peak_level_spectra
        while it_peak_level_spectra != v2.end():
           item_peak_level_spectra = MSSpectrum.__new__(MSSpectrum)
           item_peak_level_spectra.inst = shared_ptr[_MSSpectrum](new _MSSpectrum((deref(it_peak_level_spectra)).second))
           replace[<size_t> deref(it_peak_level_spectra).first] = item_peak_level_spectra
           inc(it_peak_level_spectra)
        peak_level_spectra.clear()
        peak_level_spectra.update(replace)
        del v2
        py_result = <double>_r
        return py_result
    
    def _computePScore_1(self, double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , dict peak_level_spectra , MSSpectrum theo_spectrum , double mz_window ):
        """
        _computePScore_1(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectrum: MSSpectrum , mz_window: float ) -> float
        Computes the PScore for a single theoretical spectrum
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank)
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
        """
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(peak_level_spectra, dict) and all(isinstance(k, int) and k >= 0 for k in peak_level_spectra.keys()) and all(isinstance(v, MSSpectrum) for v in peak_level_spectra.values()), 'arg peak_level_spectra wrong type'
        assert isinstance(theo_spectrum, MSSpectrum), 'arg theo_spectrum wrong type'
        assert isinstance(mz_window, float), 'arg mz_window wrong type'
    
    
        cdef libcpp_map[size_t, _MSSpectrum] * v2 = new libcpp_map[size_t, _MSSpectrum]()
        for key, value in peak_level_spectra.items():
        
        
            deref(v2)[ (<size_t>key) ] = deref((<MSSpectrum>value).inst.get())
        
        
    
    
        cdef double _r = self.inst.get().computePScore((<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), deref(v2), (deref(theo_spectrum.inst.get())), (<double>mz_window))
        replace = dict()
        cdef libcpp_map[size_t, _MSSpectrum].iterator it_peak_level_spectra = v2.begin()
        cdef MSSpectrum item_peak_level_spectra
        while it_peak_level_spectra != v2.end():
           item_peak_level_spectra = MSSpectrum.__new__(MSSpectrum)
           item_peak_level_spectra.inst = shared_ptr[_MSSpectrum](new _MSSpectrum((deref(it_peak_level_spectra)).second))
           replace[<size_t> deref(it_peak_level_spectra).first] = item_peak_level_spectra
           inc(it_peak_level_spectra)
        peak_level_spectra.clear()
        peak_level_spectra.update(replace)
        del v2
        py_result = <double>_r
        return py_result
    
    def computePScore(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: computePScore(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectra: List[MSSpectrum] , mz_window: float ) -> float
          :noindex:
        
        Computes the PScore for a vector of theoretical spectra
        
        Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
        The best score obtained by scoring all those theoretical spectra against the experimental ones is returned
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank).
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
        
        .. rubric:: Overload:
        .. py:function:: computePScore(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectrum: MSSpectrum , mz_window: float ) -> float
          :noindex:
        
        Computes the PScore for a single theoretical spectrum
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank)
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
    
        """
        if (len(args)==5) and (isinstance(args[0], float)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], dict) and all(isinstance(k, int) and k >= 0 for k in args[2].keys()) and all(isinstance(v, MSSpectrum) for v in args[2].values())) and (isinstance(args[3], list) and all(isinstance(elemt_rec, MSSpectrum) for elemt_rec in args[3])) and (isinstance(args[4], float)):
            return self._computePScore_0(*args)
        elif (len(args)==5) and (isinstance(args[0], float)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], dict) and all(isinstance(k, int) and k >= 0 for k in args[2].keys()) and all(isinstance(v, MSSpectrum) for v in args[2].values())) and (isinstance(args[3], MSSpectrum)) and (isinstance(args[4], float)):
            return self._computePScore_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ProteaseDigestion:
    """
    Cython implementation of _ProteaseDigestion

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProteaseDigestion.html>`_
      -- Inherits from ['EnzymaticDigestion']

    Class for the enzymatic digestion of proteins
    
    Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
    Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
    due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
    peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
    Thus no random selection of just n specific missed cleavage sites is performed.
    If specificity is set to semi-specific, digestion also returns semi-specific products,
    i.e. with only one end at actual cleavage sites.
    
    Usage:
    
    .. code-block:: python
    
          from pyopenms import *
          from urllib.request import urlretrieve
          #
          urlretrieve ("http://www.uniprot.org/uniprot/P02769.fasta", "bsa.fasta")
          #
          dig = ProteaseDigestion()
          dig.setEnzyme('Lys-C')
          bsa_string = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
          bsa_oms_string = String(bsa_string) # convert python string to OpenMS::String for further processing
          #
          minlen = 6
          maxlen = 30
          #
          # Using AASequence and digest
          result_digest = []
          result_digest_min_max = []
          bsa_aaseq = AASequence.fromString(bsa_oms_string)
          dig.digest(bsa_aaseq, result_digest)
          dig.digest(bsa_aaseq, result_digest_min_max, minlen, maxlen)
          print(result_digest[4].toString()) # GLVLIAFSQYLQQCPFDEHVK
          print(len(result_digest)) # 57 peptides
          print(result_digest_min_max[4].toString()) # LVNELTEFAK
          print(len(result_digest_min_max)) # 42 peptides
          #
          # Semi-specific digestion
          result_semispecific = []
          dig.setSpecificity(EnzymaticDigestion.SPEC_SEMI)
          dig.digest(bsa_aaseq, result_semispecific)
          #
          # Using digestUnmodified without the need for AASequence from the EnzymaticDigestion base class
          result_digest_unmodified = []
          dig.digestUnmodified(StringView(bsa_oms_string), result_digest_unmodified, minlen, maxlen)
          print(result_digest_unmodified[4].getString()) # LVNELTEFAK
          print(len(result_digest_unmodified)) # 42 peptides
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ProteaseDigestion rv = ProteaseDigestion.__new__(ProteaseDigestion)
       rv.inst = shared_ptr[_ProteaseDigestion](new _ProteaseDigestion(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProteaseDigestion rv = ProteaseDigestion.__new__(ProteaseDigestion)
       rv.inst = shared_ptr[_ProteaseDigestion](new _ProteaseDigestion(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProteaseDigestion](new _ProteaseDigestion())
    
    def _init_1(self, ProteaseDigestion in_0 ):
        """
        _init_1(self, in_0: ProteaseDigestion ) -> None
        """
        assert isinstance(in_0, ProteaseDigestion), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProteaseDigestion](new _ProteaseDigestion((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProteaseDigestion ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProteaseDigestion)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setEnzyme_0(self,  name ):
        """
        _setEnzyme_0(self, name: Union[bytes, str, String] ) -> None
        Sets the enzyme for the digestion (by name)
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setEnzyme(deref((convString(name)).get()))
    
    def _setEnzyme_1(self, DigestionEnzyme enzyme ):
        """
        _setEnzyme_1(self, enzyme: DigestionEnzyme ) -> None
        Sets the enzyme for the digestion
        """
        assert isinstance(enzyme, DigestionEnzyme), 'arg enzyme wrong type'
    
        self.inst.get().setEnzyme((enzyme.inst.get()))
    
    def setEnzyme(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setEnzyme(self, name: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the enzyme for the digestion (by name)

        
        .. rubric:: Overload:
        .. py:function:: setEnzyme(self, enzyme: DigestionEnzyme ) -> None
          :noindex:
        
        Sets the enzyme for the digestion
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setEnzyme_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DigestionEnzyme)):
            return self._setEnzyme_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _digest_0(self, AASequence protein , list output ):
        """
        _digest_0(self, protein: AASequence , output: List[AASequence] ) -> int
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, AASequence) for elemt_rec in output), 'arg output wrong type'
    
        cdef libcpp_vector[_AASequence] * v1 = new libcpp_vector[_AASequence]()
        cdef AASequence item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
        cdef size_t _r = self.inst.get().digest((deref(protein.inst.get())), deref(v1))
        cdef libcpp_vector[_AASequence].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = AASequence.__new__(AASequence)
            item1.inst = shared_ptr[_AASequence](new _AASequence(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
        py_result = <size_t>_r
        return py_result
    
    def _digest_1(self, AASequence protein , list output ,  min_length ,  max_length ):
        """
        _digest_1(self, protein: AASequence , output: List[AASequence] , min_length: int , max_length: int ) -> int
          Performs the enzymatic digestion of a protein.
        
        
          :param protein: Sequence to digest
          :param output: Digestion products (peptides)
          :param min_length: Minimal length of reported products
          :param max_length: Maximal length of reported products (0 = no restriction)
          :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, AASequence) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(min_length, int) and min_length >= 0, 'arg min_length wrong type'
        assert isinstance(max_length, int) and max_length >= 0, 'arg max_length wrong type'
    
        cdef libcpp_vector[_AASequence] * v1 = new libcpp_vector[_AASequence]()
        cdef AASequence item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
    
    
        cdef size_t _r = self.inst.get().digest((deref(protein.inst.get())), deref(v1), (<size_t>min_length), (<size_t>max_length))
        cdef libcpp_vector[_AASequence].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = AASequence.__new__(AASequence)
            item1.inst = shared_ptr[_AASequence](new _AASequence(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
        py_result = <size_t>_r
        return py_result
    
    def digest(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: digest(self, protein: AASequence , output: List[AASequence] ) -> int
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: digest(self, protein: AASequence , output: List[AASequence] , min_length: int , max_length: int ) -> int
          :noindex:
        
          Performs the enzymatic digestion of a protein.
        
        
          :param protein: Sequence to digest
          :param output: Digestion products (peptides)
          :param min_length: Minimal length of reported products
          :param max_length: Maximal length of reported products (0 = no restriction)
          :return: Number of discarded digestion products (which are not matching length restrictions)
    
        """
        if (len(args)==2) and (isinstance(args[0], AASequence)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, AASequence) for elemt_rec in args[1])):
            return self._digest_0(*args)
        elif (len(args)==4) and (isinstance(args[0], AASequence)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, AASequence) for elemt_rec in args[1])) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], int) and args[3] >= 0):
            return self._digest_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def peptideCount(self, AASequence protein ):
        """
        peptideCount(self, protein: AASequence ) -> int
        Returns the number of peptides a digestion of protein would yield under the current enzyme and missed cleavage settings
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
    
        cdef size_t _r = self.inst.get().peptideCount((deref(protein.inst.get())))
        py_result = <size_t>_r
        return py_result
    
    def _isValidProduct_0(self, AASequence protein ,  pep_pos ,  pep_length , bool ignore_missed_cleavages , bool methionine_cleavage ):
        """
        _isValidProduct_0(self, protein: AASequence , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool
          Variant of EnzymaticDigestion::isValidProduct() with support for n-term protein cleavage and random D|P cleavage
        
          Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the flags provided here
        
        
          :param protein: Protein sequence
          :param pep_pos: Starting index of potential peptide
          :param pep_length: Length of potential peptide
          :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
          :param allow_nterm_protein_cleavage: Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 and protein starts with 'M'
          :param allow_random_asp_pro_cleavage: Allow cleavage at D|P sites to count as n/c-terminal
          :return: True if peptide has correct n/c terminals (according to enzyme, specificity and above flags)
        """
        assert isinstance(protein, AASequence), 'arg protein wrong type'
        assert isinstance(pep_pos, int) and pep_pos >= 0, 'arg pep_pos wrong type'
        assert isinstance(pep_length, int) and pep_length >= 0, 'arg pep_length wrong type'
        assert isinstance(ignore_missed_cleavages, pybool_t), 'arg ignore_missed_cleavages wrong type'
        assert isinstance(methionine_cleavage, pybool_t), 'arg methionine_cleavage wrong type'
    
    
    
    
    
        cdef bool _r = self.inst.get().isValidProduct((deref(protein.inst.get())), (<size_t>pep_pos), (<size_t>pep_length), (<bool>ignore_missed_cleavages), (<bool>methionine_cleavage))
        py_result = <bool>_r
        return py_result
    
    def _isValidProduct_1(self,  protein ,  pep_pos ,  pep_length , bool ignore_missed_cleavages , bool methionine_cleavage ):
        """
        _isValidProduct_1(self, protein: Union[bytes, str, String] , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool
        Forwards to isValidProduct using protein.toUnmodifiedString()
        """
        assert (isinstance(protein, str) or isinstance(protein, bytes) or isinstance(protein, String)), 'arg protein wrong type'
        assert isinstance(pep_pos, int) and pep_pos >= 0, 'arg pep_pos wrong type'
        assert isinstance(pep_length, int) and pep_length >= 0, 'arg pep_length wrong type'
        assert isinstance(ignore_missed_cleavages, pybool_t), 'arg ignore_missed_cleavages wrong type'
        assert isinstance(methionine_cleavage, pybool_t), 'arg methionine_cleavage wrong type'
    
    
    
    
    
        cdef bool _r = self.inst.get().isValidProduct(deref((convString(protein)).get()), (<size_t>pep_pos), (<size_t>pep_length), (<bool>ignore_missed_cleavages), (<bool>methionine_cleavage))
        py_result = <bool>_r
        return py_result
    
    def _isValidProduct_2(self,  sequence ,  pos ,  length , bool ignore_missed_cleavages ):
        """
        _isValidProduct_2(self, sequence: Union[bytes, str, String] , pos: int , length: int , ignore_missed_cleavages: bool ) -> bool
        Boolean operator returns true if the peptide fragment starting at position `pos` with length `length` within the sequence `sequence` generated by the current enzyme\n
        Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here
        
        
        :param protein: Protein sequence
        :param pep_pos: Starting index of potential peptide
        :param pep_length: Length of potential peptide
        :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
        :return: True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
        """
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
        assert isinstance(pos, int), 'arg pos wrong type'
        assert isinstance(length, int), 'arg length wrong type'
        assert isinstance(ignore_missed_cleavages, pybool_t), 'arg ignore_missed_cleavages wrong type'
    
    
    
    
        cdef bool _r = self.inst.get().isValidProduct(deref((convString(sequence)).get()), (<int>pos), (<int>length), (<bool>ignore_missed_cleavages))
        py_result = <bool>_r
        return py_result
    
    def isValidProduct(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: isValidProduct(self, protein: AASequence , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool
          :noindex:
        
          Variant of EnzymaticDigestion::isValidProduct() with support for n-term protein cleavage and random D|P cleavage
        
          Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the flags provided here
        
        
          :param protein: Protein sequence
          :param pep_pos: Starting index of potential peptide
          :param pep_length: Length of potential peptide
          :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
          :param allow_nterm_protein_cleavage: Regard peptide as n-terminal of protein if it starts only at pos=1 or 2 and protein starts with 'M'
          :param allow_random_asp_pro_cleavage: Allow cleavage at D|P sites to count as n/c-terminal
          :return: True if peptide has correct n/c terminals (according to enzyme, specificity and above flags)
        
        .. rubric:: Overload:
        .. py:function:: isValidProduct(self, protein: Union[bytes, str, String] , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool
          :noindex:
        
        Forwards to isValidProduct using protein.toUnmodifiedString()

        
        .. rubric:: Overload:
        .. py:function:: isValidProduct(self, sequence: Union[bytes, str, String] , pos: int , length: int , ignore_missed_cleavages: bool ) -> bool
          :noindex:
        
        Boolean operator returns true if the peptide fragment starting at position `pos` with length `length` within the sequence `sequence` generated by the current enzyme\n
        Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here
        
        
        :param protein: Protein sequence
        :param pep_pos: Starting index of potential peptide
        :param pep_length: Length of potential peptide
        :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
        :return: True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
    
        """
        if (len(args)==5) and (isinstance(args[0], AASequence)) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)):
            return self._isValidProduct_0(*args)
        elif (len(args)==5) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)):
            return self._isValidProduct_1(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int)) and (isinstance(args[2], int)) and (isinstance(args[3], pybool_t)):
            return self._isValidProduct_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getMissedCleavages(self):
        """
        getMissedCleavages(self) -> int
        Returns the max. number of allowed missed cleavages for the digestion
        """
        cdef size_t _r = self.inst.get().getMissedCleavages()
        py_result = <size_t>_r
        return py_result
    
    def setMissedCleavages(self,  missed_cleavages ):
        """
        setMissedCleavages(self, missed_cleavages: int ) -> None
        Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used
        """
        assert isinstance(missed_cleavages, int) and missed_cleavages >= 0, 'arg missed_cleavages wrong type'
    
        self.inst.get().setMissedCleavages((<size_t>missed_cleavages))
    
    def countInternalCleavageSites(self,  sequence ):
        """
        countInternalCleavageSites(self, sequence: Union[bytes, str, String] ) -> int
        Returns the number of internal cleavage sites for this sequence.
        """
        assert (isinstance(sequence, str) or isinstance(sequence, bytes) or isinstance(sequence, String)), 'arg sequence wrong type'
    
        cdef size_t _r = self.inst.get().countInternalCleavageSites(deref((convString(sequence)).get()))
        py_result = <size_t>_r
        return py_result
    
    def getEnzymeName(self):
        """
        getEnzymeName(self) -> Union[bytes, str, String]
        Returns the enzyme for the digestion
        """
        cdef _String _r = self.inst.get().getEnzymeName()
        py_result = convOutputString(_r)
        return py_result
    
    def getSpecificity(self):
        """
        getSpecificity(self) -> int
        Returns the specificity for the digestion
        """
        cdef _Specificity _r = self.inst.get().getSpecificity()
        py_result = <int>_r
        return py_result
    
    def setSpecificity(self, int spec ):
        """
        setSpecificity(self, spec: int ) -> None
        Sets the specificity for the digestion (default is SPEC_FULL)
        """
        assert spec in [0, 1, 2, 3, 8, 9, 10], 'arg spec wrong type'
    
        self.inst.get().setSpecificity((<_Specificity>spec))
    
    def getSpecificityByName(self,  name ):
        """
        getSpecificityByName(self, name: Union[bytes, str, String] ) -> int
        Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef _Specificity _r = self.inst.get().getSpecificityByName(deref((convString(name)).get()))
        py_result = <int>_r
        return py_result
    
    def digestUnmodified(self, StringView sequence , list output ,  min_length ,  max_length ):
        """
        digestUnmodified(self, sequence: StringView , output: List[StringView] , min_length: int , max_length: int ) -> int
        Performs the enzymatic digestion of an unmodified sequence\n
        By returning only references into the original string this is very fast
        
        
        :param sequence: Sequence to digest
        :param output: Digestion products
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        assert isinstance(sequence, StringView), 'arg sequence wrong type'
        assert isinstance(output, list) and all(isinstance(elemt_rec, StringView) for elemt_rec in output), 'arg output wrong type'
        assert isinstance(min_length, int) and min_length >= 0, 'arg min_length wrong type'
        assert isinstance(max_length, int) and max_length >= 0, 'arg max_length wrong type'
    
        cdef libcpp_vector[_StringView] * v1 = new libcpp_vector[_StringView]()
        cdef StringView item1
        for item1 in output:
            v1.push_back(deref(item1.inst.get()))
    
    
        cdef size_t _r = self.inst.get().digestUnmodified((deref(sequence.inst.get())), deref(v1), (<size_t>min_length), (<size_t>max_length))
        cdef libcpp_vector[_StringView].iterator it_output = v1.begin()
        replace_0 = []
        while it_output != v1.end():
            item1 = StringView.__new__(StringView)
            item1.inst = shared_ptr[_StringView](new _StringView(deref(it_output)))
            replace_0.append(item1)
            inc(it_output)
        output[:] = replace_0
        del v1
        py_result = <size_t>_r
        return py_result 

cdef class Residue:
    """
    Cython implementation of _Residue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Residue.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).getName().c_str() )

    
    def __copy__(self):
       cdef Residue rv = Residue.__new__(Residue)
       rv.inst = shared_ptr[_Residue](new _Residue(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Residue rv = Residue.__new__(Residue)
       rv.inst = shared_ptr[_Residue](new _Residue(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Residue](new _Residue())
    
    def _init_1(self, Residue in_0 ):
        """
        _init_1(self, in_0: Residue ) -> None
        """
        assert isinstance(in_0, Residue), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Residue](new _Residue((deref(in_0.inst.get()))))
    
    def _init_2(self,  name ,  three_letter_code ,  one_letter_code , EmpiricalFormula formula ):
        """
        _init_2(self, name: Union[bytes, str, String] , three_letter_code: Union[bytes, str, String] , one_letter_code: Union[bytes, str, String] , formula: EmpiricalFormula ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert (isinstance(three_letter_code, str) or isinstance(three_letter_code, bytes) or isinstance(three_letter_code, String)), 'arg three_letter_code wrong type'
        assert (isinstance(one_letter_code, str) or isinstance(one_letter_code, bytes) or isinstance(one_letter_code, String)), 'arg one_letter_code wrong type'
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
    
    
    
    
        self.inst = shared_ptr[_Residue](new _Residue(deref((convString(name)).get()), deref((convString(three_letter_code)).get()), deref((convString(one_letter_code)).get()), (deref(formula.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Residue ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, name: Union[bytes, str, String] , three_letter_code: Union[bytes, str, String] , one_letter_code: Union[bytes, str, String] , formula: EmpiricalFormula ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Residue)):
             self._init_1(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and ((isinstance(args[2], str) or isinstance(args[2], bytes) or isinstance(args[2], String))) and (isinstance(args[3], EmpiricalFormula)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getInternalToFull(self):
        """
        getInternalToFull(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToFull())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToNTerm(self):
        """
        getInternalToNTerm(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToNTerm())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToCTerm(self):
        """
        getInternalToCTerm(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToCTerm())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToAIon(self):
        """
        getInternalToAIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToAIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToBIon(self):
        """
        getInternalToBIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToBIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToCIon(self):
        """
        getInternalToCIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToCIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToXIon(self):
        """
        getInternalToXIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToXIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToYIon(self):
        """
        getInternalToYIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToYIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getInternalToZIon(self):
        """
        getInternalToZIon(self) -> EmpiricalFormula
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getInternalToZIon())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getResidueTypeName(self, int res_type ):
        """
        getResidueTypeName(self, res_type: int ) -> Union[bytes, str, String]
        Returns the ion name given as a residue type
        """
        assert res_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg res_type wrong type'
    
        cdef _String _r = self.inst.get().getResidueTypeName((<_ResidueType>res_type))
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the residue
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the residue
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setSynonyms(self, set synonyms ):
        """
        setSynonyms(self, synonyms: Set[bytes] ) -> None
        Sets the synonyms
        """
        assert isinstance(synonyms, set) and all(isinstance(i, bytes) for i in synonyms), 'arg synonyms wrong type'
        cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
        cdef bytes item0
        for item0 in synonyms:
           v0.insert(_String(<char *>item0))
        self.inst.get().setSynonyms(deref(v0))
        del v0
    
    def addSynonym(self,  synonym ):
        """
        addSynonym(self, synonym: Union[bytes, str, String] ) -> None
        Adds a synonym
        """
        assert (isinstance(synonym, str) or isinstance(synonym, bytes) or isinstance(synonym, String)), 'arg synonym wrong type'
    
        self.inst.get().addSynonym(deref((convString(synonym)).get()))
    
    def getSynonyms(self):
        """
        getSynonyms(self) -> Set[bytes]
        Returns the sysnonyms
        """
        _r = self.inst.get().getSynonyms()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setThreeLetterCode(self,  three_letter_code ):
        """
        setThreeLetterCode(self, three_letter_code: Union[bytes, str, String] ) -> None
        Sets the name of the residue as three letter code
        """
        assert (isinstance(three_letter_code, str) or isinstance(three_letter_code, bytes) or isinstance(three_letter_code, String)), 'arg three_letter_code wrong type'
    
        self.inst.get().setThreeLetterCode(deref((convString(three_letter_code)).get()))
    
    def getThreeLetterCode(self):
        """
        getThreeLetterCode(self) -> Union[bytes, str, String]
        Returns the name of the residue as three letter code
        """
        cdef _String _r = self.inst.get().getThreeLetterCode()
        py_result = convOutputString(_r)
        return py_result
    
    def setOneLetterCode(self,  one_letter_code ):
        """
        setOneLetterCode(self, one_letter_code: Union[bytes, str, String] ) -> None
        Sets the name as one letter code
        """
        assert (isinstance(one_letter_code, str) or isinstance(one_letter_code, bytes) or isinstance(one_letter_code, String)), 'arg one_letter_code wrong type'
    
        self.inst.get().setOneLetterCode(deref((convString(one_letter_code)).get()))
    
    def getOneLetterCode(self):
        """
        getOneLetterCode(self) -> Union[bytes, str, String]
        Returns the name as one letter code
        """
        cdef _String _r = self.inst.get().getOneLetterCode()
        py_result = convOutputString(_r)
        return py_result
    
    def addLossFormula(self, EmpiricalFormula in_0 ):
        """
        addLossFormula(self, in_0: EmpiricalFormula ) -> None
        Adds a neutral loss formula
        """
        assert isinstance(in_0, EmpiricalFormula), 'arg in_0 wrong type'
    
        self.inst.get().addLossFormula((deref(in_0.inst.get())))
    
    def setLossFormulas(self, list in_0 ):
        """
        setLossFormulas(self, in_0: List[EmpiricalFormula] ) -> None
        Sets the neutral loss formulas
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, EmpiricalFormula) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_EmpiricalFormula] * v0 = new libcpp_vector[_EmpiricalFormula]()
        cdef EmpiricalFormula item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setLossFormulas(deref(v0))
        del v0
    
    def addNTermLossFormula(self, EmpiricalFormula in_0 ):
        """
        addNTermLossFormula(self, in_0: EmpiricalFormula ) -> None
        Adds N-terminal losses
        """
        assert isinstance(in_0, EmpiricalFormula), 'arg in_0 wrong type'
    
        self.inst.get().addNTermLossFormula((deref(in_0.inst.get())))
    
    def setNTermLossFormulas(self, list in_0 ):
        """
        setNTermLossFormulas(self, in_0: List[EmpiricalFormula] ) -> None
        Sets the N-terminal losses
        """
        assert isinstance(in_0, list) and all(isinstance(elemt_rec, EmpiricalFormula) for elemt_rec in in_0), 'arg in_0 wrong type'
        cdef libcpp_vector[_EmpiricalFormula] * v0 = new libcpp_vector[_EmpiricalFormula]()
        cdef EmpiricalFormula item0
        for item0 in in_0:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setNTermLossFormulas(deref(v0))
        del v0
    
    def getLossFormulas(self):
        """
        getLossFormulas(self) -> List[EmpiricalFormula]
        Returns the neutral loss formulas
        """
        _r = self.inst.get().getLossFormulas()
        py_result = []
        cdef libcpp_vector[_EmpiricalFormula].iterator it__r = _r.begin()
        cdef EmpiricalFormula item_py_result
        while it__r != _r.end():
           item_py_result = EmpiricalFormula.__new__(EmpiricalFormula)
           item_py_result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getNTermLossFormulas(self):
        """
        getNTermLossFormulas(self) -> List[EmpiricalFormula]
        Returns N-terminal loss formulas
        """
        _r = self.inst.get().getNTermLossFormulas()
        py_result = []
        cdef libcpp_vector[_EmpiricalFormula].iterator it__r = _r.begin()
        cdef EmpiricalFormula item_py_result
        while it__r != _r.end():
           item_py_result = EmpiricalFormula.__new__(EmpiricalFormula)
           item_py_result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setLossNames(self, list name ):
        """
        setLossNames(self, name: List[bytes] ) -> None
        Sets the neutral loss molecule name
        """
        assert isinstance(name, list) and all(isinstance(i, bytes) for i in name), 'arg name wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in name:
           v0.push_back(_String(<char *>item0))
        self.inst.get().setLossNames(deref(v0))
        del v0
    
    def setNTermLossNames(self, list name ):
        """
        setNTermLossNames(self, name: List[bytes] ) -> None
        Sets the N-terminal loss names
        """
        assert isinstance(name, list) and all(isinstance(i, bytes) for i in name), 'arg name wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in name:
           v0.push_back(_String(<char *>item0))
        self.inst.get().setNTermLossNames(deref(v0))
        del v0
    
    def addLossName(self,  name ):
        """
        addLossName(self, name: Union[bytes, str, String] ) -> None
        Adds neutral loss molecule name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().addLossName(deref((convString(name)).get()))
    
    def addNTermLossName(self,  name ):
        """
        addNTermLossName(self, name: Union[bytes, str, String] ) -> None
        Adds a N-terminal loss name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().addNTermLossName(deref((convString(name)).get()))
    
    def getLossNames(self):
        """
        getLossNames(self) -> List[bytes]
        Gets neutral loss name (if there is one, else returns an empty string)
        """
        _r = self.inst.get().getLossNames()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def getNTermLossNames(self):
        """
        getNTermLossNames(self) -> List[bytes]
        Returns the N-terminal loss names
        """
        _r = self.inst.get().getNTermLossNames()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def setFormula(self, EmpiricalFormula formula ):
        """
        setFormula(self, formula: EmpiricalFormula ) -> None
        Sets empirical formula of the residue (must be full, with N and C-terminus)
        """
        assert isinstance(formula, EmpiricalFormula), 'arg formula wrong type'
    
        self.inst.get().setFormula((deref(formula.inst.get())))
    
    def _getFormula_0(self):
        """
        _getFormula_0(self) -> EmpiricalFormula
        Returns the empirical formula of the residue
        """
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula())
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def _getFormula_1(self, int res_type ):
        """
        _getFormula_1(self, res_type: int ) -> EmpiricalFormula
        """
        assert res_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg res_type wrong type'
    
        cdef _EmpiricalFormula * _r = new _EmpiricalFormula(self.inst.get().getFormula((<_ResidueType>res_type)))
        cdef EmpiricalFormula py_result = EmpiricalFormula.__new__(EmpiricalFormula)
        py_result.inst = shared_ptr[_EmpiricalFormula](_r)
        return py_result
    
    def getFormula(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getFormula(self, ) -> EmpiricalFormula
          :noindex:
        
        Returns the empirical formula of the residue

        
        .. rubric:: Overload:
        .. py:function:: getFormula(self, res_type: int ) -> EmpiricalFormula
          :noindex:
    
        """
        if not args:
            return self._getFormula_0(*args)
        elif (len(args)==1) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]):
            return self._getFormula_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setAverageWeight(self, double weight ):
        """
        setAverageWeight(self, weight: float ) -> None
        Sets average weight of the residue (must be full, with N and C-terminus)
        """
        assert isinstance(weight, float), 'arg weight wrong type'
    
        self.inst.get().setAverageWeight((<double>weight))
    
    def _getAverageWeight_0(self):
        """
        _getAverageWeight_0(self) -> float
        Returns average weight of the residue
        """
        cdef double _r = self.inst.get().getAverageWeight()
        py_result = <double>_r
        return py_result
    
    def _getAverageWeight_1(self, int res_type ):
        """
        _getAverageWeight_1(self, res_type: int ) -> float
        """
        assert res_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg res_type wrong type'
    
        cdef double _r = self.inst.get().getAverageWeight((<_ResidueType>res_type))
        py_result = <double>_r
        return py_result
    
    def getAverageWeight(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getAverageWeight(self, ) -> float
          :noindex:
        
        Returns average weight of the residue

        
        .. rubric:: Overload:
        .. py:function:: getAverageWeight(self, res_type: int ) -> float
          :noindex:
    
        """
        if not args:
            return self._getAverageWeight_0(*args)
        elif (len(args)==1) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]):
            return self._getAverageWeight_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMonoWeight(self, double weight ):
        """
        setMonoWeight(self, weight: float ) -> None
        Sets monoisotopic weight of the residue (must be full, with N and C-terminus)
        """
        assert isinstance(weight, float), 'arg weight wrong type'
    
        self.inst.get().setMonoWeight((<double>weight))
    
    def _getMonoWeight_0(self):
        """
        _getMonoWeight_0(self) -> float
        Returns monoisotopic weight of the residue
        """
        cdef double _r = self.inst.get().getMonoWeight()
        py_result = <double>_r
        return py_result
    
    def _getMonoWeight_1(self, int res_type ):
        """
        _getMonoWeight_1(self, res_type: int ) -> float
        """
        assert res_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg res_type wrong type'
    
        cdef double _r = self.inst.get().getMonoWeight((<_ResidueType>res_type))
        py_result = <double>_r
        return py_result
    
    def getMonoWeight(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getMonoWeight(self, ) -> float
          :noindex:
        
        Returns monoisotopic weight of the residue

        
        .. rubric:: Overload:
        .. py:function:: getMonoWeight(self, res_type: int ) -> float
          :noindex:
    
        """
        if not args:
            return self._getMonoWeight_0(*args)
        elif (len(args)==1) and (args[0] in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]):
            return self._getMonoWeight_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getModification(self):
        """
        getModification(self) -> ResidueModification
        """
        cdef const _ResidueModification * __r = (self.inst.get().getModification())
        if __r == NULL:
            return None
        cdef _ResidueModification * _r = new _ResidueModification(deref(__r))
        cdef ResidueModification py_result = ResidueModification.__new__(ResidueModification)
        py_result.inst = shared_ptr[_ResidueModification](_r)
        return py_result
    
    def _setModification_0(self,  name ):
        """
        _setModification_0(self, name: Union[bytes, str, String] ) -> None
        Sets the modification by name; the mod should be present in ModificationsDB
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setModification(deref((convString(name)).get()))
    
    def _setModification_1(self, ResidueModification mod ):
        """
        _setModification_1(self, mod: ResidueModification ) -> None
        Sets the modification by a ResidueModification object; checks if present in ModificationsDB and adds if not.
        """
        assert isinstance(mod, ResidueModification), 'arg mod wrong type'
    
        self.inst.get().setModification((deref(mod.inst.get())))
    
    def setModification(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setModification(self, name: Union[bytes, str, String] ) -> None
          :noindex:
        
        Sets the modification by name; the mod should be present in ModificationsDB

        
        .. rubric:: Overload:
        .. py:function:: setModification(self, mod: ResidueModification ) -> None
          :noindex:
        
        Sets the modification by a ResidueModification object; checks if present in ModificationsDB and adds if not.
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._setModification_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ResidueModification)):
            return self._setModification_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setModificationByDiffMonoMass(self, double diffMonoMass ):
        """
        setModificationByDiffMonoMass(self, diffMonoMass: float ) -> None
        Sets the modification by monoisotopic mass difference in Da; checks if present in ModificationsDB with tolerance and adds a "user-defined" modification if not (for later lookups).
        """
        assert isinstance(diffMonoMass, float), 'arg diffMonoMass wrong type'
    
        self.inst.get().setModificationByDiffMonoMass((<double>diffMonoMass))
    
    def getModificationName(self):
        """
        getModificationName(self) -> Union[bytes, str, String]
        Returns the name of the modification to the modification
        """
        cdef _String _r = self.inst.get().getModificationName()
        py_result = convOutputString(_r)
        return py_result
    
    def setLowMassIons(self, list low_mass_ions ):
        """
        setLowMassIons(self, low_mass_ions: List[EmpiricalFormula] ) -> None
        Sets the low mass marker ions as a vector of formulas
        """
        assert isinstance(low_mass_ions, list) and all(isinstance(elemt_rec, EmpiricalFormula) for elemt_rec in low_mass_ions), 'arg low_mass_ions wrong type'
        cdef libcpp_vector[_EmpiricalFormula] * v0 = new libcpp_vector[_EmpiricalFormula]()
        cdef EmpiricalFormula item0
        for item0 in low_mass_ions:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setLowMassIons(deref(v0))
        del v0
    
    def getLowMassIons(self):
        """
        getLowMassIons(self) -> List[EmpiricalFormula]
        Returns a vector of formulas with the low mass markers of the residue
        """
        _r = self.inst.get().getLowMassIons()
        py_result = []
        cdef libcpp_vector[_EmpiricalFormula].iterator it__r = _r.begin()
        cdef EmpiricalFormula item_py_result
        while it__r != _r.end():
           item_py_result = EmpiricalFormula.__new__(EmpiricalFormula)
           item_py_result.inst = shared_ptr[_EmpiricalFormula](new _EmpiricalFormula(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setResidueSets(self, set residues_sets ):
        """
        setResidueSets(self, residues_sets: Set[bytes] ) -> None
        Sets the residue sets the amino acid is contained in
        """
        assert isinstance(residues_sets, set) and all(isinstance(i, bytes) for i in residues_sets), 'arg residues_sets wrong type'
        cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
        cdef bytes item0
        for item0 in residues_sets:
           v0.insert(_String(<char *>item0))
        self.inst.get().setResidueSets(deref(v0))
        del v0
    
    def addResidueSet(self,  residue_sets ):
        """
        addResidueSet(self, residue_sets: Union[bytes, str, String] ) -> None
        Adds a residue set to the residue sets
        """
        assert (isinstance(residue_sets, str) or isinstance(residue_sets, bytes) or isinstance(residue_sets, String)), 'arg residue_sets wrong type'
    
        self.inst.get().addResidueSet(deref((convString(residue_sets)).get()))
    
    def getResidueSets(self):
        """
        getResidueSets(self) -> Set[bytes]
        Returns the residue sets this residue is contained in
        """
        _r = self.inst.get().getResidueSets()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def hasNeutralLoss(self):
        """
        hasNeutralLoss(self) -> bool
        True if the residue has neutral loss
        """
        cdef bool _r = self.inst.get().hasNeutralLoss()
        py_result = <bool>_r
        return py_result
    
    def hasNTermNeutralLosses(self):
        """
        hasNTermNeutralLosses(self) -> bool
        True if N-terminal neutral losses are set
        """
        cdef bool _r = self.inst.get().hasNTermNeutralLosses()
        py_result = <bool>_r
        return py_result
    
    def getPka(self):
        """
        getPka(self) -> float
        Returns the pka of the residue
        """
        cdef double _r = self.inst.get().getPka()
        py_result = <double>_r
        return py_result
    
    def getPkb(self):
        """
        getPkb(self) -> float
        Returns the pkb of the residue
        """
        cdef double _r = self.inst.get().getPkb()
        py_result = <double>_r
        return py_result
    
    def getPkc(self):
        """
        getPkc(self) -> float
        Returns the pkc of the residue if it exists otherwise -1
        """
        cdef double _r = self.inst.get().getPkc()
        py_result = <double>_r
        return py_result
    
    def getPiValue(self):
        """
        getPiValue(self) -> float
        Calculates the isoelectric point using the pk values
        """
        cdef double _r = self.inst.get().getPiValue()
        py_result = <double>_r
        return py_result
    
    def setPka(self, double value ):
        """
        setPka(self, value: float ) -> None
        Sets the pka of the residue
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setPka((<double>value))
    
    def setPkb(self, double value ):
        """
        setPkb(self, value: float ) -> None
        Sets the pkb of the residue
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setPkb((<double>value))
    
    def setPkc(self, double value ):
        """
        setPkc(self, value: float ) -> None
        Sets the pkc of the residue
        """
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().setPkc((<double>value))
    
    def getSideChainBasicity(self):
        """
        getSideChainBasicity(self) -> float
        Returns the side chain basicity
        """
        cdef double _r = self.inst.get().getSideChainBasicity()
        py_result = <double>_r
        return py_result
    
    def setSideChainBasicity(self, double gb_sc ):
        """
        setSideChainBasicity(self, gb_sc: float ) -> None
        Sets the side chain basicity
        """
        assert isinstance(gb_sc, float), 'arg gb_sc wrong type'
    
        self.inst.get().setSideChainBasicity((<double>gb_sc))
    
    def getBackboneBasicityLeft(self):
        """
        getBackboneBasicityLeft(self) -> float
        Returns the backbone basicitiy if located in N-terminal direction
        """
        cdef double _r = self.inst.get().getBackboneBasicityLeft()
        py_result = <double>_r
        return py_result
    
    def setBackboneBasicityLeft(self, double gb_bb_l ):
        """
        setBackboneBasicityLeft(self, gb_bb_l: float ) -> None
        Sets the N-terminal direction backbone basicitiy
        """
        assert isinstance(gb_bb_l, float), 'arg gb_bb_l wrong type'
    
        self.inst.get().setBackboneBasicityLeft((<double>gb_bb_l))
    
    def getBackboneBasicityRight(self):
        """
        getBackboneBasicityRight(self) -> float
        Returns the C-terminal direction backbone basicitiy
        """
        cdef double _r = self.inst.get().getBackboneBasicityRight()
        py_result = <double>_r
        return py_result
    
    def setBackboneBasicityRight(self, double gb_bb_r ):
        """
        setBackboneBasicityRight(self, gb_bb_r: float ) -> None
        Sets the C-terminal direction backbone basicity
        """
        assert isinstance(gb_bb_r, float), 'arg gb_bb_r wrong type'
    
        self.inst.get().setBackboneBasicityRight((<double>gb_bb_r))
    
    def isModified(self):
        """
        isModified(self) -> bool
        True if the residue is a modified one
        """
        cdef bool _r = self.inst.get().isModified()
        py_result = <bool>_r
        return py_result
    
    def isInResidueSet(self,  residue_set ):
        """
        isInResidueSet(self, residue_set: Union[bytes, str, String] ) -> bool
        True if the residue is contained in the set
        """
        assert (isinstance(residue_set, str) or isinstance(residue_set, bytes) or isinstance(residue_set, String)), 'arg residue_set wrong type'
    
        cdef bool _r = self.inst.get().isInResidueSet(deref((convString(residue_set)).get()))
        py_result = <bool>_r
        return py_result
    
    def residueTypeToIonLetter(self, int res_type ):
        """
        residueTypeToIonLetter(self, res_type: int ) -> Union[bytes, str, String]
        Helper for mapping residue types to letters for Text annotations and labels
        """
        assert res_type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17], 'arg res_type wrong type'
    
        cdef _String _r = self.inst.get().residueTypeToIonLetter((<_ResidueType>res_type))
        py_result = convOutputString(_r)
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, Residue):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef Residue other_casted = other
        cdef Residue self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    ResidueType = __ResidueType 

cdef class RichPeak2D:
    """
    Cython implementation of _RichPeak2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RichPeak2D.html>`_
      -- Inherits from ['Peak2D', 'UniqueIdInterface', 'MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef RichPeak2D rv = RichPeak2D.__new__(RichPeak2D)
       rv.inst = shared_ptr[_RichPeak2D](new _RichPeak2D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef RichPeak2D rv = RichPeak2D.__new__(RichPeak2D)
       rv.inst = shared_ptr[_RichPeak2D](new _RichPeak2D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        A 2-dimensional raw data point or peak with meta information
        """
        self.inst = shared_ptr[_RichPeak2D](new _RichPeak2D())
    
    def _init_1(self, RichPeak2D in_0 ):
        """
        _init_1(self, in_0: RichPeak2D ) -> None
        """
        assert isinstance(in_0, RichPeak2D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_RichPeak2D](new _RichPeak2D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        A 2-dimensional raw data point or peak with meta information

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: RichPeak2D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], RichPeak2D)):
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
        if not isinstance(other, RichPeak2D):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef RichPeak2D other_casted = other
        cdef RichPeak2D self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class SavitzkyGolayFilter:
    """
    Cython implementation of _SavitzkyGolayFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SavitzkyGolayFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SavitzkyGolayFilter rv = SavitzkyGolayFilter.__new__(SavitzkyGolayFilter)
       rv.inst = shared_ptr[_SavitzkyGolayFilter](new _SavitzkyGolayFilter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SavitzkyGolayFilter rv = SavitzkyGolayFilter.__new__(SavitzkyGolayFilter)
       rv.inst = shared_ptr[_SavitzkyGolayFilter](new _SavitzkyGolayFilter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SavitzkyGolayFilter](new _SavitzkyGolayFilter())
    
    def _init_1(self, SavitzkyGolayFilter in_0 ):
        """
        _init_1(self, in_0: SavitzkyGolayFilter ) -> None
        """
        assert isinstance(in_0, SavitzkyGolayFilter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SavitzkyGolayFilter](new _SavitzkyGolayFilter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SavitzkyGolayFilter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SavitzkyGolayFilter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filter(self, MSSpectrum spectrum ):
        """
        filter(self, spectrum: MSSpectrum ) -> None
        Removed the noise from an MSSpectrum containing profile data
        """
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
        self.inst.get().filter((deref(spectrum.inst.get())))
    
    def filterExperiment(self, MSExperiment exp ):
        """
        filterExperiment(self, exp: MSExperiment ) -> None
        Removed the noise from an MSExperiment containing profile data
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

cdef class SignalToNoiseEstimatorMeanIterative:
    """
    Cython implementation of _SignalToNoiseEstimatorMeanIterative[_MSSpectrum]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMeanIterative[_MSSpectrum].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SignalToNoiseEstimatorMeanIterative rv = SignalToNoiseEstimatorMeanIterative.__new__(SignalToNoiseEstimatorMeanIterative)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMeanIterative[_MSSpectrum]](new _SignalToNoiseEstimatorMeanIterative[_MSSpectrum](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SignalToNoiseEstimatorMeanIterative rv = SignalToNoiseEstimatorMeanIterative.__new__(SignalToNoiseEstimatorMeanIterative)
       rv.inst = shared_ptr[_SignalToNoiseEstimatorMeanIterative[_MSSpectrum]](new _SignalToNoiseEstimatorMeanIterative[_MSSpectrum](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SignalToNoiseEstimatorMeanIterative[_MSSpectrum]](new _SignalToNoiseEstimatorMeanIterative[_MSSpectrum]())
    
    def _init_1(self, SignalToNoiseEstimatorMeanIterative in_0 ):
        """
        _init_1(self, in_0: SignalToNoiseEstimatorMeanIterative ) -> None
        """
        assert isinstance(in_0, SignalToNoiseEstimatorMeanIterative), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SignalToNoiseEstimatorMeanIterative[_MSSpectrum]](new _SignalToNoiseEstimatorMeanIterative[_MSSpectrum]((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SignalToNoiseEstimatorMeanIterative ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SignalToNoiseEstimatorMeanIterative)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def init(self, MSSpectrum c ):
        """
        init(self, c: MSSpectrum ) -> None
        """
        assert isinstance(c, MSSpectrum), 'arg c wrong type'
    
        self.inst.get().init((deref(c.inst.get())))
    
    def getSignalToNoise(self,  index ):
        """
        getSignalToNoise(self, index: int ) -> float
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
    
        cdef double _r = self.inst.get().getSignalToNoise((<size_t>index))
        py_result = <double>_r
        return py_result
    IntensityThresholdCalculation = __IntensityThresholdCalculation 

cdef class TMTSixPlexQuantitationMethod:
    """
    Cython implementation of _TMTSixPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTSixPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef TMTSixPlexQuantitationMethod rv = TMTSixPlexQuantitationMethod.__new__(TMTSixPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTSixPlexQuantitationMethod](new _TMTSixPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef TMTSixPlexQuantitationMethod rv = TMTSixPlexQuantitationMethod.__new__(TMTSixPlexQuantitationMethod)
       rv.inst = shared_ptr[_TMTSixPlexQuantitationMethod](new _TMTSixPlexQuantitationMethod(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_TMTSixPlexQuantitationMethod](new _TMTSixPlexQuantitationMethod())
    
    def _init_1(self, TMTSixPlexQuantitationMethod in_0 ):
        """
        _init_1(self, in_0: TMTSixPlexQuantitationMethod ) -> None
        """
        assert isinstance(in_0, TMTSixPlexQuantitationMethod), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_TMTSixPlexQuantitationMethod](new _TMTSixPlexQuantitationMethod((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: TMTSixPlexQuantitationMethod ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], TMTSixPlexQuantitationMethod)):
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
