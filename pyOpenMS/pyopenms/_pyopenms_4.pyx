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

cdef class __CombinationsLogic:
    None
    OR = 0
    AND = 1
    XOR = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class QuotingMethod:
    None
    NONE = 0
    ESCAPE = 1
    DOUBLE = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __RequirementLevel:
    None
    MUST = 0
    SHOULD = 1
    MAY = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __XFDRAlgorithm_ExitCodes:
    None
    EXECUTION_OK = 0
    ILLEGAL_PARAMETERS = 1
    UNEXPECTED_RESULT = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class CVMappingRule:
    """
    Cython implementation of _CVMappingRule

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingRule.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CVMappingRule rv = CVMappingRule.__new__(CVMappingRule)
       rv.inst = shared_ptr[_CVMappingRule](new _CVMappingRule(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CVMappingRule rv = CVMappingRule.__new__(CVMappingRule)
       rv.inst = shared_ptr[_CVMappingRule](new _CVMappingRule(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CVMappingRule](new _CVMappingRule())
    
    def _init_1(self, CVMappingRule in_0 ):
        """
        _init_1(self, in_0: CVMappingRule ) -> None
        """
        assert isinstance(in_0, CVMappingRule), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CVMappingRule](new _CVMappingRule((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CVMappingRule ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CVMappingRule)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setIdentifier(self,  identifier ):
        """
        setIdentifier(self, identifier: Union[bytes, str, String] ) -> None
        Sets the identifier of the rule
        """
        assert (isinstance(identifier, str) or isinstance(identifier, bytes) or isinstance(identifier, String)), 'arg identifier wrong type'
    
        self.inst.get().setIdentifier(deref((convString(identifier)).get()))
    
    def getIdentifier(self):
        """
        getIdentifier(self) -> Union[bytes, str, String]
        Returns the identifier of the rule
        """
        cdef _String _r = self.inst.get().getIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setElementPath(self,  element_path ):
        """
        setElementPath(self, element_path: Union[bytes, str, String] ) -> None
        Sets the path of the DOM element, where this rule is allowed
        """
        assert (isinstance(element_path, str) or isinstance(element_path, bytes) or isinstance(element_path, String)), 'arg element_path wrong type'
    
        self.inst.get().setElementPath(deref((convString(element_path)).get()))
    
    def getElementPath(self):
        """
        getElementPath(self) -> Union[bytes, str, String]
        Returns the path of the DOM element, where this rule is allowed
        """
        cdef _String _r = self.inst.get().getElementPath()
        py_result = convOutputString(_r)
        return py_result
    
    def setRequirementLevel(self, int level ):
        """
        setRequirementLevel(self, level: int ) -> None
        Sets the requirement level of this rule
        """
        assert level in [0, 1, 2], 'arg level wrong type'
    
        self.inst.get().setRequirementLevel((<_RequirementLevel>level))
    
    def getRequirementLevel(self):
        """
        getRequirementLevel(self) -> int
        Returns the requirement level of this rule
        """
        cdef _RequirementLevel _r = self.inst.get().getRequirementLevel()
        py_result = <int>_r
        return py_result
    
    def setCombinationsLogic(self, int combinations_logic ):
        """
        setCombinationsLogic(self, combinations_logic: int ) -> None
        Sets the combination operator of the rule
        """
        assert combinations_logic in [0, 1, 2], 'arg combinations_logic wrong type'
    
        self.inst.get().setCombinationsLogic((<_CombinationsLogic>combinations_logic))
    
    def getCombinationsLogic(self):
        """
        getCombinationsLogic(self) -> int
        Returns the combinations operator of the rule
        """
        cdef _CombinationsLogic _r = self.inst.get().getCombinationsLogic()
        py_result = <int>_r
        return py_result
    
    def setScopePath(self,  path ):
        """
        setScopePath(self, path: Union[bytes, str, String] ) -> None
        Sets the scope path of the rule
        """
        assert (isinstance(path, str) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'
    
        self.inst.get().setScopePath(deref((convString(path)).get()))
    
    def getScopePath(self):
        """
        getScopePath(self) -> Union[bytes, str, String]
        Returns the scope path of the rule
        """
        cdef _String _r = self.inst.get().getScopePath()
        py_result = convOutputString(_r)
        return py_result
    
    def setCVTerms(self, list cv_terms ):
        """
        setCVTerms(self, cv_terms: List[CVMappingTerm] ) -> None
        Sets the terms which are allowed
        """
        assert isinstance(cv_terms, list) and all(isinstance(elemt_rec, CVMappingTerm) for elemt_rec in cv_terms), 'arg cv_terms wrong type'
        cdef libcpp_vector[_CVMappingTerm] * v0 = new libcpp_vector[_CVMappingTerm]()
        cdef CVMappingTerm item0
        for item0 in cv_terms:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setCVTerms(deref(v0))
        del v0
    
    def getCVTerms(self):
        """
        getCVTerms(self) -> List[CVMappingTerm]
        Returns the allowed terms
        """
        _r = self.inst.get().getCVTerms()
        py_result = []
        cdef libcpp_vector[_CVMappingTerm].iterator it__r = _r.begin()
        cdef CVMappingTerm item_py_result
        while it__r != _r.end():
           item_py_result = CVMappingTerm.__new__(CVMappingTerm)
           item_py_result.inst = shared_ptr[_CVMappingTerm](new _CVMappingTerm(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def addCVTerm(self, CVMappingTerm cv_terms ):
        """
        addCVTerm(self, cv_terms: CVMappingTerm ) -> None
        Adds a term to the allowed terms
        """
        assert isinstance(cv_terms, CVMappingTerm), 'arg cv_terms wrong type'
    
        self.inst.get().addCVTerm((deref(cv_terms.inst.get())))
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, CVMappingRule):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef CVMappingRule other_casted = other
        cdef CVMappingRule self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    CombinationsLogic = __CombinationsLogic
    RequirementLevel = __RequirementLevel 

cdef class Fitter1D:
    """
    Cython implementation of _Fitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Fitter1D.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Fitter1D rv = Fitter1D.__new__(Fitter1D)
       rv.inst = shared_ptr[_Fitter1D](new _Fitter1D(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Fitter1D rv = Fitter1D.__new__(Fitter1D)
       rv.inst = shared_ptr[_Fitter1D](new _Fitter1D(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Abstract base class for all 1D-dimensional model fitter
        """
        self.inst = shared_ptr[_Fitter1D](new _Fitter1D())
    
    def _init_1(self, Fitter1D in_0 ):
        """
        _init_1(self, in_0: Fitter1D ) -> None
        """
        assert isinstance(in_0, Fitter1D), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Fitter1D](new _Fitter1D((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Abstract base class for all 1D-dimensional model fitter

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Fitter1D ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Fitter1D)):
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

cdef class InspectOutfile:
    """
    Cython implementation of _InspectOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InspectOutfile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef InspectOutfile rv = InspectOutfile.__new__(InspectOutfile)
       rv.inst = shared_ptr[_InspectOutfile](new _InspectOutfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef InspectOutfile rv = InspectOutfile.__new__(InspectOutfile)
       rv.inst = shared_ptr[_InspectOutfile](new _InspectOutfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        This class serves to read in an Inspect outfile and write an idXML file
        """
        self.inst = shared_ptr[_InspectOutfile](new _InspectOutfile())
    
    def _init_1(self, InspectOutfile in_0 ):
        """
        _init_1(self, in_0: InspectOutfile ) -> None
        """
        assert isinstance(in_0, InspectOutfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_InspectOutfile](new _InspectOutfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        This class serves to read in an Inspect outfile and write an idXML file

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: InspectOutfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], InspectOutfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  result_filename , PeptideIdentificationList peptide_identifications , ProteinIdentification protein_identification , double p_value_threshold ,  database_filename ):
        """
        load(self, result_filename: Union[bytes, str, String] , peptide_identifications: PeptideIdentificationList , protein_identification: ProteinIdentification , p_value_threshold: float , database_filename: Union[bytes, str, String] ) -> List[int]
        Load the results of an Inspect search
        
        
        :param result_filename: Input parameter which is the file name of the input file
        :param peptide_identifications: Output parameter which holds the peptide identifications from the given file
        :param protein_identification: Output parameter which holds the protein identifications from the given file
        :param p_value_threshold:
        :param database_filename:
        :raises:
          Exception: FileNotFound is thrown if the given file could not be found
        :raises:
          Exception: ParseError is thrown if the given file could not be parsed
        :raises:
          Exception: FileEmpty is thrown if the given file is empty
        """
        assert (isinstance(result_filename, str) or isinstance(result_filename, bytes) or isinstance(result_filename, String)), 'arg result_filename wrong type'
        assert isinstance(peptide_identifications, PeptideIdentificationList), 'arg peptide_identifications wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(p_value_threshold, float), 'arg p_value_threshold wrong type'
        assert (isinstance(database_filename, str) or isinstance(database_filename, bytes) or isinstance(database_filename, String)), 'arg database_filename wrong type'
    
    
    
    
    
        _r = self.inst.get().load(deref((convString(result_filename)).get()), (deref(peptide_identifications.inst.get())), (deref(protein_identification.inst.get())), (<double>p_value_threshold), deref((convString(database_filename)).get()))
        cdef list py_result = _r
        return py_result
    
    def getWantedRecords(self,  result_filename , double p_value_threshold ):
        """
        getWantedRecords(self, result_filename: Union[bytes, str, String] , p_value_threshold: float ) -> List[int]
        Loads only results which exceeds a given p-value threshold
        
        
        :param result_filename: The filename of the results file
        :param p_value_threshold: Only identifications exceeding this threshold are read
        :raises:
          Exception: FileNotFound is thrown if the given file could not be found
        :raises:
          Exception: FileEmpty is thrown if the given file is empty
        """
        assert (isinstance(result_filename, str) or isinstance(result_filename, bytes) or isinstance(result_filename, String)), 'arg result_filename wrong type'
        assert isinstance(p_value_threshold, float), 'arg p_value_threshold wrong type'
    
    
        _r = self.inst.get().getWantedRecords(deref((convString(result_filename)).get()), (<double>p_value_threshold))
        cdef list py_result = _r
        return py_result
    
    def compressTrieDB(self,  database_filename ,  index_filename , list wanted_records ,  snd_database_filename ,  snd_index_filename , bool append ):
        """
        compressTrieDB(self, database_filename: Union[bytes, str, String] , index_filename: Union[bytes, str, String] , wanted_records: List[int] , snd_database_filename: Union[bytes, str, String] , snd_index_filename: Union[bytes, str, String] , append: bool ) -> None
        Generates a trie database from another one, using the wanted records only
        """
        assert (isinstance(database_filename, str) or isinstance(database_filename, bytes) or isinstance(database_filename, String)), 'arg database_filename wrong type'
        assert (isinstance(index_filename, str) or isinstance(index_filename, bytes) or isinstance(index_filename, String)), 'arg index_filename wrong type'
        assert isinstance(wanted_records, list) and all(isinstance(elemt_rec, int) and elemt_rec >= 0 for elemt_rec in wanted_records), 'arg wanted_records wrong type'
        assert (isinstance(snd_database_filename, str) or isinstance(snd_database_filename, bytes) or isinstance(snd_database_filename, String)), 'arg snd_database_filename wrong type'
        assert (isinstance(snd_index_filename, str) or isinstance(snd_index_filename, bytes) or isinstance(snd_index_filename, String)), 'arg snd_index_filename wrong type'
        assert isinstance(append, pybool_t), 'arg append wrong type'
    
    
        cdef libcpp_vector[size_t] v2 = wanted_records
    
    
    
        self.inst.get().compressTrieDB(deref((convString(database_filename)).get()), deref((convString(index_filename)).get()), v2, deref((convString(snd_database_filename)).get()), deref((convString(snd_index_filename)).get()), (<bool>append))
        wanted_records[:] = v2
    
    def generateTrieDB(self,  source_database_filename ,  database_filename ,  index_filename , bool append ,  species ):
        """
        generateTrieDB(self, source_database_filename: Union[bytes, str, String] , database_filename: Union[bytes, str, String] , index_filename: Union[bytes, str, String] , append: bool , species: Union[bytes, str, String] ) -> None
        Generates a trie database from a given one (the type of database is determined by getLabels)
        """
        assert (isinstance(source_database_filename, str) or isinstance(source_database_filename, bytes) or isinstance(source_database_filename, String)), 'arg source_database_filename wrong type'
        assert (isinstance(database_filename, str) or isinstance(database_filename, bytes) or isinstance(database_filename, String)), 'arg database_filename wrong type'
        assert (isinstance(index_filename, str) or isinstance(index_filename, bytes) or isinstance(index_filename, String)), 'arg index_filename wrong type'
        assert isinstance(append, pybool_t), 'arg append wrong type'
        assert (isinstance(species, str) or isinstance(species, bytes) or isinstance(species, String)), 'arg species wrong type'
    
    
    
    
    
        self.inst.get().generateTrieDB(deref((convString(source_database_filename)).get()), deref((convString(database_filename)).get()), deref((convString(index_filename)).get()), (<bool>append), deref((convString(species)).get()))
    
    def getACAndACType(self,  line ,  accession ,  accession_type ):
        """
        getACAndACType(self, line: Union[bytes, str, String] , accession: String , accession_type: String ) -> None
        Retrieve the accession type and accession number from a protein description line
        """
        assert (isinstance(line, str) or isinstance(line, bytes) or isinstance(line, String)), 'arg line wrong type'
        assert isinstance(accession, String), 'arg accession wrong type'
        assert isinstance(accession_type, String), 'arg accession_type wrong type'
    
    
    
        self.inst.get().getACAndACType(deref((convString(line)).get()), deref((<String>accession).inst.get()), deref((<String>accession_type).inst.get()))
    
    def getLabels(self,  source_database_filename ,  ac_label ,  sequence_start_label ,  sequence_end_label ,  comment_label ,  species_label ):
        """
        getLabels(self, source_database_filename: Union[bytes, str, String] , ac_label: String , sequence_start_label: String , sequence_end_label: String , comment_label: String , species_label: String ) -> None
        Retrieve the labels of a given database (at the moment FASTA and Swissprot)
        """
        assert (isinstance(source_database_filename, str) or isinstance(source_database_filename, bytes) or isinstance(source_database_filename, String)), 'arg source_database_filename wrong type'
        assert isinstance(ac_label, String), 'arg ac_label wrong type'
        assert isinstance(sequence_start_label, String), 'arg sequence_start_label wrong type'
        assert isinstance(sequence_end_label, String), 'arg sequence_end_label wrong type'
        assert isinstance(comment_label, String), 'arg comment_label wrong type'
        assert isinstance(species_label, String), 'arg species_label wrong type'
    
    
    
    
    
    
        self.inst.get().getLabels(deref((convString(source_database_filename)).get()), deref((<String>ac_label).inst.get()), deref((<String>sequence_start_label).inst.get()), deref((<String>sequence_end_label).inst.get()), deref((<String>comment_label).inst.get()), deref((<String>species_label).inst.get()))
    
    def getSequences(self,  database_filename , dict wanted_records , list sequences ):
        """
        getSequences(self, database_filename: Union[bytes, str, String] , wanted_records: Dict[int, int] , sequences: List[bytes] ) -> List[int]
        Retrieve sequences from a trie database
        """
        assert (isinstance(database_filename, str) or isinstance(database_filename, bytes) or isinstance(database_filename, String)), 'arg database_filename wrong type'
        assert isinstance(wanted_records, dict) and all(isinstance(k, int) and k >= 0 for k in wanted_records.keys()) and all(isinstance(v, int) and v >= 0 for v in wanted_records.values()), 'arg wanted_records wrong type'
        assert isinstance(sequences, list) and all(isinstance(i, bytes) for i in sequences), 'arg sequences wrong type'
    
        cdef libcpp_map[size_t, size_t] * v1 = new libcpp_map[size_t, size_t]()
        for key, value in wanted_records.items():
        
        
            deref(v1)[ (<size_t>key) ] = (<size_t>value)
        
        
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in sequences:
           v2.push_back(_String(<char *>item2))
        _r = self.inst.get().getSequences(deref((convString(database_filename)).get()), deref(v1), deref(v2))
        replace = []
        cdef libcpp_vector[_String].iterator it_2 = v2.begin()
        while it_2 != v2.end():
           replace.append(<char*>deref(it_2).c_str())
           inc(it_2)
        sequences[:] = replace
        del v2
        replace = dict()
        cdef libcpp_map[size_t, size_t].iterator it_wanted_records = v1.begin()
        while it_wanted_records != v1.end():
           replace[<size_t> deref(it_wanted_records).first] = <size_t> deref(it_wanted_records).second
           inc(it_wanted_records)
        wanted_records.clear()
        wanted_records.update(replace)
        del v1
        cdef list py_result = _r
        return py_result
    
    def getExperiment(self, MSExperiment exp ,  type_ ,  in_filename ):
        """
        getExperiment(self, exp: MSExperiment , type_: String , in_filename: Union[bytes, str, String] ) -> None
        Get the experiment from a file
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert isinstance(type_, String), 'arg type_ wrong type'
        assert (isinstance(in_filename, str) or isinstance(in_filename, bytes) or isinstance(in_filename, String)), 'arg in_filename wrong type'
    
    
    
        self.inst.get().getExperiment((deref(exp.inst.get())), deref((<String>type_).inst.get()), deref((convString(in_filename)).get()))
    
    def getSearchEngineAndVersion(self,  cmd_output , ProteinIdentification protein_identification ):
        """
        getSearchEngineAndVersion(self, cmd_output: Union[bytes, str, String] , protein_identification: ProteinIdentification ) -> bool
        Get the search engine and its version from the output of the InsPecT executable without parameters. Returns true on success, false otherwise
        """
        assert (isinstance(cmd_output, str) or isinstance(cmd_output, bytes) or isinstance(cmd_output, String)), 'arg cmd_output wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
    
    
        cdef bool _r = self.inst.get().getSearchEngineAndVersion(deref((convString(cmd_output)).get()), (deref(protein_identification.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def readOutHeader(self,  filename ,  header_line ,  spectrum_file_column ,  scan_column ,  peptide_column ,  protein_column ,  charge_column ,  MQ_score_column ,  p_value_column ,  record_number_column ,  DB_file_pos_column ,  spec_file_pos_column ,  number_of_columns ):
        """
        readOutHeader(self, filename: Union[bytes, str, String] , header_line: Union[bytes, str, String] , spectrum_file_column: int , scan_column: int , peptide_column: int , protein_column: int , charge_column: int , MQ_score_column: int , p_value_column: int , record_number_column: int , DB_file_pos_column: int , spec_file_pos_column: int , number_of_columns: int ) -> None
        Read the header of an inspect output file and retrieve various information
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert (isinstance(header_line, str) or isinstance(header_line, bytes) or isinstance(header_line, String)), 'arg header_line wrong type'
        assert isinstance(spectrum_file_column, int), 'arg spectrum_file_column wrong type'
        assert isinstance(scan_column, int), 'arg scan_column wrong type'
        assert isinstance(peptide_column, int), 'arg peptide_column wrong type'
        assert isinstance(protein_column, int), 'arg protein_column wrong type'
        assert isinstance(charge_column, int), 'arg charge_column wrong type'
        assert isinstance(MQ_score_column, int), 'arg MQ_score_column wrong type'
        assert isinstance(p_value_column, int), 'arg p_value_column wrong type'
        assert isinstance(record_number_column, int), 'arg record_number_column wrong type'
        assert isinstance(DB_file_pos_column, int), 'arg DB_file_pos_column wrong type'
        assert isinstance(spec_file_pos_column, int), 'arg spec_file_pos_column wrong type'
        assert isinstance(number_of_columns, int) and number_of_columns >= 0, 'arg number_of_columns wrong type'
    
    
    
    
    
    
    
    
    
    
    
    
    
        self.inst.get().readOutHeader(deref((convString(filename)).get()), deref((convString(header_line)).get()), (<int &>spectrum_file_column), (<int &>scan_column), (<int &>peptide_column), (<int &>protein_column), (<int &>charge_column), (<int &>MQ_score_column), (<int &>p_value_column), (<int &>record_number_column), (<int &>DB_file_pos_column), (<int &>spec_file_pos_column), (<size_t &>number_of_columns))
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, InspectOutfile):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef InspectOutfile other_casted = other
        cdef InspectOutfile self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class IntegerDataArray:
    """
    Cython implementation of _IntegerDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1IntegerDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra integer data attached to a spectrum or chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IntegerDataArray rv = IntegerDataArray.__new__(IntegerDataArray)
       rv.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IntegerDataArray rv = IntegerDataArray.__new__(IntegerDataArray)
       rv.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray())
    
    def _init_1(self, IntegerDataArray in_0 ):
        """
        _init_1(self, in_0: IntegerDataArray ) -> None
        """
        assert isinstance(in_0, IntegerDataArray), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IntegerDataArray](new _IntegerDataArray((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IntegerDataArray ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IntegerDataArray)):
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
    
    def resize(self,  n ):
        """
        resize(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().resize((<size_t>n))
    
    def reserve(self,  n ):
        """
        reserve(self, n: int ) -> None
        """
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
        self.inst.get().reserve((<size_t>n))
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def push_back(self,  in_0 ):
        """
        push_back(self, in_0: int ) -> None
        """
        assert isinstance(in_0, int), 'arg in_0 wrong type'
    
        self.inst.get().push_back((<int>in_0))
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        Returns the name of the peak annotations
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result
    
    def setName(self,  name ):
        """
        setName(self, name: Union[bytes, str, String] ) -> None
        Sets the name of the peak annotations
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        self.inst.get().setName(deref((convString(name)).get()))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        Returns a reference to the description of the applied processing
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
        if not isinstance(other, IntegerDataArray):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef IntegerDataArray other_casted = other
        cdef IntegerDataArray self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    


    def __getitem__(self,  in_0 ):
        assert isinstance(in_0, int), 'arg in_0 wrong type'
        assert in_0 >= 0, 'arg in_0 cannot be negative'

        cdef unsigned int _idx = (<int>in_0)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _r = deref(self.inst.get())[(<int>in_0)]
        py_result = <int>_r
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, int), 'arg key wrong type'
        assert isinstance(value, int), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef int _v = (<int>value)
        deref(self.inst.get())[(<int>key)] = _v

    def get_data(self):
        """
        Gets the raw data for the integer data array

        Example usage: 

          idata = pyopenms.IntegerDataArray()
          data = idata.get_data()

        """

        cdef _IntegerDataArray * ida_ = self.inst.get()
        cdef libcpp_vector[int].iterator it = ida_.begin()
        cdef unsigned int n = ida_.size()
        cdef np.ndarray[int, ndim=1] data
        data = np.zeros( (n,), dtype=np.intc)

        cdef int i = 0
        while it != ida_.end():
            data[i] = deref(it)
            inc(it)
            i += 1

        return data

    def set_data(self, np.ndarray[int, ndim=1, mode="c"] data not None):
        """
        Sets the raw data for the integer data array

        Example usage: 

          idata = pyopenms.IntegerDataArray()
          data = numpy.array( [1, 2, 3, 5 ,6] ).astype(np.intc)
          idata.set_data(data)

        """

        cdef _IntegerDataArray * ida_ = self.inst.get()
        ida_.clear()

        cdef int N
        N = len(data)
        ida_.reserve(N) # allocate space for incoming data
        for i in range(N):
            ida_.push_back(data[i]) 

cdef class MRMScoring:
    """
    Cython implementation of _MRMScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1MRMScoring.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MRMScoring rv = MRMScoring.__new__(MRMScoring)
       rv.inst = shared_ptr[_MRMScoring](new _MRMScoring(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MRMScoring rv = MRMScoring.__new__(MRMScoring)
       rv.inst = shared_ptr[_MRMScoring](new _MRMScoring(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MRMScoring](new _MRMScoring())
    
    def _init_1(self, MRMScoring in_0 ):
        """
        _init_1(self, in_0: MRMScoring ) -> None
        """
        assert isinstance(in_0, MRMScoring), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MRMScoring](new _MRMScoring((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MRMScoring ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MRMScoring)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def calcXcorrCoelutionScore(self):
        """
        calcXcorrCoelutionScore(self) -> float
        Calculate the cross-correlation coelution score. The score is a distance where zero indicates perfect coelution
        """
        cdef double _r = self.inst.get().calcXcorrCoelutionScore()
        py_result = <double>_r
        return py_result
    
    def calcXcorrCoelutionWeightedScore(self, list normalized_library_intensity ):
        """
        calcXcorrCoelutionWeightedScore(self, normalized_library_intensity: List[float] ) -> float
        Calculate the weighted cross-correlation coelution score
        
        The score is a distance where zero indicates perfect coelution. The
        score is weighted by the transition intensities, non-perfect coelution
        in low-intensity transitions should thus become less important
        """
        assert isinstance(normalized_library_intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in normalized_library_intensity), 'arg normalized_library_intensity wrong type'
        cdef libcpp_vector[double] v0 = normalized_library_intensity
        cdef double _r = self.inst.get().calcXcorrCoelutionWeightedScore(v0)
        normalized_library_intensity[:] = v0
        py_result = <double>_r
        return py_result
    
    def calcSeparateXcorrContrastCoelutionScore(self):
        """
        calcSeparateXcorrContrastCoelutionScore(self) -> List[float]
        Calculate the separate cross-correlation contrast score
        """
        _r = self.inst.get().calcSeparateXcorrContrastCoelutionScore()
        cdef list py_result = _r
        return py_result
    
    def calcXcorrPrecursorContrastCoelutionScore(self):
        """
        calcXcorrPrecursorContrastCoelutionScore(self) -> float
        Calculate the precursor cross-correlation contrast score against the transitions
        
        The score is a distance where zero indicates perfect coelution
        """
        cdef double _r = self.inst.get().calcXcorrPrecursorContrastCoelutionScore()
        py_result = <double>_r
        return py_result
    
    def calcXcorrShapeScore(self):
        """
        calcXcorrShapeScore(self) -> float
        Calculate the cross-correlation shape score
        
        The score is a correlation measure where 1 indicates perfect correlation
        and 0 means no correlation.
        """
        cdef double _r = self.inst.get().calcXcorrShapeScore()
        py_result = <double>_r
        return py_result
    
    def calcXcorrShapeWeightedScore(self, list normalized_library_intensity ):
        """
        calcXcorrShapeWeightedScore(self, normalized_library_intensity: List[float] ) -> float
        Calculate the weighted cross-correlation shape score
        
        The score is a correlation measure where 1 indicates perfect correlation
        and 0 means no correlation. The score is weighted by the transition
        intensities, non-perfect coelution in low-intensity transitions should
        thus become less important
        """
        assert isinstance(normalized_library_intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in normalized_library_intensity), 'arg normalized_library_intensity wrong type'
        cdef libcpp_vector[double] v0 = normalized_library_intensity
        cdef double _r = self.inst.get().calcXcorrShapeWeightedScore(v0)
        normalized_library_intensity[:] = v0
        py_result = <double>_r
        return py_result
    
    def calcSeparateXcorrContrastShapeScore(self):
        """
        calcSeparateXcorrContrastShapeScore(self) -> List[float]
        Calculate the separate cross-correlation contrast shape score
        """
        _r = self.inst.get().calcSeparateXcorrContrastShapeScore()
        cdef list py_result = _r
        return py_result
    
    def calcXcorrPrecursorContrastShapeScore(self):
        """
        calcXcorrPrecursorContrastShapeScore(self) -> float
        Calculate the precursor cross-correlation shape score against the transitions
        """
        cdef double _r = self.inst.get().calcXcorrPrecursorContrastShapeScore()
        py_result = <double>_r
        return py_result
    
    def calcRTScore(self, LightCompound peptide , double normalized_experimental_rt ):
        """
        calcRTScore(self, peptide: LightCompound , normalized_experimental_rt: float ) -> float
        """
        assert isinstance(peptide, LightCompound), 'arg peptide wrong type'
        assert isinstance(normalized_experimental_rt, float), 'arg normalized_experimental_rt wrong type'
    
    
        cdef double _r = self.inst.get().calcRTScore((deref(peptide.inst.get())), (<double>normalized_experimental_rt))
        py_result = <double>_r
        return py_result
    
    def calcMIScore(self):
        """
        calcMIScore(self) -> float
        """
        cdef double _r = self.inst.get().calcMIScore()
        py_result = <double>_r
        return py_result
    
    def calcMIWeightedScore(self, list normalized_library_intensity ):
        """
        calcMIWeightedScore(self, normalized_library_intensity: List[float] ) -> float
        """
        assert isinstance(normalized_library_intensity, list) and all(isinstance(elemt_rec, float) for elemt_rec in normalized_library_intensity), 'arg normalized_library_intensity wrong type'
        cdef libcpp_vector[double] v0 = normalized_library_intensity
        cdef double _r = self.inst.get().calcMIWeightedScore(v0)
        
        py_result = <double>_r
        return py_result
    
    def calcMIPrecursorScore(self):
        """
        calcMIPrecursorScore(self) -> float
        """
        cdef double _r = self.inst.get().calcMIPrecursorScore()
        py_result = <double>_r
        return py_result
    
    def calcMIPrecursorContrastScore(self):
        """
        calcMIPrecursorContrastScore(self) -> float
        """
        cdef double _r = self.inst.get().calcMIPrecursorContrastScore()
        py_result = <double>_r
        return py_result
    
    def calcMIPrecursorCombinedScore(self):
        """
        calcMIPrecursorCombinedScore(self) -> float
        """
        cdef double _r = self.inst.get().calcMIPrecursorCombinedScore()
        py_result = <double>_r
        return py_result
    
    def calcSeparateMIContrastScore(self):
        """
        calcSeparateMIContrastScore(self) -> List[float]
        """
        _r = self.inst.get().calcSeparateMIContrastScore()
        cdef list py_result = _r
        return py_result
    
    def getMIMatrix(self):
        """
        getMIMatrix(self) -> MatrixDouble
        """
        cdef _Matrix[double] * _r = new _Matrix[double](self.inst.get().getMIMatrix())
        cdef MatrixDouble py_result = MatrixDouble.__new__(MatrixDouble)
        py_result.inst = shared_ptr[_Matrix[double]](_r)
        return py_result 

cdef class MapAlignmentTransformer:
    """
    Cython implementation of _MapAlignmentTransformer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentTransformer.html>`_

    This class collects functions for applying retention time transformations to data structures
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MapAlignmentTransformer rv = MapAlignmentTransformer.__new__(MapAlignmentTransformer)
       rv.inst = shared_ptr[_MapAlignmentTransformer](new _MapAlignmentTransformer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MapAlignmentTransformer rv = MapAlignmentTransformer.__new__(MapAlignmentTransformer)
       rv.inst = shared_ptr[_MapAlignmentTransformer](new _MapAlignmentTransformer(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MapAlignmentTransformer](new _MapAlignmentTransformer())
    
    def _init_1(self, MapAlignmentTransformer in_0 ):
        """
        _init_1(self, in_0: MapAlignmentTransformer ) -> None
        """
        assert isinstance(in_0, MapAlignmentTransformer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MapAlignmentTransformer](new _MapAlignmentTransformer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MapAlignmentTransformer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MapAlignmentTransformer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _transformRetentionTimes_0(self, MSExperiment in_0 , TransformationDescription in_1 , bool in_2 ):
        """
        _transformRetentionTimes_0(self, in_0: MSExperiment , in_1: TransformationDescription , in_2: bool ) -> None
        Applies the given transformation to a peak map
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
        assert isinstance(in_2, pybool_t), 'arg in_2 wrong type'
    
    
    
        self.inst.get().transformRetentionTimes((deref(in_0.inst.get())), (deref(in_1.inst.get())), (<bool>in_2))
    
    def _transformRetentionTimes_1(self, FeatureMap in_0 , TransformationDescription in_1 , bool in_2 ):
        """
        _transformRetentionTimes_1(self, in_0: FeatureMap , in_1: TransformationDescription , in_2: bool ) -> None
        Applies the given transformation to a feature map
        """
        assert isinstance(in_0, FeatureMap), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
        assert isinstance(in_2, pybool_t), 'arg in_2 wrong type'
    
    
    
        self.inst.get().transformRetentionTimes((deref(in_0.inst.get())), (deref(in_1.inst.get())), (<bool>in_2))
    
    def _transformRetentionTimes_2(self, ConsensusMap in_0 , TransformationDescription in_1 , bool in_2 ):
        """
        _transformRetentionTimes_2(self, in_0: ConsensusMap , in_1: TransformationDescription , in_2: bool ) -> None
        Applies the given transformation to a consensus map
        """
        assert isinstance(in_0, ConsensusMap), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
        assert isinstance(in_2, pybool_t), 'arg in_2 wrong type'
    
    
    
        self.inst.get().transformRetentionTimes((deref(in_0.inst.get())), (deref(in_1.inst.get())), (<bool>in_2))
    
    def _transformRetentionTimes_3(self, PeptideIdentificationList in_0 , TransformationDescription in_1 , bool in_2 ):
        """
        _transformRetentionTimes_3(self, in_0: PeptideIdentificationList , in_1: TransformationDescription , in_2: bool ) -> None
        Applies the given transformation to peptide identifications
        """
        assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
        assert isinstance(in_1, TransformationDescription), 'arg in_1 wrong type'
        assert isinstance(in_2, pybool_t), 'arg in_2 wrong type'
    
    
    
        self.inst.get().transformRetentionTimes((deref(in_0.inst.get())), (deref(in_1.inst.get())), (<bool>in_2))
    
    def transformRetentionTimes(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: transformRetentionTimes(self, in_0: MSExperiment , in_1: TransformationDescription , in_2: bool ) -> None
          :noindex:
        
        Applies the given transformation to a peak map

        
        .. rubric:: Overload:
        .. py:function:: transformRetentionTimes(self, in_0: FeatureMap , in_1: TransformationDescription , in_2: bool ) -> None
          :noindex:
        
        Applies the given transformation to a feature map

        
        .. rubric:: Overload:
        .. py:function:: transformRetentionTimes(self, in_0: ConsensusMap , in_1: TransformationDescription , in_2: bool ) -> None
          :noindex:
        
        Applies the given transformation to a consensus map

        
        .. rubric:: Overload:
        .. py:function:: transformRetentionTimes(self, in_0: PeptideIdentificationList , in_1: TransformationDescription , in_2: bool ) -> None
          :noindex:
        
        Applies the given transformation to peptide identifications
    
        """
        if (len(args)==3) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], TransformationDescription)) and (isinstance(args[2], pybool_t)):
            return self._transformRetentionTimes_0(*args)
        elif (len(args)==3) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], TransformationDescription)) and (isinstance(args[2], pybool_t)):
            return self._transformRetentionTimes_1(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], TransformationDescription)) and (isinstance(args[2], pybool_t)):
            return self._transformRetentionTimes_2(*args)
        elif (len(args)==3) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], TransformationDescription)) and (isinstance(args[2], pybool_t)):
            return self._transformRetentionTimes_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class OMSSAXMLFile:
    """
    Cython implementation of _OMSSAXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OMSSAXMLFile.html>`_
      -- Inherits from ['XMLFile']

    Used to load OMSSAXML files
    
    This class is used to load documents that implement
    the schema of OMSSAXML files
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_OMSSAXMLFile](new _OMSSAXMLFile())
    
    def load(self,  filename , ProteinIdentification protein_identification , PeptideIdentificationList id_data , bool load_proteins , bool load_empty_hits ):
        """
        load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , load_proteins: bool , load_empty_hits: bool ) -> None
        Loads data from a OMSSAXML file
        
        
        :param filename: The file to be loaded
        :param protein_identification: Protein identifications belonging to the whole experiment
        :param id_data: The identifications with m/z and RT
        :param load_proteins: If this flag is set to false, the protein identifications are not loaded
        :param load_empty_hits: Many spectra will not return a hit. Report empty peptide identifications?
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_identification, ProteinIdentification), 'arg protein_identification wrong type'
        assert isinstance(id_data, PeptideIdentificationList), 'arg id_data wrong type'
        assert isinstance(load_proteins, pybool_t), 'arg load_proteins wrong type'
        assert isinstance(load_empty_hits, pybool_t), 'arg load_empty_hits wrong type'
    
    
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(protein_identification.inst.get())), (deref(id_data.inst.get())), (<bool>load_proteins), (<bool>load_empty_hits))
    
    def setModificationDefinitionsSet(self, ModificationDefinitionsSet rhs ):
        """
        setModificationDefinitionsSet(self, rhs: ModificationDefinitionsSet ) -> None
        Sets the valid modifications
        """
        assert isinstance(rhs, ModificationDefinitionsSet), 'arg rhs wrong type'
    
        self.inst.get().setModificationDefinitionsSet((deref(rhs.inst.get())))
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
        return py_result 

cdef class OSChromatogramMeta:
    """
    Cython implementation of _OSChromatogramMeta

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSChromatogramMeta.html>`_
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
    
    def __copy__(self):
       cdef OSChromatogramMeta rv = OSChromatogramMeta.__new__(OSChromatogramMeta)
       rv.inst = shared_ptr[_OSChromatogramMeta](new _OSChromatogramMeta(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OSChromatogramMeta rv = OSChromatogramMeta.__new__(OSChromatogramMeta)
       rv.inst = shared_ptr[_OSChromatogramMeta](new _OSChromatogramMeta(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OSChromatogramMeta](new _OSChromatogramMeta())
    
    def _init_1(self, OSChromatogramMeta in_0 ):
        """
        _init_1(self, in_0: OSChromatogramMeta ) -> None
        """
        assert isinstance(in_0, OSChromatogramMeta), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OSChromatogramMeta](new _OSChromatogramMeta((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OSChromatogramMeta ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OSChromatogramMeta)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeakBoundary:
    """
    Cython implementation of _PeakBoundary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakBoundary.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property mz_min:
        def __set__(self, double mz_min):
        
            self.inst.get().mz_min = (<double>mz_min)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz_min
            py_result = <double>_r
            return py_result
    
    property mz_max:
        def __set__(self, double mz_max):
        
            self.inst.get().mz_max = (<double>mz_max)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().mz_max
            py_result = <double>_r
            return py_result
    
    def __copy__(self):
       cdef PeakBoundary rv = PeakBoundary.__new__(PeakBoundary)
       rv.inst = shared_ptr[_PeakBoundary](new _PeakBoundary(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakBoundary rv = PeakBoundary.__new__(PeakBoundary)
       rv.inst = shared_ptr[_PeakBoundary](new _PeakBoundary(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakBoundary](new _PeakBoundary())
    
    def _init_1(self, PeakBoundary in_0 ):
        """
        _init_1(self, in_0: PeakBoundary ) -> None
        """
        assert isinstance(in_0, PeakBoundary), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakBoundary](new _PeakBoundary((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakBoundary ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakBoundary)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeakPickerChromatogram:
    """
    Cython implementation of _PeakPickerChromatogram

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerChromatogram.html>`_
      -- Inherits from ['DefaultParamHandler']

    The PeakPickerChromatogram finds peaks a single chromatogram
    
    It uses the PeakPickerHiRes internally to find interesting seed candidates.
    These candidates are then expanded and a right/left border of the peak is
    searched
    Additionally, overlapping peaks can be removed
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakPickerChromatogram rv = PeakPickerChromatogram.__new__(PeakPickerChromatogram)
       rv.inst = shared_ptr[_PeakPickerChromatogram](new _PeakPickerChromatogram(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakPickerChromatogram rv = PeakPickerChromatogram.__new__(PeakPickerChromatogram)
       rv.inst = shared_ptr[_PeakPickerChromatogram](new _PeakPickerChromatogram(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakPickerChromatogram](new _PeakPickerChromatogram())
    
    def _init_1(self, PeakPickerChromatogram in_0 ):
        """
        _init_1(self, in_0: PeakPickerChromatogram ) -> None
        """
        assert isinstance(in_0, PeakPickerChromatogram), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakPickerChromatogram](new _PeakPickerChromatogram((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakPickerChromatogram ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakPickerChromatogram)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def pickChromatogram(self, MSChromatogram chromatogram , MSChromatogram picked_chrom ):
        """
        pickChromatogram(self, chromatogram: MSChromatogram , picked_chrom: MSChromatogram ) -> None
        Finds peaks in a single chromatogram and annotates left/right borders
        
        It uses a modified algorithm of the PeakPickerHiRes
        
        This function will return a picked chromatogram
        """
        assert isinstance(chromatogram, MSChromatogram), 'arg chromatogram wrong type'
        assert isinstance(picked_chrom, MSChromatogram), 'arg picked_chrom wrong type'
    
    
        self.inst.get().pickChromatogram((deref(chromatogram.inst.get())), (deref(picked_chrom.inst.get())))
    
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

cdef class PeakPickerHiRes:
    """
    Cython implementation of _PeakPickerHiRes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerHiRes.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeakPickerHiRes rv = PeakPickerHiRes.__new__(PeakPickerHiRes)
       rv.inst = shared_ptr[_PeakPickerHiRes](new _PeakPickerHiRes(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeakPickerHiRes rv = PeakPickerHiRes.__new__(PeakPickerHiRes)
       rv.inst = shared_ptr[_PeakPickerHiRes](new _PeakPickerHiRes(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeakPickerHiRes](new _PeakPickerHiRes())
    
    def _init_1(self, PeakPickerHiRes in_0 ):
        """
        _init_1(self, in_0: PeakPickerHiRes ) -> None
        """
        assert isinstance(in_0, PeakPickerHiRes), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeakPickerHiRes](new _PeakPickerHiRes((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeakPickerHiRes ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeakPickerHiRes)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _pick_0(self, MSSpectrum input , MSSpectrum output ):
        """
        _pick_0(self, input: MSSpectrum , output: MSSpectrum ) -> None
        """
        assert isinstance(input, MSSpectrum), 'arg input wrong type'
        assert isinstance(output, MSSpectrum), 'arg output wrong type'
    
    
        self.inst.get().pick((deref(input.inst.get())), (deref(output.inst.get())))
    
    def _pick_1(self, MSChromatogram input , MSChromatogram output ):
        """
        _pick_1(self, input: MSChromatogram , output: MSChromatogram ) -> None
        """
        assert isinstance(input, MSChromatogram), 'arg input wrong type'
        assert isinstance(output, MSChromatogram), 'arg output wrong type'
    
    
        self.inst.get().pick((deref(input.inst.get())), (deref(output.inst.get())))
    
    def pick(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: pick(self, input: MSSpectrum , output: MSSpectrum ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: pick(self, input: MSChromatogram , output: MSChromatogram ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], MSSpectrum)):
            return self._pick_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], MSChromatogram)):
            return self._pick_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _pickExperiment_0(self, MSExperiment input , MSExperiment output , bool check_spectrum_type ):
        """
        _pickExperiment_0(self, input: MSExperiment , output: MSExperiment , check_spectrum_type: bool ) -> None
        Applies the peak-picking algorithm to a map (MSExperiment). This method picks peaks for each scan in the map consecutively. The resulting
        picked peaks are written to the output map
        
        
        :param input: Input map in profile mode
        :param output: Output map with picked peaks
        :param check_spectrum_type: If set, checks spectrum type and throws an exception if a centroided spectrum is passed
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
        assert isinstance(output, MSExperiment), 'arg output wrong type'
        assert isinstance(check_spectrum_type, pybool_t), 'arg check_spectrum_type wrong type'
    
    
    
        self.inst.get().pickExperiment((deref(input.inst.get())), (deref(output.inst.get())), (<bool>check_spectrum_type))
    
    def _pickExperiment_1(self, MSExperiment input , MSExperiment output , list boundaries_spec , list boundaries_chrom , bool check_spectrum_type ):
        """
        _pickExperiment_1(self, input: MSExperiment , output: MSExperiment , boundaries_spec: List[List[PeakBoundary]] , boundaries_chrom: List[List[PeakBoundary]] , check_spectrum_type: bool ) -> None
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
        assert isinstance(output, MSExperiment), 'arg output wrong type'
        assert isinstance(boundaries_spec, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in boundaries_spec), 'arg boundaries_spec wrong type'
        assert isinstance(boundaries_chrom, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in boundaries_chrom), 'arg boundaries_chrom wrong type'
        assert isinstance(check_spectrum_type, pybool_t), 'arg check_spectrum_type wrong type'
    
    
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]] * v2 = new libcpp_vector[libcpp_vector[_PeakBoundary]]()
        cdef libcpp_vector[_PeakBoundary] * v2_rec = new libcpp_vector[_PeakBoundary]()
        cdef PeakBoundary item2_rec
        cdef libcpp_vector[_PeakBoundary].iterator it_boundaries_spec_rec
        for boundaries_spec_rec in boundaries_spec:
            v2_rec.clear()
            for item2_rec in boundaries_spec_rec:
                v2_rec.push_back(deref(item2_rec.inst.get()))
            v2.push_back(deref(v2_rec))
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]] * v3 = new libcpp_vector[libcpp_vector[_PeakBoundary]]()
        cdef libcpp_vector[_PeakBoundary] * v3_rec = new libcpp_vector[_PeakBoundary]()
        cdef PeakBoundary item3_rec
        cdef libcpp_vector[_PeakBoundary].iterator it_boundaries_chrom_rec
        for boundaries_chrom_rec in boundaries_chrom:
            v3_rec.clear()
            for item3_rec in boundaries_chrom_rec:
                v3_rec.push_back(deref(item3_rec.inst.get()))
            v3.push_back(deref(v3_rec))
    
        self.inst.get().pickExperiment((deref(input.inst.get())), (deref(output.inst.get())), deref(v2), deref(v3), (<bool>check_spectrum_type))
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]].iterator it_boundaries_chrom = v3.begin()
        replace_0 = []
        while it_boundaries_chrom != v3.end():
            it_boundaries_chrom_rec = deref(it_boundaries_chrom).begin()
            replace_1 = []
            while it_boundaries_chrom_rec != deref(it_boundaries_chrom).end():
                item3_rec = PeakBoundary.__new__(PeakBoundary)
                item3_rec.inst = shared_ptr[_PeakBoundary](new _PeakBoundary(deref(it_boundaries_chrom_rec)))
                replace_1.append(item3_rec)
                inc(it_boundaries_chrom_rec)
            replace_0.append(replace_1)
            inc(it_boundaries_chrom)
        boundaries_chrom[:] = replace_0
        del v3
        cdef libcpp_vector[libcpp_vector[_PeakBoundary]].iterator it_boundaries_spec = v2.begin()
        replace_0 = []
        while it_boundaries_spec != v2.end():
            it_boundaries_spec_rec = deref(it_boundaries_spec).begin()
            replace_1 = []
            while it_boundaries_spec_rec != deref(it_boundaries_spec).end():
                item2_rec = PeakBoundary.__new__(PeakBoundary)
                item2_rec.inst = shared_ptr[_PeakBoundary](new _PeakBoundary(deref(it_boundaries_spec_rec)))
                replace_1.append(item2_rec)
                inc(it_boundaries_spec_rec)
            replace_0.append(replace_1)
            inc(it_boundaries_spec)
        boundaries_spec[:] = replace_0
        del v2
    
    def pickExperiment(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: pickExperiment(self, input: MSExperiment , output: MSExperiment , check_spectrum_type: bool ) -> None
          :noindex:
        
        Applies the peak-picking algorithm to a map (MSExperiment). This method picks peaks for each scan in the map consecutively. The resulting
        picked peaks are written to the output map
        
        
        :param input: Input map in profile mode
        :param output: Output map with picked peaks
        :param check_spectrum_type: If set, checks spectrum type and throws an exception if a centroided spectrum is passed
        
        .. rubric:: Overload:
        .. py:function:: pickExperiment(self, input: MSExperiment , output: MSExperiment , boundaries_spec: List[List[PeakBoundary]] , boundaries_chrom: List[List[PeakBoundary]] , check_spectrum_type: bool ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], MSExperiment)) and (isinstance(args[2], pybool_t)):
            return self._pickExperiment_0(*args)
        elif (len(args)==5) and (isinstance(args[0], MSExperiment)) and (isinstance(args[1], MSExperiment)) and (isinstance(args[2], list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in args[2])) and (isinstance(args[3], list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, PeakBoundary) for elemt_rec_rec in elemt_rec) for elemt_rec in args[3])) and (isinstance(args[4], pybool_t)):
            return self._pickExperiment_1(*args)
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

cdef class String:
    """
    Cython implementation of _String

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1String.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()


    def __hash__(self):
      # The only required property is that objects which compare equal have
      # the same hash value:
      return hash(deref(self.inst.get()).c_str() )

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_String](new _String())
    
    def __richcmp__(self, other, op):
        if op not in (2, 3):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, String):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef String other_casted = other
        cdef String self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    
    # This goes into the class
    def _init_0(self):
        self.inst = shared_ptr[_String](new _String())
    
    # TODO this confuses typing since it overwrites the autogenerated
    #  overloads from autowrap. We need to support custom pyi files
    #  for each pyx to fix this. OR maybe we need to include wrap-ignored
    #  functions for typing.
    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==1):
            self.inst = convString(args[0])
        else:
             raise Exception('can not handle type of %s' % (args,)) 

    def toString(self):
        """Cython signature: str toString()
        -- Note: this returns a unicode string and assumes the underlying
           std::string object is UTF8 encoded (ASCII is a subset
           and works, too)
        """
        # Decodes the C string to unicode
        cdef char* c_string = _cast_const_away(self.inst.get().c_str())
        cdef Py_ssize_t length = self.inst.get().length()
        ustring = c_string[:length].decode('UTF-8')
        return ustring

    def __bytes__(self):
        return self.c_str()

    def __str__(self):
        # Since python2 is deprecated we return unicode for python3.
        return str(self.toString())

    def __unicode__(self):
        return self.toString()

    def __repr__(self):
        return f"String('{self.toString()}')"

    # TODO does this really need to be callable from python? How about cdef?
    def c_str(self):
        """Cython signature: const_char * c_str()"""
        # See https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
        #    py_string = <bytes> c_string
        # This creates a Python byte string object that holds a copy of the
        # original C string. It can be safely passed around in Python code, and
        # will be garbage collected when the last reference to it goes out of
        # scope. It is important to remember that null bytes in the string act
        # as terminator character, as generally known from C. The above will
        # therefore only work correctly for C strings that do not contain null
        # bytes.
        cdef const_char  * _r = _cast_const_away(self.inst.get().c_str())
        cdef Py_ssize_t length = self.inst.get().length()
        py_result = _r[:length] # This will work correctly also if the char array contains null bytes
        return py_result 

cdef class XFDRAlgorithm:
    """
    Cython implementation of _XFDRAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XFDRAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef XFDRAlgorithm rv = XFDRAlgorithm.__new__(XFDRAlgorithm)
       rv.inst = shared_ptr[_XFDRAlgorithm](new _XFDRAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef XFDRAlgorithm rv = XFDRAlgorithm.__new__(XFDRAlgorithm)
       rv.inst = shared_ptr[_XFDRAlgorithm](new _XFDRAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_XFDRAlgorithm](new _XFDRAlgorithm())
    
    def _init_1(self, XFDRAlgorithm in_0 ):
        """
        _init_1(self, in_0: XFDRAlgorithm ) -> None
        """
        assert isinstance(in_0, XFDRAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_XFDRAlgorithm](new _XFDRAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: XFDRAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], XFDRAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, PeptideIdentificationList peptide_ids , ProteinIdentification protein_id ):
        """
        run(self, peptide_ids: PeptideIdentificationList , protein_id: ProteinIdentification ) -> int
        """
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert isinstance(protein_id, ProteinIdentification), 'arg protein_id wrong type'
    
    
        cdef _XFDRAlgorithm_ExitCodes _r = self.inst.get().run((deref(peptide_ids.inst.get())), (deref(protein_id.inst.get())))
        py_result = <int>_r
        return py_result
    
    def validateClassArguments(self):
        """
        validateClassArguments(self) -> int
        """
        cdef _XFDRAlgorithm_ExitCodes _r = self.inst.get().validateClassArguments()
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
    XFDRAlgorithm_ExitCodes = __XFDRAlgorithm_ExitCodes 

cdef class XQuestResultXMLFile:
    """
    Cython implementation of _XQuestResultXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XQuestResultXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef XQuestResultXMLFile rv = XQuestResultXMLFile.__new__(XQuestResultXMLFile)
       rv.inst = shared_ptr[_XQuestResultXMLFile](new _XQuestResultXMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef XQuestResultXMLFile rv = XQuestResultXMLFile.__new__(XQuestResultXMLFile)
       rv.inst = shared_ptr[_XQuestResultXMLFile](new _XQuestResultXMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_XQuestResultXMLFile](new _XQuestResultXMLFile())
    
    def _init_1(self, XQuestResultXMLFile in_0 ):
        """
        _init_1(self, in_0: XQuestResultXMLFile ) -> None
        """
        assert isinstance(in_0, XQuestResultXMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_XQuestResultXMLFile](new _XQuestResultXMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: XQuestResultXMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], XQuestResultXMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , PeptideIdentificationList pep_ids , list prot_ids ):
        """
        load(self, filename: Union[bytes, str, String] , pep_ids: PeptideIdentificationList , prot_ids: List[ProteinIdentification] ) -> None
        Load the content of the xquest.xml file into the provided data structures
        
        :param filename: Filename of the file which is to be loaded
        :param pep_ids: Where the spectra with identifications of the input file will be loaded to
        :param prot_ids: Where the protein identification of the input file will be loaded to
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(pep_ids, PeptideIdentificationList), 'arg pep_ids wrong type'
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
    
    
        cdef libcpp_vector[_ProteinIdentification] * v2 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item2
        for item2 in prot_ids:
            v2.push_back(deref(item2.inst.get()))
        self.inst.get().load(deref((convString(filename)).get()), (deref(pep_ids.inst.get())), deref(v2))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v2.begin()
        replace_0 = []
        while it_prot_ids != v2.end():
            item2 = ProteinIdentification.__new__(ProteinIdentification)
            item2.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item2)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v2
    
    def store(self,  filename , list poid , PeptideIdentificationList peid ):
        """
        store(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None
        Stores the identifications in a xQuest XML file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(poid, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in poid), 'arg poid wrong type'
        assert isinstance(peid, PeptideIdentificationList), 'arg peid wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in poid:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peid.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_poid = v1.begin()
        replace_0 = []
        while it_poid != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_poid)))
            replace_0.append(item1)
            inc(it_poid)
        poid[:] = replace_0
        del v1
    
    def getNumberOfHits(self):
        """
        getNumberOfHits(self) -> int
        Returns the total number of hits in the file
        """
        cdef int _r = self.inst.get().getNumberOfHits()
        py_result = <int>_r
        return py_result
    
    def getMinScore(self):
        """
        getMinScore(self) -> float
        Returns minimum score among the hits in the file
        """
        cdef double _r = self.inst.get().getMinScore()
        py_result = <double>_r
        return py_result
    
    def getMaxScore(self):
        """
        getMaxScore(self) -> float
        Returns maximum score among the hits in the file
        """
        cdef double _r = self.inst.get().getMaxScore()
        py_result = <double>_r
        return py_result
    
    def _writeXQuestXMLSpec_0(self,  out_file ,  base_name , OPXL_PreprocessedPairSpectra preprocessed_pair_spectra , list spectrum_pairs , list all_top_csms , MSExperiment spectra , bool test_mode ):
        """
        _writeXQuestXMLSpec_0(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , preprocessed_pair_spectra: OPXL_PreprocessedPairSpectra , spectrum_pairs: List[List[int, int]] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None
        Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param preprocessed_pair_spectra: The preprocessed spectra after comparing and filtering
        :param spectrum_pairs: Indices of spectrum pairs in the input map
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map
        """
        assert (isinstance(out_file, str) or isinstance(out_file, bytes) or isinstance(out_file, String)), 'arg out_file wrong type'
        assert (isinstance(base_name, str) or isinstance(base_name, bytes) or isinstance(base_name, String)), 'arg base_name wrong type'
        assert isinstance(preprocessed_pair_spectra, OPXL_PreprocessedPairSpectra), 'arg preprocessed_pair_spectra wrong type'
        assert isinstance(spectrum_pairs, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in spectrum_pairs), 'arg spectrum_pairs wrong type'
        assert isinstance(all_top_csms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in all_top_csms), 'arg all_top_csms wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(test_mode, pybool_t), 'arg test_mode wrong type'
    
    
    
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v3 = spectrum_pairs
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]] * v4 = new libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]]()
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v4_rec = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item4_rec
        for all_top_csms_rec in all_top_csms:
            v4_rec.clear()
            for item4_rec in all_top_csms_rec:
                v4_rec.push_back(deref(item4_rec.inst.get()))
            v4.push_back(deref(v4_rec))
    
    
        self.inst.get().writeXQuestXMLSpec(deref((convString(out_file)).get()), deref((convString(base_name)).get()), (deref(preprocessed_pair_spectra.inst.get())), v3, deref(v4), (deref(spectra.inst.get())), (<const bool &>test_mode))
        del v4
        del v4_rec
        
    
    def _writeXQuestXMLSpec_1(self,  out_file ,  base_name , list all_top_csms , MSExperiment spectra , bool test_mode ):
        """
        _writeXQuestXMLSpec_1(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None
        Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap
        """
        assert (isinstance(out_file, str) or isinstance(out_file, bytes) or isinstance(out_file, String)), 'arg out_file wrong type'
        assert (isinstance(base_name, str) or isinstance(base_name, bytes) or isinstance(base_name, String)), 'arg base_name wrong type'
        assert isinstance(all_top_csms, list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in all_top_csms), 'arg all_top_csms wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(test_mode, pybool_t), 'arg test_mode wrong type'
    
    
        cdef libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]] * v2 = new libcpp_vector[libcpp_vector[_CrossLinkSpectrumMatch]]()
        cdef libcpp_vector[_CrossLinkSpectrumMatch] * v2_rec = new libcpp_vector[_CrossLinkSpectrumMatch]()
        cdef CrossLinkSpectrumMatch item2_rec
        for all_top_csms_rec in all_top_csms:
            v2_rec.clear()
            for item2_rec in all_top_csms_rec:
                v2_rec.push_back(deref(item2_rec.inst.get()))
            v2.push_back(deref(v2_rec))
    
    
        self.inst.get().writeXQuestXMLSpec(deref((convString(out_file)).get()), deref((convString(base_name)).get()), deref(v2), (deref(spectra.inst.get())), (<const bool &>test_mode))
        del v2
        del v2_rec
    
    def writeXQuestXMLSpec(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: writeXQuestXMLSpec(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , preprocessed_pair_spectra: OPXL_PreprocessedPairSpectra , spectrum_pairs: List[List[int, int]] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None
          :noindex:
        
        Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param preprocessed_pair_spectra: The preprocessed spectra after comparing and filtering
        :param spectrum_pairs: Indices of spectrum pairs in the input map
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map
        
        .. rubric:: Overload:
        .. py:function:: writeXQuestXMLSpec(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None
          :noindex:
        
        Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap
    
        """
        if (len(args)==7) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], OPXL_PreprocessedPairSpectra)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in args[3])) and (isinstance(args[4], list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in args[4])) and (isinstance(args[5], MSExperiment)) and (isinstance(args[6], pybool_t)):
            return self._writeXQuestXMLSpec_0(*args)
        elif (len(args)==5) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], list) and all(isinstance(elemt_rec, list) and all(isinstance(elemt_rec_rec, CrossLinkSpectrumMatch) for elemt_rec_rec in elemt_rec) for elemt_rec in args[2])) and (isinstance(args[3], MSExperiment)) and (isinstance(args[4], pybool_t)):
            return self._writeXQuestXMLSpec_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getVersion(self):
        """
        getVersion(self) -> Union[bytes, str, String]
        Return the version of the schema
        """
        cdef _String _r = self.inst.get().getVersion()
        py_result = convOutputString(_r)
        return py_result 

cdef class XQuestScores:
    """
    Cython implementation of _XQuestScores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XQuestScores.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef XQuestScores rv = XQuestScores.__new__(XQuestScores)
       rv.inst = shared_ptr[_XQuestScores](new _XQuestScores(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef XQuestScores rv = XQuestScores.__new__(XQuestScores)
       rv.inst = shared_ptr[_XQuestScores](new _XQuestScores(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_XQuestScores](new _XQuestScores())
    
    def _init_1(self, XQuestScores in_0 ):
        """
        _init_1(self, in_0: XQuestScores ) -> None
        """
        assert isinstance(in_0, XQuestScores), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_XQuestScores](new _XQuestScores((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: XQuestScores ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], XQuestScores)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _preScore_0(self,  matched_alpha ,  ions_alpha ,  matched_beta ,  ions_beta ):
        """
        _preScore_0(self, matched_alpha: int , ions_alpha: int , matched_beta: int , ions_beta: int ) -> float
        Compute a simple and fast to compute pre-score for a cross-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
        :param matched_beta: Number of experimental peaks matched to theoretical linear ions from the beta peptide
        :param ions_beta: Number of theoretical ions from the beta peptide
        """
        assert isinstance(matched_alpha, int) and matched_alpha >= 0, 'arg matched_alpha wrong type'
        assert isinstance(ions_alpha, int) and ions_alpha >= 0, 'arg ions_alpha wrong type'
        assert isinstance(matched_beta, int) and matched_beta >= 0, 'arg matched_beta wrong type'
        assert isinstance(ions_beta, int) and ions_beta >= 0, 'arg ions_beta wrong type'
    
    
    
    
        cdef float _r = self.inst.get().preScore((<size_t>matched_alpha), (<size_t>ions_alpha), (<size_t>matched_beta), (<size_t>ions_beta))
        py_result = <float>_r
        return py_result
    
    def _preScore_1(self,  matched_alpha ,  ions_alpha ):
        """
        _preScore_1(self, matched_alpha: int , ions_alpha: int ) -> float
        Compute a simple and fast to compute pre-score for a mono-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
        """
        assert isinstance(matched_alpha, int) and matched_alpha >= 0, 'arg matched_alpha wrong type'
        assert isinstance(ions_alpha, int) and ions_alpha >= 0, 'arg ions_alpha wrong type'
    
    
        cdef float _r = self.inst.get().preScore((<size_t>matched_alpha), (<size_t>ions_alpha))
        py_result = <float>_r
        return py_result
    
    def preScore(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: preScore(self, matched_alpha: int , ions_alpha: int , matched_beta: int , ions_beta: int ) -> float
          :noindex:
        
        Compute a simple and fast to compute pre-score for a cross-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
        :param matched_beta: Number of experimental peaks matched to theoretical linear ions from the beta peptide
        :param ions_beta: Number of theoretical ions from the beta peptide
        
        .. rubric:: Overload:
        .. py:function:: preScore(self, matched_alpha: int , ions_alpha: int ) -> float
          :noindex:
        
        Compute a simple and fast to compute pre-score for a mono-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
    
        """
        if (len(args)==4) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int) and args[2] >= 0) and (isinstance(args[3], int) and args[3] >= 0):
            return self._preScore_0(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], int) and args[1] >= 0):
            return self._preScore_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def matchOddsScore(self, MSSpectrum theoretical_spec , double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , bool is_xlink_spectrum ,  n_charges ):
        """
        matchOddsScore(self, theoretical_spec: MSSpectrum , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , is_xlink_spectrum: bool , n_charges: int ) -> float
        Compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
        
        :param theoretical_spec: Theoretical spectrum, sorted by position
        :param matched_size: Alignment between the theoretical and the experimental spectra
        :param fragment_mass_tolerance: Fragment mass tolerance of the alignment
        :param fragment_mass_tolerance_unit_ppm: Fragment mass tolerance unit of the alignment, true = ppm, false = Da
        :param is_xlink_spectrum: Type of cross-link, true = cross-link, false = mono-link
        :param n_charges: Number of considered charges in the theoretical spectrum
        """
        assert isinstance(theoretical_spec, MSSpectrum), 'arg theoretical_spec wrong type'
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(is_xlink_spectrum, pybool_t), 'arg is_xlink_spectrum wrong type'
        assert isinstance(n_charges, int) and n_charges >= 0, 'arg n_charges wrong type'
    
    
    
    
    
        cdef double _r = self.inst.get().matchOddsScore((deref(theoretical_spec.inst.get())), (<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), (<bool>is_xlink_spectrum), (<size_t>n_charges))
        py_result = <double>_r
        return py_result
    
    def logOccupancyProb(self, MSSpectrum theoretical_spec ,  matched_size , double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm ):
        """
        logOccupancyProb(self, theoretical_spec: MSSpectrum , matched_size: int , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool ) -> float
        Compute the logOccupancyProb score, similar to the match_odds, a score based on the probability of getting the given number of matched peaks by chance
        
        :param theoretical_spec: Theoretical spectrum, sorted by position
        :param matched_size: Number of matched peaks between experimental and theoretical spectra
        :param fragment_mass_tolerance: The tolerance of the alignment
        :param fragment_mass_tolerance_unit: The tolerance unit of the alignment, true = ppm, false = Da
        """
        assert isinstance(theoretical_spec, MSSpectrum), 'arg theoretical_spec wrong type'
        assert isinstance(matched_size, int) and matched_size >= 0, 'arg matched_size wrong type'
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
    
    
    
    
        cdef double _r = self.inst.get().logOccupancyProb((deref(theoretical_spec.inst.get())), (<size_t>matched_size), (<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm))
        py_result = <double>_r
        return py_result
    
    def weightedTICScoreXQuest(self,  alpha_size ,  beta_size , double intsum_alpha , double intsum_beta , double total_current , bool type_is_cross_link ):
        """
        weightedTICScoreXQuest(self, alpha_size: int , beta_size: int , intsum_alpha: float , intsum_beta: float , total_current: float , type_is_cross_link: bool ) -> float
        """
        assert isinstance(alpha_size, int) and alpha_size >= 0, 'arg alpha_size wrong type'
        assert isinstance(beta_size, int) and beta_size >= 0, 'arg beta_size wrong type'
        assert isinstance(intsum_alpha, float), 'arg intsum_alpha wrong type'
        assert isinstance(intsum_beta, float), 'arg intsum_beta wrong type'
        assert isinstance(total_current, float), 'arg total_current wrong type'
        assert isinstance(type_is_cross_link, pybool_t), 'arg type_is_cross_link wrong type'
    
    
    
    
    
    
        cdef double _r = self.inst.get().weightedTICScoreXQuest((<size_t>alpha_size), (<size_t>beta_size), (<double>intsum_alpha), (<double>intsum_beta), (<double>total_current), (<bool>type_is_cross_link))
        py_result = <double>_r
        return py_result
    
    def weightedTICScore(self,  alpha_size ,  beta_size , double intsum_alpha , double intsum_beta , double total_current , bool type_is_cross_link ):
        """
        weightedTICScore(self, alpha_size: int , beta_size: int , intsum_alpha: float , intsum_beta: float , total_current: float , type_is_cross_link: bool ) -> float
        """
        assert isinstance(alpha_size, int) and alpha_size >= 0, 'arg alpha_size wrong type'
        assert isinstance(beta_size, int) and beta_size >= 0, 'arg beta_size wrong type'
        assert isinstance(intsum_alpha, float), 'arg intsum_alpha wrong type'
        assert isinstance(intsum_beta, float), 'arg intsum_beta wrong type'
        assert isinstance(total_current, float), 'arg total_current wrong type'
        assert isinstance(type_is_cross_link, pybool_t), 'arg type_is_cross_link wrong type'
    
    
    
    
    
    
        cdef double _r = self.inst.get().weightedTICScore((<size_t>alpha_size), (<size_t>beta_size), (<double>intsum_alpha), (<double>intsum_beta), (<double>total_current), (<bool>type_is_cross_link))
        py_result = <double>_r
        return py_result
    
    def matchedCurrentChain(self, list matched_spec_common , list matched_spec_xlinks , MSSpectrum spectrum_common_peaks , MSSpectrum spectrum_xlink_peaks ):
        """
        matchedCurrentChain(self, matched_spec_common: List[List[int, int]] , matched_spec_xlinks: List[List[int, int]] , spectrum_common_peaks: MSSpectrum , spectrum_xlink_peaks: MSSpectrum ) -> float
        """
        assert isinstance(matched_spec_common, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_common), 'arg matched_spec_common wrong type'
        assert isinstance(matched_spec_xlinks, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_xlinks), 'arg matched_spec_xlinks wrong type'
        assert isinstance(spectrum_common_peaks, MSSpectrum), 'arg spectrum_common_peaks wrong type'
        assert isinstance(spectrum_xlink_peaks, MSSpectrum), 'arg spectrum_xlink_peaks wrong type'
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v0 = matched_spec_common
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v1 = matched_spec_xlinks
    
    
        cdef double _r = self.inst.get().matchedCurrentChain(v0, v1, (deref(spectrum_common_peaks.inst.get())), (deref(spectrum_xlink_peaks.inst.get())))
        matched_spec_xlinks[:] = v1
        matched_spec_common[:] = v0
        py_result = <double>_r
        return py_result
    
    def totalMatchedCurrent(self, list matched_spec_common_alpha , list matched_spec_common_beta , list matched_spec_xlinks_alpha , list matched_spec_xlinks_beta , MSSpectrum spectrum_common_peaks , MSSpectrum spectrum_xlink_peaks ):
        """
        totalMatchedCurrent(self, matched_spec_common_alpha: List[List[int, int]] , matched_spec_common_beta: List[List[int, int]] , matched_spec_xlinks_alpha: List[List[int, int]] , matched_spec_xlinks_beta: List[List[int, int]] , spectrum_common_peaks: MSSpectrum , spectrum_xlink_peaks: MSSpectrum ) -> float
        """
        assert isinstance(matched_spec_common_alpha, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_common_alpha), 'arg matched_spec_common_alpha wrong type'
        assert isinstance(matched_spec_common_beta, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_common_beta), 'arg matched_spec_common_beta wrong type'
        assert isinstance(matched_spec_xlinks_alpha, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_xlinks_alpha), 'arg matched_spec_xlinks_alpha wrong type'
        assert isinstance(matched_spec_xlinks_beta, list) and all(isinstance(elemt_rec, list) and len(elemt_rec) == 2 and isinstance(elemt_rec[0], int) and elemt_rec[0] >= 0 and isinstance(elemt_rec[1], int) and elemt_rec[1] >= 0 for elemt_rec in matched_spec_xlinks_beta), 'arg matched_spec_xlinks_beta wrong type'
        assert isinstance(spectrum_common_peaks, MSSpectrum), 'arg spectrum_common_peaks wrong type'
        assert isinstance(spectrum_xlink_peaks, MSSpectrum), 'arg spectrum_xlink_peaks wrong type'
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v0 = matched_spec_common_alpha
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v1 = matched_spec_common_beta
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v2 = matched_spec_xlinks_alpha
        cdef libcpp_vector[libcpp_pair[size_t,size_t]] v3 = matched_spec_xlinks_beta
    
    
        cdef double _r = self.inst.get().totalMatchedCurrent(v0, v1, v2, v3, (deref(spectrum_common_peaks.inst.get())), (deref(spectrum_xlink_peaks.inst.get())))
        matched_spec_xlinks_beta[:] = v3
        matched_spec_xlinks_alpha[:] = v2
        matched_spec_common_beta[:] = v1
        matched_spec_common_alpha[:] = v0
        py_result = <double>_r
        return py_result
    
    def xCorrelation(self, MSSpectrum spec1 , MSSpectrum spec2 ,  maxshift , double tolerance ):
        """
        xCorrelation(self, spec1: MSSpectrum , spec2: MSSpectrum , maxshift: int , tolerance: float ) -> List[float]
        """
        assert isinstance(spec1, MSSpectrum), 'arg spec1 wrong type'
        assert isinstance(spec2, MSSpectrum), 'arg spec2 wrong type'
        assert isinstance(maxshift, int), 'arg maxshift wrong type'
        assert isinstance(tolerance, float), 'arg tolerance wrong type'
    
    
    
    
        _r = self.inst.get().xCorrelation((deref(spec1.inst.get())), (deref(spec2.inst.get())), (<int>maxshift), (<double>tolerance))
        cdef list py_result = _r
        return py_result
    
    def xCorrelationPrescore(self, MSSpectrum spec1 , MSSpectrum spec2 , double tolerance ):
        """
        xCorrelationPrescore(self, spec1: MSSpectrum , spec2: MSSpectrum , tolerance: float ) -> float
        """
        assert isinstance(spec1, MSSpectrum), 'arg spec1 wrong type'
        assert isinstance(spec2, MSSpectrum), 'arg spec2 wrong type'
        assert isinstance(tolerance, float), 'arg tolerance wrong type'
    
    
    
        cdef double _r = self.inst.get().xCorrelationPrescore((deref(spec1.inst.get())), (deref(spec2.inst.get())), (<double>tolerance))
        py_result = <double>_r
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
