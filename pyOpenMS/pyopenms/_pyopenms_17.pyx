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
def __static_File_absolutePath( file ):
    """
    __static_File_absolutePath(file: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef _String _r = _absolutePath_File(deref((convString(file)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_basename( file ):
    """
    __static_File_basename(file: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef _String _r = _basename_File(deref((convString(file)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_empty( file ):
    """
    __static_File_empty(file: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef bool _r = _empty_File(deref((convString(file)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_exists( file ):
    """
    __static_File_exists(file: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef bool _r = _exists_File(deref((convString(file)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_fileList( dir ,  file_pattern , list output , bool full_path ):
    """
    __static_File_fileList(dir: Union[bytes, str, String] , file_pattern: Union[bytes, str, String] , output: List[bytes] , full_path: bool ) -> bool
    """
    assert (isinstance(dir, str) or isinstance(dir, bytes) or isinstance(dir, String)), 'arg dir wrong type'
    assert (isinstance(file_pattern, str) or isinstance(file_pattern, bytes) or isinstance(file_pattern, String)), 'arg file_pattern wrong type'
    assert isinstance(output, list) and all(isinstance(li, bytes) for li in output), 'arg output wrong type'
    assert isinstance(full_path, pybool_t), 'arg full_path wrong type'


    cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
    cdef bytes item2
    for item2 in output:
       v2.push_back(_String(<char *>item2))

    cdef bool _r = _fileList_File(deref((convString(dir)).get()), deref((convString(file_pattern)).get()), deref(v2), (<bool>full_path))
    del v2
    py_result = <bool>_r
    return py_result

def __static_File_find( filename , list directories ):
    """
    __static_File_find(filename: Union[bytes, str, String] , directories: List[bytes] ) -> Union[bytes, str, String]
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert isinstance(directories, list) and all(isinstance(li, bytes) for li in directories), 'arg directories wrong type'

    cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
    cdef bytes item1
    for item1 in directories:
       v1.push_back(_String(<char *>item1))
    cdef _String _r = _find_File(deref((convString(filename)).get()), deref(v1))
    del v1
    py_result = convOutputString(_r)
    return py_result

def __static_File_findDatabase( db_name ):
    """
    __static_File_findDatabase(db_name: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(db_name, str) or isinstance(db_name, bytes) or isinstance(db_name, String)), 'arg db_name wrong type'

    cdef _String _r = _findDatabase_File(deref((convString(db_name)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_findDoc( filename ):
    """
    __static_File_findDoc(filename: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'

    cdef _String _r = _findDoc_File(deref((convString(filename)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_findExecutable( toolName ):
    """
    __static_File_findExecutable(toolName: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(toolName, str) or isinstance(toolName, bytes) or isinstance(toolName, String)), 'arg toolName wrong type'

    cdef _String _r = _findExecutable_File(deref((convString(toolName)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_getExecutablePath():
    """
    __static_File_getExecutablePath() -> Union[bytes, str, String]
    """
    cdef _String _r = _getExecutablePath_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_getOpenMSDataPath():
    """
    __static_File_getOpenMSDataPath() -> Union[bytes, str, String]
    """
    cdef _String _r = _getOpenMSDataPath_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_getOpenMSHomePath():
    """
    __static_File_getOpenMSHomePath() -> Union[bytes, str, String]
    """
    cdef _String _r = _getOpenMSHomePath_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_getSystemParameters():
    """
    __static_File_getSystemParameters() -> Param
    """
    cdef _Param * _r = new _Param(_getSystemParameters_File())
    cdef Param py_result = Param.__new__(Param)
    py_result.inst = shared_ptr[_Param](_r)
    return py_result

def __static_File_getTempDirectory():
    """
    __static_File_getTempDirectory() -> Union[bytes, str, String]
    """
    cdef _String _r = _getTempDirectory_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_getTemporaryFile( alternative_file ):
    """
    __static_File_getTemporaryFile(alternative_file: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(alternative_file, str) or isinstance(alternative_file, bytes) or isinstance(alternative_file, String)), 'arg alternative_file wrong type'

    cdef _String _r = _getTemporaryFile_File(deref((convString(alternative_file)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_getUniqueName():
    """
    __static_File_getUniqueName() -> Union[bytes, str, String]
    """
    cdef _String _r = _getUniqueName_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_getUserDirectory():
    """
    __static_File_getUserDirectory() -> Union[bytes, str, String]
    """
    cdef _String _r = _getUserDirectory_File()
    py_result = convOutputString(_r)
    return py_result

def __static_File_isDirectory( path ):
    """
    __static_File_isDirectory(path: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(path, str) or isinstance(path, bytes) or isinstance(path, String)), 'arg path wrong type'

    cdef bool _r = _isDirectory_File(deref((convString(path)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_path( file ):
    """
    __static_File_path(file: Union[bytes, str, String] ) -> Union[bytes, str, String]
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef _String _r = _path_File(deref((convString(file)).get()))
    py_result = convOutputString(_r)
    return py_result

def __static_File_readable( file ):
    """
    __static_File_readable(file: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef bool _r = _readable_File(deref((convString(file)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_remove( file ):
    """
    __static_File_remove(file: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef bool _r = _remove_File(deref((convString(file)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_removeDirRecursively( dir_name ):
    """
    __static_File_removeDirRecursively(dir_name: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(dir_name, str) or isinstance(dir_name, bytes) or isinstance(dir_name, String)), 'arg dir_name wrong type'

    cdef bool _r = _removeDirRecursively_File(deref((convString(dir_name)).get()))
    py_result = <bool>_r
    return py_result

def __static_File_rename( old_filename ,  new_filename , bool overwrite_existing , bool verbose ):
    """
    __static_File_rename(old_filename: Union[bytes, str, String] , new_filename: Union[bytes, str, String] , overwrite_existing: bool , verbose: bool ) -> bool
    """
    assert (isinstance(old_filename, str) or isinstance(old_filename, bytes) or isinstance(old_filename, String)), 'arg old_filename wrong type'
    assert (isinstance(new_filename, str) or isinstance(new_filename, bytes) or isinstance(new_filename, String)), 'arg new_filename wrong type'
    assert isinstance(overwrite_existing, pybool_t), 'arg overwrite_existing wrong type'
    assert isinstance(verbose, pybool_t), 'arg verbose wrong type'




    cdef bool _r = _rename_File(deref((convString(old_filename)).get()), deref((convString(new_filename)).get()), (<bool>overwrite_existing), (<bool>verbose))
    py_result = <bool>_r
    return py_result

def __static_File_writable( file ):
    """
    __static_File_writable(file: Union[bytes, str, String] ) -> bool
    """
    assert (isinstance(file, str) or isinstance(file, bytes) or isinstance(file, String)), 'arg file wrong type'

    cdef bool _r = _writable_File(deref((convString(file)).get()))
    py_result = <bool>_r
    return py_result 

cdef class __Sorted:
    None
    INTENSITY = 0
    MASS = 1
    UNDEFINED = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class ClusteringGrid:
    """
    Cython implementation of _ClusteringGrid

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ClusteringGrid.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ClusteringGrid rv = ClusteringGrid.__new__(ClusteringGrid)
       rv.inst = shared_ptr[_ClusteringGrid](new _ClusteringGrid(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ClusteringGrid rv = ClusteringGrid.__new__(ClusteringGrid)
       rv.inst = shared_ptr[_ClusteringGrid](new _ClusteringGrid(deref(self.inst.get())))
       return rv
    
    def _init_0(self, list grid_spacing_x , list grid_spacing_y ):
        """
        _init_0(self, grid_spacing_x: List[float] , grid_spacing_y: List[float] ) -> None
        """
        assert isinstance(grid_spacing_x, list) and all(isinstance(elemt_rec, float) for elemt_rec in grid_spacing_x), 'arg grid_spacing_x wrong type'
        assert isinstance(grid_spacing_y, list) and all(isinstance(elemt_rec, float) for elemt_rec in grid_spacing_y), 'arg grid_spacing_y wrong type'
        cdef libcpp_vector[double] v0 = grid_spacing_x
        cdef libcpp_vector[double] v1 = grid_spacing_y
        self.inst = shared_ptr[_ClusteringGrid](new _ClusteringGrid(v0, v1))
        grid_spacing_y[:] = v1
        grid_spacing_x[:] = v0
    
    def _init_1(self, ClusteringGrid in_0 ):
        """
        _init_1(self, in_0: ClusteringGrid ) -> None
        """
        assert isinstance(in_0, ClusteringGrid), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ClusteringGrid](new _ClusteringGrid((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, grid_spacing_x: List[float] , grid_spacing_y: List[float] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ClusteringGrid ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ClusteringGrid)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getGridSpacingX(self):
        """
        getGridSpacingX(self) -> List[float]
        """
        _r = self.inst.get().getGridSpacingX()
        cdef list py_result = _r
        return py_result
    
    def getGridSpacingY(self):
        """
        getGridSpacingY(self) -> List[float]
        """
        _r = self.inst.get().getGridSpacingY()
        cdef list py_result = _r
        return py_result
    
    def addCluster(self, list cell_index ,  cluster_index ):
        """
        addCluster(self, cell_index: List[int, int] , cluster_index: int ) -> None
        Adds a cluster to this grid cell
        """
        assert isinstance(cell_index, list) and len(cell_index) == 2 and isinstance(cell_index[0], int) and isinstance(cell_index[1], int), 'arg cell_index wrong type'
        assert isinstance(cluster_index, int), 'arg cluster_index wrong type'
        cdef libcpp_pair[int, int] v0
        v0.first = cell_index[0]
        v0.second = cell_index[1]
    
        self.inst.get().addCluster(v0, (<int &>cluster_index))
    
    def removeCluster(self, list cell_index ,  cluster_index ):
        """
        removeCluster(self, cell_index: List[int, int] , cluster_index: int ) -> None
        Removes a cluster from this grid cell and removes the cell if no other cluster left
        """
        assert isinstance(cell_index, list) and len(cell_index) == 2 and isinstance(cell_index[0], int) and isinstance(cell_index[1], int), 'arg cell_index wrong type'
        assert isinstance(cluster_index, int), 'arg cluster_index wrong type'
        cdef libcpp_pair[int, int] v0
        v0.first = cell_index[0]
        v0.second = cell_index[1]
    
        self.inst.get().removeCluster(v0, (<int &>cluster_index))
    
    def removeAllClusters(self):
        """
        removeAllClusters(self) -> None
        Removes all clusters from this grid (and hence all cells)
        """
        self.inst.get().removeAllClusters()
    
    def getIndex(self,  position ):
        """
        getIndex(self, position: Union[Sequence[int], Sequence[float]] ) -> List[int, int]
        """
        assert len(position) == 2 and isinstance(position[0], (int, float)) and isinstance(position[1], (int, float)), 'arg position wrong type'
        cdef _DPosition2 _dp_0
        _dp_0[0] = <float>position[0]
        _dp_0[1] = <float>position[1]
        _r = self.inst.get().getIndex(_dp_0)
        cdef list py_result = [_r.first, _r.second]
        return py_result
    
    def isNonEmptyCell(self, list cell_index ):
        """
        isNonEmptyCell(self, cell_index: List[int, int] ) -> bool
        Checks if there are clusters at this cell index
        """
        assert isinstance(cell_index, list) and len(cell_index) == 2 and isinstance(cell_index[0], int) and isinstance(cell_index[1], int), 'arg cell_index wrong type'
        cdef libcpp_pair[int, int] v0
        v0.first = cell_index[0]
        v0.second = cell_index[1]
        cdef bool _r = self.inst.get().isNonEmptyCell(v0)
        py_result = <bool>_r
        return py_result
    
    def getCellCount(self):
        """
        getCellCount(self) -> int
        Returns number of grid cells occupied by one or more clusters
        """
        cdef int _r = self.inst.get().getCellCount()
        py_result = <int>_r
        return py_result 

cdef class CoarseIsotopePatternGenerator:
    """
    Cython implementation of _CoarseIsotopePatternGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CoarseIsotopePatternGenerator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_CoarseIsotopePatternGenerator](new _CoarseIsotopePatternGenerator())
    
    def _init_1(self,  max_isotope ):
        """
        _init_1(self, max_isotope: int ) -> None
        """
        assert isinstance(max_isotope, int) and max_isotope >= 0, 'arg max_isotope wrong type'
    
        self.inst = shared_ptr[_CoarseIsotopePatternGenerator](new _CoarseIsotopePatternGenerator((<size_t>max_isotope)))
    
    def _init_2(self,  max_isotope , bool round_masses ):
        """
        _init_2(self, max_isotope: int , round_masses: bool ) -> None
        """
        assert isinstance(max_isotope, int) and max_isotope >= 0, 'arg max_isotope wrong type'
        assert isinstance(round_masses, pybool_t), 'arg round_masses wrong type'
    
    
        self.inst = shared_ptr[_CoarseIsotopePatternGenerator](new _CoarseIsotopePatternGenerator((<size_t>max_isotope), (<bool>round_masses)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, max_isotope: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, max_isotope: int , round_masses: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], int) and args[0] >= 0):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], pybool_t)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def run(self, EmpiricalFormula in_0 ):
        """
        run(self, in_0: EmpiricalFormula ) -> IsotopeDistribution
        """
        assert isinstance(in_0, EmpiricalFormula), 'arg in_0 wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().run((deref(in_0.inst.get()))))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def getRoundMasses(self):
        """
        getRoundMasses(self) -> bool
        Returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)
        """
        cdef bool _r = self.inst.get().getRoundMasses()
        py_result = <bool>_r
        return py_result
    
    def setRoundMasses(self, bool round_masses_ ):
        """
        setRoundMasses(self, round_masses_: bool ) -> None
        Sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)
        """
        assert isinstance(round_masses_, pybool_t), 'arg round_masses_ wrong type'
    
        self.inst.get().setRoundMasses((<bool>round_masses_))
    
    def getMaxIsotope(self):
        """
        getMaxIsotope(self) -> int
        Returns the currently set maximum isotope
        """
        cdef size_t _r = self.inst.get().getMaxIsotope()
        py_result = <size_t>_r
        return py_result
    
    def setMaxIsotope(self,  max_isotope ):
        """
        setMaxIsotope(self, max_isotope: int ) -> None
        Sets the maximal isotope with 'max_isotope'
        """
        assert isinstance(max_isotope, int) and max_isotope >= 0, 'arg max_isotope wrong type'
    
        self.inst.get().setMaxIsotope((<size_t>max_isotope))
    
    def estimateFromPeptideWeight(self, double average_weight ):
        """
        estimateFromPeptideWeight(self, average_weight: float ) -> IsotopeDistribution
        Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromPeptideWeight((<double>average_weight)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateFromPeptideWeightAndS(self, double average_weight ,  S ):
        """
        estimateFromPeptideWeightAndS(self, average_weight: float , S: int ) -> IsotopeDistribution
        Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(S, int), 'arg S wrong type'
    
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromPeptideWeightAndS((<double>average_weight), (<unsigned int>S)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateFromRNAWeight(self, double average_weight ):
        """
        estimateFromRNAWeight(self, average_weight: float ) -> IsotopeDistribution
        Estimate Nucleotide Isotopedistribution from weight
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromRNAWeight((<double>average_weight)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateFromDNAWeight(self, double average_weight ):
        """
        estimateFromDNAWeight(self, average_weight: float ) -> IsotopeDistribution
        Estimate Nucleotide Isotopedistribution from weight
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromDNAWeight((<double>average_weight)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateFromWeightAndComp(self, double average_weight , double C , double H , double N , double O , double S , double P ):
        """
        estimateFromWeightAndComp(self, average_weight: float , C: float , H: float , N: float , O: float , S: float , P: float ) -> IsotopeDistribution
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(C, float), 'arg C wrong type'
        assert isinstance(H, float), 'arg H wrong type'
        assert isinstance(N, float), 'arg N wrong type'
        assert isinstance(O, float), 'arg O wrong type'
        assert isinstance(S, float), 'arg S wrong type'
        assert isinstance(P, float), 'arg P wrong type'
    
    
    
    
    
    
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromWeightAndComp((<double>average_weight), (<double>C), (<double>H), (<double>N), (<double>O), (<double>S), (<double>P)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateFromWeightAndCompAndS(self, double average_weight ,  S , double C , double H , double N , double O , double P ):
        """
        estimateFromWeightAndCompAndS(self, average_weight: float , S: int , C: float , H: float , N: float , O: float , P: float ) -> IsotopeDistribution
        Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition
        """
        assert isinstance(average_weight, float), 'arg average_weight wrong type'
        assert isinstance(S, int), 'arg S wrong type'
        assert isinstance(C, float), 'arg C wrong type'
        assert isinstance(H, float), 'arg H wrong type'
        assert isinstance(N, float), 'arg N wrong type'
        assert isinstance(O, float), 'arg O wrong type'
        assert isinstance(P, float), 'arg P wrong type'
    
    
    
    
    
    
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateFromWeightAndCompAndS((<double>average_weight), (<unsigned int>S), (<double>C), (<double>H), (<double>N), (<double>O), (<double>P)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateForFragmentFromPeptideWeight(self, double average_weight_precursor , double average_weight_fragment , set precursor_isotopes ):
        """
        estimateForFragmentFromPeptideWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution
        Estimate peptide fragment IsotopeDistribution from the precursor's average weight, fragment's average weight, and a set of isolated precursor isotopes
        """
        assert isinstance(average_weight_precursor, float), 'arg average_weight_precursor wrong type'
        assert isinstance(average_weight_fragment, float), 'arg average_weight_fragment wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
    
    
        cdef libcpp_set[unsigned int] v2 = precursor_isotopes
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateForFragmentFromPeptideWeight((<double>average_weight_precursor), (<double>average_weight_fragment), v2))
        precursor_isotopes.clear()
        precursor_isotopes.update(v2)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateForFragmentFromPeptideWeightAndS(self, double average_weight_precursor ,  S_precursor , double average_weight_fragment ,  S_fragment , set precursor_isotopes ):
        """
        estimateForFragmentFromPeptideWeightAndS(self, average_weight_precursor: float , S_precursor: int , average_weight_fragment: float , S_fragment: int , precursor_isotopes: Set[int] ) -> IsotopeDistribution
        Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
        and a set of isolated precursor isotopes.
        """
        assert isinstance(average_weight_precursor, float), 'arg average_weight_precursor wrong type'
        assert isinstance(S_precursor, int), 'arg S_precursor wrong type'
        assert isinstance(average_weight_fragment, float), 'arg average_weight_fragment wrong type'
        assert isinstance(S_fragment, int), 'arg S_fragment wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
    
    
    
    
        cdef libcpp_set[unsigned int] v4 = precursor_isotopes
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateForFragmentFromPeptideWeightAndS((<double>average_weight_precursor), (<unsigned int>S_precursor), (<double>average_weight_fragment), (<unsigned int>S_fragment), v4))
        precursor_isotopes.clear()
        precursor_isotopes.update(v4)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def approximateFromPeptideWeight(self, double mass ,  num_peaks ,  charge ):
        """
        approximateFromPeptideWeight(self, mass: float , num_peaks: int , charge: int ) -> IsotopeDistribution
        Roughly approximate peptide IsotopeDistribution from monoisotopic weight using Poisson distribution.
        m/z values approximated by adding one neutron mass (divided by charge) for every peak, starting at
        the given monoisotopic weight. Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
        This method is around 50 times faster than estimateFromPeptideWeight, but only an approximation.
        The following are the intensities of the first 6 peaks generated for a monoisotopic mass of 1000:
        estimateFromPeptideWeight:    0.571133000;0.306181000;0.095811100;0.022036900;0.004092170;0.000644568
        approximateFromPeptideWeight: 0.573753000;0.318752000;0.088542200;0.016396700;0.002277320;0.000253036
        KL divergences of the first 20 intensities of estimateFromPeptideWeight and this approximation range from 4.97E-5 for a
        monoisotopic mass of 20 to 0.0144 for a mass of 2500. For comparison, when comparing an observed pattern with a
        theoretical ground truth, the observed pattern is said to be an isotopic pattern if the KL between the two is below 0.05
        for 2 peaks and below 0.6 for >=6 peaks by Guo Ci Teo et al.
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(num_peaks, int), 'arg num_peaks wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
    
    
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().approximateFromPeptideWeight((<double>mass), (<unsigned int>num_peaks), (<unsigned int>charge)))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def approximateIntensities(self, double mass ,  num_peaks ):
        """
        approximateIntensities(self, mass: float , num_peaks: int ) -> List[float]
        Roughly approximate peptidic isotope pattern intensities from monoisotopic weight using Poisson distribution.
        Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
        This method is around 100 times faster than estimateFromPeptideWeight, but only an approximation, see approximateFromPeptideWeight.
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(num_peaks, int), 'arg num_peaks wrong type'
    
    
        _r = self.inst.get().approximateIntensities((<double>mass), (<unsigned int>num_peaks))
        cdef list py_result = _r
        return py_result
    
    def estimateForFragmentFromRNAWeight(self, double average_weight_precursor , double average_weight_fragment , set precursor_isotopes ):
        """
        estimateForFragmentFromRNAWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution
        Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, and a set of isolated precursor isotopes
        """
        assert isinstance(average_weight_precursor, float), 'arg average_weight_precursor wrong type'
        assert isinstance(average_weight_fragment, float), 'arg average_weight_fragment wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
    
    
        cdef libcpp_set[unsigned int] v2 = precursor_isotopes
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateForFragmentFromRNAWeight((<double>average_weight_precursor), (<double>average_weight_fragment), v2))
        precursor_isotopes.clear()
        precursor_isotopes.update(v2)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateForFragmentFromDNAWeight(self, double average_weight_precursor , double average_weight_fragment , set precursor_isotopes ):
        """
        estimateForFragmentFromDNAWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution
        Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, and a set of isolated precursor isotopes.
        """
        assert isinstance(average_weight_precursor, float), 'arg average_weight_precursor wrong type'
        assert isinstance(average_weight_fragment, float), 'arg average_weight_fragment wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
    
    
        cdef libcpp_set[unsigned int] v2 = precursor_isotopes
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateForFragmentFromDNAWeight((<double>average_weight_precursor), (<double>average_weight_fragment), v2))
        precursor_isotopes.clear()
        precursor_isotopes.update(v2)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def estimateForFragmentFromWeightAndComp(self, double average_weight_precursor , double average_weight_fragment , set precursor_isotopes , double C , double H , double N , double O , double S , double P ):
        """
        estimateForFragmentFromWeightAndComp(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] , C: float , H: float , N: float , O: float , S: float , P: float ) -> IsotopeDistribution
        Estimate fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, a set of isolated precursor isotopes, and average composition
        """
        assert isinstance(average_weight_precursor, float), 'arg average_weight_precursor wrong type'
        assert isinstance(average_weight_fragment, float), 'arg average_weight_fragment wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
        assert isinstance(C, float), 'arg C wrong type'
        assert isinstance(H, float), 'arg H wrong type'
        assert isinstance(N, float), 'arg N wrong type'
        assert isinstance(O, float), 'arg O wrong type'
        assert isinstance(S, float), 'arg S wrong type'
        assert isinstance(P, float), 'arg P wrong type'
    
    
        cdef libcpp_set[unsigned int] v2 = precursor_isotopes
    
    
    
    
    
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().estimateForFragmentFromWeightAndComp((<double>average_weight_precursor), (<double>average_weight_fragment), v2, (<double>C), (<double>H), (<double>N), (<double>O), (<double>S), (<double>P)))
        precursor_isotopes.clear()
        precursor_isotopes.update(v2)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result
    
    def calcFragmentIsotopeDist(self, IsotopeDistribution fragment_isotope_dist , IsotopeDistribution comp_fragment_isotope_dist , set precursor_isotopes , double fragment_mono_mass ):
        """
        calcFragmentIsotopeDist(self, fragment_isotope_dist: IsotopeDistribution , comp_fragment_isotope_dist: IsotopeDistribution , precursor_isotopes: Set[int] , fragment_mono_mass: float ) -> IsotopeDistribution
        Calculate isotopic distribution for a fragment molecule
        """
        assert isinstance(fragment_isotope_dist, IsotopeDistribution), 'arg fragment_isotope_dist wrong type'
        assert isinstance(comp_fragment_isotope_dist, IsotopeDistribution), 'arg comp_fragment_isotope_dist wrong type'
        assert isinstance(precursor_isotopes, set) and all(isinstance(li, int) for li in precursor_isotopes), 'arg precursor_isotopes wrong type'
        assert isinstance(fragment_mono_mass, float), 'arg fragment_mono_mass wrong type'
    
    
        cdef libcpp_set[unsigned int] v2 = precursor_isotopes
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().calcFragmentIsotopeDist((deref(fragment_isotope_dist.inst.get())), (deref(comp_fragment_isotope_dist.inst.get())), v2, (<double>fragment_mono_mass)))
        precursor_isotopes.clear()
        precursor_isotopes.update(v2)
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result 

cdef class ConsensusMapNormalizerAlgorithmThreshold:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmThreshold

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusMapNormalizerAlgorithmThreshold.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ConsensusMapNormalizerAlgorithmThreshold](new _ConsensusMapNormalizerAlgorithmThreshold())
    
    def computeCorrelation(self, ConsensusMap input_map , double ratio_threshold ,  acc_filter ,  desc_filter ):
        """
        computeCorrelation(self, input_map: ConsensusMap , ratio_threshold: float , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> List[float]
        Determines the ratio of all maps to the map with the most features
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
        assert isinstance(ratio_threshold, float), 'arg ratio_threshold wrong type'
        assert (isinstance(acc_filter, str) or isinstance(acc_filter, bytes) or isinstance(acc_filter, String)), 'arg acc_filter wrong type'
        assert (isinstance(desc_filter, str) or isinstance(desc_filter, bytes) or isinstance(desc_filter, String)), 'arg desc_filter wrong type'
    
    
    
    
        _r = self.inst.get().computeCorrelation((deref(input_map.inst.get())), (<double>ratio_threshold), deref((convString(acc_filter)).get()), deref((convString(desc_filter)).get()))
        cdef list py_result = _r
        return py_result
    
    def normalizeMaps(self, ConsensusMap input_map , list ratios ):
        """
        normalizeMaps(self, input_map: ConsensusMap , ratios: List[float] ) -> None
        Applies the given ratio to the maps of the consensusMap
        """
        assert isinstance(input_map, ConsensusMap), 'arg input_map wrong type'
        assert isinstance(ratios, list) and all(isinstance(elemt_rec, float) for elemt_rec in ratios), 'arg ratios wrong type'
    
        cdef libcpp_vector[double] v1 = ratios
        self.inst.get().normalizeMaps((deref(input_map.inst.get())), v1)
        ratios[:] = v1 

cdef class DocumentIdentifier:
    """
    Cython implementation of _DocumentIdentifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DocumentIdentifier.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DocumentIdentifier rv = DocumentIdentifier.__new__(DocumentIdentifier)
       rv.inst = shared_ptr[_DocumentIdentifier](new _DocumentIdentifier(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DocumentIdentifier rv = DocumentIdentifier.__new__(DocumentIdentifier)
       rv.inst = shared_ptr[_DocumentIdentifier](new _DocumentIdentifier(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DocumentIdentifier](new _DocumentIdentifier())
    
    def _init_1(self, DocumentIdentifier in_0 ):
        """
        _init_1(self, in_0: DocumentIdentifier ) -> None
        """
        assert isinstance(in_0, DocumentIdentifier), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DocumentIdentifier](new _DocumentIdentifier((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DocumentIdentifier ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DocumentIdentifier)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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

cdef class File:
    """
    Cython implementation of _File

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1File.html>`_
    """

    absolutePath = __static_File_absolutePath
    basename = __static_File_basename
    empty = __static_File_empty
    exists = __static_File_exists
    fileList = __static_File_fileList
    find = __static_File_find
    findDatabase = __static_File_findDatabase
    findDoc = __static_File_findDoc
    findExecutable = __static_File_findExecutable
    getExecutablePath = __static_File_getExecutablePath
    getOpenMSDataPath = __static_File_getOpenMSDataPath
    getOpenMSHomePath = __static_File_getOpenMSHomePath
    getSystemParameters = __static_File_getSystemParameters
    getTempDirectory = __static_File_getTempDirectory
    getTemporaryFile = __static_File_getTemporaryFile
    getUniqueName = __static_File_getUniqueName
    getUserDirectory = __static_File_getUserDirectory
    isDirectory = __static_File_isDirectory
    path = __static_File_path
    readable = __static_File_readable
    remove = __static_File_remove
    removeDirRecursively = __static_File_removeDirRecursively
    rename = __static_File_rename
    writable = __static_File_writable 

cdef class FineIsotopePatternGenerator:
    """
    Cython implementation of _FineIsotopePatternGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FineIsotopePatternGenerator.html>`_

    Isotope pattern generator for fine isotope distributions.
    Generates isotopes until a stop condition (threshold) is reached,
    the lower the threshold the more isotopes are generated. The
    parameter use_total_prob defines whether the stop condition is
    interpreted as the total probability that the distribution should
    cover (default) or as a threshold for individual peaks. Finally,
    the absolute parameter specifies for individual peak thresholding
    if the threshold is absolute or relative.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_FineIsotopePatternGenerator](new _FineIsotopePatternGenerator())
    
    def _init_1(self, double threshold ):
        """
        _init_1(self, threshold: float ) -> None
        """
        assert isinstance(threshold, float), 'arg threshold wrong type'
    
        self.inst = shared_ptr[_FineIsotopePatternGenerator](new _FineIsotopePatternGenerator((<double>threshold)))
    
    def _init_2(self, double threshold , bool use_total_prob ):
        """
        _init_2(self, threshold: float , use_total_prob: bool ) -> None
        """
        assert isinstance(threshold, float), 'arg threshold wrong type'
        assert isinstance(use_total_prob, pybool_t), 'arg use_total_prob wrong type'
    
    
        self.inst = shared_ptr[_FineIsotopePatternGenerator](new _FineIsotopePatternGenerator((<double>threshold), (<bool>use_total_prob)))
    
    def _init_3(self, double threshold , bool use_total_prob , bool absolute ):
        """
        _init_3(self, threshold: float , use_total_prob: bool , absolute: bool ) -> None
        """
        assert isinstance(threshold, float), 'arg threshold wrong type'
        assert isinstance(use_total_prob, pybool_t), 'arg use_total_prob wrong type'
        assert isinstance(absolute, pybool_t), 'arg absolute wrong type'
    
    
    
        self.inst = shared_ptr[_FineIsotopePatternGenerator](new _FineIsotopePatternGenerator((<double>threshold), (<bool>use_total_prob), (<bool>absolute)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, threshold: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, threshold: float , use_total_prob: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, threshold: float , use_total_prob: bool , absolute: bool ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], float)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], float)) and (isinstance(args[1], pybool_t)):
             self._init_2(*args)
        elif (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], pybool_t)):
             self._init_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setThreshold(self, double threshold ):
        """
        setThreshold(self, threshold: float ) -> None
        """
        assert isinstance(threshold, float), 'arg threshold wrong type'
    
        self.inst.get().setThreshold((<double>threshold))
    
    def getThreshold(self):
        """
        getThreshold(self) -> float
        """
        cdef double _r = self.inst.get().getThreshold()
        py_result = <double>_r
        return py_result
    
    def setAbsolute(self, bool absolute ):
        """
        setAbsolute(self, absolute: bool ) -> None
        """
        assert isinstance(absolute, pybool_t), 'arg absolute wrong type'
    
        self.inst.get().setAbsolute((<bool>absolute))
    
    def getAbsolute(self):
        """
        getAbsolute(self) -> bool
        """
        cdef bool _r = self.inst.get().getAbsolute()
        py_result = <bool>_r
        return py_result
    
    def setTotalProbability(self, bool total ):
        """
        setTotalProbability(self, total: bool ) -> None
        """
        assert isinstance(total, pybool_t), 'arg total wrong type'
    
        self.inst.get().setTotalProbability((<bool>total))
    
    def getTotalProbability(self):
        """
        getTotalProbability(self) -> bool
        """
        cdef bool _r = self.inst.get().getTotalProbability()
        py_result = <bool>_r
        return py_result
    
    def run(self, EmpiricalFormula in_0 ):
        """
        run(self, in_0: EmpiricalFormula ) -> IsotopeDistribution
        """
        assert isinstance(in_0, EmpiricalFormula), 'arg in_0 wrong type'
    
        cdef _IsotopeDistribution * _r = new _IsotopeDistribution(self.inst.get().run((deref(in_0.inst.get()))))
        cdef IsotopeDistribution py_result = IsotopeDistribution.__new__(IsotopeDistribution)
        py_result.inst = shared_ptr[_IsotopeDistribution](_r)
        return py_result 

cdef class Gradient:
    """
    Cython implementation of _Gradient

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Gradient.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef Gradient rv = Gradient.__new__(Gradient)
       rv.inst = shared_ptr[_Gradient](new _Gradient(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef Gradient rv = Gradient.__new__(Gradient)
       rv.inst = shared_ptr[_Gradient](new _Gradient(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Representation of a HPLC gradient
        """
        self.inst = shared_ptr[_Gradient](new _Gradient())
    
    def _init_1(self, Gradient in_0 ):
        """
        _init_1(self, in_0: Gradient ) -> None
        """
        assert isinstance(in_0, Gradient), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Gradient](new _Gradient((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Representation of a HPLC gradient

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Gradient ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Gradient)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def addEluent(self,  eluent ):
        """
        addEluent(self, eluent: Union[bytes, str, String] ) -> None
        Adds an eluent at the end of the eluent array
        """
        assert (isinstance(eluent, str) or isinstance(eluent, bytes) or isinstance(eluent, String)), 'arg eluent wrong type'
    
        self.inst.get().addEluent(deref((convString(eluent)).get()))
    
    def clearEluents(self):
        """
        clearEluents(self) -> None
        Removes all eluents
        """
        self.inst.get().clearEluents()
    
    def getEluents(self):
        """
        getEluents(self) -> List[bytes]
        Returns a reference to the list of eluents
        """
        _r = self.inst.get().getEluents()
        py_result = []
        cdef libcpp_vector[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.append(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def addTimepoint(self,  timepoint ):
        """
        addTimepoint(self, timepoint: int ) -> None
        Adds a timepoint at the end of the timepoint array
        """
        assert isinstance(timepoint, int), 'arg timepoint wrong type'
    
        self.inst.get().addTimepoint((<int>timepoint))
    
    def clearTimepoints(self):
        """
        clearTimepoints(self) -> None
        Removes all timepoints
        """
        self.inst.get().clearTimepoints()
    
    def getTimepoints(self):
        """
        getTimepoints(self) -> List[int]
        Returns a reference to the list of timepoints
        """
        _r = self.inst.get().getTimepoints()
        cdef list py_result = _r
        return py_result
    
    def setPercentage(self,  eluent ,  timepoint ,  percentage ):
        """
        setPercentage(self, eluent: Union[bytes, str, String] , timepoint: int , percentage: int ) -> None
        Sets the percentage of 'eluent' at 'timepoint'
        """
        assert (isinstance(eluent, str) or isinstance(eluent, bytes) or isinstance(eluent, String)), 'arg eluent wrong type'
        assert isinstance(timepoint, int), 'arg timepoint wrong type'
        assert isinstance(percentage, int), 'arg percentage wrong type'
    
    
    
        self.inst.get().setPercentage(deref((convString(eluent)).get()), (<int>timepoint), (<unsigned int>percentage))
    
    def getPercentage(self,  eluent ,  timepoint ):
        """
        getPercentage(self, eluent: Union[bytes, str, String] , timepoint: int ) -> int
        Returns a const reference to the percentages
        """
        assert (isinstance(eluent, str) or isinstance(eluent, bytes) or isinstance(eluent, String)), 'arg eluent wrong type'
        assert isinstance(timepoint, int), 'arg timepoint wrong type'
    
    
        cdef unsigned int _r = self.inst.get().getPercentage(deref((convString(eluent)).get()), (<int>timepoint))
        py_result = <unsigned int>_r
        return py_result
    
    def clearPercentages(self):
        """
        clearPercentages(self) -> None
        Sets all percentage values to 0
        """
        self.inst.get().clearPercentages()
    
    def isValid(self):
        """
        isValid(self) -> bool
        Checks if the percentages of all timepoints add up to 100%
        """
        cdef bool _r = self.inst.get().isValid()
        py_result = <bool>_r
        return py_result 

cdef class IsotopeDistribution:
    """
    Cython implementation of _IsotopeDistribution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeDistribution.html>`_

    Isotope distribution class
    
    A container that holds an isotope distribution. It consists of mass values
    and their correspondent probabilities (stored in the intensity slot)
    
    Isotope distributions can be calculated using either the
    CoarseIsotopePatternGenerator for quantized atomic masses which group
    isotopes with the same atomic number. Alternatively, the
    FineIsotopePatternGenerator can be used that calculates hyperfine isotopic
    distributions
    
    This class only describes the container that holds the isotopic
    distribution, calculations are done using classes derived from
    IsotopePatternGenerator
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsotopeDistribution rv = IsotopeDistribution.__new__(IsotopeDistribution)
       rv.inst = shared_ptr[_IsotopeDistribution](new _IsotopeDistribution(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopeDistribution rv = IsotopeDistribution.__new__(IsotopeDistribution)
       rv.inst = shared_ptr[_IsotopeDistribution](new _IsotopeDistribution(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IsotopeDistribution](new _IsotopeDistribution())
    
    def _init_1(self, IsotopeDistribution in_0 ):
        """
        _init_1(self, in_0: IsotopeDistribution ) -> None
        """
        assert isinstance(in_0, IsotopeDistribution), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopeDistribution](new _IsotopeDistribution((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopeDistribution ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopeDistribution)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def set(self, list distribution ):
        """
        set(self, distribution: List[Peak1D] ) -> None
        Overwrites the container which holds the distribution using 'distribution'
        """
        assert isinstance(distribution, list) and all(isinstance(elemt_rec, Peak1D) for elemt_rec in distribution), 'arg distribution wrong type'
        cdef libcpp_vector[_Peak1D] * v0 = new libcpp_vector[_Peak1D]()
        cdef Peak1D item0
        for item0 in distribution:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().set(deref(v0))
        cdef libcpp_vector[_Peak1D].iterator it_distribution = v0.begin()
        replace_0 = []
        while it_distribution != v0.end():
            item0 = Peak1D.__new__(Peak1D)
            item0.inst = shared_ptr[_Peak1D](new _Peak1D(deref(it_distribution)))
            replace_0.append(item0)
            inc(it_distribution)
        distribution[:] = replace_0
        del v0
    
    def insert(self, double mass , float intensity ):
        """
        insert(self, mass: float , intensity: float ) -> None
        """
        assert isinstance(mass, float), 'arg mass wrong type'
        assert isinstance(intensity, float), 'arg intensity wrong type'
    
    
        self.inst.get().insert((<double>mass), (<float>intensity))
    
    def getContainer(self):
        """
        getContainer(self) -> List[Peak1D]
        Returns the container which holds the distribution
        """
        _r = self.inst.get().getContainer()
        py_result = []
        cdef libcpp_vector[_Peak1D].iterator it__r = _r.begin()
        cdef Peak1D item_py_result
        while it__r != _r.end():
           item_py_result = Peak1D.__new__(Peak1D)
           item_py_result.inst = shared_ptr[_Peak1D](new _Peak1D(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getMax(self):
        """
        getMax(self) -> float
        Returns the maximal weight isotope which is stored in the distribution
        """
        cdef double _r = self.inst.get().getMax()
        py_result = <double>_r
        return py_result
    
    def getMin(self):
        """
        getMin(self) -> float
        Returns the minimal weight isotope which is stored in the distribution
        """
        cdef double _r = self.inst.get().getMin()
        py_result = <double>_r
        return py_result
    
    def getMostAbundant(self):
        """
        getMostAbundant(self) -> Peak1D
        Returns the most abundant isotope which is stored in the distribution
        """
        cdef _Peak1D * _r = new _Peak1D(self.inst.get().getMostAbundant())
        cdef Peak1D py_result = Peak1D.__new__(Peak1D)
        py_result.inst = shared_ptr[_Peak1D](_r)
        return py_result
    
    def size(self):
        """
        size(self) -> int
        Returns the size of the distribution which is the number of isotopes in the distribution
        """
        cdef size_t _r = self.inst.get().size()
        py_result = <size_t>_r
        return py_result
    
    def clear(self):
        """
        clear(self) -> None
        Clears the distribution and resets max isotope to 0
        """
        self.inst.get().clear()
    
    def renormalize(self):
        """
        renormalize(self) -> None
        Renormalizes the sum of the probabilities of the isotopes to 1
        """
        self.inst.get().renormalize()
    
    def trimRight(self, double cutoff ):
        """
        trimRight(self, cutoff: float ) -> None
        Trims the right side of the isotope distribution to isotopes with a significant contribution
        """
        assert isinstance(cutoff, float), 'arg cutoff wrong type'
    
        self.inst.get().trimRight((<double>cutoff))
    
    def trimLeft(self, double cutoff ):
        """
        trimLeft(self, cutoff: float ) -> None
        Trims the left side of the isotope distribution to isotopes with a significant contribution
        """
        assert isinstance(cutoff, float), 'arg cutoff wrong type'
    
        self.inst.get().trimLeft((<double>cutoff))
    
    def merge(self, double in_0 , double in_1 ):
        """
        merge(self, in_0: float , in_1: float ) -> None
        Merges distributions of arbitrary data points with constant defined resolution
        """
        assert isinstance(in_0, float), 'arg in_0 wrong type'
        assert isinstance(in_1, float), 'arg in_1 wrong type'
    
    
        self.inst.get().merge((<double>in_0), (<double>in_1))
    
    def resize(self,  size ):
        """
        resize(self, size: int ) -> None
        Resizes distribution container
        """
        assert isinstance(size, int), 'arg size wrong type'
    
        self.inst.get().resize((<unsigned int>size))
    
    def trimIntensities(self, double cutoff ):
        """
        trimIntensities(self, cutoff: float ) -> None
        Remove intensities below the cutoff
        """
        assert isinstance(cutoff, float), 'arg cutoff wrong type'
    
        self.inst.get().trimIntensities((<double>cutoff))
    
    def sortByIntensity(self):
        """
        sortByIntensity(self) -> None
        Sort isotope distribution by intensity
        """
        self.inst.get().sortByIntensity()
    
    def sortByMass(self):
        """
        sortByMass(self) -> None
        Sort isotope distribution by mass
        """
        self.inst.get().sortByMass()
    
    def averageMass(self):
        """
        averageMass(self) -> float
        Compute average mass of isotope distribution (weighted average of all isotopes)
        """
        cdef double _r = self.inst.get().averageMass()
        py_result = <double>_r
        return py_result
    
    def __iter__(self):
        it = self.inst.get().begin()
        cdef Peak1D out
        while it != self.inst.get().end():
            out = Peak1D.__new__(Peak1D)
            out.inst = shared_ptr[_Peak1D](new _Peak1D(deref(it)))
            yield out
            inc(it)
    Sorted = __Sorted 

cdef class MSDataSqlConsumer:
    """
    Cython implementation of _MSDataSqlConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataSqlConsumer.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSDataSqlConsumer rv = MSDataSqlConsumer.__new__(MSDataSqlConsumer)
       rv.inst = shared_ptr[_MSDataSqlConsumer](new _MSDataSqlConsumer(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSDataSqlConsumer rv = MSDataSqlConsumer.__new__(MSDataSqlConsumer)
       rv.inst = shared_ptr[_MSDataSqlConsumer](new _MSDataSqlConsumer(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  filename ,  run_id ,  buffer_size , bool full_meta , bool lossy_compression , double linear_mass_acc ):
        """
        _init_0(self, filename: Union[bytes, str, String] , run_id: int , buffer_size: int , full_meta: bool , lossy_compression: bool , linear_mass_acc: float ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(run_id, int) and run_id >= 0, 'arg run_id wrong type'
        assert isinstance(buffer_size, int), 'arg buffer_size wrong type'
        assert isinstance(full_meta, pybool_t), 'arg full_meta wrong type'
        assert isinstance(lossy_compression, pybool_t), 'arg lossy_compression wrong type'
        assert isinstance(linear_mass_acc, float), 'arg linear_mass_acc wrong type'
    
    
    
    
    
    
        self.inst = shared_ptr[_MSDataSqlConsumer](new _MSDataSqlConsumer(deref((convString(filename)).get()), (<uint64_t>run_id), (<int>buffer_size), (<bool>full_meta), (<bool>lossy_compression), (<double>linear_mass_acc)))
    
    def _init_1(self, MSDataSqlConsumer in_0 ):
        """
        _init_1(self, in_0: MSDataSqlConsumer ) -> None
        """
        assert isinstance(in_0, MSDataSqlConsumer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSDataSqlConsumer](new _MSDataSqlConsumer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, filename: Union[bytes, str, String] , run_id: int , buffer_size: int , full_meta: bool , lossy_compression: bool , linear_mass_acc: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSDataSqlConsumer ) -> None
          :noindex:
    
        """
        if (len(args)==6) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int)) and (isinstance(args[3], pybool_t)) and (isinstance(args[4], pybool_t)) and (isinstance(args[5], float)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSDataSqlConsumer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def flush(self):
        """
        flush(self) -> None
        Flushes the data for good
        
        After calling this function, no more data is held in the buffer but the
        class is still able to receive new data
        """
        self.inst.get().flush()
    
    def consumeSpectrum(self, MSSpectrum s ):
        """
        consumeSpectrum(self, s: MSSpectrum ) -> None
        Write a spectrum to the output file
        """
        assert isinstance(s, MSSpectrum), 'arg s wrong type'
    
        self.inst.get().consumeSpectrum((deref(s.inst.get())))
    
    def consumeChromatogram(self, MSChromatogram c ):
        """
        consumeChromatogram(self, c: MSChromatogram ) -> None
        Write a chromatogram to the output file
        """
        assert isinstance(c, MSChromatogram), 'arg c wrong type'
    
        self.inst.get().consumeChromatogram((deref(c.inst.get())))
    
    def setExpectedSize(self,  expectedSpectra ,  expectedChromatograms ):
        """
        setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None
        """
        assert isinstance(expectedSpectra, int) and expectedSpectra >= 0, 'arg expectedSpectra wrong type'
        assert isinstance(expectedChromatograms, int) and expectedChromatograms >= 0, 'arg expectedChromatograms wrong type'
    
    
        self.inst.get().setExpectedSize((<size_t>expectedSpectra), (<size_t>expectedChromatograms))
    
    def setExperimentalSettings(self, ExperimentalSettings exp ):
        """
        setExperimentalSettings(self, exp: ExperimentalSettings ) -> None
        """
        assert isinstance(exp, ExperimentalSettings), 'arg exp wrong type'
    
        self.inst.get().setExperimentalSettings((deref(exp.inst.get()))) 

cdef class MSPFile:
    """
    Cython implementation of _MSPFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSPFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MSPFile rv = MSPFile.__new__(MSPFile)
       rv.inst = shared_ptr[_MSPFile](new _MSPFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MSPFile rv = MSPFile.__new__(MSPFile)
       rv.inst = shared_ptr[_MSPFile](new _MSPFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        File adapter for MSP files (NIST spectra library)
        """
        self.inst = shared_ptr[_MSPFile](new _MSPFile())
    
    def _init_1(self, MSPFile in_0 ):
        """
        _init_1(self, in_0: MSPFile ) -> None
        """
        assert isinstance(in_0, MSPFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MSPFile](new _MSPFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        File adapter for MSP files (NIST spectra library)

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MSPFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MSPFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename , AnnotatedMSRun exp ):
        """
        store(self, filename: Union[bytes, str, String] , exp: AnnotatedMSRun ) -> None
        Stores a map in a MSPFile file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(exp, AnnotatedMSRun), 'arg exp wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(exp.inst.get())))
    
    def load(self,  filename , PeptideIdentificationList ids , MSExperiment exp ):
        """
        load(self, filename: Union[bytes, str, String] , ids: PeptideIdentificationList , exp: MSExperiment ) -> None
        Loads a map from a MSPFile file
        
        
        :param exp: PeakMap which contains the spectra after reading
        :param filename: The filename of the experiment
        :param ids: Output parameter which contains the peptide identifications from the spectra annotations
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(ids.inst.get())), (deref(exp.inst.get()))) 

cdef class MzTabFile:
    """
    Cython implementation of _MzTabFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzTabFile rv = MzTabFile.__new__(MzTabFile)
       rv.inst = shared_ptr[_MzTabFile](new _MzTabFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzTabFile rv = MzTabFile.__new__(MzTabFile)
       rv.inst = shared_ptr[_MzTabFile](new _MzTabFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzTabFile](new _MzTabFile())
    
    def _init_1(self, MzTabFile in_0 ):
        """
        _init_1(self, in_0: MzTabFile ) -> None
        """
        assert isinstance(in_0, MzTabFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzTabFile](new _MzTabFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzTabFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzTabFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename , MzTab mz_tab ):
        """
        store(self, filename: Union[bytes, str, String] , mz_tab: MzTab ) -> None
        Stores MzTab file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(mz_tab, MzTab), 'arg mz_tab wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(mz_tab.inst.get())))
    
    def load(self,  filename , MzTab mz_tab ):
        """
        load(self, filename: Union[bytes, str, String] , mz_tab: MzTab ) -> None
        Loads MzTab file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(mz_tab, MzTab), 'arg mz_tab wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(mz_tab.inst.get()))) 

cdef class ParamXMLFile:
    """
    Cython implementation of _ParamXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamXMLFile.html>`_

    The file pendant of the Param class used to load and store the param
    datastructure as paramXML
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ParamXMLFile rv = ParamXMLFile.__new__(ParamXMLFile)
       rv.inst = shared_ptr[_ParamXMLFile](new _ParamXMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ParamXMLFile rv = ParamXMLFile.__new__(ParamXMLFile)
       rv.inst = shared_ptr[_ParamXMLFile](new _ParamXMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ParamXMLFile](new _ParamXMLFile())
    
    def _init_1(self, ParamXMLFile in_0 ):
        """
        _init_1(self, in_0: ParamXMLFile ) -> None
        """
        assert isinstance(in_0, ParamXMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ParamXMLFile](new _ParamXMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ParamXMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ParamXMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , Param param ):
        """
        load(self, filename: Union[bytes, str, String] , param: Param ) -> None
        Read XML file
        
        
        :param filename: The file from where to read the Param object
        :param param: The param object where the read data should be stored
        :raises:
          Exception: FileNotFound is thrown if the file could not be found
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(param.inst.get())))
    
    def store(self,  filename , Param param ):
        """
        store(self, filename: Union[bytes, str, String] , param: Param ) -> None
        Write XML file
        
        
        :param filename: The filename where the param data structure should be stored
        :param param: The Param class that should be stored in the file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(param.inst.get()))) 

cdef class ResidueDB:
    """
    Cython implementation of _ResidueDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ResidueDB.html>`_
    """

    
    def getNumberOfResidues(self):
        """
        getNumberOfResidues(self) -> int
        Returns the number of residues stored
        """
        cdef size_t _r = self.inst.get().getNumberOfResidues()
        py_result = <size_t>_r
        return py_result
    
    def getNumberOfModifiedResidues(self):
        """
        getNumberOfModifiedResidues(self) -> int
        Returns the number of modified residues stored
        """
        cdef size_t _r = self.inst.get().getNumberOfModifiedResidues()
        py_result = <size_t>_r
        return py_result
    
    def getResidue(self,  name ):
        """
        getResidue(self, name: Union[bytes, str, String] ) -> Residue
        Returns a pointer to the residue with name, 3 letter code or 1 letter code name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef const _Residue * __r = (self.inst.get().getResidue(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _Residue * _r = new _Residue(deref(__r))
        cdef Residue py_result = Residue.__new__(Residue)
        py_result.inst = shared_ptr[_Residue](_r)
        return py_result
    
    def _getModifiedResidue_0(self,  name ):
        """
        _getModifiedResidue_0(self, name: Union[bytes, str, String] ) -> Residue
        Returns a pointer to a modified residue given a modification name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef const _Residue * __r = (self.inst.get().getModifiedResidue(deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _Residue * _r = new _Residue(deref(__r))
        cdef Residue py_result = Residue.__new__(Residue)
        py_result.inst = shared_ptr[_Residue](_r)
        return py_result
    
    def _getModifiedResidue_1(self, Residue residue ,  name ):
        """
        _getModifiedResidue_1(self, residue: Residue , name: Union[bytes, str, String] ) -> Residue
        Returns a pointer to a modified residue given a residue and a modification name
        """
        assert isinstance(residue, Residue), 'arg residue wrong type'
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
    
        cdef const _Residue * __r = (self.inst.get().getModifiedResidue((residue.inst.get()), deref((convString(name)).get())))
        if __r == NULL:
            return None
        cdef _Residue * _r = new _Residue(deref(__r))
        cdef Residue py_result = Residue.__new__(Residue)
        py_result.inst = shared_ptr[_Residue](_r)
        return py_result
    
    def getModifiedResidue(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getModifiedResidue(self, name: Union[bytes, str, String] ) -> Residue
          :noindex:
        
        Returns a pointer to a modified residue given a modification name

        
        .. rubric:: Overload:
        .. py:function:: getModifiedResidue(self, residue: Residue , name: Union[bytes, str, String] ) -> Residue
          :noindex:
        
        Returns a pointer to a modified residue given a residue and a modification name
    
        """
        if (len(args)==1) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))):
            return self._getModifiedResidue_0(*args)
        elif (len(args)==2) and (isinstance(args[0], Residue)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._getModifiedResidue_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getResidues(self,  residue_set ):
        """
        getResidues(self, residue_set: Union[bytes, str, String] ) -> Set[Residue]
        Returns a set of all residues stored in this residue db
        """
        assert (isinstance(residue_set, str) or isinstance(residue_set, bytes) or isinstance(residue_set, String)), 'arg residue_set wrong type'
    
        _r = self.inst.get().getResidues(deref((convString(residue_set)).get()))
        py_result = set()
        cdef libcpp_set[const _Residue *].iterator it__r = _r.begin()
        cdef Residue item_py_result
        while it__r != _r.end():
           item_py_result = Residue.__new__(Residue)
           item_py_result.inst = shared_ptr[_Residue](new _Residue(deref(deref(it__r))))
           py_result.add(item_py_result)
           inc(it__r)
        return py_result
    
    def getResidueSets(self):
        """
        getResidueSets(self) -> Set[bytes]
        Returns all residue sets that are registered which this instance
        """
        _r = self.inst.get().getResidueSets()
        py_result = set()
        cdef libcpp_set[_String].iterator it__r = _r.begin()
        while it__r != _r.end():
           py_result.add(<char*>deref(it__r).c_str())
           inc(it__r)
        return py_result
    
    def hasResidue(self,  name ):
        """
        hasResidue(self, name: Union[bytes, str, String] ) -> bool
        Returns true if the db contains a residue with the given name
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
    
        cdef bool _r = self.inst.get().hasResidue(deref((convString(name)).get()))
        py_result = <bool>_r
        return py_result
    
    # NOTE: using shared_ptr for a singleton will lead to segfaults, use raw ptr instead
    # cdef AutowrapPtrHolder[_ResidueDB] inst

    def __init__(self):
      self.inst = AutowrapPtrHolder[_ResidueDB](_getInstance_ResidueDB())

    def __dealloc__(self):
      # Careful here, the wrapped ptr is a single instance and we should not
      # reset it, therefore use 'wrap-manual-dealloc'
      pass 

cdef class SpectrumAccessQuadMZTransforming:
    """
    Cython implementation of _SpectrumAccessQuadMZTransforming

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessQuadMZTransforming.html>`_
      -- Inherits from ['SpectrumAccessTransforming']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SpectrumAccessQuadMZTransforming rv = SpectrumAccessQuadMZTransforming.__new__(SpectrumAccessQuadMZTransforming)
       rv.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumAccessQuadMZTransforming rv = SpectrumAccessQuadMZTransforming.__new__(SpectrumAccessQuadMZTransforming)
       rv.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        pass
    
    def _init_1(self, SpectrumAccessQuadMZTransforming in_0 ):
        """
        _init_1(self, in_0: SpectrumAccessQuadMZTransforming ) -> None
        """
        assert isinstance(in_0, SpectrumAccessQuadMZTransforming), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming((deref(in_0.inst.get()))))
    
    def _init_2(self, SpectrumAccessOpenMS in_0 , double a , double b , double c , bool ppm ):
        """
        _init_2(self, in_0: SpectrumAccessOpenMS , a: float , b: float , c: float , ppm: bool ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMS), 'arg in_0 wrong type'
        assert isinstance(a, float), 'arg a wrong type'
        assert isinstance(b, float), 'arg b wrong type'
        assert isinstance(c, float), 'arg c wrong type'
        assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMS] input_in_0 = in_0.inst
    
    
    
    
        self.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming(input_in_0, (<double>a), (<double>b), (<double>c), (<bool>ppm)))
    
    def _init_3(self, SpectrumAccessOpenMSCached in_0 , double a , double b , double c , bool ppm ):
        """
        _init_3(self, in_0: SpectrumAccessOpenMSCached , a: float , b: float , c: float , ppm: bool ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMSCached), 'arg in_0 wrong type'
        assert isinstance(a, float), 'arg a wrong type'
        assert isinstance(b, float), 'arg b wrong type'
        assert isinstance(c, float), 'arg c wrong type'
        assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMSCached] input_in_0 = in_0.inst
    
    
    
    
        self.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming(input_in_0, (<double>a), (<double>b), (<double>c), (<bool>ppm)))
    
    def _init_4(self, SpectrumAccessOpenMSInMemory in_0 , double a , double b , double c , bool ppm ):
        """
        _init_4(self, in_0: SpectrumAccessOpenMSInMemory , a: float , b: float , c: float , ppm: bool ) -> None
        """
        assert isinstance(in_0, SpectrumAccessOpenMSInMemory), 'arg in_0 wrong type'
        assert isinstance(a, float), 'arg a wrong type'
        assert isinstance(b, float), 'arg b wrong type'
        assert isinstance(c, float), 'arg c wrong type'
        assert isinstance(ppm, pybool_t), 'arg ppm wrong type'
        cdef shared_ptr[_SpectrumAccessOpenMSInMemory] input_in_0 = in_0.inst
    
    
    
    
        self.inst = shared_ptr[_SpectrumAccessQuadMZTransforming](new _SpectrumAccessQuadMZTransforming(input_in_0, (<double>a), (<double>b), (<double>c), (<bool>ppm)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessQuadMZTransforming ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMS , a: float , b: float , c: float , ppm: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMSCached , a: float , b: float , c: float , ppm: bool ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumAccessOpenMSInMemory , a: float , b: float , c: float , ppm: bool ) -> None
          :noindex:
    
        """
        if kwargs.get("__createUnsafeObject__") is True:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumAccessQuadMZTransforming)):
             self._init_1(*args)
        elif (len(args)==5) and (isinstance(args[0], SpectrumAccessOpenMS)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)) and (isinstance(args[4], pybool_t)):
             self._init_2(*args)
        elif (len(args)==5) and (isinstance(args[0], SpectrumAccessOpenMSCached)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)) and (isinstance(args[4], pybool_t)):
             self._init_3(*args)
        elif (len(args)==5) and (isinstance(args[0], SpectrumAccessOpenMSInMemory)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)) and (isinstance(args[4], pybool_t)):
             self._init_4(*args)
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

cdef class StablePairFinder:
    """
    Cython implementation of _StablePairFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1StablePairFinder.html>`_
      -- Inherits from ['BaseGroupFinder']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_StablePairFinder](new _StablePairFinder())
    
    def run(self, list input_maps , ConsensusMap result_map ):
        """
        run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None
        """
        assert isinstance(input_maps, list) and all(isinstance(elemt_rec, ConsensusMap) for elemt_rec in input_maps), 'arg input_maps wrong type'
        assert isinstance(result_map, ConsensusMap), 'arg result_map wrong type'
        cdef libcpp_vector[_ConsensusMap] * v0 = new libcpp_vector[_ConsensusMap]()
        cdef ConsensusMap item0
        for item0 in input_maps:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().run(deref(v0), (deref(result_map.inst.get())))
        cdef libcpp_vector[_ConsensusMap].iterator it_input_maps = v0.begin()
        replace_0 = []
        while it_input_maps != v0.end():
            item0 = ConsensusMap.__new__(ConsensusMap)
            item0.inst = shared_ptr[_ConsensusMap](new _ConsensusMap(deref(it_input_maps)))
            replace_0.append(item0)
            inc(it_input_maps)
        input_maps[:] = replace_0
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
