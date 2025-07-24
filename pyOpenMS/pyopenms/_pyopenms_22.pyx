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
def __static_PercolatorInfile_store( pin_file , PeptideIdentificationList peptide_ids , list feature_set , bytes in_3 ,  min_charge ,  max_charge ):
    """
    __static_PercolatorInfile_store(pin_file: Union[bytes, str, String] , peptide_ids: PeptideIdentificationList , feature_set: List[bytes] , in_3: bytes , min_charge: int , max_charge: int ) -> None
    """
    assert (isinstance(pin_file, str) or isinstance(pin_file, bytes) or isinstance(pin_file, String)), 'arg pin_file wrong type'
    assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    assert isinstance(feature_set, list) and all(isinstance(li, bytes) for li in feature_set), 'arg feature_set wrong type'
    assert isinstance(in_3, bytes), 'arg in_3 wrong type'
    assert isinstance(min_charge, int), 'arg min_charge wrong type'
    assert isinstance(max_charge, int), 'arg max_charge wrong type'


    cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
    cdef bytes item2
    for item2 in feature_set:
       v2.push_back(_String(<char *>item2))



    _store_PercolatorInfile(deref((convString(pin_file)).get()), (deref(peptide_ids.inst.get())), deref(v2), (<libcpp_string>in_3), (<int>min_charge), (<int>max_charge))
    del v2 

cdef class ChargedIndexSet:
    """
    Cython implementation of _ChargedIndexSet

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChargedIndexSet.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property charge:
        def __set__(self,  charge):
        
            self.inst.get().charge = (<int>charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().charge
            py_result = <int>_r
            return py_result
    
    def __init__(self):
        """
        __init__(self) -> None
        Index set with associated charge estimate
        """
        self.inst = shared_ptr[_ChargedIndexSet](new _ChargedIndexSet()) 

cdef class EmgGradientDescent:
    """
    Cython implementation of _EmgGradientDescent

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgGradientDescent.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef EmgGradientDescent rv = EmgGradientDescent.__new__(EmgGradientDescent)
       rv.inst = shared_ptr[_EmgGradientDescent](new _EmgGradientDescent(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef EmgGradientDescent rv = EmgGradientDescent.__new__(EmgGradientDescent)
       rv.inst = shared_ptr[_EmgGradientDescent](new _EmgGradientDescent(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Compute the area, background and shape metrics of a peak
        """
        self.inst = shared_ptr[_EmgGradientDescent](new _EmgGradientDescent())
    
    def _init_1(self, EmgGradientDescent in_0 ):
        """
        _init_1(self, in_0: EmgGradientDescent ) -> None
        """
        assert isinstance(in_0, EmgGradientDescent), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_EmgGradientDescent](new _EmgGradientDescent((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Compute the area, background and shape metrics of a peak

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: EmgGradientDescent ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], EmgGradientDescent)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getDefaultParameters(self, Param in_0 ):
        """
        getDefaultParameters(self, in_0: Param ) -> None
        """
        assert isinstance(in_0, Param), 'arg in_0 wrong type'
    
        self.inst.get().getDefaultParameters((deref(in_0.inst.get())))
    
    def _fitEMGPeakModel_0(self, MSChromatogram input_peak , MSChromatogram output_peak ):
        """
        _fitEMGPeakModel_0(self, input_peak: MSChromatogram , output_peak: MSChromatogram ) -> None
        """
        assert isinstance(input_peak, MSChromatogram), 'arg input_peak wrong type'
        assert isinstance(output_peak, MSChromatogram), 'arg output_peak wrong type'
    
    
        self.inst.get().fitEMGPeakModel((deref(input_peak.inst.get())), (deref(output_peak.inst.get())))
    
    def _fitEMGPeakModel_1(self, MSSpectrum input_peak , MSSpectrum output_peak ):
        """
        _fitEMGPeakModel_1(self, input_peak: MSSpectrum , output_peak: MSSpectrum ) -> None
        """
        assert isinstance(input_peak, MSSpectrum), 'arg input_peak wrong type'
        assert isinstance(output_peak, MSSpectrum), 'arg output_peak wrong type'
    
    
        self.inst.get().fitEMGPeakModel((deref(input_peak.inst.get())), (deref(output_peak.inst.get())))
    
    def _fitEMGPeakModel_2(self, MSChromatogram input_peak , MSChromatogram output_peak , double left_pos , double right_pos ):
        """
        _fitEMGPeakModel_2(self, input_peak: MSChromatogram , output_peak: MSChromatogram , left_pos: float , right_pos: float ) -> None
        """
        assert isinstance(input_peak, MSChromatogram), 'arg input_peak wrong type'
        assert isinstance(output_peak, MSChromatogram), 'arg output_peak wrong type'
        assert isinstance(left_pos, float), 'arg left_pos wrong type'
        assert isinstance(right_pos, float), 'arg right_pos wrong type'
    
    
    
    
        self.inst.get().fitEMGPeakModel((deref(input_peak.inst.get())), (deref(output_peak.inst.get())), (<double>left_pos), (<double>right_pos))
    
    def _fitEMGPeakModel_3(self, MSSpectrum input_peak , MSSpectrum output_peak , double left_pos , double right_pos ):
        """
        _fitEMGPeakModel_3(self, input_peak: MSSpectrum , output_peak: MSSpectrum , left_pos: float , right_pos: float ) -> None
        """
        assert isinstance(input_peak, MSSpectrum), 'arg input_peak wrong type'
        assert isinstance(output_peak, MSSpectrum), 'arg output_peak wrong type'
        assert isinstance(left_pos, float), 'arg left_pos wrong type'
        assert isinstance(right_pos, float), 'arg right_pos wrong type'
    
    
    
    
        self.inst.get().fitEMGPeakModel((deref(input_peak.inst.get())), (deref(output_peak.inst.get())), (<double>left_pos), (<double>right_pos))
    
    def fitEMGPeakModel(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: fitEMGPeakModel(self, input_peak: MSChromatogram , output_peak: MSChromatogram ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: fitEMGPeakModel(self, input_peak: MSSpectrum , output_peak: MSSpectrum ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: fitEMGPeakModel(self, input_peak: MSChromatogram , output_peak: MSChromatogram , left_pos: float , right_pos: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: fitEMGPeakModel(self, input_peak: MSSpectrum , output_peak: MSSpectrum , left_pos: float , right_pos: float ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], MSChromatogram)):
            return self._fitEMGPeakModel_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], MSSpectrum)):
            return self._fitEMGPeakModel_1(*args)
        elif (len(args)==4) and (isinstance(args[0], MSChromatogram)) and (isinstance(args[1], MSChromatogram)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._fitEMGPeakModel_2(*args)
        elif (len(args)==4) and (isinstance(args[0], MSSpectrum)) and (isinstance(args[1], MSSpectrum)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._fitEMGPeakModel_3(*args)
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

cdef class FeatureFileOptions:
    """
    Cython implementation of _FeatureFileOptions

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFileOptions.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef FeatureFileOptions rv = FeatureFileOptions.__new__(FeatureFileOptions)
       rv.inst = shared_ptr[_FeatureFileOptions](new _FeatureFileOptions(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef FeatureFileOptions rv = FeatureFileOptions.__new__(FeatureFileOptions)
       rv.inst = shared_ptr[_FeatureFileOptions](new _FeatureFileOptions(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Options for loading files containing features.
        """
        self.inst = shared_ptr[_FeatureFileOptions](new _FeatureFileOptions())
    
    def _init_1(self, FeatureFileOptions in_0 ):
        """
        _init_1(self, in_0: FeatureFileOptions ) -> None
        """
        assert isinstance(in_0, FeatureFileOptions), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_FeatureFileOptions](new _FeatureFileOptions((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Options for loading files containing features.

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: FeatureFileOptions ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], FeatureFileOptions)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setMetadataOnly(self, bool in_0 ):
        """
        setMetadataOnly(self, in_0: bool ) -> None
        Sets whether or not to load only meta data
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setMetadataOnly((<bool>in_0))
    
    def getMetadataOnly(self):
        """
        getMetadataOnly(self) -> bool
        Returns whether or not to load only meta data
        """
        cdef bool _r = self.inst.get().getMetadataOnly()
        py_result = <bool>_r
        return py_result
    
    def setSizeOnly(self, bool in_0 ):
        """
        setSizeOnly(self, in_0: bool ) -> None
        Sets whether or not to load only feature count
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setSizeOnly((<bool>in_0))
    
    def getSizeOnly(self):
        """
        getSizeOnly(self) -> bool
        Returns whether or not to load only meta data
        """
        cdef bool _r = self.inst.get().getSizeOnly()
        py_result = <bool>_r
        return py_result
    
    def setLoadConvexHull(self, bool in_0 ):
        """
        setLoadConvexHull(self, in_0: bool ) -> None
        Sets whether or not to load convex hull
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setLoadConvexHull((<bool>in_0))
    
    def getLoadConvexHull(self):
        """
        getLoadConvexHull(self) -> bool
        Returns whether or not to load convex hull
        """
        cdef bool _r = self.inst.get().getLoadConvexHull()
        py_result = <bool>_r
        return py_result
    
    def setLoadSubordinates(self, bool in_0 ):
        """
        setLoadSubordinates(self, in_0: bool ) -> None
        Sets whether or not load subordinates
        """
        assert isinstance(in_0, pybool_t), 'arg in_0 wrong type'
    
        self.inst.get().setLoadSubordinates((<bool>in_0))
    
    def getLoadSubordinates(self):
        """
        getLoadSubordinates(self) -> bool
        Returns whether or not to load subordinates
        """
        cdef bool _r = self.inst.get().getLoadSubordinates()
        py_result = <bool>_r
        return py_result
    
    def setRTRange(self, DRange1 range_ ):
        """
        setRTRange(self, range_: DRange1 ) -> None
        Restricts the range of RT values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setRTRange((deref(range_.inst.get())))
    
    def hasRTRange(self):
        """
        hasRTRange(self) -> bool
        Returns true if an RT range has been set
        """
        cdef bool _r = self.inst.get().hasRTRange()
        py_result = <bool>_r
        return py_result
    
    def getRTRange(self):
        """
        getRTRange(self) -> DRange1
        Returns the RT range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getRTRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def setMZRange(self, DRange1 range_ ):
        """
        setMZRange(self, range_: DRange1 ) -> None
        Restricts the range of MZ values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setMZRange((deref(range_.inst.get())))
    
    def hasMZRange(self):
        """
        hasMZRange(self) -> bool
        Returns true if an MZ range has been set
        """
        cdef bool _r = self.inst.get().hasMZRange()
        py_result = <bool>_r
        return py_result
    
    def getMZRange(self):
        """
        getMZRange(self) -> DRange1
        Returns the MZ range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getMZRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result
    
    def setIntensityRange(self, DRange1 range_ ):
        """
        setIntensityRange(self, range_: DRange1 ) -> None
        Restricts the range of intensity values for peaks to load
        """
        assert isinstance(range_, DRange1), 'arg range_ wrong type'
    
        self.inst.get().setIntensityRange((deref(range_.inst.get())))
    
    def hasIntensityRange(self):
        """
        hasIntensityRange(self) -> bool
        Returns true if an intensity range has been set
        """
        cdef bool _r = self.inst.get().hasIntensityRange()
        py_result = <bool>_r
        return py_result
    
    def getIntensityRange(self):
        """
        getIntensityRange(self) -> DRange1
        Returns the intensity range
        """
        cdef _DRange1 * _r = new _DRange1(self.inst.get().getIntensityRange())
        cdef DRange1 py_result = DRange1.__new__(DRange1)
        py_result.inst = shared_ptr[_DRange1](_r)
        return py_result 

cdef class FeatureFinderIdentificationAlgorithm:
    """
    Cython implementation of _FeatureFinderIdentificationAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FeatureFinderIdentificationAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']

    Algorithm class for FeatureFinderIdentification
    
    External IDs (peptides_ext, proteins_ext) may be empty,
    in which case no machine learning or FDR estimation will be performed.
    Optional seeds from e.g. untargeted FeatureFinders can be added with
    seeds.
    Results will be written to features .
    Caution: peptide IDs will be shrunk to best hit, FFid metavalues added
    and potential seed IDs added.
    
    Usage:
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_FeatureFinderIdentificationAlgorithm](new _FeatureFinderIdentificationAlgorithm())
    
    def _run_0(self, PeptideIdentificationList peptides , list proteins , PeptideIdentificationList peptides_ext , list proteins_ext , FeatureMap features ):
        """
        _run_0(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap ) -> None
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides_ext, PeptideIdentificationList), 'arg peptides_ext wrong type'
        assert isinstance(proteins_ext, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins_ext), 'arg proteins_ext wrong type'
        assert isinstance(features, FeatureMap), 'arg features wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in proteins:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_ProteinIdentification] * v3 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item3
        for item3 in proteins_ext:
            v3.push_back(deref(item3.inst.get()))
    
        self.inst.get().run((deref(peptides.inst.get())), deref(v1), (deref(peptides_ext.inst.get())), deref(v3), (deref(features.inst.get())))
        del v3
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v1.begin()
        replace_0 = []
        while it_proteins != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item1)
            inc(it_proteins)
        proteins[:] = replace_0
        del v1
    
    def _run_1(self, PeptideIdentificationList peptides , list proteins , PeptideIdentificationList peptides_ext , list proteins_ext , FeatureMap features , FeatureMap seeds ):
        """
        _run_1(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap ) -> None
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides_ext, PeptideIdentificationList), 'arg peptides_ext wrong type'
        assert isinstance(proteins_ext, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins_ext), 'arg proteins_ext wrong type'
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        assert isinstance(seeds, FeatureMap), 'arg seeds wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in proteins:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_ProteinIdentification] * v3 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item3
        for item3 in proteins_ext:
            v3.push_back(deref(item3.inst.get()))
    
    
        self.inst.get().run((deref(peptides.inst.get())), deref(v1), (deref(peptides_ext.inst.get())), deref(v3), (deref(features.inst.get())), (deref(seeds.inst.get())))
        del v3
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v1.begin()
        replace_0 = []
        while it_proteins != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item1)
            inc(it_proteins)
        proteins[:] = replace_0
        del v1
    
    def _run_2(self, PeptideIdentificationList peptides , list proteins , PeptideIdentificationList peptides_ext , list proteins_ext , FeatureMap features , FeatureMap seeds ,  spectra_file ):
        """
        _run_2(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap , spectra_file: String ) -> None
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        :param spectra_file: Path will be stored in features in case the MSExperiment has no proper primaryMSRunPath
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides_ext, PeptideIdentificationList), 'arg peptides_ext wrong type'
        assert isinstance(proteins_ext, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins_ext), 'arg proteins_ext wrong type'
        assert isinstance(features, FeatureMap), 'arg features wrong type'
        assert isinstance(seeds, FeatureMap), 'arg seeds wrong type'
        assert isinstance(spectra_file, String), 'arg spectra_file wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in proteins:
            v1.push_back(deref(item1.inst.get()))
    
        cdef libcpp_vector[_ProteinIdentification] * v3 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item3
        for item3 in proteins_ext:
            v3.push_back(deref(item3.inst.get()))
    
    
    
        self.inst.get().run((deref(peptides.inst.get())), deref(v1), (deref(peptides_ext.inst.get())), deref(v3), (deref(features.inst.get())), (deref(seeds.inst.get())), deref((<String>spectra_file).inst.get()))
        del v3
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v1.begin()
        replace_0 = []
        while it_proteins != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item1)
            inc(it_proteins)
        proteins[:] = replace_0
        del v1
    
    def run(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap ) -> None
          :noindex:
        
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        
        .. rubric:: Overload:
        .. py:function:: run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap ) -> None
          :noindex:
        
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        
        .. rubric:: Overload:
        .. py:function:: run(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , peptides_ext: PeptideIdentificationList , proteins_ext: List[ProteinIdentification] , features: FeatureMap , seeds: FeatureMap , spectra_file: String ) -> None
          :noindex:
        
        Run feature detection
        
        
        :param peptides: Vector of identified peptides
        :param proteins: Vector of identified proteins
        :param peptides_ext: Vector of external identified peptides, can be used to transfer ids from other runs
        :param proteins_ext: Vector of external identified proteins, can be used to transfer ids from other runs
        :param features: Feature detection results will be added here
        :param seeds: Optional seeds for feature detection from e.g. untargeted FeatureFinders
        :param spectra_file: Path will be stored in features in case the MSExperiment has no proper primaryMSRunPath
    
        """
        if (len(args)==5) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[3])) and (isinstance(args[4], FeatureMap)):
            return self._run_0(*args)
        elif (len(args)==6) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[3])) and (isinstance(args[4], FeatureMap)) and (isinstance(args[5], FeatureMap)):
            return self._run_1(*args)
        elif (len(args)==7) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and (isinstance(args[3], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[3])) and (isinstance(args[4], FeatureMap)) and (isinstance(args[5], FeatureMap)) and (isinstance(args[6], String)):
            return self._run_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def runOnCandidates(self, FeatureMap features ):
        """
        runOnCandidates(self, features: FeatureMap ) -> None
        Run feature detection on identified features (e.g. loaded from an IdXML file)
        """
        assert isinstance(features, FeatureMap), 'arg features wrong type'
    
        self.inst.get().runOnCandidates((deref(features.inst.get())))
    
    def setMSData(self, MSExperiment in_0 ):
        """
        setMSData(self, in_0: MSExperiment ) -> None
        Sets ms data
        """
        assert isinstance(in_0, MSExperiment), 'arg in_0 wrong type'
    
        self.inst.get().setMSData((deref(in_0.inst.get())))
    
    def getMSData(self):
        """
        getMSData(self) -> MSExperiment
        Returns ms data as MSExperiment
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getMSData())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def getChromatograms(self):
        """
        getChromatograms(self) -> MSExperiment
        Returns chromatogram data as MSExperiment
        """
        cdef _MSExperiment * _r = new _MSExperiment(self.inst.get().getChromatograms())
        cdef MSExperiment py_result = MSExperiment.__new__(MSExperiment)
        py_result.inst = shared_ptr[_MSExperiment](_r)
        return py_result
    
    def getLibrary(self):
        """
        getLibrary(self) -> TargetedExperiment
        Returns constructed assay library
        """
        cdef _TargetedExperiment * _r = new _TargetedExperiment(self.inst.get().getLibrary())
        cdef TargetedExperiment py_result = TargetedExperiment.__new__(TargetedExperiment)
        py_result.inst = shared_ptr[_TargetedExperiment](_r)
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

cdef class IBSpectraFile:
    """
    Cython implementation of _IBSpectraFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IBSpectraFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IBSpectraFile rv = IBSpectraFile.__new__(IBSpectraFile)
       rv.inst = shared_ptr[_IBSpectraFile](new _IBSpectraFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IBSpectraFile rv = IBSpectraFile.__new__(IBSpectraFile)
       rv.inst = shared_ptr[_IBSpectraFile](new _IBSpectraFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Implements the export of consensusmaps into the IBSpectra format used by isobar to load quantification results
        """
        self.inst = shared_ptr[_IBSpectraFile](new _IBSpectraFile())
    
    def _init_1(self, IBSpectraFile in_0 ):
        """
        _init_1(self, in_0: IBSpectraFile ) -> None
        """
        assert isinstance(in_0, IBSpectraFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IBSpectraFile](new _IBSpectraFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Implements the export of consensusmaps into the IBSpectra format used by isobar to load quantification results

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IBSpectraFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IBSpectraFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename , ConsensusMap cm ):
        """
        store(self, filename: Union[bytes, str, String] , cm: ConsensusMap ) -> None
        Writes the contents of the ConsensusMap cm into the file named by filename
        
        
        :param filename: The name of the file where the contents of cm should be stored
        :param cm: The ConsensusMap that should be exported to filename
        :raises:
          Exception: InvalidParameter if the ConsensusMap does not hold the result of an isobaric quantification experiment (e.g., itraq)
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(cm, ConsensusMap), 'arg cm wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(cm.inst.get()))) 

cdef class IsotopeCluster:
    """
    Cython implementation of _IsotopeCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsotopeCluster.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property peaks:
        def __set__(self, ChargedIndexSet peaks):
        
            self.inst.get().peaks = (deref(peaks.inst.get()))
        
    
        def __get__(self):
            cdef _ChargedIndexSet * _r = new _ChargedIndexSet(self.inst.get().peaks)
            cdef ChargedIndexSet py_result = ChargedIndexSet.__new__(ChargedIndexSet)
            py_result.inst = shared_ptr[_ChargedIndexSet](_r)
            return py_result
    
    property scans:
        def __set__(self, list scans):
            cdef libcpp_vector[size_t] v0 = scans
            self.inst.get().scans = v0
            
    
        def __get__(self):
            _r = self.inst.get().scans
            cdef list py_result = _r
            return py_result
    
    def __copy__(self):
       cdef IsotopeCluster rv = IsotopeCluster.__new__(IsotopeCluster)
       rv.inst = shared_ptr[_IsotopeCluster](new _IsotopeCluster(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsotopeCluster rv = IsotopeCluster.__new__(IsotopeCluster)
       rv.inst = shared_ptr[_IsotopeCluster](new _IsotopeCluster(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Stores information about an isotopic cluster (i.e. potential peptide charge variants)
        """
        self.inst = shared_ptr[_IsotopeCluster](new _IsotopeCluster())
    
    def _init_1(self, IsotopeCluster in_0 ):
        """
        _init_1(self, in_0: IsotopeCluster ) -> None
        """
        assert isinstance(in_0, IsotopeCluster), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsotopeCluster](new _IsotopeCluster((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Stores information about an isotopic cluster (i.e. potential peptide charge variants)

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsotopeCluster ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsotopeCluster)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MapAlignmentAlgorithmKD:
    """
    Cython implementation of _MapAlignmentAlgorithmKD

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmKD.html>`_

    An efficient reference-free feature map alignment algorithm for unlabeled data
    
    This algorithm uses a kd-tree to efficiently compute conflict-free connected components (CCC)
    in a compatibility graph on feature data. This graph is comprised of nodes corresponding
    to features and edges connecting features f and f' iff both are within each other's tolerance
    windows (wrt. RT and m/z difference). CCCs are those CCs that do not contain multiple features
    from the same input map, and whose features all have the same charge state
    
    All CCCs above a user-specified minimum size are considered true sets of corresponding features
    and based on these, LOWESS transformations are computed for each input map such that the average
    deviation from the mean retention time within all CCCs is minimized
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MapAlignmentAlgorithmKD rv = MapAlignmentAlgorithmKD.__new__(MapAlignmentAlgorithmKD)
       rv.inst = shared_ptr[_MapAlignmentAlgorithmKD](new _MapAlignmentAlgorithmKD(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MapAlignmentAlgorithmKD rv = MapAlignmentAlgorithmKD.__new__(MapAlignmentAlgorithmKD)
       rv.inst = shared_ptr[_MapAlignmentAlgorithmKD](new _MapAlignmentAlgorithmKD(deref(self.inst.get())))
       return rv
    
    def _init_0(self,  num_maps , Param param ):
        """
        _init_0(self, num_maps: int , param: Param ) -> None
        """
        assert isinstance(num_maps, int) and num_maps >= 0, 'arg num_maps wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
    
    
        self.inst = shared_ptr[_MapAlignmentAlgorithmKD](new _MapAlignmentAlgorithmKD((<size_t>num_maps), (deref(param.inst.get()))))
    
    def _init_1(self, MapAlignmentAlgorithmKD in_0 ):
        """
        _init_1(self, in_0: MapAlignmentAlgorithmKD ) -> None
        """
        assert isinstance(in_0, MapAlignmentAlgorithmKD), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MapAlignmentAlgorithmKD](new _MapAlignmentAlgorithmKD((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, num_maps: int , param: Param ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MapAlignmentAlgorithmKD ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], Param)):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MapAlignmentAlgorithmKD)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def addRTFitData(self, KDTreeFeatureMaps kd_data ):
        """
        addRTFitData(self, kd_data: KDTreeFeatureMaps ) -> None
        Compute data points needed for RT transformation in the current `kd_data`, add to `fit_data_`
        """
        assert isinstance(kd_data, KDTreeFeatureMaps), 'arg kd_data wrong type'
    
        self.inst.get().addRTFitData((deref(kd_data.inst.get())))
    
    def fitLOWESS(self):
        """
        fitLOWESS(self) -> None
        Fit LOWESS to fit_data_, store final models in `transformations_`
        """
        self.inst.get().fitLOWESS()
    
    def transform(self, KDTreeFeatureMaps kd_data ):
        """
        transform(self, kd_data: KDTreeFeatureMaps ) -> None
        Transform RTs for `kd_data`
        """
        assert isinstance(kd_data, KDTreeFeatureMaps), 'arg kd_data wrong type'
    
        self.inst.get().transform((deref(kd_data.inst.get()))) 

cdef class MetaInfoDescription:
    """
    Cython implementation of _MetaInfoDescription

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaInfoDescription.html>`_
      -- Inherits from ['MetaInfoInterface']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaInfoDescription rv = MetaInfoDescription.__new__(MetaInfoDescription)
       rv.inst = shared_ptr[_MetaInfoDescription](new _MetaInfoDescription(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaInfoDescription rv = MetaInfoDescription.__new__(MetaInfoDescription)
       rv.inst = shared_ptr[_MetaInfoDescription](new _MetaInfoDescription(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaInfoDescription](new _MetaInfoDescription())
    
    def _init_1(self, MetaInfoDescription in_0 ):
        """
        _init_1(self, in_0: MetaInfoDescription ) -> None
        """
        assert isinstance(in_0, MetaInfoDescription), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaInfoDescription](new _MetaInfoDescription((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaInfoDescription ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaInfoDescription)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
        if not isinstance(other, MetaInfoDescription):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef MetaInfoDescription other_casted = other
        cdef MetaInfoDescription self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class MetaboTargetedAssay:
    """
    Cython implementation of _MetaboTargetedAssay

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedAssay.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboTargetedAssay rv = MetaboTargetedAssay.__new__(MetaboTargetedAssay)
       rv.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboTargetedAssay rv = MetaboTargetedAssay.__new__(MetaboTargetedAssay)
       rv.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        This class provides methods for the extraction of targeted assays for metabolomics
        """
        self.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay())
    
    def _init_1(self, MetaboTargetedAssay in_0 ):
        """
        _init_1(self, in_0: MetaboTargetedAssay ) -> None
        """
        assert isinstance(in_0, MetaboTargetedAssay), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        This class provides methods for the extraction of targeted assays for metabolomics

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboTargetedAssay ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboTargetedAssay)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractMetaboTargetedAssay(self, MSExperiment spectra , FeatureMapping_FeatureToMs2Indices feature_ms2_index , double precursor_rt_tol , double precursor_mz_distance , double cosine_sim_threshold , double transition_threshold , double min_fragment_mz , double max_fragment_mz , bool method_consensus_spectrum , bool exclude_ms2_precursor ,  file_counter ):
        """
        extractMetaboTargetedAssay(self, spectra: MSExperiment , feature_ms2_index: FeatureMapping_FeatureToMs2Indices , precursor_rt_tol: float , precursor_mz_distance: float , cosine_sim_threshold: float , transition_threshold: float , min_fragment_mz: float , max_fragment_mz: float , method_consensus_spectrum: bool , exclude_ms2_precursor: bool , file_counter: int ) -> List[MetaboTargetedAssay]
        Extract a vector of MetaboTargetedAssays without using fragment annotation
        
        
        :param spectra: Input of MSExperiment with spectra information
        :param feature_ms2_spectra_map: FeatureMapping class with associated MS2 spectra
        :param precursor_rt_tol: Retention time tolerance of the precursor
        :param precursor_mz_distance: Max m/z distance of the precursor entries of two spectra to be merged
        :param cosine_sim_threshold: Cosine similarty threshold for the usage of SpectraMerger
        :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
        :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
        :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
        :param method_consensus_spectrum: Boolean to use consensus spectrum method
        :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
        :return: Vector of MetaboTargetedAssay
        """
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(feature_ms2_index, FeatureMapping_FeatureToMs2Indices), 'arg feature_ms2_index wrong type'
        assert isinstance(precursor_rt_tol, float), 'arg precursor_rt_tol wrong type'
        assert isinstance(precursor_mz_distance, float), 'arg precursor_mz_distance wrong type'
        assert isinstance(cosine_sim_threshold, float), 'arg cosine_sim_threshold wrong type'
        assert isinstance(transition_threshold, float), 'arg transition_threshold wrong type'
        assert isinstance(min_fragment_mz, float), 'arg min_fragment_mz wrong type'
        assert isinstance(max_fragment_mz, float), 'arg max_fragment_mz wrong type'
        assert isinstance(method_consensus_spectrum, pybool_t), 'arg method_consensus_spectrum wrong type'
        assert isinstance(exclude_ms2_precursor, pybool_t), 'arg exclude_ms2_precursor wrong type'
        assert isinstance(file_counter, int), 'arg file_counter wrong type'
    
    
    
    
    
    
    
    
    
    
    
        _r = self.inst.get().extractMetaboTargetedAssay((deref(spectra.inst.get())), (deref(feature_ms2_index.inst.get())), (<double &>precursor_rt_tol), (<double &>precursor_mz_distance), (<double &>cosine_sim_threshold), (<double &>transition_threshold), (<double &>min_fragment_mz), (<double &>max_fragment_mz), (<bool &>method_consensus_spectrum), (<bool &>exclude_ms2_precursor), (<unsigned int &>file_counter))
        py_result = []
        cdef libcpp_vector[_MetaboTargetedAssay].iterator it__r = _r.begin()
        cdef MetaboTargetedAssay item_py_result
        while it__r != _r.end():
           item_py_result = MetaboTargetedAssay.__new__(MetaboTargetedAssay)
           item_py_result.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def extractMetaboTargetedAssayFragmentAnnotation(self, list v_cmp_spec , double transition_threshold , double min_fragment_mz , double max_fragment_mz , bool use_exact_mass , bool exclude_ms2_precursor ):
        """
        extractMetaboTargetedAssayFragmentAnnotation(self, v_cmp_spec: List[MetaboTargetedAssay_CompoundTargetDecoyPair] , transition_threshold: float , min_fragment_mz: float , max_fragment_mz: float , use_exact_mass: bool , exclude_ms2_precursor: bool ) -> List[MetaboTargetedAssay]
        Extract a vector of MetaboTargetedAssays using fragment
        
        
        :param v_cmp_spec: Vector of CompoundInfo with associated fragment annotated MSspectrum
        :param transition_threshold: Intensity threshold for MS2 peak used in MetaboTargetedAssay
        :param min_fragment_mz: Minimum m/z a fragment ion has to have to be considered as a transition
        :param max_fragment_mz: Maximum m/z a fragment ion has to have to be considered as a transition
        :param use_exact_mass: Boolean if exact mass should be used as peak mass for annotated fragments
        :param exclude_ms2_precursor: Boolean to exclude MS2 precursor from MetaboTargetedAssay
        :param file_counter: Count if multiple files are used.
        :return: Vector of MetaboTargetedAssay
        """
        assert isinstance(v_cmp_spec, list) and all(isinstance(elemt_rec, MetaboTargetedAssay_CompoundTargetDecoyPair) for elemt_rec in v_cmp_spec), 'arg v_cmp_spec wrong type'
        assert isinstance(transition_threshold, float), 'arg transition_threshold wrong type'
        assert isinstance(min_fragment_mz, float), 'arg min_fragment_mz wrong type'
        assert isinstance(max_fragment_mz, float), 'arg max_fragment_mz wrong type'
        assert isinstance(use_exact_mass, pybool_t), 'arg use_exact_mass wrong type'
        assert isinstance(exclude_ms2_precursor, pybool_t), 'arg exclude_ms2_precursor wrong type'
        cdef libcpp_vector[_MetaboTargetedAssay_CompoundTargetDecoyPair] * v0 = new libcpp_vector[_MetaboTargetedAssay_CompoundTargetDecoyPair]()
        cdef MetaboTargetedAssay_CompoundTargetDecoyPair item0
        for item0 in v_cmp_spec:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
    
        _r = self.inst.get().extractMetaboTargetedAssayFragmentAnnotation(deref(v0), (<double &>transition_threshold), (<double &>min_fragment_mz), (<double &>max_fragment_mz), (<bool &>use_exact_mass), (<bool &>exclude_ms2_precursor))
        cdef libcpp_vector[_MetaboTargetedAssay_CompoundTargetDecoyPair].iterator it_v_cmp_spec = v0.begin()
        replace_0 = []
        while it_v_cmp_spec != v0.end():
            item0 = MetaboTargetedAssay_CompoundTargetDecoyPair.__new__(MetaboTargetedAssay_CompoundTargetDecoyPair)
            item0.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair(deref(it_v_cmp_spec)))
            replace_0.append(item0)
            inc(it_v_cmp_spec)
        v_cmp_spec[:] = replace_0
        del v0
        py_result = []
        cdef libcpp_vector[_MetaboTargetedAssay].iterator it__r = _r.begin()
        cdef MetaboTargetedAssay item_py_result
        while it__r != _r.end():
           item_py_result = MetaboTargetedAssay.__new__(MetaboTargetedAssay)
           item_py_result.inst = shared_ptr[_MetaboTargetedAssay](new _MetaboTargetedAssay(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def pairCompoundWithAnnotatedTDSpectraPairs(self, list v_cmpinfo , list annotated_spectra ):
        """
        pairCompoundWithAnnotatedTDSpectraPairs(self, v_cmpinfo: List[SiriusMSFile_CompoundInfo] , annotated_spectra: List[SiriusFragmentAnnotation_SiriusTargetDecoySpectra] ) -> List[MetaboTargetedAssay_CompoundTargetDecoyPair]
        Pair compound information (SiriusMSFile) with the annotated target and decoy spectrum from SIRIUS/Passatutto based on the m_id (unique identifier composed of description_filepath_native_id_k introduced in the SiriusMSConverter)
        
        
        :param v_cmpinfo: Vector of SiriusMSFile::CompoundInfo
        :param annotated_spectra: Vector of SiriusTargetDecoySpectra
        :return: Vector of MetaboTargetedAssay::CompoundTargetDecoyPair
        """
        assert isinstance(v_cmpinfo, list) and all(isinstance(elemt_rec, SiriusMSFile_CompoundInfo) for elemt_rec in v_cmpinfo), 'arg v_cmpinfo wrong type'
        assert isinstance(annotated_spectra, list) and all(isinstance(elemt_rec, SiriusFragmentAnnotation_SiriusTargetDecoySpectra) for elemt_rec in annotated_spectra), 'arg annotated_spectra wrong type'
        cdef libcpp_vector[_SiriusMSFile_CompoundInfo] * v0 = new libcpp_vector[_SiriusMSFile_CompoundInfo]()
        cdef SiriusMSFile_CompoundInfo item0
        for item0 in v_cmpinfo:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra] * v1 = new libcpp_vector[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra]()
        cdef SiriusFragmentAnnotation_SiriusTargetDecoySpectra item1
        for item1 in annotated_spectra:
            v1.push_back(deref(item1.inst.get()))
        _r = self.inst.get().pairCompoundWithAnnotatedTDSpectraPairs(deref(v0), deref(v1))
        cdef libcpp_vector[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra].iterator it_annotated_spectra = v1.begin()
        replace_0 = []
        while it_annotated_spectra != v1.end():
            item1 = SiriusFragmentAnnotation_SiriusTargetDecoySpectra.__new__(SiriusFragmentAnnotation_SiriusTargetDecoySpectra)
            item1.inst = shared_ptr[_SiriusFragmentAnnotation_SiriusTargetDecoySpectra](new _SiriusFragmentAnnotation_SiriusTargetDecoySpectra(deref(it_annotated_spectra)))
            replace_0.append(item1)
            inc(it_annotated_spectra)
        annotated_spectra[:] = replace_0
        del v1
        cdef libcpp_vector[_SiriusMSFile_CompoundInfo].iterator it_v_cmpinfo = v0.begin()
        replace_0 = []
        while it_v_cmpinfo != v0.end():
            item0 = SiriusMSFile_CompoundInfo.__new__(SiriusMSFile_CompoundInfo)
            item0.inst = shared_ptr[_SiriusMSFile_CompoundInfo](new _SiriusMSFile_CompoundInfo(deref(it_v_cmpinfo)))
            replace_0.append(item0)
            inc(it_v_cmpinfo)
        v_cmpinfo[:] = replace_0
        del v0
        py_result = []
        cdef libcpp_vector[_MetaboTargetedAssay_CompoundTargetDecoyPair].iterator it__r = _r.begin()
        cdef MetaboTargetedAssay_CompoundTargetDecoyPair item_py_result
        while it__r != _r.end():
           item_py_result = MetaboTargetedAssay_CompoundTargetDecoyPair.__new__(MetaboTargetedAssay_CompoundTargetDecoyPair)
           item_py_result.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result 

cdef class MetaboTargetedAssay_CompoundTargetDecoyPair:
    """
    Cython implementation of _MetaboTargetedAssay_CompoundTargetDecoyPair

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboTargetedAssay_CompoundTargetDecoyPair.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboTargetedAssay_CompoundTargetDecoyPair rv = MetaboTargetedAssay_CompoundTargetDecoyPair.__new__(MetaboTargetedAssay_CompoundTargetDecoyPair)
       rv.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboTargetedAssay_CompoundTargetDecoyPair rv = MetaboTargetedAssay_CompoundTargetDecoyPair.__new__(MetaboTargetedAssay_CompoundTargetDecoyPair)
       rv.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair())
    
    def _init_1(self, MetaboTargetedAssay_CompoundTargetDecoyPair in_0 ):
        """
        _init_1(self, in_0: MetaboTargetedAssay_CompoundTargetDecoyPair ) -> None
        """
        assert isinstance(in_0, MetaboTargetedAssay_CompoundTargetDecoyPair), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboTargetedAssay_CompoundTargetDecoyPair](new _MetaboTargetedAssay_CompoundTargetDecoyPair((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboTargetedAssay_CompoundTargetDecoyPair ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboTargetedAssay_CompoundTargetDecoyPair)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MorpheusScore:
    """
    Cython implementation of _MorpheusScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorpheusScore.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MorpheusScore rv = MorpheusScore.__new__(MorpheusScore)
       rv.inst = shared_ptr[_MorpheusScore](new _MorpheusScore(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MorpheusScore rv = MorpheusScore.__new__(MorpheusScore)
       rv.inst = shared_ptr[_MorpheusScore](new _MorpheusScore(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MorpheusScore](new _MorpheusScore())
    
    def _init_1(self, MorpheusScore in_0 ):
        """
        _init_1(self, in_0: MorpheusScore ) -> None
        """
        assert isinstance(in_0, MorpheusScore), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MorpheusScore](new _MorpheusScore((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MorpheusScore ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MorpheusScore)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, double fragment_mass_tolerance , bool fragment_mass_tolerance_unit_ppm , MSSpectrum exp_spectrum , MSSpectrum theo_spectrum ):
        """
        compute(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , exp_spectrum: MSSpectrum , theo_spectrum: MSSpectrum ) -> MorpheusScore_Result
        Returns Morpheus Score
        """
        assert isinstance(fragment_mass_tolerance, float), 'arg fragment_mass_tolerance wrong type'
        assert isinstance(fragment_mass_tolerance_unit_ppm, pybool_t), 'arg fragment_mass_tolerance_unit_ppm wrong type'
        assert isinstance(exp_spectrum, MSSpectrum), 'arg exp_spectrum wrong type'
        assert isinstance(theo_spectrum, MSSpectrum), 'arg theo_spectrum wrong type'
    
    
    
    
        cdef _MorpheusScore_Result * _r = new _MorpheusScore_Result(self.inst.get().compute((<double>fragment_mass_tolerance), (<bool>fragment_mass_tolerance_unit_ppm), (deref(exp_spectrum.inst.get())), (deref(theo_spectrum.inst.get()))))
        cdef MorpheusScore_Result py_result = MorpheusScore_Result.__new__(MorpheusScore_Result)
        py_result.inst = shared_ptr[_MorpheusScore_Result](_r)
        return py_result 

cdef class MorpheusScore_Result:
    """
    Cython implementation of _MorpheusScore_Result

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MorpheusScore_Result.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property matches:
        def __set__(self,  matches):
        
            self.inst.get().matches = (<size_t>matches)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().matches
            py_result = <size_t>_r
            return py_result
    
    property n_peaks:
        def __set__(self,  n_peaks):
        
            self.inst.get().n_peaks = (<size_t>n_peaks)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().n_peaks
            py_result = <size_t>_r
            return py_result
    
    property score:
        def __set__(self, float score):
        
            self.inst.get().score = (<float>score)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().score
            py_result = <float>_r
            return py_result
    
    property MIC:
        def __set__(self, float MIC):
        
            self.inst.get().MIC = (<float>MIC)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().MIC
            py_result = <float>_r
            return py_result
    
    property TIC:
        def __set__(self, float TIC):
        
            self.inst.get().TIC = (<float>TIC)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().TIC
            py_result = <float>_r
            return py_result
    
    property err:
        def __set__(self, float err):
        
            self.inst.get().err = (<float>err)
        
    
        def __get__(self):
            cdef float _r = self.inst.get().err
            py_result = <float>_r
            return py_result
    
    def __copy__(self):
       cdef MorpheusScore_Result rv = MorpheusScore_Result.__new__(MorpheusScore_Result)
       rv.inst = shared_ptr[_MorpheusScore_Result](new _MorpheusScore_Result(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MorpheusScore_Result rv = MorpheusScore_Result.__new__(MorpheusScore_Result)
       rv.inst = shared_ptr[_MorpheusScore_Result](new _MorpheusScore_Result(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MorpheusScore_Result](new _MorpheusScore_Result())
    
    def _init_1(self, MorpheusScore_Result in_0 ):
        """
        _init_1(self, in_0: MorpheusScore_Result ) -> None
        """
        assert isinstance(in_0, MorpheusScore_Result), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MorpheusScore_Result](new _MorpheusScore_Result((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MorpheusScore_Result ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MorpheusScore_Result)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class MzTabMFile:
    """
    Cython implementation of _MzTabMFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabMFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzTabMFile rv = MzTabMFile.__new__(MzTabMFile)
       rv.inst = shared_ptr[_MzTabMFile](new _MzTabMFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzTabMFile rv = MzTabMFile.__new__(MzTabMFile)
       rv.inst = shared_ptr[_MzTabMFile](new _MzTabMFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzTabMFile](new _MzTabMFile())
    
    def _init_1(self, MzTabMFile in_0 ):
        """
        _init_1(self, in_0: MzTabMFile ) -> None
        """
        assert isinstance(in_0, MzTabMFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzTabMFile](new _MzTabMFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzTabMFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzTabMFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def store(self,  filename , MzTabM mztab_m ):
        """
        store(self, filename: Union[bytes, str, String] , mztab_m: MzTabM ) -> None
        Store MzTabM file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(mztab_m, MzTabM), 'arg mztab_m wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(mztab_m.inst.get()))) 

cdef class ParamCTDFile:
    """
    Cython implementation of _ParamCTDFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamCTDFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_ParamCTDFile](new _ParamCTDFile())
    
    def store(self,  filename , Param param , ToolInfo tool_info ):
        """
        store(self, filename: Union[bytes, str] , param: Param , tool_info: ToolInfo ) -> None
        """
        assert isinstance(filename, (bytes, str)), 'arg filename wrong type'
        assert isinstance(param, Param), 'arg param wrong type'
        assert isinstance(tool_info, ToolInfo), 'arg tool_info wrong type'
        if isinstance(filename, str):
            filename = filename.encode('utf-8')
    
    
        self.inst.get().store((<libcpp_string>filename), (deref(param.inst.get())), (deref(tool_info.inst.get()))) 

cdef class PercolatorInfile:
    """
    Cython implementation of _PercolatorInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PercolatorInfile.html>`_

    Class for storing Percolator tab-delimited input files
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PercolatorInfile rv = PercolatorInfile.__new__(PercolatorInfile)
       rv.inst = shared_ptr[_PercolatorInfile](new _PercolatorInfile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PercolatorInfile rv = PercolatorInfile.__new__(PercolatorInfile)
       rv.inst = shared_ptr[_PercolatorInfile](new _PercolatorInfile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PercolatorInfile](new _PercolatorInfile())
    
    def _init_1(self, PercolatorInfile in_0 ):
        """
        _init_1(self, in_0: PercolatorInfile ) -> None
        """
        assert isinstance(in_0, PercolatorInfile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PercolatorInfile](new _PercolatorInfile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PercolatorInfile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PercolatorInfile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    store = __static_PercolatorInfile_store 

cdef class SiriusExportAlgorithm:
    """
    Cython implementation of _SiriusExportAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SiriusExportAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SiriusExportAlgorithm rv = SiriusExportAlgorithm.__new__(SiriusExportAlgorithm)
       rv.inst = shared_ptr[_SiriusExportAlgorithm](new _SiriusExportAlgorithm(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SiriusExportAlgorithm rv = SiriusExportAlgorithm.__new__(SiriusExportAlgorithm)
       rv.inst = shared_ptr[_SiriusExportAlgorithm](new _SiriusExportAlgorithm(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SiriusExportAlgorithm](new _SiriusExportAlgorithm())
    
    def _init_1(self, SiriusExportAlgorithm in_0 ):
        """
        _init_1(self, in_0: SiriusExportAlgorithm ) -> None
        """
        assert isinstance(in_0, SiriusExportAlgorithm), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SiriusExportAlgorithm](new _SiriusExportAlgorithm((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SiriusExportAlgorithm ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SiriusExportAlgorithm)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def isFeatureOnly(self):
        """
        isFeatureOnly(self) -> bool
        """
        cdef bool _r = self.inst.get().isFeatureOnly()
        py_result = <bool>_r
        return py_result
    
    def getFilterByNumMassTraces(self):
        """
        getFilterByNumMassTraces(self) -> int
        """
        cdef unsigned int _r = self.inst.get().getFilterByNumMassTraces()
        py_result = <unsigned int>_r
        return py_result
    
    def getPrecursorMzTolerance(self):
        """
        getPrecursorMzTolerance(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorMzTolerance()
        py_result = <double>_r
        return py_result
    
    def getPrecursorRtTolerance(self):
        """
        getPrecursorRtTolerance(self) -> float
        """
        cdef double _r = self.inst.get().getPrecursorRtTolerance()
        py_result = <double>_r
        return py_result
    
    def precursorMzToleranceUnitIsPPM(self):
        """
        precursorMzToleranceUnitIsPPM(self) -> bool
        """
        cdef bool _r = self.inst.get().precursorMzToleranceUnitIsPPM()
        py_result = <bool>_r
        return py_result
    
    def isNoMasstraceInfoIsotopePattern(self):
        """
        isNoMasstraceInfoIsotopePattern(self) -> bool
        """
        cdef bool _r = self.inst.get().isNoMasstraceInfoIsotopePattern()
        py_result = <bool>_r
        return py_result
    
    def getIsotopePatternIterations(self):
        """
        getIsotopePatternIterations(self) -> int
        """
        cdef int _r = self.inst.get().getIsotopePatternIterations()
        py_result = <int>_r
        return py_result
    
    def preprocessing(self,  featureXML_path , MSExperiment spectra , FeatureMapping_FeatureMappingInfo feature_mapping_info , FeatureMapping_FeatureToMs2Indices feature_ms2_indices ):
        """
        preprocessing(self, featureXML_path: Union[bytes, str, String] , spectra: MSExperiment , feature_mapping_info: FeatureMapping_FeatureMappingInfo , feature_ms2_indices: FeatureMapping_FeatureToMs2Indices ) -> None
        Preprocessing needed for SIRIUS
        
        Filter number of masstraces and perform feature mapping
        
        :param featureXML_path: Path to featureXML
        :param spectra: Input of MSExperiment with spectra information
        :param feature_mapping_info: Emtpy - stores FeatureMaps and KDTreeMaps internally
        :param feature_ms2_indices: Empty FeatureToMs2Indices
        """
        assert (isinstance(featureXML_path, str) or isinstance(featureXML_path, bytes) or isinstance(featureXML_path, String)), 'arg featureXML_path wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert isinstance(feature_mapping_info, FeatureMapping_FeatureMappingInfo), 'arg feature_mapping_info wrong type'
        assert isinstance(feature_ms2_indices, FeatureMapping_FeatureToMs2Indices), 'arg feature_ms2_indices wrong type'
    
    
    
    
        self.inst.get().preprocessing(deref((convString(featureXML_path)).get()), (deref(spectra.inst.get())), (deref(feature_mapping_info.inst.get())), (deref(feature_ms2_indices.inst.get())))
    
    def logFeatureSpectraNumber(self,  featureXML_path , FeatureMapping_FeatureToMs2Indices feature_ms2_indices , MSExperiment spectra ):
        """
        logFeatureSpectraNumber(self, featureXML_path: Union[bytes, str, String] , feature_ms2_indices: FeatureMapping_FeatureToMs2Indices , spectra: MSExperiment ) -> None
        Logs number of features and spectra used
        
        Prints the number of features and spectra used (OPENMS_LOG_INFO)
        
        :param featureXML_path: Path to featureXML
        :param feature_ms2_indices: FeatureToMs2Indices with feature mapping
        :param spectra: Input of MSExperiment with spectra information
        """
        assert (isinstance(featureXML_path, str) or isinstance(featureXML_path, bytes) or isinstance(featureXML_path, String)), 'arg featureXML_path wrong type'
        assert isinstance(feature_ms2_indices, FeatureMapping_FeatureToMs2Indices), 'arg feature_ms2_indices wrong type'
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
    
    
    
        self.inst.get().logFeatureSpectraNumber(deref((convString(featureXML_path)).get()), (deref(feature_ms2_indices.inst.get())), (deref(spectra.inst.get())))
    
    def run(self, list mzML_files , list featureXML_files ,  out_ms ,  out_compoundinfo ):
        """
        run(self, mzML_files: List[bytes] , featureXML_files: List[bytes] , out_ms: Union[bytes, str, String] , out_compoundinfo: Union[bytes, str, String] ) -> None
        Runs SiriusExport with mzML and featureXML (optional) files as input.
        
        Generates a SIRIUS .ms file and compound info table (optional).
        
        :param mzML_files: List with paths to mzML files
        :param featureXML_files: List with paths to featureXML files
        :param out_ms: Output file name for SIRIUS .ms file
        :param out_compoundinfo: Output file name for tsv file with compound info
        """
        assert isinstance(mzML_files, list) and all(isinstance(li, bytes) for li in mzML_files), 'arg mzML_files wrong type'
        assert isinstance(featureXML_files, list) and all(isinstance(li, bytes) for li in featureXML_files), 'arg featureXML_files wrong type'
        assert (isinstance(out_ms, str) or isinstance(out_ms, bytes) or isinstance(out_ms, String)), 'arg out_ms wrong type'
        assert (isinstance(out_compoundinfo, str) or isinstance(out_compoundinfo, bytes) or isinstance(out_compoundinfo, String)), 'arg out_compoundinfo wrong type'
        cdef libcpp_vector[_String] * v0 = new libcpp_vector[_String]()
        cdef bytes item0
        for item0 in mzML_files:
           v0.push_back(_String(<char *>item0))
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in featureXML_files:
           v1.push_back(_String(<char *>item1))
    
    
        self.inst.get().run(deref(v0), deref(v1), deref((convString(out_ms)).get()), deref((convString(out_compoundinfo)).get()))
        replace = []
        cdef libcpp_vector[_String].iterator it_1 = v1.begin()
        while it_1 != v1.end():
           replace.append(<char*>deref(it_1).c_str())
           inc(it_1)
        featureXML_files[:] = replace
        del v1
        replace = []
        cdef libcpp_vector[_String].iterator it_0 = v0.begin()
        while it_0 != v0.end():
           replace.append(<char*>deref(it_0).c_str())
           inc(it_0)
        mzML_files[:] = replace
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

cdef class ToolInfo:
    """
    Cython implementation of _ToolInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ToolInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, ToolInfo in_0 ):
        """
        __init__(self, in_0: ToolInfo ) -> None
        """
        assert isinstance(in_0, ToolInfo), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ToolInfo](new _ToolInfo((deref(in_0.inst.get())))) 
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
