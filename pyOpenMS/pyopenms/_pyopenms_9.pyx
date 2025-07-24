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
def __static_SpectrumMetaDataLookup_addMissingIMToPeptideIDs(PeptideIdentificationList in_0 , MSExperiment exp ):
    """
    __static_SpectrumMetaDataLookup_addMissingIMToPeptideIDs(in_0: PeptideIdentificationList , exp: MSExperiment ) -> bool
    """
    assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'


    cdef bool _r = _addMissingIMToPeptideIDs_SpectrumMetaDataLookup((deref(in_0.inst.get())), (deref(exp.inst.get())))
    py_result = <bool>_r
    return py_result

def __static_SpectrumMetaDataLookup_addMissingRTsToPeptideIDs(PeptideIdentificationList in_0 , MSExperiment exp ):
    """
    __static_SpectrumMetaDataLookup_addMissingRTsToPeptideIDs(in_0: PeptideIdentificationList , exp: MSExperiment ) -> bool
    """
    assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    assert isinstance(exp, MSExperiment), 'arg exp wrong type'


    cdef bool _r = _addMissingRTsToPeptideIDs_SpectrumMetaDataLookup((deref(in_0.inst.get())), (deref(exp.inst.get())))
    py_result = <bool>_r
    return py_result

def __static_SpectrumMetaDataLookup_addMissingSpectrumReferences(PeptideIdentificationList in_0 ,  filename , bool stop_on_error , bool override_spectra_data , bool override_spectra_references , list proteins ):
    """
    __static_SpectrumMetaDataLookup_addMissingSpectrumReferences(in_0: PeptideIdentificationList , filename: Union[bytes, str, String] , stop_on_error: bool , override_spectra_data: bool , override_spectra_references: bool , proteins: List[ProteinIdentification] ) -> bool
    """
    assert isinstance(in_0, PeptideIdentificationList), 'arg in_0 wrong type'
    assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    assert isinstance(stop_on_error, pybool_t), 'arg stop_on_error wrong type'
    assert isinstance(override_spectra_data, pybool_t), 'arg override_spectra_data wrong type'
    assert isinstance(override_spectra_references, pybool_t), 'arg override_spectra_references wrong type'
    assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'





    cdef libcpp_vector[_ProteinIdentification] * v5 = new libcpp_vector[_ProteinIdentification]()
    cdef ProteinIdentification item5
    for item5 in proteins:
        v5.push_back(deref(item5.inst.get()))
    cdef bool _r = _addMissingSpectrumReferences_SpectrumMetaDataLookup((deref(in_0.inst.get())), deref((convString(filename)).get()), (<bool>stop_on_error), (<bool>override_spectra_data), (<bool>override_spectra_references), deref(v5))
    del v5
    py_result = <bool>_r
    return py_result

def __static_SpectrumMetaDataLookup_getSpectrumMetaData(MSSpectrum spectrum , SpectrumMetaData meta ):
    """
    __static_SpectrumMetaDataLookup_getSpectrumMetaData(spectrum: MSSpectrum , meta: SpectrumMetaData ) -> None
    """
    assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    assert isinstance(meta, SpectrumMetaData), 'arg meta wrong type'


    _getSpectrumMetaData_SpectrumMetaDataLookup((deref(spectrum.inst.get())), (deref(meta.inst.get()))) 

cdef class __CHARGEMODE_MFD:
    None
    QFROMFEATURE = 0
    QHEURISTIC = 1
    QALL = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class LogType:
    None
    CMD = 0
    GUI = 1
    NONE = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __PeptidePosition:
    None
    INTERNAL = 0
    C_TERM = 1
    N_TERM = 2

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class __ProteinProteinCrossLinkType:
    None
    CROSS = 0
    MONO = 1
    LOOP = 2
    NUMBER_OF_CROSS_LINK_TYPES = 3

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class BilinearInterpolation:
    """
    Cython implementation of _BilinearInterpolation[double,double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1BilinearInterpolation[double,double].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef BilinearInterpolation rv = BilinearInterpolation.__new__(BilinearInterpolation)
       rv.inst = shared_ptr[_BilinearInterpolation[double,double]](new _BilinearInterpolation[double,double](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef BilinearInterpolation rv = BilinearInterpolation.__new__(BilinearInterpolation)
       rv.inst = shared_ptr[_BilinearInterpolation[double,double]](new _BilinearInterpolation[double,double](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_BilinearInterpolation[double,double]](new _BilinearInterpolation[double,double]())
    
    def _init_1(self, BilinearInterpolation in_0 ):
        """
        _init_1(self, in_0: BilinearInterpolation ) -> None
        """
        assert isinstance(in_0, BilinearInterpolation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_BilinearInterpolation[double,double]](new _BilinearInterpolation[double,double]((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: BilinearInterpolation ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], BilinearInterpolation)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def value(self, double arg_pos_0 , double arg_pos_1 ):
        """
        value(self, arg_pos_0: float , arg_pos_1: float ) -> float
        """
        assert isinstance(arg_pos_0, float), 'arg arg_pos_0 wrong type'
        assert isinstance(arg_pos_1, float), 'arg arg_pos_1 wrong type'
    
    
        cdef double _r = self.inst.get().value((<double>arg_pos_0), (<double>arg_pos_1))
        py_result = <double>_r
        return py_result
    
    def addValue(self, double arg_pos_0 , double arg_pos_1 , double arg_value ):
        """
        addValue(self, arg_pos_0: float , arg_pos_1: float , arg_value: float ) -> None
        Performs bilinear resampling. The arg_value is split up and added to the data points around arg_pos. ("forward resampling")
        """
        assert isinstance(arg_pos_0, float), 'arg arg_pos_0 wrong type'
        assert isinstance(arg_pos_1, float), 'arg arg_pos_1 wrong type'
        assert isinstance(arg_value, float), 'arg arg_value wrong type'
    
    
    
        self.inst.get().addValue((<double>arg_pos_0), (<double>arg_pos_1), (<double>arg_value))
    
    def getData(self):
        """
        getData(self) -> MatrixDouble
        """
        cdef _Matrix[double] * _r = new _Matrix[double](self.inst.get().getData())
        cdef MatrixDouble py_result = MatrixDouble.__new__(MatrixDouble)
        py_result.inst = shared_ptr[_Matrix[double]](_r)
        return py_result
    
    def setData(self, MatrixDouble data ):
        """
        setData(self, data: MatrixDouble ) -> None
        Assigns data to the internal random access container storing the data. SourceContainer must be assignable to ContainerType
        """
        assert isinstance(data, MatrixDouble), 'arg data wrong type'
    
        self.inst.get().setData((deref(data.inst.get())))
    
    def empty(self):
        """
        empty(self) -> bool
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def key2index_0(self, double pos ):
        """
        key2index_0(self, pos: float ) -> float
        The transformation from "outside" to "inside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().key2index_0((<double>pos))
        py_result = <double>_r
        return py_result
    
    def index2key_0(self, double pos ):
        """
        index2key_0(self, pos: float ) -> float
        The transformation from "inside" to "outside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().index2key_0((<double>pos))
        py_result = <double>_r
        return py_result
    
    def key2index_1(self, double pos ):
        """
        key2index_1(self, pos: float ) -> float
        The transformation from "outside" to "inside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().key2index_1((<double>pos))
        py_result = <double>_r
        return py_result
    
    def index2key_1(self, double pos ):
        """
        index2key_1(self, pos: float ) -> float
        The transformation from "inside" to "outside" coordinates
        """
        assert isinstance(pos, float), 'arg pos wrong type'
    
        cdef double _r = self.inst.get().index2key_1((<double>pos))
        py_result = <double>_r
        return py_result
    
    def getScale_0(self):
        """
        getScale_0(self) -> float
        """
        cdef double _r = self.inst.get().getScale_0()
        py_result = <double>_r
        return py_result
    
    def setScale_0(self, double scale ):
        """
        setScale_0(self, scale: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
    
        self.inst.get().setScale_0((<double &>scale))
    
    def getScale_1(self):
        """
        getScale_1(self) -> float
        """
        cdef double _r = self.inst.get().getScale_1()
        py_result = <double>_r
        return py_result
    
    def setScale_1(self, double scale ):
        """
        setScale_1(self, scale: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
    
        self.inst.get().setScale_1((<double &>scale))
    
    def getOffset_0(self):
        """
        getOffset_0(self) -> float
        Accessor. "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        """
        cdef double _r = self.inst.get().getOffset_0()
        py_result = <double>_r
        return py_result
    
    def setOffset_0(self, double offset ):
        """
        setOffset_0(self, offset: float ) -> None
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset_0((<double &>offset))
    
    def getOffset_1(self):
        """
        getOffset_1(self) -> float
        Accessor. "Offset" is the point (in "outside" units) which corresponds to "Data(0,0)"
        """
        cdef double _r = self.inst.get().getOffset_1()
        py_result = <double>_r
        return py_result
    
    def setOffset_1(self, double offset ):
        """
        setOffset_1(self, offset: float ) -> None
        """
        assert isinstance(offset, float), 'arg offset wrong type'
    
        self.inst.get().setOffset_1((<double &>offset))
    
    def _setMapping_0_0(self, double scale , double inside , double outside ):
        """
        _setMapping_0_0(self, scale: float , inside: float , outside: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
        assert isinstance(inside, float), 'arg inside wrong type'
        assert isinstance(outside, float), 'arg outside wrong type'
    
    
    
        self.inst.get().setMapping_0((<double &>scale), (<double &>inside), (<double &>outside))
    
    def _setMapping_0_1(self, double inside_low , double outside_low , double inside_high , double outside_high ):
        """
        _setMapping_0_1(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
        """
        assert isinstance(inside_low, float), 'arg inside_low wrong type'
        assert isinstance(outside_low, float), 'arg outside_low wrong type'
        assert isinstance(inside_high, float), 'arg inside_high wrong type'
        assert isinstance(outside_high, float), 'arg outside_high wrong type'
    
    
    
    
        self.inst.get().setMapping_0((<double &>inside_low), (<double &>outside_low), (<double &>inside_high), (<double &>outside_high))
    
    def setMapping_0(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMapping_0(self, scale: float , inside: float , outside: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMapping_0(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._setMapping_0_0(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._setMapping_0_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _setMapping_1_0(self, double scale , double inside , double outside ):
        """
        _setMapping_1_0(self, scale: float , inside: float , outside: float ) -> None
        """
        assert isinstance(scale, float), 'arg scale wrong type'
        assert isinstance(inside, float), 'arg inside wrong type'
        assert isinstance(outside, float), 'arg outside wrong type'
    
    
    
        self.inst.get().setMapping_1((<double &>scale), (<double &>inside), (<double &>outside))
    
    def _setMapping_1_1(self, double inside_low , double outside_low , double inside_high , double outside_high ):
        """
        _setMapping_1_1(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
        """
        assert isinstance(inside_low, float), 'arg inside_low wrong type'
        assert isinstance(outside_low, float), 'arg outside_low wrong type'
        assert isinstance(inside_high, float), 'arg inside_high wrong type'
        assert isinstance(outside_high, float), 'arg outside_high wrong type'
    
    
    
    
        self.inst.get().setMapping_1((<double &>inside_low), (<double &>outside_low), (<double &>inside_high), (<double &>outside_high))
    
    def setMapping_1(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: setMapping_1(self, scale: float , inside: float , outside: float ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: setMapping_1(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._setMapping_1_0(*args)
        elif (len(args)==4) and (isinstance(args[0], float)) and (isinstance(args[1], float)) and (isinstance(args[2], float)) and (isinstance(args[3], float)):
            return self._setMapping_1_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getInsideReferencePoint_0(self):
        """
        getInsideReferencePoint_0(self) -> float
        """
        cdef double _r = self.inst.get().getInsideReferencePoint_0()
        py_result = <double>_r
        return py_result
    
    def getInsideReferencePoint_1(self):
        """
        getInsideReferencePoint_1(self) -> float
        """
        cdef double _r = self.inst.get().getInsideReferencePoint_1()
        py_result = <double>_r
        return py_result
    
    def getOutsideReferencePoint_0(self):
        """
        getOutsideReferencePoint_0(self) -> float
        """
        cdef double _r = self.inst.get().getOutsideReferencePoint_0()
        py_result = <double>_r
        return py_result
    
    def getOutsideReferencePoint_1(self):
        """
        getOutsideReferencePoint_1(self) -> float
        """
        cdef double _r = self.inst.get().getOutsideReferencePoint_1()
        py_result = <double>_r
        return py_result
    
    def supportMin_0(self):
        """
        supportMin_0(self) -> float
        Lower boundary of the support, in "outside" coordinates
        """
        cdef double _r = self.inst.get().supportMin_0()
        py_result = <double>_r
        return py_result
    
    def supportMin_1(self):
        """
        supportMin_1(self) -> float
        Lower boundary of the support, in "outside" coordinates
        """
        cdef double _r = self.inst.get().supportMin_1()
        py_result = <double>_r
        return py_result
    
    def supportMax_0(self):
        """
        supportMax_0(self) -> float
        Upper boundary of the support, in "outside" coordinates
        """
        cdef double _r = self.inst.get().supportMax_0()
        py_result = <double>_r
        return py_result
    
    def supportMax_1(self):
        """
        supportMax_1(self) -> float
        Upper boundary of the support, in "outside" coordinates
        """
        cdef double _r = self.inst.get().supportMax_1()
        py_result = <double>_r
        return py_result 

cdef class CachedMzMLHandler:
    """
    Cython implementation of _CachedMzMLHandler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1CachedMzMLHandler.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CachedMzMLHandler rv = CachedMzMLHandler.__new__(CachedMzMLHandler)
       rv.inst = shared_ptr[_CachedMzMLHandler](new _CachedMzMLHandler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CachedMzMLHandler rv = CachedMzMLHandler.__new__(CachedMzMLHandler)
       rv.inst = shared_ptr[_CachedMzMLHandler](new _CachedMzMLHandler(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        An internal class that handles single spectra and chromatograms
        """
        self.inst = shared_ptr[_CachedMzMLHandler](new _CachedMzMLHandler())
    
    def _init_1(self, CachedMzMLHandler in_0 ):
        """
        _init_1(self, in_0: CachedMzMLHandler ) -> None
        """
        assert isinstance(in_0, CachedMzMLHandler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CachedMzMLHandler](new _CachedMzMLHandler((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        An internal class that handles single spectra and chromatograms

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CachedMzMLHandler ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CachedMzMLHandler)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def writeMemdump(self, MSExperiment exp ,  out ):
        """
        writeMemdump(self, exp: MSExperiment , out: Union[bytes, str, String] ) -> None
        Write complete spectra as a dump to the disk
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(out, str) or isinstance(out, bytes) or isinstance(out, String)), 'arg out wrong type'
    
    
        self.inst.get().writeMemdump((deref(exp.inst.get())), deref((convString(out)).get()))
    
    def writeMetadata(self, MSExperiment exp ,  out_meta ):
        """
        writeMetadata(self, exp: MSExperiment , out_meta: Union[bytes, str, String] ) -> None
        Write only the meta data of an MSExperiment
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(out_meta, str) or isinstance(out_meta, bytes) or isinstance(out_meta, String)), 'arg out_meta wrong type'
    
    
        self.inst.get().writeMetadata((deref(exp.inst.get())), deref((convString(out_meta)).get()))
    
    def readMemdump(self, MSExperiment exp ,  filename ):
        """
        readMemdump(self, exp: MSExperiment , filename: Union[bytes, str, String] ) -> None
        Read all spectra from a dump from the disk
        """
        assert isinstance(exp, MSExperiment), 'arg exp wrong type'
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
    
        self.inst.get().readMemdump((deref(exp.inst.get())), deref((convString(filename)).get()))
    
    def getSpectraIndex(self):
        """
        getSpectraIndex(self) -> List[streampos]
        """
        _r = self.inst.get().getSpectraIndex()
        py_result = []
        cdef libcpp_vector[_streampos].iterator it__r = _r.begin()
        cdef streampos item_py_result
        while it__r != _r.end():
           item_py_result = streampos.__new__(streampos)
           item_py_result.inst = shared_ptr[_streampos](new _streampos(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def getChromatogramIndex(self):
        """
        getChromatogramIndex(self) -> List[streampos]
        """
        _r = self.inst.get().getChromatogramIndex()
        py_result = []
        cdef libcpp_vector[_streampos].iterator it__r = _r.begin()
        cdef streampos item_py_result
        while it__r != _r.end():
           item_py_result = streampos.__new__(streampos)
           item_py_result.inst = shared_ptr[_streampos](new _streampos(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def createMemdumpIndex(self,  filename ):
        """
        createMemdumpIndex(self, filename: Union[bytes, str, String] ) -> None
        Create an index on the location of all the spectra and chromatograms
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
    
        self.inst.get().createMemdumpIndex(deref((convString(filename)).get()))
    
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

cdef class IndexedMzMLDecoder:
    """
    Cython implementation of _IndexedMzMLDecoder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IndexedMzMLDecoder.html>`_

    A class to analyze indexedmzML files and extract the offsets of individual tags
    
    Specifically, this class allows one to extract the offsets of the <indexList>
    tag and of all <spectrum> and <chromatogram> tag using the indices found at
    the end of the indexedmzML XML structure
    
    While findIndexListOffset tries extracts the offset of the indexList tag from
    the last 1024 bytes of the file, this offset allows the function parseOffsets
    to extract all elements contained in the <indexList> tag and thus get access
    to all spectra and chromatogram offsets
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IndexedMzMLDecoder rv = IndexedMzMLDecoder.__new__(IndexedMzMLDecoder)
       rv.inst = shared_ptr[_IndexedMzMLDecoder](new _IndexedMzMLDecoder(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IndexedMzMLDecoder rv = IndexedMzMLDecoder.__new__(IndexedMzMLDecoder)
       rv.inst = shared_ptr[_IndexedMzMLDecoder](new _IndexedMzMLDecoder(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IndexedMzMLDecoder](new _IndexedMzMLDecoder())
    
    def _init_1(self, IndexedMzMLDecoder in_0 ):
        """
        _init_1(self, in_0: IndexedMzMLDecoder ) -> None
        """
        assert isinstance(in_0, IndexedMzMLDecoder), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IndexedMzMLDecoder](new _IndexedMzMLDecoder((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IndexedMzMLDecoder ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IndexedMzMLDecoder)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def findIndexListOffset(self,  in_ ,  buffersize ):
        """
        findIndexListOffset(self, in_: Union[bytes, str, String] , buffersize: int ) -> streampos
        Tries to extract the indexList offset from an indexedmzML\n
        
        This function reads by default the last few (1024) bytes of the given
        input file and tries to read the content of the <indexListOffset> tag
        The idea is that somewhere in the last parts of the file specified by the
        input string, the string <indexListOffset>xxx</indexListOffset> occurs
        This function returns the xxx part converted to an integer\n
        
        Since this function cannot determine where it will start reading
        the XML, no regular XML parser can be used for this. Therefore it uses
        regex to do its job. It matches the <indexListOffset> part and any
        numerical characters that follow
        
        
        :param in: Filename of the input indexedmzML file
        :param buffersize: How many bytes of the input file should be searched for the tag
        :return: A positive integer containing the content of the indexListOffset tag, returns -1 in case of failure no tag was found (you can re-try with a larger buffersize but most likely its not an indexed mzML). Using -1 is what the reference docu recommends: http://en.cppreference.com/w/cpp/io/streamoff
        :raises:
          Exception: FileNotFound is thrown if file cannot be found
        :raises:
          Exception: ParseError if offset cannot be parsed
        """
        assert (isinstance(in_, str) or isinstance(in_, bytes) or isinstance(in_, String)), 'arg in_ wrong type'
        assert isinstance(buffersize, int), 'arg buffersize wrong type'
    
    
        cdef _streampos * _r = new _streampos(self.inst.get().findIndexListOffset(deref((convString(in_)).get()), (<int>buffersize)))
        cdef streampos py_result = streampos.__new__(streampos)
        py_result.inst = shared_ptr[_streampos](_r)
        return py_result 

cdef class IsobaricIsotopeCorrector:
    """
    Cython implementation of _IsobaricIsotopeCorrector

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricIsotopeCorrector.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IsobaricIsotopeCorrector rv = IsobaricIsotopeCorrector.__new__(IsobaricIsotopeCorrector)
       rv.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricIsotopeCorrector rv = IsobaricIsotopeCorrector.__new__(IsobaricIsotopeCorrector)
       rv.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector())
    
    def _init_1(self, IsobaricIsotopeCorrector in_0 ):
        """
        _init_1(self, in_0: IsobaricIsotopeCorrector ) -> None
        """
        assert isinstance(in_0, IsobaricIsotopeCorrector), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricIsotopeCorrector](new _IsobaricIsotopeCorrector((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricIsotopeCorrector ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsobaricIsotopeCorrector)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _correctIsotopicImpurities_0(self, ConsensusMap consensus_map_in , ConsensusMap consensus_map_out , ItraqEightPlexQuantitationMethod quant_method ):
        """
        _correctIsotopicImpurities_0(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqEightPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
        """
        assert isinstance(consensus_map_in, ConsensusMap), 'arg consensus_map_in wrong type'
        assert isinstance(consensus_map_out, ConsensusMap), 'arg consensus_map_out wrong type'
        assert isinstance(quant_method, ItraqEightPlexQuantitationMethod), 'arg quant_method wrong type'
    
    
    
        cdef _IsobaricQuantifierStatistics * _r = new _IsobaricQuantifierStatistics(self.inst.get().correctIsotopicImpurities((deref(consensus_map_in.inst.get())), (deref(consensus_map_out.inst.get())), (quant_method.inst.get())))
        cdef IsobaricQuantifierStatistics py_result = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
        py_result.inst = shared_ptr[_IsobaricQuantifierStatistics](_r)
        return py_result
    
    def _correctIsotopicImpurities_1(self, ConsensusMap consensus_map_in , ConsensusMap consensus_map_out , ItraqFourPlexQuantitationMethod quant_method ):
        """
        _correctIsotopicImpurities_1(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqFourPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
        """
        assert isinstance(consensus_map_in, ConsensusMap), 'arg consensus_map_in wrong type'
        assert isinstance(consensus_map_out, ConsensusMap), 'arg consensus_map_out wrong type'
        assert isinstance(quant_method, ItraqFourPlexQuantitationMethod), 'arg quant_method wrong type'
    
    
    
        cdef _IsobaricQuantifierStatistics * _r = new _IsobaricQuantifierStatistics(self.inst.get().correctIsotopicImpurities((deref(consensus_map_in.inst.get())), (deref(consensus_map_out.inst.get())), (quant_method.inst.get())))
        cdef IsobaricQuantifierStatistics py_result = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
        py_result.inst = shared_ptr[_IsobaricQuantifierStatistics](_r)
        return py_result
    
    def _correctIsotopicImpurities_2(self, ConsensusMap consensus_map_in , ConsensusMap consensus_map_out , TMTSixPlexQuantitationMethod quant_method ):
        """
        _correctIsotopicImpurities_2(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTSixPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
        """
        assert isinstance(consensus_map_in, ConsensusMap), 'arg consensus_map_in wrong type'
        assert isinstance(consensus_map_out, ConsensusMap), 'arg consensus_map_out wrong type'
        assert isinstance(quant_method, TMTSixPlexQuantitationMethod), 'arg quant_method wrong type'
    
    
    
        cdef _IsobaricQuantifierStatistics * _r = new _IsobaricQuantifierStatistics(self.inst.get().correctIsotopicImpurities((deref(consensus_map_in.inst.get())), (deref(consensus_map_out.inst.get())), (quant_method.inst.get())))
        cdef IsobaricQuantifierStatistics py_result = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
        py_result.inst = shared_ptr[_IsobaricQuantifierStatistics](_r)
        return py_result
    
    def _correctIsotopicImpurities_3(self, ConsensusMap consensus_map_in , ConsensusMap consensus_map_out , TMTTenPlexQuantitationMethod quant_method ):
        """
        _correctIsotopicImpurities_3(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTTenPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
        """
        assert isinstance(consensus_map_in, ConsensusMap), 'arg consensus_map_in wrong type'
        assert isinstance(consensus_map_out, ConsensusMap), 'arg consensus_map_out wrong type'
        assert isinstance(quant_method, TMTTenPlexQuantitationMethod), 'arg quant_method wrong type'
    
    
    
        cdef _IsobaricQuantifierStatistics * _r = new _IsobaricQuantifierStatistics(self.inst.get().correctIsotopicImpurities((deref(consensus_map_in.inst.get())), (deref(consensus_map_out.inst.get())), (quant_method.inst.get())))
        cdef IsobaricQuantifierStatistics py_result = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
        py_result.inst = shared_ptr[_IsobaricQuantifierStatistics](_r)
        return py_result
    
    def correctIsotopicImpurities(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqEightPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: ItraqFourPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTSixPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: correctIsotopicImpurities(self, consensus_map_in: ConsensusMap , consensus_map_out: ConsensusMap , quant_method: TMTTenPlexQuantitationMethod ) -> IsobaricQuantifierStatistics
          :noindex:
    
        """
        if (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ConsensusMap)) and (isinstance(args[2], ItraqEightPlexQuantitationMethod)):
            return self._correctIsotopicImpurities_0(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ConsensusMap)) and (isinstance(args[2], ItraqFourPlexQuantitationMethod)):
            return self._correctIsotopicImpurities_1(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ConsensusMap)) and (isinstance(args[2], TMTSixPlexQuantitationMethod)):
            return self._correctIsotopicImpurities_2(*args)
        elif (len(args)==3) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ConsensusMap)) and (isinstance(args[2], TMTTenPlexQuantitationMethod)):
            return self._correctIsotopicImpurities_3(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class IsobaricQuantifierStatistics:
    """
    Cython implementation of _IsobaricQuantifierStatistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IsobaricQuantifierStatistics.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property channel_count:
        def __set__(self,  channel_count):
        
            self.inst.get().channel_count = (<size_t>channel_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().channel_count
            py_result = <size_t>_r
            return py_result
    
    property iso_number_ms2_negative:
        def __set__(self,  iso_number_ms2_negative):
        
            self.inst.get().iso_number_ms2_negative = (<size_t>iso_number_ms2_negative)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().iso_number_ms2_negative
            py_result = <size_t>_r
            return py_result
    
    property iso_number_reporter_negative:
        def __set__(self,  iso_number_reporter_negative):
        
            self.inst.get().iso_number_reporter_negative = (<size_t>iso_number_reporter_negative)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().iso_number_reporter_negative
            py_result = <size_t>_r
            return py_result
    
    property iso_number_reporter_different:
        def __set__(self,  iso_number_reporter_different):
        
            self.inst.get().iso_number_reporter_different = (<size_t>iso_number_reporter_different)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().iso_number_reporter_different
            py_result = <size_t>_r
            return py_result
    
    property iso_solution_different_intensity:
        def __set__(self, double iso_solution_different_intensity):
        
            self.inst.get().iso_solution_different_intensity = (<double>iso_solution_different_intensity)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().iso_solution_different_intensity
            py_result = <double>_r
            return py_result
    
    property iso_total_intensity_negative:
        def __set__(self, double iso_total_intensity_negative):
        
            self.inst.get().iso_total_intensity_negative = (<double>iso_total_intensity_negative)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().iso_total_intensity_negative
            py_result = <double>_r
            return py_result
    
    property number_ms2_total:
        def __set__(self,  number_ms2_total):
        
            self.inst.get().number_ms2_total = (<size_t>number_ms2_total)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().number_ms2_total
            py_result = <size_t>_r
            return py_result
    
    property number_ms2_empty:
        def __set__(self,  number_ms2_empty):
        
            self.inst.get().number_ms2_empty = (<size_t>number_ms2_empty)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().number_ms2_empty
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef IsobaricQuantifierStatistics rv = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
       rv.inst = shared_ptr[_IsobaricQuantifierStatistics](new _IsobaricQuantifierStatistics(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IsobaricQuantifierStatistics rv = IsobaricQuantifierStatistics.__new__(IsobaricQuantifierStatistics)
       rv.inst = shared_ptr[_IsobaricQuantifierStatistics](new _IsobaricQuantifierStatistics(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IsobaricQuantifierStatistics](new _IsobaricQuantifierStatistics())
    
    def _init_1(self, IsobaricQuantifierStatistics in_0 ):
        """
        _init_1(self, in_0: IsobaricQuantifierStatistics ) -> None
        """
        assert isinstance(in_0, IsobaricQuantifierStatistics), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IsobaricQuantifierStatistics](new _IsobaricQuantifierStatistics((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IsobaricQuantifierStatistics ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IsobaricQuantifierStatistics)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def reset(self):
        """
        reset(self) -> None
        """
        self.inst.get().reset() 

cdef class MetaboliteFeatureDeconvolution:
    """
    Cython implementation of _MetaboliteFeatureDeconvolution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MetaboliteFeatureDeconvolution.html>`_
      -- Inherits from ['DefaultParamHandler']

    An algorithm to decharge small molecule features (i.e. as found by FeatureFinder)
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MetaboliteFeatureDeconvolution rv = MetaboliteFeatureDeconvolution.__new__(MetaboliteFeatureDeconvolution)
       rv.inst = shared_ptr[_MetaboliteFeatureDeconvolution](new _MetaboliteFeatureDeconvolution(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MetaboliteFeatureDeconvolution rv = MetaboliteFeatureDeconvolution.__new__(MetaboliteFeatureDeconvolution)
       rv.inst = shared_ptr[_MetaboliteFeatureDeconvolution](new _MetaboliteFeatureDeconvolution(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MetaboliteFeatureDeconvolution](new _MetaboliteFeatureDeconvolution())
    
    def _init_1(self, MetaboliteFeatureDeconvolution in_0 ):
        """
        _init_1(self, in_0: MetaboliteFeatureDeconvolution ) -> None
        """
        assert isinstance(in_0, MetaboliteFeatureDeconvolution), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MetaboliteFeatureDeconvolution](new _MetaboliteFeatureDeconvolution((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MetaboliteFeatureDeconvolution ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MetaboliteFeatureDeconvolution)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def compute(self, FeatureMap fm_in , FeatureMap fm_out , ConsensusMap cons_map , ConsensusMap cons_map_p ):
        """
        compute(self, fm_in: FeatureMap , fm_out: FeatureMap , cons_map: ConsensusMap , cons_map_p: ConsensusMap ) -> None
        Compute a zero-charge feature map from a set of charged features
        
        Find putative ChargePairs, then score them and hand over to ILP
        
        
        :param fm_in: Input feature-map
        :param fm_out: Output feature-map (sorted by position and augmented with user params)
        :param cons_map: Output of grouped features belonging to a charge group
        :param cons_map_p: Output of paired features connected by an edge
        """
        assert isinstance(fm_in, FeatureMap), 'arg fm_in wrong type'
        assert isinstance(fm_out, FeatureMap), 'arg fm_out wrong type'
        assert isinstance(cons_map, ConsensusMap), 'arg cons_map wrong type'
        assert isinstance(cons_map_p, ConsensusMap), 'arg cons_map_p wrong type'
    
    
    
    
        self.inst.get().compute((deref(fm_in.inst.get())), (deref(fm_out.inst.get())), (deref(cons_map.inst.get())), (deref(cons_map_p.inst.get())))
    
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
    CHARGEMODE_MFD = __CHARGEMODE_MFD 

cdef class MzIdentMLFile:
    """
    Cython implementation of _MzIdentMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzIdentMLFile.html>`_
      -- Inherits from ['ProgressLogger']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MzIdentMLFile rv = MzIdentMLFile.__new__(MzIdentMLFile)
       rv.inst = shared_ptr[_MzIdentMLFile](new _MzIdentMLFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MzIdentMLFile rv = MzIdentMLFile.__new__(MzIdentMLFile)
       rv.inst = shared_ptr[_MzIdentMLFile](new _MzIdentMLFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MzIdentMLFile](new _MzIdentMLFile())
    
    def _init_1(self, MzIdentMLFile in_0 ):
        """
        _init_1(self, in_0: MzIdentMLFile ) -> None
        """
        assert isinstance(in_0, MzIdentMLFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MzIdentMLFile](new _MzIdentMLFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MzIdentMLFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MzIdentMLFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , list poid , PeptideIdentificationList peid ):
        """
        load(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None
        Loads the identifications from a MzIdentML file
        
        
        :param filename: File name of the file to be checked
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsin
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(poid, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in poid), 'arg poid wrong type'
        assert isinstance(peid, PeptideIdentificationList), 'arg peid wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in poid:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peid.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_poid = v1.begin()
        replace_0 = []
        while it_poid != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_poid)))
            replace_0.append(item1)
            inc(it_poid)
        poid[:] = replace_0
        del v1
    
    def store(self,  filename , list poid , PeptideIdentificationList peid ):
        """
        store(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None
        Stores the identifications in a MzIdentML file
        
        
        :raises:
          Exception: UnableToCreateFile is thrown if the file could not be created
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
    
    def isSemanticallyValid(self,  filename , list errors , list warnings ):
        """
        isSemanticallyValid(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool
        Checks if a file is valid with respect to the mapping file and the controlled vocabulary
        
        
        :param filename: File name of the file to be checked
        :param errors: Errors during the validation are returned in this output parameter
        :param warnings: Warnings during the validation are returned in this output parameter
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(errors, list) and all(isinstance(li, bytes) for li in errors), 'arg errors wrong type'
        assert isinstance(warnings, list) and all(isinstance(li, bytes) for li in warnings), 'arg warnings wrong type'
    
        cdef libcpp_vector[_String] * v1 = new libcpp_vector[_String]()
        cdef bytes item1
        for item1 in errors:
           v1.push_back(_String(<char *>item1))
        cdef libcpp_vector[_String] * v2 = new libcpp_vector[_String]()
        cdef bytes item2
        for item2 in warnings:
           v2.push_back(_String(<char *>item2))
        cdef bool _r = self.inst.get().isSemanticallyValid(deref((convString(filename)).get()), deref(v1), deref(v2))
        del v2
        del v1
        py_result = <bool>_r
        return py_result
    
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

cdef class Normalizer:
    """
    Cython implementation of _Normalizer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Normalizer.html>`_
      -- Inherits from ['DefaultParamHandler']

    Normalizes the peak intensities spectrum-wise
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_Normalizer](new _Normalizer())
    
    def _init_1(self, Normalizer in_0 ):
        """
        _init_1(self, in_0: Normalizer ) -> None
        """
        assert isinstance(in_0, Normalizer), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_Normalizer](new _Normalizer((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: Normalizer ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], Normalizer)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def filterSpectrum(self, MSSpectrum spec ):
        """
        filterSpectrum(self, spec: MSSpectrum ) -> None
        Normalizes the spectrum
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterSpectrum((deref(spec.inst.get())))
    
    def filterPeakSpectrum(self, MSSpectrum spec ):
        """
        filterPeakSpectrum(self, spec: MSSpectrum ) -> None
        Normalizes the peak spectrum
        """
        assert isinstance(spec, MSSpectrum), 'arg spec wrong type'
    
        self.inst.get().filterPeakSpectrum((deref(spec.inst.get())))
    
    def filterPeakMap(self, MSExperiment exp ):
        """
        filterPeakMap(self, exp: MSExperiment ) -> None
        Normalizes the peak map
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

cdef class OPXLDataStructs:
    """
    Cython implementation of _OPXLDataStructs

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OPXLDataStructs.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef OPXLDataStructs rv = OPXLDataStructs.__new__(OPXLDataStructs)
       rv.inst = shared_ptr[_OPXLDataStructs](new _OPXLDataStructs(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef OPXLDataStructs rv = OPXLDataStructs.__new__(OPXLDataStructs)
       rv.inst = shared_ptr[_OPXLDataStructs](new _OPXLDataStructs(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_OPXLDataStructs](new _OPXLDataStructs())
    
    def _init_1(self, OPXLDataStructs in_0 ):
        """
        _init_1(self, in_0: OPXLDataStructs ) -> None
        """
        assert isinstance(in_0, OPXLDataStructs), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_OPXLDataStructs](new _OPXLDataStructs((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: OPXLDataStructs ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], OPXLDataStructs)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    PeptidePosition = __PeptidePosition
    ProteinProteinCrossLinkType = __ProteinProteinCrossLinkType 

cdef class PeptideAndProteinQuant:
    """
    Cython implementation of _PeptideAndProteinQuant

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef PeptideAndProteinQuant rv = PeptideAndProteinQuant.__new__(PeptideAndProteinQuant)
       rv.inst = shared_ptr[_PeptideAndProteinQuant](new _PeptideAndProteinQuant(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideAndProteinQuant rv = PeptideAndProteinQuant.__new__(PeptideAndProteinQuant)
       rv.inst = shared_ptr[_PeptideAndProteinQuant](new _PeptideAndProteinQuant(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Helper class for peptide and protein quantification based on feature data annotated with IDs
        """
        self.inst = shared_ptr[_PeptideAndProteinQuant](new _PeptideAndProteinQuant())
    
    def _init_1(self, PeptideAndProteinQuant in_0 ):
        """
        _init_1(self, in_0: PeptideAndProteinQuant ) -> None
        """
        assert isinstance(in_0, PeptideAndProteinQuant), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideAndProteinQuant](new _PeptideAndProteinQuant((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Helper class for peptide and protein quantification based on feature data annotated with IDs

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideAndProteinQuant ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideAndProteinQuant)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _readQuantData_0(self, FeatureMap map_in , ExperimentalDesign ed ):
        """
        _readQuantData_0(self, map_in: FeatureMap , ed: ExperimentalDesign ) -> None
        Read quantitative data from a feature map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        assert isinstance(map_in, FeatureMap), 'arg map_in wrong type'
        assert isinstance(ed, ExperimentalDesign), 'arg ed wrong type'
    
    
        self.inst.get().readQuantData((deref(map_in.inst.get())), (deref(ed.inst.get())))
    
    def _readQuantData_1(self, ConsensusMap map_in , ExperimentalDesign ed ):
        """
        _readQuantData_1(self, map_in: ConsensusMap , ed: ExperimentalDesign ) -> None
        Read quantitative data from a consensus map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        assert isinstance(map_in, ConsensusMap), 'arg map_in wrong type'
        assert isinstance(ed, ExperimentalDesign), 'arg ed wrong type'
    
    
        self.inst.get().readQuantData((deref(map_in.inst.get())), (deref(ed.inst.get())))
    
    def _readQuantData_2(self, list proteins , PeptideIdentificationList peptides , ExperimentalDesign ed ):
        """
        _readQuantData_2(self, proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList , ed: ExperimentalDesign ) -> None
        Read quantitative data from identification results (for quantification via spectral counting)
        
        Parameters should be set before using this method, as setting parameters will clear all results
        """
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(ed, ExperimentalDesign), 'arg ed wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in proteins:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().readQuantData(deref(v0), (deref(peptides.inst.get())), (deref(ed.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v0.begin()
        replace_0 = []
        while it_proteins != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item0)
            inc(it_proteins)
        proteins[:] = replace_0
        del v0
    
    def readQuantData(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: readQuantData(self, map_in: FeatureMap , ed: ExperimentalDesign ) -> None
          :noindex:
        
        Read quantitative data from a feature map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        
        .. rubric:: Overload:
        .. py:function:: readQuantData(self, map_in: ConsensusMap , ed: ExperimentalDesign ) -> None
          :noindex:
        
        Read quantitative data from a consensus map
        
        Parameters should be set before using this method, as setting parameters will clear all results
        
        .. rubric:: Overload:
        .. py:function:: readQuantData(self, proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList , ed: ExperimentalDesign ) -> None
          :noindex:
        
        Read quantitative data from identification results (for quantification via spectral counting)
        
        Parameters should be set before using this method, as setting parameters will clear all results
    
        """
        if (len(args)==2) and (isinstance(args[0], FeatureMap)) and (isinstance(args[1], ExperimentalDesign)):
            return self._readQuantData_0(*args)
        elif (len(args)==2) and (isinstance(args[0], ConsensusMap)) and (isinstance(args[1], ExperimentalDesign)):
            return self._readQuantData_1(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], PeptideIdentificationList)) and (isinstance(args[2], ExperimentalDesign)):
            return self._readQuantData_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def quantifyPeptides(self, PeptideIdentificationList peptides ):
        """
        quantifyPeptides(self, peptides: PeptideIdentificationList ) -> None
        Compute peptide abundances
        
        Based on quantitative data for individual charge states (in member `pep_quant_`), overall abundances for peptides are computed (and stored again in `pep_quant_`)
        Quantitative data must first be read via readQuantData()
        Optional (peptide-level) protein inference information (e.g. from Fido or ProteinProphet) can be supplied via `peptides`. In that case, peptide-to-protein associations - the basis for protein-level quantification - will also be read from `peptides`!
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().quantifyPeptides((deref(peptides.inst.get())))
    
    def quantifyProteins(self, ProteinIdentification proteins ):
        """
        quantifyProteins(self, proteins: ProteinIdentification ) -> None
        Compute protein abundances
        
        Peptide abundances must be computed first with quantifyPeptides(). Optional protein inference information (e.g. from Fido or ProteinProphet) can be supplied via `proteins`
        """
        assert isinstance(proteins, ProteinIdentification), 'arg proteins wrong type'
    
        self.inst.get().quantifyProteins((deref(proteins.inst.get())))
    
    def getStatistics(self):
        """
        getStatistics(self) -> PeptideAndProteinQuant_Statistics
        """
        cdef _PeptideAndProteinQuant_Statistics * _r = new _PeptideAndProteinQuant_Statistics(self.inst.get().getStatistics())
        cdef PeptideAndProteinQuant_Statistics py_result = PeptideAndProteinQuant_Statistics.__new__(PeptideAndProteinQuant_Statistics)
        py_result.inst = shared_ptr[_PeptideAndProteinQuant_Statistics](_r)
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

cdef class PeptideAndProteinQuant_PeptideData:
    """
    Cython implementation of _PeptideAndProteinQuant_PeptideData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_PeptideData.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property accessions:
        def __set__(self, set accessions):
            cdef libcpp_set[_String] * v0 = new libcpp_set[_String]()
            cdef bytes item0
            for item0 in accessions:
               v0.insert(_String(<char *>item0))
            self.inst.get().accessions = deref(v0)
            del v0
    
        def __get__(self):
            _r = self.inst.get().accessions
            py_result = set()
            cdef libcpp_set[_String].iterator it__r = _r.begin()
            while it__r != _r.end():
               py_result.add(<char*>deref(it__r).c_str())
               inc(it__r)
            return py_result
    
    property psm_count:
        def __set__(self,  psm_count):
        
            self.inst.get().psm_count = (<size_t>psm_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().psm_count
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef PeptideAndProteinQuant_PeptideData rv = PeptideAndProteinQuant_PeptideData.__new__(PeptideAndProteinQuant_PeptideData)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_PeptideData](new _PeptideAndProteinQuant_PeptideData(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideAndProteinQuant_PeptideData rv = PeptideAndProteinQuant_PeptideData.__new__(PeptideAndProteinQuant_PeptideData)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_PeptideData](new _PeptideAndProteinQuant_PeptideData(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideAndProteinQuant_PeptideData](new _PeptideAndProteinQuant_PeptideData())
    
    def _init_1(self, PeptideAndProteinQuant_PeptideData in_0 ):
        """
        _init_1(self, in_0: PeptideAndProteinQuant_PeptideData ) -> None
        """
        assert isinstance(in_0, PeptideAndProteinQuant_PeptideData), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideAndProteinQuant_PeptideData](new _PeptideAndProteinQuant_PeptideData((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideAndProteinQuant_PeptideData ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideAndProteinQuant_PeptideData)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeptideAndProteinQuant_ProteinData:
    """
    Cython implementation of _PeptideAndProteinQuant_ProteinData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_ProteinData.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property psm_count:
        def __set__(self,  psm_count):
        
            self.inst.get().psm_count = (<size_t>psm_count)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().psm_count
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef PeptideAndProteinQuant_ProteinData rv = PeptideAndProteinQuant_ProteinData.__new__(PeptideAndProteinQuant_ProteinData)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_ProteinData](new _PeptideAndProteinQuant_ProteinData(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideAndProteinQuant_ProteinData rv = PeptideAndProteinQuant_ProteinData.__new__(PeptideAndProteinQuant_ProteinData)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_ProteinData](new _PeptideAndProteinQuant_ProteinData(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideAndProteinQuant_ProteinData](new _PeptideAndProteinQuant_ProteinData())
    
    def _init_1(self, PeptideAndProteinQuant_ProteinData in_0 ):
        """
        _init_1(self, in_0: PeptideAndProteinQuant_ProteinData ) -> None
        """
        assert isinstance(in_0, PeptideAndProteinQuant_ProteinData), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideAndProteinQuant_ProteinData](new _PeptideAndProteinQuant_ProteinData((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideAndProteinQuant_ProteinData ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideAndProteinQuant_ProteinData)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class PeptideAndProteinQuant_Statistics:
    """
    Cython implementation of _PeptideAndProteinQuant_Statistics

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideAndProteinQuant_Statistics.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property n_samples:
        def __set__(self,  n_samples):
        
            self.inst.get().n_samples = (<size_t>n_samples)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().n_samples
            py_result = <size_t>_r
            return py_result
    
    property quant_proteins:
        def __set__(self,  quant_proteins):
        
            self.inst.get().quant_proteins = (<size_t>quant_proteins)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().quant_proteins
            py_result = <size_t>_r
            return py_result
    
    property too_few_peptides:
        def __set__(self,  too_few_peptides):
        
            self.inst.get().too_few_peptides = (<size_t>too_few_peptides)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().too_few_peptides
            py_result = <size_t>_r
            return py_result
    
    property quant_peptides:
        def __set__(self,  quant_peptides):
        
            self.inst.get().quant_peptides = (<size_t>quant_peptides)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().quant_peptides
            py_result = <size_t>_r
            return py_result
    
    property total_peptides:
        def __set__(self,  total_peptides):
        
            self.inst.get().total_peptides = (<size_t>total_peptides)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().total_peptides
            py_result = <size_t>_r
            return py_result
    
    property quant_features:
        def __set__(self,  quant_features):
        
            self.inst.get().quant_features = (<size_t>quant_features)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().quant_features
            py_result = <size_t>_r
            return py_result
    
    property total_features:
        def __set__(self,  total_features):
        
            self.inst.get().total_features = (<size_t>total_features)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().total_features
            py_result = <size_t>_r
            return py_result
    
    property blank_features:
        def __set__(self,  blank_features):
        
            self.inst.get().blank_features = (<size_t>blank_features)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().blank_features
            py_result = <size_t>_r
            return py_result
    
    property ambig_features:
        def __set__(self,  ambig_features):
        
            self.inst.get().ambig_features = (<size_t>ambig_features)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().ambig_features
            py_result = <size_t>_r
            return py_result
    
    def __copy__(self):
       cdef PeptideAndProteinQuant_Statistics rv = PeptideAndProteinQuant_Statistics.__new__(PeptideAndProteinQuant_Statistics)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_Statistics](new _PeptideAndProteinQuant_Statistics(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef PeptideAndProteinQuant_Statistics rv = PeptideAndProteinQuant_Statistics.__new__(PeptideAndProteinQuant_Statistics)
       rv.inst = shared_ptr[_PeptideAndProteinQuant_Statistics](new _PeptideAndProteinQuant_Statistics(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_PeptideAndProteinQuant_Statistics](new _PeptideAndProteinQuant_Statistics())
    
    def _init_1(self, PeptideAndProteinQuant_Statistics in_0 ):
        """
        _init_1(self, in_0: PeptideAndProteinQuant_Statistics ) -> None
        """
        assert isinstance(in_0, PeptideAndProteinQuant_Statistics), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_PeptideAndProteinQuant_Statistics](new _PeptideAndProteinQuant_Statistics((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: PeptideAndProteinQuant_Statistics ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], PeptideAndProteinQuant_Statistics)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class ProgressLogger:
    """
    Cython implementation of _ProgressLogger

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProgressLogger.html>`_

    Base class for all classes that want to report their progress
    
    Per default the progress log is disabled. Use setLogType to enable it
    
    Use startProgress, setProgress and endProgress for the actual logging
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ProgressLogger rv = ProgressLogger.__new__(ProgressLogger)
       rv.inst = shared_ptr[_ProgressLogger](new _ProgressLogger(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ProgressLogger rv = ProgressLogger.__new__(ProgressLogger)
       rv.inst = shared_ptr[_ProgressLogger](new _ProgressLogger(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ProgressLogger](new _ProgressLogger())
    
    def _init_1(self, ProgressLogger in_0 ):
        """
        _init_1(self, in_0: ProgressLogger ) -> None
        """
        assert isinstance(in_0, ProgressLogger), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ProgressLogger](new _ProgressLogger((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ProgressLogger ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ProgressLogger)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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

cdef class SpectrumMetaData:
    """
    Cython implementation of _SpectrumMetaData

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumMetaData.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property rt:
        def __set__(self, double rt):
        
            self.inst.get().rt = (<double>rt)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().rt
            py_result = <double>_r
            return py_result
    
    property precursor_rt:
        def __set__(self, double precursor_rt):
        
            self.inst.get().precursor_rt = (<double>precursor_rt)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_rt
            py_result = <double>_r
            return py_result
    
    property precursor_mz:
        def __set__(self, double precursor_mz):
        
            self.inst.get().precursor_mz = (<double>precursor_mz)
        
    
        def __get__(self):
            cdef double _r = self.inst.get().precursor_mz
            py_result = <double>_r
            return py_result
    
    property precursor_charge:
        def __set__(self,  precursor_charge):
        
            self.inst.get().precursor_charge = (<int>precursor_charge)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().precursor_charge
            py_result = <int>_r
            return py_result
    
    property ms_level:
        def __set__(self,  ms_level):
        
            self.inst.get().ms_level = (<size_t>ms_level)
        
    
        def __get__(self):
            cdef size_t _r = self.inst.get().ms_level
            py_result = <size_t>_r
            return py_result
    
    property scan_number:
        def __set__(self,  scan_number):
        
            self.inst.get().scan_number = (<int>scan_number)
        
    
        def __get__(self):
            cdef int _r = self.inst.get().scan_number
            py_result = <int>_r
            return py_result
    
    property native_id:
        def __set__(self,  native_id):
        
            self.inst.get().native_id = deref((convString(native_id)).get())
        
    
        def __get__(self):
            cdef _String _r = self.inst.get().native_id
            py_result = convOutputString(_r)
            return py_result
    
    def __copy__(self):
       cdef SpectrumMetaData rv = SpectrumMetaData.__new__(SpectrumMetaData)
       rv.inst = shared_ptr[_SpectrumMetaData](new _SpectrumMetaData(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SpectrumMetaData rv = SpectrumMetaData.__new__(SpectrumMetaData)
       rv.inst = shared_ptr[_SpectrumMetaData](new _SpectrumMetaData(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SpectrumMetaData](new _SpectrumMetaData())
    
    def _init_1(self, SpectrumMetaData in_0 ):
        """
        _init_1(self, in_0: SpectrumMetaData ) -> None
        """
        assert isinstance(in_0, SpectrumMetaData), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SpectrumMetaData](new _SpectrumMetaData((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SpectrumMetaData ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SpectrumMetaData)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,)) 

cdef class SpectrumMetaDataLookup:
    """
    Cython implementation of _SpectrumMetaDataLookup

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumMetaDataLookup.html>`_
      -- Inherits from ['SpectrumLookup']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_SpectrumMetaDataLookup](new _SpectrumMetaDataLookup())
    
    def _readSpectra_0(self, MSExperiment spectra ,  scan_regexp , bool get_precursor_rt ):
        """
        _readSpectra_0(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] , get_precursor_rt: bool ) -> None
        Read spectra and store their meta data
        
        :param SpectrumContainer: Spectrum container class, must support `size` and `operator[]`
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")
        :param get_precursor_rt: Assign precursor retention times? (This relies on all precursor spectra being present and in the right order.)
        """
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert (isinstance(scan_regexp, str) or isinstance(scan_regexp, bytes) or isinstance(scan_regexp, String)), 'arg scan_regexp wrong type'
        assert isinstance(get_precursor_rt, pybool_t), 'arg get_precursor_rt wrong type'
    
    
    
        self.inst.get().readSpectra((deref(spectra.inst.get())), deref((convString(scan_regexp)).get()), (<bool>get_precursor_rt))
    
    def _readSpectra_1(self, MSExperiment spectra ,  scan_regexp ):
        """
        _readSpectra_1(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] ) -> None
        Read and index spectra for later look-up
        
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
        """
        assert isinstance(spectra, MSExperiment), 'arg spectra wrong type'
        assert (isinstance(scan_regexp, str) or isinstance(scan_regexp, bytes) or isinstance(scan_regexp, String)), 'arg scan_regexp wrong type'
    
    
        self.inst.get().readSpectra((deref(spectra.inst.get())), deref((convString(scan_regexp)).get()))
    
    def readSpectra(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] , get_precursor_rt: bool ) -> None
          :noindex:
        
        Read spectra and store their meta data
        
        :param SpectrumContainer: Spectrum container class, must support `size` and `operator[]`
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")
        :param get_precursor_rt: Assign precursor retention times? (This relies on all precursor spectra being present and in the right order.)
        
        .. rubric:: Overload:
        .. py:function:: readSpectra(self, spectra: MSExperiment , scan_regexp: Union[bytes, str, String] ) -> None
          :noindex:
        
        Read and index spectra for later look-up
        
        :param spectra: Container of spectra
        :param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
    
        """
        if (len(args)==3) and (isinstance(args[0], MSExperiment)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))) and (isinstance(args[2], pybool_t)):
            return self._readSpectra_0(*args)
        elif (len(args)==2) and (isinstance(args[0], MSExperiment)) and ((isinstance(args[1], str) or isinstance(args[1], bytes) or isinstance(args[1], String))):
            return self._readSpectra_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getSpectrumMetaData_0(self,  index , SpectrumMetaData meta ):
        """
        _getSpectrumMetaData_0(self, index: int , meta: SpectrumMetaData ) -> None
        Look up meta data of a spectrum
        
        :param index: Index of the spectrum
        :param meta: Meta data output
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(meta, SpectrumMetaData), 'arg meta wrong type'
    
    
        self.inst.get().getSpectrumMetaData((<size_t>index), (deref(meta.inst.get())))
    
    def _getSpectrumMetaData_1(self,  spectrum_ref , SpectrumMetaData meta ):
        """
        _getSpectrumMetaData_1(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData ) -> None
        Extract meta data from a spectrum
        
        :param spectrum: Spectrum input
        :param meta: Meta data output
        :param scan_regexp: Regular expression for extracting scan number from spectrum native ID
        :param precursor_rts: RTs of potential precursor spectra of different MS levels
        """
        assert (isinstance(spectrum_ref, str) or isinstance(spectrum_ref, bytes) or isinstance(spectrum_ref, String)), 'arg spectrum_ref wrong type'
        assert isinstance(meta, SpectrumMetaData), 'arg meta wrong type'
    
    
        self.inst.get().getSpectrumMetaData(deref((convString(spectrum_ref)).get()), (deref(meta.inst.get())))
    
    def _getSpectrumMetaData_2(self,  spectrum_ref , SpectrumMetaData meta , bytes flags ):
        """
        _getSpectrumMetaData_2(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData , flags: bytes ) -> None
        Extract meta data via a spectrum reference
        
        :param spectrum_ref: Spectrum reference to parse
        :param metadata: Meta data output
        :param flags: What meta data to extract
        """
        assert (isinstance(spectrum_ref, str) or isinstance(spectrum_ref, bytes) or isinstance(spectrum_ref, String)), 'arg spectrum_ref wrong type'
        assert isinstance(meta, SpectrumMetaData), 'arg meta wrong type'
        assert isinstance(flags, bytes) and len(flags) == 1, 'arg flags wrong type'
    
    
    
        self.inst.get().getSpectrumMetaData(deref((convString(spectrum_ref)).get()), (deref(meta.inst.get())), (<char>((flags)[0])))
    
    def getSpectrumMetaData(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getSpectrumMetaData(self, index: int , meta: SpectrumMetaData ) -> None
          :noindex:
        
        Look up meta data of a spectrum
        
        :param index: Index of the spectrum
        :param meta: Meta data output
        
        .. rubric:: Overload:
        .. py:function:: getSpectrumMetaData(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData ) -> None
          :noindex:
        
        Extract meta data from a spectrum
        
        :param spectrum: Spectrum input
        :param meta: Meta data output
        :param scan_regexp: Regular expression for extracting scan number from spectrum native ID
        :param precursor_rts: RTs of potential precursor spectra of different MS levels
        
        .. rubric:: Overload:
        .. py:function:: getSpectrumMetaData(self, spectrum_ref: Union[bytes, str, String] , meta: SpectrumMetaData , flags: bytes ) -> None
          :noindex:
        
        Extract meta data via a spectrum reference
        
        :param spectrum_ref: Spectrum reference to parse
        :param metadata: Meta data output
        :param flags: What meta data to extract
    
        """
        if (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], SpectrumMetaData)):
            return self._getSpectrumMetaData_0(*args)
        elif (len(args)==2) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], SpectrumMetaData)):
            return self._getSpectrumMetaData_1(*args)
        elif (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], SpectrumMetaData)) and (isinstance(args[2], bytes) and len(args[2]) == 1):
            return self._getSpectrumMetaData_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setSpectraDataRef(self,  spectra_data ):
        """
        setSpectraDataRef(self, spectra_data: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(spectra_data, str) or isinstance(spectra_data, bytes) or isinstance(spectra_data, String)), 'arg spectra_data wrong type'
    
        self.inst.get().setSpectraDataRef(deref((convString(spectra_data)).get()))
    
    def empty(self):
        """
        empty(self) -> bool
        Check if any spectra were set
        """
        cdef bool _r = self.inst.get().empty()
        py_result = <bool>_r
        return py_result
    
    def findByRT(self, double rt ):
        """
        findByRT(self, rt: float ) -> int
        Look up spectrum by retention time (RT)
        
        :param rt: Retention time to look up
        :returns: Index of the spectrum that matched
        """
        assert isinstance(rt, float), 'arg rt wrong type'
    
        cdef size_t _r = self.inst.get().findByRT((<double>rt))
        py_result = <size_t>_r
        return py_result
    
    def findByNativeID(self,  native_id ):
        """
        findByNativeID(self, native_id: Union[bytes, str, String] ) -> int
        Look up spectrum by native ID
        
        :param native_id: Native ID to look up
        :returns: Index of the spectrum that matched
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
    
        cdef size_t _r = self.inst.get().findByNativeID(deref((convString(native_id)).get()))
        py_result = <size_t>_r
        return py_result
    
    def findByIndex(self,  index , bool count_from_one ):
        """
        findByIndex(self, index: int , count_from_one: bool ) -> int
        Look up spectrum by index (position in the vector of spectra)
        
        :param index: Index to look up
        :param count_from_one: Do indexes start counting at one (default zero)?
        :returns: Index of the spectrum that matched
        """
        assert isinstance(index, int) and index >= 0, 'arg index wrong type'
        assert isinstance(count_from_one, pybool_t), 'arg count_from_one wrong type'
    
    
        cdef size_t _r = self.inst.get().findByIndex((<size_t>index), (<bool>count_from_one))
        py_result = <size_t>_r
        return py_result
    
    def findByScanNumber(self,  scan_number ):
        """
        findByScanNumber(self, scan_number: int ) -> int
        Look up spectrum by scan number (extracted from the native ID)
        
        :param scan_number: Scan number to look up
        :returns: Index of the spectrum that matched
        """
        assert isinstance(scan_number, int) and scan_number >= 0, 'arg scan_number wrong type'
    
        cdef size_t _r = self.inst.get().findByScanNumber((<size_t>scan_number))
        py_result = <size_t>_r
        return py_result
    
    def findByReference(self,  spectrum_ref ):
        """
        findByReference(self, spectrum_ref: Union[bytes, str, String] ) -> int
        Look up spectrum by reference
        
        :param spectrum_ref: Spectrum reference to parse
        :returns: Index of the spectrum that matched
        """
        assert (isinstance(spectrum_ref, str) or isinstance(spectrum_ref, bytes) or isinstance(spectrum_ref, String)), 'arg spectrum_ref wrong type'
    
        cdef size_t _r = self.inst.get().findByReference(deref((convString(spectrum_ref)).get()))
        py_result = <size_t>_r
        return py_result
    
    def addReferenceFormat(self,  regexp ):
        """
        addReferenceFormat(self, regexp: Union[bytes, str, String] ) -> None
        Register a possible format for a spectrum reference
        
        :param regexp: Regular expression defining the format
        """
        assert (isinstance(regexp, str) or isinstance(regexp, bytes) or isinstance(regexp, String)), 'arg regexp wrong type'
    
        self.inst.get().addReferenceFormat(deref((convString(regexp)).get()))
    
    def extractScanNumber(self,  native_id ,  native_id_type_accession ):
        """
        extractScanNumber(self, native_id: Union[bytes, str, String] , native_id_type_accession: Union[bytes, str, String] ) -> int
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
        assert (isinstance(native_id_type_accession, str) or isinstance(native_id_type_accession, bytes) or isinstance(native_id_type_accession, String)), 'arg native_id_type_accession wrong type'
    
    
        cdef int _r = self.inst.get().extractScanNumber(deref((convString(native_id)).get()), deref((convString(native_id_type_accession)).get()))
        py_result = <int>_r
        return py_result
    addMissingIMToPeptideIDs = __static_SpectrumMetaDataLookup_addMissingIMToPeptideIDs
    addMissingRTsToPeptideIDs = __static_SpectrumMetaDataLookup_addMissingRTsToPeptideIDs
    addMissingSpectrumReferences = __static_SpectrumMetaDataLookup_addMissingSpectrumReferences
    getSpectrumMetaData = __static_SpectrumMetaDataLookup_getSpectrumMetaData 

cdef class StringDataArray:
    """
    Cython implementation of _StringDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1StringDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra string data attached to a spectrum or chromatogram.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef StringDataArray rv = StringDataArray.__new__(StringDataArray)
       rv.inst = shared_ptr[_StringDataArray](new _StringDataArray(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef StringDataArray rv = StringDataArray.__new__(StringDataArray)
       rv.inst = shared_ptr[_StringDataArray](new _StringDataArray(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_StringDataArray](new _StringDataArray())
    
    def _init_1(self, StringDataArray in_0 ):
        """
        _init_1(self, in_0: StringDataArray ) -> None
        """
        assert isinstance(in_0, StringDataArray), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_StringDataArray](new _StringDataArray((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: StringDataArray ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], StringDataArray)):
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
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def push_back(self,  in_0 ):
        """
        push_back(self, in_0: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(in_0, str) or isinstance(in_0, bytes) or isinstance(in_0, String)), 'arg in_0 wrong type'
    
        self.inst.get().push_back(deref((convString(in_0)).get()))
    
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
        if not isinstance(other, StringDataArray):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef StringDataArray other_casted = other
        cdef StringDataArray self_casted = self
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

        cdef _String _r = deref(self.inst.get())[(<int>in_0)]
        py_result = _cast_const_away(<char*>_r.c_str())
        return py_result

    def __setitem__(self, key, value):
        assert isinstance(key, int), 'arg key wrong type'
        assert (isinstance(value, str) or isinstance(value, unicode) or isinstance(value, bytes) or isinstance(value, String)), 'arg value wrong type'
        assert key >= 0, 'arg key cannot be negative'

        cdef unsigned int _idx = (<int>key)
        if _idx >= self.inst.get().size():
            raise IndexError("invalid index %d" % _idx)

        cdef shared_ptr[_String] _s = convString(value)
        deref(self.inst.get())[(<int>key)] = deref(_s.get()) 
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
