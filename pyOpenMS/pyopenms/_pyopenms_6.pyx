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

cdef class __ChromatogramType:
    None
    MASS_CHROMATOGRAM = 0
    TOTAL_ION_CURRENT_CHROMATOGRAM = 1
    SELECTED_ION_CURRENT_CHROMATOGRAM = 2
    BASEPEAK_CHROMATOGRAM = 3
    SELECTED_ION_MONITORING_CHROMATOGRAM = 4
    SELECTED_REACTION_MONITORING_CHROMATOGRAM = 5
    ELECTROMAGNETIC_RADIATION_CHROMATOGRAM = 6
    ABSORPTION_CHROMATOGRAM = 7
    EMISSION_CHROMATOGRAM = 8
    SIZE_OF_CHROMATOGRAM_TYPE = 9

    def getMapping(self):
        return dict([ (v, k) for k, v in self.__class__.__dict__.items() if isinstance(v, int) ]) 

cdef class AMSE_AdductInfo:
    """
    Cython implementation of _AMSE_AdductInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AMSE_AdductInfo.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self,  name , EmpiricalFormula adduct ,  charge ,  mol_multiplier ):
        """
        __init__(self, name: Union[bytes, str, String] , adduct: EmpiricalFormula , charge: int , mol_multiplier: int ) -> None
        """
        assert (isinstance(name, str) or isinstance(name, bytes) or isinstance(name, String)), 'arg name wrong type'
        assert isinstance(adduct, EmpiricalFormula), 'arg adduct wrong type'
        assert isinstance(charge, int), 'arg charge wrong type'
        assert isinstance(mol_multiplier, int), 'arg mol_multiplier wrong type'
    
    
    
    
        self.inst = shared_ptr[_AMSE_AdductInfo](new _AMSE_AdductInfo(deref((convString(name)).get()), (deref(adduct.inst.get())), (<int>charge), (<unsigned int>mol_multiplier)))
    
    def getNeutralMass(self, double observed_mz ):
        """
        getNeutralMass(self, observed_mz: float ) -> float
        """
        assert isinstance(observed_mz, float), 'arg observed_mz wrong type'
    
        cdef double _r = self.inst.get().getNeutralMass((<double>observed_mz))
        py_result = <double>_r
        return py_result
    
    def getMZ(self, double neutral_mass ):
        """
        getMZ(self, neutral_mass: float ) -> float
        """
        assert isinstance(neutral_mass, float), 'arg neutral_mass wrong type'
    
        cdef double _r = self.inst.get().getMZ((<double>neutral_mass))
        py_result = <double>_r
        return py_result
    
    def isCompatible(self, EmpiricalFormula db_entry ):
        """
        isCompatible(self, db_entry: EmpiricalFormula ) -> bool
        """
        assert isinstance(db_entry, EmpiricalFormula), 'arg db_entry wrong type'
    
        cdef bool _r = self.inst.get().isCompatible((deref(db_entry.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def getCharge(self):
        """
        getCharge(self) -> int
        """
        cdef int _r = self.inst.get().getCharge()
        py_result = <int>_r
        return py_result
    
    def getName(self):
        """
        getName(self) -> Union[bytes, str, String]
        """
        cdef _String _r = self.inst.get().getName()
        py_result = convOutputString(_r)
        return py_result 

cdef class AbsoluteQuantitation:
    """
    Cython implementation of _AbsoluteQuantitation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitation.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef AbsoluteQuantitation rv = AbsoluteQuantitation.__new__(AbsoluteQuantitation)
       rv.inst = shared_ptr[_AbsoluteQuantitation](new _AbsoluteQuantitation(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef AbsoluteQuantitation rv = AbsoluteQuantitation.__new__(AbsoluteQuantitation)
       rv.inst = shared_ptr[_AbsoluteQuantitation](new _AbsoluteQuantitation(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_AbsoluteQuantitation](new _AbsoluteQuantitation())
    
    def _init_1(self, AbsoluteQuantitation in_0 ):
        """
        _init_1(self, in_0: AbsoluteQuantitation ) -> None
        """
        assert isinstance(in_0, AbsoluteQuantitation), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_AbsoluteQuantitation](new _AbsoluteQuantitation((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: AbsoluteQuantitation ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], AbsoluteQuantitation)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def setQuantMethods(self, list quant_methods ):
        """
        setQuantMethods(self, quant_methods: List[AbsoluteQuantitationMethod] ) -> None
        """
        assert isinstance(quant_methods, list) and all(isinstance(elemt_rec, AbsoluteQuantitationMethod) for elemt_rec in quant_methods), 'arg quant_methods wrong type'
        cdef libcpp_vector[_AbsoluteQuantitationMethod] * v0 = new libcpp_vector[_AbsoluteQuantitationMethod]()
        cdef AbsoluteQuantitationMethod item0
        for item0 in quant_methods:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setQuantMethods(deref(v0))
        cdef libcpp_vector[_AbsoluteQuantitationMethod].iterator it_quant_methods = v0.begin()
        replace_0 = []
        while it_quant_methods != v0.end():
            item0 = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
            item0.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(it_quant_methods)))
            replace_0.append(item0)
            inc(it_quant_methods)
        quant_methods[:] = replace_0
        del v0
    
    def getQuantMethods(self):
        """
        getQuantMethods(self) -> List[AbsoluteQuantitationMethod]
        """
        _r = self.inst.get().getQuantMethods()
        py_result = []
        cdef libcpp_vector[_AbsoluteQuantitationMethod].iterator it__r = _r.begin()
        cdef AbsoluteQuantitationMethod item_py_result
        while it__r != _r.end():
           item_py_result = AbsoluteQuantitationMethod.__new__(AbsoluteQuantitationMethod)
           item_py_result.inst = shared_ptr[_AbsoluteQuantitationMethod](new _AbsoluteQuantitationMethod(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def calculateRatio(self, Feature component_1 , Feature component_2 ,  feature_name ):
        """
        calculateRatio(self, component_1: Feature , component_2: Feature , feature_name: Union[bytes, str, String] ) -> float
        """
        assert isinstance(component_1, Feature), 'arg component_1 wrong type'
        assert isinstance(component_2, Feature), 'arg component_2 wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
    
    
    
        cdef double _r = self.inst.get().calculateRatio((deref(component_1.inst.get())), (deref(component_2.inst.get())), deref((convString(feature_name)).get()))
        py_result = <double>_r
        return py_result
    
    def applyCalibration(self, Feature component , Feature IS_component ,  feature_name ,  transformation_model , Param transformation_model_params ):
        """
        applyCalibration(self, component: Feature , IS_component: Feature , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param ) -> float
        """
        assert isinstance(component, Feature), 'arg component wrong type'
        assert isinstance(IS_component, Feature), 'arg IS_component wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert (isinstance(transformation_model, str) or isinstance(transformation_model, bytes) or isinstance(transformation_model, String)), 'arg transformation_model wrong type'
        assert isinstance(transformation_model_params, Param), 'arg transformation_model_params wrong type'
    
    
    
    
    
        cdef double _r = self.inst.get().applyCalibration((deref(component.inst.get())), (deref(IS_component.inst.get())), deref((convString(feature_name)).get()), deref((convString(transformation_model)).get()), (deref(transformation_model_params.inst.get())))
        py_result = <double>_r
        return py_result
    
    def quantifyComponents(self, FeatureMap unknowns ):
        """
        quantifyComponents(self, unknowns: FeatureMap ) -> None
        This function applies the calibration curve, hence quantifying all the components
        """
        assert isinstance(unknowns, FeatureMap), 'arg unknowns wrong type'
    
        self.inst.get().quantifyComponents((deref(unknowns.inst.get())))
    
    def optimizeCalibrationCurveIterative(self, list component_concentrations ,  feature_name ,  transformation_model , Param transformation_model_params , Param optimized_params ):
        """
        optimizeCalibrationCurveIterative(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param , optimized_params: Param ) -> bool
        """
        assert isinstance(component_concentrations, list) and all(isinstance(elemt_rec, AQS_featureConcentration) for elemt_rec in component_concentrations), 'arg component_concentrations wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert (isinstance(transformation_model, str) or isinstance(transformation_model, bytes) or isinstance(transformation_model, String)), 'arg transformation_model wrong type'
        assert isinstance(transformation_model_params, Param), 'arg transformation_model_params wrong type'
        assert isinstance(optimized_params, Param), 'arg optimized_params wrong type'
        cdef libcpp_vector[_AQS_featureConcentration] * v0 = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration item0
        for item0 in component_concentrations:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
        cdef bool _r = self.inst.get().optimizeCalibrationCurveIterative(deref(v0), deref((convString(feature_name)).get()), deref((convString(transformation_model)).get()), (deref(transformation_model_params.inst.get())), (deref(optimized_params.inst.get())))
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_component_concentrations = v0.begin()
        replace_0 = []
        while it_component_concentrations != v0.end():
            item0 = AQS_featureConcentration.__new__(AQS_featureConcentration)
            item0.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_component_concentrations)))
            replace_0.append(item0)
            inc(it_component_concentrations)
        component_concentrations[:] = replace_0
        del v0
        py_result = <bool>_r
        return py_result
    
    def optimizeSingleCalibrationCurve(self,  component_name , list component_concentrations ):
        """
        optimizeSingleCalibrationCurve(self, component_name: Union[bytes, str, String] , component_concentrations: List[AQS_featureConcentration] ) -> None
        """
        assert (isinstance(component_name, str) or isinstance(component_name, bytes) or isinstance(component_name, String)), 'arg component_name wrong type'
        assert isinstance(component_concentrations, list) and all(isinstance(elemt_rec, AQS_featureConcentration) for elemt_rec in component_concentrations), 'arg component_concentrations wrong type'
    
        cdef libcpp_vector[_AQS_featureConcentration] * v1 = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration item1
        for item1 in component_concentrations:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().optimizeSingleCalibrationCurve(deref((convString(component_name)).get()), deref(v1))
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_component_concentrations = v1.begin()
        replace_0 = []
        while it_component_concentrations != v1.end():
            item1 = AQS_featureConcentration.__new__(AQS_featureConcentration)
            item1.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_component_concentrations)))
            replace_0.append(item1)
            inc(it_component_concentrations)
        component_concentrations[:] = replace_0
        del v1
    
    def calculateBias(self, double actual_concentration , double calculated_concentration ):
        """
        calculateBias(self, actual_concentration: float , calculated_concentration: float ) -> float
        This function calculates the bias of the calibration
        """
        assert isinstance(actual_concentration, float), 'arg actual_concentration wrong type'
        assert isinstance(calculated_concentration, float), 'arg calculated_concentration wrong type'
    
    
        cdef double _r = self.inst.get().calculateBias((<double>actual_concentration), (<double>calculated_concentration))
        py_result = <double>_r
        return py_result
    
    def fitCalibration(self, list component_concentrations ,  feature_name ,  transformation_model , Param transformation_model_params ):
        """
        fitCalibration(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param ) -> Param
        """
        assert isinstance(component_concentrations, list) and all(isinstance(elemt_rec, AQS_featureConcentration) for elemt_rec in component_concentrations), 'arg component_concentrations wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert (isinstance(transformation_model, str) or isinstance(transformation_model, bytes) or isinstance(transformation_model, String)), 'arg transformation_model wrong type'
        assert isinstance(transformation_model_params, Param), 'arg transformation_model_params wrong type'
        cdef libcpp_vector[_AQS_featureConcentration] * v0 = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration item0
        for item0 in component_concentrations:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        cdef _Param * _r = new _Param(self.inst.get().fitCalibration(deref(v0), deref((convString(feature_name)).get()), deref((convString(transformation_model)).get()), (deref(transformation_model_params.inst.get()))))
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_component_concentrations = v0.begin()
        replace_0 = []
        while it_component_concentrations != v0.end():
            item0 = AQS_featureConcentration.__new__(AQS_featureConcentration)
            item0.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_component_concentrations)))
            replace_0.append(item0)
            inc(it_component_concentrations)
        component_concentrations[:] = replace_0
        del v0
        cdef Param py_result = Param.__new__(Param)
        py_result.inst = shared_ptr[_Param](_r)
        return py_result
    
    def calculateBiasAndR(self, list component_concentrations ,  feature_name ,  transformation_model , Param transformation_model_params , list biases , double correlation_coefficient ):
        """
        calculateBiasAndR(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param , biases: List[float] , correlation_coefficient: float ) -> None
        """
        assert isinstance(component_concentrations, list) and all(isinstance(elemt_rec, AQS_featureConcentration) for elemt_rec in component_concentrations), 'arg component_concentrations wrong type'
        assert (isinstance(feature_name, str) or isinstance(feature_name, bytes) or isinstance(feature_name, String)), 'arg feature_name wrong type'
        assert (isinstance(transformation_model, str) or isinstance(transformation_model, bytes) or isinstance(transformation_model, String)), 'arg transformation_model wrong type'
        assert isinstance(transformation_model_params, Param), 'arg transformation_model_params wrong type'
        assert isinstance(biases, list) and all(isinstance(elemt_rec, float) for elemt_rec in biases), 'arg biases wrong type'
        assert isinstance(correlation_coefficient, float), 'arg correlation_coefficient wrong type'
        cdef libcpp_vector[_AQS_featureConcentration] * v0 = new libcpp_vector[_AQS_featureConcentration]()
        cdef AQS_featureConcentration item0
        for item0 in component_concentrations:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        cdef libcpp_vector[double] v4 = biases
    
        self.inst.get().calculateBiasAndR(deref(v0), deref((convString(feature_name)).get()), deref((convString(transformation_model)).get()), (deref(transformation_model_params.inst.get())), v4, (<double &>correlation_coefficient))
        biases[:] = v4
        cdef libcpp_vector[_AQS_featureConcentration].iterator it_component_concentrations = v0.begin()
        replace_0 = []
        while it_component_concentrations != v0.end():
            item0 = AQS_featureConcentration.__new__(AQS_featureConcentration)
            item0.inst = shared_ptr[_AQS_featureConcentration](new _AQS_featureConcentration(deref(it_component_concentrations)))
            replace_0.append(item0)
            inc(it_component_concentrations)
        component_concentrations[:] = replace_0
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

cdef class ChromatogramSettings:
    """
    Cython implementation of _ChromatogramSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramSettings.html>`_
      -- Inherits from ['MetaInfoInterface']

    Description of the chromatogram settings, provides meta-information
    about a single chromatogram.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ChromatogramSettings rv = ChromatogramSettings.__new__(ChromatogramSettings)
       rv.inst = shared_ptr[_ChromatogramSettings](new _ChromatogramSettings(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ChromatogramSettings rv = ChromatogramSettings.__new__(ChromatogramSettings)
       rv.inst = shared_ptr[_ChromatogramSettings](new _ChromatogramSettings(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ChromatogramSettings](new _ChromatogramSettings())
    
    def _init_1(self, ChromatogramSettings in_0 ):
        """
        _init_1(self, in_0: ChromatogramSettings ) -> None
        """
        assert isinstance(in_0, ChromatogramSettings), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ChromatogramSettings](new _ChromatogramSettings((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ChromatogramSettings ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ChromatogramSettings)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getProduct(self):
        """
        getProduct(self) -> Product
        Returns the product ion
        """
        cdef _Product * _r = new _Product(self.inst.get().getProduct())
        cdef Product py_result = Product.__new__(Product)
        py_result.inst = shared_ptr[_Product](_r)
        return py_result
    
    def setProduct(self, Product p ):
        """
        setProduct(self, p: Product ) -> None
        Sets the product ion
        """
        assert isinstance(p, Product), 'arg p wrong type'
    
        self.inst.get().setProduct((deref(p.inst.get())))
    
    def getNativeID(self):
        """
        getNativeID(self) -> Union[bytes, str, String]
        Returns the native identifier for the spectrum, used by the acquisition software.
        """
        cdef _String _r = self.inst.get().getNativeID()
        py_result = convOutputString(_r)
        return py_result
    
    def setNativeID(self,  native_id ):
        """
        setNativeID(self, native_id: Union[bytes, str, String] ) -> None
        Sets the native identifier for the spectrum, used by the acquisition software.
        """
        assert (isinstance(native_id, str) or isinstance(native_id, bytes) or isinstance(native_id, String)), 'arg native_id wrong type'
    
        self.inst.get().setNativeID(deref((convString(native_id)).get()))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the free-text comment
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  comment ):
        """
        setComment(self, comment: Union[bytes, str, String] ) -> None
        Sets the free-text comment
        """
        assert (isinstance(comment, str) or isinstance(comment, bytes) or isinstance(comment, String)), 'arg comment wrong type'
    
        self.inst.get().setComment(deref((convString(comment)).get()))
    
    def getInstrumentSettings(self):
        """
        getInstrumentSettings(self) -> InstrumentSettings
        Returns the instrument settings of the current spectrum
        """
        cdef _InstrumentSettings * _r = new _InstrumentSettings(self.inst.get().getInstrumentSettings())
        cdef InstrumentSettings py_result = InstrumentSettings.__new__(InstrumentSettings)
        py_result.inst = shared_ptr[_InstrumentSettings](_r)
        return py_result
    
    def setInstrumentSettings(self, InstrumentSettings instrument_settings ):
        """
        setInstrumentSettings(self, instrument_settings: InstrumentSettings ) -> None
        Sets the instrument settings of the current spectrum
        """
        assert isinstance(instrument_settings, InstrumentSettings), 'arg instrument_settings wrong type'
    
        self.inst.get().setInstrumentSettings((deref(instrument_settings.inst.get())))
    
    def getAcquisitionInfo(self):
        """
        getAcquisitionInfo(self) -> AcquisitionInfo
        Returns the acquisition info
        """
        cdef _AcquisitionInfo * _r = new _AcquisitionInfo(self.inst.get().getAcquisitionInfo())
        cdef AcquisitionInfo py_result = AcquisitionInfo.__new__(AcquisitionInfo)
        py_result.inst = shared_ptr[_AcquisitionInfo](_r)
        return py_result
    
    def setAcquisitionInfo(self, AcquisitionInfo acquisition_info ):
        """
        setAcquisitionInfo(self, acquisition_info: AcquisitionInfo ) -> None
        Sets the acquisition info
        """
        assert isinstance(acquisition_info, AcquisitionInfo), 'arg acquisition_info wrong type'
    
        self.inst.get().setAcquisitionInfo((deref(acquisition_info.inst.get())))
    
    def getSourceFile(self):
        """
        getSourceFile(self) -> SourceFile
        Returns the source file
        """
        cdef _SourceFile * _r = new _SourceFile(self.inst.get().getSourceFile())
        cdef SourceFile py_result = SourceFile.__new__(SourceFile)
        py_result.inst = shared_ptr[_SourceFile](_r)
        return py_result
    
    def setSourceFile(self, SourceFile source_file ):
        """
        setSourceFile(self, source_file: SourceFile ) -> None
        Sets the source file
        """
        assert isinstance(source_file, SourceFile), 'arg source_file wrong type'
    
        self.inst.get().setSourceFile((deref(source_file.inst.get())))
    
    def getPrecursor(self):
        """
        getPrecursor(self) -> Precursor
        Returns the precursors
        """
        cdef _Precursor * _r = new _Precursor(self.inst.get().getPrecursor())
        cdef Precursor py_result = Precursor.__new__(Precursor)
        py_result.inst = shared_ptr[_Precursor](_r)
        return py_result
    
    def setPrecursor(self, Precursor precursor ):
        """
        setPrecursor(self, precursor: Precursor ) -> None
        Sets the precursors
        """
        assert isinstance(precursor, Precursor), 'arg precursor wrong type'
    
        self.inst.get().setPrecursor((deref(precursor.inst.get())))
    
    def getDataProcessing(self):
        """
        getDataProcessing(self) -> List[DataProcessing]
        Returns the description of the applied processing
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
        
    
    def setChromatogramType(self, int type ):
        """
        setChromatogramType(self, type: int ) -> None
        Sets the chromatogram type
        """
        assert type in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 'arg type wrong type'
    
        self.inst.get().setChromatogramType((<_ChromatogramType>type))
    
    def getChromatogramType(self):
        """
        getChromatogramType(self) -> int
        Get the chromatogram type
        """
        cdef _ChromatogramType _r = self.inst.get().getChromatogramType()
        py_result = <int>_r
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
        if not isinstance(other, ChromatogramSettings):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ChromatogramSettings other_casted = other
        cdef ChromatogramSettings self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get())
    ChromatogramType = __ChromatogramType 

cdef class CubicSpline2d:
    """
    Cython implementation of _CubicSpline2d

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CubicSpline2d.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef CubicSpline2d rv = CubicSpline2d.__new__(CubicSpline2d)
       rv.inst = shared_ptr[_CubicSpline2d](new _CubicSpline2d(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef CubicSpline2d rv = CubicSpline2d.__new__(CubicSpline2d)
       rv.inst = shared_ptr[_CubicSpline2d](new _CubicSpline2d(deref(self.inst.get())))
       return rv
    
    def _init_0(self, list x , list y ):
        """
        _init_0(self, x: List[float] , y: List[float] ) -> None
        """
        assert isinstance(x, list) and all(isinstance(elemt_rec, float) for elemt_rec in x), 'arg x wrong type'
        assert isinstance(y, list) and all(isinstance(elemt_rec, float) for elemt_rec in y), 'arg y wrong type'
        cdef libcpp_vector[double] v0 = x
        cdef libcpp_vector[double] v1 = y
        self.inst = shared_ptr[_CubicSpline2d](new _CubicSpline2d(v0, v1))
        
        
    
    def _init_1(self, CubicSpline2d in_0 ):
        """
        _init_1(self, in_0: CubicSpline2d ) -> None
        """
        assert isinstance(in_0, CubicSpline2d), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_CubicSpline2d](new _CubicSpline2d((deref(in_0.inst.get()))))
    
    def _init_2(self, dict m ):
        """
        _init_2(self, m: Dict[float, float] ) -> None
        """
        assert isinstance(m, dict) and all(isinstance(k, float) for k in m.keys()) and all(isinstance(v, float) for v in m.values()), 'arg m wrong type'
        cdef libcpp_map[double, double] * v0 = new libcpp_map[double, double]()
        for key, value in m.items():
        
        
            deref(v0)[ (<double>key) ] = (<double>value)
        
        
        self.inst = shared_ptr[_CubicSpline2d](new _CubicSpline2d(deref(v0)))
        del v0
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, x: List[float] , y: List[float] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: CubicSpline2d ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, m: Dict[float, float] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[0])) and (isinstance(args[1], list) and all(isinstance(elemt_rec, float) for elemt_rec in args[1])):
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], CubicSpline2d)):
             self._init_1(*args)
        elif (len(args)==1) and (isinstance(args[0], dict) and all(isinstance(k, float) for k in args[0].keys()) and all(isinstance(v, float) for v in args[0].values())):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def eval(self, double x ):
        """
        eval(self, x: float ) -> float
        Evaluates the cubic spline
        """
        assert isinstance(x, float), 'arg x wrong type'
    
        cdef double _r = self.inst.get().eval((<double>x))
        py_result = <double>_r
        return py_result
    
    def derivatives(self, double x ,  order ):
        """
        derivatives(self, x: float , order: int ) -> float
        Returns first, second or third derivative of cubic spline
        """
        assert isinstance(x, float), 'arg x wrong type'
        assert isinstance(order, int), 'arg order wrong type'
    
    
        cdef double _r = self.inst.get().derivatives((<double>x), (<unsigned int>order))
        py_result = <double>_r
        return py_result 

cdef class __DigestionFilter:
    """
    Cython implementation of _DigestionFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionFilter.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    property digestion_:
        def __set__(self, ProteaseDigestion digestion_):
        
            self.inst.get().digestion_ = (deref(digestion_.inst.get()))
        
    
        def __get__(self):
            cdef _ProteaseDigestion * _r = new _ProteaseDigestion(self.inst.get().digestion_)
            cdef ProteaseDigestion py_result = ProteaseDigestion.__new__(ProteaseDigestion)
            py_result.inst = shared_ptr[_ProteaseDigestion](_r)
            return py_result
    
    property ignore_missed_cleavages_:
        def __set__(self, bool ignore_missed_cleavages_):
        
            self.inst.get().ignore_missed_cleavages_ = (<bool>ignore_missed_cleavages_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().ignore_missed_cleavages_
            py_result = <bool>_r
            return py_result
    
    property methionine_cleavage_:
        def __set__(self, bool methionine_cleavage_):
        
            self.inst.get().methionine_cleavage_ = (<bool>methionine_cleavage_)
        
    
        def __get__(self):
            cdef bool _r = self.inst.get().methionine_cleavage_
            py_result = <bool>_r
            return py_result
    
    def __init__(self, list entries , ProteaseDigestion digestion , bool ignore_missed_cleavages , bool methionine_cleavage ):
        """
        __init__(self, entries: List[FASTAEntry] , digestion: ProteaseDigestion , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> None
        """
        assert isinstance(entries, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in entries), 'arg entries wrong type'
        assert isinstance(digestion, ProteaseDigestion), 'arg digestion wrong type'
        assert isinstance(ignore_missed_cleavages, pybool_t), 'arg ignore_missed_cleavages wrong type'
        assert isinstance(methionine_cleavage, pybool_t), 'arg methionine_cleavage wrong type'
        cdef libcpp_vector[_FASTAEntry] * v0 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item0
        for item0 in entries:
            v0.push_back(deref(item0.inst.get()))
    
    
    
        self.inst = shared_ptr[_DigestionFilter](new _DigestionFilter(deref(v0), (deref(digestion.inst.get())), (<bool>ignore_missed_cleavages), (<bool>methionine_cleavage)))
        cdef libcpp_vector[_FASTAEntry].iterator it_entries = v0.begin()
        replace_0 = []
        while it_entries != v0.end():
            item0 = FASTAEntry.__new__(FASTAEntry)
            item0.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_entries)))
            replace_0.append(item0)
            inc(it_entries)
        entries[:] = replace_0
        del v0
    
    def filterPeptideEvidences(self, PeptideIdentificationList peptides ):
        """
        filterPeptideEvidences(self, peptides: PeptideIdentificationList ) -> None
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().filterPeptideEvidences((deref(peptides.inst.get()))) 

cdef class DistanceMatrix:
    """
    Cython implementation of _DistanceMatrix[float]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DistanceMatrix[float].html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef DistanceMatrix rv = DistanceMatrix.__new__(DistanceMatrix)
       rv.inst = shared_ptr[_DistanceMatrix[float]](new _DistanceMatrix[float](deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef DistanceMatrix rv = DistanceMatrix.__new__(DistanceMatrix)
       rv.inst = shared_ptr[_DistanceMatrix[float]](new _DistanceMatrix[float](deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_DistanceMatrix[float]](new _DistanceMatrix[float]())
    
    def _init_1(self, DistanceMatrix in_0 ):
        """
        _init_1(self, in_0: DistanceMatrix ) -> None
        """
        assert isinstance(in_0, DistanceMatrix), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_DistanceMatrix[float]](new _DistanceMatrix[float]((deref(in_0.inst.get()))))
    
    def _init_2(self,  dimensionsize , float value ):
        """
        _init_2(self, dimensionsize: int , value: float ) -> None
        """
        assert isinstance(dimensionsize, int) and dimensionsize >= 0, 'arg dimensionsize wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
        self.inst = shared_ptr[_DistanceMatrix[float]](new _DistanceMatrix[float]((<size_t>dimensionsize), (<float>value)))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: DistanceMatrix ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, dimensionsize: int , value: float ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], DistanceMatrix)):
             self._init_1(*args)
        elif (len(args)==2) and (isinstance(args[0], int) and args[0] >= 0) and (isinstance(args[1], float)):
             self._init_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getValue(self,  i ,  j ):
        """
        getValue(self, i: int , j: int ) -> float
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
    
    
        cdef float _r = self.inst.get().getValue((<size_t>i), (<size_t>j))
        py_result = <float>_r
        return py_result
    
    def setValue(self,  i ,  j , float value ):
        """
        setValue(self, i: int , j: int , value: float ) -> None
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
    
        self.inst.get().setValue((<size_t>i), (<size_t>j), (<float>value))
    
    def setValueQuick(self,  i ,  j , float value ):
        """
        setValueQuick(self, i: int , j: int , value: float ) -> None
        """
        assert isinstance(i, int) and i >= 0, 'arg i wrong type'
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
    
        self.inst.get().setValueQuick((<size_t>i), (<size_t>j), (<float>value))
    
    def clear(self):
        """
        clear(self) -> None
        """
        self.inst.get().clear()
    
    def resize(self,  dimensionsize , float value ):
        """
        resize(self, dimensionsize: int , value: float ) -> None
        """
        assert isinstance(dimensionsize, int) and dimensionsize >= 0, 'arg dimensionsize wrong type'
        assert isinstance(value, float), 'arg value wrong type'
    
    
        self.inst.get().resize((<size_t>dimensionsize), (<float>value))
    
    def reduce(self,  j ):
        """
        reduce(self, j: int ) -> None
        """
        assert isinstance(j, int) and j >= 0, 'arg j wrong type'
    
        self.inst.get().reduce((<size_t>j))
    
    def dimensionsize(self):
        """
        dimensionsize(self) -> int
        """
        cdef size_t _r = self.inst.get().dimensionsize()
        py_result = <size_t>_r
        return py_result
    
    def updateMinElement(self):
        """
        updateMinElement(self) -> None
        """
        self.inst.get().updateMinElement()
    
    def getMinElementCoordinates(self):
        """
        getMinElementCoordinates(self) -> List[int, int]
        """
        _r = self.inst.get().getMinElementCoordinates()
        cdef list py_result = [_r.first, _r.second]
        return py_result
    
    def __richcmp__(self, other, op):
        if op not in (2,):
           op_str = {0: '<', 2: '==', 4: '>', 1: '<=', 3: '!=', 5: '>='}[op]
           raise NotImplementedError("Comparison operator %s not implemented" % op_str)
        if not isinstance(other, DistanceMatrix):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef DistanceMatrix other_casted = other
        cdef DistanceMatrix self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get()) 

cdef class ExperimentalSettings:
    """
    Cython implementation of _ExperimentalSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalSettings.html>`_
      -- Inherits from ['DocumentIdentifier', 'MetaInfoInterface']

    Description of the experimental settings, provides meta-information
    about an LC-MS/MS injection.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ExperimentalSettings rv = ExperimentalSettings.__new__(ExperimentalSettings)
       rv.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ExperimentalSettings rv = ExperimentalSettings.__new__(ExperimentalSettings)
       rv.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings())
    
    def _init_1(self, ExperimentalSettings in_0 ):
        """
        _init_1(self, in_0: ExperimentalSettings ) -> None
        """
        assert isinstance(in_0, ExperimentalSettings), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ExperimentalSettings](new _ExperimentalSettings((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ExperimentalSettings ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ExperimentalSettings)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getSourceFiles(self):
        """
        getSourceFiles(self) -> List[SourceFile]
        Returns a reference to the source data file
        """
        _r = self.inst.get().getSourceFiles()
        py_result = []
        cdef libcpp_vector[_SourceFile].iterator it__r = _r.begin()
        cdef SourceFile item_py_result
        while it__r != _r.end():
           item_py_result = SourceFile.__new__(SourceFile)
           item_py_result.inst = shared_ptr[_SourceFile](new _SourceFile(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setSourceFiles(self, list source_files ):
        """
        setSourceFiles(self, source_files: List[SourceFile] ) -> None
        Sets the source data file
        """
        assert isinstance(source_files, list) and all(isinstance(elemt_rec, SourceFile) for elemt_rec in source_files), 'arg source_files wrong type'
        cdef libcpp_vector[_SourceFile] * v0 = new libcpp_vector[_SourceFile]()
        cdef SourceFile item0
        for item0 in source_files:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setSourceFiles(deref(v0))
        del v0
    
    def getDateTime(self):
        """
        getDateTime(self) -> DateTime
        Returns the date the experiment was performed
        """
        cdef _DateTime * _r = new _DateTime(self.inst.get().getDateTime())
        cdef DateTime py_result = DateTime.__new__(DateTime)
        py_result.inst = shared_ptr[_DateTime](_r)
        return py_result
    
    def setDateTime(self, DateTime date_time ):
        """
        setDateTime(self, date_time: DateTime ) -> None
        Sets the date the experiment was performed
        """
        assert isinstance(date_time, DateTime), 'arg date_time wrong type'
    
        self.inst.get().setDateTime((deref(date_time.inst.get())))
    
    def getSample(self):
        """
        getSample(self) -> Sample
        Returns a reference to the sample description
        """
        cdef _Sample * _r = new _Sample(self.inst.get().getSample())
        cdef Sample py_result = Sample.__new__(Sample)
        py_result.inst = shared_ptr[_Sample](_r)
        return py_result
    
    def setSample(self, Sample sample ):
        """
        setSample(self, sample: Sample ) -> None
        Sets the sample description
        """
        assert isinstance(sample, Sample), 'arg sample wrong type'
    
        self.inst.get().setSample((deref(sample.inst.get())))
    
    def getContacts(self):
        """
        getContacts(self) -> List[ContactPerson]
        Returns a reference to the list of contact persons
        """
        _r = self.inst.get().getContacts()
        py_result = []
        cdef libcpp_vector[_ContactPerson].iterator it__r = _r.begin()
        cdef ContactPerson item_py_result
        while it__r != _r.end():
           item_py_result = ContactPerson.__new__(ContactPerson)
           item_py_result.inst = shared_ptr[_ContactPerson](new _ContactPerson(deref(it__r)))
           py_result.append(item_py_result)
           inc(it__r)
        return py_result
    
    def setContacts(self, list contacts ):
        """
        setContacts(self, contacts: List[ContactPerson] ) -> None
        Sets the list of contact persons
        """
        assert isinstance(contacts, list) and all(isinstance(elemt_rec, ContactPerson) for elemt_rec in contacts), 'arg contacts wrong type'
        cdef libcpp_vector[_ContactPerson] * v0 = new libcpp_vector[_ContactPerson]()
        cdef ContactPerson item0
        for item0 in contacts:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().setContacts(deref(v0))
        del v0
    
    def getInstrument(self):
        """
        getInstrument(self) -> Instrument
        Returns a reference to the MS instrument description
        """
        cdef _Instrument * _r = new _Instrument(self.inst.get().getInstrument())
        cdef Instrument py_result = Instrument.__new__(Instrument)
        py_result.inst = shared_ptr[_Instrument](_r)
        return py_result
    
    def setInstrument(self, Instrument instrument ):
        """
        setInstrument(self, instrument: Instrument ) -> None
        Sets the MS instrument description
        """
        assert isinstance(instrument, Instrument), 'arg instrument wrong type'
    
        self.inst.get().setInstrument((deref(instrument.inst.get())))
    
    def getHPLC(self):
        """
        getHPLC(self) -> HPLC
        Returns a reference to the description of the HPLC run
        """
        cdef _HPLC * _r = new _HPLC(self.inst.get().getHPLC())
        cdef HPLC py_result = HPLC.__new__(HPLC)
        py_result.inst = shared_ptr[_HPLC](_r)
        return py_result
    
    def setHPLC(self, HPLC hplc ):
        """
        setHPLC(self, hplc: HPLC ) -> None
        Sets the description of the HPLC run
        """
        assert isinstance(hplc, HPLC), 'arg hplc wrong type'
    
        self.inst.get().setHPLC((deref(hplc.inst.get())))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the free-text comment
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  comment ):
        """
        setComment(self, comment: Union[bytes, str, String] ) -> None
        Sets the free-text comment
        """
        assert (isinstance(comment, str) or isinstance(comment, bytes) or isinstance(comment, String)), 'arg comment wrong type'
    
        self.inst.get().setComment(deref((convString(comment)).get()))
    
    def getFractionIdentifier(self):
        """
        getFractionIdentifier(self) -> Union[bytes, str, String]
        Returns fraction identifier
        """
        cdef _String _r = self.inst.get().getFractionIdentifier()
        py_result = convOutputString(_r)
        return py_result
    
    def setFractionIdentifier(self,  fraction_identifier ):
        """
        setFractionIdentifier(self, fraction_identifier: Union[bytes, str, String] ) -> None
        Sets the fraction identifier
        """
        assert (isinstance(fraction_identifier, str) or isinstance(fraction_identifier, bytes) or isinstance(fraction_identifier, String)), 'arg fraction_identifier wrong type'
    
        self.inst.get().setFractionIdentifier(deref((convString(fraction_identifier)).get()))
    
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
        if not isinstance(other, ExperimentalSettings):
           raise NotImplementedError("Comparison currently only allowed with objects of the same type. Use isinstance and define yourself.")
        cdef ExperimentalSettings other_casted = other
        cdef ExperimentalSettings self_casted = self
        if op==2:
            return deref(self_casted.inst.get()) == deref(other_casted.inst.get())
        if op==3:
            return deref(self_casted.inst.get()) != deref(other_casted.inst.get()) 

cdef class HPLC:
    """
    Cython implementation of _HPLC

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1HPLC.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef HPLC rv = HPLC.__new__(HPLC)
       rv.inst = shared_ptr[_HPLC](new _HPLC(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef HPLC rv = HPLC.__new__(HPLC)
       rv.inst = shared_ptr[_HPLC](new _HPLC(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        Representation of a HPLC experiment
        """
        self.inst = shared_ptr[_HPLC](new _HPLC())
    
    def _init_1(self, HPLC in_0 ):
        """
        _init_1(self, in_0: HPLC ) -> None
        """
        assert isinstance(in_0, HPLC), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_HPLC](new _HPLC((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        Representation of a HPLC experiment

        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: HPLC ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], HPLC)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def getInstrument(self):
        """
        getInstrument(self) -> Union[bytes, str, String]
        Returns a reference to the instument name
        """
        cdef _String _r = self.inst.get().getInstrument()
        py_result = convOutputString(_r)
        return py_result
    
    def setInstrument(self,  instrument ):
        """
        setInstrument(self, instrument: Union[bytes, str, String] ) -> None
        Sets the instument name
        """
        assert (isinstance(instrument, str) or isinstance(instrument, bytes) or isinstance(instrument, String)), 'arg instrument wrong type'
    
        self.inst.get().setInstrument(deref((convString(instrument)).get()))
    
    def getColumn(self):
        """
        getColumn(self) -> Union[bytes, str, String]
        Returns a reference to the column description
        """
        cdef _String _r = self.inst.get().getColumn()
        py_result = convOutputString(_r)
        return py_result
    
    def setColumn(self,  column ):
        """
        setColumn(self, column: Union[bytes, str, String] ) -> None
        Sets the column description
        """
        assert (isinstance(column, str) or isinstance(column, bytes) or isinstance(column, String)), 'arg column wrong type'
    
        self.inst.get().setColumn(deref((convString(column)).get()))
    
    def getTemperature(self):
        """
        getTemperature(self) -> int
        Returns the temperature (in degree C)
        """
        cdef int _r = self.inst.get().getTemperature()
        py_result = <int>_r
        return py_result
    
    def setTemperature(self,  temperature ):
        """
        setTemperature(self, temperature: int ) -> None
        Sets the temperature (in degree C)
        """
        assert isinstance(temperature, int), 'arg temperature wrong type'
    
        self.inst.get().setTemperature((<int>temperature))
    
    def getPressure(self):
        """
        getPressure(self) -> int
        Returns the pressure (in bar)
        """
        cdef unsigned int _r = self.inst.get().getPressure()
        py_result = <unsigned int>_r
        return py_result
    
    def setPressure(self,  pressure ):
        """
        setPressure(self, pressure: int ) -> None
        Sets the pressure (in bar)
        """
        assert isinstance(pressure, int), 'arg pressure wrong type'
    
        self.inst.get().setPressure((<unsigned int>pressure))
    
    def getFlux(self):
        """
        getFlux(self) -> int
        Returns the flux (in microliter/sec)
        """
        cdef unsigned int _r = self.inst.get().getFlux()
        py_result = <unsigned int>_r
        return py_result
    
    def setFlux(self,  flux ):
        """
        setFlux(self, flux: int ) -> None
        Sets the flux (in microliter/sec)
        """
        assert isinstance(flux, int), 'arg flux wrong type'
    
        self.inst.get().setFlux((<unsigned int>flux))
    
    def getComment(self):
        """
        getComment(self) -> Union[bytes, str, String]
        Returns the comments
        """
        cdef _String _r = self.inst.get().getComment()
        py_result = convOutputString(_r)
        return py_result
    
    def setComment(self,  comment ):
        """
        setComment(self, comment: Union[bytes, str, String] ) -> None
        Sets the comments
        """
        assert (isinstance(comment, str) or isinstance(comment, bytes) or isinstance(comment, String)), 'arg comment wrong type'
    
        self.inst.get().setComment(deref((convString(comment)).get()))
    
    def getGradient(self):
        """
        getGradient(self) -> Gradient
        Returns a mutable reference to the used gradient
        """
        cdef _Gradient * _r = new _Gradient(self.inst.get().getGradient())
        cdef Gradient py_result = Gradient.__new__(Gradient)
        py_result.inst = shared_ptr[_Gradient](_r)
        return py_result
    
    def setGradient(self, Gradient gradient ):
        """
        setGradient(self, gradient: Gradient ) -> None
        Sets the used gradient
        """
        assert isinstance(gradient, Gradient), 'arg gradient wrong type'
    
        self.inst.get().setGradient((deref(gradient.inst.get()))) 

cdef class IDFilter:
    """
    Cython implementation of _IDFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDFilter.html>`_

    Finds the best-scoring hit in a vector of peptide or protein identifications\n
    
    This class provides functions for filtering collections of peptide or protein identifications according to various criteria.
    It also contains helper functions and classes (functors that implement predicates) that are used in this context.\n
    
    The filter functions modify their inputs, rather than creating filtered copies.\n
    
    Most filters work on the hit level, i.e. they remove peptide or protein hits from peptide or protein identifications (IDs).
    A few filters work on the ID level instead, i.e. they remove peptide or protein IDs from vectors thereof.
    Independent of this, the inputs for all filter functions are vectors of IDs, because the data most often comes in this form.
    This design also allows many helper objects to be set up only once per vector, rather than once per ID.\n
    
    The filter functions for vectors of peptide/protein IDs do not include clean-up steps (e.g. removal of IDs without hits, reassignment of hit ranks, ...).
    They only carry out their specific filtering operations.
    This is so filters can be chained without having to repeat clean-up operations.
    The group of clean-up functions provides helpers that are useful to ensure data integrity after filters have been applied, but it is up to the individual developer to use them when necessary.\n
    
    The filter functions for MS/MS experiments do include clean-up steps, because they filter peptide and protein IDs in conjunction and potential contradictions between the two must be eliminated.
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef IDFilter rv = IDFilter.__new__(IDFilter)
       rv.inst = shared_ptr[_IDFilter](new _IDFilter(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef IDFilter rv = IDFilter.__new__(IDFilter)
       rv.inst = shared_ptr[_IDFilter](new _IDFilter(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_IDFilter](new _IDFilter())
    
    def _init_1(self, IDFilter in_0 ):
        """
        _init_1(self, in_0: IDFilter ) -> None
        """
        assert isinstance(in_0, IDFilter), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_IDFilter](new _IDFilter((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: IDFilter ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], IDFilter)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _countHits_0(self, PeptideIdentificationList identifications ):
        """
        _countHits_0(self, identifications: PeptideIdentificationList ) -> int
        Returns the total number of peptide hits in a vector of peptide identifications
        """
        assert isinstance(identifications, PeptideIdentificationList), 'arg identifications wrong type'
    
        cdef size_t _r = self.inst.get().countHits((deref(identifications.inst.get())))
        py_result = <size_t>_r
        return py_result
    
    def _countHits_1(self, list identifications ):
        """
        _countHits_1(self, identifications: List[ProteinIdentification] ) -> int
        Returns the total number of protein hits in a vector of protein identifications
        """
        assert isinstance(identifications, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in identifications), 'arg identifications wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in identifications:
            v0.push_back(deref(item0.inst.get()))
        cdef size_t _r = self.inst.get().countHits(deref(v0))
        del v0
        py_result = <size_t>_r
        return py_result
    
    def countHits(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: countHits(self, identifications: PeptideIdentificationList ) -> int
          :noindex:
        
        Returns the total number of peptide hits in a vector of peptide identifications

        
        .. rubric:: Overload:
        .. py:function:: countHits(self, identifications: List[ProteinIdentification] ) -> int
          :noindex:
        
        Returns the total number of protein hits in a vector of protein identifications
    
        """
        if (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._countHits_0(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])):
            return self._countHits_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _getBestHit_0(self, PeptideIdentificationList identifications , bool assume_sorted , PeptideHit best_hit ):
        """
        _getBestHit_0(self, identifications: PeptideIdentificationList , assume_sorted: bool , best_hit: PeptideHit ) -> bool
        Finds the best-scoring hit in a vector of peptide or protein identifications\n
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of peptide identifications
        :return: true if a hit was present, false otherwise
        """
        assert isinstance(identifications, PeptideIdentificationList), 'arg identifications wrong type'
        assert isinstance(assume_sorted, pybool_t), 'arg assume_sorted wrong type'
        assert isinstance(best_hit, PeptideHit), 'arg best_hit wrong type'
    
    
    
        cdef bool _r = self.inst.get().getBestHit((deref(identifications.inst.get())), (<bool>assume_sorted), (deref(best_hit.inst.get())))
        py_result = <bool>_r
        return py_result
    
    def _getBestHit_1(self, list identifications , bool assume_sorted , ProteinHit best_hit ):
        """
        _getBestHit_1(self, identifications: List[ProteinIdentification] , assume_sorted: bool , best_hit: ProteinHit ) -> bool
        Finds the best-scoring hit in a vector of peptide or protein identifications
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of protein identifications
        :return: true if a hit was present, false otherwise
        """
        assert isinstance(identifications, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in identifications), 'arg identifications wrong type'
        assert isinstance(assume_sorted, pybool_t), 'arg assume_sorted wrong type'
        assert isinstance(best_hit, ProteinHit), 'arg best_hit wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in identifications:
            v0.push_back(deref(item0.inst.get()))
    
    
        cdef bool _r = self.inst.get().getBestHit(deref(v0), (<bool>assume_sorted), (deref(best_hit.inst.get())))
        del v0
        py_result = <bool>_r
        return py_result
    
    def getBestHit(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: getBestHit(self, identifications: PeptideIdentificationList , assume_sorted: bool , best_hit: PeptideHit ) -> bool
          :noindex:
        
        Finds the best-scoring hit in a vector of peptide or protein identifications\n
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of peptide identifications
        :return: true if a hit was present, false otherwise
        
        .. rubric:: Overload:
        .. py:function:: getBestHit(self, identifications: List[ProteinIdentification] , assume_sorted: bool , best_hit: ProteinHit ) -> bool
          :noindex:
        
        Finds the best-scoring hit in a vector of peptide or protein identifications
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of protein identifications
        :return: true if a hit was present, false otherwise
    
        """
        if (len(args)==3) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], PeptideHit)):
            return self._getBestHit_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], pybool_t)) and (isinstance(args[2], ProteinHit)):
            return self._getBestHit_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def extractPeptideSequences(self, PeptideIdentificationList peptides , set sequences , bool ignore_mods ):
        """
        extractPeptideSequences(self, peptides: PeptideIdentificationList , sequences: Set[bytes] , ignore_mods: bool ) -> None
        Extracts all unique peptide sequences from a list of peptide IDs
        
        
        :param peptides:
        :param ignore_mods: Boolean operator default to false in case of any modifications in sequences during extraction
        :return: Sequences
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(sequences, set) and all(isinstance(i, bytes) for i in sequences), 'arg sequences wrong type'
        assert isinstance(ignore_mods, pybool_t), 'arg ignore_mods wrong type'
    
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in sequences:
           v1.insert(_String(<char *>item1))
    
        self.inst.get().extractPeptideSequences((deref(peptides.inst.get())), deref(v1), (<bool>ignore_mods))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        sequences.clear()
        sequences.update(replace)
        del v1
    
    def removeUnreferencedProteins(self, list proteins , PeptideIdentificationList peptides ):
        """
        removeUnreferencedProteins(self, proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList ) -> None
        Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list
        """
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in proteins:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().removeUnreferencedProteins(deref(v0), (deref(peptides.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v0.begin()
        replace_0 = []
        while it_proteins != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item0)
            inc(it_proteins)
        proteins[:] = replace_0
        del v0
    
    def updateProteinReferences(self, PeptideIdentificationList peptides , list proteins , bool remove_peptides_without_reference ):
        """
        updateProteinReferences(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , remove_peptides_without_reference: bool ) -> None
        Removes references to missing proteins. Only PeptideEvidence entries that reference protein hits in 'proteins' are kept in the peptide hits
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in proteins), 'arg proteins wrong type'
        assert isinstance(remove_peptides_without_reference, pybool_t), 'arg remove_peptides_without_reference wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in proteins:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().updateProteinReferences((deref(peptides.inst.get())), deref(v1), (<bool>remove_peptides_without_reference))
        cdef libcpp_vector[_ProteinIdentification].iterator it_proteins = v1.begin()
        replace_0 = []
        while it_proteins != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_proteins)))
            replace_0.append(item1)
            inc(it_proteins)
        proteins[:] = replace_0
        del v1
    
    def updateProteinGroups(self, list groups , list hits ):
        """
        updateProteinGroups(self, groups: List[ProteinGroup] , hits: List[ProteinHit] ) -> bool
        Update protein groups after protein hits were filtered
        
        
        :param groups: Input/output protein groups
        :param hits: Available protein hits (all others are removed from the groups)
        :return: Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed)
        """
        assert isinstance(groups, list) and all(isinstance(elemt_rec, ProteinGroup) for elemt_rec in groups), 'arg groups wrong type'
        assert isinstance(hits, list) and all(isinstance(elemt_rec, ProteinHit) for elemt_rec in hits), 'arg hits wrong type'
        cdef libcpp_vector[_ProteinGroup] * v0 = new libcpp_vector[_ProteinGroup]()
        cdef ProteinGroup item0
        for item0 in groups:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_vector[_ProteinHit] * v1 = new libcpp_vector[_ProteinHit]()
        cdef ProteinHit item1
        for item1 in hits:
            v1.push_back(deref(item1.inst.get()))
        cdef bool _r = self.inst.get().updateProteinGroups(deref(v0), deref(v1))
        cdef libcpp_vector[_ProteinHit].iterator it_hits = v1.begin()
        replace_0 = []
        while it_hits != v1.end():
            item1 = ProteinHit.__new__(ProteinHit)
            item1.inst = shared_ptr[_ProteinHit](new _ProteinHit(deref(it_hits)))
            replace_0.append(item1)
            inc(it_hits)
        hits[:] = replace_0
        del v1
        cdef libcpp_vector[_ProteinGroup].iterator it_groups = v0.begin()
        replace_0 = []
        while it_groups != v0.end():
            item0 = ProteinGroup.__new__(ProteinGroup)
            item0.inst = shared_ptr[_ProteinGroup](new _ProteinGroup(deref(it_groups)))
            replace_0.append(item0)
            inc(it_groups)
        groups[:] = replace_0
        del v0
        py_result = <bool>_r
        return py_result
    
    def _removeEmptyIdentifications_0(self, PeptideIdentificationList ids ):
        """
        _removeEmptyIdentifications_0(self, ids: PeptideIdentificationList ) -> None
        Removes peptide or protein identifications that have no hits in them
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
    
        self.inst.get().removeEmptyIdentifications((deref(ids.inst.get())))
    
    def _removeEmptyIdentifications_1(self, list ids ):
        """
        _removeEmptyIdentifications_1(self, ids: List[ProteinIdentification] ) -> None
        Removes peptide or protein identifications that have no hits in them
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().removeEmptyIdentifications(deref(v0))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def removeEmptyIdentifications(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeEmptyIdentifications(self, ids: PeptideIdentificationList ) -> None
          :noindex:
        
        Removes peptide or protein identifications that have no hits in them

        
        .. rubric:: Overload:
        .. py:function:: removeEmptyIdentifications(self, ids: List[ProteinIdentification] ) -> None
          :noindex:
        
        Removes peptide or protein identifications that have no hits in them
    
        """
        if (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._removeEmptyIdentifications_0(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])):
            return self._removeEmptyIdentifications_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _filterHitsByScore_0(self, PeptideIdentificationList ids , double threshold_score ):
        """
        _filterHitsByScore_0(self, ids: PeptideIdentificationList , threshold_score: float ) -> None
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(threshold_score, float), 'arg threshold_score wrong type'
    
    
        self.inst.get().filterHitsByScore((deref(ids.inst.get())), (<double>threshold_score))
    
    def _filterHitsByScore_1(self, list ids , double threshold_score ):
        """
        _filterHitsByScore_1(self, ids: List[ProteinIdentification] , threshold_score: float ) -> None
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification/ProteinIdentifiation. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(threshold_score, float), 'arg threshold_score wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().filterHitsByScore(deref(v0), (<double>threshold_score))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def _filterHitsByScore_2(self, AnnotatedMSRun experiment , double peptide_threshold_score , double protein_threshold_score ):
        """
        _filterHitsByScore_2(self, experiment: AnnotatedMSRun , peptide_threshold_score: float , protein_threshold_score: float ) -> None
        Filters an MS/MS experiment according to score thresholds
        """
        assert isinstance(experiment, AnnotatedMSRun), 'arg experiment wrong type'
        assert isinstance(peptide_threshold_score, float), 'arg peptide_threshold_score wrong type'
        assert isinstance(protein_threshold_score, float), 'arg protein_threshold_score wrong type'
    
    
    
        self.inst.get().filterHitsByScore((deref(experiment.inst.get())), (<double>peptide_threshold_score), (<double>protein_threshold_score))
    
    def filterHitsByScore(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: filterHitsByScore(self, ids: PeptideIdentificationList , threshold_score: float ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept

        
        .. rubric:: Overload:
        .. py:function:: filterHitsByScore(self, ids: List[ProteinIdentification] , threshold_score: float ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification/ProteinIdentifiation. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept

        
        .. rubric:: Overload:
        .. py:function:: filterHitsByScore(self, experiment: AnnotatedMSRun , peptide_threshold_score: float , protein_threshold_score: float ) -> None
          :noindex:
        
        Filters an MS/MS experiment according to score thresholds
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], float)):
            return self._filterHitsByScore_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], float)):
            return self._filterHitsByScore_1(*args)
        elif (len(args)==3) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], float)) and (isinstance(args[2], float)):
            return self._filterHitsByScore_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def keepNBestSpectra(self, PeptideIdentificationList peptides ,  n ):
        """
        keepNBestSpectra(self, peptides: PeptideIdentificationList , n: int ) -> None
        Filter identifications by "N best" PeptideIdentification objects (better PeptideIdentification means better [best] PeptideHit than other)
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
    
        self.inst.get().keepNBestSpectra((deref(peptides.inst.get())), (<size_t>n))
    
    def _keepNBestHits_0(self, PeptideIdentificationList ids ,  n ):
        """
        _keepNBestHits_0(self, ids: PeptideIdentificationList , n: int ) -> None
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
    
        self.inst.get().keepNBestHits((deref(ids.inst.get())), (<size_t>n))
    
    def _keepNBestHits_1(self, list ids ,  n ):
        """
        _keepNBestHits_1(self, ids: List[ProteinIdentification] , n: int ) -> None
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
    
        self.inst.get().keepNBestHits(deref(v0), (<size_t>n))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def _keepNBestHits_2(self, AnnotatedMSRun experiment ,  n ):
        """
        _keepNBestHits_2(self, experiment: AnnotatedMSRun , n: int ) -> None
        Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
        """
        assert isinstance(experiment, AnnotatedMSRun), 'arg experiment wrong type'
        assert isinstance(n, int) and n >= 0, 'arg n wrong type'
    
    
        self.inst.get().keepNBestHits((deref(experiment.inst.get())), (<size_t>n))
    
    def keepNBestHits(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: keepNBestHits(self, ids: PeptideIdentificationList , n: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: keepNBestHits(self, ids: List[ProteinIdentification] , n: int ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: keepNBestHits(self, experiment: AnnotatedMSRun , n: int ) -> None
          :noindex:
        
        Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], int) and args[1] >= 0):
            return self._keepNBestHits_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], int) and args[1] >= 0):
            return self._keepNBestHits_1(*args)
        elif (len(args)==2) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], int) and args[1] >= 0):
            return self._keepNBestHits_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _filterHitsByRank_0(self, PeptideIdentificationList ids ,  min_rank ,  max_rank ):
        """
        _filterHitsByRank_0(self, ids: PeptideIdentificationList , min_rank: int , max_rank: int ) -> None
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(min_rank, int) and min_rank >= 0, 'arg min_rank wrong type'
        assert isinstance(max_rank, int) and max_rank >= 0, 'arg max_rank wrong type'
    
    
    
        self.inst.get().filterHitsByRank((deref(ids.inst.get())), (<size_t>min_rank), (<size_t>max_rank))
    
    def _filterHitsByRank_1(self, list ids ,  min_rank ,  max_rank ):
        """
        _filterHitsByRank_1(self, ids: List[ProteinIdentification] , min_rank: int , max_rank: int ) -> None
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(min_rank, int) and min_rank >= 0, 'arg min_rank wrong type'
        assert isinstance(max_rank, int) and max_rank >= 0, 'arg max_rank wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
    
    
        self.inst.get().filterHitsByRank(deref(v0), (<size_t>min_rank), (<size_t>max_rank))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def filterHitsByRank(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: filterHitsByRank(self, ids: PeptideIdentificationList , min_rank: int , max_rank: int ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
        
        .. rubric:: Overload:
        .. py:function:: filterHitsByRank(self, ids: List[ProteinIdentification] , min_rank: int , max_rank: int ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
    
        """
        if (len(args)==3) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int) and args[2] >= 0):
            return self._filterHitsByRank_0(*args)
        elif (len(args)==3) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], int) and args[1] >= 0) and (isinstance(args[2], int) and args[2] >= 0):
            return self._filterHitsByRank_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _removeDecoyHits_0(self, PeptideIdentificationList ids ):
        """
        _removeDecoyHits_0(self, ids: PeptideIdentificationList ) -> None
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
    
        self.inst.get().removeDecoyHits((deref(ids.inst.get())))
    
    def _removeDecoyHits_1(self, list ids ):
        """
        _removeDecoyHits_1(self, ids: List[ProteinIdentification] ) -> None
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        self.inst.get().removeDecoyHits(deref(v0))
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def removeDecoyHits(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeDecoyHits(self, ids: PeptideIdentificationList ) -> None
          :noindex:
        
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively

        
        .. rubric:: Overload:
        .. py:function:: removeDecoyHits(self, ids: List[ProteinIdentification] ) -> None
          :noindex:
        
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
    
        """
        if (len(args)==1) and (isinstance(args[0], PeptideIdentificationList)):
            return self._removeDecoyHits_0(*args)
        elif (len(args)==1) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])):
            return self._removeDecoyHits_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _removeHitsMatchingProteins_0(self, PeptideIdentificationList ids , set accessions ):
        """
        _removeHitsMatchingProteins_0(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None
        Filters peptide or protein identifications according to the given proteins (negative)
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(accessions, set) and all(isinstance(i, bytes) for i in accessions), 'arg accessions wrong type'
    
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in accessions:
           v1.insert(_String(<char *>item1))
        self.inst.get().removeHitsMatchingProteins((deref(ids.inst.get())), deref(v1))
        del v1
    
    def _removeHitsMatchingProteins_1(self, list ids , set accessions ):
        """
        _removeHitsMatchingProteins_1(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None
        Filters peptide or protein identifications according to the given proteins (negative)
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(accessions, set) and all(isinstance(i, bytes) for i in accessions), 'arg accessions wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in accessions:
           v1.insert(_String(<char *>item1))
        self.inst.get().removeHitsMatchingProteins(deref(v0), deref(v1))
        del v1
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def removeHitsMatchingProteins(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: removeHitsMatchingProteins(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the given proteins (negative)

        
        .. rubric:: Overload:
        .. py:function:: removeHitsMatchingProteins(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the given proteins (negative)
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], set) and all(isinstance(i, bytes) for i in args[1])):
            return self._removeHitsMatchingProteins_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], set) and all(isinstance(i, bytes) for i in args[1])):
            return self._removeHitsMatchingProteins_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _keepHitsMatchingProteins_0(self, PeptideIdentificationList ids , set accessions ):
        """
        _keepHitsMatchingProteins_0(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None
        Filters peptide or protein identifications according to the given proteins (positive)
        """
        assert isinstance(ids, PeptideIdentificationList), 'arg ids wrong type'
        assert isinstance(accessions, set) and all(isinstance(i, bytes) for i in accessions), 'arg accessions wrong type'
    
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in accessions:
           v1.insert(_String(<char *>item1))
        self.inst.get().keepHitsMatchingProteins((deref(ids.inst.get())), deref(v1))
        del v1
    
    def _keepHitsMatchingProteins_1(self, list ids , set accessions ):
        """
        _keepHitsMatchingProteins_1(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None
        Filters peptide or protein identifications according to the given proteins (positive)
        """
        assert isinstance(ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in ids), 'arg ids wrong type'
        assert isinstance(accessions, set) and all(isinstance(i, bytes) for i in accessions), 'arg accessions wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in ids:
            v0.push_back(deref(item0.inst.get()))
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in accessions:
           v1.insert(_String(<char *>item1))
        self.inst.get().keepHitsMatchingProteins(deref(v0), deref(v1))
        del v1
        cdef libcpp_vector[_ProteinIdentification].iterator it_ids = v0.begin()
        replace_0 = []
        while it_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_ids)))
            replace_0.append(item0)
            inc(it_ids)
        ids[:] = replace_0
        del v0
    
    def _keepHitsMatchingProteins_2(self, AnnotatedMSRun experiment , list proteins ):
        """
        _keepHitsMatchingProteins_2(self, experiment: AnnotatedMSRun , proteins: List[FASTAEntry] ) -> None
        """
        assert isinstance(experiment, AnnotatedMSRun), 'arg experiment wrong type'
        assert isinstance(proteins, list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in proteins), 'arg proteins wrong type'
    
        cdef libcpp_vector[_FASTAEntry] * v1 = new libcpp_vector[_FASTAEntry]()
        cdef FASTAEntry item1
        for item1 in proteins:
            v1.push_back(deref(item1.inst.get()))
        self.inst.get().keepHitsMatchingProteins((deref(experiment.inst.get())), deref(v1))
        cdef libcpp_vector[_FASTAEntry].iterator it_proteins = v1.begin()
        replace_0 = []
        while it_proteins != v1.end():
            item1 = FASTAEntry.__new__(FASTAEntry)
            item1.inst = shared_ptr[_FASTAEntry](new _FASTAEntry(deref(it_proteins)))
            replace_0.append(item1)
            inc(it_proteins)
        proteins[:] = replace_0
        del v1
    
    def keepHitsMatchingProteins(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: keepHitsMatchingProteins(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the given proteins (positive)

        
        .. rubric:: Overload:
        .. py:function:: keepHitsMatchingProteins(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None
          :noindex:
        
        Filters peptide or protein identifications according to the given proteins (positive)

        
        .. rubric:: Overload:
        .. py:function:: keepHitsMatchingProteins(self, experiment: AnnotatedMSRun , proteins: List[FASTAEntry] ) -> None
          :noindex:
    
        """
        if (len(args)==2) and (isinstance(args[0], PeptideIdentificationList)) and (isinstance(args[1], set) and all(isinstance(i, bytes) for i in args[1])):
            return self._keepHitsMatchingProteins_0(*args)
        elif (len(args)==2) and (isinstance(args[0], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[0])) and (isinstance(args[1], set) and all(isinstance(i, bytes) for i in args[1])):
            return self._keepHitsMatchingProteins_1(*args)
        elif (len(args)==2) and (isinstance(args[0], AnnotatedMSRun)) and (isinstance(args[1], list) and all(isinstance(elemt_rec, FASTAEntry) for elemt_rec in args[1])):
            return self._keepHitsMatchingProteins_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def keepBestPeptideHits(self, PeptideIdentificationList peptides , bool strict ):
        """
        keepBestPeptideHits(self, peptides: PeptideIdentificationList , strict: bool ) -> None
        Filters peptide identifications keeping only the single best-scoring hit per ID
        
        
        :param peptides: Input/output
        :param strict: If set, keep the best hit only if its score is unique - i.e. ties are not allowed. (Otherwise all hits with the best score is kept.)
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(strict, pybool_t), 'arg strict wrong type'
    
    
        self.inst.get().keepBestPeptideHits((deref(peptides.inst.get())), (<bool>strict))
    
    def filterPeptidesByLength(self, PeptideIdentificationList peptides ,  min_length ,  max_length ):
        """
        filterPeptidesByLength(self, peptides: PeptideIdentificationList , min_length: int , max_length: int ) -> None
        Filters peptide identifications according to peptide sequence length
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(min_length, int) and min_length >= 0, 'arg min_length wrong type'
        assert isinstance(max_length, int) and max_length >= 0, 'arg max_length wrong type'
    
    
    
        self.inst.get().filterPeptidesByLength((deref(peptides.inst.get())), (<size_t>min_length), (<size_t>max_length))
    
    def filterPeptidesByCharge(self, PeptideIdentificationList peptides ,  min_charge ,  max_charge ):
        """
        filterPeptidesByCharge(self, peptides: PeptideIdentificationList , min_charge: int , max_charge: int ) -> None
        Filters peptide identifications according to charge state
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(min_charge, int) and min_charge >= 0, 'arg min_charge wrong type'
        assert isinstance(max_charge, int) and max_charge >= 0, 'arg max_charge wrong type'
    
    
    
        self.inst.get().filterPeptidesByCharge((deref(peptides.inst.get())), (<size_t>min_charge), (<size_t>max_charge))
    
    def filterPeptidesByRT(self, PeptideIdentificationList peptides ,  min_rt ,  max_rt ):
        """
        filterPeptidesByRT(self, peptides: PeptideIdentificationList , min_rt: int , max_rt: int ) -> None
        Filters peptide identifications by precursor RT, keeping only IDs in the given range
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(min_rt, int) and min_rt >= 0, 'arg min_rt wrong type'
        assert isinstance(max_rt, int) and max_rt >= 0, 'arg max_rt wrong type'
    
    
    
        self.inst.get().filterPeptidesByRT((deref(peptides.inst.get())), (<size_t>min_rt), (<size_t>max_rt))
    
    def filterPeptidesByMZ(self, PeptideIdentificationList peptides ,  min_mz ,  max_mz ):
        """
        filterPeptidesByMZ(self, peptides: PeptideIdentificationList , min_mz: int , max_mz: int ) -> None
        Filters peptide identifications by precursor m/z, keeping only IDs in the given range
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(min_mz, int) and min_mz >= 0, 'arg min_mz wrong type'
        assert isinstance(max_mz, int) and max_mz >= 0, 'arg max_mz wrong type'
    
    
    
        self.inst.get().filterPeptidesByMZ((deref(peptides.inst.get())), (<size_t>min_mz), (<size_t>max_mz))
    
    def filterPeptidesByMZError(self, PeptideIdentificationList peptides , double mass_error , bool unit_ppm ):
        """
        filterPeptidesByMZError(self, peptides: PeptideIdentificationList , mass_error: float , unit_ppm: bool ) -> None
        Filter peptide identifications according to mass deviation
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(mass_error, float), 'arg mass_error wrong type'
        assert isinstance(unit_ppm, pybool_t), 'arg unit_ppm wrong type'
    
    
    
        self.inst.get().filterPeptidesByMZError((deref(peptides.inst.get())), (<double>mass_error), (<bool>unit_ppm))
    
    def filterPeptidesByRTPredictPValue(self, PeptideIdentificationList peptides ,  metavalue_key , double threshold ):
        """
        filterPeptidesByRTPredictPValue(self, peptides: PeptideIdentificationList , metavalue_key: Union[bytes, str, String] , threshold: float ) -> None
        Filters peptide identifications according to p-values from RTPredict\n
        
        Filters the peptide hits by the probability (p-value) of a correct peptide identification having a deviation between observed and predicted RT equal to or greater than allowed
        
        
        :param peptides: Input/output
        :param metavalue_key: Name of the meta value that holds the p-value: "predicted_RT_p_value" or "predicted_RT_p_value_first_dim"
        :param threshold: P-value threshold
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert (isinstance(metavalue_key, str) or isinstance(metavalue_key, bytes) or isinstance(metavalue_key, String)), 'arg metavalue_key wrong type'
        assert isinstance(threshold, float), 'arg threshold wrong type'
    
    
    
        self.inst.get().filterPeptidesByRTPredictPValue((deref(peptides.inst.get())), deref((convString(metavalue_key)).get()), (<double>threshold))
    
    def removePeptidesWithMatchingModifications(self, PeptideIdentificationList peptides , set modifications ):
        """
        removePeptidesWithMatchingModifications(self, peptides: PeptideIdentificationList , modifications: Set[bytes] ) -> None
        Removes all peptide hits that have at least one of the given modifications
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(modifications, set) and all(isinstance(i, bytes) for i in modifications), 'arg modifications wrong type'
    
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in modifications:
           v1.insert(_String(<char *>item1))
        self.inst.get().removePeptidesWithMatchingModifications((deref(peptides.inst.get())), deref(v1))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        modifications.clear()
        modifications.update(replace)
        del v1
    
    def keepPeptidesWithMatchingModifications(self, PeptideIdentificationList peptides , set modifications ):
        """
        keepPeptidesWithMatchingModifications(self, peptides: PeptideIdentificationList , modifications: Set[bytes] ) -> None
        Keeps only peptide hits that have at least one of the given modifications
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(modifications, set) and all(isinstance(i, bytes) for i in modifications), 'arg modifications wrong type'
    
        cdef libcpp_set[_String] * v1 = new libcpp_set[_String]()
        cdef bytes item1
        for item1 in modifications:
           v1.insert(_String(<char *>item1))
        self.inst.get().keepPeptidesWithMatchingModifications((deref(peptides.inst.get())), deref(v1))
        cdef replace = set()
        cdef libcpp_set[_String].iterator it = v1.begin()
        while it != v1.end():
           replace.add(<char*>deref(it).c_str())
           inc(it)
        modifications.clear()
        modifications.update(replace)
        del v1
    
    def removePeptidesWithMatchingSequences(self, PeptideIdentificationList peptides , PeptideIdentificationList bad_peptides , bool ignore_mods ):
        """
        removePeptidesWithMatchingSequences(self, peptides: PeptideIdentificationList , bad_peptides: PeptideIdentificationList , ignore_mods: bool ) -> None
        Removes all peptide hits with a sequence that matches one in 'bad_peptides'
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(bad_peptides, PeptideIdentificationList), 'arg bad_peptides wrong type'
        assert isinstance(ignore_mods, pybool_t), 'arg ignore_mods wrong type'
    
    
    
        self.inst.get().removePeptidesWithMatchingSequences((deref(peptides.inst.get())), (deref(bad_peptides.inst.get())), (<bool>ignore_mods))
    
    def keepPeptidesWithMatchingSequences(self, PeptideIdentificationList peptides , PeptideIdentificationList bad_peptides , bool ignore_mods ):
        """
        keepPeptidesWithMatchingSequences(self, peptides: PeptideIdentificationList , bad_peptides: PeptideIdentificationList , ignore_mods: bool ) -> None
        Removes all peptide hits with a sequence that does not match one in 'good_peptides'
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(bad_peptides, PeptideIdentificationList), 'arg bad_peptides wrong type'
        assert isinstance(ignore_mods, pybool_t), 'arg ignore_mods wrong type'
    
    
    
        self.inst.get().keepPeptidesWithMatchingSequences((deref(peptides.inst.get())), (deref(bad_peptides.inst.get())), (<bool>ignore_mods))
    
    def keepUniquePeptidesPerProtein(self, PeptideIdentificationList peptides ):
        """
        keepUniquePeptidesPerProtein(self, peptides: PeptideIdentificationList ) -> None
        Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().keepUniquePeptidesPerProtein((deref(peptides.inst.get())))
    
    def removeDuplicatePeptideHits(self, PeptideIdentificationList peptides ):
        """
        removeDuplicatePeptideHits(self, peptides: PeptideIdentificationList ) -> None
        Removes duplicate peptide hits from each peptide identification, keeping only unique hits (per ID)
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
    
        self.inst.get().removeDuplicatePeptideHits((deref(peptides.inst.get())))
    
    def keepBestPerPeptide(self, PeptideIdentificationList peptides , bool ignore_mods , bool ignore_charges ,  nr_best_spectrum ):
        """
        keepBestPerPeptide(self, peptides: PeptideIdentificationList , ignore_mods: bool , ignore_charges: bool , nr_best_spectrum: int ) -> None
        Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence
        """
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(ignore_mods, pybool_t), 'arg ignore_mods wrong type'
        assert isinstance(ignore_charges, pybool_t), 'arg ignore_charges wrong type'
        assert isinstance(nr_best_spectrum, int) and nr_best_spectrum >= 0, 'arg nr_best_spectrum wrong type'
    
    
    
    
        self.inst.get().keepBestPerPeptide((deref(peptides.inst.get())), (<bool>ignore_mods), (<bool>ignore_charges), (<size_t>nr_best_spectrum))
    
    def keepBestPerPeptidePerRun(self, list prot_ids , PeptideIdentificationList peptides , bool ignore_mods , bool ignore_charges ,  nr_best_spectrum ):
        """
        keepBestPerPeptidePerRun(self, prot_ids: List[ProteinIdentification] , peptides: PeptideIdentificationList , ignore_mods: bool , ignore_charges: bool , nr_best_spectrum: int ) -> None
        Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence on a per run basis
        """
        assert isinstance(prot_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in prot_ids), 'arg prot_ids wrong type'
        assert isinstance(peptides, PeptideIdentificationList), 'arg peptides wrong type'
        assert isinstance(ignore_mods, pybool_t), 'arg ignore_mods wrong type'
        assert isinstance(ignore_charges, pybool_t), 'arg ignore_charges wrong type'
        assert isinstance(nr_best_spectrum, int) and nr_best_spectrum >= 0, 'arg nr_best_spectrum wrong type'
        cdef libcpp_vector[_ProteinIdentification] * v0 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item0
        for item0 in prot_ids:
            v0.push_back(deref(item0.inst.get()))
    
    
    
    
        self.inst.get().keepBestPerPeptidePerRun(deref(v0), (deref(peptides.inst.get())), (<bool>ignore_mods), (<bool>ignore_charges), (<size_t>nr_best_spectrum))
        cdef libcpp_vector[_ProteinIdentification].iterator it_prot_ids = v0.begin()
        replace_0 = []
        while it_prot_ids != v0.end():
            item0 = ProteinIdentification.__new__(ProteinIdentification)
            item0.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_prot_ids)))
            replace_0.append(item0)
            inc(it_prot_ids)
        prot_ids[:] = replace_0
        del v0
    DigestionFilter = __DigestionFilter 

cdef class Internal_MzMLValidator:
    """
    Cython implementation of _Internal_MzMLValidator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1Internal_MzMLValidator.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self, CVMappings mapping , ControlledVocabulary cv ):
        """
        __init__(self, mapping: CVMappings , cv: ControlledVocabulary ) -> None
        """
        assert isinstance(mapping, CVMappings), 'arg mapping wrong type'
        assert isinstance(cv, ControlledVocabulary), 'arg cv wrong type'
    
    
        self.inst = shared_ptr[_Internal_MzMLValidator](new _Internal_MzMLValidator((deref(mapping.inst.get())), (deref(cv.inst.get())))) 

cdef class LinearResampler:
    """
    Cython implementation of _LinearResampler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LinearResampler.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Annotates and filters transitions in a TargetedExperiment
    
    
    :param exp: The input, unfiltered transitions
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef LinearResampler rv = LinearResampler.__new__(LinearResampler)
       rv.inst = shared_ptr[_LinearResampler](new _LinearResampler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef LinearResampler rv = LinearResampler.__new__(LinearResampler)
       rv.inst = shared_ptr[_LinearResampler](new _LinearResampler(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_LinearResampler](new _LinearResampler())
    
    def _init_1(self, LinearResampler in_0 ):
        """
        _init_1(self, in_0: LinearResampler ) -> None
        """
        assert isinstance(in_0, LinearResampler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_LinearResampler](new _LinearResampler((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: LinearResampler ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], LinearResampler)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def raster(self, MSSpectrum input ):
        """
        raster(self, input: MSSpectrum ) -> None
        Applies the resampling algorithm to an MSSpectrum
        """
        assert isinstance(input, MSSpectrum), 'arg input wrong type'
    
        self.inst.get().raster((deref(input.inst.get())))
    
    def rasterExperiment(self, MSExperiment input ):
        """
        rasterExperiment(self, input: MSExperiment ) -> None
        Resamples the data in an MSExperiment
        """
        assert isinstance(input, MSExperiment), 'arg input wrong type'
    
        self.inst.get().rasterExperiment((deref(input.inst.get())))
    
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

cdef class MsInspectFile:
    """
    Cython implementation of _MsInspectFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MsInspectFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef MsInspectFile rv = MsInspectFile.__new__(MsInspectFile)
       rv.inst = shared_ptr[_MsInspectFile](new _MsInspectFile(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef MsInspectFile rv = MsInspectFile.__new__(MsInspectFile)
       rv.inst = shared_ptr[_MsInspectFile](new _MsInspectFile(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_MsInspectFile](new _MsInspectFile())
    
    def _init_1(self, MsInspectFile in_0 ):
        """
        _init_1(self, in_0: MsInspectFile ) -> None
        """
        assert isinstance(in_0, MsInspectFile), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_MsInspectFile](new _MsInspectFile((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: MsInspectFile ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], MsInspectFile)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def load(self,  filename , FeatureMap feature_map ):
        """
        load(self, filename: Union[bytes, str, String] , feature_map: FeatureMap ) -> None
        Loads a MsInspect file into a featureXML
        
        The content of the file is stored in `features`
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(feature_map, FeatureMap), 'arg feature_map wrong type'
    
    
        self.inst.get().load(deref((convString(filename)).get()), (deref(feature_map.inst.get())))
    
    def store(self,  filename , MSSpectrum spectrum ):
        """
        store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None
        Stores a featureXML as a MsInspect file
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(spectrum, MSSpectrum), 'arg spectrum wrong type'
    
    
        self.inst.get().store(deref((convString(filename)).get()), (deref(spectrum.inst.get()))) 

cdef class PepXMLFile:
    """
    Cython implementation of _PepXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PepXMLFile.html>`_
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __init__(self):
        """
        __init__(self) -> None
        """
        self.inst = shared_ptr[_PepXMLFile](new _PepXMLFile())
    
    def _load_0(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ):
        """
        _load_0(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _load_1(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  experiment_name ):
        """
        _load_1(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert (isinstance(experiment_name, str) or isinstance(experiment_name, bytes) or isinstance(experiment_name, String)), 'arg experiment_name wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((convString(experiment_name)).get()))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _load_2(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  experiment_name , SpectrumMetaDataLookup lookup ):
        """
        _load_2(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] , lookup: SpectrumMetaDataLookup ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert (isinstance(experiment_name, str) or isinstance(experiment_name, bytes) or isinstance(experiment_name, String)), 'arg experiment_name wrong type'
        assert isinstance(lookup, SpectrumMetaDataLookup), 'arg lookup wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
    
        self.inst.get().load(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((convString(experiment_name)).get()), (deref(lookup.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def load(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] , lookup: SpectrumMetaDataLookup ) -> None
          :noindex:
    
        """
        if (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)):
            return self._load_0(*args)
        elif (len(args)==4) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))):
            return self._load_1(*args)
        elif (len(args)==5) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and (isinstance(args[4], SpectrumMetaDataLookup)):
            return self._load_2(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def _store_0(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ):
        """
        _store_0(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def _store_1(self,  filename , list protein_ids , PeptideIdentificationList peptide_ids ,  mz_file ,  mz_name , bool peptideprophet_analyzed , double rt_tolerance ):
        """
        _store_1(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , mz_file: Union[bytes, str, String] , mz_name: Union[bytes, str, String] , peptideprophet_analyzed: bool , rt_tolerance: float ) -> None
        """
        assert (isinstance(filename, str) or isinstance(filename, bytes) or isinstance(filename, String)), 'arg filename wrong type'
        assert isinstance(protein_ids, list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in protein_ids), 'arg protein_ids wrong type'
        assert isinstance(peptide_ids, PeptideIdentificationList), 'arg peptide_ids wrong type'
        assert (isinstance(mz_file, str) or isinstance(mz_file, bytes) or isinstance(mz_file, String)), 'arg mz_file wrong type'
        assert (isinstance(mz_name, str) or isinstance(mz_name, bytes) or isinstance(mz_name, String)), 'arg mz_name wrong type'
        assert isinstance(peptideprophet_analyzed, pybool_t), 'arg peptideprophet_analyzed wrong type'
        assert isinstance(rt_tolerance, float), 'arg rt_tolerance wrong type'
    
        cdef libcpp_vector[_ProteinIdentification] * v1 = new libcpp_vector[_ProteinIdentification]()
        cdef ProteinIdentification item1
        for item1 in protein_ids:
            v1.push_back(deref(item1.inst.get()))
    
    
    
    
    
        self.inst.get().store(deref((convString(filename)).get()), deref(v1), (deref(peptide_ids.inst.get())), deref((convString(mz_file)).get()), deref((convString(mz_name)).get()), (<bool>peptideprophet_analyzed), (<double>rt_tolerance))
        cdef libcpp_vector[_ProteinIdentification].iterator it_protein_ids = v1.begin()
        replace_0 = []
        while it_protein_ids != v1.end():
            item1 = ProteinIdentification.__new__(ProteinIdentification)
            item1.inst = shared_ptr[_ProteinIdentification](new _ProteinIdentification(deref(it_protein_ids)))
            replace_0.append(item1)
            inc(it_protein_ids)
        protein_ids[:] = replace_0
        del v1
    
    def store(self, *args ):
        """
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , mz_file: Union[bytes, str, String] , mz_name: Union[bytes, str, String] , peptideprophet_analyzed: bool , rt_tolerance: float ) -> None
          :noindex:
    
        """
        if (len(args)==3) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)):
            return self._store_0(*args)
        elif (len(args)==7) and ((isinstance(args[0], str) or isinstance(args[0], bytes) or isinstance(args[0], String))) and (isinstance(args[1], list) and all(isinstance(elemt_rec, ProteinIdentification) for elemt_rec in args[1])) and (isinstance(args[2], PeptideIdentificationList)) and ((isinstance(args[3], str) or isinstance(args[3], bytes) or isinstance(args[3], String))) and ((isinstance(args[4], str) or isinstance(args[4], bytes) or isinstance(args[4], String))) and (isinstance(args[5], pybool_t)) and (isinstance(args[6], float)):
            return self._store_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def keepNativeSpectrumName(self, bool keep ):
        """
        keepNativeSpectrumName(self, keep: bool ) -> None
        """
        assert isinstance(keep, pybool_t), 'arg keep wrong type'
    
        self.inst.get().keepNativeSpectrumName((<bool>keep))
    
    def setParseUnknownScores(self, bool parse_unknown_scores ):
        """
        setParseUnknownScores(self, parse_unknown_scores: bool ) -> None
        """
        assert isinstance(parse_unknown_scores, pybool_t), 'arg parse_unknown_scores wrong type'
    
        self.inst.get().setParseUnknownScores((<bool>parse_unknown_scores)) 

cdef class SqrtScaler:
    """
    Cython implementation of _SqrtScaler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqrtScaler.html>`_
      -- Inherits from ['DefaultParamHandler']

    Scales the intensity of peaks to the sqrt
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef SqrtScaler rv = SqrtScaler.__new__(SqrtScaler)
       rv.inst = shared_ptr[_SqrtScaler](new _SqrtScaler(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef SqrtScaler rv = SqrtScaler.__new__(SqrtScaler)
       rv.inst = shared_ptr[_SqrtScaler](new _SqrtScaler(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_SqrtScaler](new _SqrtScaler())
    
    def _init_1(self, SqrtScaler in_0 ):
        """
        _init_1(self, in_0: SqrtScaler ) -> None
        """
        assert isinstance(in_0, SqrtScaler), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_SqrtScaler](new _SqrtScaler((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: SqrtScaler ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], SqrtScaler)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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

cdef class ThresholdMower:
    """
    Cython implementation of _ThresholdMower

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ThresholdMower.html>`_
      -- Inherits from ['DefaultParamHandler']
    """

    # see .pxd file for cdef of inst ptr

    def __dealloc__(self):
         self.inst.reset()

    
    def __copy__(self):
       cdef ThresholdMower rv = ThresholdMower.__new__(ThresholdMower)
       rv.inst = shared_ptr[_ThresholdMower](new _ThresholdMower(deref(self.inst.get())))
       return rv
    
    def __deepcopy__(self, memo):
       cdef ThresholdMower rv = ThresholdMower.__new__(ThresholdMower)
       rv.inst = shared_ptr[_ThresholdMower](new _ThresholdMower(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        """
        _init_0(self) -> None
        """
        self.inst = shared_ptr[_ThresholdMower](new _ThresholdMower())
    
    def _init_1(self, ThresholdMower in_0 ):
        """
        _init_1(self, in_0: ThresholdMower ) -> None
        """
        assert isinstance(in_0, ThresholdMower), 'arg in_0 wrong type'
    
        self.inst = shared_ptr[_ThresholdMower](new _ThresholdMower((deref(in_0.inst.get()))))
    
    def __init__(self, *args , **kwargs):
        """
        .. rubric:: Overload:
        .. py:function:: __init__(self, ) -> None
          :noindex:
        
        .. rubric:: Overload:
        .. py:function:: __init__(self, in_0: ThresholdMower ) -> None
          :noindex:
    
        """
        if not args:
             self._init_0(*args)
        elif (len(args)==1) and (isinstance(args[0], ThresholdMower)):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
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
