from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AMSE_AdductInfo:
    """
    Cython implementation of _AMSE_AdductInfo

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AMSE_AdductInfo.html>`_
    """
    
    def __init__(self, name: Union[bytes, str, String] , adduct: EmpiricalFormula , charge: int , mol_multiplier: int ) -> None:
        """
        Cython signature: void AMSE_AdductInfo(const String & name, EmpiricalFormula & adduct, int charge, unsigned int mol_multiplier)
        """
        ...
    
    def getNeutralMass(self, observed_mz: float ) -> float:
        """
        Cython signature: double getNeutralMass(double observed_mz)
        """
        ...
    
    def getMZ(self, neutral_mass: float ) -> float:
        """
        Cython signature: double getMZ(double neutral_mass)
        """
        ...
    
    def isCompatible(self, db_entry: EmpiricalFormula ) -> bool:
        """
        Cython signature: bool isCompatible(EmpiricalFormula db_entry)
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ... 


class AbsoluteQuantitation:
    """
    Cython implementation of _AbsoluteQuantitation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AbsoluteQuantitation.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AbsoluteQuantitation()
        """
        ...
    
    @overload
    def __init__(self, in_0: AbsoluteQuantitation ) -> None:
        """
        Cython signature: void AbsoluteQuantitation(AbsoluteQuantitation &)
        """
        ...
    
    def setQuantMethods(self, quant_methods: List[AbsoluteQuantitationMethod] ) -> None:
        """
        Cython signature: void setQuantMethods(libcpp_vector[AbsoluteQuantitationMethod] & quant_methods)
        """
        ...
    
    def getQuantMethods(self) -> List[AbsoluteQuantitationMethod]:
        """
        Cython signature: libcpp_vector[AbsoluteQuantitationMethod] getQuantMethods()
        """
        ...
    
    def calculateRatio(self, component_1: Feature , component_2: Feature , feature_name: Union[bytes, str, String] ) -> float:
        """
        Cython signature: double calculateRatio(Feature & component_1, Feature & component_2, const String & feature_name)
        """
        ...
    
    def applyCalibration(self, component: Feature , IS_component: Feature , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param ) -> float:
        """
        Cython signature: double applyCalibration(const Feature & component, const Feature & IS_component, const String & feature_name, const String & transformation_model, const Param & transformation_model_params)
        """
        ...
    
    def quantifyComponents(self, unknowns: FeatureMap ) -> None:
        """
        Cython signature: void quantifyComponents(FeatureMap & unknowns)
        This function applies the calibration curve, hence quantifying all the components
        """
        ...
    
    def optimizeCalibrationCurveIterative(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param , optimized_params: Param ) -> bool:
        """
        Cython signature: bool optimizeCalibrationCurveIterative(libcpp_vector[AQS_featureConcentration] & component_concentrations, const String & feature_name, const String & transformation_model, const Param & transformation_model_params, Param & optimized_params)
        """
        ...
    
    def optimizeSingleCalibrationCurve(self, component_name: Union[bytes, str, String] , component_concentrations: List[AQS_featureConcentration] ) -> None:
        """
        Cython signature: void optimizeSingleCalibrationCurve(const String & component_name, libcpp_vector[AQS_featureConcentration] & component_concentrations)
        """
        ...
    
    def calculateBias(self, actual_concentration: float , calculated_concentration: float ) -> float:
        """
        Cython signature: double calculateBias(double actual_concentration, double calculated_concentration)
        This function calculates the bias of the calibration
        """
        ...
    
    def fitCalibration(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param ) -> Param:
        """
        Cython signature: Param fitCalibration(libcpp_vector[AQS_featureConcentration] & component_concentrations, const String & feature_name, const String & transformation_model, Param transformation_model_params)
        """
        ...
    
    def calculateBiasAndR(self, component_concentrations: List[AQS_featureConcentration] , feature_name: Union[bytes, str, String] , transformation_model: Union[bytes, str, String] , transformation_model_params: Param , biases: List[float] , correlation_coefficient: float ) -> None:
        """
        Cython signature: void calculateBiasAndR(libcpp_vector[AQS_featureConcentration] & component_concentrations, const String & feature_name, const String & transformation_model, Param & transformation_model_params, libcpp_vector[double] & biases, double & correlation_coefficient)
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class ChromatogramSettings:
    """
    Cython implementation of _ChromatogramSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramSettings.html>`_
      -- Inherits from ['MetaInfoInterface']

    Description of the chromatogram settings, provides meta-information
    about a single chromatogram.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramSettings()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramSettings ) -> None:
        """
        Cython signature: void ChromatogramSettings(ChromatogramSettings &)
        """
        ...
    
    def getProduct(self) -> Product:
        """
        Cython signature: Product getProduct()
        Returns the product ion
        """
        ...
    
    def setProduct(self, p: Product ) -> None:
        """
        Cython signature: void setProduct(Product p)
        Sets the product ion
        """
        ...
    
    def getNativeID(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNativeID()
        Returns the native identifier for the spectrum, used by the acquisition software.
        """
        ...
    
    def setNativeID(self, native_id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNativeID(String native_id)
        Sets the native identifier for the spectrum, used by the acquisition software.
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the free-text comment
        """
        ...
    
    def setComment(self, comment: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String comment)
        Sets the free-text comment
        """
        ...
    
    def getInstrumentSettings(self) -> InstrumentSettings:
        """
        Cython signature: InstrumentSettings getInstrumentSettings()
        Returns the instrument settings of the current spectrum
        """
        ...
    
    def setInstrumentSettings(self, instrument_settings: InstrumentSettings ) -> None:
        """
        Cython signature: void setInstrumentSettings(InstrumentSettings instrument_settings)
        Sets the instrument settings of the current spectrum
        """
        ...
    
    def getAcquisitionInfo(self) -> AcquisitionInfo:
        """
        Cython signature: AcquisitionInfo getAcquisitionInfo()
        Returns the acquisition info
        """
        ...
    
    def setAcquisitionInfo(self, acquisition_info: AcquisitionInfo ) -> None:
        """
        Cython signature: void setAcquisitionInfo(AcquisitionInfo acquisition_info)
        Sets the acquisition info
        """
        ...
    
    def getSourceFile(self) -> SourceFile:
        """
        Cython signature: SourceFile getSourceFile()
        Returns the source file
        """
        ...
    
    def setSourceFile(self, source_file: SourceFile ) -> None:
        """
        Cython signature: void setSourceFile(SourceFile source_file)
        Sets the source file
        """
        ...
    
    def getPrecursor(self) -> Precursor:
        """
        Cython signature: Precursor getPrecursor()
        Returns the precursors
        """
        ...
    
    def setPrecursor(self, precursor: Precursor ) -> None:
        """
        Cython signature: void setPrecursor(Precursor precursor)
        Sets the precursors
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[shared_ptr[DataProcessing]] getDataProcessing()
        Returns the description of the applied processing
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[shared_ptr[DataProcessing]])
        Sets the description of the applied processing
        """
        ...
    
    def setChromatogramType(self, type: int ) -> None:
        """
        Cython signature: void setChromatogramType(ChromatogramType type)
        Sets the chromatogram type
        """
        ...
    
    def getChromatogramType(self) -> int:
        """
        Cython signature: ChromatogramType getChromatogramType()
        Get the chromatogram type
        """
        ...
    
    def isMetaEmpty(self) -> bool:
        """
        Cython signature: bool isMetaEmpty()
        Returns if the MetaInfo is empty
        """
        ...
    
    def clearMetaInfo(self) -> None:
        """
        Cython signature: void clearMetaInfo()
        Removes all meta values
        """
        ...
    
    def metaRegistry(self) -> MetaInfoRegistry:
        """
        Cython signature: MetaInfoRegistry metaRegistry()
        Returns a reference to the MetaInfoRegistry
        """
        ...
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        Fills the given vector with a list of all keys for which a value is set
        """
        ...
    
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        Returns the value corresponding to a string, or
        """
        ...
    
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        Sets the DataValue corresponding to a name
        """
        ...
    
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        Returns whether an entry with the given name exists
        """
        ...
    
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        Removes the DataValue corresponding to `name` if it exists
        """
        ...
    
    def __richcmp__(self, other: ChromatogramSettings, op: int) -> Any:
        ...
    ChromatogramType : __ChromatogramType 


class CubicSpline2d:
    """
    Cython implementation of _CubicSpline2d

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CubicSpline2d.html>`_
    """
    
    @overload
    def __init__(self, x: List[float] , y: List[float] ) -> None:
        """
        Cython signature: void CubicSpline2d(libcpp_vector[double] x, libcpp_vector[double] y)
        """
        ...
    
    @overload
    def __init__(self, in_0: CubicSpline2d ) -> None:
        """
        Cython signature: void CubicSpline2d(CubicSpline2d &)
        """
        ...
    
    @overload
    def __init__(self, m: Dict[float, float] ) -> None:
        """
        Cython signature: void CubicSpline2d(libcpp_map[double,double] m)
        """
        ...
    
    def eval(self, x: float ) -> float:
        """
        Cython signature: double eval(double x)
        Evaluates the cubic spline
        """
        ...
    
    def derivatives(self, x: float , order: int ) -> float:
        """
        Cython signature: double derivatives(double x, unsigned int order)
        Returns first, second or third derivative of cubic spline
        """
        ... 


class __DigestionFilter:
    """
    Cython implementation of _DigestionFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DigestionFilter.html>`_
    """
    
    digestion_: ProteaseDigestion
    
    ignore_missed_cleavages_: bool
    
    methionine_cleavage_: bool
    
    def __init__(self, entries: List[FASTAEntry] , digestion: ProteaseDigestion , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> None:
        """
        Cython signature: void DigestionFilter(libcpp_vector[FASTAEntry] & entries, ProteaseDigestion & digestion, bool ignore_missed_cleavages, bool methionine_cleavage)
        """
        ...
    
    def filterPeptideEvidences(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void filterPeptideEvidences(PeptideIdentificationList & peptides)
        """
        ... 


class DistanceMatrix:
    """
    Cython implementation of _DistanceMatrix[float]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DistanceMatrix[float].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DistanceMatrix()
        """
        ...
    
    @overload
    def __init__(self, in_0: DistanceMatrix ) -> None:
        """
        Cython signature: void DistanceMatrix(DistanceMatrix &)
        """
        ...
    
    @overload
    def __init__(self, dimensionsize: int , value: float ) -> None:
        """
        Cython signature: void DistanceMatrix(size_t dimensionsize, float value)
        """
        ...
    
    def getValue(self, i: int , j: int ) -> float:
        """
        Cython signature: float getValue(size_t i, size_t j)
        """
        ...
    
    def setValue(self, i: int , j: int , value: float ) -> None:
        """
        Cython signature: void setValue(size_t i, size_t j, float value)
        """
        ...
    
    def setValueQuick(self, i: int , j: int , value: float ) -> None:
        """
        Cython signature: void setValueQuick(size_t i, size_t j, float value)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def resize(self, dimensionsize: int , value: float ) -> None:
        """
        Cython signature: void resize(size_t dimensionsize, float value)
        """
        ...
    
    def reduce(self, j: int ) -> None:
        """
        Cython signature: void reduce(size_t j)
        """
        ...
    
    def dimensionsize(self) -> int:
        """
        Cython signature: size_t dimensionsize()
        """
        ...
    
    def updateMinElement(self) -> None:
        """
        Cython signature: void updateMinElement()
        """
        ...
    
    def getMinElementCoordinates(self) -> List[int, int]:
        """
        Cython signature: libcpp_pair[size_t,size_t] getMinElementCoordinates()
        """
        ...
    
    def __richcmp__(self, other: DistanceMatrix, op: int) -> Any:
        ... 


class ExperimentalSettings:
    """
    Cython implementation of _ExperimentalSettings

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalSettings.html>`_
      -- Inherits from ['DocumentIdentifier', 'MetaInfoInterface']

    Description of the experimental settings, provides meta-information
    about an LC-MS/MS injection.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ExperimentalSettings()
        """
        ...
    
    @overload
    def __init__(self, in_0: ExperimentalSettings ) -> None:
        """
        Cython signature: void ExperimentalSettings(ExperimentalSettings &)
        """
        ...
    
    def getSourceFiles(self) -> List[SourceFile]:
        """
        Cython signature: libcpp_vector[SourceFile] getSourceFiles()
        Returns a reference to the source data file
        """
        ...
    
    def setSourceFiles(self, source_files: List[SourceFile] ) -> None:
        """
        Cython signature: void setSourceFiles(libcpp_vector[SourceFile] source_files)
        Sets the source data file
        """
        ...
    
    def getDateTime(self) -> DateTime:
        """
        Cython signature: DateTime getDateTime()
        Returns the date the experiment was performed
        """
        ...
    
    def setDateTime(self, date_time: DateTime ) -> None:
        """
        Cython signature: void setDateTime(DateTime date_time)
        Sets the date the experiment was performed
        """
        ...
    
    def getSample(self) -> Sample:
        """
        Cython signature: Sample getSample()
        Returns a reference to the sample description
        """
        ...
    
    def setSample(self, sample: Sample ) -> None:
        """
        Cython signature: void setSample(Sample sample)
        Sets the sample description
        """
        ...
    
    def getContacts(self) -> List[ContactPerson]:
        """
        Cython signature: libcpp_vector[ContactPerson] getContacts()
        Returns a reference to the list of contact persons
        """
        ...
    
    def setContacts(self, contacts: List[ContactPerson] ) -> None:
        """
        Cython signature: void setContacts(libcpp_vector[ContactPerson] contacts)
        Sets the list of contact persons
        """
        ...
    
    def getInstrument(self) -> Instrument:
        """
        Cython signature: Instrument getInstrument()
        Returns a reference to the MS instrument description
        """
        ...
    
    def setInstrument(self, instrument: Instrument ) -> None:
        """
        Cython signature: void setInstrument(Instrument instrument)
        Sets the MS instrument description
        """
        ...
    
    def getHPLC(self) -> HPLC:
        """
        Cython signature: HPLC getHPLC()
        Returns a reference to the description of the HPLC run
        """
        ...
    
    def setHPLC(self, hplc: HPLC ) -> None:
        """
        Cython signature: void setHPLC(HPLC hplc)
        Sets the description of the HPLC run
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the free-text comment
        """
        ...
    
    def setComment(self, comment: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String comment)
        Sets the free-text comment
        """
        ...
    
    def getFractionIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getFractionIdentifier()
        Returns fraction identifier
        """
        ...
    
    def setFractionIdentifier(self, fraction_identifier: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setFractionIdentifier(String fraction_identifier)
        Sets the fraction identifier
        """
        ...
    
    def setIdentifier(self, id: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(String id)
        Sets document identifier (e.g. an LSID)
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        Retrieve document identifier (e.g. an LSID)
        """
        ...
    
    def setLoadedFileType(self, file_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLoadedFileType(String file_name)
        Sets the file_type according to the type of the file loaded from, preferably done whilst loading
        """
        ...
    
    def getLoadedFileType(self) -> int:
        """
        Cython signature: int getLoadedFileType()
        Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded
        """
        ...
    
    def setLoadedFilePath(self, file_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setLoadedFilePath(String file_name)
        Sets the file_name according to absolute path of the file loaded, preferably done whilst loading
        """
        ...
    
    def getLoadedFilePath(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getLoadedFilePath()
        Returns the file_name which is the absolute path to the file loaded
        """
        ...
    
    def isMetaEmpty(self) -> bool:
        """
        Cython signature: bool isMetaEmpty()
        Returns if the MetaInfo is empty
        """
        ...
    
    def clearMetaInfo(self) -> None:
        """
        Cython signature: void clearMetaInfo()
        Removes all meta values
        """
        ...
    
    def metaRegistry(self) -> MetaInfoRegistry:
        """
        Cython signature: MetaInfoRegistry metaRegistry()
        Returns a reference to the MetaInfoRegistry
        """
        ...
    
    def getKeys(self, keys: List[bytes] ) -> None:
        """
        Cython signature: void getKeys(libcpp_vector[String] & keys)
        Fills the given vector with a list of all keys for which a value is set
        """
        ...
    
    def getMetaValue(self, in_0: Union[bytes, str, String] ) -> Union[int, float, bytes, str, List[int], List[float], List[bytes]]:
        """
        Cython signature: DataValue getMetaValue(String)
        Returns the value corresponding to a string, or
        """
        ...
    
    def setMetaValue(self, in_0: Union[bytes, str, String] , in_1: Union[int, float, bytes, str, List[int], List[float], List[bytes]] ) -> None:
        """
        Cython signature: void setMetaValue(String, DataValue)
        Sets the DataValue corresponding to a name
        """
        ...
    
    def metaValueExists(self, in_0: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool metaValueExists(String)
        Returns whether an entry with the given name exists
        """
        ...
    
    def removeMetaValue(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void removeMetaValue(String)
        Removes the DataValue corresponding to `name` if it exists
        """
        ...
    
    def __richcmp__(self, other: ExperimentalSettings, op: int) -> Any:
        ... 


class HPLC:
    """
    Cython implementation of _HPLC

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1HPLC.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void HPLC()
        Representation of a HPLC experiment
        """
        ...
    
    @overload
    def __init__(self, in_0: HPLC ) -> None:
        """
        Cython signature: void HPLC(HPLC &)
        """
        ...
    
    def getInstrument(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getInstrument()
        Returns a reference to the instument name
        """
        ...
    
    def setInstrument(self, instrument: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setInstrument(String instrument)
        Sets the instument name
        """
        ...
    
    def getColumn(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getColumn()
        Returns a reference to the column description
        """
        ...
    
    def setColumn(self, column: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setColumn(String column)
        Sets the column description
        """
        ...
    
    def getTemperature(self) -> int:
        """
        Cython signature: int getTemperature()
        Returns the temperature (in degree C)
        """
        ...
    
    def setTemperature(self, temperature: int ) -> None:
        """
        Cython signature: void setTemperature(int temperature)
        Sets the temperature (in degree C)
        """
        ...
    
    def getPressure(self) -> int:
        """
        Cython signature: unsigned int getPressure()
        Returns the pressure (in bar)
        """
        ...
    
    def setPressure(self, pressure: int ) -> None:
        """
        Cython signature: void setPressure(unsigned int pressure)
        Sets the pressure (in bar)
        """
        ...
    
    def getFlux(self) -> int:
        """
        Cython signature: unsigned int getFlux()
        Returns the flux (in microliter/sec)
        """
        ...
    
    def setFlux(self, flux: int ) -> None:
        """
        Cython signature: void setFlux(unsigned int flux)
        Sets the flux (in microliter/sec)
        """
        ...
    
    def getComment(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getComment()
        Returns the comments
        """
        ...
    
    def setComment(self, comment: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setComment(String comment)
        Sets the comments
        """
        ...
    
    def getGradient(self) -> Gradient:
        """
        Cython signature: Gradient getGradient()
        Returns a mutable reference to the used gradient
        """
        ...
    
    def setGradient(self, gradient: Gradient ) -> None:
        """
        Cython signature: void setGradient(Gradient gradient)
        Sets the used gradient
        """
        ... 


class IDFilter:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IDFilter()
        """
        ...
    
    @overload
    def __init__(self, in_0: IDFilter ) -> None:
        """
        Cython signature: void IDFilter(IDFilter &)
        """
        ...
    
    @overload
    def countHits(self, identifications: PeptideIdentificationList ) -> int:
        """
        Cython signature: size_t countHits(PeptideIdentificationList identifications)
        Returns the total number of peptide hits in a vector of peptide identifications
        """
        ...
    
    @overload
    def countHits(self, identifications: List[ProteinIdentification] ) -> int:
        """
        Cython signature: size_t countHits(libcpp_vector[ProteinIdentification] identifications)
        Returns the total number of protein hits in a vector of protein identifications
        """
        ...
    
    @overload
    def getBestHit(self, identifications: PeptideIdentificationList , assume_sorted: bool , best_hit: PeptideHit ) -> bool:
        """
        Cython signature: bool getBestHit(PeptideIdentificationList identifications, bool assume_sorted, PeptideHit & best_hit)
        Finds the best-scoring hit in a vector of peptide or protein identifications\n
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of peptide identifications
        :return: true if a hit was present, false otherwise
        """
        ...
    
    @overload
    def getBestHit(self, identifications: List[ProteinIdentification] , assume_sorted: bool , best_hit: ProteinHit ) -> bool:
        """
        Cython signature: bool getBestHit(libcpp_vector[ProteinIdentification] identifications, bool assume_sorted, ProteinHit & best_hit)
        Finds the best-scoring hit in a vector of peptide or protein identifications
        
        If there are several hits with the best score, the first one is taken
        
        
        :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
        :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
        :param best_hit: Contains the best hit if successful in a vector of protein identifications
        :return: true if a hit was present, false otherwise
        """
        ...
    
    def extractPeptideSequences(self, peptides: PeptideIdentificationList , sequences: Set[bytes] , ignore_mods: bool ) -> None:
        """
        Cython signature: void extractPeptideSequences(PeptideIdentificationList & peptides, libcpp_set[String] & sequences, bool ignore_mods)
        Extracts all unique peptide sequences from a list of peptide IDs
        
        
        :param peptides:
        :param ignore_mods: Boolean operator default to false in case of any modifications in sequences during extraction
        :return: Sequences
        """
        ...
    
    def removeUnreferencedProteins(self, proteins: List[ProteinIdentification] , peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void removeUnreferencedProteins(libcpp_vector[ProteinIdentification] & proteins, PeptideIdentificationList & peptides)
        Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list
        """
        ...
    
    def updateProteinReferences(self, peptides: PeptideIdentificationList , proteins: List[ProteinIdentification] , remove_peptides_without_reference: bool ) -> None:
        """
        Cython signature: void updateProteinReferences(PeptideIdentificationList & peptides, libcpp_vector[ProteinIdentification] & proteins, bool remove_peptides_without_reference)
        Removes references to missing proteins. Only PeptideEvidence entries that reference protein hits in 'proteins' are kept in the peptide hits
        """
        ...
    
    def updateProteinGroups(self, groups: List[ProteinGroup] , hits: List[ProteinHit] ) -> bool:
        """
        Cython signature: bool updateProteinGroups(libcpp_vector[ProteinGroup] & groups, libcpp_vector[ProteinHit] & hits)
        Update protein groups after protein hits were filtered
        
        
        :param groups: Input/output protein groups
        :param hits: Available protein hits (all others are removed from the groups)
        :return: Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed)
        """
        ...
    
    @overload
    def removeEmptyIdentifications(self, ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void removeEmptyIdentifications(PeptideIdentificationList & ids)
        Removes peptide or protein identifications that have no hits in them
        """
        ...
    
    @overload
    def removeEmptyIdentifications(self, ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void removeEmptyIdentifications(libcpp_vector[ProteinIdentification] & ids)
        Removes peptide or protein identifications that have no hits in them
        """
        ...
    
    @overload
    def filterHitsByScore(self, ids: PeptideIdentificationList , threshold_score: float ) -> None:
        """
        Cython signature: void filterHitsByScore(PeptideIdentificationList & ids, double threshold_score)
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept
        """
        ...
    
    @overload
    def filterHitsByScore(self, ids: List[ProteinIdentification] , threshold_score: float ) -> None:
        """
        Cython signature: void filterHitsByScore(libcpp_vector[ProteinIdentification] & ids, double threshold_score)
        Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification/ProteinIdentifiation. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept
        """
        ...
    
    @overload
    def filterHitsByScore(self, experiment: AnnotatedMSRun , peptide_threshold_score: float , protein_threshold_score: float ) -> None:
        """
        Cython signature: void filterHitsByScore(AnnotatedMSRun & experiment, double peptide_threshold_score, double protein_threshold_score)
        Filters an MS/MS experiment according to score thresholds
        """
        ...
    
    def keepNBestSpectra(self, peptides: PeptideIdentificationList , n: int ) -> None:
        """
        Cython signature: void keepNBestSpectra(PeptideIdentificationList & peptides, size_t n)
        Filter identifications by "N best" PeptideIdentification objects (better PeptideIdentification means better [best] PeptideHit than other)
        """
        ...
    
    @overload
    def keepNBestHits(self, ids: PeptideIdentificationList , n: int ) -> None:
        """
        Cython signature: void keepNBestHits(PeptideIdentificationList & ids, size_t n)
        """
        ...
    
    @overload
    def keepNBestHits(self, ids: List[ProteinIdentification] , n: int ) -> None:
        """
        Cython signature: void keepNBestHits(libcpp_vector[ProteinIdentification] & ids, size_t n)
        """
        ...
    
    @overload
    def keepNBestHits(self, experiment: AnnotatedMSRun , n: int ) -> None:
        """
        Cython signature: void keepNBestHits(AnnotatedMSRun & experiment, size_t n)
        Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
        """
        ...
    
    @overload
    def filterHitsByRank(self, ids: PeptideIdentificationList , min_rank: int , max_rank: int ) -> None:
        """
        Cython signature: void filterHitsByRank(PeptideIdentificationList & ids, size_t min_rank, size_t max_rank)
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
        """
        ...
    
    @overload
    def filterHitsByRank(self, ids: List[ProteinIdentification] , min_rank: int , max_rank: int ) -> None:
        """
        Cython signature: void filterHitsByRank(libcpp_vector[ProteinIdentification] & ids, size_t min_rank, size_t max_rank)
        Filters peptide or protein identifications according to the ranking of the hits\n
        
        The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
        Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
        The ranks are (re-)computed before filtering
        'max_rank' is ignored if it is smaller than 'min_rank'
        
        
        Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis
        """
        ...
    
    @overload
    def removeDecoyHits(self, ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void removeDecoyHits(PeptideIdentificationList & ids)
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
        """
        ...
    
    @overload
    def removeDecoyHits(self, ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void removeDecoyHits(libcpp_vector[ProteinIdentification] & ids)
        Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
        """
        ...
    
    @overload
    def removeHitsMatchingProteins(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None:
        """
        Cython signature: void removeHitsMatchingProteins(PeptideIdentificationList & ids, libcpp_set[String] accessions)
        Filters peptide or protein identifications according to the given proteins (negative)
        """
        ...
    
    @overload
    def removeHitsMatchingProteins(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None:
        """
        Cython signature: void removeHitsMatchingProteins(libcpp_vector[ProteinIdentification] & ids, libcpp_set[String] accessions)
        Filters peptide or protein identifications according to the given proteins (negative)
        """
        ...
    
    @overload
    def keepHitsMatchingProteins(self, ids: PeptideIdentificationList , accessions: Set[bytes] ) -> None:
        """
        Cython signature: void keepHitsMatchingProteins(PeptideIdentificationList & ids, libcpp_set[String] accessions)
        Filters peptide or protein identifications according to the given proteins (positive)
        """
        ...
    
    @overload
    def keepHitsMatchingProteins(self, ids: List[ProteinIdentification] , accessions: Set[bytes] ) -> None:
        """
        Cython signature: void keepHitsMatchingProteins(libcpp_vector[ProteinIdentification] & ids, libcpp_set[String] accessions)
        Filters peptide or protein identifications according to the given proteins (positive)
        """
        ...
    
    @overload
    def keepHitsMatchingProteins(self, experiment: AnnotatedMSRun , proteins: List[FASTAEntry] ) -> None:
        """
        Cython signature: void keepHitsMatchingProteins(AnnotatedMSRun & experiment, libcpp_vector[FASTAEntry] & proteins)
        """
        ...
    
    def keepBestPeptideHits(self, peptides: PeptideIdentificationList , strict: bool ) -> None:
        """
        Cython signature: void keepBestPeptideHits(PeptideIdentificationList & peptides, bool strict)
        Filters peptide identifications keeping only the single best-scoring hit per ID
        
        
        :param peptides: Input/output
        :param strict: If set, keep the best hit only if its score is unique - i.e. ties are not allowed. (Otherwise all hits with the best score is kept.)
        """
        ...
    
    def filterPeptidesByLength(self, peptides: PeptideIdentificationList , min_length: int , max_length: int ) -> None:
        """
        Cython signature: void filterPeptidesByLength(PeptideIdentificationList & peptides, size_t min_length, size_t max_length)
        Filters peptide identifications according to peptide sequence length
        """
        ...
    
    def filterPeptidesByCharge(self, peptides: PeptideIdentificationList , min_charge: int , max_charge: int ) -> None:
        """
        Cython signature: void filterPeptidesByCharge(PeptideIdentificationList & peptides, size_t min_charge, size_t max_charge)
        Filters peptide identifications according to charge state
        """
        ...
    
    def filterPeptidesByRT(self, peptides: PeptideIdentificationList , min_rt: int , max_rt: int ) -> None:
        """
        Cython signature: void filterPeptidesByRT(PeptideIdentificationList & peptides, size_t min_rt, size_t max_rt)
        Filters peptide identifications by precursor RT, keeping only IDs in the given range
        """
        ...
    
    def filterPeptidesByMZ(self, peptides: PeptideIdentificationList , min_mz: int , max_mz: int ) -> None:
        """
        Cython signature: void filterPeptidesByMZ(PeptideIdentificationList & peptides, size_t min_mz, size_t max_mz)
        Filters peptide identifications by precursor m/z, keeping only IDs in the given range
        """
        ...
    
    def filterPeptidesByMZError(self, peptides: PeptideIdentificationList , mass_error: float , unit_ppm: bool ) -> None:
        """
        Cython signature: void filterPeptidesByMZError(PeptideIdentificationList & peptides, double mass_error, bool unit_ppm)
        Filter peptide identifications according to mass deviation
        """
        ...
    
    def filterPeptidesByRTPredictPValue(self, peptides: PeptideIdentificationList , metavalue_key: Union[bytes, str, String] , threshold: float ) -> None:
        """
        Cython signature: void filterPeptidesByRTPredictPValue(PeptideIdentificationList & peptides, const String & metavalue_key, double threshold)
        Filters peptide identifications according to p-values from RTPredict\n
        
        Filters the peptide hits by the probability (p-value) of a correct peptide identification having a deviation between observed and predicted RT equal to or greater than allowed
        
        
        :param peptides: Input/output
        :param metavalue_key: Name of the meta value that holds the p-value: "predicted_RT_p_value" or "predicted_RT_p_value_first_dim"
        :param threshold: P-value threshold
        """
        ...
    
    def removePeptidesWithMatchingModifications(self, peptides: PeptideIdentificationList , modifications: Set[bytes] ) -> None:
        """
        Cython signature: void removePeptidesWithMatchingModifications(PeptideIdentificationList & peptides, libcpp_set[String] & modifications)
        Removes all peptide hits that have at least one of the given modifications
        """
        ...
    
    def keepPeptidesWithMatchingModifications(self, peptides: PeptideIdentificationList , modifications: Set[bytes] ) -> None:
        """
        Cython signature: void keepPeptidesWithMatchingModifications(PeptideIdentificationList & peptides, libcpp_set[String] & modifications)
        Keeps only peptide hits that have at least one of the given modifications
        """
        ...
    
    def removePeptidesWithMatchingSequences(self, peptides: PeptideIdentificationList , bad_peptides: PeptideIdentificationList , ignore_mods: bool ) -> None:
        """
        Cython signature: void removePeptidesWithMatchingSequences(PeptideIdentificationList & peptides, PeptideIdentificationList & bad_peptides, bool ignore_mods)
        Removes all peptide hits with a sequence that matches one in 'bad_peptides'
        """
        ...
    
    def keepPeptidesWithMatchingSequences(self, peptides: PeptideIdentificationList , bad_peptides: PeptideIdentificationList , ignore_mods: bool ) -> None:
        """
        Cython signature: void keepPeptidesWithMatchingSequences(PeptideIdentificationList & peptides, PeptideIdentificationList & bad_peptides, bool ignore_mods)
        Removes all peptide hits with a sequence that does not match one in 'good_peptides'
        """
        ...
    
    def keepUniquePeptidesPerProtein(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void keepUniquePeptidesPerProtein(PeptideIdentificationList & peptides)
        Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)
        """
        ...
    
    def removeDuplicatePeptideHits(self, peptides: PeptideIdentificationList ) -> None:
        """
        Cython signature: void removeDuplicatePeptideHits(PeptideIdentificationList & peptides)
        Removes duplicate peptide hits from each peptide identification, keeping only unique hits (per ID)
        """
        ...
    
    def keepBestPerPeptide(self, peptides: PeptideIdentificationList , ignore_mods: bool , ignore_charges: bool , nr_best_spectrum: int ) -> None:
        """
        Cython signature: void keepBestPerPeptide(PeptideIdentificationList & peptides, bool ignore_mods, bool ignore_charges, size_t nr_best_spectrum)
        Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence
        """
        ...
    
    def keepBestPerPeptidePerRun(self, prot_ids: List[ProteinIdentification] , peptides: PeptideIdentificationList , ignore_mods: bool , ignore_charges: bool , nr_best_spectrum: int ) -> None:
        """
        Cython signature: void keepBestPerPeptidePerRun(libcpp_vector[ProteinIdentification] & prot_ids, PeptideIdentificationList & peptides, bool ignore_mods, bool ignore_charges, size_t nr_best_spectrum)
        Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence on a per run basis
        """
        ... 


class Internal_MzMLValidator:
    """
    Cython implementation of _Internal_MzMLValidator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1Internal_MzMLValidator.html>`_
    """
    
    def __init__(self, mapping: CVMappings , cv: ControlledVocabulary ) -> None:
        """
        Cython signature: void Internal_MzMLValidator(CVMappings & mapping, ControlledVocabulary & cv)
        """
        ... 


class LinearResampler:
    """
    Cython implementation of _LinearResampler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LinearResampler.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']

    Annotates and filters transitions in a TargetedExperiment
    
    
    :param exp: The input, unfiltered transitions
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LinearResampler()
        """
        ...
    
    @overload
    def __init__(self, in_0: LinearResampler ) -> None:
        """
        Cython signature: void LinearResampler(LinearResampler &)
        """
        ...
    
    def raster(self, input: MSSpectrum ) -> None:
        """
        Cython signature: void raster(MSSpectrum & input)
        Applies the resampling algorithm to an MSSpectrum
        """
        ...
    
    def rasterExperiment(self, input: MSExperiment ) -> None:
        """
        Cython signature: void rasterExperiment(MSExperiment & input)
        Resamples the data in an MSExperiment
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ...
    
    def setLogType(self, in_0: int ) -> None:
        """
        Cython signature: void setLogType(LogType)
        Sets the progress log that should be used. The default type is NONE!
        """
        ...
    
    def getLogType(self) -> int:
        """
        Cython signature: LogType getLogType()
        Returns the type of progress log being used
        """
        ...
    
    def startProgress(self, begin: int , end: int , label: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void startProgress(ptrdiff_t begin, ptrdiff_t end, String label)
        """
        ...
    
    def setProgress(self, value: int ) -> None:
        """
        Cython signature: void setProgress(ptrdiff_t value)
        Sets the current progress
        """
        ...
    
    def endProgress(self) -> None:
        """
        Cython signature: void endProgress()
        Ends the progress display
        """
        ...
    
    def nextProgress(self) -> None:
        """
        Cython signature: void nextProgress()
        Increment progress by 1 (according to range begin-end)
        """
        ... 


class MsInspectFile:
    """
    Cython implementation of _MsInspectFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MsInspectFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MsInspectFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MsInspectFile ) -> None:
        """
        Cython signature: void MsInspectFile(MsInspectFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , feature_map: FeatureMap ) -> None:
        """
        Cython signature: void load(const String & filename, FeatureMap & feature_map)
        Loads a MsInspect file into a featureXML
        
        The content of the file is stored in `features`
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void store(const String & filename, MSSpectrum & spectrum)
        Stores a featureXML as a MsInspect file
        """
        ... 


class PepXMLFile:
    """
    Cython implementation of _PepXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PepXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void PepXMLFile()
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids)
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String experiment_name)
        """
        ...
    
    @overload
    def load(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , experiment_name: Union[bytes, str, String] , lookup: SpectrumMetaDataLookup ) -> None:
        """
        Cython signature: void load(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String experiment_name, SpectrumMetaDataLookup lookup)
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids)
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , protein_ids: List[ProteinIdentification] , peptide_ids: PeptideIdentificationList , mz_file: Union[bytes, str, String] , mz_name: Union[bytes, str, String] , peptideprophet_analyzed: bool , rt_tolerance: float ) -> None:
        """
        Cython signature: void store(String filename, libcpp_vector[ProteinIdentification] & protein_ids, PeptideIdentificationList & peptide_ids, String mz_file, String mz_name, bool peptideprophet_analyzed, double rt_tolerance)
        """
        ...
    
    def keepNativeSpectrumName(self, keep: bool ) -> None:
        """
        Cython signature: void keepNativeSpectrumName(bool keep)
        """
        ...
    
    def setParseUnknownScores(self, parse_unknown_scores: bool ) -> None:
        """
        Cython signature: void setParseUnknownScores(bool parse_unknown_scores)
        """
        ... 


class SqrtScaler:
    """
    Cython implementation of _SqrtScaler

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SqrtScaler.html>`_
      -- Inherits from ['DefaultParamHandler']

    Scales the intensity of peaks to the sqrt
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SqrtScaler()
        """
        ...
    
    @overload
    def __init__(self, in_0: SqrtScaler ) -> None:
        """
        Cython signature: void SqrtScaler(SqrtScaler &)
        """
        ...
    
    def filterSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class ThresholdMower:
    """
    Cython implementation of _ThresholdMower

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ThresholdMower.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ThresholdMower()
        """
        ...
    
    @overload
    def __init__(self, in_0: ThresholdMower ) -> None:
        """
        Cython signature: void ThresholdMower(ThresholdMower &)
        """
        ...
    
    def filterSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void filterPeakSpectrum(MSSpectrum & spec)
        """
        ...
    
    def filterPeakMap(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterPeakMap(MSExperiment & exp)
        """
        ...
    
    def getSubsections(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getSubsections()
        """
        ...
    
    def setParameters(self, param: Param ) -> None:
        """
        Cython signature: void setParameters(Param & param)
        Sets the parameters
        """
        ...
    
    def getParameters(self) -> Param:
        """
        Cython signature: Param getParameters()
        Returns the parameters
        """
        ...
    
    def getDefaults(self) -> Param:
        """
        Cython signature: Param getDefaults()
        Returns the default parameters
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class __ChromatogramType:
    None
    MASS_CHROMATOGRAM : int
    TOTAL_ION_CURRENT_CHROMATOGRAM : int
    SELECTED_ION_CURRENT_CHROMATOGRAM : int
    BASEPEAK_CHROMATOGRAM : int
    SELECTED_ION_MONITORING_CHROMATOGRAM : int
    SELECTED_REACTION_MONITORING_CHROMATOGRAM : int
    ELECTROMAGNETIC_RADIATION_CHROMATOGRAM : int
    ABSORPTION_CHROMATOGRAM : int
    EMISSION_CHROMATOGRAM : int
    SIZE_OF_CHROMATOGRAM_TYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

