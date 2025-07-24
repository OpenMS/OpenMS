from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class BiGaussModel:
    """
    Cython implementation of _BiGaussModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BiGaussModel.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void BiGaussModel()
        """
        ...
    
    @overload
    def __init__(self, in_0: BiGaussModel ) -> None:
        """
        Cython signature: void BiGaussModel(BiGaussModel &)
        """
        ...
    
    def setOffset(self, offset: float ) -> None:
        """
        Cython signature: void setOffset(double offset)
        """
        ...
    
    def setSamples(self) -> None:
        """
        Cython signature: void setSamples()
        """
        ...
    
    def getCenter(self) -> float:
        """
        Cython signature: double getCenter()
        """
        ... 


class DataProcessing:
    """
    Cython implementation of _DataProcessing

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DataProcessing.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DataProcessing()
        """
        ...
    
    @overload
    def __init__(self, in_0: DataProcessing ) -> None:
        """
        Cython signature: void DataProcessing(DataProcessing &)
        """
        ...
    
    def setProcessingActions(self, in_0: Set[int] ) -> None:
        """
        Cython signature: void setProcessingActions(libcpp_set[ProcessingAction])
        """
        ...
    
    def getProcessingActions(self) -> Set[int]:
        """
        Cython signature: libcpp_set[ProcessingAction] getProcessingActions()
        """
        ...
    
    def getSoftware(self) -> Software:
        """
        Cython signature: Software getSoftware()
        """
        ...
    
    def setSoftware(self, s: Software ) -> None:
        """
        Cython signature: void setSoftware(Software s)
        """
        ...
    
    def getCompletionTime(self) -> DateTime:
        """
        Cython signature: DateTime getCompletionTime()
        """
        ...
    
    def setCompletionTime(self, t: DateTime ) -> None:
        """
        Cython signature: void setCompletionTime(DateTime t)
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
    
    def __richcmp__(self, other: DataProcessing, op: int) -> Any:
        ...
    ProcessingAction : __ProcessingAction 


class EmgFitter1D:
    """
    Cython implementation of _EmgFitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EmgFitter1D.html>`_
      -- Inherits from ['LevMarqFitter1D']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EmgFitter1D()
        Exponentially modified gaussian distribution fitter (1-dim.) using Levenberg-Marquardt algorithm (Eigen implementation) for parameter optimization
        """
        ...
    
    @overload
    def __init__(self, in_0: EmgFitter1D ) -> None:
        """
        Cython signature: void EmgFitter1D(EmgFitter1D &)
        """
        ... 


class GridBasedCluster:
    """
    Cython implementation of _GridBasedCluster

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1GridBasedCluster.html>`_
    """
    
    @overload
    def __init__(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] , property_A: int , properties_B: List[int] ) -> None:
        """
        Cython signature: void GridBasedCluster(DPosition2 centre, DBoundingBox2 bounding_box, libcpp_vector[int] point_indices, int property_A, libcpp_vector[int] properties_B)
        """
        ...
    
    @overload
    def __init__(self, centre: Union[Sequence[int], Sequence[float]] , bounding_box: DBoundingBox2 , point_indices: List[int] ) -> None:
        """
        Cython signature: void GridBasedCluster(DPosition2 centre, DBoundingBox2 bounding_box, libcpp_vector[int] point_indices)
        """
        ...
    
    @overload
    def __init__(self, in_0: GridBasedCluster ) -> None:
        """
        Cython signature: void GridBasedCluster(GridBasedCluster &)
        """
        ...
    
    def getCentre(self) -> Union[Sequence[int], Sequence[float]]:
        """
        Cython signature: DPosition2 getCentre()
        Returns cluster centre
        """
        ...
    
    def getBoundingBox(self) -> DBoundingBox2:
        """
        Cython signature: DBoundingBox2 getBoundingBox()
        Returns bounding box
        """
        ...
    
    def getPoints(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getPoints()
        Returns indices of points in cluster
        """
        ...
    
    def getPropertyA(self) -> int:
        """
        Cython signature: int getPropertyA()
        Returns property A
        """
        ...
    
    def getPropertiesB(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getPropertiesB()
        Returns properties B of all points
        """
        ... 


class IDConflictResolverAlgorithm:
    """
    Cython implementation of _IDConflictResolverAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDConflictResolverAlgorithm.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IDConflictResolverAlgorithm()
        Resolves ambiguous annotations of features with peptide identifications
        """
        ...
    
    @overload
    def __init__(self, in_0: IDConflictResolverAlgorithm ) -> None:
        """
        Cython signature: void IDConflictResolverAlgorithm(IDConflictResolverAlgorithm &)
        """
        ...
    
    @overload
    def resolve(self, features: FeatureMap ) -> None:
        """
        Cython signature: void resolve(FeatureMap & features)
        Resolves ambiguous annotations of features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
        """
        ...
    
    @overload
    def resolve(self, features: ConsensusMap ) -> None:
        """
        Cython signature: void resolve(ConsensusMap & features)
        Resolves ambiguous annotations of consensus features with peptide identifications\n
        
        The the filtered identifications are added to the vector of unassigned peptides
        and also reduced to a single best hit
        
        
        :param keep_matching: Keeps all IDs that match the modified sequence of the best hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
        """
        ...
    
    @overload
    def resolveBetweenFeatures(self, features: FeatureMap ) -> None:
        """
        Cython signature: void resolveBetweenFeatures(FeatureMap & features)
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
        """
        ...
    
    @overload
    def resolveBetweenFeatures(self, features: ConsensusMap ) -> None:
        """
        Cython signature: void resolveBetweenFeatures(ConsensusMap & features)
        In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear\n
        
        This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
        combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
        appears only once, i.e. no multiplicities
        """
        ... 


class MSExperiment:
    """
    Cython implementation of _MSExperiment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSExperiment.html>`_
      -- Inherits from ['ExperimentalSettings']

    In-Memory representation of a mass spectrometry experiment.
    
    Contains the data and metadata of an experiment performed with an MS (or
    HPLC and MS). This representation of an MS experiment is organized as list
    of spectra and chromatograms and provides an in-memory representation of
    popular mass-spectrometric file formats such as mzXML or mzML. The
    meta-data associated with an experiment is contained in
    ExperimentalSettings (by inheritance) while the raw data (as well as
    spectra and chromatogram level meta data) is stored in objects of type
    MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
    and getChromatogram functions.
    
    Spectra can be accessed by direct iteration or by getSpectrum(),
    while chromatograms are accessed through getChromatogram().
    See help(ExperimentalSettings) for information about meta-data.
    
    Usage:
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSExperiment()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSExperiment ) -> None:
        """
        Cython signature: void MSExperiment(MSExperiment &)
        """
        ...
    
    def getExperimentalSettings(self) -> ExperimentalSettings:
        """
        Cython signature: ExperimentalSettings getExperimentalSettings()
        """
        ...
    
    def __getitem__(self, in_0: int ) -> MSSpectrum:
        """
        Cython signature: MSSpectrum & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: MSSpectrum ) -> None:
        """Cython signature: MSSpectrum & operator[](size_t)"""
        ...
    
    def addSpectrum(self, spec: MSSpectrum ) -> None:
        """
        Cython signature: void addSpectrum(MSSpectrum spec)
        """
        ...
    
    def setSpectra(self, spectra: List[MSSpectrum] ) -> None:
        """
        Cython signature: void setSpectra(libcpp_vector[MSSpectrum] & spectra)
        """
        ...
    
    def getSpectra(self) -> List[MSSpectrum]:
        """
        Cython signature: libcpp_vector[MSSpectrum] getSpectra()
        """
        ...
    
    def aggregateFromMatrix(self, ranges: MatrixDouble , ms_level: int , mz_agg: bytes ) -> List[List[float]]:
        """
        Cython signature: libcpp_vector[libcpp_vector[double]] aggregateFromMatrix(MatrixDouble & ranges, unsigned int ms_level, libcpp_string mz_agg)
        Aggregates intensity values for multiple m/z and RT ranges specified in a matrix
        """
        ...
    
    def extractXICsFromMatrix(self, ranges: MatrixDouble , ms_level: int , mz_agg: bytes ) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] extractXICsFromMatrix(MatrixDouble & ranges, unsigned int ms_level, libcpp_string mz_agg)
        Extracts XIC chromatograms for multiple m/z and RT ranges specified in a matrix
        """
        ...
    
    def addChromatogram(self, chromatogram: MSChromatogram ) -> None:
        """
        Cython signature: void addChromatogram(MSChromatogram chromatogram)
        """
        ...
    
    def setChromatograms(self, chromatograms: List[MSChromatogram] ) -> None:
        """
        Cython signature: void setChromatograms(libcpp_vector[MSChromatogram] chromatograms)
        """
        ...
    
    def getChromatograms(self) -> List[MSChromatogram]:
        """
        Cython signature: libcpp_vector[MSChromatogram] getChromatograms()
        """
        ...
    
    def calculateTIC(self) -> MSChromatogram:
        """
        Cython signature: MSChromatogram calculateTIC()
        Returns the total ion chromatogram
        """
        ...
    
    def clear(self, clear_meta_data: bool ) -> None:
        """
        Cython signature: void clear(bool clear_meta_data)
        Clear all spectra data and meta data (if called with True)
        """
        ...
    
    def updateRanges(self) -> None:
        """
        Cython signature: void updateRanges()
        Recalculate global ranges for both spectra and chromatrograms after changes to the data has been made.
        """
        ...
    
    def reserveSpaceSpectra(self, s: int ) -> None:
        """
        Cython signature: void reserveSpaceSpectra(size_t s)
        """
        ...
    
    def reserveSpaceChromatograms(self, s: int ) -> None:
        """
        Cython signature: void reserveSpaceChromatograms(size_t s)
        """
        ...
    
    def getSize(self) -> int:
        """
        Cython signature: uint64_t getSize()
        Returns the total number of peaks
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def resize(self, s: int ) -> None:
        """
        Cython signature: void resize(size_t s)
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        """
        ...
    
    def reserve(self, s: int ) -> None:
        """
        Cython signature: void reserve(size_t s)
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        Returns the number of MS spectra
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        Returns the number of chromatograms
        """
        ...
    
    @overload
    def sortSpectra(self, sort_mz: bool ) -> None:
        """
        Cython signature: void sortSpectra(bool sort_mz)
        Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z
        """
        ...
    
    @overload
    def sortSpectra(self, ) -> None:
        """
        Cython signature: void sortSpectra()
        """
        ...
    
    @overload
    def sortChromatograms(self, sort_rt: bool ) -> None:
        """
        Cython signature: void sortChromatograms(bool sort_rt)
        Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT
        """
        ...
    
    @overload
    def sortChromatograms(self, ) -> None:
        """
        Cython signature: void sortChromatograms()
        """
        ...
    
    @overload
    def isSorted(self, check_mz: bool ) -> bool:
        """
        Cython signature: bool isSorted(bool check_mz)
        Checks if all spectra are sorted with respect to ascending RT
        """
        ...
    
    @overload
    def isSorted(self, ) -> bool:
        """
        Cython signature: bool isSorted()
        """
        ...
    
    def getPrimaryMSRunPath(self, toFill: List[bytes] ) -> None:
        """
        Cython signature: void getPrimaryMSRunPath(StringList & toFill)
        References to the first MS file(s) after conversions. Used to trace results back to original data.
        """
        ...
    
    def swap(self, in_0: MSExperiment ) -> None:
        """
        Cython signature: void swap(MSExperiment)
        """
        ...
    
    def reset(self) -> None:
        """
        Cython signature: void reset()
        """
        ...
    
    def clearMetaDataArrays(self) -> bool:
        """
        Cython signature: bool clearMetaDataArrays()
        """
        ...
    
    def getPrecursorSpectrum(self, zero_based_index: int ) -> int:
        """
        Cython signature: int getPrecursorSpectrum(int zero_based_index)
        Returns the index of the precursor spectrum for spectrum at index @p zero_based_index
        """
        ...
    
    def spectrumRanges(self) -> SpectrumRangeManager:
        """
        Cython signature: SpectrumRangeManager spectrumRanges()
        Returns a reference to the spectrum range manager
        """
        ...
    
    def chromatogramRanges(self) -> ChromatogramRangeManager:
        """
        Cython signature: ChromatogramRangeManager chromatogramRanges()
        Returns a reference to the chromatogram range manager
        """
        ...
    
    def getMinRT(self) -> float:
        """
        Cython signature: double getMinRT()
        Get the minimum RT value from the combined ranges
        """
        ...
    
    def getMaxRT(self) -> float:
        """
        Cython signature: double getMaxRT()
        Get the maximum RT value from the combined ranges
        """
        ...
    
    def getMinMZ(self) -> float:
        """
        Cython signature: double getMinMZ()
        Get the minimum m/z value from the combined ranges
        """
        ...
    
    def getMaxMZ(self) -> float:
        """
        Cython signature: double getMaxMZ()
        Get the maximum m/z value from the combined ranges
        """
        ...
    
    def getMinIntensity(self) -> float:
        """
        Cython signature: double getMinIntensity()
        Get the minimum intensity value from the combined ranges
        """
        ...
    
    def getMaxIntensity(self) -> float:
        """
        Cython signature: double getMaxIntensity()
        Get the maximum intensity value from the combined ranges
        """
        ...
    
    def getMinMobility(self) -> float:
        """
        Cython signature: double getMinMobility()
        Get the minimum mobility value from the combined ranges
        """
        ...
    
    def getMaxMobility(self) -> float:
        """
        Cython signature: double getMaxMobility()
        Get the maximum mobility value from the combined ranges
        """
        ...
    
    def clearRanges(self) -> None:
        """
        Cython signature: void clearRanges()
        Clear all ranges in all range managers
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
    
    def __richcmp__(self, other: MSExperiment, op: int) -> Any:
        ...
    
    def __iter__(self) -> MSSpectrum:
       ... 


class MzTab:
    """
    Cython implementation of _MzTab

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTab.html>`_

    Data model of MzTab files
    
    Please see the official MzTab specification at https://code.google.com/p/mztab/
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzTab()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzTab ) -> None:
        """
        Cython signature: void MzTab(MzTab &)
        """
        ... 


class OSW_ChromExtractParams:
    """
    Cython implementation of _OSW_ChromExtractParams

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OSW_ChromExtractParams.html>`_
    """
    
    min_upper_edge_dist: float
    
    mz_extraction_window: float
    
    ppm: bool
    
    extraction_function: bytes
    
    rt_extraction_window: float
    
    extra_rt_extract: float
    
    im_extraction_window: float
    
    def __init__(self, in_0: OSW_ChromExtractParams ) -> None:
        """
        Cython signature: void OSW_ChromExtractParams(OSW_ChromExtractParams &)
        """
        ... 


class OpenSwathOSWWriter:
    """
    Cython implementation of _OpenSwathOSWWriter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OpenSwathOSWWriter.html>`_
    """
    
    @overload
    def __init__(self, output_filename: Union[bytes, str, String] , run_id: int , input_filename: Union[bytes, str, String] , uis_scores: bool ) -> None:
        """
        Cython signature: void OpenSwathOSWWriter(String output_filename, uint64_t run_id, String input_filename, bool uis_scores)
        """
        ...
    
    @overload
    def __init__(self, in_0: OpenSwathOSWWriter ) -> None:
        """
        Cython signature: void OpenSwathOSWWriter(OpenSwathOSWWriter &)
        """
        ...
    
    def isActive(self) -> bool:
        """
        Cython signature: bool isActive()
        """
        ...
    
    def writeHeader(self) -> None:
        """
        Cython signature: void writeHeader()
        Initializes file by generating SQLite tables
        """
        ...
    
    def prepareLine(self, compound: LightCompound , tr: LightTransition , output: FeatureMap , id_: Union[bytes, str, String] ) -> Union[bytes, str, String]:
        """
        Cython signature: String prepareLine(LightCompound & compound, LightTransition * tr, FeatureMap & output, String id_)
        Prepare a single line (feature) for output
        
        The result can be flushed to disk using writeLines (either line by line or after collecting several lines)
        
        
        :param pep: The compound (peptide/metabolite) used for extraction
        :param transition: The transition used for extraction
        :param output: The feature map containing all features (each feature will generate one entry in the output)
        :param id: The transition group identifier (peptide/metabolite id)
        :return: A String to be written using writeLines
        """
        ...
    
    def writeLines(self, to_osw_output: List[bytes] ) -> None:
        """
        Cython signature: void writeLines(libcpp_vector[String] to_osw_output)
        Write data to disk
        
        Takes a set of pre-prepared data statements from prepareLine and flushes them to disk
        
        
        :param to_osw_output: Statements generated by prepareLine
        """
        ... 


class PeakTypeEstimator:
    """
    Cython implementation of _PeakTypeEstimator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakTypeEstimator.html>`_

    Estimates if the data of a spectrum is raw data or peak data
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakTypeEstimator()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakTypeEstimator ) -> None:
        """
        Cython signature: void PeakTypeEstimator(PeakTypeEstimator &)
        """
        ... 


class SequestInfile:
    """
    Cython implementation of _SequestInfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SequestInfile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SequestInfile()
        Sequest input file adapter
        """
        ...
    
    @overload
    def __init__(self, in_0: SequestInfile ) -> None:
        """
        Cython signature: void SequestInfile(SequestInfile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void store(const String & filename)
        Stores the experiment data in a Sequest input file that can be used as input for Sequest shell execution
        
        :param filename: the name of the file in which the infile is stored into
        """
        ...
    
    def getEnzymeInfoAsString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEnzymeInfoAsString()
        Returns the enzyme list as a string
        """
        ...
    
    def getDatabase(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getDatabase()
        Returns the used database
        """
        ...
    
    def setDatabase(self, database: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setDatabase(const String & database)
        Sets the used database
        """
        ...
    
    def getNeutralLossesForIons(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getNeutralLossesForIons()
        Returns whether neutral losses are considered for the a-, b- and y-ions
        """
        ...
    
    def setNeutralLossesForIons(self, neutral_losses_for_ions: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNeutralLossesForIons(const String & neutral_losses_for_ions)
        Sets whether neutral losses are considered for the a-, b- and y-ions
        """
        ...
    
    def getIonSeriesWeights(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIonSeriesWeights()
        Returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        """
        ...
    
    def setIonSeriesWeights(self, ion_series_weights: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIonSeriesWeights(const String & ion_series_weights)
        Sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series
        """
        ...
    
    def getPartialSequence(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getPartialSequence()
        Returns the partial sequences (space delimited) that have to occur in the theoretical spectra
        """
        ...
    
    def setPartialSequence(self, partial_sequence: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setPartialSequence(const String & partial_sequence)
        Sets the partial sequences (space delimited) that have to occur in the theoretical spectra
        """
        ...
    
    def getSequenceHeaderFilter(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getSequenceHeaderFilter()
        Returns the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        """
        ...
    
    def setSequenceHeaderFilter(self, sequence_header_filter: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setSequenceHeaderFilter(const String & sequence_header_filter)
        Sets the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered
        """
        ...
    
    def getProteinMassFilter(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getProteinMassFilter()
        Returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        """
        ...
    
    def setProteinMassFilter(self, protein_mass_filter: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setProteinMassFilter(const String & protein_mass_filter)
        Sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)
        """
        ...
    
    def getPeakMassTolerance(self) -> float:
        """
        Cython signature: float getPeakMassTolerance()
        Returns the peak mass tolerance
        """
        ...
    
    def setPeakMassTolerance(self, peak_mass_tolerance: float ) -> None:
        """
        Cython signature: void setPeakMassTolerance(float peak_mass_tolerance)
        Sets the peak mass tolerance
        """
        ...
    
    def getPrecursorMassTolerance(self) -> float:
        """
        Cython signature: float getPrecursorMassTolerance()
        Returns the precursor mass tolerance
        """
        ...
    
    def setPrecursorMassTolerance(self, precursor_mass_tolerance: float ) -> None:
        """
        Cython signature: void setPrecursorMassTolerance(float precursor_mass_tolerance)
        Sets the precursor mass tolerance
        """
        ...
    
    def getMatchPeakTolerance(self) -> float:
        """
        Cython signature: float getMatchPeakTolerance()
        Returns the match peak tolerance
        """
        ...
    
    def setMatchPeakTolerance(self, match_peak_tolerance: float ) -> None:
        """
        Cython signature: void setMatchPeakTolerance(float match_peak_tolerance)
        Sets the match peak tolerance
        """
        ...
    
    def getIonCutoffPercentage(self) -> float:
        """
        Cython signature: float getIonCutoffPercentage()
        Returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks
        """
        ...
    
    def setIonCutoffPercentage(self, ion_cutoff_percentage: float ) -> None:
        """
        Cython signature: void setIonCutoffPercentage(float ion_cutoff_percentage)
        Sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks
        """
        ...
    
    def getPeptideMassUnit(self) -> int:
        """
        Cython signature: size_t getPeptideMassUnit()
        Returns the peptide mass unit
        """
        ...
    
    def setPeptideMassUnit(self, peptide_mass_unit: int ) -> None:
        """
        Cython signature: void setPeptideMassUnit(size_t peptide_mass_unit)
        Sets the peptide mass unit
        """
        ...
    
    def getOutputLines(self) -> int:
        """
        Cython signature: size_t getOutputLines()
        Returns the number of peptides to be displayed
        """
        ...
    
    def setOutputLines(self, output_lines: int ) -> None:
        """
        Cython signature: void setOutputLines(size_t output_lines)
        Sets the number of peptides to be displayed
        """
        ...
    
    def getEnzymeNumber(self) -> int:
        """
        Cython signature: size_t getEnzymeNumber()
        Returns the enzyme used for cleavage (by means of the number from a list of enzymes)
        """
        ...
    
    def getEnzymeName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEnzymeName()
        Returns the enzyme used for cleavage
        """
        ...
    
    def setEnzyme(self, enzyme_name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t setEnzyme(String enzyme_name)
        Sets the enzyme used for cleavage (by means of the number from a list of enzymes)
        """
        ...
    
    def getMaxAAPerModPerPeptide(self) -> int:
        """
        Cython signature: size_t getMaxAAPerModPerPeptide()
        Returns the maximum number of amino acids containing the same modification in a peptide
        """
        ...
    
    def setMaxAAPerModPerPeptide(self, max_aa_per_mod_per_peptide: int ) -> None:
        """
        Cython signature: void setMaxAAPerModPerPeptide(size_t max_aa_per_mod_per_peptide)
        Sets the maximum number of amino acids containing the same modification in a peptide
        """
        ...
    
    def getMaxModsPerPeptide(self) -> int:
        """
        Cython signature: size_t getMaxModsPerPeptide()
        Returns the maximum number of modifications that are allowed in a peptide
        """
        ...
    
    def setMaxModsPerPeptide(self, max_mods_per_peptide: int ) -> None:
        """
        Cython signature: void setMaxModsPerPeptide(size_t max_mods_per_peptide)
        Sets the maximum number of modifications that are allowed in a peptide
        """
        ...
    
    def getNucleotideReadingFrame(self) -> int:
        """
        Cython signature: size_t getNucleotideReadingFrame()
        Returns the nucleotide reading frame
        """
        ...
    
    def setNucleotideReadingFrame(self, nucleotide_reading_frame: int ) -> None:
        """
        Cython signature: void setNucleotideReadingFrame(size_t nucleotide_reading_frame)
        Sets the nucleotide reading frame
        """
        ...
    
    def getMaxInternalCleavageSites(self) -> int:
        """
        Cython signature: size_t getMaxInternalCleavageSites()
        Returns the maximum number of internal cleavage sites
        """
        ...
    
    def setMaxInternalCleavageSites(self, max_internal_cleavage_sites: int ) -> None:
        """
        Cython signature: void setMaxInternalCleavageSites(size_t max_internal_cleavage_sites)
        Sets the maximum number of internal cleavage sites
        """
        ...
    
    def getMatchPeakCount(self) -> int:
        """
        Cython signature: size_t getMatchPeakCount()
        Returns the number of top abundant peaks to match with theoretical ones
        """
        ...
    
    def setMatchPeakCount(self, match_peak_count: int ) -> None:
        """
        Cython signature: void setMatchPeakCount(size_t match_peak_count)
        Sets the number of top abundant peaks to with theoretical ones
        """
        ...
    
    def getMatchPeakAllowedError(self) -> int:
        """
        Cython signature: size_t getMatchPeakAllowedError()
        Returns the number of top abundant peaks that are allowed not to match with a theoretical peak
        """
        ...
    
    def setMatchPeakAllowedError(self, match_peak_allowed_error: int ) -> None:
        """
        Cython signature: void setMatchPeakAllowedError(size_t match_peak_allowed_error)
        Sets the number of top abundant peaks that are allowed not to match with a theoretical peak
        """
        ...
    
    def getShowFragmentIons(self) -> bool:
        """
        Cython signature: bool getShowFragmentIons()
        Returns whether fragment ions shall be displayed
        """
        ...
    
    def setShowFragmentIons(self, show_fragments: bool ) -> None:
        """
        Cython signature: void setShowFragmentIons(bool show_fragments)
        Sets whether fragment ions shall be displayed
        """
        ...
    
    def getPrintDuplicateReferences(self) -> bool:
        """
        Cython signature: bool getPrintDuplicateReferences()
        Returns whether all proteins containing a found peptide should be displayed
        """
        ...
    
    def setPrintDuplicateReferences(self, print_duplicate_references: bool ) -> None:
        """
        Cython signature: void setPrintDuplicateReferences(bool print_duplicate_references)
        Sets whether all proteins containing a found peptide should be displayed
        """
        ...
    
    def getRemovePrecursorNearPeaks(self) -> bool:
        """
        Cython signature: bool getRemovePrecursorNearPeaks()
        Returns whether peaks near (15 amu) the precursor peak are removed
        """
        ...
    
    def setRemovePrecursorNearPeaks(self, remove_precursor_near_peaks: bool ) -> None:
        """
        Cython signature: void setRemovePrecursorNearPeaks(bool remove_precursor_near_peaks)
        Sets whether peaks near (15 amu) the precursor peak are removed
        """
        ...
    
    def getMassTypeParent(self) -> bool:
        """
        Cython signature: bool getMassTypeParent()
        Returns the mass type of the parent (0 - monoisotopic, 1 - average mass)
        """
        ...
    
    def setMassTypeParent(self, mass_type_parent: bool ) -> None:
        """
        Cython signature: void setMassTypeParent(bool mass_type_parent)
        Sets the mass type of the parent (0 - monoisotopic, 1 - average mass)
        """
        ...
    
    def getMassTypeFragment(self) -> bool:
        """
        Cython signature: bool getMassTypeFragment()
        Returns the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        """
        ...
    
    def setMassTypeFragment(self, mass_type_fragment: bool ) -> None:
        """
        Cython signature: void setMassTypeFragment(bool mass_type_fragment)
        Sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)
        """
        ...
    
    def getNormalizeXcorr(self) -> bool:
        """
        Cython signature: bool getNormalizeXcorr()
        Returns whether normalized xcorr values are displayed
        """
        ...
    
    def setNormalizeXcorr(self, normalize_xcorr: bool ) -> None:
        """
        Cython signature: void setNormalizeXcorr(bool normalize_xcorr)
        Sets whether normalized xcorr values are displayed
        """
        ...
    
    def getResiduesInUpperCase(self) -> bool:
        """
        Cython signature: bool getResiduesInUpperCase()
        Returns whether residues are in upper case
        """
        ...
    
    def setResiduesInUpperCase(self, residues_in_upper_case: bool ) -> None:
        """
        Cython signature: void setResiduesInUpperCase(bool residues_in_upper_case)
        Sets whether residues are in upper case
        """
        ...
    
    def addEnzymeInfo(self, enzyme_info: List[bytes] ) -> None:
        """
        Cython signature: void addEnzymeInfo(libcpp_vector[String] & enzyme_info)
        Adds an enzyme to the list and sets is as used
        """
        ...
    
    def handlePTMs(self, modification_line: Union[bytes, str, String] , modifications_filename: Union[bytes, str, String] , monoisotopic: bool ) -> None:
        """
        Cython signature: void handlePTMs(const String & modification_line, const String & modifications_filename, bool monoisotopic)
        """
        ...
    
    def __richcmp__(self, other: SequestInfile, op: int) -> Any:
        ... 


class TMTTenPlexQuantitationMethod:
    """
    Cython implementation of _TMTTenPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTTenPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TMTTenPlexQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: TMTTenPlexQuantitationMethod ) -> None:
        """
        Cython signature: void TMTTenPlexQuantitationMethod(TMTTenPlexQuantitationMethod &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        """
        ...
    
    def getChannelInformation(self) -> List[IsobaricChannelInformation]:
        """
        Cython signature: libcpp_vector[IsobaricChannelInformation] getChannelInformation()
        """
        ...
    
    def getNumberOfChannels(self) -> int:
        """
        Cython signature: size_t getNumberOfChannels()
        """
        ...
    
    def getIsotopeCorrectionMatrix(self) -> MatrixDouble:
        """
        Cython signature: MatrixDouble getIsotopeCorrectionMatrix()
        """
        ...
    
    def getReferenceChannel(self) -> int:
        """
        Cython signature: size_t getReferenceChannel()
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
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String &)
        Sets the name
        """
        ... 


class TheoreticalSpectrumGeneratorXLMS:
    """
    Cython implementation of _TheoreticalSpectrumGeneratorXLMS

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TheoreticalSpectrumGeneratorXLMS.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TheoreticalSpectrumGeneratorXLMS()
        """
        ...
    
    @overload
    def __init__(self, in_0: TheoreticalSpectrumGeneratorXLMS ) -> None:
        """
        Cython signature: void TheoreticalSpectrumGeneratorXLMS(TheoreticalSpectrumGeneratorXLMS &)
        """
        ...
    
    def getLinearIonSpectrum(self, spectrum: MSSpectrum , peptide: AASequence , link_pos: int , frag_alpha: bool , charge: int , link_pos_2: int ) -> None:
        """
        Cython signature: void getLinearIonSpectrum(MSSpectrum & spectrum, AASequence peptide, size_t link_pos, bool frag_alpha, int charge, size_t link_pos_2)
            Generates fragment ions not containing the cross-linker for one peptide
        
            B-ions are generated from the beginning of the peptide up to the first linked position,
            y-ions are generated from the second linked position up the end of the peptide.
            If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
            For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param peptide: The peptide to fragment
            :param link_pos: The position of the cross-linker on the given peptide
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation
            :param charge: The maximal charge of the ions
            :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        ...
    
    @overload
    def getXLinkIonSpectrum(self, spectrum: MSSpectrum , peptide: AASequence , link_pos: int , precursor_mass: float , frag_alpha: bool , mincharge: int , maxcharge: int , link_pos_2: int ) -> None:
        """
        Cython signature: void getXLinkIonSpectrum(MSSpectrum & spectrum, AASequence peptide, size_t link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, size_t link_pos_2)
            Generates fragment ions containing the cross-linker for one peptide
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
            For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
            Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
            a precursor mass for the two peptides and the linker is needed.
            In the case of a loop link the precursor mass is the mass of the only peptide and the linker.
            Although this function is more general, currently it is mainly used for loop-links and mono-links,
            because residues in the second, unknown peptide cannot be considered for possible neutral losses.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param peptide: The peptide to fragment
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
            :param link_pos_2: A second position for the linker, in case it is a loop link
        """
        ...
    
    @overload
    def getXLinkIonSpectrum(self, spectrum: MSSpectrum , crosslink: ProteinProteinCrossLink , frag_alpha: bool , mincharge: int , maxcharge: int ) -> None:
        """
        Cython signature: void getXLinkIonSpectrum(MSSpectrum & spectrum, ProteinProteinCrossLink crosslink, bool frag_alpha, int mincharge, int maxcharge)
            Generates fragment ions containing the cross-linker for a pair of peptides
        
            B-ions are generated from the first linked position up to the end of the peptide,
            y-ions are generated from the beginning of the peptide up to the second linked position.
            This function generates neutral loss ions by considering both linked peptides.
            Only one of the peptides, decided by @frag_alpha, is fragmented.
            This function is not suitable to generate fragments for mono-links or loop-links.
            This simplifies the function, but it has to be called twice to get all fragments of a peptide pair.
            The generated ion types and other additional settings are determined by the tool parameters
        
            :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
            :param crosslink: ProteinProteinCrossLink to be fragmented
            :param link_pos: The position of the cross-linker on the given peptide
            :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
            :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
            :param mincharge: The minimal charge of the ions
            :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
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


class __ProcessingAction:
    None
    DATA_PROCESSING : int
    CHARGE_DECONVOLUTION : int
    DEISOTOPING : int
    SMOOTHING : int
    CHARGE_CALCULATION : int
    PRECURSOR_RECALCULATION : int
    BASELINE_REDUCTION : int
    PEAK_PICKING : int
    ALIGNMENT : int
    CALIBRATION : int
    NORMALIZATION : int
    FILTERING : int
    QUANTITATION : int
    FEATURE_GROUPING : int
    IDENTIFICATION_MAPPING : int
    FORMAT_CONVERSION : int
    CONVERSION_MZDATA : int
    CONVERSION_MZML : int
    CONVERSION_MZXML : int
    CONVERSION_DTA : int
    IDENTIFICATION : int
    SIZE_OF_PROCESSINGACTION : int

    def getMapping(self) -> Dict[int, str]:
       ... 

