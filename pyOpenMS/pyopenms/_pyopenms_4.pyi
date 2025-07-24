from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class CVMappingRule:
    """
    Cython implementation of _CVMappingRule

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingRule.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVMappingRule()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVMappingRule ) -> None:
        """
        Cython signature: void CVMappingRule(CVMappingRule &)
        """
        ...
    
    def setIdentifier(self, identifier: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(String identifier)
        Sets the identifier of the rule
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        Returns the identifier of the rule
        """
        ...
    
    def setElementPath(self, element_path: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setElementPath(String element_path)
        Sets the path of the DOM element, where this rule is allowed
        """
        ...
    
    def getElementPath(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getElementPath()
        Returns the path of the DOM element, where this rule is allowed
        """
        ...
    
    def setRequirementLevel(self, level: int ) -> None:
        """
        Cython signature: void setRequirementLevel(RequirementLevel level)
        Sets the requirement level of this rule
        """
        ...
    
    def getRequirementLevel(self) -> int:
        """
        Cython signature: RequirementLevel getRequirementLevel()
        Returns the requirement level of this rule
        """
        ...
    
    def setCombinationsLogic(self, combinations_logic: int ) -> None:
        """
        Cython signature: void setCombinationsLogic(CombinationsLogic combinations_logic)
        Sets the combination operator of the rule
        """
        ...
    
    def getCombinationsLogic(self) -> int:
        """
        Cython signature: CombinationsLogic getCombinationsLogic()
        Returns the combinations operator of the rule
        """
        ...
    
    def setScopePath(self, path: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setScopePath(String path)
        Sets the scope path of the rule
        """
        ...
    
    def getScopePath(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getScopePath()
        Returns the scope path of the rule
        """
        ...
    
    def setCVTerms(self, cv_terms: List[CVMappingTerm] ) -> None:
        """
        Cython signature: void setCVTerms(libcpp_vector[CVMappingTerm] cv_terms)
        Sets the terms which are allowed
        """
        ...
    
    def getCVTerms(self) -> List[CVMappingTerm]:
        """
        Cython signature: libcpp_vector[CVMappingTerm] getCVTerms()
        Returns the allowed terms
        """
        ...
    
    def addCVTerm(self, cv_terms: CVMappingTerm ) -> None:
        """
        Cython signature: void addCVTerm(CVMappingTerm cv_terms)
        Adds a term to the allowed terms
        """
        ...
    
    def __richcmp__(self, other: CVMappingRule, op: int) -> Any:
        ...
    CombinationsLogic : __CombinationsLogic
    RequirementLevel : __RequirementLevel 


class Fitter1D:
    """
    Cython implementation of _Fitter1D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Fitter1D.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Fitter1D()
        Abstract base class for all 1D-dimensional model fitter
        """
        ...
    
    @overload
    def __init__(self, in_0: Fitter1D ) -> None:
        """
        Cython signature: void Fitter1D(Fitter1D &)
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


class InspectOutfile:
    """
    Cython implementation of _InspectOutfile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1InspectOutfile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void InspectOutfile()
        This class serves to read in an Inspect outfile and write an idXML file
        """
        ...
    
    @overload
    def __init__(self, in_0: InspectOutfile ) -> None:
        """
        Cython signature: void InspectOutfile(InspectOutfile &)
        """
        ...
    
    def load(self, result_filename: Union[bytes, str, String] , peptide_identifications: PeptideIdentificationList , protein_identification: ProteinIdentification , p_value_threshold: float , database_filename: Union[bytes, str, String] ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] load(const String & result_filename, PeptideIdentificationList & peptide_identifications, ProteinIdentification & protein_identification, double p_value_threshold, const String & database_filename)
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
        ...
    
    def getWantedRecords(self, result_filename: Union[bytes, str, String] , p_value_threshold: float ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] getWantedRecords(const String & result_filename, double p_value_threshold)
        Loads only results which exceeds a given p-value threshold
        
        
        :param result_filename: The filename of the results file
        :param p_value_threshold: Only identifications exceeding this threshold are read
        :raises:
          Exception: FileNotFound is thrown if the given file could not be found
        :raises:
          Exception: FileEmpty is thrown if the given file is empty
        """
        ...
    
    def compressTrieDB(self, database_filename: Union[bytes, str, String] , index_filename: Union[bytes, str, String] , wanted_records: List[int] , snd_database_filename: Union[bytes, str, String] , snd_index_filename: Union[bytes, str, String] , append: bool ) -> None:
        """
        Cython signature: void compressTrieDB(const String & database_filename, const String & index_filename, libcpp_vector[size_t] & wanted_records, const String & snd_database_filename, const String & snd_index_filename, bool append)
        Generates a trie database from another one, using the wanted records only
        """
        ...
    
    def generateTrieDB(self, source_database_filename: Union[bytes, str, String] , database_filename: Union[bytes, str, String] , index_filename: Union[bytes, str, String] , append: bool , species: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void generateTrieDB(const String & source_database_filename, const String & database_filename, const String & index_filename, bool append, const String species)
        Generates a trie database from a given one (the type of database is determined by getLabels)
        """
        ...
    
    def getACAndACType(self, line: Union[bytes, str, String] , accession: String , accession_type: String ) -> None:
        """
        Cython signature: void getACAndACType(String line, String & accession, String & accession_type)
        Retrieve the accession type and accession number from a protein description line
        """
        ...
    
    def getLabels(self, source_database_filename: Union[bytes, str, String] , ac_label: String , sequence_start_label: String , sequence_end_label: String , comment_label: String , species_label: String ) -> None:
        """
        Cython signature: void getLabels(const String & source_database_filename, String & ac_label, String & sequence_start_label, String & sequence_end_label, String & comment_label, String & species_label)
        Retrieve the labels of a given database (at the moment FASTA and Swissprot)
        """
        ...
    
    def getSequences(self, database_filename: Union[bytes, str, String] , wanted_records: Dict[int, int] , sequences: List[bytes] ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] getSequences(const String & database_filename, libcpp_map[size_t,size_t] & wanted_records, libcpp_vector[String] & sequences)
        Retrieve sequences from a trie database
        """
        ...
    
    def getExperiment(self, exp: MSExperiment , type_: String , in_filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void getExperiment(MSExperiment & exp, String & type_, const String & in_filename)
        Get the experiment from a file
        """
        ...
    
    def getSearchEngineAndVersion(self, cmd_output: Union[bytes, str, String] , protein_identification: ProteinIdentification ) -> bool:
        """
        Cython signature: bool getSearchEngineAndVersion(const String & cmd_output, ProteinIdentification & protein_identification)
        Get the search engine and its version from the output of the InsPecT executable without parameters. Returns true on success, false otherwise
        """
        ...
    
    def readOutHeader(self, filename: Union[bytes, str, String] , header_line: Union[bytes, str, String] , spectrum_file_column: int , scan_column: int , peptide_column: int , protein_column: int , charge_column: int , MQ_score_column: int , p_value_column: int , record_number_column: int , DB_file_pos_column: int , spec_file_pos_column: int , number_of_columns: int ) -> None:
        """
        Cython signature: void readOutHeader(const String & filename, const String & header_line, int & spectrum_file_column, int & scan_column, int & peptide_column, int & protein_column, int & charge_column, int & MQ_score_column, int & p_value_column, int & record_number_column, int & DB_file_pos_column, int & spec_file_pos_column, size_t & number_of_columns)
        Read the header of an inspect output file and retrieve various information
        """
        ...
    
    def __richcmp__(self, other: InspectOutfile, op: int) -> Any:
        ... 


class IntegerDataArray:
    """
    Cython implementation of _IntegerDataArray

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::DataArrays_1_1IntegerDataArray.html>`_
      -- Inherits from ['MetaInfoDescription']

    The representation of extra integer data attached to a spectrum or chromatogram.
    Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IntegerDataArray()
        """
        ...
    
    @overload
    def __init__(self, in_0: IntegerDataArray ) -> None:
        """
        Cython signature: void IntegerDataArray(IntegerDataArray &)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        """
        ...
    
    def resize(self, n: int ) -> None:
        """
        Cython signature: void resize(size_t n)
        """
        ...
    
    def reserve(self, n: int ) -> None:
        """
        Cython signature: void reserve(size_t n)
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        """
        ...
    
    def push_back(self, in_0: int ) -> None:
        """
        Cython signature: void push_back(int)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the peak annotations
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the peak annotations
        """
        ...
    
    def getDataProcessing(self) -> List[DataProcessing]:
        """
        Cython signature: libcpp_vector[shared_ptr[DataProcessing]] getDataProcessing()
        Returns a reference to the description of the applied processing
        """
        ...
    
    def setDataProcessing(self, in_0: List[DataProcessing] ) -> None:
        """
        Cython signature: void setDataProcessing(libcpp_vector[shared_ptr[DataProcessing]])
        Sets the description of the applied processing
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
    
    def __richcmp__(self, other: IntegerDataArray, op: int) -> Any:
        ... 


class MRMScoring:
    """
    Cython implementation of _MRMScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1MRMScoring.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MRMScoring()
        """
        ...
    
    @overload
    def __init__(self, in_0: MRMScoring ) -> None:
        """
        Cython signature: void MRMScoring(MRMScoring &)
        """
        ...
    
    def calcXcorrCoelutionScore(self) -> float:
        """
        Cython signature: double calcXcorrCoelutionScore()
        Calculate the cross-correlation coelution score. The score is a distance where zero indicates perfect coelution
        """
        ...
    
    def calcXcorrCoelutionWeightedScore(self, normalized_library_intensity: List[float] ) -> float:
        """
        Cython signature: double calcXcorrCoelutionWeightedScore(libcpp_vector[double] & normalized_library_intensity)
        Calculate the weighted cross-correlation coelution score
        
        The score is a distance where zero indicates perfect coelution. The
        score is weighted by the transition intensities, non-perfect coelution
        in low-intensity transitions should thus become less important
        """
        ...
    
    def calcSeparateXcorrContrastCoelutionScore(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] calcSeparateXcorrContrastCoelutionScore()
        Calculate the separate cross-correlation contrast score
        """
        ...
    
    def calcXcorrPrecursorContrastCoelutionScore(self) -> float:
        """
        Cython signature: double calcXcorrPrecursorContrastCoelutionScore()
        Calculate the precursor cross-correlation contrast score against the transitions
        
        The score is a distance where zero indicates perfect coelution
        """
        ...
    
    def calcXcorrShapeScore(self) -> float:
        """
        Cython signature: double calcXcorrShapeScore()
        Calculate the cross-correlation shape score
        
        The score is a correlation measure where 1 indicates perfect correlation
        and 0 means no correlation.
        """
        ...
    
    def calcXcorrShapeWeightedScore(self, normalized_library_intensity: List[float] ) -> float:
        """
        Cython signature: double calcXcorrShapeWeightedScore(libcpp_vector[double] & normalized_library_intensity)
        Calculate the weighted cross-correlation shape score
        
        The score is a correlation measure where 1 indicates perfect correlation
        and 0 means no correlation. The score is weighted by the transition
        intensities, non-perfect coelution in low-intensity transitions should
        thus become less important
        """
        ...
    
    def calcSeparateXcorrContrastShapeScore(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] calcSeparateXcorrContrastShapeScore()
        Calculate the separate cross-correlation contrast shape score
        """
        ...
    
    def calcXcorrPrecursorContrastShapeScore(self) -> float:
        """
        Cython signature: double calcXcorrPrecursorContrastShapeScore()
        Calculate the precursor cross-correlation shape score against the transitions
        """
        ...
    
    def calcRTScore(self, peptide: LightCompound , normalized_experimental_rt: float ) -> float:
        """
        Cython signature: double calcRTScore(LightCompound & peptide, double normalized_experimental_rt)
        """
        ...
    
    def calcMIScore(self) -> float:
        """
        Cython signature: double calcMIScore()
        """
        ...
    
    def calcMIWeightedScore(self, normalized_library_intensity: List[float] ) -> float:
        """
        Cython signature: double calcMIWeightedScore(const libcpp_vector[double] & normalized_library_intensity)
        """
        ...
    
    def calcMIPrecursorScore(self) -> float:
        """
        Cython signature: double calcMIPrecursorScore()
        """
        ...
    
    def calcMIPrecursorContrastScore(self) -> float:
        """
        Cython signature: double calcMIPrecursorContrastScore()
        """
        ...
    
    def calcMIPrecursorCombinedScore(self) -> float:
        """
        Cython signature: double calcMIPrecursorCombinedScore()
        """
        ...
    
    def calcSeparateMIContrastScore(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] calcSeparateMIContrastScore()
        """
        ...
    
    def getMIMatrix(self) -> MatrixDouble:
        """
        Cython signature: MatrixDouble getMIMatrix()
        """
        ... 


class MapAlignmentTransformer:
    """
    Cython implementation of _MapAlignmentTransformer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentTransformer.html>`_

    This class collects functions for applying retention time transformations to data structures
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MapAlignmentTransformer()
        """
        ...
    
    @overload
    def __init__(self, in_0: MapAlignmentTransformer ) -> None:
        """
        Cython signature: void MapAlignmentTransformer(MapAlignmentTransformer &)
        """
        ...
    
    @overload
    def transformRetentionTimes(self, in_0: MSExperiment , in_1: TransformationDescription , in_2: bool ) -> None:
        """
        Cython signature: void transformRetentionTimes(MSExperiment &, TransformationDescription &, bool)
        Applies the given transformation to a peak map
        """
        ...
    
    @overload
    def transformRetentionTimes(self, in_0: FeatureMap , in_1: TransformationDescription , in_2: bool ) -> None:
        """
        Cython signature: void transformRetentionTimes(FeatureMap &, TransformationDescription &, bool)
        Applies the given transformation to a feature map
        """
        ...
    
    @overload
    def transformRetentionTimes(self, in_0: ConsensusMap , in_1: TransformationDescription , in_2: bool ) -> None:
        """
        Cython signature: void transformRetentionTimes(ConsensusMap &, TransformationDescription &, bool)
        Applies the given transformation to a consensus map
        """
        ...
    
    @overload
    def transformRetentionTimes(self, in_0: PeptideIdentificationList , in_1: TransformationDescription , in_2: bool ) -> None:
        """
        Cython signature: void transformRetentionTimes(PeptideIdentificationList &, TransformationDescription &, bool)
        Applies the given transformation to peptide identifications
        """
        ... 


class OMSSAXMLFile:
    """
    Cython implementation of _OMSSAXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OMSSAXMLFile.html>`_
      -- Inherits from ['XMLFile']

    Used to load OMSSAXML files
    
    This class is used to load documents that implement
    the schema of OMSSAXML files
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void OMSSAXMLFile()
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , load_proteins: bool , load_empty_hits: bool ) -> None:
        """
        Cython signature: void load(const String & filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data, bool load_proteins, bool load_empty_hits)
        Loads data from a OMSSAXML file
        
        
        :param filename: The file to be loaded
        :param protein_identification: Protein identifications belonging to the whole experiment
        :param id_data: The identifications with m/z and RT
        :param load_proteins: If this flag is set to false, the protein identifications are not loaded
        :param load_empty_hits: Many spectra will not return a hit. Report empty peptide identifications?
        """
        ...
    
    def setModificationDefinitionsSet(self, rhs: ModificationDefinitionsSet ) -> None:
        """
        Cython signature: void setModificationDefinitionsSet(ModificationDefinitionsSet rhs)
        Sets the valid modifications
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
        """
        ... 


class OSChromatogramMeta:
    """
    Cython implementation of _OSChromatogramMeta

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSChromatogramMeta.html>`_
    """
    
    index: int
    
    id: bytes
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OSChromatogramMeta()
        """
        ...
    
    @overload
    def __init__(self, in_0: OSChromatogramMeta ) -> None:
        """
        Cython signature: void OSChromatogramMeta(OSChromatogramMeta &)
        """
        ... 


class PeakBoundary:
    """
    Cython implementation of _PeakBoundary

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakBoundary.html>`_
    """
    
    mz_min: float
    
    mz_max: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakBoundary()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakBoundary ) -> None:
        """
        Cython signature: void PeakBoundary(PeakBoundary &)
        """
        ... 


class PeakPickerChromatogram:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakPickerChromatogram()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakPickerChromatogram ) -> None:
        """
        Cython signature: void PeakPickerChromatogram(PeakPickerChromatogram &)
        """
        ...
    
    def pickChromatogram(self, chromatogram: MSChromatogram , picked_chrom: MSChromatogram ) -> None:
        """
        Cython signature: void pickChromatogram(MSChromatogram & chromatogram, MSChromatogram & picked_chrom)
        Finds peaks in a single chromatogram and annotates left/right borders
        
        It uses a modified algorithm of the PeakPickerHiRes
        
        This function will return a picked chromatogram
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


class PeakPickerHiRes:
    """
    Cython implementation of _PeakPickerHiRes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerHiRes.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakPickerHiRes()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakPickerHiRes ) -> None:
        """
        Cython signature: void PeakPickerHiRes(PeakPickerHiRes &)
        """
        ...
    
    @overload
    def pick(self, input: MSSpectrum , output: MSSpectrum ) -> None:
        """
        Cython signature: void pick(MSSpectrum & input, MSSpectrum & output)
        """
        ...
    
    @overload
    def pick(self, input: MSChromatogram , output: MSChromatogram ) -> None:
        """
        Cython signature: void pick(MSChromatogram & input, MSChromatogram & output)
        """
        ...
    
    @overload
    def pickExperiment(self, input: MSExperiment , output: MSExperiment , check_spectrum_type: bool ) -> None:
        """
        Cython signature: void pickExperiment(MSExperiment & input, MSExperiment & output, bool check_spectrum_type)
        Applies the peak-picking algorithm to a map (MSExperiment). This method picks peaks for each scan in the map consecutively. The resulting
        picked peaks are written to the output map
        
        
        :param input: Input map in profile mode
        :param output: Output map with picked peaks
        :param check_spectrum_type: If set, checks spectrum type and throws an exception if a centroided spectrum is passed
        """
        ...
    
    @overload
    def pickExperiment(self, input: MSExperiment , output: MSExperiment , boundaries_spec: List[List[PeakBoundary]] , boundaries_chrom: List[List[PeakBoundary]] , check_spectrum_type: bool ) -> None:
        """
        Cython signature: void pickExperiment(MSExperiment & input, MSExperiment & output, libcpp_vector[libcpp_vector[PeakBoundary]] & boundaries_spec, libcpp_vector[libcpp_vector[PeakBoundary]] & boundaries_chrom, bool check_spectrum_type)
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


class String:
    """
    Cython implementation of _String

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1String.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void String()
        """
        ...
    
    def __richcmp__(self, other: String, op: int) -> Any:
        ... 


class XFDRAlgorithm:
    """
    Cython implementation of _XFDRAlgorithm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XFDRAlgorithm.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void XFDRAlgorithm()
        """
        ...
    
    @overload
    def __init__(self, in_0: XFDRAlgorithm ) -> None:
        """
        Cython signature: void XFDRAlgorithm(XFDRAlgorithm &)
        """
        ...
    
    def run(self, peptide_ids: PeptideIdentificationList , protein_id: ProteinIdentification ) -> int:
        """
        Cython signature: XFDRAlgorithm_ExitCodes run(PeptideIdentificationList & peptide_ids, ProteinIdentification & protein_id)
        """
        ...
    
    def validateClassArguments(self) -> int:
        """
        Cython signature: XFDRAlgorithm_ExitCodes validateClassArguments()
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
    XFDRAlgorithm_ExitCodes : __XFDRAlgorithm_ExitCodes 


class XQuestResultXMLFile:
    """
    Cython implementation of _XQuestResultXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XQuestResultXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void XQuestResultXMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: XQuestResultXMLFile ) -> None:
        """
        Cython signature: void XQuestResultXMLFile(XQuestResultXMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , pep_ids: PeptideIdentificationList , prot_ids: List[ProteinIdentification] ) -> None:
        """
        Cython signature: void load(const String & filename, PeptideIdentificationList & pep_ids, libcpp_vector[ProteinIdentification] & prot_ids)
        Load the content of the xquest.xml file into the provided data structures
        
        :param filename: Filename of the file which is to be loaded
        :param pep_ids: Where the spectra with identifications of the input file will be loaded to
        :param prot_ids: Where the protein identification of the input file will be loaded to
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , poid: List[ProteinIdentification] , peid: PeptideIdentificationList ) -> None:
        """
        Cython signature: void store(const String & filename, libcpp_vector[ProteinIdentification] & poid, PeptideIdentificationList & peid)
        Stores the identifications in a xQuest XML file
        """
        ...
    
    def getNumberOfHits(self) -> int:
        """
        Cython signature: int getNumberOfHits()
        Returns the total number of hits in the file
        """
        ...
    
    def getMinScore(self) -> float:
        """
        Cython signature: double getMinScore()
        Returns minimum score among the hits in the file
        """
        ...
    
    def getMaxScore(self) -> float:
        """
        Cython signature: double getMaxScore()
        Returns maximum score among the hits in the file
        """
        ...
    
    @overload
    def writeXQuestXMLSpec(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , preprocessed_pair_spectra: OPXL_PreprocessedPairSpectra , spectrum_pairs: List[List[int, int]] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None:
        """
        Cython signature: void writeXQuestXMLSpec(const String & out_file, const String & base_name, OPXL_PreprocessedPairSpectra preprocessed_pair_spectra, libcpp_vector[libcpp_pair[size_t,size_t]] spectrum_pairs, libcpp_vector[libcpp_vector[CrossLinkSpectrumMatch]] all_top_csms, MSExperiment spectra, const bool & test_mode)
        Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param preprocessed_pair_spectra: The preprocessed spectra after comparing and filtering
        :param spectrum_pairs: Indices of spectrum pairs in the input map
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map
        """
        ...
    
    @overload
    def writeXQuestXMLSpec(self, out_file: Union[bytes, str, String] , base_name: Union[bytes, str, String] , all_top_csms: List[List[CrossLinkSpectrumMatch]] , spectra: MSExperiment , test_mode: bool ) -> None:
        """
        Cython signature: void writeXQuestXMLSpec(const String & out_file, const String & base_name, libcpp_vector[libcpp_vector[CrossLinkSpectrumMatch]] all_top_csms, MSExperiment spectra, const bool & test_mode)
        Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers
        
        :param out_file: Path and filename for the output file
        :param base_name: The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra
        :param all_top_csms: CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out
        :param spectra: The spectra, that were searched as a PeakMap
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
        """
        ... 


class XQuestScores:
    """
    Cython implementation of _XQuestScores

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1XQuestScores.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void XQuestScores()
        """
        ...
    
    @overload
    def __init__(self, in_0: XQuestScores ) -> None:
        """
        Cython signature: void XQuestScores(XQuestScores &)
        """
        ...
    
    @overload
    def preScore(self, matched_alpha: int , ions_alpha: int , matched_beta: int , ions_beta: int ) -> float:
        """
        Cython signature: float preScore(size_t matched_alpha, size_t ions_alpha, size_t matched_beta, size_t ions_beta)
        Compute a simple and fast to compute pre-score for a cross-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
        :param matched_beta: Number of experimental peaks matched to theoretical linear ions from the beta peptide
        :param ions_beta: Number of theoretical ions from the beta peptide
        """
        ...
    
    @overload
    def preScore(self, matched_alpha: int , ions_alpha: int ) -> float:
        """
        Cython signature: float preScore(size_t matched_alpha, size_t ions_alpha)
        Compute a simple and fast to compute pre-score for a mono-link spectrum match
        
        :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
        :param ions_alpha: Number of theoretical ions from the alpha peptide
        """
        ...
    
    def matchOddsScore(self, theoretical_spec: MSSpectrum , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , is_xlink_spectrum: bool , n_charges: int ) -> float:
        """
        Cython signature: double matchOddsScore(MSSpectrum & theoretical_spec, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum, size_t n_charges)
        Compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
        
        :param theoretical_spec: Theoretical spectrum, sorted by position
        :param matched_size: Alignment between the theoretical and the experimental spectra
        :param fragment_mass_tolerance: Fragment mass tolerance of the alignment
        :param fragment_mass_tolerance_unit_ppm: Fragment mass tolerance unit of the alignment, true = ppm, false = Da
        :param is_xlink_spectrum: Type of cross-link, true = cross-link, false = mono-link
        :param n_charges: Number of considered charges in the theoretical spectrum
        """
        ...
    
    def logOccupancyProb(self, theoretical_spec: MSSpectrum , matched_size: int , fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool ) -> float:
        """
        Cython signature: double logOccupancyProb(MSSpectrum theoretical_spec, size_t matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
        Compute the logOccupancyProb score, similar to the match_odds, a score based on the probability of getting the given number of matched peaks by chance
        
        :param theoretical_spec: Theoretical spectrum, sorted by position
        :param matched_size: Number of matched peaks between experimental and theoretical spectra
        :param fragment_mass_tolerance: The tolerance of the alignment
        :param fragment_mass_tolerance_unit: The tolerance unit of the alignment, true = ppm, false = Da
        """
        ...
    
    def weightedTICScoreXQuest(self, alpha_size: int , beta_size: int , intsum_alpha: float , intsum_beta: float , total_current: float , type_is_cross_link: bool ) -> float:
        """
        Cython signature: double weightedTICScoreXQuest(size_t alpha_size, size_t beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
        """
        ...
    
    def weightedTICScore(self, alpha_size: int , beta_size: int , intsum_alpha: float , intsum_beta: float , total_current: float , type_is_cross_link: bool ) -> float:
        """
        Cython signature: double weightedTICScore(size_t alpha_size, size_t beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
        """
        ...
    
    def matchedCurrentChain(self, matched_spec_common: List[List[int, int]] , matched_spec_xlinks: List[List[int, int]] , spectrum_common_peaks: MSSpectrum , spectrum_xlink_peaks: MSSpectrum ) -> float:
        """
        Cython signature: double matchedCurrentChain(libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_common, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_xlinks, MSSpectrum & spectrum_common_peaks, MSSpectrum & spectrum_xlink_peaks)
        """
        ...
    
    def totalMatchedCurrent(self, matched_spec_common_alpha: List[List[int, int]] , matched_spec_common_beta: List[List[int, int]] , matched_spec_xlinks_alpha: List[List[int, int]] , matched_spec_xlinks_beta: List[List[int, int]] , spectrum_common_peaks: MSSpectrum , spectrum_xlink_peaks: MSSpectrum ) -> float:
        """
        Cython signature: double totalMatchedCurrent(libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_common_alpha, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_common_beta, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_xlinks_alpha, libcpp_vector[libcpp_pair[size_t,size_t]] & matched_spec_xlinks_beta, MSSpectrum & spectrum_common_peaks, MSSpectrum & spectrum_xlink_peaks)
        """
        ...
    
    def xCorrelation(self, spec1: MSSpectrum , spec2: MSSpectrum , maxshift: int , tolerance: float ) -> List[float]:
        """
        Cython signature: libcpp_vector[double] xCorrelation(MSSpectrum & spec1, MSSpectrum & spec2, int maxshift, double tolerance)
        """
        ...
    
    def xCorrelationPrescore(self, spec1: MSSpectrum , spec2: MSSpectrum , tolerance: float ) -> float:
        """
        Cython signature: double xCorrelationPrescore(MSSpectrum & spec1, MSSpectrum & spec2, double tolerance)
        """
        ... 


class __CombinationsLogic:
    None
    OR : int
    AND : int
    XOR : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class QuotingMethod:
    None
    NONE : int
    ESCAPE : int
    DOUBLE : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __RequirementLevel:
    None
    MUST : int
    SHOULD : int
    MAY : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __XFDRAlgorithm_ExitCodes:
    None
    EXECUTION_OK : int
    ILLEGAL_PARAMETERS : int
    UNEXPECTED_RESULT : int

    def getMapping(self) -> Dict[int, str]:
       ... 

