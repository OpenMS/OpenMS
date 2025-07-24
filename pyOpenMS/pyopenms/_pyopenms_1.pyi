from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AScore:
    """
    Cython implementation of _AScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AScore.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AScore()
        """
        ...
    
    @overload
    def __init__(self, in_0: AScore ) -> None:
        """
        Cython signature: void AScore(AScore &)
        """
        ...
    
    def compute(self, hit: PeptideHit , real_spectrum: MSSpectrum ) -> PeptideHit:
        """
        Cython signature: PeptideHit compute(PeptideHit & hit, MSSpectrum & real_spectrum)
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


class CVTermListInterface:
    """
    Cython implementation of _CVTermListInterface

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVTermListInterface.html>`_
      -- Inherits from ['MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVTermListInterface()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVTermListInterface ) -> None:
        """
        Cython signature: void CVTermListInterface(CVTermListInterface &)
        """
        ...
    
    @overload
    def replaceCVTerms(self, cv_terms: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] & cv_terms)
        """
        ...
    
    @overload
    def replaceCVTerms(self, cv_terms: List[CVTerm] , accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void replaceCVTerms(libcpp_vector[CVTerm] & cv_terms, const String & accession)
        """
        ...
    
    def setCVTerms(self, terms: List[CVTerm] ) -> None:
        """
        Cython signature: void setCVTerms(libcpp_vector[CVTerm] & terms)
        """
        ...
    
    def replaceCVTerm(self, cv_term: CVTerm ) -> None:
        """
        Cython signature: void replaceCVTerm(CVTerm & cv_term)
        """
        ...
    
    def consumeCVTerms(self, cv_term_map: Dict[bytes,List[CVTerm]] ) -> None:
        """
        Cython signature: void consumeCVTerms(libcpp_map[String,libcpp_vector[CVTerm]] & cv_term_map)
        Merges the given map into the member map, no duplicate checking
        """
        ...
    
    def getCVTerms(self) -> Dict[bytes,List[CVTerm]]:
        """
        Cython signature: libcpp_map[String,libcpp_vector[CVTerm]] getCVTerms()
        """
        ...
    
    def addCVTerm(self, term: CVTerm ) -> None:
        """
        Cython signature: void addCVTerm(CVTerm & term)
        Adds a CV term
        """
        ...
    
    def hasCVTerm(self, accession: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasCVTerm(const String & accession)
        Checks whether the term has a value
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
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
    
    def __richcmp__(self, other: CVTermListInterface, op: int) -> Any:
        ... 


class ChromatogramTools:
    """
    Cython implementation of _ChromatogramTools

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ChromatogramTools.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ChromatogramTools()
        """
        ...
    
    @overload
    def __init__(self, in_0: ChromatogramTools ) -> None:
        """
        Cython signature: void ChromatogramTools(ChromatogramTools &)
        """
        ...
    
    def convertChromatogramsToSpectra(self, epx: MSExperiment ) -> None:
        """
        Cython signature: void convertChromatogramsToSpectra(MSExperiment & epx)
        Converts the chromatogram to a list of spectra with instrument settings
        """
        ...
    
    def convertSpectraToChromatograms(self, epx: MSExperiment , remove_spectra: bool , force_conversion: bool ) -> None:
        """
        Cython signature: void convertSpectraToChromatograms(MSExperiment & epx, bool remove_spectra, bool force_conversion)
        Converts e.g. SRM spectra to chromatograms
        """
        ... 


class EDTAFile:
    """
    Cython implementation of _EDTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1EDTAFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void EDTAFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: EDTAFile ) -> None:
        """
        Cython signature: void EDTAFile(EDTAFile &)
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , map: FeatureMap ) -> None:
        """
        Cython signature: void store(String filename, FeatureMap & map)
        """
        ...
    
    @overload
    def store(self, filename: Union[bytes, str, String] , map: ConsensusMap ) -> None:
        """
        Cython signature: void store(String filename, ConsensusMap & map)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , consensus_map: ConsensusMap ) -> None:
        """
        Cython signature: void load(String filename, ConsensusMap & consensus_map)
        """
        ... 


class FileTypes:
    """
    Cython implementation of _FileTypes

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FileTypes.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FileTypes()
        Centralizes the file types recognized by FileHandler
        """
        ...
    
    @overload
    def __init__(self, in_0: FileTypes ) -> None:
        """
        Cython signature: void FileTypes(FileTypes &)
        """
        ...
    
    def typeToName(self, t: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String typeToName(FileType t)
        Returns the name/extension of the type
        """
        ...
    
    def typeToMZML(self, t: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String typeToMZML(FileType t)
        Returns the mzML name
        """
        ...
    
    def nameToType(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: FileType nameToType(String name)
        Converts a file type name into a Type
        
        
        :param name: A case-insensitive name (e.g. FASTA or Fasta, etc.)
        """
        ... 


class IDDecoyProbability:
    """
    Cython implementation of _IDDecoyProbability

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1IDDecoyProbability.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IDDecoyProbability()
        IDDecoyProbability calculates probabilities using decoy approach
        """
        ...
    
    @overload
    def __init__(self, in_0: IDDecoyProbability ) -> None:
        """
        Cython signature: void IDDecoyProbability(IDDecoyProbability)
        """
        ...
    
    @overload
    def apply(self, prob_ids: PeptideIdentificationList , fwd_ids: PeptideIdentificationList , rev_ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void apply(PeptideIdentificationList & prob_ids, PeptideIdentificationList & fwd_ids, PeptideIdentificationList & rev_ids)
        Converts the forward and reverse identification into probabilities
        
        
        :param prob_ids: Output of the algorithm which includes identifications with probability based scores
        :param fwd_ids: Input parameter which represents the identifications of the forward search
        :param rev_ids: Input parameter which represents the identifications of the reversed search
        """
        ...
    
    @overload
    def apply(self, ids: PeptideIdentificationList ) -> None:
        """
        Cython signature: void apply(PeptideIdentificationList & ids)
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


class KroenikFile:
    """
    Cython implementation of _KroenikFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1KroenikFile.html>`_

    File adapter for Kroenik (HardKloer sibling) files
    
    The first line is the header and contains the column names:
    File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
    
    Every subsequent line is a feature
    
    All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
    metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications"
    
    The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void KroenikFile()
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void store(String filename, MSSpectrum & spectrum)
        Stores a MSExperiment into a Kroenik file
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , feature_map: FeatureMap ) -> None:
        """
        Cython signature: void load(String filename, FeatureMap & feature_map)
        Loads a Kroenik file into a featureXML
        
        The content of the file is stored in `features`
        
        :raises:
          Exception: FileNotFound is thrown if the file could not be opened
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ... 


class MSPGenericFile:
    """
    Cython implementation of _MSPGenericFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSPGenericFile.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSPGenericFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MSPGenericFile ) -> None:
        """
        Cython signature: void MSPGenericFile(MSPGenericFile &)
        """
        ...
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None:
        """
        Cython signature: void MSPGenericFile(const String & filename, MSExperiment & library)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None:
        """
        Cython signature: void load(const String & filename, MSExperiment & library)
        Load the file's data and metadata, and save it into an `MSExperiment`
        
        
        :param filename: Path to the MSP input file
        :param library: The variable into which the extracted information will be saved
        :raises:
          Exception: FileNotFound If the file could not be found
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , library: MSExperiment ) -> None:
        """
        Cython signature: void store(const String & filename, const MSExperiment & library)
        Save data and metadata into a file
        
        
        :param filename: Path to the MSP input file
        :param library: The variable from which extracted information will be saved
        :raises:
          Exception: FileNotWritable If the file is not writable
        """
        ...
    
    def getDefaultParameters(self, params: Param ) -> None:
        """
        Cython signature: void getDefaultParameters(Param & params)
        Returns the class' default parameters
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


class PeakPickerIterative:
    """
    Cython implementation of _PeakPickerIterative

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerIterative.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeakPickerIterative()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeakPickerIterative ) -> None:
        """
        Cython signature: void PeakPickerIterative(PeakPickerIterative &)
        """
        ...
    
    def pick(self, input: MSSpectrum , output: MSSpectrum ) -> None:
        """
        Cython signature: void pick(MSSpectrum & input, MSSpectrum & output)
        This will pick one single spectrum. The PeakPickerHiRes is used to
        generate seeds, these seeds are then used to re-center the mass and
        compute peak width and integrated intensity of the peak
        
        Finally, other peaks that would fall within the primary peak are
        discarded
        
        The output are the remaining peaks
        """
        ...
    
    def pickExperiment(self, input: MSExperiment , output: MSExperiment ) -> None:
        """
        Cython signature: void pickExperiment(MSExperiment & input, MSExperiment & output)
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


class PeptideIdentificationList:
    """
    Cython implementation of _PeptideIdentificationList

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideIdentificationList.html>`_

    A container for peptide identifications from multiple spectra.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideIdentificationList()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideIdentificationList ) -> None:
        """
        Cython signature: void PeptideIdentificationList(PeptideIdentificationList)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        Returns the number of peptide identifications
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Returns true if the container is empty
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Removes all peptide identifications from the container
        """
        ...
    
    def push_back(self, in_0: PeptideIdentification ) -> None:
        """
        Cython signature: void push_back(PeptideIdentification)
        Adds a peptide identification to the end of the container
        """
        ...
    
    def __getitem__(self, in_0: int ) -> PeptideIdentification:
        """
        Cython signature: PeptideIdentification & operator[](size_t)
        """
        ...
    def __setitem__(self, key: int, value: PeptideIdentification ) -> None:
        """Cython signature: PeptideIdentification & operator[](size_t)"""
        ...
    
    def at(self, in_0: int ) -> PeptideIdentification:
        """
        Cython signature: PeptideIdentification at(size_t)
        Returns the peptide identification at the given index with bounds checking
        """
        ...
    
    def back(self) -> PeptideIdentification:
        """
        Cython signature: PeptideIdentification back()
        Returns the last peptide identification in the container
        """
        ...
    
    def front(self) -> PeptideIdentification:
        """
        Cython signature: PeptideIdentification front()
        Returns the first peptide identification in the container
        """
        ...
    
    def __iter__(self) -> PeptideIdentification:
       ... 


class ProbablePhosphoSites:
    """
    Cython implementation of _ProbablePhosphoSites

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ProbablePhosphoSites.html>`_
    """
    
    first: int
    
    second: int
    
    seq_1: int
    
    seq_2: int
    
    peak_depth: int
    
    AScore: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProbablePhosphoSites()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProbablePhosphoSites ) -> None:
        """
        Cython signature: void ProbablePhosphoSites(ProbablePhosphoSites &)
        """
        ... 


class RealMassDecomposer:
    """
    Cython implementation of _RealMassDecomposer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims_1_1RealMassDecomposer.html>`_
    """
    
    @overload
    def __init__(self, in_0: RealMassDecomposer ) -> None:
        """
        Cython signature: void RealMassDecomposer(RealMassDecomposer)
        """
        ...
    
    @overload
    def __init__(self, weights: IMSWeights ) -> None:
        """
        Cython signature: void RealMassDecomposer(IMSWeights & weights)
        """
        ...
    
    def getNumberOfDecompositions(self, mass: float , error: float ) -> int:
        """
        Cython signature: uint64_t getNumberOfDecompositions(double mass, double error)
        Gets a number of all decompositions for amass with an error
        allowed. It's similar to thegetDecompositions(double,double) function
        but less space consuming, since doesn't use container to store decompositions
        
        
        :param mass: Mass to be decomposed
        :param error: Error allowed between given and result decomposition
        :return: Number of all decompositions for a given mass and error
        """
        ... 


class SemanticValidator:
    """
    Cython implementation of _SemanticValidator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1SemanticValidator.html>`_
    """
    
    def __init__(self, mapping: CVMappings , cv: ControlledVocabulary ) -> None:
        """
        Cython signature: void SemanticValidator(CVMappings mapping, ControlledVocabulary cv)
        """
        ...
    
    def validate(self, filename: Union[bytes, str, String] , errors: List[bytes] , warnings: List[bytes] ) -> bool:
        """
        Cython signature: bool validate(String filename, StringList errors, StringList warnings)
        """
        ...
    
    def locateTerm(self, path: Union[bytes, str, String] , parsed_term: SemanticValidator_CVTerm ) -> bool:
        """
        Cython signature: bool locateTerm(String path, SemanticValidator_CVTerm & parsed_term)
        Checks if a CVTerm is allowed in a given path
        """
        ...
    
    def setTag(self, tag: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setTag(String tag)
        Sets the CV parameter tag name (default 'cvParam')
        """
        ...
    
    def setAccessionAttribute(self, accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setAccessionAttribute(String accession)
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'accession')
        """
        ...
    
    def setNameAttribute(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setNameAttribute(String name)
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'name')
        """
        ...
    
    def setValueAttribute(self, value: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setValueAttribute(String value)
        Sets the name of the attribute for accessions in the CV parameter tag name (default 'value')
        """
        ...
    
    def setCheckTermValueTypes(self, check: bool ) -> None:
        """
        Cython signature: void setCheckTermValueTypes(bool check)
        Sets if CV term value types should be check (enabled by default)
        """
        ...
    
    def setCheckUnits(self, check: bool ) -> None:
        """
        Cython signature: void setCheckUnits(bool check)
        Sets if CV term units should be check (disabled by default)
        """
        ...
    
    def setUnitAccessionAttribute(self, accession: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setUnitAccessionAttribute(String accession)
        Sets the name of the unit accession attribute (default 'unitAccession')
        """
        ...
    
    def setUnitNameAttribute(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setUnitNameAttribute(String name)
        Sets the name of the unit name attribute (default 'unitName')
        """
        ... 


class SemanticValidator_CVTerm:
    """
    Cython implementation of _SemanticValidator_CVTerm

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Internal_1_1SemanticValidator_CVTerm.html>`_
    """
    
    accession: Union[bytes, str, String]
    
    name: Union[bytes, str, String]
    
    value: Union[bytes, str, String]
    
    has_value: bool
    
    unit_accession: Union[bytes, str, String]
    
    has_unit_accession: bool
    
    unit_name: Union[bytes, str, String]
    
    has_unit_name: bool
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SemanticValidator_CVTerm()
        """
        ...
    
    @overload
    def __init__(self, in_0: SemanticValidator_CVTerm ) -> None:
        """
        Cython signature: void SemanticValidator_CVTerm(SemanticValidator_CVTerm &)
        """
        ... 


class SpectrumAlignment:
    """
    Cython implementation of _SpectrumAlignment

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAlignment.html>`_
      -- Inherits from ['DefaultParamHandler']

    Aligns the peaks of two sorted spectra
    Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given
        Scoring function is the m/z distance between peaks. Intensity does not play a role!
    Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
    Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window
    
    note: A peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2
    note: Intensity is ignored
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAlignment()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAlignment ) -> None:
        """
        Cython signature: void SpectrumAlignment(SpectrumAlignment &)
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


class SpectrumAnnotator:
    """
    Cython implementation of _SpectrumAnnotator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAnnotator.html>`_
      -- Inherits from ['DefaultParamHandler']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAnnotator()
        Annotates spectra from identifications and theoretical spectra or
        identifications from spectra and theoretical spectra matching
        with various options
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAnnotator ) -> None:
        """
        Cython signature: void SpectrumAnnotator(SpectrumAnnotator &)
        """
        ...
    
    def annotateMatches(self, spec: MSSpectrum , ph: PeptideHit , tg: TheoreticalSpectrumGenerator , sa: SpectrumAlignment ) -> None:
        """
        Cython signature: void annotateMatches(MSSpectrum & spec, PeptideHit & ph, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa)
        Adds ion match annotation to the `spec` input spectrum
        
        :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
        :param ph: A spectrum identifications to be used for the annotation, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
        :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
        :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
        """
        ...
    
    def addIonMatchStatistics(self, pi: PeptideIdentification , spec: MSSpectrum , tg: TheoreticalSpectrumGenerator , sa: SpectrumAlignment ) -> None:
        """
        Cython signature: void addIonMatchStatistics(PeptideIdentification & pi, MSSpectrum & spec, TheoreticalSpectrumGenerator & tg, SpectrumAlignment & sa)
        Adds ion match statistics to `pi` PeptideIdentifcation
        
        :param pi: A spectrum identifications to be annotated, looking up matches from a spectrum and the theoretical spectrum inferred from the identifications sequence
        :param spec: A PeakSpectrum containing the peaks from which the `pi` identifications are made
        :param tg: A TheoreticalSpectrumGenerator to infer the theoretical spectrum. Its own parameters define which ion types are referred
        :param sa: A SpectrumAlignment to match the theoretical spectrum with the measured. Its own parameters define the match tolerance
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


class FileType:
    None
    UNKNOWN : int
    DTA : int
    DTA2D : int
    MZDATA : int
    MZXML : int
    FEATUREXML : int
    IDXML : int
    CONSENSUSXML : int
    MGF : int
    INI : int
    TOPPAS : int
    TRANSFORMATIONXML : int
    MZML : int
    CACHEDMZML : int
    MS2 : int
    PEPXML : int
    PROTXML : int
    MZIDENTML : int
    QCML : int
    GELML : int
    TRAML : int
    MSP : int
    OMSSAXML : int
    MASCOTXML : int
    PNG : int
    XMASS : int
    TSV : int
    PEPLIST : int
    HARDKLOER : int
    KROENIK : int
    FASTA : int
    EDTA : int
    CSV : int
    TXT : int
    OBO : int
    HTML : int
    XML : int
    ANALYSISXML : int
    XSD : int
    PSQ : int
    MRM : int
    SQMASS : int
    PQP : int
    OSW : int
    PSMS : int
    PARAMXML : int
    SIZE_OF_TYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

