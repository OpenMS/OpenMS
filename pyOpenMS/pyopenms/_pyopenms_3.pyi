from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum




class AccurateMassSearchEngine:
    """
    Cython implementation of _AccurateMassSearchEngine

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1AccurateMassSearchEngine.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void AccurateMassSearchEngine()
        """
        ...
    
    @overload
    def __init__(self, in_0: AccurateMassSearchEngine ) -> None:
        """
        Cython signature: void AccurateMassSearchEngine(AccurateMassSearchEngine &)
        """
        ...
    
    def queryByMZ(self, observed_mz: float , observed_charge: int , ion_mode: Union[bytes, str, String] , in_3: List[AccurateMassSearchResult] , observed_adduct: EmpiricalFormula ) -> None:
        """
        Cython signature: void queryByMZ(double observed_mz, int observed_charge, String ion_mode, libcpp_vector[AccurateMassSearchResult] &, EmpiricalFormula & observed_adduct)
        """
        ...
    
    def queryByFeature(self, feature: Feature , feature_index: int , ion_mode: Union[bytes, str, String] , in_3: List[AccurateMassSearchResult] ) -> None:
        """
        Cython signature: void queryByFeature(Feature feature, size_t feature_index, String ion_mode, libcpp_vector[AccurateMassSearchResult] &)
        """
        ...
    
    def queryByConsensusFeature(self, cfeat: ConsensusFeature , cf_index: int , number_of_maps: int , ion_mode: Union[bytes, str, String] , results: List[AccurateMassSearchResult] ) -> None:
        """
        Cython signature: void queryByConsensusFeature(ConsensusFeature cfeat, size_t cf_index, size_t number_of_maps, String ion_mode, libcpp_vector[AccurateMassSearchResult] & results)
        """
        ...
    
    @overload
    def run(self, in_0: FeatureMap , in_1: MzTab ) -> None:
        """
        Cython signature: void run(FeatureMap &, MzTab &)
        """
        ...
    
    @overload
    def run(self, in_0: FeatureMap , in_1: MzTabM ) -> None:
        """
        Cython signature: void run(FeatureMap &, MzTabM &)
        """
        ...
    
    @overload
    def run(self, in_0: ConsensusMap , in_1: MzTab ) -> None:
        """
        Cython signature: void run(ConsensusMap &, MzTab &)
        """
        ...
    
    def init(self) -> None:
        """
        Cython signature: void init()
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


class CVMappingFile:
    """
    Cython implementation of _CVMappingFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVMappingFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void CVMappingFile()
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , cv_mappings: CVMappings , strip_namespaces: bool ) -> None:
        """
        Cython signature: void load(const String & filename, CVMappings & cv_mappings, bool strip_namespaces)
        Loads CvMappings from the given file
        """
        ... 


class FASTAEntry:
    """
    Cython implementation of _FASTAEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FASTAEntry.html>`_
    """
    
    identifier: Union[bytes, str, String]
    
    description: Union[bytes, str, String]
    
    sequence: Union[bytes, str, String]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FASTAEntry()
        """
        ...
    
    @overload
    def __init__(self, in_0: FASTAEntry ) -> None:
        """
        Cython signature: void FASTAEntry(FASTAEntry)
        """
        ...
    
    def headerMatches(self, rhs: FASTAEntry ) -> bool:
        """
        Cython signature: bool headerMatches(const FASTAEntry & rhs)
        """
        ...
    
    def sequenceMatches(self, rhs: FASTAEntry ) -> bool:
        """
        Cython signature: bool sequenceMatches(const FASTAEntry & rhs)
        """
        ... 


class FASTAFile:
    """
    Cython implementation of _FASTAFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1FASTAFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void FASTAFile()
        This class serves for reading in and writing FASTA files
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , data: List[FASTAEntry] ) -> None:
        """
        Cython signature: void load(const String & filename, libcpp_vector[FASTAEntry] & data)
        Loads a FASTA file given by 'filename' and stores the information in 'data'
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , data: List[FASTAEntry] ) -> None:
        """
        Cython signature: void store(const String & filename, libcpp_vector[FASTAEntry] & data)
        Stores the data given by 'data' at the file 'filename'
        """
        ...
    
    def readStart(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void readStart(const String & filename)
        Prepares a FASTA file given by 'filename' for streamed reading using readNext()
        
        :raises:
            Exception:FileNotFound is thrown if the file does not exists
        :raises:
            Exception:ParseError is thrown if the file does not suit to the standard
        """
        ...
    
    def readNext(self, protein: FASTAEntry ) -> bool:
        """
        Cython signature: bool readNext(FASTAEntry & protein)
        Reads the next FASTA entry from file
        
        If you want to read all entries in one go, use load()
        
        :return: true if entry was read; false if eof was reached
        :raises:
            Exception:FileNotFound is thrown if the file does not exists
        :raises:
            Exception:ParseError is thrown if the file does not suit to the standard
        """
        ...
    
    def atEnd(self) -> bool:
        """
        Cython signature: bool atEnd()
        Boolean function to check if streams is at end of file
        """
        ...
    
    def writeStart(self, filename: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void writeStart(const String & filename)
        Prepares a FASTA file given by 'filename' for streamed writing using writeNext()
        
        :raises:
            Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        """
        ...
    
    def writeNext(self, protein: FASTAEntry ) -> None:
        """
        Cython signature: void writeNext(const FASTAEntry & protein)
        Stores the data given by `protein`. Call writeStart() once before calling writeNext()
        
        Call writeEnd() when done to close the file!
        
        :raises:
            Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        """
        ...
    
    def writeEnd(self) -> None:
        """
        Cython signature: void writeEnd()
        Closes the file (flush). Called implicitly when FASTAFile object does out of scope
        """
        ... 


class IMSAlphabet:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IMSAlphabet()
        """
        ...
    
    @overload
    def __init__(self, in_0: IMSAlphabet ) -> None:
        """
        Cython signature: void IMSAlphabet(IMSAlphabet &)
        """
        ...
    
    @overload
    def __init__(self, elements: List[IMSElement] ) -> None:
        """
        Cython signature: void IMSAlphabet(libcpp_vector[IMSElement] & elements)
        """
        ...
    
    @overload
    def getElement(self, name: bytes ) -> IMSElement:
        """
        Cython signature: IMSElement getElement(libcpp_string & name)
        Gets the element with 'index' and returns element with the given index in alphabet
        """
        ...
    
    @overload
    def getElement(self, index: int ) -> IMSElement:
        """
        Cython signature: IMSElement getElement(int index)
        Gets the element with 'index'
        """
        ...
    
    def getName(self, index: int ) -> bytes:
        """
        Cython signature: libcpp_string getName(int index)
        Gets the symbol of the element with an 'index' in alphabet
        """
        ...
    
    @overload
    def getMass(self, name: bytes ) -> float:
        """
        Cython signature: double getMass(libcpp_string & name)
        Gets mono isotopic mass of the element with the symbol 'name'
        """
        ...
    
    @overload
    def getMass(self, index: int ) -> float:
        """
        Cython signature: double getMass(int index)
        Gets mass of the element with an 'index' in alphabet
        """
        ...
    
    def getMasses(self, isotope_index: int ) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getMasses(int isotope_index)
        Gets masses of elements isotopes given by 'isotope_index'
        """
        ...
    
    def getAverageMasses(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getAverageMasses()
        Gets average masses of elements
        """
        ...
    
    def hasName(self, name: bytes ) -> bool:
        """
        Cython signature: bool hasName(libcpp_string & name)
        Returns true if there is an element with symbol 'name' in the alphabet, false - otherwise
        """
        ...
    
    @overload
    def push_back(self, name: bytes , value: float ) -> None:
        """
        Cython signature: void push_back(libcpp_string & name, double value)
        Adds a new element with 'name' and mass 'value'
        """
        ...
    
    @overload
    def push_back(self, element: IMSElement ) -> None:
        """
        Cython signature: void push_back(IMSElement & element)
        Adds a new 'element' to the alphabet
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Clears the alphabet data
        """
        ...
    
    def sortByNames(self) -> None:
        """
        Cython signature: void sortByNames()
        Sorts the alphabet by names
        """
        ...
    
    def sortByValues(self) -> None:
        """
        Cython signature: void sortByValues()
        Sorts the alphabet by mass values
        """
        ...
    
    def load(self, fname: String ) -> None:
        """
        Cython signature: void load(String & fname)
        Loads the alphabet data from the file 'fname' using the default parser. If there is no file 'fname', throws an 'IOException'
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def setElement(self, name: bytes , mass: float , forced: bool ) -> None:
        """
        Cython signature: void setElement(libcpp_string & name, double mass, bool forced)
        Overwrites an element in the alphabet with the 'name' with a new element constructed from the given 'name' and 'mass'
        """
        ...
    
    def erase(self, name: bytes ) -> bool:
        """
        Cython signature: bool erase(libcpp_string & name)
        Removes the element with 'name' from the alphabet
        """
        ... 


class MapAlignmentAlgorithmIdentification:
    """
    Cython implementation of _MapAlignmentAlgorithmIdentification

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MapAlignmentAlgorithmIdentification.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MapAlignmentAlgorithmIdentification()
        """
        ...
    
    @overload
    def align(self, in_0: List[FeatureMap] , in_1: List[TransformationDescription] , in_2: int ) -> None:
        """
        Cython signature: void align(libcpp_vector[FeatureMap] &, libcpp_vector[TransformationDescription] &, int)
        """
        ...
    
    @overload
    def align(self, in_0: List[ConsensusMap] , in_1: List[TransformationDescription] , in_2: int ) -> None:
        """
        Cython signature: void align(libcpp_vector[ConsensusMap] &, libcpp_vector[TransformationDescription] &, int)
        """
        ...
    
    @overload
    def setReference(self, in_0: FeatureMap ) -> None:
        """
        Cython signature: void setReference(FeatureMap &)
        """
        ...
    
    @overload
    def setReference(self, in_0: ConsensusMap ) -> None:
        """
        Cython signature: void setReference(ConsensusMap &)
        """
        ...
    
    @overload
    def setReference(self, in_0: PeptideIdentificationList ) -> None:
        """
        Cython signature: void setReference(PeptideIdentificationList &)
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


class OSSpectrumMeta:
    """
    Cython implementation of _OSSpectrumMeta

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenSwath_1_1OSSpectrumMeta.html>`_
    """
    
    index: int
    
    id: bytes
    
    RT: float
    
    ms_level: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void OSSpectrumMeta()
        """
        ...
    
    @overload
    def __init__(self, in_0: OSSpectrumMeta ) -> None:
        """
        Cython signature: void OSSpectrumMeta(OSSpectrumMeta &)
        """
        ... 


class PScore:
    """
    Cython implementation of _PScore

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PScore.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PScore()
        """
        ...
    
    @overload
    def __init__(self, in_0: PScore ) -> None:
        """
        Cython signature: void PScore(PScore &)
        """
        ...
    
    def calculateIntensityRankInMZWindow(self, mz: List[float] , intensities: List[float] , mz_window: float ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] calculateIntensityRankInMZWindow(libcpp_vector[double] & mz, libcpp_vector[double] & intensities, double mz_window)
        Calculate local (windowed) peak ranks
        
        The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
        The result can be used to efficiently filter spectra for top 1..n peaks in mass windows
        
        
        :param mz: The m/z positions of the peaks
        :param intensities: The intensities of the peaks
        :param mz_window: The window in Thomson centered at each peak
        """
        ...
    
    def calculateRankMap(self, peak_map: MSExperiment , mz_window: float ) -> List[List[int]]:
        """
        Cython signature: libcpp_vector[libcpp_vector[size_t]] calculateRankMap(MSExperiment & peak_map, double mz_window)
        Precalculated, windowed peak ranks for a whole experiment
        
        The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity
        
        
        :param peak_map: Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra
        :param mz_window: Window in Thomson centered at each peak
        """
        ...
    
    def calculatePeakLevelSpectra(self, spec: MSSpectrum , ranks: List[int] , min_level: int , max_level: int ) -> Dict[int, MSSpectrum]:
        """
        Cython signature: libcpp_map[size_t,MSSpectrum] calculatePeakLevelSpectra(MSSpectrum & spec, libcpp_vector[size_t] & ranks, size_t min_level, size_t max_level)
        Calculates spectra for peak level between min_level to max_level and stores them in the map
        
        A spectrum of peak level n retains the (n+1) top intensity peaks in a sliding mz_window centered at each peak
        """
        ...
    
    @overload
    def computePScore(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectra: List[MSSpectrum] , mz_window: float ) -> float:
        """
        Cython signature: double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, libcpp_map[size_t,MSSpectrum] & peak_level_spectra, libcpp_vector[MSSpectrum] & theo_spectra, double mz_window)
        Computes the PScore for a vector of theoretical spectra
        
        Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
        The best score obtained by scoring all those theoretical spectra against the experimental ones is returned
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank).
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
        """
        ...
    
    @overload
    def computePScore(self, fragment_mass_tolerance: float , fragment_mass_tolerance_unit_ppm: bool , peak_level_spectra: Dict[int, MSSpectrum] , theo_spectrum: MSSpectrum , mz_window: float ) -> float:
        """
        Cython signature: double computePScore(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, libcpp_map[size_t,MSSpectrum] & peak_level_spectra, MSSpectrum & theo_spectrum, double mz_window)
        Computes the PScore for a single theoretical spectrum
        
        
        :param fragment_mass_tolerance: Mass tolerance for matching peaks
        :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
        :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank)
        :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
        :param mz_window: Window in Thomson centered at each peak
        """
        ... 


class ProteaseDigestion:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ProteaseDigestion()
        """
        ...
    
    @overload
    def __init__(self, in_0: ProteaseDigestion ) -> None:
        """
        Cython signature: void ProteaseDigestion(ProteaseDigestion &)
        """
        ...
    
    @overload
    def setEnzyme(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setEnzyme(String name)
        Sets the enzyme for the digestion (by name)
        """
        ...
    
    @overload
    def setEnzyme(self, enzyme: DigestionEnzyme ) -> None:
        """
        Cython signature: void setEnzyme(DigestionEnzyme * enzyme)
        Sets the enzyme for the digestion
        """
        ...
    
    @overload
    def digest(self, protein: AASequence , output: List[AASequence] ) -> int:
        """
        Cython signature: size_t digest(AASequence & protein, libcpp_vector[AASequence] & output)
        """
        ...
    
    @overload
    def digest(self, protein: AASequence , output: List[AASequence] , min_length: int , max_length: int ) -> int:
        """
        Cython signature: size_t digest(AASequence & protein, libcpp_vector[AASequence] & output, size_t min_length, size_t max_length)
          Performs the enzymatic digestion of a protein.
        
        
          :param protein: Sequence to digest
          :param output: Digestion products (peptides)
          :param min_length: Minimal length of reported products
          :param max_length: Maximal length of reported products (0 = no restriction)
          :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        ...
    
    def peptideCount(self, protein: AASequence ) -> int:
        """
        Cython signature: size_t peptideCount(AASequence & protein)
        Returns the number of peptides a digestion of protein would yield under the current enzyme and missed cleavage settings
        """
        ...
    
    @overload
    def isValidProduct(self, protein: AASequence , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool:
        """
        Cython signature: bool isValidProduct(AASequence protein, size_t pep_pos, size_t pep_length, bool ignore_missed_cleavages, bool methionine_cleavage)
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
        ...
    
    @overload
    def isValidProduct(self, protein: Union[bytes, str, String] , pep_pos: int , pep_length: int , ignore_missed_cleavages: bool , methionine_cleavage: bool ) -> bool:
        """
        Cython signature: bool isValidProduct(String protein, size_t pep_pos, size_t pep_length, bool ignore_missed_cleavages, bool methionine_cleavage)
        Forwards to isValidProduct using protein.toUnmodifiedString()
        """
        ...
    
    @overload
    def isValidProduct(self, sequence: Union[bytes, str, String] , pos: int , length: int , ignore_missed_cleavages: bool ) -> bool:
        """
        Cython signature: bool isValidProduct(String sequence, int pos, int length, bool ignore_missed_cleavages)
        Boolean operator returns true if the peptide fragment starting at position `pos` with length `length` within the sequence `sequence` generated by the current enzyme\n
        Checks if peptide is a valid digestion product of the enzyme, taking into account specificity and the MC flag provided here
        
        
        :param protein: Protein sequence
        :param pep_pos: Starting index of potential peptide
        :param pep_length: Length of potential peptide
        :param ignore_missed_cleavages: Do not compare MC's of potential peptide to the maximum allowed MC's
        :return: True if peptide has correct n/c terminals (according to enzyme, specificity and missed cleavages)
        """
        ...
    
    def getMissedCleavages(self) -> int:
        """
        Cython signature: size_t getMissedCleavages()
        Returns the max. number of allowed missed cleavages for the digestion
        """
        ...
    
    def setMissedCleavages(self, missed_cleavages: int ) -> None:
        """
        Cython signature: void setMissedCleavages(size_t missed_cleavages)
        Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used
        """
        ...
    
    def countInternalCleavageSites(self, sequence: Union[bytes, str, String] ) -> int:
        """
        Cython signature: size_t countInternalCleavageSites(String sequence)
        Returns the number of internal cleavage sites for this sequence.
        """
        ...
    
    def getEnzymeName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getEnzymeName()
        Returns the enzyme for the digestion
        """
        ...
    
    def getSpecificity(self) -> int:
        """
        Cython signature: Specificity getSpecificity()
        Returns the specificity for the digestion
        """
        ...
    
    def setSpecificity(self, spec: int ) -> None:
        """
        Cython signature: void setSpecificity(Specificity spec)
        Sets the specificity for the digestion (default is SPEC_FULL)
        """
        ...
    
    def getSpecificityByName(self, name: Union[bytes, str, String] ) -> int:
        """
        Cython signature: Specificity getSpecificityByName(String name)
        Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid
        """
        ...
    
    def digestUnmodified(self, sequence: StringView , output: List[StringView] , min_length: int , max_length: int ) -> int:
        """
        Cython signature: size_t digestUnmodified(StringView sequence, libcpp_vector[StringView] & output, size_t min_length, size_t max_length)
        Performs the enzymatic digestion of an unmodified sequence\n
        By returning only references into the original string this is very fast
        
        
        :param sequence: Sequence to digest
        :param output: Digestion products
        :param min_length: Minimal length of reported products
        :param max_length: Maximal length of reported products (0 = no restriction)
        :return: Number of discarded digestion products (which are not matching length restrictions)
        """
        ... 


class Residue:
    """
    Cython implementation of _Residue

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Residue.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Residue()
        """
        ...
    
    @overload
    def __init__(self, in_0: Residue ) -> None:
        """
        Cython signature: void Residue(Residue &)
        """
        ...
    
    @overload
    def __init__(self, name: Union[bytes, str, String] , three_letter_code: Union[bytes, str, String] , one_letter_code: Union[bytes, str, String] , formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void Residue(String name, String three_letter_code, String one_letter_code, EmpiricalFormula formula)
        """
        ...
    
    def getInternalToFull(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToFull()
        """
        ...
    
    def getInternalToNTerm(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToNTerm()
        """
        ...
    
    def getInternalToCTerm(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToCTerm()
        """
        ...
    
    def getInternalToAIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToAIon()
        """
        ...
    
    def getInternalToBIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToBIon()
        """
        ...
    
    def getInternalToCIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToCIon()
        """
        ...
    
    def getInternalToXIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToXIon()
        """
        ...
    
    def getInternalToYIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToYIon()
        """
        ...
    
    def getInternalToZIon(self) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getInternalToZIon()
        """
        ...
    
    def getResidueTypeName(self, res_type: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String getResidueTypeName(ResidueType res_type)
        Returns the ion name given as a residue type
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String name)
        Sets the name of the residue
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the residue
        """
        ...
    
    def setSynonyms(self, synonyms: Set[bytes] ) -> None:
        """
        Cython signature: void setSynonyms(libcpp_set[String] synonyms)
        Sets the synonyms
        """
        ...
    
    def addSynonym(self, synonym: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addSynonym(String synonym)
        Adds a synonym
        """
        ...
    
    def getSynonyms(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getSynonyms()
        Returns the sysnonyms
        """
        ...
    
    def setThreeLetterCode(self, three_letter_code: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setThreeLetterCode(String three_letter_code)
        Sets the name of the residue as three letter code
        """
        ...
    
    def getThreeLetterCode(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getThreeLetterCode()
        Returns the name of the residue as three letter code
        """
        ...
    
    def setOneLetterCode(self, one_letter_code: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setOneLetterCode(String one_letter_code)
        Sets the name as one letter code
        """
        ...
    
    def getOneLetterCode(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getOneLetterCode()
        Returns the name as one letter code
        """
        ...
    
    def addLossFormula(self, in_0: EmpiricalFormula ) -> None:
        """
        Cython signature: void addLossFormula(EmpiricalFormula)
        Adds a neutral loss formula
        """
        ...
    
    def setLossFormulas(self, in_0: List[EmpiricalFormula] ) -> None:
        """
        Cython signature: void setLossFormulas(libcpp_vector[EmpiricalFormula])
        Sets the neutral loss formulas
        """
        ...
    
    def addNTermLossFormula(self, in_0: EmpiricalFormula ) -> None:
        """
        Cython signature: void addNTermLossFormula(EmpiricalFormula)
        Adds N-terminal losses
        """
        ...
    
    def setNTermLossFormulas(self, in_0: List[EmpiricalFormula] ) -> None:
        """
        Cython signature: void setNTermLossFormulas(libcpp_vector[EmpiricalFormula])
        Sets the N-terminal losses
        """
        ...
    
    def getLossFormulas(self) -> List[EmpiricalFormula]:
        """
        Cython signature: libcpp_vector[EmpiricalFormula] getLossFormulas()
        Returns the neutral loss formulas
        """
        ...
    
    def getNTermLossFormulas(self) -> List[EmpiricalFormula]:
        """
        Cython signature: libcpp_vector[EmpiricalFormula] getNTermLossFormulas()
        Returns N-terminal loss formulas
        """
        ...
    
    def setLossNames(self, name: List[bytes] ) -> None:
        """
        Cython signature: void setLossNames(libcpp_vector[String] name)
        Sets the neutral loss molecule name
        """
        ...
    
    def setNTermLossNames(self, name: List[bytes] ) -> None:
        """
        Cython signature: void setNTermLossNames(libcpp_vector[String] name)
        Sets the N-terminal loss names
        """
        ...
    
    def addLossName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addLossName(String name)
        Adds neutral loss molecule name
        """
        ...
    
    def addNTermLossName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addNTermLossName(String name)
        Adds a N-terminal loss name
        """
        ...
    
    def getLossNames(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getLossNames()
        Gets neutral loss name (if there is one, else returns an empty string)
        """
        ...
    
    def getNTermLossNames(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getNTermLossNames()
        Returns the N-terminal loss names
        """
        ...
    
    def setFormula(self, formula: EmpiricalFormula ) -> None:
        """
        Cython signature: void setFormula(EmpiricalFormula formula)
        Sets empirical formula of the residue (must be full, with N and C-terminus)
        """
        ...
    
    @overload
    def getFormula(self, ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula()
        Returns the empirical formula of the residue
        """
        ...
    
    @overload
    def getFormula(self, res_type: int ) -> EmpiricalFormula:
        """
        Cython signature: EmpiricalFormula getFormula(ResidueType res_type)
        """
        ...
    
    def setAverageWeight(self, weight: float ) -> None:
        """
        Cython signature: void setAverageWeight(double weight)
        Sets average weight of the residue (must be full, with N and C-terminus)
        """
        ...
    
    @overload
    def getAverageWeight(self, ) -> float:
        """
        Cython signature: double getAverageWeight()
        Returns average weight of the residue
        """
        ...
    
    @overload
    def getAverageWeight(self, res_type: int ) -> float:
        """
        Cython signature: double getAverageWeight(ResidueType res_type)
        """
        ...
    
    def setMonoWeight(self, weight: float ) -> None:
        """
        Cython signature: void setMonoWeight(double weight)
        Sets monoisotopic weight of the residue (must be full, with N and C-terminus)
        """
        ...
    
    @overload
    def getMonoWeight(self, ) -> float:
        """
        Cython signature: double getMonoWeight()
        Returns monoisotopic weight of the residue
        """
        ...
    
    @overload
    def getMonoWeight(self, res_type: int ) -> float:
        """
        Cython signature: double getMonoWeight(ResidueType res_type)
        """
        ...
    
    def getModification(self) -> ResidueModification:
        """
        Cython signature: const ResidueModification * getModification()
        """
        ...
    
    @overload
    def setModification(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setModification(String name)
        Sets the modification by name; the mod should be present in ModificationsDB
        """
        ...
    
    @overload
    def setModification(self, mod: ResidueModification ) -> None:
        """
        Cython signature: void setModification(const ResidueModification & mod)
        Sets the modification by a ResidueModification object; checks if present in ModificationsDB and adds if not.
        """
        ...
    
    def setModificationByDiffMonoMass(self, diffMonoMass: float ) -> None:
        """
        Cython signature: void setModificationByDiffMonoMass(double diffMonoMass)
        Sets the modification by monoisotopic mass difference in Da; checks if present in ModificationsDB with tolerance and adds a "user-defined" modification if not (for later lookups).
        """
        ...
    
    def getModificationName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getModificationName()
        Returns the name of the modification to the modification
        """
        ...
    
    def setLowMassIons(self, low_mass_ions: List[EmpiricalFormula] ) -> None:
        """
        Cython signature: void setLowMassIons(libcpp_vector[EmpiricalFormula] low_mass_ions)
        Sets the low mass marker ions as a vector of formulas
        """
        ...
    
    def getLowMassIons(self) -> List[EmpiricalFormula]:
        """
        Cython signature: libcpp_vector[EmpiricalFormula] getLowMassIons()
        Returns a vector of formulas with the low mass markers of the residue
        """
        ...
    
    def setResidueSets(self, residues_sets: Set[bytes] ) -> None:
        """
        Cython signature: void setResidueSets(libcpp_set[String] residues_sets)
        Sets the residue sets the amino acid is contained in
        """
        ...
    
    def addResidueSet(self, residue_sets: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addResidueSet(String residue_sets)
        Adds a residue set to the residue sets
        """
        ...
    
    def getResidueSets(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getResidueSets()
        Returns the residue sets this residue is contained in
        """
        ...
    
    def hasNeutralLoss(self) -> bool:
        """
        Cython signature: bool hasNeutralLoss()
        True if the residue has neutral loss
        """
        ...
    
    def hasNTermNeutralLosses(self) -> bool:
        """
        Cython signature: bool hasNTermNeutralLosses()
        True if N-terminal neutral losses are set
        """
        ...
    
    def getPka(self) -> float:
        """
        Cython signature: double getPka()
        Returns the pka of the residue
        """
        ...
    
    def getPkb(self) -> float:
        """
        Cython signature: double getPkb()
        Returns the pkb of the residue
        """
        ...
    
    def getPkc(self) -> float:
        """
        Cython signature: double getPkc()
        Returns the pkc of the residue if it exists otherwise -1
        """
        ...
    
    def getPiValue(self) -> float:
        """
        Cython signature: double getPiValue()
        Calculates the isoelectric point using the pk values
        """
        ...
    
    def setPka(self, value: float ) -> None:
        """
        Cython signature: void setPka(double value)
        Sets the pka of the residue
        """
        ...
    
    def setPkb(self, value: float ) -> None:
        """
        Cython signature: void setPkb(double value)
        Sets the pkb of the residue
        """
        ...
    
    def setPkc(self, value: float ) -> None:
        """
        Cython signature: void setPkc(double value)
        Sets the pkc of the residue
        """
        ...
    
    def getSideChainBasicity(self) -> float:
        """
        Cython signature: double getSideChainBasicity()
        Returns the side chain basicity
        """
        ...
    
    def setSideChainBasicity(self, gb_sc: float ) -> None:
        """
        Cython signature: void setSideChainBasicity(double gb_sc)
        Sets the side chain basicity
        """
        ...
    
    def getBackboneBasicityLeft(self) -> float:
        """
        Cython signature: double getBackboneBasicityLeft()
        Returns the backbone basicitiy if located in N-terminal direction
        """
        ...
    
    def setBackboneBasicityLeft(self, gb_bb_l: float ) -> None:
        """
        Cython signature: void setBackboneBasicityLeft(double gb_bb_l)
        Sets the N-terminal direction backbone basicitiy
        """
        ...
    
    def getBackboneBasicityRight(self) -> float:
        """
        Cython signature: double getBackboneBasicityRight()
        Returns the C-terminal direction backbone basicitiy
        """
        ...
    
    def setBackboneBasicityRight(self, gb_bb_r: float ) -> None:
        """
        Cython signature: void setBackboneBasicityRight(double gb_bb_r)
        Sets the C-terminal direction backbone basicity
        """
        ...
    
    def isModified(self) -> bool:
        """
        Cython signature: bool isModified()
        True if the residue is a modified one
        """
        ...
    
    def isInResidueSet(self, residue_set: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool isInResidueSet(String residue_set)
        True if the residue is contained in the set
        """
        ...
    
    def residueTypeToIonLetter(self, res_type: int ) -> Union[bytes, str, String]:
        """
        Cython signature: String residueTypeToIonLetter(ResidueType res_type)
        Helper for mapping residue types to letters for Text annotations and labels
        """
        ...
    
    def __richcmp__(self, other: Residue, op: int) -> Any:
        ...
    ResidueType : __ResidueType 


class RichPeak2D:
    """
    Cython implementation of _RichPeak2D

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1RichPeak2D.html>`_
      -- Inherits from ['Peak2D', 'UniqueIdInterface', 'MetaInfoInterface']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RichPeak2D()
        A 2-dimensional raw data point or peak with meta information
        """
        ...
    
    @overload
    def __init__(self, in_0: RichPeak2D ) -> None:
        """
        Cython signature: void RichPeak2D(RichPeak2D &)
        """
        ...
    
    def getIntensity(self) -> float:
        """
        Cython signature: float getIntensity()
        Returns the data point intensity (height)
        """
        ...
    
    def getMZ(self) -> float:
        """
        Cython signature: double getMZ()
        Returns the m/z coordinate (index 1)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        Returns the RT coordinate (index 0)
        """
        ...
    
    def setMZ(self, in_0: float ) -> None:
        """
        Cython signature: void setMZ(double)
        Returns the m/z coordinate (index 1)
        """
        ...
    
    def setRT(self, in_0: float ) -> None:
        """
        Cython signature: void setRT(double)
        Returns the RT coordinate (index 0)
        """
        ...
    
    def setIntensity(self, in_0: float ) -> None:
        """
        Cython signature: void setIntensity(float)
        Returns the data point intensity (height)
        """
        ...
    
    def getUniqueId(self) -> int:
        """
        Cython signature: size_t getUniqueId()
        Returns the unique id
        """
        ...
    
    def clearUniqueId(self) -> int:
        """
        Cython signature: size_t clearUniqueId()
        Clear the unique id. The new unique id will be invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        ...
    
    def hasValidUniqueId(self) -> int:
        """
        Cython signature: size_t hasValidUniqueId()
        Returns whether the unique id is valid. Returns 1 if the unique id is valid, 0 otherwise
        """
        ...
    
    def hasInvalidUniqueId(self) -> int:
        """
        Cython signature: size_t hasInvalidUniqueId()
        Returns whether the unique id is invalid. Returns 1 if the unique id is invalid, 0 otherwise
        """
        ...
    
    def setUniqueId(self, rhs: int ) -> None:
        """
        Cython signature: void setUniqueId(uint64_t rhs)
        Assigns a new, valid unique id. Always returns 1
        """
        ...
    
    def ensureUniqueId(self) -> int:
        """
        Cython signature: size_t ensureUniqueId()
        Assigns a valid unique id, but only if the present one is invalid. Returns 1 if the unique id was changed, 0 otherwise
        """
        ...
    
    def isValid(self, unique_id: int ) -> bool:
        """
        Cython signature: bool isValid(uint64_t unique_id)
        Returns true if the unique_id is valid, false otherwise
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
    
    def __richcmp__(self, other: RichPeak2D, op: int) -> Any:
        ... 


class SavitzkyGolayFilter:
    """
    Cython implementation of _SavitzkyGolayFilter

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SavitzkyGolayFilter.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SavitzkyGolayFilter()
        """
        ...
    
    @overload
    def __init__(self, in_0: SavitzkyGolayFilter ) -> None:
        """
        Cython signature: void SavitzkyGolayFilter(SavitzkyGolayFilter &)
        """
        ...
    
    def filter(self, spectrum: MSSpectrum ) -> None:
        """
        Cython signature: void filter(MSSpectrum & spectrum)
        Removed the noise from an MSSpectrum containing profile data
        """
        ...
    
    def filterExperiment(self, exp: MSExperiment ) -> None:
        """
        Cython signature: void filterExperiment(MSExperiment & exp)
        Removed the noise from an MSExperiment containing profile data
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


class SignalToNoiseEstimatorMeanIterative:
    """
    Cython implementation of _SignalToNoiseEstimatorMeanIterative[_MSSpectrum]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SignalToNoiseEstimatorMeanIterative[_MSSpectrum].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMeanIterative()
        """
        ...
    
    @overload
    def __init__(self, in_0: SignalToNoiseEstimatorMeanIterative ) -> None:
        """
        Cython signature: void SignalToNoiseEstimatorMeanIterative(SignalToNoiseEstimatorMeanIterative &)
        """
        ...
    
    def init(self, c: MSSpectrum ) -> None:
        """
        Cython signature: void init(MSSpectrum & c)
        """
        ...
    
    def getSignalToNoise(self, index: int ) -> float:
        """
        Cython signature: double getSignalToNoise(size_t index)
        """
        ...
    IntensityThresholdCalculation : __IntensityThresholdCalculation 


class TMTSixPlexQuantitationMethod:
    """
    Cython implementation of _TMTSixPlexQuantitationMethod

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TMTSixPlexQuantitationMethod.html>`_
      -- Inherits from ['IsobaricQuantitationMethod']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void TMTSixPlexQuantitationMethod()
        """
        ...
    
    @overload
    def __init__(self, in_0: TMTSixPlexQuantitationMethod ) -> None:
        """
        Cython signature: void TMTSixPlexQuantitationMethod(TMTSixPlexQuantitationMethod &)
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


class __IntensityThresholdCalculation:
    None
    MANUAL : int
    AUTOMAXBYSTDEV : int
    AUTOMAXBYPERCENT : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class __ResidueType:
    None
    Full : int
    Internal : int
    NTerminal : int
    CTerminal : int
    AIon : int
    BIon : int
    CIon : int
    XIon : int
    YIon : int
    ZIon : int
    Precursor_ion : int
    BIonMinusH20 : int
    YIonMinusH20 : int
    BIonMinusNH3 : int
    YIonMinusNH3 : int
    NonIdentified : int
    Unannotated : int
    SizeOfResidueType : int

    def getMapping(self) -> Dict[int, str]:
       ... 

