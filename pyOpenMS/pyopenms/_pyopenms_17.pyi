from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_File_absolutePath(file: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String absolutePath(String file)
    """
    ...

def __static_File_basename(file: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String basename(String file)
    """
    ...

def __static_File_empty(file: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool empty(String file)
    """
    ...

def __static_File_exists(file: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool exists(String file)
    """
    ...

def __static_File_fileList(dir: Union[bytes, str, String] , file_pattern: Union[bytes, str, String] , output: List[bytes] , full_path: bool ) -> bool:
    """
    Cython signature: bool fileList(String dir, String file_pattern, StringList output, bool full_path)
    """
    ...

def __static_File_find(filename: Union[bytes, str, String] , directories: List[bytes] ) -> Union[bytes, str, String]:
    """
    Cython signature: String find(String filename, StringList directories)
    """
    ...

def __static_File_findDatabase(db_name: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String findDatabase(String db_name)
    """
    ...

def __static_File_findDoc(filename: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String findDoc(String filename)
    """
    ...

def __static_File_findExecutable(toolName: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String findExecutable(String toolName)
    """
    ...

def __static_File_getExecutablePath() -> Union[bytes, str, String]:
    """
    Cython signature: String getExecutablePath()
    """
    ...

def __static_File_getOpenMSDataPath() -> Union[bytes, str, String]:
    """
    Cython signature: String getOpenMSDataPath()
    """
    ...

def __static_File_getOpenMSHomePath() -> Union[bytes, str, String]:
    """
    Cython signature: String getOpenMSHomePath()
    """
    ...

def __static_File_getSystemParameters() -> Param:
    """
    Cython signature: Param getSystemParameters()
    """
    ...

def __static_File_getTempDirectory() -> Union[bytes, str, String]:
    """
    Cython signature: String getTempDirectory()
    """
    ...

def __static_File_getTemporaryFile(alternative_file: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String getTemporaryFile(const String & alternative_file)
    """
    ...

def __static_File_getUniqueName() -> Union[bytes, str, String]:
    """
    Cython signature: String getUniqueName()
    """
    ...

def __static_File_getUserDirectory() -> Union[bytes, str, String]:
    """
    Cython signature: String getUserDirectory()
    """
    ...

def __static_File_isDirectory(path: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool isDirectory(String path)
    """
    ...

def __static_File_path(file: Union[bytes, str, String] ) -> Union[bytes, str, String]:
    """
    Cython signature: String path(String file)
    """
    ...

def __static_File_readable(file: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool readable(String file)
    """
    ...

def __static_File_remove(file: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool remove(String file)
    """
    ...

def __static_File_removeDirRecursively(dir_name: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool removeDirRecursively(String dir_name)
    """
    ...

def __static_File_rename(old_filename: Union[bytes, str, String] , new_filename: Union[bytes, str, String] , overwrite_existing: bool , verbose: bool ) -> bool:
    """
    Cython signature: bool rename(const String & old_filename, const String & new_filename, bool overwrite_existing, bool verbose)
    """
    ...

def __static_File_writable(file: Union[bytes, str, String] ) -> bool:
    """
    Cython signature: bool writable(String file)
    """
    ...


class ClusteringGrid:
    """
    Cython implementation of _ClusteringGrid

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ClusteringGrid.html>`_
    """
    
    @overload
    def __init__(self, grid_spacing_x: List[float] , grid_spacing_y: List[float] ) -> None:
        """
        Cython signature: void ClusteringGrid(libcpp_vector[double] & grid_spacing_x, libcpp_vector[double] & grid_spacing_y)
        """
        ...
    
    @overload
    def __init__(self, in_0: ClusteringGrid ) -> None:
        """
        Cython signature: void ClusteringGrid(ClusteringGrid &)
        """
        ...
    
    def getGridSpacingX(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getGridSpacingX()
        """
        ...
    
    def getGridSpacingY(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getGridSpacingY()
        """
        ...
    
    def addCluster(self, cell_index: List[int, int] , cluster_index: int ) -> None:
        """
        Cython signature: void addCluster(libcpp_pair[int,int] cell_index, int & cluster_index)
        Adds a cluster to this grid cell
        """
        ...
    
    def removeCluster(self, cell_index: List[int, int] , cluster_index: int ) -> None:
        """
        Cython signature: void removeCluster(libcpp_pair[int,int] cell_index, int & cluster_index)
        Removes a cluster from this grid cell and removes the cell if no other cluster left
        """
        ...
    
    def removeAllClusters(self) -> None:
        """
        Cython signature: void removeAllClusters()
        Removes all clusters from this grid (and hence all cells)
        """
        ...
    
    def getIndex(self, position: Union[Sequence[int], Sequence[float]] ) -> List[int, int]:
        """
        Cython signature: libcpp_pair[int,int] getIndex(DPosition2 position)
        """
        ...
    
    def isNonEmptyCell(self, cell_index: List[int, int] ) -> bool:
        """
        Cython signature: bool isNonEmptyCell(libcpp_pair[int,int] cell_index)
        Checks if there are clusters at this cell index
        """
        ...
    
    def getCellCount(self) -> int:
        """
        Cython signature: int getCellCount()
        Returns number of grid cells occupied by one or more clusters
        """
        ... 


class CoarseIsotopePatternGenerator:
    """
    Cython implementation of _CoarseIsotopePatternGenerator

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CoarseIsotopePatternGenerator.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CoarseIsotopePatternGenerator()
        """
        ...
    
    @overload
    def __init__(self, max_isotope: int ) -> None:
        """
        Cython signature: void CoarseIsotopePatternGenerator(size_t max_isotope)
        """
        ...
    
    @overload
    def __init__(self, max_isotope: int , round_masses: bool ) -> None:
        """
        Cython signature: void CoarseIsotopePatternGenerator(size_t max_isotope, bool round_masses)
        """
        ...
    
    def run(self, in_0: EmpiricalFormula ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution run(EmpiricalFormula)
        """
        ...
    
    def getRoundMasses(self) -> bool:
        """
        Cython signature: bool getRoundMasses()
        Returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)
        """
        ...
    
    def setRoundMasses(self, round_masses_: bool ) -> None:
        """
        Cython signature: void setRoundMasses(bool round_masses_)
        Sets the round_masses_ flag to round masses to integer values (true) or return accurate masses (false)
        """
        ...
    
    def getMaxIsotope(self) -> int:
        """
        Cython signature: size_t getMaxIsotope()
        Returns the currently set maximum isotope
        """
        ...
    
    def setMaxIsotope(self, max_isotope: int ) -> None:
        """
        Cython signature: void setMaxIsotope(size_t max_isotope)
        Sets the maximal isotope with 'max_isotope'
        """
        ...
    
    def estimateFromPeptideWeight(self, average_weight: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromPeptideWeight(double average_weight)
        Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported
        """
        ...
    
    def estimateFromPeptideWeightAndS(self, average_weight: float , S: int ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromPeptideWeightAndS(double average_weight, unsigned int S)
        Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs
        """
        ...
    
    def estimateFromRNAWeight(self, average_weight: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromRNAWeight(double average_weight)
        Estimate Nucleotide Isotopedistribution from weight
        """
        ...
    
    def estimateFromDNAWeight(self, average_weight: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromDNAWeight(double average_weight)
        Estimate Nucleotide Isotopedistribution from weight
        """
        ...
    
    def estimateFromWeightAndComp(self, average_weight: float , C: float , H: float , N: float , O: float , S: float , P: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromWeightAndComp(double average_weight, double C, double H, double N, double O, double S, double P)
        """
        ...
    
    def estimateFromWeightAndCompAndS(self, average_weight: float , S: int , C: float , H: float , N: float , O: float , P: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateFromWeightAndCompAndS(double average_weight, unsigned int S, double C, double H, double N, double O, double P)
        Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition
        """
        ...
    
    def estimateForFragmentFromPeptideWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateForFragmentFromPeptideWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[unsigned int] & precursor_isotopes)
        Estimate peptide fragment IsotopeDistribution from the precursor's average weight, fragment's average weight, and a set of isolated precursor isotopes
        """
        ...
    
    def estimateForFragmentFromPeptideWeightAndS(self, average_weight_precursor: float , S_precursor: int , average_weight_fragment: float , S_fragment: int , precursor_isotopes: Set[int] ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateForFragmentFromPeptideWeightAndS(double average_weight_precursor, unsigned int S_precursor, double average_weight_fragment, unsigned int S_fragment, libcpp_set[unsigned int] & precursor_isotopes)
        Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
        number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
        and a set of isolated precursor isotopes.
        """
        ...
    
    def approximateFromPeptideWeight(self, mass: float , num_peaks: int , charge: int ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution approximateFromPeptideWeight(double mass, unsigned int num_peaks, unsigned int charge)
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
        ...
    
    def approximateIntensities(self, mass: float , num_peaks: int ) -> List[float]:
        """
        Cython signature: libcpp_vector[double] approximateIntensities(double mass, unsigned int num_peaks)
        Roughly approximate peptidic isotope pattern intensities from monoisotopic weight using Poisson distribution.
        Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
        This method is around 100 times faster than estimateFromPeptideWeight, but only an approximation, see approximateFromPeptideWeight.
        """
        ...
    
    def estimateForFragmentFromRNAWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateForFragmentFromRNAWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[unsigned int] & precursor_isotopes)
        Estimate RNA fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, and a set of isolated precursor isotopes
        """
        ...
    
    def estimateForFragmentFromDNAWeight(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateForFragmentFromDNAWeight(double average_weight_precursor, double average_weight_fragment, libcpp_set[unsigned int] & precursor_isotopes)
        Estimate DNA fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, and a set of isolated precursor isotopes.
        """
        ...
    
    def estimateForFragmentFromWeightAndComp(self, average_weight_precursor: float , average_weight_fragment: float , precursor_isotopes: Set[int] , C: float , H: float , N: float , O: float , S: float , P: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution estimateForFragmentFromWeightAndComp(double average_weight_precursor, double average_weight_fragment, libcpp_set[unsigned int] & precursor_isotopes, double C, double H, double N, double O, double S, double P)
        Estimate fragment IsotopeDistribution from the precursor's average weight,
        fragment's average weight, a set of isolated precursor isotopes, and average composition
        """
        ...
    
    def calcFragmentIsotopeDist(self, fragment_isotope_dist: IsotopeDistribution , comp_fragment_isotope_dist: IsotopeDistribution , precursor_isotopes: Set[int] , fragment_mono_mass: float ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution calcFragmentIsotopeDist(IsotopeDistribution & fragment_isotope_dist, IsotopeDistribution & comp_fragment_isotope_dist, libcpp_set[unsigned int] & precursor_isotopes, double fragment_mono_mass)
        Calculate isotopic distribution for a fragment molecule
        """
        ... 


class ConsensusMapNormalizerAlgorithmThreshold:
    """
    Cython implementation of _ConsensusMapNormalizerAlgorithmThreshold

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ConsensusMapNormalizerAlgorithmThreshold.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void ConsensusMapNormalizerAlgorithmThreshold()
        """
        ...
    
    def computeCorrelation(self, input_map: ConsensusMap , ratio_threshold: float , acc_filter: Union[bytes, str, String] , desc_filter: Union[bytes, str, String] ) -> List[float]:
        """
        Cython signature: libcpp_vector[double] computeCorrelation(ConsensusMap & input_map, double ratio_threshold, const String & acc_filter, const String & desc_filter)
        Determines the ratio of all maps to the map with the most features
        """
        ...
    
    def normalizeMaps(self, input_map: ConsensusMap , ratios: List[float] ) -> None:
        """
        Cython signature: void normalizeMaps(ConsensusMap & input_map, libcpp_vector[double] & ratios)
        Applies the given ratio to the maps of the consensusMap
        """
        ... 


class DocumentIdentifier:
    """
    Cython implementation of _DocumentIdentifier

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1DocumentIdentifier.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void DocumentIdentifier()
        """
        ...
    
    @overload
    def __init__(self, in_0: DocumentIdentifier ) -> None:
        """
        Cython signature: void DocumentIdentifier(DocumentIdentifier &)
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


class File:
    """
    Cython implementation of _File

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1File.html>`_
    """
    
    absolutePath: __static_File_absolutePath
    
    basename: __static_File_basename
    
    empty: __static_File_empty
    
    exists: __static_File_exists
    
    fileList: __static_File_fileList
    
    find: __static_File_find
    
    findDatabase: __static_File_findDatabase
    
    findDoc: __static_File_findDoc
    
    findExecutable: __static_File_findExecutable
    
    getExecutablePath: __static_File_getExecutablePath
    
    getOpenMSDataPath: __static_File_getOpenMSDataPath
    
    getOpenMSHomePath: __static_File_getOpenMSHomePath
    
    getSystemParameters: __static_File_getSystemParameters
    
    getTempDirectory: __static_File_getTempDirectory
    
    getTemporaryFile: __static_File_getTemporaryFile
    
    getUniqueName: __static_File_getUniqueName
    
    getUserDirectory: __static_File_getUserDirectory
    
    isDirectory: __static_File_isDirectory
    
    path: __static_File_path
    
    readable: __static_File_readable
    
    remove: __static_File_remove
    
    removeDirRecursively: __static_File_removeDirRecursively
    
    rename: __static_File_rename
    
    writable: __static_File_writable 


class FineIsotopePatternGenerator:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void FineIsotopePatternGenerator()
        """
        ...
    
    @overload
    def __init__(self, threshold: float ) -> None:
        """
        Cython signature: void FineIsotopePatternGenerator(double threshold)
        """
        ...
    
    @overload
    def __init__(self, threshold: float , use_total_prob: bool ) -> None:
        """
        Cython signature: void FineIsotopePatternGenerator(double threshold, bool use_total_prob)
        """
        ...
    
    @overload
    def __init__(self, threshold: float , use_total_prob: bool , absolute: bool ) -> None:
        """
        Cython signature: void FineIsotopePatternGenerator(double threshold, bool use_total_prob, bool absolute)
        """
        ...
    
    def setThreshold(self, threshold: float ) -> None:
        """
        Cython signature: void setThreshold(double threshold)
        """
        ...
    
    def getThreshold(self) -> float:
        """
        Cython signature: double getThreshold()
        """
        ...
    
    def setAbsolute(self, absolute: bool ) -> None:
        """
        Cython signature: void setAbsolute(bool absolute)
        """
        ...
    
    def getAbsolute(self) -> bool:
        """
        Cython signature: bool getAbsolute()
        """
        ...
    
    def setTotalProbability(self, total: bool ) -> None:
        """
        Cython signature: void setTotalProbability(bool total)
        """
        ...
    
    def getTotalProbability(self) -> bool:
        """
        Cython signature: bool getTotalProbability()
        """
        ...
    
    def run(self, in_0: EmpiricalFormula ) -> IsotopeDistribution:
        """
        Cython signature: IsotopeDistribution run(EmpiricalFormula)
        """
        ... 


class Gradient:
    """
    Cython implementation of _Gradient

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Gradient.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Gradient()
        Representation of a HPLC gradient
        """
        ...
    
    @overload
    def __init__(self, in_0: Gradient ) -> None:
        """
        Cython signature: void Gradient(Gradient &)
        """
        ...
    
    def addEluent(self, eluent: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void addEluent(String eluent)
        Adds an eluent at the end of the eluent array
        """
        ...
    
    def clearEluents(self) -> None:
        """
        Cython signature: void clearEluents()
        Removes all eluents
        """
        ...
    
    def getEluents(self) -> List[bytes]:
        """
        Cython signature: libcpp_vector[String] getEluents()
        Returns a reference to the list of eluents
        """
        ...
    
    def addTimepoint(self, timepoint: int ) -> None:
        """
        Cython signature: void addTimepoint(int timepoint)
        Adds a timepoint at the end of the timepoint array
        """
        ...
    
    def clearTimepoints(self) -> None:
        """
        Cython signature: void clearTimepoints()
        Removes all timepoints
        """
        ...
    
    def getTimepoints(self) -> List[int]:
        """
        Cython signature: libcpp_vector[int] getTimepoints()
        Returns a reference to the list of timepoints
        """
        ...
    
    def setPercentage(self, eluent: Union[bytes, str, String] , timepoint: int , percentage: int ) -> None:
        """
        Cython signature: void setPercentage(String eluent, int timepoint, unsigned int percentage)
        Sets the percentage of 'eluent' at 'timepoint'
        """
        ...
    
    def getPercentage(self, eluent: Union[bytes, str, String] , timepoint: int ) -> int:
        """
        Cython signature: unsigned int getPercentage(String eluent, int timepoint)
        Returns a const reference to the percentages
        """
        ...
    
    def clearPercentages(self) -> None:
        """
        Cython signature: void clearPercentages()
        Sets all percentage values to 0
        """
        ...
    
    def isValid(self) -> bool:
        """
        Cython signature: bool isValid()
        Checks if the percentages of all timepoints add up to 100%
        """
        ... 


class IsotopeDistribution:
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
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IsotopeDistribution()
        """
        ...
    
    @overload
    def __init__(self, in_0: IsotopeDistribution ) -> None:
        """
        Cython signature: void IsotopeDistribution(IsotopeDistribution &)
        """
        ...
    
    def set(self, distribution: List[Peak1D] ) -> None:
        """
        Cython signature: void set(libcpp_vector[Peak1D] & distribution)
        Overwrites the container which holds the distribution using 'distribution'
        """
        ...
    
    def insert(self, mass: float , intensity: float ) -> None:
        """
        Cython signature: void insert(double mass, float intensity)
        """
        ...
    
    def getContainer(self) -> List[Peak1D]:
        """
        Cython signature: libcpp_vector[Peak1D] & getContainer()
        Returns the container which holds the distribution
        """
        ...
    
    def getMax(self) -> float:
        """
        Cython signature: double getMax()
        Returns the maximal weight isotope which is stored in the distribution
        """
        ...
    
    def getMin(self) -> float:
        """
        Cython signature: double getMin()
        Returns the minimal weight isotope which is stored in the distribution
        """
        ...
    
    def getMostAbundant(self) -> Peak1D:
        """
        Cython signature: Peak1D getMostAbundant()
        Returns the most abundant isotope which is stored in the distribution
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: size_t size()
        Returns the size of the distribution which is the number of isotopes in the distribution
        """
        ...
    
    def clear(self) -> None:
        """
        Cython signature: void clear()
        Clears the distribution and resets max isotope to 0
        """
        ...
    
    def renormalize(self) -> None:
        """
        Cython signature: void renormalize()
        Renormalizes the sum of the probabilities of the isotopes to 1
        """
        ...
    
    def trimRight(self, cutoff: float ) -> None:
        """
        Cython signature: void trimRight(double cutoff)
        Trims the right side of the isotope distribution to isotopes with a significant contribution
        """
        ...
    
    def trimLeft(self, cutoff: float ) -> None:
        """
        Cython signature: void trimLeft(double cutoff)
        Trims the left side of the isotope distribution to isotopes with a significant contribution
        """
        ...
    
    def merge(self, in_0: float , in_1: float ) -> None:
        """
        Cython signature: void merge(double, double)
        Merges distributions of arbitrary data points with constant defined resolution
        """
        ...
    
    def resize(self, size: int ) -> None:
        """
        Cython signature: void resize(unsigned int size)
        Resizes distribution container
        """
        ...
    
    def trimIntensities(self, cutoff: float ) -> None:
        """
        Cython signature: void trimIntensities(double cutoff)
        Remove intensities below the cutoff
        """
        ...
    
    def sortByIntensity(self) -> None:
        """
        Cython signature: void sortByIntensity()
        Sort isotope distribution by intensity
        """
        ...
    
    def sortByMass(self) -> None:
        """
        Cython signature: void sortByMass()
        Sort isotope distribution by mass
        """
        ...
    
    def averageMass(self) -> float:
        """
        Cython signature: double averageMass()
        Compute average mass of isotope distribution (weighted average of all isotopes)
        """
        ...
    
    def __iter__(self) -> Peak1D:
       ...
    Sorted : __Sorted 


class MSDataSqlConsumer:
    """
    Cython implementation of _MSDataSqlConsumer

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSDataSqlConsumer.html>`_
    """
    
    @overload
    def __init__(self, filename: Union[bytes, str, String] , run_id: int , buffer_size: int , full_meta: bool , lossy_compression: bool , linear_mass_acc: float ) -> None:
        """
        Cython signature: void MSDataSqlConsumer(String filename, uint64_t run_id, int buffer_size, bool full_meta, bool lossy_compression, double linear_mass_acc)
        """
        ...
    
    @overload
    def __init__(self, in_0: MSDataSqlConsumer ) -> None:
        """
        Cython signature: void MSDataSqlConsumer(MSDataSqlConsumer &)
        """
        ...
    
    def flush(self) -> None:
        """
        Cython signature: void flush()
        Flushes the data for good
        
        After calling this function, no more data is held in the buffer but the
        class is still able to receive new data
        """
        ...
    
    def consumeSpectrum(self, s: MSSpectrum ) -> None:
        """
        Cython signature: void consumeSpectrum(MSSpectrum & s)
        Write a spectrum to the output file
        """
        ...
    
    def consumeChromatogram(self, c: MSChromatogram ) -> None:
        """
        Cython signature: void consumeChromatogram(MSChromatogram & c)
        Write a chromatogram to the output file
        """
        ...
    
    def setExpectedSize(self, expectedSpectra: int , expectedChromatograms: int ) -> None:
        """
        Cython signature: void setExpectedSize(size_t expectedSpectra, size_t expectedChromatograms)
        """
        ...
    
    def setExperimentalSettings(self, exp: ExperimentalSettings ) -> None:
        """
        Cython signature: void setExperimentalSettings(ExperimentalSettings & exp)
        """
        ... 


class MSPFile:
    """
    Cython implementation of _MSPFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MSPFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MSPFile()
        File adapter for MSP files (NIST spectra library)
        """
        ...
    
    @overload
    def __init__(self, in_0: MSPFile ) -> None:
        """
        Cython signature: void MSPFile(MSPFile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , exp: AnnotatedMSRun ) -> None:
        """
        Cython signature: void store(String filename, AnnotatedMSRun & exp)
        Stores a map in a MSPFile file
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , ids: PeptideIdentificationList , exp: MSExperiment ) -> None:
        """
        Cython signature: void load(String filename, PeptideIdentificationList & ids, MSExperiment & exp)
        Loads a map from a MSPFile file
        
        
        :param exp: PeakMap which contains the spectra after reading
        :param filename: The filename of the experiment
        :param ids: Output parameter which contains the peptide identifications from the spectra annotations
        """
        ... 


class MzTabFile:
    """
    Cython implementation of _MzTabFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MzTabFile.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MzTabFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MzTabFile ) -> None:
        """
        Cython signature: void MzTabFile(MzTabFile &)
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , mz_tab: MzTab ) -> None:
        """
        Cython signature: void store(String filename, MzTab & mz_tab)
        Stores MzTab file
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , mz_tab: MzTab ) -> None:
        """
        Cython signature: void load(String filename, MzTab & mz_tab)
        Loads MzTab file
        """
        ... 


class ParamXMLFile:
    """
    Cython implementation of _ParamXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ParamXMLFile.html>`_

    The file pendant of the Param class used to load and store the param
    datastructure as paramXML
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ParamXMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: ParamXMLFile ) -> None:
        """
        Cython signature: void ParamXMLFile(ParamXMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , param: Param ) -> None:
        """
        Cython signature: void load(String filename, Param & param)
        Read XML file
        
        
        :param filename: The file from where to read the Param object
        :param param: The param object where the read data should be stored
        :raises:
          Exception: FileNotFound is thrown if the file could not be found
        :raises:
          Exception: ParseError is thrown if an error occurs during parsing
        """
        ...
    
    def store(self, filename: Union[bytes, str, String] , param: Param ) -> None:
        """
        Cython signature: void store(String filename, Param & param)
        Write XML file
        
        
        :param filename: The filename where the param data structure should be stored
        :param param: The Param class that should be stored in the file
        """
        ... 


class ResidueDB:
    """
    Cython implementation of _ResidueDB

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ResidueDB.html>`_
    """
    
    def getNumberOfResidues(self) -> int:
        """
        Cython signature: size_t getNumberOfResidues()
        Returns the number of residues stored
        """
        ...
    
    def getNumberOfModifiedResidues(self) -> int:
        """
        Cython signature: size_t getNumberOfModifiedResidues()
        Returns the number of modified residues stored
        """
        ...
    
    def getResidue(self, name: Union[bytes, str, String] ) -> Residue:
        """
        Cython signature: const Residue * getResidue(const String & name)
        Returns a pointer to the residue with name, 3 letter code or 1 letter code name
        """
        ...
    
    @overload
    def getModifiedResidue(self, name: Union[bytes, str, String] ) -> Residue:
        """
        Cython signature: const Residue * getModifiedResidue(const String & name)
        Returns a pointer to a modified residue given a modification name
        """
        ...
    
    @overload
    def getModifiedResidue(self, residue: Residue , name: Union[bytes, str, String] ) -> Residue:
        """
        Cython signature: const Residue * getModifiedResidue(Residue * residue, const String & name)
        Returns a pointer to a modified residue given a residue and a modification name
        """
        ...
    
    def getResidues(self, residue_set: Union[bytes, str, String] ) -> Set[Residue]:
        """
        Cython signature: libcpp_set[const Residue *] getResidues(const String & residue_set)
        Returns a set of all residues stored in this residue db
        """
        ...
    
    def getResidueSets(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getResidueSets()
        Returns all residue sets that are registered which this instance
        """
        ...
    
    def hasResidue(self, name: Union[bytes, str, String] ) -> bool:
        """
        Cython signature: bool hasResidue(const String & name)
        Returns true if the db contains a residue with the given name
        """
        ... 


class SpectrumAccessQuadMZTransforming:
    """
    Cython implementation of _SpectrumAccessQuadMZTransforming

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1SpectrumAccessQuadMZTransforming.html>`_
      -- Inherits from ['SpectrumAccessTransforming']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void SpectrumAccessQuadMZTransforming()
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessQuadMZTransforming ) -> None:
        """
        Cython signature: void SpectrumAccessQuadMZTransforming(SpectrumAccessQuadMZTransforming &)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMS , a: float , b: float , c: float , ppm: bool ) -> None:
        """
        Cython signature: void SpectrumAccessQuadMZTransforming(shared_ptr[SpectrumAccessOpenMS], double a, double b, double c, bool ppm)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMSCached , a: float , b: float , c: float , ppm: bool ) -> None:
        """
        Cython signature: void SpectrumAccessQuadMZTransforming(shared_ptr[SpectrumAccessOpenMSCached], double a, double b, double c, bool ppm)
        """
        ...
    
    @overload
    def __init__(self, in_0: SpectrumAccessOpenMSInMemory , a: float , b: float , c: float , ppm: bool ) -> None:
        """
        Cython signature: void SpectrumAccessQuadMZTransforming(shared_ptr[SpectrumAccessOpenMSInMemory], double a, double b, double c, bool ppm)
        """
        ...
    
    def getSpectrumById(self, id_: int ) -> OSSpectrum:
        """
        Cython signature: shared_ptr[OSSpectrum] getSpectrumById(int id_)
        Returns a pointer to a spectrum at the given string id
        """
        ...
    
    def getSpectraByRT(self, RT: float , deltaRT: float ) -> List[int]:
        """
        Cython signature: libcpp_vector[size_t] getSpectraByRT(double RT, double deltaRT)
        Returns a vector of ids of spectra that are within RT +/- deltaRT
        """
        ...
    
    def getNrSpectra(self) -> int:
        """
        Cython signature: size_t getNrSpectra()
        Returns the number of spectra available
        """
        ...
    
    def getChromatogramById(self, id_: int ) -> OSChromatogram:
        """
        Cython signature: shared_ptr[OSChromatogram] getChromatogramById(int id_)
        Returns a pointer to a chromatogram at the given id
        """
        ...
    
    def getNrChromatograms(self) -> int:
        """
        Cython signature: size_t getNrChromatograms()
        Returns the number of chromatograms available
        """
        ...
    
    def getChromatogramNativeID(self, id_: int ) -> str:
        """
        Cython signature: libcpp_utf8_output_string getChromatogramNativeID(int id_)
        """
        ... 


class StablePairFinder:
    """
    Cython implementation of _StablePairFinder

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1StablePairFinder.html>`_
      -- Inherits from ['BaseGroupFinder']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void StablePairFinder()
        """
        ...
    
    def run(self, input_maps: List[ConsensusMap] , result_map: ConsensusMap ) -> None:
        """
        Cython signature: void run(libcpp_vector[ConsensusMap] & input_maps, ConsensusMap & result_map)
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


class __Sorted:
    None
    INTENSITY : int
    MASS : int
    UNDEFINED : int

    def getMapping(self) -> Dict[int, str]:
       ... 

