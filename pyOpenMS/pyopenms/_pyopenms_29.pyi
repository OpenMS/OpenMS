from __future__ import annotations
from typing import overload, Any, List, Dict, Tuple, Set, Sequence, Union
from pyopenms import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
import numpy as _np

from enum import Enum as _PyEnum


def __static_MZTrafoModel_enumToName(mt: int ) -> bytes:
    """
    Cython signature: libcpp_string enumToName(MZTrafoModel_MODELTYPE mt)
    """
    ...

def __static_MZTrafoModel_findNearest(tms: List[MZTrafoModel] , rt: float ) -> int:
    """
    Cython signature: size_t findNearest(libcpp_vector[MZTrafoModel] & tms, double rt)
    """
    ...

def __static_ExperimentalDesign_fromConsensusMap(c: ConsensusMap ) -> ExperimentalDesign:
    """
    Cython signature: ExperimentalDesign fromConsensusMap(ConsensusMap c)
    """
    ...

def __static_ExperimentalDesign_fromFeatureMap(f: FeatureMap ) -> ExperimentalDesign:
    """
    Cython signature: ExperimentalDesign fromFeatureMap(FeatureMap f)
    """
    ...

def __static_ExperimentalDesign_fromIdentifications(proteins: List[ProteinIdentification] ) -> ExperimentalDesign:
    """
    Cython signature: ExperimentalDesign fromIdentifications(libcpp_vector[ProteinIdentification] & proteins)
    """
    ...

def __static_MZTrafoModel_isValidModel(trafo: MZTrafoModel ) -> bool:
    """
    Cython signature: bool isValidModel(MZTrafoModel & trafo)
    """
    ...

def __static_MZTrafoModel_nameToEnum(name: bytes ) -> int:
    """
    Cython signature: MZTrafoModel_MODELTYPE nameToEnum(libcpp_string name)
    """
    ...

def __static_MZTrafoModel_setCoefficientLimits(offset: float , scale: float , power: float ) -> None:
    """
    Cython signature: void setCoefficientLimits(double offset, double scale, double power)
    """
    ...

def __static_MZTrafoModel_setRANSACParams(p: RANSACParam ) -> None:
    """
    Cython signature: void setRANSACParams(RANSACParam p)
    """
    ...


class BSpline2d:
    """
    Cython implementation of _BSpline2d

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1BSpline2d.html>`_
    """
    
    def __init__(self, x: List[float] , y: List[float] , wave_length: float , boundary_condition: int , num_nodes: int ) -> None:
        """
        Cython signature: void BSpline2d(libcpp_vector[double] x, libcpp_vector[double] y, double wave_length, BoundaryCondition boundary_condition, size_t num_nodes)
        """
        ...
    
    def solve(self, y: List[float] ) -> bool:
        """
        Cython signature: bool solve(libcpp_vector[double] y)
        Solve the spline curve for a new set of y values. Returns false if the solution fails
        """
        ...
    
    def eval(self, x: float ) -> float:
        """
        Cython signature: double eval(double x)
        Returns the evaluation of the smoothed curve at a particular x value. If current state is not ok(), returns zero
        """
        ...
    
    def derivative(self, x: float ) -> float:
        """
        Cython signature: double derivative(double x)
        Returns the first derivative of the spline curve at the given position x. Returns zero if the current state is not ok()
        """
        ...
    
    def ok(self) -> bool:
        """
        Cython signature: bool ok()
        Returns whether the spline fit was successful
        """
        ...
    
    def debug(self, enable: bool ) -> None:
        """
        Cython signature: void debug(bool enable)
        Enable or disable debug messages from the B-spline library
        """
        ... 


class CVReference:
    """
    Cython implementation of _CVReference

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1CVReference.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void CVReference()
        """
        ...
    
    @overload
    def __init__(self, in_0: CVReference ) -> None:
        """
        Cython signature: void CVReference(CVReference &)
        """
        ...
    
    def setName(self, name: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(const String & name)
        Sets the name of the CV reference
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the CV reference
        """
        ...
    
    def setIdentifier(self, identifier: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setIdentifier(const String & identifier)
        Sets the CV identifier which is referenced
        """
        ...
    
    def getIdentifier(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getIdentifier()
        Returns the CV identifier which is referenced
        """
        ...
    
    def __richcmp__(self, other: CVReference, op: int) -> Any:
        ... 


class ExperimentalDesign:
    """
    Cython implementation of _ExperimentalDesign

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign.html>`_

    Representation of an experimental design in OpenMS. Instances can be loaded with the ExperimentalDesignFile class
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ExperimentalDesign()
        """
        ...
    
    @overload
    def __init__(self, in_0: ExperimentalDesign ) -> None:
        """
        Cython signature: void ExperimentalDesign(ExperimentalDesign &)
        """
        ...
    
    def getMSFileSection(self) -> List[ExperimentalDesign_MSFileSectionEntry]:
        """
        Cython signature: libcpp_vector[ExperimentalDesign_MSFileSectionEntry] getMSFileSection()
        """
        ...
    
    def setMSFileSection(self, msfile_section: List[ExperimentalDesign_MSFileSectionEntry] ) -> None:
        """
        Cython signature: void setMSFileSection(libcpp_vector[ExperimentalDesign_MSFileSectionEntry] msfile_section)
        """
        ...
    
    def getSampleSection(self) -> ExperimentalDesign_SampleSection:
        """
        Cython signature: ExperimentalDesign_SampleSection getSampleSection()
        Returns the Sample Section of the experimental design file
        """
        ...
    
    def setSampleSection(self, sample_section: ExperimentalDesign_SampleSection ) -> None:
        """
        Cython signature: void setSampleSection(ExperimentalDesign_SampleSection sample_section)
        Sets the Sample Section of the experimental design file
        """
        ...
    
    def getNumberOfSamples(self) -> int:
        """
        Cython signature: unsigned int getNumberOfSamples()
        Returns the number of samples measured (= highest sample index)
        """
        ...
    
    def getNumberOfFractions(self) -> int:
        """
        Cython signature: unsigned int getNumberOfFractions()
        Returns the number of fractions (= highest fraction index)
        """
        ...
    
    def getNumberOfLabels(self) -> int:
        """
        Cython signature: unsigned int getNumberOfLabels()
        Returns the number of labels per file
        """
        ...
    
    def getNumberOfMSFiles(self) -> int:
        """
        Cython signature: unsigned int getNumberOfMSFiles()
        Returns the number of MS files (= fractions * fraction_groups)
        """
        ...
    
    def getNumberOfFractionGroups(self) -> int:
        """
        Cython signature: unsigned int getNumberOfFractionGroups()
        Allows to group fraction ids and source files. Return the number of fraction_groups
        """
        ...
    
    def getSample(self, fraction_group: int , label: int ) -> int:
        """
        Cython signature: unsigned int getSample(unsigned int fraction_group, unsigned int label)
        Returns sample index (depends on fraction_group and label)
        """
        ...
    
    def isFractionated(self) -> bool:
        """
        Cython signature: bool isFractionated()
        Returns whether at least one fraction_group in this experimental design is fractionated
        """
        ...
    
    def sameNrOfMSFilesPerFraction(self) -> bool:
        """
        Cython signature: bool sameNrOfMSFilesPerFraction()
        Returns if each fraction number is associated with the same number of fraction_group
        """
        ...
    
    fromConsensusMap: __static_ExperimentalDesign_fromConsensusMap
    
    fromFeatureMap: __static_ExperimentalDesign_fromFeatureMap
    
    fromIdentifications: __static_ExperimentalDesign_fromIdentifications 


class ExperimentalDesign_MSFileSectionEntry:
    """
    Cython implementation of _ExperimentalDesign_MSFileSectionEntry

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign_MSFileSectionEntry.html>`_
    """
    
    path: bytes
    
    fraction_group: int
    
    fraction: int
    
    label: int
    
    sample: int
    
    def __init__(self) -> None:
        """
        Cython signature: void ExperimentalDesign_MSFileSectionEntry()
        """
        ... 


class ExperimentalDesign_SampleSection:
    """
    Cython implementation of _ExperimentalDesign_SampleSection

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1ExperimentalDesign_SampleSection.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void ExperimentalDesign_SampleSection()
        """
        ...
    
    @overload
    def __init__(self, in_0: ExperimentalDesign_SampleSection ) -> None:
        """
        Cython signature: void ExperimentalDesign_SampleSection(ExperimentalDesign_SampleSection)
        """
        ...
    
    def getSamples(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getSamples()
        Returns a set of all samples that are present in the sample section
        """
        ...
    
    def getFactors(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] getFactors()
        Returns a set of all factors (column names) that were defined for the sample section
        """
        ...
    
    def hasSample(self, sample: int ) -> bool:
        """
        Cython signature: bool hasSample(unsigned int sample)
        Checks whether sample section has row for a sample number
        """
        ...
    
    def hasFactor(self, factor: String ) -> bool:
        """
        Cython signature: bool hasFactor(String & factor)
        Checks whether Sample Section has a specific factor (i.e. column name)
        """
        ...
    
    def getFactorValue(self, sample: int , factor: String ) -> Union[bytes, str, String]:
        """
        Cython signature: String getFactorValue(unsigned int sample, String & factor)
        Returns value of factor for given sample and factor name
        """
        ... 


class IMSIsotopeDistribution:
    """
    Cython implementation of _IMSIsotopeDistribution

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSIsotopeDistribution_1_1IMSIsotopeDistribution.html>`_

    Represents a distribution of isotopes restricted to the first K elements
    
    Represents a distribution of isotopes of chemical elements as a list
    of peaks each as a pair of mass and abundance. 'IsotopeDistribution'
    unlike 'IsotopeSpecies' has one abundance per a nominal mass.
    Here is an example in the format (mass; abundance %)
    for molecule H2O (values are taken randomly):
    
    - IsotopeDistribution
        (18.00221; 99.03 %)
        (19.00334; 0.8 %)
        (20.00476; 0.17 %)
    
    - IsotopeSpecies
        (18.00197; 98.012 %)
        (18.00989; 1.018 %)
        (19.00312; 0.683 %)
        (19.00531; 0.117 %)
        (20.00413; 0.134 %)
        (20.00831; 0.036 %)
    
    To the sake of faster computations distribution is restricted
    to the first K elements, where K can be set by adjusting size
    'SIZE' of distribution. @note For the elements most abundant in
    living beings (CHNOPS) this restriction is negligible, since abundances
    decrease dramatically in isotopes order and are usually of no interest
    starting from +10 isotope.
    
    'IsotopeDistribution' implements folding with other distribution using an
    algorithm described in details in paper:
    Boecker et al. "Decomposing metabolic isotope patterns" WABI 2006. doi: 10.1007/11851561_2
    
    Folding with itself is done using Russian Multiplication Scheme
    """
    
    ABUNDANCES_SUM_ERROR: float
    
    SIZE: int
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution()
        """
        ...
    
    @overload
    def __init__(self, in_0: IMSIsotopeDistribution ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution(IMSIsotopeDistribution &)
        """
        ...
    
    @overload
    def __init__(self, nominalMass: int ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution(unsigned int nominalMass)
        """
        ...
    
    @overload
    def __init__(self, mass: float ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution(double mass)
        """
        ...
    
    @overload
    def __init__(self, peaks: List[IMSIsotopeDistribution_Peak] , nominalMass: int ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution(libcpp_vector[IMSIsotopeDistribution_Peak] & peaks, unsigned int nominalMass)
        """
        ...
    
    def size(self) -> int:
        """
        Cython signature: int size()
        """
        ...
    
    def getMass(self, i: int ) -> float:
        """
        Cython signature: double getMass(int i)
        """
        ...
    
    def getAbundance(self, i: int ) -> float:
        """
        Cython signature: double getAbundance(int i)
        """
        ...
    
    def getAverageMass(self) -> float:
        """
        Cython signature: double getAverageMass()
        """
        ...
    
    def getNominalMass(self) -> int:
        """
        Cython signature: unsigned int getNominalMass()
        """
        ...
    
    def setNominalMass(self, nominalMass: int ) -> None:
        """
        Cython signature: void setNominalMass(unsigned int nominalMass)
        """
        ...
    
    def getMasses(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getMasses()
        Gets a mass of isotope 'i'
        """
        ...
    
    def getAbundances(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getAbundances()
        Gets an abundance of isotope 'i'
        """
        ...
    
    def normalize(self) -> None:
        """
        Cython signature: void normalize()
        Normalizes distribution, i.e. scaling abundances to be summed up to 1 with an error
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Returns true if the distribution has no peaks, false - otherwise
        """
        ...
    
    def __richcmp__(self, other: IMSIsotopeDistribution, op: int) -> Any:
        ... 


class IMSIsotopeDistribution_Peak:
    """
    Cython implementation of _IMSIsotopeDistribution_Peak

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::ims::IMSIsotopeDistribution_1_1IMSIsotopeDistribution_Peak.html>`_
    """
    
    mass: float
    
    abundance: float
    
    def __init__(self, mass: float , abundance: float ) -> None:
        """
        Cython signature: void IMSIsotopeDistribution_Peak(double mass, double abundance)
        """
        ...
    
    def __richcmp__(self, other: IMSIsotopeDistribution_Peak, op: int) -> Any:
        ... 


class LinearInterpolation:
    """
    Cython implementation of _LinearInterpolation[double,double]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1LinearInterpolation[double,double].html>`_

    Provides access to linearly interpolated values (and
    derivatives) from discrete data points.  Values beyond the given range
    of data points are implicitly taken as zero.
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void LinearInterpolation()
        """
        ...
    
    @overload
    def __init__(self, in_0: LinearInterpolation ) -> None:
        """
        Cython signature: void LinearInterpolation(LinearInterpolation &)
        """
        ...
    
    @overload
    def __init__(self, scale: float , offset: float ) -> None:
        """
        Cython signature: void LinearInterpolation(double scale, double offset)
        """
        ...
    
    def value(self, arg_pos: float ) -> float:
        """
        Cython signature: double value(double arg_pos)
        Returns the interpolated value
        """
        ...
    
    def addValue(self, arg_pos: float , arg_value: float ) -> None:
        """
        Cython signature: void addValue(double arg_pos, double arg_value)
        Performs linear resampling. The `arg_value` is split up and added to the data points around `arg_pos`
        """
        ...
    
    def derivative(self, arg_pos: float ) -> float:
        """
        Cython signature: double derivative(double arg_pos)
        Returns the interpolated derivative
        """
        ...
    
    def getData(self) -> List[float]:
        """
        Cython signature: libcpp_vector[double] getData()
        Returns the internal random access container from which interpolated values are being sampled
        """
        ...
    
    def setData(self, data: List[float] ) -> None:
        """
        Cython signature: void setData(libcpp_vector[double] & data)
        Assigns data to the internal random access container from which interpolated values are being sampled
        """
        ...
    
    def empty(self) -> bool:
        """
        Cython signature: bool empty()
        Returns `true` if getData() is empty
        """
        ...
    
    def key2index(self, pos: float ) -> float:
        """
        Cython signature: double key2index(double pos)
        The transformation from "outside" to "inside" coordinates
        """
        ...
    
    def index2key(self, pos: float ) -> float:
        """
        Cython signature: double index2key(double pos)
        The transformation from "inside" to "outside" coordinates
        """
        ...
    
    def getScale(self) -> float:
        """
        Cython signature: double getScale()
        "Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        """
        ...
    
    def setScale(self, scale: float ) -> None:
        """
        Cython signature: void setScale(double & scale)
        "Scale" is the difference (in "outside" units) between consecutive entries in "Data"
        """
        ...
    
    def getOffset(self) -> float:
        """
        Cython signature: double getOffset()
        "Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        """
        ...
    
    def setOffset(self, offset: float ) -> None:
        """
        Cython signature: void setOffset(double & offset)
        "Offset" is the point (in "outside" units) which corresponds to "Data[0]"
        """
        ...
    
    @overload
    def setMapping(self, scale: float , inside: float , outside: float ) -> None:
        """
        Cython signature: void setMapping(double & scale, double & inside, double & outside)
        """
        ...
    
    @overload
    def setMapping(self, inside_low: float , outside_low: float , inside_high: float , outside_high: float ) -> None:
        """
        Cython signature: void setMapping(double & inside_low, double & outside_low, double & inside_high, double & outside_high)
        """
        ...
    
    def getInsideReferencePoint(self) -> float:
        """
        Cython signature: double getInsideReferencePoint()
        """
        ...
    
    def getOutsideReferencePoint(self) -> float:
        """
        Cython signature: double getOutsideReferencePoint()
        """
        ...
    
    def supportMin(self) -> float:
        """
        Cython signature: double supportMin()
        """
        ...
    
    def supportMax(self) -> float:
        """
        Cython signature: double supportMax()
        """
        ... 


class LinearResamplerAlign:
    """
    Cython implementation of _LinearResamplerAlign

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1LinearResamplerAlign.html>`_
      -- Inherits from ['LinearResampler']
    """
    
    def __init__(self, in_0: LinearResamplerAlign ) -> None:
        """
        Cython signature: void LinearResamplerAlign(LinearResamplerAlign &)
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


class MRMFeatureFinderScoring:
    """
    Cython implementation of _MRMFeatureFinderScoring

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MRMFeatureFinderScoring.html>`_
      -- Inherits from ['DefaultParamHandler', 'ProgressLogger']
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void MRMFeatureFinderScoring()
        """
        ...
    
    def pickExperiment(self, chromatograms: MSExperiment , output: FeatureMap , transition_exp_: TargetedExperiment , trafo: TransformationDescription , swath_map: MSExperiment ) -> None:
        """
        Cython signature: void pickExperiment(MSExperiment & chromatograms, FeatureMap & output, TargetedExperiment & transition_exp_, TransformationDescription trafo, MSExperiment & swath_map)
        Pick features in one experiment containing chromatogram
        
        Function for for wrapping in Python, only uses OpenMS datastructures and does not return the map
        
        
        :param chromatograms: The input chromatograms
        :param output: The output features with corresponding scores
        :param transition_exp: The transition list describing the experiment
        :param trafo: Optional transformation of the experimental retention time to the normalized retention time space used in the transition list
        :param swath_map: Optional SWATH-MS (DIA) map corresponding from which the chromatograms were extracted
        """
        ...
    
    def setStrictFlag(self, flag: bool ) -> None:
        """
        Cython signature: void setStrictFlag(bool flag)
        """
        ...
    
    @overload
    def setMS1Map(self, ms1_map: SpectrumAccessOpenMS ) -> None:
        """
        Cython signature: void setMS1Map(shared_ptr[SpectrumAccessOpenMS] ms1_map)
        """
        ...
    
    @overload
    def setMS1Map(self, ms1_map: SpectrumAccessOpenMSCached ) -> None:
        """
        Cython signature: void setMS1Map(shared_ptr[SpectrumAccessOpenMSCached] ms1_map)
        """
        ...
    
    def scorePeakgroups(self, transition_group: LightMRMTransitionGroupCP , trafo: TransformationDescription , swath_maps: List[SwathMap] , output: FeatureMap , ms1only: bool ) -> None:
        """
        Cython signature: void scorePeakgroups(LightMRMTransitionGroupCP transition_group, TransformationDescription trafo, libcpp_vector[SwathMap] swath_maps, FeatureMap & output, bool ms1only)
        Score all peak groups of a transition group
        
        Iterate through all features found along the chromatograms of the transition group and score each one individually
        
        
        :param transition_group: The MRMTransitionGroup to be scored (input)
        :param trafo: Optional transformation of the experimental retention time
            to the normalized retention time space used in thetransition list
        :param swath_maps: Optional SWATH-MS (DIA) map corresponding from which
            the chromatograms were extracted. Use empty map if no data is available
        :param output: The output features with corresponding scores (the found
            features will be added to this FeatureMap)
        :param ms1only: Whether to only do MS1 scoring and skip all MS2 scoring
        """
        ...
    
    def prepareProteinPeptideMaps_(self, transition_exp: LightTargetedExperiment ) -> None:
        """
        Cython signature: void prepareProteinPeptideMaps_(LightTargetedExperiment & transition_exp)
        Prepares the internal mappings of peptides and proteins
        
        Calling this method _is_ required before calling scorePeakgroups
        
        
        :param transition_exp: The transition list describing the experiment
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


class MZTrafoModel:
    """
    Cython implementation of _MZTrafoModel

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MZTrafoModel.html>`_

    Create and apply models of a mass recalibration function
    
    The input is a list of calibration points (ideally spanning a wide m/z range to prevent extrapolation when applying to model)
    
    Models (LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED) can be trained using CalData points (or a subset of them)
    Calibration points can have different retention time points, and a model should be build such that it captures
    the local (in time) decalibration of the instrument, i.e. choose appropriate time windows along RT to calibrate the
    spectra in this RT region
    From the available calibrant data, a model is build. Later, any uncalibrated m/z value can be fed to the model, to obtain
    a calibrated m/z
    
    The input domain can either be absolute mass differences in [Th], or relative differences in [ppm]
    The models are build based on this input
    
    Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MZTrafoModel()
        """
        ...
    
    @overload
    def __init__(self, in_0: MZTrafoModel ) -> None:
        """
        Cython signature: void MZTrafoModel(MZTrafoModel &)
        """
        ...
    
    @overload
    def __init__(self, in_0: bool ) -> None:
        """
        Cython signature: void MZTrafoModel(bool)
        """
        ...
    
    def isTrained(self) -> bool:
        """
        Cython signature: bool isTrained()
        Returns true if the model have coefficients (i.e. was trained successfully)
        """
        ...
    
    def getRT(self) -> float:
        """
        Cython signature: double getRT()
        Get RT associated with the model (training region)
        """
        ...
    
    def predict(self, mz: float ) -> float:
        """
        Cython signature: double predict(double mz)
        Apply the model to an uncalibrated m/z value
        
        Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!
        
        Applies the function y = intercept + slope*mz + power*mz^2
        and returns y
        
        
        :param mz: The uncalibrated m/z value
        :return: The calibrated m/z value
        """
        ...
    
    @overload
    def train(self, cd: CalibrationData , md: int , use_RANSAC: bool , rt_left: float , rt_right: float ) -> bool:
        """
        Cython signature: bool train(CalibrationData cd, MZTrafoModel_MODELTYPE md, bool use_RANSAC, double rt_left, double rt_right)
        Train a model using calibrant data
        
        If the CalibrationData was created using peak groups (usually corresponding to mass traces),
        the median for each group is used as a group representative. This
        is more robust, and reduces the number of data points drastically, i.e. one value per group
        
        Internally, these steps take place:
        - apply RT filter
        - [compute median per group] (only if groups were given in 'cd')
        - set Model's rt position
        - call train() (see overloaded method)
        
        
        :param cd: List of calibrants
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :param rt_left: Filter 'cd' by RT; all calibrants with RT < 'rt_left' are removed
        :param rt_right: Filter 'cd' by RT; all calibrants with RT > 'rt_right' are removed
        :return: True if model was build, false otherwise
        """
        ...
    
    @overload
    def train(self, error_mz: List[float] , theo_mz: List[float] , weights: List[float] , md: int , use_RANSAC: bool ) -> bool:
        """
        Cython signature: bool train(libcpp_vector[double] error_mz, libcpp_vector[double] theo_mz, libcpp_vector[double] weights, MZTrafoModel_MODELTYPE md, bool use_RANSAC)
        Train a model using calibrant data
        
        Given theoretical and observed mass values (and corresponding weights),
        a model (linear, quadratic, ...) is build
        Outlier removal is applied before
        The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm]
        The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
        rather on how 'predict()' works internally
        
        Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
        
        Internally, these steps take place:
        - [apply RANSAC] (depending on 'use_RANSAC')
        - build model and store its parameters internally
        
        
        :param error_mz: Observed Mass error (in ppm or Th)
        :param theo_mz: Theoretical m/z values, corresponding to 'error_mz'
        :param weights: For weighted models only: weight of calibrants; ignored otherwise
        :param md: Type of model (linear, quadratic, ...)
        :param use_RANSAC: Remove outliers before computing the model?
        :return: True if model was build, false otherwise
        """
        ...
    
    def getCoefficients(self, intercept: float , slope: float , power: float ) -> None:
        """
        Cython signature: void getCoefficients(double & intercept, double & slope, double & power)
        Get model coefficients
        
        Parameters will be filled with internal model parameters
        The model must be trained before; Exception is thrown otherwise!
        
        
        :param intercept: The intercept
        :param slope: The slope
        :param power: The coefficient for x*x (will be 0 for linear models)
        """
        ...
    
    @overload
    def setCoefficients(self, in_0: MZTrafoModel ) -> None:
        """
        Cython signature: void setCoefficients(MZTrafoModel)
        Copy model coefficients from another model
        """
        ...
    
    @overload
    def setCoefficients(self, in_0: float , in_1: float , in_2: float ) -> None:
        """
        Cython signature: void setCoefficients(double, double, double)
        Manually set model coefficients
        
        Can be used instead of train(), so manually set coefficients
        It must be exactly three values. If you want a linear model, set 'power' to zero
        If you want a constant model, set slope to zero in addition
        
        
        :param intercept: The offset
        :param slope: The slope
        :param power: The x*x coefficient (for quadratic models)
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    enumToName: __static_MZTrafoModel_enumToName
    
    findNearest: __static_MZTrafoModel_findNearest
    
    isValidModel: __static_MZTrafoModel_isValidModel
    
    nameToEnum: __static_MZTrafoModel_nameToEnum
    
    setCoefficientLimits: __static_MZTrafoModel_setCoefficientLimits
    
    setRANSACParams: __static_MZTrafoModel_setRANSACParams 


class MascotXMLFile:
    """
    Cython implementation of _MascotXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1MascotXMLFile.html>`_
      -- Inherits from ['XMLFile']
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void MascotXMLFile()
        """
        ...
    
    @overload
    def __init__(self, in_0: MascotXMLFile ) -> None:
        """
        Cython signature: void MascotXMLFile(MascotXMLFile &)
        """
        ...
    
    def load(self, filename: Union[bytes, str, String] , protein_identification: ProteinIdentification , id_data: PeptideIdentificationList , rt_mapping: SpectrumMetaDataLookup ) -> None:
        """
        Cython signature: void load(const String & filename, ProteinIdentification & protein_identification, PeptideIdentificationList & id_data, SpectrumMetaDataLookup & rt_mapping)
        Loads data from a Mascot XML file
        
        
        :param filename: The file to be loaded
        :param protein_identification: Protein identifications belonging to the whole experiment
        :param id_data: The identifications with m/z and RT
        :param lookup: Helper object for looking up spectrum meta data
        :raises:
          Exception: FileNotFound is thrown if the file does not exists
        :raises:
          Exception: ParseError is thrown if the file does not suit to the standard
        """
        ...
    
    def initializeLookup(self, lookup: SpectrumMetaDataLookup , experiment: MSExperiment , scan_regex: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void initializeLookup(SpectrumMetaDataLookup & lookup, MSExperiment & experiment, const String & scan_regex)
        Initializes a helper object for looking up spectrum meta data (RT, m/z)
        
        
        :param lookup: Helper object to initialize
        :param experiment: Experiment containing the spectra
        :param scan_regex: Optional regular expression for extracting information from references to spectra
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Return the version of the schema
        """
        ... 


class OSWFile:
    """
    Cython implementation of _OSWFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1OSWFile.html>`_

    This class serves for reading in and writing OpenSWATH OSW files
    
    See OpenSwathOSWWriter for more functionality
    
    The reader and writer returns data in a format suitable for PercolatorAdapter.
    OSW files have a flexible data structure. They contain all peptide query
    parameters of TraML/PQP files with the detected and quantified features of
    OpenSwathWorkflow (feature, feature_ms1, feature_ms2 & feature_transition)
    
    The OSWFile reader extracts the feature information from the OSW file for
    each level (MS1, MS2 & transition) separately and generates Percolator input
    files. For each of the three Percolator reports, OSWFile writer adds a table
    (score_ms1, score_ms2, score_transition) with the respective confidence metrics.
    These tables can be mapped to the corresponding feature tables, are very similar
    to PyProphet results and can thus be used interchangeably
    """
    
    @overload
    def __init__(self, filename: Union[bytes, str] ) -> None:
        """
        Cython signature: void OSWFile(const libcpp_utf8_string filename)
        """
        ...
    
    @overload
    def __init__(self, in_0: OSWFile ) -> None:
        """
        Cython signature: void OSWFile(OSWFile &)
        """
        ... 


class PeptideHit:
    """
    Cython implementation of _PeptideHit

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit.html>`_
      -- Inherits from ['MetaInfoInterface']

    Representation of a peptide hit
    
    It contains the fields score, score_type, rank, and sequence
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideHit()
        """
        ...
    
    @overload
    def __init__(self, score: float , rank: int , charge: int , sequence: AASequence ) -> None:
        """
        Cython signature: void PeptideHit(double score, unsigned int rank, int charge, AASequence sequence)
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideHit ) -> None:
        """
        Cython signature: void PeptideHit(PeptideHit &)
        """
        ...
    
    def getScore(self) -> float:
        """
        Cython signature: float getScore()
        Returns the PSM score
        """
        ...
    
    def getRank(self) -> int:
        """
        Cython signature: unsigned int getRank()
        Returns the PSM rank
        """
        ...
    
    def getSequence(self) -> AASequence:
        """
        Cython signature: AASequence getSequence()
        Returns the peptide sequence without trailing or following spaces
        """
        ...
    
    def getCharge(self) -> int:
        """
        Cython signature: int getCharge()
        Returns the charge of the peptide
        """
        ...
    
    def getPeptideEvidences(self) -> List[PeptideEvidence]:
        """
        Cython signature: libcpp_vector[PeptideEvidence] getPeptideEvidences()
        Returns information on peptides (potentially) identified by this PSM
        """
        ...
    
    def setPeptideEvidences(self, in_0: List[PeptideEvidence] ) -> None:
        """
        Cython signature: void setPeptideEvidences(libcpp_vector[PeptideEvidence])
        Sets information on peptides (potentially) identified by this PSM
        """
        ...
    
    def addPeptideEvidence(self, in_0: PeptideEvidence ) -> None:
        """
        Cython signature: void addPeptideEvidence(PeptideEvidence)
        Adds information on a peptide that is (potentially) identified by this PSM
        """
        ...
    
    def extractProteinAccessionsSet(self) -> Set[bytes]:
        """
        Cython signature: libcpp_set[String] extractProteinAccessionsSet()
        Extracts the set of non-empty protein accessions from peptide evidences
        """
        ...
    
    def setAnalysisResults(self, aresult: List[PeptideHit_AnalysisResult] ) -> None:
        """
        Cython signature: void setAnalysisResults(libcpp_vector[PeptideHit_AnalysisResult] aresult)
        Sets information on (search engine) sub scores associated with this PSM
        """
        ...
    
    def addAnalysisResults(self, aresult: PeptideHit_AnalysisResult ) -> None:
        """
        Cython signature: void addAnalysisResults(PeptideHit_AnalysisResult aresult)
        Add information on (search engine) sub scores associated with this PSM
        """
        ...
    
    def getAnalysisResults(self) -> List[PeptideHit_AnalysisResult]:
        """
        Cython signature: libcpp_vector[PeptideHit_AnalysisResult] getAnalysisResults()
        Returns information on (search engine) sub scores associated with this PSM
        """
        ...
    
    def setPeakAnnotations(self, in_0: List[PeptideHit_PeakAnnotation] ) -> None:
        """
        Cython signature: void setPeakAnnotations(libcpp_vector[PeptideHit_PeakAnnotation])
        Sets the fragment annotations
        """
        ...
    
    def getPeakAnnotations(self) -> List[PeptideHit_PeakAnnotation]:
        """
        Cython signature: libcpp_vector[PeptideHit_PeakAnnotation] getPeakAnnotations()
        Returns the fragment annotations
        """
        ...
    
    def setScore(self, in_0: float ) -> None:
        """
        Cython signature: void setScore(double)
        Sets the PSM score
        """
        ...
    
    def setRank(self, in_0: int ) -> None:
        """
        Cython signature: void setRank(unsigned int)
        Sets the PSM rank
        """
        ...
    
    def setSequence(self, in_0: AASequence ) -> None:
        """
        Cython signature: void setSequence(AASequence)
        Sets the peptide sequence
        """
        ...
    
    def setCharge(self, in_0: int ) -> None:
        """
        Cython signature: void setCharge(int)
        Sets the charge of the peptide
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
    
    def __richcmp__(self, other: PeptideHit, op: int) -> Any:
        ... 


class PeptideHit_AnalysisResult:
    """
    Cython implementation of _PeptideHit_AnalysisResult

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit_AnalysisResult.html>`_
    """
    
    score_type: Union[bytes, str, String]
    
    higher_is_better: bool
    
    main_score: float
    
    sub_scores: Dict[Union[bytes, str, String], float]
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideHit_AnalysisResult()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideHit_AnalysisResult ) -> None:
        """
        Cython signature: void PeptideHit_AnalysisResult(PeptideHit_AnalysisResult &)
        """
        ... 


class PeptideHit_PeakAnnotation:
    """
    Cython implementation of _PeptideHit_PeakAnnotation

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1PeptideHit_PeakAnnotation.html>`_
    """
    
    annotation: Union[bytes, str, String]
    
    charge: int
    
    mz: float
    
    intensity: float
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void PeptideHit_PeakAnnotation()
        """
        ...
    
    @overload
    def __init__(self, in_0: PeptideHit_PeakAnnotation ) -> None:
        """
        Cython signature: void PeptideHit_PeakAnnotation(PeptideHit_PeakAnnotation &)
        """
        ...
    
    def writePeakAnnotationsString_(self, annotation_string: String , annotations: List[PeptideHit_PeakAnnotation] ) -> None:
        """
        Cython signature: void writePeakAnnotationsString_(String & annotation_string, libcpp_vector[PeptideHit_PeakAnnotation] annotations)
        """
        ...
    
    def __richcmp__(self, other: PeptideHit_PeakAnnotation, op: int) -> Any:
        ... 


class RANSAC:
    """
    Cython implementation of _RANSAC[_RansacModelLinear]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSAC[_RansacModelLinear].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RANSAC()
        """
        ...
    
    @overload
    def __init__(self, seed: int ) -> None:
        """
        Cython signature: void RANSAC(uint64_t seed)
        """
        ...
    
    def ransac(self, pairs: List[List[float, float]] , n: int , k: int , t: float , d: int , relative_d: bool ) -> List[List[float, float]]:
        """
        Cython signature: libcpp_vector[libcpp_pair[double,double]] ransac(libcpp_vector[libcpp_pair[double,double]] pairs, size_t n, size_t k, double t, size_t d, bool relative_d)
        """
        ... 


class RANSACParam:
    """
    Cython implementation of _RANSACParam

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSACParam.html>`_
    """
    
    n: int
    
    k: int
    
    t: float
    
    d: int
    
    relative_d: bool
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RANSACParam()
        A simple struct to carry all the parameters required for a RANSAC run
        """
        ...
    
    @overload
    def __init__(self, p_n: int , p_k: int , p_t: float , p_d: int , p_relative_d: bool ) -> None:
        """
        Cython signature: void RANSACParam(size_t p_n, size_t p_k, double p_t, size_t p_d, bool p_relative_d)
        """
        ...
    
    def toString(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ...
    
    def __str__(self) -> Union[bytes, str, String]:
        """
        Cython signature: String toString()
        """
        ... 


class RANSACQuadratic:
    """
    Cython implementation of _RANSAC[_RansacModelQuadratic]

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS::Math_1_1RANSAC[_RansacModelQuadratic].html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void RANSACQuadratic()
        """
        ...
    
    @overload
    def __init__(self, seed: int ) -> None:
        """
        Cython signature: void RANSACQuadratic(uint64_t seed)
        """
        ...
    
    def ransac(self, pairs: List[List[float, float]] , n: int , k: int , t: float , d: int , relative_d: bool ) -> List[List[float, float]]:
        """
        Cython signature: libcpp_vector[libcpp_pair[double,double]] ransac(libcpp_vector[libcpp_pair[double,double]] pairs, size_t n, size_t k, double t, size_t d, bool relative_d)
        """
        ... 


class Software:
    """
    Cython implementation of _Software

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1Software.html>`_
    """
    
    @overload
    def __init__(self, ) -> None:
        """
        Cython signature: void Software()
        """
        ...
    
    @overload
    def __init__(self, in_0: Software ) -> None:
        """
        Cython signature: void Software(Software &)
        """
        ...
    
    def getName(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getName()
        Returns the name of the software
        """
        ...
    
    def getVersion(self) -> Union[bytes, str, String]:
        """
        Cython signature: String getVersion()
        Returns the software version
        """
        ...
    
    def setName(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setName(String)
        Sets the name of the software
        """
        ...
    
    def setVersion(self, in_0: Union[bytes, str, String] ) -> None:
        """
        Cython signature: void setVersion(String)
        Sets the software version
        """
        ... 


class TransformationXMLFile:
    """
    Cython implementation of _TransformationXMLFile

    Original C++ documentation is available `here <http://www.openms.de/current_doxygen/html/classOpenMS_1_1TransformationXMLFile.html>`_
    """
    
    def __init__(self) -> None:
        """
        Cython signature: void TransformationXMLFile()
        """
        ...
    
    def load(self, in_0: Union[bytes, str, String] , in_1: TransformationDescription , fit_model: bool ) -> None:
        """
        Cython signature: void load(String, TransformationDescription &, bool fit_model)
        """
        ...
    
    def store(self, in_0: Union[bytes, str, String] , in_1: TransformationDescription ) -> None:
        """
        Cython signature: void store(String, TransformationDescription)
        """
        ... 


class BoundaryCondition:
    None
    BC_ZERO_ENDPOINTS : int
    BC_ZERO_FIRST : int
    BC_ZERO_SECOND : int

    def getMapping(self) -> Dict[int, str]:
       ... 


class MZTrafoModel_MODELTYPE:
    None
    LINEAR : int
    LINEAR_WEIGHTED : int
    QUADRATIC : int
    QUADRATIC_WEIGHTED : int
    SIZE_OF_MODELTYPE : int

    def getMapping(self) -> Dict[int, str]:
       ... 

