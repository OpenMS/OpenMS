"""
Nanobind C++ Code Emitter for pyOpenMS2 (v2 - using libclang info)

This module generates nanobind C++ binding code from merged C++/pxd information.
It uses accurate C++ type information from libclang for correct method signatures.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from .cpp_parser import CppMethod, CppParameter, MergedClass, MergedMethod


def scan_caster_files_for_types(caster_dir: Optional[Path] = None) -> Set[str]:
    """Scan type caster header files to find types that have casters.

    Types with casters should NOT be bound as classes since nanobind
    will automatically convert them via the caster.

    Returns a set of unqualified type names (e.g., 'String', 'DataValue').
    """
    if caster_dir is None:
        # Default to the standard type_casters directory
        caster_dir = Path(__file__).parent.parent / "bindings" / "type_casters"

    casted_types = set()

    if not caster_dir.exists():
        return casted_types

    # Pattern to match: struct type_caster<OpenMS::TypeName>
    # or NB_TYPE_CASTER(OpenMS::TypeName, ...)
    type_caster_pattern = re.compile(r'type_caster<OpenMS::(\w+)>')
    nb_type_caster_pattern = re.compile(r'NB_TYPE_CASTER\(OpenMS::(\w+)')

    for header_file in caster_dir.glob("*.h"):
        try:
            content = header_file.read_text()
            # Find all type_caster<OpenMS::X> declarations
            for match in type_caster_pattern.finditer(content):
                casted_types.add(match.group(1))
            # Find all NB_TYPE_CASTER(OpenMS::X, ...) declarations
            for match in nb_type_caster_pattern.finditer(content):
                casted_types.add(match.group(1))
        except Exception:
            pass  # Skip files that can't be read

    return casted_types


# Auto-detect types with casters - these should be skipped from class binding
_CASTER_OWNED_TYPES: Optional[Set[str]] = None

def get_caster_owned_types() -> Set[str]:
    """Get types that have type casters (cached)."""
    global _CASTER_OWNED_TYPES
    if _CASTER_OWNED_TYPES is None:
        _CASTER_OWNED_TYPES = scan_caster_files_for_types()
    return _CASTER_OWNED_TYPES


# Methods to skip for specific classes (complex return types or signature mismatches)
SKIP_METHODS = {
    "MSSpectrum": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "getIMData",  # Returns pair
        "calculateTIC",  # Return type mismatch (double in .pxd vs float in C++)
        "findNearest",  # Return type mismatch (Size vs int in .pxd)
        "findHighestInWindow",  # Return type mismatch
        # getName - const String& works fine with type caster
        "getType",  # Complex return type - handled via SPECIAL_METHODS
        "select",  # Complex parameter types
        "getSourceFile",  # const/non-const overload
    },
    "MSChromatogram": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "getSourceFile",  # const/non-const overload
    },
    "MSExperiment": {
        "getSpectrum", "getSpectra",  # const/non-const overloads - handled via SPECIAL_METHODS
        "getChromatogram", "getChromatograms",  # const/non-const overloads
        "getSourceFiles",  # const/non-const overloads
        "setSourceFile",  # const/non-const overloads
    },
    "Feature": {
        "getSubordinates",  # const/non-const overloads
        "getPeptideIdentifications",  # const/non-const overloads
        "getConvexHulls",  # const/non-const overloads
    },
    "FeatureMap": {
        "getProteinIdentifications",  # const/non-const overloads
        "getUnassignedPeptideIdentifications",  # const/non-const overloads
        "getDataProcessing",  # const/non-const overloads
    },
    "ConsensusFeature": {
        "getPeptideIdentifications",  # const/non-const overloads
        "computeDechargeConsensus",  # Takes forward-declared FeatureMap
    },
    "ConsensusMap": {
        "getProteinIdentifications",  # const/non-const overloads
        "getUnassignedPeptideIdentifications",  # const/non-const overloads
        "getDataProcessing",  # const/non-const overloads
        "getColumnHeaders",  # const/non-const overloads
    },
    "PeptideIdentification": {
        "getHits",  # const/non-const overloads - handled via SPECIAL_METHODS
    },
    "ProteinIdentification": {
        "getHits",  # const/non-const overloads - handled via SPECIAL_METHODS
        "getProteinGroups",  # const/non-const overloads
        "getIndistinguishableProteins",  # const/non-const overloads
        "computeCoverage",  # Takes ConsensusMap/PeptideIdentificationList - complex types
        "setPrimaryMSRunPath",  # Overloaded with MSExperiment parameter
        "getPrimaryMSRunPath",  # Output parameter
        "setSearchParameters",  # SearchParameters nested type
        "getSearchParameters",  # SearchParameters nested type
        "insertProteinGroup",  # ProteinGroup nested type
        "insertIndistinguishableProteins",  # ProteinGroup nested type
    },
    "AASequence": {
        "getFormula",  # Complex return type
        "begin", "end",  # Iterator types
        "getResidue",  # const/non-const overloads
    },
    "Param": {
        "begin", "end",  # Iterator types
        "getDescription", "getTags",  # Complex return types
    },
    "MzMLFile": {
        "getOptions",  # const/non-const overloads
    },
    "InternalCalibration": {
        "fillCalibrants",  # 7 parameters + self = 8, hits nanobind limit
    },
    "ModificationsDB": {
        "getInstance",  # Singleton pattern - needs wrap-manual-memory
        "getAllSearchModifications",  # Complex return type
        "getModification",  # Overloaded
        "addModification",  # Takes unique_ptr
    },
    "ResidueModification": {
        "getDiffFormula",  # Complex return type
    },
    "EmpiricalFormula": {
        "begin", "end",  # Iterator types
        "getIsotopeDistribution",  # Takes IsotopePatternGenerator
        "getConditionalFragmentIsotopeDist",  # Takes CoarseIsotopePatternGenerator
    },
    "Residue": {
        "getModification",  # Returns pointer
    },
    "Sample": {
        "getSubsamples",  # const/non-const overloads
    },
    "Precursor": {
        "getActivationMethods",  # const/non-const overloads
    },
    "PeptideIdentificationList": {
        # Inherits from ExposedVector - most methods come from there
    },
    "DataProcessing": {
        "getProcessingActions",  # Returns set
    },
    "InstrumentSettings": {
        # Most methods are simple getters/setters
    },
    "Mobilogram": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "findNearest",  # Return type mismatch
        "select",  # Complex parameter types
    },
    "BilinearInterpolation": {
        # Template class - may need special handling
    },
    "MSNumpressCoder": {
        # Static methods with complex buffer types
        "encodeNP", "decodeNP",
    },
    "MultipleTesting": {
        # Static methods - should work
    },
    "File": {
        # Static utility methods
        "getUniqueName",  # Returns temp file path
    },
    "PepXMLFile": {
        "load",  # Complex signature with PeptideIdentificationList
        "store",
    },
    "MzIdentMLFile": {
        "load",  # Complex signature with PeptideIdentificationList
        "store",
    },
    "DRange1": {
        # Template instantiation
    },
    "DRange2": {
        # Template instantiation
    },
    "DataFilters": {
        # passes() overloads with forward-declared types - skip the ones with Feature/ConsensusFeature
        "passes",  # Multiple overloads including Feature and ConsensusFeature
    },
    "ProteinInference": {
        # Methods with forward-declared ConsensusMap
        "run",  # Takes ConsensusMap
    },
}

# Classes to skip due to incomplete type dependencies or other issues
# These will need manual handling or fixes in the C++ headers
#
# NOTE: Types with type casters (String, DataValue, ParamValue, DPosition) are
# now AUTO-DETECTED by scanning bindings/type_casters/*.h and don't need to be
# listed here. Use get_caster_owned_types() to get the auto-detected list.
SKIP_CLASSES = {
    # Incomplete type issues
    "MassExplainer",        # References Compomer which is forward-declared
    "Compomer",             # Forward-declared type, not complete when included
    "ILPDCWrapper",         # Complex ILP dependencies
    "SwathFileConsumer",    # Template complexity
    "MSDataWritingConsumer", # Template complexity
    "SpectrumAccessOpenMSInMemory",  # Complex template
    "IsobaricChannelExtractor",  # Uses forward-declared IsobaricQuantitationMethod
    "IsobaricQuantifier",        # Uses forward-declared IsobaricQuantitationMethod
    "IsobaricNormalizer",        # Uses forward-declared IsobaricQuantitationMethod
    "IsobaricIsotopeCorrector",  # Uses forward-declared IsobaricQuantitationMethod
    # Abstract classes that libclang may not detect
    "BaseGroupFinder",
    "BaseSuperimposer",
    "ConsensusIDAlgorithm",
    "ConsensusIDAlgorithmIdentity",
    "ConsensusIDAlgorithmSimilarity",  # Abstract class
    "IsobaricQuantitationMethod",  # Pure virtual methods
    # Nested classes or special cases
    "XMLHandler",
    "IndexedMzMLHandler",
    # Classes with static methods not marked in .pxd (need @staticmethod or wrap-static)
    "PercolatorFeatureSetHelper",
    "PercolatorInfile",
    "PercolatorOutfile",
    "TransformationXMLFile",
    "OpenSwathDataAccessHelper",
    # Classes with no default constructor
    "QTCluster",
    "SpectrumAccessOpenMS",
    "ProFormaParser",               # Only explicit(string_view) constructor
    # Classes with Cython template syntax issues
    "TraMLFile",
    "MZTrafoModel",
    # Private inheritance issues (inherit from std::vector but not accessible)
    "AcquisitionInfo",

    # Classes with complex constructors that reference unbound types
    "TransformationModel",          # Base class with DataPoints nested type
    "TransformationModelBSpline",
    "TransformationModelLinear",
    "TransformationModelLowess",
    "TransformationModelInterpolated",

    # Classes with copy constructors that cause overload ambiguity
    "CVMappings",
    "ConsensusMapNormalizerAlgorithmMedian",
    "ConsensusMapNormalizerAlgorithmQuantile",

    # Classes where the type is actually a template alias (DPosition2 -> DPosition<2>)
    "DPosition2",

    # Classes with overloaded methods that can't be resolved without libclang
    "Biosaur2Algorithm",  # setMSData has const ref and move versions
    # ConsensusMap - handled via SPECIAL_METHODS
    # Feature - handled via SPECIAL_METHODS
    # FeatureMap - handled via SPECIAL_METHODS
    "FeatureFinderAlgorithmMetaboIdent",  # setMSData overloads
    "FeatureFinderIdentificationAlgorithm",  # constructor parameter issues

    # Classes with overloaded static methods (need explicit overload_cast)
    "IMTypes",                      # determineIMFormat has overloads
    "XQuestScores",                 # preScore has overloads
    "ExperimentalDesignFile",       # load has overloads
    "IDConflictResolverAlgorithm",  # resolve has overloads
    "OPXLHelper",                   # addProteinPositionMetaValues has overloads

    # Classes using abstract parser types
    "IMSAlphabet",                  # Uses abstract IMSAlphabetParser
    "IMSAlphabetTextParser",        # Abstract template

    # Classes with wrong parameter types in .pxd (int instead of complex types)
    "AccurateMassSearchResult",     # setIndividualIntensities needs vector<double>
    "ControlledVocabulary",         # getAllChildTerms needs set<String>&
    "TransitionPQPFile",            # convertPQPToTargetedExperiment param issues
    "TransitionTSVFile",            # convertTSVToTargetedExperiment param issues
    "MRMFeaturePickerFile",         # load() needs nested type vectors
    "MRMTransitionGroupPicker",     # findLargestPeak needs vector<MSChromatogram>
    "PeakIntegrator",               # iterator reference issues in integratePeak

    # Classes with OpenSwath type dependencies
    "MRMFeatureFinderScoring",      # prepareProteinPeptideMaps_ needs LightTargetedExperiment
    "OpenSwathOSWWriter",           # prepareLine needs LightCompound

    # Classes with namespace parsing issues
    "TargetedExperiment",           # setInstruments namespace issue
    "TargetedExperimentHelper",     # nested type references

    # Classes with incomplete type errors in lambdas
    "IndexedMzMLFileLoader",        # OnDiscMSExperiment incomplete type
    "OnDiscMSExperiment",           # incomplete type
    "ModificationDefinitionsSet",   # parameter type issues

    # Abstract classes with pure virtual destructor
    "SpectrumAccessTransforming",   # pure virtual destructor
    "SpectrumAccessQuadMZTransforming",
    "SpectrumAccessSqMass",
    "SpectrumAccessOpenMSCached",

    # Classes with QDate/Qt dependencies
    "Date",                         # constructor needs QDate, not int

    # Classes with unknown parameter types
    "FalseDiscoveryRate",           # ScoreToTgtDecLabelPairs unknown
    "ProteaseDB",                   # parameter type issues (int instead of proper type)
    "RNaseDB",                      # similar to ProteaseDB
    "SpectralDeconvolution",        # parameter type issues

    # Classes inheriting from XMLHandler (non-copyable)
    "QcMLFile",                     # XMLHandler base is non-copyable

    # Additional problematic classes found in build
    "Instrument",         # getIonSources has const/non-const overloads
    "ExperimentalDesign", # getters with overloads
    "ExperimentalSettings",  # similar pattern
    "CVMappingRule",      # copy constructor issues

    # Classes with deleted default constructor
    "MRMBatchFeatureSelector",      # default constructor deleted
    "AAIndex",                      # default constructor deleted

    # Classes with constructor type mismatches
    "NASequence",                   # constructor needs vector<const Ribonucleotide*>
    "SemanticValidator",            # ControlledVocabulary incomplete type

    # Classes with wrong parameter types (int instead of proper types)
    "AbsoluteQuantitationMethodFile",  # load/store need vector<AQMethod>
    "AbsoluteQuantitationStandardsFile",  # load needs vector<AQStandards>
    "DeconvolvedSpectrum",          # setCharges needs vector<int>
    "FLASHDeconvAlgorithm",         # performSpectrumDeconvolution params
    "FASTAFile",                    # load/store need vector<FASTAEntry>
    "MRMScoring",                   # OpenSwath namespace, param type issues

    # Classes with template parameters or complex constructors
    "OpenSwathWorkflowSonar",
    "OpenSwathWorkflow",
    "MRMTransitionGroup",
    "MRMTransitionGroupCP",
    "LightMRMTransitionGroup",

    # Classes with private copy constructors
    "CVMappingFile",

    # Classes with private/protected constructors (singletons, utility classes)
    "UniqueIdGenerator",            # Protected default constructor, private copy
    "GridBasedCluster",             # Complex constructors with nested Point/Rectangle types
    "LogConfigHandler",             # Private default constructor
    "MultipleTesting",              # Type not found/properly bound
    "PosteriorErrorProbabilityModel",  # OpenMS::Math:: namespace causes lambda analysis issues
    "MultiplexDeltaMassesGenerator",  # Constructor type mismatch: pxd says (String,int,int), C++ expects (String,int,map)
    "PeptideAndProteinQuant",  # PeptideIdentificationList/ExperimentalDesign types incomplete
    "FeatureGroupingAlgorithmLabeled",  # ConsensusMap type incomplete in lambda
    "FeatureGroupingAlgorithmUnlabeled",  # ConsensusMap type incomplete in lambda
    "KDTreeFeatureMaps",  # Lambda analysis fails
    "MapAlignmentAlgorithmKD",  # pxd type mismatch: int instead of vector<FeatureMap>
    "MapAlignmentAlgorithmPoseClustering",  # pxd type mismatch: int instead of vector<TransformationModelLowess*>
    "ReactionMonitoringTransition",  # Lambda analysis fails
    "MRMFeature",  # Lambda analysis fails (inherits from Feature)
    "IncludeExcludeTarget",  # Lambda analysis fails
    "InstrumentSettings",  # Lambda analysis fails
    "ScanWindow",  # pxd type mismatch
    "SignalToNoiseEstimatorMedianRapid",  # pxd type mismatch and lambda analysis fails
    "LPWrapper",  # pxd type mismatch: int instead of proper vector/reference types
    "ProteaseDigestion",  # pxd type mismatch in digest method

    # Classes with unresolved overloads
    "KroenikFile",
    "MascotGenericFile",
    "MzTabFile",
    "NLargest",
    "RNaseDigestion",
    "RankScaler",
    "PeakGroup",
    "SpectrumLookup",               # extractScanNumber has overloads
    "DIAScoring",                   # parameter type issues
    "MapAlignmentTransformer",      # transformRetentionTimes has overloads
    "SpectrumMetaDataLookup",       # getSpectrumMetaData has overloads
    # Mobilogram removed from SKIP_CLASSES - now has SPECIAL_METHODS

    # Template classes (need specialized instantiation)
    "Matrix",                       # Matrix<double> template
    "MassDecomposer",               # template class
    "DistanceMatrix",               # template class
    "RANSAC",                       # template class

    # Forward-declared/incomplete types
    "MzMLSqliteHandler",            # forward declared in MSDataSqlConsumer.h
    "OSWFile",                      # SqliteConnector header parsing issues

    # More classes with overloaded static methods
    "IDFilter",                     # countHits has overloads
    "NonNegativeLeastSquaresSolver", # solve has overloads

    # Classes with non-const reference constructor parameters
    "MSPGenericFile",               # constructor takes MSExperiment&

    # Classes with type alias issues
    "OSW_ChromExtractParams",

    # TODO: Add wrap-static support in .pxd files and remove from SKIP_CLASSES

    # Classes with constructor/copy issues found in module 3
    "MultiplexDeltaMasses",         # No matching constructor for initializer list
    "PeakWidthEstimator",           # No matching constructor for initializer list
    "RibonucleotideDB",             # Deleted copy constructor (singleton-like)

    # Classes with constructor parameter issues (const correctness)
    "ModifiedPeptideGenerator",     # MapToResidueType& should be const ref
    "IonMobilityScoring",           # OpenSwath::LightTransition& should be const ref
    "AASequence",                   # char* constructor, nanobind passes const char*
    "Peak2D",                       # DPosition<2> constructor issues
    "ChromatogramPeak",             # Similar to Peak2D

    # Classes referencing Param::ParamEntry (nested type incomplete)
    "ElutionModelFitter",           # Param::ParamEntry incomplete type
    "MSSim",                        # Param::ParamEntry incomplete type
    "File",                         # isDirectory param issues + Param::ParamEntry

    # Classes with forward-declared types or nested enum issues
    "OpenSwathWorkflowBase",        # Nested enum / forward decl issues
    "KDTreeFeatureNode",            # KDTreeFeatureMaps is forward-declared

    # Classes with lambda analysis failures (nanobind can't deduce types)
    "TransformationDescription",    # Lambda analysis fails
    "MorpheusScore",                # Copy assignment operator issues
    "Normalizer",                   # Lambda analysis fails
    "SeedListGenerator",            # Lambda analysis fails
    "NucleicAcidSearchEngine",      # Lambda analysis fails
    "MSFraggerAdapter",             # Lambda analysis fails
    "CometAdapter",                 # Lambda analysis fails
    "OMSSAAdapter",                 # Lambda analysis fails
    "MyriMatchAdapter",             # Lambda analysis fails
    "PepNovoAdapter",               # Lambda analysis fails

    # Classes with pxd type mismatches (method signature issues)
    "DecoyGenerator",               # shuffle() pxd mismatch
    "GaussFitter",                  # fit() pxd type mismatch (int vs vector<DPosition<2>>)
    "IsotopeDistribution",          # set() pxd type mismatch (int vs ContainerType)
    "MobilityPeak1D",               # Constructor type mismatch (int vs DPosition)
    "FileHandler",                  # Methods use int instead of ProgressLogger::LogType enum
    "SequestInfile",                # pxd type mismatch
    "SequestOutfile",               # pxd type mismatch
    "BSpline2d",                    # Lambda analysis fails
    "CrossLinksDB",                 # Lambda analysis fails
    "CubicSpline2d",                # Lambda analysis fails
    "FileTypes",                    # Lambda analysis fails
    "IMSIsotopeDistribution",       # Constructor mismatch in ims namespace
    "LinearInterpolation",          # Template class with unresolved KeyType
    "SignalToNoiseEstimator",       # Template class needs instantiation
    "SignalToNoiseEstimatorMeanIterative",  # Template class needs instantiation
    "SignalToNoiseEstimatorMedianRapid",  # Template class needs instantiation
    "FullSwathFileConsumer",        # Abstract class (can't instantiate)
    "SpectrumAlignmentScore",       # Template class
    "SteinScottImproveScore",       # Template class
    "SignalToNoiseEstimatorMedian", # Template class
    "IMSWeights",                   # Not found in C++ headers (nested/helper class)
    "IntegerMassDecomposer",        # Uses IMSWeights which isn't bound
    "RealMassDecomposer",           # Uses IMSWeights which isn't bound
    "BilinearInterpolation",        # Template class
}

# Additional headers needed for specific classes
ADDITIONAL_INCLUDES = {
    "MSSpectrum": ["<OpenMS/KERNEL/Peak1D.h>", "<OpenMS/METADATA/DataArrays.h>", "<OpenMS/IONMOBILITY/IMTypes.h>", "<OpenMS/CONCEPT/ProgressLogger.h>", "<OpenMS/FORMAT/FileTypes.h>"],
    "MSChromatogram": ["<OpenMS/KERNEL/ChromatogramPeak.h>"],
    "MSExperiment": ["<OpenMS/KERNEL/MSSpectrum.h>", "<OpenMS/KERNEL/MSChromatogram.h>", "<OpenMS/KERNEL/MSExperiment.h>"],
    "Feature": ["<OpenMS/KERNEL/Feature.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "FeatureMap": ["<OpenMS/KERNEL/FeatureMap.h>", "<OpenMS/KERNEL/Feature.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "ConsensusFeature": ["<OpenMS/KERNEL/ConsensusFeature.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "ConsensusMap": ["<OpenMS/KERNEL/ConsensusMap.h>", "<OpenMS/KERNEL/ConsensusFeature.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "PeptideIdentification": ["<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideHit.h>"],
    "ProteinIdentification": ["<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/ProteinHit.h>"],
    "PeptideHit": ["<OpenMS/METADATA/PeptideHit.h>"],
    "ProteinHit": ["<OpenMS/METADATA/ProteinHit.h>"],
    "PeptideEvidence": ["<OpenMS/METADATA/PeptideEvidence.h>"],
    "AASequence": ["<OpenMS/CHEMISTRY/AASequence.h>"],
    "EmpiricalFormula": ["<OpenMS/CHEMISTRY/EmpiricalFormula.h>"],
    "Residue": ["<OpenMS/CHEMISTRY/Residue.h>"],
    "Param": ["<OpenMS/DATASTRUCTURES/Param.h>"],
    "MzMLFile": ["<OpenMS/FORMAT/MzMLFile.h>", "<OpenMS/KERNEL/MSExperiment.h>"],
    "IdXMLFile": ["<OpenMS/FORMAT/IdXMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "FeatureXMLFile": ["<OpenMS/FORMAT/FeatureXMLFile.h>", "<OpenMS/KERNEL/FeatureMap.h>"],
    "ConsensusXMLFile": ["<OpenMS/FORMAT/ConsensusXMLFile.h>", "<OpenMS/KERNEL/ConsensusMap.h>"],
    "SourceFile": ["<OpenMS/METADATA/SourceFile.h>"],
    "Software": ["<OpenMS/METADATA/Software.h>"],
    "ContactPerson": ["<OpenMS/METADATA/ContactPerson.h>"],
    "Sample": ["<OpenMS/METADATA/Sample.h>"],
    "Precursor": ["<OpenMS/METADATA/Precursor.h>"],
    "Product": ["<OpenMS/METADATA/Product.h>"],
    "CVTerm": ["<OpenMS/METADATA/CVTerm.h>"],
    "ModificationsDB": ["<OpenMS/CHEMISTRY/ModificationsDB.h>"],
    "ResidueModification": ["<OpenMS/CHEMISTRY/ResidueModification.h>"],
    "PeptideIdentificationList": ["<OpenMS/METADATA/PeptideIdentificationList.h>", "<OpenMS/METADATA/PeptideIdentification.h>"],
    "DataProcessing": ["<OpenMS/METADATA/DataProcessing.h>"],
    "InstrumentSettings": ["<OpenMS/METADATA/InstrumentSettings.h>"],
    "Mobilogram": ["<OpenMS/KERNEL/Mobilogram.h>", "<OpenMS/KERNEL/MobilityPeak1D.h>", "<OpenMS/METADATA/DataArrays.h>"],
    "BilinearInterpolation": ["<OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>"],
    "MSNumpressCoder": ["<OpenMS/FORMAT/MSNumpressCoder.h>"],
    "MultipleTesting": ["<OpenMS/MATH/STATISTICS/MultipleTesting.h>"],
    "File": ["<OpenMS/SYSTEM/File.h>"],
    "PepXMLFile": ["<OpenMS/FORMAT/PepXMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "MzIdentMLFile": ["<OpenMS/FORMAT/MzIdentMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "DRange1": ["<OpenMS/DATASTRUCTURES/DRange.h>"],
    "DRange2": ["<OpenMS/DATASTRUCTURES/DRange.h>"],
    "DataFilters": ["<OpenMS/PROCESSING/MISC/DataFilters.h>", "<OpenMS/KERNEL/Feature.h>", "<OpenMS/KERNEL/ConsensusFeature.h>"],
    "ProteinInference": ["<OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>", "<OpenMS/KERNEL/ConsensusMap.h>"],
}

# FALLBACK: Classes that need special __len__ support (container-like)
# NOTE: When libclang is used, this is detected automatically via size() method.
# TODO: Remove once all modes use AST-based detection
CONTAINER_CLASSES = {
    "MSSpectrum", "MSChromatogram", "MSExperiment", "FeatureMap",
    "ConsensusMap", "PeakMap", "Mobilogram",
    "FloatDataArray", "IntegerDataArray", "StringDataArray",
}

# FALLBACK: Classes that need iterator support
# NOTE: When libclang is used, this is detected automatically via begin()/end() methods.
# TODO: Remove once all modes use AST-based detection
ITERABLE_CLASSES = {
    "MSSpectrum", "MSChromatogram", "MSExperiment", "FeatureMap",
    "ConsensusMap", "PeakMap", "Mobilogram", "AASequence",
}

# Core classes to bind first - these are well-tested and have simple APIs
# All other classes require more work on type normalization
CORE_CLASSES = {
    # === BATCH 1: Already working (50 classes) ===
    # Basic peak types
    "Peak1D", "Peak2D", "ChromatogramPeak", "MobilityPeak1D",
    # Spectrum and chromatogram
    "MSSpectrum", "MSChromatogram",
    # Mobilogram (ion mobility)
    "Mobilogram",
    # Data arrays (required for MSSpectrum tests)
    "FloatDataArray", "IntegerDataArray", "StringDataArray",
    # Basic data types (DataValue has type caster, excluded)
    "DateTime",
    # MSExperiment / PeakMap
    "MSExperiment",
    # Identification classes
    "PeptideIdentification", "ProteinIdentification",
    "PeptideHit", "ProteinHit", "PeptideEvidence",
    "PeptideIdentificationList",  # Container for peptide identifications
    # Sequence classes
    "AASequence", "EmpiricalFormula", "Residue",
    # Param handling
    "Param",
    # Feature classes (handled via special methods to avoid overload issues)
    "Feature", "FeatureMap",
    "ConsensusFeature", "ConsensusMap",
    # File I/O
    "MzMLFile", "IdXMLFile", "FeatureXMLFile", "ConsensusXMLFile",
    "PepXMLFile", "MzIdentMLFile",
    # Metadata
    "SourceFile", "Software", "ContactPerson", "Sample",
    "Precursor", "Product",
    "DataProcessing", "InstrumentSettings",
    # CV terms
    "CVTerm",
    # ModificationsDB access
    "ModificationsDB", "ResidueModification",
    # Compression
    "MSNumpressCoder",
    # Statistics
    "MultipleTesting",
    # File utilities
    "File",

    # === BATCH 2: Classes A-B ===
    "AAIndex", "AASeqWithMass", "AScore",
    "AbsoluteQuantitation", "AbsoluteQuantitationMethod", "AbsoluteQuantitationMethodFile",
    "AbsoluteQuantitationStandards", "AbsoluteQuantitationStandardsFile",
    "AccurateMassSearchEngine", "AccurateMassSearchResult",
    "Acquisition", "AcquisitionInfo",
    "Adduct", "AdductInfo",
    "AnnotationStatistics",
    "Attachment",
    "AverageLinkage",
    "BSpline2d", "Base64", "BaseFeature",
    "BasicProteinInferenceAlgorithm", "BayesianProteinInferenceAlgorithm",
    "BiGaussFitter1D", "BiGaussModel",
    "BinnedSpectrum",
    "BuildInfo",

    # === BATCH 3: Classes C ===
    "CVMappingFile", "CVMappingRule", "CVMappingTerm", "CVMappings", "CVReference",
    "CVTermList", "CVTermListInterface",
    "CachedmzML", "CachedmzMLHandler",
    "CalibrationData",
    "ChargePair",
    "ChromExtractParams",
    "ChromatogramExtractor", "ChromatogramExtractorAlgorithm",
    "ChromatogramSettings", "ChromatogramTools",
    "ChromeleonFile",
    "ClusterProxyKD", "ClusteringGrid",
    "Compomer",
    "ConfidenceScoring",
    "ConsensusIDAlgorithm", "ConsensusIDAlgorithmAverage", "ConsensusIDAlgorithmBest",
    "ConsensusIDAlgorithmIdentity", "ConsensusIDAlgorithmPEPIons", "ConsensusIDAlgorithmPEPMatrix",
    "ConsensusIDAlgorithmRanks", "ConsensusIDAlgorithmSimilarity", "ConsensusIDAlgorithmWorst",
    "ConsensusMapNormalizerAlgorithmMedian", "ConsensusMapNormalizerAlgorithmQuantile",
    "ConsensusMapNormalizerAlgorithmThreshold",
    "ControlledVocabulary",
    "ConvexHull2D",
    "CrossLinkSpectrumMatch", "CrossLinksDB",
    "CsvFile",
    "CubicSpline2d",

    # === BATCH 4: Classes D ===
    "DIAScoring",
    "DTA2DFile", "DTAFile",
    "DataFilters",
    "Date",
    "DeconvolvedSpectrum",
    "DecoyGenerator",
    "DefaultParamHandler",
    "Deisotoper",
    "DigestionEnzyme", "DigestionEnzymeProtein", "DigestionEnzymeRNA",
    "DistanceMatrix",
    "DocumentIdentifier",

    # === BATCH 5: Classes E ===
    "EDTAFile",
    "Element", "ElementDB",
    "ElutionModelFitter",
    "ElutionPeakDetection",
    "EmgFitter1D", "EmgGradientDescent", "EmgModel", "EmgScoring",
    "EnzymaticDigestion",
    "ExperimentalDesign", "ExperimentalDesignFile",
    "ExperimentalSettings",

    # === BATCH 6: Classes F ===
    "FASTAFile",
    "FIAMSDataProcessor", "FIAMSScheduler",
    "FLASHDeconvAlgorithm", "FLASHDeconvFeatureFile", "FLASHDeconvSpectrumFile",
    "FalseDiscoveryRate",
    "FeatureDeconvolution",
    "FeatureDistance",
    "FeatureFileOptions",
    "FeatureFinderAlgorithmMetaboIdent", "FeatureFinderAlgorithmPicked",
    "FeatureFinderIdentificationAlgorithm", "FeatureFinderMultiplexAlgorithm",
    "FeatureFindingMetabo",
    "FeatureGroupingAlgorithm", "FeatureGroupingAlgorithmKD", "FeatureGroupingAlgorithmLabeled",
    "FeatureGroupingAlgorithmQT", "FeatureGroupingAlgorithmUnlabeled",
    "FeatureHandle",
    "FeatureMapping",
    "FeatureOverlapFilter",
    "FileHandler", "FileTypes",
    "Fitter1D",

    # === BATCH 7: Classes G-H ===
    "GNPSMGFFile", "GNPSMetaValueFile", "GNPSQuantificationFile",
    "GaussFilter", "GaussFitter", "GaussTraceFitter",
    "Gradient",
    "GridBasedCluster",
    "HPLC",
    "HyperScore",

    # === BATCH 8: Classes I ===
    "IBSpectraFile",
    "IDConflictResolverAlgorithm",
    "IDDecoyProbability",
    "IDFilter",
    "IDMapper",
    "IDRipper",
    "ILPDCWrapper",
    "IMSAlphabet", "IMSAlphabetTextParser",
    "IMSElement", "IMSIsotopeDistribution",
    "IMTypes",
    "IncludeExcludeTarget",
    "IndexedMzMLDecoder", "IndexedMzMLFile", "IndexedMzMLFileLoader",
    "InspectInfile", "InspectOutfile",
    "Instrument",
    "InternalCalibration",
    "InterpolationModel",
    "IonDetector",
    "IonIdentityMolecularNetworking",
    "IonSource",
    "IsobaricChannelExtractor", "IsobaricChannelInformation",
    "IsobaricIsotopeCorrector", "IsobaricNormalizer", "IsobaricQuantifier",
    "IsobaricQuantifierStatistics", "IsobaricQuantitationMethod",
    "IsotopeCluster", "IsotopeDistribution",
    "IsotopeFitter1D", "IsotopeLabelingMDVs", "IsotopeModel",
    "ItraqConstants", "ItraqEightPlexQuantitationMethod", "ItraqFourPlexQuantitationMethod",

    # === BATCH 9: Classes J-K ===
    "JavaInfo",
    "KDTreeFeatureMaps", "KDTreeFeatureNode",
    "KernelDensityEstimation",
    "KroenikFile",

    # === BATCH 10: Classes L ===
    "LPWrapper",
    "LabeledPairFinder",
    "LevMarqFitter1D",
    "LightTargetedExperiment",
    "LinearResampler", "LinearResamplerAlign",
    "LogConfigHandler",
    "LowessSmoothing",

    # === BATCH 11: Classes M (MRM) ===
    "MRMAssay",
    "MRMBatchFeatureSelector",
    "MRMDecoy",
    "MRMFeature", "MRMFeatureFilter", "MRMFeatureFinderScoring",
    "MRMFeaturePicker", "MRMFeaturePickerFile",
    "MRMFeatureQC", "MRMFeatureQCFile",
    "MRMFeatureSelector",
    "MRMIonSeries",
    "MRMMapping",
    "MRMRTNormalizer",
    "MRMScoring",
    "MRMTransitionGroup", "MRMTransitionGroupPicker",

    # === BATCH 12: Classes M (MS, Map, Mass) ===
    "MS2File",
    "MSDataAggregatingConsumer", "MSDataCachedConsumer",
    "MSDataSqlConsumer", "MSDataStoringConsumer", "MSDataWritingConsumer",
    "MSPFile", "MSPGenericFile",
    "MSstatsFile",
    "MZTrafoModel",
    "MapAlignmentAlgorithmIdentification", "MapAlignmentAlgorithmKD",
    "MapAlignmentAlgorithmPoseClustering",
    "MapAlignmentEvaluationAlgorithm", "MapAlignmentEvaluationAlgorithmPrecision",
    "MapAlignmentEvaluationAlgorithmRecall",
    "MapAlignmentTransformer",
    "MascotGenericFile", "MascotXMLFile",
    "MassAnalyzer",
    "MassDecomposer", "MassDecomposition", "MassDecompositionAlgorithm",
    "MassExplainer",
    "MassTrace", "MassTraceDetection",
    "MasstraceCorrelator",
    "Matrix",

    # === BATCH 13: Classes M (Meta, Mod, Mor, Mz) ===
    "MetaInfo", "MetaInfoDescription", "MetaInfoInterface", "MetaInfoRegistry",
    "MetaboTargetedAssay", "MetaboTargetedTargetDecoy",
    "MetaboliteFeatureDeconvolution", "MetaboliteSpectralMatching",
    "ModificationDefinition", "ModificationDefinitionsSet",
    "ModifiedPeptideGenerator",
    "MorpheusScore",
    "MorphologicalFilter",
    "MsInspectFile",
    "MultiplexDeltaMasses", "MultiplexDeltaMassesGenerator", "MultiplexIsotopicPeakPattern",
    "MzDataFile",
    "MzMLSpectrumDecoder", "MzMLSqliteHandler", "MzMLValidator",
    "MzQCFile",
    "MzTab", "MzTabFile", "MzTabM", "MzTabMFile",
    "MzXMLFile",

    # === BATCH 14: Classes N-O ===
    "NASequence",
    "NLargest",
    "NonNegativeLeastSquaresSolver",
    "NucleicAcidSpectrumGenerator",
    "OMSSACSVFile", "OMSSAXMLFile",
    "OPXLDataStructs", "OPXLHelper", "OPXLSpectrumProcessingAlgorithms",
    "OSChromatogramMeta", "OSSpectrumMeta",
    "OSWFile",
    "OnDiscMSExperiment",
    "OpenPepXLAlgorithm",
    "OpenSwathDataAccessHelper", "OpenSwathDataStructures",
    "OpenSwathHelper", "OpenSwathOSWWriter", "OpenSwathScoring",

    # === BATCH 15: Classes P ===
    "PScore",
    "ParamCTDFile", "ParamEntry", "ParamNode", "ParamValue", "ParamXMLFile",
    "PeakFileOptions",
    "PeakIndex",
    "PeakIntegrator",
    "PeakPickerChromatogram", "PeakPickerHiRes", "PeakPickerIM", "PeakPickerIterative",
    "PeakTypeEstimator", "PeakWidthEstimator",
    "PepXMLFileMascot",
    "PeptideAndProteinQuant",
    "PeptideIndexing",
    "PeptideProteinResolution",
    "PeptideSearchEngineFIAlgorithm",
    "PercolatorFeatureSetHelper", "PercolatorInfile", "PercolatorOutfile",
    # PosteriorErrorProbabilityModel removed - OpenMS::Math:: namespace causes issues
    "PrecursorCorrection", "PrecursorPurity",
    "PreprocessedPairSpectra",
    "ProFormaData", "ProFormaParser",
    "ProgressLogger",
    "ProtXMLFile",
    "ProteaseDB", "ProteaseDigestion",
    "ProteinInference",
    "ProteinProteinCrossLink",

    # === BATCH 16: Classes Q-R ===
    "QTCluster", "QTClusterFinder",
    "QcMLFile",
    "RANSAC", "RANSACModelLinear", "RANSACModelQuadratic",
    "RNaseDB", "RNaseDigestion",
    "RankScaler",
    "ReactionMonitoringTransition",
    "RealMassDecomposer",
    "ResidueDB",
    "Ribonucleotide", "RibonucleotideDB",
    "RichPeak2D",

    # === BATCH 17: Classes S (Sa-Sp) ===
    "SavitzkyGolayFilter",
    "ScanWindow",
    "SeedListGenerator",
    "SemanticValidator",
    "SequestInfile", "SequestOutfile",
    "SignalToNoiseEstimatorMedianRapid",
    "SimpleOpenMSSpectraFactory",
    "SimpleSearchEngineAlgorithm",
    "SimpleTSGXLMS",
    "SiriusExportAlgorithm", "SiriusFragmentAnnotation", "SiriusMSFile",
    "SpectraMerger",
    "SpectraSTSimilarityScore",
    "SpectralDeconvolution",
    "SpectrumAccessOpenMS", "SpectrumAccessOpenMSCached", "SpectrumAccessOpenMSInMemory",
    "SpectrumAccessQuadMZTransforming", "SpectrumAccessSqMass", "SpectrumAccessTransforming",
    "SpectrumAlignment", "SpectrumAlignmentScore",
    "SpectrumAnnotator",
    "SpectrumLookup", "SpectrumMetaDataLookup",
    "SpectrumSettings",
    "SplineInterpolatedPeaks", "SplinePackage",

    # === BATCH 18: Classes S (Sq-Sw) ===
    "SqMassFile",
    "SqrtScaler",
    "StablePairFinder",
    "SwathFile", "SwathFileConsumer",
    "SwathMap", "SwathMapMassCorrection", "SwathWindowLoader",

    # === BATCH 19: Classes T ===
    "TMTEighteenPlexQuantitationMethod", "TMTElevenPlexQuantitationMethod",
    "TMTSixPlexQuantitationMethod", "TMTSixteenPlexQuantitationMethod",
    "TMTTenPlexQuantitationMethod",
    "Tagger",
    "TargetedExperiment", "TargetedExperimentHelper",
    "TargetedSpectraExtractor",
    "TextFile",
    "TheoreticalSpectrumGenerator", "TheoreticalSpectrumGeneratorXLMS",
    "ThresholdMower",
    "TraMLFile",
    "TransformationDescription", "TransformationModel",
    "TransformationModelBSpline", "TransformationModelInterpolated",
    "TransformationModelLinear", "TransformationModelLowess",
    "TransformationXMLFile",
    "TransitionPQPFile", "TransitionTSVFile",

    # === BATCH 20: Classes U-X ===
    "UnimodXMLFile",
    "UniqueIdGenerator", "UniqueIdInterface",
    "VersionInfo",
    "WindowMower",
    "XFDRAlgorithm",
    "XLPrecursor", "XLPrecursorComparator",
    "XMLFile", "XMLHandler",
    "XQuestResultXMLFile", "XQuestScores",
    "XTandemInfile", "XTandemXMLFile",

    # Note: Excluded classes (need special handling):
    # - BilinearInterpolation: template class
    # - DPosition, DRange, DBoundingBox: have type casters
    # - DataValue, ParamValue: have type casters (ParamValue included for method access)
    # - DoubleList, IntList, StringList, StringView: basic types with casters
    # - LinearInterpolation: template class
    # - Normalizer: filterSpectrum is a template method
    # - String: has type caster
    # - Types, Math, ctime, smart_ptr, streampos: low-level/utility
    # - AnnotatedMSRun: MSExperiment typedef
    # - BaseGroupFinder: abstract class
    # - ChromatogramRangeManager, SpectrumRangeManager, RangeManager: templates
    # - ConversionHelper: static utility functions
    # - DataArrays: namespace-like
    # - FLASHHelperClasses: nested types
    # - IMSAlphabetParser, IMSDataConsumer, ISpectrumAccess: abstract/interfaces
    # - IntegerMassDecomposer: template
    # - InterfaceDataStructures: structs/helpers
    # - FeatureFinderAlgorithmPickedHelperStructs: helper structs
    # - PeakGroup: complex inheritance
    # - PythonMSDataConsumer: special Python handling
    # - RankData: simple struct
    # - SignalToNoiseEstimator, SignalToNoiseEstimatorMeanIterative, SignalToNoiseEstimatorMedian: templates
    # - SpectrumHelper: static helpers
    # - SpectrumNativeIDParser: static
    # - Weights: specialized struct
    # - Biosaur2Algorithm: if it has issues
}

# FALLBACK: Classes that inherit from std::vector
# NOTE: When libclang is used, this is detected automatically via canonical base classes.
# The element type is also extracted from std::vector<T> inheritance.
# TODO: Remove once all modes use AST-based detection
VECTOR_BASED_CLASSES = {
    "MSSpectrum", "MSChromatogram", "Mobilogram",
    "FloatDataArray", "IntegerDataArray", "StringDataArray",
}

# Methods to skip for vector-based classes (inherited from private std::vector base)
# We provide lambda implementations instead via CONTAINER_CLASSES/__len__
VECTOR_INHERITED_METHODS = {
    "size", "reserve", "resize", "push_back", "pop_back", "clear",
    "empty", "front", "back", "begin", "end", "cbegin", "cend",
    "at", "operator[]", "data", "capacity", "max_size", "shrink_to_fit",
    "insert", "erase", "emplace", "emplace_back", "assign", "swap",
}

# Special method implementations for critical classes
# These are C++ lambdas that implement Python methods
SPECIAL_METHODS = {
    "MSSpectrum": {
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::MSSpectrum& self) {
            // Return (mz_array, intensity_array) as numpy float64 arrays
            // Both returned as float64 for pyOpenMS backward compatibility
            const size_t n = self.size();
            double* mz_data = new double[n];
            double* int_data = new double[n];
            for (size_t i = 0; i < n; ++i) {
                mz_data[i] = self[i].getMZ();
                int_data[i] = static_cast<double>(self[i].getIntensity());
            }
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            auto mz_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(mz_data, {n}, mz_owner);
            auto int_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(int_data, {n}, int_owner);
            return nb::make_tuple(mz_arr, int_arr);
        }, "Returns a tuple of (mz_array, intensity_array) as numpy arrays")''',
        "set_peaks": '''
        .def("set_peaks", [](OpenMS::MSSpectrum& self, nb::object mz_obj, nb::object int_obj) {
            // Accept two arrays (mz, intensity) for compatibility with pyOpenMS API
            std::vector<double> mz = nb::cast<std::vector<double>>(mz_obj);
            std::vector<double> intensity = nb::cast<std::vector<double>>(int_obj);
            const size_t n = mz.size();
            if (intensity.size() != n) {
                throw std::runtime_error("mz and intensity arrays must have same length");
            }
            self.resize(n);
            for (size_t i = 0; i < n; ++i) {
                self[i].setMZ(mz[i]);
                self[i].setIntensity(intensity[i]);
            }
        }, "mz"_a, "intensity"_a, "Set peaks from mz and intensity arrays")''',
        "push_back": '''
        .def("push_back", [](OpenMS::MSSpectrum& self, const OpenMS::Peak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak to the spectrum")''',
        "size": '''
        .def("size", [](const OpenMS::MSSpectrum& self) {
            return self.size();
        }, "Returns the number of peaks")''',
        "clear": '''
        .def("clear", [](OpenMS::MSSpectrum& self, bool clear_meta_data) {
            self.clear(clear_meta_data);
        }, "clear_meta_data"_a = true, "Remove all peaks (and optionally metadata)")''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::MSSpectrum& self, size_t i) -> OpenMS::Peak1D& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "getMinMZ": '''
        .def("getMinMZ", [](const OpenMS::MSSpectrum& self) {
            return self.getMinMZ();
        }, "Returns the minimum m/z value")''',
        "getMaxMZ": '''
        .def("getMaxMZ", [](const OpenMS::MSSpectrum& self) {
            return self.getMaxMZ();
        }, "Returns the maximum m/z value")''',
        "getMinIntensity": '''
        .def("getMinIntensity", [](const OpenMS::MSSpectrum& self) {
            return self.getMinIntensity();
        }, "Returns the minimum intensity")''',
        "getMaxIntensity": '''
        .def("getMaxIntensity", [](const OpenMS::MSSpectrum& self) {
            return self.getMaxIntensity();
        }, "Returns the maximum intensity")''',
        "get_mz_array": '''
        .def("get_mz_array", [](const OpenMS::MSSpectrum& self) {
            size_t n = self.size();
            double* data = new double[n];
            for (size_t i = 0; i < n; ++i) data[i] = self[i].getMZ();
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "Returns m/z values as numpy array")''',
        "get_intensity_array": '''
        .def("get_intensity_array", [](const OpenMS::MSSpectrum& self) {
            size_t n = self.size();
            float* data = new float[n];
            for (size_t i = 0; i < n; ++i) data[i] = self[i].getIntensity();
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
        }, "Returns intensity values as numpy array")''',
        "get_drift_time_array": '''
        .def("get_drift_time_array", [](const OpenMS::MSSpectrum& self) -> std::optional<nb::ndarray<nb::numpy, float, nb::ndim<1>>> {
            // Check if IM data exists
            if (!self.containsIMData()) return std::nullopt;
            const auto& fda = self.getFloatDataArrays();
            for (const auto& arr : fda) {
                if (arr.getName() == "Ion Mobility" || arr.getMetaValue("name") == "Ion Mobility") {
                    size_t n = arr.size();
                    float* data = new float[n];
                    std::copy(arr.begin(), arr.end(), data);
                    nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
                    return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
                }
            }
            return std::nullopt;
        }, "Returns drift time array if ion mobility data exists, else None")''',
        "get_drift_time_array_mv": '''
        .def("get_drift_time_array_mv", [](OpenMS::MSSpectrum& self) -> std::optional<nb::ndarray<nb::numpy, float, nb::ndim<1>>> {
            // Memory view version - returns view into float data array (zero-copy)
            if (!self.containsIMData()) return std::nullopt;
            auto& fda = self.getFloatDataArrays();
            for (auto& arr : fda) {
                if (arr.getName() == "Ion Mobility" || arr.getMetaValue("name") == "Ion Mobility") {
                    if (arr.empty()) return std::nullopt;
                    return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                        arr.data(), {arr.size()}, nb::handle()
                    );
                }
            }
            return std::nullopt;
        }, "Returns view of drift time array if ion mobility data exists, else None")''',
        "getFloatDataArrays": '''
        .def("getFloatDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::FloatDataArray>& {
            return self.getFloatDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the float data arrays")''',
        "setFloatDataArrays": '''
        .def("setFloatDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::FloatDataArray>& arrays) {
            self.setFloatDataArrays(arrays);
        }, "arrays"_a, "Set the float data arrays")''',
        "getIntegerDataArrays": '''
        .def("getIntegerDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::IntegerDataArray>& {
            return self.getIntegerDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the integer data arrays")''',
        "setIntegerDataArrays": '''
        .def("setIntegerDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::IntegerDataArray>& arrays) {
            self.setIntegerDataArrays(arrays);
        }, "arrays"_a, "Set the integer data arrays")''',
        "getStringDataArrays": '''
        .def("getStringDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::StringDataArray>& {
            return self.getStringDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the string data arrays")''',
        "setStringDataArrays": '''
        .def("setStringDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::StringDataArray>& arrays) {
            self.setStringDataArrays(arrays);
        }, "arrays"_a, "Set the string data arrays")''',
        "getType": '''
        .def("getType", [](const OpenMS::MSSpectrum& self, bool query_data) {
            return self.getType(query_data);
        }, "query_data"_a = false, "Returns the spectrum type (centroided, profile, etc.)")''',
        "setType": '''
        .def("setType", [](OpenMS::MSSpectrum& self, OpenMS::SpectrumSettings::SpectrumType type) {
            self.setType(type);
        }, "type"_a, "Set the spectrum type")''',
        "get_drift_time_unit": '''
        .def("get_drift_time_unit", [](const OpenMS::MSSpectrum& self) -> std::optional<OpenMS::DriftTimeUnit> {
            if (!self.containsIMData()) return std::nullopt;
            return self.getDriftTimeUnit();
        }, "Returns drift time unit if ion mobility data exists, else None")''',
    },
    "MSChromatogram": {
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::MSChromatogram& self) {
            const size_t n = self.size();
            std::vector<double> rt(n);
            std::vector<double> intensity(n);
            for (size_t i = 0; i < n; ++i) {
                rt[i] = self[i].getRT();
                intensity[i] = self[i].getIntensity();
            }
            return nb::make_tuple(rt, intensity);
        }, "Returns a tuple of (rt_array, intensity_array)")''',
    },
    # Standalone enum bindings that need to be added before classes that use them
    "__enums__": {
        "DriftTimeUnit": '''
    // DriftTimeUnit enum
    nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit")
        .value("NONE", OpenMS::DriftTimeUnit::NONE)
        .value("MILLISECOND", OpenMS::DriftTimeUnit::MILLISECOND)
        .value("VSSC", OpenMS::DriftTimeUnit::VSSC)
        .export_values();''',
        "LogType": '''
    // ProgressLogger::LogType enum
    nb::enum_<OpenMS::ProgressLogger::LogType>(m, "LogType")
        .value("CMD", OpenMS::ProgressLogger::LogType::CMD)
        .value("GUI", OpenMS::ProgressLogger::LogType::GUI)
        .value("NONE", OpenMS::ProgressLogger::LogType::NONE)
        .export_values();''',
        "FileType": '''
    // FileTypes::Type enum
    nb::enum_<OpenMS::FileTypes::Type>(m, "FileType")
        .value("UNKNOWN", OpenMS::FileTypes::Type::UNKNOWN)
        .value("DTA", OpenMS::FileTypes::Type::DTA)
        .value("DTA2D", OpenMS::FileTypes::Type::DTA2D)
        .value("MZDATA", OpenMS::FileTypes::Type::MZDATA)
        .value("MZXML", OpenMS::FileTypes::Type::MZXML)
        .value("FEATUREXML", OpenMS::FileTypes::Type::FEATUREXML)
        .value("IDXML", OpenMS::FileTypes::Type::IDXML)
        .value("CONSENSUSXML", OpenMS::FileTypes::Type::CONSENSUSXML)
        .value("MGF", OpenMS::FileTypes::Type::MGF)
        .value("INI", OpenMS::FileTypes::Type::INI)
        .value("TOPPAS", OpenMS::FileTypes::Type::TOPPAS)
        .value("TRANSFORMATIONXML", OpenMS::FileTypes::Type::TRANSFORMATIONXML)
        .value("MZML", OpenMS::FileTypes::Type::MZML)
        .value("CACHEDMZML", OpenMS::FileTypes::Type::CACHEDMZML)
        .value("MS2", OpenMS::FileTypes::Type::MS2)
        .value("PEPXML", OpenMS::FileTypes::Type::PEPXML)
        .value("PROTXML", OpenMS::FileTypes::Type::PROTXML)
        .value("MZIDENTML", OpenMS::FileTypes::Type::MZIDENTML)
        .value("QCML", OpenMS::FileTypes::Type::QCML)
        .value("MZQC", OpenMS::FileTypes::Type::MZQC)
        .value("GELML", OpenMS::FileTypes::Type::GELML)
        .value("TRAML", OpenMS::FileTypes::Type::TRAML)
        .value("MSP", OpenMS::FileTypes::Type::MSP)
        .value("OMSSAXML", OpenMS::FileTypes::Type::OMSSAXML)
        .value("MASCOTXML", OpenMS::FileTypes::Type::MASCOTXML)
        .value("PNG", OpenMS::FileTypes::Type::PNG)
        .value("XMASS", OpenMS::FileTypes::Type::XMASS)
        .value("TSV", OpenMS::FileTypes::Type::TSV)
        .value("MZTAB", OpenMS::FileTypes::Type::MZTAB)
        .value("PEPLIST", OpenMS::FileTypes::Type::PEPLIST)
        .value("HARDKLOER", OpenMS::FileTypes::Type::HARDKLOER)
        .value("KROENIK", OpenMS::FileTypes::Type::KROENIK)
        .value("FASTA", OpenMS::FileTypes::Type::FASTA)
        .value("EDTA", OpenMS::FileTypes::Type::EDTA)
        .value("CSV", OpenMS::FileTypes::Type::CSV)
        .value("TXT", OpenMS::FileTypes::Type::TXT)
        .value("OBO", OpenMS::FileTypes::Type::OBO)
        .value("HTML", OpenMS::FileTypes::Type::HTML)
        .value("ANALYSISXML", OpenMS::FileTypes::Type::ANALYSISXML)
        .value("XSD", OpenMS::FileTypes::Type::XSD)
        .value("PSQ", OpenMS::FileTypes::Type::PSQ)
        .value("MRM", OpenMS::FileTypes::Type::MRM)
        .value("SQMASS", OpenMS::FileTypes::Type::SQMASS)
        .value("PQP", OpenMS::FileTypes::Type::PQP)
        .value("MS", OpenMS::FileTypes::Type::MS)
        .value("OSW", OpenMS::FileTypes::Type::OSW)
        .value("PSMS", OpenMS::FileTypes::Type::PSMS)
        .value("PIN", OpenMS::FileTypes::Type::PIN)
        .value("PARAMXML", OpenMS::FileTypes::Type::PARAMXML)
        .value("SPLIB", OpenMS::FileTypes::Type::SPLIB)
        .export_values();''',
    },
    # Nested enum bindings that must be added AFTER the containing class is bound
    "__post_class_enums__": {
        "SpectrumSettings_SpectrumType": '''
    // Bind nested SpectrumType enum to existing SpectrumSettings class
    {
        // Get the SpectrumSettings type from the module (already bound above)
        nb::handle spectrum_settings_type = m.attr("SpectrumSettings");
        nb::enum_<OpenMS::SpectrumSettings::SpectrumType>(spectrum_settings_type, "SpectrumType")
            .value("UNKNOWN", OpenMS::SpectrumSettings::SpectrumType::UNKNOWN)
            .value("CENTROID", OpenMS::SpectrumSettings::SpectrumType::CENTROID)
            .value("PROFILE", OpenMS::SpectrumSettings::SpectrumType::PROFILE)
            .export_values();
    }''',
    },
    "FloatDataArray": {
        "__init__": '''
        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::FloatDataArray&>(), "Copy constructor")''',
        "setName": '''
        .def("setName", [](OpenMS::DataArrays::FloatDataArray& self, const OpenMS::String& name) {
            self.setName(name);
        }, "name"_a, "Set the name")''',
        "getName": '''
        .def("getName", [](const OpenMS::DataArrays::FloatDataArray& self) {
            return self.getName();
        }, "Get the name")''',
        "get_data": '''
        .def("get_data", [](const OpenMS::DataArrays::FloatDataArray& self) {
            // Return as numpy array - copy data into new array
            size_t n = self.size();
            float* data = new float[n];
            std::copy(self.begin(), self.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
        }, "Returns a copy of the data as numpy array")''',
        "get_data_mv": '''
        .def("get_data_mv", [](nb::object self_obj) -> std::optional<nb::ndarray<nb::numpy, float, nb::ndim<1>>> {
            // Return a writable memory view (zero-copy)
            auto& self = nb::cast<OpenMS::DataArrays::FloatDataArray&>(self_obj);
            if (self.empty()) return std::nullopt;
            // Create view with self_obj as owner to keep the array alive and writable
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                self.data(), {self.size()}, self_obj
            );
        }, "Returns a view of the data as numpy array (fast but unsafe, None if empty)")''',
        "set_data": '''
        .def("set_data", [](OpenMS::DataArrays::FloatDataArray& self, std::vector<float> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")''',
        "push_back": '''
        .def("push_back", [](OpenMS::DataArrays::FloatDataArray& self, float val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")''',
        "resize": '''
        .def("resize", [](OpenMS::DataArrays::FloatDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")''',
        "clear": '''
        .def("clear", [](OpenMS::DataArrays::FloatDataArray& self) {
            self.clear();
        }, "Clear the array")''',
    },
    "IntegerDataArray": {
        "__init__": '''
        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::IntegerDataArray&>(), "Copy constructor")''',
        "get_data": '''
        .def("get_data", [](const OpenMS::DataArrays::IntegerDataArray& self) {
            std::vector<OpenMS::Int> arr(self.begin(), self.end());
            return arr;
        }, "Returns a copy of the data as a list")''',
        "set_data": '''
        .def("set_data", [](OpenMS::DataArrays::IntegerDataArray& self, std::vector<OpenMS::Int> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")''',
        "push_back": '''
        .def("push_back", [](OpenMS::DataArrays::IntegerDataArray& self, OpenMS::Int val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")''',
        "resize": '''
        .def("resize", [](OpenMS::DataArrays::IntegerDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")''',
        "clear": '''
        .def("clear", [](OpenMS::DataArrays::IntegerDataArray& self) {
            self.clear();
        }, "Clear the array")''',
    },
    "StringDataArray": {
        "__init__": '''
        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::StringDataArray&>(), "Copy constructor")''',
        "get_data": '''
        .def("get_data", [](const OpenMS::DataArrays::StringDataArray& self) {
            std::vector<OpenMS::String> arr(self.begin(), self.end());
            return arr;
        }, "Returns a copy of the data as a list")''',
        "set_data": '''
        .def("set_data", [](OpenMS::DataArrays::StringDataArray& self, std::vector<OpenMS::String> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")''',
        "push_back": '''
        .def("push_back", [](OpenMS::DataArrays::StringDataArray& self, const OpenMS::String& val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")''',
        "resize": '''
        .def("resize", [](OpenMS::DataArrays::StringDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")''',
        "clear": '''
        .def("clear", [](OpenMS::DataArrays::StringDataArray& self) {
            self.clear();
        }, "Clear the array")''',
    },
    "MSExperiment": {
        "getSpectra": '''
        .def("getSpectra", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSSpectrum>& {
            return self.getSpectra();
        }, nb::rv_policy::reference_internal, "Returns reference to spectra vector")''',
        "getChromatograms": '''
        .def("getChromatograms", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSChromatogram>& {
            return self.getChromatograms();
        }, nb::rv_policy::reference_internal, "Returns reference to chromatograms vector")''',
        "getSpectrum": '''
        .def("getSpectrum", [](OpenMS::MSExperiment& self, size_t index) -> OpenMS::MSSpectrum& {
            return self.getSpectrum(index);
        }, "index"_a, nb::rv_policy::reference_internal, "Returns reference to spectrum at index")''',
        "getChromatogram": '''
        .def("getChromatogram", [](OpenMS::MSExperiment& self, size_t index) -> OpenMS::MSChromatogram& {
            return self.getChromatogram(index);
        }, "index"_a, nb::rv_policy::reference_internal, "Returns reference to chromatogram at index")''',
        "addSpectrum": '''
        .def("addSpectrum", [](OpenMS::MSExperiment& self, const OpenMS::MSSpectrum& spec) {
            self.addSpectrum(spec);
        }, "spectrum"_a, "Add a spectrum to the experiment")''',
        "addChromatogram": '''
        .def("addChromatogram", [](OpenMS::MSExperiment& self, const OpenMS::MSChromatogram& chrom) {
            self.addChromatogram(chrom);
        }, "chromatogram"_a, "Add a chromatogram to the experiment")''',
        "getNrSpectra": '''
        .def("getNrSpectra", [](const OpenMS::MSExperiment& self) {
            return self.getNrSpectra();
        }, "Returns the number of spectra")''',
        "getNrChromatograms": '''
        .def("getNrChromatograms", [](const OpenMS::MSExperiment& self) {
            return self.getNrChromatograms();
        }, "Returns the number of chromatograms")''',
        "size": '''
        .def("size", [](const OpenMS::MSExperiment& self) {
            return self.size();
        }, "Returns the number of spectra")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::MSExperiment& self) {
            return self.size();
        })''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::MSExperiment& self, size_t i) -> OpenMS::MSSpectrum& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "__iter__": '''
        .def("__iter__", [](OpenMS::MSExperiment& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "MSExperiment_iter", self.begin(), self.end());
        })''',
    },
    "PeptideIdentification": {
        "getHits": '''
        .def("getHits", [](OpenMS::PeptideIdentification& self) -> std::vector<OpenMS::PeptideHit>& {
            return self.getHits();
        }, nb::rv_policy::reference_internal, "Returns reference to peptide hits")''',
        "setHits": '''
        .def("setHits", [](OpenMS::PeptideIdentification& self, const std::vector<OpenMS::PeptideHit>& hits) {
            self.setHits(hits);
        }, "hits"_a, "Set the peptide hits")''',
        "insertHit": '''
        .def("insertHit", [](OpenMS::PeptideIdentification& self, const OpenMS::PeptideHit& hit) {
            self.insertHit(hit);
        }, "hit"_a, "Insert a peptide hit")''',
    },
    "ProteinIdentification": {
        "getHits": '''
        .def("getHits", [](OpenMS::ProteinIdentification& self) -> std::vector<OpenMS::ProteinHit>& {
            return self.getHits();
        }, nb::rv_policy::reference_internal, "Returns reference to protein hits")''',
        "setHits": '''
        .def("setHits", [](OpenMS::ProteinIdentification& self, const std::vector<OpenMS::ProteinHit>& hits) {
            self.setHits(hits);
        }, "hits"_a, "Set the protein hits")''',
        "insertHit": '''
        .def("insertHit", [](OpenMS::ProteinIdentification& self, const OpenMS::ProteinHit& hit) {
            self.insertHit(hit);
        }, "hit"_a, "Insert a protein hit")''',
    },
    "Feature": {
        "getSubordinates": '''
        .def("getSubordinates", [](OpenMS::Feature& self) {
            return std::vector<OpenMS::Feature>(self.getSubordinates().begin(), self.getSubordinates().end());
        }, "Returns copy of subordinate features")''',
        "setSubordinates": '''
        .def("setSubordinates", [](OpenMS::Feature& self, const std::vector<OpenMS::Feature>& subs) {
            self.setSubordinates(subs);
        }, "subordinates"_a, "Set the subordinate features")''',
        "getPeptideIdentifications": '''
        .def("getPeptideIdentifications", [](OpenMS::Feature& self) {
            auto& list = self.getPeptideIdentifications();
            return std::vector<OpenMS::PeptideIdentification>(list.begin(), list.end());
        }, "Returns copy of peptide identifications")''',
        "setPeptideIdentifications": '''
        .def("setPeptideIdentifications", [](OpenMS::Feature& self, const std::vector<OpenMS::PeptideIdentification>& ids) {
            self.setPeptideIdentifications(ids);
        }, "ids"_a, "Set the peptide identifications")''',
    },
    "FeatureMap": {
        "getProteinIdentifications": '''
        .def("getProteinIdentifications", [](OpenMS::FeatureMap& self) {
            return std::vector<OpenMS::ProteinIdentification>(self.getProteinIdentifications().begin(), self.getProteinIdentifications().end());
        }, "Returns copy of protein identifications")''',
        "setProteinIdentifications": '''
        .def("setProteinIdentifications", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::ProteinIdentification>& ids) {
            self.setProteinIdentifications(ids);
        }, "ids"_a, "Set the protein identifications")''',
        "getUnassignedPeptideIdentifications": '''
        .def("getUnassignedPeptideIdentifications", [](OpenMS::FeatureMap& self) {
            auto& list = self.getUnassignedPeptideIdentifications();
            return std::vector<OpenMS::PeptideIdentification>(list.begin(), list.end());
        }, "Returns copy of unassigned peptide identifications")''',
        "setUnassignedPeptideIdentifications": '''
        .def("setUnassignedPeptideIdentifications", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::PeptideIdentification>& ids) {
            self.setUnassignedPeptideIdentifications(ids);
        }, "ids"_a, "Set the unassigned peptide identifications")''',
        "size": '''
        .def("size", [](const OpenMS::FeatureMap& self) {
            return self.size();
        }, "Returns the number of features")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::FeatureMap& self) {
            return self.size();
        })''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::FeatureMap& self, size_t i) -> OpenMS::Feature& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "__iter__": '''
        .def("__iter__", [](OpenMS::FeatureMap& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "FeatureMap_iter", self.begin(), self.end());
        })''',
        "push_back": '''
        .def("push_back", [](OpenMS::FeatureMap& self, const OpenMS::Feature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a feature to the map")''',
    },
    "ConsensusFeature": {
        "getPeptideIdentifications": '''
        .def("getPeptideIdentifications", [](OpenMS::ConsensusFeature& self) {
            auto& list = self.getPeptideIdentifications();
            return std::vector<OpenMS::PeptideIdentification>(list.begin(), list.end());
        }, "Returns copy of peptide identifications")''',
        "setPeptideIdentifications": '''
        .def("setPeptideIdentifications", [](OpenMS::ConsensusFeature& self, const std::vector<OpenMS::PeptideIdentification>& ids) {
            self.setPeptideIdentifications(ids);
        }, "ids"_a, "Set the peptide identifications")''',
    },
    "ConsensusMap": {
        "getProteinIdentifications": '''
        .def("getProteinIdentifications", [](OpenMS::ConsensusMap& self) {
            return std::vector<OpenMS::ProteinIdentification>(self.getProteinIdentifications().begin(), self.getProteinIdentifications().end());
        }, "Returns copy of protein identifications")''',
        "setProteinIdentifications": '''
        .def("setProteinIdentifications", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::ProteinIdentification>& ids) {
            self.setProteinIdentifications(ids);
        }, "ids"_a, "Set the protein identifications")''',
        "getUnassignedPeptideIdentifications": '''
        .def("getUnassignedPeptideIdentifications", [](OpenMS::ConsensusMap& self) {
            auto& list = self.getUnassignedPeptideIdentifications();
            return std::vector<OpenMS::PeptideIdentification>(list.begin(), list.end());
        }, "Returns copy of unassigned peptide identifications")''',
        "setUnassignedPeptideIdentifications": '''
        .def("setUnassignedPeptideIdentifications", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::PeptideIdentification>& ids) {
            self.setUnassignedPeptideIdentifications(ids);
        }, "ids"_a, "Set the unassigned peptide identifications")''',
        "size": '''
        .def("size", [](const OpenMS::ConsensusMap& self) {
            return self.size();
        }, "Returns the number of consensus features")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::ConsensusMap& self) {
            return self.size();
        })''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::ConsensusMap& self, size_t i) -> OpenMS::ConsensusFeature& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "__iter__": '''
        .def("__iter__", [](OpenMS::ConsensusMap& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "ConsensusMap_iter", self.begin(), self.end());
        })''',
        "push_back": '''
        .def("push_back", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusFeature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a consensus feature to the map")''',
    },
    "AASequence": {
        "size": '''
        .def("size", [](const OpenMS::AASequence& self) {
            return self.size();
        }, "Returns the number of residues")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::AASequence& self) {
            return self.size();
        })''',
        "toString": '''
        .def("toString", [](const OpenMS::AASequence& self) {
            return self.toString();
        }, "Returns string representation")''',
        "toUnmodifiedString": '''
        .def("toUnmodifiedString", [](const OpenMS::AASequence& self) {
            return self.toUnmodifiedString();
        }, "Returns unmodified string representation")''',
        "getMonoWeight": '''
        .def("getMonoWeight", [](const OpenMS::AASequence& self) {
            return self.getMonoWeight();
        }, "Returns monoisotopic weight")''',
        "fromString": '''
        .def_static("fromString", [](const OpenMS::String& s) {
            return OpenMS::AASequence::fromString(s);
        }, "sequence"_a, "Create AASequence from string")''',
    },
    "Param": {
        "getValue": '''
        .def("getValue", [](const OpenMS::Param& self, const OpenMS::String& key) {
            return self.getValue(key);
        }, "key"_a, "Returns the value for a key")''',
        "setValue": '''
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value) {
            self.setValue(key, value);
        }, "key"_a, "value"_a, "Set a value for a key")''',
        "exists": '''
        .def("exists", [](const OpenMS::Param& self, const OpenMS::String& key) {
            return self.exists(key);
        }, "key"_a, "Check if a key exists")''',
        "size": '''
        .def("size", [](const OpenMS::Param& self) {
            return self.size();
        }, "Returns the number of parameters")''',
        "empty": '''
        .def("empty", [](const OpenMS::Param& self) {
            return self.empty();
        }, "Check if empty")''',
        "clear": '''
        .def("clear", [](OpenMS::Param& self) {
            self.clear();
        }, "Clear all parameters")''',
    },
    "MzMLFile": {
        "load": '''
        .def("load", [](OpenMS::MzMLFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load an mzML file into an MSExperiment")''',
        "store": '''
        .def("store", [](OpenMS::MzMLFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store an MSExperiment to an mzML file")''',
    },
    "IdXMLFile": {
        "load": '''
        .def("load", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an idXML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("store", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename,
                         const std::vector<OpenMS::ProteinIdentification>& proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an idXML file")''',
    },
    "FeatureXMLFile": {
        "load": '''
        .def("load", [](OpenMS::FeatureXMLFile& self, const OpenMS::String& filename, OpenMS::FeatureMap& map) {
            self.load(filename, map);
        }, "filename"_a, "map"_a, "Load a featureXML file")''',
        "store": '''
        .def("store", [](OpenMS::FeatureXMLFile& self, const OpenMS::String& filename, const OpenMS::FeatureMap& map) {
            self.store(filename, map);
        }, "filename"_a, "map"_a, "Store a FeatureMap to featureXML")''',
    },
    "ConsensusXMLFile": {
        "load": '''
        .def("load", [](OpenMS::ConsensusXMLFile& self, const OpenMS::String& filename, OpenMS::ConsensusMap& map) {
            self.load(filename, map);
        }, "filename"_a, "map"_a, "Load a consensusXML file")''',
        "store": '''
        .def("store", [](OpenMS::ConsensusXMLFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& map) {
            self.store(filename, map);
        }, "filename"_a, "map"_a, "Store a ConsensusMap to consensusXML")''',
    },
    "PeptideIdentificationList": {
        "size": '''
        .def("size", [](const OpenMS::PeptideIdentificationList& self) {
            return self.size();
        }, "Returns the number of peptide identifications")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::PeptideIdentificationList& self) {
            return self.size();
        })''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::PeptideIdentificationList& self, size_t i) -> OpenMS::PeptideIdentification& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "__iter__": '''
        .def("__iter__", [](OpenMS::PeptideIdentificationList& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "PeptideIdentificationList_iter", self.begin(), self.end());
        })''',
        "push_back": '''
        .def("push_back", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentification& id) {
            self.push_back(id);
        }, "id"_a, "Add a peptide identification")''',
        "clear": '''
        .def("clear", [](OpenMS::PeptideIdentificationList& self) {
            self.clear();
        }, "Clear all identifications")''',
    },
    "PepXMLFile": {
        "load": '''
        .def("load", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load a pepXML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("store", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to a pepXML file")''',
    },
    "MzIdentMLFile": {
        "load": '''
        .def("load", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an mzIdentML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("store", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an mzIdentML file")''',
    },
    "Mobilogram": {
        "size": '''
        .def("size", [](const OpenMS::Mobilogram& self) {
            return self.size();
        }, "Returns the number of peaks")''',
        "__len__": '''
        .def("__len__", [](const OpenMS::Mobilogram& self) {
            return self.size();
        })''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::Mobilogram& self, size_t i) -> OpenMS::MobilityPeak1D& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        "__iter__": '''
        .def("__iter__", [](OpenMS::Mobilogram& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "Mobilogram_iter", self.begin(), self.end());
        })''',
        "push_back": '''
        .def("push_back", [](OpenMS::Mobilogram& self, const OpenMS::MobilityPeak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak")''',
        "clear": '''
        .def("clear", [](OpenMS::Mobilogram& self) {
            self.clear();
        }, "Clear all peaks")''',
    },
    "File": {
        "exists": '''
        .def_static("exists", [](const OpenMS::String& file) {
            return OpenMS::File::exists(file);
        }, "file"_a, "Check if a file exists")''',
        "readable": '''
        .def_static("readable", [](const OpenMS::String& file) {
            return OpenMS::File::readable(file);
        }, "file"_a, "Check if a file is readable")''',
        "writable": '''
        .def_static("writable", [](const OpenMS::String& file) {
            return OpenMS::File::writable(file);
        }, "file"_a, "Check if a file is writable")''',
        "isDirectory": '''
        .def_static("isDirectory", [](const OpenMS::String& path) {
            return OpenMS::File::isDirectory(path);
        }, "path"_a, "Check if a path is a directory")''',
        "basename": '''
        .def_static("basename", [](const OpenMS::String& file) {
            return OpenMS::File::basename(file);
        }, "file"_a, "Get the basename of a file path")''',
        "path": '''
        .def_static("path", [](const OpenMS::String& file) {
            return OpenMS::File::path(file);
        }, "file"_a, "Get the directory path of a file")''',
        "absolutePath": '''
        .def_static("absolutePath", [](const OpenMS::String& file) {
            return OpenMS::File::absolutePath(file);
        }, "file"_a, "Get the absolute path")''',
        "getTempDirectory": '''
        .def_static("getTempDirectory", []() {
            return OpenMS::File::getTempDirectory();
        }, "Get the temp directory")''',
        "getUserDirectory": '''
        .def_static("getUserDirectory", []() {
            return OpenMS::File::getUserDirectory();
        }, "Get the user home directory")''',
    },
    "DataProcessing": {
        "setProcessingActions": '''
        .def("setProcessingActions", [](OpenMS::DataProcessing& self, const std::set<OpenMS::DataProcessing::ProcessingAction>& actions) {
            self.setProcessingActions(actions);
        }, "actions"_a, "Set processing actions")''',
    },
}


@dataclass
class ModuleContent:
    """Content for a single C++ module file."""
    includes: Set[str]
    forward_declarations: List[str]
    class_bindings: List[str]
    enum_bindings: List[str]
    post_class_enums: List[str]  # Enums that need to be bound after classes exist


class NanobindEmitterV2:
    """
    Generates nanobind C++ binding code from merged C++/pxd declarations.
    Uses accurate C++ type information from libclang.
    """

    def __init__(self, num_modules: int = 8, core_only: bool = True):
        self.num_modules = num_modules
        self.core_only = core_only  # Only bind CORE_CLASSES for now

        # Standard includes for all modules
        self._standard_includes = {
            "<nanobind/nanobind.h>",
            "<nanobind/operators.h>",
            "<nanobind/stl/string.h>",
            "<nanobind/stl/vector.h>",
            "<nanobind/stl/map.h>",
            "<nanobind/stl/set.h>",
            "<nanobind/stl/shared_ptr.h>",
            "<nanobind/stl/optional.h>",
            "<nanobind/ndarray.h>",
            "<nanobind/make_iterator.h>",
            "<sstream>",       # For std::ostringstream in __repr__
            "<iomanip>",       # For std::setprecision in __repr__
            '"all_casters.h"',
        }

        # Operator name mappings
        self._operator_map = {
            "operator==": "__eq__",
            "operator!=": "__ne__",
            "operator<": "__lt__",
            "operator<=": "__le__",
            "operator>": "__gt__",
            "operator>=": "__ge__",
            "operator+": "__add__",
            "operator-": "__sub__",
            "operator*": "__mul__",
            "operator/": "__truediv__",
            "operator+=": "__iadd__",
            "operator-=": "__isub__",
            "operator*=": "__imul__",
            "operator/=": "__itruediv__",
            "operator[]": "__getitem__",
            "operator()": "__call__",
        }

    def emit(
        self,
        classes: Dict[str, MergedClass],
        output_dir: Path,
    ) -> None:
        """
        Generate nanobind C++ files.

        Parameters
        ----------
        classes : dict
            Merged class declarations with accurate C++ info.
        output_dir : Path
            Output directory for generated files.
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Distribute classes across modules
        modules = self._distribute_classes(classes)

        # Generate module files
        for i, module_classes in enumerate(modules, 1):
            content = self._generate_module_content(module_classes, classes)
            output_file = output_dir / f"module_{i}.cpp"
            self._write_module_file(output_file, i, content)

        # Generate main module file
        self._write_main_module(output_dir / "main_module.cpp", len(modules))

    def _distribute_classes(
        self, classes: Dict[str, MergedClass]
    ) -> List[List[MergedClass]]:
        """Distribute classes across modules for parallel compilation.

        Classes in the same inheritance chain are kept together in the same module
        to ensure proper nanobind base class registration during module init.
        """
        # Build inheritance groups using union-find
        parent: Dict[str, str] = {}  # Union-find parent pointers

        def find(x: str) -> str:
            """Find root of the set containing x."""
            if x not in parent:
                parent[x] = x
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: str, y: str) -> None:
            """Merge sets containing x and y."""
            px, py = find(x), find(y)
            if px != py:
                # Always make the alphabetically first one the parent for determinism
                if px < py:
                    parent[py] = px
                else:
                    parent[px] = py

        caster_types = get_caster_owned_types()

        # First pass: Identify abstract classes that are needed as base classes
        # These need to be bound (without constructors) before their derived classes
        abstract_base_classes_needed: Set[str] = set()
        for merged_class in classes.values():
            class_name = merged_class.name
            # Only check non-abstract CORE_CLASSES
            if self.core_only and class_name not in CORE_CLASSES:
                continue
            if class_name in SKIP_CLASSES or class_name in caster_types:
                continue
            if "::" in class_name:
                continue
            if merged_class.is_abstract:
                continue
            # Check base classes
            for base in getattr(merged_class, 'base_classes', []):
                if '::' in base:
                    base_name = base.split('::')[-1]
                else:
                    base_name = base
                if base.startswith('std::'):
                    continue
                # If base class is abstract and in CORE_CLASSES, we need it
                if base_name in CORE_CLASSES and base_name not in SKIP_CLASSES and base_name not in caster_types:
                    if base_name in classes:
                        if classes[base_name].is_abstract:
                            abstract_base_classes_needed.add(base_name)

        # Process all classes and union related inheritance chains
        for merged_class in classes.values():
            class_name = merged_class.name

            # Skip classes that won't be bound
            if self.core_only and class_name not in CORE_CLASSES:
                continue
            if class_name in SKIP_CLASSES or class_name in caster_types:
                continue
            if "::" in class_name:  # Skip nested classes
                continue
            # Skip abstract classes UNLESS they're needed as base classes
            if merged_class.is_abstract and class_name not in abstract_base_classes_needed:
                continue

            # Initialize this class in union-find
            find(class_name)

            # Union with base classes that will also be bound
            for base in getattr(merged_class, 'base_classes', []):
                if '::' in base:
                    base_name = base.split('::')[-1]
                else:
                    base_name = base

                # Skip std:: types
                if base.startswith('std::'):
                    continue

                # Check if base class will be bound
                if (base_name in CORE_CLASSES and
                    base_name not in SKIP_CLASSES and
                    base_name not in caster_types and
                    base_name in classes):
                    union(class_name, base_name)

        # Group classes by their root
        groups: Dict[str, List[MergedClass]] = {}
        for merged_class in classes.values():
            class_name = merged_class.name

            # Apply same filters as above
            if self.core_only and class_name not in CORE_CLASSES:
                continue
            if class_name in SKIP_CLASSES or class_name in caster_types:
                continue
            if "::" in class_name:
                continue
            # Skip abstract classes UNLESS they're needed as base classes
            if merged_class.is_abstract and class_name not in abstract_base_classes_needed:
                continue

            root = find(class_name)
            if root not in groups:
                groups[root] = []
            groups[root].append(merged_class)

        # Distribute groups across modules, keeping each group in one module
        modules: List[List[MergedClass]] = [[] for _ in range(self.num_modules)]
        module_sizes = [0] * self.num_modules

        # Sort groups by size (largest first) for better load balancing
        sorted_groups = sorted(groups.values(), key=lambda g: -len(g))

        for group in sorted_groups:
            # Assign to the module with the smallest current size
            min_idx = module_sizes.index(min(module_sizes))
            modules[min_idx].extend(group)
            module_sizes[min_idx] += len(group)

        # Store abstract base classes needed for use in _generate_class_bindings
        self._abstract_base_classes_needed = abstract_base_classes_needed

        # Store the set of all class names that will be bound
        self._bound_class_names: Set[str] = set()
        for module in modules:
            for cls in module:
                self._bound_class_names.add(cls.name)

        return modules

    def _topological_sort_classes(self, classes: List[MergedClass]) -> List[MergedClass]:
        """Sort classes so that base classes come before derived classes.

        This ensures that when nanobind binds a derived class with a base class
        specification, the base class is already registered.
        """
        # Build a map of class name -> MergedClass
        class_map = {c.name: c for c in classes}
        class_names = set(class_map.keys())

        # Build dependency graph: derived -> [bases it depends on]
        dependencies: Dict[str, Set[str]] = {}
        caster_types = get_caster_owned_types()

        for cls in classes:
            deps = set()
            for base in getattr(cls, 'base_classes', []):
                # Extract unqualified name
                if '::' in base:
                    base_name = base.split('::')[-1]
                else:
                    base_name = base

                # Only track dependencies on classes we're binding in this module
                if (base_name in class_names and
                    base_name in CORE_CLASSES and
                    base_name not in SKIP_CLASSES and
                    base_name not in caster_types):
                    deps.add(base_name)

            dependencies[cls.name] = deps

        # Kahn's algorithm for topological sort
        # Start with classes that have no dependencies (base classes)
        in_degree = {name: len(deps) for name, deps in dependencies.items()}
        queue = [name for name, degree in in_degree.items() if degree == 0]
        result = []

        while queue:
            # Sort queue for deterministic order
            queue.sort()
            name = queue.pop(0)
            result.append(class_map[name])

            # Reduce in-degree for classes that depend on this one
            for other_name, deps in dependencies.items():
                if name in deps:
                    in_degree[other_name] -= 1
                    if in_degree[other_name] == 0:
                        queue.append(other_name)

        # If we didn't process all classes, there's a cycle (shouldn't happen with inheritance)
        # Just append remaining classes in their original order
        processed = {c.name for c in result}
        for cls in classes:
            if cls.name not in processed:
                result.append(cls)

        return result

    def _generate_module_content(
        self,
        module_classes: List[MergedClass],
        all_classes: Dict[str, MergedClass],
    ) -> ModuleContent:
        """Generate content for a module."""
        content = ModuleContent(
            includes=self._standard_includes.copy(),
            forward_declarations=[],
            class_bindings=[],
            enum_bindings=[],
            post_class_enums=[],
        )

        # Sort classes topologically so base classes are bound before derived classes
        sorted_module_classes = self._topological_sort_classes(module_classes)

        # Check if this module contains classes that need enum bindings
        class_names_in_module = {m.name for m in module_classes}

        # Add standalone enum bindings (DriftTimeUnit) to the module with MSSpectrum
        if "MSSpectrum" in class_names_in_module and self.core_only and "MSSpectrum" in CORE_CLASSES:
            if "__enums__" in SPECIAL_METHODS:
                for enum_name, enum_code in SPECIAL_METHODS["__enums__"].items():
                    content.enum_bindings.append(enum_code)

        # Add nested enum bindings (SpectrumType) to the module with SpectrumSettings
        if "SpectrumSettings" in class_names_in_module and self.core_only and "SpectrumSettings" in CORE_CLASSES:
            content.includes.add("<OpenMS/METADATA/SpectrumSettings.h>")
            if "__post_class_enums__" in SPECIAL_METHODS:
                for enum_name, enum_code in SPECIAL_METHODS["__post_class_enums__"].items():
                    content.post_class_enums.append(enum_code)

        for merged_class in sorted_module_classes:
            class_name = merged_class.name

            # Check if this class will actually be bound (avoid including unnecessary headers)
            # In core_only mode, only CORE_CLASSES are bound
            if self.core_only and class_name not in CORE_CLASSES:
                continue

            # Skip problematic classes and types with casters
            if class_name in SKIP_CLASSES or class_name in get_caster_owned_types():
                continue

            # Skip nested classes (but allow nested namespaces)
            if "::" in class_name:
                continue

            # Skip abstract classes UNLESS they're needed as base classes
            if merged_class.is_abstract:
                if not (hasattr(self, '_abstract_base_classes_needed') and
                        class_name in self._abstract_base_classes_needed):
                    continue

            # Generate class binding first - only add headers if successful
            # Pass the set of class names in this module for base class checking
            module_class_names = {m.name for m in sorted_module_classes}
            binding = self._generate_class_binding(merged_class, module_class_names)
            if binding:
                content.class_bindings.append(binding)

                # Add OpenMS header - extract relative path from full path
                header_file = merged_class.header_file
                if header_file:
                    # Extract OpenMS-relative path (e.g., OpenMS/KERNEL/Peak1D.h)
                    # Look for the last occurrence of "OpenMS/" to handle paths like
                    # /path/to/OpenMS/src/openms/include/OpenMS/KERNEL/Peak1D.h
                    if "OpenMS/" in header_file:
                        idx = header_file.rfind("/include/")
                        if idx != -1:
                            header = header_file[idx + len("/include/"):]
                        else:
                            idx = header_file.rfind("OpenMS/")
                            header = header_file[idx:]
                    else:
                        header = header_file
                    content.includes.add(f"<{header}>")

                # Add additional includes for this class if needed
                if class_name in ADDITIONAL_INCLUDES:
                    for inc in ADDITIONAL_INCLUDES[class_name]:
                        content.includes.add(inc)

        return content

    def _generate_class_binding(
        self, merged_class: MergedClass, module_class_names: Optional[Set[str]] = None
    ) -> Optional[str]:
        """Generate nanobind binding code for a class.

        Args:
            merged_class: The class to generate bindings for
            module_class_names: Set of class names in this module (for base class filtering)
        """
        lines = []

        class_name = merged_class.name
        qualified_name = merged_class.qualified_name

        # In core_only mode, only bind CORE_CLASSES
        if self.core_only and class_name not in CORE_CLASSES:
            return None

        # Skip problematic classes and types with casters
        if class_name in SKIP_CLASSES or class_name in get_caster_owned_types():
            return None

        # Skip nested classes (but allow nested namespaces like OpenMS::DataArrays)
        # Nested classes have :: in the class name itself, not just in the namespace
        if "::" in class_name:
            return None

        # Skip abstract classes UNLESS they're needed as base classes
        is_abstract_base = (merged_class.is_abstract and
                           hasattr(self, '_abstract_base_classes_needed') and
                           class_name in self._abstract_base_classes_needed)
        if merged_class.is_abstract and not is_abstract_base:
            return None

        # Handle base classes - nanobind supports multiple inheritance
        # We can only specify base classes that are also bound to nanobind AND
        # that are in the same module (cross-module base classes aren't registered yet at init time)
        bound_bases = self._get_bound_base_classes(merged_class, module_class_names)
        if bound_bases:
            # nanobind syntax: nb::class_<Derived, Base1, Base2>(m, "Derived")
            base_spec = ", " + ", ".join(bound_bases)
        else:
            base_spec = ""

        # Start class definition
        lines.append(f'    nb::class_<{qualified_name}{base_spec}>(m, "{class_name}")')

        # Add docstring
        if merged_class.doc:
            doc = self._escape_string(merged_class.doc)
            lines[-1] = lines[-1][:-1]  # Remove trailing paren
            lines[-1] += f', "{doc}")'

        # Generate constructors (skip for abstract classes - they can't be instantiated)
        if not is_abstract_base:
            for ctor in merged_class.constructors:
                ctor_code = self._generate_constructor(ctor, qualified_name)
                if ctor_code:
                    lines.append(f"        {ctor_code}")

        # Generate methods
        for merged_method in merged_class.methods:
            if merged_method.wrap_ignore:
                continue
            method_code = self._generate_method(merged_method, qualified_name, merged_class)
            if method_code:
                lines.append(f"        {method_code}")

        # Add explicit inherited method bindings for classes where we skipped base class
        # specification due to non-virtual destructor issues
        inherited_methods = self._get_explicit_inherited_methods(merged_class, qualified_name)
        for method_line in inherited_methods:
            lines.append(f"        {method_line}")

        # Add hash support
        if merged_class.wrap_hash:
            lines.append(f'        .def("__hash__", [](const {qualified_name}& self) {{ return std::hash<{qualified_name}>{{}}(self); }})')

        # Detect container/iterator traits from AST (with fallback to hardcoded lists)
        has_iterator = (
            merged_class.wrap_iter
            or self._has_iterator_methods(merged_class)
            or class_name in ITERABLE_CLASSES  # Fallback for non-libclang modes
        )
        has_size = (
            self._has_size_method(merged_class)
            or class_name in CONTAINER_CLASSES
            or class_name in VECTOR_BASED_CLASSES
        )
        # Detect vector inheritance from canonical base classes
        vector_elem_type = self._get_vector_element_type(merged_class)
        is_vector_based = vector_elem_type is not None

        # Add iterator support (detected from begin/end methods or wrap-iter directive)
        if has_iterator:
            lines.append(f'        .def("__iter__", []({qualified_name}& self) {{ return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "{class_name}_iter", self.begin(), self.end()); }})')

        # Add __len__ for container classes (detected from size() method)
        if has_size:
            lines.append(f'        .def("__len__", []({qualified_name}& self) {{ return self.size(); }})')

        # Add __getitem__ for vector-based classes (detected from std::vector inheritance)
        if is_vector_based:
            lines.append(f'        .def("__getitem__", []({qualified_name}& self, size_t i) -> {vector_elem_type}& {{ ')
            lines.append(f'            if (i >= self.size()) throw nb::index_error();')
            lines.append(f'            return self[i];')
            lines.append(f'        }}, nb::rv_policy::reference_internal)')

        # Add special methods for this class (get_peaks, set_peaks, etc.)
        if class_name in SPECIAL_METHODS:
            for method_name, method_code in SPECIAL_METHODS[class_name].items():
                lines.append(method_code)

        # Add __repr__ for key classes
        repr_code = self._generate_repr(class_name, qualified_name)
        if repr_code:
            lines.append(repr_code)

        # Close class definition
        lines.append("        ;")

        return "\n".join(lines)

    def _has_method(self, merged_class: MergedClass, method_name: str) -> bool:
        """Check if a class has a specific method by name."""
        for m in merged_class.methods:
            if m.cpp_method.name == method_name:
                return True
        return False

    def _has_size_method(self, merged_class: MergedClass) -> bool:
        """Check if a class has a size() method (container-like)."""
        return self._has_method(merged_class, "size")

    def _has_iterator_methods(self, merged_class: MergedClass) -> bool:
        """Check if a class has begin() and end() methods (iterable)."""
        return self._has_method(merged_class, "begin") and self._has_method(merged_class, "end")

    def _get_vector_element_type(self, merged_class: MergedClass) -> Optional[str]:
        """Detect if class inherits from std::vector and extract element type.

        Parses canonical base classes looking for std::vector<T> pattern.
        Returns the element type T if found, None otherwise.
        """
        import re

        # Check canonical base classes for std::vector inheritance
        canonical_bases = getattr(merged_class.cpp_class, 'canonical_base_classes', [])
        for base in canonical_bases:
            # Match std::vector<ElementType> or std::__1::vector<ElementType> (libc++)
            match = re.match(r'std::(?:__\d+::)?vector<([^,>]+)', base)
            if match:
                elem_type = match.group(1).strip()
                # Ensure OpenMS types are qualified
                if not elem_type.startswith('OpenMS::') and not elem_type.startswith('std::'):
                    if elem_type not in ('int', 'float', 'double', 'bool', 'char', 'size_t'):
                        elem_type = f'OpenMS::{elem_type}'
                return elem_type

        # Fallback to hardcoded map for compatibility with non-libclang modes
        return self._get_element_type_fallback(merged_class.name)

    def _get_element_type_fallback(self, class_name: str) -> Optional[str]:
        """Fallback element type map for non-libclang modes."""
        element_types = {
            "MSSpectrum": "OpenMS::Peak1D",
            "MSChromatogram": "OpenMS::ChromatogramPeak",
            "Mobilogram": "OpenMS::MobilityPeak1D",
            "FloatDataArray": "float",
            "IntegerDataArray": "OpenMS::Int",
            "StringDataArray": "OpenMS::String",
        }
        return element_types.get(class_name)

    def _get_element_type(self, class_name: str) -> Optional[str]:
        """Get the element type for a vector-based class (legacy method)."""
        element_types = {
            "MSSpectrum": "OpenMS::Peak1D",
            "MSChromatogram": "OpenMS::ChromatogramPeak",
            "Mobilogram": "OpenMS::MobilityPeak1D",
            "FloatDataArray": "float",
            "IntegerDataArray": "OpenMS::Int",
            "StringDataArray": "OpenMS::String",
        }
        return element_types.get(class_name)

    def _get_bound_base_classes(
        self, merged_class: MergedClass, module_class_names: Optional[Set[str]] = None
    ) -> List[str]:
        """Get the primary base class that is also bound to nanobind.

        nanobind only supports single inheritance in the template syntax:
        nb::class_<Derived, Base>(m, "Derived")

        For classes with multiple inheritance, we only specify the first (primary)
        base class that is bound. This means we pick the base class that provides
        the core data layout (usually the first in the inheritance list).

        IMPORTANT: For classes with multiple inheritance where the first base class
        has a NON-VIRTUAL destructor, we must NOT specify any base class. This is
        because nanobind's destructor handling relies on proper virtual destructor
        chains. Without a virtual destructor, the destruction will cause segfaults.

        Only returns base classes that are:
        - In CORE_CLASSES (and not in SKIP_CLASSES or have type casters)
        - Actually being bound (in _bound_class_names)

        Cross-module inheritance is supported through NB_DOMAIN "pyopenms" as long as
        modules are imported in the correct order (base class modules before derived).

        Args:
            merged_class: The class to check base classes for
            module_class_names: Set of class names in the same module (optional, not used)

        Returns a list with at most one qualified C++ name like ['OpenMS::BaseFeature'].
        """
        caster_types = get_caster_owned_types()

        # Get base classes from C++ parsing
        base_classes = getattr(merged_class, 'base_classes', [])

        # Classes with NON-VIRTUAL destructors - specifying these as base classes
        # in nanobind when combined with multiple inheritance causes segfaults
        # because the destructor chain is not properly virtual
        NONVIRTUAL_DESTRUCTOR_CLASSES = {
            "Peak1D", "Peak2D", "MetaInfoInterface", "UniqueIdInterface",
            "DocumentIdentifier", "MetaInfoDescription"
        }

        # For multiple inheritance, check if first base has non-virtual destructor
        if len(base_classes) > 1:
            first_base = base_classes[0]
            first_base_name = first_base.split('::')[-1] if '::' in first_base else first_base
            if first_base_name in NONVIRTUAL_DESTRUCTOR_CLASSES:
                # Skip base class specification to avoid destructor issues
                return []

        for base in base_classes:
            # Extract unqualified name from qualified name like 'OpenMS::BaseFeature'
            if '::' in base:
                base_name = base.split('::')[-1]
            else:
                base_name = base

            # Skip std:: types (like std::vector) - they're not bound as classes
            if base.startswith('std::'):
                continue

            # Check if this base class is bound (in CORE_CLASSES, not skipped, no caster)
            if base_name in CORE_CLASSES and base_name not in SKIP_CLASSES and base_name not in caster_types:
                # Also check that the base class is actually in our bound classes
                # (i.e., it wasn't skipped due to pxd parsing errors or other issues)
                if hasattr(self, '_bound_class_names') and base_name not in self._bound_class_names:
                    continue  # Base class not being bound, skip it

                # Cross-module inheritance is OK with NB_DOMAIN - modules are imported in order
                # (module 1, 2, ... N) so base classes from lower-numbered modules are registered
                # before derived classes from higher-numbered modules.

                # Use fully qualified name for C++ template
                # Ensure we have OpenMS:: prefix for all types
                if base.startswith('OpenMS::'):
                    return [base]  # Already fully qualified
                elif '::' in base:
                    # Partially qualified (e.g., Internal::XMLFile) - add OpenMS prefix
                    return [f'OpenMS::{base}']
                else:
                    return [f'OpenMS::{base_name}']

        return []  # No valid base class found

    def _get_explicit_inherited_methods(
        self, merged_class: MergedClass, qualified_name: str
    ) -> List[str]:
        """Generate explicit method bindings for inherited methods.

        When we skip base class specification in nanobind due to non-virtual destructor
        issues with multiple inheritance, the derived class loses access to base class
        methods. This method adds explicit bindings for the most commonly used inherited
        methods.

        This also handles transitive inheritance - for example:
        - RichPeak2D (Peak2D + MetaInfoInterface + UniqueIdInterface) - skips base
        - BaseFeature (inherits from RichPeak2D) - needs explicit methods
        - Feature (inherits from BaseFeature) - needs explicit methods

        Returns a list of .def() method binding strings.
        """
        base_classes = getattr(merged_class, 'base_classes', [])
        class_name = merged_class.name

        # Classes with non-virtual destructors
        NONVIRTUAL_DESTRUCTOR_CLASSES = {
            "Peak1D", "Peak2D", "MetaInfoInterface", "UniqueIdInterface",
            "DocumentIdentifier", "MetaInfoDescription"
        }

        # Classes that inherit from non-virtual destructor classes and need explicit methods
        # These are classes where we skip base class specification
        CLASSES_NEEDING_EXPLICIT_METHODS = {
            "RichPeak2D",        # Peak2D + MetaInfoInterface + UniqueIdInterface
            "BaseFeature",      # -> RichPeak2D
            "Feature",          # -> BaseFeature -> RichPeak2D
            "ConsensusFeature", # -> BaseFeature -> RichPeak2D
            "FeatureHandle",    # -> Peak2D
            "FeatureMap",       # MetaInfoInterface + ...
            "ConsensusMap",     # MetaInfoInterface + ...
            "ExperimentalSettings",  # MetaInfoInterface + DocumentIdentifier
        }

        # Check if this class needs explicit inherited methods
        needs_explicit = False

        # Case 1: Class has multiple inheritance with non-virtual first base
        if len(base_classes) > 1:
            first_base = base_classes[0]
            first_base_name = first_base.split('::')[-1] if '::' in first_base else first_base
            if first_base_name in NONVIRTUAL_DESTRUCTOR_CLASSES:
                needs_explicit = True

        # Case 2: Class inherits from a class that has skipped base
        if not needs_explicit:
            for base in base_classes:
                base_name = base.split('::')[-1] if '::' in base else base
                if base_name in CLASSES_NEEDING_EXPLICIT_METHODS:
                    needs_explicit = True
                    break

        if not needs_explicit:
            return []

        methods = []

        # Determine which method sets to include based on inheritance chain
        needs_peak2d = False
        needs_metainfo = False
        needs_uniqueid = False
        needs_docid = False

        # Check direct bases
        for base in base_classes:
            base_name = base.split('::')[-1] if '::' in base else base
            if base_name in {"Peak2D", "RichPeak2D", "BaseFeature", "Feature", "ConsensusFeature", "FeatureHandle"}:
                needs_peak2d = True
            if base_name in {"MetaInfoInterface", "RichPeak2D", "BaseFeature", "Feature",
                            "ConsensusFeature", "FeatureMap", "ConsensusMap", "ExperimentalSettings"}:
                needs_metainfo = True
            if base_name in {"UniqueIdInterface", "RichPeak2D", "BaseFeature", "Feature",
                            "ConsensusFeature", "FeatureMap", "ConsensusMap"}:
                needs_uniqueid = True
            if base_name in {"DocumentIdentifier", "ExperimentalSettings", "FeatureMap", "ConsensusMap"}:
                needs_docid = True

        # Also check for the specific classes
        if class_name in {"RichPeak2D", "BaseFeature", "Feature", "ConsensusFeature", "FeatureHandle"}:
            needs_peak2d = True
        if class_name in {"RichPeak2D", "BaseFeature", "Feature", "ConsensusFeature",
                         "FeatureMap", "ConsensusMap", "ExperimentalSettings"}:
            needs_metainfo = True
        if class_name in {"RichPeak2D", "BaseFeature", "Feature", "ConsensusFeature",
                         "FeatureMap", "ConsensusMap"}:
            needs_uniqueid = True
        if class_name in {"ExperimentalSettings", "FeatureMap", "ConsensusMap"}:
            needs_docid = True

        # Peak2D methods
        if needs_peak2d:
            methods.extend([
                f'.def("getRT", []({qualified_name}& self) {{ return self.getRT(); }}, "Returns the retention time")',
                f'.def("setRT", []({qualified_name}& self, double rt) {{ self.setRT(rt); }}, "rt"_a, "Sets the retention time")',
                f'.def("getMZ", []({qualified_name}& self) {{ return self.getMZ(); }}, "Returns the m/z")',
                f'.def("setMZ", []({qualified_name}& self, double mz) {{ self.setMZ(mz); }}, "mz"_a, "Sets the m/z")',
                f'.def("getIntensity", []({qualified_name}& self) {{ return self.getIntensity(); }}, "Returns the intensity")',
                f'.def("setIntensity", []({qualified_name}& self, float intensity) {{ self.setIntensity(intensity); }}, "intensity"_a, "Sets the intensity")',
            ])

        # MetaInfoInterface methods
        if needs_metainfo:
            methods.extend([
                f'.def("getMetaValue", []({qualified_name}& self, const OpenMS::String& name) {{ return self.getMetaValue(name); }}, "name"_a, "Returns the value corresponding to a string")',
                f'.def("setMetaValue", []({qualified_name}& self, const OpenMS::String& name, const OpenMS::DataValue& value) {{ self.setMetaValue(name, value); }}, "name"_a, "value"_a, "Sets the DataValue corresponding to a name")',
                f'.def("metaValueExists", []({qualified_name}& self, const OpenMS::String& name) {{ return self.metaValueExists(name); }}, "name"_a, "Returns whether an entry with the given name exists")',
                f'.def("removeMetaValue", []({qualified_name}& self, const OpenMS::String& name) {{ self.removeMetaValue(name); }}, "name"_a, "Removes the DataValue corresponding to name")',
                f'.def("isMetaEmpty", []({qualified_name}& self) {{ return self.isMetaEmpty(); }}, "Returns if the MetaInfo is empty")',
                f'.def("clearMetaInfo", []({qualified_name}& self) {{ self.clearMetaInfo(); }}, "Removes all meta values")',
            ])

        # UniqueIdInterface methods
        if needs_uniqueid:
            methods.extend([
                f'.def("getUniqueId", []({qualified_name}& self) {{ return self.getUniqueId(); }}, "Returns the unique id")',
                f'.def("setUniqueId", []({qualified_name}& self, uint64_t id) {{ self.setUniqueId(id); }}, "id"_a, "Sets the unique id")',
                f'.def("hasValidUniqueId", []({qualified_name}& self) {{ return self.hasValidUniqueId(); }}, "Returns whether the unique id is valid")',
                f'.def("clearUniqueId", []({qualified_name}& self) {{ self.clearUniqueId(); }}, "Clears the unique id")',
                f'.def("ensureUniqueId", []({qualified_name}& self) {{ return self.ensureUniqueId(); }}, "Assigns a valid unique id if the current one is invalid")',
            ])

        # DocumentIdentifier methods
        if needs_docid:
            methods.extend([
                f'.def("getIdentifier", []({qualified_name}& self) {{ return self.getIdentifier(); }}, "Returns the document identifier")',
                f'.def("setIdentifier", []({qualified_name}& self, const OpenMS::String& id) {{ self.setIdentifier(id); }}, "id"_a, "Sets the document identifier")',
            ])

        return methods

    def _generate_repr(self, class_name: str, qualified_name: str) -> Optional[str]:
        """Generate __repr__ method for a class."""
        # Define repr formats for key classes
        repr_formats = {
            "Peak1D": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<Peak1D mz=" + std::to_string(self.getMZ()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "Peak2D": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<Peak2D rt=" + std::to_string(self.getRT()) + " mz=" + std::to_string(self.getMZ()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "ChromatogramPeak": f'''        .def("__repr__", []({qualified_name}& self) {{
            return "<ChromatogramPeak rt=" + std::to_string(self.getRT()) + " intensity=" + std::to_string(self.getIntensity()) + ">";
        }})''',
            "MSSpectrum": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSSpectrum ms_level=" << self.getMSLevel()
                << " rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " num_peaks=" << self.size() << ">";
            return oss.str();
        }})''',
            "MSChromatogram": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSChromatogram native_id='" << self.getNativeID()
                << "' num_peaks=" << self.size() << ">";
            return oss.str();
        }})''',
            "MSExperiment": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<MSExperiment num_spectra=" << self.getNrSpectra()
                << " num_chromatograms=" << self.getNrChromatograms() << ">";
            return oss.str();
        }})''',
            "Feature": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<Feature rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " mz=" << std::setprecision(4) << self.getMZ()
                << " intensity=" << std::setprecision(0) << self.getIntensity()
                << " charge=" << self.getCharge() << ">";
            return oss.str();
        }})''',
            "FeatureMap": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<FeatureMap num_features=" << self.size() << ">";
            return oss.str();
        }})''',
            "ConsensusFeature": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<ConsensusFeature rt=" << std::fixed << std::setprecision(2) << self.getRT()
                << " mz=" << std::setprecision(4) << self.getMZ()
                << " intensity=" << std::setprecision(0) << self.getIntensity()
                << " size=" << self.size() << ">";
            return oss.str();
        }})''',
            "ConsensusMap": f'''        .def("__repr__", []({qualified_name}& self) {{
            std::ostringstream oss;
            oss << "<ConsensusMap num_consensus_features=" << self.size() << ">";
            return oss.str();
        }})''',
        }
        return repr_formats.get(class_name)

    # Classes with private or deleted copy constructors - skip copy ctor binding
    PRIVATE_COPY_CTOR_CLASSES = {
        "ElementDB", "ModificationsDB", "CrossLinksDB", "EnzymesDB", "RibonucleotideDB",
        "ProteaseDB", "ResidueDB",  # Singleton databases
        "ProgressLogger",  # Non-copyable
        "UniqueIdGenerator",  # Protected/private constructors
        "IsobaricQuantitationMethod",  # Abstract class
        "GridBasedCluster",  # Complex constructors with nested types
    }

    def _generate_constructor(self, ctor: CppMethod, class_name: str = "") -> Optional[str]:
        """Generate constructor binding."""
        # Skip move constructors (&&)
        if len(ctor.parameters) == 1:
            param_type = ctor.parameters[0].type_str
            if "&&" in param_type:
                return None

            # Detect copy constructors: parameter type is the class name (with or without &)
            is_copy_ctor = (
                ctor.name == param_type or  # Exact match: ClassName(ClassName)
                ctor.name == param_type.rstrip(" &") or  # With reference: ClassName(ClassName &)
                param_type.startswith(f"const {ctor.name}")  # Const ref: ClassName(const ClassName &)
            )

            # Skip copy constructors for classes with private copy ctors
            if is_copy_ctor and class_name in self.PRIVATE_COPY_CTOR_CLASSES:
                return None

        if not ctor.parameters:
            return ".def(nb::init<>())"

        # Use canonical types from libclang when available
        param_types = []
        for p in ctor.parameters:
            normalized = self._normalize_type(
                p.type_str, canonical_type=getattr(p, 'canonical_type', '')
            )
            # Detect if this is a copy constructor parameter
            # .pxd files may use "ClassName", "ClassName &", or "const ClassName &"
            is_copy_param = (
                ctor.name == p.type_str or
                ctor.name == p.type_str.rstrip(" &") or
                p.type_str.startswith(f"const {ctor.name}")
            ) and "&&" not in p.type_str

            if is_copy_param:
                # Copy constructor - use const ref with full qualified name
                # class_name is now the qualified name (e.g., OpenMS::DataArrays::FloatDataArray)
                normalized = f"const {class_name} &"
            param_types.append(normalized)

        types_str = ", ".join(param_types)
        return f".def(nb::init<{types_str}>())"

    def _generate_method(
        self, merged_method: MergedMethod, qualified_name: str, merged_class: MergedClass
    ) -> Optional[str]:
        """Generate method binding."""
        method = merged_method.cpp_method
        method_name = merged_method.wrap_as or method.name
        class_name = merged_class.name

        # Skip methods inherited from private std::vector base
        if class_name in VECTOR_BASED_CLASSES and method.name in VECTOR_INHERITED_METHODS:
            return None

        # Skip specific methods with complex return types
        if class_name in SKIP_METHODS and method.name in SKIP_METHODS[class_name]:
            return None

        # Handle operators
        if method.name.startswith("operator"):
            py_op = self._operator_map.get(method.name)
            if py_op:
                return self._generate_operator_binding(method, qualified_name, py_op)
            return None

        # Skip static methods for now (need different binding approach)
        if method.is_static:
            return self._generate_static_method(method, qualified_name, method_name)

        # Regular methods
        return self._generate_regular_method(merged_method, qualified_name, method_name, merged_class)

    def _generate_regular_method(
        self,
        merged_method: MergedMethod,
        qualified_name: str,
        method_name: str,
        merged_class: MergedClass,
    ) -> Optional[str]:
        """Generate binding for a regular method using lambda wrapper.

        We use lambdas for all methods because:
        1. We can't reliably detect C++ overloads without libclang
        2. Lambdas let us explicitly call the method with correct signature
        3. This avoids all overload resolution issues
        """
        method = merged_method.cpp_method

        # nanobind doesn't support lambdas with more than 8 parameters
        # Since we include 'self', methods can have at most 7 parameters
        if len(method.parameters) > 7:
            return None

        # Check if method is overloaded by parameter count/types
        overloads = [m for m in merged_class.methods if m.cpp_method.name == method.name]
        is_overloaded = len(overloads) > 1

        # For const/non-const overloads with same signature, prefer const version
        if is_overloaded:
            # Check if this is a const/non-const pair (same params, different const-ness)
            same_params_overloads = [
                m for m in overloads
                if len(m.cpp_method.parameters) == len(method.parameters)
            ]
            if len(same_params_overloads) == 2:
                const_versions = [m for m in same_params_overloads if m.cpp_method.is_const]
                non_const_versions = [m for m in same_params_overloads if not m.cpp_method.is_const]
                if const_versions and non_const_versions:
                    # This is a const/non-const pair - skip the non-const version
                    if not method.is_const:
                        return None

        # Build lambda parameters
        param_decls = []
        param_names_call = []
        param_names_arg = []

        # C++ reserved keywords that can't be used as parameter names
        cpp_keywords = {
            'int', 'char', 'float', 'double', 'void', 'bool', 'long', 'short',
            'unsigned', 'signed', 'const', 'static', 'class', 'struct', 'enum',
            'union', 'virtual', 'public', 'private', 'protected', 'new', 'delete',
            'return', 'if', 'else', 'for', 'while', 'do', 'switch', 'case',
            'break', 'continue', 'default', 'goto', 'try', 'catch', 'throw',
            'namespace', 'using', 'typedef', 'typename', 'template', 'auto',
            'register', 'extern', 'volatile', 'mutable', 'inline', 'explicit',
            'export', 'true', 'false', 'nullptr', 'this', 'operator', 'sizeof',
            'alignof', 'decltype', 'constexpr', 'noexcept', 'override', 'final',
        }

        for p in method.parameters:
            # Use canonical type from libclang when available (resolves typedefs)
            ptype = self._normalize_type(
                p.type_str,
                canonical_type=getattr(p, 'canonical_type', ''),
            )
            # Check if name is valid (not empty, not "arg*", not a C++ keyword)
            valid_name = (
                p.name
                and not p.name.startswith("arg")
                and p.name not in cpp_keywords
            )
            pname = p.name if valid_name else f"p{len(param_decls)}"
            param_decls.append(f"{ptype} {pname}")
            param_names_call.append(pname)
            if valid_name:
                param_names_arg.append(f'"{p.name}"_a')

        params_decl = ", ".join(param_decls)
        params_call = ", ".join(param_names_call)

        # Determine const qualifier for lambda capture
        if method.is_const:
            self_decl = f"const {qualified_name}& self"
        else:
            self_decl = f"{qualified_name}& self"

        # Build the lambda
        if params_decl:
            lambda_sig = f"[]({self_decl}, {params_decl})"
            call_expr = f"self.{method.name}({params_call})"
        else:
            lambda_sig = f"[]({self_decl})"
            call_expr = f"self.{method.name}()"

        # Build def call
        result = f'.def("{method_name}", {lambda_sig} {{ return {call_expr}; }}'
        if param_names_arg:
            result += f", {', '.join(param_names_arg)}"
        if merged_method.doc:
            doc = self._escape_string(merged_method.doc)
            result += f', "{doc}"'
        result += ")"

        return result

    def _generate_static_method(
        self, method: CppMethod, qualified_name: str, method_name: str
    ) -> Optional[str]:
        """Generate binding for a static method."""
        method_ptr = f"&{qualified_name}::{method.name}"
        return f'.def_static("{method_name}", {method_ptr})'

    def _generate_operator_binding(
        self, method: CppMethod, qualified_name: str, py_op: str
    ) -> str:
        """Generate binding for an operator."""
        op_name = method.name.replace("operator", "")

        if op_name == "[]":
            return f'.def("__getitem__", []({qualified_name}& self, size_t i) {{ return self[i]; }})'
        elif op_name in ["==", "!=", "<", "<=", ">", ">="]:
            return f'.def(nb::self {op_name} nb::self)'
        elif op_name in ["+", "-", "*", "/"]:
            return f'.def(nb::self {op_name} nb::self)'
        elif op_name in ["+=", "-=", "*=", "/="]:
            return f'.def(nb::self {op_name} nb::self)'

        return ""

    def _is_resolved_type(self, canonical_type: str, original_type: str) -> bool:
        """Check if canonical type represents a resolved typedef.

        Returns True if the canonical type differs meaningfully from the original,
        indicating libclang resolved a typedef to its underlying type.

        Used to detect cases like:
        - Peak1D::IntensityType -> float
        - Peak1D::PositionType -> OpenMS::DPosition<1>
        - OpenMS::Int -> int
        """
        # Primitive types that libclang resolves to
        primitives = {
            'void', 'bool', 'char', 'short', 'int', 'long', 'long long',
            'unsigned char', 'unsigned short', 'unsigned int', 'unsigned long',
            'unsigned long long', 'float', 'double', 'long double',
        }

        # Check if canonical resolved to a primitive
        canonical_clean = canonical_type.replace('const ', '').replace(' &', '').replace('&', '').strip()
        if canonical_clean in primitives:
            return True

        # Check if the types differ (typedef was resolved)
        original_clean = original_type.replace('const ', '').replace(' &', '').replace('&', '').strip()
        if canonical_clean != original_clean:
            # Only use canonical if it looks like a valid resolved type
            # (not a template parameter or other incomplete type)
            if '::' in canonical_clean or canonical_clean in primitives:
                return True

        return False

    def _normalize_type(
        self,
        type_str: str,
        preserve_reference: bool = False,
        canonical_type: str = "",
    ) -> str:
        """Normalize C++ type string for use in bindings.

        Parameters
        ----------
        type_str : str
            The type string to normalize (original spelling)
        preserve_reference : bool
            If True, preserve reference qualifiers for return types in static_cast
        canonical_type : str
            Optional canonical type from libclang (typedefs resolved)
        """
        if not type_str:
            return "void"

        # If we have a canonical type from libclang that looks resolved
        # (i.e., it's a primitive or already fully qualified), use it directly
        # for nested class typedefs like Peak1D::IntensityType -> float
        if canonical_type and self._is_resolved_type(canonical_type, type_str):
            result = canonical_type
        else:
            result = type_str

        # Convert Cython template syntax to C++: libcpp_vector[T] -> std::vector<T>
        import re
        result = re.sub(r'libcpp_vector\[([^\]]+)\]', r'std::vector<\1>', result)
        result = re.sub(r'libcpp_map\[([^,]+),\s*([^\]]+)\]', r'std::map<\1, \2>', result)
        result = re.sub(r'libcpp_set\[([^\]]+)\]', r'std::set<\1>', result)
        result = re.sub(r'libcpp_pair\[([^,]+),\s*([^\]]+)\]', r'std::pair<\1, \2>', result)
        result = re.sub(r'shared_ptr\[([^\]]+)\]', r'std::shared_ptr<\1>', result)

        # OpenMS type aliases that need namespace
        openms_typedefs = {
            'Int': 'OpenMS::Int',
            'UInt': 'OpenMS::UInt',
            'Size': 'OpenMS::Size',
            'String': 'OpenMS::String',
            'StringList': 'OpenMS::StringList',
            'IntList': 'OpenMS::IntList',
            'DoubleList': 'OpenMS::DoubleList',
            'FloatDataArray': 'OpenMS::DataArrays::FloatDataArray',
            'IntegerDataArray': 'OpenMS::DataArrays::IntegerDataArray',
            'StringDataArray': 'OpenMS::DataArrays::StringDataArray',
            'FloatDataArrays': 'OpenMS::MSSpectrum::FloatDataArrays',
            'IntegerDataArrays': 'OpenMS::MSSpectrum::IntegerDataArrays',
            'StringDataArrays': 'OpenMS::MSSpectrum::StringDataArrays',
            'Param': 'OpenMS::Param',
            'DataValue': 'OpenMS::DataValue',
            'ParamValue': 'OpenMS::ParamValue',
            'DateTime': 'OpenMS::DateTime',
            'CVTerm': 'OpenMS::CVTerm',
            'AASequence': 'OpenMS::AASequence',
            'EmpiricalFormula': 'OpenMS::EmpiricalFormula',
            'ProteinIdentification': 'OpenMS::ProteinIdentification',
            'PeptideIdentification': 'OpenMS::PeptideIdentification',
            'PeptideHit': 'OpenMS::PeptideHit',
            'ProteinHit': 'OpenMS::ProteinHit',
            'FeatureHandle': 'OpenMS::FeatureHandle',
            'Precursor': 'OpenMS::Precursor',
            'Product': 'OpenMS::Product',
            'SourceFile': 'OpenMS::SourceFile',
            'ContactPerson': 'OpenMS::ContactPerson',
            'Sample': 'OpenMS::Sample',
            'Software': 'OpenMS::Software',
            'Instrument': 'OpenMS::Instrument',
            'IonSource': 'OpenMS::IonSource',
            'MassAnalyzer': 'OpenMS::MassAnalyzer',
            'IonDetector': 'OpenMS::IonDetector',
            'HPLC': 'OpenMS::HPLC',
            'Gradient': 'OpenMS::Gradient',
            'Digestion': 'OpenMS::Digestion',
            'Modification': 'OpenMS::Modification',
            'Tagging': 'OpenMS::Tagging',
            'Element': 'OpenMS::Element',
            'Residue': 'OpenMS::Residue',
            'ModificationsDB': 'OpenMS::ModificationsDB',
            'ResidueDB': 'OpenMS::ResidueDB',
            'ElementDB': 'OpenMS::ElementDB',
            'EnzymeDB': 'OpenMS::EnzymeDB',
            'DigestionEnzyme': 'OpenMS::DigestionEnzyme',
            'DigestionEnzymeProtein': 'OpenMS::DigestionEnzymeProtein',
            'DigestionEnzymeRNA': 'OpenMS::DigestionEnzymeRNA',
            'ResidueModification': 'OpenMS::ResidueModification',
            'ModificationDefinition': 'OpenMS::ModificationDefinition',
            'ModificationDefinitionsSet': 'OpenMS::ModificationDefinitionsSet',
            'ModifiedPeptideGenerator': 'OpenMS::ModifiedPeptideGenerator',
            'ProteaseDigestion': 'OpenMS::ProteaseDigestion',
            'RNaseDigestion': 'OpenMS::RNaseDigestion',
            'TheoreticalSpectrumGenerator': 'OpenMS::TheoreticalSpectrumGenerator',
            'PeakPickerHiRes': 'OpenMS::PeakPickerHiRes',
            'Normalizer': 'OpenMS::Normalizer',
            'SpectrumLookup': 'OpenMS::SpectrumLookup',
            'MRMMapping': 'OpenMS::MRMMapping',
            'MRMFeatureFinderScoring': 'OpenMS::MRMFeatureFinderScoring',
            'MRMTransitionGroupPicker': 'OpenMS::MRMTransitionGroupPicker',
            'TransitionTSVFile': 'OpenMS::TransitionTSVFile',
            'TransitionPQPFile': 'OpenMS::TransitionPQPFile',
            'TargetedExperiment': 'OpenMS::TargetedExperiment',
            'LightTargetedExperiment': 'OpenMS::LightTargetedExperiment',
            'ReactionMonitoringTransition': 'OpenMS::ReactionMonitoringTransition',
            'IncludeExcludeTarget': 'OpenMS::IncludeExcludeTarget',
            'TargetedExperimentHelper': 'OpenMS::TargetedExperimentHelper',
            # Note: SwathMap is in OpenSwath:: namespace, not OpenMS::
            'OpenSwathScoring': 'OpenMS::OpenSwathScoring',
            'ChromatogramExtractor': 'OpenMS::ChromatogramExtractor',
            'MzMLFile': 'OpenMS::MzMLFile',
            'MzXMLFile': 'OpenMS::MzXMLFile',
            'MzDataFile': 'OpenMS::MzDataFile',
            'FeatureXMLFile': 'OpenMS::FeatureXMLFile',
            'ConsensusXMLFile': 'OpenMS::ConsensusXMLFile',
            'IdXMLFile': 'OpenMS::IdXMLFile',
            'PepXMLFile': 'OpenMS::PepXMLFile',
            'ProtXMLFile': 'OpenMS::ProtXMLFile',
            'MascotXMLFile': 'OpenMS::MascotXMLFile',
            'FASTAFile': 'OpenMS::FASTAFile',
            'FASTAEntry': 'OpenMS::FASTAFile::FASTAEntry',
            # Constructor parameter types
            'CVMappingFile': 'OpenMS::CVMappingFile',
            'FeatureFinderIdentificationAlgorithm': 'OpenMS::FeatureFinderIdentificationAlgorithm',
            'FeatureFinderAlgorithmMetaboIdent': 'OpenMS::FeatureFinderAlgorithmMetaboIdent',
            'ConsensusMapNormalizerAlgorithmMedian': 'OpenMS::ConsensusMapNormalizerAlgorithmMedian',
            'ConsensusMapNormalizerAlgorithmQuantile': 'OpenMS::ConsensusMapNormalizerAlgorithmQuantile',
            # Note: IonSource, MassAnalyzer, IonDetector defined above
            'ControlledVocabulary': 'OpenMS::ControlledVocabulary',
            'CVMappingRule': 'OpenMS::CVMappingRule',
            'CVMappings': 'OpenMS::CVMappings',
            'ConsensusFeature': 'OpenMS::ConsensusFeature',
            'Peak1D': 'OpenMS::Peak1D',
            'Peak2D': 'OpenMS::Peak2D',
            'ChromatogramPeak': 'OpenMS::ChromatogramPeak',
            'MobilityPeak1D': 'OpenMS::MobilityPeak1D',
            'MSSpectrum': 'OpenMS::MSSpectrum',
            'MSChromatogram': 'OpenMS::MSChromatogram',
            'MSExperiment': 'OpenMS::MSExperiment',
            'PeakMap': 'OpenMS::PeakMap',
            # More constructor types
            'SpectrumAccessOpenMS': 'OpenMS::SpectrumAccessOpenMS',
            'SpectrumAccessOpenMSCached': 'OpenMS::SpectrumAccessOpenMSCached',
            'MapAlignmentEvaluationAlgorithmRecall': 'OpenMS::MapAlignmentEvaluationAlgorithmRecall',
            'UInt64': 'OpenMS::UInt64',
            'LightTransition': 'OpenSwath::LightTransition',
            'LightCompound': 'OpenSwath::LightCompound',
            'LightProtein': 'OpenSwath::LightProtein',
            'LightPeptide': 'OpenSwath::LightPeptide',
            'ChromExtractParams': 'OpenMS::ChromExtractParams',
        }

        # OpenMS classes that need namespace (only add if not already namespaced)
        for typedef, qualified in openms_typedefs.items():
            # Match the type when not already qualified
            pattern = r'\b(?<!OpenMS::)(?<!std::)' + typedef + r'\b'
            result = re.sub(pattern, qualified, result)

        # Convert DPosition2 to DPosition<2> (common alias in .pxd)
        result = re.sub(r'\bDPosition2\b', r'OpenMS::DPosition<2>', result)
        result = re.sub(r'\bDPosition1\b', r'OpenMS::DPosition<1>', result)

        # Handle common nested type aliases (convert to actual types)
        # NOTE: When libclang is used with canonical types, these are already resolved
        # (e.g., Peak1D::IntensityType -> float). This fallback is only needed for
        # Doxygen mode or when canonical types aren't available.
        # TODO: Remove this map once all modes use canonical types from libclang
        if not canonical_type:
            nested_types = {
                'PositionType': 'double',           # DPosition<N> reduces to coordinate type
                'CoordinateType': 'double',
                'IntensityType': 'float',           # Peak intensity is float in OpenMS
                'MassType': 'double',
            }
            for nested, actual in nested_types.items():
                # First, replace fully-qualified forms like OpenMS::Peak1D::PositionType -> double
                result = re.sub(r'OpenMS::\w+::' + nested + r'\b', actual, result)
                # Then replace any remaining unqualified forms
                result = re.sub(r'\b' + nested + r'\b', actual, result)

        # Handle enum types that need full qualification
        enum_types = {
            'SpectrumType': 'OpenMS::SpectrumSettings::SpectrumType',
            'SpectrumLevel': 'OpenMS::SpectrumLevel',
            'IMFormat': 'OpenMS::IMFormat',
            'DriftTimeUnit': 'OpenMS::DriftTimeUnit',
            'PeakType': 'OpenMS::Peak1D::PeakType',
            'DataType': 'OpenMS::DataValue::DataType',
            'Polarity': 'OpenMS::IonSource::Polarity',
            'InletType': 'OpenMS::IonSource::InletType',
            'IonizationMethod': 'OpenMS::IonSource::IonizationMethod',
            'AnalyzerType': 'OpenMS::MassAnalyzer::AnalyzerType',
            'ResolutionMethod': 'OpenMS::MassAnalyzer::ResolutionMethod',
            'ResolutionType': 'OpenMS::MassAnalyzer::ResolutionType',
            'ScanDirection': 'OpenMS::MassAnalyzer::ScanDirection',
            'ScanLaw': 'OpenMS::MassAnalyzer::ScanLaw',
            'ReflectronState': 'OpenMS::MassAnalyzer::ReflectronState',
            'Type': 'OpenMS::IonDetector::Type',
            'AcquisitionMode': 'OpenMS::IonDetector::AcquisitionMode',
            'DecoyType': 'OpenMS::TargetedExperimentHelper::RetentionTime::RTType',
            'TermSpecificity': 'OpenMS::ResidueModification::TermSpecificity',
        }
        for enum_type, qualified in enum_types.items():
            pattern = r'\b(?<!OpenMS::)(?<!::)' + enum_type + r'\b'
            result = re.sub(pattern, qualified, result)

        # Strip reference qualifiers unless we're preserving them for static_cast
        if not preserve_reference:
            result = result.replace("const ", "").replace(" const", "")
            result = result.replace("&", "").strip()

        return result

    def _escape_string(self, s: str) -> str:
        """Escape a string for use in C++ code."""
        return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")

    def _qualify_openms_types(self, type_str: str) -> str:
        """Add OpenMS:: prefix to unqualified OpenMS types.

        Handles both the main type and template parameters, e.g.:
        RangeManagerContainer< RangeMZ, RangeIntensity >
        -> OpenMS::RangeManagerContainer<OpenMS::RangeMZ, OpenMS::RangeIntensity>
        """
        import re

        # OpenMS types that need qualification when unqualified
        openms_types = {
            # Range types
            "RangeManagerContainer", "RangeManager",
            "RangeMZ", "RangeRT", "RangeIntensity", "RangeMobility",
            "RangeBase",
            # Settings types
            "SpectrumSettings", "ChromatogramSettings",
            "ExperimentalSettings", "AcquisitionInfo",
            # Data types
            "MetaInfoInterface", "MetaInfo",
            "Peak1D", "Peak2D", "ChromatogramPeak", "MobilityPeak1D",
            "DefaultParamHandler",
            # Other common types
            "String", "DataValue", "ParamValue",
        }

        result = type_str
        for t in openms_types:
            # Match the type when not already qualified with OpenMS::
            # Use word boundaries and negative lookbehind for ::
            pattern = r'(?<![:\w])' + t + r'(?![:\w])'
            result = re.sub(pattern, f'OpenMS::{t}', result)

        return result

    def _write_module_file(
        self, output_path: Path, module_num: int, content: ModuleContent
    ) -> None:
        """Write a module C++ file as a standalone NB_MODULE."""
        lines = []

        # Header comment
        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append(f"// Module {module_num}")
        lines.append("")

        # Includes
        for inc in sorted(content.includes):
            lines.append(f"#include {inc}")
        lines.append("")

        # Namespace
        lines.append("namespace nb = nanobind;")
        lines.append("using namespace nb::literals;")
        lines.append("")

        # Standalone NB_MODULE for this sub-module
        lines.append(f'NB_MODULE(_pyopenms2_{module_num}, m) {{')
        lines.append(f'    m.doc() = "pyOpenMS2 module {module_num}";')
        lines.append("")

        # Standalone enum bindings (must come before classes that use them)
        for binding in content.enum_bindings:
            lines.append(binding)
            lines.append("")

        # Class bindings
        for binding in content.class_bindings:
            lines.append(binding)
            lines.append("")

        # Post-class enum bindings (nested enums that need the class to exist first)
        for binding in content.post_class_enums:
            lines.append(binding)
            lines.append("")

        lines.append("}")

        output_path.write_text("\n".join(lines))

    def _write_main_module(self, output_path: Path, num_modules: int) -> None:
        """Write the main module file.

        Note: With standalone NB_MODULE per sub-module, this main module is
        optional. It's kept for backward compatibility but just creates an
        empty module - the actual classes are in _pyopenms2_N modules.
        """
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append("// Main module - placeholder, actual bindings in _pyopenms2_N modules")
        lines.append("")
        lines.append("#include <nanobind/nanobind.h>")
        lines.append("")
        lines.append("namespace nb = nanobind;")
        lines.append("")

        # Main module definition - empty, just for compatibility
        lines.append('NB_MODULE(_pyopenms2, m) {')
        lines.append('    m.doc() = "pyOpenMS2: Python bindings for OpenMS (nanobind)\\n"')
        lines.append('             "Note: Import pyopenms instead for all classes.";')
        lines.append("}")

        output_path.write_text("\n".join(lines))
