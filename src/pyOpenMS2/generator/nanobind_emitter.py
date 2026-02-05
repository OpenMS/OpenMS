"""
Nanobind C++ Code Emitter for pyOpenMS2 (v2 - using libclang info)

This module generates nanobind C++ binding code from merged C++/pxd information.
It uses accurate C++ type information from libclang for correct method signatures.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

from .cpp_parser import CppMethod, MergedClass, MergedMethod, _PXD_TO_CPP_TYPE

logger = logging.getLogger(__name__)


def _unqualified_name(name: str) -> str:
    """Extract unqualified name from a possibly qualified C++ name (e.g. 'OpenMS::Peak1D' -> 'Peak1D')."""
    return name.split('::')[-1] if '::' in name else name


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


# Map from OpenMS include subdirectory to domain module name.
# Classes are assigned to domain files (bind_<domain>.cpp) based on their header path.
DOMAIN_MAP = {
    "KERNEL": "kernel",
    "METADATA": "metadata",
    "CHEMISTRY": "chemistry",
    "ANALYSIS": "analysis",
    "FEATUREFINDER": "featurefinder",
    "FORMAT": "format",
    "PROCESSING": "processing",
    "DATASTRUCTURES": "datastructures",
    "MATH": "datastructures",
    "COMPARISON": "datastructures",
    "CONCEPT": "datastructures",
    "ML": "ml",
    "SYSTEM": "misc",
    "INTERFACES": "misc",
    "IONMOBILITY": "misc",
    "QC": "misc",
    "APPLICATIONS": "misc",
    "OPENSWATHALGO": "analysis",
}

# Ordered list of all domain names (determines file generation order)
DOMAIN_NAMES = [
    "kernel", "metadata", "chemistry", "analysis", "featurefinder",
    "format", "processing", "datastructures", "ml", "misc",
]

def _get_handwritten_class_domain(hw_info: dict) -> str:
    """Derive domain from a handwritten class's first include path."""
    includes = hw_info.get("includes", [])
    if includes:
        # Strip angle brackets: "<OpenMS/KERNEL/Peak1D.h>" -> "OpenMS/KERNEL/Peak1D.h"
        header = includes[0].strip("<>")
        return _get_domain_from_header(header)
    return "misc"


def _get_domain_from_header(header_file: str) -> str:
    """Extract domain name from an OpenMS header path.

    E.g. 'OpenMS/KERNEL/Peak1D.h' -> 'kernel'
         '/path/to/include/OpenMS/FORMAT/MzMLFile.h' -> 'format'
    """
    if not header_file:
        return "misc"

    # Find the OpenMS subdirectory in the path
    # Look for pattern: OpenMS/<SUBDIR>/...
    idx = header_file.rfind("OpenMS/")
    if idx == -1:
        return "misc"

    after_openms = header_file[idx + len("OpenMS/"):]
    # First path component is the subdirectory
    parts = after_openms.split("/")
    if len(parts) < 2:
        return "misc"

    subdir = parts[0]
    return DOMAIN_MAP.get(subdir, "misc")


# Methods to skip for specific classes (complex return types or signature mismatches)
#
# AUTO-DETECTION (via libclang, no need to list here):
# - Const/non-const overloads: detected and auto-resolved (prefer const version)
# - Overloaded methods: detected in overloaded_methods set
# - Methods using incomplete types: detected via uses_incomplete_type flag
#
# STILL NEED MANUAL LISTING:
# - Return type mismatches between pxd and C++ (handled via SPECIAL_METHODS)
# - Complex parameter types (pairs, nested types, etc.)
# - Iterator methods (begin/end) - special handling needed
# - Methods with > 7 parameters (nanobind limit)
SKIP_METHODS = {
    "MSSpectrum": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "getIMData",  # Returns pair
        "findHighestInWindow",  # Return type mismatch
        "select",  # Complex parameter types
    },
    "MSChromatogram": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
    },
    "AASequence": {
        "getFormula",  # Custom overload with defaults in SPECIAL_METHODS
        "begin", "end",  # Iterator types
        "getMZ", "getMonoWeight", "getAverageWeight",  # Custom overloads in SPECIAL_METHODS
        "fromString", "fromStringPermissive",  # Custom implementations in SPECIAL_METHODS
        "toBracketString",  # Need default parameter values
    },
    "Param": {
        "begin", "end",  # Iterator types
        "getDescription", "getTags",  # Complex return types
    },
    "TheoreticalSpectrumGenerator": {
        "getSpectrum",  # Need backward-compatible output parameter signature
    },
    "InternalCalibration": {
        "fillCalibrants",  # 7 parameters + self = 8, hits nanobind limit
    },
    "ModificationsDB": {
        "getAllSearchModifications",  # Complex return type
        "getModification",  # Overloaded
        "addModification",  # Takes unique_ptr
        "searchModificationsByDiffMonoMass",  # Uses TermSpecificity from module_8 - causes std::bad_cast
        "getBestModificationByDiffMonoMass",  # Uses TermSpecificity from module_8
        "searchModifications",  # Uses TermSpecificity from module_8
    },
    "ResidueModification": {
        "getDiffFormula",  # Complex return type
    },
    "EmpiricalFormula": {
        "begin", "end",  # Iterator types
        "getIsotopeDistribution",  # Takes IsotopePatternGenerator
        "getConditionalFragmentIsotopeDist",  # Takes CoarseIsotopePatternGenerator
    },
    "Residue": set(),  # getModification pointer return now auto-handled
    "Mobilogram": {
        "getFloatDataArrays", "setFloatDataArrays",
        "getIntegerDataArrays", "setIntegerDataArrays",
        "getStringDataArrays", "setStringDataArrays",
        "findNearest",  # Return type mismatch
        "select",  # Complex parameter types
    },
    "RankData": {
        "rankdata_double", "rankdata_float", "rankdata_int",  # Return numpy arrays
    },
    # MultipleTesting: Consolidated in later entry to avoid duplicate key issues
    "File": {
        "getUniqueName",  # Returns temp file path
    },
    "MzMLFile": {
        "load", "store", "loadBuffer", "storeBuffer", "isSemanticallyValid",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
    },
    "MzXMLFile": {
        "load", "store",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
    },
    "CachedmzML": {
        "load",  # Need backward-compatible signature with output parameter
    },
    "IdXMLFile": {
        "load", "store",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
    },
    "FeatureXMLFile": {
        "load", "store",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
    },
    "ConsensusXMLFile": {
        "load", "store",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
    },
    "MascotXMLFile": {
        "initializeLookup",  # Static method takes SpectrumMetaDataLookup by value, private copy ctor
        "load",  # Overload with SpectrumMetaDataLookup by value, private copy ctor
    },
    "PeakPickerHiRes": {
        "pickExperiment",  # Overload uses OnDiscMSExperiment (needs special handling)
    },
    "LinearResamplerAlign": {
        "raster",  # Inherited const overload conflicts with non-const override
    },
    "PepXMLFile": {
        "load", "store",  # Complex signature with PeptideIdentificationList
    },
    "MzIdentMLFile": {
        "load", "store",  # Complex signature with PeptideIdentificationList
    },
    "XQuestResultXMLFile": {
        "load", "store",  # Complex signature with PeptideIdentificationList
    },
    "DataFilters": {
        "passes",  # Multiple overloads including Feature/ConsensusFeature (forward-declared)
    },
    "ProteinIdentification": {
        "computeCoverage",  # Takes ConsensusMap/PeptideIdentificationList - complex types
        "setPrimaryMSRunPath",  # Overloaded with MSExperiment parameter
        "setSearchParameters", "getSearchParameters",  # SearchParameters nested type
        "insertProteinGroup", "insertIndistinguishableProteins",  # ProteinGroup nested type
    },
    "ConsensusFeature": {
        "computeDechargeConsensus",  # Takes forward-declared FeatureMap
    },
    "ProteinInference": {
        "run",  # Takes ConsensusMap (forward-declared)
    },
    # === High-priority classes being unblocked from SKIP_CLASSES ===
    "IDFilter": {
        "countHits",  # Overloaded (PeptideIdentificationList vs ProteinIdentification vector)
        "updateProteinGroups",  # ProteinGroup nested type
    },
    "FalseDiscoveryRate": {
        "applyBasic",  # ScoreToTgtDecLabelPairs unknown type
        "rocN",  # ScoreToTgtDecLabelPairs unknown type
        "applyEstimated",  # ScoreToTgtDecLabelPairs unknown type
        "trapezoidal_area_xEqy",  # ScoreToTgtDecLabelPairs unknown type
    },
    "FASTAFile": {
        "readNext",  # FASTAEntry nested type issue
        "writeNext",  # FASTAEntry nested type issue
    },
    "TraMLFile": {
        "load", "store",  # By-value params in auto-gen; SPECIAL_METHODS uses refs
        "isSemanticallyValid",  # Complex validator types
    },
    "TargetedExperiment": {
        "setInstruments",  # Namespace issue with vector<Instrument>
        "addInstrument",  # TargetedExperimentHelper::Instrument double-qualified
        "getInstruments",  # TargetedExperimentHelper namespace issue
        "setTargets",  # TargetedExperimentHelper namespace
        "getTargets",  # TargetedExperimentHelper namespace
        "addTarget",  # TargetedExperimentHelper namespace
        "setContacts",  # TargetedExperimentHelper::Contact namespace
        "getContacts",  # TargetedExperimentHelper::Contact namespace
        "addContact",  # TargetedExperimentHelper::Contact namespace
        "setPublications",  # TargetedExperimentHelper::Publication namespace
        "getPublications",  # TargetedExperimentHelper::Publication namespace
        "addPublication",  # TargetedExperimentHelper::Publication namespace
        "setInterpretations",  # TargetedExperimentHelper::Interpretation namespace
        "getInterpretations",  # TargetedExperimentHelper::Interpretation namespace
    },
    "MapAlignmentAlgorithmPoseClustering": {
        "align",  # pxd type mismatch
        "setReference",  # pxd type mismatch
    },
    "MapAlignmentTransformer": set(),  # transformRetentionTimes overload now auto-handled
    "TransformationDescription": {
        "fitModel",  # Lambda analysis / complex param types
    },
    "TransformationModelLinear": set(),  # getDefaultParameters output param now auto-handled
    "EmgGradientDescent": {"getDefaultParameters"},  # Custom overloads in SPECIAL_METHODS
    "TargetedSpectraExtractor": {"getDefaultParameters"},  # Custom overloads in SPECIAL_METHODS
    "InstrumentSettings": {
        "getScanWindows",  # ScanWindow incomplete type
        "setScanWindows",  # ScanWindow incomplete type
    },
    "Instrument": {
        "getIonSources",  # const/non-const overloads + incomplete types
        "setIonSources",  # IonSource incomplete type
        "getIonDetectors",  # const/non-const overloads + incomplete types
        "setIonDetectors",  # IonDetector incomplete type
        "getMassAnalyzers",  # const/non-const overloads + incomplete types
        "setMassAnalyzers",  # MassAnalyzer incomplete type
    },
    "ReactionMonitoringTransition": {
        "getProduct",  # Lambda analysis
        "setProduct",  # Lambda analysis
        "getIntermediateProduct",  # Lambda analysis
        "setIntermediateProduct",  # Lambda analysis
        "getPrediction",  # Lambda analysis
        "setPrediction",  # Lambda analysis
    },
    "PeakIntegrator": {
        "getDefaultParameters",  # Output parameter
    },
    "MRMFeatureQC": {
        "component_qcs",  # Nested type vector - handled in SPECIAL_METHODS
        "component_group_qcs",  # Nested type vector
        "component_group_pair_qcs",  # Nested type vector
    },
    "ProteaseDigestion": {
        "digest",  # pxd type mismatch
        "digestUnmodified",  # pxd type mismatch
        "isValidProduct",  # Need default parameter values
    },
    "ProteaseDB": {
        "getInstance",  # pxd has wrap-ignore, need manual binding
        "getAllNames",  # Need backward-compatible output param signature
        "getEnzyme",  # Returns pointer, need reference policy
    },
    "ElementDB": {
        "getElement",  # Overloaded, returns pointer - SPECIAL_METHODS uses suffixed keys
        "hasElement",  # Overloaded - SPECIAL_METHODS uses suffixed keys
        "getNames",  # Returns unordered_map of pointers
        "getSymbols",  # Returns unordered_map of pointers
        "getAtomicNumbers",  # Returns unordered_map of pointers
    },
    "RibonucleotideDB": {
        "getRibonucleotide",  # Returns shared_ptr - must use reference policy
        "getRibonucleotidePrefix",  # Returns shared_ptr - must use reference policy
    },
    "FileHandler": {
        "loadExperiment",  # Need overloads with default args
        "storeExperiment",  # Need overloads with default args
        "getType",  # pxd returns int, C++ returns FileTypes::Type enum
    },
    "MultipleTesting": {
        "pi0Est", "lfdr",  # Enum default args require enum to be registered first; handled in SPECIAL_METHODS
        # pi0MethodToString, toPi0Method, lfdrTransformToString, toLfdrTransform moved to SPECIAL_METHODS
    },
    # === Batch unblock: FLASHDeconv, FileHandler, IMTypes, etc. ===
    "IMTypes": set(),  # determineIMFormat overload now auto-handled
    "DeconvolvedSpectrum": {
        "setCharges",  # pxd type mismatch: vector<int> vs actual type
    },
    "FLASHDeconvAlgorithm": {
        "performSpectrumDeconvolution",  # Complex parameter types
        "calculateAveragine",  # Complex parameter types
    },
    "SpectralDeconvolution": set(),  # getIsotopeCosineAndIsoOffset output param now auto-handled
    "PeakGroup": {
        "getRepAbundance",  # pxd type mismatch
    },
    # === Unblocked classes: OnDiscMSExperiment, IndexedMzMLFileLoader, ExperimentalDesign, AcquisitionInfo, ExperimentalSettings ===
    "OnDiscMSExperiment": {
        "getExperimentalSettings",  # Returns shared_ptr<const ExperimentalSettings>
        "getMetaData",  # Returns shared_ptr<MSExperiment>
        "getSpectrumById",  # Returns shared_ptr<Spectrum> (interface type)
        "getChromatogramById",  # Returns shared_ptr<Chromatogram> (interface type)
    },
    "IndexedMzMLFileLoader": set(),  # store overload now auto-handled
    "ExperimentalDesign": {
        "getMSFileSection",  # Uses ExperimentalDesign_MSFileSectionEntry (unbound nested type)
        "setMSFileSection",  # Uses ExperimentalDesign_MSFileSectionEntry (unbound nested type)
        "getSampleSection",  # Uses ExperimentalDesign_SampleSection (unbound nested type)
        "setSampleSection",  # Uses ExperimentalDesign_SampleSection (unbound nested type)
        "fromConsensusMap",  # Static method
        "fromFeatureMap",  # Static method
        "fromIdentifications",  # Static method
    },
    "AcquisitionInfo": {
        "operator[]",  # Returns Acquisition& by index - needs special handling
    },
    # === Batch 2: More high-value classes unblocked from SKIP_CLASSES ===
    "MzTabFile": set(),  # store/load overloads now auto-handled
    "ControlledVocabulary": set(),  # getAllChildTerms output param now auto-handled
    "RNaseDB": set(),  # getInstance now auto-generated
    "Biosaur2Algorithm": set(),  # setMSData/getMSData overloads now auto-handled
    "FeatureFinderAlgorithmMetaboIdent": set(),  # setMSData overload now auto-handled
    "FeatureFinderIdentificationAlgorithm": {
        "getProgressLogger",  # wrap-ignored
    },
    "MascotGenericFile": {
        "updateMembers_",  # Protected method
    },
    "XQuestScores": set(),  # preScore overload now auto-handled
    "ExperimentalDesignFile": {"load"},  # TextFile overload uses unbound type
    "IDConflictResolverAlgorithm": set(),  # resolve/resolveBetweenFeatures overloads now auto-handled
    "SpectrumLookup": set(),  # extractScanNumber overload now auto-handled
    "RNaseDigestion": {
        "digest",  # Overloaded (with/without length params) + NASequence type
    },
    "CVMappings": {
        "setMappingRules",  # CVMappingRule incomplete type
        "getMappingRules",  # CVMappingRule incomplete type
        "addMappingRule",  # CVMappingRule incomplete type
    },
    "AccurateMassSearchResult": {
        "setIndividualIntensities",  # pxd type mismatch
    },
    "TransitionPQPFile": set(),  # overloads now auto-handled
    "TransitionTSVFile": set(),  # overload now auto-handled
    "NonNegativeLeastSquaresSolver": set(),  # solve overload now auto-handled
    "IncludeExcludeTarget": set(),  # replaceCVTerms overload now auto-handled
    "SeedListGenerator": {
        "generateSeedList",  # Overloaded (3 variants, one uses nested map[UInt64, ...])
        "convertSeedList",  # Overloaded (two variants)
    },
    "SequestInfile": {
        "getModifications",  # Returns nested map<String, vector<String>>
    },
    "SequestOutfile": {
        "getSequences",  # Complex nested types
    },
    "QcMLFile": {
        "map2csv",  # Nested map<String, map<String, String>>
    },
    "PeptideAndProteinQuant": {
        "readQuantData",  # Output param overloads have same Python signature after wrapping
    },
    "KDTreeFeatureMaps": {
        "addMaps",  # Overloaded (FeatureMap vs ConsensusMap)
        "applyTransformations",  # Uses TransformationModelLowess* (raw pointer)
    },
    "OPXLHelper": {
        "addProteinPositionMetaValues",  # Overloaded
        "enumerateCrossLinksAndMasses",  # Complex parameter types
        "digestDatabase",  # Complex parameter types
        "buildCandidates",  # Complex parameter types
        "buildFragmentAnnotations",  # Complex parameter types
        "buildPeptideIDs",  # Complex parameter types
        "collectPrecursorCandidates",  # Complex parameter types
        "isoPeakMeans",  # Complex parameter types
    },
    "MRMFeaturePickerFile": {
        "load",  # Uses nested type vectors (MRMFP_ComponentParams, MRMFP_ComponentGroupParams)
    },
    "MRMTransitionGroupPicker": {
        "pickTransitionGroup",  # Overloaded (LightTransition vs ReactionMonitoringTransition)
        "createMRMFeature",  # Complex template parameter
        "findLargestPeak",  # Output parameter (int& chr_idx, int& peak_idx)
        "findWidestPeakIndices",  # Output parameter
    },
    "AbsoluteQuantitationMethodFile": {
        "load",  # Uses vector<AbsoluteQuantitationMethod>
        "store",  # Uses vector<AbsoluteQuantitationMethod>
    },
    "AbsoluteQuantitationStandardsFile": {
        "load",  # Uses vector<AQS_runConcentration>
    },
    "OpenSwathDataAccessHelper": {
        "convertToSpectrumPtr",  # Returns shared_ptr<OSSpectrum>
        "convertToChromatogramPtr",  # Returns shared_ptr<OSChromatogram>
        "convertTargetedCompound",  # Uses Peptide (TargetedExperimentHelper nested type)
    },
    # === Batch 4: More classes unblocked from SKIP_CLASSES ===
    "OpenSwathHelper": {
        "checkSwathMapAndSelectTransitions",  # Complex OpenSwath parameter types
    },
    "MRMScoring": {
        "getMIMatrix",  # Returns Matrix<double> — template type
    },
    "IsobaricQuantifier": {
        "quantify",  # ConsensusMap parameter — may have forward-decl issues
    },
    "IsobaricIsotopeCorrector": {
        "correctIsotopicImpurities",  # Overloaded with pointer params
    },
    # === Batch 5: Final unblockable classes ===
    "OpenSwathOSWWriter": {
        "prepareLine",  # Uses LightCompound&, LightTransition* — complex types
    },
    # === Batch 6: Non-template classes that autowrap wraps ===
    "Compomer": {
        "getComponent",  # Returns complex nested map<String, vector<pair<...>>>
    },
    "SpectrumAlignmentScore": {
        "operator()",  # Both overloads wrap-ignored in pxd
    },
    "CrossLinksDB": {
        "getBestModificationByDiffMonoMass",  # Not in SPECIAL_METHODS for CrossLinksDB
    },
}

# Classes to skip due to incomplete type dependencies or other issues
# These will need manual handling or fixes in the C++ headers
#
# AUTO-DETECTION (via libclang, no need to list here):
# - Type casters: scanned from bindings/type_casters/*.h (use get_caster_owned_types())
# - Abstract classes: detected via pure virtual methods (is_abstract flag)
# - Deleted default constructors: detected via token analysis (has_deleted_default_constructor)
# - Private constructors: detected via access specifiers (has_private_constructor)
# - Const/non-const overloads: detected and handled in _generate_regular_method
# - Methods using incomplete types: detected via type.get_declaration().is_definition()
# - Qt base classes: auto-skipped in _get_bound_base_classes (QDate, QString, QObject, etc.)
#
# STILL NEED MANUAL LISTING:
# Pxd-to-C++ type mapping for template instance generation.
# Extends _PXD_TO_CPP_TYPE from cpp_parser with OpenMS class types needed
# for template parameter substitution (e.g., MSChromatogram -> OpenMS::MSChromatogram).
_TEMPLATE_PXD_TO_CPP = {
    **_PXD_TO_CPP_TYPE,
    "String": "OpenMS::String",
    "MSSpectrum": "OpenMS::MSSpectrum",
    "MSChromatogram": "OpenMS::MSChromatogram",
    "MRMFeature": "OpenMS::MRMFeature",
    "ReactionMonitoringTransition": "OpenMS::ReactionMonitoringTransition",
    "LightTransition": "OpenSwath::LightTransition",
    "RansacModelLinear": "OpenMS::Math::RansacModelLinear",
    "RansacModelQuadratic": "OpenMS::Math::RansacModelQuadratic",
    "Param": "OpenMS::Param",
    "Int64": "OpenMS::Int64",
    "Matrix": "OpenMS::Matrix",
    "Peak1D": "OpenMS::Peak1D",
    "Peak2D": "OpenMS::Peak2D",
    "ChromatogramPeak": "OpenMS::ChromatogramPeak",
}


# - Classes where libclang can't parse the header (missing includes, etc.)
# - Complex template instantiations
# - Parameter type mismatches between pxd and C++
# - Classes with constructor type mismatches
SKIP_CLASSES = {
    # Incomplete type issues (libclang can't fully resolve)
    "SwathFileConsumer",    # Template complexity
    "MSDataWritingConsumer", # Template complexity

    # Abstract classes (fallback - auto-detected via pxd comment or libclang, but these lack both)
    "FullSwathFileConsumer",        # Abstract but no pxd ABSTRACT comment
    "SpectrumAccessQuadMZTransforming",  # No pxd ABSTRACT comment
    "SpectrumAccessSqMass",              # No pxd ABSTRACT comment
    "SpectrumAccessOpenMSCached",        # No pxd ABSTRACT comment
    "IMSAlphabet",                  # Uses abstract IMSAlphabetParser
    "IMSAlphabetTextParser",        # Abstract template
    "BaseSuperimposer",             # No pxd ABSTRACT comment

    # Nested classes or special cases
    "XMLHandler",

    # Classes with complex constructors that reference unbound types
    "TransformationModel",          # Base class with DataPoints nested type — wrap-ignore in pxd
    "TransformationModelInterpolated",  # Copy ctor deleted (base inaccessible)

    # Classes where the type is actually a template alias (DPosition2 -> DPosition<2>)
    "DPosition2",

    # Classes with template parameters or complex constructors
    "OpenSwathWorkflowSonar",
    "OpenSwathWorkflow",
    "LightMRMTransitionGroup",

    # Template classes
    "MassDecomposer",               # template class - ABSTRACT, wrap-ignore in pxd

    # Classes with type alias issues
    "OSW_ChromExtractParams",

    # Classes with external type issues
    "IonMobilityScoring",           # Complex OpenSwath type dependencies

    # Classes provided as HANDWRITTEN_CLASSES (skip auto-generation)
    "OSSpectrum",                   # HANDWRITTEN_CLASS - OpenSwath data structure
    "OSChromatogram",               # HANDWRITTEN_CLASS - OpenSwath data structure
    "OSBinaryDataArray",            # HANDWRITTEN_CLASS - OpenSwath data structure
    "IsobaricChannelInformation",   # HANDWRITTEN_CLASS - nested type
    "TransformationModelBSpline",   # HANDWRITTEN_CLASS - dummy with static method
    "TransformationModelLowess",    # HANDWRITTEN_CLASS - dummy with static method
    "NASequence",                   # HANDWRITTEN_CLASS - minimal binding
    "SpectrumMetaDataLookup",       # HANDWRITTEN_CLASS - static methods
    "IsobaricChannelExtractor",     # HANDWRITTEN_CLASS
    "IsobaricNormalizer",           # HANDWRITTEN_CLASS
    "SpectrumAccessOpenMS",         # HANDWRITTEN_CLASS
    "SpectrumAccessOpenMSInMemory", # HANDWRITTEN_CLASS
    "SwathMap",                     # HANDWRITTEN_CLASS
    "DIAScoring",                   # HANDWRITTEN_CLASS
    "MRMFeatureFinderScoring",      # HANDWRITTEN_CLASS
    "OpenSwathScoring",             # HANDWRITTEN_CLASS
    "MRMFeature",                   # HANDWRITTEN_CLASS

    # No .pxd file in autowrap — not exposed
    "MSSim",
    "OpenSwathWorkflowBase",        # Nested enum / forward decl issues
    "NucleicAcidSearchEngine",
    "MSFraggerAdapter",
    "CometAdapter",
    "OMSSAAdapter",
    "MyriMatchAdapter",
    "PepNovoAdapter",

    # Miscellaneous
    "SignalToNoiseEstimator",       # ABSTRACT class - wrap-ignore in pxd
    "SteinScottImproveScore",       # No wrap-instances in pxd, no .pxd file
    "IMSWeights",                   # Not found in C++ headers (nested/helper class)
    "IntegerMassDecomposer",        # Uses IMSWeights which isn't bound

    # Incomplete types in members/parameters
    "CVMappingFile",                # Uses CVReference/CVMappingTerm (incomplete types)
    "CVMappingRule",                # Uses CVMappingTerm (incomplete type)
    "CVMappings",                   # Uses CVReference/CVMappingTerm (incomplete types)
    "MassExplainer",                # Uses Compomer (incomplete type in vector member)
    "ControlledVocabulary",         # Uses CVReference (incomplete type)
    "SemanticValidator",            # No default constructor (needs ControlledVocabulary&)
    "Date",                         # No default constructor (deleted)

    # Scoped enums not in OpenMS namespace
    "MzPAF",                        # MzPAFIonSeries/MzPAFDeltaUnit not in OpenMS:: namespace
    "SignalToNoiseEstimatorMedianChrom",  # Type alias, not real class

    # Const correctness issues
    "CrossLinksDB",                 # Const method mismatch

    # Constructor parameter issues (non-const lvalue ref from rvalue)
    "MSPGenericFile",               # Constructor takes MSExperiment by non-const lvalue ref

    # SQLite dependency issues
    "OSWFile",                      # Depends on SqliteConnector which needs SQLite headers
    "MzMLSqliteHandler",            # Depends on SqliteConnector which needs SQLite headers
    "MzMLSqliteSwathHandler",       # Depends on SqliteConnector which needs SQLite headers

    # Incomplete types in parameters / type caster issues
    "ILPDCWrapper",                 # Uses ChargePair (incomplete type)
    "FeatureDeconvolution",         # Uses map<String, Adduct> (type caster issue)
    "Compomer",                     # Uses map<String, Adduct> (type caster issue)
    "ChargePair",                   # Incomplete type (forward-declared)

}

# Classes with NON-VIRTUAL destructors - specifying these as base classes
# in nanobind when combined with multiple inheritance causes segfaults
NONVIRTUAL_DESTRUCTOR_CLASSES = {
    "Peak1D", "Peak2D", "MetaInfoInterface", "UniqueIdInterface",
    "DocumentIdentifier", "MetaInfoDescription"
}

# C++ reserved keywords that can't be used as parameter names
CPP_KEYWORDS = {
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

# Enums to skip (value name mismatches between pxd and C++, or other issues)
# NOTE: These are typically due to pxd files having different enum value names
# than the actual C++ headers (e.g., Precursor_ion vs Precursor)
SKIP_ENUMS = {
    "CHARGEMODE",       # In FeatureDeconvolution - pxd uses 'class' enum name
    "DriftTimeUnit",    # Already bound via SPECIAL_METHODS["__enums__"]
    "IMFormat",         # Already bound via SPECIAL_METHODS["__enums__"]
    "Method",           # Already bound via SPECIAL_METHODS["__enums__"] (RankData::Method)
    "NaNPolicy",        # Already bound via SPECIAL_METHODS["__enums__"] (RankData::NaNPolicy)
    "IntensityThresholdCalculation",  # References SignalToNoiseEstimatorMedianChrom (skipped class)
}

# Additional headers needed for specific classes
ADDITIONAL_INCLUDES = {
    "MSSpectrum": ["<OpenMS/KERNEL/Peak1D.h>", "<OpenMS/METADATA/DataArrays.h>", "<OpenMS/IONMOBILITY/IMTypes.h>", "<OpenMS/CONCEPT/ProgressLogger.h>", "<OpenMS/FORMAT/FileTypes.h>"],
    "MSChromatogram": ["<OpenMS/KERNEL/ChromatogramPeak.h>"],
    "MSExperiment": ["<OpenMS/KERNEL/MSSpectrum.h>", "<OpenMS/KERNEL/MSChromatogram.h>", "<OpenMS/KERNEL/MSExperiment.h>"],
    "Feature": ["<OpenMS/KERNEL/Feature.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "FeatureMap": ["<OpenMS/KERNEL/FeatureMap.h>", "<OpenMS/KERNEL/Feature.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "ConsensusFeature": ["<OpenMS/KERNEL/ConsensusFeature.h>", "<OpenMS/KERNEL/FeatureHandle.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "ConsensusMap": ["<OpenMS/KERNEL/ConsensusMap.h>", "<OpenMS/KERNEL/ConsensusFeature.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "PeptideIdentification": ["<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/PeptideHit.h>"],
    "ProteinIdentification": ["<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/ProteinHit.h>"],
    "PeptideHit": ["<OpenMS/METADATA/PeptideHit.h>"],
    "ProteinHit": ["<OpenMS/METADATA/ProteinHit.h>"],
    "PeptideEvidence": ["<OpenMS/METADATA/PeptideEvidence.h>"],
    "AASequence": ["<OpenMS/CHEMISTRY/AASequence.h>"],
    "EmpiricalFormula": ["<OpenMS/CHEMISTRY/EmpiricalFormula.h>", "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>", "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>"],
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
    "IDMapper": ["<OpenMS/METADATA/AnnotatedMSRun.h>"],
    "MRMTransitionGroup": ["<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>", "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "PosteriorErrorProbabilityModel": ["<OpenMS/FORMAT/TextFile.h>"],
    "IMTypes": ["<OpenMS/KERNEL/MSExperiment.h>", "<OpenMS/KERNEL/MSSpectrum.h>"],
    "Tagger": ["<OpenMS/KERNEL/MSSpectrum.h>"],
    "ModificationsDB": ["<OpenMS/CHEMISTRY/ModificationsDB.h>"],
    "ResidueModification": ["<OpenMS/CHEMISTRY/ResidueModification.h>"],
    "PeptideIdentificationList": ["<OpenMS/METADATA/PeptideIdentificationList.h>", "<OpenMS/METADATA/PeptideIdentification.h>"],
    "DataProcessing": ["<OpenMS/METADATA/DataProcessing.h>"],
    "InstrumentSettings": ["<OpenMS/METADATA/InstrumentSettings.h>"],
    "Mobilogram": ["<OpenMS/KERNEL/Mobilogram.h>", "<OpenMS/KERNEL/MobilityPeak1D.h>", "<OpenMS/METADATA/DataArrays.h>"],
    "BilinearInterpolation": ["<OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>"],
    "PrecalAveragine": ["<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>", "<OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>"],
    "LogMzPeak": ["<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>", "<OpenMS/KERNEL/Peak1D.h>"],
    "MassFeature_FDHS": ["<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>"],
    "IsobaricQuantities": ["<OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>"],
    "LightTransition": ["<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "LightModification": ["<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "LightCompound": ["<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "LightProtein": ["<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "LightTargetedExperiment": ["<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "ColumnHeader": ["<OpenMS/KERNEL/ConsensusMap.h>"],
    "ExtractionCoordinates": ["<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>"],
    "MRMFQC_ComponentQCs": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>"],
    "MRMFQC_ComponentGroupQCs": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>"],
    "MRMFQC_ComponentGroupPairQCs": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>"],
    "SelectorParameters": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>"],
    "PI_PeakArea": ["<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>"],
    "PI_PeakBackground": ["<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>"],
    "PI_PeakShapeMetrics": ["<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>"],
    "DBoundingBox2": ["<OpenMS/DATASTRUCTURES/DBoundingBox.h>"],
    "MSNumpressCoder": ["<OpenMS/FORMAT/MSNumpressCoder.h>"],
    "MultipleTesting": ["<OpenMS/MATH/STATISTICS/MultipleTesting.h>"],
    "File": ["<OpenMS/SYSTEM/File.h>"],
    "PepXMLFile": ["<OpenMS/FORMAT/PepXMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "MzIdentMLFile": ["<OpenMS/FORMAT/MzIdentMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "XQuestResultXMLFile": ["<OpenMS/FORMAT/XQuestResultXMLFile.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "DRange1": ["<OpenMS/DATASTRUCTURES/DRange.h>"],
    "DRange2": ["<OpenMS/DATASTRUCTURES/DRange.h>"],
    "DataFilters": ["<OpenMS/PROCESSING/MISC/DataFilters.h>", "<OpenMS/KERNEL/Feature.h>", "<OpenMS/KERNEL/ConsensusFeature.h>"],
    "ProteinInference": ["<OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>", "<OpenMS/KERNEL/ConsensusMap.h>"],
    "FASTAFile": ["<OpenMS/FORMAT/FASTAFile.h>"],
    "TraMLFile": ["<OpenMS/FORMAT/TraMLFile.h>", "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>"],
    "IDFilter": ["<OpenMS/PROCESSING/ID/IDFilter.h>", "<OpenMS/METADATA/PeptideIdentification.h>", "<OpenMS/METADATA/ProteinIdentification.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    "FalseDiscoveryRate": ["<OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>"],
    "TargetedExperiment": ["<OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>", "<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>"],
    "WindowMower": ["<OpenMS/PROCESSING/FILTERING/WindowMower.h>"],
    "Deisotoper": ["<OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>"],
    "ReactionMonitoringTransition": ["<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>"],
    "TransformationDescription": ["<OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>"],
    "ChromatogramExtractor": ["<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>", "<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>", "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>", "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>"],
    "ChromatogramExtractorAlgorithm": ["<OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>", "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>", "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "TransformationModelLinear": ["<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>"],
    "MapAlignmentAlgorithmPoseClustering": ["<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>"],
    "MapAlignmentTransformer": ["<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>"],
    "Instrument": ["<OpenMS/METADATA/Instrument.h>"],
    "ProteaseDigestion": ["<OpenMS/CHEMISTRY/ProteaseDigestion.h>"],
    "ProteaseDB": ["<OpenMS/CHEMISTRY/ProteaseDB.h>"],
    "RNaseDB": ["<OpenMS/CHEMISTRY/RNaseDB.h>", "<OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>"],
    "CrossLinksDB": ["<OpenMS/CHEMISTRY/CrossLinksDB.h>", "<OpenMS/CHEMISTRY/ResidueModification.h>"],
    "ProFormaParser": ["<OpenMS/CHEMISTRY/ProFormaParser.h>", "<OpenMS/CHEMISTRY/ProFormaData.h>", "<OpenMS/CHEMISTRY/AASequence.h>"],
    "IsobaricQuantitationMethod": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>"],
    "AbsoluteQuantitation": ["<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>", "<OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>"],
    "Peptide": ["<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>"],
    "ElementDB": ["<OpenMS/CHEMISTRY/ElementDB.h>", "<OpenMS/CHEMISTRY/Element.h>"],
    "ParamEntry": ["<OpenMS/DATASTRUCTURES/Param.h>"],
    "ParamNode": ["<OpenMS/DATASTRUCTURES/Param.h>"],
    "FeatureGroupingAlgorithmLabeled": ["<OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>", "<OpenMS/KERNEL/ConsensusMap.h>", "<OpenMS/KERNEL/FeatureMap.h>"],
    "OnDiscMSExperiment": ["<OpenMS/KERNEL/OnDiscMSExperiment.h>", "<OpenMS/KERNEL/MSExperiment.h>", "<OpenMS/KERNEL/MSSpectrum.h>", "<OpenMS/KERNEL/MSChromatogram.h>"],
    "IndexedMzMLFileLoader": ["<OpenMS/FORMAT/IndexedMzMLFileLoader.h>", "<OpenMS/KERNEL/OnDiscMSExperiment.h>", "<OpenMS/KERNEL/MSExperiment.h>"],
    "ExperimentalDesign": ["<OpenMS/METADATA/ExperimentalDesign.h>"],
    "ExperimentalSettings": ["<OpenMS/METADATA/ExperimentalSettings.h>", "<OpenMS/METADATA/Instrument.h>", "<OpenMS/METADATA/HPLC.h>", "<OpenMS/METADATA/Sample.h>", "<OpenMS/METADATA/ContactPerson.h>", "<OpenMS/METADATA/SourceFile.h>"],
    "AcquisitionInfo": ["<OpenMS/METADATA/AcquisitionInfo.h>", "<OpenMS/METADATA/Acquisition.h>"],
    "RibonucleotideDB": ["<OpenMS/CHEMISTRY/RibonucleotideDB.h>"],
    "PeakIntegrator": ["<OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>", "<OpenMS/KERNEL/MSChromatogram.h>", "<OpenMS/KERNEL/MSSpectrum.h>"],
    "ConvexHull2D": ["<OpenMS/DATASTRUCTURES/ConvexHull2D.h>", "<OpenMS/DATASTRUCTURES/DBoundingBox.h>"],
    "MRMFeatureQC": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>"],
    "OpenSwathOSWWriter": ["<OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>", "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>", "<OpenMS/KERNEL/FeatureMap.h>"],
    "OpenSwathHelper": ["<OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>"],
    "OSWFile": ["<OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>"],
    "MzMLSqliteHandler": ["<OpenMS/FORMAT/SqliteConnector.h>", "<OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>"],
    "MRMScoring": ["<OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>", "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    "TargetedExperimentHelper": ["<OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>"],
    "IsobaricQuantifier": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>", "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>"],
    "IsobaricIsotopeCorrector": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>", "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>"],
    "Compomer": ["<OpenMS/DATASTRUCTURES/Compomer.h>"],
    "MassExplainer": ["<OpenMS/DATASTRUCTURES/MassExplainer.h>"],
    "ILPDCWrapper": ["<OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>"],
    "KDTreeFeatureNode": ["<OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>"],
    "SpectrumAlignmentScore": ["<OpenMS/COMPARISON/SpectrumAlignmentScore.h>"],
    "IMSIsotopeDistribution": ["<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>"],
    "RealMassDecomposer": ["<OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>"],
    "DistanceMatrix": ["<OpenMS/DATASTRUCTURES/DistanceMatrix.h>"],
    "RANSAC": ["<OpenMS/ML/RANSAC/RANSAC.h>", "<OpenMS/ML/RANSAC/RANSACModelLinear.h>", "<OpenMS/ML/RANSAC/RANSACModelQuadratic.h>"],
    "LinearInterpolation": ["<OpenMS/ML/INTERPOLATION/LinearInterpolation.h>"],
    "SignalToNoiseEstimatorMedian": ["<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>", "<OpenMS/KERNEL/MSSpectrum.h>"],
    "SignalToNoiseEstimatorMeanIterative": ["<OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>", "<OpenMS/KERNEL/MSSpectrum.h>"],
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

# 
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
    "size", "push_back", "pop_back",
    "empty", "front", "back", "begin", "end", "cbegin", "cend",
    "at", "operator[]", "data", "capacity", "max_size", "shrink_to_fit",
    "insert", "erase", "emplace", "emplace_back", "assign", "swap",
}

# Methods with output parameters that need special handling.
# Maps method_name -> C++ lambda template for wrapping output parameter methods.
# The template uses {FULL_CLASS} as placeholder for the fully-qualified C++ class name.
# These replace the auto-generated binding with one that modifies the Python list in-place.
OUTPUT_PARAM_METHODS = {
    "getKeys": '''
        .def("getKeys", []({FULL_CLASS}& self, nb::list py_keys) {{
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {{
                py_keys.append(nb::cast(std::string(k)));
            }}
        }}, "keys"_a, "Fills the given list with all meta value keys")''',
}

# Primitive types for which reference returns are copied by value automatically
_PRIMITIVE_TYPES = {
    'int', 'unsigned int', 'long', 'unsigned long', 'long long',
    'unsigned long long', 'float', 'double', 'long double', 'bool',
    'char', 'unsigned char', 'short', 'unsigned short', 'size_t',
    'ptrdiff_t', 'int32_t', 'int64_t', 'uint32_t', 'uint64_t',
}

# Types with nanobind type casters that copy the value into a Python object.
# Returning a const-ref to these is equivalent to return-by-value, so
# rv_policy::reference_internal must NOT be applied (it would reference a temporary).
_TYPE_CASTER_TYPES = {
    'std::string', 'std::basic_string<char>',
    'OpenMS::String', 'String',
    'OpenMS::DataValue', 'DataValue',
    'OpenMS::ParamValue', 'ParamValue',
}

# Method name prefixes that indicate output-parameter semantics (void f(T& out)).
# Only void methods matching these prefixes will have non-const ref params
# auto-wrapped as return values. Methods like merge(T&) or update(T&) are
# input-output and must NOT be default-constructed by the wrapper.
_OUTPUT_METHOD_PREFIXES = (
    'get', 'extract', 'fill', 'compute', 'calculate', 'collect',
    'retrieve', 'fetch', 'load', 'read', 'find', 'lookup',
)

# Public field type patterns that nanobind cannot handle.
_UNBINDABLE_FIELD_PATTERNS = (
    'std::function', 'std::unique_ptr', 'std::weak_ptr',
    '(*)', 'void (*', 'void(*',
    'std::mutex', 'std::atomic', 'std::thread',
    'std::condition_variable', 'boost::',
)

# Special method implementations for critical classes
# These are C++ lambdas that implement Python methods
SPECIAL_METHODS = {
    "LightTargetedExperiment": {
        "getCompoundByRef": '''
        .def("getCompoundByRef", [](OpenSwath::LightTargetedExperiment& self, const std::string& ref) { return self.getCompoundByRef(ref); }, "ref"_a)''',
        "getPeptideByRef": '''
        .def("getPeptideByRef", [](OpenSwath::LightTargetedExperiment& self, const std::string& ref) { return self.getPeptideByRef(ref); }, "ref"_a)''',
    },
    # Compat shim: old pyOpenMS returned lists for DPosition<2>, nanobind type
    # caster returns tuples. Return lists here for backwards compatibility.
    "DBoundingBox2": {
        "minPosition": '''
        .def("minPosition", [](const OpenMS::DBoundingBox<2>& self) {
            auto p = self.minPosition();
            nb::list result;
            result.append(p[0]);
            result.append(p[1]);
            return result;
        })''',
        "maxPosition": '''
        .def("maxPosition", [](const OpenMS::DBoundingBox<2>& self) {
            auto p = self.maxPosition();
            nb::list result;
            result.append(p[0]);
            result.append(p[1]);
            return result;
        })''',
    },
    "LogMzPeak": {
        "__init__peak": '''
        .def("__init__", [](OpenMS::FLASHHelperClasses::LogMzPeak* self, const OpenMS::Peak1D& peak, bool positive) {
            new (self) OpenMS::FLASHHelperClasses::LogMzPeak(peak, positive);
        }, "peak"_a, "positive"_a)''',
    },
    "PrecalAveragine": {
        "__init__5args": '''
        .def("__init__", [](OpenMS::FLASHHelperClasses::PrecalculatedAveragine* self, double min_mass, double max_mass, double delta, OpenMS::CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine) {
            new (self) OpenMS::FLASHHelperClasses::PrecalculatedAveragine(min_mass, max_mass, delta, generator, use_RNA_averagine);
        }, "min_mass"_a, "max_mass"_a, "delta"_a, "generator"_a, "use_RNA_averagine"_a)''',
        "__init__6args": '''
        .def("__init__", [](OpenMS::FLASHHelperClasses::PrecalculatedAveragine* self, double min_mass, double max_mass, double delta, OpenMS::CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance) {
            new (self) OpenMS::FLASHHelperClasses::PrecalculatedAveragine(min_mass, max_mass, delta, generator, use_RNA_averagine, decoy_iso_distance);
        }, "min_mass"_a, "max_mass"_a, "delta"_a, "generator"_a, "use_RNA_averagine"_a, "decoy_iso_distance"_a)''',
    },
    "Software": {
        "__init__default": '''
        .def(nb::init<>(), "Default constructor")''',
    },
    "RibonucleotideDB": {
        # getInstance manually added (singleton detection not working - pxd wrap-ignore on getInstance)
        "getInstance": '''
        .def_static("getInstance", []() -> OpenMS::RibonucleotideDB* { return OpenMS::RibonucleotideDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")''',
        "getRibonucleotide": '''
        .def("getRibonucleotide", [](OpenMS::RibonucleotideDB& self, const OpenMS::String& code) -> const OpenMS::Ribonucleotide* {
            return self.getRibonucleotide(code);
        }, "code"_a, nb::rv_policy::reference, "Returns the ribonucleotide with the given code")''',
        "getRibonucleotidePrefix": '''
        .def("getRibonucleotidePrefix", [](OpenMS::RibonucleotideDB& self, const OpenMS::String& seq) -> const OpenMS::Ribonucleotide* {
            return self.getRibonucleotidePrefix(seq);
        }, "seq"_a, nb::rv_policy::reference, "Returns the ribonucleotide matching the longest prefix")''',
    },
    # NOTE: Mobilogram SPECIAL_METHODS are consolidated below (search for "Mobilogram")
    # NOTE: FileHandler SPECIAL_METHODS are consolidated below (search for "FileHandler")
    "ExperimentalDesign": {
        "fromConsensusMap": '''
        .def_static("fromConsensusMap", [](const OpenMS::ConsensusMap& c) {
            return OpenMS::ExperimentalDesign::fromConsensusMap(c);
        }, "c"_a, "Extract experimental design from consensus map")''',
        "fromFeatureMap": '''
        .def_static("fromFeatureMap", [](const OpenMS::FeatureMap& f) {
            return OpenMS::ExperimentalDesign::fromFeatureMap(f);
        }, "f"_a, "Extract experimental design from feature map")''',
        "fromIdentifications": '''
        .def_static("fromIdentifications", [](const std::vector<OpenMS::ProteinIdentification>& proteins) {
            return OpenMS::ExperimentalDesign::fromIdentifications(proteins);
        }, "proteins"_a, "Extract experimental design from identifications")''',
    },
    "IMTypes": {
        "toDriftTimeUnit": '''
        .def_static("toDriftTimeUnit", [](const OpenMS::String& dtu_string) {
            return OpenMS::toDriftTimeUnit(dtu_string);
        }, "dtu_string"_a, "Convert string to DriftTimeUnit")''',
        "driftTimeUnitToString": '''
        .def_static("driftTimeUnitToString", [](OpenMS::DriftTimeUnit value) {
            return OpenMS::driftTimeUnitToString(value);
        }, "value"_a, "Convert DriftTimeUnit to string")''',
        "toIMFormat": '''
        .def_static("toIMFormat", [](const OpenMS::String& im_format) {
            return OpenMS::toIMFormat(im_format);
        }, "im_format"_a, "Convert string to IMFormat")''',
        "imFormatToString": '''
        .def_static("imFormatToString", [](OpenMS::IMFormat value) {
            return OpenMS::imFormatToString(value);
        }, "value"_a, "Convert IMFormat to string")''',
    },
    "AcquisitionInfo": {
        "size": '''
        .def("size", [](const OpenMS::AcquisitionInfo& self) -> size_t {
            return self.size();
        }, "Returns the number of Acquisition objects")''',
        "push_back": '''
        .def("push_back", [](OpenMS::AcquisitionInfo& self, const OpenMS::Acquisition& acq) {
            self.push_back(acq);
        }, "acq"_a, "Append an Acquisition object")''',
        "resize": '''
        .def("resize", [](OpenMS::AcquisitionInfo& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the AcquisitionInfo")''',
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::AcquisitionInfo& self, size_t i) -> OpenMS::Acquisition& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, "i"_a, nb::rv_policy::reference_internal)''',
        # __len__ auto-generated from size() method
        "__setitem__": '''
        .def("__setitem__", [](OpenMS::AcquisitionInfo& self, size_t i, const OpenMS::Acquisition& acq) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = acq;
        }, "i"_a, "acq"_a)''',
    },
    "MSSpectrum": {
        "getIMData": '''
        .def("getIMData", [](const OpenMS::MSSpectrum& self) {
            auto result = self.getIMData();
            return nb::make_tuple((int)result.first, (int)result.second);
        }, "Returns (index, drift_time_unit) for ion mobility data")''',
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::MSSpectrum& self) {
            // Return (mz_array, intensity_array) as numpy arrays
            // mz as float64 (double), intensity as float32 (float) matching C++ storage
            const size_t n = self.size();
            double* mz_data = new double[n];
            float* int_data = new float[n];
            for (size_t i = 0; i < n; ++i) {
                mz_data[i] = self[i].getMZ();
                int_data[i] = self[i].getIntensity();
            }
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            auto mz_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(mz_data, {n}, mz_owner);
            auto int_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, int_owner);
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
        }, "mz"_a, "intensity"_a, "Set peaks from mz and intensity arrays")
        .def("set_peaks", [](OpenMS::MSSpectrum& self, nb::object peaks_seq) {
            // Accept a tuple or list of (mz_array, intensity_array) for compatibility with pyOpenMS API
            if (nb::len(peaks_seq) != 2) {
                throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays (mz, intensity)");
            }
            std::vector<double> mz = nb::cast<std::vector<double>>(peaks_seq[0]);
            std::vector<double> intensity = nb::cast<std::vector<double>>(peaks_seq[1]);
            const size_t n = mz.size();
            if (intensity.size() != n) {
                throw std::runtime_error("mz and intensity arrays must have same length");
            }
            self.resize(n);
            for (size_t i = 0; i < n; ++i) {
                self[i].setMZ(mz[i]);
                self[i].setIntensity(intensity[i]);
            }
        }, "peaks"_a, "Set peaks from a tuple of (mz_array, intensity_array)")''',
        "push_back": '''
        .def("push_back", [](OpenMS::MSSpectrum& self, const OpenMS::Peak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak to the spectrum")''',
        "size": '''
        .def("size", [](const OpenMS::MSSpectrum& self) {
            return self.size();
        }, "Returns the number of peaks")''',
        # __getitem__ returns copy for pyOpenMS backward compatibility
        "__getitem__": '''
        .def("__getitem__", [](const OpenMS::MSSpectrum& self, size_t i) {
            if (i >= self.size()) throw nb::index_error();
            return self[i];  // Return by value (copy)
        }, "i"_a, "Returns a copy of the peak at index i")''',
        # RangeManager methods (from RangeManagerMzInt base class)
        "getMinMZ": '''
        .def("getMinMZ", &OpenMS::MSSpectrum::getMinMZ, "Returns minimum m/z value")''',
        "getMaxMZ": '''
        .def("getMaxMZ", &OpenMS::MSSpectrum::getMaxMZ, "Returns maximum m/z value")''',
        "getMinIntensity": '''
        .def("getMinIntensity", &OpenMS::MSSpectrum::getMinIntensity, "Returns minimum intensity value")''',
        "getMaxIntensity": '''
        .def("getMaxIntensity", &OpenMS::MSSpectrum::getMaxIntensity, "Returns maximum intensity value")''',
        "clearRanges": '''
        .def("clearRanges", &OpenMS::MSSpectrum::clearRanges, "Resets all range dimensions")''',
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
        # getType, setType auto-generated
        "get_drift_time_unit": '''
        .def("get_drift_time_unit", [](const OpenMS::MSSpectrum& self) -> std::optional<OpenMS::DriftTimeUnit> {
            if (!self.containsIMData()) return std::nullopt;
            return self.getDriftTimeUnit();
        }, "Returns drift time unit if ion mobility data exists, else None")''',
        # findNearest (3 overloads) auto-generated
        "calculateTIC": '''
        .def("calculateTIC", [](const OpenMS::MSSpectrum& self) -> double {
            return static_cast<double>(self.calculateTIC());
        }, "Returns the total ion current (sum of all peak intensities)")''',
        # reserve, resize: inherited from private std::vector base, not visible to libclang
        "reserve": '''
        .def("reserve", [](OpenMS::MSSpectrum& self, size_t n) {
            self.reserve(n);
        }, "n"_a, "Reserves space for n peaks in the underlying container")''',
        "resize": '''
        .def("resize", [](OpenMS::MSSpectrum& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resizes the spectrum to contain n peaks")''',
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
        "set_peaks": '''
        .def("set_peaks", [](OpenMS::MSChromatogram& self, nb::object rt_obj, nb::object int_obj) {
            std::vector<double> rt = nb::cast<std::vector<double>>(rt_obj);
            std::vector<double> intensity = nb::cast<std::vector<double>>(int_obj);
            const size_t n = rt.size();
            if (intensity.size() != n) {
                throw std::runtime_error("rt and intensity arrays must have same length");
            }
            self.resize(n);
            for (size_t i = 0; i < n; ++i) {
                self[i].setRT(rt[i]);
                self[i].setIntensity(intensity[i]);
            }
        }, "rt"_a, "intensity"_a, "Set peaks from rt and intensity arrays")
        .def("set_peaks", [](OpenMS::MSChromatogram& self, nb::object peaks_seq) {
            if (nb::len(peaks_seq) != 2) {
                throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays (rt, intensity)");
            }
            std::vector<double> rt = nb::cast<std::vector<double>>(peaks_seq[0]);
            std::vector<double> intensity = nb::cast<std::vector<double>>(peaks_seq[1]);
            const size_t n = rt.size();
            if (intensity.size() != n) {
                throw std::runtime_error("rt and intensity arrays must have same length");
            }
            self.resize(n);
            for (size_t i = 0; i < n; ++i) {
                self[i].setRT(rt[i]);
                self[i].setIntensity(intensity[i]);
            }
        }, "peaks"_a, "Set peaks from a tuple/list of (rt_array, intensity_array)")''',
        "size": '''
        .def("size", [](const OpenMS::MSChromatogram& self) {
            return self.size();
        }, "Returns the number of peaks")''',
        "push_back": '''
        .def("push_back", [](OpenMS::MSChromatogram& self, const OpenMS::ChromatogramPeak& peak) {
            self.push_back(peak);
        }, "peak"_a, "Append a peak")''',
    },
    "Peak1D": {
        "__hash__": '''
        .def("__hash__", [](const OpenMS::Peak1D& self) {
            // Content-based hash using mz and intensity
            size_t h1 = std::hash<double>{}(self.getMZ());
            size_t h2 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1);
        })''',
    },
    "Peak2D": {
        "__hash__": '''
        .def("__hash__", [](const OpenMS::Peak2D& self) {
            // Content-based hash using mz, rt, and intensity
            size_t h1 = std::hash<double>{}(self.getMZ());
            size_t h2 = std::hash<double>{}(self.getRT());
            size_t h3 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        })''',
    },
    "ChromatogramPeak": {
        "__hash__": '''
        .def("__hash__", [](const OpenMS::ChromatogramPeak& self) {
            // Content-based hash using rt and intensity
            size_t h1 = std::hash<double>{}(self.getRT());
            size_t h2 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1);
        })''',
    },
    # Standalone enum bindings that need to be added before classes that use them
    "__enums__": {
        "DriftTimeUnit": '''
    // DriftTimeUnit enum
    nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit", nb::is_arithmetic())
        .value("NONE", OpenMS::DriftTimeUnit::NONE)
        .value("MILLISECOND", OpenMS::DriftTimeUnit::MILLISECOND)
        .value("VSSC", OpenMS::DriftTimeUnit::VSSC)
        .export_values();''',
        "IMFormat": '''
    // IMFormat enum
    nb::enum_<OpenMS::IMFormat>(m, "IMFormat", nb::is_arithmetic())
        .value("NONE", OpenMS::IMFormat::NONE)
        .value("CONCATENATED", OpenMS::IMFormat::CONCATENATED)
        .value("MULTIPLE_SPECTRA", OpenMS::IMFormat::MULTIPLE_SPECTRA)
        .value("MIXED", OpenMS::IMFormat::MIXED)
        .value("CENTROIDED", OpenMS::IMFormat::CENTROIDED)
        .value("UNKNOWN", OpenMS::IMFormat::UNKNOWN)
        .export_values();''',
        "Pi0Method": '''
    // Pi0Method enum (used by MultipleTesting)
    nb::enum_<OpenMS::Math::MultipleTesting::Pi0Method>(m, "Pi0Method", nb::is_arithmetic())
        .value("Smoother", OpenMS::Math::MultipleTesting::Pi0Method::Smoother)
        .value("Bootstrap", OpenMS::Math::MultipleTesting::Pi0Method::Bootstrap)
        .export_values();''',
        "LfdrTransform": '''
    // LfdrTransform enum (used by MultipleTesting)
    nb::enum_<OpenMS::Math::MultipleTesting::LfdrTransform>(m, "LfdrTransform", nb::is_arithmetic())
        .value("Probit", OpenMS::Math::MultipleTesting::LfdrTransform::Probit)
        .value("Logit", OpenMS::Math::MultipleTesting::LfdrTransform::Logit)
        .export_values();''',
        # LogType and FileType enums are now auto-generated from pxd files
        "Method": '''
    // RankData::Method enum
    nb::enum_<OpenMS::Math::RankData::Method>(m, "Method", nb::is_arithmetic())
        .value("Average", OpenMS::Math::RankData::Method::Average)
        .value("Min", OpenMS::Math::RankData::Method::Min)
        .value("Max", OpenMS::Math::RankData::Method::Max)
        .value("Dense", OpenMS::Math::RankData::Method::Dense)
        .value("Ordinal", OpenMS::Math::RankData::Method::Ordinal)
        .export_values();''',
        "NaNPolicy": '''
    // RankData::NaNPolicy enum
    nb::enum_<OpenMS::Math::RankData::NaNPolicy>(m, "NaNPolicy", nb::is_arithmetic())
        .value("Propagate", OpenMS::Math::RankData::NaNPolicy::Propagate)
        .value("Omit", OpenMS::Math::RankData::NaNPolicy::Omit)
        .value("Raise", OpenMS::Math::RankData::NaNPolicy::Raise)
        .export_values();''',
    },
    # Nested enum bindings that must be added AFTER the containing class is bound
    "__post_class_enums__": {
        "SpectrumSettings_SpectrumType": '''
    // Bind nested SpectrumType enum to existing SpectrumSettings class
    {
        // Get the SpectrumSettings type from the module (already bound above)
        nb::handle spectrum_settings_type = m.attr("SpectrumSettings");
        nb::enum_<OpenMS::SpectrumSettings::SpectrumType>(spectrum_settings_type, "SpectrumType", nb::is_arithmetic())
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
        # getSpectra, getChromatograms, getSpectrum, getChromatogram auto-generated
        # __getitem__ returns copy for pyOpenMS backward compatibility
        "__getitem__": '''
        .def("__getitem__", [](const OpenMS::MSExperiment& self, size_t i) {
            if (i >= self.size()) throw nb::index_error();
            return self[i];  // Return by value (copy)
        }, "i"_a, "Returns a copy of the spectrum at index i")''',
    },
    # PeptideIdentification: getHits, setHits, insertHit auto-generated
    # (mutable ref dedup prefers non-const getHits)
    # ProteinIdentification: getHits, setHits, insertHit auto-generated
    # (mutable ref dedup prefers non-const getHits)
    # Feature: getSubordinates, setSubordinates, getPeptideIdentifications,
    # setPeptideIdentifications auto-generated by regular method generator
    "FeatureMap": {
        # getProteinIdentifications, setProteinIdentifications,
        # getUnassignedPeptideIdentifications, setUnassignedPeptideIdentifications
        # auto-generated by regular method generator
        # __getitem__ auto-generated from operator[] detection (CONTAINER_CLASSES)
        "size": '''
        .def("size", [](const OpenMS::FeatureMap& self) {
            return self.size();
        }, "Returns the number of features")''',
        "push_back": '''
        .def("push_back", [](OpenMS::FeatureMap& self, const OpenMS::Feature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a feature to the map")''',
    },
    # ConsensusFeature: getPeptideIdentifications, setPeptideIdentifications,
    # getFeatureList auto-generated by regular method generator
    "ConsensusMap": {
        # getProteinIdentifications, setProteinIdentifications,
        # getUnassignedPeptideIdentifications, setUnassignedPeptideIdentifications
        # auto-generated by regular method generator
        # __getitem__ auto-generated from operator[] detection (CONTAINER_CLASSES)
        "getColumnHeaders": '''
        .def("getColumnHeaders", [](const OpenMS::ConsensusMap& self) {
            return self.getColumnHeaders();
        }, "Returns the column headers")''',
        "setColumnHeaders": '''
        .def("setColumnHeaders", [](OpenMS::ConsensusMap& self, const std::map<OpenMS::UInt64, OpenMS::ConsensusMap::ColumnHeader>& headers) {
            self.setColumnHeaders(headers);
        }, "headers"_a, "Sets the column headers")''',
        "size": '''
        .def("size", [](const OpenMS::ConsensusMap& self) {
            return self.size();
        }, "Returns the number of consensus features")''',
        "push_back": '''
        .def("push_back", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusFeature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a consensus feature to the map")''',
    },
    "AASequence": {
        "__init__": '''
        .def(nb::init<>(), "Default constructor - creates empty sequence")
        .def(nb::init<const OpenMS::AASequence&>(), "Copy constructor")
        .def("__init__", [](OpenMS::AASequence* self, const std::string& s) {
            new (self) OpenMS::AASequence(OpenMS::String(s));
        }, "sequence"_a, "Create AASequence from string (e.g., 'PEPTIDE')")''',
        "fromString": '''
        .def_static("fromString", [](const std::string& s) {
            return OpenMS::AASequence::fromString(OpenMS::String(s));
        }, "s"_a, "Create AASequence from string (deprecated - use AASequence(s) constructor)")''',
        # Default overloads using ResidueType::Full (backward compatibility with pyOpenMS)
        "getMZ": '''
        .def("getMZ", [](const OpenMS::AASequence& self, int charge) {
            return self.getMZ(charge, OpenMS::Residue::ResidueType::Full);
        }, "charge"_a, "Returns the mass-to-charge ratio of the peptide")
        .def("getMZ", [](const OpenMS::AASequence& self, int charge, OpenMS::Residue::ResidueType type) {
            return self.getMZ(charge, type);
        }, "charge"_a, "type"_a, "Returns the mass-to-charge ratio of the peptide with specified residue type")''',
        "getMonoWeight": '''
        .def("getMonoWeight", [](const OpenMS::AASequence& self) {
            return self.getMonoWeight(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the monoisotopic weight of the peptide")
        .def("getMonoWeight", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getMonoWeight(type, charge);
        }, "type"_a, "charge"_a, "Returns the monoisotopic weight of the peptide with specified residue type and charge")''',
        "getAverageWeight": '''
        .def("getAverageWeight", [](const OpenMS::AASequence& self) {
            return self.getAverageWeight(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the average weight of the peptide")
        .def("getAverageWeight", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getAverageWeight(type, charge);
        }, "type"_a, "charge"_a, "Returns the average weight of the peptide with specified residue type and charge")''',
        "fromStringPermissive": '''
        .def_static("fromStringPermissive", [](const std::string& s, bool permissive) {
            return OpenMS::AASequence::fromString(OpenMS::String(s), permissive);
        }, "s"_a, "permissive"_a, "Create AASequence from string with permissive mode")''',
        "getFormula": '''
        .def("getFormula", [](const OpenMS::AASequence& self) {
            return self.getFormula(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the empirical formula of the peptide")
        .def("getFormula", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getFormula(type, charge);
        }, "type"_a, "charge"_a, "Returns the empirical formula of the peptide with specified residue type and charge")''',
        "toBracketString": '''
        .def("toBracketString", &OpenMS::AASequence::toBracketString,
            "integer_mass"_a = true, "mass_delta"_a = false, "fixed_modifications"_a = std::vector<OpenMS::String>(),
            "Returns the bracket string representation of the peptide")''',
    },
    "Param": {
        # getValue auto-generated
        "setValue": '''
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value, const OpenMS::String& description, const std::vector<std::string>& tags) {
            self.setValue(key, value, description, tags);
        }, "key"_a, "value"_a, "description"_a, "tags"_a, "Sets a value with description and tags")
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value, const OpenMS::String& description) {
            self.setValue(key, value, description);
        }, "key"_a, "value"_a, "description"_a, "Sets a value with description")
        .def("setValue", [](OpenMS::Param& self, const OpenMS::String& key, const OpenMS::ParamValue& value) {
            self.setValue(key, value);
        }, "key"_a, "value"_a, "Set a value for a key")''',
        # exists, size, empty, clear auto-generated
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
        .def("_load_internal", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an idXML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("_store_internal", [](OpenMS::IdXMLFile& self, const OpenMS::String& filename,
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
    "MzXMLFile": {
        "load": '''
        .def("load", [](OpenMS::MzXMLFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load an mzXML file into an MSExperiment")''',
        "store": '''
        .def("store", [](OpenMS::MzXMLFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store an MSExperiment to an mzXML file")''',
    },
    "CachedmzML": {
        "load": '''
        .def_static("load", [](const OpenMS::String& filename, OpenMS::CachedmzML& cached) {
            OpenMS::CachedmzML::load(filename, cached);
        }, "filename"_a, "cached"_a, "Load a cached mzML file into a CachedmzML object (backward-compatible)")''',
    },
    "TraMLFile": {
        "load": '''
        .def("load", [](OpenMS::TraMLFile& self, const OpenMS::String& filename, OpenMS::TargetedExperiment& exp) {
            self.load(filename, exp);
        }, "filename"_a, "exp"_a, "Load a TraML file")''',
        "store": '''
        .def("store", [](OpenMS::TraMLFile& self, const OpenMS::String& filename, const OpenMS::TargetedExperiment& exp) {
            self.store(filename, exp);
        }, "filename"_a, "exp"_a, "Store to a TraML file")''',
    },
    "PeptideIdentificationList": {
        "__init__": '''
        .def(nb::init<>(), "Default constructor - creates an empty list")''',
        "size": '''
        .def("size", [](const OpenMS::PeptideIdentificationList& self) {
            return self.size();
        }, "Returns the number of peptide identifications")''',
        # __len__ auto-generated from size() method
        "__getitem__": '''
        .def("__getitem__", [](OpenMS::PeptideIdentificationList& self, size_t i) -> OpenMS::PeptideIdentification& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)''',
        # __iter__ auto-generated from begin()/end() methods
        "append": '''
        .def("append", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentification& id) {
            self.push_back(id);
        }, "id"_a, "Add a peptide identification (alias for push_back)")''',
        "extend": '''
        .def("extend", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentificationList& other) {
            for (const auto& id : other) self.push_back(id);
        }, "other"_a, "Extend with items from another list")''',
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
        .def("_load_internal", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load a pepXML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("_store_internal", [](OpenMS::PepXMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to a pepXML file")''',
    },
    "MzIdentMLFile": {
        "load": '''
        .def("_load_internal", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.load(filename, proteins, peptides);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an mzIdentML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("_store_internal", [](OpenMS::MzIdentMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an mzIdentML file")''',
    },
    "XQuestResultXMLFile": {
        "load": '''
        .def("_load_internal", [](OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename) {
            OpenMS::PeptideIdentificationList peptides;
            std::vector<OpenMS::ProteinIdentification> proteins;
            self.load(filename, peptides, proteins);
            std::vector<OpenMS::PeptideIdentification> peptide_vec(peptides.begin(), peptides.end());
            return nb::make_tuple(proteins, peptide_vec);
        }, "filename"_a, "Load an xQuest result XML file, returns tuple (proteins, peptides)")''',
        "store": '''
        .def("_store_internal", [](const OpenMS::XQuestResultXMLFile& self, const OpenMS::String& filename,
                         std::vector<OpenMS::ProteinIdentification> proteins,
                         const std::vector<OpenMS::PeptideIdentification>& peptides) {
            OpenMS::PeptideIdentificationList peptide_list(peptides);
            self.store(filename, proteins, peptide_list);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store to an xQuest result XML file")''',
    },
    "Mobilogram": {
        "size": '''
        .def("size", [](const OpenMS::Mobilogram& self) {
            return self.size();
        }, "Returns the number of peaks")''',
        "push_back": '''
        .def("push_back", [](OpenMS::Mobilogram& self, const OpenMS::MobilityPeak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak")''',
        "clear": '''
        .def("clear", [](OpenMS::Mobilogram& self) {
            self.clear();
        }, "Clear all peaks")''',
        # __len__, __getitem__, __iter__ auto-generated from vector base and begin()/end()
        "get_peaks": '''
        .def("get_peaks", [](const OpenMS::Mobilogram& self) {
            size_t n = self.size();
            double* mob_data = new double[n];
            float* int_data = new float[n];
            for (size_t i = 0; i < n; ++i) {
                mob_data[i] = self[i].getMobility();
                int_data[i] = self[i].getIntensity();
            }
            nb::capsule mob_owner(mob_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            auto mob_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(mob_data, {n}, mob_owner);
            auto int_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, int_owner);
            return nb::make_tuple(mob_arr, int_arr);
        }, "Get mobility and intensity arrays as numpy arrays")''',
        "set_peaks": '''
        .def("set_peaks", [](OpenMS::Mobilogram& self, nb::object mob_obj, nb::object int_obj) {
            auto mob_arr = nb::cast<nb::ndarray<nb::numpy, double, nb::ndim<1>>>(mob_obj);
            auto int_arr = nb::cast<nb::ndarray<nb::numpy, float, nb::ndim<1>>>(int_obj);
            size_t n = mob_arr.shape(0);
            self.resize(n);
            const double* mob_ptr = static_cast<const double*>(mob_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMobility(mob_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "mob"_a, "intensity"_a, "Set mobility and intensity from numpy arrays")
        .def("set_peaks", [](OpenMS::Mobilogram& self, nb::object peaks_seq) {
            nb::object item0, item1;
            if (nb::isinstance<nb::tuple>(peaks_seq)) {
                auto tup = nb::cast<nb::tuple>(peaks_seq);
                if (nb::len(tup) != 2) throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays");
                item0 = tup[0]; item1 = tup[1];
            } else {
                auto lst = nb::cast<nb::list>(peaks_seq);
                if (nb::len(lst) != 2) throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays");
                item0 = lst[0]; item1 = lst[1];
            }
            auto mob_arr = nb::cast<nb::ndarray<nb::numpy, double, nb::ndim<1>>>(item0);
            auto int_arr = nb::cast<nb::ndarray<nb::numpy, float, nb::ndim<1>>>(item1);
            size_t n = mob_arr.shape(0);
            self.resize(n);
            const double* mob_ptr = static_cast<const double*>(mob_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMobility(mob_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "peaks"_a, "Set peaks from [mobility_array, intensity_array]")''',
        "getFloatDataArrays": '''
        .def("getFloatDataArrays", [](const OpenMS::Mobilogram& self) {
            return self.getFloatDataArrays();
        }, "Get float data arrays")''',
        "setFloatDataArrays": '''
        .def("setFloatDataArrays", [](OpenMS::Mobilogram& self, std::vector<OpenMS::DataArrays::FloatDataArray> fda) {
            self.setFloatDataArrays(fda);
        }, "fda"_a, "Set float data arrays")''',
        "getMinMobility": '''
        .def("getMinMobility", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double min_mob = self[0].getMobility();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getMobility() < min_mob) min_mob = self[i].getMobility();
            }
            return min_mob;
        }, "Get minimum mobility value")''',
        "getMaxMobility": '''
        .def("getMaxMobility", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double max_mob = self[0].getMobility();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getMobility() > max_mob) max_mob = self[i].getMobility();
            }
            return max_mob;
        }, "Get maximum mobility value")''',
        "getMinIntensity": '''
        .def("getMinIntensity", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double min_int = self[0].getIntensity();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getIntensity() < min_int) min_int = self[i].getIntensity();
            }
            return min_int;
        }, "Get minimum intensity value")''',
        "getMaxIntensity": '''
        .def("getMaxIntensity", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double max_int = self[0].getIntensity();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getIntensity() > max_int) max_int = self[i].getIntensity();
            }
            return max_int;
        }, "Get maximum intensity value")''',
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
        "getUniqueName": '''
        .def_static("getUniqueName", [](bool include_hostname) {
            return OpenMS::File::getUniqueName(include_hostname);
        }, "include_hostname"_a = true, "Get a unique name")
        .def_static("getUniqueName", []() {
            return OpenMS::File::getUniqueName();
        }, "Get a unique name")''',
    },
    "FileHandler": {
        "loadExperiment": '''
        .def("loadExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            self.loadExperiment(filename, exp);
        }, "filename"_a, "exp"_a, "Load experiment from file")
        .def("loadExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log,
             bool rewrite_source_file, bool compute_hash) {
            self.loadExperiment(filename, exp, allowed_types, log, rewrite_source_file, compute_hash);
        }, "filename"_a, "exp"_a, "allowed_types"_a, "log"_a, "rewrite_source_file"_a, "compute_hash"_a, "Load experiment with options")''',
        "storeExperiment": '''
        .def("storeExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            self.storeExperiment(filename, exp);
        }, "filename"_a, "exp"_a, "Store experiment to file")
        .def("storeExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            self.storeExperiment(filename, exp, allowed_types, log);
        }, "filename"_a, "exp"_a, "allowed_types"_a, "log"_a, "Store experiment with options")''',
        "getType": '''
        .def_static("getType", [](const OpenMS::String& filename) {
            return OpenMS::FileHandler::getType(filename);
        }, "filename"_a, "Determine the file type from the file name")''',
    },
    "WindowMower": {
        "filterPeakSpectrumForTopNInSlidingWindow": '''
        .def("filterPeakSpectrumForTopNInSlidingWindow", [](OpenMS::WindowMower& self, OpenMS::MSSpectrum& spectrum) {
            self.filterPeakSpectrumForTopNInSlidingWindow(spectrum);
        }, "spectrum"_a, "Sliding window version (slower)")''',
    },
    # DefaultParamHandler-derived classes - these have getDefaultParameters(Param&)
    # We provide both the output-param version and a return-value version for convenience
    "EmgGradientDescent": {
        "getDefaultParameters": '''
        .def("getDefaultParameters", [](OpenMS::EmgGradientDescent& self, OpenMS::Param& params) {
            self.getDefaultParameters(params);
        }, "params"_a, "Get default parameters (fills output param)")
        .def("getDefaultParameters", [](OpenMS::EmgGradientDescent& self) {
            OpenMS::Param params;
            self.getDefaultParameters(params);
            return params;
        }, "Get default parameters (returns Param)")''',
    },
    "TargetedSpectraExtractor": {
        "getDefaultParameters": '''
        .def("getDefaultParameters", [](const OpenMS::TargetedSpectraExtractor& self, OpenMS::Param& params) {
            self.getDefaultParameters(params);
        }, "params"_a, "Get default parameters (fills output param)")
        .def("getDefaultParameters", [](const OpenMS::TargetedSpectraExtractor& self) {
            OpenMS::Param params;
            self.getDefaultParameters(params);
            return params;
        }, "Get default parameters (returns Param)")''',
    },
    "TransformationModelLinear": {
        "getDefaultParameters": '''
        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelLinear::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")''',
    },
    "Deisotoper": {
        "deisotopeAndSingleCharge": '''
        .def_static("deisotopeAndSingleCharge", [](OpenMS::MSSpectrum& spectrum,
             double fragment_tolerance, bool fragment_unit_ppm,
             int min_charge, int max_charge, bool keep_only_deisotoped,
             unsigned int min_isopeaks, unsigned int max_isopeaks,
             bool make_single_charged, bool annotate_charge,
             bool annotate_iso_peak_count, bool use_decreasing_model,
             unsigned int start_intensity_check, bool add_up_intensity,
             bool annotate_features) {
            OpenMS::Deisotoper::deisotopeAndSingleCharge(spectrum,
                fragment_tolerance, fragment_unit_ppm, min_charge, max_charge,
                keep_only_deisotoped, min_isopeaks, max_isopeaks,
                make_single_charged, annotate_charge, annotate_iso_peak_count,
                use_decreasing_model, start_intensity_check, add_up_intensity,
                annotate_features);
        }, "spectrum"_a, "fragment_tolerance"_a, "fragment_unit_ppm"_a,
           "min_charge"_a = 1, "max_charge"_a = 3, "keep_only_deisotoped"_a = false,
           "min_isopeaks"_a = 3, "max_isopeaks"_a = 10,
           "make_single_charged"_a = true, "annotate_charge"_a = false,
           "annotate_iso_peak_count"_a = false, "use_decreasing_model"_a = true,
           "start_intensity_check"_a = 2, "add_up_intensity"_a = false,
           "annotate_features"_a = false,
           "Deisotope and single charge a spectrum")''',
    },
    # NOTE: TransformationDescription SPECIAL_METHODS consolidated below
    "ModificationsDB": {
        # getInstance now auto-generated by singleton detection (Phase 6)
        "searchModifications": '''
        .def("searchModifications", [](const OpenMS::ModificationsDB& self, const OpenMS::String& mod_name, const OpenMS::String& residue, int term_spec) {
            std::set<const OpenMS::ResidueModification*> mods;
            self.searchModifications(mods, mod_name, residue, static_cast<OpenMS::ResidueModification::TermSpecificity>(term_spec));
            nb::list result;
            for (auto* m : mods) {
                result.append(nb::cast(m, nb::rv_policy::reference));
            }
            return result;
        }, "mod_name"_a, "residue"_a = "", "term_spec"_a = static_cast<int>(OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY), "Search for modifications by name")''',
        "getNumberOfModifications": '''
        .def("getNumberOfModifications", [](const OpenMS::ModificationsDB& self) {
            return self.getNumberOfModifications();
        }, "Get the number of modifications")''',
    },
    "ProteaseDB": {
        "getInstance": '''
        .def_static("getInstance", []() -> OpenMS::ProteaseDB* { return OpenMS::ProteaseDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")''',
        "getAllNames": '''
        .def("getAllNames", [](const OpenMS::ProteaseDB& self, nb::list output) {
            std::vector<OpenMS::String> all_names;
            self.getAllNames(all_names);
            for (const auto& name : all_names) {
                output.append(nb::cast(name));
            }
        }, "output"_a, "Returns all enzyme names (appends to output list)")
        .def("getAllNames", [](const OpenMS::ProteaseDB& self) {
            std::vector<OpenMS::String> all_names;
            self.getAllNames(all_names);
            return all_names;
        }, "Returns all enzyme names")''',
        "getEnzyme": '''
        .def("getEnzyme", [](const OpenMS::ProteaseDB& self, const OpenMS::String& name) -> const OpenMS::DigestionEnzymeProtein* {
            return self.getEnzyme(name);
        }, "name"_a, nb::rv_policy::reference, "Returns the enzyme with the given name (supports synonym names)")''',
        "hasEnzyme": '''
        .def("hasEnzyme", [](const OpenMS::ProteaseDB& self, const OpenMS::String& name) {
            return self.hasEnzyme(name);
        }, "name"_a, "Check if an enzyme with the given name exists")''',
    },
    "RNaseDB": {
        # getInstance now auto-generated by singleton detection (Phase 6)
        "hasEnzyme": '''
        .def("hasEnzyme", [](const OpenMS::RNaseDB& self, const OpenMS::String& name) {
            return self.hasEnzyme(name);
        }, "name"_a, "Check if an enzyme with the given name exists")''',
        "getEnzyme": '''
        .def("getEnzyme", [](const OpenMS::RNaseDB& self, const OpenMS::String& name) -> const OpenMS::DigestionEnzymeRNA* {
            return self.getEnzyme(name);
        }, "name"_a, nb::rv_policy::reference, "Get the enzyme with the given name")''',
        "getEnzymeByRegEx": '''
        .def("getEnzymeByRegEx", [](const OpenMS::RNaseDB& self, const OpenMS::String& cleavage_regex) -> const OpenMS::DigestionEnzymeRNA* {
            return self.getEnzymeByRegEx(cleavage_regex);
        }, "cleavage_regex"_a, nb::rv_policy::reference, "Get the enzyme with the given cleavage regex")''',
        "getAllNames": '''
        .def("getAllNames", [](const OpenMS::RNaseDB& self) {
            std::vector<OpenMS::String> names;
            self.getAllNames(names);
            nb::list result;
            for (const auto& n : names) result.append(nb::cast(n));
            return result;
        }, "Get all enzyme names")''',
    },
    "CrossLinksDB": {
        # getInstance now auto-generated by singleton detection (Phase 6)
        "getNumberOfModifications": '''
        .def("getNumberOfModifications", [](const OpenMS::CrossLinksDB& self) {
            return self.getNumberOfModifications();
        }, "Returns the number of modifications stored in the database")''',
        "getModification": '''
        .def("getModification", [](const OpenMS::CrossLinksDB& self, OpenMS::Size index) -> const OpenMS::ResidueModification* {
            return self.getModification(index);
        }, "index"_a, nb::rv_policy::reference, "Returns the modification with the given index")
        .def("getModification", [](const OpenMS::CrossLinksDB& self, const OpenMS::String& mod_name) -> const OpenMS::ResidueModification* {
            return self.getModification(mod_name);
        }, "mod_name"_a, nb::rv_policy::reference, "Returns the modification with the given name")''',
        "has": '''
        .def("has", [](const OpenMS::CrossLinksDB& self, const OpenMS::String& modification) {
            return self.has(modification);
        }, "modification"_a, "Returns True if the modification exists in the database")''',
        "findModificationIndex": '''
        .def("findModificationIndex", [](const OpenMS::CrossLinksDB& self, const OpenMS::String& mod_name) {
            return self.findModificationIndex(mod_name);
        }, "mod_name"_a, "Returns the index of the modification with the given name")''',
        "searchModificationsByDiffMonoMass": '''
        .def("searchModificationsByDiffMonoMass", [](const OpenMS::CrossLinksDB& self, double mass, double max_error, const OpenMS::String& residue, OpenMS::ResidueModification::TermSpecificity term_spec) {
            std::vector<OpenMS::String> mods;
            self.searchModificationsByDiffMonoMass(mods, mass, max_error, residue, term_spec);
            nb::list result;
            for (const auto& m : mods) result.append(nb::cast(m));
            return result;
        }, "mass"_a, "max_error"_a, "residue"_a, "term_spec"_a, "Collects all modifications with delta mass inside a tolerance window")''',
        "getAllSearchModifications": '''
        .def("getAllSearchModifications", [](const OpenMS::CrossLinksDB& self) {
            std::vector<OpenMS::String> mods;
            self.getAllSearchModifications(mods);
            nb::list result;
            for (const auto& m : mods) result.append(nb::cast(m));
            return result;
        }, "Collects all modifications that can be used for identification searches")''',
        "readFromOBOFile": '''
        .def("readFromOBOFile", [](OpenMS::CrossLinksDB& self, const OpenMS::String& filename) {
            self.readFromOBOFile(filename);
        }, "filename"_a, "Adds modifications from a given file in OBO format")''',
        "isInstantiated": '''
        .def("isInstantiated", [](const OpenMS::CrossLinksDB& self) {
            return self.isInstantiated();
        }, "Returns True if the database has been instantiated")''',
    },
    "ProFormaParser": {
        "parse": '''
        .def_static("parse", [](const OpenMS::String& input) {
            return OpenMS::ProFormaParser::parse(input);
        }, "input"_a, "Parse a ProForma string into a Peptidoform AST")''',
        "parseIon": '''
        .def_static("parseIon", [](const OpenMS::String& input) {
            return OpenMS::ProFormaParser::parseIon(input);
        }, "input"_a, "Parse a ProForma string into a PeptidoformIon AST")''',
        "toString": '''
        .def_static("toString", [](const OpenMS::Peptidoform& pf, OpenMS::ProFormaWriteMode mode) {
            return OpenMS::ProFormaParser::toString(pf, mode);
        }, "pf"_a, "mode"_a, "Convert a Peptidoform AST back to ProForma string notation")
        .def_static("toStringIon", [](const OpenMS::PeptidoformIon& pfi, OpenMS::ProFormaWriteMode mode) {
            return OpenMS::ProFormaParser::toString(pfi, mode);
        }, "pfi"_a, "mode"_a, "Convert a PeptidoformIon AST back to ProForma string notation")''',
        "resolveModifications": '''
        .def_static("resolveModifications", [](OpenMS::Peptidoform& pf) {
            OpenMS::ProFormaParser::resolveModifications(pf);
        }, "pf"_a, "Resolve all modifications in a Peptidoform using ModificationsDB")''',
        "toAASequence": '''
        .def_static("toAASequence", [](const OpenMS::Peptidoform& pf, OpenMS::AASequenceConversionPolicy policy) {
            return OpenMS::ProFormaParser::toAASequence(pf, policy);
        }, "pf"_a, "policy"_a, "Convert a Peptidoform to an OpenMS AASequence")''',
        "fromAASequence": '''
        .def_static("fromAASequence", [](const OpenMS::AASequence& seq) {
            return OpenMS::ProFormaParser::fromAASequence(seq);
        }, "seq"_a, "Create a Peptidoform from an OpenMS AASequence")''',
        "isRepresentableAsAASequence": '''
        .def_static("isRepresentableAsAASequence", [](const OpenMS::Peptidoform& pf) {
            return OpenMS::ProFormaParser::isRepresentableAsAASequence(pf);
        }, "pf"_a, "Check if a Peptidoform can be fully represented as an AASequence")''',
        "canCalculateMass": '''
        .def_static("canCalculateMass", [](const OpenMS::Peptidoform& pf) {
            return OpenMS::ProFormaParser::canCalculateMass(pf);
        }, "pf"_a, "Check if mass can be calculated for a Peptidoform")''',
        "getMonoWeight": '''
        .def_static("getMonoWeight", [](const OpenMS::Peptidoform& pf) {
            return OpenMS::ProFormaParser::getMonoWeight(pf);
        }, "pf"_a, "Calculate monoisotopic mass of a Peptidoform in Daltons")''',
        "getMZ": '''
        .def_static("getMZ", [](const OpenMS::PeptidoformIon& pfi) {
            return OpenMS::ProFormaParser::getMZ(pfi);
        }, "pfi"_a, "Calculate m/z for a PeptidoformIon")
        .def_static("getMZCharge", [](const OpenMS::Peptidoform& pf, int charge) {
            return OpenMS::ProFormaParser::getMZ(pf, charge);
        }, "pf"_a, "charge"_a, "Calculate m/z for a Peptidoform at a given charge state")''',
        "generateSpectrum": '''
        .def_static("generateSpectrum", [](const OpenMS::Peptidoform& pf, int min_charge, int max_charge, const std::string& ion_types, bool add_losses, bool add_metainfo) {
            return OpenMS::ProFormaParser::generateSpectrum(pf, min_charge, max_charge, ion_types, add_losses, add_metainfo);
        }, "pf"_a, "min_charge"_a, "max_charge"_a, "ion_types"_a, "add_losses"_a, "add_metainfo"_a, "Generate theoretical MS/MS spectrum for a Peptidoform")''',
        "peptidoformToJSON": '''
        .def_static("peptidoformToJSON", [](const OpenMS::Peptidoform& pf) {
            return OpenMS::toJSON(pf);
        }, "pf"_a, "Convert Peptidoform to JSON string")''',
        "peptidoformFromJSON": '''
        .def_static("peptidoformFromJSON", [](const OpenMS::String& json_str) {
            return OpenMS::peptidoformFromJSON(json_str);
        }, "json_str"_a, "Construct Peptidoform from JSON string")''',
    },
    "RankData": {
        "rankdata_double": '''
        .def_static("rankdata_double", [](std::vector<double> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_double(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)''',
        "rankdata_float": '''
        .def_static("rankdata_float", [](std::vector<float> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_float(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)''',
        "rankdata_int": '''
        .def_static("rankdata_int", [](std::vector<int> a, OpenMS::Math::RankData::Method m, OpenMS::Math::RankData::NaNPolicy p) {
            auto result = OpenMS::Math::RankData::rankdata_int(a, m, p);
            size_t n = result.size();
            double* data = new double[n];
            std::copy(result.begin(), result.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "a"_a, "method"_a, "nan_policy"_a)''',
    },
    "MultipleTesting": {
        # These methods have enum default arguments that would be evaluated at module init time
        # before the enums are registered. Provide bindings without C++ default values.
        "pi0Est": '''
        .def_static("pi0Est", [](std::vector<double> p_values,
                                  std::vector<double> lambda,
                                  OpenMS::Math::MultipleTesting::Pi0Method method,
                                  int smooth_df,
                                  bool smooth_log_pi0) {
            return OpenMS::Math::MultipleTesting::pi0Est(p_values, lambda, method, smooth_df, smooth_log_pi0);
        }, "p_values"_a, "lambda_"_a, "method"_a, "smooth_df"_a = 3, "smooth_log_pi0"_a = false,
           "Estimate pi0 (proportion of true null hypotheses)")''',
        "lfdr": '''
        .def_static("lfdr", [](std::vector<double> p_values,
                               double pi0,
                               bool trunc,
                               bool monotone,
                               OpenMS::Math::MultipleTesting::LfdrTransform transf,
                               double adj,
                               double eps,
                               size_t gridsize,
                               double cut) {
            return OpenMS::Math::MultipleTesting::lfdr(p_values, pi0, trunc, monotone, transf, adj, eps, gridsize, cut);
        }, "p_values"_a, "pi0"_a, "trunc"_a = true, "monotone"_a = true, "transf"_a, "adj"_a = 1.5, "eps"_a = 1e-8, "gridsize"_a = 100, "cut"_a = 0.05,
           "Compute local FDR values")''',
        # Enum conversion methods
        "pi0MethodToString": '''
        .def_static("pi0MethodToString", &OpenMS::Math::MultipleTesting::pi0MethodToString,
           "m"_a, "Convert Pi0Method enum to string")''',
        "toPi0Method": '''
        .def_static("toPi0Method", &OpenMS::Math::MultipleTesting::toPi0Method,
           "s"_a, "Convert string to Pi0Method enum (case insensitive)")''',
        "lfdrTransformToString": '''
        .def_static("lfdrTransformToString", &OpenMS::Math::MultipleTesting::lfdrTransformToString,
           "t"_a, "Convert LfdrTransform enum to string")''',
        "toLfdrTransform": '''
        .def_static("toLfdrTransform", &OpenMS::Math::MultipleTesting::toLfdrTransform,
           "s"_a, "Convert string to LfdrTransform enum")''',
    },
    "MSNumpressCoder": {
        "encodeNP": '''
        .def("encodeNP", [](OpenMS::MSNumpressCoder& self, std::vector<double> in, nb::object out_obj, bool zlib_compression, OpenMS::MSNumpressCoder::NumpressConfig config) {
            OpenMS::String result;
            self.encodeNP(in, result, zlib_compression, config);
            // Write result back to the output String-like object
            if (nb::hasattr(out_obj, "_value")) {
                out_obj.attr("_value") = std::string(result);
            }
        }, "in"_a, "result"_a, "zlib_compression"_a, "config"_a, "Encode vector of doubles to Base64 numpress string")''',
        "decodeNP": '''
        .def("decodeNP", [](OpenMS::MSNumpressCoder& self, nb::object in_obj, nb::list out, bool zlib_compression, OpenMS::MSNumpressCoder::NumpressConfig config) {
            std::string in_str;
            if (nb::isinstance<nb::bytes>(in_obj)) {
                auto b = nb::cast<nb::bytes>(in_obj);
                in_str = std::string(b.c_str(), b.size());
            } else {
                in_str = nb::cast<std::string>(in_obj);
            }
            std::vector<double> result;
            self.decodeNP(in_str, result, zlib_compression, config);
            for (double v : result) out.append(v);
        }, "in"_a, "out"_a, "zlib_compression"_a, "config"_a, "Decode Base64 numpress string to vector of doubles")''',
        "encodeNPRaw": '''
        .def("encodeNPRaw", [](OpenMS::MSNumpressCoder& self, std::vector<double> in, nb::object out_obj, OpenMS::MSNumpressCoder::NumpressConfig config) {
            OpenMS::String result;
            self.encodeNPRaw(in, result, config);
            if (nb::hasattr(out_obj, "_value")) {
                // Raw encoding may contain null bytes - use the full size
                out_obj.attr("_value") = nb::bytes(result.c_str(), result.size());
            }
        }, "in"_a, "result"_a, "config"_a, "Encode vector of doubles to raw numpress byte array")''',
        "decodeNPRaw": '''
        .def("decodeNPRaw", [](OpenMS::MSNumpressCoder& self, nb::object in_obj, nb::list out, OpenMS::MSNumpressCoder::NumpressConfig config) {
            std::string in_str;
            if (nb::isinstance<nb::bytes>(in_obj)) {
                auto b = nb::cast<nb::bytes>(in_obj);
                in_str = std::string(b.c_str(), b.size());
            } else {
                in_str = nb::cast<std::string>(in_obj);
            }
            std::vector<double> result;
            self.decodeNPRaw(in_str, result, config);
            for (double v : result) out.append(v);
        }, "in"_a, "out"_a, "config"_a, "Decode raw numpress byte array to vector of doubles")''',
    },
    "IsobaricQuantitationMethod": {
        "getIsotopeCorrectionMatrix": '''
        .def("getIsotopeCorrectionMatrix", [](const OpenMS::IsobaricQuantitationMethod& self) {
            return self.getIsotopeCorrectionMatrix();
        }, "Get the isotope correction matrix")''',
    },
    "AbsoluteQuantitation": {
        "setQuantMethods": '''
        .def("setQuantMethods", [](OpenMS::AbsoluteQuantitation& self, const std::vector<OpenMS::AbsoluteQuantitationMethod>& quant_methods) {
            std::vector<OpenMS::AbsoluteQuantitationMethod> methods_copy(quant_methods);
            self.setQuantMethods(methods_copy);
        }, "quant_methods"_a, "Set the quantitation methods")''',
        "getQuantMethods": '''
        .def("getQuantMethods", [](OpenMS::AbsoluteQuantitation& self) {
            return self.getQuantMethods();
        }, "Get the quantitation methods")''',
    },
    # Peptide: public fields now auto-generated from libclang (Phase 5)
    "ElementDB": {
        # getInstance now auto-generated by singleton detection (Phase 6)
        "getElement_name": '''
        .def("getElement", [](const OpenMS::ElementDB& self, const std::string& name) -> const OpenMS::Element* {
            return self.getElement(name);
        }, "name"_a, nb::rv_policy::reference, "Get element by name or symbol")''',
        "getElement_number": '''
        .def("getElement", [](const OpenMS::ElementDB& self, unsigned int atomic_number) -> const OpenMS::Element* {
            return self.getElement(atomic_number);
        }, "atomic_number"_a, nb::rv_policy::reference, "Get element by atomic number")''',
        "hasElement_name": '''
        .def("hasElement", [](const OpenMS::ElementDB& self, const std::string& name) {
            return self.hasElement(name);
        }, "name"_a, "Check if element exists by name or symbol")''',
        "hasElement_number": '''
        .def("hasElement", [](const OpenMS::ElementDB& self, unsigned int atomic_number) {
            return self.hasElement(atomic_number);
        }, "atomic_number"_a, "Check if element exists by atomic number")''',
        "addElement": '''
        .def("addElement", [](OpenMS::ElementDB& self, const std::string& name, const std::string& symbol,
                unsigned int an, const std::map<unsigned int, double>& abundance,
                const std::map<unsigned int, double>& mass, bool replace_existing) {
            self.addElement(name, symbol, an, abundance, mass, replace_existing);
        }, "name"_a, "symbol"_a, "atomic_number"_a, "abundance"_a, "mass"_a, "replace_existing"_a,
        "Add a new element to the database")''',
    },
    "FASTAFile": {
        "load": '''
        .def("load", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, nb::list& output) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> entries;
            self.load(filename, entries);
            for (const auto& e : entries) {
                output.append(nb::cast(e));
            }
        }, "filename"_a, "data"_a, "Load a FASTA file. Appends FASTAEntry objects to output list")
        .def("load", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> entries;
            self.load(filename, entries);
            return entries;
        }, "filename"_a, "Load a FASTA file. Returns list of FASTAEntry objects")''',
        "store": '''
        .def("store", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, const std::vector<OpenMS::FASTAFile::FASTAEntry>& entries) {
            self.store(filename, entries);
        }, "filename"_a, "entries"_a, "Store a FASTA file. Takes list of FASTAEntry objects")
        .def("store", [](const OpenMS::FASTAFile& self, const OpenMS::String& filename, const nb::list& entries) {
            std::vector<OpenMS::FASTAFile::FASTAEntry> fasta_entries;
            for (auto item : entries) {
                // Try to cast as FASTAEntry first, fall back to tuple
                try {
                    auto entry = nb::cast<OpenMS::FASTAFile::FASTAEntry>(item);
                    fasta_entries.push_back(entry);
                } catch (const nb::cast_error&) {
                    auto tup = nb::cast<nb::tuple>(item);
                    OpenMS::FASTAFile::FASTAEntry e;
                    e.identifier = nb::cast<std::string>(tup[0]);
                    e.description = nb::cast<std::string>(tup[1]);
                    e.sequence = nb::cast<std::string>(tup[2]);
                    fasta_entries.push_back(e);
                }
            }
            self.store(filename, fasta_entries);
        }, "filename"_a, "entries"_a, "Store a FASTA file. Takes list of FASTAEntry objects or (identifier, description, sequence) tuples")''',
    },
    "TransformationDescription": {
        "getDataPoints": '''
        .def("getDataPoints", [](const OpenMS::TransformationDescription& self) {
            return self.getDataPoints();
        }, "Get the data points used for the transformation")''',
        "setDataPoints": '''
        .def("setDataPoints", [](OpenMS::TransformationDescription& self, const std::vector<std::pair<double, double>>& data) {
            OpenMS::TransformationDescription::DataPoints dp;
            for (const auto& p : data) {
                dp.push_back(std::make_pair(p.first, p.second));
            }
            self.setDataPoints(dp);
        }, "data"_a, "Set the data points for the transformation")''',
        "apply": '''
        .def("apply", [](const OpenMS::TransformationDescription& self, double value) {
            return self.apply(value);
        }, "value"_a, "Apply the transformation to a value")''',
        "getModelTypes": '''
        .def_static("getModelTypes", [](nb::list result) {
            std::vector<OpenMS::String> types;
            OpenMS::TransformationDescription::getModelTypes(types);
            for (const auto& t : types) {
                result.append(nb::cast(t));
            }
        }, "result"_a, "Get available model types (fills list)")''',
    },
    # Pi0Result: public fields now auto-generated from libclang (Phase 5)
    # MultipleTesting: Moved to earlier in the file to avoid duplicate key overwrite
    "ProteaseDigestion": {
        "digest": '''
        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein, nb::list output) {
            // Backward-compatible: appends results to the provided Python list
            std::vector<OpenMS::AASequence> result;
            self.digest(protein, result);
            for (const auto& seq : result) {
                output.append(nb::cast(seq));
            }
        }, "protein"_a, "output"_a, "Performs the enzymatic digestion of a protein sequence (appends to output list)")
        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein) {
            std::vector<OpenMS::AASequence> output;
            self.digest(protein, output);
            return output;
        }, "protein"_a, "Performs the enzymatic digestion of a protein sequence, returns list of peptides")''',
        "isValidProduct": '''
        .def("isValidProduct", nb::overload_cast<const OpenMS::AASequence&, int, int, bool, bool, bool>(&OpenMS::ProteaseDigestion::isValidProduct, nb::const_),
            "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a = true, "allow_nterm_protein_cleavage"_a = false, "allow_random_asp_pro_cleavage"_a = false,
            "Check if peptide is a valid digestion product of protein (AASequence version)")
        .def("isValidProduct", nb::overload_cast<const OpenMS::String&, int, int, bool, bool, bool>(&OpenMS::ProteaseDigestion::isValidProduct, nb::const_),
            "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a = true, "allow_nterm_protein_cleavage"_a = false, "allow_random_asp_pro_cleavage"_a = false,
            "Check if peptide is a valid digestion product of protein (String version)")''',
    },
    "CoarseIsotopePatternGenerator": {
        "__init__default": '''
        .def(nb::init<>())
        .def(nb::init<OpenMS::Size>(), "max_isotope"_a)''',
    },
    "MRMFeatureQC": {
        "component_qcs_prop": '''
        .def_prop_rw("component_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentQCs>& { return self.component_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentQCs> v) { self.component_qcs = std::move(v); })''',
        "component_group_qcs_prop": '''
        .def_prop_rw("component_group_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentGroupQCs>& { return self.component_group_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentGroupQCs> v) { self.component_group_qcs = std::move(v); })''',
        "component_group_pair_qcs_prop": '''
        .def_prop_rw("component_group_pair_qcs",
            [](OpenMS::MRMFeatureQC& self) -> std::vector<OpenMS::MRMFeatureQC::ComponentGroupPairQCs>& { return self.component_group_pair_qcs; },
            [](OpenMS::MRMFeatureQC& self, std::vector<OpenMS::MRMFeatureQC::ComponentGroupPairQCs> v) { self.component_group_pair_qcs = std::move(v); })''',
    },
    # ConvexHull2D.getBoundingBox: auto-generated (template incomplete type fix)
    "PeakIntegrator": {
        "integratePeak_chrom": '''
        .def("integratePeak", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right) {
            return self.integratePeak(chrom, left, right);
        }, "chromatogram"_a, "left"_a, "right"_a, "Integrate peak in chromatogram")''',
        "integratePeak_spec": '''
        .def("integratePeak", [](OpenMS::PeakIntegrator& self, const OpenMS::MSSpectrum& spec, double left, double right) {
            return self.integratePeak(spec, left, right);
        }, "spectrum"_a, "left"_a, "right"_a, "Integrate peak in spectrum")''',
        "estimateBackground_chrom": '''
        .def("estimateBackground", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right, double apex) {
            return self.estimateBackground(chrom, left, right, apex);
        }, "chromatogram"_a, "left"_a, "right"_a, "peak_apex_pos"_a, "Estimate background in chromatogram")''',
        "calculatePeakShapeMetrics_chrom": '''
        .def("calculatePeakShapeMetrics", [](OpenMS::PeakIntegrator& self, const OpenMS::MSChromatogram& chrom, double left, double right, double height, double apex) {
            return self.calculatePeakShapeMetrics(chrom, left, right, height, apex);
        }, "chromatogram"_a, "left"_a, "right"_a, "peak_height"_a, "peak_apex_pos"_a, "Calculate peak shape metrics")''',
    },
    "PeakGroup": {
        "getMonoMass": '''
        .def("getMonoMass", &OpenMS::PeakGroup::getMonoMass, "Returns the monoisotopic mass")''',
    },
    "TheoreticalSpectrumGenerator": {
        "getSpectrum": '''
        .def("getSpectrum", [](OpenMS::TheoreticalSpectrumGenerator& self, OpenMS::MSSpectrum& spec, const OpenMS::AASequence& peptide, int min_charge, int max_charge) {
            self.getSpectrum(spec, peptide, min_charge, max_charge);
        }, "spec"_a, "peptide"_a, "min_charge"_a, "max_charge"_a, "Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters")''',
    },
    # DeconvolvedSpectrum: __getitem__ auto-generated from operator[] detection
    "Spectrum": {
        "setMZArray": '''
        .def("setMZArray", [](OpenMS::Interfaces::Spectrum& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setMZArray(arr);
        }, "data"_a, "Set m/z array from list")''',
        "setIntensityArray": '''
        .def("setIntensityArray", [](OpenMS::Interfaces::Spectrum& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")''',
        "getMZArray": '''
        .def("getMZArray", [](const OpenMS::Interfaces::Spectrum& self) {
            auto arr = self.getMZArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get m/z array")''',
        "getIntensityArray": '''
        .def("getIntensityArray", [](const OpenMS::Interfaces::Spectrum& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array")''',
    },
    "ChromatogramExtractor": {
        "extractChromatograms": '''
        .def("extractChromatograms", [](OpenMS::ChromatogramExtractor& self,
                std::shared_ptr<OpenMS::SpectrumAccessOpenMS> input,
                nb::list output_py,
                std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates,
                double mz_extraction_window,
                bool ppm,
                double im_extraction_window,
                const std::string& filter) {
            // Convert Python list to C++ vector
            std::vector<std::shared_ptr<OpenSwath::OSChromatogram>> output;
            for (size_t i = 0; i < nb::len(output_py); ++i) {
                output.push_back(nb::cast<std::shared_ptr<OpenSwath::OSChromatogram>>(output_py[i]));
            }
            self.extractChromatograms(input, output, extraction_coordinates, mz_extraction_window, ppm, im_extraction_window, filter);
            // Update the Python list in-place
            while (nb::len(output_py) > 0) { output_py.attr("pop")(); }
            for (auto& c : output) { output_py.append(nb::cast(c)); }
        }, "input"_a, "output"_a, "extraction_coordinates"_a, "mz_extraction_window"_a, "ppm"_a, "im_extraction_window"_a, "filter"_a)''',
        "prepare_coordinates": '''
        .def_static("prepare_coordinates", [](
                nb::list output_chromatograms_py,
                nb::list extraction_coordinates_py,
                OpenMS::TargetedExperiment& targeted,
                double rt_extraction_window,
                bool ms1,
                int ms1_isotopes) {
            std::vector<std::shared_ptr<OpenSwath::OSChromatogram>> output_chromatograms;
            std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates;
            OpenMS::ChromatogramExtractor::prepare_coordinates(output_chromatograms, extraction_coordinates, targeted, rt_extraction_window, ms1, ms1_isotopes);
            // Update Python lists in-place
            while (nb::len(output_chromatograms_py) > 0) { output_chromatograms_py.attr("pop")(); }
            for (auto& c : output_chromatograms) { output_chromatograms_py.append(nb::cast(c)); }
            while (nb::len(extraction_coordinates_py) > 0) { extraction_coordinates_py.attr("pop")(); }
            for (auto& c : extraction_coordinates) { extraction_coordinates_py.append(nb::cast(c)); }
        }, "output_chromatograms"_a, "extraction_coordinates"_a, "targeted"_a, "rt_extraction_window"_a, "ms1"_a = false, "ms1_isotopes"_a = 0, "Prepare extraction coordinates from targeted experiment")''',
    },
    "ChromatogramExtractorAlgorithm": {
        "extractChromatograms": '''
        .def("extractChromatograms", [](OpenMS::ChromatogramExtractorAlgorithm& self,
                std::shared_ptr<OpenMS::SpectrumAccessOpenMS> input,
                std::vector<std::shared_ptr<OpenSwath::OSChromatogram>>& output,
                std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> extraction_coordinates,
                double mz_extraction_window,
                bool ppm,
                double im_extraction_window,
                const std::string& filter) {
            self.extractChromatograms(input, output, extraction_coordinates, mz_extraction_window, ppm, im_extraction_window, filter);
        }, "input"_a, "output"_a, "extraction_coordinates"_a, "mz_extraction_window"_a, "ppm"_a, "im_extraction_window"_a, "filter"_a)''',
    },
    "Chromatogram": {
        "setTimeArray": '''
        .def("setTimeArray", [](OpenMS::Interfaces::Chromatogram& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setTimeArray(arr);
        }, "data"_a, "Set time array from list")''',
        "setIntensityArray": '''
        .def("setIntensityArray", [](OpenMS::Interfaces::Chromatogram& self, std::vector<double> data) {
            auto arr = std::make_shared<OpenMS::Interfaces::BinaryDataArray>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")''',
        "getTimeArray": '''
        .def("getTimeArray", [](const OpenMS::Interfaces::Chromatogram& self) {
            auto arr = self.getTimeArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get time array")''',
        "getIntensityArray": '''
        .def("getIntensityArray", [](const OpenMS::Interfaces::Chromatogram& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array")''',
    },
    # SequenceCoverage: static method taking vector<AASequence>
    "SequenceCoverage": {
        "getCoverage": '''
        .def_static("getCoverage", [](const OpenMS::AASequence& protein, const std::vector<OpenMS::AASequence>& peptides) {
            return OpenMS::SequenceCoverage::getCoverage(protein, peptides);
        }, "protein"_a, "peptides"_a, "Compute sequence coverage percentage (0-100)")''',
    },
    # IDScoreSwitcherAlgorithm: template method for findScoreType and static methods
    "IDScoreSwitcherAlgorithm": {
        "findScoreType": '''
        .def("findScoreType", [](OpenMS::IDScoreSwitcherAlgorithm& self, OpenMS::PeptideIdentification& id, OpenMS::Scores::IDType score_type) {
            return self.findScoreType(id, score_type);
        }, "id"_a, "score_type"_a, "Search for a score type in a PeptideIdentification")''',
        "switchToScoreType": '''
        .def_static("switchToScoreType", nb::overload_cast<OpenMS::PeptideIdentificationList&, OpenMS::String>(&OpenMS::IDScoreSwitcherAlgorithm::switchToScoreType),
            "pep_ids"_a, "requested_score_type_as_string"_a, "Switch the score type of peptide identifications")
        .def_static("switchToScoreType", [](OpenMS::ConsensusMap& cmap, OpenMS::String requested_score_type, bool include_unassigned) {
            return OpenMS::IDScoreSwitcherAlgorithm::switchToScoreType(cmap, requested_score_type, include_unassigned);
        }, "cmap"_a, "requested_score_type_as_string"_a, "include_unassigned"_a = true, "Switch the score type of a ConsensusMap")''',
        "switchBackScoreType": '''
        .def_static("switchBackScoreType", nb::overload_cast<OpenMS::PeptideIdentificationList&, OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult>(&OpenMS::IDScoreSwitcherAlgorithm::switchBackScoreType),
            "pep_ids"_a, "isr"_a, "Revert scores to original type")
        .def_static("switchBackScoreType", [](OpenMS::ConsensusMap& cmap, OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult isr, bool include_unassigned) {
            OpenMS::IDScoreSwitcherAlgorithm::switchBackScoreType(cmap, isr, include_unassigned);
        }, "cmap"_a, "isr"_a, "include_unassigned"_a = true, "Revert scores of a ConsensusMap to original type")''',
    },
    # Scores: static utility class for score type handling
    "Scores": {
        "isScoreType": '''
        .def_static("isScoreType", &OpenMS::Scores::isScoreType, "score_name"_a, "type"_a,
            "Check if the given score name corresponds to a specific ID score type")''',
        "parseIDType": '''
        .def_static("parseIDType", &OpenMS::Scores::parseIDType, "score_type"_a,
            "Convert a string representation of an ID score type to an IDType enum")''',
        "isHigherBetter": '''
        .def_static("isHigherBetter", &OpenMS::Scores::isHigherBetter, "type"_a,
            "Determine whether a higher score is better for the given ID score type")''',
        "getAllIDScoreNames": '''
        .def_static("getAllIDScoreNames", &OpenMS::Scores::getAllIDScoreNames,
            "Get a vector of all ID score names used in OpenMS")''',
        "normalizeScoreName": '''
        .def_static("normalizeScoreName", &OpenMS::Scores::normalizeScoreName, "score_name"_a,
            "Normalize a score name by removing the '_score' suffix if present")''',
        "isKnownScoreType": '''
        .def_static("isKnownScoreType", &OpenMS::Scores::isKnownScoreType, "score_name"_a,
            "Check if a score name is a known score type after normalization")''',
    },
    # AASequence: fromStringPermissive moved to main entry to avoid duplicate key issue
    # EmpiricalFormula: add getElementalComposition (alias for toMap)
    "EmpiricalFormula": {
        "getElementalComposition": '''
        .def("getElementalComposition", [](const OpenMS::EmpiricalFormula& self) {
            return self.toMap();
        }, "Get elemental composition as a dict {'Symbol': NrAtoms}")''',
        "getIsotopeDistribution": '''
        .def("getIsotopeDistribution", [](const OpenMS::EmpiricalFormula& self, const OpenMS::CoarseIsotopePatternGenerator& method) {
            return self.getIsotopeDistribution(method);
        }, "method"_a, "Returns isotope distribution using CoarseIsotopePatternGenerator")
        .def("getIsotopeDistribution", [](const OpenMS::EmpiricalFormula& self, const OpenMS::FineIsotopePatternGenerator& method) {
            return self.getIsotopeDistribution(method);
        }, "method"_a, "Returns isotope distribution using FineIsotopePatternGenerator")''',
        "__str__": '''
        .def("__str__", [](const OpenMS::EmpiricalFormula& self) {
            return std::string(self.toString());
        }, "Returns the formula as a string")''',
    },
}

# Nested class bindings that are emitted in the same module as a parent class.
# Maps parent_class_name -> list of (binding_code, extra_includes) tuples.
# These are standalone nb::class_<> definitions for C++ nested types
# (e.g. OpenMS::Param::ParamEntry) that can't be auto-generated.
NESTED_CLASS_BINDINGS = {
    "Param": [
        ('''
    nb::class_<OpenMS::Param::ParamEntry>(m, "ParamEntry")
        .def(nb::init<>())
        .def(nb::init<const std::string&, const OpenMS::ParamValue&, const std::string&>(), "name"_a, "value"_a, "description"_a)
        .def_rw("name", &OpenMS::Param::ParamEntry::name)
        .def_rw("description", &OpenMS::Param::ParamEntry::description)
        .def_rw("value", &OpenMS::Param::ParamEntry::value)
        .def_rw("valid_strings", &OpenMS::Param::ParamEntry::valid_strings)
        .def_rw("max_float", &OpenMS::Param::ParamEntry::max_float)
        .def_rw("min_float", &OpenMS::Param::ParamEntry::min_float)
        .def_rw("max_int", &OpenMS::Param::ParamEntry::max_int)
        .def_rw("min_int", &OpenMS::Param::ParamEntry::min_int)
        .def("isValid", [](const OpenMS::Param::ParamEntry& self) {
            std::string msg;
            bool valid = self.isValid(msg);
            return nb::make_tuple(valid, msg);
        }, "Check if value fulfills restrictions. Returns (valid, message)")
        .def("__eq__", &OpenMS::Param::ParamEntry::operator==)
        ;''',
         ["<OpenMS/DATASTRUCTURES/Param.h>"]),
        ('''
    nb::class_<OpenMS::Param::ParamNode>(m, "ParamNode")
        .def(nb::init<>())
        .def(nb::init<const std::string&, const std::string&>(), "name"_a, "description"_a)
        .def_rw("name", &OpenMS::Param::ParamNode::name)
        .def_rw("description", &OpenMS::Param::ParamNode::description)
        .def_rw("entries", &OpenMS::Param::ParamNode::entries)
        .def_rw("nodes", &OpenMS::Param::ParamNode::nodes)
        .def("size", &OpenMS::Param::ParamNode::size)
        .def("suffix", &OpenMS::Param::ParamNode::suffix, "key"_a)
        .def("__eq__", &OpenMS::Param::ParamNode::operator==)
        ;''',
         ["<OpenMS/DATASTRUCTURES/Param.h>"]),
    ],
    "PeakFileOptions": [
        ('''
    nb::class_<OpenMS::DRange<1>>(m, "DRange1")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DRange<1>&>())
        .def("__init__", [](OpenMS::DRange<1>* self, double min_val, double max_val) {
            new (self) OpenMS::DRange<1>(OpenMS::DPosition<1>(min_val), OpenMS::DPosition<1>(max_val));
        }, "min"_a, "max"_a)
        .def("__eq__", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self == other;
        })
        .def("encloses", [](const OpenMS::DRange<1>& self, double pos) {
            return self.encloses(OpenMS::DPosition<1>(pos));
        }, "position"_a, "Check if position is within this range")
        .def("united", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self.united(other);
        }, "other"_a, "Returns the union of this range with another")
        .def("isIntersected", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self.isIntersected(other);
        }, "other"_a, "Check if two ranges intersect")
        .def("isEmpty", [](const OpenMS::DRange<1>& self) { return self.isEmpty(); })
        .def("minX", [](const OpenMS::DRange<1>& self) { return self.minPosition()[0]; })
        .def("maxX", [](const OpenMS::DRange<1>& self) { return self.maxPosition()[0]; })
        .def("setMinX", [](OpenMS::DRange<1>& self, double c) {
            self.min_[0] = c;
        }, "c"_a)
        .def("setMaxX", [](OpenMS::DRange<1>& self, double c) {
            self.max_[0] = c;
        }, "c"_a)
        .def("__repr__", [](const OpenMS::DRange<1>& self) {
            std::ostringstream oss;
            oss << "DRange1(" << self.minPosition()[0] << ", " << self.maxPosition()[0] << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::DRange<1>& self) {
            std::ostringstream oss;
            oss << "DRange1(" << self.minPosition()[0] << ", " << self.maxPosition()[0] << ")";
            return oss.str();
        })
        ;''',
         ["<OpenMS/DATASTRUCTURES/DRange.h>", "<OpenMS/DATASTRUCTURES/DPosition.h>"]),
    ],
    "MSNumpressCoder": [
        ('''
    nb::class_<OpenMS::MSNumpressCoder::NumpressConfig>(m, "NumpressConfig")
        .def(nb::init<>())
        .def_rw("numpressFixedPoint", &OpenMS::MSNumpressCoder::NumpressConfig::numpressFixedPoint)
        .def_rw("numpressErrorTolerance", &OpenMS::MSNumpressCoder::NumpressConfig::numpressErrorTolerance)
        .def_rw("np_compression", &OpenMS::MSNumpressCoder::NumpressConfig::np_compression)
        .def_rw("estimate_fixed_point", &OpenMS::MSNumpressCoder::NumpressConfig::estimate_fixed_point)
        .def_rw("linear_fp_mass_acc", &OpenMS::MSNumpressCoder::NumpressConfig::linear_fp_mass_acc)
        ;''',
         ["<OpenMS/FORMAT/MSNumpressCoder.h>"]),
    ],
    "IDScoreSwitcherAlgorithm": [
        ('''
    nb::class_<OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult>(m, "ScoreSearchResult")
        .def(nb::init<>())
        .def_rw("is_main_score_type", &OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult::is_main_score_type,
                "True if the main score is already of the requested score type")
        .def_rw("score_name", &OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult::score_name,
                "Name of score to use (main score name or meta value name)")
        .def("__repr__", [](const OpenMS::IDScoreSwitcherAlgorithm::ScoreSearchResult& self) {
            std::ostringstream oss;
            oss << "ScoreSearchResult(is_main_score_type=" << (self.is_main_score_type ? "True" : "False")
                << ", score_name='" << self.score_name << "')";
            return oss.str();
        })
        ;''',
         ["<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>"]),
        ('''
    nb::class_<OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult>(m, "IDSwitchResult")
        .def(nb::init<>())
        .def_rw("original_score_name", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_name,
                "The name of the original score used before the switch")
        .def_rw("original_score_higher_better", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_higher_better,
                "Whether a higher original score is better")
        .def_rw("original_score_type", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::original_score_type,
                "The type of the original score")
        .def_rw("requested_score_higher_better", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_higher_better,
                "Whether a higher requested score is better")
        .def_rw("requested_score_type", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_type,
                "The type of the requested score")
        .def_rw("requested_score_name", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::requested_score_name,
                "The search engine score name (e.g. 'X!Tandem_score') or score category (e.g. 'PEP')")
        .def_rw("score_switched", &OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult::score_switched,
                "Flag indicating whether the main score was switched")
        .def("__repr__", [](const OpenMS::IDScoreSwitcherAlgorithm::IDSwitchResult& self) {
            std::ostringstream oss;
            oss << "IDSwitchResult(original_score_name='" << self.original_score_name
                << "', original_score_higher_better=" << (self.original_score_higher_better ? "True" : "False")
                << ", score_switched=" << (self.score_switched ? "True" : "False")
                << ")";
            return oss.str();
        })
        ;''',
         ["<OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>"]),
    ],
    "FASTAFile": [
        ('''
    nb::class_<OpenMS::FASTAFile::FASTAEntry>(m, "FASTAEntry", "Represents a single FASTA entry with identifier, description, and sequence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FASTAFile::FASTAEntry&>())
        .def_rw("identifier", &OpenMS::FASTAFile::FASTAEntry::identifier, "Protein/sequence identifier")
        .def_rw("description", &OpenMS::FASTAFile::FASTAEntry::description, "Description from FASTA header")
        .def_rw("sequence", &OpenMS::FASTAFile::FASTAEntry::sequence, "The amino acid or nucleotide sequence")
        .def("headerMatches", &OpenMS::FASTAFile::FASTAEntry::headerMatches, "rhs"_a, "Returns True if identifier and description match the given entry")
        .def("sequenceMatches", &OpenMS::FASTAFile::FASTAEntry::sequenceMatches, "rhs"_a, "Returns True if the sequence matches the given entry")
        .def("__eq__", &OpenMS::FASTAFile::FASTAEntry::operator==)
        .def("__repr__", [](const OpenMS::FASTAFile::FASTAEntry& self) {
            std::ostringstream oss;
            oss << "FASTAEntry(identifier='" << self.identifier << "', description='" << self.description.substr(0, 50) << (self.description.size() > 50 ? "..." : "") << "')";
            return oss.str();
        })
        ;''',
         ["<OpenMS/FORMAT/FASTAFile.h>"]),
    ],
}

# Fully handwritten class bindings for classes that can't go through the normal
# generation pipeline (e.g., template specializations, complex nested types).
# These are assigned to modules by hash of class name.
HANDWRITTEN_CLASSES = {
    "OSBinaryDataArray": {
        "binding": '''
    using OSBDA = OpenSwath::OSBinaryDataArray;
    nb::class_<OSBDA>(m, "OSBinaryDataArray")
        .def(nb::init<>())
        .def(nb::init<const OSBDA&>())
        .def_rw("data", &OSBDA::data)
        .def_rw("description", &OSBDA::description)
        .def("get_data", [](const OSBDA& self) {
            return self.data;
        }, "Access to a copy of the underlying data")
        .def("get_data_mv", [](nb::object self_obj) -> nb::object {
            auto& self = nb::cast<OSBDA&>(self_obj);
            size_t shape[] = {self.data.size()};
            return nb::ndarray<nb::numpy, double>(self.data.data(), 1, shape, self_obj).cast();
        }, "Access to the underlying data using a memory view")
        ;''',
        "includes": ["<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>"],
    },
    "OSSpectrum": {
        "binding": '''
    using OSSpec = OpenSwath::OSSpectrum;
    using OSBDA = OpenSwath::OSBinaryDataArray;
    using BDAPtr = std::shared_ptr<OpenSwath::BinaryDataArray>;
    nb::class_<OSSpec>(m, "OSSpectrum")
        .def(nb::init<>())
        .def(nb::init<const OSSpec&>())
        .def("getMZArray", &OSSpec::getMZArray)
        .def("getIntensityArray", &OSSpec::getIntensityArray)
        .def("getDriftTimeArray", &OSSpec::getDriftTimeArray)
        .def("setMZArray", &OSSpec::setMZArray, "data"_a)
        .def("setIntensityArray", &OSSpec::setIntensityArray, "data"_a)
        .def("setDriftTimeArray", &OSSpec::setDriftTimeArray, "data"_a)
        .def("set_mz_array", [](OSSpec& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA>();
            arr->data = std::move(data);
            self.setMZArray(arr);
        }, "data"_a, "Set m/z array from list")
        .def("set_intensity_array", [](OSSpec& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")
        .def("get_mz_array", [](const OSSpec& self) {
            auto arr = self.getMZArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get m/z array as list")
        .def("get_mz_array_mv", [](nb::object self_obj) -> nb::object {
            auto& self = nb::cast<OSSpec&>(self_obj);
            auto& data = self.getMZArray()->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, self_obj).cast();
        }, "Get m/z array as writable memory view")
        .def("get_intensity_array", [](const OSSpec& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array as list")
        .def("get_intensity_array_mv", [](nb::object self_obj) -> nb::object {
            auto& self = nb::cast<OSSpec&>(self_obj);
            auto& data = self.getIntensityArray()->data;
            size_t shape[] = {data.size()};
            return nb::ndarray<nb::numpy, double>(data.data(), 1, shape, self_obj).cast();
        }, "Get intensity array as writable memory view")
        .def("get_drift_time_array", [](const OSSpec& self) -> nb::object {
            auto arr = self.getDriftTimeArray();
            if (!arr) return nb::none();
            return nb::cast(arr->data);
        }, "Get drift time array or None")
        .def("get_data_arrays", [](OSSpec& self) {
            auto& arrays = self.getDataArrays();
            std::vector<std::shared_ptr<OSBDA>> result;
            for (auto& a : arrays) result.push_back(a);
            return result;
        }, "Get all data arrays")
        .def("set_data_arrays", [](OSSpec& self, std::vector<std::shared_ptr<OSBDA>> arrays) {
            std::vector<BDAPtr> ptrs;
            for (auto& a : arrays) ptrs.push_back(a);
            self.setDataArrays(ptrs);
        }, "arrays"_a, "Set all data arrays")
        ;''',
        "includes": ["<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>"],
    },
    "OSChromatogram": {
        "binding": '''
    using OSChrom = OpenSwath::OSChromatogram;
    using OSBDA_C = OpenSwath::OSBinaryDataArray;
    using BDAPtr_C = std::shared_ptr<OpenSwath::BinaryDataArray>;
    nb::class_<OSChrom>(m, "OSChromatogram")
        .def(nb::init<>())
        .def(nb::init<const OSChrom&>())
        .def("getTimeArray", &OSChrom::getTimeArray)
        .def("getIntensityArray", &OSChrom::getIntensityArray)
        .def("setTimeArray", &OSChrom::setTimeArray, "data"_a)
        .def("setIntensityArray", &OSChrom::setIntensityArray, "data"_a)
        .def("set_time_array", [](OSChrom& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA_C>();
            arr->data = std::move(data);
            self.setTimeArray(arr);
        }, "data"_a, "Set time array from list")
        .def("set_intensity_array", [](OSChrom& self, std::vector<double> data) {
            auto arr = std::make_shared<OSBDA_C>();
            arr->data = std::move(data);
            self.setIntensityArray(arr);
        }, "data"_a, "Set intensity array from list")
        .def("get_time_array", [](OSChrom& self) {
            auto arr = self.getTimeArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get time array as list")
        .def("get_intensity_array", [](OSChrom& self) {
            auto arr = self.getIntensityArray();
            if (!arr) return std::vector<double>();
            return arr->data;
        }, "Get intensity array as list")
        .def("get_data_arrays", [](OSChrom& self) {
            auto& arrays = self.getDataArrays();
            std::vector<std::shared_ptr<OSBDA_C>> result;
            for (auto& a : arrays) result.push_back(a);
            return result;
        }, "Get all data arrays")
        .def("set_data_arrays", [](OSChrom& self, std::vector<std::shared_ptr<OSBDA_C>> arrays) {
            std::vector<BDAPtr_C> ptrs;
            for (auto& a : arrays) ptrs.push_back(a);
            self.setDataArrays(ptrs);
        }, "arrays"_a, "Set all data arrays")
        ;''',
        "includes": ["<OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>"],
    },
    "IsobaricChannelInformation": {
        "binding": '''
    using ICI = OpenMS::IsobaricQuantitationMethod::IsobaricChannelInformation;
    nb::class_<ICI>(m, "IsobaricChannelInformation")
        .def("__init__", [](ICI* self, const std::string& name, int id, const std::string& description, double center, std::vector<int> affected_channels) {
            new (self) ICI(name, id, description, center, affected_channels);
        }, "name"_a, "id"_a, "description"_a, "center"_a, "affected_channels"_a)
        .def(nb::init<const ICI&>())
        .def_rw("name", &ICI::name)
        .def_rw("id", &ICI::id)
        .def_rw("description", &ICI::description)
        .def_rw("center", &ICI::center)
        .def_rw("affected_channels", &ICI::affected_channels)
        ;''',
        "includes": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>"],
    },
    "SpectrumHelper": {
        "binding": '''
    // SpectrumHelper is a namespace-level helper class with only static methods
    struct SpectrumHelper_Dummy {};
    nb::class_<SpectrumHelper_Dummy>(m, "SpectrumHelper")
        .def_static("removePeaks", [](OpenMS::MSSpectrum& spec, double min_mz, double max_mz) {
            OpenMS::removePeaks(spec, min_mz, max_mz);
        }, "spectrum"_a, "min_mz"_a, "max_mz"_a, "Remove peaks outside the given m/z range")
        .def_static("subtractMinimumIntensity", [](OpenMS::MSSpectrum& spec) {
            OpenMS::subtractMinimumIntensity(spec);
        }, "spectrum"_a, "Subtract the minimum intensity from all peaks")
        ;''',
        "includes": ["<OpenMS/KERNEL/SpectrumHelper.h>", "<OpenMS/KERNEL/MSSpectrum.h>"],
    },
    # TransformationModelBSpline - only static getDefaultParameters needed by tests
    "TransformationModelBSpline": {
        "binding": '''
    struct TransformationModelBSpline_Dummy {};
    nb::class_<TransformationModelBSpline_Dummy>(m, "TransformationModelBSpline")
        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelBSpline::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>", "<OpenMS/KERNEL/StandardTypes.h>"],
    },
    # TransformationModelLowess - only static getDefaultParameters needed by tests
    "TransformationModelLowess": {
        "binding": '''
    struct TransformationModelLowess_Dummy {};
    nb::class_<TransformationModelLowess_Dummy>(m, "TransformationModelLowess")
        .def_static("getDefaultParameters", [](OpenMS::Param& params) {
            OpenMS::TransformationModelLowess::getDefaultParameters(params);
        }, "params"_a, "Get default parameters")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>", "<OpenMS/KERNEL/StandardTypes.h>"],
    },
    # NASequence - minimal binding for fromString static method
    "NASequence": {
        "binding": '''
    nb::class_<OpenMS::NASequence>(m, "NASequence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::NASequence&>())
        .def("toString", &OpenMS::NASequence::toString, "Get string representation")
        .def("__str__", [](const OpenMS::NASequence& self) { return self.toString(); })
        .def("size", &OpenMS::NASequence::size, "Get number of residues")
        .def("empty", &OpenMS::NASequence::empty, "Check if empty")
        .def("getFormula", [](const OpenMS::NASequence& self) { return self.getFormula(); }, "Get empirical formula")
        .def("getMonoWeight", [](const OpenMS::NASequence& self) { return self.getMonoWeight(); }, "Get monoisotopic weight")
        .def("getAverageWeight", [](const OpenMS::NASequence& self) { return self.getAverageWeight(); }, "Get average weight")
        .def_static("fromString", [](const std::string& s) { return OpenMS::NASequence::fromString(s); }, "s"_a, "Create NASequence from string")
        .def("__eq__", [](const OpenMS::NASequence& self, const OpenMS::NASequence& other) { return self == other; }, "other"_a)
        .def("__ne__", [](const OpenMS::NASequence& self, const OpenMS::NASequence& other) { return self != other; }, "other"_a)
        ;''',
        "includes": ["<OpenMS/CHEMISTRY/NASequence.h>"],
    },
    # SpectrumMetaDataLookup - static methods for adding missing RTs/references
    "SpectrumMetaDataLookup": {
        "binding": '''
    struct SpectrumMetaDataLookup_Dummy {};

    auto smld_addMissingRTs = [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::MSExperiment& exp) -> bool {
        return OpenMS::SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptides, exp);
    };

    auto smld_addMissingRefs = [](OpenMS::PeptideIdentificationList& peptides,
            const std::string& filename, bool stop_on_error, bool override_spectra_data, bool override_spectra_references) -> bool {
        std::vector<OpenMS::ProteinIdentification> proteins;
        return OpenMS::SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, stop_on_error, override_spectra_data, override_spectra_references, proteins);
    };

    nb::class_<SpectrumMetaDataLookup_Dummy>(m, "SpectrumMetaDataLookup")
        .def_static("addMissingRTsToPeptideIDs", smld_addMissingRTs, "peptides"_a, "exp"_a, "Add missing RTs to peptide IDs")
        .def_static("addMissingSpectrumReferences", smld_addMissingRefs, "peptides"_a, "filename"_a, "stop_on_error"_a = false, "override_spectra_data"_a = false, "override_spectra_references"_a = false, "Add missing spectrum references")
        ;''',
        "includes": ["<OpenMS/METADATA/SpectrumMetaDataLookup.h>", "<OpenMS/KERNEL/MSExperiment.h>", "<OpenMS/METADATA/PeptideIdentificationList.h>"],
    },
    # IsobaricChannelExtractor
    "IsobaricChannelExtractor": {
        "binding": '''
    nb::class_<OpenMS::IsobaricChannelExtractor, OpenMS::DefaultParamHandler>(m, "IsobaricChannelExtractor")
        .def(nb::init<const OpenMS::ItraqFourPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::ItraqEightPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTSixPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTTenPlexQuantitationMethod*>(), "quant_method"_a)
        .def("extractChannels", &OpenMS::IsobaricChannelExtractor::extractChannels, "ms_exp_data"_a, "consensus_map"_a, "Extract isobaric channels")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>"],
    },
    # IsobaricNormalizer
    "IsobaricNormalizer": {
        "binding": '''
    nb::class_<OpenMS::IsobaricNormalizer>(m, "IsobaricNormalizer")
        .def(nb::init<const OpenMS::ItraqFourPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::ItraqEightPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTSixPlexQuantitationMethod*>(), "quant_method"_a)
        .def(nb::init<const OpenMS::TMTTenPlexQuantitationMethod*>(), "quant_method"_a)
        .def("normalize", &OpenMS::IsobaricNormalizer::normalize, "consensus_map"_a, "Normalize consensus map")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>",
                      "<OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>"],
    },
    # SpectrumAccessOpenMS
    "SpectrumAccessOpenMS": {
        "binding": '''
    nb::class_<OpenMS::SpectrumAccessOpenMS>(m, "SpectrumAccessOpenMS")
        .def(nb::init<std::shared_ptr<OpenMS::MSExperiment>>(), "ms_experiment"_a)
        .def("getSpectrumById", &OpenMS::SpectrumAccessOpenMS::getSpectrumById, "id"_a, "Get spectrum by index")
        .def("getChromatogramById", &OpenMS::SpectrumAccessOpenMS::getChromatogramById, "id"_a, "Get chromatogram by index")
        .def("getNrSpectra", &OpenMS::SpectrumAccessOpenMS::getNrSpectra, "Get number of spectra")
        .def("getNrChromatograms", &OpenMS::SpectrumAccessOpenMS::getNrChromatograms, "Get number of chromatograms")
        .def("getSpectrumMetaById", &OpenMS::SpectrumAccessOpenMS::getSpectrumMetaById, "id"_a, "Get spectrum metadata by index")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>",
                      "<OpenMS/KERNEL/MSExperiment.h>"],
    },
    # SpectrumAccessOpenMSInMemory
    "SpectrumAccessOpenMSInMemory": {
        "binding": '''
    nb::class_<OpenMS::SpectrumAccessOpenMSInMemory>(m, "SpectrumAccessOpenMSInMemory")
        .def("__init__", [](OpenMS::SpectrumAccessOpenMSInMemory* self, OpenMS::SpectrumAccessOpenMS& other) {
            new (self) OpenMS::SpectrumAccessOpenMSInMemory(other);
        }, "other"_a)
        .def("getSpectrumById", &OpenMS::SpectrumAccessOpenMSInMemory::getSpectrumById, "id"_a)
        .def("getChromatogramById", &OpenMS::SpectrumAccessOpenMSInMemory::getChromatogramById, "id"_a)
        .def("getNrSpectra", &OpenMS::SpectrumAccessOpenMSInMemory::getNrSpectra)
        .def("getNrChromatograms", &OpenMS::SpectrumAccessOpenMSInMemory::getNrChromatograms)
        ;''',
        "includes": ["<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>"],
    },
    # SwathMap
    "SwathMap": {
        "binding": '''
    nb::class_<OpenSwath::SwathMap>(m, "SwathMap")
        .def(nb::init<>())
        .def(nb::init<const OpenSwath::SwathMap&>())
        .def(nb::init<double, double, double, bool>(), "mz_start"_a, "mz_end"_a, "mz_center"_a, "is_ms1"_a)
        .def_rw("lower", &OpenSwath::SwathMap::lower)
        .def_rw("upper", &OpenSwath::SwathMap::upper)
        .def_rw("center", &OpenSwath::SwathMap::center)
        .def_rw("ms1", &OpenSwath::SwathMap::ms1)
        .def("setSpectrumPtr", [](OpenSwath::SwathMap& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMS> sa) {
            self.sptr = sa;
        }, "spectrum_access"_a, "Set spectrum access pointer (SpectrumAccessOpenMS)")
        .def("setSpectrumPtr", [](OpenSwath::SwathMap& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMSInMemory> sa) {
            self.sptr = sa;
        }, "spectrum_access"_a, "Set spectrum access pointer (SpectrumAccessOpenMSInMemory)")
        .def("getSpectrumPtr", [](OpenSwath::SwathMap& self) -> nb::object {
            if (!self.sptr) return nb::none();
            // Try dynamic_pointer_cast to known types
            auto sa = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOpenMS>(self.sptr);
            if (sa) return nb::cast(sa);
            auto inmem = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOpenMSInMemory>(self.sptr);
            if (inmem) return nb::cast(inmem);
            return nb::none();
        }, "Get spectrum access pointer")
        ;''',
        "includes": ["<OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>"],
    },
    # ExtractionCoordinates (nested in ChromatogramExtractorAlgorithm)
    # DIAScoring
    "DIAScoring": {
        "binding": '''
    nb::class_<OpenMS::DIAScoring, OpenMS::DefaultParamHandler>(m, "DIAScoring")
        .def(nb::init<>())
        .def("score_with_isotopes", [](OpenMS::DIAScoring& self,
                std::vector<std::shared_ptr<OpenSwath::OSSpectrum>>& spectrum,
                std::vector<OpenSwath::LightTransition>& transitions,
                OpenMS::RangeMobility& im_range,
                double dotprod, double manhattan) {
            self.score_with_isotopes(spectrum, transitions, im_range, dotprod, manhattan);
            return nb::make_tuple(dotprod, manhattan);
        }, "spectrum"_a, "transitions"_a, "im_range"_a, "dotprod"_a, "manhattan"_a)
        .def("dia_by_ion_score", [](OpenMS::DIAScoring& self,
                std::vector<std::shared_ptr<OpenSwath::OSSpectrum>>& spectrum,
                OpenMS::AASequence& sequence,
                int charge,
                OpenMS::RangeMobility& im_range,
                double bseries_score, double yseries_score) {
            self.dia_by_ion_score(spectrum, sequence, charge, im_range, bseries_score, yseries_score);
            return nb::make_tuple(bseries_score, yseries_score);
        }, "spectrum"_a, "sequence"_a, "charge"_a, "im_range"_a, "bseries_score"_a, "yseries_score"_a,
        "Score the DIA window for b/y ion series")
        ;''',
        "includes": ["<OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>",
                      "<OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>"],
    },
    # MRMFeature
    "MRMFeature": {
        "binding": '''
    nb::class_<OpenMS::MRMFeature, OpenMS::Feature>(m, "MRMFeature")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeature&>())
        .def("getScores", [](const OpenMS::MRMFeature& self) { return self.getScores(); })
        .def("setScores", &OpenMS::MRMFeature::setScores, "s"_a)
        .def("getFeature", [](const OpenMS::MRMFeature& self, const OpenMS::String& key) { return self.getFeature(key); }, "key"_a)
        .def("addFeature", [](OpenMS::MRMFeature& self, const OpenMS::Feature& f, const OpenMS::String& key) { self.addFeature(f, key); }, "f"_a, "key"_a)
        .def("getFeatures", [](const OpenMS::MRMFeature& self) { return self.getFeatures(); })
        .def("getFeatureIDs", [](OpenMS::MRMFeature& self) {
            std::vector<OpenMS::String> result;
            self.getFeatureIDs(result);
            return result;
        })
        .def("getPrecursorFeature", [](const OpenMS::MRMFeature& self, const OpenMS::String& key) { return self.getPrecursorFeature(key); }, "key"_a)
        .def("addPrecursorFeature", [](OpenMS::MRMFeature& self, const OpenMS::Feature& f, const OpenMS::String& key) { self.addPrecursorFeature(f, key); }, "f"_a, "key"_a)
        .def("getPrecursorFeatureIDs", [](OpenMS::MRMFeature& self) {
            std::vector<OpenMS::String> result;
            self.getPrecursorFeatureIDs(result);
            return result;
        })
        .def("__eq__", [](const OpenMS::MRMFeature& a, const OpenMS::MRMFeature& b) { return a == b; })
        .def("__ne__", [](const OpenMS::MRMFeature& a, const OpenMS::MRMFeature& b) { return a != b; })
        ;''',
        "includes": ["<OpenMS/KERNEL/MRMFeature.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>"],
    },
    # MRMTransitionGroupCP (template specialization)
    # MRMFeatureFinderScoring
    "MRMFeatureFinderScoring": {
        "binding": '''
    nb::class_<OpenMS::MRMFeatureFinderScoring, OpenMS::DefaultParamHandler>(m, "MRMFeatureFinderScoring")
        .def(nb::init<>())
        .def("pickExperiment", [](OpenMS::MRMFeatureFinderScoring& self,
                OpenMS::MSExperiment& chromatograms,
                OpenMS::FeatureMap& output,
                OpenMS::TargetedExperiment& transition_exp,
                OpenMS::TransformationDescription& trafo,
                OpenMS::MSExperiment& swath_map) {
            self.pickExperiment(chromatograms, output, transition_exp, trafo, swath_map);
        }, "chromatograms"_a, "output"_a, "transition_exp"_a, "trafo"_a, "swath_map"_a)
        .def("setStrictFlag", &OpenMS::MRMFeatureFinderScoring::setStrictFlag, "flag"_a)
        .def("setMS1Map", [](OpenMS::MRMFeatureFinderScoring& self, std::shared_ptr<OpenMS::SpectrumAccessOpenMS> ms1_map) {
            self.setMS1Map(ms1_map);
        }, "ms1_map"_a)
        ;''',
        "includes": ["<OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>",
                      "<OpenMS/KERNEL/MSExperiment.h>",
                      "<OpenMS/KERNEL/FeatureMap.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>"],
    },
    "OpenSwathScoring": {
        "binding": '''
    // OpenSwathScoring
    nb::class_<OpenMS::OpenSwathScoring>(m, "OpenSwathScoring")
        .def(nb::init<>())
        .def("initialize", &OpenMS::OpenSwathScoring::initialize,
            "rt_normalization_factor"_a, "add_up_spectra"_a,
            "spacing_for_spectra_resampling"_a, "merge_spectra_by_peak_width_fraction"_a,
            "drift_extra"_a, "su"_a, "spectrum_addition_method"_a,
            "spectrum_merge_method_type"_a, "use_ms1_ion_mobility"_a,
            "apply_im_peak_picking"_a)
        ;''',
        "includes": ["<OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>",
                      "<OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>"],
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


class NanobindEmitter:
    """
    Generates nanobind C++ binding code from merged C++/pxd declarations.
    Uses accurate C++ type information from libclang.
    """

    def __init__(self, single_module: bool = False):
        self.single_module = single_module

        # Standard includes for all modules
        self._standard_includes = {
            "<nanobind/nanobind.h>",
            "<nanobind/operators.h>",
            # Note: <nanobind/stl/string.h> is NOT included here because we use
            # a custom std::string caster (std_string_bytes_caster.h via all_casters.h)
            # that also accepts Python bytes, matching pyOpenMS/Cython behavior.
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

        # Distribute classes across domain-based modules
        domains = self._distribute_classes(classes)

        if self.single_module:
            # Single-module mode: emit everything into main_module.cpp
            all_classes = [cls for domain_classes in domains.values() for cls in domain_classes]
            content = self._generate_module_content(all_classes, classes, domain_name=None)
            self._write_single_module(output_dir / "main_module.cpp", content)
        else:
            # Multi-module mode: emit bind_<domain>.cpp files + placeholder main
            for domain_name in DOMAIN_NAMES:
                domain_classes = domains[domain_name]
                content = self._generate_module_content(domain_classes, classes, domain_name=domain_name)
                output_file = output_dir / f"bind_{domain_name}.cpp"
                self._write_module_file(output_file, domain_name, content)
            self._write_main_module(output_dir / "main_module.cpp")

    def _distribute_classes(
        self, classes: Dict[str, MergedClass]
    ) -> Dict[str, List[MergedClass]]:
        """Distribute classes across domain-based modules.

        Each class is assigned to a domain based on its OpenMS header path
        (e.g., OpenMS/KERNEL/Peak1D.h -> kernel). Cross-module inheritance
        is handled by nanobind's NB_DOMAIN type sharing.
        """
        caster_types = get_caster_owned_types()

        # Collect all bindable classes (skip nested, skipped, and caster-owned)
        bindable: Dict[str, MergedClass] = {}
        class_domain: Dict[str, str] = {}

        for merged_class in classes.values():
            class_name = merged_class.name
            if class_name in SKIP_CLASSES or class_name in caster_types:
                continue
            if "::" in class_name:
                continue

            bindable[class_name] = merged_class
            class_domain[class_name] = _get_domain_from_header(merged_class.header_file)

        # Identify classes that would be skipped (abstract or private ctor)
        # but are needed as base classes by other bindable classes (transitively)
        base_classes_needed: Set[str] = set()
        changed = True
        while changed:
            changed = False
            for mc in bindable.values():
                # Skip classes that would be filtered out, unless already rescued
                would_skip = mc.is_abstract or mc.has_private_constructor
                if would_skip and mc.name not in base_classes_needed:
                    continue
                for base in getattr(mc, 'base_classes', []):
                    base_name = _unqualified_name(base)
                    if base.startswith('std::'):
                        continue
                    if base_name in bindable:
                        base_mc = bindable[base_name]
                        base_would_skip = base_mc.is_abstract or base_mc.has_private_constructor
                        if base_would_skip and base_name not in base_classes_needed:
                            base_classes_needed.add(base_name)
                            changed = True

        # Remove abstract classes that aren't needed as bases.
        # Note: private-ctor classes are NOT filtered here because they may still
        # be bindable (singletons, SPECIAL_METHODS). They're filtered later in
        # _generate_class_binding with more nuanced checks.
        bindable = {
            name: mc for name, mc in bindable.items()
            if not mc.is_abstract or name in base_classes_needed
        }

        # Union-find to group inheritance chains, then move entire chain
        # to the root base class's domain. This ensures base classes are
        # registered before derived classes within the same NB_MODULE.
        parent: Dict[str, str] = {}

        def find(x: str) -> str:
            if x not in parent:
                parent[x] = x
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(x: str, y: str) -> None:
            px, py = find(x), find(y)
            if px != py:
                # Prefer the base class (y) as root so its domain wins
                parent[px] = py

        for mc in bindable.values():
            find(mc.name)
            for base in getattr(mc, 'base_classes', []):
                base_name = _unqualified_name(base)
                if base.startswith('std::'):
                    continue
                if base_name in bindable:
                    union(mc.name, base_name)  # base becomes root

        # Assign each class to the domain of its inheritance chain root
        for class_name in bindable:
            root = find(class_name)
            class_domain[class_name] = class_domain[root]

        # Build domain lists
        domains: Dict[str, List[MergedClass]] = {name: [] for name in DOMAIN_NAMES}
        for class_name, mc in bindable.items():
            domains[class_domain[class_name]].append(mc)

        # Store abstract base classes needed for use in _generate_class_bindings
        self._base_classes_needed = base_classes_needed

        # Store class->domain mapping for handwritten class assignment
        self._class_domain = class_domain

        # Store the set of all class names that will be bound
        self._bound_class_names: Set[str] = set()
        for domain_classes in domains.values():
            for cls in domain_classes:
                self._bound_class_names.add(cls.name)

        # Log domain distribution
        for domain_name in DOMAIN_NAMES:
            count = len(domains[domain_name])
            if count > 0:
                logger.info(f"  Domain '{domain_name}': {count} classes")

        return domains

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
                base_name = _unqualified_name(base)

                # Only track dependencies on classes we're binding in this module
                if (base_name in class_names and
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
        domain_name: Optional[str] = None,
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

        # Collect enums from merged classes that are in this module
        # Enums attached to classes are now generated as nested enums within the class binding
        # (see _generate_class_binding). We track their names here to avoid double-binding.
        generated_enum_names = set()
        for merged_class in module_classes:
            class_name = merged_class.name
            if class_name not in class_names_in_module:
                continue
            # Skip problematic classes
            if class_name in SKIP_CLASSES or class_name in get_caster_owned_types():
                continue

            # Mark attached enums as generated (they'll be nested under the class)
            # But only if they're not in SKIP_ENUMS (those will be skipped and need fallback)
            for enum_decl in merged_class.enums:
                if enum_decl.name not in SKIP_ENUMS:
                    generated_enum_names.add(enum_decl.name)

        # Add hard-coded enum bindings from SPECIAL_METHODS as fallback
        # (for enums not yet in pxd files or with complex definitions)
        if "__enums__" in SPECIAL_METHODS:
            # Map enum names to the class they should be bound with
            # Note: Most enums are now auto-generated from pxd files
            enum_class_map = {
                "DriftTimeUnit": "MSSpectrum",  # Not attached to class in pxd
                "IMFormat": "MSSpectrum",  # Not attached to class in pxd
                "Pi0Method": "MultipleTesting",  # Nested in MultipleTesting
                "LfdrTransform": "MultipleTesting",  # Nested in MultipleTesting
                "Method": "RankData",  # Nested in RankData
                "NaNPolicy": "RankData",  # Nested in RankData
            }
            for enum_name, enum_code in SPECIAL_METHODS["__enums__"].items():
                # Skip if already generated from pxd
                if enum_name in generated_enum_names:
                    continue
                target_class = enum_class_map.get(enum_name)
                if target_class and target_class in class_names_in_module:
                    content.enum_bindings.append(enum_code)
                    # Add includes for enums that need specific headers
                    enum_include_map = {
                        "Method": "<OpenMS/MATH/STATISTICS/RankData.h>",
                        "NaNPolicy": "<OpenMS/MATH/STATISTICS/RankData.h>",
                    }
                    if enum_name in enum_include_map:
                        content.includes.add(enum_include_map[enum_name])

        # Add nested enum bindings (SpectrumType) to the module with SpectrumSettings
        # Only add if not already auto-generated from pxd
        if "SpectrumSettings" in class_names_in_module:
            content.includes.add("<OpenMS/METADATA/SpectrumSettings.h>")
            if "__post_class_enums__" in SPECIAL_METHODS:
                for enum_name, enum_code in SPECIAL_METHODS["__post_class_enums__"].items():
                    # Skip if this enum was already auto-generated
                    simple_name = enum_name.split("_")[-1] if "_" in enum_name else enum_name
                    if simple_name in generated_enum_names:
                        logger.debug(f"Skipping hardcoded enum {enum_name} - already auto-generated")
                        continue
                    content.post_class_enums.append(enum_code)

        for merged_class in sorted_module_classes:
            class_name = merged_class.name

            # Skip problematic classes and types with casters
            if class_name in SKIP_CLASSES or class_name in get_caster_owned_types():
                continue

            # Skip nested classes (but allow nested namespaces)
            if "::" in class_name:
                continue

            # Skip abstract classes UNLESS they're needed as base classes
            if merged_class.is_abstract:
                if not (hasattr(self, '_base_classes_needed') and
                        class_name in self._base_classes_needed):
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

                # Emit nested class bindings associated with this parent class
                if class_name in NESTED_CLASS_BINDINGS:
                    for binding_code, extra_includes in NESTED_CLASS_BINDINGS[class_name]:
                        content.class_bindings.append(binding_code)
                        for inc in extra_includes:
                            content.includes.add(inc)

        # Emit handwritten full class bindings.
        # If a handwritten class inherits from a bound class, place it in
        # the same domain as the base (which may have been moved by union-find).
        _hw_base_re = re.compile(r'nb::class_<[^,>]+,\s*OpenMS::(\w+)>')
        class_domain_map = getattr(self, '_class_domain', {})
        for hw_class_name, hw_info in HANDWRITTEN_CLASSES.items():
            hw_domain = _get_handwritten_class_domain(hw_info)
            m = _hw_base_re.search(hw_info.get("binding", ""))
            if m and m.group(1) in class_domain_map:
                hw_domain = class_domain_map[m.group(1)]
            if hw_domain == domain_name:
                content.class_bindings.append(hw_info["binding"])
                for inc in hw_info.get("includes", []):
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


        # Handle template classes with wrap-instances (before SKIP_CLASSES check,
        # since the template base may be skipped but instances should be generated)
        if merged_class.template_instances and merged_class.cpp_class.is_template:
            return self._generate_template_instances(merged_class, module_class_names)

        # Skip problematic classes and types with casters
        if class_name in SKIP_CLASSES or class_name in get_caster_owned_types():
            return None

        # Skip nested classes (but allow nested namespaces like OpenMS::DataArrays)
        # Nested classes have :: in the class name itself, not just in the namespace
        if "::" in class_name:
            return None

        # Skip abstract classes UNLESS they're needed as base classes
        is_abstract_base = (merged_class.is_abstract and
                           hasattr(self, '_base_classes_needed') and
                           class_name in self._base_classes_needed)
        if merged_class.is_abstract and not is_abstract_base:
            logger.debug(f"Skipping abstract class: {class_name}")
            return None

        # Classes with deleted default constructors: don't skip the class,
        # just skip the default constructor in _generate_constructor.
        # These classes can still be useful as return types or have copy/parameterized ctors.
        if merged_class.has_deleted_default_constructor:
            logger.debug(f"Class {class_name} has deleted default constructor - will skip default ctor only")

        # Skip classes with only private/protected constructors (singletons, etc.)
        # Exception: allow if class has SPECIAL_METHODS or is a singleton (getInstance pattern)
        is_singleton = self._is_singleton(merged_class)
        is_needed_base = (hasattr(self, '_base_classes_needed') and
                          class_name in self._base_classes_needed)
        if merged_class.has_private_constructor:
            if class_name not in SPECIAL_METHODS and not is_singleton and not is_needed_base:
                logger.debug(f"Skipping class with private constructors: {class_name}")
                return None
            else:
                logger.debug(f"Class {class_name} has private constructors but has SPECIAL_METHODS/singleton - binding anyway")

        # Handle base classes - nanobind supports multiple inheritance
        # We can only specify base classes that are also bound to nanobind AND
        # that are in the same module (cross-module base classes aren't registered yet at init time)
        bound_bases = self._get_bound_base_classes(merged_class, module_class_names)
        if bound_bases:
            # nanobind syntax: nb::class_<Derived, Base1, Base2>(m, "Derived")
            base_spec = ", " + ", ".join(bound_bases)
        else:
            base_spec = ""

        # Check if class has attached enums - if so, we need to capture the class
        # binding to use as the scope for nested enum definitions
        has_nested_enums = bool(merged_class.enums)
        class_var = f"{class_name.lower()}_class" if has_nested_enums else None

        # Start class definition
        if has_nested_enums:
            lines.append(f'    auto {class_var} = nb::class_<{qualified_name}{base_spec}>(m, "{class_name}")')
        else:
            lines.append(f'    nb::class_<{qualified_name}{base_spec}>(m, "{class_name}")')

        # Add docstring
        if merged_class.doc:
            doc = self._escape_string(merged_class.doc)
            lines[-1] = lines[-1][:-1]  # Remove trailing paren
            lines[-1] += f', "{doc}")'

        # Generate constructors (skip for abstract/private-ctor base classes - they can't be instantiated)
        self._current_merged_class = merged_class
        if not is_abstract_base and not (is_needed_base and merged_class.has_private_constructor):
            for ctor in merged_class.constructors:
                ctor_code = self._generate_constructor(ctor, qualified_name)
                if ctor_code:
                    lines.append(f"        {ctor_code}")

        # Generate methods
        gen_method_count = 0
        for merged_method in merged_class.methods:
            if merged_method.wrap_ignore:
                continue
            method_code = self._generate_method(merged_method, qualified_name, merged_class)
            if method_code:
                lines.append(f"        {method_code}")
                gen_method_count += 1
        if class_name == "MSSpectrum":
            logger.info(f"MSSpectrum: {len(merged_class.methods)} merged, {gen_method_count} generated, lines so far: {len(lines)}")

        # Add explicit inherited method bindings for classes where we skipped base class
        # specification due to non-virtual destructor issues
        inherited_methods = self._get_explicit_inherited_methods(merged_class, qualified_name)
        for method_line in inherited_methods:
            lines.append(f"        {method_line}")

        # Add member variable bindings (public fields from pxd)
        for mv in getattr(merged_class, 'member_variables', []):
            lines.append(f'        .def_rw("{mv.name}", &{qualified_name}::{mv.name})')

        # Phase 5: Auto-generate public field bindings from libclang
        pxd_field_names = {mv.name for mv in getattr(merged_class, 'member_variables', [])}
        special_field_names = set(SPECIAL_METHODS.get(class_name, {}).keys())
        skip_field_names = SKIP_METHODS.get(class_name, set())
        for cpp_field in getattr(merged_class.cpp_class, 'public_fields', []):
            if cpp_field.name in pxd_field_names:
                continue  # Already emitted from pxd
            if cpp_field.name in special_field_names or cpp_field.name in skip_field_names:
                continue
            # Skip fields with types that nanobind cannot handle
            field_type = cpp_field.type_str
            canonical_field = cpp_field.canonical_type or field_type
            if any(pat in field_type or pat in canonical_field
                   for pat in _UNBINDABLE_FIELD_PATTERNS):
                continue
            if cpp_field.is_const:
                lines.append(f'        .def_ro("{cpp_field.name}", &{qualified_name}::{cpp_field.name})')
            else:
                lines.append(f'        .def_rw("{cpp_field.name}", &{qualified_name}::{cpp_field.name})')

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

        # Collect which dunders SPECIAL_METHODS will provide (to avoid duplicates)
        special_dunders = set()
        if class_name in SPECIAL_METHODS:
            for method_name in SPECIAL_METHODS[class_name]:
                if method_name.startswith("__") and method_name.endswith("__"):
                    special_dunders.add(method_name)

        # Add iterator support (detected from begin/end methods or wrap-iter directive)
        if has_iterator and "__iter__" not in special_dunders:
            lines.append(f'        .def("__iter__", []({qualified_name}& self) {{ return nb::make_iterator<nb::rv_policy::reference_internal>(nb::handle(), "{class_name}_iter", self.begin(), self.end()); }})')

        # Add __len__ for container classes (detected from size() method)
        if has_size and "__len__" not in special_dunders:
            lines.append(f'        .def("__len__", []({qualified_name}& self) {{ return self.size(); }})')

        # Add __getitem__ for vector-based classes (detected from std::vector inheritance)
        if is_vector_based and "__getitem__" not in special_dunders:
            lines.append(f'        .def("__getitem__", []({qualified_name}& self, size_t i) -> {vector_elem_type}& {{ ')
            lines.append(f'            if (i >= self.size()) throw nb::index_error();')
            lines.append(f'            return self[i];')
            lines.append(f'        }}, nb::rv_policy::reference_internal)')

        # Add __getitem__ for non-vector containers with operator[] and size()
        # operator[] may be inherited (e.g. from ExposedVector<T>) and not enumerated
        # by libclang, so we check multiple sources: cpp_class.operators, merged methods,
        # and the pxd-declared method list (which always includes operator[] if declared).
        if has_size and not is_vector_based and "__getitem__" not in special_dunders:
            # Find operator[] and its return type from multiple sources:
            # 1. Direct operators on the class (libclang)
            # 2. Merged methods from pxd (may include inherited operator[])
            # 3. CONTAINER_CLASSES fallback (ExposedVector<T> base)
            bracket_op = None
            for m in getattr(merged_class.cpp_class, 'operators', []):
                if m.name == "operator[]" and not m.is_const:
                    bracket_op = m
                    break
            if bracket_op is None:
                for m in getattr(merged_class.cpp_class, 'operators', []):
                    if m.name == "operator[]":
                        bracket_op = m
                        break
            # Also check merged methods (pxd-declared operator[] with return type)
            if bracket_op is None:
                for m in merged_class.methods:
                    if m.cpp_method.name == "operator[]" and not m.cpp_method.is_const:
                        bracket_op = m.cpp_method
                        break
                if bracket_op is None:
                    for m in merged_class.methods:
                        if m.cpp_method.name == "operator[]":
                            bracket_op = m.cpp_method
                            break
            has_bracket = (
                bracket_op is not None
                or class_name in CONTAINER_CLASSES  # ExposedVector<T> base provides operator[]
            )
            if has_bracket:
                # Determine if operator[] returns a reference (use reference_internal)
                ret = bracket_op.return_type if bracket_op else None
                # For ExposedVector-based CONTAINER_CLASSES without a detected operator[],
                # infer the element type from the class (operator[] returns T&)
                if ret is None and class_name in CONTAINER_CLASSES:
                    elem_type = self._get_exposed_vector_element_type(class_name)
                    if elem_type:
                        ret = f"{elem_type} &"
                is_ref_return = (ret and '&' in ret and '&&' not in ret
                                 and ret.replace('const ', '').replace('&', '').strip() not in _PRIMITIVE_TYPES
                                 and ret.replace('const ', '').replace('&', '').strip() not in _TYPE_CASTER_TYPES)
                if is_ref_return:
                    norm_ret = self._normalize_type(ret, preserve_reference=True,
                                                     canonical_type=getattr(bracket_op, 'canonical_return_type', ret))
                    if '&' not in norm_ret:
                        norm_ret = f"{norm_ret}&"
                    lines.append(f'        .def("__getitem__", []({qualified_name}& self, size_t i) -> {norm_ret} {{ ')
                    lines.append(f'            if (i >= self.size()) throw nb::index_error();')
                    lines.append(f'            return self[i];')
                    lines.append(f'        }}, nb::rv_policy::reference_internal)')
                else:
                    lines.append(f'        .def("__getitem__", []({qualified_name}& self, size_t i) {{ ')
                    lines.append(f'            if (i >= self.size()) throw nb::index_error();')
                    lines.append(f'            return self[i];')
                    lines.append(f'        }})')

        # Add special methods for this class (get_peaks, set_peaks, etc.)
        post_class_code = []
        if class_name in SPECIAL_METHODS:
            for method_name, method_code in SPECIAL_METHODS[class_name].items():
                if method_name == "__post_class__":
                    post_class_code.append(method_code)
                else:
                    lines.append(method_code)

        # __repr__ is now provided by Python addons (pyopenms/addons/)

        # Auto-generate conversion operators (operator bool -> __bool__, etc.)
        self._generate_conversion_operators(lines, merged_class, qualified_name)

        # Auto-generate singleton getInstance if detected and not in SPECIAL_METHODS
        if is_singleton and not (class_name in SPECIAL_METHODS and "getInstance" in SPECIAL_METHODS[class_name]):
            # Read the actual return type from the C++ method to preserve const-correctness
            singleton_ret = f"{qualified_name}*"
            for cpp_m in merged_class.cpp_class.methods:
                if cpp_m.name == "getInstance" and cpp_m.is_static:
                    if cpp_m.return_type and 'const' in cpp_m.return_type:
                        singleton_ret = f"const {qualified_name}*"
                    break
            lines.append(f'        .def_static("getInstance", []() -> {singleton_ret} {{ return {qualified_name}::getInstance(); }}, nb::rv_policy::reference, "Returns the singleton instance")')

        # Close class definition
        lines.append("        ;")

        # Generate nested enum bindings for attached enums
        if has_nested_enums and class_var:
            for enum_decl in merged_class.enums:
                enum_code = self._generate_nested_enum_binding(enum_decl, class_var, qualified_name)
                if enum_code:
                    lines.append(enum_code)

        # Emit post-class code (e.g., nested enum bindings)
        for pc in post_class_code:
            lines.append(pc)

        return "\n".join(lines)

    def _generate_conversion_operators(self, lines: list, merged_class: MergedClass, qualified_name: str):
        """Auto-generate Python dunders from C++ conversion operators."""
        conversion_fns = getattr(merged_class.cpp_class, 'conversion_functions', [])
        if not conversion_fns:
            return

        # Map C++ conversion function return types to Python dunders
        type_to_dunder = {
            "bool": "__bool__",
            "_Bool": "__bool__",
            "int": "__int__",
            "long": "__int__",
            "long long": "__int__",
            "unsigned int": "__int__",
            "size_t": "__int__",
            "double": "__float__",
            "float": "__float__",
            "std::string": "__str__",
            "std::basic_string<char>": "__str__",
        }

        for conv in conversion_fns:
            ret_type = conv.return_type.strip()
            # Also check canonical type
            canonical = getattr(conv, 'canonical_return_type', ret_type).strip()
            dunder = type_to_dunder.get(ret_type) or type_to_dunder.get(canonical)
            if dunder:
                lines.append(f'        .def("{dunder}", []({qualified_name}& self) {{ return static_cast<{ret_type}>(self); }})')

    def _is_singleton(self, merged_class: MergedClass) -> bool:
        """Detect singleton pattern: private constructor + public static getInstance.

        Checks both merged methods (pxd-filtered) and raw C++ methods (bypasses
        wrap-ignore on getInstance).
        """
        if not merged_class.has_private_constructor:
            return False
        # Check merged methods (pxd-allowed)
        for m in merged_class.methods:
            cpp = m.cpp_method
            if cpp.name == "getInstance" and cpp.is_static:
                return True
        # Also check raw C++ methods (bypasses pxd wrap-ignore on getInstance)
        for cpp_m in merged_class.cpp_class.methods:
            if cpp_m.name == "getInstance" and cpp_m.is_static:
                return True
        return False

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

    def _generate_enum_binding(self, enum_decl: Any) -> Optional[str]:
        """Generate nanobind enum binding code from an EnumDecl.

        Args:
            enum_decl: EnumDecl object from pxd_parser

        Returns:
            C++ code string for the enum binding, or None if enum has no values
        """
        if not enum_decl.values:
            return None

        enum_name = enum_decl.name
        namespace = enum_decl.namespace

        # Skip enums that have known issues (e.g., value name mismatches with C++)
        if enum_name in SKIP_ENUMS:
            logger.debug(f"Skipping enum {enum_name} - in SKIP_ENUMS list")
            return None

        # Use explicit C++ name if provided (e.g., "OpenMS::FileTypes::Type"),
        # otherwise construct from namespace and enum name
        if hasattr(enum_decl, 'cpp_name') and enum_decl.cpp_name:
            qualified_name = enum_decl.cpp_name
        else:
            qualified_name = f"{namespace}::{enum_name}"

        # Check if this is a scoped enum (C++11 enum class)
        is_scoped = getattr(enum_decl, 'is_scoped', False)

        lines = [f"    // {enum_name} enum (auto-generated from pxd)"]
        lines.append(f'    nb::enum_<{qualified_name}>(m, "{enum_name}", nb::is_arithmetic())')

        for value in enum_decl.values:
            lines.append(f'        .value("{value.name}", {qualified_name}::{value.name})')

        # Only export values for non-scoped enums (regular enums)
        # Scoped enums (enum class) should keep values scoped
        if not is_scoped:
            lines.append("        .export_values();")
        else:
            lines.append("        ;")

        return "\n".join(lines)

    def _generate_nested_enum_binding(
        self, enum_decl: Any, class_var: str, parent_qualified_name: str
    ) -> Optional[str]:
        """Generate nanobind enum binding nested under a class.

        Args:
            enum_decl: EnumDecl object from pxd_parser
            class_var: Variable name holding the parent class binding
            parent_qualified_name: Fully qualified C++ name of parent class

        Returns:
            C++ code string for the nested enum binding, or None if enum has no values
        """
        if not enum_decl.values:
            return None

        enum_name = enum_decl.name
        namespace = enum_decl.namespace

        # Skip enums that have known issues
        if enum_name in SKIP_ENUMS:
            logger.debug(f"Skipping nested enum {enum_name} - in SKIP_ENUMS list")
            return None

        # Use the enum's actual C++ qualified name from its namespace
        # (not the parent class name, which is for Python nesting only)
        if hasattr(enum_decl, 'cpp_name') and enum_decl.cpp_name:
            cpp_qualified_name = enum_decl.cpp_name
        else:
            cpp_qualified_name = f"{namespace}::{enum_name}"

        # Check if this is a scoped enum (C++11 enum class)
        is_scoped = getattr(enum_decl, 'is_scoped', False)

        parent_class_name = parent_qualified_name.split('::')[-1]
        lines = [f"    // {enum_name} enum nested under {parent_class_name}"]
        lines.append(f'    nb::enum_<{cpp_qualified_name}>({class_var}, "{enum_name}", nb::is_arithmetic())')

        for value in enum_decl.values:
            lines.append(f'        .value("{value.name}", {cpp_qualified_name}::{value.name})')

        # Only export values for non-scoped enums (regular enums)
        if not is_scoped:
            lines.append("        .export_values();")
        else:
            lines.append("        ;")

        return "\n".join(lines)

    def _get_vector_element_type(self, merged_class: MergedClass) -> Optional[str]:
        """Detect if class inherits from std::vector and extract element type.

        Parses canonical base classes looking for std::vector<T> pattern.
        Returns the element type T if found, None otherwise.
        """
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

    def _get_exposed_vector_element_type(self, class_name: str) -> Optional[str]:
        """Get element type for ExposedVector-based container classes.

        These classes use composition (ExposedVector<T>) rather than std::vector
        inheritance, so _get_vector_element_type won't detect them. Their operator[]
        returns T& but isn't visible to libclang as a direct method.
        """
        exposed_vector_types = {
            "MSExperiment": "OpenMS::MSSpectrum",
            "FeatureMap": "OpenMS::Feature",
            "ConsensusMap": "OpenMS::ConsensusFeature",
        }
        return exposed_vector_types.get(class_name)

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
        - Actually being bound (in _bound_class_names)
        - Not in SKIP_CLASSES or have type casters

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

        # Classes with private std::vector<> as first base — the memory layout
        # has std::vector before any public bases, so nanobind's static_cast
        # upcast to the first PUBLIC base gives a wrong pointer offset
        if len(base_classes) > 0 and base_classes[0].startswith('std::'):
            return []

        # For multiple inheritance, check if first base has non-virtual destructor
        if len(base_classes) > 1:
            first_base = base_classes[0]
            first_base_name = _unqualified_name(first_base)
            if first_base_name in NONVIRTUAL_DESTRUCTOR_CLASSES:
                # Skip base class specification to avoid destructor issues
                return []

        for base in base_classes:
            base_name = _unqualified_name(base)

            # Skip std:: types (like std::vector) - they're not bound as classes
            if base.startswith('std::'):
                continue

            # Skip Qt base classes (QDate, QString, QObject, etc.)
            # Qt classes follow pattern: Q + CapitalLetter (but not QC=QualityControl, QT=OpenMS)
            if self._is_qt_class(base_name):
                logger.debug(f"Skipping Qt base class: {base_name}")
                continue

            # Check if this base class is bound (not skipped, no caster)
            if base_name in self._bound_class_names and base_name not in SKIP_CLASSES and base_name not in caster_types:
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
            "Precursor",        # Peak1D + CVTermList
            "Product",          # Peak1D
            "MSSpectrum",       # private std::vector<Peak1D> + SpectrumSettings
            "MSChromatogram",   # private std::vector<ChromatogramPeak> + ChromatogramSettings
        }

        # Check if this class needs explicit inherited methods
        needs_explicit = False

        # Case 0: Class has std:: first base (private vector)
        if len(base_classes) > 0 and base_classes[0].startswith('std::'):
            needs_explicit = True

        # Case 1: Class has any base with non-virtual destructor
        if not needs_explicit:
            for base in base_classes:
                if _unqualified_name(base) in NONVIRTUAL_DESTRUCTOR_CLASSES:
                    needs_explicit = True
                    break

        # Case 2: Class inherits from a class that has skipped base
        if not needs_explicit:
            for base in base_classes:
                if _unqualified_name(base) in CLASSES_NEEDING_EXPLICIT_METHODS:
                    needs_explicit = True
                    break

        if not needs_explicit:
            return []

        methods = []

        # Determine which method sets to include based on inheritance chain
        needs_peak1d = False
        needs_peak2d = False
        needs_metainfo = False
        needs_uniqueid = False
        needs_docid = False

        # Check direct bases
        for base in base_classes:
            base_name = _unqualified_name(base)
            if base_name in {"Peak1D", "Precursor", "Product"}:
                needs_peak1d = True
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
        if class_name in {"Precursor", "Product"}:
            needs_peak1d = True
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

        # Peak1D methods (getMZ, setMZ, getIntensity, setIntensity)
        if needs_peak1d:
            methods.extend([
                f'.def("getMZ", []({qualified_name}& self) {{ return self.getMZ(); }}, "Returns the m/z")',
                f'.def("setMZ", []({qualified_name}& self, double mz) {{ self.setMZ(mz); }}, "mz"_a, "Sets the m/z")',
                f'.def("getIntensity", []({qualified_name}& self) {{ return self.getIntensity(); }}, "Returns the intensity")',
                f'.def("setIntensity", []({qualified_name}& self, float intensity) {{ self.setIntensity(intensity); }}, "intensity"_a, "Sets the intensity")',
                f'.def("getPos", []({qualified_name}& self) {{ return self.getPos(); }}, "Returns the position (m/z)")',
                f'.def("setPos", []({qualified_name}& self, double pos) {{ self.setPos(pos); }}, "pos"_a, "Sets the position (m/z)")',
            ])

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

    # Copy constructor skipping is now handled via AST detection
    # (has_deleted_copy_constructor from cpp_parser.py)

    def _generate_template_instances(
        self, merged_class: MergedClass, module_class_names: Optional[Set[str]] = None
    ) -> Optional[str]:
        """Generate bindings for template class instantiations.

        For a template class like Matrix[ValueT] with wrap-instances:
          MatrixDouble := Matrix[double]
        Generates: nb::class_<OpenMS::Matrix<double>>(m, "MatrixDouble")
        """
        lines = []
        class_name = merged_class.name
        namespace = merged_class.namespace

        # Type mapping from pxd type names to C++ types
        def _convert_template_type(type_str: str, param_map: dict) -> str:
            """Convert a pxd type to C++ type with template substitution."""
            result = type_str.strip()
            # Apply template parameter substitution using word boundaries
            # to avoid partial matches (e.g., 'Key' matching inside 'KeyType')
            for pn, ca in sorted(param_map.items(), key=lambda x: -len(x[0])):
                result = re.sub(r'\b' + re.escape(pn) + r'\b', ca, result)
            # Handle libcpp containers
            m = re.match(r'libcpp_vector\[\s*(.+)\s*\]', result)
            if m:
                inner = _convert_template_type(m.group(1), param_map)
                return f"std::vector<{inner}>"
            m = re.match(r'libcpp_pair\[\s*(.+?)\s*,\s*(.+)\s*\]', result)
            if m:
                return f"std::pair<{_convert_template_type(m.group(1), param_map)}, {_convert_template_type(m.group(2), param_map)}>"
            m = re.match(r'libcpp_set\[\s*(.+)\s*\]', result)
            if m:
                return f"std::set<{_convert_template_type(m.group(1), param_map)}>"
            if result in ("libcpp_utf8_string", "libcpp_string", "libcpp_utf8_output_string"):
                return "std::string"
            # Handle generic pxd ClassName[Args] -> C++ ClassName<Args>
            # (only non-libcpp names; libcpp_* already handled above)
            m = re.match(r'(?!libcpp_)(\w+)\[(.+)\]', result)
            if m:
                cls_name = m.group(1)
                inner_args = [_convert_template_type(a.strip(), param_map)
                              for a in m.group(2).split(',')]
                cpp_cls = _TEMPLATE_PXD_TO_CPP.get(cls_name, cls_name)
                return f"{cpp_cls}<{', '.join(inner_args)}>"
            # Direct mapping
            bare = result.replace("const ", "").rstrip(" &*").strip()
            if bare in _TEMPLATE_PXD_TO_CPP:
                result = result.replace(bare, _TEMPLATE_PXD_TO_CPP[bare])
            return result

        template_params = merged_class.cpp_class.template_params
        pxd_template_params = getattr(merged_class, 'pxd_template_params', [])

        for instance_name, instance_args in merged_class.template_instances.items():
            # Skip instances that are in SKIP_CLASSES
            if instance_name in SKIP_CLASSES:
                continue

            # Build C++ type args
            cpp_args = []
            for arg in instance_args:
                cpp_type = _TEMPLATE_PXD_TO_CPP.get(arg, arg)
                cpp_args.append(cpp_type)

            cpp_type_str = f"{namespace}::{class_name}<{', '.join(cpp_args)}>"

            # Build template param -> arg mapping for type substitution
            # Use both libclang and pxd template param names (they may differ,
            # e.g., libclang has 'Key' while pxd has 'KeyType')
            param_map = {}
            if len(template_params) == len(instance_args):
                for param, arg in zip(template_params, instance_args):
                    param_map[param] = _TEMPLATE_PXD_TO_CPP.get(arg, arg)
            if pxd_template_params and len(pxd_template_params) == len(instance_args):
                for param, arg in zip(pxd_template_params, instance_args):
                    param_map[param] = _TEMPLATE_PXD_TO_CPP.get(arg, arg)

            lines.append(f'    nb::class_<{cpp_type_str}>(m, "{instance_name}")')

            # Generate default constructor
            lines.append(f"        .def(nb::init<>())")

            # Generate copy constructor
            lines.append(f"        .def(nb::init<const {cpp_type_str}&>())")

            # Generate parameterized constructors
            for ctor in merged_class.constructors:
                params = ctor.parameters
                if not params:
                    continue  # default ctor already emitted
                # Skip copy constructors (including template forms like "MRMTransitionGroup[T1, T2] &")
                if len(params) == 1:
                    bare = params[0].type_str.replace("const ", "").rstrip(" &*").strip()
                    if bare == class_name or bare.endswith(f"::{class_name}") or bare.startswith(f"{class_name}["):
                        continue
                # Substitute template params
                ctor_types = []
                for p in params:
                    p_type = _convert_template_type(p.type_str, param_map)
                    ctor_types.append(p_type)
                type_list = ", ".join(ctor_types)
                lines.append(f"        .def(nb::init<{type_list}>())")

            # Generate methods from the template class
            for merged_method in merged_class.methods:
                method = merged_method.cpp_method
                method_name = method.name

                if merged_method.wrap_ignore:
                    continue
                if method_name in SKIP_METHODS.get(class_name, set()):
                    continue
                if method_name in SKIP_METHODS.get(instance_name, set()):
                    continue
                # Auto-skip methods handled via SPECIAL_METHODS
                if method_name in SPECIAL_METHODS.get(class_name, {}):
                    continue
                if method_name in SPECIAL_METHODS.get(instance_name, {}):
                    continue

                # Substitute template params and convert pxd types to C++
                ret_type = _convert_template_type(method.return_type or "void", param_map)

                params_code = []
                param_names = []
                for p in method.parameters:
                    p_type = _convert_template_type(p.type_str, param_map)
                    params_code.append(f"{p_type} {p.name}")
                    param_names.append(p.name)

                # Build lambda — always use non-const self for template instances
                # because pxd const inference is unreliable (e.g., getTransitionsMuteable
                # is inferred as const because it starts with "get")
                self_param = f"{cpp_type_str}& self"

                if params_code:
                    lambda_params = f"{self_param}, {', '.join(params_code)}"
                    call_args = ', '.join(param_names)
                    call = f"self.{method_name}({call_args})"
                else:
                    lambda_params = self_param
                    call = f"self.{method_name}()"

                if ret_type == "void":
                    body = call
                else:
                    body = f"return {call}"

                named_args = ''.join(f', "{p.name}"_a' for p in method.parameters)

                # Detect rv_policy for reference/pointer returns (same logic as regular methods)
                rv_policy = ""
                explicit_return_type = ""
                if ret_type != "void":
                    if '*' in ret_type and 'shared_ptr' not in ret_type:
                        rv_policy = ", nb::rv_policy::reference_internal"
                    elif '&' in ret_type and '&&' not in ret_type:
                        bare_ret = ret_type.replace('const ', '').replace('&', '').strip()
                        if bare_ret not in _PRIMITIVE_TYPES and bare_ret not in _TYPE_CASTER_TYPES:
                            rv_policy = ", nb::rv_policy::reference_internal"
                            # Lambda needs explicit return type for references
                            norm_ret = ret_type.strip()
                            if '&' not in norm_ret:
                                norm_ret = f"{norm_ret}&"
                            explicit_return_type = f" -> {norm_ret}"

                python_name = merged_method.wrap_as or method_name
                if explicit_return_type:
                    lines.append(
                        f'        .def("{python_name}", []({lambda_params}){explicit_return_type} {{ {body}; }}{named_args}{rv_policy})'
                    )
                else:
                    lines.append(
                        f'        .def("{python_name}", []({lambda_params}) {{ {body}; }}{named_args}{rv_policy})'
                    )

            # Add member variable bindings
            for mv in getattr(merged_class, 'member_variables', []):
                lines.append(f'        .def_rw("{mv.name}", &{cpp_type_str}::{mv.name})')

            # Add __len__ if size() method exists
            has_size = any(m.cpp_method.name == "size" for m in merged_class.methods if not m.wrap_ignore)
            if has_size:
                lines.append(f'        .def("__len__", []({cpp_type_str}& self) {{ return self.size(); }})')

            # Add special methods for this instance
            if instance_name in SPECIAL_METHODS:
                for method_name, method_code in SPECIAL_METHODS[instance_name].items():
                    if method_name != "__post_class__":
                        lines.append(method_code)

            lines.append("        ;")
            lines.append("")

        if not lines:
            return None

        return "\n".join(lines)

    def _generate_constructor(self, ctor: CppMethod, class_name: str = "") -> Optional[str]:
        """Generate constructor binding."""
        # Skip constructors with char* parameters (any variant)
        # These cause issues with nanobind's string conversion
        # The OpenMS::String constructors work via the type caster
        for p in ctor.parameters:
            ptype = p.type_str.strip()
            # Skip any char* constructor - use String versions instead
            if "char *" in ptype or "char*" in ptype:
                return None

        # Skip constructors with raw pointer parameters (except copy ctor pattern)
        # These often reference internal types (like BulkData) that aren't bound
        for p in ctor.parameters:
            ptype = p.type_str.strip()
            # Check for pointer parameters (but allow shared_ptr, unique_ptr, etc.)
            if "*" in ptype and "shared_ptr" not in ptype and "unique_ptr" not in ptype:
                # Skip constructors with raw pointer params (except char* handled above)
                return None

        # Skip constructors with non-const reference parameters (nb::init can't handle them)
        # These need manual lambda wrappers via SPECIAL_METHODS
        for p in ctor.parameters:
            ptype = p.type_str.strip()
            if "&" in ptype and "&&" not in ptype and "const" not in ptype:
                # Non-const ref — check it's not a copy constructor param
                bare = ptype.rstrip(" &").strip()
                if bare != ctor.name and not bare.endswith(f"::{ctor.name}"):
                    return None

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

            # Skip copy constructors for classes with deleted/private copy ctors
            # Detected from AST (has_deleted_copy_constructor) or copy ctor attributes
            if is_copy_ctor:
                merged = getattr(self, '_current_merged_class', None)
                if merged and merged.has_deleted_copy_constructor:
                    return None

        if not ctor.parameters:
            # Skip default constructor if it's deleted
            if hasattr(self, '_current_merged_class') and self._current_merged_class and self._current_merged_class.has_deleted_default_constructor:
                return None
            return ".def(nb::init<>())"

        # Use canonical types from libclang when available
        param_types = []
        for p in ctor.parameters:
            normalized = self._normalize_type(
                p.type_str, canonical_type=getattr(p, 'canonical_type', '')
            )
            # Detect if this is a copy constructor parameter
            # .pxd files may use "ClassName", "ClassName &", or "const ClassName &"
            # For fallback classes, type may be qualified (e.g., "OpenMS::ColumnHeader")
            bare_type = p.type_str.replace("const ", "").rstrip(" &*").strip()
            is_copy_param = (
                ctor.name == p.type_str or
                ctor.name == p.type_str.rstrip(" &") or
                p.type_str.startswith(f"const {ctor.name}") or
                bare_type.endswith(f"::{ctor.name}") or
                bare_type == ctor.name
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

        # Skip specific methods with complex return types (hardcoded list)
        if class_name in SKIP_METHODS and method.name in SKIP_METHODS[class_name]:
            return None

        # Auto-skip methods handled via SPECIAL_METHODS
        if method.name in SPECIAL_METHODS.get(class_name, {}):
            return None

        # Replace output parameter methods with fixed bindings (e.g., getKeys)
        if method.name in OUTPUT_PARAM_METHODS:
            return OUTPUT_PARAM_METHODS[method.name].format(FULL_CLASS=f"const {qualified_name}")

        # Auto-skip methods that use incomplete (forward-declared) types
        if getattr(method, 'uses_incomplete_type', False):
            incomplete = getattr(method, 'incomplete_types', [])
            logger.debug(f"Skipping {class_name}.{method.name} - uses incomplete types: {incomplete}")
            return None

        # Skip methods with char* parameters - use OpenMS::String versions instead
        # Also skip methods with C++ iterator parameters. nanobind has no type
        # caster for accepting iterators (e.g. __gnu_cxx::__normal_iterator) as
        # function arguments — only for *exposing* them via nb::make_iterator.
        # The non-iterator overloads (taking double/index boundaries) are the
        # Python-appropriate API. If a specific iterator overload is needed, add
        # an index-based wrapper to SPECIAL_METHODS, e.g.:
        #   .def("method", [](Class& self, Container& c, size_t left, size_t right) {
        #       return self.method(c, c.begin() + left, c.begin() + right);
        #   })
        for p in method.parameters:
            ptype = p.type_str.strip()
            if "char *" in ptype or "char*" in ptype:
                logger.debug(f"Skipping {class_name}.{method.name} - has char* parameter")
                return None
            if "iterator" in ptype.lower() or "__gnu_cxx" in ptype:
                logger.debug(f"Skipping {class_name}.{method.name} - has iterator parameter: {ptype}")
                return None

        # Auto-skip methods that have const/non-const overloads (detected via libclang)
        # These cause binding issues - prefer const version (handled below)
        if method.name in merged_class.const_overloaded_methods:
            # Will be handled in _generate_regular_method
            pass

        # Handle operators
        if method.name.startswith("operator"):
            py_op = self._operator_map.get(method.name)
            if py_op:
                return self._generate_operator_binding(method, qualified_name, py_op, merged_class)
            return None

        # Handle static methods
        if method.is_static:
            return self._generate_static_method(
                method, qualified_name, method_name, doc=merged_method.doc,
                cpp_name=merged_method.cpp_name
            )

        # Regular methods
        return self._generate_regular_method(merged_method, qualified_name, method_name, merged_class)

    def _build_lambda_params(self, method: CppMethod):
        """Build lambda parameter declarations, call names, and arg annotations for a method."""
        param_decls = []
        param_names_call = []
        param_names_arg = []

        for p in method.parameters:
            ptype = self._normalize_type(
                p.type_str,
                canonical_type=getattr(p, 'canonical_type', ''),
            )
            # Preserve references from the C++ signature:
            # - Non-const references (output parameters) so mutations
            #   are visible to the caller.
            # - Const references for non-copyable types (e.g. QcMLFile).
            #   Always preserving const& is safe and avoids copy issues.
            is_ref = getattr(p, 'is_reference', False) or '&' in p.type_str
            is_const = getattr(p, 'is_const', False) or 'const' in p.type_str
            if is_ref and '&' not in ptype:
                if is_const:
                    ptype = f"const {ptype}&"
                else:
                    ptype = ptype + '&'
            valid_name = (
                p.name
                and not p.name.startswith("arg")
                and p.name not in CPP_KEYWORDS
            )
            pname = p.name if valid_name else f"p{len(param_decls)}"
            param_decls.append(f"{ptype} {pname}")
            param_names_call.append(pname)
            if valid_name:
                param_names_arg.append(f'"{p.name}"_a')

        return (
            ", ".join(param_decls),
            ", ".join(param_names_call),
            param_names_arg,
        )

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

        # For const/non-const overloads with identical parameter types, keep only const version.
        # For overloads with different parameter types/counts, emit ALL of them — nanobind
        # handles dispatch at runtime based on argument types.
        if is_overloaded:
            # Group overloads by parameter signature (count + types)
            def _param_sig(m):
                return tuple(p.type_str for p in m.cpp_method.parameters)

            my_sig = _param_sig(merged_method)
            same_sig_overloads = [m for m in overloads if _param_sig(m) == my_sig]

            if len(same_sig_overloads) == 2:
                const_versions = [m for m in same_sig_overloads if m.cpp_method.is_const]
                non_const_versions = [m for m in same_sig_overloads if not m.cpp_method.is_const]
                if const_versions and non_const_versions:
                    # Check if non-const version returns a mutable reference
                    non_const_ret = non_const_versions[0].cpp_method.return_type or ""
                    is_mutable_ref = ('&' in non_const_ret and 'const' not in non_const_ret
                                      and '&&' not in non_const_ret)
                    if is_mutable_ref:
                        # Prefer non-const (mutable ref return) - skip const version
                        if method.is_const:
                            return None
                    else:
                        # Prefer const version (current behavior) - skip non-const
                        if not method.is_const:
                            return None
            # Otherwise: different parameter signatures — emit all, nanobind dispatches

        # Phase 2: Detect output parameter pattern
        # void f(T& out) or void f(const X& in, T& out) where non-const ref params are outputs
        output_param_result = self._try_generate_output_param_wrapper(
            method, qualified_name, method_name, merged_method.doc, merged_class
        )
        if output_param_result is not None:
            return output_param_result

        # Build lambda parameters
        params_decl, params_call, param_names_arg = self._build_lambda_params(method)

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

        # Phase 3 & 4: Detect return type policies
        ret_type = method.return_type or ""
        canonical_ret = getattr(method, 'canonical_return_type', ret_type)
        rv_policy = ""
        explicit_return_type = ""

        # Phase 3: Pointer return → rv_policy
        # For singleton DB classes (e.g. RibonucleotideDB), returned pointers are
        # globally owned and outlive `self`, so use rv_policy::reference.
        # For regular instance methods, the pointer's lifetime is tied to `self`,
        # so use reference_internal to prevent the parent from being GC'd.
        if '*' in ret_type and 'shared_ptr' not in ret_type:
            if self._is_singleton(merged_class):
                rv_policy = ", nb::rv_policy::reference"
            else:
                rv_policy = ", nb::rv_policy::reference_internal"

        # Phase 4: Reference return → rv_policy::reference_internal with explicit return type
        elif '&' in ret_type and '&&' not in ret_type:
            # Strip the return type to check if it's a primitive
            bare_ret = ret_type.replace('const ', '').replace('&', '').strip()
            canonical_bare = canonical_ret.replace('const ', '').replace('&', '').strip()
            if (bare_ret not in _PRIMITIVE_TYPES and canonical_bare not in _PRIMITIVE_TYPES
                    and bare_ret not in _TYPE_CASTER_TYPES and canonical_bare not in _TYPE_CASTER_TYPES):
                rv_policy = ", nb::rv_policy::reference_internal"
                # Lambda must have explicit return type for references
                norm_ret = self._normalize_type(ret_type, preserve_reference=True, canonical_type=canonical_ret)
                # Ensure the return type has a reference
                if '&' not in norm_ret:
                    if 'const' in ret_type:
                        norm_ret = f"const {norm_ret}&"
                    else:
                        norm_ret = f"{norm_ret}&"
                explicit_return_type = f" -> {norm_ret}"

        # Build def call
        if explicit_return_type:
            result = f'.def("{method_name}", {lambda_sig}{explicit_return_type} {{ return {call_expr}; }}'
        else:
            result = f'.def("{method_name}", {lambda_sig} {{ return {call_expr}; }}'
        if param_names_arg:
            result += f", {', '.join(param_names_arg)}"
        if rv_policy:
            result += rv_policy
        if merged_method.doc:
            doc = self._escape_string(merged_method.doc)
            result += f', "{doc}"'
        result += ")"

        return result

    def _try_generate_output_param_wrapper(
        self,
        method: CppMethod,
        qualified_name: str,
        method_name: str,
        doc: str = "",
        merged_class: Optional[MergedClass] = None,
    ) -> Optional[str]:
        """Try to generate an output parameter wrapper for void methods.

        Detects pattern: void f(T& out) or void f(const X& in, T& out)
        where non-const reference params are output-only.

        Returns the binding code if output param pattern is detected, None otherwise.
        """
        ret_type = (method.return_type or "").strip()
        canonical_ret = getattr(method, 'canonical_return_type', ret_type).strip()

        # Only applies to void methods
        if ret_type != "void" and canonical_ret != "void":
            return None

        # Only auto-wrap methods whose name strongly suggests output semantics.
        if not any(method.name.startswith(p) for p in _OUTPUT_METHOD_PREFIXES):
            return None

        # Find non-const reference parameters (output params)
        output_params = []
        input_params = []
        for p in method.parameters:
            is_ref = getattr(p, 'is_reference', False) or '&' in p.type_str
            is_const = getattr(p, 'is_const', False) or 'const' in p.type_str
            if is_ref and not is_const:
                output_params.append(p)
            else:
                input_params.append(p)

        # Must have at least one output param
        if not output_params:
            return None

        # Conservative: only auto-detect when ALL non-const ref params are outputs
        # and there's at least one such param
        # Skip if too many output params (complex case, leave to SPECIAL_METHODS)
        if len(output_params) > 2:
            return None

        # Check for ambiguous overloads: if other overloads of the same method
        # would produce the same Python signature after output param removal, skip
        if merged_class is not None:
            overloads = [m for m in merged_class.methods if m.cpp_method.name == method.name]
            if len(overloads) > 1:
                # Count how many overloads would have the same input param count
                my_input_count = len(input_params)
                same_input_count = 0
                for m in overloads:
                    other = m.cpp_method
                    other_inputs = sum(
                        1 for p in other.parameters
                        if (getattr(p, 'is_const', False) or 'const' in p.type_str)
                        or not (getattr(p, 'is_reference', False) or '&' in p.type_str)
                    )
                    if other_inputs == my_input_count:
                        same_input_count += 1
                if same_input_count > 1:
                    # Multiple overloads would have same Python signature — skip
                    return None

        # Build the wrapper
        const_qual = "const " if method.is_const else ""
        self_decl = f"{const_qual}{qualified_name}& self"

        # Build input param declarations and call args
        input_decls = []
        input_args = []
        call_args = []
        for p in method.parameters:
            is_ref = getattr(p, 'is_reference', False) or '&' in p.type_str
            is_const = getattr(p, 'is_const', False) or 'const' in p.type_str
            ptype = self._normalize_type(
                p.type_str,
                canonical_type=getattr(p, 'canonical_type', ''),
            )
            valid_name = (
                p.name and not p.name.startswith("arg") and p.name not in CPP_KEYWORDS
            )
            pname = p.name if valid_name else f"p{len(input_decls) + len(call_args)}"

            if is_ref and not is_const:
                # Output param: declare local variable, pass to method
                # Strip reference from type for local variable declaration
                local_type = ptype.replace('&', '').strip()
                call_args.append(("output", pname, local_type))
            else:
                # Input param: include in lambda signature
                if is_ref and '&' not in ptype:
                    ptype = f"const {ptype}&" if is_const else f"{ptype}&"
                input_decls.append(f"{ptype} {pname}")
                if valid_name:
                    input_args.append(f'"{p.name}"_a')
                call_args.append(("input", pname, ptype))

        # Build the lambda body
        local_vars = []
        method_call_args = []
        return_vars = []
        for kind, name, type_str in call_args:
            if kind == "output":
                local_vars.append(f"{type_str} {name}")
                method_call_args.append(name)
                return_vars.append(name)
            else:
                method_call_args.append(name)

        # Build lambda
        if input_decls:
            lambda_params = f"{self_decl}, {', '.join(input_decls)}"
        else:
            lambda_params = self_decl

        body_parts = []
        for lv in local_vars:
            body_parts.append(f"{lv};")
        body_parts.append(f"self.{method.name}({', '.join(method_call_args)});")

        if len(return_vars) == 1:
            body_parts.append(f"return {return_vars[0]};")
        else:
            # Multiple output params → return as tuple
            body_parts.append(f"return std::make_tuple({', '.join(return_vars)});")

        body = " ".join(body_parts)
        result = f'.def("{method_name}", []({lambda_params}) {{ {body} }}'
        if input_args:
            result += f", {', '.join(input_args)}"
        if doc:
            escaped_doc = self._escape_string(doc)
            result += f', "{escaped_doc}"'
        result += ")"

        return result

    def _generate_static_method(
        self,
        method: CppMethod,
        qualified_name: str,
        method_name: str,
        doc: str = "",
        cpp_name: Optional[str] = None,
    ) -> Optional[str]:
        """Generate binding for a static method using lambda wrapper.

        Static methods need lambda wrappers for:
        1. Proper type conversion (OpenMS::String <-> Python str, etc.)
        2. Handling overloads
        3. Consistent interface with regular methods

        If cpp_name is provided (for standalone namespace functions), it is used
        as the full qualified function name instead of qualified_name::method.name.
        """
        # nanobind doesn't support lambdas with more than 8 parameters
        if len(method.parameters) > 8:
            return None

        # Detect output parameter pattern for static methods
        output_result = self._try_generate_static_output_param_wrapper(
            method, qualified_name, method_name, doc
        )
        if output_result is not None:
            return output_result

        # Build lambda parameters
        params_decl, params_call, param_names_arg = self._build_lambda_params(method)

        # Determine the C++ function to call
        # For standalone functions, cpp_name contains the full qualified name
        # e.g., "OpenMS::PEFFFile::isPEFFFile"
        if cpp_name:
            cpp_func = cpp_name
        else:
            cpp_func = f"{qualified_name}::{method.name}"

        # Build the lambda (no 'self' for static methods)
        if params_decl:
            lambda_sig = f"[]({params_decl})"
            call_expr = f"{cpp_func}({params_call})"
        else:
            lambda_sig = "[]()"
            call_expr = f"{cpp_func}()"

        # Detect return type policies for static methods
        ret_type = method.return_type or ""
        rv_policy = ""

        # Static methods have no `self` to tie lifetime to, so always use
        # rv_policy::reference for pointer returns (caller must manage lifetime).
        if '*' in ret_type and 'shared_ptr' not in ret_type:
            rv_policy = ", nb::rv_policy::reference"

        # Build def_static call
        result = f'.def_static("{method_name}", {lambda_sig} {{ return {call_expr}; }}'
        # Only add arg annotations if ALL parameters have valid names
        # (nanobind requires either all or none)
        if param_names_arg and len(param_names_arg) == len(method.parameters):
            result += f", {', '.join(param_names_arg)}"
        if rv_policy:
            result += rv_policy
        if doc:
            escaped_doc = self._escape_string(doc)
            result += f', "{escaped_doc}"'
        result += ")"

        return result

    def _try_generate_static_output_param_wrapper(
        self,
        method: CppMethod,
        qualified_name: str,
        method_name: str,
        doc: str = "",
    ) -> Optional[str]:
        """Try to generate output param wrapper for static void methods."""
        ret_type = (method.return_type or "").strip()
        canonical_ret = getattr(method, 'canonical_return_type', ret_type).strip()

        if ret_type != "void" and canonical_ret != "void":
            return None

        # Only auto-wrap methods whose name strongly suggests output semantics.
        if not any(method.name.startswith(p) for p in _OUTPUT_METHOD_PREFIXES):
            return None

        output_params = []
        input_params = []
        for p in method.parameters:
            is_ref = getattr(p, 'is_reference', False) or '&' in p.type_str
            is_const = getattr(p, 'is_const', False) or 'const' in p.type_str
            if is_ref and not is_const:
                output_params.append(p)
            else:
                input_params.append(p)

        if not output_params or len(output_params) > 2:
            return None

        input_decls = []
        input_args = []
        call_args = []
        for p in method.parameters:
            is_ref = getattr(p, 'is_reference', False) or '&' in p.type_str
            is_const = getattr(p, 'is_const', False) or 'const' in p.type_str
            ptype = self._normalize_type(
                p.type_str,
                canonical_type=getattr(p, 'canonical_type', ''),
            )
            valid_name = (
                p.name and not p.name.startswith("arg") and p.name not in CPP_KEYWORDS
            )
            pname = p.name if valid_name else f"p{len(input_decls) + len(call_args)}"

            if is_ref and not is_const:
                local_type = ptype.replace('&', '').strip()
                call_args.append(("output", pname, local_type))
            else:
                if is_ref and '&' not in ptype:
                    ptype = f"const {ptype}&" if is_const else f"{ptype}&"
                input_decls.append(f"{ptype} {pname}")
                if valid_name:
                    input_args.append(f'"{p.name}"_a')
                call_args.append(("input", pname, ptype))

        local_vars = []
        method_call_args = []
        return_vars = []
        for kind, name, type_str in call_args:
            if kind == "output":
                local_vars.append(f"{type_str} {name}")
                method_call_args.append(name)
                return_vars.append(name)
            else:
                method_call_args.append(name)

        if input_decls:
            lambda_params = ', '.join(input_decls)
        else:
            lambda_params = ""

        body_parts = []
        for lv in local_vars:
            body_parts.append(f"{lv};")
        body_parts.append(f"{qualified_name}::{method.name}({', '.join(method_call_args)});")

        if len(return_vars) == 1:
            body_parts.append(f"return {return_vars[0]};")
        else:
            body_parts.append(f"return std::make_tuple({', '.join(return_vars)});")

        body = " ".join(body_parts)
        result = f'.def_static("{method_name}", []({lambda_params}) {{ {body} }}'
        if input_args:
            result += f", {', '.join(input_args)}"
        if doc:
            escaped_doc = self._escape_string(doc)
            result += f', "{escaped_doc}"'
        result += ")"

        return result

    def _generate_operator_binding(
        self, method: CppMethod, qualified_name: str, py_op: str,
        merged_class: MergedClass = None,
    ) -> str:
        """Generate binding for an operator."""
        op_name = method.name.replace("operator", "")

        if op_name == "[]":
            # Add bounds checking if class has size() method
            has_size = merged_class and (
                self._has_size_method(merged_class)
                or (merged_class.name in CONTAINER_CLASSES)
            )
            # Detect if operator[] returns a non-primitive reference
            ret_type = method.return_type or ""
            needs_ref_policy = (
                '&' in ret_type and '&&' not in ret_type
                and 'const' not in ret_type
            )
            if needs_ref_policy:
                bare_ret = ret_type.replace('&', '').strip()
                if bare_ret in _PRIMITIVE_TYPES or bare_ret in _TYPE_CASTER_TYPES:
                    needs_ref_policy = False
            rv_suffix = ", nb::rv_policy::reference_internal" if needs_ref_policy else ""
            ret_arrow = f" -> {ret_type}" if needs_ref_policy else ""
            if has_size:
                return (
                    f'.def("__getitem__", []({qualified_name}& self, size_t i){ret_arrow} {{ '
                    f'if (i >= self.size()) throw nb::index_error(); '
                    f'return self[i]; }}{rv_suffix})'
                )
            return f'.def("__getitem__", []({qualified_name}& self, size_t i){ret_arrow} {{ return self[i]; }}{rv_suffix})'
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
        result = re.sub(r'libcpp_vector\[([^\]]+)\]', r'std::vector<\1>', result)
        result = re.sub(r'libcpp_map\[([^,]+),\s*([^\]]+)\]', r'std::map<\1, \2>', result)
        result = re.sub(r'libcpp_set\[([^\]]+)\]', r'std::set<\1>', result)
        result = re.sub(r'libcpp_pair\[([^,]+),\s*([^\]]+)\]', r'std::pair<\1, \2>', result)
        result = re.sub(r'shared_ptr\[([^\]]+)\]', r'std::shared_ptr<\1>', result)

        # libclang resolves libcpp_utf8_string to std::basic_string<char>;
        # map it to OpenMS::String so nanobind uses our custom caster (accepts bytes & str).
        # Only replace when it IS the type (possibly with const/ref), not inside containers.
        result = re.sub(r'^(const\s+)?std::basic_string<char>(\s*[&*])?$',
                         lambda m: (m.group(1) or '') + 'OpenMS::String' + (m.group(2) or ''),
                         result)

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
            'LightTargetedExperiment': 'OpenSwath::LightTargetedExperiment',
            'LightModification': 'OpenSwath::LightModification',
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
            # FLASHHelperClasses types used in constructors/parameters
            'CoarseIsotopePatternGenerator': 'OpenMS::CoarseIsotopePatternGenerator',
            'IsotopeDistribution': 'OpenMS::IsotopeDistribution',
        }

        # OpenMS classes that need namespace (only add if not already namespaced)
        for typedef, qualified in openms_typedefs.items():
            # Match the type when not already qualified
            pattern = r'\b(?<!::)(?<!OpenMS::)(?<!std::)' + typedef + r'\b'
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

        # Strip reference qualifiers unless we're preserving them for static_cast
        if not preserve_reference:
            result = result.replace("const ", "").replace(" const", "")
            result = result.replace("&", "").strip()

        return result

    def _is_qt_class(self, class_name: str) -> bool:
        """Check if a class name is a Qt class (QDate, QString, QObject, etc.).

        Qt classes follow the pattern: Q + CapitalLetter (e.g., QDate, QString, QObject)
        But NOT:
        - QC* (QCBase = Quality Control in OpenMS)
        - QT* (QTCluster = OpenMS class)
        - Q followed by lowercase (not a Qt convention)
        """
        if not class_name.startswith('Q'):
            return False
        if len(class_name) < 2:
            return False
        second_char = class_name[1]
        # Exclude QC (Quality Control) and QT (OpenMS) prefixes
        if second_char in ('C', 'T'):
            return False
        # Qt classes have Q + uppercase letter
        return second_char.isupper()

    def _escape_string(self, s: str) -> str:
        """Escape a string for use in C++ code."""
        return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")

    def _qualify_openms_types(self, type_str: str) -> str:
        """Add OpenMS:: prefix to unqualified OpenMS types.

        Handles both the main type and template parameters, e.g.:
        RangeManagerContainer< RangeMZ, RangeIntensity >
        -> OpenMS::RangeManagerContainer<OpenMS::RangeMZ, OpenMS::RangeIntensity>
        """
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
        self, output_path: Path, domain_name: str, content: ModuleContent
    ) -> None:
        """Write a domain module C++ file as a standalone NB_MODULE."""
        lines = []

        # Header comment
        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append(f"// Domain: {domain_name}")
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
        lines.append(f'NB_MODULE(_pyopenms2_{domain_name}, m) {{')
        lines.append(f'    m.doc() = "pyOpenMS2 {domain_name} bindings";')
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

    def _write_single_module(
        self, output_path: Path, content: ModuleContent
    ) -> None:
        """Write all bindings into a single main_module.cpp (no sub-modules)."""
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append("// Single-module build")
        lines.append("")

        for inc in sorted(content.includes):
            lines.append(f"#include {inc}")
        lines.append("")

        lines.append("namespace nb = nanobind;")
        lines.append("using namespace nb::literals;")
        lines.append("")

        lines.append('NB_MODULE(_pyopenms2, m) {')
        lines.append('    m.doc() = "pyOpenMS2: Python bindings for OpenMS (nanobind)";')
        lines.append("")

        for binding in content.enum_bindings:
            lines.append(binding)
            lines.append("")

        for binding in content.class_bindings:
            lines.append(binding)
            lines.append("")

        for binding in content.post_class_enums:
            lines.append(binding)
            lines.append("")

        lines.append("}")

        output_path.write_text("\n".join(lines))

    def _write_main_module(self, output_path: Path) -> None:
        """Write the main module file.

        Note: With standalone NB_MODULE per sub-module, this main module is
        optional. It's kept for backward compatibility but just creates an
        empty module - the actual classes are in _pyopenms2_<domain> modules.
        """
        lines = []

        lines.append("// Auto-generated by pyOpenMS2 generator (v2 - libclang)")
        lines.append("// Main module - placeholder, actual bindings in _pyopenms2_<domain> modules")
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
