// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// TypeScript definitions of the JavaScript values of js_types.h, for openms.d.ts.

#include "js_types.h"

#include <emscripten/bind.h>

using namespace emscripten;
using namespace OpenMS::WASM::JS;

EMSCRIPTEN_BINDINGS(openms_js_types)
{
  // generic
  register_type<NumberArray>("number[]");
  register_type<NumberInput>("ArrayLike<number>");
  register_type<StringArray>("string[]");
  register_type<MetaValue>("MetaValue", "number | string | number[] | string[] | null | undefined");
  register_type<MetaValueInput>("number | string | boolean | number[] | string[]");
  register_type<MetaValues>("MetaValues", "Record<string, number | string | number[] | string[] | null>");

  // parameters
  register_type<ParamValue>("ParamValue", "number | string | boolean | number[] | string[] | null");
  register_type<ParamValues>("ParamValues", "Record<string, ParamValue>");
  register_type<ParamEntry>("ParamEntry",
    "{ name: string; value: ParamValue; type: 'string' | 'int' | 'double' | 'bool' | 'string[]' | 'int[]' | 'double[]' | 'empty'; "
    "description: string; tags: string[]; validStrings?: string[]; min?: number; max?: number }");
  register_type<ParamEntries>("ParamEntry[]");

  // spectra and chromatograms
  register_type<SpectrumPeaks>("SpectrumPeaks", "{ mz: Float64Array; intensity: Float32Array }");
  register_type<ChromatogramPeaks>("ChromatogramPeaks", "{ rt: Float64Array; intensity: Float32Array }");
  register_type<PeakObject>("Peak", "{ mz: number; intensity: number }");
  register_type<OptionalPeak>("Peak | null");
  register_type<BasePeak>("{ index: number; mz: number; intensity: number } | null");
  register_type<PrecursorObject>("Precursor",
    "{ mz: number; intensity: number; charge: number; isolationWindowLowerOffset: number; isolationWindowUpperOffset: number; "
    "activationEnergy: number; activationMethods: string[]; driftTime: number }");
  register_type<PrecursorInput>("Partial<Precursor>");
  register_type<PrecursorList>("Precursor[]");
  register_type<PrecursorInputList>("Partial<Precursor>[]");
  register_type<ProductObject>("Product", "{ mz: number; isolationWindowLowerOffset: number; isolationWindowUpperOffset: number }");
  register_type<ProductList>("Product[]");
  register_type<FloatDataArrays>("{ name: string; data: Float32Array }[]");
  register_type<IntegerDataArrays>("{ name: string; data: Int32Array }[]");
  register_type<StringDataArrays>("{ name: string; data: string[] }[]");
  register_type<SpectrumSummary>(
    "{ index: number; nativeID: string; rt: number; msLevel: number; size: number; tic: number; basePeakMZ?: number; "
    "basePeakIntensity?: number; minMZ?: number; maxMZ?: number; precursorMZ?: number; precursorCharge?: number }[]");
  register_type<PeakData2D>("PeakData2D", "{ rt: Float64Array; mz: Float64Array; intensity: Float32Array }");

  // chemistry
  register_type<ResidueObject>("Residue",
    "{ name: string; oneLetterCode: string; threeLetterCode: string; formula: string; monoWeight: number; "
    "averageWeight: number; isModified: boolean; modificationName: string | null }");
  register_type<ResidueList>("Residue[]");
  register_type<ModificationObject>("Modification",
    "{ id: string; fullId: string; fullName: string; unimodAccession: string; psiModAccession: string; origin: string; "
    "termSpecificity: string; diffMonoMass: number; monoMass: number; diffFormula: string } | null");
  register_type<TermSpecificityInput>("'any' | 'C-term' | 'N-term' | 'Protein C-term' | 'Protein N-term' | 'none' | null | undefined");
  register_type<EnzymeObject>("Enzyme", "{ name: string; regEx: string; regExDescription: string; psiId: string }");
  register_type<ElementCounts>("Record<string, number>");
  register_type<IsotopePeaks>("IsotopePeaks", "{ mz: Float64Array; intensity: Float64Array }");
  register_type<PeakObjectList>("Peak[]");
  register_type<AASequenceList>("AASequence[]");

  // processing
  register_type<Alignment>("[number, number][]");
  register_type<ThreadingInfo>("ThreadingInfo", "{ threaded: boolean; available: boolean; numThreads: number }");
  register_type<MemoryInfo>("MemoryInfo", "{ heapSize: number; heapMax: number }");

  // files
  register_type<PeakFileOptions>("PeakFileOptions",
    "{ msLevels?: number[]; rtRange?: [number, number]; mzRange?: [number, number]; intensityRange?: [number, number]; "
    "metadataOnly?: boolean; fillData?: boolean; sortSpectraByMZ?: boolean; skipXMLChecks?: boolean }");
  register_type<FASTAEntries>("{ identifier: string; description?: string; sequence: string }[]");

  // identifications
  register_type<PeptideHitObject>("PeptideHitData",
    "{ sequence: string; score: number; charge: number; rank: number; isDecoy: boolean; proteinAccessions: string[]; metaValues: MetaValues }");
  register_type<PeptideIdentificationObject>("PeptideIdentificationData",
    "{ rt: number | null; mz: number | null; scoreType: string; higherScoreBetter: boolean; identifier: string; "
    "spectrumReference: string; hits: PeptideHitData[]; metaValues: MetaValues }");
  register_type<PeptideIdentificationObjects>("PeptideIdentificationData[]");
  register_type<ProteinHitObject>("ProteinHitData",
    "{ accession: string; score: number; description: string; sequence: string; coverage: number; metaValues: MetaValues }");
  register_type<ProteinIdentificationObject>("ProteinIdentificationData",
    "{ identifier: string; searchEngine: string; searchEngineVersion: string; scoreType: string; higherScoreBetter: boolean; "
    "hits: ProteinHitData[]; metaValues: MetaValues }");
  register_type<ProteinIdentificationObjects>("ProteinIdentificationData[]");
  register_type<PeptideHitList>("PeptideHit[]");
  register_type<ProteinHitList>("ProteinHit[]");
  register_type<PeakAnnotations>("{ mz: number; intensity: number; annotation: string; charge: number }[]");

  // ProSE
  register_type<ProgressCallback>(
    "((progress: { event: 'start' | 'progress' | 'end'; label: string; depth: number; begin: number; end: number; value: number }) => void) "
    "| null | undefined");
  register_type<ProSESummary>("ProSESummary",
    "{ status: 'EXECUTION_OK' | 'INPUT_FILE_EMPTY' | 'UNEXPECTED_RESULT' | 'ILLEGAL_PARAMETERS' | 'UNKNOWN_ERROR'; "
    "openSearch: boolean; haveDecoys: boolean; decoyString: string; decoyIsPrefix: boolean; proteinFDRApplied: boolean; "
    "percolator: { applied: boolean; skipped: string | null; seconds: number } | null; "
    "run: { inputFile: string; ms2Spectra: number; matchedSpectra: number; targetPSMs: number; decoyPSMs: number; "
    "fdrApplied: boolean; achievedPsmFDR: number | null; uniquePeptides: number; uniqueProteins: number; "
    "chargeHistogram: Record<string, number>; missedCleavageHistogram: Record<string, number>; "
    "hyperscore?: { min: number; median: number; max: number }; "
    "precursorError: { median: number; mad: number; recommended: number } | null; "
    "fragmentError: { median: number; mad: number; recommended: number } | null; "
    "secondsSearch: number; secondsCalibration: number; secondsFDR: number }; "
    "search: { databaseFile: string; enzyme: string; openSearch: boolean; calibration: boolean; chunked: boolean; "
    "decoyMode: string; targetProteins: number; decoyProteins: number; indexedPeptides: number; indexedFragments: number; "
    "secondsIndexBuild: number; secondsTotal: number }; "
    "modificationAnalysis: { "
    "deltaMasses: { totalPSMs: number; modifiedPSMs: number; unmodifiedPSMs: number; meanDeltaMass: number; medianDeltaMass: number; "
    "entries: { deltaMass: number; count: number; uniquePeptides: number; chargeStates: number; percentage: number; "
    "modification: string; known: boolean }[] }; "
    "ptms: { totalModifiedPSMs: number; unknownModificationPSMs: number; uniqueModifications: number; "
    "entries: { name: string; theoreticalMass: number; observedMass: number; massDeviation: number; count: number; "
    "uniquePeptides: number; chargeStates: number; percentage: number; targetResidues: string; residueCounts: Record<string, number> }[] } "
    "} | null }");
  register_type<ProSEOptions>("ProSEOptions", "{ /** Percolator rescoring (default: true) */ rescore?: boolean } | null | undefined");
  register_type<AnnotationOptions>("AnnotationOptions",
    "{ tolerance?: number; unit?: 'Da' | 'ppm'; maxCharge?: number; ions?: ('a' | 'b' | 'c' | 'x' | 'y' | 'z')[]; "
    "losses?: boolean; precursor?: boolean } | null | undefined");
  register_type<SpectrumAnnotation>("SpectrumAnnotation",
    "{ matches: { mz: number; intensity: number; ion: string; charge: number; theoreticalMZ: number; error: number }[]; "
    "theoretical: { mz: Float64Array; ion: string[]; charge: number[] }; unit: 'Da' | 'ppm' }");

  // features
  register_type<BoundingBox>("BoundingBox", "{ minRT: number; maxRT: number; minMZ: number; maxMZ: number }");
  register_type<OptionalBoundingBox>("BoundingBox | null");
  register_type<BoundingBoxList>("BoundingBox[]");
  register_type<HullPoints>("[number, number][]");
  register_type<FeatureObject>("FeatureData",
    "{ uniqueId: string; rt: number; mz: number; intensity: number; charge: number; quality: number; width: number; "
    "boundingBox: BoundingBox | null; metaValues: MetaValues }");
  register_type<FeatureObjects>("FeatureData[]");
  register_type<FeatureTable>("FeatureTable",
    "{ rt: Float64Array; mz: Float64Array; intensity: Float32Array; charge: Int32Array; quality: Float32Array; width: Float32Array }");
  register_type<ConsensusFeatureObject>("ConsensusFeatureData",
    "{ uniqueId: string; rt: number; mz: number; intensity: number; charge: number; quality: number; "
    "elements: { mapIndex: number; uniqueId: string; rt: number; mz: number; intensity: number; charge: number }[]; metaValues: MetaValues }");
  register_type<ConsensusFeatureObjects>("ConsensusFeatureData[]");
  register_type<ColumnHeaders>("{ index: number; filename: string; label: string; size: number }[]");
  register_type<MassTraceObjects>(
    "{ label: string; mz: number; rt: number; mzSD: number; intensity: number; maxIntensity: number; fwhm: number; traceLength: number; size: number }[]");
}
