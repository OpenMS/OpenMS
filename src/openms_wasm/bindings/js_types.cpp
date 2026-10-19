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
  register_type<Bytes>("Uint8Array");

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
    "((progress: { event: 'start' | 'progress' | 'end'; label: string; depth: number; begin: number; end: number; value: number; "
    "/** position of the input file being searched (its step 'Searching <name> (<i> of <n>)...': i - 1), -1 outside a file */ "
    "file: number }) => void) | null | undefined");
  register_type<ProSESummary>("ProSESummary",
    "{ status: 'EXECUTION_OK' | 'INPUT_FILE_EMPTY' | 'UNEXPECTED_RESULT' | 'ILLEGAL_PARAMETERS' | 'UNKNOWN_ERROR'; "
    "openSearch: boolean; haveDecoys: boolean; decoyString: string; decoyIsPrefix: boolean; proteinFDRApplied: boolean; "
    "percolator: { applied: boolean; skipped: string | null; "
    "/** a caveat of an applied rescoring (Percolator's sanity check replaced the learned model by the best single feature); null: none */ "
    "warning: string | null } | null; "
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
  register_type<ProSEOptions>("ProSEOptions",
    "{ /** Percolator rescoring: sets the parameter 'rescore' for the call (default: the parameter, true) */ rescore?: boolean } "
    "| null | undefined");
  register_type<ProSEFilesOptions>("ProSEFilesOptions",
    "{ /** Percolator rescoring of each file: sets the parameter 'rescore' for the call (default: the parameter, true) */ "
    "rescore?: boolean; "
    "/** keep the peaks of the identified spectra with the results (default true) */ keepSpectra?: boolean; "
    "/** with the protein FDR, remove the decoy PSMs and proteins from the pooled results, as the ProSE tool's -out_merged "
    "does (default false: keep those that pass) */ removeDecoys?: boolean; "
    "/** search the files one at a time: called before the search of file i (its position in files) to write it at path, "
    "e.g. with allocateFile(path, size, capacity) and writeFileChunk(); only the largest file then has to fit into the "
    "module's memory. Give allocateFile() the size of the largest file as capacity: the files, which come one after the "
    "other, then reuse one block instead of growing the heap with each larger file. Has to be synchronous */ "
    "provideFile?: ((i: number, path: string) => void) | null; "
    "/** with provideFile: called after the search of file i (and after reading the spectra of its PSMs), to remove it */ "
    "releaseFile?: ((i: number, path: string) => void) | null } | null | undefined");
  register_type<ProSEFileSummary>("ProSEFileSummary",
    "{ path: string; status: 'EXECUTION_OK' | 'INPUT_FILE_EMPTY' | 'UNEXPECTED_RESULT' | 'ILLEGAL_PARAMETERS' | 'UNKNOWN_ERROR'; "
    "ms2Spectra: number; targetPSMs: number; decoyPSMs: number; "
    "percolator: { applied: boolean; skipped: string | null; warning: string | null } | null; "
    "peptdeepNce: number | null; "
    "/** the precursor mass error of the file's PSMs and the tolerance ProSE recommends for it, as in search(); null: not estimated */ "
    "precursorError: ProSESummary['run']['precursorError']; "
    "/** the fragment mass error and the tolerance ProSE recommends, as precursorError */ "
    "fragmentError: ProSESummary['run']['fragmentError']; secondsSearch: number }");
  register_type<ProSEMultiSummary>("ProSEMultiSummary",
    "{ files: ProSEFileSummary[]; "
    "/** the database, the fragment index (indexedPeptides, indexedFragments), the decoys (decoyMode) and the times, as "
    "search() reports them (summary.search). secondsTotal is the time of searchFiles() (with provideFile, the writing and removing "
    "of the files too); secondsIndexBuild that of building the index (with provideFile, of all the files together) */ "
    "search: ProSESummary['search']; proteinFdrApplied: boolean; psmFdr: number; proteinFdr: number; decoyString: string; "
    "openSearch: boolean; modificationAnalysis: ProSESummary['modificationAnalysis'] }");
  register_type<ResultTables>("ResultTables",
    "{ /** the input files, in the order of the search */ files: string[]; "
    "/** one row per identified spectrum (its best hit), in the order of the pooled PSMs */ "
    "psms: { count: number; /** position of the PSM's file in files */ file: Uint16Array; "
    "/** native ID of the spectrum */ nativeId: string[]; /** NaN: unknown */ rt: Float64Array; "
    "/** precursor m/z; NaN: unknown */ mz: Float64Array; charge: Int8Array; /** row of peptides */ peptide: Uint32Array; "
    "/** Percolator's (SVM) score if the search rescored the PSMs, else ProSE's ln(hyperscore); higher is better. A hit "
    "of a rescored file that Percolator could not rescore (e.g. without protein reference) has -100, with q-value and "
    "PEP 1, as Percolator's output gives it */ score: Float64Array; "
    "/** Percolator's q-value if the file was rescored, else the q-value of target-decoy competition of the file's PSMs, "
    "as FDR:PSM computes it (with FDR:PSM 0, which filters nothing, too). NaN: none, because the search had no decoys "
    "(or the file no decoy PSM) */ qValue: Float64Array; "
    "/** posterior error probability (Percolator); NaN: not rescored */ "
    "pep: Float64Array; /** NaN: unknown */ precursorErrorPpm: Float64Array; "
    "/** precursor isotope error, observed minus theoretical (in 13C spacings); 0: unknown */ isotopeError: Int8Array; "
    "/** mass shift of an open search (Da), not localized: peptides.proforma never has it (write '[+delta]?' + proforma, "
    "delta being the signed psms.deltaMass, e.g. '[-57.0215]?'); NaN: none */ deltaMass: Float64Array; "
    "/** 1: decoy */ decoy: Uint8Array; "
    "/** Percolator's input features, as ProSE's rescoring gives them: 'ln(hyperscore)' (its 'score'), the standard "
    "features of the .pin format that Percolator trains on (CalcMass, mass, peplen, charge<N> for N over the precursor "
    "charges of the search, enzN, enzC, enzInt, dm, absdm; computed by PercolatorInfile, also for searches without "
    "rescoring), the search engine's features (the 'extra_features' of the runs) and PeptDeep's (ms2_cosine, "
    "ms2_spectral_angle, ms2_pearson, ms2_frac_pred_found, rt_abs_error); of these, those that some PSM has */ "
    "featureNames: string[]; /** row-major, count * featureNames.length: features[i * featureNames.length + j]; NaN: absent */ "
    "features: Float64Array }; "
    "/** unique by openms */ "
    "peptides: { count: number; /** AASequence.toString(), e.g. 'AVM(Oxidation)DDFAAFVEK', '.(Acetyl)DAHKSEVAHR' */ openms: string[]; "
    "/** ProForma with Unimod names, fixed modifications included: 'AVM[Oxidation]DDFAAFVEK', '[Acetyl]-DAHKSEVAHR', "
    "'PEPTIDER-[Amidated]'; other modifications as mass deltas, e.g. '[+15.9949]'. Never the unlocalized mass shift of an "
    "open search, which is a value of each PSM (psms.deltaMass) */ proforma: string[]; "
    "/** the residues */ sequence: string[]; /** 1: decoy */ decoy: Uint8Array; "
    "/** the modifications of peptide k: modOffsets[k] .. modOffsets[k + 1] of the mod* columns (count + 1 entries) */ "
    "modOffsets: Uint32Array; /** -1: N-terminus, sequence.length: C-terminus, else the 0-based residue */ modPosition: Int16Array; "
    "/** the (Unimod) name, e.g. 'Oxidation' */ modName: string[]; /** Unimod accession number, e.g. 35; -1: none */ "
    "modUnimod: Int32Array; /** monoisotopic mass delta (Da) */ modMass: Float64Array; "
    "/** the proteins of peptide k: proteinOffsets[k] .. proteinOffsets[k + 1] of protein, start and end (count + 1 entries); "
    "one per occurrence. proteinOffsets[k] === proteinOffsets[k + 1]: the peptide has no protein, because all its proteins "
    "failed the protein FDR (the references to them are removed, the peptide and its PSMs stay: 'below protein FDR') */ "
    "proteinOffsets: Uint32Array; /** row of proteins */ protein: Uint32Array; "
    "/** the peptide in proteins.sequence: 0-based start, exclusive end. Where the sequence has the peptide; if the "
    "sequence is unknown ('') or lacks the peptide, the position the search recorded, which may exceed sequence.length; "
    "0 and 0 if none is recorded */ start: Uint32Array; end: Uint32Array }; "
    "/** the protein hits of the pooled run */ "
    "proteins: { count: number; accession: string[]; /** from the FASTA */ description: string[]; "
    "/** as searched: from the FASTA; a decoy that ProSE generated: as generated; '' if unknown */ sequence: string[]; "
    "/** the best score (psms.score) of its PSMs; NaN: none */ score: Float64Array; "
    "/** picked-protein FDR; NaN: no protein FDR */ qValue: Float64Array; /** 1: decoy */ decoy: Uint8Array; "
    "/** index of its indistinguishable-protein group; -1: none */ group: Int32Array } }");
  register_type<StoredSpectra>("StoredSpectra",
    "{ /** per PSM row of ResultTables.psms: the index of its spectrum; -1: none (keepSpectra false, a q-value above "
    "FDR:PSM, no MS2 spectrum with its native ID in its file, or several of them and not exactly one at the PSM's RT, "
    "or a file that could not be read again). The PSMs of the same spectrum share it */ "
    "psmSpectrum: Int32Array; "
    "/** the peaks of spectrum s: offsets[s] .. offsets[s + 1] of mz and intensity (n + 1 entries) */ offsets: Uint32Array; "
    "/** the peaks of the MS2 spectra as read from the files (not preprocessed), sorted by m/z */ mz: Float32Array; "
    "intensity: Float32Array; /** m/z of the first precursor; NaN: none */ precursorMz: Float64Array; "
    "/** collision energy of the first precursor in eV (MS:1000045 \"collision energy\" of the file, else its MS:1000509 \"activation energy\"); NaN: unknown */ collisionEnergy: Float32Array }");
  register_type<AnnotationOptions>("AnnotationOptions",
    "{ tolerance?: number; unit?: 'Da' | 'ppm'; maxCharge?: number; ions?: ('a' | 'b' | 'c' | 'x' | 'y' | 'z')[]; "
    "losses?: boolean; precursor?: boolean } | null | undefined");
  register_type<SpectrumAnnotation>("SpectrumAnnotation",
    "{ matches: { mz: number; intensity: number; ion: string; charge: number; theoreticalMZ: number; error: number }[]; "
    "theoretical: { mz: Float64Array; ion: string[]; charge: number[] }; unit: 'Da' | 'ppm' }");
  register_type<PeakMatch>("PeakMatch",
    "{ /** index of the peak in the mz and intensity arrays as passed (not in sorted order); in ppm several matches can "
    "have the same peak */ peak: number; "
    "/** the ion as labelled in a spectrum: 'b3'; 'y5++' (one '+' per charge, none for charge 1); 'y7-H2O' with a neutral "
    "loss; '[M+2H]++' and '[M+H-NH3]' for precursor peaks */ ion: string; "
    "/** 'a' to 'z', or 'precursor' */ type: string; /** the ion number (3 for b3); 0 for a precursor */ ordinal: number; "
    "charge: number; /** formula of the neutral loss ('H2O', 'NH3', 'H3PO4'); '' if none */ loss: string; "
    "theoreticalMz: number; /** (observed - theoretical) / theoretical * 1e6, whatever the unit of the tolerance */ "
    "errorPpm: number }");
  register_type<PeakAnnotationResult>("PeakAnnotationResult",
    "{ /** the ions that have a peak within the tolerance, by ascending m/z: one match per ion. In ppm each ion takes the "
    "closest peak, so several ions can have the same peak (the same index in peak), and of peaks of equal m/z the first in "
    "the arrays. In Da the matching is one to one and in m/z order: a peak has at most one ion (of two ions of one m/z, "
    "one is left out, unless there are two peaks), and of peaks of equal m/z one is taken; which one is unspecified */ "
    "matches: PeakMatch[]; "
    "/** every ion of the peptide for the options (charges 1 to maxCharge), by ascending m/z, matched or not */ "
    "theoretical: { ion: string; type: string; ordinal: number; charge: number; mz: number }[] }");
  register_type<PredictedIon>("PredictedIon",
    "{ /** the ion as annotatePeaks names it: 'b3'; 'y5++' (one '+' per charge, none for charge 1) */ ion: string; "
    "type: 'b' | 'y'; /** the ion number, 1 to the length of the peptide - 1 (3 for b3) */ ordinal: number; "
    "/** 2 only if the precursor has charge 2 or more */ charge: 1 | 2; "
    "/** from the theoretical spectrum, with the modifications of the peptide as given; the ions of a peptide come by "
    "ascending m/z */ mz: number; "
    "/** relative to the most intense returned b/y ion, which is 1; ions of 0.01 or less are left out. (For a precursor "
    "of charge 1 the model's own maximum may be in a channel that is not returned.) */ intensity: number }");
  register_type<PredictedIons>("PredictedIon[]");
  register_type<ModelPathInput>("string | null | undefined");

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
