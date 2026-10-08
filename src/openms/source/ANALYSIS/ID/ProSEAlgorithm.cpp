// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <limits>
#include <numeric>
#include <iomanip>
#include <iterator>
#include <locale>
#ifdef __GLIBC__
#include <malloc.h> // malloc_trim
#endif
#include <ostream>
#include <map>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  ProSEAlgorithm::ProSEAlgorithm() :
    DefaultParamHandler("ProSEAlgorithm"),
    ProgressLogger()
  {
    defaults_.setValue("precursor:mass_tolerance_lower", 10.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 10.0,
                       "Upper-side precursor-mass tolerance (positive magnitude).");
    defaults_.setMinFloat("precursor:mass_tolerance_upper", 0.0);
    defaults_.setValue("precursor:mass_tolerance_unit", "ppm", "Unit of precursor mass tolerance.");
    defaults_.setValidStrings("precursor:mass_tolerance_unit", {"ppm", "Da"});

    defaults_.setValue("precursor:min_charge", 2, "Minimum precursor charge to be considered.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered.");

    defaults_.setSectionDescription("precursor",
      "Precursor (Parent Ion) Options. mass_tolerance_lower/_upper are positive magnitudes "
      "applied as [-lower, +upper] around the precursor mass.");

    defaults_.setValue("fragment:mass_tolerance", 20.0, "Fragment mass tolerance");

    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    defaults_.setValue("fragment:deisotope", "auto", "MS2 deisotoping (single-charge deconvolution) before searching. 'auto' deisotopes only when the fragment tolerance is within the high-resolution deisotoper range (<= 0.1 Da / <= 100 ppm) and skips it for low-resolution (e.g. ion-trap CID) data; 'true' always deisotopes (requires a high-resolution fragment tolerance); 'false' never deisotopes.");
    defaults_.setValidStrings("fragment:deisotope", {"auto", "true", "false"});
    defaults_.setValue("fragment:deisotope_min_peaks", 2,
                       "Minimum number of peaks, the monoisotopic peak included, of an isotope envelope that MS2 deisotoping "
                       "collapses into its monoisotopic peak. With 2 (Sage-like), an ion whose M+2 peak is lost in the noise loses its "
                       "M+1 peak too, so the M+1 peak no longer takes a place in the per-window peak quota, and a multiply charged ion "
                       "seen with two isotope peaks is converted to charge 1; in dense spectra (e.g. timsTOF) more of the two-peak "
                       "multiply charged envelopes are chance pairs. 3 is the earlier behaviour.",
                       {"advanced"});
    defaults_.setMinInt("fragment:deisotope_min_peaks", 2);
    defaults_.setMaxInt("fragment:deisotope_min_peaks", 10);
    defaults_.setValue("fragment:deisotope_charge_cap", "precursor",
                       "Fragment charges that MS2 deisotoping tries. 'precursor' tries charges up to the precursor charge, but at most 3 "
                       "(1-3 when the precursor charge is unknown): a fragment cannot carry more charges than its precursor, so a chance "
                       "peak 1/3 Th above an ion of a 2+ precursor cannot turn the ion into a '3+' envelope at a wrong m/z. This is "
                       "Sage-like; Sage caps at the precursor charge without the limit of 3 and counts an unknown charge as 3. 'none' "
                       "tries charges 1-3 in every spectrum (the earlier behaviour).",
                       {"advanced"});
    defaults_.setValidStrings("fragment:deisotope_charge_cap", {"precursor", "none"});
    defaults_.setValue("fragment:deisotope_sum_intensity", "true",
                       "Give the monoisotopic peak of each isotope envelope the summed intensity of the envelope (Sage-like). "
                       "'false' keeps the monoisotopic peak's own intensity (the earlier behaviour).",
                       {"advanced"});
    defaults_.setValidStrings("fragment:deisotope_sum_intensity", {"true", "false"});


    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("fragment:min_ion_index", 2, "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.");
    defaults_.setMinInt("fragment:min_ion_index", 0);

    defaults_.setValue("peaks:keep_n", 0, "Maximum number of MS2 peaks kept per spectrum (NLargest) before scoring. 0 = auto: a resolution-aware cap (~400 for high-res, ~80 at 0.5 Da, floor 60) that avoids the spurious low-resolution fragment matching which otherwise inflates HyperScore for targets and decoys alike. Set an explicit value to override (e.g. 400 for the legacy behavior).", {"advanced"});
    defaults_.setMinInt("peaks:keep_n", 0);
    defaults_.setValue("peaks:window_top", 20, "Maximum number of MS2 peaks kept per 100 Da window (WindowMower) before scoring.", {"advanced"});
    defaults_.setMinInt("peaks:window_top", 1);
    defaults_.setValue(
      "fragment:query_spectrum", "auto",
      "Spectrum used for candidate retrieval. 'raw' queries the fragment index with all positive-intensity peaks, before "
      "deisotoping and local/top-N filtering, so that spectra whose matching fragments are not among the retained peaks still get "
      "candidates; it costs an additional peak list per spectrum and more fragment matches. 'processed' queries with the scoring peak "
      "list. 'auto' uses 'raw' for high-resolution fragments (<= 0.1 Da / <= 100 ppm) and 'processed' otherwise: at low resolution, "
      "the unfiltered peaks match many random fragments and crowd out the correct candidates. Scoring always uses the processed spectrum.",
      {"advanced"});
    defaults_.setValidStrings("fragment:query_spectrum", {"auto", "raw", "processed"});
    defaults_.setValue("peaks:window_type", "auto",
                       "Local peak filtering. 'auto' uses jump_full for high-resolution fragments (<= 0.1 Da / <= 100 ppm), "
                       "preserving the strongest peaks even in a short final 100 Da window, and legacy jump filtering otherwise. "
                       "'jump' scales the last window's peak quota by its observed width; 'jump_full' keeps the full quota in every window.",
                       {"advanced"});
    defaults_.setValidStrings("peaks:window_type", {"auto", "jump", "jump_full"});

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<std::string> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'. A terminus carries one modification: a variable modification of the whole terminus (e.g. 'Acetyl (Protein N-term)') is not searched where a fixed one sits on it (e.g. 'TMT6plex (N-term)').");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    vector<std::string> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));

    defaults_.setValue("decoys", "auto",
      "Decoy handling for target-decoy FDR. "
      "'auto': ensure decoys are available — reuse decoys already present in the database "
      "(marker auto-detected, prefix or suffix) or generate them if none are found. "
      "'generate': always (re)generate decoys from the target proteins (any pre-existing "
      "decoys are removed first to avoid decoy-of-decoy entries). "
      "'ignore': search the target proteins only — any decoys present in the database are removed.");
    defaults_.setValidStrings("decoys", {"auto","generate","ignore"} );

    defaults_.setValue("decoy_prefix", "DECOY_", "Accession prefix prepended to decoy proteins when ProSE generates them (modes 'auto' without existing decoys, and 'generate'). Pre-existing decoys in the database are detected automatically (any common marker, used as prefix or suffix) and do not require this setting.", {"advanced"});

    defaults_.setValue("annotate:PSM",  std::vector<std::string>{"ALL"}, "Annotations added to each PSM.");
    defaults_.setValidStrings("annotate:PSM",
      std::vector<std::string>{
        "ALL",
        Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM,
        Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
        Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
        Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION,
        Constants::UserParam::NUM_MATCHED_PEAKS,
        Constants::UserParam::MATCHED_PREFIX_IONS,
        Constants::UserParam::MATCHED_SUFFIX_IONS,
        Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE,
        Constants::UserParam::MATCHED_ION_CURRENT,
        Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM,
        Constants::UserParam::HYPERSCORE_ZSCORE,
        Constants::UserParam::LN_NUM_CANDIDATES,
        Constants::UserParam::MATCHED_ION_CURRENT_FRACTION,
        Constants::UserParam::COMPLEMENTARY_IONS_FRACTION}
      );

    defaults_.setValue("annotate:self_trained_ion_priors", "true",
      "Learn per-run fragment ion likelihoods from the confident PSMs of each file and add the Percolator features of "
      "annotate:ion_prior_features to every hit. The spectra are split into two halves "
      "by scan parity. Each half trains a model on its rank-one target hits at target-decoy competition q <= "
      "annotate:ion_prior_train_fdr of the native score (their reversed sequences on the same spectra are the noise "
      "model), and every hit is scored by the model of the other half, so no PSM is scored by a model that saw its "
      "spectrum or its label. Needs annotate:ion_prior_min_psms training PSMs in each half, otherwise the features are 0. "
      "Not applied to target-only searches (decoys=ignore), which get no features. Native scores and candidate selection "
      "are unchanged. Independent of annotate:PSM.", {"advanced"});
    defaults_.setValidStrings("annotate:self_trained_ion_priors", {"true", "false"});
    defaults_.setValue("annotate:ion_prior_model", "rich",
      "Contexts and outcomes of the ion priors. 'rich': ion series, precursor charge, fragment charge, relative position, "
      "the residues at the cleavage site (bond before P, bond after D or E) and whether the complementary singly charged "
      "ion matched; outcome: intensity rank x mass error (|error| / tolerance below 0.25, below 0.5, up to 1). 'basic': "
      "ion series, precursor charge, fragment charge and relative position; outcome: intensity rank.", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_model", {"rich", "basic"});
    defaults_.setValue("annotate:ion_prior_peaks", "all",
      "Peaks the ion priors rank and match. 'all': every peak after deisotoping, before the peaks:window_top and "
      "peaks:keep_n filters, so weak fragment ions count as evidence for the model without entering HyperScore. "
      "'scored': the peaks that are scored.", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_peaks", {"all", "scored"});
    defaults_.setValue("annotate:ion_prior_max_fragment_charge", 2,
      "Highest fragment charge of the ions of the ion priors; capped at the precursor charge - 1 (at least 1).", {"advanced"});
    defaults_.setMinInt("annotate:ion_prior_max_fragment_charge", 1);
    defaults_.setValue("annotate:ion_prior_train_fdr", 0.01,
      "Target-decoy competition q-value threshold selecting the training PSMs of annotate:self_trained_ion_priors.", {"advanced"});
    defaults_.setMinFloat("annotate:ion_prior_train_fdr", 0.0);
    defaults_.setMaxFloat("annotate:ion_prior_train_fdr", 1.0);
    defaults_.setValue("annotate:ion_prior_min_psms", 100,
      "Minimum number of training PSMs of annotate:self_trained_ion_priors in each half of the spectra; with fewer, the "
      "features are 0.", {"advanced"});
    defaults_.setMinInt("annotate:ion_prior_min_psms", 1);
    defaults_.setValue("annotate:ion_prior_features",
      std::vector<std::string>{Constants::UserParam::ION_PRIOR_LLR, Constants::UserParam::ION_PRIOR_EXPLAINED},
      "Features of annotate:self_trained_ion_priors written to every hit and listed for Percolator: "
      "ion_prior_llr (summed log-likelihood ratio of the ion outcomes, signal vs reversed noise), ion_prior_explained "
      "(share of the predicted ion presence that was observed), ion_prior_topk_observed (share of the 6 ions most likely "
      "present that were observed).", {"advanced"});
    defaults_.setValidStrings("annotate:ion_prior_features",
      std::vector<std::string>{Constants::UserParam::ION_PRIOR_LLR, Constants::UserParam::ION_PRIOR_EXPLAINED,
                               Constants::UserParam::ION_PRIOR_TOPK_OBSERVED});

    defaults_.setValue("annotate:local_fragment_evidence", "true",
      "Add chance_match_surprise and mass_competition_evidence annotations and Percolator features. "
      "Uses local peak density and alternative intact/neutral-loss fragment assignments. "
      "Retains an additional peak list before window/top-N filtering, increasing memory use. "
      "Does not change native scores or candidate selection. Independent of annotate:PSM.", {"advanced"});
    defaults_.setValidStrings("annotate:local_fragment_evidence", {"true", "false"});
    defaults_.setSectionDescription("annotate", "Annotation Options");

    defaults_.setValue("peptide:min_size", 7, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");
    defaults_.setValue("peptide:clip_nterm_methionine", "true",
                       "Also search protein N-terminal peptides after removal of their initial methionine. "
                       "The retained-M form is still searched; length, mass and missed-cleavage limits apply to each form. "
                       "Set false for the previous search space.");
    defaults_.setValidStrings("peptide:clip_nterm_methionine", {"true", "false"});
    defaults_.setValue("peptide:deduplicate", "true",
                       "Index each exact peptidoform once before candidate selection. Protein mappings still list every protein occurrence. "
                       "Across database chunks, count each peptide/charge/isotope hypothesis once; retaining queried keys adds memory per spectrum. "
                       "SNES mother indices retain their existing behavior. Set false for occurrence-based legacy candidates.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:deduplicate", {"true", "false"});
    defaults_.setValue("peptide:protein_mapping", "index",
                       "How the hits are mapped to their proteins (peptide evidences, target_decoy, protein_references, protein hits). "
                       "'index': from the digest of the fragment index, which holds every protein occurrence of a candidate; the result is "
                       "that of PeptideIndexing, which runs instead wherever the index cannot reproduce it exactly (e.g. SNES, a specificity "
                       "other than full, protein-terminal modifications, a chunked database, symbols other than letters or long stretches "
                       "of ambiguous residues in the database). 'PeptideIndexing': always search every hit in the whole database "
                       "(Aho-Corasick), as before.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:protein_mapping", {"index", "PeptideIndexing"});
    defaults_.setValue("peptide:enzyme_specificity", "full",
      "Enzyme cleavage specificity required for both peptide termini.\n"
      "  'full' : both termini must be enzyme-specific (canonical, e.g. tryptic).\n"
      "  'semi' : only one terminus needs to be enzyme-specific (semi-tryptic).\n"
      "  'none' : no enzyme constraint at either terminus; every substring of length\n"
      "           [min_size, max_size] is enumerated. Use this for immunopeptidomics\n"
      "           (e.g. HLA peptides 8..12mers). For very large search spaces consider\n"
      "           tightening 'peptide:min_size'/'peptide:max_size'.");
    defaults_.setValidStrings("peptide:enzyme_specificity", {"full", "semi", "none"});
    defaults_.setValue("peptide:motif", "", "If set, only peptides that contain this motif (provided as RegEx) will be considered.");
    defaults_.setSectionDescription("peptide", "Peptide Options");

    // SNES (Speedy Non-specific Enzyme Search): forwarded to FragmentIndex. Only
    // takes effect when peptide:enzyme_specificity=none. v1.1 opt-in (default false).
    defaults_.setValue("snes_enabled", "false",
      "[experimental, v1.1 opt-in] When peptide:enzyme_specificity=none, use mother-"
      "peptide indexing (Single-N + Single-C, one ion series per mother) instead of "
      "naïve O(L^2) sub-peptide enumeration. Much smaller index and faster search on "
      "non-specific workloads (immunopeptidomics). Ignored for specific/semi-"
      "specific enzymes. Variable modifications are supported in v1.1 via query-time "
      "subset enumeration on the realized sub-peptide.");
    defaults_.setValidStrings("snes_enabled", {"true", "false"});

    defaults_.setValue("report:top_hits", 1, "Maximum number of top scoring hits per spectrum that are reported.");
    defaults_.setValue("report:isotope_error_convention", "observed_minus_theoretical",
                       "Sign of the precursor isotope error reported for each PSM (meta value 'isotope_error', in 13C "
                       "spacings of 1.00336 Da). 'observed_minus_theoretical': +1 means the observed precursor lies one "
                       "13C spacing above the peptide (its first 13C isotope peak was selected), as in MS-GF+, pepXML and "
                       "Sage; PercolatorInfile and PercolatorAdapter (-out_pin, -percolator_executable) remove the offset "
                       "from the precursor mass difference with this sign. The search parameters of the result record it "
                       "as 'isotope_error_convention'. 'theoretical_minus_observed': the opposite sign, as written by "
                       "earlier ProSE versions (no record in the search parameters); with it the Percolator input "
                       "doubles the offset in dm/absdm instead of removing it.",
                       {"advanced"});
    defaults_.setValidStrings("report:isotope_error_convention", {"observed_minus_theoretical", "theoretical_minus_observed"});
    defaults_.setSectionDescription("report", "Reporting Options");

    defaults_.setValue("FDR:PSM", 0.0, "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, set to 0 to disable filtering and report all PSMs with q-values). Target and decoy PSMs are filtered alike by the q-value threshold; decoys that pass are kept (no decoy-specific stripping here — decoys are removed only at protein-FDR finalization). Requires '-decoys' to be set.");
    defaults_.setMinFloat("FDR:PSM", 0.0);
    defaults_.setMaxFloat("FDR:PSM", 1.0);
    defaults_.setValue("FDR:PSM_groups", "scored_charges",
                       "Target-decoy competitions of FDR:PSM. 'pooled' ranks all PSMs together. 'scored_charges' computes q-values "
                       "separately for PSMs scored with different numbers of fragment charges (scoring:fragment_charges): HyperScore "
                       "grows with the number of theoretical ions, so when multiply charged fragments are scored for precursors of "
                       "charge 3 and above but not for charge 2, a pooled competition is dominated by the higher precursor charges. "
                       "With singly charged fragments only (deisotoped, high-resolution spectra) there is one group, as with 'pooled'. "
                       "With fragment:max_charge above 2, precursor charges of 4 and above form further groups, which accept few PSMs "
                       "per run (tens), so their q-values (D/T, as FalseDiscoveryRate) are coarse. "
                       "A group without decoy or without target PSMs falls back to 'pooled'.",
                       {"advanced"});
    defaults_.setValidStrings("FDR:PSM_groups", {"scored_charges", "pooled"});
    defaults_.setValue("FDR:protein", 0.0, "Filter proteins based on picked-protein FDR q-value (e.g., 0.01 = 1% protein FDR, set to 0 to disable). Applied after PSM-level FDR on a complete protein set (single file, or the -out_merged aggregate). Setting this > 0 finalizes the result: identified decoys are removed. With 0, decoys are retained for downstream/merged FDR. Uses the picked-protein approach (Savitski et al. 2015) which pairs target and decoy proteins by accession. Requires '-decoys' to be set.");
    defaults_.setMinFloat("FDR:protein", 0.0);
    defaults_.setMaxFloat("FDR:protein", 1.0);
    defaults_.setSectionDescription("FDR", "False Discovery Rate control (requires decoys)");

    // Add parameters which are only used by FragmentIndex
    defaults_.setValue("peptide:min_mass", 100, "Minimal peptide mass for database");
    defaults_.setValue("peptide:max_mass", 9000, "Maximal peptide mass for database");

    // Fragment-level filtering
    defaults_.setValue("fragment:min_matched_ions", 5, "Minimal number of matched ions to report a PSM");

    // Precursor isotope error handling. The search convention (theoretical minus observed) is kept for these two
    // parameters; the reported isotope_error follows report:isotope_error_convention.
    defaults_.setValue("precursor:isotope_error_min", -1,
                       "Minimum precursor isotope error searched, in 13C spacings (1.00336 Da) added to the observed "
                       "precursor mass: -1 finds a peptide whose first 13C isotope peak was selected as the precursor. "
                       "The PSM meta value 'isotope_error' reports the opposite sign by default (see "
                       "report:isotope_error_convention).");
    defaults_.setValue("precursor:isotope_error_max", 1,
                       "Maximum precursor isotope error searched, with the sign of precursor:isotope_error_min: +1 finds "
                       "a peptide whose precursor was selected one 13C spacing below its monoisotopic peak.");

    // Fragment and scoring limits
    defaults_.setValue("fragment:max_charge", 2, "max fragment charge");
    defaults_.setValue("scoring:max_candidates_per_spectrum", 50, "The number of initial hits for which we calculate a score");
    defaults_.setValue(
      "scoring:fragment_charges", "auto",
      "Fragment charges of the theoretical spectra that score each candidate. 'single' scores singly charged fragments only. "
      "'multiple' adds charges up to min(precursor charge - 1, fragment:max_charge). 'auto' uses 'multiple' when the MS2 spectra "
      "are not deisotoped (fragment:deisotope=false, or a fragment tolerance above the deisotoper's limit of 0.1 Da / 100 ppm, "
      "as for ion-trap CID) and 'single' otherwise, because deisotoping converts fragments to charge 1.",
      {"advanced"});
    defaults_.setValidStrings("scoring:fragment_charges", {"auto", "single", "multiple"});
    defaults_.setValue("scoring:method", "auto",
                       "Native score. 'hyperscore' is the X!Tandem HyperScore. 'mass_accuracy' weights each matched fragment's "
                       "contribution to HyperScore's ion counts and intensity product by exp(-0.5 * (error / scoring:mass_error_sd)^2), "
                       "with the error in ppm, so that matches near the theoretical m/z count fully and matches at the edge of the "
                       "tolerance count little; exact matches score as HyperScore. The kernel is centred at 0 and suits high-resolution "
                       "fragments only. 'auto' (default) uses 'mass_accuracy' for fragment tolerances in ppm within the high-resolution "
                       "range (<= 100 ppm) and 'hyperscore' otherwise (Da tolerances, e.g. ion-trap CID). Matching tolerances, candidates "
                       "and the unweighted match annotations are unchanged. The score type stays 'ln(hyperscore)'; the search parameters "
                       "record scoring:method_resolved when the weighted score is used.",
                       {"advanced"});
    defaults_.setValidStrings("scoring:method", {"hyperscore", "mass_accuracy", "auto"});
    defaults_.setValue("scoring:mass_error_sd", 7.0,
                       "Standard deviation in ppm of the Gaussian fragment mass-error kernel of scoring:method 'mass_accuracy'. "
                       "Fixed (not fitted) and independent of the matching tolerance and of calibration.",
                       {"advanced"});
    defaults_.setMinFloat("scoring:mass_error_sd", 1e-6);
    defaults_.setSectionDescription("scoring", "Search/Scoring Limits");

    // Ion series toggles
    defaults_.setValue("ions:add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("ions:add_y_ions", {"true","false"});
    defaults_.setValue("ions:add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("ions:add_b_ions", {"true","false"});
    defaults_.setValue("ions:add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("ions:add_a_ions", {"true","false"});
    defaults_.setValue("ions:add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("ions:add_c_ions", {"true","false"});
    defaults_.setValue("ions:add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("ions:add_x_ions", {"true","false"});
    defaults_.setValue("ions:add_z_ions", "false", "Add peaks of z-ions (y - NH3) to the spectrum. These are not the z+1 ions of ETD-type spectra; see ions:by_activation and ions:add_zp1_ions.");
    defaults_.setValidStrings("ions:add_z_ions", {"true","false"});
    defaults_.setValue("ions:add_zp1_ions", "false", "Add peaks of z+1 ions (z-dot, y - NH2) to the spectrum, the main C-terminal fragments of ETD, EThcD and ETciD spectra, for all spectra. ions:by_activation adds them (with c ions) for spectra recorded as electron-activated; set this, typically with ions:add_c_ions, for ETD-type data without activation information.");
    defaults_.setValidStrings("ions:add_zp1_ions", {"true","false"});
    defaults_.setValue("ions:by_activation", "true", "Score spectra whose precursor was activated by electrons (ETD, ECD, EThcD or ETciD, as recorded in the input file) with c and z+1 ions in addition to the ion series above. Other spectra, and spectra without activation information, use the ion series above.");
    defaults_.setValidStrings("ions:by_activation", {"true","false"});
    defaults_.setSectionDescription("ions", "Theoretical ion series toggles");

    defaults_.setValue("calibration:enabled", "false",
      "If enabled, run a fast calibration pass on a subset of spectra before the main search. "
      "Estimates tighter precursor and fragment tolerances from confident PSMs. "
      "The fragment index is NOT rebuilt — only query-time tolerances are tightened. "
      "Inspired by MSFragger's calibrate_mass and OpenNuXL's autotune.");
    defaults_.setValidStrings("calibration:enabled", {"true", "false"});
    defaults_.setValue("calibration:subset_ratio", 0.1,
      "Fraction of spectra (by TIC, highest first) used for the calibration pass (0.0-1.0).");
    defaults_.setMinFloat("calibration:subset_ratio", 0.01);
    defaults_.setMaxFloat("calibration:subset_ratio", 1.0);
    defaults_.setValue("calibration:min_psms", 50,
      "Minimum number of confident PSMs required for calibration. If fewer are found, "
      "calibration is skipped and the user-configured tolerances are used.");
    defaults_.setMinInt("calibration:min_psms", 1);

    defaults_.setValue("database:chunk_size", 0,
      "Split the protein database into chunks of at most this many proteins for fragment index "
      "building. 0 = disabled (load entire database at once). Enable for very large databases "
      "(e.g. MHC-II immunopeptidomics with variants) that exceed available memory. Each chunk "
      "builds its own fragment index; results are merged across chunks before post-processing. "
      "Calibration (when enabled) runs on a strided, size-bounded sample of the full "
      "decoy-augmented database — sample size is tied to chunk_size so calibration memory "
      "respects the same budget as main-search chunks. In multi-file mode (-in a.mzML b.mzML), "
      "the chunk-major path builds each chunk's fragment index once and scores all files "
      "against it before moving to the next chunk. Note that chunk_size bounds only the "
      "fragment-index memory: the chunk-major schedule holds ALL input files' preprocessed "
      "MS2 spectra in memory for the whole search, so budget for their sum as well, or split "
      "very large cohorts across separate ProSE invocations.");
    defaults_.setMinInt("database:chunk_size", 0);

    defaults_.setSectionDescription("calibration",
      "Automatic mass accuracy calibration (two-pass search). A fast first pass on a subset of "
      "spectra estimates instrument-specific mass accuracy, then the main search uses the "
      "calibrated (typically tighter) tolerances for better discrimination.");

    defaultsToParam_();
  }

  void ProSEAlgorithm::updateMembers_()
  {
    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();

    precursor_min_charge_ = param_.getValue("precursor:min_charge");
    precursor_max_charge_ = param_.getValue("precursor:max_charge");

    peaks_keep_n_ = (Size)(int)param_.getValue("peaks:keep_n");
    peaks_window_top_ = (Int)param_.getValue("peaks:window_top");
    peaks_window_type_ = param_.getValue("peaks:window_type").toString();

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");
    // Every scoring method and the fragment index match within this tolerance (0, negative and NaN tolerances match
    // nothing, infinity everything), and HyperScore::computeMassAccuracy() and local fragment evidence require it.
    // Checked here, before the parallel scoring loops.
    if (! std::isfinite(fragment_mass_tolerance_) || fragment_mass_tolerance_ <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "fragment:mass_tolerance must be finite and positive (got " + std::to_string(fragment_mass_tolerance_) + ").");
    }

    fragment_mass_tolerance_unit_ = param_.getValue("fragment:mass_tolerance_unit").toString();

    // Resolve the MS2 deisotoping mode (fragment:deisotope = auto|true|false) ONCE.
    // The Deisotoper supports a fragment tolerance <= 100 ppm / <= 0.1 Da and throws
    // otherwise; consult its non-throwing predicate so the limit lives in one place
    // (OpenMS#9619). 'false' never deisotopes; 'true'/'auto' deisotope only when the
    // tolerance is supported -- 'true' fails fast here with a clear parameter error
    // rather than aborting later inside preprocessSpectra_'s OpenMP region, while
    // 'auto' silently skips deisotoping for low-resolution data.
    const std::string deisotope_mode = param_.getValue("fragment:deisotope").toString();
    const bool deisotope_supported =
      Deisotoper::isToleranceSupported(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm");
    // Resolved from the configured tolerance, with the deisotoper's resolution boundary; calibration does not switch it.
    const std::string query_spectrum = param_.getValue("fragment:query_spectrum").toString();
    query_raw_spectrum_ = query_spectrum == "raw" || (query_spectrum == "auto" && deisotope_supported);
    deisotope_requested_ = (deisotope_mode != "false");
    if (deisotope_mode == "true" && !deisotope_supported)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "fragment:deisotope=true requires a high-resolution fragment tolerance (<= 0.1 Da or <= 100 ppm). "
        "Use 'auto' (deisotopes only for high-resolution data) or 'false' for low-resolution data.");
    }
    if (deisotope_mode == "auto" && !deisotope_supported)
    {
      OPENMS_LOG_WARN << "[ProSE] Fragment tolerance " << fragment_mass_tolerance_
                      << (fragment_mass_tolerance_unit_ == "ppm" ? " ppm" : " Da")
                      << " exceeds the deisotoping limit (100 ppm / 0.1 Da); skipping MS2 "
                      << "deisotoping (expected for low-resolution data)." << endl;
    }
    deisotoping_.min_peaks = static_cast<unsigned int>(static_cast<int>(param_.getValue("fragment:deisotope_min_peaks")));
    deisotoping_.charge_cap_precursor = param_.getValue("fragment:deisotope_charge_cap").toString() == "precursor";
    deisotoping_.sum_intensity = param_.getValue("fragment:deisotope_sum_intensity").toBool();

    // Spectra that are not deisotoped keep multiply charged fragments at their own m/z;
    // deisotoped spectra hold charge-1 fragments only.
    const std::string fragment_charges = param_.getValue("scoring:fragment_charges").toString();
    const bool deisotoped = deisotope_requested_ && deisotope_supported;
    scoring_multiple_charges_ = fragment_charges == "multiple" || (fragment_charges == "auto" && ! deisotoped);
    scoring_max_charge_ = static_cast<int>(param_.getValue("fragment:max_charge"));

    // Mass-accuracy weighting needs fragment errors of a few ppm; 'auto' uses it for high-resolution ppm tolerances only.
    const std::string scoring_method = param_.getValue("scoring:method").toString();
    const bool high_resolution_ppm = fragment_mass_tolerance_unit_ == "ppm" && deisotope_supported;
    mass_accuracy_score_ = scoring_method == "mass_accuracy" || (scoring_method == "auto" && high_resolution_ppm);
    mass_error_sd_ppm_ = param_.getValue("scoring:mass_error_sd");
    if (mass_accuracy_score_ && (! std::isfinite(mass_error_sd_ppm_) || mass_error_sd_ppm_ <= 0.0))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "scoring:mass_error_sd must be finite and positive.");
    }
    if (scoring_method == "mass_accuracy" && ! high_resolution_ppm)
    {
      OPENMS_LOG_WARN << "[ProSE] scoring:method=mass_accuracy weights fragment matches by a " << mass_error_sd_ppm_
                      << " ppm kernel, but the fragment tolerance is " << fragment_mass_tolerance_
                      << (fragment_mass_tolerance_unit_ == "ppm" ? " ppm" : " Da")
                      << "; matches far from the theoretical m/z score little. 'auto' uses HyperScore for such tolerances." << endl;
    }

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    set<std::string> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());
    if (fixed_unique.size() != modifications_fixed_.size())
    {
      OPENMS_LOG_WARN << "Duplicate fixed modification provided. Making them unique." << endl;
      modifications_fixed_.assign(fixed_unique.begin(), fixed_unique.end());
    }
    // fixed terminal modifications the fragment index cannot restrict: fail here, not after reading the database
    FragmentIndex::checkFixedModifications(modifications_fixed_);

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<std::string> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      OPENMS_LOG_WARN << "Duplicate variable modification provided. Making them unique." << endl;
      modifications_variable_.assign(var_unique.begin(), var_unique.end());
    }
    for (const std::string& mod : FragmentIndex::shadowedVariableTerminalModifications(modifications_fixed_, modifications_variable_))
    {
      OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: a fixed modification already sits on "
                      << "that terminus, which carries one modification. To search it, specify the fixed terminal "
                      << "modification as a variable one as well." << endl;
    }

    modifications_max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    enzyme_ = param_.getValue("enzyme").toString();

    peptide_min_size_ = param_.getValue("peptide:min_size");
    peptide_max_size_ = param_.getValue("peptide:max_size");
    peptide_missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_enzyme_specificity_ = EnzymaticDigestion::getSpecificityByName(
      param_.getValue("peptide:enzyme_specificity").toString());
    peptide_motif_ = param_.getValue("peptide:motif").toString(); // TODO: remove unused parameters

    report_top_hits_ = param_.getValue("report:top_hits");
    isotope_error_observed_minus_theoretical_ =
      param_.getValue("report:isotope_error_convention").toString() == "observed_minus_theoretical";

    const std::string decoy_mode_str = param_.getValue("decoys").toString();
    if (decoy_mode_str == "generate")   { decoy_mode_ = DecoyMode_::GENERATE; }
    else if (decoy_mode_str == "ignore") { decoy_mode_ = DecoyMode_::IGNORE; }
    else                                 { decoy_mode_ = DecoyMode_::AUTO; }
    decoy_prefix_ = param_.getValue("decoy_prefix").toString();
    annotate_psm_ = ListUtils::toStringList<std::string>(param_.getValue("annotate:PSM"));
    self_trained_ion_priors_ = param_.getValue("annotate:self_trained_ion_priors").toBool();
    ion_prior_model_ = param_.getValue("annotate:ion_prior_model").toString() == "basic"
                         ? FragmentIonLikelihoodModel::ContextSet::BASIC : FragmentIonLikelihoodModel::ContextSet::RICH;
    ion_prior_scored_peaks_ = param_.getValue("annotate:ion_prior_peaks").toString() == "scored";
    ion_prior_max_fragment_charge_ = static_cast<int>(param_.getValue("annotate:ion_prior_max_fragment_charge"));
    ion_prior_train_fdr_ = param_.getValue("annotate:ion_prior_train_fdr");
    ion_prior_min_psms_ = static_cast<Size>(static_cast<int>(param_.getValue("annotate:ion_prior_min_psms")));
    {
      const StringList features = ListUtils::toStringList<std::string>(param_.getValue("annotate:ion_prior_features"));
      auto selected = [&features](const std::string& name) { return std::find(features.begin(), features.end(), name) != features.end(); };
      ion_prior_feature_llr_ = selected(Constants::UserParam::ION_PRIOR_LLR);
      ion_prior_feature_explained_ = selected(Constants::UserParam::ION_PRIOR_EXPLAINED);
      ion_prior_feature_topk_ = selected(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED);
      if (self_trained_ion_priors_ && !(ion_prior_feature_llr_ || ion_prior_feature_explained_ || ion_prior_feature_topk_))
      {
        OPENMS_LOG_WARN << "[ProSE] annotate:self_trained_ion_priors is set, but annotate:ion_prior_features is empty: "
                        << "no ion priors are learned." << endl;
        self_trained_ion_priors_ = false;
      }
    }
    if (self_trained_ion_priors_ && decoy_mode_ == DecoyMode_::IGNORE)
    {
      // A target-only search has no target-decoy competition to select training PSMs, and Percolator cannot use its
      // output: nothing is learned or written (the output equals the one with the ion priors switched off).
      OPENMS_LOG_INFO << "[ProSE] annotate:self_trained_ion_priors: not applied to a target-only search (decoys=ignore)." << endl;
      self_trained_ion_priors_ = false;
    }
    fdr_psm_ = param_.getValue("FDR:PSM");
    fdr_psm_by_scored_charges_ = param_.getValue("FDR:PSM_groups").toString() == "scored_charges";
    fdr_protein_ = param_.getValue("FDR:protein");

    // Open search mode is automatically determined based on precursor tolerance in isOpenSearchMode_()

    add_a_ions_ = param_.getValue("ions:add_a_ions").toBool();
    add_b_ions_ = param_.getValue("ions:add_b_ions").toBool();
    add_c_ions_ = param_.getValue("ions:add_c_ions").toBool();
    add_x_ions_ = param_.getValue("ions:add_x_ions").toBool();
    add_y_ions_ = param_.getValue("ions:add_y_ions").toBool();
    add_z_ions_ = param_.getValue("ions:add_z_ions").toBool();
    add_zp1_ions_ = param_.getValue("ions:add_zp1_ions").toBool();
    ions_by_activation_ = param_.getValue("ions:by_activation").toBool();
    protein_mapping_from_index_ = param_.getValue("peptide:protein_mapping").toString() == "index";

    database_chunk_size_ = param_.getValue("database:chunk_size");

    calibration_enabled_ = param_.getValue("calibration:enabled") == "true";
    calibration_subset_ratio_ = param_.getValue("calibration:subset_ratio");
    calibration_min_psms_ = param_.getValue("calibration:min_psms");
  }

  // static
  struct ProSEAlgorithm::SpectrumGenerators_
  {
    TheoreticalSpectrumGenerator standard; ///< configured ion series
    TheoreticalSpectrumGenerator electron; ///< configured ion series plus c and z+1 ions
    bool by_activation = false;

    /// True if @p spectrum is matched against, and scored with, c and z+1 ions as well: its
    /// precursor was activated by electrons
    bool electronIons(const MSSpectrum& spectrum) const
    {
      return by_activation && isElectronActivated_(spectrum);
    }

    /// The generator for @p spectrum, chosen by the activation of its precursor
    const TheoreticalSpectrumGenerator& forSpectrum(const MSSpectrum& spectrum) const
    {
      return electronIons(spectrum) ? electron : standard;
    }
  };

  ProSEAlgorithm::SpectrumGenerators_ ProSEAlgorithm::spectrumGenerators_() const
  {
    // Scoring, annotation and calibration use the ion series of the fragment index, so that
    // a candidate is scored against the ions it was retrieved by.
    auto configure = [this](TheoreticalSpectrumGenerator& tsg, bool electron)
    {
      Param p(tsg.getParameters());
      p.setValue("add_first_prefix_ion", "true");
      p.setValue("add_metainfo", "true");
      p.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
      p.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
      p.setValue("add_c_ions", (add_c_ions_ || electron) ? "true" : "false");
      p.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
      p.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
      p.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
      p.setValue("add_zp1_ions", (add_zp1_ions_ || electron) ? "true" : "false");
      tsg.setParameters(p);
    };
    SpectrumGenerators_ generators;
    configure(generators.standard, false);
    configure(generators.electron, true);
    generators.by_activation = ions_by_activation_;
    return generators;
  }

  std::vector<double> ProSEAlgorithm::localPeakDensities_(const MSSpectrum& spectrum)
  {
    // Both bounds only move forward: O(number of peaks), once per spectrum.
    constexpr double half_window = 50.0;
    std::vector<double> densities(spectrum.size());
    Size left = 0, right = 0;
    for (Size i = 0; i < spectrum.size(); ++i)
    {
      const double mz = spectrum[i].getMZ();
      while (left < spectrum.size() && spectrum[left].getMZ() < mz - half_window) ++left;
      while (right < spectrum.size() && spectrum[right].getMZ() <= mz + half_window) ++right;
      densities[i] = static_cast<double>(right - left) / (2.0 * half_window);
    }
    return densities;
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical,
      const std::vector<double>& densities, double tolerance, bool ppm)
  {
    std::vector<double> alternatives;
    return localFragmentEvidence_(spectrum, theoretical, std::numeric_limits<int>::max(), densities, tolerance, ppm, alternatives);
  }

  ProSEAlgorithm::LocalFragmentEvidence_ ProSEAlgorithm::localFragmentEvidence_(
      const MSSpectrum& spectrum, const MSSpectrum& theoretical, int max_charge,
      const std::vector<double>& densities, double tolerance, bool ppm, std::vector<double>& alternatives)
  {
    // The ions of charge <= max_charge of a TSG spectrum of charges 1..Z are the TSG spectrum of charges
    // 1..max_charge, in the same order: the generator adds each charge independently and sorts stably.
    // Annotation therefore passes its own theoretical spectrum instead of generating a second one.
    LocalFragmentEvidence_ evidence;
    if (spectrum.empty() || theoretical.empty()) return evidence;
    OPENMS_PRECONDITION(densities.size() == spectrum.size(), "One local density per experimental peak required")
    OPENMS_PRECONDITION(!theoretical.getIntegerDataArrays().empty()
      && theoretical.getIntegerDataArrays()[0].size() == theoretical.size(), "Aligned fragment charges required")

    static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
    static const double ammonia_mass = EmpiricalFormula("NH3").getMonoWeight();
    const auto& charges = theoretical.getIntegerDataArrays()[0];
    alternatives.clear();
    alternatives.reserve(3 * theoretical.size());
    for (Size i = 0; i < theoretical.size(); ++i)
    {
      OPENMS_PRECONDITION(charges[i] > 0, "Positive fragment charge required")
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      alternatives.push_back(mz);
      for (const double loss : {water_mass, ammonia_mass})
      {
        const double loss_mz = mz - loss / charges[i];
        if (loss_mz > 0.0) alternatives.push_back(loss_mz);
      }
    }
    std::sort(alternatives.begin(), alternatives.end());

    for (Size i = 0; i < theoretical.size(); ++i)
    {
      if (charges[i] > max_charge) continue;
      const double mz = theoretical[i].getMZ();
      const double width = ppm ? mz * tolerance * 1e-6 : tolerance;
      const Size nearest = spectrum.findNearest(mz);
      const double observed = spectrum[nearest].getMZ();
      if (std::abs(mz - observed) > width) continue;

      // The matched peak itself is in its density window, hence density > 0.
      const double density = densities[nearest];
      evidence.chance_match_surprise += std::max(0.0, -std::log(2.0 * width * density));
      const auto first = std::lower_bound(alternatives.begin(), alternatives.end(), observed - width);
      const auto last = std::upper_bound(first, alternatives.end(), observed + width);
      // One intact hypothesis is the ion currently being scored. All other
      // hypotheses, including coincident ions, represent alternative assignments.
      const Size count = static_cast<Size>(last - first);
      const Size competing = count > 0 ? count - 1 : 0;
      evidence.mass_competition_evidence += 1.0 / (1.0 + competing + density);
    }
    return evidence;
  }

  bool ProSEAlgorithm::isElectronActivated_(const MSSpectrum& spectrum)
  {
    if (spectrum.getPrecursors().empty()) return false;
    for (const Precursor::ActivationMethod method : spectrum.getPrecursors()[0].getActivationMethods())
    {
      if (method == Precursor::ActivationMethod::ETD || method == Precursor::ActivationMethod::ECD
          || method == Precursor::ActivationMethod::EThcD || method == Precursor::ActivationMethod::ETciD)
      {
        return true;
      }
    }
    return false;
  }

  Size ProSEAlgorithm::countElectronActivated_(const PeakMap& spectra) const
  {
    if (!ions_by_activation_) return 0;
    return static_cast<Size>(std::count_if(spectra.begin(), spectra.end(),
                                           [](const MSSpectrum& spectrum) { return isElectronActivated_(spectrum); }));
  }

  namespace
  {
    // Keeps the peaks WindowMower::filterPeakSpectrumForTopNInJumpingWindow() keeps: the peak_count most intense peaks of
    // each m/z window of window_size that starts at a peak (of the last window a share of peak_count by its width), and
    // every peak equal to a kept one. It selects them through indices instead of copies of the peaks, and finds the
    // equal peaks among their neighbours instead of searching all kept peaks for each peak. The same std::partial_sort
    // calls on the same intensities select the same peaks.
    void filterTopNInJumpingWindows(MSSpectrum& spectrum, double window_size, UInt peak_count)
    {
      if (spectrum.empty()) { return; }
      spectrum.sortByPosition();
      const Size n = spectrum.size();
      std::vector<char> kept(n, 0);
      std::vector<Size> window;
      const auto more_intense = [&spectrum](Size a, Size b) { return spectrum[b].getIntensity() < spectrum[a].getIntensity(); };
      const auto keep_most_intense = [&](Size begin, Size end, Size count)
      {
        if (end - begin > count)
        {
          window.resize(end - begin);
          std::iota(window.begin(), window.end(), begin);
          std::partial_sort(window.begin(), window.begin() + count, window.end(), more_intense);
          for (Size k = 0; k < count; ++k) { kept[window[k]] = 1; }
        }
        else
        {
          std::fill(kept.begin() + begin, kept.begin() + end, 1);
        }
      };
      double window_start = spectrum[0].getMZ();
      Size begin = 0;
      for (Size i = 0; i != n; ++i)
      {
        if (spectrum[i].getMZ() - window_start < window_size) { continue; }
        // a gap may leave windows empty: the next window starts at the next peak
        window_start = spectrum[i].getMZ();
        keep_most_intense(begin, i, peak_count);
        begin = i;
      }
      const double last_window_fraction = (spectrum[n - 1].getMZ() - window_start) / window_size;
      keep_most_intense(begin, n, static_cast<Size>(std::round(last_window_fraction * peak_count)));

      // peaks of equal m/z are neighbours (the spectrum is sorted by m/z)
      std::vector<Size> selected;
      selected.reserve(n);
      for (Size i = 0; i < n;)
      {
        Size j = i + 1;
        while (j < n && spectrum[j].getMZ() == spectrum[i].getMZ()) { ++j; }
        for (Size a = i; a < j; ++a)
        {
          for (Size b = i; b < j; ++b)
          {
            if (kept[b] && spectrum[b] == spectrum[a])
            {
              selected.push_back(a);
              break;
            }
          }
        }
        i = j;
      }
      spectrum.select(selected);
    }
  }

  void ProSEAlgorithm::filterLocalPeaks_(MSSpectrum& spectrum, Size peaks_per_window)
  {
    // Work on indices so all peak-associated data arrays survive the selection.
    // An isolated high-m/z ion is still evidence; its window's observed width
    // must not round the retention quota down to zero.
    std::vector<Size> indices(spectrum.size()), selected;
    std::iota(indices.begin(), indices.end(), 0);
    selected.reserve(spectrum.size());
    for (Size begin = 0; begin < spectrum.size();)
    {
      Size end = begin + 1;
      while (end < spectrum.size() && spectrum[end].getMZ() - spectrum[begin].getMZ() < 100.0)
      {
        ++end;
      }
      const Size keep = std::min(peaks_per_window, end - begin);
      std::partial_sort(indices.begin() + begin, indices.begin() + begin + keep, indices.begin() + end, [&spectrum](Size a, Size b) {
        if (spectrum[a].getIntensity() != spectrum[b].getIntensity()) { return spectrum[a].getIntensity() > spectrum[b].getIntensity(); }
        return a < b;
      });
      selected.insert(selected.end(), indices.begin() + begin, indices.begin() + begin + keep);
      begin = end;
    }
    std::sort(selected.begin(), selected.end());
    spectrum.select(selected);
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type)
  {
    preprocessSpectra_(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, deisotope_requested, peaks_keep_n,
                       peaks_window_top, window_type, DeisotopingSettings_{});
  }

  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp,
                                          double fragment_mass_tolerance,
                                          bool fragment_mass_tolerance_unit_ppm,
                                          bool deisotope_requested,
                                          Size peaks_keep_n,
                                          Int peaks_window_top,
                                          const std::string& window_type,
                                          const DeisotopingSettings_& deisotoping,
                                          FragmentIonLikelihoodModel::PeakLists* ion_evidence,
                                          bool ion_evidence_scored_peaks,
                                          PeakMap* evidence_spectra,
                                          PeakMap* query_spectra)
  {
    // Intensity threshold + normalization used to run here as two extra SERIAL full-map
    // passes. Both are strictly per-spectrum: ThresholdMower::filterPeakMap and
    // Normalizer::filterPeakMap are literally "for (auto& s : exp) filterSpectrum(s);" and
    // neither iterates chromatograms. They are therefore applied at the top of the parallel
    // loop below instead, which is per-spectrum equivalent and removes two full sweeps over
    // the peak data. Both objects are configured once here; like nlargest_filter below, each
    // OpenMP thread works on its own copy (firstprivate): ThresholdMower stores its 'threshold'
    // Param in a member on every call, and concurrent writes are a data race even when every
    // thread writes the same value. One copy per thread costs a few Param copies per search.
    // Peaks without intensity (zero or negative, e.g. empty centroids) would still count as
    // matched ions, so they are removed. Nothing else is (every positive float intensity,
    // subnormal ones included, passes): the ThresholdMower default of 0.05 is an absolute
    // cutoff on the raw intensities before normalization and deleted real peaks from
    // intensity-scaled input (e.g. spectra normalized to a base peak of 1).
    ThresholdMower threshold_mower_filter;
    Param threshold_param = threshold_mower_filter.getParameters();
    threshold_param.setValue("threshold", static_cast<double>(std::numeric_limits<float>::denorm_min()));
    threshold_mower_filter.setParameters(threshold_param);
    Normalizer normalizer;

    // sort by rt; done before the loop because sortSpectra(false) only permutes whole
    // spectra (it does not touch peak order, see MSExperiment::sortSpectra) and therefore
    // commutes with the per-spectrum filters that now run inside the loop.
    exp.sortSpectra(false);
    if (ion_evidence != nullptr) { ion_evidence->reset(exp.size()); }
    if (evidence_spectra != nullptr)
    {
      evidence_spectra->clear(true);
      evidence_spectra->resize(exp.size());
    }
    if (query_spectra != nullptr)
    {
      query_spectra->clear(true);
      query_spectra->resize(exp.size());
    }

    // filter settings: the most intense peaks_window_top peaks per 100 Th window (jumping windows as WindowMower's,
    // unless full_window_quota)
    const bool full_window_quota
      = window_type == "jump_full"
        || (window_type == "auto" && Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm));

    // Resolution-aware peak retention. peaks_keep_n == 0 => auto. For HIGH-resolution fragments
    // (within the deisotoper range, <= 0.1 Da / <= 100 ppm) keep the legacy cap of 400. For
    // LOW-resolution data (e.g. 0.5 Da ion-trap CID — the same regime where deisotoping is
    // skipped) a wide match window admits many spurious low-intensity matches that inflate the
    // count-sensitive HyperScore for targets AND decoys alike, collapsing the target/decoy
    // margin; retain far fewer peaks. Random matches scale with peak density x tolerance width,
    // so sqrt-dampen around the 0.02 Da reference and clamp to [60, 400]:
    //   0.2 Da -> 126, 0.5 Da -> 80, 1 Da -> 60.
    Size effective_keep_n = peaks_keep_n;
    if (effective_keep_n == 0)
    {
      if (Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm))
      {
        effective_keep_n = 400; // high-resolution: unchanged legacy behavior
      }
      else
      {
        const double eff_tol_da = fragment_mass_tolerance_unit_ppm
          ? fragment_mass_tolerance * 500.0 * 1e-6   // ppm -> Da at a 500 m/z reference
          : fragment_mass_tolerance;
        const double n = std::round(400.0 * std::sqrt(0.02 / std::max(eff_tol_da, 1e-6)));
        effective_keep_n = static_cast<Size>(std::min(400.0, std::max(60.0, n)));
      }
    }
    NLargest nlargest_filter = NLargest(effective_keep_n);

    // Deisotope only when requested (param fragment:deisotope, resolved in
    // updateMembers_) AND the tolerance is in the Deisotoper's supported range.
    // The range guard keeps the Deisotoper call below from ever throwing, so no
    // exception can escape the OpenMP region (OpenMS#9619). Mode resolution and the
    // low-resolution warning are handled once, in updateMembers_.
    const bool do_deisotope = deisotope_requested &&
      Deisotoper::isToleranceSupported(fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);

#pragma omp parallel for default(none) shared(exp, evidence_spectra, query_spectra, do_deisotope, fragment_mass_tolerance, \
                                                fragment_mass_tolerance_unit_ppm, full_window_quota, peaks_window_top, \
                                                deisotoping, ion_evidence, ion_evidence_scored_peaks) \
                                         firstprivate(threshold_mower_filter, normalizer, nlargest_filter)
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // remove 0 intensities, then normalize (formerly two serial full-map passes)
      threshold_mower_filter.filterPeakSpectrum(exp[exp_index]);
      normalizer.filterPeakSpectrum(exp[exp_index]);

      // sort by mz
      exp[exp_index].sortByPosition();

      // fragment:query_spectrum=raw retrieves candidates with this peak list (before deisotoping and local filtering).
      if (query_spectra != nullptr) { (*query_spectra)[exp_index] = exp[exp_index]; }

      // deisotope (skipped for low-resolution data; see do_deisotope above)
      // Isotope intensities must fall from the monoisotopic peak on (start_intensity_check = 1).
      // With the library default of 2, a small peak one isotope spacing below a fragment ion
      // became the envelope's monoisotopic peak and the ion itself was removed as its isotope.
      // TMT/TMTpro-labelled fragments carry such a peak (reagent isotope impurity), and dense
      // Orbitrap Astral and timsTOF spectra often hold one by chance.
      // The envelope rule is set by fragment:deisotope_min_peaks, _charge_cap and _sum_intensity; their defaults
      // are Sage-like (two peaks, charges up to the precursor charge but at most 3, summed intensity; the earlier
      // rule was three peaks, charges 1-3, own intensity). The Sage-like rule differs from Sage itself: Sage links isotope
      // pairs with a fixed 10 ppm tolerance, requires each isotope peak to be weaker than its parent and does not
      // cap the charge at 3, while the OpenMS deisotoper extends an envelope from its monoisotopic peak with the
      // search tolerance and stops at the first isotope peak that is more intense than its predecessor.
      // The deisotoper does not exclude peaks that already belong to an earlier envelope, so an isotope peak can be
      // claimed by two monoisotopic peaks and, with summing, counted in both (two-peak envelopes make this more
      // frequent; Sage likewise adds a peak to several parents of the same charge). It does not depend on labels.
      if (do_deisotope)
      {
        int max_charge = 3;
        if (deisotoping.charge_cap_precursor && ! exp[exp_index].getPrecursors().empty())
        {
          const int precursor_charge = exp[exp_index].getPrecursors()[0].getCharge();
          if (precursor_charge > 0) { max_charge = std::min(max_charge, precursor_charge); }
        }
        Deisotoper::deisotopeAndSingleCharge(exp[exp_index],
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
          1, max_charge,  // min / max charge
          false,  // keep only deisotoped
          deisotoping.min_peaks, 10,  // min / max isopeaks
          true,   // convert fragment m/z to mono-charge
          false,  // annotate charge
          false,  // annotate isotopic peak counts
          true,   // decreasing isotope intensities
          1,      // start the intensity check at the monoisotopic peak
          deisotoping.sum_intensity);  // the monoisotopic peak carries the envelope's intensity
      }

      // ion priors: every peak after deisotoping, before the window and top-N filters
      if (ion_evidence != nullptr && !ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
      // Local fragment evidence uses this peak list: after deisotoping, so that densities and
      // matching share the m/z space of the scoring spectrum, but before local and top-N filtering.
      if (evidence_spectra != nullptr) { (*evidence_spectra)[exp_index] = exp[exp_index]; }

      // remove noise
      if (full_window_quota) { filterLocalPeaks_(exp[exp_index], static_cast<Size>(peaks_window_top)); }
      else { filterTopNInJumpingWindows(exp[exp_index], 100.0, static_cast<UInt>(peaks_window_top)); }
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();

      // ion priors on the scored peaks
      if (ion_evidence != nullptr && ion_evidence_scored_peaks) { ion_evidence->assign(exp_index, exp[exp_index]); }
    }
  }

  double ProSEAlgorithm::CandidatePoolStats_::zScore() const
  {
    if (count < 2) return 0.0;
    const double n = static_cast<double>(count);
    const double mean = sum / n;
    const double var = std::max(0.0, sumsq / n - mean * mean);
    if (var <= 0.0) return 0.0; // every candidate scored the same: the best one is not an outlier
    return (best - mean) / std::sqrt(var);
  }

  namespace
  {
    /// Scratch buffers of alignAbsoluteTolerance_(), reused across the hits of one spectrum.
    struct AlignmentScratch_
    {
      std::vector<double> prev;         ///< scores of the computed cells of row i - 1
      std::vector<double> cur;          ///< scores of the computed cells of row i
      std::vector<unsigned char> dir;   ///< traceback direction of every computed cell, rows concatenated
      std::vector<Size> row_begin;      ///< first computed column of row i
      std::vector<Size> row_offset;     ///< position of row i in dir (entry i + 1 closes row i)
    };

    /**
      @brief Absolute-tolerance (Da) branch of SpectrumAlignment::getSpectrumAlignment on flat storage.

      Returns exactly the alignment of the banded dynamic programme in SpectrumAlignment.h, which keeps
      its cells in std::map<Size, std::map<Size, ...>> (two tree nodes and several lookups per cell).
      Identity: the reference computes, per row i, one contiguous run of columns (from the left border
      it carries along to the column where it leaves the band), in the same order and with the same
      expressions as below. Every other cell it reads is either a border cell (i * tolerance or
      j * tolerance, (0,0) = 0) or absent, in which case it substitutes (i + j) * tolerance: the same
      value, so one formula serves both. A traceback step onto an absent cell reads a default (0,0)
      entry there, which ends the walk.
    */
    void alignAbsoluteTolerance_(std::vector<std::pair<Size, Size>>& alignment,
                                 const MSSpectrum& s1, const MSSpectrum& s2,
                                 const double tolerance, AlignmentScratch_& scratch)
    {
      if (!s1.isSorted() || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
      }
      alignment.clear();

      const Size n1 = s1.size();
      const Size n2 = s2.size();
      // cell (i, j) outside the computed runs: border value or the reference's substitute for an absent cell
      const auto outside = [tolerance](const Size i_plus_j) { return i_plus_j == 0 ? 0.0 : i_plus_j * tolerance; };
      enum : unsigned char { DIAGONAL = 1, UP = 2, LEFT = 3 }; // predecessor (i-1,j-1), (i,j-1), (i-1,j)

      std::vector<double>& prev = scratch.prev;
      std::vector<double>& cur = scratch.cur;
      std::vector<unsigned char>& dir = scratch.dir;
      prev.clear();
      dir.clear();
      scratch.row_begin.assign(n1 + 2, 0);
      scratch.row_offset.assign(n1 + 2, 0);

      Size left_ptr(1);
      Size last_i(0), last_j(0);
      Size prev_begin(1), prev_end(0); // computed columns of row i - 1 (row 0: none)
      for (Size i = 1; i <= n1; ++i)
      {
        const double pos1(s1[i - 1].getMZ());
        const Size cur_begin = left_ptr;
        cur.clear();
        scratch.row_begin[i] = cur_begin;
        scratch.row_offset[i] = dir.size();

        for (Size j = cur_begin; j <= n2; ++j)
        {
          bool off_band(false);
          const double pos2(s2[j - 1].getMZ());
          const double diff_align = fabs(pos1 - pos2);

          // running off the right border of the band?
          if (pos2 > pos1 && diff_align > tolerance)
          {
            if (i < n1 && j < n2 && s1[i].getMZ() < pos2)
            {
              off_band = true;
            }
          }

          // can we tighten the left border of the band?
          if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
          {
            ++left_ptr;
          }

          double score_align = diff_align;
          if (j - 1 >= prev_begin && j - 1 <= prev_end)
          {
            score_align += prev[j - 1 - prev_begin];
          }
          else
          {
            score_align += outside(i - 1 + j - 1);
          }

          double score_up = tolerance;
          if (j > cur_begin)
          {
            score_up += cur.back();
          }
          else
          {
            score_up += outside(i + j - 1);
          }

          double score_left = tolerance;
          if (j >= prev_begin && j <= prev_end)
          {
            score_left += prev[j - prev_begin];
          }
          else
          {
            score_left += outside(i - 1 + j);
          }

          if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
          {
            cur.push_back(score_align);
            dir.push_back(DIAGONAL);
            last_i = i;
            last_j = j;
          }
          else if (score_up <= score_left)
          {
            cur.push_back(score_up);
            dir.push_back(UP);
          }
          else
          {
            cur.push_back(score_left);
            dir.push_back(LEFT);
          }

          if (off_band)
          {
            break;
          }
        }
        prev_begin = cur_begin;
        prev_end = cur_begin + cur.size() - 1;
        prev.swap(cur);
      }
      scratch.row_offset[n1 + 1] = dir.size();

      // do traceback
      Size i = last_i;
      Size j = last_j;
      while (i >= 1 && j >= 1)
      {
        const Size begin = scratch.row_begin[i];
        const Size width = scratch.row_offset[i + 1] - scratch.row_offset[i];
        if (j < begin || j >= begin + width)
        {
          // absent cell: the reference reads (0,0), which counts as a diagonal step only at (1,1)
          if (i == 1 && j == 1) alignment.emplace_back(0, 0);
          break;
        }
        const unsigned char d = dir[scratch.row_offset[i] + (j - begin)];
        if (d == DIAGONAL)
        {
          alignment.emplace_back(i - 1, j - 1);
          --i;
          --j;
        }
        else if (d == UP)
        {
          --j;
        }
        else
        {
          --i;
        }
      }
      std::reverse(alignment.begin(), alignment.end());
    }
  }

  void ProSEAlgorithm::postProcessHits_(const PeakMap& exp,
        std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
        const std::vector<CandidatePoolStats_>& pool_stats,
        std::vector<ProteinIdentification>& protein_ids,
        PeptideIdentificationList& peptide_ids,
        Size top_hits,
  //      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
  //      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
  //      Size max_variable_mods_per_peptide, TODO: what about this parameter?
        const StringList& modifications_fixed,
        const StringList& modifications_variable,
        Int peptide_missed_cleavages,
        double precursor_mass_tolerance,
        double fragment_mass_tolerance,
        const std::string& precursor_mass_tolerance_unit_ppm,
        const std::string& fragment_mass_tolerance_unit_ppm,
        const Int precursor_min_charge,
        const Int precursor_max_charge,
        const std::string& enzyme,
        const std::string& database_name,
        const PeakMap* evidence_spectra) const
  {
    // Candidate-pool features (delta score, z-score, candidate count) are derived
    // from @p pool_stats rather than from @p annotated_hits: scoreSpectraAgainstIndex_
    // already pruned each spectrum to max(top_hits, 2) candidates, so the hits left
    // here are the report set, not the search space. pool_stats summarises every
    // candidate that was scored, zero-scoring ones included, and is accumulated
    // across chunks by the chunked search paths.
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }

    // One value per spectrum, kept in side vectors so scoring does not have to carry
    // them on every AnnotatedHit_ candidate.
    std::vector<float> delta_scores(annotated_hits.size(), 0.0f);
    std::vector<float> hyperscore_zscores(annotated_hits.size(), 0.0f);
    std::vector<float> ln_num_candidates(annotated_hits.size(), 0.0f);
#pragma omp parallel for default(none) shared(annotated_hits, top_hits, pool_stats, delta_scores, hyperscore_zscores, ln_num_candidates)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      const CandidatePoolStats_& stats = pool_stats[scan_index];
      delta_scores[scan_index] = static_cast<float>(stats.deltaScore());
      hyperscore_zscores[scan_index] = static_cast<float>(stats.zScore());
      // log1p rather than log: a spectrum can end up with zero candidates, and
      // ln(0) is not representable. The +1 offset is order-preserving, so this
      // stays a monotone proxy for the size of the scored search space.
      ln_num_candidates[scan_index] = static_cast<float>(std::log1p(static_cast<double>(stats.count)));

      // sort and keep n best elements according to score (unchanged)
      auto& hits = annotated_hits[scan_index];
      Size topn = top_hits > hits.size() ? hits.size() : top_hits;
      std::partial_sort(hits.begin(), hits.begin() + topn, hits.end(), AnnotatedHit_::hasBetterScore);
      hits.resize(topn);
      hits.shrink_to_fit();
    }

    bool annotation_precursor_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_fragment_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_prefix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_suffix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_num_matched_peaks = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::NUM_MATCHED_PEAKS) != annotate_psm_.end();
    bool annotation_matched_prefix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS) != annotate_psm_.end();
    bool annotation_matched_suffix_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS) != annotate_psm_.end();
    bool annotation_longest_ion_run = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE) != annotate_psm_.end();
    bool annotation_matched_ion_current = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT) != annotate_psm_.end();
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end();
    bool annotation_hyperscore_zscore = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::HYPERSCORE_ZSCORE) != annotate_psm_.end();
    bool annotation_ln_num_candidates = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LN_NUM_CANDIDATES) != annotate_psm_.end();
    bool annotation_matched_ion_current_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_ION_CURRENT_FRACTION) != annotate_psm_.end();
    bool annotation_complementary_ions_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::COMPLEMENTARY_IONS_FRACTION) != annotate_psm_.end();
    const bool annotation_local_evidence = param_.getValue("annotate:local_fragment_evidence").toBool();
    const bool evidence_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    if (annotation_local_evidence && (evidence_spectra == nullptr || evidence_spectra->size() != exp.size()))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Local fragment evidence requires the aligned peak lists retained during preprocessing.");
    }

    // "ALL" adds all annotations
    if (std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end())
    {
      annotation_precursor_error_ppm = true;
      annotation_fragment_error_ppm = true;
      annotation_prefix_fraction = true;
      annotation_suffix_fraction = true;
      annotation_num_matched_peaks = true;
      annotation_matched_prefix_ions = true;
      annotation_matched_suffix_ions = true;
      annotation_longest_ion_run = true;
      annotation_matched_ion_current = true;
      annotation_fragment_annotations = true;
      annotation_hyperscore_zscore = true;
      annotation_ln_num_candidates = true;
      annotation_matched_ion_current_fraction = true;
      annotation_complementary_ions_fraction = true;
    }

    // Alignment is needed for fragment error, fragment annotations, longest ion run, MIC,
    // normalized MIC, and complementary ion pairs
    const bool need_alignment = annotation_fragment_error_ppm || annotation_fragment_annotations || annotation_longest_ion_run
      || annotation_matched_ion_current || annotation_matched_ion_current_fraction || annotation_complementary_ions_fraction;

    // Both configurations depend only on the fragment tolerance, so they are
    // shared read-only by every thread of the loop below: getSpectrum() and
    // getSpectrumAlignment() are const and hold no mutable state (the scoring
    // loop in scoreSpectraAgainstIndex_ shares its generators the same way).
    // The alignment tolerance mirrors the search's fragment tolerance so the
    // reported FRAGMENT_ERROR_MEDIAN_PPM is not polluted by far-off spurious
    // matches — SpectrumAlignment's default is 0.3 Da absolute, which is ~30×
    // looser than a typical 20 ppm search.
    // The generators mirror the ion series of the scoring, so annotation features (fragment
    // error, peak annotations, longest ion run, MIC) are computed on the same theoretical
    // spectrum the candidate was scored against.
    const SpectrumGenerators_ generators = spectrumGenerators_();
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }
    // Da tolerance: same alignment from a flat-storage copy of SpectrumAlignment's banded DP
    // (alignAbsoluteTolerance_ above); ppm tolerance keeps SpectrumAlignment's cheap matching.
    const bool sa_absolute = !sa.getParameters().getValue("is_relative_tolerance").toBool();
    const double sa_tolerance = (double)sa.getParameters().getValue("tolerance");

    // Meta value keys of the loop below, resolved to registry indices once per call: every
    // string-keyed setMetaValue()/getMetaValue() takes the registry's process-wide lock (also for
    // names that are already registered), which serialised the annotation threads.
    // Identity: MetaInfo keeps (and idXML writes) the values of a hit ordered by registry index, so
    // new names must get their indices in the order in which the loop registered them: the first
    // spectrum with hits registers "scan_index", then IM if it carries a drift time, then its first
    // hit registers the per-hit names in the order they are set below. All hits set the same names,
    // so no thread could register a later name before an earlier one and this order did not depend
    // on the thread count. registerName() only looks up names that are already known, and without
    // any hit nothing was registered, so nothing is registered here either. Names that stay
    // string-keyed (spectrum reference, IM, rank; at most one lookup per spectrum or hit) are
    // still registered inside the loop, after the ones below, as before (checked: the registry
    // contents after a search are the same as before this change, at 1 and at 16 threads).
    const int isotope_error_sign = isotope_error_observed_minus_theoretical_ ? -1 : 1;
    const bool open_search_mode = isOpenSearchMode_();
    UInt mv_scan_index{}, mv_fragment_error{}, mv_precursor_error{}, mv_prefix_fraction{}, mv_suffix_fraction{},
         mv_num_matched_peaks{}, mv_matched_prefix_ions{}, mv_matched_suffix_ions{}, mv_delta_score{},
         mv_hyperscore_zscore{}, mv_ln_num_candidates{}, mv_matched_ion_current{}, mv_matched_ion_current_fraction{},
         mv_longest_ion_run{}, mv_complementary_ions_fraction{}, mv_chance_match_surprise{},
         mv_mass_competition_evidence{}, mv_isotope_error{}, mv_delta_mass{};
    {
      const auto first_with_hits = std::find_if(annotated_hits.begin(), annotated_hits.end(),
                                                [](const std::vector<AnnotatedHit_>& hits) { return !hits.empty(); });
      if (first_with_hits != annotated_hits.end())
      {
        MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
        mv_scan_index = registry.registerName("scan_index");
        if (IMTypes::determineIMFormat(exp[first_with_hits - annotated_hits.begin()]) == IMFormat::IM_SPECTRUM)
        {
          registry.registerName(Constants::UserParam::IM);
        }
        if (annotation_fragment_error_ppm) mv_fragment_error = registry.registerName(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
        if (annotation_precursor_error_ppm) mv_precursor_error = registry.registerName(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
        if (annotation_prefix_fraction) mv_prefix_fraction = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
        if (annotation_suffix_fraction) mv_suffix_fraction = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
        if (annotation_num_matched_peaks) mv_num_matched_peaks = registry.registerName(Constants::UserParam::NUM_MATCHED_PEAKS);
        if (annotation_matched_prefix_ions) mv_matched_prefix_ions = registry.registerName(Constants::UserParam::MATCHED_PREFIX_IONS);
        if (annotation_matched_suffix_ions) mv_matched_suffix_ions = registry.registerName(Constants::UserParam::MATCHED_SUFFIX_IONS);
        mv_delta_score = registry.registerName(Constants::UserParam::DELTA_SCORE);
        if (annotation_hyperscore_zscore) mv_hyperscore_zscore = registry.registerName(Constants::UserParam::HYPERSCORE_ZSCORE);
        if (annotation_ln_num_candidates) mv_ln_num_candidates = registry.registerName(Constants::UserParam::LN_NUM_CANDIDATES);
        if (annotation_matched_ion_current) mv_matched_ion_current = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT);
        if (annotation_matched_ion_current_fraction) mv_matched_ion_current_fraction = registry.registerName(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
        if (annotation_longest_ion_run) mv_longest_ion_run = registry.registerName(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
        if (annotation_complementary_ions_fraction) mv_complementary_ions_fraction = registry.registerName(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
        if (annotation_local_evidence)
        {
          mv_chance_match_surprise = registry.registerName(Constants::UserParam::CHANCE_MATCH_SURPRISE);
          mv_mass_competition_evidence = registry.registerName(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
        }
        mv_isotope_error = registry.registerName(Constants::UserParam::ISOTOPE_ERROR);
        if (open_search_mode) mv_delta_mass = registry.registerName("DeltaMass");
      }
    }

    // With an empty list to start from, every spectrum with hits gets a slot in scan order: the identifications
    // end up in scan order without a lock and without the sort below (which orders them by scan_index, i.e. the
    // same order). Otherwise they are appended under a lock and sorted, as before.
    const bool slots = peptide_ids.empty();
    std::vector<Size> slot_of_scan;
    if (slots)
    {
      slot_of_scan.resize(annotated_hits.size());
      Size used_slots = 0;
      for (Size scan_index = 0; scan_index < annotated_hits.size(); ++scan_index)
      {
        slot_of_scan[scan_index] = used_slots;
        if (!annotated_hits[scan_index].empty()) ++used_slots;
      }
      peptide_ids.resize(used_slots);
    }

    // The work per spectrum varies with its number of hits: hand out small blocks
#pragma omp parallel for schedule(dynamic, 16)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spec = exp[scan_index];
        const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setSpectrumReference( spec.getNativeID());
        pi.setMetaValue(mv_scan_index, static_cast<unsigned int>(scan_index));
        pi.setScoreType("ln(hyperscore)"); // also for scoring:method mass_accuracy, a log-space HyperScore (recorded in the search parameters)
        pi.setHigherScoreBetter(true);
        double mz = spec.getPrecursors()[0].getMZ();
        pi.setRT(spec.getRT());
        pi.setMZ(mz);

        // Annotate ion mobility if spectrum has a single drift time (DDA-PASEF)
        if (IMTypes::determineIMFormat(spec) == IMFormat::IM_SPECTRUM)
        {
          pi.setMetaValue(Constants::UserParam::IM, spec.getDriftTime());
        }

        Size charge = spec.getPrecursors()[0].getCharge();

        // Spectrum-level quantity, identical for every candidate of this spectrum, so it is
        // computed once here rather than per hit.
        const double spectrum_tic =
          annotation_matched_ion_current_fraction ? spec.calculateTIC() : 0.0;
        const MSSpectrum& evidence_spec = annotation_local_evidence ? (*evidence_spectra)[scan_index] : spec;
        const std::vector<double> local_densities = annotation_local_evidence
          ? localPeakDensities_(evidence_spec) : std::vector<double>{};

        AlignmentScratch_ alignment_scratch; // reused by all hits of this spectrum
        std::vector<double> evidence_alternatives; // reused by all hits of this spectrum

        // create full peptide hit structure from annotated hits
        vector<PeptideHit> phs;
        for (const auto& ah : annotated_hits[scan_index])
        {
          PeptideHit ph;
          // Prefer spectrum charge; if absent (0), fall back to the charge actually used by FI for this candidate
          const Size used_charge = (charge > 0) ? charge : static_cast<Size>(ah.applied_charge);
          ph.setCharge(used_charge);
          ph.setScore(ah.score);
          ph.setSequence(ah.sequence);

          // Generate theoretical spectrum + alignment for annotations that need it.
          std::vector<std::pair<Size, Size>> alignment;
          MSSpectrum theoretical_spec;
          // Annotate the charges actually scored when higher charges are scored.
          const int max_frag_z = scoring_multiple_charges_ ? scoringMaxCharge_(static_cast<int>(used_charge))
                                                           : ((charge >= 2) ? std::min<int>(charge - 1, 2) : 1);
          if (need_alignment)
          {
            tsg.getSpectrum(theoretical_spec, ah.sequence, 1, max_frag_z);
            if (sa_absolute)
            {
              alignAbsoluteTolerance_(alignment, theoretical_spec, spec, sa_tolerance, alignment_scratch);
            }
            else
            {
              sa.getSpectrumAlignment(alignment, theoretical_spec, spec);
            }
          }

          if (annotation_fragment_error_ppm)
          {
            std::vector<double> err;
            for (const auto& match : alignment)
            {
              double fragment_error = fabs(Math::getPPM(spec[match.second].getMZ(), theoretical_spec[match.first].getMZ()));
              err.push_back(fragment_error);
            }
            double median_ppm_error(0);
            if (!err.empty()) { median_ppm_error = Math::median(err.begin(), err.end(), false); }
            ph.setMetaValue(mv_fragment_error, median_ppm_error);
          }

          if (annotation_precursor_error_ppm)
          {
            // Subtract out the isotope offset FI matched at — FragmentIndex searches
            // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs
            // + isotope_error * C13C12, and the observed-to-monoiso correction in m/z is
            //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
            // (ah.isotope_error is this search offset, whatever sign the PSM reports).
            // Without this, a ±1 Da FI match reports ~1000 ppm / charge for the Percolator
            // feature, corrupting target/decoy discrimination.
            const double corrected_mz = mz
              + static_cast<double>(ah.isotope_error) * Constants::C13C12_MASSDIFF_U / used_charge;
            double theo_mz = ah.sequence.getMZ(used_charge);
            double ppm_difference = Math::getPPM(corrected_mz, theo_mz);
            ph.setMetaValue(mv_precursor_error, ppm_difference);
          }

          if (annotation_prefix_fraction)
          {
            ph.setMetaValue(mv_prefix_fraction, ah.prefix_fraction);
          }

          if (annotation_suffix_fraction)
          {
            ph.setMetaValue(mv_suffix_fraction, ah.suffix_fraction);
          }

          // Matched ion counts (from scoring, no alignment needed)
          if (annotation_num_matched_peaks)
          {
            ph.setMetaValue(mv_num_matched_peaks, static_cast<int>(ah.matched_prefix_ions + ah.matched_suffix_ions));
          }
          if (annotation_matched_prefix_ions)
          {
            ph.setMetaValue(mv_matched_prefix_ions, static_cast<int>(ah.matched_prefix_ions));
          }
          if (annotation_matched_suffix_ions)
          {
            ph.setMetaValue(mv_matched_suffix_ions, static_cast<int>(ah.matched_suffix_ions));
          }

          ph.setMetaValue(mv_delta_score, delta_scores[scan_index]);

          if (annotation_hyperscore_zscore)
          {
            ph.setMetaValue(mv_hyperscore_zscore, hyperscore_zscores[scan_index]);
          }
          if (annotation_ln_num_candidates)
          {
            ph.setMetaValue(mv_ln_num_candidates, ln_num_candidates[scan_index]);
          }

          // Fragment annotations, longest ion run, MIC, normalized MIC, and complementary
          // ion pairs all iterate the alignment + ion names
          if (annotation_fragment_annotations || annotation_longest_ion_run || annotation_matched_ion_current
            || annotation_matched_ion_current_fraction || annotation_complementary_ions_fraction)
          {
            const auto& ion_names = theoretical_spec.getStringDataArrays()[0];
            const auto& ion_charges = theoretical_spec.getIntegerDataArrays()[0];

            // Build PeakAnnotation vector + collect ion ordinals for longest run.
            // Prefix = a/b/c (N-terminal), suffix = x/y/z (C-terminal). Ordinals
            // for different ion types at the same cleavage position (e.g. a3 + b3)
            // are merged via std::unique below — each ordinal is a backbone
            // position, not an ion-type-specific identifier.
            std::vector<PeptideHit::PeakAnnotation> peak_annotations;
            std::vector<int> prefix_ordinals, suffix_ordinals;
            double matched_ion_current = 0.0;
            const bool need_mic = annotation_matched_ion_current || annotation_matched_ion_current_fraction;
            const bool need_ordinals = annotation_longest_ion_run || annotation_complementary_ions_fraction;
            // Dedup guard for MIC: in ppm-alignment mode a single experimental
            // peak can match multiple theoretical peaks (e.g. b-ion and near
            // isotope), so we must sum each exp_idx at most once. Sized only
            // when MIC is actually requested.
            std::vector<char> counted_exp_peaks(need_mic ? spec.size() : 0, 0);
            peak_annotations.reserve(alignment.size());

            for (const auto& [theo_idx, exp_idx] : alignment)
            {
              if (annotation_fragment_annotations)
              {
                PeptideHit::PeakAnnotation pa;
                pa.mz = spec[exp_idx].getMZ();
                pa.intensity = spec[exp_idx].getIntensity();
                pa.annotation = ion_names[theo_idx];
                pa.charge = ion_charges[theo_idx];
                peak_annotations.push_back(pa);
              }

              if (need_mic && !counted_exp_peaks[exp_idx])
              {
                matched_ion_current += spec[exp_idx].getIntensity();
                counted_exp_peaks[exp_idx] = 1;
              }

              if (need_ordinals && ion_names[theo_idx].size() >= 2)
              {
                const std::string& name = ion_names[theo_idx];
                const char c = name[0];
                const bool is_prefix = (c == 'a' || c == 'b' || c == 'c');
                const bool is_suffix = (c == 'x' || c == 'y' || c == 'z');
                if (is_prefix || is_suffix)
                {
                  // Extract ordinal: "b5", "y3-H2O1+", "c12++", "z.4+" (z+1) -> 5, 3, 12, 4
                  Size pos = 1;
                  while (pos < name.size() && (name[pos] == '.' || name[pos] == '\'')) ++pos; // z. (z+1), z' (z+2)
                  const Size ordinal_begin = pos;
                  while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
                  if (pos > ordinal_begin)
                  {
                    int ordinal = StringUtils::toInt32(StringUtils::substr(name, ordinal_begin, pos - ordinal_begin));
                    (is_prefix ? prefix_ordinals : suffix_ordinals).push_back(ordinal);
                  }
                }
              }
            }

            if (annotation_fragment_annotations)
            {
              ph.setPeakAnnotations(std::move(peak_annotations));
            }

            if (annotation_matched_ion_current)
            {
              ph.setMetaValue(mv_matched_ion_current, matched_ion_current);
            }

            if (annotation_matched_ion_current_fraction)
            {
              ph.setMetaValue(mv_matched_ion_current_fraction,
                              spectrum_tic > 0 ? matched_ion_current / spectrum_tic : 0.0);
            }

            if (need_ordinals)
            {
              // Compute longest consecutive run across prefix and suffix series.
              // Sorts + deduplicates each ordinal vector in place (a backbone position
              // matched by multiple ion types, e.g. a3 and b3, counts once).
              auto longestRun = [](std::vector<int>& v) -> int {
                if (v.empty()) return 0;
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                int best = 1, run = 1;
                for (Size i = 1; i < v.size(); ++i)
                {
                  if (v[i] == v[i - 1] + 1) { ++run; if (run > best) best = run; }
                  else run = 1;
                }
                return best;
              };
              int longest_prefix = longestRun(prefix_ordinals);
              int longest_suffix = longestRun(suffix_ordinals);

              if (annotation_longest_ion_run)
              {
                ph.setMetaValue(mv_longest_ion_run, std::max(longest_prefix, longest_suffix));
              }

              if (annotation_complementary_ions_fraction)
              {
                // A prefix ion at backbone position i (a_i/b_i/c_i) is complementary to a
                // suffix ion at position (peptide_length - i) (x/y/z at the same cleavage
                // site). Andromeda-style structural signal, distinct from the independently
                // computed prefix/suffix fractions above.
                const int pep_len = static_cast<int>(ah.sequence.size());
                double complementary_fraction = 0.0;
                if (pep_len > 1)
                {
                  Size n_complementary = 0;
                  for (int p : prefix_ordinals)
                  {
                    if (std::binary_search(suffix_ordinals.begin(), suffix_ordinals.end(), pep_len - p)) { ++n_complementary; }
                  }
                  complementary_fraction = static_cast<double>(n_complementary) / static_cast<double>(pep_len - 1);
                }
                ph.setMetaValue(mv_complementary_ions_fraction, complementary_fraction);
              }
            }
          }


          if (annotation_local_evidence)
          {
            // Use the fragment charges that were scored. The annotation spectrum (charges 1..max_frag_z)
            // contains them whenever it was generated, so it is reused instead of generating a second one.
            const int evidence_z = scoringMaxCharge_(static_cast<int>(used_charge));
            LocalFragmentEvidence_ evidence;
            if (need_alignment && evidence_z <= max_frag_z)
            {
              evidence = localFragmentEvidence_(evidence_spec, theoretical_spec, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            else
            {
              MSSpectrum evidence_theory;
              tsg.getSpectrum(evidence_theory, ah.sequence, 1, evidence_z);
              evidence = localFragmentEvidence_(evidence_spec, evidence_theory, evidence_z, local_densities,
                fragment_mass_tolerance, evidence_ppm, evidence_alternatives);
            }
            ph.setMetaValue(mv_chance_match_surprise, evidence.chance_match_surprise);
            ph.setMetaValue(mv_mass_competition_evidence, evidence.mass_competition_evidence);
          }

          // Add isotope error metavalue (always; exposed as Percolator feature). ah.isotope_error is the offset
          // FragmentIndex added to the observed mass (theoretical minus observed); report:isotope_error_convention
          // selects the reported sign.
          ph.setMetaValue(mv_isotope_error, isotope_error_sign * ah.isotope_error);

          // Add delta mass metavalue for open search
          if (open_search_mode)
          {
            ph.setMetaValue(mv_delta_mass, ah.delta_mass);
          }

          // store PSM
          phs.push_back(std::move(ph));
        }
        pi.setHits(std::move(phs));
        // Ensure hits are sorted by score (best first), then assign ranks explicitly (0 = top hit)
        pi.sort();
        {
          std::vector<PeptideHit>& hits = pi.getHits();
          for (Size r = 0; r < hits.size(); ++r)
          {
            hits[r].setRank(static_cast<UInt>(r));
          }
        }

        // Debug: log spectrum-level top hit details before storing PeptideIdentification.
        // DEBUG (not INFO) because multi-file mode would emit one line per scan per input.
        if (!pi.getHits().empty())
        {
          const PeptideHit& top_hit = pi.getHits().front();
          OPENMS_LOG_DEBUG << "[ProSE] scan_index=" << scan_index
                           << " top_ln(hyperscore)=" << top_hit.getScore()
                           << " top_charge=" << top_hit.getCharge()
                           << " top_isotope_error=" << (int)top_hit.getMetaValue(mv_isotope_error)
                           << std::endl;
        }
        if (slots)
        {
          peptide_ids[slot_of_scan[scan_index]] = std::move(pi);
        }
        else
        {
#pragma omp critical (peptide_ids_access)
          {
            //clang-tidy: seems to be a false-positive in combination with omp
            peptide_ids.push_back(std::move(pi));
          }
        }
      }
    }

#ifdef _OPENMP
    // we need to sort the peptide_ids by scan_index in order to have the same output in the idXML-file
    if (omp_get_max_threads() > 1 && !slots)
    {
      // one registry lookup for the whole sort instead of two (locked) lookups per comparison;
      // getMetaValue(name) is getMetaValue(getIndex(name)), so the comparisons are unchanged
      const UInt scan_index_key = MetaInfoInterface::metaRegistry().getIndex("scan_index");
      std::sort(peptide_ids.begin(), peptide_ids.end(), [scan_index_key](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return a.getMetaValue(scan_index_key) < b.getMetaValue(scan_index_key);
      });
    }
#endif

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("ProSE");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());

    DateTime now = DateTime::now();
    std::string identifier("ProSE_" + now.get());
    protein_ids[0].setIdentifier(identifier);
    for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges =StringUtils::toStr(precursor_min_charge) + ":" + StringUtils::toStr(precursor_max_charge);

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::PeakMassType::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = modifications_fixed;
    search_parameters.variable_modifications = modifications_variable;
    search_parameters.missed_cleavages = peptide_missed_cleavages;
    search_parameters.fragment_mass_tolerance = fragment_mass_tolerance;
    search_parameters.precursor_mass_tolerance = precursor_mass_tolerance;
    search_parameters.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm == "ppm";
    search_parameters.fragment_mass_tolerance_ppm = fragment_mass_tolerance_unit_ppm == "ppm";
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);

    // add additional percolator features or post-processing
    StringList feature_set{"score"};
    if (annotation_fragment_error_ppm) feature_set.push_back(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
    if (annotation_prefix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
    if (annotation_suffix_fraction) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
    if (annotation_longest_ion_run) feature_set.push_back(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE);
    if (annotation_matched_prefix_ions) feature_set.push_back(Constants::UserParam::MATCHED_PREFIX_IONS);
    if (annotation_matched_suffix_ions) feature_set.push_back(Constants::UserParam::MATCHED_SUFFIX_IONS);
    if (annotation_matched_ion_current) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT);
    if (annotation_matched_ion_current_fraction) feature_set.push_back(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION);
    if (annotation_complementary_ions_fraction) feature_set.push_back(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION);
    if (annotation_hyperscore_zscore) feature_set.push_back(Constants::UserParam::HYPERSCORE_ZSCORE);
    if (annotation_ln_num_candidates) feature_set.push_back(Constants::UserParam::LN_NUM_CANDIDATES);
    if (annotation_local_evidence)
    {
      feature_set.push_back(Constants::UserParam::CHANCE_MATCH_SURPRISE);
      feature_set.push_back(Constants::UserParam::MASS_COMPETITION_EVIDENCE);
    }
    feature_set.push_back(Constants::UserParam::DELTA_SCORE);
    feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
    // note: precursor error is calculated by percolator itself
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    // Readers tell the sign of the PSMs' isotope_error by this record; without it (earlier ProSE versions, or
    // report:isotope_error_convention=theoretical_minus_observed) the sign is theoretical minus observed.
    if (isotope_error_observed_minus_theoretical_)
    {
      search_parameters.setMetaValue("isotope_error_convention", "observed_minus_theoretical");
    }
    // record whether open-search mode was used
    search_parameters.setMetaValue("open_search", isOpenSearchMode_() ? "true" : "false");

    search_parameters.setMetaValue("peptide:clip_nterm_methionine", param_.getValue("peptide:clip_nterm_methionine"));
    search_parameters.setMetaValue("peptide:deduplicate", param_.getValue("peptide:deduplicate"));
    search_parameters.setMetaValue("fragment:query_spectrum", param_.getValue("fragment:query_spectrum"));
    search_parameters.setMetaValue("fragment:query_spectrum_resolved", query_raw_spectrum_ ? "raw" : "processed");
    search_parameters.setMetaValue("annotate:local_fragment_evidence", param_.getValue("annotate:local_fragment_evidence"));
    search_parameters.setMetaValue("scoring:fragment_charges", param_.getValue("scoring:fragment_charges"));
    search_parameters.setMetaValue("scoring:fragment_charges_resolved", scoring_multiple_charges_ ? "multiple" : "single");
    if (mass_accuracy_score_) // recorded when it changes the score, so HyperScore searches keep their parameter list
    {
      search_parameters.setMetaValue("scoring:method", param_.getValue("scoring:method"));
      search_parameters.setMetaValue("scoring:method_resolved", "mass_accuracy");
      search_parameters.setMetaValue("scoring:mass_error_sd", mass_error_sd_ppm_);
    }
    search_parameters.setMetaValue("peaks:window_type", peaks_window_type_);
    search_parameters.setMetaValue(
      "peaks:window_type_resolved",
      peaks_window_type_ == "auto"
        ? (Deisotoper::isToleranceSupported(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm") ? "jump_full" : "jump")
        : peaks_window_type_);
    // MS2 deisotoping rule (used only when the spectra are deisotoped, see fragment:deisotope)
    search_parameters.setMetaValue("fragment:deisotope_min_peaks", static_cast<int>(deisotoping_.min_peaks));
    search_parameters.setMetaValue("fragment:deisotope_charge_cap", deisotoping_.charge_cap_precursor ? "precursor" : "none");
    search_parameters.setMetaValue("fragment:deisotope_sum_intensity", deisotoping_.sum_intensity ? "true" : "false");

    search_parameters.enzyme_term_specificity = peptide_enzyme_specificity_;
    protein_ids[0].setSearchParameters(std::move(search_parameters));

    // Annotate IM unit on ProteinIdentification if all PeptideIdentifications have IM
    if (!peptide_ids.empty())
    {
      bool all_have_im = std::all_of(
        peptide_ids.begin(), peptide_ids.end(),
        [](const PeptideIdentification& pid) { return pid.metaValueExists(Constants::UserParam::IM); });

      if (all_have_im)
      {
        protein_ids[0].setMetaValue(Constants::UserParam::IM,
          exp[0].getDriftTimeUnitAsString());
      }
    }
  }

  // =====================================================================
  // Build a SearchContext (decoys + FragmentIndex) from a FASTA database.
  // Hoisted out of the in-memory search() body so callers can build the
  // index once and reuse it across many spectrum files.
  // =====================================================================
  Param ProSEAlgorithm::fragmentIndexParameters_(bool electron_ions) const
  {
    Param p = getParameters();
    // FragmentIndex has its own boolean 'decoys' flag {true,false}; ProSE's enum
    // value (auto/generate/ignore) would fail its validation. ProSE builds the
    // decoy-augmented database itself, so the FragmentIndex must never generate.
    p.remove("decoys");
    p.setValue("decoys", "false");
    // ions:by_activation is resolved by the caller, which knows the spectra. The c and z+1 ions
    // for electron-activated spectra go into a set of their own, which the other spectra are not
    // matched against (see scoreSpectraAgainstIndex_).
    p.remove("ions:by_activation");
    p.setValue("ions:electron_ions", electron_ions ? "true" : "false");
    // annotation-only parameters of ProSE (see annotateIonPriors_)
    p.remove("annotate:self_trained_ion_priors");
    p.removeAll("annotate:ion_prior_");
    p.remove("peptide:protein_mapping"); // ProSE's own (FragmentIndex::getProteinOccurrences() serves it)
    return p;
  }

  namespace
  {
    /**
      @brief DecoyHelper::findDecoyString(@p db, quiet = true), with the OpenMP threads for large databases.

      DecoyHelper::countDecoys() sums, over the proteins, the matches of the prefix and the suffix pattern and keeps,
      per match, the spelling of its last occurrence (within a protein the suffix after the prefix). Chunks of
      consecutive proteins counted on their own and joined in their order give the same. The decision is the one of
      findDecoyString() for these statistics; it does not depend on the order in which it looks at the matches,
      since at most one prefix and one suffix can reach 80% of all prefixes or suffixes.

      At most 8 threads: the step (about 60 ms serially for 50,000 entries) gains little beyond 4 threads, and a
      larger team costs CPU time (measured at 64 threads: 2-4 CPU-s per run for no shorter wall time).
    */
    DecoyHelper::Result findDecoyStringInParallel(const std::vector<FASTAFile::FASTAEntry>& db)
    {
#ifdef _OPENMP
      constexpr int max_threads = 8;
      const int threads = omp_in_parallel() ? 1 : std::min(omp_get_max_threads(), max_threads);
#else
      const int threads = 1;
#endif
      if (threads < 2 || db.size() < 8192)
      {
        FASTAContainer<TFI_Vector> container(db);
        return DecoyHelper::findDecoyString(container, /*quiet=*/true);
      }
      struct Counts
      {
        std::map<std::string, std::pair<Size, Size>> count; ///< match -> occurrences as prefix, as suffix
        std::map<std::string, std::string> spelling;        ///< match -> spelling of its last occurrence
        Size prefixes = 0, suffixes = 0;
      };
      const SignedSize n = static_cast<SignedSize>(db.size());
      const SignedSize num_chunks = 4 * static_cast<SignedSize>(threads);
      std::vector<Counts> chunks(num_chunks);
#pragma omp parallel num_threads(threads)
      {
        const RegularExpression prefix_pattern(DecoyHelper::regexstr_prefix);
        const RegularExpression suffix_pattern(DecoyHelper::regexstr_suffix);
        std::string lower, match;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize c = 0; c < num_chunks; ++c)
        {
          Counts& chunk = chunks[c];
          for (SignedSize i = n * c / num_chunks; i < n * (c + 1) / num_chunks; ++i)
          {
            const std::string& identifier = db[i].identifier;
            lower = identifier;
            StringUtils::toLower(lower);
            if (prefix_pattern.search(lower, &match))
            {
              ++chunk.prefixes;
              ++chunk.count[match].first;
              chunk.spelling[match] = StringUtils::prefix(identifier, match.length());
            }
            if (suffix_pattern.search(lower, &match))
            {
              ++chunk.suffixes;
              ++chunk.count[match].second;
              chunk.spelling[match] = StringUtils::suffix(identifier, match.length());
            }
          }
        }
      }
      Counts all;
      for (const Counts& chunk : chunks) // in order: later chunks hold later proteins
      {
        all.prefixes += chunk.prefixes;
        all.suffixes += chunk.suffixes;
        for (const auto& [match, count] : chunk.count)
        {
          all.count[match].first += count.first;
          all.count[match].second += count.second;
        }
        for (const auto& [match, spelling] : chunk.spelling) all.spelling[match] = spelling;
      }
      // DecoyHelper::findDecoyString()'s decision
      const double proteins = static_cast<double>(db.size());
      if (static_cast<double>(all.prefixes + all.suffixes) < 0.4 * proteins || all.prefixes == all.suffixes)
      {
        return {false, "?", true};
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.first) / static_cast<double>(all.prefixes) >= 0.8 && static_cast<double>(count.first) / proteins >= 0.4)
        {
          return {true, all.spelling[match], true};
        }
      }
      for (const auto& [match, count] : all.count)
      {
        if (static_cast<double>(count.second) / static_cast<double>(all.suffixes) >= 0.8 && static_cast<double>(count.second) / proteins >= 0.4)
        {
          return {true, all.spelling[match], false};
        }
      }
      return {false, "?", true};
    }
  }

  ProSEAlgorithm::DecoyStrategy_
  ProSEAlgorithm::resolveDecoyStrategy_(const std::vector<FASTAFile::FASTAEntry>& db) const
  {
    // Detect pre-existing decoys. First the common-marker heuristic
    // (DecoyHelper: "decoy", "rev", "xxx", ... as prefix or suffix), then a
    // literal fall-back to the configured decoy_prefix so custom markers
    // outside the vocabulary are still recognised.
    // quiet=true: a target-only database is a normal case here (auto/generate
    // then synthesise decoys), so suppress DecoyHelper's "unable to determine
    // decoy string" ERROR/WARN noise — we handle the negative result ourselves.
    const DecoyHelper::Result det = findDecoyStringInParallel(db);

    bool existing = det.success;
    std::string ext_string = det.success ? det.name : decoy_prefix_;
    bool ext_is_prefix = det.success ? det.is_prefix : true;

    if (!existing && !decoy_prefix_.empty())
    {
      bool any_prefix = false, any_suffix = false;
      for (const FASTAFile::FASTAEntry& e : db)
      {
        if (StringUtils::hasPrefix(e.identifier, decoy_prefix_)) { any_prefix = true; break; }
        if (StringUtils::hasSuffix(e.identifier, decoy_prefix_)) { any_suffix = true; }
      }
      if (any_prefix || any_suffix)
      {
        existing = true;
        ext_string = decoy_prefix_;
        ext_is_prefix = any_prefix; // prefer prefix when both orientations occur
      }
    }

    DecoyStrategy_ s;
    switch (decoy_mode_)
    {
      case DecoyMode_::AUTO:
        if (existing)
        {
          // Reuse the decoys already in the database; never double-generate.
          s.generate = false; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = ext_string; s.is_prefix = ext_is_prefix;
        }
        else
        {
          // No decoys present: synthesise them so target-decoy FDR is possible.
          s.generate = true; s.strip_existing = false; s.have_decoys = true;
          s.decoy_string = decoy_prefix_; s.is_prefix = true;
        }
        break;

      case DecoyMode_::GENERATE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=generate: removing pre-existing decoys (marker '"
                          << ext_string << "') and regenerating from the target proteins.\n";
        }
        else
        {
          // Nothing detected to strip. If the database still carries some decoy-like accessions
          // with a non-standard or too-sparse marker (below DecoyHelper's threshold, and not the
          // configured decoy_prefix), they are treated as targets and reversed into decoys-of-decoys.
          // Warn so the user can set -Search:decoy_prefix to the real marker or pre-clean the FASTA.
          FASTAContainer<TFI_Vector> stats_container(db);
          const DecoyHelper::DecoyStatistics stats = DecoyHelper::countDecoys(stats_container);
          if (stats.all_prefix_occur + stats.all_suffix_occur > 0)
          {
            OPENMS_LOG_WARN << "[ProSE] decoys=generate: " << (stats.all_prefix_occur + stats.all_suffix_occur)
                            << " accession(s) carry a decoy-like marker but too few to auto-detect, so "
                            << "they are not stripped and will be treated as targets (risking "
                            << "decoys-of-decoys). Set -Search:decoy_prefix to the actual marker or "
                            << "pre-clean the database.\n";
          }
        }
        s.generate = true; s.strip_existing = existing; s.have_decoys = true;
        s.decoy_string = decoy_prefix_; s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;

      case DecoyMode_::IGNORE:
        if (existing)
        {
          OPENMS_LOG_WARN << "[ProSE] decoys=ignore: removing pre-existing decoys (marker '"
                          << ext_string << "'); searching the target proteins only.\n";
        }
        s.generate = false; s.strip_existing = existing; s.have_decoys = false;
        s.decoy_string.clear(); s.is_prefix = true;
        s.strip_string = ext_string; s.strip_is_prefix = ext_is_prefix;
        break;
    }
    return s;
  }

  // Build the searched database according to @p strategy: optionally strip
  // pre-existing decoys, optionally generate fresh decoys from the targets.
  // Produces FASTA entries only — does not build a FragmentIndex.
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    return buildDecoyAugmentedDB_(std::vector<FASTAFile::FASTAEntry>(fasta_db), strategy);
  }

  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db,
      const DecoyStrategy_& strategy) const
  {
    // The entries are filtered in place instead of being copied one by one: the same entries
    // in the same order.
    std::vector<FASTAFile::FASTAEntry> db = std::move(fasta_db);

    // decoys=auto reusing pre-existing decoys logs nothing otherwise; surface the auto-detected
    // marker so a rare DecoyHelper misdetection (a target DB whose accessions start with a decoy
    // affix) is diagnosable. Single emission: buildDecoyAugmentedDB_ runs once per search.
    const bool reuses_decoys = decoy_mode_ == DecoyMode_::AUTO && !strategy.generate && strategy.have_decoys;
    if (reuses_decoys)
    {
      OPENMS_LOG_INFO << "[ProSE] decoys=auto: reusing existing decoys detected in the database "
                      << "(marker '" << strategy.decoy_string << "', "
                      << (strategy.is_prefix ? "prefix" : "suffix") << ")." << std::endl;
    }
    // Initial-Met clipping adds N-terminal peptides of every protein that starts with M. Generated decoys keep the
    // initial Met (below); supplied decoys often do not (a reversed protein ends with it), and then the clipped
    // peptides enlarge the target space only. Count both kinds while the entries are filtered.
    const bool check_met_symmetry = reuses_decoys && param_.getValue("peptide:clip_nterm_methionine").toBool();
    Size targets = 0, targets_with_met = 0, decoys = 0, decoys_with_met = 0;

    // 1. Keep targets, dropping pre-existing decoys when requested. A stop codon ('*') that ends
    //    a sequence, as in databases translated from genomes (e.g. SGD), is not a residue: remove
    //    it so the C-terminal peptide stays searchable and decoys are built from the protein
    //    alone. FragmentIndex skips peptides that contain a stop codon inside the sequence.
    //    An entry left without residues has nothing to search, and decoy generation needs
    //    residues: drop it.
    auto kept = db.begin();
    for (FASTAFile::FASTAEntry& e : db)
    {
      const bool is_existing_decoy = strategy.strip_existing &&
        (strategy.strip_is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.strip_string)
                                  : StringUtils::hasSuffix(e.identifier, strategy.strip_string));
      if (!is_existing_decoy)
      {
        std::string& sequence = e.sequence;
        while (!sequence.empty() && sequence.back() == '*') { sequence.pop_back(); }
        if (!sequence.empty())
        {
          if (check_met_symmetry)
          {
            const bool is_decoy = strategy.is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.decoy_string)
                                                     : StringUtils::hasSuffix(e.identifier, strategy.decoy_string);
            (is_decoy ? decoys : targets) += 1;
            (is_decoy ? decoys_with_met : targets_with_met) += (sequence.size() > 1 && sequence[0] == 'M') ? 1 : 0;
          }
          if (&*kept != &e) { *kept = std::move(e); }
          ++kept;
        }
      }
    }
    db.erase(kept, db.end());
    // Warn when the share of decoys that start with M is below half the share of targets that do.
    if (check_met_symmetry && targets_with_met > 0 && 2 * decoys_with_met * targets < targets_with_met * decoys)
    {
      OPENMS_LOG_WARN << "[ProSE] peptide:clip_nterm_methionine: " << targets_with_met << " of " << targets
                      << " target proteins but only " << decoys_with_met << " of " << decoys << " decoy proteins "
                      << "in the database start with M. Removing the initial Met adds N-terminal peptides to the "
                      << "targets that have no decoy counterpart, which makes target-decoy FDR estimates slightly "
                      << "optimistic. Use '-Search:decoys generate' (generated decoys keep the initial Met) or "
                      << "'-Search:peptide:clip_nterm_methionine false' for a symmetric search space." << std::endl;
    }

    // 2. Generate decoys by reversing the (remaining) target proteins.
    if (strategy.generate)
    {
      // decoy_string is the prefix prepended to generated accessions; an empty prefix would
      // produce decoys with the same accession as their targets, silently breaking FDR.
      if (strategy.decoy_string.empty())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "decoy_prefix must be non-empty to generate decoys (decoys=auto without existing decoys, "
          "or decoys=generate).");
      }
      DecoyGenerator decoy_generator;
      const size_t old_size = db.size();
      db.reserve(old_size * 2);
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = db[i];
        // Keep an initial Met on generated decoys as well: otherwise clipping
        // expands the target N-terminal search space without its decoy counterpart.
        const bool preserve_met = param_.getValue("peptide:clip_nterm_methionine").toBool() && e.sequence.size() > 1 && e.sequence[0] == 'M';
        if (preserve_met) { e.sequence.erase(0, 1); }
        if (peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE)
          e.sequence = decoy_generator.reverseProtein(AASequence::fromString(e.sequence)).toString();
        else
          e.sequence = decoy_generator.reversePeptides(AASequence::fromString(e.sequence), enzyme_).toString();
        if (preserve_met) { e.sequence.insert(e.sequence.begin(), 'M'); }
        e.identifier = strategy.decoy_string + e.identifier;  // decoy_string is the prefix to add
        db.push_back(std::move(e));
      }
      Math::RandomShuffler shuffler(42);  // fixed seed for reproducible decoy ordering across runs/files
      shuffler.portable_random_shuffle(db.begin(), db.end());

      // Under the default decoys=auto, decoys are generated whenever the database lacks them —
      // including no-ProSE-FDR runs, because the decoys are typically consumed by a downstream or
      // external FDR step (FalseDiscoveryRate / Percolator) or a merged run. This roughly doubles
      // index build + search; for a genuinely target-only search, use decoys=ignore.
      if (decoy_mode_ == DecoyMode_::AUTO && fdr_psm_ == 0.0 && fdr_protein_ == 0.0)
      {
        OPENMS_LOG_INFO << "[ProSE] decoys=auto generated " << old_size << " decoys (the database had "
                        << "none) for downstream target-decoy FDR. Use '-Search:decoys ignore' for a "
                        << "target-only search if you do not need FDR." << std::endl;
      }
    }
    return db;
  }

  // =====================================================================
  // Strided calibration sample. Used by the chunked search paths so that
  // calibration always sees a representative, suitably-sized subset of the
  // full DB — independent of database_chunk_size_. Prevents the first-chunk
  // starvation mode where small chunk_size values (e.g. 500) produce a
  // calibration pool below calibration_min_psms_ and silently disable
  // calibration. See issue #9182.
  // =====================================================================
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildCalibrationSample_(
      const std::vector<FASTAFile::FASTAEntry>& full_db) const
  {
    // Calibration FI size is bounded by the user's memory budget — which is
    // exactly what they declared via database_chunk_size_. This keeps peak RSS
    // during calibration identical to peak RSS during the main search
    // (one chunk FI at a time), regardless of digestion mode.
    //
    // Critical for immunopeptidomics (non-specific digestion, MHC-I/II variant
    // databases): a naive fixed-size 5000-protein cal sample can generate
    // 50+ GB of fragment index because each protein produces thousands of
    // candidate peptides — defeating the memory savings chunking is designed
    // to deliver. Tying cal_target to chunk_size keeps the promise that
    // "user's memory budget = one chunk's FI."
    //
    // Strided (not contiguous / random) to give the calibration pool even
    // coverage of the full DB.
    //
    // Residual gap vs non-chunked calibration:
    //   - Trypsin / moderate chunk_size (e.g. 5000): cal sample is adequate;
    //     yield is within ~2% of non-chunked.
    //   - Small chunk_size (e.g. 100 proteins for MHC-II variants): cal pool
    //     can fall below calibration_min_psms_ → runCalibrationPass_ returns
    //     success=false, calibration silently no-ops. The main search uses
    //     user-configured tolerances.
    // Follow-up issue: pool calibration PSMs across the first K chunks of the
    // main search loop, so calibration quality becomes independent of
    // chunk_size.
    const Size cal_target = std::min(
        std::max<Size>(database_chunk_size_, Size(1)),
        full_db.size());
    if (cal_target >= full_db.size()) return full_db;

    const Size stride = std::max<Size>(1, full_db.size() / cal_target);
    std::vector<FASTAFile::FASTAEntry> cal_db;
    cal_db.reserve(cal_target);
    for (Size i = 0; i < full_db.size() && cal_db.size() < cal_target; i += stride)
    {
      cal_db.push_back(full_db[i]);
    }
    return cal_db;
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db) const
  {
    return prepareContext(fasta_db, false);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db, bool electron_ions) const
  {
    return prepareContext_(std::vector<FASTAFile::FASTAEntry>(fasta_db), electron_ions);
  }

  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext_(
      std::vector<FASTAFile::FASTAEntry>&& fasta_db, bool electron_ions,
      const std::function<const PeakMap*(Size)>& searched_spectra) const
  {
    SearchContext ctx;

    startProgress(0, 1, "Generate decoys...");
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    ctx.db = buildDecoyAugmentedDB_(std::move(fasta_db), strategy);
    ctx.decoy_string = strategy.decoy_string;
    ctx.decoy_is_prefix = strategy.is_prefix;
    ctx.have_decoys = strategy.have_decoys;
    endProgress();

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    Param this_params = fragmentIndexParameters_(electron_ions);
    ctx.fragment_index.setParameters(this_params);
    ctx.fragment_index.build(ctx.db, searched_spectra);
    ctx.electron_ions = electron_ions;
    endProgress();

    return ctx;
  }

  // =====================================================================
  // Score all spectra against a single FragmentIndex, appending results to
  // annotated_hits. Used by both chunked and non-chunked search paths.
  // =====================================================================
  void ProSEAlgorithm::scoreSpectraAgainstIndex_(
      const PeakMap& spectra,
      FragmentIndex& fi,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const SpectrumGenerators_& generators,
      double effective_fragment_tol,
      bool fragment_mass_tolerance_unit_ppm,
      bool open_search_mode,
      std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
      std::vector<CandidatePoolStats_>& pool_stats,
      const std::string& progress_label,
      const PeakMap* query_spectra) const
  {
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }
    if (query_spectra != nullptr && query_spectra->size() != spectra.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, query_spectra->size(),
        "Candidate-retrieval spectra must hold one entry per scoring spectrum.");
    }
    startProgress(0, spectra.size(), progress_label);
    size_t count_spectra{};
    const double proton_mass_u = Constants::PROTON_MASS_U;
    // Hoisted out of the omp parallel block: clang with `default(none)` forbids
    // referencing namespace-scope constants inside the loop without explicit sharing.
    const double c13c12_massdiff_u = Constants::C13C12_MASSDIFF_U;
    const Size keep = std::max(report_top_hits_, Size(2)); // keep ≥2 for delta score
    const bool deduplicate_chunks = param_.getValue("peptide:deduplicate").toBool() && database_chunk_size_ > 0 && ! fi.isSnesMode();

    // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
    // leave the parallel region, which would terminate the process: the first one is rethrown after the loop.
    std::exception_ptr scoring_error;
    std::atomic<bool> scoring_failed{false};
#pragma omp parallel for schedule(dynamic) default(none) shared(annotated_hits, pool_stats, query_spectra, count_spectra, fi, generators, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, c13c12_massdiff_u, effective_fragment_tol, keep, deduplicate_chunks, scoring_error, scoring_failed)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {
      if (scoring_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      #pragma omp atomic
      ++count_spectra;

      IF_MASTERTHREAD { setProgress(count_spectra); }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      const TheoreticalSpectrumGenerator& spectrum_generator = generators.forSpectrum(exp_spectrum);
      FragmentIndex::SpectrumMatchesTopN top_sms;
      // ions:by_activation: only electron-activated spectra are matched against the c and z+1 ions,
      // so that these ions do not change which candidates the other spectra keep
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_index] : exp_spectrum;
      fi.querySpectrum(query, db, top_sms, generators.electronIons(exp_spectrum));

      const bool snes_mode = fi.isSnesMode();
      const bool prec_tol_ppm = precursor_mass_tolerance_unit_ == "ppm";
      // SNES realization uses the asymmetric precursor tolerance — same signed
      // window the FragmentIndex candidate filter used at bin-walk time.
      // Previously collapsed to max(lower, upper), which over-admitted by ~20×
      // on calibrated asymmetric configs like [100 ppm, 5 ppm]. Review L3.
      const double snes_realize_tol_lo = precursor_mass_tolerance_lower_;
      const double snes_realize_tol_hi = precursor_mass_tolerance_upper_;

      // Reused across candidates of this spectrum. Avoids per-candidate heap
      // churn of a fresh PeakSpectrum + its DataArrays (TSG's add_metainfo fills
      // StringDataArrays with ion names — these are a notable allocation hot spot
      // when the candidate count per spectrum is in the tens/hundreds).
      PeakSpectrum theo_spectrum;

      // SNES-mode dedup: a sub-peptide [i, i+k) can be produced by both a
      // Single-N mother anchored at i and a Single-C mother ending at i+k-1.
      // Both realize to the same AASequence (same protein, start, length,
      // variable-mod subset). Without this guard, both are scored and land
      // in annotated_hits, inflating the candidate list and biasing delta
      // scores / Percolator features. Per-spectrum state — cheap, bounded
      // by max_candidates_per_spectrum. Empty for non-SNES queries.
      // Key: (protein_idx, realized_start, realized_length, subset_bitmask).
      std::set<std::tuple<UInt32, uint16_t, uint16_t, uint32_t>> seen_realizations;

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& sms_pep = fi.getPeptides()[sms.peptide_idx_];

        AASequence mod_candidate;
        if (snes_mode)
        {
          // Realize the sub-peptide at the length whose mass best matches the
          // observed precursor (iso-corrected per the FI's shifted-mass convention).
          const double exp_mz = exp_spectrum.getPrecursors()[0].getMZ();
          const double observed_mh_plus =
              exp_mz * sms.precursor_charge_ - (sms.precursor_charge_ - 1) * proton_mass_u;
          // SNES v1.1: subtract the variable-mod Σ from the realization target
          // so realizeSNESLength compares against the *unmodified* realized mass.
          // For v1 (unmodified) hits, sms.sigma_delta_ == 0 — same semantics as before.
          const double iso_shifted_target = observed_mh_plus
              + static_cast<double>(sms.isotope_error_) * c13c12_massdiff_u
              - static_cast<double>(sms.sigma_delta_);
          const int realized_len = fi.realizeSNESLength(
              sms_pep, db, iso_shifted_target,
              snes_realize_tol_lo, snes_realize_tol_hi, prec_tol_ppm);
          if (realized_len < 0) continue; // no realizable length within tolerance

          // L2 dedup: skip if this exact realization was already scored via
          // the opposite-kind mother. Same protein + start + length + subset
          // → same AASequence → same score, redundant work and inflated hits.
          const uint16_t realized_start = FragmentIndex::isSingleCMother(sms_pep.mod_bitmask_)
              ? static_cast<uint16_t>(sms_pep.sequence_.first + sms_pep.sequence_.second
                                       - static_cast<uint16_t>(realized_len))
              : sms_pep.sequence_.first;
          const auto key = std::make_tuple(sms_pep.protein_idx, realized_start,
                                            static_cast<uint16_t>(realized_len),
                                            sms.subset_bitmask_);
          if (!seen_realizations.insert(key).second) continue;

          mod_candidate = fi.reconstructRealizedSubSequence(
              sms_pep, db, static_cast<size_t>(realized_len), sms.subset_bitmask_);
        }
        else
        {
          mod_candidate = fi.reconstructModifiedSequence(sms_pep, db);
        }

        // Index construction already removes repeated occurrences within a
        // chunk. Across chunks, skip the same hypothesis before updating either
        // scores or pool statistics. Charge and isotope hypotheses stay distinct.
        if (deduplicate_chunks)
        {
          std::string key = mod_candidate.toString();
          key += "\t" + std::to_string(sms.precursor_charge_) + "\t" + std::to_string(sms.isotope_error_);
          if (! pool_stats[scan_index].seen_candidates.insert(std::move(key)).second) { continue; }
        }

        // Clear peaks + data arrays (ion names / charges) before refilling for the
        // next candidate; getSpectrum appends to whatever is there.
        theo_spectrum.clear(true);
        spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, scoringMaxCharge_(sms.precursor_charge_));
        // Note: TSG emits sorted output when add_metainfo=true (see the
        // sortByPositionPresorted() call at the tail of getSpectrum_); the extra
        // sortByPosition() pass here was a redundant O(N) scan per candidate.

        HyperScore::PSMDetail detail;
        const double score = mass_accuracy_score_
          ? HyperScore::computeMassAccuracy(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum,
                                            mass_error_sd_ppm_, detail)
          : HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);

        // Summarise the candidate before it can be dropped below or pruned at the
        // end of the loop: the pool-derived PSM features describe the whole search
        // space, so a candidate that scored 0 (no fragment matched) still counts
        // and still belongs in the null distribution the z-score is measured against.
        // Each scan_index is owned by exactly one thread, so this needs no guard.
        pool_stats[scan_index].add(score);

        if (score == 0) continue;

        AnnotatedHit_ ah;
        ah.sequence = std::move(mod_candidate);
        ah.score = score;
        // Account for the additional charge hypotheses in the ion-count fractions.
        double seq_length = static_cast<double>(ah.sequence.size()) * scoringMaxCharge_(sms.precursor_charge_);
        ah.prefix_fraction = static_cast<float>(detail.matched_prefix_ions / seq_length);
        ah.suffix_fraction = static_cast<float>(detail.matched_suffix_ions / seq_length);
        ah.mean_error = static_cast<float>(detail.mean_error);
        ah.matched_prefix_ions = static_cast<uint16_t>(detail.matched_prefix_ions);
        ah.matched_suffix_ions = static_cast<uint16_t>(detail.matched_suffix_ions);
        ah.isotope_error = sms.isotope_error_;
        ah.applied_charge = sms.precursor_charge_;
        ah.index_peptide = static_cast<UInt32>(sms.peptide_idx_);
        ah.delta_mass = 0.0;
        if (open_search_mode)
        {
          double theo_mh_plus = ah.sequence.getMZ(1);
          double exp_mz = exp_spectrum.getPrecursors()[0].getMZ();
          double exp_mh_plus = exp_mz * sms.precursor_charge_ - ((sms.precursor_charge_ - 1) * proton_mass_u);
          ah.delta_mass = exp_mh_plus - theo_mh_plus;
        }

        annotated_hits[scan_index].push_back(std::move(ah));
      }

      // Prune to top-N per spectrum to bound memory: up to
      // scoring:max_candidates_per_spectrum hits (each owning a heap AASequence)
      // would otherwise stay resident until postProcessHits_ truncates them.
      // Correct because all of this spectrum's candidates are appended above and
      // scores are independent of the chunk a candidate came from: a hit outside
      // the top-N here can never re-enter it.
      auto& hits = annotated_hits[scan_index];
      if (hits.size() > keep)
      {
        std::partial_sort(hits.begin(), hits.begin() + keep, hits.end(), AnnotatedHit_::hasBetterScore);
        hits.resize(keep);
        hits.shrink_to_fit();
      }
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_scoring_error)
        {
          if (!scoring_error) scoring_error = std::current_exception();
        }
        scoring_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (scoring_error) std::rethrow_exception(scoring_error);

    endProgress();
  }

  // =====================================================================
  // In-memory search: thin wrapper that builds a fresh SearchContext per call.
  // For repeated searches against the same database, prefer the
  // context-taking overload below.
  //
  // When database:chunk_size > 0 and the DB exceeds that size, the database
  // is split into chunks: each chunk builds its own FragmentIndex, all spectra
  // are scored against each chunk, and hits are accumulated across chunks
  // before a single postProcessHits_ + PeptideIndexing + FDR pass. This
  // trades speed (N × spectrum-scoring passes) for memory (only one chunk's
  // FragmentIndex in memory at a time).
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      PeakMap& spectra,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // Reset per-run stats here so the chunked single-file path (searchChunked_,
    // which does not route through search(spectra, ctx, ...)) starts clean too.
    last_run_stats_ = RunStatistics{};

    // Chunking disabled → take the existing single-context path (decoys built
    // lazily by prepareContext). The ctx is locally owned and not reused,
    // so opt in to eager FI release (M1) before PeptideIndexing.
    // ions:by_activation: index c and z+1 ions up front if electron-activated spectra are searched
    const bool electron_ions = countElectronActivated_(spectra) > 0;
    if (database_chunk_size_ == 0)
    {
      // single use: index only the peptides these spectra can reach
      std::function<const PeakMap*(Size)> searched_spectra;
      if (restrictIndexToSpectra_()) { searched_spectra = [&spectra](Size) { return &spectra; }; }
      SearchContext ctx = prepareContext_(std::vector<FASTAFile::FASTAEntry>(fasta_db), electron_ions, searched_spectra);
      ctx.release_fragment_index_after_scoring = true;
      return search(spectra, ctx, protein_ids, peptide_ids);
    }

    // Chunking configured: build the decoy-augmented DB once up-front, then
    // decide based on the AUGMENTED size. If the target DB is 3000 proteins
    // and chunk_size is 5000 with decoys enabled, the augmented DB (6000)
    // still exceeds chunk_size — decide-before-augment would have skipped
    // chunking and built a full FI 2× the user's declared memory budget.
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    auto full_db = buildDecoyAugmentedDB_(fasta_db, strategy);
    if (full_db.size() <= database_chunk_size_)
    {
      // Augmented DB fits in one chunk — use the single-context path but skip
      // prepareContext's internal decoy re-generation by building ctx inline.
      SearchContext ctx;
      ctx.db = std::move(full_db);
      ctx.decoy_string = strategy.decoy_string;
      ctx.decoy_is_prefix = strategy.is_prefix;
      ctx.have_decoys = strategy.have_decoys;
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      startProgress(0, 1, "Building fragment index...");
      ctx.fragment_index.setParameters(fragmentIndexParameters_(electron_ions));
      ctx.fragment_index.build(ctx.db);
      ctx.electron_ions = electron_ions;
      endProgress();
      return search(spectra, ctx, protein_ids, peptide_ids);
    }
    return searchChunked_(spectra, full_db, strategy, protein_ids, peptide_ids);
  }

  // =====================================================================
  // Chunked search implementation. Takes a pre-built decoy-augmented DB
  // (from buildDecoyAugmentedDB_) and splits it into chunks for scoring.
  // Called from the single-file search() wrapper and the multi-file wrapper.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::searchChunked_(
      PeakMap& spectra,
      std::vector<FASTAFile::FASTAEntry>& full_db,
      const DecoyStrategy_& strategy,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    const Size n_chunks = (full_db.size() + database_chunk_size_ - 1) / database_chunk_size_;
    OPENMS_LOG_INFO << "[ProSE] Database chunking enabled: " << full_db.size()
                    << " proteins (incl. decoys), chunk_size=" << database_chunk_size_
                    << " → " << n_chunks << " chunks." << std::endl;

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    bool open_search_mode = isOpenSearchMode_();
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_, peaks_keep_n_,
                       peaks_window_top_, peaks_window_type_, deisotoping_, self_trained_ion_priors_ ? &ion_evidence : nullptr,
                       ion_prior_scored_peaks_, evidence_ptr, query_ptr);

    // ions:by_activation: the chunk indices hold c and z+1 ions if electron-activated spectra are searched
    const Size n_electron_activated = countElectronActivated_(spectra);
    const bool electron_ions = n_electron_activated > 0;
    if (electron_ions)
    {
      OPENMS_LOG_INFO << "[ProSE] " << n_electron_activated << " of " << spectra.size()
                      << " spectra are electron-activated (ETD, ECD, EThcD or ETciD) and are also scored"
                      << " with c and z+1 ions." << std::endl;
    }

    // Effective tolerances — may be narrowed by calibration below. Kept asymmetric
    // (lower / upper separately) so the FragmentIndex query can exploit the full
    // calibrated window rather than the looser max() collapse.
    double effective_precursor_tol_lower = precursor_mass_tolerance_lower_;
    double effective_precursor_tol_upper = precursor_mass_tolerance_upper_;
    double effective_fragment_tol = fragment_mass_tolerance_;

    // Save originals so we can restore algo-level members after search (they may be
    // mutated below for the OSMA mod-match tolerance computation).
    const double orig_prec_tol_lower = precursor_mass_tolerance_lower_;
    const double orig_prec_tol_upper = precursor_mass_tolerance_upper_;
    bool calibration_applied = false;

    // Optional calibration on a strided sample of the full DB. The sample size is
    // bounded (see buildCalibrationSample_) so calibration memory stays O(chunk)
    // regardless of database_chunk_size_. Strided rather than first-N so small
    // chunk_size values don't starve the calibration pool (#9182).
    if (calibration_enabled_ && !open_search_mode)
    {
      std::vector<FASTAFile::FASTAEntry> cal_db = buildCalibrationSample_(full_db);
      FragmentIndex cal_fi;
      cal_fi.setParameters(fragmentIndexParameters_(electron_ions));
      cal_fi.build(cal_db);

      CalibrationResult_ cal = runCalibrationPass_(spectra, cal_fi, cal_db, query_ptr);
      if (cal.success)
      {
        effective_fragment_tol = cal.fragment_tolerance;
        if (!cal.extreme_bias)
        {
          effective_precursor_tol_lower = cal.cal_lower;
          effective_precursor_tol_upper = cal.cal_upper;
          precursor_mass_tolerance_lower_ = cal.cal_lower;
          precursor_mass_tolerance_upper_ = cal.cal_upper;
          calibration_applied = true;
          OPENMS_LOG_INFO << "[ProSE] Calibration (chunked, strided sample): shift=" << cal.precursor_shift
                          << " " << precursor_mass_tolerance_unit_
                          << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                          << " fragment=" << cal.fragment_tolerance << std::endl;
        }
        else
        {
          OPENMS_LOG_WARN << "[ProSE] Calibration: extreme bias, precursor calibration discarded. "
                          << "Fragment calibration applied." << std::endl;
        }
      }
    }

    // 3. Prepare spectrum generators (once).
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // 4. Allocate per-spectrum hit accumulator (persists across chunks).
    vector<vector<AnnotatedHit_>> annotated_hits(spectra.size());
    for (auto& a : annotated_hits) { a.reserve(report_top_hits_); }
    // Accumulates across chunks: every chunk contributes its own candidates for the
    // same spectrum, and the pool features must see all of them.
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    // 5. Chunk loop.
    const Size chunk_size = database_chunk_size_;
    Size chunk_idx = 0;
    for (Size start = 0; start < full_db.size(); start += chunk_size)
    {
      ++chunk_idx;
      const Size end = std::min(start + chunk_size, full_db.size());

      OPENMS_LOG_INFO << "[ProSE] Chunk " << chunk_idx << ": proteins "
                      << start << "–" << (end - 1) << " (" << (end - start) << " proteins)" << std::endl;

      // Build FragmentIndex for this chunk only.
      std::vector<FASTAFile::FASTAEntry> chunk_db(full_db.begin() + start, full_db.begin() + end);
      FragmentIndex chunk_fi;
      {
        Param fi_params = fragmentIndexParameters_(electron_ions);
        // Apply calibrated tolerances (if calibration succeeded above). Asymmetric
        // lower/upper preserved — collapsing to max() would re-open the tight side
        // of the calibrated window and admit spurious decoy candidates.
        fi_params.setValue("fragment:mass_tolerance", effective_fragment_tol);
        fi_params.setValue("precursor:mass_tolerance_lower", effective_precursor_tol_lower);
        fi_params.setValue("precursor:mass_tolerance_upper", effective_precursor_tol_upper);
        chunk_fi.setParameters(fi_params);
      }
      chunk_fi.build(chunk_db);

      // Score all spectra against this chunk's index.
      scoreSpectraAgainstIndex_(spectra, chunk_fi, chunk_db, generators,
                                effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                                open_search_mode, annotated_hits, pool_stats,
                                "Scoring chunk " + StringUtils::toStr(chunk_idx) + "...", query_ptr);

      // Prune to top-N per spectrum after each chunk to bound memory growth.
      // Without this, K chunks × T candidates per spectrum could accumulate
      // K*T hits before postProcessHits_ truncates — problematic for 100M+
      // peptide databases with many chunks. The pruning is correct because
      // scores are independent across chunks: a hit that fails to place in the
      // current top-N cannot improve when more chunks are added.
      const Size keep = std::max(report_top_hits_, Size(2)); // keep ≥2 for delta score
#pragma omp parallel for default(none) shared(annotated_hits, keep)
      for (SignedSize si = 0; si < (SignedSize)annotated_hits.size(); ++si)
      {
        if (annotated_hits[si].size() > keep)
        {
          std::partial_sort(annotated_hits[si].begin(),
                            annotated_hits[si].begin() + keep,
                            annotated_hits[si].end(),
                            AnnotatedHit_::hasBetterScore);
          annotated_hits[si].resize(keep);
        }
      }
    } // end chunk loop
    query_spectra.clear(true);
    for (auto& stats : pool_stats)
    {
      std::unordered_set<std::string>().swap(stats.seen_candidates);
    }

    // 6. Post-process merged hits (sort, annotate, PeptideIndexing, FDR).
    //    This runs once on ALL hits accumulated across all chunks.
    //    Set mod-match tolerance for open-search modification analysis downstream.
    last_mod_match_tolerance_used_ = computeModMatchTolerance_();

    startProgress(0, 1, "Post-processing PSMs...");
    postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      report_top_hits_,
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      std::max(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_),
      effective_fragment_tol,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      "", // no database filename for in-memory search
      evidence_ptr
      );
    evidence_spectra.clear(true); // Release the additional peak list before PeptideIndexing.
    endProgress();

    // 7. PeptideIndexing against the FULL database (not per-chunk).
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    // Use the effective decoy marker/position that was actually searched. A
    // non-empty string (decoy_prefix_ when target-only) avoids PeptideIndexing's
    // own auto-detection; missing_decoy_action=silent keeps it quiet when none match.
    param_pi.setValue("decoy_string", strategy.decoy_string.empty() ? decoy_prefix_ : strategy.decoy_string);
    param_pi.setValue("decoy_string_position", strategy.is_prefix ? "prefix" : "suffix");
    param_pi.setValue("enzyme:name", enzyme_);
    param_pi.setValue("enzyme:specificity",
                      EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
    param_pi.setValue("missing_decoy_action", "silent");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(full_db, protein_ids, peptide_ids);

    // Restore algo-level tolerance members on every return path. Calibration may have
    // mutated them above; leaving them mutated would poison subsequent searches that
    // reuse this ProSEAlgorithm instance (e.g. the multi-file wrapper).
    auto restore_tolerances = [&]()
    {
      if (calibration_applied)
      {
        precursor_mass_tolerance_lower_ = orig_prec_tol_lower;
        precursor_mass_tolerance_upper_ = orig_prec_tol_upper;
      }
    };

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      restore_tolerances();
      if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
        return ExitCodes::INPUT_FILE_EMPTY;
      else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
        return ExitCodes::UNEXPECTED_RESULT;
      else
        return ExitCodes::UNKNOWN_ERROR;
    }

    // 8. PSM-level FDR only. Protein-level FDR is the caller's responsibility
    //    (matching the non-chunked search(spectra, ctx, ...) semantics — that
    //    method also does PSM FDR only; the multi-file wrapper and file-based
    //    searchWithModificationAnalysis apply protein FDR post-call).
    const bool has_decoys = strategy.have_decoys;

    // Optional per-run ion priors need the target/decoy labels and the native scores: after
    // PeptideIndexing, before FDR overwrites the scores.
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Pre-FDR stats (target/decoy counts + HyperScore distribution).
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    if (fdr_psm_ > 0.0 && has_decoys)
    {
      // PSM-level FDR: annotate q-values and filter target+decoy PSMs alike by the threshold.
      // No decoy-specific stripping happens here (decoupled from decoy removal); the decoys that
      // pass are kept because downstream protein-level FDR and cross-file merging need them.
      // Categorical decoy removal happens only at protein-FDR finalization (file-based
      // single-file search below, or ProSE.cpp).
      StopWatch sw_fdr; sw_fdr.start();
      annotatePsmQValues(peptide_ids);
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);
      last_run_stats_.fdr_applied = true;
      last_run_stats_.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      last_run_stats_.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0 && !has_decoys)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << std::endl;
    }

    restore_tolerances();

    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // Protein mapping from the fragment index (peptide:protein_mapping = index):
  // what PeptideIndexing::run() writes, without searching the database
  // =====================================================================
  namespace
  {
    bool isAmbiguousResidue(const char c)
    {
      return c == 'B' || c == 'J' || c == 'Z' || c == 'X';
    }

    // The residues AhoCorasickAmbiguous matches an ambiguous protein residue with: B = D or N, J = I or L,
    // Z = E or Q, X = any of its unambiguous amino acids (every other letter, O and U included)
    const std::string& residuesMatchedBy(const char ambiguous)
    {
      static const std::string b = "DN", j = "IL", z = "EQ";
      static const std::string x = []
      {
        std::string all;
        for (char c = 'A'; c <= 'Z'; ++c)
        {
          if (!AA(c).isAmbiguous()) all += c;
        }
        return all;
      }();
      switch (ambiguous)
      {
        case 'B': return b;
        case 'J': return j;
        case 'Z': return z;
        default: return x;
      }
    }

    // The first @p length (<= 8) residues of @p residues as a number that sorts like them
    uint64_t residueKey(const std::string_view residues, const Size length)
    {
      uint64_t key = 0;
      for (Size i = 0; i < length; ++i) key = (key << 8) | static_cast<unsigned char>(residues[i]);
      return key;
    }

    using Occurrence = std::pair<UInt32, UInt32>; // {protein index, position}

    /**
      @brief The spans that step 3 of buildProteinMapping_() looks up for one ambiguous residue of a protein: those
      with no ambiguous residue before it and 1 to aaa_max ambiguous residues from it on (so every span over ambiguous
      residues belongs to its first one).

      With r_0 = position < r_1 < ... the ambiguous residues from the position on, a span [start, end) with start <=
      r_0 contains exactly k of them if r_{k-1} < end <= r_k (r_k = the protein end if there is none). Each span is
      looked up once per combination of the residues its ambiguous ones stand for. Used for both the count (before
      any lookup, to decide on the fallback) and the lookups, so the two agree.
    */
    class AmbiguousSpans
    {
    public:
      /// PeptideIndexing's maximum of aaa_max
      static constexpr Size MAX_AMBIGUOUS = 10;

      AmbiguousSpans(const std::string& protein, const Size position, const Size max_length, const Size aaa_max) :
        protein_(protein), first_start_(position), aaa_max_(std::min(aaa_max, MAX_AMBIGUOUS))
      {
        // the first start: after the previous ambiguous residue, and close enough to reach the position
        while (first_start_ > 0 && position - first_start_ + 1 < max_length && !isAmbiguousResidue(protein[first_start_ - 1]))
        {
          --first_start_;
        }
        // the ambiguous residues from the position on that a span can reach: up to aaa_max + 1 (the last one only
        // bounds the spans with aaa_max)
        const Size reach = std::min(protein.size(), position + max_length);
        ambiguous_[0] = position;
        for (Size i = position + 1; i < reach && count_ <= aaa_max_; ++i)
        {
          if (isAmbiguousResidue(protein[i])) ambiguous_[count_++] = i;
        }
      }

      /// The largest number of ambiguous residues in a span
      Size maxAmbiguous() const { return std::min(aaa_max_, count_); }

      /// The k-th ambiguous residue from the position on (k = 0: the position)
      Size ambiguous(const Size k) const { return ambiguous_[k]; }

      /// The starts [first, last] of the spans of @p length with exactly @p k (1..maxAmbiguous()) ambiguous residues
      /// (first > last if there are none)
      std::pair<SignedSize, SignedSize> starts(const Size length, const Size k) const
      {
        const SignedSize end_after = static_cast<SignedSize>(ambiguous_[k - 1]);
        const SignedSize end_at_most = static_cast<SignedSize>(k < count_ ? ambiguous_[k] : protein_.size());
        const SignedSize len = static_cast<SignedSize>(length);
        return {std::max(static_cast<SignedSize>(first_start_), end_after + 1 - len),
                std::min(static_cast<SignedSize>(ambiguous_[0]), end_at_most - len)};
      }

      /// The sequences looked up for the spans of @p length (combinations saturate at 2^32, far above any budget)
      Size lookups(const Size length) const
      {
        Size total = 0, combinations = 1;
        for (Size k = 1; k <= maxAmbiguous(); ++k)
        {
          combinations = std::min(combinations * residuesMatchedBy(protein_[ambiguous_[k - 1]]).size(), Size(1) << 32);
          const auto [first, last] = starts(length, k);
          if (first <= last) total += combinations * static_cast<Size>(last - first + 1);
        }
        return total;
      }

    private:
      const std::string& protein_;
      Size first_start_;
      Size aaa_max_;
      std::array<Size, MAX_AMBIGUOUS + 1> ambiguous_{};
      Size count_ = 1; ///< entries of ambiguous_
    };
  }

  ProSEAlgorithm::ProteinMapping_ ProSEAlgorithm::buildProteinMapping_(const FragmentIndex& index,
                                                                      const std::vector<FASTAFile::FASTAEntry>& db,
                                                                      std::vector<UInt32> candidates,
                                                                      const Param& indexer_parameters)
  {
    ProteinMapping_ mapping;

    // PeptideIndexing settings the mapping reproduces (see ProteinMapping_)
    if (static_cast<Int>(indexer_parameters.getValue("mismatches_max")) != 0
        || indexer_parameters.getValue("IL_equivalent").toBool()
        || !indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool()
        || indexer_parameters.getValue("write_protein_sequence").toBool()
        || indexer_parameters.getValue("write_protein_description").toBool()
        || indexer_parameters.getValue("keep_unreferenced_proteins").toBool()
        || indexer_parameters.getValue("decoy_string").toString().empty())
    {
      mapping.fallback_reason = "PeptideIndexing settings other than no mismatches, no I/L equivalence, N-terminal methionine "
                                "cleavage, no protein sequences or descriptions, no unreferenced proteins and a given decoy string";
      return mapping;
    }
    const std::string enzyme_name = indexer_parameters.getValue("enzyme:name").toString();
    const std::string specificity = indexer_parameters.getValue("enzyme:specificity").toString();
    const Param index_parameters = index.getParameters();
    if (specificity != EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL]
        || enzyme_name == EnzymaticDigestion::UnspecificCleavage || enzyme_name == EnzymaticDigestion::NoCleavage
        || index_parameters.getValue("enzyme").toString() != enzyme_name
        || index_parameters.getValue("peptide:enzyme_specificity").toString() != specificity)
    {
      mapping.fallback_reason = "not a fully specific digestion with the enzyme of the fragment index";
      return mapping;
    }
    if (!index.hasProteinOccurrences(db))
    {
      mapping.fallback_reason = "the fragment index does not list every protein occurrence of its peptides "
                                "(SNES, protein-terminal modifications, or not built from this database)";
      return mapping;
    }
    const Size aaa_max = static_cast<Size>(static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(enzyme_name);
    enzyme.setSpecificity(EnzymaticDigestion::SPEC_FULL);

    // The candidates, once each, and the lengths of their sequences. A hit with an ambiguous residue leaves the
    // mapping to PeptideIndexing; checked first, since it needs no pass over the database.
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) return mapping; // nothing to map (applyProteinMapping_() decides)
    const std::vector<FragmentIndex::Peptide>& peptides = index.getPeptides();
    if (candidates.back() >= peptides.size())
    {
      mapping.fallback_reason = "a candidate the fragment index does not hold";
      return mapping;
    }
    const auto residues_of = [&](const UInt32 candidate)
    {
      const FragmentIndex::Peptide& peptide = peptides[candidate];
      return std::string_view(db[peptide.protein_idx].sequence).substr(peptide.sequence_.first, peptide.sequence_.second);
    };
    bool unambiguous_candidates = true;
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : unambiguous_candidates)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      const std::string_view residues = residues_of(candidates[i]);
      unambiguous_candidates = unambiguous_candidates && std::none_of(residues.begin(), residues.end(), isAmbiguousResidue);
    }
    if (!unambiguous_candidates)
    {
      mapping.fallback_reason = "a hit with an ambiguous residue";
      return mapping;
    }
    Size min_length = std::numeric_limits<Size>::max(), max_length = 0;
    for (const UInt32 candidate : candidates)
    {
      min_length = std::min<Size>(min_length, peptides[candidate].sequence_.second);
      max_length = std::max<Size>(max_length, peptides[candidate].sequence_.second);
    }
    std::vector<uint8_t> length_used(max_length + 1, 0);
    for (const UInt32 candidate : candidates) length_used[peptides[candidate].sequence_.second] = 1;

    // The database: letters A-Z only (PeptideIndexing removes '*' and skips other symbols, which shifts positions),
    // no stretch of more than aaa_max X (PeptideIndexing splits proteins at such stretches), and at most a budget of
    // sequences that 3. looks up for the spans over ambiguous residues (beyond it, PeptideIndexing does this better).
    // The lookups are counted here, without looking anything up, so that a fallback costs little more than this pass.
    constexpr Size lookup_budget = Size(1) << 22;
    std::array<uint8_t, 256> residue_class{}; // 0: unambiguous letter, 1: ambiguous letter, 2: anything else
    residue_class.fill(2);
    for (char c = 'A'; c <= 'Z'; ++c) residue_class[static_cast<unsigned char>(c)] = isAmbiguousResidue(c) ? 1 : 0;
    std::vector<uint8_t> with_spans(db.size(), 0); // has spans for 3.
    bool letters_only = true, short_stretches = true;
    std::atomic<Size> lookups{0}; // exact while at most lookup_budget (then counting stops)
#pragma omp parallel for schedule(dynamic, 256) reduction(&& : letters_only, short_stretches)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      const std::string& protein = db[p].sequence;
      uint8_t classes = 0; // OR of its residue classes
      for (const char c : protein) classes |= residue_class[static_cast<unsigned char>(c)];
      letters_only = letters_only && classes < 2;
      if (classes != 1) continue;
      Size stretch = 0;
      for (const char c : protein)
      {
        stretch = (c == 'X') ? stretch + 1 : 0;
        short_stretches = short_stretches && stretch <= aaa_max;
      }
      Size protein_lookups = 0;
      for (Size position = 0; position < protein.size(); ++position)
      {
        if (!isAmbiguousResidue(protein[position])) continue;
        if (protein_lookups > lookup_budget || lookups.load(std::memory_order_relaxed) > lookup_budget) break;
        const AmbiguousSpans spans(protein, position, max_length, aaa_max);
        for (Size length = min_length; length <= max_length; ++length)
        {
          if (length_used[length]) protein_lookups += spans.lookups(length);
        }
      }
      with_spans[p] = protein_lookups > 0;
      lookups.fetch_add(std::min(protein_lookups, lookup_budget + 1), std::memory_order_relaxed);
    }
    if (!letters_only || !short_stretches)
    {
      mapping.fallback_reason = letters_only ? "the database has stretches of more than aaa_max X"
                                             : "the database has symbols other than the letters A-Z";
      return mapping;
    }
    if (lookups.load() > lookup_budget)
    {
      mapping.fallback_reason = "too many spans over ambiguous residues in the database";
      return mapping;
    }

    // 1. the spans of the digest: per candidate, then per sequence
    std::vector<std::vector<Occurrence>> candidate_occurrences(candidates.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      index.getProteinOccurrences(candidates[i], db, candidate_occurrences[i]);
    }
    // One entry per sequence (the candidates of a sequence differ in their modifications only: they have the same
    // spans). The sequences go to ProteinMapping_::SHARDS hash tables by their hash, filled in parallel; the entries
    // are numbered table by table, in the order of the candidates.
    constexpr Size shards = ProteinMapping_::SHARDS;
    std::vector<Size> shard_of(candidates.size());
#pragma omp parallel for schedule(static)
    for (SignedSize i = 0; i < static_cast<SignedSize>(candidates.size()); ++i)
    {
      shard_of[i] = std::hash<std::string_view>{}(residues_of(candidates[i])) % shards;
    }
    std::vector<std::vector<Size>> shard_candidates(shards);
    for (Size i = 0; i < candidates.size(); ++i) shard_candidates[shard_of[i]].push_back(i);
    mapping.entries.resize(shards);
    std::vector<std::vector<Size>> shard_entries(shards); // per table: the candidate of each of its entries
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      auto& table = mapping.entries[shard];
      table.reserve(shard_candidates[shard].size());
      for (const Size i : shard_candidates[shard])
      {
        const auto [entry, added] = table.try_emplace(residues_of(candidates[i]), shard_entries[shard].size());
        if (added)
        {
          shard_entries[shard].push_back(i);
        }
        else
        { // another modified form of the same residues: the same spans (merged and made unique below)
          std::vector<Occurrence>& first = candidate_occurrences[shard_entries[shard][entry->second]];
          first.insert(first.end(), candidate_occurrences[i].begin(), candidate_occurrences[i].end());
        }
      }
    }
    std::vector<Size> shard_start(shards + 1, 0);
    for (Size shard = 0; shard < shards; ++shard) shard_start[shard + 1] = shard_start[shard] + shard_entries[shard].size();
    mapping.sequences.resize(shard_start[shards]);
    mapping.occurrences.resize(shard_start[shards]);
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize shard = 0; shard < static_cast<SignedSize>(shards); ++shard)
    {
      for (auto& [sequence, entry] : mapping.entries[shard]) entry += shard_start[shard];
      for (Size k = 0; k < shard_entries[shard].size(); ++k)
      {
        const Size i = shard_entries[shard][k];
        mapping.sequences[shard_start[shard] + k] = residues_of(candidates[i]);
        mapping.occurrences[shard_start[shard] + k] = std::move(candidate_occurrences[i]);
      }
    }

    std::vector<std::tuple<Size, UInt32, UInt32>> extra; // spans of 2. and 3.: {sequence entry, protein, position}

    // 2. positions 1 and 2 of proteins that start with M: found through their first residues
    {
      // {key of the residues from the position on, protein, position}, sorted. Filled in parallel by chunks of
      // proteins, grouped by the first residue (the highest byte of the key), then each group sorted on its own.
      const Size key_length = std::min<Size>(min_length, 8);
      using Start = std::tuple<uint64_t, UInt32, UInt32>;
      const auto is_start = [&](const std::string& protein, Size position)
      {
        return position + min_length <= protein.size() && protein[0] == 'M';
      };
      constexpr Size letters = 26;
      const Size num_chunks = std::max<Size>(1, std::min<Size>(256, db.size() / 1024));
      std::vector<Size> offset(num_chunks * letters + 1, 0); // [letter][chunk]: count, then where the chunk fills it
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (is_start(protein, position)) ++offset[(protein[position] - 'A') * num_chunks + chunk + 1];
          }
        }
      }
      for (Size i = 1; i < offset.size(); ++i) offset[i] += offset[i - 1];
      std::vector<Start> starts(offset.back());
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize chunk = 0; chunk < static_cast<SignedSize>(num_chunks); ++chunk)
      {
        for (Size p = db.size() * chunk / num_chunks; p < db.size() * (chunk + 1) / num_chunks; ++p)
        {
          const std::string& protein = db[p].sequence;
          for (Size position = 1; position <= 2; ++position)
          {
            if (!is_start(protein, position)) continue;
            starts[offset[(protein[position] - 'A') * num_chunks + chunk]++] =
              Start(residueKey(std::string_view(protein).substr(position), key_length), static_cast<UInt32>(p), static_cast<UInt32>(position));
          }
        }
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize letter = 0; letter < static_cast<SignedSize>(letters); ++letter)
      {
        // after the fill, offset[letter * num_chunks + num_chunks - 1] is where the letter's group ends
        const Size begin = letter == 0 ? 0 : offset[letter * num_chunks - 1];
        const Size end = offset[(letter + 1) * num_chunks - 1];
        std::sort(starts.begin() + begin, starts.begin() + end);
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme; // as PeptideIndexing: one per thread
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
#pragma omp for schedule(dynamic, 256) nowait
        for (SignedSize entry = 0; entry < static_cast<SignedSize>(mapping.sequences.size()); ++entry)
        {
          const std::string_view sequence = mapping.sequences[entry];
          const uint64_t key = residueKey(sequence, key_length);
          for (auto s = std::lower_bound(starts.begin(), starts.end(), std::make_tuple(key, UInt32(0), UInt32(0)));
               s != starts.end() && std::get<0>(*s) == key; ++s)
          {
            const UInt32 protein = std::get<1>(*s), position = std::get<2>(*s);
            const std::string& protein_sequence = db[protein].sequence;
            if (position + sequence.size() <= protein_sequence.size() && protein_sequence.compare(position, sequence.size(), sequence) == 0
                && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(position), static_cast<int>(sequence.size()), true, true, false))
            {
              found.emplace_back(static_cast<Size>(entry), protein, position);
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // 3. spans over ambiguous residues (AmbiguousSpans; within the budget checked above): each looked up with every
    //    combination of the residues its ambiguous ones stand for
    {
      // work items: an ambiguous residue and a length with spans (a few residues with several ambiguous ones around
      // them would otherwise be the critical path); at most one per two lookups, so bounded by the budget
      struct Item { UInt32 protein, position, length; };
      std::vector<Item> items;
#pragma omp parallel
      {
        std::vector<Item> thread_items;
#pragma omp for schedule(dynamic, 16) nowait
        for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
        {
          if (!with_spans[p]) continue;
          const std::string& protein = db[p].sequence;
          for (Size position = 0; position < protein.size(); ++position)
          {
            if (!isAmbiguousResidue(protein[position])) continue;
            const AmbiguousSpans spans(protein, position, max_length, aaa_max);
            for (Size length = min_length; length <= max_length; ++length)
            {
              if (length_used[length] && spans.lookups(length) > 0)
              {
                thread_items.push_back({static_cast<UInt32>(p), static_cast<UInt32>(position), static_cast<UInt32>(length)});
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        items.insert(items.end(), thread_items.begin(), thread_items.end());
      }
#pragma omp parallel
      {
        const ProteaseDigestion thread_enzyme = enzyme;
        std::vector<std::tuple<Size, UInt32, UInt32>> found;
        std::string window;
        std::array<Size, AmbiguousSpans::MAX_AMBIGUOUS> choice{};
#pragma omp for schedule(dynamic, 1) nowait
        for (SignedSize i = 0; i < static_cast<SignedSize>(items.size()); ++i)
        {
          const Item& item = items[i];
          const std::string& protein_sequence = db[item.protein].sequence;
          const AmbiguousSpans spans(protein_sequence, item.position, max_length, aaa_max);
          for (Size k = 1; k <= spans.maxAmbiguous(); ++k)
          {
            const auto [first, last] = spans.starts(item.length, k);
            for (SignedSize start = first; start <= last; ++start)
            {
              window.assign(protein_sequence, static_cast<Size>(start), item.length);
              std::fill_n(choice.begin(), k, 0);
              while (true)
              {
                for (Size w = 0; w < k; ++w)
                {
                  window[spans.ambiguous(w) - static_cast<Size>(start)] = residuesMatchedBy(protein_sequence[spans.ambiguous(w)])[choice[w]];
                }
                const Size entry = mapping.find(window);
                if (entry < mapping.sequences.size()
                    && thread_enzyme.isValidProduct(protein_sequence, static_cast<int>(start), static_cast<int>(item.length), true, true, false))
                {
                  found.emplace_back(entry, item.protein, static_cast<UInt32>(start));
                }
                Size w = 0; // the next combination
                for (; w < k; ++w)
                {
                  if (++choice[w] < residuesMatchedBy(protein_sequence[spans.ambiguous(w)]).size()) break;
                  choice[w] = 0;
                }
                if (w == k) break;
              }
            }
          }
        }
#pragma omp critical (ProSEAlgorithm_proteinMapping)
        extra.insert(extra.end(), found.begin(), found.end());
      }
    }

    // all spans of every sequence, ascending and once (2. and 3. may find spans of 1. again)
    for (const auto& [entry, protein, position] : extra) mapping.occurrences[entry].emplace_back(protein, position);
#pragma omp parallel for schedule(dynamic, 256)
    for (SignedSize k = 0; k < static_cast<SignedSize>(mapping.occurrences.size()); ++k)
    {
      std::vector<Occurrence>& occurrences = mapping.occurrences[k];
      std::sort(occurrences.begin(), occurrences.end());
      occurrences.erase(std::unique(occurrences.begin(), occurrences.end()), occurrences.end());
    }
    return mapping;
  }

  bool ProSEAlgorithm::applyProteinMapping_(const ProteinMapping_& mapping,
                                            const std::vector<FASTAFile::FASTAEntry>& db,
                                            const Param& indexer_parameters,
                                            std::vector<ProteinIdentification>& protein_ids,
                                            PeptideIdentificationList& peptide_ids,
                                            std::string& fallback_reason)
  {
    if (!mapping.fallback_reason.empty())
    {
      fallback_reason = mapping.fallback_reason;
      return false;
    }
    if (protein_ids.size() != 1)
    {
      fallback_reason = "not exactly one identification run";
      return false;
    }
    {
      // PeptideIndexing changes its cleavage rules for results of X! Tandem and MS-GF+
      std::string engine = protein_ids[0].getOriginalSearchEngineName();
      StringUtils::toUpper(engine);
      const ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
      if (engine == "XTANDEM" || engine == "MS-GF+" || engine == "MSGFPLUS"
          || search_parameters.metaValueExists("SE:XTandem") || search_parameters.metaValueExists("SE:MS-GF+"))
      {
        fallback_reason = "results of X! Tandem or MS-GF+";
        return false;
      }
    }

    // the sequence entry of every hit (PeptideIndexing maps the unmodified sequence)
    std::vector<Size> first_hit(peptide_ids.size() + 1, 0);
    for (Size i = 0; i < peptide_ids.size(); ++i) first_hit[i + 1] = first_hit[i] + peptide_ids[i].getHits().size();
    if (first_hit.back() == 0)
    {
      fallback_reason = "no peptide hits"; // PeptideIndexing's own handling (PEPTIDE_IDS_EMPTY)
      return false;
    }
    std::vector<Size> hit_entry(first_hit.back());
    bool all_mapped = true;
#pragma omp parallel for schedule(dynamic, 64) reduction(&& : all_mapped)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      const std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const Size entry = mapping.find(hits[h].getSequence().toUnmodifiedString());
        const bool mapped = entry < mapping.sequences.size() && !mapping.occurrences[entry].empty();
        all_mapped = all_mapped && mapped;
        if (mapped) hit_entry[first_hit[i] + h] = entry;
      }
    }
    if (!all_mapped)
    {
      fallback_reason = "a hit whose sequence was not mapped from the fragment index";
      return false;
    }

    const std::string decoy_string = indexer_parameters.getValue("decoy_string").toString();
    const bool decoy_prefix = indexer_parameters.getValue("decoy_string_position").toString() == "prefix";
    std::vector<uint8_t> is_decoy(db.size());
#pragma omp parallel for schedule(static)
    for (SignedSize p = 0; p < static_cast<SignedSize>(db.size()); ++p)
    {
      is_decoy[p] = decoy_prefix ? StringUtils::hasPrefix(db[p].identifier, decoy_string) : StringUtils::hasSuffix(db[p].identifier, decoy_string);
    }
    // the proteins of the hits
    std::vector<uint8_t> entry_used(mapping.occurrences.size(), 0);
    for (const Size entry : hit_entry) entry_used[entry] = 1;
    std::vector<uint8_t> protein_used(db.size(), 0);
    for (Size entry = 0; entry < entry_used.size(); ++entry)
    {
      if (!entry_used[entry]) continue;
      for (const Occurrence& occurrence : mapping.occurrences[entry]) protein_used[occurrence.first] = 1;
    }
    std::vector<UInt32> proteins; // in database order (PeptideIndexing's std::set of protein indices)
    bool any_decoy = false;
    for (Size p = 0; p < db.size(); ++p)
    {
      if (!protein_used[p]) continue;
      proteins.push_back(static_cast<UInt32>(p));
      any_decoy = any_decoy || is_decoy[p];
    }
    if (!any_decoy && indexer_parameters.getValue("missing_decoy_action").toString() != "silent")
    {
      fallback_reason = "no hit maps to a decoy"; // PeptideIndexing's warning or error
      return false;
    }

    // The meta value names in the order in which PeptideIndexing registers them (its first hit is mapped): MetaInfo
    // keeps, and idXML writes, the values of a hit ordered by registry index
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt target_decoy = registry.registerName("target_decoy");
    const UInt protein_references = registry.registerName("protein_references");
    const DataValue target("target"), decoy("decoy"), target_and_decoy("target+decoy"), unique("unique"), non_unique("non-unique");

    // the evidences of every hit, in PeptideIndexing's order (protein, position)
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
    {
      std::vector<PeptideHit>& hits = peptide_ids[i].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        PeptideHit& hit = hits[h];
        const std::vector<Occurrence>& occurrences = mapping.occurrences[hit_entry[first_hit[i] + h]];
        const Size length = hit.getSequence().size();
        std::vector<PeptideEvidence> evidences;
        evidences.reserve(occurrences.size());
        bool in_target = false, in_decoy = false;
        Size protein_count = 0;
        UInt32 last_protein = std::numeric_limits<UInt32>::max();
        for (const auto& [protein, position] : occurrences)
        {
          const std::string& sequence = db[protein].sequence;
          const char before = (position == 0) ? PeptideEvidence::N_TERMINAL_AA : sequence[position - 1];
          const char after = (position + length >= sequence.size()) ? PeptideEvidence::C_TERMINAL_AA : sequence[position + length];
          evidences.emplace_back(db[protein].identifier, static_cast<Int>(position), static_cast<Int>(position + length) - 1, before, after);
          if (protein != last_protein)
          {
            last_protein = protein;
            ++protein_count;
          }
          (is_decoy[protein] ? in_decoy : in_target) = true;
        }
        hit.setPeptideEvidences(std::move(evidences));
        hit.setMetaValue(target_decoy, (in_target && in_decoy) ? target_and_decoy : (in_target ? target : decoy));
        hit.setMetaValue(protein_references, protein_count == 1 ? unique : non_unique);
      }
    }

    // the protein hits: the proteins of the hits (keep_unreferenced_proteins = false), in database order
    std::vector<ProteinHit> protein_hits(proteins.size());
#pragma omp parallel for schedule(static)
    for (SignedSize k = 0; k < static_cast<SignedSize>(proteins.size()); ++k)
    {
      protein_hits[k].setAccession(db[proteins[k]].identifier);
      protein_hits[k].setMetaValue(target_decoy, is_decoy[proteins[k]] ? decoy : target);
    }
    protein_ids[0].getHits() = std::move(protein_hits);

    // PeptideIndexing's settings, as it records them
    ProteaseDigestion enzyme;
    enzyme.setEnzyme(indexer_parameters.getValue("enzyme:name").toString());
    enzyme.setSpecificity(ProteaseDigestion::getSpecificityByName(indexer_parameters.getValue("enzyme:specificity").toString()));
    ProteinIdentification::SearchParameters search_parameters = protein_ids[0].getSearchParameters();
    search_parameters.setMetaValue("PeptideIndexer:decoy_string", decoy_string);
    search_parameters.setMetaValue("PeptideIndexer:decoy_string_position", decoy_prefix ? "prefix" : "suffix");
    search_parameters.setMetaValue("PeptideIndexer:enzyme", enzyme.getEnzymeName());
    search_parameters.setMetaValue("PeptideIndexer:enzyme_specificity", EnzymaticDigestion::NamesOfSpecificity[enzyme.getSpecificity()]);
    search_parameters.setMetaValue("PeptideIndexer:aaa_max", static_cast<Int>(indexer_parameters.getValue("aaa_max")));
    search_parameters.setMetaValue("PeptideIndexer:mismatches_max", static_cast<Int>(indexer_parameters.getValue("mismatches_max")));
    search_parameters.setMetaValue("PeptideIndexer:IL_equivalent", indexer_parameters.getValue("IL_equivalent").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:allow_nterm_protein_cleavage",
                                   indexer_parameters.getValue("allow_nterm_protein_cleavage").toBool() ? "true" : "false");
    search_parameters.setMetaValue("PeptideIndexer:unmatched_action", indexer_parameters.getValue("unmatched_action").toString());
    search_parameters.setMetaValue("PeptideIndexer:missing_decoy_action", indexer_parameters.getValue("missing_decoy_action").toString());
    protein_ids[0].setSearchParameters(std::move(search_parameters));
    return true;
  }

  // =====================================================================
  // In-memory search using a pre-built SearchContext (no index rebuild).
  // Takes the context by non-const reference because the underlying
  // FragmentIndex::querySpectrum() and PeptideIndexing::run() APIs are both
  // non-const, even though neither conceptually mutates shared state.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      PeakMap& spectra,
      SearchContext& ctx,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // Reset the per-run statistics bridge for this file. Callers copy it into
    // their SearchResult::stats after search() returns OK.
    last_run_stats_ = RunStatistics{};

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    bool open_search = isOpenSearchMode_();
    OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search ? "true" : "false")
                    << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] "
                    << precursor_mass_tolerance_unit_ << ")" << std::endl;

    startProgress(0, 1, "Filtering spectra...");
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_, peaks_keep_n_,
                       peaks_window_top_, peaks_window_type_, deisotoping_, self_trained_ion_priors_ ? &ion_evidence : nullptr,
                       ion_prior_scored_peaks_, evidence_ptr, query_ptr);
    endProgress();

    // ions:by_activation: electron-activated spectra are also scored with c and z+1 ions, so the
    // index must hold them. A context without them is left unchanged, as other searches may share
    // it; this call builds its own index instead.
    const Size n_electron_activated = countElectronActivated_(spectra);
    FragmentIndex electron_index;
    FragmentIndex* index = &ctx.fragment_index;
    if (n_electron_activated > 0)
    {
      OPENMS_LOG_INFO << "[ProSE] " << n_electron_activated << " of " << spectra.size()
                      << " spectra are electron-activated (ETD, ECD, EThcD or ETciD) and are also scored"
                      << " with c and z+1 ions." << std::endl;
      if (!ctx.electron_ions)
      {
        OPENMS_LOG_WARN << "[ProSE] The prepared fragment index holds no c and z+1 ions; building one for this"
                        << " search. prepareContext(fasta_db, true) prepares a context that has them." << std::endl;
        startProgress(0, 1, "Building fragment index with c and z+1 ions...");
        electron_index.setParameters(fragmentIndexParameters_(true));
        electron_index.build(ctx.db);
        endProgress();
        index = &electron_index;
      }
    }

    // Reference the prepared (decoy-augmented) database and the fragment index to search.
    std::vector<FASTAFile::FASTAEntry>& db = ctx.db;
    FragmentIndex& fragment_index_ = *index;

    // Effective tolerances: may be overridden by calibration pass below.
    // The precursor scalar passed to postProcessHits_ is the widest bound
    // (legacy API); calibration may narrow it. Fragment tolerance is a
    // single scalar throughout.
    double effective_precursor_tol = std::max(precursor_mass_tolerance_lower_,
                                              precursor_mass_tolerance_upper_);
    double effective_fragment_tol = fragment_mass_tolerance_;

    // Save original FragmentIndex parameters AND algo-level tolerance members so we
    // can restore them after the search if calibration modifies query-time tolerances.
    // This avoids persistent mutation of the shared SearchContext and the algorithm
    // instance across multi-file calls (the multi-file wrapper reuses a single
    // ProSEAlgorithm instance, so member leaks corrupt later runs).
    const Param fi_params_original = fragment_index_.getParameters();
    const double orig_precursor_mass_tolerance_lower = precursor_mass_tolerance_lower_;
    const double orig_precursor_mass_tolerance_upper = precursor_mass_tolerance_upper_;
    bool fi_params_modified = false;

    // --- Optional calibration pass ---
    if (calibration_enabled_ && !open_search)
    {
      startProgress(0, 1, "Running calibration pass...");
      StopWatch sw_cal; sw_cal.start();
      last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db, query_ptr);
      sw_cal.stop();
      last_run_stats_.seconds_calibration = sw_cal.getClockTime();
      const CalibrationResult_& cal = last_calibration_result_;
      endProgress();

      if (cal.success)
      {
        // Fragment-side calibration is independent of the precursor positive-magnitude
        // representability check — always apply it when calibration succeeds.
        Param fi_params = fi_params_original;
        fi_params.setValue("fragment:mass_tolerance", cal.fragment_tolerance);
        effective_fragment_tol = cal.fragment_tolerance;

        if (!cal.extreme_bias)
        {
          // Precursor calibration is representable in the positive-magnitude schema —
          // apply the calibrated bounds.
          fi_params.setValue("precursor:mass_tolerance_lower", cal.cal_lower);
          fi_params.setValue("precursor:mass_tolerance_upper", cal.cal_upper);

          // Refresh algo-level member copies — OpenSearchModificationAnalysis reads them
          // via computeModMatchTolerance_() and would otherwise see stale pre-calibration
          // values. Restored to originals by restore_fi_params() below on return paths.
          precursor_mass_tolerance_lower_ = cal.cal_lower;
          precursor_mass_tolerance_upper_ = cal.cal_upper;
          effective_precursor_tol = std::max(cal.cal_lower, cal.cal_upper);

          OPENMS_LOG_INFO << "[ProSE] Calibration: shift=" << cal.precursor_shift
                          << " spread=" << cal.precursor_spread << " "
                          << precursor_mass_tolerance_unit_
                          << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                          << std::endl;
        }
        else
        {
          OPENMS_LOG_WARN << "[ProSE] Calibration: |shift|=" << std::abs(cal.precursor_shift)
                          << " > spread=" << cal.precursor_spread << " "
                          << precursor_mass_tolerance_unit_ << " - precursor calibration discarded. "
                          << "The true signed window ["
                          << (cal.precursor_shift - cal.precursor_spread) << ", "
                          << (cal.precursor_shift + cal.precursor_spread)
                          << "] lies entirely on one side of zero (not representable in the "
                          << "positive-magnitude schema without loosening). Fragment calibration "
                          << "still applied. Fix external calibration, or configure "
                          << "mass_tolerance_lower/_upper manually." << std::endl;
          // Precursor bounds unchanged; fragment_tolerance applied above.
        }

        fragment_index_.setParameters(fi_params);
        fi_params_modified = true;
      }
    }

    // Capture the mod-match tolerance reflecting the current (post-calibration)
    // member state. This fires unconditionally — even for closed searches that
    // won't invoke OpenSearchModificationAnalysis — so tests can observe whether
    // calibration wiring reaches the helper regardless of search mode. search()
    // itself does not run OSMA: the searchWithModificationAnalysis wrappers run it
    // on the PSMs returned here and read this field instead of re-computing (by
    // then restore_fi_params() has reset the members to user-configured values).
    last_mod_match_tolerance_used_ = computeModMatchTolerance_();

    // spectrum generators with the ion series of the FragmentIndex
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(report_top_hits_); }
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    bool open_search_mode = open_search;

    StopWatch sw_search; sw_search.start();
    scoreSpectraAgainstIndex_(spectra, fragment_index_, db, generators,
                              effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                              open_search_mode, annotated_hits, pool_stats,
                              "Scoring peptide models against spectra...", query_ptr);
    query_spectra.clear(true);

    // The PeptideIndexing that maps the hits to their proteins below.
    // The PeptideIndexer drops peptides whose termini do not match the configured
    // specificity, so it must agree with the search-time setting — otherwise
    // semi-specific / non-specific PSMs would be silently filtered out here.
    PeptideIndexing indexer;
    {
      Param param_pi = indexer.getParameters();
      // Use the effective decoy marker/position recorded in the context (the same
      // decoys that were searched), avoiding PeptideIndexing's own auto-detection.
      param_pi.setValue("decoy_string", ctx.decoy_string.empty() ? decoy_prefix_ : ctx.decoy_string);
      param_pi.setValue("decoy_string_position", ctx.decoy_is_prefix ? "prefix" : "suffix");
      param_pi.setValue("enzyme:name", enzyme_);
      param_pi.setValue("enzyme:specificity",
                        EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
      param_pi.setValue("missing_decoy_action", "silent");
      indexer.setParameters(param_pi);
    }
    // peptide:protein_mapping = index: the protein occurrences of the hits' candidates are read from the index
    // while it is still there (the index holds them all, PeptideIndexing searches the whole database for them)
    ProteinMapping_ protein_mapping;
    if (protein_mapping_from_index_)
    {
      std::vector<UInt32> candidates;
      for (const auto& hits : annotated_hits)
      {
        for (const auto& ah : hits) candidates.push_back(ah.index_peptide);
      }
      protein_mapping = buildProteinMapping_(fragment_index_, db, std::move(candidates), indexer.getParameters());
    }

    // M1: release the fragment index eagerly when the caller opted in (single-
    // use context). All downstream work (postProcessHits_, open-search mod
    // analysis, PeptideIndexing) is FI-independent, and the subsequent
    // Aho-Corasick pass inside PeptideIndexing::run() is the RSS high-water
    // mark of the whole search on large databases. Freeing here cuts
    // hundreds of MB of steady-state peak on human-proteome runs.
    // Not unconditional: external callers building ctx via prepareContext()
    // and calling search(spectra, ctx, ...) multiple times would break.
    if (ctx.release_fragment_index_after_scoring)
    {
      fragment_index_.clear();
    }

    startProgress(0, 1, "Post-processing PSMs...");
    ProSEAlgorithm::postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      report_top_hits_,
      modifications_fixed_,
      modifications_variable_,
      peptide_missed_cleavages_,
      effective_precursor_tol,
      effective_fragment_tol,
      precursor_mass_tolerance_unit_,
      fragment_mass_tolerance_unit_,
      precursor_min_charge_,
      precursor_max_charge_,
      enzyme_,
      "", // no database filename for in-memory search
      evidence_ptr
      );
    evidence_spectra.clear(true); // Release the additional peak list before PeptideIndexing.
    endProgress();
    sw_search.stop();
    last_run_stats_.seconds_search = sw_search.getClockTime();

    // map the hits to their proteins: from the index, or with PeptideIndexing (the same result)
    PeptideIndexing::ExitCodes indexer_exit = PeptideIndexing::ExitCodes::EXECUTION_OK;
    std::string fallback_reason;
    if (protein_mapping_from_index_
        && applyProteinMapping_(protein_mapping, db, indexer.getParameters(), protein_ids, peptide_ids, fallback_reason))
    {
      OPENMS_LOG_INFO << "[ProSE] Protein mapping from the fragment index: " << protein_mapping.sequences.size()
                      << " sequences, " << protein_ids[0].getHits().size() << " proteins." << std::endl;
    }
    else
    {
      if (protein_mapping_from_index_)
      {
        OPENMS_LOG_INFO << "[ProSE] Protein mapping with PeptideIndexing (" << fallback_reason << ")." << std::endl;
      }
      indexer_exit = indexer.run(db, protein_ids, peptide_ids);
    }
    protein_mapping = ProteinMapping_(); // released before the FDR

    // Helper lambda: restore FragmentIndex parameters AND algo-level tolerance members
    // before returning if calibration modified them, so the shared SearchContext and the
    // algorithm instance are both clean for subsequent per-file searches in the
    // multi-file wrapper.
    auto restore_fi_params = [&]()
    {
      if (fi_params_modified)
      {
        fragment_index_.setParameters(fi_params_original);
        precursor_mass_tolerance_lower_ = orig_precursor_mass_tolerance_lower;
        precursor_mass_tolerance_upper_ = orig_precursor_mass_tolerance_upper;
      }
    };

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      restore_fi_params();
      if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
      {
        return ExitCodes::INPUT_FILE_EMPTY;
      }
      else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
      {
        return ExitCodes::UNEXPECTED_RESULT;
      }
      else
      {
        return ExitCodes::UNKNOWN_ERROR;
      }
    }

    // PSM-level FDR filtering. The context records whether decoys are present
    // (generated internally, or external decoys reused from the input FASTA).
    const bool has_decoys = ctx.have_decoys;

    // Optional per-run ion priors need the target/decoy labels and the native scores: after
    // PeptideIndexing, before FDR overwrites the scores.
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Capture pre-FDR stats now — BEFORE FDR, which drops pure-decoy hits
    // (FalseDiscoveryRate default add_decoy_peptides=false), may strip decoys
    // entirely, and overwrites each hit's HyperScore with its q-value.
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    if (fdr_psm_ > 0.0 && has_decoys)
    {
      // PSM-level FDR: annotate q-values and filter target+decoy PSMs alike by the threshold.
      // No decoy-specific stripping happens here (decoupled from decoy removal); the decoys that
      // pass are kept because downstream protein-level FDR and cross-file merging need them.
      // Categorical decoy removal happens only at protein-FDR finalization (file-based
      // single-file search below, or ProSE.cpp).
      StopWatch sw_fdr; sw_fdr.start();
      annotatePsmQValues(peptide_ids);
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);
      last_run_stats_.fdr_applied = true;
      last_run_stats_.achieved_psm_fdr = maxRetainedScore_(peptide_ids);
      sw_fdr.stop();
      last_run_stats_.seconds_fdr = sw_fdr.getClockTime();
    }
    else if (fdr_psm_ > 0.0 && !has_decoys)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but the search has no decoys (decoys=ignore). "
                         "Use decoys=auto to search/generate decoys and enable FDR. Skipping FDR filtering." << endl;
    }

    restore_fi_params();

    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // Shared protein-FDR finalization for a COMPLETE protein set. Single source
  // of truth for both the file-based search() below and the ProSE TOPP tool's
  // single-file path (the two previously held byte-identical copies).
  // =====================================================================
  // static
  void ProSEAlgorithm::applyCompleteSetProteinFDR(
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      const std::string& decoy_string,
      bool decoy_is_prefix,
      double protein_fdr)
  {
    // Defensive: an empty decoy marker would make every protein match the prefix/suffix guard
    // below (hasPrefix(x, "") is always true) and corrupt picked-protein FDR. Callers gate on
    // have_decoys (so decoy_string is non-empty in practice), but this is a public static helper.
    if (decoy_string.empty())
    {
      OPENMS_LOG_WARN << "[ProSE] applyCompleteSetProteinFDR called with an empty decoy marker; "
                      << "skipping protein FDR (cannot identify decoys)." << std::endl;
      return;
    }

    // Aggregate the best PSM score per peptide per protein, then apply picked-protein
    // FDR (Savitski et al. 2015) over the full target+decoy protein set.
    BasicProteinInferenceAlgorithm bpia;
    bpia.run(peptide_ids, protein_ids);

    // Picked-protein FDR needs identified decoy PROTEINS, not merely a decoy database: if PSM-level
    // filtering removed all decoy evidence, q-values would be estimated from targets only. The
    // callers gate on a DB-level / pre-inference "decoys present" check, which does not catch this
    // post-inference target-only edge -- guard it here.
    const bool has_decoy_proteins = std::any_of(
        protein_ids[0].getHits().begin(), protein_ids[0].getHits().end(),
        [&](const ProteinHit& ph) {
          return decoy_is_prefix ? StringUtils::hasPrefix(ph.getAccession(), decoy_string)
                                 : StringUtils::hasSuffix(ph.getAccession(), decoy_string);
        });
    if (!has_decoy_proteins)
    {
      OPENMS_LOG_WARN << "[ProSE] Protein FDR requested but no decoy proteins remain after inference "
                      << "(marker '" << decoy_string << "'); skipping picked-protein FDR to avoid "
                      << "target-only q-values." << std::endl;
      return;
    }

    FalseDiscoveryRate fdr;
    fdr.applyPickedProteinFDR(protein_ids[0], decoy_string, decoy_is_prefix);
    IDFilter::filterHitsByScore(protein_ids, protein_fdr);

    // Decoy removal is the finalization step. removeDecoyHits strips decoy PSMs; decoy
    // proteins then fall out as unreferenced. Repair indistinguishable-protein and
    // protein-group references and drop peptide evidence pointing at removed proteins,
    // else idXML storage fails on dangling references.
    IDFilter::removeDecoyHits(peptide_ids);
    IDFilter::removeEmptyIdentifications(peptide_ids);
    IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);
    IDFilter::updateProteinGroups(protein_ids[0].getIndistinguishableProteins(), protein_ids[0].getHits());
    IDFilter::updateProteinGroups(protein_ids[0].getProteinGroups(), protein_ids[0].getHits());
    IDFilter::removeDanglingProteinReferences(peptide_ids, protein_ids);

    OPENMS_LOG_INFO << "[ProSE] Protein inference + picked-protein FDR: "
                    << protein_ids[0].getHits().size() << " proteins at "
                    << protein_fdr * 100 << "% FDR." << std::endl;
  }

  void ProSEAlgorithm::annotatePsmQValues(PeptideIdentificationList& peptide_ids) const
  {
    FalseDiscoveryRate fdr;
    Param fdr_params = fdr.getParameters();
    fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
    fdr.setParameters(fdr_params);

    // Group the spectra by the number of fragment charges their best hit was scored with. FalseDiscoveryRate
    // competes the best hit of each spectrum (stable score order, as sort() below).
    std::map<int, std::vector<Size>> groups;
    bool separate = fdr_psm_by_scored_charges_ && scoring_multiple_charges_;
    if (separate)
    {
      std::map<int, std::pair<bool, bool>> has_target_decoy;
      std::vector<Size> without_hits;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty())
        {
          without_hits.push_back(i);
          continue;
        }
        id.sort();
        const PeptideHit& best = id.getHits()[0];
        const int group = scoringMaxCharge_(best.getCharge());
        groups[group].push_back(i);
        auto& [has_target, has_decoy] = has_target_decoy[group];
        (best.isDecoy() ? has_decoy : has_target) = true;
      }
      separate = groups.size() > 1
                 && std::all_of(has_target_decoy.begin(), has_target_decoy.end(),
                                [](const auto& entry) { return entry.second.first && entry.second.second; });
      if (groups.size() > 1 && ! separate)
      {
        OPENMS_LOG_WARN << "[ProSE] FDR:PSM_groups: a group of PSMs scored with the same number of fragment charges has no "
                        << "target or no decoy PSM; computing q-values over all PSMs together." << std::endl;
      }
      if (separate && ! without_hits.empty()) // they carry no score, but get the q-value score type with the others
      {
        auto& first = groups.begin()->second;
        first.insert(first.end(), without_hits.begin(), without_hits.end());
      }
    }
    if (! separate)
    {
      fdr.apply(peptide_ids);
      return;
    }
    for (const auto& [group, indices] : groups)
    {
      PeptideIdentificationList part;
      part.reserve(indices.size());
      for (Size i : indices) { part.push_back(std::move(peptide_ids[i])); }
      fdr.apply(part);
      for (Size k = 0; k < indices.size(); ++k) { peptide_ids[indices[k]] = std::move(part[k]); }
    }
  }

  namespace
  {
    // Starts the OpenMP threads of the calling thread, if they do not run yet, and returns their
    // number: 1 without OpenMP, with a single OpenMP thread, and inside a parallel region.
    Size startOpenMPThreads()
    {
      Size threads = 0;
#pragma omp parallel reduction(+ : threads)
      {
        ++threads;
      }
      return threads;
    }

    // The index build waits for the spectra that are read meanwhile, to index only the peptides they can reach, if
    // generating the fragments it skips takes longer than the rest of the read: with few threads (measured on
    // 8,000-spectrum files: faster for all instrument types up to 4 threads, slower for some from 6 threads on) and a
    // read that is short against the build (the read is serial, the build parallel): at most
    // MAX_SPECTRA_BYTES_PER_PEPTIDE / threads bytes of spectra file per peptide (measured with 9 M peptides: faster
    // with 11, 27 and 41 bytes per peptide at 2 threads and with 11 at 4 threads, even with 27, slower with 41).
    constexpr Size MAX_THREADS_WAITING_FOR_SPECTRA = 4;
    constexpr UInt64 MAX_SPECTRA_BYTES_PER_PEPTIDE = 100;

    // FASTAFile::load() with the OpenMP threads: the file is cut into pieces at starts of entries,
    // the threads read the pieces with FASTAFile::readNext(), and the pieces are joined in file
    // order. Same entries as load(): readNext() reads an entry from its '>' up to a line break
    // followed by '>'. A piece starts at a '>' after a line break whose line holds no '>': that
    // line is not a header (whose line break would not end an entry), so the reader of load()
    // starts an entry there as well, and the reader of the piece before stops when it arrives
    // there. Whenever this does not work out (small file, no such place, a reader that fails or
    // does not arrive at the next piece), load() reads the file and reports errors as before.
    void loadFASTA(const std::string& filename, std::vector<FASTAFile::FASTAEntry>& entries)
    {
      std::vector<std::streamoff> pieces{0}; // where each piece starts; the first one where readStart() starts
#ifdef _OPENMP
      if (!omp_in_parallel() && omp_get_max_threads() > 1)
      {
        std::ifstream in(filename, std::ios::binary);
        in.seekg(0, std::ios::end);
        const std::streamoff size = in.tellg();
        // several pieces per thread: the threads do not read equally fast
        const std::streamoff count = std::min<std::streamoff>(4 * omp_get_max_threads(), size / (std::streamoff(1) << 18));
        std::string window(Size(1) << 16, '\0');
        for (std::streamoff i = 1; i < count && in.good(); ++i)
        {
          const std::streamoff from = size / count * i;
          in.seekg(from);
          in.read(window.data(), static_cast<std::streamsize>(window.size()));
          const std::string_view text(window.data(), static_cast<Size>(in.gcount()));
          in.clear(); // reading up to the end of the file is fine
          // lines of the window that are complete: between two line breaks
          for (Size line = text.find('\n'); line != std::string_view::npos && line + 1 < text.size();)
          {
            const Size next = text.find('\n', line + 1);
            if (next == std::string_view::npos || next + 1 >= text.size()) { break; }
            if (text[next + 1] == '>' && text.substr(line + 1, next - line - 1).find('>') == std::string_view::npos)
            {
              if (from + static_cast<std::streamoff>(next) + 1 > pieces.back()) { pieces.push_back(from + static_cast<std::streamoff>(next) + 1); }
              break;
            }
            line = next;
          }
        }
      }
#endif
      if (pieces.size() < 2)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      std::vector<std::vector<FASTAFile::FASTAEntry>> read(pieces.size());
      bool complete = true;
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize i = 0; i < static_cast<SignedSize>(pieces.size()); ++i)
      {
        bool ok = false;
        try // exceptions must not leave the parallel region
        {
          FASTAFile file;
          FASTAFile::FASTAEntry entry;
          file.readStart(filename);
          ok = i == 0 || file.setPosition(pieces[i]);
          if (ok && i + 1 == static_cast<SignedSize>(pieces.size()))
          { // the last piece: up to the end of the file, as load()
            while (file.readNext(entry)) { read[i].push_back(std::move(entry)); }
          }
          else if (ok)
          {
            const std::streampos end = pieces[i + 1];
            while (ok && file.position() < end)
            {
              ok = file.readNext(entry);
              if (ok) { read[i].push_back(std::move(entry)); }
            }
            ok = ok && file.position() == end;
          }
        }
        catch (...)
        {
          ok = false;
        }
        if (!ok)
        {
#pragma omp critical (ProSEAlgorithm_loadFASTA)
          complete = false;
        }
      }
      if (!complete)
      {
        FASTAFile().load(filename, entries);
        return;
      }

      Size total = 0;
      for (const auto& piece : read) { total += piece.size(); }
      entries.clear();
      entries.reserve(total);
      for (auto& piece : read)
      {
        entries.insert(entries.end(), std::make_move_iterator(piece.begin()), std::make_move_iterator(piece.end()));
      }
    }

    // MzMLFile::load() with @p threads OpenMP threads (at least 2). The spectra are cut into chunks at
    // their <spectrum> start tags; every chunk is parsed together with the file's header and closing
    // tags by an MzMLFile of its own (MzMLFile::loadBuffer), and the spectra are joined in file order.
    // The last chunk runs up to </run>, so it also holds the chromatograms. Each thread reads its chunks
    // from the file itself, so the data held beyond the spectra is at most about 9 MB per thread,
    // whatever the size of the file. The header is parsed first on this thread: it gives the experimental
    // settings, registers its meta value names in file order and initialises the XML parser.
    // Returns false if the file does not have the layout the chunks rely on (compressed; a comment,
    // CDATA section or DOCTYPE before the end of the run; not exactly one spectrum list; a count
    // attribute other than the number of spectra found; few spectra) or a chunk fails to parse. The
    // caller then reads the file with MzMLFile::load(), which also reports its errors.
    bool loadMzMLChunked(const std::string& filename, const PeakFileOptions& options, int threads, PeakMap& exp)
    {
      std::ifstream in(filename, std::ios::binary);
      if (!in) { return false; }
      in.seekg(0, std::ios::end);
      const Size size = static_cast<Size>(in.tellg());
      char magic[2] = {0, 0};
      in.seekg(0);
      in.read(magic, 2);
      if (!in || (magic[0] == 'B' && magic[1] == 'Z') || (magic[0] == '\x1f' && magic[1] == '\x8b') || (magic[0] == 'P' && magic[1] == 'K'))
      {
        return false; // compressed (bzip2, gzip, zip) or too short
      }

      // 1. Where the markup of interest starts, found by the threads in blocks of the file.
      enum class Mark : unsigned char { spectrum, list_start, list_end, run_end, declaration };
      const auto is_tag = [](std::string_view text, std::string_view name)
      {
        return text.size() > name.size() && text.substr(0, name.size()) == name
               && std::string_view(" \t\r\n>/").find(text[name.size()]) != std::string_view::npos;
      };
      constexpr Size block = Size(1) << 20;
      constexpr Size lookahead = 16; // longer than the longest name compared below, "</spectrumList" plus one
      const Size n_blocks = (size + block - 1) / block;
      std::vector<std::vector<std::pair<Size, Mark>>> marks(n_blocks);
      bool ok = true;
#pragma omp parallel num_threads(threads)
      {
        std::ifstream file(filename, std::ios::binary);
        std::string buffer;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(n_blocks); ++sb)
        {
          const Size b = static_cast<Size>(sb);
          const Size from = b * block;
          const Size end = std::min(block, size - from); // marks starting before end belong to this block
          buffer.resize(std::min(block + lookahead, size - from));
          file.seekg(static_cast<std::streamoff>(from));
          file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
          if (!file)
          {
#pragma omp critical (ProSEAlgorithm_loadMzMLChunked)
            ok = false;
            file.clear();
            continue;
          }
          for (Size i = 0; i < end; ++i)
          {
            const char* p = static_cast<const char*>(std::memchr(buffer.data() + i, '<', end - i));
            if (p == nullptr) { break; }
            i = static_cast<Size>(p - buffer.data());
            const std::string_view text(p, std::min(lookahead, buffer.size() - i));
            if (is_tag(text, "<spectrumList")) { marks[b].emplace_back(from + i, Mark::list_start); }
            else if (is_tag(text, "<spectrum")) { marks[b].emplace_back(from + i, Mark::spectrum); }
            else if (is_tag(text, "</spectrumList")) { marks[b].emplace_back(from + i, Mark::list_end); }
            else if (is_tag(text, "</run")) { marks[b].emplace_back(from + i, Mark::run_end); }
            else if (text.substr(0, 2) == "<!") { marks[b].emplace_back(from + i, Mark::declaration); }
          }
        }
      }
      constexpr Size none = std::numeric_limits<Size>::max();
      Size list_start = none, list_end = none, run_end = none;
      std::vector<Size> starts; // of the spectra
      for (const auto& block_marks : marks)
      {
        for (const auto& [pos, mark] : block_marks)
        {
          switch (mark)
          {
            case Mark::list_start: ok = ok && list_start == none; list_start = pos; break;
            case Mark::spectrum: ok = ok && list_start != none && list_end == none; starts.push_back(pos); break;
            case Mark::list_end: ok = ok && list_start != none && list_end == none; list_end = pos; break;
            case Mark::run_end: if (run_end == none) { ok = ok && list_end != none; run_end = pos; } break;
            case Mark::declaration: ok = ok && run_end != none; break;
          }
        }
      }
      marks = {};
      if (!ok || run_end == none || starts.size() < 64) { return false; }

      // 2. The header (everything before the first spectrum) and the count attribute of its
      //    <spectrumList>, which each chunk sets to its own number of spectra.
      std::string header(starts[0], '\0');
      in.seekg(0);
      in.read(header.data(), static_cast<std::streamsize>(header.size()));
      if (!in) { return false; }
      Size count_begin = none, count_end = none;
      for (Size p = list_start + std::string_view("<spectrumList").size();;)
      {
        p = header.find_first_not_of(" \t\r\n", p);
        if (p == std::string::npos) { return false; }
        if (header[p] == '>') { break; }
        const Size name_end = header.find_first_of("= \t\r\n>", p);
        Size q = header.find_first_not_of(" \t\r\n", name_end);
        if (name_end == std::string::npos || q == std::string::npos || header[q] != '=') { return false; }
        q = header.find_first_not_of(" \t\r\n", q + 1);
        if (q == std::string::npos || (header[q] != '"' && header[q] != '\'')) { return false; }
        const Size value_end = header.find(header[q], q + 1);
        if (value_end == std::string::npos) { return false; }
        if (std::string_view(header).substr(p, name_end - p) == "count") { count_begin = q + 1; count_end = value_end; }
        p = value_end + 1;
      }
      if (count_begin == none || header.compare(count_begin, count_end - count_begin, std::to_string(starts.size())) != 0)
      {
        return false;
      }
      const std::string end_of_run = std::string("</run></mzML>") + (header.find("<indexedmzML") != std::string::npos ? "</indexedmzML>" : "");
      const auto chunk_text = [&](Size n_spectra, std::string& text)
      {
        text.assign(header, 0, count_begin);
        text += std::to_string(n_spectra);
        text.append(header, count_end, std::string::npos);
      };

      // 3. Chunks of about equal size: two per thread (the threads do not parse equally fast; every chunk
      //    costs the set-up of a parser and handler), of at most 8 MB, and of at least 16 spectra.
      constexpr Size max_chunk_bytes = Size(8) << 20;
      const Size n = starts.size();
      const Size bytes = list_end - starts[0];
      const Size n_chunks = std::min(n / 16, std::max(2 * static_cast<Size>(threads), (bytes + max_chunk_bytes - 1) / max_chunk_bytes));
      std::vector<Size> first{0}; // the first spectrum of each chunk
      for (Size c = 1; c < n_chunks; ++c)
      {
        const Size f = static_cast<Size>(std::lower_bound(starts.begin(), starts.end(), starts[0] + bytes / n_chunks * c) - starts.begin());
        if (f > first.back() && f < n) { first.push_back(f); }
      }
      first.push_back(n);

      // A chunk decodes its spectra on its thread; in batches of the default size, so that it does not
      // hold the encoded data of all its spectra at once.
      PeakFileOptions chunk_options = options;
      chunk_options.setMaxDataPoolSize(PeakFileOptions().getMaxDataPoolSize());
      try
      {
        std::string text;
        chunk_text(0, text);
        text += "</spectrumList>" + end_of_run;
        MzMLFile mzml;
        mzml.getOptions() = chunk_options;
        mzml.loadBuffer(text, exp);
      }
      catch (...)
      {
        return false;
      }
      // A slot for every spectrum of the file (MzMLHandler reserves as many): a chunk moves the spectra it
      // keeps into the first of its slots right after parsing, so that only the chunks being parsed hold
      // spectra of their own.
      const Size n_chunks_used = first.size() - 1;
      std::vector<MSSpectrum> spectra(n);
      std::vector<Size> kept(n_chunks_used, 0);
      std::vector<MSChromatogram> chromatograms;
#pragma omp parallel num_threads(threads)
      {
        std::ifstream file(filename, std::ios::binary);
        std::string chunk;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sc = 0; sc < static_cast<SignedSize>(n_chunks_used); ++sc)
        {
          const Size c = static_cast<Size>(sc);
          const bool last = c + 1 == n_chunks_used;
          const Size from = starts[first[c]];
          const Size to = last ? run_end : starts[first[c + 1]];
          chunk_text(first[c + 1] - first[c], chunk);
          const Size at = chunk.size();
          chunk.resize(at + (to - from));
          file.seekg(static_cast<std::streamoff>(from));
          file.read(chunk.data() + at, static_cast<std::streamsize>(to - from));
          bool parsed = static_cast<bool>(file);
          if (parsed)
          {
            chunk += last ? end_of_run : "</spectrumList>" + end_of_run;
            try // exceptions must not leave the parallel region
            {
              MzMLFile mzml;
              mzml.getOptions() = chunk_options;
              PeakMap part;
              mzml.loadBuffer(chunk, part);
              parsed = part.size() <= first[c + 1] - first[c];
              if (parsed)
              {
                std::move(part.begin(), part.end(), spectra.begin() + static_cast<SignedSize>(first[c]));
                kept[c] = part.size();
                if (last) { chromatograms = std::move(part.getChromatograms()); }
              }
            }
            catch (...)
            {
              parsed = false;
            }
          }
          if (!parsed)
          {
#pragma omp critical (ProSEAlgorithm_loadMzMLChunked)
            ok = false;
            file.clear();
          }
        }
      }
      if (!ok) { return false; }

      Size kept_total = 0; // the spectra in file order: each chunk's slots, without the unused ones
      for (Size c = 0; c < n_chunks_used; ++c)
      {
        for (Size i = first[c]; i < first[c] + kept[c]; ++i, ++kept_total)
        {
          if (kept_total != i) { spectra[kept_total] = std::move(spectra[i]); }
        }
      }
      spectra.erase(spectra.begin() + static_cast<SignedSize>(kept_total), spectra.end());
      exp.setSpectra(std::move(spectra));
      exp.setChromatograms(std::move(chromatograms));
      exp.setLoadedFileType(filename);
      exp.setLoadedFilePath(filename);
      exp.updateRanges();
      return true;
    }

    // The MS2 spectra of a spectrum file (mzML, Bruker .d or Thermo .raw), sorted by RT, read with
    // @p threads OpenMP threads: an mzML file in chunks parsed in parallel (loadMzMLChunked) with two
    // or more, otherwise with FileHandler.
    PeakMap loadMS2Spectra(const std::string& filename, int threads)
    {
      PeakMap spectra;
      FileHandler f;
      f.getOptions().clearMSLevels();
      f.getOptions().addMSLevel(2);
      // MzMLHandler decodes the spectra in a parallel region every maxDataPoolSize spectra (100 by
      // default): with several threads, fewer and larger regions
      if (threads > 1) { f.getOptions().setMaxDataPoolSize(1000); }
      if (threads > 1 && FileHandler::getTypeByFileName(filename) == FileTypes::MZML
          && loadMzMLChunked(filename, f.getOptions(), threads, spectra))
      {
        ChromatogramTools().convertSpectraToChromatograms<PeakMap>(spectra, true); // as FileHandler::loadExperiment()
#ifdef __GLIBC__
        // Each reader thread freed its chunk buffers and parser temporaries between the spectra it
        // keeps, and glibc holds these pages in the thread's arena, where nothing allocates again: with
        // 16 readers 109-165 MB on 8,000-24,000 spectra (4 readers: 43 MB), +0.14 GiB max RSS of the
        // whole search at 64 threads. Return them; with fewer readers the gain is small and returning
        // pages the search is about to reuse cost +0.4% wall time at 4 threads.
        if (threads >= 8) { malloc_trim(0); }
#endif
      }
      else
      {
        f.loadExperiment(filename, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW});
      }
      spectra.sortSpectra(true);
      return spectra;
    }

    // The OpenMP threads a spectrum file is read with: those of this thread, at most 16 (the chunked
    // read gets no faster with more). A read in the background, next to the index build or the search
    // of another file, gets a quarter of them, at least 2: with all of them, it slowed the index build
    // down when the read was hidden behind it anyway.
    int readerThreads(bool background)
    {
#ifdef _OPENMP
      const int threads = omp_get_max_threads();
      return background ? std::clamp(threads / 4, std::min(threads, 2), 16) : std::min(threads, 16);
#else
      return 1;
#endif
    }

    // loadMS2Spectra() on a helper thread. The helper sets its OpenMP threads (readerThreads(true) of
    // this thread): a new thread starts with the default of the process (OMP_NUM_THREADS or all cores),
    // not with what was set for this one (e.g. by TOPPBase from -threads).
    std::future<PeakMap> loadMS2SpectraAsync(const std::string& filename)
    {
      const int threads = readerThreads(true);
      return std::async(std::launch::async, [&filename, threads]()
      {
#ifdef _OPENMP
        omp_set_num_threads(threads);
#endif
        return loadMS2Spectra(filename, threads);
      });
    }

    // Whether the first MBs of an mzML file name an electron-based activation (the terms that
    // MzMLHandler reads as ECD, ETD, ETciD or EThcD). A cheap prediction of what the spectra will
    // hold, used only to decide when to build the index, not what the index holds.
    bool mzMLHeadNamesElectronActivation(const std::string& filename)
    {
      std::string head(Size(8) << 20, '\0');
      std::ifstream in(filename, std::ios::binary);
      in.read(head.data(), static_cast<std::streamsize>(head.size()));
      head.resize(static_cast<Size>(in.gcount()));
      // MS:1000250, MS:1000598, MS:1003182 and MS:1002631: one pass for what they share
      const std::string_view text(head);
      for (Size pos = text.find("MS:100"); pos != std::string_view::npos; pos = text.find("MS:100", pos + 1))
      {
        const std::string_view number = text.substr(pos + 6, 4);
        if (number == "0250" || number == "0598" || number == "3182" || number == "2631")
        {
          return true;
        }
      }
      return false;
    }
  }

  // =====================================================================
  // File-based search: thin I/O wrapper that delegates to in-memory search
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      const std::string& in_spectra, const std::string& in_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // load MS2 map
    PeakMap spectra;

    vector<FASTAFile::FASTAEntry> fasta_db;
    DecoyStrategy_ strategy; // decoys of the searched database, for protein FDR below
    bool strategy_resolved = false;
    ExitCodes ec;
    const Size threads = startOpenMPThreads();
    if (database_chunk_size_ == 0 && FileHandler::getTypeByFileName(in_spectra) == FileTypes::MZML && threads > 1)
    {
      // Unchunked, multi-threaded search of an mzML file: a helper thread reads the spectra while
      // this thread reads the FASTA file and builds the fragment index; then the search continues
      // as search(spectra, fasta_db, ...) does. The two sides share no data: the spectra are used
      // here only after get(), and the registries of OpenMS are written only by the mzML reader
      // in the meantime (FASTA reading, decoys and the index build register no meta value names),
      // so they end up as after reading the spectra first. Other readers (.d, .raw) have not
      // been checked for this and keep the order below.
      // The OpenMP threads of this thread run before the helper starts (see the condition), as
      // they did when the spectra were read first (decoding them was the first parallel region):
      // started next to a busy helper, they more often end up on another NUMA node than this
      // thread, which slows down the index build and the scoring.
      std::future<PeakMap> spectra_ready = loadMS2SpectraAsync(in_spectra);

      // load FASTA
      loadFASTA(in_db, fasta_db);

      // ions:by_activation: the index needs c and z+1 ions if electron-activated spectra are
      // searched, which is known once the spectra are read. Wait for them if they are read already
      // (or failed to be read: get() rethrows) or if the file names an electron-based activation
      // early on. Otherwise build the index without these ions meanwhile, which is what most data
      // need, and rebuild it should such spectra turn up after all (as the multi-file search does
      // when a later file has them): the index searched is the one prepareContext(fasta_db, true)
      // builds.
      const bool spectra_read = spectra_ready.wait_for(std::chrono::seconds(0)) == std::future_status::ready
                                || (ions_by_activation_ && mzMLHeadNamesElectronActivation(in_spectra));
      if (spectra_read) { spectra = spectra_ready.get(); }
      bool spectra_waited = spectra_read;
      // The single-use index holds only the peptides the spectra can reach if they are read by the time the
      // peptides are generated, or if waiting for them pays off (see MAX_THREADS_WAITING_FOR_SPECTRA).
      std::function<const PeakMap*(Size)> searched_spectra;
      if (restrictIndexToSpectra_())
      {
        const UInt64 spectra_bytes = File::fileSize(in_spectra); // UInt64(-1) if unknown: no waiting
        searched_spectra = [&, spectra_bytes](Size peptides) -> const PeakMap*
        {
          if (!spectra_waited)
          {
            const bool wait = threads <= MAX_THREADS_WAITING_FOR_SPECTRA
                              && spectra_bytes <= MAX_SPECTRA_BYTES_PER_PEPTIDE * peptides / threads;
            if (!wait && spectra_ready.wait_for(std::chrono::seconds(0)) != std::future_status::ready)
            {
              return nullptr;
            }
            spectra = spectra_ready.get();
            spectra_waited = true;
          }
          return &spectra;
        };
      }
      // The context takes the entries over instead of copying them. fasta_db is not read
      // afterwards: protein FDR below takes the decoy marker from the context, which holds what
      // resolveDecoyStrategy_(fasta_db) returned.
      SearchContext ctx = prepareContext_(std::move(fasta_db), spectra_read && countElectronActivated_(spectra) > 0, searched_spectra);
      strategy.have_decoys = ctx.have_decoys;
      strategy.decoy_string = ctx.decoy_string;
      strategy.is_prefix = ctx.decoy_is_prefix;
      strategy_resolved = true;
      if (!spectra_read)
      {
        if (!spectra_waited)
        {
          spectra = spectra_ready.get();
          spectra_waited = true;
        }
        if (countElectronActivated_(spectra) > 0)
        {
          startProgress(0, 1, "Building fragment index with c and z+1 ions...");
          ctx.fragment_index.clear();
          ctx.fragment_index.setParameters(fragmentIndexParameters_(true));
          ctx.fragment_index.build(ctx.db, searched_spectra);
          ctx.electron_ions = true;
          endProgress();
        }
      }
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      ec = search(spectra, ctx, protein_ids, peptide_ids);
    }
    else
    {
      spectra = loadMS2Spectra(in_spectra, readerThreads(false));

      // load FASTA
      loadFASTA(in_db, fasta_db);

      // delegate to in-memory search
      ec = search(spectra, fasta_db, protein_ids, peptide_ids);
    }

    if (ec != ExitCodes::EXECUTION_OK)
    {
      return ec;
    }

    // Protein inference + picked-protein FDR for single-file search.
    // Must run before decoy removal so both target and decoy proteins
    // receive aggregated scores from BPIA. Resolve the decoy strategy from the
    // same input FASTA the search used so the marker/position match.
    // The strategy is only needed for protein FDR: do not scan the accessions again without it.
    if (fdr_protein_ > 0.0)
    {
      if (!strategy_resolved) { strategy = resolveDecoyStrategy_(fasta_db); } // else: known from the context
      if (strategy.have_decoys)
      {
        // Single input file = complete experiment, so picked-protein FDR is valid. Use the resolved
        // decoy marker (prefix or suffix, detected by DecoyHelper in resolveDecoyStrategy_) so the
        // shared finalization recognises the same decoys that were searched.
        applyCompleteSetProteinFDR(protein_ids, peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_);
      }
    }

    // patch file-specific metadata
    protein_ids[0].getSearchParameters().db = in_db;
    protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // In-memory searchWithModificationAnalysis
  // =====================================================================
  ProSEAlgorithm::SearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      PeakMap& spectra,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const std::string& output_base_name) const
  {
    SearchResult result;
    result.is_open_search = isOpenSearchMode_();

    result.exit_code = search(spectra, fasta_db, result.protein_ids, result.peptide_ids);

    // Carry whatever statistics the search managed to populate, even on failure.
    result.stats = last_run_stats_;

    if (result.exit_code != ExitCodes::EXECUTION_OK)
    {
      return result;
    }

    if (result.is_open_search)
    {
      OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for open search results..." << std::endl;

      OpenSearchModificationAnalysis mod_analyzer;

      std::string output_file;
      if (!output_base_name.empty())
      {
        output_file = output_base_name + "_ModificationAnalysis.idXML";
      }

      // Read the post-calibration value captured during the internal search()
      // call. Do NOT re-compute here — by now search()'s restore_fi_params() has
      // reset the tolerance members to user-configured values.
      result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
        result.peptide_ids,
        last_mod_match_tolerance_used_,
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        output_file
      );
    }
    else
    {
      OPENMS_LOG_INFO << "[ProSE] Closed search mode - modification analysis skipped" << std::endl;
    }

    return result;
  }

  // =====================================================================
  // Multi-file searchWithModificationAnalysis (in-memory FASTA).
  //
  // Builds the FragmentIndex once via prepareContext() and reuses it across
  // all input spectrum files. Each input file produces its own SearchResult
  // (with per-file modification analysis); a final aggregate SearchResult is
  // computed by pooling all per-file PSMs and running modification analysis
  // once on the pooled set. build_pooled_aggregate == false suppresses keeping
  // that pooled PSM copy (see the header for the exact contract).
  // =====================================================================
  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      const std::vector<std::string>& in_spectra_files,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    if (!output_base_names.empty() && output_base_names.size() != in_spectra_files.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "output_base_names must be empty or have exactly one entry per spectrum file (got "
        + StringUtils::toStr(output_base_names.size()) + " entries for " + StringUtils::toStr(in_spectra_files.size())
        + " spectrum files).");
    }

    MultiFileSearchResult mfres;
    mfres.aggregate.is_open_search = isOpenSearchMode_();

    if (in_spectra_files.empty())
    {
      mfres.aggregate.exit_code = ExitCodes::INPUT_FILE_EMPTY;
      return mfres;
    }

    // Multi-threaded: a helper thread reads the spectra of an mzML file while this thread works, the
    // first file from here on while the database is prepared and the index built, and, in the
    // unchunked search below, file i + 1 while file i is searched. The OpenMP threads of this thread
    // are started first, as in search(file).
    // Registry order: idXML writes the UserParams of an object in the order in which their names were
    // registered. The mzML reader registers names (CV terms and user parameters of the file, and its
    // own, see MzMLHandler) while it runs. Preparing the database and building the index register
    // none, so with the first file the registry ends up as if the file had been read first. While
    // file i is searched, the search registers the names of what it writes (scan_index, the PSM
    // features, target_decoy, PeptideIndexer:*, spectra_data, ...); if the reader of file i + 1
    // registers names in between, these keep their order relative to each other. ProSE writes none
    // of the reader's names, unless a file carries a user parameter that is named like one of them.
    const Size threads_started = startOpenMPThreads();
    const bool read_in_background = threads_started > 1;
    const auto is_mzml = [&in_spectra_files](Size i) { return FileHandler::getTypeByFileName(in_spectra_files[i]) == FileTypes::MZML; };
    std::future<PeakMap> next_spectra; // the spectra of the next file to search, if read in the background
    const auto read_next = [&](Size i)
    {
      if (read_in_background && i < in_spectra_files.size() && is_mzml(i))
      {
        next_spectra = loadMS2SpectraAsync(in_spectra_files[i]);
      }
    };
    read_next(0);

    // Resolve decoy handling once from the shared input FASTA; reused for the
    // chunk-major path, the single-context path, and the downstream FDR steps.
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    // Surface the effective marker so a caller doing merged-PSM protein FDR
    // (e.g. ProSE -out_merged) recognises the same decoys that were searched.
    mfres.decoy_string = strategy.decoy_string;
    mfres.decoy_is_prefix = strategy.is_prefix;
    mfres.have_decoys = strategy.have_decoys;

    // -- Shared configuration recap for the end-of-search report (built once) --
    {
      SharedSearchStats& sh = mfres.shared;
      sh.enzyme = enzyme_;
      sh.precursor_tol_lower = precursor_mass_tolerance_lower_;
      sh.precursor_tol_upper = precursor_mass_tolerance_upper_;
      sh.precursor_tol_unit = precursor_mass_tolerance_unit_;
      sh.fragment_tol = fragment_mass_tolerance_;
      sh.fragment_tol_unit = fragment_mass_tolerance_unit_;
      sh.min_charge = static_cast<Int>(precursor_min_charge_);
      sh.max_charge = static_cast<Int>(precursor_max_charge_);
      sh.missed_cleavages = peptide_missed_cleavages_;
      sh.fixed_mods.assign(modifications_fixed_.begin(), modifications_fixed_.end());
      sh.variable_mods.assign(modifications_variable_.begin(), modifications_variable_.end());
      if (add_a_ions_) sh.ion_series.push_back("a");
      if (add_b_ions_) sh.ion_series.push_back("b");
      if (add_c_ions_) sh.ion_series.push_back("c");
      if (add_x_ions_) sh.ion_series.push_back("x");
      if (add_y_ions_) sh.ion_series.push_back("y");
      if (add_z_ions_) sh.ion_series.push_back("z");
      if (add_zp1_ions_) sh.ion_series.push_back("z+1");
      if (ions_by_activation_) sh.ion_series.push_back("+c/z+1 for ETD/ECD/EThcD/ETciD");
      sh.open_search = isOpenSearchMode_();
      sh.calibration_enabled = calibration_enabled_;
      sh.psm_fdr_threshold = fdr_psm_;
      sh.protein_fdr_threshold = fdr_protein_;
      // Decoy handling resolved via the auto/generate/ignore strategy (#9634):
      // target-only when no decoys are present, "generated" when ProSE synthesises
      // them, "external" when pre-existing decoys in the FASTA are reused.
      sh.decoy_mode = !strategy.have_decoys ? "none (target-only)"
                      : (strategy.generate ? "generated" : "external");
    }

    // Decide chunking on the decoy-augmented DB size (#9180): if decoy generation
    // doubles the target DB, a 3000-protein target against chunk_size=5000 should
    // still chunk because the augmented DB is 6000 — otherwise the resulting FI
    // would exceed the user's memory budget by 2×.
    std::vector<FASTAFile::FASTAEntry> full_db;
    bool use_chunked = false;
    if (database_chunk_size_ > 0)
    {
      full_db = buildDecoyAugmentedDB_(fasta_db, strategy);
      use_chunked = (full_db.size() > database_chunk_size_);
    }

    if (use_chunked)
    {
      // ================================================================
      // Chunk-major multi-file path (MSFragger-style): build each chunk's
      // FragmentIndex ONCE and score ALL files against it before moving
      // to the next chunk. This gives C FI builds instead of N×C.
      // ================================================================
      // full_db already built above.
      const bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
      const bool open_search_mode = isOpenSearchMode_();

      // Shared report stats: chunk-major builds C indices (one per chunk),
      // reused across all files — count/time them once, not per file.
      mfres.shared.chunked = true;
      for (const auto& e : full_db)
      {
        // Count by the RESOLVED decoy marker (prefix OR suffix), not the hardcoded
        // decoy_prefix_: otherwise reused external/suffix decoys (decoy_mode "external")
        // would be miscounted as targets. have_decoys is false for target-only (ignore),
        // where decoys are stripped and decoy_string is empty.
        const bool is_decoy = strategy.have_decoys &&
            accessionHasDecoyMarker_(e.identifier, strategy.decoy_string, strategy.is_prefix);
        if (is_decoy) { ++mfres.shared.db_decoy_proteins; }
        else { ++mfres.shared.db_target_proteins; }
      }

      OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search_mode ? "true" : "false")
                      << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                      << ", +" << precursor_mass_tolerance_upper_ << "] "
                      << precursor_mass_tolerance_unit_ << ")" << std::endl;

      // Phase 1: Load + preprocess all files.
      std::vector<PeakMap> all_spectra(in_spectra_files.size());
      // annotate:self_trained_ion_priors: the compact peak lists of every file, until its PSMs are annotated
      std::vector<FragmentIonLikelihoodModel::PeakLists> all_ion_evidence(in_spectra_files.size());
      const bool retain_evidence = param_.getValue("annotate:local_fragment_evidence").toBool();
      std::vector<PeakMap> all_evidence_spectra(retain_evidence ? in_spectra_files.size() : 0);
      std::vector<PeakMap> all_query_spectra(query_raw_spectrum_ ? in_spectra_files.size() : 0);
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        OPENMS_LOG_INFO << "[ProSE] Loading " << in_spectra_files[i] << std::endl;
        all_spectra[i] = next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra_files[i], readerThreads(false));
        preprocessSpectra_(all_spectra[i], fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_,
                           peaks_keep_n_, peaks_window_top_, peaks_window_type_, deisotoping_,
                           self_trained_ion_priors_ ? &all_ion_evidence[i] : nullptr, ion_prior_scored_peaks_,
                           retain_evidence ? &all_evidence_spectra[i] : nullptr,
                           query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
      }

      // ions:by_activation: the chunk indices, shared by all files, hold c and z+1 ions if any file has
      // electron-activated spectra. Only those spectra are matched against them (see
      // scoreSpectraAgainstIndex_), so the other files get the results they would get alone.
      bool electron_ions = false;
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const Size n_electron_activated = countElectronActivated_(all_spectra[i]);
        if (n_electron_activated == 0) continue;
        electron_ions = true;
        OPENMS_LOG_INFO << "[ProSE] " << in_spectra_files[i] << ": " << n_electron_activated << " of "
                        << all_spectra[i].size() << " spectra are electron-activated (ETD, ECD, EThcD or ETciD)"
                        << " and are also scored with c and z+1 ions." << std::endl;
      }

      // Per-file calibration: build a strided calibration FI once, run calibration per file.
      // Stores per-file effective tolerances (asymmetric lower/upper preserved — see #9180)
      // for use during scoring.
      const Size chunk_size = database_chunk_size_;
      const Size n_chunks = (full_db.size() + chunk_size - 1) / chunk_size;
      struct PerFileCalibration
      {
        double effective_precursor_tol_lower;
        double effective_precursor_tol_upper;
        double effective_fragment_tol;
        double mod_match_tol;  // for open-search mod analysis
      };
      std::vector<PerFileCalibration> per_file_cal(in_spectra_files.size());

      // Default: user-configured tolerances.
      for (auto& cal : per_file_cal)
      {
        cal.effective_precursor_tol_lower = precursor_mass_tolerance_lower_;
        cal.effective_precursor_tol_upper = precursor_mass_tolerance_upper_;
        cal.effective_fragment_tol = fragment_mass_tolerance_;
        cal.mod_match_tol = computeModMatchTolerance_();
      }

      if (calibration_enabled_ && !open_search_mode)
      {
        // Build a strided-sample calibration FI once, reused across files.
        std::vector<FASTAFile::FASTAEntry> cal_db = buildCalibrationSample_(full_db);
        FragmentIndex cal_fi;
        cal_fi.setParameters(fragmentIndexParameters_(electron_ions));
        StopWatch sw_cal_idx; sw_cal_idx.start();
        cal_fi.build(cal_db);
        sw_cal_idx.stop();
        mfres.shared.seconds_index_build += sw_cal_idx.getClockTime();

        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
          OPENMS_LOG_INFO << "[ProSE] Calibration for " << in_spectra_files[i]
                          << " (strided sample, " << cal_db.size() << " proteins)" << std::endl;
          CalibrationResult_ cal = runCalibrationPass_(all_spectra[i], cal_fi, cal_db, query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
          if (cal.success)
          {
            per_file_cal[i].effective_fragment_tol = cal.fragment_tolerance;
            if (!cal.extreme_bias)
            {
              per_file_cal[i].effective_precursor_tol_lower = cal.cal_lower;
              per_file_cal[i].effective_precursor_tol_upper = cal.cal_upper;
              OPENMS_LOG_INFO << "[ProSE] Calibration: shift=" << cal.precursor_shift
                              << " " << precursor_mass_tolerance_unit_
                              << " -> window [-" << cal.cal_lower << ", +" << cal.cal_upper << "]"
                              << " fragment=" << cal.fragment_tolerance << std::endl;
            }
            else
            {
              OPENMS_LOG_WARN << "[ProSE] Calibration for " << in_spectra_files[i]
                              << ": extreme bias, precursor calibration discarded. Fragment calibration applied." << std::endl;
            }
            // Recompute mod-match tolerance with calibrated values.
            // Temporarily set member variables, compute, then restore.
            const double orig_lower = precursor_mass_tolerance_lower_;
            const double orig_upper = precursor_mass_tolerance_upper_;
            if (!cal.extreme_bias)
            {
              precursor_mass_tolerance_lower_ = cal.cal_lower;
              precursor_mass_tolerance_upper_ = cal.cal_upper;
            }
            per_file_cal[i].mod_match_tol = computeModMatchTolerance_();
            precursor_mass_tolerance_lower_ = orig_lower;
            precursor_mass_tolerance_upper_ = orig_upper;
          }
          else
          {
            OPENMS_LOG_INFO << "[ProSE] Calibration failed for " << in_spectra_files[i]
                            << ", using configured tolerances." << std::endl;
          }
        }
        // cal_fi freed here.
      }
      else if (calibration_enabled_ && open_search_mode)
      {
        OPENMS_LOG_WARN << "Warning: calibration not applicable in open-search mode." << std::endl;
      }

      // Prepare spectrum generators (once).
      const SpectrumGenerators_ generators = spectrumGenerators_();

      // Per-file hit accumulators.
      std::vector<std::vector<std::vector<AnnotatedHit_>>> per_file_hits(in_spectra_files.size());
      // One pool summary per spectrum per file; accumulates across chunks.
      std::vector<std::vector<CandidatePoolStats_>> per_file_pool_stats(in_spectra_files.size());
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        per_file_hits[i].resize(all_spectra[i].size());
        for (auto& a : per_file_hits[i]) a.reserve(report_top_hits_);
        per_file_pool_stats[i].resize(all_spectra[i].size());
      }
      OPENMS_LOG_INFO << "[ProSE] Chunk-major multi-file: " << full_db.size()
                      << " proteins, " << n_chunks << " chunks, "
                      << in_spectra_files.size() << " files." << std::endl;

      Size chunk_idx = 0;
      for (Size start = 0; start < full_db.size(); start += chunk_size)
      {
        ++chunk_idx;
        const Size end = std::min(start + chunk_size, full_db.size());
        OPENMS_LOG_INFO << "[ProSE] Chunk " << chunk_idx << "/" << n_chunks
                        << " (" << (end - start) << " proteins)" << std::endl;

        std::vector<FASTAFile::FASTAEntry> chunk_db(full_db.begin() + start, full_db.begin() + end);
        FragmentIndex chunk_fi;
        chunk_fi.setParameters(fragmentIndexParameters_(electron_ions));
        StopWatch sw_chunk; sw_chunk.start();
        chunk_fi.build(chunk_db);
        sw_chunk.stop();
        mfres.shared.seconds_index_build += sw_chunk.getClockTime();
        mfres.shared.indexed_peptides += chunk_fi.getPeptides().size();
        mfres.shared.indexed_fragments += chunk_fi.getNumFragments();
        if (chunk_fi.isSnesMode()) { mfres.shared.snes_mode = true; }

        // Score ALL files against this chunk's index.
        // Each file may have different calibrated tolerances — apply per-file
        // precursor bounds to the FI before scoring, then use per-file
        // fragment tolerance for HyperScore.
        const Param base_fi_params = chunk_fi.getParameters();
        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
          // Apply per-file calibrated precursor bounds to FI query params. Asymmetric
          // lower/upper preserved — collapsing to max() would re-open the tight side
          // of the calibrated window and admit spurious decoy candidates (#9180).
          if (calibration_enabled_ && !open_search_mode)
          {
            Param fi_params = base_fi_params;
            fi_params.setValue("fragment:mass_tolerance", per_file_cal[i].effective_fragment_tol);
            fi_params.setValue("precursor:mass_tolerance_lower", per_file_cal[i].effective_precursor_tol_lower);
            fi_params.setValue("precursor:mass_tolerance_upper", per_file_cal[i].effective_precursor_tol_upper);
            chunk_fi.setParameters(fi_params);
          }
          scoreSpectraAgainstIndex_(all_spectra[i], chunk_fi, chunk_db,
                                    generators, per_file_cal[i].effective_fragment_tol,
                                    fragment_mass_tolerance_unit_ppm, open_search_mode,
                                    per_file_hits[i], per_file_pool_stats[i],
                                    "  file " + StringUtils::toStr(i + 1) + " chunk " + StringUtils::toStr(chunk_idx),
                                    query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
        }
        // Restore base FI params for next chunk (in case calibration modified them).
        if (calibration_enabled_ && !open_search_mode)
          chunk_fi.setParameters(base_fi_params);

        // Per-chunk pruning for each file.
        const Size keep = std::max(report_top_hits_, Size(2));
        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
#pragma omp parallel for default(none) shared(per_file_hits, i, keep)
          for (SignedSize si = 0; si < (SignedSize)per_file_hits[i].size(); ++si)
          {
            if (per_file_hits[i][si].size() > keep)
            {
              std::partial_sort(per_file_hits[i][si].begin(),
                                per_file_hits[i][si].begin() + keep,
                                per_file_hits[i][si].end(),
                                AnnotatedHit_::hasBetterScore);
              per_file_hits[i][si].resize(keep);
            }
          }
        }
      } // end chunk loop

      // Phase 3: Per-file postprocess with per-file calibrated tolerances.
      for (auto& file_stats : per_file_pool_stats)
      {
        for (auto& stats : file_stats)
        {
          std::unordered_set<std::string>().swap(stats.seen_candidates);
        }
      }
      mfres.per_file.reserve(in_spectra_files.size());

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");
        last_mod_match_tolerance_used_ = per_file_cal[i].mod_match_tol;

        SearchResult result;
        result.is_open_search = open_search_mode;

        if (query_raw_spectrum_) all_query_spectra[i].clear(true);

        postProcessHits_(all_spectra[i], per_file_hits[i], per_file_pool_stats[i],
          result.protein_ids, result.peptide_ids,
          report_top_hits_, modifications_fixed_, modifications_variable_,
          peptide_missed_cleavages_,
          std::max(per_file_cal[i].effective_precursor_tol_lower,
                   per_file_cal[i].effective_precursor_tol_upper),
          per_file_cal[i].effective_fragment_tol,
          precursor_mass_tolerance_unit_, fragment_mass_tolerance_unit_,
          precursor_min_charge_, precursor_max_charge_, enzyme_, "",
          retain_evidence ? &all_evidence_spectra[i] : nullptr);
        if (retain_evidence) all_evidence_spectra[i].clear(true);

        PeptideIndexing indexer;
        Param param_pi = indexer.getParameters();
        param_pi.setValue("decoy_string", strategy.decoy_string.empty() ? decoy_prefix_ : strategy.decoy_string);
        param_pi.setValue("decoy_string_position", strategy.is_prefix ? "prefix" : "suffix");
        param_pi.setValue("enzyme:name", enzyme_);
        param_pi.setValue("enzyme:specificity",
                          EnzymaticDigestion::NamesOfSpecificity[peptide_enzyme_specificity_]);
        param_pi.setValue("missing_decoy_action", "silent");
        indexer.setParameters(param_pi);
        PeptideIndexing::ExitCodes indexer_exit =
            indexer.run(full_db, result.protein_ids, result.peptide_ids);
        if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
            (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
        {
          OPENMS_LOG_WARN << "[ProSE] PeptideIndexing failed for " << in_spectra
                          << " (exit code " << static_cast<int>(indexer_exit)
                          << "). Skipping this file." << std::endl;
          if (indexer_exit == PeptideIndexing::ExitCodes::DATABASE_EMPTY)
            result.exit_code = ExitCodes::INPUT_FILE_EMPTY;
          else if (indexer_exit == PeptideIndexing::ExitCodes::UNEXPECTED_RESULT)
            result.exit_code = ExitCodes::UNEXPECTED_RESULT;
          else
            result.exit_code = ExitCodes::UNKNOWN_ERROR;
          result.stats.input_file = File::basename(in_spectra);
          mfres.per_file.push_back(std::move(result));
          all_ion_evidence[i].clear();
          continue;
        }

        // PSM-level FDR only (matching non-chunked search semantics). Target+decoy PSMs are
        // filtered alike by the q-value threshold; no decoy-specific stripping happens here
        // (decoupled from decoy removal). The decoys that pass are kept because per-file results
        // feed cross-file merging and protein-level FDR. Categorical decoy removal is a
        // protein-FDR finalization step performed by ProSE.cpp. have_decoys is the marker-aware
        // result from resolveDecoyStrategy_ (prefix or suffix, generated or external).
        const bool has_decoys = strategy.have_decoys;
        // Optional per-run ion priors, after PeptideIndexing and before FDR (see search()); each file learns its own.
        if (self_trained_ion_priors_)
        {
          annotateIonPriors_(all_spectra[i], all_ion_evidence[i], result.protein_ids, result.peptide_ids);
          all_ion_evidence[i].clear();
        }
        // Pre-FDR stats (target/decoy counts + HyperScore distribution).
        capturePreFdrStats_(result.peptide_ids, result.stats);
        if (fdr_psm_ > 0.0 && has_decoys)
        {
          StopWatch sw_fdr; sw_fdr.start();
          annotatePsmQValues(result.peptide_ids);
          IDFilter::filterHitsByScore(result.peptide_ids, fdr_psm_);
          result.stats.fdr_applied = true;
          result.stats.achieved_psm_fdr = maxRetainedScore_(result.peptide_ids);
          sw_fdr.stop();
          result.stats.seconds_fdr = sw_fdr.getClockTime();
        }

        result.exit_code = ExitCodes::EXECUTION_OK;

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, all_spectra[i]);

        // In chunk-major mode scoring is shared across files (one pass per chunk),
        // so per-file seconds_search is not separable; it is left 0 and the shared
        // index/total timing carries the cost.
        collectRunStatistics_(all_spectra[i], result.protein_ids, result.peptide_ids, result.stats);
        result.stats.input_file = File::basename(in_spectra);

        // Per-file modification analysis (uses per-file calibrated tolerance).
        if (result.is_open_search)
        {
          OpenSearchModificationAnalysis mod_analyzer;
          std::string output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, per_file_cal[i].mod_match_tol,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
        }

        mfres.per_file.push_back(std::move(result));
      }
    }
    else
    {
      // ================================================================
      // Non-chunked multi-file: shared SearchContext (existing path).
      // ================================================================
      SearchContext ctx;
      bool ctx_built = false;
      // ions:by_activation: the shared index gets c and z+1 ions once a file has electron-activated
      // spectra, by rebuilding it. Only such spectra are matched against these ions (see
      // scoreSpectraAgainstIndex_), so the results of the other files depend neither on whether
      // the index holds them nor on the input order, and no file has to be read in advance.
      // A single file is the only search of the context: its index holds only the peptides that the
      // spectra (read before) can reach.
      auto prepare_context = [&](bool electron_ions, const std::function<const PeakMap*(Size)>& searched_spectra)
      {
        if (ctx_built && (ctx.electron_ions || !electron_ions)) { return; }
        StopWatch sw_idx; sw_idx.start();
        if (ctx_built)
        {
          startProgress(0, 1, "Building fragment index with c and z+1 ions...");
          ctx.fragment_index.clear();
          ctx.fragment_index.setParameters(fragmentIndexParameters_(true));
          ctx.fragment_index.build(ctx.db);
          ctx.electron_ions = true;
          endProgress();
        }
        else
        {
          // As prepareContext(fasta_db, electron_ions), but with the decoy strategy resolved above
          // instead of detecting the decoys of fasta_db a second time. If chunk_size was set but the
          // augmented DB fits in one chunk, full_db holds it already.
          if (full_db.empty())
          {
            startProgress(0, 1, "Generate decoys...");
            full_db = buildDecoyAugmentedDB_(fasta_db, strategy);
            endProgress();
          }
          ctx.db = std::move(full_db);
          ctx.decoy_string = strategy.decoy_string;
          ctx.decoy_is_prefix = strategy.is_prefix;
          ctx.have_decoys = strategy.have_decoys;
          startProgress(0, 1, "Building fragment index...");
          ctx.fragment_index.setParameters(fragmentIndexParameters_(electron_ions));
          ctx.fragment_index.build(ctx.db, searched_spectra);
          ctx.electron_ions = electron_ions;
          endProgress();
        }
        sw_idx.stop();

        // Shared report stats: index built once (and rebuilt at most once, with c and z+1 ions)
        // and reused across all files.
        mfres.shared.chunked = false;
        mfres.shared.seconds_index_build += sw_idx.getClockTime();
        mfres.shared.indexed_fragments = ctx.fragment_index.getNumFragments();
        if (ctx_built) { return; } // a rebuild: same database and peptides
        ctx_built = true;
        mfres.shared.indexed_peptides = ctx.fragment_index.getPeptides().size();
        mfres.shared.snes_mode = ctx.fragment_index.isSnesMode();
        for (const auto& e : ctx.db)
        {
          // Count by the RESOLVED decoy marker (prefix OR suffix), not the hardcoded
          // decoy_prefix_: otherwise reused external/suffix decoys (decoy_mode "external")
          // would be miscounted as targets. have_decoys is false for target-only (ignore),
          // where decoys are stripped and decoy_string is empty.
          const bool is_decoy = strategy.have_decoys &&
              accessionHasDecoyMarker_(e.identifier, strategy.decoy_string, strategy.is_prefix);
          if (is_decoy) { ++mfres.shared.db_decoy_proteins; }
          else { ++mfres.shared.db_target_proteins; }
        }
      };

      mfres.per_file.reserve(in_spectra_files.size());

      // While the first file is read: build the index, as search(file) does without c and z+1 ions
      // unless the spectra are read already or the file names an electron-based activation early on.
      // prepare_context() below adds them if the spectra need them after all.
      // A single file is the only search of the context (r3 merge of WP7 and WP10): as in search(file), the index
      // holds only the peptides its spectra can reach if they are read by the time the peptides are generated, or
      // if waiting for them pays off (see MAX_THREADS_WAITING_FOR_SPECTRA); the spectra taken for that are searched.
      const bool restrict_index = in_spectra_files.size() == 1 && restrictIndexToSpectra_();
      PeakMap first_spectra;
      bool first_spectra_taken = false;
      if (next_spectra.valid() && next_spectra.wait_for(std::chrono::seconds(0)) != std::future_status::ready
          && !(ions_by_activation_ && mzMLHeadNamesElectronActivation(in_spectra_files[0])))
      {
        std::function<const PeakMap*(Size)> searched_spectra;
        if (restrict_index)
        {
          const UInt64 spectra_bytes = File::fileSize(in_spectra_files[0]); // UInt64(-1) if unknown: no waiting
          searched_spectra = [&, spectra_bytes](Size peptides) -> const PeakMap*
          {
            if (!first_spectra_taken)
            {
              const bool wait = threads_started <= MAX_THREADS_WAITING_FOR_SPECTRA
                                && spectra_bytes <= MAX_SPECTRA_BYTES_PER_PEPTIDE * peptides / threads_started;
              if (!wait && next_spectra.wait_for(std::chrono::seconds(0)) != std::future_status::ready)
              {
                return nullptr;
              }
              first_spectra = next_spectra.get();
              first_spectra_taken = true;
            }
            return &first_spectra;
          };
        }
        prepare_context(false, searched_spectra);
      }

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");

        OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                        << "] Searching " << in_spectra << std::endl;

        PeakMap spectra = first_spectra_taken ? std::move(first_spectra)
                          : next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra, readerThreads(false));
        first_spectra_taken = false;
        read_next(i + 1);
        std::function<const PeakMap*(Size)> searched_spectra;
        if (restrict_index) { searched_spectra = [&spectra](Size) { return &spectra; }; }
        prepare_context(countElectronActivated_(spectra) > 0, searched_spectra);

        SearchResult result;
        result.is_open_search = isOpenSearchMode_();
        // No file after this one: search() releases the index right after scoring (in parallel, and
        // before the post-processing allocates) instead of the context's destructor at the end.
        ctx.release_fragment_index_after_scoring = i + 1 == in_spectra_files.size();
        result.exit_code = search(spectra, ctx, result.protein_ids, result.peptide_ids);

        if (result.exit_code != ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_WARN << "[ProSE] Search failed for " << in_spectra
                          << " (exit code " << static_cast<int>(result.exit_code) << "). Continuing." << std::endl;
          result.stats = last_run_stats_;
          result.stats.input_file = File::basename(in_spectra);
          mfres.per_file.push_back(std::move(result));
          continue;
        }

        result.stats = last_run_stats_;
        result.stats.input_file = File::basename(in_spectra);

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);

        if (result.is_open_search)
        {
          OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for " << in_spectra << std::endl;
          OpenSearchModificationAnalysis mod_analyzer;
          std::string output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, last_mod_match_tolerance_used_,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
        }
        else
        {
          OPENMS_LOG_INFO << "[ProSE] Closed search mode - per-file modification analysis skipped" << std::endl;
        }

        mfres.per_file.push_back(std::move(result));
      }
    }

    // Build the aggregate result by pooling per-file PSMs.
    // If every per-file run failed, propagate the first non-OK exit code into the
    // aggregate (rather than overloading UNEXPECTED_RESULT, which has a specific
    // meaning - PeptideIndexing failure).
    bool any_ok = false;
    for (const auto& pf : mfres.per_file)
    {
      if (pf.exit_code == ExitCodes::EXECUTION_OK) { any_ok = true; break; }
    }

    if (!any_ok)
    {
      for (const auto& pf : mfres.per_file)
      {
        if (pf.exit_code != ExitCodes::EXECUTION_OK)
        {
          mfres.aggregate.exit_code = pf.exit_code;
          break;
        }
      }
      return mfres;
    }

    // Single-file fast path: the aggregate would just duplicate the only per-file
    // result and re-run modification analysis on the same PSMs. Skip the pooling
    // and analysis entirely. The aggregate is left with only @c is_open_search
    // and @c exit_code populated; callers should use @c per_file[0] for the
    // actual identifications. This is documented on @c MultiFileSearchResult.
    if (mfres.per_file.size() == 1 && mfres.per_file[0].exit_code == ExitCodes::EXECUTION_OK)
    {
      return mfres;
    }

    // Pooling copies every per-file PSM a second time. Skip it entirely when the caller
    // does not want the pooled set AND nothing here needs it: only the open-search
    // aggregate modification analysis below consumes it, so in closed search a
    // build_pooled_aggregate == false caller gets the same almost-empty aggregate as the
    // single-file fast path above. In open search the pooled set is built, analyzed and
    // then released again (see the release block after the analysis).
    const bool need_pooled = build_pooled_aggregate || mfres.aggregate.is_open_search;

    if (need_pooled)
    {
      // Merge per-file identifications into a single aggregate using
      // IDMergerAlgorithm — the canonical OpenMS pattern for cross-file protein
      // inference (mirrors ConsensusMapMergerAlgorithm::mergeAllIDRuns in
      // ProteomicsLFQ). This deduplicates ProteinHits by accession (union) and
      // remaps all PeptideIdentification identifiers to the merged run. No
      // PeptideIndexing re-run needed: per-file searches already linked
      // peptides to proteins.
      // The COPY overload is mandatory here: per_file is returned to the caller and must
      // stay intact, so nothing may be moved out of it.
      IDMergerAlgorithm merger;
      for (const auto& pf : mfres.per_file)
      {
        if (pf.exit_code != ExitCodes::EXECUTION_OK) { continue; }
        merger.insertRuns(pf.protein_ids, pf.peptide_ids);
      }
      ProteinIdentification merged_proteins;
      merger.returnResultsAndClear(merged_proteins, mfres.aggregate.peptide_ids);
      mfres.aggregate.protein_ids = {std::move(merged_proteins)};
      mfres.aggregate.protein_ids[0].setPrimaryMSRunPath(in_spectra_files);
    }

    // Note: protein inference + FDR on the aggregate is left to the caller
    // (e.g., the TOPP tool's -out_merged option) so that per-file outputs
    // retain run-level information without being overwritten by aggregate
    // protein lists. The aggregate here contains the merged (unfiltered)
    // proteins + pooled PSMs for downstream use.

    // Aggregate modification analysis on the pooled PSM set.
    if (mfres.aggregate.is_open_search && !mfres.aggregate.peptide_ids.empty())
    {
      OPENMS_LOG_INFO << "[ProSE] Running aggregate modification analysis on "
                      << mfres.aggregate.peptide_ids.size() << " pooled PSM(s) from "
                      << in_spectra_files.size() << " input file(s)." << std::endl;

      OpenSearchModificationAnalysis mod_analyzer;
      std::string agg_output_file;
      if (!aggregate_base_name.empty())
      {
        agg_output_file = aggregate_base_name + "_ModificationAnalysis.idXML";
      }

      mfres.aggregate.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
        mfres.aggregate.peptide_ids,
        computeModMatchTolerance_(),
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        agg_output_file
      );
    }

    // The pooled PSMs were only needed as input to the analysis above; release the second
    // copy again so a caller that asked for no pooled aggregate does not pay for it. The
    // analysis result stays. Leaves the aggregate exactly as in the single-file fast path
    // (only is_open_search / exit_code populated), as documented on MultiFileSearchResult.
    if (!build_pooled_aggregate)
    {
      mfres.aggregate.peptide_ids.clear();
      mfres.aggregate.peptide_ids.shrink_to_fit();
      mfres.aggregate.protein_ids.clear();
      mfres.aggregate.protein_ids.shrink_to_fit();
    }

    return mfres;
  }

  // =====================================================================
  // Multi-file searchWithModificationAnalysis (FASTA file path).
  //
  // Loads the FASTA from disk once, delegates to the in-memory multi-file
  // overload, and patches the database file path on each per-file (and the
  // aggregate) ProteinIdentification.
  // =====================================================================
  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      const std::vector<std::string>& in_spectra_files,
      const std::string& in_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    // load FASTA once
    vector<FASTAFile::FASTAEntry> fasta_db;
    loadFASTA(in_db, fasta_db);

    MultiFileSearchResult mfres = searchWithModificationAnalysis(
      in_spectra_files, fasta_db, output_base_names, aggregate_base_name, build_pooled_aggregate);

    mfres.shared.database_file = in_db;

    // Patch database file path into each per-file and aggregate ProteinIdentification.
    for (auto& pf : mfres.per_file)
    {
      if (pf.exit_code == ExitCodes::EXECUTION_OK && !pf.protein_ids.empty())
      {
        pf.protein_ids[0].getSearchParameters().db = in_db;
      }
    }
    if (mfres.aggregate.exit_code == ExitCodes::EXECUTION_OK && !mfres.aggregate.protein_ids.empty())
    {
      mfres.aggregate.protein_ids[0].getSearchParameters().db = in_db;
    }

    return mfres;
  }

  // =====================================================================
  // File-based single-file searchWithModificationAnalysis: thin wrapper around
  // the multi-file overload using a single-element list.
  // =====================================================================
  ProSEAlgorithm::SearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(const std::string& in_spectra,
                                                                  const std::string& in_db,
                                                                  const std::string& output_base_name) const
  {
    std::vector<std::string> in_files{in_spectra};
    std::vector<std::string> base_names;
    if (!output_base_name.empty()) { base_names.push_back(output_base_name); }

    MultiFileSearchResult mfres = searchWithModificationAnalysis(in_files, in_db, base_names, "");

    if (mfres.per_file.empty())
    {
      SearchResult empty_result;
      empty_result.exit_code = ExitCodes::INPUT_FILE_EMPTY;
      return empty_result;
    }

    return std::move(mfres.per_file[0]);
  }

  // =====================================================================
  // Helper: decoy-marker membership test (prefix/suffix), pure std::string.
  // =====================================================================
  bool ProSEAlgorithm::accessionHasDecoyMarker_(const std::string& accession,
                                                const std::string& marker, bool is_prefix)
  {
    if (marker.empty() || accession.size() < marker.size()) { return false; }
    return is_prefix ? (accession.compare(0, marker.size(), marker) == 0)
                     : (accession.compare(accession.size() - marker.size(), marker.size(), marker) == 0);
  }

  // =====================================================================
  // Helper: capture pre-FDR statistics (target/decoy counts + HyperScore
  // distribution). Must run BEFORE FalseDiscoveryRate, which overwrites the
  // hit score with the q-value and may drop decoy hits.
  // =====================================================================
  void ProSEAlgorithm::capturePreFdrStats_(const PeptideIdentificationList& peptide_ids,
                                           RunStatistics& stats)
  {
    stats.target_psms = 0;
    stats.decoy_psms = 0;
    std::vector<double> hyperscores;
    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      const PeptideHit& top = pid.getHits().front();
      hyperscores.push_back(top.getScore());

      if (!top.metaValueExists(Constants::UserParam::TARGET_DECOY)) { continue; }
      const std::string td = top.getMetaValue(Constants::UserParam::TARGET_DECOY).toString();
      // "target+decoy" counts as target (OpenMS FDR semantics).
      if (td == "decoy") { ++stats.decoy_psms; }
      else if (!td.empty()) { ++stats.target_psms; }
    }

    if (!hyperscores.empty())
    {
      auto [mn, mx] = std::minmax_element(hyperscores.begin(), hyperscores.end());
      stats.hyperscore_min = *mn;
      stats.hyperscore_max = *mx;
      stats.hyperscore_median = Math::median(hyperscores.begin(), hyperscores.end());
      stats.score_stats_valid = true;
    }
  }

  // =====================================================================
  // Recompute result-level stats from a FINAL (post-rescoring/post-FDR) PSM list.
  // =====================================================================
  void ProSEAlgorithm::updateFinalStats(RunStatistics& stats,
                                        const PeptideIdentificationList& peptide_ids,
                                        const std::string& enzyme,
                                        bool fdr_applied)
  {
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    stats.matched_spectra = static_cast<Size>(std::count_if(peptide_ids.begin(), peptide_ids.end(),
      [](const PeptideIdentification& pid) { return !pid.getHits().empty(); }));
    stats.charge_histogram.clear();
    stats.missed_cleavage_histogram.clear();
    Size n_target = 0, n_decoy = 0;

    EnzymaticDigestion digestor;
    digestor.setEnzyme(ProteaseDB::getInstance()->getEnzyme(enzyme));

    set<std::string> unique_peptides, unique_proteins;
    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      const PeptideHit& top = pid.getHits().front();
      unique_peptides.insert(top.getSequence().toString());
      for (const auto& ev : top.getPeptideEvidences()) { unique_proteins.insert(ev.getProteinAccession()); }
      if (top.getCharge() > 0) { ++stats.charge_histogram[top.getCharge()]; }
      ++stats.missed_cleavage_histogram[digestor.countInternalCleavageSites(top.getSequence().toUnmodifiedString())];
      if (top.metaValueExists(Constants::UserParam::TARGET_DECOY))
      {
        const std::string td = top.getMetaValue(Constants::UserParam::TARGET_DECOY).toString();
        if (td == "decoy") { ++n_decoy; }
        else if (!td.empty()) { ++n_target; }
      }
    }
    stats.unique_peptides = unique_peptides.size();
    stats.unique_proteins = unique_proteins.size();
    stats.target_psms = n_target;
    stats.decoy_psms = n_decoy;
    stats.fdr_applied = fdr_applied;
    stats.achieved_psm_fdr = fdr_applied ? maxRetainedScore_(peptide_ids) : -1.0;
  }

  // =====================================================================
  // Helper: maximum top-hit score among retained PSMs (== achieved FDR after
  // q-value filtering). Returns -1.0 for an empty set.
  // =====================================================================
  double ProSEAlgorithm::maxRetainedScore_(const PeptideIdentificationList& peptide_ids)
  {
    double max_score = -1.0;
    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      max_score = std::max(max_score, pid.getHits().front().getScore());
    }
    return max_score;
  }

  // =====================================================================
  // Helper: fill per-run identification statistics (silent — no logging).
  // Computes target/decoy counts from the (post-FDR) hits it is given, so
  // SearchResult::stats is consistent for any caller; leaves achieved FDR and
  // timing fields untouched (captured at well-defined points in the search paths).
  // =====================================================================
  void ProSEAlgorithm::collectRunStatistics_(
      const PeakMap& spectra,
      const std::vector<ProteinIdentification>& /*protein_ids*/,
      const PeptideIdentificationList& peptide_ids,
      RunStatistics& stats) const
  {
    stats.ms2_spectra = std::count_if(spectra.begin(), spectra.end(),
                                      [](const MSSpectrum& s) { return s.getMSLevel() == 2; });
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    stats.matched_spectra = static_cast<Size>(std::count_if(peptide_ids.begin(), peptide_ids.end(),
      [](const PeptideIdentification& pid) { return !pid.getHits().empty(); }));

    set<std::string> unique_peptides;
    set<std::string> unique_proteins;
    Size n_target = 0, n_decoy = 0;

    // Per-PSM error values for tolerance estimation (top-ranked hits only)
    vector<double> precursor_errors;
    vector<double> fragment_errors;

    EnzymaticDigestion digestor;
    digestor.setEnzyme(ProteaseDB::getInstance()->getEnzyme(enzyme_));

    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) { continue; }
      const PeptideHit& top = pid.getHits().front();
      unique_peptides.insert(top.getSequence().toString());

      for (const auto& ev : top.getPeptideEvidences())
      {
        unique_proteins.insert(ev.getProteinAccession());
      }

      if (top.getCharge() > 0) { ++stats.charge_histogram[top.getCharge()]; }

      // Precursor error: always compute inline (cheap), independent of annotate:PSM.
      // Skip hits with unresolved charge (0) — getMZ() would throw.
      if (top.getCharge() > 0)
      {
        double theo_mz = top.getSequence().getMZ(top.getCharge());
        double prec_error_ppm = Math::getPPM(pid.getMZ(), theo_mz);
        precursor_errors.push_back(fabs(prec_error_ppm));
      }

      // Fragment error: use annotation metavalue if available (computing it
      // inline would require spectrum alignment which is expensive).
      if (top.metaValueExists(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM))
      {
        fragment_errors.push_back(static_cast<double>(top.getMetaValue(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM)));
      }

      ++stats.missed_cleavage_histogram[digestor.countInternalCleavageSites(top.getSequence().toUnmodifiedString())];

      // Target/decoy of the FINAL (post-FDR) hits, so SearchResult::stats stays consistent
      // for library callers that read it without calling updateFinalStats().
      if (top.metaValueExists(Constants::UserParam::TARGET_DECOY))
      {
        const std::string td = top.getMetaValue(Constants::UserParam::TARGET_DECOY).toString();
        if (td == "decoy") { ++n_decoy; }
        else if (!td.empty()) { ++n_target; }
      }
    }

    stats.unique_peptides = unique_peptides.size();
    stats.unique_proteins = unique_proteins.size();
    stats.target_psms = n_target;
    stats.decoy_psms = n_decoy;

    // -- Per-run tolerance estimation (median + 3*MAD, matching prior behaviour) --
    const Size min_psms_for_estimation = 10;
    if (precursor_errors.size() >= min_psms_for_estimation)
    {
      double med = Math::median(precursor_errors.begin(), precursor_errors.end());
      double mad = Math::MAD(precursor_errors.begin(), precursor_errors.end(), med);
      stats.prec_err_median = med;
      stats.prec_err_mad = mad;
      stats.prec_err_recommended = std::ceil(med + 3.0 * mad);
      stats.prec_tol_valid = true;
    }
    if (fragment_errors.size() >= min_psms_for_estimation)
    {
      double med = Math::median(fragment_errors.begin(), fragment_errors.end());
      double mad = Math::MAD(fragment_errors.begin(), fragment_errors.end(), med);
      stats.frag_err_median = med;
      stats.frag_err_mad = mad;
      stats.frag_err_recommended = std::ceil(med + 3.0 * mad);
      stats.frag_tol_valid = true;
    }
  }

  // =====================================================================
  // Self-trained ion priors: learn presence, intensity rank and mass error
  // of fragment ions from the confident PSMs of this run, one model per
  // half of the spectra, and annotate every hit with the other half's model.
  // =====================================================================
  // static
  AASequence ProSEAlgorithm::reversedNoiseSequence_(const AASequence& sequence)
  {
    // Reverse all but the C-terminal residue and keep every residue's modification, so the noise
    // hypothesis shares composition, mass and enzymatic C-terminus with the peptide.
    const Size n = sequence.size();
    if (n < 3) return sequence;
    AASequence reversed;
    for (Size i = n - 1; i-- > 0;) { reversed += &sequence[i]; }
    reversed += &sequence[n - 1];
    if (sequence.hasNTerminalModification()) { reversed.setNTerminalModification(sequence.getNTerminalModification()); }
    if (sequence.hasCTerminalModification()) { reversed.setCTerminalModification(sequence.getCTerminalModification()); }
    return reversed;
  }

  void ProSEAlgorithm::annotateIonPriors_(const PeakMap& spectra,
                                          const FragmentIonLikelihoodModel::PeakLists& evidence,
                                          std::vector<ProteinIdentification>& protein_ids,
                                          PeptideIdentificationList& peptide_ids) const
  {
    if (!self_trained_ion_priors_) return;
    StopWatch sw_train, sw_annotate;
    sw_train.start();
    const bool ppm = fragment_mass_tolerance_unit_ == "ppm";
    const double tolerance = fragment_mass_tolerance_;
    const SpectrumGenerators_ generators = spectrumGenerators_();
    const bool rich = ion_prior_model_ == FragmentIonLikelihoodModel::ContextSet::RICH;
    // Fragment charges of the theoretical ions: b/y ions of a z+ precursor carry at most z - 1 charges.
    const int max_fragment_charge = ion_prior_max_fragment_charge_;
    const auto fragment_charges = [max_fragment_charge](int precursor_charge) {
      return std::max(1, std::min(max_fragment_charge, precursor_charge - 1));
    };

    // 1. The spectrum of every PSM (scan_index of postProcessHits_); -1 if it has no peaks to match.
    const Size n_ids = peptide_ids.size();
    std::vector<SignedSize> scans(n_ids, -1);
    for (Size i = 0; i < n_ids; ++i)
    {
      const PeptideIdentification& pi = peptide_ids[i];
      if (pi.getHits().empty() || !pi.metaValueExists("scan_index")) continue;
      const SignedSize scan = static_cast<int>(pi.getMetaValue("scan_index"));
      if (scan >= 0 && static_cast<Size>(scan) < spectra.size() && static_cast<Size>(scan) < evidence.size()
          && evidence.peaks(static_cast<Size>(scan)) > 0)
      {
        scans[i] = scan;
      }
    }
    // Best hit by the identification's score orientation (postProcessHits_ sorts the hits, but the
    // annotation must not depend on it).
    const auto best_hit = [](const PeptideIdentification& pi) -> const PeptideHit& {
      const std::vector<PeptideHit>& hits = pi.getHits();
      const bool higher_better = pi.isHigherScoreBetter();
      Size best = 0;
      for (Size h = 1; h < hits.size(); ++h)
      {
        if (higher_better ? hits[h].getScore() > hits[best].getScore() : hits[h].getScore() < hits[best].getScore()) best = h;
      }
      return hits[best];
    };

    // 2. Training PSMs of each half (fold = scan parity): its best target hits at TDC q <= threshold over the
    //    native score, computed within the half, so the other half's labels are not used. Decoys win exact
    //    ties, as in the native yield evaluation, which keeps the selection conservative. The search has decoys
    //    (target-only searches never get here), so (D + 1) / T is an estimate even in a half without a decoy hit.
    struct Row { double score; bool decoy; Size index; };
    std::array<std::vector<Row>, 2> rows;
    for (Size i = 0; i < n_ids; ++i)
    {
      if (scans[i] < 0) continue;
      const PeptideHit& top = best_hit(peptide_ids[i]);
      if (!top.metaValueExists(Constants::UserParam::TARGET_DECOY)) continue;
      const bool decoy = top.getMetaValue(Constants::UserParam::TARGET_DECOY).toString() == "decoy";
      const Size fold = static_cast<Size>(scans[i]) % 2;
      rows[fold].push_back({peptide_ids[i].isHigherScoreBetter() ? top.getScore() : -top.getScore(), decoy, i});
    }
    std::array<std::vector<Size>, 2> training;
    for (Size fold = 0; fold < 2; ++fold)
    {
      std::vector<Row>& fold_rows = rows[fold];
      std::sort(fold_rows.begin(), fold_rows.end(), [](const Row& a, const Row& b) {
        if (a.score != b.score) return a.score > b.score;
        if (a.decoy != b.decoy) return a.decoy;
        return a.index < b.index;
      });
      std::vector<double> q(fold_rows.size(), 1.0);
      Size targets = 0, decoys = 0;
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (fold_rows[r].decoy) { ++decoys; }
        else { ++targets; }
        q[r] = static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1));
      }
      for (Size r = fold_rows.size(); r-- > 1;) { q[r - 1] = std::min(q[r - 1], q[r]); }
      for (Size r = 0; r < fold_rows.size(); ++r)
      {
        if (!fold_rows[r].decoy && q[r] <= ion_prior_train_fdr_) training[fold].push_back(fold_rows[r].index);
      }
    }
    const bool trained = training[0].size() >= ion_prior_min_psms_ && training[1].size() >= ion_prior_min_psms_
                         && std::isfinite(tolerance) && tolerance > 0.0;

    // 3. One model per half; its signal: the ions of the training PSMs, its noise: those of their reversed sequences
    //    on the same spectra. Ions are matched in parallel; the integer counts are added up in training order.
    std::array<FragmentIonLikelihoodModel, 2> models{FragmentIonLikelihoodModel(ion_prior_model_),
                                                     FragmentIonLikelihoodModel(ion_prior_model_)};
    if (trained)
    {
      for (Size fold = 0; fold < 2; ++fold)
      {
        const std::vector<Size>& psms = training[fold];
        FragmentIonLikelihoodModel& model = models[fold];
        std::vector<std::vector<FragmentIonLikelihoodModel::Ion>> signal(psms.size()), noise(psms.size());
#pragma omp parallel for schedule(dynamic, 16)
        for (SignedSize t = 0; t < static_cast<SignedSize>(psms.size()); ++t)
        {
          const Size scan = static_cast<Size>(scans[psms[t]]);
          const PeptideHit& hit = best_hit(peptide_ids[psms[t]]);
          const int charge = hit.getCharge();
          const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spectra[scan]);
          PeakSpectrum theo;
          tsg.getSpectrum(theo, hit.getSequence(), 1, fragment_charges(charge));
          model.matchIons(evidence, scan, theo, hit.getSequence(), charge, tolerance, ppm, signal[t]);
          const AASequence reversed = reversedNoiseSequence_(hit.getSequence());
          if (reversed == hit.getSequence()) continue; // no noise hypothesis (a palindrome)
          theo.clear(true);
          tsg.getSpectrum(theo, reversed, 1, fragment_charges(charge));
          model.matchIons(evidence, scan, theo, reversed, charge, tolerance, ppm, noise[t]);
        }
        for (Size t = 0; t < psms.size(); ++t)
        {
          model.addObservations(signal[t], false);
          if (!noise[t].empty()) model.addObservations(noise[t], true);
        }
        model.finalize();
      }
    }
    else
    {
      OPENMS_LOG_WARN << "[ProSE] Ion priors: " << training[0].size() << " and " << training[1].size()
                      << " confident training PSMs in the two halves of the spectra, fewer than " << ion_prior_min_psms_
                      << "; the ion_prior_* features are 0 for this file." << std::endl;
    }
    sw_train.stop();

    // 4. Annotate every hit, scoring it with the model of the other half (cross-fitting). Untrained runs get zeros
    //    so the feature columns stay complete. The names of annotate:ion_prior_features are registered here, always in
    //    this order, before the parallel loop.
    sw_annotate.start();
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const bool write_llr = ion_prior_feature_llr_, write_explained = ion_prior_feature_explained_, write_topk = ion_prior_feature_topk_;
    const UInt mv_llr = write_llr ? registry.registerName(Constants::UserParam::ION_PRIOR_LLR) : 0;
    const UInt mv_explained = write_explained ? registry.registerName(Constants::UserParam::ION_PRIOR_EXPLAINED) : 0;
    const UInt mv_topk = write_topk ? registry.registerName(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED) : 0;
    Size annotated = 0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : annotated)
    for (SignedSize i = 0; i < static_cast<SignedSize>(n_ids); ++i)
    {
      PeptideIdentification& pi = peptide_ids[i];
      const SignedSize scan = scans[i];
      const bool scorable = trained && scan >= 0;
      const FragmentIonLikelihoodModel& model = models[scorable ? 1 - static_cast<Size>(scan) % 2 : 0];
      PeakSpectrum theo;
      std::vector<FragmentIonLikelihoodModel::Ion> ions;
      for (PeptideHit& hit : pi.getHits())
      {
        FragmentIonLikelihoodModel::Features features;
        if (scorable)
        {
          const int charge = hit.getCharge();
          theo.clear(true);
          generators.forSpectrum(spectra[scan]).getSpectrum(theo, hit.getSequence(), 1, fragment_charges(charge));
          model.matchIons(evidence, static_cast<Size>(scan), theo, hit.getSequence(), charge, tolerance, ppm, ions);
          features = model.score(ions);
        }
        if (write_llr) hit.setMetaValue(mv_llr, features.log_likelihood_ratio);
        if (write_explained) hit.setMetaValue(mv_explained, features.explained_presence);
        if (write_topk) hit.setMetaValue(mv_topk, features.top_predicted_observed);
        ++annotated;
      }
    }
    sw_annotate.stop();
    OPENMS_LOG_INFO << "[ProSE] Ion priors (" << (rich ? "rich" : "basic") << " model, "
                    << (ion_prior_scored_peaks_ ? "scored" : "all") << " peaks): "
                    << (trained ? "trained on " : "not trained; ") << training[0].size() << " + " << training[1].size()
                    << " PSMs at q <= " << ion_prior_train_fdr_ << " (halves by scan parity, each scored by the other's model) in "
                    << sw_train.getClockTime() << " s; " << annotated << " hits annotated in " << sw_annotate.getClockTime()
                    << " s; peak lists " << evidence.totalPeaks() << " peaks of " << evidence.size() << " spectra, "
                    << static_cast<double>(evidence.memoryUsage()) / (1024.0 * 1024.0) << " MiB." << std::endl;

    // 5. Percolator sees the new columns; the search parameters record every setting that determines the features
    //    (model, peaks, fragment charges, training threshold and minimum) and the training set.
    if (protein_ids.empty()) return;
    ProteinIdentification::SearchParameters& search_parameters = protein_ids[0].getSearchParameters();
    StringList features;
    if (search_parameters.metaValueExists("extra_features"))
    {
      const std::string existing = search_parameters.getMetaValue("extra_features").toString();
      if (!existing.empty()) features = ListUtils::create<std::string>(existing);
    }
    StringList written;
    if (write_llr) written.push_back(Constants::UserParam::ION_PRIOR_LLR);
    if (write_explained) written.push_back(Constants::UserParam::ION_PRIOR_EXPLAINED);
    if (write_topk) written.push_back(Constants::UserParam::ION_PRIOR_TOPK_OBSERVED);
    for (const std::string& name : written)
    {
      if (std::find(features.begin(), features.end(), name) == features.end()) features.push_back(name);
    }
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
    search_parameters.setMetaValue("annotate:self_trained_ion_priors", "true");
    search_parameters.setMetaValue("annotate:ion_prior_features", ListUtils::concatenate(written, ","));
    search_parameters.setMetaValue("annotate:ion_prior_model", rich ? "rich" : "basic");
    search_parameters.setMetaValue("annotate:ion_prior_peaks", ion_prior_scored_peaks_ ? "scored" : "all");
    search_parameters.setMetaValue("annotate:ion_prior_max_fragment_charge", ion_prior_max_fragment_charge_);
    search_parameters.setMetaValue("annotate:ion_prior_train_fdr", ion_prior_train_fdr_);
    search_parameters.setMetaValue("annotate:ion_prior_min_psms", static_cast<int>(ion_prior_min_psms_));
    search_parameters.setMetaValue("ion_prior:trained", trained ? "true" : "false");
    search_parameters.setMetaValue("ion_prior:training_psms", static_cast<int>(training[0].size() + training[1].size()));
    search_parameters.setMetaValue("ion_prior:fold_training_psms",
                                   IntList{static_cast<int>(training[0].size()), static_cast<int>(training[1].size())});
  }

  // =====================================================================
  // Helper: run calibration pass on a subset of spectra
  // =====================================================================
  ProSEAlgorithm::CalibrationResult_
  ProSEAlgorithm::runCalibrationPass_(
      PeakMap& spectra,
      FragmentIndex& fragment_index,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const PeakMap* query_spectra) const
  {
    CalibrationResult_ result;

    // Calibration queries the FI and reconstructs candidate peptides to measure the
    // precursor mass error distribution. In SNES mode, a candidate is a *mother* whose
    // realized sub-peptide depends on the observed precursor — baking in a circular
    // dependency with the calibration target. v1 disables calibration in SNES mode;
    // users who need calibrated tolerances should pre-calibrate spectra upstream
    // (e.g. with InternalCalibration) before running ProSE with SNES.
    if (fragment_index.isSnesMode())
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration is not supported in SNES mode (v1.1). "
                      << "Using configured precursor tolerance unchanged.\n";
      return result; // success=false, tolerances untouched
    }

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    // Select subset by TIC (highest-quality spectra first)
    vector<pair<double, Size>> tic_index;
    tic_index.reserve(spectra.size());
    for (Size i = 0; i < spectra.size(); ++i)
    {
      double tic = 0;
      for (const auto& p : spectra[i]) { tic += p.getIntensity(); }
      tic_index.emplace_back(tic, i);
    }
    std::sort(tic_index.rbegin(), tic_index.rend());
    Size subset_size = std::max<Size>(1, static_cast<Size>(spectra.size() * calibration_subset_ratio_));
    if (subset_size > tic_index.size()) subset_size = tic_index.size();

    OPENMS_LOG_INFO << "[ProSE] Calibration: scoring " << subset_size << " / " << spectra.size()
                    << " spectra (top TIC)..." << std::endl;

    // Score subset and collect errors from the best hit per spectrum, with the ion series
    // the fragment index was built from (e.g. c/z+1 for ETD, where b/y would match nothing)
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // Collect per-spectrum best hits with scores and errors
    struct CalHit { double score; double prec_error; double frag_error; };
    vector<CalHit> cal_hits;

    // Parallelize over the calibration subset, mirroring the main scoring loop
    // (scoreSpectraAgainstIndex_). Each iteration is independent: querySpectrum and
    // the shared TheoreticalSpectrumGenerators expose const, thread-safe methods (the
    // main loop already calls them concurrently), and every working variable below is
    // loop-local. The only cross-thread write is the push into cal_hits, guarded by a
    // critical section. Without this the calibration pass ran single-threaded: on a
    // many-core machine the 10% subset took about as long as the entire parallel
    // search, roughly doubling wall time for no extra useful work. Downstream is
    // order-independent — cal_hits is sorted by score and the error vectors are sorted
    // before quantiles — so parallel insertion order does not change the result.
    // An exception on a worker thread (a modification AASequence rejects, a scorer's parameter check, ...) must not
    // leave the parallel region, which would terminate the process: the first one is rethrown after the loop.
    std::exception_ptr calibration_error;
    std::atomic<bool> calibration_failed{false};
#pragma omp parallel for schedule(dynamic)
    for (SignedSize si = 0; si < (SignedSize)subset_size; ++si)
    {
      if (calibration_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      const Size scan_idx = tic_index[si].second;
      const MSSpectrum& spec = spectra[scan_idx];
      const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);

      FragmentIndex::SpectrumMatchesTopN top_sms;
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_idx] : spec;
      fragment_index.querySpectrum(query, db, top_sms, generators.electronIons(spec));

      // Find the best-scoring hit for this spectrum
      double best_score = 0;
      AASequence best_seq;
      int best_isotope_error = 0;
      uint16_t best_charge = 0;
      float best_mean_error = 0;

      // Reused across this spectrum's candidates — same rationale as the main
      // scoring loop: a fresh PeakSpectrum per candidate churns its DataArrays
      // (ion names / charges filled by add_metainfo) on the heap.
      PeakSpectrum theo;

      for (const auto& sms : top_sms.hits_)
      {
        AASequence seq = fragment_index.reconstructModifiedSequence(
            fragment_index.getPeptides()[sms.peptide_idx_], db);
        // Clear peaks + data arrays before refilling; getSpectrum appends to
        // whatever is there. Its output is already sorted with add_metainfo=true.
        theo.clear(true);
        tsg.getSpectrum(theo, seq, 1, scoringMaxCharge_(sms.precursor_charge_));

        // The calibration PSMs are selected with the search's own score (upstream 825c33bb). The mass-accuracy score
        // favours PSMs with small fragment errors, so it narrows the fragment tolerance estimated from them below
        // (2-10% narrower than with HyperScore on 12 of 14 ppm runs); selecting by HyperScore instead did not
        // change the yield measurably.
        HyperScore::PSMDetail detail;
        double score = mass_accuracy_score_
          ? HyperScore::computeMassAccuracy(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, spec, theo, mass_error_sd_ppm_, detail)
          : HyperScore::computeWithDetail(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, spec, theo, detail);

        if (score > best_score)
        {
          best_score = score;
          best_seq = std::move(seq);
          best_isotope_error = sms.isotope_error_;
          best_charge = sms.precursor_charge_;
          best_mean_error = static_cast<float>(detail.mean_error);
        }
      }

      if (best_score == 0 || best_seq.empty()) continue;

      // Skip PSMs matched at a non-zero isotope offset. Their precursor m/z carries
      // an extra source of uncertainty (the +1/+2 isotope peak can be ambiguous in
      // the MS1 precursor picking), which inflates the variance of the calibration
      // quantile estimate. iso_err=0 PSMs are the gold-standard "true monoisotopic
      // peak picked" subset — the right population for estimating instrument bias.
      if (best_isotope_error != 0) continue;

      // Compute precursor error (signed), isotope-corrected. FragmentIndex searches
      // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs +
      // isotope_error * C13C12; the observed-to-monoiso m/z correction is
      //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
      // Matches the sign used by postProcessHits_'s PRECURSOR_ERROR_PPM annotation.
      double exp_mz = spec.getPrecursors()[0].getMZ();
      double theo_mz = best_seq.getMZ(best_charge);
      double corrected_exp_mz = exp_mz + static_cast<double>(best_isotope_error)
                                          * Constants::C13C12_MASSDIFF_U / best_charge;
      double prec_err = (precursor_mass_tolerance_unit_ == "ppm")
                          ? Math::getPPM(corrected_exp_mz, theo_mz)
                          : (corrected_exp_mz - theo_mz);

#pragma omp critical (prose_calibration_hits)
      cal_hits.push_back({best_score, prec_err, static_cast<double>(best_mean_error)});
      }
      catch (...)
      {
#pragma omp critical (ProSEAlgorithm_calibration_error)
        {
          if (!calibration_error) calibration_error = std::current_exception();
        }
        calibration_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (calibration_error) std::rethrow_exception(calibration_error);


    // Filter to high-confidence PSMs: keep only top 50% by score (robust against
    // random matches inflating the error distribution tails). Only crop if the
    // top 50% still meets the minimum PSM threshold; otherwise leave all hits
    // and let the calibration_min_psms_ check below disable calibration.
    Size cal_hits_total = cal_hits.size();
    std::sort(cal_hits.begin(), cal_hits.end(),
              [](const CalHit& a, const CalHit& b) { return a.score > b.score; });
    Size keep = cal_hits.size() / 2;
    if (keep >= calibration_min_psms_ && keep < cal_hits.size())
    {
      cal_hits.resize(keep);
    }

    // Extract error vectors from filtered hits.
    // Additionally discard PSMs whose signed precursor error falls outside the band the
    // FragmentIndex actually delivers. The mass-space window is [observed - lower, observed + upper];
    // rewriting in terms of the signed error e = observed - theoretical:
    //   theoretical in [observed - lower, observed + upper]
    //     <=>  e in [-upper, +lower]
    // So legitimate hits from FragmentIndex have e in [-upper, +lower], NOT [-lower, +upper].
    // The bounds flip for asymmetric windows; matching on the wrong band silently drops
    // real matches when the user has compensated an instrument bias.
    vector<double> precursor_errors;
    vector<double> fragment_errors_abs;
    for (const auto& h : cal_hits)
    {
      if (h.prec_error < -precursor_mass_tolerance_upper_ ||
          h.prec_error > precursor_mass_tolerance_lower_) continue; // wrong match
      precursor_errors.push_back(h.prec_error);
      if (h.frag_error > 0) fragment_errors_abs.push_back(h.frag_error);
    }

    if (precursor_errors.size() < calibration_min_psms_)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: insufficient PSMs (" << precursor_errors.size()
                      << " < " << calibration_min_psms_
                      << "), using configured tolerances." << std::endl;
      return result; // success=false
    }

    // --- Compute calibrated tolerances ---
    const double min_tolerance = 1e-6; // avoid non-positive tolerances

    // Signed precursor errors: 0.5% and 99.5% empirical quantiles give a distribution-free
    // 99% window (same approach as NuXL autotune, OpenNuXL.cpp:4939-4941). Asymmetric by
    // construction, no Gaussian assumption required — heavy tails and biased distributions
    // are handled correctly. We previously used residual_median + 3*residual_MAD which
    // captures only ~94% coverage on a Gaussian and less on heavy-tailed data; the
    // quantile method consistently produces more realistic windows.
    std::vector<double> sorted_errs = precursor_errors;
    std::sort(sorted_errs.begin(), sorted_errs.end());
    const size_t n = sorted_errs.size();
    // Bias estimate — median of the signed errors. Kept as a diagnostic field; the
    // quantile bounds below don't use it directly (they already encode the bias via
    // the asymmetry of [lo, hi]), but callers log it and tests assert on it.
    result.precursor_shift = Math::median(sorted_errs.begin(), sorted_errs.end(), /*sorted=*/true);
    const double lo = sorted_errs[static_cast<size_t>(n * 0.005)];                   // ~most negative
    const double hi = sorted_errs[std::min(n - 1, static_cast<size_t>(n * 0.995))];  // ~most positive
    // Convention (see lines 2193-2197): signed error e = observed - theoretical lies in
    // [-cal_upper, +cal_lower]. Matching the observed [lo, hi] against that gives:
    //   -cal_upper = lo  →  cal_upper = -lo
    //   +cal_lower = hi  →  cal_lower =  hi
    const double cal_lower_raw = std::max(min_tolerance, hi);
    const double cal_upper_raw = std::max(min_tolerance, -lo);
    // "Extreme" guards the degenerate case where the signed-error distribution
    // has essentially zero spread — e.g. a pure uniform shift fixture or a
    // single-peptide corner case — and the quantile bounds can't usefully inform
    // a calibrated window. Precursor calibration is then discarded; fragment
    // calibration still applies. A distribution strictly on one side of zero is
    // NOT extreme by this definition: we emit a legal half-line window
    // (cal_upper or cal_lower clamped to min_tolerance).
    //
    // Unit-aware threshold: 1 ppm is the right floor for ppm tolerances —
    // realistic proteomics-scale fixtures never trip it. In Da mode, 1.0 Da
    // would flag every realistic calibration (spreads are <0.05 Da), so fall
    // back to min_tolerance.
    const double extreme_bias_threshold =
        (precursor_mass_tolerance_unit_ == "ppm") ? 1.0 : min_tolerance;
    result.extreme_bias = (hi - lo) < extreme_bias_threshold;
    result.precursor_spread = std::max(cal_lower_raw, cal_upper_raw);  // diagnostic only
    if (!result.extreme_bias)
    {
      // Only tighten — cap against user-configured bounds.
      result.cal_lower = std::min(cal_lower_raw, precursor_mass_tolerance_lower_);
      result.cal_upper = std::min(cal_upper_raw, precursor_mass_tolerance_upper_);
    }
    // else: cal_lower/cal_upper stay at 0; writeback block skips the calibration result.

    // Fragment: 4 × 68th percentile (following NuXL / identipy convention).
    std::sort(fragment_errors_abs.begin(), fragment_errors_abs.end());
    if (!fragment_errors_abs.empty())
    {
      double frag_68 = fragment_errors_abs[static_cast<Size>((fragment_errors_abs.size() - 1) * 0.68)];
      result.fragment_tolerance = 4.0 * frag_68;
      if (result.fragment_tolerance < min_tolerance) result.fragment_tolerance = min_tolerance;
    }
    else
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    result.fragment_shift = 0.0;

    // Don't widen beyond configured fragment tolerance (only tighten).
    if (result.fragment_tolerance > fragment_mass_tolerance_)
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    result.success = true;

    OPENMS_LOG_INFO << "[ProSE] Calibration: " << precursor_errors.size() << " PSMs used (top "
                    << static_cast<int>(100.0 * precursor_errors.size() / cal_hits_total) << "% by score)" << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor signed: shift=" << std::fixed << std::setprecision(2) << result.precursor_shift
                    << " " << precursor_mass_tolerance_unit_ << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor spread: -> " << result.precursor_spread
                    << " " << precursor_mass_tolerance_unit_
                    << (result.extreme_bias ? " (extreme bias, discarded)" : "") << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Fragment tolerance:  " << fragment_mass_tolerance_
                    << " -> " << result.fragment_tolerance << " " << fragment_mass_tolerance_unit_ << std::endl;

    return result;
  }

  // =====================================================================
  // Serialize the end-of-search report to a YAML string. Hand-rolled (no YAML
  // library), so neither the TOPP tool nor the library needs an extra dependency
  // for this. Every string scalar is double-quoted (q()) so values containing ':'
  // (e.g. Windows paths) or other metacharacters cannot break the structure, and
  // non-finite numbers are emitted as null (num()).
  // =====================================================================
  std::string ProSEAlgorithm::renderRunSummaryYaml(
      const MultiFileSearchResult& mfres,
      const std::vector<std::pair<std::string, std::vector<std::string>>>& manifest,
      Size files_failed,
      Size files_total)
  {
    // Double-quoted YAML scalar: escape the two structural chars, the common
    // whitespace escapes, and every other C0 control + DEL as \xNN. A literal
    // control character is ill-formed in a double-quoted scalar and would corrupt
    // or invalidate the document — reachable via arbitrary input file paths.
    auto q = [](const std::string& s) {
      static const char* const hex = "0123456789ABCDEF";
      std::string out;
      out.reserve(s.size() + 2);
      out += '"';
      for (char ch : s)
      {
        const unsigned char c = static_cast<unsigned char>(ch);
        if (c == '\\' || c == '"') { out += '\\'; out += static_cast<char>(c); }
        else if (c == '\n') { out += "\\n"; }
        else if (c == '\t') { out += "\\t"; }
        else if (c == '\r') { out += "\\r"; }
        else if (c < 0x20 || c == 0x7F) { out += "\\x"; out += hex[c >> 4]; out += hex[c & 0x0F]; }
        else { out += static_cast<char>(c); }
      }
      out += '"';
      return out;
    };
    // Locale-independent number; non-finite -> null. YAML 1.1 reads an exponent
    // without a mantissa dot ("8e-05") as a string, so ensure a '.' before 'e'.
    auto num = [](double v) -> std::string {
      if (!std::isfinite(v)) { return "null"; }
      std::ostringstream o;
      o.imbue(std::locale::classic());
      o << v;
      std::string s = o.str();
      const auto e = s.find_first_of("eE");
      if (e != std::string::npos && s.find('.') == std::string::npos) { s.insert(e, ".0"); }
      return s;
    };
    auto bol = [](bool b) -> std::string { return b ? "true" : "false"; };

    std::ostringstream y;
    const SharedSearchStats& sh = mfres.shared;

    y << "shared:\n";
    y << "  database_file: " << q(sh.database_file) << "\n";
    y << "  enzyme: " << q(sh.enzyme) << "\n";
    y << "  precursor_tol_lower: " << num(sh.precursor_tol_lower) << "\n";
    y << "  precursor_tol_upper: " << num(sh.precursor_tol_upper) << "\n";
    y << "  precursor_tol_unit: " << q(sh.precursor_tol_unit) << "\n";
    y << "  fragment_tol: " << num(sh.fragment_tol) << "\n";
    y << "  fragment_tol_unit: " << q(sh.fragment_tol_unit) << "\n";
    y << "  min_charge: " << sh.min_charge << "\n";
    y << "  max_charge: " << sh.max_charge << "\n";
    y << "  missed_cleavages: " << sh.missed_cleavages << "\n";
    auto str_list = [&](const char* key, const std::vector<std::string>& items) {
      y << "  " << key << ":";
      if (items.empty()) { y << " []\n"; return; }
      y << "\n";
      for (const auto& it : items) { y << "    - " << q(it) << "\n"; }
    };
    str_list("fixed_mods", std::vector<std::string>(sh.fixed_mods.begin(), sh.fixed_mods.end()));
    str_list("variable_mods", std::vector<std::string>(sh.variable_mods.begin(), sh.variable_mods.end()));
    str_list("ion_series", std::vector<std::string>(sh.ion_series.begin(), sh.ion_series.end()));
    y << "  open_search: " << bol(sh.open_search) << "\n";
    y << "  calibration_enabled: " << bol(sh.calibration_enabled) << "\n";
    y << "  snes_mode: " << bol(sh.snes_mode) << "\n";
    y << "  chunked: " << bol(sh.chunked) << "\n";
    y << "  decoy_mode: " << q(sh.decoy_mode) << "\n";
    y << "  psm_fdr_threshold: " << num(sh.psm_fdr_threshold) << "\n";
    y << "  protein_fdr_threshold: " << num(sh.protein_fdr_threshold) << "\n";
    y << "  db_target_proteins: " << sh.db_target_proteins << "\n";
    y << "  db_decoy_proteins: " << sh.db_decoy_proteins << "\n";
    y << "  indexed_peptides: " << sh.indexed_peptides << "\n";
    y << "  indexed_fragments: " << sh.indexed_fragments << "\n";
    y << "  seconds_index_build: " << num(sh.seconds_index_build) << "\n";
    y << "  seconds_total: " << num(sh.seconds_total) << "\n";

    y << "per_file:";
    if (mfres.per_file.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& pf : mfres.per_file)
      {
        const RunStatistics& st = pf.stats;
        // The first key of each list item carries the "- " marker; the rest indent to 4.
        y << "  - input_file: " << q(st.input_file) << "\n";
        y << "    exit_code: " << static_cast<int>(pf.exit_code) << "\n";
        y << "    ms2_spectra: " << st.ms2_spectra << "\n";
        y << "    matched_spectra: " << st.matched_spectra << "\n";
        y << "    target_psms: " << st.target_psms << "\n";
        y << "    decoy_psms: " << st.decoy_psms << "\n";
        y << "    fdr_applied: " << bol(st.fdr_applied) << "\n";
        y << "    achieved_psm_fdr: " << num(st.achieved_psm_fdr) << "\n";
        y << "    unique_peptides: " << st.unique_peptides << "\n";
        y << "    unique_proteins: " << st.unique_proteins << "\n";
        y << "    hyperscore:\n";
        y << "      valid: " << bol(st.score_stats_valid) << "\n";
        y << "      min: " << num(st.hyperscore_min) << "\n";
        y << "      median: " << num(st.hyperscore_median) << "\n";
        y << "      max: " << num(st.hyperscore_max) << "\n";
        y << "    charge_histogram:";
        if (st.charge_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [z, c] : st.charge_histogram) { y << "      " << q(std::to_string(z)) << ": " << c << "\n"; } }
        y << "    missed_cleavage_histogram:";
        if (st.missed_cleavage_histogram.empty()) { y << " {}\n"; }
        else { y << "\n"; for (const auto& [m, c] : st.missed_cleavage_histogram) { y << "      " << q(std::to_string(m)) << ": " << c << "\n"; } }
        y << "    precursor_error:\n";
        y << "      valid: " << bol(st.prec_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.prec_err_median) << "\n";
        y << "      mad_ppm: " << num(st.prec_err_mad) << "\n";
        y << "      recommended_ppm: " << num(st.prec_err_recommended) << "\n";
        y << "    fragment_error:\n";
        y << "      valid: " << bol(st.frag_tol_valid) << "\n";
        y << "      median_ppm: " << num(st.frag_err_median) << "\n";
        y << "      mad_ppm: " << num(st.frag_err_mad) << "\n";
        y << "      recommended_ppm: " << num(st.frag_err_recommended) << "\n";
        y << "    timing_seconds:\n";
        y << "      calibration: " << num(st.seconds_calibration) << "\n";
        y << "      search: " << num(st.seconds_search) << "\n";
        y << "      fdr: " << num(st.seconds_fdr) << "\n";
      }
    }

    y << "outputs:";
    if (manifest.empty()) { y << " []\n"; }
    else
    {
      y << "\n";
      for (const auto& [label, paths] : manifest)
      {
        y << "  - type: " << q(label) << "\n";
        y << "    paths:";
        if (paths.empty()) { y << " []\n"; }
        else { y << "\n"; for (const auto& p : paths) { y << "      - " << q(p) << "\n"; } }
      }
    }
    y << "files_failed: " << files_failed << "\n";
    y << "files_total: " << files_total << "\n";

    return y.str();
  }

  // =====================================================================
  // Render the modification-discovery section (open search) to an ostream.
  // =====================================================================
  void ProSEAlgorithm::renderModificationSummary(
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      std::ostream& os)
  {
    const auto& dm_stats = mod_analysis.delta_mass_stats;
    const auto& ptm_stats = mod_analysis.ptm_stats;

    os << "[ProSE] ---------------- Modification discovery ---------------\n";
    os << "[ProSE] Delta mass : " << dm_stats.modified_psms << " modified / "
       << dm_stats.total_psms << " PSMs ("
       << std::fixed << std::setprecision(1)
       << (dm_stats.total_psms > 0 ? (100.0 * dm_stats.modified_psms / dm_stats.total_psms) : 0.0)
       << "%), median Δ=" << std::setprecision(4) << dm_stats.median_delta_mass
       << " Da, " << dm_stats.entries.size() << " bins\n";
    os << "[ProSE] PTMs      : " << ptm_stats.total_modified_psms << " PSMs with known PTMs, "
       << ptm_stats.unknown_modification_psms << " unknown, "
       << ptm_stats.num_unique_modifications << " unique PTMs\n";

    if (!ptm_stats.entries.empty())
    {
      os << "[ProSE]   Top PTMs (rank | name | count | % | mass Da):\n";
      size_t rank = 1;
      for (const auto& ptm : ptm_stats.entries)
      {
        if (rank > 15) { break; }
        std::string name = ptm.name;
        if (name.size() > 30) { name = StringUtils::substr(name, 0, 27) + "..."; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | "
           << std::setw(31) << std::left << name << std::right << " | "
           << std::setw(6) << ptm.count << " | "
           << std::setw(5) << std::fixed << std::setprecision(1) << ptm.percentage << " | "
           << std::setw(9) << std::fixed << std::setprecision(4) << ptm.theoretical_mass << "\n";
      }
    }

    std::vector<OpenSearchModificationAnalysis::DeltaMassEntry> unknown_dm;
    for (const auto& entry : dm_stats.entries)
    {
      if (!entry.is_known_modification && entry.count >= 5) { unknown_dm.push_back(entry); }
    }
    if (!unknown_dm.empty())
    {
      os << "[ProSE]   Top unknown delta masses (potential novel PTMs):\n";
      size_t rank = 1;
      for (const auto& dm : unknown_dm)
      {
        if (rank > 10) { break; }
        os << "[ProSE]   " << std::setw(2) << rank++ << " | Δ="
           << std::setw(11) << std::fixed << std::setprecision(4) << dm.delta_mass << " Da | "
           << std::setw(6) << dm.count << " PSMs | "
           << dm.unique_peptides << " peptides\n";
      }
    }
  }

  // =====================================================================
  // Render a human-readable single-run summary block.
  // =====================================================================
  void ProSEAlgorithm::renderRunSummary(
      const RunStatistics& s,
      const SharedSearchStats& sh,
      const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
      bool is_open_search,
      std::ostream& os)
  {
    auto join = [](const std::vector<std::string>& v) -> std::string {
      if (v.empty()) { return "(none)"; }
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? ", " : "") + v[i]; }
      return out;
    };
    auto join_int = [](const std::vector<std::string>& v) -> std::string {
      std::string out;
      for (size_t i = 0; i < v.size(); ++i) { out += (i ? "," : "") + v[i]; }
      return out.empty() ? std::string("(none)") : out;
    };

    if (!s.input_file.empty())
    {
      os << "[ProSE] Input        : " << s.input_file << "  (" << s.ms2_spectra << " MS2 spectra)\n";
    }
    // -- Configuration recap --
    os << "[ProSE] Config       : " << sh.enzyme << ", " << sh.missed_cleavages << " MC | prec [-"
       << sh.precursor_tol_lower << ", +" << sh.precursor_tol_upper << "] " << sh.precursor_tol_unit
       << " | frag " << sh.fragment_tol << " " << sh.fragment_tol_unit
       << " | z " << sh.min_charge << "-" << sh.max_charge << "\n";
    os << "[ProSE]                fixed: " << join(sh.fixed_mods)
       << " | variable: " << join(sh.variable_mods) << "\n";
    os << "[ProSE]                ions: " << join_int(sh.ion_series)
       << " | calibration: " << (sh.calibration_enabled ? "on" : "off")
       << " | mode: " << (sh.open_search ? "open" : "closed")
       << (sh.snes_mode ? " | SNES" : "")
       << (sh.chunked ? " | chunked" : "") << "\n";

    // -- Database / fragment index --
    os << "[ProSE] Database     : " << (sh.database_file.empty() ? std::string("(in-memory)") : std::string(sh.database_file))
       << "  (" << sh.db_target_proteins << " target";
    if (sh.db_decoy_proteins > 0) { os << " + " << sh.db_decoy_proteins << " decoy"; }
    os << " proteins, decoys: " << (sh.decoy_mode.empty() ? std::string("n/a") : std::string(sh.decoy_mode)) << ")\n";
    os << "[ProSE] Fragment idx : " << sh.indexed_peptides << " peptides / "
       << sh.indexed_fragments << " fragments  (build "
       << std::fixed << std::setprecision(1) << sh.seconds_index_build << " s)\n";

    // -- Results --
    os << "[ProSE] --------------------------- Results ---------------------------\n";
    os << "[ProSE] Matched      : " << s.matched_spectra << " / " << s.ms2_spectra << " spectra";
    if (s.ms2_spectra > 0)
    {
      os << "   (ID rate " << std::fixed << std::setprecision(1)
         << (100.0 * s.matched_spectra / s.ms2_spectra) << "%)";
    }
    os << "\n";
    os << "[ProSE] Peptides/prot: " << s.unique_peptides << " unique peptides | "
       << s.unique_proteins << " unique proteins\n";

    // FDR / target-decoy. Decide target-only strictly from the database, not from
    // a zero PSM count (a decoy DB with no decoy PSMs is still FDR-capable).
    if (sh.decoy_mode == "none (target-only)")
    {
      os << "[ProSE] FDR          : n/a (target-only database)\n";
    }
    else
    {
      os << "[ProSE] Target/decoy : " << s.target_psms << " target / " << s.decoy_psms << " decoy PSMs";
      if (s.fdr_applied)
      {
        os << "  | PSM FDR <= " << std::fixed << std::setprecision(1) << (sh.psm_fdr_threshold * 100.0) << "%";
        if (s.achieved_psm_fdr >= 0.0)
        {
          os << ", achieved " << std::setprecision(2) << (s.achieved_psm_fdr * 100.0) << "%";
        }
        else
        {
          os << " (0 PSMs retained)";
        }
      }
      else
      {
        os << "  | PSM FDR filtering: off";
      }
      os << "\n";
    }

    // Charge distribution
    if (!s.charge_histogram.empty())
    {
      os << "[ProSE] Charges      :";
      for (const auto& [z, c] : s.charge_histogram) { os << " " << z << ":" << c; }
      os << "\n";
    }
    // Missed cleavages
    if (!s.missed_cleavage_histogram.empty())
    {
      os << "[ProSE] Missed clv   :";
      for (const auto& [mc, c] : s.missed_cleavage_histogram) { os << " " << mc << ":" << c; }
      os << "\n";
    }
    // Score distribution
    if (s.score_stats_valid)
    {
      os << "[ProSE] HyperScore   : min " << std::fixed << std::setprecision(1) << s.hyperscore_min
         << "  median " << s.hyperscore_median << "  max " << s.hyperscore_max << "\n";
    }

    // -- Tolerance estimation --
    if (s.prec_tol_valid || s.frag_tol_valid)
    {
      os << "[ProSE] --------------------- Tolerance estimate ---------------------\n";
      if (s.prec_tol_valid)
      {
        os << "[ProSE] Precursor    : median " << std::fixed << std::setprecision(2) << s.prec_err_median
           << " ppm, MAD " << s.prec_err_mad << "  -> recommend " << static_cast<int>(s.prec_err_recommended) << " ppm\n";
      }
      if (s.frag_tol_valid)
      {
        os << "[ProSE] Fragment     : median " << std::fixed << std::setprecision(2) << s.frag_err_median
           << " ppm, MAD " << s.frag_err_mad << "  -> recommend " << static_cast<int>(s.frag_err_recommended) << " ppm\n";
      }
    }

    // -- Timing --
    os << "[ProSE] ------------------------- Timing -----------------------------\n";
    os << "[ProSE] " << std::fixed << std::setprecision(1);
    if (s.seconds_calibration > 0.0) { os << "calib " << s.seconds_calibration << "s | "; }
    if (s.seconds_search > 0.0)      { os << "search " << s.seconds_search << "s | "; }
    if (s.seconds_fdr > 0.0)         { os << "fdr " << s.seconds_fdr << "s | "; }
    os << "index(shared) " << sh.seconds_index_build << "s\n";

    // -- Modification discovery (open search only) --
    if (is_open_search)
    {
      renderModificationSummary(mod_analysis, os);
    }
  }

} // namespace OpenMS
