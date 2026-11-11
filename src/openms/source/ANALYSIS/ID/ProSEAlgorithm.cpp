// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include "ProSEAlgorithmImpl.h"

#include <OpenMS/ANALYSIS/ID/CoIsolationDetector.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
// Not guarded: header-only constants with no ONNX dependency, needed to register the same
// 'peptdeep:instrument' values in either build.
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <fstream>
#include <functional>
#include <future>
#include <limits>
#include <numeric>
#include <iomanip>
#include <locale>
#include <ostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  using namespace Internal::ProSEAlgorithmHelpers;

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

    defaults_.setValue("precursor:min_charge", 2,
                       "Minimum precursor charge to be considered. Spectra without a precursor charge are searched with every "
                       "charge of the range; spectra whose recorded charge lies outside of it are searched at that charge, or "
                       "not at all with precursor:recorded_charge_outside_range 'skip'.");
    defaults_.setValue("precursor:max_charge", 5, "Maximum precursor charge to be considered (see precursor:min_charge).");
    defaults_.setValue("precursor:recorded_charge_outside_range", "skip",
                       "A spectrum whose recorded precursor charge lies outside of precursor:min_charge..max_charge: 'skip' "
                       "(the default, as OpenMS develop since #10713) does not search it; 'search' (the behaviour before #10713) "
                       "searches it at its recorded charge, up to charge 20. Spectra without a recorded charge are searched with every charge of the range.",
                       {"advanced"});
    defaults_.setValidStrings("precursor:recorded_charge_outside_range", {"search", "skip"});

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
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'. A terminus and a residue carry one modification each: a variable terminal modification (e.g. 'Acetyl (Protein N-term)', 'Gln->pyro-Glu (N-term Q)') is not searched where a fixed one sits on that terminus (e.g. 'TMT6plex (N-term)'), and a variable residue modification (e.g. 'Glutathione (C)') is not searched where a fixed one sits on that residue (e.g. 'Carbamidomethyl (C)').");
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
      "(marker auto-detected, prefix or suffix) or generate them if none are found. Supplied decoys are reused only "
      "at 0.5 to 2 decoys per target protein: target-decoy FDR needs about one decoy per target, so a database with "
      "far fewer or more is refused (use 'generate'). "
      "'generate': always (re)generate decoys from the target proteins (any pre-existing "
      "decoys are removed first to avoid decoy-of-decoy entries). "
      "'ignore': search the target proteins only — any decoys present in the database are removed.");
    defaults_.setValidStrings("decoys", {"auto","generate","ignore"} );

    defaults_.setValue("decoy_prefix", "DECOY_", "Accession prefix prepended to decoy proteins when ProSE generates them (modes 'auto' without existing decoys, and 'generate'). Pre-existing decoys in the database are detected automatically (any common marker, used as prefix or suffix) and do not require this setting.", {"advanced"});
    defaults_.setValue("decoy_method", "reverse",
      "How ProSE generates decoys (decoys 'generate', or 'auto' with a database without decoys). 'reverse' (ProSEAlgorithm's "
      "default; the ProSE tool's default is 'keep_termini'): "
      "every enzymatic segment of the protein is reversed except its last residue (the last segment entirely), an "
      "initial Met stays first with peptide:clip_nterm_methionine. 'keep_termini': the construction of OpenMS "
      "DecoyDatabase (method reverse, Decoy:keepPeptideNTerm/keepPeptideCTerm true) with the search enzyme and "
      "Decoy:non_shuffle_pattern 'KRP': in every segment of a digest without missed cleavages, the first and the last "
      "residue and every K, R and P stay in place and the rest is reversed. For enzymes whose cleavage sites are fixed by "
      "these residues (Trypsin, Trypsin/P, Lys-C, Arg-C and the like) every target peptide then has a decoy with the same "
      "start, length, mass and cleavage sites: also every protein N-terminal and Met-clipped one, so that "
      "peptide:clip_nterm_methionine and protein N-terminal modifications add the same number of target and decoy "
      "candidates (except peptides that equal their decoy, e.g. made of short segments or of K, R and P only; and "
      "modifications specific to a peptide's first residue, e.g. pyro-Glu on an N-terminal Q, apply to the target and "
      "its decoy alike only where that residue is kept in place). Semi- or "
      "non-specific digestion: as 'reverse'. Other enzymes: a warning names the decoys whose cleavage sites differ from the "
      "target's. No effect on decoys taken from the database (decoys 'auto' with decoys in it). The results differ from "
      "'reverse'.", {"advanced"});
    defaults_.setValidStrings("decoy_method", {"reverse", "keep_termini"});
    defaults_.setValue("decoy_parity_report", "false",
      "Diagnostic: log the candidate counts of the fragment index per label (target, decoy) and subclass (internal, "
      "protein N-terminal with and without initial M, Met-clipped; and their protein N-terminally acetylated occurrences, "
      "none with a fixed N-terminal modification such as TMT, which leaves no room for a variable one), and a "
      "checksum of generated decoys. Output unchanged.",
      {"advanced"});
    defaults_.setValidStrings("decoy_parity_report", {"true", "false"});

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
      "Contexts and outcomes of the ion priors. Both keep the ion series a/b, x/y, c and z apart (c and z ions are scored "
      "for electron-activated spectra, see ions:by_activation). 'rich': ion series, precursor charge, fragment charge, relative position, "
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

    defaults_.setValue("peptdeep:enable", "false", "Add PeptDeep prediction-based rescoring features (ms2_cosine, ms2_spectral_angle, ms2_pearson, ms2_frac_pred_found, rt_abs_error) to every PSM. Uses the models shipped in OpenMS' 'share/OpenMS/models' unless 'peptdeep:ms2_model'/'peptdeep:rt_model' name others. Requires an OpenMS built with ONNX support.");
    defaults_.setValidStrings("peptdeep:enable", {"true", "false"});
    defaults_.setValue("peptdeep:ms2_model", "models/peptdeep_ms2_dynamic.onnx", "PeptDeep MS2 fragment-intensity ONNX model, used when 'peptdeep:enable' is true. A relative name is resolved against OpenMS' shared-data directory; an absolute path is used as given.");
    defaults_.setValue("peptdeep:rt_model", "models/peptdeep_rt_dynamic.onnx", "PeptDeep retention-time ONNX model. See 'peptdeep:ms2_model'.");
    defaults_.setValue("peptdeep:instrument", "QE", "Instrument class passed to the PeptDeep MS2 model. An instrument that is not in this list is best approximated by the closest one that is.");
    // peptdeep's own names, in peptdeep's own order, and registered whether or not this build has
    // ONNX so that an ini file means the same thing either way.
    std::vector<std::string> peptdeep_instruments;
    peptdeep_instruments.reserve(ML::ALPHAPEPTDEEP_INSTRUMENTS.size());
    for (const std::string_view name : ML::ALPHAPEPTDEEP_INSTRUMENTS) { peptdeep_instruments.emplace_back(name); }
    defaults_.setValidStrings("peptdeep:instrument", peptdeep_instruments);
    defaults_.setValue("peptdeep:nce", -1.0, "Normalised collision energy for the PeptDeep MS2 model. Negative selects it automatically from the collision energy recorded in the spectra, refined by scoring a small grid on confident PSMs.");
    defaults_.setValue("peptdeep:rt_model_type", "b_spline", "Model mapping predicted onto observed retention time.", {"advanced"});
    defaults_.setValidStrings("peptdeep:rt_model_type", {"b_spline", "lowess", "linear"});
    defaults_.setValue("peptdeep:batch_size", 64, "Peptides per ONNX inference call. The calls run side by side, one per thread: small ones keep the threads evenly busy and the memory low (on iPRG2015 with 16 threads, 64 took 40% less time than 500 and 1 GB less memory).", {"advanced"});
    defaults_.setMinInt("peptdeep:batch_size", 1);
    defaults_.setSectionDescription("peptdeep", "PeptDeep (ONNX) prediction-based rescoring features");

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
                       "With low_memory, the repeated peptidoforms among the candidates of a spectrum are dropped instead. "
                       "Set false for occurrence-based legacy candidates.",
                       {"advanced"});
    defaults_.setValidStrings("peptide:deduplicate", {"true", "false"});
    defaults_.setValue("peptide:protein_mapping", "index",
                       "How the hits are mapped to their proteins (peptide evidences, target_decoy, protein_references, protein hits). "
                       "'index': from the digest of the fragment index, which holds every protein occurrence of a candidate (with a "
                       "protein-terminal modification: from a pass over the database for the hits' sequences); the result is "
                       "that of PeptideIndexing, which runs instead wherever the index cannot reproduce it exactly (e.g. low_memory, a specificity "
                       "other than full, a chunked database, symbols other than letters or long stretches "
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
    defaults_.setValue("peptide:motif", "", ""); // description: FragmentIndex::restrictSearchParameters()
    defaults_.setSectionDescription("peptide", "Peptide Options");

    // forwarded to FragmentIndex
    defaults_.setValue("low_memory", "false",
      "Non-specific searches (peptide:enzyme_specificity=none; ignored otherwise): index the peptide masses, as the "
      "prefixes of one mother peptide per protein position, instead of the fragments of every peptidoform, and "
      "generate the fragments of the peptidoforms in the precursor window at query time. Finds the same candidates "
      "with a fraction of the memory (e.g. for immunopeptidomics against a human database), but takes more time per "
      "spectrum, the more the wider the precursor window.");
    defaults_.setValidStrings("low_memory", {"true", "false"});

    defaults_.setValue("report:top_hits", 1,
                       "Maximum number of top scoring hits per spectrum that are reported. With 'rescore', only the best "
                       "hit per spectrum is reported after the rescoring, and the hits that Percolator rescores are set by "
                       "'rescoring:top_hits' (report:top_hits then applies to a run that is not rescored and to -out_pin). "
                       "FDR:PSM > 0 and the protein FDR (FDR:protein > 0) also keep the best hit of each spectrum only, so "
                       "values above 1 apply with 'rescore' false, FDR:PSM 0 and FDR:protein 0.");
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

    defaults_.setValue("rescore", "true",
                       "Rescore the PSMs of every search run (input file) in-process with Percolator (OpenMS' built-in "
                       "implementation, as PercolatorAdapter uses by default) before the FDR filters: Percolator q-values become "
                       "the main score (the SVM score and PEP are kept as meta values), Percolator rescores the top "
                       "'rescoring:top_hits' hits of each spectrum and only its best one per spectrum is kept, "
                       "and FDR:PSM and FDR:protein are applied to the rescored PSMs. Needs decoys (see 'decoys') and at least "
                       "100 PSMs: a run without them, or for which Percolator fails, is not rescored and keeps its search engine "
                       "scores, which the FDR filters then use (the ProSE tool warns). If that happens to some files of a "
                       "multi-file search, the files rescored on their own keep their results, and the others get the "
                       "scores of a model trained on the PSMs of those (on the scale of its training PSMs, with their PEP "
                       "at the same score; q-values by target-decoy competition of each file's PSMs on that score); if "
                       "Percolator can rescore no file on its own, all files are rescored together with one model. Only "
                       "if that fails too, the files keep their search engine scores, and the pooled evaluation of the "
                       "files (the ProSE tool's -out_merged) scores the proteins by ProSE's score 'ln(hyperscore)' of all "
                       "PSMs instead of Percolator's PEP. See percolator:max_train_psms for the training PSMs. Set to "
                       "'false' to report the search engine scores. Keep it on for intermediates (FDR:protein 0) that "
                       "ProteomicsLFQ finishes: it takes the PEP (meta value MS:1001493) as main score and stops for an input "
                       "without PEPs.");
    defaults_.setValidStrings("rescore", {"true", "false"});
    defaults_.setValue("percolator:max_train_psms", 500000,
                       "With 'rescore': the most PSMs Percolator trains a model on (0: all). Every PSM is scored, also "
                       "beyond the cap. For a run with more PSMs, Percolator samples the training PSMs by spectrum (see "
                       "'rescoring:subset_mode'); for the model of the files of a multi-file search too small for a model of "
                       "their own (trained on the files rescored on their own), each of these files contributes at most "
                       "this many divided by the number of input files (reproducible random spectra), so that none "
                       "dominates. The run log says when the training PSMs were sampled. Searches below "
                       "the cap are not affected. Until all files of such a search are searched, ProSE keeps their PSMs "
                       "before the rescoring for that model: all PSMs of the files too small, and of the others only their "
                       "share (with 0, all PSMs of every file, also with the ProSE tool's -out_pin).", {"advanced"});
    defaults_.setMinInt("percolator:max_train_psms", 0);
    defaults_.setSectionDescription("percolator", "Percolator rescoring ('rescore')");
    defaults_.setValue("rescore_solver", "newton",
                       "Solver of Percolator's SVMs with 'rescore' (Percolator's parameter 'solver'): 'newton' (default; "
                       "explicit-Hessian Newton method with exact line search and warm starts: the same SVM optimum for a "
                       "given training set as TRON, within the solver tolerance, much faster on large runs; as the earlier "
                       "iterations of the cross-validation stop at a looser tolerance, they may choose slightly different "
                       "training sets, so the models can differ from TRON's), 'tron' (Percolator's trust region Newton "
                       "method, the solver of the percolator executable and of earlier ProSE versions) or "
                       "'lda' (Fisher's linear discriminant instead of the SVM: fastest, an approximation).", {"advanced"});
    defaults_.setValidStrings("rescore_solver", {"newton", "tron", "lda"});
    defaults_.setValue("rescoring:top_hits", 2,
                       "With 'rescore': the top scoring hits per spectrum (by the search engine score) that Percolator "
                       "rescores: its training and scoring rows, of which its target-decoy competition keeps the best per "
                       "spectrum, which is reported. 2 (default) gives Percolator the runner-up of each spectrum as well; "
                       "more hits cost time and memory. 0 = 'report:top_hits' (as earlier "
                       "ProSE versions, whose Percolator rescored the reported hits). The search keeps these hits for "
                       "Percolator (with -out_pin also the reported ones): a run that Percolator does not rescore (fewer "
                       "than 100 PSMs, no decoy PSM, a failure) reports at most rescoring:top_hits (if above 0) of its "
                       "report:top_hits hits per spectrum; with decoys 'ignore' (no rescoring) report:top_hits.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:top_hits", 0);
    defaults_.setValue("rescoring:subset_mode", "out_of_fold",
                       "With a training cap ('percolator:max_train_psms' > 0) on a run above it (Percolator's parameter 'subset_mode'): 'out_of_fold' "
                       "(default) assigns every spectrum to one of Percolator's three cross-validation folds first, draws "
                       "the training sample (whole spectra, at random from the entire run, seeded), trains each fold's "
                       "model on the sampled spectra of the other folds only, and scores every PSM with the model of its "
                       "own fold, which never saw its spectrum. 'out_of_fold_avg': the same, but the PSMs of the spectra "
                       "outside the sample (which no model saw) are scored by the average of the fold models. "
                       "'in_sample' = the percolator executable's --subset-max-train: cross-validation on the sample and "
                       "all PSMs scored by the average of the fold models (the sampled PSMs by two models that trained on "
                       "them). 'learning_curve' (also without a cap): out-of-fold training on nested samples of "
                       "'rescoring:curve_start', twice, four times ... as many PSMs, up to 'percolator:max_train_psms' (0: "
                       "1,000,000) or the whole run but the validation spectra (the end of the random order, at most "
                       "250,000 PSMs or 20 %); each fold's model stops at the first doubling that adds less than "
                       "'rescoring:curve_min_gain' or 2 standard errors of targets at 1 % among the validation spectra of the "
                       "other folds (no label of its own fold is used), and is not trained further; every PSM is scored by its fold's model. The sizes and the curve are "
                       "recorded as Percolator:subset_curve.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:subset_mode", {"out_of_fold", "out_of_fold_avg", "in_sample", "learning_curve"});
    defaults_.setValue("rescoring:curve_start", 100000, "'rescoring:subset_mode' learning_curve: PSMs of the first sample.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:curve_start", 1);
    defaults_.setValue("rescoring:curve_min_gain", 0.002,
                       "'rescoring:subset_mode' learning_curve: a fold stops at the first doubling that adds less than "
                       "this fraction of validation targets.", {"advanced"});
    defaults_.setMinFloat("rescoring:curve_min_gain", 0.0);
    defaults_.setValue("rescoring:iterations", 10,
                       "Iterations of Percolator's cross-validation training (the percolator executable's --maxiter). "
                       "Fewer are faster; the results differ.",
                       {"advanced"});
    defaults_.setMinInt("rescoring:iterations", 1);
    defaults_.setValue("rescoring:tron_epsilon", 1e-7,
                       "Stopping tolerance of the SVM solver (TRON) of Percolator's training (Percolator: 1e-7). Larger "
                       "values (e.g. 1e-3) train much faster on large runs; the results differ.",
                       {"advanced"});
    defaults_.setMinFloat("rescoring:tron_epsilon", 1e-15);
    defaults_.setValue("rescoring:warm_start", "false",
                       "Each SVM training of an iteration starts from the weights of the previous iteration instead of "
                       "0 (faster; the results differ).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:warm_start", {"true", "false"});
    defaults_.setValue("rescoring:single_cneg", "false",
                       "Percolator trains with one Cneg (3 x the target/decoy ratio x Cpos) instead of selecting among "
                       "three by cross-validation (faster; the results differ).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:single_cneg", {"true", "false"});
    defaults_.setValue("rescoring:use_pi0", "true",
                       "Percolator's q-values and PEPs with its pi0 estimate (default) or with pi0 = 1, as the percolator "
                       "executable computes them with target-decoy competition (-Y). 'false' gives slightly fewer PSMs and "
                       "a lower measured false discovery proportion; the results differ. 'auto' decides as the executable (OpenMS "
                       "develop's in-process default since #10721): with ProSE's target-decoy competition that is 'false'.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:use_pi0", {"true", "false", "auto"});
    defaults_.setValue("rescoring:calc_mass_feature", "false",
                       "Percolator's training uses the CalcMass column (the theoretical m/z) as a feature. 'false' (the default): "
                       "as the percolator executable, which reads it as metadata (OpenMS develop since #10717); 'true': the "
                       "behaviour before #10717. The results differ.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:calc_mass_feature", {"true", "false"});
    defaults_.setValue("rescoring:fixed_enzyme_termini", "false",
                       "A protein terminus counts as an enzymatic site for the enzN/enzC features of the rescoring and of "
                       "'-out_pin' (Percolator's convention). Default 'false' keeps the former features, in which a peptide "
                       "at a protein N-terminus has enzN = 0 and one at a protein C-terminus that does not end in a "
                       "cleavage residue has enzC = 0; the results differ. The peptide after an initial Met (the Met-clipped protein "
                       "N-terminal peptide of peptide:clip_nterm_methionine) counts as N-terminal too. Experimental: tested "
                       "only with decoy_method 'keep_termini'; without the Met-clip rule, with the decoys of a database of "
                       "reversed proteins (which rarely start with M) enzN = 0 was left almost only to Met-clipped peptides "
                       "of targets and the FDR was no longer controlled (entrapment FDP 13-24 % at 1 % q).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:fixed_enzyme_termini", {"true", "false"});
    defaults_.setValue("rescoring:terminal_features", "false",
                       "Two more rescoring features: 'protein_nterm' (1 if the peptide occurs at a protein N-terminus) and "
                       "'met_clipped' (1 if it occurs only as the Met-clipped protein N-terminal peptide of "
                       "peptide:clip_nterm_methionine). They let the rescoring learn the different rate of correct "
                       "matches of these small peptide classes. Use them only with decoys that populate both classes like "
                       "the targets (decoy_method 'keep_termini'): with decoys that rarely start with M (reversed "
                       "proteins, reused with decoys 'auto') they would tell targets from decoys, so they are switched off "
                       "with a warning. The results differ.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:terminal_features", {"true", "false"});
    defaults_.setValue("rescoring:stratum_q_values", "false",
                       "Also q-values per stratum (meta values 'stratum' and 'stratum_q_value'): the rescored PSMs of "
                       "protein N-terminal peptides (incl. Met-clipped ones) and the others each get target-decoy "
                       "competition of their own ((decoys + 1) / targets, as mokapot's grouped confidence). For a list of "
                       "N-terminal peptides at an FDR, use the stratum q-value; the main score stays the pooled q-value. "
                       "Needs decoys that populate the N-terminal classes like the targets (decoy_method 'keep_termini'; with reused "
                       "decoys that rarely start with M it is switched off with a warning). With "
                       "rescore_method 'reset' the values are computed on RESET's scores (both decoy sets count as decoys "
                       "with reset:decoys 'two'), not as RESET's estimate.",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:stratum_q_values", {"true", "false"});
    defaults_.setValue("rescoring:peptide_level", "true",
                       "Also Percolator's peptide-level q-values (default; 'false' leaves them out), as the percolator executable reports peptides (without -U): "
                       "after the target-decoy competition of the spectra, the best PSM of each peptide (modified sequence) "
                       "and label competes, and q-values and PEPs are estimated on these peptides. Every rescored hit gets "
                       "the values of its peptide as meta values 'percolator_peptide_q_value' and 'percolator_peptide_pep'; "
                       "a list of peptides at an FDR is the target peptides with percolator_peptide_q_value at most that FDR "
                       "(the distinct peptides of the PSMs at an FDR are not FDR-controlled at the peptide level). The main "
                       "score stays the PSM-level q-value. Default Percolator rescoring only (not with rescore_method 'reset', "
                       "which reports RESET:peptide_q_value, nor with rescore_reuse models).",
                       {"advanced"});
    defaults_.setValidStrings("rescoring:peptide_level", {"true", "false"});
    defaults_.setSectionDescription("rescoring", "Options of the Percolator training of the rescoring (see 'rescore')");
    defaults_.setValue("rescore_method", "percolator",
                       "With 'rescore': 'percolator' = Percolator's PSM-level FDR (default); 'reset' = Percolator-RESET "
                       "(Freestone et al., J Proteome Res 2025): half of the decoy peptides (stems = unmodified sequences, by a "
                       "random coin per stem) train the model, the other half estimates the FDR, which is controlled at the "
                       "peptide (stem) level: every rescored hit gets the q-value of its stem as meta value "
                       "'RESET:peptide_q_value'. The main score (q-value) stays PSM-level, but under RESET it is a heuristic "
                       "(not FDR-controlled): 2 (estimating decoy PSMs + 1) / target PSMs. Experimental.", {"advanced"});
    defaults_.setValidStrings("rescore_method", {"percolator", "reset"});
    defaults_.setValue("reset:salt", "",
                       "rescore_method 'reset': salt of the decoy coins (unsigned 64-bit integer); empty = random per run. "
                       "The salt is recorded in the search parameters ('RESET:salt'). Never pick the best of several salts.",
                       {"advanced"});
    defaults_.setValue("reset:learner", "svm",
                       "rescore_method 'reset': 'svm' = linear SVM; 'mlp' = small MLP ensemble (experimental).", {"advanced"});
    defaults_.setValidStrings("reset:learner", {"svm", "mlp"});
    defaults_.setValue("reset:mlp_members", 8, "rescore_method 'reset', learner 'mlp': networks per ensemble.", {"advanced"});
    defaults_.setMinInt("reset:mlp_members", 1);
    defaults_.setValue("reset:mlp_hidden", 16, "rescore_method 'reset', learner 'mlp': units of each of the two hidden layers.", {"advanced"});
    defaults_.setMinInt("reset:mlp_hidden", 1);
    defaults_.setValue("reset:mlp_epochs", 30, "rescore_method 'reset', learner 'mlp': training epochs of the first round.", {"advanced"});
    defaults_.setMinInt("reset:mlp_epochs", 1);
    defaults_.setValue("reset:mlp_warm_epochs", 3,
                       "rescore_method 'reset', learner 'mlp': > 0 = the networks continue from the previous round for this many "
                       "epochs (warm starts); 0 = trained afresh (mlp_epochs) every round.", {"advanced"});
    defaults_.setMinInt("reset:mlp_warm_epochs", 0);
    defaults_.setValue("reset:mlp_crossfit", "true",
                       "rescore_method 'reset', learner 'mlp': score each peptide by the networks that did not train on it (3 folds "
                       "by sequence hash); 'false' = one final ensemble (memorizes its training decoys: fewer peptides).", {"advanced"});
    defaults_.setValidStrings("reset:mlp_crossfit", {"true", "false"});
    defaults_.setValue("reset:mlp_norm_min_decoys", 1000,
                       "rescore_method 'reset', learner 'mlp', cross-fitted: normalize the folds' scores (as Percolator) only if every "
                       "fold has at least this many training decoys (0 = always).", {"advanced"});
    defaults_.setMinInt("reset:mlp_norm_min_decoys", 0);
    defaults_.setValue("reset:representative", "model",
                       "rescore_method 'reset': the PSM that represents a peptide (stem): 'model' (default) = after each training "
                       "round, best score under the new model (+2 % peptides on DIA; the FDR then also assumes that null targets "
                       "and decoys are exchangeable over their whole PSM collections, including their sizes); 'search' = best "
                       "search engine score (the paper).", {"advanced"});
    defaults_.setValidStrings("reset:representative", {"search", "model"});
    defaults_.setValue("reset:iterations", 25,
                       "rescore_method 'reset': training rounds (the paper: 5; on DIA pseudo-spectra 5 rounds stop far from "
                       "convergence: 82 % of the peptides of 20 rounds on HeLa).", {"advanced"});
    defaults_.setMinInt("reset:iterations", 1);
    defaults_.setValue("reset:audit", "false",
                       "rescore_method 'reset': feature-symmetry audit before training (log).", {"advanced"});
    defaults_.setValidStrings("reset:audit", {"true", "false"});
    defaults_.setValue("reset:decoys", "single",
                       "rescore_method 'reset': 'single' = one decoy database (default); 'two' = the two-decoy mode of "
                       "Percolator-RESET: the database holds a second decoy set of the same construction (accessions starting "
                       "with reset:second_decoy_prefix, also recognized as decoys by the decoy marker), the spectra are "
                       "searched against targets and both decoy sets, and RESET estimates the FDR with the factors of two decoy "
                       "databases (fewer peptides lost to the decoy split, less variation between salts; the search is larger).",
                       {"advanced"});
    defaults_.setValidStrings("reset:decoys", {"single", "two"});
    defaults_.setValue("reset:second_decoy_prefix", "DECOY_X_",
                       "rescore_method 'reset', reset:decoys 'two': accession prefix of the proteins of the second decoy set.",
                       {"advanced"});
    defaults_.setSectionDescription("reset", "Percolator-RESET options (rescore_method 'reset')");
    defaults_.setValue("reset:stem_table", "",
                       "rescore_method 'reset': write the stems with their q-values (tab-separated) to this file, or, if it "
                       "ends with '/' or is a directory, to <dir>/<spectra file base name>.reset_stems.tsv per run.",
                       {"advanced"});
    defaults_.setValue("rescore_reuse:mode", "per_file",
                       "Which PSMs train the Percolator model of a run ('rescore'). 'per_file': the run itself, with Percolator's "
                       "three-fold cross-validation. The other modes cross-fit: the spectra are split into two halves at random "
                       "(by Percolator's seed), one model is trained per half on about 'rescore_reuse:max_train' PSMs, and every PSM is "
                       "scored by the model of the other half, so that no PSM is scored by a model that was trained on it (5 % of the "
                       "spectra of each half are not trained on: their decoys put the scores of their own half's model on one scale; they are scored by the other model). 'joint' (a "
                       "multi-file search): the training PSMs are drawn from all files, the same number from each, and the files are "
                       "rescored once all are searched (mokapot's joint model). 'static' (a multi-file search): they are drawn from the first "
                       "'rescore_reuse:static_files' files, and every later file is rescored with these models after its search, unless "
                       "the features of its PSMs differ from those of the training PSMs by more than 'rescore_reuse:max_shift': "
                       "then it is rescored as with 'per_file'. 'own': every file trains its two models on its own PSMs. In every mode the "
                       "q-values and PEPs of a file come from its own target-decoy competition (Percolator's, as with 'per_file'). 'joint' and "
                       "'static' need an unchunked closed search of several files and keep all files' PSMs in memory until they are "
                       "rescored; otherwise 'per_file' applies.",
                       {"advanced"});
    defaults_.setValidStrings("rescore_reuse:mode", {"per_file", "joint", "static", "own"});
    defaults_.setValue("rescore_reuse:max_train", 250000,
                       "Training PSMs per model of 'rescore_reuse:mode' (0: all PSMs of the half but the calibration spectra). The spectra are "
                       "drawn at random, all PSMs of a drawn spectrum together, so the number is met on average; the files share it "
                       "equally (a file with fewer PSMs gives all).",
                       {"advanced"});
    defaults_.setMinInt("rescore_reuse:max_train", 0);
    defaults_.setValue("rescore_reuse:static_files", 1, "With 'rescore_reuse:mode' static: the number of (first) input files that train the models.",
                       {"advanced"});
    defaults_.setMinInt("rescore_reuse:static_files", 1);
    defaults_.setValue("rescore_reuse:max_shift", 0.5,
                       "With 'rescore_reuse:mode' static: the largest domain shift for which a file is rescored with the static models. "
                       "The shift is max |2 AUC - 1| over the features the models use, AUC comparing the file's PSMs with the "
                       "training PSMs, targets and decoys alike (0: same distribution, 1: separated).",
                       {"advanced"});
    defaults_.setMinFloat("rescore_reuse:max_shift", 0.0);
    defaults_.setMaxFloat("rescore_reuse:max_shift", 1.0);
    defaults_.setSectionDescription("rescore_reuse", "Percolator model reuse across the files of a search ('rescore')");

    defaults_.setValue("FDR:PSM", 0.0,
                       "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR; 0 disables the filter). Target and decoy PSMs are filtered "
                       "alike by the threshold; decoys that pass are kept (no decoy-specific stripping here — decoys are removed only at "
                       "protein-FDR finalization). With 'rescore' (the default), the threshold applies to the Percolator q-values, the main "
                       "score of the rescored PSMs. Without rescoring ('rescore' false, or a run that Percolator did not rescore), 0 computes "
                       "no q-values: every PSM keeps its search engine score ('ln(hyperscore)') as main score, e.g. as input for Percolator "
                       "or a downstream FDR step. Above 0, the target-decoy q-value (see FDR:PSM_groups) replaces the main score (the search "
                       "engine score is kept as meta value 'ln(hyperscore)_score') and only the best hit of each spectrum is kept. If no PSM "
                       "matched a decoy, the q-value of the T PSMs is the conservative 1/T (target-decoy competition with +1), "
                       "so a lower threshold reports no PSM instead of lifting the filter. Needs decoys ('decoys' auto or "
                       "generate); with 'ignore', a threshold above 0 is skipped with a warning (the QPX run view of the ProSE tool's "
                       "-qpx_extended records per input whether the filter was applied: 'psm_fdr_applied'). Afterwards, the "
                       "proteins that no remaining PSM references are removed. The final PSM FDR with FDR:protein > 0 (e.g. "
                       "0.01), a loose evidence threshold of the intermediates with FDR:protein 0 (e.g. 0.1; see there).");
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
    defaults_.setValue("FDR:protein", 0.0, "Picked-protein FDR threshold (e.g. 0.01 = 1 %), and the switch between the two modes of a search. Above 0, the search is the whole experiment: its proteins are inferred on the pooled PSMs of all its inputs, filtered at this protein FDR, and the PSMs filtered finally, so every output is the final result (the decoys removed, see the ProSE tool's -keep_decoys). With 0 (the default), the search is one part of a larger experiment, e.g. one shard of files searched on one compute node: every output is a per-run intermediate (PSM q-values and PEPs per run, FDR:PSM as a loose evidence threshold, decoys and the proteins the search indexed for the run's PSMs kept), to be finished downstream on all parts together (e.g. ProteomicsLFQ, or IDMerger, ProteinInference or Epifany, FalseDiscoveryRate and IDFilter). Uses the picked-protein approach (Savitski et al. 2015) which pairs target and decoy proteins by accession. The protein inference and the picked-protein FDR use the PSMs up to FDR:protein_inference_PSM, and a protein scores with the best posterior error probability of its PSMs (Percolator's, with 'rescore'; without it the best search engine score). The proteins and the protein groups (which have q-values of their own) are filtered at the threshold, and only the proteins of the groups that pass are kept. Proteins without a reported PSM (one that passes FDR:PSM) are removed afterwards, and so are the PSMs all of whose proteins failed the protein FDR. Needs decoys ('decoys' auto or generate).");
    defaults_.setMinFloat("FDR:protein", 0.0);
    defaults_.setMaxFloat("FDR:protein", 1.0);
    defaults_.setValue("FDR:protein_inference_PSM", 0.1,
                       "q-value cutoff of the PSMs that protein inference and the protein FDR (FDR:protein > 0) use; FDR:PSM "
                       "still filters the reported PSMs. Without the protein FDR it has no effect. The q-values are those FDR:PSM filters on (Percolator's with 'rescore', else "
                       "see FDR:PSM_groups). A looser cutoff than FDR:PSM leaves more decoy proteins for the picked-protein FDR, "
                       "as ProteomicsLFQ infers proteins on PSMs at a looser q-value than its PSM FDR. Values at or below FDR:PSM "
                       "use the PSMs of FDR:PSM; with FDR:PSM 0 (no PSM filter) all PSMs are used.",
                       {"advanced"});
    defaults_.setMinFloat("FDR:protein_inference_PSM", 0.0);
    defaults_.setMaxFloat("FDR:protein_inference_PSM", 1.0);
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
      "Measure the mass error of the run and search with windows derived from it: a fast calibration pass on the most "
      "intense spectra (calibration:subset_ratio), searched with at least 20 ppm (in Da: twice the configured window, up "
      "to 0.5 Da) precursor and fragment windows (the configured ones if wider) so the error is measured and not censored, takes the offset (median) and the "
      "spread (robust standard deviation, 1.4826 * MAD) of the signed mass errors of its confident PSMs; the main search "
      "then uses the fragment tolerance |median| + k sigma (calibration:sigma_factor k) within "
      "[calibration:min_fragment_tolerance, the configured tolerance], and the precursor window [median - k sigma, "
      "median + k sigma] only with calibration:precursor_sigma_factor k > 0 (by default the precursor errors are "
      "measured and reported, the configured window kept): the configured tolerances are the widest windows searched. "
      "The pass takes its spectra per tenth of the run; a window that covers less than 98 % of the calibration PSMs "
      "(precursor; heavy tails or a mass offset that changes over the run) or 95 % of their matched fragment peaks "
      "(fragment) is not trusted, and the configured tolerance is kept with a warning. "
      "The fragment index is not rebuilt, only the query-time tolerances change. With too few PSMs "
      "(calibration:min_psms) the configured tolerances are used. 'false' searches with the configured tolerances as given.");
    defaults_.setValidStrings("calibration:enabled", {"true", "false"});
    defaults_.setValue("calibration:subset_ratio", 0.1,
      "Fraction of spectra (by TIC, highest first) used for the calibration pass (0.0-1.0).");
    defaults_.setMinFloat("calibration:subset_ratio", 0.01);
    defaults_.setMaxFloat("calibration:subset_ratio", 1.0);
    defaults_.setValue("calibration:min_psms", 50,
      "Minimum number of confident PSMs required for calibration. If fewer are found, "
      "calibration is skipped and the user-configured tolerances are used.");
    defaults_.setMinInt("calibration:min_psms", 1);
    defaults_.setValue("calibration:sigma_factor", 4.0,
      "Half-width of the calibrated fragment tolerance in robust standard deviations (1.4826 * MAD) of the signed errors "
      "of the matched fragment peaks of the calibration PSMs: |median| + k sigma. 4 covers practically all correct "
      "matches of a normal error distribution.");
    defaults_.setMinFloat("calibration:sigma_factor", 1.0);
    defaults_.setValue("calibration:precursor_sigma_factor", 0.0,
      "Half-width of the calibrated precursor window in robust standard deviations of the precursor errors of the "
      "calibration PSMs, around their median (asymmetric about zero when the run has an offset). 0 keeps the configured "
      "precursor window and only measures and reports the errors: precursor yields are flat from 5 to 50 ppm, and a "
      "tightened precursor window cost identifications on every Orbitrap HCD set of the benchmark (heavy tails), while "
      "the fragment window gains. A value such as 4 tightens it, within calibration:min_precursor_tolerance and the "
      "configured window, when the window covers 98 % of the calibration PSMs.");
    defaults_.setMinFloat("calibration:precursor_sigma_factor", 0.0);
    defaults_.setValue("calibration:min_precursor_tolerance", 0.0,
      "Smallest calibrated precursor tolerance per side, in precursor:mass_tolerance_unit (0: 2.5 ppm or 0.005 Da). "
      "The calibrated window never exceeds the configured one.");
    defaults_.setMinFloat("calibration:min_precursor_tolerance", 0.0);
    defaults_.setValue("calibration:min_fragment_tolerance", 0.0,
      "Smallest calibrated fragment tolerance, in fragment:mass_tolerance_unit (0: 10 ppm or 0.02 Da; 5 ppm or 0.01 Da "
      "with calibration:recalibrate). The calibrated tolerance never exceeds the configured one.");
    defaults_.setMinFloat("calibration:min_fragment_tolerance", 0.0);
    defaults_.setValue("calibration:recalibrate", "none",
      "Correct the m/z of the spectra by the mass errors the calibration pass measures, before the main search: "
      "'offset' subtracts the median error, 'linear' a line in m/z fitted to the errors (precursor errors against "
      "the precursor m/z, fragment errors against the fragment m/z), each per retention-time chunk "
      "(calibration:recalibrate_rt_chunks) and interpolated between the chunks. The precursor m/z and the fragment "
      "peaks of the searched MS2 spectra are corrected in place (the candidate retrieval, the scoring and the PSMs' "
      "mass errors use the corrected values; a PeakMap passed to the library is corrected too); the PSMs report the "
      "observed precursor m/z and their fragment annotations the observed peak m/z. The windows are then derived from "
      "the residual errors, with a fragment floor of 5 ppm (0.01 Da) unless calibration:min_fragment_tolerance is set. "
      "The co-isolated precursors of a chimeric search, from the same survey scans, are corrected by the precursor "
      "model. 'none' measures only.");
    defaults_.setValidStrings("calibration:recalibrate", {"none", "offset", "linear"});
    defaults_.setValue("calibration:recalibrate_rt_chunks", 1,
      "Retention-time chunks of the recalibration models (calibration:recalibrate): the calibration PSMs are split in "
      "acquisition order into this many groups of equal size (fewer when a group would hold less than "
      "calibration:min_psms PSMs), one model each, interpolated between the groups' median retention times. 1: one "
      "model for the run; 10 follows an offset that changes over the run (the pass samples each tenth of it).");
    defaults_.setMinInt("calibration:recalibrate_rt_chunks", 1);

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

    defaults_.setValue("chimeric:enabled", "false",
      "Also identify precursors co-isolated with the selected one (after the andes chimeric cascade). The survey scan of "
      "each MS2 spectrum is searched for isotope envelopes in its isolation window that fit averagine (needs centroided "
      "MS1 spectra in the input and isolation windows on the precursors). After the search, the peaks the best hit of the "
      "spectrum explains are removed, and the rest is searched for each co-isolated precursor in turn, at its m/z and "
      "charge; the best hit of each removes the peaks it explains before the next is searched. The best hit's peptide is "
      "not a hit of a co-isolated precursor, and a co-isolated precursor whose best hit repeats a PSM of the spectrum "
      "(same peptide and charge) is dropped. These PSMs get the spectrum reference of their spectrum and the PSM meta "
      "value 'chimeric_rank' 2 (the others 1). Not with database:chunk_size > 0.");
    defaults_.setValidStrings("chimeric:enabled", {"true", "false"});
    defaults_.setValue("chimeric:max_coisolated", 4, "At most this many co-isolated precursors per spectrum, the most intense first.");
    defaults_.setMinInt("chimeric:max_coisolated", 1);
    defaults_.setValue("chimeric:max_kl", 0.3,
      "Largest Kullback-Leibler divergence of an isotope envelope in the survey scan from averagine (first four isotope "
      "peaks) to accept it as a co-isolated precursor.");
    defaults_.setMinFloat("chimeric:max_kl", 0.0);
    defaults_.setValue("chimeric:top_hits", 3,
      "At most this many hits of a co-isolated precursor are reported (and at most report:top_hits). Each takes about as "
      "much memory as a hit of a spectrum.");
    defaults_.setMinInt("chimeric:top_hits", 1);
    defaults_.setValue("chimeric:mz_tolerance", 10.0, "Tolerance (ppm) for the m/z of the isotope peaks in the survey scan.");
    defaults_.setMinFloat("chimeric:mz_tolerance", 0.0);
    defaults_.setSectionDescription("chimeric", "Search of precursors co-isolated with the selected one (chimeric spectra).");

    defaults_.setSectionDescription("calibration",
      "Automatic mass tolerances (two-pass search). A fast first pass on the most intense spectra measures the offset and "
      "spread of the precursor and fragment mass errors of the run; the main search uses windows derived from them, "
      "never wider than the configured tolerances.");

    FragmentIndex::restrictSearchParameters(defaults_); // the bounds of the parameters it shares with the FragmentIndex

    defaultsToParam_();
  }

  void ProSEAlgorithm::updateMembers_()
  {
    // inverted ranges and an invalid peptide:motif, before anything is read (the FragmentIndex checks them again)
    FragmentIndex::checkSearchParameters(param_, parameter_name_prefix_);

    precursor_mass_tolerance_lower_ = param_.getValue("precursor:mass_tolerance_lower");
    precursor_mass_tolerance_upper_ = param_.getValue("precursor:mass_tolerance_upper");
    precursor_mass_tolerance_unit_ = param_.getValue("precursor:mass_tolerance_unit").toString();

    precursor_min_charge_ = param_.getValue("precursor:min_charge");
    precursor_max_charge_ = param_.getValue("precursor:max_charge");

    peptdeep_enable_ = param_.getValue("peptdeep:enable").toString() == "true";
    peptdeep_ms2_model_ = param_.getValue("peptdeep:ms2_model").toString();
    peptdeep_rt_model_ = param_.getValue("peptdeep:rt_model").toString();
    peptdeep_instrument_ = param_.getValue("peptdeep:instrument").toString();
    peptdeep_nce_ = param_.getValue("peptdeep:nce");
    peptdeep_rt_model_type_ = param_.getValue("peptdeep:rt_model_type").toString();
    peptdeep_batch_size_ = param_.getValue("peptdeep:batch_size");

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
        paramName_("fragment:mass_tolerance") + " must be finite and positive (got " + std::to_string(fragment_mass_tolerance_) + ").");
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
        paramName_("fragment:deisotope") + "=true requires a high-resolution fragment tolerance (<= 0.1 Da or <= 100 ppm). "
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
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, paramName_("scoring:mass_error_sd") + " must be finite and positive.");
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
    FragmentIndex::checkFixedModifications(modifications_fixed_, parameter_name_prefix_);

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<std::string> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      OPENMS_LOG_WARN << "Duplicate variable modification provided. Making them unique." << endl;
      modifications_variable_.assign(var_unique.begin(), var_unique.end());
    }
    // A modification that is fixed as well applies everywhere already: the remedy is to drop the fixed one.
    const auto also_fixed = [this](const std::string& mod)
    {
      return std::find(modifications_fixed_.begin(), modifications_fixed_.end(), mod) != modifications_fixed_.end();
    };
    for (const std::string& mod : FragmentIndex::shadowedVariableTerminalModifications(modifications_fixed_, modifications_variable_))
    {
      if (also_fixed(mod))
      {
        OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: it is also a fixed modification, which "
                        << "applies to every such terminus. To search the terminus with and without it, remove it from "
                        << "modifications:fixed." << endl;
        continue;
      }
      OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: a fixed modification already sits on "
                      << "that terminus, which carries one modification. To search it, specify the fixed terminal "
                      << "modification as a variable one as well." << endl;
    }
    for (const std::string& mod : FragmentIndex::shadowedVariableResidueModifications(modifications_fixed_, modifications_variable_))
    {
      if (also_fixed(mod))
      {
        OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: it is also a fixed modification, which "
                        << "applies to every such residue. To search the residue with and without it, remove it from "
                        << "modifications:fixed." << endl;
        continue;
      }
      OPENMS_LOG_WARN << "Variable modification '" << mod << "' is not searched: a fixed modification already sits on "
                      << "that residue, which carries one modification. To search it, specify the fixed residue "
                      << "modification as a variable one instead." << endl;
    }

    modifications_max_variable_mods_per_peptide_ = param_.getValue("modifications:variable_max_per_peptide");

    enzyme_ = param_.getValue("enzyme").toString();

    peptide_min_size_ = param_.getValue("peptide:min_size");
    peptide_max_size_ = param_.getValue("peptide:max_size");
    peptide_missed_cleavages_ = param_.getValue("peptide:missed_cleavages");
    peptide_enzyme_specificity_ = EnzymaticDigestion::getSpecificityByName(
      param_.getValue("peptide:enzyme_specificity").toString());

    report_top_hits_ = param_.getValue("report:top_hits");
    rescore_top_hits_ = static_cast<Size>(static_cast<int>(param_.getValue("rescoring:top_hits")));
    isotope_error_observed_minus_theoretical_ =
      param_.getValue("report:isotope_error_convention").toString() == "observed_minus_theoretical";

    const std::string decoy_mode_str = param_.getValue("decoys").toString();
    if (decoy_mode_str == "generate")   { decoy_mode_ = DecoyMode_::GENERATE; }
    else if (decoy_mode_str == "ignore") { decoy_mode_ = DecoyMode_::IGNORE; }
    else                                 { decoy_mode_ = DecoyMode_::AUTO; }
    decoy_prefix_ = param_.getValue("decoy_prefix").toString();
    decoy_keep_termini_ = param_.getValue("decoy_method").toString() == "keep_termini";
    decoy_parity_report_ = param_.getValue("decoy_parity_report").toBool();
    terminal_features_ = param_.getValue("rescoring:terminal_features").toString() == "true";
    stratum_q_values_ = param_.getValue("rescoring:stratum_q_values").toString() == "true";
    annotate_psm_ = ListUtils::toStringList<std::string>(param_.getValue("annotate:PSM"));
    if (peptdeep_enable_
        && std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") == annotate_psm_.end()
        && std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) == annotate_psm_.end())
    {
      // PeptDeepRescoring compares the predicted with the annotated fragment peaks of each PSM
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        paramName_("peptdeep:enable") + " needs the fragment annotation of the PSMs: " + paramName_("annotate:PSM") + " must contain 'ALL' or '"
        + std::string(Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) + "'.");
    }
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
    rescore_ = param_.getValue("rescore").toString() == "true";
    max_train_psms_ = static_cast<Size>(static_cast<int>(param_.getValue("percolator:max_train_psms")));
    percolator_training_.subset_max_train = static_cast<int>(std::min<Size>(max_train_psms_, std::numeric_limits<int>::max()));
    percolator_training_.subset_mode = param_.getValue("rescoring:subset_mode").toString();
    percolator_training_.curve_start = static_cast<int>(param_.getValue("rescoring:curve_start"));
    percolator_training_.curve_min_gain = static_cast<double>(param_.getValue("rescoring:curve_min_gain"));
    percolator_training_.iterations = static_cast<int>(param_.getValue("rescoring:iterations"));
    percolator_training_.tron_epsilon = static_cast<double>(param_.getValue("rescoring:tron_epsilon"));
    percolator_training_.warm_start = param_.getValue("rescoring:warm_start").toString() == "true";
    percolator_training_.single_cneg = param_.getValue("rescoring:single_cneg").toString() == "true";
    percolator_training_.solver = param_.getValue("rescore_solver").toString();
    percolator_training_.use_pi0 = param_.getValue("rescoring:use_pi0").toString();
    percolator_training_.calc_mass_feature = param_.getValue("rescoring:calc_mass_feature").toString() == "true";
    percolator_training_.fixed_enzyme_termini = param_.getValue("rescoring:fixed_enzyme_termini").toString() == "true";
    percolator_training_.peptide_level = param_.getValue("rescoring:peptide_level").toString() == "true";
    percolator_extra_ = Param();
    if (param_.getValue("rescore_method").toString() == "reset")
    {
      percolator_extra_.setValue("algorithm", "reset");
      for (const std::string key : {"salt", "learner", "audit", "stem_table", "representative", "mlp_crossfit", "decoys",
                                    "second_decoy_prefix"})
      {
        percolator_extra_.setValue("reset:" + key, param_.getValue("reset:" + key).toString());
      }
      for (const std::string key : {"iterations", "mlp_members", "mlp_hidden", "mlp_epochs", "mlp_warm_epochs", "mlp_norm_min_decoys"})
      {
        percolator_extra_.setValue("reset:" + key, static_cast<int>(param_.getValue("reset:" + key)));
      }
      if (static_cast<double>(param_.getValue("FDR:PSM")) > 0.0)
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_method 'reset' with FDR:PSM > 0: the PSM q-values that FDR:PSM filters on are "
                           "heuristic under RESET (FDR is controlled at the peptide level: RESET:peptide_q_value)." << endl;
      }
    }
    reuse_mode_ = param_.getValue("rescore_reuse:mode").toString();
    reuse_max_train_ = static_cast<Size>(static_cast<int>(param_.getValue("rescore_reuse:max_train")));
    reuse_static_files_ = static_cast<Size>(static_cast<int>(param_.getValue("rescore_reuse:static_files")));
    reuse_max_shift_ = param_.getValue("rescore_reuse:max_shift");
    if (percolator_extra_.exists("algorithm") && reuse_mode_ != "per_file")
    {
      OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode '" << reuse_mode_ << "' does not apply with rescore_method 'reset' (the "
                         "reused models are Percolator models): every file is rescored by RESET ('per_file')." << endl;
      reuse_mode_ = "per_file";
    }
    fdr_psm_ = param_.getValue("FDR:PSM");
    fdr_psm_by_scored_charges_ = param_.getValue("FDR:PSM_groups").toString() == "scored_charges";
    fdr_protein_ = param_.getValue("FDR:protein");
    fdr_protein_inference_psm_ = param_.getValue("FDR:protein_inference_PSM");

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
    chimeric_enabled_ = param_.getValue("chimeric:enabled").toBool();

    calibration_enabled_ = param_.getValue("calibration:enabled") == "true";
    calibration_subset_ratio_ = param_.getValue("calibration:subset_ratio");
    calibration_min_psms_ = param_.getValue("calibration:min_psms");
    calibration_sigma_factor_ = param_.getValue("calibration:sigma_factor");
    calibration_precursor_sigma_factor_ = param_.getValue("calibration:precursor_sigma_factor");
    calibration_min_precursor_tol_ = param_.getValue("calibration:min_precursor_tolerance");
    calibration_min_fragment_tol_ = param_.getValue("calibration:min_fragment_tolerance");
    calibration_recalibrate_ = param_.getValue("calibration:recalibrate").toString();
    calibration_recalibrate_rt_chunks_ = static_cast<Size>(static_cast<Int>(param_.getValue("calibration:recalibrate_rt_chunks")));
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
    last_percolator_ = PercolatorRescoring{};
    inference_only_peptide_ids_ = PeptideIdentificationList();

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
      reportDecoyParity_(ctx);
      return search(spectra, ctx, protein_ids, peptide_ids);
    }
    return searchChunked_(spectra, full_db, strategy, protein_ids, peptide_ids);
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
    // their SearchResult::stats after search() returns OK. The same for the pre-FDR PSMs: a copy
    // left by an earlier search must not end up in this file's result.
    last_run_stats_ = RunStatistics{};
    last_percolator_ = PercolatorRescoring{};
    pre_fdr_peptide_ids_ = PeptideIdentificationList();
    inference_only_peptide_ids_ = PeptideIdentificationList();

    // chimeric:enabled: find the precursors co-isolated with the selected one in the survey scans, which are not
    // searched themselves
    if (chimeric_enabled_ && !chimericSearch_())
    {
      OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
    }
    if (chimericSearch_())
    {
      std::vector<MSSpectrum>& all = spectra.getSpectra();
      const bool survey_scans = std::any_of(all.begin(), all.end(), [](const MSSpectrum& s) { return s.getMSLevel() == 1; });
      if (survey_scans)
      {
        const Size found = CoIsolationDetector::annotate(spectra, precursor_min_charge_, precursor_max_charge_,
                                                         param_.getValue("chimeric:mz_tolerance"), true,
                                                         param_.getValue("chimeric:max_kl"),
                                                         static_cast<Size>(static_cast<int>(param_.getValue("chimeric:max_coisolated"))));
        all.erase(std::remove_if(all.begin(), all.end(), [](const MSSpectrum& s) { return s.getMSLevel() == 1; }), all.end());
        OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << found << " co-isolated precursors in the isolation windows of "
                        << spectra.size() << " spectra." << std::endl;
      }
      else if (std::none_of(all.begin(), all.end(), [](const MSSpectrum& s) { return s.metaValueExists(CoIsolationDetector::META_MZ); }))
      {
        OPENMS_LOG_WARN << "[ProSE] chimeric:enabled needs the survey (MS1) scans, which the spectra do not include;"
                        << " no co-isolated precursors are searched." << std::endl;
      }
    }

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    bool open_search = isOpenSearchMode_();
    OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search ? "true" : "false")
                    << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] "
                    << precursor_mass_tolerance_unit_ << ")" << std::endl;

    takeBackThreadsOfEndedRead(); // here and between the other phases: see BackgroundSpectra
    startProgress(0, 1, "Filtering spectra...");
    FragmentIonLikelihoodModel::PeakLists ion_evidence; // annotate:self_trained_ion_priors
    PeakMap evidence_spectra, query_spectra;
    PeakMap* query_ptr = query_raw_spectrum_ ? &query_spectra : nullptr;
    PeakMap* evidence_ptr = param_.getValue("annotate:local_fragment_evidence").toBool() ? &evidence_spectra : nullptr;
    preprocessSearchedSpectra_(spectra, "", ion_evidence, evidence_ptr, query_ptr);
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
    // Restores the FragmentIndex parameters AND the algo-level tolerance members if calibration modified them, so the
    // shared SearchContext and the algorithm instance are both clean for subsequent per-file searches in the
    // multi-file wrapper: on the return paths below, and when an exception leaves this function (e.g. a worker
    // exception that the scoring or calibration loop rethrows) for a caller that catches it.
    auto restore_fi_params = [&]()
    {
      if (fi_params_modified)
      {
        fi_params_modified = false;
        precursor_mass_tolerance_lower_ = orig_precursor_mass_tolerance_lower;
        precursor_mass_tolerance_upper_ = orig_precursor_mass_tolerance_upper;
        fragment_index_.setParameters(fi_params_original);
      }
    };
    const OnScopeExit restore_on_exit([&restore_fi_params]() noexcept
    {
      try { restore_fi_params(); }
      catch (...) {} // only after another exception (e.g. std::bad_alloc copying the parameters): the members are restored
    });

    // --- Optional calibration pass ---
    last_calibration_result_ = CalibrationResult_(); // (of this search only; recorded in the search parameters below)
    if (calibration_enabled_ && !open_search)
    {
      startProgress(0, 1, "Running calibration pass...");
      StopWatch sw_cal; sw_cal.start();
      {
        OnScopeExit end_progress([this]() { endProgress(); }); // also when the calibration loop rethrows
        last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db, query_ptr);
      }
      sw_cal.stop();
      last_run_stats_.seconds_calibration = sw_cal.getClockTime();
      const CalibrationResult_& cal = last_calibration_result_;

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
          OPENMS_LOG_WARN << "[ProSE] Calibration: the precursor errors have (almost) zero spread (sigma="
                          << cal.precursor_sigma << " " << precursor_mass_tolerance_unit_
                          << "): the configured precursor window is kept. Fragment calibration is still applied." << std::endl;
          // Precursor bounds unchanged; fragment_tolerance applied above.
        }

        fragment_index_.setParameters(fi_params);
        fi_params_modified = true;
        applyRecalibration_(cal, spectra, query_ptr, evidence_ptr);
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

    // storage for the PSMs of every spectrum (each grows on the thread that scores the spectrum)
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    bool open_search_mode = open_search;

    StopWatch sw_search; sw_search.start();
    scoreSpectraAgainstIndex_(spectra, fragment_index_, db, generators,
                              effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                              open_search_mode, annotated_hits, pool_stats,
                              "Scoring peptide models against spectra...", query_ptr);
    freeInParallel(query_spectra);

    // chimeric:enabled: search the co-isolated precursors in what the other PSMs of their spectrum leave of it. Their
    // hits and pool statistics are appended after those of the searched spectra; their residual spectra are built
    // only after the fragment index is released (below), and removed again after the post-processing.
    const Size n_searched_spectra = spectra.size();
    std::vector<CoIsolatedPSM_> coisolated;
    std::vector<Size> coisolated_parent;
    if (chimericSearch_())
    {
      std::vector<std::vector<AnnotatedHit_>> coisolated_hits;
      std::vector<CandidatePoolStats_> coisolated_stats;
      searchCoIsolated_(spectra, annotated_hits, fragment_index_, db, generators, effective_fragment_tol,
                        fragment_mass_tolerance_unit_ppm, open_search_mode, coisolated, coisolated_hits, coisolated_stats);
      annotated_hits.reserve(annotated_hits.size() + coisolated.size());
      pool_stats.reserve(pool_stats.size() + coisolated.size());
      coisolated_parent.reserve(coisolated.size());
      for (Size k = 0; k < coisolated.size(); ++k)
      {
        annotated_hits.push_back(std::move(coisolated_hits[k]));
        pool_stats.push_back(std::move(coisolated_stats[k]));
        coisolated_parent.push_back(coisolated[k].parent);
      }
    }

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

    takeBackThreadsOfEndedRead();
    if (!coisolated.empty())
    {
      spectra.reserve(spectra.size() + coisolated.size());
      for (CoIsolatedPSM_& psm : coisolated)
      {
        MSSpectrum residual = residualSpectrum_(spectra[psm.parent], psm.removed, psm.precursor);
        spectra.addSpectrum(std::move(residual));
        psm = CoIsolatedPSM_();
      }
      std::vector<CoIsolatedPSM_>().swap(coisolated);
    }

    startProgress(0, 1, "Post-processing PSMs...");
    // The fragment annotations of a hit are needed only for the output: when the rescoring (all but the best hit of a
    // spectrum) or the PSM FDR drops hits, they are made after that, for the hits that remain. Not if anything reads
    // them before (PeptDeep features) or keeps the hits of before (the pre-FDR copy). A search that fails in the protein
    // mapping makes them before it returns its error.
    // Only for hits of this search (identifications already in @p peptide_ids keep theirs), found again by their
    // scan_index (stored as unsigned int).
    // Not with chimeric:enabled: its residual spectra, which annotate the PSMs of co-isolated precursors, are gone by then.
    const bool defer_peak_annotations = (rescore_ || fdr_psm_ > 0.0) && !peptdeep_enable_ && !keep_pre_fdr_psms_
      && !chimericSearch_() && peptide_ids.empty() && spectra.size() <= std::numeric_limits<unsigned int>::max();
    ProSEAlgorithm::postProcessHits_(spectra,
      annotated_hits,
      pool_stats,
      protein_ids,
      peptide_ids,
      searchTopHits_(),
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
      evidence_ptr,
      &coisolated_parent, // the residual spectra use the peak lists of their spectra
      defer_peak_annotations,
      last_calibration_result_.appliedFragmentRecalibration()
      );
    recordRecalibration_(last_calibration_result_, protein_ids);
    // Release the hits and the additional peak list before PeptideIndexing.
    freeOnOwnerThreads(annotated_hits, [&pool_stats](SignedSize i) { return pool_stats[i].scoring_thread; }); // (where scored)
    freeInParallel(pool_stats);
    freeInParallel(evidence_spectra);
    endProgress();

    takeBackThreadsOfEndedRead();
    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);

    // chimeric:enabled: the PSMs of co-isolated precursors refer to the spectrum they were found in (scan_index and
    // spectrum reference), next to its own PSMs, and are told apart by chimeric_rank
    if (chimericSearch_())
    {
      Size n_coisolated_psms = 0;
      for (PeptideIdentification& pi : peptide_ids)
      {
        bool coisolated = false;
        if (pi.metaValueExists("scan_index"))
        {
          const Size scan = static_cast<unsigned int>(pi.getMetaValue("scan_index"));
          if (scan >= n_searched_spectra && scan - n_searched_spectra < coisolated_parent.size())
          {
            pi.setMetaValue("scan_index", static_cast<unsigned int>(coisolated_parent[scan - n_searched_spectra]));
            coisolated = true;
            ++n_coisolated_psms;
          }
        }
        for (PeptideHit& hit : pi.getHits()) { hit.setMetaValue(CHIMERIC_RANK, coisolated ? 2 : 1); }
      }
      std::stable_sort(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return static_cast<unsigned int>(a.getMetaValue("scan_index", 0u)) < static_cast<unsigned int>(b.getMetaValue("scan_index", 0u));
      });
      std::vector<MSSpectrum>& all = spectra.getSpectra();
      all.erase(all.begin() + static_cast<SignedSize>(n_searched_spectra), all.end());
      OPENMS_LOG_INFO << "[ProSE] Chimeric search: " << n_coisolated_psms << " PSMs of co-isolated precursors." << std::endl;
    }
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

    if ((indexer_exit != PeptideIndexing::ExitCodes::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::ExitCodes::PEPTIDE_IDS_EMPTY))
    {
      // (the hits returned with the error carry their fragment annotations, as before the deferral)
      if (defer_peak_annotations) annotateFragmentPeaks_(spectra, peptide_ids, effective_fragment_tol, fragment_mass_tolerance_unit_,
                                                     last_calibration_result_.appliedFragmentRecalibration());
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
    takeBackThreadsOfEndedRead();
    if (self_trained_ion_priors_)
    {
      annotateIonPriors_(spectra, ion_evidence, protein_ids, peptide_ids);
      ion_evidence.clear();
    }

    // Capture pre-FDR stats now — BEFORE FDR, which drops pure-decoy hits
    // (FalseDiscoveryRate default add_decoy_peptides=false), may strip decoys
    // entirely, and overwrites each hit's HyperScore with its q-value.
    capturePreFdrStats_(peptide_ids, last_run_stats_);

    // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
    last_percolator_ = applyPsmFDR_(protein_ids, peptide_ids, has_decoys, last_run_stats_,
                                    keepPreRescoringPsms_() ? &pre_fdr_peptide_ids_ : nullptr, &inference_only_peptide_ids_);
    if (defer_peak_annotations)
    {
      takeBackThreadsOfEndedRead();
      // (counted as search time, as when postProcessHits_ made them)
      StopWatch sw_annotations;
      sw_annotations.start();
      annotateFragmentPeaks_(spectra, peptide_ids, effective_fragment_tol, fragment_mass_tolerance_unit_,
                             last_calibration_result_.appliedFragmentRecalibration());
      // A file that Percolator did not rescore on its own may end with its PSMs of before the rescoring
      // (rescoreFilesTogether_()): they get theirs too.
      if (pool_rescoring_ && !last_percolator_.applied)
      {
        annotateFragmentPeaks_(spectra, pre_fdr_peptide_ids_, effective_fragment_tol, fragment_mass_tolerance_unit_,
                               last_calibration_result_.appliedFragmentRecalibration());
      }
      sw_annotations.stop();
      last_run_stats_.seconds_search += sw_annotations.getClockTime();
    }

    restore_fi_params();

    takeBackThreadsOfEndedRead();
    collectRunStatistics_(spectra, protein_ids, peptide_ids, last_run_stats_);

    return ExitCodes::EXECUTION_OK;
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
    // With more threads, the single-file search always waits (and helps with the rest of the read, see
    // readerThreads()): the read then mostly ends while the peptides are generated and sorted, and the fragment passes over the third or so of the peptides that no spectrum reaches take longer
    // than the wait (measured at 16 threads on Orbitrap, Astral, timsTOF and TMT files of 8,000 spectra: 0.1-0.2 s
    // shorter per file; waiting for a read with a quarter of the threads was up to 0.1 s slower for some).
    constexpr Size MAX_THREADS_WAITING_FOR_SPECTRA = 4;
    constexpr UInt64 MAX_SPECTRA_BYTES_PER_PEPTIDE = 100;

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
      const std::string& in_spectra_location, const std::string& in_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    // The name of the spectrum file, as it is logged and stored: of a remote location (e.g. s3://...) without its
    // credentials and query (FileHandler::displayName()). in_spectra_location is what the readers are given.
    const std::string in_spectra = FileHandler::displayName(in_spectra_location);
    reset_run_name_ = in_spectra; // names the RESET stem table of the run (reset:stem_table a directory)
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
      BackgroundSpectra spectra_ready(in_spectra_location, false, chimericSearch_());

      // load FASTA
      loadFASTA(in_db, fasta_db);

      // ions:by_activation: the index needs c and z+1 ions if electron-activated spectra are
      // searched, which is known once the spectra are read. Wait for them if they are read already
      // (or failed to be read: get() rethrows) or if the file names an electron-based activation
      // early on. Otherwise build the index without these ions meanwhile, which is what most data
      // need, and rebuild it should such spectra turn up after all (as the multi-file search does
      // when a later file has them): the index searched is the one prepareContext(fasta_db, true)
      // builds.
      const bool spectra_read = spectra_ready.ready()
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
            if (!wait && !spectra_ready.ready())
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
          buildFragmentIndex_(ctx, fragmentIndexParameters_(true), searched_spectra);
          ctx.electron_ions = true;
          endProgress();
        }
      }
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      ec = search(spectra, ctx, protein_ids, peptide_ids);
    }
    else
    {
      spectra = loadMS2Spectra(in_spectra_location, readerThreads(false), nullptr, 0, chimericSearch_());

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
        applyCompleteSetProteinFDR(protein_ids, peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_,
                                   true, inference_only_peptide_ids_);
      }
    }
    inference_only_peptide_ids_ = PeptideIdentificationList();

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
    result.percolator = last_percolator_;
    result.inference_only_peptide_ids = std::move(inference_only_peptide_ids_);
    inference_only_peptide_ids_ = PeptideIdentificationList();

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

    // A single input is a complete experiment: protein FDR as in search(file).
    if (fdr_protein_ > 0.0 && !result.protein_ids.empty() && !result.peptide_ids.empty())
    {
      const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
      if (strategy.have_decoys)
      {
        applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_,
                                   true, result.inference_only_peptide_ids);
        updateFinalStats(result.stats, result.peptide_ids, enzyme_, result.stats.fdr_applied);
      }
    }
    result.inference_only_peptide_ids = PeptideIdentificationList(); // used by the protein FDR only

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
    return searchWithModificationAnalysis_(in_spectra_files, fasta_db, nullptr, output_base_names, aggregate_base_name,
                                           build_pooled_aggregate);
  }

  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis_(
      const std::vector<std::string>& in_spectra_locations,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db_given,
      const std::string* in_db,
      const std::vector<std::string>& output_base_names,
      const std::string& aggregate_base_name,
      bool build_pooled_aggregate) const
  {
    // The names of the spectrum files, as they are logged and stored: of a remote location (e.g. s3://...) without
    // its credentials and query (FileHandler::displayName()). in_spectra_locations is what the readers are given.
    std::vector<std::string> in_spectra_files;
    in_spectra_files.reserve(in_spectra_locations.size());
    for (const std::string& location : in_spectra_locations) { in_spectra_files.push_back(FileHandler::displayName(location)); }
    // The database: fasta_db_given, or read from *in_db below, once the first spectrum file is read in the background.
    std::vector<FASTAFile::FASTAEntry> fasta_db_read;
    const std::vector<FASTAFile::FASTAEntry>& fasta_db = in_db == nullptr ? fasta_db_given : fasta_db_read;
    if (!output_base_names.empty() && output_base_names.size() != in_spectra_files.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "output_base_names must be empty or have exactly one entry per spectrum file (got "
        + StringUtils::toStr(output_base_names.size()) + " entries for " + StringUtils::toStr(in_spectra_files.size())
        + " spectrum files).");
    }

    MultiFileSearchResult mfres;
    mfres.aggregate.is_open_search = isOpenSearchMode_();
    released_hits_.reset(); // the releases of a search before (see the end of the loop below), so that they do not pile up
    // 'rescore' with several files: keep each file's PSMs before the rescoring until all files are searched, for
    // rescoreFilesTogether_() (in case Percolator cannot rescore some file on its own)
    pool_rescoring_ = rescore_ && in_spectra_files.size() > 1;
    pool_files_ = in_spectra_files.size();
    const OnScopeExit end_pool_rescoring([this]() { pool_rescoring_ = false; pool_files_ = 0; });

    if (in_spectra_files.empty())
    {
      if (in_db != nullptr) { loadFASTA(*in_db, fasta_db_read); } // reports a database that cannot be read, as before
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
    BackgroundSpectra next_spectra; // the spectra of the next file to search, if read in the background
    // A single file is the only search of the context (r3 merge of WP7 and WP10): as in search(file), the index
    // holds only the peptides its spectra can reach if they are read by the time the peptides are generated, or
    // if waiting for them pays off (see MAX_THREADS_WAITING_FOR_SPECTRA); the spectra taken for that are searched.
    const bool restrict_index = in_spectra_files.size() == 1 && restrictIndexToSpectra_();
    // With more threads than MAX_THREADS_WAITING_FOR_SPECTRA, the index build waits for the spectra of such a file
    // (and helps to read them, see BackgroundSpectra::get()).
    const bool wait_for_first_spectra = restrict_index && threads_started > MAX_THREADS_WAITING_FOR_SPECTRA;
    const auto read_next = [&](Size i)
    {
      if (read_in_background && i < in_spectra_files.size() && is_mzml(i))
      {
        // File 0 next to the index build. A later one next to the search of the file before: from its rescoring on if
        // Percolator, which trains on at most 9 threads, leaves the reader's threads (a quarter) idle; otherwise from the
        // start (measured on three Orbitrap files at 16 threads without rescoring: from the end of the scoring on, 3%
        // slower, as more of the read is left for get()). The reader then registers the names of the file's meta values
        // after the protein mapping of the file before (see the registry order above).
        if (i == 0) { next_spectra.start(in_spectra_locations[i], wait_for_first_spectra, 0, chimericSearch_()); }
        else if (rescore_ && threads_started >= 12) { next_spectra.startDeferred(in_spectra_locations[i], chimericSearch_()); }
        else { next_spectra.start(in_spectra_locations[i], false, readerThreads(true, false), chimericSearch_()); }
      }
    };
    read_next(0);
    // Whether the first file names an electron-based activation early on, if build_index_while_first_file_is_read()
    // below asks: looked up meanwhile, from here on.
    std::future<bool> first_file_names_electron_activation;
    if (next_spectra.valid() && ions_by_activation_)
    {
      first_file_names_electron_activation = std::async(std::launch::async | std::launch::deferred,
        [&in_spectra_files]() { return mzMLHeadNamesElectronActivation(in_spectra_files[0]); });
    }
    // A database read here is read while the first file is: neither registers meta value names (see above).
    if (in_db != nullptr) { loadFASTA(*in_db, fasta_db_read); }

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
      sh.recalibrate = calibration_enabled_ ? calibration_recalibrate_ : "none";
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
    // A database read from *in_db is not used after it is augmented (at most once, here or in prepare_context()
    // below): its entries are taken over instead of being copied.
    const auto augment_db = [&]()
    {
      return in_db == nullptr ? buildDecoyAugmentedDB_(fasta_db, strategy) : buildDecoyAugmentedDB_(std::move(fasta_db_read), strategy);
    };
    std::vector<FASTAFile::FASTAEntry> full_db;
    bool use_chunked = false;
    if (database_chunk_size_ > 0)
    {
      full_db = augment_db();
      use_chunked = (full_db.size() > database_chunk_size_);
    }

    if (use_chunked)
    {
      // ================================================================
      // Chunk-major multi-file path (MSFragger-style): build each chunk's
      // FragmentIndex ONCE and score ALL files against it before moving
      // to the next chunk. This gives C FI builds instead of N×C.
      // No file is searched on its own here, so there is no "Searching <file>" progress step per file as in the
      // non-chunked path below; the steps of the phases and chunks are reported as before.
      // ================================================================
      // full_db already built above.
      if (rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static"))
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode " << reuse_mode_ << " needs an unchunked search: every file is "
                        << "rescored by Percolator trained on it ('per_file')." << std::endl;
      }
      const bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
      const bool open_search_mode = isOpenSearchMode_();

      // Shared report stats: chunk-major builds C indices (one per chunk),
      // reused across all files — count/time them once, not per file.
      mfres.shared.chunked = true;
      // Count by the RESOLVED decoy marker (prefix OR suffix), not the hardcoded decoy_prefix_: otherwise reused
      // external/suffix decoys (decoy_mode "external") would be miscounted as targets. have_decoys is false for
      // target-only (ignore), where decoys are stripped and decoy_string is empty.
      mfres.shared.db_decoy_proteins = strategy.have_decoys ? countDecoyEntries_(full_db, strategy.decoy_string, strategy.is_prefix) : 0;
      mfres.shared.db_target_proteins = full_db.size() - mfres.shared.db_decoy_proteins;

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
        all_spectra[i] = next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra_locations[i], readerThreads(false), nullptr, 0, chimericSearch_());
        preprocessSearchedSpectra_(all_spectra[i], in_spectra_files[i], all_ion_evidence[i],
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
      std::vector<CalibrationResult_> per_file_result(in_spectra_files.size()); // (the recalibration, recorded per file below)

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
          CalibrationResult_& cal = per_file_result[i];
          cal = runCalibrationPass_(all_spectra[i], cal_fi, cal_db, query_raw_spectrum_ ? &all_query_spectra[i] : nullptr);
          if (cal.success)
          {
            applyRecalibration_(cal, all_spectra[i], query_raw_spectrum_ ? &all_query_spectra[i] : nullptr,
                                retain_evidence ? &all_evidence_spectra[i] : nullptr);
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
        per_file_pool_stats[i].resize(all_spectra[i].size());
      }
      if (chimeric_enabled_)
      {
        OPENMS_LOG_WARN << "[ProSE] chimeric:enabled is not supported with database:chunk_size > 0 and is ignored." << std::endl;
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
        if (chunk_fi.isPeptideMassMode()) { mfres.shared.peptide_mass_mode = true; }

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
        const Size keep = std::max(searchTopHits_(), Size(2));
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
          searchTopHits_(), modifications_fixed_, modifications_variable_,
          peptide_missed_cleavages_,
          std::max(per_file_cal[i].effective_precursor_tol_lower,
                   per_file_cal[i].effective_precursor_tol_upper),
          per_file_cal[i].effective_fragment_tol,
          precursor_mass_tolerance_unit_, fragment_mass_tolerance_unit_,
          precursor_min_charge_, precursor_max_charge_, enzyme_, "",
          retain_evidence ? &all_evidence_spectra[i] : nullptr, nullptr, false,
          per_file_result[i].appliedFragmentRecalibration());
        recordRecalibration_(per_file_result[i], result.protein_ids);
        if (retain_evidence) all_evidence_spectra[i].clear(true);

        // Per input file, while each file is still its own identification run. Moving this
        // after the merge below would leave one run to calibrate NCE and RT on, which are
        // per-run quantities; see the note on PeptDeepRescoring::annotate().
        annotatePeptDeepFeatures_(all_spectra[i], result.protein_ids, result.peptide_ids);

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
        // PSM-level FDR, after Percolator rescoring with 'rescore' (see applyPsmFDR_)
        reset_run_name_ = in_spectra;
        result.percolator = applyPsmFDR_(result.protein_ids, result.peptide_ids, has_decoys, result.stats,
                                         keepPreRescoringPsms_() ? &result.pre_fdr_peptide_ids : nullptr,
                                         &result.inference_only_peptide_ids);

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
            full_db = augment_db();
            endProgress();
          }
          ctx.db = std::move(full_db);
          ctx.decoy_string = strategy.decoy_string;
          ctx.decoy_is_prefix = strategy.is_prefix;
          ctx.have_decoys = strategy.have_decoys;
          startProgress(0, 1, "Building fragment index...");
          buildFragmentIndex_(ctx, fragmentIndexParameters_(electron_ions), searched_spectra);
          ctx.electron_ions = electron_ions;
          endProgress();
          reportDecoyParity_(ctx);
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
        mfres.shared.peptide_mass_mode = ctx.fragment_index.isPeptideMassMode();
        // (counted by the resolved decoy marker, as in the chunk-major path)
        mfres.shared.db_decoy_proteins = strategy.have_decoys ? countDecoyEntries_(ctx.db, strategy.decoy_string, strategy.is_prefix) : 0;
        mfres.shared.db_target_proteins = ctx.db.size() - mfres.shared.db_decoy_proteins;
      };

      mfres.per_file.reserve(in_spectra_files.size());
      // the consumer of the files' results (setFileResultConsumer()); waited for before the results change
      BackgroundConsumer file_results;

      // rescore_reuse joint/static: the searches of the first reuse_files files defer their rescoring and FDR:PSM until
      // the models are trained on them (finishDeferredRescoring_()); with 'static', the later files are rescored with
      // these models in their search.
      const bool reuse_across_files = rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static")
                                      && in_spectra_files.size() > 1 && !isOpenSearchMode_();
      if (rescore_ && (reuse_mode_ == "joint" || reuse_mode_ == "static") && in_spectra_files.size() > 1 && isOpenSearchMode_())
      {
        OPENMS_LOG_WARN << "[ProSE] rescore_reuse:mode " << reuse_mode_ << " needs a closed search: every file is rescored "
                        << "by Percolator trained on it ('per_file')." << std::endl;
      }
      const Size reuse_files = !reuse_across_files ? 0
                               : reuse_mode_ == "joint" ? in_spectra_files.size()
                               : std::min(reuse_static_files_, in_spectra_files.size());
      std::vector<Size> reuse_deferred;
      reuse_model_.reset();
      // the state of this search is not left behind, whichever way it ends
      const OnScopeExit reset_reuse([this]() { reuse_defer_ = false; reuse_model_.reset(); reuse_file_index_ = 0; });

      // While the first file is read: build the index, as search(file) does without c and z+1 ions
      // unless the spectra are read already or the file names an electron-based activation early on.
      // prepare_context() below adds them if the spectra need them after all.
      // With restrict_index (see above), the spectra taken for the index are searched.
      PeakMap first_spectra;
      bool first_spectra_taken = false;
      // Called for the first file, inside its progress step (the build overlaps with reading it), if its spectra are
      // not read yet.
      const auto build_index_while_first_file_is_read = [&]()
      {
        if (next_spectra.valid() && !next_spectra.ready()
            && !(ions_by_activation_ && first_file_names_electron_activation.get()))
        {
          std::function<const PeakMap*(Size)> searched_spectra;
          if (restrict_index)
          {
            const UInt64 spectra_bytes = File::fileSize(in_spectra_files[0]); // UInt64(-1) if unknown: no waiting
            searched_spectra = [&, spectra_bytes](Size peptides) -> const PeakMap*
            {
              if (!first_spectra_taken)
              {
                const bool wait = wait_for_first_spectra || (threads_started <= MAX_THREADS_WAITING_FOR_SPECTRA
                                  && spectra_bytes <= MAX_SPECTRA_BYTES_PER_PEPTIDE * peptides / threads_started);
                if (!wait && !next_spectra.ready())
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
      };

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");

        // One progress step per input file, identified by its label, so that a progress display can tell which file
        // is searched. It ends with this iteration, whichever way it is left: the end, the continue of a failed
        // search, an exception. The steps of the search of the file nest below it (depth 1) only with one thread:
        // the background mzML reader shares the static step depth of the ProgressLogger.
        if (reuse_across_files && i == reuse_files)
        {
          finishDeferredRescoring_(mfres, reuse_deferred); // 'static': the models for this file and the later ones
          reuse_deferred.clear();
        }
        reuse_defer_ = reuse_across_files && i < reuse_files;
        reuse_file_index_ = i;

        startProgress(0, 1, "Searching " + File::basename(in_spectra) + " (" + StringUtils::toStr(i + 1) + " of "
                                + StringUtils::toStr(in_spectra_files.size()) + ")...");
        const OnScopeExit end_file_progress([this]() { endProgress(); });

        OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                        << "] Searching " << in_spectra << std::endl;

        if (i == 0) { build_index_while_first_file_is_read(); }

        PeakMap spectra = first_spectra_taken ? std::move(first_spectra)
                          : next_spectra.valid() ? next_spectra.get() : loadMS2Spectra(in_spectra_locations[i], readerThreads(false), nullptr, 0, chimericSearch_());
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
        reset_run_name_ = in_spectra_files[i]; // names the RESET stem table of the run (reset:stem_table a directory)
        result.exit_code = search(spectra, ctx, result.protein_ids, result.peptide_ids);

        if (result.exit_code != ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_WARN << "[ProSE] Search failed for " << in_spectra
                          << " (exit code " << static_cast<int>(result.exit_code) << "). Continuing." << std::endl;
          result.stats = last_run_stats_;
          result.percolator = last_percolator_;
          result.stats.input_file = File::basename(in_spectra);
          mfres.per_file.push_back(std::move(result));
          continue;
        }

        result.stats = last_run_stats_;
        result.percolator = last_percolator_;
        result.stats.input_file = File::basename(in_spectra);
        if (keepPreRescoringPsms_()) { result.pre_fdr_peptide_ids = std::move(pre_fdr_peptide_ids_); }
        pre_fdr_peptide_ids_ = PeptideIdentificationList();
        result.inference_only_peptide_ids = std::move(inference_only_peptide_ids_);
        inference_only_peptide_ids_ = PeptideIdentificationList();

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
        // After the last file, its spectra are released by a helper thread, which the caller then writes the results
        // next to (0.5 s for the 250,000 spectra of a timsTOF run at 64 threads, most of it in the kernel returning their
        // memory); the holder joins it like the hits that the rescoring dropped (see released_hits_). Otherwise here.
        bool released_later = false;
#if !defined(OPENMS_WINDOWSPLATFORM) && defined(_OPENMP) // Windows: the holder may go under the loader lock (a static object of a DLL), where joining deadlocks
        if (i + 1 == in_spectra_files.size() && spectra.size() >= 1024 && omp_get_max_threads() > 1 && !helpersOversubscribe())
        {
          try
          {
            auto release = std::make_shared<SpectraRelease>();
            release->before = released_hits_;
            // with a quarter of the threads of this thread, as the dropped hits (the release of the database and the
            // output after this use parallel regions of their own)
            release->thread = std::thread([to_release = std::move(spectra), threads = std::max(1, omp_get_max_threads() / 4)]() mutable
            {
              omp_set_num_threads(threads);
              freeInParallel(to_release);
            });
            released_hits_ = std::move(release);
            released_later = true;
          }
          catch (...) {} // the spectra (if not moved into the thread's function, and released with it) are released here
        }
#endif
        if (!released_later) { freeInParallel(spectra); }
        // after the last file, the database the index was built from
        if (i + 1 == in_spectra_files.size()) { freeInParallel(ctx.db); }

        if (reuse_defer_) { reuse_deferred.push_back(mfres.per_file.size()); }
        mfres.per_file.push_back(std::move(result));
        // (per_file holds all files' places, reserved above: the consumer's result stays where it is)
        // (with pool_rescoring_, a file that Percolator did not rescore on its own may still change: rescoreFilesTogether_())
        if (file_result_consumer_ && !reuse_defer_ && i + 1 < in_spectra_files.size()
            && mfres.per_file.back().exit_code == ExitCodes::EXECUTION_OK
            && (!pool_rescoring_ || mfres.per_file.back().percolator.applied))
        {
          file_results.start(file_result_consumer_, i, mfres.per_file.back());
        }
      }
      file_results.finish();
      finishDeferredRescoring_(mfres, reuse_deferred); // 'joint': all files (nothing without rescore_reuse)
    }

    // 'rescore': the files that Percolator could not rescore on their own are rescored with all others, with one model.
    if (pool_rescoring_) { rescoreFilesTogether_(mfres, in_spectra_files); }

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
      // Copying is mandatory here: per_file is returned to the caller and must
      // stay intact, so nothing may be moved out of it.
      poolSearchedFiles_(mfres.per_file, /*move=*/false, in_spectra_files,
                         mfres.aggregate.protein_ids, mfres.aggregate.peptide_ids);
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
    // load FASTA once (while the first spectrum file is read in the background)
    MultiFileSearchResult mfres = searchWithModificationAnalysis_(
      in_spectra_files, {}, &in_db, output_base_names, aggregate_base_name, build_pooled_aggregate);

    mfres.shared.database_file = in_db;

    // Patch database file path into each per-file and aggregate ProteinIdentification.
    for (auto& pf : mfres.per_file)
    {
      if (searched_(pf)) { pf.protein_ids[0].getSearchParameters().db = in_db; }
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

    // A single input is a complete experiment: protein FDR as in search(file).
    SearchResult& result = mfres.per_file[0];
    if (fdr_protein_ > 0.0 && mfres.have_decoys && result.exit_code == ExitCodes::EXECUTION_OK
        && !result.protein_ids.empty() && !result.peptide_ids.empty())
    {
      applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids, mfres.decoy_string, mfres.decoy_is_prefix, fdr_protein_,
                                 true, result.inference_only_peptide_ids);
      updateFinalStats(result.stats, result.peptide_ids, enzyme_, result.stats.fdr_applied);
    }
    result.inference_only_peptide_ids = PeptideIdentificationList(); // used by the protein FDR only
    return std::move(result);
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

  bool ProSEAlgorithm::useSearchEngineScore_(PeptideIdentificationList& ids)
  {
    // The name of the score of each PSM: none for those that have it as main score already
    const std::string main_score = SEARCH_ENGINE_SCORE, meta_score = main_score + "_score";
    std::vector<const std::string*> source(ids.size(), nullptr);
    for (Size i = 0; i < ids.size(); ++i)
    {
      const PeptideIdentification& id = ids[i];
      if (id.getHits().empty() || id.getScoreType() == main_score) { continue; }
      const PeptideHit& hit = id.getHits()[0];
      if (hit.metaValueExists(main_score)) { source[i] = &main_score; }
      else if (hit.metaValueExists(meta_score)) { source[i] = &meta_score; }
      else { return false; }
    }
    // one switcher per meta value name
    const auto switcher = [&main_score](const std::string& name)
    {
      IDScoreSwitcherAlgorithm s;
      Param p = s.getParameters();
      p.setValue("new_score", name);
      p.setValue("new_score_type", main_score);
      p.setValue("new_score_orientation", "higher_better");
      s.setParameters(p);
      return s;
    };
    IDScoreSwitcherAlgorithm from_main = switcher(main_score), from_meta = switcher(meta_score);
    Size counter = 0;
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (source[i] != nullptr) { (source[i] == &main_score ? from_main : from_meta).switchScores(ids[i], counter); }
    }
    return true;
  }

  Size ProSEAlgorithm::countDecoyEntries_(const std::vector<FASTAFile::FASTAEntry>& db, const std::string& marker, bool is_prefix)
  {
    return static_cast<Size>(std::count_if(db.begin(), db.end(), [&](const FASTAFile::FASTAEntry& e)
                                           { return accessionHasDecoyMarker_(e.identifier, marker, is_prefix); }));
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
    // The top hits of the spectra are read in parallel (the target/decoy value by name index: looking the name up takes
    // a lock on the registry, which the threads would contend for; by name if it is not registered), then their scores
    // collected and counted in spectrum order, as before (also where an error stops it).
    enum TopHit : char { NO_HIT, OTHER, TARGET, DECOY };
    const UInt td_index = MetaInfoInterface::metaRegistry().getIndex(Constants::UserParam::TARGET_DECOY);
    const SignedSize num_ids = static_cast<SignedSize>(peptide_ids.size());
    std::vector<TopHit> top_hits(peptide_ids.size(), NO_HIT);
    std::vector<double> top_scores(peptide_ids.size());
    std::vector<std::exception_ptr> errors(peptide_ids.size());
#pragma omp parallel for schedule(static, 1024)
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      try
      {
        const PeptideIdentification& pid = peptide_ids[i];
        if (pid.getHits().empty()) { continue; }
        const PeptideHit& top = pid.getHits().front();
        top_hits[i] = OTHER;
        top_scores[i] = top.getScore();

        const bool by_name = td_index == UInt(-1);
        if (by_name ? !top.metaValueExists(Constants::UserParam::TARGET_DECOY) : !top.metaValueExists(td_index)) { continue; }
        const std::string td = (by_name ? top.getMetaValue(Constants::UserParam::TARGET_DECOY) : top.getMetaValue(td_index)).toString();
        // "target+decoy" counts as target (OpenMS FDR semantics).
        if (td == "decoy") { top_hits[i] = DECOY; }
        else if (!td.empty()) { top_hits[i] = TARGET; }
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    std::vector<double> hyperscores;
    for (SignedSize i = 0; i < num_ids; ++i)
    {
      if (top_hits[i] == NO_HIT) { continue; }
      hyperscores.push_back(top_scores[i]);
      if (errors[i]) { std::rethrow_exception(errors[i]); }
      if (top_hits[i] == DECOY) { ++stats.decoy_psms; }
      else if (top_hits[i] == TARGET) { ++stats.target_psms; }
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

  namespace
  {
    // What collectRunStatistics_() and updateFinalStats() count over the top hits
    struct TopHitSummary
    {
      Size matched = 0, n_target = 0, n_decoy = 0;
      std::map<Int, Size> charges;
      std::map<Size, Size> missed_cleavages;
      Size unique_peptides = 0, unique_proteins = 0; // distinct sequence strings and protein accessions
      std::vector<double> precursor_errors, fragment_errors; // top hits in input order (with_errors only)
      Size fragment_errors_psms = 0; // the PSMs with any fragment error (with_errors only)
    };

    // One pass over the top hits, in parallel: each thread takes a contiguous range of the identifications, and the
    // ranges are merged in order, so the error lists come out in input order as from a single pass. The regular
    // expression of the missed-cleavage count and the sequence strings made this a serial step of 20 - 40 ms per file.
    // The distinct peptides and proteins are counted in parallel, too: each range notes its sequence strings (in one
    // buffer) and accessions in buckets by hash, and the buckets are counted one per thread. Merging sets of strings
    // per range and freeing them afterwards took 0.6 s on one thread for the 250,000 top hits of a timsTOF run.
    //
    // The mass errors of the tolerance estimate (with_errors), of the top hits of the target PSMs (a decoy is a wrong
    // match: its errors are noise), signed (observed minus theoretical, ppm): of the precursor, corrected by the isotope
    // error of the hit (as the search matched it; @p isotope_observed_minus_theoretical: the convention of the reported
    // 'isotope_error'), without the hits with a delta mass of an open search; and of the annotated b and y fragment peaks.
    TopHitSummary summarizeTopHits(const PeptideIdentificationList& peptide_ids, const std::string& enzyme, bool with_errors,
                                   bool isotope_observed_minus_theoretical = true)
    {
      const DigestionEnzymeProtein* protease = ProteaseDB::getInstance()->getEnzyme(enzyme);
      const MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
      const UInt target_decoy_index = registry.getIndex(Constants::UserParam::TARGET_DECOY); // UInt(-1): exists nowhere
      const UInt delta_mass_index = registry.getIndex("DeltaMass");
      const UInt isotope_error_index = registry.getIndex(Constants::UserParam::ISOTOPE_ERROR);
      static const double water_mass = EmpiricalFormula("H2O").getMonoWeight();
      const SignedSize n = static_cast<SignedSize>(peptide_ids.size());
      // Declared without OpenMP, too: with the -fopenmp-simd fallback the compiler still parses the clauses below.
      int num_parts = 1;
#ifdef _OPENMP
      num_parts = static_cast<int>(std::clamp<SignedSize>(n / 512, 1, omp_get_max_threads()));
#endif
      const Size num_buckets = static_cast<Size>(num_parts);
      struct Part
      {
        TopHitSummary summary;
        std::string sequences;                                             // the sequence strings, one after the other
        std::vector<std::vector<std::pair<Size, Size>>> peptide_buckets;   // (offset, length) in sequences, by hash
        std::vector<std::vector<const std::string*>> protein_buckets;      // accessions of the peptide evidences, by hash
      };
      std::vector<Part> parts(num_parts);
      std::vector<std::exception_ptr> errors(num_parts); // an exception must not leave the parallel region
#pragma omp parallel for schedule(static, 1) num_threads(num_parts)
      for (int part = 0; part < num_parts; ++part)
      try
      {
        Part& p = parts[part];
        TopHitSummary& s = p.summary;
        p.peptide_buckets.resize(num_buckets);
        p.protein_buckets.resize(num_buckets);
        const std::hash<std::string_view> hash;
        EnzymaticDigestion digestor;
        digestor.setEnzyme(protease);
        std::unordered_map<std::string, Size> missed_cleavages_of; // by unmodified sequence
        std::vector<double> residue_masses; // of a hit, cumulative (prefix sums, for the b and y ions)
        for (SignedSize i = n * part / num_parts; i < n * (part + 1) / num_parts; ++i)
        {
          const PeptideIdentification& pid = peptide_ids[i];
          if (pid.getHits().empty()) { continue; }
          ++s.matched;
          const PeptideHit& top = pid.getHits().front();
          const std::string sequence = top.getSequence().toString();
          p.peptide_buckets[hash(sequence) % num_buckets].emplace_back(p.sequences.size(), sequence.size());
          p.sequences += sequence;
          for (const auto& ev : top.getPeptideEvidences())
          {
            const std::string& accession = ev.getProteinAccession();
            p.protein_buckets[hash(accession) % num_buckets].push_back(&accession);
          }
          if (top.getCharge() > 0) { ++s.charges[top.getCharge()]; }
          // the target/decoy type as PeptideHit::getTargetDecoyType() (by the index of the meta value; case-insensitive;
          // another value is neither, where getTargetDecoyType() throws)
          std::string td;
          if (top.metaValueExists(target_decoy_index)) { td = StringUtils::toLowered(top.getMetaValue(target_decoy_index).toString()); }
          const bool is_target = td == "target" || td == "target+decoy";
          if (with_errors)
          {
            // Precursor error: skip hits with unresolved charge (0) - getMZ() would throw - and, in an open search, those
            // with a delta mass (a modification the sequence lacks, not a mass error; as OpenSearchModificationAnalysis
            // tells unmodified PSMs).
            const bool shifted = top.metaValueExists(delta_mass_index)
                                 && !OpenSearchModificationAnalysis::isUnmodified(static_cast<double>(top.getMetaValue(delta_mass_index)));
            if (is_target && !shifted && top.getCharge() > 0)
            {
              // the isotope offset the search matched at (theoretical minus observed)
              int offset = 0;
              if (top.metaValueExists(isotope_error_index))
              {
                const int reported = static_cast<int>(top.getMetaValue(isotope_error_index));
                offset = isotope_observed_minus_theoretical ? -reported : reported;
              }
              const double corrected_mz = monoisotopicPrecursorMZ(pid.getMZ(), offset, static_cast<double>(top.getCharge()));
              s.precursor_errors.push_back(Math::getPPM(corrected_mz, top.getSequence().getMZ(top.getCharge())));
            }
            // Fragment errors: of the annotated b and y peaks ("b3+", "y12++"; not those with a neutral loss or of other
            // ions), against their m/z from the hit's sequence: (prefix or suffix residue masses + terminal modification
            // [+ H2O] + charge protons) / charge, the masses TheoreticalSpectrumGenerator gives the annotated ions
            if (is_target && !top.getPeakAnnotations().empty())
            {
              const AASequence& sequence = top.getSequence();
              const Size length = sequence.size();
              residue_masses.assign(length + 1, 0.0);
              for (Size k = 0; k < length; ++k)
              {
                residue_masses[k + 1] = residue_masses[k] + sequence[k].getMonoWeight(Residue::Internal);
              }
              const double n_term = sequence.hasNTerminalModification() ? sequence.getNTerminalModification()->getDiffMonoMass() : 0.0;
              const double c_term = sequence.hasCTerminalModification() ? sequence.getCTerminalModification()->getDiffMonoMass() : 0.0;
              bool any = false;
              for (const PeptideHit::PeakAnnotation& pa : top.getPeakAnnotations())
              {
                const std::string& name = pa.annotation;
                bool prefix = false;
                Size ordinal = 0;
                if (name.empty() || (name[0] != 'b' && name[0] != 'y') || pa.charge < 1
                    || !FragmentIonLikelihoodModel::parseIonName(name, prefix, ordinal) || ordinal >= length
                    || name.find_first_not_of('+', name.find_first_not_of("0123456789", 1)) != std::string::npos)
                {
                  continue; // another ion, or one with a neutral loss ("y3-H2O1+")
                }
                const double neutral = prefix ? residue_masses[ordinal] + n_term
                                              : residue_masses[length] - residue_masses[length - ordinal] + c_term + water_mass;
                const double theoretical_mz = (neutral + pa.charge * Constants::PROTON_MASS_U) / pa.charge;
                s.fragment_errors.push_back(Math::getPPM(pa.mz, theoretical_mz));
                any = true;
              }
              if (any) { ++s.fragment_errors_psms; }
            }
          }
          // counted once per sequence of the range: the regular expression takes a block from boost's process-wide
          // memory cache per search, which the threads contended for (93 ms at 16 and 64 threads for a timsTOF run)
          std::string unmodified = top.getSequence().toUnmodifiedString();
          auto cleavages = missed_cleavages_of.find(unmodified);
          if (cleavages == missed_cleavages_of.end())
          {
            const Size count = digestor.countInternalCleavageSites(unmodified);
            cleavages = missed_cleavages_of.emplace(std::move(unmodified), count).first;
          }
          ++s.missed_cleavages[cleavages->second];
          if (td == "decoy") { ++s.n_decoy; }
          else if (is_target) { ++s.n_target; }
        }
      }
      catch (...)
      {
        errors[part] = std::current_exception();
      }
      for (const std::exception_ptr& error : errors)
      {
        if (error) { std::rethrow_exception(error); }
      }

      // The distinct strings of each bucket (equal strings have equal hashes, so they share a bucket), and the buffers of
      // each range freed by a thread of the team
      std::vector<Size> unique_peptides(num_buckets, 0), unique_proteins(num_buckets, 0);
      std::vector<std::exception_ptr> bucket_errors(num_buckets);
#pragma omp parallel num_threads(num_parts)
      {
        std::vector<std::string_view> strings;
#pragma omp for schedule(dynamic, 1)
        for (SignedSize sb = 0; sb < static_cast<SignedSize>(num_buckets); ++sb)
        try
        {
          const Size b = static_cast<Size>(sb);
          strings.clear();
          for (const Part& p : parts)
          {
            for (const auto& [offset, length] : p.peptide_buckets[b]) { strings.emplace_back(p.sequences.data() + offset, length); }
          }
          std::sort(strings.begin(), strings.end());
          unique_peptides[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
          strings.clear();
          for (const Part& p : parts)
          {
            for (const std::string* accession : p.protein_buckets[b]) { strings.emplace_back(*accession); }
          }
          std::sort(strings.begin(), strings.end());
          unique_proteins[b] = static_cast<Size>(std::unique(strings.begin(), strings.end()) - strings.begin());
        }
        catch (...)
        {
          bucket_errors[sb] = std::current_exception();
        }
#pragma omp for schedule(static, 1)
        for (int part = 0; part < num_parts; ++part)
        {
          std::string().swap(parts[part].sequences);
          decltype(Part::peptide_buckets)().swap(parts[part].peptide_buckets);
          decltype(Part::protein_buckets)().swap(parts[part].protein_buckets);
        }
      }

      for (const std::exception_ptr& error : bucket_errors)
      {
        if (error) { std::rethrow_exception(error); }
      }
      Size num_precursor_errors = 0, num_fragment_errors = 0;
      for (const Part& p : parts)
      {
        num_precursor_errors += p.summary.precursor_errors.size();
        num_fragment_errors += p.summary.fragment_errors.size();
      }
      TopHitSummary all = std::move(parts.front().summary);
      all.unique_peptides = std::accumulate(unique_peptides.begin(), unique_peptides.end(), Size(0));
      all.unique_proteins = std::accumulate(unique_proteins.begin(), unique_proteins.end(), Size(0));
      all.precursor_errors.reserve(num_precursor_errors);
      all.fragment_errors.reserve(num_fragment_errors);
      for (int part = 1; part < num_parts; ++part)
      {
        TopHitSummary& s = parts[part].summary;
        all.matched += s.matched;
        all.n_target += s.n_target;
        all.n_decoy += s.n_decoy;
        all.fragment_errors_psms += s.fragment_errors_psms;
        for (const auto& [charge, count] : s.charges) { all.charges[charge] += count; }
        for (const auto& [cleavages, count] : s.missed_cleavages) { all.missed_cleavages[cleavages] += count; }
        all.precursor_errors.insert(all.precursor_errors.end(), s.precursor_errors.begin(), s.precursor_errors.end());
        all.fragment_errors.insert(all.fragment_errors.end(), s.fragment_errors.begin(), s.fragment_errors.end());
      }
      return all;
    }
    // Math::median() and Math::MAD() of @p values (reordered), by selection instead of sorts: the same order statistics,
    // so the same values (the sorts of the 250,000 errors of a timsTOF run took 20-30 ms each). Unless a value is not
    // finite (a NaN is not ordered; an infinite median gives NaN deviations) or a negative zero (equal to 0, but written
    // differently), where the result may depend on the algorithm or on the order the sort left: then Math itself.
    std::pair<double, double> medianAndMADBySelection(std::vector<double>& values)
    {
      const auto ambiguous = [](double v) { return !std::isfinite(v) || (v == 0.0 && std::signbit(v)); };
      if (values.empty() || std::any_of(values.begin(), values.end(), ambiguous))
      {
        const double median = Math::median(values.begin(), values.end());
        return {median, Math::MAD(values.begin(), values.end(), median)};
      }
      // the median of finite values without a negative zero; their deviations are not negative zeros and not NaN
      // Large lists (the b/y fragment errors of a timsTOF run: millions) with the OpenMP threads: the selection by
      // libstdc++'s parallel nth_element (the same k-th smallest value), the deviations in parallel.
      constexpr Size parallel_from = 65536;
      const auto median_of = [](std::vector<double>& v)
      {
        const auto upper = v.begin() + static_cast<SignedSize>(v.size() / 2);
#ifdef PROSE_GNU_PARALLEL
        if (v.size() >= parallel_from && omp_get_max_threads() > 1) { __gnu_parallel::nth_element(v.begin(), upper, v.end()); }
        else
#endif
        { std::nth_element(v.begin(), upper, v.end()); }
        if (v.size() % 2 == 1) { return *upper; }
        return (*std::max_element(v.begin(), upper) + *upper) / 2.0;
      };
      const double median = median_of(values);
      std::vector<double> diffs(values.size());
#pragma omp parallel for schedule(static) if (values.size() >= parallel_from)
      for (SignedSize i = 0; i < static_cast<SignedSize>(values.size()); ++i) { diffs[i] = fabs(values[i] - median); }
      return {median, median_of(diffs)};
    }

  // =====================================================================
  // Recompute result-level stats from a FINAL (post-rescoring/post-FDR) PSM list.
  // =====================================================================
    // Spectra with at least one PSM. The PSMs of co-isolated precursors (chimeric:enabled) belong to the spectrum
    // (scan_index) of the selected precursor, which counts once.
    Size countMatchedSpectra(const PeptideIdentificationList& peptide_ids)
    {
      Size matched = 0;
      std::set<unsigned int> selected_scans, coisolated_scans;
      for (const PeptideIdentification& pid : peptide_ids)
      {
        if (pid.getHits().empty()) continue;
        const PeptideHit& top = pid.getHits().front();
        const bool coisolated = top.metaValueExists(ProSEAlgorithm::CHIMERIC_RANK)
                                && static_cast<int>(top.getMetaValue(ProSEAlgorithm::CHIMERIC_RANK)) == 2;
        if (!pid.metaValueExists("scan_index"))
        {
          if (!coisolated) ++matched;
          continue;
        }
        const unsigned int scan = static_cast<unsigned int>(pid.getMetaValue("scan_index"));
        (coisolated ? coisolated_scans : selected_scans).insert(scan);
        if (!coisolated) ++matched;
      }
      for (unsigned int scan : coisolated_scans)
      {
        if (selected_scans.count(scan) == 0) ++matched;
      }
      return matched;
    }

    // countMatchedSpectra(), given @p with_hits, the number of identifications with hits: the same unless a top hit is
    // of a co-isolated precursor (chimeric:enabled), which is checked first (in parallel: no set of scans per search)
    Size matchedSpectra(const PeptideIdentificationList& peptide_ids, Size with_hits)
    {
      const UInt rank_index = MetaInfoInterface::metaRegistry().getIndex(ProSEAlgorithm::CHIMERIC_RANK);
      if (rank_index == UInt(-1)) { return with_hits; } // (a name that was never registered is on no hit)
      bool any_coisolated = false;
#pragma omp parallel for schedule(static) reduction(||: any_coisolated)
      for (SignedSize i = 0; i < static_cast<SignedSize>(peptide_ids.size()); ++i)
      {
        const auto& hits = peptide_ids[i].getHits();
        if (!hits.empty() && hits.front().metaValueExists(rank_index)
            && static_cast<int>(hits.front().getMetaValue(rank_index)) == 2) { any_coisolated = true; }
      }
      return any_coisolated ? countMatchedSpectra(peptide_ids) : with_hits;
    }
  }

  void ProSEAlgorithm::updateFinalStats(RunStatistics& stats,
                                        const PeptideIdentificationList& peptide_ids,
                                        const std::string& enzyme,
                                        bool fdr_applied)
  {
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    const TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme, false);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    stats.charge_histogram = summary.charges;
    stats.missed_cleavage_histogram = summary.missed_cleavages;
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    stats.fdr_applied = fdr_applied;
    stats.achieved_psm_fdr = fdr_applied ? maxRetainedScore_(peptide_ids) : -1.0;
  }

  void ProSEAlgorithm::releaseInParallel(SearchResult& result)
  {
    for (auto& protein_id : result.protein_ids) { freeInParallel(protein_id.getHits()); }
    std::vector<ProteinIdentification>().swap(result.protein_ids);
    freeInParallel(result.peptide_ids);
    freeInParallel(result.pre_fdr_peptide_ids);
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
  ProSEAlgorithm::SpectraCounts_ ProSEAlgorithm::countSpectra_(const PeakMap& spectra) const
  {
    // in parallel: one pass over the spectra of a timsTOF run (600,000) took 0.13 s on one thread (memory bound)
    const bool search_outside = param_.getValue("precursor:recorded_charge_outside_range").toString() == "search";
    const int min_charge = static_cast<int>(precursor_min_charge_), max_charge = static_cast<int>(precursor_max_charge_);
    Size ms2 = 0, profile = 0, out_of_range = 0;
#pragma omp parallel reduction(+: ms2, profile, out_of_range) if (spectra.size() >= 16384)
    {
      std::vector<uint16_t> charges;
#pragma omp for schedule(static)
      for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
      {
        const MSSpectrum& s = spectra[i];
        if (s.getMSLevel() != 2) continue;
        ++ms2;
        if (s.getType() == SpectrumSettings::SpectrumType::PROFILE) ++profile;
        if (s.getPrecursors().empty()) continue;
        FragmentIndex::precursorCharges(s.getPrecursors()[0], min_charge, max_charge, charges, search_outside);
        if (charges.empty()) ++out_of_range;
      }
    }
    SpectraCounts_ counts;
    counts.ms2 = ms2;
    counts.profile = profile;
    counts.charge_out_of_range = out_of_range;
    return counts;
  }

  void ProSEAlgorithm::preprocessSearchedSpectra_(PeakMap& spectra, const std::string& input,
                                                   FragmentIonLikelihoodModel::PeakLists& ion_evidence,
                                                   PeakMap* evidence_spectra, PeakMap* query_spectra) const
  {
    const SpectraCounts_ counts = countSpectra_(spectra);
    const std::string of = input.empty() ? std::string() : " of '" + input + "'";
    if (counts.ms2 == 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << (input.empty() ? std::string("the input") : "'" + input + "'")
                      << " holds no MS2 spectra (e.g. an MS1-only file): "
                      << "nothing is searched and no peptides are identified." << std::endl;
    }
    if (counts.profile > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.profile << " of " << counts.ms2 << " MS2 spectra" << of
                      << " are recorded as profile spectra. ProSE expects centroided spectra and searches profile "
                      << "spectra as they are, which identifies few peptides, if any. Centroid them first (e.g. with "
                      << "PeakPickerHiRes, or peak picking when converting the raw files)." << std::endl;
    }
    if (counts.charge_out_of_range > 0)
    {
      OPENMS_LOG_WARN << "[ProSE] Warning: " << counts.charge_out_of_range << " of " << counts.ms2
                      << " MS2 spectra" << of << " have a precursor charge outside of precursor:min_charge.."
                      << "max_charge (" << precursor_min_charge_ << ".." << precursor_max_charge_
                      << ") and are not searched." << std::endl;
    }
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ == "ppm", deisotope_requested_,
                       peaks_keep_n_, peaks_window_top_, peaks_window_type_, deisotoping_,
                       self_trained_ion_priors_ ? &ion_evidence : nullptr, ion_prior_scored_peaks_, evidence_spectra,
                       query_spectra);
  }

  void ProSEAlgorithm::collectRunStatistics_(
      const PeakMap& spectra,
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids,
      RunStatistics& stats) const
  {
    const SpectraCounts_ counts = countSpectra_(spectra);
    stats.ms2_spectra = counts.ms2;
    stats.profile_spectra = counts.profile;
    stats.charge_out_of_range_spectra = counts.charge_out_of_range;
    // Count spectra that actually retained a hit. ProSE pushes one PeptideIdentification per
    // searched spectrum (incl. empty ones for non-matches), and not every caller strips empties
    // before stats are collected, so peptide_ids.size() would over-count matched spectra / ID rate.
    // Target/decoy of the FINAL (post-FDR) hits, so SearchResult::stats stays consistent
    // for library callers that read it without calling updateFinalStats().
    // Mass errors for the tolerance estimate: see summarizeTopHits() (signed, isotope-corrected, of the target PSMs).
    TopHitSummary summary = summarizeTopHits(peptide_ids, enzyme_, true, isotope_error_observed_minus_theoretical_);
    stats.matched_spectra = matchedSpectra(peptide_ids, summary.matched);
    for (const auto& [charge, count] : summary.charges) { stats.charge_histogram[charge] += count; }
    for (const auto& [cleavages, count] : summary.missed_cleavages) { stats.missed_cleavage_histogram[cleavages] += count; }
    stats.unique_peptides = summary.unique_peptides;
    stats.unique_proteins = summary.unique_proteins;
    stats.target_psms = summary.n_target;
    stats.decoy_psms = summary.n_decoy;
    vector<double>& precursor_errors = summary.precursor_errors;
    vector<double>& fragment_errors = summary.fragment_errors;
    const Size fragment_errors_psms = summary.fragment_errors_psms;

    // -- Per-run tolerance estimate: the median of the signed errors (the bias) and their MAD (the spread; unscaled).
    //    The recommended (symmetric) tolerance covers the bias plus three standard deviations of a normal distribution,
    //    estimated robustly as 1.4826 * MAD: ceil(|median| + 3 * 1.4826 * MAD). Of the fragments it is an estimate
    //    within the searched tolerance: peaks beyond it were not matched. --
    const Size min_psms_for_estimation = 10;
    const auto recommended_tolerance = [](double median, double mad) { return std::ceil(std::fabs(median) + 3.0 * 1.4826 * mad); };
    if (precursor_errors.size() >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(precursor_errors); // = Math::median() / Math::MAD()
      stats.prec_err_median = med;
      stats.prec_err_mad = mad;
      stats.prec_err_recommended = recommended_tolerance(med, mad);
      stats.prec_tol_valid = true;
    }
    if (fragment_errors_psms >= min_psms_for_estimation)
    {
      const auto [med, mad] = medianAndMADBySelection(fragment_errors);
      stats.frag_err_median = med;
      stats.frag_err_mad = mad;
      stats.frag_err_recommended = recommended_tolerance(med, mad);
      stats.frag_tol_valid = true;
    }

    // -- calibration: the records of recordRecalibration_() (measured errors; recalibration), for the run summary --
    stats.recalibration.clear();
    if (!protein_ids.empty() && protein_ids[0].getSearchParameters().metaValueExists("calibration:precursor_error_offset"))
    {
      const ProteinIdentification::SearchParameters& sp = protein_ids[0].getSearchParameters();
      std::vector<std::string> keys;
      sp.getKeys(keys);
      for (const std::string& key : keys)
      {
        if (key.compare(0, 12, "calibration:") == 0) { stats.recalibration.emplace_back(key.substr(12), sp.getMetaValue(key)); }
      }
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
    // The residues are those of the sequence: fromResidues() takes them as they are, where operator+=(const Residue*)
    // looks every residue up in the ResidueDB under a critical section (the threads of annotateIonPriors_() spent
    // a third of their time in that section in stack samples of a 27-run search).
    std::vector<const Residue*> residues;
    residues.reserve(n);
    for (Size i = n - 1; i-- > 0;) { residues.push_back(&sequence[i]); }
    residues.push_back(&sequence[n - 1]);
    AASequence reversed = AASequence::fromResidues(std::move(residues));
    if (sequence.hasNTerminalModification()) { reversed.setNTerminalModification(sequence.getNTerminalModification()); }
    if (sequence.hasCTerminalModification()) { reversed.setCTerminalModification(sequence.getCTerminalModification()); }
    return reversed;
  }

  namespace
  {
  // =====================================================================
  // Helper: run calibration pass on a subset of spectra
  // =====================================================================
    /// The signed errors (ppm, or Da) of the peaks of @p spec closest to the theoretical peaks of @p theo within
    /// @p tolerance (unmatched theoretical peaks contribute nothing); both spectra sorted by m/z. With @p mzs, the
    /// theoretical m/z of each matched peak, aligned with @p errors.
    void matchedFragmentErrors(const MSSpectrum& spec, const PeakSpectrum& theo, double tolerance, bool ppm, std::vector<float>& errors,
                               std::vector<float>* mzs = nullptr)
    {
      for (const Peak1D& t : theo)
      {
        const Int i = spec.findNearest(t.getMZ(), ppm ? t.getMZ() * tolerance * 1e-6 : tolerance);
        if (i < 0) continue;
        const double e = spec[i].getMZ() - t.getMZ();
        errors.push_back(static_cast<float>(ppm ? e / t.getMZ() * 1e6 : e));
        if (mzs != nullptr) mzs->push_back(static_cast<float>(t.getMZ()));
      }
    }
  }

  ProSEAlgorithm::CalibrationResult_
  ProSEAlgorithm::runCalibrationPass_(
      PeakMap& spectra,
      FragmentIndex& fragment_index,
      const std::vector<FASTAFile::FASTAEntry>& db,
      const PeakMap* query_spectra) const
  {
    CalibrationResult_ result;

    const bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    const bool precursor_ppm = (precursor_mass_tolerance_unit_ == "ppm");

    // The pass searches with windows wide enough to measure the error (calibrationPassWindow_(); the windows of the main
    // search are derived from the measurement below and never exceed the configured ones). Query-time parameters of the
    // index only, restored on exit.
    const double pass_lower = calibrationPassWindow_(precursor_mass_tolerance_lower_, precursor_ppm);
    const double pass_upper = calibrationPassWindow_(precursor_mass_tolerance_upper_, precursor_ppm);
    const double pass_fragment = calibrationPassWindow_(fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
    const Param fi_params_saved = fragment_index.getParameters();
    {
      Param p = fi_params_saved;
      p.setValue("precursor:mass_tolerance_lower", pass_lower);
      p.setValue("precursor:mass_tolerance_upper", pass_upper);
      p.setValue("fragment:mass_tolerance", pass_fragment);
      fragment_index.setParameters(p);
    }
    const OnScopeExit restore_fi([&fragment_index, &fi_params_saved]()
    {
      try { fragment_index.setParameters(fi_params_saved); }
      catch (...) {} // only after another exception: the caller restores the parameters of its index
    });

    // The subset: the spectra of the highest TIC (the best-measured ones), taken per tenth of the run (the spectra are in
    // acquisition order) so that a mass offset that changes over the run is represented in the sample
    vector<pair<double, Size>> tic_index(spectra.size());
#pragma omp parallel for schedule(static) if (spectra.size() >= 4096)
    for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
    {
      double tic = 0;
      for (const auto& p : spectra[i]) { tic += p.getIntensity(); }
      tic_index[i] = {tic, static_cast<Size>(i)};
    }
    vector<Size> subset;
    {
      const Size blocks = std::min<Size>(10, spectra.size());
      for (Size b = 0; b < blocks; ++b)
      {
        const Size begin = spectra.size() * b / blocks, end = spectra.size() * (b + 1) / blocks;
        const Size take = std::min<Size>(end - begin, std::max<Size>(1, static_cast<Size>(std::llround((end - begin) * calibration_subset_ratio_))));
        std::partial_sort(tic_index.begin() + begin, tic_index.begin() + begin + take, tic_index.begin() + end,
                          [](const pair<double, Size>& a, const pair<double, Size>& c) { return a.first > c.first; });
        for (Size i = begin; i < begin + take; ++i) { subset.push_back(tic_index[i].second); }
      }
    }
    const Size subset_size = subset.size();
    // Scored in the order of the precursor mass, as the main search (scoringOrder_()): neighbouring iterations probe the
    // same precursor-window slices of the fragment index, which stay in the caches. The hits are collected in the order
    // of the subset (per tenth of the run, by TIC), as one thread scored them before: the sort by score below is not
    // stable, so this order decides which of two hits with the same score at the cut survives; it depends neither on
    // the threads nor on the chunks.
    std::vector<Size> order; // positions in subset, by precursor mass
    {
      std::vector<SignedSize> position_in_subset(spectra.size(), -1);
      for (Size k = 0; k < subset_size; ++k) position_in_subset[subset[k]] = static_cast<SignedSize>(k);
      order.reserve(subset_size);
      for (SignedSize i : scoringOrder_(spectra, query_spectra))
      {
        if (position_in_subset[i] >= 0) order.push_back(static_cast<Size>(position_in_subset[i]));
      }
    }

    OPENMS_LOG_INFO << "[ProSE] Calibration: scoring " << subset_size << " / " << spectra.size()
                    << " spectra (top TIC per tenth of the run)..." << std::endl;

    // Score subset and collect errors from the best hit per spectrum, with the ion series
    // the fragment index was built from (e.g. c/z+1 for ETD, where b/y would match nothing)
    const SpectrumGenerators_ generators = spectrumGenerators_();

    // Collect per-spectrum best hits with scores and errors
    // (rt, prec_mz [theoretical] and frag_mzs [theoretical, aligned with frag_errors] feed the models of calibration:recalibrate)
    struct CalHit { double score; double prec_error; std::vector<float> frag_errors; double rt; double prec_mz; std::vector<float> frag_mzs; };
    // At most one hit per spectrum of the subset, in its slot (written by one iteration only): collected in subset order
    // after the loop
    vector<CalHit> slots(subset_size);
    std::vector<unsigned char> filled(subset_size, 0);

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
    // Runs of neighbouring masses per core, at least 8 chunks per thread; one spectrum per chunk in peptide-mass mode
    // (as the main loop)
#ifdef _OPENMP
    const int cal_chunk = fragment_index.isPeptideMassMode() ? 1 : static_cast<int>(std::clamp<SignedSize>(
      static_cast<SignedSize>(subset_size) / (SignedSize(8) * std::max(1, omp_get_max_threads())), 1, 64));
#else
    [[maybe_unused]] const int cal_chunk = 1; // (only the OpenMP directive uses it)
#endif
#pragma omp parallel for schedule(dynamic, cal_chunk)
    for (SignedSize oi = 0; oi < (SignedSize)subset_size; ++oi)
    {
      if (calibration_failed.load(std::memory_order_relaxed)) continue;
      try
      {
      const Size si = order[oi];
      const Size scan_idx = subset[si];
      const MSSpectrum& spec = spectra[scan_idx];
      const TheoreticalSpectrumGenerator& tsg = generators.forSpectrum(spec);

      // The per-thread containers of the main loop (scoreSpectraAgainstIndex_()): the candidates and the theoretical
      // spectrum keep the capacity the thread's largest spectrum needed instead of taking fresh buffers per spectrum
      ScoringScratch_& scratch = scoringScratch_();
      FragmentIndex::SpectrumMatchesTopN& top_sms = scratch.top_sms;
      top_sms.clear();
      const MSSpectrum& query = query_spectra != nullptr ? (*query_spectra)[scan_idx] : spec;
      fragment_index.querySpectrum(query, db, top_sms, generators.electronIons(spec));

      // Find the best-scoring hit for this spectrum
      double best_score = 0;
      AASequence best_seq;
      int best_isotope_error = 0;
      uint16_t best_charge = 0;

      PeakSpectrum& theo = scratch.theo_spectrum;
      std::vector<uint8_t>& theo_ion_series = scratch.theo_ion_series;

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& peptide = fragment_index.getPeptides()[sms.peptide_idx_];
        AASequence seq = fragment_index.isPeptideMassMode()
          ? fragment_index.reconstructRealizedSubSequence(peptide, db, sms.realized_length_, sms.subset_bitmask_)
          : fragment_index.reconstructModifiedSequence(peptide, db);

        // The calibration PSMs are selected with the search's own score (upstream 825c33bb). The mass-accuracy score
        // favours PSMs with small fragment errors, so it narrows the fragment tolerance estimated from them below
        // (2-10% narrower than with HyperScore on 12 of 14 ppm runs); selecting by HyperScore instead did not
        // change the yield measurably.
        // The peaks of getSpectrum() with the ion series per peak instead of the ion names (no strings, no data
        // arrays), as the main loop: the same scores. The string path only when the generator cannot (losses, ...).
        HyperScore::PSMDetail detail;
        double score;
        if (tsg.getIonSeriesSpectrum(theo, theo_ion_series, seq, 1, scoringMaxCharge_(sms.precursor_charge_)))
        {
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, theo_ion_series, mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, theo_ion_series, detail);
        }
        else
        {
          theo.clear(false); // peaks and data arrays; getSpectrum appends, its output is sorted with add_metainfo=true
          theo.getPrecursors().clear(); // (clear(false) keeps the meta data; getSpectrum() appends a precursor per call)
          tsg.getSpectrum(theo, seq, 1, scoringMaxCharge_(sms.precursor_charge_));
          score = mass_accuracy_score_
            ? HyperScore::computeMassAccuracy(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, mass_error_sd_ppm_, detail)
            : HyperScore::computeWithDetail(pass_fragment, fragment_mass_tolerance_unit_ppm, spec, theo, detail);
        }

        if (score > best_score)
        {
          best_score = score;
          best_seq = std::move(seq);
          best_isotope_error = sms.isotope_error_;
          best_charge = sms.precursor_charge_;
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
      double corrected_exp_mz = monoisotopicPrecursorMZ(exp_mz, best_isotope_error, static_cast<double>(best_charge));
      double prec_err = (precursor_mass_tolerance_unit_ == "ppm")
                          ? Math::getPPM(corrected_exp_mz, theo_mz)
                          : (corrected_exp_mz - theo_mz);

      // The signed errors of the fragment peaks the best hit explains, within the pass window
      std::vector<float> frag_errors, frag_mzs;
      if (!tsg.getIonSeriesSpectrum(theo, theo_ion_series, best_seq, 1, scoringMaxCharge_(best_charge))) // (the same peaks)
      {
        theo.clear(false);
        theo.getPrecursors().clear();
        tsg.getSpectrum(theo, best_seq, 1, scoringMaxCharge_(best_charge));
      }
      matchedFragmentErrors(spec, theo, pass_fragment, fragment_mass_tolerance_unit_ppm, frag_errors, &frag_mzs);

      slots[si] = CalHit{best_score, prec_err, std::move(frag_errors), spec.getRT(), theo_mz, std::move(frag_mzs)};
      filled[si] = 1;
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
    vector<CalHit> cal_hits;
    cal_hits.reserve(static_cast<Size>(std::count(filled.begin(), filled.end(), 1)));
    for (Size k = 0; k < subset_size; ++k) { if (filled[k]) cal_hits.push_back(std::move(slots[k])); }
    vector<CalHit>().swap(slots);

    // The confident half: the top 50 % by score (random matches inflate the tails), if that leaves calibration_min_psms_
    Size cal_hits_total = cal_hits.size();
    std::sort(cal_hits.begin(), cal_hits.end(),
              [](const CalHit& a, const CalHit& b) { return a.score > b.score; });
    Size keep = cal_hits.size() / 2;
    if (keep >= calibration_min_psms_ && keep < cal_hits.size())
    {
      cal_hits.resize(keep);
    }

    // The hits the pass window delivers: the signed error e = observed - theoretical lies in [-upper, +lower] (the
    // mass-space window is [observed - lower, observed + upper]; the bounds flip for asymmetric windows)
    vector<double> precursor_errors, precursor_mzs, precursor_rts;
    vector<double> fragment_errors, fragment_mzs, fragment_rts;
    for (const auto& h : cal_hits)
    {
      if (h.prec_error < -pass_upper || h.prec_error > pass_lower) continue; // wrong match
      precursor_errors.push_back(h.prec_error);
      precursor_mzs.push_back(h.prec_mz);
      precursor_rts.push_back(h.rt);
      fragment_errors.insert(fragment_errors.end(), h.frag_errors.begin(), h.frag_errors.end());
      fragment_mzs.insert(fragment_mzs.end(), h.frag_mzs.begin(), h.frag_mzs.end());
      fragment_rts.insert(fragment_rts.end(), h.frag_errors.size(), h.rt);
    }

    if (precursor_errors.size() < calibration_min_psms_)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: insufficient PSMs (" << precursor_errors.size()
                      << " < " << calibration_min_psms_
                      << "), using configured tolerances." << std::endl;
      return result; // success=false
    }

    // The share of the errors as measured (before a recalibration) that the measured error recordRecalibration_()
    // reports covers: precursor offset +/- 4 sigma, fragment |offset| + 4 sigma. A measurement that misses more than the
    // windows below accept (two error modes, heavy tails) is recorded as not trusted.
    result.precursor_measured_coverage = measuredCoverage_(precursor_errors, false);
    if (fragment_errors.size() >= calibration_min_psms_) result.fragment_measured_coverage = measuredCoverage_(fragment_errors, true);

    // --- calibration:recalibrate: the error models, which applyRecalibration_() subtracts from the spectra. The
    //     windows below are then derived from the residual errors, the errors the main search sees. ---
    const bool recalibrate = calibration_recalibrate_ != "none";
    if (recalibrate)
    {
      const bool linear = calibration_recalibrate_ == "linear";
      const auto robust = [](std::vector<double> values)
      {
        const auto [med, mad] = medianAndMADBySelection(values);
        return std::make_pair(med, 1.4826 * mad);
      };
      std::tie(result.precursor_raw_shift, result.precursor_raw_sigma) = robust(precursor_errors);
      result.precursor_recalibration = MzRecalibration_::fit(precursor_errors, precursor_mzs, precursor_rts,
        calibration_recalibrate_rt_chunks_, calibration_min_psms_, linear, precursor_ppm);
      for (Size i = 0; i < precursor_errors.size(); ++i)
      {
        precursor_errors[i] -= result.precursor_recalibration.error(precursor_mzs[i], precursor_rts[i]);
      }
      if (fragment_errors.size() >= calibration_min_psms_)
      {
        std::tie(result.fragment_raw_shift, result.fragment_raw_sigma) = robust(fragment_errors);
        // chunks of as many PSMs as the precursor model's: the peaks of calibration:min_psms PSMs at least
        const Size min_peaks = calibration_min_psms_ * fragment_errors.size() / precursor_errors.size();
        result.fragment_recalibration = MzRecalibration_::fit(fragment_errors, fragment_mzs, fragment_rts,
          calibration_recalibrate_rt_chunks_, min_peaks, linear, fragment_mass_tolerance_unit_ppm);
        for (Size i = 0; i < fragment_errors.size(); ++i)
        {
          fragment_errors[i] -= result.fragment_recalibration.error(fragment_mzs[i], fragment_rts[i]);
        }
      }
    }

    // --- The calibrated windows: the median (offset) and 1.4826 * MAD (spread) of the signed errors are robust against
    //     the false matches among the PSMs; quantiles of the errors are not (the accepted false matches spread over the
    //     whole window and make up the tails). Precursor: [median - k sigma, median + k sigma], asymmetric about zero
    //     when the run has an offset; fragment: |median| + k sigma, symmetric (fragment masses are not shifted). Both only
    //     tighten: within [the floor, the configured tolerance]. ---
    const double k = calibration_sigma_factor_;
    const double kp = calibration_precursor_sigma_factor_; // 0: the precursor window is kept
    const double min_tolerance = 1e-6; // avoid non-positive tolerances
    // The share of the calibration PSMs (fragment peaks) a window must cover to be trusted. Precursor errors have heavy
    // tails (5 % of the PSMs of an Astral run beyond 8 ppm at sigma 0.7), and a window that cuts them costs more
    // identifications than it gains (the yield is flat from 5 to 50 ppm): 98 %. Fragment windows gain (10 ppm on
    // Astral: +4 %): 95 %.
    constexpr double min_coverage_precursor = 0.98, min_coverage_fragment = 0.95;
    {
      std::vector<double> errs = precursor_errors; // reordered by the selection
      const auto [med, mad] = medianAndMADBySelection(errs); // = Math::median() / Math::MAD()
      const double sigma = 1.4826 * mad;
      result.precursor_shift = med;
      result.precursor_sigma = sigma;
      // Window of the precursor errors e = observed - theoretical in [-cal_upper, +cal_lower]; with kp 0 the configured
      // window: the calibration then measures and reports the precursor errors and tightens the fragment tolerance only
      // (precursor yields are flat from 5 to 50 ppm; a tightened precursor window cost identifications on every Orbitrap
      // HCD set of the 20-run benchmark at equal entrapment FDP, HF-X -1.6 %, Lumos LFQ -2.8 %)
      const double lo = kp > 0 ? med - kp * sigma : -precursor_mass_tolerance_upper_;
      const double hi = kp > 0 ? med + kp * sigma : precursor_mass_tolerance_lower_;
      const double cal_lower_raw = std::max(min_tolerance, hi);
      const double cal_upper_raw = std::max(min_tolerance, -lo);
      // Zero spread (a synthetic fixture, a single peptide): nothing to calibrate the precursor window on, it is kept;
      // the fragment calibration still applies. Unit-aware: 1 ppm, or min_tolerance in Da.
      result.extreme_bias = kp > 0 && (hi - lo) < (precursor_ppm ? 1.0 : min_tolerance);
      result.precursor_spread = std::max(cal_lower_raw, cal_upper_raw);
      // The median and MAD describe one mode. A run whose offset changes (a lock mass that drops out: half of the
      // precursors near 0 ppm, half near +10 ppm) has errors the window does not cover; then the measurement is not
      // trusted and the configured window stays.
      result.precursor_coverage = std::count_if(precursor_errors.begin(), precursor_errors.end(),
                                                [&](double e) { return e >= lo && e <= hi; }) / static_cast<double>(precursor_errors.size());
      if (kp <= 0)
      {
        result.cal_lower = precursor_mass_tolerance_lower_;
        result.cal_upper = precursor_mass_tolerance_upper_;
      }
      else if (!result.extreme_bias && result.precursor_coverage < min_coverage_precursor)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the precursor window median +/- " << kp << " sigma covers only "
                        << static_cast<int>(100.0 * result.precursor_coverage) << " % of the calibration PSMs (heavy tails, or a mass "
                        << "offset that changes over the run): the configured precursor window is kept." << std::endl;
        result.cal_lower = precursor_mass_tolerance_lower_;
        result.cal_upper = precursor_mass_tolerance_upper_;
      }
      else if (!result.extreme_bias)
      {
        const double floor_side = calibration_min_precursor_tol_ > 0 ? calibration_min_precursor_tol_ : (precursor_ppm ? 2.5 : 0.005);
        result.cal_lower = std::clamp(cal_lower_raw, std::min(floor_side, precursor_mass_tolerance_lower_), precursor_mass_tolerance_lower_);
        result.cal_upper = std::clamp(cal_upper_raw, std::min(floor_side, precursor_mass_tolerance_upper_), precursor_mass_tolerance_upper_);
      }
      // else: cal_lower/cal_upper stay at 0; the writeback skips the precursor window.
    }
    result.psms = precursor_errors.size();
    result.fragment_peaks = fragment_errors.size();
    if (fragment_errors.size() >= calibration_min_psms_)
    {
      const auto [med, mad] = medianAndMADBySelection(fragment_errors);
      const double sigma = 1.4826 * mad;
      result.fragment_shift = med;
      result.fragment_sigma = sigma;
      const double measured = std::fabs(med) + k * sigma;
      // (recalibrated: the residual errors have no offset, so the floor is half)
      const double floor = calibration_min_fragment_tol_ > 0 ? calibration_min_fragment_tol_
                         : (fragment_mass_tolerance_unit_ppm ? (recalibrate ? 5.0 : 10.0) : (recalibrate ? 0.01 : 0.02));
      result.fragment_tolerance = std::clamp(measured, std::min(floor, fragment_mass_tolerance_), fragment_mass_tolerance_);
      result.fragment_coverage = std::count_if(fragment_errors.begin(), fragment_errors.end(),
                                               [&](double e) { return std::fabs(e) <= result.fragment_tolerance; }) / static_cast<double>(fragment_errors.size());
      if (result.fragment_coverage < min_coverage_fragment)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment tolerance |median| + " << k << " sigma covers only "
                        << static_cast<int>(100.0 * result.fragment_coverage) << " % of the matched fragment peaks of the calibration PSMs: "
                        << "the configured fragment tolerance is kept." << std::endl;
        result.fragment_tolerance = fragment_mass_tolerance_;
      }
      if (!fragment_mass_tolerance_unit_ppm && measured < 0.05)
      {
        OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment errors (|offset| + " << k << " sigma = " << std::fixed << std::setprecision(3)
                        << measured << " Da) are those of high-resolution MS2 spectra, but the fragment tolerance is given in Da ("
                        << fragment_mass_tolerance_ << "); a ppm tolerance would suit the data." << std::endl;
      }
    }
    else
    {
      result.fragment_tolerance = fragment_mass_tolerance_;
    }

    // Trusted (after the precursor and fragment fits above): the measured error covers the errors (one mode, no heavy
    // tails), and it lies inside the window the calibration pass searched with. A measured range that reaches that window
    // was censored (errors beyond it were not matched, so the median and MAD come from a truncated sample) and is a lower
    // bound, not a tolerance to adopt.
    {
      const double p_off = recalibrate ? result.precursor_raw_shift : result.precursor_shift;
      const double p_sigma = recalibrate ? result.precursor_raw_sigma : result.precursor_sigma;
      const bool p_censored = p_off + 4.0 * p_sigma >= pass_lower || p_off - 4.0 * p_sigma <= -pass_upper;
      result.precursor_measurement_trusted = !p_censored && result.precursor_measured_coverage >= min_coverage_precursor;
      const double f_off = recalibrate ? result.fragment_raw_shift : result.fragment_shift;
      const double f_sigma = recalibrate ? result.fragment_raw_sigma : result.fragment_sigma;
      const bool f_censored = measuredError_(f_off, f_sigma) >= pass_fragment;
      result.fragment_measurement_trusted = fragment_errors.size() >= calibration_min_psms_ && !f_censored
                                            && result.fragment_measured_coverage >= min_coverage_fragment;
    }

    // A configured window tighter than the error of the run cuts correct matches; the calibration only tightens, so it is
    // reported (a precursor window of 4.5-5 ppm from a methods text or an SDRF is usually a post-recalibration value)
    const double prec_hi = result.precursor_shift + 4.0 * result.precursor_sigma, prec_lo = result.precursor_shift - 4.0 * result.precursor_sigma;
    if (prec_hi > precursor_mass_tolerance_lower_ || -prec_lo > precursor_mass_tolerance_upper_)
    {
      // In the parameters' terms: an error e = observed - theoretical > 0 needs mass_tolerance_lower >= e (the mass-space
      // window is [observed - lower, observed + upper]), e < 0 needs mass_tolerance_upper >= -e
      std::ostringstream sides;
      sides << std::fixed << std::setprecision(3);
      if (prec_hi > precursor_mass_tolerance_lower_)
      {
        sides << " precursor:mass_tolerance_lower " << precursor_mass_tolerance_lower_ << " " << precursor_mass_tolerance_unit_
              << " is below the measured +" << prec_hi << ";";
      }
      if (-prec_lo > precursor_mass_tolerance_upper_)
      {
        sides << " precursor:mass_tolerance_upper " << precursor_mass_tolerance_upper_ << " " << precursor_mass_tolerance_unit_
              << " is below the measured " << -prec_lo << ";";
      }
      OPENMS_LOG_WARN << "[ProSE] Calibration: the precursor window is tighter than the errors the search sees (observed - "
                      << "theoretical from " << std::fixed << std::setprecision(3) << prec_lo << " to " << prec_hi << " "
                      << precursor_mass_tolerance_unit_ << ", offset +/- 4 sigma):" << sides.str() << " this loses correct matches; "
                      << "ProSE gains from wide precursor windows, e.g. 20 ppm" << (recalibrate ? "" : " (or calibration:recalibrate offset)")
                      << "." << std::endl;
    }
    if (fragment_errors.size() >= calibration_min_psms_ && measuredError_(result.fragment_shift, result.fragment_sigma) > fragment_mass_tolerance_)
    {
      OPENMS_LOG_WARN << "[ProSE] Calibration: the fragment tolerance " << fragment_mass_tolerance_ << " " << fragment_mass_tolerance_unit_
                      << " is tighter than the errors the search sees (|offset| + 4 sigma = "
                      << measuredError_(result.fragment_shift, result.fragment_sigma) << ") and loses matched fragments." << std::endl;
    }

    if (!recalibrate)
    {
      result.precursor_raw_shift = result.precursor_shift;
      result.precursor_raw_sigma = result.precursor_sigma;
      result.fragment_raw_shift = result.fragment_shift;
      result.fragment_raw_sigma = result.fragment_sigma;
    }
    result.success = true;

    if (result.recalibrated())
    {
      OPENMS_LOG_INFO << "[ProSE]   Recalibration (" << calibration_recalibrate_ << "): precursor offset " << std::fixed
                      << std::setprecision(2) << result.precursor_raw_shift << ", sigma " << result.precursor_raw_sigma
                      << " -> residual offset " << result.precursor_shift << ", sigma " << result.precursor_sigma << " "
                      << precursor_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.precursor_raw_shift, result.precursor_raw_sigma)
                      << " -> " << measuredError_(result.precursor_shift, result.precursor_sigma) << "); model: " << result.precursor_recalibration.toString() << std::endl;
      if (result.fragment_recalibration.active())
      {
        OPENMS_LOG_INFO << "[ProSE]   Recalibration (" << calibration_recalibrate_ << "): fragment offset " << std::fixed
                        << std::setprecision(2) << result.fragment_raw_shift << ", sigma " << result.fragment_raw_sigma
                        << " -> residual offset " << result.fragment_shift << ", sigma " << result.fragment_sigma << " "
                        << fragment_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.fragment_raw_shift, result.fragment_raw_sigma)
                        << " -> " << measuredError_(result.fragment_shift, result.fragment_sigma) << "); model: " << result.fragment_recalibration.toString() << std::endl;
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Calibration: " << precursor_errors.size() << " PSMs used (top "
                    << static_cast<int>(100.0 * precursor_errors.size() / std::max<Size>(1, cal_hits_total)) << "% by score); pass windows [-"
                    << pass_lower << ", +" << pass_upper << "] " << precursor_mass_tolerance_unit_ << ", fragment "
                    << pass_fragment << " " << fragment_mass_tolerance_unit_ << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Precursor: offset " << std::fixed << std::setprecision(2) << result.precursor_shift
                    << ", sigma " << result.precursor_sigma << " " << precursor_mass_tolerance_unit_
                    << " (|offset| + 4 sigma " << measuredError_(result.precursor_shift, result.precursor_sigma) << ") -> window [-" << result.cal_lower << ", +" << result.cal_upper << "]"
                    << " (covers " << static_cast<int>(100.0 * result.precursor_coverage) << " % of the PSMs)"
                    << (kp <= 0 ? " (the configured window, calibration:precursor_sigma_factor 0)" : "")
                    << (result.extreme_bias ? " (zero spread: the configured window is kept)" : "") << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Fragment:  offset " << result.fragment_shift << ", sigma " << result.fragment_sigma
                    << " " << fragment_mass_tolerance_unit_ << " (|offset| + 4 sigma " << measuredError_(result.fragment_shift, result.fragment_sigma)
                    << "; " << result.fragment_peaks << " peaks) -> tolerance "
                    << fragment_mass_tolerance_ << " -> " << result.fragment_tolerance << " " << fragment_mass_tolerance_unit_
                    << " (covers " << static_cast<int>(100.0 * result.fragment_coverage) << " % of the peaks)" << std::endl;

    return result;
  }

  // =====================================================================
  // calibration:recalibrate: the m/z error models and their application to the spectra
  // =====================================================================
  ProSEAlgorithm::MzRecalibration_::AtRT ProSEAlgorithm::MzRecalibration_::at(double rt) const
  {
    AtRT model;
    model.ppm = ppm;
    if (chunks.empty()) return model;
    if (chunks.size() == 1 || rt <= chunks.front().rt) { model.lo = &chunks.front(); return model; }
    if (rt >= chunks.back().rt) { model.lo = &chunks.back(); return model; }
    // the chunks around rt: lo->rt <= rt < hi->rt
    const auto hi = std::upper_bound(chunks.begin(), chunks.end(), rt, [](double r, const Chunk& c) { return r < c.rt; });
    const auto lo = hi - 1;
    model.lo = &*lo;
    model.hi = &*hi;
    model.w = (rt - lo->rt) / (hi->rt - lo->rt);
    return model;
  }

  double ProSEAlgorithm::MzRecalibration_::error(double mz, double rt) const
  {
    return at(rt).error(mz);
  }

  double ProSEAlgorithm::MzRecalibration_::correct(double mz, double rt) const
  {
    return at(rt).correct(mz);
  }

  double ProSEAlgorithm::MzRecalibration_::uncorrect(double corrected, double rt) const
  {
    if (chunks.empty()) return corrected;
    // correct() maps mz to mz (1 - e(mz) 1e-6), or mz - e(mz) in Da: iterate mz = corrected / (1 - e(mz) 1e-6). Exact
    // after one step for 'offset' (its error does not depend on m/z); for a line in m/z each step shrinks the
    // remaining difference by a factor of about slope * mz * 1e-6, so three steps suffice.
    double mz = corrected;
    for (int step = 0; step < 3; ++step)
    {
      const double e = error(mz, rt);
      const double next = ppm ? corrected / (1.0 - e * 1e-6) : corrected + e;
      if (next == mz) break; // converged
      mz = next;
    }
    return mz;
  }

  std::string ProSEAlgorithm::MzRecalibration_::toString() const
  {
    if (chunks.empty()) return "";
    std::ostringstream s;
    s.imbue(std::locale::classic());
    s << (ppm ? "ppm" : "Da");
    for (const Chunk& c : chunks)
    {
      s << "; rt " << std::fixed << std::setprecision(1) << c.rt << ": " << std::setprecision(4) << c.intercept;
      if (c.slope != 0.0) { s << (c.slope < 0 ? " - " : " + ") << std::scientific << std::setprecision(3) << std::fabs(c.slope) << " mz"; }
      s << " (" << c.points << ")";
    }
    return s.str();
  }

  ProSEAlgorithm::MzRecalibration_ ProSEAlgorithm::MzRecalibration_::fit(const std::vector<double>& errors, const std::vector<double>& mzs,
    const std::vector<double>& rts, Size n_chunks, Size min_points, bool linear, bool ppm)
  {
    MzRecalibration_ model;
    model.ppm = ppm;
    const Size n = errors.size();
    if (n == 0 || mzs.size() != n || rts.size() != n) return model;
    std::vector<Size> order(n);
    std::iota(order.begin(), order.end(), Size(0));
    std::stable_sort(order.begin(), order.end(), [&rts](Size a, Size b) { return rts[a] < rts[b]; });
    n_chunks = std::max<Size>(1, std::min(n_chunks, n / std::max<Size>(1, min_points)));
    for (Size c = 0; c < n_chunks; ++c)
    {
      const Size begin = n * c / n_chunks, end = n * (c + 1) / n_chunks;
      std::vector<double> e, x;
      e.reserve(end - begin);
      x.reserve(end - begin);
      for (Size i = begin; i < end; ++i) { e.push_back(errors[order[i]]); x.push_back(mzs[order[i]]); }
      std::vector<double> sorted = e;
      const auto [med, mad] = medianAndMADBySelection(sorted);
      Chunk chunk{rts[order[(begin + end) / 2]], med, 0.0, e.size()}; // the centre: the median RT of the chunk
      if (linear)
      {
        // the least-squares line through the points within 3 sigma of the median (the false matches among the PSMs)
        const double cut = 3.0 * 1.4826 * mad;
        std::vector<double> te, tx;
        for (Size i = 0; i < e.size(); ++i) { if (std::fabs(e[i] - med) <= cut) { te.push_back(e[i]); tx.push_back(x[i]); } }
        MZTrafoModel line(ppm);
        double intercept = 0.0, slope = 0.0, power = 0.0;
        if (te.size() >= 2 && line.train(te, tx, std::vector<double>(te.size(), 1.0), MZTrafoModel::MODELTYPE::LINEAR, false))
        {
          line.getCoefficients(intercept, slope, power);
          if (std::isfinite(intercept) && std::isfinite(slope)) { chunk.intercept = intercept; chunk.slope = slope; chunk.points = te.size(); }
        }
      }
      model.chunks.push_back(chunk);
    }
    return model;
  }

  void ProSEAlgorithm::applyRecalibration_(const CalibrationResult_& cal, PeakMap& spectra, PeakMap* query_spectra,
                                           PeakMap* evidence_spectra) const
  {
    if (!cal.success || !cal.recalibrated()) return;
    const MzRecalibration_& precursor_model = cal.precursor_recalibration;
    const MzRecalibration_& fragment_model = cal.fragment_recalibration;
    MetaInfoRegistry& registry = MetaInfoInterface::metaRegistry();
    const UInt observed_key = registry.registerName(OBSERVED_PRECURSOR_MZ);
    // the co-isolated precursors of a chimeric search (CoIsolationDetector::annotate(); UInt(-1): none annotated in this
    // process), from the same survey scans as the selected precursor: corrected by its model
    const UInt coisolated_key = registry.getIndex(CoIsolationDetector::META_MZ);
    const UInt observed_coisolated_key = coisolated_key != UInt(-1) ? registry.registerName(OBSERVED_COISOLATED_MZ) : UInt(-1);
    // The model at the spectrum's RT once, per peak only its arithmetic (the chunk lookup per peak was a third of the
    // time of this function in stack samples of a 27-run search)
    const auto correct_peaks = [](MSSpectrum& s, const MzRecalibration_::AtRT& model)
    {
      for (Peak1D& p : s) { p.setMZ(model.correct(p.getMZ())); } // (monotone for any sane model: the order stays)
    };
    Size precursors = 0, coisolated = 0, peaks = 0;
#pragma omp parallel for schedule(static) reduction(+: precursors, coisolated, peaks)
    for (SignedSize i = 0; i < static_cast<SignedSize>(spectra.size()); ++i)
    {
      MSSpectrum& s = spectra[i];
      if (s.getMSLevel() != 2) continue;
      const double rt = s.getRT();
      if (precursor_model.active() && !s.getPrecursors().empty())
      {
        Precursor& p = s.getPrecursors()[0];
        if (!p.metaValueExists(observed_key)) { p.setMetaValue(observed_key, p.getMZ()); } // (a second correction keeps the first observation)
        const double corrected = precursor_model.correct(p.getMZ(), rt);
        p.setMZ(corrected);
        if (query_spectra != nullptr && !(*query_spectra)[i].getPrecursors().empty()) { (*query_spectra)[i].getPrecursors()[0].setMZ(corrected); }
        ++precursors;
      }
      if (precursor_model.active() && coisolated_key != UInt(-1) && s.metaValueExists(coisolated_key))
      {
        DoubleList mzs = s.getMetaValue(coisolated_key);
        if (!s.metaValueExists(observed_coisolated_key)) { s.setMetaValue(observed_coisolated_key, mzs); } // (as above)
        for (double& mz : mzs) { mz = precursor_model.correct(mz, rt); }
        coisolated += mzs.size();
        s.setMetaValue(coisolated_key, std::move(mzs));
      }
      if (fragment_model.active())
      {
        const MzRecalibration_::AtRT model = fragment_model.at(rt);
        correct_peaks(s, model);
        if (query_spectra != nullptr) { correct_peaks((*query_spectra)[i], model); }
        if (evidence_spectra != nullptr) { correct_peaks((*evidence_spectra)[i], model); }
        peaks += s.size();
      }
    }
    OPENMS_LOG_INFO << "[ProSE] Recalibration (" << calibration_recalibrate_ << "): corrected the precursor m/z of " << precursors
                    << " spectra" << (coisolated > 0 ? " and of " + std::to_string(coisolated) + " co-isolated precursors" : std::string())
                    << " and " << peaks << " fragment peaks in place (the PSMs report the observed m/z)" << std::endl;
  }

  double ProSEAlgorithm::measuredCoverage_(const std::vector<double>& errors, bool fragment)
  {
    if (errors.empty()) return 0.0;
    std::vector<double> values = errors; // reordered by the selection
    const std::pair<double, double> med_mad = medianAndMADBySelection(values);
    const double med = med_mad.first; // (no structured binding: clang cannot capture one in a lambda with OpenMP)
    const double sigma = 1.4826 * med_mad.second;
    const double measured = measuredError_(med, sigma);
    const auto covered = [&](double e) { return fragment ? std::fabs(e) <= measured : (e >= med - 4.0 * sigma && e <= med + 4.0 * sigma); };
    return std::count_if(errors.begin(), errors.end(), covered) / static_cast<double>(errors.size());
  }

  void ProSEAlgorithm::recordRecalibration_(const CalibrationResult_& cal, std::vector<ProteinIdentification>& protein_ids) const
  {
    if (protein_ids.empty() || !cal.success) return;
    ProteinIdentification::SearchParameters sp = protein_ids[0].getSearchParameters();
    // The mass errors the calibration pass measured: offset (median), sigma (1.4826 * MAD) and the measured error
    // |offset| + 4 sigma, the range that holds practically all correct matches (in the tolerance units)
    sp.setMetaValue("calibration:precursor_error_offset", cal.precursor_raw_shift);
    sp.setMetaValue("calibration:precursor_error_sigma", cal.precursor_raw_sigma);
    sp.setMetaValue("calibration:precursor_measured_error", measuredError_(cal.precursor_raw_shift, cal.precursor_raw_sigma));
    // (no fragment records without a fragment measurement, fewer than calibration:min_psms matched peaks: their zeros
    // would read as a perfect fragment calibration)
    const bool fragments_measured = cal.fragment_peaks >= calibration_min_psms_;
    if (fragments_measured)
    {
      sp.setMetaValue("calibration:fragment_error_offset", cal.fragment_raw_shift);
      sp.setMetaValue("calibration:fragment_error_sigma", cal.fragment_raw_sigma);
      sp.setMetaValue("calibration:fragment_measured_error", measuredError_(cal.fragment_raw_shift, cal.fragment_raw_sigma));
    }
    // Whether a measurement describes the run: it covers at least 98 % (precursor) or 95 % (fragment peaks) of the
    // calibration PSMs, lies inside the calibration pass's window (not censored), and the fragments come from
    // calibration:min_psms peaks or more. A pipeline that takes the measured error as a tolerance should not take an
    // untrusted one.
    sp.setMetaValue("calibration:psms", static_cast<Int>(cal.psms));
    sp.setMetaValue("calibration:precursor_measurement_trusted", cal.precursor_measurement_trusted ? "true" : "false");
    sp.setMetaValue("calibration:fragment_measurement_trusted", cal.fragment_measurement_trusted ? "true" : "false");
    if (cal.recalibrated())
    {
      sp.setMetaValue("calibration:recalibrate", calibration_recalibrate_);
      sp.setMetaValue("calibration:recalibrate_rt_chunks", static_cast<Int>(calibration_recalibrate_rt_chunks_));
      sp.setMetaValue("calibration:precursor_model", cal.precursor_recalibration.toString());
      sp.setMetaValue("calibration:precursor_residual_offset", cal.precursor_shift);
      sp.setMetaValue("calibration:precursor_residual_sigma", cal.precursor_sigma);
      sp.setMetaValue("calibration:precursor_residual_measured_error", measuredError_(cal.precursor_shift, cal.precursor_sigma));
      if (fragments_measured)
      {
        sp.setMetaValue("calibration:fragment_model", cal.fragment_recalibration.toString());
        sp.setMetaValue("calibration:fragment_residual_offset", cal.fragment_shift);
        sp.setMetaValue("calibration:fragment_residual_sigma", cal.fragment_sigma);
        sp.setMetaValue("calibration:fragment_residual_measured_error", measuredError_(cal.fragment_shift, cal.fragment_sigma));
      }
    }
    protein_ids[0].setSearchParameters(std::move(sp));
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
    y << "  peptide_mass_mode: " << bol(sh.peptide_mass_mode) << "\n";
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
    y << "  seconds_output: " << num(sh.seconds_output) << "\n";
    y << "  seconds_qpx_extended: " << num(sh.seconds_qpx_extended) << "\n";

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
        y << "    profile_spectra: " << st.profile_spectra << "\n";
        y << "    charge_out_of_range_spectra: " << st.charge_out_of_range_spectra << "\n";
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
        if (!st.recalibration.empty()) // calibration: the measured errors; with recalibrate the mode, models and residuals
        {
          y << "    calibration:\n";
          for (const auto& [name, value] : st.recalibration)
          {
            y << "      " << name << ": ";
            if (value.valueType() == DataValue::DOUBLE_VALUE) { y << num(static_cast<double>(value)); }
            else if (value.valueType() == DataValue::INT_VALUE) { y << static_cast<long>(value); }
            else { y << q(value.toString()); }
            y << "\n";
          }
        }
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
  std::string ProSEAlgorithm::renderTotalTime(const SharedSearchStats& sh)
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << "total " << sh.seconds_total << "s | search of the file(s) "
       << std::max(0.0, sh.seconds_total - sh.seconds_output) << "s | outputs " << sh.seconds_output << "s";
    if (sh.seconds_qpx_extended > 0.0) { os << " (of which -qpx_extended " << sh.seconds_qpx_extended << "s)"; }
    return os.str();
  }

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
    if (s.ms2_spectra == 0)
    {
      os << "[ProSE] Warning      : no MS2 spectra in the input (e.g. an MS1-only file): nothing was searched\n";
    }
    if (s.profile_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.profile_spectra << " MS2 spectra are profile spectra; centroid them first\n";
    }
    if (s.charge_out_of_range_spectra > 0)
    {
      os << "[ProSE] Warning      : " << s.charge_out_of_range_spectra << " MS2 spectra not searched: precursor charge outside of "
         << sh.min_charge << "-" << sh.max_charge << "\n";
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
       << (sh.recalibrate.empty() || sh.recalibrate == "none" ? std::string() : " (recalibrate " + sh.recalibrate + ")")
       << " | mode: " << (sh.open_search ? "open" : "closed")
       << (sh.peptide_mass_mode ? " | low_memory" : "")
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
    if (s.seconds_search > 0.0)      { os << "scoring " << s.seconds_search << "s | "; }
    if (s.seconds_fdr > 0.0)         { os << "fdr " << s.seconds_fdr << "s | "; }
    os << "index(shared) " << sh.seconds_index_build << "s\n";
    // the whole run (set by the caller): the search of the file (reading it, rescoring, FDR, ...) and the outputs
    if (sh.seconds_total > 0.0) { os << "[ProSE] " << renderTotalTime(sh) << "\n"; }

    // -- Modification discovery (open search only) --
    if (is_open_search)
    {
      renderModificationSummary(mod_analysis, os);
    }
  }
} // namespace OpenMS
