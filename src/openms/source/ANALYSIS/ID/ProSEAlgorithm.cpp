// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#endif
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
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
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <locale>
#include <ostream>
#include <map>
#include <set>
#include <sstream>
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

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0, 1};
    defaults_.setValue("precursor:isotopes", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak)");

    defaults_.setValue("fragment:mass_tolerance", 20.0, "Fragment mass tolerance");

    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    defaults_.setValue("fragment:deisotope", "auto", "MS2 deisotoping (single-charge deconvolution) before searching. 'auto' deisotopes only when the fragment tolerance is within the high-resolution deisotoper range (<= 0.1 Da / <= 100 ppm) and skips it for low-resolution (e.g. ion-trap CID) data; 'true' always deisotopes (requires a high-resolution fragment tolerance); 'false' never deisotopes.");
    defaults_.setValidStrings("fragment:deisotope", {"auto", "true", "false"});


    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("fragment:min_ion_index", 2, "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.");
    defaults_.setMinInt("fragment:min_ion_index", 0);

    defaults_.setValue("peaks:keep_n", 0, "Maximum number of MS2 peaks kept per spectrum (NLargest) before scoring. 0 = auto: a resolution-aware cap (~400 for high-res, ~80 at 0.5 Da, floor 60) that avoids the spurious low-resolution fragment matching which otherwise inflates HyperScore for targets and decoys alike. Set an explicit value to override (e.g. 400 for the legacy behavior).", {"advanced"});
    defaults_.setMinInt("peaks:keep_n", 0);
    defaults_.setValue("peaks:window_top", 20, "Maximum number of MS2 peaks kept per 100 Da window (WindowMower) before scoring.", {"advanced"});
    defaults_.setMinInt("peaks:window_top", 1);

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<std::string> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
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

    defaults_.setSectionDescription("annotate", "Annotation Options");

    defaults_.setValue("peptdeep:enable", "false", "Add PeptDeep prediction-based rescoring features (ms2_cosine, ms2_spectral_angle, ms2_pearson, ms2_frac_pred_found, rt_abs_error) to every PSM. Uses the models shipped in OpenMS' 'share/OpenMS/models' unless 'peptdeep:ms2_model'/'peptdeep:rt_model' name others. Requires an OpenMS built with ONNX support.");
    defaults_.setValidStrings("peptdeep:enable", {"true", "false"});
    defaults_.setValue("peptdeep:ms2_model", "models/peptdeep_ms2_dynamic.onnx", "PeptDeep MS2 fragment-intensity ONNX model, used when 'peptdeep:enable' is true. A relative name is resolved against OpenMS' shared-data directory; an absolute path is used as given.");
    defaults_.setValue("peptdeep:rt_model", "models/peptdeep_rt_dynamic.onnx", "PeptDeep retention-time ONNX model. See 'peptdeep:ms2_model'.");
    defaults_.setValue("peptdeep:instrument", "QE", "Instrument class passed to the PeptDeep MS2 model.");
    defaults_.setValidStrings("peptdeep:instrument", {"Lumos", "QE", "timsTOF", "SciexTOF"});
    defaults_.setValue("peptdeep:nce", -1.0, "Normalised collision energy for the PeptDeep MS2 model. Negative selects it automatically from the collision energy recorded in the spectra, refined by scoring a small grid on confident PSMs.");
    defaults_.setValue("peptdeep:rt_model_type", "b_spline", "Model mapping predicted onto observed retention time.", {"advanced"});
    defaults_.setValidStrings("peptdeep:rt_model_type", {"b_spline", "lowess", "linear"});
    defaults_.setValue("peptdeep:batch_size", 500, "Peptides per ONNX inference call.", {"advanced"});
    defaults_.setMinInt("peptdeep:batch_size", 1);
    defaults_.setSectionDescription("peptdeep", "PeptDeep (ONNX) prediction-based rescoring features");

    defaults_.setValue("peptide:min_size", 7, "Minimum size a peptide must have after digestion to be considered in the search.");
    defaults_.setValue("peptide:max_size", 40, "Maximum size a peptide must have after digestion to be considered in the search (0 = disabled).");
    defaults_.setValue("peptide:missed_cleavages", 1, "Number of missed cleavages.");
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
    defaults_.setSectionDescription("report", "Reporting Options");

    defaults_.setValue("FDR:PSM", 0.0, "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, set to 0 to disable filtering and report all PSMs with q-values). Target and decoy PSMs are filtered alike by the q-value threshold; decoys that pass are kept (no decoy-specific stripping here — decoys are removed only at protein-FDR finalization). Requires '-decoys' to be set.");
    defaults_.setMinFloat("FDR:PSM", 0.0);
    defaults_.setMaxFloat("FDR:PSM", 1.0);
    defaults_.setValue("FDR:protein", 0.0, "Filter proteins based on picked-protein FDR q-value (e.g., 0.01 = 1% protein FDR, set to 0 to disable). Applied after PSM-level FDR on a complete protein set (single file, or the -out_merged aggregate). Setting this > 0 finalizes the result: identified decoys are removed. With 0, decoys are retained for downstream/merged FDR. Uses the picked-protein approach (Savitski et al. 2015) which pairs target and decoy proteins by accession. Requires '-decoys' to be set.");
    defaults_.setMinFloat("FDR:protein", 0.0);
    defaults_.setMaxFloat("FDR:protein", 1.0);
    defaults_.setSectionDescription("FDR", "False Discovery Rate control (requires decoys)");

    // Add parameters which are only used by FragmentIndex
    defaults_.setValue("peptide:min_mass", 100, "Minimal peptide mass for database");
    defaults_.setValue("peptide:max_mass", 9000, "Maximal peptide mass for database");

    // Fragment-level filtering
    defaults_.setValue("fragment:min_matched_ions", 5, "Minimal number of matched ions to report a PSM");

    // Precursor isotope error handling
    defaults_.setValue("precursor:isotope_error_min", -1, "Minimum allowed precursor isotope error");
    defaults_.setValue("precursor:isotope_error_max", 1, "Maximum allowed precursor isotope error");

    // Fragment and scoring limits
    defaults_.setValue("fragment:max_charge", 2, "max fragment charge");
    defaults_.setValue("scoring:max_candidates_per_spectrum", 50, "The number of initial hits for which we calculate a score");
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
    defaults_.setValue("ions:add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("ions:add_z_ions", {"true","false"});
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

    peptdeep_enable_ = param_.getValue("peptdeep:enable").toString() == "true";
    peptdeep_ms2_model_ = param_.getValue("peptdeep:ms2_model").toString();
    peptdeep_rt_model_ = param_.getValue("peptdeep:rt_model").toString();
    peptdeep_instrument_ = param_.getValue("peptdeep:instrument").toString();
    peptdeep_nce_ = param_.getValue("peptdeep:nce");
    peptdeep_rt_model_type_ = param_.getValue("peptdeep:rt_model_type").toString();
    peptdeep_batch_size_ = param_.getValue("peptdeep:batch_size");

    precursor_isotopes_ = param_.getValue("precursor:isotopes");
    peaks_keep_n_ = (Size)(int)param_.getValue("peaks:keep_n");
    peaks_window_top_ = (Int)param_.getValue("peaks:window_top");

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");

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

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    set<std::string> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());
    if (fixed_unique.size() != modifications_fixed_.size())
    {
      OPENMS_LOG_WARN << "Duplicate fixed modification provided. Making them unique." << endl;
      modifications_fixed_.assign(fixed_unique.begin(), fixed_unique.end());
    }    

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<std::string> var_unique(modifications_variable_.begin(), modifications_variable_.end());
    if (var_unique.size() != modifications_variable_.size())
    {
      OPENMS_LOG_WARN << "Duplicate variable modification provided. Making them unique." << endl;
      modifications_variable_.assign(var_unique.begin(), var_unique.end());
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

    const std::string decoy_mode_str = param_.getValue("decoys").toString();
    if (decoy_mode_str == "generate")   { decoy_mode_ = DecoyMode_::GENERATE; }
    else if (decoy_mode_str == "ignore") { decoy_mode_ = DecoyMode_::IGNORE; }
    else                                 { decoy_mode_ = DecoyMode_::AUTO; }
    decoy_prefix_ = param_.getValue("decoy_prefix").toString();
    annotate_psm_ = ListUtils::toStringList<std::string>(param_.getValue("annotate:PSM"));
    fdr_psm_ = param_.getValue("FDR:PSM");
    fdr_protein_ = param_.getValue("FDR:protein");

    // Open search mode is automatically determined based on precursor tolerance in isOpenSearchMode_()

    add_a_ions_ = param_.getValue("ions:add_a_ions").toBool();
    add_b_ions_ = param_.getValue("ions:add_b_ions").toBool();
    add_c_ions_ = param_.getValue("ions:add_c_ions").toBool();
    add_x_ions_ = param_.getValue("ions:add_x_ions").toBool();
    add_y_ions_ = param_.getValue("ions:add_y_ions").toBool();
    add_z_ions_ = param_.getValue("ions:add_z_ions").toBool();

    database_chunk_size_ = param_.getValue("database:chunk_size");

    calibration_enabled_ = param_.getValue("calibration:enabled") == "true";
    calibration_subset_ratio_ = param_.getValue("calibration:subset_ratio");
    calibration_min_psms_ = param_.getValue("calibration:min_psms");
  }

  // static
  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool deisotope_requested, Size peaks_keep_n, Int peaks_window_top)
  {
    // Intensity threshold + normalization used to run here as two extra SERIAL full-map
    // passes. Both are strictly per-spectrum: ThresholdMower::filterPeakMap and
    // Normalizer::filterPeakMap are literally "for (auto& s : exp) filterSpectrum(s);" and
    // neither iterates chromatograms. They are therefore applied at the top of the parallel
    // loop below instead, which is per-spectrum equivalent and removes two full sweeps over
    // the peak data. Both objects are configured once here and shared across the OpenMP
    // threads, exactly like window_mower_filter / nlargest_filter below:
    // Normalizer::filterPeakSpectrum is const (reads the resolved 'method_' only), and
    // ThresholdMower only re-reads its 'threshold' Param into a member on every call -- the
    // same idempotent same-value write the already-shared WindowMower performs.
    ThresholdMower threshold_mower_filter;
    Normalizer normalizer;

    // sort by rt; done before the loop because sortSpectra(false) only permutes whole
    // spectra (it does not touch peak order, see MSExperiment::sortSpectra) and therefore
    // commutes with the per-spectrum filters that now run inside the loop.
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", peaks_window_top, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

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

#pragma omp parallel for default(none) shared(exp, do_deisotope, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, threshold_mower_filter, normalizer, window_mower_filter, nlargest_filter)
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // remove 0 intensities, then normalize (formerly two serial full-map passes)
      threshold_mower_filter.filterPeakSpectrum(exp[exp_index]);
      normalizer.filterPeakSpectrum(exp[exp_index]);

      // sort by mz
      exp[exp_index].sortByPosition();

      // deisotope (skipped for low-resolution data; see do_deisotope above)
      if (do_deisotope)
      {
        Deisotoper::deisotopeAndSingleCharge(exp[exp_index],
          fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,
          1, 3,   // min / max charge
          false,  // keep only deisotoped
          3, 10,  // min / max isopeaks
          true);  // convert fragment m/z to mono-charge
      }

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
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

  void ProSEAlgorithm::annotatePeptDeepFeatures_(const PeakMap& spectra,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    if (!peptdeep_enable_) { return; }
#ifdef WITH_ONNX
    startProgress(0, 1, "Adding PeptDeep rescoring features...");
    PeptDeepRescoring rescoring;
    Param p = rescoring.getParameters();
    // A bare name resolves against share/OpenMS, an absolute path is returned unchanged.
    // The models are downloaded to the build tree's share/OpenMS/models, which is not the
    // compiled-in data path, so search relative to the executable as well: '../share/OpenMS'
    // holds them both in a build tree (bin/../share) and in an install tree.
    const StringList model_dirs = {File::getExecutablePath() + "../share/OpenMS"};
    p.setValue("ms2_model", File::find(peptdeep_ms2_model_, model_dirs));
    p.setValue("rt_model", File::find(peptdeep_rt_model_, model_dirs));
    p.setValue("instrument", peptdeep_instrument_);
    p.setValue("nce", peptdeep_nce_);
    p.setValue("rt_model_type", peptdeep_rt_model_type_);
    p.setValue("batch_size", peptdeep_batch_size_);
    // Inference otherwise runs at its own default while the search around it scales with
    // the thread count the user actually asked for.
#ifdef _OPENMP
    p.setValue("threads", std::max(1, omp_get_max_threads()));
#else
    p.setValue("threads", 1);
#endif
    rescoring.setParameters(p);
    rescoring.setLogType(getLogType());
    rescoring.annotate(spectra, protein_ids, peptide_ids);
    endProgress();
#else
    (void)spectra; (void)protein_ids; (void)peptide_ids;
    OPENMS_LOG_WARN << "[ProSE] 'peptdeep:enable' is set, but this OpenMS was built without "
                       "ONNX support (WITH_ONNX=OFF); the prediction-based rescoring features "
                       "are not added." << '\n';
#endif
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
        const std::string& database_name) const
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
    // loop in scoreSpectraAgainstIndex_ shares one generator the same way).
    // The alignment tolerance mirrors the search's fragment tolerance so the
    // reported FRAGMENT_ERROR_MEDIAN_PPM is not polluted by far-off spurious
    // matches — SpectrumAlignment's default is 0.3 Da absolute, which is ~30×
    // looser than a typical 20 ppm search.
    TheoreticalSpectrumGenerator tsg;
    {
      Param tsg_param(tsg.getParameters());
      tsg_param.setValue("add_metainfo", "true");
      tsg_param.setValue("add_first_prefix_ion", "true");
      // Mirror the configured ion series so annotation features (fragment error,
      // peak annotations, longest ion run, MIC) are computed on the same
      // theoretical spectrum the candidate was scored against.
      tsg_param.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
      tsg_param.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
      tsg_param.setValue("add_c_ions", add_c_ions_ ? "true" : "false");
      tsg_param.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
      tsg_param.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
      tsg_param.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
      tsg.setParameters(tsg_param);
    }
    SpectrumAlignment sa;
    {
      Param sa_param(sa.getParameters());
      sa_param.setValue("tolerance", fragment_mass_tolerance);
      sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
      sa.setParameters(sa_param);
    }

#pragma omp parallel for
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (!annotated_hits[scan_index].empty())
      {
        const MSSpectrum& spec = exp[scan_index];
        // create empty PeptideIdentification object and fill meta data
        PeptideIdentification pi{};
        pi.setSpectrumReference( spec.getNativeID());
        pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
        pi.setScoreType("ln(hyperscore)");
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
          if (need_alignment)
          {
            const int max_frag_z = (charge >= 2) ? std::min<int>(charge - 1, 2) : 1;
            tsg.getSpectrum(theoretical_spec, ah.sequence, 1, max_frag_z);
            sa.getSpectrumAlignment(alignment, theoretical_spec, spec);
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
            ph.setMetaValue(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM, median_ppm_error);
          }

          if (annotation_precursor_error_ppm)
          {
            // Subtract out the isotope offset FI matched at — FragmentIndex searches
            // shifted_mass = precursor_mass + isotope_error * C13C12, so M_theo ≈ N_obs
            // + isotope_error * C13C12, and the observed-to-monoiso correction in m/z is
            //   corrected_mz = observed_mz + isotope_error * C13C12 / charge
            // Without this, a ±1 Da FI match reports ~1000 ppm / charge for the Percolator
            // feature, corrupting target/decoy discrimination.
            const double corrected_mz = mz
              + static_cast<double>(ah.isotope_error) * Constants::C13C12_MASSDIFF_U / used_charge;
            double theo_mz = ah.sequence.getMZ(used_charge);
            double ppm_difference = Math::getPPM(corrected_mz, theo_mz);
            ph.setMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM, ppm_difference);
          }

          if (annotation_prefix_fraction)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION, ah.prefix_fraction);
          }

          if (annotation_suffix_fraction)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION, ah.suffix_fraction);
          }

          // Matched ion counts (from scoring, no alignment needed)
          if (annotation_num_matched_peaks)
          {
            ph.setMetaValue(Constants::UserParam::NUM_MATCHED_PEAKS, static_cast<int>(ah.matched_prefix_ions + ah.matched_suffix_ions));
          }
          if (annotation_matched_prefix_ions)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS, static_cast<int>(ah.matched_prefix_ions));
          }
          if (annotation_matched_suffix_ions)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS, static_cast<int>(ah.matched_suffix_ions));
          }

          ph.setMetaValue(Constants::UserParam::DELTA_SCORE, delta_scores[scan_index]);

          if (annotation_hyperscore_zscore)
          {
            ph.setMetaValue(Constants::UserParam::HYPERSCORE_ZSCORE, hyperscore_zscores[scan_index]);
          }
          if (annotation_ln_num_candidates)
          {
            ph.setMetaValue(Constants::UserParam::LN_NUM_CANDIDATES, ln_num_candidates[scan_index]);
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
                  // Extract ordinal: "b5", "y3-H2O1+", "c12++" -> 5, 3, 12
                  Size pos = 1;
                  while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
                  if (pos > 1)
                  {
                    int ordinal = StringUtils::toInt32(StringUtils::substr(name, 1, pos - 1));
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
              ph.setMetaValue(Constants::UserParam::MATCHED_ION_CURRENT, matched_ion_current);
            }

            if (annotation_matched_ion_current_fraction)
            {
              ph.setMetaValue(Constants::UserParam::MATCHED_ION_CURRENT_FRACTION,
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
                ph.setMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE, std::max(longest_prefix, longest_suffix));
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
                ph.setMetaValue(Constants::UserParam::COMPLEMENTARY_IONS_FRACTION, complementary_fraction);
              }
            }
          }


          // Add isotope error metavalue (always; exposed as Percolator feature)
          ph.setMetaValue(Constants::UserParam::ISOTOPE_ERROR, ah.isotope_error);

          // Add delta mass metavalue for open search
          if (isOpenSearchMode_())
          {
            ph.setMetaValue("DeltaMass", ah.delta_mass);
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
                           << " top_isotope_error=" << (int)top_hit.getMetaValue(Constants::UserParam::ISOTOPE_ERROR)
                           << std::endl;
        }
#pragma omp critical (peptide_ids_access)
        {
          //clang-tidy: seems to be a false-positive in combination with omp
          peptide_ids.push_back(std::move(pi));
        }
      }
    }

#ifdef _OPENMP
    // we need to sort the peptide_ids by scan_index in order to have the same output in the idXML-file
    if (omp_get_max_threads() > 1)
    {
      std::sort(peptide_ids.begin(), peptide_ids.end(), [](const PeptideIdentification& a, const PeptideIdentification& b)
      {
        return a.getMetaValue("scan_index") < b.getMetaValue("scan_index");
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
    feature_set.push_back(Constants::UserParam::DELTA_SCORE);
    feature_set.push_back(Constants::UserParam::ISOTOPE_ERROR);
    // note: precursor error is calculated by percolator itself
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));
    // record whether open-search mode was used
    search_parameters.setMetaValue("open_search", isOpenSearchMode_() ? "true" : "false");

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
  Param ProSEAlgorithm::fragmentIndexParameters_() const
  {
    Param p = getParameters();
    // FragmentIndex has its own boolean 'decoys' flag {true,false}; ProSE's enum
    // value (auto/generate/ignore) would fail its validation. ProSE builds the
    // decoy-augmented database itself, so the FragmentIndex must never generate.
    p.remove("decoys");
    p.setValue("decoys", "false");
    return p;
  }

  ProSEAlgorithm::DecoyStrategy_
  ProSEAlgorithm::resolveDecoyStrategy_(const std::vector<FASTAFile::FASTAEntry>& db) const
  {
    // Detect pre-existing decoys. First the common-marker heuristic
    // (DecoyHelper: "decoy", "rev", "xxx", ... as prefix or suffix), then a
    // literal fall-back to the configured decoy_prefix so custom markers
    // outside the vocabulary are still recognised.
    FASTAContainer<TFI_Vector> container(db);
    // quiet=true: a target-only database is a normal case here (auto/generate
    // then synthesise decoys), so suppress DecoyHelper's "unable to determine
    // decoy string" ERROR/WARN noise — we handle the negative result ourselves.
    const DecoyHelper::Result det = DecoyHelper::findDecoyString(container, /*quiet=*/true);

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
    std::vector<FASTAFile::FASTAEntry> db;
    db.reserve(fasta_db.size() * (strategy.generate ? 2 : 1));

    // decoys=auto reusing pre-existing decoys logs nothing otherwise; surface the auto-detected
    // marker so a rare DecoyHelper misdetection (a target DB whose accessions start with a decoy
    // affix) is diagnosable. Single emission: buildDecoyAugmentedDB_ runs once per search.
    if (decoy_mode_ == DecoyMode_::AUTO && !strategy.generate && strategy.have_decoys)
    {
      OPENMS_LOG_INFO << "[ProSE] decoys=auto: reusing existing decoys detected in the database "
                      << "(marker '" << strategy.decoy_string << "', "
                      << (strategy.is_prefix ? "prefix" : "suffix") << ")." << std::endl;
    }

    // 1. Copy targets, dropping pre-existing decoys when requested.
    for (const FASTAFile::FASTAEntry& e : fasta_db)
    {
      const bool is_existing_decoy = strategy.strip_existing &&
        (strategy.strip_is_prefix ? StringUtils::hasPrefix(e.identifier, strategy.strip_string)
                                  : StringUtils::hasSuffix(e.identifier, strategy.strip_string));
      if (!is_existing_decoy) { db.push_back(e); }
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
        if (peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE)
          e.sequence = decoy_generator.reverseProtein(AASequence::fromString(e.sequence)).toString();
        else
          e.sequence = decoy_generator.reversePeptides(AASequence::fromString(e.sequence), enzyme_).toString();
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
    SearchContext ctx;

    startProgress(0, 1, "Generate decoys...");
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    ctx.db = buildDecoyAugmentedDB_(fasta_db, strategy);
    ctx.decoy_string = strategy.decoy_string;
    ctx.decoy_is_prefix = strategy.is_prefix;
    ctx.have_decoys = strategy.have_decoys;
    endProgress();

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    Param this_params = fragmentIndexParameters_();
    ctx.fragment_index.setParameters(this_params);
    ctx.fragment_index.build(ctx.db);
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
      const TheoreticalSpectrumGenerator& spectrum_generator,
      double effective_fragment_tol,
      bool fragment_mass_tolerance_unit_ppm,
      bool open_search_mode,
      std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
      std::vector<CandidatePoolStats_>& pool_stats,
      const std::string& progress_label) const
  {
    if (pool_stats.size() != annotated_hits.size())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, pool_stats.size(),
        "Candidate-pool statistics must hold one entry per spectrum.");
    }
    startProgress(0, spectra.size(), progress_label);
    size_t count_spectra{};
    const double proton_mass_u = Constants::PROTON_MASS_U;
    // Hoisted out of the omp parallel block: clang with `default(none)` forbids
    // referencing namespace-scope constants inside the loop without explicit sharing.
    const double c13c12_massdiff_u = Constants::C13C12_MASSDIFF_U;
    const Size keep = std::max(report_top_hits_, Size(2)); // keep ≥2 for delta score

#pragma omp parallel for schedule(dynamic) default(none) shared(annotated_hits, pool_stats, count_spectra, fi, spectrum_generator, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, c13c12_massdiff_u, effective_fragment_tol, keep)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {
      #pragma omp atomic
      ++count_spectra;

      IF_MASTERTHREAD { setProgress(count_spectra); }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      FragmentIndex::SpectrumMatchesTopN top_sms;
      fi.querySpectrum(exp_spectrum, db, top_sms);

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

        // Clear peaks + data arrays (ion names / charges) before refilling for the
        // next candidate; getSpectrum appends to whatever is there.
        theo_spectrum.clear(true);
        spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, 1);
        // Note: TSG emits sorted output when add_metainfo=true (see the
        // sortByPositionPresorted() call at the tail of getSpectrum_); the extra
        // sortByPosition() pass here was a redundant O(N) scan per candidate.

        HyperScore::PSMDetail detail;
        const double& score = HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);

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
        double seq_length = (double)ah.sequence.size();
        ah.prefix_fraction = static_cast<float>(detail.matched_prefix_ions / seq_length);
        ah.suffix_fraction = static_cast<float>(detail.matched_suffix_ions / seq_length);
        ah.mean_error = static_cast<float>(detail.mean_error);
        ah.matched_prefix_ions = static_cast<uint16_t>(detail.matched_prefix_ions);
        ah.matched_suffix_ions = static_cast<uint16_t>(detail.matched_suffix_ions);
        ah.isotope_error = sms.isotope_error_;
        ah.applied_charge = sms.precursor_charge_;
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
    if (database_chunk_size_ == 0)
    {
      SearchContext ctx = prepareContext(fasta_db);
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
      ctx.fragment_index.setParameters(fragmentIndexParameters_());
      ctx.fragment_index.build(ctx.db);
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
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_, peaks_keep_n_, peaks_window_top_);

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
      cal_fi.setParameters(fragmentIndexParameters_());
      cal_fi.build(cal_db);

      CalibrationResult_ cal = runCalibrationPass_(spectra, cal_fi, cal_db);
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

    // 3. Prepare spectrum generator (once).
    TheoreticalSpectrumGenerator spectrum_generator;
    {
      Param tsg_param(spectrum_generator.getParameters());
      tsg_param.setValue("add_first_prefix_ion", "true");
      tsg_param.setValue("add_metainfo", "true");
      tsg_param.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
      tsg_param.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
      tsg_param.setValue("add_c_ions", add_c_ions_ ? "true" : "false");
      tsg_param.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
      tsg_param.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
      tsg_param.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
      spectrum_generator.setParameters(tsg_param);
    }

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
        Param fi_params = fragmentIndexParameters_();
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
      scoreSpectraAgainstIndex_(spectra, chunk_fi, chunk_db, spectrum_generator,
                                effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                                open_search_mode, annotated_hits, pool_stats,
                                "Scoring chunk " + StringUtils::toStr(chunk_idx) + "...");

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
      "" // no database filename for in-memory search
      );
    endProgress();

    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);

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
      FalseDiscoveryRate fdr;
      Param fdr_params = fdr.getParameters();
      fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
      fdr.setParameters(fdr_params);
      fdr.apply(peptide_ids);
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
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_, peaks_keep_n_, peaks_window_top_);
    endProgress();

    // Reference the prepared (decoy-augmented) database and prebuilt fragment index from the context.
    std::vector<FASTAFile::FASTAEntry>& db = ctx.db;
    FragmentIndex& fragment_index_ = ctx.fragment_index;

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
      last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db);
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

    // create spectrum generator — forward the ion-series toggles so scoring
    // uses the same ion types as the FragmentIndex (which already reads them
    // via setParameters in prepareContext).
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param(spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
    param.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
    param.setValue("add_c_ions", add_c_ions_ ? "true" : "false");
    param.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
    param.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
    param.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
    spectrum_generator.setParameters(param);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(report_top_hits_); }
    std::vector<CandidatePoolStats_> pool_stats(spectra.size());

    bool open_search_mode = open_search;

    StopWatch sw_search; sw_search.start();
    scoreSpectraAgainstIndex_(spectra, fragment_index_, db, spectrum_generator,
                              effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                              open_search_mode, annotated_hits, pool_stats,
                              "Scoring peptide models against spectra...");

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
      "" // no database filename for in-memory search
      );
    endProgress();

    annotatePeptDeepFeatures_(spectra, protein_ids, peptide_ids);
    sw_search.stop();
    last_run_stats_.seconds_search = sw_search.getClockTime();

    // reindex peptides to proteins.
    // The PeptideIndexer drops peptides whose termini do not match the configured
    // specificity, so it must agree with the search-time setting — otherwise
    // semi-specific / non-specific PSMs would be silently filtered out here.
    PeptideIndexing indexer;
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

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(db, protein_ids, peptide_ids);

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
      FalseDiscoveryRate fdr;
      Param fdr_params = fdr.getParameters();
      fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
      fdr.setParameters(fdr_params);
      fdr.apply(peptide_ids);
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
    FileHandler f;
    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.loadExperiment(in_spectra, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW});
    spectra.sortSpectra(true);

    // load FASTA
    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(in_db, fasta_db);

    // delegate to in-memory search
    ExitCodes ec = search(spectra, fasta_db, protein_ids, peptide_ids);

    if (ec != ExitCodes::EXECUTION_OK)
    {
      return ec;
    }

    // Protein inference + picked-protein FDR for single-file search.
    // Must run before decoy removal so both target and decoy proteins
    // receive aggregated scores from BPIA. Resolve the decoy strategy from the
    // same input FASTA the search used so the marker/position match.
    const DecoyStrategy_ strategy = resolveDecoyStrategy_(fasta_db);
    if (fdr_protein_ > 0.0 && strategy.have_decoys)
    {
      // Single input file = complete experiment, so picked-protein FDR is valid. Use the resolved
      // decoy marker (prefix or suffix, detected by DecoyHelper in resolveDecoyStrategy_) so the
      // shared finalization recognises the same decoys that were searched.
      applyCompleteSetProteinFDR(protein_ids, peptide_ids, strategy.decoy_string, strategy.is_prefix, fdr_protein_);
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
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        OPENMS_LOG_INFO << "[ProSE] Loading " << in_spectra_files[i] << std::endl;
        FileHandler f;
        PeakFileOptions options;
        options.clearMSLevels();
        options.addMSLevel(2);
        f.getOptions() = options;
        f.loadExperiment(in_spectra_files[i], all_spectra[i], {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW});
        all_spectra[i].sortSpectra(true);
        preprocessSpectra_(all_spectra[i], fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, deisotope_requested_, peaks_keep_n_, peaks_window_top_);
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
        cal_fi.setParameters(fragmentIndexParameters_());
        StopWatch sw_cal_idx; sw_cal_idx.start();
        cal_fi.build(cal_db);
        sw_cal_idx.stop();
        mfres.shared.seconds_index_build += sw_cal_idx.getClockTime();

        for (Size i = 0; i < in_spectra_files.size(); ++i)
        {
          OPENMS_LOG_INFO << "[ProSE] Calibration for " << in_spectra_files[i]
                          << " (strided sample, " << cal_db.size() << " proteins)" << std::endl;
          CalibrationResult_ cal = runCalibrationPass_(all_spectra[i], cal_fi, cal_db);
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

      // Prepare spectrum generator (once).
      TheoreticalSpectrumGenerator spectrum_generator;
      {
        Param tsg_param(spectrum_generator.getParameters());
        tsg_param.setValue("add_first_prefix_ion", "true");
        tsg_param.setValue("add_metainfo", "true");
        tsg_param.setValue("add_a_ions", add_a_ions_ ? "true" : "false");
        tsg_param.setValue("add_b_ions", add_b_ions_ ? "true" : "false");
        tsg_param.setValue("add_c_ions", add_c_ions_ ? "true" : "false");
        tsg_param.setValue("add_x_ions", add_x_ions_ ? "true" : "false");
        tsg_param.setValue("add_y_ions", add_y_ions_ ? "true" : "false");
        tsg_param.setValue("add_z_ions", add_z_ions_ ? "true" : "false");
        spectrum_generator.setParameters(tsg_param);
      }

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
        chunk_fi.setParameters(fragmentIndexParameters_());
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
                                    spectrum_generator, per_file_cal[i].effective_fragment_tol,
                                    fragment_mass_tolerance_unit_ppm, open_search_mode,
                                    per_file_hits[i], per_file_pool_stats[i],
                                    "  file " + StringUtils::toStr(i + 1) + " chunk " + StringUtils::toStr(chunk_idx));
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
      mfres.per_file.reserve(in_spectra_files.size());

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");
        last_mod_match_tolerance_used_ = per_file_cal[i].mod_match_tol;

        SearchResult result;
        result.is_open_search = open_search_mode;

        postProcessHits_(all_spectra[i], per_file_hits[i], per_file_pool_stats[i],
          result.protein_ids, result.peptide_ids,
          report_top_hits_, modifications_fixed_, modifications_variable_,
          peptide_missed_cleavages_,
          std::max(per_file_cal[i].effective_precursor_tol_lower,
                   per_file_cal[i].effective_precursor_tol_upper),
          per_file_cal[i].effective_fragment_tol,
          precursor_mass_tolerance_unit_, fragment_mass_tolerance_unit_,
          precursor_min_charge_, precursor_max_charge_, enzyme_, "");

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
          continue;
        }

        // PSM-level FDR only (matching non-chunked search semantics). Target+decoy PSMs are
        // filtered alike by the q-value threshold; no decoy-specific stripping happens here
        // (decoupled from decoy removal). The decoys that pass are kept because per-file results
        // feed cross-file merging and protein-level FDR. Categorical decoy removal is a
        // protein-FDR finalization step performed by ProSE.cpp. have_decoys is the marker-aware
        // result from resolveDecoyStrategy_ (prefix or suffix, generated or external).
        const bool has_decoys = strategy.have_decoys;
        // Pre-FDR stats (target/decoy counts + HyperScore distribution).
        capturePreFdrStats_(result.peptide_ids, result.stats);
        if (fdr_psm_ > 0.0 && has_decoys)
        {
          StopWatch sw_fdr; sw_fdr.start();
          FalseDiscoveryRate fdr;
          Param fdr_params = fdr.getParameters();
          fdr_params.setValue("add_decoy_peptides", "true"); // keep decoys eligible (q-value filtered, but no decoy-specific stripping)
          fdr.setParameters(fdr_params);
          fdr.apply(result.peptide_ids);
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
      StopWatch sw_idx; sw_idx.start();
      if (!full_db.empty())
      {
        // chunk_size was set but augmented DB fits in one chunk — reuse the
        // already-built decoy-augmented DB instead of re-augmenting inside
        // prepareContext.
        ctx.db = std::move(full_db);
        ctx.decoy_string = strategy.decoy_string;
        ctx.decoy_is_prefix = strategy.is_prefix;
        ctx.have_decoys = strategy.have_decoys;
        startProgress(0, 1, "Building fragment index...");
        ctx.fragment_index.setParameters(fragmentIndexParameters_());
        ctx.fragment_index.build(ctx.db);
        endProgress();
      }
      else
      {
        ctx = prepareContext(fasta_db);
      }
      sw_idx.stop();

      // Shared report stats: index built once and reused across all files.
      mfres.shared.chunked = false;
      mfres.shared.seconds_index_build = sw_idx.getClockTime();
      mfres.shared.indexed_peptides = ctx.fragment_index.getPeptides().size();
      mfres.shared.indexed_fragments = ctx.fragment_index.getNumFragments();
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

      mfres.per_file.reserve(in_spectra_files.size());

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const std::string& in_spectra = in_spectra_files[i];
        const std::string per_file_base = (i < output_base_names.size()) ? output_base_names[i] : std::string("");

        OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                        << "] Searching " << in_spectra << std::endl;

        PeakMap spectra;
        {
          FileHandler f;
          PeakFileOptions options;
          options.clearMSLevels();
          options.addMSLevel(2);
          f.getOptions() = options;
          f.loadExperiment(in_spectra, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF, FileTypes::RAW});
        }
        spectra.sortSpectra(true);

        SearchResult result;
        result.is_open_search = isOpenSearchMode_();
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
    FASTAFile().load(in_db, fasta_db);

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
  // Helper: run calibration pass on a subset of spectra
  // =====================================================================
  ProSEAlgorithm::CalibrationResult_
  ProSEAlgorithm::runCalibrationPass_(
      PeakMap& spectra,
      FragmentIndex& fragment_index,
      const std::vector<FASTAFile::FASTAEntry>& db) const
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

    // Score subset and collect errors from the best hit per spectrum
    TheoreticalSpectrumGenerator tsg;
    Param tsg_param(tsg.getParameters());
    tsg_param.setValue("add_first_prefix_ion", "true");
    tsg_param.setValue("add_metainfo", "true");
    tsg.setParameters(tsg_param);

    // Collect per-spectrum best hits with scores and errors
    struct CalHit { double score; double prec_error; double frag_error; };
    vector<CalHit> cal_hits;

    // Parallelize over the calibration subset, mirroring the main scoring loop
    // (scoreSpectraAgainstIndex_). Each iteration is independent: querySpectrum and
    // the shared TheoreticalSpectrumGenerator expose const, thread-safe methods (the
    // main loop already calls them concurrently), and every working variable below is
    // loop-local. The only cross-thread write is the push into cal_hits, guarded by a
    // critical section. Without this the calibration pass ran single-threaded: on a
    // many-core machine the 10% subset took about as long as the entire parallel
    // search, roughly doubling wall time for no extra useful work. Downstream is
    // order-independent — cal_hits is sorted by score and the error vectors are sorted
    // before quantiles — so parallel insertion order does not change the result.
#pragma omp parallel for schedule(dynamic)
    for (SignedSize si = 0; si < (SignedSize)subset_size; ++si)
    {
      const Size scan_idx = tic_index[si].second;
      const MSSpectrum& spec = spectra[scan_idx];

      FragmentIndex::SpectrumMatchesTopN top_sms;
      fragment_index.querySpectrum(spec, top_sms);

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
        tsg.getSpectrum(theo, seq, 1, 1);

        HyperScore::PSMDetail detail;
        double score = HyperScore::computeWithDetail(
            fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm, spec, theo, detail);

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
