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
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
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
#include <iomanip>
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
    defaults_.setValue("precursor:mass_tolerance_lower", 20.0,
                       "Lower-side precursor-mass tolerance (positive magnitude; effective window "
                       "is [-lower, +upper] around the precursor). "
                       "When strongly asymmetric, also review precursor:isotope_error_min.");
    defaults_.setMinFloat("precursor:mass_tolerance_lower", 0.0);
    defaults_.setValue("precursor:mass_tolerance_upper", 20.0,
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

    defaults_.setValue("fragment:mass_tolerance", 10.0, "Fragment mass tolerance");

    std::vector<std::string> fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    defaults_.setValue("fragment:mass_tolerance_unit", "ppm", "Unit of fragment m");
    defaults_.setValidStrings("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);


    defaults_.setValue("fragment:min_mz", 150, "Minimal fragment mz for database");
    defaults_.setValue("fragment:max_mz", 2000, "Maximal fragment mz for database");
    defaults_.setValue("fragment:min_ion_index", 2, "Ions with index less than or equal to this value are not added to the fragment index (use 0 to include all ions; 2 skips b1/b2/y1/y2). Low-index ions are often noisy and unreliable.");
    defaults_.setMinInt("fragment:min_ion_index", 0);

    defaults_.setSectionDescription("fragment", "Fragments (Product Ion) Options");

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

    defaults_.setValue("modifications:fixed", std::vector<std::string>{"Carbamidomethyl (C)"}, "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'");
    defaults_.setValidStrings("modifications:fixed", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable", std::vector<std::string>{"Oxidation (M)"}, "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'");
    defaults_.setValidStrings("modifications:variable", ListUtils::create<std::string>(all_mods));
    defaults_.setValue("modifications:variable_max_per_peptide", 2, "Maximum number of residues carrying a variable modification per candidate peptide");
    defaults_.setSectionDescription("modifications", "Modifications Options");

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);

    defaults_.setValue("enzyme", "Trypsin", "The enzyme used for peptide digestion.");
    defaults_.setValidStrings("enzyme", ListUtils::create<std::string>(all_enzymes));

    defaults_.setValue("decoys", "false", "Should decoys be generated?");
    defaults_.setValidStrings("decoys", {"true","false"} );

    defaults_.setValue("decoy_prefix", "DECOY_", "Accession prefix used for decoy proteins. When decoy generation is enabled (-decoys), this prefix is added to generated decoy accessions. When decoy generation is disabled, the FASTA database may already contain decoy proteins with this prefix (e.g., from DecoyDatabase) and FDR can still be applied.", {"advanced"});

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
        Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM}
      );

    defaults_.setSectionDescription("annotate", "Annotation Options");

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

    defaults_.setValue("FDR:PSM", 0.0, "Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, set to 0 to disable filtering and report all PSMs with q-values). Requires '-decoys' to be set.");
    defaults_.setMinFloat("FDR:PSM", 0.0);
    defaults_.setMaxFloat("FDR:PSM", 1.0);
    defaults_.setValue("FDR:protein", 0.0, "Filter proteins based on picked-protein FDR q-value (e.g., 0.01 = 1% protein FDR, set to 0 to disable). Applied after PSM-level FDR. Uses the picked-protein approach (Savitski et al. 2015) which pairs target and decoy proteins by accession. Requires '-decoys' to be set.");
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
      "against it before moving to the next chunk.");
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

    precursor_isotopes_ = param_.getValue("precursor:isotopes");

    fragment_mass_tolerance_ = param_.getValue("fragment:mass_tolerance");

    fragment_mass_tolerance_unit_ = param_.getValue("fragment:mass_tolerance_unit").toString();

    modifications_fixed_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:fixed"));
    set<String> fixed_unique(modifications_fixed_.begin(), modifications_fixed_.end());
    if (fixed_unique.size() != modifications_fixed_.size())
    {
      OPENMS_LOG_WARN << "Duplicate fixed modification provided. Making them unique." << endl;
      modifications_fixed_.assign(fixed_unique.begin(), fixed_unique.end());
    }    

    modifications_variable_ = ListUtils::toStringList<std::string>(param_.getValue("modifications:variable"));
    set<String> var_unique(modifications_variable_.begin(), modifications_variable_.end());
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

    decoys_ = param_.getValue("decoys") == "true";
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
  void ProSEAlgorithm::preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 20, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

#pragma omp parallel for default(none) shared(exp, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, window_mower_filter, nlargest_filter)
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      // sort by mz
      exp[exp_index].sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(exp[exp_index], 
        fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, 
        1, 3,   // min / max charge 
        false,  // keep only deisotoped
        3, 10,  // min / max isopeaks 
        true);  // convert fragment m/z to mono-charge

      // remove noise
      window_mower_filter.filterPeakSpectrum(exp[exp_index]);
      nlargest_filter.filterPeakSpectrum(exp[exp_index]);

      // sort (nlargest changes order)
      exp[exp_index].sortByPosition();
    }
  }

  void ProSEAlgorithm::postProcessHits_(const PeakMap& exp,
        std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
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
        const String& precursor_mass_tolerance_unit_ppm,
        const String& fragment_mass_tolerance_unit_ppm,
        const Int precursor_min_charge,
        const Int precursor_max_charge,
        const String& enzyme,
        const String& database_name) const
  {
    // remove all but top n scoring; compute delta_score (best - second-best)
    // before truncation so it's available even when top_hits=1. Delta is stored
    // in a side vector (one float per spectrum) to avoid adding a field to every
    // AnnotatedHit_ candidate during scoring.
    std::vector<float> delta_scores(annotated_hits.size(), 0.0f);
#pragma omp parallel for default(none) shared(annotated_hits, top_hits, delta_scores)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      auto& hits = annotated_hits[scan_index];

      // O(N) pass: find top-2 scores for delta. Leaves partial_sort unchanged.
      double best_score = 0.0, second_best_score = 0.0;
      for (const auto& h : hits)
      {
        if (h.score > best_score) { second_best_score = best_score; best_score = h.score; }
        else if (h.score > second_best_score) { second_best_score = h.score; }
      }
      delta_scores[scan_index] = static_cast<float>(best_score - second_best_score);

      // sort and keep n best elements according to score (unchanged)
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
    }

    // Alignment is needed for fragment error, fragment annotations, longest ion run, and MIC
    const bool need_alignment = annotation_fragment_error_ppm || annotation_fragment_annotations || annotation_longest_ion_run || annotation_matched_ion_current;

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
          // The alignment tolerance mirrors the search's fragment tolerance so the
          // reported FRAGMENT_ERROR_MEDIAN_PPM is not polluted by far-off spurious
          // matches — SpectrumAlignment's default is 0.3 Da absolute, which is ~30×
          // looser than a typical 20 ppm search.
          std::vector<std::pair<Size, Size>> alignment;
          MSSpectrum theoretical_spec;
          if (need_alignment)
          {
            TheoreticalSpectrumGenerator tsg;
            Param tsg_param(tsg.getParameters());
            tsg_param.setValue("add_metainfo", "true");
            tsg_param.setValue("add_first_prefix_ion", "true");
            tsg.setParameters(tsg_param);

            const int max_frag_z = (charge >= 2) ? std::min<int>(charge - 1, 2) : 1;
            tsg.getSpectrum(theoretical_spec, ah.sequence, 1, max_frag_z);
            SpectrumAlignment sa;
            Param sa_param(sa.getParameters());
            sa_param.setValue("tolerance", fragment_mass_tolerance);
            sa_param.setValue("is_relative_tolerance", fragment_mass_tolerance_unit_ppm == "ppm" ? "true" : "false");
            sa.setParameters(sa_param);
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

          // Fragment annotations, longest ion run, and MIC all iterate the alignment + ion names
          if (annotation_fragment_annotations || annotation_longest_ion_run || annotation_matched_ion_current)
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
            // Dedup guard for MIC: in ppm-alignment mode a single experimental
            // peak can match multiple theoretical peaks (e.g. b-ion and near
            // isotope), so we must sum each exp_idx at most once. Sized only
            // when MIC is actually requested.
            std::vector<char> counted_exp_peaks(
                annotation_matched_ion_current ? spec.size() : 0, 0);
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

              if (annotation_matched_ion_current && !counted_exp_peaks[exp_idx])
              {
                matched_ion_current += spec[exp_idx].getIntensity();
                counted_exp_peaks[exp_idx] = 1;
              }

              if (annotation_longest_ion_run && ion_names[theo_idx].size() >= 2)
              {
                const String& name = ion_names[theo_idx];
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
                    int ordinal = String(name.substr(1, pos - 1)).toInt();
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

            if (annotation_longest_ion_run)
            {
              // Compute longest consecutive run across prefix and suffix series
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
              int longest = std::max(longestRun(prefix_ordinals), longestRun(suffix_ordinals));
              ph.setMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE, longest);
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
          phs.push_back(ph);
        }
        pi.setHits(phs);
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
    String identifier("ProSE_" + now.get());
    protein_ids[0].setIdentifier(identifier);
    for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = database_name;
    search_parameters.charges = String(precursor_min_charge) + ":" + String(precursor_max_charge);

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
  // Build a decoy-augmented copy of a FASTA database without building a
  // FragmentIndex. Used by prepareContext() and the chunked search path.
  std::vector<FASTAFile::FASTAEntry>
  ProSEAlgorithm::buildDecoyAugmentedDB_(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db) const
  {
    std::vector<FASTAFile::FASTAEntry> db = fasta_db;
    if (decoys_)
    {
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
        e.identifier = decoy_prefix_ + e.identifier;
        db.push_back(std::move(e));
      }
      Math::RandomShuffler shuffler(42);  // fixed seed for reproducible decoy ordering across runs/files
      shuffler.portable_random_shuffle(db.begin(), db.end());
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
    ctx.db = buildDecoyAugmentedDB_(fasta_db);
    endProgress();

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    auto this_params = getParameters();
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
      const String& progress_label) const
  {
    startProgress(0, spectra.size(), progress_label);
    size_t count_spectra{};
    const double proton_mass_u = Constants::PROTON_MASS_U;
    // Hoisted out of the omp parallel block: clang with `default(none)` forbids
    // referencing namespace-scope constants inside the loop without explicit sharing.
    const double c13c12_massdiff_u = Constants::C13C12_MASSDIFF_U;

#pragma omp parallel for schedule(static) default(none) shared(annotated_hits, count_spectra, fi, spectrum_generator, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, c13c12_massdiff_u, effective_fragment_tol)
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
    auto full_db = buildDecoyAugmentedDB_(fasta_db);
    if (full_db.size() <= database_chunk_size_)
    {
      // Augmented DB fits in one chunk — use the single-context path but skip
      // prepareContext's internal decoy re-generation by building ctx inline.
      SearchContext ctx;
      ctx.db = std::move(full_db);
      ctx.release_fragment_index_after_scoring = true; // single-use ctx (M1)
      startProgress(0, 1, "Building fragment index...");
      ctx.fragment_index.setParameters(getParameters());
      ctx.fragment_index.build(ctx.db);
      endProgress();
      return search(spectra, ctx, protein_ids, peptide_ids);
    }
    return searchChunked_(spectra, full_db, protein_ids, peptide_ids);
  }

  // =====================================================================
  // Chunked search implementation. Takes a pre-built decoy-augmented DB
  // (from buildDecoyAugmentedDB_) and splits it into chunks for scoring.
  // Called from the single-file search() wrapper and the multi-file wrapper.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::searchChunked_(
      PeakMap& spectra,
      std::vector<FASTAFile::FASTAEntry>& full_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    const Size n_chunks = (full_db.size() + database_chunk_size_ - 1) / database_chunk_size_;
    OPENMS_LOG_INFO << "[ProSE] Database chunking enabled: " << full_db.size()
                    << " proteins (incl. decoys), chunk_size=" << database_chunk_size_
                    << " → " << n_chunks << " chunks." << std::endl;

    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");
    bool open_search_mode = isOpenSearchMode_();
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);

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
      cal_fi.setParameters(getParameters());
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
        Param fi_params = getParameters();
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
                                open_search_mode, annotated_hits,
                                String("Scoring chunk ") + String(chunk_idx) + "...");

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

    // 7. PeptideIndexing against the FULL database (not per-chunk).
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string", decoy_prefix_);
    param_pi.setValue("decoy_string_position", "prefix");
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
    bool has_decoys = std::any_of(full_db.begin(), full_db.end(),
        [this](const FASTAFile::FASTAEntry& e) { return e.identifier.hasPrefix(decoy_prefix_); });

    if (fdr_psm_ > 0.0 && has_decoys)
    {
      FalseDiscoveryRate fdr;
      fdr.apply(peptide_ids);
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);

      if (fdr_protein_ == 0.0)
      {
        // No protein FDR → safe to remove decoys now (matches non-chunked path)
        IDFilter::removeDecoyHits(peptide_ids);
        IDFilter::removeEmptyIdentifications(peptide_ids);
        IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);
      }
    }
    else if (fdr_psm_ > 0.0 && !has_decoys)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but no decoys found (decoy_prefix='" << decoy_prefix_
                      << "'). Provide a FASTA with decoy proteins or enable '-decoys'. Skipping FDR filtering." << std::endl;
    }

    restore_tolerances();

    logSearchDiagnostics_(spectra, protein_ids, peptide_ids);

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
    bool fragment_mass_tolerance_unit_ppm = (fragment_mass_tolerance_unit_ == "ppm");

    bool open_search = isOpenSearchMode_();
    OPENMS_LOG_INFO << "[ProSE] open_search=" << (open_search ? "true" : "false")
                    << " (precursor tolerance [-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] "
                    << precursor_mass_tolerance_unit_ << ")" << std::endl;

    startProgress(0, 1, "Filtering spectra...");
    preprocessSpectra_(spectra, fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
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
      last_calibration_result_ = runCalibrationPass_(spectra, fragment_index_, db);
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
    // calibration wiring reaches the helper regardless of search mode. The OSMA
    // call sites below read from this field instead of re-computing.
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

    bool open_search_mode = open_search;

    scoreSpectraAgainstIndex_(spectra, fragment_index_, db, spectrum_generator,
                              effective_fragment_tol, fragment_mass_tolerance_unit_ppm,
                              open_search_mode, annotated_hits,
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

    // Perform modification analysis for open search results
    if (open_search)
    {
      OPENMS_LOG_INFO << "[ProSE] Performing open search modification analysis..." << std::endl;
      startProgress(0, 1, "Analyzing modification patterns...");

      OpenSearchModificationAnalysis mod_analyzer;
      // Read from last_mod_match_tolerance_used_ (captured earlier post-calibration,
      // pre-restoration) rather than re-computing: by now the restore_fi_params()
      // call at the end of search() has already reset members to user-configured
      // values, so computeModMatchTolerance_() would return the wrong value here.
      auto modification_summaries = mod_analyzer.analyzeModifications(
        peptide_ids,
        last_mod_match_tolerance_used_,
        precursor_mass_tolerance_unit_ == "ppm",
        false, // no smoothing
        ""     // no output file for in-memory search
      );

      OPENMS_LOG_INFO << "[ProSE] Found " << modification_summaries.size()
                      << " modification patterns in open search results." << std::endl;

      endProgress();
    }

    // reindex peptides to proteins.
    // The PeptideIndexer drops peptides whose termini do not match the configured
    // specificity, so it must agree with the search-time setting — otherwise
    // semi-specific / non-specific PSMs would be silently filtered out here.
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string", decoy_prefix_);
    param_pi.setValue("decoy_string_position", "prefix");
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

    // PSM-level FDR filtering.
    // Decoys may be generated internally (decoys_=true) or present in the
    // input FASTA (external, identified by decoy_prefix_).
    bool has_decoys = std::any_of(db.begin(), db.end(),
        [this](const FASTAFile::FASTAEntry& e) { return e.identifier.hasPrefix(decoy_prefix_); });

    if (fdr_psm_ > 0.0 && has_decoys)
    {
      FalseDiscoveryRate fdr;
      fdr.apply(peptide_ids);
      IDFilter::filterHitsByScore(peptide_ids, fdr_psm_);

      if (fdr_protein_ == 0.0)
      {
        // No protein FDR → safe to remove decoys now
        IDFilter::removeDecoyHits(peptide_ids);
        IDFilter::removeEmptyIdentifications(peptide_ids);
        IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);
      }
    }
    else if (fdr_psm_ > 0.0 && !has_decoys)
    {
      OPENMS_LOG_WARN << "FDR:PSM is set but no decoys found (decoy_prefix='" << decoy_prefix_
                      << "'). Provide a FASTA with decoy proteins or enable '-decoys'. Skipping FDR filtering." << endl;
    }

    restore_fi_params();

    logSearchDiagnostics_(spectra, protein_ids, peptide_ids);

    return ExitCodes::EXECUTION_OK;
  }

  // =====================================================================
  // File-based search: thin I/O wrapper that delegates to in-memory search
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      const String& in_spectra, const String& in_db,
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
    f.loadExperiment(in_spectra, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF});
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
    // receive aggregated scores from BPIA.
    bool has_decoys_single = std::any_of(fasta_db.begin(), fasta_db.end(),
        [this](const FASTAFile::FASTAEntry& e) { return e.identifier.hasPrefix(decoy_prefix_); });
    if (fdr_protein_ > 0.0 && (decoys_ || has_decoys_single))
    {
      BasicProteinInferenceAlgorithm bpia;
      bpia.run(peptide_ids, protein_ids);

      FalseDiscoveryRate fdr;
      fdr.applyPickedProteinFDR(protein_ids[0], decoy_prefix_, true);
      IDFilter::filterHitsByScore(protein_ids, fdr_protein_);

      // Now safe to remove decoys (deferred from search() above)
      IDFilter::removeDecoyHits(peptide_ids);
      IDFilter::removeEmptyIdentifications(peptide_ids);
      IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);
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
      const String& output_base_name) const
  {
    SearchResult result;
    result.is_open_search = isOpenSearchMode_();

    result.exit_code = search(spectra, fasta_db, result.protein_ids, result.peptide_ids);

    if (result.exit_code != ExitCodes::EXECUTION_OK)
    {
      return result;
    }

    if (result.is_open_search)
    {
      OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for open search results..." << std::endl;

      OpenSearchModificationAnalysis mod_analyzer;

      String output_file = "";
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

      logModificationAnalysisSummary_(result, output_base_name);
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
  // once on the pooled set.
  // =====================================================================
  ProSEAlgorithm::MultiFileSearchResult
  ProSEAlgorithm::searchWithModificationAnalysis(
      const std::vector<String>& in_spectra_files,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      const std::vector<String>& output_base_names,
      const String& aggregate_base_name) const
  {
    if (!output_base_names.empty() && output_base_names.size() != in_spectra_files.size())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "output_base_names must be empty or have exactly one entry per spectrum file (got "
        + String(output_base_names.size()) + " entries for " + String(in_spectra_files.size())
        + " spectrum files).");
    }

    MultiFileSearchResult mfres;
    mfres.aggregate.is_open_search = isOpenSearchMode_();

    if (in_spectra_files.empty())
    {
      mfres.aggregate.exit_code = ExitCodes::INPUT_FILE_EMPTY;
      return mfres;
    }

    // Decide chunking on the decoy-augmented DB size (#9180): if decoys_ doubles
    // the target DB, a 3000-protein target against chunk_size=5000 should still
    // chunk because the augmented DB is 6000 — otherwise the resulting FI would
    // exceed the user's memory budget by 2×.
    std::vector<FASTAFile::FASTAEntry> full_db;
    bool use_chunked = false;
    if (database_chunk_size_ > 0)
    {
      full_db = buildDecoyAugmentedDB_(fasta_db);
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
        f.loadExperiment(in_spectra_files[i], all_spectra[i], {FileTypes::MZML, FileTypes::BRUKER_TDF});
        all_spectra[i].sortSpectra(true);
        preprocessSpectra_(all_spectra[i], fragment_mass_tolerance_, fragment_mass_tolerance_unit_ppm);
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
        cal_fi.setParameters(getParameters());
        cal_fi.build(cal_db);

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
      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        per_file_hits[i].resize(all_spectra[i].size());
        for (auto& a : per_file_hits[i]) a.reserve(report_top_hits_);
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
        chunk_fi.setParameters(getParameters());
        chunk_fi.build(chunk_db);

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
                                    per_file_hits[i],
                                    String("  file ") + String(i + 1) + " chunk " + String(chunk_idx));
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
        const String& in_spectra = in_spectra_files[i];
        const String per_file_base = (i < output_base_names.size()) ? output_base_names[i] : String("");
        last_mod_match_tolerance_used_ = per_file_cal[i].mod_match_tol;

        SearchResult result;
        result.is_open_search = open_search_mode;

        postProcessHits_(all_spectra[i], per_file_hits[i],
          result.protein_ids, result.peptide_ids,
          report_top_hits_, modifications_fixed_, modifications_variable_,
          peptide_missed_cleavages_,
          std::max(per_file_cal[i].effective_precursor_tol_lower,
                   per_file_cal[i].effective_precursor_tol_upper),
          per_file_cal[i].effective_fragment_tol,
          precursor_mass_tolerance_unit_, fragment_mass_tolerance_unit_,
          precursor_min_charge_, precursor_max_charge_, enzyme_, "");

        PeptideIndexing indexer;
        Param param_pi = indexer.getParameters();
        param_pi.setValue("decoy_string", decoy_prefix_);
        param_pi.setValue("decoy_string_position", "prefix");
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
          mfres.per_file.push_back(std::move(result));
          continue;
        }

        // PSM-level FDR only (matching non-chunked search semantics).
        bool has_decoys = std::any_of(full_db.begin(), full_db.end(),
            [this](const FASTAFile::FASTAEntry& e) { return e.identifier.hasPrefix(decoy_prefix_); });
        if (fdr_psm_ > 0.0 && has_decoys)
        {
          FalseDiscoveryRate fdr;
          fdr.apply(result.peptide_ids);
          IDFilter::filterHitsByScore(result.peptide_ids, fdr_psm_);
          if (fdr_protein_ == 0.0)
          {
            IDFilter::removeDecoyHits(result.peptide_ids);
            IDFilter::removeEmptyIdentifications(result.peptide_ids);
            IDFilter::removeUnreferencedProteins(result.protein_ids, result.peptide_ids);
          }
        }

        result.exit_code = ExitCodes::EXECUTION_OK;

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, all_spectra[i]);

        logSearchDiagnostics_(all_spectra[i], result.protein_ids, result.peptide_ids);

        // Per-file modification analysis (uses per-file calibrated tolerance).
        if (result.is_open_search)
        {
          OpenSearchModificationAnalysis mod_analyzer;
          String output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, per_file_cal[i].mod_match_tol,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
          logModificationAnalysisSummary_(result, per_file_base);
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
      if (!full_db.empty())
      {
        // chunk_size was set but augmented DB fits in one chunk — reuse the
        // already-built decoy-augmented DB instead of re-augmenting inside
        // prepareContext.
        ctx.db = std::move(full_db);
        startProgress(0, 1, "Building fragment index...");
        ctx.fragment_index.setParameters(getParameters());
        ctx.fragment_index.build(ctx.db);
        endProgress();
      }
      else
      {
        ctx = prepareContext(fasta_db);
      }

      mfres.per_file.reserve(in_spectra_files.size());

      for (Size i = 0; i < in_spectra_files.size(); ++i)
      {
        const String& in_spectra = in_spectra_files[i];
        const String per_file_base = (i < output_base_names.size()) ? output_base_names[i] : String("");

        OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                        << "] Searching " << in_spectra << std::endl;

        PeakMap spectra;
        {
          FileHandler f;
          PeakFileOptions options;
          options.clearMSLevels();
          options.addMSLevel(2);
          f.getOptions() = options;
          f.loadExperiment(in_spectra, spectra, {FileTypes::MZML, FileTypes::BRUKER_TDF});
        }
        spectra.sortSpectra(true);

        SearchResult result;
        result.is_open_search = isOpenSearchMode_();
        result.exit_code = search(spectra, ctx, result.protein_ids, result.peptide_ids);

        if (result.exit_code != ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_WARN << "[ProSE] Search failed for " << in_spectra
                          << " (exit code " << static_cast<int>(result.exit_code) << "). Continuing." << std::endl;
          mfres.per_file.push_back(std::move(result));
          continue;
        }

        if (!result.protein_ids.empty())
          result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);

        if (result.is_open_search)
        {
          OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for " << in_spectra << std::endl;
          OpenSearchModificationAnalysis mod_analyzer;
          String output_file = per_file_base.empty() ? "" : per_file_base + "_ModificationAnalysis.idXML";
          result.modification_analysis = mod_analyzer.analyzeModificationsWithStatistics(
            result.peptide_ids, last_mod_match_tolerance_used_,
            precursor_mass_tolerance_unit_ == "ppm", false, output_file);
          logModificationAnalysisSummary_(result, per_file_base);
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

    // Merge per-file identifications into a single aggregate using
    // IDMergerAlgorithm — the canonical OpenMS pattern for cross-file protein
    // inference (mirrors ConsensusMapMergerAlgorithm::mergeAllIDRuns in
    // ProteomicsLFQ). This deduplicates ProteinHits by accession (union) and
    // remaps all PeptideIdentification identifiers to the merged run. No
    // PeptideIndexing re-run needed: per-file searches already linked
    // peptides to proteins.
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
      String agg_output_file = "";
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

      logModificationAnalysisSummary_(mfres.aggregate, aggregate_base_name);
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
      const std::vector<String>& in_spectra_files,
      const String& in_db,
      const std::vector<String>& output_base_names,
      const String& aggregate_base_name) const
  {
    // load FASTA once
    vector<FASTAFile::FASTAEntry> fasta_db;
    FASTAFile().load(in_db, fasta_db);

    MultiFileSearchResult mfres = searchWithModificationAnalysis(
      in_spectra_files, fasta_db, output_base_names, aggregate_base_name);

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
  ProSEAlgorithm::searchWithModificationAnalysis(const String& in_spectra,
                                                                  const String& in_db,
                                                                  const String& output_base_name) const
  {
    std::vector<String> in_files{in_spectra};
    std::vector<String> base_names;
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
  // Helper: log search summary statistics and per-run tolerance estimation
  // =====================================================================
  void ProSEAlgorithm::logSearchDiagnostics_(
      const PeakMap& spectra,
      const std::vector<ProteinIdentification>& /*protein_ids*/,
      const PeptideIdentificationList& peptide_ids) const
  {
    // -- Search summary --
    Size num_ms2 = std::count_if(spectra.begin(), spectra.end(),
                                 [](const MSSpectrum& s) { return s.getMSLevel() == 2; });
    Size num_identified = peptide_ids.size();

    set<String> unique_peptides;
    set<String> unique_proteins;

    // Collect per-PSM error values for tolerance estimation (top-ranked hits only)
    vector<double> precursor_errors;
    vector<double> fragment_errors;

    // Missed cleavages
    EnzymaticDigestion digestor;
    digestor.setEnzyme(ProteaseDB::getInstance()->getEnzyme(enzyme_));
    map<Size, Size> mc_counts;

    for (const auto& pid : peptide_ids)
    {
      if (pid.getHits().empty()) continue;
      const PeptideHit& top = pid.getHits().front();
      unique_peptides.insert(top.getSequence().toString());

      for (const auto& ev : top.getPeptideEvidences())
      {
        unique_proteins.insert(ev.getProteinAccession());
      }

      // Precursor error: always compute inline (cheap) so tolerance estimation
      // does not depend on the annotate:PSM setting.
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

      mc_counts[digestor.countInternalCleavageSites(top.getSequence().toUnmodifiedString())]++;
    }

    OPENMS_LOG_INFO << "\n[ProSE] ============ Search Summary ============" << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   MS2 spectra:          " << num_ms2 << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Matched spectra:      " << num_identified << std::endl;
    if (num_ms2 > 0)
    {
      OPENMS_LOG_INFO << "[ProSE]   MS2 ID rate:          "
                      << std::fixed << std::setprecision(1) << (100.0 * num_identified / num_ms2) << "%" << std::endl;
    }
    OPENMS_LOG_INFO << "[ProSE]   Unique peptides:      " << unique_peptides.size() << std::endl;
    OPENMS_LOG_INFO << "[ProSE]   Unique proteins:      " << unique_proteins.size() << std::endl;

    // Missed cleavages distribution
    if (!mc_counts.empty())
    {
      std::ostringstream mc_oss;
      for (const auto& [mc, count] : mc_counts)
      {
        if (mc_oss.tellp() > 0) mc_oss << ", ";
        mc_oss << mc << ": " << count;
      }
      OPENMS_LOG_INFO << "[ProSE]   Missed cleavages:    " << mc_oss.str() << std::endl;
    }

    // -- Per-run tolerance estimation --
    const Size min_psms_for_estimation = 10;
    if (precursor_errors.size() >= min_psms_for_estimation || fragment_errors.size() >= min_psms_for_estimation)
    {
      OPENMS_LOG_INFO << "[ProSE] -------- Tolerance Estimation --------" << std::endl;
    }

    if (precursor_errors.size() >= min_psms_for_estimation)
    {
      double med = Math::median(precursor_errors.begin(), precursor_errors.end());
      double mad = Math::MAD(precursor_errors.begin(), precursor_errors.end(), med);
      double recommended = std::ceil(med + 3.0 * mad);
      OPENMS_LOG_INFO << "[ProSE]   Precursor error: median=" << std::fixed << std::setprecision(2) << med
                      << " ppm, MAD=" << mad << " ppm"
                      << " -> recommended: " << static_cast<int>(recommended) << " ppm" << std::endl;
    }

    if (fragment_errors.size() >= min_psms_for_estimation)
    {
      double med = Math::median(fragment_errors.begin(), fragment_errors.end());
      double mad = Math::MAD(fragment_errors.begin(), fragment_errors.end(), med);
      double recommended = std::ceil(med + 3.0 * mad);
      OPENMS_LOG_INFO << "[ProSE]   Fragment error:  median=" << std::fixed << std::setprecision(2) << med
                      << " ppm, MAD=" << mad << " ppm"
                      << " -> recommended: " << static_cast<int>(recommended) << " ppm" << std::endl;
    }

    OPENMS_LOG_INFO << "[ProSE]   (configured: precursor=[-" << precursor_mass_tolerance_lower_
                    << ", +" << precursor_mass_tolerance_upper_ << "] " << precursor_mass_tolerance_unit_
                    << ", fragment=" << fragment_mass_tolerance_ << " " << fragment_mass_tolerance_unit_ << ")" << std::endl;
    OPENMS_LOG_INFO << "[ProSE] ============================================\n" << std::endl;
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

    for (Size si = 0; si < subset_size; ++si)
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

      for (const auto& sms : top_sms.hits_)
      {
        AASequence seq = fragment_index.reconstructModifiedSequence(
            fragment_index.getPeptides()[sms.peptide_idx_], db);
        PeakSpectrum theo;
        tsg.getSpectrum(theo, seq, 1, 1);
        theo.sortByPosition();

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
  // Helper: log modification analysis summary
  // =====================================================================
  void ProSEAlgorithm::logModificationAnalysisSummary_(
      const SearchResult& result,
      const String& output_base_name) const
  {
    OPENMS_LOG_INFO << "[ProSE] ============================================" << std::endl;
    OPENMS_LOG_INFO << "[ProSE] MODIFICATION DISCOVERY SUMMARY" << std::endl;
    OPENMS_LOG_INFO << "[ProSE] ============================================" << std::endl;

    const auto& dm_stats = result.modification_analysis.delta_mass_stats;
    OPENMS_LOG_INFO << "[ProSE] Delta Mass Analysis:" << std::endl;
    OPENMS_LOG_INFO << "  Total PSMs analyzed: " << dm_stats.total_psms << std::endl;
    OPENMS_LOG_INFO << "  Modified PSMs: " << dm_stats.modified_psms
                    << " (" << (dm_stats.total_psms > 0 ? (100.0 * dm_stats.modified_psms / dm_stats.total_psms) : 0.0)
                    << "%)" << std::endl;
    OPENMS_LOG_INFO << "  Unmodified PSMs: " << dm_stats.unmodified_psms << std::endl;
    OPENMS_LOG_INFO << "  Mean delta mass: " << dm_stats.mean_delta_mass << " Da" << std::endl;
    OPENMS_LOG_INFO << "  Median delta mass: " << dm_stats.median_delta_mass << " Da" << std::endl;
    OPENMS_LOG_INFO << "  Unique delta mass bins: " << dm_stats.entries.size() << std::endl;

    const auto& ptm_stats = result.modification_analysis.ptm_stats;
    OPENMS_LOG_INFO << "[ProSE] PTM Analysis:" << std::endl;
    OPENMS_LOG_INFO << "  PSMs with known PTMs: " << ptm_stats.total_modified_psms << std::endl;
    OPENMS_LOG_INFO << "  PSMs with unknown modifications: " << ptm_stats.unknown_modification_psms << std::endl;
    OPENMS_LOG_INFO << "  Unique PTMs identified: " << ptm_stats.num_unique_modifications << std::endl;

    if (!ptm_stats.entries.empty())
    {
      OPENMS_LOG_INFO << "[ProSE] Top PTMs Discovered:" << std::endl;
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;
      OPENMS_LOG_INFO << "  Rank | Name                            | Count | %     | Mass (Da)" << std::endl;
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;

      size_t rank = 1;
      for (const auto& ptm : ptm_stats.entries)
      {
        if (rank > 15) break;
        String name = ptm.name;
        if (name.size() > 30) name = name.substr(0, 27) + "...";

        std::ostringstream oss;
        oss << "  " << std::setw(4) << rank++ << " | "
            << std::setw(31) << std::left << name << " | "
            << std::setw(5) << std::right << ptm.count << " | "
            << std::setw(5) << std::fixed << std::setprecision(1) << ptm.percentage << " | "
            << std::setw(9) << std::fixed << std::setprecision(4) << ptm.theoretical_mass;
        OPENMS_LOG_INFO << oss.str() << std::endl;
      }
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;
    }

    std::vector<OpenSearchModificationAnalysis::DeltaMassEntry> unknown_dm;
    for (const auto& entry : dm_stats.entries)
    {
      if (!entry.is_known_modification && entry.count >= 5)
      {
        unknown_dm.push_back(entry);
      }
    }

    if (!unknown_dm.empty())
    {
      OPENMS_LOG_INFO << "[ProSE] Top Unknown Delta Masses (potential novel PTMs):" << std::endl;
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;
      OPENMS_LOG_INFO << "  Rank | Delta Mass (Da) | Count | Unique Peptides" << std::endl;
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;

      size_t rank = 1;
      for (const auto& dm : unknown_dm)
      {
        if (rank > 10) break;
        std::ostringstream oss;
        oss << "  " << std::setw(4) << rank++ << " | "
            << std::setw(15) << std::fixed << std::setprecision(4) << dm.delta_mass << " | "
            << std::setw(5) << dm.count << " | "
            << std::setw(15) << dm.unique_peptides;
        OPENMS_LOG_INFO << oss.str() << std::endl;
      }
      OPENMS_LOG_INFO << "  ----------------------------------------------------------------" << std::endl;
    }

    OPENMS_LOG_INFO << "[ProSE] ============================================" << std::endl;

    if (!output_base_name.empty())
    {
      OPENMS_LOG_INFO << "[ProSE] Statistics tables written to:" << std::endl;
      OPENMS_LOG_INFO << "  - " << output_base_name << "_ModificationAnalysis_DeltaMassStats.tsv" << std::endl;
      OPENMS_LOG_INFO << "  - " << output_base_name << "_ModificationAnalysis_PTMStats.tsv" << std::endl;
    }
  }

} // namespace OpenMS
