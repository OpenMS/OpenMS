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
#include <sstream>

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
        Constants::UserParam::MATCHED_B_IONS,
        Constants::UserParam::MATCHED_Y_IONS,
        Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE,
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
    // remove all but top n scoring TODO: use two parameters to distinguish between number of reported peptides and number of pre-scored peptides
#pragma omp parallel for default(none) shared(annotated_hits, top_hits)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      Size topn = top_hits > annotated_hits[scan_index].size() ? annotated_hits[scan_index].size() : top_hits;
      std::partial_sort(annotated_hits[scan_index].begin(), annotated_hits[scan_index].begin() + topn, annotated_hits[scan_index].end(), AnnotatedHit_::hasBetterScore);
      annotated_hits[scan_index].resize(topn);
      annotated_hits[scan_index].shrink_to_fit();
    }

    bool annotation_precursor_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_fragment_error_ppm = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM) != annotate_psm_.end();
    bool annotation_prefix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_suffix_fraction = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION) != annotate_psm_.end();
    bool annotation_num_matched_peaks = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::NUM_MATCHED_PEAKS) != annotate_psm_.end();
    bool annotation_matched_b_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_B_IONS) != annotate_psm_.end();
    bool annotation_matched_y_ions = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::MATCHED_Y_IONS) != annotate_psm_.end();
    bool annotation_longest_ion_run = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE) != annotate_psm_.end();
    bool annotation_fragment_annotations = std::find(annotate_psm_.begin(), annotate_psm_.end(), Constants::UserParam::FRAGMENT_ANNOTATION_USERPARAM) != annotate_psm_.end();

    // "ALL" adds all annotations
    if (std::find(annotate_psm_.begin(), annotate_psm_.end(), "ALL") != annotate_psm_.end())
    {
      annotation_precursor_error_ppm = true;
      annotation_fragment_error_ppm = true;
      annotation_prefix_fraction = true;
      annotation_suffix_fraction = true;
      annotation_num_matched_peaks = true;
      annotation_matched_b_ions = true;
      annotation_matched_y_ions = true;
      annotation_longest_ion_run = true;
      annotation_fragment_annotations = true;
    }

    // Alignment is needed for fragment error, fragment annotations, and longest ion run
    const bool need_alignment = annotation_fragment_error_ppm || annotation_fragment_annotations || annotation_longest_ion_run;

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

          // Generate theoretical spectrum + alignment for annotations that need it
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
            double theo_mz = ah.sequence.getMZ(used_charge);
            double ppm_difference = Math::getPPM(mz, theo_mz);
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
            ph.setMetaValue(Constants::UserParam::NUM_MATCHED_PEAKS, static_cast<int>(ah.matched_b_ions + ah.matched_y_ions));
          }
          if (annotation_matched_b_ions)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_B_IONS, static_cast<int>(ah.matched_b_ions));
          }
          if (annotation_matched_y_ions)
          {
            ph.setMetaValue(Constants::UserParam::MATCHED_Y_IONS, static_cast<int>(ah.matched_y_ions));
          }

          // Fragment annotations and longest ion run (both need alignment + ion names)
          if (annotation_fragment_annotations || annotation_longest_ion_run)
          {
            const auto& ion_names = theoretical_spec.getStringDataArrays()[0];
            const auto& ion_charges = theoretical_spec.getIntegerDataArrays()[0];

            // Build PeakAnnotation vector + collect ion ordinals for longest run
            std::vector<PeptideHit::PeakAnnotation> peak_annotations;
            std::vector<int> b_ordinals, y_ordinals;
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

              if (annotation_longest_ion_run)
              {
                const String& name = ion_names[theo_idx];
                if (name.size() >= 2 && (name[0] == 'b' || name[0] == 'y'))
                {
                  // Extract ordinal: "b5", "y3-H2O1+", "b12++" -> 5, 3, 12
                  Size pos = 1;
                  while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
                  if (pos > 1)
                  {
                    int ordinal = String(name.substr(1, pos - 1)).toInt();
                    if (name[0] == 'b') b_ordinals.push_back(ordinal);
                    else y_ordinals.push_back(ordinal);
                  }
                }
              }
            }

            if (annotation_fragment_annotations)
            {
              ph.setPeakAnnotations(std::move(peak_annotations));
            }

            if (annotation_longest_ion_run)
            {
              // Compute longest consecutive run across b and y series
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
              int longest = std::max(longestRun(b_ordinals), longestRun(y_ordinals));
              ph.setMetaValue(Constants::UserParam::LONGEST_PEPTIDE_ION_SEQUENCE, longest);
            }
          }

          // Add isotope error metavalue (always)
          ph.setMetaValue("isotope_error", ah.isotope_error);

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
                           << " top_isotope_error=" << (int)top_hit.getMetaValue("isotope_error")
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
  ProSEAlgorithm::SearchContext
  ProSEAlgorithm::prepareContext(
      const std::vector<FASTAFile::FASTAEntry>& fasta_db) const
  {
    SearchContext ctx;
    // Always work with a mutable local copy (PeptideIndexing::run requires non-const ref)
    ctx.db = fasta_db;

    if (decoys_)
    {
      startProgress(0, 1, "Generate decoys...");
      DecoyGenerator decoy_generator;

      const size_t old_size = ctx.db.size();
      ctx.db.reserve(ctx.db.size() * 2);
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = ctx.db[i];
        // Non-specific search: plain reverse (no enzyme boundaries to preserve).
        // Enzyme-specific search: reverse within enzymatic peptide boundaries.
        if (peptide_enzyme_specificity_ == EnzymaticDigestion::SPEC_NONE)
        {
          e.sequence = decoy_generator.reverseProtein(AASequence::fromString(e.sequence)).toString();
        }
        else
        {
          e.sequence = decoy_generator.reversePeptides(AASequence::fromString(e.sequence), enzyme_).toString();
        }
        e.identifier = decoy_prefix_ + e.identifier;
        ctx.db.push_back(std::move(e));
      }
      Math::RandomShuffler shuffler;
      shuffler.portable_random_shuffle(ctx.db.begin(), ctx.db.end());
      endProgress();
    }

    // build fragment index
    startProgress(0, 1, "Building fragment index...");
    auto this_params = getParameters();
    ctx.fragment_index.setParameters(this_params);
    ctx.fragment_index.build(ctx.db);
    endProgress();

    return ctx;
  }

  // =====================================================================
  // In-memory search: thin wrapper that builds a fresh SearchContext per call.
  // For repeated searches against the same database, prefer the
  // context-taking overload below.
  // =====================================================================
  ProSEAlgorithm::ExitCodes ProSEAlgorithm::search(
      PeakMap& spectra,
      const std::vector<FASTAFile::FASTAEntry>& fasta_db,
      vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids) const
  {
    SearchContext ctx = prepareContext(fasta_db);
    return search(spectra, ctx, protein_ids, peptide_ids);
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

    // create spectrum generator
    TheoreticalSpectrumGenerator spectrum_generator;
    Param param(spectrum_generator.getParameters());
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_metainfo", "true");
    spectrum_generator.setParameters(param);

    // preallocate storage for PSMs
    vector<vector<AnnotatedHit_> > annotated_hits(spectra.size(), vector<AnnotatedHit_>());
    for (auto & a : annotated_hits) { a.reserve(report_top_hits_); }

    startProgress(0, spectra.size(), "Scoring peptide models against spectra...");
    size_t count_spectra{};

    bool open_search_mode = open_search;
    const double proton_mass_u = Constants::PROTON_MASS_U;

#pragma omp parallel for schedule(static) default(none) shared(annotated_hits, count_spectra, fragment_index_, spectrum_generator, db, fragment_mass_tolerance_unit_ppm, spectra, open_search_mode, proton_mass_u, effective_fragment_tol)
    for (SignedSize scan_index = 0; scan_index < (SignedSize)spectra.size(); ++scan_index)
    {

      #pragma omp atomic
      ++count_spectra;

      IF_MASTERTHREAD
      {
        setProgress(count_spectra);
      }

      const MSSpectrum& exp_spectrum = spectra[scan_index];
      FragmentIndex::SpectrumMatchesTopN top_sms;
      fragment_index_.querySpectrum(exp_spectrum, top_sms);

      for (const auto& sms : top_sms.hits_)
      {
        const FragmentIndex::Peptide& sms_pep = fragment_index_.getPeptides()[sms.peptide_idx_];
        AASequence mod_candidate = fragment_index_.reconstructModifiedSequence(sms_pep, db);

        PeakSpectrum theo_spectrum;
        spectrum_generator.getSpectrum(theo_spectrum, mod_candidate, 1, 1);
        theo_spectrum.sortByPosition();

        HyperScore::PSMDetail detail;
        const double& score = HyperScore::computeWithDetail(effective_fragment_tol, fragment_mass_tolerance_unit_ppm, exp_spectrum, theo_spectrum, detail);

        if (score == 0)
        {
          continue;
        }

        AnnotatedHit_ ah;
        ah.sequence = std::move(mod_candidate);
        ah.score = score;
        double seq_length = (double)ah.sequence.size();
        ah.prefix_fraction = static_cast<float>(detail.matched_b_ions / seq_length);
        ah.suffix_fraction = static_cast<float>(detail.matched_y_ions / seq_length);
        ah.mean_error = static_cast<float>(detail.mean_error);
        ah.matched_b_ions = static_cast<uint16_t>(detail.matched_b_ions);
        ah.matched_y_ions = static_cast<uint16_t>(detail.matched_y_ions);

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

    // Build the (decoy-augmented) database and FragmentIndex exactly once.
    SearchContext ctx = prepareContext(fasta_db);

    mfres.per_file.reserve(in_spectra_files.size());

    for (Size i = 0; i < in_spectra_files.size(); ++i)
    {
      const String& in_spectra = in_spectra_files[i];
      const String per_file_base = (i < output_base_names.size()) ? output_base_names[i] : String("");

      OPENMS_LOG_INFO << "[ProSE] [" << (i + 1) << "/" << in_spectra_files.size()
                      << "] Searching " << in_spectra << std::endl;

      // load MS2 map
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

      // patch per-file metadata
      if (!result.protein_ids.empty())
      {
        result.protein_ids[0].setPrimaryMSRunPath({in_spectra}, spectra);
      }

      // per-file modification analysis (only meaningful in open-search mode)
      if (result.is_open_search)
      {
        OPENMS_LOG_INFO << "[ProSE] Running detailed modification analysis for " << in_spectra << std::endl;

        OpenSearchModificationAnalysis mod_analyzer;
        String output_file = "";
        if (!per_file_base.empty())
        {
          output_file = per_file_base + "_ModificationAnalysis.idXML";
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

        logModificationAnalysisSummary_(result, per_file_base);
      }
      else
      {
        OPENMS_LOG_INFO << "[ProSE] Closed search mode - per-file modification analysis skipped" << std::endl;
      }

      mfres.per_file.push_back(std::move(result));
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

      // Compute precursor error (signed)
      double exp_mz = spec.getPrecursors()[0].getMZ();
      double theo_mz = best_seq.getMZ(best_charge);
      double corrected_exp_mz = exp_mz - static_cast<double>(best_isotope_error)
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

    // Signed median gives the calibration bias; spread is median(|e - shift|) + 3*MAD(|e - shift|),
    // i.e., the same "typical residual + 3*scale" shape as the pre-refactor formula, just applied
    // to residuals around the signed center instead of raw errors. Zero-centered distributions
    // get identical spread; biased distributions correctly separate shift from spread.
    const double prec_shift = Math::median(precursor_errors.begin(), precursor_errors.end());

    std::vector<double> residuals;
    residuals.reserve(precursor_errors.size());
    for (double e : precursor_errors) { residuals.push_back(std::abs(e - prec_shift)); }
    const double res_median = Math::median(residuals.begin(), residuals.end());
    const double res_mad    = Math::MAD(residuals.begin(), residuals.end(), res_median);

    result.precursor_shift  = prec_shift;
    result.precursor_spread = std::max(min_tolerance, res_median + 3.0 * res_mad);
    // Strict `>`: at |shift| == spread the signed window is exactly [0, 2*spread] or
    // [-2*spread, 0] — both expressible in the positive-magnitude schema (one side == 0
    // is legal and a well-defined half-line window). Only |shift| strictly greater
    // than spread requires discarding the calibration result.
    result.extreme_bias     = std::abs(prec_shift) > result.precursor_spread;

    if (!result.extreme_bias)
    {
      // Map signed window [shift - spread, shift + spread] to the positive-magnitude
      // (lower, upper) schema. Per the convention established above (lines ~1615-1617),
      // a (lower, upper) pair means theoretical ∈ [observed - lower, observed + upper],
      // equivalently the signed error e = observed - theoretical lies in [-upper, +lower].
      // Matching [shift - spread, shift + spread] against [-upper, +lower] gives:
      //   -cal_upper = shift - spread  ->  cal_upper = spread - shift
      //   +cal_lower = shift + spread  ->  cal_lower = spread + shift
      // Both are non-negative when |shift| < spread.
      const double cal_lower_raw = result.precursor_spread + prec_shift;
      const double cal_upper_raw = result.precursor_spread - prec_shift;
      // Only tighten - cap against user-configured bounds.
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
    OPENMS_LOG_INFO << "[ProSE]   Precursor signed: shift=" << std::fixed << std::setprecision(2) << prec_shift
                    << ", residual median=" << res_median << ", residual MAD=" << res_mad
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
