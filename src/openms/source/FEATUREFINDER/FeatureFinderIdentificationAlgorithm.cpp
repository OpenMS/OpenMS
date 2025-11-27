// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/FEATUREFINDER/FFIDAlgoExternalIDHandler.h>
#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>

#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ML/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/MATH/MathFunctions.h>


#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace OpenMS::Internal;

namespace OpenMS
{

  FeatureFinderIdentificationAlgorithm::FeatureFinderIdentificationAlgorithm() :
    DefaultParamHandler("FeatureFinderIdentificationAlgorithm")
  {
    std::vector<std::string> output_file_tags;
    output_file_tags.emplace_back("output file");

    defaults_.setValue("candidates_out", "", "Optional output file with feature candidates.", output_file_tags);

    defaults_.setValue("debug", 0, "Debug level for feature detection.", {"advanced"});
    defaults_.setMinInt("debug", 0);

    defaults_.setValue("extract:batch_size", 5000, "Nr of peptides used in each batch of chromatogram extraction."
                         " Smaller values decrease memory usage but increase runtime.");
    defaults_.setMinInt("extract:batch_size", 1);
    defaults_.setValue("extract:mz_window", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)");
    defaults_.setMinFloat("extract:mz_window", 0.0);
    defaults_.setValue(
      "extract:IM_window",
      0.06,
      "Ion mobility (IM) window for chromatogram extraction in the IM dimension. "
      "Set to 0.0 to disable IM filtering (even if data contains IM information). "
      "The window is applied as +/- IM_window/2 around the median IM value of identified peptides. "
      "This parameter is automatically ignored if the input data does not contain IM information "
      "(determined via IMTypes::determineIMFormat). "
      "Currently only concatenated IM format is supported. "
      "Typical values: 0.05-0.10 for TIMS data (1/K0 units), 3-5 for FAIMS data (compensation voltage). "
      "Note: IM values are calculated per peptide/charge/RT-region, using the median of all identifications "
      "in that region for robustness. The median, min, and max IM values are propagated to output features "
      "as meta-values (IM_median, IM_min, IM_max) for quality control.");
    defaults_.setMinFloat("extract:IM_window", 0.0);

    defaults_.setValue("extract:n_isotopes", 2, "Number of isotopes to include in each peptide assay.");
    defaults_.setMinInt("extract:n_isotopes", 2);
    defaults_.setValue(
      "extract:isotope_pmin",
      0.0, 
      "Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.",
      {"advanced"});
    defaults_.setMinFloat("extract:isotope_pmin", 0.0);
    defaults_.setMaxFloat("extract:isotope_pmin", 1.0);
    defaults_.setValue(
      "extract:rt_quantile", 
      0.95, 
      "Quantile of the RT deviations between aligned internal and external IDs to use for scaling the RT extraction window",
      {"advanced"});
    defaults_.setMinFloat("extract:rt_quantile", 0.0);
    defaults_.setMaxFloat("extract:rt_quantile", 1.0);

    defaults_.setValue(
      "extract:rt_window", 
      0.0, 
      "RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.",
      {"advanced"});
    defaults_.setMinFloat("extract:rt_window", 0.0);

    defaults_.setSectionDescription("extract", "Parameters for ion chromatogram extraction");

    defaults_.setValue("detect:peak_width", 60.0, "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.");
    defaults_.setMinFloat("detect:peak_width", 0.0);
    defaults_.setValue(
      "detect:min_peak_width", 
      0.2, 
      "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.",
      {"advanced"});
    defaults_.setMinFloat("detect:min_peak_width", 0.0);

    defaults_.setValue(
      "detect:signal_to_noise", 
      0.8, 
      "Signal-to-noise threshold for OpenSWATH feature detection",
       {"advanced"});
    defaults_.setMinFloat("detect:signal_to_noise", 0.1);
    defaults_.setValue("detect:mapping_tolerance", 0.0, "RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.");
    defaults_.setMinFloat("detect:mapping_tolerance", 0.0);

    defaults_.setSectionDescription("detect", "Parameters for detecting features in extracted ion chromatograms");

    // parameters for SVM classification:
    defaults_.setValue("svm:samples", 0, "Number of observations to use for training ('0' for all)");
    defaults_.setMinInt("svm:samples", 0);
    defaults_.setValue("svm:no_selection", "false", "By default, roughly the same number of positive and negative observations, with the same intensity distribution, are selected for training. This aims to reduce biases, but also reduces the amount of training data. Set this flag to skip this procedure and consider all available observations (subject to 'svm:samples').");
    defaults_.setValidStrings("svm:no_selection", {"true","false"});
    defaults_.setValue("svm:xval_out", "", "Output file: SVM cross-validation (parameter optimization) results", output_file_tags);
    defaults_.setValidStrings("svm:xval_out", {"csv"});
    defaults_.insert("svm:", SimpleSVM().getParameters());

    defaults_.setValue("quantify_decoys", "false", "Whether decoy peptides should be quantified (true) or skipped (false).");
    defaults_.setValidStrings("quantify_decoys", {"true","false"});
    defaults_.setValue("min_psm_cutoff", "none", "Minimum score for the best PSM of a spectrum to be used as seed. Use 'none' for no cutoff.");

    defaults_.setValue("add_mass_offset_peptides", 0.0, "If for every peptide (or seed) also an offset peptide is extracted (true). Can be used to downstream to determine MBR false transfer rates. (0.0 = disabled)");
    defaults_.setMinFloat("add_mass_offset_peptides", 0.0);

    // available scores: initialPeakQuality,total_xic,peak_apices_sum,var_xcorr_coelution,var_xcorr_coelution_weighted,var_xcorr_shape,var_xcorr_shape_weighted,var_library_corr,var_library_rmsd,var_library_sangle,var_library_rootmeansquare,var_library_manhattan,var_library_dotprod,var_intensity_score,nr_peaks,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,var_massdev_score_weighted,var_bseries_score,var_yseries_score,var_dotprod_score,var_manhatt_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score
    // exclude some redundant/uninformative scores:
    // @TODO: intensity bias introduced by "peak_apices_sum"?
    // names of scores to use as SVM features
    String score_metavalues = "peak_apices_sum,var_xcorr_coelution,var_xcorr_shape,var_library_sangle,var_intensity_score,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_ms1_isotope_correlation_score,var_ms1_isotope_overlap_score,var_massdev_score,main_var_xx_swath_prelim_score";

    defaults_.setValue(
      "svm:predictors", 
      score_metavalues, 
      "Names of OpenSWATH scores to use as predictors for the SVM (comma-separated list)",
      {"advanced"});

    defaults_.setValue(
      "svm:min_prob", 
      0.0, 
      "Minimum probability of correctness, as predicted by the SVM, required to retain a feature candidate",
      {"advanced"});
    defaults_.setMinFloat("svm:min_prob", 0.0);
    defaults_.setMaxFloat("svm:min_prob", 1.0);

    defaults_.setSectionDescription("svm", "Parameters for scoring features using a support vector machine (SVM)");

    // parameters for model fitting (via ElutionModelFitter):
    std::vector<std::string> models = {"symmetric","asymmetric","none"};
    defaults_.setValue("model:type", models[0], "Type of elution model to fit to features");
    defaults_.setValidStrings("model:type", models);
    defaults_.insert("model:", ElutionModelFitter().getParameters()); // copy parameters
    defaults_.remove("model:asymmetric");

    defaults_.setSectionDescription("model", "Parameters for fitting elution models to features");

    defaults_.setValue("EMGScoring:max_iteration", 100, "Maximum number of iterations for EMG fitting.");
    defaults_.setMinInt("EMGScoring:max_iteration", 1);
    defaults_.setValue("EMGScoring:init_mom", "false", "Alternative initial parameters for fitting through method of moments.");
    defaults_.setValidStrings("EMGScoring:init_mom", {"true","false"});

    defaults_.setSectionDescription("EMGScoring", "Parameters for fitting exp. mod. Gaussians to mass traces.");

    defaultsToParam_();
  }

  PeakMap& FeatureFinderIdentificationAlgorithm::getMSData()
  {
    return ms_data_;
  }

  const PeakMap& FeatureFinderIdentificationAlgorithm::getMSData() const
  {
    return ms_data_;
  }

  void FeatureFinderIdentificationAlgorithm::setMSData(const PeakMap& ms_data)
  {
    ms_data_ = ms_data; 
    
    vector<MSSpectrum>& specs = ms_data_.getSpectra();

    // keep only MS1
    specs.erase(
      std::remove_if(specs.begin(), specs.end(),
        [](const MSSpectrum & s) { return s.getMSLevel() != 1; }),
      specs.end());    
  }

  void FeatureFinderIdentificationAlgorithm::setMSData(PeakMap&& ms_data)
  {
    ms_data_ = std::move(ms_data); 
    
    vector<MSSpectrum>& specs = ms_data_.getSpectra();

    // keep only MS1
    specs.erase(
      std::remove_if(specs.begin(), specs.end(),
        [](const MSSpectrum & s) { return s.getMSLevel() != 1; }),
      specs.end());    
  }

  PeakMap& FeatureFinderIdentificationAlgorithm::getChromatograms()
  {
    return chrom_data_;
  }

  const PeakMap& FeatureFinderIdentificationAlgorithm::getChromatograms() const
  {
    return chrom_data_;
  }

  ProgressLogger& FeatureFinderIdentificationAlgorithm::getProgressLogger()
  {
    return prog_log_;
  }

  const ProgressLogger& FeatureFinderIdentificationAlgorithm::getProgressLogger() const
  {
    return prog_log_;
  }

  TargetedExperiment& FeatureFinderIdentificationAlgorithm::getLibrary()
  {
    return output_library_;
  }

  const TargetedExperiment& FeatureFinderIdentificationAlgorithm::getLibrary() const
  {
    return output_library_;
  }


  Size FeatureFinderIdentificationAlgorithm::addOffsetPeptides_(PeptideIdentificationList& peptides, double offset)
  {
    // WARNING: Superhack! Use unique ID to distinguish seeds from real IDs. Use a mod that will never occur to
    // make them truly unique and not be converted to an actual modification.
    const String pseudo_mod_name = String(10000);
    AASequence some_seq = AASequence::fromString("XXX[" + pseudo_mod_name + "]");

    PeptideIdentificationList offset_peptides;
    offset_peptides.reserve(peptides.size());
    Size n_added{};
    for (const auto & p : peptides) // for every peptide (or seed) we add an offset peptide
    {
      /*
      // check if already a peptide in peptide_map_ that is close in RT and MZ
      // if so don't add seed
      bool peptide_already_exists = false;
      double offset_RT = p.getRT();
      double offset_MZ = p.getMZ() + offset;
      double offset_charge = p.getHits()[0].getCharge();

      for (const auto & peptide : peptides)
      {
        double peptide_RT = peptide.getRT();
        double peptide_MZ = peptide.getMZ();

        // RT or MZ values of seed match in range -> peptide already exists -> don't add seed
        // Consider up to 5th isotopic trace (e.g., because of seed misassignment)
        double th_tolerance = mz_window_ppm_ ? mz_window_ * 1e-6 * peptide_MZ : mz_window_;
        if ((fabs(offset_RT - peptide_RT) <= seed_rt_window_ / 2.0) &&
           ((fabs(offset_MZ - peptide_MZ) <= th_tolerance) ||
             fabs(offset_MZ - (1.0/offset_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(offset_MZ - (2.0/offset_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(offset_MZ - (3.0/offset_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(offset_MZ - (4.0/offset_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(offset_MZ - (5.0/offset_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance)
            )
        {
          peptide_already_exists = true;
          break;
        }
      }

      // prevent decoys to be extracted at other target peptide
      if (!peptide_already_exists)
      {
      */
        offset_peptides.emplace_back();
        PeptideHit hit;
        hit.setCharge(p.getHits()[0].getCharge());
        hit.setSequence(some_seq);
        offset_peptides.back().getHits().push_back(std::move(hit));
        offset_peptides.back().setRT(p.getRT());
        offset_peptides.back().setMZ(p.getMZ() + offset);
        offset_peptides.back().setMetaValue("FFId_category", "internal");
        offset_peptides.back().setMetaValue("OffsetPeptide", "true");  // mark as offset peptide 
        offset_peptides.back().setMetaValue("SeedFeatureID", String(UniqueIdGenerator::getUniqueId())); // also mark as seed so we can indicate that we have a mass without sequence
      //}
    }

    for (auto & p : offset_peptides) // add offset peptides
    {
      peptides.push_back(std::move(p));
      addPeptideToMap_(peptides.back(), peptide_map_);
      n_added++;
    }
    
    return n_added;
  }

  Size FeatureFinderIdentificationAlgorithm::addSeeds_(PeptideIdentificationList& peptides, const FeatureMap& seeds)
  {
    size_t seeds_added{};
    // WARNING: Superhack! Use unique ID to distinguish seeds from real IDs. Use a mod that will never occur to
    // make them truly unique and not be converted to an actual modification.
    const String pseudo_mod_name = String(10000);
    AASequence some_seq = AASequence::fromString("XXX[" + pseudo_mod_name + "]");
    for (const Feature& feat : seeds)
    {
      // check if already a peptide in peptide_map_ that is close in RT and MZ
      // if so don't add seed
      bool peptide_already_exists = false;
      for (const auto & peptide : peptides)
      {
        double seed_RT = feat.getRT();
        double seed_MZ = feat.getMZ();
        double seed_charge = feat.getCharge();
        double peptide_RT = peptide.getRT();
        double peptide_MZ = peptide.getMZ();

        // RT or MZ values of seed match in range -> peptide already exists -> don't add seed
        // Consider up to 5th isotopic trace (e.g., because of seed misassignment)
        double th_tolerance = mz_window_ppm_ ? mz_window_ * 1e-6 * peptide_MZ : mz_window_;
        if ((fabs(seed_RT - peptide_RT) <= seed_rt_window_ / 2.0) &&
           ((fabs(seed_MZ - peptide_MZ) <= th_tolerance) ||
             fabs(seed_MZ - (1.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (2.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (3.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (4.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (5.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance)
            )
        {
          peptide_already_exists = true;
          String seq = "empty";
          int chg = 0;
          if (!peptide.getHits().empty())
          {
            seq = peptide.getHits()[0].getSequence().toString();
            chg = peptide.getHits()[0].getCharge();
          }
          OPENMS_LOG_DEBUG_NOFILE << "Skipping seed from FeatureID " << String(feat.getUniqueId()) << " with CHG: " << seed_charge << "; RT: " << seed_RT << "; MZ: " << seed_MZ <<
          " due to overlap with " << seq << "/" << chg << " at MZ: " << peptide_MZ << "; RT: " << peptide_RT << endl;

          break;
        }
      }

      if (!peptide_already_exists)
      {
        // WARNING: Superhack! Store ID generated from seed in the original input peptide
        // vector to make sure that the pointers that will be added to peptide_map_
        // stay valid for the duration of the function.
        peptides.emplace_back();
        PeptideHit seed_hit;
        seed_hit.setCharge(feat.getCharge());
        seed_hit.setSequence(some_seq);
        peptides.back().getHits().push_back(std::move(seed_hit));
        peptides.back().setRT(feat.getRT());
        peptides.back().setMZ(feat.getMZ());
        peptides.back().setMetaValue("FFId_category", "internal");
        peptides.back().setMetaValue("SeedFeatureID", String(feat.getUniqueId()));

        // Copy IM meta value from feature if present (some feature finders annotate IM)
        // If not present, the seed will be extracted across the full IM range
        if (feat.metaValueExists(Constants::UserParam::IM))
        {
          peptides.back().setMetaValue(Constants::UserParam::IM, feat.getMetaValue(Constants::UserParam::IM));
        }

        addPeptideToMap_(peptides.back(), peptide_map_);
        ++seeds_added;
      }
    }
    
    return seeds_added;
  }

  // ===== Helper functions for run() =====

  void FeatureFinderIdentificationAlgorithm::validateSVMParameters_() const
  {
    if ((svm_n_samples_ > 0) && (svm_n_samples_ < 2 * svm_n_parts_))
    {
      String msg = "Sample size of " + String(svm_n_samples_) +
        " (parameter 'svm:samples') is not enough for " + String(svm_n_parts_) +
        "-fold cross-validation (parameter 'svm:xval').";
      throw Exception::InvalidParameter(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, msg);
    }
  }

  void FeatureFinderIdentificationAlgorithm::initializeFeatureFinder_()
  {
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("EMGScoring:max_iteration", param_.getValue("EMGScoring:max_iteration"));
    params.setValue("EMGScoring:init_mom", param_.getValue("EMGScoring:init_mom"));
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    params.setValue("Scores:use_ionseries_scores", "false"); // since FFID only uses MS1 spectra, this is useless
    params.setValue("Scores:use_ms2_isotope_scores", "false"); // since FFID only uses MS1 spectra, this is useless
    params.setValue("Scores:use_ms1_correlation", "false"); // this would be redundant to the "MS2" correlation and since
    // precursor transition = first product transition, additionally biased
    params.setValue("Scores:use_ms1_mi", "false"); // same as above. On MS1 level we basically only care about the "MS1 fullscan" scores
    //TODO for MS1 level scoring there is an additional parameter add_up_spectra with which we can add up spectra
    // around the apex, to complete isotopic envelopes (and therefore make this score more robust).

    if ((elution_model_ != "none") || (!candidates_out_.empty()))
    {
      params.setValue("write_convex_hull", "true");
    }
    if (min_peak_width_ < 1.0)
    {
      min_peak_width_ *= peak_width_;
    }
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width",
                    peak_width_);
    params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width_);
    // disabling the signal-to-noise threshold (setting the parameter to zero)
    // totally breaks the OpenSWATH feature detection (no features found)!
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise",
                    signal_to_noise_);
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:method",
                    "corrected");
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false"); // disabled in OpenSWATH

    feat_finder_.setParameters(params);
    feat_finder_.setLogType(ProgressLogger::NONE);
    feat_finder_.setStrictFlag(false);
    // to use MS1 Swath scores:
    feat_finder_.setMS1Map(SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(boost::make_shared<MSExperiment>(ms_data_)));
  }

  double FeatureFinderIdentificationAlgorithm::calculateRTWindow_(double rt_uncertainty) const
  {
    if (rt_window_ != 0.0)
    {
      return rt_window_; // Already set, return it
    }

    // Calculate RT window based on other parameters and alignment quality:
    double map_tol = mapping_tolerance_;
    if (map_tol < 1.0)
    {
      map_tol *= (2 * peak_width_); // relative tolerance
    }
    double calculated_window = (rt_uncertainty + 2 * peak_width_ + map_tol) * 2;
    OPENMS_LOG_INFO << "RT window size calculated as " << calculated_window << " seconds." << endl;
    return calculated_window;
  }

  bool FeatureFinderIdentificationAlgorithm::isSeedPseudoHit_(const PeptideHit& hit)
  {
    return hit.getSequence().toUnmodifiedString().hasPrefix("XXX");
  }

  void FeatureFinderIdentificationAlgorithm::removeSeedPseudoIDs_(FeatureMap& features)
  {
    // Remove all hits with pseudo ids (seeds) from features
    for (Feature& f : features)
    {
      PeptideIdentificationList& ids = f.getPeptideIdentifications();

      // if we have peptide identifications assigned and all are annotated as OffsetPeptide,
      // we mark the feature as also an OffsetPeptide
      if (!ids.empty() && std::all_of(ids.begin(), ids.end(),
          [](const PeptideIdentification& pid) { return pid.metaValueExists("OffsetPeptide"); }))
      {
        f.setMetaValue("OffsetPeptide", "true");
      }

      // remove all seed pseudo hits (PSM details)
      for (auto& pid : ids)
      {
        std::vector<PeptideHit>& hits = pid.getHits();
        auto it = remove_if(hits.begin(), hits.end(), isSeedPseudoHit_);
        hits.erase(it, hits.end());
      }

      // remove empty PeptideIdentifications
      auto it = remove_if(ids.begin(), ids.end(),
        [](const PeptideIdentification& pid) { return pid.empty(); });
      ids.erase(it, ids.end());
    }

    // clean up unassigned PeptideIdentifications
    PeptideIdentificationList& ids = features.getUnassignedPeptideIdentifications();
    for (auto& pid : ids)
    {
      std::vector<PeptideHit>& hits = pid.getHits();
      auto it = remove_if(hits.begin(), hits.end(), isSeedPseudoHit_);
      hits.erase(it, hits.end());
    }

    // remove empty PeptideIdentifications
    auto it = remove_if(ids.begin(), ids.end(),
      [](const PeptideIdentification& pid) { return pid.empty(); });
    ids.erase(it, ids.end());
  }

  std::pair<double, double> FeatureFinderIdentificationAlgorithm::calculateRTBounds_(
    double rt_min, double rt_max) const
  {
    if (mapping_tolerance_ > 0.0)
    {
      double abs_tol = mapping_tolerance_;
      if (abs_tol < 1.0)
      {
        abs_tol *= (rt_max - rt_min);
      }
      return {rt_min - abs_tol, rt_max + abs_tol};
    }
    return {rt_min, rt_max};
  }

  // ===== End of helper functions =====

  void FeatureFinderIdentificationAlgorithm::run(
    PeptideIdentificationList peptides,
    const vector<ProteinIdentification>& proteins,
    PeptideIdentificationList peptides_ext,
    vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds,
    const String& spectra_file
    )
  {
    // Check for FAIMS data
    auto faims_groups = IMDataConverter::splitByFAIMSCV(std::move(ms_data_));
    const bool has_faims = faims_groups.size() > 1 || !std::isnan(faims_groups[0].first);

    if (!has_faims)
    {
      // Non-FAIMS data: restore ms_data_ and run directly
      ms_data_ = std::move(faims_groups[0].second);
      runSingleGroup_(peptides, proteins, peptides_ext, proteins_ext, features, seeds, spectra_file);
      return;
    }

    // FAIMS data: process each CV group separately
    OPENMS_LOG_INFO << "FAIMS data detected with " << faims_groups.size() << " compensation voltage(s)." << endl;

    // Clear combined outputs
    features.clear(true);
    chrom_data_.clear(true);
    output_library_.clear(true);
    bool first_group = true;

    for (auto& [group_cv, faims_group] : faims_groups)
    {
      OPENMS_LOG_INFO << "Processing FAIMS CV group: " << group_cv << " V" << endl;

      // Filter peptide IDs for this FAIMS CV
      PeptideIdentificationList peptides_cv = FAIMSHelper::filterPeptidesByFAIMSCV(peptides, group_cv);
      PeptideIdentificationList peptides_ext_cv = FAIMSHelper::filterPeptidesByFAIMSCV(peptides_ext, group_cv);

      if (peptides_cv.empty() && peptides_ext_cv.empty())
      {
        OPENMS_LOG_INFO << "No peptide IDs for FAIMS CV " << group_cv << " V. Skipping." << endl;
        continue;
      }

      OPENMS_LOG_INFO << "  " << peptides_cv.size() << " internal IDs, "
                      << peptides_ext_cv.size() << " external IDs" << endl;

      // Create algorithm instance for this group (use same parameters)
      FeatureFinderIdentificationAlgorithm ffid_group;
      ffid_group.getProgressLogger().setLogType(prog_log_.getLogType());
      ffid_group.setParameters(this->getParameters());
      ffid_group.setMSData(std::move(faims_group));

      // Run feature detection
      FeatureMap features_cv;
      ffid_group.runSingleGroup_(peptides_cv, proteins, peptides_ext_cv, proteins_ext, features_cv, seeds, spectra_file);

      // Annotate features with FAIMS CV and add to results
      for (auto& feat : features_cv)
      {
        feat.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv);
        features.push_back(feat);
      }

      // Copy UnassignedPeptideIdentifications with FAIMS CV annotation
      for (auto& pep_id : features_cv.getUnassignedPeptideIdentifications())
      {
        pep_id.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv);
        features.getUnassignedPeptideIdentifications().push_back(std::move(pep_id));
      }

      // Copy ProteinIdentifications only from the first group
      if (first_group)
      {
        features.setProteinIdentifications(features_cv.getProteinIdentifications());
      }

      // Combine chromatograms
      for (const auto& chrom : ffid_group.getChromatograms().getChromatograms())
      {
        MSChromatogram chrom_copy = chrom;
        chrom_copy.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv);
        chrom_data_.addChromatogram(std::move(chrom_copy));
      }

      first_group = false;
    }

    // Warn about library output for FAIMS data
    OPENMS_LOG_WARN << "Warning: Library output is not available for multi-FAIMS data. "
                    << "Each FAIMS CV group has its own assay library." << endl;
    OPENMS_LOG_INFO << "Combined " << features.size() << " features from all FAIMS CV groups." << endl;

    // Set primary MS run path
    features.setPrimaryMSRunPath({spectra_file});
    features.ensureUniqueId();
  }

  void FeatureFinderIdentificationAlgorithm::runSingleGroup_(
    PeptideIdentificationList peptides,
    const vector<ProteinIdentification>& proteins,
    PeptideIdentificationList peptides_ext,
    vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds,
    const String& spectra_file
    )
  {
    // Clear output library from any previous run
    output_library_.clear(true);

    // Validate parameters
    validateSVMParameters_();

    // annotate mzML file
    features.setPrimaryMSRunPath({spectra_file}, ms_data_);

    // Check IM format
    double IM_window = param_.getValue("extract:IM_window");
    IMFormat im_format = IMTypes::determineIMFormat(ms_data_);
    bool has_IM = false;
    if (im_format == IMFormat::CONCATENATED)
    {
      has_IM = true;
    }
    else if (im_format != IMFormat::NONE) // has IM but wrong format
    {
      OPENMS_LOG_ERROR << "Wrong IM format detected. Expecting in concatenated format (float data arrays)" << std::endl;
    }

    // Initialize feature finder with appropriate parameters
    initializeFeatureFinder_();

    bool with_external_ids = !peptides_ext.empty();

    if (with_external_ids && !seeds.empty())
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Using seeds and external ids is currently not supported.");
    }

    double rt_uncertainty(0);
    if (with_external_ids)
    {
      // Use the external ID handler to align internal and external IDs
      rt_uncertainty = external_id_handler_.alignInternalAndExternalIDs(peptides, peptides_ext, rt_quantile_);
    }

    // Calculate RT window if not already set
    rt_window_ = calculateRTWindow_(rt_uncertainty);

    //-------------------------------------------------------------
    // prepare peptide map
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Preparing mapping of peptide data..." << endl;
    peptide_map_.clear();

    // Reserve enough space for all possible seeds
    {
      Size max_size = peptides.size() + seeds.size();
      if (add_mass_offset_peptides_)
      {
        max_size *= 2;
      }
      peptides.reserve(max_size);
    }

    for (auto& pep : peptides)
    {
      addPeptideToMap_(pep, peptide_map_); // stores pointer to pep in map
      pep.setMetaValue("FFId_category", "internal");
    }

    // Calculate global IM statistics BEFORE adding seeds
    // This ensures we only learn from real peptide identifications with IM data
    // and don't need to check/skip seeds during calculation
    calculateGlobalIMStats_();

    // TODO make sure that only assembled traces (more than one trace -> has a charge) if FFMetabo is used
    // see FeatureFindingMetabo: defaults_.setValue("remove_single_traces", "false", "Remove unassembled traces (single traces).");
    Size seeds_added = addSeeds_(peptides, seeds);
    OPENMS_LOG_INFO << "#Seeds without RT and m/z overlap with identified peptides added: " << seeds_added << endl;

    if (add_mass_offset_peptides_ > 0.0)
    {
      Size n_added = addOffsetPeptides_(peptides, add_mass_offset_peptides_);
      OPENMS_LOG_INFO << "#Offset peptides without RT and m/z overlap with other peptides added: " << n_added << endl;
    }

    n_internal_peps_ = peptide_map_.size();

    if (with_external_ids)
    {
      // Process and add external peptides
      for (PeptideIdentification& pep : peptides_ext)
      {
        addPeptideToMap_(pep, peptide_map_, true);
        pep.setMetaValue("FFId_category", "external");
      }
      n_external_peps_ = peptide_map_.size() - n_internal_peps_;
    }

    boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
    OpenSwath::SpectrumAccessPtr spec_temp =
        SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
    auto chunks = chunk_(peptide_map_.begin(), peptide_map_.end(), batch_size_);

    PeptideRefRTMap ref_rt_map;

    if (debug_level_ >= 668)
    {
      OPENMS_LOG_INFO << "Creating full assay library for debugging." << endl;
      // Warning: this step is pretty inefficient, since it does the whole library generation twice
      // Really use for debug only
      createAssayLibrary_(peptide_map_.begin(), peptide_map_.end(), ref_rt_map, false);
      cout << "Writing debug.traml file." << endl;
      FileHandler().storeTransitions("debug.traml", library_);
      ref_rt_map.clear();
      library_.clear(true);
    }

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    //Note: progress only works in non-debug when no logs come in-between
    getProgressLogger().startProgress(0, chunks.size(), "Creating assay library and extracting chromatograms");
    Size chunk_count = 0;
    for (auto& chunk : chunks)
    {
      //TODO since ref_rt_map is only used after chunking, we could create
      // maps per chunk and merge them in the end. Would help in parallelizing as well.
      // fills library_ (TargetedExperiment)
      createAssayLibrary_(chunk.first, chunk.second, ref_rt_map);
      OPENMS_LOG_DEBUG << "#Transitions: " << library_.getTransitions().size() << endl;

      ChromatogramExtractor extractor;
      // extractor.setLogType(ProgressLogger::NONE);
      {
        vector<OpenSwath::ChromatogramPtr> chrom_temp;
        vector<ChromatogramExtractor::ExtractionCoordinates> coords;
        // take entries in library_ and put to chrom_temp and coords
        extractor.prepare_coordinates(chrom_temp, coords, library_,
                                      numeric_limits<double>::quiet_NaN(), false);

        if (has_IM && IM_window > 0.0)
        {
          extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
                                        mz_window_ppm_, IM_window, "tophat");
        }
        else
        {
          extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
                                        mz_window_ppm_, "tophat");
        }

        extractor.return_chromatogram(chrom_temp, coords, library_, (*shared)[0],
                                      chrom_data_.getChromatograms(), false);
      }

      OPENMS_LOG_DEBUG << "Extracted " << chrom_data_.getNrChromatograms()
                       << " chromatogram(s)." << endl;

      OPENMS_LOG_DEBUG << "Detecting chromatographic peaks..." << endl;
      // suppress status output from OpenSWATH, unless in debug mode:
      if (debug_level_ < 1)
      {
        OpenMS_Log_info.remove(cout);
      }
      feat_finder_.pickExperiment(chrom_data_, features, library_,
                                  TransformationDescription(), ms_data_);
      if (debug_level_ < 1)
      {
        OpenMS_Log_info.insert(cout); // revert logging change
      }
      chrom_data_.clear(true);
      // Accumulate library entries for output before clearing
      for (const auto& pep : library_.getPeptides())
      {
        output_library_.addPeptide(pep);
      }
      for (const auto& prot : library_.getProteins())
      {
        output_library_.addProtein(prot);
      }
      for (const auto& trans : library_.getTransitions())
      {
        output_library_.addTransition(trans);
      }
      library_.clear(true);
      // since chrom_data_ here is just a container for the chromatograms and identifications will be empty,
      // pickExperiment above will only add empty ProteinIdentification runs with colliding identifiers.
      // Usually we could sanitize the identifiers or merge the runs, but since they are empty and we add the
      // "real" proteins later -> just clear them
      features.getProteinIdentifications().clear();
      getProgressLogger().setProgress(++chunk_count);
    }
    getProgressLogger().endProgress();

    OPENMS_LOG_INFO << "Found " << features.size() << " feature candidates in total."
                    << endl;

    ms_data_.reset(); // not needed anymore, free up the memory
    // complete feature annotation:
    annotateFeatures_(features, ref_rt_map);

    // sort everything:
    sort(features.getUnassignedPeptideIdentifications().begin(),
         features.getUnassignedPeptideIdentifications().end(),
         peptide_compare_);
    sort(features.begin(), features.end(), feature_compare_);

    postProcess_(features, with_external_ids);
    statistics_(features);

    features.setProteinIdentifications(proteins);
    // add external IDs (if any):
    features.getProteinIdentifications().insert(
      features.getProteinIdentifications().end(), proteins_ext.begin(),
      proteins_ext.end());
    features.getUnassignedPeptideIdentifications().insert(
      features.getUnassignedPeptideIdentifications().end(),
      peptides_ext.begin(), peptides_ext.end());

    // remove all hits with pseudo ids (seeds)
    removeSeedPseudoIDs_(features);

    // add back ignored PSMs
    features.getUnassignedPeptideIdentifications().insert(features.getUnassignedPeptideIdentifications().end(),
                                                          std::move_iterator(unassignedIDs_.begin()),
                                                          std::move_iterator(unassignedIDs_.end()));

    features.ensureUniqueId();
  }

  void FeatureFinderIdentificationAlgorithm::postProcess_(
   FeatureMap & features,
   bool with_external_ids)
  {   
    // don't do SVM stuff unless we have external data to apply the model to:
    if (with_external_ids)
    {
      external_id_handler_.classifyFeaturesWithSVM(features, param_);
    }
    // make sure proper unique ids get assigned to all features
    features.ensureUniqueId();

    // store feature candidates before filtering
    if (!candidates_out_.empty())
    {
      FileHandler().storeFeatures(candidates_out_, features);
    }

    // Use ExternalIDHandler for feature filtering
    if (with_external_ids)
    {
      external_id_handler_.filterClassifiedFeatures(features, external_id_handler_.getSVMProbsInternal().empty() ? 0.0 : double(param_.getValue("svm:min_prob")));
      OPENMS_LOG_INFO << features.size() << " features left after filtering." << endl;
    }
    else
    {
      filterFeatures_(features, with_external_ids);
      OPENMS_LOG_INFO << features.size() << " features left after filtering." << endl;
    }

    if (features.empty()) return; // elution model fit throws on empty features

    // Calculate FDR if we have external IDs
    if (with_external_ids) 
    {
      external_id_handler_.calculateFDR(features);
    }     
    
    //TODO MRMFeatureFinderScoring already does an ElutionModel scoring. It uses EMG fitting.
    // Would be nice if we could only do the fitting once, since it is one of the bottlenecks.
    // What is the intention of this post-processing here anyway? Does it filter anything?
    // If so, why not filter based on the corresponding Swath/MRM score?
    if (elution_model_ != "none")
    {
      ElutionModelFitter emf;
      Param emf_params = param_.copy("model:", true);
      emf_params.remove("type");
      emf_params.setValue("asymmetric",
                          (elution_model_ == "asymmetric") ? "true" : "false");
      emf.setParameters(emf_params);
      emf.fitElutionModels(features);
    }
    else if (!candidates_out_.empty()) // hulls not needed, remove them
    {
      for (auto & feat : features)
      {
        for (auto & sub : feat.getSubordinates())
        {
          sub.getConvexHulls().clear();
        }
      }
    }

  }

  void FeatureFinderIdentificationAlgorithm::runOnCandidates(FeatureMap & features)
  {
    if ((svm_n_samples_ > 0) && (svm_n_samples_ < 2 * svm_n_parts_))
    {
      String msg = "Sample size of " + String(svm_n_samples_) +
        " (parameter 'svm:samples') is not enough for " + String(svm_n_parts_) +
        "-fold cross-validation (parameter 'svm:xval').";
      throw Exception::InvalidParameter(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, msg);
    }

    bool with_external_ids = (!features.empty() && features[0].metaValueExists("predicted_class"));

    // extract ID information for statistics:
    peptide_map_.clear();
    set<AASequence> internal_seqs;
    for (PeptideIdentification& pep : features.getUnassignedPeptideIdentifications())
    {
      const AASequence& seq = pep.getHits()[0].getSequence();
      if (pep.getMetaValue("FFId_category") == "internal")
      {
        internal_seqs.insert(seq);
      }
      peptide_map_[seq];
    }
    for (const Feature& feat : features)
    {
      if (feat.getPeptideIdentifications().empty())
      {
        continue;
      }
      const PeptideIdentification& pep_id = feat.getPeptideIdentifications()[0];
      const AASequence& seq = pep_id.getHits()[0].getSequence();
      if (pep_id.getMetaValue("FFId_category") == "internal")
      {
        internal_seqs.insert(seq);
      }
      peptide_map_[seq];
    }
    n_internal_peps_ = internal_seqs.size();
    n_external_peps_ = peptide_map_.size() - internal_seqs.size();

    // sort everything:
    sort(features.getUnassignedPeptideIdentifications().begin(),
         features.getUnassignedPeptideIdentifications().end(),
         peptide_compare_);
    sort(features.begin(), features.end(), feature_compare_);

    postProcess_(features, with_external_ids);

    statistics_(features);

  }

  void FeatureFinderIdentificationAlgorithm::statistics_(FeatureMap const & features) const
  {
    // same peptide sequence may be quantified based on internal and external
    // IDs if charge states differ!
    set<AASequence> quantified_internal, quantified_all;
    for (const auto& f : features)
    {
      const PeptideIdentification& pep_id = f.getPeptideIdentifications()[0];
      const AASequence& seq = pep_id.getHits()[0].getSequence();
      if (f.getIntensity() > 0.0)
      {
        quantified_all.insert(seq);
        if (pep_id.getMetaValue("FFId_category") == "internal")
        {
          quantified_internal.insert(seq);
        }
      }
    }
    Size n_quant_external = quantified_all.size() - quantified_internal.size();
    // If internal and external IDs for a peptide map to different RT regions,
    // it is possible that there is a quantification from the "external" region,
    // but not from the "internal" region (no matching feature) - therefore the
    // number of "missing" external peptides can be negative!
    Int n_missing_external = Int(n_external_peps_) - n_quant_external;

    OPENMS_LOG_INFO << "\nSummary statistics (counting distinct peptides including "
      "PTMs):\n"
             << peptide_map_.size() << " peptides identified ("
             << n_internal_peps_ << " internal, " << n_external_peps_
             << " additional external)\n"
             << quantified_all.size() << " peptides with features ("
             << quantified_internal.size() << " internal, "
             << n_quant_external << " external)\n"
             << peptide_map_.size() - quantified_all.size()
             << " peptides without features ("
             << n_internal_peps_ - quantified_internal.size() << " internal, "
             << n_missing_external << " external)\n" << endl;

  }

  void FeatureFinderIdentificationAlgorithm::calculateGlobalIMStats_()
  {
    // Update ranges from MS data to get IM range from raw data
    ms_data_.updateRanges();

    // Try to get IM range from MS data (will throw if no IM data available)
    try
    {
      global_im_stats_.min = ms_data_.getMinMobility();
      global_im_stats_.max = ms_data_.getMaxMobility();
    }
    catch (const Exception::InvalidRange&)
    {
      OPENMS_LOG_DEBUG << "No IM data found in MS data." << endl;
      global_im_stats_ = IMStats(); // Empty stats
      return;
    }

    // Calculate median from peptide identifications for robust central tendency
    // (more representative of where peptides actually elute than data range center)
    std::vector<double> im_values_from_ids;

    // Collect IM values from all peptide identifications
    // Note: This is called BEFORE adding seeds, so we only see real identifications
    for (const auto& pep_entry : peptide_map_)
    {
      const ChargeMap& charge_map = pep_entry.second;
      for (const auto& charge_entry : charge_map)
      {
        const RTMap& internal_ids = charge_entry.second.first;
        for (const auto& rt_pepid : internal_ids)
        {
          const PeptideIdentification& pep_id = *rt_pepid.second;
          const double im = pep_id.getMetaValue(Constants::UserParam::IM, -1.0);
          if (im >= 0.0)  // Only collect valid IM values (>= 0.0 matches ChromatogramExtractor)
          {
            im_values_from_ids.push_back(im);
          }
        }
      }
    }

    // If we have ID-based IM values, use them for median calculation
    // Otherwise, use center of data range as fallback
    if (!im_values_from_ids.empty())
    {
      std::sort(im_values_from_ids.begin(), im_values_from_ids.end());
      Size n = im_values_from_ids.size();
      if (n % 2 == 0)
      {
        global_im_stats_.median = (im_values_from_ids[n/2 - 1] + im_values_from_ids[n/2]) / 2.0;
      }
      else
      {
        global_im_stats_.median = im_values_from_ids[n/2];
      }

      OPENMS_LOG_INFO << "Global IM statistics: median=" << global_im_stats_.median
                      << " (from " << n << " identifications), "
                      << "min=" << global_im_stats_.min << ", "
                      << "max=" << global_im_stats_.max << " (from MS data range)" << endl;
    }
    else
    {
      // No IDs with IM annotation - use center of data range
      global_im_stats_.median = (global_im_stats_.min + global_im_stats_.max) / 2.0;

      OPENMS_LOG_INFO << "Global IM statistics: median=" << global_im_stats_.median
                      << " (center of data range), "
                      << "min=" << global_im_stats_.min << ", "
                      << "max=" << global_im_stats_.max << " (from MS data range)" << endl;
    }
  }

  /**
   * @brief Calculate ion mobility statistics (median, min, max) for a given RT region
   *
   * This function computes robust statistics for ion mobility values from peptide
   * identifications within a single RT region. The statistics are used for:
   * 1. Setting the drift time on peptide assays (using median)
   * 2. Extracting chromatograms with appropriate IM windows
   * 3. Annotating features with IM QC metrics
   *
   * Implementation Strategy:
   * - Collects IM values from internal peptide identifications in the RT region
   * - Skips individual IDs that lack IM annotation (logged as debug, summary as info)
   * - Uses MEDIAN instead of mean for robustness against outliers
   * - Computes min/max to characterize the IM distribution spread
   * - Returns empty stats only if NO valid IM values are available
   *
   * Note: Seeds from untargeted feature finders may or may not have IM annotation,
   * depending on the feature finder. If IM is annotated, it is used; otherwise the
   * seed is extracted across the full IM range (ChromatogramExtractor disables IM
   * filtering when ion_mobility < 0).
   *
   * Note: RT region boundaries are determined from ALL IDs (including those without IM),
   * so skipping individual IDs for IM statistics does not affect RT extraction.
   *
   * @param r RT region containing peptide identifications (per charge state)
   * @return IMStats structure with median, min, and max IM values
   *         Returns {-1.0, -1.0, -1.0} only if no valid IM data is available
   *
   * @note The median is calculated using the standard definition:
   *       - For odd n: middle element of sorted array
   *       - For even n: average of two middle elements
   */
  FeatureFinderIdentificationAlgorithm::IMStats
  FeatureFinderIdentificationAlgorithm::getRTRegionIMStats_(const RTRegion& r)
  {
    IMStats stats;
    const ChargeMap& cm = r.ids;
    std::vector<double> im_values;
    Size n_ids_without_im = 0;

    // Collect all IM values from internal IDs
    for (const auto& e : cm)
    {
      const RTMap& internal_ids = e.second.first; // internal
      for (const auto& rt_pepidptr : internal_ids)
      {
        const PeptideIdentification& pep_id = *rt_pepidptr.second;
        const double im = pep_id.getMetaValue(Constants::UserParam::IM, -1.0);

        if (im < 0.0)
        {
          // ID without IM (negative value) - skip this ID but continue with others
          n_ids_without_im++;
          OPENMS_LOG_DEBUG << "Identification at RT " << pep_id.getRT()
                           << " lacks IM annotation - skipping for IM statistics" << endl;
        }
        else
        {
          im_values.push_back(im); // includes 0.0 as valid IM value
        }
      }
    }

    if (im_values.empty())
    {
      return stats; // Return empty stats
    }

    // Calculate statistics
    std::sort(im_values.begin(), im_values.end());

    stats.min = im_values.front();
    stats.max = im_values.back();

    // Calculate median
    Size n = im_values.size();
    if (n % 2 == 0)
    {
      stats.median = (im_values[n/2 - 1] + im_values[n/2]) / 2.0;
    }
    else
    {
      stats.median = im_values[n/2];
    }

    if (n_ids_without_im > 0)
    {
      OPENMS_LOG_INFO << "Calculated IM statistics from " << n << " IDs with IM data "
                      << "(skipped " << n_ids_without_im << " IDs without IM annotation)" << endl;
    }

    return stats;
  }

  void FeatureFinderIdentificationAlgorithm::createAssayLibrary_(
    const PeptideMap::iterator& begin, 
    const PeptideMap::iterator& end, 
    PeptideRefRTMap& ref_rt_map, 
    bool clear_IDs)
  {
    std::set<String> protein_accessions;

    Size seedcount = 0;
    for (auto pm_it = begin;
         pm_it != end; ++pm_it)
    {
      TargetedExperiment::Peptide peptide;
      const AASequence &seq = pm_it->first;


      // @NOTE: Technically, "TargetedExperiment::Peptide" stores the unmodified
      // sequence and the modifications separately. Unfortunately, creating the
      // modifications vector is complex and there is currently no convenient
      // conversion function (see "TargetedExperimentHelper::getAASequence" for
      // the reverse conversion). However, "Peptide" is later converted to
      // "OpenSwath::LightPeptide" anyway, and this is done via "AASequence"
      // (see "OpenSwathDataAccessHelper::convertTargetedPeptide"). So for our
      // purposes it works to just store the sequence including modifications in
      // "Peptide".

      // for now, seeds are stored in the same PeptideRefMap, all
      // under the same fake sequence key entry
      // TODO add own data structure for them
      if (seq.toUnmodifiedString().hasPrefix("XXX")) // seed
      {
        // This will force the SWATH scores to consider it like an unidentified peptide and e.g. use averagine isotopes
        peptide.sequence = "";
        // we do not have to aggregate their retention times, therefore just
        // iterate over the entries
        const ChargeMap& cm = pm_it->second;
        for (const auto& charge_rtmap : cm)
        {
          Int charge = charge_rtmap.first;
          // only go through internals for seeds (->first). External seeds are not supported
          for (const auto& rt_pep : charge_rtmap.second.first)
          {
            // since we don't know their IDs, seeds will all need a different grouplabel in SWATH
            // to not be combined
            seedcount++;

            double mz = rt_pep.second->getMZ();
            double rt = rt_pep.second->getRT();
            String uid = rt_pep.second->getMetaValue("SeedFeatureID");

            // UID should be enough, but let's add the seed count to be sure.
            String peptide_id = peptide.sequence + "[" + uid + "][" + String(seedcount) + "]/" + String(charge);
            peptide.setChargeState(charge);
            peptide.id = peptide_id;
            peptide.protein_refs = {"not_available"};
            peptide.setPeptideGroupLabel(peptide_id);

            //create an entry in the "output" ref_rt_map for internals
            RTMap &internal_ids = ref_rt_map[peptide_id].first;

            // get isotope distribution for peptide:
            //TODO Why 10? Document constant?
            Size n_isotopes = (isotope_pmin_ > 0.0) ? 10 : n_isotopes_;
            CoarseIsotopePatternGenerator generator(n_isotopes);
            IsotopeDistribution iso_dist = generator
                .estimateFromPeptideWeight(mz * charge - charge * Constants::PROTON_MASS_U);
            if (isotope_pmin_ > 0.0)
            {
              iso_dist.trimLeft(isotope_pmin_);
              iso_dist.trimRight(isotope_pmin_);
              iso_dist.renormalize();
            }

            double rt_tolerance = seed_rt_window_ / 2.0;

            // store beginning and end of RT region: here we only need one entry
            peptide.rts.clear();
            addPeptideRT_(peptide, rt - rt_tolerance);
            addPeptideRT_(peptide, rt + rt_tolerance);

            // Use IM from seed if annotated (some feature finders provide IM)
            // If not annotated, drift time stays at default (-1) -> full IM range extraction
            // Check >= 0.0 to match ChromatogramExtractor's IM filtering logic
            const double seed_im = rt_pep.second->getMetaValue(Constants::UserParam::IM, -1.0);
            if (seed_im >= 0.0)
            {
              peptide.setDriftTime(seed_im);
              // Store IM stats for annotation (use seed IM as median, with some tolerance for min/max)
              im_stats_[peptide.id] = {seed_im, seed_im, seed_im};
            }
            else
            {
              // Reset drift time to -1 (no IM filtering) - peptide object is reused
              peptide.setDriftTime(-1.0);
            }

            library_.addPeptide(peptide);
            generateTransitions_(peptide.id, mz, charge, iso_dist);
            internal_ids.emplace(rt_pep);
          }
        }
      }
      else
      {
        peptide.sequence = seq.toString();
        // keep track of protein accessions:
        set<String> current_accessions;
        
        const pair<RTMap, RTMap> &pair = pm_it->second.begin()->second; // internal/external pair
        const RTMap& internal_ids = pair.first;
        const RTMap& external_ids = pair.second;

        // WARNING: This assumes that at least one hit is present.
        const PeptideHit &hit = (internal_ids.empty() ?
                                 external_ids.begin()->second->getHits()[0] :
                                 internal_ids.begin()->second->getHits()[0]);
        current_accessions = hit.extractProteinAccessionsSet();
        protein_accessions.insert(current_accessions.begin(),
                                  current_accessions.end());
        // missing protein accession would crash OpenSWATH algorithms:
        if (current_accessions.empty())
        {
          current_accessions.insert("not_available");
        }

        peptide.protein_refs = vector<String>(current_accessions.begin(),
                                              current_accessions.end());
        // get regions in RT which peptide eludes (ideally only one):
        std::vector<RTRegion> rt_regions;
        getRTRegions_(pm_it->second, rt_regions, clear_IDs);

        // note: IM values are stored in the PeptideIdentifications* for the different
        // peptides, charges, and regions

        // get isotope distribution for peptide:
        Size n_isotopes = (isotope_pmin_ > 0.0) ? 10 : n_isotopes_;
        IsotopeDistribution iso_dist =
            seq.getFormula(Residue::Full, 0).getIsotopeDistribution(CoarseIsotopePatternGenerator(n_isotopes));
        if (isotope_pmin_ > 0.0)
        {
          iso_dist.trimLeft(isotope_pmin_);
          iso_dist.trimRight(isotope_pmin_);
          iso_dist.renormalize();
        }

        // go through different charge states:
        for (ChargeMap::const_iterator cm_it = pm_it->second.begin();
             cm_it != pm_it->second.end(); ++cm_it)
        {
          Int charge = cm_it->first;

          double mz = seq.getMZ(charge);
          OPENMS_LOG_DEBUG << "\nPeptide " << peptide.sequence << "/" << charge << " (m/z: " << mz << "):" << endl;
          peptide.setChargeState(charge);
          String peptide_id = peptide.sequence + "/" + String(charge);

          // we want to detect one feature per peptide and charge state - if there
          // are multiple RT regions, group them together:
          peptide.setPeptideGroupLabel(peptide_id);
          peptide.rts.clear();
          Size counter = 0;
          // accumulate IDs over multiple regions:
          RTMap &internal_ids = ref_rt_map[peptide_id].first;
          RTMap &external_ids = ref_rt_map[peptide_id].second;
          for (RTRegion& reg : rt_regions)
          {
            if (reg.ids.count(charge))
            {
              OPENMS_LOG_DEBUG_NOFILE << "Charge " << charge << ", Region# " << counter + 1 << " (RT: "
                               << float(reg.start) << "-" << float(reg.end)
                               << ", size " << float(reg.end - reg.start) << ")"
                               << endl;

              peptide.id = peptide_id;
              if (rt_regions.size() > 1)
                peptide.id += ":" + String(++counter);

              // store beginning and end of RT region:
              peptide.rts.clear();
              addPeptideRT_(peptide, reg.start);
              addPeptideRT_(peptide, reg.end);

              // determine IM statistics (median, min, max)
              // for the peptide and current charge state in the region
              // (Note: because it is the same peptide and charge state the IM should not differ that much)
              // Check >= 0.0 to match ChromatogramExtractor's IM filtering logic
              IMStats im_stats = getRTRegionIMStats_(reg);
              if (im_stats.median >= 0.0)
              {
                peptide.setDriftTime(im_stats.median); // use median (more robust than mean)
                im_stats_[peptide.id] = im_stats; // store for later annotation
              }
              else
              {
                // Reset drift time to -1 (no IM filtering) - peptide object is reused
                peptide.setDriftTime(-1.0);
              }

              library_.addPeptide(peptide);
              generateTransitions_(peptide.id, mz, charge, iso_dist);
            }
            internal_ids.insert(reg.ids[charge].first.begin(),
                                reg.ids[charge].first.end());
            external_ids.insert(reg.ids[charge].second.begin(),
                                reg.ids[charge].second.end());
          }
        }
      }
    }
    // add proteins to library:
    for (String const &acc : protein_accessions)
    {
      TargetedExperiment::Protein protein;
      protein.id = acc;
      library_.addProtein(protein);
    }
  }

  // extract RT regions of identified peptides (from charge map)
  void FeatureFinderIdentificationAlgorithm::getRTRegions_(
    ChargeMap& peptide_data,
    std::vector<RTRegion>& rt_regions,
    bool clear_IDs) const
  {
    // use RTs from all charge states of a single peptide to get a more complete picture:
    std::vector<double> rts;
    for (auto& cm : peptide_data)
    {
      // "internal" IDs:
      for (auto& rt : cm.second.first)
      {
        rts.push_back(rt.first);
      }
      // "external" IDs:
      for (auto& rt : cm.second.second)
      {
        rts.push_back(rt.first);
      }
    }
    sort(rts.begin(), rts.end());
    double rt_tolerance = rt_window_ / 2.0;

    for (auto& rt : rts)
    {
      // large gap between last RT of last region and current RT? then create a new region?
      if (rt_regions.empty() || (rt_regions.back().end < rt - rt_tolerance))
      {
        RTRegion region;
        region.start = rt - rt_tolerance;
        // TODO
        // cppcheck-suppress uninitStructMember
        rt_regions.push_back(region);
      }
      rt_regions.back().end = rt + rt_tolerance;
    }

    // sort the peptide IDs into the regions:
    for (auto& cm : peptide_data)
    {
      // regions are sorted by RT, as are IDs, so just iterate linearly:
      std::vector<RTRegion>::iterator reg_it = rt_regions.begin();
      int charge = cm.first;

      // "internal" IDs:
      for (auto& rt : cm.second.first)
      {
        // while RT larger than current region end: skip to next region (or end)
        while (rt.first > reg_it->end) 
        {
          ++reg_it;
        }
        RTMap& internal_ids = reg_it->ids[charge].first;
         // insert RT and peptide id object into multimap (for current charge of the peptide)
        internal_ids.insert(rt);
      }
      reg_it = rt_regions.begin(); // reset to start
      // "external" IDs:
      for (auto& rt : cm.second.second)
      {
        while (rt.first > reg_it->end)
        {
          ++reg_it;
        }
        RTMap& external_ids = reg_it->ids[charge].second;
        external_ids.insert(rt);
      }
      if (clear_IDs)
      {
        // ID references no longer needed (now stored in the RT regions):
        cm.second.first.clear();
        cm.second.second.clear();
      }
    }
  }

  void FeatureFinderIdentificationAlgorithm::addPeptideRT_(TargetedExperiment::Peptide& peptide, double rt) const
  {
    TargetedExperiment::RetentionTime te_rt;
    te_rt.setRT(rt);
    te_rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED;
    peptide.rts.push_back(te_rt);
  }

  /// generate transitions (isotopic traces) for a peptide ion and add them to the library:
  void FeatureFinderIdentificationAlgorithm::generateTransitions_(
    const String& peptide_id, 
    double mz, 
    Int charge,
    const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    for (const Peak1D& iso : iso_dist)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter + 1);
      String transition_name = peptide_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U * float(counter) / charge);
      transition.setLibraryIntensity(iso.getIntensity());
      transition.setMetaValue("annotation", annotation);
      transition.setPeptideRef(peptide_id);

      //TODO what about transition charge? A lot of DIA scores depend on it and default to charge 1 otherwise.
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso.getIntensity();
      ++counter;
    }
  }

  void FeatureFinderIdentificationAlgorithm::annotateFeaturesFinalizeAssay_(
    FeatureMap& features, map<Size, vector<PeptideIdentification*> >& feat_ids,
    RTMap& rt_internal)
  {
    set<PeptideIdentification*> assigned_ids;
    if (!feat_ids.empty())
    {
      // find the "best" feature (with the most IDs):
      Size best_index = 0;
      Size best_count = 0;
      for (map<Size, vector<PeptideIdentification*> >::iterator fi_it =
             feat_ids.begin(); fi_it != feat_ids.end(); ++fi_it)
      {
        Size current_index = fi_it->first;
        Size current_count = fi_it->second.size();
        if ((current_count > best_count) ||
            ((current_count == best_count) && // break ties by intensity
             (features[current_index].getIntensity() >
              features[best_index].getIntensity())))
        {
          best_count = current_count;
          best_index = current_index;
        }
      }
      // assign IDs:
      if (best_count > 0)
      {
        // we define the (one) feature with most matching IDs as correct:
        features[best_index].setMetaValue("feature_class", "positive");
        features[best_index].getPeptideIdentifications().resize(best_count);
        for (Size i = 0; i < best_count; ++i)
        {
          features[best_index].getPeptideIdentifications()[i] =
              *(feat_ids[best_index][i]);
        }
        assigned_ids.insert(feat_ids[best_index].begin(),
                            feat_ids[best_index].end());
      }
    }
    // store unassigned IDs from the current RT region:
    for (RTMap::const_iterator rt_it = rt_internal.begin();
         rt_it != rt_internal.end(); ++rt_it)
    {
      if (!assigned_ids.count(rt_it->second))
      {
        const PeptideIdentification& pep_id = *(rt_it->second);
        features.getUnassignedPeptideIdentifications().push_back(pep_id);
      }
    }
    // clean-up:
    feat_ids.clear();
    rt_internal.clear();
  }

  /// annotate identified features with m/z, isotope probabilities, etc.
  void FeatureFinderIdentificationAlgorithm::annotateFeatures_(FeatureMap& features, PeptideRefRTMap& ref_rt_map)
  {
    String previous_ref, peptide_ref;
    RTMap transformed_internal;
    Size i = 0;
    map<Size, vector<PeptideIdentification*> > feat_ids;
    for (Feature& feat : features)
    {
      feat.setMZ(feat.getMetaValue("PrecursorMZ"));
      feat.setCharge(feat.getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      ensureConvexHulls_(feat);
      // remove "fake" IDs generated by OpenSWATH (they would be removed with
      // a warning when writing output, because of missing protein
      // identification with corresponding identifier):
      feat.getPeptideIdentifications().clear();
      // annotate subordinates with theoretical isotope intensities:
      for (Feature& sub : feat.getSubordinates())
      {
        String native_id = sub.getMetaValue("native_id");
        sub.setMetaValue("isotope_probability", isotope_probs_[native_id]);
      }

      peptide_ref = feat.getMetaValue("PeptideRef");

      // Annotate feature with ion mobility statistics (if available)
      // This provides QC metrics for IM-enabled experiments:
      // - IM_median: robust central IM value used for chromatogram extraction
      // - IM_min/max: spread of IM distribution (large spread may indicate issues)
      // Note: Uses full peptide ref (with region number) as this is the key in im_stats_
      String full_peptide_ref = peptide_ref; // keep full ref with region number
      if (im_stats_.count(full_peptide_ref))
      {
        const IMStats& stats = im_stats_.at(full_peptide_ref);
        feat.setMetaValue("IM_median", stats.median);
        feat.setMetaValue("IM_min", stats.min);
        feat.setMetaValue("IM_max", stats.max);
      }

      // remove region number, if present:
      Size pos_slash = peptide_ref.rfind('/');
      Size pos_colon = peptide_ref.find(':', pos_slash + 2);
      peptide_ref = peptide_ref.substr(0, pos_colon);

      if (peptide_ref != previous_ref)
      {
        if (!previous_ref.empty())
        {
          annotateFeaturesFinalizeAssay_(
            features, feat_ids, ref_rt_map[previous_ref].first);
        }
        previous_ref = peptide_ref;
      }

      RTMap& rt_internal = ref_rt_map[peptide_ref].first;
      RTMap& rt_external = ref_rt_map[peptide_ref].second;

      if (rt_internal.empty() && rt_external.empty())
      {
        OPENMS_LOG_DEBUG << "PeptideRefs in RTMap:" << endl;
        for (const auto& rtm : ref_rt_map)
        {
          OPENMS_LOG_DEBUG << rtm.first << endl;
        }

        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RT internal and external are both empty for peptide '" + String(peptide_ref) + "' stored as '" + String(feat.getMetaValue("PeptideRef")) + "'.");
      }

      if (!rt_internal.empty()) // validate based on internal IDs
      {
        // map IDs to features (based on RT):
        double rt_min = features[i].getMetaValue("leftWidth");
        double rt_max = features[i].getMetaValue("rightWidth");
        std::tie(rt_min, rt_max) = calculateRTBounds_(rt_min, rt_max);

        RTMap::const_iterator lower = rt_internal.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_internal.upper_bound(rt_max);
        int id_count = 0;
        for (; lower != upper; ++lower)
        {
          feat_ids[i].push_back(lower->second);
          ++id_count;
        }
        // "total" only includes IDs from this RT region:
        feat.setMetaValue("n_total_ids", rt_internal.size());
        feat.setMetaValue("n_matching_ids", id_count);
        if (id_count > 0) // matching IDs -> feature may be correct
        {
          feat.setMetaValue("feature_class", "ambiguous");
        }
        else // no matching IDs -> feature is wrong
        {
          feat.setMetaValue("feature_class", "negative");
        }
      }
      else // only external IDs -> no validation possible
      {
        // Set feature class to unknown
        feat.setMetaValue("n_total_ids", 0);
        feat.setMetaValue("n_matching_ids", -1);
        feat.setMetaValue("feature_class", "unknown");
        
        // Add a dummy peptide identification from external data
        if (!rt_external.empty())
        {
          PeptideIdentification id = *(rt_external.begin()->second);
          id.clearMetaInfo();
          id.setMetaValue("FFId_category", "implied");
          id.setRT(feat.getRT());
          id.setMZ(feat.getMZ());
          // only one peptide hit per ID - see function "addPeptideToMap_":
          PeptideHit& hit = id.getHits()[0];
          hit.clearMetaInfo();
          hit.setScore(0.0);
          feat.getPeptideIdentifications().push_back(id);
        }
      }

      // distance from feature to closest peptide ID:
      if (external_id_handler_.hasRTTransformation())
      {
        // use external IDs if available, otherwise RT-transformed internal IDs
        // (but only compute the transform if necessary, once per assay!):
        if (rt_external.empty() && (transformed_internal.empty() ||
                                     (peptide_ref != previous_ref)))
        {
          transformed_internal.clear();
          for (RTMap::const_iterator it = rt_internal.begin();
               it != rt_internal.end(); ++it)
          {
            double transformed_rt = external_id_handler_.transformRT(it->first);
            RTMap::value_type pair = make_pair(transformed_rt, it->second);
            transformed_internal.insert(transformed_internal.end(), pair);
          }
        }
        const RTMap& rt_ref = (rt_external.empty() ? transformed_internal :
                                rt_external);

        double rt_min = feat.getMetaValue("leftWidth");
        double rt_max = feat.getMetaValue("rightWidth");
        std::tie(rt_min, rt_max) = calculateRTBounds_(rt_min, rt_max);

        RTMap::const_iterator lower = rt_ref.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_ref.upper_bound(rt_max);
        if (lower != upper) // there's at least one ID within the feature
        {
          feat.setMetaValue("rt_delta", 0.0);
        }
        else // check closest ID
        {
          double rt_delta1 = numeric_limits<double>::infinity();
          if (lower != rt_ref.begin())
          {
            rt_delta1 = fabs((--lower)->first - rt_min);
          }
          double rt_delta2 = numeric_limits<double>::infinity();
          if (upper != rt_ref.end())
          {
            rt_delta2 = fabs(upper->first - rt_min);
          }
          feat.setMetaValue("rt_delta", min(rt_delta1, rt_delta2));
        }
      }
      ++i;
    }
    // set of features from the last assay:
    annotateFeaturesFinalizeAssay_(features, feat_ids,
                                   ref_rt_map[peptide_ref].first);
    // store unassigned peptide IDs from assays that did not generate any
    // feature candidates:
    for (PeptideRefRTMap::iterator ref_it = ref_rt_map.begin();
         ref_it != ref_rt_map.end(); ++ref_it)
    {
      RTMap& rt_internal = ref_it->second.first;
      if (!rt_internal.empty()) // not cleared by '...FinalizeAssay()'
      {
        for (RTMap::const_iterator rt_it = rt_internal.begin();
             rt_it != rt_internal.end(); ++rt_it)
        {
          const PeptideIdentification& pep_id = *(rt_it->second);
          features.getUnassignedPeptideIdentifications().push_back(pep_id);
        }
      }
    }
  }

  void FeatureFinderIdentificationAlgorithm::ensureConvexHulls_(Feature& feature) const
  {
    if (feature.getConvexHulls().empty()) // add hulls for mass traces
    {
      double rt_min = feature.getMetaValue("leftWidth");
      double rt_max = feature.getMetaValue("rightWidth");
      for (Feature& sub : feature.getSubordinates())
      {
        double abs_mz_tol = mz_window_ / 2.0;
        if (mz_window_ppm_)
        {
          abs_mz_tol = sub.getMZ() * abs_mz_tol * 1.0e-6;
        }
        ConvexHull2D hull;
        hull.addPoint(DPosition<2>(rt_min, sub.getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_min, sub.getMZ() + abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub.getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub.getMZ() + abs_mz_tol));
        feature.getConvexHulls().push_back(hull);
      }
    }
  }

  void FeatureFinderIdentificationAlgorithm::addPeptideToMap_(PeptideIdentification& peptide, PeptideMap& peptide_map, bool external)
  {
    if (peptide.getHits().empty())
    {
      return;
    }
    peptide.sort();
    PeptideHit& hit = peptide.getHits()[0];
    peptide.getHits().resize(1);

    // if we don't quantify decoys we don't add them to the peptide list
    if (!quantify_decoys_)
    {
      if (hit.isDecoy())
      {
        unassignedIDs_.push_back(peptide);
        return;
      }
    }
    if (use_psm_cutoff_)
    {
      if ( (peptide.isHigherScoreBetter() && hit.getScore() < psm_score_cutoff_) ||
           (!peptide.isHigherScoreBetter() && hit.getScore() > psm_score_cutoff_) )
      {
        unassignedIDs_.push_back(peptide);
        return;
      }
    }

    Int charge = hit.getCharge();
    // precursor information
    double rt = peptide.getRT();
    double mz = peptide.getMZ();

    // meta value to forcefully overwrite m/z value with external one
    // to quantify this kind of data, we need to introduce a modification that matches the mass difference
    // we just start at the N-term and look for the first unmodified AA
    if (hit.metaValueExists("CalcMass"))
    {
      double diff_mz = (double)hit.getMetaValue("CalcMass") - hit.getSequence().getMZ(charge);
      double diff_mass = diff_mz * charge;
      if (fabs(diff_mass) > 0.01)
      {
        OPENMS_LOG_DEBUG_NOFILE << "Peptide m/z value and m/z of CalcMass meta value differ (" << hit.getSequence().getMZ(charge) << " / " << hit.getMetaValue("CalcMass")
                                << "Assuming unspecified/unlocalized modification." << endl;
        AASequence seq = hit.getSequence(); // TODO: add ref version to PeptideHit
        for (auto r = seq.begin(); r != seq.end(); ++r)
        {
          if (r->isModified()) continue;
          int residue_index = r - seq.begin();
          seq.setModificationByDiffMonoMass(residue_index, diff_mass);
          break;
        }
        hit.setSequence(std::move(seq)); 
      }
    }
    
    if (external)
    {
      OPENMS_LOG_DEBUG_NOFILE << "Adding peptide (external) " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      peptide_map[hit.getSequence()][charge].second.emplace(rt, &peptide);
    }
    else
    {
      if (peptide.metaValueExists("SeedFeatureID"))
      {
        OPENMS_LOG_DEBUG_NOFILE << "Adding seed (internal) from FeatureID " << peptide.getMetaValue("SeedFeatureID") << ": " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      }
      else
      {
        OPENMS_LOG_DEBUG_NOFILE << "Adding peptide (internal) " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      }

      peptide_map[hit.getSequence()][charge].first.emplace(rt, &peptide); // place into multimap
    }
  }

  void FeatureFinderIdentificationAlgorithm::updateMembers_()
  {
    peak_width_ = param_.getValue("detect:peak_width");
    min_peak_width_ = param_.getValue("detect:min_peak_width");
    signal_to_noise_ = param_.getValue("detect:signal_to_noise");

    batch_size_ = param_.getValue("extract:batch_size");
    rt_quantile_ = param_.getValue("extract:rt_quantile");
    rt_window_ = param_.getValue("extract:rt_window");
    mz_window_ = param_.getValue("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;

    isotope_pmin_ = param_.getValue("extract:isotope_pmin");
    n_isotopes_ = param_.getValue("extract:n_isotopes");

    mapping_tolerance_ = param_.getValue("detect:mapping_tolerance");

    elution_model_ = param_.getValue("model:type").toString();
    // SVM related parameters
    svm_min_prob_ = param_.getValue("svm:min_prob");
    svm_predictor_names_ = ListUtils::create<String>(param_.getValue("svm:predictors").toString());
    svm_xval_out_ = param_.getValue("svm:xval_out").toString();
    svm_quality_cutoff = param_.getValue("svm:min_prob");
    svm_n_parts_ = param_.getValue("svm:xval");
    svm_n_samples_ = param_.getValue("svm:samples");

    // debug
    debug_level_ = param_.getValue("debug");
    candidates_out_ = param_.getValue("candidates_out").toString();

    // quantification of decoys
    quantify_decoys_ = param_.getValue("quantify_decoys").toBool();
    use_psm_cutoff_ = param_.getValue("min_psm_cutoff") != "none";
    if (use_psm_cutoff_)
    {
      psm_score_cutoff_ = double(param_.getValue("min_psm_cutoff"));
    }

    add_mass_offset_peptides_ = double(param_.getValue("add_mass_offset_peptides"));
  }

  
  void FeatureFinderIdentificationAlgorithm::filterFeatures_(OpenMS::FeatureMap& features, bool classified)
  {
    if (features.empty())
    {
      return;
    }
    
    // For non-classified features, we still use the original filtering
    if (!classified)
    {
      // remove features without ID (or pseudo ID from seeds)
      features.erase(std::remove_if(features.begin(), features.end(),
                               feature_filter_peptides_), features.end());
    }
    // Note: The classified case is now handled by ExternalIDHandler::filterClassifiedFeatures
    // in the postProcess_ method
  }

}
