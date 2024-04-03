// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/MATH/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

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
    return library_;
  }

  const TargetedExperiment& FeatureFinderIdentificationAlgorithm::getLibrary() const
  {
    return library_;
  }


  Size FeatureFinderIdentificationAlgorithm::addOffsetPeptides_(vector<PeptideIdentification>& peptides, double offset)
  {
    // WARNING: Superhack! Use unique ID to distinguish seeds from real IDs. Use a mod that will never occur to
    // make them truly unique and not be converted to an actual modification.
    const String pseudo_mod_name = String(10000);
    AASequence some_seq = AASequence::fromString("XXX[" + pseudo_mod_name + "]");

    vector<PeptideIdentification> offset_peptides;
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

  Size FeatureFinderIdentificationAlgorithm::addSeeds_(vector<PeptideIdentification>& peptides, const FeatureMap& seeds)
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
        addPeptideToMap_(peptides.back(), peptide_map_);
        ++seeds_added;
      }
    }
    
    return seeds_added;
  }

  void FeatureFinderIdentificationAlgorithm::run(
    vector<PeptideIdentification> peptides,
    const vector<ProteinIdentification>& proteins,
    vector<PeptideIdentification> peptides_ext,
    vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds,
    const String& spectra_file
    )
  {
    if ((svm_n_samples_ > 0) && (svm_n_samples_ < 2 * svm_n_parts_))
    {
      String msg = "Sample size of " + String(svm_n_samples_) +
        " (parameter 'svm:samples') is not enough for " + String(svm_n_parts_) +
        "-fold cross-validation (parameter 'svm:xval').";
      throw Exception::InvalidParameter(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, msg);
    }

    // annotate mzML file
    features.setPrimaryMSRunPath({spectra_file}, ms_data_);

    // initialize algorithm classes needed later:
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

    double rt_uncertainty(0);
    bool with_external_ids = !peptides_ext.empty();

    if (with_external_ids && !seeds.empty())
    {
      throw Exception::IllegalArgument(
        __FILE__,
        __LINE__,
        OPENMS_PRETTY_FUNCTION,
        "Using seeds and external ids is currently not supported.");
    }

    if (with_external_ids)
    {
      // align internal and external IDs to estimate RT shifts:
      MapAlignmentAlgorithmIdentification aligner;
      aligner.setReference(peptides_ext); // go from internal to external scale
      vector<vector<PeptideIdentification> > aligner_peptides(1, peptides);
      vector<TransformationDescription> aligner_trafos;

      OPENMS_LOG_INFO << "Realigning internal and external IDs...";
      aligner.align(aligner_peptides, aligner_trafos);
      trafo_external_ = aligner_trafos[0];
      vector<double> aligned_diffs;
      trafo_external_.getDeviations(aligned_diffs);
      Size index = std::max(Size(0), Size(rt_quantile_ * static_cast<double>(aligned_diffs.size())) - 1);
      rt_uncertainty = aligned_diffs[index];
      try
      {
        aligner_trafos[0].fitModel("lowess");
        trafo_external_ = aligner_trafos[0];
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to align RTs of internal/external peptides. RT information will not be considered in the SVM classification. The original error message was:\n" << e.what() << endl;
      }
    }

    if (rt_window_ == 0.0)
    {
      // calculate RT window based on other parameters and alignment quality:
      double map_tol = mapping_tolerance_;
      if (map_tol < 1.0)
      {
        map_tol *= (2 * peak_width_); // relative tolerance
      }
      rt_window_ = (rt_uncertainty + 2 * peak_width_ + map_tol) * 2;
      OPENMS_LOG_INFO << "RT window size calculated as " << rt_window_ << " seconds."
               << endl;
    }

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
      addPeptideToMap_(pep, peptide_map_);
      pep.setMetaValue("FFId_category", "internal");
    }

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
    for (PeptideIdentification& pep : peptides_ext)
    {
      addPeptideToMap_(pep, peptide_map_, true);
      pep.setMetaValue("FFId_category", "external");
    }
    n_external_peps_ = peptide_map_.size() - n_internal_peps_;

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


        extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
                                       mz_window_ppm_, "tophat");
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
    for (Feature& f : features)
    {
      std::vector<PeptideIdentification>& ids = f.getPeptideIdentifications();

      // if we have peptide identifications assigned and all are annotated as OffsetPeptide, we mark the feature is also an OffsetPeptide
      if (!ids.empty() && std::all_of(ids.begin(), ids.end(), [](const PeptideIdentification & pid) { return pid.metaValueExists("OffsetPeptide"); }))
      {
        f.setMetaValue("OffsetPeptide", "true");
      }

      // remove all hits (PSM details)
      for (auto & pid : ids)
      {
        std::vector<PeptideHit>& hits = pid.getHits();
        auto it = remove_if(hits.begin(), hits.end(),
          [](const PeptideHit & ph)
          {
            return (ph.getSequence().toUnmodifiedString().hasPrefix("XXX"));
          });
        hits.erase(it, hits.end()); // remove / erase idiom
      }

      // remove empty PeptideIdentifications
      auto it = remove_if(ids.begin(), ids.end(),
        [](const PeptideIdentification & pid)
        {
          return pid.empty();
        });
      ids.erase(it, ids.end()); // remove / erase idiom
    }

    // clean up unassigned PeptideIdentifications
    std::vector<PeptideIdentification>& ids = features.getUnassignedPeptideIdentifications();
    for (auto & pid : ids)
    {
      std::vector<PeptideHit>& hits = pid.getHits();
      auto it = remove_if(hits.begin(), hits.end(), [](const PeptideHit & ph)
      {
        return (ph.getSequence().toUnmodifiedString().hasPrefix("XXX"));
      });
      hits.erase(it, hits.end());
    }

    // remove empty PeptideIdentifications
    auto it = remove_if(ids.begin(), ids.end(),
      [](const PeptideIdentification & pid)
      {
        return pid.empty();
      });
    ids.erase(it, ids.end()); // remove / erase idiom

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
      classifyFeatures_(features);
    }
    // make sure proper unique ids get assigned to all features
    features.ensureUniqueId();

    // store feature candidates before filtering
    if (!candidates_out_.empty())
    {
      FileHandler().storeFeatures(candidates_out_, features);
    }

    filterFeatures_(features, with_external_ids);
    OPENMS_LOG_INFO << features.size() << " features left after filtering." << endl;

    if (features.empty()) return; // elution model fit throws on empty features

    if (!svm_probs_internal_.empty())
    {
      calculateFDR_(features);
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

  void FeatureFinderIdentificationAlgorithm::createAssayLibrary_(const PeptideMap::iterator& begin, const PeptideMap::iterator& end, PeptideRefRTMap& ref_rt_map, bool clear_IDs)
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
        // internal/external pair
        const pair<RTMap, RTMap> &pair = pm_it->second.begin()->second;

        // WARNING: This assumes that at least one hit is present.
        const PeptideHit &hit = (pair.first.empty() ?
                                 pair.second.begin()->second->getHits()[0] :
                                 pair.first.begin()->second->getHits()[0]);
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
        // get regions in which peptide eludes (ideally only one):
        std::vector<RTRegion> rt_regions;
        getRTRegions_(pm_it->second, rt_regions, clear_IDs);

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

  void FeatureFinderIdentificationAlgorithm::getRTRegions_(
    ChargeMap& peptide_data,
    std::vector<RTRegion>& rt_regions,
    bool clear_IDs) const
  {
    // use RTs from all charge states here to get a more complete picture:
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
      // create a new region?
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
      auto reg_it = rt_regions.begin();
      // "internal" IDs:
      for (auto& rt : cm.second.first)
      {
        while (rt.first > reg_it->end)
        {
          ++reg_it;
        }
        reg_it->ids[cm.first].first.insert(rt);
      }
      reg_it = rt_regions.begin(); // reset to start
      // "external" IDs:
      for (auto& rt : cm.second.second)
      {
        while (rt.first > reg_it->end)
        {
          ++reg_it;
        }
        reg_it->ids[cm.first].second.insert(rt);
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

  void FeatureFinderIdentificationAlgorithm::checkNumObservations_(Size n_pos, Size n_neg, const String& note) const
  {
    if (n_pos < svm_n_parts_)
    {
      String msg = "Not enough positive observations for " + 
        String(svm_n_parts_) + "-fold cross-validation" + note + ".";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    if (n_neg < svm_n_parts_)
    {
      String msg = "Not enough negative observations for " + 
        String(svm_n_parts_) + "-fold cross-validation" + note + ".";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
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
        if (mapping_tolerance_ > 0.0)
        {
          double abs_tol = mapping_tolerance_;
          if (abs_tol < 1.0)
          {
            abs_tol *= (rt_max - rt_min);
          }
          rt_min -= abs_tol;
          rt_max += abs_tol;
        }
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
        feat.setMetaValue("n_total_ids", 0);
        feat.setMetaValue("n_matching_ids", -1);
        feat.setMetaValue("feature_class", "unknown");
        // add "dummy" peptide identification:
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

      // distance from feature to closest peptide ID:
      if (!trafo_external_.getDataPoints().empty())
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
            double transformed_rt = trafo_external_.apply(it->first);
            RTMap::value_type pair = make_pair(transformed_rt, it->second);
            transformed_internal.insert(transformed_internal.end(), pair);
          }
        }
        const RTMap& rt_ref = (rt_external.empty() ? transformed_internal :
                               rt_external);

        double rt_min = feat.getMetaValue("leftWidth");
        double rt_max = feat.getMetaValue("rightWidth");
        if (mapping_tolerance_ > 0.0)
        {
          double abs_tol = mapping_tolerance_;
          if (abs_tol < 1.0)
          {
            abs_tol *= (rt_max - rt_min);
          }
          rt_min -= abs_tol;
          rt_max += abs_tol;
        }
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
      if (hit.metaValueExists("target_decoy") && hit.getMetaValue("target_decoy") == "decoy")
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
    double rt = peptide.getRT();
    double mz = peptide.getMZ();
    if (!external)
    {
      if (peptide.metaValueExists("SeedFeatureID"))
      {
        OPENMS_LOG_DEBUG_NOFILE << "Adding seed (internal) from FeatureID " << peptide.getMetaValue("SeedFeatureID") << ": " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      }
      else
      {
        OPENMS_LOG_DEBUG_NOFILE << "Adding peptide (internal) " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      }
      peptide_map[hit.getSequence()][charge].first.emplace(rt, &peptide);
    }
    else
    {
      OPENMS_LOG_DEBUG_NOFILE << "Adding peptide (external) " << hit.getSequence() << "; CHG: " << charge << "; RT: " << rt << "; MZ: " << mz << endl;
      peptide_map[hit.getSequence()][charge].second.emplace(rt, &peptide);
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

  void FeatureFinderIdentificationAlgorithm::getUnbiasedSample_(const multimap<double, pair<Size, bool> >& valid_obs,
                          map<Size, double>& training_labels)
  {
    // Create an unbiased training sample:
    // - same number of pos./neg. observations (approx.),
    // - same intensity distribution of pos./neg. observations.
    // We use a sliding window over the set of observations, ordered by
    // intensity. At each step, we examine the proportion of both pos./neg.
    // observations in the window and select the middle element with according
    // probability. (We use an even window size, to cover the ideal case where
    // the two classes are balanced.)
    const Size window_size = 8;
    const Size half_win_size = window_size / 2;
    if (valid_obs.size() < half_win_size + 1)
    {
      String msg = "Not enough observations for intensity-bias filtering.";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    srand(time(nullptr)); // seed random number generator
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    Size counts[2] = {0, 0}; // pos./neg. counts in current window
    // iterators to begin, middle and past-the-end of sliding window:
    multimap<double, pair<Size, bool> >::const_iterator begin, middle, end;
    begin = middle = end = valid_obs.begin();
    // initialize ("middle" is at beginning of sequence, so no full window):
    for (Size i = 0; i <= half_win_size; ++i, ++end)
    {
      ++counts[end->second.second]; // increase counter for pos./neg. obs.
    }
    // "i" is the index of one of the two middle values of the sliding window:
    // - in the left half of the sequence, "i" is left-middle,
    // - in the right half of the sequence, "i" is right-middle.
    // The counts are updated as "i" and the sliding window move to the right.
    for (Size i = 0; i < valid_obs.size(); ++i, ++middle)
    {
      // if count for either class is zero, we don't select anything:
      if ((counts[0] > 0) && (counts[1] > 0))
      {
        // probability thresholds for neg./pos. observations:
        double thresholds[2] = {counts[1] / float(counts[0]),
                                counts[0] / float(counts[1])};
        // check middle values:
        double rnd = rand() / double(RAND_MAX); // random num. in range 0-1
        if (rnd < thresholds[middle->second.second])
        {
          training_labels[middle->second.first] = Int(middle->second.second);
          ++n_obs[middle->second.second];
        }
      }
      // update sliding window and class counts;
      // when we reach the middle of the sequence, we keep the window in place
      // for one step, to change from "left-middle" to "right-middle":
      if (i != valid_obs.size() / 2)
      {
        // only move "begin" when "middle" has advanced far enough:
        if (i > half_win_size)
        {
          --counts[begin->second.second];
          ++begin;
        }
        // don't increment "end" beyond the defined range:
        if (end != valid_obs.end())
        {
          ++counts[end->second.second];
          ++end;
        }
      }
    }
    checkNumObservations_(n_obs[1], n_obs[0], " after bias filtering");
  }


  void FeatureFinderIdentificationAlgorithm::getRandomSample_(std::map<Size, double>& training_labels) const
  {
    // @TODO: can this be done with less copying back and forth of data?
    // Pick a random subset of size "svm_n_samples_" for training: Shuffle the whole
    // sequence, then select the first "svm_n_samples_" elements.
    std::vector<Size> selection;
    selection.reserve(training_labels.size());
    for (auto it = training_labels.begin(); it != training_labels.end(); ++it)
    {
      selection.push_back(it->first);
    }
    //TODO check how often this is potentially called and move out the initialization
    Math::RandomShuffler shuffler;
    shuffler.portable_random_shuffle(selection.begin(), selection.end());
    // However, ensure that at least "svm_n_parts_" pos./neg. observations are
    // included (for cross-validation) - there must be enough, otherwise
    // "checkNumObservations_" would have thrown an error. To this end, move
    // "svm_n_parts_" pos. observations to the beginning of sequence, followed by
    // "svm_n_parts_" neg. observations (pos. first - see reason below):
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Int label = 1; label >= 0; --label)
    {
      for (Size i = n_obs[1]; i < selection.size(); ++i)
      {
        Size obs_index = selection[i];
        if (training_labels[obs_index] == label)
        {
          std::swap(selection[i], selection[n_obs[label]]);
          ++n_obs[label];
        }
        if (n_obs[label] == svm_n_parts_)
        {
          break;
        }
      }
    }
    selection.resize(svm_n_samples_);
    // copy the selected subset back:
    std::map<Size, double> temp;
    for (vector<Size>::iterator it = selection.begin(); it != selection.end();
         ++it)
    {
      temp[*it] = training_labels[*it];
    }
    training_labels.swap(temp);
  }

  void FeatureFinderIdentificationAlgorithm::classifyFeatures_(FeatureMap& features)
  {
    if (features.empty())
    {
      return;
    }
    if (features[0].metaValueExists("rt_delta")) // include RT feature
    {
      if (std::find(svm_predictor_names_.begin(), svm_predictor_names_.end(), "rt_delta") == svm_predictor_names_.end())
      {
        svm_predictor_names_.push_back("rt_delta");
      }
    }
    // values for all features per predictor (this way around to simplify scaling
    // of predictors):
    SimpleSVM::PredictorMap predictors;
    for (const String& pred : svm_predictor_names_)
    {
      predictors[pred].reserve(features.size());
      for (Feature& feat : features)
      {
        if (!feat.metaValueExists(pred))
        {
          OPENMS_LOG_ERROR << "Meta value '" << pred << "' missing for feature '"
                    << feat.getUniqueId() << "'" << endl;
          predictors.erase(pred);
          break;
        }
        predictors[pred].push_back(feat.getMetaValue(pred));
      }
    }

    // get labels for SVM:
    std::map<Size, double> training_labels;
    bool no_selection = param_.getValue("svm:no_selection") == "true";
    // mapping (for bias correction): intensity -> (index, positive?)
    std::multimap<double, pair<Size, bool> > valid_obs;
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      int label = -1;
      if (feature_class == "positive")
      {
        label = 1;
      }
      else if (feature_class == "negative")
      {
        label = 0;
      }
      if (label != -1)
      {
        ++n_obs[label];
        if (!no_selection)
        {
          double intensity = features[feat_index].getIntensity();
          valid_obs.insert(make_pair(intensity, make_pair(feat_index,
                                                          bool(label))));
        }
        else
        {
          training_labels[feat_index] = (double)label;
        }
      }
    }
    checkNumObservations_(n_obs[1], n_obs[0]);

    if (!no_selection)
    {
      getUnbiasedSample_(valid_obs, training_labels);
    }
    if (svm_n_samples_ > 0) // limited number of samples for training
    {
      if (training_labels.size() < svm_n_samples_)
      {
        OPENMS_LOG_WARN << "Warning: There are only " << training_labels.size()
                 << " valid observations for training." << endl;
      }
      else if (training_labels.size() > svm_n_samples_)
      {
        getRandomSample_(training_labels);
      }
    }

    SimpleSVM svm;
    // set (only) the relevant parameters:
    Param svm_params = svm.getParameters();
    Logger::LogStream no_log; // suppress warnings about additional parameters
    svm_params.update(param_.copy("svm:", true), false, no_log);
    svm.setParameters(svm_params);
    svm.setup(predictors, training_labels);
    if (!svm_xval_out_.empty())
    {
      svm.writeXvalResults(svm_xval_out_);
    }
    if ((debug_level_ > 0) && svm_params.getValue("kernel") == "linear")
    {
      std::map<String, double> feature_weights;
      svm.getFeatureWeights(feature_weights);
      OPENMS_LOG_DEBUG << "SVM feature weights:" << endl;
      for (std::map<String, double>::iterator it = feature_weights.begin();
           it != feature_weights.end(); ++it)
      {
        OPENMS_LOG_DEBUG << "- " << it->first << ": " << it->second << endl;
      }
    }

    std::vector<SimpleSVM::Prediction> predictions;
    svm.predict(predictions);
    OPENMS_POSTCONDITION(predictions.size() == features.size(), 
                         "SVM predictions for all features expected");
    for (Size i = 0; i < features.size(); ++i)
    {
      features[i].setMetaValue("predicted_class", predictions[i].outcome);
      double prob_positive = predictions[i].probabilities[1];
      features[i].setMetaValue("predicted_probability", prob_positive);
      // @TODO: store previous (OpenSWATH) overall quality in a meta value?
      features[i].setOverallQuality(prob_positive);
    }
  }


  void FeatureFinderIdentificationAlgorithm::filterFeaturesFinalizeAssay_(Feature& best_feature, double best_quality,
                                    const double quality_cutoff)
  {
    const String& feature_class = best_feature.getMetaValue("feature_class");
    if (feature_class == "positive") // true positive prediction
    {
      svm_probs_internal_[best_quality].first++;
    }
    else if ((feature_class == "negative")  || // false positive prediction
             (feature_class == "ambiguous")) // let's be strict about this
    {
      svm_probs_internal_[best_quality].second++;
    }
    else if (feature_class == "unknown")
    {
      svm_probs_external_.insert(best_quality);
      if (best_quality >= quality_cutoff) 
      {
        best_feature.setOverallQuality(best_quality);
        ++n_external_features_;
      }
    }
  }

  void FeatureFinderIdentificationAlgorithm::filterFeatures_(FeatureMap& features, bool classified)
  {
    if (features.empty())
    {
      return;
    }
    if (classified)
    {
      // Remove features with class "negative" or "ambiguous", keep "positive".
      // For class "unknown", for every assay (meta value "PeptideRef"), keep
      // the feature with highest "predicted_probability" (= overall quality),
      // subject to the "svm:min_prob" threshold.
      // We mark features for removal by setting their overall quality to zero.
      n_internal_features_ = n_external_features_ = 0;
      FeatureMap::Iterator best_it = features.begin();
      double best_quality = 0.0;
      String previous_ref;
      for (FeatureMap::Iterator it = features.begin(); it != features.end();
           ++it)
      {
        // features from same assay (same "PeptideRef") appear consecutively;
        // if this is a new assay, finalize the previous one:
        String peptide_ref = it->getMetaValue("PeptideRef");
        // remove region number, if present:
        Size pos_slash = peptide_ref.rfind('/');
        Size pos_colon = peptide_ref.find(':', pos_slash + 2);
        peptide_ref = peptide_ref.substr(0, pos_colon);

        if (peptide_ref != previous_ref)
        {
          if (!previous_ref.empty())
          {
            filterFeaturesFinalizeAssay_(*best_it, best_quality,
                                         svm_quality_cutoff);
            best_quality = 0.0;
          }
          previous_ref = peptide_ref;
        }

        // update qualities:
        if ((it->getOverallQuality() > best_quality) ||
            // break ties by intensity:
            ((it->getOverallQuality() == best_quality) &&
             (it->getIntensity() > best_it->getIntensity())))
        {
          best_it = it;
          best_quality = it->getOverallQuality();
        }
        if (it->getMetaValue("feature_class") == "positive")
        {
          n_internal_features_++;
        }
        else
        {
          it->setOverallQuality(0.0); // gets overwritten for "best" candidate
        }
      }
      // set of features from the last assay:
      filterFeaturesFinalizeAssay_(*best_it, best_quality, svm_quality_cutoff);

      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_quality_), features.end());
    }
    else
    {
      // remove features without ID (or pseudo ID from seeds)
      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_peptides_), features.end());
    }
  }


  void FeatureFinderIdentificationAlgorithm::calculateFDR_(FeatureMap& features)
  {
    // cumulate the true/false positive counts, in decreasing probability order:
    Size n_false = 0, n_true = 0;
    for (std::map<double, pair<Size, Size> >::reverse_iterator prob_it =
           svm_probs_internal_.rbegin(); prob_it != svm_probs_internal_.rend();
         ++prob_it)
    {
      n_true += prob_it->second.first;
      n_false += prob_it->second.second;
      prob_it->second.first = n_true;
      prob_it->second.second = n_false;
    }

    // print FDR for features that made the cut-off:
    std::map<double, pair<Size, Size> >::iterator prob_it =
      svm_probs_internal_.lower_bound(svm_min_prob_);
    if (prob_it != svm_probs_internal_.end())
    {
      float fdr = float(prob_it->second.second) / (prob_it->second.first +
                                                   prob_it->second.second);
      OPENMS_LOG_INFO << "Estimated FDR of features detected based on 'external' IDs: "
               << fdr * 100.0 << "%" << endl;
      fdr = (fdr * n_external_features_) / (n_external_features_ + 
                                            n_internal_features_);
      OPENMS_LOG_INFO << "Estimated FDR of all detected features: " << fdr * 100.0
               << "%" << endl;
    }

    // calculate q-values:
    std::vector<double> qvalues;
    qvalues.reserve(svm_probs_internal_.size());
    double min_fdr = 1.0;
    for (prob_it = svm_probs_internal_.begin();
         prob_it != svm_probs_internal_.end(); ++prob_it)
    {
      double fdr = double(prob_it->second.second) / (prob_it->second.first +
                                                     prob_it->second.second);
      if (fdr < min_fdr)
      {
        min_fdr = fdr;
      }
      qvalues.push_back(min_fdr);
    }
    // record only probabilities where q-value changes:
    std::vector<double> fdr_probs, fdr_qvalues;
    std::vector<double>::iterator qv_it = qvalues.begin();
    double previous_qvalue = -1.0;
    for (prob_it = svm_probs_internal_.begin();
         prob_it != svm_probs_internal_.end(); ++prob_it, ++qv_it)
    {
      if (*qv_it != previous_qvalue)
      {
        fdr_probs.push_back(prob_it->first);
        fdr_qvalues.push_back(*qv_it);
        previous_qvalue = *qv_it;
      }
    }
    features.setMetaValue("FDR_probabilities", fdr_probs);
    features.setMetaValue("FDR_qvalues_raw", fdr_qvalues);

    // FDRs are estimated from "internal" features, but apply only to "external"
    // ones. "Internal" features are considered "correct" by definition.
    // We need to adjust the q-values to take this into account:
    std::multiset<double>::reverse_iterator ext_it = svm_probs_external_.rbegin();
    Size external_count = 0;
    for (Int i = fdr_probs.size() - 1; i >= 0; --i)
    {
      double cutoff = fdr_probs[i];
      while ((ext_it != svm_probs_external_.rend()) && (*ext_it >= cutoff))
      {
        ++external_count;
        ++ext_it;
      }
      fdr_qvalues[i] = (fdr_qvalues[i] * external_count) /
        (external_count + n_internal_features_);
    }
    features.setMetaValue("FDR_qvalues_corrected", fdr_qvalues);

    // @TODO: should we use "1 - qvalue" as overall quality for features?
    // assign q-values to features:
    for (Feature& feat : features)
    {
      if (feat.getMetaValue("feature_class") == "positive")
      {
        feat.setMetaValue("q-value", 0.0);
      }
      else
      {
        double prob = feat.getOverallQuality();
        // find the highest FDR prob. that is less-or-equal to the feature prob.:
        std::vector<double>::iterator pos = upper_bound(fdr_probs.begin(),
                                                   fdr_probs.end(), prob);
        if (pos != fdr_probs.begin())
        {
          --pos;
        }
        Size dist = distance(fdr_probs.begin(), pos);
        feat.setMetaValue("q-value", fdr_qvalues[dist]);
      }
    }
  }
}
