// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h> // for legacy "run" method
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS
{
  FeatureFinderIdentificationAlgorithm::FeatureFinderIdentificationAlgorithm() :
    DefaultParamHandler("FeatureFinderIdentificationAlgorithm"),
    n_internal_targets_(0), n_external_targets_(0), n_seed_targets_(0)
  {
    std::vector<std::string> output_file_tags;
    output_file_tags.push_back("output file");

    defaults_.setValue("candidates_out", "", "Optional output file with feature candidates.", output_file_tags);

    defaults_.setValue("debug", 0, "Debug level for feature detection.", {"advanced"});
    defaults_.setMinInt("debug", 0);

    defaults_.setValue("extract:batch_size", 5000, "Number of target molecules to consider in each batch of chromatogram extraction."
                         " Smaller values decrease memory usage but increase runtime.");
    defaults_.setMinInt("extract:batch_size", 1);
    defaults_.setValue("extract:mz_window", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)");
    defaults_.setMinFloat("extract:mz_window", 0.0);
    defaults_.setValue("extract:n_isotopes", 2, "Number of isotopes to include in each target molecule assay.");
    defaults_.setMinInt("extract:n_isotopes", 2);
    defaults_.setValue("extract:max_isotopes", "false", "Be default, isotopes are extracted starting from the monoisotope. Set this flag to extract the most abundant isotopes instead. (Not relevant if advanced parameter 'extract:isotope_pmin' is used.)");
    defaults_.setValidStrings("extract:max_isotopes", {"true", "false"});
    defaults_.setValue("extract:isotope_pmin", 0.0, "Minimum probability for an isotope to be included in the assay for a molecule. If set, this parameter takes precedence over 'extract:n_isotopes'.", {"advanced"});
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
    defaults_.setValue("detect:mapping_tolerance", 0.0, "RT tolerance (plus/minus) for mapping molecule IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.");
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

  void FeatureFinderIdentificationAlgorithm::run(
    FeatureMap& features,
    IdentificationData& id_data,
    IdentificationData& id_data_ext,
    const String& spectra_file)
  {
    target_map_.clear();

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
    params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width",
                    peak_width_);
    params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width_);
    // disabling the signal-to-noise threshold (setting the parameter to zero)
    // totally breaks the OpenSWATH feature detection (no features found)!
    params.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise",
                    signal_to_noise_);
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method",
                    "corrected");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages", "false"); // disabled in OpenSWATH

    feat_finder_.setParameters(params);
    feat_finder_.setLogType(ProgressLogger::NONE);
    feat_finder_.setStrictFlag(false);
    // to use MS1 Swath scores:
    feat_finder_.setMS1Map(SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(boost::make_shared<MSExperiment>(ms_data_)));

    double rt_uncertainty = 0;
    bool with_external_ids = !id_data_ext.empty();
    if (with_external_ids)
    {
      // align internal and external IDs to estimate RT shifts:
      MapAlignmentAlgorithmIdentification aligner;
      aligner.setReference(id_data_ext); // go from internal to external scale
      vector<IdentificationData> aligner_ids(1, id_data);
      vector<TransformationDescription> aligner_trafos;

      OPENMS_LOG_INFO << "Realigning internal and external IDs...";
      aligner.align(aligner_ids, aligner_trafos);
      trafo_external_ = aligner_trafos[0];
      vector<double> aligned_diffs;
      trafo_external_.getDeviations(aligned_diffs);
      int index = max(0, int(rt_quantile_ * aligned_diffs.size()) - 1);
      rt_uncertainty = aligned_diffs[index];
      try
      {
        aligner_trafos[0].fitModel("lowess");
        trafo_external_ = aligner_trafos[0];
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: Failed to align RTs of internal/external IDs. RT information will not be considered in the SVM classification. The original error message was:\n" << e.what() << endl;
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
      OPENMS_LOG_INFO << "RT window size calculated as " << rt_window_ << " seconds." << endl;
    }

    //-------------------------------------------------------------
    // prepare target ion map
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Preparing data..." << endl;
    target_map_.clear();

    // @TODO: expose score choice to user via a parameter
    ID::ScoreTypeRef score_ref = id_data.pickScoreType(id_data.getObservationMatches());
    IDFilter::keepBestMatchPerObservation(id_data, score_ref);
    for (ID::ObservationMatchRef ref = id_data.getObservationMatches().begin();
         ref != id_data.getObservationMatches().end(); ++ref)
    {
      const String& cat = ref->getMetaValue("FFId_category", "");
      if (cat == "seed")
      {
        n_seed_targets_++;
      }
      else if (cat.empty()) // default to "internal"
      {
        id_data.setMetaValue(ref, "FFId_category", "internal");
      }
      addMatchToTargetMap_(ref);
    }
    n_internal_targets_ = target_map_.size() - n_seed_targets_;

    if (with_external_ids)
    {
      // find the same score type in the external data:
      String score_name = score_ref->cv_term.getName();
      score_ref = id_data_ext.findScoreType(score_name);
      if (score_ref == id_data_ext.getScoreTypes().end()) // not found
      {
        String msg = "score type '" + score_name + "' not found in external IDs";
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
      }
      IDFilter::keepBestMatchPerObservation(id_data_ext, score_ref);
      for (ID::ObservationMatchRef ref = id_data_ext.getObservationMatches().begin();
         ref != id_data_ext.getObservationMatches().end(); ++ref)
      {
        addMatchToTargetMap_(ref, true);
        id_data_ext.setMetaValue(ref, "FFId_category", "external");
      }
    }
    n_external_targets_ = target_map_.size() - n_internal_targets_ - n_seed_targets_;

    boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
    OpenSwath::SpectrumAccessPtr spec_temp = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
    auto chunks = chunk_(target_map_.begin(), target_map_.end(), batch_size_);

    if (debug_level_ >= 666)
    {
      OPENMS_LOG_INFO << "Creating full assay library for debugging." << endl;
      // Warning: this step is pretty inefficient, since it does the whole library generation twice
      // Really use for debug only
      createAssayLibrary_(target_map_.begin(), target_map_.end());
      OPENMS_LOG_DEBUG << "Writing debug.traml file..." << endl;
      TraMLFile().store("debug.traml", library_);
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
      OPENMS_LOG_INFO << "Creating assay library..." << endl;
      createAssayLibrary_(chunk.first, chunk.second);
      OPENMS_LOG_DEBUG << "#Transitions: " << library_.getTransitions().size() << endl;

      OPENMS_LOG_INFO << "Extracting chromatograms..." << endl;
      ChromatogramExtractor extractor;
      // extractor.setLogType(ProgressLogger::NONE);
      {
        vector<OpenSwath::ChromatogramPtr> chrom_temp;
        vector<ChromatogramExtractor::ExtractionCoordinates> coords;
        // take entries in library_ and write to chrom_temp and coords
        extractor.prepare_coordinates(chrom_temp, coords, library_,
                                      numeric_limits<double>::quiet_NaN(), false);
        extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
                                       mz_window_ppm_, "tophat");
        extractor.return_chromatogram(chrom_temp, coords, library_, (*shared)[0],
                                      chrom_data_.getChromatograms(), false);
      }

      OPENMS_LOG_INFO << "Extracted " << chrom_data_.getNrChromatograms()
                      << " chromatogram(s)." << endl;

      OPENMS_LOG_INFO << "Detecting chromatographic peaks..." << endl;
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
      combined_library_ += library_;
      library_.clear(true);
      // since chrom_data_ here is just a container for the chromatograms and identifications will be empty,
      // pickExperiment above will only add empty ProteinIdentification runs with colliding identifiers.
      // Usually we could sanitize the identifiers or merge the runs, but since they are empty and we add the
      // "real" proteins later -> just clear them
      features.getProteinIdentifications().clear();
      getProgressLogger().setProgress(++chunk_count);
    }
    getProgressLogger().endProgress();

    OPENMS_LOG_INFO << "Found " << features.size() << " feature candidates in total." << endl;
    // if (!candidates_out_.empty()) FeatureXMLFile().store(candidates_out_, features);

    ms_data_.reset(); // not needed anymore, free up the memory

    // complete feature annotation:
    annotateFeatures_(features);
    // sort everything:
    sort(features.begin(), features.end(), feature_compare_);

    postProcess_(features, with_external_ids);
    statistics_(features, with_external_ids);

    // move IDs:
    IdentificationData& feat_ids = features.getIdentificationData();
    feat_ids.swap(id_data);
    if (with_external_ids)
    {
      ID::RefTranslator trans = feat_ids.merge(id_data_ext);
      // references to internal IDs don't need to be translated:
      trans.allow_missing = true;
      for (Feature& feature : features)
      {
        feature.updateAllIDReferences(trans);
      }
    }
  }

  // represent seeds in IdentificationData format
  void FeatureFinderIdentificationAlgorithm::convertSeeds(
    const FeatureMap& seeds, IdentificationData& id_data,
    const String& input_file, Size n_overlap_traces)
  {
    // TODO make sure that only assembled traces (more than one trace -> has a charge) are used
    // see FeatureFindingMetabo: defaults_.setValue("remove_single_traces", "false", "Remove unassembled traces (single traces).");

    ID::ProcessingSoftware software("FeatureFinderIdentification", VersionInfo::getVersion());
    ID::ProcessingSoftwareRef sw_ref = id_data.registerProcessingSoftware(software);
    ID::InputFile input(input_file);
    if (input.name.empty())
    {
      input.name = seeds.getLoadedFilePath();
    }
    if (input.name.empty())
    {
      input.name = "UNKNOWN_SEEDS_INPUT";
    }
    ID::InputFileRef file_ref = id_data.registerInputFile(input);
    ID::ProcessingStep step(sw_ref, {file_ref});
    ID::ProcessingStepRef step_ref = id_data.registerProcessingStep(step);
    id_data.setCurrentProcessingStep(step_ref);

    Size seeds_added(0);
    vector<bool> target_already_exists(seeds.size(), false);
    Size feature_index = 0;
    double rt_tolerance = rt_window_seeds_ / 2.0;
    for (const Feature& seed : seeds)
    {
      double seed_rt = seed.getRT();
      double seed_mz = seed.getMZ();
      double seed_charge = seed.getCharge();

      // consider multiple isotopic traces when checking m/z overlap
      vector<double> isotopes_mz(n_overlap_traces, seed_mz);
      for (Size i = 1; i < n_overlap_traces; ++i)
      {
        isotopes_mz[i] -= double(i) / seed_charge * Constants::C13C12_MASSDIFF_U;
      }

      // check if there's already a target that is close in RT and MZ; if so, don't add seed
      // @TODO: this checks every observation for every seed; can we be more efficient?
      // (doing this after we have RT regions for the "real" IDs would be ideal)
      for (ID::ObservationRef ref = id_data.getObservations().begin();
           ref != id_data.getObservations().end(); ++ref)
      {
        auto pair = id_data.getMatchesForObservation(ref);
        if (pair.first == pair.second) continue; // no matches for this observation

        // RT or MZ values of seed match in range -> target already exists -> don't add seed
        double th_tolerance = mz_window_ppm_ ? mz_window_ * 1e-6 * ref->mz : mz_window_;
        if ((fabs(seed_rt - ref->rt) <= rt_tolerance) &&
            any_of(isotopes_mz.begin(), isotopes_mz.end(), [&](double isotope_mz)
            {
              return fabs(isotope_mz - ref->mz) <= th_tolerance;
            }))
        {
          target_already_exists[feature_index] = true;
          OPENMS_LOG_DEBUG_NOFILE << "Skipping seed from feature " << String(seed.getUniqueId())
                                  << " with z=" << seed_charge << ", RT=" << seed_rt << ", m/z=" << seed_mz
                                  << " due to overlap with observation " << ref->data_id
                                  << " at RT=" << ref->rt << ", m/z=" << ref->mz << endl;
          break;
        }
      }
      ++feature_index;
    }

    for (feature_index = 0; feature_index < seeds.size(); ++feature_index)
    {
      if (!target_already_exists[feature_index])
      {
        const Feature& seed = seeds[feature_index];
        ++seeds_added; // count from 1
        String id = String(seed.getUniqueId());
        ID::Observation query(id, file_ref, seed.getRT(), seed.getMZ());
        ID::ObservationRef query_ref = id_data.registerObservation(query);

        // represent seeds as compounds - no need to fake a peptide sequence:
        ID::IdentifiedCompound compound("SEED:" + String(seeds_added));
        ID::IdentifiedCompoundRef com_ref =
          id_data.registerIdentifiedCompound(compound);

        ID::ObservationMatch match(com_ref, query_ref, seed.getCharge());
        match.setMetaValue("FFId_category", "seed");
        ID::ObservationMatchRef match_ref =
          id_data.registerObservationMatch(match);
        addMatchToTargetMap_(match_ref);
      }
    }
    OPENMS_LOG_INFO << seeds_added << " seeds without RT and m/z overlap with existing IDs added" << endl;
  }

  pair<String, Int> FeatureFinderIdentificationAlgorithm::extractTargetID_(
    const Feature& feature, bool extract_charge)
  {
    String target_id = feature.getMetaValue("PeptideRef"); // e.g. "PEP:XXXXX/2#1"
    Size pos_slash = target_id.rfind('/');
    Int charge = 0;
    if (extract_charge)
    {
      Size pos_hash = target_id.find('#', pos_slash + 2);
      // second arg. to "substr" is "count", not "end_pos"!
      charge = target_id.substr(pos_slash + 1, pos_hash - pos_slash - 1).toInt();
    }
    return make_pair(target_id.substr(0, pos_slash), charge);
  }


  void FeatureFinderIdentificationAlgorithm::postProcess_(FeatureMap& features,
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

    // fill in meta data for the features that are left:
    for (auto& feature : features)
    {
      ensureConvexHulls_(feature);
      // annotate subordinates with theoretical isotope intensities:
      for (Feature& sub_feature : feature.getSubordinates())
      {
        String native_id = sub_feature.getMetaValue("native_id");
        sub_feature.setMetaValue("isotope_probability", isotope_probs_[native_id]);
      }
      // add label:
      String target_id = extractTargetID_(feature).first;
      target_id = target_id.substr(target_id.find(':') + 1); // remove type prefix
      feature.setMetaValue("label", target_id);
    }

    if (!svm_probs_internal_.empty()) calculateFDR_(features);

    // TODO: MRMFeatureFinderScoring already does an ElutionModel scoring. It uses EMG fitting.
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
      for (auto& feat : features)
      {
        for (auto& sub : feat.getSubordinates())
        {
          sub.getConvexHulls().clear();
        }
      }
    }
  }

/*
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
    molecule_map_.clear();
    set<AASequence> internal_seqs;
    for (PeptideIdentification& pep : features.getUnassignedPeptideIdentifications())
    {
      const AASequence& seq = pep.getHits()[0].getSequence();
      if (pep.getMetaValue("FFId_category") == "internal")
      {
        internal_seqs.insert(seq);
      }
      molecule_map_[seq];
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
      molecule_map_[seq];
    }
    n_internal_molecules_ = internal_seqs.size();
    n_external_molecules_ = molecule_map_.size() - internal_seqs.size();

    // sort everything:
    sort(features.getUnassignedPeptideIdentifications().begin(),
         features.getUnassignedPeptideIdentifications().end(),
         peptide_compare_);
    sort(features.begin(), features.end(), feature_compare_);

    postProcess_(features, with_external_ids);

    statistics_(features);
  }
*/

  void FeatureFinderIdentificationAlgorithm::statistics_(
    const FeatureMap& features, bool with_external_ids) const
  {
    // same molecule may be quantified based on internal and external IDs if
    // charge states differ!
    set<String> quantified_internal, quantified_seed, quantified_all;
    for (const Feature& feature : features)
    {
      String target_id = extractTargetID_(feature, false).first;
      if (feature.getIntensity() > 0.0)
      {
        quantified_all.insert(target_id);
        if (feature.getIDMatches().empty()) continue;
        String category = (*feature.getIDMatches().begin())->getMetaValue("FFId_category");
        // @TODO: we might need to introduce another category for transfers or even declare them external
        // but we do not always want to trigger the SVM due to runtime.
        if (category == "internal" || category == "transfer")
        {
          quantified_internal.insert(target_id);
        }
        else if (category == "seed")
        {
          quantified_seed.insert(target_id);
        }
      }
    }
    Size n_quant_external = quantified_all.size() - quantified_internal.size() -
      quantified_seed.size();
    // If internal and external IDs for a peptide map to different RT regions,
    // it is possible that there is a quantification from the "external" region,
    // but not from the "internal" region (no matching feature) - therefore the
    // number of "missing" external peptides can be negative!
    Int n_missing_external = Int(n_external_targets_) - n_quant_external;

    // @TODO: break targets down into peptides, compounds, RNA oligos?
    OPENMS_LOG_INFO << "\nSummary statistics (counting distinct targets including modifications):\n"
                    << target_map_.size() << " targets identified";
    if (with_external_ids || (n_seed_targets_ > 0))
    {
      OPENMS_LOG_INFO << " (" << n_internal_targets_ << " internal";
      if (with_external_ids)
      {
        OPENMS_LOG_INFO << ", " << n_external_targets_ << " additional external";
      }
      if (n_seed_targets_ > 0)
      {
        OPENMS_LOG_INFO << ", " << n_seed_targets_ << " seed-based";
      }
      OPENMS_LOG_INFO << ")";
    }
    OPENMS_LOG_INFO << "\n" << quantified_all.size() << " targets with features";
    if (with_external_ids || (n_seed_targets_ > 0))
    {
      OPENMS_LOG_INFO << " (" << quantified_internal.size() << " internal";
      if (with_external_ids)
      {
        OPENMS_LOG_INFO << ", " << n_quant_external << " external";
      }
      if (n_seed_targets_ > 0)
      {
        OPENMS_LOG_INFO << ", " << quantified_seed.size() << " seed-based";
      }
      OPENMS_LOG_INFO << ")";
    }
    OPENMS_LOG_INFO << "\n" << target_map_.size() - quantified_all.size()
                    << " targets without features";
    if (with_external_ids || (n_seed_targets_ > 0))
    {
      OPENMS_LOG_INFO << " (" << n_internal_targets_ - quantified_internal.size() << " internal";
      if (with_external_ids)
      {
        OPENMS_LOG_INFO << ", " << n_missing_external << " external";
        if (n_missing_external < 0) OPENMS_LOG_INFO << "*";
      }
      if (n_seed_targets_ > 0)
      {
        OPENMS_LOG_INFO << ", " << n_seed_targets_ - quantified_seed.size() << " seed-based";
      }
      OPENMS_LOG_INFO << ")";
    }
    if (n_missing_external < 0)
    {
      OPENMS_LOG_INFO << "\n*: At least " << -n_missing_external << " targets without valid features based on internal IDs were quantified based on external IDs (in a different charge state) instead.";
    }
    OPENMS_LOG_INFO << "\n" << endl;
  }


  void FeatureFinderIdentificationAlgorithm::createAssayLibrary_(
    TargetMap::iterator begin, TargetMap::iterator end)
  {
    Size n_isotopes = 10;
    if (!max_isotopes_ && (isotope_pmin_ == 0.0)) n_isotopes = n_isotopes_;
    CoarseIsotopePatternGenerator iso_gen(n_isotopes);
    iso_gen.setRoundMasses(false); // this is already the default, but be explicit

    for (auto target_it = begin; target_it != end; ++target_it)
    {
      const String& target_id = target_it->first;
      OPENMS_LOG_DEBUG << "\nTarget molecule: " << target_id << endl;
      const ID::IdentifiedMolecule& molecule = target_it->second.molecule;
      const ID::AdductOpt& adduct = target_it->second.adduct;
      ID::MoleculeType molecule_type = molecule.getMoleculeType();
      // special case: peptides can contain custom modifications that are
      // only defined by mass shift; these get lost during conversion to an
      // empirical formula, but are correctly included in the mass here:
      double mass = (molecule_type == ID::MoleculeType::PROTEIN) ?
        molecule.getIdentifiedPeptideRef()->sequence.getMonoWeight() : 0.0;

      // @TODO: use "TargetedExperiment::Peptide" for peptides?
      TargetedExperiment::Compound target;
      IsotopeDistribution iso_dist;
      EmpiricalFormula formula = molecule.getFormula();
      if (!formula.isEmpty())
      {
        target.theoretical_mass = (mass > 0.0) ? mass : formula.getMonoWeight();
        if (adduct)
        {
          target.theoretical_mass += (*adduct)->getMassShift();
          // when adding sum formulas, balance adduct charge with hydrogens:
          formula += (*adduct)->getEmpiricalFormula() -
            EmpiricalFormula::hydrogen((*adduct)->getCharge());
        }
        target.molecular_formula = formula.toString();
        iso_dist = formula.getIsotopeDistribution(iso_gen);
      }
      else // seed
      {
        if (!target_id.hasPrefix("SEED:"))
        {
          OPENMS_LOG_WARN
            << "Warning: no chemical formula given for compound '" << target_id
            << "'; estimating isotopic distribution based on peptide averagine" << endl;
        }
        // find m/z from an "input item" that we generated for the seed:
        ID::ObservationMatchRef match_ref = target_it->second.hits_by_charge.begin()->second.first.begin()->second;
        ID::ObservationRef obs_ref = match_ref->observation_ref;
        target.theoretical_mass = (obs_ref->mz - Constants::PROTON_MASS_U) * match_ref->charge;
        // @TODO: add support for RNA "averagine" option
        iso_dist = iso_gen.estimateFromPeptideWeight(target.theoretical_mass);
      }

      // mass difference between molecule and isotopic distribution can occur
      // when "anonymous" modifications are used in peptides (see comment above):
      double mass_diff = target.theoretical_mass - iso_dist[0].getMZ();
      double epsilon = numeric_limits<double>::epsilon() * target.theoretical_mass * 2;
      if (mass_diff > epsilon)
      {
        OPENMS_LOG_DEBUG << "Mass shift for target " << target_id << ": " << mass_diff << endl;
        for (Peak1D& pair : iso_dist)
        {
          pair.setMZ(pair.getMZ() + mass_diff);
        }
      }

      if (isotope_pmin_ > 0.0)
      {
        iso_dist.trimIntensities(isotope_pmin_);
        iso_dist.renormalize(); // @TODO: is this useful here?
      }
      else if (max_isotopes_)
      {
        iso_dist.sortByIntensity();
        iso_dist.resize(n_isotopes_);
      }

      // get regions in which target elutes (ideally only one):
      vector<RTRegion> rt_regions;
      bool is_seed = target_id.hasPrefix("SEED:");
      makeRTRegions_(target_it->second.hits_by_charge, rt_regions, is_seed);
      OPENMS_LOG_DEBUG << "Found " << rt_regions.size() << " RT region(s)." << endl;

      // go through different charge states:
      for (const auto& charge_pair : target_it->second.hits_by_charge)
      {
        Int charge = charge_pair.first; // note that charge may be negative
        target.setChargeState(charge);
        OPENMS_LOG_DEBUG << "Charge: " << charge << endl;
        String ion_id = target_id + "/" + String(charge);

        // we want to detect one feature per peptide and charge state - if there
        // are multiple RT regions, group them together:
        // peptide.setPeptideGroupLabel(peptide_id);
        target.rts.clear();
        Size counter = 0;
        for (auto& region : rt_regions)
        {
          if (region.ids.count(charge))
          {
            OPENMS_LOG_DEBUG
              << "Region " << counter + 1 << " (RT: " << float(region.start)
              << "-" << float(region.end) << ", size "
              << float(region.end - region.start) << ")" << endl;

            target.id = ion_id;
            if (rt_regions.size() > 1) target.id += "#" + String(++counter);

            // store beginning and end of RT region:
            target.rts.clear();
            addTargetRT_(target, region.start);
            addTargetRT_(target, region.end);
            library_.addCompound(target);
            generateTransitions_(target.id, target.theoretical_mass, charge, iso_dist);
          }
        }
      }
    }
  }


  void FeatureFinderIdentificationAlgorithm::makeRTRegions_(
    const ChargeMap& charge_data, vector<RTRegion>& rt_regions, bool is_seed) const
  {
    // use RTs from all charge states here to get a more complete picture:
    vector<double> rts;
    for (const auto& charge_pair : charge_data)
    {
      for (const auto& rt_pair : charge_pair.second.first) // "internal" IDs
      {
        rts.push_back(rt_pair.first);
      }
      for (const auto& rt_pair : charge_pair.second.second) // "external" IDs
      {
        rts.push_back(rt_pair.first);
      }
    }
    sort(rts.begin(), rts.end());
    double rt_tolerance = (is_seed ? rt_window_seeds_ : rt_window_) / 2.0;

    // create RT regions based on how close together the RTs are:
    for (double rt : rts)
    {
      // create a new region?
      if (rt_regions.empty() || (rt_regions.back().end < rt - rt_tolerance))
      {
        RTRegion region;
        region.start = rt - rt_tolerance;
        // cppcheck-suppress uninitStructMember
        rt_regions.push_back(region);
      }
      rt_regions.back().end = rt + rt_tolerance;
    }

    // sort the IDs into the regions:
    for (const auto& charge_pair : charge_data)
    {
      // regions are sorted by RT, as are IDs, so just iterate linearly:
      auto reg_it = rt_regions.begin();
      // "internal" IDs:
      for (auto rt_it = charge_pair.second.first.begin();
           rt_it != charge_pair.second.first.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[charge_pair.first].first.insert(*rt_it);
      }
      reg_it = rt_regions.begin(); // reset to start
      // "external" IDs:
      for (auto rt_it = charge_pair.second.second.begin();
           rt_it != charge_pair.second.second.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[charge_pair.first].second.insert(*rt_it);
      }
    }
  }


  void FeatureFinderIdentificationAlgorithm::addTargetRT_(
    TargetedExperiment::Compound& target, double rt) const
  {
    TargetedExperiment::RetentionTime te_rt;
    te_rt.setRT(rt);
    te_rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
    te_rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::LOCAL;
    target.rts.push_back(te_rt);
  }


  /// generate transitions (isotopic traces) for an ion and add them to the library:
  void FeatureFinderIdentificationAlgorithm::generateTransitions_(
    const String& target_id, double target_mass, Int charge, const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    double precursor_mz = (target_mass + charge * Constants::PROTON_MASS_U) / abs(charge);
    for (const auto& isotope : iso_dist)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter + 1);
      String transition_name = target_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(precursor_mz);
      // transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U *
      //                         float(counter) / abs(charge));
      transition.setProductMZ((isotope.getMZ() + charge * Constants::PROTON_MASS_U) / abs(charge));
      transition.setLibraryIntensity(isotope.getIntensity());
      transition.setMetaValue("annotation", annotation);
      transition.setCompoundRef(target_id);
      //TODO what about transition charge? A lot of DIA scores depend on it and default to charge 1 otherwise.
      library_.addTransition(transition);
      isotope_probs_[transition_name] = isotope.getIntensity();
      ++counter;
    }
  }


  void FeatureFinderIdentificationAlgorithm::checkNumObservations_(
    Size n_pos, Size n_neg, const String& note) const
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


  /// annotate identified features with m/z, isotope probabilities, etc.
  void FeatureFinderIdentificationAlgorithm::annotateFeatures_(FeatureMap& features)
  {
    // map indexes of features to the targets (ID, charge) that generated them:
    map<pair<String, Int>, vector<Size>> features_per_target;
    for (Size index = 0; index < features.size(); ++index)
    {
      Feature& feature = features[index];
      // to which target does this feature (candidate) belong?
      pair<String, Int> target_id_pair = extractTargetID_(feature, true);
      feature.setCharge(target_id_pair.second);
      feature.setMZ(feature.getMetaValue("PrecursorMZ"));
      // remove "fake" IDs generated by OpenSWATH (they would be removed with
      // a warning when writing output, because of missing protein
      // identification with corresponding identifier):
      feature.getPeptideIdentifications().clear();
      features_per_target[target_id_pair].push_back(index);
    }

    for (const auto& element : features_per_target)
    {
      annotateFeaturesOneTarget_(features, element.first.first,
                                 element.first.second, element.second);
    }
  }


  // process feature candidates from the same target together:
  void FeatureFinderIdentificationAlgorithm::annotateFeaturesOneTarget_(
    FeatureMap& features, const String& target_id, Int charge, const vector<Size>& indexes)
  {
    // @TODO: how much of this makes sense for features derived from seeds?
    RTMap transformed_internal;
    TargetData& target_data = target_map_.at(target_id);
    auto& rt_maps = target_data.hits_by_charge.at(charge);
    RTMap& rt_internal = rt_maps.first;
    RTMap& rt_external = rt_maps.second;

    // track which feature candidate (index) has most matching internal IDs:
    Size best_index = 0, best_count = 0;

    for (Size index : indexes)
    {
      Feature& feature = features[index];
      if (rt_internal.empty() && rt_external.empty()) // this shouldn't happen
      {
        OPENMS_LOG_DEBUG << "Valid target IDs ('PeptideRef'):" << endl;
        for (const auto& pair : target_map_)
        {
          OPENMS_LOG_DEBUG << pair.first << endl;
        }
        String msg = "No internal or external RTs for target '" + target_id +
          "', charge " + String(charge) + ", stored as '" +
          String(feature.getMetaValue("PeptideRef")) + "'";
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
      }

      // RT span of feature including tolerance for mapping IDs:
      double rt_min = 0, rt_max = 0;
      // calculate values only if needed below:
      if (!rt_internal.empty() || !trafo_external_.getDataPoints().empty())
      {
        rt_min = feature.getMetaValue("leftWidth");
        rt_max = feature.getMetaValue("rightWidth");
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
      }

      if (!rt_internal.empty()) // validate based on internal IDs
      {
        // map IDs to features based on RT:
        RTMap::const_iterator lower = rt_internal.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_internal.upper_bound(rt_max);
        Size id_count = 0;
        for (; lower != upper; ++lower)
        {
          feature.addIDMatch(lower->second);
          ++id_count;
        }
        // "total" only includes IDs from this RT region:
        feature.setMetaValue("n_total_ids", rt_internal.size());
        feature.setMetaValue("n_matching_ids", id_count);
        if (id_count > 0) // matching IDs -> feature may be "correct"
        {
          feature.setMetaValue("feature_class", "ambiguous");
          // is this the new "best" feature candidate (with most matching IDs)?
          if ((id_count > best_count) ||
              ((id_count == best_count) && // break ties by intensity
               (feature.getIntensity() > features[best_index].getIntensity())))
          {
            best_count = id_count;
            best_index = index;
          }
        }
        else // no matching IDs -> feature is wrong
        {
          feature.setMetaValue("feature_class", "negative");
        }
      }
      else // only external IDs -> no validation possible
      {
        feature.setMetaValue("n_total_ids", 0);
        feature.setMetaValue("n_matching_ids", -1);
        feature.setMetaValue("feature_class", "unknown");
      }

      // distance from feature to closest peptide ID (for training classifier):
      if (!trafo_external_.getDataPoints().empty())
      {
        // use external IDs if available, otherwise RT-transformed internal IDs
        // (but only compute the transform if necessary, once per assay!):
        if (rt_external.empty() && transformed_internal.empty())
        {
          transformed_internal.clear();
          for (const auto& rt_pair : rt_internal)
          {
            double transformed_rt = trafo_external_.apply(rt_pair.first);
            RTMap::value_type new_pair = make_pair(transformed_rt, rt_pair.second);
            transformed_internal.insert(transformed_internal.end(), new_pair);
          }
        }
        const RTMap& rt_ref = (rt_external.empty() ? transformed_internal : rt_external);
        RTMap::const_iterator lower = rt_ref.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_ref.upper_bound(rt_max);
        if (lower != upper) // there's at least one ID within the feature
        {
          feature.setMetaValue("rt_delta", 0.0);
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
          feature.setMetaValue("rt_delta", min(rt_delta1, rt_delta2));
        }
      }
    }

    // if internal IDs matched to feature candidates, select the best candidate:
    // (with the most IDs):
    if (best_count > 0)
    {
      // we define the (one) feature with most matching IDs as correct:
      features[best_index].setMetaValue("feature_class", "positive");
      features[best_index].setPrimaryID(target_data.molecule);
    }

    // clear internal IDs for this target (no longer needed):
    rt_internal.clear();
  }


  void FeatureFinderIdentificationAlgorithm::ensureConvexHulls_(Feature& feature)
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


  String makeTargetID(ID::ObservationMatchRef ref)
  {
    const ID::IdentifiedMolecule& molecule = ref->identified_molecule_var;
    const ID::AdductOpt& adduct = ref->adduct_opt;
    ID::MoleculeType molecule_type = molecule.getMoleculeType();
    String target_id;
    switch (molecule_type)
    {
      case ID::MoleculeType::PROTEIN:
        target_id = "PEP:" + molecule.getIdentifiedPeptideRef()->sequence.toString();
        break;
      case ID::MoleculeType::COMPOUND: // actual compound or seed?
        target_id = molecule.getIdentifiedCompoundRef()->identifier;
        if (!target_id.hasPrefix("SEED:")) target_id = "COM:" + target_id;
        break;
      case ID::MoleculeType::RNA:
        target_id = "RNA:" + molecule.getIdentifiedOligoRef()->sequence.toString();
        break;
      default: // avoid compiler warning
        break;
    }
    if (adduct) target_id += "+[" + (*adduct)->getName() + "]";
    return target_id;
  }


  void FeatureFinderIdentificationAlgorithm::addMatchToTargetMap_(
    ID::ObservationMatchRef ref, bool external)
  {
    String target_id = makeTargetID(ref);
    auto pos = target_map_.find(target_id);
    if (pos == target_map_.end()) // no entry for this target yet
    {
      if (!quantify_decoys_)
      {
        // check target/decoy status of molecule, skip decoys:
        const ID::IdentifiedMolecule& molecule = ref->identified_molecule_var;
        ID::MoleculeType type = molecule.getMoleculeType();
        if ((type == ID::MoleculeType::PROTEIN) &&
            (molecule.getIdentifiedPeptideRef()->allParentsAreDecoys())) return;
        if ((type == ID::MoleculeType::RNA) &&
            (molecule.getIdentifiedOligoRef()->allParentsAreDecoys())) return;
      }

      TargetData data;
      data.molecule = ref->identified_molecule_var;
      data.adduct = ref->adduct_opt;
      pos = target_map_.emplace(target_id, data).first;
    }

    double rt = ref->observation_ref->rt;
    if (!external)
    {
      OPENMS_LOG_DEBUG << "Adding " << target_id << " " << ref->charge << " " << rt << endl;
      pos->second.hits_by_charge[ref->charge].first.emplace(rt, ref);
    }
    else
    {
      pos->second.hits_by_charge[ref->charge].second.emplace(rt, ref);
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
    rt_window_seeds_ = 2.0 * peak_width_; // @TODO: add a parameter for this?

    n_isotopes_ = param_.getValue("extract:n_isotopes");
    max_isotopes_ = param_.getValue("extract:max_isotopes") == "true";
    isotope_pmin_ = param_.getValue("extract:isotope_pmin");

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
  }


  void FeatureFinderIdentificationAlgorithm::getUnbiasedSample_(
    const multimap<double, pair<Size, bool>>& valid_obs,
    map<Size, Int>& training_labels)
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


  void FeatureFinderIdentificationAlgorithm::getRandomSample_(map<Size, Int>& training_labels)
  {
    // @TODO: can this be done with less copying back and forth of data?
    // Pick a random subset of size "svm_n_samples_" for training: Shuffle the whole
    // sequence, then select the first "svm_n_samples_" elements.
    vector<Size> selection;
    selection.reserve(training_labels.size());
    for (map<Size, Int>::iterator it = training_labels.begin();
         it != training_labels.end(); ++it)
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
          swap(selection[i], selection[n_obs[label]]);
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
    map<Size, Int> temp;
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
      if (find(svm_predictor_names_.begin(), svm_predictor_names_.end(), "rt_delta") ==
          svm_predictor_names_.end())
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
    map<Size, Int> training_labels;
    bool no_selection = param_.getValue("svm:no_selection") == "true";
    // mapping (for bias correction): intensity -> (index, positive?)
    multimap<double, pair<Size, bool> > valid_obs;
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      Int label = -1;
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
          training_labels[feat_index] = label;
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
      map<String, double> feature_weights;
      svm.getFeatureWeights(feature_weights);
      OPENMS_LOG_DEBUG << "SVM feature weights:" << endl;
      for (map<String, double>::iterator it = feature_weights.begin();
           it != feature_weights.end(); ++it)
      {
        OPENMS_LOG_DEBUG << "- " << it->first << ": " << it->second << endl;
      }
    }

    vector<SimpleSVM::Prediction> predictions;
    svm.predict(predictions);
    OPENMS_POSTCONDITION(predictions.size() == features.size(),
                         "SVM predictions for all features expected");
    for (Size i = 0; i < features.size(); ++i)
    {
      features[i].setMetaValue("predicted_class", predictions[i].label);
      double prob_positive = predictions[i].probabilities[1];
      features[i].setMetaValue("predicted_probability", prob_positive);
      // @TODO: store previous (OpenSWATH) overall quality in a meta value?
      features[i].setOverallQuality(prob_positive);
    }
  }


  void FeatureFinderIdentificationAlgorithm::filterFeaturesFinalizeAssay_(
    Feature& best_feature, double best_quality, const double quality_cutoff,
    const String& target_id)
  {
    const String& feature_class = best_feature.getMetaValue("feature_class");
    if (feature_class == "positive") // true positive prediction
    {
      svm_probs_internal_[best_quality].first++;
      // feature primary ID is already set in this case
    }
    else if ((feature_class == "negative")  || // false positive prediction
             (feature_class == "ambiguous")) // let's be strict about this
    {
      svm_probs_internal_[best_quality].second++;
    }
    else if (feature_class == "unknown")
    {
      svm_probs_external_.insert(best_quality);
      if (best_quality >= quality_cutoff) // winning prediction for this assay
      {
        best_feature.setPrimaryID(target_map_.at(target_id).molecule);
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
      String previous_id;
      pair<String, Int> target_id_pair;
      for (FeatureMap::Iterator it = features.begin(); it != features.end(); ++it)
      {
        // features from same assay (same "PeptideRef") appear consecutively;
        // if this is a new assay, finalize the previous one:
        target_id_pair = extractTargetID_(*it, true);
        // target ID incl. charge:
        String target_id = target_id_pair.first + "/" + String(target_id_pair.second);
        if (target_id != previous_id)
        {
          if (!previous_id.empty())
          {
            filterFeaturesFinalizeAssay_(*best_it, best_quality, svm_quality_cutoff,
                                         target_id_pair.first);
            best_quality = 0.0;
          }
          previous_id = target_id;
        }

        // find feature with best quality (from SVM prediction):
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
      }
      // finalize set of features from the last assay:
      filterFeaturesFinalizeAssay_(*best_it, best_quality, svm_quality_cutoff,
                                   target_id_pair.first);
    }

    // remove "losing" feature candidates (without assigned primary ID):
    features.erase(remove_if(features.begin(), features.end(),
                             [](const Feature& feature)
                             { return !feature.hasPrimaryID(); }),
                   features.end());
  }


  void FeatureFinderIdentificationAlgorithm::calculateFDR_(FeatureMap& features)
  {
    // cumulate the true/false positive counts, in decreasing probability order:
    Size n_false = 0, n_true = 0;
    for (map<double, pair<Size, Size> >::reverse_iterator prob_it =
           svm_probs_internal_.rbegin(); prob_it != svm_probs_internal_.rend();
         ++prob_it)
    {
      n_true += prob_it->second.first;
      n_false += prob_it->second.second;
      prob_it->second.first = n_true;
      prob_it->second.second = n_false;
    }

    // print FDR for features that made the cut-off:
    map<double, pair<Size, Size> >::iterator prob_it =
      svm_probs_internal_.lower_bound(svm_min_prob_);
    if (prob_it != svm_probs_internal_.end())
    {
      float fdr = float(prob_it->second.second) / (prob_it->second.first +
                                                   prob_it->second.second);
      OPENMS_LOG_INFO << "Estimated FDR of features detected based on 'external' IDs: "
                      << fdr * 100.0 << "%" << endl;
      fdr = (fdr * n_external_features_) / (n_external_features_ +
                                            n_internal_features_);
      OPENMS_LOG_INFO << "Estimated FDR of all detected features: "
                      << fdr * 100.0 << "%" << endl;
    }

    // calculate q-values:
    vector<double> qvalues;
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
    vector<double> fdr_probs, fdr_qvalues;
    vector<double>::iterator qv_it = qvalues.begin();
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
    multiset<double>::reverse_iterator ext_it = svm_probs_external_.rbegin();
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
        // find highest FDR prob. that is less-or-equal to the feature prob.:
        vector<double>::iterator pos = upper_bound(fdr_probs.begin(),
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


  void FeatureFinderIdentificationAlgorithm::run(FeatureMap& features,
                                                 const vector<PeptideIdentification>& peptides,
                                                 const vector<ProteinIdentification>& proteins,
                                                 const vector<PeptideIdentification>& peptides_ext,
                                                 const vector<ProteinIdentification>& proteins_ext,
                                                 const FeatureMap& seeds,
                                                 const String& spectra_file) {
    IdentificationData id_data, id_data_ext;
    IdentificationDataConverter::importIDs(id_data, proteins, peptides);
    if (!peptides_ext.empty())
    {
      IdentificationDataConverter::importIDs(id_data_ext, proteins_ext, peptides_ext);
    }
    if (!seeds.empty())
    {
      convertSeeds(seeds, id_data, spectra_file);
    }

    run(features, id_data, id_data_ext, spectra_file);
    IdentificationDataConverter::exportFeatureIDs(features);
  }

}
