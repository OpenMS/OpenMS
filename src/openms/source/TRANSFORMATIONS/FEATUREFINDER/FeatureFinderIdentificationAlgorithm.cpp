// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  FeatureFinderIdentificationAlgorithm::FeatureFinderIdentificationAlgorithm() :
    DefaultParamHandler("FeatureFinderIdentificationAlgorithm")
  {
    StringList output_file_tags;
    output_file_tags.push_back("output file");

    defaults_.setValue("candidates_out", "", "Optional output file with feature candidates.", output_file_tags);

    defaults_.setValue("debug", 0, "Debug level for feature detection.", ListUtils::create<String>("advanced"));
    defaults_.setMinInt("debug", 0);

    defaults_.setValue("extract:mz_window", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)");
    defaults_.setMinFloat("extract:mz_window", 0.0);
    defaults_.setValue("extract:n_isotopes", 2, "Number of isotopes to include in each peptide assay.");
    defaults_.setMinInt("extract:n_isotopes", 2);
    defaults_.setValue(
      "extract:isotope_pmin",
      0.0, 
      "Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.",
      ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("extract:isotope_pmin", 0.0);
    defaults_.setMaxFloat("extract:isotope_pmin", 1.0);
    defaults_.setValue(
      "extract:rt_quantile", 
      0.95, 
      "Quantile of the RT deviations between aligned internal and external IDs to use for scaling the RT extraction window",
      ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("extract:rt_quantile", 0.0);
    defaults_.setMaxFloat("extract:rt_quantile", 1.0);

    defaults_.setValue(
      "extract:rt_window", 
      0.0, 
      "RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.",
      ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("extract:rt_window", 0.0);

    defaults_.setSectionDescription("extract", "Parameters for ion chromatogram extraction");

    defaults_.setValue("detect:peak_width", 60.0, "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.");
    defaults_.setMinFloat("detect:peak_width", 0.0);
    defaults_.setValue(
      "detect:min_peak_width", 
      0.2, 
      "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.",
      ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("detect:min_peak_width", 0.0);

    defaults_.setValue(
      "detect:signal_to_noise", 
      0.8, 
      "Signal-to-noise threshold for OpenSWATH feature detection",
       ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("detect:signal_to_noise", 0.1);
    defaults_.setValue("detect:mapping_tolerance", 0.0, "RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.");
    defaults_.setMinFloat("detect:mapping_tolerance", 0.0);

    defaults_.setSectionDescription("detect", "Parameters for detecting features in extracted ion chromatograms");

    // parameters for SVM classification:
    defaults_.setValue("svm:samples", 0, "Number of observations to use for training ('0' for all)");
    defaults_.setMinInt("svm:samples", 0);
    defaults_.setValue("svm:no_selection", "false", "By default, roughly the same number of positive and negative observations, with the same intensity distribution, are selected for training. This aims to reduce biases, but also reduces the amount of training data. Set this flag to skip this procedure and consider all available observations (subject to 'svm:samples').");
    defaults_.setValidStrings("svm:no_selection", ListUtils::create<String>("true,false"));
    defaults_.setValue("svm:xval_out", "", "Output file: SVM cross-validation (parameter optimization) results", output_file_tags);
    defaults_.setValidStrings("svm:xval_out", ListUtils::create<String>("csv"));
    defaults_.insert("svm:", SimpleSVM().getParameters());

    // available scores: initialPeakQuality,total_xic,peak_apices_sum,var_xcorr_coelution,var_xcorr_coelution_weighted,var_xcorr_shape,var_xcorr_shape_weighted,var_library_corr,var_library_rmsd,var_library_sangle,var_library_rootmeansquare,var_library_manhattan,var_library_dotprod,var_intensity_score,nr_peaks,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,var_massdev_score_weighted,var_bseries_score,var_yseries_score,var_dotprod_score,var_manhatt_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score
    // exclude some redundant/uninformative scores:
    // @TODO: intensity bias introduced by "peak_apices_sum"?
    // names of scores to use as SVM features
    String score_metavalues = "peak_apices_sum,var_xcorr_coelution,var_xcorr_shape,var_library_sangle,var_intensity_score,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,main_var_xx_swath_prelim_score";

    defaults_.setValue(
      "svm:predictors", 
      score_metavalues, 
      "Names of OpenSWATH scores to use as predictors for the SVM (comma-separated list)",
      ListUtils::create<String>("advanced"));

    defaults_.setValue(
      "svm:min_prob", 
      0.0, 
      "Minimum probability of correctness, as predicted by the SVM, required to retain a feature candidate",
      ListUtils::create<String>("advanced"));
    defaults_.setMinFloat("svm:min_prob", 0.0);
    defaults_.setMaxFloat("svm:min_prob", 1.0);

    defaults_.setSectionDescription("svm", "Parameters for scoring features using a support vector machine (SVM)");

    // parameters for model fitting (via ElutionModelFitter):
    StringList models = ListUtils::create<String>("symmetric,asymmetric,none");
    defaults_.setValue("model:type", models[0], "Type of elution model to fit to features");
    defaults_.setValidStrings("model:type", models);
    defaults_.insert("model:", ElutionModelFitter().getParameters()); // copy parameters
    defaults_.remove("model:asymmetric");

    defaults_.setSectionDescription("model", "Parameters for fitting elution models to features");

    defaultsToParam_();
  }

  void FeatureFinderIdentificationAlgorithm::run(
    vector<PeptideIdentification> peptides,
    const vector<ProteinIdentification>& proteins,
    vector<PeptideIdentification> peptides_ext,
    vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds
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

    // initialize algorithm classes needed later:
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    if ((elution_model_ != "none") || (!candidates_out_.empty()))
    {
      params.setValue("write_convex_hull", "true");
    }
    if (min_peak_width_ < 1.0) min_peak_width_ *= peak_width_;
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
    feat_finder_.setParameters(params);
    feat_finder_.setLogType(ProgressLogger::NONE);
    feat_finder_.setStrictFlag(false);

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
      int index = max(0, int(rt_quantile_ * aligned_diffs.size()) - 1);
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
      if (map_tol < 1.0) map_tol *= (2 * peak_width_); // relative tolerance
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
    peptides.reserve(peptides.size() + seeds.size());

    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      addPeptideToMap_(*pep_it, peptide_map_);
      pep_it->setMetaValue("FFId_category", "internal");
    }

    // TODO make sure that only assembled traces (more than one trace -> has a charge)
    // see FeatureFindingMetabo: defaults_.setValue("remove_single_traces", "false", "Remove unassembled traces (single traces).");
    Size seeds_added(0);
    for (FeatureMap::ConstIterator f_it = seeds.begin(); f_it != seeds.end(); ++f_it)
    {
      // check if already a peptide in peptide_map_ that is close in RT and MZ
      // if so don't add seed
      bool peptide_already_exists = false;
      for (const auto & peptide : peptides)
      {
        double seed_RT = static_cast<double>(f_it->getRT());
        double seed_MZ = static_cast<double>(f_it->getMZ());
	double seed_charge = f_it->getCharge();
        double peptide_RT = peptide.getRT();
        double peptide_MZ = peptide.getMZ();

        // RT or MZ values of seed match in range -> peptide already exists -> don't add seed
        // Consider up to 5th isotopic trace (e.g., because of seed misassignment)
        double th_tolerance = mz_window_ppm_ ? mz_window_ * 1e-6 * peptide_MZ : mz_window_;
        if ((fabs(seed_RT - peptide_RT) <= rt_window_) &&
           ((fabs(seed_MZ - peptide_MZ) <= th_tolerance) ||
             fabs(seed_MZ - (1.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (2.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (3.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (4.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance ||
             fabs(seed_MZ - (5.0/seed_charge) * Constants::C13C12_MASSDIFF_U - peptide_MZ) <= th_tolerance)
            )
        {
          peptide_already_exists = true;
          break;
        }
      }

      if (!peptide_already_exists)
      {
        peptides.emplace_back();
        PeptideHit seed_hit;
        seed_hit.setCharge(f_it->getCharge());

        const String pseudo_mod_name = String(100000 + seeds_added);

        // Check if pseudo mod is already there.
        // Multiple runs of the algorithm might have already registered it
        if (!ModificationsDB::getInstance()->has("[" + pseudo_mod_name + "]"))
        {
          ResidueModification * new_mod = new ResidueModification();
          new_mod->setFullId("[" + pseudo_mod_name + "]"); // setting FullId but not Id makes it a user-defined mod
          new_mod->setTermSpecificity(ResidueModification::ANYWHERE);
          new_mod->setUniModRecordId(100000 + seeds_added); // required for TargetedExperimentHelper
          new_mod->setOrigin('X');
          ModificationsDB::getInstance()->addModification(new_mod);
        }

        AASequence some_seq = AASequence::fromString("XXX");
        some_seq.setModification(1, "[" + pseudo_mod_name + "]");
        seed_hit.setSequence(some_seq);
        vector<PeptideHit> seed_hits;
        seed_hits.push_back(seed_hit);
        peptides.back().setHits(seed_hits);
        peptides.back().setRT(f_it->getRT());
        peptides.back().setMZ(f_it->getMZ());
        peptides.back().setMetaValue("FFId_category", "internal");
        addPeptideToMap_(peptides.back(), peptide_map_);
        ++seeds_added;
      }
    }
   OPENMS_LOG_INFO << "Seeds added: " << seeds_added << endl;


    n_internal_peps_ = peptide_map_.size();
    for (vector<PeptideIdentification>::iterator pep_it =
           peptides_ext.begin(); pep_it != peptides_ext.end(); ++pep_it)
    {
      addPeptideToMap_(*pep_it, peptide_map_, true);
      pep_it->setMetaValue("FFId_category", "external");
    }
    n_external_peps_ = peptide_map_.size() - n_internal_peps_;

    OPENMS_LOG_INFO << "Creating assay library..." << endl;
    PeptideRefRTMap ref_rt_map;
    createAssayLibrary_(peptide_map_, ref_rt_map);

    if (debug_level_ >= 666)
    {
      cout << "Writing debug.traml file." << endl; 
      TraMLFile().store("debug.traml", library_);
    }

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Extracting chromatograms..." << endl;
    ChromatogramExtractor extractor;
    // extractor.setLogType(ProgressLogger::NONE);
    {
      vector<OpenSwath::ChromatogramPtr> chrom_temp;
      vector<ChromatogramExtractor::ExtractionCoordinates> coords;
      extractor.prepare_coordinates(chrom_temp, coords, library_,
          numeric_limits<double>::quiet_NaN(), false);

      boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
      OpenSwath::SpectrumAccessPtr spec_temp =
        SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
      extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
          mz_window_ppm_, "tophat");
      extractor.return_chromatogram(chrom_temp, coords, library_, (*shared)[0],
          chrom_data_.getChromatograms(), false);
    }

    OPENMS_LOG_DEBUG << "Extracted " << chrom_data_.getNrChromatograms()
              << " chromatogram(s)." << endl;

    OPENMS_LOG_INFO << "Detecting chromatographic peaks..." << endl;
    // suppress status output from OpenSWATH, unless in debug mode:
    if (debug_level_ < 1) OpenMS_Log_info.remove(cout);
    feat_finder_.pickExperiment(chrom_data_, features, library_,
                                TransformationDescription(), ms_data_);
    if (debug_level_ < 1) OpenMS_Log_info.insert(cout); // revert logging change
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

    features.ensureUniqueId();
  }

  void FeatureFinderIdentificationAlgorithm::postProcess_(
   FeatureMap & features,
   bool with_external_ids)
  {
    // don't do SVM stuff unless we have external data to apply the model to:
    if (with_external_ids) classifyFeatures_(features);

    // store feature candidates before filtering
    if (!candidates_out_.empty())
    {
      FeatureXMLFile().store(candidates_out_, features);
    }

    filterFeatures_(features, with_external_ids);
    OPENMS_LOG_INFO << features.size() << " features left after filtering." << endl;

    if (!svm_probs_internal_.empty()) calculateFDR_(features);

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
      for (FeatureMap::Iterator feat_it = features.begin(); 
           feat_it != features.end(); ++feat_it)
      {
        for (vector<Feature>::iterator sub_it =
               feat_it->getSubordinates().begin(); sub_it != 
               feat_it->getSubordinates().end(); ++sub_it)
        {
          sub_it->getConvexHulls().clear();
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
    for (vector<PeptideIdentification>::iterator pep_it =
           features.getUnassignedPeptideIdentifications().begin(); pep_it !=
           features.getUnassignedPeptideIdentifications().end(); ++pep_it)
    {
      const AASequence& seq = pep_it->getHits()[0].getSequence();
      if (pep_it->getMetaValue("FFId_category") == "internal")
      {
        internal_seqs.insert(seq);
      }
      peptide_map_[seq];
    }
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      if (feat_it->getPeptideIdentifications().empty()) continue;
      const PeptideIdentification& pep_id =
        feat_it->getPeptideIdentifications()[0];
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

  void FeatureFinderIdentificationAlgorithm::createAssayLibrary_(PeptideMap& peptide_map, PeptideRefRTMap& ref_rt_map)
  {
    std::set<String> protein_accessions;

    for (PeptideMap::iterator pm_it = peptide_map.begin();
         pm_it != peptide_map.end(); ++pm_it)
    {
      TargetedExperiment::Peptide peptide;

      const AASequence &seq = pm_it->first;
      OPENMS_LOG_DEBUG << "\nPeptide: " << seq.toString() << std::endl;
      peptide.sequence = seq.toString();
      // @NOTE: Technically, "TargetedExperiment::Peptide" stores the unmodified
      // sequence and the modifications separately. Unfortunately, creating the
      // modifications vector is complex and there is currently no convenient
      // conversion function (see "TargetedExperimentHelper::getAASequence" for
      // the reverse conversion). However, "Peptide" is later converted to
      // "OpenSwath::LightPeptide" anyway, and this is done via "AASequence"
      // (see "OpenSwathDataAccessHelper::convertTargetedPeptide"). So for our
      // purposes it works to just store the sequence including modifications in
      // "Peptide".

      // keep track of protein accessions:
      set<String> current_accessions;
      // internal/external pair
      const pair<RTMap, RTMap> &pair = pm_it->second.begin()->second;
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
      getRTRegions_(pm_it->second, rt_regions);
      OPENMS_LOG_DEBUG << "Found " << rt_regions.size() << " RT region(s)." << std::endl;

      // go through different charge states:
      for (ChargeMap::const_iterator cm_it = pm_it->second.begin();
           cm_it != pm_it->second.end(); ++cm_it)
      {
        Int charge = cm_it->first;

        if (seq.toUnmodifiedString().hasPrefix("XXX")) // seed
        {
          //cout << peptide.sequence << " " << charge << endl;

          String peptide_id = peptide.sequence + "/" + String(charge);
          peptide.setChargeState(charge);
          peptide.id = peptide_id;
          peptide.setPeptideGroupLabel(peptide_id);
          peptide.rts.clear();

          Size counter = 0;
          // accumulate IDs over multiple regions: potentially not needed for seeds
          RTMap &internal_ids = ref_rt_map[peptide_id].first;
          RTMap &external_ids = ref_rt_map[peptide_id].second;
          for (vector<RTRegion>::iterator reg_it = rt_regions.begin();
               reg_it != rt_regions.end(); ++reg_it)
          {
            if (reg_it->ids.count(charge))
            {
              OPENMS_LOG_DEBUG << "Region " << counter + 1 << " (RT: "
                               << float(reg_it->start) << "-" << float(reg_it->end)
                               << ", size " << float(reg_it->end - reg_it->start) << ")"
                               << std::endl;

              peptide.id = peptide_id;
              if (rt_regions.size() > 1)
                peptide.id += ":" + String(++counter);

              auto &a = reg_it->ids[charge].first;
              double mz = a.begin()->second->getMZ();
              // get isotope distribution for peptide:
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

              OPENMS_LOG_DEBUG << "Seed Charge: " << charge << " (m/z: " << mz << ")" << std::endl;

              // store beginning and end of RT region:
              peptide.rts.clear();
              addPeptideRT_(peptide, reg_it->start);
              addPeptideRT_(peptide, reg_it->end);
              library_.addPeptide(peptide);
              generateTransitions_(peptide.id, mz, charge, iso_dist);
            }
            internal_ids.insert(reg_it->ids[charge].first.begin(),
                                reg_it->ids[charge].first.end());
            external_ids.insert(reg_it->ids[charge].second.begin(), // Note: empty
                                reg_it->ids[charge].second.end());
          }
        }
        else
        {
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

          double mz = seq.getMonoWeight(Residue::Full, charge) / charge;
          OPENMS_LOG_DEBUG << "Charge: " << charge << " (m/z: " << mz << ")" << std::endl;
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
          for (vector<RTRegion>::iterator reg_it = rt_regions.begin();
               reg_it != rt_regions.end(); ++reg_it)
          {
            if (reg_it->ids.count(charge))
            {
              OPENMS_LOG_DEBUG << "Region " << counter + 1 << " (RT: "
                               << float(reg_it->start) << "-" << float(reg_it->end)
                               << ", size " << float(reg_it->end - reg_it->start) << ")"
                               << std::endl;

              peptide.id = peptide_id;
              if (rt_regions.size() > 1)
                peptide.id += ":" + String(++counter);

              // store beginning and end of RT region:
              peptide.rts.clear();
              addPeptideRT_(peptide, reg_it->start);
              addPeptideRT_(peptide, reg_it->end);
              library_.addPeptide(peptide);
              generateTransitions_(peptide.id, mz, charge, iso_dist);
            }
            internal_ids.insert(reg_it->ids[charge].first.begin(),
                                reg_it->ids[charge].first.end());
            external_ids.insert(reg_it->ids[charge].second.begin(),
                                reg_it->ids[charge].second.end());
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
    std::vector<RTRegion>& rt_regions) const
  {
    // use RTs from all charge states here to get a more complete picture:
    std::vector<double> rts;
    for (ChargeMap::iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      // "internal" IDs:
      for (RTMap::iterator rt_it = cm_it->second.first.begin();
           rt_it != cm_it->second.first.end(); ++rt_it)
      {
        rts.push_back(rt_it->first);
      }
      // "external" IDs:
      for (RTMap::iterator rt_it = cm_it->second.second.begin();
           rt_it != cm_it->second.second.end(); ++rt_it)
      {
        rts.push_back(rt_it->first);
      }
    }
    sort(rts.begin(), rts.end());
    double rt_tolerance = rt_window_ / 2.0;

    for (vector<double>::iterator rt_it = rts.begin(); rt_it != rts.end();
         ++rt_it)
    {
      // create a new region?
      if (rt_regions.empty() || (rt_regions.back().end < *rt_it - rt_tolerance))
      {
        RTRegion region;
        region.start = *rt_it - rt_tolerance;
        // TODO
        // cppcheck-suppress uninitStructMember
        rt_regions.push_back(region);
      }
      rt_regions.back().end = *rt_it + rt_tolerance;
    }

    // sort the peptide IDs into the regions:
    for (ChargeMap::iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      // regions are sorted by RT, as are IDs, so just iterate linearly:
      std::vector<RTRegion>::iterator reg_it = rt_regions.begin();
      // "internal" IDs:
      for (RTMap::iterator rt_it = cm_it->second.first.begin();
           rt_it != cm_it->second.first.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[cm_it->first].first.insert(*rt_it);
      }
      reg_it = rt_regions.begin(); // reset to start
      // "external" IDs:
      for (RTMap::iterator rt_it = cm_it->second.second.begin();
           rt_it != cm_it->second.second.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[cm_it->first].second.insert(*rt_it);
      }
      // ID references no longer needed (now stored in the RT regions):
      cm_it->second.first.clear();
      cm_it->second.second.clear();
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
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin();
         iso_it != iso_dist.end(); ++iso_it, ++counter)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter + 1);
      String transition_name = peptide_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U * 
                              float(counter) / charge);
      transition.setLibraryIntensity(iso_it->getIntensity());
      transition.setMetaValue("annotation", annotation);
      transition.setPeptideRef(peptide_id);
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso_it->getIntensity();
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
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++i)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      ensureConvexHulls_(*feat_it);
      // remove "fake" IDs generated by OpenSWATH (they would be removed with
      // a warning when writing output, because of missing protein
      // identification with corresponding identifier):
      feat_it->getPeptideIdentifications().clear();
      // annotate subordinates with theoretical isotope intensities:
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        String native_id = sub_it->getMetaValue("native_id");
        sub_it->setMetaValue("isotope_probability", isotope_probs_[native_id]);
      }

      peptide_ref = feat_it->getMetaValue("PeptideRef");
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
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RT internal and external are both empty.");
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
        feat_it->setMetaValue("n_total_ids", rt_internal.size());
        feat_it->setMetaValue("n_matching_ids", id_count);
        if (id_count > 0) // matching IDs -> feature may be correct
        {
          feat_it->setMetaValue("feature_class", "ambiguous");
        }
        else // no matching IDs -> feature is wrong
        {
          feat_it->setMetaValue("feature_class", "negative");
        }
      }
      else // only external IDs -> no validation possible
      {
        feat_it->setMetaValue("n_total_ids", 0);
        feat_it->setMetaValue("n_matching_ids", -1);
        feat_it->setMetaValue("feature_class", "unknown");
        // add "dummy" peptide identification:
        PeptideIdentification id = *(rt_external.begin()->second);
        id.clearMetaInfo();
        id.setMetaValue("FFId_category", "implied");
        id.setRT(feat_it->getRT());
        id.setMZ(feat_it->getMZ());
        // only one peptide hit per ID - see function "addPeptideToMap_":
        PeptideHit& hit = id.getHits()[0];
        hit.clearMetaInfo();
        hit.setScore(0.0);
        feat_it->getPeptideIdentifications().push_back(id);
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

        double rt_min = feat_it->getMetaValue("leftWidth");
        double rt_max = feat_it->getMetaValue("rightWidth");
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
          feat_it->setMetaValue("rt_delta", 0.0);
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
          feat_it->setMetaValue("rt_delta", min(rt_delta1, rt_delta2));
        }
      }
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

  void FeatureFinderIdentificationAlgorithm::ensureConvexHulls_(Feature& feature)
  {
    if (feature.getConvexHulls().empty()) // add hulls for mass traces
    {
      double rt_min = feature.getMetaValue("leftWidth");
      double rt_max = feature.getMetaValue("rightWidth");
      for (vector<Feature>::iterator sub_it = feature.getSubordinates().begin();
           sub_it != feature.getSubordinates().end(); ++sub_it)
      {
        double abs_mz_tol = mz_window_ / 2.0;
        if (mz_window_ppm_)
        {
          abs_mz_tol = sub_it->getMZ() * abs_mz_tol * 1.0e-6;
        }
        ConvexHull2D hull;
        hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() + abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() + abs_mz_tol));
        feature.getConvexHulls().push_back(hull);
      }
    }
  }

  void FeatureFinderIdentificationAlgorithm::addPeptideToMap_(PeptideIdentification& peptide, PeptideMap& peptide_map, bool external) const
  {
    if (peptide.getHits().empty()) return;
    peptide.sort();
    PeptideHit& hit = peptide.getHits()[0];
    if (hit.metaValueExists("target_decoy") && hit.getMetaValue("target_decoy") == "decoy") { return; }
    peptide.getHits().resize(1);
    Int charge = hit.getCharge();
    double rt = peptide.getRT();
    RTMap::value_type pair = make_pair(rt, &peptide);
    if (!external)
    {
      OPENMS_LOG_DEBUG << "Adding " << hit.getSequence() << " " << charge << " " << rt << endl;
      peptide_map[hit.getSequence()][charge].first.insert(pair);
    }
    else
    {
      peptide_map[hit.getSequence()][charge].second.insert(pair);
    }
  }

  void FeatureFinderIdentificationAlgorithm::updateMembers_()
  {
    peak_width_ = param_.getValue("detect:peak_width");
    min_peak_width_ = param_.getValue("detect:min_peak_width");
    signal_to_noise_ = param_.getValue("detect:signal_to_noise");

    rt_quantile_ = param_.getValue("extract:rt_quantile");
    rt_window_ = param_.getValue("extract:rt_window");
    mz_window_ = param_.getValue("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;

    isotope_pmin_ = param_.getValue("extract:isotope_pmin");
    n_isotopes_ = param_.getValue("extract:n_isotopes");

    mapping_tolerance_ = param_.getValue("detect:mapping_tolerance");

    elution_model_ = param_.getValue("model:type");
    // SVM related parameters
    svm_min_prob_ = param_.getValue("svm:min_prob");
    svm_predictor_names_ = ListUtils::create<String>(param_.getValue("svm:predictors").toString());
    svm_xval_out_ = param_.getValue("svm:xval_out");
    svm_quality_cutoff = param_.getValue("svm:min_prob");
    svm_n_parts_ = param_.getValue("svm:xval");
    svm_n_samples_ = param_.getValue("svm:samples");

    // debug
    debug_level_ = param_.getValue("debug");
    candidates_out_ = param_.getValue("candidates_out");
  }

  void FeatureFinderIdentificationAlgorithm::getUnbiasedSample_(const multimap<double, pair<Size, bool> >& valid_obs,
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


  void FeatureFinderIdentificationAlgorithm::getRandomSample_(std::map<Size, Int>& training_labels)
  {
    // @TODO: can this be done with less copying back and forth of data?
    // Pick a random subset of size "svm_n_samples_" for training: Shuffle the whole
    // sequence, then select the first "svm_n_samples_" elements.
    std::vector<Size> selection;
    selection.reserve(training_labels.size());
    for (std::map<Size, Int>::iterator it = training_labels.begin();
         it != training_labels.end(); ++it)
    {
      selection.push_back(it->first);
    }
    random_shuffle(selection.begin(), selection.end());
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
        if (n_obs[label] == svm_n_parts_) break;
      }
    }
    selection.resize(svm_n_samples_);
    // copy the selected subset back:
    std::map<Size, Int> temp;
    for (vector<Size>::iterator it = selection.begin(); it != selection.end();
         ++it)
    {
      temp[*it] = training_labels[*it];
    }
    training_labels.swap(temp);
  }

  void FeatureFinderIdentificationAlgorithm::classifyFeatures_(FeatureMap& features)
  {
    if (features.empty()) return;

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
    for (vector<String>::iterator pred_it = svm_predictor_names_.begin();
         pred_it != svm_predictor_names_.end(); ++pred_it)
    {
      predictors[*pred_it].reserve(features.size());
      for (FeatureMap::Iterator feat_it = features.begin(); 
           feat_it < features.end(); ++feat_it)
      {
        if (!feat_it->metaValueExists(*pred_it))
        {
          OPENMS_LOG_ERROR << "Meta value '" << *pred_it << "' missing for feature '"
                    << feat_it->getUniqueId() << "'" << std::endl;
          predictors.erase(*pred_it);
          break;
        }
        predictors[*pred_it].push_back(feat_it->getMetaValue(*pred_it));
      }
    }

    // get labels for SVM:
    std::map<Size, Int> training_labels;
    bool no_selection = param_.getValue("svm:no_selection") == "true" ? true : false;
    // mapping (for bias correction): intensity -> (index, positive?)
    std::multimap<double, pair<Size, bool> > valid_obs;
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      Int label = -1;
      if (feature_class == "positive") label = 1;
      else if (feature_class == "negative") label = 0;

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

    if (!no_selection) getUnbiasedSample_(valid_obs, training_labels);

    if (svm_n_samples_ > 0) // limited number of samples for training
    {
      if (training_labels.size() < svm_n_samples_)
      {
        OPENMS_LOG_WARN << "Warning: There are only " << training_labels.size()
                 << " valid observations for training." << std::endl;
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
    if (!svm_xval_out_.empty()) svm.writeXvalResults(svm_xval_out_);
    if ((debug_level_ > 0) && String(svm_params.getValue("kernel")) == "linear")
    {
      std::map<String, double> feature_weights;
      svm.getFeatureWeights(feature_weights);
      OPENMS_LOG_DEBUG << "SVM feature weights:" << std::endl;
      for (std::map<String, double>::iterator it = feature_weights.begin();
           it != feature_weights.end(); ++it)
      {
        OPENMS_LOG_DEBUG << "- " << it->first << ": " << it->second << std::endl;
      }
    }

    std::vector<SimpleSVM::Prediction> predictions;
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
    if (features.empty()) return;

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
               << fdr * 100.0 << "%" << std::endl;
      fdr = (fdr * n_external_features_) / (n_external_features_ + 
                                            n_internal_features_);
      OPENMS_LOG_INFO << "Estimated FDR of all detected features: " << fdr * 100.0
               << "%" << std::endl;
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
      if (fdr < min_fdr) min_fdr = fdr;
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
    for (FeatureMap::iterator feat_it = features.begin(); 
         feat_it != features.end(); ++feat_it)
    {
      if (feat_it->getMetaValue("feature_class") == "positive")
      {
        feat_it->setMetaValue("q-value", 0.0);
      }
      else
      {
        double prob = feat_it->getOverallQuality();
        // find highest FDR prob. that is less-or-equal to the feature prob.:
        std::vector<double>::iterator pos = upper_bound(fdr_probs.begin(),
                                                   fdr_probs.end(), prob);
        if (pos != fdr_probs.begin()) --pos;
        Size dist = distance(fdr_probs.begin(), pos);
        feat_it->setMetaValue("q-value", fdr_qvalues[dist]);
      }
    }
  }
}
