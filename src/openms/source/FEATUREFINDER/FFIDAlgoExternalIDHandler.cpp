// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FFIDAlgoExternalIDHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ML/SVM/SimpleSVM.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <algorithm>
#include <random>

namespace OpenMS
{
  FFIDAlgoExternalIDHandler::FFIDAlgoExternalIDHandler() :
    n_external_peptides_(0),
    n_external_features_(0),
    svm_n_parts_(3),
    svm_n_samples_(0),
    svm_min_prob_(0.0),
    n_internal_features_(0)
  {
  }

  void FFIDAlgoExternalIDHandler::initSVMParameters_(const Param& param)
  {
    svm_min_prob_ = param.getValue("svm:min_prob");
    svm_n_parts_ = param.getValue("svm:xval");
    svm_n_samples_ = param.getValue("svm:samples");
    svm_xval_out_ = param.getValue("svm:xval_out").toString();
    svm_quality_cutoff = svm_min_prob_;
    svm_predictor_names_ = ListUtils::create<String>(param.getValue("svm:predictors").toString());
    debug_level_ = param.getValue("debug");
  }

  void FFIDAlgoExternalIDHandler::reset()
  {
    external_peptide_map_.clear();
    rt_transformation_ = TransformationDescription();
    n_external_peptides_ = 0;
    n_external_features_ = 0;
    svm_probs_external_.clear();
    svm_probs_internal_.clear();
    n_internal_features_ = 0;
  }
  
  void FFIDAlgoExternalIDHandler::addExternalPeptide(PeptideIdentification& peptide)
  {
    if (peptide.getHits().empty())
    {
      return;
    }
    
    peptide.sort();
    PeptideHit& hit = peptide.getHits()[0];
    peptide.getHits().resize(1);
    
    Int charge = hit.getCharge();
    double rt = peptide.getRT();
    double mz = peptide.getMZ();
    
    external_peptide_map_[hit.getSequence()][charge].emplace(rt, &peptide);
    
    OPENMS_LOG_DEBUG_NOFILE << "Adding peptide (external) " << hit.getSequence()
                         << "; CHG: " << charge << "; RT: " << rt
                         << "; MZ: " << mz << std::endl;
  }
  
  void FFIDAlgoExternalIDHandler::processExternalPeptides(std::vector<PeptideIdentification>& peptides_ext)
  {
    for (PeptideIdentification& pep : peptides_ext)
    {
      addExternalPeptide(pep);
      pep.setMetaValue("FFId_category", "external");
    }
    
    n_external_peptides_ = external_peptide_map_.size();
  }
  
  double FFIDAlgoExternalIDHandler::alignInternalAndExternalIDs(
      const std::vector<PeptideIdentification>& peptides_internal,
      const std::vector<PeptideIdentification>& peptides_external,
      double rt_quantile)
  {
    // Reset the handler state
    reset();
    
    // Align internal and external IDs to estimate RT shifts:
    MapAlignmentAlgorithmIdentification aligner;
    aligner.setReference(peptides_external); // go from internal to external scale
    std::vector<std::vector<PeptideIdentification>> aligner_peptides(1, peptides_internal);
    std::vector<TransformationDescription> aligner_trafos;

    OPENMS_LOG_INFO << "Realigning internal and external IDs...";
    aligner.align(aligner_peptides, aligner_trafos);
    rt_transformation_ = aligner_trafos[0];
    
    std::vector<double> aligned_diffs;
    rt_transformation_.getDeviations(aligned_diffs);
    
    // Calculate RT uncertainty based on quantile
    Size index = std::max(Size(0), Size(rt_quantile * static_cast<double>(aligned_diffs.size())) - 1);
    double rt_uncertainty = aligned_diffs[index];
    
    try
    {
      aligner_trafos[0].fitModel("lowess");
      rt_transformation_ = aligner_trafos[0];
    }
    catch (Exception::BaseException& e)
    {
      OPENMS_LOG_ERROR << "Error: Failed to align RTs of internal/external peptides. "
                     << "RT information will not be considered in the SVM classification. "
                     << "The original error message was:\n" << e.what() << std::endl;
    }
    
    return rt_uncertainty;
  }
  
  double FFIDAlgoExternalIDHandler::transformRT(double rt) const
  {
    return rt_transformation_.apply(rt);
  }
  
  bool FFIDAlgoExternalIDHandler::hasRTTransformation() const
  {
    return !rt_transformation_.getDataPoints().empty();
  }
  
  const TransformationDescription& FFIDAlgoExternalIDHandler::getRTTransformation() const
  {
    return rt_transformation_;
  }
  
  void FFIDAlgoExternalIDHandler::addExternalPeptideToMap_(PeptideIdentification& peptide,
                             std::map<AASequence,
                             std::map<Int, std::pair<std::multimap<double, PeptideIdentification*>,
                                                    std::multimap<double, PeptideIdentification*>>>>& peptide_map)
  {
    if (peptide.getHits().empty()) return;
    
    peptide.sort();
    PeptideHit& hit = peptide.getHits()[0];
    peptide.getHits().resize(1);
    
    Int charge = hit.getCharge();
    double rt = peptide.getRT();
    
    // Add to the external map (second in the pair)
    peptide_map[hit.getSequence()][charge].second.emplace(rt, &peptide);
  }
  
  bool FFIDAlgoExternalIDHandler::fillExternalRTMap_(const AASequence& sequence, Int charge,
                         std::multimap<double, PeptideIdentification*>& rt_map)
  {
    auto seq_it = external_peptide_map_.find(sequence);
    if (seq_it == external_peptide_map_.end()) return false;
    
    auto charge_it = seq_it->second.find(charge);
    if (charge_it == seq_it->second.end()) return false;
    
    rt_map.insert(charge_it->second.begin(), charge_it->second.end());
    return true;
  }
  
  void FFIDAlgoExternalIDHandler::annotateFeatureWithExternalIDs_(Feature& feature)
  {
    feature.setMetaValue("n_total_ids", 0);
    feature.setMetaValue("n_matching_ids", -1);
    feature.setMetaValue("feature_class", "unknown");
  }
  
  void FFIDAlgoExternalIDHandler::addDummyPeptideID_(Feature& feature, const PeptideIdentification* ext_id)
  {
    if (!ext_id) return;
    
    PeptideIdentification id = *ext_id;
    id.clearMetaInfo();
    id.setMetaValue("FFId_category", "implied");
    id.setRT(feature.getRT());
    id.setMZ(feature.getMZ());
    // Only one peptide hit per ID - see function "addPeptideToMap_":
    PeptideHit& hit = id.getHits()[0];
    hit.clearMetaInfo();
    hit.setScore(0.0);
    feature.getPeptideIdentifications().push_back(id);
  }
  
  void FFIDAlgoExternalIDHandler::handleExternalFeature_(Feature& feature, double prob_positive, double quality_cutoff)
  {
    svm_probs_external_.insert(prob_positive);
    
    if (prob_positive >= quality_cutoff)
    {
      feature.setOverallQuality(prob_positive);
      ++n_external_features_;
    }
  }
  
  void FFIDAlgoExternalIDHandler::adjustFDRForExternalFeatures_(std::vector<double>& fdr_probs,
                                   std::vector<double>& fdr_qvalues,
                                   Size n_internal_features)
  {
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
        (external_count + n_internal_features);
    }
  }

  void FFIDAlgoExternalIDHandler::checkNumObservations_(Size n_pos, Size n_neg, const String& note) const
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
void FFIDAlgoExternalIDHandler::getUnbiasedSample_(const std::multimap<double, std::pair<Size, bool> >& valid_obs,
                          std::map<Size, double>& training_labels)
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
    std::multimap<double, std::pair<Size, bool> >::const_iterator begin, middle, end;
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

  void FFIDAlgoExternalIDHandler::getRandomSample_(std::map<Size, double>& training_labels)
  {
    // Pick a random subset of size "svm_n_samples_" for training: Shuffle the whole
    // sequence, then select the first "svm_n_samples_" elements.
    std::vector<Size> selection;
    selection.reserve(training_labels.size());
    for (auto it = training_labels.begin(); it != training_labels.end(); ++it)
    {
      selection.push_back(it->first);
    }
    Math::RandomShuffler shuffler;
    shuffler.portable_random_shuffle(selection.begin(), selection.end());
    // However, ensure that at least "svm_n_parts_" pos./neg. observations are
    // included (for cross-validation) - there must be enough, otherwise
    // "checkNumObservations" would have thrown an error. To this end, move
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
    for (std::vector<Size>::iterator it = selection.begin(); it != selection.end();
          ++it)
    {
      temp[*it] = training_labels[*it];
    }
    training_labels.swap(temp);
  }

  void FFIDAlgoExternalIDHandler::classifyFeaturesWithSVM(FeatureMap& features, const Param& param)
  {
    // Initialize SVM parameters in the external ID handler
    initSVMParameters_(param);

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
                    << feat.getUniqueId() << "'" << std::endl;
          predictors.erase(pred);
          break;
        }
        predictors[pred].push_back(feat.getMetaValue(pred));
      }
    }

    // get labels for SVM:
    std::map<Size, double> training_labels;
    bool no_selection = param.getValue("svm:no_selection") == "true";
    // mapping (for bias correction): intensity -> (index, positive?)
    std::multimap<double, std::pair<Size, bool> > valid_obs;
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
          valid_obs.insert(std::make_pair(intensity, std::make_pair(feat_index,
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
    svm_params.update(param.copy("svm:", true), false, no_log);
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
      features[i].setMetaValue("predicted_class", predictions[i].outcome);
      double prob_positive = predictions[i].probabilities[1];
      features[i].setMetaValue("predicted_probability", prob_positive);
      // @TODO: store previous (OpenSWATH) overall quality in a meta value?
      features[i].setOverallQuality(prob_positive);
    }
  }

  void FFIDAlgoExternalIDHandler::finalizeAssayFeatures_(Feature& best_feature, double best_quality, double quality_cutoff)
  {
    const String& feature_class = best_feature.getMetaValue("feature_class");
    if (feature_class == "positive") // true positive prediction
    {
      svm_probs_internal_[best_quality].first++;
    }
    else if ((feature_class == "negative") || // false positive prediction
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

  void FFIDAlgoExternalIDHandler::filterClassifiedFeatures(FeatureMap& features, double quality_cutoff)
  {
    if (features.empty())
    {
      return;
    }
    
    // Remove features with class "negative" or "ambiguous", keep "positive".
    // For class "unknown", for every assay (meta value "PeptideRef"), keep
    // the feature with highest "predicted_probability" (= overall quality),
    // subject to the "svm:min_prob" threshold.
    // We mark features for removal by setting their overall quality to zero.
    n_internal_features_ = 0;
    n_external_features_ = 0;
    FeatureMap::Iterator best_it = features.begin();
    double best_quality = 0.0;
    String previous_ref;
    for (FeatureMap::Iterator it = features.begin(); it != features.end(); ++it)
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
          finalizeAssayFeatures_(*best_it, best_quality, quality_cutoff);
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
    finalizeAssayFeatures_(*best_it, best_quality, quality_cutoff);

    features.erase(std::remove_if(features.begin(), features.end(),
                             [](const Feature& f) {
                               return f.getOverallQuality() == 0.0;
                             }),
                   features.end());
  }

  void FFIDAlgoExternalIDHandler::calculateFDR(FeatureMap& features)
  {
    if (getSVMProbsInternal().empty()) return;

    // cumulate the true/false positive counts, in decreasing probability order:
    Size n_false = 0, n_true = 0;
    for (std::map<double, std::pair<Size, Size> >::reverse_iterator prob_it =
           svm_probs_internal_.rbegin(); prob_it != svm_probs_internal_.rend();
         ++prob_it)
    {
      n_true += prob_it->second.first;
      n_false += prob_it->second.second;
      prob_it->second.first = n_true;
      prob_it->second.second = n_false;
    }

    // print FDR for features that made the cut-off:
    std::map<double, std::pair<Size, Size> >::iterator prob_it =
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
    adjustFDRForExternalFeatures_(fdr_probs, fdr_qvalues, n_internal_features_);
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
        std::vector<double>::iterator pos = std::upper_bound(fdr_probs.begin(),
                                                   fdr_probs.end(), prob);
        if (pos != fdr_probs.begin())
        {
          --pos;
        }
        Size dist = std::distance(fdr_probs.begin(), pos);
        feat.setMetaValue("q-value", fdr_qvalues[dist]);
      }
    }
  }
  
  const std::map<double, std::pair<Size, Size> >& FFIDAlgoExternalIDHandler::getSVMProbsInternal() const
  {
    return svm_probs_internal_;
  }
  
} // namespace OpenMS