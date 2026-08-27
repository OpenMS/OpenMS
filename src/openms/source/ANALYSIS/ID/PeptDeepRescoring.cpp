// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>

namespace
{
  using namespace OpenMS;

  /// Slot layout of the PeptDeep MS2 output: eight ion channels per fragment position,
  /// of which the first four are the singly and doubly charged b and y ions.
  constexpr int ION_B1 = 0, ION_B2 = 1, ION_Y1 = 2, ION_Y2 = 3;
  constexpr Size DEFAULT_ION_CHANNELS = 8;

  /// (ion type, charge, ordinal) -> intensity
  using IonMap = std::map<std::tuple<char, int, int>, double>;

  /// Observed fragment intensities from a hit's peak annotations.
  IonMap observedIons_(const PeptideHit& hit)
  {
    IonMap out;
    for (const PeptideHit::PeakAnnotation& pa : hit.getPeakAnnotations())
    {
      const std::string& ann = pa.annotation;
      // Annotations look like "y5+", "b12++"; take the leading ion letter and its ordinal.
      Size i = 0;
      while (i < ann.size() && !(ann[i] >= 'a' && ann[i] <= 'z')) { ++i; }
      if (i >= ann.size()) { continue; }
      const char type = ann[i];
      if (type != 'a' && type != 'b' && type != 'c' && type != 'x' && type != 'y' && type != 'z') { continue; }
      Size j = i + 1;
      int ordinal = 0;
      while (j < ann.size() && ann[j] >= '0' && ann[j] <= '9') { ordinal = ordinal * 10 + (ann[j] - '0'); ++j; }
      if (ordinal == 0) { continue; }
      const int z = pa.charge > 0 ? pa.charge : 1;
      out[std::make_tuple(type, z, ordinal)] += pa.intensity;
    }
    return out;
  }

  /// Predicted b/y intensities, keyed like observedIons_().
  ///
  /// Fragment position i carries b_(i+1) and y_(len-1-i): the model emits fragments
  /// N-to-C, so the y series is numbered from the other end.
  IonMap predictedIons_(const std::vector<float>& flat, Size pep_len, Size channels)
  {
    IonMap out;
    if (pep_len < 2 || channels == 0) { return out; }
    for (Size pos = 0; pos + 1 < pep_len; ++pos)
    {
      const Size base = pos * channels;
      if (base + ION_Y2 >= flat.size()) { break; }
      out[std::make_tuple('b', 1, static_cast<int>(pos + 1))] = flat[base + ION_B1];
      out[std::make_tuple('b', 2, static_cast<int>(pos + 1))] = flat[base + ION_B2];
      out[std::make_tuple('y', 1, static_cast<int>(pep_len - 1 - pos))] = flat[base + ION_Y1];
      out[std::make_tuple('y', 2, static_cast<int>(pep_len - 1 - pos))] = flat[base + ION_Y2];
    }
    return out;
  }

  struct Similarity
  {
    double cosine = 0.0;
    double spectral_angle = 0.0;
    double pearson = 0.0;
    double frac_found = 0.0;
  };

  /// Agreement between predicted and observed intensities over the predicted ion set.
  /// Predicted ions that were never matched enter with an observed intensity of 0.
  Similarity similarity_(const IonMap& pred, const IonMap& obs, int min_ordinal, Size pep_len,
                         double strong_fraction)
  {
    std::vector<double> p, o;
    p.reserve(pred.size()); o.reserve(pred.size());
    for (const auto& [key, pv] : pred)
    {
      const int ordinal = std::get<2>(key);
      if (ordinal < min_ordinal || ordinal > static_cast<int>(pep_len) - 1) { continue; }
      p.push_back(pv);
      const auto it = obs.find(key);
      o.push_back(it == obs.end() ? 0.0 : it->second);
    }
    Similarity s;
    if (p.empty()) { return s; }

    const double sum_p = std::accumulate(p.begin(), p.end(), 0.0);
    const double sum_o = std::accumulate(o.begin(), o.end(), 0.0);
    if (sum_p <= 0.0 || sum_o <= 0.0) { return s; }

    double dot = 0.0, np2 = 0.0, no2 = 0.0;
    for (Size i = 0; i < p.size(); ++i) { dot += p[i] * o[i]; np2 += p[i] * p[i]; no2 += o[i] * o[i]; }
    if (np2 <= 0.0 || no2 <= 0.0) { return s; }
    s.cosine = std::clamp(dot / (std::sqrt(np2) * std::sqrt(no2)), -1.0, 1.0);
    s.spectral_angle = 1.0 - 2.0 * std::acos(s.cosine) / Constants::PI;

    const double mp = sum_p / p.size(), mo = sum_o / o.size();
    double cov = 0.0, vp = 0.0, vo = 0.0;
    for (Size i = 0; i < p.size(); ++i)
    {
      cov += (p[i] - mp) * (o[i] - mo);
      vp += (p[i] - mp) * (p[i] - mp);
      vo += (o[i] - mo) * (o[i] - mo);
    }
    s.pearson = (vp > 0.0 && vo > 0.0) ? std::clamp(cov / std::sqrt(vp * vo), -1.0, 1.0) : 0.0;

    const double peak = *std::max_element(p.begin(), p.end());
    Size strong = 0, strong_found = 0;
    for (Size i = 0; i < p.size(); ++i)
    {
      if (p[i] > strong_fraction * peak) { ++strong; if (o[i] > 0.0) { ++strong_found; } }
    }
    s.frac_found = strong > 0 ? static_cast<double>(strong_found) / static_cast<double>(strong) : 0.0;
    return s;
  }

  /// peptdeep's categorical instrument encoding.
  int64_t instrumentIndex_(const std::string& name)
  {
    if (name == "Lumos") { return 0; }
    if (name == "QE") { return 1; }
    if (name == "timsTOF") { return 2; }
    if (name == "SciexTOF") { return 3; }
    return 1;
  }
}

namespace OpenMS
{
  PeptDeepRescoring::PeptDeepRescoring() :
    DefaultParamHandler("PeptDeepRescoring")
  {
    defaults_.setValue("ms2_model", "", "Path to the PeptDeep MS2 fragment-intensity ONNX model. Required.");
    defaults_.setValue("rt_model", "", "Path to the PeptDeep retention-time ONNX model. Required.");

    defaults_.setValue("instrument", "QE", "Instrument class passed to the MS2 model.");
    defaults_.setValidStrings("instrument", {"Lumos", "QE", "timsTOF", "SciexTOF"});

    defaults_.setValue("nce", -1.0, "Normalised collision energy for the MS2 model. Negative selects it automatically: a grid centred on the collision energy recorded in the spectra is scored on a sample of confident PSMs and the best is kept. The instrument's own value is only a starting point -- the model's NCE scale does not coincide with it.");
    defaults_.setValue("nce_grid_halfwidth", 6.0, "Half-width of the automatic NCE search grid around its centre.", {"advanced"});
    defaults_.setMinFloat("nce_grid_halfwidth", 0.0);
    defaults_.setValue("nce_grid_step", 2.0, "Step size of the automatic NCE search grid.", {"advanced"});
    defaults_.setMinFloat("nce_grid_step", 0.1);
    defaults_.setValue("nce_fallback", 30.0, "NCE grid centre used when the spectra do not record a collision energy.", {"advanced"});
    defaults_.setValue("nce_sample_size", 2000, "Number of confident PSMs scored per candidate NCE.", {"advanced"});
    defaults_.setMinInt("nce_sample_size", 50);

    defaults_.setValue("calibration_quantile", 0.5, "Fraction of PSMs, ranked by search score, held out of the calibration sets. 0.5 fits on the better-scoring half. Selection uses the search score rather than spectral similarity, so the RT feature stays independent of the MS2 features.", {"advanced"});
    defaults_.setMinFloat("calibration_quantile", 0.0);
    defaults_.setMaxFloat("calibration_quantile", 0.99);

    defaults_.setValue("rt_model_type", "b_spline", "Model mapping predicted onto observed retention time.");
    defaults_.setValidStrings("rt_model_type", {"b_spline", "lowess", "linear"});
    defaults_.setValue("rt_num_nodes", 5, "Number of B-spline nodes for the retention-time calibration. Fewer means smoother; an over-flexible fit tracks the calibration set instead of the gradient.", {"advanced"});
    defaults_.setMinInt("rt_num_nodes", 2);

    defaults_.setValue("ms2_min_ordinal", 1, "Ignore predicted fragments below this ordinal.", {"advanced"});
    defaults_.setMinInt("ms2_min_ordinal", 1);
    defaults_.setValue("ms2_strong_fraction", 0.05, "An ion counts as confidently predicted, for ms2_frac_pred_found, above this fraction of the peptide's most intense predicted ion.", {"advanced"});
    defaults_.setMinFloat("ms2_strong_fraction", 0.0);
    defaults_.setMaxFloat("ms2_strong_fraction", 1.0);

    defaults_.setValue("batch_size", 500, "Peptides per ONNX inference call.", {"advanced"});
    defaults_.setMinInt("batch_size", 1);
    defaults_.setValue("threads", 4, "ONNX intra-op threads.", {"advanced"});
    defaults_.setMinInt("threads", 1);

    defaultsToParam_();
  }

  PeptDeepRescoring::~PeptDeepRescoring() = default;

  void PeptDeepRescoring::updateMembers_()
  {
    ms2_model_ = param_.getValue("ms2_model").toString();
    rt_model_ = param_.getValue("rt_model").toString();
    instrument_ = param_.getValue("instrument").toString();
    nce_ = param_.getValue("nce");
    nce_grid_halfwidth_ = param_.getValue("nce_grid_halfwidth");
    nce_grid_step_ = param_.getValue("nce_grid_step");
    nce_fallback_ = param_.getValue("nce_fallback");
    nce_sample_size_ = static_cast<Size>((int)param_.getValue("nce_sample_size"));
    calibration_quantile_ = param_.getValue("calibration_quantile");
    rt_model_type_ = param_.getValue("rt_model_type").toString();
    rt_num_nodes_ = static_cast<Size>((int)param_.getValue("rt_num_nodes"));
    ms2_min_ordinal_ = (int)param_.getValue("ms2_min_ordinal");
    ms2_strong_fraction_ = param_.getValue("ms2_strong_fraction");
    batch_size_ = static_cast<Size>((int)param_.getValue("batch_size"));
    threads_ = (int)param_.getValue("threads");
  }

  StringList PeptDeepRescoring::featureNames()
  {
    return {Constants::UserParam::MS2_COSINE,
            Constants::UserParam::MS2_SPECTRAL_ANGLE,
            Constants::UserParam::MS2_PEARSON,
            Constants::UserParam::MS2_FRAC_PRED_FOUND,
            Constants::UserParam::RT_ABS_ERROR};
  }

  double PeptDeepRescoring::getUsedNCE() const { return used_nce_; }

  double PeptDeepRescoring::getRTCalibrationError() const { return rt_calibration_error_; }

  void PeptDeepRescoring::annotate(const PeakMap& spectra,
                                   std::vector<ProteinIdentification>& protein_ids,
                                   PeptideIdentificationList& peptide_ids)
  {
    if (ms2_model_.empty() || rt_model_.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "PeptDeepRescoring needs both 'ms2_model' and 'rt_model' to be set.");
    }
    for (const std::string& m : {ms2_model_, rt_model_})
    {
      if (!File::exists(m))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, m);
      }
    }
    if (peptide_ids.empty()) { return; }

    // ---- gather the PSMs, and the distinct (sequence, charge) pairs to predict ----
    struct Row { Size pi, hit; Size len; Size key; double score; };
    std::vector<Row> rows;
    std::map<std::pair<std::string, int>, Size> seen;
    std::vector<std::string> uniq_seq;
    std::vector<float> uniq_charge;
    Size n_with_annotations = 0;

    for (Size pi = 0; pi < peptide_ids.size(); ++pi)
    {
      const std::vector<PeptideHit>& hits = peptide_ids[pi].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const AASequence& seq = hits[h].getSequence();
        if (seq.empty()) { continue; }
        if (!hits[h].getPeakAnnotations().empty()) { ++n_with_annotations; }
        const std::string str = seq.toString();
        int z = hits[h].getCharge();
        if (z <= 0) { z = 2; }
        const auto key = std::make_pair(str, z);
        auto it = seen.find(key);
        if (it == seen.end())
        {
          it = seen.emplace(key, uniq_seq.size()).first;
          uniq_seq.push_back(str);
          uniq_charge.push_back(static_cast<float>(z));
        }
        rows.push_back(Row{pi, h, seq.size(), it->second, hits[h].getScore()});
      }
    }
    if (rows.empty()) { return; }
    if (n_with_annotations == 0)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "No PSM carries peak annotations, so every MS2 feature would be zero. Re-run the "
        "search with fragment annotation enabled (ProSE: annotate:PSM contains 'ALL' or "
        "'fragment_annotation').");
    }

    // Confident PSMs, used both to score NCE candidates and to fit the RT calibration.
    // Ranking by search score keeps this independent of the MS2 features themselves.
    const bool higher_better = peptide_ids[0].isHigherScoreBetter();
    std::vector<Size> by_score(rows.size());
    std::iota(by_score.begin(), by_score.end(), 0);
    std::sort(by_score.begin(), by_score.end(), [&](Size a, Size b)
      { return higher_better ? rows[a].score > rows[b].score : rows[a].score < rows[b].score; });
    const Size n_conf = std::max<Size>(1, static_cast<Size>(rows.size() * (1.0 - calibration_quantile_)));
    std::vector<char> is_confident(rows.size(), 0);
    for (Size i = 0; i < n_conf; ++i) { is_confident[by_score[i]] = 1; }

    PeptDeepMS2Inference ms2(ms2_model_, threads_, batch_size_);

    // ---- normalised collision energy ----
    const int64_t instrument = instrumentIndex_(instrument_);
    if (nce_ > 0.0)
    {
      used_nce_ = nce_;
    }
    else
    {
      // Centre the grid on what the instrument recorded, when it recorded anything.
      std::vector<double> ces;
      for (const MSSpectrum& s : spectra)
      {
        if (s.getMSLevel() != 2 || s.getPrecursors().empty()) { continue; }
        const Precursor& prec = s.getPrecursors()[0];
        if (prec.metaValueExists("collision energy"))
        {
          ces.push_back(static_cast<double>(prec.getMetaValue("collision energy")));
        }
      }
      double centre = nce_fallback_;
      if (!ces.empty())
      {
        std::nth_element(ces.begin(), ces.begin() + ces.size() / 2, ces.end());
        centre = ces[ces.size() / 2];
      }

      std::vector<double> grid;
      for (double x = centre - nce_grid_halfwidth_; x <= centre + nce_grid_halfwidth_ + 1e-9; x += nce_grid_step_)
      {
        if (x > 0.0) { grid.push_back(x); }
      }
      if (grid.empty()) { grid.push_back(centre); }

      // Score each candidate on a sample of confident PSMs and keep the best.
      std::vector<Size> sample;
      for (Size i : by_score)
      {
        if (sample.size() >= nce_sample_size_) { break; }
        sample.push_back(i);
      }
      std::vector<std::string> s_seq;
      std::vector<float> s_charge;
      s_seq.reserve(sample.size()); s_charge.reserve(sample.size());
      for (Size i : sample) { s_seq.push_back(uniq_seq[rows[i].key]); s_charge.push_back(uniq_charge[rows[i].key]); }

      double best_median = -1.0;
      used_nce_ = centre;
      startProgress(0, grid.size(), "Selecting collision energy...");
      Size done = 0;
      for (double cand : grid)
      {
        const std::vector<float> s_nce(s_seq.size(), static_cast<float>(cand));
        const std::vector<int64_t> s_inst(s_seq.size(), instrument);
        const std::vector<std::vector<float>> pred = ms2.predictMS2(s_seq, s_charge, s_nce, s_inst);

        std::vector<double> cosines;
        cosines.reserve(sample.size());
        for (Size k = 0; k < sample.size(); ++k)
        {
          const Row& r = rows[sample[k]];
          const Size channels = r.len > 1 ? pred[k].size() / (r.len - 1) : DEFAULT_ION_CHANNELS;
          const PeptideHit& hit = peptide_ids[r.pi].getHits()[r.hit];
          cosines.push_back(similarity_(predictedIons_(pred[k], r.len, channels ? channels : DEFAULT_ION_CHANNELS),
                                        observedIons_(hit), ms2_min_ordinal_, r.len, ms2_strong_fraction_).cosine);
        }
        std::nth_element(cosines.begin(), cosines.begin() + cosines.size() / 2, cosines.end());
        const double med = cosines[cosines.size() / 2];
        OPENMS_LOG_DEBUG << "[PeptDeepRescoring] NCE " << cand << " -> median cosine " << med << std::endl;
        if (med > best_median) { best_median = med; used_nce_ = cand; }
        setProgress(++done);
      }
      endProgress();
      OPENMS_LOG_INFO << "[PeptDeepRescoring] NCE " << used_nce_ << " selected from a grid centred on "
                      << centre << " (median cosine " << best_median << ")" << std::endl;
    }

    // ---- predictions for every distinct peptide ----
    startProgress(0, 2, "Predicting fragment intensities and retention times...");
    const std::vector<float> nces(uniq_seq.size(), static_cast<float>(used_nce_));
    const std::vector<int64_t> instruments(uniq_seq.size(), instrument);
    const std::vector<std::vector<float>> ms2_pred = ms2.predictMS2(uniq_seq, uniq_charge, nces, instruments);
    setProgress(1);
    PeptDeepRTInference rt(rt_model_, threads_, batch_size_);
    const std::vector<float> rt_pred_uniq = rt.predictRT(uniq_seq);
    setProgress(2);
    endProgress();

    // ---- MS2 agreement ----
    std::vector<Similarity> sims(rows.size());
    for (Size i = 0; i < rows.size(); ++i)
    {
      const Row& r = rows[i];
      const std::vector<float>& pred = ms2_pred[r.key];
      Size channels = (r.len > 1 && !pred.empty()) ? pred.size() / (r.len - 1) : DEFAULT_ION_CHANNELS;
      if (channels == 0) { channels = DEFAULT_ION_CHANNELS; }
      const PeptideHit& hit = peptide_ids[r.pi].getHits()[r.hit];
      sims[i] = similarity_(predictedIons_(pred, r.len, channels), observedIons_(hit),
                            ms2_min_ordinal_, r.len, ms2_strong_fraction_);
    }

    // ---- retention time calibration ----
    // Predicted RT is in the model's own units; map it onto this run's gradient before
    // taking a residual, otherwise the feature measures the unit mismatch.
    TransformationModel::DataPoints points;
    points.reserve(n_conf);
    for (Size i = 0; i < rows.size(); ++i)
    {
      if (!is_confident[i]) { continue; }
      points.emplace_back(static_cast<double>(rt_pred_uniq[rows[i].key]), peptide_ids[rows[i].pi].getRT());
    }

    std::vector<double> rt_err(rows.size(), 0.0);
    std::unique_ptr<TransformationModel> model;
    if (points.size() >= 4)
    {
      Param mp;
      try
      {
        if (rt_model_type_ == "b_spline")
        {
          TransformationModelBSpline::getDefaultParameters(mp);
          mp.setValue("num_nodes", static_cast<int>(rt_num_nodes_));
          model = std::make_unique<TransformationModelBSpline>(points, mp);
        }
        else if (rt_model_type_ == "lowess")
        {
          TransformationModelLowess::getDefaultParameters(mp);
          model = std::make_unique<TransformationModelLowess>(points, mp);
        }
        else
        {
          TransformationModelLinear::getDefaultParameters(mp);
          model = std::make_unique<TransformationModelLinear>(points, mp);
        }
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "[PeptDeepRescoring] retention-time calibration (" << rt_model_type_
                        << ") failed (" << e.getName() << "); falling back to a linear fit." << std::endl;
        model.reset();
      }
      if (!model)
      {
        Param lp;
        TransformationModelLinear::getDefaultParameters(lp);
        model = std::make_unique<TransformationModelLinear>(points, lp);
      }

      std::vector<double> conf_residuals;
      conf_residuals.reserve(points.size());
      for (Size i = 0; i < rows.size(); ++i)
      {
        const double fitted = model->evaluate(static_cast<double>(rt_pred_uniq[rows[i].key]));
        rt_err[i] = std::abs(peptide_ids[rows[i].pi].getRT() - fitted);
        if (is_confident[i]) { conf_residuals.push_back(rt_err[i]); }
      }
      if (!conf_residuals.empty())
      {
        std::nth_element(conf_residuals.begin(), conf_residuals.begin() + conf_residuals.size() / 2, conf_residuals.end());
        rt_calibration_error_ = conf_residuals[conf_residuals.size() / 2];
      }
      OPENMS_LOG_INFO << "[PeptDeepRescoring] retention-time calibration (" << rt_model_type_ << ") on "
                      << points.size() << " PSMs, median residual " << rt_calibration_error_ << " s" << std::endl;
    }
    else
    {
      OPENMS_LOG_WARN << "[PeptDeepRescoring] too few PSMs to calibrate retention time; "
                      << Constants::UserParam::RT_ABS_ERROR << " will be 0 for every PSM." << std::endl;
    }

    // ---- write the features back ----
    // Every hit gets every feature: a feature present on only some PSMs of a run is
    // dropped for all of them by the Percolator feature-set check.
    std::map<std::pair<Size, Size>, Size> row_of;
    for (Size i = 0; i < rows.size(); ++i) { row_of[{rows[i].pi, rows[i].hit}] = i; }

    for (Size pi = 0; pi < peptide_ids.size(); ++pi)
    {
      std::vector<PeptideHit> hits = peptide_ids[pi].getHits();
      for (Size h = 0; h < hits.size(); ++h)
      {
        const auto it = row_of.find({pi, h});
        const Similarity s = (it == row_of.end()) ? Similarity{} : sims[it->second];
        const double e = (it == row_of.end()) ? 0.0 : rt_err[it->second];
        hits[h].setMetaValue(Constants::UserParam::MS2_COSINE, s.cosine);
        hits[h].setMetaValue(Constants::UserParam::MS2_SPECTRAL_ANGLE, s.spectral_angle);
        hits[h].setMetaValue(Constants::UserParam::MS2_PEARSON, s.pearson);
        hits[h].setMetaValue(Constants::UserParam::MS2_FRAC_PRED_FOUND, s.frac_found);
        hits[h].setMetaValue(Constants::UserParam::RT_ABS_ERROR, e);
      }
      peptide_ids[pi].setHits(std::move(hits));
    }

    if (!protein_ids.empty())
    {
      ProteinIdentification::SearchParameters sp = protein_ids[0].getSearchParameters();
      StringList features = sp.metaValueExists("extra_features")
        ? ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString(), ',')
        : StringList();
      for (const std::string& f : featureNames())
      {
        if (std::find(features.begin(), features.end(), f) == features.end()) { features.push_back(f); }
      }
      sp.setMetaValue("extra_features", ListUtils::concatenate(features, ","));
      protein_ids[0].setSearchParameters(sp);
    }
  }

} // namespace OpenMS
