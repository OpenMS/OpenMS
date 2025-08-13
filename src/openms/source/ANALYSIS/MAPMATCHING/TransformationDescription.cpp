// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

#include <iomanip>
#include <iostream>

using namespace std;

namespace OpenMS
{
  TransformationDescription::TransformationDescription() :
    data_(TransformationDescription::DataPoints()),
    model_type_("none"),
    model_(new TransformationModel())
  {
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription::DataPoints& data) :
    data_(data), model_type_("none"),
    model_(new TransformationModel())
  {
  }

  TransformationDescription::~TransformationDescription()
  {
    delete model_;
  }

  TransformationDescription::TransformationDescription(
    const TransformationDescription& rhs)
  {
    data_ = rhs.data_;
    model_type_ = "none";
    model_ = nullptr; // initialize this before the "delete" call in "fitModel"!
    Param params = rhs.getModelParameters();
    fitModel(rhs.model_type_, params);
  }

  TransformationDescription& TransformationDescription::operator=(
    const TransformationDescription& rhs)
  {
    if (this == &rhs)
      return *this;

    data_ = rhs.data_;
    model_type_ = "none";
    Param params = rhs.getModelParameters();
    fitModel(rhs.model_type_, params);

    return *this;
  }

  void TransformationDescription::fitModel(const String& model_type,
                                           const Param& params)
  {
    // if (previous) transformation is the identity, don't fit another model:
    if (model_type_ == "identity") return;

    delete model_;
    model_ = nullptr; // avoid segmentation fault in case of exception
    if ((model_type == "none") || (model_type == "identity"))
    {
      model_ = new TransformationModel();
    }
    else if (model_type == "linear")
    {
      model_ = new TransformationModelLinear(data_, params);
      // // debug output:
      // double slope, intercept;
      // TransformationModelLinear* lm = dynamic_cast<TransformationModelLinear*>(model_);
      // lm->getParameters(slope, intercept);
      // cout << "slope: " << slope << ", intercept: " << intercept << endl;
    }
    else if (model_type == "b_spline")
    {
      model_ = new TransformationModelBSpline(data_, params);
    }
    else if (model_type == "lowess")
    {
      model_ = new TransformationModelLowess(data_, params);
    }
    else if (model_type == "interpolated")
    {
      model_ = new TransformationModelInterpolated(data_, params);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unknown model type '" + model_type + "'");
    }
    model_type_ = model_type;
  }

  double TransformationDescription::apply(double value) const
  {
    return model_->evaluate(value);
  }

  const String& TransformationDescription::getModelType() const
  {
    return model_type_;
  }

  void TransformationDescription::getModelTypes(StringList& result)
  {
    result = ListUtils::create<String>("linear,b_spline,interpolated,lowess");
    // "none" and "identity" don't count
  }

  void TransformationDescription::setDataPoints(const DataPoints& data)
  {
    data_ = data;
    model_type_ = "none"; // reset the model even if it was "identity"
    delete model_;
    model_ = new TransformationModel();
  }

  void TransformationDescription::setDataPoints(const vector<pair<double, double> >& data)
  {
    data_.resize(data.size());
    for (Size i = 0; i < data.size(); ++i)
    {
      data_[i] = data[i];
    }
    model_type_ = "none"; // reset the model even if it was "identity"
    delete model_;
    model_ = new TransformationModel();
  }

  const TransformationDescription::DataPoints&
  TransformationDescription::getDataPoints() const
  {
    return data_;
  }

  const Param& TransformationDescription::getModelParameters() const
  {
    return model_->getParameters();
  }

  void TransformationDescription::invert()
  {
    for (TransformationDescription::DataPoints::iterator it = data_.begin();
         it != data_.end(); ++it)
    {
      *it = TransformationDescription::DataPoint(it->second, it->first,
                                                 it->note);
    }
    // ugly hack for linear model with explicit slope/intercept parameters:
    if ((model_type_ == "linear") && data_.empty())
    {
      TransformationModelLinear* lm =
        dynamic_cast<TransformationModelLinear*>(model_);
      lm->invert();
    }
    else
    {
      Param params = getModelParameters();
      fitModel(model_type_, params);
    }
  }

  void TransformationDescription::getDeviations(vector<double>& diffs, 
                                                bool do_apply,
                                                bool do_sort) const
  {
    diffs.clear();
    diffs.reserve(data_.size());
    for (DataPoints::const_iterator it = data_.begin(); it != data_.end(); ++it)
    {
      double x = it->first;
      if (do_apply) x = apply(x);
      diffs.push_back(abs(x - it->second));
    }
    if (do_sort) sort(diffs.begin(), diffs.end());
  }

  TransformationDescription::TransformationStatistics TransformationDescription::getStatistics() const
  {
    TransformationDescription::TransformationStatistics s;

    if (data_.empty()) return s;

    // x/y data ranges:
    double xmin, xmax, ymin, ymax;
    xmin = xmax = data_[0].first;
    ymin = ymax = data_[0].second;
    for (DataPoints::const_iterator it = ++data_.begin(); it != data_.end();
         ++it)
    {
      if (xmin > it->first) xmin = it->first;
      if (xmax < it->first) xmax = it->first;
      if (ymin > it->second) ymin = it->second;
      if (ymax < it->second) ymax = it->second;
    }

    s.xmin = xmin;
    s.xmax = xmax;
    s.ymin = ymin;
    s.ymax = ymax;

    // deviations:
    vector<double> diffs;
    getDeviations(diffs);
    bool no_model = (model_type_ == "none") || (model_type_ == "identity");

    for (const Size p : s.percents)
    {
      Size index = p / 100.0 * diffs.size() - 1;
      s.percentiles_before[p] = diffs[index];
    }

    // if we have a model, calculate deviations after applying the model
    // else set the same values
    if (!no_model) { getDeviations(diffs, true); }
  
    for (const Size p : s.percents)
    {
      Size index = p / 100.0 * diffs.size() - 1;
      s.percentiles_after[p] = diffs[index];
    }

    return s;
  }

  void TransformationDescription::printSummary(ostream& os) const
  {
    const TransformationStatistics s = getStatistics();

    os << "Number of data points (x/y pairs): " << data_.size() << "\n";
    if (data_.empty()) return;

    os << "Data range (x): " << s.xmin << " to " << s.xmax
       << "\nData range (y): " << s.ymin << " to " << s.ymax << "\n";

    // deviations:
    vector<double> diffs;
    getDeviations(diffs);
    bool no_model = (model_type_ == "none") || (model_type_ == "identity");
    os << String("Summary of x/y deviations") +
      (no_model ? "" : " before transformation") + ":\n";

    for (Size p : s.percents)
    {
      os << "- " << setw(3) << p << "% of data points within (+/-)"
         << s.percentiles_before.at(p) << "\n";
    }
    if (no_model)
    {
      os << endl;
      return;
    }
    // else:
    getDeviations(diffs, true);
    os << "Summary of x/y deviations after applying '" << model_type_
       << "' transformation:\n";

    for (Size p : s.percents)
    {
      os << "- " << setw(3) << p << "% of data points within (+/-)"
         << s.percentiles_after.at(p) << "\n";
    }
    os << endl;
  }


  double TransformationDescription::estimateWindow(double quantile,
                                                   bool invert,
                                                   bool full_window) const
  {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "[estimateWindow] enter: quantile=" << quantile
              << " invert=" << (invert ? "true" : "false")
              << " full_window=" << (full_window ? "true" : "false") << "\n";

    // Work on a copy so we don't mutate the original
    TransformationDescription tmp(*this);
    std::cout << "[estimateWindow] original: model=" << this->getModelType()
              << " n_points=" << this->getDataPoints().size() << "\n";

    if (invert)
    {
      std::cout << "[estimateWindow] action: invert() (swap datapoints + refit inverse model)\n";
      tmp.invert();
      std::cout << "[estimateWindow] after invert: model=" << tmp.getModelType()
                << " n_points=" << tmp.getDataPoints().size() << "\n";
    }

    // Helper lambdas
    auto quantile_idx = [](size_t n, double q) -> size_t {
      if (n == 0) return 0;
      size_t i = static_cast<size_t>(std::floor(q * static_cast<double>(n)));
      if (i >= n) i = n - 1;
      return i;
    };

    auto summarize = [](const std::vector<double>& v) -> std::string {
      if (v.empty()) return "empty";
      std::vector<double> s = v;
      std::sort(s.begin(), s.end());
      auto getq = [&](double q) {
        size_t i = std::min(s.size() - 1, static_cast<size_t>(std::floor(q * s.size())));
        return s[i];
      };
      std::ostringstream oss;
      oss << "n=" << s.size()
          << " p50=" << getq(0.50)
          << " p90=" << getq(0.90)
          << " p95=" << getq(0.95)
          << " p99=" << getq(0.99)
          << " max=" << s.back();
      return oss.str();
    };

    auto all_almost_zero = [](const std::vector<double>& v) -> bool {
      double mx = 0.0;
      for (double d : v) if (std::isfinite(d) && d > mx) mx = d;
      return mx <= 1e-6; // ~1 microsecond tolerance
    };

    // 1) In-sample residuals: | y - f(x) | in current tmp orientation
    std::vector<double> diffs;
    tmp.getDeviations(diffs, /*do_apply=*/true, /*do_sort=*/false);
    const size_t n0 = diffs.size();

    // Drop NaN/Inf if any
    diffs.erase(std::remove_if(diffs.begin(), diffs.end(),
                               [](double v){ return !std::isfinite(v); }),
                diffs.end());
    const size_t n_dropped = n0 - diffs.size();

    std::cout << "[estimateWindow] in-sample residuals: "
              << "n=" << diffs.size()
              << " dropped_nonfinite=" << n_dropped << "\n";
    std::cout << "[estimateWindow] in-sample summary (sec): " << summarize(diffs) << "\n";

    double full = 0.0;
    bool used_in_sample = false;
    bool used_loo = false;
    bool used_fallback = false;

    if (!diffs.empty() && !all_almost_zero(diffs))
    {
      std::sort(diffs.begin(), diffs.end());
      const size_t n = diffs.size();
      const size_t idx = quantile_idx(n, quantile);
      const double half = diffs[idx];
      full = full_window ? (2.0 * half) : half;

      used_in_sample = true;
      std::cout << "[estimateWindow] using IN-SAMPLE residuals: "
                << "idx=" << idx << " half=" << half << " -> full=" << full << " sec\n";
    }
    else
    {
      std::cout << "[estimateWindow] in-sample residuals ~zero or empty; trying LOO residuals...\n";

      // 2) LOO residuals (prediction error) in RT space
      const auto& pts = tmp.getDataPoints(); // (iRT, RT) if invert==true
      const String mt = tmp.getModelType();
      Param params = tmp.getModelParameters();

      std::vector<double> diffs_loo;
      diffs_loo.reserve(pts.size());

      for (size_t i = 0; i < pts.size(); ++i)
      {
        // Build fold without point i
        DataPoints fold;
        fold.reserve(pts.size() > 0 ? pts.size() - 1 : 0);
        for (size_t j = 0; j < pts.size(); ++j) if (j != i) fold.push_back(pts[j]);

        if (fold.size() < 2)
        {
          diffs_loo.push_back(0.0);
          continue;
        }

        TransformationDescription cv(fold);
        cv.fitModel(mt, params);
        const double y_pred = cv.apply(pts[i].first);
        const double err = std::abs(y_pred - pts[i].second);
        diffs_loo.push_back(std::isfinite(err) ? err : 0.0);
      }

      std::cout << "[estimateWindow] LOO summary (sec): " << summarize(diffs_loo) << "\n";

      if (!diffs_loo.empty() && !all_almost_zero(diffs_loo))
      {
        std::sort(diffs_loo.begin(), diffs_loo.end());
        const size_t n = diffs_loo.size();
        const double q_eff = (n >= 30 ? quantile
                                      : std::min(0.95, std::max(0.85, quantile)));
        const size_t idx = quantile_idx(n, q_eff);
        const double half = diffs_loo[idx];
        full = full_window ? (2.0 * half) : half;

        used_loo = true;
        std::cout << "[estimateWindow] using LOO residuals: "
                  << "q_eff=" << q_eff << " idx=" << idx
                  << " half=" << half << " -> full=" << full << " sec\n";
      }
      else
      {
        // 3) Final fallback: small fraction of RT span (seconds) of anchors
        double y_min =  std::numeric_limits<double>::infinity();
        double y_max = -std::numeric_limits<double>::infinity();
        for (const auto& p : pts)
        {
          if (p.second < y_min) y_min = p.second;
          if (p.second > y_max) y_max = p.second;
        }
        const double rt_span = (std::isfinite(y_min) && std::isfinite(y_max)) ? (y_max - y_min) : 0.0;

        const double alpha = 0.015;   // 1.5% of span
        const double min_full_sec = 30.0;
        const double max_full_sec = 600.0;

        double guess = (rt_span > 0.0) ? (alpha * rt_span) : 0.0;
        double unclamped = guess;

        if (!std::isfinite(guess) || guess <= 0.0) guess = min_full_sec;
        if (guess < min_full_sec)    guess = min_full_sec;
        if (guess > max_full_sec)    guess = max_full_sec;

        full = guess;
        used_fallback = true;

        std::cout << "[estimateWindow] using FALLBACK alpha*span: "
                  << "rt_span=" << rt_span << " alpha=" << alpha
                  << " raw=" << unclamped << " clamped=" << full << " sec\n";
      }
    }

    // Final sanity for invert==true (RT space): small floor & large cap
    if (invert)
    {
      const double floor_sec = 30.0;
      const double cap_sec   = 1800.0; // 30 min max
      const double before = full;
      if (!std::isfinite(full) || full <= 0.0) full = floor_sec;
      if (full < floor_sec) full = floor_sec;
      if (full > cap_sec)   full = cap_sec;

      if (before != full)
      {
        std::cout << "[estimateWindow] clamp (RT-space): " << before
                  << " -> " << full << " sec\n";
      }
    }

    std::cout << "[estimateWindow] exit: "
              << (used_in_sample ? "IN-SAMPLE" : used_loo ? "LOO" : "FALLBACK")
              << " full_window=" << full << " sec\n";

    return full;
  }


} // end of namespace OpenMS
