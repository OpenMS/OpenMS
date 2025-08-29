// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/PROCESSING/SMOOTHING/FastLowessSmoothing.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/CROSSVALIDATION/CrossValidation.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>     
#include <sstream>    

using namespace std;

namespace OpenMS
{
  bool cmpFirstDimension(const TransformationModel::DataPoint& x, const TransformationModel::DataPoint& y)
  {
    return (x.first < y.first);
  }

  TransformationModelLowess::TransformationModelLowess(
      const TransformationModel::DataPoints& data_,
      const Param& params) : model_(nullptr)
  {
    // parameter handling/checking:
    params_ = params;
    Param defaults;
    getDefaultParameters(defaults);
    params_.setDefaults(defaults);

    if (data_.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "'lowess' model requires more data");
    }

    // TODO copy ... 
    TransformationModel::DataPoints data(data_);

    // sort data
    std::sort(data.begin(), data.end(), cmpFirstDimension);

    vector<double> x(data.size()), y(data.size()), result(data.size());
    double xmin_ = data[0].first;
    double xmax_ = xmin_;
    for (Size i = 0; i < data.size(); ++i)
    {
      x[i] = data[i].first;
      y[i] = data[i].second;
      if (x[i] < xmin_) 
      {
        xmin_ = x[i];
      }
      else if (x[i] > xmax_)
      {
        xmax_ = x[i];
      }
    }

    double span = params_.getValue("span");
    int nsteps = params_.getValue("num_iterations");
    double delta = params_.getValue("delta");
    
    if (delta < 0.0)
    {
      delta = (xmax_ - xmin_) * 0.01; // automatically determine delta
    }

    // Auto-span options
    const bool auto_span_flag = params_.getValue("auto_span").toBool();
    const double span_min_param = (double)params_.getValue("auto_span_min");
    const double span_max_param = (double)params_.getValue("auto_span_max");
    const int    min_neighbors  = (int)params_.getValue("auto_min_neighbors");
    const int    k_folds_param  = (int)params_.getValue("auto_k_folds");
    const auto metric = (TransformationModelLowess::CVMetric)Helpers::indexOf(TransformationModelLowess::names_of_cvmetric, params_.getValue("auto_metric").toString());

    // Determine optimal span using CV
    const Size n = data.size();
    if (auto_span_flag)
    {
      OPENMS_LOG_INFO << "Will perform CV to determine optimal span for lowess fit..." << std::endl;

      // Build folds: LOO (n<=50) else K-fold
      const bool use_loo = (n <= 50);
      const Size K = use_loo ? n : static_cast<Size>(std::max(2, k_folds_param));
      const auto folds = OpenMS::CrossValidation::makeKFolds(n, K);

      // Build candidate grid spans
      const String grid_str = params_.getValue("auto_span_grid").toString();
      const std::vector<double> spans = buildSpanGrid(n, grid_str,
                                                      span_min_param, span_max_param,
                                                      min_neighbors);

      // Train evaluation callback
      auto train_eval = [&](double s,
                            const std::vector<std::vector<Size>>& folds_in,
                            std::vector<double>& abs_errs)
      {
        if (s * static_cast<double>(n) < static_cast<double>(min_neighbors)) return;

        Param p_cv = params_;
        p_cv.setValue("span", s);
        p_cv.setValue("auto_span", "false");

        for (Size f = 0; f < folds_in.size(); ++f)
        {
          TransformationModel::DataPoints train;
          train.reserve(n - folds_in[f].size());
          std::vector<char> held(n, 0);
          for (Size j : folds_in[f]) held[j] = 1;
          for (Size j = 0; j < n; ++j) if (!held[j]) train.push_back(data[j]);

          if (train.size() < static_cast<Size>(std::max(3, min_neighbors))) continue;

          TransformationModelLowess cv_model(train, p_cv);

          for (Size j : folds_in[f])
          {
            const double yhat = cv_model.evaluate(x[j]);
            const double e    = std::fabs(yhat - y[j]);
            if (std::isfinite(e)) abs_errs.push_back(e);
          }
        }
      };

      // Later, when scoring folds:
      auto score = [&](const std::vector<double>& errs)
      {
        return scoreResiduals(errs, metric);
      };

      // Run 1-D grid search
      const bool prefer_larger = true;
      const auto [best_span, best_score] =
        OpenMS::CrossValidation::gridSearch1D(spans.begin(), spans.end(),
                                              folds, train_eval, score,
                                              kTieTol, prefer_larger);

      span = best_span;
      OPENMS_LOG_INFO << "Optimal selected span=" << span
                << " (" << params_.getValue("auto_metric").toString() << " = " << best_score << ")" << std::endl;

      // persist for downstream and prevent re-entry
      params_.setValue("span", span);
      params_.setValue("auto_span", "false");
    }

    FastLowessSmoothing::lowess(x, y, span, nsteps, delta, result);

    TransformationModel::DataPoints data_out;
    for (Size i = 0; i < result.size(); ++i)
    {
      data_out.push_back( std::make_pair(x[i], result[i]) );
    }

    // TODO thin out data here ? we may not need that many points here to interpolate ...  it is enough if we store a few datapoints

    Param p;
    TransformationModelInterpolated::getDefaultParameters(p);
    /// p.setValue("interpolation_type", "cspline"); // linear interpolation between lowess pts
    /// p.setValue("extrapolation_type", "four-point-linear");
    p.setValue("interpolation_type", params_.getValue("interpolation_type"));
    p.setValue("extrapolation_type", params_.getValue("extrapolation_type"));

    // create new interpolation model based on the lowess data
    model_ = new TransformationModelInterpolated(data_out, p);
  }

  TransformationModelLowess::~TransformationModelLowess()
  {
    if (model_) delete model_;
  }

  void TransformationModelLowess::getDefaultParameters(Param& params)
  {
    params.clear();
    params.setValue("span", 2/3.0, "Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.");
    params.setMinFloat("span", 0.0);
    params.setMaxFloat("span", 1.0);

    // --- Auto-span tuning (optional) ---
    params.setValue("auto_span", "false",
                    "If true, or if 'span' is 0, automatically select LOWESS span by cross-validation.");
    params.setValidStrings("auto_span", {"true","false"});

    params.setValue("auto_span_min", 0.15,
                    "Lower bound for auto-selected span.");
    params.setMinFloat("auto_span_min", 0.001);

    params.setValue("auto_span_max", 0.80,
                    "Upper bound for auto-selected span.");
    params.setMaxFloat("auto_span_max", 0.99);

    params.setValue("auto_min_neighbors", 5,
                    "Minimum number of neighbors (span*n) enforced in auto mode.");
    params.setMinInt("auto_min_neighbors", 3);

    params.setValue("auto_k_folds", 5,
                    "K-folds for CV when n>50 (else LOO is used).");
    params.setMinInt("auto_k_folds", 2);

    params.setValue("auto_metric", "mae",
                    "Metric for CV selection: one of {'p90','p95','p99','rmse','mae'}.");
    params.setValidStrings("auto_metric", {"p90","p95","p99","rmse","mae"});

    params.setValue("auto_span_grid", "",
                    "Optional explicit grid of span candidates in (0,1]. Comma-separated list, e.g. '0.2,0.3,0.5'.  If empty, a default grid is used.");


    params.setValue("num_iterations", 3, "Number of robustifying iterations for lowess fitting.");
    params.setMinInt("num_iterations", 0);

    params.setValue("delta", -1.0, "Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.");

    params.setValue("interpolation_type", "cspline", "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation");
    params.setValidStrings("interpolation_type", {"linear","cspline","akima"});

    params.setValue("extrapolation_type", "four-point-linear", "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation.");
    params.setValidStrings("extrapolation_type", {"two-point-linear","four-point-linear","global-linear"});
  }

  std::vector<double> TransformationModelLowess::buildSpanGrid(Size n_pts,
                                    const String& grid_str,
                                    double span_min_param,
                                    double span_max_param,
                                    int min_neighbors)
  {
    // Parse user grid if supplied, needs to be comma-separated doubles
    std::vector<double> grid;
    if (!grid_str.empty())
    {
      grid = ListUtils::create<double>(grid_str);
    }
    else
    {
      static const double cand[] = {0.01, 0.05, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
      grid.assign(std::begin(cand), std::end(cand));
    }

    double span_min = std::max(0.01, span_min_param);
    double span_max = std::min(0.99, span_max_param);
    if (span_min > span_max) std::swap(span_min, span_max);

    const double min_span_neighbors = (n_pts > 0) ? static_cast<double>(min_neighbors) / static_cast<double>(n_pts) : 1.0;
    const double lower = std::max(span_min, std::min(0.9, min_span_neighbors)); // avoid 1–2 neighbor fits

    for (double& v : grid)
    {
      if (v < lower) v = lower;
      if (v > span_max) v = span_max;
    }
    std::sort(grid.begin(), grid.end());
    grid.erase(std::unique(grid.begin(), grid.end(),
                           [](double a, double b){ return std::fabs(a - b) < 1e-9; }),
               grid.end());

    if (grid.empty())
    {
      grid.push_back(std::min(0.95, std::max(0.01, lower)));
    }
    return grid;
  }

  const std::array<std::string, static_cast<size_t>(TransformationModelLowess::CVMetric::SIZE_OF_CVMETRIC)> TransformationModelLowess::names_of_cvmetric = { "rmse", "mae", "p90", "p95", "p99" };

  double TransformationModelLowess::scoreResiduals(const std::vector<double>& errs,
                                                   CVMetric metric)
  {
    if (errs.empty()) return std::numeric_limits<double>::infinity();

    switch (metric)
    {
      case CVMetric::RMSE:
      {
        std::vector<double> zeros(errs.size(), 0.0);
        return OpenMS::Math::rootMeanSquareError(errs.begin(), errs.end(), zeros.begin(), zeros.end());
      }
      case CVMetric::MAE:
        return OpenMS::Math::MeanAbsoluteDeviation(errs.begin(), errs.end(), 0.0);

      case CVMetric::P90:
      {
        auto tmp = errs;
        return OpenMS::Math::quantileUnsortedInput(tmp.begin(), tmp.end(), 0.90);
      }
      case CVMetric::P99:
      {
        auto tmp = errs;
        return OpenMS::Math::quantileUnsortedInput(tmp.begin(), tmp.end(), 0.99);
      }
      case CVMetric::P95:
      default:
      {
        auto tmp = errs;
        return OpenMS::Math::quantileUnsortedInput(tmp.begin(), tmp.end(), 0.95);
      }
    }
  }

  }
