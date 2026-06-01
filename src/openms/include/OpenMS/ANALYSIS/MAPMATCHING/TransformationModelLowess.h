// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // is this needed?
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

namespace OpenMS
{

  /**
    @brief Lowess (non-linear) model for transformations

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelLowess :
    public TransformationModel
  {
public:
    /**
      @brief Constructor

      @exception Exception::IllegalArgument is thrown if too few data points are provided.
    */
    TransformationModelLowess(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelLowess() override;

    /// Evaluates the model at the given value
    double evaluate(double value) const override
    {
      return model_->evaluate(value);
    }

    using TransformationModel::getParameters;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

private:
  // Hide helper utilities from public API. Defined in the .cpp.

  /// @cond INTERNAL

  /**
    @brief Build the candidate LOWESS span grid with neighbor-based flooring and min/max clamping.

    Rules applied:
      - If @p grid_str is empty, use a sensible default grid.
      - Enforce a lower bound based on both @p span_min_param and the minimum number
        of neighbors in the window: `lower = max(span_min_param, min(0.9, min_neighbors/n_pts))`.
      - Clamp candidates to `[lower, span_max_param]`, with an additional safety cap at `0.99`.
      - Sort and de-duplicate candidates with a small tolerance.
      - Guarantee a non-empty result by falling back to a single value if necessary.

    @param[in] n_pts           Number of data points available.
    @param[in] grid_str        Span candidates in (0, 1). If empty, a default span grid is used.
    @param[in] span_min_param  Minimum allowed span (will be raised to at least `0.01`).
    @param[in] span_max_param  Maximum allowed span (will be lowered to at most `0.99`).
    @param[in] min_neighbors   Minimum number of neighbors to include in each local regression.

    @return Vector of candidate spans after clamping and de-duplication.
  */
  static std::vector<double> buildSpanGrid(Size n_pts,
                                           const std::vector<double>& candidate_spans,
                                           double span_min_param,
                                           double span_max_param,
                                           int min_neighbors);

  /// Cross-validation metric for auto-span selection
  enum class CVMetric
  {
    RMSE,
    MAE,
    P90,
    P95,
    P99,
    SIZE_OF_CVMETRIC
  };
  static const std::array<std::string, (Size)CVMetric::SIZE_OF_CVMETRIC> names_of_cvmetric;

  /**
    @brief Score absolute residuals according to the selected metric.

    Supported @p metric values:
      - `"rmse"`: root mean square of errors (vs zero)
      - `"mae"` : mean absolute error
      - `"p90"`, `"p95"`, `"p99"`: corresponding empirical quantiles of |errors|

    If @p errs is empty, returns `+inf` to signal an unusable model fit.

    @param[in] errs   Absolute residuals (validation errors).
    @param[in] metric Metric name (`"rmse"`, `"mae"`, `"p90"`, `"p95"`, `"p99"`).
    @return Scalar loss; lower values indicate better fit.
  */
  static double scoreResiduals(const std::vector<double>& errs, CVMetric metric);

  /// Numerical tie tolerance for CV score comparisons (absolute difference)
  static constexpr double kTieTol = 1e-12;

  /// @endcond


protected:
    /// Pointer to the underlying interpolation
    TransformationModelInterpolated* model_;

  };
} // namespace

