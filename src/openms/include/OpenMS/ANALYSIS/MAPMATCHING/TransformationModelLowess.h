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
    @brief Parse a user-supplied span grid string into candidate values in (0, 1).

    Accepts comma and/or whitespace as separators (e.g. `"0.2, 0.3 0.5"`).
                           Values `<= 0` or `>= 1` are ignored.

    @param s Input string containing candidate spans.
    @return Vector of valid candidate spans in (0, 1).
 */
  static std::vector<double> parseSpanGridString_(const String& s);

  /**
    @brief Build the candidate LOWESS span grid with neighbor-based flooring and min/max clamping.

    Rules applied:
      - If @p grid_str is empty, use a sensible default grid.
      - Enforce a lower bound based on both @p span_min_param and the minimum number
        of neighbors in the window: `lower = max(span_min_param, min(0.9, min_neighbors/n_pts))`.
      - Clamp candidates to `[lower, span_max_param]`, with an additional safety cap at `0.99`.
      - Sort and de-duplicate candidates with a small tolerance.
      - Guarantee a non-empty result by falling back to a single value if necessary.

    @param n_pts           Number of data points available.
    @param grid_str        Optional user-specified grid string (comma/space separated).
    @param span_min_param  Minimum allowed span (will be raised to at least `0.01`).
    @param span_max_param  Maximum allowed span (will be lowered to at most `0.99`).
    @param min_neighbors   Minimum number of neighbors to include in each local regression.

    @return Vector of candidate spans after clamping and de-duplication.
  */
  static std::vector<double> buildSpanGrid(Size n_pts,
                                           const String& grid_str,
                                           double span_min_param,
                                           double span_max_param,
                                           int min_neighbors);

  /**
    @brief Score absolute residuals according to the selected metric.

    Supported @p metric values:
      - `"rmse"`: root mean square of errors (vs zero)
      - `"mae"` : mean absolute error
      - `"p90"`, `"p95"`, `"p99"`: corresponding empirical quantiles of |errors|

    If @p errs is empty, returns `+inf` to signal an unusable model fit.

    @param errs   Absolute residuals (validation errors).
    @param metric Metric name (`"rmse"`, `"mae"`, `"p90"`, `"p95"`, `"p99"`).
    @return Scalar loss; lower values indicate better fit.
  */
  static double scoreResiduals(const std::vector<double>& errs, const String& metric);

  /// Numerical tie tolerance for CV score comparisons (absolute difference)
  static constexpr double kTieTol = 1e-12;

  /// @endcond


protected:
    /// Pointer to the underlying interpolation
    TransformationModelInterpolated* model_;

  };
} // namespace

