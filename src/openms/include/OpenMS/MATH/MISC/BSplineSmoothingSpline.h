// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
     * @brief Smoothing spline implementation using B-spline basis
     * 
     * This class implements a smoothing spline that uses B-spline basis functions
     * with automatic knot selection based on a smoothing parameter, similar to
     * scipy's UnivariateSpline.
     * 
     * Unlike the BSpline2d class which requires explicit node specification,
     * this class implements the scipy-style optimization:
     *   minimize: RSS(f) + λ * ∫(f''(x))² dx
     * 
     * Or equivalently with smoothing parameter s:
     *   minimize: ∫(f''(x))² dx
     *   subject to: RSS(f) ≤ s
     * 
     * The key difference from BSpline2d:
     * - BSpline2d: Minimizes RSS for a fixed number of nodes (interpolating spline)
     * - BSplineSmoothingSpline: Balances RSS vs smoothness (smoothing spline)
     * 
     * For compatibility with PyProphet's UnivariateSpline:
     * - Uses automatic knot selection based on smoothing parameter
     * - Default s = m - sqrt(2*m) where m is number of data points
     * - Allows curve to deviate from data points (controlled smoothing)
     * 
     * Implementation strategy:
     * Since we cannot modify the eol-bspline library directly, we implement
     * the smoothing spline optimization by:
     * 1. Testing different wavelength/node combinations in BSpline2d
     * 2. Finding the combination that best satisfies the smoothing constraint
     * 3. Using the resulting spline for evaluation
     * 
     * @ingroup Math
     */
    class OPENMS_DLLAPI BSplineSmoothingSpline
    {
    public:
      /**
       * @brief Construct a smoothing spline with automatic knot selection
       * 
       * @param x Data points (must be sorted)
       * @param y Values at data points  
       * @param s Smoothing parameter (RSS target). If negative, uses automatic: s = m - sqrt(2*m)
       * @param k Spline degree (default: 3 for cubic)
       * 
       * The smoothing parameter s controls the trade-off between fit and smoothness:
       * - s = 0: Interpolating spline (passes through all points)
       * - s > 0: Allows deviation from data points
       * - s = m - sqrt(2*m): Default used by scipy.UnivariateSpline
       * - s → ∞: Approaches simple polynomial fit
       * 
       * @note For small datasets (n < 4), falls back to interpolation
       */
      BSplineSmoothingSpline(const std::vector<double>& x,
                             const std::vector<double>& y,
                             double s = -1.0,
                             int k = 3);

      /**
       * @brief Evaluate the smoothing spline at x
       * 
       * @param x Point at which to evaluate
       * @return Smoothed value, or NaN if fit failed
       */
      double eval(double x) const;

      /**
       * @brief Check if spline fit was successful
       */
      bool ok() const { return ok_; }

      /**
       * @brief Get the number of interior knots selected
       */
      int num_interior_knots() const { return num_interior_knots_; }

      /**
       * @brief Get the achieved RSS
       */
      double rss() const { return rss_; }

      /**
       * @brief Get the smoothing parameter used
       */
      double smoothing_param() const { return s_; }

    private:
      bool ok_ = false;
      int num_interior_knots_ = 0;
      double rss_ = 0.0;
      double s_ = 0.0;
      int k_ = 3;
      std::vector<double> x_;
      std::vector<double> y_;
      
      // Fit type
      enum FitType { POLYNOMIAL, BSPLINE };
      FitType fit_type_ = BSPLINE;
      
      // Polynomial coefficients [a0, a1, a2, a3] for a0 + a1*x + a2*x^2 + a3*x^3
      std::vector<double> poly_coeffs_;
      
      // The fitted BSpline (owned, only used if fit_type_ == BSPLINE)
      BSpline2d* spline_ = nullptr;

      /**
       * @brief Find optimal wavelength/nodes combination to match smoothing target
       * 
       * Since BSpline2d doesn't support direct smoothing parameter s,
       * we search for a wavelength that produces RSS close to the target s.
       * 
       * Strategy:
       * 1. Start with minimal nodes (single cubic segment)
       * 2. Gradually increase nodes if RSS > s (underfitting)
       * 3. Find the configuration where RSS ≈ s
       */
      void fit_smoothing_spline(const std::vector<double>& x,
                                const std::vector<double>& y,
                                double s_target);

      /**
       * @brief Try polynomial fit (matches scipy when 0 interior knots)
       * 
       * @return true if polynomial fit is acceptable (RSS <= s_target)
       */
      bool try_polynomial_fit(const std::vector<double>& x,
                              const std::vector<double>& y,
                              double s_target,
                              int degree = 3);

      /**
       * @brief Compute RSS for given spline
       */
      double compute_rss(const BSpline2d* spl,
                        const std::vector<double>& x,
                        const std::vector<double>& y) const;
      
      /**
       * @brief Compute RSS for polynomial fit
       */
      double compute_polynomial_rss(const std::vector<double>& x,
                                    const std::vector<double>& y) const;
      
      /**
       * @brief Evaluate polynomial at x
       */
      double eval_polynomial(double x) const;
    };

  } // namespace Math
} // namespace OpenMS
