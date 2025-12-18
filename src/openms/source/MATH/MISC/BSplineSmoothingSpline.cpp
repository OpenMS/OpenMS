// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/BSplineSmoothingSpline.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <cmath>
#include <limits>
#include <algorithm>

namespace OpenMS
{
  namespace Math
  {

    BSplineSmoothingSpline::BSplineSmoothingSpline(const std::vector<double>& x,
                                                   const std::vector<double>& y,
                                                   double s,
                                                   int k)
      : k_(k), x_(x), y_(y)
    {
      ok_ = false;
      spline_ = nullptr;

      const int n = static_cast<int>(x.size());
      if (n < 2 || static_cast<int>(y.size()) != n)
      {
        OPENMS_LOG_ERROR << "BSplineSmoothingSpline: Invalid input sizes" << std::endl;
        return;
      }

      // Check that x is sorted
      for (size_t i = 1; i < x.size(); ++i)
      {
        if (x[i] <= x[i-1])
        {
          OPENMS_LOG_ERROR << "BSplineSmoothingSpline: x values must be sorted" << std::endl;
          return;
        }
      }

      // Compute automatic smoothing parameter if needed
      if (s < 0.0)
      {
        // Default from scipy: s = m - sqrt(2*m)
        double m = static_cast<double>(n);
        s_ = m - std::sqrt(2.0 * m);
      }
      else
      {
        s_ = s;
      }

      // For very small n or s=0, use interpolation
      if (n < 4 || s_ <= 0.0)
      {
        // Pure interpolation: use BSpline with wavelength=0, auto nodes
        try
        {
          spline_ = new BSpline2d(x, y, 0.0, BSpline2d::BC_ZERO_SECOND, 0);
          if (spline_->ok())
          {
            ok_ = true;
            rss_ = compute_rss(spline_, x, y);
            num_interior_knots_ = n - 2; // Interpolation uses all points
          }
        }
        catch (...)
        {
          if (spline_) delete spline_;
          spline_ = nullptr;
        }
        return;
      }

      // Fit smoothing spline with target RSS ≈ s_
      fit_smoothing_spline(x, y, s_);
    }

    BSplineSmoothingSpline::~BSplineSmoothingSpline()
    {
      delete spline_;
    }

    bool BSplineSmoothingSpline::try_polynomial_fit(const std::vector<double>& x,
                                                      const std::vector<double>& y,
                                                      double s_target,
                                                      int degree)
    {
      const int n = static_cast<int>(x.size());
      const int m = degree + 1;  // Number of coefficients
      
      if (n < m)
      {
        return false;  // Not enough points
      }

      // Solve least squares: X^T X c = X^T y
      // where X is the Vandermonde matrix
      
      // Build X^T X (m x m symmetric matrix)
      std::vector<double> XtX(m * m, 0.0);
      std::vector<double> Xty(m, 0.0);
      
      for (int i = 0; i < n; ++i)
      {
        double xi = x[i];
        double yi = y[i];
        
        // Powers of xi
        std::vector<double> powers(m);
        powers[0] = 1.0;
        for (int p = 1; p < m; ++p)
        {
          powers[p] = powers[p-1] * xi;
        }
        
        // Accumulate X^T X and X^T y
        for (int r = 0; r < m; ++r)
        {
          Xty[r] += powers[r] * yi;
          for (int c = 0; c < m; ++c)
          {
            XtX[r * m + c] += powers[r] * powers[c];
          }
        }
      }
      
      // Solve using Gaussian elimination
      std::vector<double> coeffs(m);
      
      // Augment XtX with Xty
      for (int i = 0; i < m; ++i)
      {
        // Partial pivoting
        int pivot = i;
        double maxval = std::abs(XtX[i * m + i]);
        for (int j = i + 1; j < m; ++j)
        {
          double val = std::abs(XtX[j * m + i]);
          if (val > maxval)
          {
            maxval = val;
            pivot = j;
          }
        }
        
        if (maxval < 1e-10)
        {
          return false;  // Singular matrix
        }
        
        // Swap rows
        if (pivot != i)
        {
          for (int c = 0; c < m; ++c)
          {
            std::swap(XtX[i * m + c], XtX[pivot * m + c]);
          }
          std::swap(Xty[i], Xty[pivot]);
        }
        
        // Eliminate
        double diag = XtX[i * m + i];
        for (int j = i + 1; j < m; ++j)
        {
          double factor = XtX[j * m + i] / diag;
          for (int c = i; c < m; ++c)
          {
            XtX[j * m + c] -= factor * XtX[i * m + c];
          }
          Xty[j] -= factor * Xty[i];
        }
      }
      
      // Back substitution
      for (int i = m - 1; i >= 0; --i)
      {
        double sum = Xty[i];
        for (int j = i + 1; j < m; ++j)
        {
          sum -= XtX[i * m + j] * coeffs[j];
        }
        coeffs[i] = sum / XtX[i * m + i];
      }
      
      // Store coefficients
      poly_coeffs_ = coeffs;
      fit_type_ = POLYNOMIAL;
      
      // Compute RSS
      double rss = compute_polynomial_rss(x, y);
      
      // Check for debug output
      const char* dbg_env = std::getenv("OPENMS_TEST_VERBOSE");
      const char* dbg_pi0 = std::getenv("OPENMS_PI0_BSPLINE_DEBUG");
      bool debug = (dbg_env && dbg_env[0] != '\0') || (dbg_pi0 && dbg_pi0[0] != '\0');
      
      if (debug)
      {
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: Polynomial fit (degree=" << degree << ")" << std::endl;
        OPENMS_LOG_INFO << "  Coefficients: [";
        for (size_t i = 0; i < coeffs.size(); ++i)
        {
          if (i > 0) OPENMS_LOG_INFO << ", ";
          OPENMS_LOG_INFO << coeffs[i];
        }
        OPENMS_LOG_INFO << "]" << std::endl;
        OPENMS_LOG_INFO << "  RSS: " << rss << " (target: " << s_target << ")" << std::endl;
      }
      
      // Accept if RSS <= s_target (or close enough)
      return rss <= s_target * 1.1;  // Allow 10% margin
    }

    double BSplineSmoothingSpline::compute_polynomial_rss(const std::vector<double>& x,
                                                           const std::vector<double>& y) const
    {
      double rss = 0.0;
      for (size_t i = 0; i < x.size(); ++i)
      {
        double fitted = eval_polynomial(x[i]);
        double residual = y[i] - fitted;
        rss += residual * residual;
      }
      return rss;
    }

    double BSplineSmoothingSpline::eval_polynomial(double x) const
    {
      double result = 0.0;
      double xpow = 1.0;
      for (size_t i = 0; i < poly_coeffs_.size(); ++i)
      {
        result += poly_coeffs_[i] * xpow;
        xpow *= x;
      }
      return result;
    }

    void BSplineSmoothingSpline::fit_smoothing_spline(const std::vector<double>& x,
                                                       const std::vector<double>& y,
                                                       double s_target)
    {
      const int n = static_cast<int>(x.size());
      
      // Check for debug output
      const char* dbg_env = std::getenv("OPENMS_TEST_VERBOSE");
      const char* dbg_pi0 = std::getenv("OPENMS_PI0_BSPLINE_DEBUG");
      bool debug = (dbg_env && dbg_env[0] != '\0') || (dbg_pi0 && dbg_pi0[0] != '\0');

      if (debug)
      {
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: Fitting with s_target=" << s_target 
                        << " for n=" << n << " points" << std::endl;
      }

      // KEY INSIGHT from scipy analysis:
      // UnivariateSpline with typical smoothing (s = m - sqrt(2m)) often selects
      // ZERO interior knots, which is mathematically equivalent to a cubic polynomial.
      // 
      // For example, with m=12 and s≈7.1, scipy creates a single cubic segment
      // with knots [x_min, x_min, x_min, x_min, x_max, x_max, x_max, x_max]
      // This is exactly a cubic polynomial fit!
      //
      // Since BSpline2d uses different basis functions than scipy's FITPACK,
      // we should first try a simple polynomial fit, which will exactly match
      // scipy when it chooses 0 interior knots.

      // STEP 1: Try polynomial fit (matches scipy exactly for 0 interior knots)
      if (try_polynomial_fit(x, y, s_target, 3))  // Cubic polynomial
      {
        // Polynomial fit succeeded and RSS <= s_target
        rss_ = compute_polynomial_rss(x, y);
        num_interior_knots_ = 0;  // No interior knots (single polynomial segment)
        ok_ = true;
        
        if (debug)
        {
          OPENMS_LOG_INFO << "BSplineSmoothingSpline: Using polynomial fit" << std::endl;
          OPENMS_LOG_INFO << "  RSS=" << rss_ << " for target s=" << s_target << std::endl;
          double x_max = x.back();
          double y_max = eval_polynomial(x_max);
          OPENMS_LOG_INFO << "  eval(" << x_max << ") = " << y_max << std::endl;
        }
        return;
      }
      
      // STEP 2: Polynomial didn't fit well enough, try BSpline with more knots
      if (debug)
      {
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: Polynomial fit insufficient, trying BSpline with interior knots" << std::endl;
      }

      // Strategy: Try different numbers of interior knots
      // Start from 0 (single cubic segment) and increase until RSS drops below s_target
      // 
      // Key insight from scipy analysis:
      // - For typical pi0 data with s≈7, scipy chooses 0 interior knots (single cubic)
      // - This gives RSS ≈ 0.0024, which is less than s=7.1
      // - The smoothing comes from using minimal DOF, not from penalties
      
      struct Candidate
      {
        int num_nodes;
        int num_interior_knots;
        double rss;
        BSpline2d* spline;
        bool valid;
      };

      std::vector<Candidate> candidates;

      // Try node counts from minimal (single segment) to full interpolation
      // For cubic splines (k=3), we need at least 4 coefficients → num_nodes >= 4
      // num_nodes relates to interior knots as: num_interior = num_nodes - 2
      
      std::vector<int> node_counts;
      // Start with very few nodes (high smoothing)
      node_counts.push_back(4);  // 2 interior knots → nearly polynomial
      node_counts.push_back(6);  // 4 interior knots
      node_counts.push_back(8);  // 6 interior knots
      node_counts.push_back(std::max(4, n/2));  // medium
      node_counts.push_back(std::max(4, 3*n/4));  // higher
      node_counts.push_back(n);  // full interpolation

      // Remove duplicates and sort
      std::sort(node_counts.begin(), node_counts.end());
      node_counts.erase(std::unique(node_counts.begin(), node_counts.end()), node_counts.end());

      if (debug)
      {
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: Testing node counts: [";
        for (size_t i = 0; i < node_counts.size(); ++i)
        {
          if (i > 0) OPENMS_LOG_INFO << ", ";
          OPENMS_LOG_INFO << node_counts[i];
        }
        OPENMS_LOG_INFO << "]" << std::endl;
      }

      for (int num_nodes : node_counts)
      {
        try
        {
          // Try with wavelength=0 (no smoothing penalty in BSpline optimization)
          BSpline2d* spl = new BSpline2d(x, y, 0.0, BSpline2d::BC_ZERO_SECOND, num_nodes);
          
          if (!spl->ok())
          {
            delete spl;
            continue;
          }

          double rss = compute_rss(spl, x, y);
          int num_interior = num_nodes - 2;

          candidates.push_back(Candidate{num_nodes, num_interior, rss, spl, true});

          if (debug)
          {
            OPENMS_LOG_INFO << "BSplineSmoothingSpline: num_nodes=" << num_nodes
                            << " (interior=" << num_interior << ") → RSS=" << rss 
                            << " (target=" << s_target << ")" << std::endl;
          }
        }
        catch (...)
        {
          // Skip this configuration
          continue;
        }
      }

      if (candidates.empty())
      {
        OPENMS_LOG_ERROR << "BSplineSmoothingSpline: No valid spline configurations found" << std::endl;
        return;
      }

      // Select best candidate:
      // 1. Prefer RSS closest to s_target (but ideally RSS < s_target)
      // 2. Among similar RSS values, prefer simpler model (fewer knots)
      
      std::sort(candidates.begin(), candidates.end(), 
                [s_target](const Candidate& a, const Candidate& b) {
                  // Primary criterion: closeness to s_target
                  double err_a = std::abs(a.rss - s_target);
                  double err_b = std::abs(b.rss - s_target);
                  
                  // If both are below target and close, prefer simpler
                  if (a.rss <= s_target && b.rss <= s_target && std::abs(err_a - err_b) < 0.001)
                  {
                    return a.num_interior_knots < b.num_interior_knots;
                  }
                  
                  // Otherwise prefer closest to target
                  return err_a < err_b;
                });

      // Use the best candidate
      const Candidate& best = candidates[0];
      spline_ = best.spline;
      rss_ = best.rss;
      num_interior_knots_ = best.num_interior_knots;
      ok_ = true;

      // Clean up other candidates
      for (size_t i = 1; i < candidates.size(); ++i)
      {
        delete candidates[i].spline;
      }

      if (debug)
      {
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: Selected num_nodes=" << best.num_nodes
                        << " (interior=" << num_interior_knots_ << ") with RSS=" << rss_
                        << " for target s=" << s_target << std::endl;
        
        // Evaluate at some test points
        double x_max = x.back();
        double y_max = eval(x_max);
        OPENMS_LOG_INFO << "BSplineSmoothingSpline: eval(" << x_max << ") = " << y_max << std::endl;
      }
    }

    double BSplineSmoothingSpline::compute_rss(const BSpline2d* spl,
                                                const std::vector<double>& x,
                                                const std::vector<double>& y) const
    {
      if (!spl || !spl->ok()) return std::numeric_limits<double>::infinity();

      double rss = 0.0;
      for (size_t i = 0; i < x.size(); ++i)
      {
        double fitted = spl->eval(x[i]);
        double residual = y[i] - fitted;
        rss += residual * residual;
      }
      return rss;
    }

    double BSplineSmoothingSpline::eval(double x) const
    {
      if (!ok_) return std::numeric_limits<double>::quiet_NaN();
      
      if (fit_type_ == POLYNOMIAL)
      {
        return eval_polynomial(x);
      }
      else  // BSPLINE
      {
        if (!spline_) return std::numeric_limits<double>::quiet_NaN();
        return spline_->eval(x);
      }
    }

  } // namespace Math
} // namespace OpenMS
