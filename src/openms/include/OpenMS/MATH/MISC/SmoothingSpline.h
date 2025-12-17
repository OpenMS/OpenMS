// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
#pragma once

#include <vector>
#include <limits>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

namespace OpenMS
{
  namespace Math
  {
    /// Simple smoothing-spline utility (discrete penalty approximation).
    ///
    /// This class implements a small, robust smoothing spline fitter using a
    /// discrete second-difference penalty. It is not a full replacement for
    /// R's smooth.spline or SciPy's sophisticated handler, but provides a
    /// deterministic smoothing operator that is sufficient for reproducing
    /// the pi0 smoothing behavior used in the tests.
    class SmoothingSpline
    {
    public:
      /// Construct and fit a smoothing spline to (xs, ys).
      /// If lambda <= 0 the class will choose a smoothing parameter using
      /// generalized cross-validation (GCV) on a log-grid search.
      SmoothingSpline(const std::vector<double>& xs,
                      const std::vector<double>& ys,
                      double lambda = -1.0);

      /// Evaluate the fitted smoothing spline at x. If construction failed
      /// the function returns NaN.
      double eval(double x) const;

      /// Return the smoothed values at the original knots (xs order).
      const std::vector<double>& smoothed_y() const { return smoothed_y_; }

      /// Whether the fit succeeded
      bool ok() const { return ok_; }

    private:
      bool ok_ = false;
      std::vector<double> xs_;
      std::vector<double> smoothed_y_;
      CubicSpline2d spline_; // built from xs_ and smoothed_y_
    };
  }
}
