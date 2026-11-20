// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
//

#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace std;

namespace OpenMS::Math
{
    namespace
    {
      /// Weighted sums needed by the Gumbel score equations for a given scale b.
      /// The data are shifted by their minimum (s = x - x_min >= 0), so that exp(-s/b) <= 1 cannot overflow.
      struct GumbelScaleSums
      {
        double sum_e = 0.0;   ///< sum w * exp(-s/b)
        double sum_es = 0.0;  ///< sum w * s * exp(-s/b)
        double sum_es2 = 0.0; ///< sum w * s^2 * exp(-s/b)
      };

      GumbelScaleSums gumbelScaleSums(const std::vector<double>& s, const std::vector<double>& w, double b)
      {
        GumbelScaleSums r;
        for (size_t i = 0; i < s.size(); ++i)
        {
          const double e = w[i] * std::exp(-s[i] / b);
          r.sum_e += e;
          r.sum_es += e * s[i];
          r.sum_es2 += e * s[i] * s[i];
        }
        return r;
      }
    }

    GumbelMaxLikelihoodFitter::GumbelDistributionFitResult GumbelMaxLikelihoodFitter::fitWeighted(const std::vector<double> & x, const std::vector<double> & w)
    {
      if (x.size() != w.size())
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Number of data points and weights differ.");
      }

      // keep only finite data points with positive, finite weight
      std::vector<double> xs, ws;
      xs.reserve(x.size());
      ws.reserve(w.size());
      for (size_t i = 0; i < x.size(); ++i)
      {
        if (std::isfinite(x[i]) && std::isfinite(w[i]) && w[i] > 0.0)
        {
          xs.push_back(x[i]);
          ws.push_back(w[i]);
        }
      }

      // nothing to fit: return the initial parameters (documented behaviour)
      if (xs.empty())
      {
        return init_param_;
      }

      // Weighted maximum likelihood for the Gumbel (maximum) distribution
      //   f(x) = 1/b * exp(-(z + exp(-z))),  z = (x - a) / b.
      // The score equations reduce to a one-dimensional equation for the scale b:
      //   g(b) = b - mean_w(x) + sum(w x e^{-x/b}) / sum(w e^{-x/b}) = 0,
      // and, given b, the location follows in closed form:
      //   a = -b * log( sum(w e^{-x/b}) / sum(w) ).
      // g is strictly increasing in b (g'(b) = 1 + Var_q(x) / b^2 with q ~ w e^{-x/b}), so the root is unique.
      const double x_min = *std::min_element(xs.begin(), xs.end());
      double sum_w = 0.0, sum_ws = 0.0;
      std::vector<double> s(xs.size());
      for (size_t i = 0; i < xs.size(); ++i)
      {
        s[i] = xs[i] - x_min;
        sum_w += ws[i];
        sum_ws += ws[i] * s[i];
      }
      const double mean_s = sum_ws / sum_w;
      double var_s = 0.0;
      for (size_t i = 0; i < s.size(); ++i)
      {
        var_s += ws[i] * (s[i] - mean_s) * (s[i] - mean_s);
      }
      var_s /= sum_w;

      if (!(var_s > 0.0))
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GumbelMaxLikelihoodFitter",
                                     "Could not fit the Gumbel distribution: the (weighted) data have no spread.");
      }

      auto g = [&](double b, double* derivative) -> double
      {
        GumbelScaleSums r = gumbelScaleSums(s, ws, b);
        const double m = r.sum_es / r.sum_e;
        if (derivative != nullptr)
        {
          const double v = std::max(0.0, r.sum_es2 / r.sum_e - m * m);
          *derivative = 1.0 + v / (b * b);
        }
        return b - mean_s + m;
      };

      // Start from the initial (e.g. previous) scale parameter if available (warm start, e.g. in EM iterations),
      // otherwise from the moment estimate b0 = sqrt(6)/pi * sd, and bracket the root.
      // g(0+) = min(s) - mean(s) < 0 because the data have spread, and g(b) -> infinity for large b.
      const double b_moment = std::sqrt(6.0 * var_s) / Constants::PI;
      double lo = 0.0, hi = 0.0;
      // brackets the root by doubling/halving from b0; false if no sign change was found within 200 steps
      auto bracket = [&](double b0) -> bool
      {
        lo = b0;
        hi = b0;
        if (g(b0, nullptr) < 0.0)
        {
          for (int i = 0; i < 200; ++i)
          {
            lo = hi;
            hi *= 2.0;
            if (g(hi, nullptr) > 0.0) { return true; }
          }
          return false;
        }
        for (int i = 0; i < 200; ++i)
        {
          hi = lo;
          lo *= 0.5;
          if (g(lo, nullptr) < 0.0) { return true; }
        }
        return false;
      };
      double b_start = (std::isfinite(init_param_.b) && init_param_.b > 0.0) ? init_param_.b : b_moment;
      if (!bracket(b_start))
      {
        // a warm start far away from the solution: start again from the moment estimate
        b_start = b_moment;
        if (!bracket(b_start))
        {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GumbelMaxLikelihoodFitter",
                                       "Could not bracket the maximum likelihood estimate of the Gumbel scale parameter.");
        }
      }

      // safeguarded Newton iterations (fall back to bisection whenever Newton leaves the bracket)
      double b = std::clamp(b_start, lo, hi);
      for (int iteration = 0; iteration < 200; ++iteration)
      {
        double dg = 0.0;
        const double gb = g(b, &dg);
        if (gb > 0.0)
        {
          hi = b;
        }
        else
        {
          lo = b;
        }
        double b_new = b - gb / dg;
        if (!(b_new > lo && b_new < hi))
        {
          b_new = 0.5 * (lo + hi);
        }
        if (std::fabs(b_new - b) <= 1e-12 * std::max(1.0, b))
        {
          b = b_new;
          break;
        }
        b = b_new;
      }

      const GumbelScaleSums r = gumbelScaleSums(s, ws, b);
      const double a = x_min - b * std::log(r.sum_e / sum_w);

      if (!std::isfinite(a) || !std::isfinite(b) || b <= 0.0)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GumbelMaxLikelihoodFitter",
                                     "Could not fit the Gumbel distribution to the data.");
      }

      init_param_.a = a;
      init_param_.b = b;
      return {a, b};
    }

    double GumbelMaxLikelihoodFitter::GumbelDistributionFitResult::log_eval_no_normalize(const double x) const
    {
      // note: despite its name, this is the fully normalized log density
      double diff = (x - a)/b;
      return -log(b) - diff - exp(- diff);
    }

    double GumbelMaxLikelihoodFitter::GumbelDistributionFitResult::eval(const double x) const
    {
      return std::exp(log_eval_no_normalize(x));
    }

    double GumbelMaxLikelihoodFitter::GumbelDistributionFitResult::mode() const
    {
      return a;
    }

    GumbelMaxLikelihoodFitter::GumbelMaxLikelihoodFitter():
    init_param_({0.25, 0.1})
    {}

    GumbelMaxLikelihoodFitter::GumbelMaxLikelihoodFitter(GumbelMaxLikelihoodFitter::GumbelDistributionFitResult init):
        init_param_(init)
    {}

    GumbelMaxLikelihoodFitter::~GumbelMaxLikelihoodFitter() = default;

    void GumbelMaxLikelihoodFitter::setInitialParameters(const GumbelDistributionFitResult & param)
    {
      init_param_ = param;
    }

} // namespace OpenMS   //namespace Math
