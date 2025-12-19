// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

 // Use the evergreen umbrella header so the vendor headers are included
 // consistently inside the `evergreen` namespace and their helper
 // macros are defined. Then include the FFT submodule *inside* the
 // `evergreen` namespace so its template and type definitions live in
 // `evergreen::` as the rest of the code expects.
#include <Evergreen/evergreen.hpp>
namespace evergreen {
#include <FFT/FFT.hpp>
}

#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/STATISTICS/KernelDensityEstimation.h>
#include <OpenMS/MATH/StatisticFunctions.h>

namespace OpenMS
{
  namespace Math
  {

    // Helper: next power of two >= v
    static std::size_t next_pow2(std::size_t v)
    {
      if (v == 0) return 1;
      std::size_t p = 1;
      while (p < v) p <<= 1;
      return p;
    }

    // Fast linear binning matching statsmodels' endpoint rules (internal helper)
    static std::vector<double> fast_linbin(const std::vector<double>& x, double a, double b, std::size_t M)
    {
      // linear binning with inclusive a, inclusive b into M bins (np.linspace(a,b,M))
      std::vector<double> bins(M, 0.0);
      if (x.empty() || !(b > a)) return bins;
      const double delta = (b - a) / static_cast<double>(M - 1); // linspace with M points has M-1 intervals
      for (double v : x)
      {
        if (!std::isfinite(v)) continue;
        if (v < a || v > b) continue;
        double pos = (v - a) / delta;
        std::size_t left = static_cast<std::size_t>(std::floor(pos));
        if (left >= M - 1)
        {
          bins[M - 1] += 1.0;
          continue;
        }
        double frac = pos - static_cast<double>(left);
        bins[left] += 1.0 - frac;
        bins[left + 1] += frac;
      }
      return bins;
    }

    double bwNrd0(const std::vector<double>& x)
    {
      // filter finite
      std::vector<double> xf;
      xf.reserve(x.size());
      for (double v : x) if (std::isfinite(v)) xf.push_back(v);
      const std::size_t n = xf.size();
      if (n < 2) return 0.0;

      // sample standard deviation (ddof=1)
      double mean = 0.0;
      for (double v : xf) mean += v;
      mean /= static_cast<double>(n);
      double var = 0.0;
      for (double v : xf) { double d = v - mean; var += d * d; }
      var /= static_cast<double>(n - 1);
      double sd = std::sqrt(var);

      // compute percentiles using existing Math::quantile (Type 7, numpy-like linear interpolation)
      std::vector<double> cp = xf;
      std::sort(cp.begin(), cp.end());
      double q25 = Math::quantile(cp.begin(), cp.end(), 0.25);
      double q75 = Math::quantile(cp.begin(), cp.end(), 0.75);
      double iqr = q75 - q25;

      double lo = std::min(sd, iqr / 1.34);
      // fallbacks to match pyprophet/numpy implementation
      if (!(lo > 0.0))
      {
        if (sd > 0.0) lo = sd;
        else if (!cp.empty()) lo = std::abs(cp[0]);
        else lo = 1.0;
      }

      double bw = 0.9 * lo * std::pow(static_cast<double>(n), -0.2); // Silverman variant used by statsmodels/pyprophet
      return bw;
    }

    std::vector<double> linBin(const std::vector<double>& x, double xmin, double xmax, std::size_t nbins, const std::vector<double>* weights)
    {
      if (nbins == 0) throw std::invalid_argument("linBin: nbins must be > 0");
      std::vector<double> bins(nbins, 0.0);
      if (x.empty()) return bins;
      if (!(xmax > xmin)) throw std::invalid_argument("linBin: xmax must be > xmin");

      const double width = (xmax - xmin) / static_cast<double>(nbins);
      if (!(width > 0.0)) return bins;

      const bool use_weights = (weights != nullptr && weights->size() == x.size());

      for (std::size_t i = 0; i < x.size(); ++i)
      {
        double v = x[i];
        if (!std::isfinite(v)) continue;
        if (v < xmin || v > xmax) continue; // ignore outside range; include xmax in last bin
        std::size_t idx = static_cast<std::size_t>(std::floor((v - xmin) / width));
        if (idx >= nbins) idx = nbins - 1; // edge case when v == xmax
        double w = 1.0;
        if (use_weights) w = (*weights)[i];
        bins[idx] += w;
      }

      return bins;
    }

    std::vector<double> forRt(const std::vector<double>& X, std::size_t M)
    {
      if (M == 0) M = X.size();
      if (M == 0) return std::vector<double>();

      const std::size_t half = M / 2;
      const std::size_t out_len = M;

      // create a 1-D tensor and copy (zero-pad if needed)
      evergreen::Tensor<double> t({M});
      for (std::size_t i = 0; i < M; ++i)
      {
        double v = 0.0;
        if (i < X.size() && std::isfinite(X[i])) v = X[i];
        t[i] = v;
      }

      // use DIF-based real FFT provided by evergreen (packed real format)
      // Use FRESHLY_ZERO_PADDED=true for the single-axis case to match optimized path
      evergreen::Tensor<evergreen::cpx> packed = evergreen::real_fft<evergreen::DIF, false, false, true>(t);

      const std::size_t ny = packed.flat_size(); // should be half + 1
      std::vector<double> out(out_len, 0.0);

      // store real parts
      for (std::size_t k = 0; k < ny; ++k) out[k] = packed[k].r;
      // store imag parts for k=1..half-1
      for (std::size_t k = 1; k < ny - 1; ++k) out[half + k] = packed[k].i;

      return out;
    }

    std::vector<double> revRt(const std::vector<double>& Xp, std::size_t M)
    {
      if (M == 0) M = Xp.size();
      if (M == 0) return std::vector<double>();

      const std::size_t half = M / 2;
      const std::size_t ny = half + 1;

      if (Xp.size() < M) throw std::invalid_argument("revRt: input length must equal M");

      // build packed complex tensor
      evergreen::Tensor<evergreen::cpx> packed({ny});
      for (std::size_t k = 0; k < ny; ++k)
      {
        double re = Xp[k];
        double im = 0.0;
        if (k > 0 && k < ny - 1)
        {
          im = Xp[half + k];
        }
        packed[k] = evergreen::cpx{re, im};
      }

      // use evergreen real_ifft which returns a real tensor of length M
      evergreen::Tensor<double> rec = evergreen::real_ifft<evergreen::DIF, false, false>(packed);

      std::vector<double> out(M, 0.0);
      const std::size_t ncopy = std::min(rec.flat_size(), M);
      for (std::size_t i = 0; i < ncopy; ++i) out[i] = rec[i];
      return out;
    }

    std::vector<double> silvermanKernelFft(double bw, std::size_t M, double RANGE)
    {
      if (M == 0) return std::vector<double>();
      const std::size_t half = M / 2;
      // J = 0..M/2
      std::vector<double> J(half + 1);
      for (std::size_t k = 0; k <= half; ++k) J[k] = static_cast<double>(k);
      const double FAC1 = 2.0 * std::pow(M_PI * bw / RANGE, 2);
      std::vector<double> FAC(half + 1);
      for (std::size_t k = 0; k <= half; ++k)
      {
        double j = J[k];
        double jfac = j * j * FAC1;
        double BC = 1.0 - (j * (M_PI / static_cast<double>(M))) * (j * (M_PI / static_cast<double>(M))) / 3.0;
        if (!(BC > 0.0)) BC = std::numeric_limits<double>::min();
        FAC[k] = std::exp(-jfac) / BC;
      }

      // mirror to Munro-packed length M: [FAC[0..half], FAC[1..half-1]]
      std::vector<double> out(M, 0.0);
      // copy FAC[0..half]
      for (std::size_t k = 0; k <= half; ++k) out[k] = FAC[k];
      // mirror FAC[1..half-1]
      if (half >= 2)
      {
        for (std::size_t k = 1; k < half; ++k) out[half + k] = FAC[k];
      }
      return out;
    }

    std::pair<std::vector<double>, std::vector<double>> gridKdeFft(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
    {
      std::vector<double> xx = x;
      // ensure column vector semantics
      // n
      const std::size_t n = xx.size();

      // M = 2^ceil(log2(max(gridsize, n, 512)))
      std::size_t target = std::max<std::size_t>(gridsize, std::max<std::size_t>(n, 512));
      std::size_t M = next_pow2(target);

      // compute a,b
      double a = 0.0, b = 0.0;
      if (n > 0)
      {
        auto mm = std::minmax_element(xx.begin(), xx.end());
        a = *mm.first - cut * bw;
        b = *mm.second + cut * bw;
      }
      else
      {
        a = -cut * bw;
        b = cut * bw;
      }
      // build grid (M points linspace a..b)
      std::vector<double> grid(M);
      double delta = (b - a) / static_cast<double>(M - 1);
      for (std::size_t i = 0; i < M; ++i) grid[i] = a + static_cast<double>(i) * delta;

      double RANGE = b - a;

      // linear binning normalized by (delta * n)
      std::vector<double> binned = fast_linbin(xx, a, b, M);
      if (n > 0)
      {
        for (double &v : binned) v /= (delta * static_cast<double>(n));
      }

      // forward Munro RFFT
      std::vector<double> Y = forRt(binned, M);

      // multiply by Silverman kernel FFT
      std::vector<double> K = silvermanKernelFft(bw, M, RANGE);
      if (K.size() != Y.size()) throw std::runtime_error("gridKdeFft: internal size mismatch");
      std::vector<double> Z(M);
      for (std::size_t i = 0; i < M; ++i) Z[i] = K[i] * Y[i];

      // inverse
      std::vector<double> dens = revRt(Z, M);

      // renormalize to enforce integral 1
      double sum = 0.0;
      for (double v : dens) sum += v;
      if (sum > 0.0 && delta > 0.0)
      {
        double scale = 1.0 / (sum * delta);
        for (double &v : dens) v *= scale;
      }

      return {dens, grid};
    }

    std::vector<double> kdeFftEval(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
    {
      // build grid density
      auto pg = gridKdeFft(x, bw, gridsize, cut);
      const std::vector<double>& dens = pg.first;
      const std::vector<double>& grid = pg.second;

      // build cubic spline (CubicSpline2d expects sorted x)
      OpenMS::CubicSpline2d spline(grid, dens);

      std::vector<double> out;
      out.reserve(x.size());
      for (double xi : x)
      {
        out.push_back(spline.eval(xi));
      }
      return out;
    }

  } // namespace Math
} // namespace OpenMS
