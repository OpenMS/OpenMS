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
#include <bit>
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
      std::sort(xf.begin(), xf.end());
      double q25 = Math::quantile(xf.begin(), xf.end(), 0.25);
      double q75 = Math::quantile(xf.begin(), xf.end(), 0.75);
      double iqr = q75 - q25;

      double lo = std::min(sd, iqr / 1.34);
      // fallbacks to match pyprophet/numpy implementation
      if (!(lo > 0.0))
      {
        if (sd > 0.0) lo = sd;
        else if (!xf.empty()) lo = std::abs(xf[0]);
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

      // create a 1-D tensor with zero-padding to size M
      evergreen::Tensor<double> t(evergreen::Vector<unsigned long>{static_cast<unsigned long>(M)});
      for (std::size_t i = 0; i < M; ++i)
      {
        double v = 0.0;
        if (i < X.size() && std::isfinite(X[i])) v = X[i];
        t[i] = v;
      }

      // Compute real FFT: input size M produces (M/2 + 1) complex values
      // due to conjugate symmetry of real FFT output
      evergreen::Tensor<evergreen::cpx> fft_result = evergreen::real_fft<evergreen::DIF, false, false, true>(t);

      // Convert to Munro-packed format for efficient storage and pointwise multiplication:
      // [real_0, real_1, ..., real_M/2, imag_1, ..., imag_(M/2-1)]
      // This avoids storing redundant imaginary parts and reduces memory overhead
      const std::size_t half = M / 2;
      const std::size_t ny = fft_result.flat_size(); // should be half + 1
      std::vector<double> out(M, 0.0);

      // store real parts: real_0, real_1, ..., real_M/2
      for (std::size_t k = 0; k < ny; ++k) out[k] = fft_result[k].r;
      // store imaginary parts: imag_1, ..., imag_(M/2-1) in second half
      for (std::size_t k = 1; k < ny - 1; ++k) out[half + k] = fft_result[k].i;

      return out;
    }

    std::vector<double> revRt(const std::vector<double>& Xp, std::size_t M)
    {
      if (M == 0) M = Xp.size();
      if (M == 0) return std::vector<double>();

      const std::size_t half = M / 2;
      const std::size_t ny = half + 1;

      if (Xp.size() < M) throw std::invalid_argument("revRt: input length must equal M");

      // Unpack Munro-packed format: [real_0, real_1, ..., real_M/2, imag_1, ..., imag_(M/2-1)]
      // into complex tensor for inverse FFT
      evergreen::Tensor<evergreen::cpx> packed(evergreen::Vector<unsigned long>{static_cast<unsigned long>(ny)});
      for (std::size_t k = 0; k < ny; ++k)
      {
        double re = Xp[k];
        double im = 0.0;
        if (k > 0 && k < ny - 1)
        {
          im = Xp[half + k];  // imaginary parts stored in second half of input vector
        }
        packed[k] = evergreen::cpx{re, im};
      }

      // Compute inverse real FFT: (M/2 + 1) complex values produce M real values
      evergreen::Tensor<double> rec = evergreen::real_ifft<evergreen::DIF, false, false>(packed);

      std::vector<double> out(M, 0.0);
      const std::size_t ncopy = std::min(static_cast<std::size_t>(rec.flat_size()), static_cast<std::size_t>(M));
      for (std::size_t i = 0; i < ncopy; ++i) out[i] = rec[i];
      return out;
    }

    std::vector<double> silvermanKernelFFT(double bw, std::size_t M, double RANGE)
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

    std::pair<std::vector<double>, std::vector<double>> gridKdeFFT(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
    {
      const std::size_t n = x.size();

      // M = 2^ceil(log2(max(gridsize, n, 512)))
      std::size_t target = std::max<std::size_t>(gridsize, std::max<std::size_t>(n, 512));
      std::size_t M = std::bit_ceil(target);

      // compute a,b
      double a = 0.0, b = 0.0;
      if (n > 0)
      {
        auto mm = std::minmax_element(x.begin(), x.end());
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
      std::vector<double> binned = fast_linbin(x, a, b, M);
      if (n > 0)
      {
        for (double &v : binned) v /= (delta * static_cast<double>(n));
      }

      // FFT-based convolution using Silverman's algorithm (see gridKdeFFT() documentation)
      std::vector<double> Y = forRt(binned, M);

      // multiply by Silverman kernel FFT
      std::vector<double> K = silvermanKernelFFT(bw, M, RANGE);
      if (K.size() != Y.size()) throw std::runtime_error("gridKdeFFT: internal size mismatch");
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

    std::vector<double> kdeFFTEval(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
    {
      // build grid density
      auto pg = gridKdeFFT(x, bw, gridsize, cut);
      const std::vector<double>& dens = pg.first;
      const std::vector<double>& grid = pg.second;

      // build cubic spline (CubicSpline2d expects sorted x)
      OpenMS::CubicSpline2d spline(grid, dens);

      std::vector<double> out;
      out.reserve(x.size());
      for (double xi : x)
      {
        double val = spline.eval(xi);
        // Clamp negative densities to zero (can occur from extrapolation outside grid range)
        out.push_back(val < 0.0 ? 0.0 : val);
      }
      return out;
    }

  } // namespace Math
} // namespace OpenMS
