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

#include <cstddef>
#include <utility>
#include <vector>

#include <OpenMS/config.h>

namespace OpenMS
{
  namespace Math
  {
    /**
      @brief Kernel Density Estimation utilities using FFT-based methods.

      This class provides efficient kernel density estimation using Fast Fourier
      Transform (FFT) based convolution, following the Silverman (1982) algorithm
      as implemented in Python's statsmodels package.

      The FFT-based approach scales as O(n + M*log(M)) compared to O(n*M) for
      direct evaluation, making it highly efficient for large datasets or fine grids.

      Key features:
      - Automatic bandwidth selection using Silverman's rule-of-thumb
      - FFT-based convolution for fast computation
      - Linear binning for grid-based estimation
      - Munro-packed real FFT for memory efficiency

      References:
      Silverman BW. (1982) "Algorithm AS 176: Kernel density estimation using the
      Fast Fourier Transform." J. R. Statist. Soc. C 31(1):93-99. DOI: 10.2307/2347084

      @ingroup MathFunctionsStatistics
    */
    struct OPENMS_DLLAPI KernelDensityEstimation
    {
      /**
        @brief Bandwidth selector using the "nrd0" rule-of-thumb for kernel density estimation.

        Computes an appropriate bandwidth for Gaussian kernel density estimation using
        Silverman's rule-of-thumb (also known as "normal reference distribution" or "nrd0").
        This is a quick, robust bandwidth estimator suitable for unimodal distributions.

        The bandwidth is calculated as:
          bw = 0.9 * min(sd, IQR/1.34) * n^(-1/5)

        where sd is standard deviation, IQR is interquartile range, and n is sample size.

        Reference:
        Silverman BW. (1986) "Density Estimation for Statistics and Data Analysis."
        Chapman & Hall/CRC. ISBN 978-0412246203

        @param x Vector of data values for which bandwidth is computed
        @return Bandwidth value (positive scalar)

        @note This is equivalent to R's bw.nrd0() and Python statsmodels' bw_silverman()

        @ingroup MathFunctionsStatistics
      */
      static double bwNrd0(const std::vector<double>& x);

      /**
        @brief Linear binning of data onto an equally-spaced grid.

        Distributes data points onto a regular grid using linear interpolation, which is
        a key preprocessing step for FFT-based kernel density estimation. Each data point
        is allocated to its two nearest grid points proportionally based on distance.

        If weights are provided, weighted counts are computed; otherwise uniform weights
        of 1.0 are used.

        Reference:
        Scott DW. (1992) "Multivariate Density Estimation: Theory, Practice, and Visualization."
        Wiley. ISBN 978-0471547709 (Section on binned kernel estimation)

        @param x Vector of data values to be binned
        @param xmin Minimum value of grid (inclusive)
        @param xmax Maximum value of grid (inclusive)
        @param nbins Number of bins in the grid
        @param weights Optional vector of weights for each data point. If nullptr, uniform weights are used.
        @return Vector of length @p nbins containing (weighted) counts at each grid point

        @note Values outside [xmin, xmax] are ignored (not binned)

        @ingroup MathFunctionsStatistics
      */
      static std::vector<double> linBin(const std::vector<double>& x,
                                        double xmin,
                                        double xmax,
                                        std::size_t nbins,
                                        const std::vector<double>* weights);

      /// Convenience overload that uses uniform weights.
      static std::vector<double> linBin(const std::vector<double>& x,
                                        double xmin,
                                        double xmax,
                                        std::size_t nbins);

      /**
        @brief Forward FFT of real-valued data using Munro-packed format.

        Performs a forward Fast Fourier Transform on real-valued input data, returning
        the result in a compact Munro-packed format. This format takes advantage of the
        Hermitian symmetry property of real FFTs to store all unique frequency information
        in a real-valued vector of the same length as the input.

        The output format is: [Re(Y_0), Re(Y_1), ..., Re(Y_{M/2}), Im(Y_1), ..., Im(Y_{M/2-1})]
        where Y_k are the complex FFT coefficients. The result is UNSCALED (no normalization
        factor applied); the caller must handle any necessary scaling.

        Reference:
        Munro WD. (1976) "Efficient computation of the discrete Fourier transform."
        Also described in: Press WH et al. (2007) "Numerical Recipes" 3rd Ed. Section 12.3

        @param X Real-valued input vector
        @param M Length of FFT (if 0, uses X.size() and rounds up to next power of 2)
        @return Real-valued Munro-packed FFT coefficients of length M (unscaled)

        @see revRt() for the inverse transform (which applies scaling by multiplying by M)

        @ingroup MathFunctionsStatistics
      */
      static std::vector<double> forRt(const std::vector<double>& X, std::size_t M = 0);

      /**
        @brief Inverse FFT of Munro-packed data to real-valued output.

        Performs the inverse operation of forRt(), reconstructing a real-valued signal
        from its Munro-packed frequency domain representation. The output is scaled by
        multiplying by M (applied by the evergreen::real_ifft function).

        The input must be in Munro-packed format:
        [Re(Y_0), Re(Y_1), ..., Re(Y_{M/2}), Im(Y_1), ..., Im(Y_{M/2-1})]

        Reference:
        Munro WD. (1976) "Efficient computation of the discrete Fourier transform."
        Also described in: Press WH et al. (2007) "Numerical Recipes" 3rd Ed. Section 12.3

        @param Xp Munro-packed real-valued frequency coefficients
        @param M Length of output (if 0, uses Xp.size())
        @return Real-valued reconstructed signal of length M (scaled by M)

        @see forRt() for the forward transform (which returns unscaled coefficients)

        @ingroup MathFunctionsStatistics
      */
      static std::vector<double> revRt(const std::vector<double>& Xp, std::size_t M = 0);

      /**
        @brief Compute the FFT of a Gaussian kernel in Munro-packed format.

        Generates the frequency-domain representation of a Gaussian kernel with specified
        bandwidth, suitable for convolution-based kernel density estimation via FFT.
        The result is in Munro-packed format compatible with forRt() and revRt().

        This implementation uses the Silverman transform, which efficiently computes the
        Gaussian kernel in frequency space as: K(f) = exp(-2*pi^2*sigma^2*f^2) where sigma = bw/4.

        Reference:
        Silverman BW. (1982) "Algorithm AS 176: Kernel density estimation using the Fast
        Fourier Transform." Journal of the Royal Statistical Society. Series C (Applied
        Statistics) 31(1):93-99. DOI: 10.2307/2347084

        Also implemented in: Python statsmodels.nonparametric.kdetools.silverman_transform

        @param bw Bandwidth of the Gaussian kernel (standard deviation)
        @param M Length of the FFT grid (number of frequency bins)
        @param RANGE Range of the spatial domain over which the kernel will be applied
        @return Munro-packed FFT coefficients of the Gaussian kernel

        @ingroup MathFunctionsStatistics
      */
      static std::vector<double> silvermanKernelFFT(double bw, std::size_t M, double RANGE);

      /**
        @brief Fast kernel density estimation on a regular grid using FFT convolution.

        Computes Gaussian kernel density estimates at equally-spaced grid points using
        the FFT-based convolution algorithm. This approach scales as O(n + M*log(M))
        compared to O(n*M) for direct evaluation, making it highly efficient for large
        datasets or fine grids.

        The Silverman Algorithm:
        Uses the convolution theorem to compute KDE efficiently in the frequency domain:
        density(x) = IFFT(FFT(binned_data) * FFT(kernel))

        Silverman's analytical kernel representation:
        Instead of FFT-ing a spatial Gaussian kernel, the algorithm computes the kernel
        analytically in the frequency domain using: K(f) = exp(-2*pi^2*sigma^2*f^2) / (1 - f^2*pi^2/3)
        This is both more numerically stable and efficient than FFT-ing the spatial kernel,
        and avoids the need for an additional FFT operation.

        Implementation steps:
        1. Bin data onto a regular grid using linear binning
        2. Compute FFT of binned data using real FFT (via evergreen library)
        3. Compute Silverman's analytical Gaussian kernel in frequency domain
        4. Multiply FFT(data) pointwise by kernel in frequency domain
        5. Inverse FFT to obtain density estimates

        The grid extends from min(x) - cut*bw to max(x) + cut*bw.

        Reference:
        Silverman BW. (1982) "Algorithm AS 176: Kernel density estimation using the Fast
        Fourier Transform." J. R. Statist. Soc. C 31(1):93-99. DOI: 10.2307/2347084

        @param x Data values for which to estimate density
        @param bw Bandwidth of the Gaussian kernel
        @param gridsize Number of equally-spaced grid points (M). Should be power of 2 for FFT efficiency
        @param cut Extension factor: grid extends cut*bw beyond data range on each side
        @return Pair of vectors: (density_values, grid_points), each of length @p gridsize

        @note Returns uniform density if data range is zero (all identical values)

        @ingroup MathFunctionsStatistics
      */
      static std::pair<std::vector<double>, std::vector<double>> gridKdeFFT(const std::vector<double>& x,
                                                                            double bw,
                                                                            std::size_t gridsize = 512,
                                                                            double cut = 3.0);

      /**
        @brief Evaluate kernel density estimates at the data points themselves.

        Computes KDE values at each input data point using the FFT-grid method followed
        by cubic spline interpolation. This is more efficient than direct evaluation for
        large datasets, as the FFT-grid computation scales as O(n + M*log(M)) followed by
        O(n*log(M)) for spline interpolation, compared to O(n^2) for direct methods.

        The function first computes KDE on a regular grid using gridKdeFFT(), then
        interpolates these grid values to the query points using cubic spline interpolation.

        @param x Data values at which to evaluate the density
        @param bw Bandwidth of the Gaussian kernel
        @param gridsize Number of grid points for FFT computation (default 512)
        @param cut Extension factor for grid range (default 3.0)
        @return Vector of KDE values, one for each point in @p x

        @see gridKdeFFT() for the underlying grid-based KDE computation

        @ingroup MathFunctionsStatistics
      */
      static std::vector<double> kdeFFTEval(const std::vector<double>& x,
                                            double bw,
                                            std::size_t gridsize = 512,
                                            double cut = 3.0);
    };

  } // namespace Math
} // namespace OpenMS
