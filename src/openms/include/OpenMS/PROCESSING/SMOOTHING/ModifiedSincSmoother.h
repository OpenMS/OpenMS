// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>
#include <cmath>

namespace OpenMS
{
  /**
    @brief Computes a modified sinc smoothing filter for mass spectrometry data.

    This class implements a modified sinc-based smoothing filter optimized for mass spectrometry data.
    The filter supports two main variants:
    - MS1 mode: optimized for full scan mass spectra
    - MS2 mode: optimized for fragmentation spectra
    
    The algorithm constructs a symmetric convolution kernel based on a modified sinc function with
    exponential decay terms. The kernel shape is determined by:
    - @p degree: polynomial degree (must be even, between 2-10)
    - @p m: kernel half-width parameter (determines filter width)
    - @p isMS1: selects MS1 vs MS2 optimization variant

    @b Boundary @b Handling:
    For finite data, the filter extends boundaries using weighted linear regression
    on the nearest available data points to minimize edge artifacts.

    @b Container @b Support:
    The filter can process individual data vectors or OpenMS container types:
    - MSSpectrum: preserves m/z positions, smooths intensities
    - MSChromatogram: preserves RT positions, smooths intensities  
    - Mobilogram: preserves IM positions, smooths intensities
    - PeakMap (MSExperiment): batch processing with progress logging

    @b Usage @b Requirements:
    - Input data should represent uniform profile data for optimal results
    - Data must have sufficient length relative to kernel width (≥ 2*m+1 points recommended)
    - For container filters, positions are preserved exactly while only intensities are modified

    @b Parameter @b Helpers:
    Static utility methods provide conversions between different parameterizations:
    - bandwidthToM(): convert frequency domain bandwidth to spatial domain m
    - noiseGainToM(): convert desired noise gain to optimal m value
    - savitzkyGolayBandwidth(): compute equivalent Savitzky-Golay bandwidth

    @note This filter works optimally with uniform profile data!
    @note Data should be sorted by position (m/z, RT, or IM) before filtering.
    @note Kernel parameters must satisfy degree ∈ [2,4,6,8,10] and m ≥ degree/2 + (1 or 2).

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI ModifiedSincSmoother :
    public ProgressLogger
  {
  public:
    /**
      @brief Constructor initializing the modified sinc smoother.

      @param isMS1 true for MS1 optimization, false for MS2 optimization
      @param degree polynomial degree for sinc modification (must be even, 2-10)
      @param m kernel half-width parameter (spatial domain)

      @throws std::invalid_argument if degree is not even or outside [2,10]
      @throws std::invalid_argument if m is too small for the given degree
        (minimum: m ≥ degree/2 + 1 for MS1, m ≥ degree/2 + 2 for MS2)
    */
    ModifiedSincSmoother(bool isMS1, int degree, int m);

    /**
      @brief Smooth a vector of data values with boundary extrapolation.

      Applies the modified sinc filter to the input data. Boundaries are extended
      using weighted linear regression before filtering, then the result is 
      trimmed to match the original input size.

      @param data input intensity values to be smoothed
      @return smoothed data vector of same size as input

      @note Returns empty vector if input is empty
      @note All values in returned vector correspond to positions in original data
    */
    std::vector<double> smooth(const std::vector<double>& data);

    /**
      @brief Smooth data without boundary extension (interior points only).

      Applies the modified sinc filter only where a full symmetric neighborhood
      exists. Boundary points are left as zeros. This method is typically used
      with pre-extended data.

      @param data input data (potentially pre-extended)
      @return smoothed data vector of same size, with boundary zeros

      @note Use this method when you have already extended the data externally
      @note Boundary points (within kernel radius) will be zero in output
    */
    std::vector<double> smoothExceptBoundaries(const std::vector<double>& data);

    /**
      @brief Convert frequency domain bandwidth to spatial domain parameter m.

      @param isMS1 true for MS1 mode, false for MS2 mode  
      @param degree polynomial degree used in the filter
      @param bandwidth normalized frequency bandwidth in (0, 0.5)
      @return corresponding kernel half-width parameter m

      @throws std::invalid_argument if bandwidth ∉ (0, 0.5)

      @note Bandwidth interpretation: fraction of Nyquist frequency
      @note Smaller bandwidth → larger m → more smoothing
    */
    static int bandwidthToM(bool isMS1, int degree, double bandwidth);

    /**
      @brief Convert desired noise gain to optimal parameter m.

      Computes the kernel half-width m that achieves the specified noise gain
      factor for the given filter configuration.

      @param isMS1 true for MS1 mode, false for MS2 mode
      @param degree polynomial degree used in the filter  
      @param noiseGain desired noise amplification factor (typically < 1)
      @return corresponding kernel half-width parameter m

      @note noiseGain = 1 means noise unchanged, < 1 means noise reduction
      @note Formula based on empirical characterization of filter response
    */
    static int noiseGainToM(bool isMS1, int degree, double noiseGain);

    /**
      @brief Compute equivalent Savitzky-Golay bandwidth for comparison.

      Given modified sinc parameters, computes the frequency bandwidth that
      a Savitzky-Golay filter of the same degree would need to achieve
      similar smoothing characteristics.

      @param degree polynomial degree
      @param m kernel half-width parameter
      @return equivalent Savitzky-Golay bandwidth

      @note Useful for comparing filter characteristics across algorithms
      @note Based on empirical matching of filter frequency responses
    */
    static double savitzkyGolayBandwidth(int degree, int m);

    /**
      @brief Apply modified sinc smoothing to an MSSpectrum.

      Smooths the intensity values of a mass spectrum while preserving
      all m/z positions and metadata. Uses copy-filter-swap pattern
      for exception safety.

      @param spectrum input/output spectrum to be smoothed in-place

      @note Requires uniform profile data for optimal results
      @note m/z positions are preserved exactly; only intensities modified
      @note Negative intensities are clipped to zero after smoothing
    */
    void filter(MSSpectrum& spectrum);

    /**
      @brief Apply modified sinc smoothing to an MSChromatogram.

      Smooths the intensity values of a chromatogram while preserving
      all retention time positions and metadata.

      @param chromatogram input/output chromatogram to be smoothed in-place

      @note RT positions are preserved exactly; only intensities modified  
      @note Negative intensities are clipped to zero after smoothing
    */
    void filter(MSChromatogram& chromatogram);

    /**
      @brief Apply modified sinc smoothing to a Mobilogram.

      Smooths the intensity values of an ion mobility spectrum while preserving
      all mobility positions and metadata.

      @param mobilogram input/output mobilogram to be smoothed in-place

      @note IM positions are preserved exactly; only intensities modified
      @note Negative intensities are clipped to zero after smoothing  
    */
    void filter(Mobilogram& mobilogram);

    /**
      @brief Apply modified sinc smoothing to all spectra and chromatograms in an experiment.

      Batch processes all MSSpectra and MSChromatograms in the peak map,
      applying the filter to each container. Progress is reported via
      the ProgressLogger interface.

      @param map input/output experiment to be smoothed in-place

      @note Processes both spectra and chromatograms
      @note Progress is reported as "smoothing data" with step counts
      @note All metadata and positions are preserved across containers
    */
    void filterExperiment(PeakMap& map);

  private:
    /// MS1 vs MS2 mode selection
    bool isMS1_;
    /// Polynomial degree for sinc modification (even, 2-10)
    int degree_;
    /// Kernel half-width parameter (spatial domain)
    int m_;
    /// Symmetric convolution kernel coefficients [k_0, k_1, ..., k_m]
    std::vector<double> kernel_;
    /// Weights for boundary linear regression fitting
    std::vector<double> fit_weights_;

    /**
      @brief Construct the modified sinc convolution kernel.

      @param isMS1 MS1 vs MS2 mode selection
      @param degree polynomial degree for modifications
      @param m kernel half-width
      @param coeffs optional correction coefficients (currently unused)
      @return normalized symmetric kernel [k_0, k_1, ..., k_m]

      @note Kernel is normalized to unity gain at DC (frequency 0)
      @note Incorporates exponential decay terms for improved characteristics
    */
    std::vector<double> makeKernel(bool isMS1, int degree, int m, const std::vector<double>& coeffs);

    /**
      @brief Get correction coefficients for kernel modification.

      @param isMS1 MS1 vs MS2 mode selection  
      @param degree polynomial degree
      @param m kernel half-width
      @return correction coefficients (currently empty placeholder)

      @note Currently returns empty vector; reserved for future enhancements
    */
    std::vector<double> getCoefficients(bool isMS1, int degree, int m);

    /**
      @brief Construct weights for boundary linear regression.

      Derives fitting window size and cosine-squared weights based on
      kernel first zero location and empirical beta parameter.

      @param isMS1 MS1 vs MS2 mode selection
      @param degree polynomial degree  
      @param m kernel half-width
      @return weight vector for boundary fitting

      @note Weight derivation: first_zero → beta → fit_length → cos²(·) weights
      @note Optimized separately for MS1 and MS2 boundary characteristics
    */
    std::vector<double> makeFitWeights(bool isMS1, int degree, int m);

    /**
      @brief Extend data boundaries using weighted linear regression.

      Fits linear models to boundary regions using cosine-squared weights,
      then extrapolates to extend the data symmetrically by m points on each side.

      @param data original data vector
      @param m extension half-width (points to add on each side)
      @param degree polynomial degree (affects fit length)
      @return extended data: [left_ext..., original_data..., right_ext...]

      @note Left extension: fit forward from start, extrapolate backward
      @note Right extension: fit backward from end, extrapolate forward  
      @note Uses fit_weights_ for regression weighting
    */
    std::vector<double> extendData(const std::vector<double>& data, int m, int degree);

    /**
      @brief Helper class for weighted linear regression.

      Accumulates weighted (x,y) points and computes least-squares
      line fit on demand. Used for boundary extrapolation.
    */
    struct LinearRegression
    {
      double sum_weights = 0;  ///< Sum of weights
      double sum_x = 0;        ///< Sum of weighted x values
      double sum_y = 0;        ///< Sum of weighted y values  
      double sum_xy = 0;       ///< Sum of weighted x*y products
      double sum_x2 = 0;       ///< Sum of weighted x² values
      double sum_y2 = 0;       ///< Sum of weighted y² values
      double offset = NAN;     ///< Regression intercept (cached)
      double slope = NAN;      ///< Regression slope (cached)
      bool calculated = false; ///< Whether regression has been computed

      /**
        @brief Reset accumulator for new regression.
      */
      void clear();

      /**
        @brief Add a weighted data point to the regression.

        @param x independent variable value
        @param y dependent variable value  
        @param weight point weight (≥ 0)
      */
      void addPointW(double x, double y, double weight);

      /**
        @brief Compute least-squares line parameters.

        Solves normal equations for weighted linear regression.
        Sets calculated = true and caches offset/slope.

        @note Handles singular cases by setting slope = 0
        @note Automatically called by getOffset()/getSlope() if needed
      */
      void calculate();

      /**
        @brief Get regression intercept (y-value at x=0).

        @return regression offset parameter
        @note Triggers calculate() if not already computed
      */
      double getOffset();

      /**
        @brief Get regression slope (dy/dx).

        @return regression slope parameter  
        @note Triggers calculate() if not already computed
      */
      double getSlope();
    };

    /// Helper: square function
    static double sqr(double x) { return x * x; }
  };
} // namespace OpenMS
