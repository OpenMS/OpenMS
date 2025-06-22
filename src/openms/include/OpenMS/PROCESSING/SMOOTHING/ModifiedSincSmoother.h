#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <vector>
#include <string>

namespace OpenMS
{
  /**
   * @brief Implements modified sinc-based smoothing algorithms for data processing.
   *
   * The ModifiedSincSmoother class provides smoothing functionality using 
   * modified sinc filters with windowing. It supports both direct convolution-based
   * smoothing and polynomial-based smoothing approaches.
   */
  class OPENMS_DLLAPI ModifiedSincSmoother
  {
  public:
    ModifiedSincSmoother();
    ~ModifiedSincSmoother();

    /**
     * @brief Smooths input data using a modified sinc filter with Hamming windowing.
     *
     * @param input Input data vector to be smoothed
     * @param output Output vector to store smoothed results
     * @param window_size Size of the smoothing window (must be odd, >= 3)
     * @param cutoff_frequency Normalized cutoff frequency (0 < cutoff <= 1)
     */
    void smoothData(const std::vector<double>& input,
                    std::vector<double>& output,
                    const int window_size,
                    const double cutoff_frequency);

    /**
     * @brief Smooths data in-place using the modified sinc algorithm.
     *
     * @param data Vector of data points to be smoothed.
     * @param degree Polynomial degree for the smoothing filter.
     * @param half_kernel_size Half-size of the smoothing kernel.
     */
    void smooth(std::vector<double>& data, int degree, int half_kernel_size);

  private:
    /** @brief Extends data with boundary conditions. */
    void extendData(const std::vector<double>& data, std::vector<double>& extended, int m, const std::string& boundary);

    /** @brief Computes edge weights for boundary handling. */
    std::vector<double> edgeWeights(int degree, int m);

    /** @brief Performs weighted fitting on data points. */
    void fitWeighted(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& w, std::vector<double>& result);
  };
}
