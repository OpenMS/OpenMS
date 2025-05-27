#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <algorithm>

namespace OpenMS
{
  ModifiedSincSmoother::ModifiedSincSmoother() = default;
  ModifiedSincSmoother::~ModifiedSincSmoother() = default;

  void ModifiedSincSmoother::smoothData(
    const std::vector<double>& input,
    std::vector<double>& output,
    const int window_size,
    const double cutoff)
  {
    if (window_size % 2 == 0 || window_size < 3)
    {
      throw std::invalid_argument("Window size must be an odd number >= 3.");
    }

    const int m = window_size / 2;
    const int N = static_cast<int>(input.size());
    output.resize(N);

    std::vector<double> sinc_kernel(window_size);
    const double pi = M_PI;
    for (int i = -m; i <= m; ++i)
    {
      if (i == 0)
      {
        sinc_kernel[i + m] = 2.0 * cutoff;
      }
      else
      {
        double x = pi * i;
        sinc_kernel[i + m] = std::sin(2.0 * cutoff * x) / x;
      }
    }

    // Apply Hamming window
    for (int i = 0; i < window_size; ++i)
    {
      sinc_kernel[i] *= 0.54 - 0.46 * std::cos(2.0 * pi * i / (window_size - 1));
    }

    // Normalize kernel
    const double sum = std::accumulate(sinc_kernel.begin(), sinc_kernel.end(), 0.0);
    for (double& val : sinc_kernel)
    {
      val /= sum;
    }

    // Apply convolution
    for (int i = 0; i < N; ++i)
    {
      double smoothed = 0.0;
      for (int j = -m; j <= m; ++j)
      {
        int idx = std::max(0, std::min(i + j, N - 1)); // handle borders
        smoothed += input[idx] * sinc_kernel[j + m];
      }
      output[i] = smoothed;
    }
  }
}
