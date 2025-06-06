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
    const double cutoff_frequency)
  {
    if (window_size % 2 == 0 || window_size < 3)
    {
      throw std::invalid_argument("Window size must be an odd number >= 3.");
    }

    if (cutoff_frequency <= 0.0 || cutoff_frequency > 1.0)
    {
      throw std::invalid_argument("Cutoff frequency must be between 0 and 1.");
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
        sinc_kernel[i + m] = 2.0 * cutoff_frequency;
      }
      else
      {
        double x = pi * i;
        sinc_kernel[i + m] = std::sin(2.0 * cutoff_frequency * x) / x;
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

  void ModifiedSincSmoother::smooth(std::vector<double>& data, int degree, int half_kernel_size)
  {
    const int window_size = 2 * half_kernel_size + 1;
    const double cutoff = static_cast<double>(degree) / (degree + half_kernel_size);
    std::vector<double> result;
    smoothData(data, result, window_size, cutoff);
    data = std::move(result);
  }

  std::vector<double> ModifiedSincSmoother::generateKernel(int degree, int m)
  {
    std::vector<double> kernel(2 * m + 1);
    for (int i = -m; i <= m; ++i)
    {
      kernel[i + m] = std::exp(-std::pow(static_cast<double>(i), 2) / (2 * degree));
    }
    double sum = std::accumulate(kernel.begin(), kernel.end(), 0.0);
    for (double& val : kernel) val /= sum;
    return kernel;
  }

  void ModifiedSincSmoother::extendData(const std::vector<double>& data, std::vector<double>& extended, int m, const std::string& boundary)
  {
    const int N = static_cast<int>(data.size());
    extended.resize(N + 2 * m);

    if (boundary == "mirror")
    {
      for (int i = 0; i < m; ++i)
      {
        extended[i] = data[m - i];
        extended[N + m + i] = data[N - 2 - i];
      }
    }
    else
    {
      for (int i = 0; i < m; ++i)
      {
        extended[i] = data.front();
        extended[N + m + i] = data.back();
      }
    }
    std::copy(data.begin(), data.end(), extended.begin() + m);
  }

  std::vector<double> ModifiedSincSmoother::edgeWeights(int degree, int m)
  {
    std::vector<double> weights(2 * m + 1, 1.0);
    for (int i = 0; i < m; ++i)
    {
      weights[i] = weights[2 * m - i] = static_cast<double>(i + 1) / (m + 1);
    }
    return weights;
  }

  void ModifiedSincSmoother::fitWeighted(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& w, std::vector<double>& result)
  {
    const size_t N = x.size();
    result.resize(N);
    for (size_t i = 0; i < N; ++i)
    {
      result[i] = x[i] * w[i];
    }
  }
} // namespace OpenMS
