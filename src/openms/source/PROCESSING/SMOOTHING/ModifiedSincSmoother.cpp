#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

namespace OpenMS
{
  class ModifiedSincSmoother
  {
  public:
    ModifiedSincSmoother(bool isMS1, int degree, int m)
      : isMS1_(isMS1), degree_(degree), m_(m)
    {
      if (degree < 2 || degree > 10 || (degree % 2) != 0)
        throw std::invalid_argument("Invalid degree; must be even and between 2 and 10.");

      int mMin = isMS1 ? degree / 2 + 1 : degree / 2 + 2;
      if (m < mMin)
        throw std::invalid_argument("Invalid kernel half-width m; too small for the given degree.");

      kernel_ = makeKernel(isMS1_, degree_, m_, getCoefficients());
      fit_weights_ = makeFitWeights();
    }

    std::vector<double> smooth(const std::vector<double>& data)
    {
      int radius = kernel_.size() - 1;
      std::vector<double> extended = extendData(data, radius);
      std::vector<double> filtered = smoothExceptBoundaries(extended);
      return std::vector<double>(filtered.begin() + radius, filtered.begin() + radius + data.size());
    }

    std::vector<double> smoothExceptBoundaries(const std::vector<double>& data)
    {
      std::vector<double> out(data.size(), 0.0);
      int radius = kernel_.size() - 1;
      for (size_t i = radius; i < data.size() - radius; ++i)
      {
        double sum = kernel_[0] * data[i];
        for (size_t j = 1; j < kernel_.size(); ++j)
          sum += kernel_[j] * (data[i - j] + data[i + j]);
        out[i] = sum;
      }
      return out;
    }

    static int bandwidthToM(bool isMS1, int degree, double bandwidth)
    {
      if (bandwidth <= 0 || bandwidth >= 0.5)
        throw std::invalid_argument("Invalid bandwidth value.");
      double radius = isMS1 ?
        (0.27037 + 0.24920 * degree) / bandwidth - 1.0 :
        (0.74548 + 0.24943 * degree) / bandwidth - 1.0;
      return static_cast<int>(std::round(radius));
    }

    static int noiseGainToM(bool isMS1, int degree, double noiseGain)
    {
      double invNG2 = 1. / (noiseGain * noiseGain);
      double exp = -2.5 - 0.8 * degree;
      double m = isMS1 ?
        -1 + invNG2 * (0.543 + 0.4974 * degree) + 0.47 * std::pow(invNG2, exp) :
        -1 + invNG2 * (1.494 + 0.4965 * degree) + 0.52 * std::pow(invNG2, exp);
      return static_cast<int>(std::round(m));
    }

    static double savitzkyGolayBandwidth(int degree, int m)
    {
      return 1. / (6.352 * (m + 0.5) / (degree + 1.379) - (0.513 + 0.316 * degree) / (m + 0.5));
    }

  private:
    bool isMS1_;
    int degree_;
    int m_;
    std::vector<double> kernel_;
    std::vector<double> fit_weights_;

    std::vector<double> makeKernel(bool isMS1, int degree, int m, const std::vector<double>& coeffs)
    {
      std::vector<double> kernel(m + 1, 0.0);
      double sum = 0.0;
      for (int i = 0; i <= m; ++i)
      {
        double x = static_cast<double>(i) / (m + 1);
        double sinc_arg = M_PI * 0.5 * (isMS1 ? degree + 2 : degree + 4) * x;
        double k = (i == 0) ? 1.0 : std::sin(sinc_arg) / sinc_arg;

        for (size_t j = 0; j < coeffs.size(); ++j)
        {
          if (isMS1)
            k += coeffs[j] * x * std::sin((j + 1) * M_PI * x);
          else
          {
            int nu = ((degree / 2) % 2 == 0) ? 2 : 1;
            k += coeffs[j] * x * std::sin((2 * j + nu) * M_PI * x);
          }
        }

        double decay = isMS1 ? 2.0 : 4.0;
        k *= std::exp(-x * x * decay) + std::exp(-(x - 2) * (x - 2) * decay) +
             std::exp(-(x + 2) * (x + 2) * decay) - 2 * std::exp(-decay) - std::exp(-9 * decay);

        kernel[i] = k;
        sum += (i == 0) ? k : 2 * k;
      }

      for (auto& val : kernel)
        val /= sum;

      return kernel;
    }

    std::vector<double> getCoefficients()
    {
      return {}; // Placeholder: Static correction data would go here
    }

    std::vector<double> makeFitWeights()
    {
      double first_zero = isMS1_ ? (m_ + 1) / (1 + 0.5 * degree_) : (m_ + 1) / (1.5 + 0.5 * degree_);
      double beta = isMS1_ ? 0.65 + 0.35 * std::exp(-0.55 * (degree_ - 4)) :
                             0.70 + 0.14 * std::exp(-0.60 * (degree_ - 4));
      int fit_len = static_cast<int>(std::ceil(first_zero * beta));
      std::vector<double> weights(fit_len);
      for (int p = 0; p < fit_len; ++p)
        weights[p] = sqr(std::cos(0.5 * M_PI / (first_zero * beta) * p));
      return weights;
    }

    std::vector<double> extendData(const std::vector<double>& data, int m)
    {
      std::vector<double> extended(data.size() + 2 * m);
      std::copy(data.begin(), data.end(), extended.begin() + m);

      LinearRegression linreg;
      int fitLength = std::min(static_cast<int>(fit_weights_.size()), static_cast<int>(data.size()));

      for (int p = 0; p < fitLength; ++p)
        linreg.addPointW(p, data[p], fit_weights_[p]);

      double offset = linreg.getOffset();
      double slope = linreg.getSlope();
      for (int p = -1; p >= -m; --p)
        extended[m + p] = offset + slope * p;

      linreg.clear();
      for (int p = 0; p < fitLength; ++p)
        linreg.addPointW(p, data[data.size() - 1 - p], fit_weights_[p]);

      offset = linreg.getOffset();
      slope = linreg.getSlope();
      for (int p = -1; p >= -m; --p)
        extended[data.size() + m - 1 - p] = offset + slope * p;

      return extended;
    }

    struct LinearRegression
    {
      double sum_weights = 0, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
      double offset = NAN, slope = NAN;
      bool calculated = false;

      void clear()
      {
        sum_weights = sum_x = sum_y = sum_xy = sum_x2 = sum_y2 = 0;
        calculated = false;
      }

      void addPointW(double x, double y, double w)
      {
        sum_weights += w;
        sum_x += w * x;
        sum_y += w * y;
        sum_xy += w * x * y;
        sum_x2 += w * x * x;
        sum_y2 += w * y * y;
        calculated = false;
      }

      void calculate()
      {
        double denom = sum_x2 - (sum_x * sum_x) / sum_weights;
        slope = (sum_xy - sum_x * sum_y / sum_weights) / denom;
        if (std::isnan(slope)) slope = 0;
        offset = (sum_y - slope * sum_x) / sum_weights;
        calculated = true;
      }

      double getOffset()
      {
        if (!calculated) calculate();
        return offset;
      }

      double getSlope()
      {
        if (!calculated) calculate();
        return slope;
      }
    };

    static double sqr(double x) { return x * x; }
  };
} // namespace OpenMS