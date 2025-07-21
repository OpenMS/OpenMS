#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>

namespace OpenMS
{
  class ModifiedSincSmoother
  {
  public:
    ModifiedSincSmoother(bool isMS1, int degree, int m);

    std::vector<double> smooth(const std::vector<double>& data);
    std::vector<double> smoothExceptBoundaries(const std::vector<double>& data);
    static int bandwidthToM(bool isMS1, int degree, double bandwidth);
    static int noiseGainToM(bool isMS1, int degree, double noiseGain);
    static double savitzkyGolayBandwidth(int degree, int m);

  private:
    bool isMS1_;
    int degree_;
    int m_;
    std::vector<double> kernel_;
    std::vector<double> fit_weights_;

    std::vector<double> makeKernel(bool isMS1, int degree, int m, const std::vector<double>& coeffs);
    std::vector<double> getCoefficients(bool isMS1, int degree, int m);
    std::vector<double> makeFitWeights(bool isMS1, int degree, int m);
    std::vector<double> extendData(const std::vector<double>& data, int m, int degree);

    struct LinearRegression
    {
      double sum_weights = 0;
      double sum_x = 0;
      double sum_y = 0;
      double sum_xy = 0;
      double sum_x2 = 0;
      double sum_y2 = 0;
      double offset = NAN;
      double slope = NAN;
      bool calculated = false;

      void clear()
      {
        sum_weights = sum_x = sum_y = sum_xy = sum_x2 = sum_y2 = 0;
        calculated = false;
      }

      void addPointW(double x, double y, double weight)
      {
        sum_weights += weight;
        sum_x += weight * x;
        sum_y += weight * y;
        sum_xy += weight * x * y;
        sum_x2 += weight * x * x;
        sum_y2 += weight * y * y;
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
