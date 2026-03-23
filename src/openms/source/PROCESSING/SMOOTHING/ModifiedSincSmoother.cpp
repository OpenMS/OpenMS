#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace OpenMS
{

  void ModifiedSincSmoother::registerDefaults_()
  {
    defaults_.setValue("is_ms1", "false",
      "If true, uses the MS1 filter variant (smaller kernel, faster). "
      "If false, uses the MS variant (better stopband suppression).");
    defaults_.setValidStrings("is_ms1", {"true", "false"});
    defaults_.setValue("degree", 6,
      "Polynomial degree of the filter (even, 2-10). Higher degree gives sharper frequency cutoff.");
    defaults_.setMinInt("degree", 2);
    defaults_.setMaxInt("degree", 10);
    defaults_.setValue("m", 7,
      "Kernel half-width. Larger m means more smoothing.");
    defaults_.setMinInt("m", 2);
  }

  ModifiedSincSmoother::ModifiedSincSmoother()
    : ProgressLogger(),
      DefaultParamHandler("ModifiedSincSmoother"),
      isMS1_(false), degree_(6), m_(7)
  {
    registerDefaults_();
    defaultsToParam_(); // copies defaults to param_ and calls updateMembers_()
  }

  ModifiedSincSmoother::ModifiedSincSmoother(bool isMS1, int degree, int m)
    : ProgressLogger(),
      DefaultParamHandler("ModifiedSincSmoother"),
      isMS1_(isMS1), degree_(degree), m_(m)
  {
    if (degree < 2 || degree > 10 || (degree % 2) != 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid degree; must be even and between 2 and 10.");
    }

    const int mMin = isMS1 ? degree / 2 + 1 : degree / 2 + 2;
    if (m < mMin)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid kernel half-width m; too small for the given degree.");
    }

    registerDefaults_();
    defaultsToParam_();
    // override the default param values with the actual constructor arguments,
    // then call updateMembers_() to sync member variables and rebuild the kernel
    param_.setValue("is_ms1", isMS1 ? "true" : "false");
    param_.setValue("degree", degree);
    param_.setValue("m", m);
    updateMembers_();
  }

  void ModifiedSincSmoother::updateMembers_()
  {
    isMS1_ = param_.getValue("is_ms1").toBool();
    degree_ = static_cast<int>(param_.getValue("degree"));
    m_ = static_cast<int>(param_.getValue("m"));

    if (degree_ < 2 || degree_ > 10 || (degree_ % 2) != 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid degree; must be even and between 2 and 10.");
    }
    const int mMin = isMS1_ ? degree_ / 2 + 1 : degree_ / 2 + 2;
    if (m_ < mMin)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid kernel half-width m; too small for the given degree.");
    }

    kernel_ = makeKernel(isMS1_, degree_, m_, getCoefficients(isMS1_, degree_, m_));
    fit_weights_ = makeFitWeights(isMS1_, degree_, m_);
  }

  std::vector<double> ModifiedSincSmoother::smooth(const std::vector<double>& data)
  {
    if (data.empty()) return {};
    const int radius = static_cast<int>(kernel_.size()) - 1; // kernel_.size() == m_ + 1
    std::vector<double> extended = extendData(data, radius);
    std::vector<double> filtered = smoothExceptBoundaries(extended);
    return std::vector<double>(filtered.begin() + radius, filtered.begin() + radius + static_cast<ptrdiff_t>(data.size()));
  }

  std::vector<double> ModifiedSincSmoother::smoothExceptBoundaries(const std::vector<double>& data)
  {
    const size_t N = data.size();
    if (N == 0) return {};
    std::vector<double> out(N, 0.0);
    const int radius = static_cast<int>(kernel_.size()) - 1;

    // compute only where we have full symmetric neighborhood
    for (size_t i = static_cast<size_t>(radius); i + static_cast<size_t>(radius) < N; ++i)
    {
      double sum = kernel_[0] * data[i];
      for (size_t j = 1; j < kernel_.size(); ++j)
      {
        sum += kernel_[j] * (data[i - j] + data[i + j]);
      }
      out[i] = sum;
    }
    return out;
  }

  int ModifiedSincSmoother::bandwidthToM(bool isMS1, int degree, double bandwidth)
  {
    if (bandwidth <= 0.0 || bandwidth >= 0.5)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid bandwidth value. Expected in (0, 0.5).");
    }
    const double radius = isMS1
      ? (0.27037 + 0.24920 * degree) / bandwidth - 1.0
      : (0.74548 + 0.24943 * degree) / bandwidth - 1.0;
    return static_cast<int>(std::round(radius));
  }

  int ModifiedSincSmoother::noiseGainToM(bool isMS1, int degree, double noiseGain)
  {
    if (noiseGain <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Invalid noiseGain value; must be > 0.");
    }
    const double invNG2 = 1.0 / (noiseGain * noiseGain);
    const double e = -2.5 - 0.8 * degree;
    const double m = isMS1
      ? (-1.0 + invNG2 * (0.543 + 0.4974 * degree) + 0.47 * std::pow(invNG2, e))
      : (-1.0 + invNG2 * (1.494 + 0.4965 * degree) + 0.52 * std::pow(invNG2, e));
    return static_cast<int>(std::round(m));
  }

  double ModifiedSincSmoother::savitzkyGolayBandwidth(int degree, int m)
  {
    return 1.0 / (6.352 * (m + 0.5) / (degree + 1.379) - (0.513 + 0.316 * degree) / (m + 0.5));
  }

  std::vector<double> ModifiedSincSmoother::makeKernel(bool isMS1, int degree, int m, const std::vector<double>& coeffs)
  {
    std::vector<double> kernel(static_cast<size_t>(m) + 1u, 0.0);
    double sum = 0.0;

    for (int i = 0; i <= m; ++i)
    {
      const double x = static_cast<double>(i) / (m + 1);
      const double sinc_arg = Constants::PI * 0.5 * (isMS1 ? degree + 2 : degree + 4) * x;
      double k = (i == 0) ? 1.0 : std::sin(sinc_arg) / sinc_arg;

      // optional correction terms
      for (size_t j = 0; j < coeffs.size(); ++j)
      {
        if (isMS1)
        {
          k += coeffs[j] * x * std::sin((j + 1) * Constants::PI * x);
        }
        else
        {
          const int nu = ((degree / 2) % 2 == 0) ? 2 : 1;
          k += coeffs[j] * x * std::sin((2 * static_cast<int>(j) + nu) * Constants::PI * x);
        }
      }

      const double decay = isMS1 ? 2.0 : 4.0;
      k *= std::exp(-x * x * decay)
         + std::exp(-(x - 2) * (x - 2) * decay)
         + std::exp(-(x + 2) * (x + 2) * decay)
         - 2 * std::exp(-decay)
         - std::exp(-9 * decay);

      kernel[static_cast<size_t>(i)] = k;
      sum += (i == 0) ? k : 2.0 * k; // account for symmetric duplication
    }

    // normalize to unity gain at DC
    if (sum != 0.0)
    {
      for (double& v : kernel) v /= sum;
    }
    return kernel;
  }

  std::vector<double> ModifiedSincSmoother::getCoefficients(bool isMS1, int degree, int m)
  {
    struct ABC { double a, b, c; };
    std::vector<ABC> table;

    if (isMS1)
    {
      switch (degree)
      {
        case 2: break;
        case 4:
          table = {{0.021944195, 0.050284006, 0.765625}};
          break;
        case 6:
          table = {{0.0018977303, 0.008476806, 1.2625},
                   {0.023064667, 0.13047926, 1.2265625}};
          break;
        case 8:
          table = {{0.0065903002, 0.057929456, 1.915625},
                   {0.0023234477, 0.010298849, 2.2726562},
                   {0.021046653, 0.16646601, 1.98125}};
          break;
        case 10:
          table = {{9.749618e-4, 0.0020742896, 3.74375},
                   {0.008975366, 0.09902466, 2.7078125},
                   {0.0024195414, 0.010064855, 3.296875},
                   {0.019185117, 0.18953617, 2.784961}};
          break;
        default: break;
      }
    }
    else
    {
      switch (degree)
      {
        case 2: break;
        case 4: break;
        case 6:
          table = {{0.001717576, 0.02437382, 1.64375}};
          break;
        case 8:
          table = {{0.0043993373, 0.088211164, 2.359375},
                   {0.006146815, 0.024715371, 3.6359375}};
          break;
        case 10:
          table = {{0.0011840032, 0.04219344, 2.746875},
                   {0.0036718843, 0.12780383, 2.7703125}};
          break;
        default: break;
      }
    }

    std::vector<double> coeffs;
    coeffs.reserve(table.size());
    for (const auto& abc : table)
    {
      double cm = abc.c - m;
      coeffs.push_back(abc.a + abc.b / (cm * cm * cm));
    }
    return coeffs;
  }

  std::vector<double> ModifiedSincSmoother::makeFitWeights(bool isMS1, int degree, int m)
  {
    // derive boundary fit window and weights from first zero and empirical beta
    const double first_zero = isMS1
      ? static_cast<double>(m + 1) / (1.0 + 0.5 * degree)
      : static_cast<double>(m + 1) / (1.5 + 0.5 * degree);
    const double beta = isMS1
      ? 0.65 + 0.35 * std::exp(-0.55 * (degree - 4))
      : 0.70 + 0.14 * std::exp(-0.60 * (degree - 4));
    const int fit_len = static_cast<int>(std::ceil(first_zero * beta));

    std::vector<double> weights(static_cast<size_t>(std::max(fit_len, 0)), 0.0);
    for (int p = 0; p < fit_len; ++p)
    {
      weights[static_cast<size_t>(p)] = sqr(std::cos(0.5 * Constants::PI / (first_zero * beta) * p));
    }
    return weights;
  }

  std::vector<double> ModifiedSincSmoother::extendData(const std::vector<double>& data, int m)
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

  void ModifiedSincSmoother::filter(MSSpectrum& spectrum)
  {
    // copy data AND meta data
    MSSpectrum output = spectrum;

    // collect intensities
    std::vector<double> intensities;
    intensities.reserve(spectrum.size());
    for (auto it = spectrum.begin(); it != spectrum.end(); ++it)
    {
      intensities.push_back(it->getIntensity());
    }

    // smooth intensities
    std::vector<double> smoothed = smooth(intensities);

    // write back, preserve position
    auto out_it = output.begin();
    auto in_it = spectrum.begin();
    for (size_t i = 0; i < smoothed.size() && out_it != output.end() && in_it != spectrum.end(); ++i, ++out_it, ++in_it)
    {
      out_it->setPosition(in_it->getPosition());
      out_it->setIntensity(std::max(0.0, smoothed[i]));
    }

    // swap back
    std::swap(spectrum, output);
  }

  void ModifiedSincSmoother::filter(MSChromatogram& chromatogram)
  {
    MSChromatogram output = chromatogram;

    std::vector<double> intensities;
    intensities.reserve(chromatogram.size());
    for (auto it = chromatogram.begin(); it != chromatogram.end(); ++it)
    {
      intensities.push_back(it->getIntensity());
    }

    std::vector<double> smoothed = smooth(intensities);

    auto out_it = output.begin();
    auto in_it = chromatogram.begin();
    for (size_t i = 0; i < smoothed.size() && out_it != output.end() && in_it != chromatogram.end(); ++i, ++out_it, ++in_it)
    {
      out_it->setPosition(in_it->getPosition());
      out_it->setIntensity(std::max(0.0, smoothed[i]));
    }

    std::swap(chromatogram, output);
  }

  void ModifiedSincSmoother::filter(Mobilogram& mobilogram)
  {
    Mobilogram output = mobilogram;

    std::vector<double> intensities;
    intensities.reserve(mobilogram.size());
    for (auto it = mobilogram.begin(); it != mobilogram.end(); ++it)
    {
      intensities.push_back(it->getIntensity());
    }

    std::vector<double> smoothed = smooth(intensities);

    auto out_it = output.begin();
    auto in_it = mobilogram.begin();
    for (size_t i = 0; i < smoothed.size() && out_it != output.end() && in_it != mobilogram.end(); ++i, ++out_it, ++in_it)
    {
      out_it->setPosition(in_it->getPosition());
      out_it->setIntensity(std::max(0.0, smoothed[i]));
    }

    std::swap(mobilogram, output);
  }

  void ModifiedSincSmoother::filterExperiment(PeakMap& map)
  {
    Size progress = 0;
    startProgress(0, map.size() + map.getChromatograms().size(), "smoothing data");
    for (Size i = 0; i < map.size(); ++i)
    {
      filter(map[i]);
      setProgress(++progress);
    }
    for (Size i = 0; i < map.getChromatograms().size(); ++i)
    {
      filter(map.getChromatogram(i));
      setProgress(++progress);
    }
    endProgress();
  }

  void ModifiedSincSmoother::LinearRegression::clear()
  {
    sum_weights = sum_x = sum_y = sum_xy = sum_x2 = 0;
    calculated = false;
  }

  void ModifiedSincSmoother::LinearRegression::addPointW(double x, double y, double weight)
  {
    sum_weights += weight;
    sum_x += weight * x;
    sum_y += weight * y;
    sum_xy += weight * x * y;
    sum_x2 += weight * x * x;
    calculated = false;
  }

  void ModifiedSincSmoother::LinearRegression::calculate()
  {
    double denom = sum_x2 - (sum_x * sum_x) / sum_weights;
    slope = (sum_xy - sum_x * sum_y / sum_weights) / denom;
    if (std::isnan(slope)) slope = 0;
    offset = (sum_y - slope * sum_x) / sum_weights;
    calculated = true;
  }

  double ModifiedSincSmoother::LinearRegression::getOffset()
  {
    if (!calculated) calculate();
    return offset;
  }

  double ModifiedSincSmoother::LinearRegression::getSlope()
  {
    if (!calculated) calculate();
    return slope;
  }

} // namespace OpenMS