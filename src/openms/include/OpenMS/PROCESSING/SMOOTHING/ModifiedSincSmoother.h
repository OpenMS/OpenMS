#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI ModifiedSincSmoother
  {
  public:
    ModifiedSincSmoother();
    ~ModifiedSincSmoother();
    void smoothData(const std::vector<double>& input,
      std::vector<double>& output,
      const int window_size,
      const double cutoff);

    void smooth(std::vector<double>& data, int degree, int half_kernel_size);

  private:
    std::vector<double> generateKernel(int degree, int m);
    void extendData(const std::vector<double>& data, std::vector<double>& extended, int m, const std::string& boundary);
    std::vector<double> edgeWeights(int degree, int m);
    void fitWeighted(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& w, std::vector<double>& result);
  };
}
