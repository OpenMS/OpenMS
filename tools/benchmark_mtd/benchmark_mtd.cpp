// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS contributors $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <stdexcept>

using namespace OpenMS;

namespace
{
bool equalTraces(const std::vector<MassTrace>& a, const std::vector<MassTrace>& b)
{
  if (a.size() != b.size()) { return false; }
  for (Size i = 0; i < a.size(); ++i)
  {
    if (a[i].getSize() != b[i].getSize() || a[i].getLabel() != b[i].getLabel() || a[i].getCentroidMZ() != b[i].getCentroidMZ()
        || a[i].getCentroidRT() != b[i].getCentroidRT() || a[i].getCentroidIM() != b[i].getCentroidIM()
        || a[i].getCentroidSD() != b[i].getCentroidSD() || a[i].fwhm_mz_avg != b[i].fwhm_mz_avg || a[i].fwhm_im_avg != b[i].fwhm_im_avg
        || a[i].computePeakArea() != b[i].computePeakArea() || ! std::equal(a[i].begin(), a[i].end(), b[i].begin()))
    {
      return false;
    }
  }
  return true;
}

double median(std::vector<double> values)
{
  std::sort(values.begin(), values.end());
  const Size middle = values.size() / 2;
  return values.size() % 2 == 0 ? (values[middle - 1] + values[middle]) / 2.0 : values[middle];
}
} // namespace

int main(int argc, char** argv)
{
  if (argc < 2 || argc > 5)
  {
    OPENMS_LOG_ERROR << "Usage: benchmark_mtd input.mzML [threads=2] [repeats=9] [apex_cap_ppm=20]\n";
    return 1;
  }
  try
  {
    const int threads = argc > 2 ? StringUtils::toInt32(argv[2]) : 2;
    const int repeats = argc > 3 ? StringUtils::toInt32(argv[3]) : 9;
    const double cap = argc > 4 ? StringUtils::toDouble(argv[4]) : 20.0;
    if (threads < 2 || repeats < 1) { throw std::invalid_argument("threads must be >= 2 and repeats >= 1"); }
    omp_set_dynamic(0);
    PeakMap input;
    MzMLFile file;
    file.getOptions().setMSLevels({1});
    file.load(argv[1], input);
    input.sortSpectra(true);

    MassTraceDetection detector;
    detector.setLogType(ProgressLogger::NONE);
    Param parameters = detector.getParameters();
    parameters.setValue("max_apex_deviation_ppm", cap);
    detector.setParameters(parameters);
    auto detect = [&](int count, std::vector<MassTrace>& traces) {
      omp_set_num_threads(count);
      const auto start = std::chrono::steady_clock::now();
      detector.run(input, traces);
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    };

    std::vector<MassTrace> reference, parallel, serial;
    const double cold_serial = detect(1, reference);
    const double cold_parallel = detect(threads, parallel);
    if (! equalTraces(reference, parallel)) { throw std::runtime_error("Serial/parallel output mismatch"); }
    std::vector<double> serial_times, parallel_times;
    for (int repeat = 0; repeat < repeats; ++repeat)
    {
      if (repeat % 2 == 0)
      {
        serial_times.push_back(detect(1, serial));
        parallel_times.push_back(detect(threads, parallel));
      }
      else
      {
        parallel_times.push_back(detect(threads, parallel));
        serial_times.push_back(detect(1, serial));
      }
      if (! equalTraces(reference, serial) || ! equalTraces(reference, parallel)) { throw std::runtime_error("Output changed between repetitions"); }
    }
    double max_span = 0.0;
    const MassTrace* widest = nullptr;
    double widest_min = 0.0, widest_max = 0.0;
    for (const auto& trace : reference)
    {
      const auto [low, high] = std::minmax_element(trace.begin(), trace.end(), [](const auto& a, const auto& b) { return a.getMZ() < b.getMZ(); });
      const double span = (high->getMZ() - low->getMZ()) / trace.getCentroidMZ() * 1e6;
      if (span > max_span)
      {
        max_span = span;
        widest = &trace;
        widest_min = low->getMZ();
        widest_max = high->getMZ();
      }
    }
    const double serial_seconds = median(serial_times), parallel_seconds = median(parallel_times);
    OPENMS_LOG_INFO << std::setprecision(12) << "input=" << argv[1] << " scans=" << input.size() << " traces=" << reference.size()
                    << " cap_ppm=" << cap << " threads=" << threads << " repeats=" << repeats << " cold_serial_s=" << cold_serial
                    << " cold_parallel_s=" << cold_parallel << " serial_median_s=" << serial_seconds << " parallel_median_s=" << parallel_seconds
                    << " speedup=" << serial_seconds / parallel_seconds << " exact_equal=true\n";
    if (widest != nullptr)
    {
      OPENMS_LOG_INFO << std::setprecision(15) << "widest=" << widest->getLabel() << " points=" << widest->getSize()
                      << " rt_span_s=" << widest->getTraceLength() << " min_mz=" << widest_min << " max_mz=" << widest_max
                      << " centroid_mz=" << widest->getCentroidMZ() << " span_ppm=" << max_span << '\n';
    }
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
  return 0;
}
