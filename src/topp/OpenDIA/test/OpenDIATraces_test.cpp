// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIATraces.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace OpenMS;

namespace
{
  void require(bool value, const char* message)
  {
    if (!value)
    {
      throw std::runtime_error(message);
    }
  }

  template <typename Action> void rejects(Action action)
  {
    bool rejected = false;
    try
    {
      action();
    }
    catch (const Exception::BaseException&)
    {
      rejected = true;
    }
    require(rejected, "Invalid trace request was accepted.");
  }

  MSSpectrum scan(double rt, bool with_im = true)
  {
    MSSpectrum result;
    result.setMSLevel(2);
    result.setRT(rt);
    Precursor precursor;
    precursor.setMZ(500);
    precursor.setIsolationWindowLowerOffset(20);
    precursor.setIsolationWindowUpperOffset(20);
    result.setPrecursors({precursor});
    // Intentionally unsorted peaks, including wrong-IM and out-of-tolerance signal.
    result.emplace_back(200.0, 20.0f);
    result.emplace_back(100.0, 10.0f);
    result.emplace_back(100.005, 5.0f);
    result.emplace_back(100.0, 1000.0f);
    result.emplace_back(101.0, 2000.0f);
    if (with_im)
    {
      DataArrays::FloatDataArray values;
      IMDataConverter::setIMUnit(values, DriftTimeUnit::VSSC);
      values.assign({1.01f, 0.99f, 0.99f, 1.5f, 1.0f});
      result.getFloatDataArrays().push_back(std::move(values));
    }
    return result;
  }
}

int main(int argc, char** argv)
{
  require(argc == 2, "Expected unique test output prefix.");
  const std::string prefix(argv[1]);
  OpenSwath::LightTargetedExperiment library;
  OpenSwath::LightCompound compound;
  compound.id = "target";
  compound.drift_time = 1.0;
  library.compounds.push_back(compound);
  OpenSwath::LightTransition first;
  first.peptide_ref = "target";
  first.precursor_mz = 500;
  first.product_mz = 100;
  first.library_intensity = 1;
  library.transitions.push_back(first);
  first.product_mz = 200;
  first.library_intensity = 3;
  library.transitions.push_back(first);
  const std::vector<OpenDIA::TraceWindow> windows = {{"s0", "target", "peptide_family", 1, 10}};
  OpenDIA::TraceCollector collector(windows, library);
  auto spectrum = scan(10);
  collector.consumeSpectrum(spectrum);
  require(collector.hasIM(0), "Compatible mobility was not detected.");
  require(collector.window(0).sample_id == "s0" && collector.window(0).group_id == "peptide_family" &&
          collector.window(0).label == 1 && !collector.isDecoy(0), "Native metadata bridge changed candidate identity or label.");
  require(collector.fragmentMZ(0)[0] == 200 && collector.libraryIntensities(0)[0] == 3 &&
          collector.fragmentMZ(0)[1] == 100 && collector.libraryIntensities(0)[1] == 1 &&
          collector.fragmentMZ(0)[2] == 0 && collector.libraryIntensities(0)[2] == 0,
          "Native metadata bridge disagrees with tensor channel ordering or padding.");
  const auto& data = collector.intensities();
  require(data.size() == 6 * 32 * 8, "Tensor dimensions changed.");
  require(std::accumulate(data.begin(), data.begin() + 32 * 8, 0.0f) == 20.0f,
          "Top library fragment not first, or raw IM filtering wrong.");
  require(std::accumulate(data.begin() + 32 * 8, data.end(), 0.0f) == 15.0f,
          "Raw fragment intensity summation or m/z filtering wrong.");
  const std::size_t time_offset = 16 * 8;
  require(data[time_offset + 5] == 20.0f && data[32 * 8 + time_offset + 2] == 15.0f,
          "RT / IM bin coordinates or per-peak IM sorting alignment wrong.");
  const std::vector<float> data_first_scan = data;
  auto wrong_rt = scan(50);
  collector.consumeSpectrum(wrong_rt);
  require(std::accumulate(data.begin(), data.end(), 0.0f) == 35.0f, "Out-of-window RT contaminated tensor.");
  auto wrong_isolation = scan(10);
  wrong_isolation.getPrecursors()[0].setMZ(700);
  collector.consumeSpectrum(wrong_isolation);
  require(std::accumulate(data.begin(), data.end(), 0.0f) == 35.0f, "Wrong isolation contaminated tensor.");
  collector.store(prefix);
  std::ifstream numpy(prefix + ".traces.npy", std::ios::binary);
  char header[10];
  numpy.read(header, 10);
  require(std::string(header, 6) == std::string("\x93NUMPY", 6), "NumPy magic incorrect.");
  const unsigned length = static_cast<unsigned char>(header[8]) + 256 * static_cast<unsigned char>(header[9]);
  numpy.seekg(10 + length);
  std::vector<float> roundtrip(data.size());
  numpy.read(reinterpret_cast<char*>(roundtrip.data()), roundtrip.size() * sizeof(float));
  require(numpy.good() && roundtrip == data, "Binary tensor round trip changed raw intensities.");

  OpenDIA::TraceCollector no_im(windows, library);
  auto without_im = scan(10, false);
  no_im.consumeSpectrum(without_im);
  require(!no_im.hasIM(0), "Missing raw mobility reported as present.");
  // Unitless "Ion Mobility" arrays are ignored unless explicitly treated as 1/K0.
  for (const bool opt_in : {false, true})
  {
    OpenDIA::TraceCollector unitless(windows, library);
    unitless.setIonMobilityArrayAsVSSC(opt_in);
    auto unitless_scan = scan(10, false);
    DataArrays::FloatDataArray unitless_values;
    unitless_values.setName("Ion Mobility");
    unitless_values.assign({1.01f, 0.99f, 0.99f, 1.5f, 1.0f});
    unitless_scan.getFloatDataArrays().push_back(std::move(unitless_values));
    unitless.consumeSpectrum(unitless_scan);
    require(unitless.hasIM(0) == opt_in, "Unitless IM opt-in did not control mobility use.");
    require(opt_in ? unitless.intensities() == data_first_scan : !unitless.hasIM(0),
            "Opt-in unitless IM tensor differs from the annotated 1/K0 tensor.");
  }
  require(no_im.intensities()[32 * 8 + time_offset + 4] == 1015.0f,
          "Missing mobility did not preserve raw XIC intensity in central IM bin.");
  auto mixed = scan(10, false);
  collector.consumeSpectrum(mixed);
  require(!collector.hasIM(0), "Mixed IM availability was not marked missing.");
  require(collector.intensities()[time_offset + 4] == 40.0f &&
          collector.intensities()[32 * 8 + time_offset + 4] == 1030.0f,
          "Mixed-IM collapse did not preserve fragment intensities.");
  auto scalar = scan(10, false);
  scalar.setDriftTimeUnit(DriftTimeUnit::VSSC);
  scalar.setDriftTime(1.0);
  OpenDIA::TraceCollector scalar_im(windows, library);
  scalar_im.consumeSpectrum(scalar);
  require(scalar_im.hasIM(0), "Scalar inverse mobility was not detected.");

  OpenDIA::TraceCollector overflow(windows, library);
  auto huge = scan(10);
  for (auto& peak : huge)
  {
    peak.setIntensity(0.0f);
  }
  huge[1].setIntensity(std::numeric_limits<float>::max() * 0.75f);
  huge[2].setIntensity(std::numeric_limits<float>::max() * 0.75f);
  huge.getFloatDataArrays()[0][2] = 1.01f;
  overflow.consumeSpectrum(huge);
  rejects([&] { overflow.consumeSpectrum(without_im); });

  rejects([&] { OpenDIA::TraceCollector too_many(windows, library, 0.02, 15, 0.03, 0); });
  // A budget below one window's tensor rejects; exactly one window per sample fits.
  rejects([&] { OpenDIA::TraceCollector over_budget(windows, library, 0.02, 15, 0.03, 16384,
                                                    windows.size() * OpenDIA::TraceCollector::SAMPLE_SIZE * sizeof(float) - 1); });
  OpenDIA::TraceCollector exact_budget(windows, library, 0.02, 15, 0.03, 16384,
                                       windows.size() * OpenDIA::TraceCollector::SAMPLE_SIZE * sizeof(float));
  rejects([&] { OpenDIA::TraceCollector bad_tolerance(windows, library, 0); });
  auto bad_windows = windows;
  bad_windows[0].label = 2;
  rejects([&] { OpenDIA::TraceCollector invalid_label(bad_windows, library); });
  bad_windows = windows;
  bad_windows.push_back(windows[0]);
  rejects([&] { OpenDIA::TraceCollector duplicate(bad_windows, library); });
  bad_windows = windows;
  bad_windows[0].precursor_id = "absent";
  rejects([&] { OpenDIA::TraceCollector unknown(bad_windows, library); });
  library.transitions[0].setDecoy(true);
  rejects([&] { OpenDIA::TraceCollector inconsistent(windows, library); });
  library.transitions[1].setDecoy(true);
  rejects([&] { OpenDIA::TraceCollector positive_decoy(windows, library); });

  std::ofstream table(prefix + ".windows.tsv");
  table << "sample_id\tprecursor_id\tgroup_id\tlabel\trt_apex\n"
           "s0\ttarget\tpeptide_family\t1\t10\n";
  table.close();
  require(OpenDIA::readTraceWindows(prefix + ".windows.tsv", 1).front().group_id == "peptide_family",
          "Window parser lost explicit grouping.");
  rejects([&] { OpenDIA::readTraceWindows(prefix + ".windows.tsv", 0); });
  // Retain the tiny .npy/metadata fixture for a real NumPy interoperability check.
  return 0;
}
