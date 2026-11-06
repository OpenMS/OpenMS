// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  MSSpectrum spectrum(double rt, double isolation_mz, std::initializer_list<Peak1D> peaks)
  {
    MSSpectrum result;
    result.setMSLevel(2);
    result.setRT(rt);
    Precursor precursor;
    precursor.setMZ(isolation_mz);
    precursor.setIsolationWindowLowerOffset(50.0);
    precursor.setIsolationWindowUpperOffset(50.0);
    result.setPrecursors({precursor});
    for (const auto& peak : peaks)
    {
      result.push_back(peak);
    }
    return result;
  }

  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }
}

int main()
{
  std::vector<OpenDIA::TransitionQuery> queries = {
    {"analyte", 500.0, 100.0, 0.0, false},
    {"analyte", 500.0, 200.0, 0.0, false},
    {"other", 700.0, 300.0, 0.0, false}
  };
  OpenDIA::SearchOptions options;
  options.fragment_tolerance_da = 0.01;
  options.peak_gap_seconds = 15.0;
  options.min_transitions = 2;
  std::ostringstream output;
  OpenDIA::StreamSearch search(std::move(queries), options, output);

  MSSpectrum first = spectrum(10.0, 500.0,
                              {Peak1D(100.0, 10.0), Peak1D(100.0, 12.0), Peak1D(200.0, 20.0), Peak1D(300.0, 999.0)});
  MSSpectrum second = spectrum(20.0, 500.0, {Peak1D(100.0, 20.0), Peak1D(200.0, 30.0)});
  MSSpectrum separated = spectrum(100.0, 500.0, {Peak1D(100.0, 5.0), Peak1D(200.0, 5.0)});
  search.consumeSpectrum(first);
  search.consumeSpectrum(second);
  search.consumeSpectrum(separated);
  search.finish();

  const std::string rows = output.str();
  require(std::count(rows.begin(), rows.end(), '\n') == 2, "expected two separated peak groups");
  require(rows.find("analyte\t500\t10\t20\t20\t410\t2\t2") != std::string::npos,
          "peak-group area, apex, or unique transition support changed");
  require(rows.find("other") == std::string::npos, "precursor outside isolation window was reported");

  std::vector<OpenDIA::TransitionQuery> window_queries = { {"analyte", 500.0, 100.0, 0.0, false} };
  std::ostringstream window_output;
  OpenDIA::StreamSearch window_search(std::move(window_queries), OpenDIA::SearchOptions{}, window_output);
  MSSpectrum missing_window = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
  missing_window.getPrecursors().front().setIsolationWindowLowerOffset(0.0);
  missing_window.getPrecursors().front().setIsolationWindowUpperOffset(0.0);
  window_search.consumeSpectrum(missing_window);
  window_search.finish();
  require(window_output.str().empty(), "spectrum with no isolation-width metadata was matched");

  std::vector<OpenDIA::TransitionQuery> multi_window_queries = { {"other", 700.0, 300.0, 0.0, false} };
  std::ostringstream multi_window_output;
  OpenDIA::StreamSearch multi_window_search(std::move(multi_window_queries), OpenDIA::SearchOptions{}, multi_window_output);
  MSSpectrum multi_window = spectrum(10.0, 500.0, {Peak1D(300.0, 9.0)});
  Precursor second_window;
  second_window.setMZ(700.0);
  second_window.setIsolationWindowLowerOffset(50.0);
  second_window.setIsolationWindowUpperOffset(50.0);
  Precursor overlapping_window;
  overlapping_window.setMZ(710.0);
  overlapping_window.setIsolationWindowLowerOffset(50.0);
  overlapping_window.setIsolationWindowUpperOffset(50.0);
  multi_window.getPrecursors().push_back(second_window);
  multi_window.getPrecursors().push_back(overlapping_window);
  multi_window_search.consumeSpectrum(multi_window);
  multi_window_search.finish();
  require(multi_window_output.str().find("other\t700\t10\t10\t10\t9\t1\t1") != std::string::npos,
          "multiple or overlapping isolation windows were not handled as one scan");

  std::vector<OpenDIA::TransitionQuery> invalid_queries = {{"analyte", 500.0, 100.0, 0.0, false}};
  OpenDIA::SearchOptions invalid_options;
  invalid_options.fragment_tolerance_da = 0.0;
  std::ostringstream unused_output;
  bool rejected = false;
  try
  {
    OpenDIA::StreamSearch invalid_search(std::move(invalid_queries), invalid_options, unused_output);
  }
  catch (const Exception::InvalidParameter&)
  {
    rejected = true;
  }
  require(rejected, "non-positive m/z tolerance was accepted");
  return 0;
}
