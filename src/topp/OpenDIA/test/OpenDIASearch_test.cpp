// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
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

#ifdef _OPENMP
#include <omp.h>
#endif

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

  MSSpectrum spectrumWithIM(double rt, double isolation_mz, std::initializer_list<Peak1D> peaks,
                            std::initializer_list<float> mobilities)
  {
    MSSpectrum result = spectrum(rt, isolation_mz, peaks);
    require(peaks.size() == mobilities.size(), "test spectrum peak and IM array lengths differ");
    DataArrays::FloatDataArray im_array;
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    im_array.assign(mobilities.begin(), mobilities.end());
    result.getFloatDataArrays().push_back(std::move(im_array));
    return result;
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

  // Lane prefiltering is inclusive at both isolation-window boundaries.
  const std::vector<OpenDIA::TransitionQuery> boundary_queries = {
    {"lower_edge", 450.0, 100.0, 0.0, false},
    {"upper_edge", 550.0, 101.0, 0.0, false},
    {"below_edge", 449.999, 102.0, 0.0, false},
    {"above_edge", 550.001, 103.0, 0.0, false}
  };
  std::ostringstream boundary_output;
  OpenDIA::StreamSearch boundary_search(boundary_queries, OpenDIA::SearchOptions{}, boundary_output);
  MSSpectrum boundary_spectrum = spectrum(10.0, 500.0,
                                          {Peak1D(100.0, 10.0), Peak1D(101.0, 10.0),
                                           Peak1D(102.0, 10.0), Peak1D(103.0, 10.0)});
  boundary_search.consumeSpectrum(boundary_spectrum);
  boundary_search.finish();
  require(boundary_output.str().find("lower_edge\t450\t") != std::string::npos &&
          boundary_output.str().find("upper_edge\t550\t") != std::string::npos,
          "inclusive isolation-window boundary candidates were lost");
  require(boundary_output.str().find("below_edge") == std::string::npos &&
          boundary_output.str().find("above_edge") == std::string::npos,
          "candidate outside isolation-window boundaries was reported");

  // Interleaved lanes and skipped lanes preserve each candidate's ordered trace.
  std::vector<OpenDIA::TransitionQuery> cycle_queries = {
    {"lane_a", 500.0, 100.0, 0.0, false},
    {"lane_b", 700.0, 100.0, 0.0, false}
  };
  OpenDIA::SearchOptions cycle_options;
  cycle_options.fragment_tolerance_da = 0.01;
  cycle_options.peak_gap_seconds = 30.0;
  std::ostringstream cycle_output;
  OpenDIA::StreamSearch cycle_search(cycle_queries, cycle_options, cycle_output);
  MSSpectrum lane_a_first = spectrum(10.0, 500.0, {Peak1D(100.0, 2.0)});
  MSSpectrum lane_b_first = spectrum(20.0, 700.0, {Peak1D(100.0, 5.0)});
  MSSpectrum lane_a_second = spectrum(30.0, 500.0, {Peak1D(100.0, 4.0)});
  cycle_search.consumeSpectrum(lane_a_first);
  cycle_search.consumeSpectrum(lane_b_first);
  cycle_search.consumeSpectrum(lane_a_second);
  cycle_search.finish();
  require(cycle_output.str().find("lane_a\t500\t10\t30\t30\t60\t1\t2") != std::string::npos,
          "grid cycle boundary or sparse missing-lane trace changed");
  require(cycle_output.str().find("lane_b\t700\t20\t20\t20\t5\t1\t1") != std::string::npos,
          "window-lane query was assigned to the wrong grid row");

  const std::vector<OpenDIA::TransitionQuery> im_queries = {
    {"target", 500.0, 100.0, 0.0, false, 0.80, true},
    {"mobility_interferer", 500.0, 100.0, 0.0, false, 0.95, true}
  };
  OpenDIA::SearchOptions im_options;
  im_options.fragment_tolerance_da = 0.01;
  im_options.ion_mobility_tolerance_vssc = 0.02;
  im_options.peak_gap_seconds = 15.0;
  std::ostringstream im_output;
  OpenDIA::StreamSearch im_search(im_queries, im_options, im_output);
  MSSpectrum im_first = spectrumWithIM(10.0, 500.0, {Peak1D(100.0, 10.0)}, {0.801f});
  MSSpectrum im_second = spectrumWithIM(20.0, 500.0, {Peak1D(100.0, 20.0)}, {0.799f});
  im_search.consumeSpectrum(im_first);
  im_search.consumeSpectrum(im_second);
  im_search.finish();
  require(im_search.ionMobilityFilteringUsed(), "compatible library and raw mobility did not activate filtering");
  require(im_output.str().find("target\t500\t") != std::string::npos,
          "mobility-matched candidate was lost");
  require(im_output.str().find("mobility_interferer") == std::string::npos,
          "mobility-separated interferer was not rejected");

  // Sorting the cell by m/z must keep nonuniform mobility values paired with
  // their original centroids in the compact IM sidecar.
  std::ostringstream sorted_im_output;
  OpenDIA::StreamSearch sorted_im_search(im_queries, im_options, sorted_im_output);
  MSSpectrum unsorted_im = spectrumWithIM(10.0, 500.0,
                                          {Peak1D(100.005, 50.0), Peak1D(100.0, 10.0)},
                                          {0.95f, 0.80f});
  sorted_im_search.consumeSpectrum(unsorted_im);
  sorted_im_search.finish();
  require(sorted_im_output.str().find("target\t500\t10\t10\t10\t10\t1\t1") != std::string::npos,
          "m/z sorting detached the target mobility from its centroid");
  require(sorted_im_output.str().find("mobility_interferer\t500\t10\t10\t10\t50\t1\t1") != std::string::npos,
          "m/z sorting detached the interferer mobility from its centroid");

  // Missing raw mobility automatically falls back to the m/z-only path.
  std::ostringstream no_im_output;
  OpenDIA::StreamSearch no_im_search(im_queries, im_options, no_im_output);
  MSSpectrum no_im_spectrum = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
  no_im_search.consumeSpectrum(no_im_spectrum);
  no_im_search.finish();
  require(!no_im_search.ionMobilityFilteringUsed(), "mobility filtering reported active without raw mobility");
  require(no_im_output.str().find("target\t500\t") != std::string::npos &&
          no_im_output.str().find("mobility_interferer\t500\t") != std::string::npos,
          "m/z-only fallback changed when raw mobility was absent");

  // Raw mobility without library coordinates also keeps m/z-only matching.
  std::vector<OpenDIA::TransitionQuery> no_library_im_queries = {
    {"unannotated", 500.0, 100.0, 0.0, false}
  };
  std::ostringstream no_library_im_output;
  OpenDIA::StreamSearch no_library_im_search(no_library_im_queries, im_options, no_library_im_output);
  MSSpectrum unannotated_spectrum = spectrumWithIM(10.0, 500.0, {Peak1D(100.0, 10.0)}, {0.95f});
  no_library_im_search.consumeSpectrum(unannotated_spectrum);
  no_library_im_search.finish();
  require(!no_library_im_search.ionMobilityFilteringUsed(), "filtering reported active without library mobility");
  require(no_library_im_output.str().find("unannotated\t500\t") != std::string::npos,
          "m/z-only fallback changed when library mobility was absent");

  // A scalar VSSC drift time is also detected without a CLI mode switch.
  MSSpectrum scalar_im = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
  scalar_im.setDriftTime(0.95);
  scalar_im.setDriftTimeUnit(DriftTimeUnit::VSSC);
  std::ostringstream scalar_im_output;
  OpenDIA::StreamSearch scalar_im_search(im_queries, im_options, scalar_im_output);
  scalar_im_search.consumeSpectrum(scalar_im);
  scalar_im_search.finish();
  require(scalar_im_search.ionMobilityFilteringUsed(), "scalar VSSC mobility was not detected");
  require(scalar_im_output.str().find("mobility_interferer\t500\t") != std::string::npos &&
          scalar_im_output.str().find("target\t500\t") == std::string::npos,
          "scalar VSSC mobility did not select the matching candidate");

#ifdef _OPENMP
  // More than one candidate block and more than one scan batch must produce
  // byte-identical output when the worker count changes.
  std::vector<std::string> parallel_ids;
  parallel_ids.reserve(300);
  for (std::size_t candidate = 0; candidate < 300; ++candidate)
  {
    parallel_ids.push_back("parallel_" + std::to_string(candidate));
  }
  const auto run_candidate_blocks = [&parallel_ids](int thread_count)
  {
    omp_set_num_threads(thread_count);
    std::vector<OpenDIA::TransitionQuery> candidate_queries;
    candidate_queries.reserve(parallel_ids.size() * 2);
    for (std::size_t candidate = 0; candidate < parallel_ids.size(); ++candidate)
    {
      const double precursor_mz = 450.0 + static_cast<double>(candidate) * 0.3;
      for (std::size_t transition = 0; transition < 2; ++transition)
      {
        const double fragment_mz = 100.0 + static_cast<double>(candidate % 256) * 0.02 +
                                   static_cast<double>(transition) * 0.005;
        candidate_queries.push_back({parallel_ids[candidate], precursor_mz, fragment_mz, 0.0, false});
      }
    }

    OpenDIA::SearchOptions candidate_options;
    candidate_options.fragment_tolerance_da = 0.001;
    candidate_options.peak_gap_seconds = 15.0;
    candidate_options.min_transitions = 2;
    std::ostringstream candidate_output;
    OpenDIA::StreamSearch candidate_search(std::move(candidate_queries), candidate_options, candidate_output);
    for (std::size_t scan = 0; scan < 40; ++scan)
    {
      const double rt = scan < 20 ? 10.0 + static_cast<double>(scan) : 110.0 + static_cast<double>(scan);
      MSSpectrum scan_spectrum = spectrum(rt, 500.0, {});
      for (std::size_t candidate = 0; candidate < parallel_ids.size(); ++candidate)
      {
        const double fragment_mz = 100.0 + static_cast<double>(candidate % 256) * 0.02;
        scan_spectrum.push_back(Peak1D(fragment_mz, 1.0));
        scan_spectrum.push_back(Peak1D(fragment_mz + 0.005, 2.0));
      }
      candidate_search.consumeSpectrum(scan_spectrum);
    }
    candidate_search.finish();
    return candidate_output.str();
  };

  const int original_threads = omp_get_max_threads();
  const std::string serial_output = run_candidate_blocks(1);
  const std::string parallel_output = run_candidate_blocks(4);
  omp_set_num_threads(original_threads);
  require(serial_output == parallel_output, "candidate-block extraction changed exact TSV output");
  require(std::count(serial_output.begin(), serial_output.end(), '\n') == 600,
          "candidate-block extraction lost completed groups across scan batches");
  require(serial_output.find("parallel_0\t450\t10\t10\t29\t57\t2\t20") != std::string::npos &&
          serial_output.find("parallel_0\t450\t130\t130\t149\t57\t2\t20") != std::string::npos,
          "candidate-block extraction changed completed or active peak-group values");
#endif

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

  invalid_options = OpenDIA::SearchOptions{};
  invalid_options.ion_mobility_tolerance_vssc = 0.0;
  rejected = false;
  try
  {
    OpenDIA::StreamSearch invalid_search(std::move(invalid_queries), invalid_options, unused_output);
  }
  catch (const Exception::InvalidParameter&)
  {
    rejected = true;
  }
  require(rejected, "non-positive ion-mobility tolerance was accepted");
  return 0;
}
