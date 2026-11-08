// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAPeakMaps.h"

#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include <OpenMS/FORMAT/XIPMParquetFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  MSSpectrum ms1(double rt, std::vector<std::pair<double, double>> peaks)
  {
    MSSpectrum s;
    s.setMSLevel(1);
    s.setRT(rt);
    for (const auto& [mz, intensity] : peaks) s.push_back(Peak1D(mz, intensity));
    return s;
  }

  MSSpectrum ms2(double rt, double isolation_center, std::vector<std::pair<double, double>> peaks)
  {
    MSSpectrum s = ms1(rt, std::move(peaks));
    s.setMSLevel(2);
    Precursor precursor;
    precursor.setMZ(isolation_center);
    precursor.setIsolationWindowLowerOffset(5.0);
    precursor.setIsolationWindowUpperOffset(5.0);
    s.setPrecursors({precursor});
    return s;
  }

  struct Captured
  {
    std::string native_id;
    int ms_level = 0;
    std::vector<double> mz, rt, im, intensity;
    double rt_start = 0.0, rt_end = 0.0;
  };
}

int main()
{
  // Precursor 500.0 2+ at RT 100 with two fragments; a second precursor never inside the run's RT range.
  OpenDIA::PeakMapTarget a;
  a.precursor_id = "10";
  a.precursor_mz = 500.0;
  a.charge = 2;
  a.center_rt = 100.0;
  a.has_rt = true;
  a.fragments = {{"101", 300.0}, {"102", 400.0}};
  OpenDIA::PeakMapTarget late = a;
  late.precursor_id = "11";
  late.center_rt = 900.0;
  late.fragments = {{"111", 350.0}};

  OpenDIA::PeakMapOptions options;
  options.mz_tolerance_ppm = 20.0;
  options.rt_half_window_s = 10.0;
  options.ms1_isotopes = 1;

  std::vector<Captured> captured;
  OpenDIA::PeakMapCollector collector({a, late}, options, [&](const PeakMapExtractor::ExtractedPeakMap& map, int level)
  {
    captured.push_back({map.native_id, level, map.mz, map.rt, map.ion_mobility, map.intensity, map.rt_start, map.rt_end});
  });

  const double iso1 = 500.0 + 1.0033548 / 2.0;
  std::vector<MSSpectrum> run = {
    ms1(85.0, {{500.0, 1.0}}),                                // before the window
    ms1(95.0, {{499.9, 1.0}, {500.004, 2.0}, {iso1, 3.0}}),   // 8 ppm and isotope inside; 499.9 outside
    ms2(95.5, 500.0, {{300.003, 4.0}, {400.0, 5.0}, {450.0, 9.0}}),
    ms2(95.6, 520.0, {{300.0, 7.0}}),                         // isolation window misses the precursor
    ms1(105.0, {{500.0, 6.0}}),
    ms2(105.5, 500.0, {{300.0, 8.0}}),
    ms1(111.0, {{500.0, 99.0}}),                              // after the window: closes target 10
    ms2(111.5, 500.0, {{300.0, 99.0}}),
  };
  for (auto& spectrum : run) collector.consumeSpectrum(spectrum);
  // Target 10 was emitted when RT passed 110, before finish().
  require(captured.size() == 4, "target must be emitted once the run passes its RT window");
  collector.finish();
  collector.finish();
  require(captured.size() == 7, "finish emits the never-opened target (2 MS1 + 1 fragment) once");
  require(collector.outOfOrderSpectra() == 0, "run is in RT order");

  std::map<std::string, Captured> by_id;
  for (const auto& c : captured) by_id[c.native_id] = c;
  const Captured& mono = by_id.at("10_Precursor_i0");
  require(mono.ms_level == 1 && mono.mz.size() == 2, "monoisotope keeps 2 MS1 points inside 20 ppm and the RT window");
  require(std::abs(mono.mz[0] - 500.004) < 1e-9 && mono.rt[0] == 95.0 && mono.intensity[1] == 6.0, "MS1 point values");
  require(mono.im[0] == -1.0, "no ion mobility is recorded as -1");
  require(mono.rt_start == 90.0 && mono.rt_end == 110.0, "RT window metadata");
  require(by_id.at("10_Precursor_i1").mz.size() == 1, "first isotope point");
  const Captured& f300 = by_id.at("101");
  require(f300.ms_level == 2 && f300.mz.size() == 2 && f300.intensity[0] == 4.0 && f300.intensity[1] == 8.0,
          "fragment keeps only spectra whose isolation window holds the precursor");
  require(by_id.at("102").mz.size() == 1, "second fragment");
  require(by_id.at("11_Precursor_i0").mz.empty() && by_id.at("111").mz.empty(), "late target is written empty");

  // Ion mobility: per-peak 1/K0 is recorded and filtered by the IM window.
  {
    OpenDIA::PeakMapTarget t = a;
    t.has_im = true;
    t.precursor_im = 1.0;
    OpenDIA::PeakMapOptions im_options = options;
    im_options.im_half_window = 0.05;
    std::vector<Captured> out;
    OpenDIA::PeakMapCollector im_collector({t}, im_options, [&](const PeakMapExtractor::ExtractedPeakMap& map, int level)
    { out.push_back({map.native_id, level, map.mz, map.rt, map.ion_mobility, map.intensity, 0, 0}); });
    MSSpectrum s = ms2(100.0, 500.0, {{300.0, 1.0}, {300.001, 2.0}});
    DataArrays::FloatDataArray im;
    IMDataConverter::setIMUnit(im, DriftTimeUnit::VSSC);
    im.assign({1.02f, 1.2f});
    s.getFloatDataArrays().push_back(im);
    im_collector.consumeSpectrum(s);
    im_collector.finish();
    for (const auto& c : out)
    {
      if (c.native_id == "101")
      {
        require(c.mz.size() == 1 && std::abs(c.im[0] - 1.02) < 1e-6, "IM window keeps 1.02 and drops 1.2");
      }
    }
  }

  // Bad options fail closed.
  bool threw = false;
  try
  {
    OpenDIA::PeakMapOptions bad;
    bad.mz_tolerance_ppm = 0.0;
    OpenDIA::PeakMapCollector({a}, bad, [](const PeakMapExtractor::ExtractedPeakMap&, int) {});
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  require(threw, "zero m/z tolerance is rejected");

  // Round trip through the OpenSWATH .xipm writer and reader.
  {
    OpenSwath::LightTargetedExperiment exp;
    OpenSwath::LightCompound compound;
    compound.id = "10";
    compound.sequence = "PEPTIDEK";
    compound.charge = 2;
    compound.rt = 100.0;
    exp.compounds.push_back(compound);
    for (const auto& [name, mz] : a.fragments)
    {
      OpenSwath::LightTransition transition;
      transition.transition_name = name;
      transition.peptide_ref = "10";
      transition.precursor_mz = 500.0;
      transition.product_mz = mz;
      transition.fragment_charge = 1;
      transition.setDetectingTransition(true);
      exp.transitions.push_back(transition);
    }
    const std::string path = "OpenDIAPeakMaps_test_roundtrip.xipm";
    {
      XIPMParquetConsumer writer(path, exp);
      OpenDIA::PeakMapCollector to_file({a}, options, [&](const PeakMapExtractor::ExtractedPeakMap& map, int level)
      { writer.consumePeakMap(map, 0, "test.mzML", level); });
      for (auto& spectrum : run) to_file.consumeSpectrum(spectrum);
      to_file.finish();
      writer.finalize();
    }
    std::vector<XIPMParquetFile::XIPMPeakMap> maps;
    XIPMParquetFile(path).load(maps);
    require(maps.size() == 4, "xipm holds two MS1 and two MS2 peak maps");
    std::size_t ms1_points = 0, ms2_points = 0;
    for (const auto& map : maps)
    {
      require(map.has_precursor_id && map.precursor_id == 10 && map.modified_sequence == "PEPTIDEK", "precursor metadata");
      (map.ms_level == 1 ? ms1_points : ms2_points) += map.mz.size();
      if (map.peakmap_type == "transition" && map.transition_id == 101)
      {
        require(map.mz.size() == 2 && std::abs(map.rt[1] - 105.5) < 1e-3, "fragment points survive the round trip");
      }
    }
    require(ms1_points == 3 && ms2_points == 3, "point counts survive the round trip");
    std::remove(path.c_str());
  }

  std::cout << "OpenDIAPeakMaps_test passed\n";
  return 0;
}
