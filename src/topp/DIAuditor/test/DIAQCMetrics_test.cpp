// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, David L. Tabb $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/Precursor.h>
#include "DIAQCMetrics.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <sstream>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  MSSpectrum makeSpectrum(double rt, UInt ms_level, const std::vector<double>& intensities)
  {
    MSSpectrum s;
    s.setRT(rt);
    s.setMSLevel(ms_level);
    for (Size i = 0; i < intensities.size(); ++i) s.push_back(Peak1D(100.0 + static_cast<double>(i), intensities[i]));
    return s;
  }

  void setWindow(MSSpectrum& s, double target, double lower, double upper)
  {
    Precursor p;
    p.setMZ(target);
    p.setIsolationWindowLowerOffset(lower);
    p.setIsolationWindowUpperOffset(upper);
    s.getPrecursors().push_back(p);
  }

  void setIonMobility(MSSpectrum& s, double value, DriftTimeUnit unit)
  {
    s.setDriftTime(value);
    s.setDriftTimeUnit(unit);
  }

  /*
    Four cycles, 3 s apart, starting at 60 s. Each cycle: one MS1 spectrum (FAIMS CV -45, 3 to 6 peaks of intensity 10)
    and three MS2 spectra: 500 +- 12.5 m/z at CV -45 ("A"), 525 +- 12.5 m/z at CV -45 ("B") and 500 +- 12.5 m/z at
    CV -65 ("A'"). One more MS2 spectrum at 75 s, 800 +- 50 m/z without FAIMS ("C"), is added first.
  */
  void addRun(DIAQCMetrics& metrics)
  {
    MSSpectrum c = makeSpectrum(75.0, 2, {7.0});
    setWindow(c, 800.0, 50.0, 50.0);
    metrics.addSpectrum(c);
    for (int cycle = 0; cycle < 4; ++cycle)
    {
      const double t = 60.0 + 3.0 * cycle;
      MSSpectrum ms1 = makeSpectrum(t, 1, std::vector<double>(3 + cycle, 10.0));
      setIonMobility(ms1, -45.0, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      ms1.setMetaValue("mass resolving power", "120000");
      metrics.addSpectrum(ms1);

      MSSpectrum a = makeSpectrum(t + 0.5, 2, {5.0, 5.0});
      setWindow(a, 500.0, 12.5, 12.5);
      setIonMobility(a, -45.0, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      a.setMetaValue("total ion current", 10.0);
      metrics.addSpectrum(a);

      MSSpectrum b = makeSpectrum(t + 1.0, 2, {1.0, 2.0, 3.0});
      setWindow(b, 525.0, 12.5, 12.5);
      setIonMobility(b, -45.0, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      metrics.addSpectrum(b);

      MSSpectrum a2 = makeSpectrum(t + 1.5, 2, {0.0, 4.0});
      setWindow(a2, 500.0, 12.5, 12.5);
      setIonMobility(a2, -65.0, DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      metrics.addSpectrum(a2);
    }
  }

  const DIAQCMetrics::WindowMetrics* findWindow(const DIAQCMetrics::RunMetrics& run, double target, double faims_cv)
  {
    for (const auto& w : run.windows)
    {
      const bool same_im = std::isnan(faims_cv) ? std::isnan(w.faims_cv) : w.faims_cv == faims_cv;
      if (w.target_mz == target && same_im) return &w;
    }
    return nullptr;
  }
}

START_TEST(DIAQCMetrics, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const double NaN = std::numeric_limits<double>::quiet_NaN();

DIAQCMetrics* ptr = nullptr;
DIAQCMetrics* null_ptr = nullptr;
START_SECTION(DIAQCMetrics())
{
  ptr = new DIAQCMetrics();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->size(), 0)
}
END_SECTION

START_SECTION(~DIAQCMetrics())
{
  delete ptr;
}
END_SECTION

START_SECTION(RunMetrics compute() const)
{
  DIAQCMetrics metrics;
  addRun(metrics);
  TEST_EQUAL(metrics.size(), 17)
  const DIAQCMetrics::RunMetrics run = metrics.compute();

  TEST_EQUAL(run.ms1_count, 4)
  TEST_EQUAL(run.msn_count, 13)
  TEST_EQUAL(run.ms2_count, 13)
  TEST_REAL_SIMILAR(run.rt_min, 60.0)
  TEST_REAL_SIMILAR(run.rt_max, 75.0)

  // MS1: TIC 30, 40, 50, 60 (total 180); 25% (45) is reached at 63 s, 50% (90) at 66 s, 75% (135) at 69 s
  TEST_REAL_SIMILAR(run.ms1_total_tic, 180.0)
  TEST_REAL_SIMILAR(run.ms1_tic_quantile_rt[0], 63.0)
  TEST_REAL_SIMILAR(run.ms1_tic_quantile_rt[1], 66.0)
  TEST_REAL_SIMILAR(run.ms1_tic_quantile_rt[2], 69.0)
  TEST_REAL_SIMILAR(run.ms1_cycle_time_median, 3.0)
  TEST_REAL_SIMILAR(run.ms1_mass_resolving_power, 120000.0)
  // peak counts 3, 4, 5, 6: values at positions floor(q * 4)
  TEST_REAL_SIMILAR(run.ms1_peak_count[0], 3.0)
  TEST_REAL_SIMILAR(run.ms1_peak_count[1], 4.0)
  TEST_REAL_SIMILAR(run.ms1_peak_count[2], 5.0)
  TEST_REAL_SIMILAR(run.ms1_peak_count[3], 6.0)
  TEST_REAL_SIMILAR(run.ms1_peak_count[4], 6.0)

  // MS2 peak counts: 1 (C), 8 x 2 (A, A'), 4 x 3 (B)
  TEST_REAL_SIMILAR(run.ms2_peak_count[0], 1.0)
  TEST_REAL_SIMILAR(run.ms2_peak_count[1], 2.0)
  TEST_REAL_SIMILAR(run.ms2_peak_count[2], 2.0)
  TEST_REAL_SIMILAR(run.ms2_peak_count[3], 3.0)
  TEST_REAL_SIMILAR(run.ms2_peak_count[4], 3.0)
  TEST_REAL_SIMILAR(run.ms2_total_tic, 40.0 + 24.0 + 16.0 + 7.0)
  TEST_EQUAL(run.spectra_without_rt, 0)
  TEST_EQUAL(run.ms2_multiple_precursors, 0)
  TEST_EQUAL(run.ms2_scan_ion_mobility, 0)

  // windows in order of first acquisition, although C was added first
  TEST_EQUAL(run.window_count, 4)
  ABORT_IF(run.windows.size() != 4)
  TEST_REAL_SIMILAR(run.windows[0].target_mz, 500.0)
  TEST_REAL_SIMILAR(run.windows[0].faims_cv, -45.0)
  TEST_REAL_SIMILAR(run.windows[1].target_mz, 525.0)
  TEST_REAL_SIMILAR(run.windows[2].target_mz, 500.0)
  TEST_REAL_SIMILAR(run.windows[2].faims_cv, -65.0)
  TEST_TRUE(std::isnan(run.windows[2].ion_mobility_lower))
  TEST_REAL_SIMILAR(run.windows[3].target_mz, 800.0)
  TEST_TRUE(std::isnan(run.windows[3].faims_cv))

  const DIAQCMetrics::WindowMetrics& a = run.windows[0];
  TEST_TRUE(a.has_isolation_window)
  TEST_REAL_SIMILAR(a.lower_mz, 487.5)
  TEST_REAL_SIMILAR(a.upper_mz, 512.5)
  TEST_REAL_SIMILAR(a.width_mz, 25.0)
  TEST_EQUAL(a.spectrum_count, 4)
  TEST_REAL_SIMILAR(a.rt_min, 60.5)
  TEST_REAL_SIMILAR(a.rt_max, 69.5)
  TEST_REAL_SIMILAR(a.cycle_time_median, 3.0)
  TEST_REAL_SIMILAR(a.total_tic, 40.0)
  TEST_REAL_SIMILAR(a.tic_quantile_rt[0], 60.5)
  TEST_REAL_SIMILAR(a.tic_quantile_rt[1], 63.5)
  TEST_REAL_SIMILAR(a.tic_quantile_rt[2], 66.5)
  TEST_REAL_SIMILAR(a.peak_count[2], 2.0)
  TEST_TRUE(std::isnan(a.mass_resolving_power))

  // a window measured once has no cycle time
  const DIAQCMetrics::WindowMetrics& c = run.windows[3];
  TEST_EQUAL(c.spectrum_count, 1)
  TEST_TRUE(std::isnan(c.cycle_time_median))
  TEST_REAL_SIMILAR(c.width_mz, 100.0)

  TEST_REAL_SIMILAR(run.window_spectra_min, 1.0)
  TEST_REAL_SIMILAR(run.window_spectra_max, 4.0)
  TEST_REAL_SIMILAR(run.window_mz_min, 487.5)
  TEST_REAL_SIMILAR(run.window_mz_max, 850.0)
  TEST_REAL_SIMILAR(run.window_width_min, 25.0)
  TEST_REAL_SIMILAR(run.window_width_max, 100.0)
  TEST_REAL_SIMILAR(run.window_cycle_time_mean, 3.0)
  TEST_REAL_SIMILAR(run.window_cycle_time_median, 3.0)
  TEST_REAL_SIMILAR(run.window_half_tic_rt_min, 63.5)
  TEST_REAL_SIMILAR(run.window_half_tic_rt_max, 75.0)
  TEST_REAL_SIMILAR(run.window_total_tic_min, 7.0)
  TEST_REAL_SIMILAR(run.window_total_tic_max, 40.0)
  TEST_REAL_SIMILAR(run.window_peak_count_median_min, 1.0)
  TEST_REAL_SIMILAR(run.window_peak_count_median_max, 3.0)
  TEST_EQUAL(run.without_isolation_window.spectrum_count, 0)
  // the precursor m/z values (here the isolation targets), not the m/z range the windows cover (487.5 to 850)
  TEST_REAL_SIMILAR(run.precursor_mz_min, 500.0)
  TEST_REAL_SIMILAR(run.precursor_mz_max, 800.0)

  // no spectra
  DIAQCMetrics empty;
  const DIAQCMetrics::RunMetrics none = empty.compute();
  TEST_EQUAL(none.ms1_count, 0)
  TEST_EQUAL(none.window_count, 0)
  TEST_TRUE(std::isnan(none.rt_max))
  TEST_TRUE(std::isnan(none.ms1_cycle_time_median))
  TEST_TRUE(std::isnan(none.ms1_peak_count[2]))
  TEST_TRUE(std::isnan(none.window_mz_min))
  TEST_TRUE(std::isnan(none.precursor_mz_min))
}
END_SECTION

START_SECTION(void addSpectrum(const MSSpectrum& spectrum))
{
  // ion mobility: FAIMS voltages and ion mobility ranges of frames separate windows with the same m/z range, depending
  // on the option; the ion mobility of a single scan is a position and never does
  auto windowsWith = [](DIAQCMetrics::IonMobilityKey key)
  {
    DIAQCMetrics::Options options;
    options.ion_mobility = key;
    DIAQCMetrics metrics(options);
    addRun(metrics);
    for (double lower : {0.7, 1.0}) // two diaPASEF frames with the same m/z range
    {
      MSSpectrum s = makeSpectrum(80.0 + lower, 2, {1.0});
      setWindow(s, 900.0, 10.0, 10.0);
      s.setMetaValue("ion mobility lower limit", lower);
      s.setMetaValue("ion mobility upper limit", lower + 0.2);
      metrics.addSpectrum(s);
    }
    for (double k0 : {0.8, 1.0}) // two TIMS scans of one window
    {
      MSSpectrum s = makeSpectrum(90.0 + k0, 2, {1.0});
      setWindow(s, 950.0, 10.0, 10.0);
      setIonMobility(s, k0, DriftTimeUnit::VSSC);
      metrics.addSpectrum(s);
    }
    return metrics.compute();
  };
  const DIAQCMetrics::RunMetrics im_auto = windowsWith(DIAQCMetrics::IonMobilityKey::AUTO);
  TEST_EQUAL(im_auto.window_count, 7) // A, B, A', C, 900 @ 0.7-0.9, 900 @ 1.0-1.2, 950
  TEST_EQUAL(im_auto.ms2_scan_ion_mobility, 2)

  // a frame that also carries a single ion mobility value (e.g. its centre) is not a single scan: it has an ion
  // mobility range or an ion mobility array
  {
    DIAQCMetrics frames;
    MSSpectrum with_range = makeSpectrum(1.0, 2, {1.0});
    setWindow(with_range, 500.0, 12.5, 12.5);
    setIonMobility(with_range, 0.8, DriftTimeUnit::VSSC);
    with_range.setMetaValue("ion mobility lower limit", 0.7);
    with_range.setMetaValue("ion mobility upper limit", 0.9);
    frames.addSpectrum(with_range);
    MSSpectrum with_array = makeSpectrum(2.0, 2, {1.0, 2.0});
    setWindow(with_array, 500.0, 12.5, 12.5);
    setIonMobility(with_array, 0.8, DriftTimeUnit::VSSC);
    with_array.getFloatDataArrays().resize(1);
    with_array.getFloatDataArrays()[0].assign({0.75f, 0.85f});
    IMDataConverter::setIMUnit(with_array.getFloatDataArrays()[0], DriftTimeUnit::VSSC);
    frames.addSpectrum(with_array);
    TEST_EQUAL(frames.compute().ms2_scan_ion_mobility, 0)
  }
  ABORT_IF(im_auto.windows.size() != 7)
  TEST_REAL_SIMILAR(im_auto.windows[4].ion_mobility_lower, 0.7)
  TEST_REAL_SIMILAR(im_auto.windows[4].ion_mobility_upper, 0.9)
  TEST_REAL_SIMILAR(im_auto.windows[5].ion_mobility_lower, 1.0)
  TEST_EQUAL(im_auto.windows[6].spectrum_count, 2)
  TEST_TRUE(std::isnan(im_auto.windows[6].ion_mobility_lower))
  TEST_EQUAL(windowsWith(DIAQCMetrics::IonMobilityKey::FAIMS).window_count, 6) // the frames merge
  TEST_EQUAL(windowsWith(DIAQCMetrics::IonMobilityKey::NONE).window_count, 5)  // A and A' merge as well

  {
    DIAQCMetrics::Options options;
    options.ion_mobility = DIAQCMetrics::IonMobilityKey::NONE;
    DIAQCMetrics metrics(options);
    addRun(metrics);
    const DIAQCMetrics::RunMetrics run = metrics.compute();
    const DIAQCMetrics::WindowMetrics* merged = findWindow(run, 500.0, NaN);
    TEST_NOT_EQUAL(merged, nullptr)
    ABORT_IF(merged == nullptr)
    TEST_EQUAL(merged->spectrum_count, 8)
    // A and A' alternate at 60.5, 61.5, 63.5, 64.5, ... s: differences 1, 2, 1, 2, 1, 2, 1
    TEST_REAL_SIMILAR(merged->cycle_time_median, 1.0)
  }

  // TIC source: file value 100 vs. sum of intensities 6
  auto tic = [](DIAQCMetrics::TICSource source, bool with_file_value)
  {
    DIAQCMetrics::Options options;
    options.tic_source = source;
    DIAQCMetrics metrics(options);
    MSSpectrum s = makeSpectrum(10.0, 2, {1.0, 2.0, 3.0});
    setWindow(s, 500.0, 10.0, 10.0);
    if (with_file_value) s.setMetaValue("total ion current", 100.0);
    metrics.addSpectrum(s);
    return metrics.compute().ms2_total_tic;
  };
  TEST_REAL_SIMILAR(tic(DIAQCMetrics::TICSource::AUTO, true), 100.0)
  TEST_REAL_SIMILAR(tic(DIAQCMetrics::TICSource::AUTO, false), 6.0)
  TEST_REAL_SIMILAR(tic(DIAQCMetrics::TICSource::FILE, true), 100.0)
  TEST_REAL_SIMILAR(tic(DIAQCMetrics::TICSource::FILE, false), 0.0)
  TEST_REAL_SIMILAR(tic(DIAQCMetrics::TICSource::COMPUTED, true), 6.0)

  // the sum of intensities is accumulated in double precision (1e8 + 1 is not representable as float)
  {
    DIAQCMetrics::Options options;
    options.tic_source = DIAQCMetrics::TICSource::COMPUTED;
    DIAQCMetrics metrics(options);
    MSSpectrum s = makeSpectrum(10.0, 2, {1e8, 1.0});
    setWindow(s, 500.0, 10.0, 10.0);
    metrics.addSpectrum(s);
    TEST_TRUE(metrics.compute().ms2_total_tic == 100000001.0)
  }

  // peak counts with and without zero-intensity points (A' has one of two)
  DIAQCMetrics::Options options;
  options.peak_count = DIAQCMetrics::PeakCountMode::NONZERO;
  DIAQCMetrics nonzero(options);
  addRun(nonzero);
  const DIAQCMetrics::RunMetrics run = nonzero.compute();
  const DIAQCMetrics::WindowMetrics* a2 = findWindow(run, 500.0, -65.0);
  TEST_NOT_EQUAL(a2, nullptr)
  ABORT_IF(a2 == nullptr)
  TEST_REAL_SIMILAR(a2->peak_count[2], 1.0)

  // with a selected ion m/z, the reader moves the isolation window target into a meta value
  DIAQCMetrics selected_ion;
  for (double selected : {501.3, 507.9})
  {
    MSSpectrum s = makeSpectrum(selected, 2, {1.0});
    setWindow(s, selected, 12.5, 12.5);
    s.getPrecursors()[0].setMetaValue("isolation window target m/z", 500.0);
    selected_ion.addSpectrum(s);
  }
  const DIAQCMetrics::RunMetrics si = selected_ion.compute();
  TEST_EQUAL(si.window_count, 1)
  TEST_REAL_SIMILAR(si.windows[0].target_mz, 500.0)
  TEST_REAL_SIMILAR(si.windows[0].lower_mz, 487.5)
  TEST_REAL_SIMILAR(si.precursor_mz_min, 501.3) // the recorded precursor m/z is the selected ion
  TEST_REAL_SIMILAR(si.precursor_mz_max, 507.9)

  // MS2 spectra without an isolation window are not a DIA isolation window: they are collected in
  // without_isolation_window and left out of the window count and the window statistics
  DIAQCMetrics no_precursor;
  no_precursor.addSpectrum(makeSpectrum(1.0, 2, {1.0}));
  no_precursor.addSpectrum(makeSpectrum(2.0, 2, {1.0}));
  const DIAQCMetrics::RunMetrics np = no_precursor.compute();
  TEST_EQUAL(np.ms2_count, 2)
  TEST_EQUAL(np.window_count, 0)
  TEST_EQUAL(np.windows.size(), 0)
  TEST_EQUAL(np.windows_measured_once, 0)
  TEST_EQUAL(np.without_isolation_window.spectrum_count, 2)
  TEST_FALSE(np.without_isolation_window.has_isolation_window)
  TEST_TRUE(std::isnan(np.without_isolation_window.target_mz))
  TEST_TRUE(std::isnan(np.without_isolation_window.lower_mz))
  TEST_REAL_SIMILAR(np.without_isolation_window.cycle_time_median, 1.0)
  TEST_TRUE(std::isnan(np.window_mz_min))
  TEST_TRUE(std::isnan(np.window_spectra_max))
  TEST_TRUE(std::isnan(np.window_cycle_time_median))
  TEST_TRUE(std::isnan(np.precursor_mz_min))

  // a precursor with a selected ion but no isolation window: one group whatever the precursor m/z, no m/z range
  DIAQCMetrics no_window;
  double rt_no_window = 1.0;
  for (double mz : {600.0, 650.0, 600.0})
  {
    MSSpectrum s = makeSpectrum(rt_no_window++, 2, {1.0});
    Precursor p;
    p.setMZ(mz);
    p.setMetaValue("isolation window target m/z", 0.0); // what the reader stores for a missing isolation window
    s.getPrecursors().push_back(p);
    no_window.addSpectrum(s);
  }
  const DIAQCMetrics::RunMetrics nw = no_window.compute();
  TEST_EQUAL(nw.window_count, 0)
  TEST_EQUAL(nw.without_isolation_window.spectrum_count, 3)
  TEST_TRUE(std::isnan(nw.without_isolation_window.target_mz))
  TEST_TRUE(std::isnan(nw.window_mz_min))
  TEST_TRUE(std::isnan(nw.window_width_max))
  TEST_REAL_SIMILAR(nw.precursor_mz_min, 600.0)
  TEST_REAL_SIMILAR(nw.precursor_mz_max, 650.0)

  // a DIA run with one MS2 spectrum without an isolation window: the windows and their statistics are unchanged
  DIAQCMetrics mixed;
  addRun(mixed);
  mixed.addSpectrum(makeSpectrum(70.0, 2, {1.0})); // TIC 1, below that of every window
  const DIAQCMetrics::RunMetrics mx = mixed.compute();
  TEST_EQUAL(mx.ms2_count, 14)
  TEST_EQUAL(mx.window_count, 4)
  TEST_EQUAL(mx.windows_measured_once, 1) // C
  TEST_EQUAL(mx.without_isolation_window.spectrum_count, 1)
  TEST_REAL_SIMILAR(mx.without_isolation_window.total_tic, 1.0)
  TEST_REAL_SIMILAR(mx.window_total_tic_min, 7.0)
  TEST_REAL_SIMILAR(mx.window_half_tic_rt_max, 75.0)
  TEST_REAL_SIMILAR(mx.ms2_total_tic, 40.0 + 24.0 + 16.0 + 7.0 + 1.0) // MS2 statistics include it

  // spectra without retention time, MS3 spectra and multiplexed spectra
  DIAQCMetrics special;
  addRun(special);
  special.addSpectrum(makeSpectrum(-1.0, 1, {1.0})); // no scan start time
  MSSpectrum ms3 = makeSpectrum(61.2, 3, std::vector<double>(50, 1.0));
  setWindow(ms3, 500.0, 12.5, 12.5);
  special.addSpectrum(ms3);
  MSSpectrum msx = makeSpectrum(80.0, 2, {1.0});
  setWindow(msx, 410.0, 5.0, 5.0);
  setWindow(msx, 810.0, 5.0, 5.0);
  special.addSpectrum(msx);
  const DIAQCMetrics::RunMetrics sp = special.compute();
  TEST_EQUAL(sp.ms1_count, 5)
  TEST_EQUAL(sp.spectra_without_rt, 1)
  TEST_REAL_SIMILAR(sp.rt_min, 60.0)
  TEST_REAL_SIMILAR(sp.ms1_cycle_time_median, 3.0)
  TEST_EQUAL(sp.msn_count, 15)
  TEST_EQUAL(sp.ms2_count, 14)
  TEST_EQUAL(sp.window_count, 5) // A, B, A', C, 410 (MS3 is not in a window)
  TEST_EQUAL(sp.windows[0].spectrum_count, 4)
  TEST_REAL_SIMILAR(sp.ms2_peak_count[4], 3.0) // the 50 peaks of the MS3 spectrum are not counted
  TEST_EQUAL(sp.ms2_multiple_precursors, 1)
  TEST_EQUAL(sp.windows_measured_once, 2) // C and 410

  // a multiplexed MS2 spectrum without retention time still counts for the warning
  MSSpectrum msx_no_rt = makeSpectrum(-1.0, 2, {1.0});
  setWindow(msx_no_rt, 410.0, 5.0, 5.0);
  setWindow(msx_no_rt, 810.0, 5.0, 5.0);
  special.addSpectrum(msx_no_rt);
  TEST_EQUAL(special.compute().ms2_multiple_precursors, 2)

  // values that define a window are matched within 1e-6, as converters may write them with different last digits
  {
    DIAQCMetrics close_values;
    const std::vector<std::pair<double, double>> target_and_limit = {
      {500.0, 0.571948494756735}, {500.0 + 1e-9, 0.5719484947567351}, {500.0 - 1e-9, 0.571948494756735}, {500.001, 0.571948494756735}};
    double rt = 10.0;
    for (const auto& [target, limit] : target_and_limit)
    {
      MSSpectrum s = makeSpectrum(rt++, 2, {1.0});
      setWindow(s, target, 12.5, 12.5);
      s.setMetaValue("ion mobility lower limit", limit);
      s.setMetaValue("ion mobility upper limit", 0.9);
      close_values.addSpectrum(s);
    }
    const DIAQCMetrics::RunMetrics cv = close_values.compute();
    TEST_EQUAL(cv.window_count, 2)
    ABORT_IF(cv.windows.size() != 2)
    TEST_EQUAL(cv.windows[0].spectrum_count, 3)
    TEST_EQUAL(cv.windows_measured_once, 1)
  }

  // peak count quartiles as in DIAuditor: positions n/4, n/2 and n/4 + n/2 of the n sorted counts
  {
    DIAQCMetrics three;
    double rt = 10.0;
    for (Size peaks : {12, 5, 8})
    {
      MSSpectrum s = makeSpectrum(rt++, 2, std::vector<double>(peaks, 1.0));
      setWindow(s, 500.0, 12.5, 12.5);
      three.addSpectrum(s);
    }
    const DIAQCMetrics::PeakCountSummary pc = three.compute().windows[0].peak_count;
    TEST_REAL_SIMILAR(pc[0], 5.0)
    TEST_REAL_SIMILAR(pc[1], 5.0)
    TEST_REAL_SIMILAR(pc[2], 8.0)
    TEST_REAL_SIMILAR(pc[3], 8.0) // floor(0.75 * 3) would give 12
    TEST_REAL_SIMILAR(pc[4], 12.0)
  }
}
END_SECTION

START_SECTION(void setExperimentalSettings(const ExperimentalSettings& settings))
{
  ExperimentalSettings settings;
  settings.getInstrument().setName("Orbitrap Astral");
  settings.getInstrument().setMetaValue("instrument serial number", "SN-1");
  DateTime date;
  date.set(1, 2, 2026, 3, 4, 5);
  settings.setDateTime(date);

  DIAQCMetrics metrics;
  metrics.setExperimentalSettings(settings);
  DIAQCMetrics::RunMetrics run = metrics.compute();
  TEST_STRING_EQUAL(run.instrument, "Orbitrap Astral")
  TEST_STRING_EQUAL(run.serial_number, "SN-1")
  TEST_STRING_EQUAL(run.start_time_stamp, "2026-01-02T03:04:05")

  // the mzML reader keeps a time stamp with time zone verbatim
  settings.setMetaValue("mzml_start_time_stamp", "2026-01-02T03:04:05+01:00");
  metrics.setExperimentalSettings(settings);
  TEST_STRING_EQUAL(metrics.compute().start_time_stamp, "2026-01-02T03:04:05+01:00")

  metrics.clear();
  TEST_STRING_EQUAL(metrics.compute().instrument, "")
}
END_SECTION

START_SECTION(static void writeRunTable(const std::vector<RunMetrics>& runs, std::ostream& os))
{
  DIAQCMetrics metrics;
  addRun(metrics);
  DIAQCMetrics::RunMetrics run = metrics.compute();
  run.source_file = "run1";
  DIAQCMetrics::RunMetrics none = DIAQCMetrics().compute();
  none.source_file = "empty";

  std::ostringstream os;
  DIAQCMetrics::writeRunTable({run, none}, os);
  std::istringstream is(os.str());
  std::string header, row1, row2, extra;
  std::getline(is, header);
  std::getline(is, row1);
  std::getline(is, row2);
  TEST_FALSE(static_cast<bool>(std::getline(is, extra)))
  const std::string header_start = "SourceFile\tInstrument\tSerialNumber\tStartTimeStamp\tRTDuration\tmzMLMS1Count";
  TEST_STRING_EQUAL(header.substr(0, header_start.size()), header_start)
  TEST_TRUE(header.find("IsolationWindowWidthMax") != std::string::npos)
  TEST_TRUE(header.ends_with("MS2PkCountMax\tMedianWindowCycleTime"))
  // RTDuration 75 s = 1.25 min; resolving power without exponent
  const std::string row1_start = "run1\tNA\tNA\tNA\t1.25\t4\t13\t120000\t1.05\t1.1\t1.15\t180\t3\t3\t4\t5\t6\t6\t4\t1\t4\t487.5\t850\t25\t100\t3\t";
  TEST_STRING_EQUAL(row1.substr(0, row1_start.size()), row1_start)
  const std::string row2_start = "empty\tNA\tNA\tNA\tNA\t0\t0\tNA\tNA";
  TEST_STRING_EQUAL(row2.substr(0, row2_start.size()), row2_start)
}
END_SECTION

START_SECTION(static void writeWindowTable(const std::vector<RunMetrics>& runs, std::ostream& os))
{
  DIAQCMetrics metrics;
  addRun(metrics);
  metrics.addSpectrum(makeSpectrum(70.0, 2, {1.0})); // an MS2 spectrum without isolation window
  DIAQCMetrics::RunMetrics run = metrics.compute();
  run.source_file = "run1";

  std::ostringstream os;
  DIAQCMetrics::writeWindowTable({run}, os);
  std::istringstream is(os.str());
  std::vector<std::string> lines;
  for (std::string line; std::getline(is, line);) lines.push_back(line);
  TEST_EQUAL(lines.size(), 6)
  ABORT_IF(lines.size() != 6)
  TEST_STRING_EQUAL(lines[0], "SourceFile\tLoMZ\tHiMZ\tWidthMZ\tIonMobility\tMassResolvingPower\tMSMSCount\tRTMin\tRTMax\tCycleTimeMedian\tTIC25ileRT\tTIC50ileRT\tTIC75ileRT\tTotalTIC\tPkCountMin\tPkCount25ile\tPkCount50ile\tPkCount75ile\tPkCountMax\tIonMobilityLow\tIonMobilityHigh\tTargetMZ")
  TEST_STRING_EQUAL(lines[1], "run1\t487.5\t512.5\t25\t-45\tNA\t4\t1.00833333333\t1.15833333333\t3\t1.00833333333\t1.05833333333\t1.10833333333\t40\t2\t2\t2\t2\t2\tNA\tNA\t500")
  TEST_STRING_EQUAL(lines[4], "run1\t750\t850\t100\tNA\tNA\t1\t1.25\t1.25\tNA\t1.25\t1.25\t1.25\t7\t1\t1\t1\t1\t1\tNA\tNA\t800")
  // the MS2 spectra without isolation window follow the windows, without m/z values
  TEST_STRING_EQUAL(lines[5], "run1\tNA\tNA\tNA\tNA\tNA\t1\t1.16666666667\t1.16666666667\tNA\t1.16666666667\t1.16666666667\t1.16666666667\t1\t1\t1\t1\t1\t1\tNA\tNA\tNA")
}
END_SECTION

START_SECTION(static void writeMzQC(const std::vector<RunMetrics>& runs, std::ostream& os, const std::string& software_version, const std::string& creation_date))
{
  DIAQCMetrics metrics;
  addRun(metrics);
  DIAQCMetrics::RunMetrics run = metrics.compute();
  run.source_file = "run1";
  run.input_path = "dir with space/run#1 50%B.mzML";
  run.instrument = "Orbitrap Astral";
  DIAQCMetrics::RunMetrics none = DIAQCMetrics().compute();
  none.source_file = "empty";
  none.input_path = "empty.mzML";
  // MS2 spectra without precursor, and with a selected ion but no isolation window
  DIAQCMetrics no_precursor_metrics;
  no_precursor_metrics.addSpectrum(makeSpectrum(1.0, 2, {1.0}));
  no_precursor_metrics.addSpectrum(makeSpectrum(2.0, 2, {1.0}));
  DIAQCMetrics::RunMetrics no_precursor = no_precursor_metrics.compute();
  no_precursor.source_file = "no_precursor";
  no_precursor.input_path = "no_precursor.mzML";
  DIAQCMetrics selected_ion_metrics;
  for (double mz : {600.0, 650.0})
  {
    MSSpectrum s = makeSpectrum(mz / 100.0, 2, {1.0});
    Precursor p;
    p.setMZ(mz);
    s.getPrecursors().push_back(p);
    selected_ion_metrics.addSpectrum(s);
  }
  DIAQCMetrics::RunMetrics selected_ion = selected_ion_metrics.compute();
  selected_ion.source_file = "selected_ion";
  selected_ion.input_path = "selected_ion.mzML";

  std::ostringstream os;
  DIAQCMetrics::writeMzQC({run, none, no_precursor, selected_ion}, os, "1.2.3", "2026-01-02T03:04:05Z");
  const nlohmann::json mzqc = nlohmann::json::parse(os.str());
  TEST_STRING_EQUAL(mzqc["mzQC"]["version"].get<std::string>(), "1.0.0")
  TEST_STRING_EQUAL(mzqc["mzQC"]["creationDate"].get<std::string>(), "2026-01-02T03:04:05Z")
  TEST_EQUAL(mzqc["mzQC"]["controlledVocabularies"].size(), 1)
  const auto& run_qualities = mzqc["mzQC"]["runQualities"];
  TEST_EQUAL(run_qualities.size(), 4)
  ABORT_IF(run_qualities.size() != 4)

  const auto& metadata = run_qualities[0]["metadata"];
  TEST_STRING_EQUAL(metadata["label"].get<std::string>(), "run1")
  TEST_STRING_EQUAL(metadata["inputFiles"][0]["name"].get<std::string>(), "run#1 50%B.mzML")
  const std::string location = metadata["inputFiles"][0]["location"].get<std::string>();
  TEST_TRUE(location.starts_with("file:///"))
  TEST_TRUE(location.ends_with("/dir%20with%20space/run%231%2050%25B.mzML"))
  TEST_STRING_EQUAL(metadata["inputFiles"][0]["fileProperties"][0]["accession"].get<std::string>(), "MS:1000031")
  TEST_STRING_EQUAL(metadata["analysisSoftware"][0]["accession"].get<std::string>(), "MS:1000752")
  TEST_STRING_EQUAL(metadata["analysisSoftware"][0]["version"].get<std::string>(), "1.2.3")

  // metrics by accession; each at most once per run
  auto metricsOf = [](const nlohmann::json& run_quality)
  {
    std::map<std::string, nlohmann::json> by_accession;
    for (const auto& m : run_quality["qualityMetrics"])
    {
      TEST_EQUAL(by_accession.count(m["accession"].get<std::string>()), 0)
      by_accession[m["accession"].get<std::string>()] = m;
    }
    return by_accession;
  };
  auto m = metricsOf(run_qualities[0]);
  TEST_EQUAL(m.size(), 17)
  TEST_STRING_EQUAL(m["MS:4000194"]["name"].get<std::string>(), "DIA isolation window count")
  TEST_EQUAL(m["MS:4000194"]["value"].get<int>(), 4)
  TEST_REAL_SIMILAR(m["MS:4000067"]["value"].get<double>(), 15.0) // MS run duration: 60 s to 75 s
  TEST_STRING_EQUAL(m["MS:4000067"]["unit"]["accession"].get<std::string>(), "UO:0000010")
  TEST_TRUE(m["MS:4000061"]["value"] == nlohmann::json::array({4, 5, 6}))
  TEST_REAL_SIMILAR(m["MS:4000190"]["value"][1].get<double>(), 1.1) // minutes
  TEST_TRUE(m["MS:4000196"]["value"] == nlohmann::json::array({1, 4}))
  TEST_TRUE(m["MS:4000199"]["value"] == nlohmann::json::array({1, 3}))
  // m/z acquisition range: the precursor m/z values, not the m/z range the windows cover (487.5 to 850)
  TEST_TRUE(m["MS:4000069"]["value"] == nlohmann::json::array({500.0, 800.0}))
  TEST_TRUE(m["MS:4000190"].contains("description"))

  // undefined values are left out, never written as null
  auto e = metricsOf(run_qualities[1]);
  TEST_EQUAL(e.count("MS:4000059"), 1)
  TEST_EQUAL(e["MS:4000059"]["value"].get<int>(), 0)
  TEST_EQUAL(e.count("MS:4000067"), 0)
  TEST_EQUAL(e.count("MS:4000194"), 0)
  TEST_EQUAL(os.str().find("null"), std::string::npos)

  // MS2 spectra without isolation window: no DIA isolation window metrics, and without precursor no m/z range
  const std::vector<std::string> window_metrics = {"MS:4000193", "MS:4000194", "MS:4000195", "MS:4000196",
                                                   "MS:4000197", "MS:4000198", "MS:4000199"};
  auto np = metricsOf(run_qualities[2]);
  TEST_EQUAL(np["MS:4000060"]["value"].get<int>(), 2)
  for (const std::string& accession : window_metrics) TEST_EQUAL(np.count(accession), 0)
  TEST_EQUAL(np.count("MS:4000069"), 0)
  auto si = metricsOf(run_qualities[3]);
  for (const std::string& accession : window_metrics) TEST_EQUAL(si.count(accession), 0)
  TEST_TRUE(si["MS:4000069"]["value"] == nlohmann::json::array({600.0, 650.0}))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
