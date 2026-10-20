// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/OnDiscBrukerTimsExperiment.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <filesystem>
#include <functional>
#include <numeric>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  using Config = BrukerTimsFile::Config;

  // Whether two spectra are the same in every member but the ranges. MSSpectrum::operator== leaves out the name, the
  // ion mobility format and the ion mobility peak type; they are compared here.
  [[maybe_unused]] bool sameSpectrum(const MSSpectrum& expected, const MSSpectrum& actual)
  {
    return expected == actual &&
           expected.getName() == actual.getName() &&
           expected.getIMFormat() == actual.getIMFormat() &&
           expected.getIMPeakType() == actual.getIMPeakType();
  }

  // What differs between two spectra, in short
  [[maybe_unused]] std::string difference(const MSSpectrum& expected, const MSSpectrum& actual)
  {
    std::ostringstream d;
    d.precision(17);
    d << "expected '" << expected.getNativeID() << "' (RT " << expected.getRT() << ", MS" << expected.getMSLevel() << ", "
      << expected.size() << " peaks, " << expected.getFloatDataArrays().size() << " float arrays, IM peak type "
      << static_cast<int>(expected.getIMPeakType()) << "), got '" << actual.getNativeID() << "' (RT " << actual.getRT()
      << ", MS" << actual.getMSLevel() << ", " << actual.size() << " peaks, " << actual.getFloatDataArrays().size()
      << " float arrays, IM peak type " << static_cast<int>(actual.getIMPeakType()) << ")";
    if (expected.size() == actual.size())
    {
      for (Size i = 0; i < expected.size(); ++i)
      {
        if (expected[i] == actual[i]) continue;
        d << "; peak " << i << " differs";
        break;
      }
    }
    return d.str();
  }

  // Checks that @p exp holds the spectra of load() with the same configuration: in order, the spectra of load(), and in
  // between only MS2 spectra without peaks (those that load() leaves out). Also checks the index (native ID, retention
  // time, MS level) against the spectra, and the run metadata. Returns the number of spectra without peaks.
  [[maybe_unused]] Size checkAgainstLoad(OnDiscBrukerTimsExperiment& exp, const MSExperiment& loaded, const Config& config,
                                         const std::string& what)
  {
    Size j = 0, left_out = 0, index_mismatches = 0;
    std::string first_difference;
    for (Size i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum spectrum = exp.getSpectrum(i);
      if (spectrum.getNativeID() != exp.getSpectrumNativeID(i) || spectrum.getRT() != exp.getSpectrumRT(i) ||
          spectrum.getMSLevel() != exp.getSpectrumMSLevel(i))
      {
        ++index_mismatches;
      }
      if (j < loaded.size() && sameSpectrum(loaded[j], spectrum))
      {
        ++j;
        continue;
      }
      if (spectrum.empty() && spectrum.getMSLevel() == 2 && config.export_mode != Config::FRAME)
      {
        ++left_out;
        continue;
      }
      if (first_difference.empty())
      {
        first_difference = "spectrum " + std::to_string(i) + ": " +
                           (j < loaded.size() ? difference(loaded[j], spectrum) : std::string("load() has no more spectra"));
      }
    }
    STATUS(what << ": " << exp.size() << " spectra, " << loaded.size() << " from load(), " << left_out
                << " without peaks that load() leaves out");
    TEST_EQUAL(first_difference, "")
    TEST_EQUAL(j, loaded.size())
    TEST_EQUAL(index_mismatches, 0)
    TEST_EQUAL(exp.size(), loaded.size() + left_out)
    TEST_EQUAL(*exp.getExperimentalSettings() == static_cast<const ExperimentalSettings&>(loaded), true)
    return left_out;
  }

  [[maybe_unused]] void checkAgainstLoad(const std::string& path, const Config& config, const std::string& what)
  {
    MSExperiment loaded;
    BrukerTimsFile().load(path, loaded, config);
    OnDiscBrukerTimsExperiment exp;
    exp.openFile(path, config);
    checkAgainstLoad(exp, loaded, config, what);
  }

  // The spectra of @p exp read in the order of @p order equal those read one after the other
  [[maybe_unused]] void checkOrder(OnDiscBrukerTimsExperiment exp, const std::vector<MSSpectrum>& sequential,
                                   const std::vector<Size>& order, const std::string& what)
  {
    Size mismatches = 0;
    for (const Size i : order)
    {
      if (!sameSpectrum(sequential[i], exp.getSpectrum(i))) ++mismatches;
    }
    STATUS(what << ": " << order.size() << " spectra");
    TEST_EQUAL(mismatches, 0)
  }

  // A part of a fixture run (frames 1500-1900), with @p set applied, for the configurations that are compared
  [[maybe_unused]] Config part(const std::function<void(Config&)>& set = [](Config&) {})
  {
    Config config;
    config.frame_id_min = 1500;
    config.frame_id_max = 1900;
    set(config);
    return config;
  }
}

START_TEST(OnDiscBrukerTimsExperiment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OnDiscBrukerTimsExperiment* ptr = nullptr;
OnDiscBrukerTimsExperiment* null_ptr = nullptr;
START_SECTION(OnDiscBrukerTimsExperiment())
{
  ptr = new OnDiscBrukerTimsExperiment();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~OnDiscBrukerTimsExperiment())
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] an experiment that is not open)
{
  OnDiscBrukerTimsExperiment exp;
  TEST_EQUAL(exp.isOpen(), false)
  TEST_EQUAL(exp.getNrSpectra(), 0)
  TEST_EQUAL(exp.getNrChromatograms(), 0)
  TEST_EQUAL(exp.empty(), true)
  TEST_EQUAL(exp.getFilename(), "")
  TEST_EQUAL(exp.findFirstRT(100.0), 0)
  TEST_EQUAL(exp.getSpectrumIndexByNativeId("frame=1"), 0)
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getSpectrum(0))
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getSpectrumMetadata(0))
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getChromatogram(0))
  TEST_EXCEPTION(Exception::ElementNotFound, exp.getSpectrumByNativeId("frame=1"))
  TEST_NOT_EQUAL(exp.getExperimentalSettings().get(), nullptr)
  OnDiscBrukerTimsExperiment copy(exp);
  TEST_EQUAL(copy.isOpen(), false)
  TEST_NOT_EQUAL(exp.clone().get(), nullptr)
}
END_SECTION

START_SECTION((void openFile(const std::string& path, const Config& config = Config())))
{
  // what load() throws for a directory without Bruker data
  OnDiscBrukerTimsExperiment exp;
  const std::filesystem::path not_bruker = std::filesystem::temp_directory_path() / ("OnDiscBrukerTimsExperiment_" + File::getUniqueName() + ".d");
  std::filesystem::create_directories(not_bruker);
  MSExperiment loaded;
  TEST_EXCEPTION(Exception::FileNotReadable, BrukerTimsFile().load(not_bruker.string(), loaded))
  TEST_EXCEPTION(Exception::FileNotReadable, exp.openFile(not_bruker.string()))
  TEST_EQUAL(exp.isOpen(), false)
  std::filesystem::remove_all(not_bruker);
}
END_SECTION

#ifdef OPENTIMS_DDA_TEST_DATA
const std::string dda = OPENTIMS_DDA_TEST_DATA;
const std::string dia = OPENTIMS_DIA_TEST_DATA;
using CA = Config::CentroidAlgo;

START_SECTION([EXTRA] DDA: getSpectrum() gives the spectra of load())
{
  checkAgainstLoad(dda, Config(), "DDA, default");
  checkAgainstLoad(dda, part(), "DDA, frames 1500-1900");
  checkAgainstLoad(dda, part([](Config& c) { c.load_ms1 = false; }), "DDA, MS2 only");
  checkAgainstLoad(dda, part([](Config& c) { c.load_ms2 = false; }), "DDA, MS1 only");
  checkAgainstLoad(dda, part([](Config& c) { c.export_mode = Config::FRAME; }), "DDA, FRAME");
  checkAgainstLoad(dda, part([](Config& c) { c.export_mode = Config::FRAME; c.ms1_centroid_algo = CA::HILL_BASED; }), "DDA, FRAME, MS1 hill");
  checkAgainstLoad(dda, part([](Config& c) { c.ms1_centroid_algo = CA::HILL_BASED; c.ms2_centroid_algo = CA::HILL_BASED; c.ms2_centroid_min_hill_length = 1; }), "DDA, hill MS1 + MS2");
  checkAgainstLoad(dda, part([](Config& c) { c.ms1_centroid_im_pct = 3.0f; }), "DDA, Greedy2D MS1");
  checkAgainstLoad(dda, part([](Config& c) { c.ms1_n_neighbors = 2; c.ms1_min_support = 1; c.ms1_max_rt_distance_sec = 5.0; }), "DDA, MS1 aggregation");
  checkAgainstLoad(dda, part([](Config& c) { c.isotopic_prefilter = true; c.ms1_centroid_algo = CA::HILL_BASED; }), "DDA, isotopic prefilter");
  // DDA MS2 spectra whose ion mobility lies outside the range: load() leaves them out
  checkAgainstLoad(dda, part([](Config& c) { c.im_min = 0.9; c.im_max = 1.1; }), "DDA, ion mobility range");
  checkAgainstLoad(dda, [] { Config c; c.rt_min_sec = 120.0; c.rt_max_sec = 150.0; c.rt_selects_spectra = true; c.ms1_n_neighbors = 2; return c; }(),
                   "DDA, retention time range selecting spectra, MS1 aggregation");
  checkAgainstLoad(dda, [] { Config c; c.rt_min_sec = 120.0; c.rt_max_sec = 150.0; return c; }(), "DDA, retention time range selecting frames");
}
END_SECTION

START_SECTION([EXTRA] DIA: getSpectrum() gives the spectra of load())
{
  checkAgainstLoad(dia, Config(), "DIA, default");
  checkAgainstLoad(dia, part([](Config& c) { c.load_ms1 = false; }), "DIA, MS2 only");
  checkAgainstLoad(dia, part([](Config& c) { c.load_ms2 = false; }), "DIA, MS1 only");
  checkAgainstLoad(dia, part([](Config& c) { c.export_mode = Config::FRAME; }), "DIA, FRAME");
  checkAgainstLoad(dia, part([](Config& c) { c.export_mode = Config::SPECTRUM; }), "DIA, SPECTRUM (MS1 only: no precursors)");
  checkAgainstLoad(dia, part([](Config& c) { c.ms2_centroid_algo = CA::HILL_BASED; c.ms1_centroid_algo = CA::HILL_BASED; }), "DIA, hill MS1 + MS2");
  checkAgainstLoad(dia, part([](Config& c) { c.isotopic_prefilter = true; }), "DIA, isotopic prefilter");
  checkAgainstLoad(dia, part([](Config& c) { c.dia_ms2_n_neighbors = 1; }), "DIA, MS2 aggregation");
  checkAgainstLoad(dia, part([](Config& c) { c.dia_ms2_n_neighbors = 1; c.ms2_centroid_algo = CA::GREEDY2D; }), "DIA, MS2 aggregation, Greedy2D");
  checkAgainstLoad(dia, part([](Config& c) { c.dia_ms2_n_neighbors = 1; c.ms2_centroid_algo = CA::HILL_BASED; c.expose_hill_bounds = true; }), "DIA, MS2 aggregation, hill with bounds");
  checkAgainstLoad(dia, part([](Config& c) { c.ms1_n_neighbors = 1; c.ms1_centroid_algo = CA::HILL_BASED; }), "DIA, MS1 aggregation, hill");
  checkAgainstLoad(dia, part([](Config& c) { c.im_min = 0.9; c.im_max = 1.1; }), "DIA, ion mobility range");
  checkAgainstLoad(dia, [] { Config c; c.rt_min_sec = 120.0; c.rt_max_sec = 150.0; c.rt_selects_spectra = true; c.ms1_n_neighbors = 2; c.dia_ms2_n_neighbors = 1; return c; }(),
                   "DIA, retention time range selecting spectra, MS1 + MS2 aggregation");
}
END_SECTION

START_SECTION(([EXTRA] random and reverse order, and copies on several threads))
{
  for (const auto& [path, config, what] : {std::tuple<std::string, Config, std::string>{dda, part(), "DDA"},
                                           {dia, part(), "DIA"},
                                           {dia, part([](Config& c) { c.dia_ms2_n_neighbors = 1; c.ms1_n_neighbors = 2; }), "DIA, aggregation"},
                                           {dda, part([](Config& c) { c.ms1_n_neighbors = 2; c.ms1_min_support = 1; }), "DDA, MS1 aggregation"}})
  {
    OnDiscBrukerTimsExperiment exp;
    exp.openFile(path, config);
    std::vector<MSSpectrum> sequential;
    for (Size i = 0; i < exp.size(); ++i) sequential.push_back(exp.getSpectrum(i));

    std::vector<Size> order(exp.size());
    std::iota(order.begin(), order.end(), Size(0));
    std::reverse(order.begin(), order.end());
    checkOrder(exp, sequential, order, what + ", reverse order");
    std::shuffle(order.begin(), order.end(), std::mt19937(42));
    checkOrder(exp, sequential, order, what + ", random order");

    // clone() and the copy share the list of spectra and read with their own handle
    std::unique_ptr<IMSExperimentRandomAccess> clone = exp.clone();
    TEST_EQUAL(clone->getNrSpectra(), exp.getNrSpectra())
    TEST_EQUAL(sameSpectrum(clone->getSpectrum(exp.size() / 2), sequential[exp.size() / 2]), true)

    std::vector<MSSpectrum> parallel(exp.size());
#pragma omp parallel for firstprivate(exp) schedule(dynamic, 7)
    for (SignedSize i = static_cast<SignedSize>(order.size()) - 1; i >= 0; --i)
    {
      parallel[order[i]] = exp.getSpectrum(order[i]);
    }
    Size mismatches = 0;
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (!sameSpectrum(sequential[i], parallel[i])) ++mismatches;
    }
    STATUS(what << ", copies on several threads");
    TEST_EQUAL(mismatches, 0)
  }
}
END_SECTION

START_SECTION(MSSpectrum getSpectrumMetadata(Size index))
{
  for (const auto& [path, config] : {std::pair<std::string, Config>{dda, part()}, {dia, part()},
                                     {dia, part([](Config& c) { c.export_mode = Config::FRAME; })},
                                     {dia, part([](Config& c) { c.ms1_n_neighbors = 1; })}})
  {
    OnDiscBrukerTimsExperiment exp;
    exp.openFile(path, config);
    Size mismatches = 0;
    for (Size i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum meta = exp.getSpectrumMetadata(i);
      MSSpectrum spectrum = exp.getSpectrum(i);
      TEST_EQUAL(meta.empty(), true)
      TEST_EQUAL(meta.getFloatDataArrays().empty(), true)
      // all but the values that come from the peaks: the drift time of DDA MS2 spectra and the ion mobility peak type
      spectrum.clear(false);
      spectrum.setIMPeakType(meta.getIMPeakType());
      spectrum.setDriftTime(meta.getDriftTime());
      for (auto& precursor : spectrum.getPrecursors()) precursor.setDriftTime(IMTypes::DRIFTTIME_NOT_SET);
      if (!(meta == spectrum)) ++mismatches;
    }
    TEST_EQUAL(mismatches, 0)
  }
}
END_SECTION

START_SECTION(MSSpectrum getSpectrumByNativeId(const std::string& native_id))
{
  OnDiscBrukerTimsExperiment exp;
  exp.openFile(dda, part());
  TEST_EQUAL(exp.size() > 10, true)
  for (const Size i : {Size(0), Size(1), exp.size() / 2, exp.size() - 1})
  {
    const std::string native_id = exp.getSpectrumNativeID(i);
    TEST_EQUAL(exp.getSpectrumIndexByNativeId(native_id), i)
    TEST_EQUAL(sameSpectrum(exp.getSpectrumByNativeId(native_id), exp.getSpectrum(i)), true)
  }
  TEST_EQUAL(exp.getSpectrumIndexByNativeId("frame=0"), exp.size())
  TEST_EXCEPTION(Exception::ElementNotFound, exp.getSpectrumByNativeId("frame=0"))
}
END_SECTION

START_SECTION(Size findFirstRT(double rt) const)
{
  OnDiscBrukerTimsExperiment exp;
  exp.openFile(dia, part());
  TEST_EQUAL(exp.isSortedByRT(), true)
  TEST_EQUAL(exp.findFirstRT(-1.0), 0)
  TEST_EQUAL(exp.findFirstRT(1e9), exp.size())
  const double rt = exp.getSpectrumRT(exp.size() / 2);
  const Size first = exp.findFirstRT(rt);
  TEST_EQUAL(exp.getSpectrumRT(first), rt)
  TEST_EQUAL(first == 0 || exp.getSpectrumRT(first - 1) < rt, true)
  for (Size i = 1; i < exp.size(); ++i)
  {
    if (exp.getSpectrumRT(i - 1) > exp.getSpectrumRT(i)) TEST_EQUAL(i, 0) // sorted
  }
}
END_SECTION

START_SECTION([EXTRA] accessors of an open experiment)
{
  OnDiscBrukerTimsExperiment exp;
  const Config config = part([](Config& c) { c.load_ms1 = false; });
  exp.openFile(dia, config);
  TEST_EQUAL(exp.isOpen(), true)
  TEST_EQUAL(exp.getFilename(), dia)
  TEST_EQUAL(exp.getConfig().load_ms1, false)
  TEST_EQUAL(exp.getConfig().frame_id_min, 1500)
  TEST_EQUAL(exp.getNrChromatograms(), 0)
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getSpectrum(exp.size()))
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getSpectrumRT(exp.size()))
  TEST_EXCEPTION(Exception::IndexOverflow, exp.getChromatogram(0))
  for (Size i = 0; i < exp.size(); ++i)
  {
    if (exp.getSpectrumMSLevel(i) != 2) TEST_EQUAL(exp.getSpectrumMSLevel(i), 2)
  }
  // a copy and a moved-to experiment read the same
  OnDiscBrukerTimsExperiment copy(exp);
  OnDiscBrukerTimsExperiment assigned;
  assigned = copy;
  OnDiscBrukerTimsExperiment moved(std::move(copy));
  TEST_EQUAL(sameSpectrum(moved.getSpectrum(3), exp.getSpectrum(3)), true)
  TEST_EQUAL(sameSpectrum(assigned[4], exp[4]), true)
}
END_SECTION

START_SECTION([EXTRA] Config::calibrate with the linear m/z model is rejected)
{
  // with the MzCalibration table (and with the Bruker SDK) calibrate has no effect, and is accepted
  Config with_table = part();
  with_table.calibrate = true;
  OnDiscBrukerTimsExperiment exp;
  exp.openFile(dda, with_table);
  TEST_EQUAL(exp.isOpen(), true)
  Config linear = part();
  linear.calibrate = true;
  linear.tims_calibration_strategy = Config::TimsCalibrationStrategy::LINEAR;
  TEST_EXCEPTION(Exception::IllegalArgument, OnDiscBrukerTimsExperiment().openFile(dda, linear))
  // without MS2 there is nothing to fit
  linear.load_ms2 = false;
  checkAgainstLoad(dda, linear, "DDA, MS1 only, linear m/z model with calibrate");
  // nor with DIA data
  linear.load_ms2 = true;
  checkAgainstLoad(dia, linear, "DIA, linear m/z model with calibrate");
}
END_SECTION
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#else
#include <OpenMS/CONCEPT/ClassTest.h>
START_TEST(OnDiscBrukerTimsExperiment, "$Id$")
END_TEST
#endif // WITH_OPENTIMS
