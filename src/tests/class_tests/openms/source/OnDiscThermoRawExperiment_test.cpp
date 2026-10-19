// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#ifdef WITH_THERMO_RAW

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/OnDiscThermoRawExperiment.h>
#include <OpenMS/FORMAT/ThermoRawFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

#ifdef THERMO_RAW_TEST_DATA
namespace
{
  ThermoRawFile::Options centroided()
  {
    ThermoRawFile::Options options;
    options.centroid = true;
    options.checksum = false;
    return options;
  }

  MSExperiment loadAll(const ThermoRawFile::Options& options)
  {
    ThermoRawFile file;
    file.setLogType(ProgressLogger::NONE);
    file.setOptions(options);
    MSExperiment exp;
    file.load(THERMO_RAW_TEST_DATA, exp);
    return exp;
  }

  /// Number of indices i (in @p order) for which getSpectrum(i) differs from load()[i]
  Size countMismatches(OnDiscThermoRawExperiment& on_disc, const MSExperiment& loaded, const vector<Size>& order)
  {
    Size mismatches = 0;
    for (Size i : order)
    {
      if (!(on_disc.getSpectrum(i) == loaded[i])) { ++mismatches; }
    }
    return mismatches;
  }
}
#endif

START_TEST(OnDiscThermoRawExperiment, "$Id$")

OnDiscThermoRawExperiment* ptr = nullptr;
OnDiscThermoRawExperiment* null_ptr = nullptr;

START_SECTION(OnDiscThermoRawExperiment())
{
  ptr = new OnDiscThermoRawExperiment();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_FALSE(ptr->isOpen())
  TEST_EQUAL(ptr->getNrSpectra(), 0)
  TEST_EQUAL(ptr->getNrChromatograms(), 0)
  TEST_TRUE(ptr->empty())
  TEST_EXCEPTION(Exception::IndexOverflow, ptr->getSpectrum(0))
}
END_SECTION

START_SECTION(~OnDiscThermoRawExperiment())
{
  delete ptr;
}
END_SECTION

START_SECTION(void openFile(const std::string& path, const Options& options))
{
  OnDiscThermoRawExperiment exp;
  TEST_EXCEPTION(Exception::FileNotFound, exp.openFile("/definitely/not/present.raw"))
  TEST_FALSE(exp.isOpen())
}
END_SECTION

#ifdef THERMO_RAW_TEST_DATA

START_SECTION(MSSpectrum getSpectrum(Size index))
{
  // every spectrum equals the one of load() with the same options, for several selections
  struct Case
  {
    string name;
    ThermoRawFile::Options options;
  };
  vector<Case> cases;
  cases.push_back({"centroid", centroided()});
  ThermoRawFile::Options profile;
  profile.checksum = false;
  cases.push_back({"profile", profile});
  ThermoRawFile::Options ms2 = centroided();
  ms2.ms_levels = {2};
  cases.push_back({"ms2", ms2});
  ThermoRawFile::Options full = centroided();
  full.charge_data = full.noise_data = full.all_detectors = true;
  cases.push_back({"full", full});
  ThermoRawFile::Options lean = centroided();
  lean.preserve_trailers = lean.instrument_methods = false;
  cases.push_back({"lean", lean});
  ThermoRawFile::Options metadata = centroided();
  metadata.metadata_only = true;
  cases.push_back({"metadata", metadata});

  // RT ranges: the middle third of the run, also combined with an MS level
  const MSExperiment all = loadAll(centroided());
  TEST_EQUAL(all.empty(), false)
  const double rt_lo = all[0].getRT() + (all[all.size() - 1].getRT() - all[0].getRT()) / 3.0;
  const double rt_hi = all[0].getRT() + 2.0 * (all[all.size() - 1].getRT() - all[0].getRT()) / 3.0;
  ThermoRawFile::Options rt = centroided();
  rt.rt_min = rt_lo;
  rt.rt_max = rt_hi;
  cases.push_back({"rt", rt});
  ThermoRawFile::Options ms2_rt = rt;
  ms2_rt.ms_levels = {2};
  ms2_rt.skip_chromatograms = true;
  cases.push_back({"ms2_rt", ms2_rt});

  for (const Case& c : cases)
  {
    const MSExperiment loaded = loadAll(c.options);
    OnDiscThermoRawExperiment exp;
    exp.openFile(THERMO_RAW_TEST_DATA, c.options);
    TEST_TRUE(exp.isOpen())
    TEST_STRING_EQUAL(exp.getFilename(), THERMO_RAW_TEST_DATA)
    TEST_FALSE(exp.getOptions().checksum)
    TEST_EQUAL(exp.getNrSpectra(), loaded.size())
    TEST_EQUAL(exp.getNrChromatograms(), loaded.getNrChromatograms())
    TEST_TRUE(exp.isSortedByRT())
    vector<Size> order(exp.getNrSpectra());
    iota(order.begin(), order.end(), 0);
    TEST_EQUAL(countMismatches(exp, loaded, order), 0)
    for (Size i = 0; i < exp.getNrChromatograms(); ++i)
    {
      TEST_TRUE(exp.getChromatogram(i) == loaded.getChromatogram(i))
    }
    TEST_TRUE(*exp.getExperimentalSettings() == static_cast<const ExperimentalSettings&>(loaded))
    TEST_EXCEPTION(Exception::IndexOverflow, exp.getSpectrum(exp.getNrSpectra()))
    TEST_EXCEPTION(Exception::IndexOverflow, exp.getChromatogram(exp.getNrChromatograms()))
  }

  // MS2 only: the precursors (intensity estimated from the unindexed parent MS1) still match
  {
    const MSExperiment loaded = loadAll(ms2);
    OnDiscThermoRawExperiment exp;
    exp.openFile(THERMO_RAW_TEST_DATA, ms2);
    TEST_EQUAL(exp.getNrSpectra() > 0, true)
    for (Size i = 0; i < exp.getNrSpectra(); ++i)
    {
      TEST_EQUAL(exp.getSpectrumMSLevel(i), 2)
      TEST_TRUE(exp.getSpectrum(i).getPrecursors() == loaded[i].getPrecursors())
    }
  }
}
END_SECTION

START_SECTION([EXTRA] random access order)
{
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  vector<Size> order(exp.getNrSpectra());
  iota(order.begin(), order.end(), 0);
  reverse(order.begin(), order.end());
  TEST_EQUAL(countMismatches(exp, loaded, order), 0)
  shuffle(order.begin(), order.end(), mt19937(42));
  TEST_EQUAL(countMismatches(exp, loaded, order), 0)
}
END_SECTION

START_SECTION(MSSpectrum getSpectrumMetadata(Size index))
{
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  ThermoRawFile::Options headers_only = options;
  headers_only.metadata_only = true;
  const MSExperiment headers = loadAll(headers_only);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  Size mismatches = 0;
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    const MSSpectrum meta = exp.getSpectrumMetadata(i);
    if (!meta.empty() || !(meta == headers[i])) { ++mismatches; }
    const MSSpectrum peaks = exp.getSpectrumPeaks(i);
    if (peaks.size() != loaded[i].size() || !equal(peaks.begin(), peaks.end(), loaded[i].begin())) { ++mismatches; }
  }
  TEST_EQUAL(mismatches, 0)
}
END_SECTION

START_SECTION(double getSpectrumRT(Size index) const)
{
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    TEST_EQUAL(exp.getSpectrumRT(i), loaded[i].getRT())
    TEST_EQUAL(exp.getSpectrumMSLevel(i), loaded[i].getMSLevel())
    TEST_STRING_EQUAL(exp.getSpectrumNativeID(i), loaded[i].getNativeID())
  }
}
END_SECTION

START_SECTION(UInt getSpectrumMSLevel(Size index) const)
  NOT_TESTABLE // tested with getSpectrumRT
END_SECTION

START_SECTION(const std::string& getSpectrumNativeID(Size index) const)
  NOT_TESTABLE // tested with getSpectrumRT
END_SECTION

START_SECTION(int getScanNumber(Size index) const)
{
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, centroided());
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    const int scan = exp.getScanNumber(i);
    TEST_STRING_EQUAL(exp.getSpectrumNativeID(i).substr(exp.getSpectrumNativeID(i).rfind("scan=")), "scan=" + to_string(scan))
    TEST_EQUAL(exp.getSpectrumIndexByScanNumber(scan), i)
  }
  TEST_EQUAL(exp.getSpectrumIndexByScanNumber(-5), exp.getNrSpectra())
  TEST_EQUAL(exp.getSpectrumIndexByScanNumber(numeric_limits<int>::max()), exp.getNrSpectra())

  // MS2 only: the MS1 scans are not indexed
  ThermoRawFile::Options ms2 = centroided();
  ms2.ms_levels = {2};
  OnDiscThermoRawExperiment ms2_exp;
  ms2_exp.openFile(THERMO_RAW_TEST_DATA, ms2);
  const Size ms1_index = exp.getNrSpectra();
  Size first_ms1 = ms1_index;
  for (Size i = 0; i < exp.getNrSpectra(); ++i)
  {
    if (exp.getSpectrumMSLevel(i) == 1) { first_ms1 = i; break; }
  }
  TEST_NOT_EQUAL(first_ms1, ms1_index)
  TEST_EQUAL(ms2_exp.getSpectrumIndexByScanNumber(exp.getScanNumber(first_ms1)), ms2_exp.getNrSpectra())
}
END_SECTION

START_SECTION(Size getSpectrumIndexByScanNumber(int scan) const)
  NOT_TESTABLE // tested with getScanNumber
END_SECTION

START_SECTION(MSSpectrum getSpectrumByNativeId(const std::string& native_id))
{
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  const Size last = exp.getNrSpectra() - 1;
  TEST_TRUE(exp.getSpectrumByNativeId(loaded[last].getNativeID()) == loaded[last])
  TEST_EXCEPTION(Exception::ElementNotFound, exp.getSpectrumByNativeId("controllerType=0 controllerNumber=1 scan=0"))
}
END_SECTION

START_SECTION(Size findFirstRT(double rt) const)
{
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, centroided());
  const Size n = exp.getNrSpectra();
  TEST_EQUAL(exp.findFirstRT(-numeric_limits<double>::infinity()), 0)
  TEST_EQUAL(exp.findFirstRT(exp.getSpectrumRT(n - 1) + 1.0), n)
  for (Size i : {Size(0), n / 3, n / 2, n - 1})
  {
    const double rt = exp.getSpectrumRT(i);
    const Size first = exp.findFirstRT(rt);
    TEST_EQUAL(first <= i, true)
    TEST_EQUAL(exp.getSpectrumRT(first), rt)
    TEST_EQUAL(first == 0 || exp.getSpectrumRT(first - 1) < rt, true)
  }
}
END_SECTION

START_SECTION(OnDiscThermoRawExperiment(const OnDiscThermoRawExperiment& source))
{
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  OnDiscThermoRawExperiment copy(exp);
  TEST_TRUE(copy.isOpen())
  TEST_EQUAL(copy.getNrSpectra(), exp.getNrSpectra())
  vector<Size> order(copy.getNrSpectra());
  iota(order.begin(), order.end(), 0);
  TEST_EQUAL(countMismatches(copy, loaded, order), 0)

  OnDiscThermoRawExperiment assigned;
  assigned = exp;
  TEST_EQUAL(assigned.getNrSpectra(), exp.getNrSpectra())
  TEST_TRUE(assigned.getSpectrum(0) == loaded[0])

  unique_ptr<IMSExperimentRandomAccess> clone = exp.clone();
  TEST_EQUAL(clone->getNrSpectra(), exp.getNrSpectra())
  TEST_TRUE(clone->getSpectrum(exp.getNrSpectra() - 1) == loaded[loaded.size() - 1])
}
END_SECTION

START_SECTION(std::unique_ptr<IMSExperimentRandomAccess> clone() const)
{
  // one clone per thread, all reading concurrently; results equal the sequential ones
  const ThermoRawFile::Options options = centroided();
  const MSExperiment loaded = loadAll(options);
  OnDiscThermoRawExperiment exp;
  exp.openFile(THERMO_RAW_TEST_DATA, options);
  const Size n = exp.getNrSpectra();
  const int threads = 4;
  vector<Size> mismatches(threads, 0);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (int t = 0; t < threads; ++t)
  {
    unique_ptr<IMSExperimentRandomAccess> mine = exp.clone();
    vector<Size> order(n);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), mt19937(t));
    for (Size i : order)
    {
      if (!(mine->getSpectrum(i) == loaded[i])) { ++mismatches[t]; }
    }
  }
  TEST_EQUAL(accumulate(mismatches.begin(), mismatches.end(), Size(0)), 0)
}
END_SECTION

#endif // THERMO_RAW_TEST_DATA

END_TEST

#else  // !WITH_THERMO_RAW

// Placeholder main so the test target still links when the bridge is disabled.
int main(int /*argc*/, const char** /*argv*/)
{
  return 0;
}

#endif // WITH_THERMO_RAW
