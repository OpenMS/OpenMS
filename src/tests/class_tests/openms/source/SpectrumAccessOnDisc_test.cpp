// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOnDisc.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#ifdef WITH_THERMO_RAW
  #include <OpenMS/FORMAT/OnDiscThermoRawExperiment.h>
  #include <OpenMS/FORMAT/ThermoRawFile.h>
#endif

#include <memory>
#include <numeric>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  /// Whether spectrum @p id of @p a and @p b have the same meta data and m/z and intensity arrays
  bool sameSpectrum(OpenSwath::ISpectrumAccess& a, int id_a, OpenSwath::ISpectrumAccess& b, int id_b)
  {
    const OpenSwath::SpectrumMeta meta_a = a.getSpectrumMetaById(id_a);
    const OpenSwath::SpectrumMeta meta_b = b.getSpectrumMetaById(id_b);
    if (meta_a.RT != meta_b.RT || meta_a.ms_level != meta_b.ms_level) { return false; }
    const OpenSwath::SpectrumPtr sa = a.getSpectrumById(id_a);
    const OpenSwath::SpectrumPtr sb = b.getSpectrumById(id_b);
    return sa->getMZArray()->data == sb->getMZArray()->data &&
           sa->getIntensityArray()->data == sb->getIntensityArray()->data;
  }

  /// Number of spectra of @p on_disc that differ from those of @p reference (same id)
  Size countMismatches(OpenSwath::ISpectrumAccess& on_disc, OpenSwath::ISpectrumAccess& reference)
  {
    if (on_disc.getNrSpectra() != reference.getNrSpectra()) { return numeric_limits<Size>::max(); }
    Size mismatches = 0;
    for (int i = 0; i < static_cast<int>(on_disc.getNrSpectra()); ++i)
    {
      if (!sameSpectrum(on_disc, i, reference, i)) { ++mismatches; }
    }
    return mismatches;
  }
}

START_TEST(SpectrumAccessOnDisc, "$Id$")

const string mzml = OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML");
auto on_disc = make_shared<OnDiscMSExperiment>();
on_disc->openFile(mzml);
auto in_memory = make_shared<MSExperiment>();
MzMLFile().load(mzml, *in_memory);

SpectrumAccessOnDisc* ptr = nullptr;
SpectrumAccessOnDisc* null_ptr = nullptr;

START_SECTION(explicit SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment))
{
  ptr = new SpectrumAccessOnDisc(on_disc);
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getNrSpectra(), in_memory->size())
  TEST_EQUAL(ptr->getNrChromatograms(), in_memory->getNrChromatograms())
  const std::shared_ptr<IMSExperimentRandomAccess> no_experiment;
  TEST_EXCEPTION(Exception::NullPointer, SpectrumAccessOnDisc{no_experiment})
}
END_SECTION

START_SECTION(~SpectrumAccessOnDisc())
{
  delete ptr;
}
END_SECTION

START_SECTION(OpenSwath::SpectrumPtr getSpectrumById(int id))
{
  SpectrumAccessOnDisc access(on_disc);
  SpectrumAccessOpenMS reference(in_memory);
  TEST_EQUAL(access.getNrSpectra() > 0, true)
  TEST_EQUAL(countMismatches(access, reference), 0)
  TEST_EXCEPTION(Exception::IndexOverflow, access.getSpectrumById(-1))
  TEST_EXCEPTION(Exception::IndexOverflow, access.getSpectrumById(static_cast<int>(access.getNrSpectra())))
}
END_SECTION

START_SECTION(OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const)
{
  SpectrumAccessOnDisc access(on_disc);
  for (Size i = 0; i < access.getNrSpectra(); ++i)
  {
    const OpenSwath::SpectrumMeta meta = access.getSpectrumMetaById(static_cast<int>(i));
    TEST_REAL_SIMILAR(meta.RT, (*in_memory)[i].getRT())
    TEST_EQUAL(meta.ms_level, (*in_memory)[i].getMSLevel())
  }
}
END_SECTION

START_SECTION(std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const)
{
  SpectrumAccessOnDisc access(on_disc);
  SpectrumAccessOpenMS reference(in_memory);
  for (Size i = 0; i < access.getNrSpectra(); ++i)
  {
    const double rt = (*in_memory)[i].getRT();
    for (double delta : {0.0, 0.5, 1e6})
    {
      TEST_TRUE(access.getSpectraByRT(rt, delta) == reference.getSpectraByRT(rt, delta))
    }
  }
  TEST_TRUE(access.getSpectraByRT(-1e9, 1.0) == reference.getSpectraByRT(-1e9, 1.0))
  TEST_TRUE(access.getSpectraByRT(1e9, 1.0) == reference.getSpectraByRT(1e9, 1.0))
}
END_SECTION

START_SECTION(size_t getNrSpectra() const)
  NOT_TESTABLE // tested in the constructor
END_SECTION

START_SECTION(OpenSwath::ChromatogramPtr getChromatogramById(int id))
{
  SpectrumAccessOnDisc access(on_disc);
  SpectrumAccessOpenMS reference(in_memory);
  TEST_EQUAL(access.getNrChromatograms() > 0, true)
  for (int i = 0; i < static_cast<int>(access.getNrChromatograms()); ++i)
  {
    TEST_TRUE(access.getChromatogramById(i)->getTimeArray()->data == reference.getChromatogramById(i)->getTimeArray()->data)
    TEST_TRUE(access.getChromatogramById(i)->getIntensityArray()->data == reference.getChromatogramById(i)->getIntensityArray()->data)
    TEST_STRING_EQUAL(access.getChromatogramNativeID(i), reference.getChromatogramNativeID(i))
  }
}
END_SECTION

START_SECTION(size_t getNrChromatograms() const)
  NOT_TESTABLE // tested in the constructor
END_SECTION

START_SECTION(std::string getChromatogramNativeID(int id) const)
  NOT_TESTABLE // tested with getChromatogramById
END_SECTION

START_SECTION(SpectrumAccessOnDisc(std::shared_ptr<IMSExperimentRandomAccess> experiment, std::vector<Size> indices))
{
  // a view of the spectra in reverse order is ordered by RT again; ids map to the experiment
  vector<Size> indices(in_memory->size());
  iota(indices.begin(), indices.end(), 0);
  vector<Size> reversed(indices.rbegin(), indices.rend());
  SpectrumAccessOnDisc view(on_disc, reversed);
  SpectrumAccessOpenMS reference(in_memory);
  TEST_EQUAL(countMismatches(view, reference), 0)
  for (Size i = 0; i < view.getNrSpectra(); ++i)
  {
    TEST_EQUAL(view.getExperimentIndex(static_cast<int>(i)), i)
  }
  // a subset with only the last spectrum
  SpectrumAccessOnDisc last(on_disc, {in_memory->size() - 1});
  TEST_EQUAL(last.getNrSpectra(), 1)
  TEST_EQUAL(last.getExperimentIndex(0), in_memory->size() - 1)
  TEST_TRUE(sameSpectrum(last, 0, reference, static_cast<int>(in_memory->size() - 1)))
  TEST_EQUAL(last.getSpectraByRT((*in_memory)[0].getRT(), 0.0).size() <= 1, true)
  TEST_EXCEPTION(Exception::IndexOverflow, SpectrumAccessOnDisc(on_disc, {in_memory->size()}))
}
END_SECTION

START_SECTION(Size getExperimentIndex(int id) const)
  NOT_TESTABLE // tested with the subset constructor
END_SECTION

START_SECTION(const IMSExperimentRandomAccess& getExperiment() const)
{
  SpectrumAccessOnDisc access(on_disc);
  TEST_EQUAL(&access.getExperiment(), on_disc.get())
  SpectrumAccessOnDisc copy(access);
  TEST_NOT_EQUAL(&copy.getExperiment(), on_disc.get()) // the copy reads through its own clone
}
END_SECTION

START_SECTION(SpectrumAccessOnDisc(const SpectrumAccessOnDisc& rhs))
{
  SpectrumAccessOnDisc access(on_disc, {in_memory->size() - 1, 0});
  SpectrumAccessOnDisc copy(access);
  TEST_EQUAL(copy.getNrSpectra(), access.getNrSpectra())
  TEST_EQUAL(countMismatches(copy, access), 0)
}
END_SECTION

START_SECTION(std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const)
{
  // one clone per thread, reading concurrently, gives the sequential results
  SpectrumAccessOnDisc access(on_disc);
  SpectrumAccessOpenMS reference(in_memory);
  const int threads = 4;
  vector<Size> mismatches(threads, 0);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
  for (int t = 0; t < threads; ++t)
  {
    shared_ptr<OpenSwath::ISpectrumAccess> mine = access.lightClone();
    SpectrumAccessOpenMS mine_reference(in_memory);
    for (int repeat = 0; repeat < 20; ++repeat)
    {
      mismatches[t] += countMismatches(*mine, mine_reference);
    }
  }
  TEST_EQUAL(accumulate(mismatches.begin(), mismatches.end(), Size(0)), 0)
}
END_SECTION

#if defined(WITH_THERMO_RAW) && defined(THERMO_RAW_TEST_DATA)
START_SECTION([EXTRA] Thermo RAW views)
{
  ThermoRawFile::Options options;
  options.centroid = true;
  options.checksum = false;
  auto raw = make_shared<OnDiscThermoRawExperiment>();
  raw->openFile(THERMO_RAW_TEST_DATA, options);
  auto loaded = make_shared<MSExperiment>();
  ThermoRawFile file;
  file.setLogType(ProgressLogger::NONE);
  file.setOptions(options);
  file.load(THERMO_RAW_TEST_DATA, *loaded);

  SpectrumAccessOnDisc all(raw);
  SpectrumAccessOpenMS reference(loaded);
  TEST_EQUAL(countMismatches(all, reference), 0)

  // per-MS-level views, as SwathFile builds per-window views
  for (UInt level : {1u, 2u})
  {
    vector<Size> indices;
    auto subset = make_shared<MSExperiment>();
    for (Size i = 0; i < raw->getNrSpectra(); ++i)
    {
      if (raw->getSpectrumMSLevel(i) == level)
      {
        indices.push_back(i);
        subset->addSpectrum((*loaded)[i]);
      }
    }
    SpectrumAccessOnDisc view(raw, indices);
    SpectrumAccessOpenMS subset_reference(subset);
    TEST_EQUAL(view.getNrSpectra(), indices.size())
    TEST_EQUAL(countMismatches(view, subset_reference), 0)
    if (!subset->empty())
    {
      const double rt = (*subset)[subset->size() / 2].getRT();
      TEST_TRUE(view.getSpectraByRT(rt, 5.0) == subset_reference.getSpectraByRT(rt, 5.0))
    }
    const int threads = 4;
    vector<Size> mismatches(threads, 0);
#pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (int t = 0; t < threads; ++t)
    {
      shared_ptr<OpenSwath::ISpectrumAccess> mine = view.lightClone();
      SpectrumAccessOpenMS mine_reference(subset);
      mismatches[t] = countMismatches(*mine, mine_reference);
    }
    TEST_EQUAL(accumulate(mismatches.begin(), mismatches.end(), Size(0)), 0)
  }
}
END_SECTION
#endif

END_TEST
