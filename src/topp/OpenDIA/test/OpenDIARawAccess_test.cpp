// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARawAccess.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <filesystem>
#include <future>
#include <limits>
#include <stdexcept>

using namespace OpenDIA;
namespace
{
  constexpr ByteCount MIB = 1024 * 1024;
  void check(bool condition, const char* message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }

  template <typename F> void rejects(F operation)
  {
    bool rejected = false;
    try { operation(); }
    catch (const OpenMS::Exception::BaseException&) { rejected = true; }
    check(rejected, "Expected explicit failure");
  }

  struct IndexedFile
  {
    std::string path = (std::filesystem::current_path() / ("OpenDIARawAccess_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
    explicit IndexedFile(const OpenMS::MSExperiment& experiment)
    {
      OpenMS::MzMLFile writer;
      writer.getOptions().setWriteIndex(true);
      writer.store(path, experiment);
    }
    ~IndexedFile() { std::error_code error; std::filesystem::remove(path, error); }
    std::shared_ptr<OpenMS::SpectrumAccessOnDisc> source() const
    {
      auto experiment = std::make_shared<OpenMS::OnDiscMSExperiment>();
      check(experiment->openFile(path), "Cannot open indexed fixture");
      return std::make_shared<OpenMS::SpectrumAccessOnDisc>(experiment);
    }
  };

  void exact(const RawScan& scan, const OpenSwath::Spectrum& direct, const OpenMS::MSSpectrum& metadata)
  {
    check(scan.metadata() == metadata, "Full native acquisition metadata changed");
    check(scan.metadata().empty(), "Metadata contains peaks");
    check(scan.arrayCount() == direct.getDataArrays().size(), "Column count changed");
    for (OpenMS::Size i = 0; i < scan.arrayCount(); ++i)
    {
      check(scan.array(i).description == direct.getDataArrays()[i]->description, "Column description changed");
      check(scan.array(i).data == direct.getDataArrays()[i]->data, "Column values changed");
    }
  }

  void fixtures()
  {
    const std::vector<std::string> files{
      "topp/OpenSwathWorkflow_1_input.mzML", "topp/OpenSwathWorkflow_23_input.mzML",
      "class_tests/openms/data/IM_FAIMS_test.mzML"};
    for (const auto& file : files)
    {
      OpenMS::MSExperiment loaded;
      OpenMS::MzMLFile().load(std::string(OPENDIA_TEST_DATA_DIR) + "/" + file, loaded);
      IndexedFile indexed(loaded);
      auto direct = indexed.source();
      auto cache = std::make_shared<RawCache>(128 * MIB, 32 * MIB, 16 * MIB);
      RawAccess access(indexed.source(), cache);
      auto metadata = direct->getExperiment().clone();
      check(access.getNrSpectra() == direct->getNrSpectra(), "Scan count changed");
      for (size_t i = 0; i < access.getNrSpectra(); ++i)
      {
        const auto scan = access.acquire(static_cast<int>(i));
        exact(*scan, *direct->getSpectrumById(static_cast<int>(i)),
          metadata->getSpectrumMetadata(direct->getExperimentIndex(static_cast<int>(i))));
        check(access.getSpectrumMetaById(static_cast<int>(i)).RT == direct->getSpectrumMetaById(static_cast<int>(i)).RT,
          "RT metadata changed");
      }
      check(access.getSpectraByRT(0, 0) == direct->getSpectraByRT(0, 0), "Legacy RT semantics changed");
      check(cache->snapshot().budget.peak_total_bytes <= cache->snapshot().budget.usable_bytes, "Budget exceeded");
    }
  }

  OpenMS::MSExperiment smallRun()
  {
    OpenMS::MSExperiment experiment;
    for (int i = 0; i != 4; ++i)
    {
      OpenMS::MSSpectrum scan;
      scan.setNativeID("scan=" + std::to_string(i + 1));
      scan.setRT(10 + i);
      scan.setMSLevel(2);
      scan.getAcquisitionInfo().push_back(OpenMS::Acquisition());
      if (i == 1)
      {
        scan.setDriftTime(0);
        scan.setDriftTimeUnit(OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
      }
      if (i == 2)
      {
        scan.setDriftTime(0);
        scan.setDriftTimeUnit(OpenMS::DriftTimeUnit::MILLISECOND);
      }
      OpenMS::Precursor precursor;
      precursor.setMZ(500);
      precursor.setIsolationWindowLowerOffset(100);
      precursor.setIsolationWindowUpperOffset(100);
      scan.getPrecursors().push_back(precursor);
      if (i != 3)
      {
        for (int p = 0; p != 32; ++p)
        {
          OpenMS::Peak1D peak;
          peak.setMZ(400 + p);
          peak.setIntensity(p * 2);
          scan.push_back(peak);
        }
        OpenMS::MSSpectrum::FloatDataArray floats;
        floats.setName("Ion Mobility");
        floats.assign(32, 0.25f);
        scan.getFloatDataArrays().push_back(floats);
        OpenMS::MSSpectrum::IntegerDataArray integers;
        integers.setName("count");
        integers.assign(32, 7);
        scan.getIntegerDataArrays().push_back(integers);
      }
      experiment.addSpectrum(scan);
    }
    return experiment;
  }

  void queries(const IndexedFile& indexed)
  {
    auto cache = std::make_shared<RawCache>(4 * MIB, MIB, MIB);
    RawAccess access(indexed.source(), cache);
    auto direct = indexed.source();
    auto metadata = direct->getExperiment().clone();
    for (int i = 0; i != 4; ++i)
    {
      exact(*access.acquire(i), *direct->getSpectrumById(i), metadata->getSpectrumMetadata(direct->getExperimentIndex(i)));
    }
    check(access.acquire(0)->arrayCount() == 4 && access.acquire(0)->array(3).description == "count",
      "Float or integer column was lost");
    const QueryInterval all{0, 100, EdgePolicy::CLOSED};
    check(access.select(all).size() == 4, "Empty scan was lost");
    check(access.select(all, QueryInterval{450, 451, EdgePolicy::OPEN}).size() == 4, "Long isolation overlap was lost");
    check(access.select(all, QueryInterval{600, 601, EdgePolicy::OPEN}).empty(), "Open isolation edge included");
    check(access.select(all, QueryInterval{600, 601, EdgePolicy::LOWER_CLOSED}).size() == 4, "Closed isolation edge lost");
    check(access.select(all, QueryInterval{500, 500, EdgePolicy::OPEN}).empty(), "Empty isolation query matched");
    check(access.select({10, 10, EdgePolicy::CLOSED}).size() == 1, "Closed point RT query lost");
    check(access.select({10, 10, EdgePolicy::OPEN}).empty(), "Open point RT query matched");
    check(access.select(all, std::nullopt, 0) == std::vector<size_t>{1}, "Absent FAIMS became zero");
    check(access.acquire(0)->metadata().getDriftTimeUnit() == OpenMS::DriftTimeUnit::NONE, "Absent mobility changed");
    check(access.acquire(2)->metadata().getDriftTimeUnit() == OpenMS::DriftTimeUnit::MILLISECOND, "Zero IM unit lost");
    check(access.acquire(3)->array(0).data.empty(), "Empty scan gained points");
    rejects([&] { access.acquire(-1); });
    rejects([&] { access.acquire(4); });
    rejects([&] { access.acquire(0)->array(100); });
    rejects([&] { access.select({2, 1}); });
    rejects([&] { access.select({200, 300}, QueryInterval{2, 1}); });
    rejects([&] { access.select({0, std::numeric_limits<double>::infinity()}); });
    rejects([&] { access.select(all, std::nullopt, std::numeric_limits<double>::quiet_NaN()); });
    rejects([&] { access.select({0, 1, static_cast<EdgePolicy>(99)}); });
    rejects([&] { access.getSpectraByRT(10, -1); });
    rejects([&] { access.getChromatogramById(-1); });
    rejects([&] { RawCache invalid(100, 0, 101); });
    auto small = std::make_shared<RawCache>(MIB, 0, 1);
    RawAccess oversized(indexed.source(), small);
    rejects([&] { oversized.acquire(0); });
    check(small->snapshot().cache_entries == 0 && small->snapshot().budget.committed_bytes == 0 &&
      small->snapshot().budget.reserved_bytes == 0, "Oversize scan was partially published");
  }

  void leases(const IndexedFile& indexed)
  {
    auto probe = std::make_shared<RawCache>(4 * MIB, 0, MIB);
    RawAccess measure(indexed.source(), probe);
    measure.acquire(0);
    const auto bytes = probe->snapshot().budget.committed_bytes;
    auto cache = std::make_shared<RawCache>(bytes * 2, 0, bytes);
    auto access = std::make_shared<RawAccess>(indexed.source(), cache);
    auto retained = access->acquire(0);
    check(access->acquire(0) == retained && cache->snapshot().hits == 1, "Cache hit did not reuse scan");
    auto second = access->acquire(1);
    rejects([&] { access->acquire(2); });
    check(cache->snapshot().cache_entries == 0 && cache->snapshot().evictions == 2, "Pinned cache did not evict before failure");
    second.reset();
    check(access->acquire(2)->array(0).data == retained->array(0).data, "Eviction changed view");
    cache->clear();
    check(cache->snapshot().budget.committed_bytes == bytes, "Clear released live view accounting");
    retained.reset();
    check(cache->snapshot().budget.committed_bytes == 0, "Final scan owner leaked charge");
    cache = std::make_shared<RawCache>(4 * MIB, 0, MIB);
    access = std::make_shared<RawAccess>(indexed.source(), cache);
    auto mutable_scan = access->getSpectrumById(0);
    auto array_only = mutable_scan->getMZArray();
    mutable_scan->getMZArray()->data[0] = -1;
    check(access->acquire(0)->array(0).data[0] == 400, "Mutable legacy copy changed cache");
    mutable_scan.reset();
    cache->clear();
    check(cache->snapshot().budget.committed_bytes > 0, "Array-only owner lost charge");
    array_only.reset();
    check(cache->snapshot().budget.committed_bytes == 0, "Array-only owner leaked charge");
    auto container_only = access->getSpectrumById(0);
    cache->clear();
    container_only->getDataArrays().clear();
    check(cache->snapshot().budget.committed_bytes > 0, "Live spectrum container lost charge");
    container_only.reset();
    check(cache->snapshot().budget.committed_bytes == 0, "Spectrum container leaked charge");
    auto view = access->acquire(0);
    cache->cancel();
    check(cache->snapshot().budget.cancelled && view->array(0).data[0] == 400, "Cancel revoked view");
    rejects([&] { access->acquire(0); });
    cache->clear();
    rejects([&] { access->acquire(1); });
    access.reset();
    cache.reset();
    check(view->metadata().getNativeID() == "scan=1" && view->array(0).data[0] == 400, "Destroyed cache revoked view");
  }

  void concurrent(const IndexedFile& indexed)
  {
    auto cache = std::make_shared<RawCache>(32 * MIB, 0, MIB);
    RawAccess access(indexed.source(), cache);
    auto first = std::dynamic_pointer_cast<RawAccess>(access.lightClone());
    auto second = std::dynamic_pointer_cast<RawAccess>(access.lightClone());
    auto work = [](std::shared_ptr<RawAccess> clone)
    {
      for (int i = 0; i != 80; ++i)
      {
        const auto scan = clone->acquire(i % 4);
        check(scan->metadata().getRT() == 10 + i % 4, "Concurrent metadata changed");
        check(clone->getSpectrumById(i % 4)->getDataArrays().size() == scan->arrayCount(), "Concurrent columns changed");
      }
    };
    auto a = std::async(std::launch::async, work, first);
    auto b = std::async(std::launch::async, work, second);
    a.get();
    b.get();
    check(cache->snapshot().misses == 4 && cache->snapshot().cache_entries == 4, "Clones duplicated cache");
  }
}

int main()
{
  fixtures();
  IndexedFile indexed(smallRun());
  queries(indexed);
  leases(indexed);
  concurrent(indexed);
}
