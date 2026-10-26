// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARawAccess.h"
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/SwathFile.h>
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
      size_t points = 0;
      ByteCount largest = 0;
      check(access.getNrSpectra() == direct->getNrSpectra(), "Scan count changed");
      for (size_t i = 0; i < access.getNrSpectra(); ++i)
      {
        cache->clear();
        const auto scan = access.acquire(static_cast<int>(i));
        largest = std::max(largest, cache->snapshot().budget.committed_bytes);
        points += scan->array(0).data.size();
        exact(*scan, *direct->getSpectrumById(static_cast<int>(i)),
          metadata->getSpectrumMetadata(direct->getExperimentIndex(static_cast<int>(i))));
        check(access.getSpectrumMetaById(static_cast<int>(i)).RT == direct->getSpectrumMetaById(static_cast<int>(i)).RT,
          "RT metadata changed");
      }
      check(access.getSpectraByRT(0, 0) == direct->getSpectraByRT(0, 0), "Legacy RT semantics changed");
      check(cache->snapshot().budget.peak_total_bytes <= cache->snapshot().budget.usable_bytes, "Budget exceeded");
      // Exercise the public extractor with both scan copies and shared native clones.
      using Extractor = OpenMS::ChromatogramExtractorAlgorithm;
      const auto sample = direct->getSpectrumById(0);
      check(!sample->getMZArray()->data.empty(), "Expected a nonempty reference scan");
      const auto mobility = sample->getDriftTimeArray();
      const auto positive = std::find_if(sample->getIntensityArray()->data.begin(), sample->getIntensityArray()->data.end(),
        [](double value) { return value > 0; });
      check(positive != sample->getIntensityArray()->data.end(), "Reference scan has no signal");
      const auto peak = static_cast<size_t>(positive - sample->getIntensityArray()->data.begin());
      std::vector<RawCacheSnapshot> baselines;
      for (const size_t query_count : {size_t{1}, size_t{16}, size_t{4096}})
      {
        size_t traversal = 0;
        std::vector<Extractor::ExtractionCoordinates> coordinates(query_count);
        for (size_t q = 0; q < query_count; ++q)
        {
          coordinates[q].mz = sample->getMZArray()->data[peak] + q * 0.01;
          coordinates[q].ion_mobility = mobility && !mobility->data.empty() ? mobility->data[peak] : 0;
          coordinates[q].rt_end = -1;
        }
        const auto extract = [&](const OpenSwath::SpectrumAccessPtr& source, double im_window)
        {
          std::vector<OpenSwath::ChromatogramPtr> output;
          for (size_t q = 0; q < query_count; ++q) { output.push_back(std::make_shared<OpenSwath::Chromatogram>()); }
          Extractor().extractChromatograms(source, output, coordinates, 0.05, false, im_window, "tophat");
          return output;
        };
        for (const double im_window : {-1.0, mobility ? 0.05 : -1.0})
        {
          const auto reference = extract(direct, im_window);
          check(std::any_of(reference.front()->getIntensityArray()->data.begin(), reference.front()->getIntensityArray()->data.end(),
            [](double value) { return value > 0; }), "Reference extraction has no signal");
          for (const auto bytes : {2 * (largest + 4096), 64 * MIB})
          {
            auto bounded = std::make_shared<RawCache>(bytes, 0, largest + 4096);
            auto tiled = std::make_shared<RawAccess>(indexed.source(), bounded);
            const auto actual = extract(tiled->lightClone(), im_window);
            size_t output_elements = 0;
            for (size_t q = 0; q < query_count; ++q)
            {
              check(actual[q]->getTimeArray()->data == reference[q]->getTimeArray()->data, "Extraction RT changed");
              check(actual[q]->getIntensityArray()->data == reference[q]->getIntensityArray()->data, "Extraction intensity changed");
              output_elements += actual[q]->getTimeArray()->data.size() + actual[q]->getIntensityArray()->data.size();
            }
            const auto measured = bounded->snapshot();
            OPENMS_LOG_INFO << "WP03 " << file << " queries=" << query_count << " limit=" << bytes
              << " im_window=" << im_window << " private_peak=" << measured.budget.peak_total_bytes
              << " hits=" << measured.hits << " misses=" << measured.misses
              << " output_elements=" << output_elements << '\n';
            if (query_count == 1) { baselines.push_back(measured); }
            const auto& baseline = baselines[traversal++];
            check(measured.budget.peak_total_bytes == baseline.budget.peak_total_bytes,
              "Raw cache peak changed with query count at fixed traversal and limit");
            check(measured.misses == baseline.misses,
              "Native scan misses changed with query count at fixed traversal and limit");
            check(measured.hits == baseline.hits,
              "Raw cache hits changed with query count at fixed traversal and limit");
            check(measured.budget.peak_total_bytes <= bytes, "Public extraction exceeded its private cap");
            bounded->clear();
            check(bounded->snapshot().budget.committed_bytes == 0, "Extraction retained scan charges");
          }
        }
      }
      OPENMS_LOG_INFO << "WP03 " << file << " scans=" << access.getNrSpectra() << " points=" << points
        << " exact_metadata=1 exact_arrays=1 extraction_query_counts=1,16,4096 cache_limits=2\n";
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
        // The vendor writer uses a standard CV array, not a non-standard label.
        floats.setName(i == 0 ? "raw inverse reduced ion mobility array" : "Ion Mobility");
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
    check(access.acquire(0)->array(2).description == "raw inverse reduced ion mobility array",
      "Indexed decoder lost the standard mobility name");
    // Use the same shared drift-filter overload as the transition-list filter.
    const auto reference = direct->lightClone()->getSpectrumById(0, 0.25, 0.25);
    const auto filtered = access.lightClone()->getSpectrumById(0, 0.25, 0.25);
    check(filtered->getDriftTimeArray() && filtered->getMZArray()->data.size() == 32,
      "Indexed drift filtering lost boundary points");
    check(filtered->getMZArray()->data == reference->getMZArray()->data &&
      filtered->getIntensityArray()->data == reference->getIntensityArray()->data &&
      filtered->getDriftTimeArray()->data == reference->getDriftTimeArray()->data,
      "Indexed drift filtering changed values");
    check(access.lightClone()->getSpectrumById(0, 0.5, 1)->getMZArray()->data.empty(),
      "Indexed drift filtering kept points outside the interval");
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

  void detachedCopies(const IndexedFile& indexed)
  {
    auto cache = std::make_shared<RawCache>(MIB, 0, MIB);
    RawAccess access(indexed.source(), cache);
    auto scan = access.acquire(0);
    const auto decoded = cache->snapshot();
    auto copy = access.getSpectrumById(0);
    const auto copied = cache->snapshot();
    check(copied.cache_entries == 1 && copied.evictions == decoded.evictions,
      "Fitting detached copy evicted its cached scan");
    check(copied.budget.committed_bytes > decoded.budget.committed_bytes && copied.budget.reserved_bytes == 0,
      "Detached copy did not retain its charge");
    auto array = copy->getMZArray();
    array->data[0] = -1;
    check(scan->array(0).data[0] == 400, "Detached copy mutation changed cached scan");
    copy.reset();
    cache->clear();
    scan.reset();
    check(cache->snapshot().budget.committed_bytes > 0, "Live detached array lost its charge");
    array.reset();
    check(cache->snapshot().budget.committed_bytes == 0, "Detached array leaked its charge");

    copy = access.getSpectrumById(0);
    cache->clear();
    copy->getDataArrays().clear();
    check(cache->snapshot().budget.committed_bytes > 0, "Live detached container lost its charge");
    copy.reset();
    check(cache->snapshot().budget.committed_bytes == 0, "Detached container leaked its charge");

    // Pin real copies until the actual remaining capacity cannot admit another.
    scan = access.acquire(0);
    std::vector<OpenSwath::SpectrumPtr> pinned;
    bool rejected = false;
    for (size_t i = 0; i != 4096 && !rejected; ++i)
    {
      const auto before = cache->snapshot().budget.committed_bytes;
      try { pinned.push_back(access.getSpectrumById(0)); }
      catch (const OpenMS::Exception::BaseException&)
      {
        rejected = true;
        check(cache->snapshot().budget.committed_bytes == before, "Rejected copy leaked a charge");
      }
      const auto budget = cache->snapshot().budget;
      check(budget.reserved_bytes == 0 && budget.peak_total_bytes <= budget.usable_bytes,
        "Pinned detached copies exceeded the cap or leaked a reservation");
    }
    check(!pinned.empty() && rejected, "Pinned detached copies did not reach bounded admission failure");
    cache->clear();
    scan.reset();
    check(cache->snapshot().budget.committed_bytes > 0, "Pinned detached copies lost their charges");
    pinned.clear();
    check(cache->snapshot().budget.committed_bytes == 0, "Pinned detached copies leaked their charges");
  }

  void leases(const IndexedFile& indexed)
  {
    auto probe = std::make_shared<RawCache>(4 * MIB, 0, MIB);
    RawAccess measure(indexed.source(), probe);
    measure.acquire(0);
    const auto bytes = probe->snapshot().budget.committed_bytes;
    auto exact_cache = std::make_shared<RawCache>(bytes * 3, 0, bytes);
    RawAccess exact_access(indexed.source(), exact_cache);
    exact_access.acquire(0);
    rejects([&] { exact_access.getSpectrumById(0); });
    check(exact_cache->snapshot().budget.committed_bytes == bytes && exact_cache->snapshot().budget.reserved_bytes == 0,
      "Rejected mutable copy changed retained scan accounting");
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

  void readOnlyExtraction(const IndexedFile& indexed)
  {
    auto probe = std::make_shared<RawCache>(MIB, 0, MIB);
    RawAccess(indexed.source(), probe).acquire(0);
    const auto bytes = probe->snapshot().budget.committed_bytes;
    auto cache = std::make_shared<RawCache>(2 * bytes, 0, bytes);
    auto access = std::make_shared<RawAccess>(indexed.source(), cache);
    using Extractor = OpenMS::ChromatogramExtractorAlgorithm;
    std::vector<Extractor::ExtractionCoordinates> coordinates(1);
    coordinates[0].mz = 401;
    coordinates[0].ion_mobility = 0.25;
    coordinates[0].rt_end = -1;
    const auto extract = [&](const OpenSwath::SpectrumAccessPtr& source, double im_window)
    {
      std::vector<OpenSwath::ChromatogramPtr> output{std::make_shared<OpenSwath::Chromatogram>()};
      Extractor().extractChromatograms(source, output, coordinates, 0.1, false, im_window, "tophat");
      return output.front();
    };
    for (const double im_window : {-1.0, 0.1})
    {
      const auto expected = extract(indexed.source(), im_window);
      const auto actual = extract(access->lightClone(), im_window);
      check(actual->getTimeArray()->data == expected->getTimeArray()->data &&
        actual->getIntensityArray()->data == expected->getIntensityArray()->data,
        "Read-only extraction changed IM or non-IM evidence");
    }
    rejects([&] { access->getSpectrumById(0); });
    cache->clear();
    check(cache->snapshot().budget.committed_bytes == 0 && cache->snapshot().budget.reserved_bytes == 0,
      "Read-only extraction retained a scan or reservation");

    for (const int column : {0, 1, 2})
    {
      auto retained_cache = std::make_shared<RawCache>(2 * bytes, 0, bytes);
      auto reader = std::make_shared<RawAccess>(indexed.source(), retained_cache);
      auto scan = reader->acquire(0);
      auto view = reader->lightClone()->getSpectrumView(0);
      check(view.mz.get() == &scan->array(0) && view.intensity.get() == &scan->array(1) &&
        view.ion_mobility.get() == &scan->array(2), "Read-only view copied cached arrays");
      check(retained_cache->snapshot().budget.committed_bytes == bytes,
        "Read-only view added a numeric charge");
      auto array = column == 0 ? view.mz : column == 1 ? view.intensity : view.ion_mobility;
      view = {};
      scan.reset();
      retained_cache->cancel();
      rejects([&] { reader->getSpectrumView(0); });
      reader.reset();
      check(array->data.size() == 32 && retained_cache->snapshot().budget.committed_bytes == bytes,
        "Cancellation or reader destruction revoked a live array lease");
      array.reset();
      check(retained_cache->snapshot().budget.committed_bytes == 0,
        "Last read-only array owner leaked its scan charge");
    }

    // A transforming reader must use its existing mutable path, never the uncalibrated raw view.
    auto roomy = std::make_shared<RawCache>(4 * MIB, 0, MIB);
    auto tiled = std::make_shared<RawAccess>(indexed.source(), roomy);
    for (const bool ppm : {false, true})
    {
      auto direct = std::make_shared<OpenMS::SpectrumAccessQuadMZTransforming>(indexed.source(), 0.01, 1, 1e-8, ppm);
      auto calibrated = std::make_shared<OpenMS::SpectrumAccessQuadMZTransforming>(tiled->lightClone(), 0.01, 1, 1e-8, ppm);
      const auto expected = direct->getSpectrumById(0);
      const auto actual = calibrated->lightClone()->getSpectrumView(0);
      check(actual.mz->data == expected->getMZArray()->data &&
        actual.mz->data != tiled->acquire(0)->array(0).data,
        "Spectrum view bypassed quadratic or ppm calibration");
      const auto expected_chrom = extract(direct, 0.1);
      const auto actual_chrom = extract(calibrated->lightClone(), 0.1);
      check(actual_chrom->getTimeArray()->data == expected_chrom->getTimeArray()->data &&
        actual_chrom->getIntensityArray()->data == expected_chrom->getIntensityArray()->data,
        "Calibrated view extraction changed evidence");
    }
    roomy->clear();
    check(roomy->snapshot().budget.committed_bytes == 0, "Calibration fallback retained a view charge");
    OPENMS_LOG_INFO << "WP05 exact IM/non-IM extraction without mutable admission, zero-copy aliases, "
      "all array leases/cancellation and quadratic/ppm fallback passed\n";
  }

  void phaseCounters()
  {
    using Scope = OpenSwath::SpectrumAccessScope;
    using Purpose = Scope::Purpose;
    check(Scope::current() == Purpose::UNSPECIFIED, "Read purpose leaked into caller");
    auto experiment = smallRun();
    experiment[0].setMSLevel(1);
    IndexedFile indexed(experiment);
    auto cache = std::make_shared<RawCache>(4 * MIB, 0, MIB);
    auto access = std::make_shared<RawAccess>(indexed.source(), cache);
    access->acquire(0);
    {
      Scope scoring(Purpose::SCORING);
      auto child = std::async(std::launch::async, []
      {
        check(Scope::current() == Purpose::UNSPECIFIED, "Purpose leaked across threads");
        Scope extraction(Purpose::EXTRACTION);
        check(Scope::current() == Purpose::EXTRACTION, "Worker purpose was not set");
      });
      child.get();
      try { Scope extraction(Purpose::EXTRACTION); throw std::runtime_error("scope check"); }
      catch (const std::runtime_error&) {}
      check(Scope::current() == Purpose::SCORING, "Nested exception did not restore purpose");
      cache->setCalibration(true);
      access->getSpectrumById(0);
      cache->setCalibration(false);
      using Extractor = OpenMS::ChromatogramExtractorAlgorithm;
      std::vector<Extractor::ExtractionCoordinates> coordinates(1);
      coordinates[0].mz = 401;
      coordinates[0].rt_end = -1;
      const auto extract = [&](const OpenSwath::SpectrumAccessPtr& reader)
      {
        std::vector<OpenSwath::ChromatogramPtr> output{std::make_shared<OpenSwath::Chromatogram>()};
        Extractor().extractChromatograms(reader, output, coordinates, 0.1, false, -1, "tophat");
      };
      extract(access);
      extract(std::make_shared<OpenMS::SpectrumAccessQuadMZTransforming>(access, 0, 1, 0, false));
      check(Scope::current() == Purpose::SCORING, "Actual extractor did not restore caller purpose");
      access->getSpectrumById(1);
    }
    check(Scope::current() == Purpose::UNSPECIFIED, "Purpose leaked after scope exit");
    const auto stats = cache->snapshot();
    const auto& phase = stats.phases;
    check(stats.hits == 7 && stats.misses == 4 && stats.work.decodes == 4 && stats.work.copies == 6,
      "Read/decode/copy totals changed");
    check(phase[static_cast<size_t>(RawPhase::UNCLASSIFIED)].misses == 1 &&
      phase[static_cast<size_t>(RawPhase::CALIBRATION)].copies == 1 &&
      phase[static_cast<size_t>(RawPhase::MS1_EXTRACTION)].hits == 2 &&
      phase[static_cast<size_t>(RawPhase::MS1_EXTRACTION)].copies == 1 &&
      phase[static_cast<size_t>(RawPhase::MS2_EXTRACTION)].misses == 3 &&
      phase[static_cast<size_t>(RawPhase::MS2_EXTRACTION)].copies == 3 &&
      phase[static_cast<size_t>(RawPhase::SCORING)].copies == 1,
      "Calibration, actual extraction or calibrated-wrapper reads were misclassified");
    check(stats.work.decoded_numeric_bytes == 3072 && stats.work.copied_numeric_bytes == 5120 &&
      stats.work.read_ns >= stats.work.decode_ns && stats.work.decode_ns > 0 && stats.work.copy_ns > 0,
      "Numeric payload byte or elapsed counters changed");
    for (auto member : {&RawWorkCounters::hits, &RawWorkCounters::misses, &RawWorkCounters::evictions,
      &RawWorkCounters::decodes, &RawWorkCounters::copies, &RawWorkCounters::decoded_numeric_bytes,
      &RawWorkCounters::copied_numeric_bytes, &RawWorkCounters::read_ns, &RawWorkCounters::decode_ns, &RawWorkCounters::copy_ns})
    {
      ByteCount sum = 0;
      for (const auto& row : phase) sum += row.*member;
      check(sum == stats.work.*member, "Phase counters do not reconcile with independent totals");
    }
    auto pressure = std::make_shared<RawCache>(8192, 0, 8192);
    RawAccess reader(indexed.source(), pressure);
    reader.acquire(0);
    { Scope scoring(Purpose::SCORING); reader.acquire(1); }
    check(pressure->snapshot().evictions == 1 &&
      pressure->snapshot().phases[static_cast<size_t>(RawPhase::SCORING)].evictions == 1,
      "Eviction was not attributed to the consumer causing pressure");
    cache->clear();
    pressure->clear();
    check(cache->snapshot().budget.committed_bytes == 0 && pressure->snapshot().budget.committed_bytes == 0,
      "Phase counters retained scan ownership");
    OPENMS_LOG_INFO << "PASS phase counters: nested/thread/exception scopes, actual extractor, calibrated fallback, "
      "all ten totals, exact bytes, eviction owner and final release\n";
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

  // Optional real-input acceptance path. No conversion and no vendor fixture download.
  void nativeInput(const std::string& path)
  {
    std::shared_ptr<OpenMS::ExperimentalSettings> settings;
    auto maps = OpenMS::SwathFile().loadOnDisc(path, settings);
    check(!maps.empty() && settings, "Native input has no spectrum maps or settings");
    std::vector<OpenSwath::SwathMap> streamed;
#ifdef WITH_OPENTIMS
    if (OpenMS::FileHandler::getTypeByFileName(path) == OpenMS::FileTypes::BRUKER_TDF)
    {
      streamed = OpenMS::SwathFile().loadBrukerTdf(path, "", settings, "normal");
      check(streamed.size() == maps.size(), "Native and streaming map membership differs");
    }
#endif
    auto cache = std::make_shared<RawCache>(256 * MIB, 64 * MIB, 32 * MIB);
    OpenMS::Size checked = 0;
    for (const auto& map : maps)
    {
      auto source = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr);
      check(source && source->getNrSpectra(), "Native map is empty or not random access");
      const auto reference = std::find_if(streamed.begin(), streamed.end(), [&](const auto& other)
      { return map.isEqual(other); });
      check(streamed.empty() || reference != streamed.end(), "Native isolation/IM window differs from streaming");
      if (reference != streamed.end())
      {
        check(reference->sptr->getNrSpectra() == source->getNrSpectra(), "Native scan membership differs from streaming");
      }
      auto access = std::make_shared<RawAccess>(source, cache);
      std::vector<OpenMS::Size> indices;
      const int count = static_cast<int>(source->getNrSpectra());
      for (const int i : {0, count / 2, count - 1})
      {
        auto experiment = source->getExperiment().clone();
        const auto index = source->getExperimentIndex(i);
        const auto expected = reference == streamed.end() ? source->getSpectrumById(i) : reference->sptr->getSpectrumById(i);
        if (reference != streamed.end())
        {
          const auto a = source->getSpectrumMetaById(i), b = reference->sptr->getSpectrumMetaById(i);
          check(a.RT == b.RT && a.ms_level == b.ms_level, "Native RT/MS level differs from streaming");
        }
        exact(*access->acquire(i), *expected, experiment->getSpectrumMetadata(index));
        const auto view = access->lightClone()->getSpectrumView(i);
        check(view.mz->data == expected->getMZArray()->data &&
          view.intensity->data == expected->getIntensityArray()->data,
          "Native cached view changed numeric data");
        check(source->getSpectraByRT(source->getSpectrumMetaById(i).RT, 0) ==
          access->getSpectraByRT(source->getSpectrumMetaById(i).RT, 0), "Native RT lookup changed");
        indices.push_back(index);
        ++checked;
      }
      auto sample = std::make_shared<OpenMS::SpectrumAccessOnDisc>(
        std::shared_ptr<OpenMS::IMSExperimentRandomAccess>(source->getExperiment().clone()), indices);
      const auto first = sample->getSpectrumById(0);
      if (first->getMZArray()->data.empty()) continue;
      using Extractor = OpenMS::ChromatogramExtractorAlgorithm;
      std::vector<Extractor::ExtractionCoordinates> coordinates(1);
      const auto peak = first->getMZArray()->data.size() / 2;
      coordinates[0].mz = first->getMZArray()->data[peak];
      coordinates[0].rt_end = -1;
      const auto im = first->getDriftTimeArray();
      coordinates[0].ion_mobility = im && !im->data.empty() ? im->data[peak] : 0;
      for (const double window : {-1.0, 0.05})
      {
        if (window > 0 && (!im || im->data.empty())) continue;
        const auto extract = [&](const OpenSwath::SpectrumAccessPtr& reader)
        {
          std::vector<OpenSwath::ChromatogramPtr> output{std::make_shared<OpenSwath::Chromatogram>()};
          Extractor().extractChromatograms(reader, output, coordinates, 0.1, false, window, "tophat");
          return output.front();
        };
        const auto expected = extract(sample);
        const auto actual = extract(std::make_shared<RawAccess>(sample, cache));
        check(expected->getTimeArray()->data == actual->getTimeArray()->data &&
          expected->getIntensityArray()->data == actual->getIntensityArray()->data,
          "Native cached extraction changed IM or non-IM evidence");
      }
    }
    cache->clear();
    const auto stats = cache->snapshot();
    check(stats.budget.committed_bytes == 0 && stats.budget.reserved_bytes == 0,
      "Native input retained a scan charge");
    OPENMS_LOG_INFO << "PASS native input " << path << " maps=" << maps.size()
      << " sampled_scans=" << checked << " streaming_maps=" << streamed.size()
      << " exact_metadata_arrays_views_RT_extraction=1 final_live_bytes=0\n";
  }

  void truncated()
  {
    IndexedFile indexed(smallRun());
    // Remove the native index. Do not let the adapter treat a bad file as an empty run.
    std::filesystem::resize_file(indexed.path, std::filesystem::file_size(indexed.path) / 2);
    bool rejected = false;
    try
    {
      auto experiment = std::make_shared<OpenMS::OnDiscMSExperiment>();
      rejected = !experiment->openFile(indexed.path);
      if (!rejected)
      {
        RawAccess access(std::make_shared<OpenMS::SpectrumAccessOnDisc>(experiment),
          std::make_shared<RawCache>(4 * MIB, MIB, MIB));
        access.acquire(0);
      }
    }
    catch (const OpenMS::Exception::BaseException&) { rejected = true; }
    check(rejected, "Truncated indexed input was accepted");
    IndexedFile payload(smallRun());
    RawAccess access(payload.source(), std::make_shared<RawCache>(4 * MIB, MIB, MIB));
    // Keep the valid native index, but remove the bytes it addresses.
    std::filesystem::resize_file(payload.path, 0);
    rejects([&] { access.acquire(0); });
  }
}

int main(int argc, char** argv)
{
  try
  {
    if (argc > 1)
    {
      for (int i = 1; i < argc; ++i) nativeInput(argv[i]);
      return 0;
    }
    fixtures();
    IndexedFile indexed(smallRun());
    detachedCopies(indexed);
    queries(indexed);
    leases(indexed);
    readOnlyExtraction(indexed);
    phaseCounters();
    concurrent(indexed);
    truncated();
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
}
