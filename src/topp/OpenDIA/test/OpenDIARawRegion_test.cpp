// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include "OpenDIARawRegion.h"

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOnDisc.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <stdexcept>

using namespace OpenDIA;
namespace
{
constexpr ByteCount MIB = 1024 * 1024;
const RecordContext context {0, 71, 9, 203};
void check(bool condition, const char* message)
{
  if (! condition) { throw std::runtime_error(message); }
}
template<typename F>
void rejects(F action)
{
  bool rejected = false;
  try
  {
    action();
  }
  catch (const OpenMS::Exception::BaseException&)
  {
    rejected = true;
  }
  check(rejected, "Expected explicit rejection");
}
struct IndexedFile
{
  std::string path = (std::filesystem::current_path() / ("OpenDIARawRegion_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
  explicit IndexedFile(const OpenMS::MSExperiment& experiment)
  {
    OpenMS::MzMLFile writer;
    writer.getOptions().setWriteIndex(true);
    writer.store(path, experiment);
  }
  ~IndexedFile()
  {
    std::error_code error;
    std::filesystem::remove(path, error);
  }
  std::unique_ptr<OpenMS::IMSExperimentRandomAccess> source() const
  {
    auto result = std::make_unique<OpenMS::OnDiscMSExperiment>();
    check(result->openFile(path), "Cannot open indexed fixture");
    return result;
  }
};

// Deliberately malformed native decoder results cannot be written as valid mzML.
struct MemorySource : OpenMS::IMSExperimentRandomAccess
{
  std::vector<OpenMS::MSSpectrum> spectra;
  std::stop_source* cancel_on_read = nullptr;
  OpenMS::Size getNrSpectra() const override
  {
    return spectra.size();
  }
  OpenMS::Size getNrChromatograms() const override
  {
    return 0;
  }
  OpenMS::MSSpectrum getSpectrum(OpenMS::Size i) override
  {
    if (cancel_on_read) { cancel_on_read->request_stop(); }
    return spectra.at(i);
  }
  OpenMS::MSSpectrum getSpectrumMetadata(OpenMS::Size i) override
  {
    return spectra.at(i);
  }
  double getSpectrumRT(OpenMS::Size i) const override
  {
    return spectra.at(i).getRT();
  }
  OpenMS::UInt getSpectrumMSLevel(OpenMS::Size i) const override
  {
    return spectra.at(i).getMSLevel();
  }
  bool isSortedByRT() const override
  {
    return true;
  }
  OpenMS::MSChromatogram getChromatogram(OpenMS::Size) override
  {
    return {};
  }
  std::shared_ptr<const OpenMS::ExperimentalSettings> getExperimentalSettings() const override
  {
    return {};
  }
  std::unique_ptr<OpenMS::IMSExperimentRandomAccess> clone() const override
  {
    return std::make_unique<MemorySource>(*this);
  }
};
MemorySource small()
{
  MemorySource source;
  OpenMS::MSSpectrum scan;
  scan.setNativeID("scan=91");
  scan.setRT(1.234567890123);
  scan.setMSLevel(2);
  for (const auto& [mz, intensity] : std::array<std::pair<double, float>, 4> {{{102, 8}, {100, 1}, {101, 2}, {101, 4}}})
  {
    OpenMS::Peak1D p;
    p.setMZ(mz);
    p.setIntensity(intensity);
    scan.push_back(p);
  }
  OpenMS::MSSpectrum::FloatDataArray im;
  im.setName("raw inverse reduced ion mobility array");
  im.assign({1, 0, 0, 1});
  scan.getFloatDataArrays().push_back(im);
  scan.setDriftTime(0);
  scan.setDriftTimeUnit(OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  source.spectra.push_back(scan);
  scan.getFloatDataArrays().clear();
  scan.setDriftTimeUnit(OpenMS::DriftTimeUnit::NONE);
  source.spectra.push_back(scan);
  scan.clear(true);
  scan.setNativeID("empty");
  scan.setMSLevel(1);
  source.spectra.push_back(scan);
  return source;
}
void calibratedCreation()
{
  auto source = small();
  // A sorted source can become nonmonotonic after calibration.
  source.spectra[0].sortByPosition();
  const auto original = source.spectra[0];
  const auto bytes = RawRegion::requiredBytes({1, original.size()}, original.getNativeID().size());
  ByteBudget budget(bytes, 0);
  for (const auto calibration : {MZCalibration {10201, -202, 1, false}, MZCalibration {4.25, 0.125, 0.0003, true}, MZCalibration {0, 1, 0, false}})
  {
    auto adapter = std::make_shared<OpenMS::SpectrumAccessOnDisc>(std::shared_ptr<OpenMS::IMSExperimentRandomAccess>(source.clone()),
                                                                  std::vector<OpenMS::Size> {0});
    OpenMS::SpectrumAccessQuadMZTransforming oracle(adapter, calibration.a, calibration.b, calibration.c, calibration.ppm);
    const auto expected = oracle.getSpectrumById(0);
    const auto& values = expected->getMZArray()->data;
    std::vector<OpenMS::Size> order(values.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return values[a] < values[b]; });
    {
      auto region = RawRegion::fromSpectrum(source.spectra[0], 91, context, budget, calibration);
      auto scan = region.scan(0);
      check(region.size() == 1 && scan.acquisition().source_index == 91 && scan.context() == context && scan.nativeID() == original.getNativeID()
              && scan.acquisition().rt_seconds == original.getRT() && scan.acquisition().ms_level == original.getMSLevel()
              && scan.acquisition().faims_voltage == 0,
            "Decoded creation changed identity or metadata");
      for (OpenMS::Size p = 0; p < order.size(); ++p)
      {
        check(scan.mz()[p] == values[order[p]] && scan.intensity()[p] == original[order[p]].getIntensity()
                && scan.mobility()[p] == original.getFloatDataArrays()[0][order[p]],
              "Calibration differs from actual transform or lost stable point alignment");
      }
      check(budget.snapshot().committed_bytes == bytes, "Decoded creation is not exact-size admitted");
      rejects([&] { RawRegion::fromSpectrum(original, 0, context, budget); });
    }
    check(source.spectra[0] == original, "Calibration mutated native source");
    check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Calibration leaked final lease");
  }
  {
    auto region = RawRegion::fromSpectrum(original, 0, context, budget);
    check(region.scan(0).mz()[0] == original[0].getMZ(), "Absent calibration changed native values");
  }
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  const double max = std::numeric_limits<double>::max();
  for (const auto calibration : {MZCalibration {nan, 1, 0}, MZCalibration {0, inf, 0}, MZCalibration {0, 1, nan}, MZCalibration {-1, 0, 0},
                                 MZCalibration {0, max, max}, MZCalibration {max, 0, 0, true}})
  {
    rejects([&] { RawRegion::fromSpectrum(original, 0, context, budget, calibration); });
    check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Failed calibration leaked admission");
  }
  auto invalid_context = context;
  invalid_context.schema_version = 1;
  rejects([&] { RawRegion::fromSpectrum(original, 0, invalid_context, budget); });
  ByteBudget tiny(bytes - 1, 0);
  rejects([&] { RawRegion::fromSpectrum(original, 0, context, tiny); });
  {
    auto empty = RawRegion::fromSpectrum(source.spectra[2], 2, context, budget);
    check(empty.size() == 1 && empty.scan(0).mz().empty(), "Decoded empty row lost");
  }
  auto scalar = original;
  scalar.setDriftTimeUnit(OpenMS::DriftTimeUnit::VSSC);
  scalar.setDriftTime(0.123456789012345);
  {
    auto region = RawRegion::fromSpectrum(scalar, 7, context, budget);
    const auto mobility = region.scan(0).acquisition().acquisition_mobility;
    check(mobility && mobility->value == scalar.getDriftTime() && mobility->unit == MobilityUnit::INVERSE_REDUCED
            && ! region.scan(0).acquisition().faims_voltage,
          "Decoded scalar mobility changed");
  }
  check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Decoded creation leaked final lease");
}
void contracts()
{
  auto source = small();
  ByteBudget budget(MIB, 4096);
  const std::array<OpenMS::Size, 4> indices {1, 0, 2, 0};
  const RecordLimits limits {4, 12};
  const auto bytes = RawRegion::requiredBytes(limits, 16);
  RawRegion::Scan retained;
  {
    auto region = RawRegion::read(source, indices, context, limits, 16, budget);
    auto copy = region;
    retained = copy.scan(1).range({101, 101, EdgePolicy::CLOSED});
    check(region.size() == 4 && region.scan(0).acquisition().source_index == 1 && region.scan(3).acquisition().source_index == 0,
          "Source indices renumbered/reordered");
    check(region.scan(1).acquisition().rt_seconds == source.spectra[0].getRT() && retained.context() == context && retained.nativeID() == "scan=91",
          "Identity changed");
    check(region.scan(0).mobility().empty() && ! region.scan(0).acquisition().faims_voltage && region.scan(1).acquisition().faims_voltage == 0
            && region.scan(1).mobility()[0] == 0,
          "Absent IM/FAIMS became measured zero");
    check(region.scan(2).mz().empty(), "Empty acquisition lost");
    check(retained.mz().size() == 2 && retained.intensity()[0] == 2 && retained.intensity()[1] == 4 && retained.mobility()[0] == 0
            && retained.mobility()[1] == 1,
          "Stable duplicates or IM association changed");
    const std::array<EdgePolicy, 4> edges {EdgePolicy::OPEN, EdgePolicy::CLOSED, EdgePolicy::LOWER_CLOSED, EdgePolicy::UPPER_CLOSED};
    const std::array<double, 4> expected {6, 15, 7, 14};
    for (OpenMS::Size i = 0; i < edges.size(); ++i)
    {
      std::array<RegionSum, 3> requests {{{1, {100, 102, edges[i]}, {}},
                                          {1, {101, 101, edges[i]}, {}},
                                          {1, {99, 103}, MobilityInterval {{0, 0, EdgePolicy::CLOSED}, MobilityUnit::INVERSE_REDUCED}}}};
      std::array<double, 3> out {};
      region.sum(requests, out);
      check(out[0] == expected[i] && out[1] == (edges[i] == EdgePolicy::CLOSED ? 6 : 0) && out[2] == 3,
            "Literal m/z or physical mobility edges wrong");
    }
    const std::array<double, 4> mobility_expected {0, 15, 3, 12};
    for (OpenMS::Size i = 0; i < edges.size(); ++i)
    {
      std::array<RegionSum, 1> query {{{1, {99, 103}, MobilityInterval {{0, 1, edges[i]}, MobilityUnit::INVERSE_REDUCED}}}};
      std::array<double, 1> output {};
      region.sum(query, output);
      check(output[0] == mobility_expected[i], "Literal physical mobility boundaries wrong");
    }
    std::array<RegionSum, 2> requests {
      {{1, {102, 103, EdgePolicy::CLOSED}, {}}, {0, {99, 103}, MobilityInterval {{0, 1}, MobilityUnit::INVERSE_REDUCED}}}};
    std::array<double, 2> out {-7, -7};
    rejects([&] { region.sum(requests, out); });
    check(out[0] == -7 && out[1] == -7, "Absent mobility batch partially wrote output");
    requests[1].mobility.reset();
    requests[1].mz = {99, 101, EdgePolicy::CLOSED};
    region.sum(requests, out);
    check(out[0] == 8 && out[1] == 7, "Unsorted queries changed result");
    requests[1].acquisition = 1;
    requests[1].mobility = MobilityInterval {{0, 1}, MobilityUnit::MILLISECOND};
    out = {-7, -7};
    rejects([&] { region.sum(requests, out); });
    check(out[0] == -7 && out[1] == -7, "Invalid batch partially wrote output");
    requests[1].mobility.reset();
    requests[1].mz = {2, 1};
    rejects([&] { region.sum(requests, out); });
    requests[1].mz = {0, 1, static_cast<EdgePolicy>(99)};
    rejects([&] { region.sum(requests, out); });
    requests[1].mz = {0, std::numeric_limits<double>::infinity()};
    rejects([&] { region.sum(requests, out); });
    rejects([&] { region.sum(requests, std::span<double>(out).first(1)); });
    rejects([&] { region.scan(4); });
    auto moved = std::move(copy);
    rejects([&] { copy.size(); });
    check(moved.size() == 4 && budget.snapshot().committed_bytes == bytes, "Copies changed accounting");
  }
  check(retained.intensity()[1] == 4 && budget.snapshot().committed_bytes == bytes, "View lost immutable ownership");
  auto retained_copy = retained;
  retained = {};
  check(retained_copy.intensity()[0] == 2 && budget.snapshot().committed_bytes == bytes, "Copied scan lost lease");
  retained_copy = {};
  check(budget.snapshot().committed_bytes == 0, "Last view leaked charge");
  rejects([&] { retained.nativeID(); });
  auto read = [&](RecordLimits bounds, OpenMS::Size id = 16, std::stop_token stop = {}) {
    return RawRegion::read(source, indices, context, bounds, id, budget, stop);
  };
  rejects([&] { read({3, 12}); });
  rejects([&] { read({4, 11}); });
  rejects([&] { read(limits, 1); });
  rejects([&] { RawRegion::requiredBytes({std::numeric_limits<OpenMS::Size>::max(), 1}, 16); });
  rejects([&] { RawRegion::requiredBytes({1, std::numeric_limits<OpenMS::Size>::max()}, 16); });
  std::stop_source stop;
  stop.request_stop();
  rejects([&] { read(limits, 16, stop.get_token()); });
  std::stop_source during;
  source.cancel_on_read = &during;
  rejects([&] { read(limits, 16, during.get_token()); });
  source.cancel_on_read = nullptr;
  source.spectra[0].getFloatDataArrays()[0].pop_back();
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0][0].setMZ(std::numeric_limits<double>::quiet_NaN());
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0].getFloatDataArrays()[0][0] = std::numeric_limits<float>::infinity();
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0][0].setIntensity(std::numeric_limits<float>::infinity());
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0].setRT(std::numeric_limits<double>::quiet_NaN());
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0].setDriftTime(std::numeric_limits<double>::infinity());
  rejects([&] { read(limits); });
  source = small();
  source.spectra[0].setDriftTimeUnit(OpenMS::DriftTimeUnit::MILLISECOND);
  {
    auto region = read(limits);
    const auto scalar = region.scan(1).acquisition().acquisition_mobility;
    check(scalar && scalar->value == 0 && scalar->unit == MobilityUnit::MILLISECOND
            && region.scan(1).acquisition().point_mobility == MobilityUnit::INVERSE_REDUCED && ! region.scan(1).acquisition().faims_voltage,
          "Scalar mobility collapsed into point IM/FAIMS");
  }
  source.spectra[0].setDriftTimeUnit(OpenMS::DriftTimeUnit::VSSC);
  source.spectra[0].setDriftTime(0.123456789012345);
  {
    auto region = read(limits);
    const auto scalar = region.scan(1).acquisition().acquisition_mobility;
    check(scalar && scalar->value == source.spectra[0].getDriftTime() && scalar->unit == MobilityUnit::INVERSE_REDUCED,
          "Scalar mobility precision changed");
  }
  source.spectra[0].setDriftTime(std::numeric_limits<double>::quiet_NaN());
  rejects([&] { read(limits); });
  source.spectra[0].setDriftTime(0);
  source.spectra[0].setDriftTimeUnit(OpenMS::DriftTimeUnit::CCS);
  rejects([&] { read(limits); });
  source = small();
  ByteBudget tiny(bytes - 1, 0);
  rejects([&] { RawRegion::read(source, indices, context, limits, 16, tiny); });
  const std::array<OpenMS::Size, 1> invalid_index {3};
  rejects([&] { RawRegion::read(source, invalid_index, context, limits, 16, budget); });
  // Fixed allocation, not observed payload size, remains charged while pinned.
  ByteBudget exact(bytes, 0);
  auto pinned = RawRegion::read(source, indices, context, limits, 16, exact).scan(0);
  rejects([&] { RawRegion::read(source, indices, context, limits, 16, exact); });
  exact.cancel();
  check(pinned.intensity()[0] == 1 && exact.snapshot().committed_bytes == bytes, "Cancellation revoked a pinned view");
  pinned = {};
  check(exact.snapshot().committed_bytes == 0, "Cancelled pinned view leaked charge");
  {
    auto empty = RawRegion::read(source, {}, context, {0, 0}, 0, budget);
    check(empty.size() == 0, "Empty region is not empty");
    empty.sum({}, {});
  }
  {
    auto invalid_context = context;
    invalid_context.schema_version = 1;
    rejects([&] { RawRegion::read(source, indices, invalid_context, limits, 16, budget); });
  }
  RawRegion::Scan detached;
  {
    ByteBudget local(bytes, 0);
    detached = RawRegion::read(source, indices, context, limits, 16, local).scan(1);
  }
  check(detached.mz().size() == 4 && detached.nativeID() == "scan=91", "Budget destruction revoked storage");
  detached = {};
  // Test the other physical point unit without coercing FAIMS into mobility.
  source.spectra[0].getFloatDataArrays()[0].setName("Ion Mobility");
  auto milliseconds = RawRegion::read(source, indices, context, limits, 16, budget);
  check(milliseconds.scan(1).acquisition().point_mobility == MobilityUnit::MILLISECOND, "Millisecond point unit changed");
  std::array<RegionSum, 1> ms_query {{{1, {99, 103}, MobilityInterval {{0, 0, EdgePolicy::CLOSED}, MobilityUnit::MILLISECOND}}}};
  std::array<double, 1> ms_output {};
  milliseconds.sum(ms_query, ms_output);
  check(ms_output[0] == 3, "Millisecond zero mobility lost");
  ms_query[0].mobility->unit = static_cast<MobilityUnit>(99);
  rejects([&] { milliseconds.sum(ms_query, ms_output); });
  milliseconds = {};
  check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Failed producer leaked bytes");
  budget.cancel();
  rejects([&] { read(limits); });
}

void compare(OpenMS::IMSExperimentRandomAccess& source, const std::vector<OpenMS::Size>& indices)
{
  RecordLimits limits {indices.size(), 0};
  OpenMS::Size id_capacity = 0;
  for (auto i : indices)
  {
    auto scan = source.getSpectrum(i);
    limits.values = checkedByteAdd(limits.values, scan.size());
    id_capacity = std::max(id_capacity, scan.getNativeID().size());
  }
  ByteBudget budget(512 * MIB, 64 * MIB);
  {
    auto region = RawRegion::read(source, indices, context, limits, id_capacity, budget);
    std::array<bool, 2> extracted {false, false};
    for (OpenMS::Size r = 0; r < indices.size(); ++r)
    {
      auto direct = source.getSpectrum(indices[r]);
      auto view = region.scan(r);
      check(view.acquisition().source_index == indices[r] && view.acquisition().rt_seconds == direct.getRT()
              && view.acquisition().ms_level == direct.getMSLevel() && view.nativeID() == direct.getNativeID(),
            "Native acquisition identity changed");
      const std::optional<double> faims = direct.getDriftTimeUnit() == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE
                                            ? std::optional<double> {direct.getDriftTime()}
                                            : std::nullopt;
      check(view.acquisition().faims_voltage == faims && view.acquisition().point_mobility.has_value() == direct.containsIMData(),
            "Native FAIMS/point mobility presence changed");
      if (direct.getDriftTimeUnit() == OpenMS::DriftTimeUnit::MILLISECOND || direct.getDriftTimeUnit() == OpenMS::DriftTimeUnit::VSSC)
      {
        const auto scalar = view.acquisition().acquisition_mobility;
        const auto expected_unit
          = direct.getDriftTimeUnit() == OpenMS::DriftTimeUnit::MILLISECOND ? MobilityUnit::MILLISECOND : MobilityUnit::INVERSE_REDUCED;
        check(scalar && scalar->value == direct.getDriftTime() && scalar->unit == expected_unit, "Native scalar mobility changed");
      }
      else { check(! view.acquisition().acquisition_mobility, "Absent scalar mobility gained a value"); }
      std::vector<OpenMS::Size> order(direct.size());
      std::iota(order.begin(), order.end(), 0);
      std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return direct[a].getMZ() < direct[b].getMZ(); });
      const auto* im = direct.containsIMData() ? &direct.getFloatDataArrays()[direct.getIMData().first] : nullptr;
      for (OpenMS::Size p = 0; p < order.size(); ++p)
      {
        check(view.mz()[p] == direct[order[p]].getMZ() && view.intensity()[p] == direct[order[p]].getIntensity()
                && (! im || view.mobility()[p] == (*im)[order[p]]),
              "Native point precision or association changed");
      }
      if (view.mz().empty()) { continue; }
      std::vector<RegionSum> requests;
      std::vector<double> expected;
      for (auto p : {view.mz().size() * 3 / 4, view.mz().size() / 4, view.mz().size() / 2})
      {
        const double center = view.mz()[p];
        for (const bool with_im : {false, true})
        {
          if (with_im && ! im) { continue; }
          RegionSum query {r, {center - 0.025, center + 0.025}, {}};
          if (with_im)
          {
            query.mobility = MobilityInterval {{view.mobility()[p] - 0.025, view.mobility()[p] + 0.025}, *view.acquisition().point_mobility};
          }
          double sum = 0;
          for (OpenMS::Size j = 0; j < direct.size(); ++j)
          {
            if (direct[j].getMZ() > query.mz.lower && direct[j].getMZ() < query.mz.upper
                && (! with_im || ((*im)[j] > query.mobility->interval.lower && (*im)[j] < query.mobility->interval.upper)))
            {
              sum += direct[j].getIntensity();
            }
          }
          requests.push_back(query);
          expected.push_back(sum);
          // Real extractor oracle on one explicitly mapped native acquisition.
          // Only ordinary interior windows are compared, not legacy edge quirks.
          if (! extracted[with_im])
          {
            using Extractor = OpenMS::ChromatogramExtractorAlgorithm;
            auto adapter = std::make_shared<OpenMS::SpectrumAccessOnDisc>(std::shared_ptr<OpenMS::IMSExperimentRandomAccess>(source.clone()),
                                                                          std::vector<OpenMS::Size> {indices[r]});
            std::vector<Extractor::ExtractionCoordinates> coordinates(1);
            coordinates[0].mz = center;
            coordinates[0].ion_mobility = with_im ? view.mobility()[p] : 0;
            coordinates[0].rt_end = -1;
            std::vector<OpenSwath::ChromatogramPtr> output {std::make_shared<OpenSwath::Chromatogram>()};
            Extractor().extractChromatograms(adapter, output, coordinates, 0.05, false, with_im ? 0.05 : -1, "tophat");
            check(output[0]->getIntensityArray()->data.size() == 1 && output[0]->getTimeArray()->data[0] == direct.getRT()
                    && std::abs(output[0]->getIntensityArray()->data[0] - sum) <= 1e-6 * std::max(1.0, std::abs(sum)),
                  "Actual extractor disagrees with native column sum");
            extracted[with_im] = true;
          }
        }
      }
      std::vector<double> actual(requests.size());
      region.sum(requests, actual);
      check(actual == expected, "Batched native sums disagree with direct decoder");
      const auto start = std::chrono::steady_clock::now();
      for (int repeat = 0; repeat < 32; ++repeat)
      {
        region.sum(requests, actual);
      }
      if (r == 0)
      {
        OPENMS_LOG_INFO << "native-region rows_capacity=" << limits.records << " points_capacity=" << limits.values << " id_capacity=" << id_capacity
                        << " retained_bytes=" << budget.snapshot().committed_bytes
                        << " query_ns=" << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start).count()
                        << '\n';
      }
    }
  }
  check(budget.snapshot().committed_bytes == 0 && budget.snapshot().reserved_bytes == 0, "Reader comparison leaked region bytes");
}
void fixtures()
{
  for (const auto* file :
       {"topp/OpenSwathWorkflow_1_input.mzML", "topp/OpenSwathWorkflow_23_input.mzML", "class_tests/openms/data/IM_FAIMS_test.mzML"})
  {
    OpenMS::MSExperiment loaded;
    OpenMS::MzMLFile().load(std::string(OPENDIA_TEST_DATA_DIR) + "/" + file, loaded);
    IndexedFile indexed(loaded);
    auto source = indexed.source();
    std::vector<OpenMS::Size> indices(source->getNrSpectra());
    std::iota(indices.begin(), indices.end(), 0);
    compare(*source, indices);
  }
  auto synthetic = small();
  OpenMS::MSExperiment experiment;
  for (auto scan : synthetic.spectra)
  {
    scan.sortByPosition();
    experiment.addSpectrum(scan);
  }
  IndexedFile indexed(experiment);
  auto source = indexed.source();
  compare(*source, {1, 0, 2, 0});
}
void nativeInput(const std::string& path)
{
  std::shared_ptr<OpenMS::ExperimentalSettings> settings;
  auto maps = OpenMS::SwathFile().loadOnDisc(path, settings);
  check(! maps.empty(), "No native maps");
  for (const auto& map : maps)
  {
    auto adapter = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr);
    check(adapter && adapter->getNrSpectra(), "No native random-access acquisitions");
    auto source = adapter->getExperiment().clone();
    const auto count = adapter->getNrSpectra();
    std::vector<OpenMS::Size> indices;
    for (auto i : {OpenMS::Size {0}, count / 2, count - 1})
    {
      indices.push_back(adapter->getExperimentIndex(static_cast<int>(i)));
    }
    compare(*source, indices);
  }
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc > 1)
    {
      for (int i = 1; i < argc; ++i)
      {
        nativeInput(argv[i]);
      }
    }
    else
    {
      contracts();
      calibratedCreation();
      fixtures();
    }
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
}
