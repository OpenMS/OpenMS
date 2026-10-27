// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include "OpenDIARawAccess.h"
#include "OpenDIARegionWorkflow.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <array>
#include <filesystem>
#include <limits>
#include <stdexcept>

using namespace OpenDIA;
namespace
{
constexpr ByteCount MIB = 1024 * 1024;
using Coordinate = OpenMS::ChromatogramExtractor::ExtractionCoordinates;
using Chromatograms = std::vector<OpenSwath::ChromatogramPtr>;
void check(bool condition, const char* message)
{
  if (! condition) { throw std::runtime_error(message); }
}
template<typename F>
void rejects(F operation)
{
  bool rejected = false;
  try
  {
    operation();
  }
  catch (const OpenMS::Exception::BaseException&)
  {
    rejected = true;
  }
  check(rejected, "Expected explicit rejection");
}
Chromatograms empty(OpenMS::Size count)
{
  Chromatograms result;
  for (OpenMS::Size i = 0; i < count; ++i)
  {
    result.push_back(std::make_shared<OpenSwath::Chromatogram>());
  }
  return result;
}
OpenSwath::SpectrumAccessPtr transform(OpenSwath::SpectrumAccessPtr input, const std::optional<MZCalibration>& calibration)
{
  if (! calibration) { return input; }
  const auto& c = *calibration;
  return std::make_shared<OpenMS::SpectrumAccessQuadMZTransforming>(input, c.a, c.b, c.c, c.ppm);
}
void equal(const Chromatograms& actual, const Chromatograms& expected)
{
  check(actual.size() == expected.size(), "Chromatogram count changed");
  for (OpenMS::Size i = 0; i < actual.size(); ++i)
  {
    check(actual[i]->getTimeArray()->data == expected[i]->getTimeArray()->data, "RT arrays changed");
    check(actual[i]->getIntensityArray()->data == expected[i]->getIntensityArray()->data, "Intensity arrays changed");
  }
}
void exercise(const std::shared_ptr<OpenMS::SpectrumAccessOnDisc>& source, bool native = false)
{
  check(source->getNrSpectra() > 0, "Fixture has no acquisitions");
  const auto uncalibrated = source->getSpectrumById(0);
  if (uncalibrated->getMZArray()->data.empty()) { return; }
  const std::array<std::optional<MZCalibration>, 3> calibrations {
    {{}, MZCalibration {0.001, 1.00001, 1e-10, false}, MZCalibration {0.1, 0.00001, 1e-10, true}}};
  RegionWorkflow workflow(false, false, false, false, false, 1);
  bool signal_observed = false;
  for (const auto& calibration : calibrations)
  {
    if (native && ! calibration) { continue; }
    auto cache = std::make_shared<RawCache>(256 * MIB, 64 * MIB, 64 * MIB);
    auto access = std::make_shared<RawAccess>(source, cache);
    auto input = transform(access, calibration);
    auto reference = transform(source->lightClone(), calibration);
    const auto sample = reference->getSpectrumById(0);
    const auto mz = sample->getMZArray();
    const auto im = sample->getDriftTimeArray();
    std::vector<Coordinate> coordinates(16);
    for (OpenMS::Size q = 0; q < coordinates.size(); ++q)
    {
      const auto p = mz->data.size() * q / coordinates.size();
      coordinates[q].mz = mz->data[p];
      coordinates[q].ion_mobility = im ? im->data[p] : -1;
      coordinates[q].rt_start = source->getSpectrumMetaById(0).RT;
      coordinates[q].rt_end = q % 2 ? coordinates[q].rt_start + 1 : -1;
    }
    std::stable_sort(coordinates.begin(), coordinates.end(), [](const auto& a, const auto& b) { return a.mz < b.mz; });
    for (const bool ppm : {false, true})
    {
      for (const bool use_im : {false, true})
      {
        if (use_im && ! im) { continue; }
        OpenMS::ChromExtractParams cp {};
        cp.mz_extraction_window = ppm ? 20 : 0.05;
        cp.ppm = ppm;
        cp.im_extraction_window = use_im ? 0.05 : -1;
        cp.extraction_function = "tophat";
        auto expected = empty(coordinates.size()), actual = empty(coordinates.size());
        OpenMS::ChromatogramExtractor().extractChromatograms(reference, expected, coordinates, cp.mz_extraction_window, cp.ppm,
                                                             cp.im_extraction_window, cp.extraction_function);
        workflow.extractChromatograms_(input, actual, coordinates, cp);
        equal(actual, expected);
        for (const auto& chromatogram : expected)
        {
          signal_observed = signal_observed
                            || std::any_of(chromatogram->getIntensityArray()->data.begin(), chromatogram->getIntensityArray()->data.end(),
                                           [](double value) { return value != 0; });
        }
        const auto first = cache->snapshot();
        auto repeated = empty(coordinates.size());
        workflow.extractChromatograms_(input->lightClone(), repeated, coordinates, cp);
        equal(repeated, expected);
        const auto second = cache->snapshot();
        check(second.work.copies == 0 && second.work.copied_numeric_bytes == 0, "Region extraction made a detached legacy copy");
        // Small fixtures fit completely. Native samples can exceed the bounded LRU.
        if (! native) { check(second.work.decodes == first.work.decodes && second.hits > first.hits, "Repeated region extraction decoded again"); }
        check(second.budget.peak_total_bytes <= second.budget.usable_bytes, "Region extraction exceeded admission");
        OPENMS_LOG_INFO << "region-consumer scans=" << source->getNrSpectra() << " ppm=" << ppm << " im=" << use_im
                        << " calibrated=" << calibration.has_value() << " decodes=" << second.work.decodes << " copies=" << second.work.copies
                        << " peak=" << second.budget.peak_total_bytes << '\n';
      }
    }
    // A changed coefficient under the same namespace cannot reuse calibrated columns.
    const RecordContext context {0, 4, 9, 1};
    const auto before = cache->snapshot().work.decodes;
    auto pinned = access->acquireRegion(0, context, MZCalibration {0, 1, 0, false}).scan(0);
    auto shifted = access->acquireRegion(0, context, MZCalibration {1, 1, 0, false}).scan(0);
    check(shifted.mz()[0] == pinned.mz()[0] + 1 && cache->snapshot().work.decodes == before + 2, "Calibration cache key is incomplete");
    auto legacy = access->acquire(0);
    check(legacy->array(0).data[0] == pinned.mz()[0], "Legacy cache received calibrated masses");
    legacy.reset();
    cache->clear();
    check(cache->snapshot().budget.committed_bytes > 0 && pinned.mz().size() > 0, "Cache removal revoked a region lease");
    pinned = {};
    shifted = {};
    check(cache->snapshot().budget.committed_bytes == 0 && cache->snapshot().budget.reserved_bytes == 0, "Region consumer retained charges");
    cache->cancel();
    rejects([&] { access->acquireRegion(0, context); });
    check(source->getSpectrumById(0)->getMZArray()->data == uncalibrated->getMZArray()->data, "Calibration mutated the native source");
  }
  check(signal_observed, "Consumer comparison contained only zero signal");
}
void fixtures()
{
  for (const auto* file :
       {"topp/OpenSwathWorkflow_1_input.mzML", "topp/OpenSwathWorkflow_23_input.mzML", "class_tests/openms/data/IM_FAIMS_test.mzML"})
  {
    OpenMS::MSExperiment loaded;
    OpenMS::MzMLFile().load(std::string(OPENDIA_TEST_DATA_DIR) + "/" + file, loaded);
    const auto path = (std::filesystem::current_path() / ("OpenDIARegionWorkflow_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
    struct Remove
    {
      std::string path;
      ~Remove()
      {
        std::error_code error;
        std::filesystem::remove(path, error);
      }
    } cleanup {path};
    OpenMS::MzMLFile writer;
    writer.getOptions().setWriteIndex(true);
    writer.store(path, loaded);
    auto experiment = std::make_shared<OpenMS::OnDiscMSExperiment>();
    check(experiment->openFile(path), "Cannot open indexed fixture");
    auto source = std::make_shared<OpenMS::SpectrumAccessOnDisc>(experiment);
    exercise(source);
    auto cache = std::make_shared<RawCache>(16 * MIB, 0, 8 * MIB);
    auto input = std::make_shared<RawAccess>(source, cache);
    RegionWorkflow workflow(false, false, false, false, false, 1);
    std::vector<Coordinate> coordinates(1);
    coordinates[0].mz = 400;
    coordinates[0].rt_end = -1;
    OpenMS::ChromExtractParams cp {};
    cp.mz_extraction_window = 0.05;
    cp.im_extraction_window = -1;
    cp.extraction_function = "tophat";
    auto output = empty(1);
    auto invalid_coordinates = coordinates;
    invalid_coordinates[0].mz = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { workflow.extractChromatograms_(input, output, invalid_coordinates, cp); });
    invalid_coordinates[0].mz = std::numeric_limits<double>::max();
    auto wide = cp;
    wide.ppm = true;
    wide.mz_extraction_window = std::numeric_limits<double>::max();
    rejects([&] { workflow.extractChromatograms_(input, output, invalid_coordinates, wide); });
    auto invalid_context = RecordContext {};
    invalid_context.schema_version = 1;
    rejects([&] { input->acquireRegion(0, invalid_context); });
    rejects([&] { input->acquireRegion(0, {}, MZCalibration {std::numeric_limits<double>::quiet_NaN(), 1, 0}); });
    check(output[0]->getTimeArray()->data.empty() && cache->snapshot().work.decodes == 0, "Invalid query reached decoder/output");
    // Non-native access keeps the shared extractor, including its errors.
    auto expected = empty(1);
    OpenMS::ChromatogramExtractor().extractChromatograms(source, expected, coordinates, 0.05, false, -1, "tophat");
    workflow.extractChromatograms_(source, output, coordinates, cp);
    equal(output, expected);
    cp.extraction_function = "not-a-filter";
    rejects([&] { workflow.extractChromatograms_(input, output, coordinates, cp); });
  }
}

void edges()
{
  OpenMS::MSExperiment experiment;
  OpenMS::MSSpectrum spectrum;
  spectrum.setRT(10);
  spectrum.setMSLevel(2);
  spectrum.setNativeID("scan=1");
  const std::array<double, 6> mz {{99, 99.25, 99.5, 100, 100.5, 101}};
  const std::array<float, 6> intensity {{1000, 1, -1e20f, 1e20f, 4, 2000}};
  for (OpenMS::Size p = 0; p < mz.size(); ++p)
  {
    spectrum.emplace_back(mz[p], intensity[p]);
  }
  experiment.addSpectrum(spectrum);
  const auto path = (std::filesystem::current_path() / ("OpenDIARegionEdges_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
  struct Remove
  {
    std::string path;
    ~Remove()
    {
      std::error_code error;
      std::filesystem::remove(path, error);
    }
  } cleanup {path};
  OpenMS::MzMLFile writer;
  writer.getOptions().setWriteIndex(true);
  writer.store(path, experiment);
  auto native = std::make_shared<OpenMS::OnDiscMSExperiment>();
  check(native->openFile(path), "Cannot open endpoint fixture");
  auto source = std::make_shared<OpenMS::SpectrumAccessOnDisc>(native);
  auto cache = std::make_shared<RawCache>(MIB, 0, MIB);
  auto access = std::make_shared<RawAccess>(source, cache);
  RegionWorkflow workflow(false, false, false, false, false, 1);
  std::vector<Coordinate> coordinates(1);
  coordinates[0].mz = 100;
  coordinates[0].rt_start = 10;
  coordinates[0].rt_end = 11;
  coordinates[0].ion_mobility = -1;
  OpenMS::ChromExtractParams cp {};
  cp.mz_extraction_window = 2;
  cp.im_extraction_window = -1;
  cp.extraction_function = "tophat";
  auto actual = empty(1), expected = empty(1);
  OpenMS::ChromatogramExtractor().extractChromatograms(source, expected, coordinates, 2, false, -1, "tophat");
  workflow.extractChromatograms_(access, actual, coordinates, cp);
  equal(actual, expected);
  check(actual[0]->getIntensityArray()->data == std::vector<double> {5}, "Open endpoints or floating-point addition order changed");
  cp.im_extraction_window = 0.1;
  auto absent_im = empty(1);
  rejects([&] { workflow.extractChromatograms_(access, absent_im, coordinates, cp); });
  check(absent_im[0]->getTimeArray()->data.empty(), "Missing point IM produced output");
  auto pinned = access->acquireRegion(0, {}).scan(0);
  cache->cancel();
  check(pinned.mz().size() == mz.size(), "Cancellation revoked region lease");
  pinned = {};
  check(cache->snapshot().budget.committed_bytes == 0, "Cancelled region retained bytes");
  const auto bytes = RawRegion::requiredBytes({1, mz.size()}, spectrum.getNativeID().size());
  auto tiny_cache = std::make_shared<RawCache>(bytes - 1, 0, bytes - 1);
  RawAccess tiny(source, tiny_cache);
  rejects([&] { tiny.acquireRegion(0, {}); });
  check(tiny_cache->snapshot().budget.committed_bytes == 0 && tiny_cache->snapshot().budget.reserved_bytes == 0,
        "Failed region admission leaked bytes");
}
} // namespace
int main(int argc, char** argv)
{
  try
  {
    if (argc == 1)
    {
      fixtures();
      edges();
    }
    else
    {
      for (int argument = 1; argument < argc; ++argument)
      {
        std::shared_ptr<OpenMS::ExperimentalSettings> settings;
        auto maps = OpenMS::SwathFile().loadOnDisc(argv[argument], settings);
        check(! maps.empty(), "Native loader returned no lanes");
        for (const auto map_index : {OpenMS::Size {0}, maps.size() / 2, maps.size() - 1})
        {
          auto& map = maps[map_index];
          auto source = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr);
          check(source != nullptr, "Native loader did not return native access");
          const auto count = source->getNrSpectra();
          if (! count) { continue; }
          std::vector<OpenMS::Size> indices {source->getExperimentIndex(0), source->getExperimentIndex(count / 2),
                                             source->getExperimentIndex(count - 1)};
          auto sampled = std::make_shared<OpenMS::SpectrumAccessOnDisc>(
            std::shared_ptr<OpenMS::IMSExperimentRandomAccess>(source->getExperiment().clone()), indices);
          exercise(sampled, true);
        }
      }
    }
    OPENMS_LOG_INFO << "PASS calibrated region production consumer, exact chromatograms, cache isolation and lease accounting\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
}
