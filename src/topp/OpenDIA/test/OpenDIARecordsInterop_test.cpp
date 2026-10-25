// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARecords.h"

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <array>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
  void check(bool condition, const std::string& message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  std::string dataFile(const std::string& relative)
  {
    return std::string(OPENDIA_TEST_DATA_DIR) + "/" + relative;
  }

  OpenDIA::RecordContext context()
  {
    // Fixture-owned identity namespace, independent of vector position or hashes.
    return {0, 71, 29, 503};
  }

  OpenDIA::MobilityUnit mobilityUnit(OpenMS::DriftTimeUnit unit)
  {
    if (unit == OpenMS::DriftTimeUnit::MILLISECOND)
    {
      return OpenDIA::MobilityUnit::MILLISECOND;
    }
    check(unit == OpenMS::DriftTimeUnit::VSSC, "Unsupported continuous source mobility unit");
    return OpenDIA::MobilityUnit::INVERSE_REDUCED;
  }

  // This is a contract-boundary fixture, not a production raw adapter or extractor.
  OpenDIA::RawTile produceRaw(const std::string& filename)
  {
    OpenMS::MSExperiment source;
    OpenMS::MzMLFile().load(filename, source);
    std::vector<OpenDIA::RawPoint> points;
    for (OpenMS::Size scan = 0; scan < source.size(); ++scan)
    {
      const auto& spectrum = source[scan];
      check(spectrum.getPrecursors().size() <= 1, "Fixture requires at most one isolation window per spectrum");
      OpenDIA::Coordinates coordinates;
      coordinates.rt_seconds = spectrum.getRT();
      coordinates.native_id = spectrum.getNativeID();
      // mzML has no universal numeric frame/scan IDs. These are caller IDs only.
      coordinates.frame_id = scan + 1000;
      coordinates.scan_id = scan + 2000;
      coordinates.ms_level = spectrum.getMSLevel();
      if (!spectrum.getPrecursors().empty())
      {
        const auto& precursor = spectrum.getPrecursors().front();
        coordinates.isolation_lower_mz = precursor.getMZ() - precursor.getIsolationWindowLowerOffset();
        coordinates.isolation_upper_mz = precursor.getMZ() + precursor.getIsolationWindowUpperOffset();
      }
      const auto unit = spectrum.getDriftTimeUnit();
      if (unit == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
      {
        coordinates.faims_voltage = spectrum.getDriftTime();
      }
      else if (unit != OpenMS::DriftTimeUnit::NONE)
      {
        coordinates.mobility = OpenDIA::Mobility{spectrum.getDriftTime(), mobilityUnit(unit)};
      }
      for (OpenMS::Size peak = 0; peak < spectrum.size(); ++peak)
      {
        coordinates.mz = spectrum[peak].getMZ();
        if (spectrum.containsIMData())
        {
          const auto [index, peak_unit] = spectrum.getIMData();
          const auto& values = spectrum.getFloatDataArrays().at(index);
          check(values.size() == spectrum.size(), "Malformed source mobility array");
          coordinates.mobility = OpenDIA::Mobility{values[peak], mobilityUnit(peak_unit)};
        }
        points.push_back({coordinates, spectrum[peak].getIntensity()});
      }
    }
    const auto count = points.size();
    return OpenDIA::RawTile(context(), {count, count}, std::move(points));
  }

  void checkContext(const OpenDIA::RecordContext& actual)
  {
    const auto expected = context();
    check(actual.schema_version == expected.schema_version && actual.run_id == expected.run_id &&
          actual.lane_id == expected.lane_id && actual.calibration_snapshot_id == expected.calibration_snapshot_id,
          "Record context changed across producer lifetime");
  }

  OpenMS::Size consumeRaw(const OpenDIA::RawTile& tile, const std::string& filename,
                         bool expect_im, bool expect_faims)
  {
    // Reload independently after the producer's experiment and vectors were destroyed.
    OpenMS::MSExperiment reference;
    OpenMS::MzMLFile().load(filename, reference);
    checkContext(tile.context());
    OpenMS::Size row = 0;
    OpenMS::Size im_rows = 0;
    OpenMS::Size faims_rows = 0;
    for (OpenMS::Size scan = 0; scan < reference.size(); ++scan)
    {
      const auto& spectrum = reference[scan];
      for (OpenMS::Size peak = 0; peak < spectrum.size(); ++peak, ++row)
      {
        const auto actual = tile.point(row);
        const auto& coordinates = actual.coordinates;
        check(coordinates.rt_seconds == spectrum.getRT() && coordinates.mz == spectrum[peak].getMZ(),
              "Measured RT/mz changed");
        check(actual.intensity == spectrum[peak].getIntensity(), "Measured intensity changed");
        check(coordinates.native_id == spectrum.getNativeID() && coordinates.frame_id == scan + 1000 &&
              coordinates.scan_id == scan + 2000 && coordinates.ms_level == spectrum.getMSLevel(),
              "Source acquisition identity changed");
        if (spectrum.getPrecursors().empty())
        {
          check(coordinates.isolation_lower_mz == 0 && coordinates.isolation_upper_mz == 0,
                "Unexpected MS1 isolation geometry");
        }
        else
        {
          const auto& precursor = spectrum.getPrecursors().front();
          check(coordinates.isolation_lower_mz == precursor.getMZ() - precursor.getIsolationWindowLowerOffset() &&
                coordinates.isolation_upper_mz == precursor.getMZ() + precursor.getIsolationWindowUpperOffset(),
                "Isolation geometry changed");
        }
        auto expected_unit = spectrum.getDriftTimeUnit();
        double expected_mobility = spectrum.getDriftTime();
        if (spectrum.containsIMData())
        {
          const auto [index, unit] = spectrum.getIMData();
          expected_unit = unit;
          expected_mobility = spectrum.getFloatDataArrays().at(index).at(peak);
        }
        const bool continuous = expected_unit == OpenMS::DriftTimeUnit::MILLISECOND ||
                                expected_unit == OpenMS::DriftTimeUnit::VSSC;
        check(coordinates.mobility.has_value() == continuous, "Continuous mobility missingness changed");
        if (continuous)
        {
          ++im_rows;
          check(coordinates.mobility->value == expected_mobility, "Continuous mobility value changed");
          check((coordinates.mobility->unit == OpenDIA::MobilityUnit::MILLISECOND) ==
                (expected_unit == OpenMS::DriftTimeUnit::MILLISECOND), "Continuous mobility unit changed");
        }
        const bool faims = expected_unit == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE;
        check(coordinates.faims_voltage.has_value() == faims, "FAIMS missingness changed");
        if (faims)
        {
          ++faims_rows;
          check(*coordinates.faims_voltage == expected_mobility && !coordinates.mobility,
                "FAIMS voltage was changed or treated as continuous mobility");
        }
      }
    }
    check(row > 0 && tile.size() == row, "Raw records were dropped or added");
    check((im_rows > 0) == expect_im && (faims_rows > 0) == expect_faims, "Fixture did not exercise expected acquisition mode");
    OPENMS_LOG_INFO << "Preserved " << row << " raw points, " << im_rows << " IM, " << faims_rows
                    << " FAIMS, from " << reference.size() << " source spectra\n";
    return row;
  }

  void consumeQuery(const OpenDIA::RawTile& tile)
  {
    const auto source = tile.point(0).coordinates;
    OpenDIA::ExtractionQuerySpec request;
    request.context = tile.context();
    request.acquisitions = std::vector<std::pair<OpenDIA::RecordId, OpenDIA::RecordId>>{{source.frame_id, source.scan_id}};
    request.ms_level = source.ms_level;
    request.rt_seconds = OpenDIA::QueryInterval{source.rt_seconds, source.rt_seconds, OpenDIA::EdgePolicy::CLOSED};
    request.mz = {source.mz, source.mz, OpenDIA::EdgePolicy::CLOSED};
    if (source.ms_level > 1)
    {
      request.isolation_mz = OpenDIA::QueryInterval{source.isolation_lower_mz, source.isolation_upper_mz,
                                                  OpenDIA::EdgePolicy::LOWER_CLOSED};
    }
    if (source.mobility)
    {
      request.mobility = OpenDIA::MobilityInterval{{source.mobility->value, source.mobility->value, OpenDIA::EdgePolicy::CLOSED},
                                                   source.mobility->unit};
    }
    request.faims_voltage = source.faims_voltage;
    const OpenDIA::ExtractionQuery query(request);
    check(query.specification() == request, "Real acquisition query coordinates or identity changed");
    ++request.context.calibration_snapshot_id;
    check(query != OpenDIA::ExtractionQuery(request), "Changed calibration reused a real acquisition key");
    OPENMS_LOG_INFO << "Preserved real query identity: MS" << source.ms_level
                    << ", IM=" << source.mobility.has_value() << ", FAIMS=" << source.faims_voltage.has_value() << '\n';
  }

  void replayExtractionQuery()
  {
    // Exercise the existing public extractor, not a replacement algorithm or cache.
    auto experiment = std::make_shared<OpenMS::MSExperiment>();
    OpenMS::MzMLFile().load(dataFile("class_tests/openms/data/ChromatogramExtractor_input.mzML"), *experiment);
    auto access = std::make_shared<OpenMS::SpectrumAccessOpenMS>(experiment);
    OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates coordinates;
    coordinates.mz = 618.31;
    coordinates.rt_start = 0;
    coordinates.rt_end = -1;
    coordinates.id = "tr1";
    OpenDIA::ExtractionQuerySpec request;
    request.context = context();
    check(!experiment->empty(), "Public extractor fixture is empty");
    request.ms_level = (*experiment)[0].getMSLevel();
    for (const auto& spectrum : *experiment)
    {
      check(spectrum.getMSLevel() == request.ms_level, "Extractor fixture spans unselected MS levels");
    }
    request.mz = {coordinates.mz - 0.025, coordinates.mz + 0.025, OpenDIA::EdgePolicy::OPEN};
    request.tolerance = {0.05, OpenDIA::ToleranceUnit::TH};
    request.aggregation = OpenDIA::QueryAggregation::SUM;
    const OpenDIA::ExtractionQuery query(request);
    const auto retained = query.specification();
    auto restored_coordinates = coordinates;
    restored_coordinates.mz = (retained.mz.lower + retained.mz.upper) / 2;
    check(restored_coordinates.mz == coordinates.mz, "Query bounds changed extraction center");
    std::vector<OpenSwath::ChromatogramPtr> expected{std::make_shared<OpenSwath::Chromatogram>()};
    std::vector<OpenSwath::ChromatogramPtr> actual{std::make_shared<OpenSwath::Chromatogram>()};
    OpenMS::ChromatogramExtractorAlgorithm extractor;
    extractor.extractChromatograms(access, expected, {coordinates}, 0.05, false, -1, "tophat");
    extractor.extractChromatograms(access, actual, {restored_coordinates}, retained.tolerance.full_width,
                                  retained.tolerance.unit == OpenDIA::ToleranceUnit::PPM, -1, "tophat");
    const auto& intensities = actual.front()->getIntensityArray()->data;
    check(!intensities.empty() && actual.front()->getTimeArray()->data == expected.front()->getTimeArray()->data &&
          intensities == expected.front()->getIntensityArray()->data, "Query round-trip changed public extraction output");
    bool observed_signal = false;
    for (const auto intensity : intensities)
    {
      observed_signal = observed_signal || intensity > 0;
    }
    check(observed_signal, "Public extraction replay produced no measured signal");
    OPENMS_LOG_INFO << "Public extractor query replay preserved " << intensities.size() << " chromatogram samples with measured signal\n";
  }

  // Explicit fixture registry. Production IDs must likewise come from a stable caller namespace.
  OpenDIA::RecordId precursorId(const std::string& id)
  {
    static const std::map<std::string, OpenDIA::RecordId> ids = {
      {"PEPTIDEA", 81}, {"PEPTIDEB", 13}, {"PEPTIDEC", 55},
      {"PEPTIDED", 27}, {"PEPTIDEE", 93}, {"PEPTIDEF", 41}};
    return ids.at(id);
  }

  OpenDIA::RecordId fragmentId(const std::string& id)
  {
    if (id == "tr1")
    {
      return 901;
    }
    if (id == "tr2")
    {
      return 902;
    }
    return 900 + std::stoull(id);
  }

  OpenSwath::LightTargetedExperiment loadLibrary(const std::string& filename)
  {
    OpenSwath::LightTargetedExperiment library;
    OpenMS::TransitionTSVFile().convertTSVToTargetedExperiment(filename.c_str(), OpenMS::FileTypes::TSV, library);
    return library;
  }

  OpenDIA::HypothesisPage produceHypotheses(const std::string& filename)
  {
    const auto source = loadLibrary(filename);
    std::vector<OpenDIA::Hypothesis> hypotheses;
    std::vector<OpenDIA::Fragment> fragments;
    std::vector<OpenMS::Size> offsets{0};
    for (const auto& compound : source.getCompounds())
    {
      OpenDIA::Hypothesis hypothesis;
      hypothesis.family_id = precursorId(compound.id) + 5000;
      hypothesis.precursor_id = precursorId(compound.id);
      hypothesis.provenance.library = filename;
      bool found = false;
      for (const auto& transition : source.getTransitions())
      {
        if (transition.peptide_ref != compound.id)
        {
          continue;
        }
        hypothesis.precursor_mz = transition.precursor_mz;
        hypothesis.target_decoy = transition.getDecoy() ? OpenDIA::TargetDecoy::DECOY : OpenDIA::TargetDecoy::TARGET;
        fragments.push_back({fragmentId(transition.getNativeID()), transition.product_mz, transition.library_intensity});
        found = true;
      }
      check(found, "Fixture compound has no transitions");
      // NormalizedRetentionTime and unit-unspecified library IM are NOT raw seconds or known IM units.
      // Leave these predictions absent until a producer supplies an explicit calibrated reference.
      hypotheses.push_back(std::move(hypothesis));
      offsets.push_back(fragments.size());
    }
    const OpenDIA::RecordLimits limits{hypotheses.size(), fragments.size()};
    return OpenDIA::HypothesisPage(context(), limits, std::move(hypotheses), std::move(offsets), std::move(fragments));
  }

  void consumeHypotheses(const OpenDIA::HypothesisPage& page, const std::string& filename)
  {
    const auto reference = loadLibrary(filename);
    checkContext(page.context());
    check(page.size() == reference.getCompounds().size(), "Hypothesis count changed");
    OpenMS::Size total_fragments = 0;
    for (OpenMS::Size i = 0; i < page.size(); ++i)
    {
      const auto actual = page.hypothesis(i);
      const auto& compound = reference.getCompounds()[i];
      check(actual.precursor_id == precursorId(compound.id) && actual.family_id == precursorId(compound.id) + 5000,
            "Stable hypothesis identity changed");
      check(actual.provenance.library == filename && !actual.provenance.prediction, "Library provenance changed");
      check(!actual.predicted_rt_seconds && !actual.predicted_mobility, "Unknown prediction reference was invented");
      OpenMS::Size count = 0;
      for (const auto& transition : reference.getTransitions())
      {
        if (transition.peptide_ref != compound.id)
        {
          continue;
        }
        const auto fragment = page.fragment(i, count++);
        check(actual.precursor_mz == transition.precursor_mz &&
              (actual.target_decoy == OpenDIA::TargetDecoy::DECOY) == transition.getDecoy(),
              "Precursor coordinate or target/decoy identity changed");
        check(fragment.fragment_id == fragmentId(transition.getNativeID()) && fragment.mz == transition.product_mz &&
              fragment.predicted_intensity && *fragment.predicted_intensity == transition.library_intensity,
              "Fragment identity or measurement changed");
      }
      check(page.fragmentCount(i) == count, "Sparse fragment offsets changed");
      total_fragments += count;
    }
    check(total_fragments == reference.getTransitions().size(), "Library transitions were dropped");
    OPENMS_LOG_INFO << "Preserved " << page.size() << " hypotheses and " << total_fragments << " source transitions\n";
  }

  OpenDIA::EvidenceBlock produceEvidence(const std::string& raw_file, const std::string& library_file)
  {
    const auto tile = produceRaw(raw_file);
    const auto page = produceHypotheses(library_file);
    const auto hypothesis = page.hypothesis(0);
    const auto fragment = page.fragment(0, 0);
    std::vector<OpenDIA::EvidencePoint> points;
    // Transfer sparse observations, not scientific matching/extraction. Ownership alone is under test.
    for (OpenMS::Size i = 0; i < tile.size(); i += 2)
    {
      const auto raw = tile.point(i);
      points.push_back({hypothesis.family_id, hypothesis.precursor_id, fragment.fragment_id,
                        raw.coordinates, raw.intensity});
    }
    const auto count = points.size();
    return OpenDIA::EvidenceBlock(context(), {count, count}, std::move(points));
  }

  void consumeEvidence(const OpenDIA::EvidenceBlock& block, const std::string& raw_file,
                       const std::string& library_file)
  {
    const auto reference = produceRaw(raw_file);
    const auto page = produceHypotheses(library_file);
    checkContext(block.context());
    check(block.size() == (reference.size() + 1) / 2, "Sparse evidence support changed");
    for (OpenMS::Size i = 0; i < block.size(); ++i)
    {
      const auto actual = block.point(i);
      const auto raw = reference.point(i * 2);
      const auto hypothesis = page.hypothesis(0);
      check(actual.family_id == hypothesis.family_id && actual.precursor_id == hypothesis.precursor_id &&
            actual.fragment_id == page.fragment(0, 0).fragment_id, "Evidence identity changed");
      check(actual.value && *actual.value == raw.intensity, "Observed evidence value changed");
      const auto& a = actual.coordinates;
      const auto& b = raw.coordinates;
      check(a.rt_seconds == b.rt_seconds && a.mz == b.mz && a.frame_id == b.frame_id && a.scan_id == b.scan_id &&
            a.native_id == b.native_id && a.ms_level == b.ms_level && a.isolation_lower_mz == b.isolation_lower_mz &&
            a.isolation_upper_mz == b.isolation_upper_mz && a.faims_voltage == b.faims_voltage &&
            a.mobility.has_value() == b.mobility.has_value(), "Retained evidence coordinates changed");
      if (a.mobility)
      {
        check(a.mobility->unit == b.mobility->unit && a.mobility->value == b.mobility->value,
              "Retained evidence mobility changed");
      }
    }
    OPENMS_LOG_INFO << "Retained " << block.size() << " sparse evidence rows after raw/page destruction\n";
  }
}

int main()
{
  try
  {
    const auto library = dataFile("topp/OpenSwathWorkflow_17_input.tsv");
    const std::array<std::string, 3> files = {
      dataFile("topp/OpenSwathWorkflow_1_input.mzML"),
      dataFile("topp/OpenSwathWorkflow_23_input.mzML"),
      dataFile("class_tests/openms/data/IM_FAIMS_test.mzML")};
    for (OpenMS::Size i = 0; i < files.size(); ++i)
    {
      const auto lease = produceRaw(files[i]);
      consumeRaw(lease, files[i], i == 1, i == 2);
      consumeQuery(lease);
      const auto evidence = produceEvidence(files[i], library);
      consumeEvidence(evidence, files[i], library);
    }
    consumeHypotheses(produceHypotheses(library), library);
    replayExtractionQuery();
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
  return 0;
}
