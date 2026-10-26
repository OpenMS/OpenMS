// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIARawAccess.h"
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
#include <cmath>
#include <exception>
#include <filesystem>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
  constexpr OpenDIA::ByteCount MIB = 1024 * 1024;
  void check(bool condition, const char* message)
  {
    if (!condition) { throw std::runtime_error(message); }
  }

  struct Input
  {
    std::vector<OpenSwath::SwathMap> maps;
    std::shared_ptr<OpenDIA::RawCache> cache = std::make_shared<OpenDIA::RawCache>(4 * MIB, MIB, 1);
    std::string indexed_path;
    Input(bool fail_ms1, bool fail_ms2, bool multiple_maps = false)
    {
      std::string path = std::string(OPENDIA_TEST_DATA_DIR) + "/topp/OpenSwathWorkflow_17_input.mzML";
      if (multiple_maps)
      {
        OpenMS::MSExperiment experiment;
        OpenMS::MzMLFile().load(std::string(OPENDIA_TEST_DATA_DIR) + "/topp/OpenSwathWorkflow_1_input.mzML", experiment);
        indexed_path = (std::filesystem::current_path() / ("OpenDIARawFailure_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
        OpenMS::MzMLFile writer;
        writer.getOptions().setWriteIndex(true);
        writer.store(indexed_path, experiment);
        path = indexed_path;
      }
      std::shared_ptr<OpenMS::ExperimentalSettings> metadata;
      maps = OpenMS::SwathFile().loadOnDisc(path, metadata);
      for (auto& map : maps)
      {
        if ((map.ms1 && fail_ms1) || (!map.ms1 && fail_ms2))
        {
          auto native = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr);
          check(native != nullptr, "Expected indexed native source");
          map.sptr = std::make_shared<OpenDIA::RawAccess>(native, cache);
        }
      }
    }
    ~Input()
    {
      if (!indexed_path.empty())
      {
        std::error_code error;
        std::filesystem::remove(indexed_path, error);
      }
    }
  };

  template <typename F> void rejects(Input& input, F operation)
  {
    bool rejected = false;
    try { operation(); }
    catch (const OpenMS::Exception::InvalidValue& error)
    {
      OPENMS_LOG_INFO << "Reader error: " << error.what() << '\n';
      rejected = std::string(error.what()).find("exceeds maximum scan bytes") != std::string::npos;
    }
    check(rejected, "Native parallel path did not return the raw admission error");
    const auto state = input.cache->snapshot();
    check(state.misses > 0, "Raw reader was not called");
    check(state.cache_entries == 0 && state.budget.committed_bytes == 0 && state.budget.reserved_bytes == 0,
      "Rejected scan retained a partial result or charge");
  }

  OpenSwath::LightTargetedExperiment library()
  {
    OpenSwath::LightTargetedExperiment result;
    OpenSwath::LightCompound compound;
    compound.id = "pep";
    compound.sequence = "PEPTIDEA";
    compound.charge = 2;
    compound.rt = 25.2;
    compound.drift_time = 0.2;
    result.compounds.push_back(compound);
    OpenSwath::LightTransition transition;
    transition.precursor_mz = 412.5;
    transition.product_mz = 500;
    transition.precursor_im = 0.2;
    transition.peptide_ref = "pep";
    transition.transition_name = "tr";
    result.transitions.push_back(transition);
    return result;
  }

  void calibration(bool fail_ms1, bool fail_ms2, bool ms1_calibration)
  {
    OPENMS_LOG_INFO << "Calibration failure check: MS1=" << fail_ms1 << " MS2=" << fail_ms2
      << " calibrate_MS1=" << ms1_calibration << '\n';
    Input input(fail_ms1, fail_ms2);
    auto assay = library();
    OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType group;
    group.addTransition(assay.transitions.front(), "tr");
    OpenMS::MRMFeature feature;
    feature.setRT(25.2);
    feature.setOverallQuality(1.0);
    group.addFeature(feature);
    std::map<std::string, OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType*> groups{{"pep", &group}};
    OpenMS::SwathMapMassCorrection correction;
    auto parameters = correction.getParameters();
    parameters.setValue("mz_extraction_window", 0.05);
    parameters.setValue("im_extraction_window", 0.05);
    parameters.setValue("im_correction_function", "linear");
    parameters.setValue("ms1_im_calibration", ms1_calibration ? "true" : "false");
    correction.setParameters(parameters);
    OpenMS::TransformationDescription transformation;
    rejects(input, [&] { correction.correctIM(groups, assay, input.maps, false, transformation); });
  }

  void evidenceFilter(bool ms1, bool multiple_maps, int threads)
  {
    OPENMS_LOG_INFO << "Evidence filter failure check: MS1=" << ms1
      << " multiple_maps=" << multiple_maps << " threads=" << threads << '\n';
    Input input(ms1, !ms1, multiple_maps);
    std::vector<OpenSwath::SwathMap> maps;
    for (const auto& map : input.maps)
    {
      if (map.ms1 == ms1) { maps.push_back(map); }
    }
    check(!maps.empty(), "Expected native evidence map");
    if (!multiple_maps) { maps.resize(1); }
    else { check(maps.size() > 1, "Expected multiple native evidence maps"); }
    auto assay = library();
    // Put a precursor with a detecting fragment inside each isolation window.
    // An empty product index would skip the native read.
    assay.transitions.front().precursor_mz = ms1 ? 412.5 : (maps.front().lower + maps.front().upper) * 0.5;
    for (OpenMS::Size i = 1; i < maps.size(); ++i)
    {
      auto transition = assay.transitions.front();
      transition.precursor_mz = (maps[i].lower + maps[i].upper) * 0.5;
      transition.peptide_ref = "pep" + std::to_string(i);
      transition.transition_name = "tr" + std::to_string(i);
      auto compound = assay.compounds.front();
      compound.id = transition.peptide_ref;
      assay.compounds.push_back(compound);
      assay.transitions.push_back(transition);
    }
    OpenMS::TransitionListEvidenceFilter filter;
    auto parameters = filter.getParameters();
    parameters.setValue("enabled", "true");
    parameters.setValue("evidence_sources", ms1 ? "ms1" : "ms2");
    parameters.setValue("ms2_min_fragment_hits", 1);
    parameters.setValue("min_supported_precursors", 1);
    filter.setParameters(parameters);
    OpenMS::ChromExtractParams extraction_parameters{0, 0.05, -1, false, "tophat", -1, 0};
    rejects(input, [&] { filter.filter(maps, assay, extraction_parameters, extraction_parameters, false, threads); });
  }

  void extraction(int batch, bool load_into_memory)
  {
    OPENMS_LOG_INFO << "Extraction failure check: batch=" << batch << " load=" << load_into_memory << '\n';
    Input input(false, true);
    auto assay = library();
    OpenMS::OpenSwathWorkflow workflow(false, false, false, false, false, -1);
    OpenMS::ChromExtractParams parameters{0, 0.05, -1, false, "tophat", -1, 0};
    OpenMS::FeatureMap result;
    OpenMS::OpenSwathOSWWriter writer("");
    rejects(input, [&]
    {
      workflow.performExtraction(input.maps, OpenMS::TransformationDescription(), parameters, parameters,
        OpenMS::MRMFeatureFinderScoring().getParameters(), assay, result, false, writer, nullptr, batch, 0, load_into_memory);
    });
  }

  void scoringOwnership(int path)
  {
    const char* names[] = {"all", "ranking", "remainder", "precursor", "ID", "fetch-access", "fetch-maps", "MS1-IM-reuse", "exception"};
    OPENMS_LOG_INFO << "Native scoring ownership check: " << names[path] << '\n';
    Input input(false, false);
    auto assay = library();
    assay.transitions.front().fragment_charge = 1;
    assay.transitions.front().library_intensity = 1000;
    if (path == 7)
    {
      auto transition = assay.transitions.front();
      transition.transition_name = "tr2";
      transition.product_mz += 1.0033548;
      assay.transitions.push_back(transition);
    }
    if (path == 8) { assay.compounds.front().sequence = "?"; }
    OpenMS::MSExperiment experiment;
    for (int level : {1, 2})
    {
      OpenMS::MSSpectrum spectrum;
      spectrum.setMSLevel(level);
      spectrum.setRT(25.2);
      spectrum.setNativeID("scan=" + std::to_string(level));
      if (level == 2)
      {
        OpenMS::Precursor precursor;
        precursor.setMZ(412.5);
        precursor.setIsolationWindowLowerOffset(20);
        precursor.setIsolationWindowUpperOffset(20);
        spectrum.getPrecursors().push_back(precursor);
      }
      OpenMS::DataArrays::FloatDataArray im;
      OpenMS::IMDataConverter::setIMUnit(im, OpenMS::DriftTimeUnit::MILLISECOND);
      im.setName(OpenMS::Constants::UserParam::ION_MOBILITY);
      for (int isotope = 0; isotope < 3; ++isotope)
      {
        OpenMS::Peak1D peak;
        peak.setMZ((level == 1 ? 412.5 : 500.0) + 0.001 + isotope * 1.0033548 / (level == 1 ? 2 : 1));
        for (int point = 0; point < (path == 7 ? 3 : 1); ++point)
        {
          peak.setIntensity((point == 0 ? 250.0 : point == 1 ? 1000.0 : 500.0) / (isotope + 1));
          spectrum.push_back(peak);
          im.push_back(path == 7 ? 0.18 + point * 0.02 : 0.2);
        }
      }
      spectrum.getFloatDataArrays().push_back(im);
      experiment.addSpectrum(spectrum);
    }
    input.indexed_path = (std::filesystem::current_path() / ("OpenDIARawFailure_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
    OpenMS::MzMLFile writer;
    writer.getOptions().setWriteIndex(true);
    writer.store(input.indexed_path, experiment);
    std::shared_ptr<OpenMS::ExperimentalSettings> metadata;
    input.maps = OpenMS::SwathFile().loadOnDisc(input.indexed_path, metadata);
    check(input.maps.size() == 2, "Expected indexed scoring MS1 and MS2 maps");

    auto run = [&](const std::vector<OpenSwath::SwathMap>& maps)
    {
      OpenSwath::SpectrumAccessPtr ms1;
      std::vector<OpenSwath::SwathMap> fragments;
      for (const auto& map : maps)
      {
        if (map.ms1) { ms1 = map.sptr; }
        else { fragments.push_back(map); }
      }
      check(ms1 && fragments.size() == 1, "Missing native scoring maps");
      OpenMS::MRMFeature feature;
      feature.setRT(25.2);
      feature.setIntensity(1000);
      OpenMS::Feature fragment;
      fragment.setIntensity(1000);
      feature.addFeature(fragment, "tr");
      if (path == 7) { feature.addFeature(fragment, "tr2"); }
      OpenMS::MRMFeatureOpenMS access(feature);
      OpenMS::OpenSwathScoring scorer;
      OpenMS::OpenSwath_Scores_Usage usage;
      usage.use_im_scores = path == 7;
      usage.use_ionseries_scores = path != 8;
      scorer.initialize(1.0, 1, 0.005, 0.2, 0.0, usage, path == 5 || path == 6 ? "resample" : "simple", "fixed", path == 7, false);
      OpenMS::DIAScoring dia;
      OpenMS::OpenSwath_Scores scores;
      OpenMS::RangeMobility im_range;
      if (path == 7) { im_range = OpenMS::RangeMobility(0.15, 0.25); }
      std::vector<double> errors;
      if (path == 5 || path == 6)
      {
        // Keep the IM array. The one-scan merge can return the raw copy.
        auto spectra = path == 5 ? scorer.fetchSpectrumSwath(fragments.front().sptr, 25.2, 1, im_range)
                                 : scorer.fetchSpectrumSwath(fragments, 25.2, 1, im_range);
        check(spectra.size() == 1, "Missing fetched scoring spectrum");
        auto result = spectra.front()->getMZArray()->data;
        const auto& intensities = spectra.front()->getIntensityArray()->data;
        result.insert(result.end(), intensities.begin(), intensities.end());
        check(result.size() == 6 && result[3] > 0, "Fetched spectrum has no signal");
        return result;
      }
      if (path == 3)
      {
        scorer.calculatePrecursorDIAScores(ms1, dia, 412.5, 25.2, assay.compounds.front(), im_range, scores);
      }
      else if (path == 4)
      {
        OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType group;
        group.addTransition(assay.transitions.front(), "tr");
        scorer.calculateDIAIdScores(&access, assay.transitions.front(), group, fragments, im_range, dia, scores, -1, nullptr, -1);
      }
      else
      {
        const auto subset = path == 1 ? OpenMS::OpenSwathScoring::ScoreSubset::RANKING
          : path == 2 ? OpenMS::OpenSwathScoring::ScoreSubset::REMAINDER : OpenMS::OpenSwathScoring::ScoreSubset::ALL;
        bool rejected = false;
        try
        {
          scorer.calculateDIAScores(&access, assay.transitions, std::vector<double>(assay.transitions.size(), 1.0 / assay.transitions.size()), fragments, ms1, dia,
            assay.compounds.front(), scores, errors, 0.2, im_range, nullptr, -1, subset);
        }
        catch (const OpenMS::Exception::ParseError&)
        {
          if (path != 8) { throw; }
          rejected = true;
        }
        // Invalid peptide text must fail after both raw scans have been read.
        check(rejected == (path == 8), "Missing native precursor sequence parse error");
      }
      if (path != 2 && path != 3) { check(scores.massdev_score > 0, "Fragment mass score has no signal"); }
      if (path == 0 || path == 2 || path == 3) { check(scores.ms1_ppm_score > 0, "Precursor mass score has no signal"); }
      std::vector<double> result{scores.massdev_score, scores.weighted_massdev_score, scores.isotope_correlation,
        scores.isotope_overlap, scores.dotprod_score_dia, scores.manhatt_score_dia, scores.bseries_score,
        scores.yseries_score, scores.ms1_ppm_score, scores.ms1_isotope_correlation, scores.ms1_isotope_overlap};
      result.insert(result.end(), errors.begin(), errors.end());
      if (path == 7)
      {
        check(scores.im_drift > 0 && scores.im_ms1_drift > 0, "Missing native MS2 or MS1 IM evidence");
        result.insert(result.end(), {scores.im_drift, scores.im_drift_weighted, scores.im_delta_score,
          scores.im_ms1_drift, scores.im_ms1_delta_score, scores.im_ms1_contrast_coelution,
          scores.im_ms1_contrast_shape, scores.im_ms1_sum_contrast_coelution, scores.im_ms1_sum_contrast_shape});
      }
      return result;
    };
    const auto reference = run(input.maps);
    input.cache = std::make_shared<OpenDIA::RawCache>(4 * MIB, MIB, MIB);
    for (auto& map : input.maps)
    {
      auto native = std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr);
      check(native != nullptr, "Expected indexed scoring source");
      map.sptr = std::make_shared<OpenDIA::RawAccess>(native, input.cache);
    }
    const auto bounded = run(input.maps);
    check(reference.size() == bounded.size(), "Native scoring result size changed");
    for (OpenMS::Size i = 0; i < reference.size(); ++i)
    {
      if (!std::isfinite(bounded[i]) || std::abs(reference[i] - bounded[i]) >= 1e-10)
      {
        OPENMS_LOG_ERROR << "Scoring value path=" << names[path] << " field=" << i
          << " native=" << reference[i] << " raw=" << bounded[i] << '\n';
      }
      check(std::isfinite(bounded[i]) && std::abs(reference[i] - bounded[i]) < 1e-10, "Raw access changed native scoring values");
    }
    check(input.cache->snapshot().misses > 0, "Scoring did not read raw scans");
    input.cache->clear();
    const auto state = input.cache->snapshot();
    OPENMS_LOG_INFO << "Scoring ownership path=" << names[path] << " native_reference=equal committed="
      << state.budget.committed_bytes << " reserved=" << state.budget.reserved_bytes << '\n';
    check(state.cache_entries == 0 && state.budget.committed_bytes == 0 && state.budget.reserved_bytes == 0,
      "Native scoring retained a final charge");
  }

  void calibrationOwnership(bool mobility)
  {
    OPENMS_LOG_INFO << "Native calibration ownership check: IM=" << mobility << '\n';
    Input input(false, false);
    auto assay = library();
    assay.compounds.clear();
    assay.transitions.clear();
    using Group = OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType;
    std::vector<Group> storage(4);
    std::map<std::string, Group*> groups;
    OpenMS::MSExperiment experiment;
    for (int level : {1, 2})
    {
      OpenMS::MSSpectrum spectrum;
      spectrum.setMSLevel(level);
      spectrum.setRT(25.2);
      spectrum.setNativeID("scan=" + std::to_string(level));
      if (level == 2)
      {
        OpenMS::Precursor precursor;
        precursor.setMZ(412.5);
        precursor.setIsolationWindowLowerOffset(20);
        precursor.setIsolationWindowUpperOffset(20);
        spectrum.getPrecursors().push_back(precursor);
      }
      OpenMS::DataArrays::FloatDataArray im;
      OpenMS::IMDataConverter::setIMUnit(im, OpenMS::DriftTimeUnit::MILLISECOND);
      im.setName(OpenMS::Constants::UserParam::ION_MOBILITY);
      for (OpenMS::Size i = 0; i < storage.size(); ++i)
      {
        const auto id = "cal" + std::to_string(i);
        const double target_im = 0.5 + 0.1 * i;
        OpenMS::Peak1D peak;
        peak.setMZ((level == 1 ? 410.0 : 500.0) + i + 0.001 * (i + 1));
        peak.setIntensity(128 + i);
        spectrum.push_back(peak);
        im.push_back(target_im + (level == 1 ? 0.003 : 0.002) * (i % 2 + 1));
        if (level == 2)
        {
          auto compound = library().compounds.front();
          compound.id = id;
          compound.drift_time = target_im;
          assay.compounds.push_back(compound);
          auto transition = library().transitions.front();
          transition.peptide_ref = id;
          transition.transition_name = id;
          transition.precursor_mz = 410.0 + i;
          transition.product_mz = 500.0 + i;
          transition.precursor_im = target_im;
          assay.transitions.push_back(transition);
          storage[i].addTransition(transition, id);
          OpenMS::MRMFeature feature;
          feature.setRT(25.2);
          feature.setOverallQuality(1.0);
          storage[i].addFeature(feature);
          groups[id] = &storage[i];
        }
      }
      spectrum.getFloatDataArrays().push_back(im);
      experiment.addSpectrum(spectrum);
    }
    input.indexed_path = (std::filesystem::current_path() / ("OpenDIARawFailure_" + OpenMS::File::getUniqueName(false) + ".mzML")).string();
    OpenMS::MzMLFile writer;
    writer.getOptions().setWriteIndex(true);
    writer.store(input.indexed_path, experiment);
    std::shared_ptr<OpenMS::ExperimentalSettings> metadata;
    input.maps = OpenMS::SwathFile().loadOnDisc(input.indexed_path, metadata);
    check(input.maps.size() == 2, "Expected native MS1 and MS2 maps");

    // Measure one native scan and its charged mutable copy, without a guessed size.
    OpenDIA::ByteCount cap = 0;
    input.cache = std::make_shared<OpenDIA::RawCache>(4 * MIB, MIB, MIB);
    for (const auto& map : input.maps)
    {
      OpenDIA::RawAccess reader(std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr), input.cache);
      auto copy = reader.getSpectrumById(0);
      cap = std::max(cap, input.cache->snapshot().budget.committed_bytes);
      copy.reset();
      input.cache->clear();
    }
    auto run = [&](std::vector<OpenSwath::SwathMap> maps)
    {
      OpenMS::SwathMapMassCorrection correction;
      auto parameters = correction.getParameters();
      parameters.setValue("mz_extraction_window", 0.05);
      parameters.setValue("im_extraction_window", 0.05);
      parameters.setValue("mz_correction_function", "unweighted_regression");
      correction.setParameters(parameters);
      OpenMS::TransformationDescription transformation;
      if (mobility) { correction.correctIM(groups, assay, maps, false, transformation); }
      else { correction.correctMZ(groups, assay, maps, false); }
      std::vector<double> result;
      if (mobility)
      {
        check(transformation.getDataPoints().size() == storage.size(), "Missing native IM calibrants");
        for (double value : {0.5, 0.8}) { result.push_back(transformation.apply(value)); }
        result.push_back(correction.getFragmentImWindow());
        result.push_back(correction.getPrecursorImWindow());
      }
      else
      {
        for (auto& map : maps)
        {
          if (!map.ms1)
          {
            auto spectrum = map.sptr->getSpectrumById(0);
            result = spectrum->getMZArray()->data;
            check(result.size() == storage.size(), "Missing native m/z calibrants");
            check(std::abs(result.front() - 500.0) < 1e-6, "Native m/z fit did not correct fragments");
          }
        }
        result.push_back(correction.getFragmentMzWindow());
        result.push_back(correction.getPrecursorMzWindow());
      }
      return result;
    };
    const auto reference = run(input.maps);
    check(reference[reference.size() - 2] > 0 && reference.back() > 0, "Native calibration windows have no signal");
    input.cache = std::make_shared<OpenDIA::RawCache>(cap + MIB, MIB, cap);
    for (auto& map : input.maps)
    {
      map.sptr = std::make_shared<OpenDIA::RawAccess>(std::dynamic_pointer_cast<OpenMS::SpectrumAccessOnDisc>(map.sptr), input.cache);
      // Each scan fits with its copy, but any live copy blocks the next decode reservation.
      auto copy = map.sptr->getSpectrumById(0);
      input.cache->clear();
      check(input.cache->snapshot().budget.committed_bytes > 0, "Mutable copy has no charge");
      copy.reset();
    }
    const auto bounded = run(input.maps);
    check(reference.size() == bounded.size(), "Calibration result size changed");
    for (OpenMS::Size i = 0; i < reference.size(); ++i)
    {
      check(std::isfinite(bounded[i]) && std::abs(reference[i] - bounded[i]) < 1e-10, "Bounded calibration changed fit or window");
    }
    const auto state = input.cache->snapshot();
    check(state.misses >= 2 && state.evictions > 0, "Calibration did not exercise bounded native reads");
    check(state.budget.peak_total_bytes <= state.budget.usable_bytes, "Calibration exceeded byte cap");
    input.cache->clear();
    check(input.cache->snapshot().budget.committed_bytes == 0 && input.cache->snapshot().budget.reserved_bytes == 0,
      "Calibration retained a final charge");
    OPENMS_LOG_INFO << "Calibration ownership IM=" << mobility << " usable_bytes=" << state.budget.usable_bytes
      << " peak_bytes=" << state.budget.peak_total_bytes << " native_reference=equal final_charge=0\n";
  }
}

int main()
{
  try
  {
#ifdef _OPENMP
    const int saved_threads = omp_get_max_threads();
    const int saved_dynamic = omp_get_dynamic();
    struct RestoreThreads
    {
      int threads, dynamic;
      ~RestoreThreads() { omp_set_num_threads(threads); omp_set_dynamic(dynamic); }
    } restore_threads{saved_threads, saved_dynamic};
    omp_set_dynamic(0);
    omp_set_num_threads(2);
#endif
    // Use the actual indexed reader. Each native read exceeds its declared one-byte limit.
    calibration(false, true, false);
    calibration(true, false, false);
    calibration(true, false, true);
    evidenceFilter(true, false, 2);
    // The public single-map MS2 path stays serial even with two requested threads.
    evidenceFilter(false, false, 2);
    // The indexed multi-window fixture selects the parallel wave path.
    evidenceFilter(false, true, 2);
    evidenceFilter(true, false, 1);
    evidenceFilter(false, false, 1);
    evidenceFilter(false, true, 1);
    extraction(0, false);
    extraction(1, false);
#ifdef _OPENMP
    // One SWATH and one scoring thread select the parallel wave extraction branch.
    omp_set_num_threads(1);
#endif
    extraction(0, true);
    OPENMS_LOG_INFO << "Returned twelve native raw admission failures without parallel termination\n";
    std::exception_ptr failure;
    for (bool mobility : {true, false})
    {
      try { calibrationOwnership(mobility); }
      catch (const std::exception& error)
      {
        OPENMS_LOG_ERROR << "Calibration ownership IM=" << mobility << ": " << error.what() << '\n';
        if (!failure) { failure = std::current_exception(); }
      }
    }
    for (int path = 0; path < 9; ++path)
    {
      try { scoringOwnership(path); }
      catch (const std::exception& error)
      {
        OPENMS_LOG_ERROR << "Scoring ownership path=" << path << ": " << error.what() << '\n';
        if (!failure) { failure = std::current_exception(); }
      }
    }
    if (failure) { std::rethrow_exception(failure); }
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
}
