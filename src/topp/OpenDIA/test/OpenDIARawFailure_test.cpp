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
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
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
}

int main()
{
  try
  {
#ifdef _OPENMP
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
  }
  catch (const std::exception& error)
  {
    OPENMS_LOG_ERROR << error.what() << '\n';
    return 1;
  }
}
