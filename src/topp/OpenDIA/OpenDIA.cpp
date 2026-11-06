// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <fstream>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace OpenMS;

/**
  @page TOPP_OpenDIA OpenDIA
  @brief Stream DIA spectra and report targeted fragment-group peak candidates.

  This first implementation reads a prepared PQP transition library and streams
  MS2 spectra through OpenMS file handlers. It reports gap-segmented groups of
  co-observed library fragments to a TSV file.

  This is a minimal extraction and candidate-generation path, not a complete
  identification workflow. It does not perform library generation, machine
  learning scoring, FDR estimation, protein inference, or ion-mobility-native
  scoring. Input spectra are expected to contain precursor isolation metadata.
*/
/// @cond TOPPCLASSES
class TOPPOpenDIA final : public TOPPBase
{
public:
  TOPPOpenDIA() : TOPPBase("OpenDIA", "Stream DIA data and detect targeted fragment-group peaks.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "DIA raw data or mzML/mzXML file");
    registerInputFile_("lib", "<file>", "", "Prepared PQP transition library");
    setValidFormats_("lib", {"pqp"});
    registerOutputFile_("out", "<file>", "", "Streamed fragment-group peak table");
    setValidFormats_("out", {"tsv"});
    registerDoubleOption_("mz_tolerance", "<Da>", 0.02, "Absolute fragment m/z matching tolerance in Da", false);
    registerDoubleOption_("rt_window", "<seconds>", 0.0, "Optional window around library RT; zero disables RT filtering", false);
    registerDoubleOption_("peak_gap", "<seconds>", 30.0, "Maximum gap between matched scans in one peak group", false);
    registerIntOption_("min_transitions", "<n>", 1, "Minimum co-observed fragments required in a peak group", false);
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string input = getStringOption_("in");
    const std::string library = getStringOption_("lib");
    const std::string output_path = getStringOption_("out");

    OpenDIA::SearchOptions options;
    options.fragment_tolerance_da = getDoubleOption_("mz_tolerance");
    options.rt_window_seconds = getDoubleOption_("rt_window");
    options.peak_gap_seconds = getDoubleOption_("peak_gap");
    const Int min_transitions = getIntOption_("min_transitions");
    if (min_transitions < 1 || options.fragment_tolerance_da <= 0.0 ||
        options.rt_window_seconds < 0.0 || options.peak_gap_seconds <= 0.0)
    {
      writeLogError_("m/z tolerance and peak gap must be positive; RT window must be non-negative; min_transitions must be at least 1.");
      return ILLEGAL_PARAMETERS;
    }
    options.min_transitions = static_cast<std::size_t>(min_transitions);

    TargetedExperiment targeted;
    TransitionPQPFile().convertPQPToTargetedExperiment(library.c_str(), targeted);
    auto queries = makeQueries_(targeted);
    if (queries.empty())
    {
      writeLogError_("The PQP library contains no usable precursor-fragment transitions.");
      return ILLEGAL_PARAMETERS;
    }

    std::ofstream output(output_path);
    if (!output)
    {
      writeLogError_("Could not open output file: " + output_path);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    output << "precursor_id\tprecursor_mz\trt_start\trt_apex\trt_end\tarea_score\ttransition_support\tmatched_scans\n";

    OpenDIA::StreamSearch consumer(std::move(queries), options, output);
    FileHandler reader;
    PeakFileOptions read_options;
    read_options.setMSLevels({2});
    reader.setOptions(read_options);
    reader.transformExperiment(input, &consumer, {}, log_type_);
    consumer.finish();
    return output ? EXECUTION_OK : CANNOT_WRITE_OUTPUT_FILE;
  }

private:
  static std::vector<OpenDIA::TransitionQuery> makeQueries_(const TargetedExperiment& targeted)
  {
    struct RTValue
    {
      double seconds = 0.0;
      bool available = false;
    };
    std::unordered_map<std::string, RTValue> retention_times;
    retention_times.reserve(targeted.getPeptides().size() + targeted.getCompounds().size());

    const auto add_rt = [&retention_times](const auto& analyte)
    {
      RTValue value;
      if (analyte.hasRetentionTime())
      {
        const auto unit = analyte.getRetentionTimeUnit();
        if (unit == TargetedExperiment::RetentionTime::RTUnit::SECOND ||
            unit == TargetedExperiment::RetentionTime::RTUnit::MINUTE)
        {
          value.seconds = analyte.getRetentionTime() *
            (unit == TargetedExperiment::RetentionTime::RTUnit::MINUTE ? 60.0 : 1.0);
          value.available = true;
        }
      }
      retention_times.emplace(analyte.id, value);
    };
    for (const auto& peptide : targeted.getPeptides()) add_rt(peptide);
    for (const auto& compound : targeted.getCompounds()) add_rt(compound);

    std::vector<OpenDIA::TransitionQuery> queries;
    queries.reserve(targeted.getTransitions().size());
    for (const auto& transition : targeted.getTransitions())
    {
      const std::string& id = transition.getPeptideRef().empty() ? transition.getCompoundRef() : transition.getPeptideRef();
      const auto rt_it = retention_times.find(id);
      OpenDIA::TransitionQuery query;
      query.precursor_id = id;
      query.precursor_mz = transition.getPrecursorMZ();
      query.fragment_mz = transition.getProductMZ();
      if (rt_it != retention_times.end())
      {
        query.library_rt = rt_it->second.seconds;
        query.has_library_rt = rt_it->second.available;
      }
      queries.push_back(std::move(query));
    }
    return queries;
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIA tool;
  return tool.main(argc, argv);
}
/// @endcond
