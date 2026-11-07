// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <cmath>
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
  co-observed library fragments to a TSV file. When both the library and raw
  spectra provide inverse reduced ion mobility, it automatically uses per-peak
  mobility to reject out-of-mobility fragment matches. If either side lacks
  compatible mobility metadata, matching falls back to m/z and isolation windows.

  This is a minimal extraction and candidate-generation path, not a complete
  identification workflow. It does not perform library generation, machine
  learning scoring, FDR estimation, or protein inference. An optional feature
  table can preserve candidate-level target/decoy labels and numeric evidence
  for later rescoring. Input spectra are expected to contain precursor isolation
  metadata.
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
    registerOutputFile_("features_out", "<file>", "", "Optional candidate-level numeric feature table for rescoring", false);
    setValidFormats_("features_out", {"tsv"});
    registerDoubleOption_("mz_tolerance", "<Da>", 0.02, "Absolute fragment m/z matching tolerance in Da", false);
    registerDoubleOption_("im_tolerance", "<1/K0>", 0.03, "Automatic inverse reduced ion-mobility matching tolerance", false);
    registerDoubleOption_("rt_window", "<seconds>", 0.0, "Optional window around library RT; zero disables RT filtering", false);
    registerDoubleOption_("peak_gap", "<seconds>", 30.0, "Maximum gap between matched scans in one peak group", false);
    registerIntOption_("min_transitions", "<n>", 1, "Minimum co-observed fragments required in a peak group", false);
  }

  ExitCodes main_(int, const char**) override
  {
    const std::string input = getStringOption_("in");
    const std::string library = getStringOption_("lib");
    const std::string output_path = getStringOption_("out");
    const std::string features_path = getStringOption_("features_out");

    OpenDIA::SearchOptions options;
    options.fragment_tolerance_da = getDoubleOption_("mz_tolerance");
    options.ion_mobility_tolerance_vssc = getDoubleOption_("im_tolerance");
    options.rt_window_seconds = getDoubleOption_("rt_window");
    options.peak_gap_seconds = getDoubleOption_("peak_gap");
    const Int min_transitions = getIntOption_("min_transitions");
    if (min_transitions < 1 || !std::isfinite(options.fragment_tolerance_da) || options.fragment_tolerance_da <= 0.0 ||
        !std::isfinite(options.ion_mobility_tolerance_vssc) || options.ion_mobility_tolerance_vssc <= 0.0 ||
        options.rt_window_seconds < 0.0 || options.peak_gap_seconds <= 0.0)
    {
      writeLogError_("m/z and ion-mobility tolerances and peak gap must be positive; RT window must be non-negative; min_transitions must be at least 1.");
      return ILLEGAL_PARAMETERS;
    }
    options.min_transitions = static_cast<std::size_t>(min_transitions);

    OpenSwath::LightTargetedExperiment targeted;
    TransitionPQPFile().convertPQPToTargetedExperiment(library.c_str(), targeted);
    auto queries = makeQueries_(targeted);
    if (queries.empty())
    {
      writeLogError_("The PQP library contains no usable precursor-fragment transitions.");
      return ILLEGAL_PARAMETERS;
    }

    if (!features_path.empty() && features_path == output_path)
    {
      writeLogError_("The candidate feature table must use a different path from the peak table.");
      return ILLEGAL_PARAMETERS;
    }

    std::ofstream output(output_path);
    if (!output)
    {
      writeLogError_("Could not open output file: " + output_path);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    output << "precursor_id\tprecursor_mz\trt_start\trt_apex\trt_end\tarea_score\ttransition_support\tmatched_scans\n";

    std::ofstream features;
    std::ostream* feature_stream = nullptr;
    if (!features_path.empty())
    {
      features.open(features_path);
      if (!features)
      {
        writeLogError_("Could not open candidate feature output file: " + features_path);
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      features << "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\trt_start\trt_apex\trt_end\t"
                  "peak_width_seconds\tarea\tapex_intensity\ttransition_support\tmatched_scans\t"
                  "has_library_rt\tapex_rt_error_seconds\thas_library_im\n";
      feature_stream = &features;
    }

    OpenDIA::StreamSearch consumer(std::move(queries), options, output, feature_stream);
    FileHandler reader;
    PeakFileOptions read_options;
    read_options.setMSLevels({2});
    reader.setOptions(read_options);
    reader.transformExperiment(input, &consumer, {}, log_type_);
    consumer.finish();
    writeLogInfo_(consumer.ionMobilityFilteringUsed()
                    ? "Ion-mobility filtering: active"
                    : "Ion-mobility filtering: not available; using m/z-only matching");
    output.flush();
    if (feature_stream != nullptr)
    {
      features.flush();
    }
    return output && (feature_stream == nullptr || features) ? EXECUTION_OK : CANNOT_WRITE_OUTPUT_FILE;
  }

private:
  static std::vector<OpenDIA::TransitionQuery> makeQueries_(OpenSwath::LightTargetedExperiment& targeted)
  {
    struct LibraryCoordinates
    {
      double seconds = 0.0;
      bool has_rt = false;
      double precursor_im = 0.0;
      bool has_precursor_im = false;
    };
    std::unordered_map<std::string, LibraryCoordinates> coordinates;
    coordinates.reserve(targeted.getCompounds().size());

    for (const auto& compound : targeted.getCompounds())
    {
      LibraryCoordinates value;
      value.seconds = compound.rt;
      value.has_rt = std::isfinite(value.seconds);
      value.precursor_im = compound.drift_time;
      value.has_precursor_im = std::isfinite(value.precursor_im) && value.precursor_im > 0.0;
      coordinates.emplace(compound.id, value);
    }

    // The lightweight reader skips its global SQL sort for large libraries.
    // Restore the heavy reader's transition order before assigning candidate indices.
    std::stable_sort(targeted.transitions.begin(), targeted.transitions.end(), [&coordinates](const auto& lhs, const auto& rhs)
    {
      if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
      if (lhs.product_mz != rhs.product_mz) return lhs.product_mz < rhs.product_mz;
      const auto lhs_coordinates = coordinates.find(lhs.peptide_ref);
      const auto rhs_coordinates = coordinates.find(rhs.peptide_ref);
      const bool lhs_has_rt = lhs_coordinates != coordinates.end() && lhs_coordinates->second.has_rt;
      const bool rhs_has_rt = rhs_coordinates != coordinates.end() && rhs_coordinates->second.has_rt;
      if (lhs_has_rt != rhs_has_rt) return !lhs_has_rt;
      if (lhs_has_rt && lhs_coordinates->second.seconds != rhs_coordinates->second.seconds)
      {
        return lhs_coordinates->second.seconds < rhs_coordinates->second.seconds;
      }
      if (lhs.transition_name != rhs.transition_name) return lhs.transition_name < rhs.transition_name;
      return lhs.peptide_ref < rhs.peptide_ref;
    });

    std::vector<OpenDIA::TransitionQuery> queries;
    queries.reserve(targeted.getTransitions().size());
    for (const auto& transition : targeted.getTransitions())
    {
      const std::string& id = transition.peptide_ref;
      const auto coordinate_it = coordinates.find(id);
      OpenDIA::TransitionQuery query;
      query.precursor_id = id;
      query.precursor_mz = transition.precursor_mz;
      query.fragment_mz = transition.product_mz;
      if (coordinate_it != coordinates.end())
      {
        query.library_rt = coordinate_it->second.seconds;
        query.has_library_rt = coordinate_it->second.has_rt;
        query.precursor_im = coordinate_it->second.precursor_im;
        query.has_precursor_im = coordinate_it->second.has_precursor_im;
      }
      query.is_decoy = transition.getDecoy();
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
