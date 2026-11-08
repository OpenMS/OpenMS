// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACalibration.h"
#include "OpenDIAPeakMaps.h"
#include "OpenDIASearch.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace OpenMS;

/**
  @page TOPP_OpenDIA OpenDIA
  @brief Stream DIA spectra and report targeted fragment-group peak candidates.

  This implementation reads a prepared PQP transition library and streams
  MS2 spectra through OpenMS file handlers. By default, it reports gap-segmented
  groups of co-observed library fragments to a TSV file. An opt-in experimental
  mode can split a group after a deep valley and recovery. When both the library
  and raw spectra provide inverse reduced ion mobility, it automatically uses
  per-peak mobility to reject out-of-mobility fragment matches. If either side lacks
  compatible mobility metadata, matching falls back to m/z and isolation windows.

  This is a minimal extraction and candidate-generation path, not a complete
  identification workflow. It does not perform library generation, machine
  learning scoring, FDR estimation, or protein inference. An optional feature
  table can preserve candidate-level target/decoy labels and numeric evidence,
  including measured mobility fit when both sides provide it, for the separate
  experimental OpenDIARerank tool. Input spectra are expected
  to contain precursor isolation metadata.
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
    registerInputFile_("calibration", "<file>", "", "Experimental: OpenDIACalibrate JSON. Before extraction, library RT, fragment m/z and precursor 1/K0 are mapped to run coordinates, and the RT window, fragment ppm and 1/K0 tolerances are set from the calibration (replacing -rt_window, -mz_tolerance and -im_tolerance)", false);
    setValidFormats_("calibration", {"json"});
    registerDoubleOption_("peak_gap", "<seconds>", 30.0, "Maximum gap between matched scans in one peak group", false);
    registerIntOption_("min_transitions", "<n>", 1, "Minimum co-observed fragments required in a peak group", false);
    registerFlag_("local_peak_split", "Split a group after two low scans and a signal recovery; experimental candidate generation");
    registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0, e.g. unitless mzML 'Ion Mobility' arrays");
    registerDoubleOption_("peak_valley_fraction", "<fraction>", 0.1, "Local-peak valley/recovery threshold as a fraction of the previous apex", false);
    setMinFloat_("peak_valley_fraction", 0.001);
    setMaxFloat_("peak_valley_fraction", 0.999);
    registerIntOption_("min_peak_scans", "<n>", 3, "Minimum matched scans in an experimental local-peak candidate", false);
    setMinInt_("min_peak_scans", 2);
    registerIntOption_("max_apexes_per_group", "<n>", 1, "Report up to n separated local apexes per peak group as separate rows sharing the group's boundaries; 1 keeps one global apex", false);
    setMinInt_("max_apexes_per_group", 1);
    registerDoubleOption_("min_apex_separation", "<seconds>", 10.0, "Minimum RT distance between reported apexes of one peak group", false);
    setMinFloat_("min_apex_separation", 0.0);

    // Inspection side channel: raw targeted peak maps in the OpenSWATH .xipm Parquet schema. Off by default;
    // the search output does not depend on it.
    registerOutputFile_("out_xipm", "<file>", "", "Optional: raw MS1 and MS2 m/z, RT, 1/K0 and intensity points of selected precursors "
                        "as an OpenSWATH peak-map Parquet file (.xipm), for inspecting extractions. Does not change -out.", false);
    setValidFormats_("out_xipm", {"xipm"});
    registerInputFile_("xipm_precursors", "<file>", "", "Text file of PQP precursor IDs (one per line, '#' comments) to export with -out_xipm. "
                       "Empty exports every precursor, which can need tens of GB on a full library.", false);
    registerDoubleOption_("xipm_rt_window", "<seconds>", -1.0, "Half width of the exported RT range around the (calibrated) library RT. "
                          "-1 uses the search RT window plus 30 s; 0 exports the whole run.", false);
    setMinFloat_("xipm_rt_window", -1.0);
    registerDoubleOption_("xipm_mz_tolerance", "<ppm>", 20.0, "Half width of every exported m/z window in ppm", false);
    setMinFloat_("xipm_mz_tolerance", 0.0);
    registerDoubleOption_("xipm_im_window", "<1/K0>", 0.0, "Half width of the exported 1/K0 range around the library mobility; 0 keeps every mobility", false);
    setMinFloat_("xipm_im_window", 0.0);
    registerIntOption_("xipm_ms1_isotopes", "<n>", 2, "MS1 isotopes exported after the monoisotope", false);
    setMinInt_("xipm_ms1_isotopes", 0);
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
    options.split_local_peaks = getFlag_("local_peak_split");
    options.ion_mobility_array_as_vssc = getFlag_("im_array_as_vssc");
    options.peak_valley_fraction = getDoubleOption_("peak_valley_fraction");
    options.max_apexes_per_group = static_cast<std::size_t>(std::max(1, getIntOption_("max_apexes_per_group")));
    options.min_apex_separation_seconds = getDoubleOption_("min_apex_separation");
    const Int min_peak_scans = getIntOption_("min_peak_scans");
    const Int min_transitions = getIntOption_("min_transitions");
    if (min_transitions < 1 || !std::isfinite(options.fragment_tolerance_da) || options.fragment_tolerance_da <= 0.0 ||
        !std::isfinite(options.ion_mobility_tolerance_vssc) || options.ion_mobility_tolerance_vssc <= 0.0 ||
        options.rt_window_seconds < 0.0 || options.peak_gap_seconds <= 0.0 || min_peak_scans < 2 ||
        !std::isfinite(options.peak_valley_fraction) || options.peak_valley_fraction <= 0.0 ||
        options.peak_valley_fraction >= 1.0)
    {
      writeLogError_("m/z and ion-mobility tolerances and peak gap must be positive; RT window must be non-negative; min_transitions must be at least 1.");
      return ILLEGAL_PARAMETERS;
    }
    options.min_transitions = static_cast<std::size_t>(min_transitions);
    options.min_peak_scans = static_cast<std::size_t>(min_peak_scans);

    // Optional identification-free calibration (OpenSwath order: calibrate, then extract the whole library).
    const std::string calibration_path = getStringOption_("calibration");
    std::optional<OpenDIA::CalibrationRT> calibration;
    if (!calibration_path.empty())
    {
      try
      {
        calibration = OpenDIA::readCalibrationRT(calibration_path);
      }
      catch (const std::exception& error)
      {
        writeLogError_(std::string("Invalid calibration: ") + error.what());
        return ILLEGAL_PARAMETERS;
      }
      if (!calibration->has_ms2)
      {
        writeLogError_("The calibration has no fragment mass model; refusing a partial calibration.");
        return ILLEGAL_PARAMETERS;
      }
      options.rt_window_seconds = calibration->tolerance_s;
      options.fragment_tolerance_ppm = calibration->ms2_tolerance_ppm;
      if (calibration->has_im)
      {
        options.ion_mobility_tolerance_vssc = calibration->im_tolerance;
      }
    }

    OpenSwath::LightTargetedExperiment targeted;
    TransitionPQPFile().convertPQPToTargetedExperiment(library.c_str(), targeted);
    if (calibration)
    {
      // Library to run coordinates. The RT field then holds the predicted run RT, so the RT window and the
      // reported apex_rt_error_seconds are relative to the calibrated position.
      for (auto& compound : targeted.compounds)
      {
        if (std::isfinite(compound.rt))
        {
          compound.rt = calibration->map(compound.rt);
        }
        if (calibration->has_im && std::isfinite(compound.drift_time) && compound.drift_time > 0.0)
        {
          compound.drift_time = calibration->im_intercept + calibration->im_slope * compound.drift_time;
        }
      }
      for (auto& transition : targeted.transitions)
      {
        transition.product_mz *= 1.0 + calibration->ms2_offset_ppm * 1e-6;
      }
      writeLogInfo_("Experimental calibrated extraction: rt_window_s=" + std::to_string(options.rt_window_seconds) +
                    " ms2_offset_ppm=" + std::to_string(calibration->ms2_offset_ppm) +
                    " ms2_tolerance_ppm=" + std::to_string(options.fragment_tolerance_ppm) +
                    (calibration->has_im ? " im=" + std::to_string(calibration->im_intercept) + "+" + std::to_string(calibration->im_slope) +
                                             "x im_tolerance=" + std::to_string(options.ion_mobility_tolerance_vssc)
                                         : std::string(" im=unavailable (library 1/K0 and -im_tolerance kept)")));
    }
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
                  "has_library_rt\tapex_rt_error_seconds\thas_library_im\t"
                  "has_im_match_evidence\tmean_abs_im_error_vssc\n";
      feature_stream = &features;
    }

    // Optional peak-map export, built before the search consumes the queries.
    const std::string xipm_path = getStringOption_("out_xipm");
    std::unique_ptr<XIPMParquetConsumer> xipm_writer;
    std::unique_ptr<OpenDIA::PeakMapCollector> peak_maps;
    if (!xipm_path.empty())
    {
      if (xipm_path == output_path || xipm_path == features_path)
      {
        writeLogError_("-out_xipm must use a different path from the other outputs.");
        return ILLEGAL_PARAMETERS;
      }
      std::optional<std::unordered_set<std::string>> selected;
      const std::string selection_path = getStringOption_("xipm_precursors");
      if (!selection_path.empty())
      {
        std::ifstream selection(selection_path);
        if (!selection)
        {
          writeLogError_("Could not read -xipm_precursors: " + selection_path);
          return INPUT_FILE_NOT_READABLE;
        }
        selected.emplace();
        for (std::string line; std::getline(selection, line);)
        {
          const auto first = line.find_first_not_of(" \t\r");
          if (first == std::string::npos || line[first] == '#') continue;
          const auto last = line.find_last_not_of(" \t\r");
          selected->insert(line.substr(first, last - first + 1));
        }
      }
      double rt_half_window = getDoubleOption_("xipm_rt_window");
      if (rt_half_window < 0.0)
      {
        rt_half_window = options.rt_window_seconds > 0.0 ? options.rt_window_seconds + 30.0 : 0.0;
      }
      OpenDIA::PeakMapOptions map_options;
      map_options.mz_tolerance_ppm = getDoubleOption_("xipm_mz_tolerance");
      map_options.rt_half_window_s = rt_half_window;
      map_options.im_half_window = getDoubleOption_("xipm_im_window");
      map_options.ms1_isotopes = static_cast<std::size_t>(getIntOption_("xipm_ms1_isotopes"));
      map_options.ion_mobility_array_as_vssc = options.ion_mobility_array_as_vssc;
      if (!(map_options.mz_tolerance_ppm > 0.0))
      {
        writeLogError_("-xipm_mz_tolerance must be positive.");
        return ILLEGAL_PARAMETERS;
      }

      // Library subset for the writer's row metadata; ids and (calibrated) coordinates are copied unchanged.
      OpenSwath::LightTargetedExperiment subset;
      std::unordered_map<std::string, std::size_t> target_index;
      std::vector<OpenDIA::PeakMapTarget> targets;
      for (const auto& compound : targeted.getCompounds())
      {
        if (selected && !selected->contains(compound.id)) continue;
        OpenDIA::PeakMapTarget target;
        target.precursor_id = compound.id;
        target.charge = compound.charge;
        target.center_rt = compound.rt;
        target.has_rt = std::isfinite(compound.rt);
        target.precursor_im = compound.drift_time;
        target.has_im = std::isfinite(compound.drift_time) && compound.drift_time > 0.0;
        target_index.emplace(compound.id, targets.size());
        targets.push_back(std::move(target));
        subset.compounds.push_back(compound);
      }
      for (const auto& transition : targeted.getTransitions())
      {
        const auto it = target_index.find(transition.peptide_ref);
        if (it == target_index.end()) continue;
        targets[it->second].precursor_mz = transition.precursor_mz;
        targets[it->second].fragments.emplace_back(transition.transition_name, transition.product_mz);
        subset.transitions.push_back(transition);
      }
      if (selected && targets.size() != selected->size())
      {
        writeLogWarn_("-xipm_precursors: " + std::to_string(selected->size() - targets.size()) + " of " +
                      std::to_string(selected->size()) + " ids are not library precursors and are skipped.");
      }
      if (targets.empty())
      {
        writeLogError_("-out_xipm: no library precursor selected.");
        return ILLEGAL_PARAMETERS;
      }
      xipm_writer = std::make_unique<XIPMParquetConsumer>(xipm_path, subset);
      const std::string source_file = File::basename(input);
      XIPMParquetConsumer* writer = xipm_writer.get();
      peak_maps = std::make_unique<OpenDIA::PeakMapCollector>(
        std::move(targets), map_options,
        [writer, source_file](const PeakMapExtractor::ExtractedPeakMap& map, int ms_level)
        { writer->consumePeakMap(map, 0, source_file, ms_level); });
      writeLogInfo_("Peak-map export: " + std::to_string(subset.compounds.size()) + " precursors, rt_half_window_s=" +
                    std::to_string(rt_half_window) + " mz_tolerance_ppm=" + std::to_string(map_options.mz_tolerance_ppm));
    }

    OpenDIA::StreamSearch consumer(std::move(queries), options, output, feature_stream);
    FileHandler reader;
    PeakFileOptions read_options;
    // MS1 is read only for the peak-map export; StreamSearch ignores MS1 spectra.
    read_options.setMSLevels(peak_maps ? std::vector<Int>{1, 2} : std::vector<Int>{2});
    reader.setOptions(read_options);
    if (peak_maps)
    {
      OpenDIA::TeeConsumer tee(consumer, *peak_maps);
      reader.transformExperiment(input, &tee, {}, log_type_);
    }
    else
    {
      reader.transformExperiment(input, &consumer, {}, log_type_);
    }
    consumer.finish();
    if (peak_maps)
    {
      peak_maps->finish();
      xipm_writer->finalize();
      writeLogInfo_("Peak-map export: " + std::to_string(peak_maps->peakMapsWritten()) + " peak maps, " +
                    std::to_string(peak_maps->pointsWritten()) + " points, " +
                    std::to_string(peak_maps->outOfOrderSpectra()) + " out-of-order spectra -> " + xipm_path);
    }
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
