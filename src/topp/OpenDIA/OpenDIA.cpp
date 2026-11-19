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

#include "OpenDIALibrary.h"

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
  @brief Stream DIA spectra and report the elution peaks of every library precursor.

  Reads a prepared PQP transition library and streams MS2 spectra through OpenMS file handlers. For each precursor,
  every scan that isolated it inside its RT window around the (calibrated) library RT forms native-scan fragment
  traces (the most intense match of each fragment, zero where nothing matched). The traces are peak-picked
  (consensus maxima split at deep valleys, ranked by fragment co-elution) and up to -max_peaks peaks are reported.
  When both the library and raw spectra provide inverse reduced ion mobility, per-peak mobility rejects
  out-of-mobility fragment matches.

  An RT window is required: -calibration (OpenDIACalibrate, recommended) sets it together with the fragment ppm
  tolerance, the mobility tolerance and the boundary cap; without a calibration use -rt_window.

  This is the extraction and candidate-generation stage, not a complete identification workflow. The optional
  feature table keeps candidate-level target/decoy labels and numeric evidence for rescoring. Input spectra are
  expected to contain precursor isolation metadata.
*/
/// @cond TOPPCLASSES
class TOPPOpenDIA final : public TOPPBase
{
public:
  TOPPOpenDIA() : TOPPBase("OpenDIA", "Stream DIA data and pick targeted fragment elution peaks.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "DIA raw data or mzML/mzXML file");
    registerInputFile_("lib", "<file>", "", "Prepared PQP transition library");
    setValidFormats_("lib", {"pqp"});
    registerOutputFile_("out", "<file>", "", "Peak table");
    setValidFormats_("out", {"tsv"});
    registerOutputFile_("features_out", "<file>", "", "Optional candidate-level numeric feature table for rescoring", false);
    setValidFormats_("features_out", {"tsv"});
    registerDoubleOption_("mz_tolerance", "<Da>", 0.02, "Absolute fragment m/z matching tolerance in Da", false);
    registerDoubleOption_("im_tolerance", "<1/K0>", 0.03, "Automatic inverse reduced ion-mobility matching tolerance", false);
    registerDoubleOption_("rt_window", "<seconds>", 0.0, "Half width of the RT window around library RT; required without -calibration", false);
    registerInputFile_("calibration", "<file>", "", "OpenDIACalibrate JSON. Before extraction, library RT, fragment m/z and precursor 1/K0 are mapped to run coordinates, and the RT window, fragment ppm and 1/K0 tolerances and the boundary cap are set from the calibration (replacing -rt_window, -mz_tolerance and -im_tolerance)", false);
    setValidFormats_("calibration", {"json"});
    registerIntOption_("min_transitions", "<n>", 1, "Minimum fragments with signal at or next to a peak apex", false);
    registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0, e.g. unitless mzML 'Ion Mobility' arrays");
    registerDoubleOption_("valley_fraction", "<fraction>", 0.5, "Two maxima stay separate peaks only if the valley between them is at most this fraction of the lower one", false);
    setMinFloat_("valley_fraction", 0.001);
    setMaxFloat_("valley_fraction", 1.0);
    registerIntOption_("max_peaks", "<n>", 5, "Peaks reported per precursor, best fragment co-elution first", false);
    setMinInt_("max_peaks", 1);
    registerDoubleOption_("boundary_half_widths", "<k>", -1.0, "Cap each peak boundary at k of that side's own half widths at half "
                          "maximum. -1 uses the value estimated by OpenDIACalibrate (peak_boundary.half_widths in -calibration; no cap if "
                          "absent); 0 disables the cap", false);
    setMinFloat_("boundary_half_widths", -1.0);
    registerIntOption_("min_profile_scans", "<n>", 2, "A fragment enters the co-elution score only with signal on at least n "
                       "scans of the correlation window (1: two fragments sharing one nonzero scan correlate perfectly)", false);
    setMinInt_("min_profile_scans", 1);
    registerOutputFile_("traces_out", "<file>", "", "Inspection: native-scan fragment traces (the peak picker's input matrix) of the "
                        "-traces_precursors precursors as TSV. Does not change -out.", false);
    setValidFormats_("traces_out", {"tsv"});
    registerInputFile_("traces_precursors", "<file>", "", "Text file of PQP precursor IDs (one per line, '#' comments) for -traces_out", false);

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
    options.ion_mobility_array_as_vssc = getFlag_("im_array_as_vssc");
    options.picker.valley_fraction = getDoubleOption_("valley_fraction");
    options.picker.max_peaks = static_cast<std::size_t>(std::max(1, getIntOption_("max_peaks")));
    options.picker.min_profile_scans = static_cast<std::size_t>(std::max(1, getIntOption_("min_profile_scans")));
    const double boundary_half_widths = getDoubleOption_("boundary_half_widths");
    options.picker.boundary_half_widths = std::max(0.0, boundary_half_widths);
    const Int min_transitions = getIntOption_("min_transitions");
    if (min_transitions < 1 || !std::isfinite(options.fragment_tolerance_da) || options.fragment_tolerance_da <= 0.0 ||
        !std::isfinite(options.ion_mobility_tolerance_vssc) || options.ion_mobility_tolerance_vssc <= 0.0 ||
        options.rt_window_seconds < 0.0)
    {
      writeLogError_("m/z and ion-mobility tolerances must be positive; RT window must be non-negative; min_transitions must be at least 1.");
      return ILLEGAL_PARAMETERS;
    }
    options.min_transitions = static_cast<std::size_t>(min_transitions);

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
      if (boundary_half_widths < 0.0)
      {
        options.picker.boundary_half_widths = calibration->boundary_half_widths;
      }
      if (calibration->has_im)
      {
        options.ion_mobility_tolerance_vssc = calibration->im_tolerance;
      }
    }

    OpenDIA::Library targeted = OpenDIA::readLibrary(library);
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
      writeLogInfo_("Calibrated extraction: rt_window_s=" + std::to_string(options.rt_window_seconds) +
                    " ms2_offset_ppm=" + std::to_string(calibration->ms2_offset_ppm) +
                    " ms2_tolerance_ppm=" + std::to_string(options.fragment_tolerance_ppm) +
                    " boundary_half_widths=" + std::to_string(options.picker.boundary_half_widths) +
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
    if (!(options.rt_window_seconds > 0.0))
    {
      writeLogError_("An RT window is required (-calibration or -rt_window).");
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
      features << OpenDIA::StreamSearch::featureHeader() << '\n';
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
      OpenDIA::Library subset;
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
      // The .xipm writer (OpenMS format I/O) takes OpenSWATH types; they are built here, at the export edge only.
      OpenSwath::LightTargetedExperiment xipm_library;
      for (const auto& c : subset.compounds)
      {
        OpenSwath::LightCompound x;
        x.id = c.id;
        x.sequence = c.sequence;
        x.charge = c.charge;
        x.rt = c.rt;
        x.drift_time = c.drift_time;
        xipm_library.compounds.push_back(std::move(x));
      }
      for (const auto& t : subset.transitions)
      {
        OpenSwath::LightTransition x;
        x.transition_name = t.transition_name;
        x.peptide_ref = t.peptide_ref;
        x.library_intensity = t.library_intensity;
        x.product_mz = t.product_mz;
        x.precursor_mz = t.precursor_mz;
        x.precursor_im = t.precursor_im;
        x.fragment_charge = t.fragment_charge;
        x.fragment_nr = t.fragment_nr;
        x.setFragmentType(t.fragment_type);
        x.setDecoy(t.decoy);
        xipm_library.transitions.push_back(std::move(x));
      }
      xipm_writer = std::make_unique<XIPMParquetConsumer>(xipm_path, xipm_library);
      const std::string source_file = File::basename(input);
      XIPMParquetConsumer* writer = xipm_writer.get();
      peak_maps = std::make_unique<OpenDIA::PeakMapCollector>(
        std::move(targets), map_options,
        [writer, source_file](const OpenDIA::PeakMap& map, int ms_level)
        {
          PeakMapExtractor::ExtractedPeakMap x;
          x.native_id = map.native_id;
          x.target_mz = map.target_mz;
          x.target_rt = map.target_rt;
          x.target_ion_mobility = map.target_ion_mobility;
          x.rt_start = map.rt_start;
          x.rt_end = map.rt_end;
          x.mz = map.mz;
          x.rt = map.rt;
          x.ion_mobility = map.ion_mobility;
          x.intensity = map.intensity;
          writer->consumePeakMap(x, 0, source_file, ms_level);
        });
      writeLogInfo_("Peak-map export: " + std::to_string(subset.compounds.size()) + " precursors, rt_half_window_s=" +
                    std::to_string(rt_half_window) + " mz_tolerance_ppm=" + std::to_string(map_options.mz_tolerance_ppm));
    }

    OpenDIA::StreamSearch consumer(std::move(queries), options, output, feature_stream);
    const std::string traces_path = getStringOption_("traces_out");
    std::ofstream traces;
    if (!traces_path.empty())
    {
      const std::string traces_selection = getStringOption_("traces_precursors");
      if (traces_selection.empty())
      {
        writeLogError_("-traces_out needs -traces_precursors.");
        return ILLEGAL_PARAMETERS;
      }
      std::ifstream selection(traces_selection);
      if (!selection)
      {
        writeLogError_("Could not read -traces_precursors: " + traces_selection);
        return INPUT_FILE_NOT_READABLE;
      }
      std::vector<std::string> ids;
      for (std::string line; std::getline(selection, line);)
      {
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;
        const auto last = line.find_last_not_of(" \t\r");
        ids.push_back(line.substr(first, last - first + 1));
      }
      traces.open(traces_path);
      if (!traces)
      {
        writeLogError_("Could not open -traces_out: " + traces_path);
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      consumer.setTraceDump(traces, ids);
    }
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
  static std::vector<OpenDIA::TransitionQuery> makeQueries_(OpenDIA::Library& targeted)
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

    // The lightweight reader skips its global SQL sort for large libraries. Restore the heavy reader's transition
    // order before assigning candidate indices: sort compact keys, then move each transition once.
    struct SortKey
    {
      double precursor_mz = 0.0;
      double product_mz = 0.0;
      double seconds = 0.0;
      bool has_rt = false;
      std::size_t index = 0;
    };
    auto& transitions = targeted.transitions;
    std::vector<SortKey> keys(transitions.size());
    for (std::size_t i = 0; i < transitions.size(); ++i)
    {
      const auto coordinate_it = coordinates.find(transitions[i].peptide_ref);
      const bool has_rt = coordinate_it != coordinates.end() && coordinate_it->second.has_rt;
      keys[i] = {transitions[i].precursor_mz, transitions[i].product_mz, has_rt ? coordinate_it->second.seconds : 0.0, has_rt, i};
    }
    std::stable_sort(keys.begin(), keys.end(), [&transitions](const SortKey& lhs, const SortKey& rhs)
    {
      if (lhs.precursor_mz != rhs.precursor_mz) return lhs.precursor_mz < rhs.precursor_mz;
      if (lhs.product_mz != rhs.product_mz) return lhs.product_mz < rhs.product_mz;
      if (lhs.has_rt != rhs.has_rt) return !lhs.has_rt;
      if (lhs.has_rt && lhs.seconds != rhs.seconds) return lhs.seconds < rhs.seconds;
      const auto& left = transitions[lhs.index];
      const auto& right = transitions[rhs.index];
      if (left.transition_name != right.transition_name) return left.transition_name < right.transition_name;
      return left.peptide_ref < right.peptide_ref;
    });
    std::vector<OpenDIA::LibraryTransition> sorted;
    sorted.reserve(transitions.size());
    for (const SortKey& key : keys)
    {
      sorted.push_back(std::move(transitions[key.index]));
    }
    transitions = std::move(sorted);

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
