// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace OpenMS;

/**
@page TOPP_TransitionListEvidenceFilter TransitionListEvidenceFilter

@brief Filters a transition list by quick precursor evidence in raw data.

This tool scans one or more mzML/sqMass input runs with the
OpenMS::TransitionListEvidenceFilter module and writes a filtered transition
list for follow-up inspection or downstream DIA processing.

For multiple input files, the default behavior is a union across runs:
a target precursor is kept if it has evidence in any mzML file. Use
@p aggregation_method all to keep only targets observed in every input run.
If @p split_file_input is set, all input files are interpreted as split SWATH
windows from one run, matching OpenSwathWorkflow's split-input mode.

The evidence filter itself is target-only. If the output is meant for full
target-decoy scoring, use @p decoy_handling keep_all or @p decoy_handling
keep_matching; otherwise decoy transitions are removed from the filtered output.

@experimental This tool is experimental and intended for evaluating prefiltering
strategies before integrating filtered libraries into production workflows.

@ingroup TOPP
*/

/// @cond TOPPCLASSES
class TOPPTransitionListEvidenceFilter :
  public TOPPOpenSwathBase
{
public:
  TOPPTransitionListEvidenceFilter() :
    TOPPOpenSwathBase("TransitionListEvidenceFilter",
                      "Filter transition-list precursors by quick raw-data evidence.",
                      true)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input mzML/sqMass files. By default each file is treated as one run.");
    StringList in_formats = {"mzML", "mzXML", "sqMass"};
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
    setValidFormats_("in", in_formats);

    registerInputFile_("tr", "<file>", "", "Input transition file ('TraML', 'tsv', 'pqp', or 'oswpq').");
    StringList tr_formats = {"traML", "tsv", "pqp"};
    tr_formats.push_back("oswpq");
    setValidFormats_("tr", tr_formats);
    registerStringOption_("tr_type", "<type>", "", "Input transition file type. Default: determined from file extension or content.", false);
    setValidStrings_("tr_type", tr_formats);

    registerOutputFile_("out", "<file>", "", "Filtered output transition file.");
    StringList out_formats = {"tsv", "pqp"};
    out_formats.push_back("oswpq");
    setValidFormats_("out", out_formats);
    registerStringOption_("out_type", "<type>", "", "Output transition file type. Default: determined from file extension.", false);
    setValidStrings_("out_type", out_formats);

    registerStringOption_("aggregation_method", "<any|all>", "any",
                          "How to combine evidence across multiple input runs. 'any' keeps a target if it is supported in at least one run; 'all' requires support in every run.",
                          false, true);
    setValidStrings_("aggregation_method", {"any", "all"});

    registerStringOption_("decoy_handling", "<remove|keep_all|keep_matching>", "remove",
                          "How to handle decoys in the output. 'remove' writes target-only output. 'keep_all' keeps all original decoys. 'keep_matching' keeps decoys whose peptide_ref maps to a kept target after removing decoy_prefix.",
                          false, true);
    setValidStrings_("decoy_handling", {"remove", "keep_all", "keep_matching"});
    registerStringOption_("decoy_prefix", "<prefix>", "DECOY_", "Prefix used by decoy_handling=keep_matching to map decoy peptide_ref values back to target peptide_ref values.", false, true);

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0,
                          "Minimal distance to the upper edge of a SWATH window to still consider a precursor, in Thomson.",
                          false, true);
    setMinFloat_("min_upper_edge_dist", 0.0);
    registerFlag_("pasef", "Data is PASEF data. If omitted, the tool auto-detects ion mobility SWATH windows.");
    registerFlag_("prm", "Data is targeted DIA / PRM-like data with potentially overlapping DIA windows.", true);
    registerFlag_("force", "Override SWATH window gap/overlap sanity checks.", true);
    registerInputFile_("swath_windows_file", "<file>", "",
                       "Optional tab-separated file containing SWATH windows for extraction: lower_offset upper_offset. The first line is a header and will be skipped.",
                       false);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file.", true);

    registerFlag_("split_file_input", "Treat all input files as one split SWATH run, with each file containing one SWATH window.", true);
    registerStringOption_("readOptions", "<normal|cache>", "normal",
                          "Whether to run directly on input data or cache data to disk first. If 'cache', set tempDirectory as needed.",
                          false, true);
    setValidStrings_("readOptions", {"normal", "cache"});
    registerStringOption_("tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory for cached data.", false, true);
    registerFlag_("keep_cached_files", "If set, do not remove cached files created in tempDirectory.", false);
    registerIntOption_("outer_loop_threads", "<number>", -1,
                       "How many threads the evidence filter should use (-1 uses the OpenMP maximum).",
                       false, true);

    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 50.0,
                          "MS1 precursor m/z extraction window full width, in Thomson or ppm.",
                          false, true);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<Th|ppm>", "ppm", "Unit of mz_extraction_window_ms1.", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", {"Th", "ppm"});
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1.0,
                          "MS1 ion mobility extraction window full width. -1 disables MS1 IM filtering.",
                          false, true);

    registerDoubleOption_("mz_extraction_window_ms2", "<double>", 50.0,
                          "MS2 fragment m/z extraction window full width, in Thomson or ppm.",
                          false, true);
    setMinFloat_("mz_extraction_window_ms2", 0.0);
    registerStringOption_("mz_extraction_window_ms2_unit", "<Th|ppm>", "ppm", "Unit of mz_extraction_window_ms2.", false, true);
    setValidStrings_("mz_extraction_window_ms2_unit", {"Th", "ppm"});
    registerDoubleOption_("im_extraction_window_ms2", "<double>", -1.0,
                          "MS2 ion mobility extraction window full width. -1 disables MS2 IM filtering.",
                          false, true);

    registerSubsection_("Algorithm", "Raw-data precursor evidence filtering parameters.");
    registerSubsection_("Library", "Transition TSV/PQP parsing parameters.");
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Algorithm")
    {
      Param defaults = TransitionListEvidenceFilter().getDefaults();
      defaults.remove("enabled");
      return defaults;
    }
    if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    return Param();
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in_files = getStringList_("in");
    if (in_files.empty())
    {
      writeLogError_("Error: No input raw data files provided.");
      return ILLEGAL_PARAMETERS;
    }

    FileHandler file_handler;
    const String tr_file = getStringOption_("tr");
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = file_handler.getType(tr_file);
    }
    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine transition input file type.");
      return PARSE_ERROR;
    }

    const String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));
    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = file_handler.getTypeByFileName(out);
    }
    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine transition output file type.");
      return PARSE_ERROR;
    }

    const String decoy_handling = getStringOption_("decoy_handling");
    Param tsv_reader_param = getParam_().copy("Library:", true);
    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionListForPrefilter_(tr_type, tr_file, tsv_reader_param, decoy_handling);
    OPENMS_LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, "
                    << transition_exp.getCompounds().size() << " compounds with "
                    << transition_exp.getTransitions().size() << " transitions.\n";

    ChromExtractParams ms1_params = makeChromExtractParams_(
      getDoubleOption_("mz_extraction_window_ms1"),
      getStringOption_("mz_extraction_window_ms1_unit") == "ppm",
      getDoubleOption_("im_extraction_window_ms1"));
    ms1_params.min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    ChromExtractParams ms2_params = makeChromExtractParams_(
      getDoubleOption_("mz_extraction_window_ms2"),
      getStringOption_("mz_extraction_window_ms2_unit") == "ppm",
      getDoubleOption_("im_extraction_window_ms2"));
    ms2_params.min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");

    Param filter_params = getParam_().copy("Algorithm:", true);
    const Size min_supported_precursors = static_cast<Size>(filter_params.getValue("min_supported_precursors"));
    filter_params.setValue("enabled", "false");

    const bool split_file_input = getFlag_("split_file_input");
    const String aggregation_method = getStringOption_("aggregation_method");
    const bool require_all_runs = aggregation_method == "all";
    const String readoptions = getStringOption_("readOptions");
    const String tmp_dir = File::absolutePath(getStringOption_("tempDirectory")).ensureLastChar('/');
    const bool keep_cached_files = getFlag_("keep_cached_files");
    const bool force = getFlag_("force");
    const bool sort_swath_maps = getFlag_("sort_swath_maps");
    const bool prm = getFlag_("prm");
    const String swath_windows_file = getStringOption_("swath_windows_file");
    const int threads = static_cast<int>(getIntOption_("outer_loop_threads"));

    std::vector<StringList> run_groups;
    if (split_file_input)
    {
      run_groups.push_back(in_files);
    }
    else
    {
      for (const auto& file : in_files)
      {
        run_groups.push_back(StringList{file});
      }
    }

    std::unordered_map<std::string, Size> supported_run_counts;
    double output_precursor_im_scale = 1.0;
    bool output_precursor_im_scaled_by_charge = false;
    Size run_index = 0;
    for (const auto& run_files : run_groups)
    {
      ++run_index;
      OPENMS_LOG_INFO << "Filtering run " << run_index << "/" << run_groups.size()
                      << ": " << ListUtils::concatenate(run_files, ", ") << "\n";

      String per_run_tmp = tmp_dir;
      std::unique_ptr<File::TempDir> per_run_temp_dir;
      if (readoptions == "cache")
      {
        per_run_temp_dir = std::make_unique<File::TempDir>(tmp_dir, keep_cached_files);
        per_run_tmp = per_run_temp_dir->getPath();
      }

      std::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
      std::vector<OpenSwath::SwathMap> swath_maps;
      std::vector<String> swath_map_sources;
      if (!loadSwathFiles(run_files, exp_meta, swath_maps, swath_map_sources, split_file_input,
                          per_run_tmp, readoptions, swath_windows_file,
                          ms2_params.min_upper_edge_dist, force, sort_swath_maps, prm))
      {
        writeLogError_("Error: Failed to load SWATH files for run " + String(run_index));
        return PARSE_ERROR;
      }

      bool run_pasef = getFlag_("pasef");
      if (!run_pasef)
      {
        run_pasef = std::any_of(swath_maps.begin(), swath_maps.end(),
                                [](const OpenSwath::SwathMap& map)
                                {
                                  return !map.ms1 && map.imLower >= 0.0 && map.imUpper >= 0.0;
                                });
        if (run_pasef)
        {
          OPENMS_LOG_INFO << "Auto-detected ion mobility/PASEF windows for run "
                          << run_index << ".\n";
        }
      }

      TransitionListEvidenceFilter filter;
      filter.setParameters(filter_params);
      filter.setLogType(log_type_);
      TransitionListEvidenceFilter::Result result = filter.filter(
        swath_maps, transition_exp, ms1_params, ms2_params, run_pasef, threads);
      if (result.precursor_im_scale != 1.0 || result.precursor_im_scaled_by_charge)
      {
        if (output_precursor_im_scale == 1.0 && !output_precursor_im_scaled_by_charge)
        {
          output_precursor_im_scale = result.precursor_im_scale;
          output_precursor_im_scaled_by_charge = result.precursor_im_scaled_by_charge;
        }
        else if (std::fabs(output_precursor_im_scale - result.precursor_im_scale) > 1e-9 ||
                 output_precursor_im_scaled_by_charge != result.precursor_im_scaled_by_charge)
        {
          OPENMS_LOG_WARN << "TransitionListEvidenceFilter observed different precursor IM transforms across runs. "
                          << "Using the first detected scale factor " << output_precursor_im_scale
                          << (output_precursor_im_scaled_by_charge ? " with precursor charge multiplication" : "")
                          << " for the aggregated output.\n";
        }
      }

      std::unordered_set<std::string> run_supported;
      for (const auto& compound : result.filtered_targets.getCompounds())
      {
        run_supported.insert(compound.id);
      }
      for (const auto& compound_id : run_supported)
      {
        ++supported_run_counts[compound_id];
      }
    }

    std::unordered_set<std::string> selected_targets;
    for (const auto& item : supported_run_counts)
    {
      if ((!require_all_runs && item.second > 0) ||
          (require_all_runs && item.second == run_groups.size()))
      {
        selected_targets.insert(item.first);
      }
    }

    if (selected_targets.size() < min_supported_precursors)
    {
      writeLogError_("Error: retained only " + String(selected_targets.size()) +
                     " target precursors after aggregation, fewer than Algorithm:min_supported_precursors=" +
                     String(min_supported_precursors) + ".");
      return INCOMPATIBLE_INPUT_DATA;
    }

    OpenSwath::LightTargetedExperiment filtered_exp = buildOutputExperiment_(
      transition_exp, selected_targets, decoy_handling, getStringOption_("decoy_prefix"),
      output_precursor_im_scale, output_precursor_im_scaled_by_charge);

    writeTransitionList_(out, out_type, filtered_exp);
    OPENMS_LOG_INFO << "Wrote filtered transition list with " << filtered_exp.getProteins().size()
                    << " proteins, " << filtered_exp.getCompounds().size()
                    << " compounds, and " << filtered_exp.getTransitions().size()
                    << " transitions to " << out << ".\n";
    return EXECUTION_OK;
  }

private:
  OpenSwath::LightTargetedExperiment loadTransitionListForPrefilter_(const FileTypes::Type& tr_type,
                                                                     const String& tr_file,
                                                                     const Param& tsv_reader_param,
                                                                     const String& decoy_handling)
  {
    if (tr_type == FileTypes::PQP && decoy_handling == "keep_matching")
    {
      OPENMS_LOG_INFO << "Loading PQP with legacy TraML IDs to preserve target-decoy precursor ID matching.\n";
      OpenSwath::LightTargetedExperiment transition_exp;
      TransitionPQPFile pqp_reader;
      pqp_reader.setLogType(log_type_);
      pqp_reader.convertPQPToTargetedExperiment(tr_file.c_str(), transition_exp, true);
      return transition_exp;
    }
    return loadTransitionList(tr_type, tr_file, tsv_reader_param);
  }

  static ChromExtractParams makeChromExtractParams_(double mz_window, bool ppm, double im_window)
  {
    ChromExtractParams params;
    params.min_upper_edge_dist = 0.0;
    params.mz_extraction_window = mz_window;
    params.im_extraction_window = im_window;
    params.ppm = ppm;
    params.extraction_function = "tophat";
    params.rt_extraction_window = -1.0;
    params.extra_rt_extract = 0.0;
    return params;
  }

  static bool hasDecoyPrefix_(const std::string& id, const String& decoy_prefix)
  {
    const std::string prefix = decoy_prefix;
    if (!prefix.empty() && id.find(prefix) == 0)
    {
      return true;
    }
    return id.find("DECOY") == 0 || id.find("Decoy") == 0 || id.find("decoy") == 0;
  }

  static bool mapsToSelectedTarget_(const std::string& decoy_ref,
                                    const std::unordered_set<std::string>& selected_targets,
                                    const String& decoy_prefix)
  {
    const std::string prefix = decoy_prefix;
    if (!prefix.empty() && decoy_ref.find(prefix) == 0)
    {
      return selected_targets.find(decoy_ref.substr(prefix.size())) != selected_targets.end();
    }
    return false;
  }

  static OpenSwath::LightTargetedExperiment buildOutputExperiment_(
    const OpenSwath::LightTargetedExperiment& transition_exp,
    const std::unordered_set<std::string>& selected_targets,
    const String& decoy_handling,
    const String& decoy_prefix,
    double precursor_im_scale,
    bool precursor_im_scaled_by_charge)
  {
    OpenSwath::LightTargetedExperiment filtered_exp;
    std::unordered_map<std::string, int> charge_by_compound;
    charge_by_compound.reserve(transition_exp.getCompounds().size());
    for (const auto& compound : transition_exp.getCompounds())
    {
      charge_by_compound[compound.id] = compound.charge;
    }

    std::unordered_set<std::string> kept_compounds;
    for (const auto& transition : transition_exp.getTransitions())
    {
      const bool decoy = transition.getDecoy() || hasDecoyPrefix_(transition.getPeptideRef(), decoy_prefix);
      bool keep = false;
      if (!decoy)
      {
        keep = selected_targets.find(transition.getPeptideRef()) != selected_targets.end();
      }
      else if (decoy_handling == "keep_all")
      {
        keep = true;
      }
      else if (decoy_handling == "keep_matching")
      {
        keep = mapsToSelectedTarget_(transition.getPeptideRef(), selected_targets, decoy_prefix);
      }

      if (keep)
      {
        OpenSwath::LightTransition transition_copy = transition;
        double precursor_im_factor = precursor_im_scale;
        const auto charge_it = charge_by_compound.find(transition.getPeptideRef());
        if (precursor_im_scaled_by_charge && charge_it != charge_by_compound.end() && charge_it->second > 0)
        {
          precursor_im_factor *= charge_it->second;
        }
        if (precursor_im_factor != 1.0 && transition_copy.precursor_im >= 0.0)
        {
          transition_copy.precursor_im *= precursor_im_factor;
        }
        filtered_exp.transitions.push_back(std::move(transition_copy));
        kept_compounds.insert(transition.getPeptideRef());
      }
    }

    std::unordered_set<std::string> kept_proteins;
    for (const auto& compound : transition_exp.getCompounds())
    {
      if (kept_compounds.find(compound.id) != kept_compounds.end())
      {
        OpenSwath::LightCompound compound_copy = compound;
        double precursor_im_factor = precursor_im_scale;
        if (precursor_im_scaled_by_charge && compound_copy.charge > 0)
        {
          precursor_im_factor *= compound_copy.charge;
        }
        if (precursor_im_factor != 1.0 && compound_copy.drift_time >= 0.0)
        {
          compound_copy.drift_time *= precursor_im_factor;
        }
        filtered_exp.compounds.push_back(std::move(compound_copy));
        for (const auto& protein_ref : compound.protein_refs)
        {
          kept_proteins.insert(protein_ref);
        }
      }
    }

    for (const auto& protein : transition_exp.getProteins())
    {
      if (kept_proteins.find(protein.id) != kept_proteins.end())
      {
        filtered_exp.proteins.push_back(protein);
      }
    }

    return filtered_exp;
  }

  void writeTransitionList_(const String& out,
                            FileTypes::Type out_type,
                            const OpenSwath::LightTargetedExperiment& filtered_exp) const
  {
    if (out_type == FileTypes::TSV)
    {
      TransitionTSVFile writer;
      writer.setLogType(log_type_);
      writer.convertLightTargetedExperimentToTSV(out.c_str(), filtered_exp);
    }
    else if (out_type == FileTypes::PQP)
    {
      TransitionPQPFile writer;
      writer.setLogType(log_type_);
      writer.convertLightTargetedExperimentToPQP(out.c_str(), filtered_exp);
    }
    else if (out_type == FileTypes::OSWPQ)
    {
      TransitionParquetFile writer;
      writer.convertLightTargetedExperimentToParquet(out, filtered_exp);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Output type must be tsv, pqp, or oswpq.");
    }
  }
};

int main(int argc, const char** argv)
{
  TOPPTransitionListEvidenceFilter tool;
  return tool.main(argc, argv);
}

/// @endcond
