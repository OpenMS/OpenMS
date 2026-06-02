// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>

using namespace OpenMS;

/**
@page TOPP_OpenSwathPeakMapExtractor OpenSwathPeakMapExtractor

@brief Extract targeted raw mz/RT/IM peak maps from DIA / diaPASEF data.

This tool reuses the OpenSWATH targeting and calibration logic, but instead of
integrating chromatograms it retains the full raw point cloud inside each
targeted extraction window. The output is written as Parquet `.xipm` files
containing one row per extracted precursor or transition with compressed
parallel arrays for m/z, RT, ion mobility, and intensity.

The tool supports:
 - single-file and multi-run mzML / sqMass / Bruker TDF style OpenSWATH input
 - transition libraries in TraML / TSV / PQP / OSWPQ
 - optional RT normalization and auto-iRT calibration
 - user-provided or auto-estimated extraction windows
 - aggregated multi-run output or one `.xipm` per run

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathPeakMapExtractor.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathPeakMapExtractor.html
*/

/// @cond TOPPCLASSES
class TOPPOpenSwathPeakMapExtractor :
  public TOPPOpenSwathBase
{
public:
  TOPPOpenSwathPeakMapExtractor() :
    TOPPOpenSwathBase("OpenSwathPeakMapExtractor", "Extract targeted mz/RT/IM peak maps from DIA or diaPASEF data", true)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    StringList in_formats = {"mzML", "mzXML", "sqMass"};
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
#ifdef WITH_THERMO_RAW
    in_formats.push_back("raw");
#endif
    setValidFormats_("in", in_formats);

    registerInputFile_("tr", "<file>", "", "transition file ('TraML','tsv','pqp','oswpq')");
    StringList tr_formats = {"traML", "tsv", "pqp", "oswpq"};
    setValidFormats_("tr", tr_formats);
    registerStringOption_("tr_type", "<type>", "", "Input transition file type -- default: determined from file extension", false);
    setValidStrings_("tr_type", tr_formats);

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab-separated file containing the SWATH windows for extraction: lower_offset upper_offset. The first line is treated as a header.", false);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerOutputFile_("out", "<file>", "", "Output .xipm parquet file containing all runs. Required unless -separate_runs is set.", false);
    setValidFormats_("out", ListUtils::create<String>("xipm"));
    registerOutputDir_("out_dir", "<dir>", "", "Output directory for per-run .xipm files when -separate_runs is set.", false);
    registerFlag_("separate_runs", "Write one .xipm file per input run instead of aggregating all runs into -out", false);

    registerStringOption_("enable_ms1", "<true|false>", "true", "Extract precursor peak maps from the MS1 map if present", false, true);
    setValidStrings_("enable_ms1", ListUtils::create<String>("true,false"));
    registerIntOption_("ms1_isotopes", "<number>", 0, "The number of additional MS1 isotopes used for extraction", false, true);
    setMinInt_("ms1_isotopes", 0);

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the upper edge of a SWATH window to still consider a precursor, in Thomson", false, true);
    registerFlag_("pasef", "Data is PASEF data");
    registerStringOption_("matching_window_only", "<name>", "false", "Assume the input data is targeted / PRM-like data with potentially overlapping DIA windows. Only extract each assay from the best matching DIA window.", false, true);
    setValidStrings_("matching_window_only", ListUtils::create<String>("true,false"));

    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range). This is the full window size.", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Extract additional RT beyond the primary window for inspection.", false, true);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerDoubleOption_("ion_mobility_window", "<double>", -1.0, "Extraction window in ion mobility dimension. This is the full window size. -1 extracts the full IM range.", false);
    registerDoubleOption_("mz_extraction_window", "<double>", 50.0, "Extraction window in Thomson or ppm (see mz_extraction_window_unit)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    registerStringOption_("mz_extraction_window_unit", "<name>", "ppm", "Unit for mz extraction", false, true);
    setValidStrings_("mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));

    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 50.0, "Extraction window used in MS1 in Thomson or ppm (see mz_extraction_window_ms1_unit)", false);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<name>", "ppm", "Unit of the MS1 m/z extraction window", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", ListUtils::create<String>("ppm,Th"));
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1.0, "Extraction window in ion mobility dimension for MS1. -1 extracts the full IM range.", false);
    registerStringOption_("use_ms1_ion_mobility", "<name>", "true", "Also apply ion mobility extraction to MS1 peak-map extraction", false, true);
    setValidStrings_("use_ms1_ion_mobility", ListUtils::create<String>("true,false"));

    registerDoubleOption_("irt_mz_extraction_window", "<double>", 50.0, "Extraction window used for iRT and m/z correction in Thomson or ppm (see irt_mz_extraction_window_unit)", false, true);
    setMinFloat_("irt_mz_extraction_window", 0.0);
    registerStringOption_("irt_mz_extraction_window_unit", "<name>", "ppm", "Unit for iRT mz extraction", false, true);
    setValidStrings_("irt_mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));
    registerDoubleOption_("irt_im_extraction_window", "<double>", -1.0, "Ion mobility extraction window used for iRT calibration. -1 disables IM calibration.", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH window (alternatively: all SWATHs are in separate files)", true);
    registerStringOption_("readOptions", "<name>", "normal", "Whether to run directly on the input data, cache data to disk first, or load working sets into memory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory,workingInMemory"));
    registerStringOption_("tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory used for cached files", false, true);
    registerFlag_("keep_cached_files", "Do not remove cached files created in tempDirectory", false);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat"));
    registerIntOption_("batchSize", "<number>", 1000, "Compound batch size per SWATH window. Set 0 to process all compounds for a window in one batch.", false, true);
    setMinInt_("batchSize", 0);

    registerSubsection_("Library", "Library parameters section");
    registerSubsection_("Calibration", "Parameters for calibrant iRT peptides for RT normalization and mass / ion mobility correction.");
    registerSubsection_("Calibration:RTNormalization", "Parameters for the RTNormalization for iRT peptides.");
    registerSubsection_("Calibration:MassIMCorrection", "Parameters for the m/z and ion mobility calibration.");
    registerSubsection_("MRMMapping", "Parameters for mapping chromatograms to transitions during iRT calibration.");

    registerTOPPSubsection_("Debugging", "Debugging");
    registerOutputFile_("Debugging:irt_mzml", "<file>", "", "Chromatogram mzML containing the iRT peptides", false);
    setValidFormats_("Debugging:irt_mzml", ListUtils::create<String>("mzML"));
    registerOutputFile_("Debugging:irt_trafo", "<file>", "", "Transformation file for RT transform", false);
    setValidFormats_("Debugging:irt_trafo", ListUtils::create<String>("trafoXML"));
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    if (name == "Calibration")
    {
      CalibrationWorkflow cal_wf;
      Param p = cal_wf.getDefaults();
      p.setValue("tr_irt_priority_sampling", "", "Optional custom transition file (TSV format only) containing additional priority peptides for iRT sampling.");
      p.setValue("rt_norm", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, no automatic RT calibration is performed.");
      return p;
    }
    if (name == "Calibration:RTNormalization")
    {
      Param p;
      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points.");
      p.setValidStrings("alignmentMethod", {"linear","interpolated","lowess","b_spline"});
      p.setValue("lowess:auto_span", "true", "If true, or if 'span' is 0, automatically select LOWESS span by cross-validation.");
      p.setValidStrings("lowess:auto_span", {"true","false"});
      p.setValue("lowess:span", 0.05, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("lowess:auto_span_min", 0.15, "Lower bound for auto-selected span.");
      p.setMinFloat("lowess:auto_span_min", 0.001);
      p.setValue("lowess:auto_span_max", 0.80, "Upper bound for auto-selected span.");
      p.setMaxFloat("lowess:auto_span_max", 0.99);
      p.setValue("lowess:auto_span_grid", "0.005,0.01,0.05,0.15,0.25,0.30,0.50,0.70,0.90", "Optional explicit grid of span candidates in (0,1].");
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);
      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods.");
      p.setValidStrings("useIterativeChauvenet", {"true","false"});
      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");
      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization.");
      p.setValidStrings("estimateBestPeptides", {"true","false"});
      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored.");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation.");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage.");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to count as covered.");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered.");
      return p;
    }
    if (name == "Calibration:MassIMCorrection")
    {
      return SwathMapMassCorrection().getDefaults();
    }
    if (name == "MRMMapping")
    {
      Param p;
      p.setValue("precursor_tolerance", 0.9, "Precursor tolerance when mapping (in Th)");
      p.setValue("product_tolerance", 1.2, "Product tolerance when mapping (in Th)");
      p.setValue("irt_precursor_tolerance", 1.5, "Precursor tolerance when mapping iRT transitions (in Th)");
      p.setValue("irt_product_tolerance", 1.5, "Product tolerance when mapping iRT transitions (in Th)");
      p.setValue("map_multiple_assays", "false", "Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.");
      p.setValidStrings("map_multiple_assays", {"true","false"});
      p.setValue("error_on_unmapped", "false", "Treat remaining, unmapped chromatograms as an error");
      p.setValidStrings("error_on_unmapped", {"true","false"});
      return p;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
  }

  static Param makeCalibrationFeatureFinderDefaults_(const bool use_ms1_ion_mobility)
  {
    Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
    feature_finder_param.setValue("stop_report_after_feature", 5);
    feature_finder_param.setValue("rt_normalization_factor", 100.0);
    feature_finder_param.setValue("Scores:use_ms1_mi", "true");
    feature_finder_param.setValue("Scores:use_mi_score", "true");
    feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
    feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false");
    feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
    feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
    feature_finder_param.setValue("TransitionGroupPicker:compute_peak_shape_metrics", "false");
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false");
    feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
    feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
    feature_finder_param.setValue("Scores:use_ion_mobility_scores", "false");
    feature_finder_param.setValue("use_ms1_ion_mobility", use_ms1_ion_mobility ? "true" : "false");
    return feature_finder_param;
  }

  static void prepareExtractionCoordinates_(std::vector<PeakMapExtractor::ExtractionCoordinates>& coordinates,
                                            const OpenSwath::LightTargetedExperiment& transition_exp_used,
                                            const TransformationDescription& trafo_inverse,
                                            const ChromExtractParams& cp,
                                            const bool ms1,
                                            const int ms1_isotopes)
  {
    std::vector<OpenSwath::ChromatogramPtr> tmp_chroms;
    if (cp.rt_extraction_window < 0)
    {
      ChromatogramExtractor::prepare_coordinates(tmp_chroms, coordinates, transition_exp_used, cp.rt_extraction_window, ms1, ms1_isotopes);
    }
    else
    {
      ChromatogramExtractor::prepare_coordinates(tmp_chroms, coordinates, transition_exp_used, 0.0, ms1, ms1_isotopes);
      for (auto& coord : coordinates)
      {
        coord.rt_start = trafo_inverse.apply(coord.rt_start) - (cp.rt_extraction_window + cp.extra_rt_extract) / 2.0;
        coord.rt_end = trafo_inverse.apply(coord.rt_end) + (cp.rt_extraction_window + cp.extra_rt_extract) / 2.0;
      }
    }
  }

  static void copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
                                    const std::vector<OpenSwath::LightTransition>& all_transitions,
                                    std::vector<OpenSwath::LightTransition>& output)
  {
    std::unordered_set<std::string> selected_compounds;
    selected_compounds.reserve(used_compounds.size());
    for (const auto& compound : used_compounds)
    {
      selected_compounds.insert(compound.id);
    }

    for (const auto& transition : all_transitions)
    {
      if (selected_compounds.find(transition.getPeptideRef()) != selected_compounds.end())
      {
        output.push_back(transition);
      }
    }
  }

  static void selectCompoundsForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all,
                                       OpenSwath::LightTargetedExperiment& transition_exp_used,
                                       const int batch_size,
                                       const size_t batch_idx)
  {
    const size_t start = batch_idx * batch_size;
    size_t end = start + batch_size;
    if (end > transition_exp_used_all.compounds.size())
    {
      end = transition_exp_used_all.compounds.size();
    }

    const Size compound_count = end - start;
    transition_exp_used.proteins = transition_exp_used_all.proteins;
    transition_exp_used.compounds.reserve(compound_count);
    transition_exp_used.compounds.insert(transition_exp_used.compounds.end(),
                                         transition_exp_used_all.compounds.begin() + start,
                                         transition_exp_used_all.compounds.begin() + end);
    if (!transition_exp_used_all.compounds.empty())
    {
      const Size transitions_per_compound =
        std::max<Size>(1, transition_exp_used_all.transitions.size() / transition_exp_used_all.compounds.size());
      transition_exp_used.transitions.reserve(compound_count * transitions_per_compound);
    }
    copyBatchTransitions_(transition_exp_used.compounds, transition_exp_used_all.transitions, transition_exp_used.transitions);
  }

  static OpenSwath::SpectrumAccessPtr loadMS1Map_(const std::vector<OpenSwath::SwathMap>& swath_maps, const bool load_into_memory)
  {
    OpenSwath::SpectrumAccessPtr ms1_map;
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      if (swath_maps[i].ms1)
      {
        ms1_map = swath_maps[i].sptr;
      }
    }
    if (load_into_memory && ms1_map)
    {
      ms1_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>(new SpectrumAccessOpenMSInMemory(*ms1_map));
    }
    return ms1_map;
  }

  static void selectBestPrmTransitions_(const OpenSwath::LightTargetedExperiment& transition_exp,
                                        std::vector<int>& tr_win_map,
                                        const std::vector<OpenSwath::SwathMap>& swath_maps,
                                        const ChromExtractParams& cp)
  {
    tr_win_map.resize(transition_exp.transitions.size(), -1);
    for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
    {
      for (Size k = 0; k < transition_exp.transitions.size(); ++k)
      {
        const auto& tr = transition_exp.transitions[k];
        if (swath_maps[i].lower < tr.getPrecursorMZ() && tr.getPrecursorMZ() < swath_maps[i].upper &&
            std::fabs(swath_maps[i].upper - tr.getPrecursorMZ()) >= cp.min_upper_edge_dist)
        {
          if (tr_win_map[k] == -1 ||
              std::fabs(swath_maps[tr_win_map[k]].center - tr.getPrecursorMZ()) >
              std::fabs(swath_maps[i].center - tr.getPrecursorMZ()))
          {
            tr_win_map[k] = i;
          }
        }
      }
    }
  }

  static std::unordered_set<std::string> loadPriorityPeptideSequences_(const std::vector<String>& tsv_files,
                                                                       const Param& tsv_reader_param)
  {
    std::unordered_set<std::string> priority_sequences;
    for (const auto& tsv_file : tsv_files)
    {
      if (tsv_file.empty() || !File::exists(tsv_file))
      {
        continue;
      }

      try
      {
        FileTypes::Type file_type = FileHandler::getType(tsv_file);
        TransitionTSVFile tsv_reader;
        tsv_reader.setParameters(tsv_reader_param);
        OpenSwath::LightTargetedExperiment priority_exp;
        tsv_reader.convertTSVToTargetedExperiment(tsv_file.c_str(), file_type, priority_exp);
        for (const auto& compound : priority_exp.getCompounds())
        {
          if (!compound.sequence.empty())
          {
            priority_sequences.insert(compound.sequence);
          }
        }
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Failed to load priority peptide file " << tsv_file
                        << ": " << e.what() << std::endl;
      }
    }
    return priority_sequences;
  }

  static String fileBasenameNoExtension_(const String& path)
  {
    return FileHandler::stripExtension(File::basename(path));
  }

  static String deriveRunOutputPath_(const String& out_dir, const StringList& run_files, const Size run_index)
  {
    const String basename = run_files.empty() ? String("run_") + String(run_index + 1) :
      fileBasenameNoExtension_(run_files.front());
    return File::absolutePath(out_dir).ensureLastChar('/') + basename + ".xipm";
  }

  static String prefixOutputPath_(const String& filename, const String& prefix)
  {
    if (filename.empty())
    {
      return "";
    }

    const String directory = File::path(filename);
    const String basename = File::basename(filename);
    if (directory == ".")
    {
      return prefix + "_" + basename;
    }
    return directory + "/" + prefix + "_" + basename;
  }

  static void writeExtractedPeakMaps_(const OpenSwath::SpectrumAccessPtr& input_map,
                                      const OpenSwath::LightTargetedExperiment& transition_exp_used,
                                      const TransformationDescription& trafo_inverse,
                                      const ChromExtractParams& cp,
                                      const bool ms1,
                                      const int ms1_isotopes,
                                      const UInt64 run_id,
                                      const String& source_file,
                                      XIPMParquetConsumer& consumer,
                                      Size& written_rows)
  {
    std::vector<PeakMapExtractor::ExtractedPeakMap> peak_maps;
    std::vector<PeakMapExtractor::ExtractionCoordinates> coordinates;
    prepareExtractionCoordinates_(coordinates, transition_exp_used, trafo_inverse, cp, ms1, ms1_isotopes);

    PeakMapExtractor extractor;
    extractor.extractPeakMaps(input_map, peak_maps, coordinates, cp.mz_extraction_window, cp.ppm, cp.im_extraction_window, cp.extraction_function);

    const Int64 ms_level = ms1 ? 1 : 2;
    for (const auto& peak_map : peak_maps)
    {
      if (peak_map.mz.empty())
      {
        continue;
      }
      consumer.consumePeakMap(peak_map, run_id, source_file, ms_level);
      ++written_rows;
    }
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList file_list = getStringList_("in");
    const String tr_file = getStringOption_("tr");
    const String out = getStringOption_("out");
    const bool separate_runs = getFlag_("separate_runs");
    const String out_dir = separate_runs ? getOutputDirOption("out_dir") : String();
    const bool split_file = getFlag_("split_file_input");
    const bool sort_swath_maps = getFlag_("sort_swath_maps");
    const bool force = getFlag_("force");
    const bool user_pasef = getFlag_("pasef");
    const bool use_ms1_traces = getStringOption_("enable_ms1") == "true";
    const bool use_ms1_im = getStringOption_("use_ms1_ion_mobility") == "true";
    const bool prm = getStringOption_("matching_window_only") == "true";
    const int ms1_isotopes = static_cast<int>(getIntOption_("ms1_isotopes"));
    const int batch_size_option = static_cast<int>(getIntOption_("batchSize"));
    const String swath_windows_file = getStringOption_("swath_windows_file");
    const double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    const String readoptions_raw = getStringOption_("readOptions");
    const bool keep_cached_files = getFlag_("keep_cached_files");
    const String tmp_dir = File::absolutePath(getStringOption_("tempDirectory")).ensureLastChar('/');

    if (separate_runs)
    {
      if (out_dir.empty())
      {
        writeLogError_("Parameter error: -out_dir is required when -separate_runs is set.");
        return ILLEGAL_PARAMETERS;
      }
      if (!out.empty())
      {
        OPENMS_LOG_WARN << "Both -separate_runs and -out were provided. -out will be ignored in favor of per-run files in -out_dir." << std::endl;
      }
    }
    else if (out.empty())
    {
      writeLogError_("Parameter error: -out is required unless -separate_runs is set.");
      return ILLEGAL_PARAMETERS;
    }

    if (prm && user_pasef)
    {
      writeLogError_("Setting -pasef and -matching_window_only simultaneously is not supported.");
      return ILLEGAL_PARAMETERS;
    }

    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = FileHandler::getType(tr_file);
    }
    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type for '-tr'.");
      return PARSE_ERROR;
    }

    bool load_into_memory = false;
    String readoptions = readoptions_raw;
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }
    else if (readoptions == "workingInMemory")
    {
      readoptions = "normal";
      load_into_memory = true;
    }

    ChromExtractParams cp;
    cp.min_upper_edge_dist = min_upper_edge_dist;
    cp.mz_extraction_window = getDoubleOption_("mz_extraction_window");
    cp.ppm = getStringOption_("mz_extraction_window_unit") == "ppm";
    cp.rt_extraction_window = getDoubleOption_("rt_extraction_window");
    cp.im_extraction_window = getDoubleOption_("ion_mobility_window");
    cp.extraction_function = getStringOption_("extraction_function");
    cp.extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1;
    cp_irt.mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
    cp_irt.im_extraction_window = getDoubleOption_("irt_im_extraction_window");
    cp_irt.ppm = getStringOption_("irt_mz_extraction_window_unit") == "ppm";

    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window = getDoubleOption_("mz_extraction_window_ms1");
    cp_ms1.ppm = getStringOption_("mz_extraction_window_ms1_unit") == "ppm";
    cp_ms1.im_extraction_window = use_ms1_im ? getDoubleOption_("im_extraction_window_ms1") : -1.0;

    Param tsv_reader_param = getParam_().copy("Library:", true);
    Param irt_calibration_params = getParam_().copy("Calibration:", true);
    Param irt_detection_param = getParam_().copy("Calibration:RTNormalization:", true);
    Param calibration_param = getParam_().copy("Calibration:MassIMCorrection:", true);
    Param calibration_workflow_param = irt_calibration_params.copy("", false);
    calibration_workflow_param.remove("MassIMCorrection:");
    calibration_workflow_param.remove("RTNormalization:");
    calibration_workflow_param.remove("rt_norm");
    calibration_workflow_param.remove("tr_irt_priority_sampling");
    calibration_param.setValue("mz_extraction_window", cp_irt.mz_extraction_window);
    calibration_param.setValue("mz_extraction_window_ppm", cp_irt.ppm ? "true" : "false");
    calibration_param.setValue("im_extraction_window", cp_irt.im_extraction_window);

    Param tmp_mrm_map_param = getParam_().copy("MRMMapping:", true);
    Param irt_mrm_map_param = OpenMS::MRMMapping().getDefaults();
    irt_mrm_map_param.setValue("precursor_tolerance", tmp_mrm_map_param.getValue("irt_precursor_tolerance"));
    irt_mrm_map_param.setValue("product_tolerance", tmp_mrm_map_param.getValue("irt_product_tolerance"));
    irt_mrm_map_param.setValue("map_multiple_assays", tmp_mrm_map_param.getValue("map_multiple_assays"));
    irt_mrm_map_param.setValue("error_on_unmapped", tmp_mrm_map_param.getValue("error_on_unmapped"));

    const Param calibration_feature_finder_param = makeCalibrationFeatureFinderDefaults_(use_ms1_im);

    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(tr_type, tr_file, tsv_reader_param);
    OPENMS_LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, "
                    << transition_exp.getCompounds().size() << " compounds with "
                    << transition_exp.getTransitions().size() << " transitions." << std::endl;

    std::vector<StringList> run_groups;
    if (split_file)
    {
      run_groups.push_back(file_list);
    }
    else
    {
      for (const String& file : file_list)
      {
        run_groups.push_back({file});
      }
    }

    String trafo_in = irt_calibration_params.getValue("rt_norm").toString();
    CalibrationWorkflow calibration_workflow;
    calibration_workflow.setLogType(log_type_);
    calibration_workflow.setParameters(calibration_workflow_param);
    const IrtStrategy irt_strategy = trafo_in.empty() ?
      calibration_workflow.determineIrtStrategy(transition_exp, run_groups.size()) :
      IrtStrategy::NULL_TRANSFORMATION;
    CalibrationWorkflow::IrtExperiments cached_irts;

    std::vector<String> priority_peptides;
    if (irt_calibration_params.getValue("auto_irt:enabled").toBool() && trafo_in.empty())
    {
      String data_path = File::getOpenMSDataPath();
      std::vector<String> priority_files;
      const String irtkit_path = data_path + "/CHEMISTRY/irtkit.tsv";
      const String cirtkit_path = data_path + "/CHEMISTRY/cirtkit.tsv";
      if (File::exists(irtkit_path)) priority_files.push_back(irtkit_path);
      if (File::exists(cirtkit_path)) priority_files.push_back(cirtkit_path);

      const String custom_priority_file = irt_calibration_params.getValue("tr_irt_priority_sampling").toString();
      if (!custom_priority_file.empty())
      {
        priority_files.push_back(custom_priority_file);
      }

      auto priority_set = loadPriorityPeptideSequences_(priority_files, tsv_reader_param);
      priority_peptides.reserve(priority_set.size());
      for (const auto& sequence : priority_set)
      {
        priority_peptides.push_back(sequence);
      }
    }

    std::unique_ptr<XIPMParquetConsumer> aggregate_writer;
    if (!separate_runs)
    {
      aggregate_writer = std::make_unique<XIPMParquetConsumer>(out, transition_exp);
      const Size expected_rows = run_groups.size() *
        (transition_exp.getTransitions().size() +
         (use_ms1_traces ? transition_exp.getCompounds().size() * static_cast<Size>(1 + ms1_isotopes) : 0));
      aggregate_writer->setExpectedSize(expected_rows);
    }

    for (Size run_index = 0; run_index < run_groups.size(); ++run_index)
    {
      const StringList& current_run_files = run_groups[run_index];
      OPENMS_LOG_INFO << "Processing run " << (run_index + 1) << "/" << run_groups.size() << std::endl;

      std::unique_ptr<File::TempDir> per_run_temp_dir;
      String per_run_tmp = tmp_dir;
      if (readoptions == "cache")
      {
        per_run_temp_dir = std::make_unique<File::TempDir>(tmp_dir, keep_cached_files);
        per_run_tmp = per_run_temp_dir->getPath();
      }

      std::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
      std::vector<OpenSwath::SwathMap> swath_maps;
      std::vector<String> swath_map_sources;
      if (!loadSwathFiles(current_run_files, exp_meta, swath_maps, swath_map_sources, split_file,
                          per_run_tmp, readoptions, swath_windows_file, min_upper_edge_dist,
                          force, sort_swath_maps, prm))
      {
        return PARSE_ERROR;
      }

      const bool has_spectra_based_map = std::any_of(swath_maps.begin(), swath_maps.end(),
        [](const OpenSwath::SwathMap& map) { return map.sptr && map.sptr->getNrSpectra() > 0; });
      if (!has_spectra_based_map)
      {
        writeLogError_("OpenSwathPeakMapExtractor requires spectra-based DIA/PRM input, not chromatogram-only SRM/MRM input.");
        return INCOMPATIBLE_INPUT_DATA;
      }

      bool pasef = user_pasef;
      if (!pasef)
      {
        pasef = std::any_of(swath_maps.begin(), swath_maps.end(),
          [](const OpenSwath::SwathMap& m) { return !m.ms1 && m.imLower >= 0 && m.imUpper >= 0; });
        if (pasef)
        {
          OPENMS_LOG_INFO << "Auto-detected ion mobility (PASEF) data from SWATH windows." << std::endl;
        }
      }
      if (pasef && prm)
      {
        writeLogError_("OpenSwathPeakMapExtractor does not support using PASEF window matching together with -matching_window_only.");
        return ILLEGAL_PARAMETERS;
      }

      ChromExtractParams cp_current = cp;
      ChromExtractParams cp_ms1_current = cp_ms1;
      OpenSwath::LightTargetedExperiment transition_exp_run = transition_exp;
      TransformationDescription rt_trafo;
      if (!trafo_in.empty())
      {
        FileHandler().loadTransformations(trafo_in, rt_trafo, false, {FileTypes::TRANSFORMATIONXML});
        Param model_params;
        model_params.setValue("symmetric_regression", "false");
        model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
        model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
        rt_trafo.fitModel(irt_detection_param.getValue("alignmentMethod").toString(), model_params);
        OPENMS_LOG_WARN << "Using existing RT transformation; m/z and ion mobility calibration will be skipped." << std::endl;
      }
      else
      {
        auto irt_experiments = calibration_workflow.prepareIrtExperiments(
          irt_strategy, transition_exp_run, priority_peptides, run_index,
          cached_irts.is_prepared ? &cached_irts : nullptr);
        if (irt_strategy == IrtStrategy::SAMPLE_ONCE && irt_experiments.is_prepared && !cached_irts.is_prepared)
        {
          cached_irts = irt_experiments;
        }

        String irt_trafo_out = getStringOption_("Debugging:irt_trafo");
        String irt_mzml_out = getStringOption_("Debugging:irt_mzml");
        if (run_groups.size() > 1)
        {
          const String prefix = current_run_files.empty() ? String("run_") + String(run_index + 1) : fileBasenameNoExtension_(current_run_files.front());
          irt_trafo_out = prefixOutputPath_(irt_trafo_out, prefix);
          irt_mzml_out = prefixOutputPath_(irt_mzml_out, prefix);
        }

        auto calibration_result = calibration_workflow.performCalibration(
          swath_maps, transition_exp_run, cp_current, cp_ms1_current, irt_experiments,
          calibration_feature_finder_param, cp_irt, irt_detection_param,
          calibration_param, irt_mrm_map_param, pasef, load_into_memory,
          irt_trafo_out, irt_mzml_out, static_cast<Size>(getIntOption_("debug")));
        rt_trafo = calibration_result.rt_trafo;
      }

      if (pasef)
      {
        for (const auto& transition : transition_exp_run.getTransitions())
        {
          if (transition.getPrecursorIM() < 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Transition " + String(transition.getNativeID()) + " is missing a precursor ion mobility value required for PASEF extraction.");
          }
        }
      }

      TransformationDescription trafo_inverse = rt_trafo;
      trafo_inverse.invert();

      std::unique_ptr<XIPMParquetConsumer> run_writer;
      XIPMParquetConsumer* writer = aggregate_writer.get();
      if (separate_runs)
      {
        const String run_out = deriveRunOutputPath_(out_dir, current_run_files, run_index);
        run_writer = std::make_unique<XIPMParquetConsumer>(run_out, transition_exp_run);
        const Size expected_rows =
          transition_exp_run.getTransitions().size() +
          (use_ms1_traces ? transition_exp_run.getCompounds().size() * static_cast<Size>(1 + ms1_isotopes) : 0);
        run_writer->setExpectedSize(expected_rows);
        writer = run_writer.get();
      }

      const UInt64 run_id = static_cast<UInt64>(run_index + 1);
      Size written_rows = 0;

      if (use_ms1_traces)
      {
        OpenSwath::SpectrumAccessPtr ms1_map = loadMS1Map_(swath_maps, load_into_memory);
        if (ms1_map)
        {
          String ms1_source_file = current_run_files.empty() ? "" : ListUtils::concatenate(current_run_files, ";");
          for (Size i = 0; i < swath_maps.size() && i < swath_map_sources.size(); ++i)
          {
            if (swath_maps[i].ms1)
            {
              ms1_source_file = swath_map_sources[i];
              break;
            }
          }

          const Size n_compounds = transition_exp_run.getCompounds().size();
          const int batch_size = batch_size_option <= 0 ? static_cast<int>(n_compounds) :
            std::min(batch_size_option, static_cast<int>(n_compounds));
          const SignedSize nr_batches = batch_size > 0 ?
            static_cast<SignedSize>((n_compounds + batch_size - 1) / batch_size) : 0;
          for (SignedSize batch_idx = 0; batch_idx < nr_batches; ++batch_idx)
          {
            OpenSwath::LightTargetedExperiment transition_exp_used;
            selectCompoundsForBatch_(transition_exp_run, transition_exp_used, batch_size, static_cast<size_t>(batch_idx));
            writeExtractedPeakMaps_(ms1_map, transition_exp_used, trafo_inverse, cp_ms1_current, true,
                                    ms1_isotopes, run_id, ms1_source_file, *writer, written_rows);
          }
        }
      }

      std::vector<int> tr_win_map;
      if (pasef)
      {
        OpenSwathHelper::selectSwathTransitionsPasef(transition_exp_run, tr_win_map, cp_current.min_upper_edge_dist, swath_maps);
      }
      else if (prm)
      {
        selectBestPrmTransitions_(transition_exp_run, tr_win_map, swath_maps, cp_current);
      }

      for (SignedSize swath_index = 0; swath_index < boost::numeric_cast<SignedSize>(swath_maps.size()); ++swath_index)
      {
        if (swath_maps[swath_index].ms1)
        {
          continue;
        }

        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        if (!(prm || pasef))
        {
          OpenSwathHelper::selectSwathTransitions(transition_exp_run, transition_exp_used_all,
                                                  cp_current.min_upper_edge_dist,
                                                  swath_maps[swath_index].lower,
                                                  swath_maps[swath_index].upper);
        }
        else
        {
          std::set<std::string> matching_compounds;
          for (Size k = 0; k < tr_win_map.size(); ++k)
          {
            if (tr_win_map[k] == swath_index)
            {
              const auto& tr = transition_exp_run.transitions[k];
              transition_exp_used_all.transitions.push_back(tr);
              matching_compounds.insert(tr.getPeptideRef());
            }
          }

          std::set<std::string> matching_proteins;
          for (const auto& compound : transition_exp_run.compounds)
          {
            if (matching_compounds.find(compound.id) != matching_compounds.end())
            {
              transition_exp_used_all.compounds.push_back(compound);
              for (const auto& protein_ref : compound.protein_refs)
              {
                matching_proteins.insert(protein_ref);
              }
            }
          }

          for (const auto& protein : transition_exp_run.proteins)
          {
            if (matching_proteins.find(protein.id) != matching_proteins.end())
            {
              transition_exp_used_all.proteins.push_back(protein);
            }
          }
        }

        if (transition_exp_used_all.getTransitions().empty())
        {
          continue;
        }

        OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[swath_index].sptr;
        if (load_into_memory)
        {
          current_swath_map = std::shared_ptr<SpectrumAccessOpenMSInMemory>(new SpectrumAccessOpenMSInMemory(*current_swath_map));
        }

        const Size n_compounds = transition_exp_used_all.getCompounds().size();
        const int batch_size = batch_size_option <= 0 ? static_cast<int>(n_compounds) :
          std::min(batch_size_option, static_cast<int>(n_compounds));
        const SignedSize nr_batches = batch_size > 0 ?
          static_cast<SignedSize>((n_compounds + batch_size - 1) / batch_size) : 0;
        const String swath_source_file =
          swath_index < boost::numeric_cast<SignedSize>(swath_map_sources.size()) ?
          swath_map_sources[swath_index] :
          ListUtils::concatenate(current_run_files, ";");

        for (SignedSize batch_idx = 0; batch_idx < nr_batches; ++batch_idx)
        {
          OpenSwath::LightTargetedExperiment transition_exp_used;
          selectCompoundsForBatch_(transition_exp_used_all, transition_exp_used, batch_size, static_cast<size_t>(batch_idx));
          writeExtractedPeakMaps_(current_swath_map, transition_exp_used, trafo_inverse, cp_current, false,
                                  0, run_id, swath_source_file, *writer, written_rows);
        }
      }

      if (run_writer)
      {
        run_writer->finalize();
      }
      OPENMS_LOG_INFO << "Wrote " << written_rows << " extracted peak maps for run "
                      << (run_index + 1) << std::endl;
    }

    if (aggregate_writer)
    {
      aggregate_writer->finalize();
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenSwathPeakMapExtractor tool;
  return tool.main(argc, argv);
}
/// @endcond
