// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIACalibration.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <optional>

namespace OpenMS
{
  // Experimental companion, not part of the stable TOPP tool registry.
  class TOPPOpenDIACalibrate final : public TOPPBase
  {
  public:
    TOPPOpenDIACalibrate() :
      TOPPBase("OpenDIACalibrate", "Identification-free RT, fragment m/z and 1/K0 calibration of a library to one DIA run", {}, false)
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "DIA raw or mzML input");
      registerInputFile_("lib", "<file>", "", "Prepared PQP library with decoys named DECOY_<target TRAML_ID>");
      setValidFormats_("lib", {"pqp"});
      registerOutputFile_("out", "<file>", "", "Calibration JSON: RT map, mass offset, 1/K0 model and tolerances");
      registerOutputFile_("out_anchors", "<file>", "", "Optional anchor TSV", false);
      registerIntOption_("panel_targets", "<n>", 10000, "Targets sampled evenly over library RT; each brings its paired decoy", false);
      setMinInt_("panel_targets", 50);
      registerIntOption_("max_panel_targets", "<n>", 0, "Experimental: if calibration fails for too few anchors, retry with a 4x larger panel up to this many targets (one more raw pass per retry; memory grows ~75 KB per pair). 0 disables", false);
      setMinInt_("max_panel_targets", 0);
      registerDoubleOption_("search_ppm", "<ppm>", 15.0, "Uncalibrated fragment search half width", false);
      registerDoubleOption_("search_im", "<1/K0>", 0.05, "Uncalibrated 1/K0 search half width around the library mobility", false);
      registerDoubleOption_("anchor_fdr", "<fraction>", 0.01, "Paired target-decoy q-value cut for anchors", false);
      registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0, e.g. unitless mzML 'Ion Mobility' arrays");
    }

    ExitCodes main_(int, const char**) override
    {
      const std::string out = getStringOption_("out");
      const std::string anchors = getStringOption_("out_anchors");
      for (const auto& output : {out, anchors})
      {
        if (output.empty())
        {
          continue;
        }
        for (const auto& option : {"in", "lib"})
        {
          const std::filesystem::path input(getStringOption_(option));
          if (std::filesystem::weakly_canonical(output) == std::filesystem::weakly_canonical(input))
          {
            writeLogError_("Calibration output would overwrite an input.");
            return ILLEGAL_PARAMETERS;
          }
        }
      }
      OpenDIA::CalibrationOptions options;
      options.panel_targets = static_cast<std::size_t>(getIntOption_("panel_targets"));
      options.search_ppm = getDoubleOption_("search_ppm");
      options.search_im = getDoubleOption_("search_im");
      options.anchor_fdr = getDoubleOption_("anchor_fdr");

      const std::string library_path = getStringOption_("lib");
      OpenSwath::LightTargetedExperiment library;
      TransitionPQPFile().convertPQPToTargetedExperiment(library_path.c_str(), library);
      const auto traml_map = TransitionPQPFile().getPQPCurrentIDToTraMLIDMap(library_path.c_str(), "PRECURSOR");
      const std::vector<std::pair<std::string, std::string>> traml_ids(traml_map.begin(), traml_map.end());

      const auto start = std::chrono::steady_clock::now();
      const auto max_panel = static_cast<std::size_t>(getIntOption_("max_panel_targets"));
      if (max_panel != 0 && max_panel <= options.panel_targets)
      {
        writeLogError_("max_panel_targets must exceed panel_targets (or be 0).");
        return ILLEGAL_PARAMETERS;
      }
      std::optional<OpenDIA::CalibrationResult> calibrated;
      while (!calibrated)
      {
        OpenDIA::CalibrationCollector collector(library, traml_ids, options);
        collector.setIonMobilityArrayAsVSSC(getFlag_("im_array_as_vssc"));
        FileHandler reader;
        PeakFileOptions read_options;
        read_options.setMSLevels({2});
        reader.setOptions(read_options);
        reader.transformExperiment(getStringOption_("in"), &collector, {}, log_type_);
        try
        {
          calibrated = collector.calibrate();
        }
        catch (const std::runtime_error& error)
        {
          // Escalate only the panel size; every other setting stays fixed. The library limits the panel size.
          const std::size_t panel_pairs = collector.panelSize() / 2;
          if (max_panel == 0 || options.panel_targets >= max_panel || panel_pairs < options.panel_targets)
          {
            throw;
          }
          writeLogWarn_("Calibration with panel_targets=" + std::to_string(options.panel_targets) + " failed (" + error.what() +
                        "); retrying with a larger panel.");
          options.panel_targets = std::min(max_panel, 4 * options.panel_targets);
        }
      }
      const auto& result = *calibrated;

      std::ofstream json(out);
      OpenDIA::writeCalibrationJSON(result, options, json);
      json.close();
      if (!json)
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      if (!anchors.empty())
      {
        std::ofstream table(anchors);
        OpenDIA::writeCalibrationAnchors(result, table);
        table.close();
        if (!table)
        {
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      writeLogInfo_("Calibration panel=" + std::to_string(result.panel_targets) + "+" + std::to_string(result.panel_decoys) +
                    " stage1_anchors=" + std::to_string(result.stage1_targets) + " (decoys above cut " + std::to_string(result.stage1_decoys) + ")" +
                    " em_iterations=" + std::to_string(result.em_iterations) + (result.em_converged ? " converged" : " NOT converged") +
                    " em_peak_responsibility=" + std::to_string(result.em_presence) + " anchors=" + std::to_string(result.anchor_targets) +
                    " rt_sigma_s=" + std::to_string(result.rt_map.sigma) + " rt_tolerance_s=" + std::to_string(result.rt_tolerance_s) +
                    " ms2_offset_ppm=" + std::to_string(result.ms2_offset_ppm) + " ms2_tolerance_ppm=" + std::to_string(result.ms2_tolerance_ppm) +
                    (result.has_im ? " im_offset_at_1=" + std::to_string(result.runIM(1.0) - 1.0) + " im_tolerance=" + std::to_string(result.im_tolerance)
                                   : std::string(" im=unavailable")) +
                    " seconds=" + std::to_string(seconds) + ". Experimental, identification-free heuristic calibration.");
      return EXECUTION_OK;
    }
  };
}

int main(int argc, const char** argv)
{
  OpenMS::TOPPOpenDIACalibrate tool;
  return tool.main(argc, argv);
}
