// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
#include "OpenDIATraces.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <filesystem>

namespace OpenMS
{
  // Experimental companion, not part of the stable TOPP tool registry.
  class TOPPOpenDIATraces final : public TOPPBase
  {
  public:
    TOPPOpenDIATraces() : TOPPBase("OpenDIATraces", "Export sampled raw XIC/XIM tensors for experimental neural scoring", {}, false) {}

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "DIA raw or mzML input");
      registerInputFile_("lib", "<file>", "", "Prepared PQP library");
      setValidFormats_("lib", {"pqp"});
      registerInputFile_("windows", "<file>", "", "Explicit sample_id, precursor_id, group_id, label (-1/0/1), rt_apex TSV");
      setValidFormats_("windows", {"tsv"});
      registerStringOption_("out_prefix", "<prefix>", "", "Output prefix for .traces.npy, .samples.tsv and .geometry.json", true);
      registerIntOption_("max_samples", "<n>", 16384, "Reject larger input samples rather than silently truncating", false);
      setMinInt_("max_samples", 1);
      registerDoubleOption_("mz_tolerance", "<Da>", 0.02, "Fragment m/z half width", false);
      registerDoubleOption_("rt_half_width", "<seconds>", 15.0, "RT half width around each candidate apex", false);
      registerDoubleOption_("im_half_width", "<1/K0>", 0.03, "Mobility half width around library precursor mobility", false);
    }

    ExitCodes main_(int, const char**) override
    {
      const std::string prefix = getStringOption_("out_prefix");
      for (const auto& suffix : {".traces.npy", ".samples.tsv", ".geometry.json"})
      {
        const std::filesystem::path output(prefix + suffix);
        for (const auto& option : {"in", "lib", "windows"})
        {
          const std::filesystem::path input(getStringOption_(option));
          if (std::filesystem::weakly_canonical(output) == std::filesystem::weakly_canonical(input) ||
              (std::filesystem::exists(output) && std::filesystem::equivalent(input, output)))
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Trace output would overwrite an input file.");
          }
        }
      }
      const auto max_samples = static_cast<std::size_t>(getIntOption_("max_samples"));
      auto windows = OpenDIA::readTraceWindows(getStringOption_("windows"), max_samples);
      OpenSwath::LightTargetedExperiment library;
      TransitionPQPFile().convertPQPToTargetedExperiment(getStringOption_("lib").c_str(), library);
      OpenDIA::TraceCollector collector(std::move(windows), library, getDoubleOption_("mz_tolerance"),
                                       getDoubleOption_("rt_half_width"), getDoubleOption_("im_half_width"), max_samples);
      FileHandler reader;
      PeakFileOptions options;
      options.setMSLevels({2});
      reader.setOptions(options);
      reader.transformExperiment(getStringOption_("in"), &collector, {}, log_type_);
      collector.store(prefix);
      return EXECUTION_OK;
    }
  };
}

int main(int argc, const char** argv)
{
  OpenMS::TOPPOpenDIATraces tool;
  return tool.main(argc, argv);
}
