// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACNN.h"
#include "OpenDIATraces.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace OpenMS
{
  // Experimental companion: the default OpenDIA executable has no Torch dependency.
  class TOPPOpenDIACNN final : public TOPPBase
  {
  public:
    TOPPOpenDIACNN() : TOPPBase("OpenDIACNN", "Train or apply a native CPU CNN directly to raw XIC/XIM candidate tensors", {}, false) {}

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFile_("in", "<file>", "", "DIA raw or mzML input");
      registerInputFile_("lib", "<file>", "", "Prepared PQP library");
      setValidFormats_("lib", {"pqp"});
      registerInputFile_("windows", "<file>", "", "Explicit sample_id, precursor_id, group_id, label (-1/0/1), rt_apex TSV");
      setValidFormats_("windows", {"tsv"});
      registerStringOption_("mode", "<train|score>", "train", "Train on explicit seed labels or score with a frozen native checkpoint", false);
      setValidStrings_("mode", {"train", "score"});
      registerStringOption_("model", "<file>", "", "Native checkpoint to write in train mode or read in score mode", true);
      registerOutputFile_("out", "<file>", "", "Candidate logits and family-split provenance, not probabilities or q-values");
      setValidFormats_("out", {"tsv"});
      registerIntOption_("max_samples", "<n>", 16384, "Reject larger window manifests rather than silently truncating; raw tensor cap is 2 GiB", false);
      setMinInt_("max_samples", 1);
      registerIntOption_("threads", "<n>", 5, "LibTorch CPU threads", false);
      setMinInt_("threads", 1);
      registerIntOption_("seed", "<n>", 42, "Training and family-partition random seed", false);
      setMinInt_("seed", 0);
      registerIntOption_("epochs", "<n>", 40, "Training epoch cap", false);
      setMinInt_("epochs", 1);
      setMaxInt_("epochs", 100);
      registerIntOption_("patience", "<n>", 6, "Stop after this many validation epochs without improvement", false);
      setMinInt_("patience", 1);
      setMaxInt_("patience", 20);
      registerIntOption_("batch_size", "<n>", 32, "Training tensor batch size", false);
      setMinInt_("batch_size", 1);
      registerDoubleOption_("mz_tolerance", "<Da>", 0.02, "Fragment m/z half width", false);
      registerDoubleOption_("rt_half_width", "<seconds>", 15.0, "RT half width around each candidate apex", false);
      registerDoubleOption_("im_half_width", "<1/K0>", 0.03, "Mobility half width around library precursor mobility", false);
      registerFlag_("allow_training_overlap", "Diagnostic scoring only: permit train-family rows, retaining their train provenance");
    }

    ExitCodes main_(int, const char**) override
    {
      try
      {
        const bool training = getStringOption_("mode") == "train";
        const std::string checkpoint = getStringOption_("model");
        const std::string output = getStringOption_("out");
        const auto same_path = [](const std::filesystem::path& a, const std::filesystem::path& b)
        {
          return std::filesystem::weakly_canonical(a) == std::filesystem::weakly_canonical(b) ||
                 (std::filesystem::exists(a) && std::filesystem::exists(b) && std::filesystem::equivalent(a, b));
        };
        // Check aliases before opening any outputs, including hard links and model/output collisions.
        if (same_path(output, checkpoint)) { throw std::invalid_argument("Model and score output must be distinct."); }
        for (const auto& option : {"in", "lib", "windows"})
        {
          const std::string input = getStringOption_(option);
          const auto inside_input = [&](const std::string& path)
          {
            const auto relative = std::filesystem::weakly_canonical(path).lexically_relative(std::filesystem::weakly_canonical(input));
            return std::filesystem::is_directory(input) && !relative.empty() && *relative.begin() != "..";
          };
          if (same_path(output, input) || same_path(checkpoint, input) || inside_input(output) || inside_input(checkpoint))
          {
            throw std::invalid_argument("Model and score outputs must not overwrite raw, library or window inputs.");
          }
        }
        if (!training && !std::filesystem::is_regular_file(checkpoint))
        {
          throw std::invalid_argument("Score mode requires an existing native checkpoint.");
        }
        if (training && getFlag_("allow_training_overlap"))
        {
          throw std::invalid_argument("allow_training_overlap is only meaningful in score mode.");
        }

        const auto max_samples = static_cast<std::size_t>(getIntOption_("max_samples"));
        auto windows = OpenDIA::readTraceWindows(getStringOption_("windows"), max_samples);
        const OpenDIA::CNNGeometry geometry{getDoubleOption_("mz_tolerance"), getDoubleOption_("rt_half_width"), getDoubleOption_("im_half_width")};
        OpenSwath::LightTargetedExperiment library;
        TransitionPQPFile().convertPQPToTargetedExperiment(getStringOption_("lib").c_str(), library);
        OpenDIA::TraceCollector collector(std::move(windows), library, geometry.mz_tolerance,
                                         geometry.rt_half_width, geometry.im_half_width, max_samples);
        FileHandler reader;
        PeakFileOptions options;
        options.setMSLevels({2});
        reader.setOptions(options);
        reader.transformExperiment(getStringOption_("in"), &collector, {}, log_type_);

        // Borrow the collector's contiguous tensor. No NumPy, intermediate tensor file or Python process.
        OpenDIA::CNNData data;
        data.intensities = collector.intensities().data();
        data.intensity_count = collector.intensities().size();
        const std::size_t count = data.intensity_count / OpenDIA::TraceCollector::SAMPLE_SIZE;
        data.samples.reserve(count);
        for (std::size_t i = 0; i < count; ++i)
        {
          const auto& window = collector.window(i);
          OpenDIA::CNNSample sample;
          sample.sample_id = window.sample_id;
          sample.group_id = window.group_id;
          sample.label = window.label;
          sample.is_decoy = collector.isDecoy(i);
          sample.has_im = collector.hasIM(i);
          sample.library_intensities = collector.libraryIntensities(i);
          sample.fragment_mz = collector.fragmentMZ(i);
          data.samples.push_back(std::move(sample));
        }

        const auto start = std::chrono::steady_clock::now();
        OpenDIA::CNNResult result;
        if (training)
        {
          OpenDIA::CNNOptions settings;
          settings.seed = getIntOption_("seed");
          settings.threads = getIntOption_("threads");
          settings.epochs = getIntOption_("epochs");
          settings.patience = getIntOption_("patience");
          settings.batch_size = getIntOption_("batch_size");
          result = OpenDIA::trainTraceCNN(data, geometry, settings, checkpoint);
          writeLogInfo_("Native CNN epochs=" + std::to_string(result.epochs_run) +
                        " validation_loss=" + std::to_string(result.validation_loss) +
                        " heldout_seed_decoy_auc=" + std::to_string(result.heldout_auc));
        }
        else
        {
          result = OpenDIA::scoreTraceCNN(data, geometry, checkpoint, getIntOption_("threads"), getFlag_("allow_training_overlap"));
        }
        std::ofstream scores(output);
        scores << "sample_id\tgroup_id\tlogit\tsplit\n" << std::setprecision(std::numeric_limits<float>::max_digits10);
        for (std::size_t i = 0; i < count; ++i)
        {
          scores << data.samples[i].sample_id << '\t' << data.samples[i].group_id << '\t'
                 << result.logits.at(i) << '\t' << result.splits.at(i) << '\n';
        }
        scores.close();
        if (!scores) { throw std::runtime_error("Could not write native CNN candidate logits."); }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        writeLogInfo_("Native CNN " + std::string(training ? "training" : "inference") + " seconds=" +
                      std::to_string(seconds) + " windows=" + std::to_string(count));
        return EXECUTION_OK;
      }
      catch (const std::exception& error)
      {
        writeLogError_(error.what());
        return ILLEGAL_PARAMETERS;
      }
    }
  };
}

int main(int argc, const char** argv)
{
  OpenMS::TOPPOpenDIACNN tool;
  return tool.main(argc, argv);
}
