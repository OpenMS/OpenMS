// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACNN.h"
#include "OpenDIACNNQValues.h"
#include "OpenDIATraces.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
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
      registerStringOption_("mode", "<train|score|crossfit>", "train", "Train, apply a frozen checkpoint, or evaluate family-disjoint outer folds", false);
      setValidStrings_("mode", {"train", "score", "crossfit"});
      registerStringOption_("model", "<file>", "", "Checkpoint in train/score mode; new checkpoint prefix in crossfit mode", true);
      registerOutputFile_("out", "<file>", "", "Candidate logits and provenance; crossfit adds experimental fold-wise precursor q-values");
      setValidFormats_("out", {"tsv"});
      registerIntOption_("max_samples", "<n>", 16384, "Reject larger window manifests rather than silently truncating; see max_tensor_gib for the raw tensor cap", false);
      setMinInt_("max_samples", 1);
      registerIntOption_("max_tensor_gib", "<GiB>", 2, "Hard cap on the in-memory raw trace tensor (6 KiB per window)", false);
      setMinInt_("max_tensor_gib", 1);
      setMaxInt_("max_tensor_gib", 1024);
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
      registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0, e.g. unitless mzML 'Ion Mobility' arrays");
      registerInputFile_("seeds", "<file>", "", "Crossfit only: explicit sampled training windows; labels in -windows are ignored", false);
      setValidFormats_("seeds", {"tsv"});
      registerIntOption_("folds", "<n>", 3, "Crossfit outer family folds; each excludes fitting and early stopping families", false);
      setMinInt_("folds", 2);
      setMaxInt_("folds", 10);
      registerDoubleOption_("target_decoy_ratio", "<ratio>", 1.0, "Crossfit decoy-count correction: searched library target/decoy ratio, never observed winner ratio", false);
      setMinFloat_("target_decoy_ratio", 0.0);
    }

    ExitCodes main_(int, const char**) override
    {
      try
      {
        const std::string mode = getStringOption_("mode");
        const bool training = mode == "train";
        const bool crossfit = mode == "crossfit";
        const std::string checkpoint = getStringOption_("model");
        const std::string output = getStringOption_("out");
        const auto same_path = [](const std::filesystem::path& a, const std::filesystem::path& b)
        {
          return std::filesystem::weakly_canonical(a) == std::filesystem::weakly_canonical(b) ||
                 (std::filesystem::exists(a) && std::filesystem::exists(b) && std::filesystem::equivalent(a, b));
        };
        std::vector<std::string> model_paths{checkpoint};
        if (crossfit)
        {
          model_paths.clear();
          for (int fold = 0; fold < getIntOption_("folds"); ++fold)
          {
            model_paths.push_back(checkpoint + ".fold" + std::to_string(fold) + ".pt");
            if (std::filesystem::exists(model_paths.back()))
            {
              throw std::invalid_argument("Crossfit refuses to overwrite existing fold checkpoints. Use a new model prefix.");
            }
          }
        }
        // Validate every derived checkpoint path before opening any outputs.
        for (const auto& path : model_paths)
        {
          if (same_path(output, path)) { throw std::invalid_argument("Model and score output must be distinct."); }
          for (const auto& option : {"in", "lib", "windows", "seeds"})
          {
            const std::string input = getStringOption_(option);
            if (input.empty()) { continue; }
            const auto inside_input = [&](const std::string& candidate)
            {
              const auto relative = std::filesystem::weakly_canonical(candidate).lexically_relative(std::filesystem::weakly_canonical(input));
              return std::filesystem::is_directory(input) && !relative.empty() && *relative.begin() != "..";
            };
            if (same_path(output, input) || same_path(path, input) || inside_input(output) || inside_input(path))
            {
              throw std::invalid_argument("Model and score outputs must not overwrite raw, library or window inputs.");
            }
          }
        }
        if (!training && !crossfit && !std::filesystem::is_regular_file(checkpoint))
        {
          throw std::invalid_argument("Score mode requires an existing native checkpoint.");
        }
        if ((training || crossfit) && getFlag_("allow_training_overlap"))
        {
          throw std::invalid_argument("allow_training_overlap is only meaningful in score mode.");
        }
        if (crossfit != !getStringOption_("seeds").empty())
        {
          throw std::invalid_argument("Crossfit requires -seeds; other modes do not accept -seeds.");
        }
        const double target_decoy_ratio = getDoubleOption_("target_decoy_ratio");
        if (!std::isfinite(target_decoy_ratio) || target_decoy_ratio <= 0)
        {
          throw std::invalid_argument("target_decoy_ratio must be finite and positive.");
        }

        const auto max_samples = static_cast<std::size_t>(getIntOption_("max_samples"));
        auto windows = OpenDIA::readTraceWindows(getStringOption_("windows"), max_samples);
        if (crossfit)
        {
          std::map<std::string, OpenDIA::TraceWindow> seeds;
          for (auto& seed : OpenDIA::readTraceWindows(getStringOption_("seeds"), max_samples))
          {
            if (seed.label < 0) { throw std::invalid_argument("Crossfit seeds must have explicit positive or negative labels."); }
            seeds.emplace(seed.sample_id, std::move(seed));
          }
          for (auto& window : windows)
          {
            const auto found = seeds.find(window.sample_id);
            window.label = -1;
            if (found == seeds.end()) { continue; }
            const auto& seed = found->second;
            if (seed.precursor_id != window.precursor_id || seed.group_id != window.group_id || seed.apex_rt != window.apex_rt)
            {
              throw std::invalid_argument("Seed and scoring manifest identities, families and RT coordinates must match exactly.");
            }
            window.label = seed.label;
            seeds.erase(found);
          }
          if (!seeds.empty()) { throw std::invalid_argument("Every seed must occur in the scoring manifest."); }
        }
        const OpenDIA::CNNGeometry geometry{getDoubleOption_("mz_tolerance"), getDoubleOption_("rt_half_width"), getDoubleOption_("im_half_width")};
        OpenSwath::LightTargetedExperiment library;
        TransitionPQPFile().convertPQPToTargetedExperiment(getStringOption_("lib").c_str(), library);
        OpenDIA::TraceCollector collector(std::move(windows), library, geometry.mz_tolerance,
                                         geometry.rt_half_width, geometry.im_half_width, max_samples,
                                         static_cast<std::size_t>(getIntOption_("max_tensor_gib")) << 30);
        collector.setIonMobilityArrayAsVSSC(getFlag_("im_array_as_vssc"));
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
        OpenDIA::CNNOptions settings;
        settings.seed = getIntOption_("seed");
        settings.threads = getIntOption_("threads");
        settings.epochs = getIntOption_("epochs");
        settings.patience = getIntOption_("patience");
        settings.batch_size = getIntOption_("batch_size");
        std::vector<OpenDIA::CNNQValueInput> crossfit_rows;
        if (crossfit)
        {
          settings.outer_folds = getIntOption_("folds");
          result.logits.resize(count);
          result.splits.assign(count, "heldout");
          std::vector<bool> covered(count, false);
          for (int fold = 0; fold < settings.outer_folds; ++fold)
          {
            settings.heldout_fold = fold;
            const auto fitted = OpenDIA::trainTraceCNN(data, geometry, settings, model_paths.at(fold));
            std::size_t heldout_count = 0;
            for (std::size_t i = 0; i < count; ++i)
            {
              if (OpenDIA::traceCNNFold(data.samples[i].group_id, settings.seed, settings.outer_folds) != fold) { continue; }
              if (covered[i] || fitted.splits.at(i) != "heldout") { throw std::runtime_error("Invalid outer-fold prediction coverage."); }
              covered[i] = true;
              ++heldout_count;
              result.logits[i] = fitted.logits.at(i);
            }
            writeLogInfo_("Native CNN fold=" + std::to_string(fold) + " heldout_windows=" + std::to_string(heldout_count) +
                          " epochs=" + std::to_string(fitted.epochs_run) + " validation_loss=" + std::to_string(fitted.validation_loss) +
                          " heldout_seed_decoy_auc=" + std::to_string(fitted.heldout_auc));
          }
          for (std::size_t i = 0; i < count; ++i)
          {
            if (!covered[i]) { throw std::runtime_error("Missing outer-fold prediction."); }
            crossfit_rows.push_back({data.samples[i].sample_id, collector.window(i).precursor_id, data.samples[i].group_id,
              data.samples[i].is_decoy, OpenDIA::traceCNNFold(data.samples[i].group_id, settings.seed, settings.outer_folds), result.logits[i]});
          }
        }
        else if (training)
        {
          result = OpenDIA::trainTraceCNN(data, geometry, settings, checkpoint);
          writeLogInfo_("Native CNN epochs=" + std::to_string(result.epochs_run) +
                        " validation_loss=" + std::to_string(result.validation_loss) +
                        " heldout_seed_decoy_auc=" + std::to_string(result.heldout_auc));
        }
        else
        {
          result = OpenDIA::scoreTraceCNN(data, geometry, checkpoint, getIntOption_("threads"), getFlag_("allow_training_overlap"));
        }
        // Finish validation/statistics before opening the score output.
        const auto evaluated = crossfit ? OpenDIA::aggregateCNNQValues(crossfit_rows, target_decoy_ratio) : std::vector<OpenDIA::CNNQValueResult>{};
        std::ofstream scores(output);
        if (crossfit)
        {
          scores << "sample_id\tprecursor_id\tgroup_id\tis_decoy\tfold\tlogit\tprecursor_best\tq_value\n"
                 << std::setprecision(std::numeric_limits<double>::max_digits10);
          std::size_t targets_at_one_percent = 0;
          for (const auto& row : evaluated)
          {
            scores << row.input.sample_id << '\t' << row.input.precursor_id << '\t' << row.input.group_id << '\t'
                   << row.input.is_decoy << '\t' << row.input.fold << '\t' << row.input.score << '\t' << row.precursor_best << '\t';
            if (row.precursor_best) { scores << row.q_value; }
            else { scores << "NA"; }
            scores << '\n';
            targets_at_one_percent += row.precursor_best && !row.input.is_decoy && row.q_value <= 0.01;
          }
          writeLogInfo_("Experimental fold-wise decoy-count estimate: target_decoy_ratio=" + std::to_string(target_decoy_ratio) +
                        " target_precursors_q01=" + std::to_string(targets_at_one_percent) +
                        ". Teacher-seeded, not independently validated FDR control.");
        }
        else
        {
          scores << "sample_id\tgroup_id\tlogit\tsplit\n" << std::setprecision(std::numeric_limits<float>::max_digits10);
          for (std::size_t i = 0; i < count; ++i)
          {
            scores << data.samples[i].sample_id << '\t' << data.samples[i].group_id << '\t'
                   << result.logits.at(i) << '\t' << result.splits.at(i) << '\n';
          }
        }
        scores.close();
        if (!scores) { throw std::runtime_error("Could not write native CNN candidate logits."); }
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        writeLogInfo_("Native CNN " + std::string(crossfit ? "crossfit" : training ? "training" : "inference") + " seconds=" +
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
