// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACalibration.h"
#include "OpenDIACNN.h"
#include "OpenDIACNNQValues.h"
#include "OpenDIACNNRTWindow.h"
#include "OpenDIAStage2.h"
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
#include <optional>
#include <set>
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
      registerFlag_("paired_tdc", "Crossfit only, experimental: paired target-decoy competition. A decoy precursor TRAML_ID 'DECOY_<X>' competes with target <X>; q=(D+1)/T per fold, so target_decoy_ratio must stay 1");
      registerFlag_("pooled_tdc", "With -paired_tdc, experimental: standardise each fold's logits by the median and MAD of its decoy precursor winners and run one paired competition over all folds instead of per fold (output logits unchanged)");
      registerFlag_("rt_window_calibration", "Crossfit only, experimental: fit a library-RT to run-RT map per fold from the other folds' top 500 target precursors by CNN score (no teacher IDs) and keep only windows within max(4 sigma, 30 s) of the map before precursor q-values");
      registerInputFile_("rt_calibration", "<file>", "", "Crossfit only, experimental: OpenDIACalibrate JSON (identification-free, CNN-independent). Keep only windows whose apex lies within the calibrated RT tolerance of the calibrated map before precursor q-values", false);
      setValidFormats_("rt_calibration", {"json"});
      registerFlag_("apply_calibration_coordinates", "With -rt_calibration, experimental: shift library fragment m/z by the calibrated ppm offset and precursor 1/K0 by the calibrated line before trace extraction (tolerances unchanged)");
      registerIntOption_("ss_iterations", "<n>", 0, "Crossfit only, experimental: PyProphet-style score iterations after the -seeds round. Each retrains every fold on target precursor winners at q <= ss_fdr against each decoy precursor's best window; 0 disables", false);
      setMinInt_("ss_iterations", 0);
      setMaxInt_("ss_iterations", 10);
      registerDoubleOption_("ss_fdr", "<q>", 0.01, "Score-iteration positive q-value cut", false);
      setMinFloat_("ss_fdr", 0.0001);
      setMaxFloat_("ss_fdr", 0.1);
      registerIntOption_("ss_max_seeds", "<n>", 6000, "Score-iteration cap per class (deterministic FNV-1a order)", false);
      setMinInt_("ss_max_seeds", 3);
      registerFlag_("stage2", "Crossfit with -paired_tdc -pooled_tdc, experimental: second-stage score. Reads MS1 as well, computes fragment co-elution, "
                    "library intensity correlation, fragment and MS1 mass error, MS1 co-elution and isotope presence per window, and combines them "
                    "with the fold-standardised CNN logit by a cross-fitted semi-supervised LDA (PyProphet-style). q-values then come from that score.");
      registerIntOption_("stage2_iterations", "<n>", 3, "Second-stage semi-supervised LDA rounds", false);
      setMinInt_("stage2_iterations", 1);
      setMaxInt_("stage2_iterations", 10);
      registerDoubleOption_("stage2_initial_fdr", "<q>", 0.15, "Second-stage round-1 positive cut, ranked by the CNN score", false);
      setMinFloat_("stage2_initial_fdr", 0.001);
      setMaxFloat_("stage2_initial_fdr", 0.5);
      registerDoubleOption_("stage2_iteration_fdr", "<q>", 0.05, "Second-stage positive cut in later rounds", false);
      setMinFloat_("stage2_iteration_fdr", 0.001);
      setMaxFloat_("stage2_iteration_fdr", 0.5);
      registerDoubleOption_("stage2_core_half_width", "<seconds>", 4.0, "Half width around each apex used by the second-stage features", false);
      setMinFloat_("stage2_core_half_width", 0.5);
      registerDoubleOption_("stage2_ms1_ppm", "<ppm>", 20.0, "Half width of the MS1 m/z windows for the second-stage features", false);
      setMinFloat_("stage2_ms1_ppm", 0.1);
      registerOutputFile_("stage2_features_out", "<file>", "", "With -stage2: per-window feature table (after the RT window) for inspection", false);
      registerStringOption_("stage2_feature_set", "<groups>", "extraction", "With -stage2: '+'-separated feature groups after the CNN logit. "
                            "'extraction' = the seven tensor/MS1 features; 'search' = search peak-group features (needs -stage2_search_features); "
                            "'isotope' = MS1 M+0..M+3 averagine correlation, M-1 overlap and isotope co-elution; 'boundary' (experimental) = each window re-picked "
                            "with the consensus peak-group picker on its scan bins, with FWHM, base width, co-elution, library correlation and fragment apex agreement "
                            "over the picked boundaries; 'byion' (experimental) = numbers of singly charged b and y ions of the library sequence matched "
                            "in the isolating MS2 spectrum nearest the apex (-stage2_byion_ppm, -stage2_byion_intensity_min); 'shape' (experimental) = lagged "
                            "fragment cross-correlation, library dot product / Manhattan / RMSD / spectral angle and signal-to-noise over the boundary-group peak. "
                            "'extraction_search' = extraction+search", false);
      registerDoubleOption_("stage2_byion_ppm", "<ppm>", 10.0, "With a stage2_feature_set containing byion: matching half width of the b/y ions", false);
      setMinFloat_("stage2_byion_ppm", 0.1);
      registerDoubleOption_("stage2_byion_intensity_min", "<intensity>", 300.0, "With a stage2_feature_set containing byion: matched peaks must be more intense", false);
      setMinFloat_("stage2_byion_intensity_min", 0.0);
      registerInputFile_("stage2_search_features", "<file>", "", "With -stage2: the OpenDIA -features_out table the windows were drawn from "
                         "(candidate_row = window sample_id): transition support, matched scans, peak width and |library RT error| per window", false);
      setValidFormats_("stage2_search_features", {"tsv"});
      registerStringOption_("cnn_arch", "<name>", "grid2d", "Experimental network architecture. 'grid2d' = original 2D CNN over the fragment RT x IM grid. "
                            "'topaz' = redeem TOPAZ-style scorer: 1D conv branches over the fragment XIC, the MS1 M+0..M+2 traces and the fragment "
                            "mobilogram with co-elution heads; MS1 and mobility branches are masked out when absent (reads MS1; uses -stage2_ms1_ppm)", false);
      setValidStrings_("cnn_arch", {"grid2d", "topaz"});
      registerStringOption_("cnn_training", "<mode>", "seeds", "Crossfit only, experimental. 'seeds' = train on the -seeds windows. 'mil' (needs -cnn_arch topaz and "
                            "-rt_calibration, excludes -ss_iterations) = masked-max multiple-instance training on every training-fold precursor: "
                            "its in-RT-window candidates form one bag, labelled target 1 / decoy 0", false);
      setValidStrings_("cnn_training", {"seeds", "mil"});
      registerIntOption_("mil_bags_per_epoch", "<n>", 40000, "With -cnn_training mil: training bags drawn per epoch, also the fixed validation subsample size", false);
      setMinInt_("mil_bags_per_epoch", 1);
      registerIntOption_("topaz_width", "<w>", 32, "With -cnn_arch topaz: base width w (conv w-2w-2w, embeddings 2w, MLP 4w-2w); 32 is the TOPAZ layout", false);
      setMinInt_("topaz_width", 2);
      setMaxInt_("topaz_width", 64);
      registerDoubleOption_("weight_decay", "<value>", 0.0, "With -cnn_arch topaz: decoupled (AdamW) weight decay; 0 is plain Adam", false);
      setMinFloat_("weight_decay", 0.0);
      registerFlag_("topaz_faithful", "With -cnn_arch topaz, experimental: redeem-topaz training defaults (per-sample block max normalisation, cosine "
                    "warm-up LR, gradient clipping at 5, 0.1 x pairwise and, with MIL, 0.1 x in-bag ranking losses). Combine with -weight_decay 0.01 for TOPAZ's AdamW default");
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
        const bool paired_tdc = getFlag_("paired_tdc");
        if (paired_tdc && (!crossfit || target_decoy_ratio != 1.0))
        {
          throw std::invalid_argument("paired_tdc requires crossfit mode and target_decoy_ratio 1 (pairing balances the classes).");
        }
        const bool pooled_tdc = getFlag_("pooled_tdc");
        if (pooled_tdc && !paired_tdc)
        {
          throw std::invalid_argument("pooled_tdc requires paired_tdc.");
        }
        const int ss_iterations = getIntOption_("ss_iterations");
        if (ss_iterations > 0 && !crossfit)
        {
          throw std::invalid_argument("ss_iterations requires crossfit mode.");
        }
        const bool stage2 = getFlag_("stage2");
        if (stage2 && !pooled_tdc)
        {
          throw std::invalid_argument("stage2 requires crossfit mode with -paired_tdc and -pooled_tdc (fold LDA scores are pooled after decoy standardisation).");
        }
        if (!stage2 && !getStringOption_("stage2_features_out").empty())
        {
          throw std::invalid_argument("stage2_features_out requires -stage2.");
        }
        const auto stage2_set = OpenDIA::parseStage2FeatureSet(getStringOption_("stage2_feature_set"));
        const bool stage2_needs_search = stage2_set.search;
        if (!stage2 && (getStringOption_("stage2_feature_set") != "extraction" || !getStringOption_("stage2_search_features").empty()))
        {
          throw std::invalid_argument("stage2_feature_set and stage2_search_features require -stage2.");
        }
        if (stage2_needs_search != !getStringOption_("stage2_search_features").empty())
        {
          throw std::invalid_argument("stage2_search_features is required by, and only accepted with, a stage2_feature_set containing search.");
        }
        // Read before extraction so a missing or malformed table fails before any expensive work.
        std::unordered_map<std::string, OpenDIA::SearchPeakFeatures> search_peaks;
        if (stage2_needs_search) { search_peaks = OpenDIA::readSearchPeakFeatures(getStringOption_("stage2_search_features")); }

        const bool rt_window_calibration = getFlag_("rt_window_calibration");
        if (rt_window_calibration && !crossfit)
        {
          throw std::invalid_argument("rt_window_calibration requires crossfit mode.");
        }
        const std::string rt_calibration_path = getStringOption_("rt_calibration");
        if (!rt_calibration_path.empty() && (!crossfit || rt_window_calibration))
        {
          throw std::invalid_argument("rt_calibration requires crossfit mode and excludes rt_window_calibration.");
        }
        // Read before extraction so an invalid file fails before any expensive work.
        std::optional<OpenDIA::CalibrationRT> rt_calibration;
        if (!rt_calibration_path.empty()) { rt_calibration = OpenDIA::readCalibrationRT(rt_calibration_path); }
        const bool rt_filter = rt_window_calibration || rt_calibration.has_value();
        const bool apply_coordinates = getFlag_("apply_calibration_coordinates");
        if (apply_coordinates && !rt_calibration)
        {
          throw std::invalid_argument("apply_calibration_coordinates requires -rt_calibration.");
        }
        const bool topaz = getStringOption_("cnn_arch") == "topaz";
        const bool mil = getStringOption_("cnn_training") == "mil";
        if (mil && (!crossfit || !topaz || !rt_calibration || ss_iterations > 0))
        {
          throw std::invalid_argument("cnn_training mil requires crossfit mode, -cnn_arch topaz and -rt_calibration, and excludes -ss_iterations.");
        }
        if (!topaz && (getIntOption_("topaz_width") != 32 || getDoubleOption_("weight_decay") != 0.0 || getFlag_("topaz_faithful")))
        {
          throw std::invalid_argument("topaz_width, weight_decay and topaz_faithful require -cnn_arch topaz.");
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
        if (apply_coordinates)
        {
          // Library to run coordinates, before extraction; RT stays library RT (the window uses the map).
          if (rt_calibration->has_ms2)
          {
            for (auto& transition : library.transitions) { transition.product_mz *= 1.0 + rt_calibration->ms2_offset_ppm * 1e-6; }
          }
          if (rt_calibration->has_im)
          {
            for (auto& compound : library.compounds)
            {
              if (std::isfinite(compound.drift_time) && compound.drift_time > 0.0)
              {
                compound.drift_time = rt_calibration->im_intercept + rt_calibration->im_slope * compound.drift_time;
              }
            }
          }
          writeLogInfo_("Experimental calibrated coordinates: ms2_offset_ppm=" +
                        (rt_calibration->has_ms2 ? std::to_string(rt_calibration->ms2_offset_ppm) : std::string("unavailable")) +
                        " im=" + (rt_calibration->has_im ? std::to_string(rt_calibration->im_intercept) + "+" + std::to_string(rt_calibration->im_slope) + "x"
                                                          : std::string("unavailable")));
        }
        OpenDIA::TraceCollector collector(std::move(windows), library, geometry.mz_tolerance,
                                         geometry.rt_half_width, geometry.im_half_width, max_samples,
                                         static_cast<std::size_t>(getIntOption_("max_tensor_gib")) << 30);
        collector.setIonMobilityArrayAsVSSC(getFlag_("im_array_as_vssc"));
        if (stage2 || topaz)
        {
          collector.enableEvidence(getDoubleOption_("stage2_ms1_ppm"), getDoubleOption_("stage2_core_half_width"));
          if (stage2 && stage2_set.byion)
          {
            collector.enableByIonCounts(library, getDoubleOption_("stage2_byion_ppm"), getDoubleOption_("stage2_byion_intensity_min"));
          }
        }
        FileHandler reader;
        PeakFileOptions options;
        // MS1 is read only for the second-stage evidence or the TOPAZ MS1 branch; otherwise the collector sees MS2 only, as before.
        options.setMSLevels(stage2 || topaz ? std::vector<Int>{1, 2} : std::vector<Int>{2});
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
        // TOPAZ MS1 branch input: M+0, M+1, M+2 traces over the tensor RT bins.
        std::vector<float> ms1_traces;
        if (topaz)
        {
          ms1_traces.resize(count * 3 * OpenDIA::TraceCollector::RT_BINS);
          for (std::size_t i = 0; i < count; ++i)
          {
            const auto& e = collector.evidence(i);
            float* out = ms1_traces.data() + i * 3 * OpenDIA::TraceCollector::RT_BINS;
            std::copy(e.ms1.begin(), e.ms1.end(), out);
            std::copy(e.ms1_isotope_traces[0].begin(), e.ms1_isotope_traces[0].end(), out + OpenDIA::TraceCollector::RT_BINS);
            std::copy(e.ms1_isotope_traces[1].begin(), e.ms1_isotope_traces[1].end(), out + 2 * OpenDIA::TraceCollector::RT_BINS);
            data.samples[i].has_ms1 = e.ms1_mask != 0;
          }
          data.ms1 = ms1_traces.data();
          data.ms1_count = ms1_traces.size();
        }
        if (mil)
        {
          // Bags: a precursor's candidates inside the identification-free calibrated RT window (CNN-independent).
          std::map<std::string, double> library_rts;
          for (const auto& compound : library.getCompounds()) { library_rts[compound.id] = compound.rt; }
          std::size_t bagged = 0;
          for (std::size_t i = 0; i < count; ++i)
          {
            const auto found = library_rts.find(collector.window(i).precursor_id);
            if (found == library_rts.end() || !std::isfinite(found->second)) { continue; }
            if (std::abs(collector.window(i).apex_rt - rt_calibration->map(found->second)) <= rt_calibration->tolerance_s)
            {
              data.samples[i].bag_id = collector.window(i).precursor_id;
              ++bagged;
            }
          }
          writeLogInfo_("Experimental MIL training: bagged_windows=" + std::to_string(bagged) + "/" + std::to_string(count));
        }

        const auto start = std::chrono::steady_clock::now();
        OpenDIA::CNNResult result;
        OpenDIA::CNNOptions settings;
        settings.seed = getIntOption_("seed");
        settings.threads = getIntOption_("threads");
        settings.epochs = getIntOption_("epochs");
        settings.patience = getIntOption_("patience");
        settings.batch_size = getIntOption_("batch_size");
        settings.architecture = topaz ? OpenDIA::CNNArchitecture::TOPAZ : OpenDIA::CNNArchitecture::GRID2D;
        settings.mil = mil;
        settings.mil_bags_per_epoch = getIntOption_("mil_bags_per_epoch");
        settings.topaz_width = getIntOption_("topaz_width");
        settings.weight_decay = getDoubleOption_("weight_decay");
        settings.topaz_faithful = getFlag_("topaz_faithful");
        std::vector<OpenDIA::CNNQValueInput> crossfit_rows;
        // Pair key: target TRAML_ID, or a decoy TRAML_ID without its DECOY_ prefix.
        std::vector<std::string> pair_keys(count);
        if (crossfit && paired_tdc)
        {
          const auto traml_ids = TransitionPQPFile().getPQPCurrentIDToTraMLIDMap(getStringOption_("lib").c_str(), "PRECURSOR");
          const std::string prefix = "DECOY_";
          for (std::size_t i = 0; i < count; ++i)
          {
            const auto& precursor_id = collector.window(i).precursor_id;
            const auto found = traml_ids.find(precursor_id);
            if (found == traml_ids.end()) { throw std::invalid_argument("Precursor " + precursor_id + " is missing from the PQP PRECURSOR table."); }
            const std::string& name = found->second;
            if (!data.samples[i].is_decoy) { pair_keys[i] = name; }
            else if (name.starts_with(prefix)) { pair_keys[i] = name.substr(prefix.size()); }
          }
        }
        std::vector<bool> in_rt_window;
        std::vector<OpenDIA::CNNQValueResult> evaluated;
        const double ss_fdr = getDoubleOption_("ss_fdr");
        const auto ss_max_seeds = static_cast<std::size_t>(getIntOption_("ss_max_seeds"));
        // Round 0 trains on -seeds. Each score iteration (PyProphet-style semi-supervised learning) relabels from the
        // previous round's held-out scores and retrains every fold; outputs and checkpoints are from the last round.
        for (int iteration = 0; iteration <= (crossfit ? ss_iterations : 0); ++iteration)
        {
        if (iteration > 0)
        {
          // Positives: target precursor winners at q <= ss_fdr. Negatives: every decoy precursor's best window.
          std::vector<std::size_t> positives, negatives;
          for (std::size_t i = 0; i < count; ++i)
          {
            data.samples[i].label = -1;
            if (!evaluated[i].precursor_best) { continue; }
            if (data.samples[i].is_decoy) { negatives.push_back(i); }
            else if (evaluated[i].q_value <= ss_fdr) { positives.push_back(i); }
          }
          // Deterministic, label-blind subset order: FNV-1a64 over precursor and sample identity.
          const auto order = [&](std::size_t i)
          {
            std::uint64_t hash = UINT64_C(14695981039346656037);
            for (unsigned char c : collector.window(i).precursor_id + ":" + data.samples[i].sample_id) { hash = (hash ^ c) * UINT64_C(1099511628211); }
            return std::make_pair(hash, i);
          };
          const auto by_order = [&](std::size_t a, std::size_t b) { return order(a) < order(b); };
          std::sort(positives.begin(), positives.end(), by_order);
          std::sort(negatives.begin(), negatives.end(), by_order);
          const std::size_t n = std::min({positives.size(), negatives.size(), ss_max_seeds});
          if (n < 3) { throw std::runtime_error("Score iteration found fewer than three positive or negative seeds."); }
          // Every fold trains on the other folds' families; each needs positives outside its own held-out fold.
          for (int fold = 0; fold < getIntOption_("folds"); ++fold)
          {
            const bool usable = std::any_of(positives.begin(), positives.begin() + static_cast<std::ptrdiff_t>(n), [&](std::size_t i)
            {
              return OpenDIA::traceCNNFold(data.samples[i].group_id, settings.seed, getIntOption_("folds")) != fold;
            });
            if (!usable)
            {
              throw std::runtime_error("Score iteration has no positive seeds outside held-out fold " + std::to_string(fold) +
                                       ", so that fold cannot be retrained; all positives come from one fold (see -pooled_tdc).");
            }
          }
          for (std::size_t j = 0; j < n; ++j)
          {
            data.samples[positives[j]].label = 1;
            data.samples[negatives[j]].label = 0;
          }
          writeLogInfo_("Experimental score iteration=" + std::to_string(iteration) + " positives_available=" + std::to_string(positives.size()) +
                        " negatives_available=" + std::to_string(negatives.size()) + " seeds_per_class=" + std::to_string(n));
        }
        crossfit_rows.clear();
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
              data.samples[i].is_decoy, OpenDIA::traceCNNFold(data.samples[i].group_id, settings.seed, settings.outer_folds), result.logits[i], pair_keys[i]});
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
        // Optional RT window: precursor winners and q-values use only in-window rows.
        in_rt_window.assign(crossfit_rows.size(), true);
        if (rt_filter)
        {
          std::map<std::string, double> library_rts;
          for (const auto& compound : library.getCompounds()) { library_rts[compound.id] = compound.rt; }
          std::vector<double> apex_rt, library_rt;
          for (std::size_t i = 0; i < crossfit_rows.size(); ++i)
          {
            const auto found = library_rts.find(crossfit_rows[i].precursor_id);
            if (found == library_rts.end() || !std::isfinite(found->second))
            {
              throw std::invalid_argument("The RT window needs a finite library RT for precursor " + crossfit_rows[i].precursor_id);
            }
            apex_rt.push_back(collector.window(i).apex_rt);
            library_rt.push_back(found->second);
          }
          if (rt_calibration)
          {
            std::set<std::string> precursors, kept;
            for (std::size_t i = 0; i < crossfit_rows.size(); ++i)
            {
              in_rt_window[i] = std::abs(apex_rt[i] - rt_calibration->map(library_rt[i])) <= rt_calibration->tolerance_s;
              precursors.insert(crossfit_rows[i].precursor_id);
              if (in_rt_window[i]) { kept.insert(crossfit_rows[i].precursor_id); }
            }
            writeLogInfo_("Experimental identification-free RT calibration: anchors=" + std::to_string(rt_calibration->map.anchors) +
                          " sigma_s=" + std::to_string(rt_calibration->map.sigma) + " half_window_s=" + std::to_string(rt_calibration->tolerance_s) +
                          " precursors_in_window=" + std::to_string(kept.size()) + "/" + std::to_string(precursors.size()));
          }
          std::vector<OpenDIA::RTWindowFold> report;
          if (rt_window_calibration) { in_rt_window = OpenDIA::crossfitRTWindowMask(crossfit_rows, apex_rt, library_rt, report); }
          for (const auto& fold : report)
          {
            writeLogInfo_("Experimental RT window fold=" + std::to_string(fold.fold) + " anchors=" + std::to_string(fold.map.anchors) +
                          " sigma_s=" + std::to_string(fold.map.sigma) + " half_window_s=" + std::to_string(fold.half_window) +
                          " precursors_in_window=" + std::to_string(fold.kept_precursors) + "/" + std::to_string(fold.precursors));
          }
        }
        std::vector<OpenDIA::CNNQValueInput> kept_rows;
        std::vector<std::size_t> kept_index(crossfit_rows.size(), SIZE_MAX);
        for (std::size_t i = 0; i < crossfit_rows.size(); ++i)
        {
          if (!in_rt_window[i]) { continue; }
          kept_index[i] = kept_rows.size();
          kept_rows.push_back(crossfit_rows[i]);
        }
        const auto kept_evaluated = !crossfit ? std::vector<OpenDIA::CNNQValueResult>{} :
          paired_tdc ? OpenDIA::aggregatePairedCNNQValues(kept_rows, pooled_tdc) : OpenDIA::aggregateCNNQValues(kept_rows, target_decoy_ratio);
        evaluated.clear();
        for (std::size_t i = 0; i < crossfit_rows.size(); ++i)
        {
          evaluated.push_back(kept_index[i] == SIZE_MAX ? OpenDIA::CNNQValueResult{crossfit_rows[i], false, 1.0} : kept_evaluated[kept_index[i]]);
        }
        if (ss_iterations > 0)
        {
          std::size_t targets = 0;
          for (const auto& row : evaluated) { targets += row.precursor_best && !row.input.is_decoy && row.q_value <= 0.01; }
          writeLogInfo_("Experimental score iteration=" + std::to_string(iteration) + " target_precursors_q01=" + std::to_string(targets));
        }
        }
        std::vector<double> stage2_score;
        if (stage2)
        {
          // CNN logit standardised per fold by that fold's decoy precursor winners, as in -pooled_tdc.
          std::map<int, std::vector<double>> decoy_logits;
          for (const auto& row : evaluated)
          {
            if (row.precursor_best && row.input.is_decoy) { decoy_logits[row.input.fold].push_back(row.input.score); }
          }
          std::map<int, std::pair<double, double>> scale;
          for (auto& [fold, values] : decoy_logits)
          {
            const auto median_of = [](std::vector<double> v)
            {
              std::sort(v.begin(), v.end());
              const std::size_t n = v.size();
              return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
            };
            const double centre = median_of(values);
            std::vector<double> deviation;
            for (double v : values) { deviation.push_back(std::abs(v - centre)); }
            const double spread = 1.4826 * median_of(deviation);
            if (!(spread > 0.0)) { throw std::invalid_argument("stage2: zero decoy logit spread in fold " + std::to_string(fold)); }
            scale[fold] = {centre, spread};
          }
          // Run mass offsets, identification-free: intensity-weighted median signed error of the positive seed windows.
          const auto seed_offset = [&](bool ms1)
          {
            std::vector<std::pair<double, double>> values;
            for (std::size_t i = 0; i < count; ++i)
            {
              if (collector.window(i).label != 1) { continue; }
              const auto& e = collector.evidence(i);
              const double w = ms1 ? e.ms1_weight : e.ms2_weight;
              if (w > 0.0) { values.emplace_back((ms1 ? e.ms1_ppm_sum : e.ms2_ppm_sum) / w, w); }
            }
            if (values.empty()) { return 0.0; }
            std::sort(values.begin(), values.end());
            double total = 0.0;
            for (const auto& v : values) { total += v.second; }
            double cumulative = 0.0;
            for (const auto& v : values)
            {
              cumulative += v.second;
              if (cumulative >= 0.5 * total) { return v.first; }
            }
            return values.back().first;
          };
          OpenDIA::Stage2ExtractionOptions extraction;
          extraction.rt_half_width = geometry.rt_half_width;
          extraction.core_half_width = getDoubleOption_("stage2_core_half_width");
          extraction.ms1_offset_ppm = seed_offset(true);
          extraction.ms2_offset_ppm = seed_offset(false);
          std::vector<OpenDIA::Stage2Row> rows;
          std::vector<std::size_t> row_window;
          for (std::size_t i = 0; i < count; ++i)
          {
            if (!in_rt_window[i]) { continue; }
            OpenDIA::Stage2Row row;
            row.precursor_id = crossfit_rows[i].precursor_id;
            row.is_decoy = crossfit_rows[i].is_decoy;
            row.fold = crossfit_rows[i].fold;
            const auto& s = scale.at(row.fold);
            const double cnn_z = (crossfit_rows[i].score - s.first) / s.second;
            if (!stage2_set.extraction)
            {
              row.features.push_back(cnn_z);
            }
            else
            {
              const auto extracted = OpenDIA::stage2Features(collector.intensities().data() + i * OpenDIA::TraceCollector::SAMPLE_SIZE,
                                                             collector.libraryIntensities(i), collector.fragmentMZ(i), collector.evidence(i), extraction);
              row.features.assign(extracted.begin(), extracted.end());
              row.features[0] = cnn_z;
            }
            if (stage2_needs_search)
            {
              // Search peak-group evidence of the candidate this window was drawn from; identity and apex must agree.
              const auto& window = collector.window(i);
              const auto found = search_peaks.find(window.sample_id);
              if (found == search_peaks.end() || found->second.precursor_id != window.precursor_id ||
                  std::abs(found->second.rt_apex - window.apex_rt) > 1e-3)
              {
                throw std::invalid_argument("stage2_search_features has no matching candidate_row, precursor and apex for window " + window.sample_id);
              }
              const auto search = OpenDIA::stage2SearchFeatures(found->second);
              row.features.insert(row.features.end(), search.begin(), search.end());
            }
            if (stage2_set.isotope)
            {
              const auto isotope = OpenDIA::stage2IsotopeFeatures(collector.evidence(i), collector.precursorMZ(i), collector.precursorCharge(i), extraction);
              row.features.insert(row.features.end(), isotope.begin(), isotope.end());
            }
            if (stage2_set.boundary)
            {
              const auto boundary = OpenDIA::stage2BoundaryFeatures(collector.intensities().data() + i * OpenDIA::TraceCollector::SAMPLE_SIZE,
                                                                    collector.libraryIntensities(i), collector.fragmentMZ(i), collector.evidence(i), extraction);
              row.features.insert(row.features.end(), boundary.begin(), boundary.end());
            }
            if (stage2_set.byion)
            {
              const auto byion = OpenDIA::stage2ByIonFeatures(collector.evidence(i));
              row.features.insert(row.features.end(), byion.begin(), byion.end());
            }
            if (stage2_set.shape)
            {
              const auto shape = OpenDIA::stage2ShapeFeatures(collector.intensities().data() + i * OpenDIA::TraceCollector::SAMPLE_SIZE,
                                                              collector.libraryIntensities(i), collector.fragmentMZ(i), collector.evidence(i));
              row.features.insert(row.features.end(), shape.begin(), shape.end());
            }
            rows.push_back(row);
            row_window.push_back(i);
          }
          OpenDIA::Stage2Options stage2_options;
          stage2_options.iterations = getIntOption_("stage2_iterations");
          stage2_options.initial_fdr = getDoubleOption_("stage2_initial_fdr");
          stage2_options.iteration_fdr = getDoubleOption_("stage2_iteration_fdr");
          std::vector<OpenDIA::Stage2FoldReport> report;
          const auto lda = OpenDIA::crossfitStage2LDA(rows, getIntOption_("folds"), stage2_options, report);
          const auto names = OpenDIA::stage2FeatureSetNames(stage2_set);
          writeLogInfo_("Experimental second stage: ms1_offset_ppm=" + std::to_string(extraction.ms1_offset_ppm) +
                        " ms2_offset_ppm=" + std::to_string(extraction.ms2_offset_ppm) + " windows=" + std::to_string(rows.size()));
          for (const auto& fold : report)
          {
            std::string weights;
            for (std::size_t p = 0; p < names.size(); ++p) { weights += " " + names[p] + "=" + std::to_string(fold.weights[p]); }
            writeLogInfo_("Experimental second stage fold=" + std::to_string(fold.fold) + " iterations=" + std::to_string(fold.iterations) +
                          " positives=" + std::to_string(fold.positives) + " negatives=" + std::to_string(fold.negatives) + weights);
          }
          // Paired, pooled q-values from the second-stage score over the in-window rows.
          std::vector<OpenDIA::CNNQValueInput> stage2_rows;
          for (std::size_t r = 0; r < rows.size(); ++r)
          {
            auto input = crossfit_rows[row_window[r]];
            input.score = lda[r];
            stage2_rows.push_back(std::move(input));
          }
          const auto stage2_evaluated = OpenDIA::aggregatePairedCNNQValues(stage2_rows, true);
          stage2_score.assign(count, std::numeric_limits<double>::quiet_NaN());
          for (std::size_t i = 0; i < count; ++i) { evaluated[i] = OpenDIA::CNNQValueResult{crossfit_rows[i], false, 1.0}; }
          for (std::size_t r = 0; r < rows.size(); ++r)
          {
            const std::size_t i = row_window[r];
            stage2_score[i] = lda[r];
            evaluated[i] = stage2_evaluated[r];
            evaluated[i].input.score = crossfit_rows[i].score;  // the logit column keeps the CNN logit
          }
          const std::string features_path = getStringOption_("stage2_features_out");
          if (!features_path.empty())
          {
            std::ofstream features(features_path);
            features << "sample_id\tprecursor_id\tis_decoy\tfold";
            for (const auto& name : names) { features << '\t' << name; }
            features << "\tstage2_score\n" << std::setprecision(8);
            for (std::size_t r = 0; r < rows.size(); ++r)
            {
              const auto& input = crossfit_rows[row_window[r]];
              features << input.sample_id << '\t' << input.precursor_id << '\t' << input.is_decoy << '\t' << input.fold;
              for (double v : rows[r].features) { features << '\t' << v; }
              features << '\t' << lda[r] << '\n';
            }
            features.close();
            if (!features) { throw std::runtime_error("Could not write -stage2_features_out."); }
          }
        }
        std::ofstream scores(output);
        if (crossfit)
        {
          scores << "sample_id\tprecursor_id\tgroup_id\tis_decoy\tfold\tlogit\tprecursor_best\tq_value"
                 << (rt_filter ? "\tin_rt_window" : "") << (stage2 ? "\tstage2_score\n" : "\n")
                 << std::setprecision(std::numeric_limits<double>::max_digits10);
          std::size_t targets_at_one_percent = 0;
          for (std::size_t i = 0; i < evaluated.size(); ++i)
          {
            const auto& row = evaluated[i];
            scores << row.input.sample_id << '\t' << row.input.precursor_id << '\t' << row.input.group_id << '\t'
                   << row.input.is_decoy << '\t' << row.input.fold << '\t' << row.input.score << '\t' << row.precursor_best << '\t';
            if (row.precursor_best) { scores << row.q_value; }
            else { scores << "NA"; }
            if (rt_filter) { scores << '\t' << in_rt_window[i]; }
            if (stage2)
            {
              scores << '\t';
              if (std::isfinite(stage2_score[i])) { scores << stage2_score[i]; }
              else { scores << "NA"; }
            }
            scores << '\n';
            targets_at_one_percent += row.precursor_best && !row.input.is_decoy && row.q_value <= 0.01;
          }
          writeLogInfo_(std::string(stage2 ? "Experimental second-stage pooled paired target-decoy estimate (cross-fitted LDA over CNN logit and extraction features):"
                                    : pooled_tdc ? "Experimental pooled paired target-decoy estimate (per-fold decoy-standardised logits):"
                                    : paired_tdc ? "Experimental fold-wise paired target-decoy estimate:" : "Experimental fold-wise decoy-count estimate:") +
                        " target_decoy_ratio=" + std::to_string(target_decoy_ratio) +
                        " target_precursors_q01=" + std::to_string(targets_at_one_percent) +
                        ". Seed-derived, not independently validated FDR control.");
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
