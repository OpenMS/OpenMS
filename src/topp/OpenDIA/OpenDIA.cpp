// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAWorkflow.h"

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <unistd.h>

using namespace OpenMS;

/**
  @page TOPP_OpenDIA OpenDIA
  @brief Complete OpenDIA analysis of DIA runs: calibration, candidate search, features, scoring and q-values.

  For every input run: identification-free calibration of library RT, fragment m/z and 1/K0; candidate search with the
  calibrated library; one more pass over the spectra for the scoring evidence (6-fragment traces and apex b/y ions,
  12-fragment grids, shared-signal interference); a cross-fitted, self-trained scorer per run; and q-values by pooled
  paired target-decoy competition. Across runs, precursors get a global q-value from their best run score.

  Output directory: <run>.precursors.tsv (one row per precursor winner), <run>.calibration.json, precursors.tsv
  (precursors at the run and global q cut, all runs) and optionally <run>.features.tsv (every scored window with its
  features).
*/
class TOPPOpenDIA final : public TOPPBase
{
public:
  TOPPOpenDIA() : TOPPBase("OpenDIA", "Complete OpenDIA analysis of DIA runs: raw spectra to precursor q-values.", {}, false) {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", {}, "DIA runs (Thermo raw, mzML, ...)");
    registerInputFile_("lib", "<file>", "", "PQP library with decoys named DECOY_<target TRAML_ID> (search fragments, up to 6 per precursor)");
    setValidFormats_("lib", {"pqp"});
    registerInputFile_("scoring_lib", "<file>", "", "Optional PQP with up to 12 scoring fragments per precursor (same precursor IDs); "
                       "default: the search fragments", false);
    setValidFormats_("scoring_lib", {"pqp"});
    registerOutputDir_("out_dir", "<dir>", "", "Output directory");
    registerDoubleOption_("q", "<q>", 0.01, "q-value cut of the combined precursor table", false);
    setMinFloat_("q", 0.0);
    setMaxFloat_("q", 1.0);
    registerFlag_("features_out", "Also write <run>.features.tsv with every scored window and its features");
    registerIntOption_("seed", "<n>", 0, "Scorer seed", false);
    setMinInt_("seed", 0);
    registerFlag_("im_array_as_vssc", "Interpret a spectrum ion-mobility array whose unit is not CCS as 1/K0");
    registerIntOption_("max_peaks", "<n>", 5, "Advanced: candidate peaks kept per precursor by the search, best fragment co-elution first", false, true);
    setMinInt_("max_peaks", 1);
    registerDoubleOption_("valley_fraction", "<fraction>", 0.5, "Advanced: in the search, two maxima stay separate candidate peaks only if the valley "
                          "between them is at most this fraction of the lower one (1 keeps every smoothed maximum)", false, true);
    setMinFloat_("valley_fraction", 0.0);
    setMaxFloat_("valley_fraction", 1.0);
    registerDoubleOption_("im_tolerance_scale", "<factor>", 0.6, "Advanced: multiply the calibrated 1/K0 half width used for search and extraction "
                          "(PLAN_IMNARROW: 0.6; 1 = as calibrated)", false, true);
    setMinFloat_("im_tolerance_scale", 0.05);
    setMaxFloat_("im_tolerance_scale", 4.0);
    registerDoubleOption_("trace_tolerance_ppm", "<ppm>", 10.0, "Advanced: fragment m/z half width of the 6-fragment trace tensor in ppm "
                          "(PLAN_TRPPM: 10, +18 % IDs at 1 % true FDP on the r ~ 7 HeLa entrapment)", false, true);
    setMinFloat_("trace_tolerance_ppm", 0.0);
    registerIntOption_("trace_tensor_gb", "<GiB>", 64, "Advanced: hard cap on the 6-fragment trace tensor in GiB (about 6.1 KiB per candidate window with ion mobility); "
                       "raise it for very large libraries, memory permitting", false, true);
    setMinInt_("trace_tensor_gb", 1);
    registerFlag_("no_rt_refine", "Advanced: skip the per-run RT refinement (PLAN_RTR): by default a ridge RT correction on sequence composition is fitted from the run's own q<=0.01 IDs after a first pass, then search, extraction and scoring run again with corrected RTs and a narrower RT window (about 1.8-1.9x wall time)", true);
    registerIntOption_("rt_refine_min_ids", "<n>", 500, "Advanced: skip the RT refinement with fewer training IDs", false, true);
    setMinInt_("rt_refine_min_ids", 1);
    registerStringOption_("decoys", "<mode>", "library",
                          "library: use the library's decoys. predict: replace them by pseudo-reversed decoys predicted like the "
                          "targets (own RT, ion mobility and fragment intensities, PeptDeep; DESIGN_OPENDIA_DECOYS), written to "
                          "<out_dir>/<lib>.predicted_decoys.pqp. Needs a WITH_ONNX build; for predicted libraries", false);
    setValidStrings_("decoys", {"library", "predict"});
    registerStringOption_("decoy_model_dir", "<dir>", "",
                          "decoys predict: directory with peptdeep_{rt,ms2,ccs}_dynamic.onnx; empty uses the installed models", false, true);
    registerStringOption_("decoy_cache", "<dir>", "",
                          "decoys predict: directory of the predicted library (default: out_dir). A library predicted earlier from the "
                          "same input library, models and OpenMS revision is reused instead of predicted again, so a study pays the "
                          "prediction once even when every run is a separate call", false, true);
  }

  static std::string stem_(const std::string& path)
  {
    return std::filesystem::path(path).stem().string();
  }

  static void number_(std::string& row, double v)
  {
    char buffer[32];
    const auto r = std::to_chars(buffer, buffer + sizeof(buffer), v);
    row.append(1, '\t').append(buffer, r.ptr);
  }

  /// -decoys predict: drop the library's decoys and give every target a pseudo-reversed decoy predicted like a target
  /// (OpenSwathLibraryPreparation::generateDecoys with predict_decoys and predict_keep_short, settings estimated from the
  /// library). Returns the path of the target/decoy PQP it writes to @p cache_dir. The library is predicted once per
  /// input: <stem>.predicted_decoys.pqp.key holds the SHA-1 of the input library, the model directory and the OpenMS
  /// revision, and a library with a matching key is reused. Prediction writes to a temporary name and renames, so a
  /// concurrent or interrupted call never leaves a partial library under the final name.
  /// ponytail: the only remaining OpenSWATH dependency of OpenDIA; to be removed or replaced (decision pending).
  std::string predictDecoys_(const std::string& lib, const std::string& cache_dir)
  {
    const std::string stem = stem_(lib);
    const std::string out = cache_dir + "/" + stem + ".predicted_decoys.pqp";
    const std::string key_file = out + ".key";
    const std::string key = "lib_sha1=" + FileHandler::computeFileHash(lib) + "\nmodel_dir=" + getStringOption_("decoy_model_dir") +
                            "\nopenms=" + VersionInfo::getVersion() + " " + VersionInfo::getRevision() + "\n";
    {
      std::ifstream in(key_file);
      const std::string cached((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      if (cached == key && std::filesystem::exists(out))
      {
        writeLogInfo_("Decoy prediction: reusing " + out + " (same input library, models and OpenMS revision)");
        return out;
      }
    }
    std::filesystem::create_directories(cache_dir);
    OpenSwath::LightTargetedExperiment all;
    // The written library has canonical precursor IDs (the output's precursor_id refers to it); sequence, charge and
    // proteins are those of the input.
    TransitionPQPFile().convertPQPToTargetedExperiment(lib.c_str(), all);
    OpenSwath::LightTargetedExperiment targets;
    std::set<std::string> decoy_ids;
    for (const auto& t : all.transitions)
    {
      if (t.getDecoy()) decoy_ids.insert(t.peptide_ref);
    }
    for (auto& c : all.compounds)
    {
      if (!decoy_ids.contains(c.id)) targets.compounds.push_back(std::move(c));
    }
    for (auto& t : all.transitions)
    {
      if (!t.getDecoy()) targets.transitions.push_back(std::move(t));
    }
    std::set<std::string> used;
    for (const auto& c : targets.compounds) used.insert(c.protein_refs.begin(), c.protein_refs.end());
    for (auto& p : all.proteins)
    {
      if (used.contains(p.id)) targets.proteins.push_back(std::move(p));
    }
    all = OpenSwath::LightTargetedExperiment();
    writeLogInfo_("Predicting decoys for " + std::to_string(targets.compounds.size()) + " target precursors (" +
                  std::to_string(decoy_ids.size()) + " library decoys dropped)");

    const std::string tag = "." + std::to_string(::getpid()) + ".tmp.pqp";
    const std::string targets_pqp = cache_dir + "/" + stem + ".targets" + tag;
    const std::string partial = cache_dir + "/" + stem + ".predicted_decoys" + tag;
    std::filesystem::remove(targets_pqp);
    std::filesystem::remove(partial);
    std::filesystem::remove(key_file);
    // Remove the temporaries also when prediction throws (e.g. too few annotated targets).
    struct Cleanup
    {
      std::vector<std::string> paths;
      ~Cleanup()
      {
        std::error_code ignored;
        for (const auto& p : paths) std::filesystem::remove(p, ignored);
      }
    } cleanup{{targets_pqp, partial}};
    TransitionPQPFile().convertLightTargetedExperimentToPQP(targets_pqp.c_str(), targets);
    targets = OpenSwath::LightTargetedExperiment();

    OpenSwathLibraryPreparation::DecoyGeneratorParameters parameters;
    parameters.method = "pseudo-reverse";
    parameters.predict_decoys = true;
    parameters.predict_keep_short = true;
    parameters.predict_model_directory = getStringOption_("decoy_model_dir");
    parameters.predict_threads = std::max(1, getIntOption_("threads"));
    OpenSwathLibraryPreparation preparation;
    preparation.setLogType(log_type_);
    preparation.generateDecoys(targets_pqp, FileTypes::PQP, partial, FileTypes::PQP, parameters);
    std::filesystem::rename(partial, out);
    std::ofstream(key_file) << key;
    return out;
  }

  ExitCodes main_(int, const char**) override
  {
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<std::string> inputs = getStringList_("in");
    const std::string out_dir = getOutputDirOption("out_dir");
    const double q_cut = getDoubleOption_("q");
    if (inputs.empty())
    {
      writeLogError_("No input runs.");
      return ILLEGAL_PARAMETERS;
    }
    std::filesystem::create_directories(out_dir);
    std::map<std::string, int> stems;
    for (const auto& in : inputs)
    {
      if (++stems[stem_(in)] > 1)
      {
        writeLogError_("Two input runs share the file name " + stem_(in) + "; outputs would collide.");
        return ILLEGAL_PARAMETERS;
      }
    }

    std::string lib_path = getStringOption_("lib");
    if (getStringOption_("decoys") == "predict")
    {
      if (!getStringOption_("scoring_lib").empty())
      {
        writeLogError_("-decoys predict does not support -scoring_lib (its decoys would not match the predicted ones).");
        return ILLEGAL_PARAMETERS;
      }
      const std::string cache = getStringOption_("decoy_cache");
      lib_path = predictDecoys_(lib_path, cache.empty() ? out_dir : cache);
    }
    writeLogInfo_("Reading library " + lib_path);
    OpenDIA::AnalysisLibrary library = OpenDIA::AnalysisLibrary::read(lib_path, getStringOption_("scoring_lib"));
    std::size_t library_rss = 0;
    SysInfo::getProcessMemoryConsumption(library_rss);
    writeLogInfo_("Library: precursors=" + std::to_string(library.precursors.size()) + " pairs=" + std::to_string(library.pairs) +
                  " rss_mb=" + std::to_string(library_rss / 1024));

    OpenDIA::RunOptions options;
    options.threads = std::max(1, getIntOption_("threads"));
    options.seed = static_cast<std::uint64_t>(getIntOption_("seed"));
    options.im_array_as_vssc = getFlag_("im_array_as_vssc");
    options.max_peaks = static_cast<std::size_t>(getIntOption_("max_peaks"));
    options.search_valley_fraction = getDoubleOption_("valley_fraction");
    options.im_tolerance_scale = getDoubleOption_("im_tolerance_scale");
    options.trace_tolerance_ppm = getDoubleOption_("trace_tolerance_ppm");
    if (!(options.search_valley_fraction > 0.0))
    {
      writeLogError_("-valley_fraction must be in (0, 1].");
      return ILLEGAL_PARAMETERS;
    }
    if (!(options.trace_tolerance_ppm > 0.0))
    {
      writeLogError_("-trace_tolerance_ppm must be positive.");
      return ILLEGAL_PARAMETERS;
    }
    options.trace_tensor_gb = static_cast<std::size_t>(getIntOption_("trace_tensor_gb"));
    options.rt_refine = !getFlag_("no_rt_refine");
    options.rt_refine_min_ids = static_cast<std::size_t>(getIntOption_("rt_refine_min_ids"));
    options.calibration.panel_targets = 40000;

    // Best (calibrated) score of each precursor per run, for the global q-value.
    struct Hit
    {
      std::size_t run;
      std::uint32_t window;
      float score;
      double q;
    };
    std::vector<std::vector<Hit>> hits(library.precursors.size());
    std::vector<OpenDIA::RunResult> results;
    for (std::size_t r = 0; r < inputs.size(); ++r)
    {
      const std::string name = stem_(inputs[r]);
      writeLogInfo_("Run " + std::to_string(r + 1) + "/" + std::to_string(inputs.size()) + ": " + inputs[r]);
      OpenDIA::RunAnalysis analysis(library, options, [&](const std::string& line) { writeLogInfo_("  " + line); });
      OpenDIA::RunResult result = analysis.run(inputs[r]);
      {
        std::ofstream json(out_dir + "/" + name + ".calibration.json");
        json << result.calibration_json;
      }
      if (!writeRun_(out_dir + "/" + name + ".precursors.tsv", name, library, result))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      if (getFlag_("features_out") && !writeFeatures_(out_dir + "/" + name + ".features.tsv", library, result))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      for (std::uint32_t i = 0; i < result.windows.size(); ++i)
      {
        const auto& w = result.windows[i];
        if (w.score > OpenDIA::CrossFitScorer::NOT_KEPT) hits[w.precursor].push_back({r, i, w.score, w.q_value});
      }
      // Keep only what the combined table needs.
      std::vector<float>().swap(result.features);
      results.push_back(std::move(result));
    }

    // Global precursor q: best calibrated score over runs, paired target-decoy competition (ties to the decoy).
    std::vector<double> global_q(library.precursors.size(), 1.0);
    {
      std::vector<double> best(library.precursors.size(), -std::numeric_limits<double>::infinity());
      for (std::size_t p = 0; p < hits.size(); ++p)
        for (const auto& h : hits[p]) best[p] = std::max(best[p], static_cast<double>(h.score));
      std::vector<std::int64_t> pair_winner(library.pairs, -1);
      for (std::uint32_t p = 0; p < library.precursors.size(); ++p)
      {
        if (!std::isfinite(best[p])) continue;
        auto& w = pair_winner[library.precursors[p].pair];
        if (w < 0 || best[p] > best[w] || (best[p] == best[w] && library.precursors[p].decoy)) w = p;
      }
      std::vector<std::uint32_t> ranked;
      for (auto w : pair_winner)
        if (w >= 0) ranked.push_back(static_cast<std::uint32_t>(w));
      std::sort(ranked.begin(), ranked.end(), [&](std::uint32_t a, std::uint32_t b)
      {
        if (best[a] != best[b]) return best[a] > best[b];
        if (library.precursors[a].decoy != library.precursors[b].decoy) return library.precursors[a].decoy;
        return a < b;
      });
      std::vector<double> q(ranked.size());
      double t = 0.0, d = 0.0;
      for (std::size_t i = 0; i < ranked.size(); ++i)
      {
        (library.precursors[ranked[i]].decoy ? d : t) += 1.0;
        q[i] = std::min(1.0, (d + 1.0) / std::max(t, 1.0));
      }
      for (std::size_t i = q.size(); i-- > 1;) q[i - 1] = std::min(q[i - 1], q[i]);
      for (std::size_t i = 0; i < ranked.size(); ++i) global_q[ranked[i]] = q[i];
    }

    std::ofstream all(out_dir + "/precursors.tsv");
    all << "run\tprecursor_id\tsequence\tcharge\tproteins\trt_apex\trt_start\trt_end\tarea\tscore\tq_value\tglobal_q_value\n";
    std::size_t rows = 0, global_targets = 0;
    for (std::uint32_t p = 0; p < library.precursors.size(); ++p)
    {
      const auto& pre = library.precursors[p];
      if (pre.decoy) continue;
      if (global_q[p] <= q_cut) ++global_targets;
      for (const auto& h : hits[p])
      {
        if (h.q > q_cut || global_q[p] > q_cut) continue;
        const auto& w = results[h.run].windows[h.window];
        std::string proteins;
        for (const auto& s : pre.proteins) proteins += (proteins.empty() ? "" : ";") + s;
        std::string row = stem_(inputs[h.run]) + '\t' + pre.id + '\t' + pre.sequence + '\t' + std::to_string(pre.charge) + '\t' + proteins;
        for (double v : {w.peak.rt_apex, w.peak.rt_start, w.peak.rt_end, w.peak.area, static_cast<double>(w.score), w.q_value, global_q[p]}) number_(row, v);
        all << row << '\n';
        ++rows;
      }
    }
    all.close();
    if (!all) return CANNOT_WRITE_OUTPUT_FILE;
    writeLogInfo_("Combined: " + std::to_string(global_targets) + " target precursors at global q <= " + std::to_string(q_cut) + ", " +
                  std::to_string(rows) + " run-level rows; seconds=" +
                  std::to_string(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()));
    return EXECUTION_OK;
  }

  /// One row per precursor winner (the kept window), with its run q-value.
  static bool writeRun_(const std::string& path, const std::string& run, const OpenDIA::AnalysisLibrary& library, const OpenDIA::RunResult& result)
  {
    std::ofstream out(path);
    out << "run\tcandidate_row\tprecursor_id\tis_decoy\tsequence\tcharge\trt_apex\trt_start\trt_end\tarea\ttransition_support\tscore\tq_value\n";
    for (const auto& w : result.windows)
    {
      if (!(w.score > OpenDIA::CrossFitScorer::NOT_KEPT)) continue;
      const auto& p = library.precursors[w.precursor];
      std::string row = run + '\t' + std::to_string(w.candidate) + '\t' + p.id + '\t' + (p.decoy ? "1" : "0") + '\t' + p.sequence + '\t' + std::to_string(p.charge);
      for (double v : {w.peak.rt_apex, w.peak.rt_start, w.peak.rt_end, w.peak.area}) number_(row, v);
      row += '\t' + std::to_string(w.peak.transition_support);
      for (double v : {static_cast<double>(w.score), w.q_value}) number_(row, v);
      out << row << '\n';
    }
    out.close();
    return static_cast<bool>(out);
  }

  /// Every scored window with its features (candidate_row = search feature table row).
  static bool writeFeatures_(const std::string& path, const OpenDIA::AnalysisLibrary& library, const OpenDIA::RunResult& result)
  {
    std::ofstream out(path);
    out << "candidate_row\tprecursor_id\tis_decoy\tfold\trt_apex";
    for (const auto& n : result.feature_names) out << '\t' << n;
    out << "\tscore\tq_value\n";
    const std::size_t d = result.feature_names.size();
    for (std::size_t i = 0; i < result.windows.size(); ++i)
    {
      const auto& w = result.windows[i];
      const auto& p = library.precursors[w.precursor];
      std::string row = std::to_string(w.candidate) + '\t' + p.id + '\t' + (p.decoy ? "1" : "0") + '\t' + std::to_string(p.fold);
      number_(row, w.peak.rt_apex);
      for (std::size_t k = 0; k < d; ++k) number_(row, result.features[i * d + k]);
      number_(row, w.score);
      number_(row, w.q_value);
      out << row << '\n';
    }
    out.close();
    return static_cast<bool>(out);
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIA tool;
  return tool.main(argc, argv);
}
