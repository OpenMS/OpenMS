// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAWorkflow.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <map>

using namespace OpenMS;

/**
  @page TOPP_OpenDIAWorkflow OpenDIAWorkflow
  @brief Complete OpenDIA analysis of DIA runs: calibration, candidate search, features, scoring and q-values.

  For every input run: identification-free calibration of library RT, fragment m/z and 1/K0; candidate search with the
  calibrated library; one more pass over the spectra for the scoring evidence (6-fragment traces and apex b/y ions,
  12-fragment grids, shared-signal interference); a cross-fitted, self-trained scorer per run; and q-values by pooled
  paired target-decoy competition. Across runs, precursors get a global q-value from their best run score.

  Output directory: <run>.precursors.tsv (one row per precursor winner), <run>.calibration.json, precursors.tsv
  (precursors at the run and global q cut, all runs) and optionally <run>.features.tsv (every scored window with its
  features).
*/
class TOPPOpenDIAWorkflow final : public TOPPBase
{
public:
  TOPPOpenDIAWorkflow() : TOPPBase("OpenDIAWorkflow", "Complete OpenDIA analysis of DIA runs: raw spectra to precursor q-values.", {}, false) {}

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
    registerStringOption_("calibration", "<files>", "", "Advanced: comma-separated OpenDIACalibrate JSON per run to reuse instead of calibrating", false, true);
    registerStringOption_("anchors", "<files>", "", "Advanced: comma-separated anchor TSV per run, with -calibration", false, true);
  }

  static std::string stem_(const std::string& path)
  {
    return std::filesystem::path(path).stem().string();
  }

  static std::vector<std::string> list_(const std::string& s)
  {
    std::vector<std::string> out;
    if (s.empty()) return out;
    StringUtils::split(s, ',', out);
    return out;
  }

  static void number_(std::string& row, double v)
  {
    char buffer[32];
    const auto r = std::to_chars(buffer, buffer + sizeof(buffer), v);
    row.append(1, '\t').append(buffer, r.ptr);
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
    const auto calibrations = list_(getStringOption_("calibration"));
    const auto anchor_files = list_(getStringOption_("anchors"));
    if ((!calibrations.empty() || !anchor_files.empty()) && (calibrations.size() != inputs.size() || anchor_files.size() != inputs.size()))
    {
      writeLogError_("-calibration and -anchors need one file per input run.");
      return ILLEGAL_PARAMETERS;
    }
    std::map<std::string, int> stems;
    for (const auto& in : inputs)
    {
      if (++stems[stem_(in)] > 1)
      {
        writeLogError_("Two input runs share the file name " + stem_(in) + "; outputs would collide.");
        return ILLEGAL_PARAMETERS;
      }
    }

    writeLogInfo_("Reading library " + getStringOption_("lib"));
    OpenDIA::AnalysisLibrary library = OpenDIA::AnalysisLibrary::read(getStringOption_("lib"), getStringOption_("scoring_lib"));
    writeLogInfo_("Library: precursors=" + std::to_string(library.precursors.size()) + " pairs=" + std::to_string(library.pairs));

    OpenDIA::RunOptions options;
    options.threads = std::max(1, getIntOption_("threads"));
    options.seed = static_cast<std::uint64_t>(getIntOption_("seed"));
    options.im_array_as_vssc = getFlag_("im_array_as_vssc");
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
      if (!calibrations.empty())
      {
        options.calibration_json = calibrations[r];
        options.anchors_tsv = anchor_files[r];
      }
      const std::string name = stem_(inputs[r]);
      writeLogInfo_("Run " + std::to_string(r + 1) + "/" + std::to_string(inputs.size()) + ": " + inputs[r]);
      OpenDIA::RunAnalysis analysis(library, options, [&](const std::string& line) { writeLogInfo_("  " + line); });
      OpenDIA::RunResult result = analysis.run(inputs[r]);
      if (!result.calibration_json.empty())
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
  TOPPOpenDIAWorkflow tool;
  return tool.main(argc, argv);
}
