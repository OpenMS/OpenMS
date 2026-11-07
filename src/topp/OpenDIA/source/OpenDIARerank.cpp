// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIARerank.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <istream>
#include <limits>
#include <map>
#include <ostream>
#include <unordered_map>
#include <utility>

namespace OpenMS::OpenDIA
{
  namespace
  {
    using AssayKey = std::pair<std::string, double>;

    std::vector<std::string> splitTabs_(const std::string& line)
    {
      std::vector<std::string> fields;
      std::size_t start = 0;
      while (true)
      {
        const std::size_t end = line.find('\t', start);
        fields.push_back(line.substr(start, end == std::string::npos ? end : end - start));
        if (end == std::string::npos) break;
        start = end + 1;
      }
      return fields;
    }

    double parseDouble_(const std::string& value, const std::string& column)
    {
      std::size_t parsed = 0;
      const double result = std::stod(value, &parsed);
      if (parsed != value.size() || !std::isfinite(result))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected a finite numeric candidate feature in column '" + column + "'.", value);
      }
      return result;
    }

    Size parseSize_(const std::string& value, const std::string& column)
    {
      if (value.empty() || value.front() == '-')
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected a non-negative integer candidate feature in column '" + column + "'.", value);
      }
      std::size_t parsed = 0;
      const unsigned long long result = std::stoull(value, &parsed);
      if (parsed != value.size() || result > std::numeric_limits<Size>::max())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected a non-negative integer candidate feature in column '" + column + "'.", value);
      }
      return static_cast<Size>(result);
    }

    void validateCandidate_(const CandidateFeatures& candidate)
    {
      if (candidate.precursor_id.empty() || !std::isfinite(candidate.precursor_mz) || !(candidate.precursor_mz > 0.0))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Candidate precursor identity and m/z must be present and valid.",
                                      candidate.precursor_id);
      }
      if (!std::isfinite(candidate.area) || candidate.area < 0.0 ||
          !std::isfinite(candidate.apex_intensity) || candidate.apex_intensity < 0.0 ||
          !std::isfinite(candidate.peak_width_seconds) || candidate.peak_width_seconds < 0.0 ||
          !std::isfinite(candidate.apex_rt_error_seconds) ||
          !std::isfinite(candidate.mean_abs_im_error_vssc) || candidate.mean_abs_im_error_vssc < 0.0 ||
          (!candidate.has_im_match_evidence && candidate.mean_abs_im_error_vssc != 0.0))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Area, apex intensity, and peak width must be finite and non-negative; "
                                      "RT error and mobility residual must be finite and valid.",
                                      candidate.precursor_id);
      }
    }

    bool parseBool_(const std::string& value, const std::string& column)
    {
      if (value == "0") return false;
      if (value == "1") return true;
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Expected 0 or 1 in candidate feature column '" + column + "'.", value);
    }

    std::vector<double> toModelFeatures_(const CandidateFeatures& candidate)
    {
      validateCandidate_(candidate);
      return {
        std::log1p(candidate.area),
        std::log1p(candidate.apex_intensity),
        std::log1p(candidate.peak_width_seconds),
        static_cast<double>(candidate.transition_support),
        std::log1p(static_cast<double>(candidate.matched_scans)),
        candidate.has_library_rt ? 1.0 : 0.0,
        candidate.has_library_rt ? candidate.apex_rt_error_seconds : 0.0,
        candidate.has_library_im ? 1.0 : 0.0,
        candidate.has_im_match_evidence ? 1.0 : 0.0,
        candidate.has_im_match_evidence ? candidate.mean_abs_im_error_vssc : 0.0
      };
    }

    std::map<AssayKey, std::vector<std::size_t>> groupAssays_(const std::vector<CandidateFeatures>& rows)
    {
      std::map<AssayKey, std::vector<std::size_t>> groups;
      std::map<AssayKey, bool> assay_labels;
      for (std::size_t i = 0; i < rows.size(); ++i)
      {
        const auto key = AssayKey{rows[i].precursor_id, rows[i].precursor_mz};
        const auto [label_it, inserted] = assay_labels.emplace(key, rows[i].is_decoy);
        if (!inserted && label_it->second != rows[i].is_decoy)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Candidate windows for one assay have inconsistent target/decoy labels.",
                                        rows[i].precursor_id);
        }
        groups[key].push_back(i);
      }
      return groups;
    }

    std::vector<std::size_t> bestByArea_(const std::vector<CandidateFeatures>& rows,
                                          const std::map<AssayKey, std::vector<std::size_t>>& groups)
    {
      std::vector<std::size_t> best;
      best.reserve(groups.size());
      for (const auto& [key, indexes] : groups)
      {
        (void)key;
        std::size_t best_index = indexes.front();
        for (const std::size_t index : indexes)
        {
          if (rows[index].area > rows[best_index].area) best_index = index;
        }
        best.push_back(best_index);
      }
      return best;
    }

    void assignAssayWinnersAndQValues_(std::vector<ScoredCandidate>& rows,
                                       const std::map<AssayKey, std::vector<std::size_t>>& groups)
    {
      std::vector<std::size_t> winners;
      winners.reserve(groups.size());
      for (const auto& [key, indexes] : groups)
      {
        (void)key;
        std::size_t winner = indexes.front();
        for (const std::size_t index : indexes)
        {
          if (rows[index].score > rows[winner].score) winner = index;
        }
        rows[winner].assay_best = true;
        winners.push_back(winner);
      }

      std::size_t target_total = 0;
      std::size_t decoy_total = 0;
      for (const std::size_t index : winners)
      {
        rows[index].candidate.is_decoy ? ++decoy_total : ++target_total;
      }
      if (target_total == 0 || decoy_total == 0)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Held-out assay set must contain at least one target and one decoy.",
                                      "target/decoy counts");
      }

      std::stable_sort(winners.begin(), winners.end(), [&rows](std::size_t lhs, std::size_t rhs)
      {
        return rows[lhs].score > rows[rhs].score;
      });

      struct ScoreGroup
      {
        std::size_t begin = 0;
        std::size_t end = 0;
        double fdr = 1.0;
        double q_value = 1.0;
      };
      std::vector<ScoreGroup> score_groups;
      std::size_t targets = 0;
      std::size_t decoys = 0;
      for (std::size_t begin = 0; begin < winners.size();)
      {
        std::size_t end = begin + 1;
        const double score = rows[winners[begin]].score;
        while (end < winners.size() && rows[winners[end]].score == score) ++end;
        for (std::size_t i = begin; i < end; ++i)
        {
          rows[winners[i]].candidate.is_decoy ? ++decoys : ++targets;
        }
        const double fdr = targets == 0 ? 1.0 : std::min(1.0,
          ((static_cast<double>(decoys) + 1.0) / static_cast<double>(targets)) *
          (static_cast<double>(target_total) / static_cast<double>(decoy_total)));
        score_groups.push_back({begin, end, fdr, 1.0});
        begin = end;
      }
      double q_value = 1.0;
      for (auto it = score_groups.rbegin(); it != score_groups.rend(); ++it)
      {
        q_value = std::min(q_value, it->fdr);
        it->q_value = q_value;
        for (std::size_t i = it->begin; i < it->end; ++i)
        {
          rows[winners[i]].assay_q_value = q_value;
        }
      }
    }
  }

  std::vector<CandidateFeatures> readCandidateFeatures(std::istream& input)
  {
    std::string line;
    if (!std::getline(input, line))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Candidate feature table is empty.", "header");
    }
    const auto header = splitTabs_(line);
    std::unordered_map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header.size(); ++i) columns.emplace(header[i], i);
    const std::vector<std::string> required = {
      "candidate_row", "precursor_id", "precursor_mz", "is_decoy", "peak_width_seconds", "area",
      "apex_intensity", "transition_support", "matched_scans", "has_library_rt",
      "apex_rt_error_seconds", "has_library_im"
    };
    for (const auto& name : required)
    {
      if (columns.find(name) == columns.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Required candidate feature column is missing: " + name, name);
      }
    }
    const bool has_im_evidence_column = columns.find("has_im_match_evidence") != columns.end();
    const bool has_im_error_column = columns.find("mean_abs_im_error_vssc") != columns.end();
    if (has_im_evidence_column != has_im_error_column)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Candidate mobility evidence and residual columns must be present together.",
                                    "mobility feature columns");
    }

    std::vector<CandidateFeatures> candidates;
    while (std::getline(input, line))
    {
      if (line.empty()) continue;
      const auto fields = splitTabs_(line);
      if (fields.size() != header.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Candidate feature row has a different number of fields than the header.", line);
      }
      const auto get = [&fields, &columns](const std::string& name) -> const std::string&
      {
        return fields[columns.at(name)];
      };
      CandidateFeatures candidate;
      candidate.row = parseSize_(get("candidate_row"), "candidate_row");
      candidate.precursor_id = get("precursor_id");
      candidate.precursor_mz = parseDouble_(get("precursor_mz"), "precursor_mz");
      candidate.is_decoy = parseBool_(get("is_decoy"), "is_decoy");
      candidate.peak_width_seconds = parseDouble_(get("peak_width_seconds"), "peak_width_seconds");
      candidate.area = parseDouble_(get("area"), "area");
      candidate.apex_intensity = parseDouble_(get("apex_intensity"), "apex_intensity");
      candidate.transition_support = parseSize_(get("transition_support"), "transition_support");
      candidate.matched_scans = parseSize_(get("matched_scans"), "matched_scans");
      candidate.has_library_rt = parseBool_(get("has_library_rt"), "has_library_rt");
      candidate.apex_rt_error_seconds = parseDouble_(get("apex_rt_error_seconds"), "apex_rt_error_seconds");
      candidate.has_library_im = parseBool_(get("has_library_im"), "has_library_im");
      if (has_im_evidence_column)
      {
        candidate.has_im_match_evidence = parseBool_(get("has_im_match_evidence"), "has_im_match_evidence");
        candidate.mean_abs_im_error_vssc = parseDouble_(get("mean_abs_im_error_vssc"), "mean_abs_im_error_vssc");
      }
      validateCandidate_(candidate);
      candidates.push_back(std::move(candidate));
    }
    if (candidates.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Candidate feature table contains no candidate rows.", "rows");
    }
    return candidates;
  }

  std::vector<ScoredCandidate> rescoreCandidates(const std::vector<CandidateFeatures>& training,
                                                  const std::vector<CandidateFeatures>& scoring,
                                                  int seed,
                                                  int threads)
  {
    if (training.empty() || scoring.empty() || threads < 1)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Training and scoring rows must be non-empty and threads must be positive.",
                                    "candidate rows");
    }
    if (seed < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Percolator seed must be non-negative.", "seed");
    }

    const auto training_groups = groupAssays_(training);
    const auto scoring_groups = groupAssays_(scoring);
    const auto training_best = bestByArea_(training, training_groups);
    RescoreInput train_input;
    train_input.feature_names = StringList{
      "log_area", "log_apex_intensity", "log_peak_width", "transition_support", "log_matched_scans",
      "has_library_rt", "apex_rt_error_seconds", "has_library_im", "has_im_match_evidence", "mean_abs_im_error_vssc"
    };
    train_input.features.reserve(training_best.size());
    train_input.is_decoy.reserve(training_best.size());
    for (const std::size_t index : training_best)
    {
      train_input.features.push_back(toModelFeatures_(training[index]));
      train_input.is_decoy.push_back(training[index].is_decoy);
    }

    RescoreInput score_input;
    score_input.feature_names = train_input.feature_names;
    score_input.features.reserve(scoring.size());
    score_input.is_decoy.reserve(scoring.size());
    for (const auto& candidate : scoring)
    {
      score_input.features.push_back(toModelFeatures_(candidate));
      score_input.is_decoy.push_back(candidate.is_decoy);
    }

    Percolator percolator;
    Param parameters = percolator.getDefaults();
    parameters.setValue("seed", seed);
    parameters.setValue("num_threads", threads);
    percolator.setParameters(parameters);
    const PercolatorModel model = percolator.train(train_input);
    const RescoreOutput score_output = percolator.score(score_input, model, train_input);
    if (score_output.scores.size() != scoring.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Percolator returned a score count different from the held-out candidates.",
                                    "score count");
    }

    std::vector<ScoredCandidate> result;
    result.reserve(scoring.size());
    for (std::size_t i = 0; i < scoring.size(); ++i)
    {
      result.push_back({scoring[i], score_output.scores[i], false, 1.0});
    }
    assignAssayWinnersAndQValues_(result, scoring_groups);
    return result;
  }

  void writeScoredCandidates(std::ostream& output, const std::vector<ScoredCandidate>& candidates)
  {
    output << "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\tpercolator_score\tassay_best\tassay_q_value\n";
    output << std::setprecision(17);
    for (const auto& row : candidates)
    {
      output << row.candidate.row << '\t' << row.candidate.precursor_id << '\t' << row.candidate.precursor_mz << '\t'
             << (row.candidate.is_decoy ? 1 : 0) << '\t' << row.score << '\t' << (row.assay_best ? 1 : 0) << '\t';
      if (row.assay_best) output << row.assay_q_value;
      else output << "NA";
      output << '\n';
    }
  }
}
