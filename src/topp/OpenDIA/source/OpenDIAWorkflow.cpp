// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAWorkflow.h>

#include <OpenDIAFragmentGrid.h>
#include <OpenDIAInterference.h>
#include <OpenDIAPeakMaps.h>
#include <OpenDIAStage2.h>
#include <OpenDIATraces.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <unordered_set>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr std::size_t K12 = FragmentGrid::FRAGMENTS;

    std::string annotation(const OpenSwath::LightTransition& t) { return t.getAnnotation(); }

    /// " rss_mb=<current> peak_mb=<peak>" for stage log lines (diagnostic only).
    std::string memoryNote()
    {
      std::size_t now = 0, peak = 0;
      SysInfo::getProcessMemoryConsumption(now);
      SysInfo::getProcessPeakMemoryConsumption(peak);
      return " rss_mb=" + std::to_string(now / 1024) + " peak_mb=" + std::to_string(peak / 1024);
    }

    std::string unmodified(const std::string& sequence)
    {
      try
      {
        return AASequence::fromString(sequence).toUnmodifiedString();
      }
      catch (const Exception::BaseException&)
      {
        std::string out;
        int depth = 0;
        for (char c : sequence)
        {
          if (c == '(' || c == '[') ++depth;
          else if (c == ')' || c == ']') --depth;
          else if (depth == 0 && c >= 'A' && c <= 'Z') out += c;
        }
        return out;
      }
    }

    double seconds(std::chrono::steady_clock::time_point since)
    {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - since).count();
    }

    /// Query list for StreamSearch from a (calibrated) library, in the heavy reader's transition order.
    std::vector<TransitionQuery> makeQueries(OpenSwath::LightTargetedExperiment& targeted)
    {
      struct Coordinates
      {
        double seconds = 0.0, im = 0.0;
        bool has_rt = false, has_im = false;
      };
      std::unordered_map<std::string, Coordinates> coordinates;
      coordinates.reserve(targeted.getCompounds().size());
      for (const auto& c : targeted.getCompounds())
      {
        coordinates.emplace(c.id, Coordinates{c.rt, c.drift_time, std::isfinite(c.rt), std::isfinite(c.drift_time) && c.drift_time > 0.0});
      }
      struct Key
      {
        double precursor_mz, product_mz, seconds;
        bool has_rt;
        std::size_t index;
      };
      auto& transitions = targeted.transitions;
      std::vector<Key> keys(transitions.size());
      for (std::size_t i = 0; i < transitions.size(); ++i)
      {
        const auto it = coordinates.find(transitions[i].peptide_ref);
        const bool has_rt = it != coordinates.end() && it->second.has_rt;
        keys[i] = {transitions[i].precursor_mz, transitions[i].product_mz, has_rt ? it->second.seconds : 0.0, has_rt, i};
      }
      std::stable_sort(keys.begin(), keys.end(), [&transitions](const Key& l, const Key& r)
      {
        if (l.precursor_mz != r.precursor_mz) return l.precursor_mz < r.precursor_mz;
        if (l.product_mz != r.product_mz) return l.product_mz < r.product_mz;
        if (l.has_rt != r.has_rt) return !l.has_rt;
        if (l.has_rt && l.seconds != r.seconds) return l.seconds < r.seconds;
        const auto& a = transitions[l.index];
        const auto& b = transitions[r.index];
        if (a.transition_name != b.transition_name) return a.transition_name < b.transition_name;
        return a.peptide_ref < b.peptide_ref;
      });
      std::vector<OpenSwath::LightTransition> sorted;
      sorted.reserve(transitions.size());
      for (const Key& k : keys) sorted.push_back(std::move(transitions[k.index]));
      transitions = std::move(sorted);
      std::vector<TransitionQuery> queries;
      queries.reserve(transitions.size());
      for (const auto& t : transitions)
      {
        TransitionQuery q;
        q.precursor_id = t.peptide_ref;
        q.precursor_mz = t.precursor_mz;
        q.fragment_mz = t.product_mz;
        const auto it = coordinates.find(t.peptide_ref);
        if (it != coordinates.end())
        {
          q.library_rt = it->second.seconds;
          q.has_library_rt = it->second.has_rt;
          q.precursor_im = it->second.im;
          q.has_precursor_im = it->second.has_im;
        }
        q.is_decoy = t.getDecoy();
        queries.push_back(q);
      }
      return queries;
    }

    void readSpectra(const std::string& path, Interfaces::IMSDataConsumer& consumer, std::vector<Int> levels)
    {
      FileHandler reader;
      PeakFileOptions options;
      options.setMSLevels(levels);
      reader.setOptions(options);
      reader.transformExperiment(path, &consumer, {}, ProgressLogger::NONE);
    }
  }

  const std::vector<std::string>& runFeatureNames()
  {
    static const std::vector<std::string> names = []
    {
      std::vector<std::string> n;
      for (const auto& s : stage2SearchFeatureNames()) n.push_back(s);
      for (const auto& s : stage2BoundaryFeatureNames()) n.push_back(s);
      for (const auto& s : stage2ByIonFeatureNames()) n.push_back(s);
      for (const char* s : FragmentGrid::evidenceNames()) n.emplace_back(s);
      for (const char* s : FragmentGrid::deficitNames()) n.emplace_back(s);
      for (const char* s : InterferenceResolver::featureNames()) n.emplace_back(s);
      return n;
    }();
    return names;
  }

  // ------------------------------------------------------------------------------------------------ library

  AnalysisLibrary AnalysisLibrary::read(const std::string& search_pqp, const std::string& scoring_pqp)
  {
    AnalysisLibrary lib;
    TransitionPQPFile().convertPQPToTargetedExperiment(search_pqp.c_str(), lib.search);
    const auto traml = TransitionPQPFile().getPQPCurrentIDToTraMLIDMap(search_pqp.c_str(), "PRECURSOR");
    lib.traml_ids.assign(traml.begin(), traml.end());
    std::unordered_map<std::string, std::string> traml_to_id;
    for (const auto& [id, name] : traml) traml_to_id.emplace(name, id);
    // Precursors from the compounds.
    for (const auto& c : lib.search.getCompounds())
    {
      Precursor p;
      p.id = c.id;
      const auto t = traml.find(c.id);
      if (t == traml.end()) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Precursor without TRAML_ID", c.id);
      p.traml_id = t->second;
      p.sequence = c.sequence;
      p.charge = c.charge;
      p.proteins = c.protein_refs;
      p.library_rt = c.rt;
      p.has_library_rt = std::isfinite(c.rt);
      lib.index.emplace(p.id, static_cast<std::uint32_t>(lib.precursors.size()));
      lib.precursors.push_back(std::move(p));
    }
    // Search fragments: m/z, decoy flag, annotations.
    std::vector<std::unordered_set<std::string>> search_annotations(lib.precursors.size());
    std::vector<std::vector<const OpenSwath::LightTransition*>> fragments(lib.precursors.size());
    for (const auto& t : lib.search.getTransitions())
    {
      const auto it = lib.index.find(t.peptide_ref);
      if (it == lib.index.end()) continue;
      auto& p = lib.precursors[it->second];
      p.mz = t.precursor_mz;
      p.decoy = t.getDecoy();
      search_annotations[it->second].insert(annotation(t));
      fragments[it->second].push_back(&t);
    }
    // Pairs and families: a decoy belongs to the target whose TRAML_ID it carries after DECOY_.
    std::unordered_map<std::string, std::uint32_t> pair_index;
    for (auto& p : lib.precursors)
    {
      p.pair_key = p.traml_id;
      if (p.decoy)
      {
        if (p.traml_id.rfind("DECOY_", 0) != 0)
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decoy TRAML_ID must be DECOY_<target TRAML_ID>", p.traml_id);
        p.pair_key = p.traml_id.substr(6);
      }
    }
    for (auto& p : lib.precursors)
    {
      std::string family_sequence = p.sequence;
      if (p.decoy)
      {
        const auto target = traml_to_id.find(p.pair_key);
        if (target != traml_to_id.end()) family_sequence = lib.precursors[lib.index.at(target->second)].sequence;
      }
      p.family = unmodified(family_sequence);
      p.fold = static_cast<std::uint8_t>(CrossFitScorer::fnv1a(p.family) % 3);
      p.pair = pair_index.try_emplace(p.pair_key, static_cast<std::uint32_t>(pair_index.size())).first->second;
    }
    lib.pairs = pair_index.size();
    // A pair shares the target's fold.
    {
      std::vector<int> pair_fold(lib.pairs, -1);
      for (const auto& p : lib.precursors)
        if (!p.decoy) pair_fold[p.pair] = p.fold;
      for (auto& p : lib.precursors)
        if (pair_fold[p.pair] >= 0) p.fold = static_cast<std::uint8_t>(pair_fold[p.pair]);
    }
    // Scoring fragments.
    const auto assign = [&](std::uint32_t pi, std::vector<std::pair<long, const OpenSwath::LightTransition*>>& list)
    {
      std::sort(list.begin(), list.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      auto& p = lib.precursors[pi];
      p.scoring_count = static_cast<std::uint8_t>(std::min(K12, list.size()));
      double maximum = 0.0;
      for (std::size_t s = 0; s < p.scoring_count; ++s) maximum = std::max(maximum, list[s].second->library_intensity);
      for (std::size_t s = 0; s < p.scoring_count; ++s)
      {
        p.scoring_mz[s] = list[s].second->product_mz;
        p.scoring_library[s] = list[s].second->library_intensity / std::max(maximum, 1e-12);
        p.scoring_in_search[s] = search_annotations[pi].contains(annotation(*list[s].second));
      }
    };
    if (scoring_pqp.empty())
    {
      for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
      {
        std::vector<std::pair<long, const OpenSwath::LightTransition*>> list;
        for (const auto* t : fragments[pi]) list.emplace_back(std::stol(t->transition_name), t);
        assign(pi, list);
      }
    }
    else
    {
      OpenSwath::LightTargetedExperiment scoring;
      TransitionPQPFile().convertPQPToTargetedExperiment(scoring_pqp.c_str(), scoring);
      std::vector<std::vector<std::pair<long, const OpenSwath::LightTransition*>>> lists(lib.precursors.size());
      for (const auto& t : scoring.getTransitions())
      {
        const auto it = lib.index.find(t.peptide_ref);
        if (it != lib.index.end()) lists[it->second].emplace_back(std::stol(t.transition_name), &t);
      }
      for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi) assign(pi, lists[pi]);
    }
    return lib;
  }

  // ------------------------------------------------------------------------------------------------ run

  RunAnalysis::RunAnalysis(AnalysisLibrary& library, RunOptions options, Log log) :
    library_(library), options_(std::move(options)), log_(std::move(log))
  {
  }

  RunResult RunAnalysis::run(const std::string& path)
  {
    RunResult result;
    result.run = path;
    const auto say = [&](const std::string& line)
    {
      result.log.push_back(line + memoryNote());
      if (log_) log_(result.log.back());
    };
    auto& lib = library_;
    const auto t_run = std::chrono::steady_clock::now();

    // 1. Calibration (uncalibrated library), or a given calibration and its anchors.
    OpenSwath::LightTargetedExperiment targeted = lib.search;
    std::vector<std::pair<std::string, double>> anchors;  // used anchors: precursor id, run RT
    if (!options_.calibration_json.empty())
    {
      result.calibration = readCalibrationRT(options_.calibration_json);
      std::ifstream in(options_.anchors_tsv);
      std::string line;
      if (!in || !std::getline(in, line)) throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, options_.anchors_tsv);
      std::vector<std::string> header;
      StringUtils::split(line, '\t', header);
      const auto col = [&](const std::string& name)
      {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end()) throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "anchors: column " + name);
        return static_cast<std::size_t>(it - header.begin());
      };
      const std::size_t c_id = col("precursor_id"), c_rt = col("run_rt"), c_used = col("used");
      while (std::getline(in, line))
      {
        if (line.empty()) continue;
        std::vector<std::string> v;
        StringUtils::split(line, '\t', v);
        if (v[c_used] == "1") anchors.emplace_back(v[c_id], StringUtils::toDouble(v[c_rt]));
      }
      say("Calibration: reused " + options_.calibration_json + ", used anchors=" + std::to_string(anchors.size()));
    }
    else
    {
      CalibrationOptions copt = options_.calibration;
      std::optional<CalibrationResult> calibrated;
      std::optional<CalibrationWindow> window;
      std::size_t passes = 0;
      while (!calibrated)
      {
        CalibrationCollector collector(targeted, lib.traml_ids, copt, window);
        collector.setIonMobilityArrayAsVSSC(options_.im_array_as_vssc);
        readSpectra(path, collector, {2});
        ++passes;
        say("Calibration: pass " + std::to_string(passes) + " collected panel=" + std::to_string(collector.panelSize()));
        CalibrationDiagnostics diagnostics;
        try
        {
          calibrated = collector.calibrate(options_.window_panel_targets ? &diagnostics : nullptr);
        }
        catch (const std::runtime_error& error)
        {
          if (!window && options_.window_panel_targets > copt.panel_targets && diagnostics.window)
          {
            say(std::string("Calibration: whole-run pass failed (") + error.what() + "); windowed pass with up to " +
                std::to_string(options_.window_panel_targets) + " targets.");
            window = diagnostics.window;
            copt.panel_targets = options_.window_panel_targets;
            continue;
          }
          throw;
        }
      }
      std::ostringstream json;
      writeCalibrationJSON(*calibrated, copt, json);
      result.calibration_json = json.str();
      std::istringstream in(result.calibration_json);
      result.calibration = readCalibrationRT(in, path);
      for (const auto& a : calibrated->anchors)
        if (a.used) anchors.emplace_back(a.precursor_id, a.run_rt);
      say("Calibration: passes=" + std::to_string(passes) + " used anchors=" + std::to_string(anchors.size()));
    }
    const CalibrationRT& cal = result.calibration;
    if (!cal.has_ms2) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Calibration has no fragment mass model", path);
    say("Calibration: rt_tolerance_s=" + std::to_string(cal.tolerance_s) + " ms2_offset_ppm=" + std::to_string(cal.ms2_offset_ppm) +
        " ms2_tolerance_ppm=" + std::to_string(cal.ms2_tolerance_ppm) + " boundary_half_widths=" + std::to_string(cal.boundary_half_widths) +
        " seconds=" + std::to_string(seconds(t_run)));

    // 2. Search with the calibrated library.
    const auto t_search = std::chrono::steady_clock::now();
    for (auto& c : targeted.compounds)
    {
      if (std::isfinite(c.rt)) c.rt = cal.map(c.rt);
      if (cal.has_im && std::isfinite(c.drift_time) && c.drift_time > 0.0) c.drift_time = cal.im_intercept + cal.im_slope * c.drift_time;
    }
    for (auto& t : targeted.transitions) t.product_mz *= 1.0 + cal.ms2_offset_ppm * 1e-6;
    std::vector<std::uint32_t> search_to_library;
    std::vector<SearchCandidate> candidates;
    {
      SearchOptions so;
      so.rt_window_seconds = cal.tolerance_s;
      so.fragment_tolerance_ppm = cal.ms2_tolerance_ppm;
      if (cal.has_im) so.ion_mobility_tolerance_vssc = cal.im_tolerance;
      so.ion_mobility_array_as_vssc = options_.im_array_as_vssc;
      so.picker.max_peaks = options_.max_peaks;
      so.picker.boundary_half_widths = cal.boundary_half_widths;
      so.min_transitions = options_.min_transitions;
      // The peak table text is not needed: candidates are kept in memory. A stream without a buffer discards it.
      std::ostream sink(nullptr);
      StreamSearch search(makeQueries(targeted), so, sink, nullptr);
      search.keepCandidates();
      readSpectra(path, search, {2});
      search.finish();
      candidates = search.candidates();
      search_to_library.resize(search.precursorCount());
      for (std::size_t p = 0; p < search.precursorCount(); ++p) search_to_library[p] = lib.index.at(search.precursorId(p));
    }
    say("Search: candidates=" + std::to_string(candidates.size()) + " seconds=" + std::to_string(seconds(t_search)));

    // 3. Population: every candidate of pairs whose best candidate has enough supporting fragments.
    std::vector<std::uint32_t> pair_support(lib.pairs, 0);
    for (const auto& c : candidates)
    {
      auto& s = pair_support[lib.precursors[search_to_library[c.precursor]].pair];
      s = std::max(s, c.transition_support);
    }
    for (std::uint32_t row = 0; row < candidates.size(); ++row)
    {
      const auto& c = candidates[row];
      const std::uint32_t pi = search_to_library[c.precursor];
      if (pair_support[lib.precursors[pi].pair] < options_.population_min_support) continue;
      RunWindow w;
      w.precursor = pi;
      w.candidate = row;
      w.peak = c;
      result.windows.push_back(w);
    }
    std::vector<SearchCandidate>().swap(candidates);
    auto& windows = result.windows;
    const std::size_t n = windows.size();

    // 4. One pass: 6-fragment traces + b/y ions (TraceCollector) and the 12-fragment grid (FragmentGridCollector).
    const auto t_extract = std::chrono::steady_clock::now();
    std::vector<std::uint32_t> local_of(lib.precursors.size(), UINT32_MAX);
    std::vector<std::uint32_t> local_precursors;
    for (const auto& w : windows)
    {
      if (local_of[w.precursor] == UINT32_MAX)
      {
        local_of[w.precursor] = static_cast<std::uint32_t>(local_precursors.size());
        local_precursors.push_back(w.precursor);
      }
    }
    std::vector<FragmentGridCollector::Precursor> grid_precursors(local_precursors.size());
    std::vector<InterferenceResolver::Precursor> resolver_precursors(local_precursors.size());
    for (std::size_t l = 0; l < local_precursors.size(); ++l)
    {
      const auto& p = lib.precursors[local_precursors[l]];
      auto& g = grid_precursors[l];
      g.mz = p.mz;
      for (std::size_t s = 0; s < K12; ++s)
      {
        g.extraction_mz[s] = p.scoring_mz[s] > 0.0 ? p.scoring_mz[s] * (1.0 + cal.ms2_offset_ppm * 1e-6) : 0.0;
        g.library[s] = p.scoring_library[s];
        g.search[s] = p.scoring_in_search[s];
        resolver_precursors[l].fragment_mz[s] = p.scoring_mz[s];
        resolver_precursors[l].library[s] = p.scoring_library[s];
      }
    }
    std::vector<InterferenceResolver::Window> grid_windows(n);
    std::vector<TraceWindow> trace_windows(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      const auto& w = windows[i];
      grid_windows[i].precursor = local_of[w.precursor];
      grid_windows[i].rt_start = w.peak.rt_start;
      grid_windows[i].rt_apex = w.peak.rt_apex;
      grid_windows[i].rt_end = w.peak.rt_end;
      trace_windows[i] = {std::to_string(w.candidate), lib.precursors[w.precursor].id, lib.precursors[w.precursor].family, -1, w.peak.rt_apex};
    }
    FragmentGridCollector grids(grid_precursors, grid_windows, options_.grid_half_width_s, options_.grid_extraction_ppm);
    // The 6-fragment tensor uses the calibrated search library (as OpenDIACNN -apply_calibration_coordinates).
    OpenSwath::LightTargetedExperiment trace_library = lib.search;
    for (auto& t : trace_library.transitions) t.product_mz *= 1.0 + cal.ms2_offset_ppm * 1e-6;
    TraceCollector traces(std::move(trace_windows), trace_library, 0.02, options_.grid_half_width_s, 0.03, n, std::size_t{64} << 30);
    traces.setIonMobilityArrayAsVSSC(options_.im_array_as_vssc);
    traces.enableEvidence(20.0, options_.core_half_width_s);
    traces.enableByIonCounts(trace_library, options_.byion_ppm, options_.byion_intensity_min);
    {
      TeeConsumer both(grids, traces);
      readSpectra(path, both, {1, 2});
    }
    grids.finish();
    const auto isolation = grids.isolationIndices();
    for (std::size_t l = 0; l < resolver_precursors.size(); ++l) resolver_precursors[l].isolation = isolation[l];
    say("Extraction: windows=" + std::to_string(n) + " precursors=" + std::to_string(local_precursors.size()) +
        " isolation_windows=" + std::to_string(grids.isolationWindowCount()) + " peak_open_grids=" + std::to_string(grids.peakOpenGrids()) +
        " seconds=" + std::to_string(seconds(t_extract)));

    // 5. Features.
    const auto t_features = std::chrono::steady_clock::now();
    const auto& names = runFeatureNames();
    const std::size_t d = names.size();
    result.feature_names = names;
    result.features.assign(n * d, 0.0f);
    const double noise = FragmentGrid::noiseLevel(grid_windows, resolver_precursors);
    const auto interference = InterferenceResolver::compute(resolver_precursors, grid_windows, cal.ms2_tolerance_ppm,
                                                            2.0 * options_.grid_half_width_s / InterferenceResolver::BINS, options_.threads);
    Stage2ExtractionOptions extraction;
    extraction.rt_half_width = options_.grid_half_width_s;
    extraction.core_half_width = options_.core_half_width_s;
#pragma omp parallel num_threads(options_.threads)
    {
      std::array<float, TraceCollector::SAMPLE_SIZE> view;
#pragma omp for schedule(dynamic, 1024)
      for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii)
      {
        const std::size_t i = static_cast<std::size_t>(ii);
        const auto& w = windows[i];
        float* f = result.features.data() + i * d;
        std::size_t o = 0;
        SearchPeakFeatures peak;
        peak.transition_support = w.peak.transition_support;
        peak.matched_scans = w.peak.matched_scans;
        peak.peak_width_seconds = w.peak.rt_end - w.peak.rt_start;
        peak.abs_rt_error_seconds = lib.precursors[w.precursor].has_library_rt ? std::abs(w.peak.apex_rt_error_seconds) : 0.0;
        for (double v : stage2SearchFeatures(peak)) f[o++] = static_cast<float>(v);
        const float* tensor = traces.windowTensor(i, view);
        for (double v : stage2BoundaryFeatures(tensor, traces.libraryIntensities(i), traces.fragmentMZ(i), traces.evidence(i), extraction)) f[o++] = static_cast<float>(v);
        for (double v : stage2ByIonFeatures(traces.evidence(i))) f[o++] = static_cast<float>(v);
        for (float v : grids.evidence()[i]) f[o++] = v;
        const auto& gp = grid_precursors[grid_windows[i].precursor];
        for (float v : FragmentGrid::deficit(grid_windows[i].core_area, gp.library, gp.search, noise)) f[o++] = v;
        for (std::size_t k = 0; k < InterferenceResolver::FEATURES; ++k) f[o++] = interference[i * InterferenceResolver::FEATURES + k];
      }
    }
    say("Features: " + std::to_string(d) + " per window, noise=" + std::to_string(noise) + " seconds=" + std::to_string(seconds(t_features)));

    // 6. Scorer.
    const auto t_score = std::chrono::steady_clock::now();
    CrossFitScorer::Input input;
    input.dimension = d;
    input.features = result.features;
    std::vector<std::uint32_t> dense_precursor(lib.precursors.size(), UINT32_MAX), dense_pair(lib.pairs, UINT32_MAX);
    std::uint32_t np = 0, npair = 0;
    for (const auto& w : windows)
    {
      const auto& p = lib.precursors[w.precursor];
      if (dense_precursor[w.precursor] == UINT32_MAX) dense_precursor[w.precursor] = np++;
      if (dense_pair[p.pair] == UINT32_MAX) dense_pair[p.pair] = npair++;
      input.precursor.push_back(dense_precursor[w.precursor]);
      input.pair.push_back(dense_pair[p.pair]);
      input.decoy.push_back(p.decoy);
      input.fold.push_back(p.fold);
      input.pair_hash.push_back(CrossFitScorer::fnv1a("v" + p.pair_key));
    }
    // Seeds: the window nearest each used anchor (|apex - anchor run RT| <= anchor_apex_s) is positive; the decoy partner's
    // best search window (most supporting fragments, then most scans, then lowest row) and a sample of other decoy windows
    // are negative.
    std::unordered_map<std::uint32_t, double> anchor_rt;
    for (const auto& [id, rt] : anchors)
    {
      const auto it = lib.index.find(id);
      if (it != lib.index.end()) anchor_rt.emplace(it->second, rt);
    }
    std::vector<std::uint32_t> best_of(lib.precursors.size(), UINT32_MAX), anchor_window(lib.precursors.size(), UINT32_MAX);
    for (std::uint32_t i = 0; i < n; ++i)
    {
      const auto& w = windows[i];
      auto& b = best_of[w.precursor];
      if (b == UINT32_MAX || w.peak.transition_support > windows[b].peak.transition_support ||
          (w.peak.transition_support == windows[b].peak.transition_support && w.peak.matched_scans > windows[b].peak.matched_scans))
      {
        b = i;
      }
      const auto a = anchor_rt.find(w.precursor);
      if (a == anchor_rt.end()) continue;
      const double dt = std::abs(w.peak.rt_apex - a->second);
      auto& aw = anchor_window[w.precursor];
      if (dt <= options_.anchor_apex_s && (aw == UINT32_MAX || dt < std::abs(windows[aw].peak.rt_apex - a->second))) aw = i;
    }
    std::unordered_map<std::string, std::uint32_t> by_traml;
    for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi) by_traml.emplace(lib.precursors[pi].traml_id, pi);
    std::unordered_set<std::uint32_t> negative_rows;
    std::size_t positives = 0;
    for (std::uint32_t pi = 0; pi < lib.precursors.size(); ++pi)
    {
      if (anchor_window[pi] == UINT32_MAX) continue;
      const auto& p = lib.precursors[pi];
      input.seeds.push_back({anchor_window[pi], true, CrossFitScorer::fnv1a("v" + p.family)});
      ++positives;
      const auto partner = by_traml.find("DECOY_" + p.traml_id);
      if (partner != by_traml.end() && best_of[partner->second] != UINT32_MAX && negative_rows.insert(best_of[partner->second]).second)
      {
        input.seeds.push_back({best_of[partner->second], false, CrossFitScorer::fnv1a("v" + lib.precursors[partner->second].family)});
      }
    }
    {
      std::vector<std::uint32_t> decoys;
      for (std::uint32_t i = 0; i < n; ++i)
        if (lib.precursors[windows[i].precursor].decoy && !negative_rows.contains(i)) decoys.push_back(i);
      // Deterministic sample (splitmix64 partial shuffle).
      std::uint64_t state = options_.seed + 13;
      const auto next = [&state]()
      {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
      };
      const std::size_t take = std::min(options_.decoy_seeds, decoys.size());
      for (std::size_t k = 0; k < take; ++k)
      {
        std::swap(decoys[k], decoys[k + next() % (decoys.size() - k)]);
        const auto& p = lib.precursors[windows[decoys[k]].precursor];
        input.seeds.push_back({decoys[k], false, CrossFitScorer::fnv1a("v" + p.family)});
      }
    }
    say("Seeds: positives=" + std::to_string(positives) + " negatives=" + std::to_string(input.seeds.size() - positives));
    CrossFitScorer::Options so;
    so.seed = options_.seed;
    so.threads = options_.threads;
    const auto scored = CrossFitScorer::run(input, so, [&](const std::string& l) { say("Scorer " + l); });
    for (std::size_t i = 0; i < n; ++i)
    {
      windows[i].score = scored.score[i];
      windows[i].q_value = scored.q_value[i];
    }
    std::size_t at1 = 0;
    for (std::size_t i = 0; i < n; ++i) at1 += !lib.precursors[windows[i].precursor].decoy && windows[i].q_value <= 0.01;
    say("Scoring: target precursors at q <= 0.01: " + std::to_string(at1) + " seconds=" + std::to_string(seconds(t_score)) +
        "; run total seconds=" + std::to_string(seconds(t_run)));
    return result;
  }
}
