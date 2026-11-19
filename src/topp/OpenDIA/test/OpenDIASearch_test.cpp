// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      std::fprintf(stderr, "FAILED: %s\n", message);
      throw std::runtime_error(message);
    }
  }

  MSSpectrum spectrum(double rt, double isolation_mz, std::vector<Peak1D> peaks, double half_width = 50.0)
  {
    MSSpectrum result;
    result.setMSLevel(2);
    result.setRT(rt);
    Precursor precursor;
    precursor.setMZ(isolation_mz);
    precursor.setIsolationWindowLowerOffset(half_width);
    precursor.setIsolationWindowUpperOffset(half_width);
    result.setPrecursors({precursor});
    for (const auto& peak : peaks)
    {
      result.push_back(peak);
    }
    return result;
  }

  MSSpectrum spectrumWithIM(double rt, double isolation_mz, std::vector<Peak1D> peaks, std::initializer_list<float> mobilities)
  {
    MSSpectrum result = spectrum(rt, isolation_mz, std::move(peaks));
    require(result.size() == mobilities.size(), "test spectrum peak and IM array lengths differ");
    DataArrays::FloatDataArray im_array;
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    im_array.assign(mobilities.begin(), mobilities.end());
    result.getFloatDataArrays().push_back(std::move(im_array));
    return result;
  }

  struct Row
  {
    std::string id;
    double mz = 0.0, start = 0.0, apex = 0.0, end = 0.0, score = 0.0;
    std::size_t support = 0, scans = 0;
  };

  std::vector<Row> parse(const std::string& text)
  {
    std::vector<Row> rows;
    std::istringstream lines(text);
    for (std::string line; std::getline(lines, line);)
    {
      std::istringstream fields(line);
      Row row;
      fields >> row.id >> row.mz >> row.start >> row.apex >> row.end >> row.score >> row.support >> row.scans;
      rows.push_back(row);
    }
    return rows;
  }

  std::size_t count(const std::vector<Row>& rows, const std::string& id)
  {
    return static_cast<std::size_t>(std::count_if(rows.begin(), rows.end(), [&id](const Row& row) { return row.id == id; }));
  }

  OpenDIA::SearchOptions defaults()
  {
    OpenDIA::SearchOptions options;
    options.fragment_tolerance_da = 0.01;
    options.rt_window_seconds = 1000.0;
    return options;
  }

  /// Runs a search over @p spectra and returns the peak table (and the feature table in @p features, if given).
  std::string run(const std::vector<OpenDIA::TransitionQuery>& queries, const OpenDIA::SearchOptions& options,
                  std::vector<MSSpectrum> spectra, std::string* features = nullptr, bool* im_used = nullptr)
  {
    std::ostringstream output, feature_output;
    OpenDIA::StreamSearch search(queries, options, output, features != nullptr ? &feature_output : nullptr);
    for (auto& s : spectra)
    {
      search.consumeSpectrum(s);
    }
    search.finish();
    if (features != nullptr)
    {
      *features = feature_output.str();
    }
    if (im_used != nullptr)
    {
      *im_used = search.ionMobilityFilteringUsed();
    }
    return output.str();
  }

  template <typename Action> bool throws(Action action)
  {
    try
    {
      action();
    }
    catch (const Exception::BaseException&)
    {
      return true;
    }
    return false;
  }

  double bell(double rt, double centre, double height = 1000.0)
  {
    return height * std::exp(-0.5 * (rt - centre) * (rt - centre) / 4.0);
  }
}

int main()
{
  // One eluting precursor inside the isolation window and one outside it.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {
      {"analyte", 500.0, 100.0, 20.0, true}, {"analyte", 500.0, 200.0, 20.0, true}, {"other", 700.0, 300.0, 20.0, true}};
    OpenDIA::SearchOptions options = defaults();
    options.min_transitions = 2;
    std::vector<MSSpectrum> spectra;
    for (int i = 0; i < 20; ++i)
    {
      const double rt = 10.0 + i;
      spectra.push_back(spectrum(rt, 500.0, {Peak1D(100.0, bell(rt, 20.0)), Peak1D(200.0, 0.5 * bell(rt, 20.0)), Peak1D(300.0, 999.0)}));
    }
    const auto rows = parse(run(queries, options, spectra));
    require(rows.size() == 1 && rows[0].id == "analyte", "expected one peak of the isolated precursor only");
    require(rows[0].apex == 20.0 && rows[0].support == 2 && rows[0].start < 20.0 && rows[0].end > 20.0, "apex, support or boundaries off");
  }

  // Transitions of one precursor must agree on the decoy label.
  require(throws([] {
            std::ostringstream sink;
            OpenDIA::StreamSearch search({{"mixed", 500.0, 100.0, 0.0, true, 0.0, false, false}, {"mixed", 500.0, 101.0, 0.0, true, 0.0, false, true}},
                                         defaults(), sink);
          }),
          "transitions with inconsistent target/decoy labels were accepted");

  // Invalid settings.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {{"analyte", 500.0, 100.0, 0.0, true}};
    for (const auto& change : std::vector<std::function<void(OpenDIA::SearchOptions&)>>{
           [](auto& o) { o.fragment_tolerance_da = 0.0; }, [](auto& o) { o.ion_mobility_tolerance_vssc = 0.0; },
           [](auto& o) { o.rt_window_seconds = 0.0; }, [](auto& o) { o.min_transitions = 0; }, [](auto& o) { o.picker.max_peaks = 0; }})
    {
      OpenDIA::SearchOptions options = defaults();
      change(options);
      require(throws([&] { std::ostringstream sink; OpenDIA::StreamSearch search(queries, options, sink); }), "invalid search option accepted");
    }
  }

  // Isolation metadata: no width is ignored, several and overlapping windows count once, edges are inclusive.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {{"analyte", 500.0, 100.0, 0.0, false}};
    MSSpectrum no_width = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
    no_width.getPrecursors().front().setIsolationWindowLowerOffset(0.0);
    no_width.getPrecursors().front().setIsolationWindowUpperOffset(0.0);
    require(run(queries, defaults(), {no_width}).empty(), "spectrum with no isolation-width metadata was matched");

    MSSpectrum multi = spectrum(10.0, 500.0, {Peak1D(300.0, 9.0)});
    for (const double centre : {700.0, 710.0})
    {
      Precursor window;
      window.setMZ(centre);
      window.setIsolationWindowLowerOffset(50.0);
      window.setIsolationWindowUpperOffset(50.0);
      multi.getPrecursors().push_back(window);
    }
    const auto multi_rows = parse(run({{"other", 700.0, 300.0, 0.0, false}}, defaults(), {multi}));
    require(multi_rows.size() == 1 && multi_rows[0].apex == 10.0 && multi_rows[0].scans == 1, "multiple or overlapping isolation windows were not one scan");

    const std::vector<OpenDIA::TransitionQuery> edges = {
      {"lower_edge", 450.0, 100.0, 0.0, false}, {"upper_edge", 550.0, 101.0, 0.0, false},
      {"below_edge", 449.999, 102.0, 0.0, false}, {"above_edge", 550.001, 103.0, 0.0, false}};
    const auto edge_rows = parse(run(edges, defaults(), {spectrum(10.0, 500.0, {Peak1D(100.0, 10.0), Peak1D(101.0, 10.0), Peak1D(102.0, 10.0), Peak1D(103.0, 10.0)})}));
    require(count(edge_rows, "lower_edge") == 1 && count(edge_rows, "upper_edge") == 1, "inclusive isolation-window edge candidates were lost");
    require(count(edge_rows, "below_edge") == 0 && count(edge_rows, "above_edge") == 0, "candidate outside the isolation window was reported");
  }

  // Interleaved lanes: each precursor's trace holds only the scans that isolated it.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {{"lane_a", 500.0, 100.0, 0.0, false}, {"lane_b", 700.0, 100.0, 0.0, false}};
    const auto rows = parse(run(queries, defaults(), {spectrum(10.0, 500.0, {Peak1D(100.0, 2.0)}), spectrum(20.0, 700.0, {Peak1D(100.0, 5.0)}),
                                                      spectrum(30.0, 500.0, {Peak1D(100.0, 4.0)})}));
    require(rows.size() == 2, "expected one peak per lane");
    for (const Row& row : rows)
    {
      if (row.id == "lane_a")
      {
        require(row.apex == 30.0 && row.scans == 2 && row.start == 10.0 && row.end == 30.0, "lane_a trace mixed lanes");
      }
      else
      {
        require(row.apex == 20.0 && row.scans == 1, "lane_b trace mixed lanes");
      }
    }
  }

  // RT window: a scan outside a precursor's window never enters its trace; precursors without a library RT are
  // searched over the whole run and written at the end.
  {
    OpenDIA::SearchOptions options = defaults();
    options.rt_window_seconds = 10.0;
    const std::vector<OpenDIA::TransitionQuery> queries = {{"early", 500.0, 100.0, 20.0, true}, {"untimed", 500.0, 100.0, 0.0, false}};
    std::vector<MSSpectrum> spectra;
    for (int i = 0; i < 100; ++i)
    {
      const double rt = 1.0 + i;
      spectra.push_back(spectrum(rt, 500.0, {Peak1D(100.0, bell(rt, 20.0) + bell(rt, 80.0, 3000.0))}));
    }
    const std::string text = run(queries, options, spectra);
    const auto rows = parse(text);
    require(count(rows, "early") == 1 && count(rows, "untimed") == 2, "RT window gating or untimed search changed");
    for (const Row& row : rows)
    {
      if (row.id == "early")
      {
        require(row.apex == 20.0 && row.start >= 10.0 && row.end <= 30.0, "scan outside the RT window entered a trace");
      }
    }
    require(text.rfind("early", 0) == 0, "rows must be written in order of RT-window end, untimed precursors last");
  }

  // Ion mobility: applied only when both the library and the spectrum provide 1/K0.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {{"target", 500.0, 100.0, 0.0, false, 0.80, true}, {"interferer", 500.0, 100.0, 0.0, false, 0.95, true}};
    OpenDIA::SearchOptions options = defaults();
    options.ion_mobility_tolerance_vssc = 0.02;
    bool used = false;
    auto rows = parse(run(queries, options, {spectrumWithIM(10.0, 500.0, {Peak1D(100.0, 10.0)}, {0.801f}), spectrumWithIM(20.0, 500.0, {Peak1D(100.0, 20.0)}, {0.799f})},
                          nullptr, &used));
    require(used && count(rows, "target") == 1 && count(rows, "interferer") == 0, "mobility did not separate target and interferer");

    // Sorting a cell by m/z keeps each mobility with its centroid.
    rows = parse(run(queries, options, {spectrumWithIM(10.0, 500.0, {Peak1D(100.005, 50.0), Peak1D(100.0, 10.0)}, {0.95f, 0.80f})}));
    require(rows.size() == 2, "both mobility-separated candidates should match one centroid each");
    for (const Row& row : rows)
    {
      require(row.score == (row.id == "target" ? 10.0 : 50.0), "m/z sorting detached a mobility from its centroid");
    }

    // A unitless "Ion Mobility" array parses as milliseconds: ignored by default, 1/K0 only with the opt-in.
    for (const bool opt_in : {false, true})
    {
      MSSpectrum unitless = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
      DataArrays::FloatDataArray im_array;
      im_array.setName("Ion Mobility");
      im_array.push_back(0.80f);
      unitless.getFloatDataArrays().push_back(std::move(im_array));
      require(unitless.getIMData().second == DriftTimeUnit::MILLISECOND, "unitless Ion Mobility array no longer parses as ms");
      OpenDIA::SearchOptions unitless_options = options;
      unitless_options.ion_mobility_array_as_vssc = opt_in;
      rows = parse(run(queries, unitless_options, {unitless}, nullptr, &used));
      require(used == opt_in && (count(rows, "interferer") == 0) == opt_in, "unitless IM opt-in did not control filtering");
    }

    rows = parse(run(queries, options, {spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)})}, nullptr, &used));
    require(!used && rows.size() == 2, "m/z-only fallback changed without raw mobility");
    rows = parse(run({{"unannotated", 500.0, 100.0, 0.0, false}}, options, {spectrumWithIM(10.0, 500.0, {Peak1D(100.0, 10.0)}, {0.95f})}, nullptr, &used));
    require(!used && rows.size() == 1, "m/z-only fallback changed without library mobility");

    MSSpectrum scalar = spectrum(10.0, 500.0, {Peak1D(100.0, 10.0)});
    scalar.setDriftTime(0.95);
    scalar.setDriftTimeUnit(DriftTimeUnit::VSSC);
    rows = parse(run(queries, options, {scalar}, nullptr, &used));
    require(used && count(rows, "interferer") == 1 && count(rows, "target") == 0, "scalar VSSC mobility did not select the matching candidate");
  }

  // RT x mobility picking (im_bins): a weak target at its library mobility, overlapping in RT with a strong precursor at
  // another mobility inside the window, is merged into the strong peak in the pooled trace; the mobility band finds it.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {{"target", 500.0, 100.0, 0.0, false, 0.80, true},
                                                           {"target", 500.0, 200.0, 0.0, false, 0.80, true},
                                                           {"target", 500.0, 300.0, 0.0, false, 0.80, true}};
    std::vector<MSSpectrum> spectra;
    // Scans every 1 s, 0..40 s. Interferer at 1/K0 0.84 (inside +-0.05), apex 18 s, height 1000; target at 0.80, apex
    // 23 s, height 60. With a pooled maximum the target sits on the interferer's tail and is no separate maximum.
    for (int s = 0; s <= 40; ++s)
    {
      const double rt = s;
      const double inter = 1000.0 * std::exp(-0.5 * std::pow((rt - 18.0) / 2.5, 2));
      const double targ = 60.0 * std::exp(-0.5 * std::pow((rt - 23.0) / 1.5, 2));
      std::vector<Peak1D> peaks;
      std::vector<float> ims;
      for (const double mz : {100.0, 200.0, 300.0})
      {
        if (inter > 1.0) { peaks.emplace_back(mz, inter); ims.push_back(0.84f); }
        if (targ > 1.0) { peaks.emplace_back(mz + 0.001, targ); ims.push_back(0.80f); }
      }
      MSSpectrum spec = spectrum(rt, 500.0, peaks);
      DataArrays::FloatDataArray im_array;
      IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
      im_array.assign(ims.begin(), ims.end());
      spec.getFloatDataArrays().push_back(std::move(im_array));
      spectra.push_back(std::move(spec));
    }
    OpenDIA::SearchOptions options = defaults();
    options.ion_mobility_tolerance_vssc = 0.05;
    const auto pooled = parse(run(queries, options, spectra));
    const auto has_apex = [](const std::vector<Row>& rows, double rt) { return std::any_of(rows.begin(), rows.end(), [rt](const Row& r) { return std::abs(r.apex - rt) <= 1.0; }); };
    require(has_apex(pooled, 18.0) && !has_apex(pooled, 23.0), "test premise: the pooled trace hides the target peak");
    options.im_bins = 8;
    const auto binned = parse(run(queries, options, spectra));
    require(has_apex(binned, 23.0), "RT x mobility picking did not recover the masked target peak");
    require(has_apex(binned, 18.0), "RT x mobility picking dropped the pooled peak");
    // The pooled peaks are unchanged: every pooled row is reproduced first.
    for (std::size_t i = 0; i < pooled.size(); ++i)
    {
      require(binned[i].apex == pooled[i].apex && binned[i].score == pooled[i].score, "im_bins changed a pooled candidate");
    }
    // Without raw mobility, im_bins has no effect.
    std::vector<MSSpectrum> plain;
    for (const auto& s : spectra) { MSSpectrum p = s; p.getFloatDataArrays().clear(); plain.push_back(std::move(p)); }
    OpenDIA::SearchOptions off = defaults();
    off.ion_mobility_tolerance_vssc = 0.05;
    require(run(queries, off, plain) == run(queries, options, plain), "im_bins changed a search without raw mobility");
    require(throws([&] { OpenDIA::SearchOptions bad = defaults(); bad.im_bins = 2; std::ostringstream o; OpenDIA::StreamSearch(queries, bad, o); }),
            "im_bins 2 must be rejected");
  }

  // A ppm fragment tolerance scales with m/z: 10 ppm keeps a 0.008 Da error at 1000 but rejects it at 100.
  {
    OpenDIA::SearchOptions options = defaults();
    options.fragment_tolerance_ppm = 10.0;
    const auto rows = parse(run({{"low_mz", 500.0, 100.0, 0.0, false}, {"high_mz", 500.0, 1000.0, 0.0, false}}, options,
                                {spectrum(10.0, 500.0, {Peak1D(100.008, 10.0), Peak1D(1000.008, 10.0)})}));
    require(count(rows, "high_mz") == 1 && count(rows, "low_mz") == 0, "ppm tolerance wrong");
  }

  // Two elutions of one precursor separated by a silent gap, plus a second precursor; feature table and trace dump.
  {
    const std::vector<OpenDIA::TransitionQuery> queries = {
      {"twice", 500.0, 100.0, 200.0, true}, {"twice", 500.0, 200.0, 200.0, true},
      {"late", 520.0, 300.0, 400.0, true, 0.0, false, true}, {"late", 520.0, 310.0, 400.0, true, 0.0, false, true}};
    OpenDIA::SearchOptions options = defaults();
    options.rt_window_seconds = 120.0;
    options.min_transitions = 2;
    const auto scans = [](bool with_twice)
    {
      std::vector<MSSpectrum> spectra;
      for (int i = 0; i < 300; ++i)
      {
        const double rt = 100.0 + 1.5 * i;
        const double a = with_twice ? (rt < 210.0 ? bell(rt, 170.0) : 0.4 * bell(rt, 250.0)) : 0.0;
        const double b = 2.0 * bell(rt, 400.0);
        std::vector<Peak1D> peaks;
        if (a > 1.0) { peaks.emplace_back(100.0, a); peaks.emplace_back(200.0, 0.5 * a); }
        if (b > 1.0) { peaks.emplace_back(with_twice ? 300.0 : 300.001, b); peaks.emplace_back(310.0, 0.7 * b); }
        spectra.push_back(spectrum(rt, 510.0, peaks, 20.0));
      }
      return spectra;
    };
    std::string features;
    const auto rows = parse(run(queries, options, scans(true), &features));
    require(rows.size() == 3, "expected two peaks of 'twice' and one of 'late'");
    std::vector<double> twice_apexes;
    for (const Row& row : rows)
    {
      require(row.end - row.start < 30.0, "boundaries should follow the elution profile, not the RT window");
      if (row.id == "twice") { twice_apexes.push_back(row.apex); }
      else { require(std::abs(row.apex - 400.0) < 1.6, "apex of 'late' is off"); }
    }
    std::sort(twice_apexes.begin(), twice_apexes.end());
    require(twice_apexes.size() == 2 && std::abs(twice_apexes[0] - 170.0) < 1.6 && std::abs(twice_apexes[1] - 250.0) < 1.6,
            "the two elutions of 'twice' were not separated");

    std::istringstream feature_lines(features);
    std::vector<std::vector<std::string>> table;
    for (std::string line; std::getline(feature_lines, line);)
    {
      std::istringstream fields(line);
      table.emplace_back();
      for (std::string field; std::getline(fields, field, '\t');) { table.back().push_back(field); }
    }
    std::size_t header_columns = 0;
    {
      std::istringstream header(OpenDIA::StreamSearch::featureHeader());
      for (std::string field; std::getline(header, field, '\t');) { ++header_columns; }
    }
    require(table.size() == 3 && header_columns == 19, "one feature row per peak, 19 columns");
    for (std::size_t r = 0; r < table.size(); ++r)
    {
      const auto& t = table[r];
      require(t.size() == header_columns && t[0] == std::to_string(r), "feature row width or running row number wrong");
      require(t[3] == (t[1] == "late" ? "1" : "0") && t[12] == "1", "decoy label or library-RT flag wrong");
      require(std::stod(t[15]) > 3.0 && std::stod(t[15]) < 8.0, "FWHM of a sigma-2 s peak should be about 4.7 s");
      require(std::abs(std::stod(t[13]) - (std::stod(t[5]) - (t[1] == "late" ? 400.0 : 200.0))) < 1e-9, "apex RT error wrong");
    }

    // The trace dump writes the picker's matrix of the selected precursor only and leaves the peak output unchanged.
    std::ostringstream dump_output, dump_traces;
    OpenDIA::StreamSearch dump_search(queries, options, dump_output);
    dump_search.setTraceDump(dump_traces, {"late"});
    for (auto& s : scans(false)) { dump_search.consumeSpectrum(s); }
    dump_search.finish();
    std::istringstream dump_lines(dump_traces.str());
    std::string header;
    std::getline(dump_lines, header);
    require(header == "precursor_id\tscan\trt\tslot\tfragment_mz\tintensity\tobserved_mz", "trace dump header changed");
    std::size_t dumped = 0, matched = 0;
    double apex_intensity = 0.0, apex_observed = 0.0;
    for (std::string line; std::getline(dump_lines, line);)
    {
      std::istringstream fields(line);
      std::string id, observed;
      double scan_id = 0.0, rt = 0.0, slot = 0.0, fragment = 0.0, intensity = 0.0;
      fields >> id >> scan_id >> rt >> slot >> fragment >> intensity >> observed;
      require(id == "late" && rt >= 280.0 && rt <= 520.0, "trace dump row of another precursor or outside the RT window");
      ++dumped;
      if (intensity > 0.0)
      {
        ++matched;
        require(observed != "nan", "matched trace point without an observed m/z");
      }
      if (slot == 0.0 && intensity > apex_intensity) { apex_intensity = intensity; apex_observed = std::stod(observed); }
    }
    // 161 scans in [280, 520] times 2 fragments, the zero filled ones included.
    require(dumped == 2 * 161 && matched > 10 && matched < dumped, "trace dump should hold every scan of the RT window");
    require(std::abs(apex_intensity - 2000.0) < 1.0 && std::abs(apex_observed - 300.001) < 1e-9, "trace dump apex value or m/z is off");
    const auto dump_rows = parse(dump_output.str());
    require(dump_rows.size() == 1 && dump_rows[0].id == "late", "trace dump changed the peak output");
  }

#ifdef _OPENMP
  // More than one chunk per lane, several lanes and several batches: output must not depend on the thread count.
  {
    std::vector<std::string> ids;
    for (std::size_t c = 0; c < 900; ++c) { ids.push_back("p" + std::to_string(c)); }
    std::vector<OpenDIA::TransitionQuery> queries;
    for (std::size_t c = 0; c < ids.size(); ++c)
    {
      const double precursor_mz = (c % 3 == 0 ? 450.0 : 650.0) + static_cast<double>(c) * 0.01;
      const double library_rt = 20.0 + static_cast<double>(c % 97);
      for (std::size_t t = 0; t < 3; ++t)
      {
        queries.push_back({ids[c], precursor_mz, 100.0 + static_cast<double>(c % 256) * 0.02 + static_cast<double>(t) * 0.005, library_rt, c % 7 != 0});
      }
    }
    std::vector<MSSpectrum> spectra;
    for (std::size_t scan = 0; scan < 240; ++scan)
    {
      const double rt = 1.0 + 0.5 * static_cast<double>(scan);
      std::vector<Peak1D> peaks;
      for (std::size_t c = 0; c < 256; ++c)
      {
        const double fragment_mz = 100.0 + static_cast<double>(c) * 0.02;
        const double height = 100.0 + static_cast<double>(c);
        for (std::size_t t = 0; t < 3; ++t) { peaks.emplace_back(fragment_mz + static_cast<double>(t) * 0.005, bell(rt, 20.0 + static_cast<double>(c % 97), height) + 1.0); }
      }
      spectra.push_back(spectrum(rt, scan % 2 == 0 ? 500.0 : 700.0, peaks));
    }
    OpenDIA::SearchOptions options = defaults();
    options.rt_window_seconds = 15.0;
    const int original = omp_get_max_threads();
    std::string serial_features, parallel_features;
    omp_set_num_threads(1);
    const std::string serial = run(queries, options, spectra, &serial_features);
    omp_set_num_threads(4);
    const std::string parallel = run(queries, options, spectra, &parallel_features);
    omp_set_num_threads(original);
    require(serial == parallel && serial_features == parallel_features, "output changed with the thread count");
    require(parse(serial).size() >= ids.size(), "expected at least one peak per precursor");
  }
#endif
  return 0;
}
