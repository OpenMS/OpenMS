// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIACalibration.h"
#include "OpenDIAInterference.h"

#include "OpenDIALibrary.h"

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <unordered_map>

using namespace OpenMS;

/**
  @page TOPP_OpenDIAInterference OpenDIAInterference
  @brief Shared-fragment interference features for OpenDIA candidate peaks (explaining away between co-eluting candidates).

  For every candidate window (one peak of one library precursor, from OpenDIA -features_out), the up to 12 library
  fragments are extracted from the raw MS2 spectra that isolate the precursor, on a 32-bin grid over apex +- 15 s. Windows
  of different precursors in one isolation window whose peaks overlap and share fragment m/z are linked; each claimant's
  abundance (from its own unshared fragments) and elution profile predict how much of the shared signal it explains, and
  that part is removed. The output gives, per window, the explained share of its central signal, the library cosine of
  the residual, its change, and the number of claimants (see OpenDIA::InterferenceResolver). Targets and decoys are
  treated identically; no labels are used.

  Input windows: a TSV with sample_id (= candidate_row of -candidates) and precursor_id columns, e.g. the scoring
  population. Output rows follow the input order.
*/
class TOPPOpenDIAInterference final : public TOPPBase
{
public:
  TOPPOpenDIAInterference() : TOPPBase("OpenDIAInterference", "Shared-fragment interference features for OpenDIA candidate peaks.", {}, false) {}

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "DIA raw or mzML input");
    registerInputFile_("lib", "<file>", "", "PQP library with up to 12 fragments per precursor (precursor IDs as in -windows)");
    setValidFormats_("lib", {"pqp"});
    registerInputFile_("calibration", "<file>", "", "OpenDIACalibrate JSON: fragment m/z offset and fitted fragment tolerance");
    setValidFormats_("calibration", {"json"});
    registerInputFile_("candidates", "<file>", "", "OpenDIA -features_out table (candidate_row, precursor_id, rt_start, rt_apex, rt_end)");
    setValidFormats_("candidates", {"tsv"});
    registerInputFile_("windows", "<file>", "", "Windows to score: TSV with sample_id (= candidate_row) and precursor_id");
    setValidFormats_("windows", {"tsv"});
    registerOutputFile_("out", "<file>", "", "Features per window: sample_id, precursor_id, i_explained, i_libcos_resid, i_dlibcos, i_nclaim");
    setValidFormats_("out", {"tsv"});
    registerOutputFile_("out_npy", "<file>", "", "Optional: the four features as float32 [windows, 4] NumPy array, input row order", false);
    registerDoubleOption_("extraction_ppm", "<ppm>", 17.5, "Half width of each fragment extraction around the calibrated m/z", false);
    setMinFloat_("extraction_ppm", 0.1);
    registerDoubleOption_("rt_half_width", "<s>", 15.0, "Half width of the 32-bin RT grid around each window apex", false);
    setMinFloat_("rt_half_width", 1.0);
    registerIntOption_("threads", "<n>", 1, "Threads for the interference graph", false);
    setMinInt_("threads", 1);
  }

  static std::vector<std::string> splitTabs_(const std::string& line)
  {
    std::vector<std::string> out;
    StringUtils::split(line, '\t', out);
    return out;
  }

  static std::size_t column_(const std::vector<std::string>& header, const std::string& name, const std::string& file)
  {
    const auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end()) throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Column " + name + " missing in " + file);
    return static_cast<std::size_t>(it - header.begin());
  }

  static double toDouble_(const std::string& s)
  {
    double v = 0.0;
    const auto r = std::from_chars(s.data(), s.data() + s.size(), v);
    if (r.ec != std::errc() || r.ptr != s.data() + s.size()) return StringUtils::toDouble(s);
    return v;
  }

  ExitCodes main_(int, const char**) override
  {
    const auto t0 = std::chrono::steady_clock::now();
    const std::string out_path = getStringOption_("out");
    const std::string npy_path = getStringOption_("out_npy");
    for (const auto& input : {"in", "lib", "calibration", "candidates", "windows"})
    {
      if (getStringOption_(input) == out_path || (!npy_path.empty() && getStringOption_(input) == npy_path))
      {
        writeLogError_("Outputs must not overwrite an input.");
        return ILLEGAL_PARAMETERS;
      }
    }
    const auto calibration = OpenDIA::readCalibrationRT(getStringOption_("calibration"));
    if (!calibration.has_ms2 || !(calibration.ms2_tolerance_ppm > 0.0))
    {
      writeLogError_("The calibration has no fragment mass model.");
      return ILLEGAL_PARAMETERS;
    }
    const double half = getDoubleOption_("rt_half_width");
    const double ppm = getDoubleOption_("extraction_ppm");
    const int threads = getIntOption_("threads");

    // Library: precursor m/z and the up to 12 fragments per precursor in transition ID order (the slot order of build12.py).
    const OpenDIA::Library library = OpenDIA::readLibrary(getStringOption_("lib"));
    struct LibraryFragment
    {
      long id;
      double mz;
      double intensity;
    };
    std::unordered_map<std::string, std::vector<LibraryFragment>> fragments;
    std::unordered_map<std::string, double> library_mz;
    fragments.reserve(library.getCompounds().size());
    for (const auto& t : library.getTransitions())
    {
      fragments[t.peptide_ref].push_back({std::stol(t.transition_name), t.product_mz, t.library_intensity});
      library_mz[t.peptide_ref] = t.precursor_mz;
    }

    // Candidate peaks.
    const std::string candidates_path = getStringOption_("candidates");
    std::ifstream candidates(candidates_path);
    std::string line;
    if (!candidates || !std::getline(candidates, line))
    {
      writeLogError_("Could not read -candidates: " + candidates_path);
      return INPUT_FILE_NOT_READABLE;
    }
    const auto ch = splitTabs_(line);
    const std::size_t c_row = column_(ch, "candidate_row", candidates_path), c_start = column_(ch, "rt_start", candidates_path),
                      c_apex = column_(ch, "rt_apex", candidates_path), c_end = column_(ch, "rt_end", candidates_path),
                      c_pid = column_(ch, "precursor_id", candidates_path);
    struct Peak
    {
      std::string precursor_id;
      double start, apex, end;
    };
    std::unordered_map<std::string, Peak> peaks;
    while (std::getline(candidates, line))
    {
      if (line.empty()) continue;
      const auto v = splitTabs_(line);
      peaks.emplace(v[c_row], Peak{v[c_pid], toDouble_(v[c_start]), toDouble_(v[c_apex]), toDouble_(v[c_end])});
    }

    // Windows, in input order.
    const std::string windows_path = getStringOption_("windows");
    std::ifstream windows_in(windows_path);
    if (!windows_in || !std::getline(windows_in, line))
    {
      writeLogError_("Could not read -windows: " + windows_path);
      return INPUT_FILE_NOT_READABLE;
    }
    const auto wh = splitTabs_(line);
    const std::size_t w_sid = column_(wh, "sample_id", windows_path), w_pid = column_(wh, "precursor_id", windows_path);
    std::vector<std::string> sample_ids, precursor_ids;
    std::unordered_map<std::string, std::uint32_t> precursor_index;
    std::vector<OpenDIA::InterferenceResolver::Precursor> precursors;
    std::vector<double> precursor_mz;
    std::vector<std::array<double, OpenDIA::InterferenceResolver::FRAGMENTS>> extraction_mz;
    std::vector<OpenDIA::InterferenceResolver::Window> windows;
    constexpr std::size_t K = OpenDIA::InterferenceResolver::FRAGMENTS;
    while (std::getline(windows_in, line))
    {
      if (line.empty()) continue;
      const auto v = splitTabs_(line);
      const std::string& sid = v[w_sid];
      const std::string& pid = v[w_pid];
      const auto peak = peaks.find(sid);
      if (peak == peaks.end() || peak->second.precursor_id != pid)
      {
        writeLogError_("Window " + sid + " has no candidate row with precursor " + pid + " in -candidates.");
        return INCOMPATIBLE_INPUT_DATA;
      }
      auto [it, inserted] = precursor_index.try_emplace(pid, static_cast<std::uint32_t>(precursors.size()));
      if (inserted)
      {
        auto frag = fragments.find(pid);
        if (frag == fragments.end())
        {
          writeLogError_("Precursor " + pid + " is not in -lib.");
          return INCOMPATIBLE_INPUT_DATA;
        }
        auto& f = frag->second;
        std::sort(f.begin(), f.end(), [](const LibraryFragment& a, const LibraryFragment& b) { return a.id < b.id; });
        OpenDIA::InterferenceResolver::Precursor p;
        std::array<double, K> extraction{};
        double maximum = 0.0;
        for (std::size_t s = 0; s < std::min(K, f.size()); ++s) maximum = std::max(maximum, f[s].intensity);
        for (std::size_t s = 0; s < std::min(K, f.size()); ++s)
        {
          p.fragment_mz[s] = f[s].mz;
          p.library[s] = f[s].intensity / std::max(maximum, 1e-12);
          extraction[s] = f[s].mz * (1.0 + calibration.ms2_offset_ppm * 1e-6);
        }
        precursors.push_back(p);
        precursor_mz.push_back(library_mz.at(pid));
        extraction_mz.push_back(extraction);
      }
      OpenDIA::InterferenceResolver::Window w;
      w.precursor = it->second;
      w.rt_start = peak->second.start;
      w.rt_apex = peak->second.apex;
      w.rt_end = peak->second.end;
      windows.push_back(w);
      sample_ids.push_back(sid);
      precursor_ids.push_back(pid);
    }
    peaks.clear();
    fragments.clear();
    writeLogInfo_("Interference: windows=" + std::to_string(windows.size()) + " precursors=" + std::to_string(precursors.size()) +
                  " fragment_tolerance_ppm=" + std::to_string(calibration.ms2_tolerance_ppm) + " ms2_offset_ppm=" +
                  std::to_string(calibration.ms2_offset_ppm));

    // Extraction.
    OpenDIA::InterferenceCollector collector(precursor_mz, std::move(extraction_mz), windows, half, ppm);
    FileHandler reader;
    PeakFileOptions options;
    options.setMSLevels({2});
    reader.setOptions(options);
    reader.transformExperiment(getStringOption_("in"), &collector, {}, log_type_);
    const auto isolation = collector.isolationIndices();
    std::size_t isolated = 0;
    for (std::size_t p = 0; p < precursors.size(); ++p)
    {
      precursors[p].isolation = isolation[p];
      isolated += isolation[p] >= 0;
    }
    const auto t1 = std::chrono::steady_clock::now();
    writeLogInfo_("Interference extraction: spectra_used=" + std::to_string(collector.spectraUsed()) + " isolation_windows=" +
                  std::to_string(collector.isolationWindows().size()) + " precursors_isolated=" + std::to_string(isolated) +
                  " out_of_order_spectra=" + std::to_string(collector.outOfOrderSpectra()) + " seconds=" +
                  std::to_string(std::chrono::duration<double>(t1 - t0).count()));

    OpenDIA::InterferenceResolver::Statistics stats;
    const auto features = OpenDIA::InterferenceResolver::compute(precursors, windows, calibration.ms2_tolerance_ppm,
                                                                 2.0 * half / OpenDIA::InterferenceResolver::BINS, threads, &stats);
    const auto t2 = std::chrono::steady_clock::now();
    writeLogInfo_("Interference graph: pairs=" + std::to_string(stats.pairs_tested) + " links=" + std::to_string(stats.links) +
                  " used_links=" + std::to_string(stats.used_links) + " windows_with_claimant=" + std::to_string(stats.windows_with_claimant) +
                  " seconds=" + std::to_string(std::chrono::duration<double>(t2 - t1).count()));

    std::ofstream out(out_path);
    if (!out)
    {
      writeLogError_("Could not open -out: " + out_path);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    out << "sample_id\tprecursor_id";
    for (const char* name : OpenDIA::InterferenceResolver::featureNames()) out << '\t' << name;
    out << '\n';
    std::string row;
    char buffer[32];
    for (std::size_t i = 0; i < windows.size(); ++i)
    {
      row.assign(sample_ids[i]).append(1, '\t').append(precursor_ids[i]);
      for (std::size_t f = 0; f < OpenDIA::InterferenceResolver::FEATURES; ++f)
      {
        const auto r = std::to_chars(buffer, buffer + sizeof(buffer), features[i * OpenDIA::InterferenceResolver::FEATURES + f]);
        row.append(1, '\t').append(buffer, r.ptr);
      }
      row.append(1, '\n');
      out << row;
    }
    out.close();
    if (!out)
    {
      writeLogError_("Could not write -out: " + out_path);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    if (!npy_path.empty() && !writeNpy_(npy_path, features, windows.size(), OpenDIA::InterferenceResolver::FEATURES))
    {
      writeLogError_("Could not write -out_npy: " + npy_path);
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    return EXECUTION_OK;
  }

  /// float32 C-order NumPy v1.0 file of shape [rows, columns].
  static bool writeNpy_(const std::string& path, const std::vector<float>& values, std::size_t rows, std::size_t columns)
  {
    std::ofstream npy(path, std::ios::binary);
    std::string header = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + std::to_string(rows) + ", " + std::to_string(columns) + "), }";
    header.append((16 - (10 + header.size() + 1) % 16) % 16, ' ');
    header += '\n';
    const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0, static_cast<char>(header.size() & 255), static_cast<char>((header.size() >> 8) & 255)};
    npy.write(magic, sizeof(magic));
    npy << header;
    npy.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(float)));
    npy.close();
    return static_cast<bool>(npy);
  }
};

int main(int argc, const char** argv)
{
  TOPPOpenDIAInterference tool;
  return tool.main(argc, argv);
}
