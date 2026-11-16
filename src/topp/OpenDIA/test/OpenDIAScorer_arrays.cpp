// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

// Validation driver (not a ctest): runs OpenDIA::CrossFitScorer on the PLAN_GAP scoring arrays of one run, the inputs
// of train_gap.py + make_w.py + calib_tail.py, and writes the calibrated scores (scores_<arm>WT.npy layout).
// Usage: OpenDIAScorer_arrays <run dir> <out.npy> <seed> <threads> <feature .npy> [...]
// The run dir holds population.tsv and train_windows.tsv; the feature arrays are float32 [rows, k], concatenated.

#include <OpenDIAScorer.h>

#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenMS::OpenDIA;

namespace
{
  std::vector<float> readNpy(const std::string& path, std::size_t& rows, std::size_t& cols)
  {
    std::ifstream in(path, std::ios::binary);
    char magic[10];
    in.read(magic, 10);
    if (!in || std::memcmp(magic, "\x93NUMPY", 6) != 0) throw std::runtime_error("not an npy file: " + path);
    std::size_t header_len = static_cast<unsigned char>(magic[8]) | (static_cast<unsigned char>(magic[9]) << 8);
    if (magic[6] != 1)
    {
      char extra[2];
      in.read(extra, 2);
      header_len = (header_len) | (static_cast<std::size_t>(static_cast<unsigned char>(extra[0])) << 16) |
                   (static_cast<std::size_t>(static_cast<unsigned char>(extra[1])) << 24);
    }
    std::string header(header_len, ' ');
    in.read(header.data(), static_cast<std::streamsize>(header_len));
    if (header.find("'<f4'") == std::string::npos) throw std::runtime_error("expected float32: " + path);
    const bool fortran = header.find("'fortran_order': True") != std::string::npos;
    const auto open = header.find('('), close = header.find(')');
    std::string shape = header.substr(open + 1, close - open - 1);
    std::stringstream ss(shape);
    std::string a, b;
    std::getline(ss, a, ',');
    std::getline(ss, b, ',');
    rows = std::stoull(a);
    cols = b.find_first_not_of(' ') == std::string::npos ? 1 : std::stoull(b);
    std::vector<float> v(rows * cols);
    in.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
    if (!in) throw std::runtime_error("short npy file: " + path);
    if (!fortran || cols == 1) return v;
    std::vector<float> c(rows * cols);
    for (std::size_t j = 0; j < cols; ++j)
      for (std::size_t r = 0; r < rows; ++r) c[r * cols + j] = v[j * rows + r];
    return c;
  }

  void writeNpy(const std::string& path, const std::vector<float>& v)
  {
    std::ofstream out(path, std::ios::binary);
    std::string header = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + std::to_string(v.size()) + ",), }";
    header.append((16 - (10 + header.size() + 1) % 16) % 16, ' ');
    header += '\n';
    const char magic[] = {'\x93', 'N', 'U', 'M', 'P', 'Y', 1, 0, static_cast<char>(header.size() & 255), static_cast<char>(header.size() >> 8)};
    out.write(magic, sizeof(magic));
    out << header;
    out.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(float)));
  }

  std::vector<std::string> split(const std::string& line)
  {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string f;
    while (std::getline(ss, f, '\t')) out.push_back(f);
    return out;
  }

  std::size_t col(const std::vector<std::string>& h, const std::string& name)
  {
    for (std::size_t i = 0; i < h.size(); ++i)
      if (h[i] == name) return i;
    throw std::runtime_error("missing column " + name);
  }
}

int main(int argc, char** argv)
{
  if (argc < 6)
  {
    std::cerr << "usage: OpenDIAScorer_arrays <run dir> <out.npy> <seed> <threads> <feature .npy> [...]\n";
    return 1;
  }
  const auto t0 = std::chrono::steady_clock::now();
  const std::string dir = argv[1], out_path = argv[2];
  CrossFitScorer::Options options;
  options.seed = std::stoull(argv[3]);
  options.threads = std::stoi(argv[4]);

  CrossFitScorer::Input input;
  std::unordered_map<std::string, std::uint32_t> precursor_index, pair_index, row_of;
  {
    std::ifstream pop(dir + "/population.tsv");
    std::string line;
    std::getline(pop, line);
    const auto h = split(line);
    const auto c_sid = col(h, "sample_id"), c_pid = col(h, "precursor_id"), c_dec = col(h, "is_decoy"), c_fold = col(h, "fold"),
               c_key = col(h, "pair_key");
    while (std::getline(pop, line))
    {
      if (line.empty()) continue;
      const auto v = split(line);
      const auto r = static_cast<std::uint32_t>(input.precursor.size());
      row_of.emplace(v[c_sid], r);
      input.precursor.push_back(precursor_index.try_emplace(v[c_pid], static_cast<std::uint32_t>(precursor_index.size())).first->second);
      input.pair.push_back(pair_index.try_emplace(v[c_key], static_cast<std::uint32_t>(pair_index.size())).first->second);
      input.decoy.push_back(static_cast<std::uint8_t>(std::stoi(v[c_dec])));
      input.fold.push_back(static_cast<std::uint8_t>(std::stoi(v[c_fold])));
      input.pair_hash.push_back(CrossFitScorer::fnv1a("v" + v[c_key]));
    }
  }
  const std::size_t n = input.rows();
  {
    std::ifstream tw(dir + "/train_windows.tsv");
    std::string line;
    std::getline(tw, line);
    const auto h = split(line);
    const auto c_sid = col(h, "sample_id"), c_group = col(h, "group"), c_label = col(h, "label"), c_family = col(h, "family");
    while (std::getline(tw, line))
    {
      if (line.empty()) continue;
      const auto v = split(line);
      if (v[c_group] != "anchor" && v[c_group] != "anchor_decoy" && v[c_group] != "decoy") continue;
      input.seeds.push_back({row_of.at(v[c_sid]), v[c_label] == "1", CrossFitScorer::fnv1a("v" + v[c_family])});
    }
  }
  std::vector<std::vector<float>> blocks;
  std::vector<std::size_t> widths;
  for (int a = 5; a < argc; ++a)
  {
    std::size_t rows = 0, cols = 0;
    blocks.push_back(readNpy(argv[a], rows, cols));
    if (rows != n) throw std::runtime_error(std::string("row count mismatch: ") + argv[a]);
    widths.push_back(cols);
    input.dimension += cols;
  }
  input.features.resize(n * input.dimension);
  for (std::size_t r = 0; r < n; ++r)
  {
    std::size_t o = 0;
    for (std::size_t b = 0; b < blocks.size(); ++b)
    {
      for (std::size_t c = 0; c < widths[b]; ++c) input.features[r * input.dimension + o + c] = blocks[b][r * widths[b] + c];
      o += widths[b];
    }
  }
  blocks.clear();
  std::cerr << "rows " << n << " dimension " << input.dimension << " seeds " << input.seeds.size() << " precursors " << precursor_index.size()
            << " pairs " << pair_index.size() << '\n';
  const auto result = CrossFitScorer::run(input, options, [](const std::string& l) { std::cerr << l << '\n'; });
  writeNpy(out_path, result.score);
  // Optional diagnostics: every window's fold-model score and linear part (the W rule's input).
  if (std::getenv("OPENDIA_WRITE_RAW"))
  {
    writeNpy(out_path + ".raw.npy", result.raw_score);
    writeNpy(out_path + ".linear.npy", result.linear);
  }
  std::size_t at1 = 0;
  for (std::size_t r = 0; r < n; ++r) at1 += !input.decoy[r] && result.q_value[r] <= 0.01;
  std::cerr << "targets at q <= 0.01: " << at1 << ", seconds "
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << '\n';
  return 0;
}
