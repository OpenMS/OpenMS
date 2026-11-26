// FDR-20261125 A/B replay harness (diagnosis only). Rebuilds the gbmRescore input from a run's features.tsv
// (feature columns, fold, precursor/pair from the library TRAML keys, pair_hash = fnv1a("v" + key)) and re-runs the tool's
// own gbmRescore with optional variants. Default settings must reproduce the exported score column exactly.
// Usage: gbm_replay <features.tsv> <keys.tsv: id \t traml \t decoy> <seed> <threads> <out.tsv> [exclude=a,b] [mono=f:+1,g:-1] [fold_local]
#include <OpenDIAGBM.h>

#include <charconv>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenMS::OpenDIA;

static std::vector<std::string_view> split(std::string_view s, char c)
{
  std::vector<std::string_view> out;
  std::size_t a = 0;
  for (std::size_t i = 0; i <= s.size(); ++i)
    if (i == s.size() || s[i] == c) { out.push_back(s.substr(a, i - a)); a = i + 1; }
  return out;
}

static double num(std::string_view s)
{
  double v = 0.0;
  if (s == "nan" || s == "-nan") return std::numeric_limits<double>::quiet_NaN();
  if (s == "inf") return std::numeric_limits<double>::infinity();
  if (s == "-inf") return -std::numeric_limits<double>::infinity();
  auto r = std::from_chars(s.data(), s.data() + s.size(), v);
  if (r.ec != std::errc()) throw std::runtime_error("bad number: " + std::string(s));
  return v;
}

int main(int argc, char** argv)
{
  if (argc < 6) { std::cerr << "usage\n"; return 2; }
  const std::string feat_path = argv[1], key_path = argv[2], out_path = argv[5];
  const std::uint64_t seed = std::stoull(argv[3]);
  const int threads = std::stoi(argv[4]);
  GBMVariant variant;
  std::vector<std::string> exclude_names;
  std::vector<std::pair<std::string, int>> mono_names;
  for (int a = 6; a < argc; ++a)
  {
    std::string arg = argv[a];
    if (arg.rfind("exclude=", 0) == 0)
    {
      const std::string list = arg.substr(8);  // split returns views: keep the string alive
      for (auto v : split(list, ',')) exclude_names.emplace_back(v);
    }
    else if (arg.rfind("mono=", 0) == 0)
    {
      const std::string list = arg.substr(5);
      for (auto v : split(list, ','))
      {
        const std::string item(v);
        auto kv = split(item, ':');
        mono_names.emplace_back(std::string(kv[0]), std::stoi(std::string(kv[1])));
      }
    }
    else if (arg == "fold_local") variant.fold_local = true;
    else { std::cerr << "bad arg " << arg << "\n"; return 2; }
  }
  std::unordered_map<std::string, std::string> key_of;
  {
    std::ifstream k(key_path);
    std::string line;
    while (std::getline(k, line))
    {
      auto f = split(line, '\t');
      std::string t(f[1]);
      if (t.rfind("DECOY_", 0) == 0) t = t.substr(6);
      key_of.emplace(std::string(f[0]), t);
    }
  }
  std::ifstream in(feat_path);
  std::string header;
  std::getline(in, header);
  auto cols = split(header, '\t');
  std::size_t c_score = 0, c_net = 0, c_fold = 3;
  std::vector<std::size_t> fcols;
  std::vector<std::string> fnames;
  for (std::size_t c = 0; c < cols.size(); ++c)
  {
    if (cols[c] == "score") c_score = c;
    else if (cols[c] == "network_score") c_net = c;
    else if (c >= 5 && cols[c] != "q_value" && cols[c] != "network_q_value") { fcols.push_back(c); fnames.emplace_back(cols[c]); }
  }
  const std::size_t d = fcols.size();
  variant.use_feature.assign(d, 1);
  variant.monotone.assign(d, 0);
  for (const auto& e : exclude_names)
  {
    bool found = false;
    for (std::size_t j = 0; j < d; ++j) if (fnames[j] == e) { variant.use_feature[j] = 0; found = true; }
    if (!found) { std::cerr << "no feature " << e << "\n"; return 2; }
  }
  for (const auto& [e, s] : mono_names)
  {
    bool found = false;
    for (std::size_t j = 0; j < d; ++j) if (fnames[j] == e) { variant.monotone[j] = static_cast<std::int8_t>(s); found = true; }
    if (!found) { std::cerr << "no feature " << e << "\n"; return 2; }
  }
  const bool fold_local = variant.fold_local;
  if (exclude_names.empty()) variant.use_feature.clear();
  if (mono_names.empty()) variant.monotone.clear();
  variant.fold_local = fold_local;
  CrossFitScorer::Input input;
  input.dimension = d;
  std::vector<float> net, tool;
  std::unordered_map<std::string, std::uint32_t> dprec, dpair;
  std::vector<std::string> ids;  // library precursor id per dense precursor index
  std::string line;
  std::size_t rows = 0;
  while (std::getline(in, line))
  {
    auto f = split(line, '\t');
    const std::string id(f[1]);
    const auto& key = key_of.at(id);
    const std::uint8_t decoy = f[2] == "1";
    auto ip = dprec.try_emplace(id, static_cast<std::uint32_t>(dprec.size())).first->second;
    auto ia = dpair.try_emplace(key, static_cast<std::uint32_t>(dpair.size())).first->second;
    if (ip == ids.size()) ids.push_back(id);
    input.precursor.push_back(ip);
    input.pair.push_back(ia);
    input.decoy.push_back(decoy);
    input.fold.push_back(static_cast<std::uint8_t>(num(f[c_fold])));
    input.pair_hash.push_back(CrossFitScorer::fnv1a("v" + key));
    for (std::size_t c : fcols) input.features.push_back(static_cast<float>(num(f[c])));
    net.push_back(static_cast<float>(num(f[c_net])));
    tool.push_back(static_cast<float>(num(f[c_score])));
    ++rows;
  }
  auto log = [](const std::string& l) { std::cerr << l << '\n'; };
  const auto res = gbmRescore(input, net, 0.02, seed, threads, log, 3, 0.01, variant);
  std::size_t same = 0, kept = 0;
  for (std::size_t r = 0; r < rows; ++r)
  {
    same += res.score[r] == tool[r];
    kept += res.score[r] > -1e29f;
  }
  std::cerr << "rows " << rows << " d " << d << " kept " << kept << " rows equal to exported score " << same << '\n';
  // Kept rows only (one window per precursor), so the result can be fetched cheaply: the replay scores and, as a
  // reproduction check, the run's own exported scores.
  char buf[64];
  auto dump = [&](const std::string& path, const std::vector<float>& sc)
  {
    std::ofstream out(path);
    out << "row\tprecursor_id\tis_decoy\tscore\n";
    for (std::size_t r = 0; r < rows; ++r)
    {
      if (!(sc[r] > -1e29f)) continue;
      auto e = std::to_chars(buf, buf + sizeof(buf), static_cast<double>(sc[r]));
      out << r << '\t' << ids[input.precursor[r]] << '\t' << int(input.decoy[r]) << '\t' << std::string_view(buf, e.ptr - buf) << '\n';
    }
  };
  dump(out_path, res.score);
  dump(out_path + ".self.tsv", tool);
  return 0;
}
