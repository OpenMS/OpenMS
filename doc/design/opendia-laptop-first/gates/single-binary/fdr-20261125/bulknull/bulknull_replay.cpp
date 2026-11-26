// FDR-20261126 bulk-null replay harness (dev screen only). Rebuilds the scorer input from a run's features.tsv, exactly as
// gbm_replay.cpp does: feature columns, fold, precursor and pair from the library TRAML keys, pair_hash = fnv1a("v" + key).
// It then takes the run's exported final `score` and calls the tool's bulkNullQ on it.
// Usage: bulknull_replay <features.tsv> <keys.tsv: id \t traml \t decoy ...> <threads> <out.tsv>
// The output has kept rows only: precursor_id, is_decoy, score, exported q_value, replayed plain TDC q, bulk-null q.
#include <OpenDIAGBM.h>

#include <charconv>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
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
  if (argc != 5) { std::cerr << "usage: bulknull_replay <features.tsv> <keys.tsv> <threads> <out.tsv>\n"; return 2; }
  const int threads = std::stoi(argv[3]);
  std::unordered_map<std::string, std::string> key_of;
  {
    std::ifstream k(argv[2]);
    std::string line;
    while (std::getline(k, line))
    {
      auto f = split(line, '\t');
      std::string t(f[1]);
      if (t.rfind("DECOY_", 0) == 0) t = t.substr(6);
      key_of.emplace(std::string(f[0]), t);
    }
  }
  std::ifstream in(argv[1]);
  std::string header;
  std::getline(in, header);
  auto cols = split(header, '\t');
  std::size_t c_score = 0, c_q = 0, c_fold = 3;
  std::vector<std::size_t> fcols;
  for (std::size_t c = 0; c < cols.size(); ++c)
  {
    if (cols[c] == "score") c_score = c;
    else if (cols[c] == "q_value") c_q = c;
    else if (c >= 5 && cols[c] != "network_score" && cols[c] != "network_q_value") fcols.push_back(c);
  }
  CrossFitScorer::Input input;
  input.dimension = fcols.size();
  std::vector<float> score;
  std::vector<double> tool_q;
  std::unordered_map<std::string, std::uint32_t> dprec, dpair;
  std::vector<std::string> ids;
  std::string line;
  while (std::getline(in, line))
  {
    auto f = split(line, '\t');
    const std::string id(f[1]);
    const auto& key = key_of.at(id);
    auto ip = dprec.try_emplace(id, static_cast<std::uint32_t>(dprec.size())).first->second;
    auto ia = dpair.try_emplace(key, static_cast<std::uint32_t>(dpair.size())).first->second;
    if (ip == ids.size()) ids.push_back(id);
    input.precursor.push_back(ip);
    input.pair.push_back(ia);
    input.decoy.push_back(f[2] == "1");
    input.fold.push_back(static_cast<std::uint8_t>(num(f[c_fold])));
    input.pair_hash.push_back(CrossFitScorer::fnv1a("v" + key));
    for (std::size_t c : fcols) input.features.push_back(static_cast<float>(num(f[c])));
    score.push_back(static_cast<float>(num(f[c_score])));
    tool_q.push_back(num(f[c_q]));
  }
  const std::size_t n = score.size();
  std::vector<std::uint32_t> kept;
  for (std::uint32_t r = 0; r < n; ++r)
    if (score[r] > -1e29f) kept.push_back(r);
  const auto comp = CrossFitScorer::pooledPairedTDC(input, score, kept);
  std::vector<double> plain(n, 1.0);
  for (std::size_t i = 0; i < comp.ranked.size(); ++i) plain[comp.ranked[i]] = comp.q[i];
  const auto bn = bulkNullQ(input, score, threads, [](const std::string& l) { std::cerr << l << '\n'; });
  std::size_t same = 0, pt = 0, bt = 0;
  for (std::uint32_t r : kept)
  {
    same += static_cast<float>(plain[r]) == static_cast<float>(tool_q[r]) || std::abs(plain[r] - tool_q[r]) < 1e-9;
    pt += !input.decoy[r] && plain[r] <= 0.01;
    bt += !input.decoy[r] && bn[r] <= 0.01;
  }
  std::cerr << "rows " << n << " d " << input.dimension << " kept " << kept.size() << " kept rows whose replayed plain q equals the exported q "
            << same << "; targets at q<=0.01 plain " << pt << " bulk-null " << bt << '\n';
  std::ofstream out(argv[4]);
  out << "precursor_id\tis_decoy\tscore\ttool_q\tplain_q\tbn_q\n";
  char buf[64];
  const auto put = [&](double v) { auto e = std::to_chars(buf, buf + sizeof(buf), v); out << '\t' << std::string_view(buf, e.ptr - buf); };
  for (std::uint32_t r : kept)
  {
    out << ids[input.precursor[r]] << '\t' << int(input.decoy[r]);
    put(static_cast<double>(score[r]));
    put(tool_q[r]);
    put(plain[r]);
    put(bn[r]);
    out << '\n';
  }
  return 0;
}
