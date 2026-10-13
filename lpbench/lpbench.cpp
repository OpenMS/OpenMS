// LP/MIP solver benchmark for OpenMS::LPWrapper (COIN-OR vs GLPK vs HiGHS).
// Build the same source against each OpenMS build; the backend is fixed at compile time.
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/config.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace OpenMS;
using Clock = std::chrono::steady_clock;

static const char* backend()
{
#if defined(OPENMS_HAS_COINOR)
  return "COIN";
#elif defined(OPENMS_HAS_HIGHS)
  return "HIGHS";
#else
  return "GLPK";
#endif
}

struct Result { double build_s, solve_s, obj; int status; };

// Decharging-like ILP: binary edges between features, each feature in at most one selected edge.
static Result edgeSelection(int nodes, int deg, unsigned seed)
{
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> pick(0, nodes - 1);
  std::uniform_real_distribution<double> w(0.1, 10.0);
  auto t0 = Clock::now();
  LPWrapper lp;
  lp.setObjectiveSense(LPWrapper::MAX);
  std::vector<std::vector<Int>> incident(nodes);
  int edges = nodes * deg / 2;
  for (int e = 0; e < edges; ++e)
  {
    int a = pick(rng), b = pick(rng);
    if (a == b) continue;
    Int c = lp.addColumn();
    lp.setColumnBounds(c, 0, 1, LPWrapper::DOUBLE_BOUNDED);
    lp.setColumnType(c, LPWrapper::INTEGER);
    lp.setObjective(c, w(rng));
    incident[a].push_back(c);
    incident[b].push_back(c);
  }
  for (int n = 0; n < nodes; ++n)
  {
    if (incident[n].size() < 2) continue;
    std::vector<double> ones(incident[n].size(), 1.0);
    lp.addRow(incident[n], ones, "n" + std::to_string(n), 0, 1, LPWrapper::DOUBLE_BOUNDED);
  }
  auto t1 = Clock::now();
  LPWrapper::SolverParam p;
  p.message_level = 0;
  p.time_limit = 120000;
  lp.solve(p);
  auto t2 = Clock::now();
  return {std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count(),
          lp.getObjectiveValue(), (int)lp.getStatus()};
}

// Generalized assignment problem: jobs -> agents with capacities (harder MIP).
static Result gap(int agents, int jobs, unsigned seed)
{
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> cost(10, 50), wt(5, 25);
  auto t0 = Clock::now();
  LPWrapper lp;
  lp.setObjectiveSense(LPWrapper::MIN);
  std::vector<std::vector<Int>> x(agents, std::vector<Int>(jobs));
  std::vector<std::vector<int>> wgt(agents, std::vector<int>(jobs));
  for (int i = 0; i < agents; ++i)
    for (int j = 0; j < jobs; ++j)
    {
      Int c = lp.addColumn();
      lp.setColumnBounds(c, 0, 1, LPWrapper::DOUBLE_BOUNDED);
      lp.setColumnType(c, LPWrapper::INTEGER);
      lp.setObjective(c, cost(rng));
      x[i][j] = c;
      wgt[i][j] = wt(rng);
    }
  for (int j = 0; j < jobs; ++j)
  {
    std::vector<Int> idx; std::vector<double> v;
    for (int i = 0; i < agents; ++i) { idx.push_back(x[i][j]); v.push_back(1.0); }
    lp.addRow(idx, v, "j" + std::to_string(j), 1, 1, LPWrapper::FIXED);
  }
  for (int i = 0; i < agents; ++i)
  {
    std::vector<Int> idx; std::vector<double> v; double tot = 0;
    for (int j = 0; j < jobs; ++j) { idx.push_back(x[i][j]); v.push_back(wgt[i][j]); tot += wgt[i][j]; }
    lp.addRow(idx, v, "a" + std::to_string(i), 0, 0.8 * tot / agents, LPWrapper::UPPER_BOUND_ONLY);
  }
  auto t1 = Clock::now();
  LPWrapper::SolverParam p;
  p.message_level = 0;
  p.time_limit = 120000;
  lp.solve(p);
  auto t2 = Clock::now();
  return {std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count(),
          lp.getObjectiveValue(), (int)lp.getStatus()};
}

// Continuous transportation LP.
static Result transport(int n, unsigned seed)
{
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> cost(1, 100), amt(10, 100);
  auto t0 = Clock::now();
  LPWrapper lp;
  lp.setObjectiveSense(LPWrapper::MIN);
  std::vector<double> supply(n), demand(n);
  double s = 0;
  for (auto& v : supply) { v = amt(rng); s += v; }
  for (auto& v : demand) v = s / n * 0.95;
  std::vector<std::vector<Int>> x(n, std::vector<Int>(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    {
      Int c = lp.addColumn();
      lp.setColumnBounds(c, 0, 0, LPWrapper::LOWER_BOUND_ONLY);
      lp.setColumnType(c, LPWrapper::CONTINUOUS);
      lp.setObjective(c, cost(rng));
      x[i][j] = c;
    }
  for (int i = 0; i < n; ++i)
  {
    std::vector<Int> idx(x[i]); std::vector<double> v(n, 1.0);
    lp.addRow(idx, v, "s" + std::to_string(i), 0, supply[i], LPWrapper::UPPER_BOUND_ONLY);
  }
  for (int j = 0; j < n; ++j)
  {
    std::vector<Int> idx; std::vector<double> v(n, 1.0);
    for (int i = 0; i < n; ++i) idx.push_back(x[i][j]);
    lp.addRow(idx, v, "d" + std::to_string(j), demand[j], 0, LPWrapper::LOWER_BOUND_ONLY);
  }
  auto t1 = Clock::now();
  LPWrapper::SolverParam p;
  p.message_level = 0;
  lp.solve(p);
  auto t2 = Clock::now();
  return {std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count(),
          lp.getObjectiveValue(), (int)lp.getStatus()};
}

int main()
{
  const int reps = 3;
  auto run = [&](const std::string& name, auto f)
  {
    std::vector<double> solve;
    Result last{};
    for (int r = 0; r < reps; ++r) { last = f(); solve.push_back(last.solve_s); }
    std::sort(solve.begin(), solve.end());
    std::printf("%s,%s,%.4f,%.4f,%.6f,%d\n", backend(), name.c_str(), last.build_s, solve[reps / 2], last.obj, last.status);
    std::fflush(stdout);
  };
  std::printf("backend,problem,build_s,solve_median_s,objective,status\n");
  for (int n : {1000, 5000, 20000, 50000})
    run("edgesel_n" + std::to_string(n), [&] { return edgeSelection(n, 6, 42 + n); });
  for (auto [a, j] : std::vector<std::pair<int, int>>{{5, 30}, {10, 60}, {15, 100}})
    run("gap_" + std::to_string(a) + "x" + std::to_string(j), [&] { return gap(a, j, 7 + a * j); });
  for (int n : {50, 150, 300})
    run("transport_" + std::to_string(n), [&] { return transport(n, 3 + n); });
}
