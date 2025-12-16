// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace OpenMS;
using namespace OpenMS::Math;

START_TEST(MultipleTesting, "$Id$")

// loosen tolerances to match the reference comparison used in pyprophet (4 decimal places)
TOLERANCE_ABSOLUTE(1e-4)
TOLERANCE_RELATIVE(1.0 + 1e-4)

START_SECTION(template<class T> std::vector<double> compute_model_fdr(const std::vector<T>&))
{
  std::vector<double> data = {0.1, 0.2, 0.3};
  auto f = compute_model_fdr<double>(data);
  TEST_EQUAL(f.size(), data.size())
  TEST_REAL_SIMILAR(f[0], 0.1)
  TEST_REAL_SIMILAR(f[1], 0.15)
  TEST_REAL_SIMILAR(f[2], 0.2)
}
END_SECTION

START_SECTION("lfdr CSV reference (pyprophet) regression")
{
  // load CSV
  const std::string path = OPENMS_GET_TEST_DATA_PATH("test_lfdr_ref_data.csv");
  std::ifstream in(path.c_str());
  TEST_EQUAL(bool(in), true)

  std::string header;
  std::getline(in, header);
  struct Row { double p; double lfdr_default; double lfdr_monotone_false; double lfdr_transf_logit; double lfdr_eps; };
  std::vector<Row> rows;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ','))
    {
      // strip optional quotes
      if (!item.empty() && item.front() == '"' && item.back() == '"') item = item.substr(1, item.size()-2);
      parts.push_back(item);
    }
    if (parts.size() < 5) continue;
    Row r;
    r.p = std::stod(parts[0]);
    r.lfdr_default = std::stod(parts[1]);
    r.lfdr_monotone_false = std::stod(parts[2]);
    r.lfdr_transf_logit = std::stod(parts[3]);
    r.lfdr_eps = std::stod(parts[4]);
    rows.push_back(r);
  }

  TEST_EQUAL(rows.empty(), false)

  // sort by p ascending
  std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b){ return a.p < b.p; });

  std::vector<double> pvec, ref_default, ref_monotone_false, ref_logit, ref_eps;
  pvec.reserve(rows.size()); ref_default.reserve(rows.size()); ref_monotone_false.reserve(rows.size()); ref_logit.reserve(rows.size()); ref_eps.reserve(rows.size());
  for (auto &r : rows) { pvec.push_back(r.p); ref_default.push_back(r.lfdr_default); ref_monotone_false.push_back(r.lfdr_monotone_false); ref_logit.push_back(r.lfdr_transf_logit); ref_eps.push_back(r.lfdr_eps); }

  const double pi0 = 0.669926026474838;
  auto out_def = lfdr(pvec, pi0);
  auto out_monofalse = lfdr(pvec, pi0, true, false);
  auto out_logit = lfdr(pvec, pi0, true, true, "logit");
  auto out_eps = lfdr(pvec, pi0, true, true, "probit", 1.5, std::pow(10.0, -2));

  TEST_EQUAL(out_def.size(), ref_default.size())
  TEST_EQUAL(out_monofalse.size(), ref_monotone_false.size())
  TEST_EQUAL(out_logit.size(), ref_logit.size())
  TEST_EQUAL(out_eps.size(), ref_eps.size())

  const double tol = 1e-2; // match python test decimal=2
  for (std::size_t i = 0; i < ref_default.size(); ++i)
  {
    TEST_EQUAL(std::fabs(out_def[i] - ref_default[i]) <= tol, true);
    TEST_EQUAL(std::fabs(out_monofalse[i] - ref_monotone_false[i]) <= tol, true);
    TEST_EQUAL(std::fabs(out_logit[i] - ref_logit[i]) <= tol, true);
    TEST_EQUAL(std::fabs(out_eps[i] - ref_eps[i]) <= tol, true);
  }

  // print results in three lines for regression capture
  auto print_vec = [&](const std::vector<double>& v)
  {
    std::cout << "[";
    for (std::size_t i = 0; i < v.size(); ++i)
    {
      if (i) std::cout << ", ";
      std::cout << v[i];
    }
    std::cout << "]\n";
  };

  print_vec(out_def);
  print_vec(out_monofalse);
  print_vec(out_logit);
  print_vec(out_eps);
}
END_SECTION

START_SECTION("lfdr basic checks (probit & logit)")
{
  std::vector<double> p = {0.001, 0.01, 0.05, 0.2, 0.8, 0.95};
  double pi0 = 0.8;

  auto out_probit = lfdr(p, pi0, true, true, "probit", 1.5, 1e-8, 256, 3.0);
  TEST_EQUAL(out_probit.size(), p.size())
  // values in [0,1] and finite
  for (std::size_t i = 0; i < out_probit.size(); ++i)
  {
    TEST_EQUAL(std::isfinite(out_probit[i]), true);
    TEST_EQUAL(out_probit[i] >= 0.0, true);
    TEST_EQUAL(out_probit[i] <= 1.0, true);
  }
  // monotonic non-decreasing in p
  std::vector<std::size_t> idx(p.size()); for (std::size_t i = 0; i < p.size(); ++i) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b){ return p[a] < p[b]; });
  for (std::size_t i = 1; i < idx.size(); ++i) TEST_EQUAL(out_probit[idx[i]] >= out_probit[idx[i-1]], true);

  auto out_logit = lfdr(p, pi0, true, true, "logit", 1.5, 1e-8, 256, 3.0);
  TEST_EQUAL(out_logit.size(), p.size())
  for (double v : out_logit) { TEST_EQUAL(std::isfinite(v), true); TEST_EQUAL(v >= 0.0, true); TEST_EQUAL(v <= 1.0, true); }
}
END_SECTION

START_SECTION("forrt/revrt roundtrip small")
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
  auto Xp = forrt(x, x.size());
  auto xr = revrt(Xp, x.size());
  TEST_EQUAL(xr.size(), x.size())
  for (std::size_t i = 0; i < x.size(); ++i)
  {
    TEST_REAL_SIMILAR(xr[i], x[i]);
  }
}
END_SECTION

START_SECTION(double bw_nrd0(const std::vector<double>&))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  double bw = bw_nrd0(x);
  // Expected value matching numpy.percentile + statsmodels/pyprophet convention (0.9 factor)
  TEST_REAL_SIMILAR(bw, 0.9735846228506357)
}
END_SECTION

START_SECTION(std::vector<double> linbin(const std::vector<double>&, double, double, std::size_t, const std::vector<double>*))
{
  std::vector<double> x = {0.1, 0.4, 0.9};
  auto bins = linbin(x, 0.0, 1.0, 5, nullptr);
  TEST_EQUAL(bins.size(), 5)
  // width = 0.2 -> indices: 0,2,4
  TEST_REAL_SIMILAR(bins[0], 1.0)
  TEST_REAL_SIMILAR(bins[1], 0.0)
  TEST_REAL_SIMILAR(bins[2], 1.0)
  TEST_REAL_SIMILAR(bins[3], 0.0)
  TEST_REAL_SIMILAR(bins[4], 1.0)

  std::vector<double> w = {2.0, 3.0, 5.0};
  auto bins_w = linbin(x, 0.0, 1.0, 5, &w);
  TEST_REAL_SIMILAR(bins_w[0], 2.0)
  TEST_REAL_SIMILAR(bins_w[2], 3.0)
  TEST_REAL_SIMILAR(bins_w[4], 5.0)
}
END_SECTION

START_SECTION(std::vector<double> pnorm(const std::vector<double>&, const std::vector<double>&))
{
  std::vector<double> stat0 = {0.0, 1.0}; // mu=0.5, sample sd = sqrt(0.5)
  std::vector<double> stat = {0.5, 1.5, -0.5, std::numeric_limits<double>::quiet_NaN()};

  auto out = pnorm(stat, stat0);
  TEST_EQUAL(out.size(), stat.size());

  // compute expected values
  double mu = 0.5;
  double var = 0.5; // (0-0.5)^2 + (1-0.5)^2 = 0.5, ddof=1 => /1 = 0.5
  double sigma = std::sqrt(var);
  const double sqrt2 = std::sqrt(2.0);

  // first value equals mu -> tail = 0.5
  double z0 = (stat[0] - mu) / sigma;
  double expected0 = 1.0 - 0.5 * (1.0 + std::erf(z0 / sqrt2));
  TEST_REAL_SIMILAR(out[0], expected0);

  double z1 = (stat[1] - mu) / sigma;
  double expected1 = 1.0 - 0.5 * (1.0 + std::erf(z1 / sqrt2));
  TEST_REAL_SIMILAR(out[1], expected1);

  double z2 = (stat[2] - mu) / sigma;
  double expected2 = 1.0 - 0.5 * (1.0 + std::erf(z2 / sqrt2));
  TEST_REAL_SIMILAR(out[2], expected2);

  // NaN propagated
  TEST_EQUAL(std::isnan(out[3]), true);
}
END_SECTION

START_SECTION(template<class T> std::vector<double> pemp(const std::vector<T>&, const std::vector<T>&))
{
  std::vector<double> stat = {0.9, 0.8};
  std::vector<double> stat0 = {0.1, 0.2};
  auto p = pemp<double>(stat, stat0);
  TEST_EQUAL(p.size(), stat.size())
  TEST_REAL_SIMILAR(p[0], 0.5)
  TEST_REAL_SIMILAR(p[1], 0.5)
}
END_SECTION

START_SECTION(std::vector<double> qvalue(const std::vector<double>&, double, bool))
{
  std::vector<double> p = {0.01, 0.02, 0.03};
  double pi0 = 1.0;
  auto q = qvalue(p, pi0, false);
  TEST_EQUAL(q.size(), p.size())
  TEST_REAL_SIMILAR(q[0], 0.03)
  TEST_REAL_SIMILAR(q[1], 0.03)
  TEST_REAL_SIMILAR(q[2], 0.03)
}
END_SECTION

START_SECTION(Pi0Result pi0est(const std::vector<double>&, const std::vector<double>&))
{
  std::vector<double> p = {0.9, 0.8, 0.85, 0.95};
  // Choose lambda grid such that all p >= lambda -> estimator becomes clipped to 1
  std::vector<double> lambda = {0.5, 0.6};
  auto res = pi0est(p, lambda);
  TEST_EQUAL(res.pi0, 1.0)
}
END_SECTION

START_SECTION("R reference qvalue CSV matches C++ implementation")
{
  // load CSV
  const std::string path = OPENMS_GET_TEST_DATA_PATH("test_qvalue_ref_data.csv");
  std::ifstream in(path.c_str());
  TEST_EQUAL(bool(in), true)

  std::string header;
  std::getline(in, header);
  // read rows
  struct Row { double p; double qd; double qpfdr; };
  std::vector<Row> rows;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    // simple CSV split by comma, remove optional quotes
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ','))
    {
      // strip quotes
      if (!item.empty() && item.front() == '"' && item.back() == '"') item = item.substr(1, item.size()-2);
      parts.push_back(item);
    }
    if (parts.size() < 3) continue;
    Row r;
    r.p = std::stod(parts[0]);
    r.qd = std::stod(parts[1]);
    r.qpfdr = std::stod(parts[2]);
    rows.push_back(r);
  }

  TEST_EQUAL(rows.empty(), false)

  // sort by p ascending to match pyprophet behavior
  std::sort(rows.begin(), rows.end(), [](const Row &a, const Row &b){ return a.p < b.p; });

  std::vector<double> pvec, qd_ref, qpfdr_ref;
  pvec.reserve(rows.size()); qd_ref.reserve(rows.size()); qpfdr_ref.reserve(rows.size());
  for (auto &r : rows) { pvec.push_back(r.p); qd_ref.push_back(r.qd); qpfdr_ref.push_back(r.qpfdr); }

  const double pi0 = 0.669926026474838;
  auto q_def = qvalue(pvec, pi0, false);
  auto q_pf = qvalue(pvec, pi0, true);

  TEST_EQUAL(q_def.size(), qd_ref.size())
  TEST_EQUAL(q_pf.size(), qpfdr_ref.size())

  for (std::size_t i = 0; i < q_def.size(); ++i)
  {
    TEST_REAL_SIMILAR(q_def[i], qd_ref[i]);
    TEST_REAL_SIMILAR(q_pf[i], qpfdr_ref[i]);
  }
}
END_SECTION

START_SECTION("silverman_kernel_fft properties")
{
  // basic sanity checks
  std::size_t M = 8;
  double bw = 1.0;
  double RANGE = 10.0;
  auto K = silverman_kernel_fft(bw, M, RANGE);
  TEST_EQUAL(K.size(), M)
  // K[0] should equal 1 (exp(0)/BC with BC==1)
  TEST_REAL_SIMILAR(K[0], 1.0)
  // positive and decreasing for first half
  std::size_t half = M / 2;
  for (std::size_t k = 0; k <= half; ++k)
  {
    TEST_EQUAL(std::isfinite(K[k]), true)
    TEST_EQUAL(K[k] > 0.0, true)
    if (k > 0) TEST_EQUAL(K[k] <= K[k-1], true);
  }
  // symmetry in Munro-packed layout: K[1..half-1] mirrored at K[half+1..]
  for (std::size_t k = 1; k < half; ++k)
  {
    TEST_REAL_SIMILAR(K[k], K[half + k]);
  }
}
END_SECTION

START_SECTION("grid_kde_fft basic integration and sizes")
{
  std::vector<double> x = {0.0, 1.0, 2.0};
  double bw = bw_nrd0(x);
  auto res = grid_kde_fft(x, bw, 64, 3.0);
  auto dens = res.first;
  auto grid = res.second;
  TEST_EQUAL(dens.size(), grid.size())
  TEST_EQUAL(dens.empty(), false)
  // integral should be ~1.0
  double delta = grid.size() > 1 ? (grid[1] - grid[0]) : 1.0;
  double sum = 0.0;
  for (double v : dens) { TEST_EQUAL(std::isfinite(v), true); TEST_EQUAL(v >= 0.0, true); sum += v; }
  TEST_REAL_SIMILAR(sum * delta, 1.0)
}
END_SECTION

START_SECTION("kde_fft_eval matches spline-interpolated grid density at sample points")
{
  std::vector<double> x = {0.0, 1.0, 2.0};
  double bw = bw_nrd0(x);
  auto res = grid_kde_fft(x, bw, 64, 3.0);
  auto dens = res.first;
  auto grid = res.second;

  // expected values: spline(grid, dens) evaluated at original sample points x
  OpenMS::CubicSpline2d spline(grid, dens);
  std::vector<double> expected;
  for (double xi : x) expected.push_back(spline.eval(xi));

  // evaluate KDE at sample points using kde_fft_eval (x used as both data and query)
  auto y = kde_fft_eval(x, bw, 64, 3.0);
  TEST_EQUAL(y.size(), expected.size())
  for (std::size_t i = 0; i < expected.size(); ++i) TEST_REAL_SIMILAR(y[i], expected[i]);
}
END_SECTION

END_TEST
