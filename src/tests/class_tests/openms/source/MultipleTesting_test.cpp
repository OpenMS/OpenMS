// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
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

END_TEST
