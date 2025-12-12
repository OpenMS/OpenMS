// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

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
