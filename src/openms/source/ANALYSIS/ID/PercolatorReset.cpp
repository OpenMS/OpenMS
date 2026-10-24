// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PercolatorReset.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

// vendored Percolator: the TRON solver of the linear L2-SVM (wired PRIVATE to libOpenMS)
#include "ssl.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <exception>
#include <map>
#include <numeric>
#include <limits>
#include <random>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
namespace P = OpenMS::Internal::Percolator;

namespace
{
  // ---------------------------------------------------------------------------------------------
  // Hashing: the coin is a function of (salt, sequence); tie breaks and fold assignment of (salt, sequence, row
  // index). None depends on labels, threads or scores; tie breaks (and so results with tied scores) depend on the
  // order of the input rows.
  // ---------------------------------------------------------------------------------------------
  std::uint64_t splitmix64(std::uint64_t x)
  {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }

  std::uint64_t fnv1a64(const std::string& s)
  {
    std::uint64_t h = 0xCBF29CE484222325ULL;
    for (const unsigned char c : s)
    {
      h ^= c;
      h *= 0x100000001B3ULL;
    }
    return h;
  }

  // purposes (distinct streams of the same salt)
  constexpr std::uint64_t COIN = 0x434F494E00000001ULL;
  constexpr std::uint64_t TIE = 0x5449450000000002ULL;
  constexpr std::uint64_t FOLD = 0x464F4C4400000003ULL;
  constexpr std::uint64_t MLP = 0x4D4C500000000004ULL;
  constexpr std::uint64_t AUDIT = 0x4155440000000005ULL;

  std::uint64_t keyedHash(std::uint64_t salt, std::uint64_t purpose, std::uint64_t value)
  {
    return splitmix64(splitmix64(salt ^ purpose) ^ value);
  }

  // Ranks @p items (indices) by score, best first; ties by @p tie (ascending).
  void rankByScore(std::vector<Size>& items, const std::vector<double>& score, const std::vector<std::uint64_t>& tie)
  {
    std::sort(items.begin(), items.end(), [&](Size a, Size b)
    {
      if (score[a] != score[b]) return score[a] > score[b];
      if (tie[a] != tie[b]) return tie[a] < tie[b];
      return a < b;
    });
  }

  // ---------------------------------------------------------------------------------------------
  // Learners on standardized feature rows
  // ---------------------------------------------------------------------------------------------
  struct Matrix
  {
    Size rows = 0, cols = 0;
    std::vector<double> v; // row-major
    const double* row(Size i) const { return v.data() + i * cols; }
  };

  // ---------------------------------------------------------------------------------------------
  // MLP learner (experimental): ensemble of networks d -> h -> h -> 1 (tanh, tanh, linear) in float32, trained by
  // mini-batch Adam with the batch's forward and backward passes as small matrix products (Eigen). Every network is
  // trained by one thread from its own random stream, so results do not depend on the number of threads.
  // ---------------------------------------------------------------------------------------------
  using MatF = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using VecF = Eigen::Matrix<float, Eigen::Dynamic, 1>;
  using RowVecF = Eigen::Matrix<float, 1, Eigen::Dynamic>;

  struct MlpNet
  {
    Eigen::Index d = 0, h = 0;
    std::vector<float> p;      // W1 (d x h, row-major), b1 (h), W2 (h x h), b2 (h), w3 (h), b3
    std::vector<float> m1, m2; // Adam moments (kept for warm starts)
    std::uint64_t step = 0;    // Adam steps so far

    Eigen::Index oB1() const { return d * h; }
    Eigen::Index oW2() const { return d * h + h; }
    Eigen::Index oB2() const { return d * h + h + h * h; }
    Eigen::Index oW3() const { return d * h + 2 * h + h * h; }
    Eigen::Index oB3() const { return d * h + 3 * h + h * h; }
  };

  // Work buffers of the forward pass (per thread)
  struct MlpWork
  {
    MatF a1, a2;
    VecF o;
  };

  // Outputs of @p net for the @p n rows (row-major, net.d columns) at @p x, into work.o
  void mlpForward(const MlpNet& net, const float* x, Eigen::Index n, MlpWork& work)
  {
    const Eigen::Index d = net.d, h = net.h;
    if (work.a1.rows() < n || work.a1.cols() != h)
    {
      work.a1.resize(n, h);
      work.a2.resize(n, h);
      work.o.resize(n);
    }
    Eigen::Map<const MatF> X(x, n, d);
    Eigen::Map<const MatF> W1(net.p.data(), d, h), W2(net.p.data() + net.oW2(), h, h);
    Eigen::Map<const RowVecF> b1(net.p.data() + net.oB1(), h), b2(net.p.data() + net.oB2(), h);
    Eigen::Map<const VecF> w3(net.p.data() + net.oW3(), h);
    auto A1 = work.a1.topRows(n);
    auto A2 = work.a2.topRows(n);
    auto o = work.o.head(n);
    A1.noalias() = X * W1;
    A1.rowwise() += b1;
    A1 = A1.array().tanh();
    A2.noalias() = A1 * W2;
    A2.rowwise() += b2;
    A2 = A2.array().tanh();
    o.noalias() = A2 * w3;
    o.array() += net.p[static_cast<Size>(net.oB3())];
  }

  // A trained model: linear (w, bias last) or an MLP ensemble.
  struct Model
  {
    std::vector<double> w; // linear: cols + 1
    std::vector<MlpNet> nets;                // MLP ensemble
    std::vector<std::vector<MlpNet>> folds;  // cross-fitted MLP: one ensemble per fold (nets empty)
    std::vector<double> fold_center, fold_scale; // cross-fitted MLP: score = (output - center) / scale per fold

    bool finite() const
    {
      auto ok = [](const auto& v) { return std::all_of(v.begin(), v.end(), [](auto x) { return std::isfinite(x); }); };
      auto ok_nets = [&](const std::vector<MlpNet>& v) { return std::all_of(v.begin(), v.end(), [&](const MlpNet& n) { return ok(n.p); }); };
      if (!folds.empty()) return std::all_of(folds.begin(), folds.end(), ok_nets);
      if (nets.empty()) return !w.empty() && ok(w);
      return ok_nets(nets);
    }

    // linear model: score of one standardized row
    double score(const double* x) const
    {
      double s = w.back();
      for (Size j = 0; j + 1 < w.size(); ++j) s += w[j] * x[j];
      return s;
    }

    // mean output of @p members (summed in order) for @p n standardized float rows at @p x
    static void scoreEnsemble(const std::vector<MlpNet>& members, const float* x, Eigen::Index n, double* out, MlpWork& work)
    {
      std::fill(out, out + n, 0.0);
      for (const MlpNet& net : members)
      {
        mlpForward(net, x, n, work);
        for (Eigen::Index i = 0; i < n; ++i) out[i] += static_cast<double>(work.o(i));
      }
      for (Eigen::Index i = 0; i < n; ++i) out[i] /= static_cast<double>(members.size());
    }
  };

  // Linear L2-SVM (squared hinge, regularized bias) with per-class costs, via Percolator's TRON solver:
  // minimizes 0.5 ||w||^2 + sum_i C_{y_i} max(0, 1 - y_i w.x_i)^2, the bias (a constant 1.0 column that the solver
  // appends itself) regularized as in liblinear / LinearSVC with intercept_scaling 1.
  Model trainSvm(const Matrix& X, const std::vector<Size>& rows, const std::vector<double>& y, double c_pos, double c_neg)
  {
    const int m = static_cast<int>(rows.size());
    const int n = static_cast<int>(X.cols) + 1;
    P::AlgIn in(static_cast<unsigned int>(m), n);
    in.m = m;
    for (int i = 0; i < m; ++i)
    {
      in.vals[i] = const_cast<double*>(X.row(rows[i]));
      in.Y[i] = y[i];
      if (y[i] > 0) ++in.positives; else ++in.negatives;
    }
    P::vector_double w;
    w.d = n;
    w.vec = new double[n]();
    P::vector_double out;
    out.d = m;
    out.vec = new double[std::max(m, 1)]();
    P::options opt;
    opt.lambda = 1.0;
    opt.lambda_u = 1.0;
    opt.epsilon = P::EPSILON;
    opt.cgitermax = P::CGITERMAX;
    opt.mfnitermax = P::MFNITERMAX;
    // TRON's per-example cost is c / (2 lambda)
    P::tron(in, opt, w, out, 2.0 * c_pos, 2.0 * c_neg);
    Model model;
    model.w.assign(w.vec, w.vec + n);
    return model;
  }

  // Trains one MLP network for @p epochs epochs of mini-batch Adam (batch 256, learning rate 1e-3) on the rows @p rows of
  // @p X with labels @p y: weighted squared hinge (per-class costs C+ / C-, normalized to a total weight of rows.size())
  // plus weight decay 1 / (mean C * rows), as the SVM's regularization. An empty @p net is first initialized (normal,
  // scaled by 1 / sqrt(fan-in)); otherwise training continues from its state (warm start). @p rng draws the
  // initialization and the batch order.
  void trainMlpNet(MlpNet& net, const MatF& X, const std::vector<Size>& rows, const std::vector<double>& y, double c_pos,
                   double c_neg, int hidden, int epochs, std::mt19937_64& rng)
  {
    const Eigen::Index d = X.cols(), h = std::max(hidden, 1);
    const Size m = rows.size();
    if (net.p.empty())
    {
      net.d = d;
      net.h = h;
      net.p.assign(static_cast<Size>(net.oB3() + 1), 0.0f);
      std::normal_distribution<double> nd(0.0, 1.0);
      const double s1 = 1.0 / std::sqrt(static_cast<double>(std::max<Eigen::Index>(d, 1))), s2 = 1.0 / std::sqrt(static_cast<double>(h));
      for (Eigen::Index k = 0; k < d * h; ++k) net.p[static_cast<Size>(k)] = static_cast<float>(nd(rng) * s1);
      for (Eigen::Index k = 0; k < h * h; ++k) net.p[static_cast<Size>(net.oW2() + k)] = static_cast<float>(nd(rng) * s2);
      for (Eigen::Index k = 0; k < h; ++k) net.p[static_cast<Size>(net.oW3() + k)] = static_cast<float>(nd(rng) * s2);
      net.m1.assign(net.p.size(), 0.0f);
      net.m2.assign(net.p.size(), 0.0f);
      net.step = 0;
    }
    if (m == 0) return;
    const double c_mean = (c_pos + c_neg) / 2.0;
    const float decay = static_cast<float>(1.0 / (c_mean * static_cast<double>(m)));
    // class weights normalized so that the total weight is m
    double wsum = 0.0;
    for (Size i = 0; i < m; ++i) wsum += (y[i] > 0 ? c_pos : c_neg);
    const double wnorm = wsum > 0 ? static_cast<double>(m) / wsum : 1.0;
    const float cw_pos = static_cast<float>(c_pos * wnorm), cw_neg = static_cast<float>(c_neg * wnorm);
    constexpr Eigen::Index batch = 256;
    MatF xb(batch, d), a1(batch, h), a2(batch, h), d1(batch, h), d2(batch, h);
    VecF o(batch), g(batch);
    std::vector<float> grad(net.p.size());
    Eigen::Map<const MatF> W1(net.p.data(), d, h), W2(net.p.data() + net.oW2(), h, h);
    Eigen::Map<const RowVecF> b1(net.p.data() + net.oB1(), h), b2(net.p.data() + net.oB2(), h);
    Eigen::Map<const VecF> w3(net.p.data() + net.oW3(), h);
    Eigen::Map<MatF> gW1(grad.data(), d, h), gW2(grad.data() + net.oW2(), h, h);
    Eigen::Map<RowVecF> gb1(grad.data() + net.oB1(), h), gb2(grad.data() + net.oB2(), h);
    Eigen::Map<VecF> gw3(grad.data() + net.oW3(), h);
    const float lr = 1e-3f, beta1 = 0.9f, beta2 = 0.999f, eps = 1e-8f;
    std::vector<Size> order(m);
    std::iota(order.begin(), order.end(), 0);
    for (int epoch = 0; epoch < epochs; ++epoch)
    {
      std::shuffle(order.begin(), order.end(), rng);
      for (Size start = 0; start < m; start += static_cast<Size>(batch))
      {
        const Eigen::Index nb = static_cast<Eigen::Index>(std::min(m - start, static_cast<Size>(batch)));
        for (Eigen::Index b = 0; b < nb; ++b) xb.row(b) = X.row(static_cast<Eigen::Index>(rows[order[start + static_cast<Size>(b)]]));
        auto Xb = xb.topRows(nb);
        auto A1 = a1.topRows(nb);
        auto A2 = a2.topRows(nb);
        auto D1 = d1.topRows(nb);
        auto D2 = d2.topRows(nb);
        auto O = o.head(nb);
        auto G = g.head(nb);
        // forward
        A1.noalias() = Xb * W1;
        A1.rowwise() += b1;
        A1 = A1.array().tanh();
        A2.noalias() = A1 * W2;
        A2.rowwise() += b2;
        A2 = A2.array().tanh();
        O.noalias() = A2 * w3;
        O.array() += net.p[static_cast<Size>(net.oB3())];
        // d loss / d output (0 for rows beyond the margin)
        for (Eigen::Index b = 0; b < nb; ++b)
        {
          const double yi = y[order[start + static_cast<Size>(b)]];
          const float margin = 1.0f - static_cast<float>(yi) * O(b);
          G(b) = margin > 0.0f ? -2.0f * (yi > 0 ? cw_pos : cw_neg) * margin * static_cast<float>(yi) : 0.0f;
        }
        // backward
        gw3.noalias() = A2.transpose() * G;
        grad[static_cast<Size>(net.oB3())] = G.sum();
        D2 = (G * w3.transpose()).array() * (1.0f - A2.array().square());
        gW2.noalias() = A1.transpose() * D2;
        gb2 = D2.colwise().sum();
        D1.noalias() = D2 * W2.transpose();
        D1.array() *= (1.0f - A1.array().square());
        gW1.noalias() = Xb.transpose() * D1;
        gb1 = D1.colwise().sum();
        // Adam
        ++net.step;
        const float inv_nb = 1.0f / static_cast<float>(nb);
        const float c1 = static_cast<float>(1.0 - std::pow(static_cast<double>(beta1), static_cast<double>(net.step)));
        const float c2 = static_cast<float>(1.0 - std::pow(static_cast<double>(beta2), static_cast<double>(net.step)));
        for (Size k = 0; k < net.p.size(); ++k)
        {
          const float gk = grad[k] * inv_nb + decay * net.p[k];
          net.m1[k] = beta1 * net.m1[k] + (1.0f - beta1) * gk;
          net.m2[k] = beta2 * net.m2[k] + (1.0f - beta2) * gk * gk;
          net.p[k] -= lr * (net.m1[k] / c1) / (std::sqrt(net.m2[k] / c2) + eps);
        }
      }
    }
  }

  // AUC = P(a > b) + P(a = b) / 2, by ranks
  double aucOf(const std::vector<double>& a, const std::vector<double>& b)
  {
    if (a.empty() || b.empty()) return 0.5;
    std::vector<std::pair<double, int>> all;
    all.reserve(a.size() + b.size());
    for (double v : a) all.emplace_back(v, 0);
    for (double v : b) all.emplace_back(v, 1);
    std::sort(all.begin(), all.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    double rank_sum_a = 0.0;
    for (Size i = 0; i < all.size();)
    {
      Size j = i;
      while (j < all.size() && all[j].first == all[i].first) ++j;
      const double avg = (static_cast<double>(i) + static_cast<double>(j - 1)) / 2.0 + 1.0;
      for (Size k = i; k < j; ++k)
      {
        if (all[k].second == 0) rank_sum_a += avg;
      }
      i = j;
    }
    const double na = static_cast<double>(a.size()), nb = static_cast<double>(b.size());
    return (rank_sum_a - na * (na + 1.0) / 2.0) / (na * nb);
  }

  double ksOf(std::vector<double> a, std::vector<double> b)
  {
    if (a.empty() || b.empty()) return 0.0;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    Size i = 0, j = 0;
    double best = 0.0;
    while (i < a.size() || j < b.size())
    {
      double v;
      if (j >= b.size() || (i < a.size() && a[i] <= b[j])) v = a[i]; else v = b[j];
      while (i < a.size() && a[i] == v) ++i;
      while (j < b.size() && b[j] == v) ++j;
      const double fa = static_cast<double>(i) / static_cast<double>(a.size());
      const double fb = static_cast<double>(j) / static_cast<double>(b.size());
      best = std::max(best, std::fabs(fa - fb));
    }
    return best;
  }

  // Ridge logistic regression by Newton steps; returns w (cols + 1, bias last)
  std::vector<double> fitLogistic(const Matrix& X, const std::vector<Size>& rows, const std::vector<double>& y01)
  {
    const Size d = X.cols + 1;
    Eigen::VectorXd w = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(d));
    const double ridge = 1e-3 * static_cast<double>(std::max<Size>(rows.size(), 1));
    for (int it = 0; it < 25; ++it)
    {
      Eigen::MatrixXd H = Eigen::MatrixXd::Identity(static_cast<Eigen::Index>(d), static_cast<Eigen::Index>(d)) * ridge;
      H(static_cast<Eigen::Index>(d - 1), static_cast<Eigen::Index>(d - 1)) = 1e-9;
      Eigen::VectorXd g = ridge * w;
      g(static_cast<Eigen::Index>(d - 1)) = 0.0;
      Eigen::VectorXd x(static_cast<Eigen::Index>(d));
      for (Size r = 0; r < rows.size(); ++r)
      {
        const double* xr = X.row(rows[r]);
        for (Size j = 0; j + 1 < d; ++j) x(static_cast<Eigen::Index>(j)) = xr[j];
        x(static_cast<Eigen::Index>(d - 1)) = 1.0;
        const double p = 1.0 / (1.0 + std::exp(-x.dot(w)));
        g += (p - y01[r]) * x;
        H.selfadjointView<Eigen::Lower>().rankUpdate(x, std::max(p * (1.0 - p), 1e-12));
      }
      H.triangularView<Eigen::StrictlyUpper>() = H.transpose();
      const Eigen::VectorXd step = H.ldlt().solve(g);
      if (!step.allFinite()) break;
      w -= step;
      if (step.norm() < 1e-8 * (1.0 + w.norm())) break;
    }
    return std::vector<double>(w.data(), w.data() + w.size());
  }
} // namespace

bool PercolatorReset::isTrainingDecoy(std::uint64_t salt, const std::string& stem)
{
  return (keyedHash(salt, COIN, fnv1a64(stem)) & 1ULL) == 1ULL;
}

std::vector<double> PercolatorReset::seqStepQValues(const std::vector<bool>& is_decoy, double factor, double plus_one)
{
  std::vector<double> q(is_decoy.size(), 1.0);
  double decoys = 0.0, targets = 0.0;
  for (Size k = 0; k < is_decoy.size(); ++k)
  {
    if (is_decoy[k]) decoys += 1.0; else targets += 1.0;
    q[k] = std::min(1.0, factor * (decoys + plus_one) / std::max(targets, 1.0));
  }
  for (Size k = q.size(); k-- > 1;)
  {
    q[k - 1] = std::min(q[k - 1], q[k]);
  }
  return q;
}

PercolatorReset::Output PercolatorReset::run(const Input& in, const Options& opt)
{
  const Size n_rows = in.features.size();
  const Size n_feat = in.feature_names.size();
  if (n_rows == 0 || in.is_decoy.size() != n_rows || in.search_scores.size() != n_rows || in.spectrum.size() != n_rows
      || in.stems.size() != n_rows)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RESET: input vectors are empty or of different lengths.", StringUtils::toStr(n_rows));
  }
  for (Size i = 0; i < n_rows; ++i)
  {
    if (in.features[i].size() != n_feat)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: a feature row has the wrong length.", StringUtils::toStr(in.features[i].size()));
    }
    if (!std::isfinite(in.search_scores[i]))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: non-finite search engine score.", StringUtils::toStr(i));
    }
    for (double v : in.features[i])
    {
      if (!std::isfinite(v))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "RESET: non-finite feature value in row " + StringUtils::toStr(i) + ".", StringUtils::toStr(v));
      }
    }
  }
  if (opt.representative != "search" && opt.representative != "model")
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: unknown representative.", opt.representative);
  }
  if (opt.learner != "svm" && opt.learner != "mlp")
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: unknown learner.", opt.learner);
  }
  auto bad_option = [](const std::string& what, const std::string& value)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: invalid option " + what + ".", value);
  };
  if (!(opt.train_fdr > 0.0 && opt.train_fdr < 1.0)) bad_option("train_fdr (0 < train_fdr < 1)", StringUtils::toStr(opt.train_fdr));
  if (opt.iterations < 1) bad_option("iterations (>= 1)", StringUtils::toStr(opt.iterations));
  if (opt.c_values.empty()) bad_option("c_values (empty)", "");
  for (double c : opt.c_values)
  {
    if (!(std::isfinite(c) && c > 0.0)) bad_option("c_values (finite, > 0)", StringUtils::toStr(c));
  }
  if (opt.learner == "mlp" && (opt.mlp_members < 1 || opt.mlp_hidden < 1 || opt.mlp_epochs < 1 || opt.mlp_warm_epochs < 0
                                || opt.mlp_norm_min_decoys < 0))
  {
    bad_option("mlp_members / mlp_hidden / mlp_epochs (>= 1), mlp_warm_epochs / mlp_norm_min_decoys (>= 0)", StringUtils::toStr(opt.mlp_members));
  }
  if (n_rows >= static_cast<Size>(std::numeric_limits<int>::max()) || n_feat + 1 >= static_cast<Size>(std::numeric_limits<int>::max()))
  {
    bad_option("input size (rows and features must fit an int for the solver)", StringUtils::toStr(n_rows));
  }
  std::uint64_t salt;
  if (opt.salt.has_value())
  {
    salt = *opt.salt;
  }
  else
  {
    std::random_device rd;
    salt = (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());
  }
  Output out;
  out.salt = salt;
  out.learner = opt.learner;
  if (opt.learner == "mlp")
  {
    out.learner_settings = StringUtils::toStr(opt.mlp_members) + " members, hidden " + StringUtils::toStr(opt.mlp_hidden) + ", epochs "
                           + StringUtils::toStr(opt.mlp_epochs) + ", warm epochs " + StringUtils::toStr(opt.mlp_warm_epochs)
                           + (opt.mlp_crossfit ? ", cross-fitted, fold normalization from " + StringUtils::toStr(opt.mlp_norm_min_decoys)
                                                 + " training decoys per fold" : "");
  }

  // per row: stem hash and tie key (random order among equal scores, fixed by the salt)
  std::vector<std::uint64_t> stem_hash(n_rows), row_tie(n_rows);
  for (Size i = 0; i < n_rows; ++i)
  {
    // label-blind (masking): nothing that orders or splits the stems for training may depend on the label
    stem_hash[i] = fnv1a64(in.stems[i]);
    row_tie[i] = keyedHash(salt, TIE, keyedHash(stem_hash[i], TIE, static_cast<std::uint64_t>(i)));
  }
  auto better = [&](Size a, Size b, const std::vector<double>& s)
  {
    if (s[a] != s[b]) return s[a] > s[b];
    return row_tie[a] < row_tie[b];
  };

  // 1. PSM-level competition: the best hit per spectrum by search engine score
  std::unordered_map<Size, Size> best_of_spectrum;
  best_of_spectrum.reserve(n_rows);
  for (Size i = 0; i < n_rows; ++i)
  {
    auto [it, inserted] = best_of_spectrum.try_emplace(in.spectrum[i], i);
    if (!inserted && better(i, it->second, in.search_scores)) it->second = i;
  }
  std::vector<Size> spectrum_winners;
  spectrum_winners.reserve(best_of_spectrum.size());
  for (const auto& [spec, row] : best_of_spectrum) spectrum_winners.push_back(row);
  std::sort(spectrum_winners.begin(), spectrum_winners.end());

  // 2. Stem reduction: best winning PSM per stem (sequence) by search engine score. Label-blind (masking): the stem's
  //    label is the label of that PSM; a sequence with target and decoy PSMs is one stem.
  std::unordered_map<std::string, Size> stem_index; // stem -> index into out.stems
  for (Size row : spectrum_winners)
  {
    auto [it, inserted] = stem_index.try_emplace(in.stems[row], out.stems.size());
    if (inserted)
    {
      Stem s;
      s.stem = in.stems[row];
      s.representative = row;
      out.stems.push_back(std::move(s));
    }
    else if (better(row, out.stems[it->second].representative, in.search_scores))
    {
      out.stems[it->second].representative = row;
    }
  }
  for (Stem& st : out.stems)
  {
    st.decoy = in.is_decoy[st.representative];
    st.search_score = in.search_scores[st.representative];
  }
  // stems in a defined order independent of hash map iteration: by sequence
  {
    std::sort(out.stems.begin(), out.stems.end(), [](const Stem& a, const Stem& b) { return a.stem < b.stem; });
    stem_index.clear();
    for (Size s = 0; s < out.stems.size(); ++s) stem_index[out.stems[s].stem] = s;
  }
  // the PSMs (spectrum winners) of each stem (Options::representative "model"; label-blind: all PSMs of the sequence)
  std::vector<std::vector<Size>> stem_rows(opt.representative == "model" ? out.stems.size() : 0);
  if (opt.representative == "model")
  {
    for (Size row : spectrum_winners) stem_rows[stem_index.at(in.stems[row])].push_back(row);
  }
  // roles: the coin of each decoy stem (all decoy stems, also those that are not considered)
  for (Stem& s : out.stems)
  {
    s.role = !s.decoy ? Role::TARGET : (isTrainingDecoy(salt, s.stem) ? Role::TRAINING_DECOY : Role::ESTIMATING_DECOY);
    s.considered = true;
  }
  // considered stems; counts
  std::vector<Size> cs; // considered stem indices
  for (Size s = 0; s < out.stems.size(); ++s)
  {
    const Stem& st = out.stems[s];
    if (!st.considered) continue;
    cs.push_back(s);
    if (!st.decoy) ++out.target_stems;
    else
    {
      ++out.decoy_stems;
      if (st.role == Role::TRAINING_DECOY) ++out.training_decoys; else ++out.estimating_decoys;
    }
  }
  if (out.training_decoys == 0 || out.estimating_decoys + out.target_stems == 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RESET needs decoy stems for training and pseudo-target stems.", StringUtils::toStr(out.decoy_stems));
  }
  const Size m = cs.size();

  // 3. features of the considered stems: zeroed features, standardization over all stems (no labels), constant
  //    columns dropped
  std::unordered_set<std::string> zero(opt.zero_features.begin(), opt.zero_features.end());
  std::vector<Size> used;
  std::vector<double> mean(n_feat, 0.0), sd(n_feat, 0.0);
  for (Size j = 0; j < n_feat; ++j)
  {
    if (zero.contains(in.feature_names[j])) continue;
    double s1 = 0.0;
    for (Size k = 0; k < m; ++k) s1 += in.features[out.stems[cs[k]].representative][j];
    const double mu = s1 / static_cast<double>(m);
    double s2 = 0.0;
    for (Size k = 0; k < m; ++k)
    {
      const double dv = in.features[out.stems[cs[k]].representative][j] - mu;
      s2 += dv * dv;
    }
    const double sigma = std::sqrt(s2 / static_cast<double>(m));
    if (!std::isfinite(mu) || !std::isfinite(sigma))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: the values of feature '" + in.feature_names[j] + "' are too large to standardize.", StringUtils::toStr(mu));
    }
    if (!(sigma > 1e-12 * std::max(1.0, std::fabs(mu)))) continue;
    mean[j] = mu;
    sd[j] = sigma;
    used.push_back(j);
  }
  if (used.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RESET: no feature varies over the stems.", StringUtils::toStr(m));
  }
  for (Size j : used) out.used_features.push_back(in.feature_names[j]);
  auto standardize = [&](const std::vector<double>& raw, double* dst)
  {
    for (Size u = 0; u < used.size(); ++u) dst[u] = (raw[used[u]] - mean[used[u]]) / sd[used[u]];
  };
  Matrix X;
  X.rows = m;
  X.cols = used.size();
  X.v.resize(m * used.size());
  for (Size k = 0; k < m; ++k) standardize(in.features[out.stems[cs[k]].representative], X.v.data() + k * X.cols);
  // learner "mlp": float32 copies of the stem rows (updated with X) and of all spectrum winners (representative "model")
  const bool mlp = opt.learner == "mlp";
  const bool crossfit = mlp && opt.mlp_crossfit;
  const Eigen::Index mlp_chunk = 1024; // rows per scoring block
  const int mlp_threads = std::max(1, opt.threads);
  MatF Xf;
  MatF Xw;                             // standardized spectrum winners, in the order of spectrum_winners
  std::vector<Size> winner_pos;        // row -> position in spectrum_winners (representative "model")
  // fixed fold of a sequence (warm starts, cross-fitting): a hash of salt and sequence (label-blind)
  auto foldOf = [&](std::uint64_t seq_hash) { return static_cast<int>(keyedHash(salt, FOLD, seq_hash) % 3); };
  std::vector<int> stem_fold(m, 0);
  // rows by ensemble (cross-fitting: by fold; else one group) for the stems, the spectrum winners and all rows
  std::vector<std::vector<Size>> stem_groups, winner_groups, row_groups;
  auto standardizeF = [&](const std::vector<double>& raw, float* dst)
  {
    for (Size u = 0; u < used.size(); ++u) dst[u] = static_cast<float>((raw[used[u]] - mean[used[u]]) / sd[used[u]]);
  };
  if (mlp)
  {
    Xf = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
           X.v.data(), static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(X.cols)).cast<float>();
    const Size n_groups = crossfit ? 3 : 1;
    stem_groups.resize(n_groups);
    row_groups.resize(n_groups);
    for (Size k = 0; k < m; ++k)
    {
      stem_fold[k] = foldOf(stem_hash[out.stems[cs[k]].representative]);
      stem_groups[crossfit ? static_cast<Size>(stem_fold[k]) : 0].push_back(k);
    }
    for (Size i = 0; i < n_rows; ++i) row_groups[crossfit ? static_cast<Size>(foldOf(stem_hash[i])) : 0].push_back(i);
    if (opt.representative == "model")
    {
      winner_groups.resize(n_groups);
      winner_pos.assign(n_rows, 0);
      Xw.resize(static_cast<Eigen::Index>(spectrum_winners.size()), static_cast<Eigen::Index>(X.cols));
      for (Size w = 0; w < spectrum_winners.size(); ++w)
      {
        winner_pos[spectrum_winners[w]] = w;
        standardizeF(in.features[spectrum_winners[w]], Xw.data() + static_cast<Eigen::Index>(w) * Xw.cols());
        winner_groups[crossfit ? static_cast<Size>(foldOf(stem_hash[spectrum_winners[w]])) : 0].push_back(w);
      }
    }
  }
  // MLP scores: dst[i] = ensemble output for item i of @p groups (group g: fold g's ensemble if the model is
  // cross-fitted, else the model's ensemble); fill(i, dst) writes item i's standardized float row. Parallel over
  // blocks of mlp_chunk items of one group.
  auto scoreGrouped = [&](const Model& md, const std::vector<std::vector<Size>>& groups, Size n, const auto& fill,
                          std::vector<double>& dst)
  {
    dst.assign(n, 0.0);
    std::vector<std::pair<Size, Size>> blocks; // (group, first position)
    for (Size g = 0; g < groups.size(); ++g)
    {
      for (Size s0 = 0; s0 < groups[g].size(); s0 += static_cast<Size>(mlp_chunk)) blocks.emplace_back(g, s0);
    }
#pragma omp parallel num_threads(mlp_threads)
    {
      MlpWork work;
      MatF xc(mlp_chunk, static_cast<Eigen::Index>(used.size()));
      std::vector<double> sc(static_cast<Size>(mlp_chunk));
#pragma omp for schedule(static)
      for (SignedSize b = 0; b < static_cast<SignedSize>(blocks.size()); ++b)
      {
        const auto [g, s0] = blocks[static_cast<Size>(b)];
        const Size nr = std::min(static_cast<Size>(mlp_chunk), groups[g].size() - s0);
        for (Size r = 0; r < nr; ++r) fill(groups[g][s0 + r], xc.data() + static_cast<Eigen::Index>(r) * xc.cols());
        Model::scoreEnsemble(md.folds.empty() ? md.nets : md.folds[g], xc.data(), static_cast<Eigen::Index>(nr), sc.data(), work);
        const bool norm = !md.folds.empty() && !md.fold_scale.empty();
        for (Size r = 0; r < nr; ++r) dst[groups[g][s0 + r]] = norm ? (sc[r] - md.fold_center[g]) / md.fold_scale[g] : sc[r];
      }
    }
  };
  auto fillStem = [&](Size k, float* dst) { std::copy_n(Xf.data() + static_cast<Eigen::Index>(k) * Xf.cols(), Xf.cols(), dst); };
  auto fillWinner = [&](Size w, float* dst) { std::copy_n(Xw.data() + static_cast<Eigen::Index>(w) * Xw.cols(), Xw.cols(), dst); };

  // pseudo-labels (masking): training decoys -1, pseudo-targets (targets, estimating decoys) +1
  std::vector<bool> training_decoy(m);
  std::vector<std::uint64_t> tie(m);
  for (Size k = 0; k < m; ++k)
  {
    training_decoy[k] = out.stems[cs[k]].role == Role::TRAINING_DECOY;
    tie[k] = row_tie[out.stems[cs[k]].representative];
  }

  // 4. feature-symmetry audit (training decoys vs low-scoring pseudo-targets, same score range)
  if (opt.audit)
  {
    // null proxy: the stems below the median search engine score of the training decoys; under the null, training
    // decoys and null pseudo-targets share one distribution, so few correct targets fall below that score
    std::vector<double> td_scores;
    for (Size k = 0; k < m; ++k)
    {
      if (training_decoy[k]) td_scores.push_back(out.stems[cs[k]].search_score);
    }
    std::vector<Size> grp_a, grp_b; // training decoys, null proxy (pseudo-targets)
    {
      const Size half = td_scores.size() / 2;
      std::nth_element(td_scores.begin(), td_scores.begin() + static_cast<std::ptrdiff_t>(half), td_scores.end());
      const double cutoff = td_scores[half];
      for (Size k = 0; k < m; ++k)
      {
        if (out.stems[cs[k]].search_score > cutoff) continue;
        (training_decoy[k] ? grp_a : grp_b).push_back(k);
      }
    }
    // about 3 null pseudo-targets per training decoy (2 incorrect targets, 1 estimating decoy): the excess over 3 is
    // the share of correct targets in the proxy
    out.audit_decoys = grp_a.size();
    out.audit_proxy = grp_b.size();
    out.audit_proxy_ratio = grp_a.empty() ? 0.0 : static_cast<double>(grp_b.size()) / static_cast<double>(grp_a.size());
    if (grp_a.size() >= 20 && grp_b.size() >= 20)
    {
      for (Size j = 0; j < n_feat; ++j)
      {
        std::vector<double> va, vb;
        for (Size k : grp_a) va.push_back(in.features[out.stems[cs[k]].representative][j]);
        for (Size k : grp_b) vb.push_back(in.features[out.stems[cs[k]].representative][j]);
        AuditFeature af;
        af.name = in.feature_names[j];
        af.auc = aucOf(va, vb);
        af.ks = ksOf(va, vb);
        af.flagged = std::fabs(af.auc - 0.5) > 0.05;
        out.audit.push_back(af);
      }
      // cross-fitted logistic regression on the standardized used features (2 folds by stem hash)
      std::vector<Size> fold_rows[2];
      std::vector<double> fold_y[2];
      for (Size k : grp_a)
      {
        const int f = static_cast<int>(keyedHash(salt, AUDIT, stem_hash[out.stems[cs[k]].representative]) & 1ULL);
        fold_rows[f].push_back(k);
        fold_y[f].push_back(1.0);
      }
      for (Size k : grp_b)
      {
        const int f = static_cast<int>(keyedHash(salt, AUDIT, stem_hash[out.stems[cs[k]].representative]) & 1ULL);
        fold_rows[f].push_back(k);
        fold_y[f].push_back(0.0);
      }
      std::vector<double> pa, pb;
      for (int f = 0; f < 2; ++f)
      {
        if (fold_rows[1 - f].empty()) continue;
        const std::vector<double> w = fitLogistic(X, fold_rows[1 - f], fold_y[1 - f]);
        for (Size r = 0; r < fold_rows[f].size(); ++r)
        {
          const double* x = X.row(fold_rows[f][r]);
          double s = w.back();
          for (Size j = 0; j < X.cols; ++j) s += w[j] * x[j];
          (fold_y[f][r] > 0.5 ? pa : pb).push_back(s);
        }
      }
      out.audit_multivariate_auc = aucOf(pa, pb);
      out.audit_flagged = out.audit_multivariate_auc > 0.55;
    }
  }

  // 5. training rounds
  std::vector<double> cur(m);
  for (Size k = 0; k < m; ++k) cur[k] = out.stems[cs[k]].search_score; // initial direction: search engine score
  std::vector<std::pair<double, double>> grid;
  for (double cp : opt.c_values)
  {
    for (double cn : opt.c_values) grid.emplace_back(cp, cn);
  }
  if (grid.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: empty C grid.", "");
  }
  const int threads = std::max(1, opt.threads);
  Model model;
  // positives: pseudo-targets passing @p alpha against the training decoys (SeqStep(+), c = 3/4)
  auto positiveSet = [&](const std::vector<Size>& ranked, double alpha, double plus_one)
  {
    std::vector<bool> dec(ranked.size());
    for (Size r = 0; r < ranked.size(); ++r) dec[r] = training_decoy[ranked[r]];
    const std::vector<double> q = seqStepQValues(dec, 3.0, plus_one);
    std::vector<Size> pos;
    for (Size r = 0; r < ranked.size(); ++r)
    {
      if (!dec[r] && q[r] <= alpha) pos.push_back(ranked[r]);
    }
    return pos;
  };
  std::vector<Size> ranked(m);
  // learner "mlp": networks per cross-validation task (fold, C pair) and of the final model; with warm starts
  // (Options::mlp_warm_epochs > 0) they carry over to the next round
  const bool warm = mlp && opt.mlp_warm_epochs > 0;
  std::vector<std::vector<MlpNet>> cv_nets(3 * grid.size());
  std::vector<MlpNet> final_nets;
  for (int round = 0; round < std::max(1, opt.iterations); ++round)
  {
    std::iota(ranked.begin(), ranked.end(), 0);
    rankByScore(ranked, cur, tie);
    double alpha = opt.train_fdr;
    std::vector<Size> pos = positiveSet(ranked, alpha, 1.0);
    if (pos.empty()) pos = positiveSet(ranked, alpha, 0.0);
    while (pos.empty() && alpha < 1.0)
    {
      alpha += 0.005;
      pos = positiveSet(ranked, alpha, 0.0);
    }
    if (round > 0 && pos.size() < 3)
    {
      // a later round whose model leaves too few positives (weak data): keep the previous round's model and scores
      out.rounds.push_back("round " + StringUtils::toStr(round + 1) + ": " + StringUtils::toStr(pos.size())
                           + " positive stems: training stopped, the model of round " + StringUtils::toStr(round) + " is kept");
      break;
    }
    if (pos.empty())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: no positive training stems.", StringUtils::toStr(round));
    }
    std::vector<Size> neg;
    for (Size k = 0; k < m; ++k)
    {
      if (training_decoy[k]) neg.push_back(k);
    }
    // stratified folds: positives and negatives each in a random (salted) order, dealt round robin
    auto foldOrder = [&](std::vector<Size> v)
    {
      std::vector<std::uint64_t> key(m);
      for (Size k : v) key[k] = keyedHash(salt, FOLD + static_cast<std::uint64_t>(round), tie[k]);
      std::sort(v.begin(), v.end(), [&](Size a, Size b) { return key[a] != key[b] ? key[a] < key[b] : a < b; });
      return v;
    };
    std::vector<Size> train_rows;
    std::vector<double> train_y;
    std::vector<int> train_fold;
    {
      const std::vector<Size> po = foldOrder(pos), no = foldOrder(neg);
      for (Size r = 0; r < po.size(); ++r) { train_rows.push_back(po[r]); train_y.push_back(1.0); train_fold.push_back(static_cast<int>(r % 3)); }
      for (Size r = 0; r < no.size(); ++r) { train_rows.push_back(no[r]); train_y.push_back(-1.0); train_fold.push_back(static_cast<int>(r % 3)); }
    }
    // warm starts and cross-fitting: a stem keeps its fold over the rounds (a network never trained on its held-out
    // stems in an earlier round; cross-fitting scores a stem by the networks of its fold); round robin as above (and no
    // cross-fitting in this round) if a fold would lack positives or negatives
    bool fixed_folds = false;
    if (warm || crossfit)
    {
      std::vector<int> fixed(train_rows.size());
      Size per_fold[3][2] = {{0, 0}, {0, 0}, {0, 0}};
      for (Size r = 0; r < train_rows.size(); ++r)
      {
        fixed[r] = stem_fold[train_rows[r]];
        ++per_fold[fixed[r]][train_y[r] > 0 ? 1 : 0];
      }
      fixed_folds = true;
      for (const auto& pf : per_fold) fixed_folds = fixed_folds && pf[0] > 0 && pf[1] > 0;
      if (fixed_folds) train_fold = std::move(fixed);
    }
    const bool crossfit_round = crossfit && fixed_folds;
    if (pos.size() < 3 || neg.size() < 3)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: too few positive stems or training decoys for 3-fold cross-validation (at least 3 each).",
        StringUtils::toStr(pos.size()) + " / " + StringUtils::toStr(neg.size()));
    }
    // learner "svm"
    auto fit = [&](const std::vector<Size>& rows, const std::vector<double>& y, double cp, double cn)
    {
      Model fm = trainSvm(X, rows, y, cp, cn);
      if (!fm.finite())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "RESET: the learner returned non-finite parameters.", opt.learner);
      }
      return fm;
    };
    // cross-validation: per (fold, pair) the held-out scores
    const Size tasks = 3 * grid.size();
    std::vector<std::vector<double>> held_scores(tasks);
    std::vector<std::vector<Size>> held_index(3); // positions in train_rows
    for (Size r = 0; r < train_rows.size(); ++r) held_index[train_fold[r]].push_back(r);
    std::vector<std::exception_ptr> errors(tasks);
    const int members = opt.mlp_members;
    const int epochs = warm && round > 0 ? opt.mlp_warm_epochs : opt.mlp_epochs;
    if (mlp)
    {
      // every network of every task is one parallel job; then the held-out scores of each task's ensemble
      std::vector<std::vector<Size>> fold_rows(3);
      std::vector<std::vector<double>> fold_y(3);
      for (int f = 0; f < 3; ++f)
      {
        for (Size r = 0; r < train_rows.size(); ++r)
        {
          if (train_fold[r] == f) continue;
          fold_rows[f].push_back(train_rows[r]);
          fold_y[f].push_back(train_y[r]);
        }
      }
      for (std::vector<MlpNet>& v : cv_nets)
      {
        if (!warm) v.clear();
        v.resize(static_cast<Size>(members));
      }
      const SignedSize jobs = static_cast<SignedSize>(tasks) * members;
      std::vector<std::exception_ptr> job_errors(static_cast<Size>(jobs));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
      for (SignedSize j = 0; j < jobs; ++j)
      {
        try
        {
          const Size t = static_cast<Size>(j / members);
          const int f = static_cast<int>(t / grid.size());
          const auto [cp, cn] = grid[t % grid.size()];
          const std::uint64_t task_seed = keyedHash(salt, MLP, static_cast<std::uint64_t>(round) * 100 + static_cast<std::uint64_t>(f));
          std::mt19937_64 rng(keyedHash(task_seed, MLP, static_cast<std::uint64_t>(j % members)));
          trainMlpNet(cv_nets[t][static_cast<Size>(j % members)], Xf, fold_rows[f], fold_y[f], cp, cn, opt.mlp_hidden, epochs, rng);
        }
        catch (...)
        {
          job_errors[static_cast<Size>(j)] = std::current_exception();
        }
      }
      for (const std::exception_ptr& e : job_errors)
      {
        if (e) std::rethrow_exception(e);
      }
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
      for (SignedSize t = 0; t < static_cast<SignedSize>(tasks); ++t)
      {
        try
        {
          const std::vector<MlpNet>& nets = cv_nets[static_cast<Size>(t)];
          const int f = static_cast<int>(t / static_cast<SignedSize>(grid.size()));
          if (!std::all_of(nets.begin(), nets.end(), [](const MlpNet& n) { return std::all_of(n.p.begin(), n.p.end(), [](float v) { return std::isfinite(v); }); }))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: the learner returned non-finite parameters.", opt.learner);
          }
          MatF xh(static_cast<Eigen::Index>(held_index[f].size()), Xf.cols());
          for (Size r = 0; r < held_index[f].size(); ++r) xh.row(static_cast<Eigen::Index>(r)) = Xf.row(static_cast<Eigen::Index>(train_rows[held_index[f][r]]));
          std::vector<double>& hs = held_scores[static_cast<Size>(t)];
          hs.assign(held_index[f].size(), 0.0);
          MlpWork work;
          for (Eigen::Index r0 = 0; r0 < xh.rows(); r0 += mlp_chunk)
          {
            Model::scoreEnsemble(nets, xh.data() + r0 * xh.cols(), std::min(mlp_chunk, xh.rows() - r0), hs.data() + r0, work);
          }
          for (double v : hs)
          {
            if (!std::isfinite(v))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: non-finite held-out score.", StringUtils::toStr(v));
            }
          }
        }
        catch (...)
        {
          errors[static_cast<Size>(t)] = std::current_exception();
        }
      }
    }
    else
    {
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
    for (SignedSize t = 0; t < static_cast<SignedSize>(tasks); ++t)
    {
      try
      {
        const int f = static_cast<int>(t / static_cast<SignedSize>(grid.size()));
        const auto [cp, cn] = grid[static_cast<Size>(t) % grid.size()];
        std::vector<Size> rows;
        std::vector<double> y;
        for (Size r = 0; r < train_rows.size(); ++r)
        {
          if (train_fold[r] == f) continue;
          rows.push_back(train_rows[r]);
          y.push_back(train_y[r]);
        }
        const Model fm = fit(rows, y, cp, cn);
        std::vector<double>& hs = held_scores[static_cast<Size>(t)];
        for (Size r : held_index[f])
        {
          hs.push_back(fm.score(X.row(train_rows[r])));
          if (!std::isfinite(hs.back()))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: non-finite held-out score.",
                                          StringUtils::toStr(hs.back()));
          }
        }
      }
      catch (...)
      {
        errors[static_cast<Size>(t)] = std::current_exception();
      }
    }
    }
    for (const std::exception_ptr& e : errors)
    {
      if (e) std::rethrow_exception(e);
    }
    // discoveries per pair: SeqStep+ (c = 3/4) on each held-out fold, positives as targets, negatives as decoys
    auto countPair = [&](Size p, double a, double plus_one)
    {
      Size total = 0;
      for (int f = 0; f < 3; ++f)
      {
        const std::vector<double>& hs = held_scores[static_cast<Size>(f) * grid.size() + p];
        std::vector<Size> order(held_index[f].size());
        std::iota(order.begin(), order.end(), 0);
        std::vector<std::uint64_t> otie(order.size());
        for (Size r = 0; r < order.size(); ++r) otie[r] = tie[train_rows[held_index[f][r]]];
        rankByScore(order, hs, otie);
        std::vector<bool> dec(order.size());
        for (Size r = 0; r < order.size(); ++r) dec[r] = train_y[held_index[f][order[r]]] < 0;
        const std::vector<double> q = seqStepQValues(dec, 3.0, plus_one);
        for (Size r = 0; r < order.size(); ++r)
        {
          if (!dec[r] && q[r] <= a) ++total;
        }
      }
      return total;
    };
    double cv_alpha = opt.train_fdr, cv_plus = 1.0;
    std::vector<Size> counts(grid.size());
    Size best = 0;
    for (int attempt = 0;; ++attempt)
    {
      for (Size p = 0; p < grid.size(); ++p) counts[p] = countPair(p, cv_alpha, cv_plus);
      best = static_cast<Size>(std::max_element(counts.begin(), counts.end()) - counts.begin()); // first max
      if (counts[best] > 0 || cv_alpha >= 1.0) break;
      if (cv_plus > 0.0) cv_plus = 0.0; else cv_alpha += 0.005;
    }
    if (crossfit_round)
    {
      // cross-fitted model: the ensembles of the chosen C pair, one per fold (normalized below)
      model.nets.clear();
      model.folds.clear();
      model.fold_center.clear();
      model.fold_scale.clear();
      for (int f = 0; f < 3; ++f) model.folds.push_back(cv_nets[static_cast<Size>(f) * grid.size() + best]);
    }
    else if (mlp)
    {
      model.folds.clear();
      if (!warm) final_nets.clear();
      final_nets.resize(static_cast<Size>(members));
      const std::uint64_t final_seed = keyedHash(salt, MLP, static_cast<std::uint64_t>(round) * 100 + 99);
      std::vector<std::exception_ptr> job_errors(static_cast<Size>(members));
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
      for (int j = 0; j < members; ++j)
      {
        try
        {
          std::mt19937_64 rng(keyedHash(final_seed, MLP, static_cast<std::uint64_t>(j)));
          trainMlpNet(final_nets[static_cast<Size>(j)], Xf, train_rows, train_y, grid[best].first, grid[best].second, opt.mlp_hidden, epochs, rng);
        }
        catch (...)
        {
          job_errors[static_cast<Size>(j)] = std::current_exception();
        }
      }
      for (const std::exception_ptr& e : job_errors)
      {
        if (e) std::rethrow_exception(e);
      }
      model.nets = final_nets;
      if (!model.finite())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: the learner returned non-finite parameters.", opt.learner);
      }
    }
    else
    {
      model = fit(train_rows, train_y, grid[best].first, grid[best].second);
    }
    if (opt.representative == "model")
    {
      // the stem's PSM with the best score under the new model represents it from now on (label-blind: a function of
      // the stem's PSMs and the model)
      std::vector<std::exception_ptr> rep_errors(m);
      std::vector<double> winner_score; // learner "mlp": all spectrum winners scored in blocks first
      if (mlp) scoreGrouped(model, winner_groups, spectrum_winners.size(), fillWinner, winner_score);
#pragma omp parallel num_threads(threads)
      {
        std::vector<double> x(X.cols);
#pragma omp for schedule(dynamic, 256)
        for (SignedSize kk = 0; kk < static_cast<SignedSize>(m); ++kk)
        {
          try
          {
            const Size k = static_cast<Size>(kk);
            Stem& st = out.stems[cs[k]];
            double best_score = -std::numeric_limits<double>::infinity();
            Size best_row = st.representative;
            for (Size row : stem_rows[cs[k]])
            {
              double sc;
              if (mlp)
              {
                sc = winner_score[winner_pos[row]];
              }
              else
              {
                standardize(in.features[row], x.data());
                sc = model.score(x.data());
              }
              if (sc > best_score || (sc == best_score && row_tie[row] < row_tie[best_row]))
              {
                best_score = sc;
                best_row = row;
              }
            }
            if (best_row != st.representative)
            {
              st.representative = best_row;
              st.search_score = in.search_scores[best_row];
              standardize(in.features[best_row], X.v.data() + k * X.cols);
              if (mlp) Xf.row(static_cast<Eigen::Index>(k)) = Xw.row(static_cast<Eigen::Index>(winner_pos[best_row]));
            }
          }
          catch (...)
          {
            rep_errors[static_cast<Size>(kk)] = std::current_exception();
          }
        }
      }
      for (const std::exception_ptr& e : rep_errors)
      {
        if (e) std::rethrow_exception(e);
      }
    }
    if (mlp) scoreGrouped(model, stem_groups, m, fillStem, cur);
    if (crossfit_round)
    {
      // the three fold ensembles are trained separately, so their score scales differ, but the stems are ranked
      // together: as Percolator normalizes the scores of its cross-validation folds, each fold's scores are shifted and
      // scaled so that the score at which its pseudo-targets pass Options::train_fdr against its training decoys
      // (SeqStep+, c = 3/4, as the positive set) becomes 0 and the median score of its training decoys -1 (labels of
      // training decoys only: masking holds); a fold without such a threshold keeps its scores
      model.fold_center.assign(3, 0.0);
      model.fold_scale.assign(3, 1.0);
      // only with enough training decoys in every fold (Options::mlp_norm_min_decoys): with few, the anchors are set by
      // a handful of decoys and add noise; all folds or none, so that raw and normalized folds are never mixed
      Size fold_decoys[3] = {0, 0, 0};
      for (Size k = 0; k < m; ++k) fold_decoys[stem_fold[k]] += training_decoy[k] ? 1 : 0;
      const bool normalize = std::min({fold_decoys[0], fold_decoys[1], fold_decoys[2]}) >= static_cast<Size>(opt.mlp_norm_min_decoys);
      for (int f = 0; normalize && f < 3; ++f)
      {
        std::vector<Size> fk;
        std::vector<double> dec_scores;
        for (Size k = 0; k < m; ++k)
        {
          if (stem_fold[k] != f) continue;
          fk.push_back(k);
          if (training_decoy[k]) dec_scores.push_back(cur[k]);
        }
        if (dec_scores.empty()) continue;
        rankByScore(fk, cur, tie);
        std::vector<bool> dec(fk.size());
        for (Size r = 0; r < fk.size(); ++r) dec[r] = training_decoy[fk[r]];
        const std::vector<double> q = seqStepQValues(dec, 3.0, 1.0);
        double thr = std::numeric_limits<double>::quiet_NaN();
        for (Size r = 0; r < fk.size(); ++r)
        {
          if (!dec[r] && q[r] <= opt.train_fdr) thr = cur[fk[r]];
        }
        std::nth_element(dec_scores.begin(), dec_scores.begin() + static_cast<std::ptrdiff_t>(dec_scores.size() / 2), dec_scores.end());
        const double med = dec_scores[dec_scores.size() / 2];
        if (!std::isfinite(thr) || !(thr - med > 1e-12)) continue;
        model.fold_center[static_cast<Size>(f)] = thr;
        model.fold_scale[static_cast<Size>(f)] = thr - med;
      }
      for (Size k = 0; k < m; ++k) cur[k] = (cur[k] - model.fold_center[static_cast<Size>(stem_fold[k])]) / model.fold_scale[static_cast<Size>(stem_fold[k])];
    }
    for (Size k = 0; k < m; ++k)
    {
      if (!mlp) cur[k] = model.score(X.row(k));
      if (!std::isfinite(cur[k]))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: non-finite stem score.",
                                      StringUtils::toStr(cur[k]));
      }
    }
    std::ostringstream os;
    os << "round " << (round + 1) << ": positives " << pos.size() << " (train FDR " << alpha << "), training decoys "
       << neg.size() << ", C+ " << grid[best].first << ", C- " << grid[best].second << ", held-out discoveries "
       << counts[best];
    out.rounds.push_back(os.str());
  }

  // 6. final ranking: targets and estimating decoys, FDP = 2 (D_est + 1) / max(T, 1)
  {
    std::vector<Size> fin;
    for (Size k = 0; k < m; ++k)
    {
      out.stems[cs[k]].score = cur[k];
      if (!training_decoy[k]) fin.push_back(k);
    }
    rankByScore(fin, cur, tie);
    std::vector<bool> dec(fin.size());
    for (Size r = 0; r < fin.size(); ++r) dec[r] = out.stems[cs[fin[r]]].decoy;
    const std::vector<double> q = seqStepQValues(dec, 2.0, 1.0);
    for (Size r = 0; r < fin.size(); ++r)
    {
      Stem& st = out.stems[cs[fin[r]]];
      st.q_value = q[r];
      if (!st.decoy && q[r] <= 0.01) ++out.discoveries_1pct;
    }
  }
  if (opt.learner == "svm")
  {
    // weights on raw features (bias last)
    out.weights.assign(n_feat + 1, 0.0);
    double bias = model.w.back();
    for (Size u = 0; u < used.size(); ++u)
    {
      out.weights[used[u]] = model.w[u] / sd[used[u]];
      bias -= model.w[u] * mean[used[u]] / sd[used[u]];
    }
    out.weights[n_feat] = bias;
  }

  // 7. per row: final scores, stem q-values and roles; heuristic PSM-level q-values and PEPs
  out.scores.assign(n_rows, 0.0);
  out.psm_q_values.assign(n_rows, 1.0);
  out.psm_peps.assign(n_rows, 1.0);
  out.stem_q_values.assign(n_rows, 1.0);
  out.roles.assign(n_rows, Role::TARGET);
  if (mlp)
  {
    scoreGrouped(model, row_groups, n_rows, [&](Size i, float* dst) { standardizeF(in.features[i], dst); }, out.scores);
    for (Size i = 0; i < n_rows; ++i)
    {
      if (!std::isfinite(out.scores[i]))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: non-finite PSM score.", StringUtils::toStr(out.scores[i]));
      }
    }
  }
  else
  {
    std::vector<std::exception_ptr> errors(n_rows);
#pragma omp parallel num_threads(threads)
    {
      std::vector<double> x(used.size());
#pragma omp for schedule(static)
      for (SignedSize i = 0; i < static_cast<SignedSize>(n_rows); ++i)
      {
        try
        {
          standardize(in.features[static_cast<Size>(i)], x.data());
          out.scores[static_cast<Size>(i)] = model.score(x.data());
          if (!std::isfinite(out.scores[static_cast<Size>(i)]))
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "RESET: non-finite PSM score.",
                                          StringUtils::toStr(out.scores[static_cast<Size>(i)]));
          }
        }
        catch (...)
        {
          errors[static_cast<Size>(i)] = std::current_exception();
        }
      }
    }
    for (const std::exception_ptr& e : errors)
    {
      if (e) std::rethrow_exception(e);
    }
  }
  for (Size i = 0; i < n_rows; ++i)
  {
    // the row's role follows its own label (decoy rows: the coin of their sequence); the stem's q-value is exported
    // only to rows of the stem's label (review r2: a target row of a sequence whose stem is a decoy must not inherit
    // that decoy's q-value)
    out.roles[i] = !in.is_decoy[i] ? Role::TARGET
                   : (isTrainingDecoy(salt, in.stems[i]) ? Role::TRAINING_DECOY : Role::ESTIMATING_DECOY);
    auto it = stem_index.find(in.stems[i]);
    if (it != stem_index.end())
    {
      const Stem& st = out.stems[it->second];
      out.stem_q_values[i] = st.considered && st.decoy == in.is_decoy[i] ? st.q_value : 1.0;
    }
  }
  {
    // the PSM of a spectrum: its hit with the best search engine score (the PSM-level competition of step 1, so that
    // the PSM that represents a stem is reported); the PSMs of training-decoy stems are dropped
    std::vector<Size> psms;
    for (Size row : spectrum_winners)
    {
      if (out.roles[row] != Role::TRAINING_DECOY) psms.push_back(row);
    }
    std::sort(psms.begin(), psms.end());
    rankByScore(psms, out.scores, row_tie);
    std::vector<bool> dec(psms.size());
    for (Size r = 0; r < psms.size(); ++r) dec[r] = in.is_decoy[psms[r]];
    const std::vector<double> q = seqStepQValues(dec, 2.0, 1.0);
    // PEP: antitonic regression (pool adjacent violators) of the estimating-decoy indicator over the ranks
    std::vector<double> blk_val, blk_w;
    std::vector<Size> blk_len;
    // one initial block per distinct score (review r3: tied scores get one PEP, whatever their tie-broken order)
    for (Size r = 0; r < dec.size();)
    {
      Size e = r;
      double decoys = 0.0;
      while (e < dec.size() && out.scores[psms[e]] == out.scores[psms[r]])
      {
        decoys += dec[e] ? 1.0 : 0.0;
        ++e;
      }
      blk_val.push_back(decoys / static_cast<double>(e - r));
      blk_w.push_back(static_cast<double>(e - r));
      blk_len.push_back(e - r);
      r = e;
      while (blk_val.size() > 1 && blk_val[blk_val.size() - 2] > blk_val.back())
      {
        const Size b = blk_val.size() - 1;
        const double wsum = blk_w[b - 1] + blk_w[b];
        blk_val[b - 1] = (blk_val[b - 1] * blk_w[b - 1] + blk_val[b] * blk_w[b]) / wsum;
        blk_w[b - 1] = wsum;
        blk_len[b - 1] += blk_len[b];
        blk_val.pop_back();
        blk_w.pop_back();
        blk_len.pop_back();
      }
    }
    Size r = 0;
    for (Size b = 0; b < blk_val.size(); ++b)
    {
      const double p = blk_val[b];
      const double pep = p >= 1.0 / 3.0 ? 1.0 : 2.0 * p / (1.0 - p);
      for (Size k = 0; k < blk_len[b]; ++k, ++r)
      {
        out.psm_q_values[psms[r]] = q[r];
        out.psm_peps[psms[r]] = pep;
      }
    }
  }
  return out;
}
} // namespace OpenMS
