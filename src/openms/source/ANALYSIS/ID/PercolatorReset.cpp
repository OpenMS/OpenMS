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

  // A trained model: linear (w, bias last) or an MLP ensemble.
  struct Model
  {
    std::vector<double> w; // linear: cols + 1
    struct Net
    {
      Size d = 0, h = 0;
      std::vector<double> w1, b1, w2, b2, w3; // d x h, h, h x h, h, h
      double b3 = 0.0;
    };
    std::vector<Net> nets;

    bool finite() const
    {
      auto ok = [](const std::vector<double>& v) { return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); }); };
      if (nets.empty()) return !w.empty() && ok(w);
      for (const Net& n : nets)
      {
        if (!ok(n.w1) || !ok(n.b1) || !ok(n.w2) || !ok(n.b2) || !ok(n.w3) || !std::isfinite(n.b3)) return false;
      }
      return true;
    }

    double score(const double* x) const
    {
      if (nets.empty())
      {
        double s = w.back();
        for (Size j = 0; j + 1 < w.size(); ++j) s += w[j] * x[j];
        return s;
      }
      double sum = 0.0;
      std::vector<double> a1, a2;
      for (const Net& n : nets)
      {
        a1.assign(n.h, 0.0);
        a2.assign(n.h, 0.0);
        for (Size k = 0; k < n.h; ++k)
        {
          double s = n.b1[k];
          for (Size j = 0; j < n.d; ++j) s += n.w1[j * n.h + k] * x[j];
          a1[k] = std::tanh(s);
        }
        for (Size k = 0; k < n.h; ++k)
        {
          double s = n.b2[k];
          for (Size j = 0; j < n.h; ++j) s += n.w2[j * n.h + k] * a1[j];
          a2[k] = std::tanh(s);
        }
        double o = n.b3;
        for (Size k = 0; k < n.h; ++k) o += n.w3[k] * a2[k];
        sum += o;
      }
      return sum / static_cast<double>(nets.size());
    }
  };

  // Linear L2-SVM (squared hinge, regularized bias) with per-class costs, via Percolator's TRON solver:
  // minimizes 0.5 ||w||^2 + sum_i C_{y_i} max(0, 1 - y_i w.x_i)^2 (the objective of LinearSVC with class weights).
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

  // Small MLP ensemble (two tanh hidden layers, linear output), weighted squared-hinge loss + weight decay
  // (decay = 1 / (C * rows), as the SVM's regularization), mini-batch Adam, fixed epochs; deterministic seeds
  // from the salt. The members differ by their initialization and batch order only.
  Model trainMlp(const Matrix& X, const std::vector<Size>& rows, const std::vector<double>& y, double c_pos, double c_neg,
                 std::uint64_t seed, int members, int hidden, int epochs)
  {
    Model model;
    const Size d = X.cols, h = static_cast<Size>(std::max(hidden, 1)), m = rows.size();
    const double c_mean = (c_pos + c_neg) / 2.0;
    const double decay = 1.0 / (c_mean * std::max<double>(static_cast<double>(m), 1.0));
    for (int member = 0; member < members; ++member)
    {
      std::mt19937_64 rng(keyedHash(seed, MLP, static_cast<std::uint64_t>(member)));
      Model::Net net;
      net.d = d;
      net.h = h;
      std::normal_distribution<double> nd(0.0, 1.0);
      net.w1.resize(d * h);
      net.b1.assign(h, 0.0);
      net.w2.resize(h * h);
      net.b2.assign(h, 0.0);
      net.w3.resize(h);
      for (double& v : net.w1) v = nd(rng) / std::sqrt(static_cast<double>(std::max<Size>(d, 1)));
      for (double& v : net.w2) v = nd(rng) / std::sqrt(static_cast<double>(h));
      for (double& v : net.w3) v = nd(rng) / std::sqrt(static_cast<double>(h));
      // Adam state over the flat parameter vector
      std::vector<double*> params;
      for (double& v : net.w1) params.push_back(&v);
      for (double& v : net.b1) params.push_back(&v);
      for (double& v : net.w2) params.push_back(&v);
      for (double& v : net.b2) params.push_back(&v);
      for (double& v : net.w3) params.push_back(&v);
      params.push_back(&net.b3);
      const Size np = params.size();
      std::vector<double> mom(np, 0.0), vel(np, 0.0), grad(np, 0.0);
      const double lr = 1e-3, beta1 = 0.9, beta2 = 0.999, eps = 1e-8;
      Size step = 0;
      std::vector<Size> order(m);
      std::iota(order.begin(), order.end(), 0);
      const Size batch = 256;
      std::vector<double> a1(h), a2(h), d2(h), d1(h);
      // class weights normalized so that the total weight is m
      double wsum = 0.0;
      for (Size i = 0; i < m; ++i) wsum += (y[i] > 0 ? c_pos : c_neg);
      const double wnorm = wsum > 0 ? static_cast<double>(m) / wsum : 1.0;
      for (int epoch = 0; epoch < epochs; ++epoch)
      {
        std::shuffle(order.begin(), order.end(), rng);
        for (Size start = 0; start < m; start += batch)
        {
          const Size end = std::min(m, start + batch);
          std::fill(grad.begin(), grad.end(), 0.0);
          for (Size bi = start; bi < end; ++bi)
          {
            const Size i = order[bi];
            const double* x = X.row(rows[i]);
            for (Size k = 0; k < h; ++k)
            {
              double s = net.b1[k];
              for (Size j = 0; j < d; ++j) s += net.w1[j * h + k] * x[j];
              a1[k] = std::tanh(s);
            }
            for (Size k = 0; k < h; ++k)
            {
              double s = net.b2[k];
              for (Size j = 0; j < h; ++j) s += net.w2[j * h + k] * a1[j];
              a2[k] = std::tanh(s);
            }
            double o = net.b3;
            for (Size k = 0; k < h; ++k) o += net.w3[k] * a2[k];
            const double margin = 1.0 - y[i] * o;
            if (margin <= 0.0) continue;
            const double cw = (y[i] > 0 ? c_pos : c_neg) * wnorm;
            const double dout = -2.0 * cw * margin * y[i]; // d loss / d o
            // parameter offsets
            const Size o_w1 = 0, o_b1 = d * h, o_w2 = o_b1 + h, o_b2 = o_w2 + h * h, o_w3 = o_b2 + h, o_b3 = o_w3 + h;
            for (Size k = 0; k < h; ++k)
            {
              grad[o_w3 + k] += dout * a2[k];
              d2[k] = dout * net.w3[k] * (1.0 - a2[k] * a2[k]);
            }
            grad[o_b3] += dout;
            for (Size j = 0; j < h; ++j)
            {
              double s = 0.0;
              for (Size k = 0; k < h; ++k)
              {
                grad[o_w2 + j * h + k] += a1[j] * d2[k];
                s += net.w2[j * h + k] * d2[k];
              }
              d1[j] = s * (1.0 - a1[j] * a1[j]);
            }
            for (Size k = 0; k < h; ++k) grad[o_b2 + k] += d2[k];
            for (Size j = 0; j < d; ++j)
            {
              for (Size k = 0; k < h; ++k) grad[o_w1 + j * h + k] += x[j] * d1[k];
            }
            for (Size k = 0; k < h; ++k) grad[o_b1 + k] += d1[k];
          }
          ++step;
          const double nb = static_cast<double>(end - start);
          const double c1 = 1.0 - std::pow(beta1, static_cast<double>(step));
          const double c2 = 1.0 - std::pow(beta2, static_cast<double>(step));
          for (Size p = 0; p < np; ++p)
          {
            const double g = grad[p] / nb + decay * *params[p];
            mom[p] = beta1 * mom[p] + (1.0 - beta1) * g;
            vel[p] = beta2 * vel[p] + (1.0 - beta2) * g * g;
            *params[p] -= lr * (mom[p] / c1) / (std::sqrt(vel[p] / c2) + eps);
          }
        }
      }
      model.nets.push_back(std::move(net));
    }
    return model;
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
  if (opt.learner == "mlp" && (opt.mlp_members < 1 || opt.mlp_hidden < 1 || opt.mlp_epochs < 1))
  {
    bad_option("mlp_members / mlp_hidden / mlp_epochs (>= 1)", StringUtils::toStr(opt.mlp_members));
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
    if (pos.size() < 3 || neg.size() < 3)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "RESET: too few positive stems or training decoys for 3-fold cross-validation (at least 3 each).",
        StringUtils::toStr(pos.size()) + " / " + StringUtils::toStr(neg.size()));
    }
    auto fit = [&](const std::vector<Size>& rows, const std::vector<double>& y, double cp, double cn, std::uint64_t seed)
    {
      Model fm = opt.learner == "mlp" ? trainMlp(X, rows, y, cp, cn, seed, opt.mlp_members, opt.mlp_hidden, opt.mlp_epochs)
                                      : trainSvm(X, rows, y, cp, cn);
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
        const Model fm = fit(rows, y, cp, cn, keyedHash(salt, MLP, static_cast<std::uint64_t>(round * 100 + f)));
        std::vector<double>& hs = held_scores[static_cast<Size>(t)];
        for (Size r : held_index[f]) hs.push_back(fm.score(X.row(train_rows[r])));
      }
      catch (...)
      {
        errors[static_cast<Size>(t)] = std::current_exception();
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
    model = fit(train_rows, train_y, grid[best].first, grid[best].second,
                keyedHash(salt, MLP, static_cast<std::uint64_t>(round * 100 + 99)));
    if (opt.representative == "model")
    {
      // the stem's PSM with the best score under the new model represents it from now on (label-blind: a function of
      // the stem's PSMs and the model)
      std::vector<std::exception_ptr> rep_errors(m);
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
              standardize(in.features[row], x.data());
              const double sc = model.score(x.data());
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
    for (Size k = 0; k < m; ++k)
    {
      cur[k] = model.score(X.row(k));
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
    for (Size r = 0; r < dec.size(); ++r)
    {
      blk_val.push_back(dec[r] ? 1.0 : 0.0);
      blk_w.push_back(1.0);
      blk_len.push_back(1);
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
