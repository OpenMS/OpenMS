// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAScorer.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    /// splitmix64: deterministic on every platform (std distributions are not).
    class Random
    {
    public:
      explicit Random(std::uint64_t seed) : state_(seed) {}
      std::uint64_t next()
      {
        std::uint64_t z = (state_ += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
      }
      /// Uniform in [0, n) without modulo bias.
      std::uint64_t below(std::uint64_t n)
      {
        const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max() - std::numeric_limits<std::uint64_t>::max() % n;
        std::uint64_t x;
        do { x = next(); } while (x >= limit);
        return x % n;
      }
      double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }

    private:
      std::uint64_t state_;
    };

    /// Draws k distinct entries (or k with replacement if the pool is smaller) by a partial Fisher-Yates shuffle of the
    /// pool in place; the pool stays a permutation, so every draw is uniform.
    void draw(std::vector<std::uint32_t>& pool, std::size_t k, Random& rng, std::vector<std::uint32_t>& out)
    {
      out.resize(k);
      const std::size_t n = pool.size();
      if (n < k)
      {
        for (std::size_t i = 0; i < k; ++i) out[i] = pool[rng.below(n)];
        return;
      }
      for (std::size_t i = 0; i < k; ++i)
      {
        std::swap(pool[i], pool[i + rng.below(n - i)]);
        out[i] = pool[i];
      }
    }

    /// Cholesky solve of the SPD system A x = b (A k x k row-major, overwritten).
    std::vector<double> solveSPD(std::vector<double>& a, std::vector<double> b, std::size_t k)
    {
      for (std::size_t j = 0; j < k; ++j)
      {
        double d = a[j * k + j];
        for (std::size_t p = 0; p < j; ++p) d -= a[j * k + p] * a[j * k + p];
        if (!(d > 0.0)) throw std::runtime_error("WhitenedNetwork: the decoy covariance is not positive definite.");
        d = std::sqrt(d);
        a[j * k + j] = d;
        for (std::size_t i = j + 1; i < k; ++i)
        {
          double s = a[i * k + j];
          for (std::size_t p = 0; p < j; ++p) s -= a[i * k + p] * a[j * k + p];
          a[i * k + j] = s / d;
        }
      }
      for (std::size_t i = 0; i < k; ++i)
      {
        for (std::size_t p = 0; p < i; ++p) b[i] -= a[i * k + p] * b[p];
        b[i] /= a[i * k + i];
      }
      for (std::size_t i = k; i-- > 0;)
      {
        for (std::size_t p = i + 1; p < k; ++p) b[i] -= a[p * k + i] * b[p];
        b[i] /= a[i * k + i];
      }
      return b;
    }

    double median(std::vector<double> v)
    {
      const std::size_t n = v.size(), h = n / 2;
      std::nth_element(v.begin(), v.begin() + h, v.end());
      const double upper = v[h];
      if (n % 2) return upper;
      return 0.5 * (upper + *std::max_element(v.begin(), v.begin() + h));
    }

    constexpr std::size_t CHUNK = 32; ///< Rows per gradient chunk; chunks are summed in order (thread independence).
  }

  // ---------------------------------------------------------------------------------------------------------------
  // WhitenedNetwork

  WhitenedNetwork::WhitenedNetwork(std::size_t dimension, std::size_t hidden, std::size_t embed, std::uint64_t seed) :
    dimension_(dimension), hidden_(hidden), embed_(embed)
  {
    if (dimension == 0 || hidden == 0 || embed == 0) throw std::invalid_argument("WhitenedNetwork: sizes must be positive.");
    parameters_.resize(hidden * dimension + hidden + embed * hidden + embed);
    Random rng(seed);
    const double b1 = 1.0 / std::sqrt(static_cast<double>(dimension)), b2 = 1.0 / std::sqrt(static_cast<double>(hidden));
    std::size_t p = 0;
    for (std::size_t i = 0; i < hidden * dimension + hidden; ++i) parameters_[p++] = (2.0 * rng.uniform() - 1.0) * b1;
    for (std::size_t i = 0; i < embed * hidden + embed; ++i) parameters_[p++] = (2.0 * rng.uniform() - 1.0) * b2;
  }

  void WhitenedNetwork::embed(const float* z, double* h) const
  {
    const double* w1 = parameters_.data();
    const double* c1 = w1 + hidden_ * dimension_;
    const double* w2 = c1 + hidden_;
    const double* c2 = w2 + embed_ * hidden_;
    double a[256];
    std::vector<double> big;
    double* act = a;
    if (hidden_ > 256)
    {
      big.resize(hidden_);
      act = big.data();
    }
    for (std::size_t j = 0; j < hidden_; ++j)
    {
      double s = c1[j];
      for (std::size_t i = 0; i < dimension_; ++i) s += w1[j * dimension_ + i] * z[i];
      act[j] = s > 0.0 ? s : 0.0;
    }
    for (std::size_t i = 0; i < dimension_; ++i) h[i] = z[i];
    for (std::size_t e = 0; e < embed_; ++e)
    {
      double s = c2[e];
      for (std::size_t j = 0; j < hidden_; ++j) s += w2[e * hidden_ + j] * act[j];
      h[dimension_ + e] = std::tanh(s);
    }
  }

  std::vector<double> WhitenedNetwork::readout(const std::vector<double>& positive, const std::vector<double>& negative, std::size_t k,
                                               double lambda, double& objective)
  {
    const std::size_t np = positive.size() / k, nn = negative.size() / k;
    if (np < 1 || nn < 2) throw std::invalid_argument("WhitenedNetwork::readout: needs >= 1 positive and >= 2 negative rows.");
    std::vector<double> mp(k, 0.0), mn(k, 0.0);
    for (std::size_t r = 0; r < np; ++r)
      for (std::size_t i = 0; i < k; ++i) mp[i] += positive[r * k + i];
    for (std::size_t r = 0; r < nn; ++r)
      for (std::size_t i = 0; i < k; ++i) mn[i] += negative[r * k + i];
    for (std::size_t i = 0; i < k; ++i)
    {
      mp[i] /= static_cast<double>(np);
      mn[i] /= static_cast<double>(nn);
    }
    std::vector<double> s(k * k, 0.0), c(k);
    for (std::size_t r = 0; r < nn; ++r)
    {
      for (std::size_t i = 0; i < k; ++i) c[i] = negative[r * k + i] - mn[i];
      for (std::size_t i = 0; i < k; ++i)
        for (std::size_t j = 0; j <= i; ++j) s[i * k + j] += c[i] * c[j];
    }
    for (std::size_t i = 0; i < k; ++i)
    {
      for (std::size_t j = 0; j <= i; ++j)
      {
        s[i * k + j] /= static_cast<double>(nn - 1);
        s[j * k + i] = s[i * k + j];
      }
      s[i * k + i] += lambda;
    }
    std::vector<double> dmu(k);
    for (std::size_t i = 0; i < k; ++i) dmu[i] = mp[i] - mn[i];
    auto w = solveSPD(s, dmu, k);
    objective = 0.0;
    for (std::size_t i = 0; i < k; ++i) objective += dmu[i] * w[i];
    return w;
  }

  double WhitenedNetwork::objective(const std::vector<const float*>& positive, const std::vector<const float*>& negative, double lambda,
                                    std::vector<double>* gradient, int threads) const
  {
    const std::size_t k = outputs(), np = positive.size(), nn = negative.size(), n = np + nn;
    std::vector<double> h(n * k), act(gradient ? n * hidden_ : 0);
    const double* w1 = parameters_.data();
    const double* c1 = w1 + hidden_ * dimension_;
    const double* w2 = c1 + hidden_;
    const double* c2 = w2 + embed_ * hidden_;
    const auto row = [&](std::size_t r) { return r < np ? positive[r] : negative[r - np]; };

#pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(n); ++r)
    {
      const float* z = row(static_cast<std::size_t>(r));
      double* hr = h.data() + r * k;
      double local[256];
      std::vector<double> big;
      double* a = gradient ? act.data() + r * hidden_ : local;
      if (!gradient && hidden_ > 256)
      {
        big.resize(hidden_);
        a = big.data();
      }
      for (std::size_t j = 0; j < hidden_; ++j)
      {
        double s = c1[j];
        for (std::size_t i = 0; i < dimension_; ++i) s += w1[j * dimension_ + i] * z[i];
        a[j] = s > 0.0 ? s : 0.0;
      }
      for (std::size_t i = 0; i < dimension_; ++i) hr[i] = z[i];
      for (std::size_t e = 0; e < embed_; ++e)
      {
        double s = c2[e];
        for (std::size_t j = 0; j < hidden_; ++j) s += w2[e * hidden_ + j] * a[j];
        hr[dimension_ + e] = std::tanh(s);
      }
    }

    // Readout on the batch.
    const std::vector<double> hp(h.begin(), h.begin() + np * k), hn(h.begin() + np * k, h.end());
    double obj = 0.0;
    const auto w = readout(hp, hn, k, lambda, obj);
    if (!gradient) return obj;

    // d obj / d h: positives 2 w / P; negatives -2 w / N - 2 w (w' c_j) / (N - 1), c_j the centred negative.
    std::vector<double> mn(k, 0.0);
    for (std::size_t r = 0; r < nn; ++r)
      for (std::size_t i = 0; i < k; ++i) mn[i] += hn[r * k + i];
    for (auto& v : mn) v /= static_cast<double>(nn);
    std::vector<double> coef(n);
    for (std::size_t r = 0; r < np; ++r) coef[r] = 2.0 / static_cast<double>(np);
    for (std::size_t r = 0; r < nn; ++r)
    {
      double wc = 0.0;
      for (std::size_t i = 0; i < k; ++i) wc += w[i] * (hn[r * k + i] - mn[i]);
      coef[np + r] = -2.0 / static_cast<double>(nn) - 2.0 * wc / static_cast<double>(nn - 1);
    }

    // Back-propagate per chunk of rows into its own buffer, then sum the chunks in order.
    const std::size_t chunks = (n + CHUNK - 1) / CHUNK, size = parameters_.size();
    std::vector<double> partial(chunks * size, 0.0);
#pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t ch = 0; ch < static_cast<std::ptrdiff_t>(chunks); ++ch)
    {
      double* g = partial.data() + ch * size;
      double* g_w1 = g;
      double* g_c1 = g_w1 + hidden_ * dimension_;
      double* g_w2 = g_c1 + hidden_;
      double* g_c2 = g_w2 + embed_ * hidden_;
      std::vector<double> gu2(embed_), ga(hidden_);
      for (std::size_t r = ch * CHUNK; r < std::min(n, (ch + 1) * CHUNK); ++r)
      {
        const float* z = row(r);
        const double* hr = h.data() + r * k;
        const double* a = act.data() + r * hidden_;
        for (std::size_t e = 0; e < embed_; ++e)
        {
          const double t = hr[dimension_ + e];
          gu2[e] = coef[r] * w[dimension_ + e] * (1.0 - t * t);
          g_c2[e] += gu2[e];
        }
        std::fill(ga.begin(), ga.end(), 0.0);
        for (std::size_t e = 0; e < embed_; ++e)
        {
          if (gu2[e] == 0.0) continue;
          for (std::size_t j = 0; j < hidden_; ++j)
          {
            g_w2[e * hidden_ + j] += gu2[e] * a[j];
            ga[j] += w2[e * hidden_ + j] * gu2[e];
          }
        }
        for (std::size_t j = 0; j < hidden_; ++j)
        {
          if (!(a[j] > 0.0)) continue;
          g_c1[j] += ga[j];
          for (std::size_t i = 0; i < dimension_; ++i) g_w1[j * dimension_ + i] += ga[j] * z[i];
        }
      }
    }
    gradient->assign(size, 0.0);
    for (std::size_t ch = 0; ch < chunks; ++ch)
      for (std::size_t p = 0; p < size; ++p) (*gradient)[p] += partial[ch * size + p];
    return obj;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // CrossFitScorer

  std::uint64_t CrossFitScorer::fnv1a(const std::string& s)
  {
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char b : s) h = (h ^ b) * 0x100000001b3ULL;
    return h;
  }

  CrossFitScorer::Competition CrossFitScorer::pooledPairedTDC(const Input& input, const std::vector<float>& score,
                                                              const std::vector<std::uint32_t>& rows)
  {
    std::uint32_t n_precursors = 0, n_pairs = 0;
    for (std::uint32_t r : rows)
    {
      n_precursors = std::max(n_precursors, input.precursor[r] + 1);
      n_pairs = std::max(n_pairs, input.pair[r] + 1);
    }
    // Precursor winners: highest score, ties to the lower row.
    constexpr std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> best(n_precursors, NONE);
    for (std::uint32_t r : rows)
    {
      auto& b = best[input.precursor[r]];
      if (b == NONE || score[r] > score[b] || (score[r] == score[b] && r < b)) b = r;
    }
    Competition out;
    for (std::uint32_t b : best)
      if (b != NONE) out.winners.push_back(b);
    std::sort(out.winners.begin(), out.winners.end());

    // Per-fold standardisation by the decoy winners' median and 1.4826 MAD.
    std::uint8_t n_folds = 0;
    for (std::uint32_t r : out.winners) n_folds = std::max<std::uint8_t>(n_folds, input.fold[r] + 1);
    std::vector<double> centre(n_folds, 0.0), spread(n_folds, 0.0);
    std::vector<bool> present(n_folds, false);
    for (std::uint8_t f = 0; f < n_folds; ++f)
    {
      std::vector<double> d;
      for (std::uint32_t r : out.winners)
      {
        if (input.fold[r] != f) continue;
        present[f] = true;
        if (input.decoy[r]) d.push_back(score[r]);
      }
      if (!present[f]) continue;
      if (d.empty()) throw std::invalid_argument("pooledPairedTDC: a fold has no decoy winners.");
      centre[f] = median(d);
      for (auto& v : d) v = std::abs(v - centre[f]);
      spread[f] = 1.4826 * median(d);
      if (!(spread[f] > 0.0)) throw std::invalid_argument("pooledPairedTDC: a fold has zero decoy spread.");
    }
    const auto standard = [&](std::uint32_t r) { return (score[r] - centre[input.fold[r]]) / spread[input.fold[r]]; };

    // Best member of each pair (ties to the decoy).
    std::vector<std::uint32_t> pair_best(n_pairs, NONE);
    for (std::uint32_t r : out.winners)
    {
      auto& b = pair_best[input.pair[r]];
      if (b == NONE) { b = r; continue; }
      const double sr = standard(r), sb = standard(b);
      if (sr > sb || (sr == sb && input.decoy[r] && !input.decoy[b])) b = r;
    }
    std::vector<std::pair<double, std::uint32_t>> ranked;
    for (std::uint32_t b : pair_best)
      if (b != NONE) ranked.emplace_back(standard(b), b);
    std::sort(ranked.begin(), ranked.end(), [&](const auto& x, const auto& y)
    {
      if (x.first != y.first) return x.first > y.first;
      if (input.decoy[x.second] != input.decoy[y.second]) return input.decoy[x.second] > input.decoy[y.second];
      return x.second < y.second;
    });
    out.ranked.resize(ranked.size());
    out.q.resize(ranked.size());
    double targets = 0.0, decoys = 0.0;
    for (std::size_t i = 0; i < ranked.size(); ++i)
    {
      out.ranked[i] = ranked[i].second;
      if (input.decoy[ranked[i].second]) decoys += 1.0;
      else targets += 1.0;
      out.q[i] = std::min(1.0, (decoys + 1.0) / std::max(targets, 1.0));
    }
    for (std::size_t i = out.q.size(); i-- > 1;) out.q[i - 1] = std::min(out.q[i - 1], out.q[i]);
    return out;
  }

  std::vector<float> CrossFitScorer::calibrateTail(const Input& input, const std::vector<float>& score)
  {
    const std::size_t n = input.rows();
    std::vector<std::uint32_t> kept;
    for (std::uint32_t r = 0; r < n; ++r)
      if (score[r] > -1e29f) kept.push_back(r);
    const auto winners = pooledPairedTDC(input, score, kept).winners;
    std::vector<float> out(n, NOT_KEPT);
    std::uint8_t n_folds = 0;
    for (std::uint32_t r : kept) n_folds = std::max<std::uint8_t>(n_folds, input.fold[r] + 1);
    for (std::uint8_t f = 0; f < n_folds; ++f)
    {
      std::vector<double> d;
      for (std::uint32_t r : winners)
        if (input.fold[r] == f && input.decoy[r]) d.push_back(score[r]);
      const std::size_t m = d.size();
      if (m == 0) continue;
      if (m < 10) throw std::invalid_argument("calibrateTail: a fold has fewer than 10 decoy winners.");
      std::sort(d.begin(), d.end());
      const std::size_t top = std::max<std::size_t>(m / 100, 10);
      const double t = d[m - top];
      double beta = 0.0;
      for (std::size_t i = m - top; i < m; ++i) beta += d[i] - t;
      beta = std::max(beta / static_cast<double>(top), 1e-9);
      const auto neg_log_survival = [&](std::size_t i) { return -std::log(static_cast<double>(m - i + 1) / static_cast<double>(m + 1)); };
      // Knots: unique decoy scores up to t, each with the survival of its first (lowest) index.
      std::vector<double> xs, ys;
      for (std::size_t i = 0; i <= m - top; ++i)
      {
        if (!xs.empty() && d[i] == xs.back()) continue;
        xs.push_back(d[i]);
        ys.push_back(neg_log_survival(i));
      }
      const double y_t = neg_log_survival(m - top);
      const double slope = xs.size() > 1 ? (ys[1] - ys[0]) / std::max(xs[1] - xs[0], 1e-12) : 0.0;
      for (std::uint32_t r : kept)
      {
        if (input.fold[r] != f) continue;
        const double v = score[r];
        double c;
        if (v > t) c = y_t + (v - t) / beta;
        else if (v < xs.front()) c = xs.size() > 1 ? ys[0] + (v - xs[0]) * slope : ys[0];
        else
        {
          const std::size_t j = static_cast<std::size_t>(std::upper_bound(xs.begin(), xs.end(), v) - xs.begin());
          if (j >= xs.size()) c = ys.back();
          else c = ys[j - 1] + (ys[j] - ys[j - 1]) * (v - xs[j - 1]) / (xs[j] - xs[j - 1]);
        }
        out[r] = static_cast<float>(c);
      }
    }
    return out;
  }

  namespace
  {
    struct Fit
    {
      WhitenedNetwork net;
      std::vector<double> w;
      int epochs = 0;
      double criterion = 0.0;
    };

    /// Adam on -objective with early stopping on the validation objective; final readout on the training rows.
    Fit fit(const std::vector<float>& z, std::size_t d, const std::vector<std::uint32_t>& rows, const std::vector<std::uint8_t>& label,
            const std::vector<bool>& valid, std::uint64_t seed, const CrossFitScorer::Options& o)
    {
      std::vector<std::uint32_t> tp, tn, vp, vn;
      for (std::size_t i = 0; i < rows.size(); ++i)
      {
        auto& target = valid[i] ? (label[i] ? vp : vn) : (label[i] ? tp : tn);
        target.push_back(rows[i]);
      }
      if (tp.empty() || tn.size() < 2 || vp.empty() || vn.size() < 2)
      {
        throw std::invalid_argument("CrossFitScorer: too few training or validation rows (" + std::to_string(tp.size()) + "/" +
                                    std::to_string(tn.size()) + "/" + std::to_string(vp.size()) + "/" + std::to_string(vn.size()) + ").");
      }
      Random rng(seed * 0x9e3779b97f4a7c15ULL + 1);
      Fit f{WhitenedNetwork(d, o.hidden, o.embed, rng.next()), {}, 0, 0.0};
      const std::size_t size = f.net.parameters().size(), k = f.net.outputs();
      std::vector<double> m(size, 0.0), v(size, 0.0), g, best_state;
      const double b1 = 0.9, b2 = 0.999, eps = 1e-8;
      double best = -std::numeric_limits<double>::infinity();
      int stale = 0;
      long step = 0;
      std::vector<std::uint32_t> bp, bn;
      std::vector<const float*> pp, pn;
      const auto pointers = [&](const std::vector<std::uint32_t>& r, std::vector<const float*>& out)
      {
        out.resize(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) out[i] = z.data() + static_cast<std::size_t>(r[i]) * d;
      };
      std::vector<const float*> vpp, vnp;
      pointers(vp, vpp);
      pointers(vn, vnp);
      for (int epoch = 0; epoch < o.max_epochs; ++epoch)
      {
        for (int s = 0; s < o.steps_per_epoch; ++s)
        {
          draw(tp, o.batch_positive, rng, bp);
          draw(tn, o.batch_negative, rng, bn);
          pointers(bp, pp);
          pointers(bn, pn);
          f.net.objective(pp, pn, o.lambda, &g, o.threads);
          ++step;
          const double c1 = 1.0 - std::pow(b1, static_cast<double>(step)), c2 = 1.0 - std::pow(b2, static_cast<double>(step));
          auto& p = f.net.parameters();
          for (std::size_t i = 0; i < size; ++i)
          {
            const double gi = -g[i]; // loss = -objective
            m[i] = b1 * m[i] + (1.0 - b1) * gi;
            v[i] = b2 * v[i] + (1.0 - b2) * gi * gi;
            p[i] -= o.learning_rate * (m[i] / c1) / (std::sqrt(v[i] / c2) + eps);
          }
        }
        f.epochs = epoch + 1;
        const double crit = f.net.objective(vpp, vnp, o.lambda, nullptr, o.threads);
        if (crit > best + 1e-4)
        {
          best = crit;
          stale = 0;
          best_state = f.net.parameters();
        }
        else if (++stale >= o.patience) break;
      }
      f.net.parameters() = best_state;
      f.criterion = best;
      // Final readout on the training rows.
      std::vector<double> hp(tp.size() * k), hn(tn.size() * k);
#pragma omp parallel for num_threads(o.threads) schedule(static)
      for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(tp.size()); ++i) f.net.embed(z.data() + static_cast<std::size_t>(tp[i]) * d, hp.data() + i * k);
#pragma omp parallel for num_threads(o.threads) schedule(static)
      for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(tn.size()); ++i) f.net.embed(z.data() + static_cast<std::size_t>(tn[i]) * d, hn.data() + i * k);
      double obj = 0.0;
      f.w = WhitenedNetwork::readout(hp, hn, k, o.lambda, obj);
      return f;
    }

    /// Full and linear scores of the given rows.
    void scoreRows(const Fit& f, const std::vector<float>& z, std::size_t d, const std::vector<std::uint32_t>& rows,
                   std::vector<float>& full, std::vector<float>* linear, int threads)
    {
      const std::size_t k = f.net.outputs();
#pragma omp parallel num_threads(threads)
      {
        std::vector<double> h(k);
#pragma omp for schedule(static)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(rows.size()); ++i)
        {
          const std::uint32_t r = rows[i];
          f.net.embed(z.data() + static_cast<std::size_t>(r) * d, h.data());
          double lin = 0.0, s = 0.0;
          for (std::size_t j = 0; j < d; ++j) lin += h[j] * f.w[j];
          s = lin;
          for (std::size_t j = d; j < k; ++j) s += h[j] * f.w[j];
          full[r] = static_cast<float>(s);
          if (linear) (*linear)[r] = static_cast<float>(lin);
        }
      }
    }
  }

  CrossFitScorer::Result CrossFitScorer::run(const Input& in, const Options& o, const std::function<void(const std::string&)>& log)
  {
    const std::size_t n = in.rows(), d = in.dimension;
    if (d == 0 || in.features.size() != n * d || in.pair.size() != n || in.decoy.size() != n || in.fold.size() != n ||
        in.pair_hash.size() != n)
    {
      throw std::invalid_argument("CrossFitScorer: inconsistent input sizes.");
    }
    if (o.threads < 1) throw std::invalid_argument("CrossFitScorer: threads must be >= 1.");
    std::uint8_t n_folds = 0;
    for (std::uint8_t f : in.fold) n_folds = std::max<std::uint8_t>(n_folds, f + 1);
    if (n_folds < 2) throw std::invalid_argument("CrossFitScorer: needs at least 2 folds.");
    {
      std::vector<int> pair_fold(*std::max_element(in.pair.begin(), in.pair.end()) + 1, -1);
      for (std::size_t r = 0; r < n; ++r)
      {
        int& f = pair_fold[in.pair[r]];
        if (f >= 0 && f != in.fold[r]) throw std::invalid_argument("CrossFitScorer: a pair spans two folds.");
        f = in.fold[r];
      }
    }
    for (const auto& s : in.seeds)
      if (s.row >= n) throw std::invalid_argument("CrossFitScorer: seed row out of range.");

    Result res;
    const auto say = [&](const std::string& line)
    {
      res.log.push_back(line);
      if (log) log(line);
    };
    std::vector<float> raw(n, 0.0f), linear(n, 0.0f), round_score(n, 0.0f);
    std::vector<float> z(n * d);
    for (std::uint8_t k = 0; k < n_folds; ++k)
    {
      std::vector<std::uint32_t> train, test;
      for (std::uint32_t r = 0; r < n; ++r) (in.fold[r] == k ? test : train).push_back(r);
      // Standardise by the training folds (population SD; constant columns keep SD 1).
      std::vector<double> mu(d, 0.0), sd(d, 0.0);
      for (std::uint32_t r : train)
        for (std::size_t j = 0; j < d; ++j) mu[j] += in.features[r * d + j];
      for (auto& v : mu) v /= static_cast<double>(train.size());
      for (std::uint32_t r : train)
        for (std::size_t j = 0; j < d; ++j)
        {
          const double x = in.features[r * d + j] - mu[j];
          sd[j] += x * x;
        }
      for (auto& v : sd)
      {
        v = std::sqrt(v / static_cast<double>(train.size()));
        if (v == 0.0) v = 1.0;
      }
#pragma omp parallel for num_threads(o.threads) schedule(static)
      for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(n); ++r)
        for (std::size_t j = 0; j < d; ++j) z[r * d + j] = static_cast<float>((in.features[r * d + j] - mu[j]) / sd[j]);

      // Round 0: seed rows of the training folds.
      std::vector<std::uint32_t> rows;
      std::vector<std::uint8_t> label;
      std::vector<bool> valid;
      for (const auto& s : in.seeds)
      {
        if (in.fold[s.row] == k) continue;
        rows.push_back(s.row);
        label.push_back(s.positive);
        valid.push_back(s.group_hash % 100 < o.validation_percent);
      }
      Fit model = fit(z, d, rows, label, valid, o.seed + 100 + k, o);
      std::string line = "fold " + std::to_string(k) + ": r0 " + std::to_string(rows.size()) + " " + std::to_string(model.epochs) +
                         " ep crit " + std::to_string(model.criterion);
      for (std::size_t r = 0; r < o.self_q.size(); ++r)
      {
        scoreRows(model, z, d, train, round_score, nullptr, o.threads);
        const auto comp = pooledPairedTDC(in, round_score, train);
        rows.clear();
        label.clear();
        valid.clear();
        std::size_t positives = 0;
        for (std::size_t i = 0; i < comp.ranked.size(); ++i)
        {
          if (comp.q[i] <= o.self_q[r] && !in.decoy[comp.ranked[i]])
          {
            rows.push_back(comp.ranked[i]);
            ++positives;
          }
        }
        for (std::uint32_t w : comp.winners)
          if (in.decoy[w]) rows.push_back(w);
        label.assign(rows.size(), 0);
        std::fill(label.begin(), label.begin() + positives, 1);
        for (std::uint32_t row : rows) valid.push_back(in.pair_hash[row] % 100 < o.validation_percent);
        model = fit(z, d, rows, label, valid, o.seed + 200 + 10 * k + r + 1, o);
        line += " | r" + std::to_string(r + 1) + " pos " + std::to_string(positives) + " neg " + std::to_string(rows.size() - positives) +
                " " + std::to_string(model.epochs) + " ep crit " + std::to_string(model.criterion);
      }
      scoreRows(model, z, d, test, raw, &linear, o.threads);
      say(line);
    }

    // W rule: per precursor the window with the largest linear part (ties: lower row) keeps its full score.
    std::uint32_t n_precursors = 0;
    for (std::uint32_t p : in.precursor) n_precursors = std::max(n_precursors, p + 1);
    constexpr std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> keep(n_precursors, NONE);
    for (std::uint32_t r = 0; r < n; ++r)
    {
      auto& b = keep[in.precursor[r]];
      if (b == NONE || linear[r] > linear[b]) b = r;
    }
    std::vector<float> kept(n, NOT_KEPT);
    for (std::uint32_t b : keep)
      if (b != NONE) kept[b] = raw[b];

    res.score = calibrateTail(in, kept);
    res.raw_score = std::move(raw);
    res.linear = std::move(linear);
    res.q_value.assign(n, 1.0);
    std::vector<std::uint32_t> all(n);
    std::iota(all.begin(), all.end(), 0u);
    const auto comp = pooledPairedTDC(in, res.score, all);
    for (std::size_t i = 0; i < comp.ranked.size(); ++i) res.q_value[comp.ranked[i]] = comp.q[i];
    return res;
  }
}
