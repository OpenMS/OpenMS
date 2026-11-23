// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAGBM.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    constexpr double EPS = 1e-12;
    constexpr double MIN_HESSIAN = 1e-3;

    struct Bin
    {
      double g = 0.0, h = 0.0;
      std::uint32_t n = 0;
    };

    struct Split
    {
      double gain = 0.0;
      int feature = -1;
      int bin = -1;
    };

    struct Leaf
    {
      std::size_t begin = 0, end = 0;  ///< Segment of the row index array.
      int node = 0;                    ///< Node id in the tree.
      double g = 0.0, h = 0.0;
      std::vector<Bin> hist;           ///< d x 256.
      Split split;
    };

    double score(double g, double h, double l2) { return g * g / (h + l2 + EPS); }

    std::uint64_t splitmix(std::uint64_t& state)
    {
      std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      return z ^ (z >> 31);
    }
  }

  void HistGBM::fit(const float* x, std::size_t d, const std::vector<std::uint32_t>& train, const std::vector<std::uint32_t>& valid,
                    const std::vector<std::uint8_t>& label, const std::vector<double>& weight, const Options& o)
  {
    trees.clear();
    const std::size_t nt = train.size();
    if (nt == 0) throw std::invalid_argument("HistGBM::fit: no training rows.");
    const int max_bins = std::clamp(o.max_bins, 2, 255);
    const int threads = std::max(1, o.threads);

    // Quantile thresholds per feature on (a deterministic stride sample of) the training rows. A value's bin is the
    // number of thresholds <= it; NaN takes the top bin. Left of split bin b means x < threshold[b].
    std::vector<std::vector<float>> thresholds(d);
    const std::size_t stride = std::max<std::size_t>(1, nt / 200000);
#pragma omp parallel for num_threads(threads) schedule(dynamic)
    for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(d); ++jj)
    {
      const std::size_t j = static_cast<std::size_t>(jj);
      std::vector<float> v;
      for (std::size_t i = 0; i < nt; i += stride)
      {
        const float val = x[static_cast<std::size_t>(train[i]) * d + j];
        if (!std::isnan(val)) v.push_back(val);
      }
      std::sort(v.begin(), v.end());
      std::vector<float> distinct;
      for (float val : v)
        if (distinct.empty() || val != distinct.back()) distinct.push_back(val);
      auto& t = thresholds[j];
      const std::size_t nv = max_bins - 1;  // value bins; the last one is shared with NaN
      if (distinct.size() <= nv)
      {
        for (std::size_t k = 1; k < distinct.size(); ++k) t.push_back(0.5f * distinct[k - 1] + 0.5f * distinct[k]);
      }
      else
      {
        for (std::size_t k = 1; k < nv; ++k)
        {
          const std::size_t p = std::max<std::size_t>(1, k * v.size() / nv);
          if (v[p - 1] == v[p]) continue;
          const float th = 0.5f * v[p - 1] + 0.5f * v[p];
          if (t.empty() || th > t.back()) t.push_back(th);
        }
      }
    }
    std::vector<std::uint8_t> bins(nt * d);
#pragma omp parallel for num_threads(threads) schedule(static)
    for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(d); ++jj)
    {
      const std::size_t j = static_cast<std::size_t>(jj);
      const auto& t = thresholds[j];
      for (std::size_t i = 0; i < nt; ++i)
      {
        const float val = x[static_cast<std::size_t>(train[i]) * d + j];
        bins[j * nt + i] = static_cast<std::uint8_t>(std::isnan(val) ? t.size() : std::upper_bound(t.begin(), t.end(), val) - t.begin());
      }
    }

    // Init: weighted log-odds.
    double wp = 0.0, wn = 0.0;
    for (std::uint32_t r : train) (label[r] ? wp : wn) += weight[r];
    init = std::log(std::max(wp, EPS) / std::max(wn, EPS));
    std::vector<double> raw(nt, init), g(nt), h(nt), raw_valid(valid.size(), init), loss_valid(valid.size());
    std::vector<std::uint32_t> idx(nt), tmp(nt);

    const auto validLoss = [&]()
    {
#pragma omp parallel for num_threads(threads) schedule(static)
      for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(valid.size()); ++i)
      {
        const std::uint32_t r = valid[i];
        const double z = raw_valid[i];
        // log(1 + exp(-y' z)) with y' in {-1, 1}
        const double m = label[r] ? -z : z;
        loss_valid[i] = weight[r] * (m > 0 ? m + std::log1p(std::exp(-m)) : std::log1p(std::exp(m)));
      }
      double s = 0.0;
      for (double l : loss_valid) s += l;
      return s;
    };
    double best_loss = valid.empty() ? 0.0 : validLoss();
    int since_best = 0;

    const auto buildHist = [&](Leaf& leaf)
    {
      leaf.hist.assign(d * 256, Bin());
#pragma omp parallel for num_threads(threads) schedule(dynamic)
      for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(d); ++jj)
      {
        const std::size_t j = static_cast<std::size_t>(jj);
        Bin* hb = leaf.hist.data() + j * 256;
        const std::uint8_t* col = bins.data() + j * nt;
        for (std::size_t k = leaf.begin; k < leaf.end; ++k)
        {
          const std::uint32_t i = idx[k];
          Bin& b = hb[col[i]];
          b.g += g[i];
          b.h += h[i];
          ++b.n;
        }
      }
    };
    const auto findSplit = [&](Leaf& leaf)
    {
      std::vector<Split> per(d);
      const std::size_t count = leaf.end - leaf.begin;
      const double parent = score(leaf.g, leaf.h, o.l2);
#pragma omp parallel for num_threads(threads) schedule(static)
      for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(d); ++jj)
      {
        const std::size_t j = static_cast<std::size_t>(jj);
        const Bin* hb = leaf.hist.data() + j * 256;
        const int nb = static_cast<int>(thresholds[j].size()) + 1;
        double gl = 0.0, hl = 0.0;
        std::size_t nl = 0;
        for (int b = 0; b + 1 < nb; ++b)
        {
          gl += hb[b].g;
          hl += hb[b].h;
          nl += hb[b].n;
          const std::size_t nr = count - nl;
          if (nl < o.min_samples_leaf) continue;
          if (nr < o.min_samples_leaf) break;
          const double hr = leaf.h - hl;
          if (hl < MIN_HESSIAN || hr < MIN_HESSIAN) continue;
          const double gain = score(gl, hl, o.l2) + score(leaf.g - gl, hr, o.l2) - parent;
          if (gain > per[j].gain) per[j] = {gain, static_cast<int>(j), b};
        }
      }
      leaf.split = Split();
      for (const auto& s : per)
        if (s.feature >= 0 && s.gain > leaf.split.gain) leaf.split = s;
    };

    for (int it = 0; it < o.max_iter; ++it)
    {
      for (std::size_t i = 0; i < nt; ++i)
      {
        const std::uint32_t r = train[i];
        const double p = 1.0 / (1.0 + std::exp(-raw[i]));
        g[i] = weight[r] * (p - (label[r] ? 1.0 : 0.0));
        h[i] = weight[r] * std::max(p * (1.0 - p), 1e-16);
      }
      std::iota(idx.begin(), idx.end(), 0u);
      std::vector<Node> tree(1);
      std::vector<Leaf> leaves(1);
      leaves[0].begin = 0;
      leaves[0].end = nt;
      for (std::size_t i = 0; i < nt; ++i)
      {
        leaves[0].g += g[i];
        leaves[0].h += h[i];
      }
      buildHist(leaves[0]);
      findSplit(leaves[0]);
      while (static_cast<int>(leaves.size()) < o.max_leaf_nodes)
      {
        int pick = -1;
        for (std::size_t l = 0; l < leaves.size(); ++l)
          if (leaves[l].split.feature >= 0 && (pick < 0 || leaves[l].split.gain > leaves[pick].split.gain)) pick = static_cast<int>(l);
        if (pick < 0) break;
        Leaf parent = std::move(leaves[pick]);
        const int j = parent.split.feature, sb = parent.split.bin;
        const std::uint8_t* col = bins.data() + static_cast<std::size_t>(j) * nt;
        // Stable partition of the segment.
        std::size_t nl = 0, nr = 0;
        for (std::size_t k = parent.begin; k < parent.end; ++k)
        {
          if (col[idx[k]] <= sb) idx[parent.begin + nl++] = idx[k];
          else tmp[nr++] = idx[k];
        }
        std::copy_n(tmp.begin(), nr, idx.begin() + parent.begin + nl);
        Leaf left, right;
        left.begin = parent.begin;
        left.end = right.begin = parent.begin + nl;
        right.end = parent.end;
        for (int b = 0; b <= sb; ++b)
        {
          left.g += parent.hist[static_cast<std::size_t>(j) * 256 + b].g;
          left.h += parent.hist[static_cast<std::size_t>(j) * 256 + b].h;
        }
        right.g = parent.g - left.g;
        right.h = parent.h - left.h;
        Leaf& small = nl <= nr ? left : right;
        Leaf& large = nl <= nr ? right : left;
        buildHist(small);
        large.hist = std::move(parent.hist);
        for (std::size_t k = 0; k < large.hist.size(); ++k)
        {
          large.hist[k].g -= small.hist[k].g;
          large.hist[k].h -= small.hist[k].h;
          large.hist[k].n -= small.hist[k].n;
        }
        Node& node = tree[parent.node];
        node.feature = j;
        node.threshold = thresholds[j][sb];
        node.left = static_cast<int>(tree.size());
        node.right = node.left + 1;
        left.node = node.left;
        right.node = node.right;
        tree.emplace_back();
        tree.emplace_back();
        findSplit(left);
        findSplit(right);
        leaves[pick] = std::move(left);
        leaves.insert(leaves.begin() + pick + 1, std::move(right));
      }
      for (auto& leaf : leaves)
      {
        const double v = -o.learning_rate * leaf.g / (leaf.h + o.l2 + EPS);
        tree[leaf.node].value = v;
        for (std::size_t k = leaf.begin; k < leaf.end; ++k) raw[idx[k]] += v;
      }
      trees.push_back(std::move(tree));
      if (leaves.size() == 1) break;  // nothing left to split
      if (!valid.empty())
      {
        const auto& t = trees.back();
#pragma omp parallel for num_threads(threads) schedule(static)
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(valid.size()); ++i)
        {
          const float* row = x + static_cast<std::size_t>(valid[i]) * d;
          int k = 0;
          while (t[k].feature >= 0) k = row[t[k].feature] < t[k].threshold ? t[k].left : t[k].right;
          raw_valid[i] += t[k].value;
        }
        const double loss = validLoss();
        if (loss < best_loss - o.tolerance * std::max(1.0, std::abs(best_loss)))
        {
          best_loss = loss;
          since_best = 0;
        }
        else if (++since_best >= o.patience) break;
      }
    }
  }

  double HistGBM::predict(const float* row) const
  {
    double z = init;
    for (const auto& t : trees)
    {
      int k = 0;
      while (t[k].feature >= 0) k = row[t[k].feature] < t[k].threshold ? t[k].left : t[k].right;
      z += t[k].value;
    }
    return z;
  }

  CrossFitScorer::Result gbmRescore(const CrossFitScorer::Input& in, const std::vector<float>& score, double hard_fraction,
                                    std::uint64_t seed, int threads, const std::function<void(const std::string&)>& log,
                                    int rounds, double positive_q)
  {
    const std::size_t n = in.rows(), d = in.dimension;
    if (score.size() != n || in.features.size() != n * d) throw std::invalid_argument("gbmRescore: inconsistent input sizes.");
    if (rounds < 1 || !(positive_q > 0.0) || positive_q > 1.0) throw std::invalid_argument("gbmRescore: rounds must be >= 1 and positive_q in (0, 1].");
    threads = std::max(1, threads);
    constexpr std::size_t MAX_NEG = 1500000;
    std::uint8_t n_folds = 0;
    for (std::uint8_t f : in.fold) n_folds = std::max<std::uint8_t>(n_folds, f + 1);
    std::uint32_t n_precursors = 0;
    for (std::uint32_t p : in.precursor) n_precursors = std::max(n_precursors, p + 1);
    std::vector<std::uint32_t> all(n);
    std::iota(all.begin(), all.end(), 0u);

    std::vector<float> s = score, next(n), prev;
    std::vector<std::uint8_t> label(n, 0);
    std::vector<double> weight(n, 1.0);
    HistGBM::Options go;
    go.threads = threads;
    const int ROUNDS = rounds;
    int round = 0;
    for (; round < ROUNDS; ++round)
    {
      CrossFitScorer::Competition comp;
      try
      {
        comp = CrossFitScorer::pooledPairedTDC(in, s, all);
      }
      catch (const std::invalid_argument& e)
      {
        // Degenerate margins (e.g. tied decoy winners on tiny inputs): keep the previous round's scores.
        if (round == 0) throw;
        if (log) log("GBM rescore r" + std::to_string(round + 1) + ": stopped (" + e.what() + "); keeping round " + std::to_string(round) + " scores");
        s = prev;
        break;
      }
      std::vector<std::uint32_t> pos, neg;
      for (std::size_t i = 0; i < comp.ranked.size(); ++i)
        if (comp.q[i] <= positive_q && !in.decoy[comp.ranked[i]]) pos.push_back(comp.ranked[i]);
      std::size_t hard_precursors = 0;
      std::vector<std::uint32_t> dw;
      for (std::uint32_t w : comp.winners)
        if (in.decoy[w]) dw.push_back(w);
      if (hard_fraction > 0.0)
      {
        std::sort(dw.begin(), dw.end(), [&](std::uint32_t a, std::uint32_t b) { return s[a] != s[b] ? s[a] > s[b] : a < b; });
        const std::size_t take = std::min(dw.size(), std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(hard_fraction * dw.size()))));
        std::vector<std::uint8_t> hard(n_precursors, 0);
        for (std::size_t i = 0; i < take; ++i) hard[in.precursor[dw[i]]] = 1;
        hard_precursors = take;
        for (std::uint32_t r = 0; r < n; ++r)
          if (in.decoy[r] && hard[in.precursor[r]]) neg.push_back(r);
      }
      else
      {
        neg = dw;
        hard_precursors = dw.size();
      }
      std::sort(pos.begin(), pos.end());
      std::string line = "GBM rescore r" + std::to_string(round + 1) + ": pos " + std::to_string(pos.size()) + " neg " +
                         std::to_string(neg.size()) + " (hard precursors " + std::to_string(hard_precursors) + ") trees ";
      if (pos.empty() || neg.empty())
      {
        if (log) log(line + "none (no positives or negatives); stopped");
        break;
      }
      for (std::uint8_t f = 0; f < n_folds; ++f)
      {
        std::vector<std::uint32_t> fp, fn;
        for (std::uint32_t r : pos)
          if (in.fold[r] != f) fp.push_back(r);
        for (std::uint32_t r : neg)
          if (in.fold[r] != f) fn.push_back(r);
        if (fn.size() > MAX_NEG)
        {
          std::uint64_t state = seed + 1000 * (round + 1) + f;
          for (std::size_t k = 0; k < MAX_NEG; ++k) std::swap(fn[k], fn[k + splitmix(state) % (fn.size() - k)]);
          fn.resize(MAX_NEG);
          std::sort(fn.begin(), fn.end());
        }
        const double wpos = fp.empty() ? 1.0 : static_cast<double>(fn.size()) / static_cast<double>(fp.size());
        std::vector<std::uint32_t> rows;
        rows.reserve(fp.size() + fn.size());
        std::merge(fp.begin(), fp.end(), fn.begin(), fn.end(), std::back_inserter(rows));
        std::vector<std::uint32_t> tr, va;
        for (std::uint32_t r : rows)
        {
          label[r] = !in.decoy[r];
          weight[r] = label[r] ? wpos : 1.0;
          (in.pair_hash[r] % 100 < 10 ? va : tr).push_back(r);
        }
        HistGBM model;
        if (!fp.empty() && !fn.empty() && !tr.empty()) model.fit(in.features.data(), d, tr, va, label, weight, go);
        else model.init = 0.0;
        line += (f ? "/" : "") + std::to_string(model.trees.size());
#pragma omp parallel for num_threads(threads) schedule(static)
        for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(n); ++r)
          if (in.fold[r] == f) next[r] = static_cast<float>(model.predict(in.features.data() + static_cast<std::size_t>(r) * d));
      }
      prev = std::move(s);
      s = next;
      if (log) log(line);
    }
    // The last round gets the same check as the others: a degenerate round (e.g. a fold whose model kept no tree, so
    // its scores are constant) is discarded.
    if (round == ROUNDS)
    {
      try
      {
        CrossFitScorer::pooledPairedTDC(in, s, all);
      }
      catch (const std::invalid_argument& e)
      {
        if (log) log("GBM rescore r" + std::to_string(ROUNDS) + ": discarded (" + e.what() + "); keeping round " + std::to_string(ROUNDS - 1) + " scores");
        s = prev;
      }
    }

    // Per precursor the best window (ties: lower row) keeps its score.
    constexpr std::uint32_t NONE = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> keep(n_precursors, NONE);
    for (std::uint32_t r = 0; r < n; ++r)
    {
      auto& b = keep[in.precursor[r]];
      if (b == NONE || s[r] > s[b]) b = r;
    }
    std::vector<float> kept(n, CrossFitScorer::NOT_KEPT);
    for (std::uint32_t b : keep)
      if (b != NONE) kept[b] = s[b];
    CrossFitScorer::Result res;
    res.score = CrossFitScorer::calibrateTail(in, kept);
    res.q_value.assign(n, 1.0);
    const auto comp = CrossFitScorer::pooledPairedTDC(in, res.score, all);
    for (std::size_t i = 0; i < comp.ranked.size(); ++i) res.q_value[comp.ranked[i]] = comp.q[i];
    return res;
  }
}
