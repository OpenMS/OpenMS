// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIARTRefine.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <set>

namespace OpenMS::OpenDIA
{
  std::vector<std::string> RTRefineModel::tokenize(const std::string& sequence)
  {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i < sequence.size() && sequence[i] == '.') ++i;
    // Length of a "(UniMod:...)" group starting at j, 0 if none.
    const auto mod_length = [&](std::size_t j) -> std::size_t
    {
      if (j >= sequence.size() || sequence[j] != '(') return 0;
      const auto close = sequence.find(')', j);
      return close == std::string::npos ? 0 : close - j + 1;
    };
    while (i < sequence.size())
    {
      const char c = sequence[i];
      if (std::isupper(static_cast<unsigned char>(c)))
      {
        const std::size_t m = mod_length(i + 1);
        out.push_back(sequence.substr(i, 1 + m));
        i += 1 + m;
      }
      else if (const std::size_t m = mod_length(i); m > 0)
      {
        out.push_back(sequence.substr(i, m));
        i += m;
      }
      else
      {
        ++i;  // other characters ('.', lowercase) carry no token
      }
    }
    return out;
  }

  RTRefineModel::RTRefineModel(const std::vector<std::string>& sequences)
  {
    std::set<std::string> vocabulary;
    std::set<char> alphabet;
    for (const auto& s : sequences)
    {
      for (auto& t : tokenize(s))
      {
        if (std::isupper(static_cast<unsigned char>(t[0]))) alphabet.insert(t[0]);
        vocabulary.insert(std::move(t));
      }
    }
    vocabulary_.assign(vocabulary.begin(), vocabulary.end());
    alphabet_.assign(alphabet.begin(), alphabet.end());
  }

  std::vector<double> RTRefineModel::features(const std::string& sequence, int charge, double mz) const
  {
    std::vector<double> f(dimension(), 0.0);
    const std::size_t nv = vocabulary_.size(), na = alphabet_.size();
    const auto letter = [&](char c) -> std::ptrdiff_t
    {
      const auto it = std::lower_bound(alphabet_.begin(), alphabet_.end(), c);
      return it != alphabet_.end() && *it == c ? it - alphabet_.begin() : -1;
    };
    const auto tokens = tokenize(sequence);
    std::size_t residues = 0;
    char first = 0, last = 0;
    for (const auto& t : tokens)
    {
      const auto it = std::lower_bound(vocabulary_.begin(), vocabulary_.end(), t);
      if (it != vocabulary_.end() && *it == t) f[static_cast<std::size_t>(it - vocabulary_.begin())] += 1.0;
      if (std::isupper(static_cast<unsigned char>(t[0])))
      {
        if (residues == 0) first = t[0];
        last = t[0];
        ++residues;
      }
    }
    if (residues > 0)
    {
      if (const auto a = letter(first); a >= 0) f[nv + static_cast<std::size_t>(a)] = 1.0;
      if (const auto a = letter(last); a >= 0) f[nv + na + static_cast<std::size_t>(a)] = 1.0;
    }
    const std::size_t o = nv + 2 * na;
    f[o] = static_cast<double>(residues);
    f[o + 1] = std::log(static_cast<double>(std::max<std::size_t>(residues, 1)));
    f[o + 2 + static_cast<std::size_t>(std::clamp(charge, 1, 4) - 1)] = 1.0;
    f[o + 6] = mz;
    return f;
  }

  bool choleskySolve(std::vector<double> a, std::vector<double> b, std::size_t d, std::vector<double>& x)
  {
    // In place: lower triangle of a becomes L with A = L L^T.
    for (std::size_t j = 0; j < d; ++j)
    {
      double s = a[j * d + j];
      for (std::size_t k = 0; k < j; ++k) s -= a[j * d + k] * a[j * d + k];
      if (!(s > 0.0)) return false;
      const double l = std::sqrt(s);
      a[j * d + j] = l;
      for (std::size_t i = j + 1; i < d; ++i)
      {
        double t = a[i * d + j];
        for (std::size_t k = 0; k < j; ++k) t -= a[i * d + k] * a[j * d + k];
        a[i * d + j] = t / l;
      }
    }
    // Forward (L z = b), then backward (L^T x = z).
    for (std::size_t i = 0; i < d; ++i)
    {
      double t = b[i];
      for (std::size_t k = 0; k < i; ++k) t -= a[i * d + k] * b[k];
      b[i] = t / a[i * d + i];
    }
    for (std::size_t ii = d; ii-- > 0;)
    {
      double t = b[ii];
      for (std::size_t k = ii + 1; k < d; ++k) t -= a[k * d + ii] * b[k];
      b[ii] = t / a[ii * d + ii];
    }
    x = std::move(b);
    return true;
  }

  double RidgeFit::predict(const double* row) const
  {
    double v = intercept;
    for (std::size_t k = 0; k < weights.size(); ++k)
      if (sd[k] > 0.0) v += weights[k] * (row[k] - mean[k]) / sd[k];
    return v;
  }

  RidgeFit fitRidge(const std::vector<double>& x, std::size_t d, const std::vector<std::size_t>& rows, const std::vector<double>& y, double alpha)
  {
    RidgeFit fit;
    fit.mean.assign(d, 0.0);
    fit.sd.assign(d, 0.0);
    fit.weights.assign(d, 0.0);
    const std::size_t m = rows.size();
    if (m == 0) return fit;
    double ymean = 0.0;
    for (std::size_t r = 0; r < m; ++r)
    {
      ymean += y[r];
      for (std::size_t k = 0; k < d; ++k) fit.mean[k] += x[rows[r] * d + k];
    }
    ymean /= static_cast<double>(m);
    for (auto& v : fit.mean) v /= static_cast<double>(m);
    for (std::size_t r = 0; r < m; ++r)
      for (std::size_t k = 0; k < d; ++k)
      {
        const double c = x[rows[r] * d + k] - fit.mean[k];
        fit.sd[k] += c * c;
      }
    for (auto& v : fit.sd) v = std::sqrt(v / static_cast<double>(m));  // population SD, as StandardScaler
    fit.intercept = ymean;
    // Normal equations on standardised features; a zero-SD feature stays 0 (its row and column only carry alpha).
    std::vector<double> a(d * d, 0.0), b(d, 0.0), z(d);
    for (std::size_t r = 0; r < m; ++r)
    {
      for (std::size_t k = 0; k < d; ++k) z[k] = fit.sd[k] > 0.0 ? (x[rows[r] * d + k] - fit.mean[k]) / fit.sd[k] : 0.0;
      const double yc = y[r] - ymean;
      for (std::size_t i = 0; i < d; ++i)
      {
        if (z[i] == 0.0) continue;
        b[i] += z[i] * yc;
        for (std::size_t j = 0; j <= i; ++j) a[i * d + j] += z[i] * z[j];
      }
    }
    for (std::size_t i = 0; i < d; ++i)
    {
      a[i * d + i] += alpha;
      for (std::size_t j = 0; j < i; ++j) a[j * d + i] = a[i * d + j];
    }
    if (!choleskySolve(std::move(a), std::move(b), d, fit.weights)) fit.weights.assign(d, 0.0);
    return fit;
  }

  RTRefineResult fitRTRefine(const std::vector<double>& x, std::size_t d, const std::vector<std::size_t>& train, const std::vector<double>& y,
                             const std::vector<std::uint32_t>& folds)
  {
    constexpr std::uint32_t K = 5;
    RTRefineResult result;
    result.training = train.size();
    // Out-of-fold predictions for one alpha.
    const auto cross_fit = [&](double alpha)
    {
      std::vector<double> pred(train.size(), 0.0);
      for (std::uint32_t f = 0; f < K; ++f)
      {
        std::vector<std::size_t> rows;
        std::vector<double> ys;
        for (std::size_t r = 0; r < train.size(); ++r)
          if (folds[r] != f) { rows.push_back(train[r]); ys.push_back(y[r]); }
        const auto fit = fitRidge(x, d, rows, ys, alpha);
        for (std::size_t r = 0; r < train.size(); ++r)
          if (folds[r] == f) pred[r] = fit.predict(x.data() + train[r] * d);
      }
      return pred;
    };
    double best = std::numeric_limits<double>::infinity();
    for (int g = 0; g < 13; ++g)
    {
      const double alpha = std::pow(10.0, -2.0 + 0.5 * g);
      auto pred = cross_fit(alpha);
      double mse = 0.0;
      for (std::size_t r = 0; r < train.size(); ++r) mse += (y[r] - pred[r]) * (y[r] - pred[r]);
      if (mse < best)
      {
        best = mse;
        result.alpha = alpha;
        result.out_of_fold = std::move(pred);
      }
    }
    const auto fit = fitRidge(x, d, train, y, result.alpha);
    const std::size_t n = d == 0 ? 0 : x.size() / d;
    result.correction.resize(n);
    for (std::size_t i = 0; i < n; ++i) result.correction[i] = fit.predict(x.data() + i * d);
    return result;
  }
}
