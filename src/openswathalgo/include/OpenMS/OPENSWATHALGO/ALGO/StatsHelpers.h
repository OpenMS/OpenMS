// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg  $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <vector>
#include <cstddef>

namespace OpenSwath
{

  /**
    @brief Normalize intensities in vector by normalization_factor
  */
  OPENSWATHALGO_DLLAPI void normalize(const std::vector<double>& intensities, double normalization_factor, std::vector<double>& normalized_intensities);

  /**
  @brief compute the Euclidean norm of the vector
  */
  template <typename T>
  double norm(T beg, T end)
  {
    double res = 0.0;
    for (; beg != end; ++beg)
    {
      double tmp = *beg;
      res += tmp * tmp;
    }
    return sqrt(res);
  }

  /**
  @brief compute dotprod of vectors
  */
  template <typename Texp, typename Ttheo>
  double dotProd(Texp intExpBeg, Texp intExpEnd, Ttheo intTheo)
  {
    std::vector<double> res(std::distance(intExpBeg, intExpEnd));
    std::transform(intExpBeg, intExpEnd, intTheo, res.begin(), std::multiplies<double>());
    double sum = std::accumulate(res.begin(), res.end(), 0.);
    return sum;
  }

  /**
    @brief the dot product scoring

    sqrt data,
    normalize by vector norm
    compute dotprod
  */
  OPENSWATHALGO_DLLAPI double dotprodScoring(std::vector<double> intExp, std::vector<double> theorint);

  /**
    @brief compute manhattan distance between Exp and Theo
  */
  template <typename Texp, typename Ttheo>
  double manhattanDist(Texp itExpBeg, Texp itExpEnd, Ttheo itTheo)
  {
    double sum = 0.0;
    for (std::size_t i = 0; itExpBeg < itExpEnd; ++itExpBeg, ++itTheo, ++i)
    {
      double x = *itExpBeg - *itTheo;
      x = fabs(x);
      sum += x;
    }
    return sum;
  }

  /**
    @brief manhattan scoring

    sqrt intensities
    normalize vector by TIC
    compute manhattan score
   */
  OPENSWATHALGO_DLLAPI double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint);


/**
  @brief compute pearson correlation of vector x and y
*/
  template <typename TInputIterator, typename TInputIteratorY>
  typename std::iterator_traits<TInputIterator>::value_type cor_pearson(
    TInputIterator xBeg,
    TInputIterator xEnd,
    TInputIteratorY yBeg
    )
  {
    typedef typename std::iterator_traits<TInputIterator>::value_type value_type;
    value_type   m1, m2;
    value_type   s1, s2;
    value_type   corr;
    m1 = m2 = s1 = s2 = 0.0;
    corr = 0.0;
    ptrdiff_t n = std::distance(xBeg, xEnd);
    value_type nd = static_cast<value_type>(n);
    for (; xBeg != xEnd; ++xBeg, ++yBeg)
    {
      corr += *xBeg * *yBeg;
      m1 += *xBeg;
      m2 += *yBeg;
      s1 += *xBeg * *xBeg;
      s2 += *yBeg * *yBeg;
    }
    m1 /= nd;
    m2 /= nd;
    s1 -= m1 * m1 * nd;
    s2 -= m2 * m2 * nd;

    if (s1 < 1.0e-12 || s2 < 1.0e-12)
      return 0.0;
    else
    {
      corr -= m1 * m2 * (double)n;
      corr /= sqrt(s1 * s2);
      return corr;
    }
  }

  /**
    @brief functor to compute the mean and stddev of sequence using the std::foreach algorithm
  */
  class OPENSWATHALGO_DLLAPI mean_and_stddev
  {
    double m_, q_;
    unsigned long c_;
public:
    typedef double argument_type, result_type;
    mean_and_stddev() :
      m_(0.0), q_(0.0), c_(0u)
    {
    }

    void operator()(double sample)
    {
      double const delta = sample - m_;
      m_ += delta / ++c_;
      q_ += delta * (sample - m_);
    }

    double sample_variance() const
    {
      return (c_ > 1u) ? (q_ / (c_ - 1)) : 0;
    }

    double standard_variance() const
    {
      return (c_ > 1u) ? (q_ / c_) : 0;
    }

    double sample_stddev() const
    {
      return std::sqrt(sample_variance());
    }

    double standard_stddev() const
    {
      return std::sqrt(standard_variance());
    }

    double mean() const
    {
      return m_;
    }

    unsigned long count() const
    {
      return c_;
    }

    double variance() const
    {
      return sample_variance();
    }

    double stddev() const
    {
      return sample_stddev();
    }

    double operator()() const
    {
      return stddev();
    }

  };

} //end namespace OpenSwath

