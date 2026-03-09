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
  double norm(T beg, T end);

  /// @cond

  // Explicit template instantiation declarations
  extern template double norm<std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator);

  extern template double norm<std::vector<double>::iterator>(
    std::vector<double>::iterator, std::vector<double>::iterator);

  /// @endcond

  /**
  @brief compute dotprod of vectors
  */
  template <typename Texp, typename Ttheo>
  double dotProd(Texp intExpBeg, Texp intExpEnd, Ttheo intTheo);

  /// @cond
  
  // Explicit template instantiation declarations (tell the compiler these exist)
  extern template double dotProd<std::vector<double>::const_iterator, std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator, std::vector<double>::const_iterator);
  
  extern template double dotProd<std::vector<float>::const_iterator, std::vector<float>::const_iterator>(
    std::vector<float>::const_iterator, std::vector<float>::const_iterator, std::vector<float>::const_iterator);
  
  extern template double dotProd<std::vector<int>::const_iterator, std::vector<int>::const_iterator>(
    std::vector<int>::const_iterator, std::vector<int>::const_iterator, std::vector<int>::const_iterator);

  /// @endcond
  
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
  double manhattanDist(Texp itExpBeg, Texp itExpEnd, Ttheo itTheo);

  /// @cond

  // Explicit template instantiation declarations
  extern template double manhattanDist<std::vector<double>::iterator, std::vector<double>::iterator>(
    std::vector<double>::iterator, std::vector<double>::iterator, std::vector<double>::iterator);

  extern template double manhattanDist<std::vector<double>::const_iterator, std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator, std::vector<double>::const_iterator);

  /// @endcond

  /**
    @brief manhattan scoring

    sqrt intensities
    normalize vector by TIC
    compute manhattan score
   */
  OPENSWATHALGO_DLLAPI double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint);


/**
  @brief compute pearson correlation of vector x and y
  @deprecated Use Math::pearsonCorrelationCoefficient instead
*/
  template <typename TInputIterator, typename TInputIteratorY>
  [[deprecated("Use Math::pearsonCorrelationCoefficient instead")]]
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

