// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski  $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>

#include <algorithm>
#include <numeric>
#include <functional>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace OpenSwath
{
  void normalize(
    const std::vector<double> & intensities,
    double normalizer,
    std::vector<double> & normalized_intensities)
  {
    //compute total intensities
    //normalize intensities
    normalized_intensities.resize(intensities.size());
    if (normalizer > 0)
    {
      std::transform(intensities.cbegin(), intensities.cend(), normalized_intensities.begin(),
                     [&normalizer](double val)
                     {
                      return val / normalizer;
                     });
    }
  }


  double dotprodScoring(std::vector<double> intExp, std::vector<double> theorint)
  {
    // Vectorized square roots
    Eigen::Map<Eigen::VectorXd> mapExp(intExp.data(), intExp.size());
    Eigen::Map<Eigen::VectorXd> mapTheor(theorint.data(), theorint.size());

    mapExp = mapExp.cwiseSqrt();
    mapTheor = mapTheor.cwiseSqrt();

    double intExptotal = norm(intExp.begin(), intExp.end());
    double intTheorTotal = norm(theorint.begin(), theorint.end());
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

  double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint)
  {
    // Vectorized square roots
    Eigen::Map<Eigen::VectorXd> mapExp(intExp.data(), intExp.size());
    Eigen::Map<Eigen::VectorXd> mapTheor(theorint.data(), theorint.size());

    mapExp = mapExp.cwiseSqrt();
    mapTheor = mapTheor.cwiseSqrt();

    double intExptotal = mapExp.sum();
    double intTheorTotal = mapTheor.sum();

    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

  // Template implementation for norm (only compiled in this .cpp file)
  template <typename T>
  double norm(T beg, T end)
  {
    if (beg == end) return 0.0;
    size_t size = std::distance(beg, end);
    using ValueType = typename std::iterator_traits<T>::value_type;
    Eigen::Map<const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>> v(&(*beg), size);
    return static_cast<double>(v.norm());
  }

  // Explicit template instantiation definitions for norm
  template OPENSWATHALGO_DLLAPI double norm<std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator);

  template OPENSWATHALGO_DLLAPI double norm<std::vector<double>::iterator>(
    std::vector<double>::iterator, std::vector<double>::iterator);

  // Template implementation for manhattanDist (only compiled in this .cpp file)
  template <typename Texp, typename Ttheo>
  double manhattanDist(Texp itExpBeg, Texp itExpEnd, Ttheo itTheo)
  {
    if (itExpBeg == itExpEnd) return 0.0;
    size_t size = std::distance(itExpBeg, itExpEnd);
    using ExpType = typename std::iterator_traits<Texp>::value_type;
    using TheoType = typename std::iterator_traits<Ttheo>::value_type;
    Eigen::Map<const Eigen::Matrix<ExpType, Eigen::Dynamic, 1>> v1(&(*itExpBeg), size);
    Eigen::Map<const Eigen::Matrix<TheoType, Eigen::Dynamic, 1>> v2(&(*itTheo), size);
    return static_cast<double>((v1.template cast<double>() - v2.template cast<double>()).cwiseAbs().sum());
  }

  // Explicit template instantiation definitions for manhattanDist
  template OPENSWATHALGO_DLLAPI double manhattanDist<std::vector<double>::iterator, std::vector<double>::iterator>(
    std::vector<double>::iterator, std::vector<double>::iterator, std::vector<double>::iterator);

  template OPENSWATHALGO_DLLAPI double manhattanDist<std::vector<double>::const_iterator, std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator, std::vector<double>::const_iterator);

  // Template implementation (only compiled in this .cpp file)
  template <typename Texp, typename Ttheo>
  double dotProd(Texp intExpBeg, Texp intExpEnd, Ttheo intTheo)
  {
    size_t size = std::distance(intExpBeg, intExpEnd);

    // Get the value types
    using ExpType = typename std::iterator_traits<Texp>::value_type;
    using TheoType = typename std::iterator_traits<Ttheo>::value_type;

    // Create appropriate Eigen maps based on the actual data types
    Eigen::Map<const Eigen::Matrix<ExpType, Eigen::Dynamic, 1>> vec1(&(*intExpBeg), size);
    Eigen::Map<const Eigen::Matrix<TheoType, Eigen::Dynamic, 1>> vec2(&(*intTheo), size);

    // Compute dot product and cast result to double
    return static_cast<double>(vec1.dot(vec2));
  }

  // Explicit template instantiation definitions (compile these specific versions)
  template OPENSWATHALGO_DLLAPI double dotProd<std::vector<double>::const_iterator, std::vector<double>::const_iterator>(
    std::vector<double>::const_iterator, std::vector<double>::const_iterator, std::vector<double>::const_iterator);

  template OPENSWATHALGO_DLLAPI double dotProd<std::vector<float>::const_iterator, std::vector<float>::const_iterator>(
    std::vector<float>::const_iterator, std::vector<float>::const_iterator, std::vector<float>::const_iterator);

  template OPENSWATHALGO_DLLAPI double dotProd<std::vector<int>::const_iterator, std::vector<int>::const_iterator>(
    std::vector<int>::const_iterator, std::vector<int>::const_iterator, std::vector<int>::const_iterator);

}
