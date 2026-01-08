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
    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
    }

    double intExptotal = norm(intExp.begin(), intExp.end());
    double intTheorTotal = norm(theorint.begin(), theorint.end());
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::dotProd(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

  double manhattanScoring(std::vector<double> intExp, std::vector<double> theorint)
  {

    for (unsigned int i = 0; i < intExp.size(); ++i)
    {
      intExp[i] = sqrt(intExp[i]);
      theorint[i] = sqrt(theorint[i]);
      //std::transform(intExp.begin(), intExp.end(), intExp.begin(), sqrt);
      //std::transform(theorint.begin(), theorint.end(), theorint.begin(), sqrt);
    }

    double intExptotal = std::accumulate(intExp.begin(), intExp.end(), 0.0);
    double intTheorTotal = std::accumulate(theorint.begin(), theorint.end(), 0.0);
    OpenSwath::normalize(intExp, intExptotal, intExp);
    OpenSwath::normalize(theorint, intTheorTotal, theorint);
    double score2 = OpenSwath::manhattanDist(intExp.begin(), intExp.end(), theorint.begin());
    return score2;
  }

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
