// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/SMOOTHING/GaussFilterAlgorithm.h>

namespace OpenMS
{

  GaussFilterAlgorithm::GaussFilterAlgorithm()  :
    coeffs_(),
    sigma_(0.1),
    spacing_(0.01),
    use_ppm_tolerance_(false),
    ppm_tolerance_(10.0)
  {
    initialize(sigma_ * 8, spacing_, ppm_tolerance_, use_ppm_tolerance_);
  }

  GaussFilterAlgorithm::~GaussFilterAlgorithm() = default;

  void GaussFilterAlgorithm::initialize(double gaussian_width, double spacing, double ppm_tolerance, bool use_ppm_tolerance)
  {
    spacing_ = spacing;
    use_ppm_tolerance_ = use_ppm_tolerance;
    ppm_tolerance_ = ppm_tolerance;
    sigma_ = gaussian_width / 8.0;
    Size number_of_points_right = (Size)(ceil(4 * sigma_ / spacing_)) + 1;
    coeffs_.resize(number_of_points_right);
    coeffs_[0] = 1.0 / (sigma_ * sqrt(2.0 * Constants::PI));

    for (Size i = 1; i < number_of_points_right; i++)
    {
      coeffs_[i] = 1.0 / (sigma_ * sqrt(2.0 * Constants::PI)) * exp(-((i * spacing_) * (i * spacing_)) / (2 * sigma_ * sigma_));
    }
#ifdef DEBUG_FILTERING
    std::cout << "Coeffs: " << std::endl;
    for (Size i = 0; i < number_of_points_right; i++)
    {
      std::cout << i * spacing_ << ' ' << coeffs_[i] << std::endl;
    }
#endif

  }

}
