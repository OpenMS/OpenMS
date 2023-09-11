// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

namespace OpenMS
{
  double ContinuousWaveletTransformNumIntegration::integrate_
    (const std::vector<double> & processed_input,
    double spacing_data,
    int index)
  {
    double v = 0.;
    int half_width = (int)wavelet_.size();
    int index_in_data = (int)floor((half_width * spacing_) / spacing_data);
    int offset_data_left = ((index - index_in_data) < 0) ? 0 : (index - index_in_data);
    int offset_data_right = ((index + index_in_data) > (int)processed_input.size() - 1) ? (int)processed_input.size() - 2 : (index + index_in_data);

    // integrate from i until offset_data_left
    {
      int index_w_r = 0;
      for (int i = index; i > offset_data_left; --i)
      {
        int index_w_l = (int)Math::round(((index - (i - 1)) * spacing_data) / spacing_);
        // we could also use:
        // v += spacing_data / 2. * (...), but this can be factored out (see below) for faster computation
        v += (processed_input[i] * wavelet_[index_w_r] + processed_input[i - 1] * wavelet_[index_w_l]);
        index_w_r = index_w_l;
      }
    }

    // integrate from i+1 until offset_data_right
    {
      int index_w_l = 0;
      for (int i = index; i < offset_data_right; ++i)
      {
        int index_w_r = (int)Math::round((((i + 1) - index) * spacing_data) / spacing_);
        v += (processed_input[i + 1] * wavelet_[index_w_r] + processed_input[i] * wavelet_[index_w_l]);
        index_w_l = index_w_r;
      }
    }

    // multiply by (spacing_data / 2.), but change order for better numerical stability
    return v / 2./ sqrt(scale_) * spacing_data;
  }

  void ContinuousWaveletTransformNumIntegration::init(double scale, double spacing)
  {
    // will set members for scale_ and spacing_
    ContinuousWaveletTransform::init(scale, spacing);
    int number_of_points = (int)(ceil(5 * scale_ / spacing_)) + 1;
    wavelet_.reserve(number_of_points);
    wavelet_.push_back(1.);

    const double spacing_scale = spacing_ / scale_;
    for (int i = 1; i < number_of_points; ++i)
    {
      wavelet_.push_back(marr_(i * spacing_scale));
    }

#ifdef DEBUG_PEAK_PICKING
    std::cout << "WAVELET" << std::endl;
    for (int i = 0; i < number_of_points; i++)
    {
      std::cout << i << ' '  <<   i * spacing_ << " " << wavelet_[i] << std::endl;
    }
#endif

  }

}
