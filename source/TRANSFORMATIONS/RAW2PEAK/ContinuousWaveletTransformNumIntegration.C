// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>

namespace OpenMS
{
  double ContinuousWaveletTransformNumIntegration::integrate_
  (const std::vector<double>& processed_input,
   double spacing_data,
   int index)
  {
    double v = 0.;
    int half_width = (int)wavelet_.size();
    int index_in_data = (int)floor((half_width*spacing_) / spacing_data);
    int offset_data_left = ((index - index_in_data) < 0) ? 0 : (index-index_in_data);
    int offset_data_right = ((index + index_in_data) > (int)processed_input.size()-1) ? (int)processed_input.size()-2 : (index+index_in_data);

    // integrate from i until offset_data_left
    for (int i = index; i > offset_data_left; --i)
    {
			int index_w_r = (int)Math::round(((index-i)*spacing_data)/spacing_);
			int index_w_l = (int)Math::round(((index-(i-1))*spacing_data)/spacing_);

      v += spacing_data / 2.*( processed_input[i]*wavelet_[index_w_r] + processed_input[i-1]*wavelet_[index_w_l] );
    }

    // integrate from i+1 until offset_data_right
    for (int i = index; i < offset_data_right; ++i)
    {
			int index_w_r = (int)Math::round((((i+1)-index)*spacing_data)/spacing_);
			int index_w_l = (int)Math::round(((i-index)*spacing_data)/spacing_);
			
      v += spacing_data / 2.*( processed_input[i+1]*wavelet_[index_w_r] + processed_input[i]*wavelet_[index_w_l]);
    }

    return v / sqrt(scale_);
  }


  void ContinuousWaveletTransformNumIntegration::init(double scale, double spacing)
  {
    ContinuousWaveletTransform::init(scale, spacing);
    int number_of_points_right = (int)(ceil(5*scale_/spacing_));
    int number_of_points = number_of_points_right + 1;
    wavelet_.resize(number_of_points);
    wavelet_[0] = 1.;

    for (int i=1; i<number_of_points; i++)
    {
      wavelet_[i] = marr_(i*spacing_/scale_ );
    }

#ifdef DEBUG_PEAK_PICKING
    std::cout << "WAVELET" << std::endl;
    for (int i=0; i<number_of_points; i++)
    {
      std::cout << i << ' '  <<   i*spacing_ << " " << wavelet_[i] << std::endl;
    }
#endif

  }


}
