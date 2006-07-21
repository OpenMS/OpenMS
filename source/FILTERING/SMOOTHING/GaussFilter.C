// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

namespace OpenMS
{
GaussFilter::GaussFilter(const Param& param) throw (Exception::InvalidValue)
{
    param_ = param;

    // if a smoothing parameter is missed in the param object the value should be substituted by a dv value
    DataValue dv;
    float kernel_width = 0.;

    dv = param_.getValue("GaussianWidth");
    if (dv.isEmpty() || dv.toString() == "")
        kernel_width = .8;
    else
        kernel_width = (float)dv;

    if (kernel_width <= 0)
    {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The kernel width should be greater than zero!",String(kernel_width));
    }

    // The standard deviation corresponds approximately kernel_width / 8
    sigma_ = kernel_width / 8.;

    //compute the filter kernel coefficients with at least 50 data points
    spacing_= 4*sigma_ / 50;
    init(sigma_,spacing_);
}

void GaussFilter::init(float sigma, float spacing)
{
    sigma_= sigma;
    spacing_=spacing;

    int number_of_points_right = (int)(ceil(4*sigma_ / spacing_))+1;
    coeffs_.resize(number_of_points_right);
    coeffs_[0] = 1.0/(sigma_ * sqrt(2.0 * M_PI));

    for (int i=1; i < number_of_points_right; i++)
    {
        coeffs_[i] = gauss_(i*spacing_);
    }
#ifdef DEBUG_FILTERING
    std::cout << "Coeffs: " << std::endl;
    for (int i=0; i < number_of_points_right; i++)
    {
        std::cout << i*spacing_ << ' ' << coeffs_[i] << std::endl;
    }
#endif

}
}
