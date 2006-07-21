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

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>

namespace OpenMS
{
double ContinuousWaveletTransform::getInterpolatedValue_(double x, RawDataPointConstIterator it_left)
{
    // Interpolate between the point to the left and the point to the right.
    double left_position = it_left->getPosition()[mz_dim_];
    double right_position = (it_left+1)->getPosition()[mz_dim_];
    double d=(x-left_position)/(right_position-left_position);

    return ((it_left+1)->getIntensity()*d+it_left->getIntensity()*(1-d));
}

void ContinuousWaveletTransform::init(double scale, double spacing)
{
    scale_ = scale;
    spacing_=spacing;
}
}
