// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

	void GaussFilter::init(DoubleReal sigma, DoubleReal spacing)
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
