// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

namespace OpenMS
{

   GaussFilter::GaussFilter()
	  : ProgressLogger(),
	    DefaultParamHandler("GaussFilter"),
			coeffs_(),
	    sigma_(0.1),
	    spacing_(0.01)
  {
  	//Parameter settings
  	defaults_.setValue("gaussian_width",0.2, "Use a gaussian filter width which has approximately the same width as your mass peaks (FWHM in m/z).");
		defaults_.setValue("ppm_tolerance", 10.0 , "Gaussian width, depending on the m/z position.\nThe higher the value, the wider the peak and therefore the wider the gaussian.");
		defaults_.setValue("use_ppm_tolerance", "false", "If true, instead of the gaussian_width value, the ppm_tolerance is used. The gaussian is calculated in each step anew, so this is much slower.");
		defaults_.setValidStrings("use_ppm_tolerance", StringList::create("true,false"));
    defaultsToParam_();
  }

  GaussFilter::~GaussFilter()
  {
  }

  void GaussFilter::updateMembers_() 
  {
    sigma_ = (DoubleReal)param_.getValue("gaussian_width") / 8.0;
		Size number_of_points_right = (Size)(ceil(4*sigma_ / spacing_))+1;
    coeffs_.resize(number_of_points_right);
    coeffs_[0] = 1.0/(sigma_ * sqrt(2.0 * Constants::PI));

    for (Size i=1; i < number_of_points_right; i++)
    {
    	coeffs_[i] = 1.0/(sigma_ * sqrt(2.0 * Constants::PI)) * exp(-((i*spacing_)*(i*spacing_)) / (2 * sigma_ * sigma_));
    }
#ifdef DEBUG_FILTERING
    std::cout << "Coeffs: " << std::endl;
    for (Size i=0; i < number_of_points_right; i++)
    {
        std::cout << i*spacing_ << ' ' << coeffs_[i] << std::endl;
    }
#endif
  }

}
