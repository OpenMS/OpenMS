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
// $Id: IsotopeDiffFilter.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  // Bern have a different one  
  //const String IsotopeDiffFilter::info_ = "returns total intensity of peak pairs that could result from isotope peaks";
    

  IsotopeDiffFilter::IsotopeDiffFilter()
    : FilterFunctor()
  { 
		name_ = IsotopeDiffFilter::getName();
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37);
		param_ = defaults_;
  }
  
  IsotopeDiffFilter::IsotopeDiffFilter(const IsotopeDiffFilter& source)
    : FilterFunctor(source)
  {
  }
  
  IsotopeDiffFilter& IsotopeDiffFilter::operator=(const IsotopeDiffFilter& source)
  {
    FilterFunctor::operator=(source);
    return *this;
  }
  
  IsotopeDiffFilter::~IsotopeDiffFilter()
  {
  }

/*
  String IsotopeDiffFilter::info() const
  {
    return info_;
  }

  vector<double> IsotopeDiffFilter::operator() (const ClusterSpectrum& cspec)
  {
    double tolerance = (double)param_.getValue("tolerance");;
    double isodiff = 0;
    //iterate over all peaks
    for (int i = 0; i < (int)cspec.getSpec().size(); ++i)
    {
      for (uint j = 1; i+j < cspec.getSpec().size() ; ++j)
      {
        if ( fabs( cspec.getSpec().getContainer()[i+j].getPosition()[0] - cspec.getSpec().getContainer()[i].getPosition()[0] + 1 ) < tolerance && cspec.getSpec().getContainer()[i-j].getIntensity() < cspec.getSpec().getContainer()[i].getIntensity())
        {
          isodiff+= cspec.getSpec().getContainer()[i].getIntensity() +
	  cspec.getSpec().getContainer()[i+j].getIntensity();
        }
        else if ( fabs( cspec.getSpec().getContainer()[i+j].getPosition()[0] - cspec.getSpec().getContainer()[i].getPosition()[0] ) > 1 + tolerance )
        {
          break;
        }
      }
    }
    vector<double> result;
    result.push_back(isodiff);
    return result;
  }
	*/
}
