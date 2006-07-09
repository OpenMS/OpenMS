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
// $Id: ComplementFilter.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <cmath>

using namespace std;

namespace OpenMS
{
  // todo charge state 2?

  ComplementFilter::ComplementFilter()
		:	FilterFunctor()
  { 
		name_ = ComplementFilter::getName();
    //value from Bioinformatics, Bern 2004
		defaults_.setValue("tolerance", 0.37f);
		param_ = defaults_;
  }

	ComplementFilter::ComplementFilter(const ComplementFilter& source)
		:	FilterFunctor(source)
	{

	}

	ComplementFilter& ComplementFilter::operator = (const ComplementFilter& source)
	{
		FilterFunctor::operator = (source);
		return *this;
	}

  ComplementFilter::~ComplementFilter()
  {
  }

/*
  vector<double> ComplementFilter::operator() (const ClusterSpectrum& cspec)
  {
    double tolerance = (double)param_.getValue("tolerance");
    double parentmass = cspec.getParentMass();
    double result = 0;
    uint j = cspec.getSpec().size() - 1;
    for (uint i = 0; i < cspec.getSpec().size() && i <= j; ++i)
    {
      double sum = cspec.getSpec().getContainer()[i].getPosition()[0] + cspec.getSpec().getContainer()[j].getPosition()[0];
      if (fabs( sum-parentmass ) < tolerance )
      {
        result += cspec.getSpec().getContainer()[i].getIntensity() + cspec.getSpec().getContainer()[j].getIntensity();
      }
      else if ( sum < parentmass )
      {
        ++i;
      }
      else if ( sum > parentmass )
      {
        --j;
      }
      else 
      {
				// @todo another Exception?!?  Which one? ????
        //throw CanNotHappen(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
    }
    vector<double> res;
    res.push_back(result);
    return res;
  }
	*/
}
