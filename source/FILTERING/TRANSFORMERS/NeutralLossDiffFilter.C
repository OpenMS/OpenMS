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
// $Id: NeutralLossDiffFilter.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  // Bern have a different one  

  NeutralLossDiffFilter::NeutralLossDiffFilter()
    : FilterFunctor()
  { 
		name_ = NeutralLossDiffFilter::getName();
    //value from Bioinformatics, Bern 2004
    defaults_.setValue("tolerance", 0.37f);
		param_ = defaults_;
  }
  
  NeutralLossDiffFilter::NeutralLossDiffFilter(const NeutralLossDiffFilter& source)
    : FilterFunctor(source)
  { 
		name_ = source.getName();
  }
  
  NeutralLossDiffFilter& NeutralLossDiffFilter::operator=(const NeutralLossDiffFilter& source)
  { 
    FilterFunctor::operator=(source);
    return *this;
  }

  NeutralLossDiffFilter::~NeutralLossDiffFilter()
  {
  }
/*
  vector<double> NeutralLossDiffFilter::operator() (const ClusterSpectrum& cspec)
  {
    double tolerance = (double)param_.getValue("tolerance");
    double isodiff = 0;
    //iterate over all peaks
    for (int i = 0; i < (int)cspec.getSpec().size(); ++i)
    {
      for (int j = 1; i-j >= 0; ++j)
      {
        if ( fabs( cspec.getSpec().getContainer()[i-j].getPosition()[0] - 
	cspec.getSpec().getContainer()[i].getPosition()[0] - 18 ) < tolerance || // water 
	fabs( cspec.getSpec().getContainer()[i-j].getPosition()[0] - 
	cspec.getSpec().getContainer()[i].getPosition()[0] - 17 ) < tolerance ) //ammoniom
        {
          isodiff += cspec.getSpec().getContainer()[i-j].getIntensity()+
	  cspec.getSpec().getContainer()[i].getIntensity();
        }
        else if ( fabs( cspec.getSpec().getContainer()[i-j].getPosition()[0] - cspec.getSpec().getContainer()[i].getPosition()[0]  ) > 18 + tolerance )
        {
          break;
        }
      }
    }
    vector<double> result;
    result.push_back(isodiff);
    return result;
  }*/
}
