// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>

using namespace std;

namespace OpenMS
{

  NeutralLossMarker::NeutralLossMarker()
    : PeakMarker()
  {
		setName(NeutralLossMarker::getProductName());
		defaults_.setValue("marks", 1, "How often a peak must be marked to be reported", false);
    defaults_.setValue("tolerance", 0.2, "Tolerance in m/z direction", false);
		defaultsToParam_();
  }

  NeutralLossMarker::NeutralLossMarker(const NeutralLossMarker& source)
    : PeakMarker(source)
  {
  }
  
  NeutralLossMarker::~NeutralLossMarker()
  {
  }

  NeutralLossMarker& NeutralLossMarker::operator=(const NeutralLossMarker& source)
  {
		if (this != &source)
		{
    	PeakMarker::operator=(source);
		}
    return *this;
  }

/*
  map<double,bool> NeutralLossMarker::operator()( MSSpectrum<>& spec) const
  {
    // how often a peak needs to be marked to be returned
    double marks = 1;
    double tolerance = (double)param_.getValue("tolerance");
    map<double,int> ions_w_neutrallosses;
    spec.getContainer().sortByPosition();
    for (uint i = 0; i < spec.size(); ++i)
    {
      double mz = spec.getContainer()[i].getPosition()[0];
      double intensity = spec.getContainer()[i].getIntensity();
      int j = i - 1;
      while ( j >= 0 )
      {
        double curmz = spec.getContainer()[j].getPosition()[0];
        double curIntensity = spec.getContainer()[j].getIntensity();
       
        // check for peak thats a a water or ammonia away
        if ( fabs( mz-curmz-17 ) < tolerance || fabs( mz-curmz-18 ) < tolerance )
        {
          // neutral loss peak should be smaller
          if ( curIntensity < intensity )
          {
            ions_w_neutrallosses[mz]++;
            // neutral loss peak not marked
            //ions_w_neutrallosses[curmz]++;
          }
        }
        else if ( mz - curmz > 18.3 )
        {
          break;
        }
        --j;
      }
    }
    map<double,bool> result;
    for ( map<double,int>::const_iterator cmit = ions_w_neutrallosses.begin(); cmit != ions_w_neutrallosses.end(); ++cmit )
    {
      if ( cmit->second >= marks )
      {
        result.insert(make_pair(cmit->first,1));
      }
    }
    return result;
  }
*/
}
