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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>

using namespace std;

namespace OpenMS
{
  ComplementMarker::ComplementMarker()
    : PeakMarker()
  {
		setName(ComplementMarker::getProductName());
    defaults_.setValue("tolerance", 1);
    defaults_.setValue("marks", 1);
		defaultsToParam_();
  }

  ComplementMarker::ComplementMarker(const ComplementMarker& source)
    : PeakMarker(source)
  {
  }
  
  ComplementMarker::~ComplementMarker()
  {
  }

  ComplementMarker& ComplementMarker::operator = (const ComplementMarker& source)
  {
		if (this != &source)
		{
    	PeakMarker::operator=(source);
		}
    return *this;
  }

/*
  map<double,bool> ComplementMarker::operator()( MSSpectrum< DPeak<1> >& spec) const
  {
    // how often a peak needs to be marked to be returned
    double marks = (double)param_.getValue("marks");
    double parentmass = spec.getPrecursorPeak().getPosition()[0];
    double tolerance = (double)param_.getValue("tolerance");
    map<double,int> matching_b_y_ions;
    spec.getContainer().sortByPosition();
    int j = spec.size() -1;
    for (uint i = 0; i < spec.size(); ++i)
    {
      while ( j >= 0 && spec.getContainer()[j].getPosition()[0] > (parentmass-spec.getContainer()[i].getPosition()[0]) + tolerance ){
        j--;
      }
      // just takes the first matching ion; todo take all
      if ( j >= 0 && fabs(spec.getContainer()[i].getPosition()[0] + spec.getContainer()[j].getPosition()[0] - parentmass ) < tolerance ) 
      {
        matching_b_y_ions[spec.getContainer()[i].getPosition()[0]]++;
        matching_b_y_ions[spec.getContainer()[j].getPosition()[0]]++;
        j--;
      }
    }
    map<double,bool> result;
    for ( map<double,int>::const_iterator cmit = matching_b_y_ions.begin(); cmit != matching_b_y_ions.end(); ++cmit )
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
