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
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>

using namespace std;

namespace OpenMS
{
 
  IsotopeMarker::IsotopeMarker()
    : PeakMarker()
  {
		setName(IsotopeMarker::getProductName());
    defaults_.setValue("marks", 1, "How often a peak must be marked to be reported", false);
    defaults_.setValue("mz_variation", 0.1, "variation in m/z direction", false);
    defaults_.setValue("in_variation", 0.5, "variation in intensity", false);
		defaultsToParam_();
  }

  IsotopeMarker::IsotopeMarker(const IsotopeMarker& source)
    : PeakMarker(source)
  {
  }
  
  IsotopeMarker::~IsotopeMarker()
  {
  }

  IsotopeMarker& IsotopeMarker::operator = (const IsotopeMarker& source)
  {
		if (this != &source)
		{
    	PeakMarker::operator=(source);
		}
    return *this;
  }

/*
  map<double,bool> IsotopeMarker::operator()(MSSpectrum<>& spec)const
  {
    double mzvariation = (double)param_.getValue("mz_variation");
    double invariation = (double)param_.getValue("in_variation");
    uint marks = (unsigned int)param_.getValue("marks");
    spec.getContainer().sortByPosition();
    map<double,uint> isotopemarks ; // possible isotopes
    for (uint i = 0; i < spec.size(); ++i)
    {
      double mz = spec.getContainer()[i].getPosition()[0];
      double intensity = spec.getContainer()[i].getIntensity();
      uint j = i+1;
      vector<pair<double,double> > isotopes = SpectrumGenerator::instance()->isotopepeaks(mz,intensity);
      while ( j < spec.getContainer().size() && spec.getContainer()[j].getPosition()[0] <= mz + 3 + mzvariation )
      {
        double curmz = spec.getContainer()[j].getPosition()[0];
        double curIntensity = spec.getContainer()[j].getIntensity();
        uint iso = (uint)(curmz - mz + 0.499999 );
        if ( iso > 0 && curmz - mz - iso > mzvariation ) 
        {
          ++j;
          continue;
        }
        if ( fabs(isotopes[iso].second-curIntensity) < invariation*isotopes[iso].second )
        {
          isotopemarks[mz]++;
          isotopemarks[curmz]++;
        }
        ++j;
      }
    }
    map<double,bool> result;
    for ( map<double,uint>::const_iterator cmit = isotopemarks.begin(); cmit != isotopemarks.end(); ++cmit )
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
