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
// $Id: MarkerMower.C,v 1.4 2006/04/05 11:18:24 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>

using namespace std;

namespace OpenMS
{

  const String MarkerMower::info_ = "Removes all peaks that are not marked by inserted PeakMarkers";
  
  /**
  IsotopeMarker, ComplementMarker, NeutralLossMarker are used per default<br>
  */
  MarkerMower::MarkerMower()
    : PreprocessingFunctor()
  {
		name_ = MarkerMower::getName();
    //todo remove
    insertmarker(new IsotopeMarker());
    insertmarker(new ComplementMarker());
    insertmarker(new NeutralLossMarker());
  }

  MarkerMower::MarkerMower(const MarkerMower& source)
    : PreprocessingFunctor(source)
  {
		name_ = source.getName();
  }

  MarkerMower::~MarkerMower()
  {
  }

  MarkerMower& MarkerMower::operator=(const MarkerMower& source)
  {
    PreprocessingFunctor::operator=(source);
    return *this;
  }

  String MarkerMower::info() const
  {
    return info_;
  }

  void MarkerMower::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    map<double,int> marks;
    for ( vector<PeakMarker*>::const_iterator cvit = markers_.begin(); cvit != markers_.end(); ++cvit )
    {
      map<double,bool> marked = (**cvit)(spec);
      for ( map<double,bool>::const_iterator cmit = marked.begin(); cmit != marked.end(); ++cmit )
      {
        if ( cmit->second ) marks[cmit->first]++;
      }
    }
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); )
    {
      if ( marks[it->getPosition()[0]] > 0 )
      {
        ++it;
      }
      else 
      {
        it = spec.getContainer().erase(it);
      }
    }
    
  }
  
  /**
  violates FactoryProduct interface
  */
  void MarkerMower::insertmarker(PeakMarker* pm)
  {
    markers_.push_back(pm);
  }

}
