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
// $Maintainer: Florian Zeller $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSHCTRL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSHCTRL_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <vector>
#include <map>

namespace OpenMS
{
	/** 
		@brief A facade for various Superhirn FeatureFinder classes. Use FeatureFinderAlgorithmSH instead.
	
		@ingroup FeatureFinder
	*/
  typedef std::map<double, RawData*> MyMap;
  typedef std::vector<MyMap> Vec;
  
  class OPENMS_DLLAPI FeatureFinderAlgorithmSHCtrl
  {
    
  public:
    
    FeatureFinderAlgorithmSHCtrl() { }
    
    std::vector<Feature> extractPeaks(Vec datavec);
    
    void initParams(Param param);
    
  protected:
    
  };
  
}

#endif
