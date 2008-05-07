// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>

namespace OpenMS
{

	FeatureGroupingAlgorithmLabeled::FeatureGroupingAlgorithmLabeled()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmLabeled");
		
		defaults_.insert("",PairMatcher().getParameters());
		
		defaultsToParam_();
	}

	FeatureGroupingAlgorithmLabeled::~FeatureGroupingAlgorithmLabeled()
	{
	}

	void FeatureGroupingAlgorithmLabeled::group(const std::vector< FeatureMap<> >& maps, ConsensusMap<>& out)
	{
		//initialize PairMatcher
    PairMatcher pm;
    
    pm.setParameters(param_.copy("",true));
    
    //run it
    pm.run(maps[0]);
    
    //store the result
    const PairMatcher::PairVectorType& pairs = pm.getBestPairs();
    for (UInt i=0; i<pairs.size(); ++i)
    {
    	UInt i1 = pairs[i].getFirst().getMetaValue(11);
    	UInt i2 = pairs[i].getSecond().getMetaValue(11);
    	ConsensusMap<>::ConsensusElementType c(1,i1,pairs[i].getFirst(),2,i2,pairs[i].getSecond());
    	out.push_back(c);
    }
    
	}

} //namespace 
