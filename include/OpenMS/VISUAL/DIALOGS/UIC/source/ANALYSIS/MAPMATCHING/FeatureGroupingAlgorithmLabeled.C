// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>

namespace OpenMS
{

	FeatureGroupingAlgorithmLabeled::FeatureGroupingAlgorithmLabeled()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmLabeled");
		
		defaults_.insert("",LabeledPairFinder().getParameters());
		
		defaultsToParam_();
	}

	FeatureGroupingAlgorithmLabeled::~FeatureGroupingAlgorithmLabeled()
	{
	}

	void FeatureGroupingAlgorithmLabeled::group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out)
	{
		//check that the number of maps is ok
		if (maps.size()!=1) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Exactly one map must be given!");
		if (out.getFileDescriptions().size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Two file descriptions must be set in 'out'!");
		
		//initialize LabeledPairFinder
    LabeledPairFinder pm;
    pm.setParameters(param_.copy("",true));
    
    //convert to consensus map
		std::vector<ConsensusMap> input(1);
		ConsensusMap::convert(0,maps[0],input[0]);
		
		//run
		pm.run(input,out);
	}

} //namespace 
