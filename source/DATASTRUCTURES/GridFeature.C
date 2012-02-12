// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/GridFeature.h>

using namespace std;

namespace OpenMS
{

	GridFeature::GridFeature(const BaseFeature& feature, Size map_index,
													 Size feature_index) : 
		feature_(feature),
		map_index_(map_index), feature_index_(feature_index), annotations_()
	{
		const vector<PeptideIdentification>& peptides = 
			feature.getPeptideIdentifications();
		for (vector<PeptideIdentification>::const_iterator pep_it = 
					 peptides.begin(); pep_it != peptides.end(); ++pep_it)
		{
			if (pep_it->getHits().empty()) continue; // shouldn't be the case
			annotations_.insert(pep_it->getHits()[0].getSequence());
		}
	}

	GridFeature::~GridFeature()
	{
	}

	const BaseFeature& GridFeature::getFeature() const
	{
		return feature_;
	}

	Size GridFeature::getMapIndex() const
	{
		return map_index_;
	}

	Size GridFeature::getFeatureIndex() const
	{
		return feature_index_;
	}

	Int GridFeature::getID() const
	{
		return (Int)feature_index_;
	}

	const set<AASequence>& GridFeature::getAnnotations() const
	{
		return annotations_;
	}

	DoubleReal GridFeature::getRT() const
	{ 
		return feature_.getRT();
}

	DoubleReal GridFeature::getMZ() const
	{ 
		return feature_.getMZ();
	}

}
