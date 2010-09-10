// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/GridFeature.h>

namespace OpenMS
{

GridFeature::GridFeature(const std::vector<std::vector<BaseFeature> >& input_maps,Size map_index,Size feature_index) : GridElement(input_maps[map_index][feature_index].getRT(),input_maps[map_index][feature_index].getMZ()),map_index_(map_index),feature_index_(feature_index), input_maps_(input_maps) {
	// TODO Auto-generated constructor stub

}

GridFeature::~GridFeature() {
	// TODO Auto-generated destructor stub
}

BaseFeature GridFeature::getFeature() const
{
	return input_maps_[map_index_][feature_index_];
}

Size GridFeature::getMapIndex()
{
	return map_index_;
}
Size GridFeature::getFeatureIndex()
{
	return feature_index_;
}

Int GridFeature::getID()
{
	return feature_index_;
}

}
