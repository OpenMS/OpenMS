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

#ifndef OPENMS_DATASTRUCTURES_GRIDFEATURE_H
#define OPENMS_DATASTRUCTURES_GRIDFEATURE_H

#include <OpenMS/DATASTRUCTURES/GridElement.h>

namespace OpenMS
{
/**
 * @brief Used for QT feature linking based on geometric hashing.
 * A GridFeature can be stored in a HashGrid and points to a Feature.
 * @see HashGrid
 */
class OPENMS_DLLAPI GridFeature : public GridElement {
private:
	/**
	 * @brief Index of the Feature- or ConsensusMap
	 */
	Size map_index_;
	/**
	 * @brief Index of the Feature in the map
	 */
	Size feature_index_;
	/**
	 * @brief Reference to all FeatureMaps
	 */
	const std::vector<std::vector<BaseFeature> >& input_maps_;
public:
	/**
	 * @brief Detailed constructor
	 * @param input_maps reference to all FeatureMaps
	 * @param map_index index of the Feature- or ConsensusMap
	 * @param feature_index index of the Feature in the map
	 */
	GridFeature(const std::vector<std::vector<BaseFeature> >& input_maps,Size map_index,Size feature_index);
	/**
	 *@brief Returns the feature
	 */
	BaseFeature getFeature() const;
	/**
	 * @brief Destructor
	 */
	virtual ~GridFeature();
	/**
	 * @brief Returns the map index
	 */
	Size getMapIndex();
	/**
	 * @brief Returns the feature index
	 */
	Size getFeatureIndex();
	/**
	 * @brief Returns the ID of the GridFeature
	 */
	Int getID();
};
}

#endif /* GRIDFEATURE_H_ */
