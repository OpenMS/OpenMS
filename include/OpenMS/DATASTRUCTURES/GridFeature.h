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

#ifndef OPENMS_DATASTRUCTURES_GRIDFEATURE_H
#define OPENMS_DATASTRUCTURES_GRIDFEATURE_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>

namespace OpenMS
{

/**
 * @brief Representation of a feature in a hash grid.
 * 
 * A GridFeature can be stored in a HashGrid and points to a BaseFeature (Feature or ConsensusFeature). Used for QT feature grouping (see QTClusterFinder).
 */

	class OPENMS_DLLAPI GridFeature 
	{
		private:
		/// Reference to the contained feature
		const BaseFeature& feature_;

		/// Index of the feature map or consensus map
		Size map_index_;

		/// Index of the feature in the map
		Size feature_index_;

		/// Set of peptide sequences annotated to the feature
		std::set<AASequence> annotations_;

	public:
		/**
		 * @brief Detailed constructor
		 * @param feature Reference to the contained feature
		 * @param map_index Index of the feature map or consensus map
		 * @param feature_index Index of the feature in the map
		 */
		GridFeature(const BaseFeature& feature, Size map_index, Size feature_index);

		/// Returns the feature
		const BaseFeature& getFeature() const;

		/// Destructor
		virtual ~GridFeature();

		/// Returns the map index
		Size getMapIndex() const;
		
		/// Returns the feature index
		Size getFeatureIndex() const;

		/// Returns the ID of the GridFeature (same as the feature index)
		Int getID() const;

		/// Returns the set of peptide sequences annotated to the cluster center
		const std::set<AASequence>& getAnnotations() const;

		/// Returns the feature RT
    DoubleReal getRT() const;

		/// Returns the feature m/z
    DoubleReal getMZ() const;
	};
}

#endif // OPENMS_DATASTRUCTURES_GRIDFEATURE_H
