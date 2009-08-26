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
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on pose clustering.

		Pose clustering analyses pair distances to find the most probable transformation of retention times.

		For further details see:
		@n Eva Lange et.al
		@n A Geometric Approach for the Alignment of Liquid Chromatography-Mass Spectrometry Data
		@n ISMB/ECCB 2007

	  @htmlinclude OpenMS_MapAlignmentAlgorithmPoseClustering.parameters

		@ingroup MapAlignment

		@todo write test, work out the TODOs (Clemens)
	*/
	class OPENMS_DLLAPI MapAlignmentAlgorithmPoseClustering
		: public MapAlignmentAlgorithm
	{
	 public:
		/// Default constructor
		MapAlignmentAlgorithmPoseClustering();

		/// Destructor
		virtual ~MapAlignmentAlgorithmPoseClustering();

		// Docu in base class
		virtual void alignPeakMaps(std::vector< MSExperiment<> >&, std::vector<TransformationDescription>&);

		// Docu in base class
		virtual void alignFeatureMaps(std::vector< FeatureMap<> >&, std::vector<TransformationDescription>&);

		/// Creates a new instance of this class (for Factory)
		static MapAlignmentAlgorithm* create()
		{
			return new MapAlignmentAlgorithmPoseClustering();
		}

		/// Returns the product name (for the Factory)
		static String getProductName()
		{
			return "pose_clustering_affine";
		}

	 protected:

		/**
			@brief This will compute a linear regression based on all consensus
			features which contain feature handles from the x and the y map.  If there
			is only one pair, assume slope is one and set intercept accordingly.  If
			there is no pair at all, throw an exception.

			@throw Exception::Precondition if no consensus feature contains feature
			handles from both maps
		*/
		TransformationDescription calculateRegression_(Size const index_x_map, Size const index_y_map, ConsensusMap const& consensus_map, bool symmetric_regression) const;

	 private:

		/// Copy constructor intentionally not implemented -> private
		MapAlignmentAlgorithmPoseClustering(const MapAlignmentAlgorithmPoseClustering& );
		///Assignment operator intentionally not implemented -> private
		MapAlignmentAlgorithmPoseClustering& operator=(const MapAlignmentAlgorithmPoseClustering& );

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMPOSECLUSTERING_H
