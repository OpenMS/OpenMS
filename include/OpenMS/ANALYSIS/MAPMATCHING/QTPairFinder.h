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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_QTPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_QTPAIRFINDER_H

#include <OpenMS/DATASTRUCTURES/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS {

/**
   @brief This class implements an element pair finding algorithm.

   It offers a method to determine element pairs across element maps.

   Motivation and definition of the distance measure:

   The similarity value should express our confidence that one element might
   possibly be matched to the other.  Larger quality values are better, the
   maximal similarity is one.  Let \f$\Delta_\textit{RT}\f$ and
   \f$\Delta_\textit{MZ}\f$ be the absolute values of the RT and MZ differences
   in the data.
   Then the distance value is
   \f[
   \big(
   \frac
     {\Delta_\textit{RT}}
     {\big( \textit{max\_pair\_distance:RT} \big)}
   ^\textit{diff\_exponent\_RT}
   +
   \frac
     {\Delta_\textit{MZ}}
     {\big( \textit{max\_pair\_distance:MZ} \big)}
   ^\textit{diff\_exponent\_MZ}
   \big)


   \f]

   Elements of different charge states are not linked.

   Choosing <i>diff_exponent</i>: This parameter controls the growth rate of the
   penalty for differences.  It is for example possible to search for pairs
   using the absolute distance in RT (which should not be very susceptible to
   outliers) and the squared distance in MZ (where small difference occur
   frequently, but large differences indicate a mismatch).

   Stability criterion:  The linking depends on the quality method of the the QT clusters.

   @htmlinclude OpenMS_QTPairFinder.parameters

   @ingroup FeatureGrouping
   */

class OPENMS_DLLAPI QTPairFinder  : public BaseGroupFinder{
private:

	/**
	 * maximal distance in RT direction.
	 */
	DoubleReal max_pair_distance_rt;
	/**
	 * maximal distance in m/z direction.
	 */
	DoubleReal max_pair_distance_mz;
	/**
	 * Distances wrt RT is raised to this power.
	 */
	DoubleReal diff_exponent_rt;
	/**
	 * Distances wrt mz is raised to this power.
	 */
	DoubleReal diff_exponent_mz;
	/**
	 * @brief Distance map.
	 * To store every distance only once, distance between to features are accessible by searching for a pair where the first position is the smaller pointer value.
	 */
	std::map<std::pair<GridFeature*,GridFeature*>,DoubleReal> distances;
	/**
	 * @brief calculates the distance between two grid features by looking it up in the distance map.
	 * Returns -1 if the distance is not contained in the distance map
	 */
	DoubleReal distance_(GridFeature* left,GridFeature* right) const;
	/**
	 * @brief calculates the distance between two features.
	 * Returns -1 if the two features have different charge.
	 */
	DoubleReal distance_(BaseFeature const & left,BaseFeature const & right ) const;
	/**
	 * @brief calculates the distance of RT difference and m/z difference.
	 * Uses m/z and RT exponent values.
	 */
	DoubleReal distance_(DoubleReal position_difference_rt,DoubleReal position_difference_mz) const;
	/**
	 * @brief number of input maps
	 */
	Size maps_size;
	/**
	 * @brief Checks if the feature identifiers are considered for linking
	 */
	bool use_IDs;
	/**
	 * @brief Checks the output mode of consensusMap input.
	 *  If set, the feature information of the consensusXMLs remains, while it gets lost if it is not set.
	 *  In this case, consensusFeatures are treated as features.
	 */
	bool keep_subelements;
	/**
	 * @brief Stores a pointer of the input map if the input is a consensusMap
	 */
	std::vector<ConsensusMap> *input_maps;
	/**
	 * @brief Produces the output.
	 * Is used if the keep_subelements flag is not set
	 */
	void makeConsensus(HashGrid& grid,ProgressLogger& logger,ConsensusMap& result_map);
	/**
	 * @brief Produces the output.
	 * Is used if the keep_subelements flag is set
	 */
	void makeConsensus(HashGrid& grid,ProgressLogger& logger,ConsensusMap& result_map,const std::vector<ConsensusMap>& input_maps);
	/**
	 * @brief List of grid features
	 */
	std::list<GridFeature> grid_features;
	/**
	 * @brief Calculates the distance in the hash grid
	 */
	void calculateDistances(HashGrid& grid);

protected:
	enum
	{
		RT = Peak2D::RT,
		MZ = Peak2D::MZ
	};

public:
	/**
	 * @brief Constructor
	 */
	QTPairFinder();
	/**
	 * @brief Destructor
	 */
	virtual ~QTPairFinder();
	/**
	 * @brief QT clustering method
	 */
	QTCluster QTClust(HashGrid& act_grid);

	/**
	 * @brief Returns the name of the product
	 */
	static const String
	getProductName()
	{
		return "qt";
	}

	/**
	       @brief Run the algorithm for FeatureMaps

	       @note Exactly two @em input maps must be provided.

	       @exception Exception::IllegalArgument is thrown if the input data is not valid.
	       */

	void run(const std::vector<FeatureMap<> >& input_maps, ConsensusMap& result_map);
	/**
	       @brief Run the algorithm for ConsensusMaps

	       @note Exactly two @em input maps must be provided.

	       @exception Exception::IllegalArgument is thrown if the input data is not valid.
	       */
	void run(const std::vector<ConsensusMap>& input_maps,ConsensusMap &result_map );


	/**
		@brief Checks if the peptide IDs of two features are compatible.

		A feature without identification is always compatible. Otherwise, two features are compatible if the best peptide hits of their identifications have the same sequence.
	 */
	bool compatibleIDs_(const BaseFeature& feat1, const BaseFeature& feat2) const;

	/**
	 * Returns an instance of this class
	 */
	static BaseGroupFinder*	create()
	{
		return new QTPairFinder();
	}
};
}
#endif /* QTCLUSTERING_H_ */
