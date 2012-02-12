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


#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{

	QTClusterFinder::QTClusterFinder(): 
		BaseGroupFinder(), feature_distance_(FeatureDistance())
	{
		setName(getProductName());

		defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
		defaults_.setValidStrings("use_identifications", StringList::create("true,false"));

		defaults_.insert("", feature_distance_.getDefaults());

		defaultsToParam_();
	}


	void QTClusterFinder::setParameters_(DoubleReal max_intensity)
	{
		use_IDs_ = String(param_.getValue("use_identifications")) == "true";
		max_diff_rt_ = param_.getValue("distance_RT:max_difference");
		max_diff_mz_ = param_.getValue("distance_MZ:max_difference");
		Param distance_params = param_.copy("");
		distance_params.remove("use_identifications");
		feature_distance_ = FeatureDistance(max_intensity, true);
		feature_distance_.setParameters(distance_params);
	}


	template <typename MapType>
	void QTClusterFinder::run_(const vector<MapType>& input_maps, 
														 ConsensusMap& result_map)
	{
		num_maps_ = input_maps.size();
		if (num_maps_ < 2)
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
																			 "At least two input maps required");
		}

		// set up the distance functor (and set other parameters):
		DoubleReal max_intensity = input_maps[0].getMaxInt();
		for (Size map_index = 1; map_index < num_maps_; ++map_index)
		{
			max_intensity = max(max_intensity, input_maps[map_index].getMaxInt());
		}
		setParameters_(max_intensity);

		// create the hash grid and fill it with features:
		// cout << "Hashing..." << endl;
		list<GridFeature> grid_features;
		Grid grid(Grid::ClusterCenter(max_diff_rt_, max_diff_mz_));
		for (Size map_index = 0; map_index < num_maps_; ++map_index)
		{
			for (Size feature_index = 0; feature_index < input_maps[map_index].size();
					 ++feature_index)
			{
				grid_features.push_back(
					GridFeature(input_maps[map_index][feature_index], map_index, 
											feature_index));
        GridFeature& gfeature = grid_features.back();
				// sort peptide hits once now, instead of multiple times later:
				BaseFeature& feature = const_cast<BaseFeature&>(
					grid_features.back().getFeature());
				for (vector<PeptideIdentification>::iterator pep_it = 
							 feature.getPeptideIdentifications().begin(); pep_it != 
							 feature.getPeptideIdentifications().end(); ++pep_it)
				{
					pep_it->sort();
				}
				grid.insert(std::make_pair(Grid::ClusterCenter(gfeature.getRT(), gfeature.getMZ()), &gfeature));
			}
		}

		// compute QT clustering:
		// cout << "Clustering..." << endl;
		list<QTCluster> clustering;
		computeClustering_(grid, clustering);
		// number of clusters == number of data points:
		Size size = clustering.size();

		ProgressLogger logger;
		logger.setLogType(ProgressLogger::CMD);
		logger.startProgress(0, size, "linking features");

		result_map.clear(false);
		while (!clustering.empty())
		{
			// cout << "Clusters: " << clustering.size() << endl;
			ConsensusFeature consensus_feature;
			makeConsensusFeature_(clustering, consensus_feature);
			result_map.push_back(consensus_feature);
			logger.setProgress(size - clustering.size());
		}

		logger.endProgress();
	}


	void QTClusterFinder::makeConsensusFeature_(list<QTCluster>& clustering,
																							ConsensusFeature& feature)
	{
		// find the best cluster:
		clustering.sort();
		list<QTCluster>::reverse_iterator best = clustering.rbegin();
		map<Size, GridFeature*> elements;
		best->getElements(elements);
		// cout << "Elements: " << elements.size() << endl;
	
		// create consensus feature:
		feature.setQuality(best->getQuality());
		for (map<Size, GridFeature*>::const_iterator it = elements.begin(); 
				 it != elements.end(); ++it)
		{
			feature.insert(it->first, it->second->getFeature());
		}
		feature.computeConsensus();
	
		// update the clustering:
		clustering.pop_back();
		for (list<QTCluster>::iterator it = clustering.begin(); 
				 it != clustering.end(); )
		{
			if (!it->update(elements)) // cluster is invalid (center point removed):
			{
				it = clustering.erase(it);
			}
			else ++it;
		}
	}


	void QTClusterFinder::run(const vector<ConsensusMap>& input_maps, 
														ConsensusMap& result_map)
	{
		run_(input_maps, result_map);
	}


	void QTClusterFinder::run(const std::vector<FeatureMap<> >& input_maps, 
														ConsensusMap& result_map)
	{
		run_(input_maps, result_map);
	}


	void QTClusterFinder::computeClustering_(Grid& grid, 
																					 list<QTCluster>& clustering)
	{
		clustering.clear();
		distances_.clear();
		// FeatureDistance produces normalized distances (between 0 and 1):
		const DoubleReal max_distance = 1.0;

		// iterate over all grid cells:
		for (Grid::iterator it = grid.begin(); it != grid.end(); ++it)
		{
      const Grid::CellIndex act_coords = it.index();
      const Int x = act_coords[0], y = act_coords[1];

      GridFeature* center_feature = it->second;
      QTCluster cluster(center_feature, num_maps_, max_distance, use_IDs_);

      // iterate over neighboring grid cells (1st dimension):
      for (int i = x - 1; i <= x + 1; ++i)
      {
        // iterate over neighboring grid cells (2nd dimension):
        for (int j = y - 1; j <= y + 1; ++j)
        {
          try
          { 
            const Grid::CellContent act_pos = grid.grid_at(Grid::CellIndex(i, j));

            for (Grid::const_cell_iterator it_cell = act_pos.begin(); it_cell != act_pos.end(); ++it_cell)
						{
							GridFeature* neighbor_feature = it_cell->second;
							// consider only "real" neighbors, not the element itself:
							if (center_feature != neighbor_feature)
							{
								DoubleReal dist = getDistance_(center_feature, 
																							 neighbor_feature);
								if (dist == FeatureDistance::infinity) 
								{
									continue; // conditions not satisfied
								}
								// if neighbor point is a possible cluster point, add it:
								if (!use_IDs_ || compatibleIDs_(cluster, neighbor_feature))
								{
									cluster.add(neighbor_feature, dist);
								}
							}
						}
					}
          catch (std::out_of_range &)
          { }
				}
			}
      clustering.push_back(cluster);
		}
	}


	DoubleReal QTClusterFinder::getDistance_(GridFeature* left,
																					 GridFeature* right)
	{
		// look-up in the distance map:
		pair<GridFeature*, GridFeature*> key = make_pair(min(left, right), 
																										 max(left, right));
		PairDistances::const_iterator pos = distances_.find(key);
		if (pos != distances_.end()) // distance computed before
		{
			return pos->second;
		}
		else // compute distance now and store it for later
		{
			DoubleReal dist = feature_distance_(left->getFeature(), 
																					right->getFeature()).second;
			distances_[key] = dist;
			return dist;
		}
	}


	bool QTClusterFinder::compatibleIDs_(QTCluster& cluster, 
																			 const GridFeature* neighbor)
	{
		if (cluster.getAnnotations().empty()) return true;
		if (neighbor->getAnnotations().empty()) return true;
		return (cluster.getAnnotations() == neighbor->getAnnotations());
	}


	QTClusterFinder::~QTClusterFinder()
	{
	}

} // namespace OpenMS
