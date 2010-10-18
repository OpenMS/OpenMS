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
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS;
using namespace std;


QTClusterFinder::QTClusterFinder(): BaseGroupFinder()
{
	BaseGroupFinder::setName(getProductName());

	defaults_.setValue("keep_subelements", "false", "For consensusXML input: Keep the sub-features of the input in the output (by default the consensus centroids are used)");
	defaults_.setValidStrings("keep_subelements", StringList::create("true,false"));

	defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
	defaults_.setValidStrings("use_identifications", StringList::create("true,false"));

	defaults_.setValue("diff_exponent:RT", 1.0, "RT differences are raised to this power", StringList::create("advanced"));
	defaults_.setMinFloat("diff_exponent:RT", 0.0);

	defaults_.setValue("diff_exponent:MZ", 2.0, "m/z differences are raised to this power", StringList::create("advanced"));
	defaults_.setMinFloat("diff_exponent:MZ", 0.0);

	defaults_.setSectionDescription("diff_exponent", "Absolute position differences are raised to this power (e.g. 1 for 'linear' distance, 2 for 'quadratic' distance).");

	defaults_.setValue("max_distance:RT", 100.0, "Maximum allowed distance in RT for a pair");
	defaults_.setMinFloat("max_distance:RT", 0.0);

	defaults_.setValue("max_distance:MZ", 0.3, "Maximum allowed distance in MZ for a pair");
	defaults_.setMinFloat("max_distance:MZ", 0.0);

	BaseGroupFinder::defaultsToParam_();
}


void QTClusterFinder::setParameters_()
{
	diff_exp_rt_ = (DoubleReal) param_.getValue("diff_exponent:RT");
	diff_exp_mz_ = (DoubleReal) param_.getValue("diff_exponent:MZ");
	max_dist_mz_ = (DoubleReal) param_.getValue("max_distance:MZ");
	max_dist_rt_ = (DoubleReal) param_.getValue("max_distance:RT");
	use_IDs_ = String(param_.getValue("use_identifications")) == "true";
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
	setParameters_();

	// create the hash grid and fill it with features:
	list<GridFeature> grid_features;
	HashGrid grid(max_dist_rt_, max_dist_mz_);
	for (Size map_index = 0; map_index < num_maps_; ++map_index)
	{
		for (Size feature_index = 0; feature_index < input_maps[map_index].size(); 
				 ++feature_index)
		{
			grid_features.push_back(GridFeature(input_maps[map_index][feature_index],
																					map_index, feature_index));
			// sort peptide hits once now, instead of multiple times later:
			BaseFeature& feature = const_cast<BaseFeature&>(
				grid_features.back().getFeature());
			for (vector<PeptideIdentification>::iterator pep_it = 
						 feature.getPeptideIdentifications().begin(); pep_it != 
						 feature.getPeptideIdentifications().end(); ++pep_it)
			{
				pep_it->sort();
			}
			grid.insert(&(grid_features.back()));
		}
	}

	// compute QT clustering:
	list<QTCluster> clustering;
	computeClustering_(grid, clustering);
	Size size = clustering.size(); // number of clusters == number of data points

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


void QTClusterFinder::computeClustering_(HashGrid& grid, 
																				 list<QTCluster>& clustering)
{
	clustering.clear();
	distances_.clear();
	DoubleReal max_distance = getDistance_(max_dist_rt_, max_dist_mz_);

	// iterate over all grid cells:
	for (GridCells::iterator it = grid.begin(); it != grid.end(); ++it)
	{
		pair<int, int> coords = it->first;
		list<GridElement*>& elements = it->second;
		int x = coords.first, y = coords.second;
		// iterate over all elements in the grid cell:
		for (list<GridElement*>::iterator center_element = elements.begin();
				 center_element != elements.end(); ++center_element)
		{
			// create a cluster for every element (this element is cluster center):
			GridFeature* center_feature = 
				dynamic_cast<GridFeature*>(*center_element);
			QTCluster cluster(center_feature, num_maps_, max_distance, use_IDs_);

			// iterate over neighboring grid cells (1st dimension):
			for (int i = x - 1; i <= x + 1; ++i)
			{
				// check for border cells:
				if ((i < 0) || (i > grid.getGridSizeX())) continue;

				// iterate over neighboring grid cells (2nd dimension):
				for (int j = y - 1; j <= y + 1; ++j)
				{
					// check for border cells:
					if ((j < 0) || (j > grid.getGridSizeY())) continue;

					// find all neighbors of the cluster:
					GridCells::iterator pos = grid.find(make_pair(i, j));
					if (pos == grid.end()) continue;
					list<GridElement*>& neighbor_elements = pos->second;

					// iterate over the neighbors:
					for (list<GridElement*>::iterator neighbor_element = 
								 neighbor_elements.begin(); neighbor_element != 
								 neighbor_elements.end(); ++neighbor_element)
					{
						GridFeature* neighbor_feature = 
							dynamic_cast<GridFeature*>(*neighbor_element);
						// consider only "real" neighbors, not the element itself:
						if (center_feature != neighbor_feature)
						{
							DoubleReal dist = getDistance_(center_feature, neighbor_feature);
							if (dist < 0) continue; // charge states don't match
							DoubleReal dist_mz = abs(neighbor_feature->mz - 
																			 center_feature->mz);
							DoubleReal dist_rt = abs(neighbor_feature->rt - 
																			 center_feature->rt);
							// if neighbor point is a possible cluster point, add it:
							if ((dist_mz <= max_dist_mz_) && (dist_rt <= max_dist_rt_) && 
									(!use_IDs_ || compatibleIDs_(cluster, neighbor_feature)))
							{
								cluster.add(neighbor_feature, dist);
							}
						}
					}
				}
			}
			clustering.push_back(cluster);
		}
	}
}


DoubleReal QTClusterFinder::getDistance_(GridFeature* left, GridFeature* right)
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
		DoubleReal dist = getDistance_(left->getFeature(), right->getFeature());
		distances_[key] = dist;
		return dist;
	}
}


DoubleReal QTClusterFinder::getDistance_(const BaseFeature& left, 
																				 const BaseFeature& right) const
{
	if (left.getCharge() != right.getCharge()) return -1; // charges don't match
	DPosition<2> pos_diff = left.getPosition() - right.getPosition();
	return getDistance_(pos_diff[RT], pos_diff[MZ]);
}


DoubleReal QTClusterFinder::getDistance_(DoubleReal pos_diff_rt, 
																				 DoubleReal pos_diff_mz) const
{
  // RT:
	pos_diff_rt = abs(pos_diff_rt) / max_dist_rt_;
  pos_diff_rt = pow(pos_diff_rt, diff_exp_rt_);

  // MZ:
	pos_diff_mz = abs(pos_diff_mz) / max_dist_mz_;
  pos_diff_mz = pow(pos_diff_mz, diff_exp_mz_);

  DoubleReal result = pos_diff_rt + pos_diff_mz;
  return result;
}


bool QTClusterFinder::compatibleIDs_(const QTCluster& cluster, 
																		 const GridFeature* neighbor) const
{
	const vector<PeptideIdentification>& peptides = 
		neighbor->getFeature().getPeptideIdentifications();
	// a neighbor feature without identifications always matches:
	if (peptides.empty()) return true;
	set<AASequence> sequences;
	for (vector<PeptideIdentification>::const_iterator pep_it = peptides.begin(); 
			 pep_it != peptides.end(); ++pep_it)
	{
		if (pep_it->getHits().empty()) continue; // shouldn't be the case
		sequences.insert(pep_it->getHits()[0].getSequence());
	}
	return (cluster.getAnnotations() == sequences);
}


QTClusterFinder::~QTClusterFinder()
{
}
