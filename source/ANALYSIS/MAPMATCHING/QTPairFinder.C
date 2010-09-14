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


#include <OpenMS/ANALYSIS/MAPMATCHING/QTPairFinder.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS;
using namespace std;

QTPairFinder::QTPairFinder() : BaseGroupFinder()
{
	BaseGroupFinder::setName(getProductName());

	defaults_.setValue("keep_subelements", "true",
		      "Use subelements of consensusXML for output");
	defaults_.setValidStrings("keep_subelements",StringList::create("true,false"));

	defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
			defaults_.setValidStrings("use_identifications", StringList::create("true,false"));

	defaults_.setValue("diff_exponent:RT", 1.0,
	      "RT differences are raised to this power", StringList::create("advanced"));
	    defaults_.setMinFloat("diff_exponent:RT",0.1);

	defaults_.setValue("diff_exponent:MZ", 2.0,
			"MZ differences are raised to this power", StringList::create("advanced"));
	defaults_.setMinFloat("diff_exponent:MZ",0.1);

	defaults_.setSectionDescription("diff_exponent",
			"Absolute position differences are raised to this power. "
			"E.g. 1 for 'linear' distance, 2 for 'quadratic' distance");

	defaults_.setValue("max_pair_distance:RT", 100.0,"Maximal allowed distance in RT for a pair");
	defaults_.setMinFloat("max_pair_distance:RT",0.);

	defaults_.setValue("max_pair_distance:MZ", 0.3, "Maximal allowed distance in MZ for a pair [Unit defined by 'mz_unit']");
	defaults_.setMinFloat("max_pair_distance:MZ",0.);

	BaseGroupFinder::defaultsToParam_();
}

void QTPairFinder::run(const std::vector<ConsensusMap>& input_maps,ConsensusMap &result_map)
{
	if ( input_maps.size() < 2 )
	{
		throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"more than on input maps required");
	}
	diff_exponent_rt = (DoubleReal) param_.getValue("diff_exponent:RT");
	diff_exponent_mz = (DoubleReal) param_.getValue("diff_exponent:MZ");
	max_pair_distance_mz = (DoubleReal) param_.getValue("max_pair_distance:MZ");
	max_pair_distance_rt = (DoubleReal) param_.getValue("max_pair_distance:RT");
	use_IDs = String(param_.getValue("use_identifications")) == "true";
	keep_subelements = String(param_.getValue("keep_subelements")) == "true";

	//Convert Features into BaseFeatures
	vector<vector<BaseFeature> > input(input_maps.size());
	for (Size i=0;i<input_maps.size();++i)
	{
		for (ConsensusMap::const_iterator feature_it=input_maps[i].begin();feature_it!=input_maps[i].end();++feature_it)
		{
			input[i].push_back(*feature_it);
		}
	}
	maps_size=input_maps.size();

	//Create the HashGrid an fill it with the BaseFeatures
	HashGrid grid(max_pair_distance_rt,max_pair_distance_mz);
	for (Size map_index=0; map_index < input_maps.size(); ++map_index)
	{
		for (Size feature_index=0; feature_index < input_maps[map_index].size(); ++ feature_index)
		{
			grid_features.push_back(GridFeature(input, map_index, feature_index));
			grid.insert(&(grid_features.back()));
		}
	}
	calculateDistances(grid);
	ProgressLogger logger;
	logger.setLogType(ProgressLogger::CMD);
	logger.startProgress(0,grid.getNumberOfElements(),"linking features");
	result_map.clear(false);
	//Check if subelements should be kept and run the specific method
	if (keep_subelements)
		makeConsensus(grid,logger,result_map,input_maps);
	else
		makeConsensus(grid,logger,result_map);
	logger.endProgress();
}

void QTPairFinder::run(const std::vector<FeatureMap<> >& input_maps, ConsensusMap& result_map)
{
	if ( input_maps.size() < 2 )
	{
		throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"more than on input maps required");
	}
	diff_exponent_rt = (DoubleReal) param_.getValue("diff_exponent:RT");
	diff_exponent_mz = (DoubleReal) param_.getValue("diff_exponent:MZ");
	max_pair_distance_mz = (DoubleReal) param_.getValue("max_pair_distance:MZ");
	max_pair_distance_rt = (DoubleReal) param_.getValue("max_pair_distance:RT");
	use_IDs = String(param_.getValue("use_identifications")) == "true";
	keep_subelements = String(param_.getValue("keep_subelements")) == "true";

	//Convert Features into BaseFeatures
	vector<vector<BaseFeature> > input(input_maps.size());
	for (Size i=0;i<input_maps.size();++i)
	{
		for (FeatureMap<>::const_iterator feature_it=input_maps[i].begin();feature_it!=input_maps[i].end();++feature_it)
		{
			input[i].push_back(*feature_it);
		}
	}
	maps_size=input_maps.size();

	//Create the HashGrid an fill it with the BaseFeatures
	HashGrid grid(max_pair_distance_rt,max_pair_distance_mz);
	for (Size map_index=0; map_index < input_maps.size(); ++map_index)
	{
		for (Size feature_index=0; feature_index < input_maps[map_index].size(); ++ feature_index)
		{
			grid_features.push_back(GridFeature(input, map_index, feature_index));
			grid.insert(&(grid_features.back()));
		}
	}
	calculateDistances(grid);
	ProgressLogger logger;
	logger.setLogType(ProgressLogger::CMD);
	logger.startProgress(0,grid.getNumberOfElements(),"linking features");
	result_map.clear(false);
	makeConsensus(grid,logger,result_map);
	logger.endProgress();
}

void QTPairFinder::calculateDistances(HashGrid& grid)
{
	//Iterate over all neighbor cells of every data point and calculate all possible distances
	for (GridCells::iterator it=grid.begin();it!=grid.end();++it)
	{
		pair<int,int> act_coords=it->first;
		list<GridElement*>& elements=it->second;
		int x=act_coords.first;
		int y=act_coords.second;
		for (list<GridElement*>::iterator center_element=elements.begin();center_element!=elements.end();++center_element)
		{
			GridFeature* center_feature_ptr = dynamic_cast<GridFeature*> (*center_element);
			const BaseFeature& center_feature=center_feature_ptr->getFeature();
			for (int i=x-1;i<=x+1;++i)
			{
				//Check the border cells
				if (i<0 || i>grid.getGridSizeX())
					continue;
				for (int j=y-1;j<=y+1;++j)
				{
					if (j<0 || j>grid.getGridSizeY())
						continue;
					GridCells::iterator act_pos=grid.find(make_pair(i,j));
					if (act_pos==grid.end())
						continue;

					list<GridElement*>& neighbor_elements=act_pos->second;
					for (list<GridElement*>::iterator neighbor_element=neighbor_elements.begin();neighbor_element!=neighbor_elements.end();++neighbor_element)
					{
						GridFeature* neighbor_feature_ptr = dynamic_cast<GridFeature*> (*neighbor_element);
						if (center_feature_ptr!=neighbor_feature_ptr)
						{
							const BaseFeature& neighbor_feature=neighbor_feature_ptr->getFeature();
							distances[make_pair(min(center_feature_ptr,neighbor_feature_ptr),max(center_feature_ptr,neighbor_feature_ptr))]=distance_(center_feature,neighbor_feature);
						}
					}
				}
			}
		}
	}
}

QTCluster QTPairFinder::QTClust(HashGrid& act_grid)
{
	DoubleReal max_distance=distance_(max_pair_distance_rt,max_pair_distance_mz);
	if (act_grid.getNumberOfElements() <= 1)
	{
			GridFeature* element_ptr = dynamic_cast<GridFeature*> (*(act_grid.begin()->second.begin()));
			return QTCluster(element_ptr,maps_size,max_distance);
	}
	//Create a multiset in which the clusters are ordered depending on their quality (lowest first)
	multiset<QTCluster> cluster_set;

	//Iterate over all neighbor cells of every data point and check if they belong together in one cluster
	for (GridCells::iterator it=act_grid.begin();it!=act_grid.end();++it)
	{
		pair<int,int> act_coords=it->first;
		list<GridElement*>& elements=it->second;
		int x=act_coords.first;
		int y=act_coords.second;
		for (list<GridElement*>::iterator center_element=elements.begin();center_element!=elements.end();++center_element)
		{
			GridFeature* center_feature_ptr = dynamic_cast<GridFeature*> (*center_element);
			//Create a cluster for every element (this element is cluster center)
			QTCluster act_cluster(center_feature_ptr,maps_size,max_distance);
			for (int i=x-1;i<=x+1;++i)
			{
				//Check the border cells
				if (i<0 || i>act_grid.getGridSizeX())
					continue;
				for (int j=y-1;j<=y+1;++j)
				{
					if (j<0 || j>act_grid.getGridSizeY())
						continue;

					//Find all neighbors of the cluster
					GridCells::iterator act_pos=act_grid.find(make_pair(i,j));
					if (act_pos==act_grid.end())
						continue;
					list<GridElement*>& neighbor_elements=act_pos->second;

					//Iterate over the neighbors
					for (list<GridElement*>::iterator neighbor_element=neighbor_elements.begin();neighbor_element!=neighbor_elements.end();++neighbor_element)
					{
						GridFeature* neighbor_feature_ptr = dynamic_cast<GridFeature*> (*neighbor_element);
						//Consider only "real" neighbors, not the element itself
						if (center_feature_ptr!=neighbor_feature_ptr)
						{
							DoubleReal distance=distance_(center_feature_ptr,neighbor_feature_ptr);
							DoubleReal mz_distance=abs(neighbor_feature_ptr->mz - center_feature_ptr->mz);
							DoubleReal rt_distance=abs(neighbor_feature_ptr->rt - center_feature_ptr->rt);
							//Check if the neighbor point is a possible cluster point
							if (distance>0 && mz_distance <= max_pair_distance_mz && rt_distance <= max_pair_distance_rt && (!use_IDs || compatibleIDs_(neighbor_feature_ptr->getFeature(), center_feature_ptr->getFeature())))
								act_cluster.add(neighbor_feature_ptr,distance);
						}
					}
				}
			}
			cluster_set.insert(act_cluster);
		}
	}
	//Return the cluster with the highest quality
	return *cluster_set.rbegin();
}

void QTPairFinder::makeConsensus(HashGrid& grid,ProgressLogger& logger,ConsensusMap& result_map)
{
	Int grid_size=grid.getNumberOfElements();
	while(grid.getNumberOfElements() > 0)
	{
		QTCluster act_cluster=QTClust(grid);
		const map<Size,GridFeature*>& cluster_members=act_cluster.getClusterMembers();

		ConsensusFeature consensus_feature;
		consensus_feature.setQuality(act_cluster.getQuality());
		for (map<Size,GridFeature*>::const_iterator it=cluster_members.begin();it!=cluster_members.end();++it)
		{
			BaseFeature act_feature=it->second->getFeature();
			consensus_feature.insert(it->first,act_feature);
			grid.removeElement(it->second);
		}
		consensus_feature.computeConsensus();
		result_map.push_back(consensus_feature);
		logger.setProgress(grid_size - grid.getNumberOfElements());
	}
}

void QTPairFinder::makeConsensus(HashGrid& grid,ProgressLogger& logger,ConsensusMap& result_map,const std::vector<ConsensusMap>& input_maps)
{
	Int grid_size=grid.getNumberOfElements();
	while(grid.getNumberOfElements() > 0)
	{
		QTCluster act_cluster=QTClust(grid);
		const map<Size,GridFeature*>& cluster_members=act_cluster.getClusterMembers();

		ConsensusFeature consensus_feature;
		consensus_feature.setQuality(act_cluster.getQuality());
		for (map<Size,GridFeature*>::const_iterator it=cluster_members.begin();it!=cluster_members.end();++it)
		{
			//Calculate the correct map index
			Size file_description_offset=0;
			for (std::vector<ConsensusMap>::const_iterator map_it=input_maps.begin();map_it!=input_maps.end()-1 && *map_it != input_maps[it->second->getMapIndex()];++map_it)
			{
				file_description_offset+=map_it->getFileDescriptions().size();
			}
			ConsensusFeature act_consensus_feature=input_maps[it->second->getMapIndex()][it->second->getFeatureIndex()];
			const ConsensusFeature::HandleSetType& act_features=act_consensus_feature.getFeatures();
			for (ConsensusFeature::HandleSetType::const_iterator feature_it=act_features.begin();feature_it!=act_features.end();++feature_it)
			{
				FeatureHandle act_feature=*feature_it;
				//Assign the correct map index
				act_feature.setMapIndex(file_description_offset+act_feature.getMapIndex());
				consensus_feature.insert(act_feature);
			}
			grid.removeElement(it->second);
		}
		consensus_feature.computeConsensus();
		result_map.push_back(consensus_feature);
		logger.setProgress(grid_size - grid.getNumberOfElements());
	}
}


DoubleReal QTPairFinder::distance_(GridFeature* left,GridFeature* right) const
{
	//Look up the distance map
	map<pair<GridFeature*,GridFeature*>,DoubleReal>::const_iterator distance_position=distances.find(make_pair(min(left,right),max(left,right)));
	if (distance_position!=distances.end())
		return distance_position->second;
	else
		return -1;
}

DoubleReal QTPairFinder::distance_(BaseFeature const & left,BaseFeature const & right ) const
{
	  // distance from charge
	  if ( left.getCharge() != right.getCharge() )
	  {
	    return -1;
	  }
	  DPosition<2> position_difference = left.getPosition() - right.getPosition();
	  return distance_(position_difference[RT],position_difference[MZ]);
}

DoubleReal QTPairFinder::distance_(DoubleReal position_difference_rt,DoubleReal position_difference_mz) const
{
  // .. in RT
  if ( position_difference_rt < 0 )
  {
	  position_difference_rt = -position_difference_rt;
  }
  position_difference_rt *= 1/max_pair_distance_rt;
  position_difference_rt = pow(position_difference_rt,diff_exponent_rt);

  // .. in MZ
  if ( position_difference_mz < 0 )
  {
	  position_difference_mz = -position_difference_mz;
  }

  DoubleReal result = position_difference_rt + position_difference_mz;

  return result;
}

bool QTPairFinder::compatibleIDs_(const BaseFeature& feat1, const BaseFeature& feat2) const
{
	vector<PeptideIdentification> pep1 = feat1.getPeptideIdentifications(),pep2 = feat2.getPeptideIdentifications();
	// a feature without identifications always matches:
	if (pep1.empty() || pep2.empty()) return true;
	set<AASequence> best1, best2;
	for (vector<PeptideIdentification>::iterator pep_it = pep1.begin(); pep_it != pep1.end(); ++pep_it)
	{
		if (pep_it->getHits().empty()) continue; // shouldn't be the case
		pep_it->sort();
		best1.insert(pep_it->getHits()[0].getSequence());
	}
	for (vector<PeptideIdentification>::iterator pep_it = pep2.begin(); pep_it != pep2.end(); ++pep_it)
	{
		if (pep_it->getHits().empty()) continue; // shouldn't be the case
		pep_it->sort();
		best2.insert(pep_it->getHits()[0].getSequence());
	}
	return (best1 == best2);
}


QTPairFinder::~QTPairFinder()
{
	// TODO Auto-generated destructor stub
}
