// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>

//using namespace std;

namespace OpenMS
{
	ClusterAnalyzer::ClusterAnalyzer()
	{
	}

	ClusterAnalyzer::~ClusterAnalyzer()
	{
	}

	ClusterAnalyzer& ClusterAnalyzer::operator = (const ClusterAnalyzer& source)
	{
		//ALWAYS CHECK FOR SELF ASSIGNEMT!
		if (this == &source)
		{
			return *this;
		}
		//...
		return *this;
	}

	std::vector< Real > ClusterAnalyzer::averageSilhouetteWidth(std::vector<BinaryTreeNode>& tree, DistanceMatrix<Real>& original)
	{
		//throw exception if cannot be legal clustering
		if(tree.size() < 1)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "tree is empty but minimal clustering hirachy has at least one level");
		}

		std::vector< Real > average_silhouette_widths; //for each step from the average silhouette widths of the clusters
		std::vector< Real > interdist_i(original.dimensionsize(),std::numeric_limits<Real>::max()); //for each element i holds the min. average intercluster distance in cluster containing i
		std::vector< UInt > cluster_with_interdist(original.dimensionsize(),0); //for each element i holds which cluster originated the min. intercluster distance
		std::vector< Real > intradist_i(original.dimensionsize(),0.0); //for each element i holds the average intracluster distance in [i]

		//inital values for interdis_i and cluster_with_interdist
		for(UInt i = 1; i < original.dimensionsize(); ++i)
		{
			for(UInt j = 0; j < i; ++j)
			{
				if(original.getValue(i,j)<interdist_i[i])
				{
					interdist_i[i] = original.getValue(i,j);
					cluster_with_interdist[i] = j;
				}
				if(original.getValue(i,j)<interdist_i[j])
				{
					interdist_i[j] = original.getValue(i,j);
					cluster_with_interdist[j] = i;
				}
			}
		}

		/* to manually retrace
		for(UInt i = 0; i < original.dimensionsize(); ++i)
		{
			std::cout << interdist_i[i] << " | " << cluster_with_interdist[i] << " | " << intradist_i[i] << std::endl;
		}
		*/

		//initial cluster state
		std::vector< std::vector < UInt > > clusters(original.dimensionsize());
		for(UInt i = 0; i < clusters.size(); ++i)
		{
			clusters[i].push_back(i);
		}

		for(UInt t = 0; t < tree.size()-1; ++t) //last steps silhouettes would be all 0 respectively not defined
		{

			for(UInt i = 0; i < original.dimensionsize(); ++i)
			{
				std::vector<UInt>::iterator in_left = std::find(clusters[tree[t].left_child].begin(),clusters[tree[t].left_child].end(),i);
				std::vector<UInt>::iterator in_right = std::find(clusters[tree[t].right_child].begin(),clusters[tree[t].right_child].end(),i);

				if(in_left==clusters[tree[t].left_child].end() && in_right==clusters[tree[t].right_child].end()) //i (!element_of) left or right
				{
					//intradist_i is always kept
					//handle interdist:
					if(tree[t].left_child != cluster_with_interdist[i] && tree[t].right_child != cluster_with_interdist[i]) //s(i)_nr (!element_of) left or right
					{
						Real interdist_merged(0);
						for(UInt j = 0; j < clusters[tree[t].left_child].size(); ++j)
						{
							interdist_merged += original.getValue(i,clusters[tree[t].left_child][j]);
						}
						for(UInt j = 0; j < clusters[tree[t].right_child].size(); ++j)
						{
							interdist_merged += original.getValue(i,clusters[tree[t].right_child][j]);
						}
						interdist_merged /= (Real)(clusters[tree[t].left_child].size()+clusters[tree[t].right_child].size());
						if (interdist_merged < interdist_i[i])
						{
							interdist_i[i]=interdist_merged;
							cluster_with_interdist[i] = tree[t].left_child;
						}
					}
					else //s(i)_nr (element_of) left or right
					{
						//calculate interdist_i to merged
						UInt k; //the one cluster of the two merged which does NOT contain s(i)_nr
						if(tree[t].right_child!=cluster_with_interdist[i] )
						{
							k=tree[t].right_child;
						}
						else
						{
							k=tree[t].left_child;
						}
						Real interdist_merged(0);
						for(UInt j = 0; j < clusters[k].size(); ++j)
						{
							interdist_merged += original.getValue(i,clusters[k][j]);
						}
						interdist_merged += (clusters[cluster_with_interdist[i]].size()*interdist_i[i]);
						interdist_merged /= (Real)(clusters[k].size()+clusters[cluster_with_interdist[i]].size());
						//if new inderdist is smaller that old min. nothing else has to be done
						if (interdist_merged <= interdist_i[i])
						{
							interdist_i[i]=interdist_merged;
							cluster_with_interdist[i] = tree[t].left_child;
						}
						// else find min av. dist from other clusters to i
						else
						{
							interdist_i[i]=interdist_merged;
							cluster_with_interdist[i] = tree[t].left_child;

							for(UInt u = 0; u < clusters.size(); ++u)
							{
								if(u != tree[t].left_child && u != tree[t].right_child && clusters[u].end() == std::find(clusters[u].begin(),clusters[u].end(),i))
								{
									Real min_interdist_i(0);
									for(UInt v = 0; v < clusters[u].size(); ++v)
									{
										min_interdist_i += original.getValue(clusters[u][v],i);
									}
									min_interdist_i /= (Real)clusters[u].size();
									if (min_interdist_i < interdist_i[i])
									{
										interdist_i[i]=min_interdist_i;
										cluster_with_interdist[i] = u;
									}
								}
							}
						}
					}

				}
				else //i (element_of) left or right
				{
					UInt k,l; //k is the cluster that is one of the merged but not the one containing i, l the cluster containing i
					if(in_left==clusters[tree[t].left_child].end())
					{
						l = tree[t].right_child;
						k = tree[t].left_child;
					}
					else
					{
						l = tree[t].left_child;
						k = tree[t].right_child;
					}

					if(k != cluster_with_interdist[i]) //s(i)_nr (!element_of) left or right cluster
					{
						//interdist_i is kept
						//but intradist_i has to be updated
						intradist_i[i] *= clusters[l].size()-1;
						for(UInt j = 0; j < clusters[k].size(); ++j)
						{
							intradist_i[i] += original.getValue(i,clusters[k][j]);
						}
						intradist_i[i] /= (Real)(clusters[k].size()+(clusters[l].size()-1));
					}
					else //s(i)_nr (element_of) left or right
					{
						//intradist_i has to be updated
						intradist_i[i] *= clusters[l].size()-1;
						intradist_i[i] += (clusters[k].size()*interdist_i[i]);
						intradist_i[i] /= (Real)(clusters[k].size()+(clusters[l].size()-1));
						//find new min av. interdist_i
						interdist_i[i]=std::numeric_limits<Real>::max();
						for(UInt u = 0; u < clusters.size(); ++u)
						{
							if(u!=l && u!=k)
							{
								Real av_interdist_i(0);
								for(UInt v = 0; v < clusters[u].size(); ++v)
								{
									av_interdist_i += original.getValue(clusters[u][v],i);
								}
								av_interdist_i /= (Real)clusters[u].size();
								if (av_interdist_i < interdist_i[i])
								{
									interdist_i[i]=av_interdist_i;
									cluster_with_interdist[i] = u;
								}
							}
						}
					}
				}
			}
			//redo clustering following tree
			//pushback elements of right_child to left_child (and then erase second)
			clusters[tree[t].left_child].insert(clusters[tree[t].left_child].end(),clusters[tree[t].right_child].begin(),clusters[tree[t].right_child].end());

			//erase second one
			clusters.erase(clusters.begin()+tree[t].right_child);

			//adept the cluster indices in clusters with interdist
			for(UInt x = 0; x < cluster_with_interdist.size();++x)
			{
				if(cluster_with_interdist[x]>tree[t].right_child)
				{
					--cluster_with_interdist[x];
				}
			}

			/* to manually retrace
			for(UInt x = 0; x < clusters.size();++x)
			{
				for(UInt y = 0; y < clusters[x].size();++y)
				{
					std::cout << clusters[x][y] << " ";
				}
				std::cout << " | ";
			}
			std::cout << std::endl;
			std::cout << "---------" << std::endl;
			for(UInt z = 0; z < original.dimensionsize(); ++z)
			{
				std::cout << interdist_i[z] << " , " << intradist_i[z] << " , " << cluster_with_interdist[z] << " , ";
				std::cout << ((interdist_i[z] - intradist_i[z]) / std::max(interdist_i[z],intradist_i[z])) << std::endl;

			}
			std::cout << "---------" << std::endl;
			*/

			//calculate average silhouette width for clusters and then overall average silhouette width for cluster step
			Real average_overall_silhouette(0); // from cluster step
			std::vector<Real> silhouettes(original.dimensionsize(),0.0);
			for(UInt g = 0; g < clusters.size(); ++g)
			{
				if(clusters[g].size()>1)
				{
					//collect silhouettes clusterwise so that average cluster silhouettes will be easily accessible
					for(UInt h = 0; h < clusters[g].size(); ++h)
					{
						if(interdist_i[clusters[g][h]]!=0)
						{
							average_overall_silhouette += (interdist_i[clusters[g][h]] - intradist_i[clusters[g][h]]) / std::max(interdist_i[clusters[g][h]],intradist_i[clusters[g][h]]);
							/* to manually retrace
								silhouettes[clusters[g][h]] = (interdist_i[clusters[g][h]] - intradist_i[clusters[g][h]]) / std::max(interdist_i[clusters[g][h]],intradist_i[clusters[g][h]]);
							*/
						}
					}
				}
			}
			/* to manually retrace
				for(UInt i = 0; i < silhouettes.size(); ++i)
				{
					std::cout << "s(" <<  (i) << ") = " << silhouettes[i] << std::endl;
				}
				std::cout << "---------" << std::endl;
			*/
			average_silhouette_widths.push_back(average_overall_silhouette/(Real)original.dimensionsize());
		}
		average_silhouette_widths.push_back(0.0);
		return average_silhouette_widths;
	}

	std::vector< Real > ClusterAnalyzer::dunnIndices(std::vector<BinaryTreeNode>& tree, DistanceMatrix<Real>& original, bool tree_from_singlelinkage)
	{
		//throw exception if cannot be legal clustering
		if(tree.size() < 1)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "tree is empty but minimal clustering hirachy has at least one level");
		}

		std::vector< std::vector < Real > > clusters(original.dimensionsize());
		std::vector< Real > all_dunn_indices;

		//initial cluster state
		for(UInt i = 0; i < clusters.size(); ++i)
		{
			clusters[i].push_back(i);
		}

		//initial state for min inter and max intra distances
		Real min_intercluster_distance(tree[0].distance), max_intracluster_distance(0);
		std::pair<UInt,UInt> clusters_with_min_intercluster_dist(tree[0].left_child,tree[0].right_child);
		for(UInt cluster_step = 0; cluster_step < tree.size()-1; ++cluster_step)
		{

			//max intracluster distance changed?
			for(UInt x = 0; x < clusters[tree[cluster_step].left_child].size();++x)
			{
				for(UInt y = 0; y < clusters[tree[cluster_step].right_child].size();++y)
				{
					if(original.getValue(clusters[tree[cluster_step].left_child][x],clusters[tree[cluster_step].right_child][y]) > max_intracluster_distance)
					{
						max_intracluster_distance = original.getValue(clusters[tree[cluster_step].left_child][x],clusters[tree[cluster_step].right_child][y]);
					}
				}
			}

			//redo clustering following tree
			//pushback elements of right_child to left_child (and then erase second)
			clusters[tree[cluster_step].left_child].insert(clusters[tree[cluster_step].left_child].end(),clusters[tree[cluster_step].right_child].begin(),clusters[tree[cluster_step].right_child].end());

			//erase second one
			clusters.erase(clusters.begin()+tree[cluster_step].right_child);

			//shortcut for single linkage generated hirachy as merging criterion is min intercluster distance
			if(tree_from_singlelinkage)
			{
				min_intercluster_distance = tree[cluster_step+1].distance;
			}
			else
			{
				if(tree[cluster_step].left_child == clusters_with_min_intercluster_dist.first || tree[cluster_step].right_child == clusters_with_min_intercluster_dist.second)
				{
					//find new min intercluster distance
					min_intercluster_distance = std::numeric_limits<Real>::max();
					for(UInt k = 1; k < clusters.size(); ++k)
					{
						for(UInt l = 0; l < k; ++l)
						{
							for(UInt i = 0; i < clusters[k].size(); ++i)
							{
								for(UInt j = 0; j < clusters[l].size(); ++j)
								{
									if(original.getValue(clusters[k][i],clusters[l][j]) < min_intercluster_distance)
									{
										min_intercluster_distance=original.getValue(clusters[k][i],clusters[l][j]);
										clusters_with_min_intercluster_dist.first = k;
										clusters_with_min_intercluster_dist.second = l;
									}
								}
							}
						}
					}
					if(clusters_with_min_intercluster_dist.first > clusters_with_min_intercluster_dist.second)
					{
						std::swap(clusters_with_min_intercluster_dist.first , clusters_with_min_intercluster_dist.second);
					}
				}
			}

			/* to manually retrace
				std::cout << min_intercluster_distance << std::endl;
				std::cout << clusters_with_min_intercluster_dist.first << " , " << clusters_with_min_intercluster_dist.second << std::endl;
				std::cout << max_intracluster_distance << std::endl;
			*/

			//find max dunn index and deduct the corresponding cluster step
			Real dunn_index;
			if(max_intracluster_distance>0)
			{
				dunn_index=min_intercluster_distance/max_intracluster_distance;
			}
			else
			{
				dunn_index=0;
			}
			all_dunn_indices.push_back(dunn_index);
		}
		all_dunn_indices.push_back(0.0); //last one is clearly 0
		return all_dunn_indices;
	}

	void ClusterAnalyzer::cut(size_t cluster_quantity, std::vector< std::vector<UInt> >& clusters, std::vector<BinaryTreeNode>& tree)
	{
		if(cluster_quantity==0)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "minimal partition contains one cluster, not zero");
		}
		if(cluster_quantity>=tree.size()+1)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "maximal partition contains singleton clusters, further separation is not possible");
		}
		clusters.clear();
		clusters.reserve(tree.size()+1);
		for(UInt i = 0; i < tree.size()+1; ++i)
		{
			clusters.push_back(std::vector<UInt>(1,i));
		}
		//redo clustering till step (original.dimensionsize()-cluster_quantity)
		for(UInt cluster_step = 0; cluster_step < tree.size()+1-cluster_quantity; ++cluster_step)
		{
			//pushback elements of right_child to left_child (and then erase second)
			clusters[tree[cluster_step].left_child].insert(clusters[tree[cluster_step].left_child].end(),clusters[tree[cluster_step].right_child].begin(),clusters[tree[cluster_step].right_child].end());

			// erase second one
			clusters.erase(clusters.begin()+tree[cluster_step].right_child);
		}
	}

	Real ClusterAnalyzer::averagePopulationAberration(UInt cluster_quantity, std::vector<BinaryTreeNode>& tree)
	{
		if(cluster_quantity==0)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "minimal partition contains one cluster, not zero");
		}
		if(cluster_quantity>=tree.size()+1)
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "maximal partition contains singleton clusters, further separation is not possible");
		}
		std::vector<Real> average_sizes;
		average_sizes.reserve(tree.size()+1);
		std::vector< std::vector<UInt> > clusters;
		clusters.reserve(tree.size()+1);

		clusters.clear();
		clusters.reserve(tree.size()+1);
		for(UInt i = 0; i < tree.size()+1; ++i)
		{
			clusters.push_back(std::vector<UInt>(1,i));
		}
		//redo clustering till step (original.dimensionsize()-cluster_quantity)
		for(UInt cluster_step = 0; cluster_step < tree.size()+1-cluster_quantity; ++cluster_step)
		{
			//pushback elements of right_child to left_child (and then erase second)
			clusters[tree[cluster_step].left_child].insert(clusters[tree[cluster_step].left_child].end(),clusters[tree[cluster_step].right_child].begin(),clusters[tree[cluster_step].right_child].end());

			// erase second one
			clusters.erase(clusters.begin()+tree[cluster_step].right_child);
		}

		Real average = (Real)(tree.size()+1)/(Real)cluster_quantity;
		Real aberration(0);
		for(UInt i = 0; i < clusters.size(); ++i)
		{
			aberration += std::fabs((Real)clusters[i].size()-average);
		}
		aberration /= (Real)clusters.size();

		return aberration;
	}

	std::vector< Real > ClusterAnalyzer::cohesion(std::vector< std::vector<UInt> >& clusters, DistanceMatrix<Real>& original)
	{
		if(clusters.size()==0 || clusters.size() > original.dimensionsize())
		{
			throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "invalid clustering");
		}

		std::vector< Real > cohesions;
		cohesions.reserve(clusters.size());
		for(UInt i = 0; i < clusters.size(); ++i)
		{
			std::vector< Real > se; // calculate squared errors for cluster i
			se.reserve(clusters[i].size());
			Real av_dist(0); // all pairwise distances in cluster i
			for(UInt j = 0; j < clusters[i].size(); ++j)
			{
				for(UInt k = 0; k < j; ++k)
				{
					se.push_back(original.getValue(clusters[i][j],clusters[i][k]));
					av_dist += original.getValue(clusters[i][j],clusters[i][k]);
				}
			}
			av_dist /= (((Real)clusters[i].size()*(Real)(clusters[i].size()-1.0))/2.0); //now av. intra cluster distance
			Real sse(0);
			for(UInt d = 0; d < se.size(); ++d)
			{
				sse += ((se[d]-av_dist) * (se[d]-av_dist));
			}
			cohesions.push_back(sse);
		}
		return cohesions;
	}

	String ClusterAnalyzer::newickTree(std::vector<BinaryTreeNode>& tree, bool include_distance)
	{

		std::vector<String> clusters;
		clusters.reserve(tree.size()+1);
		for(UInt i = 0; i < tree.size()+1; ++i)
		{
			clusters.push_back(String(i));
		}
		//redo clustering till step (original.dimensionsize()-1)
		for(UInt cluster_step = 0; cluster_step < tree.size(); ++cluster_step)
		{
			//append string right_child to left_child (and then erase second)
			if(include_distance)
			{
				clusters[tree[cluster_step].left_child].insert(0,"( ");
				clusters[tree[cluster_step].left_child] += ":" ;
				clusters[tree[cluster_step].left_child] += tree[cluster_step].distance ;
				clusters[tree[cluster_step].left_child] += " , ";
				clusters[tree[cluster_step].left_child] += clusters[tree[cluster_step].right_child] ;
				clusters[tree[cluster_step].left_child] += ":" ;
				clusters[tree[cluster_step].left_child] += String(tree[cluster_step].distance) ;
				clusters[tree[cluster_step].left_child] += " )";
			}
			else
			{
				clusters[tree[cluster_step].left_child].insert(0,"( ");
				clusters[tree[cluster_step].left_child] += " , ";
				clusters[tree[cluster_step].left_child] += clusters[tree[cluster_step].right_child] ;
				clusters[tree[cluster_step].left_child] += " )";
			}

			// erase second one
			clusters.erase(clusters.begin()+tree[cluster_step].right_child);
		}

		return clusters[0];
	}

	bool compareBinaryTreeNode(const BinaryTreeNode& x, const BinaryTreeNode& y)
	{
		return (x.distance < y.distance);
	}

	BinaryTreeNode::BinaryTreeNode(size_t i, size_t j, Real x) : left_child(i), right_child(j), distance(x)
	{
	}

	BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& source) : left_child(source.left_child), right_child(source.right_child), distance(source.distance)
	{
	}

	BinaryTreeNode::~BinaryTreeNode()
	{
	}

	BinaryTreeNode& BinaryTreeNode::operator = (const BinaryTreeNode& source)
	{
		if (this != &source)
		{
			left_child = source.left_child;
			right_child = source.right_child;
			distance = source.distance;
		}
		return *this;
	}

}
