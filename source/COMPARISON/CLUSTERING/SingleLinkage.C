// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
	SingleLinkage::SingleLinkage()
	  : ClusterFunctor(), ProgressLogger()
	{
	}

	SingleLinkage::SingleLinkage(const SingleLinkage& source)
	  : ClusterFunctor(source), ProgressLogger()
	{
	}

	SingleLinkage::~SingleLinkage()
	{
	}

	SingleLinkage& SingleLinkage::operator = (const SingleLinkage& source)
	{
		if (this != &source)
		{
			ClusterFunctor::operator = (source);
			ProgressLogger::operator = (source);
		}
		return *this;
	}

	void SingleLinkage::operator () (DistanceMatrix<Real>& original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold /*=1*/) const
	{
		// input MUST have >= 2 elements!
		if(original_distance.dimensionsize()<2)
		{
			throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from only contains one element");
		}

		cluster_tree.clear();
		if(threshold<1)
		{
      LOG_ERROR << "You tried to use Single Linkage clustering with a threshold. This is currently not supported!" << std::endl;
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}

	//SLINK
		std::vector<Size> pi;
		pi.reserve(original_distance.dimensionsize());
		std::vector<Real> lambda;
		lambda.reserve(original_distance.dimensionsize());

		startProgress(0,original_distance.dimensionsize(),"clustering data");

		//initialize first pointer values
		pi.push_back(0);
		lambda.push_back(std::numeric_limits<Real>::max());

		for (Size k = 1; k < original_distance.dimensionsize(); ++k)
		{
			std::vector<Real> row_k;
			row_k.reserve(k);

			//initialize pointer values for element to cluster
			pi.push_back(k);
			lambda.push_back(std::numeric_limits<Real>::max());

			// get the right distances
			for (Size i = 0; i < k; ++i)
			{
				row_k.push_back(original_distance.getValue(i,k));
			}

			//calculate pointer values for element k
			for (Size i = 0; i < k; ++i)
			{
				if(lambda[i] >= row_k[i])
				{
					row_k[pi[i]] = std::min(row_k[pi[i]], lambda[i]);
					lambda[i] = row_k[i];
					pi[i] = k;
				}
				else
				{
					row_k[pi[i]] = std::min(row_k[pi[i]], row_k[i]);
				}
			}

			//update clustering if neccessary
			for (Size i = 0; i < k; ++i)
			{
				if(lambda[i] >= lambda[pi[i]])
				{
					pi[i] = k;
				}
			}
			setProgress(k);
		}

		for (Size i = 0; i < pi.size()-1; ++i)
		{
			//strict order is always kept in algorithm: i < pi[i]
			cluster_tree.push_back(BinaryTreeNode(i,pi[i],lambda[i]));
			//~ std::cout << i << '\n' << pi[i] << '\n' << lambda[i] << std::endl;
		}

		//sort pre-tree
		std::sort(cluster_tree.begin(),cluster_tree.end(),compareBinaryTreeNode);

		// convert -pre-tree to correct format
		for (Size i = 0; i < cluster_tree.size(); ++i)
		{
			if(cluster_tree[i].right_child < cluster_tree[i].left_child)
			{
				std::swap(cluster_tree[i].left_child,cluster_tree[i].right_child);
			}
			for (Size k = i+1; k < cluster_tree.size(); ++k)
			{
				if(cluster_tree[k].left_child == cluster_tree[i].right_child)
				{
					cluster_tree[k].left_child = cluster_tree[i].left_child;
				}
				else if(cluster_tree[k].left_child > cluster_tree[i].right_child)
				{
					--cluster_tree[k].left_child;
				}
				if(cluster_tree[k].right_child == cluster_tree[i].right_child)
				{
					cluster_tree[k].right_child = cluster_tree[i].left_child;
				}
				else if(cluster_tree[k].right_child > cluster_tree[i].right_child)
				{
					--cluster_tree[k].right_child;
				}
			}

		}
		//~ prepare to redo clustering to get all indices for binarytree in min index element representation
		std::vector< std::set<Size> >clusters(original_distance.dimensionsize());
		for (Size i = 0; i < original_distance.dimensionsize(); ++i)
		{
			clusters[i].insert(i);
		}
		for (Size cluster_step = 0; cluster_step < cluster_tree.size(); ++cluster_step)
		{
			Size new_left_child = *(clusters[cluster_tree[cluster_step].left_child].begin());
			Size new_right_child = *(clusters[cluster_tree[cluster_step].right_child].begin());
			clusters[cluster_tree[cluster_step].left_child].insert(clusters[cluster_tree[cluster_step].right_child].begin(),clusters[cluster_tree[cluster_step].right_child].end());
			clusters.erase(clusters.begin()+cluster_tree[cluster_step].right_child);
			std::swap(cluster_tree[cluster_step].left_child, new_left_child);
			std::swap(cluster_tree[cluster_step].right_child, new_right_child);
			if(cluster_tree[cluster_step].left_child > cluster_tree[cluster_step].right_child)
			{
				std::swap(cluster_tree[cluster_step].left_child , cluster_tree[cluster_step].right_child);
			}
		}

		endProgress();
	}

}
