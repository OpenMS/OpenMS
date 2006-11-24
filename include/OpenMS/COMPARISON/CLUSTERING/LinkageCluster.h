// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_LINKAGECLUSTER_H
#define OPENMS_COMPARISON_CLUSTERING_LINKAGECLUSTER_H

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>

#include <vector>

namespace OpenMS
{

  /**
  implements the hierarchical Clustering Algorithms
  Complete, Average and Single Linkage <br>
  
    \param linkage: <br>
      0 = complete Linkage<br>
      1 = average Linkage<br>
      2 = single Linkage<br>
    \param cutoff
      this is where the tree is cut into Clusters<br>
    \param masdiff
      maximum mass difference two spectra can have and still be compared, pairs
      with a larger mass difference automatically get similarity zero
    \param split
      zero if all charge states should be clustered in one step<br>
      one if each charge state should be clustered individually (more efficient)<br>
    \param stepsize
      since only spectra with similar parentmass are considered similar the clustering
      can be done in steps of size <i>stepsize (nr of spectra in each step)</i>.
      large stepsizes minimize needed recalculation of overlapping clusters, small
      stepsizes use less ram.<br>
      if the stepsize is smaller than the largest overlap, the clustering doesnt
      succeed
    */
  class LinkageCluster :public ClusterFunctor
  {
  public:

    /** @brief standard constructor <br> */
    LinkageCluster();

    /** @brief copy constructor <br> */
    LinkageCluster(const LinkageCluster& source);

    /** @brief destructor <br> */
    ~LinkageCluster() {}

    /**@brief assignment operator <br> */
    LinkageCluster& operator=(const LinkageCluster& source);

    /** @brief create instance <br> */
    static ClusterFunctor* create() {return new LinkageCluster();}

    /** @brief function call operator <br> */
    std::map<int,ClusterNode*> operator()(const ClusterExperiment::ClusterRun* );

    String info() const;

		static const String getName()
		{
			return "LinkageCluster";
		}

  };

}

#endif // OPENMS_COMPARISON_CLUSTERING_LINKAGECLUSTER_H
