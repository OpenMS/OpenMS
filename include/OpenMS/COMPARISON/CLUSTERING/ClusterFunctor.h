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

#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H

#include <map>

#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>

namespace OpenMS
{
  /**
  ClusterFunctor classes partition a set of ClusterSpectrum into natural groups<br>
  */
  class ClusterFunctor : public FactoryProduct
  {
  public:
    /** @brief sparse similarity matrix <br> */
    typedef std::vector<SparseVector* > simmatrix;

    /** @brief standard constructor <br> */
    ClusterFunctor();

    /** @brief copy constructor <br> */
    ClusterFunctor(const ClusterFunctor& source);

    /** @brief destructor <br> */
    virtual ~ClusterFunctor() {}

    /** @brief copy constructor <br> */
    ClusterFunctor& operator=(const ClusterFunctor& source);

		static void registerChildren();
		
    /** @brief function call operator <br> */
    /**
    \param clusterrunp the parent ClusterRun, it computes the similarity score and provides the ClusterSpectrum objects
    \return the finished clusters ( the int is just a side-effect of the clustering-process )
    */
    virtual std::map<int, ClusterNode*> operator()(const ClusterExperiment::ClusterRun* clusterrunp) = 0;
  };

}

#endif // OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H
