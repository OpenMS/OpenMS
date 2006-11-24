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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERCOMPAREFUNCTOR_H

#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>

namespace OpenMS
{
  /**
  ClusterCompareFunctor compares two clusterings<br>
  the result is similar to a contingency table<br>
  a is the number of pairs that are clustered together in both clusterings<br>
  b and c are the numbers of pairs where the pair is coclustered in one clustering but in different clusters in the other<br>
  d is the number of pairs that are in different clusters in both clusterings <br> <br>
  from these values some cluster similarity values are computed<br>
  */
  class ClusterCompareFunctor : public AnalysisFunctor
  {
  public:
    /** @brief standard constructor <br> */
    ClusterCompareFunctor();

    /** @brief copy constructor <br> */
    ClusterCompareFunctor(const ClusterCompareFunctor& source);

    /** @brief destructor <br> */
    ~ClusterCompareFunctor();

    /** @brief assignment operator <br> */
    ClusterCompareFunctor& operator=(const ClusterCompareFunctor& source);

    /** @brief create instance <br> */
    static AnalysisFunctor* create() {return new ClusterCompareFunctor();}

    std::map<String,double> operator()(const std::map<int,ClusterNode*>&);

		static const String getName()
		{
			return "ClusterCompareFunctor";
		}
	
  };
}
#endif  // OPENMS_COMPARISON_CLUSTERING_CLUSTERCOMPAREFUNCTOR_H
