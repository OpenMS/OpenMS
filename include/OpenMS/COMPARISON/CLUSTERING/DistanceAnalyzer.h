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
#ifndef OPENMS_COMPARISON_CLUSTERING_DISTANCEANALYZER_H
#define OPENMS_COMPARISON_CLUSTERING_DISTANCEANALYZER_H

#include <OpenMS/COMPARISON/CLUSTERING/AnalysisFunctor.h>

namespace OpenMS
{
  /**
  DistanceAnalyzer uses the preprocessing and compare settings in the <br>
  other ClusterRun to analyze the clustering given to operator() <br>
  
  \param clustersize_cutoff minimum size a cluster has to have to be considered
  \param mass_cutoff maximum mass difference for intercluster pairs to be considered
  \param resolution of the histogram
  \return
    ROC curve for intracluster vs intercluster similarities with AUC <br>
    histograms of inter and intracluster similarities<br>
  */
  class DistanceAnalyzer : public AnalysisFunctor
  {
  public:
    /** @brief standard constructor <br>*/
    DistanceAnalyzer();

    /** @brief copy constructor <br> */
    DistanceAnalyzer(const DistanceAnalyzer& source);

    /** @brief destructor <br>*/
    ~DistanceAnalyzer();

    /** @brief assignment operator <br> */
    DistanceAnalyzer& operator=(const DistanceAnalyzer& source);

    static FactoryProduct* create() {return new DistanceAnalyzer();}

    std::map<String,double> operator()(const std::map<int,ClusterNode*>&);

		static const String getName()
		{
			return "DistanceAnalyzer";
		}

  };
}
#endif  // OPENMS_COMPARISON_CLUSTERING_CLUSTERCOMPAREFUNCTOR_H
