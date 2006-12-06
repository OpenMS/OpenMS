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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERIDENTIFICATOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERIDENTIFICATOR_H

#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>
#include <OpenMS/METADATA/Identification.h>

#include <map>

namespace OpenMS
{
	class DBAdapter;
	class Cluster;
	
  /**
  compares spectra against cluster representatives <br>
  masstollerance is the allowed difference between the spectrum and min and mass masses in the clusters<br>
  there are two ways to use this: <br>
    -quickid is for single queries since the candidate clusters are feteched from database on each call <br>
    -id is to be used in combination with setClusterRun, all clusters are created prior to searching so this is faster for many queries<br>
    - existsAsReference returns if a given peptide can be identified (needs setClusterRun)<br>
  */
  class ClusterIdentificator
  {
  public:

    ClusterIdentificator(DBAdapter* adapterp);

    ClusterIdentificator( const ClusterIdentificator& source );

    ClusterIdentificator& operator=(const ClusterIdentificator& source);

    ~ClusterIdentificator();

    /**
    sets the reference dataset against which to compare
    */
    void setClusterRun( const ClusterExperiment::ClusterRun* crp );

    /**
    sets the maximum mass difference to a cluster when comparing
    */
    void setMassTolerance(double tol);

    /**
    assign sequences of similar clusters <br>
    if setClusterRun was called
    */
    Identification id(const ClusterSpectrum& cspec);

    /**
    assign sequences of similar clusters <br>
    if setClusterRun was not called
    */
    Identification quickid(const ClusterExperiment::ClusterRun& cr ,const ClusterSpectrum& cspec, uint clusterrunid);

    /**
    for evaluation purposes <br>
    check if a particular sequence can be identified by the reference set <br>
    */
    uint existsAsReference(String sequence,uint charge);

    /**
    set if to use the median (consensus)<br>
    or centroid <br>
    */
    void setUseMedian(bool b){usemedian_ = b;}
  private:
    ClusterIdentificator();
    std::vector<const Cluster*> getCandidateClusters_(double mass, uint charge);
    DBAdapter* adapterp_;
    const ClusterExperiment::ClusterRun* clusterrunp_;
    double masstolerance_;
    std::multimap<double,int> clustermass2Clusterid_;
    std::map<int,Cluster*> clusterid2Clusterp_;
    bool usemedian_;

  };
}

#endif // OPENMS_COMPARISON_CLUSTERING_ClusterIdentificator_H
