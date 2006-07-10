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
// $Id: Cluster.h,v 1.15 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTER_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTER_H

#include <OpenMS/COMPARISON/CLUSTERING/ClusterNode.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterExperiment.h>

#include <vector>
#include <map>

namespace OpenMS
{
	class DBAdapter;
  
  /**
  a Cluster of MSMS spectra<br>
  consists of a ClusterNode, median and centroid object <br>
  */
  class Cluster: public PersistentObject
  {
    public:

    /** @brief constructor <br>*/
    Cluster( ClusterNode* clusterp, const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp );
    /** @brief constructor <br> */
    Cluster( const ClusterNode& clusterr, const ClusterExperiment::ClusterRun* crp , DBAdapter* adapterp );
    /** @brief constructor <br> */
    Cluster( const std::vector<int>& dataset, MSSpectrum< DPeak<1> > median, MSSpectrum< DPeak<1> > centroid);

    /** @brief copy constructor <br> */
    Cluster( const Cluster& source);

    /** @brief destructor <br> */
    ~Cluster();

    /** @brief readonly access ClusterNode <br> */
    const ClusterNode* datasetp() const { return datasetp_;}

    /** @brief number of spectra <br> */
    uint size() const;

    /** @brief number of spectra with most represented sequence <br> */
    uint sequencecount() const;

    /** @brief sequence most spectra have <br> */
    const String& sequence() const;

    /** @brief access to ClusterNode <br> */
    ClusterNode* datasetp(){return datasetp_;};

    /** @brief median is constructed using SpectrumCheapDPCorr <br> */
    const MSSpectrum< DPeak<1> >& median() const;

    /** @brief centroid is the spectrum that has the least distance to all <br> */
    const MSSpectrum< DPeak<1> >& centroid() const;

    /** @brief write access to median*/
    MSSpectrum< DPeak<1> >& median() ; 
    /** @brief write access to centroid*/
    MSSpectrum< DPeak<1> >& centroid() ;

    /** @brief assignment operator <br> */
    Cluster& operator=(const Cluster& source);

    /** @brief for speed reasons most of ClusterNode is usually not loaded from DB */
    void loadClusterNode(DBAdapter* adapterp);

    /** @brief sequences of objects in Cluster <br> */
    std::map<int,String> sequences(DBAdapter* adapterp);

    /** @brief minimum parent mass of spectra in cluster <br> */
    double getMinParentMass() const { return minmass_;}

    /** @brief maximum parent mass of spectra in cluster <br> */
    double getMaxParentMass() const { return minmass_;}

		// Docu in base class
		virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base);
		
		// Docu in base class
		virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base);

    /** @brief standard constructor */
    Cluster(); //TODO Persistence : make private again

	protected:
		// Docu in base class
    virtual void clearChildIds_()
    {
    	//TODO Persistence	
    };		

  private:



    /** @brief find centroid <br> */
    MSSpectrum< DPeak<1> >* findcentroid_( const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp );

    /** @brief construct median <br> */
    MSSpectrum< DPeak<1> >* findmedian_( const ClusterExperiment::ClusterRun* crp, DBAdapter* adapterp );

    void updatesizes_(DBAdapter* adapterp);

    /**
    holds ids of members
    */
    ClusterNode* datasetp_;

    /**
    median spectrum
    */
    MSSpectrum< DPeak<1> >* medianp_;

    /**
    centroid spectrum
    */
    MSSpectrum< DPeak<1> >* centroidp_;

    double minmass_;

    double maxmass_;

    uint size_;

    uint sequencecount_;

    String sequence_;

  };
}

#endif // OPENMS_COMPARISON_CLUSTERING_CLUSTER_H
