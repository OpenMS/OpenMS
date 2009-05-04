// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_FEATUREDECHARGER_H
#define OPENMS_ANALYSIS_DECHARGING_FEATUREDECHARGER_H

// OpenMS
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/COMPARISON/CLUSTERING/HierarchicalClustering.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

// STL
#include <vector>

namespace OpenMS
{
  /** 
    @brief An algorithm to decharge features (i.e. as found by FeatureFinder).
    
    The Decharger uses an hierarchical clustering (complete linkage) to group charge variants of the same peptide, which
    usually occur in ESI ionization mode. The resulting zero-charge peptides, which are defined by RT and mass
    are written in a featureXML file. Intensities of charge variants are summed up. The position of the zero charge
    variant is the average of all clustered peptides in each dimension.
    If several peptides with the same charge variant are grouped (which is clearly not allowed), a heuristic is used:
    <ul>
    <li>cluster consists of only one charge variant (but several peptides) -> split cluster in single elements</li>
    <li>cluster consists of several charge variants -> dispose cluster</li>
    </ul>
		 
		@htmlinclude OpenMS_FeatureDecharger.parameters
    
    @ingroup Analysis
  */
  
  class OPENMS_DLLAPI FeatureDecharger : public DefaultParamHandler
  {
    public:
    
      typedef FeatureMap<> FeatureMapType;
      typedef Feature FeatureType;
      typedef DPosition<2> ClusterPointType;
      using DefaultParamHandler::param_;
      using DefaultParamHandler::defaults_;
          
      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      FeatureDecharger()
      : DefaultParamHandler("FeatureDecharger")
      {
        defaults_.setValue("cluster_rt_mz_relation", 100.0, "Multiplication factor for m/z coordinates used to balance the dimension differences of RT and m/z.");
        
        HierarchicalClustering<> hc;
        defaults_.insert("hierarchical_clustering:",hc.getParameters());
        
        defaultsToParam_();
      }

      /// Copy constructor
      inline FeatureDecharger(const FeatureDecharger& source)
          : DefaultParamHandler(source),
          featuremap_dc_(source.featuremap_dc_)
      {}

      /// Assignment operator
      inline FeatureDecharger& operator=(const FeatureDecharger& source)
      {
        if (&source==this)
        {
          return *this;
        }

        DefaultParamHandler::operator=(source);
        featuremap_dc_ = source.featuremap_dc_;

        return *this;          
      };

        

      /// destructor
      virtual ~FeatureDecharger() 
      {};
      //@}    

      /** @name Acessors
       */
      //@{

      /// retrieve computed zero-charge feature map
      const FeatureMapType& getFeatureMap() const
      {
        return featuremap_dc_; 
      }

      //@}
      
      /// compute a zero-charge feature map from a set of charged features (@p map)
      void compute(FeatureMapType &map) 
      {

        std::vector<ClusterPointType> feature_stripped;

        // remove charge
        
        double RT_MZ_relation = param_.getValue("cluster_rt_mz_relation");
        
        //std::cout << "creating initial clusters ... ";
        for (FeatureMapType::iterator iter = map.begin(); iter!=map.end(); ++iter)
        {
          ClusterPointType point;
          double mass = iter->getMZ() * iter->getCharge() - iter->getCharge();   //TODO: rectify by real Proton weight (1.07...)
          iter->setMZ(mass);
          point.setX( iter->getRT() );
          point.setY( mass*RT_MZ_relation );
          //std::cout << "x:y " << point.getX() << " " << point.getY() << "\n";
          feature_stripped.push_back(point);
        }

        // cluster
        
        HierarchicalClustering<ClusterPointType> hierclust;
        Param mod_param = param_.copy("hierarchical_clustering:",true);
        if (mod_param.empty()) 
        {
          std::cout << "HierarchicalClustering: param is emtpy. Using defaults!\n";
        }
        else
        {
          hierclust.setParameters(mod_param);
        }
        
        hierclust.compute(feature_stripped);
    
        HierarchicalClustering<ClusterPointType>::ClusterIdxVectorType clusters = hierclust.getClusters();

        hierclust.printStatistics(std::cout);
        
        // combine all features which belong to the same cluster
        FeatureType feature;
        featuremap_dc_.clear();
        featuremap_dc_.assign(clusters.size(), feature);

        uint idx_validCluster = 0;
        bool is_bad_cluster = false;
        uint bad_clusters = 0;
        uint bad_clusters_resolved = 0;                
        for (uint i=0; i<clusters.size(); ++i)
        {
          
          double rt_avg = 0;
          double m_avg = 0;
          double int_sum = 0;
          
          is_bad_cluster = false;
          // 
          std::vector < int > charge_variants(clusters[i].size());
          //std::cout << "cluster size: " << clusters[i].size() << "\n";
          for (uint j=0; j<clusters[i].size() ;++j)
          {
            rt_avg += map[clusters[i][j]].getRT();
            m_avg += map[clusters[i][j]].getMZ();
            int_sum += map[clusters[i][j]].getIntensity();
            
            //std::cout << "  adding " << map[clusters[i][j]].getRT() << "   " << map[clusters[i][j]].getMZ() << "   "<< map[clusters[i][j]].getIntensity() << "\n";
            
            // store charge of current feature
            charge_variants[j] = map[clusters[i][j]].getCharge();
          }
          
          // check if charges within current cluster are unique
          std::sort(charge_variants.begin(), charge_variants.end());
          for (uint j = 1; j<clusters[i].size() ;++j)
          {
            // check if a charge variant appears more than once
            if (charge_variants[j-1]==charge_variants[j])
            {
              is_bad_cluster = true;
              ++bad_clusters;
              // if the cluster has only one charge variant ...
              if (charge_variants[0]==charge_variants[clusters[i].size()-1])
              {
                ++bad_clusters_resolved;
                // ... we can split it
                FeatureType feature;
                featuremap_dc_.insert(featuremap_dc_.end(), clusters[i].size(), feature);
                
                // append single elements of current cluster
                // --> does not work? why
                //std::vector < ClusterPointType > newCluster(1); 
                //clusters.insert(clusters.end(), clusters[i].size(), newCluster);
                
                for (uint new_cl = 0; new_cl<clusters[i].size(); ++new_cl)
                {
                  featuremap_dc_[idx_validCluster].setRT (map[clusters[i][new_cl]].getRT() );
                  featuremap_dc_[idx_validCluster].setMZ ( map[clusters[i][new_cl]].getMZ() );
                  featuremap_dc_[idx_validCluster].setIntensity ( map[clusters[i][new_cl]].getIntensity() );
                  featuremap_dc_[idx_validCluster].setCharge ( 0 );        
                  ++idx_validCluster;          
                }
              }  
            }               
          }
          
          if (is_bad_cluster == false)
          {
            featuremap_dc_[idx_validCluster].setRT (rt_avg / clusters[i].size() );
            featuremap_dc_[idx_validCluster].setMZ ( m_avg / clusters[i].size() );
            featuremap_dc_[idx_validCluster].setIntensity ( int_sum );
            featuremap_dc_[idx_validCluster].setCharge ( 0 );
            
            //TODO average over quality as well?
            //feature.setQuality(0,1); // override default
            //feature.setQuality(1,1); // override default
            //feature.setOverallQuality(1); // override default
            ++idx_validCluster;
          }
          
        }

        // erase all elements past idx_validCluster-1
        featuremap_dc_.erase(featuremap_dc_.begin()+idx_validCluster, featuremap_dc_.end());    
        
  
        std::cout << "STATISTICS:\n  #valid cluster (incl. recovered):" << idx_validCluster << "\n  #badCluster:" << bad_clusters << "\n  #Cluster recovered from bad:" << bad_clusters_resolved << "\n";
        
        featuremap_dc_.updateRanges();
        
        return;
      }

    protected:
      /// result feature map (zero-charge)
      FeatureMapType featuremap_dc_;     
  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DECHARGING_FEATUREDECHARGER_H
