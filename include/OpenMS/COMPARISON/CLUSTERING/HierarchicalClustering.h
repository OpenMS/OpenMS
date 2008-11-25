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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H
#define OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H

// OpenMS
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

// STL
#include <vector>
#include <math.h>


namespace OpenMS
{
  /** 
    @brief Complete or single linkage clustering based on DPosition<>.
    
    General clustering in arbitrary dimensions using euclidean distance up until a certain
    cutoff, i.e. the maximal allowed distance for merging two clusters.
        
    The result is an STL vector of vectors. Each inner vector respresents one cluster, holding the indices of the input elements.
		 
		@htmlinclude OpenMS_HierarchicalClustering.parameters
        
    @ingroup Comparison
  */
  
  template <typename ClusterPoint = DPosition<2> >
  class HierarchicalClustering : public DefaultParamHandler
  {
    public:
    
      typedef ClusterPoint ClusterPointType;
      typedef std::vector<uint> ClusterIdxType;
      typedef std::vector< ClusterIdxType > ClusterIdxVectorType;
      using DefaultParamHandler::param_;
      using DefaultParamHandler::defaults_;
      
      typedef std::vector <double> DistanceLineType;
      typedef std::vector < DistanceLineType > DistanceMatrixType;
      
      /// clustering method
      enum LINKAGE_TYPE {COMPLETE_LINKAGE=0, SINGLE_LINKAGE};    
          
      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      HierarchicalClustering()
      : DefaultParamHandler("HierarchicalClustering")
      {
        defaults_.setValue("cluster_cutoff", 40.0, "maximal distance allowed between two clusters before merging them");           // maximal distance allowed between two clusters before merging them
        defaults_.setValue("linkage_type", COMPLETE_LINKAGE,"clustering method: 0=complete linkage 1=single_linkage");
    
        defaultsToParam_();
      }

      /// Copy constructor
      inline HierarchicalClustering(const HierarchicalClustering& source)
          : DefaultParamHandler(source),
          clustermap_(source.clustermap_)
      {}

      /// Assignment operator
      inline HierarchicalClustering& operator=(const HierarchicalClustering& source)
      {
        if (&source==this)
        {
          return *this;
        }

        DefaultParamHandler::operator=(source);
        clustermap_ = source.clustermap_;

        return *this;          
      };

        

      /// destructor
      virtual ~HierarchicalClustering() 
      {
      };
      //@}    

      /** @name Acessors
       */
      //@{

      /// retrieve computed zero-charge feature map
      const ClusterIdxVectorType & getClusters() const
      {
        return clustermap_; 
      }

      //@}
      
      /** @brief compute clusters from a vector of DPositions up until given cutoff
			 *  
			 *	@exception Throws Exception::NotImplemented if linkage_type Param is invalid
			 */
      void compute(const std::vector<ClusterPointType>& points)
      {
        double cutoff = param_.getValue("cluster_cutoff");
        int linkage_method = param_.getValue("linkage_type");

                
        // initial result: every point represents a cluster
        clustermap_.clear();
        ClusterIdxVectorType clustermap_tmp(points.size(), ClusterIdxType(1));
        clustermap_.assign(clustermap_tmp.begin(),clustermap_tmp.end());
        
        for (uint i=0; i<points.size(); ++i)
        {
          clustermap_[i][0] = i;
        }
        
        
        // initial distances between points
        DistanceMatrixType distanceMatrix(points.size(), DistanceLineType(1));
                
        // compute distance matrix
        for (uint i=1; i<points.size(); ++i)
        {
          distanceMatrix[i].resize(i);
          for (uint j=0; j<i; ++j)
          {
            distanceMatrix[i][j] = getDistance_(points[i], points[j]);
            //std::cout << getDistance_(points[i], points[j]) << " ";
          }
          //std::cout << "\n";
        }
        
        for (uint i=0; i<points.size(); ++i)
        {
          //std::cout << "size: " <<  distanceMatrix[i].size();
        }
        
        
        double number_of_clusters = points.size()-1;
        
        double minDistance;
        uint minRow=0, minColumn=0;
        
        getMinDistance_(distanceMatrix, minDistance, minRow, minColumn);
        
        // iterate until cutoff reached
        while ((minDistance < cutoff) && (number_of_clusters>=2))
        {
           //std::cout << "starting ...";

           // find closest two clusters
           getMinDistance_(distanceMatrix, minDistance, minRow, minColumn);
           
           //std::cout << "combining clusters " << (minRow) << " and " << minColumn << "\n";
           
           // combine clusters
           // IDEA? append to larger cluster (and swap afterwards if neccessary)
           while (!clustermap_[minColumn].empty())
           {
             //std::cout << "putting " << clustermap[minColumn].back() << " from " << minColumn << "->" << minRow << "\n";
             clustermap_[minRow].push_back( clustermap_[minColumn].back() );
             clustermap_[minColumn].pop_back();
           }
           
           // partially recompute distance matrix
           if (linkage_method==COMPLETE_LINKAGE)
           {
             linkageComplete_(distanceMatrix, minRow, minColumn);
           } 
           else if (linkage_method==SINGLE_LINKAGE)
           {
             linkageSingle_(distanceMatrix, minRow, minColumn);
           } 
           else
           {
              std::cerr << "Linkage method not available" << linkage_method << "\n";
              throw Exception::NotImplemented(__FILE__,__LINE__,__FUNCTION__);
           }
           
           //delete 2nd clusters row from matrix (don´t bother about columns - they won´t save memory)
           distanceMatrix[minColumn].clear();

           --number_of_clusters;
        }
        
        //std::cout << "HC: remove dead clusters ... \n";
                
        // remove empty indices 
        ClusterIdxVectorType cl_tmp;
        ClusterIdxVectorType::iterator iter = clustermap_.begin();
        while ( iter != clustermap_.end())
        {
          if (iter->empty())
          {
            iter = clustermap_.erase(iter, iter+1);
          } else
          {
            ++iter;
          }       
        }
       
        return;
      }

     
      void printStatistics (std::ostream& os)
      {
        if (clustermap_.empty())
        {
          os << "no clusters defined! call 'compute()' first!\n";
        }
        else
        {
          os << "\nCluster size statistics:\n";
          uint maxsize = 0;
          for (uint i=0; i<clustermap_.size(); ++i)
          {
            if (clustermap_[i].size() > maxsize) maxsize = clustermap_[i].size();
          }

          std::vector<uint> size_dist(maxsize+1);
          for (uint i=0; i<clustermap_.size(); ++i)
          {
            ++size_dist[clustermap_[i].size()];  
          }
          
          for (uint i=0; i<maxsize+1; ++i)
          {
            if (size_dist[i]>0)
              os << size_dist[i] <<  " of size " << i << "\n";  
          }

        }
        return;
      }

    protected:
    
      /// find entry in distance matrix with MINIMAL value and return value and position
      /// @note minRow will always be greater than minColumn
      inline void getMinDistance_(const DistanceMatrixType& distanceMatrix, double& minDistance, uint& minRow, uint& minColumn)
      {
        minDistance = std::numeric_limits<double>::max();

        for (uint i=1; i<distanceMatrix.size(); ++i)
        {
          if (!distanceMatrix[i].empty())
          {
            //std::cout << "line length: " << distanceMatrix[i].size() << " -- ";
            for (uint j=0; j<i; ++j)
            {
              //std::cout << "i,j " << i << "," << j << "- " << distanceMatrix[i][j] << "#  ";
              if (!distanceMatrix[j].empty())
              {
                if (minDistance > distanceMatrix[i][j])
                {
                  minDistance = distanceMatrix[i][j];
                  minRow = i;
                  minColumn = j;
                }
              }
            }  // for j
            //std::cout << "\n";
          } 
        } // for i 
        
        //std::cout << "Distance_min " << minDistance << " row, column " << minRow << "," << minColumn << "\n";
        
        return;        
      }
    
      /// compute euclidean distance between two points
      inline double getDistance_(const ClusterPointType& a, const ClusterPointType& b)
      {
         // (insert new distance calculation here)
         double distance = 0;
         for (uint i=0; i<a.size(); ++i)
         {
            distance += pow(a[i]-b[i], 2);
         }
         
         return sqrt(distance);
      }

      /// complete linkage clustering
      /// we combine clusters @p minrow and @p mincolumn (@p minrow is the new cluster index)
      inline void linkageComplete_(DistanceMatrixType& distanceMatrix, const uint&  minRow, const uint& minColumn)
      {
        // recompute data matrix entries
        for (uint i=0; i<distanceMatrix.size(); ++i)
        {
          // we only have a triangular matrix ... so be careful with indices
          if (i>minRow) // && (i>minColumn) <-- can be assured, because minRow>minColumn!
          {
            distanceMatrix[i][minRow] = std::max(distanceMatrix[i][minRow], distanceMatrix[i][minColumn]); 
          }
          else if (i<minColumn) // && (i<minRow) <-- can be assured, because minRow>minColumn!
          {
            distanceMatrix[minRow][i] = std::max(distanceMatrix[minRow][i], distanceMatrix[minColumn][i]);  
          }
          else if ((i<minRow) && (i>minColumn))
          {  
            distanceMatrix[minRow][i] = std::max(distanceMatrix[minRow][i], distanceMatrix[i][minColumn]);  
          }
          
          
        }
        
        return;
      };

      /// single linkage clustering
      /// we combine clusters @p minrow and @p mincolumn (@p minrow is the new cluster index)
      inline void linkageSingle_(DistanceMatrixType& distanceMatrix, const uint& minRow, const uint& minColumn)
      {
        // recompute data matrix entries
        for (uint i=0; i<distanceMatrix.size(); ++i)
        {
          // we only have a triangular matrix ... so be careful with indices
          if (i>minRow) // && (i>minColumn) <-- can be assured, because minRow>minColumn!
          {
            distanceMatrix[i][minRow] = std::min(distanceMatrix[i][minRow], distanceMatrix[i][minColumn]); 
          }
          else if (i<minColumn) // && (i<minRow) <-- can be assured, because minRow>minColumn!
          {
            distanceMatrix[minRow][i] = std::min(distanceMatrix[minRow][i], distanceMatrix[minColumn][i]);  
          }
          else if ((i<minRow) && (i>minColumn))
          {  
            distanceMatrix[minRow][i] = std::min(distanceMatrix[minRow][i], distanceMatrix[i][minColumn]);  
          }
          
        }
        
        return;
      };    
    
      /// result cluster map with indices to given feature map
      ClusterIdxVectorType clustermap_;
  };
} // namespace OpenMS

#endif // OPENMS_COMPARISON_CLUSTERING_HIERARCHICALCLUSTERING_H   
