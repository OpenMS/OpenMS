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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEALIGNMENT_H

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/GeomHashPairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher_registerChildren.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
  /**
     @brief Base alignment class.
     
     This class is the base class for the alignment of multiple element maps (e.g. feature maps). 
     Corresponding elements (e.g. features) are grouped together and stored as a consensus
     element (e.g. consensus feature). 
     
     @todo default values e.g. for the pairfinder...
   
  **/
  template < typename ConsensusElementT >
  class BaseAlignment
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** @name Type definitions
       */
      //@{
      /// Consensuselement
      typedef ConsensusElementT ConsensusElementType;

      /// ContainerType
      typedef typename ConsensusElementType::ElementType ElementType;
      typedef typename ConsensusElementType::ElementContainerType ElementContainerType;

      typedef std::vector< ConsensusElementType > ConsensusMapType;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      BaseAlignment()
          : param_(),
          reference_map_index_(0)
      {}

      /// Copy constructor
      BaseAlignment(const BaseAlignment& source)
      {
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusMapType > >::create((source.pairwise_matcher_)->getName());
        pair_finder_ = Factory<BasePairFinder< ConsensusMapType > >::create((source.pair_finder_)->getName());
      }

      ///  Assignment operator
      virtual BaseAlignment& operator = (const BaseAlignment& source)
      {
        if (&source==this)
          return *this;

        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusMapType> >::create((source.pairwise_matcher_)->getName());
        pair_finder_ = Factory<BasePairFinder< ConsensusMapType > >::create((source.pair_finder_)->getName());
        return *this;
      }

      /// Destructor
      virtual ~BaseAlignment()
    {}
      //@}

      /** @name Accesssor methods
       */
      //@{

      /// Mutable access to the param object
      void setParam(const Param& param)
      {
        param_ = param;

        DataValue data_value = param_.getValue("matchingAlgorithm");
        pairwise_matcher_ = Factory<BasePairwiseMapMatcher< ConsensusMapType > >::create(data_value);

        data_value = param_.getValue("consensusAlgorithm");
        pair_finder_ = Factory<BasePairFinder< ConsensusMapType > >::create(data_value);

        map_type_ = param_.getValue("mapType");
      }
      /// Non-mutable access to the param object
      const Param& getParam() const
      {
        return param_;
      }

      /// Mutable access to the vector of maps
      void setElementMapVector(const std::vector< ElementContainerType* >& element_map_vector)
      {
        element_map_vector_ = element_map_vector;
      }
      /// Mutable access to the vector of maps
      std::vector< ElementContainerType* >& getElementMapVector()
      {
        return element_map_vector_;
      }
      /// Non-mutable access to the vector of maps
      const std::vector< ElementContainerType* >& getElementMapVector() const
      {
        return element_map_vector_;
      }

      /// Mutable access to the index of the reference map
      void setReferenceMapIndex(UnsignedInt index) throw (Exception::InvalidValue)
      {
        if (index > element_map_vector_.size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The index is not contained in the vector of element containers.","") ;
        }
        else
        {
          reference_map_index_ = index;
        }
      }
      /// Non-mutable access to the index of the reference map
      const UnsignedInt& getReferenceMapIndex() const
      {
        return reference_map_index_;
      }

      /// Mutable access to the map vector
      void setElementMap(UnsignedInt index, const ElementContainerType& element_map ) throw (Exception::InvalidValue)
      {
        if (index > element_map_vector_.size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The index is not contained in the vector of element containers.","") ;
        }
        else
        {
          element_map_vector_[index] = &element_map;
        }
      }
      /// Non-mutable access to the map vector
      const ElementContainerType& getElementMap(UnsignedInt index) const throw (Exception::InvalidValue)
      {
        if (index > element_map_vector_.size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The index is not contained in the vector of element containers.","") ;
        }
        else
        {
          return *(element_map_vector_[index]);
        }
      }

      /// Mutable access to the vector of file names
      void setFileNames(const std::vector< String >& file_names)
      {
        file_names_ = file_names;
      }
      /// Mutable access to the vector of file names
      std::vector< String >& getFileNames()
      {
        return file_names_;
      }
      /// Non-mutable access to the vector of file names
      const std::vector< String >& getFileNames() const
      {
        return file_names_;
      }

      /// Mutable access to the map type
      void setMapType(const String& map_type)
      {
        map_type_ = map_type;
      }
      /// Mutable access to the map type
      String& getMapType()
      {
        return map_type_;
      }
      /// Non-mutable access to the map type
      const String& getMapType() const
      {
        return map_type_;
      }
      
      /// Mutable access to the final consensus map
      void setFinalConsensusMap(const std::vector < ConsensusElementType >& final_consensus)
      {
        final_consensus_map_ = final_consensus;
      }
      /// Mutable access to the final consensus map
      std::vector < ConsensusElementType >& getFinalConsensusMap()
      {
        return final_consensus_map_;
      }
      /// Non-mutable access to the final consensus map
      const std::vector < ConsensusElementType >& getFinalConsensusMap() const
      {
        return final_consensus_map_;
      }

      /// Start the alignment
      virtual void run() throw (Exception::InvalidValue) = 0;

      /// Return the alignment tree in Newick Tree format
      virtual String getAlignmentTree() const = 0;

    protected:
      /** @name Data members
       */
      //@{
      /// Parameter
      Param param_;

      /// Final consensus map
      std::vector < ConsensusElementType > final_consensus_map_;

      /// File names of the maps to be aligned
      std::vector< String > file_names_;

      /// The maps to be aligned
      std::vector< ElementContainerType* > element_map_vector_;

      /// The map type
      String map_type_;

      /// Index of the reference map
      UnsignedInt reference_map_index_;

      /// Pairwise map matcher
      BasePairwiseMapMatcher< ConsensusMapType >* pairwise_matcher_;

      /// Pairfinder used to find corresponding elements (consensus)
      BasePairFinder< ConsensusMapType >* pair_finder_;
      //@}

      // take the map with the most elements as a reference map
      void assignReferenceMap_()
      {
        UnsignedInt n = element_map_vector_.size();
        UnsignedInt ref_index = 0;
        UnsignedInt max_number = element_map_vector_[ref_index]->size();

        for (UnsignedInt i = 1; i < n; ++i)
        {
          if (n > max_number)
          {
            max_number = n;
            ref_index = i;
          }
        }
        reference_map_index_ = ref_index;
      }

      // Build a consensus map of the map with index map_index. Take each element and
      // create a consenus element of it (the set of grouped elements contains only the element itself).
      void buildConsensusMapTypeInsertInGroup_(UnsignedInt map_index, std::vector< ConsensusElementType >& cons_map)
      {
        const ElementContainerType& map = *(element_map_vector_[map_index]);
        UnsignedInt n = map.size();
        for (UnsignedInt i=0; i < n; ++i)
        {
          ConsensusElementType c(map_index,i,map[i]);
          cons_map.push_back(c);
        }
      }


      // Build a consensus map of the map with index map_index. Take each element and
      // create a consenus element of it (don't insert the element into the group).
      void buildConsensusMapType_(UnsignedInt map_index, std::vector< ConsensusElementType >& cons_map)
      {
        const ElementContainerType& map = *(element_map_vector_[map_index]);
        UnsignedInt n = map.size();
        for (UnsignedInt i=0; i < n; ++i)
        {
          ConsensusElementType c(map[i].getPosition(),map[i].getIntensity());
          cons_map.push_back(c);
        }
      }

      //       int dumpConsensusMapTypes(const String& filename, const std::vector< ConsensusMapType >& pairwise_consensus_maps_)
      //       {
      //         // for all pairwise consensus maps
      //         for ( Size m = 0; m < pairwise_consensus_maps_.size(); ++m )
      //         {
      //           //           std::cout << "Write consensus map " << m << std::endl;
      //           String name = filename + '_' + String(m+1);
      //           std::ofstream dump_file(name.c_str());
      //           dump_file << "# " << name<< " generated " << Date::now() << std::endl;
      //           dump_file << "# 1:mapnumber 2:firstRT 3:firstMZ 4:firstIT 5:secondRT 6:secondMZ 7:secondIT\n";
      //           const ConsensusMapType& act_map = pairwise_consensus_maps_[m];
      //           // write each consensus with groupsize > 1
      //           for ( Size c = 0; c < act_map.size(); ++c)
      //           {
      //             //           std::cout << "write " << act_map[c].getPosition() << " Count " << act_map[c].count() << std::endl;
      //             if (act_map[c].count() > 1)
      //             {
      //               const ConsensusElementType& cons = act_map[c];
      //               typename ConsensusElementType::Group::const_iterator it = cons.begin();
      //               dump_file << m+1 << ' '
      //               << (it->getElement()).getPosition()[RT] << ' '
      //               << (it->getElement()).getPosition()[MZ] << ' '
      //               << (it->getElement()).getIntensity() << ' ';
      //               std::cout << m+1 << ' '
      //               << (it->getElement()).getPosition()[RT] << ' '
      //               << (it->getElement()).getPosition()[MZ] << ' '
      //               << (it->getElement()).getIntensity() << ' ';
      //               ++it;
      //               dump_file << (it->getElement()).getPosition()[RT] << ' '
      //               << (it->getElement()).getPosition()[MZ] << ' '
      //               << (it->getElement()).getIntensity() << ' '
      //               << std::endl;
      //               std::cout << (it->getElement()).getPosition()[RT] << ' '
      //               << (it->getElement()).getPosition()[MZ] << ' '
      //               << (it->getElement()).getIntensity() << ' '
      //               << std::endl;
      //             }
      //           }
      //           dump_file << "# " << filename << " EOF " << Date::now() << std::endl;
      //         }
      //
      //
      //         return 0;
      //       }
  }
  ; // BaseAlignment
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BaseAlignment_H
