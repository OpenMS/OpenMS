// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>

#include <deque>

namespace OpenMS
{

	#define DEBUG_MERGING
	#undef DEBUG_MERGING

	/**
    @brief A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature or ConsensusPeak)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements and have basically the same interface
    as an STL vector (model of Random Access Container and Back Insertion Sequence).
 
    @improvement use STL list instead of vector (because of insertion and deletion of elements) (Eva)
    
    @ingroup Kernel
  */
template < typename ConsensusElementT = ConsensusFeature < FeatureMap< > > >
class ConsensusMap : public DPeakArray<ConsensusElementT >
  {
  public:
    /// Consensus element type
    typedef ConsensusElementT ConsensusElementType;
    /// Consensus element iterator
    typedef typename ConsensusElementType::Group::const_iterator ConsensusElementIterator;
    /// Base class type
    typedef DPeakArray<ConsensusElementType > Base;
    /// Mutable iterator
    typedef typename Base::iterator Iterator;
    /// Non-mutable iterator
    typedef typename Base::const_iterator ConstIterator;
    /// Mutable reverse iterator
    typedef typename Base::reverse_iterator ReverseIterator;
    /// Non-mutable reverse iterator
    typedef typename Base::const_reverse_iterator ConstReverseIterator;
    
    /// Constructor
    inline ConsensusMap()
        : Base()
    {}

    /// Copy constructor
    inline ConsensusMap(const ConsensusMap& source)
        : Base(source),
        map_vector_(source.map_vector_),
        filenames_(source.filenames_)
    {    }

    /// Destructor
    inline ~ConsensusMap()
    {}

    /// Creates a ConsensusMap with n elements
    ConsensusMap(typename Base::size_type n) : Base(n)
    {}

    /// Assignment operator
    ConsensusMap& operator = (const ConsensusMap& source)
    {
      if (this==&source)
        return *this;

      Base::operator=(source);
      map_vector_ = source.map_vector_;
      filenames_ = source.filenames_;
      return *this;
    }

    /// Non-mutable access to the maps
    inline const std::vector < typename ConsensusElementType::ElementContainerType* >& getMapVector() const
    {
      return map_vector_;
    }
    /// Mutable access to the maps
    inline std::vector < typename ConsensusElementType::ElementContainerType* >& getMapVector()
    {
      return map_vector_;
    }
    /// Mutable access to the maps
    inline void setMapVector(const std::vector < typename ConsensusElementType::ElementContainerType* >& map_vector)
    {
      map_vector_ = map_vector;
    }

    /// Non-mutable access to the filenames
    inline const  std::vector < String >& getFilenames() const
    {
      return filenames_;
    }
    /// Mutable access to the filenames

    inline std::vector < String >& getFilenames()
    {
      return filenames_;
    }
    /// Mutable access to filenames
    inline void setFilenames(const std::vector < String >& filenames)
    {
      filenames_ = filenames;
    }
    
      /// Merge overlapping consensus elements
    void merge(ConsensusMap<ConsensusElementT>& new_map)
    {
#ifdef DEBUG_MERGING  	
      std::cout << "Number of elements " << this->size() << std::endl;
#endif
  
      std::vector<UInt> remove_indices;     	
      UInt n = this->size();
      std::vector < std::pair < DoubleReal, UInt > > rt_start_end_points(2*n);
      UInt i,j;
      for (i = 0,j = 0; i < n; ++i, j+=2)
      {
        rt_start_end_points[j] = std::pair< DoubleReal, UInt> (this->operator[](i).getPositionRange().min()[0], i);
        rt_start_end_points[j+1] = std::pair< DoubleReal,	 UInt> (this->operator[](i).getPositionRange().max()[0], i);
      }
      sort(rt_start_end_points.begin(),rt_start_end_points.end());
      
     
      // scan line for overlapping consensus ranges in retention time
      std::map < UInt, std::vector< UInt> > active_intervalls;
      std::map < UInt, std::vector< UInt> > intersecting_intervalls;
        
#ifdef DEBUG_MERGING  	
      std::cout << "Merge " << std::endl;
#endif
      n = rt_start_end_points.size();
      // find consensus elements, which intersect in retention time
      for (UInt i = 0; i < n; ++i)        
        {
          UInt current_index = rt_start_end_points[i].second;
            
          // check if the intervall with index curr_index is already active
          // if not
          std::map < UInt, std::vector< UInt> >::iterator it = active_intervalls.find(current_index);
          if (it == active_intervalls.end())
          {
            // iterate over all active elements and push the new intersecting intervall
            it = active_intervalls.begin();
            while (it != active_intervalls.end())
            {
              it->second.push_back(current_index);
              ++it;
            }
            
            // insert this intervall to the active
            active_intervalls.insert(std::pair<UInt, std::vector< UInt> > (current_index, std::vector< UInt >() ));
          }
          else
          {
            // store 
            if (it->second.size() > 0)
            {
              intersecting_intervalls.insert(*it);
            }
            // reactivate the appropriate intervall	
            active_intervalls.erase(it);
          }  
        }
        
      
        // find within the groups of intersecting elements in rt, elements that overlap in m/z
        std::map < UInt, std::vector< UInt> >::iterator it = intersecting_intervalls.begin();
        while(it != intersecting_intervalls.end())
        {
          // test if the intervall overlaps in m/z with one of the other intervalls
          DoubleReal c_mz_min = this->operator[](it->first).getPositionRange().min()[1];
          DoubleReal c_mz_max = this->operator[](it->first).getPositionRange().max()[1];
#ifdef DEBUG_MERGING  	
          std::cout << " RT " <<  this->operator[](it->first).getPositionRange().min()[0] 
                        << " - "  <<  this->operator[](it->first).getPositionRange().max()[0]
                        << " MZ  " <<  c_mz_min
                        << " - "  <<  c_mz_max
                        << "\n Indizes " << std::endl;
#endif    
           
          ConsensusElementType& cons_elem = this->operator[](it->first);
          
#ifdef DEBUG_MERGING  	
          for (ConsensusElementIterator it_cons = cons_elem.begin(); it_cons != cons_elem.end(); ++it_cons)
          {
            std::cout << it_cons->getMapIndex() << ' ';
          }
          std::cout << "----- " << std::endl;
#endif
          
          for (UInt i = 0; i < it->second.size(); ++i)
          {
            DoubleReal mz_min = this->operator[](it->second[i]).getPositionRange().min()[1];
            DoubleReal mz_max = this->operator[](it->second[i]).getPositionRange().max()[1];

            UInt duplicates = 0;
            if (((mz_min < c_mz_min) && (c_mz_min < mz_max)) || ((c_mz_min < mz_min) && (mz_min < c_mz_max)))
            {
#ifdef DEBUG_MERGING  	
              std::cout << "+++ RT " <<  this->operator[](it->second[i]).getPositionRange().min()[0] 
                        << " - "  <<  this->operator[](it->second[i]).getPositionRange().max()[0]
                        << " MZ  " <<  mz_min
                        << " - "  <<  mz_max 
                        << "\n Indizes " << std::endl;
#endif                      
              /*for (ConsensusElementIterator it_cons = (this->operator[](it->second[i])).begin(); it_cons != (this->operator[](it->second[i])).end(), duplicates < 1; ++it_cons)
              */
              for (ConsensusElementIterator it_cons = (this->operator[](it->second[i])).begin(); it_cons != (this->operator[](it->second[i])).end(); ++it_cons)
              {
#ifdef DEBUG_MERGING  	
                std::cout << it_cons->getMapIndex() << ' ';
#endif
                IndexTuple<typename ConsensusElementType::ElementContainerType> i;
                i.setMapIndex(it_cons->getMapIndex());
                if (cons_elem.find(i) != cons_elem.end())
                  {
                    ++duplicates;
                  } 
              }
#ifdef DEBUG_MERGING  	
              std::cout << "\n number of duplicates " << duplicates << std::endl;
#endif              	

              //               if (duplicates < (std::min(cons_elem.size(),this->operator[](it->second[i]).size())*0.2))
              if (duplicates > 0)
              {
                continue;
              }
              else
              {
#ifdef DEBUG_MERGING  	              	
                std::cout << "MERGE" << std::endl;
#endif                	
                // insert
                for (ConsensusElementIterator it_cons = (this->operator[](it->second[i])).begin(); it_cons != (this->operator[](it->second[i])).end(); ++it_cons)
                {
                  cons_elem.insert(*it_cons);
#ifdef DEBUG_MERGING  	                  
                  std::cout << "INSERT " << it_cons->getMapIndex() << std::endl;
#endif                  	
                }
                
                // remove 
#ifdef DEBUG_MERGING  	                
                std::cout << "Remove " << *(this->begin() + it->second[i]) << std::endl;
#endif          
                std::cout << "Remove " << *(this->begin() + it->second[i]) << std::endl;
                remove_indices.push_back(it->second[i]);      	
               }
            }
          }
#ifdef DEBUG_MERGING  	          
          std::cout << "========================" << std::endl;
#endif          	
          ++it;
        }
        
       
        sort(remove_indices.begin(),remove_indices.end());
                     
        j=0;
        n=this->size();
        UInt m=remove_indices.size();
        for (UInt i = 0; i < n; ++i)
        {
          if ((j<m) && (i == remove_indices[j]))
            {
              ++j;
              continue;
            }
            else
              {
                new_map.push_back(this->operator[](i));
              }
        }
        
#ifdef DEBUG_MERGING  	
       

        std::cout << "Number of elements " << this->size() << std::endl;
#endif
    }

  protected:
    /// Vector of element maps
    std::vector < typename ConsensusElementType::ElementContainerType* > map_vector_;

    /// Vector of element map filenames
    std::vector < String > filenames_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  template < typename ConsensusElementT >
  std::ostream& operator << (std::ostream& os, const ConsensusMap<ConsensusElementT>& cons_map)
  {
    for (UInt i = 0; i < cons_map.size(); ++i)
    {
      os << cons_map[i] << std::endl;
    }

    return os;
  }
} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
