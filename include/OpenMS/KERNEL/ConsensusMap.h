// -*- Mode: C++; tab-width: 2; -*-
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
    
    A ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements and have basically the same interface
    as an STL vector (model of Random Access Container and Back Insertion Sequence).
 
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
      filenames_ = source.filenames_;
      return *this;
    }

    /// Non-mutable access to the filenames
    inline const  std::vector < String >& getFileNames() const
    {
      return filenames_;
    }
    /// Mutable access to the filenames

    inline std::vector < String >& getFileNames()
    {
      return filenames_;
    }
    /// Mutable access to filenames
    inline void setFileNames(const std::vector < String >& filenames)
    {
      filenames_ = filenames;
    }
    
      /// Merge overlapping consensus elements
    void merge(ConsensusMap<ConsensusElementT>& new_map)
    {
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
           
          ConsensusElementType& cons_elem = this->operator[](it->first);
          
          for (UInt i = 0; i < it->second.size(); ++i)
          {
            DoubleReal mz_min = this->operator[](it->second[i]).getPositionRange().min()[1];
            DoubleReal mz_max = this->operator[](it->second[i]).getPositionRange().max()[1];

            UInt duplicates = 0;
            if (((mz_min < c_mz_min) && (c_mz_min < mz_max)) || ((c_mz_min < mz_min) && (mz_min < c_mz_max)))
            {
              for (ConsensusElementIterator it_cons = (this->operator[](it->second[i])).begin(); it_cons != (this->operator[](it->second[i])).end(); ++it_cons)
              {
                IndexTuple i;
                i.setMapIndex(it_cons->getMapIndex());
                if (cons_elem.find(i) != cons_elem.end())
                  {
                    ++duplicates;
                  } 
              }

              //               if (duplicates < (std::min(cons_elem.size(),this->operator[](it->second[i]).size())*0.2))
              if (duplicates > 0)
              {
                continue;
              }
              else
              {
                // insert
                for (ConsensusElementIterator it_cons = (this->operator[](it->second[i])).begin(); it_cons != (this->operator[](it->second[i])).end(); ++it_cons)
                {
                  cons_elem.insert(*it_cons);
                }
                
                // remove 
                std::cout << "Remove " << *(this->begin() + it->second[i]) << std::endl;
                remove_indices.push_back(it->second[i]);      	
               }
            }
          }
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
              UInt curr = remove_indices[j];
              ++j;
              while (remove_indices[j] == curr)
              {
                ++j;
              }
              continue;
            }
            else
              {
                new_map.push_back(this->operator[](i));
              }
        }
        
    }

  protected:
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
