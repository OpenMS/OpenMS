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

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/DPeakArray.h>

namespace OpenMS
{
  /**
    @brief A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements (e.g. consensus features, peaks or raw data points)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as vectors of elements and have basically the same interface
    as an STL vector has (model of Random Access Container and Back Insertion Sequence).
    
    @ingroup Kernel, Serialization
  */
  template < typename ConsensusElementT = ConsensusFeature <DFeatureMap< 2, DFeature<2, KernelTraits > > > >
  class ConsensusMap : public DPeakArray< 2, ConsensusElementT >
  {
    public:
      /// Consensus element type
      typedef ConsensusElementT ConsensusElementType;
      /// Base class type
      typedef DPeakArray< 2, ConsensusElementType > Base;
      /// Mutable iterator
      typedef typename Base::iterator Iterator;
      /// Non-mutable iterator
      typedef typename Base::const_iterator ConstIterator;
      /// Mutable reverse iterator
      typedef typename Base::reverse_iterator ReverseIterator;
      /// Non-mutable reverse iterator
      typedef typename Base::const_reverse_iterator ConstReverseIterator;

      /** @name Constructors and Destructor
      */
      //@{
      /// See DelementArray documentation.
      inline ConsensusMap() : Base()
      {}
      /// See DelementArray documentation.
      inline ConsensusMap(const ConsensusMap& p) : Base(p)
      {
        map_vector_ = p.map_vector_;
        filenames_ = p.filenames_;

      }
      /// See DelementArray documentation.
      inline ~ConsensusMap()
      {}
      /// See DelementArray documentation.
      ConsensusMap(typename Base::size_type n) : Base(n)
      {}
      /// See DelementArray documentation.
      ConsensusMap(typename Base::size_type n, const ConsensusElementType& element) : Base(n,element)
      {}
      /// See DelementArray documentation.
      template <class InputIterator>
      ConsensusMap(InputIterator f, InputIterator l) : Base(f,l)
      {}
      //@}

      /// See DelementArray documentation.
      ConsensusMap& operator = (const ConsensusMap& rhs)
      {
        if (this==&rhs)
          return *this;

        Base::operator=(rhs);
        map_vector_ = rhs.map_vector_;
        filenames_ = rhs.filenames_;
        return *this;
      }


      /// Non-mutable access to the maps
      inline const std::vector < typename ConsensusElementType::ElementContainerType >& getMapVector() const
      {
        return map_vector_;
      }
      /// Mutable access to the maps
      inline std::vector < typename ConsensusElementType::ElementContainerType >& getMapVector()
      {
        return map_vector_;
      }
      /// Mutable access to the maps
      inline void setMapVector(const std::vector < typename ConsensusElementType::ElementContainerType >& map_vector)
      {
        return map_vector_ = map_vector;
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
        return filenames_ = filenames;
      }

    protected:
      /// Vector of aligned element maps
      std::vector < typename ConsensusElementType::ElementContainerType > map_vector_;

      /// Vector of aligned element maps
      std::vector < String > filenames_;
  };

  ///Print the contents to a stream.
  template < typename ConsensusElementT >
  std::ostream& operator << (std::ostream& os, const ConsensusMap<ConsensusElementT>& cons_map)
  {
    for (UnsignedInt i = 0; i < cons_map.size(); ++i)
    {
      os << cons_map[i] << std::endl;
    }

    return os;
  }
} // namespace OpenMS

#endif // OPENMS_KERNEL_ConsensusMap_H
