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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_UNIQUEIDINDEXER_H
#define OPENMS_CONCEPT_UNIQUEIDINDEXER_H

#include <OpenMS/CONCEPT/UniqueIdInterface.h>

#ifdef _MSC_VER // disable some BOOST warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4396 )
#endif

#include <boost/unordered_map.hpp>

#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif


namespace OpenMS
{

/**@brief A base class for random access containers for classes derived from UniqueIdInterface
 * that adds functionality to convert a unique id into an index into the container.
 *
 * See FeatureMap<> and ConsensusMap for living examples.
 * The RandomAccessContainer must support operator[], at(), and size().
 */
template < typename RandomAccessContainer >
  class OPENMS_DLLAPI UniqueIdIndexer
  {
    public:

      typedef boost::unordered_map<UInt64, Size> UniqueIdMap;

      /**@brief Returns the index of the feature with the given unique id, or Size(-1) if none exists in this random access container.

       The complexity is expected constant upon success, linear upon failure.

       The lookup actually performs the following steps:
       - consult the internal hash map (i.e., uniqueid_to_index_)
       - if an index is found and that element indeed has this unique id, then return the index
       - if no index is found or the unique ids do not match, then update the hash map (using updateUniqueIdToIndex()) and lookup the index again
       - if an index is found this time, return it, otherwise return Size(-1)
       .
       Note that subordinate elements are not considered here.
       */
      Size
      uniqueIdToIndex(UInt64 unique_id) const
      {
        Size index;
        try
        {
          index = uniqueid_to_index_.at(unique_id);
          if ( getBase_().at(index).getUniqueId() != unique_id )
          {
            throw std::out_of_range("unique_id_to_index_");
          }
        }
        catch ( std::out_of_range& )
        {
          try
          {
            this->updateUniqueIdToIndex();
            index = uniqueid_to_index_.at(unique_id);
          }
          catch ( std::out_of_range& )
          {
            index = -1; // which means: invalid
          }
        }
        return index;
      }

      /**@brief Updates the hash map from unique id to index.

       @sa uniqueIdToIndex()
       */
      void
      updateUniqueIdToIndex() const
      {
        Size num_valid_unique_id = 0;
        // add or update unique id of existing features
        for ( Size index = 0; index < getBase_().size(); ++index )
        {
          UInt64 unique_id = getBase_()[index].getUniqueId();
          if ( UniqueIdInterface::isValid(unique_id) )
          {
            uniqueid_to_index_[unique_id] = index;
            ++num_valid_unique_id;
          }
        }
        // remove invalid or outdated entries
        uniqueid_to_index_.erase(UniqueIdInterface::INVALID);
        for ( UniqueIdMap::iterator iter = uniqueid_to_index_.begin(); iter != uniqueid_to_index_.end(); /* see loop */)
        {
          if ( iter->second >= getBase_().size() || getBase_()[iter->second].getUniqueId() != iter->first )
          {
            iter = uniqueid_to_index_.erase(iter);
          }
          else
          {
            ++iter;
          }
        }
        if ( uniqueid_to_index_.size() != num_valid_unique_id )
        {
          std::stringstream ss;
          ss << "Duplicate valid unique ids detected!   RandomAccessContainer has size()==" << getBase_().size();
          ss << ", num_valid_unique_id==" << num_valid_unique_id;
          ss << ", uniqueid_to_index_.size()==" << uniqueid_to_index_.size();
          throw Exception::Postcondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, ss.str());
        }
        return;
      }

      /**@brief Swap.
       *
       *  Note: obviously we can swap only with indices for the same type.
       */
      void
      swap(UniqueIdIndexer& rhs)
      {
        std::swap(uniqueid_to_index_, rhs.uniqueid_to_index_);
        return;
      }

    protected:

      /**@brief A little helper to get access to the base (!) class RandomAccessContainer.
       *
       * This is just a static_cast and probably not interesting elsewhere, so we make it a protected member.
       */
      const RandomAccessContainer&
      getBase_() const
      {
        return *static_cast<const RandomAccessContainer*> (this);
      }

      /**@brief hash map from unique id to index of features
       *
       * This is mutable because the hash map is updated on demand, even if the underlying container is const.
       */
      mutable UniqueIdMap uniqueid_to_index_;

  };

} //namespace OpenMS

#endif // OPENMS_CONCEPT_UNIQUEIDINDEXER_H
