// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_DATASUBSET_H
#define OPENMS_DATASTRUCTURES_DATASUBSET_H

#include <OpenMS/DATASTRUCTURES/DataPoint.h>
#include <OpenMS/DATASTRUCTURES/SILACTreeNode.h>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>

namespace OpenMS
{
  class DataSubset;

  /**
   * @brief structure for Distance
   */
   struct Dist
   {

   };

  /**
   * @brief Element of a DistanceSet
   * @see HashClustering
   * @ingroup Datastructures
   */
  struct DistanceEntry
  {
  /**
	 * @brief DataSubset to which the distance points
	 */
   DataSubset* data_point;

  /**
	 * @brief DataSubset from which the distance points, i.e which holds an iterator for the distance entry
	 */
   DataSubset* owner;

  /**
	 * @brief The distance
	 */
   DoubleReal distance;

  /**
	 * @brief detailed constructor
	 * @param owner_ DataSubset to which the distance points
   * @param data_point_ DataSubset from which the distance points, i.e which holds an iterator for the distance entry
	 * @param distance_ the distance
	 */
   DistanceEntry(DataSubset* owner_, DataSubset* data_point_, DoubleReal distance_)
    {
      owner = owner_;
      data_point = data_point_;
      distance = distance_;
    }

    Int operator<(const DistanceEntry &i) const
    {
      return distance < i.distance;
    }

  };

  typedef boost::multi_index::multi_index_container<
      DistanceEntry,
      boost::multi_index::indexed_by<
      boost::multi_index::hashed_unique<
      boost::multi_index::composite_key<
      DistanceEntry,
      boost::multi_index::member<DistanceEntry, DataSubset*,&DistanceEntry::owner>,
      boost::multi_index::member<DistanceEntry, DataSubset*,&DistanceEntry::data_point> > >,
      boost::multi_index::ordered_non_unique< boost::multi_index::tag<Dist>,
      boost::multi_index::member<DistanceEntry, DoubleReal,&DistanceEntry::distance> > >
      > DistanceSet;

  /**
   * @brief A DataSubset is a data structure used for hierarchical clustering based on geometric hashing.
   *
   *
   * A DataSubset represents a subset of DataPoints arranged in a HashGrid, as well as an subtree of the hierarchical clustering tree.
   * @image html DataSubset.png
   *
   * @see HashClustering
   * @ingroup Datastructures
   */

  class OPENMS_DLLAPI DataSubset : public GridElement
  {

  public:
  /**
	 * @brief Map of iterators to the distance set
	 * @see HashClustering
	 */
   std::map<GridElement*, DistanceSet::iterator> distance_iterators;

  /**
   * All data points contained in the subset
   */
   std::list<DataPoint*> data_points;

  /**
	 * Subtree of the hierarchical clustering tree representing the data points
	 */
   std::vector<SILACTreeNode> tree;

  /**
   * @brief default constructor
   */
   DataSubset();

  /**
   * @brief destructor
   */
   ~DataSubset();

  /**
	 * @brief detailed constructor
	 * @param data_point initial data point
	 */
   DataSubset(DataPoint& data_point);

  /**
	 * @brief copy constructor
	 * @param copy this DataSubset will be copied
	 */    
   DataSubset(const DataSubset& copy);

  /**
	 * @brief copy constructor
	 * @param copy_ptr the DataSubset, to which the pointer points, will be copied
	 */    
   DataSubset(const DataSubset* copy_ptr);

   Int operator<(const DataSubset &el) const;

  /**
	 * @brief gets the number of DataPoints in the subset
	 */    
   Size size();

  /**
	 * @brief gets the id of the DataSubset
	 */
  Int getID() const;

    bool operator != (const DataSubset &el) const;
    bool operator == (const DataSubset &el) const;

  };
}

#endif /* DATASUBSET_H_ */
