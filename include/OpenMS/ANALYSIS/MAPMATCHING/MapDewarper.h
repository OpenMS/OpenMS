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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPDEWARPER_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPDEWARPER_H

#include<OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/Grid.h>

#include<OpenMS/KERNEL/FeatureMap.h>
#include<vector>

namespace OpenMS
{

  /**
    @brief This class applies a transformation as computed by
    class BaseMapMatcher to a set of features.
        
  */
  template < typename MapType = FeatureMap< Feature > >
  class MapDewarper
  {
  public:

    /// Mappings
    typedef std::vector< BaseMapping* > MappingVector;

    /// Constructor
    MapDewarper()
        : grid_(), elements_()
    {}

    /// Copy constructor
    MapDewarper(const MapDewarper& source)
        : grid_(source.grid_),
        elements_(source.elements_)
    {}

    ///  Assignment operator
    MapDewarper& operator = (const MapDewarper& source);

    /// Destructor
    virtual ~MapDewarper()
    {}

    /// Dewarps the element map.
    void dewarp();

    /** @name Accesssor methods
    */
    //@{
    /// Set grid
    inline void setGrid(Grid& g)
    {
      grid_ = g;
    }
    /// Get grid
    inline Grid& getGrid()
    {
      return grid_;
    }
    /// Get grid (non-mutable)
    inline const Grid& getGrid() const
    {
      return grid_;
    }
    /// Set features
    inline void setMap(MapType& elem)
    {
      elements_ = elem;
    }
    /// Get map
    inline MapType& getMap()
    {
      return elements_;
    }
    /// Get map (non-mutable)
    inline const MapType& getMap() const
    {
      return elements_;
    }
    //@}

  protected:
    /// Vector of DRange<2> instances defining a grid over the map
    Grid grid_;

    /// The feature that we want to dewarp
    MapType elements_;
  }
  ; // end of class MapDewarper

  template < typename MapType >
  MapDewarper<MapType>& MapDewarper<MapType>::operator= (const MapDewarper& source)
  {
    if (&source == this)
      return *this;

    grid_     = source.grid_;
    elements_ = source.elements_;
    return *this;
  }

  template < typename MapType >
  void MapDewarper<MapType>::dewarp()
  {
    // iterate over all elements...
    typename MapType::iterator feat_iter = elements_.begin();
    while (feat_iter != elements_.end() )
    {
      // Test in which cell this feature is included
      // and apply the corresponding transformation
      typename Grid::iterator grid_iter = grid_.begin();
      while (grid_iter != grid_.end() )
      {
        if (grid_iter->encloses(feat_iter->getPosition() ) )
        {
          DPosition<2> pos = feat_iter->getPosition();
          // apply transform for every coordinate
          for (unsigned int i=0; i < 2; i++)
          {
            DoubleReal temp;
            temp = pos[i];
            //std::cout << "Retrieving mapping " << i << std::endl;
            grid_iter->getMappings().at(i)->apply(temp);
            pos[i] = temp;
          }
          feat_iter->setPosition(pos);
        }
        grid_iter++;

      } // end while (grid)

      feat_iter++;
    }
  }// end while (features)
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_MAPDEWARPER_H
