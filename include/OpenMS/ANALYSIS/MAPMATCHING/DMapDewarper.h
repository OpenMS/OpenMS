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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DMAPDEWARPER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DMAPDEWARPER_H

#include<OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>

#include<OpenMS/KERNEL/DFeatureMap.h>

#include<vector>

namespace OpenMS
{

  /**
    @brief This class applies a transformation as computed by
    class DBaseMapMatcher to a set of features.
        
  */
  template < typename MapT = DFeatureMap<2, DFeature<2> > >
  class DMapDewarper
  {
    public:

      /** @name Type definitions
      */
      //@{
      /// The grid is simply a vector of cells.
      typedef DGrid<2> Grid;
      ///
      typedef MapT MapType;
      ///
      typedef DBaseMapping<1> MappingType;
      ///
      typedef std::vector< MappingType* > MappingVector;
      //@}

      /// Constructor
      DMapDewarper()
          : grid_(), elements_()
      {}

      /// Copy constructor
      DMapDewarper(const DMapDewarper& source)
          : grid_(source.grid_),
          elements_(source.elements_)
      {}

      ///  Assignment operator
      DMapDewarper& operator = (const DMapDewarper& source)
      {
        if (&source == this)
          return *this;

        grid_     = source.grid_;
        elements_ = source.elements_;
        return *this;
      }

      /// equality operator
      bool operator == (const DMapDewarper& rhs)
      {
        return (grid_     == rhs.grid_ &&
                elements_ == rhs.elements_);
      }

      /// Destructor
      virtual ~DMapDewarper()
      {}

      /**
			
        @brief Dewarps the feature map. 
        
        This is a bit ugly. The mapping is D-dimensional function of
        DPosition<D>. But we work with two one dimensional functions
        (one for rt and one for m/z). Therefore we take the DPosition<2>
        of each feature and split it into two instances of DPosition<1> ....
        
      */
      void dewarp()
      {
        // iterate over all features...
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
              DPosition<2> pos         = feat_iter->getPosition();
              // apply transform for every coordinate
              for (unsigned int i=0; i<2; i++)
              {
                DPosition<1> temp;
                temp[0] = pos[i];
								//std::cout << "Retrieving mapping " << i << std::endl;
                grid_iter->getMappings().at(i)->apply(temp);
                pos[i] = temp[0];
              }
              feat_iter->setPosition(pos);
            }
            grid_iter++;

          } // end while (grid)

          feat_iter++;
        } // end while (features)

      }

      /** @name Accesssor methods
      */
      //@{
      /// Set grid
      void setGrid(Grid& g) 
      {
        grid_ = g;
      }
      /// Get grid
      Grid& getGrid()
      {
        return grid_;
      }
      /// Get grid (non-mutable)
      const Grid& getGrid() const
      {
        return grid_;
      }
      /// Set features
      void setMap(MapType& elem)
      {
        elements_ = elem;
      }
      /// Get map
      MapType& getMap()
      {
        return elements_;
      }
      /// Get map (non-mutable)
      const MapType& getMap() const
      {
        return elements_;
      }
      //@}

    protected:
      /// Vector of DRange instances defining a grid over the map
      Grid grid_;

      /// The feature that we want to dewarp
      MapType elements_;

  }
  ; // end of class DMapDewarper

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DMAPDEWARPER_H
