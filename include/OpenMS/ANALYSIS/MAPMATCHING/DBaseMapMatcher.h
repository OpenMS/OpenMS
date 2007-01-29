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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DBASEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DBASEMAPMATCHER_H

#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>

#include <utility>

namespace OpenMS
{

  /**
    @brief The base class of the map matching algorithm.
    
    This class defines the basic interface for all mapmatching
    algorithms. It expects a list of feature pairs together
    with a quality value for each pair and the coordinates
    of a grid covering the map.
    
    It then estimates the parameters of a transformation
    describing the shift in retention time and m/z between
    the two maps.  
  */
  template <typename ElementT = DFeature<2> >
  class DBaseMapMatcher
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < LCMS_Tag >::RT,
        MZ = DimensionDescription < LCMS_Tag >::MZ
    };
      /// Element type
      typedef ElementT ElementType;
      /// Traits type
      typedef typename ElementType::TraitsType TraitsType;
      /// The grid is simply a vector of cells.
      typedef DGrid<2> Grid;
      /// The feature pairs are computed by the feature matching class
      typedef DFeaturePairVector< 2, ElementType > FeaturePairVector;
      /// Quality type
      typedef typename TraitsType::QualityType QualityType;

      /// Constructor
      DBaseMapMatcher() : min_quality_(-1)
      {}

      /// Copy constructor
      DBaseMapMatcher(const DBaseMapMatcher& source)
          : grid_(source.grid_),
          feature_pairs_(source.feature_pairs_),
          min_quality_(source.min_quality_)
      {}

      ///  Assignment operator
      DBaseMapMatcher& operator = (const DBaseMapMatcher& source)
      {
        if (&source==this)
          return *this;

        grid_          = source.grid_;
        feature_pairs_ = source.feature_pairs_;
        min_quality_   = source.min_quality_;

        return *this;
      }

      /// Equality operator
      bool operator == (const DBaseMapMatcher& rhs)
      {
        return (grid_          == rhs.grid_ &&
                feature_pairs_ == rhs.feature_pairs_ &&
                min_quality_   == rhs.min_quality_);
      }

      /// Destructor
      virtual ~DBaseMapMatcher()
      {}

      /// Set grid
      void setGrid(const Grid& g)
      {
        grid_ = g;
      }
      /// Get grid
      Grid& getGrid()
      {
        return grid_;
      }
      /// Get grid (non-mutable)
      const Grid & getGrid() const
      {
        return grid_;
      }
      
      /// Set feature pair list
      void setFeaturePairs(const FeaturePairVector& plist)
      {
        feature_pairs_ = plist;
      }
      /// Get feature pair list
      FeaturePairVector& getFeaturePairs()
      {
        return feature_pairs_;
      }
      /// Get feature pair list (non-mutable)
      const FeaturePairVector& getFeaturePairs() const
      {
        return feature_pairs_;
      }
      
      /// Set quality
      void setMinQuality(const QualityType& qu)
      {
        min_quality_ = qu;
      }
      /// Get quality
      QualityType& getMinQuality()
      {
        return min_quality_;
      }
      /// Get quality
      const QualityType& getMinQuality() const
      {
        return min_quality_;
      }

      /// Estimates the transformation for each grid cell
      virtual void estimateTransform() = 0;

    protected:
      /// Vector of DRange instances defining a grid over the map
      Grid grid_;

      /// Vector of pairs of features that have been identified by the feature matcher
      FeaturePairVector feature_pairs_;

      /// Minimum quality that we accept for feature pairs, defined in param class
      QualityType min_quality_;

  }
  ; // end of class BaseMapMatcher

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DBASEMAPMATCHER_H
