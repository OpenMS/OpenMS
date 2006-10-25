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
// $Maintainer: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASESUPERIMPOSER_H

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
     @brief The base class of all superimposer algorithms.

     This class defines the basic interface for all superimposer
     algorithms. It works on two feature maps and estimates the transformation
     which maps one map (scene map) onto the other (model map).

     The matching shall minimize the shift in retention time and m/z between
     the two maps after a suitable transformation, yet have large cardinality.

  **/

  template <typename MapT = DFeatureMap<2> >
  class BaseSuperimposer
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** Symbolic names for indices of feature maps etc.
      This should make things more understandable and maintainable. 
      */
      enum Maps
      {
        MODEL = 0,
        SCENE = 1
    };

      /** @name Type definitions
       */
      //@{
      /// Container for input features
      typedef MapT PointMapType;

      /// Type of features considered here
      typedef typename PointMapType::value_type PointType;

      /// Traits type
      typedef typename PointType::TraitsType TraitsType;

      /// Quality type
      typedef typename TraitsType::QualityType QualityType;

      /// Position type
      typedef DPosition < 2, TraitsType > PositionType;

      //// Intensity type
      typedef typename TraitsType::IntensityType IntensityType;

      /// Type of estimated transformation
      typedef DBaseMapping< 1, TraitsType > TransformationType;

      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      BaseSuperimposer()
          : param_(),
          final_transformation_()
      {
        feature_map_[MODEL] = 0;
        feature_map_[SCENE] = 0;
        final_transformation_[RT] = 0;
        final_transformation_[MZ] = 0;
      }

      /** @brief Copy constructor

      The final transformation member has to be copied in the derived classes!
      */
      BaseSuperimposer(const BaseSuperimposer& source)
          : param_(source.param_)
      {
        feature_map_[MODEL] = source.feature_map_[MODEL];
        feature_map_[SCENE] = source.feature_map_[SCENE];
        final_transformation_[RT] = source.final_transformation_[RT];
        final_transformation_[MZ] = source.final_transformation_[MZ];
      }

      /** @brief  Assignment operator

      The final transformation member has to be assigned in the derived classes!
      */
      virtual BaseSuperimposer& operator = (BaseSuperimposer source)
      {
        param_ = source.param_;
        feature_map_[MODEL] = source.feature_map_[MODEL];
        feature_map_[SCENE] = source.feature_map_[SCENE];
        final_transformation_[RT] = source.final_transformation_[RT];
        final_transformation_[MZ] = source.final_transformation_[MZ];
        return *this;
      }

      /// Destructor
      virtual ~BaseSuperimposer()
      {
        for ( Size dim = 0; dim < 2; ++dim )
        {
          if (final_transformation_[dim] != 0)
          {
            delete final_transformation_[dim];
            final_transformation_[dim] = 0;
          }
        }
      }
      //@}

      /** @name Accesssor methods
       */
      //@{


      /// Set param class
      void setParam(const Param& param)
      {
        param_ = param;
      }

      /// Get param class
      Param& getParam()
      {
        return param_;
      }

      /// Get param class (non-mutable)
      const Param& getParam() const
      {
        return param_;
      }


      /// Set feature map
      void setFeatureMap(Size const index, const PointMapType& feature_map)
      {
        feature_map_[index] = &feature_map;
      }

      /// Get feature map
      const PointMapType& getFeatureMap(Size index)
      {
        return *feature_map_[index];
      }

      /// Get feature maps (non-mutable)
      const PointMapType& getFeatureMap(Size index) const
      {
        return *feature_map_[index];
      }

      //@}

      /// Estimates the transformation for each grid cell
      virtual void run()
      {}
      ;


    protected:

      /** @name Data members
       */
      //@{

      /// Param class containing the parameters for the map matching phase
      Param param_;

      /// Two maps of features to be matched
      PointMapType const * feature_map_[2];

      /// Final transformation
      TransformationType* final_transformation_[2];

      //@}

  }
  ; // BaseSuperimposer

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BaseSuperimposer_H
