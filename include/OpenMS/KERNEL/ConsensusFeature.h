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

#ifndef OPENMS_KERNEL_CONSENSUSFEATURE_H
#define OPENMS_KERNEL_CONSENSUSFEATURE_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Group.h>

namespace OpenMS
{

  /**
    @brief A 2-dimensional consensus feature.
    
    A consensus feature represents corresponding features in multiple featuremaps.
    
    @ingroup Kernel
  */
  class ConsensusFeature
  	: public Feature,
    	public Group
  {
    public:
      /**
        @name Type definitions
      */
      //@{
      typedef DRange<2> PositionBoundingBoxType;
      typedef DRange<1> IntensityBoundingBoxType;
      //@}


      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      inline ConsensusFeature()
      	: Feature(),
          Group(),
          position_range_(),
          intensity_range_()
      {
      }
      
      /// Copy constructor
      inline ConsensusFeature(const ConsensusFeature& rhs)
      	: Feature(rhs),
          Group(rhs),
          position_range_(rhs.position_range_),
          intensity_range_(rhs.intensity_range_)
      {
      }
      
      ///Constructor from raw data point
      inline ConsensusFeature(const RawDataPoint2D& point)
      	: Feature(),
          Group(),
          position_range_(),
          intensity_range_()
      {
      	this->RawDataPoint2D::operator=(point);
      }

      /// Constructor with map and element index for a singleton consensus feature group
      inline ConsensusFeature(UInt m,  UInt e, const Feature& feature)
    	{
        this->insert(m,e,feature);
      }

      /// Assignement operator
      ConsensusFeature& operator=(const ConsensusFeature& source)
      {
        if (&source==this) return *this;

        Group::operator=(source);
        Feature::operator=(source);
        position_range_=source.position_range_;
        intensity_range_=source.intensity_range_;

        return *this;
      }

      /// Destructor
      virtual ~ConsensusFeature()
   		{
   		}
      //@}


      inline void insert(const IndexTuple& tuple, bool recalculate = true)
      {
        Group::insert(tuple);

        if (recalculate) computeConsensus_();
      }

      inline void insert(UInt map_index, UInt feature_index, const Feature& feature, bool recalculate = true)
      {
        insert(IndexTuple(map_index,feature_index,feature),recalculate);
      }

      /// Non-mutable access to the position range
      inline const PositionBoundingBoxType& getPositionRange() const
      {
        return position_range_;
      }
      /// Mutable access to the position range
      inline PositionBoundingBoxType& getPositionRange()
      {
        return position_range_;
      }
      /// Set the position range
      inline void setPositionRange(const PositionBoundingBoxType& p)
      {
        position_range_= p;
      }

      /// Non-mutable access to the intensity range
      inline const IntensityBoundingBoxType& getIntensityRange() const
      {
        return intensity_range_;
      }
      /// Mutable access to the intensity range
      inline IntensityBoundingBoxType& getIntensityRange()
      {
        return intensity_range_;
      }
      /// Set the intensity range
      inline void setIntensityRange(const IntensityBoundingBoxType& i)
      {
        intensity_range_ = i;
      }

      /// Non-mutable access to the combined features
      inline const Group& getFeatures() const
      {
        return *this;
      }
      /// Mutable access to the combined features
      inline Group& getFeatures() //TODO remove
      {
        return *this;
      }
      /// Set the combined features
      inline void setFeatures(const Group& g)
      {
        Group::operator=(g);
      }

    protected:
      PositionBoundingBoxType position_range_;
      IntensityBoundingBoxType intensity_range_;

      // compute the consensus attributes like intensity and position as well as the position and intensity range given by the group elements
      void computeConsensus_() //TODO do not compute this automatically; only update ranges automatically
      {
        unsigned int n = Group::size();
        DPosition<2> sum_position;
        DPosition<2> pos_min(std::numeric_limits<DoubleReal>::max());
        DPosition<2> pos_max(std::numeric_limits<DoubleReal>::min());
        DPosition<1> sum_intensities = 0;
        DPosition<1> int_min(std::numeric_limits<DoubleReal>::max());
        DPosition<1> int_max(std::numeric_limits<DoubleReal>::min());
        for (Group::const_iterator it = Group::begin(); it != Group::end(); ++it)
        {
          DPosition<1> act_int = it->getIntensity();
          DPosition<2> act_pos = it->getPosition();

          if (int_min > act_int)
          {
            int_min = act_int;
          }
          if (int_max < act_int)
          {
            int_max = act_int;
          }

          for (UInt dim=0; dim < 2; ++dim)
          {
            if (act_pos[dim] > pos_max[dim]) pos_max[dim] = act_pos[dim];
            if (act_pos[dim] < pos_min[dim]) pos_min[dim] = act_pos[dim];
          }

          sum_intensities += act_int;
          sum_position += act_pos;
        }

        for (UInt dim = 0; dim< 2 ; ++dim)
        {
          this->position_[dim] = sum_position[dim] / n;
        }
        this->intensity_ = sum_intensities[0] / n;

        intensity_range_.setMinMax(int_min,int_max);
        position_range_.setMinMax(pos_min,pos_max);
      }
  };

  ///Print the contents of a ConsensusFeature to a stream
  std::ostream& operator << (std::ostream& os, const ConsensusFeature& cons);

} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATURE_H
