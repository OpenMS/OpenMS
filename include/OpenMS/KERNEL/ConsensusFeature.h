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

#include <OpenMS/KERNEL/RawDataPoint2D.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/Feature.h>

#include <set>

namespace OpenMS
{
  /**
    @brief A 2-dimensional consensus feature.
    
    A consensus feature represents corresponding features in multiple featuremaps.
    
    @todo do not call computeConsensus_() automatically. only calculate it on demand (Marc)
    
    @ingroup Kernel
  */
  class ConsensusFeature
  	: public RawDataPoint2D,
    	public std::set<FeatureHandle, FeatureHandle::IndexLess>
  {
    public:
      /**
        @name Type definitions
      */
      //@{
      typedef std::set<FeatureHandle, FeatureHandle::IndexLess> HandleSetType;
      typedef DRange<2> PositionBoundingBoxType;
      typedef DRange<1> IntensityBoundingBoxType;
      //@}


      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      inline ConsensusFeature()
      	: RawDataPoint2D(),
          HandleSetType(),
          position_range_(),
          intensity_range_(),
          quality_(0.0)
      {
      }
      
      /// Copy constructor
      inline ConsensusFeature(const ConsensusFeature& rhs)
      	: RawDataPoint2D(rhs),
          HandleSetType(rhs),
          position_range_(rhs.position_range_),
          intensity_range_(rhs.intensity_range_),
          quality_(rhs.quality_)
      {
      }
      
      ///Constructor from raw data point
      inline ConsensusFeature(const RawDataPoint2D& point)
      	: RawDataPoint2D(),
          HandleSetType(),
          position_range_(),
          intensity_range_(),
          quality_()
      {
      	this->RawDataPoint2D::operator=(point);
      }

      /// Constructor with map and element index for a singleton consensus feature group
      inline ConsensusFeature(UInt map_index,  UInt element_index, const Feature& feature)
    	{
        this->insert(map_index,element_index,feature);
      }

      /// Assignement operator
      ConsensusFeature& operator=(const ConsensusFeature& source)
      {
        if (&source==this) return *this;

        HandleSetType::operator=(source);
        RawDataPoint2D::operator=(source);
        position_range_=source.position_range_;
        intensity_range_=source.intensity_range_;

        return *this;
      }

      /// Destructor
      virtual ~ConsensusFeature()
   		{
   		}
      //@}

			/**
				@brief Adds an feature handle into the consensus feature
      	
      	@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
      */
      inline void insert(const FeatureHandle& handle, bool recalculate = true)
      {
        if (!(HandleSetType::insert(handle).second))
        {
        	String key = String("map") + handle.getMapIndex() + "/feature" + handle.getElementIndex();
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",key) ;
        }
        
        if (recalculate) computeConsensus_();
      }
			
			/**
				@brief Creates an FeatureHandle and adds it
      	
      	@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
      */
      inline void insert(UInt map_index, UInt feature_index, const Feature& feature, bool recalculate = true)
      {
        insert(FeatureHandle(map_index,feature_index,feature),recalculate);
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
      inline const HandleSetType& getFeatures() const
      {
        return *this;
      }

      /// Returns the quality
      inline DoubleReal getQuality() const
      {
        return quality_;
      }

      /// Sets the quality
      inline void setQuality(DoubleReal quality)
      {
        quality_ = quality;
      }

    protected:
    	///Position range
      PositionBoundingBoxType position_range_;
      ///Intensity range
      IntensityBoundingBoxType intensity_range_;
      ///Quality of the consensus feature
      DoubleReal quality_;

      // compute the consensus attributes like intensity and position as well as the position and intensity range given by the group elements
      void computeConsensus_()
      {
        unsigned int n = HandleSetType::size();
        DPosition<2> sum_position;
        DPosition<2> pos_min(std::numeric_limits<DoubleReal>::max());
        DPosition<2> pos_max(std::numeric_limits<DoubleReal>::min());
        DPosition<1> sum_intensities = 0;
        DPosition<1> int_min(std::numeric_limits<DoubleReal>::max());
        DPosition<1> int_max(std::numeric_limits<DoubleReal>::min());
        for (HandleSetType::const_iterator it = HandleSetType::begin(); it != HandleSetType::end(); ++it)
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
