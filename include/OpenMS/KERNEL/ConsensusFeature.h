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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/RawDataPoint2D.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/Feature.h>

#include <set>

namespace OpenMS
{
  /**
    @brief A 2-dimensional consensus feature.
    
    A consensus feature represents corresponding features in multiple feature maps.
    
    @improvement In order to speed up multiple calls of getPositionRange() and getIntensitRange(), one could store them in a mutable member variable (Marc, Clemens)
    
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
      //@}


      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      inline ConsensusFeature()
      	: RawDataPoint2D(),
          HandleSetType(),
          quality_(0.0)
      {
      }
      
      /// Copy constructor
      inline ConsensusFeature(const ConsensusFeature& rhs)
      	: RawDataPoint2D(rhs),
          HandleSetType(rhs),
          quality_(rhs.quality_)
      {
      }
      
      ///Constructor from raw data point
      inline ConsensusFeature(const RawDataPoint2D& point)
      	: RawDataPoint2D(point),
          HandleSetType(),
          quality_()
      {
      }

      /// Constructor with map and element index for a singleton consensus feature group. Sets the consensus feature position and intensity to the values of @p feature as well.
      inline ConsensusFeature(UInt map_index,  UInt element_index, const Feature& feature)
      	: RawDataPoint2D(),
          HandleSetType(),
          quality_()
    	{
        RawDataPoint2D::operator=(feature);
        insert(map_index,element_index,feature);
      }

      /// Assignement operator
      ConsensusFeature& operator=(const ConsensusFeature& rhs)
      {
        if (&rhs==this) return *this;

        HandleSetType::operator=(rhs);
        RawDataPoint2D::operator=(rhs);
        quality_ = rhs.quality_;

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
      inline void insert(const FeatureHandle& handle)
      {
        if (!(HandleSetType::insert(handle).second))
        {
        	String key = String("map") + handle.getMapIndex() + "/feature" + handle.getElementIndex();
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",key) ;
        }
      }
			
			/**
				@brief Creates an FeatureHandle and adds it
      	
      	@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
      */
      inline void insert(UInt map_index, UInt feature_index, const Feature& feature)
      {
        insert(FeatureHandle(map_index,feature_index,feature));
      }

      /// Non-mutable access to the contained feature handles
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

      /// Computes a consensus position and intensity from the contained feature handles uses it for this consensus feature
      void computeConsensus();
			
			/// Returns the position range of the contained elements
			inline DRange<2> getPositionRange() const
			{
		  	DPosition<2> min = DPosition<2>::max;
		  	DPosition<2> max = DPosition<2>::min;
		    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		    {
		    	if (it->getRT()<min[0]) min[0]=it->getRT();
		    	if (it->getRT()>max[0]) max[0]=it->getRT();
		    	if (it->getMZ()<min[1]) min[1]=it->getMZ();
		    	if (it->getMZ()>max[1]) max[1]=it->getMZ();
		    }
				return DRange<2>(min,max);
			}
			
			/// Returns the intensity range of the contained elements
			inline DRange<1> getIntensityRange() const
			{
		  	DPosition<1> min = DPosition<1>::max;
		  	DPosition<1> max = DPosition<1>::min;
		    for (ConsensusFeature::HandleSetType::const_iterator it = begin(); it != end(); ++it)
		    {
		    	if (it->getIntensity()<min[0]) min[0]=it->getIntensity();
		    	if (it->getIntensity()>max[0]) max[0]=it->getIntensity();
		    }
				return DRange<1>(min,max);
			}


    protected:
      ///Quality of the consensus feature
      DoubleReal quality_;
  };

  ///Print the contents of a ConsensusFeature to a stream
  std::ostream& operator << (std::ostream& os, const ConsensusFeature& cons);

} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATURE_H
