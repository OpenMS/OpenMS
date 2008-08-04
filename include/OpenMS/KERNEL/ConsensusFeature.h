// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSFEATURE_H
#define OPENMS_KERNEL_CONSENSUSFEATURE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

#include <set>

namespace OpenMS
{
  /**
		@brief A 2-dimensional consensus feature.
	    
		A consensus feature represents corresponding features in multiple feature
		maps. The corresponding features are represented a set of @ref FeatureHandle instances.
	
		This implements a variant of the <i>composite</i> design pattern: Each
		ConsensusFeature "contains" zero or more FeatureHandles.
	  
	  @see ConsensusMap
	  
		@ingroup Kernel
  */
  class ConsensusFeature
  	: public Peak2D,
    	public std::set<FeatureHandle, FeatureHandle::IndexLess>
  {
	 public:
		/**
		@name Type definitions
		*/
		//@{
		typedef std::set<FeatureHandle, FeatureHandle::IndexLess> HandleSetType;
		//@}

			/// Compare by getQuality()
			struct QualityLess
			: std::binary_function < ConsensusFeature, ConsensusFeature, bool >
			{
				inline bool operator () ( ConsensusFeature const & left, ConsensusFeature const & right ) const
				{
					return ( left.getQuality() < right.getQuality() );
				}
				inline bool operator () ( ConsensusFeature const & left, DoubleReal const & right ) const
				{
					return ( left.getQuality() < right );
				}
				inline bool operator () ( DoubleReal const & left, ConsensusFeature const & right ) const
				{
					return ( left< right.getQuality() );
				}
				inline bool operator () ( DoubleReal const & left, DoubleReal const & right ) const
				{
					return ( left < right );
				}
			};

		///@name Constructors and Destructor
		//@{
		/// Default constructor
		ConsensusFeature()
			: Peak2D(),
				HandleSetType(),
				quality_(0.0),
				charge_(0)
		{
		}
      
		/// Copy constructor
		ConsensusFeature(const ConsensusFeature& rhs)
			: Peak2D(rhs),
				HandleSetType(rhs),
				quality_(rhs.quality_),
				charge_(rhs.charge_)
		{
		}
      
		///Constructor from raw data point
		ConsensusFeature(const Peak2D& point)
			: Peak2D(point),
				HandleSetType(),
				quality_(0.0),
				charge_(0)
		{
		}

		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position and intensity to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt map_index,  UInt element_index, const Peak2D& element)
			: Peak2D(element),
				HandleSetType(),
				quality_(0.0),
				charge_(0)
		{
			insert(map_index,element_index,element);
		}


		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position, intensity, charge and quality to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt map_index,  UInt element_index, const Feature& element)
			: Peak2D(element),
				HandleSetType(),
				quality_(element.getOverallQuality()),
				charge_(element.getCharge())
		{
			insert(map_index,element_index,element);
		}

		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position, intensity, charge and quality to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt map_index,  UInt element_index, const ConsensusFeature& element)
			: Peak2D(element),
				HandleSetType(),
				quality_(element.getQuality()),
				charge_(element.getCharge())
		{
			insert(map_index,element_index,element);
		}

		
		/// Assignment operator
		ConsensusFeature& operator=(const ConsensusFeature& rhs)
		{
			if (&rhs==this) return *this;

			HandleSetType::operator=(rhs);
			Peak2D::operator=(rhs);
			quality_ = rhs.quality_;
			charge_ = rhs.charge_;
			
			return *this;
		}

		/// Destructor
		virtual ~ConsensusFeature()
		{
		}
		//@}


		///@name Management of feature handles
		//@{
		/**
			@brief Adds an feature handle into the consensus feature
	      	
			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(FeatureHandle const & handle)
		{
			if (!(HandleSetType::insert(handle).second))
			{
				String key = String("map") + handle.getMapIndex() + "/feature" + handle.getElementIndex();
				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The set already contained an element with this key.",key) ;
			}
			return;
		}
		
		///Adds all feature handles in @p handle_set to this consensus feature.
		void insert(HandleSetType const & handle_set)
		{
			for (ConsensusFeature::HandleSetType::const_iterator it = handle_set.begin(); it != handle_set.end(); ++it)
			{
				insert(*it);
			}
			return;
		}

		/**
			@brief Creates an FeatureHandle and adds it
					
			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt map_index, UInt element_index, const Peak2D& element)
		{
			insert(FeatureHandle(map_index,element_index,element));
		}

		/**
			@brief Creates an FeatureHandle and adds it
					
			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt map_index, UInt element_index, const Feature& element)
		{
			insert(FeatureHandle(map_index,element_index,element));
		}
		
		/**
			@brief Creates an FeatureHandle and adds it
					
			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt map_index, UInt element_index, const ConsensusFeature& element)
		{
			insert(FeatureHandle(map_index,element_index,element));
		}

		/// Non-mutable access to the contained feature handles
		const HandleSetType& getFeatures() const
		{
			return *this;
		}
		//@}
		
		///@name Accessors
		//@{
		/// Returns the quality
		DoubleReal getQuality() const
		{
			return quality_;
		}
		/// Sets the quality
		void setQuality(DoubleReal quality)
		{
			quality_ = quality;
		}
		/// Sets the charge
		void setCharge(Int charge)
		{
			charge_ = charge;
		}
		/// Returns the charge
		Int getCharge() const
		{
			return charge_;
		}
		/// Returns the position range of the contained elements
		DRange<2> getPositionRange() const
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
		DRange<1> getIntensityRange() const
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
		//@}
		
		/**
			@brief Computes a consensus position and intensity and assigns them to this consensus feature.
			
			The position and intensity of the contained feature handles is simply averaged and assigned to
			this consensus feature.
			
			@note This method has to be called explicitly, after adding all feature handles.
		*/
		void computeConsensus();

	 protected:
		///Quality of the consensus feature
		DoubleReal quality_;
		///Charge of the consensus feature
		Int charge_;
  };

  ///Print the contents of a ConsensusFeature to a stream
  std::ostream& operator << (std::ostream& os, const ConsensusFeature& cons);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSFEATURE_H
