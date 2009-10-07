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

#ifndef OPENMS_KERNEL_CONSENSUSFEATURE_H
#define OPENMS_KERNEL_CONSENSUSFEATURE_H

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <set>

namespace OpenMS
{
	class PeptideIdentification;

	/**
		@brief A 2-dimensional consensus feature.

		A consensus feature represents corresponding features in multiple feature
		maps. The corresponding features are represented a set of @ref
		FeatureHandle instances.  Each ConsensusFeature "contains" zero or more
		FeatureHandles.

		@see ConsensusMap

		@ingroup Kernel
	*/
	class OPENMS_DLLAPI ConsensusFeature
		: public RichPeak2D,
			public std::set<FeatureHandle, FeatureHandle::IndexLess>
	{
	 public:
		///Type definitions
		//@{
		typedef DoubleReal QualityType;
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
			inline bool operator () ( ConsensusFeature const & left, QualityType const & right ) const
			{
				return ( left.getQuality() < right );
			}
			inline bool operator () ( QualityType const & left, ConsensusFeature const & right ) const
			{
				return ( left< right.getQuality() );
			}
			inline bool operator () ( QualityType const & left, QualityType const & right ) const
			{
				return ( left < right );
			}
		};

		/// Compare by size(), the number of consensus elements
		struct SizeLess
			: std::binary_function < ConsensusFeature, ConsensusFeature, bool >
		{
			inline bool operator () ( ConsensusFeature const & left, ConsensusFeature const & right ) const
			{
				return ( left.size() < right.size() );
			}
			inline bool operator () ( ConsensusFeature const & left, UInt64 const & right ) const
			{
				return ( left.size() < right );
			}
			inline bool operator () ( UInt64 const & left, ConsensusFeature const & right ) const
			{
				return ( left< right.size() );
			}
			inline bool operator () ( const UInt64 & left,	const UInt64 & right ) const
			{
				return ( left < right );
			}
		};

		/// Compare by the sets of consensus elements (lexicographically)
		struct MapsLess
			: std::binary_function < ConsensusFeature, ConsensusFeature, bool >
		{
			inline bool operator () ( ConsensusFeature const & left, ConsensusFeature const & right ) const
			{
				return std::lexicographical_compare( left.begin(), left.end(), right.begin(), right.end(), FeatureHandle::IndexLess() );
			}
		};


		///@name Constructors and Destructor
		//@{
		/// Default constructor
		ConsensusFeature();

		/// Copy constructor
		ConsensusFeature(const ConsensusFeature& rhs);

		/// Constructor from raw data point
		ConsensusFeature(const RichPeak2D& point);

		///Constructor from raw data point
		ConsensusFeature(const Peak2D& point);

		///Constructor from raw data point
		ConsensusFeature(const Feature& feature);

		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position and intensity to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt64 map_index,	UInt64 element_index, const Peak2D& element);


		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position, intensity, charge and quality to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt64 map_index,	UInt64 element_index, const Feature& element);

		/**
			@brief Constructor with map and element index for a singleton consensus
			feature. Sets the consensus feature position, intensity, charge and quality to the values
			of @p element as well.
		*/
		ConsensusFeature(UInt64 map_index,	UInt64 element_index, const ConsensusFeature& element);


		/// Assignment operator
		ConsensusFeature& operator=(const ConsensusFeature& rhs);

		/// Destructor
		virtual ~ConsensusFeature();
		//@}


		///@name Management of feature handles
		//@{
		/**
		@brief Adds an feature handle into the consensus feature

		@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(const FeatureHandle& handle);

		/// Adds all feature handles in @p handle_set to this consensus feature.
		void insert(const HandleSetType& handle_set);

		/**
			@brief Creates a FeatureHandle and adds it

			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt64 map_index, UInt64 element_index, const Peak2D& element);

		/**
			@brief Creates a FeatureHandle and adds it

			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt64 map_index, UInt64 element_index, const Feature& element);

		/**
			@brief Creates a FeatureHandle and adds it

			@exception Exception::InvalidValue is thrown if a handle with the same map and element index already exists.
		*/
		void insert(UInt64 map_index, UInt64 element_index, const ConsensusFeature& element);

		/// Non-mutable access to the contained feature handles
		const HandleSetType& getFeatures() const;
		//@}

		///@name Accessors
		//@{
		/// Returns the quality
		QualityType getQuality() const;
		/// Sets the quality
		void setQuality(QualityType quality);
		/// Sets the charge
		void setCharge(Int charge);
		/// Returns the charge
		Int getCharge() const;
		/// Returns the position range of the contained elements
		DRange<2> getPositionRange() const;
		/// Returns the intensity range of the contained elements
		DRange<1> getIntensityRange() const;
		
		//@}

		/**
       @brief Computes and updates the consensus position, intensity, and charge.

       The position and intensity of the contained feature handles is averaged.
       The most frequent charge state wins, while the tie breaking prefers
       smaller (absolute) charges.

       @note This method has to be called explicitly, <i>after</i> adding the feature handles.
       */
		void computeConsensus();

		/**
       @brief Computes the uncharged parent RT & mass, assuming the handles are charge variants.

       The position of the feature handles (decharged) is averaged.
       Intensities are summed up.
       Charge is set to 0.
       Mass calculation: If the given features contain a metavalue "dc_charge_adduct_mass" then this will be used as adduct mass instead of
       weight(H+) * charge.

       @note This method has to be called explicitly, <i>after</i> adding the feature handles.
       */
		void computeDechargeConsensus(const FeatureMap<>& fm);

		/// returns a const reference to the PeptideIdentification vector
		const std::vector<PeptideIdentification>& getPeptideIdentifications() const;

		/// returns a mutable reference to the PeptideIdentification vector
		std::vector<PeptideIdentification>& getPeptideIdentifications();

		/// sets the PeptideIdentification vector
		void setPeptideIdentifications( const std::vector<PeptideIdentification>& peptide_identifications );

	 protected:
		/// Quality of the consensus feature
		QualityType quality_;

		/// Charge of the consensus feature
		Int charge_;

		/// Peptide identifications belonging to the consensus feature
		std::vector<PeptideIdentification> peptide_identifications_;

	};

	///Print the contents of a ConsensusFeature to a stream
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ConsensusFeature& cons);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSFEATURE_H
