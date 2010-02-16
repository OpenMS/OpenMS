// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_FEATURE_H
#define OPENMS_KERNEL_FEATURE_H

#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>

namespace OpenMS
{

	/** @brief An LC-MS feature.

	The Feature class is used to describe the two-dimensional signal caused by a
	peptide.	It can store a charge state and a list of peptide identifications.
	The area occupied by the Feature in the LC-MS data set is represented by a
	list of convex hulls (one for each isotopic peak).	There is also a convex
	hull for the entire Feature.	The model description can store the parameters
	of a two-dimensional theoretical model of the underlying signal in LC-MS.
	Currently, non-peptidic compounds are also represented as features.

	By convention in %OpenMS, the position of a feature is defined as maximum
	position of the model for the retention time dimension and the mass of the
	monoisotopic peak for the m/z dimension.	The intensity of a feature is
	(proportional to) its total ion count.

	Feature is derived from RichPeak2D.	 Also inherited is a MetaInfoInterface.
	Features as usually contained in a FeatureMap.	See also FeatureHandle and
	ConsensusFeature.

	@ingroup Kernel
	*/
	class OPENMS_DLLAPI Feature
		: public RichPeak2D
	{
	 public:
		///@name Type definitions
		//@{
		///Type of the quality values
		typedef	 DoubleReal QualityType;
		///Charge type
		typedef Int ChargeType;
		//@}

		/** @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		Feature()
			: RichPeak2D(),
				overall_quality_(),
				convex_hulls_(),
				convex_hulls_modified_(true),
				convex_hull_(),
				charge_( 0 ),
				subordinates_()
		{
			std::fill( qualities_, qualities_ + 2, 0 );
		}

		/// Copy constructor
		Feature( const Feature& feature )
			: RichPeak2D( feature ),
				overall_quality_( feature.overall_quality_ ),
				model_desc_( feature.model_desc_ ),
				convex_hulls_( feature.convex_hulls_ ),
				convex_hulls_modified_(feature.convex_hulls_modified_),
				convex_hull_( feature.convex_hull_ ),
				charge_( feature.charge_ ),
				identifications_( feature.identifications_ ),
				subordinates_( feature.subordinates_ )
		{
			std::copy( feature.qualities_, feature.qualities_ + 2, qualities_ );
		}

		/// Destructor
		~Feature()
		{
		}
		//@}

		/// @name Model and Quality methods
		//@{
		/// Non-mutable access to the overall quality
		QualityType getOverallQuality() const
		{
			return overall_quality_;
		}
		/// Set the overall quality
		void setOverallQuality( QualityType q )
		{
			overall_quality_ = q;
		}

		/// Non-mutable access to the quality in dimension c
		QualityType getQuality( Size index ) const
		{
			OPENMS_PRECONDITION( index < 2, "Feature<2>:getQuality(Size): index overflow!" );
			return qualities_[ index ];
		}
		/// Set the quality in dimension c
		void setQuality( Size index, QualityType q )
		{
			OPENMS_PRECONDITION( index < 2, "Feature<2>:setQuality(Size): index overflow!" );
			qualities_[ index ] = q;
		}

		/// Non-mutable access to the model description
		const ModelDescription<2>& getModelDescription() const
		{
			return model_desc_;
		}
		/// Mutable access to the model description
		ModelDescription<2>& getModelDescription()
		{
			return model_desc_;
		}
		/// Set the model description
		void setModelDescription( const ModelDescription<2>& q )
		{
			model_desc_ = q;
		}
		//@}

		/// Non-mutable access to charge state
		const ChargeType& getCharge() const
		{
			return charge_;
		}
		/// Set charge state
		void setCharge( const ChargeType& ch )
		{
			charge_ = ch;
		}

		///@name Convex hulls and bounding box
		//@{
		/// Non-mutable access to the convex hulls
		const std::vector<ConvexHull2D>& getConvexHulls() const
		{
			return convex_hulls_;
		}
		/// Mutable access to the convex hulls of single mass traces
		std::vector<ConvexHull2D>& getConvexHulls()
		{
			convex_hulls_modified_ = true;
			return convex_hulls_;
		}
		/// Set the convex hulls of single mass traces
		void setConvexHulls( const std::vector<ConvexHull2D>& hulls )
		{
			convex_hulls_modified_ = true;
			convex_hulls_ = hulls;
		}
		/**
		@brief Returns the overall convex hull of the feature (calculated from the convex hulls of the mass traces)

		@note the bounding box of the feature can be accessed through the returned convex hull
		*/
		ConvexHull2D& getConvexHull() const;

		/// Returns if the mass trace convex hulls of the feature enclose the position specified by @p rt and @p mz
		bool encloses(DoubleReal rt, DoubleReal mz) const;
		//@}

		/// Assignment operator
		Feature& operator = ( const Feature& rhs );

		/// Equality operator
		bool operator == ( const Feature& rhs ) const;

		/// Compare by getOverallQuality()
		struct OverallQualityLess
			: std::binary_function < Feature, Feature, bool >
		{
			bool operator () ( Feature const & left, Feature const & right ) const
			{
				return ( left.getOverallQuality() < right.getOverallQuality() );
			}
			bool operator () ( Feature const & left, QualityType right ) const
			{
				return ( left.getOverallQuality() < right );
			}
			bool operator () ( QualityType left, Feature const & right ) const
			{
				return ( left < right.getOverallQuality() );
			}
			bool operator () ( QualityType left, QualityType right ) const
			{
				return ( left < right );
			}
		};

		/// returns a const reference to the PeptideIdentification vector
		const std::vector<PeptideIdentification>& getPeptideIdentifications() const
		{
			return identifications_;
		};

		/// returns a mutable reference to the PeptideIdentification vector
		std::vector<PeptideIdentification>& getPeptideIdentifications()
		{
			return identifications_;
		};

		/// sets the PeptideIdentification vector
		void setPeptideIdentifications( const std::vector<PeptideIdentification>& identifications )
		{
			identifications_ = identifications;
		};

		/// immutable access to subordinate features
		const std::vector<Feature>& getSubordinates() const
		{
			return subordinates_;
		}

		/// mutable access to subordinate features
		std::vector<Feature>& getSubordinates()
		{
			return subordinates_;
		}

		/// mutable access to subordinate features
		void setSubordinates(const std::vector<Feature>& rhs)
		{
			subordinates_ = rhs;
		}


    /**@brief Applies a member function of Type to the feature (including subordinates).
       The returned values are accumulated.

       <b>Example:</b>  The following will print the number of features (parent feature and subordinates) with invalid unique ids:
       @code
       Feature f;
       (...)
       std::cout << f.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
       @endcode
       See e.g. UniqueIdInterface for what else can be done this way.
    */
    template < typename Type >
    Size applyMemberFunction( Size (Type::*member_function)() )
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for ( std::vector<Feature>::iterator iter = subordinates_.begin(); iter != subordinates_.end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    /// The "const" variant.
    template < typename Type >
    Size applyMemberFunction( Size (Type::*member_function)() const ) const
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for ( std::vector<Feature>::const_iterator iter = subordinates_.begin(); iter != subordinates_.end(); ++iter)
      {
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }


	 protected:

		/// Overall quality measure of the feature
		QualityType overall_quality_;

		/// Quality measures for each dimension
		QualityType qualities_[ 2 ];

		/// Description of the theoretical model the feature was constructed with
		ModelDescription<2> model_desc_;

		/// Array of convex hulls (one for each mass trace)
		std::vector<ConvexHull2D> convex_hulls_;

		/// Flag that indicates if the overall convex hull needs to be recomputed (i.e. mass trace convex hulls were modified)
		mutable bool convex_hulls_modified_;

		/// Overall convex hull of the feature
		mutable ConvexHull2D convex_hull_;

		/// Charge of the peptide represented by this feature.  The default value is 0, which represents an unknown charge state.
		ChargeType charge_;

		/// Peptide PeptideIdentifications belonging to the feature
		std::vector<PeptideIdentification> identifications_;

		/// subordinate features (e.g. features that the ModelFitter discarded due to inferior quality)
		std::vector<Feature> subordinates_;

	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_FEATURE_H
