// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_FEATURE_H
#define OPENMS_KERNEL_FEATURE_H

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
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
        Features as usually are contained in a FeatureMap. See also FeatureHandle and
	ConsensusFeature.

	@ingroup Kernel
	*/
	class OPENMS_DLLAPI Feature
		: public BaseFeature
	{
	 public:
		/** @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		Feature();

		/// Copy constructor
		Feature(const Feature& feature);

		/// Destructor
		~Feature();
		//@}

		/// @name Model and quality methods
		//@{
		/// Non-mutable access to the overall quality
		QualityType getOverallQuality() const;

		/// Set the overall quality
		void setOverallQuality(QualityType q);

		/// Non-mutable access to the quality in dimension c
		QualityType getQuality( Size index ) const;
		/// Set the quality in dimension c
		void setQuality( Size index, QualityType q );

		/// Compare by quality
		typedef QualityLess OverallQualityLess;

		/// Non-mutable access to the model description
		const ModelDescription<2>& getModelDescription() const;

		/// Mutable access to the model description
		ModelDescription<2>& getModelDescription();

		/// Set the model description
		void setModelDescription( const ModelDescription<2>& q );
		//@}

		///@name Convex hulls and bounding box
		//@{
		/// Non-mutable access to the convex hulls
		const std::vector<ConvexHull2D>& getConvexHulls() const;
		/// Mutable access to the convex hulls of single mass traces
		std::vector<ConvexHull2D>& getConvexHulls();
		/// Set the convex hulls of single mass traces
		void setConvexHulls( const std::vector<ConvexHull2D>& hulls );
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

		/// immutable access to subordinate features
		const std::vector<Feature>& getSubordinates() const;

		/// mutable access to subordinate features
		std::vector<Feature>& getSubordinates();

		/// mutable access to subordinate features
		void setSubordinates(const std::vector<Feature>& rhs);

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

		/// Quality measures for each dimension
		QualityType qualities_[2];

		/// Description of the theoretical model the feature was constructed with
		ModelDescription<2> model_desc_;

		/// Array of convex hulls (one for each mass trace)
		std::vector<ConvexHull2D> convex_hulls_;

		/// Flag that indicates if the overall convex hull needs to be recomputed (i.e. mass trace convex hulls were modified)
		mutable bool convex_hulls_modified_;

		/// Overall convex hull of the feature
		mutable ConvexHull2D convex_hull_;

		/// subordinate features (e.g. features that the ModelFitter discarded due to inferior quality)
		std::vector<Feature> subordinates_;
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_FEATURE_H
