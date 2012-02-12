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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_BASEFEATURE_H
#define OPENMS_KERNEL_BASEFEATURE_H

#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{

	/** @brief A basic LC-MS feature.

	This class represents a "minimal" feature, defined by a position in RT and m/z, intensity, charge, quality, and annotated peptides. Most classes dealing with features will use the subclasses Feature or ConsensusFeature directly. However, algorithms that rely on very general characteristics of features can use this class to provide a unified solution for both "normal" features and consensus features.

	@ingroup Kernel
	*/
	class OPENMS_DLLAPI BaseFeature
		: public RichPeak2D
	{
	 public:
		///@name Type definitions
		//@{
		/// Type of quality values
		typedef		Real QualityType;
		/// Type of charge values
		typedef		Int ChargeType;
		/// Type of feature width/FWHM (RT)
		typedef		Real WidthType;
		//@}

		/** @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		BaseFeature();

		/// Copy constructor
		BaseFeature(const BaseFeature& feature);

		///Constructor from raw data point
		BaseFeature(const Peak2D& point);

		/// Constructor from raw data point with meta information
		BaseFeature(const RichPeak2D& point);

		/// Destructor
		~BaseFeature();
		//@}

		/// @name Quality methods
		//@{
		/// Non-mutable access to the overall quality
		QualityType getQuality() const;
		/// Set the overall quality
		void setQuality(QualityType q);
		/// Compare by quality
		struct QualityLess: std::binary_function<BaseFeature, BaseFeature, bool>
		{
			bool operator()(const BaseFeature& left, const BaseFeature& right) const
			{
				return (left.getQuality() < right.getQuality());
			}
			bool operator()(const BaseFeature& left, const QualityType& right) const
			{
				return (left.getQuality() < right);
			}
			bool operator()(const QualityType& left, const BaseFeature& right) const
			{
				return (left < right.getQuality());
			}
			bool operator()(const QualityType& left, const QualityType& right) const
			{
				return (left < right);
			}
		};
		//@}

		/// Non-mutable access to the features width (full width at half max, FWHM)
		WidthType getWidth() const;
		/// Set the width of the feature (FWHM)
		void setWidth(WidthType fwhm);

		/// Non-mutable access to charge state
		const ChargeType& getCharge() const;

		/// Set charge state
		void setCharge(const ChargeType& ch);

		/// Assignment operator
		BaseFeature& operator=(const BaseFeature& rhs);

		/// Equality operator
		bool operator==(const BaseFeature& rhs) const;

		/// Unequality operator
		bool operator!=(const BaseFeature& rhs) const;

		/// returns a const reference to the PeptideIdentification vector
		const std::vector<PeptideIdentification>& getPeptideIdentifications() const;

		/// returns a mutable reference to the PeptideIdentification vector
		std::vector<PeptideIdentification>& getPeptideIdentifications();

		/// sets the PeptideIdentification vector
		void setPeptideIdentifications(const std::vector<PeptideIdentification>& peptides);

	 protected:

		/// Overall quality measure of the feature
		QualityType quality_;

		/// Charge of the peptide represented by this feature.  The default value is 0, which represents an unknown charge state.
		ChargeType charge_;

		/// Width (FWHM) for the feature. The default value is 0.0, a feature finding algorithm can compute this form the model.
		WidthType width_;

		/// Peptide PeptideIdentifications belonging to the feature
		std::vector<PeptideIdentification> peptides_;
	};

} // namespace OpenMS

#endif // OPENMS_KERNEL_BASEFEATURE_H
