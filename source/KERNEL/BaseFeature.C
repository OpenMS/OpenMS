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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/BaseFeature.h>

using namespace std;

namespace OpenMS
{
	BaseFeature::BaseFeature()
		: RichPeak2D(), quality_(0.0), charge_(0), peptides_()
	{
	}

	BaseFeature::BaseFeature(const BaseFeature& rhs)
		: RichPeak2D(rhs), quality_(rhs.quality_), charge_(rhs.charge_),
			peptides_(rhs.peptides_)
	{
	}

	BaseFeature::BaseFeature(const RichPeak2D& point)
		: RichPeak2D(point), quality_(0.0), charge_(0), peptides_()
	{
	}

	BaseFeature::BaseFeature(const Peak2D& point)
		: RichPeak2D(point), quality_(0.0), charge_(0), peptides_()
	{
	}

	BaseFeature& BaseFeature::operator=(const BaseFeature& rhs)
	{
		if (&rhs == this) return *this;

		RichPeak2D::operator=(rhs);
		quality_ = rhs.quality_;
		charge_ = rhs.charge_;
		peptides_ =  rhs.peptides_;

		return *this;
	}

	bool BaseFeature::operator==(const BaseFeature& rhs) const
	{
		return (RichPeak2D::operator==(rhs) && (quality_ == rhs.quality_) &&
						(charge_ == rhs.charge_) && (peptides_ == rhs.peptides_));
	}

	bool BaseFeature::operator!=(const BaseFeature& rhs) const
	{
		return !operator==(rhs);
	}

	BaseFeature::~BaseFeature()
	{
	}
		
	const BaseFeature::QualityType& BaseFeature::getQuality() const
	{
		return quality_;
	}

	void BaseFeature::setQuality(const BaseFeature::QualityType& quality)
	{
		quality_ = quality;
	}
	
	const BaseFeature::ChargeType& BaseFeature::getCharge() const
	{
		return charge_;
	}
	
	void BaseFeature::setCharge(const BaseFeature::ChargeType& charge)
	{
		charge_ = charge;
	}
	
	const vector<PeptideIdentification>& BaseFeature::getPeptideIdentifications()
		const
	{
		return peptides_;
	}

	vector<PeptideIdentification>& BaseFeature::getPeptideIdentifications()
	{
		return peptides_;
	}

	void BaseFeature::setPeptideIdentifications(
		const vector<PeptideIdentification>& peptides)
	{
		peptides_ = peptides;
	}

} // namespace OpenMS
