// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MRM/TransitionInterpretation.h>

#include <algorithm>

using namespace std;


namespace OpenMS
{
	TransitionInterpretation::TransitionInterpretation()
		:	CVTermList(),
			mz_delta_(numeric_limits<DoubleReal>::max()),
			is_primary_(false),
			product_ordinal_(numeric_limits<Size>::max())
	{
	}

  TransitionInterpretation::TransitionInterpretation(const TransitionInterpretation& rhs)
		:	CVTermList(rhs),
			mz_delta_(rhs.mz_delta_),
      is_primary_(rhs.is_primary_),
      product_adjustment_(rhs.product_adjustment_),
      product_ordinal_(rhs.product_ordinal_),
      product_series_(rhs.product_series_)
	{
	}

	TransitionInterpretation::~TransitionInterpretation()
	{
	}

	TransitionInterpretation& TransitionInterpretation::operator = (const TransitionInterpretation& rhs)
	{
		if (&rhs != this)
		{
			mz_delta_ = rhs.mz_delta_;
			is_primary_ = rhs.is_primary_;
			product_adjustment_ = rhs.product_adjustment_;
			product_ordinal_ = rhs.product_ordinal_;
			product_series_ = rhs.product_series_;
		}
		return *this;
	}

	void TransitionInterpretation::setMZDelta(DoubleReal delta)
	{
		mz_delta_ = delta;
	}

	DoubleReal TransitionInterpretation::getMZDelta() const
	{
		return mz_delta_;
	}

	void TransitionInterpretation::setPrimary(bool primary)
	{
		is_primary_ = primary;
	}

	bool TransitionInterpretation::getPrimary() const
	{
		return is_primary_;
	}

	void TransitionInterpretation::setProductAdjustment(const String& adjustment)
	{
		product_adjustment_ = adjustment;
	}

	const String& TransitionInterpretation::getProductAdjustment() const
	{
		return product_adjustment_;
	}

	void TransitionInterpretation::setProductOrdinal(Size ordinal)
	{
		product_ordinal_ = ordinal;
	}

	Size TransitionInterpretation::getProductOrdinal() const
	{
		return product_ordinal_;
	}

	void TransitionInterpretation::setProductSeries(const String& series)
	{
		product_series_ = series;
	}

	const String& TransitionInterpretation::getProductSeries() const
	{
		return product_series_;
	}
}


