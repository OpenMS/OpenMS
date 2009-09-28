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

#include <OpenMS/ANALYSIS/MRM/TransitionPrediction.h>

#include <algorithm>

using namespace std;


namespace OpenMS
{
	TransitionPrediction::TransitionPrediction()
		:	MetaInfoInterface(),
			CVTermList(),
			intensity_rank_(0),
			recommended_transition_rank_(0),
			relative_intensity_(0),
			transition_source_("")
	{
	}

  TransitionPrediction::TransitionPrediction(const TransitionPrediction& rhs)
		:	MetaInfoInterface(rhs),
			CVTermList(rhs),
			intensity_rank_(rhs.intensity_rank_),
      recommended_transition_rank_(rhs.recommended_transition_rank_),
      relative_intensity_(rhs.relative_intensity_),
      transition_source_(rhs.transition_source_)
	{
	}

	TransitionPrediction::~TransitionPrediction()
	{
	}

	TransitionPrediction& TransitionPrediction::operator = (const TransitionPrediction& rhs)
	{
		if (&rhs != this)
		{
			MetaInfoInterface::operator = (rhs);
			CVTermList::operator = (rhs);
			intensity_rank_ = rhs.intensity_rank_;
			recommended_transition_rank_ = rhs.recommended_transition_rank_;
			relative_intensity_ = rhs.relative_intensity_;
			transition_source_ = rhs.transition_source_;
		}
		return *this;
	}

	void TransitionPrediction::setIntensityRank(Size rank)
	{
		intensity_rank_ = rank;
	}

	Size TransitionPrediction::getIntensityRank() const
	{
		return intensity_rank_;
	}

	void TransitionPrediction::setRecommendedTransitionRank(Size rank)
	{
		recommended_transition_rank_ = rank;
	}

	Size TransitionPrediction::getRecommendedTransitionRank() const
	{
		return recommended_transition_rank_;
	}

	void TransitionPrediction::setRelativeIntensity(DoubleReal intensity)
	{
		relative_intensity_ = intensity;
	}

	DoubleReal TransitionPrediction::getRelativeIntensity() const
	{
		return relative_intensity_;
	}

	void TransitionPrediction::setTransitionSource(const String& transition_source)
	{
		transition_source_ = transition_source;
	}	

	const String& TransitionPrediction::getTransitionSource() const
	{
		return transition_source_;
	}	

}


