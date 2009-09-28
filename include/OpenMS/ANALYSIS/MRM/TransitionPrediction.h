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

#ifndef OPENMS_ANALYSIS_MRM_TRANSITIONPREDICTION_H
#define OPENMS_ANALYSIS_MRM_TRANSITIONPREDICTION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
	/**
		@brief This class stores an prediction of an SRM/MRM transition


	*/
	class OPENMS_DLLAPI TransitionPrediction 
		: public MetaInfoInterface,
			public CVTermList
	{

		public:

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		TransitionPrediction();

		/// copy constructor
		TransitionPrediction(const TransitionPrediction& rhs);

		/// destructor
		virtual ~TransitionPrediction();
		//@}

		/// assignment operator 
		TransitionPrediction& operator = (const TransitionPrediction& rhs);

		/** @name Accessors
		*/
		//@{
		void setIntensityRank(Size rank);

		Size getIntensityRank() const;

		void setRecommendedTransitionRank(Size rank);

		Size getRecommendedTransitionRank() const;

		void setRelativeIntensity(DoubleReal intensity);

		DoubleReal getRelativeIntensity() const;
		
		void setTransitionSource(const String& transition_source);

		const String& getTransitionSource() const;
		//@}

		protected:

		Size intensity_rank_;
		Size recommended_transition_rank_;
		DoubleReal relative_intensity_;
		String transition_source_;
	};
}

#endif // OPENMS_ANALYSIS_MRM_TRANSITIONINTERPRETATION_H

