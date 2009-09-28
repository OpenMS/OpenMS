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

#ifndef OPENMS_ANALYSIS_MRM_TRANSITIONINTERPRETATION_H
#define OPENMS_ANALYSIS_MRM_TRANSITIONINTERPRETATION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
	/**
		@brief This class stores an interpretation of an SRM/MRM transition

		This class is intended to be used to store intepretations of MRM
		transitions, in the class ReactionMonitoringTransition. It can store
		the type of ion and different other values.

	*/
	class OPENMS_DLLAPI TransitionInterpretation
		: public CVTermList
	{

		public:

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		TransitionInterpretation();

		/// copy constructor
		TransitionInterpretation(const TransitionInterpretation& rhs);

		/// destructor
		virtual ~TransitionInterpretation();
		//@}

		/// assignment operator 
		TransitionInterpretation& operator = (const TransitionInterpretation& rhs);

		/** @name Accessors
		*/
		//@{
		void setMZDelta(DoubleReal mz_delta);

		DoubleReal getMZDelta() const;

		void setPrimary(bool primary);

		bool getPrimary() const;

		void setProductAdjustment(const String& adjustment);

		const String& getProductAdjustment() const;

		void setProductOrdinal(Size ordinal);
	
		Size getProductOrdinal() const;

		void setProductSeries(const String& series);

		const String& getProductSeries() const;
		//@}

		protected:

		DoubleReal mz_delta_;
		
		bool is_primary_;
		
		String product_adjustment_;
		
		Size product_ordinal_;
		
		String product_series_;
	};
}

#endif // OPENMS_ANALYSIS_MRM_TRANSITIONINTERPRETATION_H

