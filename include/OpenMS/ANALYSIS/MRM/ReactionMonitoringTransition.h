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

#ifndef OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H
#define OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/ANALYSIS/MRM/TransitionInterpretation.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class stores a SRM/MRM transition

		bla
	*/
	class OPENMS_DLLAPI ReactionMonitoringTransition : public MetaInfoInterface
	{

		public:

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		ReactionMonitoringTransition();

		/// copy constructor
		ReactionMonitoringTransition(const ReactionMonitoringTransition& rhs);

		/// destructor
		virtual ~ReactionMonitoringTransition();
		//@}

		/// assignment operator 
		ReactionMonitoringTransition& operator = (const ReactionMonitoringTransition& rhs);

		/** @name Accessors
		*/
		//@{
		/// sets the precursor mz (Q1 value)
		void setPrecursorMZ(DoubleReal mz);

		DoubleReal getPrecursorMZ() const;

		void setPrecursorCharge(Int mz);

		Int getPrecursorCharge() const;

		void setProductMZ(DoubleReal mz);

		DoubleReal getProductMZ() const;

		void setProductCharge(Int mz);

		Int getProductCharge() const;

		void setInterpretationList(const std::vector<TransitionInterpretation>& interpretations);

		const std::vector<TransitionInterpretation>& getInterpretationList() const;

		void addInterpretation(const TransitionInterpretation& interpretation);
		//@}

		protected:

		void updateMembers_();

		DoubleReal precursor_mz_;

		Int precursor_charge_;

		DoubleReal product_mz_;

		Int product_charge_;

		std::vector<TransitionInterpretation> interpretation_list_;
	
		//vector<Configurations> configurations_list_; why multiple configuration per transition????
	};
}

#endif // OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H

