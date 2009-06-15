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

		The default values for precursor and product m/z values are
		set to numeric_limits<DoubleReal>::max(). Default values for 
		precursor an product charge is set to numeric_limits<Int>::max().
	*/
	class OPENMS_DLLAPI ReactionMonitoringTransition : public MetaInfoInterface
	{

		public:

		struct Validation
		{
			String transition_source;
			DoubleReal relative_intensity;
			Size recommended_transition_rank;
			Size intensity_rank;
			std::vector<MetaInfoInterface> cvs;

			Validation& operator = (const Validation& rhs)
			{
				if (this != &rhs)
				{
					transition_source = rhs.transition_source;
					relative_intensity = rhs.relative_intensity;
					recommended_transition_rank = rhs.recommended_transition_rank;
					intensity_rank = rhs.intensity_rank;
					cvs = rhs.cvs;
				}
				return *this;
			}
		};

		struct Configuration
		{
			String contact_ref;
			String instrument_ref;
			std::vector<Validation> validations;
			std::vector<MetaInfoInterface> cvs;

			Configuration& operator = (const Configuration& rhs)
			{
				if (this != &rhs)
				{
					contact_ref = rhs.contact_ref;
					instrument_ref = rhs.instrument_ref;
					validations = rhs.validations;
					cvs = rhs.cvs;
				}
				return *this;
			}
		};

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
		void setName(const String& name);

		const String& getName() const;

		void setPeptideRef(const String& peptide_ref);

		const String& getPeptideRef() const;

		void setCompoundRef(const String& compound_ref);

		const String& getCompoundRef() const;

		/// sets the precursor mz (Q1 value)
		void setPrecursorMZ(DoubleReal mz);

		DoubleReal getPrecursorMZ() const;

		void setPrecursorCharge(Int mz);

		Int getPrecursorCharge() const;

		void setProductMZ(DoubleReal mz);

		DoubleReal getProductMZ() const;

		void setProductCharge(Int mz);

		Int getProductCharge() const;

		void setInterpretations(const std::vector<TransitionInterpretation>& interpretations);

		const std::vector<TransitionInterpretation>& getInterpretations() const;

		void addInterpretation(const TransitionInterpretation& interpretation);

		void setConfigurations(const std::vector<Configuration>& configuration);
		
		const std::vector<Configuration>& getConfigurations() const;

		void addConfiguration(const Configuration& configuration);
		//@}

		protected:

		void updateMembers_();

		String name_;

		DoubleReal precursor_mz_;

		Int precursor_charge_;

		DoubleReal product_mz_;

		Int product_charge_;

		std::vector<TransitionInterpretation> interpretation_list_;
	
		//vector<Configuration> configurations_list_; why multiple configuration per transition????

		String peptide_ref_;

		String compound_ref_;

		std::vector<Configuration> configurations_;
	};
}

#endif // OPENMS_ANALYSIS_MRM_REACTIONMONITORINGTRANSITION_H

