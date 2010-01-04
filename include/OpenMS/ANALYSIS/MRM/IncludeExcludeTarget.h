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

#ifndef OPENMS_ANALYSIS_MRM_INCLUDEEXCLUDETARGET_H
#define OPENMS_ANALYSIS_MRM_INCLUDEEXCLUDETARGET_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class stores a SRM/MRM transition

		The default values for precursor and product m/z values are
		set to numeric_limits<DoubleReal>::max(). Default values for 
		precursor an product charge is set to numeric_limits<Int>::max().
	*/
	class OPENMS_DLLAPI IncludeExcludeTarget 
		: public CVTermList
	{

		public:

		struct Configuration
			: public CVTermList
		{
			String contact_ref;
			String instrument_ref;
			std::vector<CVTermList> validations;

			Configuration& operator = (const Configuration& rhs)
			{
				if (this != &rhs)
				{
					CVTermList::operator = (rhs);
					contact_ref = rhs.contact_ref;
					instrument_ref = rhs.instrument_ref;
					validations = rhs.validations;
				}
				return *this;
			}
		};

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		IncludeExcludeTarget();

		/// copy constructor
		IncludeExcludeTarget(const IncludeExcludeTarget& rhs);

		/// destructor
		virtual ~IncludeExcludeTarget();
		//@}

		/// assignment operator 
		IncludeExcludeTarget& operator = (const IncludeExcludeTarget& rhs);

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

		void setPrecursorCVTermList(const CVTermList& list);

		void addPrecursorCVTerm(const CVTerm& cv_term);

		const CVTermList& getPrecursorCVTermList() const;
		
		void setProductMZ(DoubleReal mz);

		DoubleReal getProductMZ() const;

		void setProductCVTermList(const CVTermList& list);

		void addProductCVTerm(const CVTerm& cv_term);

		const CVTermList& getProductCVTermList() const;

		void setInterpretations(const std::vector<CVTermList>& interpretations);

		const std::vector<CVTermList>& getInterpretations() const;

		void addInterpretation(const CVTermList& interpretation);

		void setConfigurations(const std::vector<Configuration>& configuration);
		
		const std::vector<Configuration>& getConfigurations() const;

		void addConfiguration(const Configuration& configuration);

		void setPrediction(const CVTermList& prediction);

		void addPredictionTerm(const CVTerm& prediction);

		const CVTermList& getPrediction() const;
		//@}

		/** @name Predicates
		*/
		//@{
		/// equality operator
		bool operator == (const IncludeExcludeTarget& rhs) const;
		
		/// inequality operator
		bool operator != (const IncludeExcludeTarget& rhs) const;
		//@}

		protected:

		void updateMembers_();

		String name_;

		DoubleReal precursor_mz_;

		CVTermList precursor_cv_terms_;

		DoubleReal product_mz_;

		CVTermList product_cv_terms_;

		std::vector<CVTermList> interpretation_list_;
	
		String peptide_ref_;

		String compound_ref_;

		std::vector<Configuration> configurations_;

		CVTermList prediction_;
	};
}

#endif // OPENMS_ANALYSIS_MRM_INCLUDEEXCLUDETARGET_H

