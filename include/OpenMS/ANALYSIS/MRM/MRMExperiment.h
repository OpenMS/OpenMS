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

#ifndef OPENMS_ANALYSIS_MRM_MRMEXPERIMENT_H
#define OPENMS_ANALYSIS_MRM_MRMEXPERIMENT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

namespace OpenMS
{
	/**
		@brief This class stores an prediction of an SRM/MRM transition


	*/
	class MRMExperiment
	{

		public:

		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		MRMExperiment();

		/// copy constructor
		MRMExperiment(const MRMExperiment& rhs);

		/// destructor
		virtual ~MRMExperiment();
		//@}

		/// assignment operator 
		MRMExperiment& operator = (const MRMExperiment& rhs);

		/** @name Accessors
		*/
		//@{
		// cv list

		// publications list

		// instrument list

		// software list

		// protein list

		// compound list

		void setTransitions(const std::vector<ReactionMonitoringTransition>& transitions);
		
		const std::vector<ReactionMonitoringTransition>& getTransitions() const;

		void addTransition(const ReactionMonitoringTransition& transition);
		//@}

		protected:

		std::vector<ReactionMonitoringTransition> transitions_;
	};
}

#endif // OPENMS_ANALYSIS_MRM_MRMEXPERIMENT_H

