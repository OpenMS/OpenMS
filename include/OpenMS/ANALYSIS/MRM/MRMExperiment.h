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
#include <OpenMS/METADATA/Software.h>

namespace OpenMS
{
	/**
		@brief This class stores an prediction of an SRM/MRM transition


	*/
	class OPENMS_DLLAPI MRMExperiment
	{
		public:

		struct CV
		{
			String id;
			String fullname;
			String version;
			String URI;
		};

		struct Protein
		{
			String id;
			String accession;
			String name;
			String description;
			String comment;
			String sequence;
		};

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

		/** @name Predicates
		*/
		//@{
		bool operator == (const MRMExperiment& rhs) const;
		//@}

		/** @name Accessors
		*/
		//@{
		// cv list
		void setCVs(const std::vector<CV>& cvs);
	
		const std::vector<CV>& getCVs() const;

		void addCV(const CV& cv);

		// contact list
		void setContacts(const std::vector<MetaInfoInterface>& contacts);

		const std::vector<MetaInfoInterface>& getContacts() const;

		void addContact(const MetaInfoInterface& contact);

		// publication list
    void setPublications(const std::vector<MetaInfoInterface>& publications);

    const std::vector<MetaInfoInterface>& getPublications() const;

    void addPublication(const MetaInfoInterface& publication);

		// instrument list
		void setInstruments(const std::vector<MetaInfoInterface>& instruments);

		const std::vector<MetaInfoInterface>& getInstruments() const;

		void addInstrument(const MetaInfoInterface& instrument);

		// software list
		void setSoftware(const std::vector<Software>& software);

		const std::vector<Software>& getSoftware() const;

		void addSoftware(const Software& software);

		// protein list
	  void setProteins(const std::vector<Protein>& proteins);

    const std::vector<Protein>& getProteins() const;

    void addProtein(const Protein& protein);

		// compound list

		/// set transition list
		void setTransitions(const std::vector<ReactionMonitoringTransition>& transitions);
		
		/// returns the transition list
		const std::vector<ReactionMonitoringTransition>& getTransitions() const;

		/// adds a transition to the list
		void addTransition(const ReactionMonitoringTransition& transition);
		//@}

		protected:

		std::vector<CV> cvs_;

		std::vector<ContactPerson> contacts_;

		std::vector<MetaInfoInterface> publications_;

		std::vector<MetaInfoInterface> instruments_;

		std::vector<Software> software_;

		std::vector<Protein> proteins_;
		//std::vector<> compounds_;

		std::vector<ReactionMonitoringTransition> transitions_;

	};
}

#endif // OPENMS_ANALYSIS_MRM_MRMEXPERIMENT_H

