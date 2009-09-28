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
			CV(const String& new_id, const String& new_fullname, const String& new_version, const String& new_URI)
				:	id(new_id),
					fullname(new_fullname),
					version(new_version),
					URI(new_URI)
			{

			}
			String id;
			String fullname;
			String version;
			String URI;
			
			bool operator == (const CV& cv) const
			{
				return 	id == cv.id &&
								fullname == cv.fullname &&
								version == cv.version &&
								URI == cv.URI;
			}

		};

		struct Protein
			: public CVTermList
		{
			Protein()
				: CVTermList()
			{
			}

			String id;
			String accession;
			String name;
			String description;
			String comment;
			String sequence;

			bool operator == (const Protein& rhs) const
			{
				return  CVTermList::operator == (rhs) &&
								id == rhs.id &&
								accession == rhs.accession &&
								name == rhs.name &&
								description == rhs.description &&
								comment == rhs.comment && 
								sequence == rhs.sequence;
			}

      Protein& operator = (const Protein& rhs)
      {
				if (&rhs != this)
				{
					CVTermList::operator = (rhs);
        	id = rhs.id;
          accession = rhs.accession;
          name = rhs.name;
          description = rhs.description;
          comment = rhs.comment;
          sequence = rhs.sequence;
				}
				return *this;
      }

		};

		class OPENMS_DLLAPI RetentionTime
			: public CVTermList
		{
			public: 

			RetentionTime()
				: CVTermList(),
					local_retention_time(0),
					normalized_retention_time(0),
					predicted_retention_time(0)
			{
			}

			RetentionTime(const RetentionTime& rhs)
        : CVTermList(rhs),
					local_retention_time(rhs.local_retention_time),
					normalization_standard(rhs.normalization_standard),
					normalized_retention_time(rhs.normalized_retention_time),
					predicted_retention_time(rhs.predicted_retention_time),
					predicted_retention_time_software_ref(rhs.predicted_retention_time_software_ref)
      {
      }

			virtual ~RetentionTime()
			{
			}

			RetentionTime& operator = (const RetentionTime& rhs)
			{
				if (&rhs != this)
				{
					CVTermList::operator = (rhs);
					local_retention_time = rhs.local_retention_time;
					normalization_standard = rhs.normalization_standard;
					normalized_retention_time = rhs.normalized_retention_time;
					predicted_retention_time = rhs.predicted_retention_time;
					predicted_retention_time_software_ref = rhs.predicted_retention_time_software_ref;
				}
				return *this;
			}

      bool operator == (const RetentionTime& rhs) const
      {
				return	CVTermList::operator == (rhs) && 
								local_retention_time == rhs.local_retention_time &&
			          normalization_standard == rhs.normalization_standard &&
          			normalized_retention_time == rhs.normalized_retention_time &&
          			predicted_retention_time == rhs.predicted_retention_time &&
          			predicted_retention_time_software_ref == rhs.predicted_retention_time_software_ref;
      }


			DoubleReal local_retention_time;
			String normalization_standard;
			DoubleReal normalized_retention_time;
			DoubleReal predicted_retention_time;
			String predicted_retention_time_software_ref;
		};

		class OPENMS_DLLAPI Compound
			: public CVTermList
		{
			public:
				
			Compound()
				: CVTermList()
			{
			}

			Compound(const Compound& rhs)
				:	CVTermList(rhs),
					id(rhs.id),
					rts(rhs.rts)
			{
			}
			
			Compound& operator = (const Compound& rhs)
			{
				if (this != &rhs)
				{
					CVTermList::operator = (rhs);
					id = rhs.id;
					rts = rhs.rts;
				}
				return *this;
			}

      bool operator == (const Compound& rhs) const
      {
				return	CVTermList::operator == (rhs) && 
								id == rhs.id &&
			         	rts == rhs.rts;
      }

			String id;
			std::vector<RetentionTime> rts;
		};
		

    class OPENMS_DLLAPI Peptide
			: public CVTermList
    {
      public:

      Peptide()
				: CVTermList()
      {
      }

      Peptide(const Peptide& rhs)
        : CVTermList(rhs),
					rts(rhs.rts),
					id(rhs.id),
					group_label(rhs.group_label),
					labeling_category(rhs.labeling_category),
					modified_sequence(rhs.modified_sequence),
					unmodified_sequence(rhs.unmodified_sequence),
					protein_ref(rhs.protein_ref),
					evidence(rhs.evidence)
      {
      }

      Peptide& operator = (const Peptide& rhs)
      {
        if (this != &rhs)
        {
					CVTermList::operator = (rhs);
          rts = rhs.rts;
					id = rhs.id;
					group_label = rhs.group_label;
					labeling_category = rhs.labeling_category;
					modified_sequence = rhs.modified_sequence;
					unmodified_sequence = rhs.unmodified_sequence;
					protein_ref = rhs.protein_ref;
					evidence = rhs.evidence;
        }
        return *this;
      }

      bool operator == (const Peptide& rhs) const
      {
				return	CVTermList::operator == (rhs) && 
								rts == rhs.rts &&
          			id == rhs.id &&
          			group_label == rhs.group_label &&
          			labeling_category == rhs.labeling_category &&
          			modified_sequence == rhs.modified_sequence &&
          			unmodified_sequence == rhs.unmodified_sequence &&
          			protein_ref == rhs.protein_ref &&
          			evidence == rhs.evidence;
      }



      std::vector<RetentionTime> rts;
			String id;
			String group_label;
			String labeling_category;
			String modified_sequence;
			String unmodified_sequence;
			String protein_ref;
			CVTermList evidence;
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
		void setContacts(const std::vector<CVTermList>& contacts);

		const std::vector<CVTermList>& getContacts() const;

		void addContact(const CVTermList& contact);

		// publication list
    void setPublications(const std::vector<CVTermList>& publications);

    const std::vector<CVTermList>& getPublications() const;

    void addPublication(const CVTermList& publication);

		// instrument list
		void setInstruments(const std::vector<CVTermList>& instruments);

		const std::vector<CVTermList>& getInstruments() const;

		void addInstrument(const CVTermList& instrument);

		// software list
		void setSoftware(const std::vector<Software>& software);

		const std::vector<Software>& getSoftware() const;

		void addSoftware(const Software& software);

		// protein list
	  void setProteins(const std::vector<Protein>& proteins);

    const std::vector<Protein>& getProteins() const;

    void addProtein(const Protein& protein);

		// compound list
		void setCompounds(const std::vector<Compound>& rhs);

		const std::vector<Compound>& getCompounds() const;

		void addCompound(const Compound& rhs);

    void setPeptides(const std::vector<Peptide>& rhs);

    const std::vector<Peptide>& getPeptides() const;

    void addPeptide(const Peptide& rhs);

		/// set transition list
		void setTransitions(const std::vector<ReactionMonitoringTransition>& transitions);
		
		/// returns the transition list
		const std::vector<ReactionMonitoringTransition>& getTransitions() const;

		/// adds a transition to the list
		void addTransition(const ReactionMonitoringTransition& transition);
		//@}

		protected:

		std::vector<CV> cvs_;

		std::vector<CVTermList> contacts_;

		std::vector<CVTermList> publications_;

		std::vector<CVTermList> instruments_;

		std::vector<Software> software_;

		std::vector<Protein> proteins_;
		
		std::vector<Compound> compounds_;

		std::vector<Peptide> peptides_;

		std::vector<ReactionMonitoringTransition> transitions_;

	};
}

#endif // OPENMS_ANALYSIS_MRM_MRMEXPERIMENT_H

