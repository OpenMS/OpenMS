// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H
#define OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/Software.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class stores an prediction of an SRM/MRM transition


	*/
	class OPENMS_DLLAPI TargetedExperiment
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
			String sequence;

			bool operator == (const Protein& rhs) const
			{
				return  CVTermList::operator == (rhs) &&
								id == rhs.id &&
								sequence == rhs.sequence;
			}

      Protein& operator = (const Protein& rhs)
      {
				if (&rhs != this)
				{
					CVTermList::operator = (rhs);
        	id = rhs.id;
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
				: CVTermList()
			{
			}

			RetentionTime(const RetentionTime& rhs)
        : CVTermList(rhs),
					software_ref(rhs.software_ref)
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
					software_ref = rhs.software_ref;
				}
				return *this;
			}

      bool operator == (const RetentionTime& rhs) const
      {
				return	CVTermList::operator == (rhs) && 
          			software_ref == rhs.software_ref;
      }


			String software_ref;
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

			struct Modification
				: public CVTermList
			{
				DoubleReal avg_mass_delta;
				Size location;
				DoubleReal mono_mass_delta;
			};

      Peptide()
				: CVTermList()
      {
      }

      Peptide(const Peptide& rhs)
        : CVTermList(rhs),
					rts(rhs.rts),
					id(rhs.id),
					protein_refs(rhs.protein_refs),
					evidence(rhs.evidence),
					sequence(rhs.sequence),
					mods(rhs.mods)
      {
      }

      Peptide& operator = (const Peptide& rhs)
      {
        if (this != &rhs)
        {
					CVTermList::operator = (rhs);
          rts = rhs.rts;
					id = rhs.id;
					protein_refs = rhs.protein_refs;
					evidence = rhs.evidence;
					sequence = rhs.sequence;
					mods = rhs.mods;
        }
        return *this;
      }

      bool operator == (const Peptide& rhs) const
      {
				return	CVTermList::operator == (rhs) && 
								rts == rhs.rts &&
          			id == rhs.id &&
          			protein_refs == rhs.protein_refs &&
          			evidence == rhs.evidence &&
								sequence == rhs.sequence &&
								mods == rhs.mods;
      }



      std::vector<RetentionTime> rts;
			String id;
			std::vector<String> protein_refs;
			CVTermList evidence;
			String sequence;
			std::vector<Modification> mods;
    };


		/** @name Constructors and destructors
		*/
		//@{
		/// default constructor
		TargetedExperiment();

		/// copy constructor
		TargetedExperiment(const TargetedExperiment& rhs);

		/// destructor
		virtual ~TargetedExperiment();
		//@}

		/// assignment operator 
		TargetedExperiment& operator = (const TargetedExperiment& rhs);

		/** @name Predicates
		*/
		//@{
		bool operator == (const TargetedExperiment& rhs) const;
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

		void setIncludeTargets(const std::vector<IncludeExcludeTarget>& targets);

		const std::vector<IncludeExcludeTarget>& getIncludeTargets() const;

		void addIncludeTarget(const IncludeExcludeTarget& target);

		void setExcludeTargets(const std::vector<IncludeExcludeTarget>& targets);

		const std::vector<IncludeExcludeTarget>& getExcludeTargets() const;
	
		void addExcludeTarget(const IncludeExcludeTarget& target);

		/// sets the source files
		void setSourceFiles(const std::vector<SourceFile>& source_files);

		/// returns the source file list
		const std::vector<SourceFile>& getSourceFiles() const;

		/// adds a source file to the list
		void addSourceFile(const SourceFile& source_file);
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

		std::vector<IncludeExcludeTarget> include_targets_;

		std::vector<IncludeExcludeTarget> exclude_targets_;

		std::vector<SourceFile> source_files_;

	};
}

#endif // OPENMS_ANALYSIS_TARGETED_TARGETEDEXPERIMENT_H

