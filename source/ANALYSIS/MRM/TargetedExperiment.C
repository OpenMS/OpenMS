// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/MRM/TargetedExperiment.h>
#include <algorithm>

using namespace std;


namespace OpenMS
{
	TargetedExperiment::TargetedExperiment()
	{
	}

  TargetedExperiment::TargetedExperiment(const TargetedExperiment& rhs)
		: cvs_(rhs.cvs_),
			contacts_(rhs.contacts_),
			publications_(rhs.publications_),
			instruments_(rhs.instruments_),
			software_(rhs.software_),
			proteins_(rhs.proteins_),
			compounds_(rhs.compounds_),
			peptides_(rhs.peptides_),
			transitions_(rhs.transitions_)
	{
	}

	TargetedExperiment::~TargetedExperiment()
	{
	}

	TargetedExperiment& TargetedExperiment::operator = (const TargetedExperiment& rhs)
	{
		if (&rhs != this)
		{
			cvs_ = rhs.cvs_;
      contacts_ = rhs.contacts_;
      publications_ = rhs.publications_;
      instruments_ = rhs.instruments_;
      software_ = rhs.software_;
      proteins_ = rhs.proteins_;
      compounds_ = rhs.compounds_;
      peptides_ = rhs.peptides_;
			transitions_ = rhs.transitions_;
		}
		return *this;
	}


	bool TargetedExperiment::operator == (const TargetedExperiment& rhs) const
	{
		return 	cvs_ == rhs.cvs_ &&
      			contacts_ == rhs.contacts_ &&
			      publications_ == rhs.publications_ &&
			      instruments_ == rhs.instruments_ &&
      			software_ == rhs.software_ &&
      			proteins_ == rhs.proteins_ &&
      			compounds_ == rhs.compounds_ &&
      			peptides_ == rhs.peptides_ &&
      			transitions_ == rhs.transitions_;
	}

	void TargetedExperiment::setCVs(const vector<CV>& cvs)
	{
		cvs_ = cvs;
	}

	const vector<TargetedExperiment::CV>& TargetedExperiment::getCVs() const
	{
		return cvs_;
	}

	void TargetedExperiment::addCV(const CV& cv)
	{
		cvs_.push_back(cv);
	}

  void TargetedExperiment::setContacts(const vector<CVTermList>& contacts)
  {
    contacts_ = contacts;
  }

  const vector<CVTermList>& TargetedExperiment::getContacts() const
  {
    return contacts_;
  }

  void TargetedExperiment::addContact(const CVTermList& contact)
  {
    contacts_.push_back(contact);
  }

	void TargetedExperiment::setPublications(const vector<CVTermList>& publications)
	{
		publications_ = publications;
	}

	const vector<CVTermList>& TargetedExperiment::getPublications() const
	{
		return publications_;
	}

	void TargetedExperiment::addPublication(const CVTermList& publication)
	{
		publications_.push_back(publication);
	}

  void TargetedExperiment::setInstruments(const vector<CVTermList>& instruments)
  {
    instruments_ = instruments;
  }

  const vector<CVTermList>& TargetedExperiment::getInstruments() const
  {
    return instruments_;
  }

  void TargetedExperiment::addInstrument(const CVTermList& instrument)
  {
    instruments_.push_back(instrument);
  }

	void TargetedExperiment::setSoftware(const vector<Software>& software)
  {
    software_ = software;
  }

  const vector<Software>& TargetedExperiment::getSoftware() const
  {
    return software_;
  }

  void TargetedExperiment::addSoftware(const Software& software)
  {
    software_.push_back(software);
  }

  void TargetedExperiment::setProteins(const vector<Protein>& proteins)
  {
    proteins_ = proteins;
  }

  const vector<TargetedExperiment::Protein>& TargetedExperiment::getProteins() const
  {
    return proteins_;
  }

  void TargetedExperiment::addProtein(const Protein& protein)
  {
    proteins_.push_back(protein);
  }

	void TargetedExperiment::setCompounds(const vector<Compound>& compounds)
	{
		compounds_ = compounds;
	}

	const vector<TargetedExperiment::Compound>& TargetedExperiment::getCompounds() const
	{
		return compounds_;
	}

	void TargetedExperiment::addCompound(const Compound& rhs)
	{
		compounds_.push_back(rhs);
	}

  void TargetedExperiment::setPeptides(const vector<Peptide>& peptides)
  {
    peptides_ = peptides;
  }

  const vector<TargetedExperiment::Peptide>& TargetedExperiment::getPeptides() const
  {
    return peptides_;
  }

  void TargetedExperiment::addPeptide(const Peptide& rhs)
  {
    peptides_.push_back(rhs);
  }

	void TargetedExperiment::setTransitions(const vector<ReactionMonitoringTransition>& transitions)
	{
		transitions_ = transitions;
	}

	const vector<ReactionMonitoringTransition>& TargetedExperiment::getTransitions() const
	{
		return transitions_;
	}

	void TargetedExperiment::addTransition(const ReactionMonitoringTransition& transition)
	{
		transitions_.push_back(transition);
	}

	void TargetedExperiment::setIncludeTargets(const vector<IncludeExcludeTarget>& targets)
	{
		include_targets_ = targets;
	}

	const vector<IncludeExcludeTarget>& TargetedExperiment::getIncludeTargets() const
	{
		return include_targets_;
	}

	void TargetedExperiment::addIncludeTarget(const IncludeExcludeTarget& target)
	{
		include_targets_.push_back(target);
	}

  void TargetedExperiment::setExcludeTargets(const vector<IncludeExcludeTarget>& targets)
  {
    exclude_targets_ = targets;
  }

  const vector<IncludeExcludeTarget>& TargetedExperiment::getExcludeTargets() const
  {
    return exclude_targets_;
  }

  void TargetedExperiment::addExcludeTarget(const IncludeExcludeTarget& target)
  {
    exclude_targets_.push_back(target);
  }


	void TargetedExperiment::setSourceFiles(const vector<SourceFile>& source_files)
	{
		source_files_ = source_files;
	}

	const vector<SourceFile>& TargetedExperiment::getSourceFiles() const
	{
		return source_files_;
	}

	void TargetedExperiment::addSourceFile(const SourceFile& source_file)
	{
		source_files_.push_back(source_file);
	}
}


