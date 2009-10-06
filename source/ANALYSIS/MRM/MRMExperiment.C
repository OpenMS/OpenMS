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

#include <OpenMS/ANALYSIS/MRM/MRMExperiment.h>
#include <algorithm>

using namespace std;


namespace OpenMS
{
	MRMExperiment::MRMExperiment()
	{
	}

  MRMExperiment::MRMExperiment(const MRMExperiment& rhs)
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

	MRMExperiment::~MRMExperiment()
	{
	}

	MRMExperiment& MRMExperiment::operator = (const MRMExperiment& rhs)
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


	bool MRMExperiment::operator == (const MRMExperiment& rhs) const
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

	void MRMExperiment::setCVs(const vector<CV>& cvs)
	{
		cvs_ = cvs;
	}

	const vector<MRMExperiment::CV>& MRMExperiment::getCVs() const
	{
		return cvs_;
	}

	void MRMExperiment::addCV(const CV& cv)
	{
		cvs_.push_back(cv);
	}

  void MRMExperiment::setContacts(const vector<CVTermList>& contacts)
  {
    contacts_ = contacts;
  }

  const vector<CVTermList>& MRMExperiment::getContacts() const
  {
    return contacts_;
  }

  void MRMExperiment::addContact(const CVTermList& contact)
  {
    contacts_.push_back(contact);
  }

	void MRMExperiment::setPublications(const vector<CVTermList>& publications)
	{
		publications_ = publications;
	}

	const vector<CVTermList>& MRMExperiment::getPublications() const
	{
		return publications_;
	}

	void MRMExperiment::addPublication(const CVTermList& publication)
	{
		publications_.push_back(publication);
	}

  void MRMExperiment::setInstruments(const vector<CVTermList>& instruments)
  {
    instruments_ = instruments;
  }

  const vector<CVTermList>& MRMExperiment::getInstruments() const
  {
    return instruments_;
  }

  void MRMExperiment::addInstrument(const CVTermList& instrument)
  {
    instruments_.push_back(instrument);
  }

	void MRMExperiment::setSoftware(const vector<Software>& software)
  {
    software_ = software;
  }

  const vector<Software>& MRMExperiment::getSoftware() const
  {
    return software_;
  }

  void MRMExperiment::addSoftware(const Software& software)
  {
    software_.push_back(software);
  }

  void MRMExperiment::setProteins(const vector<Protein>& proteins)
  {
    proteins_ = proteins;
  }

  const vector<MRMExperiment::Protein>& MRMExperiment::getProteins() const
  {
    return proteins_;
  }

  void MRMExperiment::addProtein(const Protein& protein)
  {
    proteins_.push_back(protein);
  }

	void MRMExperiment::setCompounds(const vector<Compound>& compounds)
	{
		compounds_ = compounds;
	}

	const vector<MRMExperiment::Compound>& MRMExperiment::getCompounds() const
	{
		return compounds_;
	}

	void MRMExperiment::addCompound(const Compound& rhs)
	{
		compounds_.push_back(rhs);
	}

  void MRMExperiment::setPeptides(const vector<Peptide>& peptides)
  {
    peptides_ = peptides;
  }

  const vector<MRMExperiment::Peptide>& MRMExperiment::getPeptides() const
  {
    return peptides_;
  }

  void MRMExperiment::addPeptide(const Peptide& rhs)
  {
    peptides_.push_back(rhs);
  }

	void MRMExperiment::setTransitions(const vector<ReactionMonitoringTransition>& transitions)
	{
		transitions_ = transitions;
	}

	const vector<ReactionMonitoringTransition>& MRMExperiment::getTransitions() const
	{
		return transitions_;
	}

	void MRMExperiment::addTransition(const ReactionMonitoringTransition& transition)
	{
		transitions_.push_back(transition);
	}

	void MRMExperiment::setSourceFiles(const vector<SourceFile>& source_files)
	{
		source_files_ = source_files;
	}

	const vector<SourceFile>& MRMExperiment::getSourceFiles() const
	{
		return source_files_;
	}

	void MRMExperiment::addSourceFile(const SourceFile& source_file)
	{
		source_files_.push_back(source_file);
	}
}


