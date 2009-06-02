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
		:	transitions_(rhs.transitions_)
	{
	}

	MRMExperiment::~MRMExperiment()
	{
	}

	MRMExperiment& MRMExperiment::operator = (const MRMExperiment& rhs)
	{
		if (&rhs != this)
		{
			transitions_ = rhs.transitions_;
		}
		return *this;
	}

	bool MRMExperiment::operator == (const MRMExperiment& rhs) const
	{
		return transitions_ == rhs.transitions_;
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

	void MRMExperiment::setPublications(const vector<MetaInfoInterface>& publications)
	{
		publications_ = publications;
	}

	const vector<MetaInfoInterface>& MRMExperiment::getPublications() const
	{
		return publications_;
	}

	void MRMExperiment::addPublication(const MetaInfoInterface& publication)
	{
		publications_.push_back(publication);
	}

  void MRMExperiment::setInstruments(const vector<MetaInfoInterface>& instruments)
  {
    instruments_ = instruments;
  }

  const vector<MetaInfoInterface>& MRMExperiment::getInstruments() const
  {
    return instruments_;
  }

  void MRMExperiment::addInstrument(const MetaInfoInterface& instrument)
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
}


