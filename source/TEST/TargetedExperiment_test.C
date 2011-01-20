// -*- mode: C++; tab-width: 2; -*-
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TargetedExperiment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TargetedExperiment* ptr = 0;
START_SECTION(TargetedExperiment())
{
	ptr = new TargetedExperiment();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(virtual ~TargetedExperiment())
{
	delete ptr;
}
END_SECTION

START_SECTION((TargetedExperiment(const TargetedExperiment &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const TargetedExperiment &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCVs(const std::vector< CV > &cvs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CV>& getCVs() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addCV(const CV &cv)))
{
  // TODO
}
END_SECTION

START_SECTION((void setContacts(const std::vector< CVTermList > &contacts)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getContacts() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addContact(const CVTermList &contact)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPublications(const std::vector< CVTermList > &publications)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getPublications() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addPublication(const CVTermList &publication)))
{
  // TODO
}
END_SECTION

START_SECTION((void setInstruments(const std::vector< CVTermList > &instruments)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getInstruments() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInstrument(const CVTermList &instrument)))
{
  // TODO
}
END_SECTION

START_SECTION((void setSoftware(const std::vector< Software > &software)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Software>& getSoftware() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addSoftware(const Software &software)))
{
  // TODO
}
END_SECTION

START_SECTION((void setProteins(const std::vector< Protein > &proteins)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Protein>& getProteins() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addProtein(const Protein &protein)))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompounds(const std::vector< Compound > &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Compound>& getCompounds() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addCompound(const Compound &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptides(const std::vector< Peptide > &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Peptide>& getPeptides() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addPeptide(const Peptide &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setTransitions(const std::vector< ReactionMonitoringTransition > &transitions)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<ReactionMonitoringTransition>& getTransitions() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addTransition(const ReactionMonitoringTransition &transition)))
{
  // TODO
}
END_SECTION

START_SECTION((TargetedExperiment& operator=(const TargetedExperiment &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



