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
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ReactionMonitoringTransition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ReactionMonitoringTransition* ptr = 0;
ReactionMonitoringTransition* nullPointer = 0;
START_SECTION(ReactionMonitoringTransition())
{
	ptr = new ReactionMonitoringTransition();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~ReactionMonitoringTransition())
{
	delete ptr;
}
END_SECTION

START_SECTION((ReactionMonitoringTransition(const ReactionMonitoringTransition &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorMZ(DoubleReal mz)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getPrecursorMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductMZ(DoubleReal mz)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getProductMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getProductCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setInterpretations(const std::vector< TransitionInterpretation > &interpretations)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<TransitionInterpretation>& getInterpretations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInterpretation(const TransitionInterpretation &interpretation)))
{
  // TODO
}
END_SECTION

START_SECTION((void setConfigurations(const std::vector< Configuration > &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Configuration>& getConfigurations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addConfiguration(const Configuration &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const ReactionMonitoringTransition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const ReactionMonitoringTransition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ReactionMonitoringTransition& operator=(const ReactionMonitoringTransition &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



