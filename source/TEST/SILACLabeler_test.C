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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/SILACLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SILACLabeler* ptr = 0;
START_SECTION(SILACLabeler())
{
	ptr = new SILACLabeler();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~SILACLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
  NOT_TESTABLE // empty hook
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
  FeatureMapSimVector illegal;
  illegal.push_back(FeatureMapSim());

  SILACLabeler exception_labeler;
  Param empty_param;
  exception_labeler.preCheck(empty_param);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, exception_labeler.setUpHook(illegal), "We currently support only 2-channel SILAC")

  illegal.push_back(FeatureMapSim());
  illegal.push_back(FeatureMapSim());

  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, exception_labeler.setUpHook(illegal), "We currently support only 2-channel SILAC")

  // TODO .. test the real functionality
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE // empty hook
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE // empty hook
}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  NOT_TESTABLE // empty hook
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = SILACLabeler::create();
  TEST_NOT_EQUAL(labeler, 0)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(SILACLabeler::getProductName(), "SILAC")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



