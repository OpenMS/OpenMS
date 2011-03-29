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
SILACLabeler* nullPointer = 0;
BaseLabeler*      base_nullPointer = 0;

START_SECTION(SILACLabeler())
{
	ptr = new SILACLabeler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SILACLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

// just to call the methods once
SILACLabeler dummyLabeler;
FeatureMapSimVector empty;

START_SECTION((void preCheck(Param &param) const ))
{
  Param p;
  dummyLabeler.preCheck(p);

  // preCheck has no content
  NOT_TESTABLE
}
END_SECTION


START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postRTHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postDetectabilityHook(empty);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
  // we do not modify the map in this step
  dummyLabeler.postIonizationHook(empty);
  NOT_TESTABLE
}
END_SECTION

MSSimExperiment exp;
START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  // we do not modify the map in this step
  dummyLabeler.postRawTandemMSHook(empty,exp);
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = SILACLabeler::create();
  TEST_NOT_EQUAL(labeler, base_nullPointer)
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



