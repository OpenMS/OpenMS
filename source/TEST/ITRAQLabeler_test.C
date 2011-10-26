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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ITRAQLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ITRAQLabeler* ptr = 0;
ITRAQLabeler* null_ptr = 0;
START_SECTION(ITRAQLabeler())
{
	ptr = new ITRAQLabeler();
	TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->getParameters().getValue("iTRAQ"), "4plex")
}
END_SECTION

START_SECTION(virtual ~ITRAQLabeler())
{
	delete ptr;
}
END_SECTION

START_SECTION((void preCheck(Param &param) const ))
{
  ITRAQLabeler i;
  Param p;
  p.setValue("RawTandemSignal:status", "MS^E");
  TEST_EXCEPTION(Exception::InvalidParameter, i.preCheck(p));
  
  p.setValue("RawTandemSignal:status", "precursor");
  i.preCheck(p); // should work
}
END_SECTION

START_SECTION((void setUpHook(FeatureMapSimVector &)))
{
  ITRAQLabeler i;
  // check for correct number of channels
  FeatureMapSimVector f_maps;
  f_maps.push_back(FeatureMap<>());
  i.setUpHook(f_maps);

  // add another map
  Param p = i.getParameters();
  p.setValue("channel_active_4plex", StringList::create("114:myReference, 117:blabla"), "Four-plex only: Each channel that was used in the experiment and its description (114-117) in format <channel>:<name>, e.g. \"114:myref\",\"115:liver\"."); 
  i.setParameters(p);
  f_maps.push_back(FeatureMap<>());
  i.setUpHook(f_maps);

  // if no Exception until here, all is good

  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDigestHook(FeatureMapSimVector &)))
{
  // TODO
}
END_SECTION

START_SECTION((void postRTHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postDetectabilityHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postIonizationHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawMSHook(FeatureMapSimVector &)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  // TODO
}
END_SECTION

START_SECTION((static BaseLabeler* create()))
{
  BaseLabeler* labeler = ITRAQLabeler::create();
  TEST_NOT_EQUAL(labeler, 0)
  delete labeler;
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  ITRAQLabeler i;
  TEST_EQUAL(i.getProductName(), "itraq")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



