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
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BaseLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BaseLabeler* ptr = 0;
BaseLabeler* nullPointer = 0;
START_SECTION(BaseLabeler())
{
	ptr = new BaseLabeler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~BaseLabeler())
{
	delete ptr;
}
END_SECTION

BaseLabeler labeler;
FeatureMapSimVector empty_fmsv;
MSSimExperiment empty_experiment;

START_SECTION((virtual void setUpHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.setUpHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postDigestHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postDigestHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRTHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRTHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postDetectabilityHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postDetectabilityHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postIonizationHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postIonizationHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRawMSHook(FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRawMSHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRawTandemMSHook(empty_fmsv, empty_experiment))
}
END_SECTION

START_SECTION((virtual Param getDefaultParameters() const ))
{
  Param p; // empty parameters
  TEST_EQUAL(labeler.getDefaultParameters(), p) // BaseLabeler should not have any parameters
}
END_SECTION

START_SECTION((virtual void setRnd(const SimRandomNumberGenerator &rng)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void preCheck(Param &) const ))
{
  Param p;
  TEST_EXCEPTION(Exception::NotImplemented, labeler.preCheck(p))
}
END_SECTION

START_SECTION((const ConsensusMap& getConsensus() const ))
{
  ConsensusMap cm;
  TEST_EQUAL(labeler.getConsensus(), cm) // Consensus should be empty
}
END_SECTION

START_SECTION((String getChannelIntensityName(const Size channel_index) const ))
{
  TEST_STRING_EQUAL(labeler.getChannelIntensityName(1), "channel_1_intensity")
  TEST_STRING_EQUAL(labeler.getChannelIntensityName(100), "channel_100_intensity")
}
END_SECTION

START_SECTION((void registerChildren()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



