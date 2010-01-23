// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RawMSSignalSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RawMSSignalSimulation* ptr = 0;
START_SECTION((RawMSSignalSimulation(const gsl_rng *random_generator)))
{
	ptr = new RawMSSignalSimulation(NULL);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~RawMSSignalSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((RawMSSignalSimulation(const RawMSSignalSimulation &source)))
{
  RawMSSignalSimulation source(NULL);
  Param p = source.getParameters();
  p.setValue("peak_fwhm",0.3);
  source.setParameters(p);
  
  RawMSSignalSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((RawMSSignalSimulation& operator=(const RawMSSignalSimulation &source)))
{
  RawMSSignalSimulation source(NULL);
  RawMSSignalSimulation target(source);
  
  Param p = source.getParameters();
  p.setValue("peak_fwhm",0.3);
  source.setParameters(p);
  TEST_NOT_EQUAL(source.getParameters(), target.getParameters())
  
  target = source;

  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((void generateRawSignals(FeatureMapSim &, MSSimExperiment &)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



