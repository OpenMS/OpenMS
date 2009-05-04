// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/MSSim.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MSSim, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSim* ptr = 0;
START_SECTION(MSSim())
{
	ptr = new MSSim();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~MSSim())
{
	delete ptr;
}
END_SECTION

START_SECTION((MSSim(const MSSim &source)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~MSSim()))
{
  // TODO
}
END_SECTION

START_SECTION((MSSim& operator=(const MSSim &source)))
{
  // TODO
}
END_SECTION

START_SECTION((void simulate(const gsl_rng *rnd_gen, const SamplePeptides &peptides)))
{
  // TODO
}
END_SECTION

START_SECTION((MSSimExperiment const& getExperiment() const ))
{
  // TODO
}
END_SECTION

START_SECTION((FeatureMapSim const& getSimulatedFeatures() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



