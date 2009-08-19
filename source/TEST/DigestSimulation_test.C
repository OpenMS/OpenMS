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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/DigestSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DigestSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DigestSimulation* ptr = 0;
START_SECTION(DigestSimulation())
{
	ptr = new DigestSimulation();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~DigestSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((DigestSimulation(const DigestSimulation &source)))
{
  DigestSimulation a;
	Param p = a.getParameters();
	p.setValue("enzyme","none");
	a.setParameters(p);
	DigestSimulation b(a);

	TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION

START_SECTION((DigestSimulation& operator=(const DigestSimulation &source)))
{
  DigestSimulation a,b;
	Param p = a.getParameters();
	p.setValue("enzyme","none");
	a.setParameters(p);

	TEST_NOT_EQUAL(b.getParameters(),a.getParameters());
	b = a;
	TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION


START_SECTION((void digest(FeatureMapSim & feature_map)))
{
  SampleProteins in,out;
  // TODO: re-enable test code  
/* 	in["ACDKDDLDDFRLNN"] = 100;
 * 	in["ACDKDDLASSRL"] = 50;
 * 
 * 	// test if "out" is deleted
 * 	out["AALA"] = 10000;
 * 
 * 	DigestSimulation a;
 * 	a.digest(in,out);
 * 
 * 	TEST_REAL_SIMILAR(out["ACDK"], 107.143);
 * 	TEST_REAL_SIMILAR(out["ACDKDDLASSR"], 35.7143);
 * 	TEST_REAL_SIMILAR(out["ACDKDDLDDFR"], 71.4286);
 * 	TEST_REAL_SIMILAR(out["DDLASSR"], 35.7143);
 * 	TEST_REAL_SIMILAR(out["DDLASSRL"], 35.7143);
 * 	TEST_REAL_SIMILAR(out["DDLDDFR"], 71.4286);
 * 	TEST_REAL_SIMILAR(out["DDLDDFRLNN"], 71.4286);
 * 	TEST_REAL_SIMILAR(out["LNN"], 71.4286);
 */

 //for (SampleProteins::const_iterator protein = out.begin();
 //        protein != out.end();
 //        ++protein)
 //{
 //  std::cout << protein->first.toString() << " " << protein->second << "\n";
 //}

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



