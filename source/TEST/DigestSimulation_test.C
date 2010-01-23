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
	FeatureMapSim fm;
  ProteinIdentification protIdent;
	// add new ProteinHit to ProteinIdentification
	{
  ProteinHit protHit(0.0, 1, "Hit1", "ACDKDDLDDFRLNN");
  protHit.setMetaValue("description", "desc 1");
	protHit.setMetaValue("intensity", 100);
  protIdent.insertHit(protHit);
	}
	{
  ProteinHit protHit(0.0, 1, "Hit2", "ACDKDDLASSRL");
  protHit.setMetaValue("description", "desc 1");
	protHit.setMetaValue("intensity", 50);
  protIdent.insertHit(protHit);
	}
  std::vector<ProteinIdentification> vec_protIdent;
  vec_protIdent.push_back(protIdent);
  fm.setProteinIdentifications(vec_protIdent);
     
	DigestSimulation a;
	a.digest(fm);
	  
	TEST_EQUAL(fm.size(), 8) 
	
	TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), "ACDK")
	TEST_EQUAL(fm[0].getIntensity(), 108)

	TEST_EQUAL(fm[1].getPeptideIdentifications()[0].getHits()[0].getSequence(), "ACDKDDLASSR")
	TEST_EQUAL(fm[1].getIntensity(), 36)

	TEST_EQUAL(fm[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), "ACDKDDLDDFR")
	TEST_EQUAL(fm[2].getIntensity(), 72)

	TEST_EQUAL(fm[3].getPeptideIdentifications()[0].getHits()[0].getSequence(), "DDLASSR")
	TEST_EQUAL(fm[3].getIntensity(), 36)

	TEST_EQUAL(fm[4].getPeptideIdentifications()[0].getHits()[0].getSequence(), "DDLASSRL")
	TEST_EQUAL(fm[4].getIntensity(), 36)

	TEST_EQUAL(fm[5].getPeptideIdentifications()[0].getHits()[0].getSequence(), "DDLDDFR")
	TEST_EQUAL(fm[5].getIntensity(), 72)

	TEST_EQUAL(fm[6].getPeptideIdentifications()[0].getHits()[0].getSequence(), "DDLDDFRLNN")
	TEST_EQUAL(fm[6].getIntensity(), 72)

	TEST_EQUAL(fm[7].getPeptideIdentifications()[0].getHits()[0].getSequence(), "LNN")
	TEST_EQUAL(fm[7].getIntensity(), 72)


 //for (FeatureMapSim::const_iterator f = fm.begin();
 //        f != fm.end();
 //        ++f)
 //{
 //  std::cout << f->getPeptideIdentifications()[0].getHits()[0].getSequence() << " " << f->getIntensity() << "\n";
 //}

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



