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
#include <OpenMS/SIMULATION/PTMSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PTMSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PTMSimulation* ptr = 0;
START_SECTION(PTMSimulation(const gsl_rng *rnd_gen))
{
	ptr = new PTMSimulation(NULL);
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~PTMSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((PTMSimulation(const PTMSimulation &source)))
{
  PTMSimulation a(NULL);

	Param p = a.getParameters();
	p.setValue("potential_modifications", StringList::create("MOD:00071|0.100003"));
	a.setParameters(p);
	PTMSimulation b(a);
	TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION

START_SECTION((virtual ~PTMSimulation()))
{
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((PTMSimulation& operator=(const PTMSimulation &source)))
{
  PTMSimulation a(NULL);
	PTMSimulation b(a);

	Param p = a.getParameters();
	p.setValue("potential_modifications", StringList::create("MOD:00071|0.13"));
	a.setParameters(p);
	TEST_NOT_EQUAL(b.getParameters(),a.getParameters());
	b = a;
	TEST_EQUAL(b.getParameters(),a.getParameters());
}
END_SECTION

START_SECTION((void predictPTMs(FeatureMapSim &map)))
{
	gsl_rng* rnd_gen = gsl_rng_alloc (gsl_rng_taus);

  PTMSimulation a(rnd_gen);
	
	FeatureMapSim map;
	map.reserve(3);
	StringList peps = StringList::create("ACHKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKHHACAC,AAAAAAAAHTKLRTTIPPEFG,RRRRRRRRRYCNHKTUIKL");
	for (StringList::const_iterator it=peps.begin(); it!=peps.end(); ++it)
	{
		Feature f;
		PeptideIdentification pep_id;
		pep_id.insertHit(PeptideHit(1.0, 1, 1, *it));
		f.getPeptideIdentifications().push_back(pep_id);
		f.setIntensity(1000);
		map.push_back(f);
	}

	a.predictPTMs(map);

	TEST_EQUAL(map[0].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACHKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKHHACAC")
	TEST_EQUAL(map[1].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "AAAAAAAAHTKLRTTIPPEFG")
	TEST_EQUAL(map[2].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "RRRRRRRRRYCNHKTUIKL")
	TEST_EQUAL(map[3].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "ACHK(Lys->Allysine (K))K(Lys->Allysine (K))K(Lys->Allysine (K))KKKKKKKKKKKKKKKKKKKKKKKKKKKHHACAC")
	TEST_EQUAL(map[4].getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "R(Methyl (R))RRRRRRRRYCNHK(Lys->Allysine (K))TUIKL")
	TEST_EQUAL(map.size(),5);
	//for (FeatureMapSim::const_iterator it = map.begin(); it!=map.end(); ++it)
	//{
	//	std::cout << it->getPeptideIdentifications()[0].getHits()[0].getSequence().toString() << "\n";
	//}


	gsl_rng_free (rnd_gen);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



