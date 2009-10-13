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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(EnzymaticDigestion, "$Id$")

/////////////////////////////////////////////////////////////

EnzymaticDigestion* e_ptr = 0;
START_SECTION((EnzymaticDigestion()))
	e_ptr = new EnzymaticDigestion;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION([EXTRA] ~EnzymaticDigestion())
	delete e_ptr;
END_SECTION

START_SECTION((Size getMissedCleavages() const ))
	TEST_EQUAL(EnzymaticDigestion().getMissedCleavages(),0)
END_SECTION

START_SECTION((Enzyme getEnzyme() const))
	TEST_EQUAL(EnzymaticDigestion().getEnzyme(),EnzymaticDigestion::TRYPSIN)
END_SECTION

START_SECTION((void setMissedCleavages(Size missed_cleavages)))
	EnzymaticDigestion ed;
	ed.setMissedCleavages(5);
	TEST_EQUAL(ed.getMissedCleavages(),5)
END_SECTION

START_SECTION((void setEnzyme(Enzyme enzyme)))
	//can be tested as soon as there is a second enzyme
	NOT_TESTABLE
END_SECTION

START_SECTION((Enzyme getEnzymeByName(const String& name)))
	EnzymaticDigestion ed;
	TEST_EQUAL(ed.getEnzymeByName("Trypsin"), EnzymaticDigestion::TRYPSIN);
	TEST_EQUAL(ed.getEnzymeByName("DoesNotExist"), EnzymaticDigestion::SIZE_OF_ENZYMES);
END_SECTION

START_SECTION((Size peptideCount(const AASequence &protein)))
	EnzymaticDigestion ed;
	Size tmp = ed.peptideCount(String("ACDE"));
	TEST_EQUAL(tmp,1)
	tmp = ed.peptideCount(String("ACKDE"));
	TEST_EQUAL(tmp,2)
	tmp = ed.peptideCount(String("ACRDE"));
	TEST_EQUAL(tmp,2)
	tmp = ed.peptideCount(String("ACKPDE"));
	TEST_EQUAL(tmp,1)
	tmp = ed.peptideCount(String("ACRPDE"));
	TEST_EQUAL(tmp,1)
	tmp = ed.peptideCount(String("ARCRDRE"));
	TEST_EQUAL(tmp,4)
	tmp = ed.peptideCount(String("RKR"));
	TEST_EQUAL(tmp,3)
	ed.setMissedCleavages(1);
	TEST_EQUAL(ed.peptideCount(String("ACDE")),1)
	TEST_EQUAL(ed.peptideCount(String("ACRDE")),3)
	TEST_EQUAL(ed.peptideCount(String("ARCDRE")),5)
	TEST_EQUAL(ed.peptideCount(String("RKR")),5)
	ed.setMissedCleavages(3);
	TEST_EQUAL(ed.peptideCount(String("ACDE")),1)
	TEST_EQUAL(ed.peptideCount(String("ACRDE")),3)
	TEST_EQUAL(ed.peptideCount(String("ARCDRE")),6)
	TEST_EQUAL(ed.peptideCount(String("RKR")),6)
END_SECTION

START_SECTION((void digest(const AASequence& protein, std::vector<AASequence>& output)))
	EnzymaticDigestion ed;
	vector<AASequence> out;
	
	ed.digest(String("ACDE"),out);
	TEST_EQUAL(out.size(),1)
	TEST_EQUAL(out[0].toString(),"ACDE")

	ed.digest(String("ACKDE"),out);
	TEST_EQUAL(out.size(),2)
	TEST_EQUAL(out[0].toString(),"ACK")
	TEST_EQUAL(out[1].toString(),"DE")
	
	ed.digest(String("ACRDE"),out);
	TEST_EQUAL(out.size(),2)
	TEST_EQUAL(out[0].toString(),"ACR")
	TEST_EQUAL(out[1].toString(),"DE")

	ed.digest(String("ACKPDE"),out);
	TEST_EQUAL(out.size(),1)
	TEST_EQUAL(out[0].toString(),"ACKPDE")
	
	ed.digest(String("ACRPDE"),out);
	TEST_EQUAL(out.size(),1)
	TEST_EQUAL(out[0].toString(),"ACRPDE")
	
	ed.digest(String("ARCRDRE"),out);
	TEST_EQUAL(out.size(),4)
	TEST_EQUAL(out[0].toString(),"AR")
	TEST_EQUAL(out[1].toString(),"CR")
	TEST_EQUAL(out[2].toString(),"DR")
	TEST_EQUAL(out[3].toString(),"E")

	ed.digest(String("RKR"),out);
	TEST_EQUAL(out.size(),3)
	TEST_EQUAL(out[0].toString(),"R")
	TEST_EQUAL(out[1].toString(),"K")
	TEST_EQUAL(out[2].toString(),"R")

	ed.setMissedCleavages(1);
	
	ed.digest(String("ACDE"),out);
	TEST_EQUAL(out.size(),1)
	TEST_EQUAL(out[0].toString(),"ACDE")
	
	ed.digest(String("ACRDE"),out);
	TEST_EQUAL(out.size(),3)
	TEST_EQUAL(out[0].toString(),"ACR")
	TEST_EQUAL(out[1].toString(),"DE")
	TEST_EQUAL(out[2].toString(),"ACRDE")
	
	ed.digest(String("ARCDRE"),out);
	TEST_EQUAL(out.size(),5)
	TEST_EQUAL(out[0].toString(),"AR")
	TEST_EQUAL(out[1].toString(),"CDR")
	TEST_EQUAL(out[2].toString(),"E")
	TEST_EQUAL(out[3].toString(),"ARCDR")
	TEST_EQUAL(out[4].toString(),"CDRE")
	
	ed.digest(String("RKR"),out);
	TEST_EQUAL(out.size(),5)
	TEST_EQUAL(out[0].toString(),"R")
	TEST_EQUAL(out[1].toString(),"K")
	TEST_EQUAL(out[2].toString(),"R")
	TEST_EQUAL(out[3].toString(),"RK")
	TEST_EQUAL(out[4].toString(),"KR")
	
//	ed.setMissedCleavages(3);
//	
//	ed.digest(String("ACDE"),out);
//	TEST_EQUAL(out.size(),1)
//	TEST_EQUAL(out[0],"ACDE")
//	
//	ed.digest(String("ACRDE"),out);
//	TEST_EQUAL(out.size(),3)
//	TEST_EQUAL(out[0],"ACR")
//	TEST_EQUAL(out[1],"DE")
//	TEST_EQUAL(out[2],"ACRDE")
//	
//	ed.digest(String("ARCDRE"),out);
//	TEST_EQUAL(out.size(),6)
//	TEST_EQUAL(out[0],"AR")
//	TEST_EQUAL(out[1],"CDR")
//	TEST_EQUAL(out[2],"E")
//	TEST_EQUAL(out[3],"ARCDR")
//	TEST_EQUAL(out[4],"CDRE")
//	TEST_EQUAL(out[5],"ARCDRE")
//	
//	ed.digest(String("RKR"),out);
//	TEST_EQUAL(out.size(),6)
//	TEST_EQUAL(out[0],"R")
//	TEST_EQUAL(out[1],"K")
//	TEST_EQUAL(out[2],"R")
//	TEST_EQUAL(out[3],"RK")
//	TEST_EQUAL(out[4],"KR")
//	TEST_EQUAL(out[5],"RKR")


END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
