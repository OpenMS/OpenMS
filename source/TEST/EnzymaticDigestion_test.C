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
// $Authors: Marc Sturm, Chris Bielow $
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

START_SECTION((SignedSize getMissedCleavages() const ))
	TEST_EQUAL(EnzymaticDigestion().getMissedCleavages(),0)
END_SECTION

START_SECTION((Enzyme getEnzyme() const))
	TEST_EQUAL(EnzymaticDigestion().getEnzyme(),EnzymaticDigestion::TRYPSIN)
END_SECTION

START_SECTION((void setMissedCleavages(SignedSize missed_cleavages)))
	EnzymaticDigestion ed;
	ed.setMissedCleavages(5);
	TEST_EQUAL(ed.getMissedCleavages(),5)
END_SECTION

START_SECTION((void setEnzyme(Enzyme enzyme)))
	EnzymaticDigestion ed;
  ed.setEnzyme(EnzymaticDigestion::TRYPSIN);
  TEST_EQUAL(ed.getEnzyme(), EnzymaticDigestion::TRYPSIN);
  ed.setEnzyme(EnzymaticDigestion::SIZE_OF_ENZYMES);
  TEST_EQUAL(ed.getEnzyme(), EnzymaticDigestion::SIZE_OF_ENZYMES);
END_SECTION

START_SECTION((Enzyme getEnzymeByName(const String& name)))
	EnzymaticDigestion ed;
	TEST_EQUAL(ed.getEnzymeByName("Trypsin"), EnzymaticDigestion::TRYPSIN);
	TEST_EQUAL(ed.getEnzymeByName("DoesNotExist"), EnzymaticDigestion::SIZE_OF_ENZYMES);
END_SECTION

START_SECTION((bool isLogModelEnabled() const))
	EnzymaticDigestion ed;
  TEST_EQUAL(ed.isLogModelEnabled(), false);
END_SECTION

START_SECTION((void setLogModelEnabled(bool enabled)))
	EnzymaticDigestion ed;
  ed.setLogModelEnabled(true);
  TEST_EQUAL(ed.isLogModelEnabled(), true);
  ed.setLogModelEnabled(false);
  TEST_EQUAL(ed.isLogModelEnabled(), false);
END_SECTION

START_SECTION((DoubleReal getLogThreshold() const))
	EnzymaticDigestion ed;
  ed.setLogThreshold(1.234);
  TEST_EQUAL(ed.getLogThreshold(), 1.234);
END_SECTION

START_SECTION((void setLogThreshold(DoubleReal threshold)))
  // TESTED ABOVE
  NOT_TESTABLE
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

  // with log L model:
  ed.setLogModelEnabled(true);
  TEST_EQUAL(ed.peptideCount(String("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA")),9+1+1) // K R + 1

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
	

  // with log L model:
  ed.setLogModelEnabled(true);
  ed.digest(String("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA"), out);
	TEST_EQUAL(out.size(), 11)
	TEST_EQUAL(out[0].toString(),"MKWVTFISLLLLFSSAYSRGVFRRDTHK")
	TEST_EQUAL(out[1].toString(),"SEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLK")
	TEST_EQUAL(out[2].toString(),"PDPNTLCDEFKADEKK")
	TEST_EQUAL(out[3].toString(),"FWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQK")
	TEST_EQUAL(out[4].toString(),"FPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDK")
	TEST_EQUAL(out[5].toString(),"PLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRK")
	TEST_EQUAL(out[6].toString(),"VPQVSTPTLVEVSRSLGK")
	TEST_EQUAL(out[7].toString(),"VGTRCCTK")
	TEST_EQUAL(out[8].toString(),"PESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRR")
	TEST_EQUAL(out[9].toString(),"PCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHK")
	TEST_EQUAL(out[10].toString(),"PKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA")

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
