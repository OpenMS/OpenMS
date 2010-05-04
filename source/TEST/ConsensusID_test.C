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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm, Andreas Bertsch, Sven Nahnsen $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/ConsensusID.h>
#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ResidueDB, "$Id$")

/////////////////////////////////////////////////////////////

ConsensusID* ptr = 0;
START_SECTION(ConsensusID())
	ptr = new ConsensusID();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~ConsensusID())
	delete(ptr);
END_SECTION

// PeptideIdentification with 3 id runs is created
vector<PeptideIdentification> ids(3);
vector<PeptideHit> hits;
cout<<"HELLO"<<endl;
// the first ID has 5 hits
hits.resize(5);
hits[0].setRank(1);
hits[0].setSequence("A");
hits[0].setScore(31);
hits[1].setRank(2);
hits[1].setSequence("B");	
hits[1].setScore(28);
hits[2].setRank(3);
hits[2].setSequence("C");
hits[2].setScore(17);
hits[3].setRank(4);
hits[3].setSequence("D");
hits[3].setScore(7);
hits[4].setRank(5);
hits[4].setSequence("E");
hits[4].setScore(3);
ids[0].setHits(hits);
// the second ID has 3 hits
hits.resize(3);
hits[0].setRank(1);
hits[0].setSequence("C");
hits[0].setScore(32);
hits[1].setRank(2);
hits[1].setSequence("A");
hits[1].setScore(30);
hits[2].setRank(3);
hits[2].setSequence("B");
hits[2].setScore(29);
ids[1].setHits(hits);
// the third ID has 10 hits
hits.resize(10);
hits[0].setRank(1);
hits[0].setSequence("F");
hits[0].setScore(81);
hits[1].setRank(2);
hits[1].setSequence("C");	
hits[1].setScore(60);
hits[2].setRank(3);
hits[2].setSequence("G");
hits[2].setScore(50);
hits[3].setRank(4);
hits[3].setSequence("D");
hits[3].setScore(40);
hits[4].setRank(5);
hits[4].setSequence("B");
hits[4].setScore(25);
hits[5].setRank(6);
hits[5].setSequence("E");
hits[5].setScore(5);
hits[6].setRank(7);
hits[6].setSequence("H");	
hits[6].setScore(4);
hits[7].setRank(8);
hits[7].setSequence("I");
hits[7].setScore(3);
hits[8].setRank(9);
hits[8].setSequence("J");
hits[8].setScore(2);
hits[9].setRank(10);
hits[9].setSequence("K");
hits[9].setScore(1);
ids[2].setHits(hits);

START_SECTION(void apply(std::vector<PeptideIdentification>& ids))
	TOLERANCE_ABSOLUTE(0.01)
	
	// ***** Ranked ********
	
	ConsensusID consensus;
	//define parameters
	Param param;
	param.setValue("algorithm","ranked");
	param.setValue("considered_hits",5);
	consensus.setParameters(param);
	//apply
	vector<PeptideIdentification> f = ids;
	consensus.apply(f);
	
	TEST_EQUAL(f.size(),1);
	hits = f[0].getHits();
	TEST_EQUAL(hits.size(),7);
	
	TEST_EQUAL(hits[0].getRank(),1);
	TEST_EQUAL(hits[0].getSequence(),"C");
	TEST_REAL_SIMILAR(hits[0].getScore(),80.0f);
	
	TEST_EQUAL(hits[1].getRank(),2);
	TEST_EQUAL(hits[1].getSequence(),"A");
	TEST_REAL_SIMILAR(hits[1].getScore(),60.0f);
	
	TEST_EQUAL(hits[2].getRank(),3);
	TEST_EQUAL(hits[2].getSequence(),"B");
	TEST_REAL_SIMILAR(hits[2].getScore(),53.33f);
	
	TEST_EQUAL(hits[3].getRank(),4);
	TEST_EQUAL(hits[3].getSequence(),"F");
	TEST_REAL_SIMILAR(hits[3].getScore(),33.333f);
	
	TEST_EQUAL(hits[4].getRank(),5);
	TEST_EQUAL(hits[4].getSequence(),"D");
	TEST_REAL_SIMILAR(hits[4].getScore(),26.666f);
	
	TEST_EQUAL(hits[5].getRank(),6);
	TEST_EQUAL(hits[5].getSequence(),"G");
	TEST_REAL_SIMILAR(hits[5].getScore(),20.0f);
	
	TEST_EQUAL(hits[6].getRank(),7);
	TEST_EQUAL(hits[6].getSequence(),"E");
	TEST_REAL_SIMILAR(hits[6].getScore(),6.666f);

	// ***** Merge ********

	//define parameters
	param.clear();
	param.setValue("algorithm","merge");
	param.setValue("considered_hits",6);
	consensus.setParameters(param);
	//apply
	f = ids;
	consensus.apply(f);
	
	TEST_EQUAL(f.size(),1);
	hits = f[0].getHits();
	TEST_EQUAL(hits.size(),7);
	
	TEST_EQUAL(hits[0].getRank(),1);
	TEST_EQUAL(hits[0].getSequence(),"F");
	TEST_REAL_SIMILAR(hits[0].getScore(),81.0f);
	
	TEST_EQUAL(hits[1].getRank(),2);
	TEST_EQUAL(hits[1].getSequence(),"C");
	TEST_REAL_SIMILAR(hits[1].getScore(),60.0f);
	
	TEST_EQUAL(hits[2].getRank(),3);
	TEST_EQUAL(hits[2].getSequence(),"G");
	TEST_REAL_SIMILAR(hits[2].getScore(),50.0f);
	
	TEST_EQUAL(hits[3].getRank(),4);
	TEST_EQUAL(hits[3].getSequence(),"D");
	TEST_REAL_SIMILAR(hits[3].getScore(),40.0f);
	
	TEST_EQUAL(hits[4].getRank(),5);
	TEST_EQUAL(hits[4].getSequence(),"A");
	TEST_REAL_SIMILAR(hits[4].getScore(),31.0f);
	
	TEST_EQUAL(hits[5].getRank(),6);
	TEST_EQUAL(hits[5].getSequence(),"B");
	TEST_REAL_SIMILAR(hits[5].getScore(),29.0f);
	
	TEST_EQUAL(hits[6].getRank(),7);
	TEST_EQUAL(hits[6].getSequence(),"E");
	TEST_REAL_SIMILAR(hits[6].getScore(),5.0f);

	// ***** Average ********

	//define parameters
	param.clear();
	param.setValue("algorithm","average");
	param.setValue("considered_hits",4);
	consensus.setParameters(param);
	//apply
	f = ids;
	consensus.apply(f);
	
	TEST_EQUAL(f.size(),1);
	hits = f[0].getHits();
	TEST_EQUAL(hits.size(),6);
	
	TEST_EQUAL(hits[0].getRank(),1);
	TEST_EQUAL(hits[0].getSequence(),"C");
	TEST_REAL_SIMILAR(hits[0].getScore(),36.333f);
	
	TEST_EQUAL(hits[1].getRank(),2);
	TEST_EQUAL(hits[1].getSequence(),"F");
	TEST_REAL_SIMILAR(hits[1].getScore(),27.0f);
	
	TEST_EQUAL(hits[2].getRank(),3);
	TEST_EQUAL(hits[2].getSequence(),"A");
	TEST_REAL_SIMILAR(hits[2].getScore(),20.333f);
	
	TEST_EQUAL(hits[3].getRank(),4);
	TEST_EQUAL(hits[3].getSequence(),"B");
	TEST_REAL_SIMILAR(hits[3].getScore(),19.0f);
	
	TEST_EQUAL(hits[4].getRank(),5);
	TEST_EQUAL(hits[4].getSequence(),"G");
	TEST_REAL_SIMILAR(hits[4].getScore(),16.666f);
	
	TEST_EQUAL(hits[5].getRank(),6);
	TEST_EQUAL(hits[5].getSequence(),"D");
	TEST_REAL_SIMILAR(hits[5].getScore(),15.666f);

	// ***** Average, Inverse Order ********

	//define parameters
	param.clear();
	param.setValue("algorithm","average");
	param.setValue("considered_hits",1);
	consensus.setParameters(param);
	f = ids;
	for (Size i = 0; i < f.size(); ++i)
	{
		f[i].setHigherScoreBetter(false);
	}
	//apply
	consensus.apply(f);
	
	TEST_EQUAL(f.size(),1);
	hits = f[0].getHits();
	TEST_EQUAL(hits.size(),3);
	
	TEST_EQUAL(hits[0].getRank(),1);
	TEST_EQUAL(hits[0].getSequence(),"K");
	TEST_REAL_SIMILAR(hits[0].getScore(),0.333f);
	
	TEST_EQUAL(hits[1].getRank(),2);
	TEST_EQUAL(hits[1].getSequence(),"E");
	TEST_REAL_SIMILAR(hits[1].getScore(),1.0f);
	
	TEST_EQUAL(hits[2].getRank(),3);
	TEST_EQUAL(hits[2].getSequence(),"B");
	TEST_REAL_SIMILAR(hits[2].getScore(),9.666f);
	
	// ***** probability ********

/*	//define parameters
	param.clear();
	param.setValue("algorithm","probability");
	param.setValue("considered_hits",5);
	consensus.setParameters(param);
	//apply
	f = ids;
	consensus.apply(f);
	
	TEST_EQUAL(f.size(),1);
	hits = f[0].getHits();
	TEST_EQUAL(hits.size(),7);
	
	TEST_EQUAL(hits[0].getRank(),1);
	TEST_EQUAL(hits[0].getSequence(),"C");
	TEST_REAL_SIMILAR(hits[0].getScore(),32640.0f);
	
	TEST_EQUAL(hits[1].getRank(),2);
	TEST_EQUAL(hits[1].getSequence(),"B");
	TEST_REAL_SIMILAR(hits[1].getScore(),20300.0f);
	
	TEST_EQUAL(hits[2].getRank(),3);
	TEST_EQUAL(hits[2].getSequence(),"A");
	TEST_REAL_SIMILAR(hits[2].getScore(),930.0f);
	
	TEST_EQUAL(hits[3].getRank(),4);
	TEST_EQUAL(hits[3].getSequence(),"D");
	TEST_REAL_SIMILAR(hits[3].getScore(),280.0f);
	
	TEST_EQUAL(hits[4].getRank(),5);
	TEST_EQUAL(hits[4].getSequence(),"F");
	TEST_REAL_SIMILAR(hits[4].getScore(),81.0f);
	
	TEST_EQUAL(hits[5].getRank(),6);
	TEST_EQUAL(hits[5].getSequence(),"G");
	TEST_REAL_SIMILAR(hits[5].getScore(),50.0f);
	
*/
	// ***** Exception ********
	param.setValue("algorithm","Bla4711");
	TEST_EXCEPTION(Exception::InvalidParameter,consensus.setParameters(param));
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
