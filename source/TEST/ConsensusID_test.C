// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Sven Nahnsen $
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
ConsensusID* nullPointer = 0;
START_SECTION(ConsensusID())
	ptr = new ConsensusID();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
hits[0].setSequence(AASequence("A"));
hits[0].setScore(31);
hits[1].setRank(2);
hits[1].setSequence(AASequence("B"));
hits[1].setScore(28);
hits[2].setRank(3);
hits[2].setSequence(AASequence("C"));
hits[2].setScore(17);
hits[3].setRank(4);
hits[3].setSequence(AASequence("D"));
hits[3].setScore(7);
hits[4].setRank(5);
hits[4].setSequence(AASequence("E"));
hits[4].setScore(3);
ids[0].setHits(hits);
// the second ID has 3 hits
hits.resize(3);
hits[0].setRank(1);
hits[0].setSequence(AASequence("C"));
hits[0].setScore(32);
hits[1].setRank(2);
hits[1].setSequence(AASequence("A"));
hits[1].setScore(30);
hits[2].setRank(3);
hits[2].setSequence(AASequence("B"));
hits[2].setScore(29);
ids[1].setHits(hits);
// the third ID has 10 hits
hits.resize(10);
hits[0].setRank(1);
hits[0].setSequence(AASequence("F"));
hits[0].setScore(81);
hits[1].setRank(2);
hits[1].setSequence(AASequence("C"));
hits[1].setScore(60);
hits[2].setRank(3);
hits[2].setSequence(AASequence("G"));
hits[2].setScore(50);
hits[3].setRank(4);
hits[3].setSequence(AASequence("D"));
hits[3].setScore(40);
hits[4].setRank(5);
hits[4].setSequence(AASequence("B"));
hits[4].setScore(25);
hits[5].setRank(6);
hits[5].setSequence(AASequence("E"));
hits[5].setScore(5);
hits[6].setRank(7);
hits[6].setSequence(AASequence("H"));
hits[6].setScore(4);
hits[7].setRank(8);
hits[7].setSequence(AASequence("I"));
hits[7].setScore(3);
hits[8].setRank(9);
hits[8].setSequence(AASequence("J"));
hits[8].setScore(2);
hits[9].setRank(10);
hits[9].setSequence(AASequence("K"));
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

	//~ // ***** Merge ********

	//~ //define parameters
	//~ param.clear();
	//~ param.setValue("algorithm","merge");
	//~ param.setValue("considered_hits",6);
	//~ consensus.setParameters(param);
	//~ //apply
	//~ f = ids;
	//~ consensus.apply(f);

	//~ TEST_EQUAL(f.size(),1);
	//~ hits = f[0].getHits();
	//~ TEST_EQUAL(hits.size(),7);

	//~ TEST_EQUAL(hits[0].getRank(),1);
	//~ TEST_EQUAL(hits[0].getSequence(),"F");
	//~ TEST_REAL_SIMILAR(hits[0].getScore(),81.0f);

	//~ TEST_EQUAL(hits[1].getRank(),2);
	//~ TEST_EQUAL(hits[1].getSequence(),"C");
	//~ TEST_REAL_SIMILAR(hits[1].getScore(),60.0f);

	//~ TEST_EQUAL(hits[2].getRank(),3);
	//~ TEST_EQUAL(hits[2].getSequence(),"G");
	//~ TEST_REAL_SIMILAR(hits[2].getScore(),50.0f);

	//~ TEST_EQUAL(hits[3].getRank(),4);
	//~ TEST_EQUAL(hits[3].getSequence(),"D");
	//~ TEST_REAL_SIMILAR(hits[3].getScore(),40.0f);

	//~ TEST_EQUAL(hits[4].getRank(),5);
	//~ TEST_EQUAL(hits[4].getSequence(),"A");
	//~ TEST_REAL_SIMILAR(hits[4].getScore(),31.0f);

	//~ TEST_EQUAL(hits[5].getRank(),6);
	//~ TEST_EQUAL(hits[5].getSequence(),"B");
	//~ TEST_REAL_SIMILAR(hits[5].getScore(),29.0f);

	//~ TEST_EQUAL(hits[6].getRank(),7);
	//~ TEST_EQUAL(hits[6].getSequence(),"E");
	//~ TEST_REAL_SIMILAR(hits[6].getScore(),5.0f);

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
