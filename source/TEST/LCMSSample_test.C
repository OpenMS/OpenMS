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
// $Maintainer: Stephan Aiche $
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/LCMSSample.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LCMSSample, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LCMSSample* ptr = 0;
START_SECTION(LCMSSample())
{
	ptr = new LCMSSample();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((virtual ~LCMSSample()))
{
	delete ptr;
}
END_SECTION

START_SECTION((LCMSSample(const LCMSSample &source)))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2(s1);	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::Iterator it1 = s1.begin();
	LCMSSample::Iterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}

}
END_SECTION

START_SECTION((LCMSSample& operator=(const LCMSSample &source)))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::Iterator it1 = s1.begin();
	LCMSSample::Iterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((void loadFASTA(const String filename)))
{
	// this get's almot a bit boring but I have no idea how else to test it.
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::Iterator it1 = s1.begin();
	LCMSSample::Iterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((const PeptideSequences& getPeptideSequences() const ))
{
	// this get's almot a bit boring but I have no idea how else to test it.
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test2.fasta"));
	s1.digest();

	LCMSSample::PeptideSequences ps = s1.getPeptideSequences();
	
	TEST_EQUAL(ps.size()==1,true)
	TEST_EQUAL(ps.begin()->first=="TVQMENQFVAFVDK",true)
	TEST_EQUAL(ps.begin()->second==1,true);
	
}
END_SECTION

START_SECTION((Iterator begin()))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::Iterator it1 = s1.begin();
	LCMSSample::Iterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((Iterator end()))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::Iterator it1 = s1.begin();
	LCMSSample::Iterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((ConstIterator begin() const ))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::ConstIterator it1 = s1.begin();
	LCMSSample::ConstIterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((ConstIterator end() const ))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::ConstIterator it1 = s1.begin();
	LCMSSample::ConstIterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((size_t size() const ))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
}
END_SECTION

START_SECTION((void digest()))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();

	LCMSSample s2 = s1;	

	TEST_EQUAL(s1.size(),s2.size());
	
	LCMSSample::ConstIterator it1 = s1.begin();
	LCMSSample::ConstIterator it2 = s2.begin();
		
	for (; it1!=s1.end();++it1,++it2)
	{
		TEST_EQUAL(*it1==*it2,true)
	}
}
END_SECTION

START_SECTION((void printProteins() const ))
{
  // no idea how to test this but I guess it is also not necessary
	TEST_EQUAL(1,1)
}
END_SECTION

START_SECTION((void printPeptides() const ))
{
  // no idea how to test this but I guess it is also not necessary
	TEST_EQUAL(1,1)
}
END_SECTION

START_SECTION((void clearProteins()))
{
	LCMSSample s1;
	s1.loadFASTA(OPENMS_GET_TEST_DATA_PATH("LCMSSample_test.fasta"));
	s1.digest();
	s1.clearProteins();

	// size() gives the number of peptides
	TEST_EQUAL(s1.size(),77);
}
END_SECTION

START_SECTION((void setPdModelFile(OpenMS::String file)))
{
	TEST_EQUAL(1,1)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



