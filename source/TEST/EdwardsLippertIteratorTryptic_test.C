// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// --------------------------------------------------------------------------

/*
The only member from EdwardsLippertIterator that is overwritten is the
isDigestingEnd method. So this is the only one that needs testing. 
Besides we test if every result is the product of a tryptic cleavage.
*/

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FastaIterator.h>
///////////////////////////
#include <OpenMS/CHEMISTRY/EdwardsLippertIteratorTryptic.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EdwardsLippertIteratorTryptic, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EdwardsLippertIteratorTryptic* ptr = 0;
CHECK(EdwardsLippertIteratorTryptic())
	ptr = new EdwardsLippertIteratorTryptic();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~EdwardsLippertIteratorTryptic())
	delete ptr;
RESULT

CHECK (virtual bool isDigestingEnd(char aa1,char aa2))
	ptr = new EdwardsLippertIteratorTryptic();
	TEST_EQUAL(ptr->isDigestingEnd('R','S'),1)
	TEST_EQUAL(ptr->isDigestingEnd('K','S'),1)
	TEST_EQUAL(ptr->isDigestingEnd('R','P'),0)
	TEST_EQUAL(ptr->isDigestingEnd('K','P'),0)
	TEST_EQUAL(ptr->isDigestingEnd('S','S'),0)
RESULT

CHECK (FASTAEntry operator*())
				/* TODO
	vector<float> spec;
	spec.push_back(178.1864);
	spec.push_back(441.4806);
	const vector<float> specc (spec);

	ptr = new EdwardsLippertIteratorTryptic();
	FastaIterator * fit = new FastaIterator();
	ptr->setFastaFile("data/EdwardsLippertIterator_test_2.fasta");
	fit->setFastaFile("data/EdwardsLippertIterator_test_2.fasta");
	ptr->setSpectrum(specc);
	fit->begin();
	ptr->begin();
	float tol = 0.2;
	ptr->setTolerance(tol);
	while (!ptr->isAtEnd())
	{
		while ((**fit).first != (**ptr).first) ++*fit;
		String seq = (**ptr).second;
		String realSeq = (**fit).second;
		bool isCorrect = false;
		if (seq == realSeq.substr(0,seq.length()))
		{
			if (realSeq[seq.length()-1]=='R'||realSeq[seq.length()-1]=='K') isCorrect = true;
		} 
		else if (seq == realSeq.substr(realSeq.length()-seq.length()-1,seq.length()))
		{
			if (realSeq[realSeq.length()-seq.length()-2]=='R'||realSeq[realSeq.length()-seq.length()-2]=='K') isCorrect = true;
		} 
		else 
		{
			for (unsigned int i = 1; i<realSeq.length()-1;++i)
			{
				if (realSeq.substr(i,seq.length())==seq)
				{
					if ((realSeq[i-1]=='R'||realSeq[i-1]=='K')&&seq[0]!='P' && (seq[seq.length()-1]=='R'||seq[seq.length()-1]=='K') ) isCorrect=true;
				}
			}
		}
		TEST_EQUAL (isCorrect,1);
		++*ptr;
	}*/
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



