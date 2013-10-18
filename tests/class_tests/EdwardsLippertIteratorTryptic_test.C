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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
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
EdwardsLippertIteratorTryptic* nullPointer = 0;
START_SECTION(EdwardsLippertIteratorTryptic())
	ptr = new EdwardsLippertIteratorTryptic();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~EdwardsLippertIteratorTryptic())
	delete ptr;
END_SECTION

START_SECTION(EdwardsLippertIteratorTryptic(const EdwardsLippertIteratorTryptic& rhs))
  ptr = new EdwardsLippertIteratorTryptic();
  ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("EdwardsLippertIterator_test_2.fasta"));
	vector<DoubleReal> specc;
	specc.push_back(178.1864);
	specc.push_back(441.4806);
  ptr->setSpectrum(specc);
  ptr->begin();
  ++*ptr;
  EdwardsLippertIteratorTryptic copy (*ptr);
  TEST_EQUAL ((*ptr).getFastaFile(),(copy).getFastaFile());
  TEST_EQUAL ((*ptr).getTolerance(),(copy).getTolerance());
  TEST_EQUAL ((**ptr).first,(*copy).first);
  TEST_EQUAL ((**ptr).second,(*copy).second);
END_SECTION

START_SECTION(EdwardsLippertIteratorTryptic& operator=(const EdwardsLippertIteratorTryptic &rhs))
	ptr = new EdwardsLippertIteratorTryptic();
  ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("EdwardsLippertIterator_test_2.fasta"));
  vector<DoubleReal> specc;
  specc.push_back(178.1864);
  specc.push_back(441.4806);
  ptr->setSpectrum(specc);
  ptr->begin();
  ++*ptr;
  EdwardsLippertIteratorTryptic copy;
	copy = *ptr;
  TEST_EQUAL ((*ptr).getFastaFile(),(copy).getFastaFile());
  TEST_EQUAL ((*ptr).getTolerance(),(copy).getTolerance());
  TEST_EQUAL ((**ptr).first,(*copy).first);
  TEST_EQUAL ((**ptr).second,(*copy).second);				
END_SECTION

START_SECTION(virtual bool isDigestingEnd(char aa1,char aa2))
	ptr = new EdwardsLippertIteratorTryptic();
	TEST_EQUAL(ptr->isDigestingEnd('R','S'),true)
	TEST_EQUAL(ptr->isDigestingEnd('K','S'),true)
	TEST_EQUAL(ptr->isDigestingEnd('R','P'),false)
	TEST_EQUAL(ptr->isDigestingEnd('K','P'),false)
	TEST_EQUAL(ptr->isDigestingEnd('S','S'),false)
END_SECTION

START_SECTION([EXTRA] FASTAEntry operator*())
	vector<DoubleReal> spec;
	spec.push_back(178.1864);
	spec.push_back(441.4806);
	const vector<DoubleReal> specc (spec);

	ptr = new EdwardsLippertIteratorTryptic();
	FastaIterator * fit = new FastaIterator();
	ptr->setFastaFile(OPENMS_GET_TEST_DATA_PATH("EdwardsLippertIterator_test_2.fasta"));
	fit->setFastaFile(OPENMS_GET_TEST_DATA_PATH("EdwardsLippertIterator_test_2.fasta"));
	ptr->setSpectrum(specc);
	fit->begin();
	ptr->begin();
	
	DoubleReal tol = 0.2;
	ptr->setTolerance(tol);
	
	while (!ptr->isAtEnd())
	{
		while ((**fit).first != (**ptr).first) 
		{
			++*fit;
		}
		
		String seq = (**ptr).second;
		String realSeq = (**fit).second;
		
		bool isCorrect = false;
		if (seq == realSeq.substr(0,seq.length()))
		{
			if (realSeq[seq.length() - 1] == 'R' || realSeq[seq.length() - 1] == 'K') isCorrect = true;
		} 
		else if (seq == realSeq.substr(realSeq.length() - seq.length() - 1, seq.length()))
		{
			if (realSeq[realSeq.length() - seq.length()-2] == 'R' || realSeq[realSeq.length() - seq.length() - 2] == 'K') isCorrect = true;
		} 
		else 
		{
			for (Size i = 1; i < realSeq.length() - 1; ++i)
			{
				if (realSeq.substr(i, seq.length()) == seq)
				{
					if ((realSeq[i-1] == 'R' || realSeq[i-1] == 'K') && seq[0] != 'P' && (seq[seq.length() - 1] == 'R' || seq[seq.length() - 1] == 'K')) isCorrect = true;
				}
			}
		}
		if (realSeq.hasSuffix(seq))
		{
			String suffix = realSeq.suffix(seq.size() + 1);
			if ((suffix[0] == 'K' || suffix[0] == 'R') && seq[0] != 'P')
			{
				isCorrect = true;
			}
		}
		
		TEST_EQUAL (isCorrect, true);
		++*ptr;
	}
END_SECTION

START_SECTION(static const String getProductName())
	TEST_STRING_EQUAL(EdwardsLippertIteratorTryptic::getProductName(), "EdwardsLippertIteratorTryptic")
END_SECTION
			
START_SECTION(static PepIterator* create())
  TEST_NOT_EQUAL(EdwardsLippertIterator::create(), nullPointer)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



