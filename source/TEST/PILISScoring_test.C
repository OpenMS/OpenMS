// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

///////////////////////////

START_TEST(PILISScoring_test.C, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISScoring* ptr = 0;
PILISScoring* nullPointer = 0;
String filename(OPENMS_GET_TEST_DATA_PATH("IDFilter_test2.idXML"));
START_SECTION(PILISScoring())
	ptr = new PILISScoring();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PILISScoring())
	delete ptr;
END_SECTION

ptr = new PILISScoring();

START_SECTION(PILISScoring(const PILISScoring& source))
	PILISScoring copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PILISScoring& operator = (const PILISScoring& source))
	PILISScoring copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void getScores(std::vector<PeptideIdentification>& ids))
	vector<PeptideIdentification> ids;
	vector<ProteinIdentification> prot_ids;
	String document_id;
	IdXMLFile().load(filename, prot_ids, ids, document_id);
	ptr->getScores(ids);
	for (vector<PeptideIdentification>::const_iterator it = ids.begin(); it != ids.end(); ++it)
	{
		TEST_EQUAL(it->getScoreType(), "PILIS-E-value")
	}
END_SECTION

START_SECTION(void getScore(PeptideIdentification& id))
	vector<PeptideIdentification> ids;
	vector<ProteinIdentification> prot_ids;
	String document_id;
	IdXMLFile().load(filename, prot_ids, ids, document_id);
	ptr->getScore(ids[0]);
	TEST_REAL_SIMILAR(ids[0].getHits().begin()->getScore(), 33.85)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
