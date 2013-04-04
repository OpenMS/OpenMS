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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(PILISIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISIdentification* ptr = 0;
PILISIdentification* nullPointer = 0;
RichPeakSpectrum spec;
DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), spec);
spec.setMSLevel(2);
START_SECTION(PILISIdentification())
	ptr = new PILISIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PILISIdentification())
	delete ptr;
END_SECTION

ptr = new PILISIdentification();

START_SECTION(PILISIdentification(const PILISIdentification& source))
	PILISIdentification copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(PILISIdentification& operator = (const PILISIdentification& source))
	PILISIdentification copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void setModel(PILISModel* hmm_model))
	PILISModel* model = new PILISModel();
	model->readFromFile("PILIS/PILIS_default_model.dat");
	ptr->setModel(model);
END_SECTION

START_SECTION(void getIdentification(const std::map<String, UInt>& candidates, PeptideIdentification& id, const RichPeakSpectrum& spectrum))
	map<String, UInt> candidates;
	candidates["DDFPIVIVGNKADIENQR"] = 2;
	candidates["DFPIANGER"] = 1;
	candidates["DFPIADGER"] = 1;
	PeptideIdentification id;
	ptr->getIdentification(candidates, id, spec);
	TEST_EQUAL(id.getHits().size(), 3)
	TEST_EQUAL(id.getHits().begin()->getSequence(), "DFPIANGER")
END_SECTION

START_SECTION(void getIdentifications(const std::vector<std::map<String, UInt> >& candidates, std::vector<PeptideIdentification>& ids, const RichPeakMap& exp))

	map<String, UInt> cand;
	cand["DDFPIVIVGNKADIENQR"] = 2;
	cand["DFPIANGER"] = 1;
	cand["DFPIADGER"] = 1;
	vector<map<String, UInt> > candidates;
	candidates.push_back(cand);

	vector<PeptideIdentification> ids;
	RichPeakMap map;
	map.push_back(spec);
	ptr->getIdentifications(candidates, ids, map);
	TEST_EQUAL(ids.size(), map.size())
	TEST_EQUAL(ids.begin()->getHits().size(), 3)
	TEST_EQUAL(ids.begin()->getHits().begin()->getSequence(), "DFPIANGER")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
