// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPFile* ptr = nullptr;
MSPFile* nullPointer = nullptr;
START_SECTION((MSPFile()))
	ptr = new MSPFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MSPFile()))
	delete ptr;
END_SECTION

START_SECTION(MSPFile(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	MSPFile f3(f1);
	TEST_EQUAL(f1.getParameters() == f3.getParameters(), true)
END_SECTION

START_SECTION(MSPFile& operator=(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	f2 = f1;
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), true)
END_SECTION

START_SECTION(void load(const String &filename, std::vector< PeptideIdentification > &ids, PeakMap &exp))
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
	PeakMap exp;
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(ids.size(), 5)

		//test DocumentIdentifier addition
	TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"));
  TEST_STRING_EQUAL(FileTypes::typeToName(exp.getLoadedFileType()),"msp");


	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[3].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[4].getNativeID(), "index=4")

	Param p(msp_file.getParameters());
	p.setValue("instrument", "qtof");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear(true);
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 2)
	TEST_EQUAL(ids.size(), 2)

	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")

	p.setValue("instrument", "it");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear(true);
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 3)
	TEST_EQUAL(ids.size(), 3)

	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=4")

END_SECTION

START_SECTION(void store(const String& filename, const PeakMap& exp) const)
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
  PeakMap exp;
  msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	for (Size i = 0; i != ids.size(); ++i)
	{
		exp[i].getPeptideIdentifications().push_back(ids[i]);
	}
	String filename;
	NEW_TMP_FILE(filename)
	msp_file.store(filename, exp);

	exp.clear(true);
	ids.clear();
	msp_file.load(filename, ids, exp);
	TEST_EQUAL(ids.size(), 5)
	TEST_EQUAL(exp.size(), 5)

	TEST_EQUAL(ids[0].getHits().size(), 1)
	TEST_EQUAL(ids[1].getHits().size(), 1)
	TEST_EQUAL(ids[2].getHits().size(), 1)
	TEST_EQUAL(ids[3].getHits().size(), 1)
	TEST_EQUAL(ids[4].getHits().size(), 1)
	TEST_EQUAL(ids[0].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[1].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[2].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[3].getHits().begin()->getSequence().isModified(), true)
	TEST_EQUAL(ids[4].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[0].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[1].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[2].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[3].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[4].getHits().begin()->getCharge(), 3)
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
