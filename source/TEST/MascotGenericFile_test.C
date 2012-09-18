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
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <sstream>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MascotGenericFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MascotGenericFile* ptr = 0;
MascotGenericFile* nullPointer = 0;
START_SECTION(MascotGenericFile())
{
	ptr = new MascotGenericFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MascotGenericFile())
{
	delete ptr;
}
END_SECTION

ptr = new MascotGenericFile();

START_SECTION((template < typename MapType > void load(const String &filename, MapType &exp)))
{
  PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
  TEST_EQUAL(exp.size(), 1)

  TEST_EQUAL(exp.begin()->size(), 9)
}
END_SECTION

START_SECTION((void store(std::ostream &os, const String &filename, const PeakMap &experiment)))
{
	PeakMap exp;
	ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);
	
	stringstream ss;
	ptr->store(ss, "bla", exp);

	// BEGIN IONS
	// TITLE=Testtitle
	// PEPMASS=1998
	// RTINSECONDS=25.379
	// 1 1
	// 2 4
	// 3 9
	// 4 16
	// 5 25
	// 6 36
	// 7 49
	// 8 64
	// 9 81
	// END IONS
	vector<String> strings;
	strings.push_back("BEGIN IONS");
	strings.push_back("TITLE=Testtitle");
	strings.push_back("PEPMASS=1998");
	strings.push_back("RTINSECONDS=25.37");
	strings.push_back("1 1");
	strings.push_back("2 4");
	strings.push_back("3 9");
	strings.push_back("4 16");
	strings.push_back("5 25");
	strings.push_back("6 36");
	strings.push_back("7 49");
	strings.push_back("8 64");
	strings.push_back("9 81");
	strings.push_back("END IONS");

	String mgf_file(ss.str());
	for (Size i = 0; i != strings.size(); ++i)
	{
		TEST_EQUAL(mgf_file.hasSubstring(mgf_file[i]), true)
	}	
}
END_SECTION

START_SECTION((void store(const String &filename, const PeakMap &experiment)))
{
  String tmp_name("MascotGenericFile_1.tmp");
	NEW_TMP_FILE(tmp_name)
	PeakMap exp;
  ptr->load(OPENMS_GET_TEST_DATA_PATH("MascotInfile_test.mascot_in"), exp);


	ptr->store(tmp_name, exp);

	PeakMap exp2;
	ptr->load(tmp_name, exp2);
	TEST_EQUAL(exp.size() == exp2.size(), true)
	TEST_EQUAL(exp.begin()->size() == exp2.begin()->size(), true)
	TEST_REAL_SIMILAR(exp.begin()->getRT(), exp2.begin()->getRT())
	TEST_REAL_SIMILAR(exp.begin()->getPrecursors().begin()->getMZ(), exp2.begin()->getPrecursors().begin()->getMZ())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



