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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

CsvFile* ptr = nullptr;
CsvFile* nullPointer = nullptr;
START_SECTION(CsvFile())
	ptr = new CsvFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~CsvFile())
	delete ptr;
END_SECTION

#if 0

// Something is terribly wrong here ... looks like an unintentional commit?

START_SECTION(CsvFile(const String& filename, char is = ',',bool ie = false, Int first_n = -1))
//tested in getRow
TEST_EXCEPTION(Exception::FileNotFound, CsvFile("CsvFile_1.csv"))
END_SECTION

START_SECTION(void load(const String& filename, char is = ',', bool ie = false, Int first_n = -1))
//tested in getRow
TEST_EXCEPTION(Exception::FileNotFound, f1.load("CsvFile_1.csv"))


END_SECTION

#endif

START_SECTION(bool getRow(Size row,StringList &list))
	TOLERANCE_ABSOLUTE(0.01)
	CsvFile f1,f3,f4;

	CsvFile f2(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"), '\t');
	StringList list;
	f2.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f2.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f2.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

	f3.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),'\t');
	f3.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f3.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f3.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

	f4.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),'\t',true);
	f4.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f4.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f4.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

END_SECTION

START_SECTION(void store(const String& filename))
	CsvFile f1,f2;
	StringList list;

	f1.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),'\t',true); // load from a file
	f1.store(OPENMS_GET_TEST_DATA_PATH("CsvFile_3.csv"));          // store into a new one
	f2.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_3.csv"),'\t',true); // load the new one
	f2.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f2.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f2.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))
END_SECTION

START_SECTION(void addRow(const StringList& list))
	CsvFile f1, f2;
	StringList list;

	f1.addRow(ListUtils::create<String>("first,second,third"));
	f1.addRow(ListUtils::create<String>("4,5,6"));

	f1.store(OPENMS_GET_TEST_DATA_PATH("CsvFile_4.csv"));
	f2.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_4.csv"), ',', false);
	f2.getRow(0,list);
	TEST_EQUAL(list, ListUtils::create<String>("first,second,third"))
	f2.getRow(1,list);
	TEST_EQUAL(list, ListUtils::create<String>("4,5,6"))
END_SECTION

START_SECTION(void clear())
	CsvFile f1;
	StringList list;

	f1.addRow(ListUtils::create<String>("hello,world"));
	f1.getRow(0, list);
	TEST_EQUAL(list, ListUtils::create<String>("hello,world"))
	f1.clear();
	TEST_EXCEPTION(Exception::InvalidIterator, f1.getRow(0, list))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
