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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <cstdio>
#include <cstdlib>
///////////////////////////
#include <OpenMS/METADATA/DocumentIDTagger.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DocumentIDTagger, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DocumentIDTagger* ptr = nullptr;
DocumentIDTagger* nullPointer = nullptr;
START_SECTION(DocumentIDTagger())
{
	ptr = new DocumentIDTagger("someTOPPTool");
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DocumentIDTagger())
{
	delete ptr;
}
END_SECTION

START_SECTION((DocumentIDTagger(String toolname)))
{
	DocumentIDTagger tagme("SomeTOPPTool");
  NOT_TESTABLE
}
END_SECTION

START_SECTION((DocumentIDTagger(const DocumentIDTagger &source)))
{
  DocumentIDTagger tagme("SomeTOPPTool");
  DocumentIDTagger tagme2(tagme);

  TEST_EQUAL(tagme==tagme2, true)
}
END_SECTION

START_SECTION((DocumentIDTagger& operator=(const DocumentIDTagger &source)))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2 = tagme;
	TEST_EQUAL(tagme==tagme2,true)
}
END_SECTION

START_SECTION((bool operator==(const DocumentIDTagger &source) const ))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2 = tagme;
	TEST_EQUAL(tagme==tagme2, true)
	DocumentIDTagger tagme3(tagme);
	TEST_EQUAL(tagme==tagme3, true)
}
END_SECTION

START_SECTION((bool operator!=(const DocumentIDTagger &source) const ))
{
  DocumentIDTagger tagme("SomeTOPPTool");
	DocumentIDTagger tagme2("SomeOtherTOPPTool");
	TEST_EQUAL(tagme!=tagme2, true)
}
END_SECTION

START_SECTION((String getPoolFile() const))
{
	NOT_TESTABLE; // tested below
}
END_SECTION
			
START_SECTION((void setPoolFile(const String& file)))
{
	String tmp_pool;
	NEW_TMP_FILE(tmp_pool);
	DocumentIDTagger tagme("SomeTOPPTool");
	//use custom pool file
	tagme.setPoolFile(tmp_pool);
	TEST_EQUAL(tagme.getPoolFile(), tmp_pool)
}
END_SECTION			

String tmp_pool;
NEW_TMP_FILE(tmp_pool);
ofstream outfile;
outfile.open (tmp_pool.c_str(), ofstream::out);
outfile << "ID1\nIDNew\nIDsecondtoLast\nIDLast\n";
outfile.close();

START_SECTION((bool tag(DocumentIdentifier &map) const ))
{
	DocumentIdentifier myD;
	myD.setIdentifier("");
	DocumentIDTagger tagme("SomeTOPPTool");
	//use custom pool file
	tagme.setPoolFile(tmp_pool);
	Int cnt(0);
	tagme.countFreeIDs(cnt);
	TEST_EQUAL(cnt, 4);
	tagme.tag(myD);
	TEST_EQUAL(myD.getIdentifier(), "ID1");
	tagme.tag(myD);
	TEST_EQUAL(myD.getIdentifier(), "IDNew");
	tagme.countFreeIDs(cnt);
	TEST_EQUAL(cnt, 2);
	// 2 left
	TEST_EQUAL(tagme.tag(myD), true);
	// 1 left
	TEST_EQUAL(tagme.tag(myD), true);
	//0 left, expect it to go wrong
	TEST_EXCEPTION(Exception::DepletedIDPool, tagme.tag(myD));
	// confirm 0 left
	TEST_EQUAL(tagme.countFreeIDs(cnt), true);
	TEST_EQUAL(cnt, 0);
}
END_SECTION

START_SECTION((bool countFreeIDs(Int &free) const ))
{
  NOT_TESTABLE //done above
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



