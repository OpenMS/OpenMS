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
#include <OpenMS/DATASTRUCTURES/CVReference.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVReference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVReference* ptr = nullptr;
CVReference* nullPointer = nullptr;
START_SECTION(CVReference())
{
	ptr = new CVReference();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVReference())
{
	delete ptr;
}
END_SECTION

ptr = new CVReference();

START_SECTION((CVReference(const CVReference &rhs)))
{
 	CVReference cvr;
	TEST_STRING_EQUAL(CVReference(cvr).getName(), cvr.getName())
	TEST_STRING_EQUAL(CVReference(cvr).getIdentifier(), cvr.getIdentifier())

	cvr.setName("my_test_name");
	TEST_STRING_EQUAL(CVReference(cvr).getName(), "my_test_name")

	cvr.setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(CVReference(cvr).getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((CVReference& operator=(const CVReference &rhs)))
{
  CVReference cvr, cvr_copy;
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getName(), "")
	TEST_STRING_EQUAL(cvr_copy.getIdentifier(), "")

	cvr.setName("my_test_name");
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getName(), "my_test_name")

	cvr.setIdentifier("my_test_identifier");
	cvr_copy = cvr;
	TEST_STRING_EQUAL(cvr_copy.getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((bool operator == (const CVReference& rhs) const))
{
  CVReference cvr, cvr_copy;
	TEST_EQUAL(cvr == cvr_copy, true)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)

  cvr.setName("my_test_name");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)

  cvr.setIdentifier("my_test_identifier");
	TEST_EQUAL(cvr == cvr_copy, false)
  cvr_copy = cvr;
	TEST_EQUAL(cvr == cvr_copy, true)
}
END_SECTION

START_SECTION((bool operator != (const CVReference& rhs) const))
{
  CVReference cvr, cvr_copy;
  TEST_EQUAL(cvr != cvr_copy, false)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setName("my_test_name");
  TEST_EQUAL(cvr != cvr_copy, true)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)

  cvr.setIdentifier("my_test_identifier");
  TEST_EQUAL(cvr != cvr_copy, true)
  cvr_copy = cvr;
  TEST_EQUAL(cvr != cvr_copy, false)
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  ptr->setName("my_test_name");
	TEST_STRING_EQUAL(ptr->getName(), "my_test_name")
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIdentifier(const String &identifier)))
{
  ptr->setIdentifier("my_test_identifier");
	TEST_STRING_EQUAL(ptr->getIdentifier(), "my_test_identifier")
}
END_SECTION

START_SECTION((const String& getIdentifier() const ))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



