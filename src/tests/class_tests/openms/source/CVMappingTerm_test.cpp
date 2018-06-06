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
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVMappingTerm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappingTerm* ptr = nullptr;
CVMappingTerm* nullPointer = nullptr;
START_SECTION(CVMappingTerm())
{
	ptr = new CVMappingTerm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVMappingTerm())
{
	delete ptr;
}
END_SECTION

ptr = new CVMappingTerm();

START_SECTION((CVMappingTerm(const CVMappingTerm &rhs)))
{
 	CVMappingTerm cvmt;
	cvmt.setAccession("my_test_accession");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getAccession(), "my_test_accession")

	cvmt.setUseTermName(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTermName(), true)
	cvmt.setUseTermName(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTermName(), false)

	cvmt.setUseTerm(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTerm(), true)
	cvmt.setUseTerm(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getUseTerm(), false)

	cvmt.setTermName("my_test_termname");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getTermName(), "my_test_termname")

	cvmt.setIsRepeatable(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getIsRepeatable(), true)
	cvmt.setIsRepeatable(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getIsRepeatable(), false)

	cvmt.setAllowChildren(true);
	TEST_EQUAL(CVMappingTerm(cvmt).getAllowChildren(), true)
	cvmt.setAllowChildren(false);
	TEST_EQUAL(CVMappingTerm(cvmt).getAllowChildren(), false)

	cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	TEST_STRING_EQUAL(CVMappingTerm(cvmt).getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION

START_SECTION((CVMappingTerm& operator=(const CVMappingTerm &rhs)))
{
	CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getAccession(), "my_test_accession")

  cvmt.setUseTermName(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTermName(), true)
  cvmt.setUseTermName(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTermName(), false)

  cvmt.setUseTerm(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTerm(), true)
  cvmt.setUseTerm(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getUseTerm(), false)

  cvmt.setTermName("my_test_termname");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getTermName(), "my_test_termname")

  cvmt.setIsRepeatable(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getIsRepeatable(), true)
  cvmt.setIsRepeatable(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getIsRepeatable(), false)

  cvmt.setAllowChildren(true);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getAllowChildren(), true)
  cvmt.setAllowChildren(false);
	cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy.getAllowChildren(), false)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	cvmt_copy = cvmt;
  TEST_STRING_EQUAL(cvmt_copy.getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION


START_SECTION((bool operator == (const CVMappingTerm& rhs) const))
{
  CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTermName(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTermName(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTerm(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTerm(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setTermName("my_test_termname");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  
	cvmt.setIsRepeatable(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setIsRepeatable(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setAllowChildren(true);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setAllowChildren(false);
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
	TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
	TEST_EQUAL(cvmt_copy == cvmt, true)
}
END_SECTION

START_SECTION((bool operator != (const CVMappingTerm& rhs) const))
{
  CVMappingTerm cvmt, cvmt_copy;

  cvmt.setAccession("my_test_accession");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTermName(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTermName(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setUseTerm(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setUseTerm(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setTermName("my_test_termname");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setIsRepeatable(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setIsRepeatable(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setAllowChildren(true);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
  cvmt.setAllowChildren(false);
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)

  cvmt.setCVIdentifierRef("my_test_cvidentifierref");
  TEST_EQUAL(cvmt_copy == cvmt, false)
  cvmt_copy = cvmt;
  TEST_EQUAL(cvmt_copy == cvmt, true)
}
END_SECTION

START_SECTION((void setAccession(const String &accession)))
{
  ptr->setAccession("my_test_accession");
	TEST_STRING_EQUAL(ptr->getAccession(), "my_test_accession")
}
END_SECTION

START_SECTION((const String& getAccession() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setUseTermName(bool use_term_name)))
{
  ptr->setUseTermName(true);
	TEST_EQUAL(ptr->getUseTermName(), true)
	ptr->setUseTermName(false);
	TEST_EQUAL(ptr->getUseTermName(), false)
}
END_SECTION

START_SECTION((bool getUseTermName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setUseTerm(bool use_term)))
{
  ptr->setUseTerm(true);
	TEST_EQUAL(ptr->getUseTerm(), true)
	ptr->setUseTerm(false);
	TEST_EQUAL(ptr->getUseTerm(), false)
}
END_SECTION

START_SECTION((bool getUseTerm() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setTermName(const String &term_name)))
{
  ptr->setTermName("my_test_termname");
	TEST_EQUAL(ptr->getTermName(), "my_test_termname")
}
END_SECTION

START_SECTION((const String& getTermName() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setIsRepeatable(bool is_repeatable)))
{
  ptr->setIsRepeatable(true);
	TEST_EQUAL(ptr->getIsRepeatable(), true)
	ptr->setIsRepeatable(false);
	TEST_EQUAL(ptr->getIsRepeatable(), false)
}
END_SECTION

START_SECTION((bool getIsRepeatable() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setAllowChildren(bool allow_children)))
{
  ptr->setAllowChildren(true);
	TEST_EQUAL(ptr->getAllowChildren(), true)
	ptr->setAllowChildren(false);
	TEST_EQUAL(ptr->getAllowChildren(), false)
}
END_SECTION

START_SECTION((bool getAllowChildren() const ))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void setCVIdentifierRef(const String &cv_identifier_ref)))
{
  ptr->setCVIdentifierRef("my_test_cvidentifierref");
	TEST_EQUAL(ptr->getCVIdentifierRef(), "my_test_cvidentifierref")
}
END_SECTION

START_SECTION((const String& getCVIdentifierRef() const ))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



