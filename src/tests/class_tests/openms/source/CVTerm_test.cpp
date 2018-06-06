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
#include <OpenMS/METADATA/CVTerm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVTerm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVTerm* ptr = nullptr;
CVTerm* nullPointer = nullptr;
START_SECTION(CVTerm())
{
	ptr = new CVTerm();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVTerm())
{
	delete ptr;
}
END_SECTION

START_SECTION((bool operator==(const CVTerm &rhs) const ))
{
	CVTerm term1, term2;
	TEST_EQUAL(term1 == term2, true)
	
	term1.setAccession("acc");
	TEST_EQUAL(term1 == term2, false)
	term2.setAccession("acc");
	TEST_EQUAL(term1 == term2, true)

	term1.setName("name");
	TEST_EQUAL(term1 == term2, false)
	term2.setName("name");
	TEST_EQUAL(term1 == term2, true)

	term1.setCVIdentifierRef("cv_id_ref");
	TEST_EQUAL(term1 == term2, false)
	term2.setCVIdentifierRef("cv_id_ref");
	TEST_EQUAL(term1 == term2, true)

	term1.setValue(DataValue(0.4));
	TEST_EQUAL(term1 == term2, false)
	term2.setValue(DataValue(0.4));
	TEST_EQUAL(term1 == term2, true)

	term1.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
	TEST_EQUAL(term1 == term2, false)
	term2.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
	TEST_EQUAL(term1 == term2, true)
}
END_SECTION

START_SECTION((bool operator!=(const CVTerm &rhs) const ))
{
	CVTerm term1, term2;
  TEST_EQUAL(term1 != term2, false)

  term1.setAccession("acc");
  TEST_EQUAL(term1 != term2, true)
  term2.setAccession("acc");
  TEST_EQUAL(term1 != term2, false)

  term1.setName("name");
  TEST_EQUAL(term1 != term2, true)
  term2.setName("name");
  TEST_EQUAL(term1 != term2, false)

  term1.setCVIdentifierRef("cv_id_ref");
  TEST_EQUAL(term1 != term2, true)
  term2.setCVIdentifierRef("cv_id_ref");
  TEST_EQUAL(term1 != term2, false)

  term1.setValue(DataValue(0.4));
  TEST_EQUAL(term1 != term2, true)
  term2.setValue(DataValue(0.4));
  TEST_EQUAL(term1 != term2, false)

  term1.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
  TEST_EQUAL(term1 != term2, true)
  term2.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
  TEST_EQUAL(term1 != term2, false)
}
END_SECTION

START_SECTION((bool hasValue() const ))
{
  CVTerm term;
	TEST_EQUAL(term.hasValue(), false)
	term.setValue(DataValue(0.5));
	TEST_EQUAL(term.hasValue(), true)
}
END_SECTION

START_SECTION((bool hasUnit() const ))
{
	CVTerm term;
	TEST_EQUAL(term.hasUnit(), false)
	term.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
	TEST_EQUAL(term.hasUnit(), true)
}
END_SECTION

START_SECTION((void setAccession(const String &accession)))
{
	CVTerm term;
	TEST_STRING_EQUAL(term.getAccession(), "")
	term.setAccession("acc");
	TEST_STRING_EQUAL(term.getAccession(), "acc")
}
END_SECTION

START_SECTION((const String& getAccession() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  CVTerm term;
	TEST_STRING_EQUAL(term.getName(), "")
	term.setName("name");
	TEST_STRING_EQUAL(term.getName(), "name")
}
END_SECTION

START_SECTION((const String& getName() const ))
{
	NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void setCVIdentifierRef(const String &cv_identifier_ref)))
{
  CVTerm term;
	TEST_STRING_EQUAL(term.getCVIdentifierRef(), "")
	term.setCVIdentifierRef("cv_id_ref");
	TEST_STRING_EQUAL(term.getCVIdentifierRef(), "cv_id_ref")
}
END_SECTION

START_SECTION((const String& getCVIdentifierRef() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void setValue(const DataValue &value)))
{
  CVTerm term;
	TEST_EQUAL(term.getValue() == DataValue::EMPTY, true)
	DataValue value(300.0);
	term.setValue(value);
	TEST_REAL_SIMILAR((double)term.getValue(), 300.0)
	DataValue value2("bla");
	term.setValue(value2);
	TEST_STRING_EQUAL(term.getValue().toString(), "bla")
}
END_SECTION

START_SECTION((const DataValue& getValue() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((void setUnit(const Unit &unit)))
{
	CVTerm term;
  TEST_STRING_EQUAL(term.getUnit().accession, "")
	TEST_STRING_EQUAL(term.getUnit().name, "")
	TEST_STRING_EQUAL(term.getUnit().cv_ref, "")
	CVTerm::Unit unit("u_acc", "u_name", "u_cv_ref");
	term.setUnit(unit);
  TEST_STRING_EQUAL(term.getUnit().accession, "u_acc")
  TEST_STRING_EQUAL(term.getUnit().name, "u_name")
  TEST_STRING_EQUAL(term.getUnit().cv_ref, "u_cv_ref")
}
END_SECTION

START_SECTION((const Unit& getUnit() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((CVTerm(const String &accession, const String &name, const String &cv_identifier_ref, const String &value, const Unit &unit)))
{
	CVTerm::Unit unit("u_acc", "u_name", "u_cv_ref");
  CVTerm term("acc", "name", "cv_id_ref", "value", unit);
	TEST_STRING_EQUAL(term.getAccession(), "acc")
	TEST_STRING_EQUAL(term.getName(), "name")
	TEST_STRING_EQUAL(term.getCVIdentifierRef(), "cv_id_ref")
	TEST_STRING_EQUAL(term.getValue(), "value")
	TEST_STRING_EQUAL(term.getUnit().accession, "u_acc")
	TEST_STRING_EQUAL(term.getUnit().name, "u_name")
	TEST_STRING_EQUAL(term.getUnit().cv_ref, "u_cv_ref")
}
END_SECTION

START_SECTION((CVTerm(const CVTerm &rhs)))
{
  CVTerm term1;
  term1.setAccession("acc");
  term1.setName("name");
  term1.setCVIdentifierRef("cv_id_ref");
  term1.setValue(DataValue(0.4));
  term1.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));

	CVTerm term2(term1);
	TEST_STRING_EQUAL(term2.getAccession(), "acc")
	TEST_STRING_EQUAL(term2.getName(), "name")
	TEST_STRING_EQUAL(term2.getCVIdentifierRef(), "cv_id_ref")
	TEST_EQUAL(term2.getValue() == DataValue(0.4), true)
	TEST_EQUAL(term2.getUnit() == CVTerm::Unit("u_acc", "u_name", "u_cv_ref"), true)
}
END_SECTION

START_SECTION((CVTerm& operator=(const CVTerm &rhs)))
{
  CVTerm term1, term2;
  TEST_EQUAL(term1 == term2, true)
  
  term1.setAccession("acc");
  TEST_EQUAL(term1 == term2, false)
	term2 = term1;
  TEST_EQUAL(term1 == term2, true)
  
  term1.setName("name");
  TEST_EQUAL(term1 == term2, false)
	term2 = term1;
  TEST_EQUAL(term1 == term2, true)
  
  term1.setCVIdentifierRef("cv_id_ref");
  TEST_EQUAL(term1 == term2, false)
	term2 = term1;
  TEST_EQUAL(term1 == term2, true)
  
  term1.setValue(DataValue(0.4));
  TEST_EQUAL(term1 == term2, false)
	term2 = term1;
  TEST_EQUAL(term1 == term2, true)
  
  term1.setUnit(CVTerm::Unit("u_acc", "u_name", "u_cv_ref"));
  TEST_EQUAL(term1 == term2, false)
	term2 = term1;
  TEST_EQUAL(term1 == term2, true)
}
END_SECTION

CVTerm::Unit* ptr_unit = nullptr;
CVTerm::Unit* nullPointer_unit = nullptr;
START_SECTION(([CVTerm::Unit] Unit()))
{
  ptr_unit = new CVTerm::Unit();
  TEST_NOT_EQUAL(ptr_unit, nullPointer_unit)
}
END_SECTION

START_SECTION(([CVTerm::Unit] Unit(const String &p_accession, const String &p_name, const String &p_cv_ref)))
{
  CVTerm::Unit u("ACCESSION", "p_name", "p_cv_ref");
  TEST_EQUAL(u.accession, "ACCESSION")
  TEST_EQUAL(u.cv_ref, "p_cv_ref")
  TEST_EQUAL(u.name, "p_name")
}
END_SECTION

START_SECTION(([CVTerm::Unit] Unit(const Unit &rhs)))
{
  CVTerm::Unit u("ACCESSION", "p_name", "p_cv_ref");
  TEST_EQUAL(u.accession, "ACCESSION")
  TEST_EQUAL(u.cv_ref, "p_cv_ref")
  TEST_EQUAL(u.name, "p_name")

  CVTerm::Unit cu(u);
  TEST_EQUAL(cu.accession, u.accession)
  TEST_EQUAL(cu.cv_ref, u.cv_ref)
  TEST_EQUAL(cu.name, u.name)
}
END_SECTION

START_SECTION(([CVTerm::Unit] virtual ~Unit()))
{
  delete ptr_unit;
}
END_SECTION

START_SECTION(([CVTerm::Unit] Unit& operator=(const Unit &rhs)))
{
  CVTerm::Unit u("ACCESSION", "p_name", "p_cv_ref");
  TEST_EQUAL(u.accession, "ACCESSION")
  TEST_EQUAL(u.cv_ref, "p_cv_ref")
  TEST_EQUAL(u.name, "p_name")

  CVTerm::Unit cu;
  cu = u;
  TEST_EQUAL(cu.accession, u.accession)
  TEST_EQUAL(cu.cv_ref, u.cv_ref)
  TEST_EQUAL(cu.name, u.name)
}
END_SECTION

START_SECTION(([CVTerm::Unit] bool operator==(const Unit &rhs) const ))
{
  CVTerm::Unit u("ACCESSION", "p_name", "p_cv_ref");
  CVTerm::Unit cu("ACCESSION", "p_name", "p_cv_ref");
  CVTerm::Unit nu("ACCESSION2", "p_name", "p_cv_ref");
  CVTerm::Unit nu2("ACCESSION", "p_name2", "p_cv_ref");
  CVTerm::Unit nu3("ACCESSION", "p_name", "p_cv_ref2");

  TEST_EQUAL(u==cu, true)
  TEST_EQUAL(u==u, true)
  TEST_EQUAL(u==nu, false)
  TEST_EQUAL(u==nu2, false)
  TEST_EQUAL(u==nu3, false)
  TEST_EQUAL(cu==nu, false)
}
END_SECTION

START_SECTION(([CVTerm::Unit] bool operator!=(const Unit &rhs) const ))
{
  CVTerm::Unit u("ACCESSION", "p_name", "p_cv_ref");
  CVTerm::Unit cu("ACCESSION", "p_name", "p_cv_ref");
  CVTerm::Unit nu("ACCESSION2", "p_name", "p_cv_ref");
  CVTerm::Unit nu2("ACCESSION", "p_name2", "p_cv_ref");
  CVTerm::Unit nu3("ACCESSION", "p_name", "p_cv_ref2");

  TEST_EQUAL(u!=cu, false)
  TEST_EQUAL(u!=u, false)
  TEST_EQUAL(u!=nu, true)
  TEST_EQUAL(u!=nu2, true)
  TEST_EQUAL(u!=nu3, true)
  TEST_EQUAL(cu!=nu, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



