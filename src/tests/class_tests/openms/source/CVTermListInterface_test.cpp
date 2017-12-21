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
#include <OpenMS/METADATA/CVTermListInterface.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVTermListInterface, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVTermListInterface* ptr = nullptr;
CVTermListInterface* nullPointer = nullptr;
START_SECTION(CVTermListInterface())
{
	ptr = new CVTermListInterface();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVTermListInterface())
{
	delete ptr;
}
END_SECTION

START_SECTION((bool operator==(const CVTermListInterface &cv_term_list) const ))
{
  CVTermListInterface cv_term_list, cv_term_list2;
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
	cv_term_list.setMetaValue("blubb", "blubber");
	TEST_EQUAL(cv_term_list == cv_term_list2, false)
	cv_term_list2.setMetaValue("blubb", "blubber");
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
	CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
	CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	cv_term_list.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list == cv_term_list2, false)
	cv_term_list2.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list == cv_term_list2, true)
}
END_SECTION

START_SECTION((bool operator!=(const CVTermListInterface &cv_term_list) const ))
{
  CVTermListInterface cv_term_list, cv_term_list2;
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  cv_term_list.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.setMetaValue("blubb", "blubber");
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, false)
  cv_term_list2.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
}
END_SECTION

START_SECTION((bool hasCVTerm(const String &accession) const ))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
	CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	CVTermListInterface cv_term_list;
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
	cv_term_list.addCVTerm(cv_term);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
}
END_SECTION

/*
START_SECTION((bool checkCVTerms(const CVMappingRule &rule, const ControlledVocabulary &cv) const ))
{
}
END_SECTION
*/

START_SECTION((void setCVTerms(const std::vector< CVTerm > &terms)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
	CVTerm cv_term2("my_accession2", "my_name2", "my_cv_identifier_ref2", "4.0", unit);
	CVTermListInterface cv_term_list;
	vector<CVTerm> cv_terms;
	cv_terms.push_back(cv_term);
	cv_terms.push_back(cv_term2);
	cv_term_list.setCVTerms(cv_terms);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true);
	TEST_EQUAL(cv_term_list.hasCVTerm("my_accession2"), true);
}
END_SECTION

START_SECTION((const Map<String, std::vector<CVTerm> >& getCVTerms() const ))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession2", "my_name2", "my_cv_identifier_ref2", "4.0", unit);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.getCVTerms().size(), 0); // test empty one
  CVTermListInterface cv_term_list2;
  TEST_EQUAL(cv_term_list2.getCVTerms().size(), 0); // test empty one

  vector<CVTerm> cv_terms;
  cv_terms.push_back(cv_term);
  cv_terms.push_back(cv_term2);
  cv_term_list.setCVTerms(cv_terms);
  TEST_EQUAL(cv_term_list.getCVTerms().size(), 2);
  TEST_EQUAL(cv_term_list.getCVTerms().has("my_accession"), true);
  TEST_EQUAL(cv_term_list.getCVTerms().has("my_accession2"), true);
}
END_SECTION

START_SECTION((void addCVTerm(const CVTerm &term)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
}
END_SECTION

START_SECTION((void replaceCVTerm(const CVTerm &cv_term)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"].size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"][0].getValue(), "3.0")
  CVTerm cv_term2("my_accession", "my_name", "my_cv_identifier_ref", "2.0", unit);
  cv_term_list.replaceCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"].size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"][0].getValue(), "2.0")
}
END_SECTION

START_SECTION((void replaceCVTerms(const std::vector<CVTerm> &cv_terms)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession", "my_name", "my_cv_identifier_ref", "2.0", unit);
  std::vector<CVTerm> tmp;
  tmp.push_back(cv_term);
  tmp.push_back(cv_term2);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerms(tmp, "my_accession");
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"].size(), 2)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"][0].getValue(), "3.0")
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"][1].getValue(), "2.0")
  cv_term_list.replaceCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"].size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"][0].getValue(), "2.0")
}
END_SECTION

START_SECTION((void replaceCVTerms(const Map<String, vector<CVTerm> >& cv_term_map)))
{
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  CVTerm cv_term2("my_accession2", "my_name", "my_cv_identifier_ref", "2.0", unit);
  std::vector<CVTerm> tmp;
  tmp.push_back(cv_term);
  std::vector<CVTerm> tmp2;
  tmp2.push_back(cv_term2);
  Map<String, std::vector<CVTerm> >new_terms;
  new_terms["my_accession2"] = tmp2;
  TEST_EQUAL(new_terms.has("my_accession2"), true);

  // create CVTermListInterface with old "my_accession"
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.replaceCVTerms(tmp, "my_accession");
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession"].size(), 1)

  // replace the terms, delete "my_accession" and introduce "my_accession2"
  cv_term_list.replaceCVTerms(new_terms);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession2"), true)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession2"].size(), 1)
  TEST_EQUAL(cv_term_list.getCVTerms()["my_accession2"][0].getValue(), "2.0")
}
END_SECTION

/*
START_SECTION(([EXTRA] bool checkCVTerms(const ControlledVocabulary &cv) const ))
{
// [Term]
// id: MS:1000132
// name: percent of base peak
// def: "The magnitude of a peak or measurement element expressed in terms of the percentage of the magnitude of the base peak intensity." [PSI:MS]
// xref: value-type:xsd\:float "The allowed value-type for this CV term."
// is_a: MS:1000043 ! intensity unit
  CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)

	ControlledVocabulary cv;
	cv.loadFromOBO("MS", "CV/psi-ms.obo");

	TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)

  CVTerm cv_term2("MS:1000132", "percent of base peaks wrong", "MS", "3.0", unit);
	cv_term_list.addCVTerm(cv_term2);
	TEST_EQUAL(cv_term_list.checkCVTerms(cv), false)
}
END_SECTION

START_SECTION(([EXTRA] void correctCVTermNames()))
{
  CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermListInterface cv_term_list;
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), true)

  ControlledVocabulary cv;
  cv.loadFromOBO("MS", "CV/psi-ms.obo");

  TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)

  CVTerm cv_term2("MS:1000132", "percent of base peaks wrong", "MS", "3.0", unit);
  cv_term_list.addCVTerm(cv_term2);
  TEST_EQUAL(cv_term_list.checkCVTerms(cv), false)

	cv_term_list.correctCVTermNames();
	TEST_EQUAL(cv_term_list.checkCVTerms(cv), true)	
}
END_SECTION
*/

START_SECTION((CVTermListInterface(const CVTermListInterface &rhs)))
{
  CVTermListInterface cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
	CVTermListInterface cv_term_list2(cv_term_list);
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermListInterface cv_term_list3(cv_term_list);
  TEST_EQUAL(cv_term_list == cv_term_list3, true)
}
END_SECTION

START_SECTION((CVTermListInterface& operator=(const CVTermListInterface &rhs)))
{
  CVTermListInterface cv_term_list;
  cv_term_list.setMetaValue("blubb", "blubber");
  CVTermListInterface cv_term_list2;
	cv_term_list2 = cv_term_list;
  TEST_EQUAL(cv_term_list == cv_term_list2, true)
  CVTerm::Unit unit("my_unit_accession", "my_unit_name", "my_unit_ontology_name");
  CVTerm cv_term("my_accession", "my_name", "my_cv_identifier_ref", "3.0", unit);
  cv_term_list.addCVTerm(cv_term);
  CVTermListInterface cv_term_list3;
	cv_term_list3 = cv_term_list;
  TEST_EQUAL(cv_term_list == cv_term_list3, true)
}
END_SECTION

START_SECTION((bool empty() const))
{
	CVTerm::Unit unit("MS:1000043", "intensity unit", "MS");
  CVTerm cv_term("MS:1000132", "percent of base peak", "MS", "3.0", unit);
  CVTermListInterface cv_term_list;
	TEST_EQUAL(cv_term_list.empty(), true)
  TEST_EQUAL(cv_term_list.hasCVTerm("my_accession"), false)
  cv_term_list.addCVTerm(cv_term);
  TEST_EQUAL(cv_term_list.hasCVTerm("MS:1000132"), true)
	TEST_EQUAL(cv_term_list.empty(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



