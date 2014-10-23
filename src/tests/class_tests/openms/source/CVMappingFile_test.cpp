// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/CVMappingFile.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>

using namespace OpenMS;
using namespace std;

START_TEST(CVMappingFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappingFile* ptr = 0;
CVMappingFile* nullPointer = 0;
START_SECTION(CVMappingFile())
{
	ptr = new CVMappingFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CVMappingFile())
{
	delete ptr;
}
END_SECTION

START_SECTION((void load(const String &filename, CVMappings &cv_mappings, bool strip_namespaces=false)))
{
  CVMappings mappings;
	CVMappingFile().load(OPENMS_GET_TEST_DATA_PATH("cv_mapping_test_file.xml"), mappings);

	TEST_EQUAL(mappings.getMappingRules().size(), 9)

	vector<CVMappingRule> rules = mappings.getMappingRules();

	TEST_STRING_EQUAL(rules[0].getIdentifier(), "0")
	TEST_STRING_EQUAL(rules[1].getIdentifier(), "1")
	TEST_STRING_EQUAL(rules[2].getIdentifier(), "2")
	TEST_STRING_EQUAL(rules[3].getIdentifier(), "3")
	TEST_STRING_EQUAL(rules[4].getIdentifier(), "4")
	TEST_STRING_EQUAL(rules[5].getIdentifier(), "5")
	TEST_STRING_EQUAL(rules[6].getIdentifier(), "6")
	TEST_STRING_EQUAL(rules[7].getIdentifier(), "7")
	TEST_STRING_EQUAL(rules[8].getIdentifier(), "8")

	TEST_EQUAL(rules[0].getCVTerms().size(), 14)
	TEST_STRING_EQUAL(rules[0].getElementPath(), "/mzData/description/admin/sampleDescription/cvParam/@accession")
	TEST_EQUAL(rules[0].getRequirementLevel(), CVMappingRule::MAY)
	TEST_STRING_EQUAL(rules[0].getScopePath(), "/mzData/description/admin/sampleDescription")
	TEST_EQUAL(rules[0].getCombinationsLogic(), CVMappingRule::OR)

  TEST_EQUAL(rules[1].getCVTerms().size(), 32)
  TEST_STRING_EQUAL(rules[1].getElementPath(), "/mzData/description/instrument/source/cvParam/@accession")
  TEST_EQUAL(rules[1].getRequirementLevel(), CVMappingRule::SHOULD)
  TEST_STRING_EQUAL(rules[1].getScopePath(), "/mzData/description/instrument/source")
  TEST_EQUAL(rules[1].getCombinationsLogic(), CVMappingRule::XOR)

  TEST_EQUAL(rules[2].getCVTerms().size(), 46)
  TEST_STRING_EQUAL(rules[2].getElementPath(), "/mzData/description/instrument/analyzerList/analyzer/cvParam/@accession")
  TEST_EQUAL(rules[2].getRequirementLevel(), CVMappingRule::MUST)
  TEST_STRING_EQUAL(rules[2].getScopePath(), "/mzData/description/instrument/analyzerList/analyzer")
  TEST_EQUAL(rules[2].getCombinationsLogic(), CVMappingRule::AND)



	TEST_STRING_EQUAL(rules[0].getCVTerms()[0].getAccession(), "PSI:1000001")
	TEST_EQUAL(rules[0].getCVTerms()[0].getUseTermName(), false)
	TEST_EQUAL(rules[0].getCVTerms()[0].getUseTerm(), true)
	TEST_STRING_EQUAL(rules[0].getCVTerms()[0].getTermName(), "Sample Number")
	TEST_EQUAL(rules[0].getCVTerms()[0].getIsRepeatable(), true)
	TEST_EQUAL(rules[0].getCVTerms()[0].getAllowChildren(), true)
	TEST_STRING_EQUAL(rules[0].getCVTerms()[0].getCVIdentifierRef(), "PSI")

  TEST_STRING_EQUAL(rules[0].getCVTerms()[1].getAccession(), "PSI:1000002")
  TEST_EQUAL(rules[0].getCVTerms()[1].getUseTermName(), true)
  TEST_EQUAL(rules[0].getCVTerms()[1].getUseTerm(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[1].getTermName(), "Sample Name")
  TEST_EQUAL(rules[0].getCVTerms()[1].getIsRepeatable(), true)
  TEST_EQUAL(rules[0].getCVTerms()[1].getAllowChildren(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[1].getCVIdentifierRef(), "PSI")

  TEST_STRING_EQUAL(rules[0].getCVTerms()[2].getAccession(), "PSI:1000003")
  TEST_EQUAL(rules[0].getCVTerms()[2].getUseTermName(), true)
  TEST_EQUAL(rules[0].getCVTerms()[2].getUseTerm(), false)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[2].getTermName(), "Sample State")
  TEST_EQUAL(rules[0].getCVTerms()[2].getIsRepeatable(), true)
  TEST_EQUAL(rules[0].getCVTerms()[2].getAllowChildren(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[2].getCVIdentifierRef(), "PSI")

  TEST_STRING_EQUAL(rules[0].getCVTerms()[3].getAccession(), "PSI:1000004")
  TEST_EQUAL(rules[0].getCVTerms()[3].getUseTermName(), true)
  TEST_EQUAL(rules[0].getCVTerms()[3].getUseTerm(), false)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[3].getTermName(), "Sample Mass")
  TEST_EQUAL(rules[0].getCVTerms()[3].getIsRepeatable(), false)
  TEST_EQUAL(rules[0].getCVTerms()[3].getAllowChildren(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[3].getCVIdentifierRef(), "PSI")

  TEST_STRING_EQUAL(rules[0].getCVTerms()[4].getAccession(), "PSI:1000005")
  TEST_EQUAL(rules[0].getCVTerms()[4].getUseTermName(), true)
  TEST_EQUAL(rules[0].getCVTerms()[4].getUseTerm(), false)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[4].getTermName(), "Sample Volume")
  TEST_EQUAL(rules[0].getCVTerms()[4].getIsRepeatable(), false)
  TEST_EQUAL(rules[0].getCVTerms()[4].getAllowChildren(), false)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[4].getCVIdentifierRef(), "PSI")

  TEST_STRING_EQUAL(rules[0].getCVTerms()[5].getAccession(), "PSI:1000006")
  TEST_EQUAL(rules[0].getCVTerms()[5].getUseTermName(), false)
  TEST_EQUAL(rules[0].getCVTerms()[5].getUseTerm(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[5].getTermName(), "Sample Concentration")
  TEST_EQUAL(rules[0].getCVTerms()[5].getIsRepeatable(), true)
  TEST_EQUAL(rules[0].getCVTerms()[5].getAllowChildren(), true)
  TEST_STRING_EQUAL(rules[0].getCVTerms()[5].getCVIdentifierRef(), "PSI")


	TEST_EQUAL(mappings.getCVReferences().size(), 1)
	TEST_STRING_EQUAL(mappings.getCVReferences().begin()->getName(), "mzData CV")
	TEST_STRING_EQUAL(mappings.getCVReferences().begin()->getIdentifier(), "PSI")
	TEST_EQUAL(mappings.hasCVReference("PSI"), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



