// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/CVMappingFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CVMappingFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappingFile* ptr = 0;
START_SECTION(CVMappingFile())
{
	ptr = new CVMappingFile();
	TEST_NOT_EQUAL(ptr, 0)
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



