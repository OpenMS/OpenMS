// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricQuantifier* ptr = nullptr;
IsobaricQuantifier* null_ptr = nullptr;

ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricQuantifier(const IsobaricQuantitationMethod *const quant_method)))
{
	ptr = new IsobaricQuantifier(&quant_meth);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricQuantifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((IsobaricQuantifier(const IsobaricQuantifier &other)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier * quantifier2 = new IsobaricQuantifier(quantifier);

	TEST_NOT_EQUAL(quantifier2, null_ptr)
  delete quantifier2;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IsobaricQuantifier& operator=(const IsobaricQuantifier &rhs)))
{
  IsobaricQuantifier quantifier(&quant_meth);
  IsobaricQuantifier quantifier2(&quant_meth);

  quantifier2 = quantifier;

  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void quantify(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier.consensusXML"),cm_in);

	IsobaricQuantifier iq(&quant_meth);
	Param p;
	p.setValue("normalization", "true");
  p.setValue("isotope_correction", "true");
	iq.setParameters(p);
	iq.quantify(cm_in,cm_out);

	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);

	WHITELIST("<?xml-stylesheet");
	// WHITELIST("<?xml-stylesheet,consensusElement id=");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricQuantifier_out.consensusXML"));
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
