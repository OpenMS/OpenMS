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
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricNormalizer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricNormalizer* ptr = nullptr;
IsobaricNormalizer* null_ptr = nullptr;

// 
ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricNormalizer(const IsobaricQuantitationMethod *const quant_method)))
{
	ptr = new IsobaricNormalizer(&quant_meth);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricNormalizer())
{
	delete ptr;
}
END_SECTION

START_SECTION((IsobaricNormalizer(const IsobaricNormalizer &other)))
{
  IsobaricNormalizer normalizer(&quant_meth);
  IsobaricNormalizer * normalizer2 = new IsobaricNormalizer(normalizer);
  
	TEST_NOT_EQUAL(normalizer2, null_ptr)  
  delete normalizer2;
	
  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((IsobaricNormalizer& operator=(const IsobaricNormalizer &rhs)))
{
  IsobaricNormalizer normalizer(&quant_meth);
  IsobaricNormalizer normalizer2(&quant_meth);
  
  normalizer2 = normalizer;
  
  // equality cannot be checked
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void normalize(ConsensusMap &consensus_map)))
{
  IsobaricNormalizer normalizer(&quant_meth);
    
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer.consensusXML"),cm_in);

  normalizer.normalize(cm_in);

	String cm_file_out;
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_in);

	WHITELIST("<?xml-stylesheet");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricNormalizer_out.consensusXML"));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
