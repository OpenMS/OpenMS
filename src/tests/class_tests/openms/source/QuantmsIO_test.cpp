// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/QuantmsIO.h>
///////////////////////////

#ifdef WITH_PARQUET

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace OpenMS;
using namespace std;

START_TEST(QuantmsIO, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QuantmsIO* ptr = nullptr;
QuantmsIO* null_ptr = nullptr;

START_SECTION(QuantmsIO())
{
  ptr = new QuantmsIO();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuantmsIO())
{
  delete ptr;
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<ProteinIdentification>& protein_identifications, const PeptideIdentificationList& peptide_identifications)))
{
  QuantmsIO file;
  
  // Create test data
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  ProteinIdentification protein_id;
  protein_id.setIdentifier("test_search");
  protein_id.setSearchEngine("TestEngine");
  protein_id.setScoreType("TestScore", true);
  protein_ids.push_back(protein_id);
  
  PeptideIdentification peptide_id;
  peptide_id.setIdentifier("test_search");
  peptide_id.setRT(1234.5);
  peptide_id.setMZ(500.25);
  peptide_id.setScoreType("TestScore");
  
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  hit.setScore(0.95);
  hit.setCharge(2);
  
  PeptideEvidence evidence;
  evidence.setProteinAccession("TEST_PROTEIN");
  hit.setPeptideEvidences(vector<PeptideEvidence>{evidence});
  
  peptide_id.setHits(vector<PeptideHit>{hit});
  peptide_ids.push_back(peptide_id);
  
  String output_file = OPENMS_GET_TEST_DATA_PATH("") + "QuantmsIO_output_test.parquet";
  
  TEST_EXCEPTION_WITH_MESSAGE(Exception::NotImplemented, file.store(output_file, protein_ids, peptide_ids), ""); // Should throw when parquet not available
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

#else // WITH_PARQUET

START_TEST(QuantmsIO, "$Id$")

QuantmsIO* ptr = nullptr;
QuantmsIO* null_ptr = nullptr;

START_SECTION(QuantmsIO())
{
  ptr = new QuantmsIO();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuantmsIO())
{
  delete ptr;
}
END_SECTION

START_SECTION((void store(const String& filename, const std::vector<ProteinIdentification>& protein_identifications, const PeptideIdentificationList& peptide_identifications)))
{
  QuantmsIO file;
  
  vector<ProteinIdentification> protein_ids;
  PeptideIdentificationList peptide_ids;
  
  String output_file = OPENMS_GET_TEST_DATA_PATH("") + "QuantmsIO_output_test.parquet";
  
  TEST_EXCEPTION(Exception::NotImplemented, file.store(output_file, protein_ids, peptide_ids));
}
END_SECTION

END_TEST

#endif // WITH_PARQUET