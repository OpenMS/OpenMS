// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

PepNovoOutfile* ptr = nullptr;
PepNovoOutfile* nullPointer = nullptr;
START_SECTION(PepNovoOutfile())
  ptr = new PepNovoOutfile();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~PepNovoOutfile())
  delete ptr;
END_SECTION

START_SECTION((PepNovoOutfile& operator=(const PepNovoOutfile &pepnovo_outfile)))
  PepNovoOutfile pepnovo_outfile1;
  PepNovoOutfile pepnovo_outfile2;
  pepnovo_outfile2 = pepnovo_outfile1;
  PepNovoOutfile pepnovo_outfile3;
  pepnovo_outfile1 = PepNovoOutfile();
  TEST_EQUAL(( pepnovo_outfile2 == pepnovo_outfile3 ), true)
END_SECTION

START_SECTION((PepNovoOutfile(const PepNovoOutfile &pepnovo_outfile)))
  PepNovoOutfile pepnovo_outfile1;
  PepNovoOutfile pepnovo_outfile2(pepnovo_outfile1);
  PepNovoOutfile pepnovo_outfile3;
  pepnovo_outfile1 = PepNovoOutfile();
  TEST_EQUAL(( pepnovo_outfile2 == pepnovo_outfile3 ), true)
END_SECTION

START_SECTION((bool operator==(const PepNovoOutfile &pepnovo_outfile) const))
  PepNovoOutfile pepnovo_outfile1;
  PepNovoOutfile pepnovo_outfile2;
  TEST_EQUAL(( pepnovo_outfile1 == pepnovo_outfile2 ), true)
END_SECTION

PepNovoOutfile file;


START_SECTION((void load(const std::string &result_filename, std::vector< PeptideIdentification > &peptide_identifications, ProteinIdentification &protein_identification, const double &score_threshold, const IndexPosMappingType &id_rt_mz, const std::map< String, String > &mod_id_map)))
  std::vector< PeptideIdentification > peptide_identifications;
  ProteinIdentification protein_identification;
  map< String, double > filenames_and_precursor_retention_times;

  // test exceptions
  //TEST_EXCEPTION_WITH_MESSAGE(Exception::FileNotFound, file.load("a", peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), "the file 'a' could not be found")

  //TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out1"), peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), OPENMS_GET_TEST_DATA_PATH_MESSAGE("", "PepNovoOutfile.out1", " in: Not enough columns in file in line 2 (should be 8)!"))

  //TEST_EXCEPTION_WITH_MESSAGE(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out2"), peptide_identifications, protein_identification, 0.915f, filenames_and_precursor_retention_times), OPENMS_GET_TEST_DATA_PATH_MESSAGE("", "PepNovoOutfile.out2", " in: Not enough columns in file in line 7 (should be 8)!" ))

  peptide_identifications.clear();
  protein_identification.setHits(vector< ProteinHit >());


  // test the actual program
  map<String, String> key_to_mod;
  key_to_mod["K+42"]="Acetyl (K)";
  key_to_mod["Y+42"]="Acetyl (Y)";

  PepNovoOutfile::IndexPosMappingType rt_and_index, rt_and_index2;


  rt_and_index[0] = std::make_pair(1510.5732421875, 747.761901855469);
  rt_and_index[1] = std::make_pair(1530.11535644531, 549.856262207031);

  // check missing index-key ( rt_and_index[2] )
  TEST_EXCEPTION(Exception::ParseError, file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), peptide_identifications, protein_identification, -2.000f, rt_and_index, key_to_mod));
  rt_and_index[2] = std::make_pair(1533.16589355469, 358.174530029297);
  rt_and_index[3] = std::make_pair(1111, 2222);

  for (Size i=0; i<2; ++i)
  {
    if (i==0) rt_and_index2 = rt_and_index; // use explicit mapping
    else rt_and_index2 = PepNovoOutfile::IndexPosMappingType(); // try to reconstruct from title in pepnovo file


    file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), peptide_identifications, protein_identification, -2.000f, rt_and_index2, key_to_mod);

    TEST_EQUAL(peptide_identifications.size(), 4)
    ABORT_IF(peptide_identifications.size() != 4)

    TEST_EQUAL(peptide_identifications[0].getHits().size(), 5)
    TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), -2.0)
    TEST_REAL_SIMILAR(peptide_identifications[0].getMZ(), 747.761901855469)
    TEST_REAL_SIMILAR(peptide_identifications[0].getRT(), 1510.5732421875)

    TEST_EQUAL(peptide_identifications[1].getHits().size(), 14)
    TEST_REAL_SIMILAR(peptide_identifications[1].getSignificanceThreshold(), -2.0)
    TEST_REAL_SIMILAR(peptide_identifications[1].getMZ(), 549.856262207031)
    TEST_REAL_SIMILAR(peptide_identifications[1].getRT(), 1530.11535644531)

    TEST_EQUAL(peptide_identifications[2].getHits().size(), 20)
    TEST_EQUAL(peptide_identifications[3].getHits().size(), 20)

    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), -1.412)
    TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), AASequence::fromString("ADYGVTR"))
    TEST_EQUAL(peptide_identifications[0].getHits()[0].getRank(), 1)
    TEST_EQUAL(peptide_identifications[0].getHits()[0].getCharge(), 2)
    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getMetaValue("PnvScr"), 21.144)

    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), -1.483)
    TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), AASequence::fromString("SDYGVTR"))
    TEST_EQUAL(peptide_identifications[0].getHits()[1].getRank(), 2)
    TEST_EQUAL(peptide_identifications[0].getHits()[1].getCharge(), 2)
    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getMetaValue("PnvScr"), 18.239)




    file.load(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), peptide_identifications, protein_identification, -4.000f, rt_and_index2, key_to_mod);

    TEST_EQUAL(peptide_identifications.size(), 4)
    ABORT_IF( peptide_identifications.size() != 4 )
    TEST_EQUAL(peptide_identifications[3].getHits().size(), 20)
    TEST_REAL_SIMILAR(peptide_identifications[3].getSignificanceThreshold(), -4.0)
    TEST_REAL_SIMILAR(peptide_identifications[3].getHits()[11].getScore(),8.045)
    TEST_EQUAL(peptide_identifications[3].getHits()[11].getSequence(),  AASequence::fromString("GK(Acetyl)EAMAPK"))
    TEST_EQUAL(peptide_identifications[3].getHits()[11].getRank(), 12)
    TEST_EQUAL(peptide_identifications[3].getHits()[0].getCharge(), 2)

  }


END_SECTION


START_SECTION(void getSearchEngineAndVersion(const String& pepnovo_output_without_parameters_filename, ProteinIdentification& protein_identification))
  ProteinIdentification protein_identification;

  // test the actual program
  file.getSearchEngineAndVersion(OPENMS_GET_TEST_DATA_PATH("PepNovoOutfile.out"), protein_identification);
  TEST_EQUAL(protein_identification.getSearchEngine(), "PepNovo+");
  TEST_EQUAL(protein_identification.getSearchEngineVersion(), "Build 20081230");
  TEST_EQUAL(protein_identification.getSearchParameters().fragment_mass_tolerance, 0.5);
  TEST_REAL_SIMILAR(protein_identification.getSearchParameters().fragment_mass_tolerance, 0.5);
  TEST_REAL_SIMILAR(protein_identification.getSearchParameters().precursor_mass_tolerance, 2.5);
  TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications.size(), 2);
  if(protein_identification.getSearchParameters().variable_modifications.size()== 2)
  {
    TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications[0], "K+42");
    TEST_EQUAL(protein_identification.getSearchParameters().variable_modifications[1], "Y+42");
  }
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
