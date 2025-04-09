// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;

START_TEST(IDMergerAlgorithm, "$Id$")

    START_SECTION(insertRun())
    {

      ProteinIdentification pr1;
      pr1.setIdentifier("PR1");
      pr1.setPrimaryMSRunPath({"f1r1.mzML","f2r1.mzML"});
      pr1.getHits().emplace_back(0,1,"A","");
      pr1.getHits().emplace_back(0,1,"B","");
      pr1.getHits().emplace_back(0,1,"C","");
      ProteinIdentification pr2;
      pr2.setIdentifier("PR2");
      pr2.setPrimaryMSRunPath({"f1r2.mzML","f2r2.mzML"});
      pr2.getHits().emplace_back(0,1,"A","");
      pr2.getHits().emplace_back(0,1,"D","");
      pr2.getHits().emplace_back(0,1,"E","");
      ProteinIdentification pr3;
      pr3.setIdentifier("PR3");
      pr3.setPrimaryMSRunPath({"f1r3.mzML","f2r3.mzML"});
      pr3.getHits().emplace_back(0,1,"D","");
      pr3.getHits().emplace_back(0,1,"E","");
      pr3.getHits().emplace_back(0,1,"F","");
      pr3.getHits().emplace_back(0,1,"G","");
      ProteinIdentification pr4; //empty
      pr4.setIdentifier("PR4");
      pr4.setPrimaryMSRunPath({"control.mzML"});

      PeptideHit ph0{0,1,1,AASequence::fromString("AA")};
      ph0.addPeptideEvidence(PeptideEvidence{"A",1,3,'A','A'});
      ph0.addPeptideEvidence(PeptideEvidence{"B",1,3,'A','A'});
      PeptideHit ph1{0,1,1,AASequence::fromString("AAA")};
      ph1.addPeptideEvidence(PeptideEvidence{"A",1,4,'A','A'});
      PeptideHit ph11{0,1,1,AASequence::fromString("AAC")};
      ph11.addPeptideEvidence(PeptideEvidence{"C",1,4,'A','A'});
      PeptideHit ph2{0,1,1,AASequence::fromString("AAAA")};
      ph2.addPeptideEvidence(PeptideEvidence{"D",1,5,'A','A'});
      ph2.addPeptideEvidence(PeptideEvidence{"E",1,5,'A','A'});
      PeptideHit ph3{0,1,1,AASequence::fromString("AAAAA")};
      ph3.addPeptideEvidence(PeptideEvidence{"D",1,6,'A','A'});
      ph3.addPeptideEvidence(PeptideEvidence{"E",1,6,'A','A'});
      PeptideHit ph4{0,1,1,AASequence::fromString("AAAAAA")};
      ph4.addPeptideEvidence(PeptideEvidence{"F",1,7,'A','A'});
      // ph5 same pep sequence but different proteins -> this actually means that there was an error or a different
      // protein database was used or something was filtered. But we cannot recover it here during merging.
      // you need to re-index. We can think about warning but this requires additional checks & datastructures.
      PeptideHit ph5{0.1,1,1,AASequence::fromString("AAA")};
      ph5.addPeptideEvidence(PeptideEvidence{"F",1,4,'A','A'});

      PeptideIdentification pe1;
      pe1.setIdentifier("PR1");
      pe1.getHits().push_back(ph0);
      pe1.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0);

      PeptideIdentification pe2;
      pe2.setIdentifier("PR1");
      pe2.getHits().push_back(ph1);
      pe2.getHits().push_back(ph11);
      pe2.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0);

      PeptideIdentification pe3;
      pe3.setIdentifier("PR2");
      pe3.getHits().push_back(ph0); // how to handle accessions that are not in the corresponding list of proteins?
      //currently ignore and add nonetheless.
      pe3.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0);

      PeptideIdentification pe4;
      pe4.setIdentifier("PR2");
      pe4.getHits().push_back(ph2);
      pe4.getHits().push_back(ph3);
      pe4.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,1);

      PeptideIdentification pe5;
      pe5.setIdentifier("PR3");
      pe5.getHits().push_back(ph2);
      pe5.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0);

      PeptideIdentification pe6;
      pe6.setIdentifier("PR3");
      pe6.getHits().push_back(ph3);
      pe6.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,1);

      PeptideIdentification pe7;
      pe7.setIdentifier("PR3");
      pe7.getHits().push_back(ph4);
      pe7.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0);

      PeptideIdentification pe8;
      pe8.setIdentifier("PR3");
      pe8.getHits().push_back(ph5);
      pe8.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,0); //can happen if second file had no IDs

      PeptideIdentification pe9;
      pe9.setIdentifier("PR5"); // non-existent run: this will be ignored
      pe9.getHits().push_back(ph5);
      pe9.setMetaValue(Constants::UserParam::ID_MERGE_INDEX,564);

      vector<PeptideIdentification> pes{pe1,pe2,pe3,pe4,pe5,pe6,pe7,pe8,pe9};
      vector<ProteinIdentification> prs{pr1,pr2,pr3,pr4};
      IDMergerAlgorithm ima("mymerge");
      ima.insertRuns(prs, pes);
      ProteinIdentification prres;
      vector<PeptideIdentification> peres;
      ima.returnResultsAndClear(prres,peres);

      TEST_EQUAL(pes.size(), 9)
      TEST_EQUAL(peres.size(), 8)
      TEST_EQUAL(prres.getHits().size(), 7)
      StringList toFill; prres.getPrimaryMSRunPath(toFill);
      TEST_EQUAL(toFill.size(), 7)
      TEST_EQUAL(static_cast<int>(peres[2].getMetaValue(Constants::UserParam::ID_MERGE_INDEX)), 2)

    }
    END_SECTION

    START_SECTION(insertRun())
    {
      IdXMLFile f;
      vector<ProteinIdentification> pr1;
      vector<PeptideIdentification> pe1;
      f.load(OPENMS_GET_TEST_DATA_PATH("newIDMergerTest1.idXML"),pr1,pe1);
      Size pe1size = pe1.size();

      vector<ProteinIdentification> pr2;
      vector<PeptideIdentification> pe2;
      f.load(OPENMS_GET_TEST_DATA_PATH("newIDMergerTest2.idXML"),pr2,pe2);
      Size pe2size = pe2.size();

      IDMergerAlgorithm ima("mymerge");
      ima.insertRuns(std::move(pr1), std::move(pe1));
      ima.insertRuns(std::move(pr2), std::move(pe2));

      TEST_EXCEPTION(Exception::MissingInformation, ima.insertRuns({ProteinIdentification{}}, {PeptideIdentification{}}))

      ProteinIdentification prres;
      vector<PeptideIdentification> peres;
      ima.returnResultsAndClear(prres,peres);

      TEST_EQUAL(prres.getHits().size(),6)
      TEST_EQUAL(peres.size(),pe1size + pe2size)
      TEST_EQUAL(pr1.size(),0)
      TEST_EQUAL(pr2.size(),0)
      TEST_EQUAL(pe1.size(),0)
      TEST_EQUAL(pe2.size(),0)
      TEST_EQUAL(prres.getIdentifier().hasPrefix("mymerge"), true)
      StringList toFill; prres.getPrimaryMSRunPath(toFill);
      TEST_EQUAL(toFill.size(), 2)

    }
    END_SECTION

    START_SECTION(check search setting consistency)
        {
          IdXMLFile f;
          vector<ProteinIdentification> pr1;
          vector<PeptideIdentification> pe1;
          f.load(OPENMS_GET_TEST_DATA_PATH("newIDMergerTest1.idXML"),pr1,pe1);

          vector<ProteinIdentification> pr2;
          vector<PeptideIdentification> pe2;
          f.load(OPENMS_GET_TEST_DATA_PATH("newIDMergerTest2.idXML"),pr2,pe2);
          // fail with different db filename
          pr2[0].getSearchParameters().db = "baz";

          IDMergerAlgorithm ima("mymerge");
          ima.insertRuns(std::move(pr1), std::move(pe1));
          TEST_EXCEPTION(Exception::MissingInformation,ima.insertRuns(pr2, pe2))

          // check windows path with correct filename
          String fn = "C:\\foo\\s_pyo_sf370_potato_human_target_decoy_with_contaminants.fasta";
          pr2[0].getSearchParameters().db = fn;

          ima.insertRuns(std::move(pr2), std::move(pe2));

          ProteinIdentification prres;
          vector<PeptideIdentification> peres;
          ima.returnResultsAndClear(prres,peres);

          TEST_EQUAL(prres.getHits().size(),6)
        }
    END_SECTION

END_TEST
