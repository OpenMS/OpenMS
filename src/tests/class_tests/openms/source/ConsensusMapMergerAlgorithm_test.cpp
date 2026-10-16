// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMapMergerAlgorithm, "$Id$")

    START_SECTION(mergeAllIDRuns)
      {
        ConsensusXMLFile cf;
        ConsensusMap cmap;
        cf.load(OPENMS_GET_TEST_DATA_PATH("BSA.consensusXML"), cmap);
        ConsensusMapMergerAlgorithm cmerge;
        cmerge.mergeAllIDRuns(cmap);
        TEST_EQUAL(cmap.getProteinIdentifications().size(), 1)
      }
    END_SECTION

    START_SECTION(mergeProteinsAcrossFractionsAndReplicates (no Design))
      {
        ConsensusXMLFile cf;
        ConsensusMap cmap;
        cf.load(OPENMS_GET_TEST_DATA_PATH("BSA.consensusXML"), cmap);
        ConsensusMapMergerAlgorithm cmerge;
        ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(cmap);
        cmerge.mergeProteinsAcrossFractionsAndReplicates(cmap, ed);
        //without a special experimental design on sample level, runs are treated like replicates
        // or fractions and all are merged
        TEST_EQUAL(cmap.getProteinIdentifications().size(), 1)
        StringList toFill; cmap.getProteinIdentifications()[0].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 6)
      }
    END_SECTION

    START_SECTION(mergeProteinsAcrossFractionsAndReplicates)
      {
        ConsensusXMLFile cf;
        ConsensusMap cmap;
        cf.load(OPENMS_GET_TEST_DATA_PATH("BSA.consensusXML"), cmap);
        ConsensusMapMergerAlgorithm cmerge;
        ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(cmap);
        ExperimentalDesign::SampleSection ss{
            {{"1","C1"},{"2","C2"},{"3","C3"}},
            {{"1",0},{"2",1},{"3",2}},
            {{"Sample",0},{"Condition",1}}
        };
        ed.setSampleSection(ss);
        // The file each PSM came from, by its original run (one file per run in this map)
        std::map<std::string, std::string> file_of_run;
        for (const auto& run : cmap.getProteinIdentifications())
        {
          StringList files;
          run.getPrimaryMSRunPath(files);
          file_of_run[run.getIdentifier()] = files.at(0);
        }
        std::vector<std::string> expected_files;
        cmap.applyFunctionOnPeptideIDs([&](PeptideIdentification& pid) { expected_files.push_back(file_of_run.at(pid.getIdentifier())); });
        cmerge.mergeProteinsAcrossFractionsAndReplicates(cmap, ed);
        TEST_EQUAL(cmap.getProteinIdentifications().size(), 3)
        StringList toFill; cmap.getProteinIdentifications()[0].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA1_F2.mzML")
        toFill.clear(); cmap.getProteinIdentifications()[1].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA2_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA2_F2.mzML")
        toFill.clear(); cmap.getProteinIdentifications()[2].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA3_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA3_F2.mzML")
        TEST_EQUAL(cmap.getProteinIdentifications().size(), 3)
        // id_merge_index points into the merged run's file list at the file each PSM came from
        std::map<std::string, StringList> files_of_run;
        for (const auto& run : cmap.getProteinIdentifications())
          run.getPrimaryMSRunPath(files_of_run[run.getIdentifier()]);
        std::vector<std::string> merged_files;
        cmap.applyFunctionOnPeptideIDs([&](PeptideIdentification& pid) {
          merged_files.push_back(files_of_run.at(pid.getIdentifier()).at(static_cast<Size>(static_cast<Int>(pid.getMetaValue("id_merge_index")))));
        });
        TEST_EQUAL(merged_files.size(), expected_files.size())
        TEST_EQUAL(merged_files == expected_files, true)
      }
    END_SECTION

END_TEST
