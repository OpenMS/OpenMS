// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
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
        TEST_EQUAL(IdentificationDataConverter::proteinIdentifications(cmap).size(), 1)
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
        const auto merged = IdentificationDataConverter::proteinIdentifications(cmap);
        TEST_EQUAL(merged.size(), 1)
        StringList toFill; merged[0].getPrimaryMSRunPath(toFill);
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
        const auto before = IdentificationDataConverter::exportMapIdentifications(cmap);
        for (const auto& run : before.proteins)
        {
          StringList files;
          run.getPrimaryMSRunPath(files);
          file_of_run[run.getIdentifier()] = files.at(0);
        }
        std::vector<std::string> expected_files;
        for (const auto* pid : before.peptides()) expected_files.push_back(file_of_run.at(pid->getIdentifier()));
        cmerge.mergeProteinsAcrossFractionsAndReplicates(cmap, ed);
        const auto after = IdentificationDataConverter::exportMapIdentifications(cmap);
        TEST_EQUAL(after.proteins.size(), 3)
        StringList toFill; after.proteins[0].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA1_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA1_F2.mzML")
        toFill.clear(); after.proteins[1].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA2_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA2_F2.mzML")
        toFill.clear(); after.proteins[2].getPrimaryMSRunPath(toFill);
        TEST_EQUAL(toFill.size(), 2)
        TEST_EQUAL(toFill[0], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA3_F1.mzML")
        TEST_EQUAL(toFill[1], "/Users/pfeuffer/git/OpenMS-inference-src/share/OpenMS/examples/FRACTIONS/BSA3_F2.mzML")
        // id_merge_index points into the merged run's file list at the file each PSM came from
        std::map<std::string, StringList> files_of_run;
        for (const auto& run : after.proteins)
          run.getPrimaryMSRunPath(files_of_run[run.getIdentifier()]);
        std::vector<std::string> merged_files;
        for (const auto* pid : after.peptides())
          merged_files.push_back(files_of_run.at(pid->getIdentifier()).at(static_cast<Size>(static_cast<Int>(pid->getMetaValue("id_merge_index")))));
        TEST_EQUAL(merged_files.size(), expected_files.size())
        TEST_EQUAL(merged_files == expected_files, true)
      }
    END_SECTION

    START_SECTION([EXTRA] mergeAllIDRuns with identification data)
      {
        ConsensusXMLFile cf;
        ConsensusMap cmap;
        cf.load(OPENMS_GET_TEST_DATA_PATH("BSA.consensusXML"), cmap);
        const auto data = cmap.getIdentificationData();
        ConsensusMapMergerAlgorithm cmerge;
        cmerge.mergeAllIDRuns(cmap);
        // the runs stay; an inference result pools them
        const auto& merged = cmap.getIdentificationData();
        TEST_EQUAL(merged.getRuns().size(), data.getRuns().size())
        TEST_EQUAL(merged.getRuns()[0] == data.getRuns()[0], true)
        TEST_EQUAL(merged.getInferenceResults().size(), 1)
        ABORT_IF(merged.getInferenceResults().size() != 1)
        TEST_EQUAL(merged.getInferenceResults()[0].identifier, "merged")
        TEST_EQUAL(merged.getInferenceResults()[0].inputs.size(), data.getRuns().size())
        // export writes the merged protein run
        const auto exported = IdentificationDataConverter::proteinIdentifications(cmap);
        TEST_EQUAL(exported.size(), 1)
        ABORT_IF(exported.size() != 1)
        TEST_EQUAL(exported[0].getIdentifier(), "merged")
        StringList files;
        exported[0].getPrimaryMSRunPath(files);
        TEST_EQUAL(files.size(), 6)
      }
    END_SECTION

    START_SECTION([EXTRA] mergeProteinsAcrossFractionsAndReplicates with identification data)
      {
        ConsensusXMLFile cf;
        ConsensusMap cmap;
        cf.load(OPENMS_GET_TEST_DATA_PATH("BSA.consensusXML"), cmap);
        ExperimentalDesign ed = ExperimentalDesign::fromConsensusMap(cmap);
        ExperimentalDesign::SampleSection ss{
            {{"1","C1"},{"2","C2"},{"3","C3"}},
            {{"1",0},{"2",1},{"3",2}},
            {{"Sample",0},{"Condition",1}}
        };
        ed.setSampleSection(ss);
        ConsensusMapMergerAlgorithm cmerge;
        cmerge.mergeProteinsAcrossFractionsAndReplicates(cmap, ed);
        // an inference result pools the runs of each condition
        const auto& results = cmap.getIdentificationData().getInferenceResults();
        TEST_EQUAL(results.size(), 3)
        ABORT_IF(results.size() != 3)
        TEST_EQUAL(results[0].identifier, "condition0")
        TEST_EQUAL(results[0].inputs.size(), 2)
        // and export writes the merged protein runs, one per condition with the files of its runs
        const auto exported = IdentificationDataConverter::proteinIdentifications(cmap);
        TEST_EQUAL(exported.size(), 3)
        ABORT_IF(exported.size() != 3)
        for (Size i = 0; i < 3; ++i)
        {
          TEST_EQUAL(exported[i].getIdentifier(), results[i].identifier)
          StringList files;
          exported[i].getPrimaryMSRunPath(files);
          TEST_EQUAL(files.size(), 2)
        }
      }
    END_SECTION

END_TEST
