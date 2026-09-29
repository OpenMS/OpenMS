// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: MATHIAS WALZER$
// $Authors: MATHIAS WALZER$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <set>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

bool check_pepids(const PeptideIdentificationList& check, const PeptideIdentificationList& against)
{
    std::vector<std::string> upk, upkc;
    TEST_EQUAL(check.size(), against.size())
    if (check.size() != against.size())
        return false;
    for (size_t i = 0; i < check.size(); ++i)
    {
      TEST_EQUAL(check[i].getHits().size(), against[i].getHits().size())
      for (size_t j = 0; j < check[i].getHits().size(); ++j)
      {
        check [i].getHits()[j].getKeys(upkc);
        against[i].getHits()[j].getKeys(upk);
        TEST_EQUAL(upkc.size(), upk.size())
        if (upkc.size() != upk.size())
            return false;
        for (size_t k = 0; k < upk.size(); ++k)
           TEST_STRING_EQUAL(upkc[k],upk[k])
      }
    }
    return true;
}

bool check_proids(const vector<ProteinIdentification>& check, const vector<ProteinIdentification>& against, const vector<std::string>& fs)
{
    TEST_EQUAL(check.size(), against.size())
    if (check.size()!= against.size())
        return false;
    for (size_t i = 0; i < check.size(); ++i)
      TEST_EQUAL(check[i].getHits().size(), against[i].getHits().size())

    std::string efc = check.front().getSearchParameters().getMetaValue("extra_features");
    TEST_STRING_EQUAL(efc, ListUtils::concatenate(fs, ","))
    return true;
}

START_TEST(PercolatorFeatureSetHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

STATUS("Preparing test inputs.")

PeptideIdentificationList comet_check_pids;
PeptideIdentificationList msgf_check_pids;
PeptideIdentificationList xtandem_check_pids;
std::vector< ProteinIdentification > comet_check_pods;
std::vector< ProteinIdentification > msgf_check_pods;
std::vector< ProteinIdentification > xtandem_check_pods;

IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc_check.idXML"), comet_check_pods, comet_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc_check.idXML"), msgf_check_pods, msgf_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("xtandem.topperc_check.idXML"), xtandem_check_pods, xtandem_check_pids);

START_SECTION((static void mergeMULTISEProteinIds(std::vector< ProteinIdentification > &all_protein_ids, std::vector< ProteinIdentification > &new_protein_ids)))
{
    PeptideIdentificationList comet_pids;
    std::vector< ProteinIdentification > comet_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc.idXML"), comet_pods, comet_pids);

    PeptideIdentificationList msgf_pids;
    std::vector< ProteinIdentification > msgf_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc.idXML"), msgf_pods, msgf_pids);

    // expected result: the union of both runs' ProteinHits, keyed by accession
    std::set<std::string> expected_accessions;
    for (const ProteinHit& h : msgf_pods.front().getHits()) expected_accessions.insert(h.getAccession());
    for (const ProteinHit& h : comet_pods.front().getHits()) expected_accessions.insert(h.getAccession());
    const std::string expected_db = msgf_pods.front().getSearchParameters().db;

    std::vector< ProteinIdentification > merge_pods;
    PercolatorFeatureSetHelper::mergeMULTISEProteinIds(merge_pods, msgf_pods);

    // merging into an empty target creates a single run that adopts the incoming run's settings
    TEST_EQUAL(merge_pods.size(), 1)
    TEST_STRING_EQUAL(merge_pods.front().getSearchEngine(), "MS-GF+")

    PercolatorFeatureSetHelper::mergeMULTISEProteinIds(merge_pods, comet_pods);

    // a second run from a different engine collapses the search engine label, never the run
    TEST_EQUAL(merge_pods.size(), 1)
    TEST_STRING_EQUAL(merge_pods.front().getSearchEngine(), "multiple")
    // search parameters stay those of the first merged run
    TEST_STRING_EQUAL(merge_pods.front().getSearchParameters().db, expected_db)

    // ProteinHits are unioned by accession, so shared accessions must not be duplicated
    TEST_EQUAL(merge_pods.front().getHits().size(), expected_accessions.size())
    std::set<std::string> merged_accessions;
    for (const ProteinHit& h : merge_pods.front().getHits()) merged_accessions.insert(h.getAccession());
    TEST_EQUAL(merged_accessions.size(), expected_accessions.size())
    TEST_TRUE(merged_accessions == expected_accessions)
}
END_SECTION

START_SECTION((static void addMSGFFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
    StringList fs;
    PeptideIdentificationList msgf_pids;
    std::vector< ProteinIdentification > msgf_pods;

    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc.idXML"), msgf_pods, msgf_pids);
    PercolatorFeatureSetHelper::addMSGFFeatures(msgf_pids,fs);

    //check completeness of feature construction
    ABORT_IF(!check_pepids(msgf_check_pids, msgf_pids));

    //check registration of percolator features for adapter
    ABORT_IF(!check_proids(msgf_check_pods, msgf_pods, fs));
}
END_SECTION

START_SECTION((static void addXTANDEMFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
    StringList fs;
    PeptideIdentificationList xtandem_pids;
    std::vector< ProteinIdentification > xtandem_pods;

    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("xtandem.topperc.idXML"), xtandem_pods, xtandem_pids);
    PercolatorFeatureSetHelper::addXTANDEMFeatures(xtandem_pids, fs);

    //check completeness of feature construction
    ABORT_IF(!check_pepids(xtandem_check_pids, xtandem_pids));

    //check registration of percolator features for adapter
    ABORT_IF(!check_proids(xtandem_check_pods, xtandem_pods, fs));
}
END_SECTION

START_SECTION((static void addCOMETFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
    StringList fs;
    PeptideIdentificationList comet_pids;
    std::vector< ProteinIdentification > comet_pods;

    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc.idXML"), comet_pods, comet_pids);
    PercolatorFeatureSetHelper::addCOMETFeatures(comet_pids, fs);

    //check completeness of feature construction
    ABORT_IF(!check_pepids(comet_check_pids, comet_pids));

    //check registration of percolator features for adapter
    ABORT_IF(!check_proids(comet_check_pods, comet_pods, fs));
}
END_SECTION

START_SECTION((static void addMASCOTFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
  NOT_TESTABLE  // yet
}
END_SECTION

START_SECTION((static void addANDESFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
    // andes annotates hits with numeric "andes:"-prefixed MetaValues; these (and only these,
    // when numeric) should be collected as Percolator features. Built in-memory to avoid a
    // dependency on an andes-produced test-data file.
    PeptideHit hit;
    hit.setMetaValue("andes:RankScore", 12.5);          // double -> collected
    hit.setMetaValue("andes:NumMatchedMainIons", 7);    // int -> collected
    hit.setMetaValue("andes:DeltaRankScore", 3.25);     // double -> collected
    hit.setMetaValue("andes:flag", "yes");              // string -> NOT collected
    hit.setMetaValue("target_decoy", "target");         // non-andes -> NOT collected

    PeptideIdentification pid;
    pid.insertHit(hit);
    PeptideIdentificationList andes_pids;
    andes_pids.push_back(pid);

    StringList fs;
    PercolatorFeatureSetHelper::addANDESFeatures(andes_pids, fs);

    // only the three numeric andes-prefixed features are registered (sorted/unique via std::set)
    TEST_EQUAL(fs.size(), 3)
    TEST_EQUAL(ListUtils::contains(fs, std::string("andes:RankScore")), true)
    TEST_EQUAL(ListUtils::contains(fs, std::string("andes:NumMatchedMainIons")), true)
    TEST_EQUAL(ListUtils::contains(fs, std::string("andes:DeltaRankScore")), true)
    TEST_EQUAL(ListUtils::contains(fs, std::string("andes:flag")), false)
    TEST_EQUAL(ListUtils::contains(fs, std::string("target_decoy")), false)
}
END_SECTION

START_SECTION((static void checkExtraFeatures(const std::vector< PeptideHit > &psms, StringList &extra_features)))
{
    // a requested extra feature survives only if it is present on every PSM
    PeptideHit hit_a;
    hit_a.setMetaValue("on_all", 1.0);
    hit_a.setMetaValue("on_first_only", 2.0);

    PeptideHit hit_b;
    hit_b.setMetaValue("on_all", 3.0);

    std::vector<PeptideHit> psms {hit_a, hit_b};

    StringList extra = ListUtils::create<std::string>("on_all,on_first_only,on_none");
    PercolatorFeatureSetHelper::checkExtraFeatures(psms, extra);

    TEST_EQUAL(extra.size(), 1)
    TEST_TRUE(ListUtils::contains(extra, std::string("on_all")))
    TEST_FALSE(ListUtils::contains(extra, std::string("on_first_only")))
    TEST_FALSE(ListUtils::contains(extra, std::string("on_none")))

    // an empty request stays empty and must not touch the PSMs
    StringList none;
    PercolatorFeatureSetHelper::checkExtraFeatures(psms, none);
    TEST_EQUAL(none.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
