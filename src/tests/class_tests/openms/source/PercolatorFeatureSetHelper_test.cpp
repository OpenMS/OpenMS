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

///////////////////////////
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

bool check_pepids(const vector<PeptideIdentification>& check, const vector<PeptideIdentification>& against)
{
    std::vector<String> upk, upkc;
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

bool check_proids(const vector<ProteinIdentification>& check, const vector<ProteinIdentification>& against, const vector<String>& fs)
{
    TEST_EQUAL(check.size(), against.size())
    if (check.size()!= against.size())
        return false;
    for (size_t i = 0; i < check.size(); ++i)
      TEST_EQUAL(check[i].getHits().size(), against[i].getHits().size())

    String efc = check.front().getSearchParameters().getMetaValue("extra_features");
    TEST_STRING_EQUAL(efc, ListUtils::concatenate(fs, ","))
    return true;
}

START_TEST(PercolatorFeatureSetHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

STATUS("Preparing test inputs.")

std::vector< PeptideIdentification > comet_check_pids;
std::vector< PeptideIdentification > msgf_check_pids;
std::vector< PeptideIdentification > xtandem_check_pids;
std::vector< PeptideIdentification > merge_check_pids;
std::vector< PeptideIdentification > concat_check_pids;
std::vector< ProteinIdentification > comet_check_pods;
std::vector< ProteinIdentification > msgf_check_pods;
std::vector< ProteinIdentification > xtandem_check_pods;
std::vector< ProteinIdentification > concat_check_pods;
std::vector< ProteinIdentification > merge_check_pods;

IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc_check.idXML"), comet_check_pods, comet_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc_check.idXML"), msgf_check_pods, msgf_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("xtandem.topperc_check.idXML"), xtandem_check_pods, xtandem_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("combined.merge.perco.in.idXML"), merge_check_pods, merge_check_pids);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("combined.concat.perco.in.idXML"), concat_check_pods, concat_check_pids);

START_SECTION((static void concatMULTISEPeptideIds(std::vector< PeptideIdentification > &all_peptide_ids, std::vector< PeptideIdentification > &new_peptide_ids, String search_engine)))
{
    StringList fs;
    std::vector< PeptideIdentification > comet_pids;
    std::vector< ProteinIdentification > comet_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc.idXML"), comet_pods, comet_pids);

    std::vector< PeptideIdentification > msgf_pids;
    std::vector< ProteinIdentification > msgf_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc.idXML"), msgf_pods, msgf_pids);

    StringList ses = ListUtils::create<String>("MS-GF+,Comet");
    std::vector< PeptideIdentification > concat_pids;
    PercolatorFeatureSetHelper::concatMULTISEPeptideIds(concat_pids, msgf_pids, "MS-GF+");
    PercolatorFeatureSetHelper::concatMULTISEPeptideIds(concat_pids, comet_pids, "Comet");
    PercolatorFeatureSetHelper::addCONCATSEFeatures(concat_pids, ses, fs);

    //check completeness of feature construction
    ABORT_IF(!check_pepids(concat_check_pids, concat_pids));
}
END_SECTION

START_SECTION((static void mergeMULTISEPeptideIds(std::vector< PeptideIdentification > &all_peptide_ids, std::vector< PeptideIdentification > &new_peptide_ids, String search_engine)))
{
    std::vector< PeptideIdentification > comet_pids;
    std::vector< ProteinIdentification > comet_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc.idXML"), comet_pods, comet_pids);

    std::vector< PeptideIdentification > msgf_pids;
    std::vector< ProteinIdentification > msgf_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc.idXML"), msgf_pods, msgf_pids);

    std::vector< PeptideIdentification > merge_pids;
    StringList ses = ListUtils::create<String>("MS-GF+,Comet");
    PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(merge_pids, msgf_pids, "MS-GF+");
    PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(merge_pids, comet_pids, "Comet");
    StringList empty_extra;
    PercolatorFeatureSetHelper::addMULTISEFeatures(merge_pids, ses, empty_extra, true);
    TEST_EQUAL(merge_pids.size(),4)
    for (size_t i = merge_pids.size()-1; i > 0; --i)
    {
      PercolatorFeatureSetHelper::checkExtraFeatures(merge_pids[i].getHits(), empty_extra);  // also check against empty extra features list and inconsistency removal
      merge_pids.erase(merge_pids.begin()+i);  //erase to be able to use completeness check function below
    }
    TEST_EQUAL(merge_pids.size(),1)
    //check completeness of feature construction
    ABORT_IF(!check_pepids(merge_check_pids, merge_pids));
}
END_SECTION

START_SECTION((static void mergeMULTISEProteinIds(std::vector< ProteinIdentification > &all_protein_ids, std::vector< ProteinIdentification > &new_protein_ids)))
{
    StringList fs;
    std::vector< PeptideIdentification > comet_pids;
    std::vector< ProteinIdentification > comet_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("comet.topperc.idXML"), comet_pods, comet_pids);

    std::vector< PeptideIdentification > msgf_pids;
    std::vector< ProteinIdentification > msgf_pods;
    IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("msgf.topperc.idXML"), msgf_pods, msgf_pids);

    std::vector< ProteinIdentification > merge_pods;
    PercolatorFeatureSetHelper::mergeMULTISEProteinIds(merge_pods, msgf_pods);
    PercolatorFeatureSetHelper::mergeMULTISEProteinIds(merge_pods, comet_pods);

    std::vector< PeptideIdentification > merge_pids;
    StringList ses = ListUtils::create<String>("MS-GF+,Comet");
    PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(merge_pids, msgf_pids, "MS-GF+");
    PercolatorFeatureSetHelper::mergeMULTISEPeptideIds(merge_pids, comet_pids, "Comet");
    PercolatorFeatureSetHelper::addMULTISEFeatures(merge_pids, ses, fs, true);

    //check completeness of feature construction
    ABORT_IF(!check_proids(merge_check_pods, merge_pods, fs));
}
END_SECTION

START_SECTION((static void addMSGFFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &feature_set)))
{
    StringList fs;
    std::vector< PeptideIdentification > msgf_pids;
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
    std::vector< PeptideIdentification > xtandem_pids;
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
    std::vector< PeptideIdentification > comet_pids;
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

START_SECTION((static void addMULTISEFeatures(std::vector< PeptideIdentification > &peptide_ids, StringList &search_engines_used, StringList &feature_set, bool complete_only=true, bool limits_imputation=false)))
{
  NOT_TESTABLE  // actually tested in combination with mergeMULTISEPeptideIds
}
END_SECTION

START_SECTION((static void addCONCATSEFeatures(std::vector< PeptideIdentification > &peptide_id_list, StringList &search_engines_used, StringList &feature_set)))
{
  NOT_TESTABLE  // actually tested in combination with concatMULTISEPeptideIds
}
END_SECTION

START_SECTION((static void checkExtraFeatures(const std::vector< PeptideHit > &psms, StringList &extra_features)))
{
  NOT_TESTABLE  // actually tested in combination with mergeMULTISEPeptideIds
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
