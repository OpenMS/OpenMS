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
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>
///////////////////////////

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifierStatistics.h>

using namespace OpenMS;
using namespace std;

ConsensusFeature getCFWithIntensites(double v[])
{
  ConsensusFeature cf;
  BaseFeature bf0, bf1, bf2, bf3;
  bf0.setIntensity(v[0]);
  bf1.setIntensity(v[1]);
  bf2.setIntensity(v[2]);
  bf3.setIntensity(v[3]);
  cf.insert(0, bf0);cf.insert(1, bf1);cf.insert(2, bf2);cf.insert(3, bf3);
  cf.setIntensity(v[0]+v[1]+v[2]+v[3]);
  return cf;
}

START_TEST(IsobaricIsotopeCorrector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//
ItraqFourPlexQuantitationMethod quant_meth;

START_SECTION((IsobaricQuantifierStatistics correctIsotopicImpurities(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
  { // check the run including output
    ConsensusXMLFile cm_file;
    ConsensusMap cm_in, cm_out;
    cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricIsotopeCorrector.consensusXML"),cm_in);

    // copy in/output
    cm_out = cm_in;

    //
    IsobaricQuantifierStatistics stats = IsobaricIsotopeCorrector::correctIsotopicImpurities(cm_in, cm_out, &quant_meth);

    // 1. check the actual result
    String cm_file_out;
    NEW_TMP_FILE(cm_file_out);
    cm_file.store(cm_file_out,cm_out);

    WHITELIST("<?xml-stylesheet,id=\",href=\"file:////");
    TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricIsotopeCorrector_out.consensusXML"));

    // 2. check the returned stats -> values are based on the org. impl.
    TEST_EQUAL(stats.channel_count, 4)
    TEST_EQUAL(stats.iso_number_ms2_negative, 8)
    TEST_EQUAL(stats.iso_number_reporter_negative, 9)
    TEST_EQUAL(stats.iso_number_reporter_different, 1)
    TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0.02489106611347)
    TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 559.034896850586)
    TEST_EQUAL(stats.number_ms2_total, cm_in.size())
    // TEST_EQUAL(stats.number_ms2_empty, 0)
    // TEST_EQUAL(stats.empty_channels[114], 0)
    // TEST_EQUAL(stats.empty_channels[115], 0)
    // TEST_EQUAL(stats.empty_channels[116], 0)
    // TEST_EQUAL(stats.empty_channels[117], 0)
  }

  // 3. check stats in detail
  {
    ConsensusXMLFile cm_file;
    ConsensusMap cm_in, cm_out;
    cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricIsotopeCorrector.consensusXML"),cm_in);
    cm_in.clear(false);

    // copy in/output
    cm_out = cm_in;


    // first run (empty):
    IsobaricQuantifierStatistics stats = IsobaricIsotopeCorrector::correctIsotopicImpurities(cm_in, cm_out, &quant_meth);
    TEST_EQUAL(stats.channel_count, 4)
    TEST_EQUAL(stats.iso_number_ms2_negative, 0)
    TEST_EQUAL(stats.iso_number_reporter_negative, 0)
    TEST_EQUAL(stats.iso_number_reporter_different, 0)
    TEST_REAL_SIMILAR(stats.iso_solution_different_intensity,  0)
    TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 0)
    TEST_EQUAL(stats.number_ms2_total, cm_in.size())
    // TEST_EQUAL(stats.number_ms2_empty, 0)
    // TEST_EQUAL(stats.empty_channels[114], 0)
    // TEST_EQUAL(stats.empty_channels[115], 0)
    // TEST_EQUAL(stats.empty_channels[116], 0)
    // TEST_EQUAL(stats.empty_channels[117], 0)


    // add some target results
    double v1[4] = {1.071,  95.341,  101.998,  96.900}; // naive yields: {-1,100,100,100};  NNLS: {0.00000  99.91414 100.00375  99.99990}
    cm_in.push_back(getCFWithIntensites(v1));
    cm_out = cm_in;

    stats = IsobaricIsotopeCorrector::correctIsotopicImpurities(cm_in, cm_out, &quant_meth);

    // check the corrected intensities
    ABORT_IF(cm_out[0].getFeatures().size() != 4)
    ConsensusFeature::HandleSetType::const_iterator it = cm_out[0].getFeatures().begin();
    TEST_REAL_SIMILAR(it->getIntensity(), 0.00000)
    ++it;
    TEST_REAL_SIMILAR(it->getIntensity(), 99.91414)
    ++it;
    TEST_REAL_SIMILAR(it->getIntensity(), 100.00375)
    ++it;
    TEST_REAL_SIMILAR(it->getIntensity(), 99.99990)

    // test the stats
    TEST_EQUAL(stats.channel_count, 4)
    TEST_EQUAL(stats.iso_number_ms2_negative, 1)
    TEST_EQUAL(stats.iso_number_reporter_negative, 1)
    TEST_EQUAL(stats.iso_number_reporter_different, 0)
    TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0)
    TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 299.9178)
    TEST_EQUAL(stats.number_ms2_total, cm_in.size())
    // TEST_EQUAL(stats.number_ms2_empty, 0)
    // TEST_EQUAL(stats.empty_channels[114], 1)
    // TEST_EQUAL(stats.empty_channels[115], 0)
    // TEST_EQUAL(stats.empty_channels[116], 0)
    // TEST_EQUAL(stats.empty_channels[117], 0)

    // change some more... (second run)
    double v2[4] = {0,0,0,0};
    cm_in.push_back(getCFWithIntensites(v2));
    cm_out = cm_in;
    stats = IsobaricIsotopeCorrector::correctIsotopicImpurities(cm_in, cm_out, &quant_meth);

    TEST_EQUAL(stats.channel_count, 4)
    TEST_EQUAL(stats.iso_number_ms2_negative, 1)
    TEST_EQUAL(stats.iso_number_reporter_negative, 1)
    TEST_EQUAL(stats.iso_number_reporter_different, 0)
    TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0)
    TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 299.9178)
    TEST_EQUAL(stats.number_ms2_total, cm_in.size())
    // TEST_EQUAL(stats.number_ms2_empty, 1)
    // TEST_EQUAL(stats.empty_channels[114], 2)
    // TEST_EQUAL(stats.empty_channels[115], 1)
    // TEST_EQUAL(stats.empty_channels[116], 1)
    // TEST_EQUAL(stats.empty_channels[117], 1)
  }

  // 4. test precondition
  {
    ConsensusXMLFile cm_file;
    ConsensusMap cm_in, cm_out;
    cm_file.load(OPENMS_GET_TEST_DATA_PATH("IsobaricIsotopeCorrector.consensusXML"),cm_in);

    TEST_PRECONDITION_VIOLATED(IsobaricIsotopeCorrector::correctIsotopicImpurities(cm_in,cm_out, &quant_meth))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
