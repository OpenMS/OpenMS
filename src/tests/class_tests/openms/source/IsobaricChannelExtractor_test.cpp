// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricChannelExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricChannelExtractor * ptr = nullptr;
IsobaricChannelExtractor* null_ptr = nullptr;
IsobaricQuantitationMethod* q_method = new ItraqFourPlexQuantitationMethod();

START_SECTION((IsobaricChannelExtractor(const IsobaricQuantitationMethod * const quant_method)))
{
  ptr = new IsobaricChannelExtractor(q_method);
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricChannelExtractor())
{
  delete ptr;
}

END_SECTION

START_SECTION((IsobaricChannelExtractor(const IsobaricChannelExtractor &other)))
{
  IsobaricChannelExtractor ice(q_method);
  Param p = ice.getParameters();
  p.setValue("select_activation", "any");

  ice.setParameters(p);

  IsobaricChannelExtractor ice2(ice);
  TEST_EQUAL(ice2.getParameters(), p)
}

END_SECTION

START_SECTION((IsobaricChannelExtractor& operator=(const IsobaricChannelExtractor& rhs)))
{
  IsobaricChannelExtractor ice(q_method);
  Param p = ice.getParameters();
  p.setValue("reporter_mass_shift", 0.3);
  ice.setParameters(p);

  IsobaricChannelExtractor ice2(q_method);
  ice2 = ice;
  TEST_EQUAL(ice2.getParameters(), p)
}

END_SECTION

START_SECTION((void extractChannels(const PeakMap&ms_exp_data, ConsensusMap & consensus_map)))
{
  {
    // load test data
    PeakMap exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    // check channel meta information
    TEST_EQUAL(cm_out.getColumnHeaders().size(), 4)
    ABORT_IF(cm_out.getColumnHeaders().size() != 4)

    TEST_EQUAL(cm_out.getColumnHeaders()[0].label, "itraq4plex_114")
    TEST_EQUAL(cm_out.getColumnHeaders()[0].getMetaValue("channel_name"), "114")
    TEST_EQUAL(cm_out.getColumnHeaders()[0].getMetaValue("channel_id"), 0)
    TEST_EQUAL(cm_out.getColumnHeaders()[0].getMetaValue("channel_description"), "ref")
    TEST_EQUAL(cm_out.getColumnHeaders()[0].getMetaValue("channel_center"), 114.1112)

    TEST_EQUAL(cm_out.getColumnHeaders()[1].label, "itraq4plex_115")
    TEST_EQUAL(cm_out.getColumnHeaders()[1].getMetaValue("channel_name"), "115")
    TEST_EQUAL(cm_out.getColumnHeaders()[1].getMetaValue("channel_id"), 1)
    TEST_EQUAL(cm_out.getColumnHeaders()[1].getMetaValue("channel_description"), "something")
    TEST_EQUAL(cm_out.getColumnHeaders()[1].getMetaValue("channel_center"), 115.1082)

    TEST_EQUAL(cm_out.getColumnHeaders()[2].label, "itraq4plex_116")
    TEST_EQUAL(cm_out.getColumnHeaders()[2].getMetaValue("channel_name"), "116")
    TEST_EQUAL(cm_out.getColumnHeaders()[2].getMetaValue("channel_id"), 2)
    TEST_EQUAL(cm_out.getColumnHeaders()[2].getMetaValue("channel_description"), "else")
    TEST_EQUAL(cm_out.getColumnHeaders()[2].getMetaValue("channel_center"), 116.1116)

    TEST_EQUAL(cm_out.getColumnHeaders()[3].label, "itraq4plex_117")
    TEST_EQUAL(cm_out.getColumnHeaders()[3].getMetaValue("channel_name"), "117")
    TEST_EQUAL(cm_out.getColumnHeaders()[3].getMetaValue("channel_id"), 3)
    TEST_EQUAL(cm_out.getColumnHeaders()[3].getMetaValue("channel_description"), "")
    TEST_EQUAL(cm_out.getColumnHeaders()[3].getMetaValue("channel_center"), 117.1149)

    // compare results
    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)
    ConsensusFeature::iterator cf_it;

    TEST_EQUAL(cm_out[0].size(), 4)
    TEST_EQUAL(cm_out[0].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=2")
    TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_intensity"), 5251952.5)
    TEST_EQUAL(cm_out[0].getCharge(), 2)
    TEST_REAL_SIMILAR(cm_out[0].getIntensity(), 1490501.21)
    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 643005.56)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 458708.97)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 182238.38)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 206543.3)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())


    TEST_EQUAL(cm_out[1].size(), 4)
    TEST_EQUAL(cm_out[1].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=4")
    TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_intensity"), 7365030)
    TEST_EQUAL(cm_out[1].getCharge(), 3)
    TEST_REAL_SIMILAR(cm_out[1].getIntensity(), 2329603)
    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 847251)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 861806)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 311899)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    TEST_EQUAL(cm_out[2].size(), 4)
    TEST_EQUAL(cm_out[2].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=6")
    TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_intensity"), 6835636)
    TEST_EQUAL(cm_out[2].getCharge(), 3)
    TEST_REAL_SIMILAR(cm_out[2].getIntensity(), 2520967)
    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 894414)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 958965)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 326443)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())

    TEST_EQUAL(cm_out[3].size(), 4)
    TEST_EQUAL(cm_out[3].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=8")
    TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_intensity"), 6762358)
    TEST_EQUAL(cm_out[3].getCharge(), 3)
    TEST_REAL_SIMILAR(cm_out[3].getIntensity(), 1585286)
    cf_it = cm_out[3].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 581601)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 623851)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 191352)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 188482)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[3].end())

    TEST_EQUAL(cm_out[4].size(), 4)
    TEST_EQUAL(cm_out[4].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=10")
    TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_intensity"), 5464634.5)
    TEST_EQUAL(cm_out[4].getCharge(), 2)
    TEST_REAL_SIMILAR(cm_out[4].getIntensity(), 1746368)
    cf_it = cm_out[4].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[4].end())
  }
  { // test -> keep_unannotated_precursor
    // load test data
    PeakMap exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_7.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    p.setValue("keep_unannotated_precursor", "false");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 4)
    ABORT_IF(cm_out.size() != 4)
    TEST_EQUAL(((double)cm_out[0].getMetaValue("precursor_intensity")) != 0.0, true)

    p.setValue("keep_unannotated_precursor", "true");
    ice.setParameters(p);

    ConsensusMap cm_out_w_unannotated;
    ice.extractChannels(exp, cm_out_w_unannotated);

    TEST_EQUAL(cm_out_w_unannotated.size(), 5)
    ABORT_IF(cm_out_w_unannotated.size() != 5)
    TEST_REAL_SIMILAR(cm_out_w_unannotated[0].getMetaValue("precursor_intensity"), 0.0)
  }
  {
    // load test data
    PeakMap exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    p.setValue("min_precursor_intensity", 5300000.0);
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    // compare results
    TEST_EQUAL(cm_out.size(), 4)
    ABORT_IF(cm_out.size() != 4)
    for(ConsensusMap::Iterator cf = cm_out.begin(); cf != cm_out.end(); ++cf)
    {
      double prec_intensity = cf->getMetaValue("precursor_intensity");
      TEST_EQUAL(prec_intensity > 5300000.0, true)
    }
  }
  {
    // load test data
    PeakMap exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    p.setValue("min_reporter_intensity", 200000.0);
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ConsensusFeature::iterator cf_it;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)

    TEST_EQUAL(cm_out[0].size(), 4)
    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 643005.56)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 458708.97)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0) // is 182238.38 < 200.000
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 206543.3)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())

    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 847251)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 861806)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 311899)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 894414)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 958965)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 326443)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())

    cf_it = cm_out[3].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 581601)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 623851)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0) // is 191352 < 200.000
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0) // is 188482 < 200.000
    ++cf_it;
    ABORT_IF(cf_it != cm_out[3].end())

    cf_it = cm_out[4].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[4].end())
  }
  {
    // load test data
    PeakMap exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    p.setValue("min_reporter_intensity", 200000.0);
    p.setValue("discard_low_intensity_quantifications", "true");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ConsensusFeature::iterator cf_it;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 3)
    ABORT_IF(cm_out.size() != 3)

    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 847251)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 861806)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 311899)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())

    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 894414)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 958965)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 326443)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())
  }
  {
    // check precursor purity computation
    // - tested purities were validated manually
    // - dataset contains 2 ms1 and 5 ms2 spectra
    //   with the purity values listed below

    PeakMap exp_purity;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp_purity);

    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "any");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp_purity, cm_out);

    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)

    // check results
    TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_purity"), 1.0)
    TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_purity"), 0.692434)
    TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_purity"), 0.824561)
    TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_purity"), 0.731295)
    TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_purity"), 1.0)

    // now filter by purity
    p.setValue("min_precursor_purity", 0.75);
    ice.setParameters(p);

    ConsensusMap cm_filtered;
    ice.extractChannels(exp_purity, cm_filtered);

    TEST_EQUAL(cm_filtered.size(), 3)
    ABORT_IF(cm_filtered.size() != 3)

    // check results
    TEST_REAL_SIMILAR(cm_filtered[0].getMetaValue("precursor_purity"), 1.0)
    TEST_REAL_SIMILAR(cm_filtered[1].getMetaValue("precursor_purity"), 0.824561)
    TEST_REAL_SIMILAR(cm_filtered[2].getMetaValue("precursor_purity"), 1.0)
  }
}
END_SECTION

START_SECTION(([EXTRA] purity computation without interpolation))
{
  // check precursor purity computation
  // - tested purities were validated manually
  // - dataset contains 2 ms1 and 5 ms2 spectra
  //   with the purity values listed below

  PeakMap exp_purity;
  MzMLFile mzmlfile;
  mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp_purity);

  Param pItraq = q_method->getParameters();
  pItraq.setValue("channel_114_description", "ref");
  pItraq.setValue("channel_115_description", "something");
  pItraq.setValue("channel_116_description", "else");
  q_method->setParameters(pItraq);

  IsobaricChannelExtractor ice(q_method);

  // disable activation filtering
  Param p = ice.getParameters();
  p.setValue("select_activation", "any");
  p.setValue("purity_interpolation", "false");

  ice.setParameters(p);

  // extract channels
  ConsensusMap cm_out;
  ice.extractChannels(exp_purity, cm_out);

  TEST_EQUAL(cm_out.size(), 5)
  ABORT_IF(cm_out.size() != 5)

  // check results
  TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_purity"), 1.0)
  TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_purity"), 0.65472)
  TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_purity"), 0.775739)
  TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_purity"), 0.72009)
  TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_purity"), 1.0)

  // now filter by purity
  p.setValue("min_precursor_purity", 0.75);
  ice.setParameters(p);

  ConsensusMap cm_filtered;
  ice.extractChannels(exp_purity, cm_filtered);

  TEST_EQUAL(cm_filtered.size(), 3)
  ABORT_IF(cm_filtered.size() != 3)

  // check results
  TEST_REAL_SIMILAR(cm_filtered[0].getMetaValue("precursor_purity"), 1.0)
  TEST_REAL_SIMILAR(cm_filtered[1].getMetaValue("precursor_purity"), 0.775739)
  TEST_REAL_SIMILAR(cm_filtered[2].getMetaValue("precursor_purity"), 1.0)
}
END_SECTION

// extra test for tmt10plex to ensure high-res extraction works
START_SECTION(([EXTRA] TMT 10plex support))
{
  PeakMap tmt10plex_exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_8.mzML"), tmt10plex_exp);

  TMTTenPlexQuantitationMethod tmt10plex;
  IsobaricChannelExtractor ice(&tmt10plex);

  // disable activation filtering
  Param p = ice.getParameters();
  p.setValue("reporter_mass_shift", 0.003);
  ice.setParameters(p);

  // extract channels
  ConsensusMap cm_out;
  ice.extractChannels(tmt10plex_exp, cm_out);

  TEST_EQUAL(cm_out.size(), 5)
  ABORT_IF(cm_out.size() != 5)

  ConsensusMap::iterator cm_it = cm_out.begin();
  ConsensusFeature::iterator cf_it;

  TEST_EQUAL(cm_it->size(), 10)
  ABORT_IF(cm_it->size() != 10)
  TEST_EQUAL(cm_it->getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=7811")

  // test the extracted intensities
  cf_it = cm_it->begin();

  TEST_REAL_SIMILAR(cf_it->getIntensity(), 7759.65)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 6637.34)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 9147.74)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 8026)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 9454.86)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 21048.8)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 27783.1)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 27442.5)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 15765.4)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 17543.5)
  ++cf_it;
  ABORT_IF(cf_it != cm_it->end())

  // next scan
  ++cm_it;
  TEST_EQUAL(cm_it->size(), 10)
  ABORT_IF(cm_it->size() != 10)
  TEST_EQUAL(cm_it->getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=7812")

  // test the extracted intensities .. acutally no reporter in this scan
  cf_it = cm_it->begin();

  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  ABORT_IF(cf_it != cm_it->end())

  // next scan
  ++cm_it;
  TEST_EQUAL(cm_it->size(), 10)
  ABORT_IF(cm_it->size() != 10)
  TEST_EQUAL(cm_it->getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=7813")

  // test the extracted intensities
  cf_it = cm_it->begin();

  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1888.23)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1692.61)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1902.28)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1234.26)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1961.36)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 1560.74)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 2308.15)
  ++cf_it;
  ABORT_IF(cf_it != cm_it->end())

  // next scan
  ++cm_it;
  TEST_EQUAL(cm_it->size(), 10)
  ABORT_IF(cm_it->size() != 10)
  TEST_EQUAL(cm_it->getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=7814")

  // test the extracted intensities
  cf_it = cm_it->begin();

  TEST_REAL_SIMILAR(cf_it->getIntensity(), 26266.6)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 20802.2)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 36053.4)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 30815.4)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 34762.3)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 27767.8)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 45284.8)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 51015.2)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 29435.1)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 40080.7)
  ++cf_it;
  ABORT_IF(cf_it != cm_it->end())

  // next scan
  ++cm_it;
  TEST_EQUAL(cm_it->size(), 10)
  ABORT_IF(cm_it->size() != 10)
  TEST_EQUAL(cm_it->getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=7815")

  // test the extracted intensities
  cf_it = cm_it->begin();

  TEST_REAL_SIMILAR(cf_it->getIntensity(), 30760.9)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 17172)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 19647.1)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 24401.9)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 32279.3)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 19115.6)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 35027.3)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 34874.2)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 24060.4)
  ++cf_it;
  TEST_REAL_SIMILAR(cf_it->getIntensity(), 30866.5)
  ++cf_it;
  ABORT_IF(cf_it != cm_it->end())
}
END_SECTION

delete q_method;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
