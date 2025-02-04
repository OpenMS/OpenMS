// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MapAlignmentTransformer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentTransformer* ptr = nullptr;
MapAlignmentTransformer* null_ptr = nullptr;

TransformationDescription::DataPoints data;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(1.0, 3.0));

TransformationDescription td(data);
Param params;
td.fitModel("linear", params);

START_SECTION(MapAlignmentTransformer())
{
	ptr = new MapAlignmentTransformer();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MapAlignmentTransformer())
{
	delete ptr;
}
END_SECTION

START_SECTION((static void transformRetentionTimes(PeakMap& msexp, const TransformationDescription& trafo, bool store_original_rt = false)))
{
  PeakMap exp;
  PeakMap::SpectrumType spec;

  // first spectrum (MS)
  spec.setRT(11.1);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // second spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(11.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  // third spectrum (MS)
  spec.clear(true);
  spec.setRT(12.2);
  spec.setMSLevel(1);
  exp.addSpectrum(spec);

  // forth spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(12.5);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  MapAlignmentTransformer::transformRetentionTimes(exp, td);

  // check the spectra:
  TEST_EQUAL(exp[0].getRT(), 23.2)
  TEST_EQUAL(exp[1].getRT(), 24.0)
  TEST_EQUAL(exp[2].getRT(), 25.4)
  TEST_EQUAL(exp[3].getRT(), 26.0)

  // check storing of original RTs:
  for (Size i = 0; i < 4; ++i)
  {
    TEST_EQUAL(exp[i].metaValueExists("original_RT"), false);
  }

  MapAlignmentTransformer::transformRetentionTimes(exp, td, true);
  TEST_EQUAL(exp[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(exp[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(exp[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(exp[3].getMetaValue("original_RT"), 26.0);

  // applying a transform again doesn't overwrite the original RTs:
  MapAlignmentTransformer::transformRetentionTimes(exp, td, true);
  TEST_EQUAL(exp[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(exp[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(exp[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(exp[3].getMetaValue("original_RT"), 26.0);
}
END_SECTION

START_SECTION((static void transformRetentionTimes(FeatureMap& fmap, const TransformationDescription& trafo, bool store_original_rt = false)))
{
  Feature f;
  FeatureMap featmap;

  f.setRT(11.1);
  featmap.push_back(f);

  f.setRT(11.5);
  featmap.push_back(f);

  f.setRT(12.2);
  featmap.push_back(f);

  f.setRT(12.5);
  featmap.push_back(f);


  MapAlignmentTransformer::transformRetentionTimes(featmap, td);

  // check the features:
  TEST_EQUAL(featmap[0].getRT(), 23.2)
  TEST_EQUAL(featmap[1].getRT(), 24.0)
  TEST_EQUAL(featmap[2].getRT(), 25.4)
  TEST_EQUAL(featmap[3].getRT(), 26.0)

  // check storing of original RTs:
  for (Size i = 0; i < 4; ++i)
  {
    TEST_EQUAL(featmap[i].metaValueExists("original_RT"), false);
  }

  MapAlignmentTransformer::transformRetentionTimes(featmap, td, true);
  TEST_EQUAL(featmap[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(featmap[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(featmap[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(featmap[3].getMetaValue("original_RT"), 26.0);

  // applying a transform again doesn't overwrite the original RTs:
  MapAlignmentTransformer::transformRetentionTimes(featmap, td, true);
  TEST_EQUAL(featmap[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(featmap[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(featmap[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(featmap[3].getMetaValue("original_RT"), 26.0);
}
END_SECTION

START_SECTION((static void transformRetentionTimes(ConsensusMap& cmap, const TransformationDescription& trafo, bool store_original_rt = false)))
{
  ConsensusFeature cf;
  ConsensusMap consensusmap;

  cf.setRT(11.1);
  consensusmap.push_back(cf);

  cf.setRT(11.5);
  consensusmap.push_back(cf);

  cf.setRT(12.2);
  consensusmap.push_back(cf);

  cf.setRT(12.5);
  consensusmap.push_back(cf);

  MapAlignmentTransformer::transformRetentionTimes(consensusmap, td);

  // check the consensus features:
  TEST_EQUAL(consensusmap[0].getRT(), 23.2)
  TEST_EQUAL(consensusmap[1].getRT(), 24.0)
  TEST_EQUAL(consensusmap[2].getRT(), 25.4)
  TEST_EQUAL(consensusmap[3].getRT(), 26.0)

  // check storing of original RTs:
  for (Size i = 0; i < 4; ++i)
  {
    TEST_EQUAL(consensusmap[i].metaValueExists("original_RT"), false);
  }

  MapAlignmentTransformer::transformRetentionTimes(consensusmap, td, true);
  TEST_EQUAL(consensusmap[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(consensusmap[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(consensusmap[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(consensusmap[3].getMetaValue("original_RT"), 26.0);

  // applying a transform again doesn't overwrite the original RTs:
  MapAlignmentTransformer::transformRetentionTimes(consensusmap, td, true);
  TEST_EQUAL(consensusmap[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(consensusmap[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(consensusmap[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(consensusmap[3].getMetaValue("original_RT"), 26.0);
}
END_SECTION

START_SECTION((static void transformRetentionTimes(std::vector<PeptideIdentification>& pep_ids, const TransformationDescription& trafo, bool store_original_rt = false)))
{
  PeptideIdentification pi;
  vector<PeptideIdentification> pis;

  pi.setRT(11.1);
  pis.push_back(pi);

  pi.setRT(11.5);
  pis.push_back(pi);

  pi.setRT(12.2);
  pis.push_back(pi);

  pi.setRT(12.5);
  pis.push_back(pi);

  MapAlignmentTransformer::transformRetentionTimes(pis, td);

  // check the peptide IDs:
  TEST_EQUAL(pis[0].getRT(), 23.2)
  TEST_EQUAL(pis[1].getRT(), 24.0)
  TEST_EQUAL(pis[2].getRT(), 25.4)
  TEST_EQUAL(pis[3].getRT(), 26.0)

  // check storing of original RTs:
  for (Size i = 0; i < 4; ++i)
  {
    TEST_EQUAL(pis[i].metaValueExists("original_RT"), false);
  }

  MapAlignmentTransformer::transformRetentionTimes(pis, td, true);
  TEST_EQUAL(pis[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(pis[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(pis[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(pis[3].getMetaValue("original_RT"), 26.0);

  // applying a transform again doesn't overwrite the original RTs:
  MapAlignmentTransformer::transformRetentionTimes(pis, td, true);
  TEST_EQUAL(pis[0].getMetaValue("original_RT"), 23.2);
  TEST_EQUAL(pis[1].getMetaValue("original_RT"), 24.0);
  TEST_EQUAL(pis[2].getMetaValue("original_RT"), 25.4);
  TEST_EQUAL(pis[3].getMetaValue("original_RT"), 26.0);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



