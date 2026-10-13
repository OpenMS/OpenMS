// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/NativeIdentificationTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include "RTAlignment.h"
///////////////////////////

START_TEST(RTAlignment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

RTAlignment* ptr = nullptr;
RTAlignment* nullPointer = nullptr;

START_SECTION(RTAlignment())
{
  ptr = new RTAlignment;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(~RTAlignment())
{
  delete ptr;
}
END_SECTION

RTAlignment rtA;
START_SECTION(QCBase::Status requirements() const override)
{
  TEST_EQUAL(rtA.requirements() == (QCBase::Status() | QCBase::Requires::TRAFOALIGN | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(const std::string& getName() const override) {TEST_EQUAL(rtA.getName(), "RTAlignment")} END_SECTION

  START_SECTION((compute(FeatureMap & features, TransformationDescription& trafo)))
{
  // Valid FeatureMap
  FeatureMap fmap;
  PeptideIdentification peptide_ID;
  PeptideIdentificationList identifications;
  PeptideIdentificationList unassignedIDs;
  Feature feature1, feature2;
  peptide_ID.setRT(0);
  identifications.push_back(peptide_ID);
  peptide_ID.setRT(1);
  identifications.push_back(peptide_ID);
  MapPeptideIdentifications ids;
  ids.features[{0}] = identifications;
  identifications.clear();
  fmap.push_back(feature1);
  peptide_ID.setRT(10);
  identifications.push_back(peptide_ID);
  peptide_ID.setRT(12);
  identifications.push_back(peptide_ID);
  ids.features[{1}] = identifications;
  fmap.push_back(feature2);
  // unassigned PeptideHits
  peptide_ID.setRT(0.5);
  unassignedIDs.push_back(peptide_ID);
  peptide_ID.setRT(2.5);
  unassignedIDs.push_back(peptide_ID);
  ids.unassigned = unassignedIDs;

  // Transformation
  TransformationDescription td;
  td.fitModel("identity", Param());
  td.setDataPoints(vector<pair<double, double>> {{0.0, 1.0}, {0.25, 1.5}, {0.5, 2.0}, {1.0, 3.0}});
  td.fitModel("linear");
  RTAlignment rtA;
  Internal::ClassTest::addSearchRun(ids);
  Internal::ClassTest::setIdentifications(fmap, ids);
  rtA.compute(fmap, td);
  const auto computed = Internal::ClassTest::identificationsOf(fmap);
  // test features
  TEST_REAL_SIMILAR(computed.at({0})[0].getMetaValue("rt_align"), 1);
  TEST_REAL_SIMILAR(computed.at({0})[0].getMetaValue("rt_raw"), 0);
  TEST_REAL_SIMILAR(computed.at({0})[1].getMetaValue("rt_align"), 3);
  TEST_REAL_SIMILAR(computed.at({0})[1].getMetaValue("rt_raw"), 1);
  TEST_REAL_SIMILAR(computed.at({1})[0].getMetaValue("rt_align"), 21);
  TEST_REAL_SIMILAR(computed.at({1})[0].getMetaValue("rt_raw"), 10);
  // test unassigned
  TEST_REAL_SIMILAR(computed.unassigned[0].getMetaValue("rt_align"), 2);
  TEST_REAL_SIMILAR(computed.unassigned[0].getMetaValue("rt_raw"), 0.5);
  TEST_REAL_SIMILAR(computed.unassigned[1].getMetaValue("rt_align"), 6);
  TEST_REAL_SIMILAR(computed.unassigned[1].getMetaValue("rt_raw"), 2.5);

  // empty FeatureMap
  FeatureMap fmap_empty {};
  rtA.compute(fmap_empty, td);
  // empty feature
  Feature feature_empty {};
  fmap_empty.push_back(feature_empty);
  rtA.compute(fmap_empty, td);
  // empty PeptideIdentifications
  fmap_empty.push_back(feature1);
  rtA.compute(fmap_empty, td);
  // empty PeptideIdentification
  PeptideIdentification peptide_ID_empty {};
  fmap_empty.push_back(feature1);
  MapPeptideIdentifications empty_ids;
  empty_ids.features[{fmap_empty.size() - 1}] = {peptide_ID_empty};
  Internal::ClassTest::addSearchRun(empty_ids);
  Internal::ClassTest::setIdentifications(fmap_empty, empty_ids);
  rtA.compute(fmap_empty, td);

  // data processing: after alignment
  DataProcessing processing_method {};
  std::set<DataProcessing::ProcessingAction> dp {DataProcessing::ProcessingAction::ALIGNMENT};
  processing_method.setProcessingActions(dp);
  fmap.setDataProcessing({processing_method});
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, rtA.compute(fmap, td), "Metric RTAlignment received a featureXML AFTER map alignment, but needs a featureXML BEFORE map alignment!");
}
END_SECTION

START_SECTION((compute(std::vector<IdentificationData::Observation>& observations, const TransformationDescription& trafo)))
{
  TransformationDescription td;
  td.setDataPoints(vector<pair<double, double>> {{0.0, 1.0}, {0.25, 1.5}, {0.5, 2.0}, {1.0, 3.0}});
  td.fitModel("linear");
  std::vector<IdentificationData::Observation> observations(2);
  observations[0].rt = 0.5;
  // an observation without retention time gets no values
  RTAlignment().compute(observations, td);
  TEST_REAL_SIMILAR(observations[0].getMetaValue("rt_align"), 2);
  TEST_REAL_SIMILAR(observations[0].getMetaValue("rt_raw"), 0.5);
  TEST_EQUAL(observations[1].metaValueExists("rt_align"), false)
  TEST_EQUAL(observations[1].metaValueExists("rt_raw"), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
