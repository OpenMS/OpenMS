// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <set>

using namespace OpenMS;

START_TEST(FeatureFinderAlgorithmMetaboIdent, "$Id$")

START_SECTION(([EXTRA] void run(const std::vector<FeatureFinderMetaboIdentCompound>& metaboIdentTable, FeatureMap& features, const std::string& spectra_file) with FAIMS data))
{
  // the data of TOPP_FeatureFinderMetaboIdent_1, its spectra alternately at two compensation voltages
  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("../../../topp/FeatureFinderMetaboIdent_1_input.mzML"), exp);
  for (Size i = 0; i < exp.size(); ++i)
  {
    exp[i].setDriftTime(i % 2 ? -65.0 : -45.0);
    exp[i].setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  }

  FeatureFinderAlgorithmMetaboIdent ff;
  Param p = ff.getParameters();
  p.setValue("extract:mz_window", 5.0);
  p.setValue("extract:rt_window", 20.0);
  p.setValue("detect:peak_width", 3.0);
  ff.setParameters(p);
  ff.setMSData(std::move(exp));

  using Compound = FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound;
  std::vector<Compound> targets;
  targets.emplace_back("inosine", "C10H12N4O5", 0.0, std::vector<int> {1}, std::vector<double> {264.0}, std::vector<double> {0.0},
                       std::vector<double> {0.0});
  // a mass without signal
  targets.emplace_back("missing", "C30H50N4O5", 0.0, std::vector<int> {1}, std::vector<double> {264.0}, std::vector<double> {0.0},
                       std::vector<double> {0.0});
  FeatureMap features;
  ff.run(targets, features, "test.mzML");

  // the targets without features of both groups: identifications without matches of one run "id", which no feature
  // links, each with its compensation voltage
  const auto& data = features.getIdentificationData();
  ABORT_IF(data.getRuns().size() != 1)
  TEST_EQUAL(data.getRuns()[0].getIdentifier(), "id")
  TEST_EQUAL(data.getRuns()[0].getNumberOfMatches(), 0)
  std::set<double> cvs;
  for (const auto& target : features.getUnassignedIdentifications())
  {
    if (target.query->getMetaValue("label") != DataValue("missing")) continue;
    TEST_EQUAL(target.query->getMetaValue("PeptideRef").toString().empty(), false)
    TEST_EQUAL(bool(target.query->rt), true)
    TEST_EQUAL(bool(target.query->mz), true)
    cvs.insert(double(target.query->getMetaValue(Constants::UserParam::FAIMS_CV)));
  }
  TEST_EQUAL(cvs.size(), 2)
  TEST_EQUAL(cvs.contains(-45.0), true)
  TEST_EQUAL(cvs.contains(-65.0), true)
}
END_SECTION

END_TEST
