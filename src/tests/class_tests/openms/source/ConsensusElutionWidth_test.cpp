// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ConsensusElutionWidth.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <cmath>

using namespace OpenMS;
using namespace std;

namespace
{
  const double SQRT_2PI = std::sqrt(2.0 * Constants::PI);

  Feature makeFeature(UInt64 uid, double height, double sigma, double center, bool identified, bool fit_ok = true)
  {
    Feature f;
    f.setUniqueId(uid);
    f.setRT(center);
    f.setIntensity(static_cast<float>(height * sigma * SQRT_2PI));
    f.setMetaValue("model_height", height);
    f.setMetaValue("model_Gauss_sigma", sigma);
    f.setMetaValue("model_center", center);
    f.setMetaValue("model_status", fit_ok ? "0 (valid)" : "1 (invalid area)");
    if (identified) f.getPeptideIdentifications().resize(1);
    return f;
  }

  TransformationDescription linearTrafo(double slope)
  {
    TransformationDescription::DataPoints points;
    points.emplace_back(0.0, 0.0);
    points.emplace_back(1000.0, 1000.0 * slope);
    TransformationDescription trafo(points);
    trafo.fitModel("linear");
    return trafo;
  }
}

START_TEST(ConsensusElutionWidth, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static double localSlope(const TransformationDescription& trafo, double rt)))
{
  TEST_REAL_SIMILAR(ConsensusElutionWidth::localSlope(TransformationDescription(), 100.0), 1.0)
  TEST_REAL_SIMILAR(ConsensusElutionWidth::localSlope(linearTrafo(0.8), 500.0), 0.8)
}
END_SECTION

START_SECTION((static Stats apply(ConsensusMap& consensus, const std::vector<FeatureMap>& feature_maps, const std::vector<TransformationDescription>& transformations)))
{
  // One precursor in three runs with the same true abundance (height 1e6, width 5 s on the consensus scale).
  // Run 1's RT is stretched by 1.25 relative to the consensus scale (trafo slope 0.8, true width 6.25 s, height 0.8e6).
  // The width fits are off by factors 0.7, 1.4 (identified) and the fit of run 2 (unidentified) by 1.2.
  std::vector<FeatureMap> maps(3);
  maps[0].push_back(makeFeature(10, 1e6, 5.0 * 0.7, 500.0, true));
  maps[1].push_back(makeFeature(11, 0.8e6, 6.25 * 1.4, 625.0, true));
  maps[2].push_back(makeFeature(12, 1e6, 5.0 * 1.2, 500.0, false));
  // a second precursor without identified fits: unchanged
  maps[0].push_back(makeFeature(20, 2e6, 4.0, 900.0, false));
  maps[1].push_back(makeFeature(21, 2e6, 3.0, 900.0, true, false));
  std::vector<TransformationDescription> trafos = {TransformationDescription(), linearTrafo(0.8)};

  ConsensusMap cmap;
  ConsensusFeature cf1;
  for (Size m = 0; m < 3; ++m) cf1.insert(m, maps[m][0]);
  cmap.push_back(cf1);
  ConsensusFeature cf2;
  cf2.insert(0, maps[0][1]);
  cf2.insert(1, maps[1][1]);
  cmap.push_back(cf2);

  auto stats = ConsensusElutionWidth::apply(cmap, maps, trafos);
  TEST_EQUAL(stats.consensus_features, 1)
  TEST_EQUAL(stats.handles_updated, 3)
  TEST_EQUAL(stats.handles_kept, 2)

  // consensus width: median of 3.5 * 1 and 8.75 * 0.8 = 7 -> 5.25; each run reports height x 5.25 / slope x sqrt(2 pi)
  const double w_c = 0.5 * (3.5 + 7.0);
  std::vector<double> q;
  for (const auto& fh : cmap[0].getFeatures()) q.push_back(fh.getIntensity());
  TOLERANCE_RELATIVE(1.0001)
  TEST_REAL_SIMILAR(q[0], 1e6 * w_c * SQRT_2PI)
  TEST_REAL_SIMILAR(q[1], 0.8e6 * (w_c / 0.8) * SQRT_2PI)
  TEST_REAL_SIMILAR(q[2], 1e6 * w_c * SQRT_2PI)
  // the run-to-run difference of the width fits is gone
  TEST_REAL_SIMILAR(q[0], q[1])

  // second precursor: original intensities
  std::vector<double> q2;
  for (const auto& fh : cmap[1].getFeatures()) q2.push_back(fh.getIntensity());
  TEST_REAL_SIMILAR(q2[0], 2e6 * 4.0 * SQRT_2PI)
  TEST_REAL_SIMILAR(q2[1], 2e6 * 3.0 * SQRT_2PI)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
