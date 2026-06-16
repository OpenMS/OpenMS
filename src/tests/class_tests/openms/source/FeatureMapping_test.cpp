// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <vector>

///////////////////////////

START_TEST(FeatureMapping, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// Helper: build an MS2 spectrum with a given RT and (optional) precursor m/z.
auto mkMS2 = [](double rt, double prec_mz, bool has_prec = true) -> MSSpectrum
{
  MSSpectrum s;
  s.setMSLevel(2);
  s.setRT(rt);
  if (has_prec)
  {
    Precursor p;
    p.setMZ(prec_mz);
    s.setPrecursors(std::vector<Precursor>{p});
  }
  return s;
};

START_SECTION(([FeatureMapping::FeatureMappingInfo] default state))
{
  FeatureMapping::FeatureMappingInfo info;
  TEST_EQUAL(info.feature_maps.empty(), true)
}
END_SECTION

START_SECTION(([FeatureMapping::FeatureToMs2Indices] default state))
{
  FeatureMapping::FeatureToMs2Indices res;
  TEST_EQUAL(res.assignedMS2.empty(), true)
  TEST_EQUAL(res.unassignedMS2.empty(), true)
}
END_SECTION

START_SECTION((static FeatureToMs2Indices assignMS2IndexToFeature(const MSExperiment& spectra, const FeatureMappingInfo& fm_info, const double& precursor_mz_tolerance, const double& precursor_rt_tolerance, bool ppm)))
{
  // Two features at known (RT, m/z)
  FeatureMap fmap;
  { Feature f; f.setRT(100.0); f.setMZ(500.0); fmap.push_back(f); } // F0
  { Feature f; f.setRT(200.0); f.setMZ(600.0); fmap.push_back(f); } // F1

  FeatureMapping::FeatureMappingInfo info;
  info.feature_maps.push_back(fmap);
  info.kd_tree.addMaps(info.feature_maps); // references features in info.feature_maps

  const BaseFeature* F0 = &info.feature_maps[0][0];
  const BaseFeature* F1 = &info.feature_maps[0][1];

  // --- basic assignment, MS1 skipped, precursor-less MS2 dropped ---
  {
    MSExperiment exp;
    exp.addSpectrum(mkMS2(100.0, 500.0));         // idx0 -> F0
    exp.addSpectrum(mkMS2(200.0, 600.0));         // idx1 -> F1
    exp.addSpectrum(mkMS2(300.0, 999.0));         // idx2 -> unassigned (no feature in window)
    { MSSpectrum s; s.setMSLevel(1); s.setRT(150.0); exp.addSpectrum(s); } // idx3 MS1 -> skipped
    exp.addSpectrum(mkMS2(100.0, 0.0, false));    // idx4 no precursor -> dropped

    FeatureMapping::FeatureToMs2Indices r =
      FeatureMapping::assignMS2IndexToFeature(exp, info, 0.5, 5.0, false);

    TEST_EQUAL(r.assignedMS2.size(), 2)
    TEST_EQUAL(r.assignedMS2[F0].size(), 1)
    TEST_EQUAL(r.assignedMS2[F0][0], 0)
    TEST_EQUAL(r.assignedMS2[F1].size(), 1)
    TEST_EQUAL(r.assignedMS2[F1][0], 1)

    // idx2 had a precursor but no feature in window
    TEST_EQUAL(r.unassignedMS2.size(), 1)
    TEST_EQUAL(r.unassignedMS2[0], 2)
  }

  // --- RT outside tolerance -> unassigned ---
  {
    MSExperiment exp;
    exp.addSpectrum(mkMS2(120.0, 500.0)); // RT 120 vs F0 RT 100, rt_tol 5 -> out of window
    FeatureMapping::FeatureToMs2Indices r =
      FeatureMapping::assignMS2IndexToFeature(exp, info, 0.5, 5.0, false);
    TEST_EQUAL(r.assignedMS2.size(), 0)
    TEST_EQUAL(r.unassignedMS2.size(), 1)
    TEST_EQUAL(r.unassignedMS2[0], 0)
  }

  // --- ppm tolerance: in-window vs out-of-window ---
  {
    MSExperiment in_win;
    in_win.addSpectrum(mkMS2(200.0, 600.05)); // ~83 ppm from 600.0, inside 200 ppm
    FeatureMapping::FeatureToMs2Indices r =
      FeatureMapping::assignMS2IndexToFeature(in_win, info, 200.0, 5.0, true);
    TEST_EQUAL(r.assignedMS2.size(), 1)
    TEST_EQUAL(r.assignedMS2[F1].size(), 1)

    MSExperiment out_win;
    out_win.addSpectrum(mkMS2(200.0, 600.5)); // ~833 ppm from 600.0, outside 200 ppm
    FeatureMapping::FeatureToMs2Indices r2 =
      FeatureMapping::assignMS2IndexToFeature(out_win, info, 200.0, 5.0, true);
    TEST_EQUAL(r2.assignedMS2.size(), 0)
    TEST_EQUAL(r2.unassignedMS2.size(), 1)
  }
}
END_SECTION

START_SECTION([EXTRA] (multiple features in window: closest m/z wins))
{
  // Two features close in m/z; the precursor is nearer to G0
  FeatureMap fmap;
  { Feature f; f.setRT(100.0); f.setMZ(500.0); fmap.push_back(f); } // G0
  { Feature f; f.setRT(100.0); f.setMZ(500.3); fmap.push_back(f); } // G1

  FeatureMapping::FeatureMappingInfo info;
  info.feature_maps.push_back(fmap);
  info.kd_tree.addMaps(info.feature_maps);
  const BaseFeature* G0 = &info.feature_maps[0][0];
  const BaseFeature* G1 = &info.feature_maps[0][1];

  MSExperiment exp;
  exp.addSpectrum(mkMS2(100.0, 500.1)); // |500.1-500.0|=0.1 < |500.1-500.3|=0.2

  FeatureMapping::FeatureToMs2Indices r =
    FeatureMapping::assignMS2IndexToFeature(exp, info, 0.5, 5.0, false);

  TEST_EQUAL(r.assignedMS2.size(), 1)
  TEST_EQUAL(r.assignedMS2.count(G0), 1) // closest in m/z wins
  TEST_EQUAL(r.assignedMS2.count(G1), 0)
  TEST_EQUAL(r.unassignedMS2.empty(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
