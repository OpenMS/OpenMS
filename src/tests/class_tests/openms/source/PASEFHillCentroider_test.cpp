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

#include <OpenMS/FORMAT/HANDLERS/PASEFHillCentroider.h>

#include <algorithm>
#include <cmath>
#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal::PASEFHillCentroider;
using namespace std;

namespace
{
  // Build flat (mz, intensity, im) arrays from a list of per-IM-scan peak
  // lists. Each scan is a vector of (mz, intensity); im_value is the same for
  // all peaks within one scan.
  struct ScanPeaks
  {
    double im;
    vector<pair<double, double>> peaks; // (mz, intensity)
  };

  void buildFlat(const vector<ScanPeaks>& scans,
                 vector<double>& mz, vector<double>& intensity, vector<double>& im)
  {
    mz.clear(); intensity.clear(); im.clear();
    for (const auto& s : scans)
    {
      for (const auto& [m, i] : s.peaks)
      {
        mz.push_back(m);
        intensity.push_back(i);
        im.push_back(s.im);
      }
    }
  }

  Params defaultParams()
  {
    Params p;
    p.mz_tol_ppm        = 20.0;  // generous: synthetic data uses round numbers
    p.valley_factor     = 1.3;
    p.min_hill_length   = 2;
    p.min_intensity     = 0.0;
    p.im_group_tolerance = 0.0;
    return p;
  }
}

START_TEST(PASEFHillCentroider, "$Id$")

START_SECTION([EXTRA] empty input returns empty output)
{
  auto result = centroidFrame(nullptr, nullptr, nullptr, 0, defaultParams());
  TEST_EQUAL(result.empty(), true);
}
END_SECTION

START_SECTION([EXTRA] isolated peak across N IM scans yields one centroid)
{
  // One ion at m/z=500, present in IM scans 1.0, 1.05, 1.10, 1.15, 1.20
  // with a Gaussian-ish intensity profile peaking at the middle scan.
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}}},
    {1.05, {{500.0, 400.0}}},
    {1.10, {{500.0, 800.0}}}, // apex
    {1.15, {{500.0, 400.0}}},
    {1.20, {{500.0, 100.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  TEST_EQUAL(result.size(), 1);
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
  TEST_REAL_SIMILAR(result[0].intensity, 100.0 + 400.0 + 800.0 + 400.0 + 100.0);
  TEST_REAL_SIMILAR(result[0].im_apex, 1.10);
  TEST_REAL_SIMILAR(result[0].im_lower, 1.00);
  TEST_REAL_SIMILAR(result[0].im_upper, 1.20);
  TEST_EQUAL(result[0].length, 5);
}
END_SECTION

START_SECTION([EXTRA] two co-m/z ions separated by IM valley are split)
{
  // Same m/z=500 across 9 IM scans, two peaks separated by a deep valley.
  // Intensities: 200, 700, 200, 80 (valley), 200, 700, 200 — but with two
  // hills, each long enough to survive min_hill_length=2 after splitting.
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 200.0}}},
    {1.05, {{500.0, 700.0}}}, // first apex
    {1.10, {{500.0, 200.0}}},
    {1.15, {{500.0, 50.0}}},  // valley (700/50 = 14 >> hvf=1.3)
    {1.20, {{500.0, 200.0}}},
    {1.25, {{500.0, 700.0}}}, // second apex
    {1.30, {{500.0, 200.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  // Expect exactly two centroids — the headline win for hill detection.
  TEST_EQUAL(result.size(), 2);
  // Both centroids share the same m/z (intensity-weighted mean is 500).
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
  TEST_REAL_SIMILAR(result[1].mz, 500.0);
  // Apex IM values flank the valley.
  TEST_EQUAL(result[0].im_apex < 1.15, true);
  TEST_EQUAL(result[1].im_apex > 1.15, true);
}
END_SECTION

START_SECTION([EXTRA] single-IM-scan blip is rejected by min_hill_length=2)
{
  // One real ion + one single-scan noise blip at a different m/z.
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}}},
    {1.05, {{500.0, 400.0}}},
    {1.10, {{500.0, 800.0}}},
    {1.15, {{500.0, 400.0}}},
    {1.20, {{500.0, 100.0}, {600.0, 999.0}}}, // 600 only in this scan
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  TEST_EQUAL(result.size(), 1);          // 600 hill is length=1, dropped
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
}
END_SECTION

START_SECTION([EXTRA] min_hill_length=1 keeps single-scan peaks (sparse-data mode))
{
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}}},
    {1.05, {{500.0, 400.0}}},
    {1.10, {{500.0, 800.0}, {600.0, 999.0}}}, // 600 only here
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  Params p = defaultParams();
  p.min_hill_length = 1;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p);

  TEST_EQUAL(result.size(), 2);
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
  TEST_REAL_SIMILAR(result[1].mz, 600.0);
}
END_SECTION

START_SECTION([EXTRA] max_scan_gap=0 splits a hill broken by one empty scan)
{
  // Same ion at m/z=500 in scans {0,1,3,4} — scan 2 is missing.
  // With strict consecutive linking the hill breaks at the gap and we get
  // two length-2 hills; both survive min_hill_length=2.
  vector<double> mz =        {500.0, 500.0, 500.0, 500.0};
  vector<double> intensity = {100.0, 400.0, 400.0, 100.0};
  vector<double> im =        {1.00,  1.05,  1.15,  1.20 };
  vector<uint32_t> scan_id = {0,     1,     3,     4    };

  Params p = defaultParams();
  p.max_scan_gap = 0;  // default behavior
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p, scan_id.data());

  TEST_EQUAL(result.size(), 2);  // hill split by the empty scan
}
END_SECTION

START_SECTION([EXTRA] max_scan_gap=1 bridges a single empty IM scan)
{
  // Same input as above; with gap-1 tolerance the linker should keep the
  // scan-3 peak chained to the scan-1 tail (skipping the empty scan 2).
  vector<double> mz =        {500.0, 500.0, 500.0, 500.0};
  vector<double> intensity = {100.0, 400.0, 400.0, 100.0};
  vector<double> im =        {1.00,  1.05,  1.15,  1.20 };
  vector<uint32_t> scan_id = {0,     1,     3,     4    };

  Params p = defaultParams();
  p.max_scan_gap = 1;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p, scan_id.data());

  TEST_EQUAL(result.size(), 1);
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
  // All four observed scans are part of the hill; the gap doesn't count
  // toward length.
  TEST_EQUAL(result[0].length, 4);
}
END_SECTION

START_SECTION([EXTRA] max_scan_gap=1 does not let a length-1 hill bridge a gap)
{
  // Single observation at scan 0 then empty at scan 1, another single
  // observation at scan 2. Under the "consecutive link before any gap"
  // rule, the length-1 tail from scan 0 must be pruned when scan 1 has no
  // peak at its m/z. Then scan 2's peak starts a SEPARATE length-1 hill.
  // With min_hill_length=2 both are dropped — preventing a false length-2
  // hill from two unrelated noise blips.
  vector<double> mz =        {500.0, 500.0};
  vector<double> intensity = {100.0, 200.0};
  vector<double> im =        {1.00,  1.10 };
  vector<uint32_t> scan_id = {0,     2    };

  Params p = defaultParams();
  p.max_scan_gap   = 1;
  p.min_hill_length = 2;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p, scan_id.data());

  TEST_EQUAL(result.size(), 0);  // both length-1, gap not bridgeable
}
END_SECTION

START_SECTION([EXTRA] max_scan_gap=1 lets a paired hill bridge a subsequent gap)
{
  // Real-ion pattern: consecutive pair at scans 0,1 (hill reaches length 2),
  // empty scan 2, observation at scan 3. The paired hill is allowed to
  // bridge the gap because it's no longer length-1.
  vector<double> mz =        {500.0, 500.0, 500.0};
  vector<double> intensity = {100.0, 400.0, 200.0};
  vector<double> im =        {1.00,  1.05,  1.15 };
  vector<uint32_t> scan_id = {0,     1,     3    };

  Params p = defaultParams();
  p.max_scan_gap = 1;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p, scan_id.data());

  TEST_EQUAL(result.size(), 1);
  TEST_EQUAL(result[0].length, 3);  // gap doesn't count toward length
}
END_SECTION

START_SECTION([EXTRA] max_scan_gap=1 does not bridge two consecutive empty scans)
{
  // Two empty scans in a row exceed gap-1; the linker should not bridge.
  vector<double> mz =        {500.0, 500.0, 500.0, 500.0};
  vector<double> intensity = {100.0, 400.0, 400.0, 100.0};
  vector<double> im =        {1.00,  1.05,  1.20,  1.25 };
  vector<uint32_t> scan_id = {0,     1,     4,     5    };  // gap of 2 scans

  Params p = defaultParams();
  p.max_scan_gap = 1;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p, scan_id.data());

  TEST_EQUAL(result.size(), 2);  // still split — gap=2 exceeds tolerance
}
END_SECTION

START_SECTION([EXTRA] two close-m/z ions outside tolerance are kept separate)
{
  // m/z 500.0 and 500.1 (~200 ppm apart); with mz_tol_ppm=20 they should
  // not link into one hill.
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 200.0}, {500.1, 150.0}}},
    {1.05, {{500.0, 700.0}, {500.1, 500.0}}},
    {1.10, {{500.0, 200.0}, {500.1, 150.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  TEST_EQUAL(result.size(), 2);
  // Output is sorted by m/z ascending.
  TEST_EQUAL(result[0].mz < result[1].mz, true);
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
  TEST_REAL_SIMILAR(result[1].mz, 500.1);
}
END_SECTION

START_SECTION([EXTRA] determinism: same input produces byte-identical output)
{
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}, {600.0, 50.0}}},
    {1.05, {{500.0, 400.0}, {600.0, 200.0}}},
    {1.10, {{500.0, 800.0}, {600.0, 400.0}}},
    {1.15, {{500.0, 400.0}, {600.0, 200.0}}},
    {1.20, {{500.0, 100.0}, {600.0, 50.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto r1 = centroidFrame(mz.data(), intensity.data(), im.data(),
                          mz.size(), defaultParams());
  auto r2 = centroidFrame(mz.data(), intensity.data(), im.data(),
                          mz.size(), defaultParams());

  TEST_EQUAL(r1.size(), r2.size());
  for (Size i = 0; i < r1.size(); ++i)
  {
    TEST_REAL_SIMILAR(r1[i].mz,        r2[i].mz);
    TEST_REAL_SIMILAR(r1[i].intensity, r2[i].intensity);
    TEST_REAL_SIMILAR(r1[i].im_apex,   r2[i].im_apex);
    TEST_EQUAL(r1[i].length,           r2[i].length);
  }
}
END_SECTION

START_SECTION([EXTRA] min_intensity filter drops low-intensity input peaks)
{
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}, {600.0, 0.5}}}, // 600 below threshold
    {1.05, {{500.0, 400.0}, {600.0, 0.5}}},
    {1.10, {{500.0, 800.0}, {600.0, 0.5}}},
    {1.15, {{500.0, 400.0}}},
    {1.20, {{500.0, 100.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  Params p = defaultParams();
  p.min_intensity = 1.0;
  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), p);

  TEST_EQUAL(result.size(), 1);
  TEST_REAL_SIMILAR(result[0].mz, 500.0);
}
END_SECTION

START_SECTION([EXTRA] no valley split when intensities form a smooth single peak)
{
  // Strictly increasing then decreasing — should NOT split.
  vector<ScanPeaks> scans = {
    {1.00, {{500.0, 100.0}}},
    {1.05, {{500.0, 300.0}}},
    {1.10, {{500.0, 600.0}}},
    {1.15, {{500.0, 1000.0}}},
    {1.20, {{500.0, 600.0}}},
    {1.25, {{500.0, 300.0}}},
    {1.30, {{500.0, 100.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  TEST_EQUAL(result.size(), 1);
  TEST_EQUAL(result[0].length, 7);
}
END_SECTION

START_SECTION([EXTRA] output is sorted by m/z ascending)
{
  // Build input with peaks intentionally not in m/z order across IM scans.
  vector<ScanPeaks> scans = {
    {1.00, {{700.0, 50.0}, {500.0, 100.0}, {600.0, 80.0}}},
    {1.05, {{700.0, 50.0}, {500.0, 400.0}, {600.0, 80.0}}},
    {1.10, {{700.0, 50.0}, {500.0, 800.0}, {600.0, 80.0}}},
  };
  vector<double> mz, intensity, im;
  buildFlat(scans, mz, intensity, im);

  auto result = centroidFrame(mz.data(), intensity.data(), im.data(),
                               mz.size(), defaultParams());

  TEST_EQUAL(result.size(), 3);
  for (Size i = 1; i < result.size(); ++i)
    TEST_EQUAL(result[i].mz >= result[i-1].mz, true);
}
END_SECTION

END_TEST
