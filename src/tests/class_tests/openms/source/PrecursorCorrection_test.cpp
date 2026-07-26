// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka$
// $Authors: Oliver Alka$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <cmath>
#include <limits>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PrecursorCorrection, "$Id$")

/////////////////////////////////////////////////////////////

// Prepare dummy data
MSExperiment exp;
vector<Precursor> v_precursor_1, v_precursor_2, v_precursor_3;
Precursor precursor_1, precursor_2, precursor_3;
MSSpectrum ms1_spectrum_1, ms1_spectrum_2, ms1_spectrum_3, ms2_spectrum_1, ms2_spectrum_2, ms2_spectrum_3;
vector<MSSpectrum> v_spectra;

// precursor
precursor_1.setIntensity(240.0f);
precursor_1.setMZ(509.9999);
v_precursor_1.push_back(precursor_1);

precursor_2.setIntensity(230.0f);
precursor_2.setMZ(610.0001);
precursor_2.setCharge(1);
v_precursor_2.push_back(precursor_2);

precursor_3.setIntensity(220.0f);
precursor_3.setMZ(611.0035);
precursor_3.setCharge(1);
v_precursor_3.push_back(precursor_3);

// peaks
Peak1D p1;
p1.setIntensity(200.0f);
p1.setMZ(509.9994);

Peak1D p2;
p2.setIntensity(250.0f);
p2.setMZ(510.0000);

Peak1D p3;
p3.setIntensity(150.0f);
p3.setMZ(510.0001);

Peak1D p4;
p4.setIntensity(250.0f);
p4.setMZ(609.9998);

Peak1D p5;
p5.setIntensity(200.0f);
p5.setMZ(610.0000);

Peak1D p6;
p6.setIntensity(180.0f);
p6.setMZ(610.0005);

Peak1D p7;
p7.setIntensity(250.0f);
p7.setMZ(611.0031);

Peak1D p8;
p8.setIntensity(200.0f);
p8.setMZ(611.0033);

Peak1D p9;
p9.setIntensity(180.0f);
p9.setMZ(611.0038);

vector<Peak1D> peaks_1{p1,p2,p3};
vector<Peak1D> peaks_2{p4,p5,p6};
vector<Peak1D> peaks_3{p7,p8,p9};
vector<Peak1D> empty_peaks{};

// ms1
ms1_spectrum_1.insert(ms1_spectrum_1.begin(), peaks_1.begin(), peaks_1.end());
ms1_spectrum_1.setMSLevel(1);
ms1_spectrum_1.setNativeID("scan=1");
ms1_spectrum_1.setRT(100.0);
ms1_spectrum_2.insert(ms1_spectrum_2.begin(), peaks_2.begin(), peaks_2.end());
ms1_spectrum_2.setMSLevel(1);
ms1_spectrum_2.setNativeID("scan=3");
ms1_spectrum_2.setRT(180.85);
ms1_spectrum_3.insert(ms1_spectrum_3.begin(), peaks_3.begin(), peaks_3.end());
ms1_spectrum_2.setNativeID("scan=5");
ms1_spectrum_3.setMSLevel(1);
ms1_spectrum_3.setRT(183.85);

// ms2
ms2_spectrum_1.insert(ms2_spectrum_1.begin(), empty_peaks.begin(), empty_peaks.end());
ms2_spectrum_1.setMSLevel(2);
ms2_spectrum_1.setNativeID("scan=2");
ms2_spectrum_1.setRT(100.1);
ms2_spectrum_2.insert(ms2_spectrum_2.begin(), empty_peaks.begin(), empty_peaks.end());
ms2_spectrum_2.setMSLevel(2);
ms2_spectrum_2.setNativeID("scan=4");
ms2_spectrum_2.setRT(180.90);
ms2_spectrum_3.insert(ms2_spectrum_3.begin(), empty_peaks.begin(), empty_peaks.end());
ms2_spectrum_3.setMSLevel(2);
ms2_spectrum_3.setNativeID("scan=6");
ms2_spectrum_3.setRT(183.92);

// ms2 precursor information
ms2_spectrum_1.setPrecursors(v_precursor_1);
ms2_spectrum_2.setPrecursors(v_precursor_2);
ms2_spectrum_3.setPrecursors(v_precursor_3);

v_spectra.push_back(ms1_spectrum_1);
v_spectra.push_back(ms2_spectrum_1);
v_spectra.push_back(ms1_spectrum_2);
v_spectra.push_back(ms2_spectrum_2);
v_spectra.push_back(ms1_spectrum_3);
v_spectra.push_back(ms2_spectrum_3);

// MSExperiment
exp.setSpectra(v_spectra);
exp.sortSpectra();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PrecursorCorrection* ptr = nullptr;
PrecursorCorrection* null_ptr = nullptr;
START_SECTION(PrecursorCorrection())
{
	ptr = new PrecursorCorrection();
	TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~PrecursorCorrection())
{
	delete ptr;
}
END_SECTION

START_SECTION((static void getPrecursors(const MSExperiment &exp, std::vector< Precursor > &precursors, std::vector< double > &precursors_rt, std::vector< Size > &precursor_scan_index)))
{
  MSExperiment getP_exp = exp;
  vector<Precursor> precursor;
  vector<double> rt;
  vector<Size> index;
  PrecursorCorrection::getPrecursors(getP_exp, precursor, rt, index);

  TEST_EQUAL(precursor.size(), 3);
  TEST_EQUAL(rt.size(), 3);
  TEST_EQUAL(index.size(), 3);
  TEST_REAL_SIMILAR(precursor[0].getMZ(), 509.9999)
  TEST_REAL_SIMILAR(precursor[0].getIntensity(), 240.0);
  TEST_REAL_SIMILAR(rt[0], 100.1);
}
END_SECTION

FuzzyStringComparator fsc;
fsc.setAcceptableAbsolute(1e-8);

START_SECTION((static void writeHist(const std::string &out_csv, const std::vector< double > &deltaMZs, const std::vector< double > &mzs, const std::vector< double > &rts)))
{
  MSExperiment write_exp = exp;

  std::string csv_tmp;
  NEW_TMP_FILE(csv_tmp);
  vector<double> dmz;
  vector<double> mz;
  vector<double> rt;

  PrecursorCorrection::correctToNearestMS1Peak(write_exp, 2, true, dmz, mz, rt);
  PrecursorCorrection::writeHist(csv_tmp, dmz, mz,rt);

  TEST_EQUAL(fsc.compareFiles(csv_tmp, OPENMS_GET_TEST_DATA_PATH("PrecursorCorrection_out.csv")),true);
}
END_SECTION

START_SECTION((static std::set<Size> correctToNearestMS1Peak(MSExperiment &exp, double mz_tolerance, bool ppm, std::vector< double > &deltaMZs, std::vector< double > &mzs, std::vector< double > &rts)))
{
  // test with 1 ppm (1)
  MSExperiment nearest_exp_1 = exp;
  vector<double> dmz_1;
  vector<double> mz_1;
  vector<double> rt_1;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  // corrected precursor_3: none
  PrecursorCorrection::correctToNearestMS1Peak(nearest_exp_1, 1, true, dmz_1, mz_1, rt_1);

  TEST_REAL_SIMILAR(dmz_1[0], 0.0001);
  TEST_REAL_SIMILAR(dmz_1[1], -0.0001);

  // test with 5 ppm (2)
  MSExperiment nearest_exp_2 = exp;
  vector<double> dmz_2;
  vector<double> mz_2;
  vector<double> rt_2;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  // corrected precursor_3: 611.0033
  PrecursorCorrection::correctToNearestMS1Peak(nearest_exp_2, 5, true, dmz_2, mz_2, rt_2);

  TEST_REAL_SIMILAR(dmz_2[0], 0.0001);
  TEST_REAL_SIMILAR(dmz_2[1], -0.0001);
  TEST_REAL_SIMILAR(dmz_2[2], -0.0002);
}
END_SECTION

START_SECTION((static std::set<Size> correctToHighestIntensityMS1Peak(MSExperiment &exp, double mz_tolerance, bool ppm, std::vector< double > &deltaMZs, std::vector< double > &mzs, std::vector< double > &rts)))
{
  // test with 0.0001 Da (1)
  MSExperiment highest_exp_1 = exp;
  vector<double> dmz_1;
  vector<double> mz_1;
  vector<double> rt_1;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  // corrected precursor_3: none
  PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_1, 0.0001, false, dmz_1, mz_1, rt_1);

  TEST_REAL_SIMILAR(dmz_1[0], 0.0001);
  TEST_REAL_SIMILAR(dmz_1[1], -0.0001);

  // test with 0.0005 Da (2)
  MSExperiment highest_exp_2 = exp;
  vector<double> dmz_2;
  vector<double> mz_2;
  vector<double> rt_2;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 609.9998
  // corrected precursor_3: 611.0031
  PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_2, 0.0005, false, dmz_2, mz_2, rt_2);

  TEST_REAL_SIMILAR(dmz_2[0], 0.0001);
  TEST_REAL_SIMILAR(dmz_2[1], -0.0003);
  TEST_REAL_SIMILAR(dmz_2[2], -0.0004);
}
END_SECTION

// check ppm
START_SECTION((static std::set<Size> correctToHighestIntensityMS1Peak(MSExperiment &exp, double mz_tolerance, bool ppm, std::vector< double > &deltaMZs, std::vector< double > &mzs, std::vector< double > &rts)))
{
// test with 1 ppm (1)
MSExperiment highest_exp_1 = exp;
vector<double> dmz_1;
vector<double> mz_1;
vector<double> rt_1;

// corrected precursor_1: 510.0000
// corrected precursor_2: 609.9998 (1 ppm of 610.0001 = +/- 0.000610)
// corrected precursor_3: none
PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_1, 2, true, dmz_1, mz_1, rt_1);


std::cout << dmz_1[0] << std::endl;
std::cout << dmz_1[1] << std::endl;

TEST_REAL_SIMILAR(dmz_1[0], 0.0001);
TEST_REAL_SIMILAR(dmz_1[1], -0.0003);

// test with 5 ppm Da
MSExperiment highest_exp_2 = exp;
vector<double> dmz_2;
vector<double> mz_2;
vector<double> rt_2;

// corrected precursor_1: 510.0000
// corrected precursor_2: 609.9998
// corrected precursor_3: 611.0031
PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_2, 5, true, dmz_2, mz_2, rt_2);

std::cout << dmz_2[0] << std::endl;
std::cout << dmz_2[1] << std::endl;
std::cout << dmz_2[2] << std::endl;

TEST_REAL_SIMILAR(dmz_2[0], 0.0001);
TEST_REAL_SIMILAR(dmz_2[1], -0.0003);
TEST_REAL_SIMILAR(dmz_2[2], -0.0004);
}
END_SECTION

// FeatureMap
DPosition<2> position_1(175.0, 609.9100);
DPosition<2> position_2(185.0, 611.9300);
vector<DPosition<2> > vec;
vec.push_back(position_1);
vec.push_back(position_2);

ConvexHull2D hull;
hull.setHullPoints(vec);
hull.expandToBoundingBox();
vector<ConvexHull2D> hulls;
hulls.push_back(hull);

FeatureMap fmap;
Feature feature;
feature.setMZ(610.0000);
feature.setRT(180.0);
feature.setCharge(1);
feature.setConvexHulls(hulls);
fmap.push_back(feature);

START_SECTION((static std::set<Size> correctToNearestFeature(const FeatureMap &features, MSExperiment &exp, double rt_tolerance_s=0.0, double mz_tolerance=0.0, bool ppm=true, bool believe_charge=false, bool keep_original=false, bool all_matching_features=false, int max_trace=2, int debug_level=0)))
{
  MSExperiment f_exp = exp;
  double rt_tolerance = 5;
  double mz_tolerance = 5;
  bool ppm = true;

  vector<Precursor> precursor_before_correction;
  vector<Precursor> precursor_after_correction;

  vector<MSSpectrum> f_spectra_before = f_exp.getSpectra();
  for (const auto& it : f_spectra_before)
  {
    if (it.getNativeID() == "scan=6")
    {
      precursor_before_correction = it.getPrecursors();
    }
  }

  // the precursor of the ms2 with nativeID 6 should be corrected
  PrecursorCorrection::correctToNearestFeature(fmap, f_exp, rt_tolerance, mz_tolerance, ppm);

  vector<MSSpectrum> f_spectra_after = f_exp.getSpectra();
  for (const auto& it : f_spectra_after)
  {
    if (it.getNativeID() == "scan=6")
    {
      precursor_after_correction = it.getPrecursors();
    }
  }

  TEST_EQUAL(precursor_before_correction.size(), 1);
  TEST_EQUAL(precursor_after_correction.size(), 1);
  TEST_REAL_SIMILAR(precursor_before_correction[0].getPos(), 611.0035);
  TEST_REAL_SIMILAR(precursor_after_correction[0].getPos(), 610.0000);

  // Matching must use the feature's convex-hull bounding box, not its centroid (issue #4787).
  // The MS2 below sits on the lower edge of the hull, 9.07 s away from the feature centroid,
  // while rt_tolerance is 0 - so a centroid-based lookup would find nothing here.
  // Values mirror the TOPP fixture HighResPrecursorMassCorrector_2538_1091_2.
  {
    vector<DPosition<2> > hull_vec;
    hull_vec.push_back(DPosition<2>(2529.7605, 1091.53929624078));
    hull_vec.push_back(DPosition<2>(2572.4499, 1091.54503004844));
    ConvexHull2D wide_hull;
    wide_hull.setHullPoints(hull_vec);
    wide_hull.expandToBoundingBox();
    vector<ConvexHull2D> wide_hulls;
    wide_hulls.push_back(wide_hull);

    Feature wide_feature;
    wide_feature.setMZ(1091.54301896188);
    wide_feature.setRT(2538.83463437242); // centroid, 9.07 s away from the MS2 below
    wide_feature.setCharge(2);
    wide_feature.setConvexHulls(wide_hulls);

    FeatureMap wide_fmap;
    wide_fmap.push_back(wide_feature);

    Precursor edge_precursor;
    edge_precursor.setMZ(1091.5400);
    edge_precursor.setCharge(2);
    vector<Precursor> edge_precursors;
    edge_precursors.push_back(edge_precursor);

    MSSpectrum edge_ms2;
    edge_ms2.setMSLevel(2);
    edge_ms2.setRT(2529.7605); // on the hull edge, nowhere near the centroid
    edge_ms2.setPrecursors(edge_precursors);

    MSExperiment edge_exp;
    edge_exp.addSpectrum(edge_ms2);

    // rt_tolerance = 0: only the hull extent can bring this MS2 and this feature together
    set<Size> corrected = PrecursorCorrection::correctToNearestFeature(wide_fmap, edge_exp, 0.0, 10.0, true);

    TEST_EQUAL(corrected.size(), 1);
    TEST_REAL_SIMILAR(edge_exp[0].getPrecursors()[0].getMZ(), 1091.54301896188);
  }

  // Exercises the m/z sweep used to find candidate features (issue #4787). Three features share
  // one RT range but differ in m/z: two narrow ones at 500 and 700, plus a wide one spanning
  // 500-700 whose charge (3) does not match the first and last precursor. The MS2 spectra are
  // stored out of m/z order on purpose.
  //   - the spectra must be matched independently of their storage order
  //   - the narrow 500 feature must not still be considered once m/z has moved past it
  //   - the wide feature must survive being skipped by believe_charge at 500 and still match at 600
  //   - a non-finite precursor m/z must not disturb any of the above
  {
    FeatureMap sweep_fmap;
    const double f_mz[3]    = {500.0, 700.0, 600.0};
    const double f_lo[3]    = {499.995, 699.995, 500.0};
    const double f_hi[3]    = {500.005, 700.005, 700.0};
    const int    f_charge[3] = {2, 2, 3};
    for (Size f = 0; f < 3; ++f)
    {
      vector<DPosition<2> > pts;
      pts.push_back(DPosition<2>(100.0, f_lo[f]));
      pts.push_back(DPosition<2>(200.0, f_hi[f]));
      ConvexHull2D h;
      h.setHullPoints(pts);
      h.expandToBoundingBox();
      vector<ConvexHull2D> hs;
      hs.push_back(h);

      Feature feat;
      feat.setRT(150.0);
      feat.setMZ(f_mz[f]);
      feat.setCharge(f_charge[f]);
      feat.setConvexHulls(hs);
      sweep_fmap.push_back(feat);
    }

    // stored order is 700, 500, 600, NaN - i.e. deliberately not ascending in m/z
    const double pc_mz[4]    = {700.0020, 500.0020, 600.0020, std::numeric_limits<double>::quiet_NaN()};
    const int    pc_charge[4] = {2, 2, 3, 2};
    MSExperiment sweep_exp;
    for (Size s = 0; s < 4; ++s)
    {
      Precursor pc;
      pc.setMZ(pc_mz[s]);
      pc.setCharge(pc_charge[s]);
      vector<Precursor> pcs;
      pcs.push_back(pc);

      MSSpectrum ms2;
      ms2.setMSLevel(2);
      ms2.setRT(150.0);
      ms2.setPrecursors(pcs);
      sweep_exp.addSpectrum(ms2);
    }

    set<Size> sweep_corrected = PrecursorCorrection::correctToNearestFeature(
        sweep_fmap, sweep_exp, 0.0, 10.0, true, /* believe_charge */ true);

    TEST_EQUAL(sweep_corrected.size(), 3);
    TEST_REAL_SIMILAR(sweep_exp[0].getPrecursors()[0].getMZ(), 700.0); // narrow, charge 2
    TEST_REAL_SIMILAR(sweep_exp[1].getPrecursors()[0].getMZ(), 500.0); // narrow, charge 2
    TEST_REAL_SIMILAR(sweep_exp[2].getPrecursors()[0].getMZ(), 600.0); // wide, charge 3
    TEST_EQUAL(sweep_exp[2].getPrecursors()[0].getCharge(), 3);
    // the non-finite precursor is left untouched
    TEST_TRUE(std::isnan(sweep_exp[3].getPrecursors()[0].getMZ()));
  }

  // Covers the other fallback: a feature that cannot be placed on the sweep axis at all.
  // DBoundingBox only ever absorbs a coordinate through '<' / '>' comparisons, which a NaN never
  // satisfies, so a NaN hull point cannot reach the box - an infinite one is the only way to end up
  // with a non-finite m/z interval. Such a feature is excluded from the sorted sweep and has to be
  // compared against every precursor instead; if that fallback is dropped, the match below is lost.
  {
    vector<DPosition<2> > inf_pts;
    inf_pts.push_back(DPosition<2>(100.0, 600.0));
    inf_pts.push_back(DPosition<2>(200.0, std::numeric_limits<double>::infinity()));
    ConvexHull2D inf_hull;
    inf_hull.setHullPoints(inf_pts);
    inf_hull.expandToBoundingBox();
    vector<ConvexHull2D> inf_hulls;
    inf_hulls.push_back(inf_hull);

    Feature inf_feature;
    inf_feature.setRT(150.0);
    inf_feature.setMZ(600.0);
    inf_feature.setCharge(2);
    inf_feature.setConvexHulls(inf_hulls);

    // the m/z interval runs to +inf, so this feature must not be sorted onto the sweep axis
    TEST_FALSE(std::isfinite(inf_feature.getConvexHull().getBoundingBox().maxPosition()[1]));

    FeatureMap inf_fmap;
    inf_fmap.push_back(inf_feature);

    Precursor inf_pc;
    inf_pc.setMZ(600.0020);
    inf_pc.setCharge(2);
    vector<Precursor> inf_pcs;
    inf_pcs.push_back(inf_pc);

    MSSpectrum inf_ms2;
    inf_ms2.setMSLevel(2);
    inf_ms2.setRT(150.0);
    inf_ms2.setPrecursors(inf_pcs);

    MSExperiment inf_exp;
    inf_exp.addSpectrum(inf_ms2);

    set<Size> inf_corrected = PrecursorCorrection::correctToNearestFeature(inf_fmap, inf_exp, 0.0, 10.0, true);

    TEST_EQUAL(inf_corrected.size(), 1);
    TEST_REAL_SIMILAR(inf_exp[0].getPrecursors()[0].getMZ(), 600.0);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



