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

START_SECTION((static void writeHist(const String &out_csv, const std::vector< double > &deltaMZs, const std::vector< double > &mzs, const std::vector< double > &rts)))
{
  MSExperiment write_exp = exp;

  String csv_tmp;
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
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



