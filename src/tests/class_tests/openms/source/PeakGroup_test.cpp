// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FLASHHelperClasses::LogMzPeak LogMzPeak;
LogMzPeak fillPeak(double mz, float it, int cs, int iso_idx)
{
  Peak1D p;
  p.setIntensity(it);
  p.setMZ(mz);
  LogMzPeak lmp(p, true);
  lmp.abs_charge = cs;
  lmp.isotopeIndex = iso_idx;
  return lmp;
}

START_TEST(PeakGroup, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakGroup* ptr = 0;
PeakGroup* null_ptr = 0;
START_SECTION(PeakGroup())
{
  ptr = new PeakGroup();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeakGroup())
{
  delete ptr;
}
END_SECTION

/// test data
PeakGroup sample_pg(1, 2, true);
sample_pg.setScanNumber(3);

LogMzPeak tmp_peak0 = fillPeak(1125.5118055019082, 443505.625, 2, 0);
sample_pg.push_back(tmp_peak0);

LogMzPeak tmp_peak1 = fillPeak(1126.0134829208082, 11212854, 2, 1);
sample_pg.push_back(tmp_peak1);

LogMzPeak tmp_peak2 = fillPeak(1126.515160339708, 1214510.5, 2, 2);
sample_pg.push_back(tmp_peak2);

LogMzPeak tmp_peak3 = fillPeak(1127.0168377586081, 7506.6767578125, 2, 3);
sample_pg.push_back(tmp_peak3);
sample_pg.updateMonoMassAndIsotopeIntensities(1e-5);

/// detailed constructor test
START_SECTION((PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive)))
{
  PeakGroup tmp_pg(1, 2, true);
  TEST_EQUAL(std::get<0>(tmp_pg.getAbsChargeRange()), std::get<0>(sample_pg.getAbsChargeRange()));
  TEST_EQUAL(std::get<1>(tmp_pg.getAbsChargeRange()), std::get<1>(sample_pg.getAbsChargeRange()));
  TEST_EQUAL(tmp_pg.isPositive(), tmp_pg.isPositive());
}
END_SECTION

/// copy constructor test
START_SECTION((PeakGroup(const PeakGroup &)))
{
  PeakGroup copy_pg(sample_pg);

  TEST_EQUAL(std::get<0>(sample_pg.getAbsChargeRange()), std::get<0>(copy_pg.getAbsChargeRange()));
  TEST_EQUAL(sample_pg.size(), copy_pg.size());
  TEST_REAL_SIMILAR(sample_pg[0].intensity, copy_pg[0].intensity);
  TEST_REAL_SIMILAR(sample_pg[1].mz, copy_pg[1].mz);
}
END_SECTION

/// assignment constructor test
START_SECTION((PeakGroup& operator=(const PeakGroup &t)))
{
  PeakGroup tmp_pg = sample_pg;

  TEST_EQUAL(std::get<0>(sample_pg.getAbsChargeRange()), std::get<0>(tmp_pg.getAbsChargeRange()));
  TEST_EQUAL(sample_pg.size(), tmp_pg.size());
  TEST_REAL_SIMILAR(sample_pg[0].intensity, tmp_pg[0].intensity);
  TEST_REAL_SIMILAR(sample_pg[1].mz, tmp_pg[1].mz);
}
END_SECTION


/////////////////////////////////////////////////////////////
// accessor method tests
/////////////////////////////////////////////////////////////

START_SECTION((std::tuple<double, double> getMzRange(int abs_charge) const))
{
  std::tuple<double, double> temp_range = sample_pg.getMzRange(2);
  TEST_REAL_SIMILAR(std::get<0>(temp_range), 1125.5118055019082);
  TEST_REAL_SIMILAR(std::get<1>(temp_range), 1127.0168377586081);
}
END_SECTION

START_SECTION((bool isPositive() const))
{
  bool test_positive = sample_pg.isPositive();
  TEST_EQUAL(test_positive, true);
}
END_SECTION

START_SECTION((int getScanNumber() const))
{
  int test_scan_num = sample_pg.getScanNumber();
  TEST_EQUAL(test_scan_num, 3);
}
END_SECTION

START_SECTION((void setScanNumber(const int scan_number)))
{
  sample_pg.setScanNumber(5);
  int test_scan_num = sample_pg.getScanNumber();
  TEST_EQUAL(test_scan_num, 5);
}
END_SECTION

/// not testable : setChargePower_, setChargeSignalPower_ - no getter for private variable

START_SECTION((void setChargeIsotopeCosine(const int abs_charge, const float cos)))
{
  sample_pg.setChargeIsotopeCosine(2, 0.4);
  TEST_REAL_SIMILAR(sample_pg.getChargeIsotopeCosine(2), 0.4);
}
END_SECTION

START_SECTION((float getChargeIsotopeCosine(const int abs_charge) const))
{
  TEST_REAL_SIMILAR(sample_pg.getChargeIsotopeCosine(0), .0);
  TEST_REAL_SIMILAR(sample_pg.getChargeIsotopeCosine(2), 0.4);
}
END_SECTION


START_SECTION((float getChargeIntensity(const int abs_charge) const))
{
  TEST_REAL_SIMILAR(sample_pg.getChargeIntensity(2), .0);
}
END_SECTION

START_SECTION((void setRepAbsCharge(const int max_qscore_charge)))
{
  sample_pg.setRepAbsCharge(2);
  int temp_abs = sample_pg.getRepAbsCharge();
  TEST_EQUAL(temp_abs, 2);
}
END_SECTION


START_SECTION((std::tuple<double, double> getRepMzRange() const))
{
  std::tuple<double, double> tmp_range = sample_pg.getRepMzRange();
  TEST_REAL_SIMILAR(std::get<0>(tmp_range), 1125.5118055019082);
  TEST_REAL_SIMILAR(std::get<1>(tmp_range), 1127.0168377586081);
}
END_SECTION


START_SECTION((std::tuple<int, int> getAbsChargeRange() const))
{
  std::tuple<int, int> test_cs_range = sample_pg.getAbsChargeRange();
  TEST_EQUAL(std::get<0>(test_cs_range), 1);
  TEST_EQUAL(std::get<1>(test_cs_range), 2);
}
END_SECTION

START_SECTION((void setAbsChargeRange(const int min_abs_charge, const int max_abs_charge)))
{
  PeakGroup sample_pg2(4, 9, true); // for operator test
  std::tuple<int, int> test_cs_range = sample_pg2.getAbsChargeRange();
  TEST_EQUAL(std::get<0>(test_cs_range), 4);
  TEST_EQUAL(std::get<1>(test_cs_range), 9);
}
END_SECTION


START_SECTION((void setIsotopeCosine(const float cos)))
{
  sample_pg.setIsotopeCosine(0.3);
  double temp_iso = sample_pg.getIsotopeCosine();
  TEST_REAL_SIMILAR(temp_iso, 0.3);
}
END_SECTION

START_SECTION((float getIsotopeCosine() const))
{
  double temp_iso = sample_pg.getIsotopeCosine();
  TEST_REAL_SIMILAR(temp_iso, 0.3);
}
END_SECTION


START_SECTION((int getRepAbsCharge() const))
{
  int temp_abs = sample_pg.getRepAbsCharge();
  TEST_EQUAL(temp_abs, 2);
}
END_SECTION


START_SECTION((void setQscore(const float qscore)))
{
  sample_pg.setQscore(0.1);
  double temp_score = sample_pg.getQscore();
  TEST_REAL_SIMILAR(temp_score, 0.1);
}
END_SECTION

START_SECTION((float getQscore() const))
{
  double temp_score = sample_pg.getQscore();
  TEST_REAL_SIMILAR(temp_score, 0.1);
}
END_SECTION


START_SECTION((void setChargeScore(const float charge_score)))
{
  sample_pg.setChargeScore(0.2);
  double temp_score = sample_pg.getChargeScore();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION

START_SECTION((float getChargeScore() const))
{
  double temp_score = sample_pg.getChargeScore();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION


START_SECTION((void setAvgPPMError(const float error)))
{
  sample_pg.setAvgPPMError(0.2);
  double temp_score = sample_pg.getAvgPPMError();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION

START_SECTION((float getAvgPPMError() const))
{
  double temp_score = sample_pg.getAvgPPMError();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION


START_SECTION((void setSNR(const float snr)))
{
  sample_pg.setSNR(0.2);
  double temp_score = sample_pg.getSNR();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION

START_SECTION((float getSNR() const))
{
  double temp_score = sample_pg.getSNR();
  TEST_REAL_SIMILAR(temp_score, 0.2);
}
END_SECTION


START_SECTION((void setChargeSNR(const int abs_charge, const float c_snr)))
{
  sample_pg.setChargeSNR(2, 0.2);
  TEST_REAL_SIMILAR(sample_pg.getChargeSNR(2), 0.2);
}
END_SECTION

START_SECTION((float getChargeSNR(const int abs_charge) const))
{
  TEST_REAL_SIMILAR(sample_pg.getChargeSNR(0), .0);
  TEST_REAL_SIMILAR(sample_pg.getChargeSNR(2), 0.2);
}
END_SECTION

START_SECTION((double getMonoMass() const))
{
  double tmp_mass = sample_pg.getMonoMass();
  TEST_REAL_SIMILAR(tmp_mass, 2249.0101019557173);
}
END_SECTION

START_SECTION((double getIntensity() const))
{
  double tmp_inty = sample_pg.getIntensity();
  TEST_REAL_SIMILAR(tmp_inty, 12878380.801757813)
}
END_SECTION


PeakGroup sample_pg2(sample_pg);
LogMzPeak tmp_peak4 = fillPeak(1127.5185151766082, 2504.3433, 2, 4);
sample_pg2.push_back(tmp_peak4);
sample_pg2.updateMonoMassAndIsotopeIntensities(1e-5);
START_SECTION((void updateMonom assAndIsotopeIntensities()))
{
  double temp_mass = sample_pg2.getMonoMass();
  double temp_inty = sample_pg2.getIntensity();
  TEST_REAL_SIMILAR(temp_mass, 2249.0101025181098);
  TEST_REAL_SIMILAR(temp_inty, 12880881);
}
END_SECTION


/// operator constructor test
START_SECTION((bool operator<(const PeakGroup &a) const))
{
  // between two masses with different monoisotopic masses
  bool is_pg2_bigger = sample_pg < sample_pg2;
  TEST_EQUAL(is_pg2_bigger, true);
}
END_SECTION

START_SECTION((bool operator>(const PeakGroup &a) const))
{
  // between two masses with different monoisotopic masses
  bool is_pg2_bigger = sample_pg2 > sample_pg;
  TEST_EQUAL(is_pg2_bigger, true);
}
END_SECTION

START_SECTION((bool operator==(const PeakGroup &a) const))
{
  PeakGroup sample_pg4(sample_pg);

  bool are_two_pgs_same = sample_pg == sample_pg4;
  TEST_EQUAL(are_two_pgs_same, true);
}
END_SECTION


/// TODOs
/// - updateIsotopeCosineAndQscore, recruitAllPeaksInSpectrum, isSignalMZ, setTargeted, getIsotopeIntensities
/// - isTargeted, getTargetDecoyType, setTargetDecoyType, getQvalue, setQvalue, getQvalueWithChargeDecoyOnly, setQvalueWithChargeDecoyOnly


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
