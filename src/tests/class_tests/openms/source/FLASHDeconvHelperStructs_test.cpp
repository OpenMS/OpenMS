// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FLASHHelperClasses::LogMzPeak LogMzPeak;
typedef FLASHHelperClasses::PrecalculatedAveragine PrecalculatedAveragine;
typedef FLASHHelperClasses::MassFeature MassFeature;
typedef FLASHHelperClasses::IsobaricQuantities IsobaricQuantities;

START_TEST(FLASHDeconvHelperStructs, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHHelperClasses* ptr = 0;
FLASHHelperClasses* null_ptr = 0;
START_SECTION(FLASHHelperClasses())
{
  ptr = new FLASHHelperClasses();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHHelperClasses())
{
  delete ptr;
}
END_SECTION


START_SECTION((static double getLogMz(const double mz, const bool positive)))
{
  double mz = 1300;
  double tmp_lmz1 = OpenMS::FLASHHelperClasses::getLogMz(mz, true);
  double tmp_lmz2 = OpenMS::FLASHHelperClasses::getLogMz(mz, false);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_lmz1, 7.169344415063863);
  TEST_REAL_SIMILAR(tmp_lmz2, 7.170119121465);
}
END_SECTION

START_SECTION((static double getChargeMass(const bool positive)))
{
  double temp_pos = OpenMS::FLASHHelperClasses::getChargeMass(true);
  double temp_neg = OpenMS::FLASHHelperClasses::getChargeMass(false);
  TEST_REAL_SIMILAR(temp_pos, Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(temp_neg, -Constants::PROTON_MASS_U);
}
END_SECTION


/// testing LogMzPeak
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak()=default))
{
  LogMzPeak* lmp_ptr = new LogMzPeak();
  LogMzPeak* lmp_null_ptr = 0;

  TEST_NOT_EQUAL(lmp_ptr, lmp_null_ptr);
  delete lmp_ptr;
}
END_SECTION

// test data
Peak1D tmp_p1;
tmp_p1.setIntensity(443505.625);
tmp_p1.setMZ(1125.5118055019082);
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak(const Peak1D &peak, const bool positive)))
{
  // Test positive mode
  LogMzPeak tmp_peak(tmp_p1, true);
  TEST_REAL_SIMILAR(tmp_peak.mz, 1125.5118055019082);
  TEST_REAL_SIMILAR(tmp_peak.intensity, 443505.625);
  TEST_REAL_SIMILAR(tmp_peak.logMz, 7.0250977989903145);
  TEST_EQUAL(tmp_peak.is_positive, true);
  TEST_EQUAL(tmp_peak.abs_charge, 0);
  TEST_EQUAL(tmp_peak.isotopeIndex, 0);

  // Test negative mode
  LogMzPeak tmp_peak_neg(tmp_p1, false);
  TEST_REAL_SIMILAR(tmp_peak_neg.mz, 1125.5118055019082);
  TEST_REAL_SIMILAR(tmp_peak_neg.intensity, 443505.625);
  TEST_EQUAL(tmp_peak_neg.is_positive, false);
  // logMz should be different for negative mode due to charge mass difference
  TEST_NOT_EQUAL(tmp_peak_neg.logMz, tmp_peak.logMz);
}
END_SECTION

LogMzPeak test_peak(tmp_p1, true);
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] LogMzPeak(const LogMzPeak &)))
{
  LogMzPeak tmp_p(test_peak);
  TEST_REAL_SIMILAR(test_peak.mz, tmp_p.mz);
  TEST_REAL_SIMILAR(test_peak.intensity, tmp_p.intensity);
  TEST_REAL_SIMILAR(test_peak.logMz, tmp_p.logMz);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] double getUnchargedMass()))
{
  // Test with charge set
  test_peak.abs_charge = 2;
  TEST_REAL_SIMILAR(test_peak.getUnchargedMass(), 2249.0090580702745);

  // Test when abs_charge is 0 - should return 0.0
  LogMzPeak zero_charge_peak(tmp_p1, true);
  zero_charge_peak.abs_charge = 0;
  TEST_REAL_SIMILAR(zero_charge_peak.getUnchargedMass(), 0.0);

  // Test when mass is already set (mass > 0) - should return mass directly
  LogMzPeak preset_mass_peak(tmp_p1, true);
  preset_mass_peak.abs_charge = 2;
  preset_mass_peak.mass = 5000.0;
  TEST_REAL_SIMILAR(preset_mass_peak.getUnchargedMass(), 5000.0);

  // Test negative mode
  LogMzPeak neg_peak(tmp_p1, false);
  neg_peak.abs_charge = 2;
  double neg_mass = neg_peak.getUnchargedMass();
  // Negative mode should give different mass due to -PROTON_MASS_U
  TEST_NOT_EQUAL(neg_mass, test_peak.getUnchargedMass());
}
END_SECTION

LogMzPeak test_peak2(test_peak);
test_peak2.logMz = 8.;
START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator<(const LogMzPeak &a) const))
{
  // Test when logMz values are different
  bool is_p2_larger = test_peak < test_peak2;
  TEST_EQUAL(is_p2_larger, true);

  // Test when logMz values are equal - falls back to intensity comparison
  LogMzPeak peak_a(tmp_p1, true);
  LogMzPeak peak_b(tmp_p1, true);
  peak_a.logMz = 5.0;
  peak_b.logMz = 5.0;
  peak_a.intensity = 100.0f;
  peak_b.intensity = 200.0f;
  TEST_EQUAL(peak_a < peak_b, true);
  TEST_EQUAL(peak_b < peak_a, false);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator>(const LogMzPeak &a) const))
{
  // Test when logMz values are different
  bool is_p2_larger = test_peak2 > test_peak;
  TEST_EQUAL(is_p2_larger, true);

  // Test when logMz values are equal - falls back to intensity comparison
  LogMzPeak peak_a(tmp_p1, true);
  LogMzPeak peak_b(tmp_p1, true);
  peak_a.logMz = 5.0;
  peak_b.logMz = 5.0;
  peak_a.intensity = 200.0f;
  peak_b.intensity = 100.0f;
  TEST_EQUAL(peak_a > peak_b, true);
  TEST_EQUAL(peak_b > peak_a, false);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::LogMzPeak] bool operator==(const LogMzPeak &other) const))
{
  test_peak2 = test_peak;
  bool are_two_ps_same = test_peak2 == test_peak;
  TEST_EQUAL(are_two_ps_same, true);

  // Test inequality when logMz same but intensity different
  LogMzPeak peak_a(tmp_p1, true);
  LogMzPeak peak_b(tmp_p1, true);
  peak_a.logMz = 5.0;
  peak_b.logMz = 5.0;
  peak_a.intensity = 100.0f;
  peak_b.intensity = 200.0f;
  TEST_EQUAL(peak_a == peak_b, false);

  // Test inequality when intensity same but logMz different
  peak_b.intensity = 100.0f;
  peak_b.logMz = 6.0;
  TEST_EQUAL(peak_a == peak_b, false);
}
END_SECTION
///


/// testing PrecalculatedAveragine
START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine()))
{
  PrecalculatedAveragine* p_avg_ptr = new PrecalculatedAveragine();
  PrecalculatedAveragine* p_avg_null_ptr = 0;

  TEST_NOT_EQUAL(p_avg_ptr, p_avg_null_ptr)
  delete p_avg_ptr;
}
END_SECTION

// test data
CoarseIsotopePatternGenerator generator = CoarseIsotopePatternGenerator();
PrecalculatedAveragine p_avg_test;
START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine(const double min_mass, const double max_mass, const double delta, CoarseIsotopePatternGenerator& generator, const bool use_RNA_averagine)))
{
  p_avg_test = PrecalculatedAveragine(50, 100, 25, generator, false, -1);
  Size temp_a_idx = p_avg_test.getApexIndex(75);
  double temp_m_diff = p_avg_test.getAverageMassDelta(75);
  TEST_EQUAL(temp_a_idx, 0);
  TOLERANCE_ABSOLUTE(0.3);
  TEST_REAL_SIMILAR(temp_m_diff, 0.04);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] IsotopeDistribution get(const double mass) const))
{
  IsotopeDistribution tmp_iso = p_avg_test.get(60);
  TOLERANCE_ABSOLUTE(2);
  TEST_REAL_SIMILAR(tmp_iso.getMin(), 53.);
  TEST_REAL_SIMILAR(tmp_iso.getMax(), 59.);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] void setMaxIsotopeIndex(const int index)))
{
  p_avg_test.setMaxIsotopeIndex(4);
  TEST_EQUAL(p_avg_test.getMaxIsotopeIndex(), 4);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] int getMaxIsotopeIndex() const))
{
  int tmp_max_idx = p_avg_test.getMaxIsotopeIndex();
  TEST_EQUAL(tmp_max_idx, 4);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getLeftCountFromApex(const double mass) const))
{
  Size tmp_left = p_avg_test.getLeftCountFromApex(75);
  TEST_EQUAL(tmp_left, 2);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getRightCountFromApex(const double mass) const))
{
  Size temp_right = p_avg_test.getRightCountFromApex(75);
  TEST_EQUAL(temp_right, 2);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getApexIndex(const double mass) const))
{
  Size tmp_apex = p_avg_test.getApexIndex(75);
  TEST_EQUAL(tmp_apex, 0);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] double getAverageMassDelta(const double mass) const))
{
  double tmp_m_delta = p_avg_test.getAverageMassDelta(50);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_m_delta, 0.025);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] double getMostAbundantMassDelta(const double mass) const))
{
  double tmp_m_delta = p_avg_test.getMostAbundantMassDelta(1000);
  TOLERANCE_ABSOLUTE(0.1);
  TEST_REAL_SIMILAR(tmp_m_delta, 0);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] Size getLastIndex(const double mass) const))
{
  double last_index = p_avg_test.getLastIndex(50);
  TEST_EQUAL(last_index, 2);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine(const PrecalculatedAveragine&)))
{
  // Test copy constructor
  PrecalculatedAveragine p_avg_copy(p_avg_test);
  TEST_EQUAL(p_avg_copy.getApexIndex(75), p_avg_test.getApexIndex(75));
  TEST_EQUAL(p_avg_copy.getMaxIsotopeIndex(), p_avg_test.getMaxIsotopeIndex());
  TEST_REAL_SIMILAR(p_avg_copy.getAverageMassDelta(75), p_avg_test.getAverageMassDelta(75));
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine(PrecalculatedAveragine&&)))
{
  // Test move constructor
  CoarseIsotopePatternGenerator gen;
  PrecalculatedAveragine p_avg_temp(50, 100, 25, gen, false, -1);
  Size original_apex = p_avg_temp.getApexIndex(75);
  PrecalculatedAveragine p_avg_moved(std::move(p_avg_temp));
  TEST_EQUAL(p_avg_moved.getApexIndex(75), original_apex);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine& operator=(const PrecalculatedAveragine&)))
{
  // Test copy assignment operator
  PrecalculatedAveragine p_avg_assigned;
  p_avg_assigned = p_avg_test;
  TEST_EQUAL(p_avg_assigned.getApexIndex(75), p_avg_test.getApexIndex(75));
  TEST_EQUAL(p_avg_assigned.getMaxIsotopeIndex(), p_avg_test.getMaxIsotopeIndex());
  TEST_REAL_SIMILAR(p_avg_assigned.getAverageMassDelta(75), p_avg_test.getAverageMassDelta(75));
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine& operator=(PrecalculatedAveragine&&)))
{
  // Test move assignment operator
  CoarseIsotopePatternGenerator gen;
  PrecalculatedAveragine p_avg_temp(50, 100, 25, gen, false, -1);
  Size original_apex = p_avg_temp.getApexIndex(75);
  PrecalculatedAveragine p_avg_move_assigned;
  p_avg_move_assigned = std::move(p_avg_temp);
  TEST_EQUAL(p_avg_move_assigned.getApexIndex(75), original_apex);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] double getSNRMultiplicationFactor(const double mass) const))
{
  double snr_factor = p_avg_test.getSNRMultiplicationFactor(75);
  // SNR multiplication factor should be positive
  TEST_EQUAL(snr_factor > 0, true);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine with RNA averagine))
{
  // Test RNA averagine mode
  CoarseIsotopePatternGenerator gen_rna;
  PrecalculatedAveragine p_avg_rna(50, 100, 25, gen_rna, true, -1);
  // RNA averagine should produce different isotope patterns than peptide
  IsotopeDistribution iso_rna = p_avg_rna.get(75);
  IsotopeDistribution iso_peptide = p_avg_test.get(75);
  // The distributions should exist and be valid
  TEST_EQUAL(iso_rna.size() > 0, true);
  TEST_EQUAL(iso_peptide.size() > 0, true);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::PrecalculatedAveragine] PrecalculatedAveragine with decoy isotope distance))
{
  // Test decoy mode with positive decoy_iso_distance
  CoarseIsotopePatternGenerator gen_decoy;
  PrecalculatedAveragine p_avg_decoy(50, 100, 25, gen_decoy, false, 1.5);
  // Decoy isotope patterns should be generated
  IsotopeDistribution iso_decoy = p_avg_decoy.get(75);
  TEST_EQUAL(iso_decoy.size() > 0, true);
  // The decoy distribution should differ from normal due to scaled m/z values
  IsotopeDistribution iso_normal = p_avg_test.get(75);
  // Decoy patterns have scaled isotope distances
  TEST_NOT_EQUAL(iso_decoy.getMin(), iso_normal.getMin());
}
END_SECTION
///


/// testing MassFeature
START_SECTION(([FLASHDeconvHelperStructs::MassFeature] MassFeature default values))
{
  MassFeature mf;
  // Test that MassFeature can be default constructed
  MassFeature* mf_ptr = new MassFeature();
  MassFeature* mf_null_ptr = nullptr;
  TEST_NOT_EQUAL(mf_ptr, mf_null_ptr);
  delete mf_ptr;
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::MassFeature] bool operator<(const MassFeature &a) const))
{
  MassFeature mf1;
  MassFeature mf2;
  mf1.avg_mass = 1000.0;
  mf2.avg_mass = 2000.0;
  TEST_EQUAL(mf1 < mf2, true);
  TEST_EQUAL(mf2 < mf1, false);

  // Test equal masses
  mf2.avg_mass = 1000.0;
  TEST_EQUAL(mf1 < mf2, false);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::MassFeature] bool operator>(const MassFeature &a) const))
{
  MassFeature mf1;
  MassFeature mf2;
  mf1.avg_mass = 2000.0;
  mf2.avg_mass = 1000.0;
  TEST_EQUAL(mf1 > mf2, true);
  TEST_EQUAL(mf2 > mf1, false);

  // Test equal masses
  mf2.avg_mass = 2000.0;
  TEST_EQUAL(mf1 > mf2, false);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::MassFeature] bool operator==(const MassFeature &other) const))
{
  MassFeature mf1;
  MassFeature mf2;
  mf1.avg_mass = 1000.0;
  mf2.avg_mass = 1000.0;
  TEST_EQUAL(mf1 == mf2, true);

  mf2.avg_mass = 2000.0;
  TEST_EQUAL(mf1 == mf2, false);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::MassFeature] member variables))
{
  MassFeature mf;
  // Test setting and getting member variables
  mf.index = 42;
  mf.iso_offset = 1;
  mf.scan_number = 100;
  mf.min_scan_number = 90;
  mf.max_scan_number = 110;
  mf.rep_charge = 5;
  mf.avg_mass = 15000.0;
  mf.min_charge = 3;
  mf.max_charge = 10;
  mf.charge_count = 8;
  mf.isotope_score = 0.95;
  mf.qscore = 0.88;
  mf.rep_mz = 1500.5;
  mf.is_decoy = false;
  mf.ms_level = 1;

  TEST_EQUAL(mf.index, 42);
  TEST_EQUAL(mf.iso_offset, 1);
  TEST_EQUAL(mf.scan_number, 100);
  TEST_EQUAL(mf.min_scan_number, 90);
  TEST_EQUAL(mf.max_scan_number, 110);
  TEST_EQUAL(mf.rep_charge, 5);
  TEST_REAL_SIMILAR(mf.avg_mass, 15000.0);
  TEST_EQUAL(mf.min_charge, 3);
  TEST_EQUAL(mf.max_charge, 10);
  TEST_EQUAL(mf.charge_count, 8);
  TEST_REAL_SIMILAR(mf.isotope_score, 0.95);
  TEST_REAL_SIMILAR(mf.qscore, 0.88);
  TEST_REAL_SIMILAR(mf.rep_mz, 1500.5);
  TEST_EQUAL(mf.is_decoy, false);
  TEST_EQUAL(mf.ms_level, 1);

  // Test per_charge_intensity and per_isotope_intensity vectors
  mf.per_charge_intensity = {100.0f, 200.0f, 300.0f};
  mf.per_isotope_intensity = {50.0f, 100.0f, 75.0f, 25.0f};
  TEST_EQUAL(mf.per_charge_intensity.size(), 3);
  TEST_EQUAL(mf.per_isotope_intensity.size(), 4);
  TEST_REAL_SIMILAR(mf.per_charge_intensity[1], 200.0f);
  TEST_REAL_SIMILAR(mf.per_isotope_intensity[2], 75.0f);
}
END_SECTION
///


/// testing IsobaricQuantities
START_SECTION(([FLASHDeconvHelperStructs::IsobaricQuantities] IsobaricQuantities default values))
{
  IsobaricQuantities iq;
  // Test that IsobaricQuantities can be default constructed
  IsobaricQuantities* iq_ptr = new IsobaricQuantities();
  IsobaricQuantities* iq_null_ptr = nullptr;
  TEST_NOT_EQUAL(iq_ptr, iq_null_ptr);
  delete iq_ptr;
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::IsobaricQuantities] bool empty() const))
{
  IsobaricQuantities iq;
  // Should be empty initially (quantities vector is empty)
  TEST_EQUAL(iq.empty(), true);

  // Add quantities
  iq.quantities.push_back(100.0);
  TEST_EQUAL(iq.empty(), false);

  // Clear quantities
  iq.quantities.clear();
  TEST_EQUAL(iq.empty(), true);
}
END_SECTION

START_SECTION(([FLASHDeconvHelperStructs::IsobaricQuantities] member variables))
{
  IsobaricQuantities iq;
  // Test setting and getting member variables
  iq.scan = 500;
  iq.rt = 120.5;
  iq.precursor_mz = 750.25;
  iq.precursor_mass = 1498.48;
  iq.quantities = {100.0, 200.0, 150.0, 175.0};
  iq.merged_quantities = {450.0, 375.0};

  TEST_EQUAL(iq.scan, 500);
  TEST_REAL_SIMILAR(iq.rt, 120.5);
  TEST_REAL_SIMILAR(iq.precursor_mz, 750.25);
  TEST_REAL_SIMILAR(iq.precursor_mass, 1498.48);
  TEST_EQUAL(iq.quantities.size(), 4);
  TEST_EQUAL(iq.merged_quantities.size(), 2);
  TEST_REAL_SIMILAR(iq.quantities[0], 100.0);
  TEST_REAL_SIMILAR(iq.quantities[3], 175.0);
  TEST_REAL_SIMILAR(iq.merged_quantities[0], 450.0);
}
END_SECTION
///

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST