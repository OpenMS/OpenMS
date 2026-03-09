// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <cmath>

///////////////////////////

START_TEST(FLASHHelperClasses, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef FLASHHelperClasses::LogMzPeak LogMzPeak;
typedef FLASHHelperClasses::MassFeature MassFeature;
typedef FLASHHelperClasses::IsobaricQuantities IsobaricQuantities;
typedef FLASHHelperClasses::PrecalculatedAveragine PrecalculatedAveragine;

/////////////////////////////////////////////////////////////
// Static utility function tests
/////////////////////////////////////////////////////////////

START_SECTION(static float getChargeMass(bool positive_ioniziation_mode))
{
  // Positive mode should return proton mass
  float pos_mass = FLASHHelperClasses::getChargeMass(true);
  TEST_REAL_SIMILAR(pos_mass, Constants::PROTON_MASS_U)

  // Negative mode should return negative proton mass
  float neg_mass = FLASHHelperClasses::getChargeMass(false);
  TEST_REAL_SIMILAR(neg_mass, -Constants::PROTON_MASS_U)
}
END_SECTION

START_SECTION(static double getLogMz(double mz, bool positive))
{
  // Test positive mode
  double mz = 500.0;
  double log_mz_pos = FLASHHelperClasses::getLogMz(mz, true);
  double expected_pos = std::log(mz - Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(log_mz_pos, expected_pos)

  // Test negative mode
  double log_mz_neg = FLASHHelperClasses::getLogMz(mz, false);
  double expected_neg = std::log(mz + Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(log_mz_neg, expected_neg)

  // Test with different m/z values
  double mz2 = 1000.0;
  double log_mz2 = FLASHHelperClasses::getLogMz(mz2, true);
  TEST_REAL_SIMILAR(log_mz2, std::log(mz2 - Constants::PROTON_MASS_U))
}
END_SECTION

/////////////////////////////////////////////////////////////
// LogMzPeak tests
/////////////////////////////////////////////////////////////

LogMzPeak* lmp_ptr = nullptr;
LogMzPeak* lmp_null_ptr = nullptr;

START_SECTION(LogMzPeak())
{
  lmp_ptr = new LogMzPeak();
  TEST_NOT_EQUAL(lmp_ptr, lmp_null_ptr)

  // Check default values
  TEST_REAL_SIMILAR(lmp_ptr->mz, 0.0)
  TEST_REAL_SIMILAR(lmp_ptr->intensity, 0.0f)
  TEST_REAL_SIMILAR(lmp_ptr->logMz, -1000.0)
  TEST_REAL_SIMILAR(lmp_ptr->mass, 0.0)
  TEST_EQUAL(lmp_ptr->abs_charge, 0)
  TEST_EQUAL(lmp_ptr->is_positive, true)
  TEST_EQUAL(lmp_ptr->isotopeIndex, -1)
}
END_SECTION

START_SECTION(~LogMzPeak())
{
  delete lmp_ptr;
}
END_SECTION

START_SECTION(explicit LogMzPeak(const Peak1D& peak, bool positive))
{
  // Test with positive mode
  Peak1D peak;
  peak.setMZ(500.5);
  peak.setIntensity(10000.0f);

  LogMzPeak lmp_pos(peak, true);
  TEST_REAL_SIMILAR(lmp_pos.mz, 500.5)
  TEST_REAL_SIMILAR(lmp_pos.intensity, 10000.0f)
  TEST_EQUAL(lmp_pos.is_positive, true)
  TEST_EQUAL(lmp_pos.abs_charge, 0)
  TEST_EQUAL(lmp_pos.isotopeIndex, 0)
  double expected_logMz = std::log(500.5 - Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(lmp_pos.logMz, expected_logMz)

  // Test with negative mode
  LogMzPeak lmp_neg(peak, false);
  TEST_REAL_SIMILAR(lmp_neg.mz, 500.5)
  TEST_EQUAL(lmp_neg.is_positive, false)
  double expected_logMz_neg = std::log(500.5 + Constants::PROTON_MASS_U);
  TEST_REAL_SIMILAR(lmp_neg.logMz, expected_logMz_neg)
}
END_SECTION

START_SECTION(LogMzPeak(const LogMzPeak&))
{
  Peak1D peak;
  peak.setMZ(750.25);
  peak.setIntensity(5000.0f);
  LogMzPeak original(peak, true);
  original.abs_charge = 3;
  original.isotopeIndex = 2;
  original.mass = 2248.0;

  LogMzPeak copy(original);
  TEST_REAL_SIMILAR(copy.mz, 750.25)
  TEST_REAL_SIMILAR(copy.intensity, 5000.0f)
  TEST_EQUAL(copy.abs_charge, 3)
  TEST_EQUAL(copy.isotopeIndex, 2)
  TEST_REAL_SIMILAR(copy.mass, 2248.0)
  TEST_EQUAL(copy.is_positive, true)
  TEST_REAL_SIMILAR(copy.logMz, original.logMz)
}
END_SECTION

START_SECTION(double getUnchargedMass() const)
{
  Peak1D peak;
  peak.setMZ(500.5);
  peak.setIntensity(1000.0f);
  LogMzPeak lmp(peak, true);

  // With no charge set, should return 0
  TEST_REAL_SIMILAR(lmp.getUnchargedMass(), 0.0)

  // Set charge and mass = 0 (will calculate from mz)
  lmp.abs_charge = 2;
  lmp.mass = 0;
  double expected_mass = (500.5 - Constants::PROTON_MASS_U) * 2;
  TEST_REAL_SIMILAR(lmp.getUnchargedMass(), expected_mass)

  // With mass already set > 0, should return that mass
  lmp.mass = 999.5;
  TEST_REAL_SIMILAR(lmp.getUnchargedMass(), 999.5)

  // Test negative mode
  Peak1D peak2;
  peak2.setMZ(500.5);
  peak2.setIntensity(1000.0f);
  LogMzPeak lmp_neg(peak2, false);
  lmp_neg.abs_charge = 2;
  lmp_neg.mass = 0;
  double expected_neg = (500.5 + Constants::PROTON_MASS_U) * 2;
  TEST_REAL_SIMILAR(lmp_neg.getUnchargedMass(), expected_neg)
}
END_SECTION

START_SECTION(bool operator<(const LogMzPeak& a) const)
{
  LogMzPeak lmp1, lmp2, lmp3;

  // Different logMz values
  lmp1.logMz = 5.0;
  lmp1.intensity = 100.0f;
  lmp2.logMz = 6.0;
  lmp2.intensity = 100.0f;

  TEST_EQUAL(lmp1 < lmp2, true)
  TEST_EQUAL(lmp2 < lmp1, false)

  // Same logMz, different intensity
  lmp3.logMz = 5.0;
  lmp3.intensity = 200.0f;
  TEST_EQUAL(lmp1 < lmp3, true)   // same logMz, but lmp1 has lower intensity
  TEST_EQUAL(lmp3 < lmp1, false)
}
END_SECTION

START_SECTION(bool operator>(const LogMzPeak& a) const)
{
  LogMzPeak lmp1, lmp2, lmp3;

  // Different logMz values
  lmp1.logMz = 6.0;
  lmp1.intensity = 100.0f;
  lmp2.logMz = 5.0;
  lmp2.intensity = 100.0f;

  TEST_EQUAL(lmp1 > lmp2, true)
  TEST_EQUAL(lmp2 > lmp1, false)

  // Same logMz, different intensity
  lmp3.logMz = 6.0;
  lmp3.intensity = 50.0f;
  TEST_EQUAL(lmp1 > lmp3, true)   // same logMz, but lmp1 has higher intensity
  TEST_EQUAL(lmp3 > lmp1, false)
}
END_SECTION

START_SECTION(bool operator==(const LogMzPeak& other) const)
{
  LogMzPeak lmp1, lmp2, lmp3;

  lmp1.logMz = 5.5;
  lmp1.intensity = 100.0f;
  lmp2.logMz = 5.5;
  lmp2.intensity = 100.0f;
  lmp3.logMz = 5.5;
  lmp3.intensity = 200.0f;

  TEST_EQUAL(lmp1 == lmp2, true)
  TEST_EQUAL(lmp1 == lmp3, false)  // same logMz but different intensity
}
END_SECTION

/////////////////////////////////////////////////////////////
// MassFeature tests
/////////////////////////////////////////////////////////////

START_SECTION([MassFeature] bool operator<(const MassFeature& a) const)
{
  MassFeature mf1, mf2;
  mf1.avg_mass = 10000.0;
  mf2.avg_mass = 15000.0;

  TEST_EQUAL(mf1 < mf2, true)
  TEST_EQUAL(mf2 < mf1, false)
}
END_SECTION

START_SECTION([MassFeature] bool operator>(const MassFeature& a) const)
{
  MassFeature mf1, mf2;
  mf1.avg_mass = 15000.0;
  mf2.avg_mass = 10000.0;

  TEST_EQUAL(mf1 > mf2, true)
  TEST_EQUAL(mf2 > mf1, false)
}
END_SECTION

START_SECTION([MassFeature] bool operator==(const MassFeature& other) const)
{
  MassFeature mf1, mf2, mf3;
  mf1.avg_mass = 10000.0;
  mf2.avg_mass = 10000.0;
  mf3.avg_mass = 10001.0;

  TEST_EQUAL(mf1 == mf2, true)
  TEST_EQUAL(mf1 == mf3, false)
}
END_SECTION

START_SECTION([MassFeature] member variables)
{
  MassFeature mf;

  // Test setting and accessing all member variables
  mf.index = 42;
  mf.iso_offset = 1;
  mf.scan_number = 100;
  mf.min_scan_number = 95;
  mf.max_scan_number = 105;
  mf.rep_charge = 10;
  mf.avg_mass = 15000.5;
  mf.min_charge = 5;
  mf.max_charge = 15;
  mf.charge_count = 11;
  mf.isotope_score = 0.95;
  mf.qscore = 0.85;
  mf.rep_mz = 1500.05;
  mf.is_decoy = false;
  mf.ms_level = 1;

  TEST_EQUAL(mf.index, 42)
  TEST_EQUAL(mf.iso_offset, 1)
  TEST_EQUAL(mf.scan_number, 100)
  TEST_EQUAL(mf.min_scan_number, 95)
  TEST_EQUAL(mf.max_scan_number, 105)
  TEST_EQUAL(mf.rep_charge, 10)
  TEST_REAL_SIMILAR(mf.avg_mass, 15000.5)
  TEST_EQUAL(mf.min_charge, 5)
  TEST_EQUAL(mf.max_charge, 15)
  TEST_EQUAL(mf.charge_count, 11)
  TEST_REAL_SIMILAR(mf.isotope_score, 0.95)
  TEST_REAL_SIMILAR(mf.qscore, 0.85)
  TEST_REAL_SIMILAR(mf.rep_mz, 1500.05)
  TEST_EQUAL(mf.is_decoy, false)
  TEST_EQUAL(mf.ms_level, 1)

  // Test per_charge_intensity and per_isotope_intensity vectors
  mf.per_charge_intensity.resize(20, 0.0f);
  mf.per_isotope_intensity.resize(10, 0.0f);
  mf.per_charge_intensity[10] = 5000.0f;
  mf.per_isotope_intensity[0] = 10000.0f;

  TEST_REAL_SIMILAR(mf.per_charge_intensity[10], 5000.0f)
  TEST_REAL_SIMILAR(mf.per_isotope_intensity[0], 10000.0f)
}
END_SECTION

/////////////////////////////////////////////////////////////
// IsobaricQuantities tests
/////////////////////////////////////////////////////////////

START_SECTION([IsobaricQuantities] bool empty() const)
{
  IsobaricQuantities iq;

  // Initially empty (no quantities)
  TEST_EQUAL(iq.empty(), true)

  // Add a quantity
  iq.quantities.push_back(100.0);
  TEST_EQUAL(iq.empty(), false)

  // Clear and verify empty again
  iq.quantities.clear();
  TEST_EQUAL(iq.empty(), true)
}
END_SECTION

START_SECTION([IsobaricQuantities] member variables)
{
  IsobaricQuantities iq;

  // Test setting and accessing all member variables
  iq.scan = 150;
  iq.rt = 600.5;
  iq.precursor_mz = 800.4;
  iq.precursor_mass = 7996.0;

  TEST_EQUAL(iq.scan, 150)
  TEST_REAL_SIMILAR(iq.rt, 600.5)
  TEST_REAL_SIMILAR(iq.precursor_mz, 800.4)
  TEST_REAL_SIMILAR(iq.precursor_mass, 7996.0)

  // Test quantities vectors
  iq.quantities = {100.0, 200.0, 300.0, 400.0};
  iq.merged_quantities = {150.0, 350.0};

  TEST_EQUAL(iq.quantities.size(), 4)
  TEST_REAL_SIMILAR(iq.quantities[0], 100.0)
  TEST_REAL_SIMILAR(iq.quantities[3], 400.0)
  TEST_EQUAL(iq.merged_quantities.size(), 2)
  TEST_REAL_SIMILAR(iq.merged_quantities[0], 150.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// PrecalculatedAveragine tests
/////////////////////////////////////////////////////////////

START_SECTION(PrecalculatedAveragine())
{
  PrecalculatedAveragine* pa_ptr = new PrecalculatedAveragine();
  TEST_NOT_EQUAL(pa_ptr, nullptr)
  delete pa_ptr;
}
END_SECTION

START_SECTION(PrecalculatedAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance = -1))
{
  CoarseIsotopePatternGenerator generator(100);

  // Create PrecalculatedAveragine for peptide
  PrecalculatedAveragine pa(500.0, 5000.0, 100.0, generator, false);

  // Test that it was created and can return isotope distributions
  IsotopeDistribution iso = pa.get(1000.0);
  TEST_EQUAL(iso.size() > 0, true)

  // Test RNA averagine
  PrecalculatedAveragine pa_rna(500.0, 5000.0, 100.0, generator, true);
  IsotopeDistribution iso_rna = pa_rna.get(1000.0);
  TEST_EQUAL(iso_rna.size() > 0, true)
}
END_SECTION

START_SECTION(PrecalculatedAveragine(const PrecalculatedAveragine&))
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine original(500.0, 5000.0, 100.0, generator, false);

  PrecalculatedAveragine copy(original);
  IsotopeDistribution iso_orig = original.get(2000.0);
  IsotopeDistribution iso_copy = copy.get(2000.0);

  TEST_EQUAL(iso_orig.size(), iso_copy.size())
  if (iso_orig.size() > 0 && iso_copy.size() > 0)
  {
    TEST_REAL_SIMILAR(iso_orig[0].getIntensity(), iso_copy[0].getIntensity())
  }
}
END_SECTION

START_SECTION(PrecalculatedAveragine(PrecalculatedAveragine&& other) noexcept)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine original(500.0, 5000.0, 100.0, generator, false);
  Size orig_apex = original.getApexIndex(2000.0);

  PrecalculatedAveragine moved(std::move(original));
  TEST_EQUAL(moved.getApexIndex(2000.0), orig_apex)
}
END_SECTION

START_SECTION(PrecalculatedAveragine& operator=(const PrecalculatedAveragine& pc))
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine original(500.0, 5000.0, 100.0, generator, false);
  PrecalculatedAveragine assigned;

  assigned = original;
  IsotopeDistribution iso_orig = original.get(2000.0);
  IsotopeDistribution iso_assigned = assigned.get(2000.0);

  TEST_EQUAL(iso_orig.size(), iso_assigned.size())
}
END_SECTION

START_SECTION(PrecalculatedAveragine& operator=(PrecalculatedAveragine&& pc) noexcept)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine original(500.0, 5000.0, 100.0, generator, false);
  Size orig_apex = original.getApexIndex(2000.0);

  PrecalculatedAveragine moved_assigned;
  moved_assigned = std::move(original);
  TEST_EQUAL(moved_assigned.getApexIndex(2000.0), orig_apex)
}
END_SECTION

START_SECTION(IsotopeDistribution get(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  // Test various masses
  IsotopeDistribution iso1 = pa.get(1000.0);
  IsotopeDistribution iso2 = pa.get(5000.0);
  IsotopeDistribution iso3 = pa.get(9000.0);

  TEST_EQUAL(iso1.size() > 0, true)
  TEST_EQUAL(iso2.size() > 0, true)
  TEST_EQUAL(iso3.size() > 0, true)

  // Higher mass should have more isotopes with significant intensity
  // (apex moves to higher index)
  Size apex1 = pa.getApexIndex(1000.0);
  Size apex2 = pa.getApexIndex(5000.0);
  Size apex3 = pa.getApexIndex(9000.0);

  TEST_EQUAL(apex1 <= apex2, true)
  TEST_EQUAL(apex2 <= apex3, true)

  // Test mass exceeding max - should return distribution for max mass
  IsotopeDistribution iso_exceed = pa.get(20000.0);
  TEST_EQUAL(iso_exceed.size() > 0, true)
}
END_SECTION

START_SECTION(size_t getMaxIsotopeIndex() const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 5000.0, 100.0, generator, false);
  pa.setMaxIsotopeIndex(50);

  TEST_EQUAL(pa.getMaxIsotopeIndex(), 50)
}
END_SECTION

START_SECTION(void setMaxIsotopeIndex(int index))
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 5000.0, 100.0, generator, false);

  pa.setMaxIsotopeIndex(30);
  TEST_EQUAL(pa.getMaxIsotopeIndex(), 30)

  pa.setMaxIsotopeIndex(100);
  TEST_EQUAL(pa.getMaxIsotopeIndex(), 100)
}
END_SECTION

START_SECTION(Size getLeftCountFromApex(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  Size left1 = pa.getLeftCountFromApex(1000.0);
  Size left2 = pa.getLeftCountFromApex(5000.0);

  // Should return reasonable values (at least minimum of 2)
  TEST_EQUAL(left1 >= 2, true)
  TEST_EQUAL(left2 >= 2, true)
}
END_SECTION

START_SECTION(Size getRightCountFromApex(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  Size right1 = pa.getRightCountFromApex(1000.0);
  Size right2 = pa.getRightCountFromApex(5000.0);

  // Should return reasonable values (at least minimum of 2)
  TEST_EQUAL(right1 >= 2, true)
  TEST_EQUAL(right2 >= 2, true)
}
END_SECTION

START_SECTION(Size getApexIndex(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  Size apex1 = pa.getApexIndex(1000.0);
  Size apex5 = pa.getApexIndex(5000.0);
  Size apex9 = pa.getApexIndex(9000.0);

  // Apex should shift to higher indices for higher masses
  TEST_EQUAL(apex1 <= apex5, true)
  TEST_EQUAL(apex5 <= apex9, true)

  // For small mass, apex should be at or near 0
  Size apex_small = pa.getApexIndex(600.0);
  TEST_EQUAL(apex_small <= 2, true)
}
END_SECTION

START_SECTION(Size getLastIndex(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  Size last1 = pa.getLastIndex(1000.0);
  Size last5 = pa.getLastIndex(5000.0);

  // Last index should be apex + right count
  Size apex1 = pa.getApexIndex(1000.0);
  Size right1 = pa.getRightCountFromApex(1000.0);
  TEST_EQUAL(last1, apex1 + right1)

  Size apex5 = pa.getApexIndex(5000.0);
  Size right5 = pa.getRightCountFromApex(5000.0);
  TEST_EQUAL(last5, apex5 + right5)
}
END_SECTION

START_SECTION(double getAverageMassDelta(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  double delta1 = pa.getAverageMassDelta(1000.0);
  double delta5 = pa.getAverageMassDelta(5000.0);
  double delta9 = pa.getAverageMassDelta(9000.0);

  // Average mass should be greater than monoisotopic mass (delta > 0)
  TEST_EQUAL(delta1 > 0, true)
  TEST_EQUAL(delta5 > 0, true)
  TEST_EQUAL(delta9 > 0, true)

  // Delta should increase with mass
  TEST_EQUAL(delta1 < delta5, true)
  TEST_EQUAL(delta5 < delta9, true)
}
END_SECTION

START_SECTION(double getMostAbundantMassDelta(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  double delta1 = pa.getMostAbundantMassDelta(1000.0);
  double delta5 = pa.getMostAbundantMassDelta(5000.0);
  double delta9 = pa.getMostAbundantMassDelta(9000.0);

  // Most abundant mass delta should be non-negative
  TEST_EQUAL(delta1 >= 0, true)
  TEST_EQUAL(delta5 >= 0, true)
  TEST_EQUAL(delta9 >= 0, true)

  // For small masses, monoisotopic is often most abundant (delta ~ 0)
  // For larger masses, delta should increase
  TEST_EQUAL(delta1 <= delta5, true)
}
END_SECTION

START_SECTION(double getSNRMultiplicationFactor(double mass) const)
{
  CoarseIsotopePatternGenerator generator(100);
  PrecalculatedAveragine pa(500.0, 10000.0, 100.0, generator, false);

  double snr1 = pa.getSNRMultiplicationFactor(1000.0);
  double snr5 = pa.getSNRMultiplicationFactor(5000.0);

  // SNR factor should be positive
  TEST_EQUAL(snr1 > 0, true)
  TEST_EQUAL(snr5 > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test decoy isotope pattern generation
/////////////////////////////////////////////////////////////

START_SECTION([PrecalculatedAveragine] decoy isotope patterns)
{
  CoarseIsotopePatternGenerator generator(100);

  // Create normal and decoy averagine patterns
  PrecalculatedAveragine pa_normal(500.0, 5000.0, 100.0, generator, false, -1);
  PrecalculatedAveragine pa_decoy(500.0, 5000.0, 100.0, generator, false, 1.5);

  IsotopeDistribution iso_normal = pa_normal.get(2000.0);
  IsotopeDistribution iso_decoy = pa_decoy.get(2000.0);

  // Both should have distributions
  TEST_EQUAL(iso_normal.size() > 0, true)
  TEST_EQUAL(iso_decoy.size() > 0, true)

  // Decoy pattern should be different (scaled isotope distances)
  // The m/z values will be different
  if (iso_normal.size() > 1 && iso_decoy.size() > 1)
  {
    double normal_spacing = iso_normal[1].getMZ() - iso_normal[0].getMZ();
    double decoy_spacing = iso_decoy[1].getMZ() - iso_decoy[0].getMZ();
    // Decoy spacing should be approximately 1.5x normal spacing
    TEST_REAL_SIMILAR(decoy_spacing / normal_spacing, 1.5)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Integration tests: combining multiple classes
/////////////////////////////////////////////////////////////

START_SECTION([Integration] LogMzPeak with charge and mass calculations)
{
  // Create a peak and test full workflow
  Peak1D peak;
  peak.setMZ(857.0789);  // m/z for a +3 charged peptide
  peak.setIntensity(50000.0f);

  LogMzPeak lmp(peak, true);
  lmp.abs_charge = 3;
  lmp.isotopeIndex = 0;

  // Calculate uncharged mass
  double mass = lmp.getUnchargedMass();
  // Expected: (857.0789 - 1.007276) * 3 = 2568.2148
  double expected_mass = (857.0789 - Constants::PROTON_MASS_U) * 3;
  TEST_REAL_SIMILAR(mass, expected_mass)

  // Set mass directly and verify
  lmp.mass = 2568.0;
  TEST_REAL_SIMILAR(lmp.getUnchargedMass(), 2568.0)
}
END_SECTION

START_SECTION([Integration] MassFeature sorting)
{
  // Create multiple MassFeatures and test sorting
  std::vector<MassFeature> features;

  MassFeature mf1, mf2, mf3;
  mf1.avg_mass = 15000.0;
  mf1.index = 1;
  mf2.avg_mass = 10000.0;
  mf2.index = 2;
  mf3.avg_mass = 20000.0;
  mf3.index = 3;

  features.push_back(mf1);
  features.push_back(mf2);
  features.push_back(mf3);

  // Sort by mass using operator<
  std::sort(features.begin(), features.end());

  TEST_REAL_SIMILAR(features[0].avg_mass, 10000.0)
  TEST_REAL_SIMILAR(features[1].avg_mass, 15000.0)
  TEST_REAL_SIMILAR(features[2].avg_mass, 20000.0)
}
END_SECTION

START_SECTION([Integration] LogMzPeak sorting)
{
  // Create multiple LogMzPeaks and test sorting
  std::vector<LogMzPeak> peaks;

  LogMzPeak lmp1, lmp2, lmp3;
  lmp1.logMz = 6.5;
  lmp1.intensity = 1000.0f;
  lmp2.logMz = 6.2;
  lmp2.intensity = 2000.0f;
  lmp3.logMz = 6.8;
  lmp3.intensity = 1500.0f;

  peaks.push_back(lmp1);
  peaks.push_back(lmp2);
  peaks.push_back(lmp3);

  // Sort by logMz using operator<
  std::sort(peaks.begin(), peaks.end());

  TEST_REAL_SIMILAR(peaks[0].logMz, 6.2)
  TEST_REAL_SIMILAR(peaks[1].logMz, 6.5)
  TEST_REAL_SIMILAR(peaks[2].logMz, 6.8)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
