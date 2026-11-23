// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <cmath>

using namespace OpenMS;

START_TEST(MultiplexFiltering, "$Id$")

// read data
MSExperiment exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexFiltering.mzML"), exp);
exp.updateRanges();

// pick data
PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);
std::vector<PeakPickerHiRes::PeakBoundary> boundaries;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;
MSExperiment exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

// set parameters
int charge_min = 1;
int charge_max = 6;
int isotopes_per_peptide_min = 3;
int isotopes_per_peptide_max = 6;
double intensity_cutoff = 10.0;
double rt_band = 2.0;
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double peptide_similarity = 0.8;
double averagine_similarity = 0.75;
double averagine_similarity_scaling = 0.75;
std::string averagine_type="peptide";

// construct list of peak patterns
MultiplexDeltaMasses shifts1;
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(8.0443702794, "Arg8"));
MultiplexDeltaMasses shifts2;
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
MultiplexDeltaMasses::LabelSet label_set;
label_set.insert("Arg8");
label_set.insert("Arg8");
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(2*8.0443702794, label_set));
std::vector<MultiplexIsotopicPeakPattern> patterns;
for (int c = charge_max; c >= charge_min; --c)
{
    MultiplexIsotopicPeakPattern pattern1(c, isotopes_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    MultiplexIsotopicPeakPattern pattern2(c, isotopes_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
}

MultiplexFiltering* nullPointer = nullptr;
MultiplexFiltering* ptr;

START_SECTION(MultiplexFiltering(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, std::string averagine_type="peptide"))
  MultiplexFiltering filtering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  ptr = new MultiplexFiltering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

START_SECTION((static double deuteriumLowerFraction(unsigned deuterium, double enrichment, const IsotopeDistribution& isotopes)))
{
  IsotopeDistribution isotopes;
  isotopes.set(IsotopeDistribution::ContainerType{Peak1D(1000.0, 0.5f), Peak1D(1001.0, 0.3f), Peak1D(1002.0, 0.2f)});
  // four 2H at 96 %: q_1 = 4 (0.04/0.96), q_2 = 6 (0.04/0.96)^2, q_3 = 4 (0.04/0.96)^3;
  // (q_1 I_0 + q_2 I_1 + q_3 I_2) / (I_0 + q_1 I_1 + q_2 I_2)
  TEST_REAL_SIMILAR(MultiplexFiltering::deuteriumLowerFraction(4, 0.96, isotopes), 0.156708595)
  // about r / (1 + r I_1/I_0) with r = n (1 - p) / p
  TEST_EQUAL(std::fabs(MultiplexFiltering::deuteriumLowerFraction(4, 0.96, isotopes) - (1.0 / 6) / (1 + 0.6 / 6)) < 0.01, true)
  // more 2H atoms (two Lys4), more of it
  TEST_EQUAL(MultiplexFiltering::deuteriumLowerFraction(8, 0.96, isotopes) > MultiplexFiltering::deuteriumLowerFraction(4, 0.96, isotopes), true)
  // none without 2H, at complete enrichment, or without an isotope distribution to place it
  TEST_EQUAL(MultiplexFiltering::deuteriumLowerFraction(0, 0.96, isotopes), 0.0)
  TEST_EQUAL(MultiplexFiltering::deuteriumLowerFraction(4, 1.0, isotopes), 0.0)
  IsotopeDistribution single;
  single.set(IsotopeDistribution::ContainerType{Peak1D(1000.0, 1.0f)});
  TEST_EQUAL(MultiplexFiltering::deuteriumLowerFraction(4, 0.96, single), 0.0)
}
END_SECTION

START_SECTION((static double deuteriumAllowance(unsigned deuterium, double enrichment, const IsotopeDistribution& isotopes)))
{
  IsotopeDistribution isotopes;
  isotopes.set(IsotopeDistribution::ContainerType{Peak1D(1000.0, 0.5f), Peak1D(1001.0, 0.3f), Peak1D(1002.0, 0.2f)});
  // below the cap, the expected lower isotopologue
  TEST_REAL_SIMILAR(MultiplexFiltering::deuteriumAllowance(4, 0.96, isotopes), MultiplexFiltering::deuteriumLowerFraction(4, 0.96, isotopes))
  // at most half the way from 0.3 to I_0/I_1 = 5/3, which a pattern one isotope too high would show
  TEST_EQUAL(MultiplexFiltering::deuteriumLowerFraction(18, 0.6, isotopes) > 0.5 * (0.5 / 0.3 - 0.3), true)
  TEST_REAL_SIMILAR(MultiplexFiltering::deuteriumAllowance(18, 0.6, isotopes), 0.5 * (0.5 / 0.3 - 0.3))
  // none where I_0/I_1 is not above 0.3, without 2H or at complete enrichment
  IsotopeDistribution heavy;
  heavy.set(IsotopeDistribution::ContainerType{Peak1D(1000.0, 0.1f), Peak1D(1001.0, 0.5f), Peak1D(1002.0, 0.4f)});
  TEST_EQUAL(MultiplexFiltering::deuteriumAllowance(4, 0.96, heavy), 0.0)
  TEST_EQUAL(MultiplexFiltering::deuteriumAllowance(0, 0.96, isotopes), 0.0)
  TEST_EQUAL(MultiplexFiltering::deuteriumAllowance(4, 1.0, isotopes), 0.0)
}
END_SECTION

START_SECTION((void setDeuteriumLabels(const std::map<std::string, DeuteriumLabel>& labels)))
{
  MultiplexFiltering filtering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  filtering.setDeuteriumLabels({{"Arg8", {8, 0.97}}});
  filtering.setDeuteriumLabels({});
  TEST_EXCEPTION(Exception::InvalidValue, filtering.setDeuteriumLabels({{"Arg8", {8, 0.0}}}))
  TEST_EXCEPTION(Exception::InvalidValue, filtering.setDeuteriumLabels({{"Arg8", {8, 1.1}}}))
}
END_SECTION

END_TEST
