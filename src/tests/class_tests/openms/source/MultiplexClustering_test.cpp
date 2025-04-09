// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>

#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

START_TEST(MultiplexFilteringProfile, "$Id$")

// read data
MSExperiment exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexClustering.mzML"), exp);
exp.updateRanges();

// pick data
PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;
MSExperiment exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

// set parameters
int charge_min = 1;
int charge_max = 4;
int isotopes_per_peptide_min = 3;
int isotopes_per_peptide_max = 6;
double intensity_cutoff = 10.0;
double rt_band = 3;
double rt_typical = 90;
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double peptide_similarity = 0.8;
double averagine_similarity = 0.75;
double averagine_similarity_scaling = 0.75;
String averagine_type="peptide";

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

MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
std::vector<MultiplexFilteredMSExperiment> filter_results = filtering.filter();

MultiplexClustering* nullPointer = nullptr;
MultiplexClustering* ptr;

START_SECTION(MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical))
    MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical);
    std::vector<std::map<int,GridBasedCluster> > cluster_results = clustering.cluster(filter_results);
    ptr = new MultiplexClustering(exp, exp_picked, boundaries_exp_s, rt_typical);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical);

START_SECTION(cluster(const std::vector<MultiplexFilteredMSExperiment>& filter_results))
    std::vector<std::map<int,GridBasedCluster> > cluster_results = clustering.cluster(filter_results);
    TEST_EQUAL(cluster_results[0].size(), 0);
    TEST_EQUAL(cluster_results[1].size(), 0);
    TEST_EQUAL(cluster_results[2].size(), 0);
    TEST_EQUAL(cluster_results[3].size(), 0);
    TEST_EQUAL(cluster_results[4].size(), 2);
    TEST_EQUAL(cluster_results[5].size(), 0);
    TEST_EQUAL(cluster_results[6].size(), 0);
    TEST_EQUAL(cluster_results[7].size(), 0);
END_SECTION

END_TEST
