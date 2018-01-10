// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>

using namespace OpenMS;

START_TEST(MultiplexFilteringProfile, "$Id$")

// read data
PeakMap exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexClustering.mzML"), exp);
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
PeakMap exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

// set parameters
int charge_min = 1;
int charge_max = 4;
int peaks_per_peptide_min = 3;
int peaks_per_peptide_max = 6;
bool missing_peaks = false;
double intensity_cutoff = 10.0;
double peptide_similarity = 0.8;
double averagine_similarity = 0.75;
double averagine_similarity_scaling = 0.75;
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double rt_typical = 90;
double rt_minimum = 5;

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
    MultiplexIsotopicPeakPattern pattern1(c, peaks_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    MultiplexIsotopicPeakPattern pattern2(c, peaks_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
}

MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, peaks_per_peptide_min, peaks_per_peptide_max, missing_peaks, intensity_cutoff, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling);
std::vector<MultiplexFilterResult> filter_results = filtering.filter();

MultiplexClustering* nullPointer = nullptr;
MultiplexClustering* ptr;

START_SECTION(MultiplexClustering(const PeakMap& exp_profile, const PeakMap& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical, double rt_minimum))
    MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical, rt_minimum);
    std::vector<std::map<int,GridBasedCluster> > cluster_results = clustering.cluster(filter_results);
    ptr = new MultiplexClustering(exp, exp_picked, boundaries_exp_s, rt_typical, rt_minimum);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical, rt_minimum);

START_SECTION(std::vector<std::map<int GridBasedCluster> > cluster(const std::vector<MultiplexFilterResult>& filter_results))
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
