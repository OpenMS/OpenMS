// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexClustering.h>

using namespace OpenMS;

START_TEST(MultiplexFiltering, "$Id$")

// read data
MSExperiment<Peak1D> exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexClustering.mzML"), exp);

// pick data
PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms1_only", DataValue("true"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);
std::vector<PeakPickerHiRes::PeakBoundary> boundaries;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;
MSExperiment<Peak1D> exp_picked;
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
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double rt_typical = 90;
double rt_minimum = 5;
bool debug = false;

// construct list of peak patterns
std::vector<PeakPattern> patterns;
std::vector<double> shifts1;
shifts1.push_back(0);
shifts1.push_back(8.0443702794);
std::vector<double> shifts2;
shifts2.push_back(0);
shifts2.push_back(2*8.0443702794);
for (int c = charge_max; c >= charge_min; --c)
{
    PeakPattern pattern1(c, peaks_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    PeakPattern pattern2(c, peaks_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
} 

/*MultiplexFiltering filtering(exp, exp_picked, boundaries_exp_s, patterns, peaks_per_peptide_min, peaks_per_peptide_max, missing_peaks, intensity_cutoff, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, debug);
// The above line throws the following exception:
// 243: Error: Caught unexpected OpenMS exception of type 'Exception' thrown in line 83 of file '/home/lars/Code/git/OpenMS/src/openms/source/FILTERING/DATAREDUCTION/MultiplexFiltering.cpp' in function 'OpenMS::MultiplexFiltering::MultiplexFiltering(OpenMS::MSExperiment<OpenMS::Peak1D>, OpenMS::MSExperiment<OpenMS::Peak1D>, std::vector<std::vector<OpenMS::PeakPickerHiRes::PeakBoundary> >, std::vector<OpenMS::PeakPattern>, int, int, bool, double, double, bool, double, double, bool)' - Message: Centroided data and the corresponding list of peak boundaries do not contain same number of spectra. (72 != 72)
// But 72 == 72! Why the exception?
std::vector<FilterResult> filter_results = filtering.filter();

MultiplexClustering* nullPointer = 0;
MultiplexClustering* ptr;

START_SECTION(MultiplexClustering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, double rt_typical, double rt_minimum, bool debug))
    MultiplexClustering clustering(exp, exp_picked, boundaries_exp_s, rt_typical, rt_minimum, debug);
    std::vector<std::map<int,Cluster> > cluster_results = clustering.cluster(filter_results);
    ptr = new MultiplexClustering(exp, exp_picked, boundaries_exp_s, rt_typical, rt_minimum, debug);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION*/

END_TEST
