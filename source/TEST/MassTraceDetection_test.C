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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassTraceDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassTraceDetection* ptr = 0;
MassTraceDetection* null_ptr = 0;
START_SECTION(MassTraceDetection())
{
    ptr = new MassTraceDetection();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MassTraceDetection())
{
    delete ptr;
}
END_SECTION

MassTraceDetection test_mtd;

START_SECTION((void updateIterativeWeightedMeanMZ(const DoubleReal &, const DoubleReal &, DoubleReal &, DoubleReal &, DoubleReal &)))
{
    DoubleReal centroid_mz(150.22), centroid_int(25000000);
    DoubleReal new_mz1(150.34), new_int1(23043030);
    DoubleReal new_mz2(150.11), new_int2(1932392);

    std::vector<DoubleReal> mzs, ints;
    mzs.push_back(centroid_mz);
    mzs.push_back(new_mz1);
    mzs.push_back(new_mz2);
    ints.push_back(centroid_int);
    ints.push_back(new_int1);
    ints.push_back(new_int2);

    DoubleReal total_weight1(centroid_int + new_int1);
    DoubleReal total_weight2(centroid_int + new_int1 + new_int2);

    DoubleReal wmean1((centroid_mz * centroid_int + new_mz1 * new_int1)/total_weight1);
    DoubleReal wmean2((centroid_mz * centroid_int + new_mz1 * new_int1 + new_mz2 * new_int2)/total_weight2);

    DoubleReal prev_count(centroid_mz * centroid_int);
    DoubleReal prev_denom(centroid_int);

    test_mtd.updateIterativeWeightedMeanMZ(new_mz1, new_int1, centroid_mz, prev_count, prev_denom);

    TEST_REAL_SIMILAR(centroid_mz, wmean1);

    test_mtd.updateIterativeWeightedMeanMZ(new_mz2, new_int2, centroid_mz, prev_count, prev_denom);

    TEST_REAL_SIMILAR(centroid_mz, wmean2);

}
END_SECTION

// load a mzML file for testing the algorithm
MSExperiment<Peak1D> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input);

Size exp_mt_lengths[3] = {86, 31, 16};
DoubleReal exp_mt_rts[3] = {347.778, 346.881, 346.836};
DoubleReal exp_mt_mzs[3] = {437.26675, 438.27241, 439.27594};
DoubleReal exp_mt_ints[3] = {3124.765, 631.45, 116.966};

std::vector<MassTrace> output_mt;

Param p_mtd = MassTraceDetection().getDefaults();
p_mtd.setValue("min_trace_length", 3.0);

START_SECTION((void run(const MSExperiment< Peak1D > &, std::vector< MassTrace > &)))
{
    test_mtd.run(input, output_mt);

    // with default parameters, only 2 of 3 traces will be found
    TEST_EQUAL(output_mt.size(), 2);

    // if min_trace_length is set to 3 seconds, another mass trace is detected
    test_mtd.setParameters(p_mtd);
    output_mt.clear();

    test_mtd.run(input, output_mt);

    TEST_EQUAL(output_mt.size(), 3);

    for (Size i = 0; i < output_mt.size(); ++i)
    {
        TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
        TEST_REAL_SIMILAR(output_mt[i].getCentroidMZ(), exp_mt_mzs[i]);
        TEST_REAL_SIMILAR(output_mt[i].computePeakArea(), exp_mt_ints[i]);
    }
}
END_SECTION

std::vector<MassTrace> filt;

//START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
//{
//    test_mtd.filterByPeakWidth(output_mt, filt);

//    TEST_EQUAL(output_mt.size(), filt.size());

////    for (Size i = 0; i < output_mt.size(); ++i)
////    {
////        TEST_EQUAL(output_mt[i].getFWHMScansNum(), filt[i].getFWHMScansNum());
////    }
//}
//END_SECTION

MSExperiment<Peak1D>::ConstAreaIterator mt_it1 = input.areaBeginConst(335.0, 385.0, 437.1, 437.4);
MSExperiment<Peak1D>::ConstAreaIterator mt_it2 = input.areaBeginConst(335.0, 385.0, 438.2, 438.4);
MSExperiment<Peak1D>::ConstAreaIterator mt_it3 = input.areaBeginConst(335.0, 385.0, 439.2, 439.4);

std::vector<MassTrace> found_mtraces;

MSExperiment<Peak1D>::ConstAreaIterator mt_end = input.areaEndConst();

START_SECTION((void run(MSExperiment< Peak1D >::ConstAreaIterator &begin, MSExperiment< Peak1D >::ConstAreaIterator &end, std::vector< MassTrace > &found_masstraces)))
{

    NOT_TESTABLE
//    test_mtd.run(mt_it1, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[0]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[0]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[0]);

//    found_mtraces.clear();


//    test_mtd.run(mt_it2, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[1]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[1]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[1]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[1]);

//    found_mtraces.clear();


//    test_mtd.run(mt_it3, mt_end, found_mtraces);
//    TEST_EQUAL(found_mtraces.size(), 1);
//    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[2]);
//    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[2]);
//    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[2]);

//    found_mtraces.clear();
}
END_SECTION




/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



