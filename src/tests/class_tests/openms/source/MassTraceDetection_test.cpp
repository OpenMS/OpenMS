// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Tristan Aretz, Manuel Zschaebitz$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassTraceDetection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassTraceDetection* ptr = nullptr;
MassTraceDetection* null_ptr = nullptr;
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

START_SECTION((void updateIterativeWeightedMeanMZ(const double, const double , double&, double&, double&)))
{
  double centroid_mz(150.22), centroid_int(25000000);
  double new_mz1(150.34), new_int1(23043030);
  double new_mz2(150.11), new_int2(1932392);

  std::vector<double> mzs, ints;
  mzs.push_back(centroid_mz);
  mzs.push_back(new_mz1);
  mzs.push_back(new_mz2);
  ints.push_back(centroid_int);
  ints.push_back(new_int1);
  ints.push_back(new_int2);

  double total_weight1(centroid_int + new_int1);
  double total_weight2(centroid_int + new_int1 + new_int2);

  double wmean1((centroid_mz * centroid_int + new_mz1 * new_int1)/total_weight1);
  double wmean2((centroid_mz * centroid_int + new_mz1 * new_int1 + new_mz2 * new_int2)/total_weight2);

  double prev_nominator(centroid_mz * centroid_int);
  double prev_denomominator(centroid_int);

  test_mtd.updateIterativeWeightedMeanMZ(new_mz1, new_int1, centroid_mz, prev_nominator, prev_denomominator);

  TEST_REAL_SIMILAR(centroid_mz, wmean1);

  test_mtd.updateIterativeWeightedMeanMZ(new_mz2, new_int2, centroid_mz, prev_nominator, prev_denomominator);

  TEST_REAL_SIMILAR(centroid_mz, wmean2);

}
END_SECTION

// load a mzML file for testing the algorithm
// do it with threads so it will be faster
#ifdef _OPENMP
  omp_set_num_threads(8);
#endif
PeakMap input1;
PeakMap input2;
PeakMap input3;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input1);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_testinput1.mzML"),input2);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_testinput2.mzML"),input3);
//end thread
#ifdef _OPENMP
  omp_set_num_threads(1);
#endif


Size exp_mt_lengths[3] = {86, 31, 16};
double exp_mt_rts[3] = {341.063314463158, 339.314891947562, 350.698987241276};
double exp_mt_mzs[3] = {437.26675, 438.27241, 439.27594};
double exp_mt_ints[3] = {3381.72226139326, 664.763828332733, 109.490108620676};

std::vector<MassTrace> output_mt;

Param p_mtd = MassTraceDetection().getDefaults();
p_mtd.setValue("min_trace_length", 3.0);

  START_SECTION((void run(const PeakMap &, std::vector< MassTrace > &)))
  {
    output_mt.clear();
    // with default parameters, only 2 of 3 traces will be found
    test_mtd.run(input1, output_mt);
    TEST_EQUAL(output_mt.size(), 2);
    output_mt.clear();


    // if min_trace_length is set to 3 seconds, another mass trace is detected
    test_mtd.setParameters(p_mtd);
    test_mtd.run(input1, output_mt);
    TEST_EQUAL(output_mt.size(), 3);

    for (Size i = 0; i < output_mt.size(); ++i)
    {
      TEST_EQUAL(output_mt[i].getSize(), exp_mt_lengths[i]);
      TEST_REAL_SIMILAR(output_mt[i].getCentroidRT(), exp_mt_rts[i]);
      TEST_REAL_SIMILAR(output_mt[i].getCentroidMZ(), exp_mt_mzs[i]);
      TEST_REAL_SIMILAR(output_mt[i].computePeakArea(), exp_mt_ints[i]);
    }

    // Regression test for bug #1633
    // Test by adding MS2 spectra to the input
    
    PeakMap input_new;
    MSSpectrum s;
    s.setMSLevel(2);
    // {
    //   Peak1D p;
    //   p.setMZ( 500 );
    //   p.setIntensity( 6000 );
    //   s.push_back(p);
    // }

    // add a few additional MS2 spectra in front
    for (Size i = 0; i < input1.size(); ++i)
    {
      input_new.addSpectrum(s);
    }
    // now add the "real" spectra at the end
    for (Size i = 0; i < input1.size(); ++i)
    {
      input_new.addSpectrum(input1[i]);
    }
    output_mt.clear();
    test_mtd.run(input_new, output_mt);
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

// never implemented, but could be fun
// START_SECTION((void filterByPeakWidth(std::vector< MassTrace > &, std::vector< MassTrace > &)))
// {
//   //  test_mtd.filterByPeakWidth(output_mt, filt);

//   //  TEST_EQUAL(output_mt.size(), filt.size());

//    for (Size i = 0; i < output_mt.size(); ++i)
//    {
//        TEST_EQUAL(output_mt[i].getFWHMScansNum(), filt[i].getFWHMScansNum());
//    }
// }
// END_SECTION

PeakMap::ConstAreaIterator mt_it1 = input1.areaBeginConst(335.0, 385.0, 437.1, 437.4);
PeakMap::ConstAreaIterator mt_it2 = input1.areaBeginConst(335.0, 385.0, 438.2, 438.4);
PeakMap::ConstAreaIterator mt_it3 = input1.areaBeginConst(335.0, 385.0, 439.2, 439.4);

std::vector<MassTrace> found_mtraces;

PeakMap::ConstAreaIterator mt_end = input1.areaEndConst();

  START_SECTION((void run(PeakMap::ConstAreaIterator &begin, PeakMap::ConstAreaIterator &end, std::vector< MassTrace > &found_masstraces)))
  {

//     Daten zum testen 
// Anzahl an chrom apices: 5514
//  Anzahl an Traces: 234
//  Einzelne Peaks zum ueberpruefen RT|MZ|IN|scan_id|peak_id 

// Peak index: 0 1587.41|827.93|1.59927e+07|244|3
// Peak index: 2757 1349.09|822.429|220230|37|7
// Peak index: 5513 1329.79|991.592|76292.3|23|23
//  Einzelne Traces zum ueberpruefen RT|MZ|SD|Size|Label 

// Peak index: 0 1604.4|827.93|0.000205495|32|T1
// Peak index: 117 1495.72|844.565|0.000771067|15|T118
// Peak index: 233 1138.87|850.494|0.0168011|5|T234

// Daten zum testen 
// Anzahl an chrom apices: 5514
//  Anzahl an Traces: 234
//  Einzelne Peaks zum ueberpruefen RT|MZ|IN|scan_id|peak_id 

// Peak index: 0 1587.41|827.93|1.59927e+07|244|3
// Peak index: 2757 1349.09|822.429|220230|37|7
// Peak index: 5513 1329.79|991.592|76292.3|23|23
//  Einzelne Traces zum ueberpruefen RT|MZ|SD|Size|Label 

// Peak index: 0 1604.4|827.93|0.000205495|32|T1
// Peak index: 117 1557.22|804.446|0.000905521|15|T118
// Peak index: 233 1138.87|850.494|0.0168011|5|T234
    found_mtraces.clear();

    test_mtd.run(mt_it1, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[0]);

    found_mtraces.clear();

    test_mtd.run(mt_it2, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[1]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[1]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[1]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[1]);

    found_mtraces.clear();

    test_mtd.run(mt_it3, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[2]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[2]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[2]);

    found_mtraces.clear();


    //running test1 again with threads
    #ifdef _OPENMP
        omp_set_num_threads(8);
    #endif
    test_mtd.run(mt_it1, mt_end, found_mtraces);
    TEST_EQUAL(found_mtraces.size(), 1);
    TEST_EQUAL(found_mtraces[0].getSize(), exp_mt_lengths[0]);

    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidRT(), exp_mt_rts[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].getCentroidMZ(), exp_mt_mzs[0]);
    TEST_REAL_SIMILAR(found_mtraces[0].computePeakArea(), exp_mt_ints[0]);

    found_mtraces.clear();
  }
  END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



