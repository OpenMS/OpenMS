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


#include <fstream>
#include <omp.h>
#include <sstream>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
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



START_SECTION((void updateIterativeWeightedMeanMZ(const double, const double , double&, double&, double&)))
{
  MassTraceDetection test_mtd;
  double centroid_mz(150.22), centroid_int(25000000);
  double mzs[2] = {150.34, 150.11}; 
  double ints[2] = {23043030, 1932392};
  
  double wmean1((centroid_mz * centroid_int + mzs[0] * ints[0]) / (centroid_int + ints[0]));
  double wmean2((centroid_mz * centroid_int + mzs[0] * ints[0] + mzs[1] * ints[1])/ (centroid_int + ints[0] + ints[1]));

  double prev_nominator(centroid_mz * centroid_int);
  double prev_denomominator(centroid_int);

  test_mtd.updateIterativeWeightedMeanMZ(mzs[0], ints[0], centroid_mz, prev_nominator, prev_denomominator);

  TEST_REAL_SIMILAR(centroid_mz, wmean1);

  test_mtd.updateIterativeWeightedMeanMZ(mzs[1], ints[1], centroid_mz, prev_nominator, prev_denomominator);

  TEST_REAL_SIMILAR(centroid_mz, wmean2);
}
END_SECTION

///later
START_SECTION(void updateWeightedSDEstimateRobust(const PeakType&, const double&, double&, double& ))
{
  NOT_TESTABLE;
}
END_SECTION



START_SECTION((void run(const PeakMap &, std::vector< MassTrace > &)))
{
  MassTraceDetection test_mtd;
  std::vector<MassTrace> output_mtd;

  //input1_mtd is for testing the basic functionality
  PeakMap input1_mtd;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_input1.mzML"),input1_mtd);

  //control values
  Size exp_mt_lengths[3]{86, 31, 16};
  double exp_mt_rts[3]{341.063314463158, 339.314891947562, 350.698987241276};
  double exp_mt_mzs[3]{437.26675, 438.27241, 439.27594};
  double exp_mt_ints[3]{3381.72226139326, 664.763828332733, 109.490108620676};


  // with default parameters, only 2 of 3 traces will be found
  output_mtd.clear();
  test_mtd.run(input1_mtd, output_mtd);
  TEST_EQUAL(output_mtd.size(), 2);


  // if min_trace_length is set to 3 seconds, another mass trace is detected
  Param p_mtd = MassTraceDetection().getDefaults();
  p_mtd.setValue("min_trace_length", 3.0);
  test_mtd.setParameters(p_mtd);

  output_mtd.clear();
  test_mtd.run(input1_mtd, output_mtd);
  TEST_EQUAL(output_mtd.size(), 3);

  std::sort(output_mtd.begin(), output_mtd.end(),
    [&output_mtd](const MassTrace & a, const MassTrace & b) -> bool
    { 
    return a.getCentroidMZ() < b.getCentroidMZ();
    });
  
  for (Size i = 0; i < output_mtd.size(); ++i)
  {
    TEST_EQUAL(output_mtd[i].getSize(), exp_mt_lengths[i]);
    TEST_REAL_SIMILAR(output_mtd[i].getCentroidRT(), exp_mt_rts[i]);
    TEST_REAL_SIMILAR(output_mtd[i].getCentroidMZ(), exp_mt_mzs[i]);
    TEST_REAL_SIMILAR(output_mtd[i].computePeakArea(), exp_mt_ints[i]);
  }


  // Regression test for bug #1633
  // Test by adding MS2 spectra to the input
  PeakMap input_new;
  MSSpectrum s;
  s.setMSLevel(2);
  {
    Peak1D p;
    p.setMZ( 500 );
    p.setIntensity( 6000 );
    s.push_back(p);
  }

  // add a few additional MS2 spectra in front
  for (Size i = 0; i < input1_mtd.size(); ++i)
  {
    input_new.addSpectrum(s);
  }
  // now add the "real" spectra at the end
  for (Size i = 0; i < input1_mtd.size(); ++i)
  {
    input_new.addSpectrum(input1_mtd[i]);
  }

  output_mtd.clear();
  test_mtd.run(input_new, output_mtd);
  TEST_EQUAL(output_mtd.size(), 3);

  std::sort(output_mtd.begin(), output_mtd.end(),
    [&output_mtd](const MassTrace & a, const MassTrace & b) -> bool
    { 
    return a.getCentroidMZ() < b.getCentroidMZ();
    });

  for (Size i = 0; i < output_mtd.size(); ++i)
  {
    TEST_EQUAL(output_mtd[i].getSize(), exp_mt_lengths[i]);
    TEST_REAL_SIMILAR(output_mtd[i].getCentroidRT(), exp_mt_rts[i]);
    TEST_REAL_SIMILAR(output_mtd[i].getCentroidMZ(), exp_mt_mzs[i]);
    TEST_REAL_SIMILAR(output_mtd[i].computePeakArea(), exp_mt_ints[i]);
  }
}
END_SECTION


START_SECTION((void run(PeakMap::ConstAreaIterator &begin, PeakMap::ConstAreaIterator &end, std::vector< MassTrace > &found_masstraces)))
{
  MassTraceDetection test_mtd;
  std::vector<MassTrace> output_mtd;
  std::string ground_truth = OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_test_output1.tsv");

  /// load data with threads, so it will a bit faster
  #ifdef _OPENMP
    omp_set_num_threads(8);
  #endif

    //input2 is for testing correctness of dense data points
    PeakMap input2;
    MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz"),input2);

  #ifdef _OPENMP
    omp_set_num_threads(1);
  #endif

  PeakMap::ConstAreaIterator cai_begin = input2.areaBeginConst(2000, 4000, 500, 750, 1);
  PeakMap::ConstAreaIterator cai_end = input2.areaEndConst();


  //first run: single thread
  output_mtd.clear();
  test_mtd.run(cai_begin, cai_end, output_mtd);
  
  //has to be sorted, because of parallelization the order is different
  std::sort(output_mtd.begin(), output_mtd.end(), 
    [&output_mtd] (const MassTrace& a, const MassTrace& b) -> bool
    {
      return a.getCentroidRT() < b.getCentroidRT();
    });
  
  std::string file_single = OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_test_output2.tsv");
  NEW_TMP_FILE(file_single);
  std::ofstream os_single;

  os_single.open(file_single);
  os_single << "CentroidRT\tCentroidMZ\tCentroidSD\tTraceLength\tSize\n";
  for (const MassTrace& i : output_mtd)
  {
    os_single 
      << i.getCentroidRT() << "\t"
      << i.getCentroidMZ() << "\t"
      << i.getCentroidSD() << "\t"
      // << i.computePeakArea() << "\t" //<-- takes a bit longer, can be tested with MassTraceDetection_test_output2
      << i.getTraceLength() << "\t"
      << i.getSize() << "\n";
  }
  os_single.close();

  TEST_FILE_SIMILAR(file_single, ground_truth);

  //second run: parallel
  #ifdef _OPENMP
    omp_set_num_threads(8);
  #endif

    //where const?
    cai_begin = input2.areaBeginConst(2000, 4000, 500, 750, 1);
    cai_end = input2.areaEndConst();

    output_mtd.clear();
    test_mtd.run(cai_begin, cai_end, output_mtd);

    //has to be sorted, because of parallelization the order is different
    std::sort(output_mtd.begin(), output_mtd.end(),
      [&output_mtd](const MassTrace & a, const MassTrace & b) -> bool
      { 
      return a.getCentroidRT() < b.getCentroidRT();
      });

    std::string file_parallel = OPENMS_GET_TEST_DATA_PATH("MassTraceDetection_test_output1.tsv");
    NEW_TMP_FILE(file_parallel);
    std::ofstream os_parallel;

    os_parallel.open(file_parallel);
    os_parallel << "CentroidRT\tCentroidMZ\tCentroidSD\tTraceLength\tSize\n";
    for (const MassTrace& i : output_mtd)
    {
      os_parallel
      << i.getCentroidRT() << "\t"
      << i.getCentroidMZ() << "\t"
      << i.getCentroidSD() << "\t"
      // << i.computePeakArea() << "\t"
      << i.getTraceLength() << "\t"
      << i.getSize() << "\n";
    }
    os_parallel.close();

    TEST_FILE_SIMILAR(file_parallel, ground_truth);

  #ifdef _OPENMP
    omp_set_num_threads(1);
  #endif
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST