// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Oliver Alka$
// $Authors: Oliver Alka$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PrecursorCorrection, "$Id$")

/////////////////////////////////////////////////////////////

// Prepare dummy data
MSExperiment exp;
vector<Precursor> v_precursor_1, v_precursor_2;
Precursor precursor_1, precursor_2;
MSSpectrum ms1_spectrum_1, ms1_spectrum_2, ms2_spectrum_1, ms2_spectrum_2;
vector<MSSpectrum> v_spectra;

// precursor
precursor_1.setIntensity(240.0f);
precursor_1.setMZ(509.9999);
v_precursor_1.push_back(precursor_1);

precursor_2.setIntensity(230.0f);
precursor_2.setMZ(610.0001);
v_precursor_2.push_back(precursor_2);

// peaks
Peak1D p1;
p1.setIntensity(200.0f);
p1.setMZ(509.9994);

Peak1D p2;
p2.setIntensity(250.0f);
p2.setMZ(510.0000);

Peak1D p3;
p3.setIntensity(150.0f);
p3.setMZ(510.0001);

Peak1D p4;
p4.setIntensity(250.0f);
p4.setMZ(609.9998);

Peak1D p5;
p5.setIntensity(200.0f);
p5.setMZ(610.0000);

Peak1D p6;
p6.setIntensity(180.0f);
p6.setMZ(610.0005);

vector<Peak1D> peaks_1{p1,p2,p3};
vector<Peak1D> peaks_2{p4,p5,p6};
vector<Peak1D> empty_peaks{};

// ms1
ms1_spectrum_1.insert(ms1_spectrum_1.begin(), peaks_1.begin(), peaks_1.end());
ms1_spectrum_1.setMSLevel(1);
ms1_spectrum_1.setNativeID(1);
ms1_spectrum_1.setRT(100.0);
ms1_spectrum_2.insert(ms1_spectrum_2.begin(), peaks_2.begin(), peaks_2.end());
ms1_spectrum_2.setMSLevel(1);
ms1_spectrum_2.setNativeID(3);
ms1_spectrum_2.setRT(180.0);

// ms2
ms2_spectrum_1.insert(ms2_spectrum_1.begin(), empty_peaks.begin(), empty_peaks.end());
ms2_spectrum_1.setMSLevel(2);
ms2_spectrum_1.setNativeID(2);
ms2_spectrum_1.setRT(100.1);
ms2_spectrum_2.insert(ms2_spectrum_2.begin(), empty_peaks.begin(), empty_peaks.end());
ms2_spectrum_2.setMSLevel(2);
ms2_spectrum_2.setNativeID(4);
ms2_spectrum_2.setRT(180.1);

// ms2 precursor information
ms2_spectrum_1.setPrecursors(v_precursor_1);
ms2_spectrum_2.setPrecursors(v_precursor_2);

v_spectra.push_back(ms1_spectrum_1);
v_spectra.push_back(ms2_spectrum_1);
v_spectra.push_back(ms1_spectrum_2);
v_spectra.push_back(ms2_spectrum_2);

// MSExperiment
exp.setSpectra(v_spectra);
exp.sortSpectra();

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PrecursorCorrection* ptr = nullptr;
PrecursorCorrection* null_ptr = nullptr;
START_SECTION(PrecursorCorrection())
{
	ptr = new PrecursorCorrection();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PrecursorCorrection())
{
	delete ptr;
}
END_SECTION

START_SECTION((static void getPrecursors(const MSExperiment &exp, std::vector< Precursor > &precursors, std::vector< double > &precursors_rt, std::vector< Size > precursor_scan_index)))
{
  MSExperiment getP_exp = exp;
  vector<Precursor> precursor;
  vector<double> rt;
  vector<Size> index;
  PrecursorCorrection::getPrecursors(getP_exp, precursor, rt, index);

  TEST_EQUAL(precursor.size(), 2);
  TEST_EQUAL(rt.size(), 2);
  TEST_EQUAL(index.size(), 0);
  TEST_REAL_SIMILAR(precursor[0].getMZ(), 509.9999)
  TEST_REAL_SIMILAR(precursor[0].getIntensity(), 240.0);
  TEST_REAL_SIMILAR(rt[0], 100.1);
}
END_SECTION

FuzzyStringComparator fsc;
fsc.setAcceptableAbsolute(1e-8);

START_SECTION((static void writeHist(const String &out_csv, const std::vector< double > &deltaMZs, const std::vector< double > &mzs, const std::vector< double > &rts)))
{
  MSExperiment write_exp = exp;

  String csv_tmp;
  NEW_TMP_FILE(csv_tmp);
  vector<double> dmz;
  vector<double> mz;
  vector<double> rt;

  PrecursorCorrection::correctToNearestMS1Peak(write_exp, 1, true, dmz, mz, rt);
  PrecursorCorrection::writeHist(csv_tmp, dmz, mz,rt);

  TEST_EQUAL(fsc.compareFiles(csv_tmp, OPENMS_GET_TEST_DATA_PATH("PrecursorCorrection_out.csv")),true);
}
END_SECTION

START_SECTION((static std::set<Size> correctToNearestMS1Peak(MSExperiment &exp, double mz_tolerance, bool ppm, std::vector< double > &deltaMZs, std::vector< double > &mzs, std::vector< double > &rts)))
{
  // test with 1 ppm (1)
  MSExperiment nearest_exp_1 = exp;
  vector<double> dmz_1;
  vector<double> mz_1;
  vector<double> rt_1;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  PrecursorCorrection::correctToNearestMS1Peak(nearest_exp_1, 1, true, dmz_1, mz_1, rt_1);

  TEST_REAL_SIMILAR(dmz_1[0], 0.0001)
  TEST_REAL_SIMILAR(dmz_1[1], -0.0001)

  // test with 5 ppm (2)
  MSExperiment nearest_exp_2 = exp;
  vector<double> dmz_2;
  vector<double> mz_2;
  vector<double> rt_2;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  PrecursorCorrection::correctToNearestMS1Peak(nearest_exp_2, 5, true, dmz_2, mz_2, rt_2);

  TEST_REAL_SIMILAR(dmz_2[0], 0.0001)
  TEST_REAL_SIMILAR(dmz_2[1], -0.0001)
}
END_SECTION

START_SECTION((static std::set<Size> correctToHighestIntensityMS1Peak(MSExperiment &exp, double mz_tolerance, std::vector< double > &deltaMZs, std::vector< double > &mzs, std::vector< double > &rts)))
{
  // test with 0.0001 Da (1)
  MSExperiment highest_exp_1 = exp;
  vector<double> dmz_1;
  vector<double> mz_1;
  vector<double> rt_1;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 610.0000
  PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_1, 0.0001, dmz_1, mz_1, rt_1);

  TEST_REAL_SIMILAR(dmz_1[0], 0.0001)
  TEST_REAL_SIMILAR(dmz_1[1], -0.0001)

  // test with 0.0005 Da (2)
  MSExperiment highest_exp_2 = exp;
  vector<double> dmz_2;
  vector<double> mz_2;
  vector<double> rt_2;

  // corrected precursor_1: 510.0000
  // corrected precursor_2: 609.9998
  PrecursorCorrection::correctToHighestIntensityMS1Peak(highest_exp_2, 0.0005, dmz_2, mz_2, rt_2);

  TEST_REAL_SIMILAR(dmz_2[0], 0.0001)
  TEST_REAL_SIMILAR(dmz_2[1], -0.0003)
}
END_SECTION

// featureMap

START_SECTION((static std::set<Size> correctToNearestFeature(const FeatureMap &features, MSExperiment &exp, double rt_tolerance_s=0.0, double mz_tolerance=0.0, bool ppm=true, bool believe_charge=false, bool keep_original=false, bool all_matching_features=false, int max_trace=2, int debug_level=0)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



