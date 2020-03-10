// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumAddition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static OpenSwath::SpectrumPtr addUpSpectra(std::vector< OpenSwath::SpectrumPtr > all_spectra, double sampling_rate, double filter_zeros)) )
{
  OpenSwath::SpectrumPtr spec1(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr mass1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity1(new OpenSwath::BinaryDataArray);

  OpenSwath::SpectrumPtr spec2(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr mass2(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity2(new OpenSwath::BinaryDataArray);

  // Intensity
  static const double arr1[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  };
  std::vector<double> intensity1_ (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  intensity1->data = intensity1_;
  static const double arr2[] = {
    1, 3, 5, 7, 9, 11, 9, 7, 5, 3, 1
  };
  std::vector<double> intensity2_ (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  intensity2->data = intensity2_;

  // Mass
  static const double arr3[] = {
    100, 101.5, 101.9, 102.0, 102.1, 102.11, 102.2, 102.25, 102.3, 102.4, 102.45
  };
  std::vector<double> mass1_ (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
  mass1->data = mass1_;
  static const double arr4[] = {
    100, 101.6, 101.95, 102.0, 102.05, 102.1, 102.12, 102.15, 102.2, 102.25, 102.30
  };
  std::vector<double> mass2_ (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );
  mass2->data = mass2_;


  spec1->setMZArray( mass1);
  spec1->setIntensityArray(intensity1);

  spec2->setMZArray( mass2);
  spec2->setIntensityArray(intensity2);

  std::vector<OpenSwath::SpectrumPtr> all_spectra;

  OpenSwath::SpectrumPtr empty_result = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(empty_result->getMZArray()->data.size(), 0);

  all_spectra.clear();
  OpenSwath::SpectrumPtr a(new OpenSwath::Spectrum());
  OpenSwath::SpectrumPtr b(new OpenSwath::Spectrum());
  all_spectra.push_back(a);
  all_spectra.push_back(b);
  OpenSwath::SpectrumPtr empty2 = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(empty2->getMZArray()->data.size(), 0);

  all_spectra.clear();
  all_spectra.push_back(spec1);
  all_spectra.push_back(spec2);
  OpenSwath::SpectrumPtr result = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(result->getMZArray()->data.size(), 25);

  OpenSwath::SpectrumPtr result_filtered = SpectrumAddition::addUpSpectra(all_spectra, 0.1, true);
  TEST_EQUAL(result_filtered->getMZArray()->data.size(), 9);
  TEST_REAL_SIMILAR(result_filtered->getMZArray()->data[0], 100.0);
  TEST_REAL_SIMILAR(result_filtered->getIntensityArray()->data[0], 2);
  TEST_REAL_SIMILAR(result_filtered->getMZArray()->data[3], 101.9);
  TEST_REAL_SIMILAR(result_filtered->getIntensityArray()->data[3], 3 + 5/2.0); // 3 @ 101.9 and 5 @ 101.95

  std::cout << " result size " << result->getMZArray()->data.size() << " and result m/z" << std::endl;

  std::copy(result_filtered->getMZArray()->data.begin(), result_filtered->getMZArray()->data.end(),
      std::ostream_iterator<double>(std::cout, " "));

  std::cout << std::endl << "and result intensity " << std::endl;

  std::copy(result_filtered->getIntensityArray()->data.begin(), result_filtered->getIntensityArray()->data.end(),
      std::ostream_iterator<double>(std::cout, " "));
}
END_SECTION

START_SECTION((static OpenMS::MSSpectrum addUpSpectra(std::vector< OpenMS::Spectrum<> all_spectra, double sampling_rate, bool filter_zeros) ))
{
  // Intensity
  static const double arr1[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  };
  std::vector<double> intensity1_ (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    1, 3, 5, 7, 9, 11, 9, 7, 5, 3, 1
  };
  std::vector<double> intensity2_ (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  // Mass
  static const double arr3[] = {
    100, 101.5, 101.9, 102.0, 102.1, 102.11, 102.2, 102.25, 102.3, 102.4, 102.45
  };
  std::vector<double> mass1_ (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
  static const double arr4[] = {
    100, 101.6, 101.95, 102.0, 102.05, 102.1, 102.12, 102.15, 102.2, 102.25, 102.30
  };
  std::vector<double> mass2_ (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );

  OpenMS::MSSpectrum s1;
  for (Size k = 0; k < mass1_.size(); k++)
  {
    s1.push_back(Peak1D(mass1_[k], intensity1_[k]));
  }

  OpenMS::MSSpectrum s2;
  for (Size k = 0; k < mass2_.size(); k++)
  {
    s2.push_back(Peak1D(mass2_[k], intensity2_[k]));
  }

  std::vector<MSSpectrum> all_spectra;
  MSSpectrum empty_result = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(empty_result.empty(), true);

  all_spectra.clear();
  all_spectra.push_back( MSSpectrum() );
  all_spectra.push_back( MSSpectrum() );
  std::cout << " to do here " << std::endl;
  MSSpectrum empty2 = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(empty2.size(), 0);

  all_spectra.clear();
  all_spectra.push_back(s1);
  all_spectra.push_back(s2);
  MSSpectrum result = SpectrumAddition::addUpSpectra(all_spectra, 0.1, false);
  TEST_EQUAL(result.size(), 25);

  MSSpectrum result_filtered = SpectrumAddition::addUpSpectra(all_spectra, 0.1, true);
  TEST_EQUAL(result_filtered.size(), 9);
  TEST_REAL_SIMILAR(result_filtered[0].getMZ(), 100.0);
  TEST_REAL_SIMILAR(result_filtered[0].getIntensity(), 2);
  TEST_REAL_SIMILAR(result_filtered[3].getMZ(), 101.9);
  TEST_REAL_SIMILAR(result_filtered[3].getIntensity(), 3 + 5/2.0); // 3 @ 101.9 and 5 @ 101.95

  // automatic spacing should be the min distance found in the data in each
  // spectrum individually, i.e. it should not decrease the resolution
  result_filtered = SpectrumAddition::addUpSpectra(all_spectra, 0.01, true);
  MSSpectrum result_filtered_auto = SpectrumAddition::addUpSpectra(all_spectra, -1, true);
  // this has some numerical stability issues
  // TEST_EQUAL(result_filtered, result_filtered_auto)

  TEST_EQUAL(result_filtered.size(), 16);
  TEST_REAL_SIMILAR(result_filtered[0].getMZ(), 100.0);
  TEST_REAL_SIMILAR(result_filtered[0].getIntensity(), 2);
  TEST_REAL_SIMILAR(result_filtered[3].getMZ(), 101.9);
  TEST_REAL_SIMILAR(result_filtered[3].getIntensity(), 3);

  TEST_EQUAL(result_filtered_auto.size(), 28);
  TEST_REAL_SIMILAR(result_filtered[0].getMZ(), 100.0);
  TEST_REAL_SIMILAR(result_filtered[0].getIntensity(), 2);
  TEST_REAL_SIMILAR(result_filtered[3].getMZ(), 101.9);
  TEST_REAL_SIMILAR(result_filtered[3].getIntensity(), 3);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



