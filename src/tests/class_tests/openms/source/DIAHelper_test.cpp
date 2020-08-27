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
// $Maintainer: Witold Wolski, Hannes Roest $
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#ifdef USE_BOOST_UNIT_TEST

// include boost unit test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
// macros for boost
#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5)
#define TEST_REAL_SIMILAR(val1, val2) \
  BOOST_CHECK ( boost::test_tools::check_is_close(val1, val2, EPS_05 ));
#define TEST_EQUAL(val1, val2) BOOST_CHECK_EQUAL(val1, val2);
#define END_SECTION
#define START_TEST(var1, var2)
#define END_TEST

#else
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#endif

#include <iterator>
#include <iomanip>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

using namespace std;
using namespace OpenMS;
//using namespace OpenMS::OpenSWATH;

///////////////////////////

START_TEST(DIAHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(testIntegrateWindows_test)
{
  OpenSwath::SpectrumPtr spec(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);

  mass->data.push_back(100.);
  mass->data.push_back(101.);
  mass->data.push_back(102.);
  mass->data.push_back(103.);
  mass->data.push_back(104.);
  mass->data.push_back(105.);
  mass->data.push_back(106.);

  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);
  intensity->data.push_back(1.);

  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  ion_mobility->data.push_back(1.);
  ion_mobility->data.push_back(2.);
  ion_mobility->data.push_back(3.);
  ion_mobility->data.push_back(4.);
  ion_mobility->data.push_back(5.);
  ion_mobility->data.push_back(6.);
  ion_mobility->data.push_back(7.);
  ion_mobility->description = "Ion Mobility";

  spec->setMZArray( mass);
  spec->setIntensityArray( intensity);
  spec->getDataArrays().push_back( ion_mobility );

  double mz, intens;
  DIAHelpers::integrateWindow(spec, 101., 103., mz, intens);
  // std::cout << "mz : " << mz << " int : " << intens << std::endl;
  TEST_REAL_SIMILAR (mz, 101.5);
  TEST_REAL_SIMILAR (intens, 2)

  std::vector<double> windows, intInt, intMz;
  windows.push_back(101.);
  windows.push_back(103.);
  windows.push_back(105.);
  DIAHelpers::integrateWindows(spec, windows, 2, intInt, intMz);
  TEST_EQUAL (intInt.size(), intMz.size() )
  TEST_EQUAL (intInt.size(), 3)
  TEST_REAL_SIMILAR (intInt[0], 2)
  TEST_REAL_SIMILAR (intMz[0], 100.5);

  // std::cout << "print Int" << std::endl;
  // std::copy(intInt.begin(), intInt.end(),
  //     std::ostream_iterator<double>(std::cout, " "));
  // std::cout << std::endl << "print mz" << intMz.size() << std::endl;
  // std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
  // std::copy(intMz.begin(), intMz.end(),
  //     std::ostream_iterator<double>(std::cout, " "));

  double im(0.0), im_intens(0.0);
  DIAHelpers::integrateDriftSpectrum(spec, 101., 109., im, im_intens, 2, 5);
  TEST_REAL_SIMILAR (im, 3.5);
  TEST_REAL_SIMILAR (im_intens, 4)

  double im2(0.0), im_intens2(0.0);
  DIAHelpers::integrateDriftSpectrum(spec, 101., 103., im2, im_intens2, 2, 5);
  TEST_REAL_SIMILAR (im2, 2.5);
  TEST_REAL_SIMILAR (im_intens2, 2)
}
END_SECTION

START_SECTION([EXTRA] void adjustExtractionWindow(double& right, double& left, const double& mz_extract_window, const bool& mz_extraction_ppm))
{
  // test absolute
  {
    double left(500.0), right(500.0);
    OpenMS::DIAHelpers::adjustExtractionWindow(right, left, 0.5, false);
    TEST_REAL_SIMILAR(left, 500 - 0.25);
    TEST_REAL_SIMILAR(right, 500 + 0.25);
  }
  // test ppm
  {
    double left(500.0), right(500.0);
    OpenMS::DIAHelpers::adjustExtractionWindow(right, left, 10, true);
    TEST_REAL_SIMILAR(left, 500 - 500 * 5 /1e6);
    TEST_REAL_SIMILAR(right, 500 + 500 * 5 /1e6);
  }
}
END_SECTION

START_SECTION([EXTRA] getBYSeries_test)
{
  TheoreticalSpectrumGenerator generator;
  Param p;
  p.setValue("add_metainfo", "true",
      "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
  generator.setParameters(p);

  String sequence = "SYVAWDR";
  std::vector<double> bseries, yseries;
  OpenMS::AASequence a = OpenMS::AASequence::fromString(sequence);
  OpenMS::DIAHelpers::getBYSeries(a, bseries, yseries, &generator);
  bseries.clear();
  OpenMS::DIAHelpers::getTheorMasses(a, bseries, &generator);

}
END_SECTION

#if 0
START_SECTION([EXTRA] getAveragineIsotopeDistribution_test)
{

  std::vector<std::pair<double, double> > tmp;
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(100., tmp);
  TEST_EQUAL(tmp.size() == 4, true);

  double mass1[] = { 100, 101.00048, 102.00096, 103.00144 };
  double int1[] =
      { 0.9496341, 0.0473560, 0.0029034, 0.0001064 };

  double * mm = &mass1[0];
  double * ii = &int1[0];
  for (unsigned int i = 0; i < tmp.size(); ++i, ++mm, ++ii) {

    std::cout << "mass :" << std::setprecision(10) << tmp[i].first
        << "intensity :" << tmp[i].second << std::endl;
    TEST_REAL_SIMILAR(tmp[i].first, *mm);
    TEST_REAL_SIMILAR(tmp[i].second, *ii);
  }

  tmp.clear();
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(30., tmp);
  double mass2[] = { 30, 31.0005, 32.001, 33.0014 };
  double int2[] = { 0.987254, 0.012721, 2.41038e-05, 2.28364e-08 };
  mm = &mass2[0];
  ii = &int2[0];
  for (unsigned int i = 0; i < tmp.size(); ++i, ++mm, ++ii) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
    std::cout << "mass :" << std::setprecision(10) << tmp[i].first
        << "intensity :" << tmp[i].second << std::endl;
    std::cout << i << "dm" <<  *mm - tmp[i].first << " di " << *ii - tmp[i].second << std::endl;
    TEST_REAL_SIMILAR(tmp[i].first, *mm)
    TEST_REAL_SIMILAR(tmp[i].second, *ii)
  }

  tmp.clear();
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(110., tmp);
  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }

  tmp.clear();
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(120., tmp);
  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }

  tmp.clear();
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(300., tmp);
  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }

  tmp.clear();
  OpenMS::DIAHelpers::getAveragineIsotopeDistribution(500., tmp);
  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }

}
END_SECTION

START_SECTION([EXTRA] simulateSpectrumFromAASequence_test)
{
  TheoreticalSpectrumGenerator generator;
  Param p;
  p.setValue("add_metainfo", "false",
             "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
  p.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
  generator.setParameters(p);

  String sequence = "SYVAWDR";
  OpenMS::AASequence a = OpenMS::AASequence::fromString(sequence);
  std::vector<double> masses1;
  std::vector<std::pair<double, double> > tmp, out;
  OpenMS::DIAHelpers::simulateSpectrumFromAASequence(a, masses1, tmp, &generator);

  std::copy(masses1.begin(), masses1.end(),
      std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }
  OpenMS::DIAHelpers::modifyMassesByCharge(tmp, out, 2.);
  OpenMS::DIAHelpers::addPreisotopeWeights(masses1, tmp);
  std::cout << "preisotope weights added" << std::endl;

  for (unsigned int i = 0; i < tmp.size(); ++i) {
    std::cout << "mass :" << tmp[i].first << "intensity :" << tmp[i].second
        << std::endl;
  }

}
END_SECTION

START_SECTION([EXTRA] addIsotopesToSpec_test)
{
  std::vector<std::pair<double, double> > tmp_, out;
  tmp_.push_back(std::make_pair(100., 100.));
  tmp_.push_back(std::make_pair(200., 300.));
  tmp_.push_back(std::make_pair(300., 200.));

  OpenMS::DIAHelpers::addIsotopes2Spec(tmp_, out);
  std::cout << "addIsotopesToSpec_test" << std::endl;
  for (unsigned int i = 0; i < out.size(); ++i) {
    std::cout << out[i].first << " " << out[i].second << std::endl;
  }

}
END_SECTION
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
