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
OpenSwath::SpectrumPtr imSpec(new OpenSwath::Spectrum());
{
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

    imSpec->setMZArray( mass );
    imSpec->setIntensityArray( intensity );
    imSpec->getDataArrays().push_back( ion_mobility );
}

OpenSwath::SpectrumPtr spec(new OpenSwath::Spectrum());
{
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

    spec->setMZArray( mass );
    spec->setIntensityArray( intensity );
}


START_SECTION(bool integrateWindow(const OpenSwath::SpectrumPtr& spectrum, double mz_start, double mz_end, double & mz, double & im, double & intensity, double drift_start, double drift_end, bool centroided))
{

  {
    //Test integration of empty spectrum
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    double mz(0), intens(0), im(0);

    DIAHelpers::integrateWindow(emptySpec, 101., 103., mz, im, intens, -1, -1);

    TEST_REAL_SIMILAR(mz, -1);
    TEST_REAL_SIMILAR(im, -1);
    TEST_REAL_SIMILAR(intens, 0);
  }

  {
    // Test spectrum without ion mobility while asking for ion mobility filtering, should throw an exception
    double mz(0), intens(0), im(0);
    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, DIAHelpers::integrateWindow(spec, 101., 103., mz, im, intens, 2, 5), "Cannot integrate with drift time if no drift time is available");
  }

  {
    // Test ion mobility enhanced array with no ion mobility windows, although IM is present it should be ignored
    double mz(0), intens(0), im(0);

    DIAHelpers::integrateWindow(imSpec, 101., 103., mz, im, intens, -1, -1);
    TEST_REAL_SIMILAR (mz, 101.5);
    TEST_REAL_SIMILAR (intens, 2);
    TEST_REAL_SIMILAR (im, -1); // since no IM, this value should be -1
  }



  {
    // Test With Ion Mobility (Condition 1/2)
    double mz(0), intens(0), im(0);

    DIAHelpers::integrateWindow(imSpec, 101., 109., mz, im, intens, 2, 5);
    TEST_REAL_SIMILAR (im, 3.5);
    TEST_REAL_SIMILAR (intens, 4);
  }

  {
    // Test with Ion Mobility (Condition 2/2)
    double mz(0), intens(0), im(0);

    DIAHelpers::integrateWindow(imSpec, 101., 103., mz, im, intens, 2, 5);
    TEST_REAL_SIMILAR (im, 2.5);
    TEST_REAL_SIMILAR (intens, 2);
  }
}
END_SECTION


START_SECTION(bool integrateWindow(const std::vector<OpenSwath::SpectrumPtr>& spectra, double mz_start, double mz_end, double & mz, double & im, double & intensity, double drift_start, double drift_end, bool centroided))
{
  {
    // Test integration of empty array
    std::vector<OpenSwath::SpectrumPtr> emptySpecArr;
    double mz(0), intens(0), im(0);

    DIAHelpers::integrateWindow(emptySpecArr, 101., 103., mz, im, intens, -1, -1);
    TEST_REAL_SIMILAR(mz, -1);
    TEST_REAL_SIMILAR(im, -1);
    TEST_REAL_SIMILAR(intens, 0);
  }

  {
    //Test integration of empty spectrum
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<OpenSwath::SpectrumPtr> specArrEmptySpectrum;
    double mz(0), intens(0), im(0);

    specArrEmptySpectrum.push_back(emptySpec);
    DIAHelpers::integrateWindow(specArrEmptySpectrum, 101., 103., mz, im, intens, -1, -1);

    TEST_REAL_SIMILAR(mz, -1);
    TEST_REAL_SIMILAR(im, -1);
    TEST_REAL_SIMILAR(intens, 0);
  }

  {
    // Test ion mobility enhanced array with no ion mobility windows, although IM is present it should be ignored
    std::vector<OpenSwath::SpectrumPtr> specArr;
    double mz(0), intens(0), im(0);
    specArr.push_back(imSpec);

    DIAHelpers::integrateWindow(specArr, 101., 103., mz, im, intens, -1, -1);
    TEST_REAL_SIMILAR (mz, 101.5);
    TEST_REAL_SIMILAR (intens, 2);
    TEST_REAL_SIMILAR (im, -1); // since no IM, this value should be -1
  }

  {
    // Test With Ion Mobility (Condition 1/2)
    std::vector<OpenSwath::SpectrumPtr> specArr;
    double mz(0), intens(0), im(0);

    specArr.push_back(imSpec);

    DIAHelpers::integrateWindow(specArr, 101., 109., mz, im, intens, 2, 5);
    TEST_REAL_SIMILAR (im, 3.5);
    TEST_REAL_SIMILAR (intens, 4);
  }

  {
    // Test with Ion Mobility (Condition 2/2)
    std::vector<OpenSwath::SpectrumPtr> specArr;
    double mz(0), intens(0), im(0);

    specArr.push_back(imSpec);

    DIAHelpers::integrateWindow(specArr, 101., 103., mz, im, intens, 2, 5);
    TEST_REAL_SIMILAR (im, 2.5);
    TEST_REAL_SIMILAR (intens, 2);
  }
}
END_SECTION

START_SECTION(void integrateWindows(const OpenSwath::SpectrumPtr& spectra, const std::vector<double> & windowsCenter, double width, std::vector<double> & integratedWindowsIntensity, std::vector<double> & integratedWindowsMZ, std::vector<double> & integratedWindowsIm, double drift_start, double drift_end, bool remZero))
{
  {
    // Test empty spectrum (with non empty windows) - remove zeros
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<double> windows, intInt, intMz, intIm;

    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);

    DIAHelpers::integrateWindows(emptySpec, windows, 2, intInt, intMz, intIm, -1, -1, true);
    TEST_EQUAL (intInt.empty(), true);
    TEST_EQUAL (intIm.empty(), true);
    TEST_EQUAL (intMz.empty(), true);
  }

  {
    // Test empty spectrum (with non empty windows) - Don't remove zeros
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;

    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);

    DIAHelpers::integrateWindows(emptySpec, windows, 2, intInt, intMz, intIm, -1, -1, false);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 3)
    TEST_REAL_SIMILAR (intInt[0], 0)
    TEST_REAL_SIMILAR (intInt[1], 0)
    TEST_REAL_SIMILAR (intInt[2], 0)
    TEST_REAL_SIMILAR (intMz[0], 101.) // should be middle of window
    TEST_REAL_SIMILAR (intMz[1], 103.) // should be middle of window
    TEST_REAL_SIMILAR (intMz[2], 105.) // should be middle of window
    TEST_REAL_SIMILAR (intIm[0], -1) // should be avg. drift
    TEST_REAL_SIMILAR (intIm[1], -1) // should be avg. drift
    TEST_REAL_SIMILAR (intIm[2], -1) // should be avg. drift
  }

  {
    // Test non empty spectrum with no im
    std::vector<double> windows, intInt, intMz, intIm;
    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);

    DIAHelpers::integrateWindows(spec, windows, 2, intInt, intMz, intIm, -1, -1);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 3)
    TEST_REAL_SIMILAR (intInt[0], 2)
    TEST_REAL_SIMILAR (intMz[0], 100.5);
    TEST_REAL_SIMILAR (intIm[0], -1);

    // std::cout << "print Int" << std::endl;
    // std::copy(intInt.begin(), intInt.end(),
    //     std::ostream_iterator<double>(std::cout, " "));
    // std::cout << std::endl << "print mz" << intMz.size() << std::endl;
    // std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
    // std::copy(intMz.begin(), intMz.end(),
    //     std::ostream_iterator<double>(std::cout, " "));
  }

  {
    // Test non empty spectrum with ion mobility
    std::vector<double> windows, intInt, intMz, intIm;
    windows.push_back(102.);

    DIAHelpers::integrateWindows(imSpec, windows, 2, intInt, intMz, intIm, 2, 5);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 1)
    TEST_REAL_SIMILAR (intInt[0], 2)
    TEST_REAL_SIMILAR (intMz[0], 101.5);
    TEST_REAL_SIMILAR (intIm[0], 2.5);

    /*
    std::cout << "print Int" << std::endl;
    std::copy(intInt.begin(), intInt.end(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl << "print mz" << intMz.size() << std::endl;
    std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
    std::copy(intMz.begin(), intMz.end(),
        std::ostream_iterator<double>(std::cout, " "));
    */
  }
}
END_SECTION


START_SECTION(void integrateWindows(const std::vector<OpenSwath::SpectrumPtr>& spectrum, const std::vector<double> & windowsCenter, double width, std::vector<double> & integratedWindowsIntensity, std::vector<double> & integratedWindowsMZ, std::vector<double> & integratedWindowsIm, double drift_start, double drift_end, bool remZero))
{
  {
    // Test empty windows
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;
    specArr.push_back(emptySpec);

    TEST_EXCEPTION_WITH_MESSAGE(Exception::MissingInformation, DIAHelpers::integrateWindows(specArr, windows, 2, intInt, intMz, intIm, -1, -1), "No windows supplied!");
  }


  {
    // Test empty spectrum (with non empty windows) - remove zeros
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;
    specArr.push_back(emptySpec);

    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);


    DIAHelpers::integrateWindows(specArr, windows, 2, intInt, intMz, intIm, -1, -1, true);
    TEST_EQUAL (intInt.empty(), true);
    TEST_EQUAL (intIm.empty(), true);
    TEST_EQUAL (intMz.empty(), true);
  }

  {
    // Test empty spectrum (with non empty windows) - Don't remove zeros
    OpenSwath::SpectrumPtr emptySpec(new OpenSwath::Spectrum());
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;
    specArr.push_back(emptySpec);

    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);

    DIAHelpers::integrateWindows(specArr, windows, 2, intInt, intMz, intIm, -1, -1, false);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 3)
    TEST_REAL_SIMILAR (intInt[0], 0)
    TEST_REAL_SIMILAR (intInt[1], 0)
    TEST_REAL_SIMILAR (intInt[2], 0)
    TEST_REAL_SIMILAR (intMz[0], 101.) // should be middle of window
    TEST_REAL_SIMILAR (intMz[1], 103.) // should be middle of window
    TEST_REAL_SIMILAR (intMz[2], 105.) // should be middle of window
    TEST_REAL_SIMILAR (intIm[0], -1) // should be avg. drift
    TEST_REAL_SIMILAR (intIm[1], -1) // should be avg. drift
    TEST_REAL_SIMILAR (intIm[2], -1) // should be avg. drift
  }


  {
    // Test non empty spectrum with no im
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;
    windows.push_back(101.);
    windows.push_back(103.);
    windows.push_back(105.);

    specArr.push_back(spec);

    DIAHelpers::integrateWindows(specArr, windows, 2, intInt, intMz, intIm, -1, -1);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 3)
    TEST_REAL_SIMILAR (intInt[0], 2)
    TEST_REAL_SIMILAR (intMz[0], 100.5);
    TEST_REAL_SIMILAR (intIm[0], -1);

    // std::cout << "print Int" << std::endl;
    // std::copy(intInt.begin(), intInt.end(),
    //     std::ostream_iterator<double>(std::cout, " "));
    // std::cout << std::endl << "print mz" << intMz.size() << std::endl;
    // std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
    // std::copy(intMz.begin(), intMz.end(),
    //     std::ostream_iterator<double>(std::cout, " "));
  }

  {
    // Test non empty spectrum with ion mobility
    std::vector<OpenSwath::SpectrumPtr> specArr;
    std::vector<double> windows, intInt, intMz, intIm;
    windows.push_back(102.);

    specArr.push_back(imSpec);

    DIAHelpers::integrateWindows(specArr, windows, 2, intInt, intMz, intIm, 2, 5);
    TEST_EQUAL (intInt.size(), intMz.size() )
    TEST_EQUAL (intInt.size(), intIm.size() )
    TEST_EQUAL (intInt.size(), 1)
    TEST_REAL_SIMILAR (intInt[0], 2)
    TEST_REAL_SIMILAR (intMz[0], 101.5);
    TEST_REAL_SIMILAR (intIm[0], 2.5);

    /*
    std::cout << "print Int" << std::endl;
    std::copy(intInt.begin(), intInt.end(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl << "print mz" << intMz.size() << std::endl;
    std::cout << intMz[0] << " " << intMz[1] << " " << intMz[2] << std::endl;
    std::copy(intMz.begin(), intMz.end(),
        std::ostream_iterator<double>(std::cout, " "));
    */
  }
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

  // test case where extraction leads to a negative number, should correct this to the 0 bounds (no ppm)
  // Note since this is very unlikely with ppm, this functionality is currently not implemented in ppm
  {
    double left(500.0), right(500.0);
    OpenMS::DIAHelpers::adjustExtractionWindow(right, left, 1000, false);
    TEST_REAL_SIMILAR(left, 0)
    TEST_REAL_SIMILAR(right, 500 + 500)
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
