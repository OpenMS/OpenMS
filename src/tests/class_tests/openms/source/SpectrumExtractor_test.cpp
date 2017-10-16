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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumExtractor.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumExtractor, "$Id$")

/////////////////////////////////////////////////////////////
// Raw spectrum data acquired in DDA mode (i.e., product ion full spectrum scan)
// measured on a QTRAP 5500 corresponding to C-Aconitate
// taken from E. coli grown on glucose M9 during steady-state
// for flux analysis.

MSSpectrum spectrum;
spectrum.resize(73);
MSSpectrum::Iterator it = spectrum.begin();

it->setMZ(61.92);
it->setIntensity(6705.41660838088f);

it++->setMZ(68.88);
it->setIntensity(1676.35415209522f);

it++->setMZ(71.4);
it->setIntensity(1676.35415209522f);

it++->setMZ(79.56);
it->setIntensity(1676.35415209522f);

it++->setMZ(84.6);
it->setIntensity(3352.70830419044f);

it++->setMZ(84.72);
it->setIntensity(5029.06245628566f);

it++->setMZ(84.84);
it->setIntensity(8381.7707604761f);

it++->setMZ(84.96);
it->setIntensity(53643.332867047f);

it++->setMZ(85.08);
it->setIntensity(51966.9787149518f);

it++->setMZ(85.2);
it->setIntensity(6705.41660838088f);

it++->setMZ(85.32);
it->setIntensity(8381.7707604761f);

it++->setMZ(85.44);
it->setIntensity(1676.35415209522f);

it++->setMZ(85.68);
it->setIntensity(11734.4790646665f);

it++->setMZ(85.8);
it->setIntensity(25145.3122814283f);

it++->setMZ(85.92);
it->setIntensity(68730.520235904f);

it++->setMZ(86.04);
it->setIntensity(112315.72819038f);

it++->setMZ(86.16);
it->setIntensity(6705.41660838088f);

it++->setMZ(86.28);
it->setIntensity(6705.41660838088f);

it++->setMZ(86.4);
it->setIntensity(3352.70830419044f);

it++->setMZ(87.72);
it->setIntensity(1676.35415209522f);

it++->setMZ(87.96);
it->setIntensity(1676.35415209522f);

it++->setMZ(88.08);
it->setIntensity(1676.35415209522f);

it++->setMZ(90.36);
it->setIntensity(3352.70830419044f);

it++->setMZ(94.44);
it->setIntensity(1676.35415209522f);

it++->setMZ(99.84);
it->setIntensity(1676.35415209522f);

it++->setMZ(100.8);
it->setIntensity(1676.35415209522f);

it++->setMZ(101.04);
it->setIntensity(5029.06245628566f);

it++->setMZ(101.88);
it->setIntensity(3352.70830419044f);

it++->setMZ(102);
it->setIntensity(3352.70830419044f);

it++->setMZ(102.96);
it->setIntensity(3352.70830419044f);

it++->setMZ(110.16);
it->setIntensity(1676.35415209522f);

it++->setMZ(110.88);
it->setIntensity(5029.06245628566f);

it++->setMZ(111);
it->setIntensity(3352.70830419044f);

it++->setMZ(111.12);
it->setIntensity(5029.06245628566f);

it++->setMZ(111.24);
it->setIntensity(3352.70830419044f);

it++->setMZ(111.84);
it->setIntensity(5029.06245628566f);

it++->setMZ(111.96);
it->setIntensity(18439.8956730474f);

it++->setMZ(112.08);
it->setIntensity(20116.2498251426f);

it++->setMZ(112.2);
it->setIntensity(5029.06245628566f);

it++->setMZ(112.32);
it->setIntensity(1676.35415209522f);

it++->setMZ(112.44);
it->setIntensity(1676.35415209522f);

it++->setMZ(112.56);
it->setIntensity(3352.70830419044f);

it++->setMZ(112.68);
it->setIntensity(3352.70830419044f);

it++->setMZ(114);
it->setIntensity(3352.70830419044f);

it++->setMZ(128.16);
it->setIntensity(6705.41660838088f);

it++->setMZ(128.4);
it->setIntensity(1676.35415209522f);

it++->setMZ(128.88);
it->setIntensity(3352.70830419044f);

it++->setMZ(129);
it->setIntensity(3352.70830419044f);

it++->setMZ(129.12);
it->setIntensity(6705.41660838088f);

it++->setMZ(129.84);
it->setIntensity(5029.06245628566f);

it++->setMZ(129.96);
it->setIntensity(10058.1249125713f);

it++->setMZ(130.08);
it->setIntensity(31850.7288898092f);

it++->setMZ(130.2);
it->setIntensity(10058.1249125713f);

it++->setMZ(130.32);
it->setIntensity(1676.35415209522f);

it++->setMZ(130.44);
it->setIntensity(1676.35415209522f);

it++->setMZ(130.56);
it->setIntensity(3352.70830419044f);

it++->setMZ(132.12);
it->setIntensity(1676.35415209522f);

it++->setMZ(138);
it->setIntensity(1676.35415209522f);

it++->setMZ(139.08);
it->setIntensity(1676.35415209522f);

it++->setMZ(140.16);
it->setIntensity(3352.70830419044f);

it++->setMZ(144.12);
it->setIntensity(1676.35415209522f);

it++->setMZ(146.04);
it->setIntensity(3352.70830419044f);

it++->setMZ(146.16);
it->setIntensity(1676.35415209522f);

it++->setMZ(156);
it->setIntensity(1676.35415209522f);

it++->setMZ(156.12);
it->setIntensity(5029.06245628566f);

it++->setMZ(156.36);
it->setIntensity(1676.35415209522f);

it++->setMZ(173.76);
it->setIntensity(1676.35415209522f);

it++->setMZ(174);
it->setIntensity(1676.35415209522f);

it++->setMZ(174.12);
it->setIntensity(6705.41660838088f);

it++->setMZ(174.24);
it->setIntensity(11734.4790646665f);

it++->setMZ(174.36);
it->setIntensity(6705.41660838088f);

it++->setMZ(174.6);
it->setIntensity(1676.35415209522f);

it++->setMZ(175.08);
it->setIntensity(1676.35415209522f);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumExtractor* ptr = 0;
SpectrumExtractor* null_ptr = 0;
START_SECTION(SpectrumExtractor())
{
  ptr = new SpectrumExtractor();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~SpectrumExtractor())
{
  delete ptr;
}
END_SECTION

ptr = new SpectrumExtractor();

START_SECTION(getParameters())
{
  Param params = ptr->getParameters();
  TEST_EQUAL(params.getValue("rt_window"), 30)
  TEST_EQUAL(params.getValue("min_score"), 0.7)
  TEST_EQUAL(params.getValue("min_forward_match"), 0.9)
  TEST_EQUAL(params.getValue("min_reverse_match"), 0.9)
  TEST_EQUAL(params.getValue("mz_tolerance"), 0.1)
  TEST_EQUAL(params.getValue("mz_tolerance_units"), "Da")
  TEST_EQUAL(params.getValue("sgolay_frame_length"), 15)
  TEST_EQUAL(params.getValue("sgolay_polynomial_order"), 3)
  TEST_EQUAL(params.getValue("gauss_width"), 0.2)
  TEST_EQUAL(params.getValue("use_gauss"), "true")
  TEST_EQUAL(params.getValue("signal_to_noise"), 1.0)
}
END_SECTION

START_SECTION(setRTWindow())
{
  TEST_EQUAL(ptr->getRTWindow(), 30)
  ptr->setRTWindow(50);
  TEST_EQUAL(ptr->getRTWindow(), 50)
}
END_SECTION

START_SECTION(setMinScore())
{
  TEST_EQUAL(ptr->getMinScore(), 0.7)
  ptr->setMinScore(0.5);
  TEST_EQUAL(ptr->getMinScore(), 0.5)
}
END_SECTION

START_SECTION(setMinForwardMatch())
{
  TEST_EQUAL(ptr->getMinForwardMatch(), 0.9)
  ptr->setMinForwardMatch(0.5);
  TEST_EQUAL(ptr->getMinForwardMatch(), 0.5)
}
END_SECTION

START_SECTION(setMinReverseMatch())
{
  TEST_EQUAL(ptr->getMinReverseMatch(), 0.9)
  ptr->setMinReverseMatch(0.5);
  TEST_EQUAL(ptr->getMinReverseMatch(), 0.5)
}
END_SECTION

START_SECTION(setMZTolerance())
{
  TEST_EQUAL(ptr->getMZTolerance(), 0.1)
  ptr->setMZTolerance(0.5);
  TEST_EQUAL(ptr->getMZTolerance(), 0.5)
}
END_SECTION

START_SECTION(setMZToleranceUnits())
{
  TEST_EQUAL(ptr->getMZToleranceUnits(), "Da")
  TEST_NOT_EQUAL(ptr->getMZToleranceUnits(), "ppm")
  ptr->setMZToleranceUnits("ppm");
  TEST_EQUAL(ptr->getMZToleranceUnits(), "ppm")
}
END_SECTION

START_SECTION(setSGolayFrameLength())
{
  TEST_EQUAL(ptr->getSGolayFrameLength(), 15)
  ptr->setSGolayFrameLength(7);
  TEST_EQUAL(ptr->getSGolayFrameLength(), 7)
}
END_SECTION

START_SECTION(setSGolayPolynomialOrder())
{
  TEST_EQUAL(ptr->getSGolayPolynomialOrder(), 3)
  ptr->setSGolayPolynomialOrder(2);
  TEST_EQUAL(ptr->getSGolayPolynomialOrder(), 2)
}
END_SECTION

START_SECTION(setGaussWidth())
{
  TEST_EQUAL(ptr->getGaussWidth(), 0.2)
  ptr->setGaussWidth(0.5);
  TEST_EQUAL(ptr->getGaussWidth(), 0.5)
}
END_SECTION

START_SECTION(setUseGauss())
{
  TEST_EQUAL(ptr->getUseGauss(), true)
  ptr->setUseGauss(false);
  TEST_EQUAL(ptr->getUseGauss(), false)
}
END_SECTION

START_SECTION(setSignalToNoise())
{
  TEST_EQUAL(ptr->getSignalToNoise(), 1.0)
  ptr->setSignalToNoise(0.6);
  TEST_EQUAL(ptr->getSignalToNoise(), 0.6)
}
END_SECTION

START_SECTION(getParameters().getValue("rt_window"))
{
  TEST_EQUAL(ptr->getParameters().getValue("rt_window"), 50)
}
END_SECTION

START_SECTION(getParameters().getDescription("rt_window"))
{
  TEST_EQUAL(ptr->getParameters().getDescription("rt_window"), "Retention time window in seconds.")
}
END_SECTION

START_SECTION(pickSpectrum())
{
  // TODO replace and improve this test
  // at the moment the output is a text file to use
  // as plotly input data (x and y of a scatter plot)
  MSSpectrum picked;
  spectrum.sortByPosition();
  ptr->setGaussWidth(0.2);
  ptr->setUseGauss(true);
  ptr->pickSpectrum(spectrum, picked);
  TEST_NOT_EQUAL(spectrum.size(), picked.size())
  ofstream outfile;
  outfile.open("plot_output.txt", ios::out | ios::trunc);
  if (outfile.is_open()) {
    outfile << "Fill plotly data with this:" << std::endl << "x: [";
    for (UInt i=0; i<picked.size(); ++i) {
      outfile << picked[i].getMZ() << ", ";
    }
    outfile << "]," << std::endl << "y: [";
    for (UInt i=0; i<picked.size(); ++i) {
      outfile << picked[i].getIntensity() << ", ";
    }
    outfile << "],";
    outfile.close();
  }
  else {
    cout << "Unable to open file to write the spectrum";
  }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
