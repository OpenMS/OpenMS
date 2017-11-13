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
it++->setIntensity(6705.41660838088f);

it->setMZ(68.88);
it++->setIntensity(1676.35415209522f);

it->setMZ(71.4);
it++->setIntensity(1676.35415209522f);

it->setMZ(79.56);
it++->setIntensity(1676.35415209522f);

it->setMZ(84.6);
it++->setIntensity(3352.70830419044f);

it->setMZ(84.72);
it++->setIntensity(5029.06245628566f);

it->setMZ(84.84);
it++->setIntensity(8381.7707604761f);

it->setMZ(84.96);
it++->setIntensity(53643.332867047f);

it->setMZ(85.08);
it++->setIntensity(51966.9787149518f);

it->setMZ(85.2);
it++->setIntensity(6705.41660838088f);

it->setMZ(85.32);
it++->setIntensity(8381.7707604761f);

it->setMZ(85.44);
it++->setIntensity(1676.35415209522f);

it->setMZ(85.68);
it++->setIntensity(11734.4790646665f);

it->setMZ(85.8);
it++->setIntensity(25145.3122814283f);

it->setMZ(85.92);
it++->setIntensity(68730.520235904f);

it->setMZ(86.04);
it++->setIntensity(112315.72819038f);

it->setMZ(86.16);
it++->setIntensity(6705.41660838088f);

it->setMZ(86.28);
it++->setIntensity(6705.41660838088f);

it->setMZ(86.4);
it++->setIntensity(3352.70830419044f);

it->setMZ(87.72);
it++->setIntensity(1676.35415209522f);

it->setMZ(87.96);
it++->setIntensity(1676.35415209522f);

it->setMZ(88.08);
it++->setIntensity(1676.35415209522f);

it->setMZ(90.36);
it++->setIntensity(3352.70830419044f);

it->setMZ(94.44);
it++->setIntensity(1676.35415209522f);

it->setMZ(99.84);
it++->setIntensity(1676.35415209522f);

it->setMZ(100.8);
it++->setIntensity(1676.35415209522f);

it->setMZ(101.04);
it++->setIntensity(5029.06245628566f);

it->setMZ(101.88);
it++->setIntensity(3352.70830419044f);

it->setMZ(102);
it++->setIntensity(3352.70830419044f);

it->setMZ(102.96);
it++->setIntensity(3352.70830419044f);

it->setMZ(110.16);
it++->setIntensity(1676.35415209522f);

it->setMZ(110.88);
it++->setIntensity(5029.06245628566f);

it->setMZ(111);
it++->setIntensity(3352.70830419044f);

it->setMZ(111.12);
it++->setIntensity(5029.06245628566f);

it->setMZ(111.24);
it++->setIntensity(3352.70830419044f);

it->setMZ(111.84);
it++->setIntensity(5029.06245628566f);

it->setMZ(111.96);
it++->setIntensity(18439.8956730474f);

it->setMZ(112.08);
it++->setIntensity(20116.2498251426f);

it->setMZ(112.2);
it++->setIntensity(5029.06245628566f);

it->setMZ(112.32);
it++->setIntensity(1676.35415209522f);

it->setMZ(112.44);
it++->setIntensity(1676.35415209522f);

it->setMZ(112.56);
it++->setIntensity(3352.70830419044f);

it->setMZ(112.68);
it++->setIntensity(3352.70830419044f);

it->setMZ(114);
it++->setIntensity(3352.70830419044f);

it->setMZ(128.16);
it++->setIntensity(6705.41660838088f);

it->setMZ(128.4);
it++->setIntensity(1676.35415209522f);

it->setMZ(128.88);
it++->setIntensity(3352.70830419044f);

it->setMZ(129);
it++->setIntensity(3352.70830419044f);

it->setMZ(129.12);
it++->setIntensity(6705.41660838088f);

it->setMZ(129.84);
it++->setIntensity(5029.06245628566f);

it->setMZ(129.96);
it++->setIntensity(10058.1249125713f);

it->setMZ(130.08);
it++->setIntensity(31850.7288898092f);

it->setMZ(130.2);
it++->setIntensity(10058.1249125713f);

it->setMZ(130.32);
it++->setIntensity(1676.35415209522f);

it->setMZ(130.44);
it++->setIntensity(1676.35415209522f);

it->setMZ(130.56);
it++->setIntensity(3352.70830419044f);

it->setMZ(132.12);
it++->setIntensity(1676.35415209522f);

it->setMZ(138);
it++->setIntensity(1676.35415209522f);

it->setMZ(139.08);
it++->setIntensity(1676.35415209522f);

it->setMZ(140.16);
it++->setIntensity(3352.70830419044f);

it->setMZ(144.12);
it++->setIntensity(1676.35415209522f);

it->setMZ(146.04);
it++->setIntensity(3352.70830419044f);

it->setMZ(146.16);
it++->setIntensity(1676.35415209522f);

it->setMZ(156);
it++->setIntensity(1676.35415209522f);

it->setMZ(156.12);
it++->setIntensity(5029.06245628566f);

it->setMZ(156.36);
it++->setIntensity(1676.35415209522f);

it->setMZ(173.76);
it++->setIntensity(1676.35415209522f);

it->setMZ(174);
it++->setIntensity(1676.35415209522f);

it->setMZ(174.12);
it++->setIntensity(6705.41660838088f);

it->setMZ(174.24);
it++->setIntensity(11734.4790646665f);

it->setMZ(174.36);
it++->setIntensity(6705.41660838088f);

it->setMZ(174.6);
it++->setIntensity(1676.35415209522f);

it->setMZ(175.08);
it->setIntensity(1676.35415209522f);

START_SECTION(getMZ())
{
  TEST_EQUAL(spectrum[0].getMZ(), 61.92)
  TEST_EQUAL(spectrum[0].getIntensity(), 6705.41660838088f)
  TEST_EQUAL(spectrum[1].getMZ(), 68.88)
  TEST_EQUAL(spectrum[1].getIntensity(), 1676.35415209522f)
  TEST_EQUAL(spectrum[6].getMZ(), 84.84)
  TEST_EQUAL(spectrum[6].getIntensity(), 8381.7707604761f)
  TEST_EQUAL(spectrum[71].getMZ(), 174.6)
  TEST_EQUAL(spectrum[71].getIntensity(), 1676.35415209522f)
  TEST_EQUAL(spectrum[72].getMZ(), 175.08)
  TEST_EQUAL(spectrum[72].getIntensity(), 1676.35415209522f)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SpectrumExtractor* ptr = 0;
SpectrumExtractor* null_ptr = 0;
const String experiment_path = OPENMS_GET_TEST_DATA_PATH("SpectrumExtractor_13C1_spectra0to100.mzML");
const String target_list_path = OPENMS_GET_TEST_DATA_PATH("SpectrumExtractor_13CFlux_TraML.csv");
MzMLFile mzml;
PeakMap experiment;
TransitionTSVReader tsv_reader;
TargetedExperiment targeted_exp;
mzml.load(experiment_path, experiment);
tsv_reader.convertTSVToTargetedExperiment(target_list_path.c_str(), FileTypes::CSV, targeted_exp);

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
  TEST_EQUAL(params.getValue("rt_window"), 30.0)
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
  TEST_EQUAL(params.getValue("peak_height_min"), 0.0)
  TEST_EQUAL(params.getValue("peak_height_max"), 4e6)
  TEST_EQUAL(params.getValue("fwhm_threshold"), 0.0)
  TEST_EQUAL(params.getValue("tic_weight"), 1.0)
  TEST_EQUAL(params.getValue("fwhm_weight"), 1.0)
  TEST_EQUAL(params.getValue("snr_weight"), 1.0)
}
END_SECTION

START_SECTION(setRTWindow())
{
  TEST_EQUAL(ptr->getRTWindow(), 30.0)
  ptr->setRTWindow(50.0);
  TEST_EQUAL(ptr->getRTWindow(), 50.0)
}
END_SECTION

START_SECTION(setMinScore())
{
  TEST_EQUAL(ptr->getMinScore(), 0.7)
  ptr->setMinScore(2.5);
  TEST_EQUAL(ptr->getMinScore(), 2.5)
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

START_SECTION(getPeakHeightMin())
{
  TEST_EQUAL(ptr->getPeakHeightMin(), 0.0)
  ptr->setPeakHeightMin(0.6);
  TEST_EQUAL(ptr->getPeakHeightMin(), 0.6)
}
END_SECTION

START_SECTION(getPeakHeightMax())
{
  TEST_EQUAL(ptr->getPeakHeightMax(), 4e6)
  ptr->setPeakHeightMax(150000.0);
  TEST_EQUAL(ptr->getPeakHeightMax(), 150000.0)
}
END_SECTION

START_SECTION(getFWHMThreshold())
{
  TEST_EQUAL(ptr->getFWHMThreshold(), 0.0)
  ptr->setFWHMThreshold(0.23);
  TEST_EQUAL(ptr->getFWHMThreshold(), 0.23)
}
END_SECTION

START_SECTION(getParameters().getDescription("rt_window"))
{
  TEST_EQUAL(ptr->getParameters().getDescription("rt_window"), "Retention time window in seconds.")
}
END_SECTION

START_SECTION(setTICWeight())
{
  TEST_EQUAL(ptr->getTICWeight(), 1.0)
  ptr->setTICWeight(2.0);
  TEST_EQUAL(ptr->getTICWeight(), 2.0)
}
END_SECTION

START_SECTION(setFWHMWeight())
{
  TEST_EQUAL(ptr->getFWHMWeight(), 1.0)
  ptr->setFWHMWeight(2.0);
  TEST_EQUAL(ptr->getFWHMWeight(), 2.0)
}
END_SECTION

START_SECTION(setSNRWeight())
{
  TEST_EQUAL(ptr->getSNRWeight(), 1.0)
  ptr->setSNRWeight(2.0);
  TEST_EQUAL(ptr->getSNRWeight(), 2.0)
}
END_SECTION

START_SECTION(pickSpectrum())
{
  MSSpectrum picked_spectrum;
  spectrum.sortByPosition();

  ptr->setUseGauss(true);
  ptr->setGaussWidth(0.25);
  ptr->setPeakHeightMin(0.0);
  ptr->setPeakHeightMax(200000.0);
  ptr->setFWHMThreshold(0.0);
  ptr->pickSpectrum(spectrum, picked_spectrum);
  TEST_NOT_EQUAL(spectrum.size(), picked_spectrum.size())
  TEST_EQUAL(picked_spectrum.size(), 6)
  MSSpectrum::Iterator it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 86.0196)
  TEST_REAL_SIMILAR(it->getIntensity(), 116036.0)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 129.396)
  TEST_REAL_SIMILAR(it->getIntensity(), 10575.5)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 130.081)
  TEST_REAL_SIMILAR(it->getIntensity(), 31838.1)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 174.24)
  TEST_REAL_SIMILAR(it->getIntensity(), 11731.3)

  ptr->setPeakHeightMin(15000.0);
  ptr->setPeakHeightMax(110000.0);
  ptr->pickSpectrum(spectrum, picked_spectrum);
  // With the new filters on peaks' heights, less peaks get picked.
  TEST_EQUAL(picked_spectrum.size(), 3)
  it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 130.081)
  TEST_REAL_SIMILAR(it->getIntensity(), 31838.1)

  ptr->setFWHMThreshold(0.23);
  ptr->pickSpectrum(spectrum, picked_spectrum);
  // Filtering also on fwhm, even less peaks get picked.
  TEST_EQUAL(picked_spectrum.size(), 2)
  it = picked_spectrum.begin();
  TEST_REAL_SIMILAR(it->getMZ(), 85.014)
  TEST_REAL_SIMILAR(it->getIntensity(), 60754.7)
  ++it;
  TEST_REAL_SIMILAR(it->getMZ(), 112.033)
  TEST_REAL_SIMILAR(it->getIntensity(), 21941.9)
}
END_SECTION

START_SECTION(annotateSpectra())
{
  ptr->setUseGauss(true);
  ptr->setGaussWidth(0.25);
  ptr->setRTWindow(30.0);
  ptr->setMZTolerance(0.1);
  ptr->setPeakHeightMin(15000.0);
  ptr->setPeakHeightMax(110000.0);
  ptr->setFWHMThreshold(0.23);

  vector<MSSpectrum> spectra = experiment.getSpectra();
  vector<MSSpectrum> annotated_spectra;
  FeatureMap features;

  ptr->annotateSpectra(spectra, targeted_exp, annotated_spectra, features);

  TEST_NOT_EQUAL(annotated_spectra.size(), 0)
  TEST_EQUAL(annotated_spectra.size(), features.size())

  TEST_EQUAL(annotated_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_EQUAL(annotated_spectra[0].size(), 121)
  TEST_EQUAL(annotated_spectra[4].getName(), "glu-L.glu-L_m4-4")
  TEST_EQUAL(annotated_spectra[4].size(), 98)
  TEST_EQUAL(annotated_spectra[8].getName(), "asp-L.asp-L_m0-0")
  TEST_EQUAL(annotated_spectra[8].size(), 61)
  TEST_EQUAL(annotated_spectra[12].getName(), "asp-L.asp-L_m3-2")
  TEST_EQUAL(annotated_spectra[12].size(), 432)
  TEST_EQUAL(annotated_spectra[16].getName(), "glu-L.glu-L_m1-1")
  TEST_EQUAL(annotated_spectra[16].size(), 149)
  TEST_EQUAL(annotated_spectra[20].getName(), "skm.skm_m4-3")
  TEST_EQUAL(annotated_spectra[20].size(), 552)

  TEST_EQUAL(features[0].getMetaValue("transition_name"), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(features[0].getRT(), 1.337016666667)
  TEST_REAL_SIMILAR(features[0].getMZ(), 148.052001953125)
  TEST_EQUAL(features[4].getMetaValue("transition_name"), "glu-L.glu-L_m4-4")
  TEST_REAL_SIMILAR(features[4].getRT(), 2.050366666667)
  TEST_REAL_SIMILAR(features[4].getMZ(), 150.059005737305)
  TEST_EQUAL(features[8].getMetaValue("transition_name"), "asp-L.asp-L_m0-0")
  TEST_REAL_SIMILAR(features[8].getRT(), 2.162216666667)
  TEST_REAL_SIMILAR(features[8].getMZ(), 132.029998779297)
  TEST_EQUAL(features[12].getMetaValue("transition_name"), "asp-L.asp-L_m3-2")
  TEST_REAL_SIMILAR(features[12].getRT(), 2.339083333333)
  TEST_REAL_SIMILAR(features[12].getMZ(), 135.039993286133)
  TEST_EQUAL(features[16].getMetaValue("transition_name"), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(features[16].getRT(), 2.452616666667)
  TEST_REAL_SIMILAR(features[16].getMZ(), 147.04899597168)
  TEST_EQUAL(features[20].getMetaValue("transition_name"), "skm.skm_m4-3")
  TEST_REAL_SIMILAR(features[20].getRT(), 2.782566666667)
  TEST_REAL_SIMILAR(features[20].getMZ(), 177.057998657227)

  cout << "Annotated spectra from annotateSpectra():" << endl;
  for (auto s : annotated_spectra)
  {
    cout << "name: " << s.getName() << "\t peaks: " << s.size() << endl;
  }
  cout << "Features from annotateSpectra():" << endl;
  for (auto f : features)
  {
    cout << "name: " << f.getMetaValue("transition_name") << "\t RT: " << f.getRT() << "\t MZ: " << f.getMZ() << endl;
  }
}
END_SECTION

START_SECTION(scoreSpectra())
{
  ptr->setUseGauss(true);
  ptr->setGaussWidth(0.25);
  ptr->setRTWindow(30.0);
  ptr->setMZTolerance(0.1);
  ptr->setPeakHeightMin(15000.0);
  ptr->setPeakHeightMax(110000.0);
  ptr->setFWHMThreshold(0.23);
  ptr->setTICWeight(1.0);
  ptr->setFWHMWeight(1.0);
  ptr->setSNRWeight(1.0);

  vector<MSSpectrum> annotated_spectra;
  FeatureMap features;
  const vector<MSSpectrum> spectra = experiment.getSpectra();

  ptr->annotateSpectra(spectra, targeted_exp, annotated_spectra, features);
  TEST_EQUAL(annotated_spectra.size(), features.size())

  vector<MSSpectrum> picked_spectra(annotated_spectra.size());
  for (UInt i=0; i<annotated_spectra.size(); ++i)
  {
    ptr->pickSpectrum(annotated_spectra[i], picked_spectra[i]);
  }

  for (Int i=annotated_spectra.size()-1; i>=0; --i)
  {
    if (!picked_spectra[i].size())
    {
      annotated_spectra.erase(annotated_spectra.begin() + i);
      picked_spectra.erase(picked_spectra.begin() + i);
      features.erase(features.begin() + i);
    }
  }
  TEST_EQUAL(annotated_spectra.size(), features.size())
  TEST_EQUAL(picked_spectra.size(), features.size())

  vector<MSSpectrum> scored_spectra;
  ptr->scoreSpectra(annotated_spectra, picked_spectra, features, scored_spectra);

  TEST_NOT_EQUAL(scored_spectra.size(), 0)
  TEST_EQUAL(scored_spectra.size(), annotated_spectra.size())
  TEST_EQUAL(scored_spectra.size(), features.size())

  TEST_EQUAL(scored_spectra[0].getName(), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[1][0], 15.2046270370483) // score
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[2][0], 5.3508939743042)  // total tic
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[3][0], 3.96267318725586) // inverse average fwhm
  TEST_REAL_SIMILAR(scored_spectra[0].getFloatDataArrays()[4][0], 5.89106035232544) // average snr

  TEST_EQUAL(scored_spectra[4].getName(), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[1][0], 10.90163230896)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[2][0], 6.50192594528198)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[3][0], 2.14086890220642)
  TEST_REAL_SIMILAR(scored_spectra[4].getFloatDataArrays()[4][0], 2.25883746147156)

  TEST_EQUAL(scored_spectra[8].getName(), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[1][0], 13.7276296615601)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[2][0], 5.51675566061136)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[3][0], 3.46319246830875)
  TEST_REAL_SIMILAR(scored_spectra[8].getFloatDataArrays()[4][0], 4.74768113612061)

  TEST_EQUAL(scored_spectra[11].getName(), "skm.skm_m4-3")
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[1][0], 10.5747480392456)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[2][0], 6.60354130105922)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[3][0], 2.0288507938385)
  TEST_REAL_SIMILAR(scored_spectra[11].getFloatDataArrays()[4][0], 1.94235549504842)

  TEST_EQUAL(features[0].getMetaValue("transition_name"), "met-L.met-L_m0-0")
  TEST_REAL_SIMILAR(features[0].getIntensity(), 15.2046270370483)                  // score
  TEST_REAL_SIMILAR(features[0].getMetaValue("log10_total_tic"), 5.3508939743042)  // total tic
  TEST_REAL_SIMILAR(features[0].getMetaValue("inverse_avgFWHM"), 3.96267318725586) // inverse average fwhm
  TEST_REAL_SIMILAR(features[0].getMetaValue("avgSNR"), 5.89106035232544)          // average snr
  TEST_REAL_SIMILAR(features[0].getMetaValue("avgFWHM"), 0.252354895075162)        // average fwhm

  TEST_EQUAL(features[4].getMetaValue("transition_name"), "asp-L.asp-L_m1-0")
  TEST_REAL_SIMILAR(features[4].getIntensity(), 10.90163230896)
  TEST_REAL_SIMILAR(features[4].getMetaValue("log10_total_tic"), 6.50192594528198)
  TEST_REAL_SIMILAR(features[4].getMetaValue("inverse_avgFWHM"), 2.14086890220642)
  TEST_REAL_SIMILAR(features[4].getMetaValue("avgSNR"), 2.25883746147156)
  TEST_REAL_SIMILAR(features[4].getMetaValue("avgFWHM"), 0.467100044855705)

  TEST_EQUAL(features[8].getMetaValue("transition_name"), "glu-L.glu-L_m1-1")
  TEST_REAL_SIMILAR(features[8].getIntensity(), 13.7276296615601)
  TEST_REAL_SIMILAR(features[8].getMetaValue("log10_total_tic"), 5.51675566061136)
  TEST_REAL_SIMILAR(features[8].getMetaValue("inverse_avgFWHM"), 3.46319246830875)
  TEST_REAL_SIMILAR(features[8].getMetaValue("avgSNR"), 4.74768113612061)
  TEST_REAL_SIMILAR(features[8].getMetaValue("avgFWHM"), 0.288750916719437)

  TEST_EQUAL(features[11].getMetaValue("transition_name"), "skm.skm_m4-3")
  TEST_REAL_SIMILAR(features[11].getIntensity(), 10.5747480392456)
  TEST_REAL_SIMILAR(features[11].getMetaValue("log10_total_tic"), 6.60354130105922)
  TEST_REAL_SIMILAR(features[11].getMetaValue("inverse_avgFWHM"), 2.02885079241748)
  TEST_REAL_SIMILAR(features[11].getMetaValue("avgSNR"), 1.94235549504842)
  TEST_REAL_SIMILAR(features[11].getMetaValue("avgFWHM"), 0.492889868361609)

  // sort(scored_spectra.begin(), scored_spectra.end(), [](MSSpectrum a, MSSpectrum b)
  // {
  //   return a.getFloatDataArrays()[1][0] > b.getFloatDataArrays()[1][0];
  // });
  // sort(scored_spectra.begin(), scored_spectra.end(), [](MSSpectrum a, MSSpectrum b)
  // {
  //   return a.getName().compare(b.getName()) < 0;
  // });
  // cout <<  endl << "Scored spectra have been sorted by name." << endl;

  cout << endl << "Info from scored spectra:" << endl;
  for (auto s : scored_spectra)
  {
    cout << s.getName()
    << "\t score: " << s.getFloatDataArrays()[1][0]
    << "\t log10_tic: " << s.getFloatDataArrays()[2][0]
    << "\t 1/fwhm: " << s.getFloatDataArrays()[3][0]
    << "\t SNR: " << s.getFloatDataArrays()[4][0] << endl;
  }

  cout << endl << "Info from FeatureMap:" << endl;
  for (auto f : features)
  {
    cout << f.getMetaValue("transition_name")
    << "\t score: " << f.getIntensity()
    << "\t log10_tic: " << f.getMetaValue("log10_total_tic")
    << "\t 1/fwhm: " << f.getMetaValue("inverse_avgFWHM")
    << "\t SNR: " << f.getMetaValue("avgSNR")
    << "\t fwhm: " << f.getMetaValue("avgFWHM") << endl;
  }
}
END_SECTION

START_SECTION(extractSpectra())
{
  ptr->setUseGauss(true);
  ptr->setGaussWidth(0.25);
  ptr->setRTWindow(30.0);
  ptr->setMZTolerance(0.1);
  ptr->setPeakHeightMin(15000.0);
  ptr->setPeakHeightMax(110000.0);
  ptr->setFWHMThreshold(0.23);
  ptr->setTICWeight(1.0);
  ptr->setFWHMWeight(1.0);
  ptr->setSNRWeight(1.0);
  ptr->setMinScore(15.0);

  vector<MSSpectrum> extracted_spectra;
  FeatureMap extracted_features;
  ptr->extractSpectra(experiment, targeted_exp, extracted_spectra, extracted_features);

  TEST_EQUAL(extracted_spectra.size(), extracted_features.size())

  cout << endl << "Printing mapping of transition -> best spectrum:" << endl;
  for (UInt i=0; i<extracted_spectra.size(); ++i)
  {
    cout << extracted_spectra[i].getName() << "\t" << extracted_features[i].getIntensity() << endl;
  }
}
END_SECTION

START_SECTION(selectSpectra())
{
  ptr->setUseGauss(true);
  ptr->setGaussWidth(0.25);
  ptr->setRTWindow(30.0);
  ptr->setMZTolerance(0.1);
  ptr->setPeakHeightMin(15000.0);
  ptr->setPeakHeightMax(110000.0);
  ptr->setFWHMThreshold(0.23);
  ptr->setTICWeight(1.0);
  ptr->setFWHMWeight(1.0);
  ptr->setSNRWeight(1.0);
  ptr->setMinScore(15.0);

  const std::vector<MSSpectrum> spectra = experiment.getSpectra();
  std::vector<MSSpectrum> annotated;
  FeatureMap features;
  ptr->annotateSpectra(spectra, targeted_exp, annotated, features);
  std::vector<MSSpectrum> picked(annotated.size());
  for (UInt i=0; i<annotated.size(); ++i)
  {
    ptr->pickSpectrum(annotated[i], picked[i]);
  }
  for (Int i=annotated.size()-1; i>=0; --i)
  {
    if (!picked[i].size())
    {
      annotated.erase(annotated.begin() + i);
      picked.erase(picked.begin() + i);
      features.erase(features.begin() + i);
    }
  }
  std::vector<MSSpectrum> scored;
  ptr->scoreSpectra(annotated, picked, features, scored);

  cout << endl << "Scored spectra and their score:" << endl;
  for (auto s : scored)
  {
    cout << s.getName() << "\t" << s.getFloatDataArrays()[1][0] << endl;
  }

  std::vector<MSSpectrum> selected_spectra;
  FeatureMap selected_features;

  ptr->selectSpectra(scored, features, selected_spectra, selected_features);
  TEST_EQUAL(selected_spectra.size(), 2)
  TEST_EQUAL(selected_spectra.size(), selected_features.size())
  for (UInt i=0; i<selected_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(selected_spectra[i].getName(), "")
    TEST_EQUAL(selected_spectra[i].getName(), selected_features[i].getMetaValue("transition_name"))
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0], selected_features[i].getIntensity())
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0] >= ptr->getMinScore(), true)
  }

  ptr->selectSpectra(scored, selected_spectra);
  TEST_NOT_EQUAL(selected_spectra.size(), 0)
  for (UInt i=0; i<selected_spectra.size(); ++i)
  {
    TEST_NOT_EQUAL(selected_spectra[i].getName(), "")
    TEST_EQUAL(selected_spectra[i].getFloatDataArrays()[1][0] >= ptr->getMinScore(), true)
  }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
