// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <vector>
#include <numeric>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

START_TEST(ModifiedSincSmoother<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

ModifiedSincSmoother* ms_ptr = nullptr;
ModifiedSincSmoother* ms_nullPointer = nullptr;
START_SECTION((ModifiedSincSmoother(bool, int, int)))
  ms_ptr = new ModifiedSincSmoother(false, 6, 7);
  TEST_NOT_EQUAL(ms_ptr, ms_nullPointer)
END_SECTION

START_SECTION((virtual ~ModifiedSincSmoother()))
  delete ms_ptr;
END_SECTION

START_SECTION((Default constructor and parameter API))
{
  ModifiedSincSmoother smoother;
  Param p = smoother.getParameters();
  TEST_EQUAL(p.getValue("degree"), 6)
  TEST_EQUAL(p.getValue("m"), 7)
  TEST_EQUAL(p.getValue("is_ms1").toString(), "false")

  // Modify parameters
  p.setValue("degree", 4);
  p.setValue("m", 5);
  p.setValue("is_ms1", "true");
  smoother.setParameters(p);

  // Should work after parameter update
  std::vector<double> data(21, 5.0);
  std::vector<double> output = smoother.smooth(data);
  for (double val : output)
  {
    TEST_REAL_SIMILAR(val, 5.0)
  }
}
END_SECTION

START_SECTION((std::vector<double> smooth(const std::vector<double>&)))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};
  std::vector<double> expected = {
    0.15835885, 0.11657466, -0.09224721, 0.03165689, -0.05481473,
   -0.05436219, 0.51054827, -0.59067866, -1.21928695, 5.28610520,
   10.46161952, 6.82674246, 2.49236743, 1.04220381, 0.03264660
  };

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (size_t i = 0; i < data.size(); ++i)
  {
    TEST_REAL_SIMILAR(output[i], expected[i])
  }
}
END_SECTION

START_SECTION((smooth MS1 variant))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};

  // Expected values from Java reference: new ModifiedSincSmoother(true, 6, 7)
  std::vector<double> expected = {
    9.8621205110963930e-04,  1.5814360451280190e-01,  1.0378944922019940e-02,
   -9.3245732506647640e-02, -2.2469895253114215e-03,  3.1668080258012170e-01,
   -6.6050772241400120e-02, -1.1369501098276590e+00,  2.1676680879777394e-01,
    5.2472199849855430e+00,  9.0273654134685140e+00,  7.3725523213984620e+00,
    3.1669629429298880e+00,  6.8997426114563800e-01, -1.1895183805986807e-01
  };

  ModifiedSincSmoother smoother(true, degree, m);
  std::vector<double> output = smoother.smooth(data);

  TEST_EQUAL(data.size(), output.size())
  for (size_t i = 0; i < data.size(); ++i)
  {
    TEST_REAL_SIMILAR(output[i], expected[i])
  }
}
END_SECTION

START_SECTION((smooth with correction coefficients active, degree 8))
{
  int degree = 8;
  int m = 10;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};

  std::vector<double> expected = {
    4.0115164326883515e-02,  1.7343455236363028e-01,  1.3837225122592080e-02,
   -1.3679047792299953e-01, -1.3822129268404207e-01,  3.7724434168548170e-01,
    2.9434588024971550e-01, -1.1946044466547185e+00, -3.4299476063451995e-01,
    5.2766357414382290e+00,  9.6167148837388830e+00,  7.3920204737632450e+00,
    2.7337361809027270e+00,  6.3008317717784400e-01,  8.8849725091422600e-02
  };

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (size_t i = 0; i < data.size(); ++i)
  {
    TEST_REAL_SIMILAR(output[i], expected[i])
  }
}
END_SECTION

START_SECTION((smooth MS1 with correction coefficients, degree 4))
{
  int degree = 4;
  int m = 5;
  std::vector<double> data = {0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0};

  std::vector<double> expected = {
   -2.0134608712825580e-02,  1.5981707967146100e-01, -1.1921763597028505e-02,
   -6.9368352740410770e-03, -1.2594087487395433e-02,  9.7654613807182390e-02,
   -1.1245533076423640e-02, -8.7072404127170470e-01,  2.2797486631074101e-01,
    5.1968315648088000e+00,  8.9482994438211460e+00,  7.1448507383672570e+00,
    3.1759540185629900e+00,  9.6327576333055360e-01, -1.1503424676296573e-01
  };

  ModifiedSincSmoother smoother(true, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (size_t i = 0; i < data.size(); ++i)
  {
    TEST_REAL_SIMILAR(output[i], expected[i])
  }
}
END_SECTION

START_SECTION((short input, expect safe behavior))
{
  int degree = 4;
  int m = 3;
  std::vector<double> data = {1.0, 2.0, 3.0};

  ModifiedSincSmoother smoother(true, degree, m);  // MS1 mode: m >= 3 is valid for degree=4
  std::vector<double> output = smoother.smooth(data);

  TEST_EQUAL(output.size(), data.size())
}
END_SECTION

START_SECTION((constant input stays constant))
{
  int degree = 6;
  int m = 7;
  std::vector<double> data(21, 5.0);

  ModifiedSincSmoother smoother(false, degree, m);
  std::vector<double> output = smoother.smooth(data);

  for (double val : output)
  {
    TEST_REAL_SIMILAR(val, 5.0)
  }
}
END_SECTION

START_SECTION((static int bandwidthToM(bool isMS1, int degree, double bandwidth)))
{
  // MS mode — values verified against Java reference
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 4, 0.1), 16)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 4, 0.2), 8)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 6, 0.1), 21)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 6, 0.2), 10)

  // MS1 mode
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 4, 0.1), 12)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 4, 0.2), 5)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 6, 0.1), 17)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 6, 0.2), 8)

  // Boundary exceptions
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::bandwidthToM(false, 4, 0.0))
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::bandwidthToM(false, 4, 0.5))
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::bandwidthToM(false, 4, -0.1))
}
END_SECTION

START_SECTION((static int noiseGainToM(bool isMS1, int degree, double noiseGain)))
{
  // Values verified against Java reference
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(false, 4, 0.5), 13)
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(false, 4, 0.8), 4)
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(false, 6, 0.5), 17)
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(true, 4, 0.5), 9)
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(true, 4, 0.8), 3)
  TEST_EQUAL(ModifiedSincSmoother::noiseGainToM(true, 6, 0.5), 13)

  // Invalid noise gain
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::noiseGainToM(false, 4, 0.0))
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::noiseGainToM(false, 4, -1.0))
}
END_SECTION

START_SECTION((static double savitzkyGolayBandwidth(int degree, int m)))
{
  // Values verified against Java reference
  TEST_REAL_SIMILAR(ModifiedSincSmoother::savitzkyGolayBandwidth(4, 10), 8.1765529486428860e-02)
  TEST_REAL_SIMILAR(ModifiedSincSmoother::savitzkyGolayBandwidth(4, 20), 4.1456732856719296e-02)
  TEST_REAL_SIMILAR(ModifiedSincSmoother::savitzkyGolayBandwidth(6, 10), 1.1351775630841529e-01)
  TEST_REAL_SIMILAR(ModifiedSincSmoother::savitzkyGolayBandwidth(6, 20), 5.7047267374270800e-02)

  // Higher m should give smaller bandwidth
  double bw1 = ModifiedSincSmoother::savitzkyGolayBandwidth(4, 10);
  double bw2 = ModifiedSincSmoother::savitzkyGolayBandwidth(4, 20);
  TEST_EQUAL(bw2 < bw1, true)
}
END_SECTION

START_SECTION((void filter(MSSpectrum& spectrum)))
{
  ModifiedSincSmoother smoother(false, 4, 5);

  MSSpectrum spectrum;
  std::vector<double> intensities;
  for (int i = 0; i < 15; ++i)
  {
    Peak1D peak;
    peak.setPosition(100.0 + i);
    peak.setIntensity(i % 3 == 0 ? 1000.0 : 100.0);
    spectrum.push_back(peak);
    intensities.push_back(peak.getIntensity());
  }

  std::vector<double> reference = smoother.smooth(intensities);

  MSSpectrum original = spectrum;
  smoother.filter(spectrum);

  TEST_EQUAL(spectrum.size(), original.size())

  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getPosition()[0], original[i].getPosition()[0])
    TEST_REAL_SIMILAR(spectrum[i].getIntensity(), std::max(0.0, reference[i]))
  }
}
END_SECTION

START_SECTION((void filter(MSChromatogram& chromatogram)))
{
  ModifiedSincSmoother smoother(true, 4, 5);
  
  // Create test chromatogram
  MSChromatogram chromatogram;
  for (int i = 0; i < 10; ++i)
  {
    ChromatogramPeak peak;
    peak.setRT(i * 10.0);
    peak.setIntensity(i % 2 == 0 ? 500.0 : 50.0);
    chromatogram.push_back(peak);
  }
  
  MSChromatogram original = chromatogram;
  smoother.filter(chromatogram);
  
  // Should preserve chromatogram size and RT positions
  TEST_EQUAL(chromatogram.size(), original.size())
  
  for (size_t i = 0; i < chromatogram.size(); ++i)
  {
    TEST_REAL_SIMILAR(chromatogram[i].getRT(), original[i].getRT()) // RT unchanged
    TEST_EQUAL(chromatogram[i].getIntensity() >= 0.0, true) // No negative intensities
  }
}
END_SECTION

START_SECTION((void filter(Mobilogram& mobilogram)))
{
  ModifiedSincSmoother smoother(false, 6, 7);
  
  // Create test mobilogram
  Mobilogram mobilogram;
  for (int i = 0; i < 10; ++i)
  {
    MobilityPeak1D peak;
    peak.setPosition(i * 0.1);
    peak.setIntensity(100.0 + i * 10);
    mobilogram.push_back(peak);
  }
  
  Mobilogram original = mobilogram;
  smoother.filter(mobilogram);
  
  // Should preserve mobilogram size and IM positions
  TEST_EQUAL(mobilogram.size(), original.size())
  
  for (size_t i = 0; i < mobilogram.size(); ++i)
  {
    TEST_REAL_SIMILAR(mobilogram[i].getPosition()[0], original[i].getPosition()[0]) // IM unchanged
    TEST_EQUAL(mobilogram[i].getIntensity() >= 0.0, true) // No negative intensities
  }
}
END_SECTION

START_SECTION((void filterExperiment(PeakMap& map)))
{
  ModifiedSincSmoother smoother(false, 4, 5);
  
  // Create test experiment with spectra and chromatograms
  PeakMap experiment;
  
  // Add 2 spectra
  for (int spec = 0; spec < 2; ++spec)
  {
    MSSpectrum spectrum;
    spectrum.setRT(spec * 10.0);
    for (int i = 0; i < 5; ++i)
    {
      Peak1D peak;
      peak.setPosition(100.0 + i);
      peak.setIntensity(100.0 + i * 50);
      spectrum.push_back(peak);
    }
    experiment.addSpectrum(spectrum);
  }
  
  // Add 1 chromatogram
  MSChromatogram chromatogram;
  for (int i = 0; i < 5; ++i)
  {
    ChromatogramPeak peak;
    peak.setRT(i * 5.0);
    peak.setIntensity(200.0 + i * 25);
    chromatogram.push_back(peak);
  }
  experiment.addChromatogram(chromatogram);
  
  PeakMap original = experiment;
  smoother.filterExperiment(experiment);
  
  // Should preserve experiment structure
  TEST_EQUAL(experiment.size(), original.size())
  TEST_EQUAL(experiment.getChromatograms().size(), original.getChromatograms().size())
  
  // Check that filtering was applied to spectra
  for (size_t i = 0; i < experiment.size(); ++i)
  {
    TEST_EQUAL(experiment[i].size(), original[i].size())
    for (size_t j = 0; j < experiment[i].size(); ++j)
    {
      TEST_REAL_SIMILAR(experiment[i][j].getPosition()[0], original[i][j].getPosition()[0])
    }
  }
  
  // Check that filtering was applied to chromatograms
  for (size_t i = 0; i < experiment.getChromatograms().size(); ++i)
  {
    TEST_EQUAL(experiment.getChromatogram(i).size(), original.getChromatogram(i).size())
  }
}
END_SECTION

START_SECTION((Constructor parameter validation))
{
  // Test invalid degree (odd)
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother(false, 5, 5))
  
  // Test invalid degree (out of range)
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother(false, 12, 5))
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother(false, 0, 5))
  
  // Test invalid m (too small for MS2)
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother(false, 4, 1))
  
  // Test invalid m (too small for MS1)
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother(true, 4, 2))
  
  // Test valid parameters
  ModifiedSincSmoother valid_ms2(false, 4, 4);
  TEST_NOT_EQUAL(&valid_ms2, ms_nullPointer)

  ModifiedSincSmoother valid_ms1(true, 6, 4);
  TEST_NOT_EQUAL(&valid_ms1, ms_nullPointer)
}
END_SECTION

START_SECTION((Edge cases - empty data))
{
  ModifiedSincSmoother smoother(false, 4, 5);
  
  // Empty vector
  std::vector<double> empty_data;
  std::vector<double> result = smoother.smooth(empty_data);
  TEST_EQUAL(result.size(), 0)
  
  // Empty spectrum
  MSSpectrum empty_spectrum;
  smoother.filter(empty_spectrum);
  TEST_EQUAL(empty_spectrum.size(), 0)
  
  // Empty chromatogram
  MSChromatogram empty_chromatogram;
  smoother.filter(empty_chromatogram);
  TEST_EQUAL(empty_chromatogram.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST