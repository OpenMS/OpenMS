// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>
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

  ModifiedSincSmoother smoother(true, degree, m); // MS1 variant
  std::vector<double> output = smoother.smooth(data);

  double sum_input = std::accumulate(data.begin(), data.end(), 0.0);
  double sum_output = std::accumulate(output.begin(), output.end(), 0.0);
  TEST_REAL_SIMILAR(sum_input, sum_output)

  TEST_EQUAL(data.size(), output.size())
}
END_SECTION

START_SECTION((short input, expect safe behavior))
{
  int degree = 4;
  int m = 3;
  std::vector<double> data = {1.0, 2.0, 3.0};

  ModifiedSincSmoother smoother(false, degree, m);
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
  // Test bandwidth conversion for MS2
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 4, 0.1), 30)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(false, 6, 0.1), 32)
  
  // Test bandwidth conversion for MS1
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 4, 0.1), 10)
  TEST_EQUAL(ModifiedSincSmoother::bandwidthToM(true, 6, 0.1), 12)
  
  // Test invalid bandwidth - should throw
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::bandwidthToM(false, 4, 0.0))
  TEST_EXCEPTION(Exception::InvalidParameter, ModifiedSincSmoother::bandwidthToM(false, 4, 0.5))
}
END_SECTION

START_SECTION((static int noiseGainToM(bool isMS1, int degree, double noiseGain)))
{
  // Test noise gain conversion
  int m_ms2 = ModifiedSincSmoother::noiseGainToM(false, 4, 0.8);
  int m_ms1 = ModifiedSincSmoother::noiseGainToM(true, 4, 0.8);
  
  TEST(m_ms2 > 0)
  TEST(m_ms1 > 0)
  TEST(m_ms2 != m_ms1) // Should be different for MS1 vs MS2
}
END_SECTION

START_SECTION((static double savitzkyGolayBandwidth(int degree, int m)))
{
  // Test Savitzky-Golay bandwidth calculation
  double bw = ModifiedSincSmoother::savitzkyGolayBandwidth(4, 10);
  TEST(bw > 0.0)
  TEST(bw < 0.5)
  
  // Higher m should give smaller bandwidth
  double bw2 = ModifiedSincSmoother::savitzkyGolayBandwidth(4, 20);
  TEST(bw2 < bw)
}
END_SECTION

START_SECTION((void filter(MSSpectrum& spectrum)))
{
  ModifiedSincSmoother smoother(false, 4, 5);
  
  // Create test spectrum
  MSSpectrum spectrum;
  for (int i = 0; i < 10; ++i)
  {
    Peak1D peak;
    peak.setPosition(100.0 + i);
    peak.setIntensity(i % 3 == 0 ? 1000.0 : 100.0); // Some peaks higher
    spectrum.push_back(peak);
  }
  
  MSSpectrum original = spectrum;
  smoother.filter(spectrum);
  
  // Should preserve spectrum size and positions
  TEST_EQUAL(spectrum.size(), original.size())
  
  for (size_t i = 0; i < spectrum.size(); ++i)
  {
    TEST_REAL_SIMILAR(spectrum[i].getPosition(), original[i].getPosition()) // Positions unchanged
    TEST(spectrum[i].getIntensity() >= 0.0) // No negative intensities
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
    TEST(chromatogram[i].getIntensity() >= 0.0) // No negative intensities
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
    TEST_REAL_SIMILAR(mobilogram[i].getPosition(), original[i].getPosition()) // IM unchanged
    TEST(mobilogram[i].getIntensity() >= 0.0) // No negative intensities
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
      TEST_REAL_SIMILAR(experiment[i][j].getPosition(), original[i][j].getPosition())
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
  TEST_NOT_EXCEPTION(ModifiedSincSmoother(false, 4, 4))
  TEST_NOT_EXCEPTION(ModifiedSincSmoother(true, 6, 4))
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