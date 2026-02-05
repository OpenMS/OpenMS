// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

MSExperiment createTestExperiment()
{
  MSExperiment exp;
  
  // Create 3 spectra with known RT and m/z values
  for (Size i = 0; i < 3; ++i)
  {
    MSSpectrum spec;
    spec.setRT(100.0 + i * 50.0);  // RT: 100, 150, 200
    spec.setMSLevel(1);
    
    // Add 3 peaks per spectrum
    Peak1D peak;
    peak.setMZ(500.0);
    peak.setIntensity(1000.0);
    spec.push_back(peak);
    
    peak.setMZ(550.0);
    peak.setIntensity(2000.0);
    spec.push_back(peak);
    
    peak.setMZ(600.0);
    peak.setIntensity(3000.0);
    spec.push_back(peak);
    
    exp.addSpectrum(spec);
  }
  
  exp.updateRanges();
  return exp;
}

START_TEST(MSExperiment_rasterizeRTMZ, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((void rasterizeRTMZ(float* output, Size rt_bins, Size mz_bins, CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, UInt ms_level, RasterAggregation aggregation) const))
{
  MSExperiment exp = createTestExperiment();
  
  // Test basic rasterization with SUM aggregation
  {
    Size rt_bins = 10;
    Size mz_bins = 10;
    vector<float> output(rt_bins * mz_bins, 0.0f);
    
    exp.rasterizeRTMZ(output.data(), rt_bins, mz_bins, 
                    100.0, 200.0,    // RT range
                    500.0, 600.0,    // m/z range
                    1,               // MS level
                    MSExperiment::RasterAggregation::SUM);
    
    // Total intensity across all bins should equal sum of all peak intensities
    float total_sum = 0.0f;
    for (float val : output) total_sum += val;
    // 3 spectra * (1000 + 2000 + 3000) = 18000
    TEST_REAL_SIMILAR(total_sum, 18000.0f)
  }
  
  // Test with MAX aggregation
  {
    Size rt_bins = 10;
    Size mz_bins = 10;
    vector<float> output(rt_bins * mz_bins, 0.0f);
    
    exp.rasterizeRTMZ(output.data(), rt_bins, mz_bins, 
                    100.0, 200.0,
                    500.0, 600.0,
                    1,
                    MSExperiment::RasterAggregation::MAX);
    
    // The maximum bin value should equal the largest single peak intensity (3000)
    float max_val = 0.0f;
    for (float val : output) if (val > max_val) max_val = val;
    TEST_REAL_SIMILAR(max_val, 3000.0f)
  }
  
  // Test with empty experiment
  {
    MSExperiment empty_exp;
    Size rt_bins = 5;
    Size mz_bins = 5;
    vector<float> output(rt_bins * mz_bins, 1.0f);  // Initialize with non-zero
    
    empty_exp.rasterizeRTMZ(output.data(), rt_bins, mz_bins,
                         0.0, 100.0,
                         0.0, 1000.0,
                         1,
                         MSExperiment::RasterAggregation::SUM);
    
    // Output should be all zeros
    for (float val : output)
    {
      TEST_EQUAL(val, 0.0f)
    }
  }
  
  // Test with single spectrum
  {
    MSExperiment single_exp;
    MSSpectrum spec;
    spec.setRT(50.0);
    spec.setMSLevel(1);
    Peak1D peak;
    peak.setMZ(500.0);
    peak.setIntensity(100.0);
    spec.push_back(peak);
    single_exp.addSpectrum(spec);
    
    Size rt_bins = 10;
    Size mz_bins = 10;
    vector<float> output(rt_bins * mz_bins, 0.0f);
    
    single_exp.rasterizeRTMZ(output.data(), rt_bins, mz_bins,
                          0.0, 100.0,
                          400.0, 600.0,
                          1,
                          MSExperiment::RasterAggregation::SUM);
    
    // Should have exactly one bin with intensity equal to the peak intensity (100)
    int nonzero_count = 0;
    float found_val = 0.0f;
    for (float val : output)
    {
      if (val > 0.0f)
      {
        nonzero_count++;
        found_val = val;
      }
    }
    TEST_EQUAL(nonzero_count, 1)
    TEST_REAL_SIMILAR(found_val, 100.0f)
  }
  
  // Test MS level filtering
  {
    MSExperiment mixed_exp;
    for (Size i = 0; i < 5; ++i)
    {
      MSSpectrum spec;
      spec.setRT(100.0 + i * 10.0);
      spec.setMSLevel(i % 2 == 0 ? 1 : 2);  // Alternate MS1/MS2
      Peak1D peak;
      peak.setMZ(500.0);
      peak.setIntensity(1000.0);
      spec.push_back(peak);
      mixed_exp.addSpectrum(spec);
    }
    
    Size rt_bins = 10;
    Size mz_bins = 10;
    vector<float> output_ms1(rt_bins * mz_bins, 0.0f);
    vector<float> output_ms2(rt_bins * mz_bins, 0.0f);
    
    // Rasterize MS1 only
    mixed_exp.rasterizeRTMZ(output_ms1.data(), rt_bins, mz_bins,
                         100.0, 150.0,
                         400.0, 600.0,
                         1,
                         MSExperiment::RasterAggregation::SUM);
    
    // Rasterize MS2 only
    mixed_exp.rasterizeRTMZ(output_ms2.data(), rt_bins, mz_bins,
                         100.0, 150.0,
                         400.0, 600.0,
                         2,
                         MSExperiment::RasterAggregation::SUM);
    
    // Both should have total intensities matching the number of spectra per MS level
    float sum_ms1 = 0.0f, sum_ms2 = 0.0f;
    for (Size i = 0; i < output_ms1.size(); ++i)
    {
      sum_ms1 += output_ms1[i];
      sum_ms2 += output_ms2[i];
    }
    // There are 3 MS1 spectra and 2 MS2 spectra in the RT range, each with intensity 1000
    TEST_REAL_SIMILAR(sum_ms1, 3000.0f)
    TEST_REAL_SIMILAR(sum_ms2, 2000.0f)
  }
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Exception: null pointer))
{
  MSExperiment exp = createTestExperiment();
  
  // Test null pointer exception
  TEST_EXCEPTION(Exception::NullPointer,
    exp.rasterizeRTMZ(nullptr, 10, 10, 0.0, 100.0, 0.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Exception: zero rt_bins))
{
  MSExperiment exp = createTestExperiment();
  vector<float> output(100, 0.0f);
  
  // Test zero rt_bins
  TEST_EXCEPTION(Exception::InvalidValue,
    exp.rasterizeRTMZ(output.data(), 0, 10, 0.0, 100.0, 0.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Exception: zero mz_bins))
{
  MSExperiment exp = createTestExperiment();
  vector<float> output(100, 0.0f);
  
  // Test zero mz_bins
  TEST_EXCEPTION(Exception::InvalidValue,
    exp.rasterizeRTMZ(output.data(), 10, 0, 0.0, 100.0, 0.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Exception: inverted RT range))
{
  MSExperiment exp = createTestExperiment();
  vector<float> output(100, 0.0f);
  
  // Test inverted RT range (min >= max)
  TEST_EXCEPTION(Exception::InvalidRange,
    exp.rasterizeRTMZ(output.data(), 10, 10, 100.0, 100.0, 0.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
  
  TEST_EXCEPTION(Exception::InvalidRange,
    exp.rasterizeRTMZ(output.data(), 10, 10, 100.0, 50.0, 0.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Exception: inverted m/z range))
{
  MSExperiment exp = createTestExperiment();
  vector<float> output(100, 0.0f);
  
  // Test inverted m/z range (min >= max)
  TEST_EXCEPTION(Exception::InvalidRange,
    exp.rasterizeRTMZ(output.data(), 10, 10, 0.0, 100.0, 1000.0, 1000.0, 1, MSExperiment::RasterAggregation::SUM)
  )
  
  TEST_EXCEPTION(Exception::InvalidRange,
    exp.rasterizeRTMZ(output.data(), 10, 10, 0.0, 100.0, 1000.0, 500.0, 1, MSExperiment::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeRTMZ - Edge case: out of range))
{
  MSExperiment exp = createTestExperiment();
  Size rt_bins = 10;
  Size mz_bins = 10;
  vector<float> output(rt_bins * mz_bins, 0.0f);
  
  // Request range outside experiment data
  exp.rasterizeRTMZ(output.data(), rt_bins, mz_bins,
                  1000.0, 2000.0,  // RT outside data
                  500.0, 600.0,
                  1,
                  MSExperiment::RasterAggregation::SUM);
  
  // Output should be all zeros
  for (float val : output)
  {
    TEST_EQUAL(val, 0.0f)
  }
}
END_SECTION

START_SECTION((void rasterizeRTMZ - SUM vs MAX aggregation))
{
  // Create experiment with overlapping peaks
  MSExperiment exp;
  
  for (Size i = 0; i < 5; ++i)
  {
    MSSpectrum spec;
    spec.setRT(100.0 + i * 2.0);  // Close together to get same bin
    spec.setMSLevel(1);
    
    Peak1D peak;
    peak.setMZ(500.0);  // Same m/z
    peak.setIntensity(1000.0);
    spec.push_back(peak);
    
    exp.addSpectrum(spec);
  }
  
  Size rt_bins = 2;
  Size mz_bins = 2;
  vector<float> output_sum(rt_bins * mz_bins, 0.0f);
  vector<float> output_max(rt_bins * mz_bins, 0.0f);
  
  exp.rasterizeRTMZ(output_sum.data(), rt_bins, mz_bins,
                  100.0, 110.0,
                  490.0, 510.0,
                  1,
                  MSExperiment::RasterAggregation::SUM);
  
  exp.rasterizeRTMZ(output_max.data(), rt_bins, mz_bins,
                  100.0, 110.0,
                  490.0, 510.0,
                  1,
                  MSExperiment::RasterAggregation::MAX);
  
  // Find non-zero bin
  float sum_value = 0.0f;
  float max_value = 0.0f;
  for (Size i = 0; i < output_sum.size(); ++i)
  {
    if (output_sum[i] > sum_value) sum_value = output_sum[i];
    if (output_max[i] > max_value) max_value = output_max[i];
  }
  
  // SUM should accumulate per-RT-bin: bin with most entries gets 3 * 1000 = 3000
  // MAX should only keep max single intensity (1000)
  TEST_REAL_SIMILAR(sum_value, 3000.0f)
  TEST_REAL_SIMILAR(max_value, 1000.0f)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
