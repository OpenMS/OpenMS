// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

/// Creates a test spectrum with IM data in IM_PEAK format (timsTOF-style)
MSSpectrum createTestIMFrame()
{
  MSSpectrum spec;
  spec.setRT(100.0);
  spec.setMSLevel(1);

  // Add peaks with corresponding IM values
  // Note: IM frame data is typically NOT sorted by m/z
  vector<double> mz_values = {500.0, 550.0, 600.0, 500.5, 550.5, 600.5};
  vector<double> int_values = {1000.0, 2000.0, 3000.0, 1500.0, 2500.0, 3500.0};
  vector<float> im_values = {1.0f, 1.2f, 1.4f, 1.1f, 1.3f, 1.5f};

  for (Size i = 0; i < mz_values.size(); ++i)
  {
    Peak1D peak;
    peak.setMZ(mz_values[i]);
    peak.setIntensity(int_values[i]);
    spec.push_back(peak);
  }

  // Add IM float data array with CV term for ion mobility
  MSSpectrum::FloatDataArray im_array;
  im_array.setName("Ion Mobility");
  // Set CV term to mark this as IM data (child of MS:1002893)
  im_array.setMetaValue("unit_accession", "MS:1002815"); // inverse reduced ion mobility
  for (float im : im_values)
  {
    im_array.push_back(im);
  }
  spec.getFloatDataArrays().push_back(im_array);

  return spec;
}

START_TEST(MSSpectrum_rasterizeIMFrame, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((void rasterizeIMFrame(float* output, Size im_bins, Size mz_bins, CoordinateType min_im, CoordinateType max_im, CoordinateType min_mz, CoordinateType max_mz, RasterAggregation aggregation) const))
{
  MSSpectrum spec = createTestIMFrame();

  // Test basic rasterization with SUM aggregation
  {
    Size im_bins = 10;
    Size mz_bins = 10;
    vector<float> output(im_bins * mz_bins, 0.0f);

    spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                          1.0, 1.6,     // IM range (includes all peaks 1.0-1.5)
                          500.0, 601.0, // m/z range (includes all peaks 500.0-600.5)
                          MSSpectrum::RasterAggregation::SUM);

    // Total intensity should be sum of all peak intensities
    float total_sum = 0.0f;
    for (float val : output) total_sum += val;
    // 1000 + 2000 + 3000 + 1500 + 2500 + 3500 = 13500
    TEST_REAL_SIMILAR(total_sum, 13500.0f)
  }

  // Test with MAX aggregation
  {
    Size im_bins = 10;
    Size mz_bins = 10;
    vector<float> output(im_bins * mz_bins, 0.0f);

    spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                          1.0, 1.6,
                          500.0, 601.0,
                          MSSpectrum::RasterAggregation::MAX);

    // The maximum bin value should be at most the largest single peak intensity (3500)
    float max_val = 0.0f;
    for (float val : output) if (val > max_val) max_val = val;
    TEST_REAL_SIMILAR(max_val, 3500.0f)
  }

  // Test with empty spectrum (but with IM array)
  {
    MSSpectrum empty_spec;
    MSSpectrum::FloatDataArray im_array;
    im_array.setName("Ion Mobility");
    im_array.setMetaValue("unit_accession", "MS:1002815");
    empty_spec.getFloatDataArrays().push_back(im_array);

    Size im_bins = 5;
    Size mz_bins = 5;
    vector<float> output(im_bins * mz_bins, 1.0f);  // Initialize with non-zero

    empty_spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                                0.0, 2.0,
                                0.0, 1000.0,
                                MSSpectrum::RasterAggregation::SUM);

    // Output should be all zeros
    for (float val : output)
    {
      TEST_EQUAL(val, 0.0f)
    }
  }

  // Test with single peak
  {
    MSSpectrum single_spec;
    Peak1D peak;
    peak.setMZ(500.0);
    peak.setIntensity(100.0);
    single_spec.push_back(peak);

    MSSpectrum::FloatDataArray im_array;
    im_array.setName("Ion Mobility");
    im_array.setMetaValue("unit_accession", "MS:1002815");
    im_array.push_back(1.5f);  // Single IM value
    single_spec.getFloatDataArrays().push_back(im_array);

    Size im_bins = 10;
    Size mz_bins = 10;
    vector<float> output(im_bins * mz_bins, 0.0f);

    single_spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                                 1.0, 2.0,
                                 400.0, 600.0,
                                 MSSpectrum::RasterAggregation::SUM);

    // Should have exactly one bin with intensity 100
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

  // Test peaks at exact boundary values (should be included in last bin)
  {
    MSSpectrum boundary_spec;
    Peak1D peak;
    peak.setMZ(600.0);  // Exactly at max_mz
    peak.setIntensity(1000.0);
    boundary_spec.push_back(peak);

    MSSpectrum::FloatDataArray im_array;
    im_array.setName("Ion Mobility");
    im_array.setMetaValue("unit_accession", "MS:1002815");
    im_array.push_back(2.0f);  // Exactly at max_im
    boundary_spec.getFloatDataArrays().push_back(im_array);

    Size im_bins = 10;
    Size mz_bins = 10;
    vector<float> output(im_bins * mz_bins, 0.0f);

    boundary_spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                                   1.0, 2.0,    // IM range
                                   500.0, 600.0, // m/z range
                                   MSSpectrum::RasterAggregation::SUM);

    // The peak at exact boundary should be in the last bin
    float total_sum = 0.0f;
    for (float val : output) total_sum += val;
    TEST_REAL_SIMILAR(total_sum, 1000.0f)
  }
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: null pointer))
{
  MSSpectrum spec = createTestIMFrame();

  TEST_EXCEPTION(Exception::NullPointer,
    spec.rasterizeIMFrame(nullptr, 10, 10, 0.0, 2.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: zero im_bins))
{
  MSSpectrum spec = createTestIMFrame();
  vector<float> output(100, 0.0f);

  TEST_EXCEPTION(Exception::InvalidValue,
    spec.rasterizeIMFrame(output.data(), 0, 10, 0.0, 2.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: zero mz_bins))
{
  MSSpectrum spec = createTestIMFrame();
  vector<float> output(100, 0.0f);

  TEST_EXCEPTION(Exception::InvalidValue,
    spec.rasterizeIMFrame(output.data(), 10, 0, 0.0, 2.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: inverted IM range))
{
  MSSpectrum spec = createTestIMFrame();
  vector<float> output(100, 0.0f);

  // Test inverted IM range (min >= max)
  TEST_EXCEPTION(Exception::InvalidRange,
    spec.rasterizeIMFrame(output.data(), 10, 10, 2.0, 2.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )

  TEST_EXCEPTION(Exception::InvalidRange,
    spec.rasterizeIMFrame(output.data(), 10, 10, 2.0, 1.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: inverted m/z range))
{
  MSSpectrum spec = createTestIMFrame();
  vector<float> output(100, 0.0f);

  // Test inverted m/z range (min >= max)
  TEST_EXCEPTION(Exception::InvalidRange,
    spec.rasterizeIMFrame(output.data(), 10, 10, 0.0, 2.0, 1000.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )

  TEST_EXCEPTION(Exception::InvalidRange,
    spec.rasterizeIMFrame(output.data(), 10, 10, 0.0, 2.0, 1000.0, 500.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Exception: missing IM data))
{
  // Spectrum without IM data array
  MSSpectrum spec;
  Peak1D peak;
  peak.setMZ(500.0);
  peak.setIntensity(100.0);
  spec.push_back(peak);

  vector<float> output(100, 0.0f);

  TEST_EXCEPTION(Exception::MissingInformation,
    spec.rasterizeIMFrame(output.data(), 10, 10, 0.0, 2.0, 0.0, 1000.0, MSSpectrum::RasterAggregation::SUM)
  )
}
END_SECTION

START_SECTION((void rasterizeIMFrame - Edge case: out of range))
{
  MSSpectrum spec = createTestIMFrame();
  Size im_bins = 10;
  Size mz_bins = 10;
  vector<float> output(im_bins * mz_bins, 0.0f);

  // Request IM range outside spectrum data
  spec.rasterizeIMFrame(output.data(), im_bins, mz_bins,
                        10.0, 20.0,  // IM outside data (data is 1.0-1.5)
                        500.0, 600.0,
                        MSSpectrum::RasterAggregation::SUM);

  // Output should be all zeros
  for (float val : output)
  {
    TEST_EQUAL(val, 0.0f)
  }
}
END_SECTION

START_SECTION((void rasterizeIMFrame - SUM vs MAX aggregation))
{
  // Create spectrum with overlapping peaks (same IM and m/z bin)
  MSSpectrum spec;

  for (Size i = 0; i < 5; ++i)
  {
    Peak1D peak;
    peak.setMZ(500.0 + i * 0.01);  // Very close m/z values
    peak.setIntensity(1000.0);
    spec.push_back(peak);
  }

  // All peaks have similar IM values to fall in the same bin
  MSSpectrum::FloatDataArray im_array;
  im_array.setName("Ion Mobility");
  im_array.setMetaValue("unit_accession", "MS:1002815");
  for (Size i = 0; i < 5; ++i)
  {
    im_array.push_back(1.5f + i * 0.01f);  // Very close IM values
  }
  spec.getFloatDataArrays().push_back(im_array);

  Size im_bins = 2;
  Size mz_bins = 2;
  vector<float> output_sum(im_bins * mz_bins, 0.0f);
  vector<float> output_max(im_bins * mz_bins, 0.0f);

  spec.rasterizeIMFrame(output_sum.data(), im_bins, mz_bins,
                        1.0, 2.0,
                        490.0, 510.0,
                        MSSpectrum::RasterAggregation::SUM);

  spec.rasterizeIMFrame(output_max.data(), im_bins, mz_bins,
                        1.0, 2.0,
                        490.0, 510.0,
                        MSSpectrum::RasterAggregation::MAX);

  // Find non-zero bin values
  float sum_value = 0.0f;
  float max_value = 0.0f;
  for (Size i = 0; i < output_sum.size(); ++i)
  {
    if (output_sum[i] > sum_value) sum_value = output_sum[i];
    if (output_max[i] > max_value) max_value = output_max[i];
  }

  // SUM should accumulate all 5 peaks: 5 * 1000 = 5000
  // MAX should only keep max single intensity: 1000
  TEST_REAL_SIMILAR(sum_value, 5000.0f)
  TEST_REAL_SIMILAR(max_value, 1000.0f)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
