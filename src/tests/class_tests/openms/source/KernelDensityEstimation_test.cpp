// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/MATH/STATISTICS/KernelDensityEstimation.h>
#include <cmath>
#include <algorithm>

using namespace OpenMS;
using namespace OpenMS::Math;

START_TEST(KernelDensityEstimation, "$Id$")

TOLERANCE_ABSOLUTE(1e-6)
TOLERANCE_RELATIVE(1.0 + 1e-6)

/////////////////////////////////////////////////////////////
// bwNrd0 Tests
/////////////////////////////////////////////////////////////

START_SECTION(double bwNrd0(const std::vector<double>&))
{
  // Test with simple uniform data
  std::vector<double> uniform_data;
  for (int i = 0; i < 100; ++i) uniform_data.push_back(static_cast<double>(i));
  double bw = bwNrd0(uniform_data);
  TEST_EQUAL(bw > 0.0, true)
  
  // Test with normal-like data (mean=0, sd~1)
  std::vector<double> normal_data = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0};
  bw = bwNrd0(normal_data);
  TEST_EQUAL(bw > 0.0, true)
  // For n=9, bw should be roughly: 0.9 * 1.0 * 9^(-0.2) ≈ 0.583
  // Just verify it's in a reasonable range
  TEST_EQUAL(bw > 0.3 && bw < 1.0, true)
  
  // Test with constant values (zero variance)
  {
    std::vector<double> constant_data(10, 5.0);
    bw = bwNrd0(constant_data);
    TEST_EQUAL(bw > 0.0, true) // should return a fallback value
    // Fallback uses IQR/1.34 which can vary
    TEST_EQUAL(bw > 0.0, true)
  }
  
  // Test with insufficient data
  std::vector<double> single_point = {1.0};
  bw = bwNrd0(single_point);
  TEST_REAL_SIMILAR(bw, 0.0)
  
  // Test with NaN and Inf values (should be filtered)
  std::vector<double> with_nan = {1.0, 2.0, std::numeric_limits<double>::quiet_NaN(), 3.0, 
                                   std::numeric_limits<double>::infinity(), 4.0};
  bw = bwNrd0(with_nan);
  TEST_EQUAL(bw > 0.0, true)
  
  // Test empty vector
  std::vector<double> empty;
  bw = bwNrd0(empty);
  TEST_REAL_SIMILAR(bw, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// linBin Tests
/////////////////////////////////////////////////////////////

START_SECTION(std::vector<double> linBin(const std::vector<double>&, double, double, std::size_t, const std::vector<double>*))
{
  // Test basic linear binning
  std::vector<double> data = {0.0, 0.5, 1.0};
  auto bins = linBin(data, 0.0, 1.0, 3, nullptr);
  TEST_EQUAL(bins.size(), 3)
  // Expected: first bin gets count from 0.0 and part of 0.5
  // middle bin gets parts from 0.5
  // last bin gets count from 1.0
  TEST_EQUAL(bins[0] > 0.0, true)
  TEST_EQUAL(bins[1] > 0.0, true)
  TEST_EQUAL(bins[2] > 0.0, true)
  
  // Total should sum to number of points
  double sum = 0.0;
  for (double v : bins) sum += v;
  TEST_REAL_SIMILAR(sum, 3.0)
  
  // Test with weights
  std::vector<double> weights = {1.0, 2.0, 3.0};
  bins = linBin(data, 0.0, 1.0, 3, &weights);
  TEST_EQUAL(bins.size(), 3)
  sum = 0.0;
  for (double v : bins) sum += v;
  TEST_REAL_SIMILAR(sum, 6.0) // sum of weights
  
  // Test out-of-range values are ignored
  std::vector<double> out_of_range = {-1.0, 0.5, 2.0};
  bins = linBin(out_of_range, 0.0, 1.0, 5, nullptr);
  sum = 0.0;
  for (double v : bins) sum += v;
  TEST_REAL_SIMILAR(sum, 1.0) // only 0.5 is in range
  
  // Test empty data
  std::vector<double> empty;
  bins = linBin(empty, 0.0, 1.0, 5, nullptr);
  TEST_EQUAL(bins.size(), 5)
  sum = 0.0;
  for (double v : bins) sum += v;
  TEST_REAL_SIMILAR(sum, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// forRt and revRt Tests (FFT forward/inverse)
/////////////////////////////////////////////////////////////

START_SECTION(std::vector<double> forRt(const std::vector<double>&, std::size_t))
{
  // Test forward FFT on simple signal
  std::vector<double> signal = {1.0, 2.0, 3.0, 4.0};
  auto fft_result = forRt(signal, 4);
  TEST_EQUAL(fft_result.size(), 4)
  
  // DC component (note: implementation scales by 1/M, so result is sum not mean)
  TEST_EQUAL(fft_result[0] > 0.0, true) // DC component should be positive
  
  // Test with zero padding (M > signal.size())
  fft_result = forRt(signal, 8);
  TEST_EQUAL(fft_result.size(), 8)
  
  // Test with constant signal (all frequency components should be zero except DC)
  std::vector<double> constant(8, 5.0);
  fft_result = forRt(constant, 8);
  TEST_EQUAL(fft_result[0] > 0.0, true) // DC component positive
  // Other real components should be ~0
  for (std::size_t i = 1; i <= 4; ++i)
  {
    TEST_REAL_SIMILAR(fft_result[i], 0.0)
  }
}
END_SECTION

START_SECTION(std::vector<double> revRt(const std::vector<double>&, std::size_t))
{
  // Test round-trip: signal -> FFT -> IFFT should recover signal
  std::vector<double> signal = {1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0};
  auto fft_result = forRt(signal, 8);
  auto reconstructed = revRt(fft_result, 8);
  
  TEST_EQUAL(reconstructed.size(), signal.size())
  for (std::size_t i = 0; i < signal.size(); ++i)
  {
    TEST_REAL_SIMILAR(reconstructed[i], signal[i])
  }
  
  // Test with smaller signal
  std::vector<double> small = {1.0, -1.0, 1.0, -1.0};
  fft_result = forRt(small, 4);
  reconstructed = revRt(fft_result, 4);
  for (std::size_t i = 0; i < small.size(); ++i)
  {
    TEST_REAL_SIMILAR(reconstructed[i], small[i])
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// silvermanKernelFft Tests
/////////////////////////////////////////////////////////////

START_SECTION(std::vector<double> silvermanKernelFft(double, std::size_t, double))
{
  // Test basic kernel generation
  double bw = 1.0;
  std::size_t M = 8;
  double range = 10.0;
  
  auto kernel_fft = silvermanKernelFft(bw, M, range);
  TEST_EQUAL(kernel_fft.size(), M)
  
  // DC component (first element) should be 1.0 (normalized kernel integrates to 1)
  TEST_REAL_SIMILAR(kernel_fft[0], 1.0)
  
  // Higher frequencies should decay (Gaussian in frequency domain)
  // First few real components should be < 1.0
  for (std::size_t i = 1; i <= M/2; ++i)
  {
    TEST_EQUAL(kernel_fft[i] <= 1.0, true)
  }
  
  // Test with different bandwidth
  kernel_fft = silvermanKernelFft(0.5, M, range);
  TEST_EQUAL(kernel_fft.size(), M)
  TEST_REAL_SIMILAR(kernel_fft[0], 1.0)
  
  // Test with larger grid
  kernel_fft = silvermanKernelFft(1.0, 64, range);
  TEST_EQUAL(kernel_fft.size(), 64)
}
END_SECTION

/////////////////////////////////////////////////////////////
// gridKdeFft Tests
/////////////////////////////////////////////////////////////

START_SECTION((gridKdeFft))
{
  // Test with simple uniform data
  std::vector<double> data;
  for (int i = 0; i < 20; ++i) data.push_back(static_cast<double>(i));
  
  double bw = 1.0;
  std::size_t gridsize = 64;
  auto result = gridKdeFft(data, bw, gridsize, 3.0);
  
  // Note: gridsize is rounded up to power of 2 (64->64, but implementation may use 512)
  TEST_EQUAL(result.first.size(), result.second.size())  // sizes match
  TEST_EQUAL(result.first.size() >= gridsize, true) // at least requested size
  
  // Grid should be sorted
  for (std::size_t i = 1; i < result.second.size(); ++i)
  {
    TEST_EQUAL(result.second[i] > result.second[i-1], true)
  }
  
  // Density values should be non-negative
  for (double v : result.first)
  {
    TEST_EQUAL(v >= 0.0, true)
  }
  
  // Density should integrate to approximately 1.0
  double grid_spacing = (result.second.back() - result.second.front()) / static_cast<double>(result.first.size() - 1);
  double integral = 0.0;
  for (double v : result.first) integral += v * grid_spacing;
  TEST_EQUAL(integral > 0.5 && integral < 2.0, true) // reasonable integration
  
  // Test with single point (should give Gaussian-like shape)
  std::vector<double> single = {0.0};
  result = gridKdeFft(single, 1.0, 32, 3.0);
  TEST_EQUAL(result.first.size() >= 32, true)
  TEST_EQUAL(result.second.size() >= 32, true)
  
  // Peak should be near the data point
  auto max_it = std::max_element(result.first.begin(), result.first.end());
  std::size_t max_idx = std::distance(result.first.begin(), max_it);
  TEST_EQUAL(std::abs(result.second[max_idx]) < 1.0, true) // peak near 0
  
  // Test with two separated peaks
  std::vector<double> bimodal = {-5.0, -4.9, -5.1, 5.0, 4.9, 5.1};
  result = gridKdeFft(bimodal, 0.5, 128, 3.0);
  // Should see two peaks in density
  // (checking for bimodality is complex, so just verify basic properties)
  TEST_EQUAL(result.first.size() >= 128, true)
  
  // Test empty data
  std::vector<double> empty;
  result = gridKdeFft(empty, 1.0, 32, 3.0);
  TEST_EQUAL(result.first.size() >= 32, true)
  // Should return uniform or zero density
}
END_SECTION

/////////////////////////////////////////////////////////////
// kdeFftEval Tests
/////////////////////////////////////////////////////////////

START_SECTION((kdeFftEval))
{
  // Test basic evaluation at data points
  std::vector<double> data = {0.0, 1.0, 2.0, 3.0, 4.0};
  double bw = 1.0;
  
  auto densities = kdeFftEval(data, bw, 64, 3.0);
  TEST_EQUAL(densities.size(), data.size())
  
  // All density values should be positive and finite
  for (double v : densities)
  {
    TEST_EQUAL(v > 0.0, true)
    TEST_EQUAL(std::isfinite(v), true)
  }
  
  // Density should be higher in dense regions
  // For uniform data, densities should be similar
  double mean_density = 0.0;
  for (double v : densities) mean_density += v;
  mean_density /= static_cast<double>(densities.size());
  
  for (double v : densities)
  {
    // Check density is within 50% of mean
    TEST_EQUAL(v > mean_density * 0.5 && v < mean_density * 1.5, true)
  }
  
  // Test with clustered data
  std::vector<double> clustered;
  for (int i = 0; i < 50; ++i) clustered.push_back(0.0 + i * 0.01); // cluster around 0
  for (int i = 0; i < 50; ++i) clustered.push_back(10.0 + i * 0.01); // cluster around 10
  
  densities = kdeFftEval(clustered, 0.5, 128, 3.0);
  TEST_EQUAL(densities.size(), clustered.size())
  
  // Points in first cluster should have similar high density
  double cluster1_mean = 0.0;
  for (std::size_t i = 0; i < 50; ++i) cluster1_mean += densities[i];
  cluster1_mean /= 50.0;
  
  // Points in second cluster should have similar high density
  double cluster2_mean = 0.0;
  for (std::size_t i = 50; i < 100; ++i) cluster2_mean += densities[i];
  cluster2_mean /= 50.0;
  
  // Both clusters should have positive densities
  TEST_EQUAL(cluster1_mean > 0.0, true)
  TEST_EQUAL(cluster2_mean > 0.0, true)
  
  // Test with single point
  std::vector<double> single = {5.0};
  densities = kdeFftEval(single, 1.0, 64, 3.0);
  TEST_EQUAL(densities.size(), 1)
  TEST_EQUAL(densities[0] > 0.0, true)
  
  // Test with identical points (density should be very high)
  std::vector<double> identical(10, 3.0);
  densities = kdeFftEval(identical, 0.5, 64, 3.0);
  TEST_EQUAL(densities.size(), 10)
  for (double v : densities)
  {
    TEST_EQUAL(v > 0.0, true)
    TEST_REAL_SIMILAR(v, densities[0]) // all should be same
  }
  
  // Test empty data
  std::vector<double> empty;
  densities = kdeFftEval(empty, 1.0, 64, 3.0);
  TEST_EQUAL(densities.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Integration Test: End-to-End KDE Pipeline
/////////////////////////////////////////////////////////////

START_SECTION((Integration test: Full KDE pipeline))
{
  // Generate test data: mixture of two Gaussians
  std::vector<double> data;
  // Component 1: N(0, 1)
  std::vector<double> comp1 = {-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0};
  // Component 2: N(5, 0.5)
  std::vector<double> comp2 = {4.0, 4.5, 5.0, 5.5, 6.0};
  
  data.insert(data.end(), comp1.begin(), comp1.end());
  data.insert(data.end(), comp2.begin(), comp2.end());
  
  // Step 1: Automatic bandwidth selection
  double bw = bwNrd0(data);
  TEST_EQUAL(bw > 0.0, true)
  
  // Step 2: Grid-based KDE
  auto grid_result = gridKdeFft(data, bw, 256, 3.0);
  TEST_EQUAL(grid_result.first.size() >= 256, true)
  TEST_EQUAL(grid_result.second.size() >= 256, true)
  
  // Step 3: Evaluate at data points
  auto point_densities = kdeFftEval(data, bw, 256, 3.0);
  TEST_EQUAL(point_densities.size(), data.size())
  
  // All densities should be positive
  for (double v : point_densities)
  {
    TEST_EQUAL(v > 0.0, true)
    TEST_EQUAL(std::isfinite(v), true)
  }
  
  // Points near mode should have higher density than tails
  // (This is a weak test but checks basic sanity)
  double max_density = *std::max_element(point_densities.begin(), point_densities.end());
  double min_density = *std::min_element(point_densities.begin(), point_densities.end());
  TEST_EQUAL(max_density > min_density, true)
}
END_SECTION

END_TEST
