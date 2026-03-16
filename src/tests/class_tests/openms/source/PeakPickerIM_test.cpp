// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/test_config.h>
#include <algorithm>
#include <cmath>

///////////////////////////
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerIM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerIM* ptr = nullptr;
PeakPickerIM* nullPointer = nullptr;

START_SECTION((PeakPickerIM()))
ptr = new PeakPickerIM();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerIM()))
delete ptr;
END_SECTION


MSSpectrum input;

// Synthetic 2D Gaussian peak around m/z 300 with two IM peaks (0.8, 1.0)
std::vector<double> mzs;
std::vector<double> IMs;
std::vector<double> Intensities;
// --- Parameters for m/z Gaussian ---
const double mz_center = 300.0;
const double mz_sigma = 0.0025;
const double mz_min = 299.99;
const double mz_max = 300.01;
const double mz_step = 0.0005;
// --- Parameters for IM Gaussians ---
const double im_center1 = 0.8;
const double im_center2 = 1.0;
const double im_sigma = 0.03;
const double im_min = 0.7;
const double im_max = 1.1;
const double im_step = 0.01;

// Build the 2D grid: for each m/z, duplicate across all IM scans
for (double mz = mz_min; mz <= mz_max + 1e-9; mz += mz_step)
{
  // 1D Gaussian in m/z
  const double mz_int = std::exp(-0.5 * std::pow((mz - mz_center) / mz_sigma, 2.0));
  for (double im = im_min; im <= im_max + 1e-9; im += im_step)
  {
    // Two Gaussians in IM dimension
    const double im_int1 = std::exp(-0.5 * std::pow((im - im_center1) / im_sigma, 2.0));
    const double im_int2 = std::exp(-0.5 * std::pow((im - im_center2) / im_sigma, 2.0));
    const double im_int = im_int1 + im_int2;
    // Final separable 2D intensity
    const double intensity = mz_int * im_int;
    mzs.push_back(mz);
    IMs.push_back(im);
    Intensities.push_back(intensity);
  }
}

for (size_t i = 0; i < mzs.size(); ++i)
{
  input.emplace_back(mzs[i], Intensities[i]);
}
input.getFloatDataArrays().resize(1);
input.getFloatDataArrays()[0].setName(Constants::UserParam::ION_MOBILITY);
for (size_t i = 0; i < IMs.size(); ++i)
{
  input.getFloatDataArrays()[0].push_back(IMs[i]);
}
input.setIMFormat(IMFormat::CONCATENATED);

START_SECTION(void pickIMTraces(MSSpectrum& spectrum))
{
  PeakPickerIM pp_im;
  MSSpectrum spec = input;

  pp_im.pickIMTraces(spec);

  TEST_EQUAL(spec.size(), 2)
  TEST_REAL_SIMILAR(spec[0].getIntensity(), 5.70646)
  TEST_REAL_SIMILAR(spec[0].getMZ(), 300.0)

  TEST_EQUAL(spec.getFloatDataArrays().size(), 1)
  TEST_EQUAL(spec.getFloatDataArrays()[0].getName(), Constants::UserParam::ION_MOBILITY_CENTROID)
  TEST_REAL_SIMILAR(spec.getFloatDataArrays()[0][0], 0.8)
}
END_SECTION

START_SECTION(void pickIMSage(MSSpectrum& spectrum) const)
{
  PeakPickerIM pp_im;
  Param p;
  p.setValue("pickIMSage:im_tolerance_sage", 0.15);
  p.setValue("pickIMSage:ppm_tolerance_sage", 200.0);
  pp_im.setParameters(p);

  MSSpectrum local_input = input;
  pp_im.pickIMSage(local_input);

  TEST_EQUAL(local_input.size(), 2)

  TEST_REAL_SIMILAR(local_input[0].getMZ(), 300.0)
  TEST_REAL_SIMILAR(local_input[1].getMZ(), 300.0)

  TEST_EQUAL(local_input.getFloatDataArrays().size(), 1)
  TEST_EQUAL(local_input.getFloatDataArrays()[0].getName(), Constants::UserParam::ION_MOBILITY_CENTROID)

  double im1 = local_input.getFloatDataArrays()[0][0];
  double im2 = local_input.getFloatDataArrays()[0][1];
  if (im1 > im2) std::swap(im1, im2);

  TEST_REAL_SIMILAR(im1, 0.8)
  TEST_REAL_SIMILAR(im2, 1.0)
}
END_SECTION

START_SECTION(void pickIMCluster(MSSpectrum& spectrum) const)
{
  PeakPickerIM pp_im;
  Param p;
  p.setValue("pickIMCluster:im_tolerance_cluster", 0.15);
  p.setValue("pickIMCluster:ppm_tolerance_cluster", 200.0);
  pp_im.setParameters(p);

  MSSpectrum local_input = input;
  pp_im.pickIMCluster(local_input);

  // Cluster might produce extra peaks in the tails, but should have the main ones
  TEST_TRUE(local_input.size() >= 2)

  // Find the two most intense peaks
  std::vector<size_t> indices(local_input.size());
  for (size_t i = 0; i < local_input.size(); ++i)
    indices[i] = i;

  std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return local_input[a].getIntensity() > local_input[b].getIntensity(); });

  double im1 = local_input.getFloatDataArrays()[0][indices[0]];
  double im2 = local_input.getFloatDataArrays()[0][indices[1]];

  if (im1 > im2) std::swap(im1, im2);

  TOLERANCE_ABSOLUTE(0.001)
  TEST_REAL_SIMILAR(im1, 0.8)
  TEST_REAL_SIMILAR(im2, 1.0)
}
END_SECTION

END_TEST