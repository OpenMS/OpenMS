// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>

using namespace OpenMS;
using namespace std;

template <class SpectrumT>
void check_results(SpectrumT spec)
{
  double sum = 0.0;
  for (Size i = 0; i < spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 3 + 2);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 4 + 2.0 / 3 * 8);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 1.0 / 3 *8 +2 + 1.0 / 3);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 2.0 / 3);
}

///////////////////////////

START_TEST(LinearResamplerAlign, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum input_spectrum;
input_spectrum.resize(5);
input_spectrum[0].setMZ(0);
input_spectrum[0].setIntensity(3.0f);
input_spectrum[1].setMZ(0.5);
input_spectrum[1].setIntensity(6.0f);
input_spectrum[2].setMZ(1.);
input_spectrum[2].setIntensity(8.0f);
input_spectrum[3].setMZ(1.6);
input_spectrum[3].setIntensity(2.0f);
input_spectrum[4].setMZ(1.8);
input_spectrum[4].setIntensity(1.0f);

// A spacing of 0.75 will lead to a recalculation of intensities, each
// resampled point gets intensities from raw data points that are at most +/-
// spacing away.
double default_spacing = 0.75;
#if 1

START_SECTION(( template < template< typename > class SpecT, typename PeakType > void raster(SpecT< PeakType > &spectrum)))
{

  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", default_spacing);
  lr.setParameters(param);
  lr.raster(spec);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 3+2);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 4+2.0/3*8);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 1.0/3*8+2+1.0/3);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 2.0 / 3);
}
END_SECTION

// it should also work with chromatograms
START_SECTION([EXTRA] test_linear_res_chromat)
{
  MSChromatogram chrom;
  chrom.resize(5);
  chrom[0].setRT(0);
  chrom[0].setIntensity(3.0f);
  chrom[1].setRT(0.5);
  chrom[1].setIntensity(6.0f);
  chrom[2].setRT(1.);
  chrom[2].setIntensity(8.0f);
  chrom[3].setRT(1.6);
  chrom[3].setIntensity(2.0f);
  chrom[4].setRT(1.8);
  chrom[4].setIntensity(1.0f);

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", default_spacing);
  lr.setParameters(param);
  lr.raster(chrom);

  check_results(chrom);
}
END_SECTION

START_SECTION(( void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end)))
{

  MSSpectrum spec = input_spectrum;
  MSSpectrum output_spectrum;
  output_spectrum.resize(4);

  // We want to resample the input spectrum at these m/z positions: 0, 0.75, 1.5 and 2.25
  std::vector<double> mz_res_data(4);
  std::vector<double> int_res_data(4);
  output_spectrum[0].setMZ(0);
  output_spectrum[1].setMZ(0.75);
  output_spectrum[2].setMZ(1.5);
  output_spectrum[3].setMZ(2.25);

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", default_spacing);
  lr.setParameters(param);
  lr.raster(spec.begin(), spec.end(), output_spectrum.begin(), output_spectrum.end());

  check_results(output_spectrum);
}

END_SECTION

// it should also work with data vectors
START_SECTION( ( template <typename PeakTypeIterator, typename ConstPeakTypeIterator>
     void raster(ConstPeakTypeIterator mz_raw_it, ConstPeakTypeIterator mz_raw_end,
       ConstPeakTypeIterator int_raw_it, ConstPeakTypeIterator int_raw_end,
       PeakTypeIterator mz_resample_it, PeakTypeIterator mz_resample_end,
       PeakTypeIterator int_resample_it, PeakTypeIterator int_resample_end
       )
    ))
{
  MSChromatogram spec;

  std::vector<double> mz_data(5);
  std::vector<double> int_data(5);

  mz_data[0] = 0;
  mz_data[1] = 0.5;
  mz_data[2] = 1.;
  mz_data[3] = 1.6;
  mz_data[4] = 1.8;
  int_data[0] = 3.0f;
  int_data[1] = 6.0f;
  int_data[2] = 8.0f;
  int_data[3] = 2.0f;
  int_data[4] = 1.0f;

  // We want to resample the input spectrum at these m/z positions: 0, 0.75, 1.5 and 2.25
  std::vector<double> mz_res_data(4);
  std::vector<double> int_res_data(4);
  mz_res_data[0] = 0;
  mz_res_data[1] = 0.75;
  mz_res_data[2] = 1.5;
  mz_res_data[3] = 2.25;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", default_spacing);
  lr.setParameters(param);

  lr.raster(mz_data.begin(), mz_data.end(), int_data.begin(), int_data.end(),
  mz_res_data.begin(), mz_res_data.end(), int_res_data.begin(), int_res_data.end());

  // check_results(spec);
  double sum = 0.0;
  for (Size i=0; i<int_res_data.size(); ++i)
  {
    sum += int_res_data[i];
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(int_res_data[0], 3+2);
  TEST_REAL_SIMILAR(int_res_data[1], 4+2.0/3*8);
  TEST_REAL_SIMILAR(int_res_data[2], 1.0/3*8+2+1.0/3);
  TEST_REAL_SIMILAR(int_res_data[3], 2.0 / 3);
}
END_SECTION

// it should work with alignment to 0, 1.8 and give the same result
START_SECTION((template < template< typename > class SpecT, typename PeakType > void raster_align(SpecT< PeakType > &spectrum, double start_pos, double end_pos)))
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.75);
  lr.setParameters(param);

  lr.raster_align(spec, 0, 1.8);
  check_results(spec);
}
END_SECTION

// it should work with alignment to -0.25, 1.8
START_SECTION([EXTRA] test_linear_res_align_3)
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
  lr.raster_align(spec, -0.25, 1.8);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 1.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 1.5+3);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 3+4);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 4.0+ 0.6);
  TEST_REAL_SIMILAR(spec[4].getIntensity(), 1.4 + 0.9 );
  TEST_REAL_SIMILAR(spec[5].getIntensity(), 0.1 );
}
END_SECTION

// it should work with alignment to -2.25, 1.8
START_SECTION([EXTRA] test_linear_res_align_4)
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.75);
  lr.setParameters(param);
  lr.raster_align(spec, -2.25, 1.8);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 0);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 0);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 0);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 3+2);
  TEST_REAL_SIMILAR(spec[4].getIntensity(), 4+2.0/3*8);
  TEST_REAL_SIMILAR(spec[5].getIntensity(), 1.0/3*8+2+1.0/3);
  TEST_REAL_SIMILAR(spec[6].getIntensity(), 2.0 / 3);
}
END_SECTION

// it should work with alignment to -0.25, 1.25
START_SECTION([EXTRA] test_linear_res_align_5)
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_align(spec, -0.25, 1.25);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20 - 2.4 -0.6); // missing points 1.75 and 2.25 which have intensity 2.4 together

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 1.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 1.5+3);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 3+4);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 4.0);//+ 0.6);
}
END_SECTION

// it should work with alignment to 0.25, 1.8
START_SECTION([EXTRA] test_linear_res_align_6)
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_align(spec, 0.25, 1.8);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20 - 1.5 -1.5 ); // we loose 1.5 on the left

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 3); //+1.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 3+4);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 4.0+ 0.6);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 1.4 + 0.9 );
}
END_SECTION

// it should also work when we scale the m/z
START_SECTION([EXTRA] test_linear_res_align_scaling)
{
  MSSpectrum spec = input_spectrum;
  for (Size i = 0; i < spec.size(); i++)
  {
    spec[i].setMZ( spec[i].getMZ()*10 );
  }

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", 5.0);
  lr.setParameters(param);
	lr.raster_align(spec, -2.5, 12.5);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20 - 2.4 -0.6); // missing points 1.75 and 2.25 which have intensity 2.4 together

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 1.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 1.5+3);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 3+4);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 4.0); //+ 0.6);
}
END_SECTION
#endif

// it should work with ppm scaling
START_SECTION([EXTRA] test_linear_res_align_7)
{
  MSSpectrum spec = input_spectrum;

  // int = [3,6,8,2,1]
  // mz = [100, 101, 102, 103, 104]
  spec[0].setMZ(99 + 0.99/2.0);
  spec[1].setMZ(99.99 + 0.5);
  spec[2].setMZ(100.99 + 1.01/2.0);
  spec[3].setMZ(102 + 1.02 / 2.0);
  spec[4].setMZ(103.02 + 1.03 / 2.0);

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", 10000.0);
  param.setValue("ppm", "true");
  lr.setParameters(param);
	lr.raster_align(spec, 99, 105);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
    std::cout << spec[i] << std::endl;
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 1.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 4.4997); // 3 + 1.5
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 6.99911); // 3 + 4
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 5.0008); // 4 + 1
  TEST_REAL_SIMILAR(spec[5].getIntensity(), 0.500101);
}
END_SECTION

START_SECTION([EXTRA] test_linear_res_align_8)
{
  MSSpectrum spec = input_spectrum;

  // int = [3,6,8,2,1]
  // mz = [100, 101, 102, 103, 104]
  for (Size i=0; i<spec.size(); ++i)
  {
    spec[i].setMZ(100 + i);
  }

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", 10000.0);
  param.setValue("ppm", "true");
  lr.setParameters(param);
	lr.raster_align(spec, 99, 105);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
    std::cout << spec[i] << std::endl;
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[1].getIntensity(), 2.97);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 5.97);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 8.09725);
  TEST_REAL_SIMILAR(spec[4].getIntensity(), 2.01129);
  TEST_REAL_SIMILAR(spec[5].getIntensity(), 0.951471);
}
END_SECTION

#if 1

// also the interpolation should work
START_SECTION((template < typename PeakTypeIterator > void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator it, PeakTypeIterator resampled_end)))
{
  MSSpectrum spec = input_spectrum;
  MSSpectrum resampled;

  int i = 0;
  double start_pos = 0.25;
  double end_pos = 2.0;
  double spacing = 0.5;
  int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing + 1));
  resampled.resize(number_resampled_points);
  for (MSSpectrum::iterator it = resampled.begin(); it != resampled.end(); it++)
  {
      it->setMZ( start_pos + i*spacing);
      ++i;
  }

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_interpolate(spec.begin(), spec.end(), resampled.begin(), resampled.end() );

  spec = resampled;

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 4.5);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 7);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 5.5);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 1.25);
}
END_SECTION

START_SECTION(( template < typename PeakTypeIterator, typename ConstPeakTypeIterator > void raster(ConstPeakTypeIterator raw_it, ConstPeakTypeIterator raw_end, PeakTypeIterator resample_it, PeakTypeIterator resample_end)))
{

  MSSpectrum spec = input_spectrum;
  MSSpectrum resampled;

  int i = 0;
  double start_pos = 0;
  double end_pos = 2.25;
  double spacing = 0.75;
  int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing + 1));
  resampled.resize(number_resampled_points);
  for (MSSpectrum::iterator it = resampled.begin(); it != resampled.end(); it++)
  {
      it->setMZ( start_pos + i*spacing);
      ++i;
  }

  // A spacing of 0.75 will lead to a recalculation of intensities, each
  // resampled point gets intensities from raw data points that are at most +/-
  // spacing away.

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.75);
  lr.setParameters(param);
	lr.raster(spec.begin(), spec.end(), resampled.begin(), resampled.end() );

  spec = resampled;

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);

  TEST_REAL_SIMILAR(spec[0].getIntensity(), 3+2);
  TEST_REAL_SIMILAR(spec[1].getIntensity(), 4+2.0/3*8);
  TEST_REAL_SIMILAR(spec[2].getIntensity(), 1.0/3*8+2+1.0/3);
  TEST_REAL_SIMILAR(spec[3].getIntensity(), 2.0 / 3);
}
END_SECTION

// it should accept nonsense input values
START_SECTION([EXTRA] test_linear_res_align_input)
{
  MSSpectrum spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);

  lr.raster_align(spec, 2.25, 1.8);
  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 0);

  spec = input_spectrum;
  lr.raster_align(spec, 0.25, -1.8);
  sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 0);

  spec = input_spectrum;
  lr.raster_align(spec, 2.25, 5.8);
  sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 0);

  spec = input_spectrum;
  lr.raster_align(spec, -2.25, -2.0);
  sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 0);

}
END_SECTION

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


