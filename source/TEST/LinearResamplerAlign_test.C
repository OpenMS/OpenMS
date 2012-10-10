// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>

using namespace OpenMS;
using namespace std;

template <template <typename> class SpectrumT, typename PeakT>
void check_results(SpectrumT<PeakT> spec)
{
  DoubleReal sum = 0.0;
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

///////////////////////////

START_TEST(LinearResamplerAlign, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSSpectrum< Peak1D > input_spectrum;
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

START_SECTION((template < template< typename > class MSSpectrum, typename PeakType > void raster(MSSpectrum< PeakType > &spectrum)))
{

  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing", default_spacing);
  lr.setParameters(param);
	lr.raster(spec);

  DoubleReal sum = 0.0;
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
  MSChromatogram< Peak1D > spec;
  spec.resize(5);
  spec[0].setMZ(0);
  spec[0].setIntensity(3.0f);
  spec[1].setMZ(0.5);
  spec[1].setIntensity(6.0f);
  spec[2].setMZ(1.);
  spec[2].setIntensity(8.0f);
  spec[3].setMZ(1.6);
  spec[3].setIntensity(2.0f);
  spec[4].setMZ(1.8);
  spec[4].setIntensity(1.0f);

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",default_spacing);
  lr.setParameters(param);
	lr.raster(spec);

  check_results(spec);
}
END_SECTION

// it should work with alignment to 0, 1.8 and give the same result
START_SECTION((template < template< typename > class MSSpectrum, typename PeakType > void raster_align(MSSpectrum< PeakType > &spectrum, double start_pos, double end_pos)))
{
  MSSpectrum< Peak1D > spec = input_spectrum;

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
  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_align(spec, -0.25, 1.8);

  DoubleReal sum = 0.0;
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
  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.75);
  lr.setParameters(param);
	lr.raster_align(spec, -2.25, 1.8);

  DoubleReal sum = 0.0;
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
  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_align(spec, -0.25, 1.25);

  DoubleReal sum = 0.0;
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
  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
	lr.raster_align(spec, 0.25, 1.8);

  DoubleReal sum = 0.0;
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
  MSSpectrum< Peak1D > spec = input_spectrum;
  for (Size i = 0; i < spec.size(); i++)
  {
    spec[i].setMZ( spec[i].getMZ()*10 );
  }

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",5.0);
  lr.setParameters(param);
	lr.raster_align(spec, -2.5, 12.5);

  DoubleReal sum = 0.0;
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

// also the interpolation should work
START_SECTION((template < typename PeakTypeIterator > void raster_interpolate(PeakTypeIterator raw_it, PeakTypeIterator raw_end, PeakTypeIterator it, PeakTypeIterator resampled_end)))
{
  MSSpectrum< Peak1D > spec = input_spectrum;
  MSSpectrum< Peak1D > resampled;

  int i = 0;
  double start_pos = 0.25;
  double end_pos = 2.0;
  double spacing = 0.5;
  int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing + 1));
  resampled.resize(number_resampled_points);
  for (MSSpectrum<Peak1D>::iterator it = resampled.begin(); it != resampled.end(); it++)
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

  DoubleReal sum = 0.0;
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

  MSSpectrum< Peak1D > spec = input_spectrum;
  MSSpectrum< Peak1D > resampled;

  int i = 0;
  double start_pos = 0;
  double end_pos = 2.25;
  double spacing = 0.75;
  int number_resampled_points = (int)(ceil((end_pos -start_pos) / spacing + 1));
  resampled.resize(number_resampled_points);
  for (MSSpectrum<Peak1D>::iterator it = resampled.begin(); it != resampled.end(); it++)
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

  DoubleReal sum = 0.0;
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
  MSSpectrum< Peak1D > spec = input_spectrum;

  LinearResamplerAlign lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);

	lr.raster_align(spec, 2.25, 1.8);
  DoubleReal sum = 0.0;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


