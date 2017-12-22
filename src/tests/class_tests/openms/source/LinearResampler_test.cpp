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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>


#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////

START_TEST(LinearResampler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

LinearResampler* lr_ptr = nullptr;
LinearResampler* lr_nullPointer = nullptr;

START_SECTION((LinearResampler()))
  lr_ptr = new LinearResampler;
  TEST_NOT_EQUAL(lr_ptr, lr_nullPointer);
END_SECTION

START_SECTION((~LinearResampler()))
  delete lr_ptr;
END_SECTION

START_SECTION((template<typename PeakType> void raster(MSSpectrum& spectrum)))
{
  MSSpectrum spec;
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

  LinearResampler lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
  lr.raster(spec);

  double sum = 0.0;
  for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);
}

/////////// test raster with a spacing of 0.75
{
  MSSpectrum spec;
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

  // A spacing of 0.75 will lead to a recalculation of intensities, each
  // resampled point gets intensities from raw data points that are at most +/-
  // spacing away.
  LinearResampler lr;
  Param param;
  param.setValue("spacing",0.75);
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

START_SECTION(( template <typename PeakType > void rasterExperiment(MSExperiment<PeakType>& exp)))
{
  MSSpectrum spec;
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

  PeakMap exp;
  exp.addSpectrum(spec);
  exp.addSpectrum(spec);

  LinearResampler lr;
  Param param;
  param.setValue("spacing",0.5);
  lr.setParameters(param);
  lr.rasterExperiment(exp);


  for (Size s=0; s<exp.size(); ++s)
  {
    double sum = 0.0;
    for (Size i=0; i<exp[s].size(); ++i)
    {
      sum += exp[s][i].getIntensity();
    }
    TEST_REAL_SIMILAR(sum, 20);
  }

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
