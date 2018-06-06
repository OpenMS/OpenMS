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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/SMOOTHING/GaussFilterAlgorithm.h>

///////////////////////////

START_TEST(GaussFilterAlgorithm<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

GaussFilterAlgorithm* dgauss_ptr = nullptr;
GaussFilterAlgorithm* dgauss_nullPointer = nullptr;

START_SECTION((GaussFilterAlgorithm()))
  dgauss_ptr = new GaussFilterAlgorithm;
  TEST_NOT_EQUAL(dgauss_ptr, dgauss_nullPointer)
END_SECTION

START_SECTION((virtual ~GaussFilterAlgorithm()))
    delete dgauss_ptr;
END_SECTION

START_SECTION((void initialize(double gaussian_width, double spacing, double ppm_tolerance, bool use_ppm_tolerance)))
  // We cannot really test that the variables are correctly set since we dont
  // have access to them.
  dgauss_ptr = new GaussFilterAlgorithm;
  dgauss_ptr->initialize(0.5, 0.5, 10, false);
  dgauss_ptr->initialize(0.5, 0.5, 10, true);
  TEST_NOT_EQUAL(dgauss_ptr, dgauss_nullPointer)
  delete dgauss_ptr;
END_SECTION

START_SECTION((template <typename ConstIterT, typename IterT> bool filter(ConstIterT mz_in_start, ConstIterT mz_in_end, ConstIterT int_in_start, IterT mz_out, IterT int_out)))
  std::vector<double> mz;
  std::vector<double> intensities;
  std::vector<double> mz_out(5);
  std::vector<double> intensities_out(5);

  for (Size i=0; i<5; ++i)
  {
    intensities.push_back(1.0f);
    mz.push_back(500.0+0.2*i);
  }

  GaussFilterAlgorithm gauss;
  gauss.initialize(1.0 * 8 /* gaussian_width */, 0.01 /* spacing */, 10.0 /* ppm_tolerance */, false /* use_ppm_tolerance */);
  gauss.filter(mz.begin(), mz.end(), intensities.begin(), mz_out.begin(), intensities_out.begin());

  std::vector<double>::iterator it=intensities_out.begin();
  TEST_REAL_SIMILAR(*it,1.0)
  ++it;
  TEST_REAL_SIMILAR(*it,1.0)
  ++it;
  TEST_REAL_SIMILAR(*it,1.0)
  ++it;
  TEST_REAL_SIMILAR(*it,1.0)
  ++it;
  TEST_REAL_SIMILAR(*it,1.0)
END_SECTION 

START_SECTION((bool filter(OpenMS::Interfaces::SpectrumPtr spectrum)))

  OpenMS::Interfaces::SpectrumPtr spectrum(new OpenMS::Interfaces::Spectrum);

  spectrum->getMZArray()->data.resize(9);
  spectrum->getIntensityArray()->data.resize(9);

  for (Size i=0; i<9; ++i)
  {
    spectrum->getIntensityArray()->data[i] = 0.0f;
    spectrum->getMZArray()->data[i] = 500.0+0.03*i; 
    if (i==3)
    {
      spectrum->getIntensityArray()->data[i] = 1.0f;
    }
    if (i==4)
    {
      spectrum->getIntensityArray()->data[i] = 0.8f;
    }
    if (i==5)
    {
      spectrum->getIntensityArray()->data[i] = 1.2f;
    }
  }

  TOLERANCE_ABSOLUTE(0.01)

  GaussFilterAlgorithm gauss;
  TEST_EQUAL(spectrum->getIntensityArray()->data.size(), 9) 
  gauss.initialize(0.2, 0.01, 1.0, false);
  gauss.filter(spectrum);

  TEST_EQUAL(spectrum->getIntensityArray()->data.size(), 9) 
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[0],0.000734827)  
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[1],0.0543746)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[2],0.298025)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[3],0.707691)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[4],0.8963)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[5],0.799397)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[6],0.352416)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[7],0.065132)
  TEST_REAL_SIMILAR(spectrum->getIntensityArray()->data[8],0.000881793)
END_SECTION 


START_SECTION((bool filter(OpenMS::Interfaces::ChromatogramPtr chromatogram)))

  OpenMS::Interfaces::ChromatogramPtr chromatogram(new OpenMS::Interfaces::Chromatogram);

  chromatogram->getTimeArray()->data.resize(9);
  chromatogram->getIntensityArray()->data.resize(9);

  for (Size i=0; i<9; ++i)
  {
    chromatogram->getIntensityArray()->data[i] = 0.0f;
    chromatogram->getTimeArray()->data[i] = 500.0+0.03*i; 
    if (i==3)
    {
      chromatogram->getIntensityArray()->data[i] = 1.0f;
    }
    if (i==4)
    {
      chromatogram->getIntensityArray()->data[i] = 0.8f;
    }
    if (i==5)
    {
      chromatogram->getIntensityArray()->data[i] = 1.2f;
    }
  }

  TOLERANCE_ABSOLUTE(0.01)

  GaussFilterAlgorithm gauss;
  TEST_EQUAL(chromatogram->getIntensityArray()->data.size(), 9) 
  gauss.initialize(0.2, 0.01, 1.0, false);
  gauss.filter(chromatogram);

  TEST_EQUAL(chromatogram->getIntensityArray()->data.size(), 9) 
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[0],0.000734827)  
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[1],0.0543746)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[2],0.298025)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[3],0.707691)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[4],0.8963)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[5],0.799397)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[6],0.352416)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[7],0.065132)
  TEST_REAL_SIMILAR(chromatogram->getIntensityArray()->data[8],0.000881793)
END_SECTION 

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
