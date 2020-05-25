// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransformNumIntegration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransformNumIntegration* ptr = nullptr;
ContinuousWaveletTransformNumIntegration* nullPointer = nullptr;
START_SECTION((ContinuousWaveletTransformNumIntegration()))
  ptr = new ContinuousWaveletTransformNumIntegration();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ContinuousWaveletTransformNumIntegration()))
  delete ptr;
END_SECTION

START_SECTION((virtual void init(double scale, double spacing)))
  ContinuousWaveletTransformNumIntegration transformer;
  float scale = 0.5f;
  float spacing = 0.1f;
  
  transformer.init(scale,spacing);
  TEST_REAL_SIMILAR(transformer.getWavelet()[0],1.)
  TEST_REAL_SIMILAR(transformer.getScale(),scale)
  TEST_REAL_SIMILAR(transformer.getSpacing(),spacing)
END_SECTION

START_SECTION((template <typename InputPeakIterator> void transform(InputPeakIterator begin_input, InputPeakIterator end_input, float resolution, unsigned int zeros=0)))
  ContinuousWaveletTransformNumIntegration transformer;
  float scale = 0.5f;
  float spacing = 0.1f;
  
  transformer.init(scale,spacing);
  std::vector<Peak1D > raw_data(9);
  raw_data[4].setIntensity(1.0f);
  transformer.transform(raw_data.begin(),raw_data.end(),1.);
  TEST_REAL_SIMILAR(transformer[4],0)
  TEST_REAL_SIMILAR(transformer.getWavelet()[0],1.)
  TEST_REAL_SIMILAR(transformer.getScale(),scale)
  TEST_REAL_SIMILAR(transformer.getSpacing(),spacing)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



