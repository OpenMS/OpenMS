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

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransform, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransform* ptr = nullptr;
ContinuousWaveletTransform* nullPointer = nullptr;
START_SECTION((ContinuousWaveletTransform()))
  ptr = new ContinuousWaveletTransform();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ContinuousWaveletTransform()))
  delete ptr;
END_SECTION


START_SECTION((std::vector<Peak1D >& getSignal()))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(2);
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
END_SECTION

START_SECTION((const std::vector<Peak1D >& getSignal() const))
 ContinuousWaveletTransform cwt;
 
 TEST_EQUAL(cwt.getSignal().empty(), true)
END_SECTION

START_SECTION((double getScale() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0)
END_SECTION

START_SECTION((double getSpacing() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0)
END_SECTION

START_SECTION((double operator[](unsigned int i) const ))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(1);
  cwt.getSignal() = signal;
  
  ContinuousWaveletTransform const cwt_const(cwt);
  
  TEST_REAL_SIMILAR(cwt_const[0],0)
END_SECTION

START_SECTION((SignedSize getLeftPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 0)
END_SECTION

START_SECTION((SignedSize getRightPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 0)
END_SECTION

START_SECTION((SignedSize getSignalLength() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getSignalLength(), 0)
END_SECTION

START_SECTION((int getSize() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getSize(), 0)
END_SECTION

START_SECTION((const std::vector<double>& getWavelet() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getWavelet().size(), 0)
END_SECTION

START_SECTION((double& getScale()))
  ContinuousWaveletTransform cwt;
  cwt.getScale() = 0.2;
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0.2)
END_SECTION

START_SECTION((double& getSpacing()))
  ContinuousWaveletTransform cwt;
  cwt.getSpacing() = 0.2;
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0.2)
END_SECTION

START_SECTION((double operator[](unsigned int i)))
  std::vector<Peak1D > signal;
  Peak1D rp;
  rp.setIntensity(100.0f);
  signal.push_back(rp);
  
  ContinuousWaveletTransform cwt;
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt[0],100)
END_SECTION

START_SECTION((SignedSize& getLeftPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getLeftPaddingIndex() = 2;
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 2)
END_SECTION

START_SECTION((SignedSize& getRightPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getRightPaddingIndex() = 2;
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
END_SECTION

START_SECTION((SignedSize& getSignalLength()))
  ContinuousWaveletTransform cwt;
  cwt.getSignalLength() = 2;
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
END_SECTION

START_SECTION((std::vector<double>& getWavelet()))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
END_SECTION



START_SECTION((virtual void init(double scale, double spacing)))
  ContinuousWaveletTransform cwt;
  double scale = 0.2;
  double spacing = 2.3;
  cwt.init(scale,spacing);
  
  TEST_REAL_SIMILAR(cwt.getSpacing(),spacing)
  TEST_REAL_SIMILAR(cwt.getScale(),scale)
END_SECTION

START_SECTION((void setLeftPaddingIndex(const SignedSize end_left_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setLeftPaddingIndex(2);
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 2)
END_SECTION

START_SECTION((void setRightPaddingIndex(const SignedSize begin_right_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setRightPaddingIndex(2);
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
END_SECTION

START_SECTION((void setScale(double scale)))
  ContinuousWaveletTransform cwt;
  cwt.setScale(0.2);
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0.2)
END_SECTION

START_SECTION((void setSignal(const std::vector<Peak1D >& signal)))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(2);
  cwt.setSignal(signal);
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
END_SECTION

START_SECTION((void setSignalLength(const SignedSize signal_length)))
  ContinuousWaveletTransform cwt;
  cwt.setSignalLength(2);
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
END_SECTION

START_SECTION((void setSpacing(double spacing)))
  ContinuousWaveletTransform cwt;
  cwt.setSpacing(0.2);
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0.2)
END_SECTION

START_SECTION((void setWavelet(const std::vector<double>& wavelet)))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



