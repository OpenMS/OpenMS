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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSignalToNoiseEstimator
  : public SignalToNoiseEstimator< >
{
  public:
  TestSignalToNoiseEstimator()
    : SignalToNoiseEstimator< >()
  {
  }

  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf)
  : SignalToNoiseEstimator< >(bpf)
  {
  }

  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
    if (&bpf==this) return *this;

    SignalToNoiseEstimator< >::operator=(bpf);

    return *this;
  }

  protected:

  void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_)
      throw() override
  {
    if (scan_first_ == scan_last_)
    {
      std::cout << "bla";
    }
    // do nothing here...
  }

};

START_TEST(SignalToNoiseEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = nullptr;
TestSignalToNoiseEstimator* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimator()))
	ptr = new TestSignalToNoiseEstimator();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((SignalToNoiseEstimator(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec.begin(), spec.end());
  TestSignalToNoiseEstimator sne_copy(sne);
	NOT_TESTABLE
END_SECTION


START_SECTION((SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec.begin(), spec.end());
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
	NOT_TESTABLE
END_SECTION


START_SECTION((virtual ~SignalToNoiseEstimator()))
	delete ptr;
END_SECTION


START_SECTION((virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec.begin(), spec.end());
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual void init(const Container& c)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum spec;
  sne.init(spec);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual double getSignalToNoise(const PeakIterator& data_point)))
  // hard to do without implementing computeSTN_ properly
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual double getSignalToNoise(const PeakType &data_point)))
  // hard to do without implementing computeSTN_ properly
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


