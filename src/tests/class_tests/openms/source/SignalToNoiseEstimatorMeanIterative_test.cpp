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
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

// #define DEBUG_TEST
// #undef  DEBUG_TEST

START_TEST(SignalToNoiseEstimatorMeanIterative, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMeanIterative< >* ptr = nullptr;
SignalToNoiseEstimatorMeanIterative< >* nullPointer = nullptr;
START_SECTION((SignalToNoiseEstimatorMeanIterative()))
        ptr = new SignalToNoiseEstimatorMeanIterative<>;
	TEST_NOT_EQUAL(ptr, nullPointer)
        SignalToNoiseEstimatorMeanIterative<> sne;
END_SECTION



START_SECTION((SignalToNoiseEstimatorMeanIterative& operator=(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2 = sne;
	NOT_TESTABLE
END_SECTION

START_SECTION((SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2(sne);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual ~SignalToNoiseEstimatorMeanIterative()))
  delete ptr;
END_SECTION


START_SECTION([EXTRA](virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)))

	TOLERANCE_ABSOLUTE(0.5)

  MSSpectrum raw_data;
  MSSpectrum::const_iterator it;
  DTAFile dta_file;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimator_test.dta"), raw_data);
  
    
  SignalToNoiseEstimatorMeanIterative<> sne;  
	Param p;
	p.setValue("win_len", 40.1);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(raw_data.begin(),raw_data.end());

  MSSpectrum stn_data;
  
#ifdef DEBUG_TEST
  MSSpectrum stn_data__;
#endif
  
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.out"), stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_SIMILAR (stn_data[i].getIntensity(), sne.getSignalToNoise(it));
#ifdef DEBUG_TEST    
    Peak1D peak = (*it);
    peak.setIntensity(stn_data[i].getIntensity() / sne.getSignalToNoise(it));
    stn_data__.push_back(peak);
#endif    
    ++i;
  }
  
#ifdef DEBUG_TEST
  dta_file.store(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.debug"), stn_data__);
#endif  
  
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


