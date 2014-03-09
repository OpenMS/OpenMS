// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// Rapid
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
// Median
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(SignalToNoiseEstimatorMedianRapidComparison, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Comparisons between NoiseEstimatorRapid and  NoiseEstimatorMedian //

// Compare the two NoiseEstimator -> there should be less than 50% difference
START_SECTION([EXTRA] compare )
{

  MSExperiment<> raw_data;
  double window_length = 20.0;

  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzML"),raw_data);
  MSSpectrum<> spec = raw_data[0];
  // copy spectrum to container
  std::vector<double> mz(spec.size()), intensity(spec.size());
  for (Size p = 0; p < spec.size(); ++p)
  {
    mz[p] = spec[p].getMZ();
    intensity[p] = spec[p].getIntensity();
  }

  SignalToNoiseEstimatorMedian< MSSpectrum < > > sne;  
	Param p;
	p.setValue("win_len", window_length);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(spec);

  SignalToNoiseEstimatorMedianRapid rapid_sne(window_length);
  SignalToNoiseEstimatorMedianRapid::NoiseEstimator rapid_ne = rapid_sne.estimateNoise(mz, intensity);

  // allow for a 50% difference between the two
  TOLERANCE_RELATIVE(1.5)
  MSSpectrum<>::iterator it = spec.begin();
  for (;it!=spec.end(); ++it)
  {
    double val1 = it->getIntensity() / rapid_ne.get_noise_value( it->getMZ() );
    double val2 = sne.getSignalToNoise(it);

    TEST_REAL_SIMILAR(val1, val2);
  }

}
END_SECTION

// Compare a file with rather a lot of noise -> there should be less than 20% difference
START_SECTION([EXTRA] compare_noisy )
{
  MSExperiment<> raw_data;
  double window_length = 20.0;

  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_5_long.mzML"),raw_data);
  MSSpectrum<> spec = raw_data[0];
  // copy spectrum to container
  std::vector<double> mz(spec.size()), intensity(spec.size());
  for (Size p = 0; p < spec.size(); ++p)
  {
    mz[p] = spec[p].getMZ();
    intensity[p] = spec[p].getIntensity();
  }

  SignalToNoiseEstimatorMedian< MSSpectrum < > > sne;  
	Param p;
	p.setValue("win_len", window_length);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(spec);

  SignalToNoiseEstimatorMedianRapid rapid_sne(window_length);
  SignalToNoiseEstimatorMedianRapid::NoiseEstimator rapid_ne = rapid_sne.estimateNoise(mz, intensity);

  TOLERANCE_RELATIVE(1.20)
  MSSpectrum<>::iterator it = spec.begin();
  for (;it!=spec.end(); ++it)
  {
    double val1 = it->getIntensity() / rapid_ne.get_noise_value( it->getMZ() );
    double val2 = sne.getSignalToNoise(it);

    TEST_REAL_SIMILAR(val1, val2);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


