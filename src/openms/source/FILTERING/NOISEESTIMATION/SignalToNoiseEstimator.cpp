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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <random>

using namespace std;

namespace OpenMS
{

  float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans, const double percentile)
  {
    vector<Size> spec_indices;
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (exp[i].getMSLevel() == ms_level && !exp[i].empty())
      {
        spec_indices.push_back(i);
      }
    }

    if (spec_indices.empty()) return 0.0f;

    std::default_random_engine generator(time(nullptr));
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    float noise = 0.0;
    UInt count = 0;
    vector<float> tmp;
    while (count++ < n_scans)
    {
      UInt scan = (UInt)(distribution(generator) * (spec_indices.size() - 1));
      tmp.clear();
      for (const auto& peak : exp[scan])
      {
        tmp.push_back(peak.getIntensity());
      }
      Size idx = tmp.size() * percentile / 100.0;
      std::nth_element(tmp.begin(), tmp.begin() + idx, tmp.end());
      noise += tmp[idx];
    }
    return noise / n_scans;
  }

}
