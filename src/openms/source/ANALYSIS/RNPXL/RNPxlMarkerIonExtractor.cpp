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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

using namespace std;

namespace OpenMS
{


RNPxlMarkerIonExtractor::MarkerIonsType RNPxlMarkerIonExtractor::extractMarkerIons(const PeakSpectrum& s, const double marker_tolerance)
{
  MarkerIonsType marker_ions;
  marker_ions["A"].push_back(make_pair(136.06231, 0.0));
  marker_ions["A"].push_back(make_pair(330.06033, 0.0));
  marker_ions["C"].push_back(make_pair(112.05108, 0.0));
  marker_ions["C"].push_back(make_pair(306.04910, 0.0));
  marker_ions["G"].push_back(make_pair(152.05723, 0.0));
  marker_ions["G"].push_back(make_pair(346.05525, 0.0));
  marker_ions["U"].push_back(make_pair(113.03509, 0.0));
  marker_ions["U"].push_back(make_pair(307.03311, 0.0));

  PeakSpectrum spec(s);
  Normalizer normalizer;
  normalizer.filterSpectrum(spec);
  spec.sortByPosition();

  // for each nucleotide with marker ions
  for (MarkerIonsType::iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
  {
    // for each marker ion of the current nucleotide
    for (Size i = 0; i != it->second.size(); ++i)
    {
      double mz = it->second[i].first;
      double max_intensity = 0;
      for (PeakSpectrum::ConstIterator sit = spec.begin(); sit != spec.end(); ++sit)
      {
        if (sit->getMZ() + marker_tolerance < mz)
        {
          continue;
        }
        if (mz < sit->getMZ() - marker_tolerance)
        {
          break;
        }
        if (fabs(mz - sit->getMZ()) < marker_tolerance)
        {
          if (max_intensity < sit->getIntensity())
          {
            max_intensity = sit->getIntensity();
          }
        }
      }
      it->second[i].second = max_intensity;
    }
  }
  return marker_ions;
}
}

