// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>

using namespace std;
using namespace Eigen;

namespace OpenMS
{
  SpectraSTSimilarityScore::SpectraSTSimilarityScore() :
    PeakSpectrumCompareFunctor()
  {
    setName(SpectraSTSimilarityScore::getProductName());
  }

  SpectraSTSimilarityScore::SpectraSTSimilarityScore(const SpectraSTSimilarityScore & source) :
    PeakSpectrumCompareFunctor(source)
  {
  }

  SpectraSTSimilarityScore::~SpectraSTSimilarityScore()
  {
  }

  SpectraSTSimilarityScore & SpectraSTSimilarityScore::operator=(const SpectraSTSimilarityScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double SpectraSTSimilarityScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double SpectraSTSimilarityScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    // TODO: check if this operator makes sense (as it doesn't allow to fine tune resolution)
    BinnedSpectrum bin1(s1, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
    BinnedSpectrum bin2(s2, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);

    // normalized dot product
    bin1.getBins() /= bin1.getBins().norm();
    bin2.getBins() /= bin2.getBins().norm();
    return bin1.getBins().dot(bin2.getBins());
  }

  double SpectraSTSimilarityScore::operator()(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2) const
  {
    return bin1.getBins().dot(bin2.getBins());
  }

  bool SpectraSTSimilarityScore::preprocess(PeakSpectrum & spec,
                                            float remove_peak_intensity_threshold,
                                            UInt cut_peaks_below,
                                            Size min_peak_number,
                                            Size max_peak_number)
  {
    double min_high_intensity = 0.;
    if (!spec.empty())
    {
      double max_el = std::max_element(spec.begin(),spec.end(),Peak1D::IntensityLess())->getIntensity();
      min_high_intensity = (1.0 / cut_peaks_below) * max_el;
    }

    spec.sortByPosition();

    PeakSpectrum tmp;
    Size s = 0;
    for (PeakSpectrum::iterator k = spec.begin(); k < spec.end() && s < max_peak_number; ++k, ++s)
    {
      Peak1D peak;
      if (k->getIntensity() >  remove_peak_intensity_threshold && k->getIntensity() > min_high_intensity)
      {
        peak.setIntensity(sqrt(k->getIntensity()));
        peak.setMZ(k->getMZ());
        peak.setPosition(k->getPosition());
        tmp.push_back(peak);
      }

    }
    spec = tmp;
    //if not enough peaks in the spectrum pass that one out
    return spec.size() >= min_peak_number;
  }

  BinnedSpectrum SpectraSTSimilarityScore::transform(const PeakSpectrum & spec)
  {
    // TODO: resolution seems rather low. Check with current original implementations.
    BinnedSpectrum bin(spec, 1, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
    bin.getBins() /= bin.getBins().norm();
    return bin;
  }

  double SpectraSTSimilarityScore::dot_bias(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2, double dot_product) const
  {
    double numerator = (bin1.getBins().cwiseProduct(bin2.getBins())).norm();
    
    if (dot_product != 0)
    {
      return (double)numerator / dot_product;
    }
    else
    {
      return (double)numerator / (*this)(bin1, bin2);
    }
  }

  double SpectraSTSimilarityScore::delta_D(double top_hit, double runner_up)
  {
    if (top_hit == 0)
    {
      throw Exception::DivisionByZero(__FILE__, __LINE__, __FUNCTION__);
    }
    else
    {
      return (double)(top_hit - runner_up) / top_hit;
    }
  }

  double SpectraSTSimilarityScore::compute_F(double dot_product, double delta_D, double dot_bias)
  {
    double b(0);
    if (dot_bias < 0.1 || (0.35 < dot_bias && dot_bias <= 0.4))
    {
      b = 0.12;
    }
    else if (0.4 < dot_bias && dot_bias <= 0.45)
    {
      b = 0.18;
    }
    else if (dot_bias > 0.45)
    {
      b  = 0.24;
    }
    return 0.6 * dot_product + 0.4 * delta_D - b;
  }

}
