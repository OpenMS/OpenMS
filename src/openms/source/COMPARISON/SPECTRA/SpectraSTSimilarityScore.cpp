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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <cmath>
using namespace std;

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
    double score(0);
    BinnedSpectrum bin1(1, 1, s1);
    BinnedSpectrum bin2(1, 1, s2);

    //normalize bins

    //magnitude of the spectral vector
    float magnitude1(0);
    float magnitude2(0);
    for (SparseVector<float>::SparseVectorIterator iter1 = bin1.getBins().begin();
          iter1 < bin1.getBins().end(); ++iter1)
    {
      magnitude1 += pow((double) * iter1, 2);
    }
    magnitude1 = sqrt(magnitude1);

    //normalize bins of bin1
    for (SparseVector<float>::SparseVectorIterator iter1 = bin1.getBins().begin();
          iter1 < bin1.getBins().end(); ++iter1)
    {
      *iter1 = (float) * iter1 / magnitude1;
    }

    for (SparseVector<float>::SparseVectorIterator iter2 = bin2.getBins().begin();
          iter2 < bin2.getBins().end(); ++iter2)
    {
      magnitude2 += pow((double) * iter2, 2);
    }
    magnitude2 = sqrt(magnitude2);

    //normalize bins of bin2
    for (SparseVector<float>::SparseVectorIterator iter2 = bin2.getBins().begin();
          iter2 < bin2.getBins().end(); ++iter2)
    {
      *iter2 = (float) * iter2 / magnitude2;
    }

    Size shared_bins = min(bin1.getBinNumber(), bin2.getBinNumber());
    for (Size s = 0; s < shared_bins; ++s)
    {
      if ((double)bin1.getBins()[s] > 0.0 && (double)bin2.getBins()[s] > 0.0)
      {
        score += ((double)bin1.getBins()[s] * (double)bin2.getBins()[s]);
      }
    }

    return score;

  }

  double SpectraSTSimilarityScore::operator()(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2) const
  {
    double score(0);

    Size shared_bins = min(bin1.getBinNumber(), bin2.getBinNumber());
    for (Size s = 0; s < shared_bins; ++s)
    {
      if (bin1.getBins()[s] > 0 && bin2.getBins()[s] > 0)
      {
        score += (bin1.getBins()[s] * bin2.getBins()[s]);
      }
    }

    return score;
  }

  bool SpectraSTSimilarityScore::preprocess(PeakSpectrum & spec,
                                            float remove_peak_intensity_threshold,
                                            UInt cut_peaks_below,
                                            Size min_peak_number,
                                            Size max_peak_number)
  {
    spec.sortByIntensity(true);
    double min_high_intensity = 0;
    if (!spec.empty())
    {
      min_high_intensity = (1 / cut_peaks_below) * spec[0].getIntensity();
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
    BinnedSpectrum bin(1, 1, spec);
    float magnitude(0);
    for (SparseVector<float>::SparseVectorIterator iter = bin.getBins().begin(); iter < bin.getBins().end(); ++iter)
    {
      magnitude += pow((double) * iter, 2);
    }
    magnitude = sqrt(magnitude);
    //normalize bins
    for (SparseVector<float>::SparseVectorIterator iter = bin.getBins().begin(); iter < bin.getBins().end(); ++iter)
    {
      *iter = (float) * iter / magnitude;
    }
    return bin;
  }

  double SpectraSTSimilarityScore::dot_bias(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2, double dot_product) const
  {
    double numerator(0);

    Size shared_bins = min(bin1.getBinNumber(), bin2.getBinNumber());
    for (Size s = 0; s < shared_bins; ++s)
    {
      if (bin1.getBins()[s] > 0 && bin2.getBins()[s] > 0)
      {
        numerator += (pow(bin1.getBins()[s], 2) * pow(bin2.getBins()[s], 2));
      }
    }
    numerator = sqrt(numerator);

    if (dot_product)
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
