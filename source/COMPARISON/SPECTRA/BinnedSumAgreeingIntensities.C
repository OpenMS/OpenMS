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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>

using namespace std;

namespace OpenMS
{
  BinnedSumAgreeingIntensities::BinnedSumAgreeingIntensities() :
    BinnedSpectrumCompareFunctor()
  {
    setName(BinnedSumAgreeingIntensities::getProductName());
    defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
    defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
    defaultsToParam_();
  }

  BinnedSumAgreeingIntensities::BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities & source) :
    BinnedSpectrumCompareFunctor(source)
  {
  }

  BinnedSumAgreeingIntensities::~BinnedSumAgreeingIntensities()
  {
  }

  BinnedSumAgreeingIntensities & BinnedSumAgreeingIntensities::operator=(const BinnedSumAgreeingIntensities & source)
  {
    if (this != &source)
    {
      BinnedSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum & spec1, const BinnedSpectrum & spec2) const
  {
    // avoid crash while comparing
    if (!spec1.checkCompliance(spec2))
    {
      throw IncompatibleBinning(__FILE__, __LINE__, __PRETTY_FUNCTION__, "");
    }

    // shortcut similarity calculation by comparing PrecursorPeaks (PrecursorPeaks more than delta away from each other are supposed to be from another peptide)
    DoubleReal pre_mz1 = 0.0;
    if (!spec1.getPrecursors().empty())
      pre_mz1 = spec1.getPrecursors()[0].getMZ();
    DoubleReal pre_mz2 = 0.0;
    if (!spec1.getPrecursors().empty())
      pre_mz2 = spec2.getPrecursors()[0].getMZ();
    if (fabs(pre_mz1 - pre_mz2) > (double)param_.getValue("precursor_mass_tolerance"))
    {
      return 0;
    }

    double score(0), sharedBins(min(spec1.getBinNumber(), spec2.getBinNumber())), sum1(0), sum2(0), summax(0);

    // all bins at equal position and similar intensities contribute positively to score
    for (Size i = 0; i < sharedBins; ++i)
    {
      sum1 += spec1.getBins()[i];
      sum2 += spec2.getBins()[i];
      summax += max((float)0, ((spec1.getBins()[i] + spec2.getBins()[i]) / 2) - fabs(spec1.getBins()[i] - spec2.getBins()[i]));
    }

    // resulting score normalized to interval [0,1]
    score = summax * (2 / (sum1 + sum2));

    return score;

  }

}
