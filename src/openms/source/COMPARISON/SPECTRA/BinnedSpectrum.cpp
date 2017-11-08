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
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

using namespace std;

namespace OpenMS
{

  BinnedSpectrum::BinnedSpectrum(const float size, const UInt spread, const PeakSpectrum & ps) :
    bin_spread_(spread), 
    bin_size_(size), 
    bins_()
  {
    precursors_ = ps.getPrecursors();
    binSpectrum_(ps);
  }

  BinnedSpectrum::~BinnedSpectrum()
  {
  }

  const SparseVector<float>& BinnedSpectrum::getBins() const
  {
    return bins_;
  }

  SparseVector<float>& BinnedSpectrum::getBins()
  {
    return bins_;
  }

  const std::vector<Precursor>& BinnedSpectrum::getPrecursors() const
  {
    return precursors_;
  }

  std::vector<Precursor>& BinnedSpectrum::getPrecursors()
  {
    return precursors_;
  }

  void BinnedSpectrum::binSpectrum_(const PeakSpectrum& ps)
  {
    OPENMS_PRECONDITION(ps.isSorted(), "Spectrum needs to be sorted by m/z.");
    
    if (ps.empty()) { return; }

    const size_t highest_index = getBinIndex(ps.back().getMZ()) + bin_spread_;
    bins_ = SparseVector<float>(highest_index + 1, 0, 0);

    // put all peaks into bins
    for (auto const & p : ps)
    {
      // e.g.: bin_size_ = 1.5: first bin covers range [0, 1.5) so peak at 1.5 falls in second bin (index 1)
      const size_t idx = getBinIndex(p.getMZ());

      // add peak to corresponding bin
      bins_[idx] = bins_.at(idx) + p.getIntensity();

      // add peak to neighboring bins
      for (Size j = 0; j < bin_spread_; ++j)
      {
        bins_[idx + j + 1] = bins_.at(idx + j + 1) + p.getIntensity();
        
        // prevent spreading over left boundaries
        if (static_cast<int>(idx - j - 1) >= 0)
        {
          bins_[idx - j - 1] = bins_.at(idx - j - 1) + p.getIntensity();
        }
      }
    }
  }

  // static
  bool BinnedSpectrum::isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b)
  {
    // check if bin size and spread are equal
    return std::tie(a.bin_size_, a.bin_spread_) 
        == std::tie(b.bin_size_, b.bin_spread_);
  }

  bool BinnedSpectrum::operator!=(const BinnedSpectrum& rhs) const
  {
    return !(operator==(rhs));
  }

}

