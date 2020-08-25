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
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

using namespace std;

namespace OpenMS
{
  // Empty vector initialized to maximum supported dimensionality.
  const BinnedSpectrum::SparseVectorType BinnedSpectrum::EmptySparseVector(numeric_limits<BinnedSpectrum::SparseVectorIndexType>::max());

  BinnedSpectrum::BinnedSpectrum(const PeakSpectrum& ps, float size, bool unit_ppm, UInt spread, float offset) :
    bin_spread_(spread), 
    bin_size_(size),
    unit_ppm_(unit_ppm),
    offset_(offset),
    bins_()
  {
    precursors_ = ps.getPrecursors();
    binSpectrum_(ps);
  }

  BinnedSpectrum::~BinnedSpectrum()
  {
  }

  const BinnedSpectrum::SparseVectorType& BinnedSpectrum::getBins() const
  {
    return bins_;
  }

  BinnedSpectrum::SparseVectorType& BinnedSpectrum::getBins()
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

    bins_ = EmptySparseVector;

    for (auto const & p : ps)
    {
      // if bin size is in relative units (ppm), check if minimum value is >= 1 (otherwise we might get numerical problems with the negative log)
      OPENMS_PRECONDITION(!unit_ppm_ || p.getMZ() >= BinnedSpectrum::MIN_MZ_, "Spectrum with relative bin size contains peaks with m/z < 1");

      // e.g.: bin_size_ = 1.5: first bin covers range [0, 1.5) so peak at 1.5 falls in second bin (index 1)
      const size_t idx = getBinIndex(p.getMZ());

      // add peak to corresponding bin
      bins_.coeffRef(idx) += p.getIntensity();

      // add peak to neighboring bins
      for (Size j = 0; j < bin_spread_; ++j)
      {
         bins_.coeffRef(idx + j + 1) +=  p.getIntensity();
        
        // prevent spreading over left boundaries
        if (static_cast<int>(idx - j - 1) >= 0)
        {
          bins_.coeffRef(idx - j - 1) += p.getIntensity();
        }
      }
    }
  }

  bool BinnedSpectrum::operator==(const BinnedSpectrum& rhs) const
  {
    // first compare bin layout and precursors
    if (std::tie(unit_ppm_, bin_size_, bin_spread_, precursors_)
        != std::tie(rhs.unit_ppm_, rhs.bin_size_, rhs.bin_spread_, rhs.precursors_)) 
    {
      return false; 
    }

    // efficient look-up of number of non-zero entries, so we use this as-well
    if (bins_.nonZeros() != rhs.bins_.nonZeros()) { return false; }

    // test non-sparse (non-zero) elements for equality
    SparseVectorIteratorType it(bins_);
    SparseVectorIteratorType rhs_it(rhs.bins_);  
    while (it)      
    {
      if (it.index() != rhs_it.index()
       || it.value() != rhs_it.value()) { return false; }
      ++it;
      ++rhs_it;
    }
    return true;
  }

  // static
  bool BinnedSpectrum::isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b)
  {
    // check if bin size, offset (and unit) are equal
    return std::tie(a.unit_ppm_, a.bin_size_, a.offset_) 
        == std::tie(b.unit_ppm_, b.bin_size_, b.offset_);
  }

  bool BinnedSpectrum::operator!=(const BinnedSpectrum& rhs) const
  {
    return !(operator==(rhs));
  }

}

