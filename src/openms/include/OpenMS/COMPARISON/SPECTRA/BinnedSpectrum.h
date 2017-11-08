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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>

namespace OpenMS
{

  /**
    @brief This is a binned representation of a PeakSpectrum

    @param sz the size of the bins and
    @param sp number of neighboring bins to both sides affected by a peak contribution
    @param ps the PeakSpectrum, used to calculate the bins

    sz denotes the size of a bin in @p Th, thereby deciding the number of bins (all of size sz) the spectrum is discretized to.
    Each bin will represent a certain @p Th range and the peaks will be put in the respective bins and sum up inside.
    sp denotes the number of neighboring bins to the left and the number of neighboring bins to the right a peak is also added to.
    E.g. a BinnedSpectrum with binsize of 0.5 @p Th will have a peak at 100 @p Th in bin no. 200, a peak at 100.1 @p Th will be in bin no. 201.
    If the binspread is 1, the peak at 100 Th will be added to bin no. 199, 200 and 201.
    If the binspread is 2, the peak at 100 @p Th will also be added to bin no. 198 and 202, and so on.

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSpectrum
  {

private:
    /// the spread to left or right
    UInt bin_spread_;

    /// the size of each bin
    float bin_size_;

    /// bins
    SparseVector<float> bins_;

    void binSpectrum_(const PeakSpectrum& ps);

    std::vector<Precursor> precursors_;
public:

    typedef SparseVector<float>::const_iterator const_bin_iterator;
    typedef SparseVector<float>::iterator bin_iterator;

    /// default constructor
    BinnedSpectrum() = delete;

    /// detailed constructor
    BinnedSpectrum(float size, UInt spread, const PeakSpectrum& ps);

    /// copy constructor
    BinnedSpectrum(const BinnedSpectrum&) = default;

    /// destructor
    virtual ~BinnedSpectrum();

    /// assignment operator
    BinnedSpectrum& operator=(const BinnedSpectrum&) = default;

    /// equality operator
    bool operator==(const BinnedSpectrum& rhs) const
    {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
      return std::tie(bin_size_, bin_spread_, precursors_, bins_) 
          == std::tie(rhs.bin_size_, rhs.bin_spread_, rhs.precursors_, rhs.bins_);
#pragma clang diagnostic pop
    }

    /// inequality operator
    bool operator!=(const BinnedSpectrum& rhs) const;

    inline size_t getBinIndex(double mz) { return static_cast<int>(floor(mz / bin_size_)); }

    /// get the bin size
    inline double getBinSize() const { return bin_size_; }

    /// get the bin spread
    inline size_t getBinSpread() const { return bin_spread_; }

    /// get the number of bins
    inline size_t getBinNumber() const { return bins_.size(); }

    /// get the number of filled bins
    inline size_t getFilledBinNumber() const { return bins_.nonzero_size(); }

    /// immutable access to the bin container
    const SparseVector<float>& getBins() const;

    /// mutable access to the bin container
    SparseVector<float>& getBins();

    // inmutable access to precursors
    const std::vector<Precursor>& getPrecursors() const;

    /// mutable access to precursors
    std::vector<Precursor>& getPrecursors();

    /// returns the const begin iterator of the container
    inline const_bin_iterator begin() const { return bins_.begin(); }

    /// returns the const end iterator of the container
    inline const_bin_iterator end() const { return bins_.end(); }

    /// returns the begin iterator of the container
    inline bin_iterator begin() { return bins_.begin(); }

    /// returns the end iterator of the container
    inline bin_iterator end() { return bins_.end(); }

    /// Function to check comparability of two BinnedSpectrum objects,
    /// That is, if they have equal bin size and spread
    static bool isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b);
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

