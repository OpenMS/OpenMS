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
#include <OpenMS/CONCEPT/Exception.h>

#include <Eigen/Sparse>

#include <cmath>

namespace OpenMS
{

  /**
    @brief This is a binned representation of a PeakSpectrum

    @param sz the size of the bins and
    @param sp number of neighboring bins to both sides affected by a peak contribution
    @param ps the PeakSpectrum, used to calculate the binned spectrum

    sz denotes the size of a bin in @p Th, thereby deciding the number of bins (all of size sz) the spectrum is discretized to.
    Each bin will represent a certain @p Th range and the peaks will be put in the respective bins and sum up inside.
    sp denotes the number of neighboring bins to the left and the number of neighboring bins to the right a peak is also added to.
    E.g. a BinnedSpectrum with binsize of 0.5 @p Th will have a peak at 100 @p Th in bin no. 200, a peak at 100.1 @p Th will be in bin no. 201.
    If the binspread is 1, the peak at 100 Th will be added to bin no. 199, 200 and 201.
    If the binspread is 2, the peak at 100 @p Th will also be added to bin no. 198 and 202, and so on.

    Many operations are provided by the underlying SparseVector implementation:
    - bin-wise addition (e.g.: c = a.getBins() + b.getBins())
    - bin-wise scaling  (e.g.: c = a.getBins() * 5f)
    - to get the number of filled bins, call: getBins().nonZeros()
    - many more...
    See the Eigen SparseVector implementation for details.

    Implementation detail: Eigen SparseVectors need to have the same dimensionality. EmptySparseVector provides an empty SparseVector
                           with compatible dimension to perform all supported operations.

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSpectrum
  {

public:
    /// typedef for the underlying sparse vector
    using SparseVectorType = Eigen::SparseVector<float>;

    /// typedef for the index into the sparse vector
    using SparseVectorIndexType = Eigen::SparseVector<float>::Index;
 
    /// typedef for the index into the sparse vector
    using SparseVectorIteratorType = Eigen::SparseVector<float>::InnerIterator;

    /// the empty SparseVector
    static const SparseVectorType EmptySparseVector;

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
    bool operator==(const BinnedSpectrum& rhs) const;

    /// inequality operator
    bool operator!=(const BinnedSpectrum& rhs) const;

    /// returns the bin intensity at a given m/z position 
    inline float getBinIntensity(double mz) { return bins_.coeffRef(getBinIndex(mz)); }

    /// return the bin index of a given m/z position. mainly for internal use
    inline SparseVectorIndexType getBinIndex(double mz) { return static_cast<SparseVectorIndexType>(floor(mz / bin_size_)); }

    /// get the bin size
    inline double getBinSize() const { return bin_size_; }

    /// get the bin spread
    inline size_t getBinSpread() const { return bin_spread_; }

    /// immutable access to the bin container
    const SparseVectorType& getBins() const;

    /// mutable access to the bin container
    SparseVectorType& getBins();

    // inmutable access to precursors
    const std::vector<Precursor>& getPrecursors() const;

    /// mutable access to precursors
    std::vector<Precursor>& getPrecursors();

    /// Function to check comparability of two BinnedSpectrum objects,
    /// That is, if they have equal bin size and spread
    static bool isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b);

private:
    /// the spread to left or right
    UInt bin_spread_;

    /// the size of each bin
    float bin_size_;

    /// bins
    SparseVectorType bins_;

    /// calculate binnning of peak spectrum
    void binSpectrum_(const PeakSpectrum& ps);

    /// precursor information
    std::vector<Precursor> precursors_;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUM_H

