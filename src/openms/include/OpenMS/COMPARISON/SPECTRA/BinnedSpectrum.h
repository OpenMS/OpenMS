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
// $Authors: Timo Sachsenberg, Mathias Walzer $
// --------------------------------------------------------------------------
//
#pragma once

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
    // smallest possible m/z value (needs to be >= 1)
    static constexpr const float MIN_MZ_ = 1.0;

public:
    /** Sensible default values and notes from doi:10.1007/s13361-015-1179-x
      * Low-resolution MS/MS data:
      *   bin width = 1.0005     
      *   offset = 0.4
      *   spread should be 0
      *
      * High-resolution MS/MS data:
      *   bin width = 0.02   
      *   offset = 0.0
      *   spread should be 0 
      *   Note: in sum scores, intensities from neighboring bins should be considered with half intensity of each flanking bin.
      *         @TODO: Weighted intensity spread is currently not implemented(but could replace the spread parameter).
      */ 
 
    // default bin width for low-resolution data (adapted from doi:10.1007/s13361-015-1179-x)
    static constexpr const float DEFAULT_BIN_WIDTH_LOWRES = 1.0005f;

    // default bin width for high-resolution data (adapted from doi:10.1007/s13361-015-1179-x)
    static constexpr const float DEFAULT_BIN_WIDTH_HIRES = 0.02f;

    /// default bin offset for high-resolution data (adapted from doi:10.1007/s13361-015-1179-x)
    static constexpr const float DEFAULT_BIN_OFFSET_HIRES = 0.0f;

    /// default bin offset for low-resolution data (adapted from doi:10.1007/s13361-015-1179-x)
    static constexpr const float DEFAULT_BIN_OFFSET_LOWRES = 0.4f;

    /// typedef for the underlying sparse vector
    using SparseVectorType = Eigen::SparseVector<float>;

    /// typedef for the index into the sparse vector
    using SparseVectorIndexType = Eigen::SparseVector<float>::Index;
 
    /// typedef for the index into the sparse vector
    using SparseVectorIteratorType = Eigen::SparseVector<float>::InnerIterator;

    /// the empty SparseVector
    static const SparseVectorType EmptySparseVector;

    /// default constructor
    // BinnedSpectrum() = delete;
    BinnedSpectrum() {}

    /// detailed constructor
    BinnedSpectrum(const PeakSpectrum& ps, float size, bool unit_ppm, UInt spread, float offset);

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

    /// return the bin index of a given m/z position
    inline SparseVectorIndexType getBinIndex(float mz) const 
    {
      if (unit_ppm_)
      {
        /*
         * By solving:    mz = MIN_MZ_ * (1.0 + bin_size_ * 1e-6)^index for index
         *     we get: index = floor(log(mz/MIN_MZ_)/log(1.0 + bin_size_ * 1e-6))
         * Note: for ppm we don't need to consider an offset_.
         */  
        return static_cast<SparseVectorIndexType>(floor(log(mz/MIN_MZ_)/log1p(bin_size_ * 1e-6)));
      }
      else 
      { // implemented as described in PMC4607604
        // Note: Consider a peak offset (important for low-resolution data, where most peak boundaries
        //       may fall on the mass peak apex. See publication for details.).
        return static_cast<SparseVectorIndexType>(floor(mz / bin_size_ + offset_)); 
      }
    }

    /// return the lower m/z of a bin given its index
    inline float getBinLowerMZ(size_t i) const
    {
      if (unit_ppm_)
      {
        // mz = MIN_MZ_ * (1.0 + bin_size_)^index for index
        return (MIN_MZ_ * pow(1.0 + bin_size_ * 1e-6, i));
      }
      else 
      { 
        return ((static_cast<float>(i) - offset_) * bin_size_);
      }
    }

    /// get the bin size
    inline float getBinSize() const { return bin_size_; }

    /// get the bin spread
    inline size_t getBinSpread() const { return bin_spread_; }

    /// immutable access to the bin container
    const SparseVectorType& getBins() const;

    /// mutable access to the bin container
    SparseVectorType& getBins();

    /// return offset
    inline float getOffset() const { return offset_; }

    /// immutable access to precursors
    const std::vector<Precursor>& getPrecursors() const;

    /// mutable access to precursors
    std::vector<Precursor>& getPrecursors();

    /// Check if two BinnedSpectrum objects have equally sized bins and offset.
    //  returns true if bin size, unit and offset are equal, otherwise false
    static bool isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b);

private:
    /// the spread to left or right
    UInt bin_spread_;

    /// the size of each bin
    float bin_size_;

    /// absolute bin size or relative bin size
    bool unit_ppm_;

    /// offset of bin start
    float offset_;

    /// bins
    SparseVectorType bins_;

    /// calculate binning of peak spectrum
    void binSpectrum_(const PeakSpectrum& ps);

    /// precursor information
    std::vector<Precursor> precursors_;
  };

}

