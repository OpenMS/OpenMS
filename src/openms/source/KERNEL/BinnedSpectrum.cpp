// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//

#include <OpenMS/KERNEL/BinnedSpectrum.h>

#include <Eigen/Sparse>

using namespace std;

/// typedef for the index into the sparse vector
using SparseVectorIndexType = Eigen::SparseVector<float>::Index;

/// typedef for the index into the sparse vector
using SparseVectorIteratorType = Eigen::SparseVector<float>::InnerIterator;

namespace OpenMS
{

  BinnedSpectrum::BinnedSpectrum(const PeakSpectrum& ps, float size, bool unit_ppm, UInt spread, float offset) :
    bin_spread_(spread), 
    bin_size_(size),
    unit_ppm_(unit_ppm),
    offset_(offset)
  {
    bins_ = new SparseVectorType(numeric_limits<SparseVectorIndexType>::max());
    precursors_ = ps.getPrecursors();
    binSpectrum_(ps);
  }

  BinnedSpectrum::BinnedSpectrum(const BinnedSpectrum& rhs) :
    bin_spread_(rhs.bin_spread_), 
    bin_size_(rhs.bin_size_),
    unit_ppm_(rhs.unit_ppm_),
    offset_(rhs.offset_),
    bins_(new SparseVectorType(*rhs.bins_)),
    precursors_(rhs.precursors_)
  {
  }

  BinnedSpectrum& BinnedSpectrum::operator=(const BinnedSpectrum& rhs)
  {
    if (this != &rhs)
    {
      bin_spread_ = rhs.bin_spread_;
      bin_size_ = rhs.bin_size_;
      unit_ppm_ = rhs.unit_ppm_;
      offset_ = rhs.offset_;
      precursors_ = rhs.precursors_;

      delete bins_;
      bins_ = new SparseVectorType(*rhs.bins_);
    }

    return *this;
  }

  BinnedSpectrum::~BinnedSpectrum()
  {
    delete bins_;
  }

  const BinnedSpectrum::SparseVectorType* BinnedSpectrum::getBins() const
  {
    return bins_;
  }

  BinnedSpectrum::SparseVectorType* BinnedSpectrum::getBins()
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

    if (ps.empty())
    { 
      return;
    }

    // delete bins_;
    // bins_ = new SparseVectorType(numeric_limits<SparseVectorIndexType>::max());

    for (auto const & p : ps)
    {
      // if bin size is in relative units (ppm), check if minimum value is >= 1 (otherwise we might get numerical problems with the negative log)
      OPENMS_PRECONDITION(!unit_ppm_ || p.getMZ() >= BinnedSpectrum::MIN_MZ_, "Spectrum with relative bin size contains peaks with m/z < 1");

      // e.g.: bin_size_ = 1.5: first bin covers range [0, 1.5) so peak at 1.5 falls in second bin (index 1)
      const size_t idx = getBinIndex(p.getMZ());

      // add peak to corresponding bin
      bins_->coeffRef(idx) += p.getIntensity();

      // add peak to neighboring bins
      for (Size j = 0; j < bin_spread_; ++j)
      {
         bins_->coeffRef(idx + j + 1) +=  p.getIntensity();
        
        // prevent spreading over left boundaries
        if (static_cast<int>(idx - j - 1) >= 0)
        {
          bins_->coeffRef(idx - j - 1) += p.getIntensity();
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
    if (bins_->nonZeros() != rhs.bins_->nonZeros())
    {
      return false;
    }

    // test non-sparse (non-zero) elements for equality
    SparseVectorIteratorType it(*bins_);
    SparseVectorIteratorType rhs_it(*rhs.bins_);  
    while (it)      
    {
      if (it.index() != rhs_it.index() || it.value() != rhs_it.value())
      {
        return false;
      }
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

  size_t BinnedSpectrum::getBinIndex(float mz) const 
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

  float BinnedSpectrum::getBinIntensity(double mz)
  {
    return bins_->coeffRef(getBinIndex(mz));
  }

}

