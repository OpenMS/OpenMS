// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h> // for OPENMS_ASSERTIONS
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/KERNEL/Mobilogram.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <algorithm>
#include <numeric>

#include "PeakContainerHelpers.h"

namespace OpenMS
{
  bool Mobilogram::RTLess::operator()(const Mobilogram& a, const Mobilogram& b) const
  {
    return a.getRT() < b.getRT();
  }

  bool Mobilogram::operator==(const Mobilogram& rhs) const
  {
    return data_ == rhs.data_ && retention_time_ == rhs.retention_time_ && drift_time_unit_ == rhs.drift_time_unit_;
  }
  void Mobilogram::updateRanges()
  {
    #ifdef OPENMS_ASSERTIONS
      double im_min = RangeMobility::isEmpty() ? 0 : getMinMobility();
      double im_max = RangeMobility::isEmpty() ? 0 : getMaxMobility();
      double int_min = RangeIntensity::isEmpty() ? 0 : getMinIntensity();
      double int_max = RangeIntensity::isEmpty() ? 0 : getMaxIntensity();
    #endif

    clearRanges();
    for (const auto& peak : data_)
    {
      extendMobility(peak.getMobility());
      extendIntensity(peak.getIntensity());
    }

    #ifdef OPENMS_ASSERTIONS
      // check if updateRanges() was necessary and find places where it was not
      double im_min_new = RangeMobility::isEmpty() ? 0 : getMinMobility();
      double im_max_new = RangeMobility::isEmpty() ? 0 : getMaxMobility();
      double int_min_new =  RangeIntensity::isEmpty() ? 0 : getMinIntensity();
      double int_max_new =  RangeIntensity::isEmpty() ? 0 : getMaxIntensity();

      if (im_min_new == im_min && im_max_new == im_max 
        && int_min_new == int_min && int_max_new == int_max)
      {
        OPENMS_LOG_WARN << "Update ranges was called but ranges were already up-to-date" << std::endl;
      } 
    #endif     
  }

  const Mobilogram::FloatDataArrays &Mobilogram::getFloatDataArrays() const
  {
    return float_data_arrays_;
  }

  Mobilogram::FloatDataArrays &Mobilogram::getFloatDataArrays()
  {
    return float_data_arrays_;
  }

  const Mobilogram::StringDataArrays &Mobilogram::getStringDataArrays() const
  {
    return string_data_arrays_;
  }

  Mobilogram::StringDataArrays &Mobilogram::getStringDataArrays()
  {
    return string_data_arrays_;
  }

  const Mobilogram::IntegerDataArrays &Mobilogram::getIntegerDataArrays() const
  {
    return integer_data_arrays_;
  }

  Mobilogram::IntegerDataArrays &Mobilogram::getIntegerDataArrays()
  {
    return integer_data_arrays_;
  }

  std::string Mobilogram::getDriftTimeUnitAsString() const
  {
    return NamesOfDriftTimeUnit[(size_t)drift_time_unit_];
  }

  void Mobilogram::setDriftTimeUnit(DriftTimeUnit dt) noexcept
  {
    drift_time_unit_ = dt;
  }
  void Mobilogram::sortByIntensity(bool reverse)
  {
    if (reverse && std::is_sorted(cbegin(), cend(), [](auto& left, auto& right) {
          PeakType::IntensityLess cmp;
          return cmp(right, left);
        }))
    {
      return;
    }
    else if (!reverse && std::is_sorted(cbegin(), cend(), PeakType::IntensityLess()))
    {
      return;
    }

    // fast path: nothing to keep in sync with the peaks
    if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
    {
      if (reverse)
      {
        std::stable_sort(begin(), end(), [](const auto& left, const auto& right) {
          PeakType::IntensityLess cmp;
          return cmp(right, left);
        });
      }
      else
      {
        std::stable_sort(begin(), end(), PeakType::IntensityLess());
      }
      return;
    }

    // permute peaks and *all* data arrays together
    if (reverse)
    {
      sort([this](const Size i1, const Size i2) { return data_[i2].getIntensity() < data_[i1].getIntensity(); });
    }
    else
    {
      sort([this](const Size i1, const Size i2) { return data_[i1].getIntensity() < data_[i2].getIntensity(); });
    }
  }

  void Mobilogram::sortByPosition()
  {
    if (isSorted())
    {
      return;
    }

    // fast path: nothing to keep in sync with the peaks
    if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
    {
      std::stable_sort(begin(), end(), PeakType::PositionLess());
      return;
    }

    // permute peaks and *all* data arrays together
    sort([this](const Size i1, const Size i2) { return data_[i1].getPosition() < data_[i2].getPosition(); });
  }

  bool Mobilogram::isSorted() const
  {
    return std::is_sorted(begin(), end(), PeakType::PositionLess());
  }

  Size Mobilogram::findNearest(CoordinateType mb) const
  {
    // no peak => no search
    if (empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one peak to determine the nearest peak!");
    }
    // search for position for inserting, then pick the closer of its two neighbours
    return Internal::PeakContainerHelpers::nearestIndex(cbegin(), cend(), MBBegin(mb), mb, [](const PeakType& p) { return p.getMobility(); });
  }

  Int Mobilogram::findNearest(CoordinateType mb, CoordinateType tolerance) const
  {
    return Internal::PeakContainerHelpers::nearestIndexInTolerance(*this, mb, tolerance, [](const PeakType& p) { return p.getMobility(); });
  }

  Int Mobilogram::findNearest(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const
  {
    return Internal::PeakContainerHelpers::nearestIndexInWindow(*this, mb, tolerance_left, tolerance_right,
                                                                [](const PeakType& p) { return p.getMobility(); });
  }

  Int Mobilogram::findHighestInWindow(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const
  {
    if (empty())
    {
      return -1;
    }
    // most intense peak in [mb - tolerance_left, mb + tolerance_right], -1 if there is none
    return Internal::PeakContainerHelpers::highestIndexInRange(cbegin(), MBBegin(mb - tolerance_left), MBEnd(mb + tolerance_right));
  }

  Mobilogram::Iterator Mobilogram::MBBegin(CoordinateType mb)
  {
    PeakType p;
    p.setPosition(mb);
    return lower_bound(begin(), end(), p, PeakType::PositionLess());
  }

  Mobilogram::Iterator Mobilogram::MBBegin(Iterator begin, CoordinateType mb, Iterator end)
  {
    PeakType p;
    p.setPosition(mb);
    return lower_bound(begin, end, p, PeakType::PositionLess());
  }

  Mobilogram::Iterator Mobilogram::MBEnd(CoordinateType mb)
  {
    PeakType p;
    p.setPosition(mb);
    return upper_bound(begin(), end(), p, PeakType::PositionLess());
  }

  Mobilogram::Iterator Mobilogram::MBEnd(Iterator begin, CoordinateType mb, Iterator end)
  {
    PeakType p;
    p.setPosition(mb);
    return upper_bound(begin, end, p, PeakType::PositionLess());
  }

  Mobilogram::ConstIterator Mobilogram::MBBegin(CoordinateType mb) const
  {
    PeakType p;
    p.setPosition(mb);
    return lower_bound(cbegin(), cend(), p, PeakType::PositionLess());
  }

  Mobilogram::ConstIterator Mobilogram::MBBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const
  {
    PeakType p;
    p.setPosition(mb);
    return lower_bound(begin, end, p, PeakType::PositionLess());
  }

  Mobilogram::ConstIterator Mobilogram::MBEnd(CoordinateType mb) const
  {
    PeakType p;
    p.setPosition(mb);
    return upper_bound(cbegin(), cend(), p, PeakType::PositionLess());
  }

  Mobilogram::ConstIterator Mobilogram::MBEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const
  {
    PeakType p;
    p.setPosition(mb);
    return upper_bound(begin, end, p, PeakType::PositionLess());
  }

  Mobilogram::Iterator Mobilogram::PosBegin(CoordinateType mz)
  {
    return MBBegin(mz);
  }

  Mobilogram::Iterator Mobilogram::PosBegin(Iterator begin, CoordinateType mb, Iterator end)
  {
    return MBBegin(begin, mb, end);
  }

  Mobilogram::ConstIterator Mobilogram::PosBegin(CoordinateType mb) const
  {
    return MBBegin(mb);
  }

  Mobilogram::ConstIterator Mobilogram::PosBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const
  {
    return MBBegin(begin, mb, end);
  }

  Mobilogram::Iterator Mobilogram::PosEnd(CoordinateType mb)
  {
    return MBEnd(mb);
  }

  Mobilogram::Iterator Mobilogram::PosEnd(Iterator begin, CoordinateType mb, Iterator end)
  {
    return MBEnd(begin, mb, end);
  }

  Mobilogram::ConstIterator Mobilogram::PosEnd(CoordinateType mb) const
  {
    return MBEnd(mb);
  }

  Mobilogram::ConstIterator Mobilogram::PosEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const
  {
    return MBEnd(begin, mb, end);
  }

  void Mobilogram::clear() noexcept
  {
    data_.clear();
    // The data arrays are parallel to the peaks: keeping them while dropping the peaks
    // would leave the mobilogram internally inconsistent (and any later sort()/select()
    // would throw on the size mismatch).
    float_data_arrays_.clear();
    string_data_arrays_.clear();
    integer_data_arrays_.clear();
    RangeManager::clearRanges();
  }

  void Mobilogram::checkDataArraySizes_() const
  {
    Internal::PeakContainerHelpers::checkDataArraySizes(data_.size(), float_data_arrays_, string_data_arrays_, integer_data_arrays_,
                                                        "mobilogram", OPENMS_PRETTY_FUNCTION);
  }

  Mobilogram& Mobilogram::select(const std::vector<Size>& indices)
  {
    // The range check is the guarantee of this entry point: an out-of-range index is rejected before
    // any storage is touched, so the mobilogram is left unchanged. selectUnchecked() then performs the
    // actual reordering.
    Internal::PeakContainerHelpers::checkSelectIndices(indices, data_.size(), "mobilogram", OPENMS_PRETTY_FUNCTION);
    return selectUnchecked(indices);
  }

  Mobilogram& Mobilogram::selectUnchecked(const std::vector<Size>& indices)
  {
    // Cheap and independent of the indices: guards a mis-sized data array from being indexed out of
    // bounds during the reorder below.
    checkDataArraySizes_();

    // Duplicate indices are undefined behaviour (asserted in debug builds only).
    Internal::PeakContainerHelpers::selectByIndices(data_, float_data_arrays_, string_data_arrays_, integer_data_arrays_, indices);
    return *this;
  }

  Mobilogram::ConstIterator Mobilogram::getBasePeak() const
  {
    return Internal::PeakContainerHelpers::basePeak(cbegin(), cend());
  }

  Mobilogram::Iterator Mobilogram::getBasePeak()
  {
    ConstIterator largest = const_cast<const Mobilogram&>(*this).getBasePeak();
    return begin() + std::distance(cbegin(), largest);
  }

  MobilityPeak1D::IntensityType Mobilogram::calculateTIC() const
  {
    return std::accumulate(cbegin(), cend(), 0.0f, [](PeakType::IntensityType sum, const PeakType& p) { return sum + p.getIntensity(); });
  }

  std::ostream& operator<<(std::ostream& os, const Mobilogram& mb)
  {
    os << "-- MOBILOGRAM BEGIN --\n";

    // peaklist
    for (const auto& peak : mb)
    {
      os << peak << '\n';
    }

    os << "-- MOBILOGRAM END --\n";
    return os;
  }
} // namespace OpenMS
