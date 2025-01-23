// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/KERNEL/Mobilogram.h>

#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <numeric>

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
    clearRanges();
    for (const auto& peak : data_)
    {
      extendMobility(peak.getMobility());
      extendIntensity(peak.getIntensity());
    }
  }

  String Mobilogram::getDriftTimeUnitAsString() const
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
    if (reverse)
    {
      std::stable_sort(begin(), end(), [](auto& left, auto& right) {
        PeakType::IntensityLess cmp;
        return cmp(right, left);
      });
    }
    else
    {
      std::stable_sort(begin(), end(), PeakType::IntensityLess());
    }
  }

  void Mobilogram::sortByPosition()
  {
    if (isSorted())
    {
      return;
    }
    std::stable_sort(begin(), end(), PeakType::PositionLess());
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
    // search for position for inserting
    ConstIterator it = MBBegin(mb);
    // border cases
    if (it == cbegin())
    {
      return 0;
    }
    if (it == cend())
    {
      return size() - 1;
    }
    // the peak before or the current peak are closest
    auto it2 = it;
    --it2;
    if (std::fabs(it->getMobility() - mb) < std::fabs(it2->getMobility() - mb))
    {
      return Size(it - cbegin());
    }
    else
    {
      return Size(it2 - cbegin());
    }
  }

  Int Mobilogram::findNearest(CoordinateType mb, CoordinateType tolerance) const
  {
    if (empty())
    {
      return -1;
    }
    Size i = findNearest(mb);
    const double found_mb = this->operator[](i).getMobility();
    if (found_mb >= mb - tolerance && found_mb <= mb + tolerance)
    {
      return static_cast<Int>(i);
    }
    else
    {
      return -1;
    }
  }

  Int Mobilogram::findNearest(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const
  {
    if (empty())
    {
      return -1;
    }
    // do a binary search for nearest peak first
    Size i = findNearest(mb);

    const double nearest_mb = data_[i].getMobility();

    if (nearest_mb < mb)
    {
      if (nearest_mb >= mb - tolerance_left)
      {
        return i; // success: nearest peak is in left tolerance window
      }
      else
      {
        if (i == this->size() - 1)
        {
          return -1; // we are at the last peak which is too far left
        }
        // Nearest peak is too far left so there can't be a closer peak in the left window.
        // There still might be a peak to the right of mz that falls in the right window
        ++i; // now we are at a peak exactly on or to the right of mz
        const double next_mz = data_[i].getMobility();
        if (next_mz <= mb + tolerance_right)
        {
          return i;
        }
      }
    }
    else
    {
      if (nearest_mb <= mb + tolerance_right)
      {
        return i; // success: nearest peak is in right tolerance window
      }
      else
      {
        if (i == 0)
        {
          return -1; // we are at the first peak which is too far right
        }
        --i; // now we are at a peak exactly on or to the right of mz
        const double next_mz = data_[i].getMobility();
        if (next_mz >= mb - tolerance_left)
        {
          return i;
        }
      }
    }

    // neither in the left nor the right tolerance window
    return -1;
  }

  Int Mobilogram::findHighestInWindow(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_right) const
  {
    if (empty())
    {
      return -1;
    }
    // get left/right iterator
    auto left = this->MBBegin(mb - tolerance_left);
    auto right = this->MBEnd(mb + tolerance_right);

    if (left == right)
    {
      return -1;
    }

    auto max_intensity_it = std::max_element(left, right, MobilityPeak1D::IntensityLess());

    // find peak (index) with highest intensity to expected position
    return max_intensity_it - this->begin();
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
    RangeManager::clearRanges();
  }

  Mobilogram::ConstIterator Mobilogram::getBasePeak() const
  {
    ConstIterator largest = cbegin();
    if (empty())
    {
      return largest;
    }
    ConstIterator current = cbegin();
    ++current;
    for (; current != cend(); ++current)
    {
      if (largest->getIntensity() < current->getIntensity())
      {
        largest = current;
      }
    }
    return largest;
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
