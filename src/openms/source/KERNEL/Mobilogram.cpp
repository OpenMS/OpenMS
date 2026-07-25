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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <algorithm>
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
    // The data arrays are parallel to the peaks: keeping them while dropping the peaks
    // would leave the mobilogram internally inconsistent (and any later sort()/select()
    // would throw on the size mismatch).
    float_data_arrays_.clear();
    string_data_arrays_.clear();
    integer_data_arrays_.clear();
    RangeManager::clearRanges();
  }

  Mobilogram& Mobilogram::select(const std::vector<Size>& indices)
  {
    const Size snew = indices.size();
    const Size peaks_old = data_.size();

    std::vector<MobilityPeak1D> tmp;
    tmp.reserve(snew);
    for (Size i = 0; i < snew; ++i)
    {
      tmp.push_back(std::move(data_[indices[i]]));
    }
    data_.swap(tmp);

    std::vector<float> mda_tmp_float;
    for (Size i = 0; i < float_data_arrays_.size(); ++i)
    {
      if (float_data_arrays_[i].empty())
      {
        continue;
      }
      if (float_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "FloatDataArray[" + StringUtils::toStr(i) + "] size (" + StringUtils::toStr(float_data_arrays_[i].size()) +
                                        ") does not match mobilogram size (" + StringUtils::toStr(peaks_old) + ")");
      }
      mda_tmp_float.clear();
      mda_tmp_float.reserve(snew);
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp_float.push_back(std::move(float_data_arrays_[i][indices[j]]));
      }
      std::swap(float_data_arrays_[i], mda_tmp_float);
    }

    std::vector<std::string> mda_tmp_str;
    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      if (string_data_arrays_[i].empty())
      {
        continue;
      }
      if (string_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "StringDataArray[" + StringUtils::toStr(i) + "] size (" + StringUtils::toStr(string_data_arrays_[i].size()) +
                                        ") does not match mobilogram size (" + StringUtils::toStr(peaks_old) + ")");
      }
      mda_tmp_str.clear();
      mda_tmp_str.reserve(snew);
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp_str.push_back(std::move(string_data_arrays_[i][indices[j]]));
      }
      std::swap(string_data_arrays_[i], mda_tmp_str);
    }

    std::vector<Int> mda_tmp_int;
    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      if (integer_data_arrays_[i].empty())
      {
        continue;
      }
      if (integer_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "IntegerDataArray[" + StringUtils::toStr(i) + "] size (" + StringUtils::toStr(integer_data_arrays_[i].size()) +
                                        ") does not match mobilogram size (" + StringUtils::toStr(peaks_old) + ")");
      }
      mda_tmp_int.clear();
      mda_tmp_int.reserve(snew);
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp_int.push_back(std::move(integer_data_arrays_[i][indices[j]]));
      }
      std::swap(integer_data_arrays_[i], mda_tmp_int);
    }

    return *this;
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
