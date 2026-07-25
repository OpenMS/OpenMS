// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/KERNEL/MSChromatogram.h>

#include <algorithm>

using namespace OpenMS;

std::ostream& OpenMS::operator<<(std::ostream& os, const MSChromatogram& chrom)
{
  os << "-- MSCHROMATOGRAM BEGIN --" << std::endl;

  //chromatogram settings
  os << static_cast<const ChromatogramSettings&>(chrom);

  //data list
  for (const ChromatogramPeak& pe : chrom)
  {
    os << pe << std::endl;
  }

  os << "-- MSCHROMATOGRAM END --" << std::endl;

  return os;
}

bool MSChromatogram::MZLess::operator()(const MSChromatogram &a, const MSChromatogram &b) const
{
  return a.getMZ() < b.getMZ();
}

MSChromatogram &MSChromatogram::operator=(const MSChromatogram &source)
{
  if (&source == this)
  {
    return *this;
  }

  RangeManagerType::operator=(source);
  ContainerType::operator=(source);
  ChromatogramSettings::operator=(source);

  name_ = source.name_;
  float_data_arrays_ = source.float_data_arrays_;
  string_data_arrays_ = source.string_data_arrays_;
  integer_data_arrays_ = source.integer_data_arrays_;

  return *this;
}

bool MSChromatogram::operator==(const MSChromatogram &rhs) const
{
  //name_ can differ => it is not checked; also ranges are not checked
  return std::operator==(*this, rhs) &&
         ChromatogramSettings::operator==(rhs) &&
         float_data_arrays_ == rhs.float_data_arrays_ &&
         string_data_arrays_ == rhs.string_data_arrays_ &&
         integer_data_arrays_ == rhs.integer_data_arrays_;
}

const std::string &MSChromatogram::getName() const
{
  return name_;
}

void MSChromatogram::setName(const std::string &name)
{
  name_ = name;
}

double MSChromatogram::getMZ() const
{
  return getProduct().getMZ();
}

const MSChromatogram::FloatDataArrays &MSChromatogram::getFloatDataArrays() const
{
  return float_data_arrays_;
}

MSChromatogram::FloatDataArrays &MSChromatogram::getFloatDataArrays()
{
  return float_data_arrays_;
}

const MSChromatogram::StringDataArrays &MSChromatogram::getStringDataArrays() const
{
  return string_data_arrays_;
}

MSChromatogram::StringDataArrays &MSChromatogram::getStringDataArrays()
{
  return string_data_arrays_;
}

const MSChromatogram::IntegerDataArrays &MSChromatogram::getIntegerDataArrays() const
{
  return integer_data_arrays_;
}

MSChromatogram::IntegerDataArrays &MSChromatogram::getIntegerDataArrays()
{
  return integer_data_arrays_;
}

void MSChromatogram::sortByIntensity(bool reverse)
{
  // fast path: nothing to keep in sync with the peaks
  if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
  {
    if (reverse)
    {
      std::stable_sort(ContainerType::begin(), ContainerType::end(), [](const auto& left, const auto& right) {
        PeakType::IntensityLess cmp;
        return cmp(right, left);
      });
    }
    else
    {
      std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::IntensityLess());
    }
    return;
  }

  // permute peaks and *all* data arrays together
  if (reverse)
  {
    sort([this](const Size i1, const Size i2) {
      return this->ContainerType::operator[](i2).getIntensity() < this->ContainerType::operator[](i1).getIntensity();
    });
  }
  else
  {
    sort([this](const Size i1, const Size i2) {
      return this->ContainerType::operator[](i1).getIntensity() < this->ContainerType::operator[](i2).getIntensity();
    });
  }
}

void MSChromatogram::sortByPosition()
{
  if (isSorted())
  {
    return;
  }

  // fast path: nothing to keep in sync with the peaks
  if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
  {
    std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::PositionLess());
    return;
  }

  // permute peaks and *all* data arrays together
  sort([this](const Size i1, const Size i2) {
    return this->ContainerType::operator[](i1).getPosition() < this->ContainerType::operator[](i2).getPosition();
  });
}

bool MSChromatogram::isSorted() const
{
  for (Size i = 1; i < this->size(); ++i)
  {
    if (this->operator[](i - 1).getRT() > this->operator[](i).getRT())
    {
      return false;
    }
  }
  return true;
}

Size MSChromatogram::findNearest(MSChromatogram::CoordinateType rt) const
{
  //no peak => no search
  if (empty())
  {
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one peak to determine the nearest peak!");
  }
  //search for position for inserting
  ConstIterator it = RTBegin(rt);
  //border cases
  if (it == ContainerType::begin())
  {
    return 0;
  }
  if (it == ContainerType::end())
  {
    return ContainerType::size() - 1;
  }
  //the peak before or the current peak are closest
  ConstIterator it2 = it;
  --it2;
  if (std::fabs(it->getRT() - rt) < std::fabs(it2->getRT() - rt))
  {
    return Size(it - ContainerType::begin());
  }
  else
  {
    return Size(it2 - ContainerType::begin());
  }
}

MSChromatogram::Iterator MSChromatogram::RTBegin(MSChromatogram::CoordinateType rt)
{
  PeakType p;
  p.setPosition(rt);
  return lower_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
}

MSChromatogram::Iterator
MSChromatogram::RTBegin(MSChromatogram::Iterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::Iterator end)
{
  PeakType p;
  p.setPosition(rt);
  return lower_bound(begin, end, p, PeakType::PositionLess());
}

MSChromatogram::Iterator MSChromatogram::RTEnd(MSChromatogram::CoordinateType rt)
{
  PeakType p;
  p.setPosition(rt);
  return upper_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
}

MSChromatogram::Iterator
MSChromatogram::RTEnd(MSChromatogram::Iterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::Iterator end)
{
  PeakType p;
  p.setPosition(rt);
  return upper_bound(begin, end, p, PeakType::PositionLess());
}

MSChromatogram::ConstIterator MSChromatogram::RTBegin(MSChromatogram::CoordinateType rt) const
{
  PeakType p;
  p.setPosition(rt);
  return lower_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
}

MSChromatogram::ConstIterator
MSChromatogram::RTBegin(MSChromatogram::ConstIterator begin, MSChromatogram::CoordinateType rt,
                        MSChromatogram::ConstIterator end) const
{
  PeakType p;
  p.setPosition(rt);
  return lower_bound(begin, end, p, PeakType::PositionLess());
}

MSChromatogram::ConstIterator MSChromatogram::RTEnd(MSChromatogram::CoordinateType rt) const
{
  PeakType p;
  p.setPosition(rt);
  return upper_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
}

MSChromatogram::ConstIterator
MSChromatogram::RTEnd(MSChromatogram::ConstIterator begin, MSChromatogram::CoordinateType rt,
                      MSChromatogram::ConstIterator end) const
{
  PeakType p;
  p.setPosition(rt);
  return upper_bound(begin, end, p, PeakType::PositionLess());
}

MSChromatogram::Iterator MSChromatogram::PosBegin(MSChromatogram::CoordinateType rt)
{
  return RTBegin(rt);
}

MSChromatogram::Iterator
MSChromatogram::PosBegin(MSChromatogram::Iterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::Iterator end)
{
  return RTBegin(begin, rt, end);
}


MSChromatogram::Iterator
MSChromatogram::PosEnd(MSChromatogram::Iterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::Iterator end)
{
  return RTEnd(begin, rt, end);
}

MSChromatogram::ConstIterator MSChromatogram::PosBegin(MSChromatogram::CoordinateType rt) const
{
  return RTBegin(rt);
}

MSChromatogram::ConstIterator
MSChromatogram::PosBegin(MSChromatogram::ConstIterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::ConstIterator end) const
{
  return RTBegin(begin, rt, end);
}

MSChromatogram::ConstIterator
MSChromatogram::PosEnd(MSChromatogram::ConstIterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::ConstIterator end) const
{
  return RTEnd(begin, rt, end);
}

MSChromatogram::Iterator MSChromatogram::PosEnd(MSChromatogram::CoordinateType rt)
{
  return RTEnd(rt);
}

MSChromatogram::ConstIterator MSChromatogram::PosEnd(MSChromatogram::CoordinateType rt) const 
{
  return RTEnd(rt);
}

void MSChromatogram::clear(bool clear_meta_data)
{
  ContainerType::clear();

  // The data arrays are parallel to the peaks: keeping them while dropping the peaks
  // would leave the chromatogram internally inconsistent (and any later sort()/select()
  // would throw on the size mismatch). They are therefore always cleared, matching
  // MSSpectrum::clear(). 'clear_meta_data' only governs the descriptive metadata.
  clearRanges();
  float_data_arrays_.clear();
  string_data_arrays_.clear();
  integer_data_arrays_.clear();

  if (clear_meta_data)
  {
    this->ChromatogramSettings::operator=(ChromatogramSettings()); // no "clear" method
    name_.clear();
  }
}

MSChromatogram& MSChromatogram::select(const std::vector<Size>& indices)
{
  Size snew = indices.size();
  ContainerType tmp;
  tmp.reserve(indices.size());

  const Size peaks_old = size();

  // Validate everything *before* touching any storage, so a rejected call leaves
  // the chromatogram exactly as it was instead of half-permuted.
  for (Size i = 0; i < snew; ++i)
  {
    if (indices[i] >= peaks_old)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Index " + StringUtils::toStr(indices[i]) + " is out of range for a chromatogram of size " +
                                      StringUtils::toStr(peaks_old));
    }
  }
  {
    auto check_sizes = [peaks_old](const auto& arrays, const char* what) {
      for (Size i = 0; i < arrays.size(); ++i)
      {
        if (!arrays[i].empty() && arrays[i].size() != peaks_old)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string(what) + "[" + StringUtils::toStr(i) + "] size (" + StringUtils::toStr(arrays[i].size()) +
                                          ") does not match chromatogram size (" + StringUtils::toStr(peaks_old) + ")");
        }
      }
    };
    check_sizes(float_data_arrays_, "FloatDataArray");
    check_sizes(string_data_arrays_, "StringDataArray");
    check_sizes(integer_data_arrays_, "IntegerDataArray");
  }

  for (Size i = 0; i < snew; ++i)
  {
    tmp.push_back(std::move(ContainerType::operator[](indices[i])));
  }
  ContainerType::swap(tmp);

  std::vector<float> mda_tmp_float;
  for (Size i = 0; i < float_data_arrays_.size(); ++i)
  {
    if (float_data_arrays_[i].empty())
    {
      continue;
    }
    if (float_data_arrays_[i].size() != peaks_old)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FloatDataArray[" + StringUtils::toStr(i) + "] size (" +
                                                                                StringUtils::toStr(float_data_arrays_[i].size()) + ") does not match chromatogram size (" + StringUtils::toStr(peaks_old) + ")");
    }

    mda_tmp_float.clear();
    mda_tmp_float.reserve(float_data_arrays_[i].size());
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
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StringDataArray[" + StringUtils::toStr(i) + "] size (" +
                                                                                StringUtils::toStr(string_data_arrays_[i].size()) + ") does not match chromatogram size (" + StringUtils::toStr(peaks_old) + ")");
    }

    mda_tmp_str.clear();
    mda_tmp_str.reserve(string_data_arrays_[i].size());
    for (Size j = 0; j < snew; ++j)
    {
      // copy, not move: 'indices' may name the same source entry twice, and a
      // moved-from std::string would yield an empty second copy
      mda_tmp_str.push_back(string_data_arrays_[i][indices[j]]);
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
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray[" + StringUtils::toStr(i) + "] size (" +
                                                                                StringUtils::toStr(integer_data_arrays_[i].size()) + ") does not match chromatogram size (" + StringUtils::toStr(peaks_old) + ")");
    }

    mda_tmp_int.clear();
    mda_tmp_int.reserve(integer_data_arrays_[i].size());
    for (Size j = 0; j < snew; ++j)
    {
      mda_tmp_int.push_back(std::move(integer_data_arrays_[i][indices[j]]));
    }
    std::swap(integer_data_arrays_[i], mda_tmp_int);
  }

  return *this;
}

// This helper function is based on the cstd::set_union implementation. It is different in that it has a separate concept of "close enough to merge"
// This is defined as having retention times of within 1/1000 seconds
// Note: We assume that RTs are distinct in each of the two Chromatograms but may be the same between Chromatograms.
OpenMS::MSChromatogram::Iterator setSumSimilarUnion(OpenMS::MSChromatogram::Iterator first1,
                    OpenMS::MSChromatogram::Iterator last1,
                    OpenMS::MSChromatogram::Iterator first2,
                    OpenMS::MSChromatogram::Iterator last2,
                    OpenMS::MSChromatogram::Iterator result)
{
  while (true)
  {
    if (first1 == last1)
    {
      return std::copy(first2,last2,result);
    }
    if (first2 == last2)
    {
      return std::copy(first1,last1,result);
    }
    auto smaller_RT = [](OpenMS::MSChromatogram::Iterator a, OpenMS::MSChromatogram::Iterator b)->bool
    {
      return round(a->getRT() * 1000.0) < round(b->getRT() * 1000.0);
    };

    if (smaller_RT(first1, first2))
    {
      *result = *first1; ++first1;
    }
    else if (smaller_RT(first2, first1))
    {
      *result = *first2; ++first2;
    }
    else
    { // approx. equal
      *result = *first1;
      result->setIntensity(result->getIntensity() + first2->getIntensity());
      ++first1;
      ++first2;
    }
    ++result;
  }
}

void MSChromatogram::updateRanges()
{
  #ifdef OPENMS_ASSERTIONS
    double rt_min = RangeRT::isEmpty() ? 0 : getMinRT();
    double rt_max = RangeRT::isEmpty() ? 0 : getMaxRT();
    double int_min = RangeIntensity::isEmpty() ? 0 : getMinIntensity();
    double int_max = RangeIntensity::isEmpty() ? 0 : getMaxIntensity();
  #endif

  clearRanges();
  for (const auto& peak : (ContainerType&) *this)
  {
    extendRT(peak.getRT());
    extendIntensity(peak.getIntensity());
  }

  #ifdef OPENMS_ASSERTIONS
    double rt_min_new = RangeRT::isEmpty() ? 0 : getMinRT();
    double rt_max_new = RangeRT::isEmpty() ? 0 : getMaxRT();
    double int_min_new = RangeIntensity::isEmpty() ? 0 : getMinIntensity();
    double int_max_new = RangeIntensity::isEmpty() ? 0 : getMaxIntensity();

    // check if all are equal and no update range was necessary
    if (rt_min_new == rt_min && rt_max_new == rt_max
      && int_min_new == int_min && int_max_new == int_max)
    {
      OPENMS_LOG_WARN << "Update ranges was called but ranges were already up-to-date" << std::endl;
    }
  #endif
}

void MSChromatogram::mergePeaks(MSChromatogram& other, bool add_meta)
{
  vector<ChromatogramPeak> temp;
  temp.resize(size() + other.size());
  auto new_end = setSumSimilarUnion(begin(), end(), other.begin(), other.end(), temp.begin());
  ContainerType::assign(temp.begin(), new_end);

  if (add_meta)
  {
    DoubleList ls;
    if (metaValueExists(Constants::UserParam::MERGED_CHROMATOGRAM_MZS))
    {
      ls = getMetaValue(Constants::UserParam::MERGED_CHROMATOGRAM_MZS).toDoubleList();
    }
    ls.push_back(other.getMZ());
    setMetaValue(Constants::UserParam::MERGED_CHROMATOGRAM_MZS, ls);
  }
}
