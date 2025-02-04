// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSChromatogram.h>

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

  ContainerType::operator=(source);
  RangeManagerType::operator=(source);
  ChromatogramSettings::operator=(source);

  name_ = source.name_;
  float_data_arrays_ = source.float_data_arrays_;
  string_data_arrays_ = source.string_data_arrays_;
  integer_data_arrays_ = source.integer_data_arrays_;

  return *this;
}

bool MSChromatogram::operator==(const MSChromatogram &rhs) const
{
  //name_ can differ => it is not checked
  return std::operator==(*this, rhs) &&
         RangeManagerType::operator==(rhs) &&
         ChromatogramSettings::operator==(rhs)  &&
         float_data_arrays_ == rhs.float_data_arrays_ &&
         string_data_arrays_ == rhs.string_data_arrays_ &&
         integer_data_arrays_ == rhs.integer_data_arrays_;
}

const String &MSChromatogram::getName() const
{
  return name_;
}

void MSChromatogram::setName(const String &name)
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

void MSChromatogram::sortByIntensity(bool reverse) {
  if (float_data_arrays_.empty() && !string_data_arrays_.empty() && !integer_data_arrays_.empty())
  {
    if (reverse)
    {
      std::sort(ContainerType::begin(), ContainerType::end(), [](auto &left, auto &right) {PeakType::IntensityLess cmp; return cmp(right, left);});
    }
    else
    {
      std::sort(ContainerType::begin(), ContainerType::end(), PeakType::IntensityLess());
    }
  }
  else
  {
    //sort index list
    std::vector<std::pair<PeakType::IntensityType, Size> > sorted_indices;
    sorted_indices.reserve(ContainerType::size());
    for (Size i = 0; i < ContainerType::size(); ++i)
    {
      sorted_indices.emplace_back(ContainerType::operator[](i).getIntensity(), i);
    }

    if (reverse)
    {
      std::sort(sorted_indices.begin(), sorted_indices.end(),  [](auto& left, auto& right){return left > right;});
    }
    else
    {
      std::sort(sorted_indices.begin(), sorted_indices.end());
    }

    //apply sorting to ContainerType and to meta data arrays
    ContainerType tmp;
    for (Size i = 0; i < sorted_indices.size(); ++i)
    {
      tmp.push_back(*(ContainerType::begin() + (sorted_indices[i].second)));
    }
    ContainerType::swap(tmp);

    for (Size i = 0; i < float_data_arrays_.size(); ++i)
    {
      std::vector<float> mda_tmp;
      for (Size j = 0; j < float_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(float_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      mda_tmp.swap(float_data_arrays_[i]);
    }

    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      std::vector<String> mda_tmp;
      for (Size j = 0; j < string_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(string_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      mda_tmp.swap(string_data_arrays_[i]);
    }

    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      std::vector<Int> mda_tmp;
      for (Size j = 0; j < integer_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(integer_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      mda_tmp.swap(integer_data_arrays_[i]);
    }
  }
}

void MSChromatogram::sortByPosition()
{
  if (float_data_arrays_.empty())
  {
    std::sort(ContainerType::begin(), ContainerType::end(), PeakType::PositionLess());
  }
  else
  {
    //sort index list
    std::vector<std::pair<PeakType::PositionType, Size> > sorted_indices;
    sorted_indices.reserve(ContainerType::size());
    for (Size i = 0; i < ContainerType::size(); ++i)
    {
      sorted_indices.emplace_back(ContainerType::operator[](i).getPosition(), i);
    }
    std::sort(sorted_indices.begin(), sorted_indices.end());

    //apply sorting to ContainerType and to metadataarrays
    ContainerType tmp;
    for (Size i = 0; i < sorted_indices.size(); ++i)
    {
      tmp.push_back(*(ContainerType::begin() + (sorted_indices[i].second)));
    }
    ContainerType::swap(tmp);

    for (Size i = 0; i < float_data_arrays_.size(); ++i)
    {
      std::vector<float> mda_tmp;
      for (Size j = 0; j < float_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(float_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      std::swap(float_data_arrays_[i], mda_tmp);
    }

    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      std::vector<String> mda_tmp;
      for (Size j = 0; j < string_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(string_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      std::swap(string_data_arrays_[i], mda_tmp);
    }

    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      std::vector<Int> mda_tmp;
      for (Size j = 0; j < integer_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(integer_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      std::swap(integer_data_arrays_[i], mda_tmp);
    }
  }
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

  if (clear_meta_data)
  {
    clearRanges();
    this->ChromatogramSettings::operator=(ChromatogramSettings()); // no "clear" method
    name_.clear();
    float_data_arrays_.clear();
    string_data_arrays_.clear();
    integer_data_arrays_.clear();
  }
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
