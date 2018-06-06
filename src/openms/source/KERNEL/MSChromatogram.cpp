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
  for (MSChromatogram::ConstIterator it = chrom.begin(); it != chrom.end(); ++it)
  {
    os << *it << std::endl;
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
  if (&source == this) { return *this;}

  ContainerType::operator=(source);
  RangeManager<1>::operator=(source);
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
         RangeManager<1>::operator==(rhs) &&
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
  if (float_data_arrays_.empty() && string_data_arrays_.size() && integer_data_arrays_.size())
  {
    if (reverse)
    {
      std::sort(ContainerType::begin(), ContainerType::end(), reverseComparator(PeakType::IntensityLess()));
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
      sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(), i));
    }

    if (reverse)
    {
      std::sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement<std::pair<PeakType::IntensityType, Size> >()));
    }
    else
    {
      std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<PeakType::IntensityType, Size> >());
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
      float_data_arrays_[i].swap(mda_tmp);
    }

    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      std::vector<String> mda_tmp;
      for (Size j = 0; j < string_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(string_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      string_data_arrays_[i].swap(mda_tmp);
    }

    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      std::vector<Int> mda_tmp;
      for (Size j = 0; j < integer_data_arrays_[i].size(); ++j)
      {
        mda_tmp.push_back(*(integer_data_arrays_[i].begin() + (sorted_indices[j].second)));
      }
      integer_data_arrays_[i].swap(mda_tmp);
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
      sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getPosition(), i));
    }
    std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<PeakType::PositionType, Size> >());

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
    if (this->operator[](i - 1).getRT() > this->operator[](i).getRT()) return false;
  }
  return true;
}

Size MSChromatogram::findNearest(MSChromatogram::CoordinateType rt) const
{
  //no peak => no search
  if (ContainerType::size() == 0) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one peak to determine the nearest peak!");

  //search for position for inserting
  ConstIterator it = RTBegin(rt);
  //border cases
  if (it == ContainerType::begin()) return 0;

  if (it == ContainerType::end()) return ContainerType::size() - 1;

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

MSChromatogram::Iterator MSChromatogram::PosEnd(MSChromatogram::CoordinateType rt)
{
  return RTEnd(rt);
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

MSChromatogram::ConstIterator MSChromatogram::PosEnd(MSChromatogram::CoordinateType rt) const
{
  return RTEnd(rt);
}

MSChromatogram::ConstIterator
MSChromatogram::PosEnd(MSChromatogram::ConstIterator begin, MSChromatogram::CoordinateType rt, MSChromatogram::ConstIterator end) const
{
  return RTEnd(begin, rt, end);
}

MSChromatogram::ConstIterator MSChromatogram::MZEnd(MSChromatogram::CoordinateType rt) const {return RTEnd(rt);}

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
