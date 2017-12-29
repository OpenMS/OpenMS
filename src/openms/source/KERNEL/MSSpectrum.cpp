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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <ostream>

namespace OpenMS
{
  MSSpectrum &MSSpectrum::select(const std::vector<Size> &indices)
  {
    Size snew = indices.size();
    ContainerType tmp;
    tmp.reserve(indices.size());

    const Size peaks_old = size();

    for (Size i = 0; i < snew; ++i)
    {
      tmp.push_back(*(ContainerType::begin() + indices[i]));
    }
    ContainerType::swap(tmp);

    for (Size i = 0; i < float_data_arrays_.size(); ++i)
    {
      if (float_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FloatDataArray[" + String(i) + "] size (" +
                                                                                  String(float_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
      }

      std::vector<float> mda_tmp;
      mda_tmp.reserve(float_data_arrays_[i].size());
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp.push_back(*(float_data_arrays_[i].begin() + indices[j]));
      }
      std::swap(float_data_arrays_[i], mda_tmp);
    }

    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      if (string_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StringDataArray[" + String(i) + "] size (" +
                                                                                  String(string_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
      }
      std::vector<String> mda_tmp;
      mda_tmp.reserve(string_data_arrays_[i].size());
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp.push_back(*(string_data_arrays_[i].begin() + indices[j]));
      }
      std::swap(string_data_arrays_[i], mda_tmp);
    }

    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      if (integer_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray[" + String(i) + "] size (" +
                                                                                  String(integer_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
      }
      std::vector<Int> mda_tmp;
      mda_tmp.reserve(integer_data_arrays_[i].size());
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp.push_back(*(integer_data_arrays_[i].begin() + indices[j]));
      }
      std::swap(integer_data_arrays_[i], mda_tmp);
    }

    return *this;
  }

  void MSSpectrum::clear(bool clear_meta_data)
  {
    ContainerType::clear();

    if (clear_meta_data)
    {
      clearRanges();
      this->SpectrumSettings::operator=(SpectrumSettings()); // no "clear" method
      retention_time_ = -1.0;
      drift_time_ = -1.0;
      ms_level_ = 1;
      name_.clear();
      float_data_arrays_.clear();
      string_data_arrays_.clear();
      integer_data_arrays_.clear();
    }
  }

  MSSpectrum::ConstIterator
  MSSpectrum::MZEnd(MSSpectrum::ConstIterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::ConstIterator end) const
  {
    PeakType p;
    p.setPosition(mz);
    return upper_bound(begin, end, p, PeakType::PositionLess());
  }

  MSSpectrum::ConstIterator MSSpectrum::MZEnd(MSSpectrum::CoordinateType mz) const
  {
    PeakType p;
    p.setPosition(mz);
    return upper_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
  }

  MSSpectrum::ConstIterator MSSpectrum::MZBegin(MSSpectrum::ConstIterator begin, MSSpectrum::CoordinateType mz,
                                                MSSpectrum::ConstIterator end) const
  {
    PeakType p;
    p.setPosition(mz);
    return lower_bound(begin, end, p, PeakType::PositionLess());
  }

  Int MSSpectrum::findNearest(MSSpectrum::CoordinateType mz, MSSpectrum::CoordinateType tolerance_left,
                              MSSpectrum::CoordinateType tolerance_right) const
  {
    if (ContainerType::empty()) return -1;

    // do a binary search for nearest peak first
    Size i = findNearest(mz);

    const double nearest_mz = this->operator[](i).getMZ();

    if (nearest_mz < mz)
    {
      if (nearest_mz >= mz - tolerance_left)
      {
        return i; // success: nearest peak is in left tolerance window
      }
      else
      {
        if (i == this->size() - 1) return -1; // we are at the last peak which is too far left
        // Nearest peak is too far left so there can't be a closer peak in the left window.
        // There still might be a peak to the right of mz that falls in the right window
        ++i;  // now we are at a peak exactly on or to the right of mz
        const double next_mz = this->operator[](i).getMZ();
        if (next_mz <= mz + tolerance_right) return i;
      }
    }
    else
    {
      if (nearest_mz <= mz + tolerance_right)
      {
        return i; // success: nearest peak is in right tolerance window
      }
      else
      {
        if (i == 0) return -1; // we are at the first peak which is too far right
        --i;  // now we are at a peak exactly on or to the right of mz
        const double next_mz = this->operator[](i).getMZ();
        if (next_mz >= mz - tolerance_left) return i;
      }
    }

    // neither in the left nor the right tolerance window
    return -1;
  }

  Int MSSpectrum::findNearest(MSSpectrum::CoordinateType mz, MSSpectrum::CoordinateType tolerance) const
  {
    if (ContainerType::empty()) return -1;
    Size i = findNearest(mz);
    const double found_mz = this->operator[](i).getMZ();
    if (found_mz >= mz - tolerance && found_mz <= mz + tolerance)
    {
      return static_cast<Int>(i);
    }
    else
    {
      return -1;
    }
  }

  Size MSSpectrum::findNearest(MSSpectrum::CoordinateType mz) const
  {
    // no peak => no search
    if (ContainerType::size() == 0) throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "There must be at least one peak to determine the nearest peak!");

    // search for position for inserting
    ConstIterator it = MZBegin(mz);
    // border cases
    if (it == ContainerType::begin()) return 0;

    if (it == ContainerType::end()) return ContainerType::size() - 1;

    // the peak before or the current peak are closest
    ConstIterator it2 = it;
    --it2;
    if (std::fabs(it->getMZ() - mz) < std::fabs(it2->getMZ() - mz))
    {
      return Size(it - ContainerType::begin());
    }
    else
    {
      return Size(it2 - ContainerType::begin());
    }
  }

  void MSSpectrum::sortByPosition()
  {
    if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
    {
      std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::PositionLess());
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
      std::stable_sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<PeakType::PositionType, Size> >());

      // extract list of indices
      std::vector<Size> select_indices;
      select_indices.reserve(sorted_indices.size());
      for (Size i = 0; i < sorted_indices.size(); ++i)
      {
        select_indices.push_back(sorted_indices[i].second);
      }
      select(select_indices);
    }
  }

  void MSSpectrum::sortByIntensity(bool reverse)
  {
    if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
    {
      if (reverse)
      {
        std::stable_sort(ContainerType::begin(), ContainerType::end(), reverseComparator(PeakType::IntensityLess()));
      }
      else
      {
        std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::IntensityLess());
      }
    }
    else
    {
      // sort index list
      std::vector<std::pair<PeakType::IntensityType, Size> > sorted_indices;
      sorted_indices.reserve(ContainerType::size());
      for (Size i = 0; i < ContainerType::size(); ++i)
      {
        sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(), i));
      }

      if (reverse)
      {
        std::stable_sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement<std::pair<PeakType::IntensityType, Size> >()));
      }
      else
      {
        std::stable_sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<PeakType::IntensityType, Size> >());
      }

      // extract list of indices
      std::vector<Size> select_indices;
      select_indices.reserve(sorted_indices.size());
      for (Size i = 0; i < sorted_indices.size(); ++i)
      {
        select_indices.push_back(sorted_indices[i].second);
      }
      select(select_indices);
    }
  }

  bool MSSpectrum::isSorted() const
  {
    if (this->size() < 2) return true;

    for (Size i = 1; i < this->size(); ++i)
    {
      if (this->operator[](i - 1).getMZ() > this->operator[](i).getMZ()) return false;
    }
    return true;
  }

  bool MSSpectrum::operator==(const MSSpectrum &rhs) const
  {
    //name_ can differ => it is not checked
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
    return std::operator==(*this, rhs) &&
           RangeManager<1>::operator==(rhs) &&
           SpectrumSettings::operator==(rhs) &&
           retention_time_ == rhs.retention_time_ &&
           drift_time_ == rhs.drift_time_ &&
           ms_level_ == rhs.ms_level_ &&
           float_data_arrays_ == rhs.float_data_arrays_ &&
           string_data_arrays_ == rhs.string_data_arrays_ &&
           integer_data_arrays_ == rhs.integer_data_arrays_;

#pragma clang diagnostic pop
  }

  MSSpectrum &MSSpectrum::operator=(const MSSpectrum &source)
  {
    if (&source == this) return *this;

    ContainerType::operator=(source);
    RangeManager<1>::operator=(source);
    SpectrumSettings::operator=(source);

    retention_time_ = source.retention_time_;
    drift_time_ = source.drift_time_;
    ms_level_ = source.ms_level_;
    name_ = source.name_;
    float_data_arrays_ = source.float_data_arrays_;
    string_data_arrays_ = source.string_data_arrays_;
    integer_data_arrays_ = source.integer_data_arrays_;

    return *this;
  }

  MSSpectrum::MSSpectrum() :
    ContainerType(),
    RangeManager<1>(),
    SpectrumSettings(),
    retention_time_(-1),
    drift_time_(-1),
    ms_level_(1),
    name_(),
    float_data_arrays_(),
    string_data_arrays_(),
    integer_data_arrays_()
  {}

  MSSpectrum::MSSpectrum(const MSSpectrum &source) :
    ContainerType(source),
    RangeManager<1>(source),
    SpectrumSettings(source),
    retention_time_(source.retention_time_),
    drift_time_(source.drift_time_),
    ms_level_(source.ms_level_),
    name_(source.name_),
    float_data_arrays_(source.float_data_arrays_),
    string_data_arrays_(source.string_data_arrays_),
    integer_data_arrays_(source.integer_data_arrays_)
  {}

  MSSpectrum &MSSpectrum::operator=(const SpectrumSettings &source)
  {
    SpectrumSettings::operator=(source);
    return *this;
  }

  void MSSpectrum::updateRanges()
  {
    this->clearRanges();
    updateRanges_(ContainerType::begin(), ContainerType::end());
  }

  double MSSpectrum::getRT() const
  {
    return retention_time_;
  }

  void MSSpectrum::setRT(double rt)
  {
    retention_time_ = rt;
  }

  double MSSpectrum::getDriftTime() const
  {
    return drift_time_;
  }

  void MSSpectrum::setDriftTime(double dt)
  {
    drift_time_ = dt;
  }

  UInt MSSpectrum::getMSLevel() const
  {
    return ms_level_;
  }

  void MSSpectrum::setMSLevel(UInt ms_level)
  {
    ms_level_ = ms_level;
  }

  const String &MSSpectrum::getName() const
  {
    return name_;
  }

  void MSSpectrum::setName(const String &name)
  {
    name_ = name;
  }

  const MSSpectrum::FloatDataArrays &MSSpectrum::getFloatDataArrays() const
  {
    return float_data_arrays_;
  }

  void MSSpectrum::setFloatDataArrays(const MSSpectrum::FloatDataArrays &fda)
  {
    float_data_arrays_ = fda;
  }

  const MSSpectrum::StringDataArrays &MSSpectrum::getStringDataArrays() const
  {
    return string_data_arrays_;
  }

  void MSSpectrum::setStringDataArrays(const MSSpectrum::StringDataArrays &sda)
  {
    string_data_arrays_ = sda;
  }

  MSSpectrum::StringDataArrays &MSSpectrum::getStringDataArrays()
  {
    return string_data_arrays_;
  }

  const MSSpectrum::IntegerDataArrays &MSSpectrum::getIntegerDataArrays() const
  {
    return integer_data_arrays_;
  }

  MSSpectrum::IntegerDataArrays &MSSpectrum::getIntegerDataArrays()
  {
    return integer_data_arrays_;
  }

  void MSSpectrum::setIntegerDataArrays(const MSSpectrum::IntegerDataArrays &ida)
  {
    integer_data_arrays_ = ida;
  }

  MSSpectrum::Iterator MSSpectrum::MZBegin(MSSpectrum::CoordinateType mz)
  {
    PeakType p;
    p.setPosition(mz);
    return lower_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
  }

  MSSpectrum::Iterator
  MSSpectrum::MZBegin(MSSpectrum::Iterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::Iterator end)
  {
    PeakType p;
    p.setPosition(mz);
    return lower_bound(begin, end, p, PeakType::PositionLess());
  }

  MSSpectrum::Iterator MSSpectrum::MZEnd(MSSpectrum::CoordinateType mz)
  {
    PeakType p;
    p.setPosition(mz);
    return upper_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
  }

  MSSpectrum::Iterator
  MSSpectrum::MZEnd(MSSpectrum::Iterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::Iterator end)
  {
    PeakType p;
    p.setPosition(mz);
    return upper_bound(begin, end, p, PeakType::PositionLess());
  }

  MSSpectrum::ConstIterator MSSpectrum::MZBegin(MSSpectrum::CoordinateType mz) const
  {
    PeakType p;
    p.setPosition(mz);
    return lower_bound(ContainerType::begin(), ContainerType::end(), p, PeakType::PositionLess());
  }

  MSSpectrum::Iterator MSSpectrum::PosBegin(MSSpectrum::CoordinateType mz)
  {
    return MZBegin(mz);
  }

  MSSpectrum::Iterator
  MSSpectrum::PosBegin(MSSpectrum::Iterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::Iterator end)
  {
    return MZBegin(begin, mz, end);
  }

  MSSpectrum::Iterator MSSpectrum::PosEnd(MSSpectrum::CoordinateType mz)
  {
    return MZEnd(mz);
  }

  MSSpectrum::Iterator
  MSSpectrum::PosEnd(MSSpectrum::Iterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::Iterator end)
  {
    return MZEnd(begin, mz, end);
  }

  MSSpectrum::ConstIterator MSSpectrum::PosBegin(MSSpectrum::CoordinateType mz) const
  {
    return MZBegin(mz);
  }

  MSSpectrum::ConstIterator
  MSSpectrum::PosBegin(MSSpectrum::ConstIterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::ConstIterator end) const
  {
    return MZBegin(begin, mz, end);
  }

  MSSpectrum::ConstIterator MSSpectrum::PosEnd(MSSpectrum::CoordinateType mz) const
  {
    return MZEnd(mz);
  }

  MSSpectrum::ConstIterator
  MSSpectrum::PosEnd(MSSpectrum::ConstIterator begin, MSSpectrum::CoordinateType mz, MSSpectrum::ConstIterator end) const
  {
    return MZEnd(begin, mz, end);
  }

  bool MSSpectrum::RTLess::operator()(const MSSpectrum &a, const MSSpectrum &b) const {
    return a.getRT() < b.getRT();
  }
}
