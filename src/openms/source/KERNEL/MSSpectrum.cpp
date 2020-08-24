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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/FORMAT/PeakTypeEstimator.h>

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
      tmp.push_back(std::move(ContainerType::operator[](indices[i])));
    }
    ContainerType::swap(tmp);

    std::vector<float> mda_tmp_float;
    for (Size i = 0; i < float_data_arrays_.size(); ++i)
    {
      if (float_data_arrays_[i].empty()) continue;
      if (float_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FloatDataArray[" + String(i) + "] size (" +
                                                                                  String(float_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
      }

      mda_tmp_float.clear();
      mda_tmp_float.reserve(float_data_arrays_[i].size());
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp_float.push_back(std::move(float_data_arrays_[i][indices[j]]));
      }
      std::swap(float_data_arrays_[i], mda_tmp_float);
    }

    std::vector<String> mda_tmp_str;
    for (Size i = 0; i < string_data_arrays_.size(); ++i)
    {
      if (string_data_arrays_[i].empty()) continue;
      if (string_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "StringDataArray[" + String(i) + "] size (" +
                                                                                  String(string_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
      }

      mda_tmp_str.clear();
      mda_tmp_str.reserve(string_data_arrays_[i].size());
      for (Size j = 0; j < snew; ++j)
      {
        mda_tmp_str.push_back(std::move(string_data_arrays_[i][indices[j]]));
      }
      std::swap(string_data_arrays_[i], mda_tmp_str);
    }

    std::vector<Int> mda_tmp_int;
    for (Size i = 0; i < integer_data_arrays_.size(); ++i)
    {
      if (integer_data_arrays_[i].empty()) continue;
      if (integer_data_arrays_[i].size() != peaks_old)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IntegerDataArray[" + String(i) + "] size (" +
                                                                                  String(integer_data_arrays_[i].size()) + ") does not match spectrum size (" + String(peaks_old) + ")");
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

  SpectrumSettings::SpectrumType MSSpectrum::getType(const bool query_data) const
  {
    SpectrumSettings::SpectrumType t = SpectrumSettings::getType();
    // easy case: type is known
    if (t != SpectrumSettings::UNKNOWN) return t;

    // Some conversion software only annotate "MS:1000525 spectrum representation" leading to an UNKNOWN type
    // Fortunately, some store a data processing item that indicates that the data has been picked
    for (auto& dp : getDataProcessing())
    {
      if (dp->getProcessingActions().count(DataProcessing::PEAK_PICKING) == 1)
      {
        return SpectrumSettings::CENTROID;
      }
    }

    if (query_data)
    {
      return PeakTypeEstimator::estimateType(begin(), end());
    }
    return SpectrumSettings::UNKNOWN;
  }

  MSSpectrum::ConstIterator MSSpectrum::getBasePeak() const
  {
    ConstIterator largest = cbegin();
    if (empty()) return largest;
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

  MSSpectrum::Iterator MSSpectrum::getBasePeak()
  {
    ConstIterator largest = const_cast<const MSSpectrum&>(*this).getBasePeak();
    return begin() + std::distance(cbegin(), largest);
  }

  MSSpectrum::PeakType::IntensityType MSSpectrum::getTIC() const
  {
    return std::accumulate(cbegin(),
                           cend(),
                           0.0,
                           [](MSSpectrum::PeakType::IntensityType sum, const PeakType& p)
                              {
                                return sum + p.getIntensity();
                              });
  }

  void MSSpectrum::clear(bool clear_meta_data)
  {
    ContainerType::clear();

    if (clear_meta_data)
    {
      ContainerType::shrink_to_fit();

      clearRanges();
      this->SpectrumSettings::operator=(SpectrumSettings()); // no "clear" method
      retention_time_ = -1.0;
      drift_time_ = -1.0;
      drift_time_unit_ = MSSpectrum::DriftTimeUnit::NONE;
      ms_level_ = 1;
      name_.clear();
      name_.shrink_to_fit();
      float_data_arrays_.clear();
      float_data_arrays_.shrink_to_fit();
      string_data_arrays_.clear();
      string_data_arrays_.shrink_to_fit();
      integer_data_arrays_.clear();
      integer_data_arrays_.shrink_to_fit();
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

  Int MSSpectrum::findHighestInWindow(MSSpectrum::CoordinateType mz, MSSpectrum::CoordinateType tolerance_left,
                              MSSpectrum::CoordinateType tolerance_right) const
  {
    if (ContainerType::empty()) return -1;

    // get left/right iterator
    auto left = this->MZBegin(mz - tolerance_left);
    auto right = this->MZEnd(mz + tolerance_right);

    // no MS1 precursor peak in +- tolerance window found
    if  (left == right)
    {
      return -1;
    }

    auto max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

    // find peak (index) with highest intensity to expected position
    return (max_intensity_it - this->begin());
  }


  void MSSpectrum::sortByPositionPresorted(const std::vector<Chunk>& chunks)
  {
    if (chunks.size() == 1 && chunks[0].is_sorted) return;

    if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
    {
      std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::PositionLess());
    }
    else
    {
      std::vector<Size> select_indices(this->size());
      std::iota(select_indices.begin(), select_indices.end(), 0);

      auto comparePos = [this] (Size a, Size b) { return this->ContainerType::operator[](a).getPos() < this->ContainerType::operator[](b).getPos(); };

      // sort all chunks, that haven't been sorted yet
      for (Size i = 0; i < chunks.size(); ++i)
      {
        if (!chunks[i].is_sorted)
        {
          std::stable_sort(select_indices.begin() + chunks[i].start, select_indices.begin() + chunks[i].end, comparePos);
        }
      }

      // now we can recursively merge all chunks, which is faster than using stable_sort in the first place
      std::function<void(Size,Size)> rec;
      rec = [&chunks, &select_indices, &rec, &comparePos] (Size first, Size last)->void {
        if (last > first)
        {
          Size mid = first + (last - first) / 2;
          rec(first, mid);
          rec(mid + 1, last);
          std::inplace_merge(select_indices.begin() + chunks[first].start, select_indices.begin() + chunks[mid].end, select_indices.begin() + chunks[last].end, comparePos);
        }
      };

      rec(0, chunks.size() - 1);

      select(select_indices);
    }
  }

  void MSSpectrum::sortByPosition()
  {
    if (isSorted()) return;

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
    if (reverse && std::is_sorted(ContainerType::begin(), ContainerType::end(), reverseComparator(PeakType::IntensityLess()))) return;
    else if (!reverse && std::is_sorted(ContainerType::begin(), ContainerType::end(), PeakType::IntensityLess())) return;

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
    return std::is_sorted(ContainerType::begin(), ContainerType::end(), PeakType::PositionLess());
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
           drift_time_unit_ == rhs.drift_time_unit_ &&
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
    drift_time_unit_ = source.drift_time_unit_;
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
    drift_time_unit_(MSSpectrum::DriftTimeUnit::NONE),
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
    drift_time_unit_(source.drift_time_unit_),
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

  MSSpectrum::DriftTimeUnit MSSpectrum::getDriftTimeUnit() const
  {
    return drift_time_unit_;
  }

  void MSSpectrum::setDriftTimeUnit(DriftTimeUnit dt)
  {
    drift_time_unit_ = dt;
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

  bool MSSpectrum::containsIMData() const
  {
    const auto& s = *this;
    return (!s.getFloatDataArrays().empty() &&
      ( s.getFloatDataArrays()[0].getName().hasPrefix("Ion Mobility") ||
        s.getFloatDataArrays()[0].getName() == "ion mobility array" ||
        s.getFloatDataArrays()[0].getName() == "mean inverse reduced ion mobility array" ||
        s.getFloatDataArrays()[0].getName() == "ion mobility drift time")
      );
  }
}
