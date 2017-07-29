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

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/METADATA/DataArrays.h>

namespace OpenMS
{
  /**
    @brief The representation of a 1D spectrum.

    It contains peak data and metadata about specific instrument settings,
    acquisition settings, description of the meta values used in the peaks and precursor info
    (SpectrumSettings).

    Several MSSpectrum instances are contained in a peak map (MSExperiment), which is essentially
    a vector of spectra with additional information about the experiment.

    Precursor info from SpectrumSettings should only be used if this spectrum is a tandem-MS
    spectrum. The precursor spectrum is the first spectrum in MSExperiment, that has a lower
    MS-level than the current spectrum.

    @note For range operations, see \ref RangeUtils "RangeUtils module"!

    @ingroup Kernel
  */
  template <typename PeakT>
  class MSSpectrum :
    private std::vector<PeakT>,
    public RangeManager<1>,
    public SpectrumSettings
  {
public:

    /// Comparator for the retention time.
    struct RTLess :
      public std::binary_function<MSSpectrum, MSSpectrum, bool>
    {
      inline bool operator()(const MSSpectrum& a, const MSSpectrum& b) const
      {
        return a.getRT() < b.getRT();
      }

    };

    ///@name Base type definitions
    //@{
    /// Peak type
    typedef PeakT PeakType;
    /// Coordinate (m/z) type
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Spectrum base type
    typedef std::vector<PeakType> ContainerType;
    /// Float data array vector type
    typedef OpenMS::DataArrays::FloatDataArray FloatDataArray ;
    typedef std::vector<FloatDataArray> FloatDataArrays;
    /// String data array vector type
    typedef OpenMS::DataArrays::StringDataArray StringDataArray ;
    typedef std::vector<StringDataArray> StringDataArrays;
    /// Integer data array vector type
    typedef OpenMS::DataArrays::IntegerDataArray IntegerDataArray ;
    typedef std::vector<IntegerDataArray> IntegerDataArrays;
    //@}

    ///@name Peak container iterator type definitions
    //@{
    /// Mutable iterator
    typedef typename ContainerType::iterator Iterator;
    /// Non-mutable iterator
    typedef typename ContainerType::const_iterator ConstIterator;
    /// Mutable reverse iterator
    typedef typename ContainerType::reverse_iterator ReverseIterator;
    /// Non-mutable reverse iterator
    typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
    //@}

    ///@name Export methods from std::vector<PeakT>
    //@{
    using ContainerType::operator[];
    using ContainerType::begin;
    using ContainerType::rbegin;
    using ContainerType::end;
    using ContainerType::rend;
    using ContainerType::resize;
    using ContainerType::size;
    using ContainerType::push_back;
    using ContainerType::pop_back;
    using ContainerType::empty;
    using ContainerType::front;
    using ContainerType::back;
    using ContainerType::reserve;
    using ContainerType::insert;
    using ContainerType::erase;
    using ContainerType::swap;

    using typename ContainerType::iterator;
    using typename ContainerType::const_iterator;
    using typename ContainerType::size_type;
    using typename ContainerType::value_type;
    using typename ContainerType::reference;
    using typename ContainerType::const_reference;
    using typename ContainerType::pointer;
    using typename ContainerType::difference_type;
    //@}


    /// Constructor
    MSSpectrum() :
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

    /// Copy constructor
    MSSpectrum(const MSSpectrum& source) :
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

    /// Destructor
    ~MSSpectrum()
    {}

    /// Assignment operator
    MSSpectrum& operator=(const MSSpectrum& source)
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

    /// Assignment operator
    MSSpectrum& operator=(const SpectrumSettings & source)
    {
      SpectrumSettings::operator=(source);
      return *this;
    }

    /// Equality operator
    bool operator==(const MSSpectrum& rhs) const
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

    /// Equality operator
    bool operator!=(const MSSpectrum& rhs) const
    {
      return !(operator==(rhs));
    }

    // Docu in base class (RangeManager)
    virtual void updateRanges()
    {
      this->clearRanges();
      updateRanges_(ContainerType::begin(), ContainerType::end());
    }

    ///@name Accessors for meta information
    ///@{
    /// Returns the absolute retention time (in seconds)
    inline double getRT() const
    {
      return retention_time_;
    }

    /// Sets the absolute retention time (in seconds)
    inline void setRT(double rt)
    {
      retention_time_ = rt;
    }

    /**
      @brief Returns the ion mobility drift time in milliseconds (-1 means it is not set)

      @note Drift times may be stored directly as an attribute of the spectrum
      (if they relate to the spectrum as a whole). In case of ion mobility
      spectra, the drift time of the spectrum will always be set here while the
      drift times attribute in the Precursor class may often be unpopulated.
    */
    inline double getDriftTime() const
    {
      return drift_time_;
    }

    /**
      @brief Returns the ion mobility drift time in milliseconds
    */
    inline void setDriftTime(double dt)
    {
      drift_time_ = dt;
    }

    /**
      @brief Returns the MS level.

      For survey scans this is 1, for MS/MS scans 2, ...
    */
    inline UInt getMSLevel() const
    {
      return ms_level_;
    }

    /// Sets the MS level.
    inline void setMSLevel(UInt ms_level)
    {
      ms_level_ = ms_level;
    }

    /// Returns the name
    inline const String& getName() const
    {
      return name_;
    }

    /// Sets the name
    inline void setName(const String& name)
    {
      name_ = name;
    }

    //@}

    /**
      @name Peak data array methods

      These methods are used to annotate each peak in a spectrum with meta information.
      It is an intermediate way between storing the information in the peak's MetaInfoInterface
      and deriving a new peak type with members for this information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    //@{
    /// Returns a const reference to the float meta data arrays
    inline const FloatDataArrays& getFloatDataArrays() const
    {
      return float_data_arrays_;
    }

    /// Returns a mutable reference to the float meta data arrays
    inline FloatDataArrays& getFloatDataArrays()
    {
      return float_data_arrays_;
    }

    /// Sets the float meta data arrays
    inline void setFloatDataArrays(const FloatDataArrays& fda)
    {
      float_data_arrays_ = fda;
    }

    /// Returns a const reference to the string meta data arrays
    inline const StringDataArrays& getStringDataArrays() const
    {
      return string_data_arrays_;
    }

    /// Returns a mutable reference to the string meta data arrays
    inline StringDataArrays& getStringDataArrays()
    {
      return string_data_arrays_;
    }

    /// Sets the string meta data arrays
    inline void setStringDataArrays(const StringDataArrays& sda)
    {
      string_data_arrays_ = sda;
    }

    /// Returns a const reference to the integer meta data arrays
    inline const IntegerDataArrays& getIntegerDataArrays() const
    {
      return integer_data_arrays_;
    }

    /// Returns a mutable reference to the integer meta data arrays
    inline IntegerDataArrays& getIntegerDataArrays()
    {
      return integer_data_arrays_;
    }

    /// Sets the integer meta data arrays
    inline void setIntegerDataArrays(const IntegerDataArrays& ida)
    {
      integer_data_arrays_ = ida;
    }

    //@}

    ///@name Sorting peaks
    //@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
    */
    void sortByIntensity(bool reverse = false)
    {
      if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
      {
        if (reverse)
        {
          std::stable_sort(ContainerType::begin(), ContainerType::end(), reverseComparator(typename PeakType::IntensityLess()));
        }
        else
        {
          std::stable_sort(ContainerType::begin(), ContainerType::end(), typename PeakType::IntensityLess());
        }
      }
      else
      {
        // sort index list
        std::vector<std::pair<typename PeakType::IntensityType, Size> > sorted_indices;
        sorted_indices.reserve(ContainerType::size());
        for (Size i = 0; i < ContainerType::size(); ++i)
        {
          sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(), i));
        }

        if (reverse)
        {
          std::stable_sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement<std::pair<typename PeakType::IntensityType, Size> >()));
        }
        else
        {
          std::stable_sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<typename PeakType::IntensityType, Size> >());
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

    /**
      @brief Lexicographically sorts the peaks by their position.

      The spectrum is sorted with respect to position. Meta data arrays will be sorted accordingly.
    */
    void sortByPosition()
    {
      if (float_data_arrays_.empty() && string_data_arrays_.empty() && integer_data_arrays_.empty())
      {
        std::stable_sort(ContainerType::begin(), ContainerType::end(), typename PeakType::PositionLess());
      }
      else
      {
        //sort index list
        std::vector<std::pair<typename PeakType::PositionType, Size> > sorted_indices;
        sorted_indices.reserve(ContainerType::size());
        for (Size i = 0; i < ContainerType::size(); ++i)
        {
          sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getPosition(), i));
        }
        std::stable_sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<typename PeakType::PositionType, Size> >());

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

    /// Checks if all peaks are sorted with respect to ascending m/z
    bool isSorted() const
    {
      if (this->size() < 2) return true;

      for (Size i = 1; i < this->size(); ++i)
      {
        if (this->operator[](i - 1).getMZ() > this->operator[](i).getMZ()) return false;
      }
      return true;
    }

    //@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific m/z

      @param mz The searched for mass-to-charge ratio searched
      @return Returns the index of the peak.

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the spectrum is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType mz) const
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

    /**
      @brief Binary search for the peak nearest to a specific m/z given a +/- tolerance windows in Th

      @param mz The searched for mass-to-charge ratio searched
      @param tolerance The non-negative tolerance applied to both sides of mz

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
    */
    Int findNearest(CoordinateType mz, CoordinateType tolerance) const
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

    /**
      @brief Search for the peak nearest to a specific m/z given two +/- tolerance windows in Th

      @param mz The searched for mass-to-charge ratio searched
      @param tolerance_left The non-negative tolerance applied left of mz
      @param tolerance_right The non-negative tolerance applied right of mz

      @return Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
      @note Peaks exactly on borders are considered in tolerance window.
      @note Search for the left border is done using a binary search followed by a linear scan
    */
    Int findNearest(CoordinateType mz, CoordinateType tolerance_left, CoordinateType tolerance_right) const
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

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator MZBegin(CoordinateType mz)
    {
      PeakType p;
      p.setPosition(mz);
      return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    Iterator MZBegin(Iterator begin, CoordinateType mz, Iterator end)
    {
      PeakType p;
      p.setPosition(mz);
      return lower_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator MZEnd(CoordinateType mz)
    {
      PeakType p;
      p.setPosition(mz);
      return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    Iterator MZEnd(Iterator begin, CoordinateType mz, Iterator end)
    {
      PeakType p;
      p.setPosition(mz);
      return upper_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator MZBegin(CoordinateType mz) const
    {
      PeakType p;
      p.setPosition(mz);
      return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the spectrum is sorted with respect to m/z! Otherwise the result is undefined.
    */
    ConstIterator MZBegin(ConstIterator begin, CoordinateType mz, ConstIterator end) const
    {
      PeakType p;
      p.setPosition(mz);
      return lower_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator MZEnd(CoordinateType mz) const
    {
      PeakType p;
      p.setPosition(mz);
      return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the spectrum is sorted with respect to m/z. Otherwise the result is undefined.
    */
    ConstIterator MZEnd(ConstIterator begin, CoordinateType mz, ConstIterator end) const
    {
      PeakType p;
      p.setPosition(mz);
      return upper_bound(begin, end, p, typename PeakType::PositionLess());
    }

    //@}


    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data)
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

    /*
      @brief Select a (subset of) spectrum and its data_arrays, only retaining the indices given in @p indices

      @param indices Vector of indices to keep
      @return Reference to this MSSpectrum

    */
    MSSpectrum& select(const std::vector<Size>& indices)
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

protected:
   
    /// Retention time
    double retention_time_;

    /// Drift time
    double drift_time_;

    /// MS level
    UInt ms_level_;

    /// Name
    String name_;

    /// Float data arrays
    FloatDataArrays float_data_arrays_;

    /// String data arrays
    StringDataArrays string_data_arrays_;

    /// Integer data arrays
    IntegerDataArrays integer_data_arrays_;
  };

  /// Print the contents to a stream.
  template <typename PeakT>
  std::ostream& operator<<(std::ostream& os, const MSSpectrum<PeakT>& spec)
  {
    os << "-- MSSPECTRUM BEGIN --" << std::endl;

    //spectrum settings
    os << static_cast<const SpectrumSettings&>(spec);

    //peaklist
    for (typename MSSpectrum<PeakT>::ConstIterator it = spec.begin(); it != spec.end(); ++it)
    {
      os << *it << std::endl;
    }

    os << "-- MSSPECTRUM END --" << std::endl;
    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSSPECTRUM_H

