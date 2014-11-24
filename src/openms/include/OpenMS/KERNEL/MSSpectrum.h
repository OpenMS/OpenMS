// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>

namespace OpenMS
{
  class Peak1D;

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
  template <typename PeakT = Peak1D>
  class MSSpectrum :
    private std::vector<PeakT>,
    public RangeManager<1>,
    public SpectrumSettings
  {
public:

    ///Float data array class
    class FloatDataArray :
      public MetaInfoDescription,
      public std::vector<float>
    {};

    ///Integer data array class
    class IntegerDataArray :
      public MetaInfoDescription,
      public std::vector<Int>
    {};

    ///String data array class
    class StringDataArray :
      public MetaInfoDescription,
      public std::vector<String>
    {};

    ///Comparator for the retention time.
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
    typedef std::vector<FloatDataArray> FloatDataArrays;
    /// String data array vector type
    typedef std::vector<StringDataArray> StringDataArrays;
    /// Integer data array vector type
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
      ms_level_ = source.ms_level_;
      name_ = source.name_;
      float_data_arrays_ = source.float_data_arrays_;
      string_data_arrays_ = source.string_data_arrays_;
      integer_data_arrays_ = source.integer_data_arrays_;

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
             SpectrumSettings::operator==(rhs)  &&
             retention_time_ == rhs.retention_time_ &&
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
    /// Returns the absolute retention time (is seconds)
    inline double getRT() const
    {
      return retention_time_;
    }

    /// Sets the absolute retention time (is seconds)
    inline void setRT(double rt)
    {
      retention_time_ = rt;
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
          std::sort(ContainerType::begin(), ContainerType::end(), reverseComparator(typename PeakType::IntensityLess()));
        }
        else
        {
          std::sort(ContainerType::begin(), ContainerType::end(), typename PeakType::IntensityLess());
        }
      }
      else
      {
        //sort index list
        std::vector<std::pair<typename PeakType::IntensityType, Size> > sorted_indices;
        sorted_indices.reserve(ContainerType::size());
        for (Size i = 0; i < ContainerType::size(); ++i)
        {
          sorted_indices.push_back(std::make_pair(ContainerType::operator[](i).getIntensity(), i));
        }

        if (reverse)
        {
          std::sort(sorted_indices.begin(), sorted_indices.end(), reverseComparator(PairComparatorFirstElement<std::pair<typename PeakType::IntensityType, Size> >()));
        }
        else
        {
          std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<typename PeakType::IntensityType, Size> >());
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

    /**
      @brief Lexicographically sorts the peaks by their position.

      The spectrum is sorted with respect to position. Meta data arrays will be sorted accordingly.
    */
    void sortByPosition()
    {
      if (float_data_arrays_.empty())
      {
        std::sort(ContainerType::begin(), ContainerType::end(), typename PeakType::PositionLess());
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
        std::sort(sorted_indices.begin(), sorted_indices.end(), PairComparatorFirstElement<std::pair<typename PeakType::PositionType, Size> >());

        //apply sorting to ContainerType and to metadataarrays
        ContainerType tmp;
        tmp.reserve(sorted_indices.size());
        for (Size i = 0; i < sorted_indices.size(); ++i)
        {
          tmp.push_back(*(ContainerType::begin() + (sorted_indices[i].second)));
        }
        ContainerType::swap(tmp);

        for (Size i = 0; i < float_data_arrays_.size(); ++i)
        {
          std::vector<float> mda_tmp;
          mda_tmp.reserve(float_data_arrays_[i].size());
          for (Size j = 0; j < float_data_arrays_[i].size(); ++j)
          {
            mda_tmp.push_back(*(float_data_arrays_[i].begin() + (sorted_indices[j].second)));
          }
          std::swap(float_data_arrays_[i], mda_tmp);
        }

        for (Size i = 0; i < string_data_arrays_.size(); ++i)
        {
          std::vector<String> mda_tmp;
          mda_tmp.reserve(string_data_arrays_[i].size());
          for (Size j = 0; j < string_data_arrays_[i].size(); ++j)
          {
            mda_tmp.push_back(*(string_data_arrays_[i].begin() + (sorted_indices[j].second)));
          }
          std::swap(string_data_arrays_[i], mda_tmp);
        }

        for (Size i = 0; i < integer_data_arrays_.size(); ++i)
        {
          std::vector<Int> mda_tmp;
          mda_tmp.reserve(integer_data_arrays_[i].size());
          for (Size j = 0; j < integer_data_arrays_[i].size(); ++j)
          {
            mda_tmp.push_back(*(integer_data_arrays_[i].begin() + (sorted_indices[j].second)));
          }
          std::swap(integer_data_arrays_[i], mda_tmp);
        }
      }
    }

    /// Checks if all peaks are sorted with respect to ascending m/z
    bool isSorted() const
    {
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
      if (ContainerType::size() == 0) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There must be at least one peak to determine the nearest peak!");

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
        ms_level_ = 1;
        name_.clear();
        float_data_arrays_.clear();
        string_data_arrays_.clear();
        integer_data_arrays_.clear();
      }
    }

protected:

    /// Retention time
    double retention_time_;

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
