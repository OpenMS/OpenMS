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

#ifndef OPENMS_KERNEL_MSCHROMATOGRAM_H
#define OPENMS_KERNEL_MSCHROMATOGRAM_H

#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/METADATA/ChromatogramSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/METADATA/DataArrays.h>

namespace OpenMS
{
  class ChromatogramPeak;

  /**
    @brief The representation of a chromatogram.
    @ingroup Kernel
  */
  template <typename PeakT>
  class MSChromatogram :
    private std::vector<PeakT>,
    public RangeManager<1>,
    public ChromatogramSettings
  {

public:

    /// Comparator for the retention time.
    struct MZLess :
      public std::binary_function<MSChromatogram, MSChromatogram, bool>
    {
      inline bool operator()(const MSChromatogram& a, const MSChromatogram& b) const
      {
        return a.getMZ() < b.getMZ();
      }

    };

    ///@name Base type definitions
    ///@{
    /// Peak type
    typedef PeakT PeakType;
    /// Coordinate (RT) type
    typedef typename PeakType::CoordinateType CoordinateType;
    /// Chromatogram base type
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
    MSChromatogram() :
      ContainerType(),
      RangeManager<1>(),
      ChromatogramSettings(),
      name_(),
      float_data_arrays_(),
      string_data_arrays_(),
      integer_data_arrays_()
    {}

    /// Copy constructor
    MSChromatogram(const MSChromatogram& source) :
      ContainerType(source),
      RangeManager<1>(source),
      ChromatogramSettings(source),
      name_(source.name_),
      float_data_arrays_(source.float_data_arrays_),
      string_data_arrays_(source.string_data_arrays_),
      integer_data_arrays_(source.integer_data_arrays_)
    {}

    /// Destructor
    virtual ~MSChromatogram()
    {}

    /// Assignment operator
    MSChromatogram& operator=(const MSChromatogram& source)
    {
      if (&source == this) return *this;

      ContainerType::operator=(source);
      RangeManager<1>::operator=(source);
      ChromatogramSettings::operator=(source);

      name_ = source.name_;
      float_data_arrays_ = source.float_data_arrays_;
      string_data_arrays_ = source.string_data_arrays_;
      integer_data_arrays_ = source.integer_data_arrays_;

      return *this;
    }

    /// Equality operator
    bool operator==(const MSChromatogram& rhs) const
    {
      //name_ can differ => it is not checked
      return std::operator==(*this, rhs) &&
             RangeManager<1>::operator==(rhs) &&
             ChromatogramSettings::operator==(rhs)  &&
             float_data_arrays_ == rhs.float_data_arrays_ &&
             string_data_arrays_ == rhs.string_data_arrays_ &&
             integer_data_arrays_ == rhs.integer_data_arrays_;
    }

    /// Equality operator
    bool operator!=(const MSChromatogram& rhs) const
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

    ///@}

    /// returns the mz of the product entry, makes sense especially for MRM scans
    inline double getMZ() const
    {
      return getProduct().getMZ();
    }

    /**
      @name Peak data array methods

      These methods are used to annotate each peak in a chromatogram with meta information.
      It is an intermediate way between storing the information in the peak's MetaInfoInterface
      and deriving a new peak type with members for this information.

      These statements should help you chose which approach to use
        - Access to meta info arrays is slower than to a member variable
        - Access to meta info arrays is faster than to a %MetaInfoInterface
        - Meta info arrays are stored when using mzML format for storing
    */
    ///@{
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

    ///@}

    ///@name Sorting peaks
    ///@{
    /**
      @brief Lexicographically sorts the peaks by their intensity.

      Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly.
    */
    void sortByIntensity(bool reverse = false)
    {
      if (float_data_arrays_.empty() && string_data_arrays_.size() && integer_data_arrays_.size())
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

      The chromatogram is sorted with respect to position. Meta data arrays will be sorted
      accordingly.
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

    ///Checks if all peaks are sorted with respect to ascending RT
    bool isSorted() const
    {
      for (Size i = 1; i < this->size(); ++i)
      {
        if (this->operator[](i - 1).getRT() > this->operator[](i).getRT()) return false;
      }
      return true;
    }

    ///@}

    ///@name Searching a peak or peak range
    ///@{
    /**
      @brief Binary search for the peak nearest to a specific RT

      @param rt The searched for mass-to-charge ratio searched
      @return Returns the index of the peak.

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is undefined.

      @exception Exception::Precondition is thrown if the chromatogram is empty (not only in debug mode)
    */
    Size findNearest(CoordinateType rt) const
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

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to retention time! Otherwise the
      result is undefined.
    */
    Iterator RTBegin(CoordinateType rt)
    {
      PeakType p;
      p.setPosition(rt);
      return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    Iterator RTBegin(Iterator begin, CoordinateType rt, Iterator end)
    {
      PeakType p;
      p.setPosition(rt);
      return lower_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator RTEnd(CoordinateType rt)
    {
      PeakType p;
      p.setPosition(rt);
      return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    Iterator RTEnd(Iterator begin, CoordinateType rt, Iterator end)
    {
      PeakType p;
      p.setPosition(rt);
      return upper_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator RTBegin(CoordinateType rt) const
    {
      PeakType p;
      p.setPosition(rt);
      return lower_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range begin

      @note Make sure the chromatogram is sorted with respect to RT! Otherwise the result is
      undefined.
    */
    ConstIterator RTBegin(ConstIterator begin, CoordinateType rt, ConstIterator end) const
    {
      PeakType p;
      p.setPosition(rt);
      return lower_bound(begin, end, p, typename PeakType::PositionLess());
    }

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator RTEnd(CoordinateType rt) const
    {
      PeakType p;
      p.setPosition(rt);
      return upper_bound(ContainerType::begin(), ContainerType::end(), p, typename PeakType::PositionLess());
    }

    ConstIterator MZEnd(CoordinateType rt) const {return RTEnd(rt);}

    /**
      @brief Binary search for peak range end (returns the past-the-end iterator)

      @note Make sure the chromatogram is sorted with respect to RT. Otherwise the result is
      undefined.
    */
    ConstIterator RTEnd(ConstIterator begin, CoordinateType rt, ConstIterator end) const
    {
      PeakType p;
      p.setPosition(rt);
      return upper_bound(begin, end, p, typename PeakType::PositionLess());
    }

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
        this->ChromatogramSettings::operator=(ChromatogramSettings()); // no "clear" method
        name_.clear();
        float_data_arrays_.clear();
        string_data_arrays_.clear();
        integer_data_arrays_.clear();
      }
    }

    ///@}

protected:

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
  std::ostream& operator<<(std::ostream& os, const MSChromatogram<PeakT>& chrom)
  {
    os << "-- MSCHROMATOGRAM BEGIN --" << std::endl;

    //chromatogram settings
    os << static_cast<const ChromatogramSettings&>(chrom);

    //data list
    for (typename MSChromatogram<PeakT>::ConstIterator it = chrom.begin(); it != chrom.end(); ++it)
    {
      os << *it << std::endl;
    }

    os << "-- MSCHROMATOGRAM END --" << std::endl;

    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSCHROMATOGRAM_H
