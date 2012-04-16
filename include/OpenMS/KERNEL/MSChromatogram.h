// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSCHROMATOGRAM_H
#define OPENMS_KERNEL_MSCHROMATOGRAM_H

#include <OpenMS/METADATA/ChromatogramSettings.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/FORMAT/DB/PersistentObject.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>


namespace OpenMS
{
  /**
    @brief The representation of a chromatogram.
    @ingroup Kernel
  */
  template <typename PeakT = ChromatogramPeak>
  class MSChromatogram :
    public std::vector<PeakT>,
    public RangeManager<1>,
    public ChromatogramSettings,
    public PersistentObject
  {

public:

    /// Float data array class
    class OPENMS_DLLAPI FloatDataArray :
      public MetaInfoDescription,
      public std::vector<Real>
    {};

    /// String data array class
    class OPENMS_DLLAPI StringDataArray :
      public MetaInfoDescription,
      public std::vector<String>
    {};

    /// Float data array class
    class OPENMS_DLLAPI IntegerDataArray :
      public MetaInfoDescription,
      public std::vector<Int>
    {};

    /// Comparator for the retention time.
    struct MZLess :
      public std::binary_function<MSChromatogram, MSChromatogram, bool>
    {
      inline bool operator()(const MSChromatogram & a, const MSChromatogram & b) const
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


    /// Constructor
    MSChromatogram() :
      ContainerType(),
      RangeManager<1>(),
      ChromatogramSettings(),
      PersistentObject(),
      name_(),
      float_data_arrays_(),
      string_data_arrays_(),
      integer_data_arrays_()
    {}

    /// Copy constructor
    MSChromatogram(const MSChromatogram & source) :
      ContainerType(source),
      RangeManager<1>(source),
      ChromatogramSettings(source),
      PersistentObject(source),
      name_(source.name_),
      float_data_arrays_(source.float_data_arrays_),
      string_data_arrays_(source.string_data_arrays_),
      integer_data_arrays_(source.integer_data_arrays_)
    {}

    /// Destructor
    virtual ~MSChromatogram()
    {}

    /// Assignment operator
    MSChromatogram & operator=(const MSChromatogram & source)
    {
      if (&source == this) return *this;

      ContainerType::operator=(source);
      RangeManager<1>::operator=(source);
      ChromatogramSettings::operator=(source);
      PersistentObject::operator=(source);

      name_ = source.name_;
      float_data_arrays_ = source.float_data_arrays_;
      string_data_arrays_ = source.string_data_arrays_;
      integer_data_arrays_ = source.integer_data_arrays_;

      return *this;
    }

    /// Equality operator
    bool operator==(const MSChromatogram & rhs) const
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
    bool operator!=(const MSChromatogram & rhs) const
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
    inline const String & getName() const
    {
      return name_;
    }

    /// Sets the name
    inline void setName(const String & name)
    {
      name_ = name;
    }

    ///@}

    /// returns the mz of the product entry, makes sense especially for MRM scans
    inline DoubleReal getMZ() const
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
    inline const FloatDataArrays & getFloatDataArrays() const
    {
      return float_data_arrays_;
    }

    /// Returns a mutable reference to the float meta data arrays
    inline FloatDataArrays & getFloatDataArrays()
    {
      return float_data_arrays_;
    }

    /// Returns a const reference to the string meta data arrays
    inline const StringDataArrays & getStringDataArrays() const
    {
      return string_data_arrays_;
    }

    /// Returns a mutable reference to the string meta data arrays
    inline StringDataArrays & getStringDataArrays()
    {
      return string_data_arrays_;
    }

    /// Returns a const reference to the integer meta data arrays
    inline const IntegerDataArrays & getIntegerDataArrays() const
    {
      return integer_data_arrays_;
    }

    /// Returns a mutable reference to the integer meta data arrays
    inline IntegerDataArrays & getIntegerDataArrays()
    {
      return integer_data_arrays_;
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
          std::vector<Real> mda_tmp;
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
          std::vector<Real> mda_tmp;
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
      if (ContainerType::size() == 0) throw Exception::Precondition(__FILE__, __LINE__, __PRETTY_FUNCTION__, "There must be at least one peak to determine the nearest peak!");

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
        clearId();
        this->ChromatogramSettings::operator=(ChromatogramSettings());   // no "clear" method
        name_.clear();
        float_data_arrays_.clear();
        string_data_arrays_.clear();
        integer_data_arrays_.clear();
      }
    }

    ///@}

protected:
    // Docu in base class
    virtual void clearChildIds_()
    {}

    /// Name
    String name_;

    /// Float data arrays
    FloatDataArrays float_data_arrays_;

    /// String data arrays
    StringDataArrays string_data_arrays_;

    /// Intager data arrays
    IntegerDataArrays integer_data_arrays_;
  };

  /// Print the contents to a stream.
  template <typename PeakT>
  std::ostream & operator<<(std::ostream & os, const MSChromatogram<PeakT> & spec)
  {
    os << "-- MSSPECTRUM BEGIN --" << std::endl;

    //chromatogram settings
    os << static_cast<const ChromatogramSettings &>(spec);

    //peaklist
    os << static_cast<const typename MSChromatogram<PeakT>::ContainerType &>(spec);

    os << "-- MSSPECTRUM END --" << std::endl;

    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSCHROMATOGRAM_H
