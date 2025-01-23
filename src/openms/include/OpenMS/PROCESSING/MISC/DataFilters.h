// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

namespace OpenMS
{
  class Feature;
  class ConsensusFeature;
  /**
      @brief DataFilter array providing some convenience functions

      @note For features the meta data filtering works on the MetaDataInterface of the Feature.
      For peaks it works on the FloatDataArrays defined in MSSpectrum.
  */
  class OPENMS_DLLAPI DataFilters
  {
public:
    DataFilters() = default;

    ///Information to filter
    enum FilterType
    {
      INTENSITY,      ///< Filter the intensity value
      QUALITY,        ///< Filter the overall quality value
      CHARGE,         ///< Filter the charge value
      SIZE,           ///< Filter the number of subordinates/elements
      META_DATA       ///< Filter meta data
    };
    ///Filter operation
    enum FilterOperation
    {
      GREATER_EQUAL,  ///< Greater than the value or equal to the value
      EQUAL,          ///< Equal to the value
      LESS_EQUAL,     ///< Less than the value or equal to the value
      EXISTS          ///< Only for META_DATA filter type, tests if meta data exists
    };

    /// Representation of a peak/feature filter combining FilterType, FilterOperation and a value (either double or String)
    struct OPENMS_DLLAPI DataFilter
    {
      DataFilter(){};
      /// ctor for common case of numerical filter
      DataFilter(const FilterType type, const FilterOperation op, const double val, const String& meta_name = "")
        : field(type), op(op), value(val), value_string(), meta_name(meta_name), value_is_numerical(true)
      {};
      /// ctor for common case of string filter
      DataFilter(const FilterType type, const FilterOperation op, const String& val, const String& meta_name = "")
        : field(type), op(op), value(0.0), value_string(val), meta_name(meta_name), value_is_numerical(false)
      {};
      /// Field to filter
      FilterType field{ DataFilters::INTENSITY };
      /// Filter operation
      FilterOperation op{ DataFilters::GREATER_EQUAL} ;
      /// Value for comparison
      double value{ 0.0 };
      /// String value for comparison (for meta data)
      String value_string;
      /// Name of the considered meta information (key)
      String meta_name;
      /// use @p value or @p value_string ?
      bool value_is_numerical{ false };

      /// Returns a string representation of the filter
      String toString() const;

      /**
          @brief Parses @p filter and sets the filter properties accordingly

          This method accepts the format provided by toString().

          @exception Exception::InvalidValue is thrown when the filter is not formatted properly
      */
      void fromString(const String & filter);

      ///Equality operator
      bool operator==(const DataFilter & rhs) const;

      ///Inequality operator
      bool operator!=(const DataFilter & rhs) const;

    };

    /// Filter count
    Size size() const;

    /**
      @brief Filter accessor

      @exception Exception::IndexOverflow is thrown for invalid indices
    */
    const DataFilter & operator[](Size index) const;

    /// Adds a filter
    void add(const DataFilter & filter);

    /**
      @brief Removes the filter corresponding to @p index

      @exception Exception::IndexOverflow is thrown for invalid indices
    */
    void remove(Size index);

    /**
      @brief Replaces the filter corresponding to @p index

      @exception Exception::IndexOverflow is thrown for invalid indices
    */
    void replace(Size index, const DataFilter & filter);

    /// Removes all filters
    void clear();

    /// Enables/disables the all the filters
    void setActive(bool is_active);

    /**
        @brief Returns if the filters are enabled

        They are automatically enabled when a filter is added and
        automatically disabled when the last filter is removed
    */
    inline bool isActive() const
    {
      return is_active_;
    }

    /// Returns if the @p feature fulfills the current filter criteria
    bool passes(const Feature& feature) const;

    /// Returns if the @p consensus_feature fulfills the current filter criteria
    bool passes(const ConsensusFeature& consensus_feature) const;

    /// Returns if the a peak in a @p spectrum at @p peak_index fulfills the current filter criteria
    inline bool passes(const MSSpectrum& spectrum, Size peak_index) const
    {
      if (!is_active_) return true;

      for (Size i = 0; i < filters_.size(); i++)
      {
        const DataFilters::DataFilter & filter = filters_[i];
        if (filter.field == INTENSITY)
        {
          switch (filter.op)
          {
          case GREATER_EQUAL:
            if (spectrum[peak_index].getIntensity() < filter.value) return false;

            break;

          case EQUAL:
            if (spectrum[peak_index].getIntensity() != filter.value) return false;

            break;

          case LESS_EQUAL:
            if (spectrum[peak_index].getIntensity() > filter.value) return false;

            break;

          default:
            break;
          }
        }
        else if (filter.field == META_DATA)
        {
          const auto& f_arrays = spectrum.getFloatDataArrays();
          //find the right meta data array
          SignedSize f_index = -1;
          for (Size j = 0; j < f_arrays.size(); ++j)
          {
            if (f_arrays[j].getName() == filter.meta_name)
            {
              f_index = j;
              break;
            }
          }
          //if it is present, compare it
          if (f_index != -1)
          {
            if (filter.op == EQUAL && f_arrays[f_index][peak_index] != filter.value) return false;
            else if (filter.op == LESS_EQUAL && f_arrays[f_index][peak_index] > filter.value) return false;
            else if (filter.op == GREATER_EQUAL && f_arrays[f_index][peak_index] < filter.value) return false;
          }

          //if float array not found, search in integer arrays
          const typename MSSpectrum::IntegerDataArrays & i_arrays = spectrum.getIntegerDataArrays();
          //find the right meta data array
          SignedSize i_index = -1;
          for (Size j = 0; j < i_arrays.size(); ++j)
          {
            if (i_arrays[j].getName() == filter.meta_name)
            {
              i_index = j;
              break;
            }
          }
          //if it is present, compare it
          if (i_index != -1)
          {
            if (filter.op == EQUAL && i_arrays[i_index][peak_index] != filter.value) return false;
            else if (filter.op == LESS_EQUAL && i_arrays[i_index][peak_index] > filter.value) return false;
            else if (filter.op == GREATER_EQUAL && i_arrays[i_index][peak_index] < filter.value) return false;
          }

          //if it is not present, abort
          if (f_index == -1 && i_index == -1) return false;
        }
      }
      return true;
    }

    /// Returns if the a peak in a @p chrom at @p peak_index fulfills the current filter criteria
    inline bool passes(const MSChromatogram& chrom, Size peak_index) const
    {
      if (!is_active_) return true;

      for (Size i = 0; i < filters_.size(); i++)
      {
        const DataFilters::DataFilter& filter = filters_[i];
        if (filter.field == INTENSITY)
        {
          switch (filter.op)
          {
            case GREATER_EQUAL:
              if (chrom[peak_index].getIntensity() < filter.value)
                return false;

              break;

            case EQUAL:
              if (chrom[peak_index].getIntensity() != filter.value)
                return false;

              break;

            case LESS_EQUAL:
              if (chrom[peak_index].getIntensity() > filter.value)
                return false;

              break;

            default:
              break;
          }
        }
        else if (filter.field == META_DATA)
        {
          const auto& f_arrays = chrom.getFloatDataArrays();
          // find the right meta data array
          SignedSize f_index = -1;
          for (Size j = 0; j < f_arrays.size(); ++j)
          {
            if (f_arrays[j].getName() == filter.meta_name)
            {
              f_index = j;
              break;
            }
          }
          // if it is present, compare it
          if (f_index != -1)
          {
            if (filter.op == EQUAL && f_arrays[f_index][peak_index] != filter.value) return false;
            else if (filter.op == LESS_EQUAL && f_arrays[f_index][peak_index] > filter.value) return false;
            else if (filter.op == GREATER_EQUAL && f_arrays[f_index][peak_index] < filter.value) return false;
          }

          // if float array not found, search in integer arrays
          const typename MSSpectrum::IntegerDataArrays& i_arrays = chrom.getIntegerDataArrays();
          // find the right meta data array
          SignedSize i_index = -1;
          for (Size j = 0; j < i_arrays.size(); ++j)
          {
            if (i_arrays[j].getName() == filter.meta_name)
            {
              i_index = j;
              break;
            }
          }
          // if it is present, compare it
          if (i_index != -1)
          {
            if (filter.op == EQUAL && i_arrays[i_index][peak_index] != filter.value) return false;
            else if (filter.op == LESS_EQUAL && i_arrays[i_index][peak_index] > filter.value) return false;
            else if (filter.op == GREATER_EQUAL && i_arrays[i_index][peak_index] < filter.value) return false;
          }

          // if it is not present, abort
          if (f_index == -1 && i_index == -1) return false;
        }
      }
      return true;
    }

    /// Returns if the a peak in a @p mobilogram at @p peak_index fulfills the current filter criteria
    inline bool passes(const Mobilogram& mobilogram, Size peak_index) const
    {
      if (!is_active_) {
        return true;
      }
        

      for (Size i = 0; i < filters_.size(); i++)
      {
        const DataFilters::DataFilter& filter = filters_[i];
        if (filter.field == INTENSITY)
        {
          switch (filter.op)
          {
            case GREATER_EQUAL:
              if (mobilogram[peak_index].getIntensity() < filter.value)
                return false;

              break;

            case EQUAL:
              if (mobilogram[peak_index].getIntensity() != filter.value)
                return false;

              break;

            case LESS_EQUAL:
              if (mobilogram[peak_index].getIntensity() > filter.value)
                return false;

              break;

            default:
              break;
          }
        }
        else if (filter.field == META_DATA)
        { // no metadata arrays so far...
          return false;
        }
      }
      return true;
    }

  protected:
    ///Array of DataFilters
    std::vector<DataFilter> filters_;
    ///Vector of meta indices acting as index cache
    std::vector<Size> meta_indices_;

    ///Determines if the filters are activated
    bool is_active_ = false;

    ///Returns if the meta value at @p index of @p meta_interface (a peak or feature) passes the @p filter
    bool metaPasses_(const MetaInfoInterface& meta_interface, const DataFilters::DataFilter& filter, Size index) const;
  };

} //namespace

