// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
  class String;
  /**
    @brief Helper functions for MSSpectrum and MSChromatogram.

    @ingroup Kernel
  */

  /// Returns an iterator to the first data array with the given name. 
  /// The end iterator is returned in case no data array with given name exists.
  template <class DataArrayT>
  typename DataArrayT::iterator getDataArrayByName(DataArrayT& a, const String& name)
  {
    typename DataArrayT::iterator it = a.begin();
    for (; it != a.end(); ++it)
    {
      if (it->getName() == name) return it;
    }
    return it;
  }
  
  template <class DataArrayT>
  typename DataArrayT::const_iterator getDataArrayByName(const DataArrayT& a, const String& name)
  {
    typename DataArrayT::const_iterator it = a.begin();
    for (; it != a.end(); ++it)
    {
      if (it->getName() == name) return it;
    }
    return it;
  }

  /// remove all peaks EXCEPT in the given range
  template <typename PeakContainerT>
  void removePeaks(
    PeakContainerT& p,
    const double pos_start,
    const double pos_end,
    const bool ignore_data_arrays = false
  )
  {
    typename PeakContainerT::iterator it_start = p.PosBegin(pos_start);
    typename PeakContainerT::iterator it_end = p.PosEnd(pos_end);
    if (!ignore_data_arrays)
    {
      Size hops_left = std::distance(p.begin(), it_start);
      Size n_elems = std::distance(it_start, it_end);

      typename PeakContainerT::StringDataArrays& SDAs = p.getStringDataArrays();
      for (DataArrays::StringDataArray& sda : SDAs)
      {
        if (sda.size() == p.size())
        {
          sda.erase(sda.begin() + hops_left + n_elems, sda.end());
          sda.erase(sda.begin(), sda.begin() + hops_left);
        }
      }

      typename PeakContainerT::FloatDataArrays& FDAs = p.getFloatDataArrays();
      for (DataArrays::FloatDataArray& fda : FDAs)
      {
        if (fda.size() == p.size())
        {
          fda.erase(fda.begin() + hops_left + n_elems, fda.end());
          fda.erase(fda.begin(), fda.begin() + hops_left);
        }
      }

      typename PeakContainerT::IntegerDataArrays& IDAs = p.getIntegerDataArrays();
      for (DataArrays::IntegerDataArray& ida : IDAs)
      {
        if (ida.size() == p.size())
        {
          ida.erase(ida.begin() + hops_left + n_elems, ida.end());
          ida.erase(ida.begin(), ida.begin() + hops_left);
        }
      }
    }
    if (it_start == it_end)
    { // no elements left
      p.resize(0);
    }
    else
    { // if it_end != it_start, the second erase operation is safe
      p.erase(it_end, p.end());
      p.erase(p.begin(), it_start);
    }
  }

  template <typename PeakContainerT>
  void subtractMinimumIntensity(PeakContainerT& p)
  {
    if (p.empty()) return;

    typename PeakContainerT::iterator it = std::min_element(p.begin(), p.end(),
      [](typename PeakContainerT::PeakType& a, typename PeakContainerT::PeakType& b)
      {
        return a.getIntensity() < b.getIntensity();
      });

    const double rebase = - it->getIntensity();
    for (typename PeakContainerT::PeakType& peak : p)
    {
      peak.setIntensity(peak.getIntensity() + rebase);
    }
    // Note: data arrays are not updated
  }
  
  /**
   * @brief Possible methods for merging peak intensities.
   * 
   * @see makePeakPositionUnique() 
   **/
  enum class IntensityAveragingMethod : int { MEDIAN, MEAN, SUM, MIN, MAX };
  
  /**
   * @brief Make peak positions unique.
   * 
   * A peak container may contain multiple peaks with the same position, i.e. either
   * an MSSpectrum containing peaks with the same m/z position, or an MSChromatogram
   * containing peaks with identical RT position. One scenario where this might happen
   * is when multiple spectra are merged to a single one.
   * 
   * The method combines peaks with the same position to a single one with the
   * intensity determined by method m.
   *
   * @param[in] p The peak container to be manipulated.
   * @param[in] m The method for determining peak intensity from peaks with same position (median, mean, sum, min, max).
   **/
  template <typename PeakContainerT>
  void makePeakPositionUnique(PeakContainerT& p, const IntensityAveragingMethod m = IntensityAveragingMethod::MEDIAN)
  {
    if (!p.getFloatDataArrays().empty() || !p.getStringDataArrays().empty() || !p.getIntegerDataArrays().empty())
    {
      OPENMS_LOG_WARN << "Warning: data arrays are being ignored in the method SpectrumHelper::makePeakPositionUnique().\n";
    }
    
    if (p.empty()) return;
    
    p.sortByPosition();
    
    double current_position = p.begin()->getPos();
    PeakContainerT p_new;
    double intensity_new(0);
    std::vector<double> intensities_at_same_position;
    for (typename PeakContainerT::PeakType& peak : p)
    {
      if (peak.getPos() > current_position)
      {
        // add a peak to the new peak container
        switch(m)
        {
          case IntensityAveragingMethod::MEDIAN: intensity_new = Math::median(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
          case IntensityAveragingMethod::MEAN: intensity_new = Math::mean(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
          case IntensityAveragingMethod::SUM: intensity_new = Math::sum(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
          case IntensityAveragingMethod::MIN: intensity_new = *std::min_element(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
          case IntensityAveragingMethod::MAX: intensity_new = *std::max_element(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
        }
        typename PeakContainerT::PeakType peak_new(current_position, intensity_new);
        p_new.push_back(peak_new);
        
        current_position = peak.getPos();
        intensities_at_same_position.clear();
      }
      
      intensities_at_same_position.push_back(peak.getIntensity());
    }
    
    // add the very last peak to the new peak container
    switch(m)
    {
      case IntensityAveragingMethod::MEDIAN : intensity_new = Math::median(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
      case IntensityAveragingMethod::MEAN : intensity_new = Math::mean(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
      case IntensityAveragingMethod::SUM : intensity_new = Math::sum(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
      case IntensityAveragingMethod::MIN : intensity_new = *std::min_element(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
      case IntensityAveragingMethod::MAX : intensity_new = *std::max_element(intensities_at_same_position.begin(), intensities_at_same_position.end()); break;
    }
    typename PeakContainerT::PeakType peak_new(current_position, intensity_new);
    p_new.push_back(peak_new);

    std::swap(p_new, p);
  }

  /**
   * @brief Copies only the meta data contained in the input spectrum to the output spectrum.
   * 
   * @note Actual data is not copied.
   *
   * @param[in] input The input spectrum.
   * @param[out] output The output spectrum (will be cleared and will contain all metadata of the input spectrum).
   * @param clear_spectrum Whether the output spectrum should be cleared first (all raw data and data arrays will be deleted)
   **/
  OPENMS_DLLAPI void copySpectrumMeta(const MSSpectrum & input, MSSpectrum & output, bool clear_spectrum = true);
  
} // namespace OpenMS

