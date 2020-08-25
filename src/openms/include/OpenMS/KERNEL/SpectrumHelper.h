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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/LogStream.h>

#pragma once

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
  
} // namespace OpenMS


