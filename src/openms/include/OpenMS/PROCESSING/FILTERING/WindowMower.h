// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Timo Sachsenberg$
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <set>

namespace OpenMS
{

  /**
    @brief Retains the highest peaks in a sliding or jumping window

    @htmlinclude OpenMS_WindowMower.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI WindowMower :
    public DefaultParamHandler
  {
public:

    // @name Constructors, destructors and assignment operators
    // @{
    /// default constructor
    WindowMower();
    /// destructor
    ~WindowMower() override;

    /// copy constructor
    WindowMower(const WindowMower& source);
    /// assignment operator
    WindowMower& operator=(const WindowMower& source);
    // @}

    /// sliding window version (slower)
    template <typename SpectrumType>
    void filterPeakSpectrumForTopNInSlidingWindow(SpectrumType& spectrum)
    {
      typedef typename SpectrumType::ConstIterator ConstIterator;

      windowsize_ = (double)param_.getValue("windowsize");
      peakcount_ = (UInt)param_.getValue("peakcount");

      //copy spectrum
      SpectrumType old_spectrum = spectrum;
      old_spectrum.sortByPosition();

      //find high peak positions
      bool end  = false;
      std::set<double> positions;
      for (ConstIterator it = old_spectrum.begin(); it != old_spectrum.end(); ++it)
      {
        // copy the window from the spectrum
        SpectrumType window;
        for (ConstIterator it2 = it; (it2->getPosition() - it->getPosition() < windowsize_); )
        {
          window.push_back(*it2);
          if (++it2 == old_spectrum.end())
          {
            end = true;
            break;
          }
        }

        //extract peakcount most intense peaks
        window.sortByIntensity(true);
        for (Size i = 0; i < peakcount_; ++i)
        {
          if (i < window.size())
          {
            positions.insert(window[i].getMZ());
          }
        }
        //abort at the end of the spectrum
        if (end) break;
      }

      // select peaks that were retained
      std::vector<Size> indices;
      for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        if (positions.find(it->getMZ()) != positions.end())
        {
          Size index(it - spectrum.begin());
          indices.push_back(index);
        }
      }
      spectrum.select(indices);
    }

    void filterPeakSpectrum(PeakSpectrum& spectrum);

    void filterPeakMap(PeakMap& exp);

    // jumping window version (faster)
    template <typename SpectrumType>
    void filterPeakSpectrumForTopNInJumpingWindow(SpectrumType& spectrum)
    {
      if (spectrum.empty())
      {
        return;
      }

      spectrum.sortByPosition();

      windowsize_ = static_cast<double>(param_.getValue("windowsize"));
      peakcount_ = static_cast<UInt>(param_.getValue("peakcount"));

      // copy meta data
      SpectrumType out = spectrum;
      out.clear(false);

      SpectrumType peaks_in_window;
      double window_start = spectrum[0].getMZ();
      for (Size i = 0; i != spectrum.size(); ++i)
      {
        if (spectrum[i].getMZ() - window_start < windowsize_) // collect peaks in window
        {
          peaks_in_window.push_back(spectrum[i]);
        }
        else // step over window boundaries
        {
          window_start = spectrum[i].getMZ(); // as there might be large gaps between peaks resulting in empty windows, set new window start to next peak

          // copy N highest peaks to out
          if (peaks_in_window.size() > peakcount_)
          {
            std::partial_sort(peaks_in_window.begin(), peaks_in_window.begin() + peakcount_, peaks_in_window.end(), [](auto &left, auto &right) {typename SpectrumType::PeakType::IntensityLess cmp; return cmp(right, left);});
            copy(peaks_in_window.begin(), peaks_in_window.begin() + peakcount_, back_inserter(out));
          }
          else
          {
            std::sort(peaks_in_window.begin(), peaks_in_window.end(), [](auto &left, auto &right) {typename SpectrumType::PeakType::IntensityLess cmp; return cmp(right, left);});
            copy(peaks_in_window.begin(), peaks_in_window.end(), back_inserter(out));
          }

          peaks_in_window.clear(false);
          peaks_in_window.push_back(spectrum[i]);
        }
      }

      if (!peaks_in_window.empty()) // last window is not empty
      {
        // Note that the last window might be much smaller than windowsize.
        // Therefore the number of peaks copied from this window should be adapted accordingly.
        // Otherwise a lot of noise peaks are copied from each end of a spectrum.

        double last_window_size = peaks_in_window.back().getMZ() - window_start;
        double last_window_size_fraction = last_window_size / windowsize_;
        Size last_window_peakcount = static_cast<Size>(std::round(last_window_size_fraction * peakcount_));

        if (peaks_in_window.size() > last_window_peakcount)
        { // sort for last_window_peakcount highest peaks
          std::partial_sort(peaks_in_window.begin(), peaks_in_window.begin() + last_window_peakcount, peaks_in_window.end(), 
                            [](auto &left, auto &right) {typename SpectrumType::PeakType::IntensityLess cmp; return cmp(right, left);});
          std::copy(peaks_in_window.begin(), peaks_in_window.begin() + last_window_peakcount, back_inserter(out));
        }
        else
        {
          std::copy(peaks_in_window.begin(), peaks_in_window.end(), std::back_inserter(out));
        }
      }

      // select peaks that were retained
      std::vector<Size> indices;
      for (typename SpectrumType::ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        if (std::find(out.begin(), out.end(), *it) != out.end())
        {
          Size index(it - spectrum.begin());
          indices.push_back(index);
        }
      }
      spectrum.select(indices);

      return;
    }

    //TODO reimplement DefaultParamHandler::updateMembers_()

private:
    double windowsize_;
    UInt peakcount_;
  };

}


