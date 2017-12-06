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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Timo Sachsenberg$
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <set>

namespace OpenMS
{

  /**
    @brief WindowMower augments the highest peaks in a sliding or jumping window

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
            std::partial_sort(peaks_in_window.begin(), peaks_in_window.begin() + peakcount_, peaks_in_window.end(), reverseComparator(typename SpectrumType::PeakType::IntensityLess()));
            copy(peaks_in_window.begin(), peaks_in_window.begin() + peakcount_, back_inserter(out));
          }
          else
          {
            std::sort(peaks_in_window.begin(), peaks_in_window.end(), reverseComparator(typename SpectrumType::PeakType::IntensityLess()));
            copy(peaks_in_window.begin(), peaks_in_window.end(), back_inserter(out));
          }

          peaks_in_window.clear(false);
          peaks_in_window.push_back(spectrum[i]);
        }
      }

      if (peaks_in_window.empty()) // last window is empty -> no special handling needed
      {
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

      // Note that the last window might be much smaller than windowsize.
      // Therefor the number of peaks copied from this window should be adapted accordingly.
      // Otherwise a lot of noise peaks are copied from each end of a spectrum.

      double last_window_size = peaks_in_window.back().getMZ() - window_start;
      double last_window_size_fraction = last_window_size / windowsize_;
      Size last_window_peakcount = last_window_size_fraction * peakcount_;

      if (last_window_peakcount) // handle single peak in last window (will produce no proper fraction)
      {
        last_window_peakcount = 1;
      }

      // sort for last_window_peakcount highest peaks
      std::partial_sort(peaks_in_window.begin(), peaks_in_window.begin() + last_window_peakcount, peaks_in_window.end(), reverseComparator(typename SpectrumType::PeakType::IntensityLess()));

      if (peaks_in_window.size() > last_window_peakcount)
      {
        std::copy(peaks_in_window.begin(), peaks_in_window.begin() + last_window_peakcount, back_inserter(out));
      }
      else
      {
        std::copy(peaks_in_window.begin(), peaks_in_window.end(), std::back_inserter(out));
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

#endif //OPENMS_FILTERING_TRANSFORMERS_WINDOWMOWER_H

