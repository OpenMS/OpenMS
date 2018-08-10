// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

  template <typename PeakContainerT>
  void removePeaks(
    PeakContainerT& p,
    const double pos_start,
    const double pos_end,
    const bool ignoreDataArrays = false
  )
  {
    typename PeakContainerT::iterator it_start = p.PosBegin(pos_start);
    typename PeakContainerT::iterator it_end = p.PosEnd(pos_end);
    if (!ignoreDataArrays)
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
    p.erase(it_end, p.end());
    p.erase(p.begin(), it_start);
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
} // namespace OpenMS


