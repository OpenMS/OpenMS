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
// $Maintainer: Hannes Roest, Witold Wolski $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/Macros.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace OpenSwath
{

  void integrateWindows(const OpenSwath::SpectrumPtr spectrum,
                        const std::vector<double> & windowsCenter, double width,
                        std::vector<double> & integratedWindowsIntensity,
                        std::vector<double> & integratedWindowsMZ,
                        bool remZero)
  {
    std::vector<double>::const_iterator beg = windowsCenter.begin();
    std::vector<double>::const_iterator end = windowsCenter.end();
    double mz, intensity;
    for (; beg != end; ++beg)
    {
      double left = *beg - width / 2.0;
      double right = *beg + width / 2.0;
      if (integrateWindow(spectrum, left, right, mz, intensity, false))
      {
        integratedWindowsIntensity.push_back(intensity);
        integratedWindowsMZ.push_back(mz);
      }
      else if (!remZero)
      {
        integratedWindowsIntensity.push_back(0.);
        integratedWindowsMZ.push_back(*beg);
      }
      else
      {
      }
    }
  }

  /// integrate all masses in window
  bool integrateWindow(const OpenSwath::SpectrumPtr spectrum, double mz_start, double mz_end,
                       double & mz, double & intensity, bool centroided)
  {
    OPENSWATH_PRECONDITION( std::adjacent_find(spectrum->getMZArray()->data.begin(),
            spectrum->getMZArray()->data.end(), std::greater<double>()) == spectrum->getMZArray()->data.end(),
          "Precondition violated: m/z vector needs to be sorted!" )

    intensity = 0;
    if (!centroided)
    {
      // get the weighted average for noncentroided data.
      // TODO this is not optimal if there are two peaks in this window (e.g. if the window is too large)
      typedef std::vector<double>::const_iterator itType;
      mz = 0;
      intensity = 0;

      itType mz_arr_end = spectrum->getMZArray()->data.end();
      itType int_it = spectrum->getIntensityArray()->data.begin();

      // this assumes that the spectra are sorted!
      itType mz_it = std::lower_bound(spectrum->getMZArray()->data.begin(),
        spectrum->getMZArray()->data.end(), mz_start);
      itType mz_it_end = std::lower_bound(mz_it, mz_arr_end, mz_end);

      // also advance intensity iterator now
      std::iterator_traits< itType >::difference_type iterator_pos = std::distance((itType)spectrum->getMZArray()->data.begin(), mz_it);
      std::advance(int_it, iterator_pos);

      for (; mz_it != mz_it_end; ++mz_it, ++int_it)
      {
        intensity += (*int_it);
        mz += (*int_it) * (*mz_it);
      }

      if (intensity > 0.)
      {
        mz /= intensity;
        return true;
      }
      else
      {
        mz = -1;
        intensity = 0;
        return false;
      }

    }
    else
    {
      // not implemented
      throw "Not implemented";
    }
  }

}
