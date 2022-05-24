// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>

namespace OpenMS
{
  DeconvolvedSpectrum::DeconvolvedSpectrum(const MSSpectrum& spectrum, const int scan_number) :
      scan_number_(scan_number)
  {
    spec_ = spectrum;
  }

  MSSpectrum DeconvolvedSpectrum::toSpectrum(const int to_charge)
  {
    auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);
    if(spec_.getMSLevel() > 1 && precursor_peak_group_.empty())
    {
      return out_spec;
    }

    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      out_spec.emplace_back(pg.getMonoMass(), pg.getIntensity());
    }
    if (!precursor_peak_group_.empty() && !spec_.getPrecursors().empty())
    {
      Precursor precursor(spec_.getPrecursors()[0]);
      //precursor.setCharge((precursor_peak_group_.isPositive() ?
      //                     precursor_peak_group_.getRepAbsCharge() :
      //                     -precursor_peak_group_.getRepAbsCharge()));//getChargeMass
      precursor.setCharge(to_charge);
      precursor.setMZ(precursor_peak_group_.getMonoMass() +
                      to_charge * FLASHDeconvHelperStructs::getChargeMass(to_charge >= 0));
      precursor.setIntensity(precursor_peak_group_.getIntensity());
      out_spec.getPrecursors().clear();
      out_spec.getPrecursors().emplace_back(precursor);
    }
    return out_spec;
  }


  const MSSpectrum& DeconvolvedSpectrum::getOriginalSpectrum() const
  {
    return spec_;
  }

  const PeakGroup& DeconvolvedSpectrum::getPrecursorPeakGroup() const
  {
    if (precursor_peak_group_.empty())
    {
      return PeakGroup();
    }
    return precursor_peak_group_;
  }

  int DeconvolvedSpectrum::getPrecursorCharge() const
  {
    return precursor_peak_.getCharge();
  }

  double DeconvolvedSpectrum::getCurrentMaxMass(const double max_mass) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_mass;
    }
    return precursor_peak_group_.getMonoMass();
  }

  double DeconvolvedSpectrum::getCurrentMinMass(const double min_mass) const
  {
    if (spec_.getMSLevel() == 1)
    {
      return min_mass;
    }
    return 50.0;
  }

  int DeconvolvedSpectrum::getCurrentMaxAbsCharge(const int max_abs_charge) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_abs_charge;
    }
    return abs(precursor_peak_.getCharge());
  }

  const Precursor& DeconvolvedSpectrum::getPrecursor() const
  {
    return precursor_peak_;
  }

  int DeconvolvedSpectrum::getScanNumber() const
  {
    return scan_number_;
  }

  int DeconvolvedSpectrum::getPrecursorScanNumber() const
  {
    return precursor_scan_number_;
  }


  const String& DeconvolvedSpectrum::getActivationMethod() const
  {
    return activation_method_;
  }


  void DeconvolvedSpectrum::setPrecursor(const Precursor& precursor)
  {
    precursor_peak_ = precursor;
  }

  void DeconvolvedSpectrum::setPrecursorIntensity(const double i)
  {
    precursor_peak_.setIntensity(i);
  }

  void DeconvolvedSpectrum::setActivationMethod(const String& method)
  {
    activation_method_ = method;
  }

  void DeconvolvedSpectrum::setPrecursorPeakGroup(const PeakGroup& pg)
  {
    precursor_peak_group_ = pg;
  }

  void DeconvolvedSpectrum::setPrecursorScanNumber(const int scan_number)
  {
    precursor_scan_number_ = scan_number;
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::begin() const noexcept
  {
    return peak_groups.begin();
  }
  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::end() const noexcept
  {
    return peak_groups.end();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::begin() noexcept
  {
    return peak_groups.begin();
  }
  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::end() noexcept
  {
    return peak_groups.end();
  }

  const PeakGroup& DeconvolvedSpectrum::operator[](const Size i) const
  {
    return peak_groups[i];
  }

  void DeconvolvedSpectrum::push_back(const PeakGroup& pg)
  {
    peak_groups.push_back(pg);
  }
  Size DeconvolvedSpectrum::size() const noexcept
  {
    return peak_groups.size();
  }
  void DeconvolvedSpectrum::clear()
  {
    peak_groups.clear();
  }
  void DeconvolvedSpectrum::reserve(Size n)
  {
    peak_groups.reserve(n);
  }
  bool DeconvolvedSpectrum::empty() const
  {
    return peak_groups.empty();
  }
  void DeconvolvedSpectrum::swap (std::vector<PeakGroup>& x)
  {
    peak_groups.swap(x);
  }

  void DeconvolvedSpectrum::sort()
  {
    std::sort(peak_groups.begin(), peak_groups.end());
  }

  void DeconvolvedSpectrum::sortByQScore()
  {
    std::sort(peak_groups.begin(), peak_groups.end(), [](const PeakGroup & p1, const PeakGroup & p2){return p1.getQScore() > p2.getQScore();});
  }


}
