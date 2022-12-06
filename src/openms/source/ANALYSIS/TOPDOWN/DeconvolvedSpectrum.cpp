// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
  DeconvolvedSpectrum::DeconvolvedSpectrum(const int scan_number) : scan_number_(scan_number)
  {
  }

  MSSpectrum DeconvolvedSpectrum::toSpectrum(const int to_charge, double tol, bool retain_undeconvolved)
  {
    auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);
    if ((spec_.getMSLevel() > 1 && precursor_peak_group_.empty()) || empty())
    {
      return out_spec;
    }
    double charge_mass_offset = (double)abs(to_charge) * FLASHDeconvHelperStructs::getChargeMass(to_charge >= 0);
    std::unordered_set<double> deconvolved_mzs;
    std::stringstream val {};

    val << "tol=" << tol << ";massoffset=" << std::to_string(charge_mass_offset) << ";chargemass=" << std::to_string(FLASHDeconvHelperStructs::getChargeMass(peak_groups_[0].isPositive()))
        << ";peaks=";
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }

      out_spec.emplace_back(pg.getMonoMass() + charge_mass_offset, pg.getIntensity());
      auto [z1, z2] = pg.getAbsChargeRange();
      int min_iso = -1, max_iso = 0;

      for (auto& p : pg)
      {
        min_iso = min_iso < 0 ? p.isotopeIndex : std::min(min_iso, p.isotopeIndex);
        max_iso = std::max(max_iso, p.isotopeIndex);
      }
      val << z1 << ":" << z2 << "," << min_iso << ":" << max_iso << ";";

      if (retain_undeconvolved)
      {
        for (auto& p : pg)
        {
          deconvolved_mzs.insert(p.mz);
        }
      }
    }
    out_spec.setMetaValue("DeconvMassInfo", val.str());

    if (retain_undeconvolved)
    {
      for (auto& p : spec_)
      {
        if (deconvolved_mzs.find(p.getMZ()) != deconvolved_mzs.end()) // if p is deconvolved
        {
          continue;
        }
        out_spec.emplace_back(p.getMZ() + charge_mass_offset - FLASHDeconvHelperStructs::getChargeMass(to_charge >= 0), p.getIntensity());
      }
    }
    out_spec.sortByPosition();
    if (!precursor_peak_group_.empty() && !precursor_peak_.empty())
    {
      Precursor precursor(spec_.getPrecursors()[0]);
      precursor.setCharge(to_charge);
      precursor.setMZ(precursor_peak_group_.getMonoMass() + charge_mass_offset);
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
      return *(new PeakGroup());
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

  const Precursor::ActivationMethod& DeconvolvedSpectrum::getActivationMethod() const
  {
    return activation_method_;
  }

  void DeconvolvedSpectrum::setPrecursor(const Precursor& precursor)
  {
    precursor_peak_ = precursor;
  }

  void DeconvolvedSpectrum::setPrecursorIntensity(const float i)
  {
    precursor_peak_.setIntensity(i);
  }

  void DeconvolvedSpectrum::setActivationMethod(const Precursor::ActivationMethod& method)
  {
    activation_method_ = method;
  }

  void DeconvolvedSpectrum::setPrecursorPeakGroup(const PeakGroup& pg)
  {
    precursor_peak_group_ = pg;
  }

  void DeconvolvedSpectrum::setOriginalSpectrum(const MSSpectrum& spec)
  {
    spec_ = spec;
  }


  void DeconvolvedSpectrum::setPrecursorScanNumber(const int scan_number)
  {
    precursor_scan_number_ = scan_number;
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::begin() const noexcept
  {
    return peak_groups_.begin();
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::end() const noexcept
  {
    return peak_groups_.end();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::begin() noexcept
  {
    return peak_groups_.begin();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::end() noexcept
  {
    return peak_groups_.end();
  }

  const PeakGroup& DeconvolvedSpectrum::operator[](const Size i) const
  {
    return peak_groups_[i];
  }

  PeakGroup& DeconvolvedSpectrum::operator[](const Size i)
  {
    return peak_groups_[i];
  }

  void DeconvolvedSpectrum::push_back(const PeakGroup& pg)
  {
    peak_groups_.push_back(pg);
  }

  Size DeconvolvedSpectrum::size() const noexcept
  {
    return peak_groups_.size();
  }

  void DeconvolvedSpectrum::clear()
  {
    std::vector<PeakGroup>().swap(peak_groups_);
  }

  void DeconvolvedSpectrum::reserve(Size n)
  {
    peak_groups_.reserve(n);
  }

  bool DeconvolvedSpectrum::empty() const
  {
    return peak_groups_.empty();
  }

  void DeconvolvedSpectrum::setPeakGroups(std::vector<PeakGroup>& x)
  {
    std::vector<PeakGroup>().swap(peak_groups_);
    peak_groups_ = x;
  }

  void DeconvolvedSpectrum::sort()
  {
    std::sort(peak_groups_.begin(), peak_groups_.end());
  }

  void DeconvolvedSpectrum::sortByQScore()
  {
    std::sort(peak_groups_.begin(), peak_groups_.end(), [](const PeakGroup& p1, const PeakGroup& p2) { return p1.getQScore() > p2.getQScore(); });
  }

  void DeconvolvedSpectrum::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                   std::vector<DeconvolvedSpectrum>& deconvolved_decoy_spectra) // per ms level + precursor update as well.
  {
    std::map<uint, std::vector<float>> tscore_map; // per ms level

    std::map<uint, std::vector<float>> dscore_iso_decoy_map;
    std::map<uint, std::map<float, float>> qscore_iso_decoy_map; // maps for isotope decoy only qvalues

    std::map<uint, std::vector<float>> dscore_noise_decoy_map;
    std::map<uint, std::map<float, float>> qscore_noise_decoy_map; // maps for noise decoy only qvalues

    std::map<uint, std::vector<float>> dscore_charge_decoy_map;
    std::map<uint, std::map<float, float>> qscore_charge_decoy_map; // maps for charge decoy only qvalues

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      if (tscore_map.find(ms_level) == tscore_map.end())
      {
        tscore_map[ms_level] = std::vector<float>();
        dscore_noise_decoy_map[ms_level] = std::vector<float>();
        qscore_noise_decoy_map[ms_level] = std::map<float, float>();
        dscore_iso_decoy_map[ms_level] = std::vector<float>();
        qscore_iso_decoy_map[ms_level] = std::map<float, float>();
        dscore_charge_decoy_map[ms_level] = std::vector<float>();
        qscore_charge_decoy_map[ms_level] = std::map<float, float>();
      }
      for (auto& pg : deconvolved_spectrum)
      {
        tscore_map[ms_level].push_back(pg.getQScore());
      }
    }
    for (auto& decoy_deconvolved_spectrum : deconvolved_decoy_spectra)
    {
      uint ms_level = decoy_deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : decoy_deconvolved_spectrum)
      {
        if (pg.getDecoyFlag() == PeakGroup::DecoyFlag::isotope_decoy)
        {
          dscore_iso_decoy_map[ms_level].push_back(pg.getQScore());
        }
        else if (pg.getDecoyFlag() == PeakGroup::DecoyFlag::noise_decoy)
        {
          dscore_noise_decoy_map[ms_level].push_back(pg.getQScore());
        }
        else if (pg.getDecoyFlag() == PeakGroup::DecoyFlag::charge_decoy)
        {
          dscore_charge_decoy_map[ms_level].push_back(pg.getQScore());
        }
      }
    }

    for (auto& titem : tscore_map)
    {
      uint ms_level = titem.first;
      auto& tscore = titem.second;
      auto& dscore_iso = dscore_iso_decoy_map[ms_level];
      auto& dscore_charge = dscore_charge_decoy_map[ms_level];
      auto& dscore_noise = dscore_noise_decoy_map[ms_level];

      std::sort(tscore.begin(), tscore.end());
      std::sort(dscore_iso.begin(), dscore_iso.end());
      std::sort(dscore_noise.begin(), dscore_noise.end());
      std::sort(dscore_charge.begin(), dscore_charge.end());

      auto& map_charge = qscore_charge_decoy_map[ms_level];
      float tmp_q_charge = 1;
      for (size_t i = 0; i < tscore.size(); i++)
      {
        float ts = tscore[i];
        size_t dindex = dscore_charge.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_charge.begin(), dscore_charge.end(), ts), dscore_charge.end());
        size_t tindex = tscore.size() - i;

        tmp_q_charge = std::min(tmp_q_charge, ((float)dindex / (float)(dindex + tindex)));
        map_charge[ts] = tmp_q_charge;
      }

      auto& map_iso = qscore_iso_decoy_map[ms_level];
      float tmp_q_iso = 1;
      for (size_t i = 0; i < tscore.size(); i++)
      {
        float ts = tscore[i];
        size_t dindex = dscore_iso.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_iso.begin(), dscore_iso.end(), ts), dscore_iso.end());
        size_t tindex = tscore.size() - i;

        tmp_q_iso = std::min(tmp_q_iso, ((float)dindex / (float)(dindex + tindex) * (1 - map_charge[ts])));
        map_iso[ts] = tmp_q_iso;
      }

      auto& map_noise = qscore_noise_decoy_map[ms_level];
      float tmp_q_noise = 1;

      for (size_t i = 0; i < tscore.size(); i++)
      {
        float ts = tscore[i];
        size_t dindex = dscore_noise.size() == 0 ? 0 : std::distance(std::lower_bound(dscore_noise.begin(), dscore_noise.end(), ts), dscore_noise.end());
        size_t tindex = tscore.size() - i;

        tmp_q_noise = std::min(tmp_q_noise, ((float)dindex / (float)(dindex + tindex)));
        map_noise[ts] = tmp_q_noise;
      }
    }

    for (auto& titem : tscore_map)
    {
      uint ms_level = titem.first;
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }
        auto& map_iso = qscore_iso_decoy_map[ms_level];
        auto& map_noise = qscore_noise_decoy_map[ms_level];
        auto& map_charge = qscore_charge_decoy_map[ms_level];

        for (auto& pg : deconvolved_spectrum)
        {
          //
          pg.setQvalue(map_charge[pg.getQScore()], PeakGroup::DecoyFlag::charge_decoy);
          pg.setQvalue(map_noise[pg.getQScore()], PeakGroup::DecoyFlag::noise_decoy);
          pg.setQvalue(map_iso[pg.getQScore()], PeakGroup::DecoyFlag::isotope_decoy);

          if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() > 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
          {
            float qs = deconvolved_spectrum.getPrecursorPeakGroup().getQScore();
            auto& pmap_iso = qscore_iso_decoy_map[ms_level - 1];
            auto& pmap_charge = qscore_charge_decoy_map[ms_level - 1];
            auto& pmap_noise = qscore_noise_decoy_map[ms_level - 1];
            deconvolved_spectrum.precursor_peak_group_.setQvalue(pmap_iso[qs], PeakGroup::DecoyFlag::isotope_decoy);
            deconvolved_spectrum.precursor_peak_group_.setQvalue(pmap_noise[qs], PeakGroup::DecoyFlag::noise_decoy);
            deconvolved_spectrum.precursor_peak_group_.setQvalue(pmap_charge[qs], PeakGroup::DecoyFlag::charge_decoy);
          }
        }
      }
    }
  }
}