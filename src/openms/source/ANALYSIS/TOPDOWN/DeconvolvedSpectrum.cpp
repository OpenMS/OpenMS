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

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h"

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

  void DeconvolvedSpectrum::writeDeconvolvedMasses(std::fstream& fs,
                                                   const String& file_name,
                                                   const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                                   const bool write_detail)
  {
    if (empty())
    {
      return;
    }

    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      const double mono_mass = pg.getMonoMass();
      const double avg_mass = pg.getMonoMass() + avg.getAverageMassDelta(mono_mass);
      const double intensity = pg.getIntensity();

      auto charge_range = pg.getAbsChargeRange();
      int min_charge = pg.isPositive() ? std::get<0>(charge_range) : std::get<1>(charge_range);
      int max_charge = pg.isPositive() ? std::get<1>(charge_range) : std::get<0>(charge_range);

      fs << file_name << "\t" << pg.getScanNumber() << "\t"
         << std::to_string(spec_.getRT()) << "\t"
         << size() << "\t"
         << std::to_string(avg_mass) << "\t" << std::to_string(mono_mass) << "\t" << intensity << "\t"
         << min_charge << "\t" << max_charge << "\t"
         << pg.size() << "\t";

      if (write_detail)
      {
        fs << std::fixed << std::setprecision(2);
        for (auto& p : pg)
        {
          fs << p.mz << " ";
        }
        fs << ";\t";

        fs << std::fixed << std::setprecision(1);
        for (auto& p : pg)
        {
          fs << p.intensity << " ";
        }
        fs << ";\t";
        fs << std::setprecision(-1);

        for (auto& p : pg)
        {
          fs << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }
        fs << ";\t";
        for (auto& p : pg)
        {
          fs << p.getUnchargedMass() << " ";
        }
        fs << ";\t";
        for (auto& p : pg)
        {
          fs << p.isotopeIndex << " ";
        }
        fs << ";\t";

        for (auto& p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          double mass_error =
              (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) -
               p.mz) /
              p.mz;
          fs << 1e6 * mass_error << " ";
        }
        fs << ";\t";
      }
      if (spec_.getMSLevel() > 1)
      {
        //PrecursorScanNum	PrecursorMz	PrecursorIntensity PrecursorCharge	PrecursorMonoMass		PrecursorQScore
        fs << precursor_scan_number_ << "\t" << std::to_string(precursor_peak_.getMZ()) << "\t"
           << precursor_peak_.getIntensity() << "\t"
           << precursor_peak_.getCharge()
           << "\t";

        if (precursor_peak_group_.empty())
        {
          fs << "nan\tnan\tnan\t";
        }
        else
        {
          fs << precursor_peak_group_.getChargeSNR(precursor_peak_.getCharge()) << "\t"
             << std::to_string(precursor_peak_group_.getMonoMass()) << "\t"
             << precursor_peak_group_.getQScore() << "\t";
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeScore() << "\t";

      auto max_qscore_mz_range = pg.getMaxQScoreMzRange();
      fs << pg.getSNR() << "\t" << pg.getChargeSNR(pg.getRepAbsCharge()) << "\t"
         << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t"
         << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t"
         << pg.getQScore() << "\t" << std::setprecision(-1);//

      for (int i = std::get<0>(charge_range); i <= std::get<1>(charge_range); i++)
      {

        fs << pg.getChargeIntensity(i);

        if (i < std::get<1>(charge_range))
        {
          fs << ";";
        }
      }
      fs << "\t";
      int isotope_end_index = 0;

      for (auto& p : pg)
      {
        isotope_end_index = isotope_end_index < p.isotopeIndex ? p.isotopeIndex : isotope_end_index;
      }
      auto per_isotope_intensity = std::vector<double>(isotope_end_index + 1, .0);
      for (auto& p : pg)
      {
        per_isotope_intensity[p.isotopeIndex] += p.intensity;
      }

      for (int i = 0; i <= isotope_end_index; i++)
      {
        fs << per_isotope_intensity[i];
        if (i < isotope_end_index)
        {
          fs << ";";
        }
      }

      fs << "\n";
    }
  }


  void DeconvolvedSpectrum::writeDeconvolvedMassesHeader(std::fstream& fs, const int ms_level, const bool detail)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
    else
    {
      if (ms_level == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               //"PeakMZs\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\tPeakIntensities\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorSNR\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tChargeSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
  }

  void DeconvolvedSpectrum::writeTopFD(std::fstream& fs,
                                       const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                       const double snr_threshold,
                                       const double decoy_harmonic_factor,
                                       const double decoy_precursor_offset)//, fstream& fsm, fstream& fsp)
  {
    UInt ms_level = spec_.getMSLevel();

    if (ms_level > 1)
    {
      if (precursor_peak_group_.empty() ||
          precursor_peak_group_.getChargeSNR(precursor_peak_.getCharge()) < snr_threshold)
      {
        return;
      }
    }

    if (size() < topFD_min_peak_count_)
    {
      return;
    }

    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << scan_number_ << "\n"
       << "SCANS=" << scan_number_ << "\n"
       << "RETENTION_TIME=" << spec_.getRT() << "\n";

    if (ms_level > 1)
    {
      fs << "ACTIVATION=" << activation_method_ << "\n";
      if (!precursor_peak_group_.empty())
      {
        fs << "MS_ONE_ID=" << precursor_scan_number_ << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number_ << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak_.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << (int) (precursor_peak_.getCharge() * decoy_harmonic_factor) << "\n"
           << "PRECURSOR_MASS="
           << std::to_string(precursor_peak_group_.getMonoMass() * decoy_harmonic_factor + decoy_precursor_offset) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak_.getIntensity() << "\n";
      }
      else
      {
        double average_mass =
            (precursor_peak_.getMZ() -
             FLASHDeconvHelperStructs::getChargeMass(precursor_peak_.getCharge() > 0)) *
            abs(precursor_peak_.getCharge() * decoy_harmonic_factor);
        double mono_mass = average_mass - avg.getAverageMassDelta(average_mass) + decoy_precursor_offset;
        fs << "MS_ONE_ID=" << precursor_scan_number_ << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number_ << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak_.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << (int) (precursor_peak_.getCharge() * decoy_harmonic_factor) << "\n"
           << "PRECURSOR_MASS=" << std::to_string(mono_mass) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak_.getIntensity() << "\n";
      }
    }
    fs << std::setprecision(-1);

    double isotope_score_threshold = 0;
    std::vector<double> isotope_scores;

    if (size() > topFD_max_peak_count_)// max peak count for TopPic = 500
    {
      isotope_scores.reserve(size());
      for (auto& pg : *this)
      {
        isotope_scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(isotope_scores.begin(), isotope_scores.end());
      isotope_score_threshold = isotope_scores[isotope_scores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(isotope_scores);
    }

    int size = 0;
    for (auto& pg : *this)
    {
      if (pg.getIsotopeCosine() < isotope_score_threshold)
      {
        continue;
      }
      if (size >= topFD_max_peak_count_)
      {
        break;
      }
      size++;
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t"
         << (pg.isPositive() ? std::get<1>(pg.getAbsChargeRange()) : -std::get<1>(pg.getAbsChargeRange()))
         << "\n";
      fs << std::setprecision(-1);
      if (size >= topFD_max_peak_count_)
      {
        break;
      }

      // peaks of different charges are separately recorded even if they represent the same mass..
      /*std::set<int> charges;
      for (auto& peaks : pg)
      {
        charges.insert(peaks.abs_charge);
      }
      for (int charge : charges)
      {
        if (pg.getChargeIntensity(charge) <= 0)
        {
          continue;
        }
        size++;
        fs << std::fixed << std::setprecision(2);
        fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getChargeIntensity(charge) << "\t"
           << (pg.isPositive() ? charge : -charge)
           << "\n";
        fs << std::setprecision(-1);
        if (size >= topFD_max_peak_count_)
        {
          break;
        }
      }*/
    }

    fs << "END IONS\n\n";
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
    std::sort(peak_groups.begin(), peak_groups.end(), [](const PeakGroup & p1, const PeakGroup & p2){return p1.getQScore() < p2.getQScore();});
  }


}
