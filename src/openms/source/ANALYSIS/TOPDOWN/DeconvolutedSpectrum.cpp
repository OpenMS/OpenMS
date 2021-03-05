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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"
#include "include/OpenMS/ANALYSIS/TOPDOWN/QScore.h"
namespace OpenMS
{
  DeconvolutedSpectrum::DeconvolutedSpectrum(const MSSpectrum &spectrum, const int scan_number) :
      scan_number_(scan_number)
  {
    spec_ = spectrum;
  }

  MSSpectrum DeconvolutedSpectrum::toSpectrum(const int mass_charge)
  {
    auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);
    for (auto &pg : *this)
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
      precursor.setCharge((precursor_peak_group_.isPositive() ?
                           precursor_peak_group_.getRepAbsCharge() :
                           -precursor_peak_group_.getRepAbsCharge()));//getChargeMass
      precursor.setMZ(precursor_peak_group_.getMonoMass() +
                      mass_charge * FLASHDeconvHelperStructs::getChargeMass(mass_charge >= 0));
      precursor.setIntensity(precursor_peak_group_.getIntensity());
      out_spec.getPrecursors().clear();
      out_spec.getPrecursors().emplace_back(precursor);
    }
    return out_spec;
  }

  void DeconvolutedSpectrum::writeDeconvolutedMasses(std::fstream &fs,
                                                     const String &file_name,
                                                     const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                                     const bool write_detail)//, fstream& fsm, fstream& fsp)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg : *this)
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
        for (auto &p : pg)
        {
          fs << p.mz << " ";
        }
        fs << ";\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p : pg)
        {
          fs << p.intensity << " ";
        }
        fs << ";\t";
        fs << std::setprecision(-1);

        for (auto &p : pg)
        {
          fs << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }
        fs << ";\t";
        for (auto &p : pg)
        {
          fs << p.getUnchargedMass() << " ";
        }
        fs << ";\t";
        for (auto &p : pg)
        {
          fs << p.isotopeIndex << " ";
        }
        fs << ";\t";

        for (auto &p : pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          double mass_error =
              (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) - p.mz) / p.mz;
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
          fs << "nan\tnan\t";
        }
        else
        {
          fs << std::to_string(precursor_peak_group_.getMonoMass()) << "\t"
             << precursor_peak_group_.getQScore() << "\t";
        }
      }
      fs << pg.getIsotopeCosine() << "\t" << pg.getChargeScore() << "\t";

      auto max_qscore_mz_range = pg.getMaxQScoreMzRange();
      fs << pg.getSNR() << "\t"
         << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge()) << "\t"
         << std::to_string(std::get<0>(max_qscore_mz_range)) << "\t"
         << std::to_string(std::get<1>(max_qscore_mz_range)) << "\t"
         << pg.getQScore() << "\t" << std::setprecision(-1); //

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

      for (auto &p : pg)
      {
        isotope_end_index = isotope_end_index < p.isotopeIndex ? p.isotopeIndex : isotope_end_index;
      }
      auto per_isotope_intensity = std::vector<double>(isotope_end_index + 1, .0);
      for (auto &p : pg)
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


  void DeconvolutedSpectrum::writeDeconvolutedMassesHeader(std::fstream &fs, const int ms_level, const bool detail)
  {
    if (detail)
    {
      if (ms_level == 1)
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\tPeakMZs\tPeakIntensities\tPeakCharges\tPeakMasses\tPeakIsotopeIndices\tPeakPPMErrors\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
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
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";

      }
      else
      {
        fs
            << "FileName\tScanNum\tRetentionTime\tMassCountInSpec\tAverageMass\tMonoisotopicMass\t"
               "SumIntensity\tMinCharge\tMaxCharge\t"
               "PeakCount\t"
               "PrecursorScanNum\tPrecursorMz\tPrecursorIntensity\tPrecursorCharge\tPrecursorMonoisotopicMass\tPrecursorQScore\t"
               "IsotopeCosine\tChargeScore\tMassSNR\tRepresentativeCharge\tRepresentativeMzStart\tRepresentativeMzEnd\tQScore\tPerChargeIntensity\tPerIsotopeIntensity\n";
      }
    }
  }

  /*  void DeconvolutedSpectrum::writeAttCsvHeader(std::fstream& fs)
    {
      fs
          << "ScanNumber,RetentionTime,PrecursorScanNumber,Charge,ChargeSNR,PeakIntensity,EnvIntensity,EnvIsotopeCosine,PeakMz,"
             "MonoMass,MassSNR,IsotopeCosine,MassIntensity,QScore,Class\n";
    }*/

  void DeconvolutedSpectrum::writeThermoInclusionHeader(std::fstream &fs)
  {
    fs << "Compound,Formula,Adduct,m/z,z,t start (min),t stop (min),Isolation Window (m/z),Normalized AGC Target (%)\n";
  }


  void DeconvolutedSpectrum::clearPeakGroupsChargeInfo()
  {
    for (auto &pg: *this)
    {
      pg.clearChargeInfo();
    }
  }

  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs,
                                        const int id,
                                        const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                        const double harmonic_factor,
                                        const double precursor_offset)//, fstream& fsm, fstream& fsp)
  {

    UInt ms_level = spec_.getMSLevel();

    if (ms_level > 1)
    {
      if (precursor_peak_group_.empty())
      {
        return;
      }
    }

    fs << std::fixed << std::setprecision(2);
    fs << "BEGIN IONS\n"
       << "ID=" << id << "\n"
       << "SCANS=" << scan_number_ << "\n"
       << "RETENTION_TIME=" << spec_.getRT() << "\n";

    fs << "ACTIVATION=" << activation_method_ << "\n";

    if (ms_level > 1)
    {
      if (!precursor_peak_group_.empty())
      {
        fs << "MS_ONE_ID=" << precursor_scan_number_ << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number_ << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak_.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << (int) (precursor_peak_.getCharge() * harmonic_factor) << "\n"
           << "PRECURSOR_MASS="
           << std::to_string(precursor_peak_group_.getMonoMass() * harmonic_factor + precursor_offset) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak_.getIntensity() << "\n";
      }
      else
      {
        double average_mass =
            (precursor_peak_.getMZ() - FLASHDeconvHelperStructs::getChargeMass(precursor_peak_.getCharge() > 0)) *
            abs(precursor_peak_.getCharge() * harmonic_factor);
        double mono_mass = average_mass - avg.getAverageMassDelta(average_mass) + precursor_offset;
        fs << "MS_ONE_ID=" << precursor_scan_number_ << "\n"
           << "MS_ONE_SCAN=" << precursor_scan_number_ << "\n"
           << "PRECURSOR_MZ="
           << std::to_string(precursor_peak_.getMZ()) << "\n"
           << "PRECURSOR_CHARGE=" << (int) (precursor_peak_.getCharge() * harmonic_factor) << "\n"
           << "PRECURSOR_MASS=" << std::to_string(mono_mass) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak_.getIntensity() << "\n";
      }
    }
    fs << std::setprecision(-1);

    double isotope_score_threshold = 0;
    std::vector<double> isotope_scores;

    if (size() > 500)// max peak count for TopPic
    {
      isotope_scores.reserve(size());
      for (auto &pg : *this)
      {
        isotope_scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(isotope_scores.begin(), isotope_scores.end());
      isotope_score_threshold = isotope_scores[isotope_scores.size() - 500];
      std::vector<double>().swap(isotope_scores);
    }

    int size = 0;
    for (auto &pg : *this)
    {
      if (pg.getIsotopeCosine() < isotope_score_threshold)
      {
        continue;
      }
      if (size >= 500)
      {
        break;
      }

      size++;
      fs << std::fixed << std::setprecision(2);
      fs << std::to_string(pg.getMonoMass()) << "\t" << pg.getIntensity() << "\t"
         << (pg.isPositive() ? pg.getRepAbsCharge() : -pg.getRepAbsCharge())
         //  << "\t" << log10(pg.precursorSNR+1e-10) << "\t" << log10(pg.precursorTotalSNR+1e-10)
         //  << "\t" << log10(pg.maxSNR + 1e-10) << "\t" << log10(pg.total_snr_ + 1e-10)
         << "\n";
      fs << std::setprecision(-1);
    }

    fs << "END IONS\n\n";
  }

  bool DeconvolutedSpectrum::registerPrecursor(const std::vector<DeconvolutedSpectrum> &survey_scans,
                                               const std::map<int, std::vector<std::vector<double>>> &precursor_map_for_real_time_acquisition)
  {
    bool is_positive = true; // TODO update..
    precursor_peak_.setIntensity(.0);

    if (!precursor_map_for_real_time_acquisition.empty())
    {
      for (auto &precursor: spec_.getPrecursors())
      {
        //double max_precurosr_intensity = .0;
        for (auto &activation_method :  precursor.getActivationMethods())
        {
          activation_method_ = Precursor::NamesOfActivationMethodShort[activation_method];
          break;
        }
        if (precursor_peak_.getIntensity() > precursor.getIntensity())
        {
          continue;
        }
        precursor_peak_ = precursor;
      }

      double start_mz = precursor_peak_.getIsolationWindowLowerOffset() > 100.0 ?
                        precursor_peak_.getIsolationWindowLowerOffset() :
                        -precursor_peak_.getIsolationWindowLowerOffset() + precursor_peak_.getMZ();
      double end_mz = precursor_peak_.getIsolationWindowUpperOffset() > 100.0 ?
                      precursor_peak_.getIsolationWindowUpperOffset() :
                      precursor_peak_.getIsolationWindowUpperOffset() + precursor_peak_.getMZ();
      // ms1 scan -> mass, charge ,score, mz range
      double center_mz = (start_mz + end_mz) / 2.0;
      double mz_tol = (end_mz - center_mz) + .1;

      for (int i = survey_scans.size() - 1; i >= 0; i--)
      {
        auto precursor_spectrum = survey_scans[i];
        if (precursor_spectrum.empty())
        {
          continue;
        }
        int index = precursor_spectrum.spec_.findNearest(center_mz, mz_tol);
        if (index < 0)
        {
          continue;
        }
        if (precursor_peak_.getIntensity() > precursor_spectrum.spec_[index].getIntensity())
        {
          continue;
        }
        precursor_peak_.setIntensity(precursor_spectrum.spec_[index].getIntensity());
        precursor_peak_.setMZ(precursor_spectrum.spec_[index].getMZ());
      }

      for (auto map = precursor_map_for_real_time_acquisition.upper_bound(scan_number_);
           map != precursor_map_for_real_time_acquisition.begin();
           map--)
      {
        precursor_scan_number_ = map->first;
        //std::cout<<scan_number_ << " * " << precursor_scan_number_ <<std::endl;
        if (map != precursor_map_for_real_time_acquisition.end())
        {
          for (auto &smap : map->second)
          {
            //
            if (abs(start_mz - smap[3]) < .1 || abs(end_mz - smap[4]) < .1)
            {
              LogMzPeak precursor_log_mz_peak(precursor_peak_, is_positive);
              precursor_log_mz_peak.abs_charge = (int) smap[1];
              precursor_log_mz_peak.isotopeIndex = 0;
              precursor_log_mz_peak.mass = smap[0];
              precursor_peak_.setCharge(precursor_log_mz_peak.abs_charge);
              precursor_peak_group_.push_back(precursor_log_mz_peak);
              precursor_peak_group_.setQScore(smap[2]);
              precursor_peak_group_.setRepAbsCharge((int) smap[1]);
              precursor_peak_group_.updateMassesAndIntensity();
              return true;
            }
          }
        }
      }

      for (auto map = precursor_map_for_real_time_acquisition.upper_bound(scan_number_);
           map != precursor_map_for_real_time_acquisition.end();
           map++)
      {
        precursor_scan_number_ = map->first;
        //std::cout<<scan_number_ << " * " << precursor_scan_number_ <<std::endl;
        if (map != precursor_map_for_real_time_acquisition.end())
        {
          for (auto &smap : map->second)
          {
            //
            if (abs(start_mz - smap[3]) < .1 || abs(end_mz - smap[4]) < .1)
            {
              LogMzPeak precursor_log_mz_peak(precursor_peak_, is_positive);
              precursor_log_mz_peak.abs_charge = (int) smap[1];
              precursor_log_mz_peak.isotopeIndex = 0;
              precursor_log_mz_peak.mass = smap[0];
              precursor_peak_.setCharge(precursor_log_mz_peak.abs_charge);
              precursor_peak_group_.push_back(precursor_log_mz_peak);
              precursor_peak_group_.setQScore(smap[2]);
              precursor_peak_group_.setRepAbsCharge((int) smap[1]);
              precursor_peak_group_.updateMassesAndIntensity();
              return true;
            }
          }
        }
      }


      //std::cout<<scan_number_ << " " << precursor_peak_.getMZ() <<std::endl;
      return false;
    }

    //precursor_spectrum.updatePeakGroupMap();
    for (int i = survey_scans.size() - 1; i >= 0; i--)
    {
      auto precursor_spectrum = survey_scans[i];
      if (precursor_spectrum.empty())
      {
        continue;
      }
      for (auto &precursor: spec_.getPrecursors())
      {
        for (auto &activation_method :  precursor.getActivationMethods())
        {
          activation_method_ = Precursor::NamesOfActivationMethodShort[activation_method];
          break;
        }
        precursor_peak_ = precursor;
        precursor_scan_number_ = precursor_spectrum.scan_number_;
        double start_mz = precursor.getIsolationWindowLowerOffset() > 100.0 ?
                          precursor.getIsolationWindowLowerOffset() :
                          -precursor.getIsolationWindowLowerOffset() + precursor.getMZ();
        double end_mz = precursor.getIsolationWindowUpperOffset() > 100.0 ?
                        precursor.getIsolationWindowUpperOffset() :
                        precursor.getIsolationWindowUpperOffset() + precursor.getMZ();

        double max_isotope_cosine = -3.0;
        for (auto &pg: precursor_spectrum)
        {
          //std::sort(pg.begin(), pg.end());
          if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
          {
            continue;
          }

          //double sum_intensity = .0;
          double max_intensity = .0;
          const LogMzPeak *tmp_precursor = nullptr;
          for (auto &tmp_peak:pg)
          {
            if (tmp_peak.mz < start_mz)
            {
              continue;
            }
            if (tmp_peak.mz > end_mz)
            {
              break;
            }
            //sum_intensity += tmp_peak.intensity;

            if (tmp_peak.intensity < max_intensity)
            {
              continue;
            }
            max_intensity = tmp_peak.intensity;
            tmp_precursor = &tmp_peak;
          }

          if (pg.getIsotopeCosine() < max_isotope_cosine || tmp_precursor == nullptr)
          {
            continue;
          }

          //precursor_peak_.setMZ(tmp_precursor->mz);
          //precursor_peak_.setIntensity(tmp_precursor->intensity);
          precursor_peak_
              .setCharge(tmp_precursor->is_positive ? tmp_precursor->abs_charge : -tmp_precursor->abs_charge);
          max_isotope_cosine = pg.getIsotopeCosine();
          precursor_peak_group_ = pg;
        }
        if (!precursor_peak_group_.empty())
        {
          break;
        }
      }
      if (!precursor_peak_group_.empty())
      {
        break;
      }
    }
    return precursor_peak_group_.empty();
  }

  const MSSpectrum &DeconvolutedSpectrum::getOriginalSpectrum() const
  {
    return spec_;
  }

  PeakGroup DeconvolutedSpectrum::getPrecursorPeakGroup() const
  {
    if (precursor_peak_group_.empty())
    {
      return PeakGroup();
    }
    return precursor_peak_group_;
  }

  int DeconvolutedSpectrum::getPrecursorCharge() const
  {
    return precursor_peak_.getCharge();
  }

  double DeconvolutedSpectrum::getCurrentMaxMass(const double max_mass) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_mass;
    }
    return precursor_peak_group_.getMonoMass();
  }

  int DeconvolutedSpectrum::getCurrentMaxAbsCharge(const int max_abs_charge) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_abs_charge;
    }
    return abs(precursor_peak_.getCharge());
  }

  const Peak1D DeconvolutedSpectrum::getPrecursor() const
  {
    return precursor_peak_;
  }
}
