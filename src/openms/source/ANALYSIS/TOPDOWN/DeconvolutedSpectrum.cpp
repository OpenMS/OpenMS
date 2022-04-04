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

#include "include/OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h"

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
    for (auto &pg: *this)
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
      precursor.setCharge(mass_charge);
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
                                                     const bool write_detail)
  {
    if (empty())
    {
      return;
    }

    for (auto &pg: *this)
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
        for (auto &p: pg)
        {
          fs << p.mz << " ";
        }
        fs << ";\t";

        fs << std::fixed << std::setprecision(1);
        for (auto &p: pg)
        {
          fs << p.intensity << " ";
        }
        fs << ";\t";
        fs << std::setprecision(-1);

        for (auto &p: pg)
        {
          fs << (p.is_positive ? p.abs_charge : -p.abs_charge) << " ";
        }
        fs << ";\t";
        for (auto &p: pg)
        {
          fs << p.getUnchargedMass() << " ";
        }
        fs << ";\t";
        for (auto &p: pg)
        {
          fs << p.isotopeIndex << " ";
        }
        fs << ";\t";

        for (auto &p: pg)
        {
          double average_mass = pg.getMonoMass() + p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U;
          double mass_error =
              (average_mass / p.abs_charge + FLASHDeconvHelperStructs::getChargeMass(p.is_positive) -
               p.mz) / p.mz;
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

      for (auto &p: pg)
      {
        isotope_end_index = isotope_end_index < p.isotopeIndex ? p.isotopeIndex : isotope_end_index;
      }
      auto per_isotope_intensity = std::vector<double>(isotope_end_index + 1, .0);
      for (auto &p: pg)
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

  void DeconvolutedSpectrum::writeTopFD(std::fstream &fs,
                                        const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                                        const double snr_threshold,
                                        const double harmonic_factor,
                                        const double precursor_offset)//, fstream& fsm, fstream& fsp)
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
           << "PRECURSOR_CHARGE=" << (int) (precursor_peak_.getCharge() * harmonic_factor) << "\n"
           << "PRECURSOR_MASS="
           << std::to_string(precursor_peak_group_.getMonoMass() * harmonic_factor + precursor_offset) << "\n"
           << "PRECURSOR_INTENSITY=" << precursor_peak_.getIntensity() << "\n";
      }
      else
      {
        double average_mass =
            (precursor_peak_.getMZ() -
             FLASHDeconvHelperStructs::getChargeMass(precursor_peak_.getCharge() > 0)) *
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

    if (size() > topFD_max_peak_count_)// max peak count for TopPic = 500
    {
      isotope_scores.reserve(size());
      for (auto &pg: *this)
      {
        isotope_scores.push_back(pg.getIsotopeCosine());
      }
      std::sort(isotope_scores.begin(), isotope_scores.end());
      isotope_score_threshold = isotope_scores[isotope_scores.size() - topFD_max_peak_count_];
      std::vector<double>().swap(isotope_scores);
    }

    int size = 0;
    for (auto &pg: *this)
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
      for (auto &peaks : pg)
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


  bool DeconvolutedSpectrum::registerPrecursor(const std::vector<DeconvolutedSpectrum> &survey_scans,
                                               const bool is_positive, double isolation_window_size_,
                                               const std::map<int, std::vector<std::vector<double>>> &precursor_map_for_real_time_acquisition)
  {
    precursor_peak_.setIntensity(.0);
    double start_mz = 0;
    double end_mz = 0;
    //int target_precursor_scan = -1;
    for (auto &precursor: spec_.getPrecursors())
    {
      for (auto &activation_method: precursor.getActivationMethods())
      {
        activation_method_ = Precursor::NamesOfActivationMethodShort[activation_method];
        if (!activation_method_.compare("HCID"))
        {
          activation_method_ = "HCD";
        }
        break;
      }
      precursor_peak_ = precursor;

      double loffset = precursor.getIsolationWindowLowerOffset();
      double uoffest = precursor.getIsolationWindowUpperOffset();
      loffset = loffset <=0? isolation_window_size_/2.0 : loffset;
      uoffest = uoffest <=0? isolation_window_size_/2.0 : uoffest;


      start_mz = loffset > 100.0 ?
                     loffset :
                 -loffset + precursor.getMZ();
      end_mz = uoffest > 100.0 ?
                   uoffest:
                   uoffest + precursor.getMZ();
      //std::cout<<  precursor.getIsolationWindowLowerOffset()  << " " << precursor.getIsolationWindowUpperOffset() << std::endl;
    }
    if (!precursor_map_for_real_time_acquisition.empty() && precursor_peak_group_.empty())
    {
      auto map = precursor_map_for_real_time_acquisition.lower_bound(scan_number_);
      while (map != precursor_map_for_real_time_acquisition.begin())
      {
        --map;
        if (map->first < scan_number_ - 50)
        {
          return false;
        }

        if (map != precursor_map_for_real_time_acquisition.end())
        {
          for (auto &smap: map->second)
          {
            if (abs(start_mz - smap[3]) < .001 && abs(end_mz - smap[4]) < .001)
            {
              LogMzPeak precursor_log_mz_peak(precursor_peak_, is_positive);
              precursor_log_mz_peak.abs_charge = (int) smap[1];
              precursor_log_mz_peak.isotopeIndex = 0;
              precursor_log_mz_peak.mass = smap[0];
              precursor_log_mz_peak.intensity = smap[6];
              //precursor_peak_.setMetaValue("color", smap[7]);
              precursor_peak_group_.push_back(precursor_log_mz_peak);
              precursor_peak_.setCharge(precursor_log_mz_peak.abs_charge);
              precursor_peak_.setIntensity(smap[5]);
              precursor_peak_group_.setAbsChargeRange(smap[7], smap[8]);
              precursor_peak_group_.setChargeIsotopeCosine(precursor_log_mz_peak.abs_charge, smap[9]);
              precursor_peak_group_.setChargeSNR(precursor_log_mz_peak.abs_charge, smap[10]);//cnsr
              precursor_peak_group_.setIsotopeCosine(smap[11]);
              precursor_peak_group_.setSNR(smap[12]);
              precursor_peak_group_.setChargeScore(smap[13]);
              precursor_peak_group_.setAvgPPMError(smap[14]);
              precursor_peak_group_.setQScore(smap[2]);
              precursor_peak_group_.setRepAbsCharge((int) smap[1]);
              precursor_peak_group_.updateMassesAndIntensity();
              //precursor_peak_group_.setScanNumber()
              precursor_scan_number_ = map->first;
              //std::cout<<precursor_scan_number_<<" " << precursor_peak_group_.getMonoMass()<<std::endl;

              for (int i = survey_scans.size() - 1; i >= 0; i--)
              {
                auto precursor_spectrum = survey_scans[i];
                if (precursor_spectrum.getScanNumber() != precursor_scan_number_)
                {
                  continue;
                }

                auto o_spec = precursor_spectrum.getOriginalSpectrum();
                auto spec_iterator = o_spec.PosBegin(start_mz);
                while (spec_iterator != o_spec.end() && spec_iterator->getMZ() < end_mz)
                {
                  double intensity = spec_iterator->getIntensity();
                  spec_iterator++;
                }

                double max_score = -1.0;
                for (auto &pg: precursor_spectrum)
                {
                  if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
                  {
                    continue;
                  }
                  if (abs(pg.getMonoMass() - smap[0]) > 5.0)
                  {
                    continue;
                  }
                  double max_intensity = .0;
                  const LogMzPeak *tmp_precursor = nullptr;
                  // double power = .0;
                  for (auto &tmp_peak: pg)
                  {
                    if (tmp_peak.mz < start_mz)
                    {
                      continue;
                    }
                    if (tmp_peak.mz > end_mz)
                    {
                      break;
                    }
                    //power += tmp_peak.intensity * tmp_peak.intensity;

                    if (tmp_peak.intensity < max_intensity)
                    {
                      continue;
                    }
                    max_intensity = tmp_peak.intensity;
                    tmp_precursor = &tmp_peak;
                    //sum_int += tmp_peak.intensity * tmp_peak.intensity;
                  }

                  if (tmp_precursor == nullptr)
                  {
                    continue;
                  }
                  //double t_cos = pg.getIsotopeCosine();
                  //auto precursor_snr = pg.getChargeSNR(tmp_precursor->abs_charge);//t_cos * t_cos * power
                  //                    // / (total_power - power + (1 - t_cos) * (1 - t_cos) * power);

                  auto score = pg.getChargeSNR(tmp_precursor->abs_charge);
                  // // most intense one should determine the mass
                  if (score < max_score)
                  {
                    continue;
                  }

                  max_score = score;
                  precursor_peak_group_ = pg;
                }
              }

              return true;
            }
          }
        }
      }

      return false;
    }

    double max_score = -1.0;

    for (int i = survey_scans.size() - 1; i >= 0; i--)
    {
      auto precursor_spectrum = survey_scans[i];

      if (precursor_spectrum.empty())
      {
        continue;
      }
      auto o_spec = precursor_spectrum.getOriginalSpectrum();
      auto spec_iterator = o_spec.PosBegin(start_mz);

      // double total_power = .0;
      while (spec_iterator != o_spec.end() && spec_iterator->getMZ() < end_mz)
      {
        //total_power += spec_iterator->getIntensity() * spec_iterator->getIntensity();
        spec_iterator++;
      }
      for (auto &pg: precursor_spectrum)
      {
        //std::sort(pg.begin(),pg.end());
        if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz)
        {
          continue;
        }

        //double sum_int = 0;
        double max_intensity = .0;
        const LogMzPeak *tmp_precursor = nullptr;

        int c = int(.5 + pg.getMonoMass() / start_mz);
        //bool contained = true;
        //double power = 0;
        //double i_sum = 0;
        for (auto &tmp_peak: pg)
        {
          if (tmp_peak.abs_charge != c)
          {
            continue;
          }
          //sum += tmp_peak.mz * tmp_peak.intensity;

          if (tmp_peak.mz < start_mz || tmp_peak.mz > end_mz)
          {
            continue;
          }
          // power += tmp_peak.intensity * tmp_peak.intensity;

          if (tmp_peak.intensity < max_intensity)
          {
            continue;
          }
          max_intensity = tmp_peak.intensity;
          tmp_precursor = &tmp_peak;
        }

        if (tmp_precursor == nullptr)// || i_sum <= 0 || sum < start_mz || sum > end_mz)
        {
          continue;
        }
        //cos2θ ||V||2/(N1+sin2θ ||V||2).
        //auto t_cos = pg.getIsotopeCosine();
        //auto precursor_snr = t_cos * t_cos * power
        //                     / (total_power - power + (1 - t_cos) * (1 - t_cos) * power);

        auto score = pg.getChargeSNR(tmp_precursor->abs_charge); // most intense one should determine the mass
        if (score < max_score)
        {
          continue;
        }
        precursor_peak_
            .setCharge(tmp_precursor->is_positive ? tmp_precursor->abs_charge
                                                  : -tmp_precursor->abs_charge);

        precursor_peak_.setIntensity(tmp_precursor->intensity);
        max_score = score;
        //pg.setChargeSNR(tmp_precursor->abs_charge, precursor_snr);
        precursor_peak_group_ = pg;
        precursor_scan_number_ = precursor_spectrum.scan_number_;
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

  double DeconvolutedSpectrum::getCurrentMinMass(const double min_mass) const
  {
    if (spec_.getMSLevel() == 1)
    {
      return min_mass;
    }
    return 50.0;
  }

  int DeconvolutedSpectrum::getCurrentMaxAbsCharge(const int max_abs_charge) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_abs_charge;
    }
    return abs(precursor_peak_.getCharge());
  }

  const Precursor DeconvolutedSpectrum::getPrecursor() const
  {
    return precursor_peak_;
  }

  int DeconvolutedSpectrum::getScanNumber() const
  {
    return scan_number_;
  }

  int DeconvolutedSpectrum::getPrecursorScanNumber() const
  {
    return precursor_scan_number_;
  }
}
