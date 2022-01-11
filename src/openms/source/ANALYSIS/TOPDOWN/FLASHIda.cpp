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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <sstream>

namespace OpenMS
{
  // constructor
  FLASHIda::FLASHIda(char *arg)
  {

    std::unordered_map<std::string, std::vector<double>> inputs;
    std::vector<String> log_files;
    std::vector<String> out_files;
    char *token = std::strtok(arg, " ");
    std::string key;
    while (token != nullptr)
    {
      String token_string = std::string(token);

      if (token_string.hasSuffix(".log"))
      {
        log_files.push_back(token_string);
        std::cout << token_string << " is used for global targeting\n";
      }
      else if (token_string.hasSuffix(".out"))
      {
        out_files.push_back(token_string);
        std::cout << token_string << " is used for global targeting\n";
      }
      else
      {
        double num = atof(token_string.c_str());

        if (num == 0 && !isdigit(token_string[token_string.size() - 1]))
        {
          key = token_string;
          inputs[key] = DoubleList();
        }
        else
        {
          inputs[key].push_back(num);
        }
      }
      token = std::strtok(nullptr, " ");
    }
    rt_window_ = inputs["RT_window"][0];
    qscore_threshold_ = inputs["score_threshold"][0];
    snr_threshold_ = 1;
    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", (int) inputs["min_charge"][0]);
    fd_defaults.setValue("max_charge", (int) inputs["max_charge"][0]);
    fd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    fd_defaults.setValue("max_mass", inputs["max_mass"][0]);
    fd_defaults.setValue("min_isotope_cosine", DoubleList{.9, .9});

    fd_defaults.setValue("min_qscore", .0);
    fd_defaults.setValue("tol", inputs["tol"]);
    tol_ = inputs["tol"];
    fd_defaults.setValue("rt_window", rt_window_);
    fd_defaults.setValue("min_peaks", IntList{4, 3}); // more sensitive

    auto mass_count_double = inputs["max_mass_count"];

    for (double j: mass_count_double)
    {
      mass_count_.push_back((int) j);
    }

    for (auto &log_file: log_files)
    {
      std::ifstream instream(log_file);
      if (instream.good())
      {
        String line;
        int scan;
        double mass;
        while (std::getline(instream, line))
        {
          if (line.find("0 targets") != line.npos)
          {
            continue;
          }
          if (line.hasPrefix("MS1"))
          {
            Size st = line.find("Aceess ID ") + 10;
            Size ed = line.find(')');
            String n = line.substr(st, ed);
            scan = atoi(n.c_str());
          }
          if (line.hasPrefix("Mass"))
          {
            Size st = 5;
            Size ed = line.find('\t');
            String n = line.substr(st, ed);
            mass = atof(n.c_str());
            target_masses_.insert(mass);
            int nmass = FLASHDeconvAlgorithm::getNominalMass(mass);
            if (target_nominal_masses_.find(nmass) == target_nominal_masses_.end())
            {
              target_nominal_masses_[nmass] = std::vector<double>();
            }
            target_nominal_masses_[nmass].push_back((double) (scan * 5400.0 / 25000.0));
            //precursor_map_for_real_time_acquisition[scan].push_back(e);
          }
        }
        instream.close();
      }
    }

    for (auto &log_file: out_files)
    {
      std::ifstream instream(log_file);

      if (instream.good())
      {
        String line;
        double mass;
        while (std::getline(instream, line))
        {
          if (line.hasPrefix("rt"))
          {
            continue;
          }
          std::stringstream tmp_stream(line);
          String str;
          std::vector<String> results;
          while (getline(tmp_stream, str, '\t'))
          {
            results.push_back(str);
          }
          mass = atof(results[5].c_str());
          target_masses_.insert(mass);
          int nmass = FLASHDeconvAlgorithm::getNominalMass(mass);
          if (target_nominal_masses_.find(nmass) == target_nominal_masses_.end())
          {
            target_nominal_masses_[nmass] = std::vector<double>();
          }
          target_nominal_masses_[nmass].push_back(60.0 * atof(results[0].c_str()));
        }
        instream.close();
      }
    }

    fd_.setParameters(fd_defaults);
    fd_.calculateAveragine(false);
    fd_.setTargetMasses(target_masses_, 1);

    std::cout << "QScore threshold: " << qscore_threshold_ << std::endl;
    std::cout << fd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(const double *mzs,
                              const double *ints,
                              const int length,
                              const double rt,
                              const int ms_level,
                              const char *name)
  {

    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);
    //deconvoluted_spectrum_ = DeconvolutedSpectrum(spec, 1);
    if (ms_level == 1)
    {
      //current_max_mass_ = max_mass;
      //currentChargeRange = chargeRange;
    }
    else
    {
      return 0;
      //TODO precursor infor here
    }

    std::vector<DeconvolutedSpectrum> tmp;
    std::map<int, std::vector<std::vector<double>>> empty;
    deconvoluted_spectrum_ = fd_.getDeconvolutedSpectrum(spec, tmp, 0, empty);

    // per spec deconvolution
    FLASHIda::filterPeakGroupsUsingMassExclusion_(ms_level, rt);

    // spec.clear(true);

    return deconvoluted_spectrum_.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    std::vector<PeakGroup> filtered_peakgroups;

    std::sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end(), QscoreComparator_);
    int mass_count = mass_count_[ms_level - 1];
    trigger_charges.clear();
    trigger_charges.reserve(mass_count);

    filtered_peakgroups.reserve(mass_count_.size());
    std::set<int> current_selected_mzs; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_qscore_map_;

    for (auto&[m, r]: tqscore_exceeding_mz_rt_map_)
    {
      if (rt - r > rt_window_)
      {
        continue;
      }
      new_mz_rt_map_[m] = r;
    }
    new_mz_rt_map_.swap(tqscore_exceeding_mz_rt_map_);
    std::unordered_map<int, double>().swap(new_mz_rt_map_);

    for (auto&[m, r]: tqscore_exceeding_mass_rt_map_)
    {
      if (rt - r > rt_window_)
      {
        continue;
      }
      new_mass_rt_map_[m] = r;
    }
    new_mass_rt_map_.swap(tqscore_exceeding_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_mass_rt_map_);

    for (auto &item: all_mass_rt_map_)
    {
      if (rt - item.second > rt_window_)
      {
        continue;
      }
      new_all_mass_rt_map_[item.first] = item.second;
      new_mass_qscore_map_[item.first] = mass_qscore_map_[item.first];
    }
    new_all_mass_rt_map_.swap(all_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_all_mass_rt_map_);

    new_mass_qscore_map_.swap(mass_qscore_map_);
    std::unordered_map<int, double>().swap(new_mass_qscore_map_);

    //When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold.
    //when selection_phase == 1, consider all other masses for selection
    for (int selection_phase = 0; selection_phase < 2; selection_phase++)
    {
      for (auto &pg: deconvoluted_spectrum_)
      {

        if (filtered_peakgroups.size() >= mass_count)
        {
          break;
        }

        if (pg.getQScore() < qscore_threshold_)
        {
          break;
        }

        if (pg.getChargeSNR(pg.getRepAbsCharge()) < snr_threshold_)
        {
          continue;
        }

        int charge = pg.getRepAbsCharge();
        double qscore = pg.getQScore();

        int mz = (int) round(
            (std::get<0>(pg.getMaxQScoreMzRange()) + std::get<1>(pg.getMaxQScoreMzRange())) / 2.0);

        int nominal_mass = FLASHDeconvAlgorithm::getNominalMass(pg.getMonoMass());

        if (current_selected_mzs.find(mz) != current_selected_mzs.end())
        {
          continue;
        }

        if (selection_phase == 0)
        { // first, select masses or m/zs outside exclusion list
          if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end()
              ||
              tqscore_exceeding_mz_rt_map_.find(mz) != tqscore_exceeding_mz_rt_map_.end()
              )
          {
            continue;
          }
        }

        // here crawling isolation windows max_isolation_window_half_
        auto ospec = deconvoluted_spectrum_.getOriginalSpectrum();
        Size index = ospec.findNearest(mz);
        int tindexl = index;
        int tindexr = index + 1;
        double lmz = -1.0, rmz = -1.0;

        double sig_pwr = .0;
        double noise_pwr = .0;
        bool lkeep = true;
        bool rkeep = false;
        auto pmz_range = pg.getMzRange(charge);
        double pmzl = std::get<0>(pmz_range);
        double pmzr = std::get<1>(pmz_range);

        while (lkeep || rkeep)
        {
          if (tindexl < 0 || ospec[tindexl].getMZ() < mz - max_isolation_window_half_)// left side done
          {
            lkeep = false;
          }
          else
          {
            if (pg.isSignalMZ(ospec[tindexl].getMZ(), tol_[ospec.getMSLevel() - 1])) //
            {
              sig_pwr += pg.getIntensity() * pg.getIntensity();
            }
            else
            {
              noise_pwr += pg.getIntensity() * pg.getIntensity();
            }
          }
          if (tindexr >= ospec.size() || ospec[tindexr].getMZ() > mz + max_isolation_window_half_)// right side done
          {
            rkeep = false;
          }
          else
          {
            if (pg.isSignalMZ(ospec[tindexr].getMZ(), tol_[ospec.getMSLevel() - 1])) //
            {
              sig_pwr += pg.getIntensity() * pg.getIntensity();
            }
            else
            {
              noise_pwr += pg.getIntensity() * pg.getIntensity();
            }
          }

          tindexl--;
          tindexr++;

          if (ospec[tindexr].getMZ() < pmzr + min_isolation_window_half_)
          {
            continue;
          }

          if (ospec[tindexl].getMZ() > pmzl - min_isolation_window_half_)
          {
            continue;
          }

          if (sig_pwr / noise_pwr < snr_threshold_)
          {
            break;
          }
          else
          {
            lmz = ospec[tindexl + 1].getMZ();
            rmz = ospec[tindexr - 1].getMZ();
          }
        }

        if (lmz < 0 || rmz < 0)
        {
          continue;
        }

        all_mass_rt_map_[nominal_mass] = rt;
        auto inter = mass_qscore_map_.find(nominal_mass);
        if (inter == mass_qscore_map_.end())
        {
          mass_qscore_map_[nominal_mass] = 1 - qscore;
        }
        else
        {
          mass_qscore_map_[nominal_mass] *= 1 - qscore;
        }

        if (1 - mass_qscore_map_[nominal_mass] > tqscore_threshold)
        {
          tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
          tqscore_exceeding_mz_rt_map_[mz] = rt;
        }
        filtered_peakgroups.push_back(pg);
        trigger_charges.push_back(charge);
        trigger_left_isolation_mzs_.push_back(lmz);
        trigger_right_isolation_mzs_.push_back(rmz);

        current_selected_mzs.insert(mz);
      }
    }

    deconvoluted_spectrum_.swap(filtered_peakgroups);
    std::vector<PeakGroup>().swap(filtered_peakgroups);
  }

  void FLASHIda::getIsolationWindows(double *wstart,
                                     double *wend,
                                     double *qscores,
                                     int *charges,
                                     int *min_charges,
                                     int *max_charges,
                                     double *mono_masses,
                                     double *chare_cos,
                                     double *charge_snrs,
                                     double *iso_cos,
                                     double *snrs, double *charge_scores,
                                     double *ppm_errors,
                                     double *precursor_intensities,
                                     double *peakgroup_intensities)
  {
    //std::sort(deconvoluted_spectrum_.begin(), deconvoluted_spectrum_.end(), QscoreComparator_);

    for (int i = 0; i < deconvoluted_spectrum_.size(); i++)
    {
      if (trigger_charges[i] == 0)
      {
        continue;
      }
      auto peakgroup = deconvoluted_spectrum_[i];
      charges[i] = trigger_charges[i];
      auto cr = peakgroup.getAbsChargeRange();
      min_charges[i] = std::get<0>(cr);
      max_charges[i] = std::get<1>(cr);

      wstart[i] = trigger_left_isolation_mzs_[i];  //std::get<0>(mz_range) - min_isolation_window_half_;
      wend[i] = trigger_right_isolation_mzs_[i]; //std::get<1>(mz_range) + min_isolation_window_half_;

      qscores[i] = QScore::getQScore(&peakgroup, charges[i]);
      mono_masses[i] = peakgroup.getMonoMass();
      chare_cos[i] = peakgroup.getChargeIsotopeCosine(charges[i]);
      charge_snrs[i] = peakgroup.getChargeSNR(charges[i]);
      iso_cos[i] = peakgroup.getIsotopeCosine();
      snrs[i] = peakgroup.getSNR();
      charge_scores[i] = peakgroup.getChargeScore();
      ppm_errors[i] = peakgroup.getAvgPPMError();
      peakgroup_intensities[i] = peakgroup.getIntensity();
      precursor_intensities[i] = peakgroup.getChargeIntensity(charges[i]);
    }
    std::vector<PeakGroup> empty;
    deconvoluted_spectrum_.swap(empty);
  }

  MSSpectrum FLASHIda::makeMSSpectrum_(const double *mzs, const double *ints, const int length, const double rt,
                                       const int ms_level, const char *name)
  {
    auto spec = MSSpectrum();
    for (int i = 0; i < length; i++)
    {
      if (ints[i] <= 0)
      {
        continue;
      }
      spec.emplace_back(mzs[i], ints[i]);
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);
    return spec;
  }
}
