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
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
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

        if (num == 0&&  !isdigit(token_string[token_string.size() - 1]))
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
    fd_defaults.setValue("min_isotope_cosine", DoubleList{.85, .85});

    fd_defaults.setValue("min_qscore", .0);
    fd_defaults.setValue("tol", inputs["tol"]);
    tol_ = inputs["tol"];
    fd_defaults.setValue("rt_window", rt_window_);
    fd_defaults.setValue("min_peaks", IntList{3, 3});//

    auto mass_count_double = inputs["max_mass_count"];

    for (double j: mass_count_double)
    {
      mass_count_.push_back((int) j);
    }

    for (auto& log_file: log_files)
    {
      std::ifstream instream(log_file);
      if (instream.good())
      {
        String line;
        double rt = .0;
        double mass;
        while (std::getline(instream, line))
        {
          if (line.find("0 targets") != line.npos)
          {
            continue;
          }
          if (line.hasPrefix("MS1"))
          {
            Size st = line.find("RT ") + 3;
            Size ed = line.find('(') - 2;
            String n = line.substr(st, ed - st + 1);
            rt = atof(n.c_str());
            //std::cout << " rt " << n << std::endl ;
            //precursor_map_for_real_time_acquisition[scan] = std::vector<std::vector<double>>();//// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
          }
          if (line.hasPrefix("Mass"))
          {
            Size st = 5;
            Size ed = line.find('\t');
            String n = line.substr(st, ed - st + 1);
            mass = atof(n.c_str());

            if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
            {
              target_mass_rt_map_[mass] = std::vector<double>();
            }
            target_mass_rt_map_[mass].push_back(rt * 60.0);
            //precursor_map_for_real_time_acquisition[scan].push_back(e);
          }
        }
        instream.close();
      }
    }

    for (auto& log_file: out_files)
    {
      std::ifstream instream(log_file);

      if (instream.good())
      {
        String line;
        double mass;
        double mz;
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
          mz = (atof(results[1].c_str()) + atof(results[2].c_str())) / 2.0;
          if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
          {
            target_mass_rt_map_[mass] = std::vector<double>();
          }
          target_mass_rt_map_[mass].push_back(60.0 * atof(results[0].c_str()));

          if (target_mz_rt_map_.find(mz) == target_mz_rt_map_.end())
          {
            target_mz_rt_map_[mz] = std::vector<double>();
          }
          target_mz_rt_map_[mz].push_back(60.0 * atof(results[0].c_str()));
        }
        instream.close();
      }
    }

    fd_.setParameters(fd_defaults);
    fd_.calculateAveragine(false);

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
    //deconvolved_spectrum_ = DeconvolvedSpectrum(spec, 1);
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

    std::vector<DeconvolvedSpectrum> tmp;
    std::map<int, std::vector<std::vector<double>>> empty;

    target_masses_.clear();
    for (auto&[mass, rts]: target_mass_rt_map_)
    {
      for (double prt: rts)
      {
        if (std::abs(rt - prt) < 300)
        {
          target_masses_.push_back(mass);
          break;
        }
      }
    }
    std::sort(target_masses_.begin(), target_masses_.end());
    fd_.setTargetMasses(target_masses_);
    fd_.performSpectrumDeconvolution(spec, tmp, 0, false, empty);
    deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();
    // per spec deconvolution
    FLASHIda::filterPeakGroupsUsingMassExclusion_(ms_level, rt);
    // spec.clear(true);

    return deconvolved_spectrum_.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    std::vector<PeakGroup> filtered_peakgroups;
    deconvolved_spectrum_.sortByQScore();
    Size mass_count = (Size) mass_count_[ms_level - 1];
    trigger_charges.clear();
    trigger_charges.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);

    filtered_peakgroups.reserve(mass_count_.size());
    std::set<int> current_selected_mzs;// current selected mzs

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

    for (auto& item: all_mass_rt_map_)
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

    int selection_phase_start = 0;
    int selection_phase_end = target_masses_.size() > 0 ? 1 : 2;
    //When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold.
    //when selection_phase == 1, consider all other masses for selection
    for (int selection_phase = selection_phase_start; selection_phase < selection_phase_end; selection_phase++)
    {
      for (auto& pg: deconvolved_spectrum_)
      {

        if (filtered_peakgroups.size() >= mass_count)
        {
          break;
        }

        int charge = pg.getRepAbsCharge();
        double qscore = pg.getQScore();
        double mass = pg.getMonoMass();
        double center_mz =
            (std::get<0>(pg.getMaxQScoreMzRange()) + std::get<1>(pg.getMaxQScoreMzRange())) / 2.0;

        int nominal_mass = FLASHDeconvAlgorithm::getNominalMass(mass);
        bool include = false;
        double snr_threshold = snr_threshold_;

        if (target_masses_.size() > 0)
        {
          double delta = 2 * tol_[0] * mass * 1e-6;
          auto ub = std::upper_bound(target_masses_.begin(), target_masses_.end(), mass + delta);

          while (ub != target_masses_.begin()&&  !include&&  *ub > mass - delta)
          {
            --ub;
            if (std::abs(*ub - mass) < delta)
            {
              include = true;
              snr_threshold = 0.0;
              break;
            }
          }


          if (!include)
          {
            continue;
          }
        }

        if (!include&&  qscore < qscore_threshold_)
        {
          break;
        }

        if (pg.getChargeSNR(charge) < snr_threshold)
        {
          continue;
        }


        int integer_mz = (int) round(center_mz);

        if (current_selected_mzs.find(integer_mz) != current_selected_mzs.end())
        {
          continue;
        }

        if (selection_phase == 0)
        {// first, select masses or m/zs outside exclusion list
          if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end() ||
              tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end())
          {
            continue;
          }
        }

        // here crawling isolation windows max_isolation_window_half_
        auto ospec = deconvolved_spectrum_.getOriginalSpectrum();
        Size index = ospec.findNearest(center_mz);
        int tindexl = index;
        Size tindexr = index + 1;
        double lmz = -1.0, rmz = -1.0;

        double sig_pwr = .0;
        double noise_pwr = .0;
        bool lkeep = true;
        bool rkeep = true;
        auto pmz_range = pg.getMzRange(charge);
        double pmzl = std::get<0>(pmz_range);
        double pmzr = std::get<1>(pmz_range);
        //double final_snr = .0;

        while (lkeep || rkeep)
        {
          if (tindexl < 0 || ospec[tindexl].getMZ() < center_mz - max_isolation_window_half_)// left side done
          {
            lkeep = false;
          }
          else
          {
            double power = pow(ospec[tindexl].getIntensity(), 2);
            if (pg.isSignalMZ(ospec[tindexl].getMZ(), tol_[ospec.getMSLevel() - 1]))//
            {
              sig_pwr += power;
            }
            else
            {
              noise_pwr += power;
            }
          }
          if (tindexr >= ospec.size() || ospec[tindexr].getMZ() > center_mz + max_isolation_window_half_)// right side done
          {
            rkeep = false;
          }
          else
          {
            double power = pow(ospec[tindexr].getIntensity(), 2);
            if (pg.isSignalMZ(ospec[tindexr].getMZ(), tol_[ospec.getMSLevel() - 1]))//
            {
              sig_pwr += power;
            }
            else
            {
              noise_pwr += power;
            }
          }

          double tlmz = ospec[tindexl].getMZ();
          double trmz = ospec[tindexr].getMZ();

          tindexl--;
          tindexr++;

          if (trmz < pmzr + min_isolation_window_half_)
          {
            continue;
          }

          if (tlmz > pmzl - min_isolation_window_half_)
          {
            continue;
          }

          if (sig_pwr / noise_pwr < snr_threshold)
          {
            break;
          }
          else
          {
            if (tindexl >=0 && tindexl < (int) ospec.size() - 1)
            {
              lmz = std::max(center_mz - max_isolation_window_half_, tlmz);
            }
            if (tindexr > 0 && tindexr < ospec.size())
            {
              rmz = std::min(center_mz + max_isolation_window_half_, trmz);
            }
            //final_snr = sig_pwr / noise_pwr;
          }
        }


        if (lmz < ospec[0].getMZ() - max_isolation_window_half_ || rmz > ospec.back().getMZ() + max_isolation_window_half_ || lmz + 2 * min_isolation_window_half_ > rmz || rmz - lmz > 2 * max_isolation_window_half_)
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
          tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
        }
        filtered_peakgroups.push_back(pg);
        trigger_charges.push_back(charge);

        trigger_left_isolation_mzs_.push_back(lmz);
        trigger_right_isolation_mzs_.push_back(rmz);

        current_selected_mzs.insert(integer_mz);
      }
    }

    deconvolved_spectrum_.swap(filtered_peakgroups);
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
    //std::sort(deconvolved_spectrum_.begin(), deconvolved_spectrum_.end(), QscoreComparator_);

    for (Size i = 0; i < deconvolved_spectrum_.size(); i++)
    {
      if (trigger_charges[i] == 0)
      {
        continue;
      }
      auto peakgroup = deconvolved_spectrum_[i];
      charges[i] = trigger_charges[i];
      auto cr = peakgroup.getAbsChargeRange();
      min_charges[i] = std::get<0>(cr);
      max_charges[i] = std::get<1>(cr);

      wstart[i] = trigger_left_isolation_mzs_[i];//std::get<0>(mz_range) - min_isolation_window_half_;
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
    deconvolved_spectrum_.swap(empty);
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
}// namespace OpenMS
