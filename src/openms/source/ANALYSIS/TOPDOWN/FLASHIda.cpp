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
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  /// minimum isolation window width divided by two
  inline const double min_isolation_window_half_ = .6;
  /// maximum isolation window width divided by two
  inline const double max_isolation_window_half_ = 3.0;
  /// constructor
  FLASHIda::FLASHIda(char* arg)
  {
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    std::unordered_map<std::string, std::vector<double>> inputs;
    std::vector<String> log_files;
    std::vector<String> out_files;
    char* token = std::strtok(arg, " ");
    std::string key;
    std::stringstream ss {};

    while (token != nullptr)
    {
      String token_string = std::string(token);

      if (token_string.hasSuffix(".log"))
      {
        log_files.push_back(token_string);
        ss << token_string << " ";
      }
      else if (token_string.hasSuffix(".out"))
      {
        out_files.push_back(token_string);
        ss << token_string << " ";
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
    targeting_mode_ = (int)(inputs["target_mode"][0]);
    if (targeting_mode_ == 1)
    {
      std::cout << ss.str() << "file(s) is(are) used to generate inclusion list\n";
    }
    else if (targeting_mode_ == 2)
    {
      std::cout << ss.str() << "file(s) is(are) used to generate exclusion list\n";
    }
    Param fd_defaults = FLASHDeconvAlgorithm().getDefaults();

    fd_defaults.setValue("min_charge", (int)inputs["min_charge"][0]);
    fd_defaults.setValue("max_charge", (int)inputs["max_charge"][0]);
    fd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    fd_defaults.setValue("max_mass", inputs["max_mass"][0]);
    fd_defaults.setValue("min_isotope_cosine", DoubleList {.85, .85});
    // fd_defaults.setValue("min_qscore", .0);
    fd_defaults.setValue("tol", inputs["tol"]);
    tol_ = std::vector<double>(inputs["tol"]);
    // fd_defaults.setValue("rt_window", rt_window_);
    // fd_defaults.setValue("min_peaks", IntList{3, 3});//

    auto mass_count_double = inputs["max_mass_count"];

    for (double j : mass_count_double)
    {
      mass_count_.push_back((int)j);
    }

    for (auto& log_file : log_files)
    {
      std::ifstream instream(log_file);
      if (instream.good())
      {
        String line;
        double rt = .0;
        double mass;
        double qscore;
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
            // precursor_map_for_real_time_acquisition[scan] = std::vector<std::vector<double>>();//// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
          }
          if (line.hasPrefix("Mass"))
          {
            Size st = 5;
            Size ed = line.find('\t');
            String n = line.substr(st, ed - st + 1);
            mass = atof(n.c_str());

            st = line.find("Score=") + 6;
            ed = line.find('\t', st);
            n = line.substr(st, ed - st + 1);
            qscore = atof(n.c_str());

            if (targeting_mode_ > 0)
            {
              if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
              {
                target_mass_rt_map_[mass] = std::vector<double>();
              }
              target_mass_rt_map_[mass].push_back(rt * 60.0);
              if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end())
              {
                target_mass_qscore_map_[mass] = std::vector<double>();
              }
              target_mass_qscore_map_[mass].push_back(qscore);
            }
            // precursor_map_for_real_time_acquisition[scan].push_back(e);
          }
        }
        instream.close();
      }
    }

    for (auto& log_file : out_files)
    {
      std::ifstream instream(log_file);

      if (instream.good())
      {
        String line;
        double mass;
        double qscore;

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
          qscore = atof(results[3].c_str());

          if (targeting_mode_ > 0)
          {
            if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
            {
              target_mass_rt_map_[mass] = std::vector<double>();
            }
            target_mass_rt_map_[mass].push_back(60.0 * atof(results[0].c_str()));
            if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end())
            {
              target_mass_qscore_map_[mass] = std::vector<double>();
            }
            target_mass_qscore_map_[mass].push_back(qscore);
          }
        }
        instream.close();
      }
    }

    fd_.setParameters(fd_defaults);
    fd_.calculateAveragine(false);

    std::cout << "QScore threshold: " << qscore_threshold_ << std::endl;
    std::cout << fd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(const double* mzs, const double* ints, const int length, const double rt, const int ms_level, const char* name)
  {
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);
    // deconvolved_spectrum_ = DeconvolvedSpectrum(spec, 1);
    if (ms_level == 1)
    {
      // current_max_mass_ = max_mass;
      // currentChargeRange = chargeRange;
    }
    else
    {
      return 0;
      // TODO precursor infor here
    }

    std::vector<DeconvolvedSpectrum> tmp;
    std::map<int, std::vector<std::vector<float>>> empty;

    target_masses_.clear();
    if (targeting_mode_ > 0)
    {
      for (auto& [mass, rts] : target_mass_rt_map_)
      {
        for (double prt : rts)
        {
          if (std::abs(rt - prt) < rt_window_)
          {
            target_masses_.push_back(mass);
            break;
          }
        }
      }
      std::sort(target_masses_.begin(), target_masses_.end());
      fd_.setTargetMasses(target_masses_, false);
    }

    fd_.performSpectrumDeconvolution(spec, tmp, 0, empty);
    deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();
    // per spec deconvolution
    FLASHIda::filterPeakGroupsUsingMassExclusion_(ms_level, rt);
    // spec.clear(true);

    return (int)deconvolved_spectrum_.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    std::vector<PeakGroup> filtered_peakgroups;
    deconvolved_spectrum_.sortByQScore();
    Size mass_count = (Size)mass_count_[ms_level - 1];
    trigger_charges.clear();
    trigger_charges.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);

    filtered_peakgroups.reserve(mass_count_.size());
    std::set<int> current_selected_mzs; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_qscore_map_;
    std::unordered_map<int, double> t_mass_qscore_map_;

    if (targeting_mode_ == 2)
    {
      for (auto& [mass, rts] : target_mass_rt_map_)
      {
        int nominal_mass = FLASHDeconvAlgorithm::getNominalMass(mass);
        auto qscores = target_mass_qscore_map_[mass];
        for (int i = 0; i < rts.size(); i++)
        {
          double prt = rts[i];
          double qscore = qscores[i];
          if (std::abs(rt - prt) < rt_window_)
          {
            auto inter = t_mass_qscore_map_.find(nominal_mass);
            if (inter == t_mass_qscore_map_.end())
            {
              t_mass_qscore_map_[nominal_mass] = 1 - qscore;
            }
            else
            {
              t_mass_qscore_map_[nominal_mass] *= 1 - qscore;
            }
          }
        }
      }
    }

    for (auto& [m, r] : tqscore_exceeding_mz_rt_map_)
    {
      if (rt - r > rt_window_)
      {
        continue;
      }
      new_mz_rt_map_[m] = r;
    }
    new_mz_rt_map_.swap(tqscore_exceeding_mz_rt_map_);
    std::unordered_map<int, double>().swap(new_mz_rt_map_);

    for (auto& [m, r] : tqscore_exceeding_mass_rt_map_)
    {
      if (rt - r > rt_window_)
      {
        continue;
      }
      new_mass_rt_map_[m] = r;
    }
    new_mass_rt_map_.swap(tqscore_exceeding_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_mass_rt_map_);

    for (auto& item : all_mass_rt_map_)
    {
      if (rt - item.second > rt_window_)
      {
        continue;
      }
      new_all_mass_rt_map_[item.first] = item.second;

      auto inter = new_mass_qscore_map_.find(item.first);

      new_mass_qscore_map_[item.first] = mass_qscore_map_[item.first];

    }
    new_all_mass_rt_map_.swap(all_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_all_mass_rt_map_);

    new_mass_qscore_map_.swap(mass_qscore_map_);
    std::unordered_map<int, double>().swap(new_mass_qscore_map_);

    int selection_phase_start = 0;
    int selection_phase_end = 1; // inclusive
    // When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold.
    // when selection_phase == 1, consider all other masses for selection
    // for target inclusive masses, qscore precursor snr threshold is not applied.
    // In all phase, for target exclusive mode, all the exclusive masses are excluded. For target inclusive mode, only the target masses are considered.

    for (int selection_phase = selection_phase_start; selection_phase <= selection_phase_end; selection_phase++)
    {
      for (auto& pg : deconvolved_spectrum_)
      {
        if (filtered_peakgroups.size() >= mass_count)
        {
          break;
        }

        int charge = pg.getRepAbsCharge();
        double qscore = pg.getQScore();
        double mass = pg.getMonoMass();
        auto [mz1, mz2] = pg.getRepMzRange();
        double center_mz = (mz1 + mz2) / 2.0;

        int nominal_mass = FLASHDeconvAlgorithm::getNominalMass(mass);
        bool target_matched = false;
        double snr_threshold = snr_threshold_;
        double qscore_threshold = qscore_threshold_;

        double tqscore_factor_for_exclusion = 1.0;
        if (targeting_mode_ == 2)
        {
          auto inter = t_mass_qscore_map_.find(nominal_mass);
          if (inter != t_mass_qscore_map_.end())
          {
            tqscore_factor_for_exclusion = t_mass_qscore_map_[nominal_mass];
          }
        }

        if (targeting_mode_ == 1 && target_masses_.size() > 0)  // inclusive mode
        {
          double delta = 2 * tol_[0] * mass * 1e-6;
          auto ub = std::upper_bound(target_masses_.begin(), target_masses_.end(), mass + delta);

          while (!target_matched)
          {
            if (ub != target_masses_.end())
            {
              if (std::abs(*ub - mass) < delta) // target is detected.
              {
                target_matched = true;
              }
              if (mass - *ub > delta)
              {
                break;
              }
            }
            if (ub == target_masses_.begin())
            {
              break;
            }
            ub--;
          }

          if (target_matched)
          {
            snr_threshold = 0.0;
            qscore_threshold = 0.0;
          }
          else
          {
            continue;
          }
        }

        if (qscore < qscore_threshold)
        {
          break;
        }

        if (pg.getChargeSNR(charge) < snr_threshold)
        {
          continue;
        }

        int integer_mz = (int)round(center_mz);

        if (current_selected_mzs.find(integer_mz) != current_selected_mzs.end())
        {
          continue;
        }

        if (selection_phase == 0)
        { // first, select masses under tqscore threshold
          if (1 - tqscore_factor_for_exclusion  > tqscore_threshold ||  tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end() || tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end())
          {
            continue;
          }
        }

        // here crawling isolation windows max_isolation_window_half_
        auto ospec = deconvolved_spectrum_.getOriginalSpectrum();
        if (ospec.size() > 2)
        {
          Size index = ospec.findNearest(center_mz);
          Size tindexl = index == 0 ? index : index - 1;
          Size tindexr = index == 0 ? index + 1 : index;
          double lmz = ospec[tindexl].getMZ(), rmz = ospec[tindexr].getMZ();
          double sig_pwr = .0;
          double noise_pwr = .0;

          bool goleft = tindexl > 0 && (center_mz - lmz >= rmz - center_mz);
          bool goright = tindexr < ospec.size() - 1 && (center_mz - lmz <= rmz - center_mz);

          while (goleft || goright)
          {
            if (goleft)
            {
              double power = pow(ospec[tindexl].getIntensity(), 2);
              if (pg.isSignalMZ(ospec[tindexl].getMZ(), tol_[ospec.getMSLevel() - 1])) //
              {
                sig_pwr += power;
              }
              else
              {
                noise_pwr += power;
              }
              tindexl--;
              lmz = ospec[tindexl].getMZ();
            }
            if (goright)
            {
              double power = pow(ospec[tindexr].getIntensity(), 2);
              if (pg.isSignalMZ(ospec[tindexr].getMZ(), tol_[ospec.getMSLevel() - 1])) //
              {
                sig_pwr += power;
              }
              else
              {
                noise_pwr += power;
              }
              tindexr++;
              rmz = ospec[tindexr].getMZ();
            }

            goleft = tindexl > 0 && (center_mz - lmz >= rmz - center_mz);
            goright = tindexr < ospec.size() - 1 && (center_mz - lmz <= rmz - center_mz);

            if (lmz > mz1 || rmz < mz2 || rmz - lmz < min_isolation_window_half_ * 2)
            {
              continue;
            }

            if (sig_pwr / noise_pwr < snr_threshold)
            {
              break;
            }
          }
          mz1 = std::max(center_mz - max_isolation_window_half_, lmz);
          mz2 = std::min(center_mz + max_isolation_window_half_, rmz);
        }

        if (mz1 < ospec[0].getMZ() - max_isolation_window_half_ || mz2 > ospec.back().getMZ() + max_isolation_window_half_ ||
            mz1 + 2 * min_isolation_window_half_ > mz2 ||
            mz2 - mz1 > 2 * max_isolation_window_half_)
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

        if (1 - mass_qscore_map_[nominal_mass] * tqscore_factor_for_exclusion > tqscore_threshold)
        {
          tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
          tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
        }
        filtered_peakgroups.push_back(pg);
        trigger_charges.push_back(charge);

        trigger_left_isolation_mzs_.push_back(mz1);
        trigger_right_isolation_mzs_.push_back(mz2);

        current_selected_mzs.insert(integer_mz);
      }
    }

    deconvolved_spectrum_.setPeakGroups(filtered_peakgroups);
    std::vector<PeakGroup>().swap(filtered_peakgroups);
  }

  void FLASHIda::getIsolationWindows(double* wstart, double* wend, double* qscores, int* charges, int* min_charges, int* max_charges, double* mono_masses, double* chare_cos, double* charge_snrs,
                                     double* iso_cos, double* snrs, double* charge_scores, double* ppm_errors, double* precursor_intensities, double* peakgroup_intensities)
  {
    // std::sort(deconvolved_spectrum_.begin(), deconvolved_spectrum_.end(), QscoreComparator_);

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

      wstart[i] = trigger_left_isolation_mzs_[i]; // std::get<0>(mz_range) - min_isolation_window_half_;
      wend[i] = trigger_right_isolation_mzs_[i];  // std::get<1>(mz_range) + min_isolation_window_half_;

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
    deconvolved_spectrum_.setPeakGroups(empty);
  }

  MSSpectrum FLASHIda::makeMSSpectrum_(const double* mzs, const double* ints, const int length, const double rt, const int ms_level, const char* name)
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


  std::map<int, std::vector<std::vector<float>>> FLASHIda::parseFLASHIdaLog(const String& in_log_file)
  {
    std::map<int, std::vector<std::vector<float>>> precursor_map_for_real_time_acquisition; // ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color


    if (in_log_file.empty())
    {
      return precursor_map_for_real_time_acquisition;
    }

    std::ifstream f(in_log_file.c_str());
    if (!f.good())
    {
      std::cout << "FLASHIda log file " << in_log_file << " is NOT found. FLASHIda support is not active." << std::endl;
      return precursor_map_for_real_time_acquisition;
    }


    std::cout << "FLASHIda log file used: " << in_log_file << std::endl;
    std::ifstream instream(in_log_file);
    if (instream.good())
    {
      String line;
      int scan;
      float mass, charge, w1, w2, qscore, pint, mint, z1, z2;
      float features[6];
      while (std::getline(instream, line))
      {
        if (line.find("0 targets") != line.npos)
        {
          continue;
        }
        if (line.hasPrefix("MS1"))
        {
          Size st = line.find("MS1 Scan# ") + 10;
          Size ed = line.find(' ', st);
          String n = line.substr(st, ed);
          scan = atoi(n.c_str());
          precursor_map_for_real_time_acquisition[scan] = std::vector<std::vector<float>>(); //// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
        }
        if (line.hasPrefix("Mass"))
        {
          Size st = 5;
          Size ed = line.find('\t');
          String n = line.substr(st, ed);
          mass = (float)atof(n.c_str());

          st = line.find("Z=") + 2;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          charge = (float)atof(n.c_str());

          st = line.find("Score=") + 6;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          qscore = (float)atof(n.c_str());

          st = line.find("[") + 1;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          w1 = (float)atof(n.c_str());

          st = line.find('-', ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          w2 = (float)atof(n.c_str());

          st = line.find("PrecursorIntensity=", ed) + 19;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          pint = (float)atof(n.c_str());

          st = line.find("PrecursorMassIntensity=", ed) + 23;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          mint = (float)atof(n.c_str());

          st = line.find("Features=", ed) + 9;
          // ed = line.find(' ', st);

          st = line.find('[', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[0] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[1] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[2] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[3] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[4] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          features[5] = (float)atof(n.c_str());

          st = line.find("ChargeRange=[", ed) + 13;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          z1 = (float)atof(n.c_str());

          st = line.find("-", ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          z2 = (float)atof(n.c_str());
          std::vector<float> e(15);
          e[0] = mass;
          e[1] = charge;
          e[2] = qscore;
          e[3] = w1;
          e[4] = w2;
          e[5] = pint;
          e[6] = mint;
          e[7] = z1;
          e[8] = z2;
          for (int i = 9; i < 15; i++)
          {
            e[i] = features[i - 9];
          }
          precursor_map_for_real_time_acquisition[scan].push_back(e);
        }
      }
      instream.close();
    }
    else
    {
      std::cout << in_log_file << " not found\n";
    }
    std::cout << "Used precursor size : " << precursor_map_for_real_time_acquisition.size() << std::endl;

    return precursor_map_for_real_time_acquisition;
  }

} // namespace OpenMS
