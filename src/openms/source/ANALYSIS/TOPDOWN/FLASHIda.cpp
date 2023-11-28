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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
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
  inline const double max_isolation_window_half_ = 2.5;
  /// constructor
  FLASHIda::FLASHIda(char* arg)
  {
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
    std::unordered_map<std::string, std::vector<double>> inputs;
    std::vector<String> log_files;
    std::vector<String> out_files; /// add tsv for exclusion list in the future.
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
      std::cout << ss.str() << "file(s) is(are) used for inclusion mode\n";
    }
    else if (targeting_mode_ == 2)
    {
      std::cout << ss.str() << "file(s) is(are) used for in-depth mode\n";
    }else if (targeting_mode_ == 3)
    {
      std::cout << ss.str() << "file(s) is(are) used for exclusion mode\n";
    }
    Param sd_defaults = SpectralDeconvolution().getDefaults();

    sd_defaults.setValue("min_charge", (int)inputs["min_charge"][0]);
    sd_defaults.setValue("max_charge", (int)inputs["max_charge"][0]);
    sd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    sd_defaults.setValue("max_mass", inputs["max_mass"][0]);
    sd_defaults.setValue("tol", inputs["tol"]);
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

            if (targeting_mode_ == 1 || targeting_mode_ == 2)
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
          }
          else if (line.hasPrefix("AllMass"))
          {
            if (targeting_mode_ == 3)
            {
              Size st = 8;
              Size ed = line.size();
              String n = line.substr(st, ed - st + 1);

              std::stringstream tmp_stream(n);
              String str;
              std::vector<double> results;
              while (getline(tmp_stream, str, ' '))
              {
                results.push_back(atof(str.c_str()));
              }
              if (exclusion_rt_masses_map_.find(rt * 60.0) == exclusion_rt_masses_map_.end())
              {
                exclusion_rt_masses_map_[rt * 60.0] = std::vector<double>();
              }
              for(double m : results)
              {
                exclusion_rt_masses_map_[rt * 60.0].push_back(m);
              }
            }
          }
        }
        instream.close();
      }
    }

    for (auto& log_file : out_files)
    {
      std::ifstream instream(log_file);
      double rt = .0;
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

          if (targeting_mode_ == 1 || targeting_mode_ == 2)
          {
            if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
            {
              target_mass_rt_map_[mass] = std::vector<double>();
            }
            rt = atof(results[0].c_str());
            target_mass_rt_map_[mass].push_back(60.0 * rt);
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

    fd_.setParameters(sd_defaults);
    fd_.calculateAveragine(false);

    std::cout << sd_defaults << std::endl;
  }

  int FLASHIda::getPeakGroups(const double* mzs, const double* ints,
                              const int length, const double rt, const int ms_level, const char* name)
  {
    //int ret[2] = {0,0};
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);
    // selected_peak_groups_ = DeconvolvedSpectrum(spec, 1);
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
    PeakGroup empty;

    target_masses_.clear();
    excluded_masses_.clear();
    if (targeting_mode_ == 1)
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
    }else if (targeting_mode_ == 3)
    {
      for (auto& [prt, masses] : exclusion_rt_masses_map_)
      {
        if (std::abs(rt - prt) >= rt_window_ && prt != 0)
          continue;
        for (double mass : masses)
        {
          excluded_masses_.push_back(mass);
        }
      }
      std::sort(excluded_masses_.begin(), excluded_masses_.end());
    }

    selected_peak_groups_.clear();
    deconvolved_spectrum_.clear();

    fd_.performSpectrumDeconvolution(spec, 0, empty);
    deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();
    // per spec deconvolution
    FLASHIda::filterPeakGroupsUsingMassExclusion_(ms_level, rt);
    // spec.clear(true);
    return (int)selected_peak_groups_.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    deconvolved_spectrum_.sortByQscore();
    Size mass_count = (Size)mass_count_[ms_level - 1];
    trigger_charges.clear();
    trigger_charges.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);

    selected_peak_groups_.reserve(mass_count_.size());
    std::set<double> current_selected_mzs; // current selected mzs
    std::set<double> current_selected_masses; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_qscore_map_;
    std::unordered_map<int, double> t_mass_qscore_map_;

    if (targeting_mode_ == 2)
    {
      for (auto& [mass, rts] : target_mass_rt_map_)
      {
        int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
        auto qscores = target_mass_qscore_map_[mass];
        for (uint i = 0; i < rts.size(); i++)
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

      //auto inter = new_mass_qscore_map_.find(item.first);

      new_mass_qscore_map_[item.first] = mass_qscore_map_[item.first];

    }
    new_all_mass_rt_map_.swap(all_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_all_mass_rt_map_);

    new_mass_qscore_map_.swap(mass_qscore_map_);
    std::unordered_map<int, double>().swap(new_mass_qscore_map_);

    const int selection_phase_start = 0;
    const int selection_phase_end = 2; // inclusive
    // When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold.
    // when selection_phase == 1, consider all other masses for selection but the same m/z is avoided
    // when selection_phase == 2, consider all.
    // for target inclusive masses, qscore precursor snr threshold is not applied.
    // In all phase, for target exclusive mode, all the exclusive masses are excluded. For target inclusive mode, only the target masses are considered.

    for (int iteration = targeting_mode_ == 2? 0 : 1; iteration < 2 ;iteration++) // for mass exclusion, first collect masses with exclusion list. Then collect without exclusion. This works the best
    {
      for (int selection_phase = selection_phase_start; selection_phase <= selection_phase_end; selection_phase++)
      {
        for (auto& pg : deconvolved_spectrum_)
        {
          if (selected_peak_groups_.size() >= mass_count)
          {
            break;
          }

          int charge = pg.getRepAbsCharge();
          double qscore = pg.getQscore();
          double mass = pg.getMonoMass();
          auto [mz1, mz2] = pg.getRepMzRange();
          double center_mz = (mz1 + mz2) / 2.0;

          int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
          bool target_matched = false;
          double snr_threshold = snr_threshold_;
          double qscore_threshold = qscore_threshold_;
          double tqscore_factor_for_exclusion = 1.0;
          int integer_mz = (int)round(center_mz);

          if (iteration == 0)
          {
            auto inter = t_mass_qscore_map_.find(nominal_mass);
            if (inter != t_mass_qscore_map_.end())
            {
              tqscore_factor_for_exclusion = t_mass_qscore_map_[nominal_mass];
            }
            if (1 - tqscore_factor_for_exclusion > tqscore_threshold)
            {
              continue;
            }
          }

          if (targeting_mode_ == 1 && target_masses_.size() > 0) // inclusive mode
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
              qscore_threshold = 0.0; // stop exclusion for targets. TODO tqscore lowest first? charge change.
            }
            else
            {
              continue;
            }
          }else if (targeting_mode_ == 3 && excluded_masses_.size() > 0) // inclusive mode
          {
            bool to_exclude = false;
            double delta = 2 * tol_[0] * mass * 1e-6;
            auto ub = std::upper_bound(excluded_masses_.begin(), excluded_masses_.end(), mass + delta);

            while (!to_exclude)
            {
              if (ub != excluded_masses_.end())
              {
                if (std::abs(*ub - mass) < delta) // target is detected.
                {
                  to_exclude = true;
                }
                if (mass - *ub > delta)
                {
                  break;
                }
              }
              if (ub == excluded_masses_.begin())
              {
                break;
              }
              ub--;
            }

            if (to_exclude)
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


          if (current_selected_mzs.find(center_mz) != current_selected_mzs.end()) // mz has been triggered
          {
            if (selection_phase < selection_phase_end)
            {
              continue;
            }
            if (!target_matched && current_selected_masses.find(pg.getMonoMass()) == current_selected_masses.end()) // but mass is different
            {
              continue;
            }
          }

          if (selection_phase < selection_phase_end - 1)
          { // first, select masses under tqscore threshold
            if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end() || tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end())
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
            double lmz = center_mz, rmz = center_mz;
            double sig_int = pg.getChargeIntensity(charge);
            double noise_int = sqrt(sig_int * sig_int / (.01 + pg.getChargeSNR(charge)));

            bool goleft = tindexl > 0 && (center_mz - lmz <= rmz - center_mz);
            bool goright = tindexr < ospec.size() - 1 && (center_mz - lmz >= rmz - center_mz);

            while (goleft || goright)
            {
              if (goleft)
              {
                tindexl--;
                lmz = ospec[tindexl].getMZ();
                double intensity = ospec[tindexl].getIntensity();
                if (lmz < mz1)
                {
                  noise_int += intensity;
                }
              }
              if (goright)
              {
                tindexr++;
                rmz = ospec[tindexr].getMZ();
                double intensity = ospec[tindexr].getIntensity();
                if (rmz > mz2)
                {
                  noise_int += intensity;
                }
              }

              goleft = tindexl > 0 && (center_mz - lmz <= rmz - center_mz);
              goright = tindexr < ospec.size() - 1 && (center_mz - lmz >= rmz - center_mz);

              if (lmz > mz1 || rmz < mz2 || rmz - lmz < min_isolation_window_half_ * 2)
              {
                continue;
              }

              if (sig_int / noise_int < sqrt(snr_threshold))
              {
                break;
              }
            }
            mz1 = std::max(center_mz - max_isolation_window_half_, lmz);
            mz2 = std::min(center_mz + max_isolation_window_half_, rmz);

          }


          if (mz1 < ospec[0].getMZ() - max_isolation_window_half_ || mz2 > ospec.back().getMZ() + max_isolation_window_half_ || mz1 + 2 * min_isolation_window_half_ - .01 > mz2 ||
              mz2 - mz1 > 2 * max_isolation_window_half_ + .01)
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

          selected_peak_groups_.push_back(pg);
          trigger_charges.push_back(charge);

          trigger_left_isolation_mzs_.push_back(mz1);
          trigger_right_isolation_mzs_.push_back(mz2);
          current_selected_masses.insert(pg.getMonoMass());
          current_selected_mzs.insert(center_mz);
        }
      }
    }
  }

  void FLASHIda::getAllMonoisotopicMasses(double *masses, int length)
  {
    int len = std::min(length, (int)deconvolved_spectrum_.size());
    for(int i=0;i<len;i++)
    {
      masses[i] = deconvolved_spectrum_[i].getMonoMass();
    }
  }

  int FLASHIda::GetAllPeakGroupSize()
  {
    return deconvolved_spectrum_.size();
  }

  void FLASHIda::getIsolationWindows(double* wstart, double* wend, double* qscores, int* charges, int* min_charges, int* max_charges, double* mono_masses, double* chare_cos, double* charge_snrs,
                                     double* iso_cos, double* snrs, double* charge_scores, double* ppm_errors, double* precursor_intensities, double* peakgroup_intensities)
  {
    // std::sort(selected_peak_groups_.begin(), selected_peak_groups_.end(), QscoreComparator_);

    for (Size i = 0; i < selected_peak_groups_.size(); i++)
    {
      if (trigger_charges[i] == 0)
      {
        continue;
      }
      auto peakgroup = selected_peak_groups_[i];
      charges[i] = trigger_charges[i];
      auto cr = peakgroup.getAbsChargeRange();
      min_charges[i] = std::get<0>(cr);
      max_charges[i] = std::get<1>(cr);

      wstart[i] = trigger_left_isolation_mzs_[i]; // std::get<0>(mz_range) - min_isolation_window_half_;
      wend[i] = trigger_right_isolation_mzs_[i];  // std::get<1>(mz_range) + min_isolation_window_half_;

      qscores[i] = Qscore::getQscore(&peakgroup);
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
    int mass_cntr = 0;
    for(auto& v : precursor_map_for_real_time_acquisition)
    {
      mass_cntr+=v.second.size();
    }

    std::cout << "Used precursor size : " << precursor_map_for_real_time_acquisition.size() << " precursor masses : " << mass_cntr <<  std::endl;

    return precursor_map_for_real_time_acquisition;
  }
} // namespace OpenMS
