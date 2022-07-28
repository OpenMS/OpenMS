//--------------------------------------------------------------------------
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

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>

namespace OpenMS
{
  PeakGroup::PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive) :
      min_abs_charge_(min_abs_charge),
      max_abs_charge_(max_abs_charge),
      is_positive_(is_positive)
  {
  }

  PeakGroup::~PeakGroup()
  {
  }

  bool PeakGroup::operator<(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ < a.intensity_;
    }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
    //}
    //return this->spec->getRT() < a.spec->getRT();
  }

  bool PeakGroup::operator>(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ > a.intensity_;
    }
    return this->monoisotopic_mass_ > a.monoisotopic_mass_;
  }

  bool PeakGroup::operator==(const PeakGroup& a) const
  {
    return
        this->monoisotopic_mass_ == a.monoisotopic_mass_
&&        this->intensity_ == a.intensity_;
  }



  void PeakGroup::updateAvgPPMError_(double iso_da_distance)
  {
    std::vector<float> diffs;
    std::vector<std::vector<float>> per_isotope_masses;
    int isotope_end_index = 0;

    for (auto& p : *this)
    {
      if(p.isotopeIndex < 0)
      {
        continue;
      }
      isotope_end_index = isotope_end_index < p.isotopeIndex ? p.isotopeIndex : isotope_end_index;
    }
    per_isotope_masses = std::vector<std::vector<float>>(isotope_end_index + 1, std::vector<float>());
    for (auto& p : *this)
    {
      if(p.isotopeIndex < 0)
      {
        continue;
      }
      per_isotope_masses[p.isotopeIndex].push_back(p.getUnchargedMass());
    }
    diffs.reserve(size());
    for (Size i = 0; i < per_isotope_masses.size(); i++)
    {
      auto& v = per_isotope_masses[i];
      Size n = v.size();
      double average = n >= 2 ? accumulate(v.begin(), v.end(), 0.0) / n : getMonoMass() + i * iso_da_distance; //
      for (float& t : v)
      {
        diffs.push_back(pow(1e6 * (t - average) / average, 2.0));
      }
    }
    Size n = diffs.size();
    avg_ppm_error_ =  n == 0 ? .0 : sqrt(accumulate(diffs.begin(), diffs.end(), 0.0) / n);
  }


  void PeakGroup::updateIsotopeCosineAndQScore(const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double min_cos, double iso_da_distance)
  {
   if(empty())
    {
      return;
    }

    updateChargeFitScoreAndChargeIntensities_();
    if(charge_score_ < .7)
    {
      qscore_ = 0;
      return;
    }

    updateMonomassAndIsotopeIntensities(); //
    int h_offset;
    int second_best_offset= 0;
    isotope_cosine_score_ = FLASHDeconvAlgorithm::getIsotopeCosineAndDetermineIsotopeIndex(monoisotopic_mass_, per_isotope_int_, h_offset, second_best_offset, avg, 0);

    if(isotope_cosine_score_ < min_cos)
    {
      return;
    }

    auto iso_dist = avg.get(monoisotopic_mass_);
    int iso_size = (int)iso_dist.size();
    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeIntensity(abs_charge) <= 0)
      {
        continue;
      }
      auto current_per_isotope_intensities = std::vector<float>(getIsotopeIntensities().size(), .0f);
      int min_isotope_index = current_per_isotope_intensities.size();
      int max_isotope_index = -1; // this is inclusive!!

      for (auto& peak : logMzpeaks_)
      {
        if (peak.abs_charge != abs_charge)
        {
          continue;
        }

        if (peak.isotopeIndex >= current_per_isotope_intensities.size())
        {
          continue;
        }

        if (peak.isotopeIndex < 0)
        {
          continue;
        }
        current_per_isotope_intensities[peak.isotopeIndex] += peak.intensity;
        min_isotope_index = min_isotope_index < peak.isotopeIndex ? min_isotope_index : peak.isotopeIndex;
        max_isotope_index = max_isotope_index < peak.isotopeIndex ? peak.isotopeIndex : max_isotope_index;
      }

      float cos_score = FLASHDeconvAlgorithm::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index + 1, iso_dist, iso_size, 0);

      //cos_score = (float)min_isotope_index / (float)max_isotope_index;//current_per_isotope_intensities[max_isotope_index]/ std::accumulate(current_per_isotope_intensities.begin(), current_per_isotope_intensities.end(), 0);
      setChargeIsotopeCosine(abs_charge, cos_score);//
    }

    updateAvgPPMError_(iso_da_distance);
    updateSNR();
    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) <= 0 || getChargeIsotopeCosine(abs_charge) <= 0)
      {
        continue;
      }

      double q_score = QScore::getQScore(this, abs_charge);
      if(qscore_ > q_score)
      {
        continue;
      }
      max_qscore_abs_charge_ = abs_charge;
      qscore_ = q_score;
    }
    return;
  }

/*
  MSSpectrum PeakGroup::getSubspectrumForMass(const MSSpectrum& spec, const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg, double mono_mass)
  {
    MSSpectrum sub_spec; // is this shallow copy?
    int max_isotope = avg.getLastIndex(mono_mass);
    int left = avg.getLeftCountFromApex(mono_mass);
    int right = avg.getRightCountFromApex(mono_mass);

    sub_spec.reserve(spec.size());

    for(int c = max_abs_charge_;c >= min_abs_charge_;c--)
    {
      if (c <= 0)
      {
        break;
      }
      double cmz = (mono_mass - left) / c + FLASHDeconvHelperStructs::getChargeMass(is_positive_);
      Size index = spec.findNearest(cmz);
      double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / c;
      for (; index < spec.size(); index++)
      {
        double pint = spec[index].getIntensity();
        if (pint <= 0)
        {
          continue;
        }

        int iso_index = (int)round((spec[index].getMZ() - cmz) / iso_delta);
        sub_spec.push_back(spec[index]);

        if (iso_index > max_isotope + right + left)
        {
          break;
        }
      }
    }
    return sub_spec;
  }
*/
  void PeakGroup::recruitAllPeaksInSpectrum(const MSSpectrum& spec, const double tol,
                                            const FLASHDeconvHelperStructs::PrecalculatedAveragine& avg,
                                            double mono_mass, double mass_offset)
  {
    if(mono_mass + mass_offset < 0)
    {
      return;
    }

    int max_isotope = avg.getLastIndex(mono_mass + mass_offset);

    clear();
    reserve((max_isotope) * (max_abs_charge_ - min_abs_charge_ + 1) * 2);
    //noisy_peaks.reserve(max_isotope * (max_abs_charge_ - min_abs_charge_ + 1)* 2);

    int nmax_abs_charge = -1;
    int nmin_abs_charge = min_abs_charge_;
    int max_sig_charge = 0;
    float max_sig = 0;
    for(int c = max_abs_charge_;c >= min_abs_charge_;c--){
      if(c <= 0)
      {
        break;
      }
      double charge_noise_pwr = .0, charge_sig_pwr = .0;
      double charge_intensity = .0;
      double cmz = (mono_mass) /c + FLASHDeconvHelperStructs::getChargeMass(is_positive_);
      Size index = spec.findNearest((cmz + mass_offset/c)* (1 - tol));
      double iso_delta = Constants::ISOTOPE_MASSDIFF_55K_U / c;
      for(;index < spec.size(); index++)
      {
        double pint = spec[index].getIntensity();
        if(pint <= 0)
        {
          continue;
        }
        double pmz = spec[index].getMZ();
        int iso_index = (int)round((pmz - cmz)/iso_delta);
        int iso_index2 = mass_offset !=0? (int)round((pmz - (cmz + mass_offset/c))/iso_delta) : iso_index;

        if(iso_index2 < 0)
        {
          continue;
        }

        if(iso_index2 > max_isotope)
        {
          break;
        }

        if(iso_index < 0)
        {
          continue;
        }

        if(iso_index > max_isotope)
        {
          break;
        }

        double peak_pwr = pint * pint;

        if(abs(pmz - cmz - iso_index * iso_delta) <= std::min(.2, pmz * tol)) //
        {
          auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          push_back(p);

          charge_sig_pwr += peak_pwr;
          charge_intensity += pint;

        }
        else
        {
          /*auto p = LogMzPeak(spec[index], is_positive_);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          noisy_peaks.push_back(p); //tmp
          */
          charge_noise_pwr += peak_pwr;
        }
      }
      if(charge_sig_pwr > 0)
      {
        if(nmax_abs_charge < 0)
        {
          nmax_abs_charge = c;
        }
        nmin_abs_charge = c;

        setChargePowers_(c, charge_sig_pwr, charge_noise_pwr, charge_intensity);

        if(max_sig < charge_sig_pwr)
        {
          max_sig = charge_sig_pwr;
          max_sig_charge = c;
        }
      }

      if(index >= spec.size())
      {
        break;
      }
    }

    // determine charge ranges
    if(max_sig_charge > 0)
    {
      int t_nmax_abs_charge = nmax_abs_charge;
      int t_nmin_abs_charge = nmin_abs_charge;
      int c_zero = 0;
      for(int z=max_sig_charge;z<=t_nmax_abs_charge;z++)
      {
        if(per_charge_int_[z] <= 0)
        {
          c_zero ++;
        }else{
          c_zero = 0;
        }
        if(c_zero > 2)
        {
          break;
        }
        nmax_abs_charge = z;
      }
      c_zero = 0;
      for(int z=max_sig_charge;z>=t_nmin_abs_charge;z--)
      {
        if(per_charge_int_[z] <= 0)
        {
          c_zero ++;
        }else{
          c_zero = 0;
        }
        if(c_zero > 2)
        {
          break;
        }
        nmin_abs_charge = z;
      }
    }
    sort();
    max_abs_charge_ = nmax_abs_charge;
    min_abs_charge_ = nmin_abs_charge;

    if(min_abs_charge_>max_abs_charge_)
    {
      clear();
    }

  }

  void PeakGroup::updateChargeFitScoreAndChargeIntensities_()
  {
    double max_per_charge_intensity = .0;
    double summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = -1;

    for (int c = min_abs_charge_;c<=max_abs_charge_;c++)
    {
      summed_intensity += per_charge_int_[c];
      if (per_charge_int_[c] > 0)
      {
        if (first_index < 0)
        {
          first_index = c;
        }
        last_index = c;
      }

      if (max_per_charge_intensity > per_charge_int_[c])
      {
        continue;
      }
      max_per_charge_intensity = per_charge_int_[c];
      max_index = c;
    }
    if(max_index <0 || summed_intensity <= 0)
    {
      charge_score_ = 0;
      return;
    }

    first_index = first_index < 0 ? 0 : first_index;
    min_abs_charge_ = first_index;
    max_abs_charge_ = last_index;
    double p = .0;
    for (int c = max_index; c < last_index; c++) //
    {
      double diff = per_charge_int_[c + 1] - per_charge_int_[c];
      //double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i + 1]);
      if (diff <= 0)
      {
        continue;
      }
      p += diff;
    }
    for (int c = max_index; c > first_index; c--)
    {
      double diff = per_charge_int_[c - 1] - per_charge_int_[c];
      //      double ratio = per_charge_intensity[i] / (.1 + per_charge_intensity[i - 1]);

      if (diff <= 0)
      {
        continue;
      }
      p += diff;
    }
    charge_score_ = std::max(.0, 1.0 - p / summed_intensity);
  }

  void PeakGroup::updateMonomassAndIsotopeIntensities()
  {
    int max_isotope_index = 0;
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
    for (auto& p: logMzpeaks_)
    {
      max_isotope_index = max_isotope_index < p.isotopeIndex ? p.isotopeIndex : max_isotope_index;
    }

    per_isotope_int_ = std::vector<float>(max_isotope_index + 1, .0f);
    intensity_ = .0;
    double nominator = .0;

    //std::vector<LogMzPeak> new_logMzpeaks_;
    //new_logMzpeaks_.reserve(logMzpeaks_.size());

    for (auto& p: logMzpeaks_)
    {
      double pi = p.intensity + 1;
      if (p.isotopeIndex < 0)
      {
        continue;
      }

      per_isotope_int_[p.isotopeIndex] += p.intensity;
     // new_logMzpeaks_.push_back(p);
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * Constants::ISOTOPE_MASSDIFF_55K_U);
      intensity_ += pi;
    }
    //logMzpeaks_.swap(new_logMzpeaks_);

    monoisotopic_mass_ = nominator / intensity_;
  }

  bool PeakGroup::isSignalMZ(const double mz, const double tol) const
  {
    for (auto& p: logMzpeaks_)
    {
      if (abs(p.mz - mz) < p.mz * tol * 1e-6)
      {
        return true;
      }
    }
    return false;
  }

  void PeakGroup::setSNR(const float snr)
  {
    snr_ = snr;
  }

  void PeakGroup::setChargeSNR(const int abs_charge, const float c_snr)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_snr_.empty())
    {
      per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_snr_[abs_charge] = c_snr;
  }

  void PeakGroup::setTargeted()
  {
    is_targeted_ = true;
  }

  void PeakGroup::setChargePowers_(const int abs_charge, const double signal_pwr, const double noise_pwr, const double intensity)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_pwr_.empty())
    {
      per_charge_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    if (per_charge_signal_pwr_.empty())
    {
      per_charge_signal_pwr_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    if (per_charge_int_.empty())
    {
      per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_int_[abs_charge] = intensity;
    per_charge_pwr_[abs_charge] = signal_pwr + noise_pwr;
    per_charge_signal_pwr_[abs_charge] = signal_pwr;
  }


  void PeakGroup::setChargeIsotopeCosine(const int abs_charge, const float cos)
  {
    if (max_abs_charge_ < abs_charge)
    {
      return;
    }
    if (per_charge_cos_.empty())
    {
      per_charge_cos_ = std::vector<float>(1 + max_abs_charge_, .0);
    }
    per_charge_cos_[abs_charge] = cos;
  }


  void PeakGroup::setMaxQScoreMzRange(const double min, const double max)
  {
    max_qscore_mz_start_ = min;
    max_qscore_mz_end_ = max;
  }

  void PeakGroup::setAbsChargeRange(const int min_abs_charge, const int max_abs_charge)
  {
    min_abs_charge_ = min_abs_charge;
    max_abs_charge_ = max_abs_charge;
  }

  void PeakGroup::setScanNumber(const int sn)
  {
    scan_number_ = sn;
  }

  void PeakGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
  }

  void PeakGroup::setRepAbsCharge(const int max_qscore_charge)
  {
    max_qscore_abs_charge_ = max_qscore_charge;
  }

  void PeakGroup::setChargeScore(const float score)
  {
    charge_score_ = score;
  }

  void PeakGroup::setAvgPPMError(const float error)
  {
    avg_ppm_error_ = error;
  }

  void PeakGroup::setQScore(const float q)
  {
    qscore_ = q;
  }

  std::tuple<double, double> PeakGroup::getMaxQScoreMzRange() const
  {
    return std::tuple<double, double>{max_qscore_mz_start_, max_qscore_mz_end_};
  }

  std::tuple<double, double> PeakGroup::getMzRange(int abs_charge) const
  {
    double mz_start = -1;
    double mz_end = -10;
    if (!(abs_charge > max_abs_charge_ || abs_charge < min_abs_charge_))
    {
      for (auto& tmp_p: logMzpeaks_)
      {
        if (tmp_p.abs_charge != abs_charge)
        {
          continue;
        }
        if (mz_start < 0)
        {
          mz_start = tmp_p.mz;
        }
        else
        {
          mz_start = mz_start < tmp_p.mz ? mz_start : tmp_p.mz;
        }
        mz_end = mz_end > tmp_p.mz ? mz_end : tmp_p.mz;
      }
    }
    return std::tuple<double, double>{mz_start, mz_end};
  }

  std::tuple<int, int> PeakGroup::getAbsChargeRange() const
  {
    return std::tuple<int, int>{min_abs_charge_, max_abs_charge_};
  }

  std::vector<float> PeakGroup::getIsotopeIntensities() const
  {
    return per_isotope_int_;
  }

  int PeakGroup::getScanNumber() const
  {
    return scan_number_;
  }

  double PeakGroup::getMonoMass() const
  {
    return monoisotopic_mass_;
  }

  double PeakGroup::getIntensity() const
  {
    return intensity_;
  }

  float PeakGroup::getIsotopeCosine() const
  {
    return isotope_cosine_score_;
  }

  int PeakGroup::getRepAbsCharge() const
  {
    return max_qscore_abs_charge_;
  }

  float PeakGroup::getQScore() const
  {
    return qscore_;
  }

  bool PeakGroup::isTargeted() const
  {
    return is_targeted_;
  }

  void PeakGroup::updateSNR()
  {
    float cos_squred = isotope_cosine_score_ * isotope_cosine_score_;
    float signal = 0, noise = 0;
    per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0);

    for (int c = min_abs_charge_; c <= std::min((int) per_charge_signal_pwr_.size() - 1, max_abs_charge_); ++c)
    {
      float signal_pwr =
          per_charge_signal_pwr_[c] < per_charge_pwr_[c] ? per_charge_signal_pwr_[c] : per_charge_pwr_[c];
      if(per_charge_cos_.size() > c)
      {
        float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
        float nom = per_charge_cos_squared * signal_pwr;
        float denom = 1 + per_charge_pwr_[c] - signal_pwr + (1 - per_charge_cos_squared) * signal_pwr;

        per_charge_snr_[c] = denom <= 0 ? .0 : (nom / denom);
      }
      signal += per_charge_signal_pwr_[c];
      noise += per_charge_pwr_[c] - signal_pwr;
    }

    float t_nom = cos_squred * signal;
    float t_denom = 1 + noise
                    + (1 - cos_squred) * signal;

    snr_ = t_denom <= 0 ? .0 : (t_nom / t_denom);
  }


  float PeakGroup::getSNR() const
  {
    return snr_;
  }

  float PeakGroup::getChargeScore() const
  {
    return charge_score_;
  }

  float PeakGroup::getAvgPPMError() const
  {
    return avg_ppm_error_;
  }


  float PeakGroup::getChargeSNR(const int abs_charge) const
  {
    if (abs_charge<0 || (int) per_charge_snr_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_snr_[abs_charge];
  }

  float PeakGroup::getChargeIsotopeCosine(const int abs_charge) const
  {
    if (abs_charge<0 || (int) per_charge_cos_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_cos_[abs_charge];
  }

  float PeakGroup::getChargeIntensity(const int abs_charge) const
  {
    if (abs_charge<0 || (int) per_charge_int_.size() <= abs_charge)
    {
      return 0;
    }
    return per_charge_int_[abs_charge];
  }

  bool PeakGroup::isPositive() const
  {
    return is_positive_;
  }



  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator PeakGroup::begin() const noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::const_iterator PeakGroup::end() const noexcept
  {
    return logMzpeaks_.end();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator PeakGroup::begin() noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHDeconvHelperStructs::LogMzPeak>::iterator PeakGroup::end() noexcept
  {
    return logMzpeaks_.end();
  }

  const FLASHDeconvHelperStructs::LogMzPeak& PeakGroup::operator[](const Size i) const
  {
    return logMzpeaks_[i];
  }

  void PeakGroup::push_back (const FLASHDeconvHelperStructs::LogMzPeak& pg)
  {
    logMzpeaks_.push_back(pg);
  }

  Size PeakGroup::size() const noexcept
  {
    return logMzpeaks_.size();
  }

  void PeakGroup::clear()
  {
    logMzpeaks_.clear();
  }

  void PeakGroup::reserve (Size n)
  {
    logMzpeaks_.reserve(n);
  }

  bool PeakGroup::empty() const
  {
    return logMzpeaks_.empty();
  }

  void PeakGroup::swap (std::vector<FLASHDeconvHelperStructs::LogMzPeak>& x)
  {
    logMzpeaks_.swap(x);
  }

  void PeakGroup::shrink_to_fit()
  {
    logMzpeaks_.shrink_to_fit();
  }

  void PeakGroup::sort()
  {
    std::sort(logMzpeaks_.begin(), logMzpeaks_.end());
  }
}
