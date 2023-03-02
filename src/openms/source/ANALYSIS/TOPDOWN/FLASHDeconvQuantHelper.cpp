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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvQuantHelper.h>

namespace FLASHDeconvQuantHelper
{
  /// comparison operators (using monoisotopic_mass_)
  bool FeatureGroup::operator<(const FeatureGroup &a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ < a.intensity_;
    }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
  }

  bool FeatureGroup::operator>(const FeatureGroup &a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
    {
      return this->intensity_ > a.intensity_;
    }
    return this->monoisotopic_mass_ > a.monoisotopic_mass_;
  }

  bool FeatureGroup::operator==(const FeatureGroup &a) const
  {
    return
        this->monoisotopic_mass_ == a.monoisotopic_mass_
        && this->intensity_ == a.intensity_;
  }

  /// iterator related functions
  std::vector<FeatureSeed>::const_iterator FeatureGroup::begin() const noexcept
  {
    return feature_seeds_.cbegin();
  }

  std::vector<FeatureSeed>::const_iterator FeatureGroup::end() const noexcept
  {
    return feature_seeds_.cend();
  }

  std::vector<FeatureSeed>::iterator FeatureGroup::begin() noexcept
  {
    return feature_seeds_.begin();
  }

  std::vector<FeatureSeed>::iterator FeatureGroup::end() noexcept
  {
    return feature_seeds_.end();
  }

  const FeatureSeed& FeatureGroup::operator[](const Size i) const
  {
    return feature_seeds_[i];
  }

  void FeatureGroup::push_back(const FeatureSeed &new_feature)
  {
    feature_seeds_.push_back(new_feature);
  }

  Size FeatureGroup::size() const noexcept
  {
    return feature_seeds_.size();
  }

  void FeatureGroup::reserve (Size n)
  {
    feature_seeds_.reserve(n);
  }

  void FeatureGroup::clear()
  {
    feature_seeds_.clear();
  }

  std::vector<FeatureSeed>::iterator FeatureGroup::erase(std::vector<FeatureSeed>::iterator pos)
  {
    return feature_seeds_.erase(pos);
  }

  bool FeatureGroup::empty() const
  {
    return feature_seeds_.empty();
  }

  void FeatureGroup::swap (std::vector<FeatureSeed>& x)
  {
    feature_seeds_.swap(x);
  }

  void FeatureGroup::sort()
  {
    std::sort(feature_seeds_.begin(), feature_seeds_.end());
  }

  /// default getter and setters

  double FeatureGroup::getMonoisotopicMass() const
  {
    return monoisotopic_mass_;
  }

  int FeatureGroup::getMinCharge() const
  {
    return min_abs_charge_;
  }

  int FeatureGroup::getMaxCharge() const
  {
    return max_abs_charge_;
  }

  Size FeatureGroup::getMaxIsotopeIndex() const
  {
    return max_isotope_index_;
  }

  double FeatureGroup::getIntensity() const
  {
    return intensity_;
  }

  double FeatureGroup::getRtOfMostAbundantMT() const
  {
    return centroid_rt_of_most_abundant_mt_;
  }

  float FeatureGroup::getIsotopeCosine() const
  {
    return isotope_cosine_score_;
  }

  float FeatureGroup::getFeatureGroupScore() const
  {
    return total_score_;
  }

  const std::set<int>& FeatureGroup::getChargeSet() const
  {
    return charges_;
  }

  const std::pair<double, double>& FeatureGroup::getFwhmRange() const
  {
    return fwhm_range_;
  }

  const std::vector<Size>& FeatureGroup::getTraceIndices() const
  {
    return ltrace_indices_;
  }

  const std::vector<float>& FeatureGroup::getIsotopeIntensities() const
  {
    return per_isotope_int_;
  }

  const std::vector<float>& FeatureGroup::getChargeIntensities() const
  {
    return per_charge_int_;
  }

  float FeatureGroup::getIntensityOfCharge(const int &abs_charge) const
  {
    return per_charge_int_[abs_charge];
  }

  float FeatureGroup::getIsotopeCosineOfCharge(const int &abs_charge) const
  {
    return per_charge_cos_[abs_charge];
  }

  double FeatureGroup::getAverageMass() const
  {
    return average_mass_;
  }

  void FeatureGroup::setMonoisotopicMass(const double mass)
  {
    monoisotopic_mass_ = mass;
  }

  void FeatureGroup::setChargeRange(const int min_c, const int max_c)
  {
    min_abs_charge_ = min_c;
    max_abs_charge_ = max_c;
  }

  void FeatureGroup::setMaxIsotopeIndex(const Size index)
  {
    max_isotope_index_ = index;
  }

  void FeatureGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
  }

  void FeatureGroup::setFeatureGroupScore(const float score)
  {
    total_score_ = score;
  }

  void FeatureGroup::setPerChargeIntensities(std::vector<float> const &perChargeInt)
  {
    per_charge_int_ = perChargeInt;
  }

  void FeatureGroup::setPerChargeCosineScore(std::vector<float> const &perChargeCos)
  {
    per_charge_cos_ = perChargeCos;
  }

  void FeatureGroup::setAverageMass(double averageMass)
  {
    average_mass_ = averageMass;
  }

  void FeatureGroup::updateMembers()
  {
    /// --- Excluded members for updates ----
    // monoisotopic_mass_ & intensity_ & per_isotope_int_ & max_isotope_index_ & charges_ : used for scoring & filtering, thus should be treated separately
    // all score related members

    // TODO: check if sorting is needed here
    std::sort(feature_seeds_.begin(), feature_seeds_.end());

    // charge range
    min_abs_charge_ = *charges_.begin();
    max_abs_charge_ = *charges_.rbegin();

    // members to be changed or trackers
    FeatureSeed *most_abundant_seed = &feature_seeds_[0];
    double min_fwhm(std::numeric_limits<double>::max());
    double max_fwhm(0.0);
    ltrace_indices_.clear();
    ltrace_indices_.reserve(feature_seeds_.size());

    for(auto &s : feature_seeds_)
    {
      // find the most abundant seed
      if (most_abundant_seed->getIntensity() < s.getIntensity())
      {
        most_abundant_seed = &s;
      }

      // fwhm
      std::pair<double, double> tmp_fwhm(s.getFwhmStart(), s.getFwhmEnd());
      if (tmp_fwhm.first < min_fwhm)
        min_fwhm = tmp_fwhm.first;
      if (tmp_fwhm.second > max_fwhm)
        max_fwhm = tmp_fwhm.second;

      // update ltrace_indices_
      ltrace_indices_.push_back(s.getTraceIndex());
    }

    // find the most abundant peak
    Size max_peak_idx = most_abundant_seed->getMassTrace().findMaxByIntPeak(false);
    auto tmp_max_peak = most_abundant_seed->getMassTrace()[max_peak_idx];
    centroid_rt_of_most_abundant_mt_ = tmp_max_peak.getRT();

    fwhm_range_ = std::make_pair(min_fwhm, max_fwhm);

    // for fast searching later
    std::sort(ltrace_indices_.begin(), ltrace_indices_.end());
  }

  void FeatureGroup::updateMembersForScoring()
  {
    // based on PeakGroup::updateMonomassAndIsotopeIntensities()
    /// update 5 members: monoisotopic_mass_, max_isotope_index_, per_isotope_int_, intensity_, charges_

    // calculate max_isotope_index_
    int max_isotope_index = 0;
    for (auto& f: feature_seeds_)
    {
      max_isotope_index = max_isotope_index < f.getIsotopeIndex() ? f.getIsotopeIndex() : max_isotope_index;
    }
    max_isotope_index_ = max_isotope_index;

    per_isotope_int_ = std::vector<float>(max_isotope_index_ + 1, .0f);
    intensity_ = .0;
    double nominator = .0;

    // update per_isotope_int_, intensity_ and charge_
    for (auto &f: feature_seeds_)
    {
      if (f.getIsotopeIndex() < 0)
      {
        continue;
      }
      per_isotope_int_[f.getIsotopeIndex()] += f.getIntensity();
      charges_.insert(f.getCharge());

      double pi = f.getIntensity() + 1;
      intensity_ += pi;
      nominator += pi * (f.getUnchargedMass() - f.getIsotopeIndex() * Constants::ISOTOPE_MASSDIFF_55K_U);
    }
    // update monoisotopic mass
    monoisotopic_mass_ = nominator / intensity_;
  }

  void FeatureGroup::updateIsotopeIndices(const int offset)
  {
    for (auto &seed : feature_seeds_)
    {
      seed.setIsotopeIndex(seed.getIsotopeIndex() - offset);
    }
  }

  // TODO: need to find a smarter way
  bool FeatureGroup::doesThisChargeExist(int charge) const
  {
    bool exist = false;
    for (auto lmt_iter = feature_seeds_.begin(); lmt_iter != feature_seeds_.end(); ++lmt_iter)
    {
      if (lmt_iter->getCharge() == charge)
      {
        exist = true;
      }
    }
    return exist;
  }

  FeatureSeed* FeatureGroup::getApexLMT() const
  {
    FeatureSeed* apex_lmt = nullptr;
    double max_intensity = 0.0;

    for (auto lmt_iter = feature_seeds_.begin(); lmt_iter != feature_seeds_.end(); ++lmt_iter)
    {
      if (lmt_iter->getIntensity() > max_intensity)
      {
        max_intensity = lmt_iter->getIntensity();
        apex_lmt = (FeatureSeed*) &(*lmt_iter);
      }
    }
    return apex_lmt;
  }
}
