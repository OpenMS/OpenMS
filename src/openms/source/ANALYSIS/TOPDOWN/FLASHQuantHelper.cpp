// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantHelper.h>

namespace FLASHQuantHelper
{
  /// getter & setter
  const MassTrace& FeatureSeed::getMassTrace() const
  {
    return mass_trace_;
  }

  double FeatureSeed::getCentroidMz() const
  {
    return centroid_mz_;
  }

  int FeatureSeed::getCharge() const
  {
    return charge_;
  }

  double FeatureSeed::getFwhmStart() const
  {
    return fwhm_start_;
  }

  double FeatureSeed::getFwhmEnd() const
  {
    return fwhm_end_;
  }

  double FeatureSeed::getIntensity() const
  {
    return intensity_;
  }

  int FeatureSeed::getIsotopeIndex() const
  {
    return isotope_index_;
  }

  Size FeatureSeed::getTraceIndex() const
  {
    return trace_index_;
  }

  double FeatureSeed::getMass() const
  {
    return mass_;
  }

  void FeatureSeed::setMassTrace(MassTrace &mt)
  {
    mass_trace_ = mt;
  }

  void FeatureSeed::setCentroidMz(double &mz)
  {
    centroid_mz_ = mz;
  }

  void FeatureSeed::setCharge(int cs)
  {
    charge_ = cs;
  }

  void FeatureSeed::setFwhmStart(double fwhm_s)
  {
    fwhm_start_ = fwhm_s;
  }

  void FeatureSeed::setFwhmEnd(double fwhm_e)
  {
    fwhm_end_ = fwhm_e;
  }

  void FeatureSeed::setIntensity(double inty)
  {
    intensity_ = inty;
  }

  void FeatureSeed::setIsotopeIndex(int idx)
  {
    isotope_index_ = idx;
  }

  void FeatureSeed::setTraceIndex(Size i)
  {
    trace_index_ = i;
  }

  void FeatureSeed::setMass(double mass)
  {
    mass_ = mass;
  }

  double FeatureSeed::getUnchargedMass()
  {
    if (charge_ == 0)
    {
      return .0;
    }
    if (mass_ <= 0)
    {
      mass_ = (centroid_mz_ - Constants::PROTON_MASS_U) * charge_;
    }
    return mass_;
  }

  /// referenced: MassTrace::estimateFWHM
  std::pair<Size, Size> FeatureSeed::computeBulkRetentionTimeRange(bool use_smoothed_ints) const
  {
    /// calculating retention time of 10% of maximum (Apex)
    Size max_idx(mass_trace_.findMaxByIntPeak(use_smoothed_ints));

    std::vector<double> tmp_ints;
    for (Size vec_idx = 0; vec_idx < mass_trace_.getSize(); ++vec_idx)
    {
      tmp_ints.push_back(mass_trace_[vec_idx].getIntensity());
    }

    double inty_threshold(tmp_ints[max_idx] * 0.1); // 10 % of maximum

    // mass trace is empty OR no points left of apex in mass trace OR no points right of apex in mass trace
    if (tmp_ints.empty() || max_idx == 0 || max_idx == tmp_ints.size() - 1)
    {
      return std::make_pair(0, 0);
    }

    Size left_border(max_idx), right_border(max_idx);

    while (left_border > 0 && tmp_ints[left_border] >= inty_threshold)
    {
      --left_border;
    }

    while (right_border + 1 < tmp_ints.size() && tmp_ints[right_border] >= inty_threshold)
    {
      ++right_border;
    }

    return std::make_pair(left_border, right_border);
  }

  /// referenced: MassTrace::computeFwhmArea()
  double FeatureSeed::computeBulkPeakArea(bool use_smoothed_ints) const
  {
    if (mass_trace_.getSize()==0) // if empty
    {
      return 0.0;
    }

    /// calculating retention time of 10% of maximum (Apex)
    std::pair<double, double> rt_index_pair = computeBulkRetentionTimeRange(use_smoothed_ints);

    /// area-under-the-curve until 10% of maximum
    double peak_area(0.0);

    if (use_smoothed_ints)
    {
      auto smoothed_intensities = mass_trace_.getSmoothedIntensities();
      if (smoothed_intensities.size()==0) // if empty
      {
        return peak_area;
      }
      double int_before = smoothed_intensities[rt_index_pair.first];
      double rt_before = mass_trace_[rt_index_pair.first].getRT();
      // note: '<=' operator, since rt_index_pair are all inclusive!
      for (Size i = rt_index_pair.first + 1; i <= rt_index_pair.second; ++i)
      {
        peak_area += (int_before + smoothed_intensities[i])/2 * (mass_trace_[i].getRT() - rt_before);
        int_before = smoothed_intensities[i];
        rt_before = mass_trace_[i].getRT();
      }
    }
    else
    {
      double int_before = mass_trace_[rt_index_pair.first].getIntensity();
      double rt_before = mass_trace_[rt_index_pair.first].getRT();
      // note: '<=' operator, since rt_index_pair are all inclusive!
      for (Size i = rt_index_pair.first + 1; i <= rt_index_pair.second; ++i)
      {
        peak_area += (int_before + mass_trace_[i].getIntensity())/2 * (mass_trace_[i].getRT() - rt_before);
        int_before = mass_trace_[i].getIntensity();
        rt_before = mass_trace_[i].getRT();
      }
    }

    return peak_area;
  }

  double FeatureSeed::getCentroidRT() const
  {
    return mass_trace_.getCentroidRT();
  }

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
  const std::vector<FeatureSeed>& FeatureGroup::getSeeds() const
  {
    return feature_seeds_;
  }

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

  std::vector<FeatureSeed> FeatureGroup::getTheoreticalShapes() const
  {
    return theoretical_shapes_;
  }

  void FeatureGroup::setMonoisotopicMass(const double mass)
  {
    monoisotopic_mass_ = mass;
  }

  void FeatureGroup::setMaxIsotopeIndex(const Size index)
  {
    max_isotope_index_ = index;
  }

  void FeatureGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
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

  void FeatureGroup::updateTheoreticalShapes(std::vector<FeatureSeed> const &shapes)
  {
    if (theoretical_shapes_.empty())
    {
      theoretical_shapes_.reserve(this->size());
    }
    theoretical_shapes_.insert(theoretical_shapes_.end(), shapes.begin(), shapes.end());
  }

  void FeatureGroup::updateMembers(bool use_smoothed_ints)
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
    Size max_peak_idx = most_abundant_seed->getMassTrace().findMaxByIntPeak(use_smoothed_ints);  // smoothed intensity
    auto tmp_max_peak = most_abundant_seed->getMassTrace()[max_peak_idx];
    centroid_rt_of_most_abundant_mt_ = tmp_max_peak.getRT();

    fwhm_range_ = std::make_pair(min_fwhm, max_fwhm);

    // for fast searching later
    std::sort(ltrace_indices_.begin(), ltrace_indices_.end());
  }

  void FeatureGroup::updateMembersForScoring()
  {
    // based on PeakGroup::updateMonoMassAndIsotopeIntensities()
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
      double fi = f.getIntensity();
      per_isotope_int_[f.getIsotopeIndex()] += fi;
      charges_.insert(f.getCharge());
      intensity_ += fi;
      nominator += fi * (f.getUnchargedMass() - f.getIsotopeIndex() * Constants::ISOTOPE_MASSDIFF_55K_U);
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
