// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MassTrace.h>

using namespace OpenMS;
namespace FLASHQuantHelper
{
  /**
   * @brief Internal structure to store MassTrace and its additional information
   * */
  class OPENMS_DLLAPI FeatureSeed
  {
  public:
    /// default constructor
    FeatureSeed() = default;

    /// default destructor
    ~FeatureSeed() = default;

    /// copy constructor
    FeatureSeed(const FeatureSeed &seed) = default;

    /// move constructor
    FeatureSeed(FeatureSeed &&other) = default;

    /// comparison operator (ascending order)
    bool operator<(const FeatureSeed &fs) const
    {
      return (centroid_mz_ < fs.centroid_mz_);
    }

    /// assignment operator
    FeatureSeed &operator=(const FeatureSeed &seed) = default;

    /// constructor from MassTrace
    FeatureSeed(const MassTrace &mt) :
        mass_trace_(mt),
        centroid_mz_(mt.getCentroidMZ()),
        charge_(-1),
        intensity_(mt.computePeakArea()),
        isotope_index_(-1),
        trace_index_(0), // index of current trace (out of all input mass traces), thus not set here but after this construction
        mass_(0) // determined mass after deconvolution. NOT monoisotopic but only decharged
      {
        auto fwhm = mt.getFWHMborders();
        fwhm_start_ = mt[fwhm.first].getRT();
        fwhm_end_ = mt[fwhm.second].getRT();
      }

    /// default getter
    const MassTrace& getMassTrace() const;
    double getCentroidMz() const;
    int getCharge() const;
    double getFwhmStart() const;
    double getFwhmEnd() const;
    double getIntensity() const;
    int getIsotopeIndex() const;
    Size getTraceIndex() const;
    double getMass() const;

    /// default setter
    void setMassTrace(MassTrace &mt);
    void setCentroidMz(double &mz);
    void setCharge(int cs);
    void setFwhmStart(double fwhm_s);
    void setFwhmEnd(double fwhm_e);
    void setIntensity(double inty);
    void setIsotopeIndex(int idx);
    void setTraceIndex(Size i);
    void setMass(double mass);

    /// calculating and setting uncharged mass
    double getUnchargedMass();

    /// calculating retention time of 10% maximum (Apex) and area-under-the-curve until that point  (for FeatureGroupQuantity)
    std::pair<Size, Size> computeBulkRetentionTimeRange() const;
    double computeBulkPeakArea() const;

  private:
    MassTrace mass_trace_;

    double centroid_mz_; // centroid mz from mass trace
    int charge_;
    double fwhm_start_; // from mass trace, in sec
    double fwhm_end_;
    double intensity_;
    int isotope_index_;
    Size trace_index_;
    // determined mass after deconvolution. NOT monoisotopic but only decharged
    double mass_;
  };

  /**
   * @brief vector class for mass traces from same molecule, different charges and isotope indices
   * */
  class OPENMS_DLLAPI FeatureGroup
  {
  public :
    /// default constructor
    FeatureGroup() = default;

    /// default destructor
    ~FeatureGroup() = default;

    /// copy constructor
    FeatureGroup(const FeatureGroup &) = default;

    /// move constructor
    FeatureGroup(FeatureGroup &&other) = default;

    /// comparison operators (using monoisotopic_mass_)
    bool operator<(const FeatureGroup &a) const;
    bool operator>(const FeatureGroup &a) const;
    bool operator==(const FeatureGroup &a) const;

    /// assignment operator
    FeatureGroup &operator=(const FeatureGroup &t) = default;

    /// constructor with PeakGroup
    FeatureGroup(PeakGroup &pgroup)
    {
      monoisotopic_mass_ = pgroup.getMonoMass();
      auto cs_range = pgroup.getAbsChargeRange();
      min_abs_charge_ =  std::get<0>(cs_range);
      max_abs_charge_ = std::get<1>(cs_range);
      intensity_ = pgroup.getIntensity();
      isotope_cosine_score_ = pgroup.getIsotopeCosine();
    }

    /// explicit constructor for lower_bound / upper_bound search
    explicit FeatureGroup(const double mass) :
        monoisotopic_mass_(mass), intensity_(0)
    {
    }

    /** default getters **/
    double getMonoisotopicMass() const;
    int getMinCharge() const;
    int getMaxCharge() const;
    Size getMaxIsotopeIndex() const;
    double getIntensity() const;
    double getRtOfMostAbundantMT() const;
    float getIsotopeCosine() const;
    float getFeatureGroupScore() const;

    const std::set<int> &getChargeSet() const;
    const std::pair<double, double>& getFwhmRange() const;
    const std::vector<Size>& getTraceIndices() const;
    const std::vector<float>& getIsotopeIntensities() const;
    const std::vector<float>& getChargeIntensities() const;
    float getIntensityOfCharge(const int &abs_charge) const;
    float getIsotopeCosineOfCharge(const int &abs_charge) const;
    double getAverageMass() const;

    /** default setters **/
    void setMonoisotopicMass(const double mass);
    void setChargeRange(const int min_c, const int max_c);
    void setMaxIsotopeIndex(const Size index);
    void setIsotopeCosine(const float cos);
    void setFeatureGroupScore(const float score);
    void setPerChargeIntensities(std::vector<float> const &perChargeInt);
    void setPerChargeCosineScore(std::vector<float> const &perChargeCos);
    void setAverageMass(double averageMass);

    /// update multiple variables in one function
    void updateMembers(); // update after feature_seeds_ is changed
    void updateMembersForScoring(); // update primitively for scoring
    void updateIsotopeIndices(const int offset);

    /// checking the information within FeatureGroup
    bool doesThisChargeExist(int charge) const;
    FeatureSeed* getApexLMT() const;

    /// iterator settings
    std::vector<FeatureSeed>::const_iterator begin() const noexcept;
    std::vector<FeatureSeed>::const_iterator end() const noexcept;
    std::vector<FeatureSeed>::iterator begin() noexcept;
    std::vector<FeatureSeed>::iterator end() noexcept;

    const FeatureSeed& operator[](const Size i) const;

    void push_back(const FeatureSeed& seed);
    Size size() const noexcept;
    void reserve (Size n);
    void clear();
    std::vector<FeatureSeed>::iterator erase(std::vector<FeatureSeed>::iterator pos);
    bool empty() const;
    void swap (std::vector<FeatureSeed>& seed);
    void sort();

  private:
    /// features to be grouped
    std::vector<FeatureSeed> feature_seeds_;

    /// information of the deconvolved mass
    double monoisotopic_mass_;
    /// charge range
    int min_abs_charge_, max_abs_charge_; // absolute charge states.
    /// the largest isotope index
    int max_isotope_index_;
    /// summed intensities of feature_seeds_
    double intensity_;
    /// RT value from most abundant MassTrace.
    double centroid_rt_of_most_abundant_mt_;

    /// scores
    float isotope_cosine_score_;
    float total_score_;

    /// list of charges from feature_seeds
    std::set<int> charges_;
    /// min(fwhm_start) and max(fwhm_start) of MassTraces
    std::pair<double, double> fwhm_range_;
    /// Index to MassTraces included in feature_seeds
    std::vector<Size> ltrace_indices_;

    /// intensities per isotope index
    std::vector<float> per_isotope_int_;

    /// variables for writing results only
    std::vector<float> per_charge_int_;
    std::vector<float> per_charge_cos_;
    double average_mass_;
  };

  class OPENMS_DLLAPI CmpFeatureSeedByRT
  {
  public:
    bool operator()(const FeatureSeed& x, const FeatureSeed& y) const
    {
      if(x.getFwhmStart() == y.getFwhmStart()){
        return x.getFwhmEnd() < y.getFwhmEnd();
      }
      return x.getFwhmStart() < y.getFwhmStart();
    }
  };

  class OPENMS_DLLAPI CmpFeatureSeedByMZ
  {
  public:
    bool operator()(const FeatureSeed* x, const FeatureSeed* y) const
    {
      return x->getCentroidMz() < y->getCentroidMz();
    }
  };

  class OPENMS_DLLAPI CmpFeatureSeedByIntensity
  {
  public:
    bool operator()(const FeatureSeed* x, const FeatureSeed* y) const
    {
      // descending order
      return x->getIntensity() > y->getIntensity();
    }
  };

  class OPENMS_DLLAPI CmpFeatureGroupByScore
  {
  public:
    bool operator()(const FeatureGroup* x, const FeatureGroup* y) const
    {
      // intensity
      if(x->getIntensity() == y->getIntensity()){
        return x->getIsotopeCosine() < y->getIsotopeCosine();
      }

      return x->getIntensity() < y->getIntensity();
    }
  };

  class OPENMS_DLLAPI CmpFeatureGroupPointersByMass
  {
  public:
    bool operator()(const FeatureGroup* x, const FeatureGroup* y) const
    {
      return x->getMonoisotopicMass() < y->getMonoisotopicMass();
    }
  };

  struct OPENMS_DLLAPI Feature
  {
  public:
    std::vector<FeatureSeed> unique_traces;
    std::vector<FeatureSeed> shared_traces;
    std::vector<Size> unique_trace_indices; // index to input shared_m_traces_indices
    std::vector<Size> shared_trace_indices; // index to input shared_m_traces_indices
    std::vector<double> isotope_probabilities; // used as weights to EGHTraceFitter. Index of this vec = same index as unique's
//    LogMassTrace* most_abundant_mt_in_fg = nullptr;

    int charge;
    Size feature_group_index;

    /// default constructor
    Feature() = default;

    /// default destructor
    ~Feature() = default;

    Size getPeakSizes() const
    {
      Size total_peaks_size = 0;
      for (auto &lmt : unique_traces)
      {
        total_peaks_size += lmt.getMassTrace().getSize();
      }
      return total_peaks_size;
    }

    void prepareVectors(const Size n)
    {
      unique_traces.reserve(n);
      shared_traces.reserve(n);
      unique_trace_indices.reserve(n);
      shared_trace_indices.reserve(n);
    }

    void shrinkVectors()
    {
      unique_traces.shrink_to_fit();
      shared_traces.shrink_to_fit();
      unique_trace_indices.shrink_to_fit();
      shared_trace_indices.shrink_to_fit();
    }
  };
}