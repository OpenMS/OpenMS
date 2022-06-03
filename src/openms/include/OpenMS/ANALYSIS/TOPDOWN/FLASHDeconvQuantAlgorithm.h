// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <algorithm>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include "boost/dynamic_bitset.hpp"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <fstream>

using namespace std;
namespace OpenMS
{
  /**
   * @brief Internal structure to store MassTrace and its additional information
   */
  class OPENMS_DLLAPI FeatureSeed
  {
  public:
    /// default constructor
    FeatureSeed() :
        mass_trace_(),
        centroid_mz_(),
        charge_(),
        fwhm_start_(),
        fwhm_end_(),
        intensity_(),
        isotope_index_(),
        mass_(),
        trace_index_()
    {
    }

    /// default destructor
    ~FeatureSeed()
    {
    }

    /// copy constructor
    FeatureSeed(const FeatureSeed &seed)
    {
      mass_trace_ = seed.mass_trace_;
      centroid_mz_ = seed.centroid_mz_;
      charge_ = seed.charge_;
      fwhm_start_ = seed.fwhm_start_;
      fwhm_end_ = seed.fwhm_end_;
      intensity_ = seed.intensity_;
      isotope_index_ = seed.isotope_index_;
      mass_ = seed.mass_;
      trace_index_ = seed.trace_index_;
    }

    /// assignment operator
    FeatureSeed &operator=(const FeatureSeed &seed)
    {
      if (this == &seed)
        return *this;

      mass_trace_ = seed.mass_trace_;
      centroid_mz_ = seed.centroid_mz_;
      charge_ = seed.charge_;
      fwhm_start_ = seed.fwhm_start_;
      fwhm_end_ = seed.fwhm_end_;
      intensity_ = seed.intensity_;
      isotope_index_ = seed.isotope_index_;
      mass_ = seed.mass_;
      trace_index_ = seed.trace_index_;

      return *this;
    }

    /// constructor with mass trace
    FeatureSeed(MassTrace &mt)
    {
      mass_trace_ = &mt;
      centroid_mz_ = mt.getCentroidMZ();
      charge_ = 0;
      auto fwhm = mt.getFWHMborders();
      fwhm_start_ = mt[fwhm.first].getRT();
      fwhm_end_ = mt[fwhm.second].getRT();
      intensity_ = mt.computePeakArea();
      isotope_index_ = -1;

      /// deteremined mass after deconvolution. NOT monoisotopic but only decharged
      mass_ = 0;

      // index of current trace (out of all input mass traces)
      trace_index_ = 0;
    }

    /// comparison operator (ascending order)
    bool operator<(const FeatureSeed &fs) const
    {
      return (centroid_mz_ < fs.centroid_mz_);
    }

    /// getter & setter
    MassTrace* getMassTrace() const
    {
      return mass_trace_;
    }

    double getCentroidMz() const
    {
      return centroid_mz_;
    }

    int getCharge() const
    {
      return charge_;
    }

    double getFwhmStart() const
    {
      return fwhm_start_;
    }

    double getFwhmEnd() const
    {
      return fwhm_end_;
    }

    double getIntensity() const
    {
      return intensity_;
    }

    int getIsotopeIndex() const
    {
      return isotope_index_;
    }

    Size getTraceIndex() const
    {
      return trace_index_;
    }

    double getMass() const
    {
      return mass_;
    }

    void setMassTrace(MassTrace *mt)
    {
      mass_trace_ = mt;
    }

    void setCentroidMz(double &mz)
    {
      centroid_mz_ = mz;
    }

    void setCharge(int cs)
    {
      charge_ = cs;
    }

    void setFwhmStart(double fwhm_s)
    {
      fwhm_start_ = fwhm_s;
    }

    void setFwhmEnd(double fwhm_e)
    {
      fwhm_end_ = fwhm_e;
    }

    void setIntensity(double inty)
    {
      intensity_ = inty;
    }

    void setIsotopeIndex(double idx)
    {
      isotope_index_ = idx;
    }

    void setTraceIndex(Size i)
    {
      trace_index_ = i;
    }

    void setMass(double mass)
    {
      mass_ = mass;
    }

    double getUnchargedMass()
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

  private:
    MassTrace *mass_trace_;

    // log transformed centroid mz
    double centroid_mz_;
    int charge_;
    double fwhm_start_; // from mass trace, in sec
    double fwhm_end_;
    double intensity_;
    int isotope_index_;
    double mass_;
    Size trace_index_;
  };

  // vector class for mass traces from same molecule, different charges and isotope indices
  class OPENMS_DLLAPI FeatureGroup :
      private std::vector<FeatureSeed>
  {
  public :
    using std::vector<FeatureSeed>::push_back;
    using std::vector<FeatureSeed>::empty;
    using std::vector<FeatureSeed>::begin;
    using std::vector<FeatureSeed>::end;
    using std::vector<FeatureSeed>::size;
    using std::vector<FeatureSeed>::reserve;
    using std::vector<FeatureSeed>::swap;
    using std::vector<FeatureSeed>::shrink_to_fit;
    using std::vector<FeatureSeed>::erase;
    using std::vector<FeatureSeed>::operator[];

    /// default constructor
    FeatureGroup() = default;

    /// default destructor
    ~FeatureGroup() = default;

    /// comparison operators
    bool operator<(const FeatureGroup &a) const
    {
      if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
      {
        return this->intensity_ < a.intensity_;
      }
      return this->monoisotopic_mass_ < a.monoisotopic_mass_;
    }

    bool operator>(const FeatureGroup &a) const
    {
      if (this->monoisotopic_mass_ == a.monoisotopic_mass_)
      {
        return this->intensity_ > a.intensity_;
      }
      return this->monoisotopic_mass_ > a.monoisotopic_mass_;
    }

    bool operator==(const FeatureGroup &a) const
    {
      return
          this->monoisotopic_mass_ == a.monoisotopic_mass_
          && this->intensity_ == a.intensity_;
    }

    /// assignment operator
    FeatureGroup &operator=(const FeatureGroup &t) = default;

    /// copy constructor
    FeatureGroup(const FeatureGroup &) = default;

    /// move constructor
    FeatureGroup(FeatureGroup &&other) = default;

    /// constructor with DeconvolutedSpectrum
    FeatureGroup(PeakGroup &pgroup)
    {
      monoisotopic_mass_ = pgroup.getMonoMass();
      auto cs_range = pgroup.getAbsChargeRange();
      min_abs_charge_ =  std::get<0>(cs_range);
      max_abs_charge_ = std::get<1>(cs_range);
      charge_score_ = pgroup.getChargeScore();
      isotope_cosine_score_ = pgroup.getIsotopeCosine();
      intensity_ = pgroup.getIntensity();
    }

    // for lower_bound / upper_bound search
    FeatureGroup(const double mass) :
        monoisotopic_mass_(mass)
    {
      intensity_ = 0;
    }

    void updateMassesAndIntensity(const int offset = 0)
    {
      if (offset != 0)
      {
        std::vector<FeatureSeed> tmpPeaks;
        tmpPeaks.swap(*this);
        reserve(tmpPeaks.size());

        for (auto &p: tmpPeaks)
        {
          p.setIsotopeIndex(p.getIsotopeIndex() - offset);
          if (p.getIsotopeIndex() < 0 || p.getIsotopeIndex() >= max_isotope_index_)
          {
            continue;
          }
          push_back(p);
        }
      }

      intensity_ = .0;
      double nominator = .0;

      for (auto &p: *this)
      {
        double pi = p.getIntensity() + 1;
        intensity_ += pi;
        nominator += pi * (p.getUnchargedMass() - p.getIsotopeIndex() * Constants::ISOTOPE_MASSDIFF_55K_U);
      }
      monoisotopic_mass_ = nominator / intensity_;
    }

    double getMonoisotopicMass() const
    {
      return monoisotopic_mass_;
    }

    int getMinCharge() const
    {
      return min_abs_charge_;
    }

    int getMaxCharge() const
    {
      return max_abs_charge_;
    }

    std::vector<int> getChargeVector() const
    {
      return charges_;
    }

    Size getMaxIsotopeIndex() const
    {
      return max_isotope_index_;
    }

    double getIntensity() const
    {
      return intensity_;
    }

    float getIsotopeCosine() const
    {
      return isotope_cosine_score_;
    }

    float getChargeScore() const
    {
      return charge_score_;
    }

    float getFeatureGroupScore() const
    {
      return total_score_;
    }

    std::pair<double, double> getFwhmRange() const
    {
      return fwhm_range_;
    }

    std::vector<Size> getTraceIndices() const
    {
      return ltrace_indices_;
    }

    float getIsotopeCosineOfCharge(const int &abs_charge) const
    {
      return per_charge_cos_[abs_charge];
    }

    float getIntensityOfCharge(const int &abs_charge) const
    {
      return per_charge_int_[abs_charge];
    }

    bool hasPerChargeVector() const
    {
      if (per_charge_int_.empty())
      {
        return false;
      }
      return true;
    }

    double getCentroidRtOfApices() const
    {
      return centroid_rt_of_apices;
    }

    double getRtOfApex() const
    {
      return rt_of_apex;
    }

    FeatureSeed* getApexLMTofCharge(int charge) const
    {
      FeatureSeed* apex_lmt = nullptr;
      double max_intensity = 0.0;

      for (auto lmt_iter = this->begin(); lmt_iter != this->end(); ++lmt_iter)
      {
        if (lmt_iter->getCharge() == charge && lmt_iter->getIntensity() > max_intensity)
        {
          max_intensity = lmt_iter->getIntensity();
          apex_lmt = (FeatureSeed*) &(*lmt_iter);
        }
      }
      return apex_lmt;
    }

//    std::pair<double, double> getSummedIntensityOfMostAbundantMTperCS() const
//    {
//      auto per_cs_max_area = std::vector<double>(1 + max_abs_charge_, .0);
//      auto per_cs_max_inty = std::vector<double>(1 + max_abs_charge_, .0);
//
//      for (auto &lmt: *this)
//      {
//        double sum_inty = 0;
//        for (auto &p : *(lmt.getMassTrace()))
//        {
//          sum_inty += p.getIntensity();
//        }
//        if(per_cs_max_inty[lmt.getCharge()] < sum_inty)
//        {
//          per_cs_max_inty[lmt.getCharge()] = sum_inty;
//        }
//
//        if (per_cs_max_area[lmt.getCharge()] < lmt.getIntensity())
//        {
//          per_cs_max_area[lmt.getCharge()] = lmt.getIntensity();
//        }
//      }
//
//      return std::make_pair(std::accumulate(per_cs_max_area.begin(), per_cs_max_area.end(), .0),
//                            std::accumulate(per_cs_max_inty.begin(), per_cs_max_inty.end(), .0));
//    }

//    double getAvgFwhmLength() const
//    {
//      std::vector<double> fwhm_len_arr;
//      fwhm_len_arr.reserve(this->size());
//
//      for (auto &l_trace: *this)
//      {
//        double tmp_fwhm_len = l_trace.getFwhmEnd() - l_trace.getFwhmStart();
//        fwhm_len_arr.push_back(tmp_fwhm_len);
//      }
//      return accumulate(fwhm_len_arr.begin(), fwhm_len_arr.end(), 0.0) / (this->size());
//    }

    void setChargeRange(const int min_c, const int max_c)
    {
      min_abs_charge_ = min_c;
      max_abs_charge_ = max_c;
    }

    void setMaxIsotopeIndex(const Size index)
    {
      max_isotope_index_ = index;
    }

    void setChargeScore(const float score)
    {
      charge_score_ = score;
    }

    void setIsotopeCosine(const float cos)
    {
      isotope_cosine_score_ = cos;
    }

    void setFeatureGroupScore(const float score)
    {
      total_score_ = score;
    }

    void setChargeVector()
    {
      std::set<int> cs_set;
      for(auto &lmt : *this)
      {
        cs_set.insert(lmt.getCharge());
      }
      charges_ = std::vector<int>(cs_set.size());
      std::copy(cs_set.begin(), cs_set.end(), charges_.begin());
      std::sort(charges_.begin(), charges_.end());
    }

    void setChargeIsotopeCosine(const int abs_charge, const float cos)
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

    void setChargeIntensity(const int abs_charge, const float intensity)
    {
      if (max_abs_charge_ < abs_charge)
      {
        return;
      }
      if (per_charge_int_.empty())
      {
        per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);
      }
      per_charge_int_[abs_charge] = intensity;
    }

    void setFwhmRange()
    {
      double min_fwhm(numeric_limits<double>::max());
      double max_fwhm(0.0);

      for (auto &l_trace: *this)
      {
        std::pair<double, double> tmp_fwhm(l_trace.getFwhmStart(), l_trace.getFwhmEnd());

        if (tmp_fwhm.first < min_fwhm)
          min_fwhm = tmp_fwhm.first;
        if (tmp_fwhm.second > max_fwhm)
          max_fwhm = tmp_fwhm.second;
      }
      fwhm_range_ = std::make_pair(min_fwhm, max_fwhm);
    }

//    void initializePerChargeVectors()
//    {
//      per_charge_cos_.clear();
//      per_charge_int_.clear();
//      per_charge_cos_ = std::vector<float>(1 + max_abs_charge_, .0);
//      per_charge_int_ = std::vector<float>(1 + max_abs_charge_, .0);
//      // other vectors (per_charge_pwr_, per_charge_signal_pwr_, per_charge_snr_) don't need to be initialized
//    }

    void setTraceIndices()
    {
      ltrace_indices_.clear();
      ltrace_indices_.reserve(this->size());
      for (auto &l_trace: *this)
      {
        ltrace_indices_.push_back(l_trace.getTraceIndex());
      }
      std::sort(ltrace_indices_.begin(), ltrace_indices_.end());
    }

    void setCentroidRtOfApices()
    {
      double tmp_rt = .0;
      double max_inty = .0;
      double apex_rt = .0;

      for (const auto &lmt: *this)
      {
        Size max_idx = lmt.getMassTrace()->findMaxByIntPeak(false);
        auto tmp_max_peak = (*lmt.getMassTrace())[max_idx];
        tmp_rt += tmp_max_peak.getRT();
        if (max_inty < tmp_max_peak.getIntensity())
        {
          max_inty = tmp_max_peak.getIntensity();
          apex_rt = tmp_max_peak.getRT();
        }
      }
      tmp_rt /= this->size();
      centroid_rt_of_apices = tmp_rt;
      rt_of_apex = apex_rt;
    }


  private:
    /// information on the deconvouted mass
    double monoisotopic_mass_;
    /// charge range
    int min_abs_charge_, max_abs_charge_; // absolute charge states.
    int max_isotope_index_;
    double intensity_;
    float charge_score_;
    float isotope_cosine_score_;
    float total_score_;

    std::pair<double, double> fwhm_range_;
    std::vector<Size> ltrace_indices_;
    std::vector<int> charges_;

    std::vector<float> per_charge_cos_;
    std::vector<float> per_charge_int_;

    double centroid_rt_of_apices;
    double rt_of_apex; // Apex : apex of most abundant mass traces
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
    bool operator()(const FeatureGroup& x, const FeatureGroup& y) const
    {
      // intensity
      if(x.getIntensity() == y.getIntensity()){
        return x.getIsotopeCosine() < y.getIsotopeCosine();
      }

      return x.getIntensity() < y.getIntensity();
    }
  };

  struct OPENMS_DLLAPI FeatureElement
  {
  public:
    std::vector<FeatureSeed*> mass_traces;
    std::vector<Size> mass_trace_indices; // index to input shared_m_traces_indices
    std::vector<double> isotope_probabilities;
//    LogMassTrace* most_abundant_mt_in_fg = nullptr;

    int charge;
    Size feature_group_index;

    /// default constructor
    FeatureElement() = default;

    Size getPeakSizes() const
    {
      Size total_peaks_size = 0;
      for (auto &lmt : mass_traces)
      {
        total_peaks_size += lmt->getMassTrace()->getSize();
      }
      return total_peaks_size;
    }
   };

  class OPENMS_DLLAPI FLASHDeconvQuantAlgorithm :
      public ProgressLogger,
      public DefaultParamHandler
{
  public:
    /// Default constructor
    FLASHDeconvQuantAlgorithm();

    /// Default destructor
    ~FLASHDeconvQuantAlgorithm() override;

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace> &input_mtraces, FeatureMap &output_featmap);

    String outfile_path;

  protected:
    void updateMembers_() override;

  private:
    Param getFLASHDeconvParams_();

    // equivalent to FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_
    void getFeatureFromSpectrum_(std::vector<FeatureSeed*> &local_traces, std::vector<FeatureGroup> &local_fgroup, const double &rt);

    void buildMassTraceGroups_(std::vector<FeatureSeed> &in_seeds, std::vector<FeatureGroup> &features);

    bool scoreFeatureGroup_(FeatureGroup& fg) const;

    void calculatePerChargeIsotopeIntensity_(std::vector<double> &per_isotope_intensity,
                                                        std::vector<double> &per_charge_intensity,
                                                        FeatureGroup &fg) const;

    void refineFeatureGroups_(std::vector<FeatureGroup>& features);

    bool rescoreFeatureGroup_(FeatureGroup& fg, bool score_always = false) const;

    void setFeatureGroupScore_(FeatureGroup &fg) const;

    static double getCosine_(const std::vector<double> &a,
                             const int &a_start,
                             const int &a_end,
                             const IsotopeDistribution &b,
                             const int &b_size,
                             const int offset);

    double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const;

    double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const;

    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    bool doFWHMbordersOverlap(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    bool doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2) const;

    void clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups,
                               std::vector<MassTrace>& input_mtraces) const;

    void resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                   const std::vector<MassTrace> & input_masstraces,
                                   const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   const std::set<Size>& hypo_indices,
                                   std::vector<FeatureGroup>& out_features) const;

    void writeMassTracesOfFeatureGroup(const std::vector<FeatureGroup>& featgroups,
                                       const std::vector<std::vector<Size> >& shared_m_traces_indices) const;

    void writeFeatureGroupsInFile(std::vector<FeatureGroup>& feat, std::vector<MassTrace>& input_mtraces) const;

    void writeOutputInFeatureXML_(const std::vector<FeatureGroup> &feature_groups,
                                  const std::vector<std::vector<Size>> &shared_m_traces_indices) const;

    void storeFeatureGroupInOpenMSFeature(std::vector<FeatureGroup> &feature_groups,
                                          FeatureMap &out_featmap) const;

    void resolveConflictRegion_(std::vector<FeatureElement> &feat,
                                const std::vector<Size> &feature_idx_in_current_conflict_region,
                                const std::vector<const MassTrace*> &conflicting_mts) const;

    void runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces, EGHTraceFitter* fitter) const;

    void getMostAbundantMassTraceFromFeatureGroup(const FeatureGroup &fgroup,
                                                  const int &ignore_this_charge,
                                                  FeatureSeed* &most_abundant_mt_ptr,
                                                  const std::vector<std::vector<Size>>& shared_m_traces) const;

    void getFLASHDeconvConsensusResult();

    bool isThisMassOneOfTargets(const double &candi_mass, const double &candi_rt) const;

    void makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const;

    /// parameter stuff
//    double local_rt_range_;
    double local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    int charge_range_;
    double min_mass_;
    double max_mass_;
    double mz_tolerance_; // ppm

    double total_intensity_;

    const double mass_tolerance_da_ = 3; // Da, for feature mass collection
//    const double mass_tolerance_ppm_ = 20;

    // advanced parameter?
    Size min_nr_mtraces_ = 3; // minimum number of consecutive bridges among mass traces to support feature
    bool use_smoothed_intensities_;
    double rt_window_ = 1; // TODO : remove?

    /// variables for internal use (not for user input)
    FLASHDeconvHelperStructs::PrecalculatedAveragine iso_model_;
    Size max_nr_traces_; // calculated from iso_model_ (setAveragineModel())

    /// cosine threshold between observed and theoretical isotope patterns for MS1
    double min_isotope_cosine_ = 0.90;

    /// FLASHDeconvAlgorithm class for deconvolution
    FLASHDeconvAlgorithm fd_;

    // loop up table
    std::vector<std::pair<double, double>> target_masses_; // mass and rt
    bool with_target_masses_ = false;

  };
}