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

using namespace std;
namespace OpenMS
{
  /**
   * FeatureFindingIntact: quantification algorithm for intact proteins
   * */
  class OPENMS_DLLAPI FeatureFindingIntact :
      public ProgressLogger,
      public DefaultParamHandler
  {

  public:
    /// Default constructor
    FeatureFindingIntact();

    /// Default destructor
    ~FeatureFindingIntact();

    /**
      @brief Internal structure used in @ref FeatureFindingIntact that keeps
      track of a feature hypothesis (isotope group hypothesis).

      @ingroup Quantitation
    */
    class OPENMS_DLLAPI FeatureHypothesis
    {
    public:
      /// default constructor
      FeatureHypothesis():
        charge_(),
        feat_score_(),
        feature_mass_(),
        iso_pattern_traces_(),
        iso_mt_index_pairs_(),
        charge_score_(),
        rt_start_(),
        rt_end_(),
        fwhm_range_(),
        mz_start_(),
        mz_end_(),
        mz_score_(),
        rt_score_(),
        inty_score_(),
        scores_per_mt_(),
        quant_value_()
      {
      }

      /// default destructor
      ~FeatureHypothesis()
      {}

      /// copy constructor
      FeatureHypothesis(const FeatureHypothesis& fh):
          charge_(fh.charge_),
          feat_score_(fh.feat_score_),
          feature_mass_(fh.feature_mass_),
          iso_pattern_traces_(fh.iso_pattern_traces_),
          iso_mt_index_pairs_(fh.iso_mt_index_pairs_),
          charge_score_(fh.charge_score_),
          rt_start_(fh.rt_start_),
          rt_end_(fh.rt_end_),
          fwhm_range_(fh.fwhm_range_),
          mz_start_(fh.mz_start_),
          mz_end_(fh.mz_end_),
          mz_score_(fh.mz_score_),
          rt_score_(fh.rt_score_),
          inty_score_(fh.inty_score_),
          scores_per_mt_(fh.scores_per_mt_),
          quant_value_(fh.quant_value_)
      {
      }

      /// assignment operator
      FeatureHypothesis& operator=(const FeatureHypothesis& fh)
      {
        if (this == &fh)
          return *this;

        charge_ = fh.charge_;
        feat_score_ = fh.feat_score_;
        feature_mass_ = fh.feature_mass_;
        iso_pattern_traces_ = fh.iso_pattern_traces_;
        iso_mt_index_pairs_ = fh.iso_mt_index_pairs_;
        charge_score_ = fh.charge_score_;
        rt_start_ = fh.rt_start_;
        rt_end_ = fh.rt_end_;
        fwhm_range_ = fh.fwhm_range_;
        mz_start_ = fh.mz_start_;
        mz_end_ = fh.mz_end_;
        mz_score_ = fh.mz_score_;
        rt_score_ = fh.rt_score_;
        inty_score_ = fh.inty_score_;
        scores_per_mt_ = fh.scores_per_mt_;
        quant_value_ = fh.quant_value_;
        return *this;
      }

      /// comparison operator (descending order)
      bool operator > (const FeatureHypothesis &fh) const
      {
        return (feat_score_ > fh.feat_score_);
      }

      /// comparison operator (ascending order)
      bool operator < (const FeatureHypothesis &fh) const
      {
        return (feat_score_ < fh.feat_score_);
      }

      /// getter & setter
      int getCharge() const
      {
        return charge_;
      }

      Size getSize() const
      {
        return iso_pattern_traces_.size();
      }

      double getScore() const
      {
        return feat_score_;
      }

      double getFeatureMass() const
      {
        return feature_mass_;
      }

      std::vector<std::pair<Size, Size>> getIndicesOfMassTraces() const
      {
        return iso_mt_index_pairs_;
      }

      double getChargeScore() const
      {
        return charge_score_;
      }

      double getMZScore() const
      {
        return mz_score_;
      }

      double getRTScore() const
      {
        return rt_score_;
      }

      double getIntyScore() const
      {
        return inty_score_;
      }

      double getQuant() const
      {
        return quant_value_;
      }

      std::pair<double, double> getRTRange() const
      {
        return std::make_pair(rt_start_, rt_end_);
      }

      std::pair<double, double> getMZRange() const
      {
        return std::make_pair(mz_start_, mz_end_);
      }

      std::pair<double, double> getFwhmRange() const
      {
        return fwhm_range_;
      }

      void setCharge(const int& ch)
      {
        charge_ = ch;
      }

      void setScore(const double& score)
      {
        feat_score_ = score;
      }

      void setFeatureMass(const double & mass)
      {
        feature_mass_ = mass;
      }

      void setIndicesOfMassTraces(const std::vector<std::pair<Size, Size>>& index_pairs)
      {
        iso_mt_index_pairs_ = index_pairs;
      }

      void setChargeScore(const double cscore)
      {
        charge_score_ = cscore;
      }

      void setMZScore(const double score)
      {
        mz_score_ = score;
      }

      void setRTScore(const double score)
      {
        rt_score_ = score;
      }

      void setIntyScore(const double score)
      {
        inty_score_ = score;
      }

      void setFwhmRange()
      {
        double min_fwhm(numeric_limits<double>::max());
        double max_fwhm(0.0);
        for (auto& iso : iso_pattern_traces_)
        {
          std::pair<Size, Size> fwhm_idx(iso->getFWHMborders());
          std::pair<double, double> tmp_fwhm( (*iso)[fwhm_idx.first].getRT(), (*iso)[fwhm_idx.second].getRT());

          if (tmp_fwhm.first < min_fwhm)
            min_fwhm = tmp_fwhm.first;
          if (tmp_fwhm.second > max_fwhm)
            max_fwhm = tmp_fwhm.second;
        }
        fwhm_range_ = std::make_pair(min_fwhm, max_fwhm);
      }

      vector<double> getAllIntensities(bool smoothed = false) const
      {
        vector<double> tmp;
        for (Size i = 0; i < iso_pattern_traces_.size(); ++i)
        {
          tmp.push_back(iso_pattern_traces_[i]->getIntensity(smoothed));
        }
        return tmp;
      }

      vector<const MassTrace*> getMassTraces() const
      {
        return iso_pattern_traces_;
      }

      void removeMassTrace(const Size index)
      {
        iso_pattern_traces_.erase(iso_pattern_traces_.begin() + index);
        scores_per_mt_.erase(scores_per_mt_.begin() + index);
        iso_mt_index_pairs_.erase(iso_mt_index_pairs_.begin() + index);
      }

      /// adding mass trace
      void addMassTrace(const MassTrace& mt_ptr)
      {
        iso_pattern_traces_.push_back(&mt_ptr);
      }

      // adding scores per mass traces
      void addMassTraceScore(const double& mt_score)
      {
        scores_per_mt_.push_back(mt_score);
      }

      void updateFeatureMass()
      {
        // average mass
        double mono_mz = iso_pattern_traces_[0]->getCentroidMZ()
            - (iso_mt_index_pairs_[0].first * Constants::C13C12_MASSDIFF_U / charge_);
        feature_mass_ = (mono_mz - Constants::PROTON_MASS_U) * charge_;
      }

      void setRTRange()
      {
        double rt_lower_limit = std::numeric_limits<double>::max();
        double rt_upper_limit = 0.0;

        for (const auto& mt : iso_pattern_traces_)
        {
          const auto& bounding_box = mt->getConvexhull().getBoundingBox();

          if(bounding_box.minX() < rt_lower_limit)
          {
            rt_lower_limit = bounding_box.minX();
          }
          if(bounding_box.maxX() > rt_upper_limit)
          {
            rt_upper_limit = bounding_box.maxX();
          }
        }
        rt_start_ = rt_lower_limit;
        rt_end_ = rt_upper_limit;
      }

      void setMZRange(double mz_range)
      {
        mz_start_ = feature_mass_/charge_ + Constants::PROTON_MASS_U;
        mz_end_ = mz_start_ + mz_range;
      }

      void reset_feature_score()
      {
        double score(0.0);
        for (const auto& s : scores_per_mt_)
        {
          score += s;
        }
        feat_score_ = score;
      }

      double computeQuant()
      {
        // from mass traces
        double q(0.0);
        for (auto& mt : iso_pattern_traces_)
        {
          q += mt->computePeakArea();
        }
        quant_value_ = q;
        return q;
      }

    private:
      int charge_;
      double feat_score_;
      double feature_mass_;
      std::vector<const MassTrace*> iso_pattern_traces_;
      // first: iso index of current feature, second : masstrace index of final masstraces
      std::vector<std::pair<Size, Size>> iso_mt_index_pairs_;
      double charge_score_;
      // calculated based on iso_pattern_traces
      double rt_start_;
      double rt_end_;
      std::pair<double, double> fwhm_range_;
      double mz_start_;
      double mz_end_;
      double mz_score_;
      double rt_score_;
      double inty_score_;
      std::vector<double> scores_per_mt_;
      double quant_value_;
    };

    class OPENMS_DLLAPI CmpMassTraceByMZ
    {
    public:
      bool operator()(const MassTrace& x, const MassTrace& y) const
      {
        return x.getCentroidMZ() < y.getCentroidMZ();
      }
    };

    /// TODO: replace this with FLASHDeconv
    /// This struct contains the averagine patterns precalulated for speed up. Other variables are also calculated for fast cosine calculation
    struct OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes;
      /// L2 norms for masses
      std::vector<double> norms;
      /// mass differences between average mass and monoisotopic mass
//      std::vector<double> averageMassDelta;
      /// Isotope start indices: isotopes of the indices less than them have very low intensities
//      std::vector<Size> isotopeStartIndices;
      /// Isotope end indices: isotopes of the indices larger than them have very low intensities
//      std::vector<Size> isotopeEndIndices;
      /// max isotope index
      int maxIsotopeIndex;

      /// mass interval for calculation
      double massInterval;
      /// min mass for calculation
      double minMass;
    public:
      /// default constructor
      PrecalculatedAveragine() = default;

      /**
       @brief constructor with parameters such as mass ranges and interval ( delta ).
       @param m minMass
       @param M maxMass
       @param delta mass interval between m and M
       @param generator generator by which the calulation is done
    */
      PrecalculatedAveragine(double m,
                             double M,
                             double delta,
                             CoarseIsotopePatternGenerator *generator)
          :
          massInterval(delta), minMass(m)
      {
        int i = 0;
        while (true)
        {
          double a = i * massInterval;
          i++;
          if (a < m)
          {
            continue;
          }
          if (a > M)
          {
            break;
          }
          auto iso = generator->estimateFromPeptideWeight(a);
          auto factor = .01;
          iso.trimRight(factor * iso.getMostAbundant().getIntensity());

          double norm = .0;
          Size mostAbundantIndex = 0;
          double mostAbundantInt = 0;

          for (Size k = 0; k < iso.size(); k++)
          {
            norm += iso[k].getIntensity() * iso[k].getIntensity();
            if (mostAbundantInt >= iso[k].getIntensity())
            {
              continue;
            }
            mostAbundantInt = iso[k].getIntensity();
            mostAbundantIndex = k;
          }

//          Size leftIndex = mostAbundantIndex;
          for (Size k = 0; k <= mostAbundantIndex; k++)
          {
            if (iso[k].getIntensity() > mostAbundantInt * factor)
            {
              break;
            }
            norm -= iso[k].getIntensity() * iso[k].getIntensity();
//            leftIndex--;
            iso[k].setIntensity(0);
          }

//          Size rightIndex = iso.size() - 1 - mostAbundantIndex;
//          for (Size k = iso.size() - 1; k >= mostAbundantIndex; k--)
//          {
//            if (iso[k].getIntensity() > mostAbundantInt * factor)
//            {
//              break;
//            }
//            norm -= iso[k].getIntensity() * iso[k].getIntensity();
//            rightIndex--;
//            iso[k].setIntensity(0);
//          }

//          iso.renormalize(); // to keep isotopes?
          norms.push_back(norm);
          isotopes.push_back(iso);
        }
      }

      /// get distribution for input mass
      IsotopeDistribution get(double mass) const
      {
        Size i = (Size) (.5 + (mass - minMass) / massInterval);
        i = i >= isotopes.size() ? isotopes.size() - 1 : i;
        return isotopes[i];
      }

      /// get max isotope index
      int getMaxIsotopeIndex() const
      {
        return maxIsotopeIndex;
      }

      /// get max isotope index
      void setMaxIsotopeIndex(int index)
      {
        maxIsotopeIndex = index;
      }

      /// get norm
      double getNorm(double mass) const
      {
        Size i = (Size) (.5 + (mass - minMass) / massInterval);
        i = i >= isotopes.size() ? isotopes.size() - 1 : i;
        return norms[i];
      }

      /// get isotope start index
//      Size getIsotopeStartIndex(double mass) const
//      {
//        Size i = (Size) (.5 + (mass - minMass) / massInterval);
//        i = i >= isotopes.size() ? isotopes.size() - 1 : i;
//        return isotopeStartIndices[i];
//      }

      /// get isotope end index
//      Size getIsotopeEndIndex(double mass) const
//      {
//        Size i = (Size) (.5 + (mass - minMass) / massInterval);
//        i = i >= isotopes.size() ? isotopes.size() - 1 : i;
//        return isotopeEndIndices[i];
//      }

      /// get mass difference between avg and mono masses
//      double getAverageMassDelta(double mass) const
//      {
//        Size i = (Size) (.5 + (mass - minMass) / massInterval);
//        i = i >= isotopes.size() ? isotopes.size() - 1 : i;
//        return averageMassDelta[i];
//      }

    };

    // main method of FeatureFindingIntact
    void run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap, String in_file_path);

    struct OPENMS_DLLAPI DeconvMassStruct
    {
    public:
      /// default constructor
      DeconvMassStruct() = default;

      double deconv_mass = 0.0; // median mass of current
      std::vector<Size> feature_idx;
      std::vector<double> feature_masses; // sorted
      std::set<int> charges;
      std::pair<double, double> fwhm_border;
      double combined_score = 0.0;
      double quant_values = 0.0;

      void initialize(double mass, int cs, Size f_idx, std::pair<double, double> fwhm, double score)
      {
        deconv_mass = mass;
        charges.insert(cs);
        feature_idx.push_back(f_idx);
        feature_masses.push_back(mass);
        fwhm_border = fwhm;
        combined_score = score;
      }

      bool operator<(const DeconvMassStruct& other) const
      {
        return combined_score < other.combined_score;
      }

      bool operator==(const DeconvMassStruct& other) const
      {
        return combined_score == other.combined_score;
      }

      void addFeatureHypothesis(double mass, int cs, int f_idx, std::pair<double, double> fwhm, double s)
      {
        feature_masses.push_back(mass);
        charges.insert(cs);
        feature_idx.push_back(f_idx);
        updateFwhmBorder(fwhm);
        combined_score += s;
      }

      void removeFeatureHypothesis(double mass, double score)
      {
        combined_score -= score;

        auto iter = std::find(feature_masses.begin(), feature_masses.end(), mass);
        if (iter != feature_masses.end())
        {
          feature_masses.erase(iter);
        }
      }

      // update deconv_mass from calculating median of masses
      // return true if deconv_mass is changed - this step is needed to reduce unnecessary step in DeconvMassStruct set update
      bool updateDeconvMass()
      {
        if (feature_masses.size() == 1)
        {
          deconv_mass = feature_masses[0];
          return false;
        }
        std::sort(feature_masses.begin(), feature_masses.end());

        Size m_size = feature_masses.size();
        Size mid = static_cast<Size>(m_size / 2.0);
        double new_mass(0.0);
        if ((m_size % 2) == 0)
        {
          new_mass = (feature_masses[mid - 1] +  feature_masses[mid]) / 2;
        }
        else
        {
          new_mass = feature_masses[mid];
        }

        if (new_mass != deconv_mass)
        {
          deconv_mass = new_mass;
          return true;
        }
        return false;
      }

      bool hasContinuousCharges() const
      {
        // at least three charges should be continuously appearing
        for (auto c_itr = charges.cbegin(); c_itr != charges.cend(); ++c_itr)
        {
          auto next_itr = std::next(c_itr);
          auto next_itr2 = std::next(c_itr, 2);
          if (next_itr == charges.cend() || next_itr2 == charges.cend())
          {
            break;
          }
          if ( (*next_itr2)-(*c_itr)==2 )
          {
            return true;
          }
        }
        return false;
      }

      void updateFwhmBorder(std::pair<double, double> new_fwhm)
      {
        // it is guaranteed new_fwhm overlaps with current fwhm_border;
        if (new_fwhm.first < fwhm_border.first)
        {
          fwhm_border.first = new_fwhm.first;
        }
        if (new_fwhm.second > fwhm_border.second)
        {
          fwhm_border.second = new_fwhm.second;
        }
      }
    };

  protected:
    void updateMembers_() override;

  private:
    /// method for builiding Feature Hypotheses
    void buildFeatureHypotheses_(std::vector<MassTrace>& input_mtraces,
                                 std::vector<FeatureHypothesis>& output_hypotheses,
                                 std::vector<std::vector<Size>>& shared_m_traces_indices,
                                 std::map<double, DeconvMassStruct>& deconv_masses) const;

    void findLocalFeatures_(const std::vector<std::pair<const MassTrace*, Size>>& candidates,
                            const double& total_intensity,
                            std::vector<FeatureHypothesis>& output_hypotheses,
                            std::map<double, DeconvMassStruct>& deconv_masses) const;

    double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const;

    double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, int charge) const;

    bool doFWHMbordersOverlap(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    // TODO: based on getCosine from FLASHDeconvAlgorithm
    double computeCosineSimOfDiffSizedVector_(const std::vector<double>& a,
                                      const IsotopeDistribution& b,
                                      const int& b_size,
                                      const double& b_norm,
                                      const int offset) const;

    double computeAveragineCosineSimScore_(const std::vector<double>& hypo_ints,
                                           const IsotopeDistribution& iso_dist,
                                           const Size& iso_size,
                                           const double& iso_norm,
                                           int& offset) const;

    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    void setAveragineModel_();

    void removeMassArtifacts_(const std::vector<FeatureHypothesis>& feat_hypotheses,
                              std::map<double, DeconvMassStruct>& deconv_masses) const;

    void setChargeScoreForFeatureHypothesis(std::vector<FeatureHypothesis>& candidate_hypotheses,
                                            std::vector<std::pair<double, int>>& feat_and_charges) const;

    void clusterFeatureHypotheses_(std::vector<FeatureHypothesis>& output_hypotheses,
                                   const std::vector<std::vector<Size>>& shared_m_traces_indices,
                                   std::map<double, DeconvMassStruct>& deconv_masses) const;

    void resolveConflictInCluster_(const std::vector<FeatureHypothesis>& feat_hypo,
                                   const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   const std::set<Size>& hypo_indices,
                                   std::vector<FeatureHypothesis>& out_features,
                                   std::vector<Size>& out_feature_idx,
                                   ofstream& outs,
                                   String& cluster_name) const;

    void addFeature2DeconvMassStruct(FeatureHypothesis &in_feature,
                                     Size feature_idx,
                                     std::map<double, DeconvMassStruct> &deconv_masses) const;

    void filterDeconvMassStruct(std::map<double, DeconvMassStruct> &deconv_masses,
                                const vector<FeatureHypothesis> &feat_hypotheses,
                                std::map<double, DeconvMassStruct>::iterator &curr_it,
                                const bool struct_needs_update = false ) const;

    /// parameters
    double local_rt_range_;
    double local_mz_range_;
    int charge_lower_bound_;
    int charge_upper_bound_;
    double min_mass_;
    double max_mass_;
    bool use_smoothed_intensities_;
    const double mass_tolerance_ = 1.5; // Da, for feature mass collection

    Size max_nr_traces_; // calculated from iso_model_ (FeatureFindingIntact::setAveragineModel())
    PrecalculatedAveragine iso_model_;
  };
}
