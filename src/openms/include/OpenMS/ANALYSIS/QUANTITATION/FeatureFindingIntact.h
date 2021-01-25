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

using namespace std;
namespace OpenMS
{
  /**
   * FeatureFindingIntact: quantification algorithm for intact proteins
   * */
  class OPENMS_DLLAPI FeatureFindingIntact :
      public ProgressLogger
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
        charge_score_()
      {
      }

      /// default destructor
      virtual ~FeatureHypothesis()
      {

      }

      /// copy constructor
      FeatureHypothesis(const FeatureHypothesis& fh):
          charge_(fh.charge_),
          feat_score_(fh.feat_score_),
          feature_mass_(fh.feature_mass_),
          iso_pattern_traces_(fh.iso_pattern_traces_),
          iso_mt_index_pairs_(fh.iso_mt_index_pairs_),
          charge_score_(fh.charge_score_)
      {
      }

      /// assignment operator
      FeatureHypothesis& operator=(const FeatureHypothesis& fh)
      {
        if (this == &fh)
          return *this;

        charge_ = fh.charge_;
        iso_pattern_traces_ = fh.iso_pattern_traces_;
        feat_score_ = fh.feat_score_;
        feature_mass_ = fh.feature_mass_;
        iso_mt_index_pairs_ = fh.iso_mt_index_pairs_;
        charge_score_ = fh.charge_score_;
        return *this;
      }

      /// getter & setter
      SignedSize getCharge() const
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

      void setCharge(const SignedSize& ch)
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

      vector<double> getAllIntensities(bool smoothed = false) const
      {
        vector<double> tmp;
        for (Size i = 0; i < iso_pattern_traces_.size(); ++i)
        {
          tmp.push_back(iso_pattern_traces_[i]->getIntensity(smoothed));
        }
        return tmp;
      }

      double getCentroidMZofMonoisotopicFeature() const
      {
        if (iso_pattern_traces_.empty())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "FeatureHypothesis is empty, no centroid MZ!",
                                        String(iso_pattern_traces_.size()));
        }
        return iso_pattern_traces_[0]->getCentroidMZ();
      }

      vector<const MassTrace*> getMassTraces() const
      {
        return iso_pattern_traces_;
      }

      /// adding mass trace
      void addMassTrace(const MassTrace& mt_ptr)
      {
        iso_pattern_traces_.push_back(&mt_ptr);
      }

      void updateFeatureMass()
      {
        double mono_mz = iso_pattern_traces_[0]->getCentroidMZ()
            - (iso_mt_index_pairs_[0].first * Constants::C13C12_MASSDIFF_U / charge_);
        feature_mass_ = (mono_mz - Constants::PROTON_MASS_U) * charge_;
      }

    private:
      SignedSize charge_;
      double feat_score_;
      double feature_mass_;
      std::vector<const MassTrace*> iso_pattern_traces_;
      // first: iso index of current feature, second : masstrace index of final masstraces
      std::vector<std::pair<Size, Size>> iso_mt_index_pairs_;
      double charge_score_;
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
    void run(std::vector<MassTrace>& input_mtraces, FeatureMap& output_featmap);

    // TODO : change this according to the defaultHandler;
    void updateMembers_()
    {
      local_mz_range_ = 6.5 ; // MZ range where to look for isotopic mass traces (-> decides size of isotopes =(local_mz_range_ * lowest_charge))
      local_rt_range_ = 15.0 ; // RT range where to look for coeluting mass traces
      charge_lower_bound_ = 7;
      charge_upper_bound_ = 30;
      mass_upper_bound = 20000;
      use_smoothed_intensities_ = true; // for intensity of a mass trace
    }

  private:
    /// method for builiding Feature Hypotheses
    void buildFeatureHypotheses_(std::vector<MassTrace>& input_mtraces,
                                 std::vector<FeatureHypothesis>& output_hypotheses,
                                 std::vector<std::vector<Size>>& shared_m_traces_indices) const;

    void findLocalFeatures_(const std::vector<std::pair<const MassTrace*, Size>>& candidates,
                            const double total_intensity,
                            std::vector<FeatureHypothesis>& output_hypotheses,
                            std::vector<std::vector<Size>>& shared_m_traces_indices) const;

    double scoreRT_(const MassTrace& tr1, const MassTrace& tr2) const;

    double scoreMZ_(const MassTrace& tr1, const MassTrace& tr2, Size iso_pos, Size charge) const;

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

    void setAveragineModel();

    void clusterFeatureHypotheses_(std::vector<FeatureHypothesis>& output_hypotheses,
                                  const std::vector<std::vector<Size>>& shared_m_traces_indices) const;

    void resolveConflictInCluster_(const std::vector<FeatureHypothesis>& feat_hypo,
                                   const std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   const std::set<Size>& hypo_indices,
                                   std::vector<FeatureHypothesis>& out_features
                                   ) const;

    /// parameters
    double local_rt_range_;
    double local_mz_range_;
    double charge_lower_bound_;
    double charge_upper_bound_;
    double mass_upper_bound;
    Size max_nr_traces_; // calculated from iso_model_ (FeatureFindingIntact::setAveragineModel())
    bool use_smoothed_intensities_;

    PrecalculatedAveragine iso_model_;
  };
}
