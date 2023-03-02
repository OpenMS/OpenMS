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

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MassTrace.h>

using namespace OpenMS;
namespace FLASHDeconvQuantHelper
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
        mass_(0), // determined mass after deconvolution. NOT monoisotopic but only decharged
        trace_index_(0) // index of current trace (out of all input mass traces), thus not set here but after this construction
      {
        auto fwhm = mt.getFWHMborders();
        fwhm_start_ = mt[fwhm.first].getRT();
        fwhm_end_ = mt[fwhm.second].getRT();
      }

    /// getter & setter
    const MassTrace& getMassTrace() const
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

    void setMassTrace(MassTrace &mt)
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

    void setIsotopeIndex(int idx)
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
    MassTrace mass_trace_;

    double centroid_mz_; // centroid mz from mass trace
    int charge_;
    double fwhm_start_; // from mass trace, in sec
    double fwhm_end_;
    double intensity_;
    int isotope_index_;
    // determined mass after deconvolution. NOT monoisotopic but only decharged
    double mass_;
    Size trace_index_;
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