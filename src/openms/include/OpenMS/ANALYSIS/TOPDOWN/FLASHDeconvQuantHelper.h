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
    FeatureSeed(MassTrace &mt)
    {
      mass_trace_ = &mt;
      centroid_mz_ = mt.getCentroidMZ();
      auto fwhm = mt.getFWHMborders();
      fwhm_start_ = mt[fwhm.first].getRT();
      fwhm_end_ = mt[fwhm.second].getRT();
      intensity_ = mt.computePeakArea();

      /// initialize with -1
      charge_ = -1;
      isotope_index_ = -1;

      /// determined mass after deconvolution. NOT monoisotopic but only decharged
      mass_ = 0;

      /// index of current trace (out of all input mass traces), thus not set here but after this construction
      trace_index_ = 0;
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
      charge_score_ = pgroup.getChargeScore();
      isotope_cosine_score_ = pgroup.getIsotopeCosine();
      intensity_ = pgroup.getIntensity();
    }

    /// explicit constructor for lower_bound / upper_bound search
    explicit FeatureGroup(const double mass) :
    monoisotopic_mass_(mass)
    {
      intensity_ = 0;
    }

    /** default getters **/
    double getMonoisotopicMass() const;
    int getMinCharge() const;
    int getMaxCharge() const;
    Size getMaxIsotopeIndex() const;
    double getIntensity() const;
    double getCentroidRtOfApices() const;
    double getRtOfMostAbundantMT() const;
    float getChargeScore() const;
    float getIsotopeCosine() const;
    float getFeatureGroupScore() const;

    std::set<int> getChargeSet() const;
    std::pair<double, double> getFwhmRange() const;
    std::vector<Size> getTraceIndices() const;
    std::vector<float> getIsotopeIntensities() const;
//    float getIsotopeCosineOfCharge(const int &abs_charge) const;
//    float getIntensityOfCharge(const int &abs_charge) const;

    /** default setters **/
    void setChargeRange(const int min_c, const int max_c);
    void setMaxIsotopeIndex(const Size index);
    void setChargeScore(const float score);
    void setIsotopeCosine(const float cos);
    void setFeatureGroupScore(const float score);

//    void setChargeIsotopeCosine(const int abs_charge, const float cos);
//    void setChargeIntensity(const int abs_charge, const float intensity);

    /// update multiple variables in one function
    void updateMembers(); // update after feature_seeds_ is changed
    void updateMembersForScoring(); // update primitively for scoring

    /// find the Apex seed with specific charge
    FeatureSeed* getApexLMTofCharge(int charge) const;
//    bool hasPerChargeVector() const;

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
    /// centroid RT value among apice MT from each charges
//    double centroid_rt_of_apices_;
    /// RT value from most abundant MassTrace.
    double centroid_rt_of_most_abundant_mt_;

    /// scores
    float charge_score_;
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
    /// per charge values
//    std::vector<float> per_charge_cos_;
//    std::vector<float> per_charge_int_;
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
    std::vector<double> isotope_probabilities; // used as weights to EGHTraceFitter
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
}