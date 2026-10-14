// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// These classes are minimal implementations of the interfaces defined in ITransition:
//  - IFeature
//  - IMRMFeature
//  - ITransitionGroup
//  - ISignalToNoise

namespace OpenMS
{
  /**
    @brief An implementation of the OpenSWATH Feature Access interface using OpenMS

  */
  class OPENMS_DLLAPI FeatureOpenMS :
    public OpenSwath::IFeature
  {
public:

    explicit FeatureOpenMS(const Feature& feature);

    ~FeatureOpenMS() override;

    void getRT(std::vector<double>& rt) const override;

    void getIntensity(std::vector<double>& intens) const override;

    float getIntensity() const override;

    double getRT() const override;

private:
    const Feature* feature_;
  };

  /**
    @brief An implementation of the OpenSWATH MRM Feature Access interface using OpenMS

  */
  class OPENMS_DLLAPI MRMFeatureOpenMS :
    public OpenSwath::IMRMFeature
  {
public:

    /// Build access lanes from the native-id maps stored in @p mrmfeature.
    explicit MRMFeatureOpenMS(MRMFeature& mrmfeature);

    /// Build pointer-based access lanes from explicit fragment and precursor native-id lists.
    MRMFeatureOpenMS(MRMFeature& mrmfeature,
                     const std::vector<std::string>& feature_ids,
                     const std::vector<std::string>& precursor_feature_ids);

    /**
      @brief Build access lanes with externally aligned native-id order.

      The feature id vectors define the aligned order used for index-hint
      fast paths. Fragment features are read directly from the underlying
      MRMFeature storage, so the caller must ensure that @p feature_ids matches
      that storage order. The lookup-id vectors are only used to validate and
      initialize access to the underlying fragment and precursor features.
    */
    MRMFeatureOpenMS(MRMFeature& mrmfeature,
                     const std::vector<std::string>& feature_ids,
                     const std::vector<std::string>& precursor_feature_ids,
                     const std::vector<std::string>& feature_lookup_ids,
                     const std::vector<std::string>& precursor_feature_lookup_ids);

    ~MRMFeatureOpenMS() override;

    std::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) override;

    std::shared_ptr<OpenSwath::IFeature> getPrecursorFeature(std::string nativeID) override;

    std::vector<std::string> getNativeIDs() const override;

    std::vector<std::string> getPrecursorIDs() const override;

    float getIntensity() const override;

    double getRT() const override;

    double getMetaValue(std::string name) const;

    size_t size() const override;

    void getFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const;

    void getPrecursorFeatureIntensities(const std::vector<std::string>& native_ids, std::vector<std::vector<double>>& intensities) const;

    /**
      @brief Return the fragment-feature intensity for @p native_id.

      The @p expected_index parameter is used as an optional fast-path hint.
      It is only trusted when this object was constructed with externally
      aligned feature ids. The default MRMFeature native-id order is map-based
      and must not be assumed to match the underlying feature storage order.
    */
    float getFeatureIntensity(const std::string& native_id, Size expected_index) const;

private:
    /**
      @brief Internal access lane for one set of features stored in an MRMFeature peak group.

      A lane encapsulates the lookup and wrapper state for either:

      - fragment-transition features, i.e. the per-transition peak features
        stored in MRMFeature and addressed by transition native ids, or
      - precursor features, i.e. the precursor/isotope peak features stored
        alongside the fragment-transition features and addressed by precursor
        native ids.

      In other words, one MRMFeature peak group owns multiple individual
      Feature objects, and a lane manages access to one of those two groups of
      child features. The implementation lives in the .cpp to keep this header
      compact, but conceptually each lane owns:

      - optional aligned native-id storage used for index-hint fast paths
      - either direct access to Feature storage or pointer-based lookup access
      - native-id to index lookup state
      - lazily materialized FeatureOpenMS wrappers

      Fragment-transition and precursor-feature handling share the same
      mechanics, so MRMFeatureOpenMS keeps one lane for fragment features and
      one lane for precursor features instead of duplicating all lookup and
      caching members in the outer class.
    */
    struct FeatureLane_;

    const MRMFeature& mrmfeature_;
    /// Fragment-feature access lane. May use aligned ids when the caller provides storage-aligned feature ids.
    std::unique_ptr<FeatureLane_> feature_lane_;
    /// Precursor-feature access lane for precursor chromatogram features.
    std::unique_ptr<FeatureLane_> precursor_feature_lane_;
  };

  /**
    @brief An implementation of the OpenSWATH Transition Group Access interface using OpenMS

  */
  template <typename SpectrumT, typename TransitionT>
  class TransitionGroupOpenMS :
    public OpenSwath::ITransitionGroup
  {
public:

    TransitionGroupOpenMS(MRMTransitionGroup<SpectrumT, TransitionT>& trgroup) :
      trgroup_(trgroup)
    {
    }

    ~TransitionGroupOpenMS() override
    {
    }

    std::size_t size() const override
    {
      return trgroup_.size();
    }

    std::vector<std::string> getNativeIDs() const override
    {
      std::vector<std::string> result;
      for (std::size_t i = 0; i < this->size(); i++)
      {
        result.push_back(trgroup_.getChromatograms()[i].getNativeID());
      }
      return result;
    }

    void getLibraryIntensities(std::vector<double>& intensities) const override
    {
      trgroup_.getLibraryIntensity(intensities);
    }

private:
    const MRMTransitionGroup<SpectrumT, TransitionT>& trgroup_;
  };

  /**
    @brief An implementation of the OpenSWATH SignalToNoise Access interface using OpenMS

  */
  template <typename ContainerT>
  class SignalToNoiseOpenMS :
    public OpenSwath::ISignalToNoise
  {
public:

    SignalToNoiseOpenMS(const ContainerT& chromat,
                        double sn_win_len_, unsigned int sn_bin_count_, bool write_log_messages) :
      chromatogram_(chromat)
    {
      // One estimator is needed per chromatogram, millions of them in a DIA run. Constructing one
      // builds and validates its whole parameter set, which costs more than the estimate itself:
      // take a used one from this thread's pool instead. init() below resets everything an
      // estimator keeps from the previous chromatogram.
      std::vector<std::unique_ptr<Estimator_>>* pool = pool_();
      if (pool != nullptr && !pool->empty())
      {
        estimator_ = std::move(pool->back());
        pool->pop_back();
      }
      else
      {
        estimator_ = std::make_unique<Estimator_>();
      }
      if (!estimator_->configured || estimator_->win_len != sn_win_len_ ||
          estimator_->bin_count != sn_bin_count_ || estimator_->write_log_messages != write_log_messages)
      {
        OpenMS::Param snt_parameters = estimator_->sn.getParameters();
        snt_parameters.setValue("win_len", sn_win_len_);
        snt_parameters.setValue("bin_count", sn_bin_count_);

        if (write_log_messages)
        {
          snt_parameters.setValue("write_log_messages", "true");
        }
        else
        {
          snt_parameters.setValue("write_log_messages", "false");
        }

        estimator_->sn.setParameters(snt_parameters);
        estimator_->win_len = sn_win_len_;
        estimator_->bin_count = sn_bin_count_;
        estimator_->write_log_messages = write_log_messages;
        estimator_->configured = true;
      }
      estimator_->sn.init(chromatogram_);
    }

    ~SignalToNoiseOpenMS() override
    {
      std::vector<std::unique_ptr<Estimator_>>* pool = pool_();
      if (pool != nullptr && estimator_ && pool->size() < max_pooled_estimators_)
      {
        pool->push_back(std::move(estimator_));
      }
    }

    SignalToNoiseOpenMS(const SignalToNoiseOpenMS&) = delete;
    SignalToNoiseOpenMS& operator=(const SignalToNoiseOpenMS&) = delete;

    double getValueAtRT(double RT) override
    {
      if (chromatogram_.empty()) {return -1;}

      // Note that MZBegin does not seem to return the same iterator on
      // different setups, see https://github.com/OpenMS/OpenMS/issues/1163
      auto iter = chromatogram_.PosEnd(RT);

      // ensure that iter is valid
      if (iter == chromatogram_.end()) 
      {
        iter--;
      }

      auto prev = iter;
      if (prev != chromatogram_.begin()) 
      {
        prev--;
      }

      if (std::fabs(prev->getPos() - RT) < std::fabs(iter->getPos() - RT) )
      {
        // prev is closer to the apex
        return estimator_->sn.getSignalToNoise((Size) distance(chromatogram_.begin(),prev));
      }
      else
      {
        // iter is closer to the apex
        return estimator_->sn.getSignalToNoise((Size) distance(chromatogram_.begin(),iter));
      }
    }

private:

    /// An estimator together with the parameters it was configured with.
    struct Estimator_
    {
      OpenMS::SignalToNoiseEstimatorMedian< ContainerT > sn;
      double win_len = 0.0;
      unsigned int bin_count = 0;
      bool write_log_messages = false;
      bool configured = false;
    };

    /// The pool of a thread; clears the thread's flag when the thread ends.
    struct Pool_
    {
      std::vector<std::unique_ptr<Estimator_>> estimators;
      ~Pool_() { poolAlive_() = false; }
    };

    /// Whether the pool of this thread still exists. Objects with static storage duration (there
    /// are some in MRMFeatureAccessOpenMS.cpp) are destroyed after the thread-local pool of the
    /// main thread; the flag is a plain bool and stays readable until then.
    static bool& poolAlive_()
    {
      thread_local bool alive = true;
      return alive;
    }

    /// Estimators that are not in use, per thread; nullptr once the thread has destroyed its pool.
    static std::vector<std::unique_ptr<Estimator_>>* pool_()
    {
      if (!poolAlive_())
      {
        return nullptr;
      }
      thread_local Pool_ pool;
      return &pool.estimators;
    }

    static constexpr Size max_pooled_estimators_ = 64;

    const ContainerT& chromatogram_;
    std::unique_ptr<Estimator_> estimator_;

  };

}
