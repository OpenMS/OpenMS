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
#include <OpenMS/DATASTRUCTURES/String.h>

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

    explicit MRMFeatureOpenMS(MRMFeature& mrmfeature);

    MRMFeatureOpenMS(MRMFeature& mrmfeature,
                     const std::vector<std::string>& feature_ids,
                     const std::vector<std::string>& precursor_feature_ids);

    MRMFeatureOpenMS(MRMFeature& mrmfeature,
                     const std::vector<std::string>& feature_ids,
                     const std::vector<std::string>& precursor_feature_ids,
                     const std::vector<String>& feature_lookup_ids,
                     const std::vector<String>& precursor_feature_lookup_ids);

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

    float getFeatureIntensity(const std::string& native_id) const;

    float getFeatureIntensity(const std::string& native_id, std::size_t expected_index) const;

    float getFeatureIntensity(std::size_t index) const;

    bool hasCachedFeatureIds() const;

private:
    using FeatureWrapperList_ = std::vector<FeatureOpenMS>;
    using FeaturePointerList_ = std::vector<const Feature*>;
    using FeatureIndexList_ = std::vector<std::pair<std::string, std::size_t>>;

    std::size_t findFeatureIndex_(const std::string& native_id) const;

    std::size_t findPrecursorFeatureIndex_(const std::string& native_id) const;

    const Feature& getFeatureByIndex_(std::size_t index) const;

    void ensureFeatureWrappers_();

    void ensurePrecursorFeatureWrappers_();

    const MRMFeature& mrmfeature_;
    std::shared_ptr<FeatureWrapperList_> features_;
    std::shared_ptr<FeatureWrapperList_> precursor_features_;
    const std::vector<Feature>* direct_features_ = nullptr;
    FeaturePointerList_ feature_ptrs_;
    FeaturePointerList_ precursor_feature_ptrs_;
    FeatureIndexList_ feature_index_;
    FeatureIndexList_ precursor_feature_index_;
    std::vector<std::string> feature_ids_storage_;
    std::vector<std::string> precursor_feature_ids_storage_;
    const std::vector<std::string>* feature_ids_ = nullptr;
    const std::vector<std::string>* precursor_feature_ids_ = nullptr;
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
      chromatogram_(chromat), sn_()
    {
      OpenMS::Param snt_parameters = sn_.getParameters();
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

      sn_.setParameters(snt_parameters);
      sn_.init(chromatogram_);
    }

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
        return sn_.getSignalToNoise((Size) distance(chromatogram_.begin(),prev));
      }
      else
      {
        // iter is closer to the apex
        return sn_.getSignalToNoise((Size) distance(chromatogram_.begin(),iter));
      }
    }

private:

    const ContainerT& chromatogram_;
    OpenMS::SignalToNoiseEstimatorMedian< ContainerT > sn_;

  };

}
