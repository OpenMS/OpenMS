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

#include <boost/shared_ptr.hpp>

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

    explicit FeatureOpenMS(Feature& feature);

    ~FeatureOpenMS() override;

    void getRT(std::vector<double>& rt) const override;

    void getIntensity(std::vector<double>& intens) const override;

    float getIntensity() const override;

    double getRT() const override;

private:
    Feature* feature_;
  };

  /**
    @brief An implementation of the OpenSWATH MRM Feature Access interface using OpenMS

  */
  class OPENMS_DLLAPI MRMFeatureOpenMS :
    public OpenSwath::IMRMFeature
  {
public:

    explicit MRMFeatureOpenMS(MRMFeature& mrmfeature);

    ~MRMFeatureOpenMS() override;

    boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) override;

    boost::shared_ptr<OpenSwath::IFeature> getPrecursorFeature(std::string nativeID) override;

    std::vector<std::string> getNativeIDs() const override;

    std::vector<std::string> getPrecursorIDs() const override;

    float getIntensity() const override;

    double getRT() const override;

    size_t size() const override;

private:
    const MRMFeature& mrmfeature_;
    std::map<std::string, boost::shared_ptr<FeatureOpenMS> > features_;
    std::map<std::string, boost::shared_ptr<FeatureOpenMS> > precursor_features_;
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

    SignalToNoiseOpenMS(ContainerT& chromat,
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


