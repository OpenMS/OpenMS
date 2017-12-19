// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_MRMFEATUREACCESSOPENMS_H
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_MRMFEATUREACCESSOPENMS_H

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

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

    void getRT(std::vector<double>& rt) override;

    void getIntensity(std::vector<double>& intens) override;

    float getIntensity() override;

    double getRT() override;

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

    float getIntensity() override;

    double getRT() override;

    size_t size() override;

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

    std::size_t size() override
    {
      return trgroup_.size();
    }

    std::vector<std::string> getNativeIDs() override
    {
      std::vector<std::string> result;
      for (std::size_t i = 0; i < this->size(); i++)
      {
        result.push_back(trgroup_.getChromatograms()[i].getNativeID());
      }
      return result;
    }

    void getLibraryIntensities(std::vector<double>& intensities) override
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
      typename ContainerT::const_iterator iter = chromatogram_.MZEnd(RT);

      // ensure that iter is valid
      if (iter == chromatogram_.end()) 
      {
        iter--;
      }

      typename ContainerT::const_iterator prev = iter;
      if (prev != chromatogram_.begin() ) 
      {
        prev--;
      }

      if (std::fabs(prev->getMZ() - RT) < std::fabs(iter->getMZ() - RT) )
      {
        // prev is closer to the apex
        return sn_.getSignalToNoise(*prev);
      }
      else
      {
        // iter is closer to the apex
        return sn_.getSignalToNoise(*iter);
      }
    }

private:

    const ContainerT& chromatogram_;
    OpenMS::SignalToNoiseEstimatorMedian< ContainerT > sn_;

  };

}

#endif // OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_MRMFEATUREACCESSOPENMS_H

