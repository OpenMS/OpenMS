// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_MRMFEATUREACCESSIMPL_H_
#define OPENMS_ANALYSIS_OPENSWATH_DATAACCESS_MRMFEATUREACCESSIMPL_H_

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ITransition.h>

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include "boost/shared_ptr.hpp"

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

    FeatureOpenMS(Feature & feature)
    {
      feature_ = &feature;   // store raw ptr to the feature
    }

    ~FeatureOpenMS()
    {
    }

    void getRT(std::vector<double> & rt)
    {
      OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
      ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
      for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); it++)
      {
        rt.push_back(it->getX());
      }
    }

    void getIntensity(std::vector<double> & intens)
    {
      OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
      ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
      for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); it++)
      {
        intens.push_back(it->getY());
      }
    }

    float getIntensity()
    {
      return feature_->getIntensity();
    }

    double getRT()
    {
      return feature_->getRT();
    }

private:
    Feature * feature_;
  };

  /**
    @brief An implementation of the OpenSWATH MRM Feature Access interface using OpenMS

  */
  class OPENMS_DLLAPI MRMFeatureOpenMS :
    public OpenSwath::IMRMFeature
  {
public:

    MRMFeatureOpenMS(MRMFeature & mrmfeature);

    ~MRMFeatureOpenMS()
    {
    }

    boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID)
    {
      return boost::static_pointer_cast<OpenSwath::IFeature>(features_[nativeID]);
    }

    float getIntensity()
    {
      return mrmfeature_.getIntensity();
    }

    double getRT()
    {
      return mrmfeature_.getRT();
    }

private:
    const MRMFeature & mrmfeature_;
    std::map<std::string, boost::shared_ptr<FeatureOpenMS> > features_;
  };

  /**
    @brief An implementation of the OpenSWATH Transition Group Access interface using OpenMS

  */
  template <template <typename> class SpectrumT, typename PeakT, typename TransitionT>
  class OPENMS_DLLAPI TransitionGroupOpenMS :
    public OpenSwath::ITransitionGroup
  {
public:

    TransitionGroupOpenMS(MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & trgroup) :
      trgroup_(trgroup)
    {
    }

    ~TransitionGroupOpenMS()
    {
    }

    std::size_t size()
    {
      return trgroup_.size();
    }

    std::vector<std::string> getNativeIDs()
    {
      std::vector<std::string> result;
      for (std::size_t i = 0; i < this->size(); i++)
      {
        result.push_back(trgroup_.getChromatograms()[i].getNativeID());
      }
      return result;
    }

    void getLibraryIntensities(std::vector<double> & intensities)
    {
      trgroup_.getLibraryIntensity(intensities);
    }

private:
    const MRMTransitionGroup<SpectrumT, PeakT, TransitionT> & trgroup_;
  };

  /**
    @brief An implementation of the OpenSWATH SignalToNoise Access interface using OpenMS

  */
  template <typename PeakT>
  class OPENMS_DLLAPI SignalToNoiseOpenMS :
    public OpenSwath::ISignalToNoise
  {
public:

    SignalToNoiseOpenMS(OpenMS::MSSpectrum<PeakT> & chromat,
                        double sn_win_len_, unsigned int sn_bin_count_) :
      chromatogram_(chromat), sn_()
    {
      OpenMS::Param snt_parameters = sn_.getParameters();
      snt_parameters.setValue("win_len", sn_win_len_);
      snt_parameters.setValue("bin_count", sn_bin_count_);
      sn_.setParameters(snt_parameters);
      sn_.init(chromatogram_);
    }

    double getValueAtRT(double RT)
    {
      typename OpenMS::MSSpectrum<PeakT>::const_iterator it = chromatogram_.MZBegin(RT);
      return sn_.getSignalToNoise(*it);
    }

private:
    const OpenMS::MSSpectrum<PeakT> & chromatogram_;
    OpenMS::SignalToNoiseEstimatorMedian<OpenMS::MSSpectrum<PeakT> > sn_;
  };

}

#endif
