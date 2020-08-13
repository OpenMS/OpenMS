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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>

namespace OpenMS
{

  MRMFeatureOpenMS::MRMFeatureOpenMS(MRMFeature& mrmfeature) :
    mrmfeature_(mrmfeature)
  {
    std::vector<String> ids;
    mrmfeature.getFeatureIDs(ids);
    for (std::vector<String>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      boost::shared_ptr<FeatureOpenMS> ptr = boost::shared_ptr<FeatureOpenMS>(new FeatureOpenMS(mrmfeature.getFeature(*it)));
      features_[*it] = ptr;
    }

    std::vector<String> p_ids;
    mrmfeature.getPrecursorFeatureIDs(p_ids);
    for (std::vector<String>::iterator it = p_ids.begin(); it != p_ids.end(); ++it)
    {
      boost::shared_ptr<FeatureOpenMS> ptr = boost::shared_ptr<FeatureOpenMS>(new FeatureOpenMS(mrmfeature.getPrecursorFeature(*it)));
      precursor_features_[*it] = ptr;
    }
  }

  FeatureOpenMS::FeatureOpenMS(Feature& feature)
  {
    feature_ = &feature; // store raw ptr to the feature
  }

  FeatureOpenMS::~FeatureOpenMS()
  {
  }

  void FeatureOpenMS::getRT(std::vector<double>& rt) const
  {
    OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
    ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
    for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); ++it)
    {
      rt.push_back(it->getX());
    }
  }

  void FeatureOpenMS::getIntensity(std::vector<double>& intens) const
  {
    OPENMS_PRECONDITION(feature_->getConvexHulls().size() == 1, "There needs to exactly one convex hull per feature.");
    ConvexHull2D::PointArrayType data_points = feature_->getConvexHulls()[0].getHullPoints();
    for (ConvexHull2D::PointArrayType::iterator it = data_points.begin(); it != data_points.end(); ++it)
    {
      intens.push_back(it->getY());
    }
  }

  float FeatureOpenMS::getIntensity() const
  {
    return feature_->getIntensity();
  }

  double FeatureOpenMS::getRT() const
  {
    return feature_->getRT();
  }

  MRMFeatureOpenMS::~MRMFeatureOpenMS()
  {
  }

  boost::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getFeature(std::string nativeID)
  {
    OPENMS_PRECONDITION(features_.find(nativeID) != features_.end(), "Feature needs to exist");
    return boost::static_pointer_cast<OpenSwath::IFeature>(features_[nativeID]);
  }

  boost::shared_ptr<OpenSwath::IFeature> MRMFeatureOpenMS::getPrecursorFeature(std::string nativeID)
  {
    OPENMS_PRECONDITION(precursor_features_.find(nativeID) != precursor_features_.end(), "Precursor feature needs to exist");
    return boost::static_pointer_cast<OpenSwath::IFeature>(precursor_features_[nativeID]);
  }

  std::vector<std::string> MRMFeatureOpenMS::getNativeIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<FeatureOpenMS> >::const_iterator it = features_.begin(); it != features_.end(); ++it)
    {
      v.push_back(it->first);
    }
    return v;
  }

  std::vector<std::string> MRMFeatureOpenMS::getPrecursorIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<FeatureOpenMS> >::const_iterator it = precursor_features_.begin(); it != precursor_features_.end(); ++it) 
    {
      v.push_back(it->first);
    }
    return v;
  }

  float MRMFeatureOpenMS::getIntensity() const
  {
    return mrmfeature_.getIntensity();
  }

  double MRMFeatureOpenMS::getRT() const
  {
    return mrmfeature_.getRT();
  }

  size_t MRMFeatureOpenMS::size() const
  {
    return features_.size();
  }

  // default instances
  MSSpectrum spec;
  MSChromatogram chrom;
  SignalToNoiseOpenMS< MSSpectrum> spec_signal_to_noise_openms(spec, 1.0, 3, true);
  SignalToNoiseOpenMS< MSChromatogram > chrom_signal_to_noise_openms(chrom, 1.0, 3, true);

  MRMTransitionGroup<MSSpectrum, ReactionMonitoringTransition> trgroup;
  TransitionGroupOpenMS<MSSpectrum, ReactionMonitoringTransition> default_transition_group_openms(trgroup);

} //end namespace OpenMS

