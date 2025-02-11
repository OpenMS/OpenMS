// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h>


namespace OpenSwath
{

  MockFeature::MockFeature()
  {
  }

  MockFeature::~MockFeature()
  {
  }

  void MockFeature::getRT(std::vector<double>& rt) const
  {
    rt = m_rt_vec;
  }

  void MockFeature::getIntensity(std::vector<double>& intens) const
  {
    intens = m_intensity_vec;
  }

  float MockFeature::getIntensity() const
  {
    return m_intensity;
  }

  double MockFeature::getRT() const
  {
    return m_rt;
  }

  MockMRMFeature::MockMRMFeature()
  {
  }

  MockMRMFeature::~MockMRMFeature()
  {
  }

  boost::shared_ptr<OpenSwath::IFeature> MockMRMFeature::getFeature(std::string nativeID)
  {
    return boost::static_pointer_cast<OpenSwath::IFeature>(m_features[nativeID]);
  }

  boost::shared_ptr<OpenSwath::IFeature> MockMRMFeature::getPrecursorFeature(std::string nativeID)
  {
    return boost::static_pointer_cast<OpenSwath::IFeature>(m_precursor_features[nativeID]);
  }

  std::vector<std::string> MockMRMFeature::getNativeIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<MockFeature> >::const_iterator
         it = m_features.begin(); it != m_features.end(); ++it)
    {
      v.push_back(it->first);
    }
    return v;
  }

  std::vector<std::string> MockMRMFeature::getPrecursorIDs() const
  {
    std::vector<std::string> v;
    for (std::map<std::string, boost::shared_ptr<MockFeature> >::const_iterator
         it = m_precursor_features.begin(); it != m_precursor_features.end(); ++it)
    {
      v.push_back(it->first);
    }
    return v;
  }

  float MockMRMFeature::getIntensity() const
  {
    return m_intensity;
  }

  double MockMRMFeature::getRT() const
  {
    return m_rt;
  }

  size_t MockMRMFeature::size() const
  {
    return m_features.size();
  }

  MockTransitionGroup::MockTransitionGroup()
  {
  }

  MockTransitionGroup::~MockTransitionGroup()
  {
  }

  std::size_t MockTransitionGroup::size() const
  {
    return m_size;
  }

  std::vector<std::string> MockTransitionGroup::getNativeIDs() const
  {
    return m_native_ids;
  }

  void MockTransitionGroup::getLibraryIntensities(std::vector<double>& intensities) const
  {
    intensities = m_library_intensities;
  }

  MockSignalToNoise::MockSignalToNoise()
  {
  }

  double MockSignalToNoise::getValueAtRT(double /* RT */)
  {
    return m_sn_value;
  }

}
