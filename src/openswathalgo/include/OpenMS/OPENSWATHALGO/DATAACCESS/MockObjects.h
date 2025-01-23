// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ITransition.h>

#include <boost/shared_ptr.hpp>
#include <map>
#include <vector>
#include <string>

namespace OpenSwath
{

  /**
    @brief Mock object implementing IFeature
  */
  class OPENSWATHALGO_DLLAPI MockFeature :
    public OpenSwath::IFeature
  {
public:

    MockFeature();

    ~MockFeature() override;

    void getRT(std::vector<double>& rt) const override;

    void getIntensity(std::vector<double>& intens) const override;

    float getIntensity() const override;

    double getRT() const override;

    std::vector<double> m_rt_vec;
    std::vector<double> m_intensity_vec;
    float m_intensity;
    double m_rt;
  };

  /**
    @brief Mock object implementing IMRMFeature
  */
  class OPENSWATHALGO_DLLAPI MockMRMFeature :
    public OpenSwath::IMRMFeature
  {
public:

    MockMRMFeature();

    ~MockMRMFeature() override;

    boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) override;

    boost::shared_ptr<OpenSwath::IFeature> getPrecursorFeature(std::string nativeID) override;

    std::vector<std::string> getNativeIDs() const override;

    std::vector<std::string> getPrecursorIDs() const override;

    float getIntensity() const override;

    double getRT() const override;

    size_t size() const override;

    std::map<std::string, boost::shared_ptr<MockFeature> > m_features;
    std::map<std::string, boost::shared_ptr<MockFeature> > m_precursor_features;
    float m_intensity;
    double m_rt;
  };

  /**
    @brief Mock object implementing ITransitionGroup
  */
  class OPENSWATHALGO_DLLAPI MockTransitionGroup :
    public OpenSwath::ITransitionGroup
  {
public:

    MockTransitionGroup();

    ~MockTransitionGroup() override;

    std::size_t size() const override;

    std::vector<std::string> getNativeIDs() const override;

    void getLibraryIntensities(std::vector<double>& intensities) const override;

    std::size_t m_size;
    std::vector<std::string> m_native_ids;
    std::vector<double> m_library_intensities;
  };


  /**
    @brief Mock object implementing ISignalToNoise
  */
  class OPENSWATHALGO_DLLAPI MockSignalToNoise :
    public OpenSwath::ISignalToNoise
  {
public:
    MockSignalToNoise();

    double getValueAtRT(double /* RT */) override;

    double m_sn_value;
  };

} //end namespace OpenMS

