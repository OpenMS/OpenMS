// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

namespace OpenSwath
{
  // Datastructures for Scoring
  class OPENSWATHALGO_DLLAPI IFeature
  {
public:
    virtual ~IFeature(){}
    virtual void getRT(std::vector<double>& rt) const = 0;
    virtual void getIntensity(std::vector<double>& intens) const = 0;
    virtual float getIntensity() const = 0;
    virtual double getRT() const = 0;
  };

  class OPENSWATHALGO_DLLAPI IMRMFeature
  {
public:
    virtual ~IMRMFeature(){}
    virtual boost::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID) = 0;
    virtual boost::shared_ptr<OpenSwath::IFeature> getPrecursorFeature(std::string nativeID) = 0;
    virtual std::vector<std::string> getNativeIDs() const = 0;
    virtual std::vector<std::string> getPrecursorIDs() const = 0;
    virtual float getIntensity() const = 0;
    virtual double getRT() const = 0;
    virtual size_t size() const = 0;
  };

  struct OPENSWATHALGO_DLLAPI ITransitionGroup
  {
    virtual ~ITransitionGroup() {}
    virtual std::size_t size() const = 0;
    virtual std::vector<std::string> getNativeIDs() const = 0;
    virtual void getLibraryIntensities(std::vector<double>& intensities) const = 0;
  };

  struct OPENSWATHALGO_DLLAPI ISignalToNoise
  {
    virtual ~ISignalToNoise() {}
    virtual double getValueAtRT(double RT) = 0; // cannot be const due to OpenMS implementation
  };
  typedef boost::shared_ptr<ISignalToNoise> ISignalToNoisePtr;


} //end Namespace OpenSwath

