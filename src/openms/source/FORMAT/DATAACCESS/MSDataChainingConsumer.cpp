// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataChainingConsumer.h>

#include <utility>

namespace OpenMS
{

  MSDataChainingConsumer::MSDataChainingConsumer() = default;

  MSDataChainingConsumer::MSDataChainingConsumer(std::vector<Interfaces::IMSDataConsumer *> consumers) :
    consumers_(std::move(consumers))
  {}

  MSDataChainingConsumer::~MSDataChainingConsumer() = default;

  void MSDataChainingConsumer::appendConsumer(Interfaces::IMSDataConsumer * consumer)
  {
    consumers_.push_back(consumer);
  }

  void MSDataChainingConsumer::setExperimentalSettings(const ExperimentalSettings & settings)
  {
    for (Size i = 0; i < consumers_.size(); i++)
    {
      consumers_[i]->setExperimentalSettings(settings);
    }
  }

  void MSDataChainingConsumer::setExpectedSize(Size s_size, Size c_size) 
  {
    for (Size i = 0; i < consumers_.size(); i++)
    {
      consumers_[i]->setExpectedSize(s_size, c_size);
    }
  }

  void MSDataChainingConsumer::consumeSpectrum(SpectrumType & s)
  {
    for (Size i = 0; i < consumers_.size(); i++)
    {
      consumers_[i]->consumeSpectrum(s);
    }
  }

  void MSDataChainingConsumer::consumeChromatogram(ChromatogramType & c)
  {
    for (Size i = 0; i < consumers_.size(); i++)
    {
      consumers_[i]->consumeChromatogram(c);
    }
  }

} //end namespace OpenMS

