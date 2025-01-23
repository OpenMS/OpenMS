// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerDataBase.h>

namespace OpenMS
{
  class TOPPViewBase;

  /**
  @brief Base behavior for different visualizaton modules in TOPPView.
  */
  class TVControllerBase
    : public QObject
  {
    Q_OBJECT

public:
    ///@name Type definitions
    //@{
    /// Feature map type
    typedef LayerDataBase::FeatureMapType FeatureMapType;
    /// Feature map managed type
    typedef LayerDataBase::FeatureMapSharedPtrType FeatureMapSharedPtrType;

    /// Consensus feature map type
    typedef LayerDataBase::ConsensusMapType ConsensusMapType;
    /// Consensus  map managed type
    typedef LayerDataBase::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

    /// Peak map type
    typedef LayerDataBase::ExperimentType ExperimentType;
    /// Main managed data type (experiment)
    typedef LayerDataBase::ExperimentSharedPtrType ExperimentSharedPtrType;
    /// Peak spectrum type
    typedef ExperimentType::SpectrumType SpectrumType;
    //@}
    TVControllerBase() = delete;

    ~TVControllerBase() override = default;
public slots:
    /// Slot for behavior activation. The default behaviour does nothing. Override in child class if desired.
    virtual void activateBehavior();

    /// Slot for behavior deactivation. The default behaviour does nothing. Override in child class if desired.
    virtual void deactivateBehavior();
protected:
    /// Construct the behaviour
    TVControllerBase(TOPPViewBase* parent);

    TOPPViewBase* tv_;
  };
}
