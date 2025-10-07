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
#include <OpenMS/VISUAL/TVControllerBase.h>
#include <vector>

namespace OpenMS
{
  class TOPPViewBase;

  /**
  @brief Behavior of TOPPView in spectra view mode.
  */
  class TVSpectraViewController
    : public TVControllerBase
  {
    Q_OBJECT

public:
    /// Construct the behaviour with its parent
    TVSpectraViewController(TOPPViewBase* parent);

public slots:
    /// Behavior for showSpectrumAsNew1D
    virtual void showSpectrumAsNew1D(int index);

    /// Behavior for showChromatogramsAsNew1D
    virtual void showChromatogramsAsNew1D(const std::vector<int>& indices);

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(int index);

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(const std::vector<int>& indices);

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index);
  };
}
