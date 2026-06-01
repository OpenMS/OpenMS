// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/TVControllerBase.h>
#include <vector>

namespace OpenMS
{
  class TOPPViewBase;
  class Plot1DWidget;
  struct MiniLayer;
  struct OSWIndexTrace;

  /**
  @brief Behavior of TOPPView in spectra view mode.
  */
  class TVDIATreeTabController
    : public TVControllerBase
  {
    Q_OBJECT

public:
    /// Construct the behaviour with its parent
  TVDIATreeTabController(TOPPViewBase* parent);

public slots:
    /// shows all transitions from the given subtree
    virtual void showChromatograms(const OSWIndexTrace& trace);

    /// shows all transitions from the given subtree in a new 1D canvas
    virtual void showChromatogramsAsNew1D(const OSWIndexTrace& trace);
    
private:    
    bool showChromatogramsInCanvas_(const OSWIndexTrace& trace, MiniLayer& ml, Plot1DWidget* w);
  };
}
