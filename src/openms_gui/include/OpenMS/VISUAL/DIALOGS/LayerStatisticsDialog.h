// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <QtWidgets/QDialog>

#include <memory> // for unique_ptr

namespace Ui
{
  class LayerStatisticsDialogTemplate;
}

namespace OpenMS
{
  class LayerStatistics;
  class PlotWidget;
  class PlotCanvas;
  /**
      @brief Dialog showing statistics about the data of the current layer

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /// Constructor not implemented
    LayerStatisticsDialog() = delete;
    /// Custom constructor
    LayerStatisticsDialog(PlotWidget* parent, std::unique_ptr<LayerStatistics>&& stats);
    /// D'tor
    ~LayerStatisticsDialog() override;

protected:
    /// The statistics of the layer
    std::unique_ptr<LayerStatistics> stats_;

private:
    Ui::LayerStatisticsDialogTemplate* ui_;
  };
}
