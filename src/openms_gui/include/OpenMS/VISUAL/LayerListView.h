// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <QtWidgets/QListWidget>

namespace OpenMS
{
  class PlotWidget;
  /**
    @brief Pimped QListView for Layers of a Canvas


  */
  class OPENMS_GUI_DLLAPI LayerListView
    : public QListWidget
  {
    Q_OBJECT

  public:
    /// Default constructor
    LayerListView(QWidget* parent);

    /// rebuild list of layers and remember current widget (for context menu etc)
    void update(PlotWidget* active_widget);

  signals:
    /// emitted whenever a change to a layer happened, e.g. its name was changed, it was removed, or a new layer was selected
    void layerDataChanged();

  private:
    /// active row was changed by user to new row @p i
    void currentRowChangedAction_(int i);

    void itemChangedAction_(QListWidgetItem* item);

    void contextMenuEvent(QContextMenuEvent* event) override;

    /// show preferences dialog
    void itemDoubleClickedAction_(QListWidgetItem*);

    PlotWidget* spectrum_widget_ = nullptr; ///< holds the actual data. Might be nullptr.
  };

} //namespace

