// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>

#include <QtWidgets/QGraphicsView>

namespace OpenMS
{
  class TOPPASScene;
  class Param;

  /**
    @brief Widget visualizing and allowing to edit TOPP pipelines.

    This class is a subclass of QGraphicsView and visualizes a TOPPASScene.
    Several TOPPASWidgets can be opened in TOPPAS at the same time,
    managed by a QWorkspace.

        @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASWidget :
    public QGraphicsView,
    public EnhancedTabBarWidgetInterface
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASWidget(const Param & preferences, QWidget * parent = nullptr, const String & tmp_path = "");

    /// Destructor
    ~TOPPASWidget() override;

    /// Returns the scene
    TOPPASScene * getScene();
    /// Zooms in or out, depending on @p zoom_in
    void zoom(bool zoom_in);

signals:

    /// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
    void sendStatusMessage(std::string message, OpenMS::UInt time);
    /// Emitted when the cursor position changes (for displaying e.g. in status bar)
    void sendCursorStatus(double x = 0.0, double y = 0.0);
    /// Emitted when a drop event occurs
    void toolDroppedOnWidget(double x = 0.0, double y = 0.0);
    /// Emitted when a drop event occurs
    void pipelineDroppedOnWidget(const String & filename, bool new_window);

protected:

    /// The scene visualized by this widget
    TOPPASScene * scene_;

    ///@name reimplemented QT events
    //@{
    void wheelEvent(QWheelEvent * event) override;
    void keyPressEvent(QKeyEvent * e) override;
    void keyReleaseEvent(QKeyEvent * e) override;
    void leaveEvent(QEvent * e) override;
    void enterEvent(QEnterEvent * e) override;
    void dragEnterEvent(QDragEnterEvent * event) override;
    void dragMoveEvent(QDragMoveEvent * event) override;
    void dropEvent(QDropEvent * event) override;
    void resizeEvent(QResizeEvent * event) override;
    void closeEvent(QCloseEvent * e) override;
    //@}
  };
}

