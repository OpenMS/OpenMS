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

#include <QMdiArea>

class QMimeData;
class QDragEnterEvent;
class QDragMoveEvent;
class QDropEvent;

namespace OpenMS
{
  class EnhancedTabBarWidgetInterface;

  class OPENMS_GUI_DLLAPI EnhancedWorkspace :
    public QMdiArea
  {
    Q_OBJECT

public:
    /// Constructor
    EnhancedWorkspace(QWidget * parent);

    /// Destructor
    ~EnhancedWorkspace() override;

    /// Adds a subwindow to the QMdiArea
    /// Qt will add a System menu (which shows when you right-click on the subwindows' local top bar).
    /// This menu will contain a shortcut for Ctrl-W, which makes our custom Ctrl-W in TOPPView's menu ambiguous.
    /// Upon pressing it, you will get a `QAction::event: Ambiguous shortcut overload: Ctrl+W` on the console and no triggered signal.
    /// To prevent this we simply set the SystemMenu to null (no menu will show when you right-click).
    /// If you do not want that, call Qt's overload of addSubWindow
    QMdiSubWindow* addSubWindow(QWidget* widget);

    /// arrange all windows horizontally
    void tileHorizontal();

    /// arrange all windows vertically
    void tileVertical();

    /// get the subwindow with the given id (for all subwindows which inherit from EnhancedTabBarWidgetInterface)
    /// Returns nullptr if window is not present
    EnhancedTabBarWidgetInterface* getWidget(int id) const;

signals:

    /// Signal that is emitted, when a drag-and-drop action ends on this widget
    void dropReceived(const QMimeData * data, QWidget * source, int id);

protected:

    ///@name Reimplemented Qt events
    //@{
    void dragEnterEvent(QDragEnterEvent * event) override;
    void dragMoveEvent(QDragMoveEvent * event) override;
    void dropEvent(QDropEvent * event) override;
    //@}
  };
}

