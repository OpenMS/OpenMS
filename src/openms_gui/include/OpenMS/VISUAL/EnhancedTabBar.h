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

//QT
#include <QTabBar>
class QMouseEvent;
class QMimeData;

namespace OpenMS
{
  class String;

  /**
      @brief Convenience tab bar implementation

      This tab bar differs in several ways form the QTabBar:
      - you can close a tab by double-clicking it or through its context menu.
      - it works based on tab identifiers (a fixed id stored in tab data) rather than on tab indices, which might
        change when inserting or removing a tab.
      - it accepts all drag-and-drop actions and emits signals to handle them.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI EnhancedTabBar :
    public QTabBar
  {
    Q_OBJECT
public:
    /// Constructor
    EnhancedTabBar(QWidget * parent = nullptr);

    /// Destructor
    ~EnhancedTabBar() override;

    /// sets the text of the current tab
    void setTabText(const QString& text);

    /// Adds a new tab with the name @p text and the identifier @p id
    int addTab(const String & text, int id);

    /// Selects the tab with identifier @p id
    void show(int id);

public slots:
    /// Remove the tab with identifier @p id
    void removeId(int id);

signals:
    /// Signal that indicates that the current tab changed, giving the @p id of the Tab
    void currentIdChanged(int id);

    /// Signal that indicates that the tab with identifier @p id is requested to be removed (double click or context menu)
    void closeRequested(int id);

    /// Signal that is emitted, when a drag-and-drop action ends on a tab
    void dropOnTab(const QMimeData * data, QWidget * source, int id);

    /// Signal that is emitted, when a drag-and-drop action ends on the unused space on the right side of the tabs.
    void dropOnWidget(const QMimeData * data, QWidget * source);

protected:
    ///@name Reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QMouseEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    void dragEnterEvent(QDragEnterEvent * e) override;
    void dropEvent(QDropEvent * e) override;
    //@}

    /// Returns the QTabBar index of the tab at position @p pos. If there is no tab at that position -1 is returned.
    int tabAt_(const QPoint & pos);

protected slots:
    /// Slot that translates the currentChanged(int) signal to currentIdChanged(int)
    void currentChanged_(int id);
  };

}
