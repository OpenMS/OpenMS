// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedWorkspace.h>

#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>

#include <QtCore/QMimeData>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDropEvent>

#include <QMdiSubWindow>

namespace OpenMS
{

  EnhancedWorkspace::EnhancedWorkspace(QWidget* parent) :
    QMdiArea(parent)
  {
    setAcceptDrops(true);
  }

  EnhancedWorkspace::~EnhancedWorkspace() = default;

  QMdiSubWindow* EnhancedWorkspace::addSubWindow(QWidget* widget)
  {
    auto subwindow = QMdiArea::addSubWindow(widget);
    if (subwindow)
    {
      subwindow->setSystemMenu(nullptr);
    }
    return subwindow;
  }

  void EnhancedWorkspace::tileHorizontal()
  {
    // primitive horizontal tiling
    QList<QMdiSubWindow*> windows = this->subWindowList();

    if (!windows.count())
    {
      return;
    }

    int heightForEach = this->height() / windows.count();
    int y = 0;
    for (int i = 0; i < int(windows.count()); ++i)
    {
      QMdiSubWindow* window = windows.at(i);
      if (window->isMaximized() || window->isMinimized() || window->isFullScreen())
      {
        // prevent flicker
        window->hide();
        window->showNormal();
      }
      int preferredHeight = window->widget()->minimumHeight() + window->baseSize().height();
      int actHeight = std::max(heightForEach, preferredHeight);

      window->setGeometry(0, y, this->width(), actHeight);
      y += actHeight;
      window->setVisible(true);
      window->show();
    }
  }

  void EnhancedWorkspace::tileVertical()
  {
    // primitive vertical tiling
    QList<QMdiSubWindow*> windows = this->subWindowList();
    if (!windows.count())
    {
      return;
    }

    int widthForEach = this->width() / windows.count();
    int y = 0;
    for (int i = 0; i < int(windows.count()); ++i)
    {
      QMdiSubWindow* window = windows.at(i);
      if (window->isMaximized() || window->isMinimized() || window->isFullScreen())
      {
        // prevent flicker
        window->hide();
        window->showNormal();
      }
      int preferredWidth = window->widget()->minimumWidth() + window->baseSize().width();
      int actWidth = std::max(widthForEach, preferredWidth);

      window->setGeometry(y, 0, actWidth, this->height());
      y += actWidth;
      window->setVisible(true);
      window->show();
    }
  }

  /// get the subwindow with the given id (for all subwindows which inherit from EnhancedTabBarWidgetInterface)
  /// Returns nullptr if window is not present

  EnhancedTabBarWidgetInterface* EnhancedWorkspace::getWidget(int id) const
  {
    for (const auto& sub_window : this->subWindowList())
    {
      EnhancedTabBarWidgetInterface* w = dynamic_cast<EnhancedTabBarWidgetInterface*>(sub_window->widget());
      //cout << "  Tab " << i << ": " << w->window_id << endl;
      if (w != nullptr && w->getWindowId() == id)
      {
        return w;
      }
    }
    return nullptr;
  }

  void EnhancedWorkspace::dragEnterEvent(QDragEnterEvent * event)
  {
    if (event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
    }
  }

  void EnhancedWorkspace::dragMoveEvent(QDragMoveEvent * event)
  {
    if (event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
    }
  }

  void EnhancedWorkspace::dropEvent(QDropEvent * event)
  {
    emit dropReceived(event->mimeData(), dynamic_cast<QWidget*>(event->source()), -1);
    event->acceptProposedAction();
  }

} //namespace
