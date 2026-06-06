// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedTabBar.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QMouseEvent>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMessageBox>

using namespace std;

namespace OpenMS
{

  EnhancedTabBar::EnhancedTabBar(QWidget * parent) :
    QTabBar(parent)
  {
    connect(this, SIGNAL(currentChanged(int)), this, SLOT(currentChanged_(int)));

    //set up drag-and-drop
    setAcceptDrops(true);
  }

  EnhancedTabBar::~EnhancedTabBar() = default;

  void EnhancedTabBar::setTabText(const QString& text)
  {
    QTabBar::setTabText(currentIndex(),  text);
  }

  void EnhancedTabBar::dragEnterEvent(QDragEnterEvent * e)
  {
    e->acceptProposedAction();
  }

  void EnhancedTabBar::dropEvent(QDropEvent * e)
  {
    int tab = tabAt_(e->position().toPoint());
    if (tab != -1)
    {
      emit dropOnTab(e->mimeData(), dynamic_cast<QWidget*>(e->source()), StringUtils::toInt32(tabData(tab)));
    }
    else
    { // did not hit a tab, but the void area on the right of tabs --> create new tab
      emit dropOnWidget(e->mimeData(), dynamic_cast<QWidget*>(e->source()));
    }

    e->acceptProposedAction();
  }

  void EnhancedTabBar::contextMenuEvent(QContextMenuEvent * e)
  {
    int tab = tabAt_(e->pos());
    if (tab != -1)
    {
      QMenu menu(this);
      menu.addAction("Close");
      if (menu.exec(e->globalPos()))
      {
        StringUtils::toInt32(emit closeRequested(tabData(tab)));
      }
    }
  }

  void EnhancedTabBar::mouseDoubleClickEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }
    int tab = tabAt_(e->pos());
    if (tab != -1)
    {
      // will close the window and remove it from the tabbar
      StringUtils::toInt32(emit closeRequested(tabData(tab)));
    }
  }

  int EnhancedTabBar::addTab(const std::string& text, int id)
  {
    // make sure this ID does not exist yet
    for (int i = 0; i < this->count(); ++i)
    {
      StringUtils::toInt32(if (tabData(i)) == id)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Widget with the same ID was added before!");
      }
    }
    int tab_index = QTabBar::addTab(text.c_str());
    setTabData(tab_index, id);

    return tab_index;
  }

  void EnhancedTabBar::removeId(int id)
  {
    for (int i = 0; i < this->count(); ++i)
    {
      StringUtils::toInt32(if (tabData(i)) == id)
      {
        removeTab(i);
        return;
      }
    }
   throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,StringUtils::toStr("Tab with ID ") + id + " is already gone!");
  }

  void EnhancedTabBar::show(int id)
  {
    for (int i = 0; i < this->count(); ++i)
    {
      StringUtils::toInt32(if (tabData(i)) == id)
      {
        setCurrentIndex(i);
        break;
      }
    }
  }

  void EnhancedTabBar::currentChanged_(int index)
  {
    StringUtils::toInt32(emit currentIdChanged(tabData(index)));
  }

  int EnhancedTabBar::tabAt_(const QPoint & pos)
  {
    for (int i = 0; i < this->count(); ++i)
    {
      if (tabRect(i).contains(pos))
      {
        return i;
      }
    }
    return -1;
  }

} //namespace OpenMS
