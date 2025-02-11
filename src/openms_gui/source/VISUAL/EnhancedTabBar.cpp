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
    int tab = tabAt_(e->pos());
    if (tab != -1)
    {
      emit dropOnTab(e->mimeData(), dynamic_cast<QWidget*>(e->source()), tabData(tab).toInt());
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
        emit closeRequested(tabData(tab).toInt());
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
      emit closeRequested(tabData(tab).toInt());
    }
  }

  int EnhancedTabBar::addTab(const String& text, int id)
  {
    // make sure this ID does not exist yet
    for (int i = 0; i < this->count(); ++i)
    {
      if (tabData(i).toInt() == id)
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
      if (tabData(i).toInt() == id)
      {
        removeTab(i);
        return;
      }
    }
   throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Tab with ID ") + id + " is already gone!");
  }

  void EnhancedTabBar::show(int id)
  {
    for (int i = 0; i < this->count(); ++i)
    {
      if (tabData(i).toInt() == id)
      {
        setCurrentIndex(i);
        break;
      }
    }
  }

  void EnhancedTabBar::currentChanged_(int index)
  {
    emit currentIdChanged(tabData(index).toInt());
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
