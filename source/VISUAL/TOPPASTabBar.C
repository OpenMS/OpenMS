// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASTabBar.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QMouseEvent>
#include <QtGui/QMenu>
#include <QtGui/QMessageBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  TOPPASTabBar::TOPPASTabBar(QWidget * parent) :
    QTabBar(parent)
  {
    connect(this, SIGNAL(currentChanged(int)), this, SLOT(currentChanged_(int)));

    //set up drag-and-drop TODO
    //setAcceptDrops(true);
  }

  TOPPASTabBar::~TOPPASTabBar()
  {

  }

//  void TOPPASTabBar::dragEnterEvent(QDragEnterEvent* e)
//  {
//      e->acceptProposedAction();
//  }

//  void TOPPASTabBar::dropEvent(QDropEvent* e)
//  {
//      int tab = tabAt_(e->pos());
//      if (tab!=-1)
//      {
//          emit dropOnTab(e->mimeData(), e->source(), tabData(tab).toInt());
//      }
//      else
//      {
//        emit dropOnWidget(e->mimeData(), e->source());
//      }
//
//      e->acceptProposedAction();
//  }

  void TOPPASTabBar::contextMenuEvent(QContextMenuEvent * e)
  {
    int tab = tabAt_(e->pos());
    if (tab != -1)
    {
      QMenu menu(this);
      menu.addAction("Close");
      if (menu.exec(e->globalPos()))
      {
        emit aboutToCloseId(tabData(tab).toInt());
      }
    }
  }

  void TOPPASTabBar::mouseDoubleClickEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }
    int tab = tabAt_(e->pos());
    if (tab != -1)
    {
      emit aboutToCloseId(tabData(tab).toInt());
    }
  }

  int TOPPASTabBar::addTab(const String & text, int id)
  {
    int tab_index = QTabBar::addTab(text.c_str());
    setTabData(tab_index, id);

    return tab_index;
  }

  void TOPPASTabBar::removeId(int id)
  {
    for (int i = 0; i < this->count(); ++i)
    {
      if (tabData(i).toInt() == id)
      {
        removeTab(i);
        break;
      }
    }
  }

  void TOPPASTabBar::setCurrentId(int id)
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

  void TOPPASTabBar::currentChanged_(int id)
  {
    emit currentIdChanged(tabData(id).toInt());
  }

  int TOPPASTabBar::tabAt_(const QPoint & pos)
  {
    int tab = -1;

    for (int i = 0; i < this->count(); ++i)
    {
      if (tabRect(i).contains(pos))
      {
        tab = i;
        break;
      }
    }

    return tab;
  }

} //namespace OpenMS
