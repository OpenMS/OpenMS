// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASTreeView.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QDrag>
#include <QApplication>
#include <QtCore/QMimeData>

#include <iostream>

using namespace std;

namespace OpenMS
{

  TOPPASTreeView::TOPPASTreeView(QWidget * parent) :
    QTreeWidget(parent)
  {
    // we drag by ourselves:
    setDragEnabled(false);
  }

  TOPPASTreeView::~TOPPASTreeView() = default;

  void TOPPASTreeView::mousePressEvent(QMouseEvent * event)
  {
    QTreeWidget::mousePressEvent(event);

    if (event->button() == Qt::LeftButton)
    {
      drag_start_pos_ = event->pos();
    }
  }

  void TOPPASTreeView::mouseMoveEvent(QMouseEvent * event)
  {
    QTreeWidget::mouseMoveEvent(event);

    if (!(event->buttons() & Qt::LeftButton))
    {
      return;
    }
    if ((event->pos() - drag_start_pos_).manhattanLength() < QApplication::startDragDistance())
    {
      return;
    }
    if (currentItem() && currentItem()->childCount() > 0)
    {
      // drag item is a category or a tool with types - one of the types must be selected
      return;
    }

    QDrag * drag = new QDrag(this);
    QMimeData * mime_data = new QMimeData;

    mime_data->setText(currentItem()->text(0));
    drag->setMimeData(mime_data);

    // start drag
    drag->exec(Qt::CopyAction);
  }

  void TOPPASTreeView::keyPressEvent(QKeyEvent * e)
  {
    QTreeWidget::keyPressEvent(e);
    if (currentItem() && e->key() == Qt::Key_Return)
    {
      e->accept();
      emit itemDoubleClicked(currentItem(), 0);
    }
    else
    {
      e->ignore();
    }
  }

  void TOPPASTreeView::enterEvent(QEvent * /*e*/)
  {
    setFocus();
  }

  void TOPPASTreeView::leaveEvent(QEvent * /*e*/)
  {

  }

} //namespace OpenMS
