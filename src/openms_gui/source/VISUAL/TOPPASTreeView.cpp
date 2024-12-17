// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

  void TOPPASTreeView::filter(const QString& must_match)
  {
    // hide all
    QTreeWidgetItemIterator it(this);
    while (*it)
    {
      (*it)->setHidden(true);
      (*it)->setExpanded(false);
      ++it;
    }

    // recursive lambda: show items and its subchildren (e.g. when a category matches)
    function<void(QTreeWidgetItem*)> show_sub_tree = [&](QTreeWidgetItem* item) {
      item->setHidden(false);
      for (int i = 0; i < item->childCount(); i++)
      {
        QTreeWidgetItem* child = item->child(i);
        child->setHidden(false);
        child->setExpanded(true); // technically not required, since our tree is only 2 layers deep, but maybe in the future...
        show_sub_tree(child);
      }
    };

    // show stuff that matches
    auto items = this->findItems(must_match, Qt::MatchContains | Qt::MatchRecursive);
    for (auto& it : items)
    {
      // show parent (if any) -- otherwise the children will not be displayed
      if (it->parent())
      {
        it->parent()->setHidden(false);
        it->parent()->setExpanded(true);
      }
      show_sub_tree(it); // also show all children
      it->setExpanded(true);
    }
  }

  void TOPPASTreeView::expandAll()
  {
    QTreeWidgetItemIterator it(this);
    while (*it)
    {
      (*it)->setExpanded(true);
      ++it;
    }
  }

  void TOPPASTreeView::collapseAll()
  {
    QTreeWidgetItemIterator it(this);
    while (*it)
    {
      (*it)->setExpanded(false);
      ++it;
    }
  }

  void TOPPASTreeView::mousePressEvent(QMouseEvent* event)
  {
    QTreeWidget::mousePressEvent(event);

    if (event->button() == Qt::LeftButton)
    {
      drag_start_pos_ = event->pos();
    }
  }

  void TOPPASTreeView::mouseMoveEvent(QMouseEvent* event)
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

  void TOPPASTreeView::keyPressEvent(QKeyEvent* e)
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

  void TOPPASTreeView::enterEvent(QEnterEvent* /*e*/)
  {
    setFocus();
  }

  void TOPPASTreeView::leaveEvent(QEvent* /*e*/)
  {
  }

} //namespace OpenMS
