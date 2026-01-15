// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/CONCEPT/Types.h>


// Qt
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDropEvent>
#include <QtCore/QMimeData>
#include <QUrl>

using namespace std;

namespace OpenMS
{
  TOPPASWidget::TOPPASWidget(const Param & /*preferences*/, QWidget * parent, const String & tmp_path) :
    QGraphicsView(parent),
    EnhancedTabBarWidgetInterface(),
    scene_(new TOPPASScene(this, tmp_path.toQString()))
  {
    setAttribute(Qt::WA_DeleteOnClose);
    setAttribute(Qt::WA_AlwaysShowToolTips);
    setRenderHint(QPainter::Antialiasing);
    setScene(scene_);
    setAcceptDrops(true);
    setDragMode(QGraphicsView::ScrollHandDrag);
    setFocusPolicy(Qt::StrongFocus);
  }

  TOPPASWidget::~TOPPASWidget() = default;

  TOPPASScene * TOPPASWidget::getScene()
  {
    return scene_;
  }

  void TOPPASWidget::zoom(bool zoom_in)
  {
    qreal factor = 1.1;
    if (zoom_in)
    {
      factor = 1.0 / factor;
    }
    scale(factor, factor);

    QRectF items_rect = scene_->itemsBoundingRect();
    QRectF new_scene_rect = items_rect.united(mapToScene(rect()).boundingRect());
    qreal top_left_x = new_scene_rect.topLeft().x();
    qreal top_left_y = new_scene_rect.topLeft().y();
    qreal bottom_right_x = new_scene_rect.bottomRight().x();
    qreal bottom_right_y = new_scene_rect.bottomRight().y();
    qreal width = new_scene_rect.width();
    qreal height = new_scene_rect.height();
    new_scene_rect.setTopLeft(QPointF(top_left_x - width / 2.0, top_left_y - height / 2.0));
    new_scene_rect.setBottomRight(QPointF(bottom_right_x + width / 2.0, bottom_right_y + height / 2.0));
    scene_->setSceneRect(new_scene_rect);
  }

  void TOPPASWidget::wheelEvent(QWheelEvent * event)
  {
    zoom(event->angleDelta().y() < 0);
  }

  void TOPPASWidget::dragEnterEvent(QDragEnterEvent * event)
  {
    // TODO: test mime type/source? where?
    event->acceptProposedAction();
  }

  void TOPPASWidget::dragMoveEvent(QDragMoveEvent * event)
  {
    // TODO: test mime type/source? where?
    event->acceptProposedAction();
  }

  void TOPPASWidget::dropEvent(QDropEvent * event)
  {
    // TODO: test mime type/source? where?
    //std::cerr << "Drop Event with data:\n  " << String( event->mimeData()->formats().join("\n  ")) << "\n\n";

    if (event->mimeData()->hasUrls())
    {
      String filename = String(event->mimeData()->urls().front().toLocalFile());
      emit sendStatusMessage("loading drop file '" + filename + "' (press CRTL while dropping to insert into current window)", 0);
      // open pipeline in new window (or in current if CTRL is pressed)
      emit pipelineDroppedOnWidget(filename, event->keyboardModifiers() != Qt::ControlModifier);
    }
    else
    {
      QPointF scene_pos = mapToScene(event->pos());
      emit toolDroppedOnWidget(scene_pos.x(), scene_pos.y());
    }
    event->acceptProposedAction();
  }

  void TOPPASWidget::keyPressEvent(QKeyEvent * e)
  {
    if (e->key() == Qt::Key_C && e->modifiers() == Qt::ControlModifier)
    {
      scene_->copySelected();
      e->accept();
    }
    else if (e->key() == Qt::Key_X && e->modifiers() == Qt::ControlModifier)
    {
      scene_->copySelected();
      scene_->removeSelected();
      e->accept();
    }
    else if (e->key() == Qt::Key_V && e->modifiers() == Qt::ControlModifier)
    {
      scene_->paste();
      e->accept();
    }
    else if (e->key() == Qt::Key_Control)
    {
      setDragMode(QGraphicsView::RubberBandDrag);
      //color of hovering edge may change
      TOPPASEdge* hover_edge = scene_->getHoveringEdge();
      if (hover_edge)
      {
        hover_edge->update();
      }
      e->accept();
    }
    else if (e->key() == Qt::Key_Delete || e->key() == Qt::Key_Backspace)
    {
      scene_->removeSelected();
      e->accept();
    }
    else if (e->key() == Qt::Key_Plus)
    {
      zoom(false);
      e->accept();
    }
    else if (e->key() == Qt::Key_Minus)
    {
      zoom(true);
      e->accept();
    }
    else
    {
      e->ignore();
    }
  }

  void TOPPASWidget::keyReleaseEvent(QKeyEvent * e)
  {
    if (e->key() == Qt::Key_Control)
    {
      setDragMode(QGraphicsView::ScrollHandDrag);
      //color of hovering edge may change
      TOPPASEdge* hover_edge = scene_->getHoveringEdge();
      if (hover_edge)
      {
        hover_edge->update();
      }
      e->accept();
    }
  }

  void TOPPASWidget::leaveEvent(QEvent * /*e*/)
  {

  }

  void TOPPASWidget::enterEvent(QEnterEvent* /*e*/)
  {
#ifndef Q_WS_MAC
    setFocus();
#endif
  }

  void TOPPASWidget::resizeEvent(QResizeEvent * /*event*/)
  {
// QGraphicsView::resizeEvent(event);
// if (scene_)
// {
// QRectF items_rect = scene_->itemsBoundingRect();
// scene_->setSceneRect(items_rect.united(mapToScene(viewport()->rect()).boundingRect()));
// }
  }

  void TOPPASWidget::closeEvent(QCloseEvent * e)
  {
    bool close = scene_->saveIfChanged();
    if (close)
    {
      e->accept();
    }
    else
    {
      e->ignore();
    }
  }

} //Namespace
