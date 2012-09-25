// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/CONCEPT/Types.h>


// Qt
#include <QtGui/QDragEnterEvent>
#include <QtGui/QDragMoveEvent>
#include <QtGui/QDropEvent>
#include <QtCore/QMimeData>
#include <QUrl>

using namespace std;

namespace OpenMS
{
  TOPPASWidget::TOPPASWidget(const Param & /*preferences*/, QWidget * parent, const String & tmp_path) :
    QGraphicsView(parent),
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

  TOPPASWidget::~TOPPASWidget()
  {
    emit aboutToBeDestroyed(window_id_);
  }

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
    zoom(event->delta() < 0);
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
      e->accept();
    }
  }

  void TOPPASWidget::leaveEvent(QEvent * /*e*/)
  {

  }

  void TOPPASWidget::enterEvent(QEvent * /*e*/)
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

  // from EnhancedTabBarWidgetInterface
  void TOPPASWidget::setWindowId(Int window_id)
  {
    window_id_ = window_id;
  }

  // from EnhancedTabBarWidgetInterface
  Int TOPPASWidget::getWindowId()
  {
    return window_id_;
  }

} //Namespace
