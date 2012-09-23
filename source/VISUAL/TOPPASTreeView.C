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

#include <OpenMS/VISUAL/TOPPASTreeView.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QDrag>
#include <QtGui/QApplication>
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

  TOPPASTreeView::~TOPPASTreeView()
  {

  }

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
