// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/EnhancedWorkspace.h>
#include <QtCore/QMimeData>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDropEvent>

#include <QMdiSubWindow>
#include <QtCore/QStringList>

namespace OpenMS
{

  EnhancedWorkspace::EnhancedWorkspace(QWidget* parent) :
    QMdiArea(parent)
  {
    setAcceptDrops(true);
  }

  EnhancedWorkspace::~EnhancedWorkspace()
  {
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
