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

  EnhancedTabBar::~EnhancedTabBar()
  {

  }

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
