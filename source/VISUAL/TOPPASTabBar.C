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
