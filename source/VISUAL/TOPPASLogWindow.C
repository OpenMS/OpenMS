// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/VISUAL/TOPPASLogWindow.h>

#include <iostream>

#include <QtGui/QContextMenuEvent>
#include <QtGui/QMenu>
#include <QAction>

using namespace std;

namespace OpenMS
{

  TOPPASLogWindow::TOPPASLogWindow(QWidget * parent) :
    QTextEdit(parent),
    max_length_(-1)
  {
    // trim if required
    connect (this, SIGNAL(textChanged()), this, SLOT(trimText_())); 
  }

  TOPPASLogWindow::~TOPPASLogWindow()
  {
  }

  void TOPPASLogWindow::contextMenuEvent(QContextMenuEvent * e)
  {
    QMenu * menu = createStandardContextMenu();
    menu->addAction("Clear");
    QAction * selected = menu->exec(e->globalPos());
    if (selected && selected->text() == "Clear")
    {
      clear();
    }
    delete menu;
  }


  void TOPPASLogWindow::trimText_()
  {
     if (max_length_ <= 0) return;

     if (this->toPlainText().size() > max_length_)
     {
       this->setPlainText(this->toPlainText().right(max_length_/2));
       //std::cerr << "cut text to " << this->toPlainText().size() << "\n";
     }
  }
  int TOPPASLogWindow::maxLength() const
  {
     return (max_length_);
  }
  

  void TOPPASLogWindow::setMaxLength(int max_length)
  {
    max_length_ = max_length;
  }


} //namespace OpenMS
