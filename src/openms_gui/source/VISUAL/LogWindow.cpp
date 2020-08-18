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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LogWindow.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QContextMenuEvent>
#include <QMenu>
#include <QTextEdit>

using namespace std;

namespace OpenMS
{
  LogWindow::LogWindow(QWidget* parent)
    : QTextEdit(parent)
  {
    const auto help = "Log Window<BR>"
      "<BR>Output from TOPP tools and other status information is shown here";
    this->setWhatsThis(help);
    this->setToolTip(help);

    setReadOnly(true);

    connect(this, SIGNAL(textChanged()), this, SLOT(trimText_()));
  }

  void LogWindow::contextMenuEvent(QContextMenuEvent* event)
  {
    QMenu context_menu;
    context_menu.addAction("Clear", [&]() {
      this->clear();
    });
    context_menu.exec(this->mapToGlobal(event->pos()));
  }

  void LogWindow::appendText(const QString& text)
  {
    moveCursor(QTextCursor::End, QTextCursor::MoveAnchor); // move cursor to end, since text is inserted at cursor
    insertPlainText(text);
    // show log window
    qobject_cast<QWidget*>(this->parent())->show();
  }

  void LogWindow::appendNewHeader(const LogWindow::LogState state, const String& heading, const String& body)
  {
    String state_string;
    switch (state)
    {
      case NOTICE: state_string = "NOTICE"; break;
      case WARNING: state_string = "WARNING"; break;
      case CRITICAL: state_string = "ERROR"; break;
    }

    // update log
    append("==============================================================================");
    append((DateTime::now().getTime() + " " + state_string + ": " + heading).toQString());
    append(body.toQString());

    //show log tool window
    qobject_cast<QWidget*>(parent())->show();
  }

  void LogWindow::addNewline()
  {
    append("");
  }


  void LogWindow::trimText_()
  {
    if (max_length_ <= 0) return;

    if (this->toPlainText().size() > max_length_)
    {
      this->setPlainText(this->toPlainText().right(max_length_ / 2));
      //std::cerr << "cut text to " << this->toPlainText().size() << "\n";
    }
  }
  int LogWindow::maxLength() const
  {
    return max_length_;
  }


  void LogWindow::setMaxLength(int max_length)
  {
    max_length_ = max_length;
  }

} //Namespace
