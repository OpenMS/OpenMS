// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
