// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/KERNEL/StandardTypes.h>

#include <QTextEdit>

namespace OpenMS
{
  class String;

  /**
    @brief A log window (QTextEdit) with convenience functions


  */
  class LogWindow
    : public QTextEdit
  {
    Q_OBJECT
    Q_PROPERTY(int max_length READ maxLength WRITE setMaxLength)

  public:
    ///Log message states
    enum LogState
    {
      NOTICE, ///< Notice
      WARNING, ///< Warning
      CRITICAL ///< Fatal error
    };

    /// Default constructor
    LogWindow(QWidget* parent);

    /// appends text without adding line breaks and shows the log-window
    void appendText(const QString& text);

    /// appends a new block with @p heading and a @p body
    void appendNewHeader(const LogState state, const String& heading, const String& body);
    
    /// appends a line break (same as append(""))
    void addNewline();

    /// read max_length
    int maxLength() const;
    /// set max_length
    void setMaxLength(int max_length);

  signals:

  protected slots:
    /// if text length reached max_length_, then delete prefix until length of text is 1/2 of max_length_
    void trimText_();

  private:
    void contextMenuEvent(QContextMenuEvent* event) override;
    int max_length_ { -1 };  ///< -1 by default, which means there is no maximum length
  };

} //namespace

