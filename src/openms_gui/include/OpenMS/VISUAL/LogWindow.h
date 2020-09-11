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

    /// appends text without adding linebreaks and shows the log-window
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

