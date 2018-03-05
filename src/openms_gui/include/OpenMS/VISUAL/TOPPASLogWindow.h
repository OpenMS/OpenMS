// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASLOGWINDOW_H
#define OPENMS_VISUAL_TOPPASLOGWINDOW_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//QT
#include <QtGui/QTextEdit>
class QContextMenuEvent;

namespace OpenMS
{
  /**
      @brief QTextEdit implementation with a "clear" button in the context menu
      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI TOPPASLogWindow :
    public QTextEdit
  {
    Q_OBJECT
    Q_PROPERTY(int max_length READ maxLength WRITE setMaxLength)

  public:
    /// Constructor
    TOPPASLogWindow(QWidget * parent = nullptr);
    /// Destructor
    ~TOPPASLogWindow() override;

    /// read max_length
    int maxLength() const;
    /// set max_length
    void setMaxLength(int max_length);

  protected:
    ///@name Reimplemented Qt events
    //@{
    void contextMenuEvent(QContextMenuEvent * e) override;
    //@}

    /// Members:
    int max_length_;  /// -1 by default, which means there is no maximum length

  protected slots:
    /// if text length reached max_length_, then delete prefix until length of text is 1/2 of max_length_
    void trimText_();

  };
}
#endif // OPENMS_VISUAL_TOPPASLOGWINDOW_H
