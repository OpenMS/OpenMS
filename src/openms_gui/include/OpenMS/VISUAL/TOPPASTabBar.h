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

#ifndef OPENMS_VISUAL_TOPPASTABBAR_H
#define OPENMS_VISUAL_TOPPASTABBAR_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//QT
#include <QtGui/QTabBar>
class QMouseEvent;
class QMimeData;

namespace OpenMS
{
  class String;

  /**
      @brief Convenience tab bar implementation

      This tab bar differs in several ways from the QTabBar:
      - you can close a tab by double-clicking it or through its context menu.
      - it works based on tab identifiers (a fixed id stored in tab data) rather than on tab indices, which might
        change when inserting or removing a tab.

      @ingroup Visual
  */
  class OPENMS_GUI_DLLAPI TOPPASTabBar :
    public QTabBar
  {
    Q_OBJECT
public:
    /// Constructor
    TOPPASTabBar(QWidget * parent = nullptr);
    /// Destructor
    ~TOPPASTabBar() override;

    /// Adds a new tab with the name @p text and the identifier @p id
    int addTab(const String & text, int id);
    /// Selects the tab with identifier @p id
    void setCurrentId(int id);

public slots:
    /// Remove the tab with identifier @p id
    void removeId(int id);

signals:
    /// Signal that indicates that the current tab changed
    void currentIdChanged(int id);
    /// Signal that indicates that the tab with identifier @p id is about to be removed (double click or context menu)
    void aboutToCloseId(int id);

    // /// Signal that is emitted, when a drag-and-drop action ends on a tab
    // void dropOnTab(const QMimeData* data, QWidget* source, int id);
    // /// Signal that is emitted, when a drag-and-drop action ends on the unused space on the right side of the tabs.
    // void dropOnWidget(const QMimeData* data, QWidget* source);

protected:
    ///@name Reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QMouseEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    //void dragEnterEvent(QDragEnterEvent* e);
    //void dropEvent(QDropEvent* e);
    //@}

    ///Returns the QTabBar index of the tab at position @p pos. If there is no tab at that position -1 is returned.
    int tabAt_(const QPoint & pos);

protected slots:
    /// Slot that translates the currentChanged(int) signal to currentIdChanged(int)
    void currentChanged_(int id);
  };

}
#endif // OPENMS_VISUAL_TOPPASTABBAR_H
