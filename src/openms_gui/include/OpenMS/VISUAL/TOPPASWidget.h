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

#ifndef OPENMS_VISUAL_TOPPASWIDGET_H
#define OPENMS_VISUAL_TOPPASWIDGET_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/EnhancedTabBarWidgetInterface.h>

#include <QtGui/QGraphicsView>

namespace OpenMS
{
  class TOPPASScene;
  class Param;

  /**
    @brief Widget visualizing and allowing to edit TOPP pipelines.

    This class is a subclass of QGraphicsView and visualizes a TOPPASScene.
    Several TOPPASWidgets can be opened in TOPPAS at the same time,
    managed by a QWorkspace.

        @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASWidget :
    public QGraphicsView,
    public EnhancedTabBarWidgetInterface
  {
    Q_OBJECT

public:

    /// Default constructor
    TOPPASWidget(const Param & preferences, QWidget * parent = nullptr, const String & tmp_path = "");

    /// Destructor
    ~TOPPASWidget() override;

    /// setter from EnhancedTabBarWidgetInterface
    void setWindowId(Int id) override;

    /// getter from EnhancedTabBarWidgetInterface
    Int getWindowId() override;

    /// Returns the scene
    TOPPASScene * getScene();
    /// Zooms in or out, depending on @p zoom_in
    void zoom(bool zoom_in);

signals:

    /// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
    void sendStatusMessage(std::string message, OpenMS::UInt time);
    /// Emitted when the cursor position changes (for displaying e.g. in status bar)
    void sendCursorStatus(double x = 0.0, double y = 0.0);
    /// Message about the destruction of this widget
    void aboutToBeDestroyed(int w_id);
    /// Emitted when a drop event occurs
    void toolDroppedOnWidget(double x = 0.0, double y = 0.0);
    /// Emitted when a drop event occurs
    void pipelineDroppedOnWidget(const String & filename, bool new_window);

protected:

    /// The scene visualized by this widget
    TOPPASScene * scene_;

    ///@name reimplemented QT events
    //@{
    void wheelEvent(QWheelEvent * event) override;
    void keyPressEvent(QKeyEvent * e) override;
    void keyReleaseEvent(QKeyEvent * e) override;
    void leaveEvent(QEvent * e) override;
    void enterEvent(QEvent * e) override;
    void dragEnterEvent(QDragEnterEvent * event) override;
    void dragMoveEvent(QDragMoveEvent * event) override;
    void dropEvent(QDropEvent * event) override;
    void resizeEvent(QResizeEvent * event) override;
    void closeEvent(QCloseEvent * e) override;
    //@}

    /// Widget id used as identifier
    Int window_id_;
  };
}

#endif
