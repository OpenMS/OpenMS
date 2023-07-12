// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

// QT
#include <QtCore/QProcess>
#include <QtNetwork/QNetworkReply>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMdiArea>
#include <QtWidgets/QSplashScreen>

class QToolBar;
class QListWidget;
class QTextEdit;
class QMdiArea;
class QLabel;
class QWidget;
class QTreeWidget;
class QTreeWidgetItem;
class QWebView;
class QNetworkAccessManager;

namespace Ui
{
  class FLASHDeconvWizardBase;
}

namespace OpenMS
{
  /**
    @brief Main window of the FLASHDeconvWizard tool

  */
  class OPENMS_GUI_DLLAPI FLASHDeconvWizardBase : public QMainWindow, public DefaultParamHandler
  {
    Q_OBJECT

  public:
    /// Constructor
    FLASHDeconvWizardBase(QWidget* parent = nullptr);
    /// Destructor
    ~FLASHDeconvWizardBase() override;

    void showAboutDialog();

  protected:
    /// The current path (used for loading and storing).
    /// Depending on the preferences this is static or changes with the current window/layer.
    String current_path_;

    /// The path for temporary files
    String tmp_path_;

  private slots:
    // names created by QtCreator. Do not change them.
    void on_actionExit_triggered();
    void on_actionVisit_FLASHDeconv_homepage_triggered();
    void on_actionReport_new_issue_triggered();

  private:
    Ui::FLASHDeconvWizardBase* ui;
  }; // class

} // namespace OpenMS
