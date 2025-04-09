// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
