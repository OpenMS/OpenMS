// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

//QT
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMdiArea>
#include <QtWidgets/QButtonGroup>
#include <QtCore/QProcess>
#include <QtWidgets/QSplashScreen>
#include <QtNetwork/QNetworkReply>

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
  class FLASHQuantWizardBase;
}

namespace OpenMS
{
  /**
    @brief Main window of the FLASHDeconvWizard tool

  */
  class OPENMS_GUI_DLLAPI FLASHQuantWizardBase:
    public QMainWindow,
    public DefaultParamHandler
  {
    Q_OBJECT

public:
    /// Constructor
    FLASHQuantWizardBase(QWidget* parent = nullptr);
    /// Destructor
    ~FLASHQuantWizardBase() override;

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
    void on_actionVisit_FLASHQuant_homepage_triggered();
    void on_actionReport_new_issue_triggered();

  private:
    Ui::FLASHQuantWizardBase* ui;
  }; //class

} //namespace
