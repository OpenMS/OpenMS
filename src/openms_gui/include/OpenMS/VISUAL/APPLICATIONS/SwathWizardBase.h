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
  class SwathWizardBase;
}

namespace OpenMS
{
  /**
    @brief Main window of the SwathWizard tool

  */
  class OPENMS_GUI_DLLAPI SwathWizardBase :
    public QMainWindow,
    public DefaultParamHandler
  {
    Q_OBJECT

public:
    /// Constructor
    SwathWizardBase(QWidget* parent = nullptr);
    /// Destructor
    ~SwathWizardBase() override;
 
    void showAboutDialog();

protected slots:


protected:
    /// Log output window
    //TOPPASLogWindow* log_;

    /// The current path (used for loading and storing).
    /// Depending on the preferences this is static or changes with the current window/layer.
    String current_path_;

    /// The path for temporary files
    String tmp_path_;

  private slots:
    // names created by QtCreator. Do not change them.
    void on_actionExit_triggered();
    void on_actionVisit_OpenSwath_homepage_triggered();
    void on_actionReport_new_issue_triggered();

  private:
    Ui::SwathWizardBase* ui;
    
  }; //class

} //namespace
