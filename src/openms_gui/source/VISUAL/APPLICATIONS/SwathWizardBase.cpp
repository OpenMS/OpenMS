// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/APPLICATIONS/SwathWizardBase.h>
#include <ui_SwathWizardBase.h>

#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/DIALOGS/SwathTabWidget.h>

//Qt
#include <QtCore/QDir>
#include <QDesktopServices>
#include <QMessageBox>
#include <QSettings>

using namespace std;
using namespace OpenMS;

namespace OpenMS
{
  using namespace Internal;

  SwathWizardBase::SwathWizardBase(QWidget* parent) :
    QMainWindow(parent),
    DefaultParamHandler("SwathWizardBase"),
    //clipboard_scene_(nullptr),
    ui(new Ui::SwathWizardBase)
  {
    ui->setupUi(this);
    QSettings settings("OpenMS", "SwathWizard");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    setWindowTitle("SwathWizard");
    setWindowIcon(QIcon(":/SwathWizard.png"));

    SwathTabWidget* cw_swath = new SwathTabWidget(this);
    setCentralWidget(cw_swath);
  }


  SwathWizardBase::~SwathWizardBase()
  {
    delete ui;
  }


  void SwathWizardBase::showAboutDialog()
  {
    QApplicationTOPP::showAboutDialog(this, "SwathWizard");
  }

  void OpenMS::SwathWizardBase::on_actionExit_triggered()
  {
      QApplicationTOPP::exit();
  }

  void OpenMS::SwathWizardBase::on_actionVisit_OpenSwath_homepage_triggered()
  {
    const char* url = "http://openswath.org";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

  void OpenMS::SwathWizardBase::on_actionReport_new_issue_triggered()
  {
    const char* url = "https://github.com/OpenMS/OpenMS/issues";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

} //namespace OpenMS





