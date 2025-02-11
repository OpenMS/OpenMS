// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/VISUAL/APPLICATIONS/FLASHDeconvWizardBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/VISUAL/DIALOGS/FLASHDeconvTabWidget.h>
#include <ui_FLASHDeconvWizardBase.h>

// Qt
#include <QDesktopServices>
#include <QMessageBox>
#include <QSettings>

using namespace std;
using namespace OpenMS;

namespace OpenMS
{
  using namespace Internal;

  FLASHDeconvWizardBase::FLASHDeconvWizardBase(QWidget* parent) :
      QMainWindow(parent), DefaultParamHandler("FLASHDeconvWizardBase"),
      ui(new Ui::FLASHDeconvWizardBase)
  {
    ui->setupUi(this);
    QSettings settings("OpenMS", "FLASHDeconvWizard");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("windowState").toByteArray());
    setWindowTitle("FLASHDeconvWizard");
    setWindowIcon(QIcon(":/FLASHDeconvWizard.png"));

    FLASHDeconvTabWidget* cwidget = new FLASHDeconvTabWidget(this);
    setCentralWidget(cwidget);
  }


  FLASHDeconvWizardBase::~FLASHDeconvWizardBase()
  {
    delete ui;
  }


  void FLASHDeconvWizardBase::showAboutDialog()
  {
    QApplicationTOPP::showAboutDialog(this, "FLASHDeconvWizard");
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionExit_triggered()
  {
    QApplicationTOPP::exit();
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionVisit_FLASHDeconv_homepage_triggered()
  {
    const char* url = "https://openms.de/application/flashdeconv/";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

  void OpenMS::FLASHDeconvWizardBase::on_actionReport_new_issue_triggered()
  {
    const char* url = "https://github.com/OpenMS/OpenMS/issues";
    if (!QDesktopServices::openUrl(QUrl(url)))
    {
      QMessageBox::warning(nullptr, "Cannot open browser. Please check your default browser settings.", QString(url));
    }
  }

} // namespace OpenMS
